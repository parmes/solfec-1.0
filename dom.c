/*
 * dom.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * a domain gathers bodies and constraints
 */

/* This file is part of Solfec.
 * Solfec is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * Solfec is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with Solfec. If not, see <http://www.gnu.org/licenses/>. */

#include <string.h>
#include <float.h>
#include "alg.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "set.h"
#include "dom.h"
#include "goc.h"
#include "tmr.h"
#include "err.h"

#define CONBLK 512 /* constraints memory block size */
#define MAPBLK 128 /* map items memory block size */
#define SETBLK 512 /* set items memory block size */
#define SIZE (HASH3D+1) /* aabb timing tables size */
#define CON(ptr) ((CON*)(ptr))

typedef struct private_data DATA;

struct private_data
{
  double aabb_timings [SIZE],
	 aabb_limits [SIZE+1];

  int aabb_counter;

  BOXALG aabb_algo;
};

/* create private data */
static DATA* data_create (void)
{
  double part;
  DATA *data;

  ERRMEM (data = malloc (sizeof (DATA)));

  part = 1.0 / (double) SIZE;
  data->aabb_limits [0] = 0.0;
  data->aabb_counter = 0;

  return data;
}

/* free private data */
static void data_destroy (DATA *data)
{
  free (data);
}

/* calculate orthonormal
 * base storing it in column-wise manner
 * in 'loc' with 'n' as the last column */
static void localbase (double *n, double *loc)
{
  double  len,
    e [2][3] = { {1., 0., 0.},
                 {0., 1., 0.} };
  loc [6] = n [0];
  loc [7] = n [1];
  loc [8] = n [2];
  PRODUCT (e [0], n, loc);
  if ((len = LEN (loc)) < GEOMETRIC_EPSILON)
  {
    PRODUCT (e [1], n, loc);
    len = LEN (loc);
  }
  loc [0] /= len;
  loc [1] /= len;
  loc [2] /= len;
  PRODUCT (loc, n, (loc + 3));
}

/* insert a new constrait between two bodies */
static CON* insert (DOM *dom, BODY *master, BODY *slave)
{
  CON *con;

  /* make sure one of the body pointers was passed, tollerate
   * one of them to be NULL, as this indicates a single body constraint */
  ASSERT_DEBUG (master || slave, "At least one body pointer must be passed");

  ERRMEM (con = MEM_Alloc (&dom->conmem));
  con->master = master;
  con->slave = slave;

  /* add to the body constraint adjacency => NOTE this dupliciates
   * the dual graph structure, which will be of use in the distributed
   * memory parallel implementation (where partitionings of bodies and
   * contact points are separate) */
  SET_Insert (&dom->setmem, &master->con, con, NULL);
  if (slave) SET_Insert (&dom->setmem, &slave->con, con, NULL);
 
  /* add to list */
  con->next = dom->con;
  if (dom->con)
    dom->con->prev = con;
  dom->con = con;
  dom->ncon ++;

  /* insert into local dynamics */
  con->dia = LOCDYN_Insert (dom->ldy, con, master, slave);

  /* constraint identifier */
  con->id = dom->cid ++;
  MAP_Insert (&dom->mapmem, &dom->idc, (void*)con->id, con, NULL);

  return con;
}

/* fastest box overlap algorithm */
static BOXALG aabb_algorithm (DOM *dom)
{
  DATA *data = dom->data;
  double num, *tim, *lim;
  int i;

  if (data->aabb_counter < SIZE)
    data->aabb_algo = data->aabb_counter ++; /* at first test all algorithms */
  else /* when all tested */
  {
    tim = data->aabb_timings;
    lim = data->aabb_limits;

    for (i = 0; i < SIZE; i ++) lim [i+1] = lim [i] + tim [i]; /* sum up */
    for (i = 1; i <= SIZE; i ++) lim [i] /= lim [SIZE]; /* normalize */

    num = DRAND(); /* random in [0, 1] */

    for (i = 0; i < SIZE; i ++)
    {
      if (num >= lim [i] && num < lim [i+1])
      {
	data->aabb_algo = i;
	return i;
      }
    }
  }

  return (data->aabb_algo = HYBRID);
}

/* update timing related data */
static void aabb_timing (DOM *dom, double timing)
{
  DATA *data = dom->data;

  data->aabb_timings [data->aabb_algo] = timing;
}

/* insert a contact into the constraints set */
static CON* insert_contact (DOM *dom, BOX *one, BOX *two, BODY *master, BODY *slave,
  void *mgobj, GOBJ mkind, SHAPE *mshp, void *sgobj, GOBJ skind, SHAPE *sshp, double *spampnt,
  double *spaspnt, double *normal, double gap, double area, SURFACE_MATERIAL *mat, short paircode)
{
  CON *con;

  con = insert (dom, master, slave);
  con->kind = CONTACT;
  con->mgobj = mgobj;
  con->mkind = mkind;
  con->mshp = mshp;
  con->sgobj = sgobj;
  con->skind = skind;
  con->sshp = sshp;
  con->one = one;
  con->two = two;
  COPY (spampnt, con->point);
  BODY_Ref_Point (master, mshp, mgobj, spampnt, con->mpnt); /* referential image */
  BODY_Ref_Point (slave, sshp, sgobj, spaspnt, con->spnt); /* ... */
  localbase (normal, con->base);
  con->gap = gap;
  con->area = area;
  con->paircode = paircode;
  con->state |= SURFACE_MATERIAL_Transfer (dom->time, mat, &con->mat); /* transfer surface pair data from the database to the local variable */
  con->state |= CON_NEW;  /* mark as newly created */
  return con;
}

/* box overlap creation callback */
static void* overlap_create (DOM *dom, BOX *one, BOX *two)
{
  double onepnt [3], twopnt [3], normal [3], gap, area;
  int state, spair [2];
  short paircode;
  SURFACE_MATERIAL *mat;

  state = gobjcontact (
    CONTACT_DETECT, GOBJ_Pair_Code (one, two),
    one->sgp->shp, one->sgp->gobj,
    two->sgp->shp, two->sgp->gobj,
    onepnt, twopnt, normal, &gap, &area, spair);

  if (state)
  {
    ASSERT (gap > dom->depth, ERR_DOM_DEPTH);

    /* set surface pair data if there was a contact */
    mat = SPSET_Find (dom->sps, spair [0], spair [1]);
  }

  switch (state)
  {
    case 1: /* first body is the master */
    {
      paircode = GOBJ_Pair_Code (one, two);
      return insert_contact (dom, one, two, one->body, two->body, one->sgp->gobj, one->kind, one->sgp->shp,
	two->sgp->gobj, two->kind, two->sgp->shp, onepnt, twopnt, normal, gap, area, mat, paircode);
    }
    case 2: /* second body is the master */
    {
      paircode = GOBJ_Pair_Code (two, one);
      return insert_contact (dom, one, two, two->body, one->body, two->sgp->gobj, two->kind, two->sgp->shp,
	one->sgp->gobj, one->kind, one->sgp->shp, twopnt, onepnt, normal, gap, area, mat, paircode);
    }
  }

  return NULL; /* no contact found */
}

/* box verlap release callback */
static void overlap_release (DOM *dom, BOX *one, BOX *two, CON *con)
{
  DOM_Remove_Constraint (dom, con);
}

/* update contact data */
void update_contact (DOM *dom, CON *con)
{
  double mpnt [3], spnt [3], normal [3];
  int state, spair [2] = {con->mat.surf1,
                          con->mat.surf2};

  if (con->state & CON_NEW)
  {
    con->state &= ~CON_NEW;
    return; /* new contacts are up to date */
  }

  /* current spatial points and normal */
  BODY_Cur_Point (con->master, con->mshp, con->mgobj, con->mpnt, mpnt);
  BODY_Cur_Point (con->slave, con->sshp, con->sgobj, con->spnt, spnt);
  COPY (con->base+6, normal);

  /* update contact data => during an update 'master' and 'slave' relation does not change */
  state = gobjcontact (
    CONTACT_UPDATE, con->paircode,
    con->mshp, con->mgobj, con->sshp, con->sgobj,
    mpnt, spnt, normal, &con->gap, /* 'mpnt' and 'spnt' are updated here */
    &con->area, spair); /* surface pair might change though */

  if (state == 0) /* remove contact */
  {
    AABB_Break_Adjacency (dom->aabb, con->one, con->two); /* box overlap will be re-detected */
    DOM_Remove_Constraint (dom, con);
  }
  else
  {
    ASSERT (con->gap > dom->depth, ERR_DOM_DEPTH);

    COPY (mpnt, con->point);
    BODY_Ref_Point (con->master, con->mshp, con->mgobj, mpnt, con->mpnt);
    BODY_Ref_Point (con->slave, con->sshp, con->sgobj, spnt, con->spnt);
    localbase (normal, con->base);
    if (state > 1) /* surface pair has changed */
    {
      SURFACE_MATERIAL *mat = SPSET_Find (dom->sps, spair [0], spair [1]); /* find new surface pair description */
      con->state |= SURFACE_MATERIAL_Transfer (dom->time, mat, &con->mat); /* transfer surface pair data from the database to the local variable */
    }
  }
}

/* update fixed point data */
void update_fixpnt (DOM *dom, CON *con)
{
  BODY_Cur_Point (con->master, con->mshp, con->mgobj, con->mpnt, con->point);
}

/* update fixed direction data */
void update_fixdir (DOM *dom, CON *con)
{
  BODY_Cur_Point (con->master, con->mshp, con->mgobj, con->mpnt, con->point);
}

/* update velocity direction data */
void update_velodir (DOM *dom, CON *con)
{
  VELODIR (con->Z) = TMS_Value (con->tms, dom->time + dom->step);
  BODY_Cur_Point (con->master, con->mshp, con->mgobj, con->mpnt, con->point);
}

/* update rigid link data */
void update_riglnk (DOM *dom, CON *con)
{
  double n [3],
	 m [3],
	 s [3],
	 len;

  if (con->master && con->slave)
  {
    BODY_Cur_Point (con->master, con->mshp, con->mgobj, con->mpnt, m);
    BODY_Cur_Point (con->slave, con->sshp, con->sgobj, con->spnt, s);
  }
  else /* master point to a spatial point link */
  {
    BODY_Cur_Point (con->master, con->mshp, con->mgobj, con->mpnt, m);
    COPY (con->spnt, s);
  }

  COPY (m, con->point);
  SUB (s, m, RIGLNK_VEC (con->Z));
  SUB (s, m, n);
  len = LEN (n);
  con->gap = len - RIGLNK_LEN(con->Z);
  len = 1.0 / len;
  SCALE (n, len);
  localbase (n, con->base);
}

/* tell whether the geometric objects are topologically adjacent */
static int gobj_adjacent (short paircode, void *aobj, void *bobj)
{
  switch (paircode)
  {
    case AABB_ELEMENT_ELEMENT: return ELEMENT_Adjacent (aobj, bobj);
    case AABB_CONVEX_CONVEX: return CONVEX_Adjacent (aobj, bobj);
    case AABB_SPHERE_SPHERE: return SPHERE_Adjacent (aobj, bobj);
  }

  return 0;
}

/* go over contact points and remove those whose corresponding
 * areas are much smaller than those of other points related to
 * objects directly topologically adjacent in their shape definitions */
static void sparsify_contacts (DOM *dom)
{
  double threshold = dom->threshold;
  SET *del, *itm;
  CON *con, *adj;
  OFFB *blk;
  MEM mem;
  int n;

  MEM_Init (&mem, sizeof (SET), SETBLK);

  for (del = NULL, con = dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT && con->state & CON_NEW) /* walk over all new contacts */
    {
      for (blk = con->dia->adj; blk; blk = blk->n)
      {
	adj = blk->dia->con;
	
	if (con->area < threshold * adj->area) /* check whether the area of the diagonal element is too small (this test is cheaper => let it go first) */
	{
	  if (con->master == adj->master && con->slave == adj->slave) /* identify contacts pair sharing the same pairs of bodies */
	  {
	    if (gobj_adjacent (GOBJ_Pair_Code_Ext (con->mkind, adj->mkind), con->mgobj, adj->mgobj)) /* check whether the geometric objects are topologically adjacent */
	       SET_Insert (&mem, &del, con, NULL); /* if so schedule the current contact for deletion */
	  }
	  else if (con->master == adj->slave && con->slave == adj->master)
	  {
	    if (gobj_adjacent (GOBJ_Pair_Code_Ext (con->mkind, adj->skind), con->mgobj, adj->sgobj))
	      SET_Insert (&mem, &del, con, NULL);
	  }
	  else if (con->slave == adj->master && con->master == adj->slave)
	  {
	    if (gobj_adjacent (GOBJ_Pair_Code_Ext (con->skind, adj->mkind), con->sgobj, adj->mgobj))
	      SET_Insert (&mem, &del, con, NULL);
	  }
	  else if (con->slave == adj->slave && con->master == adj->master)
	  {
	    if (gobj_adjacent (GOBJ_Pair_Code_Ext (con->skind, adj->skind), con->sgobj, adj->sgobj))
	      SET_Insert (&mem, &del, con, NULL);
	  }
	}
      }
    }
  }

  /* now remove unwanted contacts */
  for (itm = SET_First (del), n = 0; itm; itm = SET_Next (itm), n ++)
  {
    con = itm->data;
    /* remove first from the box adjacency structure => otherwise box engine would try
     * to release this contact at a later point and that would cose memory corruption */
    AABB_Break_Adjacency (dom->aabb, con->one, con->two); /* box overlap will be re-detected */
    DOM_Remove_Constraint (dom, con); /* now remove from the domain */
  }

  /* report */
  if (dom->verbose) printf ("CONSTRAINTS: %d\nSPARSIFIED CONTACTS: %d\n",  dom->ncon, n);

  /* clean up */
  MEM_Release (&mem);

  /* TODO: possibly record recently sparsified contacts (gobj pairs) and do not insert them later for a while */
}

/* get geometrical object kind */
static int gobj_kind (SGP *sgp)
{
  switch (sgp->shp->kind)
  {
  case SHAPE_MESH: return GOBJ_ELEMENT;
  case SHAPE_CONVEX: return GOBJ_CONVEX;
  case SHAPE_SPHERE: return GOBJ_SPHERE;
  }

  ASSERT_DEBUG (0, "Invalid shape kind in gobj_kind");

  return 0;
}

/* get geometrical object extents update callback */
static BOX_Extents_Update extents_update (SGP *sgp)
{
  switch (sgp->shp->kind)
  {
  case SHAPE_MESH: return (BOX_Extents_Update) ELEMENT_Extents;
  case SHAPE_CONVEX: return (BOX_Extents_Update) CONVEX_Extents;
  case SHAPE_SPHERE: return (BOX_Extents_Update) SPHERE_Extents;
  }

  ASSERT_DEBUG (0, "Invalid shape kind in extents_update");

  return 0;
}

/* constraint kind string */
char* CON_Kind (CON *con)
{
  switch (con->kind)
  {
  case CONTACT: return "CONTACT";
  case FIXPNT: return "FIXPNT";
  case FIXDIR: return "FIXDIR";
  case VELODIR: return "VELODIR";
  case RIGLNK: return "RIGLNK";
  }

  return NULL;
}

/* create a domain */
DOM* DOM_Create (AABB *aabb, SPSET *sps, short dynamic, double step)
{
  DOM *dom;

  ERRMEM (dom = calloc (1, sizeof (DOM)));
  dom->aabb = aabb;
  dom->sps = sps;
  dom->dynamic = (dynamic == 1 ? 1 : 0);
  dom->step = step;
  dom->time = 0.0;

  MEM_Init (&dom->conmem, sizeof (CON), CONBLK);
  MEM_Init (&dom->mapmem, sizeof (MAP), MAPBLK);
  MEM_Init (&dom->setmem, sizeof (SET), SETBLK);
  dom->bid = 1;
  dom->lab = NULL;
  dom->idb = NULL;
  dom->bod = NULL;
  dom->nbod = 0;
  dom->delb = NULL;
  dom->newb = NULL;
  dom->cid = 1;
  dom->idc= NULL;
  dom->con = NULL;
  dom->ncon = 0;
  dom->prev = dom->next = NULL;
  dom->flags = 0;
  dom->threshold = 0.01;
  dom->depth = -DBL_MAX;
  ERRMEM (dom->ldy = LOCDYN_Create (dom));

  SET (dom->gravdir, 0);
  dom->gravval = NULL;

  SET (dom->extents, -DBL_MAX);
  SET (dom->extents + 3, DBL_MAX);

  dom->data = data_create ();

  dom->verbose = 0;

  return dom;
}

/* insert a body into the domain */
void DOM_Insert_Body (DOM *dom, BODY *bod)
{
  SGP *sgp, *sgpe;

  for (sgp = bod->sgp, sgpe = sgp + bod->nsgp; sgp < sgpe; sgp ++)
  {
    sgp->box = AABB_Insert (dom->aabb, bod, gobj_kind (sgp), sgp,
                            sgp->shp->data, extents_update (sgp));
  }

  if (bod->label) /* map labeled bodies */
    MAP_Insert (&dom->mapmem, &dom->lab, bod->label, bod, (MAP_Compare)strcmp);

  bod->id = dom->bid ++;
  MAP_Insert (&dom->mapmem, &dom->idb, (void*)bod->id, bod, NULL);

  bod->dom = dom;
  dom->nbod ++;

  bod->next = dom->bod;
  if (dom->bod) dom->bod->prev = bod;
  dom->bod = bod;

  if (dom->time > 0) SET_Insert (&dom->setmem, &dom->newb, bod, NULL);
}

/* remove a body from the domain */
void DOM_Remove_Body (DOM *dom, BODY *bod)
{
  SGP *sgp, *sgpe;

  for (sgp = bod->sgp, sgpe = sgp + bod->nsgp; sgp < sgpe; sgp ++)
  {
    AABB_Delete (dom->aabb, sgp->box);
  }

  /* remove all body related constraints */
  {
    SET *con = bod->con;
    bod->con = NULL; /* DOM_Remove_Constraint will try to remove the constraint
			from body constraints set, which is not nice if we try
			to iterate over the set at the same time => make it empty */
    
    for (SET *item = SET_First (con); item; item = SET_Next (item))
      DOM_Remove_Constraint (dom, item->data);

    SET_Free (&dom->setmem, &con); /* free body's constraint set */
  }

  if (bod->label) /* delete labeled body */
    MAP_Delete (&dom->mapmem, &dom->lab, bod->label, (MAP_Compare)strcmp);

  /* delete from id based map */
  MAP_Delete (&dom->mapmem, &dom->idb, (void*)bod->id, NULL);

  dom->nbod --;

  if (bod->prev) bod->prev->next = bod->next;
  else dom->bod = bod->next;
  if (bod->next) bod->next->prev = bod->prev;


  if (dom->time > 0) SET_Insert (&dom->setmem, &dom->delb, (void*) bod->id, NULL);
}

/* find labeled body */
BODY* DOM_Find_Body (DOM *dom, char *label)
{
  return MAP_Find (dom->lab, label, (MAP_Compare)strcmp);
}

/* fix a referential point of the body along all directions */
CON* DOM_Fix_Point (DOM *dom, BODY *bod, double *pnt)
{
  CON *con;

  con = insert (dom, bod, NULL);
  con->kind = FIXPNT;
  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  con->mgobj = SHAPE_Gobj (bod->shape, pnt, &con->mshp);
  IDENTITY (con->base);
  return con;
}

/* fix a referential point of the body along the spatial direction */
CON* DOM_Fix_Direction (DOM *dom, BODY *bod, double *pnt, double *dir)
{
  CON *con;

  con = insert (dom, bod, NULL);
  con->kind = FIXDIR;
  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  con->mgobj = SHAPE_Gobj (bod->shape, pnt, &con->mshp);
  localbase (dir, con->base);
  return con;
}

/* prescribe a velocity of the referential point along the spatial direction */
CON* DOM_Set_Velocity (DOM *dom, BODY *bod, double *pnt, double *dir, TMS *vel)
{
  CON *con;

  con = insert (dom, bod, NULL);
  con->kind = VELODIR;
  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  con->mgobj = SHAPE_Gobj (bod->shape, pnt, &con->mshp);
  localbase (dir, con->base);
  con->tms = vel;
  return con;
}

/* insert rigid link constraint between two (referential) points of bodies; if one of the body
 * pointers is NULL then the link acts between the other body and the fixed (spatial) point */
CON* DOM_Put_Rigid_Link (DOM *dom, BODY *master, BODY *slave, double *mpnt, double *spnt)
{
  double v [3], d;
  CON*con;

  if (!master) master = slave;
  ASSERT_DEBUG (master, "At least one body pointer must not be NULL");

  SUB (mpnt, spnt, v);
  d = LEN (v);
  
  if (d < GEOMETRIC_EPSILON) /* no point in keeping very short links */
  {
    con = insert (dom, master, slave);
    con->kind = FIXPNT;
    COPY (mpnt, con->point);
    COPY (mpnt, con->mpnt);
    COPY (spnt, con->spnt);
    con->mgobj = SHAPE_Gobj (master->shape, mpnt, &con->mshp);
    IDENTITY (con->base);

    if (master && slave) /* no contact between this pair */
    {
      int msgp, ssgp;
      SGP *sgp;

      msgp = SHAPE_Sgp (master->sgp, master->nsgp, mpnt);
      ssgp = SHAPE_Sgp (slave->sgp, slave->nsgp, spnt);
      sgp = &slave->sgp [ssgp];
      con->sgobj = sgp->gobj;
      con->sshp = sgp->shp;

      AABB_Exclude_Gobj_Pair (dom->aabb, master->id, msgp, slave->id, ssgp);
    }
  }
  else
  {
    con = insert (dom, master, slave);
    con->kind = RIGLNK;
    COPY (mpnt, con->point);
    COPY (mpnt, con->mpnt);
    COPY (spnt, con->spnt);
    RIGLNK_LEN (con->Z) = d; /* initial distance */
    con->mgobj = SHAPE_Gobj (master->shape, mpnt, &con->mshp);
    if (slave) con->sgobj = SHAPE_Gobj (slave->shape, spnt, &con->sshp);
    update_riglnk (dom, con); /* initial update */
  }
  
  return con;
}

/* remove a constraint from the domain */
void DOM_Remove_Constraint (DOM *dom, CON *con)
{
  /* remove from the body constraint adjacency  */
  SET_Delete (&dom->setmem, &con->master->con, con, NULL);
  if (con->slave) SET_Delete (&dom->setmem, &con->slave->con, con, NULL);

  /* remove from id-based map */
  MAP_Delete (&dom->mapmem, &dom->idc, (void*)con->id, NULL);

  /* remove from list */
  if (con->prev)
    con->prev->next = con->next;
  else dom->con = con->next;
  if (con->next)
    con->next->prev = con->prev;
  dom->ncon --;

  /* remove from local dynamics */
  LOCDYN_Remove (dom->ldy, con->dia);

  /* destroy passed data */
  MEM_Free (&dom->conmem, con);
}

/* set simulation scene extents */
void DOM_Extents (DOM *dom, double *extents)
{
  COPY6 (extents, dom->extents);
}

/* domain update initial half-step => bodies and constraints are
 * updated and the current local dynamic problem is returned */
LOCDYN* DOM_Update_Begin (DOM *dom)
{
  double time, step;
  TIMING timing;
  BODY *bod;
  CON *con;

  /* time and step */
  time = dom->time;
  step = dom->step;

  /* initialize bodies */
  if (time == 0.0)
  {
    if (dom->dynamic > 0)
      for (bod = dom->bod; bod; bod = bod->next)
	BODY_Dynamic_Init (bod, bod->scheme); /* integration scheme is set externally */
    else
      for (bod = dom->bod; bod; bod = bod->next)
	BODY_Static_Init (bod);
  }

  /* begin time integration */
  if (dom->dynamic)
    for (bod = dom->bod; bod; bod = bod->next)
      BODY_Dynamic_Step_Begin (bod, time, step);
  else
    for (bod = dom->bod; bod; bod = bod->next)
      BODY_Static_Step_Begin (bod, time, step);

  /* detect contacts */
  timerstart (&timing);

  AABB_Update (dom->aabb, aabb_algorithm (dom),
    dom, (BOX_Overlap_Create) overlap_create,
    (BOX_Overlap_Release) overlap_release);

  aabb_timing (dom, timerend (&timing));

  sparsify_contacts (dom);

  /* update all constraints */
  for (con = dom->con; con; con = con->next)
  {
    switch (con->kind)
    {
      case CONTACT: update_contact (dom, con); break;
      case FIXPNT:  update_fixpnt  (dom, con); break;
      case FIXDIR:  update_fixdir  (dom, con); break;
      case VELODIR: update_velodir (dom, con); break;
      case RIGLNK:  update_riglnk  (dom, con); break;
    }
  }

  /* output local dynamics */
  return dom->ldy;
}

/* domain update final half-step => once the local dynamic
 * problem has been solved (externally), motion of bodies
 * is updated with the help of new constraint reactions */
void DOM_Update_End (DOM *dom)
{
  double time, step, *de, *be;
  SET *del, *item;
  BODY *bod;

  /* time and step */
  time = dom->time;
  step = dom->step;

  /* end time integration */
  if (dom->dynamic)
    for (bod = dom->bod; bod; bod = bod->next)
      BODY_Dynamic_Step_End (bod, time, step);
  else
    for (bod = dom->bod; bod; bod = bod->next)
      BODY_Static_Step_End (bod, time, step);

  /* advance time */
  dom->time += step;
  dom->step = step;

  /* erase bodies outside of scene extents */
  for (bod = dom->bod, de = dom->extents, del = NULL; bod; bod = bod->next)
  {
    be = bod->extents;

    if (be [3] < de [0] ||
	be [4] < de [1] ||
	be [5] < de [2] ||
	be [0] > de [3] ||
	be [1] > de [4] ||
	be [2] > de [5]) SET_Insert (&dom->setmem, &del, bod, NULL); /* insert into deletion set */
  }

  for (item = SET_First (del); item; item = SET_Next (item)) /* remove bodies falling out of the scene extents */
  {
    DOM_Remove_Body (dom, item->data);
    BODY_Destroy (item->data);
  }

  SET_Free (&dom->setmem, &del); /* free up deletion set */
}

/* write constraint state */
static void write_constraint (CON *con, PBF *bf)
{
  int kind = con->kind;

  PBF_Uint (bf, &con->id, 1);
  PBF_Int (bf, &kind, 1);

  PBF_Double (bf, con->R, 3);
  PBF_Double (bf, con->point, 3);
  PBF_Double (bf, con->base, 9);
  PBF_Double (bf, &con->area, 1);
  PBF_Double (bf, &con->gap, 1);

  int count = 1;

  if (con->slave) count = 2;

  PBF_Int (bf, &count, 1);
  PBF_Uint (bf, &con->master->id, 1);
  if (con->slave) PBF_Uint (bf, &con->slave->id, 1);

  if (kind == CONTACT) PBF_String (bf, &con->mat.label);

  if (kind == RIGLNK || kind == VELODIR) PBF_Double (bf, con->Z, Z_SIZE);
}

/* read constraint state */
static CON* read_constraint (DOM *dom, PBF *bf)
{
  CON *con;
  int kind;

  ERRMEM (con = MEM_Alloc (&dom->conmem));

  PBF_Uint (bf, &con->id, 1);
  PBF_Int (bf, &kind, 1);
  con->kind = kind;

  PBF_Double (bf, con->R, 3);
  PBF_Double (bf, con->point, 3);
  PBF_Double (bf, con->base, 9);
  PBF_Double (bf, &con->area, 1);
  PBF_Double (bf, &con->gap, 1);

  unsigned int id;
  int count;

  PBF_Int (bf, &count, 1);
  PBF_Uint (bf, &id, 1);
  ASSERT_DEBUG_EXT (con->master = MAP_Find (dom->idb, (void*)id, NULL), "Invalid master id");
  if (count == 2)
  {
    PBF_Uint (bf, &id, 1);
    ASSERT_DEBUG_EXT (con->slave = MAP_Find (dom->idb, (void*)id, NULL), "Invalid slave id");
  }

  if (kind == CONTACT)
  {
    char buf [PBF_MAXSTRING], *label = buf;
    SURFACE_MATERIAL *mat;

    PBF_String (bf, &label);
    ASSERT_DEBUG_EXT (mat = SPSET_Find_Label (dom->sps, label), "Invalid material label");
    con->state |= SURFACE_MATERIAL_Transfer (dom->time, mat, &con->mat);
  }

  if (kind == RIGLNK || kind == VELODIR) PBF_Double (bf, con->Z, Z_SIZE);

  return con;
}

/* write domain state */
void DOM_Write_State (DOM *dom, PBF *bf)
{
  /* mark domain output */

  PBF_Label (bf, "DOM");

  /* write time step */

  PBF_Label (bf, "STEP");

  PBF_Double (bf, &dom->step, 1);

  /* write contacts */

  PBF_Label (bf, "CONS");
 
  PBF_Int (bf, &dom->ncon, 1);

  for (CON *con = dom->con; con; con = con->next)
  {
    write_constraint (con, bf);
  }

  /* write ids of bodies that have been deleted and empty the deleted bodies ids set */

  unsigned int id;
  SET *item;
  int size;

  PBF_Label (bf, "DELBODS");

  size = SET_Size (dom->delb);

  PBF_Int (bf, &size, 1);

  for (item = SET_First (dom->delb); item; item = SET_Next (item))
  {
    id = (int)item->data;
    PBF_Uint (bf, &id, 1);
  }

  SET_Free (&dom->setmem, &dom->delb);

  /* write complete data of newly created bodies and empty the newly created bodies set */

  int dsize = 0, doubles;
  double *d = NULL;

  int isize = 0, ints, *i = NULL;

  PBF_Label (bf, "NEWBODS");

  size = SET_Size (dom->newb);

  PBF_Int (bf, &size, 1);

  for (item = SET_First (dom->newb); item; item = SET_Next (item))
  {
    doubles = ints = 0;

    BODY_Pack (item->data, &dsize, &d, &doubles, &isize, &i, &ints);

    PBF_Int (bf, &doubles, 1);
    PBF_Double (bf, d, doubles);

    PBF_Int (bf, &ints, 1);
    PBF_Int (bf, i, ints);
  }

  free (d);
  free (i);

  SET_Free (&dom->setmem, &dom->newb);

  /* write regular bodies (this also includes states of newly created ones) */

  PBF_Label (bf, "BODS");

  PBF_Int (bf, &dom->nbod, 1);

  for (BODY *bod = dom->bod; bod; bod = bod->next)
  {
    PBF_Uint (bf, &bod->id, 1);
    BODY_Write_State (bod, bf);
  }
}

/* read domain state */
void DOM_Read_State (DOM *dom, PBF *bf)
{
  if (PBF_Label (bf, "DOM"))
  {
    /* read time step */

    ASSERT (PBF_Label (bf, "STEP"), ERR_FILE_FORMAT);

    PBF_Double (bf, &dom->step, 1);

    /* read contacts */

    ASSERT (PBF_Label (bf, "CONS"), ERR_FILE_FORMAT);

    MAP_Free (&dom->mapmem, &dom->idc);
    MEM_Release (&dom->conmem);
    dom->con = NULL;
   
    PBF_Int (bf, &dom->ncon, 1);

    for (int n = 0; n < dom->ncon; n ++)
    {
      CON *con;
      
      con = read_constraint (dom, bf);
      MAP_Insert (&dom->mapmem, &dom->idc, (void*)con->id, con, NULL);
      con->next = dom->con;
      if (dom->con) dom->con->prev = con;
      dom->con = con;
    }

    /* read ids of bodies that need to be deleted and remove them from all containers */

    unsigned int id;
    int size, n;

    ASSERT (PBF_Label (bf, "DELBODS"), ERR_FILE_FORMAT);

    PBF_Int (bf, &size, 1);

    for (n = 0; n < size; n ++)
    {
      BODY *bod;

      PBF_Uint (bf, &id, 1);

      ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*)id, NULL), "Invalid body id");

      if (bod->label) MAP_Delete (&dom->mapmem, &dom->lab, bod->label, (MAP_Compare) strcmp);
      MAP_Delete (&dom->mapmem, &dom->idb, (void*)id, NULL);
      if (bod->next) bod->next->prev = bod->prev;
      if (bod->prev) bod->prev->next = bod->next;
      else dom->bod = bod->next;
      dom->nbod --;
      BODY_Destroy (bod);
    }

    /* read complete data of newly created bodies and insert them into all containers */

    int dpos, doubles;
    double *d = NULL;

    int ipos, ints, *i = NULL;

    ASSERT (PBF_Label (bf, "NEWBODS"), ERR_FILE_FORMAT);

    PBF_Int (bf, &size, 1);

    for (n = 0; n < size; n ++)
    {
      BODY *bod;

      PBF_Int (bf, &doubles, 1);
      ERRMEM (d  = realloc (d, doubles * sizeof (double)));
      PBF_Double (bf, d, doubles);

      PBF_Int (bf, &ints, 1);
      ERRMEM (i  = realloc (i, ints * sizeof (int)));
      PBF_Int (bf, i, ints);

      dpos = ipos = 0;

      bod = BODY_Unpack (dom->owner, &dpos, d, doubles, &ipos, i, ints);

      if (bod->label) MAP_Insert (&dom->mapmem, &dom->lab, bod->label, bod, (MAP_Compare) strcmp);
      MAP_Insert (&dom->mapmem, &dom->idb, (void*)bod->id, bod, NULL);
      bod->next = dom->bod;
      if (dom->bod) dom->bod->prev = bod;
      dom->bod = bod;
      dom->nbod ++;
    }

    free (d);
    free (i);

    /* read regular bodies */

    ASSERT (PBF_Label (bf, "BODS"), ERR_FILE_FORMAT);

    int nbod;

    PBF_Int (bf, &nbod, 1);

    for (int n = 0; n < nbod; n ++)
    {
      unsigned int id;
      BODY *bod;

      PBF_Uint (bf, &id, 1);
      ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*)id, NULL), "Body id invalid");
      BODY_Read_State (bod, bf);
    }
  }
}

/* release memory */
void DOM_Destroy (DOM *dom)
{
  BODY *bod, *next;
 
  for (bod = dom->bod; bod; bod = next)
  {
    next = bod->next;
    BODY_Destroy (bod);
  }

  LOCDYN_Destroy (dom->ldy);

  MEM_Release (&dom->conmem);
  MEM_Release (&dom->setmem);

  if (dom->gravval)
    TMS_Destroy (dom->gravval);

  data_destroy (dom->data);

  free (dom);
}
