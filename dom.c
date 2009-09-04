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
#include <limits.h>
#include <float.h>
#include "sol.h"
#include "alg.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "set.h"
#include "dom.h"
#include "goc.h"
#include "tmr.h"
#include "err.h"

#if MPI
#include "put.h"
#include "pck.h"
#include "com.h"
#include "tag.h"
#endif

#define CONBLK 512 /* constraints memory block size */
#define MAPBLK 128 /* map items memory block size */
#define SETBLK 512 /* set items memory block size */
#define SIZE (HASH3D+1) /* aabb timing tables size */

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
  data->aabb_algo = 0;

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
static CON* insert (DOM *dom, BODY *master, BODY *slave, short locdyn)
{
  CON *con;

  /* make sure one of the body pointers was passed, tollerate
   * one of them to be NULL, as this indicates a single body constraint */
  ASSERT_DEBUG (master || slave, "At least one body pointer must be passed");

  ERRMEM (con = MEM_Alloc (&dom->conmem));
  con->master = master;
  con->slave = slave;

  /* add to the body constraint adjacency */
  SET_Insert (&dom->setmem, &master->con, con, NULL);
  if (slave) SET_Insert (&dom->setmem, &slave->con, con, NULL);
 
  /* add to list */
  con->next = dom->con;
  if (dom->con)
    dom->con->prev = con;
  dom->con = con;
  dom->ncon ++;

  /* constraint identifier */
#if MPI
  if (!dom->noid)
  {
    if (SET_Size (dom->sparecid))
    {
      SET *item;
     
      item = SET_First (dom->sparecid);
      con->id = (unsigned int) (long) item->data; /* use a previously freed id */
      SET_Delete (&dom->setmem, &dom->sparecid, item->data, NULL);
    }
    else con->id = dom->cid, dom->cid += dom->ncpu; /* every ncpu number */
  }
#else
  con->id = dom->cid ++;
#endif

  MAP_Insert (&dom->mapmem, &dom->idc, (void*) (long) con->id, con, NULL);

  /* insert into local dynamics */
  if (locdyn) con->dia = LOCDYN_Insert (dom->ldy, con, master, slave);

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

  return data->aabb_algo;
}

/* update timing related data */
static void aabb_timing (DOM *dom, double timing)
{
  DATA *data = dom->data;

  data->aabb_timings [data->aabb_algo] = timing;
}

/* insert a contact into the constraints set */
static CON* insert_contact (DOM *dom, BOX *mbox, BOX *sbox, BODY *master, BODY *slave,
  void *mgobj, GOBJ mkind, SHAPE *mshp, void *sgobj, GOBJ skind, SHAPE *sshp, double *spampnt,
  double *spaspnt, double *normal, double gap, double area, SURFACE_MATERIAL *mat, short paircode)
{
  CON *con;

  con = insert (dom, master, slave, 0); /* 0: do not insert into LOCDYN yet, only after sparsification */
  con->kind = CONTACT;
  con->mgobj = mgobj;
  con->mkind = mkind;
  con->mshp = mshp;
  con->sgobj = sgobj;
  con->skind = skind;
  con->sshp = sshp;
  con->mbox = mbox;
  con->sbox = sbox;
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
#if MPI
    int proc;

    ASSERT (Zoltan_LB_Point_Assign (dom->aabb->zol, state == 1 ? onepnt : twopnt, &proc) == ZOLTAN_OK, ERR_ZOLTAN);
    if (proc != dom->rank) return NULL; /* insert contacts located in this AABB partition */
#endif

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
      return insert_contact (dom, two, one, two->body, one->body, two->sgp->gobj, two->kind, two->sgp->shp,
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
    /* insert into local dynamics */
    con->dia = LOCDYN_Insert (dom->ldy, con, con->master, con->slave);
    /* invalidate newness */
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

  if (state)
  {
#if MPI
    int proc;

    ASSERT (Zoltan_LB_Point_Assign (dom->aabb->zol, mpnt, &proc) == ZOLTAN_OK, ERR_ZOLTAN);
    if (proc != dom->rank) goto del; /* delete contacts located outside of this AABB partition */
#endif

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
  else /* remove contact */
  {
#if MPI
del:
#endif
    AABB_Break_Adjacency (dom->aabb, con->mbox, con->sbox); /* box overlap will be re-detected */
    DOM_Remove_Constraint (dom, con);
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
  MEM mem;
  int n;

  MEM_Init (&mem, sizeof (SET), SETBLK);

  for (del = NULL, con = dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT && con->state & CON_NEW) /* walk over all new contacts */
    {
      SET *set [2] = {con->master->con, con->slave->con};

      for (n = 0; n < 2; n ++) for (itm = SET_First (set [n]); itm; itm = SET_Next (itm))
      {
	adj = itm->data;

	if (adj == con) continue;
	
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
    AABB_Break_Adjacency (dom->aabb, con->mbox, con->sbox); /* box overlap will be re-detected */
    DOM_Remove_Constraint (dom, con); /* now remove from the domain */
  }

  dom->nspa = n; /* record the number of sparsified contacts */

  /* clean up */
  MEM_Release (&mem);
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

  if (kind == RIGLNK || kind == VELODIR) PBF_Double (bf, con->Z, DOM_Z_SIZE);
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
  ASSERT_DEBUG_EXT (con->master = MAP_Find (dom->idb, (void*) (long) id, NULL), "Invalid master id");
  if (count == 2)
  {
    PBF_Uint (bf, &id, 1);
    ASSERT_DEBUG_EXT (con->slave = MAP_Find (dom->idb, (void*) (long) id, NULL), "Invalid slave id");
  }

  if (kind == CONTACT)
  {
    char buf [PBF_MAXSTRING], *label = buf;
    SURFACE_MATERIAL *mat;

    PBF_String (bf, &label);
    ASSERT_DEBUG_EXT (mat = SPSET_Find_Label (dom->sps, label), "Invalid material label");
    con->state |= SURFACE_MATERIAL_Transfer (dom->time, mat, &con->mat);
  }

  if (kind == RIGLNK || kind == VELODIR) PBF_Double (bf, con->Z, DOM_Z_SIZE);

  return con;
}

#if MPI
typedef struct conaux CONAUX;

/* auxiliary
 * constraint */
struct conaux
{
  int id,
      kind,
      master,
      slave;
  double R [3];
  double vec [2][3];
  TMS *tms;
};

/* pack non-contact constraint */
static void pack_constraint (CON *con, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, con->id);
  pack_int (isize, i, ints, con->kind);
  pack_doubles (dsize, d, doubles, con->R, 3);
  pack_int (isize, i, ints, con->master->id);

  switch (con->kind)
  {
  case FIXPNT:
    pack_doubles (dsize, d, doubles, con->mpnt, 3);
    break;
  case FIXDIR:
    pack_doubles (dsize, d, doubles, con->mpnt, 3);
    pack_doubles (dsize, d, doubles, con->base + 6, 3);
    break;
  case VELODIR:
    pack_doubles (dsize, d, doubles, con->mpnt, 3);
    pack_doubles (dsize, d, doubles, con->base + 6, 3);
    TMS_Pack (con->tms, dsize, d, doubles, isize, i, ints);
    break;
  case RIGLNK:
    pack_int (isize, i, ints, con->slave->id);
    pack_doubles (dsize, d, doubles, con->mpnt, 3);
    pack_doubles (dsize, d, doubles, con->spnt, 3);
    break;
  case CONTACT:
    ASSERT_DEBUG (0, "Trying to pack a contact constraint");
    break;
  }
}

/* unpack non-contact constraint */
CONAUX* unpack_constraint (MEM *mem, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  CONAUX *aux;

  ERRMEM (aux = MEM_Alloc (mem));

  aux->id = unpack_int (ipos, i, ints);
  aux->kind = unpack_int (ipos, i, ints);
  unpack_doubles (dpos, d, doubles, aux->R, 3);
  aux->master = unpack_int (ipos, i, ints);

  switch (aux->kind)
  {
  case FIXPNT:
    unpack_doubles (dpos, d, doubles, aux->vec [0], 3);
    break;
  case FIXDIR:
    unpack_doubles (dpos, d, doubles, aux->vec [0], 3);
    unpack_doubles (dpos, d, doubles, aux->vec [1], 3);
    break;
  case VELODIR:
    unpack_doubles (dpos, d, doubles, aux->vec [0], 3);
    unpack_doubles (dpos, d, doubles, aux->vec [1], 3);
    aux->tms = TMS_Unpack (dpos, d, doubles, ipos, i, ints);
    break;
  case RIGLNK:
    aux->slave = unpack_int (ipos, i, ints);
    unpack_doubles (dpos, d, doubles, aux->vec [0], 3);
    unpack_doubles (dpos, d, doubles, aux->vec [1], 3);
    break;
  }

  return aux;
}

/* number of bodies */
static int body_count (DOM *dom, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return dom->nbod;
}

/* list of body identifiers */
static void body_list (DOM *dom, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr)
{
  BODY *bod;
  int i;
  
  for (bod = dom->bod, i = 0; bod; i ++, bod = bod->next)
  {
    global_ids [i * num_gid_entries] = bod->id;
    obj_wgts [i * wgt_dim] = bod->dofs;
  }

  *ierr = ZOLTAN_OK;
}

/* number of spatial dimensions */
static int dimensions (DOM *dom, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}

/* list of body extent midpoints */
static void midpoints (DOM *dom, int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int num_dim, double *geom_vec, int *ierr)
{
  unsigned int id;
  double *e, *v;
  BODY *bod;
  int i;

  for (i = 0, bod = dom->bod; i < num_obj; i ++, bod = bod->next)
  {
    id = global_ids [i * num_gid_entries];

    if (bod && bod->id != id) /* Zoltan changed the order of bodies in the list */
    {
      ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) id, NULL), "Invalid body id");
    }

    e = bod->extents;
    v = &geom_vec [i* num_dim];
    MID (e, e+3, v);
  }

  *ierr = ZOLTAN_OK;
}

/* insert a child body into the domain */
static void insert_child (DOM *dom, BODY *bod)
{
  /* insert into id based map */
  MAP_Insert (&dom->mapmem, &dom->children, (void*) (long) bod->id, bod, NULL);

  bod->dom = dom;

  /* note, that the bounding boxes related to this
   * body migrate independently in the AABB module */
}

/* remove a child body from the domain */
static void remove_child (DOM *dom, BODY *bod)
{
  AABB_Detach_Body (dom->aabb, bod);

  /* remove child related constraints */
  {
    SET *con = bod->con;
    bod->con = NULL; /* DOM_Remove_Constraint will try to remove the constraint
			from body constraints set, which is not nice if we try
			to iterate over the set at the same time => make it empty */
    
    for (SET *item = SET_First (con); item; item = SET_Next (item))
      DOM_Remove_Constraint (dom, item->data);

    SET_Free (&dom->setmem, &con); /* free body's constraint set */
  }

  /* remove external constrains */
  {
    for (MAP *item = MAP_First (bod->conext); item; item = MAP_Next (item))
    {
      CONEXT *ext = item->data;
      ext->bod = NULL; /* see conext_remove_all to understand why */
    }

    MAP_Free (&dom->mapmem, &bod->conext);
  }

  /* delete from id based map */
  MAP_Delete (&dom->mapmem, &dom->children, (void*) (long) bod->id, NULL);
}

/* insert migrated body into the domain */
void insert_migrated_body (DOM *dom, BODY *bod)
{
  if (bod->label) /* map labeled bodies */
    MAP_Insert (&dom->mapmem, &dom->lab, bod->label, bod, (MAP_Compare)strcmp);

  MAP_Insert (&dom->mapmem, &dom->idb, (void*) (long) bod->id, bod, NULL);

  bod->dom = dom;
  dom->nbod ++;

  bod->next = dom->bod;
  if (dom->bod) dom->bod->prev = bod;
  dom->bod = bod;
}

/* remove migrated body from the domain */
static void remove_migrated_body (DOM *dom, BODY *bod)
{
  AABB_Detach_Body (dom->aabb, bod);

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

  /* remove external constrains */
  {
    for (MAP *item = MAP_First (bod->conext); item; item = MAP_Next (item))
    {
      CONEXT *ext = item->data;
      ext->bod = NULL; /* see conext_remove_all to understand why */
    }

    MAP_Free (&dom->mapmem, &bod->conext);
  }

  if (bod->label) /* delete labeled body */
    MAP_Delete (&dom->mapmem, &dom->lab, bod->label, (MAP_Compare)strcmp);

  /* delete from id based map */
  MAP_Delete (&dom->mapmem, &dom->idb, (void*) (long) bod->id, NULL);

  dom->nbod --;

  if (bod->prev) bod->prev->next = bod->next;
  else dom->bod = bod->next;
  if (bod->next) bod->next->prev = bod->prev;
}

/* balance bodies */
static int domain_balance (DOM *dom)
{
  int changes,
      num_gid_entries,
      num_lid_entries,
      num_import,
      *import_procs,
      num_export,
      *export_procs;

  ZOLTAN_ID_PTR import_global_ids,
		import_local_ids,
		export_global_ids,
		export_local_ids;

  COMOBJ *bodsend, *bodrecv, *consend, *conrecv, *ptr;

  int nbodsend, nbodrecv, nconsend, nconrecv;

  COMDATA *send, *recv, *qtr;

  int nsend, nrecv;

  MAP *export_bod = NULL,
      *export_con = NULL,
      *item;

  SET *adj, **del;

  MEM mem;

  unsigned int id;

  int i, *j;

  ERRMEM (send = calloc (dom->ncpu, sizeof (COMDATA)));

  /* before rebalancing - delete orhpans (children of removed bodies) */
  for (i = 0, del = dom->delch, qtr = send; i < dom->ncpu; i ++, del ++)
  {
    if (*del)
    {
      qtr->rank = i; /* send it there */
      qtr->ints = SET_Size (*del);
      ERRMEM (qtr->i = j = malloc (sizeof (int [qtr->ints])));
      for (SET *x = SET_First (*del); x; x = SET_Next (x), j ++) *j = (int) (long) x->data; /* ids of children to be deleted */
      SET_Free (&dom->setmem, del); /* empty this set to be reused next time */
      qtr ++;
    }
  }

  nsend = qtr - send;

  COM (MPI_COMM_WORLD, TAG_ORPHANS, send, nsend, &recv, &nrecv); /* send and receive orphans ids */

  for (i = 0, qtr = recv; i < nrecv; i ++, qtr ++)
  {
    for (j = qtr->i; j < qtr->i + qtr->ints; j ++) /* for each orphan id */
    {
      BODY *bod;

      ASSERT_DEBUG_EXT (bod = MAP_Find (dom->children, (void*) (long) (*j), NULL), "Invalid orphan id"); /* find child */
      remove_child (dom, bod);
      BODY_Destroy (bod);
    }
  }

  for (i = 0; i < nsend; i ++) free (send [i].i);
  free (send);
  free (recv);

  /* update body partitioning using mid points of their extents */
  ASSERT (Zoltan_LB_Balance (dom->zol, &changes, &num_gid_entries, &num_lid_entries,
	  &num_import, &import_global_ids, &import_local_ids, &import_procs,
	  &num_export, &export_global_ids, &export_local_ids, &export_procs) == ZOLTAN_OK, ERR_ZOLTAN);


  dom->nexpbod = num_export; /* record this figure for later output statistics */

  /* SUMMARY: After partitioning update some bodies will be exported to other partitions.
   *          We need to maintain user prescribed non-contact constraints attached to those
   *          bodies. For this reason 'bod->con' lists of exported bodies are first scanned
   *          and the 'export_con' set is created from all such constraints. At the same
   *          time 'export_bod' set is created, comprising all bodies to be exported; then
   *          all RIGLNK constraints are found and:
   *
   *       1. If a master body of RIGLNK constraint is exported, then the slave body is
   *          appended to the 'export_bod' set (so that both are exporet into the same location)
   *
   *       2. If only a slave body of RIGLNK constraint is exported, then this body
   *          is removed from the 'export_bod' set and the constraint is removed from 'export_con'
   *          set (neither the body nor the constraint is exported)
   *
   *          Constraint communication is done first, so that we have access to memory of needed
   *          constraints, before 'remove_migrated_body' would be invoked (which would delete those
   *          constrints as well). Then bodies from 'export_bod' are sent and some bodies
   *          are received. The sent bodies are removed from the domain while the reveived are
   *          insrted. Finally, the received constraints are inserted into the domain.
   *
   *          The last step is sending and receiving the 'child' bodies; these are copies
   *          of the currently stored 'parent' bodies sent into all the partitions which
   *          are overalpped by their bounding boxes.
   *          -----------------------------------------------------------------------  */

  for (i = 0; i < num_export; i ++) /* for each exported body */
  {
    BODY *bod;

    id = export_global_ids [i * num_gid_entries]; /* get id */

    ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) id, NULL), "Invalid body id"); /* identify body object */

    MAP_Insert (&dom->mapmem, &export_bod, bod, (void*) (long) export_procs [i], NULL); /* map this body to its export rank */

    for (adj = SET_First (bod->con); adj; adj = SET_Next (adj)) /* search adjacent constraints */
    {
      CON *con = adj->data;

      if (con->kind != CONTACT)
      {
	MAP_Insert (&dom->mapmem, &export_con, con,
	            (void*) (long) export_procs [i], NULL); /* map non-contact constraint to its export rank */

	con->id = 0; /* this way, constraint's id will not be freed in DOM_Remove_Constraint;
			non-contact constrints should have cluster-wide unique ids, as users
			could store some global variables in Python input and might later like
			to acces those constraints; this will not be implemented for contacts */
      }
    }
  }

  Zoltan_LB_Free_Data (&import_global_ids, &import_local_ids, &import_procs,
		       &export_global_ids, &export_local_ids, &export_procs); /* not needed */

  for (item = MAP_First (export_con); item; ) /* for each exported constraint */
  {
    CON *con = item->data;

    if (con->kind == RIGLNK) /* if a rigid link */
    {
      if (! MAP_Find (export_bod, con->master, NULL)) /* if the master body is not exported */
      {
        MAP_Delete (&dom->mapmem, &export_bod, con->slave, NULL); /* do not export the slave */
	item = MAP_Delete_Node (&dom->mapmem, &export_con, item); /* do not export the constraint */
	continue;
      }
      else if (! MAP_Find (export_bod, con->slave, NULL)) /* or if the slave body is not exported */
      {
        MAP_Insert (&dom->mapmem, &export_bod, con->slave, item->data, NULL); /* export slave */ 	
      }
    }

    item = MAP_Next (item);
  }

  /* communicate constraints */

  MEM_Init (&mem, sizeof (CONAUX), CONBLK);

  if ((nconsend = MAP_Size (export_con)))
  {
    ERRMEM (consend = malloc (sizeof (COMOBJ [nconsend])));

    for (item = MAP_First (export_con), ptr = consend; item; item = MAP_Next (item), ptr ++)
    {
      ptr->rank = (int) (long) item->data;
      ptr->o = item->key;
    }
  }
  else consend = NULL;

  COMOBJS (MPI_COMM_WORLD, TAG_CONAUX, (OBJ_Pack)pack_constraint, &mem, (OBJ_Unpack)unpack_constraint, consend, nconsend, &conrecv, &nconrecv);

  /* communicate parent bodies */

  if ((nbodsend = MAP_Size (export_bod)))
  {
    ERRMEM (bodsend = malloc (sizeof (COMOBJ [nbodsend])));

    for (item = MAP_First (export_bod), ptr = bodsend; item; item = MAP_Next (item), ptr ++)
    {
      ptr->rank = (int) (long) item->data;
      ptr->o = item->key;
    }
  }
  else bodsend = NULL;

  COMOBJS (MPI_COMM_WORLD, TAG_PARENTS, (OBJ_Pack)BODY_Parent_Pack, dom->owner, (OBJ_Unpack)BODY_Parent_Unpack, bodsend, nbodsend, &bodrecv, &nbodrecv);

  /* remove exported bodies */

  for (i = 0, ptr = bodsend; i < nbodsend; i ++, ptr ++)
  {
    remove_migrated_body (dom, ptr->o);
    BODY_Destroy (ptr->o);
  }

  /* insert imported bodies */

  for (i = 0, ptr = bodrecv; i < nbodrecv; i ++, ptr ++)
  {
    insert_migrated_body (dom, ptr->o);
  }

  /* insert imported constraints */

  dom->noid = 1; /* disable constraint ids generation */

  for (i = 0, ptr = conrecv; i < nconrecv; i ++, ptr ++)
  {
    CONAUX *aux;
    BODY *bod;
    CON *con;

    aux = ptr->o;

    ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) aux->master, NULL), "Invalid body id");

    switch (aux->kind)
    {
    case FIXPNT: con = DOM_Fix_Point (dom, bod, aux->vec [0]); break;
    case FIXDIR: con = DOM_Fix_Direction (dom, bod, aux->vec [0], aux->vec [1]); break;
    case VELODIR: con = DOM_Set_Velocity (dom, bod, aux->vec [0], aux->vec [1], aux->tms); break;
    case RIGLNK:
    {
      BODY *other;

      ASSERT_DEBUG_EXT (other = MAP_Find (dom->idb, (void*) (long) aux->slave, NULL), "Invalid body id");
      con = DOM_Put_Rigid_Link (dom, bod, other, aux->vec [0], aux->vec [1]);
    }
    break;
    }

    con->id = aux->id;
  }

  dom->noid = 0; /* enable constraint ids generation */

  MAP_Free (&dom->mapmem, &export_bod);
  MAP_Free (&dom->mapmem, &export_con);
  MEM_Release (&mem);
  free (bodsend);
  free (bodrecv);
  free (consend);
  free (conrecv);

  return changes;
}

/* try domein balancing */
static void domain_try_balance (DOM *dom)
{
  int sum, min, avg, max;

  /* compute inbalance of bodies in partitions */
  PUT_int_stats (1, &dom->nbod, &sum, &min, &avg, &max);

  double ratio = (double) max / (double) MAX (min, 1);

  if (ratio > dom->imbalance_tolerance) domain_balance (dom);
  else dom->nexpbod = 0;
}

/* return next pointer and realloc send data memory if needed */
inline static COMDATA* sendnextdata (int nsend, int *size, COMDATA **send)
{
  if (nsend >= *size)
  {
    (*size) *= 2;

    ERRMEM (*send = realloc (*send, sizeof (COMDATA [*size])));
  }

  return &(*send)[nsend];
}

/* return next pointer and realloc send object memory if needed */
inline static COMOBJ* sendnextobj (int nsend, int *size, COMOBJ **send)
{
  if (nsend >= *size)
  {
    (*size) *= 2;

    ERRMEM (*send = realloc (*send, sizeof (COMOBJ [*size])));
  }

  return &(*send)[nsend];
}

/* update children configurations and parent ranks */
static void domain_gossip (DOM *dom)
{
  COMOBJ *send, *recv, *ptr;
  int size, nsend, nrecv;
  BODY *bod;
  SET *item;

  size = MAX (dom->nbod, 128);
  ERRMEM (send = malloc (sizeof (COMOBJ [size])));

  /* gather are (body, rank) send paris */
  for (nsend = 0, ptr = send, bod = dom->bod; bod; bod = bod->next)
  {
    for (item = SET_First (bod->my.children); item; item = SET_Next (item))
    {
      ptr->rank = (int) (long) item->data;
      ptr->o = bod;
      ptr = sendnextobj (++ nsend, &size, &send);
    }
  }

  /* send configuration and velocity from parent to child bodies and
   * update children shapes while unpacking the sent data; also update parent
   * ranks as those might change when parents migrate while their children don't */
  COMOBJS (MPI_COMM_WORLD, TAG_CHILDREN_UPDATE, (OBJ_Pack)BODY_Child_Pack_State, dom,
           (OBJ_Unpack)BODY_Child_Unpack_State, send, nsend, &recv, &nrecv);

  free (send);
  free (recv);
}

/* remove all external constraints */
static void conext_remove_all (DOM *dom)
{
  CONEXT *ext;

  /* empty body external constraint sets */
  for (ext = dom->conext; ext; ext = ext->next)
  {
    if (ext->bod && ext->bod->conext) MAP_Free (&dom->mapmem, &ext->bod->conext);
  }

  /* erase in the domain */
  MEM_Release (&dom->extmem);
  dom->conext = NULL;
}

/* insert an external constraint */
static void conext_insert (DOM *dom, int rank, CONEXT *ext)
{
  /* set child rank */
  ext->rank = rank;

  /* append list */
  ext->next = dom->conext;
  dom->conext = ext;

  /* add to the body constraint adjacency */
  MAP_Insert (&dom->mapmem, &ext->bod->conext, (void*) (long) ext->id, ext, NULL);
}

/* create an external constraint */
static CONEXT* conext_create (DOM *dom, CON *con, BODY *bod)
{
  CONEXT *ext;

  ERRMEM (ext = MEM_Alloc (&dom->extmem));

  COPY (con->R, ext->R);
  NNCOPY (con->base, ext->base);
  ext->id = con->id;
  ext->bod = bod;

  if (bod == con->master)
  {
    ext->isma = 1;
    ext->sgp = con->mbox->sgp;
    COPY (con->mpnt, ext->point);
  }
  else
  {
    ext->isma = 0;
    ext->sgp = con->sbox->sgp;
    COPY (con->spnt, ext->point);
  }

  return ext;
}

/* pack an external constraint */
static void conext_pack (CONEXT *ext, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_doubles (dsize, d, doubles, ext->R, 3);
  pack_doubles (dsize, d, doubles, ext->point, 3);
  pack_doubles (dsize, d, doubles, ext->base, 9);
  pack_int (isize, i, ints, ext->id);
  pack_int (isize, i, ints, ext->bod->id);
  pack_int (isize, i, ints, ext->isma);
  pack_int (isize, i, ints, ext->sgp - ext->bod->sgp);
}

/* unpack an external constraint (do not set body pointer, but an id) */
CONEXT* conext_unpack (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  CONEXT *ext;

  ERRMEM (ext = MEM_Alloc (&dom->extmem));

  unpack_doubles (dpos, d, doubles, ext->R, 3);
  unpack_doubles (dpos, d, doubles, ext->point, 3);
  unpack_doubles (dpos, d, doubles, ext->base, 9);
  ext->id  = unpack_int (ipos, i, ints);
  ext->bod = (BODY*) (long) unpack_int (ipos, i, ints);
  ext->isma = unpack_int (ipos, i, ints);
  ext->sgp = (SGP*) (long) unpack_int (ipos, i, ints);

  return ext;
}

/* unpack an external constraint for a parent body */
CONEXT* conext_unpack_parent (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  CONEXT *ext = conext_unpack (dom, dpos, d, doubles, ipos, i, ints);

  ASSERT_DEBUG_EXT (ext->bod = MAP_Find (dom->idb, ext->bod, NULL), "Invalid body id");

  ext->sgp = &ext->bod->sgp [(int) (long) ext->sgp];

  return ext;
}

/* unpack an external constraint for a child body */
CONEXT* conext_unpack_child (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  CONEXT *ext = conext_unpack (dom, dpos, d, doubles, ipos, i, ints);

  ASSERT_DEBUG_EXT (ext->bod = MAP_Find (dom->children, ext->bod, NULL), "Invalid body id");

  ext->sgp = &ext->bod->sgp [(int) (long) ext->sgp];

  return ext;
}

/* migrate external constraints */
static void domain_glue_begin (DOM *dom)
{
  COMOBJ *send, *recv, *ptr;
  int k, size, nsend, nrecv;
  SET *item, *jtem, *del;
  CONEXT *ext;
  BODY *bod;
  CON *con;

  /* Remove all currrently stored external constraints */
  conext_remove_all (dom);

  size = MAX (dom->ncon, 128);

  ERRMEM (send = malloc (sizeof (COMOBJ [size])));

  /* Walk over all constraints and identify (con, rank) pairs to be exported */
  for (ptr = send, nsend = 0, con = dom->con; con; con = con->next)
  {
    BODY *bod [] = {con->master, con->slave};

    for (k = 0; k < 2; k ++)
    {
      if (bod [k])
      {
	if (bod [k]->flags & BODY_CHILD)
	{
	  ptr->rank = bod [k]->my.parent;
	  ptr->o = conext_create (dom, con, bod [k]);
	  ptr = sendnextobj (++ nsend, &size, &send);
	}
      }
    }
  }

  /* Send from children and receive in parents */
  COMOBJS (MPI_COMM_WORLD, TAG_CONEXT_TO_PARENTS, (OBJ_Pack)conext_pack, dom, (OBJ_Unpack)conext_unpack_parent, send, nsend, &recv, &nrecv);

  /* Insert received external constraints */
  for (k = 0, ptr = recv; k < nrecv; k ++, ptr ++) conext_insert (dom, ptr->rank, ptr->o);

  /* Clean up */
  for (k = 0, ptr = send; k < nsend; k ++, ptr ++) MEM_Free (&dom->extmem, ptr->o);
  free (recv);

  /* Now parent bodies contain complete sets of adjacent constraints;
   * send those sets back to children, ommiting their own constraints */
  for (del = NULL, ptr = send, nsend = 0, bod = dom->bod; bod; bod = bod->next)
  {
    for (item = SET_First (bod->my.children); item; item = SET_Next (item))
    {
      for (jtem = SET_First (bod->con); jtem; jtem = SET_Next (jtem)) /* local constraints */
      {
	ptr->rank = (int) (long) item->data; /* child rank */
	ptr->o = conext_create (dom, jtem->data, bod);
	SET_Insert (&dom->setmem, &del, ptr->o, NULL); /* to delete when done */
	ptr = sendnextobj (++ nsend, &size, &send);
      }

      for (MAP *jtem = MAP_First (bod->conext); jtem; jtem = MAP_Next (jtem)) /* external constraints */
      {
        ext = jtem->data;

	if (ext->rank != (int) (long) item->data) /* if not coming from this child */
	{
	  ptr->rank = (int) (long) item->data; /* child rank */
	  ptr->o = ext;
	  ptr = sendnextobj (++ nsend, &size, &send);
	}
      }
    }
  }

  /* Send from parents and receive in children */
  COMOBJS (MPI_COMM_WORLD, TAG_CONEXT_TO_CHILDREN, (OBJ_Pack)conext_pack, dom, (OBJ_Unpack)conext_unpack_child, send, nsend, &recv, &nrecv);

  /* Insert received new external constraints */
  for (k = 0, ptr = recv; k < nrecv; k ++, ptr ++) conext_insert (dom, -1, ptr->o); /* invalid -1 rank here */

  /* Clean up */
  for (item = SET_First (del); item; item = SET_Next (item)) MEM_Free (&dom->extmem, item->data);
  free (send);
  free (recv);
}

/* update external constraint reactions */
static void domain_glue_end (DOM *dom)
{
  COMDATA *send, *recv, *ptr;
  int k, *j, *l, size, nsend, nrecv;
  MEM memintpair;
  CONEXT *ext;
  double *R;
  BODY *bod;
  CON *con;

  size = MAX (dom->ncon, 128);

  ERRMEM (send = malloc (sizeof (COMDATA [size])));

  MEM_Init (&memintpair, sizeof (int [2]), size);

  /* Walk over all constraints and identify (con, rank) pairs to be exported */
  for (ptr = send, nsend = 0, con = dom->con; con; con = con->next)
  {
    BODY *bod [] = {con->master, con->slave};

    for (k = 0; k < 2; k ++)
    {
      if (bod [k])
      {
	if (bod [k]->flags & BODY_CHILD)
	{
	  ptr->rank = bod [k]->my.parent;
	  ptr->ints = 2;
	  ptr->doubles = 3;
	  ERRMEM (ptr->i = MEM_Alloc (&memintpair));
	  ptr->i[0] = bod[k]->id;
	  ptr->i[1] = con->id;
	  ptr->d = con->R;
	  ptr = sendnextdata (++ nsend, &size, &send);
	}
      }
    }
  }

  /* Send updated reactions of external constraints from children to parents */
  COM (MPI_COMM_WORLD, TAG_CONEXT_UPDATE_PARENTS, send, nsend, &recv, &nrecv);

  /* Insert received external constraints */
  for (k = 0, ptr = recv; k < nrecv; k ++, ptr ++)
  {
    for (j = ptr->i, l = j + ptr->ints, R = ptr->d; j < l; j += 2, R += 3)
    {
      ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) j [0], NULL), "Invalid body id");
      ASSERT_DEBUG_EXT (ext = MAP_Find (bod->conext, (void*) (long) j [1], NULL), "Invalid external constraint id");
      COPY (R, ext->R);
    }
  }

  /* Clean up */
  MEM_Release (&memintpair);
  free (send);
  free (recv);
}

/* create MPI related data */
static void create_mpi (DOM *dom)
{
  MPI_Comm_rank (MPI_COMM_WORLD, &dom->rank); /* store rank */

  MPI_Comm_size (MPI_COMM_WORLD, &dom->ncpu); /* store size */

  dom->cid = (dom->rank + 1); /* overwrite */

  dom->sparecid = NULL; /* initialize spare ids sets */
  dom->sparebid = NULL;

  dom->noid = 0; /* assign constraint ids in 'insert' routine (turned off when importing non-contacts) */

  dom->children = NULL; /* initially empty */

  ERRMEM (dom->delch = calloc (dom->ncpu, sizeof (SET*)));

  dom->conext = NULL;

  MEM_Init (&dom->extmem, sizeof (CONEXT), CONBLK);

  ASSERT (dom->zol = Zoltan_Create (MPI_COMM_WORLD), ERR_ZOLTAN); /* zoltan context for body partitioning */

  dom->imbalance_tolerance = 1.3;

  /* general parameters */
  Zoltan_Set_Param (dom->zol, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param (dom->zol, "DEBUG_MEMORY", "0");
  Zoltan_Set_Param (dom->zol, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param (dom->zol, "NUM_LID_ENTRIES", "0");
  Zoltan_Set_Param (dom->zol, "OBJ_WEIGHT_DIM", "1"); /* bod->ndof */
 
  /* load balancing parameters */
  Zoltan_Set_Param (dom->zol, "LB_METHOD", "RCB");
  Zoltan_Set_Param (dom->zol, "IMBALANCE_TOL", "1.2");
  Zoltan_Set_Param (dom->zol, "AUTO_MIGRATE", "FALSE"); /* we shall use COMOBJS */
  Zoltan_Set_Param (dom->zol, "RETURN_LISTS", "EXPORT"); /* the rest will be done by COMOBJS */

  /* RCB parameters */
  Zoltan_Set_Param (dom->zol, "RCB_OVERALLOC", "1.3");
  Zoltan_Set_Param (dom->zol, "RCB_REUSE", "1");
  Zoltan_Set_Param (dom->zol, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param (dom->zol, "CHECK_GEOM", "1");
  Zoltan_Set_Param (dom->zol, "KEEP_CUTS", "1");

  /* callbacks */
  Zoltan_Set_Fn (dom->zol, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) body_count, dom);
  Zoltan_Set_Fn (dom->zol, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)()) body_list, dom);
  Zoltan_Set_Fn (dom->zol, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) dimensions, dom);
  Zoltan_Set_Fn (dom->zol, ZOLTAN_GEOM_MULTI_FN_TYPE, (void (*)()) midpoints, dom);
}

/* destroy MPI related data */
static void destroy_mpi (DOM *dom)
{
  free (dom->delch);

  MEM_Release (&dom->extmem);

  Zoltan_Destroy (&dom->zol);
}
#endif

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
  aabb->dom = dom;
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

#if MPI
  create_mpi (dom);
#endif

  return dom;
}

/* insert a body into the domain */
void DOM_Insert_Body (DOM *dom, BODY *bod)
{
  if (bod->label) /* map labeled bodies */
    MAP_Insert (&dom->mapmem, &dom->lab, bod->label, bod, (MAP_Compare)strcmp);

#if MPI
  ASSERT_DEBUG (dom->rank == 0, "Bodies can only be inserted at rank 0");

  if (SET_Size (dom->sparebid))
  {
    SET *item;
   
    item = SET_First (dom->sparebid);
    bod->id = (unsigned int) (long) item->data; /* use a previously freed id */
    SET_Delete (&dom->setmem, &dom->sparebid, item->data, NULL);
  }
  else
#endif
  bod->id = dom->bid ++; /* due to the above assertion this is fine in parallel */
  MAP_Insert (&dom->mapmem, &dom->idb, (void*) (long) bod->id, bod, NULL);

  bod->dom = dom;
  dom->nbod ++;

  bod->next = dom->bod;
  if (dom->bod) dom->bod->prev = bod;
  dom->bod = bod;

  if (dom->time > 0) SET_Insert (&dom->setmem, &dom->newb, bod, NULL);

  AABB_Insert_Body (dom->aabb, bod);
}

/* remove a body from the domain */
void DOM_Remove_Body (DOM *dom, BODY *bod)
{
  AABB_Delete_Body (dom->aabb, bod);

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

#if MPI
  /* remove external constrains */
  {
    for (MAP *item = MAP_First (bod->conext); item; item = MAP_Next (item))
    {
      CONEXT *ext = item->data;
      ext->bod = NULL; /* see conext_remove_all to understand why */
    }

    MAP_Free (&dom->mapmem, &bod->conext);
  }
#endif

  if (bod->label) /* delete labeled body */
    MAP_Delete (&dom->mapmem, &dom->lab, bod->label, (MAP_Compare)strcmp);

  /* delete from id based map */
  MAP_Delete (&dom->mapmem, &dom->idb, (void*) (long) bod->id, NULL);

  dom->nbod --;

  if (bod->prev) bod->prev->next = bod->next;
  else dom->bod = bod->next;
  if (bod->next) bod->next->prev = bod->prev;


  if (dom->time > 0) SET_Insert (&dom->setmem, &dom->delb, (void*) (long) bod->id, NULL);

#if MPI
  /* free body id */
  SET_Insert (&dom->setmem, &dom->sparebid, (void*) (long) bod->id, NULL);

  /* gather children ids to be deleted during balancing */
  for (SET *item = SET_First (bod->my.children); item; item = SET_Next (item))
    SET_Insert (&dom->setmem, &dom->delch [(int) (long) item->data], (void*) (long) bod->id, NULL);
#endif
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
  SGP *sgp;
  int n;

  if ((n = SHAPE_Sgp (bod->sgp, bod->nsgp, pnt)) < 0) return NULL;

  sgp = &bod->sgp [n];
  con = insert (dom, bod, NULL, 1);
  con->kind = FIXPNT;
  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  IDENTITY (con->base);
  con->mgobj = sgp->gobj;
  con->mshp = sgp->shp;
  con->mbox = sgp->box;

  return con;
}

/* fix a referential point of the body along the spatial direction */
CON* DOM_Fix_Direction (DOM *dom, BODY *bod, double *pnt, double *dir)
{
  CON *con;
  SGP *sgp;
  int n;

  if ((n = SHAPE_Sgp (bod->sgp, bod->nsgp, pnt)) < 0) return NULL;

  sgp = &bod->sgp [n];
  con = insert (dom, bod, NULL, 1);
  con->kind = FIXDIR;
  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  localbase (dir, con->base);
  con->mgobj = sgp->gobj;
  con->mshp = sgp->shp;
  con->mbox = sgp->box;

  return con;
}

/* prescribe a velocity of the referential point along the spatial direction */
CON* DOM_Set_Velocity (DOM *dom, BODY *bod, double *pnt, double *dir, TMS *vel)
{
  CON *con;
  SGP *sgp;
  int n;

  if ((n = SHAPE_Sgp (bod->sgp, bod->nsgp, pnt)) < 0) return NULL;

  sgp = &bod->sgp [n];
  con = insert (dom, bod, NULL, 1);
  con->kind = VELODIR;
  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  localbase (dir, con->base);
  con->tms = vel;
  con->mgobj = sgp->gobj;
  con->mshp = sgp->shp;
  con->mbox = sgp->box;

  return con;
}

/* insert rigid link constraint between two (referential) points of bodies; if one of the body
 * pointers is NULL then the link acts between the other body and the fixed (spatial) point */
CON* DOM_Put_Rigid_Link (DOM *dom, BODY *master, BODY *slave, double *mpnt, double *spnt)
{
  double v [3], d;
  CON *con;
  SGP *msgp, *ssgp;
  int m, s;

  if (!master)
  {
    master = slave;
    mpnt = spnt;
    slave = NULL;
  }

  ASSERT_DEBUG (master, "At least one body pointer must not be NULL");

  if (master && (m = SHAPE_Sgp (master->sgp, master->nsgp, mpnt)) < 0) return NULL;

  if (slave && (s = SHAPE_Sgp (slave->sgp, slave->nsgp, spnt)) < 0) return NULL;

  msgp = &master->sgp [m];
  if (slave) ssgp = &slave->sgp [s];
  else ssgp = NULL;

  SUB (mpnt, spnt, v);
  d = LEN (v);
  
  if (d < GEOMETRIC_EPSILON) /* no point in keeping very short links */
  {
    con = insert (dom, master, slave, 1);
    con->kind = FIXPNT;
    COPY (mpnt, con->point);
    COPY (mpnt, con->mpnt);
    COPY (spnt, con->spnt);
    IDENTITY (con->base);
    con->mgobj = msgp->gobj;
    con->mshp = msgp->shp;
    con->mbox = msgp->box;

    if (slave)
    {
      con->sgobj = ssgp->gobj;
      con->sshp = ssgp->shp;
      con->sbox = ssgp->box;

      AABB_Exclude_Gobj_Pair (dom->aabb, master->id, m, slave->id, s); /* no contact between this pair */
    }
  }
  else
  {
    con = insert (dom, master, slave, 1);
    con->kind = RIGLNK;
    COPY (mpnt, con->point);
    COPY (mpnt, con->mpnt);
    COPY (spnt, con->spnt);
    RIGLNK_LEN (con->Z) = d; /* initial distance */
    con->mgobj = msgp->gobj;
    con->mshp = msgp->shp;
    con->mbox = msgp->box;

    if (slave)
    { 
      con->sgobj = ssgp->gobj;
      con->sshp = ssgp->shp;
      con->sbox = ssgp->box;
    }

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
  MAP_Delete (&dom->mapmem, &dom->idc, (void*) (long) con->id, NULL);

  /* remove from list */
  if (con->prev)
    con->prev->next = con->next;
  else dom->con = con->next;
  if (con->next)
    con->next->prev = con->prev;
  dom->ncon --;

  /* remove from local dynamics */
  if (con->dia) LOCDYN_Remove (dom->ldy, con->dia);

#if MPI
  /* free constraint id if needed */
  if (con->id) SET_Insert (&dom->setmem, &dom->sparecid, (void*) (long) con->id, NULL);
#endif

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

#if MPI
  if (dom->rank == 0)
#endif
  if (dom->verbose) printf ("DOMAIN ... "), fflush (stdout);

#if MPI
  SOLFEC_Timer_Start (dom->owner, "TIMBAL");

  domain_try_balance (dom);

  SOLFEC_Timer_End (dom->owner, "TIMBAL");
#endif

  SOLFEC_Timer_Start (dom->owner, "TIMINT");

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

  SOLFEC_Timer_End (dom->owner, "TIMINT");

  /* detect contacts */
  timerstart (&timing);

  AABB_Update (dom->aabb, aabb_algorithm (dom),
    dom, (BOX_Overlap_Create) overlap_create,
    (BOX_Overlap_Release) overlap_release);

  aabb_timing (dom, timerend (&timing));

  SOLFEC_Timer_Start (dom->owner, "CONDET");

  /* sparsify new contacts */
  sparsify_contacts (dom);

  SOLFEC_Timer_End (dom->owner, "CONDET");

  SOLFEC_Timer_Start (dom->owner, "TIMINT");

  /* update old constraints */
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

  SOLFEC_Timer_End (dom->owner, "TIMINT");

#if MPI
  SOLFEC_Timer_Start (dom->owner, "TIMBAL");

  domain_glue_begin (dom);

  SOLFEC_Timer_End (dom->owner, "TIMBAL");
#endif

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

#if MPI
  SOLFEC_Timer_Start (dom->owner, "TIMBAL");

  domain_glue_end (dom);

  SOLFEC_Timer_End (dom->owner, "TIMBAL");
#endif

  SOLFEC_Timer_Start (dom->owner, "TIMINT");

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

  SOLFEC_Timer_End (dom->owner, "TIMINT");
}

#if MPI
/* balance children according to a given geometric partitioning */
void DOM_Balance_Children (DOM *dom, struct Zoltan_Struct *zol)
{
  COMOBJ *bodsend, *bodrecv, *ptr;
  int *procs, numprocs, i, id;
  int nbodsend, nbodrecv;
  COMDATA *send, *recv;
  int nsend, nrecv;
  MEM mem;

  struct pair
  {
    union { int id; BODY *ptr; } bod;
    int rank;
  };

  ERRMEM (procs = malloc (sizeof (int [dom->ncpu])));

  MEM_Init (&mem, sizeof (struct pair), CONBLK);

  SET *procset,
      *delset,
      *sndset,
      *item;

  delset = sndset = NULL;

  /* 1. create an empty set 'delset' of pairs (id, rank) of children
   *    that will need to be removed from their partitions
   * 2. create an empty set 'sndset' of pairs (body, rank) of bodies
   *    to be sent as children to other partitions */

  procset = NULL;

  for (BODY *bod = dom->bod; bod; bod = bod->next)
  {
    double *e = bod->extents;

    Zoltan_LB_Box_Assign (zol, e[0], e[1], e[2], e[3], e[4], e[5], procs, &numprocs);

    /* 3. create a set copy 'procset' of 'procs' table */
    for (i = 0; i < numprocs; i ++) 
    {
      if (dom->rank != procs [i]) /* skip current rank */
	SET_Insert (&dom->setmem, &procset, (void*) (long) procs [i], NULL);
    }

    /* 3.5. if this body just migrated here, a child copy
     *      of it could still be here; schedule it for deletion */
    if (SET_Contains (bod->my.children, (void*) (long) dom->rank, NULL))
    {
      struct pair *p;

      ERRMEM (p = MEM_Alloc (&mem));
      p->bod.id = bod->id;
      p->rank = dom->rank;
      SET_Insert (&dom->setmem, &delset, p, NULL); /* schedule for deletion with other children */
      SET_Delete (NULL, &bod->my.children, (void*) (long) dom->rank, NULL); /* delete from bodies children set */
    }

     /* 4. for each x in bod->my.children, if not in procset,
      *    add (bod->id, x) to delset, remove x from bod->my.children */
    for (SET *x = SET_First (bod->my.children); x; )
    {
      if (!SET_Contains (procset, x->data, NULL))
      {
	struct pair *p;

	ERRMEM (p = MEM_Alloc (&mem));
	p->bod.id = bod->id;
	p->rank = (int) (long) x->data;
	SET_Insert (&dom->setmem, &delset, p, NULL);
	x = SET_Delete_Node (NULL, &bod->my.children, x);
      }
      else x = SET_Next (x);
    }

     /* 5. for each x in procset, if not in bod->my.children,
      *    add (bod, x) to sndset, insert x into bod->my.children */
    for (SET *x = SET_First (procset); x; x = SET_Next (x))
    {
      if (!SET_Contains (bod->my.children, x->data, NULL))
      {
	struct pair *p;

	ERRMEM (p = MEM_Alloc (&mem));
	p->bod.ptr = bod;
	p->rank = (int) (long) x->data;
	SET_Insert (&dom->setmem, &sndset, p, NULL);
	SET_Insert (NULL, &bod->my.children, x->data, NULL);
      }
    }

    /* reset processors set */
    SET_Free (&dom->setmem, &procset);
  }

  /* 6. communicate delset and delete unwanted children */

  ERRMEM (send = calloc (dom->ncpu, sizeof (COMDATA))); /* one for each processor */

  /* pack ids into specific rank data */
  for (SET *x = SET_First (delset); x; x = SET_Next (x))
  {
    struct pair *p = x->data;
    COMDATA *c = &send [p->rank];

    c->ints ++;
    ERRMEM (c->i = realloc (c->i, sizeof (int [c->ints])));
    c->i [c->ints - 1] = p->bod.id;
  }

  /* compress the send buffer */
  for (nsend = i = 0; i < dom->ncpu; i ++)
  {
    if (send [i].ints)
    {
      send [nsend] = send [i];
      send [nsend].rank = i;
      nsend ++;
    }
  }

  /* send and receive */
  COM (MPI_COMM_WORLD, TAG_CHILDREN_DELETE, send, nsend, &recv, &nrecv);

  /* delete children */
  for (nrecv --; nrecv >= 0; nrecv --)
  {
    for (i = 0; i < recv[nrecv].ints; i ++)
    {
      BODY *bod;

      id = recv[nrecv].i [i];

      ASSERT_DEBUG_EXT (bod = MAP_Find (dom->children, (void*) (long) id, NULL), "Invalid body id");
      remove_child (dom, bod);
      BODY_Destroy (bod);
    }
  }

  /* 7. communicate sndset and insert new children */

  if ((nbodsend = SET_Size (sndset)))
  {
    ERRMEM (bodsend = malloc (sizeof (COMOBJ [nbodsend])));

    /* set up send buffer */
    for (item = SET_First (sndset), ptr = bodsend; item; item = SET_Next (item), ptr ++)
    {
      struct pair *p = item->data;

      ptr->rank = p->rank;
      ptr->o = p->bod.ptr;
    }
  }
  else bodsend = NULL;

  dom->nexpchild = nbodsend; /* record for later statistics */

  /* send and receive children */
  COMOBJS (MPI_COMM_WORLD, TAG_CHILDREN_INSERT, (OBJ_Pack)BODY_Child_Pack, dom->owner, (OBJ_Unpack)BODY_Child_Unpack, bodsend, nbodsend, &bodrecv, &nbodrecv);

  /* insert children */
  for (i = 0, ptr = bodrecv; i < nbodrecv; i ++, ptr ++)
  {
    insert_child (dom, ptr->o);
  }

  /* clean up */
  for (i = 0; i < nsend; i ++) free (send [i].i);
  SET_Free (&dom->setmem, &delset);
  SET_Free (&dom->setmem, &sndset);
  MEM_Release (&mem);
  free (bodsend);
  free (bodrecv);
  free (procs);
  free (send);
  free (recv);
}

/* update children shapes */
void DOM_Update_Children (DOM *dom)
{
  /* update children configurations and parent ranks */
  domain_gossip (dom);
}
#endif

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
    id = (int) (long) item->data;
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

    if (bod->label) PBF_Label (bf, bod->label); /* label body record for fast access */

    BODY_Write_State (bod, bf);
  }
}

/* read domain state */
void DOM_Read_State (DOM *dom, PBF *bf)
{
  int ncon;

  /* clear contacts */
  MAP_Free (&dom->mapmem, &dom->idc);
  MEM_Release (&dom->conmem);
  dom->con = NULL;
  dom->ncon = 0;

  for (; bf; bf = bf->next)
  {
    if (PBF_Label (bf, "DOM"))
    {
      /* read time step */

      ASSERT (PBF_Label (bf, "STEP"), ERR_FILE_FORMAT);

      PBF_Double (bf, &dom->step, 1);

      /* read contacts */

      ASSERT (PBF_Label (bf, "CONS"), ERR_FILE_FORMAT);
    
      PBF_Int (bf, &ncon, 1);

      for (int n = 0; n < ncon; n ++)
      {
	CON *con;
	
	con = read_constraint (dom, bf);
	MAP_Insert (&dom->mapmem, &dom->idc, (void*) (long) con->id, con, NULL);
	con->next = dom->con;
	if (dom->con) dom->con->prev = con;
	dom->con = con;
      }

      dom->ncon += ncon;

      /* read ids of bodies that need to be deleted and remove them from all containers */

      unsigned int id;
      int size, n;

      ASSERT (PBF_Label (bf, "DELBODS"), ERR_FILE_FORMAT);

      PBF_Int (bf, &size, 1);

      for (n = 0; n < size; n ++)
      {
	BODY *bod;

	PBF_Uint (bf, &id, 1);

	ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) id, NULL), "Invalid body id");

	if (bod->label) MAP_Delete (&dom->mapmem, &dom->lab, bod->label, (MAP_Compare) strcmp);
	MAP_Delete (&dom->mapmem, &dom->idb, (void*) (long) id, NULL);
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
	MAP_Insert (&dom->mapmem, &dom->idb, (void*) (long) bod->id, bod, NULL);
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
	ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) id, NULL), "Body id invalid");
	BODY_Read_State (bod, bf);
      }
    }
  }
}

/* read state of an individual body */
int  DOM_Read_Body (DOM *dom, PBF *bf, BODY *bod)
{
  if (bod->label)
  {
    for (; bf; bf = bf->next)
    {
      if (PBF_Label (bf, bod->label))
      {
	BODY_Read_State (bod, bf);
	return 1;
      }
    }
  }
  else
  {
    for (; bf; bf = bf->next)
    {
      if (PBF_Label (bf, "BODS"))
      {
	int nbod;

	PBF_Int (bf, &nbod, 1);

	for (int n = 0; n < nbod; n ++)
	{
	  unsigned int id;
	  BODY *obj;

	  PBF_Uint (bf, &id, 1);
	  ASSERT_DEBUG_EXT (obj = MAP_Find (dom->idb, (void*) (long) id, NULL), "Body id invalid");
	  if (bod->id == obj->id) 
	  {
	    BODY_Read_State (bod, bf);
	    return 1;
	  }
	  else /* skip body and continue */
	  {
	    BODY fake;

	    ERRMEM (fake.conf = malloc (sizeof (double [BODY_Conf_Size (obj)])));
	    ERRMEM (fake.velo = malloc (sizeof (double [obj->dofs])));
	    fake.shape = NULL;

	    BODY_Read_State (&fake, bf);

	    free (fake.conf);
	    free (fake.velo);
	  }
	}
      }
    }
  }

  return 0;
}

/* read state of an individual constraint */
int  DOM_Read_Constraint (DOM *dom, PBF *bf, CON *con)
{
  for (; bf; bf = bf->next)
  {
    int ncon;

    if (PBF_Label (bf, "CONS"))
    {
      PBF_Int (bf, &ncon, 1);

      for (int n = 0; n < ncon; n ++)
      {
	CON *obj = read_constraint (dom, bf);

	if (con->id == obj->id)
	{
	  *con = *obj;
          MEM_Free (&dom->conmem, obj); /* not needed */
	  return 1;
	}
	else MEM_Free (&dom->conmem, obj); /* skip and continue */
      }
    }
  }

  return 0;
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

#if MPI
  destroy_mpi (dom);
#endif

  free (dom);
}
