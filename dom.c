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
#include "set.h"
#include "dom.h"
#include "goc.h"
#include "tmr.h"
#include "pck.h"
#include "err.h"
#include "dio.h"
#include "cra.h"

#if MPI
#include "put.h"
#include "com.h"
#endif

#define CONBLK 128 /* constraints memory block size */
#define MAPBLK 128 /* map items memory block size */
#define SETBLK 128 /* set items memory block size */

/* excluded surface pairs comparison */
static int pair_compare (int *a, int *b)
{
  if (a[0] < b[0]) return -1;
  else if (a[0] == b[0] && a[1] < b[1]) return -1;
  else if (a[0] == b[0] && a[1] == b[1]) return 0;
  else return 1;
}

/* create aabb data */
static AABB_DATA* aabb_create_data (void)
{
  AABB_DATA *data;
  double part;

  ERRMEM (data = malloc (sizeof (AABB_DATA)));

  part = 1.0 / (double) BOXALG_COUNT;
  data->aabb_limits [0] = 0.0;
  data->aabb_counter = 0;
  data->aabb_algo = 0;

  return data;
}

/* free aabb data */
static void aabb_destroy_data (AABB_DATA *data)
{
  free (data);
}

/* fastest box overlap algorithm */
static BOXALG aabb_algorithm (DOM *dom)
{
#if 0
  AABB_DATA *data = dom->aabb_data;
  double num, *tim, *lim;
  int i;

  if (data->aabb_counter < BOXALG_COUNT)
  {
    data->aabb_algo = data->aabb_counter ++; /* at first test all algorithms */
  }
  else /* when all have been tested */
  {
    tim = data->aabb_timings;
    lim = data->aabb_limits;

    for (i = 0; i < BOXALG_COUNT; i ++) lim [i+1] = lim [i] + tim [i]; /* sum up */
    for (i = 1; i <= BOXALG_COUNT; i ++) lim [i] /= lim [BOXALG_COUNT]; /* normalize */

    num = DRAND(); /* random in [0, 1] */

    for (i = 0; i < BOXALG_COUNT; i ++)
    {
      if (num >= lim [i] && num < lim [i+1])
      {
	data->aabb_algo = i;
	return i;
      }
    }
  }

  return data->aabb_algo;
#else
  return HYBRID; /* TODO: remove when the above proves robust */
#endif
}

/* update aabb timing related data */
static void aabb_timing (DOM *dom, double timing)
{
  AABB_DATA *data = dom->aabb_data;

  data->aabb_timings [data->aabb_algo] = timing;
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

/* compare two constraints */
#define CONCMP ((SET_Compare) constraint_compare)
static int constraint_compare (CON *one, CON *two)
{
  short oc = (one->kind == CONTACT), tc = (two->kind == CONTACT);
  BODY *onebod [2], *twobod [2];
  SGP *onesgp [2], *twosgp [2];

  if (oc > tc) return -1; /* contacts are smaller than non-contacts */
  else if (oc < tc) return 1; /* -||- */
  else if ((oc + tc) == 0) return (one->id < two->id ? -1 : (one->id == two->id ? 0 : 1)); /* non-contacts: compare IDs */

  if (one->master < one->slave) /* contacts: order pointers before comparing */
  {
    onebod [0] = one->master;
    onesgp [0] = one->msgp;
    onebod [1] = one->slave;
    onesgp [1] = one->ssgp;
  }
  else if (one->master > one->slave)
  {
    onebod [0] = one->slave;
    onesgp [0] = one->ssgp;
    onebod [1] = one->master;
    onesgp [1] = one->msgp;
  }
  else
  {
    onebod [0] = onebod [1] = one->master;

    if (one->msgp < one->ssgp)
    {
      onesgp [0] = one->msgp;
      onesgp [1] = one->ssgp;
    }
    else
    {
      onesgp [0] = one->ssgp;
      onesgp [1] = one->msgp;
    }
  }

  if (two->master < two->slave)
  {
    twobod [0] = two->master;
    twosgp [0] = two->msgp;
    twobod [1] = two->slave;
    twosgp [1] = two->ssgp;
  }
  else if (two->master > two->slave)
  {
    twobod [0] = two->slave;
    twosgp [0] = two->ssgp;
    twobod [1] = two->master;
    twosgp [1] = two->msgp;
  }
  else
  {
    twobod [0] = twobod [1] = two->master;

    if (two->msgp < two->ssgp)
    {
      twosgp [0] = two->msgp;
      twosgp [1] = two->ssgp;
    }
    else
    {
      twosgp [0] = two->ssgp;
      twosgp [1] = two->msgp;
    }
  }

  if (onebod [0] < twobod [0]) return -1; /* compare lexicographically by body pointers and sgps */
  else if (onebod [0] == twobod [0])
  {
    if (onebod [1] < twobod [1]) return -1;
    else if (onebod [1] == twobod [1])
    {
      if (onesgp [0] < twosgp [0]) return -1;
      else if (onesgp [0] == twosgp [0])
      {
	if (onesgp [1] < twosgp [1]) return -1;
	else if (onesgp [1] == twosgp [1]) return 0;
      }
    }
  }

  return 1;
}

/* insert a new constrait between two bodies */
static CON* insert (DOM *dom, BODY *master, BODY *slave, SGP *msgp, SGP *ssgp, short kind)
{
  CON *con;

  /* ensure that master pointers were passed;
   * tollerate NULL slave pointers, as this indicates a single body constraint */
  ASSERT_DEBUG (master && msgp, "At least master body pointers must be passed");

  ERRMEM (con = MEM_Alloc (&dom->conmem));
  con->kind = kind;
  con->master = master;
  con->slave = slave;
  con->msgp = msgp;
  con->ssgp = ssgp;

  /* insert into list */
  con->next = dom->con;
  if (dom->con) dom->con->prev = con;
  dom->con = con;
  dom->ncon ++;

#if MPI
  if (dom->noid == 0) /* if id generation is enabled */
  {
#endif
    SET *item = SET_First (dom->sparecid);

    if (item)
    {
      con->id = (unsigned int) (long) item->data; /* use a previously freed id */
      SET_Delete (&dom->setmem, &dom->sparecid, item->data, NULL);
    }
    else
    {
      con->id = dom->cid;

#if MPI
      ASSERT (((unsigned long long) dom->cid) + ((unsigned long long) dom->ncpu) < UINT_MAX, ERR_DOM_TOO_MANY_CONSTRAINTS);
      dom->cid += dom->ncpu; /* every ncpu number */
#else
      ASSERT (((unsigned long long) dom->cid) + ((unsigned long long) 1) < UINT_MAX, ERR_DOM_TOO_MANY_CONSTRAINTS);
      con->id = dom->cid ++;
#endif
    }
#if MPI
  }
  else con->id = dom->noid; /* assign the 'noid' as it was imported with a constraint */
#endif

#if PARDEBUG
  ASSERT_DEBUG (!SET_Contains (master->con, con, CONCMP), "Constraint %s with id %d already in body list", CON_Kind (con), con->id);
  ASSERT_DEBUG (!slave || (slave && !SET_Contains (slave->con, con, CONCMP)), "Constraint %s with id %d already in body list", CON_Kind (con), con->id);
#endif

  /* insert into body constraint adjacency */
  SET_Insert (&dom->setmem, &master->con, con, CONCMP);
  if (slave) SET_Insert (&dom->setmem, &slave->con, con, CONCMP);

#if PARDEBUG
  ASSERT_DEBUG (SET_Contains (master->con, con, CONCMP), "Failed to insert constraint %s with id %d into body list", CON_Kind (con), con->id);
  ASSERT_DEBUG (!slave || (slave && SET_Contains (slave->con, con, CONCMP)), "Failed to insert constraint %s with id %d already into body list", CON_Kind (con), con->id);
#endif

  /* map by id */
  MAP_Insert (&dom->mapmem, &dom->idc, (void*) (long) con->id, con, NULL);

  return con;
}

/* insert a contact into the constraints set */
static CON* insert_contact (DOM *dom, BODY *master, BODY *slave, SGP *msgp, SGP *ssgp, double *mpntspa,
       double *spntspa, double *normal, double area , double gap, SURFACE_MATERIAL *mat, short paircode)
{
  CON *con;

  con = insert (dom, master, slave, msgp, ssgp, CONTACT); /* do not insert into LOCDYN yet, only after sparsification */
  COPY (mpntspa, con->point);
  BODY_Ref_Point (master, msgp, mpntspa, con->mpnt); /* referential image */
  BODY_Ref_Point (slave, ssgp, spntspa, con->spnt);
  localbase (normal, con->base);
  con->area = area;
  con->gap = gap;
  con->paircode = paircode;
  con->state |= SURFACE_MATERIAL_Transfer (dom->time, mat, &con->mat); /* transfer surface pair data from the database to the local variable */
  con->state |= CON_NEW;  /* mark as newly created */
  return con;
}

/* does a potential contact already exists ? */
static int contact_exists (BOX *one, BOX *two)
{
  CON aux;

  aux.kind = CONTACT;
  aux.master = one->body;
  aux.msgp = one->sgp;
  aux.slave = two->body;
  aux.ssgp = two->sgp;

  return SET_Contains (one->body->con, &aux, CONCMP);
}

/* box overlap creation callback */
static void overlap_create (DOM *dom, BOX *one, BOX *two)
{
  double onepnt [3], twopnt [3], normal [3], gap, area;
  int state, spair [2], pair [2];
  SURFACE_MATERIAL *mat;
  short paircode;
  CON *con;

  if (contact_exists (one, two)) return;

  state = gobjcontact (
    CONTACT_DETECT, GOBJ_Pair_Code (one, two),
    one->sgp->shp, one->sgp->gobj,
    two->sgp->shp, two->sgp->gobj,
    onepnt, twopnt, normal, &gap, &area, spair);

  if (state)
  {
    ASSERT_DEBUG (gap <= 0, "A contact with positive gap (%g) was detected which indicates a bug in goc.c", gap);

    if (gap <= dom->depth) dom->flags |= DOM_DEPTH_VIOLATED;

    /* set surface pair data if there was a contact */
    mat = SPSET_Find (dom->sps, spair [0], spair [1]);

    if (dom->excluded)
    {
      if (spair [0] <= spair [1]) { pair [0] = spair [0]; pair [1] = spair [1]; }
      else { pair [0] = spair [1]; pair [1] = spair [0]; }

      if (SET_Contains (dom->excluded, pair, (SET_Compare) pair_compare)) return; /* exluded pair */
    }
  }

  switch (state)
  {
    case 1: /* first body has outward normal => second body is the master */
    {
      paircode = GOBJ_Pair_Code (one, two);
      con = insert_contact (dom, two->body, one->body, two->sgp, one->sgp, twopnt, onepnt, normal, area, gap, mat, paircode);
      con->spair [0] = spair [1];
      con->spair [1] = spair [0];
    }
    break;
    case 2:  /* second body has outward normal => first body is the master */
    {
      paircode = GOBJ_Pair_Code (two, one);
      con = insert_contact (dom, one->body, two->body, one->sgp, two->sgp, onepnt, twopnt, normal, area, gap, mat, paircode);
      con->spair [0] = spair [0];
      con->spair [1] = spair [1];
    }
    break;
  }
}

#if MPI
/* schedule remote deletion of external constraints */
static void ext_to_remove (DOM *dom, CON *con)
{
  for (SET *item = SET_First (con->ext); item; item = SET_Next (item))
  {
    SET_Insert (&dom->setmem, &dom->dbd [(int) (long) item->data].remove, (void*) (long) con->id, NULL);
  }
}
#endif

/* update contact data */
static void update_contact (DOM *dom, CON *con)
{
  double mpnt [3], spnt [3], normal [3];
  void *mgobj = mgobj(con),
       *sgobj = sgobj(con);
  SHAPE *mshp = mshp(con),
	*sshp = sshp(con);
  int state;

  /* current spatial points and normal */
  BODY_Cur_Point (con->master, con->msgp, con->mpnt, mpnt);
  BODY_Cur_Point (con->slave, con->ssgp, con->spnt, spnt);
  COPY (con->base+6, normal);

  /* update contact data => during an update 'master' and 'slave' relation does not change */
  state = gobjcontact (
    CONTACT_UPDATE, con->paircode,
    sshp, sgobj, mshp, mgobj, /* the slave body holds the outward normal */
    spnt, mpnt, normal, &con->gap, /* 'mpnt' and 'spnt' are updated here */
    &con->area, con->spair); /* surface pair might change though */

  if (state || (con->state & CON_COHESIVE))
  {
    if (con->state & CON_COHESIVE) /* reuse original point */
    {
      BODY_Cur_Point (con->master, con->msgp, con->mpnt, con->point);
      localbase (normal, con->base); /* but update normal (the interface might rotate) */
    }
    else
    {
      if (con->gap <= dom->depth) dom->flags |= DOM_DEPTH_VIOLATED;

      COPY (mpnt, con->point);
      BODY_Ref_Point (con->master, con->msgp, mpnt, con->mpnt);
      BODY_Ref_Point (con->slave, con->ssgp, spnt, con->spnt);
      localbase (normal, con->base);
      if (state > 1) /* surface pair has changed */
      {
	SURFACE_MATERIAL *mat = SPSET_Find (dom->sps, con->spair [0], con->spair [1]); /* find new surface pair description */
	con->state |= SURFACE_MATERIAL_Transfer (dom->time, mat, &con->mat); /* transfer surface pair data from the database to the local variable */
      }
    }
  }
  else
  {
#if MPI
    ext_to_remove (dom, con); /* schedule remote deletion of external constraints */
#endif
    DOM_Remove_Constraint (dom, con); /* remove from the domain */
  }
}

/* update fixed point data */
static void update_fixpnt (DOM *dom, CON *con)
{
  BODY_Cur_Point (con->master, con->msgp, con->mpnt, con->point);
}

/* update fixed direction data */
static void update_fixdir (DOM *dom, CON *con)
{
  BODY_Cur_Point (con->master, con->msgp, con->mpnt, con->point);
}

/* update velocity direction data */
static void update_velodir (DOM *dom, CON *con)
{
  VELODIR (con->Z) = TMS_Value (con->tms, dom->time + dom->step);
  BODY_Cur_Point (con->master, con->msgp, con->mpnt, con->point);
}

/* update rigid link data */
static void update_riglnk (DOM *dom, CON *con)
{
  double n [3],
	 m [3],
	 s [3],
	 len;

  if (con->master && con->slave)
  {
    BODY_Cur_Point (con->master, con->msgp, con->mpnt, m);
    BODY_Cur_Point (con->slave, con->ssgp, con->spnt, s);
  }
  else /* master point to a spatial point link */
  {
    BODY_Cur_Point (con->master, con->msgp, con->mpnt, m);
    COPY (con->spnt, s);
  }

  COPY (m, con->point);
  SUB (m, s, RIGLNK_VEC (con->Z));
  SUB (m, s, n);
  len = LEN (n);
  con->gap = len - RIGLNK_LEN(con->Z);
  len = 1.0 / len;
  SCALE (n, len);
  localbase (n, con->base);
}

/* update glue data */
static void update_glue (DOM *dom, CON *con)
{
  BODY_Cur_Point (con->master, con->msgp, con->mpnt, con->point);
}

/* tell whether the geometric objects are topologically adjacent */
static int gobj_adjacent (short paircode, void *aobj, void *bobj)
{
  switch (paircode)
  {
    case AABB_ELEMENT_ELEMENT: return ELEMENT_Adjacent (aobj, bobj);
    case AABB_CONVEX_CONVEX: return CONVEX_Adjacent (aobj, bobj);
  }

  return 0;
}

#if MPI
/* return an SGP index:
 * a) semi-positive index indicates regular surface SGP
 * b) negative index indicates a node of a non-surface FEM mesh element */
static int SGP_index (BODY *bod, SGP *sgp)
{
  long n = sgp - bod->sgp;

  if (n < 0 || n >= bod->nsgp) /* non-surface ELEMENT */
  {
    ASSERT_DEBUG (bod->kind == FEM && !bod->msh, "Regular FEM body expected");
    n = - (int) (long) sgp->box; /* GLUE-ed mesh node index (1-based) */
  }

  return n;
}
#endif

/* recreate an SGP from an index returned by SGP_index */
static SGP* SGP_from_index (DOM *dom, BODY *bod, int n)
{
  SGP *sgp;

  if (n >= 0 && n < bod->nsgp) sgp = &bod->sgp [n];
  else
  {
    ASSERT_DEBUG (bod->kind == FEM && !bod->msh, "Regular FEM body expected");
    ASSERT_DEBUG (n < 0, "Negative index expected");
    ERRMEM (sgp = MEM_Alloc (&dom->sgpmem));
    sgp->shp = bod->shape;
    ASSERT_DEBUG_EXT (sgp->gobj = MESH_Element_With_Node (bod->shape->data, -n-1), "Element with given node number not found");
    sgp->box = (BOX*) (long) (-n); /* GLUE-ed mesh node index (1-based) */
  }

  return sgp;
}

#if MPI
/* constraint weight */
static int constraint_weight (CON *con)
{
  int wgt0 = 1, /* default weight */
      wgt1 = 0; /* weight of local dynamics row */
  DOM *dom = con->master->dom;

  if (con->dia)
  {
    OFFB *blk;

    for (blk = con->dia->adjext; blk; blk = blk->n) wgt1 ++;
    for (blk = con->dia->adj; blk; blk = blk->n) wgt1 ++;
  }

  wgt0 += (int) (dom->weight_factor * (double) wgt1);

  return wgt0;
}

/* body weight */
static int body_weight (BODY *bod)
{
  return  bod->dofs + bod->nsgp; /* XXX: this is meant to represent the time integration and contact detection together */
}

/* domain weight */
static int domain_weight (DOM *dom)
{
  int weight = 0;
  BODY *bod;
  CON *con;

  for (bod = dom->bod; bod; bod = bod->next)
  {
    weight += body_weight (bod);
  }

  for (con = dom->con; con; con = con->next)
  {
    weight += constraint_weight (con);
  }

  return weight;
}

/* number of objects for balacing */
static int obj_count (DOM *dom, int *ierr)
{
  *ierr = ZOLTAN_OK;

  if (dom->ncon + dom->nbod == 0) return 1; /* XXX: Zoltan fails for 0 count */
  else return dom->ncon + dom->nbod;
}

/* list of object identifiers for load balancing */
static void obj_list (DOM *dom, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr)
{
  BODY *bod;
  CON *con;
  int i;
  
  for (con = dom->con, i = 0; con; i ++, con = con->next)
  {
    global_ids [i * num_gid_entries] = con->id;
    obj_wgts [i * wgt_dim] = constraint_weight (con);
  }

  for (bod = dom->bod; bod; i ++, bod = bod->next)
  {
    global_ids [i * num_gid_entries] = UINT_MAX - bod->id; /* XXX: vournable */
    obj_wgts [i * wgt_dim] = body_weight (bod);
  }

  if (i == 0) /* XXX: Zoltan workaround */
  {
    global_ids [0] = UINT_MAX;
    obj_wgts [0] = 1.0;
  }

  *ierr = ZOLTAN_OK;
}

/* number of spatial dimensions */
static int dimensions (DOM *dom, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}

/* list of object points exploited during load balancing */
static void obj_points (DOM *dom, int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int num_dim, double *geom_vec, int *ierr)
{
  unsigned int id;
  double *e, *v;
  BODY *bod;
  CON *con;
  int i;

  if (num_obj == 1 && global_ids [0] == UINT_MAX) /* XXX: Zoltan workaround */
  {
    SET (geom_vec, 0.0);
  }
  else for (i = 0; i < num_obj; i ++)
  {
    id = global_ids [i * num_gid_entries];
    v = &geom_vec [i* num_dim];

    con = MAP_Find (dom->idc, (void*) (long) id, NULL);

    if (con)
    {
      COPY (con->point, v);
    }
    else
    {
      ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) (UINT_MAX - id), NULL), "Invalid body id");
      e = bod->extents;
      MID (e, e+3, v);
    }
  }

  *ierr = ZOLTAN_OK;
}

/* pack constraint migrating out during */
static void pack_constraint (CON *con, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, con->id);

  /* external ranks */
  pack_int (isize, i, ints, SET_Size (con->ext));
  for (SET *item = SET_First (con->ext); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int) (long) item->data);

  pack_int (isize, i, ints, con->kind);
  pack_int (isize, i, ints, con->master->id);
  if (con->slave) pack_int (isize, i, ints, con->slave->id);
  else pack_int (isize, i, ints, 0);

  pack_int (isize, i, ints, SGP_index (con->master, con->msgp));
  if (con->slave) pack_int (isize, i, ints, SGP_index (con->slave, con->ssgp));

  pack_doubles (dsize, d, doubles, con->mpnt, 3);
  if (con->slave) pack_doubles (dsize, d, doubles, con->spnt, 3);

  pack_doubles (dsize, d, doubles, con->R, 3);
  pack_doubles (dsize, d, doubles, con->U, 3);
  pack_doubles (dsize, d, doubles, con->point, 3);
  pack_doubles (dsize, d, doubles, con->base, 9);
  pack_double  (dsize, d, doubles, con->gap);

  switch ((int) con->kind)
  {
    case CONTACT:
    pack_double  (dsize, d, doubles, con->area);
    pack_int (isize, i, ints, con->paircode);
    SURFACE_MATERIAL_Pack_State (&con->mat, dsize, d, doubles, isize, i, ints);
    break;
    case VELODIR:
    TMS_Pack (con->tms, dsize, d, doubles, isize, i, ints);
    pack_doubles (dsize, d, doubles, con->Z, DOM_Z_SIZE);
    break;
    case RIGLNK:
    pack_doubles (dsize, d, doubles, con->Z, DOM_Z_SIZE);
    break;
  }

  con->state |= CON_IDLOCK; /* prevent id deletion */

  DOM_Remove_Constraint (con->master->dom, con); /* remove from the domain */
}

/* unpack constraint migrated in during load balancing */
static void unpack_constraint (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int kind, cid, mid, sid, n, j, k;
  BODY *master, *slave;
  SGP *msgp, *ssgp;
  SET *ext;
  CON *con;

  cid = unpack_int (ipos, i, ints);

  /* external ranks */
  j = unpack_int (ipos, i, ints);
  for (ext = NULL, n = 0; n < j; n ++)
  {
    k = unpack_int (ipos, i, ints);

    if (k == dom->rank) /* migrated in plance of a previous external constraint */
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) cid, NULL), "Invalid constraint id");
      DOM_Remove_Constraint (dom, con); /* ranks packing and this removal are before constraint creation not to mess up bod->con set comparisons */
    }
    else SET_Insert (&dom->setmem, &ext, (void*) (long) k, NULL);
  }

  kind = unpack_int (ipos, i, ints);
  mid = unpack_int (ipos, i, ints);
  sid = unpack_int (ipos, i, ints);

  ASSERT_DEBUG_EXT (master = MAP_Find (dom->allbodies, (void*) (long) mid, NULL), "Invalid body id");
  if (sid) ASSERT_DEBUG_EXT (slave = MAP_Find (dom->allbodies, (void*) (long) sid, NULL), "Invalid body id"); else slave = NULL;

  n = unpack_int (ipos, i, ints);
  msgp = SGP_from_index (dom, master, n);

  if (slave) n = unpack_int (ipos, i, ints), ssgp = SGP_from_index (dom, slave, n); else ssgp = NULL;

  dom->noid = cid; /* disable constraint ids generation and use 'noid' instead */

  con = insert (dom, master, slave, msgp, ssgp, kind);

  dom->noid = 0; /* enable constraint ids generation */

  con->ext = ext; /* external ranks */

  unpack_doubles (dpos, d, doubles, con->mpnt, 3);
  if (slave) unpack_doubles (dpos, d, doubles, con->spnt, 3);

  unpack_doubles (dpos, d, doubles, con->R, 3);
  unpack_doubles (dpos, d, doubles, con->U, 3);
  unpack_doubles (dpos, d, doubles, con->point, 3);
  unpack_doubles (dpos, d, doubles, con->base, 9);
  con->gap = unpack_double  (dpos, d, doubles);

  switch (kind)
  {
    case CONTACT:
    con->area = unpack_double  (dpos, d, doubles);
    con->paircode = unpack_int (ipos, i, ints);
    SURFACE_MATERIAL_Unpack_State (dom->sps, &con->mat, dpos, d, doubles, ipos, i, ints);
    break;
    case VELODIR:
    con->tms = TMS_Unpack (dpos, d, doubles, ipos, i, ints);
    unpack_doubles (dpos, d, doubles, con->Z, DOM_Z_SIZE);
    break;
    case RIGLNK:
    unpack_doubles (dpos, d, doubles, con->Z, DOM_Z_SIZE);
    break;
  }

  con->dia = LOCDYN_Insert (dom->ldy, con, con->master, con->slave); /* insert into local dynamics */
}

/* insert a new external constraint migrated in during domain gluing */
static CON* insert_external_constraint (DOM *dom, BODY *master, BODY *slave, SGP *msgp, SGP *ssgp, short kind, unsigned int cid)
{
  CON *con;

  ERRMEM (con = MEM_Alloc (&dom->conmem));
  con->kind = kind;
  con->master = master;
  con->slave = slave;
  con->msgp = msgp;
  con->ssgp = ssgp;

  /* constraint identifier */
  con->id = cid;

  /* state */
  con->state |= CON_EXTERNAL;

#if PARDEBUG
  ASSERT_DEBUG (!SET_Contains (master->con, con, CONCMP), "Constraint %s with id %d already in body list", CON_Kind (con), con->id);
  ASSERT_DEBUG (!slave || (slave && !SET_Contains (slave->con, con, CONCMP)), "Constraint %s with id %d already in body list", CON_Kind (con), con->id);
#endif

  /* add to the body constraint adjacency */
  SET_Insert (&dom->setmem, &master->con, con, CONCMP);
  if (slave) SET_Insert (&dom->setmem, &slave->con, con, CONCMP);

#if PARDEBUG
  ASSERT_DEBUG (SET_Contains (master->con, con, CONCMP), "Failed to insert constraint %s with id %d into body list", CON_Kind (con), con->id);
  ASSERT_DEBUG (!slave || (slave && SET_Contains (slave->con, con, CONCMP)), "Failed to insert constraint %s with id %d already into body list", CON_Kind (con), con->id);
#endif
 
  /* insert into external map */
  MAP_Insert (&dom->mapmem, &dom->conext, (void*) (long) cid, con, NULL);

  return con;
}

/* pack boundary constraint migrating out during domain gluing */
static void pack_boundary_constraint (CON *con, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, con->kind);

  pack_int (isize, i, ints, con->id);
  pack_int (isize, i, ints, con->master->id);

  if (con->slave) pack_int (isize, i, ints, con->slave->id);
  else pack_int (isize, i, ints, 0);

  pack_int (isize, i, ints, SGP_index (con->master, con->msgp));
  if (con->slave) pack_int (isize, i, ints, SGP_index (con->slave, con->ssgp));

  pack_doubles (dsize, d, doubles, con->mpnt, 3);
  if (con->slave) pack_doubles (dsize, d, doubles, con->spnt, 3);

  pack_doubles (dsize, d, doubles, con->R, 3);
  pack_doubles (dsize, d, doubles, con->U, 3);
  pack_doubles (dsize, d, doubles, con->point, 3);
  pack_doubles (dsize, d, doubles, con->base, 9);

  if (con->kind == CONTACT) /* sparsification || BBS need below data  */
  {
    pack_double (dsize, d, doubles, con->area);
    pack_double (dsize, d, doubles, con->gap);
  }
}

/* unpack external constraint migrated in during domain gluing */
static CON* unpack_external_constraint (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int kind, cid, mid, sid, n;
  BODY *master, *slave;
  SGP *msgp, *ssgp;
  CON *con;

  kind = unpack_int (ipos, i, ints);

  cid = unpack_int (ipos, i, ints);
  mid = unpack_int (ipos, i, ints);
  sid = unpack_int (ipos, i, ints);

  ASSERT_DEBUG_EXT (master = MAP_Find (dom->allbodies, (void*) (long) mid, NULL), "Invalid body id");
  if (sid) ASSERT_DEBUG_EXT (slave = MAP_Find (dom->allbodies, (void*) (long) sid, NULL), "Invalid body id"); else slave = NULL;

  n = unpack_int (ipos, i, ints);
  msgp = SGP_from_index (dom, master, n);

  if (slave)
  {
    n = unpack_int (ipos, i, ints);
    ssgp = SGP_from_index (dom, slave, n);
  }
  else ssgp = NULL;

  con = insert_external_constraint (dom, master, slave, msgp, ssgp, kind, cid);

  unpack_doubles (dpos, d, doubles, con->mpnt, 3);
  if (slave) unpack_doubles (dpos, d, doubles, con->spnt, 3);

  unpack_doubles (dpos, d, doubles, con->R, 3);
  unpack_doubles (dpos, d, doubles, con->U, 3);
  unpack_doubles (dpos, d, doubles, con->point, 3);
  unpack_doubles (dpos, d, doubles, con->base, 9);

  if (kind == CONTACT) /* sparsification || BBS need below data */
  {
    con->area = unpack_double (dpos, d, doubles);
    con->gap = unpack_double (dpos, d, doubles);
  }

  return con;
}

/* pack boundary constraint update */
static void pack_boundary_constraint_update (CON *con, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, con->id);
  pack_doubles (dsize, d, doubles, con->R, 3);
  pack_doubles (dsize, d, doubles, con->U, 3);
  pack_doubles (dsize, d, doubles, con->point, 3);

  if (con->kind == CONTACT) /* sparsification || BBS need below data */
  {
    pack_doubles (dsize, d, doubles, con->mpnt, 3);
    pack_doubles (dsize, d, doubles, con->spnt, 3);
    pack_doubles (dsize, d, doubles, con->base, 9);
    pack_double (dsize, d, doubles, con->area);
    pack_double (dsize, d, doubles, con->gap);
  }
  else if (con->kind == RIGLNK)
  {
    pack_doubles (dsize, d, doubles, con->base, 9);
  }
}

/* unpack external constraint update */
static CON* unpack_external_constraint_update (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  CON *con;
  int id;

  id = unpack_int (ipos, i, ints);
  ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid constraint id");
  unpack_doubles (dpos, d, doubles, con->R, 3);
  unpack_doubles (dpos, d, doubles, con->U, 3);
  unpack_doubles (dpos, d, doubles, con->point, 3);

  if (con->kind == CONTACT) /* sparsification || BBS need below data */
  {
    unpack_doubles (dpos, d, doubles, con->mpnt, 3);
    unpack_doubles (dpos, d, doubles, con->spnt, 3);
    unpack_doubles (dpos, d, doubles, con->base, 9);
    con->area = unpack_double (dpos, d, doubles);
    con->gap = unpack_double (dpos, d, doubles);
  }
  else if (con->kind == RIGLNK)
  {
    unpack_doubles (dpos, d, doubles, con->base, 9);
  }

  return con;
}

/* pack migrating out parent body */
static void pack_parent (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  DOM *dom;

  /* must be parent */
  ASSERT_DEBUG (bod->flags & BODY_PARENT, "Not a parent");

  /* set domain */
  dom = bod->dom;

  /* pack id */
  pack_int (isize, i, ints, bod->id);

  /* pack state */
#if LOCAL_BODIES
  BODY_Pack (bod, dsize, d, doubles, isize, i, ints);
#endif
  BODY_Parent_Pack (bod, dsize, d, doubles, isize, i, ints);

  /* delete from label map */
  if (bod->label) MAP_Delete (&dom->mapmem, &dom->lab, bod->label, (MAP_Compare)strcmp);

  /* delete from id based map */
  MAP_Delete (&dom->mapmem, &dom->idb, (void*) (long) bod->id, NULL);

  /* remove from list */
  if (bod->prev) bod->prev->next = bod->next;
  else dom->bod = bod->next;
  if (bod->next) bod->next->prev = bod->prev;

  /* decrement */
  dom->nbod --;

  /* unmark parent */
  bod->flags &= ~BODY_PARENT;
}

/* unpack migrated in parent body */
static void unpack_parent (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  BODY *bod;
  int id;

  /* unpack id */
  id = unpack_int (ipos, i, ints);

  /* find body */
#if LOCAL_BODIES
  bod = MAP_Find (dom->allbodies, (void*) (long) id, NULL);

  if (bod == NULL)
  {
    bod = BODY_Unpack (dom->solfec, dpos, d, doubles, ipos, i, ints);
    bod->dom = dom;
    MAP_Insert (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, bod, NULL);
  }
  else
  {
    BODY *tmp = BODY_Unpack (dom->solfec, dpos, d, doubles, ipos, i, ints);
    BODY_Destroy (tmp);
  }
#else
  ASSERT_DEBUG_EXT (bod = MAP_Find (dom->allbodies, (void*) (long) id, NULL), "Invalid body id");
#endif

  /* must be child or dummy */
  ASSERT_DEBUG ((bod->flags & BODY_PARENT) == 0, "Neither child nor dummy");

  /* unpack state */
  BODY_Parent_Unpack (bod, dpos, d, doubles, ipos, i, ints);

  /* if it was a child */
  if (bod->flags & BODY_CHILD)
  {
    /* unmark child */
    bod->flags &= ~BODY_CHILD;

    /* delete from children set */
    SET_Delete (&dom->setmem, &dom->children, bod, NULL);
  }

  /* insert into label map */
  if (bod->label) MAP_Insert (&dom->mapmem, &dom->lab, bod->label, bod, (MAP_Compare)strcmp);

  /* insert into id based map */
  MAP_Insert (&dom->mapmem, &dom->idb, (void*) (long) bod->id, bod, NULL);

  /* insert into list */
  bod->prev = NULL;
  bod->next = dom->bod;
  if (dom->bod) dom->bod->prev = bod;
  dom->bod = bod;

  /* increment */
  dom->nbod ++;

  /* update rank */
  bod->rank = dom->rank;

  /* mark as parent */
  bod->flags |= BODY_PARENT;
}

/* pack migrating out child body */
static void pack_child (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  /* must be an exported or an existing parent */
  ASSERT_DEBUG (((bod->flags & (BODY_PARENT|BODY_CHILD)) == 0 && bod->rank != bod->dom->rank) /* just migrating out parent after being packed (hence unmarked) */
                || (bod->flags & BODY_PARENT), "Not a parent"); /* or an existing parent */

  /* pack id */
  pack_int (isize, i, ints, bod->id);

  /* pack state */
#if LOCAL_BODIES
  BODY_Pack (bod, dsize, d, doubles, isize, i, ints);
#endif
  BODY_Child_Pack (bod, dsize, d, doubles, isize, i, ints);
}

/* unpack child body */
static void unpack_child (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  BODY *bod;
  int id;

  /* unpack id */
  id = unpack_int (ipos, i, ints);

  /* find body */
#if LOCAL_BODIES
  bod = MAP_Find (dom->allbodies, (void*) (long) id, NULL);

  if (bod == NULL)
  {
    bod = BODY_Unpack (dom->solfec, dpos, d, doubles, ipos, i, ints);
    bod->dom = dom;
    MAP_Insert (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, bod, NULL);
  }
  else
  {
    BODY *tmp = BODY_Unpack (dom->solfec, dpos, d, doubles, ipos, i, ints);
    BODY_Destroy (tmp);
  }
#else
  ASSERT_DEBUG_EXT (bod = MAP_Find (dom->allbodies, (void*) (long) id, NULL), "Invalid body id");
#endif

  /* must be child or dummy */
  ASSERT_DEBUG ((bod->flags & BODY_PARENT) == 0, "Neither child nor dummy");

  /* unpack state */
  BODY_Child_Unpack (bod, dpos, d, doubles, ipos, i, ints);

  /* if it was a dummy */
  if ((bod->flags & BODY_CHILD) == 0)
  {
    /* insert into children set */
    SET_Insert (&dom->setmem, &dom->children, bod, NULL);

    /* mark as child */
    bod->flags |= BODY_CHILD;
  }

  /* mark as updated */
  bod->flags |= BODY_CHILD_UPDATED;
}

/* pack child update */
static void pack_child_update (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  /* must be an existing parent */
  ASSERT_DEBUG (bod->flags & BODY_PARENT, "Not a parent");

  /* pack id */
  pack_int (isize, i, ints, bod->id);

  /* pack state */
  BODY_Child_Update_Pack (bod, dsize, d, doubles, isize, i, ints);
}

/* unpack child update */
static void unpack_child_update (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  BODY *bod;
  int id;

  /* unpack id */
  id = unpack_int (ipos, i, ints);

  /* find body */
  ASSERT_DEBUG_EXT (bod = MAP_Find (dom->allbodies, (void*) (long) id, NULL), "Invalid body id");

  /* must be a child */
  ASSERT_DEBUG (bod->flags & BODY_CHILD, "Not a child");

  /* unpack state */
  BODY_Child_Update_Unpack (bod, dpos, d, doubles, ipos, i, ints);
}

/* pack statistics */
static void pack_stats (DOM *dom, int rank, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  DOMSTATS *s, *e;

  pack_int (isize, i, ints, dom->nbod);
  pack_int (isize, i, ints, dom->aabb->boxnum);
  pack_int (isize, i, ints, dom->ncon);
  pack_int (isize, i, ints, MAP_Size (dom->conext));
  pack_int (isize, i, ints, dom->nspa);
  pack_int (isize, i, ints, dom->bytes);
  pack_int (isize, i, ints, dom->weight);

  if (rank == (dom->ncpu-1)) /* last set was packed => zero current statistics record */
  {
    for (s = dom->stats, e = s + dom->nstats; s < e; s ++)
    {
      s->sum = s->max = 0;
      s->min = INT_MAX;
    }
  }
}

/* unpack statistics */
static void unpack_stats (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  DOMSTATS *s, *e;
  int val;

  for (s = dom->stats, e = s + dom->nstats; s < e; s ++)
  {
    val = unpack_int (ipos, i, ints);
    s->sum += val;
    s->min = MIN (s->min, val);
    s->max = MAX (s->max, val);
  }
}

/* create statistics */
static void stats_create (DOM *dom)
{
  dom->nstats = 7;

  ERRMEM (dom->stats = MEM_CALLOC (sizeof (DOMSTATS [dom->nstats])));
  
  dom->stats [0].name = "BODIES";
  dom->stats [1].name = "BOXES";
  dom->stats [2].name = "CONSTRAINTS";
  dom->stats [3].name = "EXTERNAL";
  dom->stats [4].name = "SPARSIFIED";
  dom->stats [5].name = "BYTES SENT";
  dom->stats [6].name = "DOM WEIGHT";
}

/* compute statistics */
static void stats_compute (DOM *dom)
{
  DOMSTATS *s, *e;

  for (s = dom->stats, e = s + dom->nstats; s < e; s ++)
  {
    s->avg = s->sum / dom->ncpu;
  }
}

/* destroy statistics */
static void stats_destroy (DOM *dom)
{
  free (dom->stats);
}

/* update con->point members of external constraints */
static void update_external_constraint_points (DOM *dom)
{
  int nsend, nrecv, *isize, *dsize, i, *j, *k;
  COMDATA *send, *recv, *ptr;
  double *p;
  SET *item;
  CON *con;

  nsend = dom->ncpu;
  ERRMEM (send = MEM_CALLOC (sizeof (COMDATA [nsend])));
  ERRMEM (isize = MEM_CALLOC (sizeof (int [nsend])));
  ERRMEM (dsize = MEM_CALLOC (sizeof (int [nsend])));

  for (i = 0; i < nsend; i ++) send [i].rank = i;

  for (con = dom->con; con; con = con->next)
  {
    for (item = SET_First (con->ext); item; item = SET_Next (item))
    {
      i = (int) (long) item->data;
      ptr = &send [i];
      pack_int (&isize [i], &ptr->i, &ptr->ints, con->id);
      pack_doubles (&dsize [i], &ptr->d, &ptr->doubles, con->point, 3);
    }
  }

  dom->bytes += COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);

  for (i = 0; i < nrecv; i ++)
  {
    ptr = &recv [i];
    for (j = ptr->i, k = j + ptr->ints, p = ptr->d; j < k; j ++, p += 3)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*j), NULL), "Invalid constraint id");
      COPY (p, con->point);
    }
  }

  for (i = 0; i < nsend; i ++)
  {
    ptr = &send [i];
    free (ptr->i);
    free (ptr->d);
  }
  free (send);
  free (isize);
  free (dsize);
  free (recv); /* includes recv[]->i and recv[]->d memory */
}

/* pack children update data */
static void update_children_pack (DBD *dbd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;

  pack_int (isize, i, ints, SET_Size (dbd->children));
  for (item = SET_First (dbd->children); item; item = SET_Next (item))
    pack_child_update (item->data, dsize, d, doubles, isize, i, ints);
}

/* unpack children udate data */
static void* update_children_unpack (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j;

  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    unpack_child_update (dom, dpos, d, doubles, ipos, i, ints);
  }

  return NULL;
}

/* update children shapes */
static void update_children (DOM *dom)
{
  COMOBJ *send, *recv;
  int i, nrecv;
  BODY *bod;
  SET *item;
  DBD *dbd;

  dbd = dom->dbd;

  for (bod = dom->bod; bod; bod = bod->next)
  {
    for (item = SET_First (bod->children); item; item = SET_Next (item))
    {
      SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].children, bod, NULL); /* map bodies to child rank sets */
    }
  }

  ERRMEM (send = malloc (sizeof (COMOBJ [dom->ncpu])));

  for (i = 0; i < dom->ncpu; i ++)
  {
    send [i].o = &dbd [i];
    send [i].rank = i;
  }

  /* send children updates; since this is the first communication in a sequence, we have here dom->bytes = ... rather than dom->bytes += ... */
  dom->bytes = COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)update_children_pack, dom, (OBJ_Unpack)update_children_unpack, send, dom->ncpu, &recv, &nrecv);

  for (i = 0; i < dom->ncpu; i ++) SET_Free (&dom->setmem, &dbd [i].children);
  free (send);
  free (recv);
}

/* compute ranks of migrating children */
static void children_migration_begin (DOM *dom, DBD *dbd)
{
  int *procs, numprocs, i;
  BODY *bod;

  ERRMEM (procs = malloc (sizeof (int [dom->ncpu])));

  for (bod = dom->bod; bod; bod = bod->next)
  {
    /* must be a parent */
    ASSERT_DEBUG (bod->flags & BODY_PARENT, "Not a parent");

    double *e = bod->extents;

    SET_Free (&dom->setmem, &bod->children); /* empty children set */

    Zoltan_LB_Box_Assign (dom->zol, e[0], e[1], e[2], e[3], e[4], e[5], procs, &numprocs);

    for (i = 0; i < numprocs; i ++)
    {
      if (bod->rank != procs [i]) /* if this is neither current nor the new body rank */
      {
	SET_Insert (&dom->setmem, &dbd [procs [i]].children, bod, NULL); /* schedule for sending a child */
	SET_Insert (&dom->setmem, &bod->children, (void*) (long) procs [i], NULL); /* extend parent's children set */
      }
    }
  }

  free (procs);
}

/* delete migrated out children */
static void children_migration_end (DOM *dom)
{
  SET *delset, *item;

  delset = NULL;

  for (item = SET_First (dom->children); item; item = SET_Next (item))
  {
    BODY *bod = item->data;

    /* must be a child */
    ASSERT_DEBUG (bod->flags & BODY_CHILD, "Not a child");

    if ((bod->flags & BODY_CHILD_UPDATED) == 0) /* migrated out as it wasn't updated by a parent */
    {
      bod->flags &= ~BODY_CHILD; /* unmark child */
      
      SET_Insert (&dom->setmem, &delset, bod, NULL); /* schedule deletion from dom->children */
    }
    else bod->flags &= ~BODY_CHILD_UPDATED; /* invalidate update flag */
  }

  /* subtract deleted children from domain children set */
  for (item = SET_First (delset); item; item = SET_Next (item))
  {
    SET_Delete (&dom->setmem, &dom->children, item->data, NULL);
  }

  SET_Free (&dom->setmem, &delset);
}

#if LOCAL_BODIES
/* pack parents and children */
static void pack_parents_and_children (DBD *dbd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;

  /* pack exported bodies */
  pack_int (isize, i, ints, SET_Size (dbd->bodies));
  for (item = SET_First (dbd->bodies); item; item = SET_Next (item))
    pack_parent (item->data, dsize, d, doubles, isize, i, ints);

  /* pack exported children */
  pack_int (isize, i, ints, SET_Size (dbd->children));
  for (item = SET_First (dbd->children); item; item = SET_Next (item))
    pack_child (item->data, dsize, d, doubles, isize, i, ints);
}

/* unpack parents and children */
static void* unpack_parents_and_children (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j;

  /* unpack imported bodies */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    unpack_parent (dom, dpos, d, doubles, ipos, i, ints);
  }

  /* unpack imported children */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    unpack_child (dom, dpos, d, doubles, ipos, i, ints);
  }

  return NULL;
}
#endif

/* pack domain balancing data */
static void domain_balancing_pack (DBD *dbd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;

#if LOCAL_BODIES == 0
  /* pack exported bodies */
  pack_int (isize, i, ints, SET_Size (dbd->bodies));
  for (item = SET_First (dbd->bodies); item; item = SET_Next (item))
    pack_parent (item->data, dsize, d, doubles, isize, i, ints);

  /* pack exported children */
  pack_int (isize, i, ints, SET_Size (dbd->children));
  for (item = SET_First (dbd->children); item; item = SET_Next (item))
    pack_child (item->data, dsize, d, doubles, isize, i, ints);
#endif

  /* pack exported constraints */
  pack_int (isize, i, ints, SET_Size (dbd->constraints));
  for (item = SET_First (dbd->constraints); item; item = SET_Next (item))
    pack_constraint (item->data, dsize, d, doubles, isize, i, ints);

  /* pack ids of deleted external constraints */
  pack_int (isize, i, ints, SET_Size (dbd->remove));
  for (item = SET_First (dbd->remove); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int) (long) item->data);
}

/* unpack domain balancing data */
static void* domain_balancing_unpack (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j, id;
  CON *con;

#if LOCAL_BODIES == 0
  /* unpack imported bodies */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    unpack_parent (dom, dpos, d, doubles, ipos, i, ints);
  }

  /* unpack imported children */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    unpack_child (dom, dpos, d, doubles, ipos, i, ints);
  }
#endif

  /* unpack imporeted constraints */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    unpack_constraint (dom, dpos, d, doubles, ipos, i, ints);
  }

  /* unpack deleted external constraint ids */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    id = unpack_int (ipos, i, ints);
    ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid constraint id");
    DOM_Remove_Constraint (dom, con);
  }

  return NULL;
}

/* compute migration of old boundary constraints (those remaining after constraints update) */
static void old_boundary_constraints_migration (DOM *dom, DBD *dbd)
{
  SET *item, *ext;
  BODY *bod;
  CON *con;
  int i;

  /* compute migration sets */
  for (con = dom->con; con; con = con->next)
  {
    BODY *bodies [] = {con->master, con->slave};

    ext = NULL;

    for (i = 0; i < 2; i ++)
    {
      bod = bodies [i];

      if (bod)
      {
	for (item = SET_First (bod->children); item; item = SET_Next (item))
	{
	  SET_Insert (&dom->setmem, &ext, item->data, NULL);
	}
      
	if (bod->flags & BODY_CHILD)
	{
	  SET_Insert (&dom->setmem, &ext, (void*) (long) bod->rank, NULL);
	}
      }
    }

    for (item = SET_First (ext); item; item = SET_Next (item))
    {
      if (!SET_Contains (con->ext, item->data, NULL)) /* not exported yet */
      {
	SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].glue, con, NULL); /* schedule export */
	SET_Insert (&dom->setmem, &con->ext, item->data, NULL);
      }
      else /* already exported */
      {
	SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].update, con, NULL); /* schedule update */
      }

#if LOCAL_BODIES
      for (i = 0; i < 2; i ++)
      {
	bod = bodies [i];

	if (bod)
	{
	  if (bod->flags & BODY_CHILD)
	  {
	    if (!SET_Contains (bod->children, item->data, NULL) && bod->rank != (int) (long) item->data)
	      SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].dummies, bod, NULL); /* schedule dummy */
	  }
	  else if (!SET_Contains (bod->children, item->data, NULL)) /* parent */
	    SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].dummies, bod, NULL); /* schedule dummy */
	}
      }
#endif
    }

    for (item = SET_First (con->ext); item; item = SET_Next (item))
    {
      if (!SET_Contains (ext, item->data, NULL)) /* not needed any more */
      {
	SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].remove, (void*) (long) con->id, NULL); /* schedule remote deletion */
      }
    }

    SET_Free (&dom->setmem, &con->ext);
    con->ext = ext;
  }
}

/* pack old boundary constraints */
static void old_boundary_constraints_pack (DBD *dbd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;

#if LOCAL_BODIES
  /* pack dummies */
  pack_int (isize, i, ints, SET_Size (dbd->dummies));
  for (item = SET_First (dbd->dummies); item; item = SET_Next (item))
    BODY_Pack (item->data, dsize, d, doubles, isize, i, ints);
#endif

  /* pack exported boundary constraints */
  pack_int (isize, i, ints, SET_Size (dbd->glue));
  for (item = SET_First (dbd->glue); item; item = SET_Next (item))
    pack_boundary_constraint (item->data, dsize, d, doubles, isize, i, ints);

  /* pack updated boundary constraints */
  pack_int (isize, i, ints, SET_Size (dbd->update));
  for (item = SET_First (dbd->update); item; item = SET_Next (item))
    pack_boundary_constraint_update (item->data, dsize, d, doubles, isize, i, ints);

  /* pack ids of deleted external constraints */
  pack_int (isize, i, ints, SET_Size (dbd->remove));
  for (item = SET_First (dbd->remove); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int) (long) item->data);
}

/* unpack old external constraints */
static void* old_external_constraints_unpack (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j, id;
  SET **pp;
  CON *con;

  ERRMEM (pp = MEM_CALLOC (sizeof (SET*)));

#if LOCAL_BODIES
  /* unpack dummies */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    BODY *tmp = BODY_Unpack (dom->solfec, dpos, d, doubles, ipos, i, ints),
         *bod = MAP_Find (dom->allbodies, (void*) (long) tmp->id, NULL);

    if (bod == NULL)
    {
      bod = tmp;
      bod->dom = dom;
      MAP_Insert (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, bod, NULL);
    }
    else BODY_Destroy (tmp);
  }
#endif

  /* unpack imporeted external constraints */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    con = unpack_external_constraint (dom, dpos, d, doubles, ipos, i, ints);
    SET_Insert (&dom->setmem, pp, con, NULL);
  }

  /* unpack updated external constraints */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    con = unpack_external_constraint_update (dom, dpos, d, doubles, ipos, i, ints);
    SET_Insert (&dom->setmem, pp, con, NULL);
  }

  /* unpack deleted external constraint ids */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    id = unpack_int (ipos, i, ints);
    ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid constraint id");
    DOM_Remove_Constraint (dom, con);
  }

  return pp;
}

/* insert or delete pending constraints */
static void insert_pending_constraints (DOM *dom)
{
  double d [3], *e, *p;
  PNDCON *pnd;
  SET *item;
  CON *con;

  for (item = SET_First (dom->pendingcons); item; item = SET_Next (item))
  {
    pnd = item->data;

    /* make sure that all attached body constraint points are within the body extents so that the newly
     * inserted constraints will migrate to partitions where bodies have representation (child/parent) */

    if (pnd->master->flags & BODY_PARENT)
    {
      e = pnd->master->extents;
      p = pnd->mpnt;

      if (p [0] < e [0]) e [0] = p [0];
      if (p [1] < e [1]) e [1] = p [1];
      if (p [2] < e [2]) e [2] = p [2];
      if (p [0] > e [3]) e [3] = p [0];
      if (p [1] > e [4]) e [4] = p [1];
      if (p [2] > e [5]) e [5] = p [2];
    }

    if (pnd->slave && pnd->slave->flags & BODY_PARENT)
    {
      double *pp [] = {pnd->mpnt, pnd->spnt};
      e = pnd->slave->extents;

      for (int i = 0; i < 2; i ++)
      {
        p = pp [i]; /* make sure that slave knows about both points since it can be on a processor different that
		       the parent and can have no other means of knowing where to migrate the needed child */

	if (p [0] < e [0]) e [0] = p [0];
	if (p [1] < e [1]) e [1] = p [1];
	if (p [2] < e [2]) e [2] = p [2];
	if (p [0] > e [3]) e [3] = p [0];
	if (p [1] > e [4]) e [4] = p [1];
	if (p [2] > e [5]) e [5] = p [2];
      }
    }
  }

  for (item = SET_First (dom->pendingcons); item; item = SET_Next (item))
  {
    pnd = item->data;

    /* extend extents of these bodies that might have been altered above;
     * note that we are not using GEOMETRIC_EPSILON here as this can be set
     * by the user which in turn may cause migration consitency problems */

    if (pnd->master->flags & BODY_PARENT)
    {
      e = pnd->master->extents;
      SUB (e+3, e, d);
      SCALE (d, PUT_GEOMEPS);
      SUB (e, d, e);
      ADD (e+3, d, e+3);
    }

    if (pnd->slave && pnd->slave->flags & BODY_PARENT)
    {
      e = pnd->slave->extents;
      SUB (e+3, e, d);
      SCALE (d, PUT_GEOMEPS);
      SUB (e, d, e);
      ADD (e+3, d, e+3);
    }
  }

  for (item = SET_First (dom->pendingcons); item; item = SET_Next (item))
  {
    pnd = item->data;

    if (pnd->master->flags & BODY_PARENT) /* insert only those having parent master */
    {
      switch (pnd->kind)
      {
      case FIXPNT:
	con = DOM_Fix_Point (dom, pnd->master, pnd->mpnt);
	break;
      case FIXDIR:
	con = DOM_Fix_Direction (dom, pnd->master, pnd->mpnt, pnd->dir);
	break;
      case VELODIR:
	con = DOM_Set_Velocity (dom, pnd->master, pnd->mpnt, pnd->dir, pnd->val);
	break;
      case RIGLNK:
	con = DOM_Put_Rigid_Link (dom, pnd->master, pnd->slave, pnd->mpnt, pnd->spnt);
	break;
      case GLUE:
	con = DOM_Glue_Nodes (dom, pnd->master, pnd->slave, pnd->mnode, pnd->snode);
	break;
      }

      ASSERT_TEXT (con, "Failed to insert a pending constraint.\n"
	                "Please report this bug!\n");

      switch ((int)con->kind)
      {
      case FIXPNT: update_fixpnt (dom, con); break;
      case FIXDIR: update_fixdir (dom, con); break;
      case VELODIR: update_velodir (dom, con); break;
      case RIGLNK: update_riglnk (dom, con); break;
      case GLUE: update_glue (dom, con); break;
      }
    }

    free (pnd); /* they were MEM_ALLOC-ed */
  }

  /* empty pending constraints set */
  SET_Free (&dom->setmem, &dom->pendingcons);
}

/* domain balancing */
static void domain_balancing (DOM *dom)
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

  COMOBJ *send, *recv;
  char str [128];
  int i, rank;
  int nrecv;
  SET *item;
  BODY *bod;
  DBD *dbd;
  CON *con;

  /* pending constraints */
  insert_pending_constraints (dom);

#if PARDEBUG
  /* test whether constraints attached to bodies have their points inside of body extents */
  for (bod = dom->bod; bod; bod = bod->next)
  {
    for (item = SET_First (bod->con); item; item = SET_Next (item))
    {
      con = item->data;
      if (con->point [0] < bod->extents [0] || con->point [0] > bod->extents [3] ||
	  con->point [1] < bod->extents [1] || con->point [1] > bod->extents [4] ||
	  con->point [2] < bod->extents [2] || con->point [2] > bod->extents [5])
      {
	ASSERT_DEBUG (0, "Constraint point outside of the attached body extents"); /* a debugger catchable assertion */
      }
    }
  }
#endif

  /* domain weight for statistics */
  dom->weight = domain_weight (dom);

  /* load balancing migration sets */
  dbd = dom->dbd;

#if MPI && DEBUG
  for (con = dom->con, i = 0; con; con = con->next, i ++);
  ASSERT_DEBUG (i == dom->ncon, "Inconsistent constraints count");
  for (bod = dom->bod, i = 0; bod; bod = bod->next, i ++);
  ASSERT_DEBUG (i == dom->nbod, "Inconsistent bodies count");
#endif

  if (dom->rebalanced % dom->updatefreq == 0)
  {
#if 0
    Zoltan_Generate_Files (dom->zol, "kdd", 1, 1, 0, 0);
#endif

    /* update RCB parameters */
    snprintf (str, 128, "%g", dom->imbalance_tolerance);
    Zoltan_Set_Param (dom->zol, "IMBALANCE_TOL", str);

    /* update body partitioning */
    ASSERT (Zoltan_LB_Balance (dom->zol, &changes, &num_gid_entries, &num_lid_entries,
	    &num_import, &import_global_ids, &import_local_ids, &import_procs,
	    &num_export, &export_global_ids, &export_local_ids, &export_procs) == ZOLTAN_OK, ERR_ZOLTAN);

    for (i = 0; (dom->bod || dom->con) && i < num_export; i ++)
    {
      rank = export_procs [i];

      if (rank != dom->rank)
      {
	con = MAP_Find (dom->idc, (void*) (long) export_global_ids [i], NULL);
	if (con) SET_Insert (&dom->setmem, &dbd [rank].constraints, con, NULL); /* map this constraint to its export rank */
	else
	{
	  ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) (UINT_MAX - export_global_ids [i]), NULL), "Invalid body id");
	  bod->rank = rank; /* set the new rank */
	  SET_Insert (&dom->setmem, &dbd [rank].bodies, bod, NULL); /* map this body to its export rank */
	}
      }
    }

    Zoltan_LB_Free_Data (&import_global_ids, &import_local_ids, &import_procs,
			 &export_global_ids, &export_local_ids, &export_procs);

    dom->rebalanced = 1;
  }
  else /* reuse previous geometrical partitioning */
  {
    for (bod = dom->bod; bod; bod = bod->next)
    {
      double *e = bod->extents, v [3];
      MID (e, e+3, v);
      Zoltan_LB_Point_Assign (dom->zol, v, &rank);
      if (rank != dom->rank)
      {
	bod->rank = rank;
	SET_Insert (&dom->setmem, &dbd [rank].bodies, bod, NULL);
      }
    }

    for (con = dom->con; con; con = con->next)
    {
      Zoltan_LB_Point_Assign (dom->zol, con->point, &rank);
      if (rank != dom->rank)
	SET_Insert (&dom->setmem, &dbd [rank].constraints, con, NULL);
    }

    dom->rebalanced ++;
  }

  /* --- */

  ERRMEM (send = malloc (sizeof (COMOBJ [dom->ncpu])));

  for (i = 0; i < dom->ncpu; i ++)
  {
    send [i].rank = i;
    send [i].o = &dbd [i];
  }

  /* compute chidren migration sets */
  children_migration_begin (dom, dbd);

#if LOCAL_BODIES
  /* communication */
  dom->bytes += COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)pack_parents_and_children, dom, (OBJ_Unpack)unpack_parents_and_children, send, dom->ncpu, &recv, &nrecv);
#endif

  /* communication */
  dom->bytes += COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)domain_balancing_pack, dom, (OBJ_Unpack)domain_balancing_unpack, send, dom->ncpu, &recv, &nrecv);

  /* delete migrated out children */
  children_migration_end (dom);

  /* free auxiliary sets */
  for (i = 0; i < dom->ncpu; i ++)
  {
    SET_Free (&dom->setmem, &dbd [i].bodies);
    SET_Free (&dom->setmem, &dbd [i].children);
    SET_Free (&dom->setmem, &dbd [i].constraints);
    SET_Free (&dom->setmem, &dbd [i].remove);
  }

  /* clean */
  free (recv);

#if DEBUG
  for (con = dom->con; con; con = con->next)
  {
    if ((con->master->flags & (BODY_PARENT|BODY_CHILD)) == 0 ||
	(con->slave && (con->slave->flags & (BODY_PARENT|BODY_CHILD)) == 0))
    {
      ASSERT_DEBUG (0, "Regular constraint (id = %d, %s) attached to a dummy (id = %d, %s)", /* a debugger catchable assertion */
        con->id, CON_Kind (con),
	(con->master->flags & (BODY_PARENT|BODY_CHILD)) == 0  ? con->master->id : con->slave->id,
	(con->master->flags & (BODY_PARENT|BODY_CHILD)) == 0  ? "master" : "slave");
    }
  }
#endif

  /* compute old boundary constraints migration sets */
  old_boundary_constraints_migration (dom, dbd);

  /* after this step all bodies contain sets of all old constraints (including contacts); this way during contact detection all existing contacts will get filtered out;
   * note that children sets of bodies can extend and hence constraints might need to be sent to new ranks */
  dom->bytes += COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)old_boundary_constraints_pack, dom, (OBJ_Unpack)old_external_constraints_unpack, send, dom->ncpu, &recv, &nrecv);

  /* assign external ranks */
  for (i = 0; i < nrecv; i ++)
  {
    for (item = SET_First (*(SET**) recv [i].o); item; item = SET_Next (item))
    {
      con = item->data;
      con->rank = recv [i].rank;
    }

    SET_Free (&dom->setmem, (SET**) recv [i].o);
    free (recv [i].o);
  }

  /* free auxiliary sets */
  for (i = 0; i < dom->ncpu; i ++)
  {
    SET_Free (&dom->setmem, &dbd [i].remove);
    SET_Free (&dom->setmem, &dbd [i].update);
    SET_Free (&dom->setmem, &dbd [i].glue);
#if LOCAL_BODIES
    SET_Free (&dom->setmem, &dbd [i].dummies);
#endif
  }

  /* clean */
  free (send);
  free (recv);
}

/* compute new boundary constraints migration */
static void new_boundary_constraints_migration (DOM *dom, DBD *dbd)
{
  SET *item;
  BODY *bod;
  CON *con;
  int i;

  /* compute additional boundary sets */
  for (con = dom->con; con; con = con->next)
  {
    if (con->state & CON_NEW) /* only newly created constraints */
    {
      BODY *bodies [] = {con->master, con->slave};

      for (i = 0; i < 2; i ++)
      {
	bod = bodies [i];

	for (item = SET_First (bod->children); item; item = SET_Next (item))
	{
	  SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].glue, con, NULL); /* schedule for sending to children */
	  SET_Insert (&dom->setmem, &con->ext, item->data, NULL); /* record as sent to rank */
	}

	if (bod->flags & BODY_CHILD)
	{
	  SET_Insert (&dom->setmem, &dbd [bod->rank].glue, con, NULL); /* schedule for sending to parent */
	  SET_Insert (&dom->setmem, &con->ext, (void*) (long) bod->rank, NULL); /* record as sent to rank */
	}
      }

#if LOCAL_BODIES
      for (item = SET_First (con->ext); item; item = SET_Next (item))
      {
	for (i = 0; i < 2; i ++)
	{
	  bod = bodies [i];

	  if (bod->flags & BODY_CHILD)
	  {
	    if (!SET_Contains (bod->children, item->data, NULL) && bod->rank != (int) (long) item->data)
	      SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].dummies, bod, NULL); /* schedule dummy */
	  }
	  else if (!SET_Contains (bod->children, item->data, NULL)) /* parent */
	    SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].dummies, bod, NULL); /* schedule dummy */
	}
      }
#endif
    }
  }
}

/* pack domain gluing data */
static void domain_gluing_begin_pack (DBD *dbd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;

#if LOCAL_BODIES
  /* pack dummies */
  pack_int (isize, i, ints, SET_Size (dbd->dummies));
  for (item = SET_First (dbd->dummies); item; item = SET_Next (item))
    BODY_Pack (item->data, dsize, d, doubles, isize, i, ints);
#endif

  /* pack penetration depth flag */
  pack_int (isize, i, ints, dbd->dom->flags & DOM_DEPTH_VIOLATED);

  /* pack exported boundary constraints */
  pack_int (isize, i, ints, SET_Size (dbd->glue));
  for (item = SET_First (dbd->glue); item; item = SET_Next (item))
    pack_boundary_constraint (item->data, dsize, d, doubles, isize, i, ints);
}

/* unpack domain gluing data */
static void* domain_gluing_begin_unpack (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j;
  SET **pp;
  CON *con;

  ERRMEM (pp = MEM_CALLOC (sizeof (SET*)));

#if LOCAL_BODIES
  /* unpack dummies */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    BODY *tmp = BODY_Unpack (dom->solfec, dpos, d, doubles, ipos, i, ints),
         *bod = MAP_Find (dom->allbodies, (void*) (long) tmp->id, NULL);

    if (bod == NULL)
    {
      bod = tmp;
      bod->dom = dom;
      MAP_Insert (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, bod, NULL);
    }
    else BODY_Destroy (tmp);
  }
#endif

  /* unpack penetration depth flag */
  n = unpack_int (ipos, i, ints); ASSERT (!n, ERR_DOM_DEPTH);

  /* unpack imporeted external constraints */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    con = unpack_external_constraint (dom, dpos, d, doubles, ipos, i, ints);
    SET_Insert (&dom->setmem, pp, con, NULL);
  }

  return pp;
}

/* domain gluing begin */
static void domain_gluing_begin (DOM *dom)
{
  COMOBJ *send, *recv;
  SET *item;
  int nrecv;
  CON *con;
  DBD *dbd;
  int i;

  ERRMEM (send = malloc (sizeof (COMOBJ [dom->ncpu])));
  dbd = dom->dbd;

  for (i = 0; i < dom->ncpu; i ++)
  {
    send [i].rank = i;
    send [i].o = &dbd [i];
  }

  /* compute migration sets */
  new_boundary_constraints_migration (dom, dbd);

  /* communication */
  dom->bytes += COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)domain_gluing_begin_pack, dom, (OBJ_Unpack)domain_gluing_begin_unpack, send, dom->ncpu, &recv, &nrecv);

  /* assign external ranks */
  for (i = 0; i < nrecv; i ++)
  {
    for (item = SET_First (*(SET**) recv [i].o); item; item = SET_Next (item))
    {
      con = item->data;
      con->rank = recv [i].rank;
    }

    SET_Free (&dom->setmem, (SET**) recv [i].o);
    free (recv [i].o);
  }

  /* clean */
  for (i = 0; i < dom->ncpu; i ++)
  {
    SET_Free (&dom->setmem, &dbd [i].glue);
#if LOCAL_BODIES
    SET_Free (&dom->setmem, &dbd [i].dummies);
#endif
  }
  free (send);
  free (recv);
}

/* pack domain gluing data */
static void domain_gluing_end_pack (DBD *dbd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;

  /* pack statistics */
  pack_stats (dbd->dom, dbd->rank, dsize, d, doubles, isize, i, ints);

  /* pack ids of external contacts to be deleted */
  pack_int (isize, i, ints, SET_Size (dbd->remove));
  for (item = SET_First (dbd->remove); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int) (long) item->data);
}

/* unpack domain gluing data */
static void* domain_gluing_end_unpack (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j, id;
  CON *con;

  /* unpack statistics */
  unpack_stats (dom, dpos, d, doubles, ipos, i, ints);

  /* unpack deleted external constraint ids */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    id = unpack_int (ipos, i, ints);
    ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid constraint id");
    DOM_Remove_Constraint (dom, con);
  }

  return NULL;
}

/* prepare dom->dbd[]->ext sets */
static void prepare_reaction_update_sets (DOM *dom)
{
  SET *item;
  DBD *dbd;
  CON *con;
  int i;

  dbd = dom->dbd;

  for (i = 0; i < dom->ncpu; i ++) SET_Free (&dom->setmem, &dbd [i].ext);

  for (con = dom->con; con; con = con->next)
  {
    for (item = SET_First (con->ext); item; item = SET_Next (item))
    {
      SET_Insert (&dom->setmem, &dbd [(int) (long) item->data].ext, con, NULL);
    }
  }
}

/* domain gluing end */
static void domain_gluing_end (DOM *dom)
{
  COMOBJ *send, *recv;
  int nrecv;
  DBD *dbd;
  int i;

  ERRMEM (send = malloc (sizeof (COMOBJ [dom->ncpu])));
  dbd = dom->dbd;

  for (i = 0; i < dom->ncpu; i ++)
  {
    send [i].rank = i;
    send [i].o = &dbd [i];
  }

  /* send ids of external constraints to be deleted and deleted external constraints of received ids */
  dom->bytes += COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)domain_gluing_end_pack, dom, (OBJ_Unpack)domain_gluing_end_unpack, send, dom->ncpu, &recv, &nrecv);

  /* compute statistics */
  stats_compute (dom);

  /* clean */
  for (i = 0; i < dom->ncpu; i ++) SET_Free (&dom->setmem, &dbd [i].remove);
  free (send);
  free (recv);

  /* reaction update sets */
  prepare_reaction_update_sets (dom);
}

/* pack managed bodies */
static void manage_bodies_pack (DBD *dbd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;

  /* pack spare body ids */
  pack_int (isize, i, ints, SET_Size (dbd->dom->sparebid));
  for (item = SET_First (dbd->dom->sparebid); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int) (long) item->data);

  /* pack rank */
  pack_int (isize, i, ints, dbd->dom->rank);

  if (dbd->rank != dbd->dom->rank)
  {
    /* pack pending bodies */
    pack_int (isize, i, ints, SET_Size (dbd->dom->pendingbods));
#if LOCAL_BODIES == 0
    for (item = SET_First (dbd->dom->pendingbods); item; item = SET_Next (item))
      BODY_Pack (item->data, dsize, d, doubles, isize, i, ints);
#endif
  }
}

/* unpack children udate data => NOTE, that this routine is called in the sequenc of ranks */
static void* manage_bodies_unpack (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j, k, rank;
  SET *item;

  /* unpack spare body ids */
  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    k = unpack_int (ipos, i, ints);
    SET_Insert (&dom->setmem, &dom->sparebid, (void*) (long) k, NULL); /* creates union across all ranks */
  }

  /* unpack rank */
  rank = unpack_int (ipos, i, ints);

  if (rank != dom->rank)
  {
    /* unpack pending bodies */
    j = unpack_int (ipos, i, ints);
#if LOCAL_BODIES
    dom->bid += j; /* account for the global body ID increments */
    ASSERT (dom->bid < UINT_MAX, ERR_DOM_TOO_MANY_BODIES); /* make sure we do not run out of ids */
#else
    BODY *bod;

    for (n = 0; n < j; n ++)
    {
      bod = BODY_Unpack (dom->solfec, dpos, d, doubles, ipos, i, ints);
      dom->insertbodymode = NEVER;
      DOM_Insert_Body (dom, bod); /* XXX: rely on rank ordering in communication and body
					  ordering during packing in order to get
					  the same sequence of bodies on all processors;
					  this guarantees the righ ID assignment */
    }
#endif
  }
  else
  {
    for (item = SET_First (dom->pendingbods); item; item = SET_Next (item))
    {
      dom->insertbodymode = ALWAYS;
      DOM_Insert_Body (dom, item->data); /* XXX: as above */
    }
  }

  return NULL;
}

/* deleted unwanted and insert pending bodies */
static void manage_bodies (DOM *dom)
{
  COMOBJ *send, *recv;
  int i, nrecv;
  BODY *bod;
  SET *item;
  DBD *dbd;

  dbd = dom->dbd;

  ERRMEM (send = malloc (sizeof (COMOBJ [dom->ncpu])));

  for (i = 0; i < dom->ncpu; i ++)
  {
    send [i].o = &dbd [i];
    send [i].rank = i;
  }

  /* since this is the first communication in a sequence, we have here dom->bytes = ... rather than dom->bytes += ... */
  dom->bytes = COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)manage_bodies_pack, dom, (OBJ_Unpack)manage_bodies_unpack, send, dom->ncpu, &recv, &nrecv);

  free (send);
  free (recv);

  /* empty pending bodies set */
  SET_Free (&dom->setmem, &dom->pendingbods);

  /* delete bodies associated with spare ids */
  for (item = SET_First (dom->sparebid); item; item = SET_Next (item))
  {
    if ((bod = MAP_Find (dom->allbodies, item->data, NULL)))
    {
      DOM_Remove_Body (dom, bod); /* loop over 'sparebid' => look there (***) */
      BODY_Destroy (bod);
    }
  }

  /* empty body ids set */
  SET_Free (&dom->setmem, &dom->sparebid);

  /* restore body insertion mode */
  dom->insertbodymode = EVERYNCPU;

#if LOCAL_BODIES
  /* delete local bodies which are
   * neither parent, nor child not dummy */

  SET *todel = NULL;
  MAP *jtem;

  for (jtem = MAP_First (dom->allbodies); jtem; jtem = MAP_Next (jtem))
  {
    bod = jtem->data;
    if (!((bod->flags & (BODY_PARENT|BODY_CHILD)) || bod->con))
    {
      SET_Insert (&dom->setmem, &todel, bod, NULL);
    }
  }

  for (item = SET_First (todel); item; item = SET_Next (item))
  {
    bod = item->data;
    MAP_Delete (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, NULL);
    BODY_Destroy (bod);
  }

  SET_Free (&dom->setmem, &todel);
#endif
}

/* pack normal reaction components (only for contacts) of boundary constraints */
static void pack_normal_reactions (SET *set, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;
  CON *con;

  pack_int (isize, i, ints, SET_Size (set));
  for (item = SET_First (set); item; item = SET_Next (item))
  {
    con = item->data;
    pack_int (isize, i, ints, con->id);
    if (con->kind == CONTACT) pack_double (dsize, d, doubles, con->R [2]);
    else pack_doubles (dsize, d, doubles, con->R, 3);
  }
}

/* unpack normal reaction components (only for contacts) of external constraints */
static void* unpack_normal_reactions (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j, id;
  CON *con;

  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    id = unpack_int (ipos, i, ints);
    ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid contact id");
    if (con->kind == CONTACT) con->R [2] = unpack_double (dpos, d, doubles);
    else unpack_doubles (dpos, d, doubles, con->R, 3);
  }

  return NULL;
}

/* pack boundary reactions */
static void pack_reactions (SET *set, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;
  CON *con;

  pack_int (isize, i, ints, SET_Size (set));
  for (item = SET_First (set); item; item = SET_Next (item))
  {
    con = item->data;
    pack_int (isize, i, ints, con->id);
    pack_doubles (dsize, d, doubles, con->R, 3);
  }
}

/* unpack external reactions */
static void* unpack_reactions (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j, id;
  CON *con;

  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    id = unpack_int (ipos, i, ints);
    ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid contact id");
    unpack_doubles (dpos, d, doubles, con->R, 3);
  }

  return NULL;
}

/* pack boundary R, U, V */
static void pack_RUV (SET *set, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;
  CON *con;

  pack_int (isize, i, ints, SET_Size (set));
  for (item = SET_First (set); item; item = SET_Next (item))
  {
    con = item->data;
    pack_int (isize, i, ints, con->id);
    pack_doubles (dsize, d, doubles, con->R, 3);
    pack_doubles (dsize, d, doubles, con->U, 3);
    pack_doubles (dsize, d, doubles, con->V, 3);
  }
}

/* unpack external R, U, V */
static void* unpack_RUV (DOM *dom, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int n, j, id;
  CON *con;

  j = unpack_int (ipos, i, ints);
  for (n = 0; n < j; n ++)
  {
    id = unpack_int (ipos, i, ints);
    ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid contact id");
    unpack_doubles (dpos, d, doubles, con->R, 3);
    unpack_doubles (dpos, d, doubles, con->U, 3);
    unpack_doubles (dpos, d, doubles, con->V, 3);
  }

  return NULL;
}

/* send boundary R, U, V to their external receivers; U and V are needed in bod.c: compute_contact_work;
 * note that external constraints have con->dia == NULL, hence V has been moved inside of CON structure */
void update_external_RUV (DOM *dom)
{
  COMOBJ *send, *recv;
  int i, nrecv;

  ERRMEM (send = malloc (sizeof (COMOBJ [dom->ncpu])));

  for (i = 0; i < dom->ncpu; i ++)
  {
    send [i].o = dom->dbd [i].ext;
    send [i].rank = i;
  }

  dom->bytes += COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)pack_RUV, dom,
    (OBJ_Unpack)unpack_RUV, send, dom->ncpu, &recv, &nrecv);

  free (send);
  free (recv);
}

/* create MPI related data */
static void create_mpi (DOM *dom)
{
  dom->rebalanced = 0;

  dom->updatefreq = 10;

  dom->insertbodymode = EVERYNCPU;

  dom->sparebid = NULL;

  dom->children = NULL;

  dom->conext = NULL;

  MPI_Comm_rank (MPI_COMM_WORLD, &dom->rank); /* store rank */

  MPI_Comm_size (MPI_COMM_WORLD, &dom->ncpu); /* store size */

  ERRMEM (dom->dbd  = MEM_CALLOC (sizeof (DBD [dom->ncpu])));

  for (int i = 0; i < dom->ncpu; i ++)
  {
    dom->dbd [i].dom = dom;
    dom->dbd [i].rank = i;
  }

  dom->cid = (dom->rank + 1); /* overwrite */

  dom->noid = 0; /* assign constraint ids in 'insert' routine (turned off when importing non-contacts) */

  dom->bytes = 0;

  stats_create (dom);

  ASSERT (dom->zol = Zoltan_Create (MPI_COMM_WORLD), ERR_ZOLTAN); /* zoltan context domain partitioning */

  dom->imbalance_tolerance = 1.3;
  dom->weight_factor = 1.0;

  /* general parameters */
  Zoltan_Set_Param (dom->zol, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param (dom->zol, "DEBUG_MEMORY", "0");
  Zoltan_Set_Param (dom->zol, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param (dom->zol, "NUM_LID_ENTRIES", "0");
  Zoltan_Set_Param (dom->zol, "OBJ_WEIGHT_DIM", "1");
 
  /* load balancing parameters */
  Zoltan_Set_Param (dom->zol, "LB_METHOD", "RCB");
  Zoltan_Set_Param (dom->zol, "IMBALANCE_TOL", "1.3");
  Zoltan_Set_Param (dom->zol, "AUTO_MIGRATE", "FALSE");
  Zoltan_Set_Param (dom->zol, "RETURN_LISTS", "EXPORT");

  /* RCB parameters */
  Zoltan_Set_Param (dom->zol, "RCB_OVERALLOC", "1.3");
  Zoltan_Set_Param (dom->zol, "RCB_REUSE", "1");
  Zoltan_Set_Param (dom->zol, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param (dom->zol, "CHECK_GEOM", "1");
  Zoltan_Set_Param (dom->zol, "KEEP_CUTS", "1");
  Zoltan_Set_Param (dom->zol, "REDUCE_DIMENSIONS", "0"); /* FIXME: when enabled runs into Zoltan bug with the inconsitency between
							   Zoltan_LB_Box_Assign and Zoltan_LB_Point_Assign for < 3 dimensional sets */

  /* body callbacks */
  Zoltan_Set_Fn (dom->zol, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) obj_count, dom);
  Zoltan_Set_Fn (dom->zol, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)()) obj_list, dom);
  Zoltan_Set_Fn (dom->zol, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) dimensions, dom);
  Zoltan_Set_Fn (dom->zol, ZOLTAN_GEOM_MULTI_FN_TYPE, (void (*)()) obj_points, dom);
}

/* destroy MPI related data */
static void destroy_mpi (DOM *dom)
{
  free (dom->dbd);

  stats_destroy (dom);

  Zoltan_Destroy (&dom->zol);
}
#endif

/* go over contact points and remove those whose corresponding
 * areas are much smaller than those of other points related to
 * objects directly topologically adjacent in their shape definitions */
static void sparsify_contacts (DOM *dom)
{
  double threshold = dom->threshold,
	 minarea = dom->minarea,
	 margin = 2.0 * dom->mindist, d [3], c;
  SET *del, *itm;
  CON *con, *adj;
  MEM mem;
  int n;

  MEM_Init (&mem, sizeof (SET), SETBLK);

  for (del = NULL, con = dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT && (con->state & CON_NEW)) /* walk over all new contacts */
    {
      SET *set [2] = {con->master->con, con->slave->con};

      if (con->area < minarea) /* simple criterion first */
      {
	con->state |= CON_DONE;
	SET_Insert (&mem, &del, con, NULL); /* schedule for deletion */
      }
      else /* threshold based adjacency test next */
      {
	for (n = 0; n < 2; n ++) for (itm = SET_First (set [n]); itm; itm = SET_Next (itm))
	{
	  adj = itm->data;

	  if (adj == con || adj->kind != CONTACT || (con->state & CON_DONE)) continue;

	  if (con->area < threshold * adj->area) /* check whether the area of the diagonal element is too small (this test is cheaper => let it go first) */
	  {
	    if (con->master == adj->master && con->slave == adj->slave) /* identify contacts pair sharing the same pairs of bodies */
	    {
	      if (gobj_adjacent (GOBJ_Pair_Code_Ext (mkind(con), mkind(adj)), mgobj(con), mgobj (adj))) /* check whether the geometric objects are topologically adjacent */
	      {
	        con->state |= CON_DONE;
		SET_Insert (&mem, &del, con, NULL); /* if so schedule the current contact for deletion */
	      }
	    }
	    else if (con->master == adj->slave && con->slave == adj->master)
	    {
	      if (gobj_adjacent (GOBJ_Pair_Code_Ext (mkind(con), skind(adj)), mgobj(con), sgobj(adj)))
	      {
	        con->state |= CON_DONE;
		SET_Insert (&mem, &del, con, NULL);
	      }
	    }
	    else if (con->slave == adj->master && con->master == adj->slave)
	    {
	      if (gobj_adjacent (GOBJ_Pair_Code_Ext (skind(con), mkind(adj)), sgobj(con), mgobj(adj)))
	      {
	        con->state |= CON_DONE;
		SET_Insert (&mem, &del, con, NULL);
	      }
	    }
	    else if (con->slave == adj->slave && con->master == adj->master)
	    {
	      if (gobj_adjacent (GOBJ_Pair_Code_Ext (skind(con), skind(adj)), sgobj(con), sgobj(adj)))
	      {
	        con->state |= CON_DONE;
		SET_Insert (&mem, &del, con, NULL);
	      }
	    }
	  }
	  else
	  {
	    SUB (con->point, adj->point, d);

	    MAXABS (d, c);

	    if (c < margin && con->id < adj->id) /* eliminate duplicated contact points (compare ids not to eliminate both) */
	    {
	      con->state |= CON_DONE;
	      SET_Insert (&mem, &del, con, NULL);
	    }
	  }
	}
      }
    }
  }

  /* now remove unwanted contacts */
  for (itm = SET_First (del), n = 0; itm; itm = SET_Next (itm), n ++)
  {
    con = itm->data;
#if MPI
    ext_to_remove (dom, con); /* schedule remote deletion of external constraints */
#endif
    DOM_Remove_Constraint (dom, con); /* remove from the domain */
  }

  dom->nspa = n; /* record the number of sparsified contacts */

  /* clean up */
  MEM_Release (&mem);
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
  case GLUE: return "GLUE";
  }

  return NULL;
}

/* create a domain */
DOM* DOM_Create (AABB *aabb, SPSET *sps, short dynamic, double step)
{
  DOM *dom;

  ERRMEM (dom = MEM_CALLOC (sizeof (DOM)));
  dom->aabb = aabb;
  aabb->dom = dom;
  dom->sps = sps;
  dom->dynamic = (dynamic == 1 ? 1 : 0);
  dom->step = step;
  dom->time = 0.0;

  MEM_Init (&dom->conmem, sizeof (CON), CONBLK);
  MEM_Init (&dom->mapmem, sizeof (MAP), MAPBLK);
  MEM_Init (&dom->setmem, sizeof (SET), SETBLK);
  MEM_Init (&dom->sgpmem, sizeof (SGP), CONBLK);
  MEM_Init (&dom->excmem, sizeof (int [2]), SETBLK);
  dom->bid = 1;
  dom->lab = NULL;
  dom->idb = NULL;
  dom->bod = NULL;
  dom->nbod = 0;
  dom->newb = NULL;
  dom->allbodies = NULL;
  dom->allbodiesread = 0;
  dom->sparecid = NULL;
  dom->excluded = NULL;
  dom->cid = 1;
  dom->idc= NULL;
  dom->con = NULL;
  dom->ncon = 0;
  dom->prev = dom->next = NULL;
  dom->flags = 0;
  dom->threshold = 0.01;
  dom->minarea = 0.0;
  dom->mindist = GEOMETRIC_EPSILON;
  dom->depth = -DBL_MAX;
  ERRMEM (dom->ldy = LOCDYN_Create (dom));

  dom->gravity [0] = NULL;
  dom->gravity [1] = NULL;
  dom->gravity [2] = NULL;

  SET (dom->extents, -DBL_MAX);
  SET (dom->extents + 3, DBL_MAX);

  dom->aabb_data = aabb_create_data ();

  dom->verbose = 0;

#if MPI
  create_mpi (dom);
#endif

  return dom;
}

/* insert a body into the domain */
void DOM_Insert_Body (DOM *dom, BODY *bod)
{
  /* take next id */
  bod->id = dom->bid ++; 

  /* make sure we do not run out of ids */
  ASSERT (dom->bid < UINT_MAX, ERR_DOM_TOO_MANY_BODIES);

  /* assign domain */
  bod->dom = dom;

  /* insert into the set of all created bodies */
  MAP_Insert (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, bod, NULL);

#if MPI
  /* insert every 'rank' body into this domain */
  if (dom->insertbodymode == ALWAYS ||
     (dom->insertbodymode == EVERYNCPU &&
      bod->id % (unsigned) dom->ncpu == (unsigned) dom->rank))
  {
    /* mark as parent */
    bod->flags |= BODY_PARENT;
#endif
    /* insert into overlap engine */
    AABB_Insert_Body (dom->aabb, bod);

    /* insert into label based map */
    if (bod->label) MAP_Insert (&dom->mapmem, &dom->lab, bod->label, bod, (MAP_Compare)strcmp);

    /* insert into id based map */
    MAP_Insert (&dom->mapmem, &dom->idb, (void*) (long) bod->id, bod, NULL);

    /* insert into list */
    bod->next = dom->bod;
    if (dom->bod) dom->bod->prev = bod;
    dom->bod = bod;

    /* increment */
    dom->nbod ++;

    /* schedule body insertion in the output */
    if (dom->time > 0) SET_Insert (&dom->setmem, &dom->newb, bod, NULL);

    /* detailed stats */
    switch (bod->kind)
    {
    case OBS: dom->nobs ++; break;
    case RIG: dom->nrig ++; break;
    case PRB: dom->nprb ++; break;
    case FEM: dom->nfem ++; break;
    }
    dom->dofs += bod->dofs;
#if MPI
  }
#if LOCAL_BODIES
  else
  {
    /* do not store this body */
    MAP_Delete (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, NULL);
    BODY_Destroy (bod);
  }
#endif
#endif
}

/* remove a body from the domain */
void DOM_Remove_Body (DOM *dom, BODY *bod)
{
  /* remove from overlap engine */
  AABB_Delete_Body (dom->aabb, bod);

  SET *con = NULL, *item;

  /* DOM_Remove_Constraint will try to remove the constraint
     from body constraints set, which is not nice if we try
     to iterate over the set at the same time => make a copy */
  for (item = SET_First (bod->con); item; item = SET_Next (item))
    SET_Insert (&dom->setmem, &con, item->data, NULL);
 
  /* remove all body related constraints */
  for (item = SET_First (con); item; item = SET_Next (item))
    DOM_Remove_Constraint (dom, item->data);

  /* free constraint set */
  SET_Free (&dom->setmem, &con);

#if MPI
  if (bod->flags & BODY_PARENT) /* only parent bodies are in the list and label/id maps */
  {
#endif
    /* delete from label based map */
    if (bod->label) MAP_Delete (&dom->mapmem, &dom->lab, bod->label, (MAP_Compare)strcmp);

    /* delete from id based map */
    MAP_Delete (&dom->mapmem, &dom->idb, (void*) (long) bod->id, NULL);

    /* remove from list */
    if (bod->prev) bod->prev->next = bod->next;
    else dom->bod = bod->next;
    if (bod->next) bod->next->prev = bod->prev;

    /* decrement */
    dom->nbod --;

    /* detailed stats */
    switch (bod->kind)
    {
    case OBS: dom->nobs --; break;
    case RIG: dom->nrig --; break;
    case PRB: dom->nprb --; break;
    case FEM: dom->nfem --; break;
    }
    dom->dofs -= bod->dofs;

#if MPI
  }

  /* remove from the domain childeren set if needed */
  if (bod->flags & BODY_CHILD) SET_Delete (&dom->setmem, &dom->children, bod, NULL);

  /* free children set */
  SET_Free (&dom->setmem, &bod->children);

  /* free body id => union of spare ids from all ranks is created during balancing */
  SET_Insert (&dom->setmem, &dom->sparebid, (void*) (long) bod->id, NULL); /* when called from (***) this will do nothing */
#endif

  /* delete from the set of all created bodies */
  MAP_Delete (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, NULL);

  /* make sure that it is not in the new bodies set (could have been
   * inserted and then deleted in the course of fragmentation and cracking) */
  if (dom->time > 0) SET_Delete (&dom->setmem, &dom->newb, bod, NULL);
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
  con = insert (dom, bod, NULL, sgp, NULL, FIXPNT);
  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  IDENTITY (con->base);

  /* insert into local dynamics */
  con->dia = LOCDYN_Insert (dom->ldy, con, bod, NULL);

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
  con = insert (dom, bod, NULL, sgp, NULL, FIXDIR);
  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  localbase (dir, con->base);

  /* insert into local dynamics */
  con->dia = LOCDYN_Insert (dom->ldy, con, bod, NULL);

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
  con = insert (dom, bod, NULL, sgp, NULL, VELODIR);
  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  localbase (dir, con->base);
  con->tms = vel;

  /* insert into local dynamics */
  con->dia = LOCDYN_Insert (dom->ldy, con, bod, NULL);

  return con;
}

/* insert rigid link constraint between two (referential) points of bodies; if one of the body
 * pointers is NULL then the link acts between the other body and the fixed (spatial) point;
 * if the points coincide then a gluing FIXPNT constraint is inserted instead */
CON* DOM_Put_Rigid_Link (DOM *dom, BODY *master, BODY *slave, double *mpnt, double *spnt)
{
  double v [3], d, *tmp;
  SGP *msgp, *ssgp;
  int m, s;
  CON *con;

  if (!master)
  {
    master = slave;
    slave = NULL;
    tmp = mpnt;
    mpnt = spnt;
    spnt = tmp;
  }

  ASSERT_DEBUG (master, "At least one body pointer must not be NULL");

  if ((m = SHAPE_Sgp (master->sgp, master->nsgp, mpnt)) < 0) return NULL;
  msgp = &master->sgp [m];

  if (slave)
  { 
    if ((s = SHAPE_Sgp (slave->sgp, slave->nsgp, spnt)) < 0) return NULL;
    ssgp = &slave->sgp [s];
  }
  else ssgp = NULL;

  SUB (mpnt, spnt, v);
  d = LEN (v);
  
  if (d < GEOMETRIC_EPSILON) /* glue points */
  {
    con = insert (dom, master, slave, msgp, ssgp, FIXPNT);
    COPY (mpnt, con->point);
    COPY (mpnt, con->mpnt);
    COPY (spnt, con->spnt);
    IDENTITY (con->base);
  }
  else
  {
    con = insert (dom, master, slave, msgp, ssgp, RIGLNK);
    COPY (mpnt, con->point);
    COPY (mpnt, con->mpnt);
    COPY (spnt, con->spnt);
    RIGLNK_LEN (con->Z) = d; /* initial distance */
    update_riglnk (dom, con); /* initial update */
  }
  
  /* insert into local dynamics */
  con->dia = LOCDYN_Insert (dom->ldy, con, master, slave);

  return con;
}

/* insert gluging constraint between nodes of regular FEM bodies */
CON* DOM_Glue_Nodes (DOM *dom, BODY *master, BODY *slave, int mnode, int snode)
{
  MESH *mmsh, *smsh;
  double *mpnt, *spnt;
  SGP *msgp, *ssgp;
  CON *con;

  ASSERT_DEBUG (master && slave, "Both body pointer must be valid");

  msgp = SGP_from_index (dom, master, -mnode-1);
  ssgp = SGP_from_index (dom, slave, -snode-1);
  mmsh = msgp->shp->data;
  smsh = ssgp->shp->data;
  mpnt = &mmsh->ref_nodes [mnode][0];
  spnt = &smsh->ref_nodes [snode][0];

  con = insert (dom, master, slave, msgp, ssgp, GLUE);
  COPY (mpnt, con->point);
  COPY (mpnt, con->mpnt);
  COPY (spnt, con->spnt);
  IDENTITY (con->base);

  /* insert into local dynamics */
  con->dia = LOCDYN_Insert (dom->ldy, con, master, slave);

  return con;
}

/* remove a constraint from the domain */
void DOM_Remove_Constraint (DOM *dom, CON *con)
{
  long n = con->msgp - con->master->sgp;
  if (n < 0 || n >= con->master->nsgp) MEM_Free (&dom->sgpmem, con->msgp);
  if (con->slave)
  {
    n = con->ssgp - con->slave->sgp;
    if (n < 0 || n >= con->slave->nsgp) MEM_Free (&dom->sgpmem, con->ssgp);
  }

#if DEBUG
  ASSERT_DEBUG (SET_Contains (con->master->con, con, CONCMP), "Constraint %s with id %d not present in body list", CON_Kind (con), con->id);
  ASSERT_DEBUG (!con->slave || (con->slave && SET_Contains (con->slave->con, con, CONCMP)), "Constraint %s with id %d not present in body list", CON_Kind (con), con->id);
#endif

  /* remove from the body constraint adjacency  */
  SET_Delete (&dom->setmem, &con->master->con, con, CONCMP);
  if (con->slave) SET_Delete (&dom->setmem, &con->slave->con, con, CONCMP);

#if DEBUG
  ASSERT_DEBUG (!SET_Contains (con->master->con, con, CONCMP), "Failed to delete constraint %s with id %d from body list", CON_Kind (con), con->id);
  ASSERT_DEBUG (!con->slave || (con->slave && !SET_Contains (con->slave->con, con, CONCMP)), "Failed to delete constraint %s with id %d from body list", CON_Kind (con), con->id);
#endif

#if MPI
  if (con->state & CON_EXTERNAL)
  {
    /* remove from map */
    MAP_Delete (&dom->mapmem, &dom->conext, (void*) (long) con->id, NULL);

    /* destroy passed data */
    MEM_Free (&dom->conmem, con);
  }
  else
  {
#endif
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
    /* free external ranks */
    SET_Free (&dom->setmem, &con->ext);

    /* free constraint id if possible */
    if (!(con->state & CON_IDLOCK))
#endif
    SET_Insert (&dom->setmem, &dom->sparecid, (void*) (long) con->id, NULL);

    if (con->kind == CONTACT) SURFACE_MATERIAL_Destroy_State (&con->mat); /* free contact material state */
    /* free velocity constraint time history */
    else if (con->kind == VELODIR) TMS_Destroy (con->tms);

    /* destroy passed data */
    MEM_Free (&dom->conmem, con);
#if MPI
  }
#endif
}

/* transfer constraint from the source to the destination body */
void DOM_Transfer_Constraint (DOM *dom, CON *con, BODY *src, BODY *dst)
{
  double point [3];
  int n;

  if (con->kind == CONTACT) return; /* let contacts be re-detected => left in the body will get deleted */

#if MPI
  if (con->state & CON_EXTERNAL) return; /* do not transfer external constraints => let them get deleted with bodies */
#endif

  LOCDYN_Remove (dom->ldy, con->dia);

  if (con->kind == RIGLNK && src == con->slave)
  {
    double *z = RIGLNK_VEC (con->Z);
    SUB (con->point, z, point);
  }
  else
  {
    COPY (con->point, point);
  }

  n = SHAPE_Closest_Sgp (dst->sgp, dst->nsgp, point, NULL);

  if (con->master == src)
  {
    SET_Delete (&dom->setmem, &con->master->con, con, CONCMP);
    con->msgp = &dst->sgp [n];
    con->master = dst;
    SET_Insert (&dom->setmem, &con->master->con, con, CONCMP);
  }
  else
  {
    ASSERT_DEBUG (con->slave == src, "Inconsistent constraint structure: invalid slave pointer");
    SET_Delete (&dom->setmem, &con->slave->con, con, CONCMP);
    con->ssgp = &dst->sgp [n];
    con->slave = dst;
    SET_Insert (&dom->setmem, &con->slave->con, con, CONCMP);
  }

  con->dia = LOCDYN_Insert (dom->ldy, con, con->master, con->slave);

#if MPI
  SET_Free (&dom->setmem, &con->ext);

  WARNING (con->kind != RIGLNK, "Rigid link constraint has not been tested with fragmenting bodies in parallel.\n"
                                "Errorneous results might be produced due to an incorrect handling of it.\n");

  /* TODO/FIXME: parallel rigid link constraint */
#endif
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
  CON *con, *next;
  TIMING timing;
  BOXALG alg;
  BODY *bod;

#if MPI
#if LOCAL_BODIES
  if (dom->time == 0.0) domain_balancing (dom); /* initially balance bodies */
#endif

  if (dom->rank == 0)
#endif
  if (dom->verbose) printf ("DOMAIN ... "), fflush (stdout);

  SOLFEC_Timer_Start (dom->solfec, "TIMINT");

  /* time and step */
  time = dom->time;
  step = dom->step;

  /* initialize bodies */
  if (dom->dynamic > 0)
  {
    for (bod = dom->bod; bod; bod = bod->next)
    {
      BODY_Dynamic_Init (bod); /* integration scheme is set externally */

      double h = BODY_Dynamic_Critical_Step (bod);

      if (h < step) step = h;
    }
  }
  else
  {
    for (bod = dom->bod; bod; bod = bod->next) BODY_Static_Init (bod);
  }

#if MPI
  dom->step = step = PUT_double_min (step);

  if (dom->rank == 0)
#else
  dom->step = step;
#endif
  if (dom->verbose) printf (" (STEP: %.3g) ", step), fflush (stdout);

  /* begin time integration */
  if (dom->dynamic)
    for (bod = dom->bod; bod; bod = bod->next)
      BODY_Dynamic_Step_Begin (bod, time, step);
  else
    for (bod = dom->bod; bod; bod = bod->next)
      BODY_Static_Step_Begin (bod, time, step);

  SOLFEC_Timer_End (dom->solfec, "TIMINT");

#if MPI
  SOLFEC_Timer_Start (dom->solfec, "PARBAL");

  update_children (dom); /* children need to be updated before old constraints (whose update depends on children) */

  SOLFEC_Timer_End (dom->solfec, "PARBAL");
#endif

  SOLFEC_Timer_Start (dom->solfec, "CONUPD");

  /* update old constraints */
  for (con = dom->con; con; con = next)
  {
    next = con->next; /* contact update can delete the current iterate */

    switch (con->kind)
    {
      case CONTACT: update_contact (dom, con); break;
      case FIXPNT:  update_fixpnt  (dom, con); break;
      case FIXDIR:  update_fixdir  (dom, con); break;
      case VELODIR: update_velodir (dom, con); break;
      case RIGLNK:  update_riglnk  (dom, con); break;
      case GLUE:    update_glue (dom, con); break;
    }
  }

#if MPI
  /* external con->point coordinates are need to be updated before the update of body extents;
   * this way slave bodies suitably update their extents and maintain children on the constraint owner processor */
  update_external_constraint_points (dom);
#endif

  /* update body extents after constraints update so that constraint points can be incorporated if needed */
  for (bod = dom->bod; bod; bod = bod->next) BODY_Update_Extents (bod);

  SOLFEC_Timer_End (dom->solfec, "CONUPD");

#if MPI
  SOLFEC_Timer_Start (dom->solfec, "PARBAL");

  domain_balancing (dom); /* migrate bodies (parents and children) and constraints */

  SOLFEC_Timer_End (dom->solfec, "PARBAL");
#endif

  /* detect contacts */
  alg = aabb_algorithm (dom);

#if MPI
  if (dom->rank == 0)
#endif
  if (dom->verbose) printf ("CONDET (%s) ... ", AABB_Algorithm_Name (alg)), fflush (stdout);
  
  SOLFEC_Timer_Start (dom->solfec, "CONDET");

  timerstart (&timing);

  AABB_Update (dom->aabb, alg, dom, (BOX_Overlap_Create) overlap_create);

  aabb_timing (dom, timerend (&timing));

  SOLFEC_Timer_End (dom->solfec, "CONDET");

#if MPI
  SOLFEC_Timer_Start (dom->solfec, "PARBAL");

  domain_gluing_begin (dom); /* migrate new external contacts */

  SOLFEC_Timer_End (dom->solfec, "PARBAL");
#endif

  SOLFEC_Timer_Start (dom->solfec, "CONDET");

  sparsify_contacts (dom); /* once all contacts have been collected sparsify them (aimed at greater consitency of serial and parallel runs) */

  SOLFEC_Timer_End (dom->solfec, "CONDET");

#if MPI
  SOLFEC_Timer_Start (dom->solfec, "PARBAL");

  domain_gluing_end (dom); /* migrate external constraint deletions after sparsification */

  SOLFEC_Timer_End (dom->solfec, "PARBAL");
#else
  ASSERT (!(dom->flags & DOM_DEPTH_VIOLATED), ERR_DOM_DEPTH);
#endif

  SOLFEC_Timer_Start (dom->solfec, "CONUPD");

  for (con = dom->con; con; con = con->next) /* update new constraints */
  {
    if (con->state & CON_NEW)
    {
      if (con->kind == CONTACT) /* new contacts are inserted into LOCDYN only after sparsification */
      {
	con->dia = LOCDYN_Insert (dom->ldy, con, con->master, con->slave); /* insert into local dynamics */
      }
      con->state &= ~CON_NEW; /* invalidate newness */
    }
  }

  SOLFEC_Timer_End (dom->solfec, "CONUPD");

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
  SOLFEC_Timer_Start (dom->solfec, "PARBAL");

  /* update external R, U, V after solution has completed */
  update_external_RUV (dom);

  SOLFEC_Timer_End (dom->solfec, "PARBAL");
#endif

  SOLFEC_Timer_Start (dom->solfec, "TIMINT");

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

  Propagate_Cracks (dom); /* do cracking */

#if MPI
  manage_bodies (dom); /* delete unwanted and insert pending bodies */
#endif

  SOLFEC_Timer_End (dom->solfec, "TIMINT");
}

#if MPI
/* send boundary reactions to their external receivers;
 * if 'normal' is > 0 only normal components are sent */
void DOM_Update_External_Reactions (DOM *dom, short normal)
{
  COMOBJ *send, *recv;
  int i, nrecv;

  ERRMEM (send = malloc (sizeof (COMOBJ [dom->ncpu])));

  for (i = 0; i < dom->ncpu; i ++)
  {
    send [i].o = dom->dbd [i].ext;
    send [i].rank = i;
  }

  if (normal > 0)
  {
    dom->bytes += COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)pack_normal_reactions, dom,
      (OBJ_Unpack)unpack_normal_reactions, send, dom->ncpu, &recv, &nrecv);
  }
  else
  {
    dom->bytes += COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)pack_reactions, dom,
      (OBJ_Unpack)unpack_reactions, send, dom->ncpu, &recv, &nrecv);
  }

  free (send);
  free (recv);
}

/* schedule parallel insertion of a constraint (to be called on all processors) */
int DOM_Pending_Constraint (DOM *dom, short kind, BODY *master, BODY *slave,
    double *mpnt, double *spnt, double *dir, TMS *val, int mnode, int snode)
{
  PNDCON *pnd;

  ERRMEM (pnd = MEM_CALLOC (sizeof (PNDCON)));
  pnd->kind = kind;
  pnd->master = master;
  pnd->slave = slave;
  if (mpnt) COPY (mpnt, pnd->mpnt);
  if (spnt) COPY (spnt, pnd->spnt);
  if (dir) COPY (dir, pnd->dir);
  pnd->val = val;
  pnd->mnode = mnode;
  pnd->snode = snode;

  if (master->kind == FEM && !master->msh)
  {
    if (mnode >= 0) pnd->mele = MESH_Element_With_Node (master->shape->data, mnode);
    else pnd->mele = MESH_Element_Containing_Point (master->shape->data, mpnt, 1);

    if (!pnd->mele) return 0;
  }
  else if (SHAPE_Sgp (master->sgp, master->nsgp, mpnt) < 0) return 0;

  if (slave)
  {
    if (slave->kind == FEM && !slave->msh)
    {
      if (snode >= 0) pnd->sele = MESH_Element_With_Node (slave->shape->data, snode);
      else pnd->sele = MESH_Element_Containing_Point (slave->shape->data, spnt, 1);

      if (!pnd->sele) return 0;
    }
    else if (SHAPE_Sgp (slave->sgp, slave->nsgp, spnt) < 0) return 0;
  }

  SET_Insert (&dom->setmem, &dom->pendingcons, pnd, NULL); /* they will be inserted or deleted during load balancing */

  return 1;
}

/* schedule ASAP insertion of a body in parallel */
void DOM_Pending_Body (DOM *dom, BODY *bod)
{
  SET_Insert (&dom->setmem, &dom->pendingbods, bod, NULL);
}
#endif

/* write domain state */
void DOM_Write_State (DOM *dom, PBF *bf)
{
  dom_write_state (dom, bf);
}

/* read domain state */
void DOM_Read_State (DOM *dom, PBF *bf)
{
  dom_read_state (dom, bf);
}

/* read state of an individual body */
int DOM_Read_Body (DOM *dom, PBF *bf, BODY *bod)
{
  return dom_read_body (dom, bf, bod);
}

/* read state of an individual constraint */
int DOM_Read_Constraint (DOM *dom, PBF *bf, CON *con)
{
  return dom_read_constraint (dom, bf, con);
}

/* exclude contact between a pair of surfaces */
void DOM_Exclude_Contact (DOM *dom, int surf1, int surf2)
{
  int *pair;

  ERRMEM (pair = MEM_Alloc (&dom->excmem));
  if (surf1 <= surf2)
  {
    pair [0] = surf1;
    pair [1] = surf2;
  }
  else
  {
    pair [0] = surf2;
    pair [1] = surf1;
  }
  SET_Insert (&dom->setmem, &dom->excluded, pair, (SET_Compare)pair_compare);
}

/* release memory */
void DOM_Destroy (DOM *dom)
{
  CON *con;
  MAP *item;
 
#if MPI
  destroy_mpi (dom);
#endif

  for (item = MAP_First (dom->allbodies); item; item = MAP_Next (item))
  {
    BODY_Destroy (item->data);
  }

  for (con = dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT) SURFACE_MATERIAL_Destroy_State (&con->mat);
    else if (con->kind == VELODIR && con->tms) TMS_Destroy (con->tms);
  }

  LOCDYN_Destroy (dom->ldy);

  MEM_Release (&dom->conmem);
  MEM_Release (&dom->setmem);
  MEM_Release (&dom->mapmem);
  MEM_Release (&dom->sgpmem);
  MEM_Release (&dom->excmem);

  if (dom->gravity [0]) TMS_Destroy (dom->gravity [0]);
  if (dom->gravity [1]) TMS_Destroy (dom->gravity [1]);
  if (dom->gravity [2]) TMS_Destroy (dom->gravity [2]);

  aabb_destroy_data (dom->aabb_data);

  free (dom);
}

/* export MBFCP definition */
void DOM_2_MBFCP (DOM *dom, FILE *out)
{
  BODY *bod;
  CON *con;
  int n;

  if (dom->gravity [0])
  {
    fprintf (out, "GRAVITY:\n");
    TMS_2_MBFCP (dom->gravity [0], out);
    TMS_2_MBFCP (dom->gravity [1], out);
    TMS_2_MBFCP (dom->gravity [2], out);
    fprintf (out, "\n");
  }

  fprintf (out, "BODIES:\t%d\n\n", dom->nbod);

  for (bod = dom->bod; bod; bod = bod->next)
  {
    BODY_2_MBFCP (bod, out);
  }

  for (con = dom->con, n = 0; con; con = con->next)
  {
    if (con->kind != CONTACT && con->kind != GLUE) n ++;
  }

  fprintf (out, "CONSTRAINTS:\t%d\n\n", n);

  for (con = dom->con; con; con = con->next)
  {
    fprintf (out, "ID:\t%d\n", con->id);

    switch (con->kind)
    {
    case FIXPNT:
      fprintf (out, "KIND:\tFIXPNT\n");
      fprintf (out, "BODY:\t%d\n", con->master->id);
      fprintf (out, "POINT:\t%g  %g  %g\n", con->mpnt [0], con->mpnt [1], con->mpnt [2]);
      break;
    case FIXDIR:
      fprintf (out, "KIND:\tFIXDIR\n");
      fprintf (out, "BODY:\t%d\n", con->master->id);
      fprintf (out, "POINT:\t%g  %g  %g\n", con->mpnt [0], con->mpnt [1], con->mpnt [2]);
      fprintf (out, "DIRECTION:\t%g  %g  %g\n", con->base [6], con->base [7], con->base [8]);
      break;
    case VELODIR:
      fprintf (out, "KIND:\tVELODIR\n");
      fprintf (out, "BODY:\t%d\n", con->master->id);
      fprintf (out, "POINT:\t%g  %g  %g\n", con->mpnt [0], con->mpnt [1], con->mpnt [2]);
      fprintf (out, "DIRECTION:\t%g  %g  %g\n", con->base [6], con->base [7], con->base [8]);
      TMS_2_MBFCP (con->tms, out);
      break;
    case RIGLNK:
      fprintf (out, "KIND:\tRIGLNK\n");
      fprintf (out, "BODY1:\t%d\n", con->master->id);
      if (con->slave) fprintf (out, "BODY2:\t%d\n", con->slave->id);
      else fprintf (out, "BODY2:\t%s\n", "NONE");
      fprintf (out, "POINT1:\t%g  %g  %g\n", con->mpnt [0], con->mpnt [1], con->mpnt [2]);
      fprintf (out, "POINT2:\t%g  %g  %g\n", con->spnt [0], con->spnt [1], con->spnt [2]);
      break;
    default:
      break;
    }

    fprintf (out, "\n");
  }
}
