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
#include "fra.h"
#include "psc.h"
#include "lng.h"

#if MPI
#include "put.h"
#include "com.h"
#endif

#if OMP
#include <omp.h>
#include "ompu.h"
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

  ERRMEM (data = malloc (sizeof (AABB_DATA)));

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

  /* ensure that master pointer was passed;
   * tollerate NULL slave pointers, as this indicates a single body constraint;
   * note that msgp may be NULL for dummy bodies created in get_body function */
  ASSERT_DEBUG (master, "At least master body pointers must be passed");

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
  int state, spair [2], pair [2], ntri;
  SURFACE_MATERIAL *mat;
  short paircode;
  CON *con;
  TRI *tri;

  if (contact_exists (one, two)) return;

  state = gobjcontact (
    CONTACT_DETECT, GOBJ_Pair_Code (one, two),
    one->sgp->shp, one->sgp->gobj,
    two->sgp->shp, two->sgp->gobj,
    onepnt, twopnt, normal,
    &gap, &area, spair, &tri, &ntri);

  if (state)
  {
#if OMP
    omp_set_lock (&dom->lock);
#endif

    ASSERT_DEBUG (gap <= 0, "A contact with positive gap (%g) was detected which indicates a bug in goc.c", gap);

    if (gap <= dom->depth) dom->flags |= DOM_DEPTH_VIOLATED;

    /* set surface pair data if there was a contact */
    mat = SPSET_Find (dom->sps, spair [0], spair [1]);

    if (dom->excluded)
    {
      if (spair [0] <= spair [1]) { pair [0] = spair [0]; pair [1] = spair [1]; }
      else { pair [0] = spair [1]; pair [1] = spair [0]; }

      if (SET_Contains (dom->excluded, pair, (SET_Compare) pair_compare))
      {
	free (tri);
	return; /* exluded pair */
      }
    }

    switch (state)
    {
      case 1: /* first body has outward normal => second body is the master */
      {
	paircode = GOBJ_Pair_Code (one, two);
	con = insert_contact (dom, two->body, one->body, two->sgp, one->sgp, twopnt, onepnt, normal, area, gap, mat, paircode);
	con->spair [0] = spair [0];
	con->spair [1] = spair [1];
      }
      break;
      case 2:  /* second body has outward normal => first body is the master */
      {
	paircode = GOBJ_Pair_Code (two, one);
	con = insert_contact (dom, one->body, two->body, one->sgp, two->sgp, onepnt, twopnt, normal, area, gap, mat, paircode);
	con->spair [0] = spair [1];
	con->spair [1] = spair [0];
      }
      break;
    }
#if OMP
    omp_unset_lock (&dom->lock);
#endif
  }

  if (tri) free (tri);
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
  int state, ntri;
  TRI *tri;

  /* current spatial points and normal */
  BODY_Cur_Point (con->master, con->msgp, con->mpnt, mpnt);
  BODY_Cur_Point (con->slave, con->ssgp, con->spnt, spnt);
  COPY (con->base+6, normal);

  /* update contact data => during an update 'master' and 'slave' relation does not change */
  state = gobjcontact (
    CONTACT_UPDATE, con->paircode,
    sshp, sgobj, mshp, mgobj, /* the slave body holds the outward normal */
    spnt, mpnt, normal, &con->gap, /* 'mpnt' and 'spnt' are updated here */
    &con->area, con->spair, &tri, &ntri); /* surface pair might change though */

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

  if (tri) free (tri);
}

/* update fixed point data */
static void update_fixpnt (DOM *dom, CON *con)
{
  BODY_Cur_Point (con->master, con->msgp, con->mpnt, con->point);
}

/* update fixed direction data */
static void update_fixdir (DOM *dom, CON *con)
{
  if (con->slave)
  {
    double n [3];
    
    BODY_Cur_Point (con->slave, con->ssgp, con->spnt, con->point);
    BODY_Ref_Point (con->master, con->msgp, con->point, con->mpnt);

    BODY_Cur_Vector (con->master, con->msgp->gobj, con->mpnt, con->Z, n);
    localbase (n, con->base);
  }
  else BODY_Cur_Point (con->master, con->msgp, con->mpnt, con->point);
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
static void update_spring (DOM *dom, CON *con)
{
  double v [3],
         n [3],
	 m [3],
	 s [3],
	 len,
	 inv;

  BODY_Cur_Point (con->master, con->msgp, con->mpnt, m);
  BODY_Cur_Point (con->slave, con->ssgp, con->spnt, s);

  COPY (m, con->point);
  SUB (m, s, v);

  switch (con->spair[0])
  {
  case SPRING_FOLLOW:
  {
    COPY (v, n);
    len = LEN (v);
    if (len != 0.0)
    {
      inv = 1.0 / len;
      SCALE (n, inv);
    }
    else /* avoid singularity */
    {
      VECTOR (n, 0, 0, 1);
    }
  }
  break;
  case SPRING_FIXED:
  {
    n[0] = con->Z[4];
    n[1] = con->Z[5];
    n[2] = con->Z[6];
    len = DOT (v, n);
  }
  break;
  case SPRING_CONV_MASTER:
  {
    BODY_Cur_Vector (con->master, con->msgp->gobj, con->mpnt, con->Z+4, n);
    len = DOT (v, n);
  }
  break;
  case SPRING_CONV_SLAVE:
  {
    BODY_Cur_Vector (con->slave, con->ssgp->gobj, con->spnt, con->Z+4, n);
    len = DOT (v, n);
  }
  break;
  }

  con->Z[3] = len; /* useful for visualization */
  con->gap = len - con->Z[2];
  localbase (n, con->base);
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
 * semi-positive index indicates regular surface SGP */
static int SGP_index (BODY *bod, SGP *sgp)
{
  long n = sgp - bod->sgp;

  ASSERT_DEBUG (n >= 0 && n < bod->nsgp, "Error in SGP index");

  return n;
}

/* recreate an SGP from an index returned by SGP_index */
static SGP* SGP_from_index (DOM *dom, BODY *bod, int n)
{
  SGP *sgp;

  ASSERT_DEBUG (n >= 0 && n < bod->nsgp, "Error in SGP index");

  sgp = &bod->sgp [n];

  return sgp;
}

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

#if ZOLTAN
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
  
  for (con = dom->con, i = 0; con; con = con->next)
  {
    global_ids [i * num_gid_entries] = con->id;
    obj_wgts [i * wgt_dim] = constraint_weight (con);
    i ++;
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
#endif

/* create dummy body if needed */
static BODY* get_body (DOM *dom, int id)
{
  BODY *bod;

  bod = MAP_Find (dom->allbodies, (void*) (long) id, NULL);

  if (!bod) /* create dummy */
  {
    ERRMEM (bod = MEM_CALLOC (sizeof (BODY)));
    bod->nsgp = INT_MAX;
    bod->id = id;
    bod->dom = dom;
    MAP_Insert (&dom->mapmem, &dom->allbodies, (void*) (long) id, bod, NULL);
  }

  return bod;
}

/* pack body if the child of this body was not previously on given rank */
static void pack_body (BODY *bod, int rank, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  if (!SET_Contains (bod->prevchildren, (void*) (long) rank, NULL))
  {
    pack_int (isize, i, ints, 1);
    BODY_Pack (bod, dsize, d, doubles, isize, i, ints);
  }
  else pack_int (isize, i, ints, 0);
}

/* unpack body in place of another body or a dummy */
static BODY* unpack_body (DOM *dom, int id, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  BODY *bod, *q;
  MAP *item;
  SET *jtem;
  CON *con;
  int flg;

  flg = unpack_int (ipos, i, ints);

  item = MAP_Find_Node (dom->allbodies, (void*) (long) id, NULL);

  if (flg == 0) /* previous child is expected to exist here */
  {
    ASSERT_DEBUG (item, "Invalid body id");
    return item->data;
  }

  bod = BODY_Unpack (dom->solfec, dpos, d, doubles, ipos, i, ints);
  bod->dom = dom;

  if (item == NULL) /* if not in the map simply insert it */
  {
    MAP_Insert (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, bod, NULL);
  }
  else /* already in the map */
  {
    q = item->data;
   
    if (q->nsgp == INT_MAX) /* q is a dummy */
    {
      for (jtem = SET_First (q->con); jtem; jtem = SET_Next (jtem)) /* for every constraint */
      {
	con = jtem->data;

	if (con->master == q) /* if master was the dummy */
	{
	  if (con->slave) SET_Delete (&dom->setmem, &con->slave->con, con, CONCMP); /* remove constraint from slave's set */

	  con->master = bod; /* remap constraint to the new body */
	  con->msgp = SGP_from_index (dom, bod, con->msgp - q->sgp); /* remap SGP into the new body */

	  if (con->slave) SET_Insert (&dom->setmem, &con->slave->con, con, CONCMP); /* insert updated constraint into slave's set */
	}
	else /* slave was the dummy */
	{
	  SET_Delete (&dom->setmem, &con->master->con, con, CONCMP); /* remove from master's set */

	  con->slave = bod; /* remap body */
	  con->ssgp = SGP_from_index (dom, bod, con->ssgp - q->sgp); /* remap SGP */

	  SET_Insert (&dom->setmem, &con->master->con, con, CONCMP); /* insert updated constraint into master's set */
	}

	SET_Insert (&dom->setmem, &bod->con, con, CONCMP); /* insert updated constraint into body set */
      }

      SET_Free (&dom->setmem, &q->con); /* free dummies constraint set */
      free (q); /* free dummy */
      item->data = bod; /* map new body */
    }
    else /* q is a valid body (formerly a child or a parent on this processor) */
    {
      BODY_Destroy (bod); /* destroy unpacked body */
      bod = q; /* return q */
    }
  }

  return bod;
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
    case FIXPNT:
    case FIXDIR:
    case RIGLNK:
    pack_doubles (dsize, d, doubles, con->Z, DOM_Z_SIZE);
    break;
    case SPRING:
    {
      int fid = lngcallback_id (NULL, con->tms);
      ASSERT_TEXT (fid, "failed to obtain SPRING callback function ID");
      pack_int (isize, i, ints, fid);
      pack_int (isize, i, ints, con->spair[0]); /* direction update kind */
      pack_doubles (dsize, d, doubles, con->Z, DOM_Z_SIZE);
    }
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

  master = get_body (dom, mid);
  if (sid) slave = get_body (dom, sid); else slave = NULL;

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
    case FIXPNT:
    case FIXDIR:
    case RIGLNK:
    unpack_doubles (dpos, d, doubles, con->Z, DOM_Z_SIZE);
    break;
    case SPRING:
    {
      int fid = unpack_int (ipos, i, ints); /* callback */
      fid = lngcallback_set (fid, NULL, (void**) &con->tms);
      ASSERT_TEXT (fid, "failed to set SPRING callback based on ID");
      con->spair[0] = unpack_int (ipos, i, ints); /* direction update kind */
      unpack_doubles (dpos, d, doubles, con->Z, DOM_Z_SIZE);
    }
    break;
  }
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

  master = get_body (dom, mid);
  if (sid) slave = get_body (dom, sid); else slave = NULL;

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
static void pack_parent (BODY *bod, int rank, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  DOM *dom;

  /* must be parent */
  ASSERT_DEBUG (bod->flags & BODY_PARENT, "Not a parent");

  /* set domain */
  dom = bod->dom;

  /* pack id */
  pack_int (isize, i, ints, bod->id);

  /* pack state */
  pack_body (bod, rank, dsize, d, doubles, isize, i, ints);

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
  bod = unpack_body (dom, id, dpos, d, doubles, ipos, i, ints);

  /* must be child or dummy */
  ASSERT_DEBUG ((bod->flags & BODY_PARENT) == 0, "Neither child nor dummy");

  /* unpack state */
  BODY_Parent_Unpack (bod, dpos, d, doubles, ipos, i, ints);

  /* if it was a child */
  if (bod->flags & BODY_CHILD)
  {
    /* unmark child */
    bod->flags &= ~BODY_CHILD;

    /* delete from children map */
    MAP_Delete (&dom->mapmem, &dom->children, (void*) (long) bod->id, NULL);
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

#if PSCTEST
  PSC_Test_Body (bod);
#endif
}

/* pack migrating out child body */
static void pack_child (BODY *bod, int rank, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  /* must be an exported or an existing parent */
  ASSERT_DEBUG (((bod->flags & (BODY_PARENT|BODY_CHILD)) == 0 && bod->rank != bod->dom->rank) /* just migrating out parent after being packed (hence unmarked) */
                || (bod->flags & BODY_PARENT), "Not a parent"); /* or an existing parent */

  /* pack id */
  pack_int (isize, i, ints, bod->id);

  /* pack state */
  pack_body (bod, rank, dsize, d, doubles, isize, i, ints);

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
  bod = unpack_body (dom, id, dpos, d, doubles, ipos, i, ints);

  /* must be child or dummy */
  ASSERT_DEBUG ((bod->flags & BODY_PARENT) == 0, "Neither child nor dummy");

  /* unpack state */
  BODY_Child_Unpack (bod, dpos, d, doubles, ipos, i, ints);

  /* if it was a dummy */
  if ((bod->flags & BODY_CHILD) == 0)
  {
    /* insert into children map */
    MAP_Insert (&dom->mapmem, &dom->children, (void*) (long) bod->id, bod, NULL);

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

    SET_Free (&dom->setmem, &bod->prevchildren); /* empty previous children set */
    bod->prevchildren = bod->children;
    bod->children = NULL;

#if ZOLTAN
    Zoltan_LB_Box_Assign (dom->zol, e[0], e[1], e[2], e[3], e[4], e[5], procs, &numprocs);
#else
    numprocs = dynlb_box_assign (dom->lb, e, e+3, procs);
#endif

    for (i = 0; i < numprocs; i ++) /* child migration based on geometrical extents */
    {
      if (bod->rank != procs [i]) /* if this is neither current nor the new body rank */
      {
	SET_Insert (&dom->setmem, &dbd [procs [i]].children, bod, NULL); /* schedule for sending a child */
	SET_Insert (&dom->setmem, &bod->children, (void*) (long) procs [i], NULL); /* extend parent's children set */
      }
    }

    SET *pairedup = MAP_Find (dom->pairedup, (void*) (long) bod->id, NULL);
    for (SET *item = SET_First (pairedup); item; item = SET_Next (item)) /* child migration based on involvement in two-body bilateral constraints */
    {
      MAP *node = MAP_Find_Node (dom->idtorank, item->data, NULL);
      ASSERT_TEXT (node, "Inconsistent DOM->idtorank mapping");
      int pairedup_rank = (int) (long) node->data;
      for (i = 0; i < numprocs; i ++)
      {
	if (procs [i] == pairedup_rank) break; /* if already migrating to pairedup_rank based on gometrical extents */
      }
      if (i == numprocs && bod->rank != pairedup_rank) /* can migrate there */
      {
	SET_Insert (&dom->setmem, &dbd [pairedup_rank].children, bod, NULL); /* schedule for sending a child */
	SET_Insert (&dom->setmem, &bod->children, (void*) (long) pairedup_rank, NULL); /* extend parent's children set */
      }
    }
  }

  free (procs);
}

/* delete migrated out children */
static void children_migration_end (DOM *dom)
{
  SET *delset, *item;
  MAP *jtem;

  delset = NULL;

  for (jtem = MAP_First (dom->children); jtem; jtem = MAP_Next (jtem))
  {
    BODY *bod = jtem->data;

    /* must be a child */
    ASSERT_DEBUG (bod->flags & BODY_CHILD, "Not a child");

    if ((bod->flags & BODY_CHILD_UPDATED) == 0) /* migrated out as it wasn't updated by a parent */
    {
      bod->flags &= ~BODY_CHILD; /* unmark child */
      
      SET_Insert (&dom->setmem, &delset, (void*) (long)bod->id, NULL); /* schedule deletion from dom->children */
    }
    else bod->flags &= ~BODY_CHILD_UPDATED; /* invalidate update flag */
  }

  /* subtract deleted children from domain children map */
  for (item = SET_First (delset); item; item = SET_Next (item))
  {
    MAP_Delete (&dom->mapmem, &dom->children, item->data, NULL);
  }

  SET_Free (&dom->setmem, &delset);
}

/* pack domain balancing data */
static void domain_balancing_pack (DBD *dbd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;

  /* pack exported bodies */
  pack_int (isize, i, ints, SET_Size (dbd->bodies));
  for (item = SET_First (dbd->bodies); item; item = SET_Next (item))
    pack_parent (item->data, dbd->rank, dsize, d, doubles, isize, i, ints);

  /* pack exported children */
  pack_int (isize, i, ints, SET_Size (dbd->children));
  for (item = SET_First (dbd->children); item; item = SET_Next (item))
    pack_child (item->data, dbd->rank, dsize, d, doubles, isize, i, ints);

  /* pack exported constraints */
  pack_int (isize, i, ints, SET_Size (dbd->constraints));
  for (item = SET_First (dbd->constraints); item; item = SET_Next (item))
  {
    CON *con = item->data;
    pack_constraint (con, dsize, d, doubles, isize, i, ints);
  }

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
  PNDCON *pnd;
  SET *item;
  CON *con;

  for (item = SET_First (dom->pendingcons); item; item = SET_Next (item))
  {
    pnd = item->data;

    if (pnd->master == NULL && pnd->mid) pnd->master = MAP_Find (dom->idb, (void*) (long) pnd->mid, NULL);

    if (pnd->master && pnd->master->flags & BODY_PARENT) /* insert only those having parent master */
    {
      if (pnd->slave == NULL && pnd->sid)
      {
	if (!(pnd->slave = MAP_Find (dom->idb, (void*) (long) pnd->sid, NULL)))
	{
	  ASSERT_TEXT (pnd->slave = MAP_Find (dom->children, (void*) (long) pnd->sid, NULL),
				    "Missing slave body for a two-body pending constraint");
	}
      }

      switch (pnd->kind)
      {
      case FIXPNT:
	con = DOM_Fix_Point (dom, pnd->master, pnd->mpnt, pnd->strength);
	break;
      case FIXDIR:
	con = DOM_Fix_Direction (dom, pnd->master, pnd->mpnt, pnd->dir, pnd->slave, pnd->spnt);
	break;
      case VELODIR:
	con = DOM_Set_Velocity (dom, pnd->master, pnd->mpnt, pnd->dir, pnd->val);
	break;
      case RIGLNK:
	con = DOM_Put_Rigid_Link (dom, pnd->master, pnd->slave, pnd->mpnt, pnd->spnt, pnd->strength);
	break;
      case SPRING:
	con = DOM_Put_Spring (dom, pnd->master, pnd->mpnt, pnd->slave, pnd->spnt, pnd->val, pnd->lim, pnd->dir, pnd->update);
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
      case SPRING: update_spring (dom, con); break;
      }
    }

    free (pnd); /* they were MEM_ALLOC-ed */
  }

  /* empty pending constraints set */
  SET_Free (&dom->setmem, &dom->pendingcons);
}

#if 0
/* update global body id to rank mapping */
static void update_bidsets (DOM *dom)
{
  int i, j, n, k [2], size [dom->ncpu][2], count [dom->ncpu], disp [dom->ncpu], *send, *recv;
  MAP *jtem;
  SET *item;

  /* count modifications */
  for (item = SET_First (dom->bidset [dom->rank]), k[0] = 0; item; item = SET_Next (item))
  {
    if (!MAP_Find (dom->allbodies, item->data, NULL)) k[0] ++; /* deletions */
  }

  for (jtem = MAP_First (dom->allbodies), k[1] = 0; jtem; jtem = MAP_Next (jtem))
  {
    if (!SET_Find (dom->bidset [dom->rank], jtem->key, NULL)) k[1] ++; /* insertions */
  }

  MPI_Allgather (k, 2, MPI_INT, (int*)size, 2, MPI_INT, MPI_COMM_WORLD); /* size contains deletion/insertion counts from each rank */

  j = size[dom->rank][0] + size[dom->rank][1];
  for (i = n = 0; i < dom->ncpu; i ++)
  {
    count [i] = size[i][0] + size[i][1];
    disp [i] = n;
    n += count[i];
  }

  ERRMEM (send = malloc ((n+j) * sizeof(int)));
  recv = send + j;

  /* prepare send buffer */
  for (item = SET_First (dom->bidset [dom->rank]), j = 0; item; item = SET_Next (item))
  {
    if (!MAP_Find (dom->allbodies, item->data, NULL)) send [j ++] = (int) (long) item->data; /* deletion ids */
  }

  for (jtem = MAP_First (dom->allbodies); jtem; jtem = MAP_Next (jtem))
  {
    if (!SET_Find (dom->bidset [dom->rank], jtem->key, NULL)) send [j ++] = (int) (long) jtem->key; /* insertion ids */
  }

  MPI_Allgatherv (send, j, MPI_INT, recv, count, disp, MPI_INT, MPI_COMM_WORLD); /* recv contains deletion/insertion ids from each rank */

  for (i = n = 0; i < dom->ncpu; i ++)
  {
    for (j = 0; j < size[i][0]; j ++, n ++) /* for all deletions from rank i */
    {
      SET_Delete (&dom->setmem, &dom->bidset [i], (void*) (long) recv[n], NULL);
    }

    for (j = 0; j < size[i][1]; j ++, n ++) /* for all insertions from rank i */
    {
      SET_Insert (&dom->setmem, &dom->bidset [i], (void*) (long) recv[n], NULL);
    }
  }

  free (send);

#if DEBUG
  for (item = SET_First (dom->bidset [dom->rank]); item; item = SET_Next (item))
  {
    ASSERT_DEBUG (MAP_Find (dom->allbodies, item->data, NULL), "bidset[rank] != dom->allbodies!\n");
  }

  for (jtem = MAP_First (dom->allbodies); jtem; jtem = MAP_Next (jtem))
  {
    ASSERT_DEBUG (SET_Find (dom->bidset [dom->rank], jtem->key, NULL), "bidset[rank] != dom->allbodies!\n");
  }
#endif
}
#endif

/* reset dom->idtorank mapping */
static void reset_idtorank (DOM *dom)
{
  int rank, size, *sendbuf, sendcount, *recvbuf, *recvcounts, *displs, i, recvsize, *ptr, j;
  BODY *bod;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  sendcount = 2 * dom->nbod;
  ERRMEM (recvcounts = malloc(size * sizeof(int)));
  ERRMEM (displs = malloc(size * sizeof(int)));

  MPI_Allgather (&sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

  for (displs[0] = recvsize = i = 0; i < size; i ++)
  {
    if (i) displs[i] = recvsize;
    recvsize += recvcounts[i];
  }

  ERRMEM (recvbuf = malloc(recvsize * sizeof(int)));
  ERRMEM (sendbuf = malloc(sendcount * sizeof(int)));

  for (bod = dom->bod, ptr = sendbuf; bod; bod = bod->next, ptr += 2)
  {
    ptr[0] = bod->id;
    ptr[1] = bod->rank;
  }

  MPI_Allgatherv (sendbuf, sendcount, MPI_INT, recvbuf, recvcounts, displs, MPI_INT, MPI_COMM_WORLD);

  MAP_Free (&dom->mapmem, &dom->idtorank); /* free previous mapping */

  for (i = 0; i < size; i ++)
  {
    ASSERT_TEXT (recvcounts[i] % 2 == 0, "Inconsistent receive count in dom.c:reset_idtorank");

    for (j = 0; j < recvcounts[i]/2; j ++)
    {
      long id = recvbuf[displs[i]+2*j];
      long rank = recvbuf[displs[i]+2*j+1];
      MAP_Insert (&dom->mapmem, &dom->idtorank, (void*) id, (void*) rank, NULL); /* map current ranks */
    }
  }

  free (sendbuf);
  free (recvbuf);
  free (displs);
  free (recvcounts);
}

/* domain balancing */
static void domain_balancing (DOM *dom)
{
#if ZOLTAN
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
  char str [128];
#endif

  COMOBJ *send, *recv;
  int i, rank;
  int nrecv;
  SET *item;
  BODY *bod;
  DBD *dbd;
  CON *con;

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

#if ZOLTAN
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
#if PSCTEST
	  PSC_Write_Body (bod);
#endif
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
#if PSCTEST
	PSC_Write_Body (bod);
#endif
      }
    }

    for (con = dom->con; con; con = con->next)
    {
      switch (con->kind)
      {
      case RIGLNK:
      case SPRING:
	ASSERT_DEBUG (con->master->flags & BODY_PARENT, "Inconsistency: two-body bilateral constraints migrate with parents");
	rank = con->master->rank;  
	break;
      default:
        Zoltan_LB_Point_Assign (dom->zol, con->point, &rank);
	break;
      }
      if (rank != dom->rank) SET_Insert (&dom->setmem, &dbd [rank].constraints, con, NULL);
    }

    dom->rebalanced ++;
  }
#else
  /* processing constraints dominates computational time -->
   * we favour the use of contact points to guide load balancing */
  int npoint = dom->ncon > dom->nbod ? dom->ncon : dom->nbod;
  double *point[3];

  ERRMEM (point[0] = malloc (npoint * sizeof (double)));
  ERRMEM (point[1] = malloc (npoint * sizeof (double)));
  ERRMEM (point[2] = malloc (npoint * sizeof (double)));

  if (dom->ncon > dom->nbod)
  {
    for (con = dom->con, i = 0; con; con = con->next, i ++)
    {
      point[0][i] = con->point[0];
      point[1][i] = con->point[1];
      point[2][i] = con->point[2];
    }
  }
  else
  {
    for (bod = dom->bod, i = 0; bod; bod = bod->next, i ++)
    {
      double *e = bod->extents, v[3];
      MID (e, e+3, v);
      point[0][i] = v[0];
      point[1][i] = v[1];
      point[2][i] = v[2];
    }
  }

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (dom->lb == NULL)
  {
    dom->lb = dynlb_create (0, npoint, point, 0, dom->imbalance_tolerance - 1.0, DYNLB_RCB_TREE);
  } 
  else if (dom->rebalanced % dom->updatefreq == 0)
  {
    dom->lb->epsilon = dom->imbalance_tolerance - 1.0;

    dynlb_update (dom->lb, npoint, point);
  }

  dom->rebalanced ++;

  free (point[0]);
  free (point[1]);
  free (point[2]);

  for (bod = dom->bod; bod; bod = bod->next)
  {
    double *e = bod->extents, v [3];
    MID (e, e+3, v);
    rank = dynlb_point_assign (dom->lb, v);
    if (rank != dom->rank)
    {
      bod->rank = rank;
      SET_Insert (&dom->setmem, &dbd [rank].bodies, bod, NULL);
#if PSCTEST
      PSC_Write_Body (bod);
#endif
    }
  }

  for (con = dom->con; con; con = con->next)
  {
    switch (con->kind)
    {
    case RIGLNK:
    case SPRING:
      ASSERT_DEBUG (con->master->flags & BODY_PARENT, "Inconsistency: two-body bilateral constraints migrate with parents");
      rank = con->master->rank;  
      break;
    default:
      rank = dynlb_point_assign (dom->lb, con->point);
      break;
    }
    if (rank != dom->rank) SET_Insert (&dom->setmem, &dbd [rank].constraints, con, NULL);
  }
#endif

  /* --- */

  ERRMEM (send = malloc (sizeof (COMOBJ [dom->ncpu])));

  for (i = 0; i < dom->ncpu; i ++)
  {
    send [i].rank = i;
    send [i].o = &dbd [i];
  }

  /* reset idtorank map */
  reset_idtorank (dom);

  /* compute chidren migration sets */
  children_migration_begin (dom, dbd);

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

  for (con = dom->con; con; con = con->next) /* insert into local dynamics */
  {
    if (!con->dia) con->dia = LOCDYN_Insert (dom->ldy, con, con->master, con->slave);
  }

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
    }
  }
}

/* pack domain gluing data */
static void domain_gluing_begin_pack (DBD *dbd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;

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
  }
}

/* unpack children udate data => NOTE, that this routine is called in the sequence of ranks */
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

    dom->bid += j; /* account for the global body ID increments */

    ASSERT (dom->bid < UINT_MAX, ERR_DOM_TOO_MANY_BODIES); /* make sure we do not run out of ids */
  }
  else
  {
    for (item = SET_First (dom->pendingbods); item; item = SET_Next (item))
    {
      BODY *bod = item->data;
      dom->insertbodymode = ALWAYS;
      DOM_Insert_Body (dom, bod); /* insert pending bodies */
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
      if (bod->nsgp == INT_MAX) /* remove dummy created in get_body */
      {
	/* DOM_Remove_Constraint will remove the constraint from body constraints set,
	 * which is not nice if we try to iterate over the set at the same time => make a copy */
	SET *copy = SET_Copy (&dom->setmem, bod->con), *jtem;

	/* remove all dummy related constraints */
	for (jtem = SET_First (copy); jtem; jtem = SET_Next (jtem)) DOM_Remove_Constraint (dom, jtem->data);

	/* free copy */
	SET_Free (&dom->setmem, &copy);

	/* remove dummy from all bodies map */
	MAP_Delete (&dom->mapmem, &dom->allbodies, item->data, NULL);

	/* free dummy */
	free (bod);
      }
      else /* regular body */
      {
	DOM_Remove_Body (dom, bod); /* while looping over 'sparebid' => look there (***) */
	BODY_Destroy (bod);
      }
    }
  }

  /* empty body ids set */
  SET_Free (&dom->setmem, &dom->sparebid);

  /* restore body insertion mode */
  dom->insertbodymode = EVERYNCPU;

  /* delete unused bodies with empty constraint sets */
  SET *todel = NULL;
  MAP *jtem;

  for (jtem = MAP_First (dom->allbodies); jtem; jtem = MAP_Next (jtem))
  {
    bod = jtem->data;
    if ((bod->flags & (BODY_PARENT|BODY_CHILD)) == 0 && SET_Size (bod->con) == 0)
    {
      if (!SET_Contains (dom->newb, bod, NULL)) /* after t > 0 insertion bodies are in the 'newb' set until written in dio.c:dom_write_state */
      {                                         /* since writing happens every-skip-steps, bodies in the 'newb' step should be preserved until then */
        SET_Insert (&dom->setmem, &todel, bod, NULL);
      }
    }
  }

  for (item = SET_First (todel); item; item = SET_Next (item))
  {
    bod = item->data;
    MAP_Delete (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, NULL);
    if (bod->nsgp == INT_MAX) free (bod); /* dummy created in get_body */
    else BODY_Destroy (bod); /* regular body */
  }

  SET_Free (&dom->setmem, &todel);
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

  dom->pendingcons = NULL;

  dom->pendingbods = NULL;

  dom->pendingremovebods = NULL;

  dom->pairedup = NULL;

  dom->idtorank = NULL;

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

  dom->imbalance_tolerance = 1.1;
  dom->weight_factor = 1.0;

#if ZOLTAN
  ASSERT (dom->zol = Zoltan_Create (MPI_COMM_WORLD), ERR_ZOLTAN); /* zoltan context domain partitioning */
#else
    dom->lb = NULL;
#endif

#if ZOLTAN
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
#endif
}

/* destroy MPI related data */
static void destroy_mpi (DOM *dom)
{
  free (dom->dbd);

  stats_destroy (dom);

#if ZOLTAN
  Zoltan_Destroy (&dom->zol);
#else
  dynlb_destroy (dom->lb);
#endif
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
    if (con->kind == CONTACT && (con->state & CON_NEW)) /* walk over all new primary contacts */
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
  case SPRING: return "SPRING";
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

#if OMP
  omp_init_lock (&dom->lock);
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
     (dom->insertbodymode == EVERYNCPU && dom->rank == 0)) /* all bodies created on rank 0 */
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
    
    if (dom->time > 0.0) /* initialize body */
    {
      if (dom->dynamic) BODY_Dynamic_Init (bod);
      else BODY_Static_Init (bod);
    }
#if MPI
  }
  else
  {
    /* do not store this body */
    MAP_Delete (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, NULL);
    BODY_Destroy (bod);
  }
#endif
}

/* remove a body from the domain */
void DOM_Remove_Body (DOM *dom, BODY *bod)
{
  /* remove from overlap engine */
  AABB_Delete_Body (dom->aabb, bod);

  /* DOM_Remove_Constraint will remove the constraint from body constraints set,
   * which is not nice if we try to iterate over the set at the same time => make a copy */
  SET *copy = SET_Copy (&dom->setmem, bod->con), *item;

  /* remove all body related constraints */
  for (item = SET_First (copy); item; item = SET_Next (item)) DOM_Remove_Constraint (dom, item->data);

  /* free copy */
  SET_Free (&dom->setmem, &copy);

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

  /* remove from the domain childeren map if needed */
  if (bod->flags & BODY_CHILD) MAP_Delete (&dom->mapmem, &dom->children, (void*) (long)bod->id, NULL);

  /* free children sets */
  SET_Free (&dom->setmem, &bod->children);
  SET_Free (&dom->setmem, &bod->prevchildren);

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
CON* DOM_Fix_Point (DOM *dom, BODY *bod, double *pnt, double strength)
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

  STRENGTH (con->Z) = strength; /* set up strength */

  return con;
}

/* fix a referential point of the body along the spatial direction */
CON* DOM_Fix_Direction (DOM *dom, BODY *bod, double *pnt, double *dir, BODY *bod2, double *pnt2)
{
  SGP *sgp, *sgp2;
  CON *con;
  int n, m;

  if ((n = SHAPE_Sgp (bod->sgp, bod->nsgp, pnt)) < 0) return NULL;
  sgp = &bod->sgp [n];

  if (bod2)
  {
    if ((m = SHAPE_Sgp (bod2->sgp, bod2->nsgp, pnt2)) < 0) return NULL;
    sgp2 = &bod2->sgp [m];
  }
  else sgp2 = NULL;

  con = insert (dom, bod, bod2, sgp, sgp2, FIXDIR);

  COPY (pnt, con->point);
  COPY (pnt, con->mpnt);
  localbase (dir, con->base);

  if (bod2)
  {
    COPY (pnt2, con->point); /* slider direction attached to the slave body */
    COPY (pnt2, con->spnt);
    COPY (dir, con->Z); /* normal will be taken from the master body though */
  }

  /* insert into local dynamics */
  con->dia = LOCDYN_Insert (dom->ldy, con, bod, bod2);

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
CON* DOM_Put_Rigid_Link (DOM *dom, BODY *master, BODY *slave, double *mpnt, double *spnt, double strength)
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

  STRENGTH (con->Z) = strength; /* set up strength */
  
  /* insert into local dynamics */
  con->dia = LOCDYN_Insert (dom->ldy, con, master, slave);

  return con;
}

/* create user spring constraint */
CON* DOM_Put_Spring (DOM *dom, BODY *master, double *mpnt, BODY *slave, double *spnt, void *function, double *lim, double *direction, int update)
{
  SGP *msgp, *ssgp;
  double v[3], w[3];
  int m, s;
  CON *con;

  ASSERT_DEBUG (master && slave, "Both bodies needs to be passed");

  if ((m = SHAPE_Sgp (master->sgp, master->nsgp, mpnt)) < 0) return NULL;
  msgp = &master->sgp [m];

  if ((s = SHAPE_Sgp (slave->sgp, slave->nsgp, spnt)) < 0) return NULL;
  ssgp = &slave->sgp [s];

  con = insert (dom, master, slave, msgp, ssgp, SPRING);
  COPY (mpnt, con->point);
  COPY (mpnt, con->mpnt);
  COPY (spnt, con->spnt);
  con->Z[0] = lim[0];
  con->Z[1] = lim[1];
  SUB (mpnt, spnt, v);
  COPY (direction, w);
  NORMALIZE (w);
  if (update == SPRING_FOLLOW) con->Z[2] = LEN(v); /* initial distance */
  else con->Z[2] = DOT (w, v); /* initial distance */
  /* con->Z[3] will store length of spring */
  con->Z[4] = w[0];
  con->Z[5] = w[1];
  con->Z[6] = w[2];
  con->spair[0] = update; /* store update in spair[0] */
  con->tms = function;
  update_spring (dom, con); /* initial update */

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

/* initialize domain at t == 0.0 */
void DOM_Initialize (DOM *dom)
{
  BODY *bod;

#if MPI
  SOLFEC_Timer_Start (dom->solfec, "PARBAL");

  domain_balancing (dom); /* initially balance bodies */
  
  insert_pending_constraints (dom); /* insert pending constraints */

  SOLFEC_Timer_End (dom->solfec, "PARBAL");
#endif

  SOLFEC_Timer_Start (dom->solfec, "TIMINT");

  /* initialize bodies */
  if (dom->dynamic > 0)
  {
    for (bod = dom->bod; bod; bod = bod->next)
    {
      BODY_Dynamic_Init (bod); /* integration scheme is set externally */
    }
  }
  else
  {
    for (bod = dom->bod; bod; bod = bod->next)
    {
      BODY_Static_Init (bod);
    }
  }

  SOLFEC_Timer_End (dom->solfec, "TIMINT");
}

/* domain update initial half-step => bodies and constraints are
 * updated and the current local dynamic problem is returned */
LOCDYN* DOM_Update_Begin (DOM *dom)
{
  double time, step;
  TIMING timing;
  BOXALG alg;
  BODY *bod;
  CON *con;

#if PSCTEST
    for (bod = dom->bod; bod; bod = bod->next)
    {
#if MPI
      if (dom->rank == 0)
#endif
      if (bod->kind != FEM)
      {
	fprintf (stderr, "WARNING: only FEM bodies are supported in the Parallel self-consitency test mode (PSCTEST).\n");
	fprintf (stderr, "WARNING: the simulation is aborted.\n");
	EXIT (1);
      }
    }
#endif

#if MPI
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
      double h = BODY_Dynamic_Critical_Step (bod);

      if (h < step) step = h;
    }
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
  CON *next;
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
      case SPRING:  update_spring (dom, con); break;
    }
  }

#if MPI
  /* external con->point coordinates need to be updated before the update of body extents;
   * this is way slave bodies suitably update their extents and maintain children on the constraint owner processor */
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
  CON *con, *next;
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

  Fracture_Check (dom);

  Propagate_Cracks (dom); /* do cracking */

#if MPI
  /* remove pending bodies */
  for (item = SET_First (dom->pendingremovebods); item; item = SET_Next (item))
  {
    DOM_Remove_Body (dom, item->data);
    BODY_Destroy (item->data);
  }
  SET_Free (&dom->setmem, &dom->pendingremovebods); /* empty set */

  manage_bodies (dom); /* delete unwanted and insert pending bodies */
#endif

  /* remove failed rigid links */
  for (con = dom->con; con; con = next)
  {
    next = con->next; /* contact update can delete the current iterate */

    if (con->kind == RIGLNK)
    {
      double strength = STRENGTH (con->Z);

      if (strength != DBL_MAX)
      {
	if (con->R[2] > strength) /* tensile strength */
	{
#if MPI
	  ext_to_remove (dom, con); /* schedule remote deletion of external constraints */
#endif
	  DOM_Remove_Constraint (dom, con); /* remove from the domain */
	}
      }
    }
    else if (con->kind == FIXPNT)
    {
      double strength = STRENGTH (con->Z);

      if (strength != DBL_MAX)
      {
	double norm = LEN (con->R);

	if (norm > strength)
	{
#if MPI
	  ext_to_remove (dom, con); /* schedule remote deletion of external constraints */
#endif
	  DOM_Remove_Constraint (dom, con); /* remove from the domain */
	}
      }
    }
  }

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
int DOM_Pending_Constraint (DOM *dom, short kind, BODY *master, unsigned int mid, BODY *slave,
    unsigned int sid, double *mpnt, double *spnt, double *dir, TMS *val, int mnode, int snode,
    double strength, double *lim, int update)
{
  PNDCON *pnd;

  ERRMEM (pnd = MEM_CALLOC (sizeof (PNDCON)));
  pnd->kind = kind;
  pnd->master = master;
  pnd->slave = slave;
  pnd->mid = mid;
  pnd->sid = sid;
  if (mpnt) COPY (mpnt, pnd->mpnt);
  if (spnt) COPY (spnt, pnd->spnt);
  if (dir) COPY (dir, pnd->dir);
  pnd->val = val;
  pnd->mnode = mnode;
  pnd->snode = snode;
  pnd->strength = strength;
  if (lim) COPY2 (lim, pnd->lim);
  pnd->update = update;

  if (master)
  {
    if (master->kind == FEM && !master->msh)
    {
      if (mnode >= 0) pnd->mele = MESH_Element_With_Node (master->shape->data, mnode);
      else pnd->mele = MESH_Element_Containing_Point (master->shape->data, mpnt, 1);

      if (!pnd->mele) return 0;
    }
    else if (SHAPE_Sgp (master->sgp, master->nsgp, mpnt) < 0) return 0;
  }

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

  if (mid && sid)
  {
    MAP *node;

    node = MAP_Find_Node (dom->pairedup, (void*) (long) mid, NULL);
    if (node == NULL)
    {
      node = MAP_Insert (&dom->mapmem, &dom->pairedup, (void*) (long) mid, NULL, NULL);
    }
    SET_Insert (&dom->setmem, (SET**) &node->data, (void*) (long) sid, NULL);

    node = MAP_Find_Node (dom->pairedup, (void*) (long) sid, NULL);
    if (node == NULL)
    {
      node = MAP_Insert (&dom->mapmem, &dom->pairedup, (void*) (long) sid, NULL, NULL);
    }
    SET_Insert (&dom->setmem, (SET**) &node->data, (void*) (long) mid, NULL);
  }

  return 1;
}

/* schedule ASAP insertion of a body in parallel */
void DOM_Pending_Body_Insert (DOM *dom, BODY *bod)
{
  SET_Insert (&dom->setmem, &dom->pendingbods, bod, NULL);
}

/* schedule ASAP removal of a body in parallel (to be called on one processor) */
void DOM_Pending_Body_Remove (DOM *dom, BODY *bod)
{
  SET_Insert (&dom->setmem, &dom->pendingremovebods, bod, NULL);
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

#if OMP
  omp_destroy_lock (&dom->lock);
#endif
 
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
    if (con->kind != CONTACT) n ++;
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
