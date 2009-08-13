/*
 * box.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * axis aligned bounding box overlap detection
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

#include "box.h"
#include "hyb.h"
#include "alg.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "bod.h"
#include "swp.h"
#include "hsh.h"
#include "err.h"

#if MPI
#include "dom.h"
#include "pck.h"
#include "com.h"
#include "tag.h"
#endif

#define SIZE 512 /* mempool size */

/* auxiliary data */
struct auxdata
{
  SET *nobody;
  SET *nogobj;
  MEM *mapmem;
  void *data;
  BOX_Overlap_Create create;
};

/* body pair comparison */
static int bodcmp (OPR *a, OPR *b)
{
  if (a->bod1 < b->bod1) return -1;
  else if (a->bod1 == b->bod1 && a->bod2 < b->bod2) return -1;
  else if (a->bod1 == b->bod1 && a->bod2 == b->bod2) return 0;
  else return 1;
}

/* geomtric object pair comparison */
static int gobjcmp (OPR *a, OPR *b)
{
  int cmp = bodcmp (a, b);

  if (cmp == 0)
  {
    if (a->sgp1 < b->sgp1) return -1;
    else if (a->sgp1 == b->sgp1 && a->sgp2 < b->sgp2) return -1;
    else if (a->sgp1 == b->sgp1 && a->sgp2 == b->sgp2) return 0;
    else return 1;
  }
  else return cmp;
}

/* check for a lack of overlap */
static int no_overlap (double *a, double *b)
{
  if (a[0] > b[3] || b[0] > a[3]) return 1;
  if (a[1] > b[4] || b[1] > a[4]) return 1;
  if (a[2] > b[5] || b[2] > a[5]) return 1;
  return 0;
}

/* local overlap creation callback => filters our unwnated adjacency */
static void* local_create (struct auxdata *aux, BOX *one, BOX *two)
{
  void *user;

  /* check if these are two obstacles => no need to report overlap */
  if (BODY(one->body)->kind == OBS && BODY(two->body)->kind == OBS) return NULL;

  /* boxes already adjacent ? */
  if (MAP_Find (one->adj, two, NULL)) return NULL;

  /* self-contact ? <=> one->body == two->body */
  if (one->body == two->body)
  {
    if (!(((BODY*)one->body)->flags & BODY_DETECT_SELF_CONTACT)) return NULL;
    /* it's much more efficient to handle this case separately with an on-body
     * flag rathert than use the exclusion set within the AABB structure */
  }

  int id1 = BODY (one->body)->id,
      id2 = BODY (two->body)->id,
      no1 = one->sgp - BODY (one->body)->sgp,
      no2 = two->sgp - BODY (two->body)->sgp;

  OPR pair = {MIN (id1, id2), MAX (id1, id2), MIN (no1, no2), MAX (no1, no2)};

  /* an excluded body pair ? */
  if (SET_Contains (aux->nobody, &pair, (SET_Compare) bodcmp)) return NULL;

  /* an excluded object pair ? */
  if (SET_Contains (aux->nogobj, &pair, (SET_Compare) gobjcmp)) return NULL;

  /* test topological adjacency */
  switch (GOBJ_Pair_Code (one, two))
  {
    case AABB_ELEMENT_ELEMENT:
    {
      if (one->body == two->body) /* assuming only one mesh per body is possible,
				     the two elements are within the same mesh;
				     exclude topologically adjacent elements */
	if (ELEMENT_Adjacent (one->sgp->gobj, two->sgp->gobj)) return NULL;
    }
    break;
    case AABB_CONVEX_CONVEX:
    {
      if (one->body == two->body) /* exclude topologically adjacent convices */
	if (CONVEX_Adjacent (one->sgp->gobj, two->sgp->gobj)) return NULL;
    }
    break;
    case AABB_SPHERE_SPHERE:
    {
      if (one->body == two->body) /* exclude topologically adjacent spheres */
	if (SPHERE_Adjacent (one->sgp->gobj, two->sgp->gobj)) return NULL;
    }
  }
 
  /* report overlap creation and get the user pointer */
  user = aux->create (aux->data, one, two);

  if (user) /* if the user pointer is NULL keep re-detecting the overlap */
  {
    MAP_Insert (aux->mapmem, &one->adj, two, user, NULL); /* map overlap symmetrically */
    MAP_Insert (aux->mapmem, &two->adj, one, user, NULL);
    /* one might think that only one adjacency (say box with bigger pointer
     * stores box with smaller pointer) could be used. nevertheless, for box
     * removal one needs to release all overlaps, and then the complete adjacency
     * informaltion is necessary; thus the symmetrical mapping */
  }

  return user;
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

#if MPI
/* number of boxes */
static int box_count (AABB *aabb, int *ierr)
{
  *ierr = ZOLTAN_OK;

  return aabb->nin + aabb->nlst;
}

/* list of box identifiers */
static void box_list (AABB *aabb, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr)
{
  BOX *box, **aux;
  BODY *bod;
  int i;

  free (aabb->aux);
  /* realloc auxiliary table */
  ERRMEM (aabb->aux = malloc ((aabb->nin + aabb->nlst) * sizeof (BOX*)));

  /* gather inserted boxes */
  for (aux = aabb->aux, i = 0, box = aabb->in; box; aux ++, i ++, box = box->next)
  {
    bod = box->body;

    global_ids [i * num_gid_entries] = bod->id;
    global_ids [i * num_gid_entries + 1] = box->sgp - bod->sgp; /* local SGP index */
    local_ids [i * num_lid_entries] = i;
    *aux = box;
  }

  /* gather current boxes */
  for (box = aabb->lst; box; aux ++, i ++, box = box->next)
  {
    bod = box->body;

    global_ids [i * num_gid_entries] = bod->id;
    global_ids [i * num_gid_entries + 1] = box->sgp - bod->sgp;
    local_ids [i * num_lid_entries] = i;
    *aux = box;
  }

  *ierr = ZOLTAN_OK;
}

/* number of spatial dimensions */
static int dimensions (AABB *aabb, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}

/* list of body extent low points */
static void lopoints (AABB *aabb, int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int num_dim, double *geom_vec, int *ierr)
{
  BOX **aux, *box;
  int i;

  aux = aabb->aux;

  for (i = 0; i < num_obj; i ++, geom_vec += num_dim)
  {
    box = aux [local_ids [i * num_lid_entries]];
    COPY (box->extents, geom_vec);
  }

  *ierr = ZOLTAN_OK;
}

/* pack migration data */
static void box_pack (BOX *box, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, box->bid);
  pack_int (isize, i, ints, box->kind & ~GOBJ_NEW);
  if (box->body) pack_int (isize, i, ints, box->sgp - BODY(box->body)->sgp);
  else pack_int (isize, i, ints, (int) box->sgp); /* detached box */
  pack_doubles (dsize, d, doubles, box->extents, 6);
}

/* unpack migration data */
static void* box_unpack (AABB *aabb, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  DOM *dom = aabb->dom;
  int bid, kind, isgp;
  double extents [6];
  BODY *bod;
  SGP *sgp;
  BOX *box;

  bid = unpack_int (ipos, i, ints);
  kind = unpack_int (ipos, i, ints);
  isgp = unpack_int (ipos, i, ints);
  unpack_doubles (dpos, d, doubles, extents, 6);

  bod = MAP_Find (dom->children, (void*) bid, NULL);
  if (!bod) ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) bid, NULL), "Invalid body id");

  sgp = &bod->sgp [isgp];

  if (sgp->box == NULL) /* do not insert boxes that are already attached */
  {
    switch (kind)
    {
    case GOBJ_ELEMENT: sgp->box = AABB_Insert (aabb, bod, kind, sgp, sgp->shp->data, (BOX_Extents_Update)ELEMENT_Extents); break;
    case GOBJ_CONVEX:  sgp->box = AABB_Insert (aabb, bod, kind, sgp, sgp->shp->data, (BOX_Extents_Update)CONVEX_Extents); break;
    case GOBJ_SPHERE:  sgp->box = AABB_Insert (aabb, bod, kind, sgp, sgp->shp->data, (BOX_Extents_Update)SPHERE_Extents); break;
    }

    box = sgp->box;

    ASSERT_DEBUG (box, "Invalid box kind");

    COPY6 (extents, box->extents);
  }

  return NULL;
}

/* return next pointer and realloc send memory if needed */
inline static COMOBJ* sendnext (int nsend, int *size, COMOBJ **send)
{
  if (nsend >= *size)
  {
    (*size) *= 2;
    ERRMEM (*send = realloc (*send, sizeof (COMOBJ [*size])));
  }

  return &(*send)[nsend];
}

/* synchronise exclusion set */
static void sync (SET **eset, MEM *oprmem, MEM *setmem, int nobody)
{
  int i, s, n, ncpu, step, *sizes, *disp, *set, *all, *ptr, *end;
  SET_Compare cmp;
  SET *item;

  if (nobody) 
  {
    step = 2;
    cmp = (SET_Compare)bodcmp;
  }
  else
  {
    step = 4;
    cmp = (SET_Compare)gobjcmp;
  }

  MPI_Comm_size (MPI_COMM_WORLD, &ncpu);

  s = SET_Size (*eset);
  ERRMEM (sizes = malloc (sizeof (int [ncpu])));
  ERRMEM (disp = malloc (sizeof (int [ncpu])));
  ERRMEM (set = malloc (sizeof (int [step * s])));

  /* sizes of all sets */
  MPI_Allgather (&s, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);

  for (disp [0] = 0, i = 0; i < ncpu; i ++)
    if (i < ncpu - 1) disp [i+1] = disp [i] + sizes [i];

  n = disp [ncpu-1] + sizes [ncpu-1]; /* total size */

  ERRMEM (all = malloc (sizeof (int [step * n])));

  /* copy the local set into a buffer */
  for (item = SET_First (*eset), ptr = set; item; item = SET_Next (item), ptr += step)
  {
    OPR *opr = item->data;

    ptr [0] = opr->bod1;
    ptr [1] = opr->bod2;

    if (step == 4)
    {
      ptr [2] = opr->sgp1;
      ptr [3] = opr->sgp2;
    }
  }

  /* gather all local sets into the common 'all' space */
  MPI_Allgatherv (set, step * s, MPI_INT, all, sizes, disp, MPI_INT, MPI_COMM_WORLD);

  /* complete the local set */
  for (ptr = all, end = all + step * n; ptr < end; ptr += step)
  {
    OPR pair = {ptr [0], ptr [1], 0, 0};

    if (step == 4)
    {
      pair.sgp1 = ptr [2];
      pair.sgp2 = ptr [3];
    }

    if (!SET_Contains (*eset, &pair, cmp))  /* if not found */
    {
      OPR *opr;
      
      ERRMEM (opr = MEM_Alloc (oprmem));
      *opr = pair;
      SET_Insert (setmem, eset, opr, cmp); /* insert new item */
    }
  }

  free (all);
  free (set);
  free (disp);
  free (sizes);
}

/* synchronise set */
static void syncset (SET **theset, MEM *setmem)
{
  int i, s, n, ncpu, *sizes, *disp, *set, *all, *ptr, *end;
  SET *item;

  MPI_Comm_size (MPI_COMM_WORLD, &ncpu);

  s = SET_Size (*theset);
  ERRMEM (sizes = malloc (sizeof (int [ncpu])));
  ERRMEM (disp = malloc (sizeof (int [ncpu])));
  ERRMEM (set = malloc (sizeof (int [s])));

  /* sizes of all sets */
  MPI_Allgather (&s, 1, MPI_INT, sizes, 1, MPI_INT, MPI_COMM_WORLD);

  for (disp [0] = 0, i = 0; i < ncpu; i ++)
    if (i < ncpu - 1) disp [i+1] = disp [i] + sizes [i];

  n = disp [ncpu-1] + sizes [ncpu-1]; /* total size */

  ERRMEM (all = malloc (sizeof (int [n])));

  /* copy the local set into a buffer */
  for (item = SET_First (*theset), ptr = set; item; item = SET_Next (item), ptr ++)
  {
    ptr [0] = (int) item->data;
  }

  /* gather all local sets into the common 'all' space */
  MPI_Allgatherv (set, s, MPI_INT, all, sizes, disp, MPI_INT, MPI_COMM_WORLD);

  /* complete the local set */
  for (ptr = all, end = all + n; ptr < end; ptr ++)
  {
    if (!SET_Contains (*theset, (void*) ptr [0], NULL))  /* if not found */
      SET_Insert (setmem, theset, (void*) ptr [0], NULL); /* insert new item */
  }

  free (all);
  free (set);
  free (disp);
  free (sizes);
}

/* attach deteached boxes */
static void attach_detached (AABB *aabb)
{
  BOX *box, *next;
  BODY *bod;
  DOM *dom;

  dom = aabb->dom;

  for (box = aabb->in; box; box = next)
  {
    next = box->next;

    if (box->body == NULL)
    {
      bod = MAP_Find (dom->children, (void*) box->bid, NULL);
      if (!bod) bod = MAP_Find (dom->idb, (void*) box->bid, NULL);

      if (bod)
      {
	box->body = bod;
	box->sgp = &bod->sgp [(int) box->sgp];
      }
      else AABB_Delete (aabb, box);
    }
  }

  for (box = aabb->lst; box; box = next)
  {
    next = box->next;

    if (box->body == NULL)
    {
      bod = MAP_Find (dom->children, (void*) box->bid, NULL);
      if (!bod) bod = MAP_Find (dom->idb, (void*) box->bid, NULL);

      if (bod)
      {
	box->body = bod;
	box->sgp = &bod->sgp [(int) box->sgp];
      }
      else AABB_Delete (aabb, box);
    }
  }
}

/* balance boxes */
static void aabb_balance (AABB *aabb)
{
  int flag, size;
  BOX *box;

  /* synchronise body pair exclusion set */
  MPI_Allreduce (&aabb->nobody_modified, &flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (flag) sync (&aabb->nobody, &aabb->setmem, &aabb->oprmem, 1), aabb->nobody_modified = 0;

  /* synchronise geometric object pair exclusion set */
  MPI_Allreduce (&aabb->nogobj_modified, &flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (flag) sync (&aabb->nogobj, &aabb->setmem, &aabb->oprmem, 0), aabb->nogobj_modified = 0;

  size = SET_Size (aabb->delbod);
  MPI_Allreduce (&size, &flag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  if (flag)
  {
    SET *item, *jtem;
    BOX *next;

    /* synchronise body deletion set */
    syncset (&aabb->delbod, &aabb->setmem);

    /* remove related items from exclusion sets */
    for (item = SET_First (aabb->delbod); item; item = SET_Next (item))
    {
      OPR *op;

      for (jtem = SET_First (aabb->nobody); jtem; )
      {
	op = jtem->data;

	if (op->bod1 == (unsigned int) item->data ||
	    op->bod1 == (unsigned int) item->data)
	  jtem = SET_Delete_Node (&aabb->setmem, &aabb->nobody, jtem);
	else jtem = SET_Next (jtem);
      }

      for (jtem = SET_First (aabb->nogobj); jtem; )
      {
	op = jtem->data;

	if (op->bod1 == (unsigned int) item->data ||
	    op->bod1 == (unsigned int) item->data)
	  jtem = SET_Delete_Node (&aabb->setmem, &aabb->nobody, jtem);
	else jtem = SET_Next (jtem);
      }
    }

    /* delete boxes attached to deleted bodies */
    for (box = aabb->in; box; box = next)
    {
      next = box->next;
      if (SET_Contains (aabb->delbod, (void*) box->bid, NULL)) AABB_Delete (aabb, box);
    }
    for (box = aabb->lst; box; box = next)
    {
      next = box->next;
      if (SET_Contains (aabb->delbod, (void*) box->bid, NULL)) AABB_Delete (aabb, box);
    }

    /* empty deletion set */
    SET_Free (&aabb->setmem, &aabb->delbod);
  }

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

  /* update box partitioning using low points of their extents */
  ASSERT (Zoltan_LB_Balance (aabb->zol, &changes, &num_gid_entries, &num_lid_entries,
	  &num_import, &import_global_ids, &import_local_ids, &import_procs,
	  &num_export, &export_global_ids, &export_local_ids, &export_procs) == ZOLTAN_OK, ERR_ZOLTAN);

  Zoltan_LB_Free_Data (&import_global_ids, &import_local_ids, &import_procs,
                       &export_global_ids, &export_local_ids, &export_procs);

  /* balance children bodies in the domin */
  DOM_Balance_Children (aabb->dom, aabb->zol);

  int nsend, nrecv, i, j, naux;
  int *procs, numprocs, rank;
  COMOBJ *send, *recv, *ptr;
  SET *del, *item;
  BOX **aux;
  DOM *dom;

  nsend = 0;
  del = NULL;
  dom = aabb->dom;
  rank = dom->rank;
  naux = box_count (aabb, &i);
  size = MIN (naux, 128);
  ERRMEM (send = malloc (sizeof (COMOBJ [size])));
  ERRMEM (procs = malloc (sizeof (int [DOM(aabb->dom)->size])));

  for (aux = aabb->aux, ptr = send, i = 0; i < naux; i ++)
  {
    double *e = aux [i]->extents;

    ASSERT (Zoltan_LB_Box_Assign (aabb->zol, e[0], e[1], e[2], e[3], e[4], e[5], procs, &numprocs) == ZOLTAN_OK, ERR_ZOLTAN);

    for (flag = 1, j = 0; j < numprocs; j ++)
    {
      if (procs [j] != rank) /* exported */
      {
	ptr->rank = procs [j];
	ptr->o = aux [i];
	ptr = sendnext (++ nsend, &size, &send);
      }
      else flag = 0; /* should stay here */
    }

    if (flag) SET_Insert (&aabb->setmem, &del, aux [i], NULL);
  }

  /* communicate migration data (unpacking inserts the imported boxes) */
  COMOBJS (MPI_COMM_WORLD, TAG_AABB_BALANCE, (OBJ_Pack)box_pack, aabb, (OBJ_Unpack)box_unpack, send, nsend, &recv, &nrecv);

  /* delete exported boxes */
  for (item = SET_First (del); item; item = SET_Next (item)) AABB_Delete (aabb, item->data);

  /* there should be no detached bodies left */
  for (box = aabb->in; box; box = box->next) ASSERT_DEBUG (box->body, "A deteached box found");
  for (box = aabb->lst; box; box = box->next) ASSERT_DEBUG (box->body, "A deteached box found");

  SET_Free (&aabb->setmem, &del);
  free (send);
  free (recv);
}

/* create MPI related data */
static void create_mpi (AABB *aabb)
{
  aabb->delbod = NULL;
  aabb->nobody_modified = 0;
  aabb->nogobj_modified = 0;
  aabb->aux = NULL;

  /* zoltan context for body partitioning */
  ASSERT (aabb->zol = Zoltan_Create (MPI_COMM_WORLD), ERR_ZOLTAN);

  /* general parameters */
  Zoltan_Set_Param (aabb->zol, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param (aabb->zol, "DEBUG_MEMORY", "0");
  Zoltan_Set_Param (aabb->zol, "NUM_GID_ENTRIES", "2"); /* body id, sgp index */
  Zoltan_Set_Param (aabb->zol, "NUM_LID_ENTRIES", "1"); /* indices in aux table */
 
  /* load balancing parameters */
  Zoltan_Set_Param (aabb->zol, "LB_METHOD", "RCB");
  Zoltan_Set_Param (aabb->zol, "IMBALANCE_TOL", "1.2");
  Zoltan_Set_Param (aabb->zol, "AUTO_MIGRATE", "FALSE"); /* we shall use COMOBJS */
  Zoltan_Set_Param (aabb->zol, "RETURN_LISTS", "NONE");

  /* RCB parameters */
  Zoltan_Set_Param (aabb->zol, "RCB_OVERALLOC", "1.3");
  Zoltan_Set_Param (aabb->zol, "RCB_REUSE", "1");
  Zoltan_Set_Param (aabb->zol, "RCB_OUTPUT_LEVEL", "0");
  Zoltan_Set_Param (aabb->zol, "CHECK_GEOM", "1");
  Zoltan_Set_Param (aabb->zol, "KEEP_CUTS", "1");

  /* callbacks */
  Zoltan_Set_Fn (aabb->zol, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) box_count, aabb);
  Zoltan_Set_Fn (aabb->zol, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)()) box_list, aabb);
  Zoltan_Set_Fn (aabb->zol, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) dimensions, aabb);
  Zoltan_Set_Fn (aabb->zol, ZOLTAN_GEOM_MULTI_FN_TYPE, (void (*)()) lopoints, aabb);
}

/* destroy MPI related data */
static void destroy_mpi (AABB *aabb)
{
  free (aabb->aux);

  Zoltan_Destroy (&aabb->zol);
}
#endif

/* create box overlap driver data */
AABB* AABB_Create (int size)
{
  AABB *aabb;

  ERRMEM (aabb = malloc (sizeof (AABB)));
  aabb->lst = aabb->in = aabb->out= NULL;
  aabb->tab = NULL;

  MEM_Init (&aabb->boxmem, sizeof (BOX), MIN (size, SIZE));
  MEM_Init (&aabb->mapmem, sizeof (MAP), MIN (size, SIZE));
  MEM_Init (&aabb->setmem, sizeof (SET), MIN (size, SIZE));
  MEM_Init (&aabb->oprmem, sizeof (OPR), MIN (size, SIZE));

  aabb->nobody = aabb->nogobj= NULL;

  aabb->nin = 0;
  aabb->nlst = 0;
  aabb->ntab = 0;
  aabb->modified = 0;
  aabb->swp = aabb->hsh = NULL;

#if MPI
  create_mpi (aabb);
#endif

  return aabb;
}

/* insert geometrical object */
BOX* AABB_Insert (AABB *aabb, void *body, GOBJ kind, SGP *sgp, void *data, BOX_Extents_Update update)
{
  BOX *box;

  ERRMEM (box = MEM_Alloc (&aabb->boxmem));
  box->data = data;
  box->update = update;
  box->bid = BODY(body)->id;
  box->kind = kind | GOBJ_NEW; /* mark as new */
  box->body = body;
  box->sgp = sgp;

  /* include into the
   * insertion list */
  box->next = aabb->in;
  if (aabb->in) aabb->in->prev= box;
  aabb->in = box;
  aabb->nin ++;

  /* set modified */
  aabb->modified = 1;

  return box;
}

/* delete an object */
void AABB_Delete (AABB *aabb, BOX *box)
{
#if MPI
  if (!box) return; /* possible for a child body or a parent body whose box has migrated away */
  else if (box->body) box->sgp->box = NULL; /* invalidate box pointer */
#endif

  if (box->kind & GOBJ_NEW) /* still in the insertion list */
  {
    /* remove box from the insertion list */
    if (box->prev) box->prev->next = box->next;
    else aabb->in = box->next;
    if (box->next) box->next->prev = box->prev;
    aabb->nin --;

    /* free it here */
    MEM_Free (&aabb->boxmem, box); 
  }
  else
  {
    /* remove box from the current list */
    if (box->prev) box->prev->next = box->next;
    else aabb->lst = box->next;
    if (box->next) box->next->prev = box->prev;
    aabb->nlst --; /* decrease count */

    /* insert it into
     * the deletion list */
    box->next = aabb->out;
    aabb->out = box;
  }

  /* set modified */
  aabb->modified = 1;
}

/* insert a body */
void AABB_Insert_Body (AABB *aabb, void *body)
{
  SGP *sgp, *sgpe;
  BODY *bod;

  for (bod = body, sgp = bod->sgp, sgpe = sgp + bod->nsgp; sgp < sgpe; sgp ++)
  {
    sgp->box = AABB_Insert (aabb, bod, gobj_kind (sgp), sgp,
                            sgp->shp->data, extents_update (sgp));
  }
}

/* delete a body */
void AABB_Delete_Body (AABB *aabb, void *body)
{
  SGP *sgp, *sgpe;
  BODY *bod;

  for (bod = body, sgp = bod->sgp, sgpe = sgp + bod->nsgp; sgp < sgpe; sgp ++)
  {
    AABB_Delete (aabb, sgp->box);
  }

#if MPI
  SET_Insert (&aabb->setmem, &aabb->delbod, (void*) bod->id, NULL);
#endif
}

#if MPI
/* detach a body */
void AABB_Detach_Body (AABB *aabb, void *body)
{
  SGP *sgp, *sgpe;
  BODY *bod;
  BOX *box;

  for (bod = body, sgp = bod->sgp, sgpe = sgp + bod->nsgp; sgp < sgpe; sgp ++)
  {
    box = sgp->box;

    if (box)
    {
      box->sgp = (SGP*) (box->sgp - bod->sgp);
      box->body = NULL;
    }
  }
}
#endif

/* update state => detect created and released overlaps */
void AABB_Update (AABB *aabb, BOXALG alg, void *data, BOX_Overlap_Create create, BOX_Overlap_Release release)
{
  struct auxdata aux = {aabb->nobody, aabb->nogobj, &aabb->mapmem, data, create};
  BOX *box, *next, *adj;
  MAP *item, *jtem;

#if MPI
  attach_detached (aabb);
#endif

  /* update extents */
  for (box = aabb->lst; box; box = box->next) /* for each current box */
    box->update (box->data, box->sgp->gobj, box->extents);

  for (box = aabb->in; box; box = box->next) /* for each inserted box */
    box->update (box->data, box->sgp->gobj, box->extents);

#if MPI
  aabb_balance (aabb);
#endif

  for (box = aabb->in; box; box = box->next) /* for each inserted box */
    box->kind &= ~GOBJ_NEW; /* not new any more */

  /* check for released overlaps */
  for (box = aabb->out; box; box = next) /* for each deleted box */
  {
    next = box->next;

    for (item = MAP_First (box->adj); item; item = MAP_Next (item)) /* for each adjacent box */
    {
      adj = item->key;
#if MPI
      if (box->body) /* not detached */
#endif
      release (data, box, adj, item->data); /* release overlap */
      MAP_Delete (&aabb->mapmem, &adj->adj, box, NULL); /* remove 'box' from 'adj's adjacency */
    }

    MAP_Free (&aabb->mapmem, &box->adj); /* free adjacency */
    MEM_Free (&aabb->boxmem, box); /* free box */
  }
  aabb->out = NULL; /* list emptied */

  for (box = aabb->lst; box; box = box->next) /* for each current box */
  {
    for (item = MAP_First (box->adj); item; item = jtem) /* for each adjacent box */
    {
      adj = item->key;
      if (no_overlap (box->extents, adj->extents))
      {
	release (data, box, adj, item->data); /* release overlap */
	MAP_Delete (&aabb->mapmem, &adj->adj, box, NULL); /* remove 'box' from 'adj's adjacency */
	jtem = MAP_Delete_Node (&aabb->mapmem, &box->adj, item); /* remove 'adj' from 'box's adjacency */
      }
      else jtem = MAP_Next (item);
    }
  }

  if (aabb->modified) /* merge insertion and curent lists, update pointer table */
  {
    BOX **b;

    if (aabb->in)
    {
      for (box = aabb->in; box->next; box = box->next); /* rewind till the last inserted one */
      box->next = aabb->lst; /* link it with the current list */
      if (aabb->lst) aabb->lst->prev = box; /* set previous link */
      aabb->lst = aabb->in; /* list appended */
      aabb->nlst += aabb->nin; /* increase number of current boxes */
      aabb->in = NULL; /* list emptied */
      aabb->nin = 0;
      aabb->ntab = aabb->nlst; /* resize the pointer table */
      ERRMEM (aabb->tab = realloc (aabb->tab, sizeof (BOX*) * aabb->ntab));
    }

    for (box = aabb->lst, b = aabb->tab; box; box = box->next, b ++) *b = box; /* overwrite box pointers */
  }

  /* the algorithm
   * specific part */
  switch (alg)
  {
    case HYBRID:
    {
      hybrid (aabb->tab, aabb->ntab, &aux, (BOX_Overlap_Create)local_create);
    }
    break;
    case HASH3D:
    {
      if (!aabb->hsh) aabb->hsh = HASH_Create (aabb->ntab);

      HASH_Do (aabb->hsh, aabb->ntab, aabb->tab,
	       &aux, (BOX_Overlap_Create)local_create); 
    }
    break;
    case SWEEP_HASH2D_LIST:
    case SWEEP_HASH1D_XYTREE:
    case SWEEP_HASH2D_XYTREE:
    case SWEEP_XYTREE:
    {
      if (!aabb->swp) aabb->swp = SWEEP_Create (aabb->ntab, (DRALG)alg);

      SWEEP_Do (aabb->swp, (DRALG)alg, aabb->modified, aabb->ntab,
	        aabb->tab, &aux, (BOX_Overlap_Create)local_create); 
    }
    break;
  }

  /* set unmodified */
  aabb->modified = 0;
}

/* never report overlaps betweem this pair of bodies (given by identifiers) */
void AABB_Exclude_Body_Pair (AABB *aabb, unsigned int id1, unsigned int id2)
{
  OPR *opr;
  
  ERRMEM (opr = MEM_Alloc (&aabb->oprmem));
  opr->bod1 = MIN (id1, id2);
  opr->bod2 = MAX (id1, id2);
  SET_Insert (&aabb->setmem, &aabb->nobody, opr, (SET_Compare) bodcmp);

#if MPI
  aabb->nobody_modified = 1;
#endif
}

/* never report overlaps betweem this pair of objects (bod1, sgp1), (bod1, sgp2) */
void AABB_Exclude_Gobj_Pair (AABB *aabb, unsigned int bod1, int sgp1, unsigned int bod2, int sgp2)
{
  OPR *opr;
  
  ERRMEM (opr = MEM_Alloc (&aabb->oprmem));
  opr->bod1 = MIN (bod1, bod2);
  opr->bod2 = MAX (bod1, bod2);
  opr->sgp1 = MIN (sgp1, sgp2);
  opr->sgp2 = MAX (sgp1, sgp2);
  SET_Insert (&aabb->setmem, &aabb->nogobj, opr, (SET_Compare) gobjcmp);

#if MPI
  aabb->nogobj_modified = 1;
#endif
}

/* break box adjacency if boxes associated with the objects are adjacent */
void AABB_Break_Adjacency (AABB *aabb, BOX *one, BOX *two)
{
  /* mutually delete from adjacency lists */
  MAP_Delete (&aabb->mapmem, &one->adj, two, NULL);
  MAP_Delete (&aabb->mapmem, &two->adj, one, NULL);
}

/* release memory */
void AABB_Destroy (AABB *aabb)
{
  free (aabb->tab);

  MEM_Release (&aabb->boxmem);
  MEM_Release (&aabb->mapmem);
  MEM_Release (&aabb->setmem);
  MEM_Release (&aabb->oprmem);

  if (aabb->swp) SWEEP_Destroy (aabb->swp);
  if (aabb->hsh) HASH_Destroy (aabb->hsh);

#if MPI
  destroy_mpi (aabb);
#endif

  free (aabb);
}
