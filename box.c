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

#include "sol.h"
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
#include "put.h"
#include "pck.h"
#include "com.h"
#endif

#define SIZE 128 /* mempool size */

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
  BODY *onebod = one->body, *twobod = two->body;

#if MPI
  if ((onebod->flags & BODY_CHILD) && (twobod->flags & BODY_CHILD)) return NULL; /* only parent-parent and parent-child contacts are valid */
  else if (onebod->flags & BODY_CHILD)
  {
    if (SET_Contains (twobod->children, (void*) (long) onebod->rank, NULL)) /* check whether the parent's body child and child's parent are on the same processor */
    {
      if (onebod->dom->rank < onebod->rank) return NULL; /* if so, use processor ordering in order to avoid duplicated contact detection */
    }
  }
  else if (twobod->flags & BODY_CHILD)
  {
    if (SET_Contains (onebod->children, (void*) (long) twobod->rank, NULL))
    {
      if (twobod->dom->rank < twobod->rank) return NULL;
    }
  }
#endif

  /* check if these are two obstacles => no need to report overlap */
  if (onebod->kind == OBS && twobod->kind == OBS) return NULL;

  /* boxes already adjacent ? */
  if (MAP_Find (one->adj, two, NULL)) return NULL;

  /* self-contact ? <=> one->body == two->body */
  if (onebod == twobod && !(onebod->flags & BODY_DETECT_SELF_CONTACT)) return NULL;

  int id1 = onebod->id,
      id2 = twobod->id,
      no1 = one->sgp - onebod->sgp,
      no2 = two->sgp - twobod->sgp;

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
      if (onebod == twobod) /* assuming only one mesh per body is possible,
			       the two elements are within the same mesh;
			       exclude topologically adjacent elements */
	if (ELEMENT_Adjacent (one->sgp->gobj, two->sgp->gobj)) return NULL;
    }
    break;
    case AABB_CONVEX_CONVEX:
    {
      if (onebod == twobod) /* exclude topologically adjacent convices */
	if (CONVEX_Adjacent (one->sgp->gobj, two->sgp->gobj)) return NULL;
    }
    break;
    case AABB_SPHERE_SPHERE:
    {
      if (onebod == twobod) /* exclude topologically adjacent spheres */
	if (SPHERE_Adjacent (one->sgp->gobj, two->sgp->gobj)) return NULL;
    }
  }
 
  /* report overlap creation and get the user pointer */
  void *user = aux->create (aux->data, one, two);

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

#if MPI
/* weight of a box */
static int box_weight (BOX *box)
{
  return 1 + MAP_Size (box->adj); /* using adjacency size attempts to balance contacts together with boxes */
}

/* number of boxes */
static int box_count (AABB *aabb, int *ierr)
{
  *ierr = ZOLTAN_OK;

  return aabb->boxnum;
}

/* list of box identifiers */
static void box_list (AABB *aabb, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr)
{
  BOX *box, **aux;
  BODY *bod;
  int i;

  /* realloc auxiliary table */
  free (aabb->aux);
  ERRMEM (aabb->aux = malloc (aabb->boxnum * sizeof (BOX*)));

  /* gather current boxes */
  for (aux = aabb->aux, box = aabb->lst, i = 0; box; aux ++, i ++, box = box->next)
  {
    bod = box->body;

    global_ids [i * num_gid_entries] = bod->id;
    global_ids [i * num_gid_entries + 1] = box->sgp - bod->sgp;
    local_ids [i * num_lid_entries] = i;
    obj_wgts [i * wgt_dim] = box_weight (box);
    box->update (box->data, box->sgp->gobj, box->extents);
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
static void boxpoints (AABB *aabb, int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int num_dim, double *geom_vec, int *ierr)
{
  BOX **aux, *box;
  int i;

  aux = aabb->aux;

  for (i = 0; i < num_obj; i ++, geom_vec += num_dim)
  {
    box = aux [local_ids [i * num_lid_entries]];
    MID (box->extents, box->extents + 3, geom_vec);
  }

  *ierr = ZOLTAN_OK;
}

/* pack migration data */
static void box_pack (BOX *box, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, box->body->id);
  pack_int (isize, i, ints, box->sgp - box->body->sgp);
}

/* unpack migration data */
static void* box_unpack (AABB *aabb, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  DOM *dom = aabb->dom;
  int id, kind, nsgp;
  BOX *box;
  BODY *bod;
  SGP *sgp;

  id = unpack_int (ipos, i, ints);
  nsgp = unpack_int (ipos, i, ints);

  ASSERT_DEBUG_EXT (bod = MAP_Find (dom->allbodies, (void*) (long) id, NULL), "Invalid body id");
  ASSERT_DEBUG (bod->flags & (BODY_PARENT|BODY_CHILD), "Neither child nor parent");
  sgp = &bod->sgp [nsgp];
  kind = gobj_kind  (sgp);

  if (sgp->box == NULL) /* do not insert a boxe that already exists */
  {
    switch (kind)
    {
    case GOBJ_ELEMENT: box = AABB_Insert (aabb, bod, kind, sgp, sgp->shp->data, (BOX_Extents_Update)ELEMENT_Extents); break;
    case GOBJ_CONVEX:  box = AABB_Insert (aabb, bod, kind, sgp, sgp->shp->data, (BOX_Extents_Update)CONVEX_Extents); break;
    case GOBJ_SPHERE:  box = AABB_Insert (aabb, bod, kind, sgp, sgp->shp->data, (BOX_Extents_Update)SPHERE_Extents); break;
    }

    box->update (box->data, box->sgp->gobj, box->extents);

    return box;
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

/* balance boxes */
static void aabb_balance (AABB *aabb, void *data, BOX_Overlap_Release release)
{
  int *procs, numprocs, rank, size;
  int nsend, nrecv, i, j, flag;
  SET *delcon, *delbox, *item;
  COMOBJ *send, *recv, *ptr;
  BOX **aux, *box;
  MAP *jtem;
  double *e;
  DOM *dom;

  nsend = 0;
  size = 256;
  delbox = NULL;
  delcon = NULL;
  dom = aabb->dom;
  rank = dom->rank;
  ERRMEM (send = malloc (sizeof (COMOBJ [size])));
  ERRMEM (procs = malloc (sizeof (int [aabb->dom->ncpu])));

  for (aux = aabb->aux, ptr = send, i = 0; i < aabb->boxnum; i ++)
  {
    box = aux [i];
    e = box->extents;

    ASSERT (Zoltan_LB_Box_Assign (aabb->zol, e[0], e[1], e[2], e[3], e[4], e[5], procs, &numprocs) == ZOLTAN_OK, ERR_ZOLTAN);

    for (flag = 1, j = 0; j < numprocs; j ++)
    {
      if (procs [j] != rank) /* exported */
      {
        if (!SET_Contains (box->ranks, (void*) (long) procs [j], NULL)) /* wasn't yet sent there */
	{
	  ptr->rank = procs [j];
	  ptr->o = box;
	  ptr = sendnext (++ nsend, &size, &send);
	  SET_Insert (&aabb->setmem, &box->ranks, (void*) (long) procs [j], NULL);
	}
      }
      else flag = 0; /* should stay here */
    }

    if (flag) SET_Insert (&aabb->setmem, &delbox, box, NULL); /* schedule for deletion */
  }

  /* communicate boxes */
  dom->bytes += COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack)box_pack, aabb, (OBJ_Unpack)box_unpack, send, nsend, &recv, &nrecv);

  /* delete exported boxes */
  for (item = SET_First (delbox); item; item = SET_Next (item))
  {
    box = item->data;

    for (jtem = MAP_First (box->adj); jtem; jtem = MAP_Next (jtem))
      SET_Insert (&aabb->setmem, &delcon, jtem->data, NULL); /* gather contacts to be released */

    AABB_Delete (aabb, box); /* delete box */
  }

  /* release abandoned contacts */
  for (item = SET_First (delcon); item; item = SET_Next (item))
  {
    release (data, item->data);
  }

  SET_Free (&aabb->setmem, &delbox);
  SET_Free (&aabb->setmem, &delcon);
  free (procs);
  free (send);
  free (recv);
}

/* create MPI related data */
static void create_mpi (AABB *aabb)
{
  aabb->aux = NULL;

  /* zoltan context for body partitioning */
  ASSERT (aabb->zol = Zoltan_Create (MPI_COMM_WORLD), ERR_ZOLTAN);

  /* general parameters */
  Zoltan_Set_Param (aabb->zol, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param (aabb->zol, "DEBUG_MEMORY", "0");
  Zoltan_Set_Param (aabb->zol, "NUM_GID_ENTRIES", "2"); /* body id, sgp index */
  Zoltan_Set_Param (aabb->zol, "NUM_LID_ENTRIES", "1"); /* indices in aux table */
  Zoltan_Set_Param (aabb->zol, "OBJ_WEIGHT_DIM", "1");
 
  /* load balancing parameters */
  Zoltan_Set_Param (aabb->zol, "LB_METHOD", "RCB");
  Zoltan_Set_Param (aabb->zol, "IMBALANCE_TOL", "1.3");
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
  Zoltan_Set_Fn (aabb->zol, ZOLTAN_GEOM_MULTI_FN_TYPE, (void (*)()) boxpoints, aabb);
}

/* destroy MPI related data */
static void destroy_mpi (AABB *aabb)
{
  free (aabb->aux);

  Zoltan_Destroy (&aabb->zol);
}
#endif

/* algorithm name */
char* AABB_Algorithm_Name (BOXALG alg)
{
  switch (alg)
  {
  case SWEEP_HASH2D_LIST: return "SWEEP_HASH2D_LIST";
  case SWEEP_HASH2D_XYTREE: return "SWEEP_HASH2D_XYTREE";
  case SWEEP_XYTREE: return "SWEEP_XYTREE";
  case SWEEP_HASH1D_XYTREE: return "SWEEP_HASH1D_XYTREE";
  case HYBRID: return "HYBRID";
  case HASH3D: return "HASH3D";
  }

  return NULL;
}

/* create box overlap driver data */
AABB* AABB_Create (int size)
{
  AABB *aabb;

  ERRMEM (aabb = malloc (sizeof (AABB)));
  MEM_Init (&aabb->boxmem, sizeof (BOX), MIN (size, SIZE));
  MEM_Init (&aabb->mapmem, sizeof (MAP), MIN (size, SIZE));
  MEM_Init (&aabb->setmem, sizeof (SET), MIN (size, SIZE));
  MEM_Init (&aabb->oprmem, sizeof (OPR), MIN (size, SIZE));
  aabb->lst = NULL;
  aabb->tab = NULL;
  aabb->nobody = NULL;
  aabb->nogobj= NULL;
  aabb->boxnum = 0;
  aabb->tabsize = 0;
  aabb->modified = 0;
  aabb->swp = NULL;
  aabb->hsh = NULL;

#if MPI
  create_mpi (aabb);
#endif

  return aabb;
}

/* insert geometrical object */
BOX* AABB_Insert (AABB *aabb, BODY *body, GOBJ kind, SGP *sgp, void *data, BOX_Extents_Update update)
{
  BOX *box;

  ERRMEM (box = MEM_Alloc (&aabb->boxmem));
  box->data = data;
  box->update = update;
  box->kind = kind;
  box->body = body;
  box->sgp = sgp;
  sgp->box = box;

  /* insert into list */
  box->next = aabb->lst;
  if (aabb->lst) aabb->lst->prev= box;
  aabb->lst = box;

  /* incrememt */
  aabb->boxnum ++;

  /* set modified */
  aabb->modified = 1;

  return box;
}

/* delete an object */
void AABB_Delete (AABB *aabb, BOX *box)
{
#if MPI
  if (box == NULL) return; /* possible for a body whose box has migrated away */

  box->sgp->box = NULL; /* invalidate pointer */

  SET_Free (&aabb->setmem, &box->ranks); /* free ranks set */
#endif

  MAP *item;
  BOX *adj;

  for (item = MAP_First (box->adj); item; item = MAP_Next (item)) /* for each adjacent box */
  {
    adj = item->key;
    MAP_Delete (&aabb->mapmem, &adj->adj, box, NULL); /* remove 'box' from 'adj's adjacency */
  }

  /* free adjacency */
  MAP_Free (&aabb->mapmem, &box->adj);

  /* remove from list */
  if (box->prev) box->prev->next = box->next;
  else aabb->lst = box->next;
  if (box->next) box->next->prev = box->prev;

  /* decrement */
  aabb->boxnum --;

  /* free box */
  MEM_Free (&aabb->boxmem, box);

  /* set modified */
  aabb->modified = 1;
}

/* insert a body */
void AABB_Insert_Body (AABB *aabb, BODY *body)
{
  SGP *sgp, *sgpe;

  for (sgp = body->sgp, sgpe = sgp + body->nsgp; sgp < sgpe; sgp ++)
  {
    AABB_Insert (aabb, body, gobj_kind (sgp), sgp, sgp->shp->data, SGP_Extents_Update (sgp));
  }
}

/* delete a body */
void AABB_Delete_Body (AABB *aabb, BODY *body)
{
  SGP *sgp, *sgpe;

  for (sgp = body->sgp, sgpe = sgp + body->nsgp; sgp < sgpe; sgp ++)
  {
    AABB_Delete (aabb, sgp->box);
  }
}

/* update state => detect created and released overlaps */
void AABB_Update (AABB *aabb, BOXALG alg, void *data, BOX_Overlap_Create create, BOX_Overlap_Release release)
{
  struct auxdata aux = {aabb->nobody, aabb->nogobj, &aabb->mapmem, data, create};
  SET *del, *jtem;
  BOX *box, *adj;
  MAP *item;

#if MPI
  if (aabb->dom->rank == 0)
#endif
  if (aabb->dom && aabb->dom->verbose) printf ("CONDET (%s) ... ", AABB_Algorithm_Name (alg)), fflush (stdout);
  
#if MPI
  if (aabb->dom) SOLFEC_Timer_Start (aabb->dom->solfec, "PARBAL");

  /* rebalance boxes */
  aabb_balance (aabb, data, release);

  if (aabb->dom) SOLFEC_Timer_End (aabb->dom->solfec, "PARBAL"), SOLFEC_Timer_Start (aabb->dom->solfec, "CONDET");
#else
  if (aabb->dom) SOLFEC_Timer_Start (aabb->dom->solfec, "CONDET");

  /* update extents */
  for (box = aabb->lst; box; box = box->next) /* for each current box */
    box->update (box->data, box->sgp->gobj, box->extents);
#endif

  for (box = aabb->lst; box; box = box->next) /* for each current box */
  {
    for (del = NULL, item = MAP_First (box->adj); item; item = MAP_Next (item)) /* for each adjacent box */
    {
      adj = item->key;
      if (no_overlap (box->extents, adj->extents))
      {
	release (data, item->data); /* release overlap */
	MAP_Delete (&aabb->mapmem, &adj->adj, box, NULL); /* remove 'box' from 'adj's adjacency */
	SET_Insert (&aabb->setmem, &del, adj, NULL); /* gather adjacent boxes to be released */
      }
    }

    for (jtem = SET_First (del); jtem; jtem = SET_Next (jtem))
      MAP_Delete (&aabb->mapmem, &box->adj, jtem->data, NULL); /* release symmetric overlaps */
  }

  if (aabb->modified) /* merge insertion and curent lists, update pointer table */
  {
    BOX **b;

    if (aabb->boxnum > aabb->tabsize)
    {
      free (aabb->tab);
      aabb->tabsize = 2 * aabb->boxnum;
      ERRMEM (aabb->tab = malloc (sizeof (BOX*) * aabb->tabsize));
    }

    for (box = aabb->lst, b = aabb->tab; box; box = box->next, b ++) *b = box; /* overwrite box pointers */
  }

#if DEBUG && MPI
  for (box = aabb->lst; box; box = box->next)
  {
    ASSERT_DEBUG (box->body->flags & (BODY_PARENT|BODY_CHILD), "BOX is attached to a DUMMY");
  }
#endif

  /* regardless of the current algorithm notify sweep-plane about the change */
  if (aabb->modified && aabb->swp) SWEEP_Changed (aabb->swp);

  /* the algorithm
   * specific part */
  switch (alg)
  {
    case HYBRID:
    {
      hybrid (aabb->tab, aabb->boxnum, &aux, (BOX_Overlap_Create)local_create);
    }
    break;
    case HASH3D:
    {
      if (!aabb->hsh) aabb->hsh = HASH_Create (aabb->boxnum);

      HASH_Do (aabb->hsh, aabb->boxnum, aabb->tab,
	       &aux, (BOX_Overlap_Create)local_create); 
    }
    break;
    case SWEEP_HASH2D_LIST:
    case SWEEP_HASH1D_XYTREE:
    case SWEEP_HASH2D_XYTREE:
    case SWEEP_XYTREE:
    {
      if (!aabb->swp) aabb->swp = SWEEP_Create (aabb->boxnum, (DRALG)alg);

      SWEEP_Do (aabb->swp, (DRALG)alg, aabb->boxnum, aabb->tab, &aux, (BOX_Overlap_Create)local_create); 
    }
    break;
  }

  /* set unmodified */
  aabb->modified = 0;

  if (aabb->dom)
  {
    DOM_Sparsify_Contacts (aabb->dom);

    SOLFEC_Timer_End (aabb->dom->solfec, "CONDET");
  }
}

/* never report overlaps betweem this pair of bodies (given by identifiers) */
void AABB_Exclude_Body_Pair (AABB *aabb, unsigned int id1, unsigned int id2)
{
  OPR *opr;
  
  ERRMEM (opr = MEM_Alloc (&aabb->oprmem));
  opr->bod1 = MIN (id1, id2);
  opr->bod2 = MAX (id1, id2);
  SET_Insert (&aabb->setmem, &aabb->nobody, opr, (SET_Compare) bodcmp);
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

#if MPI
/* geomtric partitioning */
void AABB_Partition (AABB *aabb)
{
  char tol [128];

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

  /* update imbalance tolerance */
  snprintf (tol, 128, "%g", aabb->dom->imbalance_tolerance);
  Zoltan_Set_Param (aabb->zol, "IMBALANCE_TOL", tol);

  /* update box partitioning */
  ASSERT (Zoltan_LB_Balance (aabb->zol, &changes, &num_gid_entries, &num_lid_entries,
	  &num_import, &import_global_ids, &import_local_ids, &import_procs,
	  &num_export, &export_global_ids, &export_local_ids, &export_procs) == ZOLTAN_OK, ERR_ZOLTAN);

  
  /* free auxiliary data *//* free auxiliary data */
  Zoltan_LB_Free_Data (&import_global_ids, &import_local_ids, &import_procs,
		       &export_global_ids, &export_local_ids, &export_procs);
}
#endif

/* get geometrical object extents update callback */
BOX_Extents_Update SGP_Extents_Update (SGP *sgp)
{
  switch (sgp->shp->kind)
  {
  case SHAPE_MESH: return (BOX_Extents_Update) ELEMENT_Extents;
  case SHAPE_CONVEX: return (BOX_Extents_Update) CONVEX_Extents;
  case SHAPE_SPHERE: return (BOX_Extents_Update) SPHERE_Extents;
  }

  ASSERT_DEBUG (0, "Invalid shape kind in SGP_Extents_Update");

  return 0;
}
