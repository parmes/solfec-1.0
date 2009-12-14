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
#if MPI
  int rank;
#endif
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

/* local overlap creation callback => filters our unwnated adjacency */
static void local_create (struct auxdata *aux, BOX *one, BOX *two)
{
  BODY *onebod = one->body, *twobod = two->body;

#if MPI
  SET *item = SET_First (one->ranks), *jtem = SET_First (two->ranks);

  while (item->data != jtem->data) /* find the lowest common rank of two boxes */
  {
    if (item->data < jtem->data) item = SET_Next (item);
    else jtem = SET_Next (jtem);

    ASSERT_DEBUG (item && jtem, "Inconsistent lowest rank search for a box pair");
  }

  if (aux->rank > (int) (long) item->data) return; /* filter out overlaps from oter than the lowest common rank of two boxes */
#endif

  /* check if these are two obstacles => no need to report overlap */
  if (onebod->kind == OBS && twobod->kind == OBS) return;

  /* self-contact ? <=> one->body == two->body */
  if (onebod == twobod && !(onebod->flags & BODY_DETECT_SELF_CONTACT)) return;

  int id1 = onebod->id,
      id2 = twobod->id,
      no1 = one->sgp - onebod->sgp,
      no2 = two->sgp - twobod->sgp;

  OPR pair = {MIN (id1, id2), MAX (id1, id2), MIN (no1, no2), MAX (no1, no2)};

  /* an excluded body pair ? */
  if (SET_Contains (aux->nobody, &pair, (SET_Compare) bodcmp)) return;

  /* an excluded object pair ? */
  if (SET_Contains (aux->nogobj, &pair, (SET_Compare) gobjcmp)) return;

  /* test topological adjacency */
  switch (GOBJ_Pair_Code (one, two))
  {
    case AABB_ELEMENT_ELEMENT:
    {
      if (onebod == twobod) /* assuming only one mesh per body is possible,
			       the two elements are within the same mesh;
			       exclude topologically adjacent elements */
	if (ELEMENT_Adjacent (one->sgp->gobj, two->sgp->gobj)) return;
    }
    break;
    case AABB_CONVEX_CONVEX:
    {
      if (onebod == twobod) /* exclude topologically adjacent convices */
	if (CONVEX_Adjacent (one->sgp->gobj, two->sgp->gobj)) return;
    }
    break;
    case AABB_SPHERE_SPHERE:
    {
      if (onebod == twobod) /* exclude topologically adjacent spheres */
	if (SPHERE_Adjacent (one->sgp->gobj, two->sgp->gobj)) return;
    }
  }
 
  /* report overlap creation */
  aux->create (aux->data, one, two);
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
/* detach boxes from outside of the domain and attach new incoming boxes */
static void detach_and_attach (AABB *aabb)
{
  int *procs, numprocs, rank;
  BOX_Extents_Update update;
  SET *delbox, *item;
  SGP *sgp, *sge;
  double e [6];
  int j, ncpu;
  BODY *bod;
  BOX *box;
  DOM *dom;

  delbox = NULL;
  dom = aabb->dom;
  rank = dom->rank;
  ncpu = dom->ncpu;
  ERRMEM (procs = malloc (sizeof (int [ncpu])));

  /* test existing boxes for overlap with this partition */
  for (box = aabb->lst; box; box = box->next)
  {
    COPY6 (box->extents, e);

    ASSERT (Zoltan_LB_Box_Assign (dom->zol, e[0], e[1], e[2], e[3], e[4], e[5], procs, &numprocs) == ZOLTAN_OK, ERR_ZOLTAN);

    for (j = 0; j < numprocs; j ++)
    {
      if (procs [j] == rank) break;
    }

    if (j == numprocs) SET_Insert (&aabb->setmem, &delbox, box, NULL); /* schedule for deletion */
    else
    {
      SET_Free (&aabb->setmem, &box->ranks);
      for (j = 0; j < numprocs; j ++) SET_Insert (&aabb->setmem, &box->ranks, (void*) (long) procs [j], NULL); /* refresh box ranks */
    }
  }

  /* attach parents */
  for (bod = dom->bod; bod; bod = bod->next)
  {
    for (sgp = bod->sgp, sge = sgp + bod->nsgp; sgp < sge; sgp ++)
    {
      if (sgp->box == NULL)
      {
	update = SGP_Extents_Update (sgp);

	update (sgp->shp->data, sgp->gobj, e);

	ASSERT (Zoltan_LB_Box_Assign (dom->zol, e[0], e[1], e[2], e[3], e[4], e[5], procs, &numprocs) == ZOLTAN_OK, ERR_ZOLTAN);

	for (j = 0; j < numprocs; j ++)
	{
	  if (procs [j] == rank) break;
	}

	if (j < numprocs)
	{
	  box = AABB_Insert (aabb, bod, gobj_kind (sgp), sgp, update);

	  COPY6 (e, box->extents);

	  for (j = 0; j < numprocs; j ++) SET_Insert (&aabb->setmem, &box->ranks, (void*) (long) procs [j], NULL); /* set box ranks */
	}
      }
    }
  }

  /* attach children */
  for (item = SET_First (dom->children); item; item = SET_Next (item))
  {
    for (bod = item->data, sgp = bod->sgp, sge = sgp + bod->nsgp; sgp < sge; sgp ++)
    {
      if (sgp->box == NULL)
      {
	update = SGP_Extents_Update (sgp);

	update (sgp->shp->data, sgp->gobj, e);

	ASSERT (Zoltan_LB_Box_Assign (dom->zol, e[0], e[1], e[2], e[3], e[4], e[5], procs, &numprocs) == ZOLTAN_OK, ERR_ZOLTAN);

	for (j = 0; j < numprocs; j ++)
	{
	  if (procs [j] == rank) break;
	}

	if (j < numprocs)
	{
	  box = AABB_Insert (aabb, bod, gobj_kind (sgp), sgp, update);

	  COPY6 (e, box->extents);

	  for (j = 0; j < numprocs; j ++) SET_Insert (&aabb->setmem, &box->ranks, (void*) (long) procs [j], NULL); /* set box ranks */
	}
      }
    }
  }

  /* delete exported boxes */
  for (item = SET_First (delbox); item; item = SET_Next (item))
  {
    AABB_Delete (aabb, item->data);
  }

  /* clean */
  SET_Free (&aabb->setmem, &delbox);
  free (procs);
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

  return aabb;
}

/* insert geometrical object */
BOX* AABB_Insert (AABB *aabb, BODY *body, GOBJ kind, SGP *sgp, BOX_Extents_Update update)
{
  BOX *box;

  ERRMEM (box = MEM_Alloc (&aabb->boxmem));
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
    AABB_Insert (aabb, body, gobj_kind (sgp), sgp, SGP_Extents_Update (sgp));
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
void AABB_Update (AABB *aabb, BOXALG alg, void *data, BOX_Overlap_Create create)
{
#if MPI
  struct auxdata aux = {aabb->dom->rank, aabb->nobody, aabb->nogobj, &aabb->mapmem, data, create};
#else
  struct auxdata aux = {aabb->nobody, aabb->nogobj, &aabb->mapmem, data, create};
#endif
  BOX *box;

#if MPI
  if (aabb->dom->rank == 0)
#endif
  if (aabb->dom && aabb->dom->verbose) printf ("CONDET (%s) ... ", AABB_Algorithm_Name (alg)), fflush (stdout);
  
  if (aabb->dom) SOLFEC_Timer_Start (aabb->dom->solfec, "CONDET");

  for (box = aabb->lst; box; box = box->next) box->update (box->sgp->shp->data, box->sgp->gobj, box->extents); /* update box extents */

#if MPI
  detach_and_attach (aabb); /* detach boxes from outside of the domain and attach new incoming boxes */
#endif

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

  free (aabb);
}

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
