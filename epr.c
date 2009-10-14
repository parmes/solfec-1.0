/*
 * epr.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * extended pseudo-rigid model
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

#include "epr.h"
#include "hyb.h"
#include "goc.h"
#include "alg.h"
#include "err.h"

typedef struct epr_element EPR_ELEMENT;
typedef struct epr_mesh EPR_MESH;

/* extended pseudo-rigid element */
struct epr_element
{
  SHAPE *head, *tail; /* beginning and end of shapes within this element */

  EPR_ELEMENT **adj;

  int nadj;

  double A [3]; /* a referential auxiliary point */
};

/* extended pseudo-rigid mesh */
struct epr_mesh
{
  EPR_ELEMENT *ele;

  int nele;
};

/* overlap callback for element adjacency */
static void* overlap (void *data, BOX *one, BOX *two)
{
  double p [3], q [3];

  if (gobjdistance (GOBJ_Pair_Code (one, two), one->sgp, two->sgp, p, q) < GEOMETRIC_EPSILON) /* if they touch */
  {
    EPR_ELEMENT *epr = one->sgp->shp->epr,
		*epq = two->sgp->shp->epr, **x, **y;

    for (x = epr->adj, y = x + epr->nadj; x < y; x ++)
    {
      if (*x == epq) return NULL; /* already adjacent */
    }

    ERRMEM (epr->adj = realloc (epr->adj, (++epr->nadj) * sizeof (EPR_ELEMENT*)));  /* extend adjacency */
    epr->adj [epr->nadj-1] = epq;
    ERRMEM (epq->adj = realloc (epq->adj, (++epq->nadj) * sizeof (EPR_ELEMENT*))); 
    epq->adj [epq->nadj-1] = epr;
  }

  return NULL;
}

/* set up adjacency pointers */
void set_up_adjacency (SHAPE *shp)
{
  BOX_Extents_Update update;
  SGP *sgp, *sge, *sgx;
  BOX **boxes;
  MEM mem;
  int num;

  sgp = SGP_Create (shp, &num);
  MEM_Init (&mem, sizeof (BOX), num);
  ERRMEM (boxes = malloc (sizeof (AABB*) * num));
  for (sgx = sgp, sge = sgp + num, num = 0; sgx < sge; sgx ++, num ++)
  {
    ERRMEM (boxes [num] = MEM_Alloc (&mem));
    update = SGP_Extents_Update (sgx);
    update (sgx->shp->data, sgx->gobj, boxes [num]->extents); /* set up extents */
    boxes [num]->sgp = sgx;
  }

  hybrid (boxes, num, NULL, overlap); /* detect boxoverlaps => set adjacency inside the callback */

  MEM_Release (&mem); /* done */
  free (boxes);
}

/* create EPR internals for a body */
void EPR_Create (SHAPE *shp, BULK_MATERIAL *mat, BODY *bod)
{
  EPR_ELEMENT *ele;
  EPR_MESH *msh;
  SHAPE *shq;

  ERRMEM (msh = malloc (sizeof (EPR_MESH)));

  for (shq = shp, msh->nele = 1; shq; shq = shq->next) /* initialize with 1 for the sake of the last NULL-terimed element */
  {
    if (shq->epr) msh->nele ++; /* extended element list end marker */
  }

  ERRMEM (msh->ele = calloc (msh->nele, sizeof (EPR_ELEMENT)));

  for (shq = shp, ele = msh->ele, ele->head = shq; shq; shq = shq->next)
  {
    shq->epr = ele; /* current shape points to current element */

    if (shq->epr) /* end of element list marker (inclusive) */
    {
      ele->tail = shq->next; /* record list tail (exclusive) */

      ele ++; /* iterate to next element */

      if (shq->next) ele->head = shq->next; /* record list head */
    }
  }

  set_up_adjacency (shp);

  bod->dofs = 9 * msh->nele + 3;

  ERRMEM (bod->conf = calloc (2 * bod->dofs, sizeof (double)));
  bod->velo = bod->conf + bod->dofs;

  for (double *F = bod->conf, *E = F + 9 * msh->nele; F < E; F += 9) { IDENTITY (F); }

  bod->priv = msh;
}

/* return configuration size */
int EPR_Conf_Size (BODY *bod)
{
  return bod->dofs;
}

/* overwrite state */
void EPR_Overwrite_State (BODY *bod, double *q, double *u)
{
  double *x, *y;

  for (x = bod->conf, y = x + bod->dofs; x < y; x ++, q ++) *x = *q;
  for (x = bod->velo, y = x + bod->dofs; x < y; x ++, u ++) *x = *u;
}

/* set initial rigid motion velocity */
void EPR_Initial_Velocity (BODY *bod, double *linear, double *angular)
{
  double *u = bod->velo, *e = u + bod->dofs - 3;

  if (angular)
  {
    for (; u < e; u += 9) { VECSKEW (angular, u); }
  }

  if (linear) { COPY (linear, e); }
}

/* initialise dynamic time stepping */
void EPR_Dynamic_Init (BODY *bod, SCHEME scheme)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* estimate critical step for the dynamic scheme */
double EPR_Dynamic_Critical_Step (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
  return 0.0;
}

/* perform the initial half-step of the dynamic scheme */
void EPR_Dynamic_Step_Begin (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* perform the final half-step of the dynamic scheme */
void EPR_Dynamic_Step_End (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* initialise static time stepping */
void EPR_Static_Init (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* perform the initial half-step of the static scheme */
void EPR_Static_Step_Begin (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* perform the final half-step of the static scheme */
void EPR_Static_Step_End (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* motion x = x (X, state) */
void EPR_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* inverse motion X = X (x, state) */
void EPR_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* obtain velocity at (element, point), expressed in the local 'base' => all entities are spatial */
void EPR_Local_Velo (BODY *bod, VELOTIME time, SHAPE *shp, void *gobj, double *point, double *base, double *velo)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* return transformation operator from the generalised to the local velocity space at (element, point, base) */
MX* EPR_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
  return NULL;
}

/* compute current kinetic energy */
double EPR_Kinetic_Energy (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
  return 0.0;
}

/* get some values at a node of a geometrical object */
void EPR_Nodal_Values (BODY *bod, SHAPE *shp, void *gobj, int node, VALUE_KIND kind, double *values)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* get some values at a referential point */
void EPR_Point_Values (BODY *bod, double *point, VALUE_KIND kind, double *values)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* release EPR memory */
void EPR_Destroy (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}
