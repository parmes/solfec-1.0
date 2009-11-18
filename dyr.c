/*
 * dyr.c
 * Copyright (C) 2005, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------------
 * dynamic rectangle structure
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


#include <stdlib.h>
#include <stdio.h>
#include "alg.h"
#include "mem.h"
#include "xyt.h"
#include "dyr.h"
#include "err.h"

#define MEMCHUNK 512

typedef struct rect RECT;
typedef struct range RANGE;
typedef struct xyaux XYAUX;
typedef struct dynrect DR;

struct rect
{
  double xmin, ymin, xmax, ymax;
  RECT *next;
  BOX *box;
};

struct range
{
  double min, max;
  BOX *box;
};

struct xyaux
{
  double min, max;
  BOX *box;
  void *data;
  BOX_Overlap_Create report;
};

struct dynrect
{
  int hsize;
  DRALG algo;
  double avglen [2];

  MEM rectpool;
  MEM xytreepool;
  MEM intpool;
  MEM rngpool;

  BOX **hashmarks;

  /* HASH2D_LIST */
  RECT **hash_list;

  /* HASH1D_XYTREE */
  /* HASH2D_XYTREE */
  struct xytree  **hash_xytree;
  
  /* XYTREE */
  struct xytree *just_xytree;

  double *rkey;
};

#define xyaux_init(aux_, min_, max_, box_, data_, report_) \
  aux_.min = min_; \
  aux_.max = max_; \
  aux_.box = box_; \
  aux_.data = data_; \
  aux_.report = report_

inline static int overlap (RECT *r, double xmin, double ymin, double xmax, double ymax)
{
  if (r->xmax < xmin || r->xmin > xmax) return 0;
  else if (r->ymax < ymin || r->ymin > ymax) return 0;
  else return 1;
}

static void hash2d_list_insert (DR *dyn, double xmin, double ymin, double xmax,
                                double ymax, BOX *box, void *data, BOX_Overlap_Create report)
{
  BOX **hmarks = dyn->hashmarks;
  int hsize = dyn->hsize;
  int i, j, imin, jmin, imax, jmax, key;
  RECT **htable = dyn->hash_list, *r;
  MEM *pool = &dyn->rectpool;
  double *avglen = dyn->avglen;
  char found;

  imin = INT (xmin, avglen [0]);
  jmin = INT (ymin, avglen [1]);
  imax = INT (xmax, avglen [0]);
  jmax = INT (ymax, avglen [1]);

  for (i = imin; i <= imax; i ++)
  for (j = jmin; j <= jmax; j ++)
  {
    key = HASH2 (i, j, hsize);

    if (hmarks [key] == box) continue;

    found = 0;

    for (r = htable [key]; r; r = r->next)
    {
      if (r->box == box) found = 1;
      else if (r->box->mark != box->sgp && overlap (r, xmin, ymin, xmax, ymax))
      {
        report (data, box, r->box);
	r->box->mark = box->sgp;
      }
    }

    if (found == 0)
    {
      ERRMEM (r = MEM_Alloc (pool));
      r->xmin = xmin; r->ymin = ymin; r->xmax = xmax; r->ymax = ymax;
      r->box = box;
      r->box->mark = NULL;
      r->next = htable [key];
      htable [key] = r;
    }

    hmarks [key] = box;
  }
}

static void xycallback (XYAUX *aux, RANGE *r)
{
  if (r->box->mark == aux->box->sgp) return;
  else if (aux->max < r->min || aux->min > r->max) return;

  aux->report (aux->data, aux->box, r->box);
  r->box->mark = aux->box->sgp;
}

static void hash1d_xytree_insert (DR *dyn, double xmin, double ymin, double xmax,
                                  double ymax, BOX *box, void *data, BOX_Overlap_Create report)
{
  BOX **hmarks = dyn->hashmarks;
  int hsize = dyn->hsize, i, imin, imax, key;
  XYAUX aux;
  struct xytree **htable = dyn->hash_xytree;
  MEM *pool = &dyn->xytreepool;
  double *avglen = dyn->avglen;
  RANGE *r;

  xyaux_init (aux, xmin, xmax, box, data, report);

  imin = INT (xmin, avglen [0]);
  imax = INT (xmax, avglen [0]);

  r = MEM_Alloc (&dyn->rngpool);
  r->min = xmin;
  r->max = xmax;
  r->box = box;
  box->mark = NULL;

  for (i = imin; i <= imax; i ++)
  {
    key = HASH1 (i, hsize);

    if (hmarks [key] == box) continue;

    XY_Query (htable [key], ymin, -ymax, &aux, (XY_Callback) xycallback);

    XY_Insert (pool, &htable [key], ymax, -ymin, box, r);

    hmarks [key] = box;
  }
}

static void hash2d_xytree_insert (DR *dyn, double xmin, double ymin, double xmax,
                                  double ymax, BOX *box, void *data, BOX_Overlap_Create report)
{
  BOX **hmarks = dyn->hashmarks;
  int hsize = dyn->hsize, i, j, imin, jmin, imax, jmax, key;
  XYAUX aux;
  struct xytree **htable = dyn->hash_xytree;
  MEM *pool = &dyn->xytreepool;
  double *avglen = dyn->avglen;
  RANGE *r;

  xyaux_init (aux, ymin, ymax, box, data, report);

  imin = INT (xmin, avglen [0]);
  jmin = INT (ymin, avglen [1]);
  imax = INT (xmax, avglen [0]);
  jmax = INT (ymax, avglen [0]);

  r = MEM_Alloc (&dyn->rngpool);
  r->min = ymin;
  r->max = ymax;
  r->box = box;
  box->mark = NULL;

  for (i = imin; i <= imax; i ++)
  for (j = jmin; j <= jmax; j ++)
  {
    key = HASH2 (i, j, hsize);

    if (hmarks [key] == box) continue;

    XY_Query (htable [key], xmin, -xmax, &aux, (XY_Callback) xycallback);

    XY_Insert (pool, &htable [key], xmax, -xmin, box, r);

    hmarks [key] = box;
  }
}

static void xytree_insert (DR *dyn, double xmin, double ymin, double xmax,
                           double ymax, BOX *box, void *data, BOX_Overlap_Create report)
{
  XYAUX aux;
  MEM *pool = &dyn->xytreepool;
  struct xytree **tree = &dyn->just_xytree;
  RANGE *r;

  xyaux_init (aux, ymin, ymax, box, data, report);

  XY_Query (*tree, xmin, -xmax, &aux, (XY_Callback) xycallback);

  r = MEM_Alloc (&dyn->rngpool);
  r->min = ymin;
  r->max = ymax;
  r->box = box;
  box->mark = NULL;

  XY_Insert (pool, tree, xmax, -xmin, box, r);
}

static void (*insert [HASH1D_XYTREE+1])
  (DR*, double, double, double, double, BOX*, void*, BOX_Overlap_Create)
  = { hash2d_list_insert, hash2d_xytree_insert, xytree_insert, hash1d_xytree_insert };

static void hash2d_list_delete (DR *dyn, double xmin, double ymin, double xmax, double ymax, BOX *box)
{
  BOX **hmarks = dyn->hashmarks;
  int hsize = dyn->hsize;
  int i, j, imin, jmin, imax, jmax, key;
  RECT **htable = dyn->hash_list, *r, *prev;
  MEM *pool = &dyn->rectpool;
  double *avglen = dyn->avglen;

  imin = INT (xmin, avglen [0]);
  jmin = INT (ymin, avglen [1]);
  imax = INT (xmax, avglen [0]);
  jmax = INT (ymax, avglen [1]);

  for (i = imin; i <= imax; i ++)
  for (j = jmin; j <= jmax; j ++)
  {
    key = HASH2 (i, j, hsize);
    hmarks [key] = NULL;
  } 

  for (i = imin; i <= imax; i ++)
  for (j = jmin; j <= jmax; j ++)
  {
    key = HASH2 (i, j, hsize);

    if (hmarks [key] == box) continue;


    for (prev = NULL, r = htable [key]; r; prev = r, r = r->next)
    {
      if (r->box == box)
      {
        if (prev) prev->next = r->next;
	else htable [key] = r->next;

	MEM_Free (pool, r);
	break;
      }
    }

    hmarks [key] = box;
  }
}

static void hash1d_xytree_delete (DR *dyn, double xmin, double ymin, double xmax, double ymax, BOX *box)
{
  BOX **hmarks = dyn->hashmarks;
  int hsize = dyn->hsize, i, imin, imax, key;
  struct xytree **htable = dyn->hash_xytree;
  MEM *pool = &dyn->xytreepool;
  double *avglen = dyn->avglen;

  imin = INT (xmin, avglen [0]);
  imax = INT (xmax, avglen [0]);

  for (i = imin; i <= imax; i ++)
  {
    key = HASH1 (i, hsize);
    hmarks [key] = NULL;
  }

  for (i = imin; i <= imax; i ++)
  {
    key = HASH1 (i, hsize);

    if (hmarks [key] == box) continue;

    XY_Delete (pool, &htable [key], ymax, -ymin, box);

    hmarks [key] = box;
  }
}

static void hash2d_xytree_delete (DR *dyn, double xmin, double ymin, double xmax, double ymax, BOX *box)
{
  BOX **hmarks = dyn->hashmarks;
  int hsize = dyn->hsize, i, j, imin, jmin, imax, jmax, key;
  struct xytree **htable = dyn->hash_xytree;
  MEM *pool = &dyn->xytreepool;
  double *avglen = dyn->avglen;

  imin = INT (xmin, avglen [0]);
  jmin = INT (ymin, avglen [1]);
  imax = INT (xmax, avglen [0]);
  jmax = INT (ymax, avglen [0]);


  for (i = imin; i <= imax; i ++)
  for (j = jmin; j <= jmax; j ++)
  {
    key = HASH2 (i, j, hsize);
    hmarks [key] = NULL;
  } 

  for (i = imin; i <= imax; i ++)
  for (j = jmin; j <= jmax; j ++)
  {
    key = HASH2 (i, j, hsize);

    if (hmarks [key] == box) continue;

    XY_Delete (pool, &htable [key], xmax, -xmin, box);

    hmarks [key] = box;
  }
}

static void xytree_delete (DR *dyn, double xmin, double ymin, double xmax, double ymax, BOX *box)
{
  MEM *pool = &dyn->xytreepool;
  struct xytree **tree = &dyn->just_xytree;

  XY_Delete (pool, tree, xmax, -xmin, box);
}

static void (*delete [HASH1D_XYTREE+1]) (DR*, double, double, double, double, BOX*)
  = { hash2d_list_delete, hash2d_xytree_delete, xytree_delete, hash1d_xytree_delete };

/* Create dynamic rectangle (DR) structure and return 'dr' pointer */
void* DR_Create (int boxnum, DRALG algo)
{
  DR *d;

  ERRMEM (d = malloc (sizeof (DR)));

  switch (algo)
  {
    case HASH2D_LIST:
      ERRMEM (d->hash_list = MEM_CALLOC (sizeof (RECT*) * boxnum));
      MEM_Init (&d->rectpool, sizeof (RECT), MEMCHUNK);
      break;
    case HASH1D_XYTREE:
    case HASH2D_XYTREE:
      ERRMEM (d->hash_xytree = MEM_CALLOC (sizeof (struct xytree*) * boxnum));
      MEM_Init (&d->xytreepool, sizeof (struct xytree), MEMCHUNK);
      MEM_Init (&d->rngpool, sizeof (RANGE), boxnum);
      break;
    case XYTREE_ALONE:
      d->just_xytree = NULL;
      MEM_Init (&d->xytreepool, sizeof (struct xytree), MEMCHUNK);
      MEM_Init (&d->rngpool, sizeof (RANGE), boxnum);
      break;
  }

  d->hsize = boxnum;
  d->algo = algo;

  ERRMEM (d->hashmarks = MEM_CALLOC (sizeof (BOX*) * d->hsize));

  return d;
}

/* Set up avarage edge lengths for hashing approaches */
void DR_Params (void *dr, double avglen [2])
{
  DR *d = dr;

  d->avglen [0] = avglen [0];
  d->avglen [1] = avglen [1];
}

/* Insert a rectnagle into DR, performing the overlap query before */
void DR_Insert (void *dr, double xmin, double ymin, double xmax,
                double ymax, BOX *box, void *data, BOX_Overlap_Create report)
{
  DR *d = dr;

  insert [d->algo] (d, xmin, ymin, xmax, ymax, box, data, report);
}

/* Delete a rectangle from the DR */
void DR_Delete (void *dr, double xmin, double ymin, double xmax, double ymax, BOX *box)
{
  DR *d = dr;

  delete [d->algo] (d, xmin, ymin, xmax, ymax, box);
}

/* Clean up */
void DR_Destroy (void *dr)
{
  DR *d = dr;

  switch (d->algo)
  {
    case HASH2D_LIST:
      free (d->hash_list);
      MEM_Release (&d->rectpool);
      break;
    case HASH1D_XYTREE:
    case HASH2D_XYTREE:
      free (d->hash_xytree);
      MEM_Release (&d->xytreepool);
      MEM_Release (&d->rngpool);
      break;
    case XYTREE_ALONE:
      MEM_Release (&d->xytreepool);
      MEM_Release (&d->rngpool);
      break;
  }

  free (d->hashmarks);
  free (d);
}
