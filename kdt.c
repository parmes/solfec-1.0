/*
 * kdt.c
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * kd-tree based map
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

#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "kdt.h"
#include "mem.h"
#include "alg.h"
#include "err.h"
#include "hyb.h"

static int compare0 (double **x, double **y)
{
  if (x [0][0] < y [0][0]) return -1;
  else if (x [0][0] > y [0][0]) return 1;
  else return 0;
}

static int compare1 (double **x, double **y)
{
  if (x [0][1] < y [0][1]) return -1;
  else if (x [0][1] > y [0][1]) return 1;
  else return 0;
}

static int compare2 (double **x, double **y)
{
  if (x [0][2] < y [0][2]) return -1;
  else if (x [0][2] > y [0][2]) return 1;
  else return 0;
}

static int (*qcmp [3]) (const void*, const void*) =
{(int (*) (const void*, const void*)) compare0,
 (int (*) (const void*, const void*)) compare1,
 (int (*) (const void*, const void*)) compare2};

static void overlap (void *data, BOX *one, BOX *two)
{
  if (one->sgp < two->sgp) two->mark = (void*) 1; /* mark higher */
  else one->mark = (void*) 1;
}

/* filter out epsilon margin separated points */
static int separate (int n, double **q, double epsilon)
{
  if (epsilon <= 0.0) return n;

  BOX *b, *g, **pb;
  int i, m;

  ERRMEM (b = malloc (n * sizeof (BOX)));
  ERRMEM (pb = malloc (n * sizeof (BOX*)));

  for (i = 0, g = b; i < n; i ++, g ++)
  {
    double *e = g->extents;
    e [0] = q [i][0] - epsilon;
    e [1] = q [i][1] - epsilon;
    e [2] = q [i][2] - epsilon;
    e [3] = e [0] + 2 * epsilon;
    e [4] = e [1] + 2 * epsilon;
    e [5] = e [2] + 2 * epsilon;
    g->sgp = (SGP*) q [i];
    g->mark = NULL;
    pb [i] = g;
  }

  hybrid (pb, n, NULL, overlap);

  for (i = m = 0, g = b; i < n; i ++, g ++)
  {
    if (!g->mark) q [m ++] = (double*) g->sgp;
  }

  free (pb);
  free (b);

  return m;
}

/* split points along most elongated direction */
static int split (int n, double **q, double *p, int *d)
{
  double extents [6] = {q[0][0], q[0][1], q[0][2], q[0][0], q[0][1], q[0][2]};
  double **x, **y;
  int i, k, t [2];

  for (x = q+1, y = q+n; x != y; x ++)
  {
    if (x [0][0] < extents [0]) extents [0] = x [0][0];
    if (x [0][1] < extents [1]) extents [1] = x [0][1];
    if (x [0][2] < extents [2]) extents [2] = x [0][2];
    if (x [0][0] > extents [3]) extents [3] = x [0][0];
    if (x [0][1] > extents [4]) extents [4] = x [0][1];
    if (x [0][2] > extents [5]) extents [5] = x [0][2];
  }

  extents [3] -= extents [0];
  extents [4] -= extents [1];
  extents [5] -= extents [2];

  if (extents [3] >= extents [4] && extents [3] >= extents [5]) *d = 0;
  else if (extents [4] >= extents [3] && extents [4] >= extents [5]) *d = 1;
  else *d = 2;

  for (k = i = 0; k < 3; k ++)
    if (k != *d) t [i ++] = k;
  i = 0;

back:
  qsort (q, n, sizeof (double*), qcmp [*d]);
  k = n / 2;
  x = &q[k];
  while (x < (y-1) && x [0][*d] == (x+1) [0][*d]) { k ++; x ++; } /* q [3*l+d] <= q [3*k+d] for l <= k */
  while (k == n && i < 2) /* no good splitting along current dimension */
  {
    *d = t [i ++];
    goto back;
  }
  if (k == n && i == 2) return -1; /* coincident points */
  COPY (x [0], p);

  return k;
}

/* recursive create */
static KDT* create (int n, double **q)
{
  KDT *kd;
  int k;

  ERRMEM (kd = MEM_CALLOC (sizeof (KDT)));

  if (n == 1)
  {
    COPY (q [0], kd->p);
    kd->d = -1; /* leaf */
    return kd;
  }

  k = split (n, q, kd->p, &kd->d);
  if (k < 0) { free (kd); return NULL; } /* coincident points */
  kd->l = create (k, q);
  kd->r = create (n-k, q+k);

  return kd;
}

/* recursive pick */
static void pick (KDT *kd, double *p, void ***data, int *n)
{
  if (kd->d < 0) /* leaf */
  {
    int k = *n, i;
    (*n) += kd->n;
    ERRMEM (*data = realloc (*data, (*n) * sizeof (void*)));
    for (i = 0; i < kd->n; i ++) (*data) [k+i] = kd->data [i];
  }
  else if (p [kd->d] <= kd->p [kd->d]) pick (kd->l, p, data, n);
  else pick (kd->r, p, data, n);
}

/* create kd-tree for n points separated by epsilon
 * margin; returns NULL for coincident points */
KDT* KDT_Create (int n, double *p, double epsilon)
{
  double **q;
  KDT *kd;
  int i;

  ERRMEM (q = malloc (n * sizeof (double*)));
  for (i = 0; i < n; i ++) q [i] = &p [3*i];
  n = separate (n, q, epsilon);
  kd = create (n, q);
  free (q);

  return kd;
}

/* drop data down the kd-tree */
void KDT_Drop (KDT *kd, double *extents, void *data)
{
  if (kd->d < 0) /* leaf */
  {
    kd->n ++;
    ERRMEM (kd->data = realloc (kd->data, kd->n * sizeof (void*)));
    kd->data [kd->n-1] = data;
  }
  else if (extents [kd->d+3] <= kd->p [kd->d]) KDT_Drop (kd->l, extents, data);
  else if (extents [kd->d] > kd->p [kd->d]) KDT_Drop (kd->r, extents, data);
  else
  {
    KDT_Drop (kd->l, extents, data);
    KDT_Drop (kd->r, extents, data);
  }
}

/* pick data for a point; free buffer after use */
void KDT_Pick (KDT *kd, double *p, void ***data, int *n)
{
  *n = 0;
  *data = NULL;
  pick (kd, p, data, n);
}

/* return nearest point in kd-tree */
double* KDT_Nearest (KDT *kd, double *p)
{
  double a [3], d, dmax, *c;

  dmax = DBL_MAX;
  c = NULL;

  while (kd)
  {
    SUB (p, kd->p, a);
    d = DOT (a, a);
    if (d < dmax)
    {
      dmax = d;
      c = kd->p;
    }

    if (kd->d < 0) break; /* leaf */
    else if (p [kd->d] <= kd->p [kd->d]) kd = kd->l;
    else kd = kd->r;
  }
 
  return c; 
}

/* destroy kd-tree */
void KDT_Destroy (KDT *kd)
{
  if (kd == NULL) return;
  KDT_Destroy (kd->l);
  KDT_Destroy (kd->r);
  free (kd->data);
  free (kd);
}
