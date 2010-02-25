/*
 * hyb.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------------------
 * box overlap detection using the hybrid streamed segment tree described in:
 * Afran Zomorodian, Herbert Edelsbrunner, "Fast software for box intersection",
 * 2002, International Journal of Computer Geometry and Applications, Vol. 12,
 * Numbers 1-2, pp. 142-172.
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
#include <string.h>
#include <float.h>
#include <math.h>
#include "err.h"
#include "hyb.h"

#define CUTOFF 1024
#define PLT(b1, b2, d) ((b1)->extents[d] < (b2)->extents[d] ? 1 : ((b1)->extents[d] == (b2)->extents[d] && (b1)->sgp < (b2)->sgp ? 1 : 0))
typedef int (*QCMP) (const void*, const void*); /* qsort comparison type */

/* compare for qsort */
static int boxcmp (BOX **a, BOX **b)
{
  if (PLT (*a,*b,0)) return -1;
  else return 1;
}

/* scan intervals with points along dimensions */
static void onewayscan (BOX **Ib, BOX **Ie, BOX **Pb, BOX **Pe,
                        int d, void *data, BOX_Overlap_Create create)
{
  BOX **P;
  int n;

  /* sort intervals and points along dimension 0 */
  qsort (Ib, Ie-Ib, sizeof (BOX*), (QCMP)boxcmp);
  qsort (Pb, Pe-Pb, sizeof (BOX*), (QCMP)boxcmp);

  for (; Ib < Ie; Ib ++)
  {
    for (; Pb < Pe && PLT (*Pb, *Ib, 0); Pb ++); /* skip points before 'Ib' */

    for (P = Pb; P < Pe && (*P)->extents [0] < (*Ib)->extents [3+d]; P ++) /* until 'P' is inside 'Ib' */
    {
      for (n = 1; n <= d; n ++) /* zero overlaps already, check 1,...,d */
      {
	if ((*Ib)->extents [n] > (*P)->extents [3+n] ||
	    (*Ib)->extents [3+n] < (*P)->extents [n]) break; /* break if 'P' does not overlap 'Ib' along dimension 'n' */
      }

      if (n > d && (*Ib)->sgp != (*P)->sgp) create (data, *Ib, *P); /* report overlap, if (P, Ib) overlap along all dimensions and P != Ib */
    }
  }
}

/* scan interchanging roles of points and intervals */
static void twowayscan (BOX **Ib, BOX **Ie, BOX **Pb, BOX **Pe,
                        int d, void *data, BOX_Overlap_Create create)
{
  BOX **I, **i, **P, **p;
  int n;

  /* sort intervals and points along dimension 0 */
  qsort (Ib, Ie-Ib, sizeof (BOX*), (QCMP)boxcmp);
  qsort (Pb, Pe-Pb, sizeof (BOX*), (QCMP)boxcmp);

  while (Ib < Ie && Pb < Pe)
  {
    if (PLT (*Ib, *Pb, 0)) /* 'Ib' begins before 'Pb' => 'Ib' is the interval */
    {
      I = Ib ++; /* iterate 'Ib', while 'I' is the current interval */
      for (p = Pb; p < Pe && (*p)->extents [0] < (*I)->extents [3]; p ++) /* for all points 'p' inside 'I' */
      {
	for (n = 1; n <= d; n ++) /* do the remaining check ... (see onewayscan) */
	{
	  if ((*I)->extents [n] > (*p)->extents [3+n] ||
	      (*I)->extents [3+n] < (*p)->extents [n]) break;
	}

	if (n > d && (*I)->sgp != (*p)->sgp) create (data, *I, *p);
      }
    }
    else /* here, 'Pb' is the interval */
    {
      P = Pb ++; /* iterate 'Pb', while 'P' is the current interval */
      for (i = Ib; i < Ie && (*i)->extents [0] < (*P)->extents [3]; i ++) /* for all points 'i' inside 'P' */
      {
	for (n = 1; n <= d; n ++) /* as above ... */
	{
	  if ((*P)->extents [n] > (*i)->extents [3+n] ||
	      (*P)->extents [3+n] < (*i)->extents [n]) break;
	}

	if (n > d && (*P)->sgp != (*i)->sgp) create (data, *P, *i);
      }
    }
  }
}

/* [Ib, lo_hi_inside) collects d-intervals containig [lo, hi) */
static BOX** lo_hi_inside (BOX **Ib, BOX **Ie, double lo, double hi, int d)
{
  BOX **i, **j, *k;

  for (i = Ib, j = Ie; i < j; )
  {
    while (i < j && (*i)->extents [d] <= lo && (*i)->extents [3+d] >= hi) i ++;
    do { j --; } while (i < j && ((*j)->extents [d] > lo || (*j)->extents [3+d] < hi));
    if (i < j)
    {
       k = *i;
      *i = *j;
      *j =  k;
    }
  }

  return i;
}

/* ternary tree height */
inline static int height (int n)
{ return (int) (1.0 + 0.5 * log ((double)n) / 1.0986122); /* log (n) / log (3) */ }

/* median of three numbers (lower extents coords along 'd') */
static BOX* median_of_three (BOX *a, BOX *b, BOX *c, int d)
{
  if (PLT (a, b, d))
  {
    if (PLT (c, a, d)) return a;
    else if (PLT (c, b, d)) return c;
    else return b;
  }
  else
  {
    if (PLT (c, b, d)) return b;
    else if (PLT (c, a, d)) return c;
    else return a;
  }
}

/* approximate median of a point set */
static BOX* median (BOX **Pb, BOX **Pe, int d, int h)
{
  if (h == 0) return *(Pb + rand () % (Pe-Pb));
  else return median_of_three (median (Pb, Pe, d, h-1),
                               median (Pb, Pe, d, h-1),
                               median (Pb, Pe, d, h-1), d);
}

/* [Pb, split) are in [-inf, mi); [split, Pe) are in [mi, +inf) */
static BOX** split (BOX **Pb, BOX **Pe, double mi, int d)
{
  BOX **i, **j, *k;

  for (i = Pb, j = Pe; i < j; )
  {
    while (i < j && (*i)->extents [d] <  mi) i ++; /* skip points already in [-inf, mi) */
    do { j --; } while (i < j && (*j)->extents [d] >= mi); /* skip points already in [mi, +inf) */
    if (i < j) /* swap points if the two ends did not meet */
    {
       k = *i;
      *i = *j;
      *j =  k;
    }
  }

  return i; /* the end of the lower range */
}

/* intervals [Ib, overlaps) overlap [lo, hi) */
static BOX** overlaps (BOX **Ib, BOX **Ie, double lo, double hi, int d)
{
  BOX **i, **j, *k;

  for (i = Ib, j = Ie; i < j; )
  {
    while (i < j && !((*i)->extents [d] >= hi || (*i)->extents [3+d] < lo)) i ++;
    do { j --; } while (i < j && ((*j)->extents [d] >= hi || (*j)->extents [3+d] < lo));
    if (i < j)
    {
       k = *i;
      *i = *j;
      *j =  k;
    }
  }

  return i;
}

/* streamed segment tree */
static void stream (BOX **Ib, BOX **Ie, BOX **Pb, BOX **Pe,
  double lo, double hi, int d, void *data, BOX_Overlap_Create create)
{
  if (Ib >= Ie || Pb >= Pe) return;
  else if (d == 0) onewayscan (Ib, Ie, Pb, Pe, d, data, create);
  else if ((Ie-Ib) < CUTOFF || (Pe-Pb) < CUTOFF) twowayscan (Ib, Ie, Pb, Pe, d, data, create);
  else
  {
    BOX **Im, **Pm, *P;
    double mi;

    Im = lo_hi_inside (Ib, Ie, lo, hi, d); /* [Ib, Im) collects intervals containig [lo, hi)
                                              [Im, Ie) enumerates the remaining intervals */

    /* recurse along lower dimensions */
    stream (Ib, Im, Pb, Pe, -DBL_MAX, DBL_MAX, d-1, data, create);
    stream (Pb, Pe, Ib, Im, -DBL_MAX, DBL_MAX, d-1, data, create);

    /* continue down the tree
     * along 'd'th dimension */

    P = median (Pb, Pe, d, height (Pe-Pb)); /* approximate median of points */
    mi = P->extents [d];
    
    Pm = split (Pb, Pe, mi, d); /* [Pb, Pm) are in [lo, mi); [Pm, Pe) are in [mi, hi) */

    Im = overlaps (Ib, Ie, lo, mi, d); /* intervals [Ib, Im) overlap [lo, mi) */
    stream (Ib, Im, Pb, Pm, lo, mi, d, data, create);

    Im = overlaps (Ib, Ie, mi, hi, d); /* intervals [Ib, Im) overlap [mi, hi) */
    stream (Ib, Im, Pm, Pe, mi, hi, d, data, create);
  }
}

/* streamed segment tree without twowayscan */
static void stream_ext (BOX **Ib, BOX **Ie, BOX **Pb, BOX **Pe,
  double lo, double hi, int d, void *data, BOX_Overlap_Create create)
{
  if (Ib >= Ie || Pb >= Pe) return;
  else if (d == 0 || (Ie-Ib) < CUTOFF || (Pe-Pb) < CUTOFF) onewayscan (Ib, Ie, Pb, Pe, d, data, create);
  else
  {
    BOX **Im, **Pm, *P;
    double mi;

    Im = lo_hi_inside (Ib, Ie, lo, hi, d); /* [Ib, Im) collects intervals containig [lo, hi)
                                              [Im, Ie) enumerates the remaining intervals */

    /* recurse along lower dimensions */
    stream_ext (Ib, Im, Pb, Pe, -DBL_MAX, DBL_MAX, d-1, data, create);
    stream_ext (Pb, Pe, Ib, Im, -DBL_MAX, DBL_MAX, d-1, data, create);

    /* continue down the tree
     * along 'd'th dimension */

    P = median (Pb, Pe, d, height (Pe-Pb)); /* approximate median of points */
    mi = P->extents [d];
    
    Pm = split (Pb, Pe, mi, d); /* [Pb, Pm) are in [lo, mi); [Pm, Pe) are in [mi, hi) */

    Im = overlaps (Ib, Ie, lo, mi, d); /* intervals [Ib, Im) overlap [lo, mi) */
    stream_ext (Ib, Im, Pb, Pm, lo, mi, d, data, create);

    Im = overlaps (Ib, Ie, mi, hi, d); /* intervals [Ib, Im) overlap [mi, hi) */
    stream_ext (Ib, Im, Pm, Pe, mi, hi, d, data, create);
  }
}

/* hybrid overlap detection driver */
void hybrid (BOX **boxes, int n, void *data, BOX_Overlap_Create create)
{
  BOX **copy;

  ERRMEM (copy = malloc (sizeof (BOX*) * n));
  memcpy (copy, boxes, sizeof (BOX*) * n); /* copy of the pointers table */
  
  stream (boxes, boxes + n, /* these are intervals */
          copy, copy + n, /* these are points */
	 -DBL_MAX, DBL_MAX, /* the top level interval is [-inf, inf] */
	  2, data, create); /* we go from the thrid (2) dimension down to one (0) */

  free (copy);
}

/* report overlaps between two sets of boxes */
void hybrid_ext (BOX **seta, int na, BOX **setb, int nb, void *data, BOX_Overlap_Create create)
{
  stream_ext (seta, seta + na, /* these are intervals */
              setb, setb + nb, /* these are points */
	     -DBL_MAX, DBL_MAX, /* the top level interval is [-inf, inf] */
	      2, data, create); /* we go from the thrid (2) dimension down to one (0) */

  stream_ext (setb, setb + nb,
              seta, seta + na,
	     -DBL_MAX, DBL_MAX,
	      2, data, create);
}
