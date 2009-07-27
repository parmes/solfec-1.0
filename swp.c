/*
 * swp.c
 * Copyright (C) 2005, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------------
 * sweep plane box intersection
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
#include "swp.h"
#include "alg.h"
#include "err.h"

typedef struct point POINT;
typedef struct sweep SWEEP;

struct point
{
  BOX *box;
  char high;  /* equal 0 or 3 */
};

struct sweep
{
  int boxnum;
  POINT *points;

  DRALG algo;
  void *dynrect;

  int time;
  int dims [3];
  double avglen [3];
};

#define BOX(pnt) ((pnt)->box)
#define ISLOW(pnt) ((pnt)->high == 0)
#define ISHIGH(pnt) ((pnt)->high == 3)

static int pntcmp0 (const POINT *a, const POINT *b)
{
  if (BOX (a)->extents [(int) a->high] < BOX (b)->extents [(int) b->high]) return -1;
  else if (BOX (a)->extents [(int) a->high] > BOX (b)->extents [(int) b->high]) return 1;
  else if (ISHIGH (a) && ISLOW (b)) return -1;
  else if (ISHIGH (b) && ISLOW (a)) return 1;
  else return 0;
}

static int pntcmp1 (const POINT *a, const POINT *b)
{
  if (BOX (a)->extents [(int) a->high + 1] < BOX (b)->extents [(int) b->high + 1]) return -1;
  else if (BOX (a)->extents [(int) a->high + 1] > BOX (b)->extents [(int) b->high + 1]) return 1;
  else if (ISHIGH (a) && ISLOW (b)) return -1;
  else if (ISHIGH (b) && ISLOW (a)) return 1;
  else return 0;
}

static int pntcmp2 (const POINT *a, const POINT *b)
{
  if (BOX (a)->extents [(int) a->high + 2] < BOX (b)->extents [(int) b->high + 2]) return -1;
  else if (BOX (a)->extents [(int) a->high + 2] > BOX (b)->extents [(int) b->high + 2]) return 1;
  else if (ISHIGH (a) && ISLOW (b)) return -1;
  else if (ISHIGH (b) && ISLOW (a)) return 1;
  else return 0;
}

typedef int (*qcmp) (const void*, const void*);

static int (*pntcmp [3]) (const void*, const void*) = { (qcmp) pntcmp0, (qcmp) pntcmp1, (qcmp) pntcmp2 };

static void getcenters (BOX **begin, BOX **end, double centers [3])
{
  BOX **box;

  centers [0] = centers [1] = centers [2] = 0.;

  for (box = begin; box != end; box ++)
  {
    centers [0] += (*box)->extents [0] + (*box)->extents [3];
    centers [1] += (*box)->extents [1] + (*box)->extents [4];
    centers [2] += (*box)->extents [2] + (*box)->extents [5];
  }

  centers [0] /= 2. * (double)(end - begin);
  centers [1] /= 2. * (double)(end - begin);
  centers [2] /= 2. * (double)(end - begin);
}

static void getavglens (BOX **begin, BOX **end, double avglen [3])
{
  BOX **box;

  avglen [0] = avglen [1] = avglen [2] = 0.;

  for (box = begin; box != end; box ++)
  {
    avglen [0] += (*box)->extents [3] - (*box)->extents [0];
    avglen [1] += (*box)->extents [4] - (*box)->extents [1];
    avglen [2] += (*box)->extents [5] - (*box)->extents [2];
  }

  avglen [0] /= (double)(end - begin);
  avglen [1] /= (double)(end - begin);
  avglen [2] /= (double)(end - begin);
}

static void getdivers (BOX **begin, BOX **end, double centers [3], double divers [3])
{
  BOX **box;
  double val [6];

  divers [0] = divers [1] = divers [2] = 0.;

  for (box = begin; box != end; box ++)
  {
    val [0] = (*box)->extents [0] - centers [0];
    val [3] = (*box)->extents [3] - centers [0];

    val [1] = (*box)->extents [1] - centers [1];
    val [4] = (*box)->extents [4] - centers [1];
    
    val [2] = (*box)->extents [2] - centers [2];
    val [5] = (*box)->extents [5] - centers [2];

    divers [0] += val [0] * val [0] + val [3] * val [3];
    divers [1] += val [1] * val [1] + val [4] * val [4];
    divers [2] += val [2] * val [2] + val [5] * val [5];
  }

  divers [0] /= 2. * (double)(end - begin);
  divers [1] /= 2. * (double)(end - begin);
  divers [2] /= 2. * (double)(end - begin);
}

static void sortdims (double divers [3], int dims [3])
{
  int tmp, i, j;

  dims [0] = 0;
  dims [1] = 1;
  dims [2] = 2;

  for (i = 0; i < 3; i ++)
    for (j = i + 1; j < 3; j ++)
      if (divers [dims [i]] < divers [dims [j]])
      {
        tmp = dims [i];
        dims [i] = dims [j];
        dims [j] = tmp;
      }
}

static void setpoints (BOX **begin, BOX **end, POINT *points)
{
  BOX **box;
  POINT *lo, *hi;

  for (box = begin, lo = points, hi = lo + 1; box != end; box ++, lo += 2, hi += 2)
  {
    lo->box = hi->box = *box;
    lo->high = 0;
    hi->high = 3;
  }
}

static void isort (POINT *begin, POINT *end, int d)
{
  POINT *p, *q, tmp;

  p = q = --end;

  do
  {
    while (q < end && pntcmp [d] (q, q+1) > 0)
    {
      tmp = *(q + 1);
      *(q + 1) = *q;
      *q = tmp;

      q ++;
    }

    q = --p;
    
  } while (p >= begin);
}

inline static void enter_sweep (void *dynrect, BOX *box, int *dims,
                                void *data, BOX_Overlap_Create report)
{
  DR_Insert (dynrect, box->extents [dims [1]], box->extents [dims [2]],
             box->extents [dims [1] + 3], box->extents [dims [2] + 3],
	     box, data, report);
}

inline static void leave_sweep (void *dynrect, BOX *box, int *dims)
{
  DR_Delete (dynrect, box->extents [dims [1]], box->extents [dims [2]],
             box->extents [dims [1] + 3], box->extents [dims [2] + 3],
	     box);
}

static void reinit (SWEEP *s, int boxnum)
{
  free (s->points);

  ERRMEM (s->points = malloc (sizeof (POINT) * boxnum * 2));

  s->boxnum = boxnum;

  DR_Destroy (s->dynrect);

  s->dynrect = DR_Create (boxnum, s->algo);
}

static void reinitdr (SWEEP *s, int boxnum, DRALG algo)
{
  DR_Destroy (s->dynrect);

  s->dynrect = DR_Create (boxnum, algo);

  s->algo = algo;
}

void* SWEEP_Create (int boxnum, DRALG algo)
{
  SWEEP *s;

  ERRMEM (s = malloc (sizeof (SWEEP)));

  ERRMEM (s->points = malloc (sizeof (POINT) * boxnum * 2));

  s->boxnum = boxnum;

  s->dynrect = DR_Create (boxnum, algo);
  s->algo = algo;
  s->time = 0;

  return  s;
}

void SWEEP_Do (void *context, DRALG algo, int changed, int boxnum, BOX **boxes, void *data, BOX_Overlap_Create report)
{
  SWEEP *s = context;
  BOX **end = boxes + boxnum;
  double centers [3], divers [3], *avglen = s->avglen;
  POINT *endp, *p;
  int *dims = s->dims;
  void *dynrect;

  if (changed) s->time = 0; /* force full reinitialisation (qsort'ed) */

  if (s->boxnum < boxnum) reinit (s, boxnum);
  else if (s->algo != algo) reinitdr (s, boxnum, algo);

  endp = s->points + boxnum * 2;
  dynrect = s->dynrect;

  if (s->time == 0)
  {
    getavglens (boxes, end, avglen);
    getcenters (boxes, end, centers);
    getdivers (boxes, end, centers, divers);
    sortdims (divers, dims);
    setpoints (boxes, end, s->points);
    qsort (s->points, endp - s->points, sizeof (POINT), pntcmp [dims [0]]);
    DR_Params (dynrect, avglen);
  }
  else
  {
    isort (s->points, endp, dims [0]); /* Insertion sort for consecutive runs */
  }

  for (p = s->points; p != endp; p ++)
  {
    if (p->high == 0)
    {
      enter_sweep (dynrect, p->box, dims, data, report);
    }
    else
    {
      leave_sweep (dynrect, p->box, dims);
    }
  }

  s->time = 1; /* no point to increment it */
}

void SWEEP_Destroy (void *context)
{
  SWEEP *s = context;

  free (s->points);
  DR_Destroy (s->dynrect);
  free (s);
}
