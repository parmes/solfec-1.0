/*
 * hsh.c
 * Copyright (C) 2005, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * ------------------------------------------------------------------------------
 * overlap detection  based on hashing box volumes and checking resulting list
 * of colliding hash table entries for overlaps.
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
#include <float.h>

#include "alg.h"
#include "mem.h"
#include "hsh.h"
#include "lis.h"
#include "err.h"

typedef struct link LINK;
typedef struct mybox MYBOX;
typedef struct hash HASH;

struct link
{
  double *extents;
  MYBOX *box;
  LINK *next, *snext;
};

struct mybox
{
  BOX *box;
  int num, mark;
};

struct hash
{
  int hsize;
  LINK **htable;

  MEM linkpool;
  
  int boxnum;
  MYBOX *boxes;
};

#define NKEYS 1024

#define LE0(q, p) (q)->extents [0] <= (p)->extents [0]
#define LE1(q, p) (q)->extents [1] <= (p)->extents [1]
#define LE2(q, p) (q)->extents [2] <= (p)->extents [2]

IMPLEMENT_LIST_SORT (SINGLE_LINKED, sort0, LINK, prev, next, LE0)
IMPLEMENT_LIST_SORT (SINGLE_LINKED, sort1, LINK, prev, next, LE1)
IMPLEMENT_LIST_SORT (SINGLE_LINKED, sort2, LINK, prev, next, LE2)

/* sorting along three dimensions */
static LINK* (*sort [3]) (LINK*) = {sort0, sort1, sort2};

static int keyscmp (const int *key1, const int *key2)
{
  return *key1 < *key2 ? -1 : 1;
}

typedef int (*qcmp) (const void*, const void*);

#define OUT(ext1, ext2, d) (GT (ext1 [d], ext2 [3 + d]) || GT (ext2 [d], ext1 [3 + d]))
#define MYBOX(lnk) ((lnk)->box)
#define BOXO(lnk) (MYBOX (lnk)->box)

static void onewayscan (LINK *list, int d,
                        void *data, BOX_Overlap_Create report)
{
  LINK *jp;
  int n;

  for (jp = list->next; jp && LE (jp->extents [d], list->extents [3 + d]); jp = jp->next)
  {
    if (MYBOX (jp)->mark >= MYBOX (list)->num) continue; /* hash list links are created in order of "my" boxes vector and then scanned
							    in the order of their consicutive creation; hence lower number boxes find
							    intersections with higher number ones; the reverse case can be skipped */

    for (n = 0; n < 3; n ++)
      if (d != n)  
      {
        if (OUT (list->extents, jp->extents, n)) break;
      }

    if (n == 3)
    {
      report (data, BOXO (list), BOXO (jp));
      MYBOX (jp)->mark = MYBOX (list)->num;
    }
  }
}

static LINK* linkalloc (HASH *h)
{
  LINK *lnk;

  ERRMEM (lnk = MEM_Alloc (&h->linkpool));

  return lnk;
}

#define linksclear(h) MEM_Release (&(h)->linkpool)

static void reinit (HASH *h, int boxnum)
{
  free (h->htable);

  ERRMEM (h->htable = malloc (sizeof (LINK*) * boxnum));

  h->hsize = boxnum; 

  free (h->boxes);
  
  ERRMEM (h->boxes = malloc (sizeof (MYBOX) * boxnum));

  h->boxnum = boxnum; 

  MEM_Release (&h->linkpool);
  MEM_Init (&h->linkpool, sizeof (LINK), boxnum);
}

void*  HASH_Create (int boxnum)
{
  HASH *h;

  ERRMEM (h = malloc (sizeof (HASH)));

  ERRMEM (h->htable = malloc (sizeof (LINK*) * boxnum));

  h->hsize = boxnum; 
  
  ERRMEM (h->boxes = malloc (sizeof (MYBOX) * boxnum));

  h->boxnum = boxnum; 

  MEM_Init (&h->linkpool, sizeof (LINK), boxnum);

  return h;
}

void HASH_Do (void *context, int boxnum, BOX **boxes, void *data, BOX_Overlap_Create report)
{
  HASH *h = context;
  double avsize = 0., a [2][3] = { {DBL_MAX,DBL_MAX,DBL_MAX} , {-DBL_MAX,-DBL_MAX,-DBL_MAX} };
  int keys [NKEYS], i, j, k, n, m, p, q [3][2], dmax, hsize;
  LINK **htable, *tail = NULL, *list = NULL, *lnk;
  BOX **ip, **end = boxes + boxnum;
  MYBOX *jp;

  if (h->boxnum < boxnum) reinit (h, boxnum);

  htable = h->htable;
  hsize = h->hsize;

  /* compute average extents and scene ranges */
  for (ip = boxes, jp = h->boxes; ip != end; ip ++, jp ++)
  {
    avsize += (*ip)->extents [3] - (*ip)->extents [0];
    avsize += (*ip)->extents [4] - (*ip)->extents [1];
    avsize += (*ip)->extents [5] - (*ip)->extents [2];

    a [0][0] = MIN (a [0][0],(*ip)->extents [0]);
    a [1][0] = MAX (a [1][0],(*ip)->extents [3]);

    a [0][1] = MIN (a [0][1],(*ip)->extents [1]);
    a [1][1] = MAX (a [1][1],(*ip)->extents [4]);

    a [0][2] = MIN (a [0][2],(*ip)->extents [2]);
    a [1][2] = MAX (a [1][2],(*ip)->extents [5]);

    jp->box = *ip;
    jp->num = (ip - boxes);
    jp->mark = 0;
  }

  a [1][0] -= a [0][0];
  a [1][1] -= a [0][1];
  a [1][2] -= a [0][2];

  /* maximally elongated dimension */
  if (a [1][0] > a [1][1]) dmax = 0;
  else dmax = 1;
  if (a [1][2] > a [1][dmax]) dmax = 2;
    
  avsize /= (double) (3 * boxnum);

  for (i = 0; i < hsize; i++) htable [i] = NULL;

  /* compute integer voxel ranges */
  for (ip = boxes, jp = h->boxes; ip != end; ip ++, jp ++)
  {
    q [0][0] = INTEGER ((*ip)->extents [0], avsize);
    q [0][1] = INTEGER ((*ip)->extents [3], avsize);
  
    q [1][0] = INTEGER ((*ip)->extents [1], avsize);
    q [1][1] = INTEGER ((*ip)->extents [4], avsize);

    q [2][0] = INTEGER ((*ip)->extents [2], avsize);
    q [2][1] = INTEGER ((*ip)->extents [5], avsize);

    m = 0;
    
    for (i = q [0][0]; i <= q [0][1]; i ++)
    for (j = q [1][0]; j <= q [1][1]; j ++)
    for (k = q [2][0]; k <= q [2][1]; k ++)
    {
      keys [m ++] = HASH3 (i, j, k, hsize);

loop:
      if (m == NKEYS) goto out;
    }

out:
    qsort (keys, m, sizeof (int), (qcmp) keyscmp); /* keys can repeat */

    for (n = 0; n < m; n ++)
    {
      p = keys [n];

      if (n && p == keys [n-1]) continue; /* skip repeated keys */

      lnk = linkalloc (h);

      lnk->extents = (*ip)->extents;
      lnk->box = jp;
      lnk->next = htable [p]; /* hash entry list */
      htable [p] = lnk;

      if (tail) tail->snext = lnk; /* list of successively created links */
      else list = lnk;
      tail = lnk;
    }
    
    if (m == NKEYS)
    {
      m = 0;
      goto loop;
    }
  }

  for (i = 0; i < hsize; i ++)
    if (htable [i])
    {
      htable [i] = sort [dmax] (htable [i]); /* short hash entry lists */
    }

  for (; list; list = list->snext) /* in the order of successively created links */
    onewayscan (list, dmax, data, report); /* one way scan from their position onward (reduces work) */

  linksclear (h);
}

void HASH_Destroy (void *context)
{
  HASH *h = context;

  free (h->htable);
  free (h->boxes);
  MEM_Release (&h->linkpool);
  free (h);
}
