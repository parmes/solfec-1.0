/*
 * cvx.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * collection of convex polyhedrons
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

#include <limits.h>
#include <float.h>
#include <string.h>
#include "mem.h"
#include "sol.h"
#include "hyb.h"
#include "gjk.h"
#include "spx.h"
#include "err.h"
#include "alg.h"
#include "msh.h"
#include "cvx.h"
#include "hul.h"
#include "pck.h"
#include "fem.h"
#include "kdt.h"

#define MEMINC 64 /* memory increment size used in convex manipulations */

/* point definition used
 * for convex splitting */
struct symbolic_point
{
  /* plane */
  double *plane;

  /* segment beginning */
  double *begin;

  /* segment end */
  double *end;

  /* resultant point
   * coordinates */
  double coord [3];
};

/* type for 'qsort' compare routine casting */
typedef int (*qcmp) (const void*, const void*);

/* point behind all planes query */
static int point_inside (int npla, double *pla, double *point)
{
  for (; npla > 0; pla += 4, npla --)
    if (PLANE (pla, point) > GEOMETRIC_EPSILON) return 0;

  return 1;
}

/* compute planes - CCW oriented faces */
inline static void computeplanes (CONVEX *cvx)
{
  double max, *a, *b, *c, *nl;
  int n, m;

  for (n = m = 0; n < cvx->nfac; n ++, m += (cvx->fac [m] + 1))
  {
    a = &cvx->cur [cvx->fac [m + 1]];
    b = &cvx->cur [cvx->fac [m + 2]];
    c = &cvx->cur [cvx->fac [m + 3]];
    nl =  &cvx->pla [n * 4];
    NORMAL (a, b, c, nl);
    MAXABS (nl, max);
    ASSERT (max > 0, ERR_CVX_ZERO_NORMAL);
    NORMALIZE (nl);
    nl [3] = - DOT (a, nl);
  }
}

/* update planes (normals unchanged) */
inline static void updateplanes (CONVEX *cvx)
{
  double *a, *nl;
  int n, m;

  for (n = m = 0; n < cvx->nfac; n ++, m += (cvx->fac [m] + 1))
  {
    a = &cvx->cur [cvx->fac [m + 1]];
    nl =  &cvx->pla [n * 4];
    nl [3] = - DOT (a, nl);
  }
}

/* compute convex volume */
static double volume (CONVEX *cvx)
{
  double *zero, *a, *b, *c;
  double J, volume;
  int i, j, *f;
 
  zero = cvx->ref;
  volume = 0.; 
  for (i = 0, f = cvx->fac; i < cvx->nfac; i ++, f = &f [f [0] + 1])
  {
    a = &cvx->ref [f [1]];
    for (j = 2; j < f [0]; j ++)
    {
      b = &cvx->ref [f [j]];
      c = &cvx->ref [f [j + 1]];

      J = simplex_J (zero, a, b, c);
      volume += simplex_1 (J, zero, a, b, c);
    }
  }

  return volume;
}

typedef struct edge EDGE;

struct edge
{
  double *v, *w; /* v < w */
};

/* create new edge (appropriately sort vertices) */
static EDGE* edge_new (MEM *mem, double *v, double *w)
{
  EDGE *e;

  ERRMEM (e = MEM_Alloc (mem));

  if (v < w)
  {
    e->v = v;
    e->w = w;
  }
  else
  {
    e->v = w;
    e->w = v;
  }

  return e;
}

static int edge_compare (EDGE *a, EDGE *b)
{
  if (a->v < b->v) return -1;
  else if (a->v == b->v && a->w < b->w) return -1;
  else if (a->v == b->v && a->w == b->w) return 0;
  else return 1;
}

/* size of a single convex */
static int size (CONVEX *cvx)
{
  int n, m, l, *fac = cvx->fac;

  /* memory occupied by face definitions */
  for (n = m = l = 0; n < cvx->nfac; n ++, m += (fac [m] + 1))
    l += (fac [m] + 1);

  return sizeof (CONVEX) + sizeof (double [6]) * cvx->nver +
    sizeof (double [4]) * cvx->nfac + sizeof (int) * (l + cvx->nfac);
}

/* copy single convex */
static CONVEX* copycvx (CONVEX *cvx)
{
  CONVEX *twin;

  ERRMEM (twin = malloc (size (cvx)));
  memcpy (twin, cvx, size (cvx));
  /* FEM specific => does not get inherited */
  twin->epn = NULL;
  twin->ele = NULL;
  twin->nele = 0;
  /* -------------------------------------- */
  twin->ref = (double*)(twin + 1);
  twin->cur = twin->ref + twin->nver * 3;
  twin->pla = twin->cur + twin->nver * 3;
  twin->surface = (int*) (twin->pla + twin->nfac * 4);
  twin->fac = twin->surface + twin->nfac;
  if (cvx->nadj)
  {
    ERRMEM (twin->adj = malloc (cvx->nadj * sizeof (CONVEX*)));
    for (int i = 0; i < cvx->nadj; i ++) twin->adj [i] = cvx->adj [i]; /* to be re-mapped outside */
    twin->nadj = cvx->nadj;
  }
  else
  {
    twin->adj = NULL;
    twin->nadj = 0;
  }
  twin->next = NULL;

  return twin;
}

/* cut single convex with a plane */
TRI* cut (CONVEX *cvx, double *point, double *normal, int *m)
{
  double val [3], pla [4], dpl [4], dp0, *ver, *hv, *vv, *ov, *pl;
  TRI *t, *q, *s, *hul, *out;
  MEM edgemem, setmem;
  SET *edges, *item;
  int k, n, *f, *g;
  EDGE *e;

  /* initialize pools */
  MEM_Init (&edgemem, sizeof (EDGE), 64);
  MEM_Init (&setmem, sizeof (SET), 64);
  out = hul = NULL;
  ver = cvx->cur;
  edges = NULL;
  hv = NULL;

  /* 4-plane */
  COPY (normal, pla);
  NORMALIZE (pla);
  pla [3] = - DOT (pla, point);

  /* create edges */
  for (f = cvx->fac, pl = cvx->pla, n = 0;
       n < cvx->nfac; f += f[0]+1, pl += 4, n ++)
  {
    ADD4 (pla, pl, dpl);
    dp0 = LEN4 (dpl);
    if (dp0 < GEOMETRIC_EPSILON) /* coplanar oposit face => skip this cut */
    {
      *m = 0;
      goto out; /* useful for pairs of internally adjacent convices; negligible on the boundary */
    }

    for (g = f+1; g < f+f[0]; g ++) /* create edges from consecutive vertex pairs */
    {
      e = edge_new (&edgemem, &ver [*g], &ver [*(g+1)]);
      if (!SET_Insert (&setmem, &edges, e, (SET_Compare) edge_compare))
	MEM_Free (&edgemem, e);
    }

    /* now the edge created by the last and the first face vertex */
    e = edge_new (&edgemem, &ver [*g], &ver [*(f+1)]);
    if (!SET_Insert (&setmem, &edges, e, (SET_Compare) edge_compare))
      MEM_Free (&edgemem, e);
  }

  /* hull vertices */
  ERRMEM (hv = malloc ((1 + SET_Size (edges)) * sizeof (double [3])));
  vv = hv;
  n = 0;

  /* split edges */
  for (item = SET_First (edges); item; item = SET_Next (item))
  {
    e = item->data;
    val [0] = PLANE (pla, e->v);
    val [1] = PLANE (pla, e->w);

    if (val [0] * val [1] < - GEOMETRIC_EPSILON * GEOMETRIC_EPSILON)
    {
      PLANESEG (pla, e->v, e->w, val);
      COPY (val, vv);
      vv += 3;
      n ++;
    }
    else if (fabs (val [0]) <= GEOMETRIC_EPSILON)
    {
      COPY (e->v, vv);
      vv += 3;
      n ++;
    }
    else if (fabs (val [1]) <= GEOMETRIC_EPSILON)
    {
      COPY (e->w, vv);
      vv += 3;
      n ++;
    }
  }

  if (n == 0)
  {
    k = 0;
    goto out;
  }

  /* add third vertex */
  SUB (hv, pla, vv);

  /* compute hull */
  hul = hull (hv, n+1, m);

  if (!hul)
  {
    k = 0;
    goto out;
  }

  /* select output triangles */
  for (k = 0, t = hul, q = t + (*m); t != q; t ++)
  {
    if (t->ver [0] != vv &&
	t->ver [1] != vv &&
	t->ver [2] != vv)
    {
      k ++;
    }
  }
  ERRMEM (out = MEM_CALLOC (k * sizeof (TRI) + n * sizeof (double [3])));
  ov = (double*) (out + k);
  memcpy (ov, hv, n * sizeof (double [3]));
  for (t = hul, q = t + (*m), s = out; t != q; t ++)
  {
    if (t->ver [0] != vv &&
	t->ver [1] != vv &&
	t->ver [2] != vv)
    {
      s->ver [0] = &ov [t->ver [0] - hv];
      s->ver [1] = &ov [t->ver [1] - hv];
      s->ver [2] = &ov [t->ver [2] - hv];
      s->ptr = cvx;
      COPY (t->out, s->out);
      s ++;
    }
  }

out:
  /* clean */
  if (hv) free (hv);
  if (hul) free (hul);
  MEM_Release (&setmem);
  MEM_Release (&edgemem);

  *m = k;
  return out;
}

/* convex triangulation contact test */
static int convex_touches_triangulation (CONVEX *cvx, KDT *kd)
{
  double extents [6], b [12], p [3], q [3], deps;
  SET *leaves = NULL, *item;
  KDT *leaf;
  TRI *t;
  int i;

  CONVEX_Extents (NULL, cvx, extents);
  KDT_Pick_Extents (kd, extents, &leaves);
  deps = 10 * GEOMETRIC_EPSILON;

  for (item = SET_First (leaves); item; item = SET_Next (item))
  {
    leaf = item->data;
    for (i = 0; i < leaf->n; i ++)
    {
      t = leaf->data [i];
      COPY (t->ver [0], b);
      COPY (t->ver [1], b+3);
      COPY (t->ver [2], b+6);
      ADDMUL (b, deps, t->out, b+9);
      if (gjk (cvx->cur, cvx->nver, b, 4, p, q) < GEOMETRIC_EPSILON)
      {
	SET_Free (NULL, &leaves);
	return 1;
      }
    }
  }

  SET_Free (NULL, &leaves);
  return 0;
}

/* split convex in two halves, using plane (point 'pt, normal 'nl') */
static int split (CONVEX *cvx, double *pt, double *nl, int surfid[2], CONVEX **one, CONVEX **two, MEM *extmem, MAP **below, MAP **above)
{
  double val [3], pla [4], dpla [4], dl, *ver, *ov, *tv, *vv, *oo, *tt;
  MEM edgemem, setmem, mapmem;
  int n, *f, *g, nov, ntv;
  MAP *normals, *jtem;
  SET *edges, *item;
  CONVEX *out;
  EDGE *e;

  /* in case convex is too small do not split it */
  if (volume (cvx) < pow (GEOMETRIC_EPSILON, 3)) return 0;

  /* initialize pools */
  MEM_Init (&edgemem, sizeof (EDGE), 64);
  MEM_Init (&setmem, sizeof (SET), 64);
  MEM_Init (&mapmem, sizeof (MAP), 64);
  normals = NULL;
  edges = NULL;
  ver = cvx->cur;

  for (f = cvx->fac, n = 0; n < cvx->nfac; f += f[0]+1, n ++)
  {
    for (g = f+1; g < f+f[0]; g ++) /* create edges from consecutive vertex pairs */
    {
      e = edge_new (&edgemem, &ver [*g], &ver [*(g+1)]);
      if (!SET_Insert (&setmem, &edges, e, (SET_Compare) edge_compare))
	MEM_Free (&edgemem, e);
    }

    /* now the edge created by the last and the first face vertex */
    e = edge_new (&edgemem, &ver [*g], &ver [*(f+1)]);
    if (!SET_Insert (&setmem, &edges, e, (SET_Compare) edge_compare))
      MEM_Free (&edgemem, e);
  }

  n = cvx->nver + SET_Size (edges); /* a safe size of split convex vertices */
  ERRMEM (ov = malloc (2 * n * sizeof (double [3])));
  tv = ov + 3 * n;

  COPY (nl, pla);
  NORMALIZE (pla);
  pla [3] = - DOT (pla, pt);

  /* split vertices (include epsilon margin near plane) */
  for (vv = ver, oo = ov, tt = tv, n = nov = ntv = 0;
       n < cvx->nver; vv += 3, n ++)
  {
    val [0] = PLANE (pla, vv);

    if (val [0] <= GEOMETRIC_EPSILON)
    {
      COPY (vv, oo);
      oo += 3;
      nov ++;
    }

    if (val [0] >= -GEOMETRIC_EPSILON)
    {
      COPY (vv, tt);
      tt += 3;
      ntv ++;
    }
  }

  /* split edges (exclude epsilon margin near plane) */
  for (item = SET_First (edges); item; item = SET_Next (item))
  {
    e = item->data;
    val [0] = PLANE (pla, e->v);
    val [1] = PLANE (pla, e->w);

    if (val [0] * val [1] < - GEOMETRIC_EPSILON * GEOMETRIC_EPSILON)
    {
      PLANESEG (pla, e->v, e->w, val);
      COPY (val, oo);
      COPY (val, tt);
      oo += 3;
      tt += 3;
      nov ++;
      ntv ++;
    }
  }

  if (nov && !ntv) /* copy input convex into the first list */
  {
    out = copycvx (cvx);
    MAP_Insert (extmem, below, cvx, out, NULL);
    *one = CONVEX_Glue (out, *one);
  }
  else if (ntv && !nov) /* copy input convex into the second list */
  {
    out = copycvx (cvx);
    MAP_Insert (extmem, above, cvx, out, NULL);
    *two = CONVEX_Glue (out, *two);
  }
  else /* create output vertices as convex hulls of point sets */
  {
    if (nov >= 3)
    {
      out = CONVEX_Hull (*one, ov, nov, cvx->surface [0], cvx->volume);
      if (cvx->nadj)
      {
	ERRMEM (out->adj = malloc (cvx->nadj * sizeof (CONVEX*)));
	for (int i = 0; i < cvx->nadj; i ++) out->adj [i] = cvx->adj [i]; /* to be re-mapped outside */
	out->nadj = cvx->nadj;
      }
      MAP_Insert (extmem, below, cvx, out, NULL);
      if (out) *one = out;
      else nov = 0; /* failed to create hull */
    }
    if (ntv >= 3)
    {
      out = CONVEX_Hull (*two, tv, ntv, cvx->surface [0], cvx->volume);
      if (cvx->nadj)
      {
	ERRMEM (out->adj = malloc (cvx->nadj * sizeof (CONVEX*)));
	for (int i = 0; i < cvx->nadj; i ++) out->adj [i] = cvx->adj [i]; /* to be re-mapped outside */
	out->nadj = cvx->nadj;
      }
      MAP_Insert (extmem, above, cvx, out, NULL);
      if (out) *two = out;
      else ntv = 0; /* failed to create hull */
    }
  }

  /* map normals to surface identifiers */
  for (oo = cvx->pla, f = cvx->surface, n = 0; n < cvx->nfac; oo += 4, f ++, n ++)
    MAP_Insert (&mapmem, &normals, oo, (void*) (long) (*f), (MAP_Compare) POINTS_COMPARE);

  /* try to associate normals of split
   * convices with surface identifiers */
  if (nov >= 3)
  {
    for (oo = (*one)->pla, f = (*one)->surface, n = 0; n < (*one)->nfac; oo += 4, f ++, n ++)
    {
      SUB4 (oo, pla, dpla); /* below => same plane orientation */
      dl = LEN4 (dpla);
      if (dl < GEOMETRIC_EPSILON) /* cutting plane */
      {
	*f = surfid[0];
      }
      else if ((jtem = MAP_Find_Node (normals, oo, (MAP_Compare) POINTS_COMPARE)))
      {
	*f = (int) (long) jtem->data;
      }
      else *f = cvx->surface [0]; /* if a normal could not be found use the first surface identifier */
    }
  }
  if (ntv >= 3)
  {
    for (oo = (*two)->pla, f = (*two)->surface, n = 0; n < (*two)->nfac; oo += 4, f ++, n ++)
    {
      ADD4 (oo, pla, dpla); /* above => reversed plane orientation */
      dl = LEN4 (dpla);
      if (dl < GEOMETRIC_EPSILON) /* cutting plane */
      {
	*f = surfid[1];
      }
      else if ((jtem = MAP_Find_Node (normals, oo, (MAP_Compare) POINTS_COMPARE)))
      {
	*f = (int) (long) jtem->data;
      }
      else *f = cvx->surface [0]; /* if a normal could not be found use the first surface identifier */
    }
  }

  /* clean up */
  free (ov);
  MEM_Release (&mapmem);
  MEM_Release (&setmem);
  MEM_Release (&edgemem);

  return 1;
}

/* size of faces definition */
static int facsize (CONVEX *cvx)
{
  int n, m, l, *fac = cvx->fac;

  /* memory occupied by face definitions */
  for (n = m = l = 0; n < cvx->nfac; n ++, m += (fac [m] + 1))
    l += (fac [m] + 1);

  return l;
}

/* overlap callback for convex adjacency */
static void overlap (void *data, BOX *one, BOX *two)
{
  double p [3], q [3];
  CONVEX *cvx = (CONVEX*)one->sgp,
	 *cvy = (CONVEX*)two->sgp;

  if (gjk (cvx->cur, cvx->nver, cvy->cur, cvy->nver, p, q) < GEOMETRIC_EPSILON) /* if they touch */
  {
    ERRMEM (cvx->adj = realloc (cvx->adj, (++cvx->nadj) * sizeof (CONVEX*)));  /* extend adjacency */
    cvx->adj [cvx->nadj-1] = cvy;
    ERRMEM (cvy->adj = realloc (cvy->adj, (++cvy->nadj) * sizeof (CONVEX*))); 
    cvy->adj [cvy->nadj-1] = cvx;
  }
}

/* recursive neighbour marking */
static void mark_neighs (CONVEX *x, short flag)
{
  if (x->flag != flag)
  {
    x->flag = flag;

    for (int i = 0; i < x->nadj; i ++)
    {
      mark_neighs (x->adj [i], flag);
    }
  }
}

/* add new convex to the 'cvx' list and return the appended list */
CONVEX* CONVEX_Create (CONVEX *cvx, double *ver, int nver, int *fac, int nfac, int *surfaces, int volume)
{
  int n, m, l;
  CONVEX *cvy;

  /* memory occupied by face definitions */
  for (n = m = l = 0; n < nfac; n ++, m += (fac [m] + 1))
    l += (fac [m] + 1);

  /* sum up remaining memory */
  m = sizeof (CONVEX) + sizeof (double [6]) * nver +
    sizeof (double [4]) * nfac + sizeof (int) * (l + nfac);
  ERRMEM (cvy = malloc (m));

  cvy->ref = (double*)(cvy + 1);
  cvy->cur = cvy->ref + nver * 3;
  cvy->pla = cvy->cur + nver * 3;
  cvy->surface = (int*) (cvy->pla + nfac * 4); /* space for surface identifiers */
  cvy->fac = cvy->surface + nfac;
  cvy->nver = nver;
  cvy->nfac = nfac;
  cvy->epn = NULL;
  cvy->adj = NULL;
  cvy->ele = NULL;
  cvy->nadj = 0;
  cvy->nele = 0;
  cvy->volume = volume; /* volume identifier */
  cvy->mat = NULL;

  memcpy (cvy->ref, ver, sizeof (double [3]) * nver);
  memcpy (cvy->cur, ver, sizeof (double [3]) * nver);
  if (surfaces) memcpy (cvy->surface, surfaces, sizeof (int) * nfac);
  else memset (cvy->surface, 0, sizeof (int) * nfac);
  memcpy (cvy->fac, fac, sizeof (int) * l);

  /* remap vertex pointers */
  for (n = m = 0; n < nfac; n ++,m += (fac [m] + 1))
    for (l = 1; l <= fac [m]; l ++)
      cvy->fac [m + l] = fac [m + l] * 3;

  /* calculate planes */
  computeplanes (cvy);
 
  /* append list */
  cvy->next = cvx;

  return cvy;
}

/* append convices list 'cvx' with a convex hull of the input point set */
CONVEX* CONVEX_Hull (CONVEX *cvx, double *pnt, int npnt, int surface, int volume)
{
  int i, j, n, m, *fac, *f;
  double *ver, *u, *v;
  MAP *map, *item;
  int *surfaces;
  CONVEX *cvy;
  TRI *tri;
  MEM mem;

  /* create convex hull */
  if (!(tri = hull (pnt, npnt, &m))) return NULL;

  /* convex hull is triangulated and its vertices
   * only point to the points in 'pnt' storage;
   * we now map nodes used in the hull */

  MEM_Init (&mem, sizeof (MAP), MEMINC);

  for (map = NULL, n = i = 0; i < m; i ++)
  {
    for (j = 0; j < 3; j ++)
    {
      if (MAP_Insert (&mem, &map, tri [i].ver [j], (void*) (long) n, NULL)) n ++;
    }
  }

  ERRMEM (ver = malloc (n * sizeof (double [3])));
  ERRMEM (fac = malloc (m * sizeof (int [4])));
  ERRMEM (surfaces = malloc (m * sizeof (int)));

  for (item = MAP_First (map); item; item = MAP_Next (item))
  {
    i = (int) (long) item->data;
    u = item->key;
    v = ver + i * 3;
    COPY (u, v); /* copy hull vertices */
  }

  for (f = fac, i = 0; i < m; f += 4, i ++)
  {
    f [0] = 3;

    for (j = 0; j < 3; j ++) f [j+1] = (int) (long) MAP_Find (map, tri [i].ver [j], NULL); /* map hull vertices */

    surfaces [i] = surface;
  }

  /* TODO: merge adjacent coplanar faces */

  cvy = CONVEX_Create (cvx, ver, n, fac, m, surfaces, volume);

  free (tri);
  free (ver);
  free (fac);
  free (surfaces);
  MEM_Release (&mem);

  return cvy;
}

/* glue to convex lists */
CONVEX* CONVEX_Glue (CONVEX *cvx, CONVEX *cvy)
{
  CONVEX *cvz = cvx;

  for (; cvx->next; cvx = cvx->next);
  cvx->next = cvy;
  return cvz;
}

/* comput adjacency pointers */
void CONVEX_Compute_Adjacency (CONVEX *cvx)
{
  MEM mem;
  BOX **boxes;
  CONVEX *cvy;
  int num;
  
  for (cvy = cvx, num = 0; cvy; cvy = cvy->next)
  {
    cvy->nadj = 0;
    num ++;
  }

  if (num < 2) return;

  MEM_Init (&mem, sizeof (BOX), num);
  ERRMEM (boxes = malloc (sizeof (AABB*) * num));
  for (cvy = cvx, num = 0; cvy; cvy = cvy->next, num ++)
  {
    ERRMEM (boxes [num] = MEM_Alloc (&mem));
    CONVEX_Extents (NULL, cvy, boxes [num]->extents); /* set up extents */
    boxes [num]->sgp = (SGP*)cvy;
  }

  hybrid (boxes, num, NULL, overlap); /* detect boxoverlaps => set adjacency inside the callback */

  MEM_Release (&mem); /* done */
  free (boxes);
}

/* copy convex list */
CONVEX* CONVEX_Copy (CONVEX *cvx)
{
  CONVEX *twin, *tail;
  MAP *map = NULL; /* adjacency map */
  MEM mem;
  int i;

  MEM_Init (&mem, sizeof (MAP), MEMINC);

  for (tail = NULL; cvx; cvx = cvx->next)
  {
    twin = copycvx (cvx);
    MAP_Insert (&mem, &map, cvx, twin, NULL);
    twin->next = tail;
    tail = twin;
  }

  /* map adjacency */
  for (cvx = twin; cvx; cvx = cvx->next)
  {
    for (i = 0; i < cvx->nadj; i ++)
    {
      ASSERT_DEBUG_EXT (cvx->adj [i] = MAP_Find (map, cvx->adj [i], NULL), "Inconsistent adjacency mapping");
    }
  }

  MEM_Release (&mem);

  return twin;
}

/* scale convex list */
void CONVEX_Scale (CONVEX *cvx, double *vector)
{
  double *v, *w;
  int n;

  for (; cvx; cvx = cvx->next)
  {
    for (n = 0, v = cvx->ref, w = cvx->cur; n < cvx->nver; n ++, v += 3, w += 3)
    {
      v [0] *= vector [0];
      v [1] *= vector [1];
      v [2] *= vector [2];
      w [0] *= vector [0];
      w [1] *= vector [1];
      w [2] *= vector [2];
    }

    updateplanes (cvx);
  }
}

/* translate convex list */
void CONVEX_Translate (CONVEX *cvx, double *vector)
{
  double *v, *w;
  int n;

  for (; cvx; cvx = cvx->next)
  {
    for (n = 0, v = cvx->ref, w = cvx->cur; n < cvx->nver; n ++, v += 3, w += 3)
    {
      v [0] += vector [0];
      v [1] += vector [1];
      v [2] += vector [2];
      w [0] += vector [0];
      w [1] += vector [1];
      w [2] += vector [2];
    }

    updateplanes (cvx);
  }
}

/* rotate convex list */
void CONVEX_Rotate (CONVEX *cvx, double *point, double *vector, double angle)
{
  double R [9], omega [3], *v, *w;
  int n;

  angle *=  ALG_PI / 180.0;
  COPY (vector, omega); 
  NORMALIZE (omega); 
  SCALE (omega, angle);
  EXPMAP (omega, R);

  for (; cvx; cvx = cvx->next)
  {
    for (n = 0, v = cvx->ref, w = cvx->cur; n < cvx->nver; n ++, v += 3, w += 3)
    {
      SUB (v, point, omega);
      NVADDMUL (point, R, omega, v);
      SUB (w, point, omega);
      NVADDMUL (point, R, omega, w);
    }

    computeplanes (cvx);
  }
}

/* cut through convices with a plane; return triangulated cross-section; vertices in the triangles
 * point to the memory allocated after the triangles memory; adjacency is not maintained;
 * TRI->ptr stores a pointer to the geometrical object that has been cut by the triangle */
TRI* CONVEX_Cut (CONVEX *cvx, double *point, double *normal, int *m)
{
  TRI **tri, *out, *t, *e, *q;
  int *ntr, n, s, i, j, k;
  SET *vertices, *item;
  double *v, *w, *z;
  KDT *kdtree, *kd;
  MEM setmem;

  n = 0;
  s = 128;
  (*m) = 0;
  z = NULL;
  out = NULL;
  kdtree = NULL;
  vertices = NULL;
  MEM_Init (&setmem, sizeof (SET), 128);
  ERRMEM (tri = malloc (s * sizeof (TRI*)));
  ERRMEM (ntr = malloc (s * sizeof (int)));

  for (; cvx; cvx = cvx->next)
  {
    tri [n] = cut (cvx, point, normal, &ntr [n]); /* find intersection */

    if (tri [n]) /* found */
    {
      for (t = tri [n], e = t + ntr [n]; t != e; t ++)
      {
	for (i = 0; i < 3; i++)
	{
	  SET_Insert (&setmem, &vertices, t->ver [i], NULL); /* collect vertices */
	}
      }
      (*m) += ntr [n];
      n ++;
    }

    if (n == s)
    {
      s *= 2;
      ERRMEM (tri = realloc (tri, s * sizeof (TRI*)));
      ERRMEM (ntr = realloc (ntr, s * sizeof (int)));
    }
  }

  if (n == 0) goto out;

  /* create kd-tree */
  i = SET_Size (vertices);
  ERRMEM (z = malloc (i * sizeof (double [3])));
  for (item = SET_First (vertices), v = z; item; item = SET_Next (item), v += 3)
  {
    w = item->data;
    COPY (w, v);
  }
  kdtree = KDT_Create (i, z, GEOMETRIC_EPSILON);
  i = KDT_Size (kdtree);

  /* output memory */
  ERRMEM (out = malloc ((*m) * sizeof (TRI) + i * sizeof (double [3])));
  v = (double*) (out + (*m));

  /* copy vertices */
  for (kd = KDT_First (kdtree); kd; kd = KDT_Next (kd))
  {
    double *a = kd->p, *b = &v [3 * kd->n];
    COPY (a, b);
  }

  /* compile output */
  for (i = 0, q = out; i < n; i ++)
  {
    for (t = tri [i], e = t + ntr [i]; t != e; t ++, q ++)
    {
      q->ptr = t->ptr;
      COPY (t->out, q->out);

      for (k = 0; k < 3; k ++)
      {
	kd = KDT_Nearest (kdtree, t->ver [k], GEOMETRIC_EPSILON);
	ASSERT_DEBUG (kd, "Kd-tree nearest point query failed");
	j = 3 * kd->n;
	q->ver [k] = &v [j]; /* map to new storage */
      }
    }
  }
out:
  for (i = 0; i < n; i ++) free (tri [i]);
  MEM_Release (&setmem);
  KDT_Destroy (kdtree);
  free (tri);
  free (ntr);
  free (z);

  return out;
}

/* split convices in two lists with plane defined by (point, normal); adjacencies between
 * the split lists elements need to be recomputed; surfid corresponds to the new surface;
 * topoadj != 0 implies cutting from the point and through the topological adjacency only */
void CONVEX_Split (CONVEX *cvx, double *point, double *normal, short topoadj, int surfid[2], CONVEX **one, CONVEX **two)
{
  MAP *below, *above, *other;
  CONVEX *x, *y, *three;
  double p [3], q [3];
  int mc, mb, i, j;
  TRI *c, *b;
  MEM mem;
  KDT *kd;

  below = above = other = NULL;
  *one = *two = three = NULL;
  kd = NULL;
  c = NULL;

  MEM_Init (&mem, sizeof (MAP), MEMINC);

  if (topoadj)
  {
    c = CONVEX_Cut (cvx, point, normal, &mc); /* compute cut surface triangulation */
    TRI_Compadj (c, mc); /* compute adjacency structure of triangles */
    b = TRI_Topoadj (c, mc, point, &mb); /* part of the triangulation adjacent to the point */
    ASSERT_DEBUG (b, "Input point is too far from mesh section by the splitting plane");
    if (mb < mc) /* two or more separated intersection surfaces */
      kd = TRI_Kdtree (b, mb);
  }

  for (x = cvx; x; x = x->next)
  {
    if (kd)
    {
      if (convex_touches_triangulation (x, kd))
	split (x, point, normal, surfid, one, two, &mem, &below, &above); /* split only if adjacent to the selected intersection patch */
      else
      {
	y = copycvx (x);
	MAP_Insert (&mem, &other, x, y, NULL);
	three = CONVEX_Glue (y, three);
      }
    }
    else split (x, point, normal, surfid, one, two, &mem, &below, &above); /* split all */
  }

  /* map adjacency */
  CONVEX *vcvx [] = {*one, *two};
  MAP *vmap [] = {below, above, other};

  for (j = 0; j < 2; j ++)
  {
    for (y = vcvx [j]; y; y = y->next)
    {
      for (i = 0; i < y->nadj; i ++)
      {
	x = MAP_Find (vmap [j], y->adj [i], NULL);

	if (x) /* neighbour was mapped on the same side of the splitting plane */
	{
	  y->adj [i] = x; /* re-map old pointers */
	}
	else
	{
	  x = MAP_Find (other, y->adj [i], NULL);

	  if (x && gjk (x->cur, x->nver, y->cur, y->nver, p, q) < GEOMETRIC_EPSILON)
	  {
	    y->adj [i] = x; /* re-map old pointers */
	  }
	  else /* broken adjacency */
	  {
	    y->adj [i] = y->adj [-- y->nadj]; /* remove from adjacency */
	    i --; /* start again */
	  }
	}
      }
    }
  }

  if (kd)
  {
    for (y = three; y; y = y->next) /* map remaining adjacency */
    {
      for (i = 0; i < y->nadj; i ++)
      {
	for (j = 0; j < 3; j ++)
	{
	  x = MAP_Find (vmap [j], y->adj [i], NULL);

	  if (x && gjk (x->cur, x->nver, y->cur, y->nver, p, q) < GEOMETRIC_EPSILON)
	  {
	    y->adj [i] = x;
	    break;
	  }
	}

	ASSERT_DEBUG (j < 3, "Inconsistent adjacency mapping");
      }
    }

    *one = CONVEX_Glue (*one, *two); /* no fragmentation => the result goes into 'one' */
    *one = CONVEX_Glue (*one, three);
    *two = NULL;

    KDT_Destroy (kd);
    free (c);
  }
  
  MEM_Release (&mem);
}

/* is convex set separable into disjoint parts */
int CONVEX_Separable (CONVEX *cvx)
{
  CONVEX *x;

  for (x = cvx; x; x = x->next) x->flag = 0;

  int flag = 1;

  mark_neighs (cvx, flag);

  for (x = cvx; x; x = x->next)
  {
    if (x->flag == 0)
    {
      flag ++;
      mark_neighs (x, flag);
    }
  }

  return (flag == 1 ? 0 : flag);
}

/* separate convex set into disjoint parts */
CONVEX** CONVEX_Separate (CONVEX *cvx, int *m)
{
  CONVEX *x, *y, **out;
  MAP *map = NULL; /* adjacency map */
  MEM mem;
  int i, j;

  MEM_Init (&mem, sizeof (MAP), MEMINC);
 
  *m = CONVEX_Separable (cvx);

  if ((*m) == 0) return NULL;

  ERRMEM (out = MEM_CALLOC ((*m) * sizeof (MESH*)));

  for (i = 1; i <= (*m); i ++)
  {
    for (x = cvx; x; x = x->next)
    {
      if (x->flag == i)
      {
	y = copycvx (x);
	MAP_Insert (&mem, &map, x, y, NULL);
	out [i-1] = CONVEX_Glue (y, out [i-1]); /* copy this convex into the output list */
      }
    }
  }

  /* map adjacency */
  for (j = 0; j < *m; j ++)
  {
    for (y = out [j]; y; y = y->next)
    {
      for (i = 0; i < y->nadj; i ++)
      {
	ASSERT_DEBUG_EXT (y->adj [i] = MAP_Find (map, y->adj [i], NULL), "Inconsistent adjacency mapping");
      }
    }
  }

  MEM_Release (&mem);

  return out;
}

/* compute partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void CONVEX_Char_Partial (CONVEX *cvx, int ref, double *vo, double *sx, double *sy, double *sz, double *eul)
{
  double zero [3] = {0, 0, 0},
	 J, *a, *b, *c, *ver;
  int i, j, *f;

  for (; cvx; cvx = cvx->next)
  {
    ver = ref ? cvx->ref : cvx->cur;

    for (i = 0, f = cvx->fac; i < cvx->nfac; i ++, f = &f [f [0] + 1])
    {
      a = &ver [f [1]];
      for (j = 2; j < f [0]; j ++)
      {
	b = &ver [f [j]];
	c = &ver [f [j + 1]];

	J = simplex_J (zero, a, b, c);
	*vo += simplex_1 (J, zero, a, b, c);
	*sx += simplex_x (J, zero, a, b, c);
	*sy += simplex_y (J, zero, a, b, c);
	*sz += simplex_z (J, zero, a, b, c);
	eul [0] += simplex_xx (J, zero, a, b, c);
	eul [3] += simplex_xy (J, zero, a, b, c);
	eul [4] += simplex_yy (J, zero, a, b, c);
	eul [6] += simplex_xz (J, zero, a, b, c);
	eul [7] += simplex_yz (J, zero, a, b, c);
	eul [8] += simplex_zz (J, zero, a, b, c);
      }
    }
  }
}

/* compute volume characteristics */
void CONVEX_Char (CONVEX *cvx, int ref, double *volume, double *center, double *euler)
{
  double vo, sx, sy, sz,
	 cen [3], eul [9];

  vo = sx = sy = sz = 0.0;
  SET9 (eul, 0.0);

  CONVEX_Char_Partial (cvx, ref, &vo, &sx, &sy, &sz, eul);

  cen [0] = sx / vo;
  cen [1] = sy / vo;
  cen [2] = sz / vo;

  eul [0] -= (2*sx - cen[0]*vo)*cen[0];
  eul [4] -= (2*sy - cen[1]*vo)*cen[1];
  eul [8] -= (2*sz - cen[2]*vo)*cen[2];
  eul [3] -= cen[0]*sy + cen[1]*sx - cen[0]*cen[1]*vo;
  eul [6] -= cen[0]*sz + cen[2]*sx - cen[0]*cen[2]*vo;
  eul [7] -= cen[1]*sz + cen[2]*sy - cen[1]*cen[2]*vo;
  eul [1] = eul[3];
  eul [2] = eul[6];
  eul [5] = eul[7];

  if (volume) *volume = vo;
  if (center) COPY (cen, center);
  if (euler) NNCOPY (eul, euler);
}

/* compute convex colume (referential
 * if ref == 1, or current otherwise) */
double CONVEX_Volume (CONVEX *cvx, int ref)
{
  double zero [3] = {0, 0, 0}, volume = 0.0,
	 *ver = ref ? cvx->ref : cvx->cur,
	 J, *a, *b, *c;
  int i, j, *f;

  for (; cvx; cvx = cvx->next)
  {
    for (i = 0, f = cvx->fac; i < cvx->nfac; i ++, f = &f [f [0] + 1])
    {
      a = &ver [f [1]];
      for (j = 2; j < f [0]; j ++)
      {
	b = &ver [f [j]];
	c = &ver [f [j + 1]];

	J = simplex_J (zero, a, b, c);
	volume += simplex_1 (J, zero, a, b, c);
      }
    }
  }

  return volume;
}


/* update bounding boxe of an individual convex */
void CONVEX_Extents (void *data, CONVEX *cvx, double *extents)
{
  double *cur;
  int n;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;
    
  for (cur = cvx->cur, n = 0; n < cvx->nver; cur += 3, n ++)
  {
    if (cur [0] < extents [0]) extents [0] = cur [0];
    if (cur [1] < extents [1]) extents [1] = cur [1];
    if (cur [2] < extents [2]) extents [2] = cur [2];
    if (cur [0] > extents [3]) extents [3] = cur [0];
    if (cur [1] > extents [4]) extents [4] = cur [1];
    if (cur [2] > extents [5]) extents [5] = cur [2];
  }

  extents [0] -= GEOMETRIC_EPSILON;
  extents [1] -= GEOMETRIC_EPSILON;
  extents [2] -= GEOMETRIC_EPSILON;
  extents [3] += GEOMETRIC_EPSILON;
  extents [4] += GEOMETRIC_EPSILON;
  extents [5] += GEOMETRIC_EPSILON;
}

/* update oriented extents of an individual convex */
void CONVEX_Oriented_Extents (CONVEX *cvx, double *vx, double *vy, double *vz, double *extents)
{
  double *cur;
  double e [3];
  int n;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;
    
  for (cur = cvx->cur, n = 0; n < cvx->nver; cur += 3, n ++)
  {
    e [0] = DOT (vx, cur);
    e [1] = DOT (vy, cur);
    e [2] = DOT (vz, cur);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [0] > extents [3]) extents [3] = e [0];
    if (e [1] > extents [4]) extents [4] = e [1];
    if (e [2] > extents [5]) extents [5] = e [2];
  }
}

/* compute extents of convex list */
void CONVEX_List_Extents (CONVEX *cvx, double *extents)
{
  double e [6];

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;
    
  for (; cvx; cvx = cvx->next)
  {
    CONVEX_Extents (NULL, cvx, e);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [3] > extents [3]) extents [3] = e [3];
    if (e [4] > extents [4]) extents [4] = e [4];
    if (e [5] > extents [5]) extents [5] = e [5];
  }
}

/* compute oriented extents of convex list */
void CONVEX_List_Oriented_Extents (CONVEX *cvx, double *vx, double *vy, double *vz, double *extents)
{
  double e [6];

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;
    
  for (; cvx; cvx = cvx->next)
  {
    CONVEX_Oriented_Extents (cvx, vx, vy, vz, e);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [3] > extents [3]) extents [3] = e [3];
    if (e [4] > extents [4]) extents [4] = e [4];
    if (e [5] > extents [5]) extents [5] = e [5];
  }
}

/* return first not NULL bulk material for a convex list */
void* CONVEX_First_Bulk_Material (CONVEX *cvx)
{
  for (; cvx; cvx = cvx->next)
    if (cvx->mat) return cvx->mat;

  return NULL;
}

/* return convex containing the point */
CONVEX* CONVEX_Containing_Point (CONVEX *cvx, double *point)
{
  for (; cvx; cvx = cvx->next)
    if (point_inside (cvx->nfac, cvx->pla, point)) return cvx;

  return NULL;
}

/* does this convex (not a list) contain the point? */
int CONVEX_Contains_Point (void *dummy, CONVEX *cvx, double *point)
{
  return point_inside (cvx->nfac, cvx->pla, point);
}

/* return distance of a spatial point to the convex */
double CONVEX_Spatial_Point_Distance (void *dummy, CONVEX *cvx, double *point)
{
  double q [3];

  return gjk_convex_point (cvx->cur, cvx->nver, point, q);
}

/* update convex list according to the given motion */
void CONVEX_Update (CONVEX *cvx, void *body, void *shp, MOTION motion)
{
  double *ref, *cur;
  ELEPNT *epn;
  int n;

  for (; cvx; cvx = cvx->next)
  {	
    if (cvx->epn && motion == (MOTION) BODY_Cur_Point) /* use fast update for FEM bodies with rough mesh */
    {
      for (epn = cvx->epn, ref = cvx->ref, cur = cvx->cur, n = 0; n < cvx->nver; epn ++, ref += 3, cur += 3, n ++)
	FEM_Cur_Point_Ext (body, epn->ele, ref, epn->pnt, cur);
    }
    else if (cvx->epn && motion)
    {
      for (epn = cvx->epn, ref = cvx->ref, cur = cvx->cur, n = 0; n < cvx->nver; epn ++, ref += 3, cur += 3, n ++)
	motion (body, epn->ele, epn->pnt, cur); /* modal_motion in fem.c */
    }
    else /* use regular update for other body types */
    {
      for (ref = cvx->ref, cur = cvx->cur, n = 0; n < cvx->nver; ref += 3, cur += 3, n ++)
      {
	if (motion) 
	{
	  SGP sgp = {shp, cvx, GOBJ_CONVEX, NULL};
	  motion (body, &sgp, ref, cur); /* move current nodes */
	}
	else { COPY (ref, cur); } /* restore reference configuration */
      }
    }

    /* calculate planes */
    computeplanes (cvx);
  }
}

/* test wether two convices are adjacent */
int CONVEX_Adjacent (CONVEX *one, CONVEX *two)
{
  int n;

  /* note that convex adjacency involves
   * face, edge, and vertex adjacent cases */
  for (n = 0; n < two->nadj; n ++)
    if (two->adj [n] == one) return 1; /* enough to compare one way (adjacency lists are symmetric) */

  return 0;
}

/* return 6-vector (normal, point) planes of convex faces */
double* CONVEX_Planes (CONVEX *cvx)
{
  double *cur, *pla, *p, *q;
  int i, *f;

  ERRMEM (p = malloc (cvx->nfac * sizeof (double [6])));

  for (i = 0, cur = cvx->cur, pla = cvx->pla, q = p, f = cvx->fac;
       i < cvx->nfac; i ++, pla += 4, q += 6, f += f[0]+1)
  {
    COPY (pla, q);
    COPY (&cur [f[1]], q+3);
  }

  return p;
}

/* free list of convices */
void CONVEX_Destroy (CONVEX *cvx)
{
  CONVEX *nxt;

  for (;cvx; cvx = nxt)
  {
    nxt = cvx->next;
    free (cvx->adj);
    free (cvx->ele);
    free (cvx->epn);
    free (cvx);
  }
}

void CONVEX_Pack (CONVEX *cvx, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  CONVEX *ptr;
  int count, n;
  MAP *map;

  for (count = 0, ptr = cvx; ptr; ptr = ptr->next) count ++; /* number of convices */

  for (map = NULL, n = 0, ptr = cvx; ptr; ptr = ptr->next, n ++)
    MAP_Insert (NULL, &map, ptr, (void*) (long) n, NULL); /* map pointers to table indices */

  pack_int (isize, i, ints, count); /* number of convices */

  for (; cvx; cvx = cvx->next)
  {
    pack_int (isize, i, ints, cvx->nver);
    pack_int (isize, i, ints, cvx->nfac);
    pack_int (isize, i, ints, cvx->nadj);
    pack_int (isize, i, ints, cvx->volume);
    pack_int (isize, i, ints, (n = facsize (cvx)));

    pack_ints (isize, i, ints, cvx->fac, n);
    pack_ints (isize, i, ints, cvx->surface, cvx->nfac);

    pack_doubles (dsize, d, doubles, cvx->pla, cvx->nfac * 4);
    pack_doubles (dsize, d, doubles, cvx->cur, cvx->nver * 3);
    pack_doubles (dsize, d, doubles, cvx->ref, cvx->nver * 3);

    /* rather than adjacency pack indices of neighbours in the output sequence */
    for (n = 0; n < cvx->nadj; n ++) pack_int (isize, i, ints, (int) (long) MAP_Find (map, cvx->adj [n], NULL));

    pack_int (isize, i, ints, cvx->mat ? 1 : 0); /* pack material existence flag */
    if (cvx->mat) pack_string (isize, i, ints, cvx->mat->label);
  }

  MAP_Free (NULL, &map);
}

CONVEX* CONVEX_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int count, k, nver, nfac, nadj, volume, facsi, size, n, j;
  CONVEX *ptr, **tab, *tail = NULL, *head;
  
  count = unpack_int (ipos, i, ints); /* number of convices */

  ERRMEM (tab = malloc (count * sizeof (CONVEX*)));
  
  for (k = 0; k < count; k ++) /* unpack convices */
  {
    nver = unpack_int (ipos, i, ints);
    nfac = unpack_int (ipos, i, ints);
    nadj = unpack_int (ipos, i, ints);
    volume = unpack_int (ipos, i, ints);
    facsi = unpack_int (ipos, i, ints);

    size = nver * sizeof (double [6]) +
           nfac * (sizeof (int) + sizeof (double [4])) +
	   facsi * sizeof (int) +
	   sizeof (CONVEX); /* total size of contiguous storage */

    ERRMEM (ptr = MEM_CALLOC (size));
    ptr->ref = (double*)(ptr + 1);
    ptr->cur = ptr->ref + nver * 3;
    ptr->pla = ptr->cur + nver * 3;
    ptr->surface = (int*) (ptr->pla + nfac * 4);
    ptr->fac = ptr->surface + nfac;
    ERRMEM (ptr->adj = malloc (nadj * sizeof (CONVEX*)));
    ptr->nver = nver;
    ptr->nfac = nfac;
    ptr->nadj = 0;
    ptr->volume = volume;
    tab [k] = ptr;

    if (tail)
    {
      tail->next = ptr;
      tail = ptr;
    }
    else head = tail = ptr;

    unpack_ints (ipos, i, ints, ptr->fac, facsi);
    unpack_ints (ipos, i, ints, ptr->surface, nfac);

    unpack_doubles (dpos, d, doubles, ptr->pla, nfac * 4);
    unpack_doubles (dpos, d, doubles, ptr->cur, nver * 3);
    unpack_doubles (dpos, d, doubles, ptr->ref, nver * 3);

    for (n = 0; n < nadj; n ++)
    {
      j = unpack_int (ipos, i, ints);
      ptr->adj [ptr->nadj ++] = (CONVEX*) (long) j; /* store index in 'tab' for the moment */
    }

    j = unpack_int (ipos, i, ints); /* unpack material existence flag */

    if (j)
    {
      SOLFEC *sol = solfec;
      char *label = unpack_string (ipos, i, ints);
      ASSERT_DEBUG_EXT (ptr->mat = MATSET_Find (sol->mat, label), "Failed to find material when unpacking a convex");
      free (label);
    }
    else ptr->mat = NULL;
  }

  /* now map adjacency */
  for (k = 0; k < count; k ++)
  {
    ptr = tab [k];

    for (n = 0; n < ptr->nadj; n ++)
    {
      ASSERT_DEBUG (0 <= (int) (long) ptr->adj [n] && (int) (long) ptr->adj [n] < count, "Adjacent convex index out of bounds");
      ptr->adj [n] = tab [(int) (long) ptr->adj [n]];
    }
  }

  free (tab);

  return head;
}

/* export MBFCP definition */
void CONVEX_2_MBFCP (CONVEX *cvx, FILE *out)
{
  int i, n, *f, *s;
  CONVEX *ptr;
  double *d;

  for (ptr = cvx, n = 0; ptr; ptr = ptr->next, n ++);

  fprintf (out, "CONVEXES:\t%d\n", n);

  for (ptr = cvx; ptr; ptr = ptr->next)
  {
    fprintf (out, "VERTEXES:\t%d\n", ptr->nver);
    for (d = ptr->ref, n = 0; n < ptr->nver; d += 3, n ++)
    {
      fprintf (out, "%g  %g  %g\n", d [0], d [1], d [2]);
    }
    fprintf (out, "FACES:\t%d\n", ptr->nfac);
    for (f = ptr->fac, s = ptr->surface, n = 0; n < ptr->nfac; f += f[0]+1, n ++, s ++)
    {
      fprintf (out, "%d  ", f [0]);
      for (i = 0; i < f [0]; i ++) fprintf (out, "%d  ", f [i+1]);
      fprintf (out, "%d\n", s [0]);
    }
  }
}
