/*
 * tri.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * three-dimensional triangle
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
#include "tri.h"
#include "mem.h"
#include "map.h"
#include "set.h"
#include "err.h"
#include "alg.h"
#include "spx.h"

/* auxiliary pointer pair */
struct pp { void *a, *b; };

/* pointer pair comparison */
static int ppcmp (struct pp *x, struct pp *y)
{
  if (x->a < y->a) return -1;
  else if (x->a == y->a && x->b < y->b) return -1;
  else if (x->a == y->a && x->b == y->b) return 0;
  else return 1;
}

/* does triangle 'tri' have edge (v,w) ? */
inline static int hasedge (TRI *tri, double *v, double *w)
{
  double **ver = tri->ver;
  int i, j;

  for (i = j = 0; i < 3; i ++)
    if (ver [i] == v || ver [i] == w) j ++;

  return j == 2 ? 1 : 0;
}

/* return the next after 'tri' CCW triangle around 'v' */
inline static TRI* nextaround (TRI *tri, double *v)
{
  if (v == tri->ver [0]) return tri->adj [2];
  else if (v == tri->ver [1]) return tri->adj [0];
  else return tri->adj [1];
}

/* sort adjacency */
void TRI_Sortadj (TRI *tri)
{
  double **ver = tri->ver,
	 *edg [3][2] = {{ver[0], ver[1]},
	                {ver[1], ver[2]},
			{ver[2], ver[0]}};
  TRI *out [3] = {NULL, NULL, NULL},
      **adj = tri->adj;
  int i, j;

  for (i = 0; i < 3; i ++)
  {
    for (j = 0; j < 3; j ++)
    {
      if (adj [j] && hasedge (adj [j], edg [i][0], edg [i][1]))
      {
	out [i] = adj [j]; /* select adj[j] if it is adjacent to the current edge */
	adj [j] = NULL; /* mark as already selected */
      }
    }
  }

  adj [0] = out [0];
  adj [1] = out [1];
  adj [2] = out [2];
}

/* make 'p' and 'q' adjacent */
int TRI_Addadj (TRI *p, TRI *q)
{
  int i, j;

  for (i = 0; i < 3 && p->adj[i]; i ++);
  for (j = 0; j < 3 && q->adj[j]; j ++);
  if (i >= 3 || j >= 3) return 0;
  p->adj [i] = q;
  q->adj [j] = p;
  return 1;
}

/* compy into a compact memory block */
TRI* TRI_Copy (TRI *tri, int n)
{
  int vcnt; /* number of vertices */
  MAP *vm, *im; /* map of visited vertices; iterator 'im' */
  MEM mem; /* memory for map items */
  TRI *t, *s, *e, *o; /* triangle iterators 't' and 's', table end 'e' and output 'o' */
  double *v;
  int i;

  MEM_Init (&mem, sizeof (MAP), n * 3);
  e = tri + n;
  vm = NULL;
  vcnt = 0;

  for (t = tri; t < e; t ++) /* for each triangle */
  {
    for (i = 0; i < 3; i ++) /* for each vertex */
    {
      if (! MAP_Find_Node (vm, t->ver [i], NULL)) /* if the vertex was not mapped before */
      {
	vcnt ++; /* make account of it */
	MAP_Insert (&mem, &vm, t->ver [i], NULL, NULL); /* map vertex 'v' */
      }
    }
  }

  /* alloc output memory */
  ERRMEM (o = malloc (sizeof (TRI)*n + sizeof (double [3]) * vcnt));
  v = (double*)(o + n);

  /* map old vertices into new placeholders */
  for (i = 0, im = MAP_First (vm); im; im = MAP_Next (im), i ++) im->data = (v + i*3);

  for (t = tri, s = o; t < e; t ++, s ++) /* again => for each triangle */
  {
    COPY (t->out, s->out);
    s->flg = t->flg;
    
    for (i = 0; i < 3; i ++)
    {
      if (t->adj [i]) s->adj [i] = o + (t->adj[i] - tri); /* map adjacency */
      s->ver [i] = MAP_Find (vm, t->ver [i], NULL); /* map new vertex */
    }
  }

  MEM_Release (&mem);

  return o;
}

/* compute adjacency structure */
void TRI_Compadj (TRI *tri, int n)
{
  MAP *vm; /* map of vertex pairs */
  TRI *t, *s, *e;
  MEM mm, mp;
  int i;

  MEM_Init (&mm, sizeof (MAP), n * 3);
  MEM_Init (&mp, sizeof (struct pp), n * 3);
  e = tri + n;
  vm = NULL;

  /* zero adjacency */ 
  for (t = tri; t < e; t ++)
    t->adj [0] = t->adj [1] =
    t->adj [2] = NULL;

  for (t = tri; t < e; t ++) /* for each triangle */
  {
    for (i = 0; i < 3; i ++) /* for each vertex */
    {
      struct pp x = {t->ver [i], t->ver[i < 2 ? i+1 : 0]};
      if (x.a > x.b) { void *c = x.b; x.b = x.a; x.a = c; }

      if (!(s =  MAP_Find (vm, &x, (MAP_Compare)ppcmp))) /* if the edge was not mapped before */
      {
	struct pp *y;

	ERRMEM (y = MEM_Alloc (&mp)); *y = x;
	MAP_Insert (&mm, &vm, y, t, (MAP_Compare)ppcmp); /* map edge to 't' */
      }
      else TRI_Addadj (t, s); /* the edge was mapped before => adjacency found */
    }
  }

  MEM_Release (&mm);
  MEM_Release (&mp);
}

/* compute polar polyhedron of (tri, n) */
PFV* TRI_Polarise (TRI *tri, int n, int *m)
{
  int pfvcnt; /* number of vertices of all polar faces */
  PFV *pfv, *p, *q; /* first 'pfcnt' entries are polar face vertex list heads, the rest is list memory of size (pfvcnt - pfcnt); and iterator 'p' */
  MAP *vm, *im; /* map of visited vertices of (tri, n) (polar faces); iterator 'im' */
  MEM mem; /* memory for map items */
  TRI *t, *s, *e; /* triangle iterators 't' and 's', and table end 'e' */
  double *v, *w, d, x;
  int i, j;

  MEM_Init (&mem, sizeof (MAP), n * 3);
  e = tri + n;
  pfvcnt = 0;
  vm = NULL;
  pfv = NULL;

  for (t = tri; t < e; t ++) /* for each triangle */
  {
    for (i = 0; i < 3; i ++) /* for each vertex */
    {
      if (! MAP_Find (vm, t->ver [i], NULL)) /* if the vertex was not mapped before */
      {
	v = t->ver [i]; /* this is new vertex => new polar face */
	pfvcnt ++; /* 't' is the first vertex of the new face */
	
        for (s = nextaround (t, v), j = 0; s && s != t && j < n; s = nextaround (s, v), j ++) pfvcnt ++; /* add more vertices */

#if GEOMDEBUG	
	ASSERT_DEBUG (s && j < n, "Topological inconsistency in the input triangle mesh (a hole was found)"); /* a hole was found */
#else
	if (!s || j >= n) goto error;
#endif

	MAP_Insert (&mem, &vm, v, t, NULL); /* map vertex 'v' to triangle 't' */
      }
    }
  }

  /* alloc output memory => PFVs and 'n' vertices */
  ERRMEM (pfv = malloc (sizeof (PFV) * pfvcnt + sizeof (double [3]) * n));
  w = (double*) (pfv + pfvcnt);

  /* compute coordinates */
  for (t = tri; t < e; t ++)
  {
    v = w + 3*(t-tri);
    for (i = 0, d = 0; i < 3; i ++) { x = - DOT (t->out, t->ver [i]); if (x < d) d = x; } /* find smallest */

#if GEOMDEBUG
    ASSERT_DEBUG (d <= GEOMETRIC_EPSILON, /* near zero accepted but not too much on the positive side */
      "Geometric inconsistency during triangle mesh polarisation (half-planes must include zero)");
#else
    if (d > GEOMETRIC_EPSILON) goto error;
#endif

    if (d > -GEOMETRIC_EPSILON) d = -GEOMETRIC_EPSILON; /* ugly regularisation */
    DIV (t->out, -d, v); /* <n,x> <= 1 */
  }

  /* now go again and create polar face lists */ 
  for (i = 0, im = MAP_First (vm); im; im = MAP_Next (im), i ++)
  {
    v = im->key;
    t = im->data;
    pfv [i].coord = w + (t-tri)*3; /* map coord to the memory placed at the end of 'pfv' block */
    pfv [i].nl = v; /* common normal */
    q = &pfv [i]; /* list tail */

    for (s = t, t = nextaround (t, v);
	 t != s; t = nextaround (t, v))
    {
      p = &pfv [--pfvcnt]; /* get list item from the end => leave the beginning for the return table */
      p->coord = w + (t-tri)*3; /* map coord ... */
      p->nl = v; /* common normal */
      q->n = p; /* add to tail */
      q = p; /* update tail */
    }

    q->n = &pfv [i]; /* circular list */
  }

  goto done;

error:
  if (pfv) free (pfv);
  pfv = NULL;
  i = 0;

done:
  MEM_Release (&mem);

  (*m) = i;
  return pfv;
}

/* copute vertices */
double* TRI_Vertices (TRI *tri, int n, int *m)
{
  int vcnt; /* number of vertices */
  MAP *vm, *im; /* map of visited vertices; iterator 'im' */
  MEM mem; /* memory for map items */
  double *v, *w, *z;
  TRI *t, *e;
  int i;

  MEM_Init (&mem, sizeof (MAP), n * 3);
  e = tri + n;
  vm = NULL;
  vcnt = 0;

  for (t = tri; t < e; t ++)
  {
    for (i = 0; i < 3; i ++)
    {
      if (! MAP_Find_Node (vm, t->ver [i], NULL)) /* vertex not mapped before */
      {
	MAP_Insert (&mem, &vm, t->ver [i], (void*) (long) vcnt, NULL); /* map vertex 'v' to a consecutive index */
	vcnt ++;
      }
    }
  }

  /* alloc output memory */
  ERRMEM (v = malloc (sizeof (double [3]) * vcnt));

  /* map old vertices into new placeholders */
  for (im = MAP_First (vm); im; im = MAP_Next (im), i ++)
  {
    w = (double*) im->key; /* vertex source */
    i = (int) (long) im->data; /* vertex index */
    z = v + i*3; /* vertex destination */
    COPY (w, z);
  }

  MEM_Release (&mem);

  *m = vcnt;
  return v;
}

/* compute planes */
double* TRI_Planes (TRI *tri, int n, int *m)
{
  double *v, *w;
  TRI *t, *e;

  ERRMEM (v = malloc (sizeof (double [6]) * n));
  e = tri + n;
  
  for (t = tri, w = v; t < e; t ++, w += 6)
  {
    COPY (t->out, w);
    COPY (t->ver[0], w+3);
  }

  *m = n;
  return v;
}

/* compute mass center and volume of triangulated solid */
double TRI_Char (TRI *tri, int n, double *center)
{
  double zero [3] = {0, 0, 0},
	 volume, sx, sy, sz,
	 J, *a, *b, *c;
  TRI *t, *e;

  volume = sx = sy = sz = 0.0;
  for (t = tri, e = t+n; t != e; t++)
  {
    a = t->ver [0];
    b = t->ver [1];
    c = t->ver [2];

    J = simplex_J (zero, a, b, c);
    volume += simplex_1 (J, zero, a, b, c);

    if (center)
    {
      sx += simplex_x (J, zero, a, b, c);
      sy += simplex_y (J, zero, a, b, c);
      sz += simplex_z (J, zero, a, b, c);
    }
  }

  if (center)
  {
    if (volume > 0.0)
    {
      center [0]  = sx / volume;
      center [1]  = sy / volume;
      center [2]  = sz / volume;
    }
    else if (n)
    {
      SET (center, 0.0);

      for (t = tri, e = t+n; t != e; t++)
      {
	a = t->ver [0];
	b = t->ver [1];
	c = t->ver [2];
	MID3 (a, b, c, zero);
        ADD (center, zero, center);
      }

      DIV (center, (double)n, center);
    }
  }

  return volume;
}
