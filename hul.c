/*
 * hul.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * convex hull in three dimensions according to the algorithm by
 * Barber et al. "The Quickhull Algorithm for Convex Hulls"
 * ACM Transactions on Mathematical Software, 1996
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
#include <float.h>
#include "mem.h"
#include "err.h"
#include "alg.h"
#include "lis.h"
#include "hul.h"

typedef struct vertex vertex;
typedef struct edge edge;
typedef struct face face;

struct vertex
{
  double *v; /* input vertex */
  vertex *n; /* next in a list */
};

struct edge
{
  double *v [2]; /* edge vertices in CCW order */
  face *f; /* a neighbouring face through this edge */
  edge *n; /* next edge in a list */
};

struct face
{
  double pla [4];
  vertex *v, *w; /* list of facial vertices 'v' and the furthest vertex 'w' */
  edge *e; /* list of edges (and implicitly, the list of neighbours */
  face *n; /* next face in a list */
  char marked; /* marker used for visible faces */
  TRI *tri; /* auxiliary adjacent triangle (used to create the output table) */
};

/* set up facial plane */
static int setplane (face *f)
{
  double ba [3], cb [3];
  edge *e1, *e2;

  e1 = f->e; e2 = e1->n;
  SUB (e1->v[1], e1->v[0], ba);
  SUB (e2->v[1], e2->v[0], cb);
  PRODUCT (ba, cb, f->pla);
  f->pla [3] = - DOT (e1->v[0], f->pla);

  MAXABS (f->pla, ba [0]);
  if (ba [0] == 0.0) return 0; /* degenerate case (colinear vertices) */
  else return 1;
}

/* return edge of f pointing to g */
inline static edge* othere (face *f, face *g)
{
  for (edge *e = f->e; e; e = e->n)
  {
    if (e->f == g) return e;
  }

#if GEOMDEBUG
  ASSERT_DEBUG (0, "Inconsitent topology (edge not found)");
#else
  return NULL;
#endif
}

/* return vertex of f other than p and q */
inline static double* otherv (face *f, double *p, double *q)
{
  for (edge *e = f->e; e; e = e->n)
  {
    if (e->v[0] != p && e->v[0] != q) return e->v [0];
  }

#if GEOMDEBUG
  ASSERT_DEBUG (0, "Inconsitent topology (vertex not found)");
#else
  return NULL;
#endif
}

/* return edge starting at v */
inline static edge* edge_0 (face *f, double *v)
{
  for (edge *e = f->e; e; e = e->n)
  {
    if (e->v [0] == v) return e;
  }

#if GEOMDEBUG
  ASSERT_DEBUG (0, "Inconsitent topology (edge not found)");
#else
  return NULL;
#endif
}

/* return edge ending at v */
inline static edge* edge_1 (face *f, double *v)
{
  for (edge *e = f->e; e; e = e->n)
  {
    if (e->v [1] == v) return e;
  }

#if GEOMDEBUG
  ASSERT_DEBUG (0, "Inconsitent topology (edge not found)");
#else
  return NULL;
#endif
}

/* mend face with colinear vertices */
static int mendface (face *f)
{
  edge *e, *h, *x, *y, *z, *w, *o;
  double u [3], v [3], d, l [2];
  double *a, *b, *c;
  face *g;

  for (e = f->e; e; e = e->n)
  {
    a = e->v [0];
    b = e->v [1];
    if (!(c = otherv (f, a, b))) return 0;

    SUB (b, a, u);
    SUB (c, a, v);

    l[0] = LEN (u);
    l[1] = LEN (v);

    d = DOT (u, v);

    if (d >= 0 && l[1] <= l[0]) break; /* c in [a, b] */
  }

  if (e)
  {
    a = c;
    g = e->f;
    if (!(h = othere (g, f))) return 0;
    if (!(b = otherv (g, h->v[0], h->v[1]))) return 0;

    if (!(x = edge_0 (f, a))) return 0;
    if (!(y = edge_1 (g, b))) return 0;
    if (!(z = edge_0 (g, b))) return 0;
    if (!(w = edge_1 (f, a))) return 0;

    if (!(o = othere (y->f, g))) return 0;
    o->f = f;
    if (!(o = othere (w->f, f))) return 0;
    o->f = g;

    e->v [0] = b;
    e->v [1] = a;
    h->v [0] = a;
    h->v [1] = b;

    e->n = NULL;
    y->n = e;
    x->n = y;
    f->e = x;

    h->n = NULL;
    w->n = h;
    z->n = w;
    g->e = z;

#if GEOMDEBUG
    ASSERT_DEBUG (setplane (f), "Zero normal (when mending face)");
    ASSERT_DEBUG (setplane (g), "Zero normal (when mending face)");
#else
    if (!(setplane (f) && setplane (g))) return 0;
#endif
  }
#if GEOMDEBUG
  else
  {
    ASSERT_DEBUG (0, "Face mending failed");
  }
#else
  else return 0;
#endif

  return 1;
}


/* vertex least equal comparison */
static int vertex_le (vertex *a, vertex  *b)
{
  if ((a->v[0]+GEOMETRIC_EPSILON) < (b->v[0]-GEOMETRIC_EPSILON)) return 1;
  else if (fabs (a->v[0]-b->v[0]) <= (2.0*GEOMETRIC_EPSILON))
  { 
    if ((a->v[1]+GEOMETRIC_EPSILON) < (b->v[1]-GEOMETRIC_EPSILON)) return 1;
    else if (fabs (a->v[1]-b->v[1]) <= (2.0*GEOMETRIC_EPSILON))
    {
     if ((a->v[2]+GEOMETRIC_EPSILON) < (b->v[2]-GEOMETRIC_EPSILON)) return 1;
     else if (fabs (a->v[2]-b->v[2]) <= (2.0*GEOMETRIC_EPSILON)) return 1;
    }
  }

  return 0;
}

/* vertex equal comparison */
static int vertex_eq (vertex *a, vertex  *b)
{
  if ((fabs (a->v[0]-b->v[0]) <= (2.0*GEOMETRIC_EPSILON)) &&
      (fabs (a->v[1]-b->v[1]) <= (2.0*GEOMETRIC_EPSILON)) &&
      (fabs (a->v[2]-b->v[2]) <= (2.0*GEOMETRIC_EPSILON))) return 1;

  return 0;
}

/* lexicographical list sort of verex list */
IMPLEMENT_LIST_SORT (SINGLE_LINKED, sort_vertices, vertex, p, n, vertex_le)

/* select vertices of an initial simplex and output the list of remaining vertices */
static vertex* simplex_vertices (double *v, int n, MEM *mv, double *sv [4])
{
  double *e = v+3*n, *w, a [3], b [3], c [3], u [3], d;
  vertex *o, *x, *y;
  int i, j;

  /* find well separated vertices */
  for (w = v, j = 0; w < e && j < 4; w += 3)
  {
    for (i = 0; i < j; i ++)
    {
      SUB (sv [i], w, u);
      MAXABS (u, d);
      if (d < 2.0 * GEOMETRIC_EPSILON) break; /* too close */
    }

    if (i == j) /* not too close, but maybe ... */
    {
      switch (j)
      {
	case 2:
	{
	  SUB (sv [1], sv [0], a);
	  PRODUCT (u, a, b);
	  MAXABS (b, d);
          if (d < GEOMETRIC_EPSILON) continue; /* ... colinear */
	}
	break;
	case 3:
	{
	  SUB (sv [1], sv [0], a);
	  SUB (sv [2], sv [1], b);
	  PRODUCT (a, b, c);
	  d = DOT (u, c);
          if (ABS (d) < GEOMETRIC_EPSILON) continue; /* ... coplanar */
	}
	break;
      }
      
      sv [j++] = w; /* add vertex */
    }
  }

#if GEOMDEBUG
  ASSERT_DEBUG (j == 4, "All input points coincide");
#else
  if (j != 4) return NULL;
#endif

  /* copy the remaining vertices */
  for (o = NULL, w = v; w < e; w += 3)
  {
    for (i = 0; i < 4; i ++)
    { if (sv [i] == w) break; } /* skip already selected */

    if (i == 4) /* this one was not selected */
    {
      ERRMEM (x = MEM_Alloc (mv));
      x->v = w;
      x->n = o; /* add to list */
      o = x; /* ... */
    }
  }

  o = sort_vertices (o); /* lexicographicaly sort vertices */

  /* remove duplicates */
  for (x = o, o = NULL; x; x = y)
  {
    for (y = x->n; y; y = y->n) /* skip duplicates within the list */
    {
      if (!vertex_eq (x, y)) break;
    }

    for (i = 0; i < 4; i ++) /* skip x if it dupilcates one of the output vertices */
    {
      SUB (x->v, sv [i], u);
      MAXABS (u, d);
      if (d <= 2.0 * GEOMETRIC_EPSILON) break;
    }

    if (i == 4) /* does not duplicate output vertices */
    {
      x->n = o;
      o = x;
    }
  }

  return o;
}

/* create a simplex and return the corresponding face list */
static face* simplex (MEM *me, MEM *mf, double *a, double *b, double *c, double *d)
{
  edge *e [12] = 
  { MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me),
    MEM_Alloc (me) }, *edg;
  
  double *o [4] = {d, a, b, c}, *u;

  int i;

  face *f [4] = 
  { MEM_Alloc (mf),
    MEM_Alloc (mf),
    MEM_Alloc (mf),
    MEM_Alloc (mf) };

  ERRMEM (e[0] && e[1] && e[2]  && e[3]  &&
          e[4] && e[5] && e[6]  && e[7]  &&
	  e[8] && e[9] && e[10] && e[11] &&
          f[0] && f[1] && f[2]  && f[3]);

  e [0]->v [0] = a;
  e [0]->v [1] = b;
  e [1]->v [0] = b;
  e [1]->v [1] = c;
  e [2]->v [0] = c;
  e [2]->v [1] = a;
  f [0]->e = e [0];
  e [0]->n = e [1];
  e [1]->n = e [2];

  e [3]->v [0] = b;
  e [3]->v [1] = d;
  e [4]->v [0] = d;
  e [4]->v [1] = c;
  e [5]->v [0] = c;
  e [5]->v [1] = b;
  f [1]->e = e [3];
  e [3]->n = e [4];
  e [4]->n = e [5];

  e [6]->v [0] = d;
  e [6]->v [1] = a;
  e [7]->v [0] = a;
  e [7]->v [1] = c;
  e [8]->v [0] = c;
  e [8]->v [1] = d;
  f [2]->e = e [6];
  e [6]->n = e [7];
  e [7]->n = e [8];

  e [9 ]->v [0] = d;
  e [9 ]->v [1] = b;
  e [10]->v [0] = b;
  e [10]->v [1] = a;
  e [11]->v [0] = a;
  e [11]->v [1] = d;
  f [3]->e = e [9 ];
  e [9]->n = e [10];
  e[10]->n = e [11];

  e [0]->f = f [3];
  e [1]->f = f [1];
  e [2]->f = f [2];

  e [3]->f = f [3];
  e [4]->f = f [2];
  e [5]->f = f [0];

  e [6]->f = f [3];
  e [7]->f = f [0];
  e [8]->f = f [1];

  e [9 ]->f = f [1];
  e [10]->f = f [0];
  e [11]->f = f [2];

  f [0]->n = f [1];
  f [1]->n = f [2];
  f [2]->n = f [3];

  for (i = 0; i < 4; i ++)
  {
#if GEOMDEBUG
    ASSERT_DEBUG (setplane (f[i]), "Zero normal (when creating initial simplex)"); /* set face plane */
#else
    if (!setplane (f[i])) return NULL; /* set face plane */
#endif
    if (PLANE (f[i]->pla, o [i]) > GEOMETRIC_EPSILON) /* the other vertex should be behind => reorient the face */
    {
      for (edg = f[i]->e; edg; edg = edg->n) /* reverse order of edge vertices */
      { u = edg->v [0]; edg->v [0] = edg->v [1]; edg->v [1] = u; }
      e [0] = f [i]->e; /* reverse edge list */
      e [1] = e [0]->n;
      e [2] = e [1]->n;
      e [2]->n = e [1];
      e [1]->n = e [0];
      e [0]->n = NULL;
      f [i]->e = e [2];
      SCALE (f[i]->pla, -1.0); /* reverse plane normal */
      f[i]->pla [3] *= -1.0; 
    }
  }

  return f [0];
}

/* mark faces visible from 'v'ertex */
static void mark (face *f, double *v, face **g)
{
  double d = PLANE (f->pla, v);
    
  if (!f->marked && d > 0.0)
  {
    f->marked = 1;
    for (edge *e = f->e; e; e = e->n) mark (e->f, v, g);
  }
  else if (d <= 0.0) *g = f;
}

/* return next CCW face after f around vertx v */
inline static face* nextaround (face *f, double *v)
{
  edge *e;

  for (e = f->e; e; e = e->n)
    if (e->v [1] == v) return e->f;

#if GEOMDEBUG
  ASSERT_DEBUG (0, "Inconsitent adjacency (in nextaround)");
#endif

  return NULL;
}

/* walk behind the horizon (unvisible side)
 * ridges and return consecutive CCW edges */
inline static edge* nextonridge (edge* e, face **g)
{
  if (g)
  {
    double *v = e->v[1];
    face *f = e->f;

    for (f = nextaround (f, v); f && f->marked; f = nextaround (f, v)); /* walk around e->v[0] until unmarked face is found */

#if GEOMDEBUG
    ASSERT_DEBUG (f, "Inconsitent topology => first edge on the ridge not found (1 in nextonridge)");
#else
    if (!f) return NULL;
#endif

    for (e = f->e; e && e->v [0] != v; e = e->n); /* find an edge adjacent to the marked region */

    *g = f; /* record new unmarked face */
  }
  else
  {
    for (; e && (!e->f->marked); e = e->n); /* should be there for the first call */
#if GEOMDEBUG
    ASSERT_DEBUG (e, "Inconsitent topology => first edge on the ridge not found (2 in nextonridge)");
#endif
  }

  return e;
}

static int testsimplex (face *h)
{
  edge *e;
  face *f;
  int i, n; 

  for (; h; h = h->n)
  {
    for (e = h->e; e; e = e->n)
    {
      for (i = 0; i < 2; i ++)
      {
	for (n = 1, f = nextaround (h, e->v[i]); f && f != h && n < 4; f = nextaround (f, e->v[i])) n ++;
#if GEOMDEBUG
	ASSERT_DEBUG (f && n == 3, "Incorrect adjacency in the initial simplex");
#else
	if (!f || n != 3) return 0;
#endif
      }
    }
  }

  return 1;
}

/* compute convex hull */
TRI* hull (double *v, int n, int *m)
{
  face *f, *g, *h, *u, *head, *cur, *tail;
  edge *e, *k, *i, *j, *ehead, *etail;
  double d, dmax, *sv [4];
  vertex *x, *y, *z, *l;
  MEM mv, me, mf;
  TRI *tri, *t;

  MEM_Init (&mv, sizeof (vertex), n);
  MEM_Init (&me, sizeof (edge), n);
  MEM_Init (&mf, sizeof (face), n);
  tri = NULL;

  /* select vertices of an initial simplex into 'sv' */
  if (!(l = simplex_vertices (v, n, &mv, sv))) goto error;
 
  /* create the initial simplex */ 
  if (!(h = simplex (&me, &mf, sv[0], sv[1], sv[2], sv[3]))) goto error;

  if (!(testsimplex (h))) goto error;

  /* initialise outside vertex lists */
  for (f = h; f; f = f->n)
  {
    for (z = NULL, x = l, dmax = 0.0; x; x = y)
    {
      y = x->n; 
      d = PLANE (f->pla, x->v);
      if (d > GEOMETRIC_EPSILON)
      {
	if (d > dmax) /* and select maximal elements */
	{
	  if (f->w)
	  { f->w->n = f->v;
	    f->v = f->w; } /* move to the regular list */
	  f->w = x; /* set as maximal */
	}
	else /* insert into the regular list */
	{ x->n = f->v;
	  f->v = x; }

        /* update 'l' list */	  
	if (z) z->n = y; /* skip one */
	else l = y; /* update head */
      }
      else z = x; /* previous element staying in the list */
    }
  }

  for (f = h; f && (!f->w); f = f->n); /* find first face with a nonempty vertex list */
  while (f)
  {
    /* mark visible faces */
    mark (f, f->w->v, &g);

    /* loop over the ridge edges */
    if (!(k = e = nextonridge (g->e, NULL))) goto error;
    ehead = etail = NULL;
    head = tail = NULL;
    do
    {
      /* create new face */
      ERRMEM (cur = MEM_Alloc (&mf));
      if (!tail) tail = cur; /* record last face */
      ERRMEM (i = MEM_Alloc (&me));
      i->v [0] = e->v [1]; /* first new edge is adjacent to 'e' => reversed */
      i->v [1] = e->v [0];
      i->f = g; /* first new edge is the neighbour of 'g' */
      cur->e = i; /* include the edge into the new face's edge list */
      ERRMEM (j = MEM_Alloc (&me));
      j->v [0] = f->w->v; /* this is the top vertes */
      j->v [1] = e->v [1];
      j->n = cur->e; cur->e = j; /* maintain edge list */
      ERRMEM (i = MEM_Alloc (&me));
      i->v [0] = e->v [0];
      i->v [1] = f->w->v; /* top vertex */
      i->n = cur->e; cur->e = i; /* maintain edge list */

      if (head) /* if there are already new faces in the list */
      { i->f = head; /* this edge's neighbour is the list head */
        ehead->f = cur; } /* and head's edge neighbour is the current face */
      else etail = i; /* or => set up tail's edge (to be connected at the end) */

      cur->n = head; head = cur; /* maintain face list */
      ehead = j; /* this is the head edge */
      
      j = e; /* back up current outer edge => 'nextonridge' needs an old 'e->f' */
      if (!(e = nextonridge (j, &g))) goto error; /* next outer edge along the visible set ridge */
      j->f = cur; /* set up new adjacency (once the old 'e->f' was utilised) */

    } while (e != k);
    /* link tail and head */
    ehead->f = tail;
    etail->f = head;

    /* free top vertex */
    MEM_Free (&mv, f->w);
    f->w = NULL;

    /* for each new face */
    for (g = head; g; g = g->n)
    {
      if (!setplane (g)) /* set up g->pla (returnes 0 if a degenerate triangle was found)  */
      {
        if (!mendface (g)) goto error; /* vertices are colinear but not coincident (that case was eliminated by sorting and filtering) */
      }
    }

    /* for each marked face f */
    for (f = h; f; f = f->n)
    {
      if (f->marked)
      {
	if (f->w)
	{ f->w->n = f->v; /* put the furthest vertex 'w' back into the 'v' list */
	  f->v = f->w; }

	if (f->v)
	{
	  /* for each new face g */
	  for (g = head; g; g = g->n)
	  {
	    /* for each v in f->v */
	    for (dmax = 0.0, l = NULL, x = f->v; x; x = y)
	    {
	      y = x->n;

	      d = PLANE (g->pla, x->v);
	      if (d > GEOMETRIC_EPSILON) /* x is above g->pla */
	      {
		if (d > dmax) /* and is maximal */
		{
		  if (g->w)
		  { g->w->n = g->v;
		    g->v = g->w; } /* move to the regular list */
		  g->w = x; /* set as maximal */
		}
		else /* insert into the regular list */
		{ x->n = g->v;
		  g->v = x; }

		if (l) l->n = y; /* x is removed from f->v */
		else f->v = y;
	      }
	      else l = x; /* last not moved vertex */
	    }
	  }
	}
      }
    }

    /* for each marked face f */
    for (u = NULL, f = h; f; f = g)
    {
      g = f->n;
      if (f->marked)
      {
        /* delete all v in f->v */
	for (x = f->v; x; x = y)
	{ y = x->n; MEM_Free (&mv, x); }

        /* delete all e in f->e */
	for (e = f->e; e; e = i)
	{ i = e->n; MEM_Free (&me, e); }

        /* delete f */
	MEM_Free (&mf, f);

	/* update 'h' list */
	if (u) u->n= g; /* skip */
	else h = g; /* new head */
      }
      else u = f; /* previous unmarked face */
    }

    /* append h with
     * the new faces */
    tail->n = h;
    h = head;
   
    /* select next face with nonempty vertex list */ 
    for (f = h; f && (!f->w); f = f->n);
  }

  /* h contains faces of the convex hull;
   * it be now translated into a table TRI[] */

  for ((*m) = 0, f = h; f; f = f->n) (*m) ++; /* count output faces */
  ERRMEM (tri = calloc ((*m), sizeof (TRI))); /* output memory (faces are triangular) */
  for (t = tri, f = h; f; f = f->n, t ++) /* translate each face into a triangle */
  {
    e = f->e; k = e->n; i = k->n;
    COPY (f->pla, t->out); /* same normal */
    t->ver [0] = e->v [0]; /* CCW ordered vertices */
    t->ver [1] = k->v [0];
    t->ver [2] = i->v [0];
#if GEOMDEBUG
    if (e->f->tri) { ASSERT_DEBUG (TRI_Addadj (e->f->tri, t), "Too many triangle neighbours"); } /* called only once for each pair => after (***) ... */
    if (k->f->tri) { ASSERT_DEBUG (TRI_Addadj (k->f->tri, t), "Too many triangle neighbours"); }
    if (i->f->tri) { ASSERT_DEBUG (TRI_Addadj (i->f->tri, t), "Too many triangle neighbours"); }
#else
    if (e->f->tri) if (!TRI_Addadj (e->f->tri, t)) goto error; /* called only once for each pair => after (***) ... */
    if (k->f->tri) if (!TRI_Addadj (k->f->tri, t)) goto error;
    if (i->f->tri) if (!TRI_Addadj (i->f->tri, t)) goto error;
#endif
    f->tri = t; /* ... (***) has been executed for the first of neighbours */
  }

  for (t --; t >= tri; t --) TRI_Sortadj (t); /* sort adjacency lists */
  
  goto done; /* skip error handling */

error:
 if (tri) free (tri); 
 tri = NULL;

done:
  /* clean up */
  MEM_Release (&mv);
  MEM_Release (&me);
  MEM_Release (&mf);

  return tri;
}
