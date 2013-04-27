/*
 * msh.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * mesh based shape
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
#include "ext/tetgen/tetsol.h"
#include "sol.h"
#include "spx.h"
#include "err.h"
#include "mem.h"
#include "map.h"
#include "alg.h"
#include "msh.h"
#include "pck.h"
#include "kdt.h"
#include "gjk.h"
#include "pbf.h"

/* flag used during MESH_Split redesign */
#define OLD_SPLIT 1 /* FIXME */

/* used in some pools */
#define MEMCHUNK 128

/* linear shape functions for the hexahedron */
#define HEX0(x,y,z) (0.125*(1.0-(x))*(1.0-(y))*(1.0-(z)))
#define HEX1(x,y,z) (0.125*(1.0+(x))*(1.0-(y))*(1.0-(z)))
#define HEX2(x,y,z) (0.125*(1.0+(x))*(1.0+(y))*(1.0-(z)))
#define HEX3(x,y,z) (0.125*(1.0-(x))*(1.0+(y))*(1.0-(z)))
#define HEX4(x,y,z) (0.125*(1.0-(x))*(1.0-(y))*(1.0+(z)))
#define HEX5(x,y,z) (0.125*(1.0+(x))*(1.0-(y))*(1.0+(z)))
#define HEX6(x,y,z) (0.125*(1.0+(x))*(1.0+(y))*(1.0+(z)))
#define HEX7(x,y,z) (0.125*(1.0-(x))*(1.0+(y))*(1.0+(z)))

static int tet [][4] =  /* 1-based indexing as the element lists start with the number of nodes */
  {{1,3,2,0},
   {1,2,4,0},
   {2,3,4,0},
   {3,1,4,0}}, pyr [][4] = 
  {{1,4,3,2},
   {1,2,5,0},
   {2,3,5,0},
   {3,4,5,0},
   {4,1,5,0}}, wed [][4] =
  {{1,3,2,0},
   {4,5,6,0},
   {1,2,5,4},
   {2,3,6,5},
   {3,1,4,6}}, hex [][4] =
  {{1,4,3,2},
   {1,2,6,5},
   {2,3,7,6},
   {3,4,8,7},
   {1,5,8,4},
   {5,6,7,8}};

/* maximal number of neighbours */
inline static int neighs (int type)
{
  switch (type)
  {
    case 4: return 4;
    case 5: return 5;
    case 6: return 5;
    case 8: return 6;
  }
  return 0;
}

inline static void swap (int *a, int *b)
{
  int c = *a; *a = *b; *b = c;
}

/* quick sort of ints */
static void sort (int* begin, int* end)
{
  int *lower = begin, *upper = end,
    bound = *(begin+(end-begin)/2);
  
  while (lower <= upper)
  {
    while (*lower < bound) lower++;
    while (bound < *upper) upper--;

    if (lower < upper) swap (lower++, upper--);
    else lower++;
  }

  if (begin < upper) sort (begin, upper);
  if (upper < end) sort (upper+1, end);
}

static int lexcmp (int *a, int *b, int m)
{
  int n;
  for (n = 0; n < m; n ++)
  {
    if (a [n] < b [n]) return -1;
    else if (a [n] > b [n]) return 1;
  }
  return 0;
}

/* comparison used in face mapping */
static int face_compare (FACE *one, FACE *two)
{
  if (one->type < two->type) return -1;
  else if (one->type > two->type) return 1;
  else return lexcmp (one->nodes, two->nodes, one->type);
}

static ELEMENT* create_element (MEM *elemem, int *element)
{
  ELEMENT *ele;
  int n;

  ele = MEM_Alloc (elemem);
  ele->type = element [0];
  for (n = 1; n <= element [0]; n ++)
    ele->nodes [n-1] = element [n];
  ele->volume = element [n]; /* volume identifier */

  return ele;
}

static void setup_face (ELEMENT *ele, int n, FACE *fac, int dosort)
{
  switch (ele->type)
  {
    case 4:
      fac->type = 3;
      fac->nodes [0] = ele->nodes [tet [n][0]-1]; /* shift due to the 1-based indexing */
      fac->nodes [1] = ele->nodes [tet [n][1]-1];
      fac->nodes [2] = ele->nodes [tet [n][2]-1];
      if (dosort) sort (fac->nodes, fac->nodes+2);
    break;
    case 5:
    if (n == 0)
    { fac->type = 4;
      fac->nodes [0] = ele->nodes [pyr [n][0]-1];
      fac->nodes [1] = ele->nodes [pyr [n][1]-1];
      fac->nodes [2] = ele->nodes [pyr [n][2]-1];
      fac->nodes [3] = ele->nodes [pyr [n][3]-1];
      if (dosort) sort (fac->nodes, fac->nodes+3); }
    else
    { fac->type = 3;
      fac->nodes [0] = ele->nodes [pyr [n][0]-1];
      fac->nodes [1] = ele->nodes [pyr [n][1]-1];
      fac->nodes [2] = ele->nodes [pyr [n][2]-1];
      if (dosort) sort (fac->nodes, fac->nodes+2); }
    break;
    case 6:
    if (n < 2)
    { fac->type = 3;
      fac->nodes [0] = ele->nodes [wed [n][0]-1];
      fac->nodes [1] = ele->nodes [wed [n][1]-1];
      fac->nodes [2] = ele->nodes [wed [n][2]-1];
      if (dosort) sort (fac->nodes, fac->nodes+2); }
    else
    { fac->type = 4;
      fac->nodes [0] = ele->nodes [wed [n][0]-1];
      fac->nodes [1] = ele->nodes [wed [n][1]-1];
      fac->nodes [2] = ele->nodes [wed [n][2]-1];
      fac->nodes [3] = ele->nodes [wed [n][3]-1];
      if (dosort) sort (fac->nodes, fac->nodes+3); }
    break;
    case 8:
      fac->type = 4;
      fac->nodes [0] = ele->nodes [hex [n][0]-1];
      fac->nodes [1] = ele->nodes [hex [n][1]-1];
      fac->nodes [2] = ele->nodes [hex [n][2]-1];
      fac->nodes [3] = ele->nodes [hex [n][3]-1];
      if (dosort) sort (fac->nodes, fac->nodes+3);
    break;
  }
}

/* copy face[n] vertex indices into 'fac'; return shifted 'fac' pointer */
static int* setup_face_vertices (ELEMENT *ele, int n, int *fac)
{
  switch (ele->type)
  {
    case 4:
      fac [0] = 3; /* number of vertices */
      fac [1] = tet [n][0]-1; /* shift due to the 1-based indexing */
      fac [2] = tet [n][1]-1;
      fac [3] = tet [n][2]-1;
      return (fac + 4);
    break;
    case 5:
    if (n == 0)
    { fac [0] = 4;
      fac [1] = pyr [n][0]-1;
      fac [2] = pyr [n][1]-1;
      fac [3] = pyr [n][2]-1;
      fac [4] = pyr [n][3]-1;
      return (fac + 5); }
    else
    { fac [0] = 3;
      fac [1] = pyr [n][0]-1;
      fac [2] = pyr [n][1]-1;
      fac [3] = pyr [n][2]-1;
      return (fac + 4); }
    break;
    case 6:
    if (n < 2)
    { fac [0] = 3;
      fac [1] = wed [n][0]-1;
      fac [2] = wed [n][1]-1;
      fac [3] = wed [n][2]-1;
      return (fac + 4); }
    else
    { fac [0] = 4;
      fac [1] = wed [n][0]-1;
      fac [2] = wed [n][1]-1;
      fac [3] = wed [n][2]-1;
      fac [4] = wed [n][3]-1;
      return (fac + 5); }
    break;
    case 8:
      fac [0] = 4;
      fac [1] = hex [n][0]-1;
      fac [2] = hex [n][1]-1;
      fac [3] = hex [n][2]-1;
      fac [4] = hex [n][3]-1;
      return (fac + 5);
    break;
  }
  return NULL;
}


static void setup_normal (double (*nodes) [3], FACE *fac)
{
  int n0, n1, n2;
  double *normal;

  n0 = fac->nodes [0];
  n1 = fac->nodes [1];
  n2 = fac->nodes [2];
  normal = fac->normal;
  NORMAL (nodes [n0], nodes [n1], nodes [n2], normal);
  NORMALIZE (normal);
}

/* get new element, old face list and create the element faces => return the new face list */
static FACE* create_faces (MEM *facmem, MEM *mapmem, MAP **faces, ELEMENT *ele, FACE *list)
{
  FACE *fac, tmp;
  int n, m;

  m = neighs (ele->type); 
      
  for (n = 0; n < m; n ++)
  {
    /* set up temporary face for to
     * be used as the map search key */
    setup_face (ele, n, &tmp, 1); /* nodes sorted for map key comparisons */
    fac = MAP_Find (*faces, &tmp, (MAP_Compare) face_compare);  /* is it there ? */

    if (fac) /* was mapped already */
    {
      /* set up element adjacency */
      ele->adj [ele->neighs] = fac->ele;
      fac->ele->adj [fac->ele->neighs] = ele;
      fac->ele->neighs ++;
      ele->neighs ++;

      fac->ele = NULL; /* mark as the inner face (***) */
    }
    else
    {
      fac = MEM_Alloc (facmem);
      fac->ele = ele;
      setup_face (ele, n, fac, 1);
      fac->index = n; /* local index */
      MAP_Insert (mapmem, faces, fac, /* map by the type/nodes key */
	fac, (MAP_Compare) face_compare);
      fac->next = list;
      list = fac;
    }
  }

  return list;
}

static int maximal (int *element)
{
  int n, ret = 0;
  for (n = 1; n <= element [0]; n++)
    if (element [n] > ret) ret = element [n];
  return ret;
}

static int minimal (int *element)
{
  int n, ret = INT_MAX;
  for (n = 1; n <= element [0]; n++)
    if (element [n] < ret) ret = element [n];
  return ret;
}

/* copy node coordinates into a local table => for cache efficiency */
static void load_nodes (double (*heap) [3], int type, int *nodes, double (*stack) [3])
{
  int n;

  for (n = 0; n < type; n ++)
  { COPY (heap [nodes [n]], stack [n]); }
}

#if 0
/* create face planes from the current shape of an element */
static int create_element_planes (double (*node) [3], ELEMENT *ele, double *pla)
{
  double a [3], b [3], c [3];
  int *n = ele->nodes - 1, /* due to one based indexing in tet, pyr, wed, hex */
      (*f) [4], i, m;

  switch (ele->type)
  { case 4: f = tet; break; /* set up face definitions */
    case 5: f = pyr; break;
    case 6: f = wed; break;
    case 8: f = hex; break; }

  m = neighs (ele->type);
  for (i = 0; i < m; i ++, pla += 4)
  {
    COPY (node [n[f[i][0]]], a);
    COPY (node [n[f[i][1]]], b);
    COPY (node [n[f[i][2]]], c);
    NORMAL (a, b, c, pla);
    NORMALIZE (pla);
    pla [3] = - DOT (a, pla);
  }

  return m;
}

/* point behind all planes query */
static int point_inside (int npla, double *pla, double *point)
{
  for (; npla > 0; pla += 4, npla --)
    if (PLANE (pla, point) > GEOMETRIC_EPSILON) return 0;

  return 1;
}

/* maximal plane distance query */
static double point_distance (int npla, double *pla, double *point)
{
  double d, dist;

  for (dist = -DBL_MAX; npla > 0; pla += 4, npla --)
  {
    d = PLANE (pla, point);
    if (d > dist) dist = d;
  }

  return dist;
}

/* check if edge [nod1, nod2] belongs to the element */
static int element_has_edge (ELEMENT *ele, int nod1, int nod2)
{
  int *nodes = ele->nodes,
      type = ele->type, n, j;

  for (n = j = 0; n < type; n ++)
    if (nodes [n] == nod1 ||
	nodes [n] == nod2) j ++;

  return (j == 2 ? 1 : 0);
}
#endif

/* compute element extents */
inline static void element_extents (MESH *msh, ELEMENT *ele, int ref, double *extents)
{
  double nodes [8][3];
  int n;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;
  
  load_nodes (ref ? msh->ref_nodes : msh->cur_nodes, ele->type, ele->nodes, nodes);

  for (n = 0; n < ele->type; n ++)
  {
    if (nodes [n][0] < extents [0]) extents [0] = nodes [n][0];
    if (nodes [n][1] < extents [1]) extents [1] = nodes [n][1];
    if (nodes [n][2] < extents [2]) extents [2] = nodes [n][2];
    if (nodes [n][0] > extents [3]) extents [3] = nodes [n][0];
    if (nodes [n][1] > extents [4]) extents [4] = nodes [n][1];
    if (nodes [n][2] > extents [5]) extents [5] = nodes [n][2];
  }

  extents [0] -= GEOMETRIC_EPSILON;
  extents [1] -= GEOMETRIC_EPSILON;
  extents [2] -= GEOMETRIC_EPSILON;
  extents [3] += GEOMETRIC_EPSILON;
  extents [4] += GEOMETRIC_EPSILON;
  extents [5] += GEOMETRIC_EPSILON;
}


/* check if element has nodes */
static int element_has_nodes (ELEMENT *ele, int n, int *nodes)
{
  int *e = ele->nodes, m = ele->type, k = 0, i, j;

  for (i = 0; i < m; i ++)
    for (j = 0; j < n; j ++)
      if (e [i] == nodes [j]) k ++;

  return (k == n ? 1 : 0); 
}

#if OLD_SPLIT
/* face triangulation contact test */
static int face_touches_triangulation (double *v [4], int n, double *normal, KDT *kd)
{
  double extents [6], a [16], b [12], p [3], q [3], deps;
  SET *leaves = NULL, *item;
  KDT *leaf;
  TRI *t;
  int i;

  COPY (v [0], extents);
  COPY (extents, extents+3);
  for (i = 1; i < n; i ++)
  {
    if (v [i][0] < extents [0]) extents [0] = v [i][0];
    else if (v [i][0] > extents [3]) extents [3] = v [i][0];
    if (v [i][1] < extents [1]) extents [1] = v [i][1];
    else if (v [i][1] > extents [4]) extents [4] = v [i][1];
    if (v [i][2] < extents [2]) extents [2] = v [i][2];
    else if (v [i][2] > extents [5]) extents [5] = v [i][2];
  }

  KDT_Pick_Extents (kd, extents, &leaves);
  deps = 10 * GEOMETRIC_EPSILON;

  COPY (v [0], a);
  ADDMUL (a, deps, normal, a+3);
  COPY (v [1], a+6);
  COPY (v [2], a+9);
  if (n == 4) { COPY (v [3], a+12); }

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
      if (gjk (a, n, b, 4, p, q) < GEOMETRIC_EPSILON)
      {
	SET_Free (NULL, &leaves);
	return 1;
      }
    }
  }

  SET_Free (NULL, &leaves);
  return 0;
}

/* element triangulation contact test */
static int element_touches_triangulation (MESH *msh, ELEMENT *ele, KDT *kd)
{
  double extents [6], a [8][3], b [12], p [3], q [3], deps;
  SET *leaves = NULL, *item;
  KDT *leaf;
  TRI *t;
  int i;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;

  load_nodes (msh->cur_nodes, ele->type, ele->nodes, a);

  for (i = 0; i < ele->type; i ++)
  {
    if (a [i][0] < extents [0]) extents [0] = a [i][0];
    if (a [i][1] < extents [1]) extents [1] = a [i][1];
    if (a [i][2] < extents [2]) extents [2] = a [i][2];
    if (a [i][0] > extents [3]) extents [3] = a [i][0];
    if (a [i][1] > extents [4]) extents [4] = a [i][1];
    if (a [i][2] > extents [5]) extents [5] = a [i][2];
  }

  extents [0] -= GEOMETRIC_EPSILON;
  extents [1] -= GEOMETRIC_EPSILON;
  extents [2] -= GEOMETRIC_EPSILON;
  extents [3] += GEOMETRIC_EPSILON;
  extents [4] += GEOMETRIC_EPSILON;
  extents [5] += GEOMETRIC_EPSILON;

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
      if (gjk ((double*)a, ele->type, b, 4, p, q) < GEOMETRIC_EPSILON)
      {
	SET_Free (NULL, &leaves);
	return 1;
      }
    }
  }

  SET_Free (NULL, &leaves);
  return 0;
}

/* trim faces adjacent to the input triangulation */
static TRI* trimadj (MESH *msh, KDT *kdtri, double *point, double *normal, int *m)
{
  double pla [4], val [2], *ver [5], *av [6], *bv [6], *pnt;
  double *v [4], (*cur) [3] = msh->cur_nodes;
  MEM setmem, pntmem, trimem;
  SET *verts, *tset, *item;
  KDT *kdt2, *kd;
  TRI *t, *s, *tri;
  int nav, nbv, i;
  FACE *fac;

  MEM_Init (&pntmem, sizeof (double [3]), 128);
  MEM_Init (&trimem, sizeof (TRI), 128);
  MEM_Init (&setmem, sizeof (SET), 128);
  verts = tset = NULL;

  COPY (normal, pla);
  NORMALIZE (pla);
  pla [3] = - DOT (pla, point);

  for (fac = msh->faces; fac; fac = fac->n)
  {
    for (i = 0; i < fac->type; i ++)
    {
      v [i] = cur [fac->nodes [i]];
    }

    if (face_touches_triangulation (v, fac->type, fac->normal, kdtri))
    {
      ver [0] = v [fac->type-1];
      for (i = 0; i < fac->type; i ++) ver [i+1] = v [i];

      val [0] = PLANE (pla, ver [0]);
      for (nav = nbv = 0, i = 1; i <= fac->type; i ++)
      {
	val [1] = PLANE (pla, ver [i]);

	if (val [0] * val [1] < - GEOMETRIC_EPSILON * GEOMETRIC_EPSILON)
	{
	  ERRMEM (pnt = MEM_Alloc (&pntmem));
	  PLANESEG (pla, ver [i-1], ver [i], pnt);
	  bv [nbv ++] = pnt;
	  av [nav ++] = pnt;
	  SET_Insert (&setmem, &verts, pnt, NULL);
	}

	if (val [1] <= GEOMETRIC_EPSILON)
	{
	  bv [nbv ++] = ver [i];
	}

	if (val [1] >= -GEOMETRIC_EPSILON)
	{
	  av [nav ++] = ver [i];
	}

	val [0] = val [1];
      }

      ASSERT_DEBUG (nbv <= 6 && nav <= 6, "Too many vertices of a trimmed face");

      for (i = 2; i < nbv; i ++)
      {
	ERRMEM (t = MEM_Alloc (&trimem));
	t->ver [0] = bv [0];
	t->ver [1] = bv [i-1];
	t->ver [2] = bv [i];
	t->flg = fac->surface;
	COPY (fac->normal, t->out);
	SET_Insert (&trimem, &tset, t, NULL);
      }

      for (i = 2; i < nav; i ++)
      {
	ERRMEM (t = MEM_Alloc (&trimem));
	t->ver [0] = av [0];
	t->ver [1] = av [i-1];
	t->ver [2] = av [i];
	t->flg = fac->surface;
	COPY (fac->normal, t->out);
	SET_Insert (&trimem, &tset, t, NULL);
      }
    }
    else
    {
      for (i = 2; i < fac->type; i ++)
      {
	ERRMEM (t = MEM_Alloc (&trimem));
	t->ver [0] = v [0];
	t->ver [1] = v [i-1];
	t->ver [2] = v [i];
	t->flg = fac->surface;
	COPY (fac->normal, t->out);
	SET_Insert (&trimem, &tset, t, NULL);
      }
    }
  }

  i = SET_Size (verts);
  ERRMEM (pnt = malloc (i * sizeof (double [3])));
  for (item = SET_First (verts), ver [1] = pnt; item; item = SET_Next (item), ver [1] += 3)
  {
    ver [0] = item->data;
    COPY (ver [0], ver [1]);
  }
  kdt2 = KDT_Create (i, pnt, GEOMETRIC_EPSILON);
  free (pnt);

  i = KDT_Size (kdt2); /* index kd-tree nodes => KDT->n */
  *m = SET_Size (tset);
  ERRMEM (tri = MEM_CALLOC (i * sizeof (double [3]) + (*m) * sizeof (TRI)));
  pnt = (double*) (tri + *m);

  for (kd = KDT_First (kdt2); kd; kd = KDT_Next (kd))
  {
    ver [0] = kd->p;
    ver [1] = &pnt [3*kd->n];
    COPY (ver[0], ver[1]);
  }

  for (item = SET_First (tset), t = tri; item; item = SET_Next (item), t ++)
  {
    s = item->data;
    t->flg = s->flg;
    COPY (s->out, t->out);
    for (i = 0; i < 3; i ++)
    {
      if (SET_Contains (verts, s->ver [i], NULL))
      {
	kd = KDT_Nearest (kdt2, s->ver [i], GEOMETRIC_EPSILON);
	ASSERT_DEBUG (kd, "Kd-tree point query failed");
	t->ver [i] = &pnt [3*kd->n];
      }
      else t->ver [i] = s->ver [i];
    }
  }

  KDT_Destroy (kdt2);
  MEM_Release (&pntmem);
  MEM_Release (&setmem);
  MEM_Release (&trimem);

  return tri;
}

/* trim the surface by the plane; output triangle vertices points either to mesh vertices or is allocated after the triangles */
static void trim (MESH *msh, double *point, double *normal, TRI **below, int *mbelow, TRI **above, int *mabove)
{
  double pla [4], val [2], *ver [5], *av [6], *bv [6], *pnt;
  double (*cur) [3] = msh->cur_nodes;
  SET *verts, *at, *bt, *item;
  MEM setmem, pntmem, trimem;
  KDT *kdtree, *kd;
  int nav, nbv, i;
  TRI *t, *s;
  FACE *fac;

  MEM_Init (&pntmem, sizeof (double [3]), 128);
  MEM_Init (&trimem, sizeof (TRI), 128);
  MEM_Init (&setmem, sizeof (SET), 128);
  verts = at = bt = NULL;

  COPY (normal, pla);
  NORMALIZE (pla);
  pla [3] = - DOT (pla, point);

  for (fac = msh->faces; fac; fac = fac->n)
  {
    ver [0] = cur [fac->nodes [fac->type-1]];
    for (i = 0; i < fac->type; i ++) ver [i+1] = cur [fac->nodes [i]];

    val [0] = PLANE (pla, ver [0]);
    for (nav = nbv = 0, i = 1; i <= fac->type; i ++)
    {
      val [1] = PLANE (pla, ver [i]);

      if (val [0] * val [1] < - GEOMETRIC_EPSILON * GEOMETRIC_EPSILON)
      {
	ERRMEM (pnt = MEM_Alloc (&pntmem));
	PLANESEG (pla, ver [i-1], ver [i], pnt);
	bv [nbv ++] = pnt;
	av [nav ++] = pnt;
	SET_Insert (&setmem, &verts, pnt, NULL);
      }

      if (val [1] <= GEOMETRIC_EPSILON)
      {
	bv [nbv ++] = ver [i];
      }

      if (val [1] >= -GEOMETRIC_EPSILON)
      {
	av [nav ++] = ver [i];
      }

      val [0] = val [1];
    }

    ASSERT_DEBUG (nbv <= 6 && nav <= 6, "Too many vertices of a trimmed face");

    for (i = 2; i < nbv; i ++)
    {
      ERRMEM (t = MEM_Alloc (&trimem));
      t->ver [0] = bv [0];
      t->ver [1] = bv [i-1];
      t->ver [2] = bv [i];
      t->flg = fac->surface;
      COPY (fac->normal, t->out);
      SET_Insert (&trimem, &bt, t, NULL);
    }

    for (i = 2; i < nav; i ++)
    {
      ERRMEM (t = MEM_Alloc (&trimem));
      t->ver [0] = av [0];
      t->ver [1] = av [i-1];
      t->ver [2] = av [i];
      t->flg = fac->surface;
      COPY (fac->normal, t->out);
      SET_Insert (&trimem, &at, t, NULL);
    }
  }

  i = SET_Size (verts);
  ERRMEM (pnt = malloc (i * sizeof (double [3])));
  for (item = SET_First (verts), ver [1] = pnt; item; item = SET_Next (item), ver [1] += 3)
  {
    ver [0] = item->data;
    COPY (ver [0], ver [1]);
  }
  kdtree = KDT_Create (i, pnt, GEOMETRIC_EPSILON);
  free (pnt);

  i = KDT_Size (kdtree); /* index kd-tree nodes => KDT->n */
  *mbelow = SET_Size (bt);
  *mabove = SET_Size (at);
  ERRMEM (*below = MEM_CALLOC (i * sizeof (double [3]) + (*mbelow) * sizeof (TRI)));
  ERRMEM (*above = MEM_CALLOC ((*mabove) * sizeof (TRI)));
  pnt = (double*) ((*below) + *mbelow);

  for (kd = KDT_First (kdtree); kd; kd = KDT_Next (kd))
  {
    ver [0] = kd->p;
    ver [1] = &pnt [3*kd->n];
    COPY (ver[0], ver[1]);
  }

  for (item = SET_First (bt), t = *below; item; item = SET_Next (item), t ++)
  {
    s = item->data;
    t->flg = s->flg;
    COPY (s->out, t->out);
    for (i = 0; i < 3; i ++)
    {
      if (SET_Contains (verts, s->ver [i], NULL))
      {
	kd = KDT_Nearest (kdtree, s->ver [i], GEOMETRIC_EPSILON);
	ASSERT_DEBUG (kd, "Kd-tree point query failed");
	t->ver [i] = &pnt [3*kd->n];
      }
      else t->ver [i] = s->ver [i];
    }
  }

  for (item = SET_First (at), t = *above; item; item = SET_Next (item), t ++)
  {
    s = item->data;
    t->flg = s->flg;
    COPY (s->out, t->out);
    for (i = 0; i < 3; i ++)
    {
      if (SET_Contains (verts, s->ver [i], NULL))
      {
	kd = KDT_Nearest (kdtree, s->ver [i], GEOMETRIC_EPSILON);
	ASSERT_DEBUG (kd, "Kd-tree point query failed");
	t->ver [i] = &pnt [3*kd->n];
      }
      else t->ver [i] = s->ver [i];
    }
  }

  MEM_Release (&setmem);
  MEM_Release (&trimem);
  MEM_Release (&pntmem);
  KDT_Destroy (kdtree);
}

/* produce a split mesh out of element set and cut faces set; 'newels' and 'newnod' are used in case of partial topological split;
 * the 'newels' contain new elements adjacent to the splitting surface, 'newnod' maps nodes which need to doubled in these elements */
static MESH* produce_split_mesh (MESH *msh, SET *els, int *cutfaces, int surfid, SET *newels, MAP *newnod)
{
  double (*nodes) [3], (*curno) [3] = msh->cur_nodes;
  int *elements, *surfaces, *ptr, i, j;
  MAP *nodmap, *jtem;
  ELEMENT *ele;
  MEM mapmem;
  FACE *fac;
  SET *item;
  MESH *out;

  MEM_Init (&mapmem, sizeof (MAP), MEMCHUNK);

  for (nodmap = NULL, item = SET_First (els), j = 0; item; item = SET_Next (item))
  {
    ele = item->data;
    for (i = 0; i < ele->type; i ++)
    {
      if (!MAP_Find_Node (nodmap, (void*) (long) ele->nodes [i], NULL))
      {
        ASSERT_DEBUG_EXT (jtem = MAP_Insert (&mapmem, &nodmap, /* map subset of used nodes */
          (void*) (long) ele->nodes [i], (void*) (long) j ++, NULL),
	  "Map insertion failed");
      }
    }
  }
  for (item = SET_First (newels); item; item = SET_Next (item)) /* add unmapped new element nodes */
  {
    ele = item->data;
    for (i = 0; i < ele->type; i ++)
    {
      if (!MAP_Find_Node (newnod, (void*) (long) ele->nodes [i], NULL) &&
	  !MAP_Find_Node (nodmap, (void*) (long) ele->nodes [i], NULL))
      {
        ASSERT_DEBUG_EXT (jtem = MAP_Insert (&mapmem, &nodmap, /* map subset of used nodes */
          (void*) (long) ele->nodes [i], (void*) (long) j ++, NULL),
	  "Map insertion failed");
      }
    }
  }

  j += MAP_Size (newnod); /* increase the node set size by |newnod| */

  ERRMEM (nodes = malloc (j * sizeof (double [3])));

  for (jtem = MAP_First (nodmap); jtem; jtem = MAP_Next (jtem)) /* copy a subset of old nodes onto new nodes */
  {
    i = (int) (long) jtem->key; /* old index */
    j = (int) (long) jtem->data; /* new index */
    double *a = curno [i],
	   *b = nodes [j];
    COPY (a, b); /* copy coordinates */
  }

  for (jtem = MAP_First (newnod); jtem; jtem = MAP_Next (jtem)) /* now for 'newnod' set => copy a subset of old nodes onto new nodes */
  {
    i = (int) (long) jtem->key; /* old index */
    j = (int) (long) jtem->data; /* new index */
    double *a = curno [i],
	   *b = nodes [j];
    COPY (a, b); /* copy coordinates */
  }

  j = SET_Size (els); /* regular elements */
  j += SET_Size (newels); /* elements one-sidedly adjacent to the splitting surface */
  ERRMEM (elements = malloc ((j + 1) * sizeof (int [10]))); /* overestimate */

  for (ptr = elements, item = SET_First (els); item; item = SET_Next (item)) /* for each regular element */
  {
    ele = item->data;
    *ptr = ele->type; ptr ++;
    for (i = 0; i < ele->type; i ++, ptr ++)
    {
      ASSERT_DEBUG_EXT (jtem = MAP_Find_Node (nodmap, (void*) (long) ele->nodes [i], NULL), "Node mapping failed");
      *ptr = (int) (long) jtem->data;
    }
    *ptr = ele->volume; ptr ++;
  }
  for (item = SET_First (newels); item; item = SET_Next (item)) /* for each element adjacent to the splitting */
  {
    ele = item->data;
    *ptr = ele->type; ptr ++;
    for (i = 0; i < ele->type; i ++, ptr ++)
    {
      if (!(jtem = MAP_Find_Node (newnod, (void*) (long) ele->nodes [i], NULL))) /* if externally given new node mapping was not found */
      {
        ASSERT_DEBUG_EXT (jtem = MAP_Find_Node (nodmap, (void*) (long) ele->nodes [i], NULL), "Node mapping failed"); /* use internal node mapping */
      }
      *ptr = (int) (long) jtem->data;
    }
    *ptr = ele->volume; ptr ++;
  }
  *ptr = 0; /* mark end */

  ERRMEM (surfaces = malloc (sizeof (int [6]) * 8 * (j + 1))); /* overestimate */

  surfaces [0] = surfid; /* global id will not be used (all faces are mapped) */

  for (ptr = surfaces + 1, item = SET_First (els); item; item = SET_Next (item)) /* map existing faces */
  {
    ele = item->data;
    for (fac = ele->faces; fac; fac = fac->next)
    {
      *ptr = fac->type; ptr ++;
      for (i = 0; i < fac->type; i ++, ptr ++)
      {
	ASSERT_DEBUG_EXT (jtem = MAP_Find_Node (nodmap, (void*) (long) fac->nodes [i], NULL), "Node mapping failed");
	*ptr = (int) (long) jtem->data;
      }
      *ptr = fac->surface; ptr ++;
    }
  }
  for (item = SET_First (newels); item; item = SET_Next (item)) /* map existing faces */
  {
    ele = item->data;
    for (fac = ele->faces; fac; fac = fac->next)
    {
      *ptr = fac->type; ptr ++;
      for (i = 0; i < fac->type; i ++, ptr ++)
      {
	if (!(jtem = MAP_Find_Node (newnod, (void*) (long) fac->nodes [i], NULL))) /* if externally given new node mapping was not found */
	{
	  ASSERT_DEBUG_EXT (jtem = MAP_Find_Node (nodmap, (void*) (long) fac->nodes [i], NULL), "Node mapping failed"); /* use internal mapping */
	}
	*ptr = (int) (long) jtem->data;
      }
      *ptr = fac->surface; ptr ++;
    }
  }
  for (; *cutfaces; cutfaces += ABS (cutfaces [0]) + 1) /* map newly created faces */
  {
    *ptr = ABS (cutfaces [0]); ptr ++;
    for (i = 1; i <= ABS (cutfaces [0]); i ++, ptr ++)
    {
      if (cutfaces [0] < 0) *ptr = cutfaces [i]; /* newly split nodes */
      else
      {
	ASSERT_DEBUG_EXT (jtem = MAP_Find_Node (nodmap, (void*) (long) cutfaces [i], NULL), "Node mapping failed"); /* locally remapped old nodes */
        *ptr = (int) (long) jtem->data;
      }
    }
    *ptr = surfid; ptr ++;
  }
  *ptr = 0; /* mark end */

  out = MESH_Create (nodes, elements, surfaces); /* this sorts out current nodes */

  /* copy reference nodes now */

  double  (*ref0) [3] = msh->ref_nodes,
	  (*ref1) [3] = out->ref_nodes;

  for (jtem = MAP_First (nodmap); jtem; jtem = MAP_Next (jtem)) /* copy a subset of old reference nodes onto new nodes */
  {
    i = (int) (long) jtem->key; /* old index */
    j = (int) (long) jtem->data; /* new index */
    double *a = ref0 [i],
	   *b = ref1 [j];
    COPY (a, b);
  }

  for (jtem = MAP_First (newnod); jtem; jtem = MAP_Next (jtem)) /* now for 'newnod' set => copy a subset of old reference nodes onto new nodes */
  {
    i = (int) (long) jtem->key; /* old index */
    j = (int) (long) jtem->data; /* new index */
    double *a = ref0 [i],
	   *b = ref1 [j];
    COPY (a, b);
  }

  /* TODO: take care of material mapping, if element materials were presecribed */

  free (nodes);
  free (elements);
  free (surfaces);
  MEM_Release (&mapmem);

  return out;
}

/* try to split mesh globally using only inter-element boundaries */
static int inter_element_global_split (MESH *msh, double *point, double *normal, int surfid, MESH **one, MESH **two)
{
  int code, bulk, i, j, onpla [4], *cutfaces, *cut;
  double (*nod) [3] = msh->cur_nodes, nn [3];
  ELEMENT *ele;
  MEM setmem;
  SET *below,
      *above;

  MEM_Init (&setmem, sizeof (SET), MEMCHUNK);
  ERRMEM (cutfaces = malloc (sizeof (int [5]) * 8 * (msh->surfeles_count + msh->bulkeles_count))); /* overestimate */

  COPY (normal, nn);
  NORMALIZE (nn);

  below = above = NULL;
  cut = cutfaces;

  for (bulk = 0, ele = msh->surfeles; ele;)
  {
    for (code = i = j = 0; i < ele->type; i ++)
    {
      double a [3], dot;

      SUB (nod [ele->nodes [i]], point, a);
      dot = DOT (a, nn);
      if (dot < -GEOMETRIC_EPSILON)
      {
	if (!code) code = -1;
	else if (code > 0) { code = 0; break; }
      }
      else if (dot > GEOMETRIC_EPSILON)
      {
	if (!code) code = 1;
	else if (code < 0) { code = 0; break; }
      }
      else /* on plane */
      {
	onpla [j ++] = ele->nodes [i];
	ASSERT_DEBUG (j <= 4, "Inconsitent on plane nodes count while splitting mesh");
      }
    }

    if (j >= 3) /* at least a trinagle */
    {
      cut [0] = j; cut ++;
      for (i = 0; i < j; i ++) cut [i] = onpla [i]; /* output cut face */
      cut += j;
    }

    if (code < 0) SET_Insert (&setmem, &below, ele, NULL);
    else if (code > 0) SET_Insert (&setmem, &above, ele, NULL);
    else /* element crossing the plane */
    {
      free (cutfaces);
      MEM_Release (&setmem);
      return 0; /* report failure */
    }

    if (!bulk && !ele->next) bulk = 1, ele = msh->bulkeles;
    else ele = ele->next;
  }
  cut [0] = 0; /* mark end */

  /* here we are => inter-element cut succeeded;
   * now we re-map nodes and faces and create two meshes */

  *one = produce_split_mesh (msh, below, cutfaces, surfid, NULL, NULL);
  *two = produce_split_mesh (msh, above, cutfaces, surfid, NULL, NULL);

  free (cutfaces);
  MEM_Release (&setmem);

  return 1;
}

/* try to split mesh locally using only inter-element boundaries and topological adjacency */
static int inter_element_local_split (MESH *msh, KDT *kdtri, double *point, double *normal, int surfid, MESH **one)
{
  int code, bulk, i, j, k, flg, onpla [4], *cutfaces, *cut;
  double (*nod) [3] = msh->cur_nodes, nn [3];
  MEM setmem, mapmem;
  SET *below, *above;
  ELEMENT *ele;
  MAP *newnod;

  MEM_Init (&setmem, sizeof (SET), MEMCHUNK);
  MEM_Init (&mapmem, sizeof (MAP), MEMCHUNK);
  ERRMEM (cutfaces = malloc (sizeof (int [5]) * 8 * (msh->surfeles_count + msh->bulkeles_count))); /* overestimate */

  COPY (normal, nn);
  NORMALIZE (nn);

  k = msh->nodes_count;
  below = above = NULL;
  newnod = NULL;
  cut = cutfaces;

  for (bulk = 0, ele = msh->surfeles; ele;)
  {
    for (code = i = j = 0; i < ele->type; i ++)
    {
      double a [3], dot;

      SUB (nod [ele->nodes [i]], point, a);
      dot = DOT (a, nn);
      if (dot < -GEOMETRIC_EPSILON)
      {
	if (!code) code = -1;
	else if (code > 0) { code = 0; break; }
      }
      else if (dot > GEOMETRIC_EPSILON)
      {
	if (!code) code = 1;
	else if (code < 0) { code = 0; break; }
      }
      else /* on plane */
      {
	onpla [j ++] = ele->nodes [i];
	ASSERT_DEBUG (j <= 4, "Inconsitent on plane nodes count while splitting mesh");
      }
    }

    flg = element_touches_triangulation (msh, ele, kdtri);

    if (flg && code == 0) /* element crossing the plane */
    {
      free (cutfaces);
      MEM_Release (&setmem);
      MEM_Release (&mapmem);
      return 0; /* report failure */
    }

    if (flg)
    {
      if (code < 0)
      {
	if (j >= 3) /* at least a trinagle */
	{
	  cut [0] = -j; cut ++; /* negative facet flag => use 'newnod' map in 'produce_split_mesh' */
	  for (i = 0; i < j; i ++)
	  {
	    MAP *item = MAP_Find_Node (newnod, (void*) (long) onpla [i], NULL);
	    if (item) cut [i] = (int) (long) item->data; /* output remapped cut face */
	    else
	    {
	      MAP_Insert (&mapmem, &newnod, (void*) (long) onpla [i], (void*) (long) k, NULL);
	      cut [i] = k ++;
	    }
	  }
	  cut += j;
	}

	SET_Insert (&setmem, &below, ele, NULL); /* below adjacent elements */
      }
      else /* code > 0 */
      {
	if (j >= 3) /* at least a trinagle */
	{
	  cut [0] = j; cut ++; /* positive facet flag => use old mesh nodes in 'produce_split_mesh' */
	  for (i = 0; i < j; i ++) cut [i] = onpla [i]; /* output cut face */
	  cut += j;
	}

	SET_Insert (&setmem, &above, ele, NULL); /* other elements */
      }
    }
    else SET_Insert (&setmem, &above, ele, NULL); /* other elements */

    if (!bulk && !ele->next) bulk = 1, ele = msh->bulkeles;
    else ele = ele->next;
  }
  cut [0] = 0; /* mark end */

  /* here we are => inter-element cut succeeded;
   * now we re-map nodes and faces and create new mesh */

  *one = produce_split_mesh (msh, above, cutfaces, surfid, below, newnod);

  free (cutfaces);
  MEM_Release (&setmem);
  MEM_Release (&mapmem);

  return 1;
}
#else

/* produce mesh from a subset of elements */
static MESH* produce_subset_mesh (MESH *msh, SET *els, int surfid)
{
  double (*nodes) [3], (*curno) [3] = msh->cur_nodes;
  int *elements, *surfaces, *ptr, i, j;
  MAP *nodmap, *jtem;
  ELEMENT *ele;
  MEM mapmem;
  FACE *fac;
  SET *item;
  MESH *out;

  MEM_Init (&mapmem, sizeof (MAP), MEMCHUNK);

  for (nodmap = NULL, item = SET_First (els), j = 0; item; item = SET_Next (item))
  {
    ele = item->data;
    for (i = 0; i < ele->type; i ++)
    {
      if (!MAP_Find_Node (nodmap, (void*) (long) ele->nodes [i], NULL))
      {
        ASSERT_DEBUG_EXT (jtem = MAP_Insert (&mapmem, &nodmap, /* map subset of used nodes */
          (void*) (long) ele->nodes [i], (void*) (long) j ++, NULL),
	  "Map insertion failed");
      }
    }
  }

  ERRMEM (nodes = malloc (j * sizeof (double [3])));

  for (jtem = MAP_First (nodmap); jtem; jtem = MAP_Next (jtem)) /* copy a subset of old nodes onto new nodes */
  {
    i = (int) (long) jtem->key; /* old index */
    j = (int) (long) jtem->data; /* new index */
    double *a = curno [i],
	   *b = nodes [j];
    COPY (a, b); /* copy coordinates */
  }

  j = SET_Size (els); /* regular elements */
  ERRMEM (elements = malloc ((j + 1) * sizeof (int [10]))); /* overestimate */

  for (ptr = elements, item = SET_First (els); item; item = SET_Next (item)) /* for each regular element */
  {
    ele = item->data;
    *ptr = ele->type; ptr ++;
    for (i = 0; i < ele->type; i ++, ptr ++)
    {
      ASSERT_DEBUG_EXT (jtem = MAP_Find_Node (nodmap, (void*) (long) ele->nodes [i], NULL), "Node mapping failed");
      *ptr = (int) (long) jtem->data;
    }
    *ptr = ele->volume; ptr ++;
  }
  *ptr = 0; /* mark end */

  ERRMEM (surfaces = malloc (sizeof (int [6]) * 8 * (j + 1))); /* overestimate */

  surfaces [0] = surfid; /* global id used for not mapped faces */

  for (ptr = surfaces + 1, item = SET_First (els); item; item = SET_Next (item)) /* map existing faces */
  {
    ele = item->data;
    for (fac = ele->faces; fac; fac = fac->next)
    {
      *ptr = fac->type; ptr ++;
      for (i = 0; i < fac->type; i ++, ptr ++)
      {
	ASSERT_DEBUG_EXT (jtem = MAP_Find_Node (nodmap, (void*) (long) fac->nodes [i], NULL), "Node mapping failed");
	*ptr = (int) (long) jtem->data;
      }
      *ptr = fac->surface; ptr ++;
    }
  }
  *ptr = 0; /* mark end */

  out = MESH_Create (nodes, elements, surfaces); /* this sorts out current nodes */

  /* copy reference nodes now */

  double  (*ref0) [3] = msh->ref_nodes,
	  (*ref1) [3] = out->ref_nodes;

  for (jtem = MAP_First (nodmap); jtem; jtem = MAP_Next (jtem)) /* copy a subset of old reference nodes onto new nodes */
  {
    i = (int) (long) jtem->key; /* old index */
    j = (int) (long) jtem->data; /* new index */
    double *a = ref0 [i],
	   *b = ref1 [j];
    COPY (a, b);
  }

  /* TODO: take care of material mapping, if element materials were presecribed */

  free (nodes);
  free (elements);
  free (surfaces);
  MEM_Release (&mapmem);

  return out;
}
#endif

/* recursive neighbour marking */
static void mark_neighs (ELEMENT *ele, short flag)
{
  if (ele->flag != flag)
  {
    ele->flag = flag;

    for (int i = 0; i < ele->neighs; i ++)
    {
      mark_neighs (ele->adj [i], flag);
    }
  }
}

/* create mesh from vector of nodes, element list in format =>
 * {nuber of nodes, node0, node1, ...}, {REPEAT}, ..., 0 (end of list); and surface kinds in format =>
 * global surface, {number of nodes, node0, node1, ..., surface}, {REPEAT}, ..., 0 (end of list); */
MESH* MESH_Create (double (*nodes) [3], int *elements, int *surfaces)
{
  int maximal_node,
      minimal_node,
      elements_count,
      faces_count,
      temp, *eleptr, n;
  double (*ref) [3],
	 (*cur) [3];
  MEM *elemem,
      facmem,
      mapmem;
  ELEMENT *ele, *enx, *elist;
  FACE *fac, *cac, *gac, *flist;
  MAP *faces, *smap;
  MESH *msh;
  
  maximal_node = 0;
  minimal_node = INT_MAX;
  elements_count = 0;
  faces_count = 0;

  /* create mesh storage */
  ERRMEM (msh = MEM_CALLOC (sizeof (MESH)));
  elemem = &msh->elemem;
 
  /* calculate elements */ 
  for (eleptr = elements; eleptr [0]; eleptr += (eleptr [0]+2)) elements_count ++;

  MEM_Init (elemem, sizeof (ELEMENT), elements_count);
  MEM_Init (&facmem, sizeof (FACE), MEMCHUNK);
  MEM_Init (&mapmem, sizeof (MAP), MEMCHUNK);
  MEM_Init (&msh->mapmem, sizeof (MAP), MIN (elements_count, MEMCHUNK));
  msh->map = NULL;

  elist = NULL;
  flist = NULL;
  faces = NULL;

  /* create elements list & face adjacency map */
  for (eleptr = elements; eleptr [0]; eleptr += (eleptr [0]+2))
  {
    ASSERT (
      eleptr [0] == 4 || /* tetrahedron */
      eleptr [0] == 5 || /* pyramid */
      eleptr [0] == 6 || /* wedge */
      eleptr [0] == 8,   /* hexahedron */
      ERR_MSH_UNSUPPORTED_ELEMENT);

    ele = create_element (elemem, eleptr);
    flist = create_faces (&facmem, &mapmem, &faces, ele, flist);
    ele->next = elist;
    elist = ele;

    /* node number extrema */
    temp = maximal (eleptr);
    if (temp > maximal_node)
      maximal_node = temp;
    temp = minimal (eleptr);
    if (temp < minimal_node)
      minimal_node = temp;
  }

  /* calculate faces */
  for (fac = flist; fac; fac = fac->next)
    if (fac->ele) faces_count ++;

  /* alocate additional storage */
  MEM_Init (&msh->facmem, sizeof (FACE), faces_count);
  msh->nodes_count = (maximal_node - minimal_node + 1);
  ERRMEM (msh->ref_nodes = malloc (sizeof (double [3]) * (msh->nodes_count * 2)));
  msh->cur_nodes = msh->ref_nodes + msh->nodes_count;
  msh->surfeles_count = msh->bulkeles_count = 0;
  msh->surfeles = msh->bulkeles = NULL;

  /* set up elements */
  for (ele = elist; ele; ele = enx)
  {
    enx = ele->next;

    if (minimal_node > 0) /* impose 0-based indexing */
    {
      for (temp = 0; temp < ele->type; temp ++)
	ele->nodes [temp] -= minimal_node;
    }

    ele->prev = NULL;
   
    if (ele->neighs < neighs (ele->type)) /* surface element */
    {
      msh->surfeles_count ++;
      ele->next = msh->surfeles;
      if (msh->surfeles) msh->surfeles->prev = ele;
      msh->surfeles = ele;
    }
    else /* bulk element */
    {
      msh->bulkeles_count ++;
      ele->next = msh->bulkeles;
      if (msh->bulkeles) msh->bulkeles->prev = ele;
      msh->bulkeles = ele;
    }
  }

  /* create surfaces map => skip first element of 'surfaces' == the global surface kind */
  for (eleptr = (surfaces + 1), smap = NULL, temp = 0;
    eleptr [0]; eleptr += (eleptr [0]+2), temp ++)
  {
    fac = MEM_Alloc (&facmem);
    
    ASSERT (
      eleptr [0] == 3 || /* triangle */
      eleptr [0] == 4,   /* quad */
      ERR_MSH_UNSUPPORTED_FACE);

    fac->type = eleptr [0];
    for (n = 0; n < eleptr [0]; n ++)
      fac->nodes [n] = eleptr [n+1];
    sort (fac->nodes, fac->nodes+fac->type-1);

    fac->surface = eleptr [eleptr [0] + 1];
    MAP_Insert (&mapmem, &smap, fac, /* map by the type/nodes key */
      fac, (MAP_Compare) face_compare);
  }

  /* set up nodes */
  for (temp = minimal_node,
       ref = msh->ref_nodes,
       cur = msh->cur_nodes;
       temp <= maximal_node;
       temp ++, ref ++, cur ++)
  {
    COPY (nodes [temp], *ref);
    COPY (nodes [temp], *cur);
  }

  /* set up faces */
  for (fac = flist; fac; fac = fac->next)
  {
    if (fac->ele) /* see (***) */
    {
      ele = fac->ele;

      cac = MEM_Alloc (&msh->facmem);
      setup_face (ele, fac->index, cac, 0); /* setup face nodes without sorting them */
      cac->index = fac->index;
      cac->ele = fac->ele;
      setup_normal (msh->cur_nodes, cac); /* calculate outer spatial normal */
      cac->next = ele->faces; /* append element face list */
      ele->faces = cac;

      /* set the mapped surface kind if possible => otherwise the global one */
      gac = MAP_Find (smap, fac, (MAP_Compare) face_compare); 
      cac->surface = (gac ? gac->surface : surfaces [0]);
    }
  }

  /* create mesh face list */
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    for (fac = ele->faces; fac; fac = fac->next)
    {
      fac->n = msh->faces;
      msh->faces = fac;
    }
  }

  /* clean up */
  MEM_Release (&facmem);
  MEM_Release (&mapmem);

  return msh;
}

/* create a meshed hexahedron by specifying its eight nodes and
 * division numbers along three edges adjacent to the 1st node */
MESH* MESH_Hex (double (*nodes) [3], int i, int j, int k, int *surfaces, int volume, double *dx, double *dy, double *dz)
{
  int *ele, *ee, *ss, *sur, n, m, mx, my, mz, ii, jj, kk, nn;
  double (*nod) [3], x, y, z, ddx, ddy, ddz;
  double vx [3], vy [3], vz [3], vxvy [3];
  MESH *msh;


  SUB (nodes [1], nodes [0], vx);
  SUB (nodes [2], nodes [1], vy);
  SUB (nodes [4], nodes [0], vz);
  PRODUCT (vx, vy, vxvy);
  if (DOT (vxvy, vz) < 0) /* change orientation: swap upper and lowe nodes */
  {
    double copy [4][3];

    for (ii = 0; ii < 4; ii ++) { COPY (nodes [ii], copy [ii]); }
    for (ii = 0; ii < 4; ii ++) { COPY (nodes [ii+4], nodes [ii]); }
    for (ii = 0; ii < 4; ii ++) { COPY (copy [ii], nodes [ii+4]); }
  }

  ERRMEM (nod = malloc (sizeof (double [3])*(i+1)*(j+1)*(k+1)));
  ERRMEM (ele = malloc (sizeof (int [10])*i*j*k + sizeof (int)));
  ERRMEM (sur = malloc (sizeof(int [2]) + sizeof (int [6])*((2*i*j)+(2*j*k)+(2*i*k))));

  /* create the unit cube mesh */

  ddx =  2.0 / ((double)i);
  ddy =  2.0 / ((double)j);
  ddz =  2.0 / ((double)k);
  mx = my = mz = 0;

  if (dx == NULL)
  {
    ERRMEM (dx = malloc (i*sizeof (double)));
    for (mx = 0; mx < i; mx ++) dx[mx] = ddx;
  }
  else
  {
    for (x = 0, ii = 0; ii < i; ii ++)
    {
      ASSERT_DEBUG (dx [ii] > 0.0, "There must hold => dx[i] > 0 for all i");
      x += dx [ii];
    }
    for (ii = 0; ii < i; ii ++) dx [ii] *= 2.0 / x;
  }

  if (dy == NULL)
  {
    ERRMEM (dy = malloc (j*sizeof (double)));
    for (my = 0; my < j; my ++) dy[my] = ddy;
  }
  else
  {
    for (y = 0, jj = 0; jj < j; jj ++)
    {
      ASSERT_DEBUG (dy [jj] > 0.0, "There must hold => dy[j] > 0 for all j");
      y += dy [jj];
    }
    for (jj = 0; jj < j; jj ++) dy [jj] *= 2.0 / y;
  }

  if (dz == NULL)
  {
    ERRMEM (dz = malloc (k*sizeof (double)));
    for (mz = 0; mz < k; mz ++) dz[mz] = ddz;
  }
  else
  {
    for (z = 0, kk = 0; kk < k; kk ++)
    {
      ASSERT_DEBUG (dz [kk] > 0.0, "There must hold => dz[k] > 0 for all k");
      z += dz [kk];
    }
    for (kk = 0; kk < k; kk ++) dz [kk] *= 2.0 / z;
  }

  /* nodes */
  for (z = -1, kk = nn = 0; kk <= k; z += dz[MIN(kk,k-1)], kk ++)
  for (y = -1, jj = 0; jj <= j; y += dy[MIN(jj,j-1)], jj ++)
  for (x = -1, ii = 0; ii <= i; x += dx[MIN(ii,i-1)], ii ++, nn ++)
  {
    nod [nn][0] = nodes[0][0]*HEX0(x,y,z) + nodes[1][0]*HEX1(x,y,z) + nodes[2][0]*HEX2(x,y,z) + nodes[3][0]*HEX3(x,y,z)+
                  nodes[4][0]*HEX4(x,y,z) + nodes[5][0]*HEX5(x,y,z) + nodes[6][0]*HEX6(x,y,z) + nodes[7][0]*HEX7(x,y,z);
    nod [nn][1] = nodes[0][1]*HEX0(x,y,z) + nodes[1][1]*HEX1(x,y,z) + nodes[2][1]*HEX2(x,y,z) + nodes[3][1]*HEX3(x,y,z)+
                  nodes[4][1]*HEX4(x,y,z) + nodes[5][1]*HEX5(x,y,z) + nodes[6][1]*HEX6(x,y,z) + nodes[7][1]*HEX7(x,y,z);
    nod [nn][2] = nodes[0][2]*HEX0(x,y,z) + nodes[1][2]*HEX1(x,y,z) + nodes[2][2]*HEX2(x,y,z) + nodes[3][2]*HEX3(x,y,z)+
                  nodes[4][2]*HEX4(x,y,z) + nodes[5][2]*HEX5(x,y,z) + nodes[6][2]*HEX6(x,y,z) + nodes[7][2]*HEX7(x,y,z);
  }

  if (mx) free (dx);
  if (my) free (dy);
  if (mz) free (dz);

  /* elements */ 
  n = (i+1)*(j+1); 
  sur [0] = surfaces [0]; /* will not be used by let be set sensibly */
  for (kk = 0, ee = ele, ss = sur + 1; kk < k; kk ++)
  for (jj = 0; jj < j; jj ++)
  for (ii = 0; ii < i; ii ++, ee += 10)
  {
    nn = ((kk*n)+(jj*(i+1))+ii);

    ee [0] = 8;
    ee [1] = nn;
    ee [2] = nn+1;
    ee [3] = nn+i+2;
    ee [4] = nn+i+1;
    ee [5] = ee[1]+n;
    ee [6] = ee[2]+n;
    ee [7] = ee[3]+n;
    ee [8] = ee[4]+n;
    ee [9] = volume; /* volume identifier */

    if (kk == 0) /* face 1 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[0][m]];
      ss [5] = surfaces [0];
      ss += 6;
    }
    if (jj == 0) /* face 2 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[1][m]];
      ss [5] = surfaces [4];
      ss += 6;
    }
    if (ii == (i-1)) /* face 3 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[2][m]];
      ss [5] = surfaces [3];
      ss += 6;
    }
    if (jj == (j-1)) /* face 4 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[3][m]];
      ss [5] = surfaces [2];
      ss += 6;
    }
    if (ii == 0) /* face 5 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[4][m]];
      ss [5] = surfaces [1];
      ss += 6;
    }
    if (kk == (k-1)) /* face 6 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[5][m]];
      ss [5] = surfaces [5];
      ss += 6;
    }
  }
  ee [0] = 0;
  ss [0] = 0;

  /* create the mesh object */
  msh = MESH_Create (nod, ele, sur);

  /* clean */
  free (nod);
  free (ele);
  free (sur);

  return msh;
}

/* create pipe like mesh using a point, a direction, an inner radius,
 * a thickness and subdivison counts along the direction, radius and thickness */
MESH* MESH_Pipe (double *pnt, double *dir, double rin, double thi,
                 int ndir, int nrad, int nthi, int *surfaces, int volume) /* surfaces: bottom, top, inner, outer */
{
  double (*nodes) [3], (*nn) [3];
  int *elems, *surfs, *ee, *ss;

  int nodes_count = (ndir+1)*nrad*(nthi+1),
      elems_count = ndir*nrad*nthi;

  double ksi_d, ksi_r, ksi_t, d_ksi_d, d_ksi_r, d_ksi_t;
  double normal [3], a [3], b [3];
  int i, j, k;
  MESH *msh;

  COPY (dir, normal);
  NORMALIZE (normal);
  d_ksi_d = LEN (dir) / (double)ndir;
  d_ksi_r = 2.0*ALG_PI / (double)nrad;
  d_ksi_t = thi / (double)nthi;
  a [0] = 1; a [1] = a [2] = 0;
  PRODUCT (normal, a, b);
  if (LEN (b) < 1E-10)
  {
    a [1] = 1; a [0] = a [2] = 0;
    PRODUCT (normal, a, b);
    if (LEN (b) < 1E-10)
    {
      a [2] = 1; a [0] = a [1] = 0;
      PRODUCT (normal, a, b);
      if (LEN (b) < 1E-10)
      {
	ASSERT_DEBUG (0, "Failed to select another direction");
	return NULL;
      }
    }
  }
  PRODUCT (normal, b, a); /* a, b PERP normal */

  ERRMEM (nodes = malloc (nodes_count * sizeof (double [3])));

  for (nn = nodes, i = 0, ksi_d = 0.0; i < ndir + 1; i ++, ksi_d += d_ksi_d)
  {
    for (j = 0, ksi_r = 0; j < nrad; j ++, ksi_r += d_ksi_r)
    {
      for (k = 0, ksi_t = rin; k < nthi + 1; k ++, nn ++, ksi_t += d_ksi_t)
      {
	nn [0][0] = pnt [0] + normal [0] * ksi_d + a [0] * ksi_t * sin (ksi_r) + b [0] * ksi_t * cos (ksi_r);
	nn [0][1] = pnt [1] + normal [1] * ksi_d + a [1] * ksi_t * sin (ksi_r) + b [1] * ksi_t * cos (ksi_r);
	nn [0][2] = pnt [2] + normal [2] * ksi_d + a [2] * ksi_t * sin (ksi_r) + b [2] * ksi_t * cos (ksi_r);
      }
    }
  }

  ERRMEM (elems = malloc (elems_count * sizeof (int [10]) + sizeof (int)));
  ERRMEM (surfs = malloc (elems_count * sizeof (int [36]) + sizeof (int [2]))); /* too many */
  surfs [0] = surfaces [0];

  for (ee = elems, ss = surfs + 1, i = 0; i < ndir; i ++)
  {
    for (j = 0; j < nrad; j ++)
    {
      for (k = 0; k < nthi; k ++, ee += 10)
      {
	if (j < nrad - 1)
	{
	  ee [0] = 8;
	  ee [1] = i*nrad*(nthi+1) + j*(nthi+1) + k;
	  ee [2] = ee [1] + 1;
	  ee [3] = i*nrad*(nthi+1) + (j+1)*(nthi+1) + k+1;
	  ee [4] = ee [3] - 1;
	  ee [5] = ee [1] + nrad*(nthi+1);
	  ee [6] = ee [2] + nrad*(nthi+1);
	  ee [7] = ee [3] + nrad*(nthi+1);
	  ee [8] = ee [4] + nrad*(nthi+1);
	  ee [9] = volume;
	}
	else
	{
	  ee [0] = 8;
	  ee [1] = i*nrad*(nthi+1) + j*(nthi+1) + k;
	  ee [2] = ee [1] + 1;
	  ee [3] = i*nrad*(nthi+1) + k+1;
	  ee [4] = ee [3] - 1;
	  ee [5] = ee [1] + nrad*(nthi+1);
	  ee [6] = ee [2] + nrad*(nthi+1);
	  ee [7] = ee [3] + nrad*(nthi+1);
	  ee [8] = ee [4] + nrad*(nthi+1);
	  ee [9] = volume;
	}

	if (i == 0)
	{
	  ss [0] = 4;
	  ss [1] = ee [1];
	  ss [2] = ee [2];
	  ss [3] = ee [3];
	  ss [4] = ee [4];
	  ss [5] = surfaces [0];
	  ss += 6;
	}
	if (i == ndir - 1)
	{
          ss [0] = 4;
	  ss [1] = ee [5];
	  ss [2] = ee [6];
	  ss [3] = ee [7];
	  ss [4] = ee [8];
	  ss [5] = surfaces [1];
	  ss += 6;
	}
	if (k == 0)
	{
          ss [0] = 4;
	  ss [1] = ee [1];
	  ss [2] = ee [4];
	  ss [3] = ee [5];
	  ss [4] = ee [8];
	  ss [5] = surfaces [2];
	  ss += 6;
	}
        if (k == nthi - 1)
	{
          ss [0] = 4;
	  ss [1] = ee [2];
	  ss [2] = ee [3];
	  ss [3] = ee [6];
	  ss [4] = ee [7];
	  ss [5] = surfaces [3];
	  ss += 6;
	}
      }
    }
  }
  ee [0] = 0;
  ss [0] = 0;

  msh = MESH_Create (nodes, elems, surfs);

  free (nodes);
  free (elems);
  free (surfs);

  return msh;
}

/* create a copy of a mesh */
MESH* MESH_Copy (MESH *msh)
{
  ELEMENT *ele, *cpy;
  FACE *fac, *gac;
  MEM mapmem;
  MESH *ret;
  MAP *map;
  int n;

  /* initialise space */
  MEM_Init (&mapmem, sizeof (MAP), msh->surfeles_count + msh->bulkeles_count);
  ERRMEM (ret = MEM_CALLOC (sizeof (MESH)));
  MEM_Init (&ret->elemem, sizeof (ELEMENT), msh->surfeles_count + msh->bulkeles_count);
  ret->surfeles_count = msh->surfeles_count;
  ret->bulkeles_count = msh->bulkeles_count;
  for (n = 0, ele = msh->surfeles; ele; ele = ele->next)
    for (fac = ele->faces; fac; fac = fac->next) n ++; /* count surface faces */
  MEM_Init (&ret->facmem, sizeof (FACE), n);
  MEM_Init (&ret->mapmem, sizeof (MAP), MIN (msh->surfeles_count + msh->bulkeles_count, MEMCHUNK));
  ERRMEM (ret->ref_nodes = malloc (sizeof (double [3]) * (msh->nodes_count * 2)));
  ret->cur_nodes = ret->ref_nodes + msh->nodes_count;
  ret->nodes_count = msh->nodes_count;
  ret->surfeles = ret->bulkeles = NULL;
  
  /* copy vertices */
  memcpy (ret->ref_nodes, msh->ref_nodes, sizeof (double [3]) * (msh->nodes_count * 2));

  /* create, copy and map elements */
  for (ele = msh->surfeles, map = NULL; ele; ele = ele->next)
  {
    ERRMEM (cpy = MEM_Alloc (&ret->elemem));
    *cpy = *ele;

    /* element subdivision
     * does not get coppied */
    cpy->domnum = 0;
    cpy->dom = NULL;

    /* maintain list */
    cpy->prev = NULL;
    cpy->next = ret->surfeles;
    if (ret->surfeles) ret->surfeles->prev = cpy;
    ret->surfeles = cpy;
   
    /* copy surface faces */
    cpy->faces = NULL;
    for (fac = ele->faces; fac; fac = fac->next)
    {
      ERRMEM (gac = MEM_Alloc (&ret->facmem));
      *gac = *fac;
      gac->ele = cpy; /* overwrite element */

      /* integration data
       * does not get coppied */
      gac->idata = NULL;

      /* maintain list */
      gac->next = cpy->faces;
      cpy->faces = gac;
    }

    /* map new and old elements */
    MAP_Insert (&mapmem, &map, ele, cpy, NULL);
  }
  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    ERRMEM (cpy = MEM_Alloc (&ret->elemem));
    *cpy = *ele;

    /* element subdivision
     * does not get coppied */
    cpy->domnum = 0;
    cpy->dom = NULL;

    /* maintain list */
    cpy->prev = NULL;
    cpy->next = ret->bulkeles;
    if (ret->bulkeles) ret->bulkeles->prev = cpy;
    ret->bulkeles = cpy;
   
    /* map new and old elements */
    MAP_Insert (&mapmem, &map, ele, cpy, NULL);
  }

  /* maintain adjacency => use the element map */
  for (ele = ret->surfeles; ele; ele = ele->next)
  {
    for (n = 0; n < ele->neighs; n ++)
    {
      ele->adj [n] = MAP_Find (map, ele->adj [n], NULL); /* find a new pointer corresponding to the old one */
    }
  }
  for (ele = ret->bulkeles; ele; ele = ele->next)
  {
    for (n = 0; n < ele->neighs; n ++)
    {
      ele->adj [n] = MAP_Find (map, ele->adj [n], NULL);
    }
  }

  /* create mesh face list */
  for (ele = ret->surfeles; ele; ele = ele->next)
  {
    for (fac = ele->faces; fac; fac = fac->next)
    {
      fac->n = ret->faces;
      ret->faces = fac;
    }
  }

  MEM_Release (&mapmem);
  return ret;
}

/* scaling of a mesh */
void MESH_Scale (MESH *msh, double *vector)
{
  double (*ref) [3] = msh->ref_nodes,
	 (*cur) [3] = msh->cur_nodes,
	 (*end) [3] = ref + msh->nodes_count;
  ELEMENT *ele;
  FACE *fac;

  for (; ref < end; ref ++, cur ++)
  {
    ref [0][0] *= vector [0];
    ref [0][1] *= vector [1];
    ref [0][2] *= vector [2];
    cur [0][0] *= vector [0];
    cur [0][1] *= vector [1];
    cur [0][2] *= vector [2];
  }

  cur = msh->cur_nodes;
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    for (fac = ele->faces; fac; fac = fac->next)
    {
      setup_normal (cur, fac);
    }
  }
}

/* translation of a mesh */
void MESH_Translate (MESH *msh, double *vector)
{
  double (*ref) [3] = msh->ref_nodes,
	 (*cur) [3] = msh->cur_nodes,
	 (*end) [3] = ref + msh->nodes_count;

  for (; ref < end; ref ++, cur ++)
  {
    ref [0][0] += vector [0];
    ref [0][1] += vector [1];
    ref [0][2] += vector [2];
    cur [0][0] += vector [0];
    cur [0][1] += vector [1];
    cur [0][2] += vector [2];
  }
}

/* rotation of a mesh */
void MESH_Rotate (MESH *msh, double *point, double *vector, double angle)
{
  double R [9], omega [3];

  angle *=  ALG_PI / 180.0;
  COPY (vector, omega); 
  NORMALIZE (omega); 
  SCALE (omega, angle);
  EXPMAP (omega, R);

  double (*ref) [3] = msh->ref_nodes,
	 (*cur) [3] = msh->cur_nodes,
	 (*end) [3] = ref + msh->nodes_count;
  ELEMENT *ele;
  FACE *fac;

  for (; ref < end; ref ++, cur ++)
  {
    SUB (ref[0], point, omega);
    NVADDMUL (point, R, omega, ref[0]);
    SUB (cur[0], point, omega);
    NVADDMUL (point, R, omega, cur[0]);
  }

  cur = msh->cur_nodes;
  for (ele = msh->surfeles; ele; ele = ele->next)
    for (fac = ele->faces; fac; fac = fac->next)
    {
      setup_normal (cur, fac);
    }
}

/* cut through mesh with a plane; return triangulated cross-section; vertices in the triangles
 * point to the memory allocated after the triangles memory; adjacency is not maintained;
 * TRI->ptr stores a pointer to the geometrical object that has been cut by the triangle */
TRI* MESH_Cut (MESH *msh, double *point, double *normal, int *m)
{
  TRI *out, *t, *e;
  CONVEX *cvx, *q;

  cvx = MESH_Convex (msh, 0, 1);
  out = CONVEX_Cut (cvx, point, normal, m);
  for (t = out, e = t + (*m); t != e; t ++)
  {
    q = t->ptr;
    t->ptr = (TRI*) q->ele [0]; /* see (&&&) */
  }
  CONVEX_Destroy (cvx);

  return out;
}

/* as above but this time the plane and the cut are in the reference configuration */
TRI* MESH_Ref_Cut (MESH *msh, double *point, double *normal, int *m)
{
  TRI *out, *t, *e;
  CONVEX *cvx, *q;

  cvx = MESH_Convex (msh, 1, 1);
  out = CONVEX_Cut (cvx, point, normal, m);
  for (t = out, e = t + (*m); t != e; t ++)
  {
    q = t->ptr;
    t->ptr = (TRI*) q->ele [0]; /* see (&&&) */
  }
  CONVEX_Destroy (cvx);

  return out;
}

/* split mesh in two with plane defined by (point, normal); output meshes are tetrahedral if some
 * elements are crossed; if only element boundaries are crossed then the original mesh is used;
 * topoadj != 0 implies cutting from the point and through the topological adjacency only */
#if OLD_SPLIT
void MESH_Split (MESH *msh, double *point, double *normal, short topoadj, int surfid, MESH **one, MESH **two)
{
  TRI *c, *b, *a, *t, *e, *q;
  int mc, mb, ma, mq;
  short onepart;
  MESH *tmp;
  KDT *kd;

  c = MESH_Cut (msh, point, normal, &mc);

  if (topoadj)
  {
    TRI_Compadj (c, mc); /* compute adjacency structure of triangles */
    b = TRI_Topoadj (c, mc, point, &mb); /* part of the triangulation adjacent to the point */
    ASSERT_DEBUG (b, "Input point is too far from mesh section by the splitting plane");
    if (mb < mc) onepart = 1;
    else onepart = 0;
    mc = mb;
    c = b;
  }

  for (t = c, e = t + mc; t != e; t ++) t->flg = surfid;

  if (c)
  {
    if (topoadj && onepart)
    {
      kd = TRI_Kdtree (c, mc);

      *two = NULL; /* this is certain since TRI_Topoadj returend less triangles then passed */

      if (!inter_element_local_split (msh, kd, point, normal, surfid, one)) /* if regular mesh interl-element split failed */
      {
	b = trimadj (msh, kd, point, normal, &mb); /* trim faces adjacent to the input triangulation */

	q = TRI_Merge (b, mb, c, mc, &mq);

	tmp = tetrahedralize3 (q, mq, 0.0, 2.0, msh->surfeles->volume); /* the internal surface is included */

	/* TODO: map volume materials */

	ASSERT_DEBUG_EXT (
	    inter_element_local_split (tmp, kd, point, normal, surfid, one),
	    "Local inter-element splitting failed"); /* this will be an inter-element split now */

        MESH_Destroy (tmp);
        free (q);
      }

      KDT_Destroy (kd);
    }
    else if (inter_element_global_split (msh, point, normal, surfid, one, two)) /* global inter-element split succeeded */
    {
      /* all done in the condition */
    }
    else
    {
      trim (msh, point, normal, &b, &mb, &a, &ma);

      q = TRI_Merge (b, mb, c, mc, &mq);
#if DEBUG
      TRI_Compadj (q, mq); /* XXX: tests triangulation consistency */
#endif
      *one = tetrahedralize3 (q, mq, 0.0, 2.0, msh->surfeles->volume);
      free (q);

      /* TODO: map volume materials */
   
      q = TRI_Merge (a, ma, c, mc, &mq);
#if DEBUG
      TRI_Compadj (q, mq); /* XXX: tests triangulation consistency */
#endif
      *two = tetrahedralize3 (q, mq, 0.0, 2.0, msh->surfeles->volume);
      free (q);

      /* TODO: map volume materials */

      free (b);
      free (a);
    }

    free (c);
  }
  else
  {
    double d [3];

    SUB (msh->cur_nodes [0], point, d);
    if (DOT (d, normal) <= 0.0)
    {
      *one = MESH_Copy (msh);
      *two = NULL;
    }
    else
    {
      *two = MESH_Copy (msh);
      *one = NULL;
    }
  }
}
#else
static void mark_adjacent_in_set (ELEMENT *ele, SET *set)
{
  ele->flag = 1;

  for (int i = 0; i < ele->neighs; i ++)
  {
    if (SET_Contains (set, ele->adj[i], NULL) && ele->adj[i]->flag == 0)
    {
      mark_adjacent_in_set (ele->adj[i], set);
    }
  }
}

void MESH_Split (MESH *msh, double *point, double *normal, short topoadj, int surfid, MESH **one, MESH **two)
{
  double (*nod) [3] = msh->cur_nodes, nn [3];
  SET *below, *above, *onpla, *item;
  int code, bulk, on, i, j, k, l;
  MEM mapmem, setmem;
  MAP *map, *jtem;
  ELEMENT *ele;
  FACE *fac;

  MEM_Init (&setmem, sizeof (SET), MEMCHUNK);
  MEM_Init (&mapmem, sizeof (MAP), MEMCHUNK);

  COPY (normal, nn);
  NORMALIZE (nn);

  below = above = onpla = NULL;
  map = NULL;

  for (bulk = 0, ele = msh->surfeles; ele;)
  {
    for (code = on = i = 0; i < ele->type; i ++)
    {
      double a [3], dot;

      SUB (nod [ele->nodes [i]], point, a);
      dot = DOT (a, nn);
      if (dot < -GEOMETRIC_EPSILON)
      {
	if (!code) code = -1;
	else if (code > 0) { code = 0; break; }
      }
      else if (dot > GEOMETRIC_EPSILON)
      {
	if (!code) code = 1;
	else if (code < 0) { code = 0; break; }
      }
      else on = 1;
    }


    if (code <= 0) /* element crossing the plane or below */
    {
      SET_Insert (&setmem, &below, ele, NULL);
    }
    else /* element above the plane  */
    {
      SET_Insert (&setmem, &above, ele, NULL);
      if (i == ele->type && on) SET_Insert (&setmem, &onpla, ele, NULL); /* on plane */
    }

    if (!bulk && !ele->next) bulk = 1, ele = msh->bulkeles;
    else ele = ele->next;
  }

  if (topoadj == 0) /* produce two new meshes from two sets */
  {
    *one = produce_subset_mesh (msh, below, surfid); /* the undefined, new faces will have surfid */
    *two = produce_subset_mesh (msh, above, surfid);
  }
  else /* produce new single mesh by doubling the nodes along the faces splitting the two sets */
  {
    for (ele = msh->surfeles; ele; ele = ele->next) ele->flag = 0;
    for (ele = msh->bulkeles; ele; ele = ele->next) ele->flag = 0; /* unmark all elements */

    for (item = SET_First (onpla); item; item = SET_Next (item))
    {
      ele = item->data;
      if (ELEMENT_Contains_Spatial_Point (msh, ele, point)) break; /* find first element adjacent to the input point */
    }

    if (item) /* found */
    {
      mark_adjacent_in_set (ele, onpla); /* mark all on plane elements topologically adjacent to that first element */

      double (*nodes) [3]; /* new mesh definition nodes */
      int *elements, *e; /* new mesh definition elements, current pointer */
      int *surfaces, *s; /* new mesh definition faces, current pointer */

      ERRMEM (nodes = malloc (2 * msh->nodes_count * (sizeof (double [3])))); /* overestimate twice */
      ERRMEM (elements = malloc ((msh->surfeles_count + msh->bulkeles_count + 1) * sizeof (int [10]))); /* overestimate a bit */
      ERRMEM (surfaces = malloc (8 * (msh->surfeles_count + msh->bulkeles_count) * sizeof (int [6]))); /* overestimate considerabely */
      e = elements;
      s = surfaces;
      *s = surfid; /* new faces will be left undefined and upon rediscovery in MESH_Create will have this value */
      s ++;
      k = msh->nodes_count; /* current new free node number at the end of old range */

      for (bulk = 0, ele = msh->surfeles; ele;)
      {
	if (ele->flag) /* if topologically adjacent to the input point */
	{
	  for (i = 0; i < ele->neighs; i ++) /* for all neighbours */
	  {
	    if (SET_Contains (below, ele->adj[i], NULL)) /* if neighbour in the other set (on plane is in above) */
	    {
	      for (j = 0; j < ele->type; j ++)
	      {
		for (l = 0; l < ele->adj[i]->type; l ++)
		{
		  if (ele->nodes [j] == ele->adj[i]->nodes [l]) /* common node needs to be doubled */
		  {
		    if (MAP_Find_Node (map, (void*) (long) ele->nodes [j], NULL) == NULL) /* if new node was not mapped yet */
		    {
		      MAP_Insert (&mapmem, &map, (void*) (long) ele->nodes [j], (void*) (long) k, NULL); /* map new node number */
		      k ++; /* increase new free node number */
		    }
		  }
		}
	      }
	    }
	  }

	  *e = ele->type; e ++;
	  for (i = 0; i < ele->type; i ++, e ++)
	  {
	    jtem = MAP_Find_Node (map, (void*) (long) ele->nodes [i], NULL); /* try finding new node mapping */
	    if (jtem) *e = (int) (long) jtem->data; /* use new node */
	    else *e = ele->nodes [i]; /* use old node */
	  }
	  *e = ele->volume; e ++;

          for (fac = ele->faces; fac; fac = fac->next) /* copy existing faces */
	  {
	    *s = fac->type; s ++;
	    for (i = 0; i < fac->type; i ++, s ++)
	    {
	      jtem = MAP_Find_Node (map, (void*) (long) fac->nodes [i], NULL); /* try finding new node mapping */
	      if (jtem) *s = (int) (long) jtem->data; /* use new node */
	      else *s = fac->nodes [i]; /* use old node */
	    }
	    *s = fac->surface; s ++;
	  }
	}
	else /* regular element */
	{
	  *e = ele->type; e ++;
	  for (i = 0; i < ele->type; i ++, e ++) *e = ele->nodes [i];
	  *e = ele->volume; e ++;

	  for (fac = ele->faces; fac; fac = fac->next) /* copy existing faces */
	  {
	    *s = fac->type; s ++;
	    for (i = 0; i < fac->type; i ++, s ++) *s = fac->nodes [i];
	    *s = fac->surface; s ++;
	  }
	}

	if (!bulk && !ele->next) bulk = 1, ele = msh->bulkeles;
	else ele = ele->next;
      }
      *e = 0;
      *s = 0;

      for (i = 0; i < msh->nodes_count; i ++) /* copy old nodes */
      {
	COPY (msh->cur_nodes [i], nodes [i]);
      }
      for (jtem = MAP_First (map); jtem; jtem = MAP_Next (jtem)) /* copy new nodes */
      {
	i = (int) (long) jtem->key;
	j = (int) (long) jtem->data;

	COPY (msh->cur_nodes [i], nodes [j]);
      }

      *one = MESH_Create (nodes, elements, surfaces); /* create new mesh */
      *two = NULL;

      free (nodes);
      free (elements);
      free (surfaces);
    }
    else *one = *two = NULL; /* topologically adjacent split has failed */
  }

  MEM_Release (&setmem);
  MEM_Release (&mapmem);
}
#endif

/* is mesh separable into disjoint parts */
int MESH_Separable (MESH *msh)
{
  ELEMENT *ele;

  for (ele = msh->surfeles; ele; ele = ele->next) ele->flag = 0;
  for (ele = msh->bulkeles; ele; ele = ele->next) ele->flag = 0;

  int flag = 1;

  mark_neighs (msh->surfeles, flag);

  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    if (ele->flag == 0)
    {
      flag ++;
      mark_neighs (ele, flag);
    }
  }

  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    if (ele->flag == 0)
    {
      flag ++;
      mark_neighs (ele, flag);
    }
  }

  return (flag == 1 ? 0 : flag);
}

/* separate mesh into disjoint parts */
MESH** MESH_Separate (MESH *msh, int *m)
{
  ELEMENT *ele;
  MESH **out;
  SET *els;
  MEM mem;
  int i;
 
  *m = MESH_Separable (msh);

  if ((*m) == 0) return NULL;

  ERRMEM (out = malloc ((*m) * sizeof (MESH*)));
  MEM_Init (&mem, sizeof (SET), MEMCHUNK);

  for (i = 1, els = NULL; i <= (*m); i ++)
  {
    SET_Free (&mem, &els);

    for (ele = msh->surfeles; ele; ele = ele->next)
      if (ele->flag == i) SET_Insert (&mem, &els, ele, NULL);

    for (ele = msh->bulkeles; ele; ele = ele->next)
      if (ele->flag == i) SET_Insert (&mem, &els, ele, NULL);

#if OLD_SPLIT
    int j = 0;
    out [i-1] = produce_split_mesh (msh, els, &j, 0, NULL, NULL);
#else
    out [i-1] = produce_subset_mesh (msh, els, 0);
#endif
  }

  MEM_Release (&mem);

  return out;
}

/* compute partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void MESH_Char_Partial (MESH *msh, int ref, double *vo, double *sx, double *sy, double *sz, double *eul)
{
  ELEMENT *ele;

  /* previously we looped over the surface faces and used simplex integration,
   * but since partitioned meshes do not have all surface elements (in order to
   * decrease contact detection effort) we now iterate over all elements */

  for (ele = msh->surfeles; ele; ele = ele->next)
    ELEMENT_Char_Partial (msh, ele, ref, vo, sx, sy, sz, eul);

  for (ele = msh->bulkeles; ele; ele = ele->next)
    ELEMENT_Char_Partial (msh, ele, ref, vo, sx, sy, sz, eul);
}

/* get characteristics of the meshed shape:
 * volume, mass center, and Euler tensor (centered) */
void MESH_Char (MESH *msh, int ref, double *volume, double *center, double *euler)
{
  double vo, sx, sy, sz,
	 cen [3], eul [9];

  vo = sx = sy = sz = 0.0;
  SET9 (eul, 0.0);

  MESH_Char_Partial (msh, ref, &vo, &sx, &sy, &sz, eul);

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

/* find an element containing a spatial or referential point */
ELEMENT* MESH_Element_Containing_Point (MESH *msh, double *point, int ref)
{
  /* XXX:
   * the linear search should be at some point optimised with a proper
     spatial search => move point to the reference configuration and
     query a search tree based on the reference configuration geometry */

  ELEMENT *ele;

  /* first search surface elements */
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    if (ELEMENT_Contains_Point (msh, ele, point, ref)) return ele;
  }

  /* then the bulk elements */
  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    if (ELEMENT_Contains_Point (msh, ele, point, ref)) return ele;
  }

  return NULL;
}

/* find an element containing a spatial point */
ELEMENT* MESH_Element_Containing_Spatial_Point (MESH *msh, double *point)
{
  return MESH_Element_Containing_Point (msh, point, 0);
}

/* find an element with a given node */
ELEMENT* MESH_Element_With_Node (MESH *msh, int node)
{
  ELEMENT *ele = MAP_Find (msh->map, (void*) (long) node, NULL),
	  *start [2] = {msh->surfeles, msh->bulkeles};
  int i, j;

  if (ele == NULL) /* if this node has not been yet mapped to an element */
  {
    for(j = 0; j < 2; j ++)
    {
      for (ele = start [j]; ele; ele = ele->next) /* make a costly linear search for an element */
      {
	for (i = 0; i < ele->type; i ++)
	{
	  if (ele->nodes [i] == node)
	  {
	    MAP_Insert (&msh->mapmem, &msh->map, (void*) (long) node, ele, NULL); /* and map it */
	    return ele;
	  }
	}
      }
    }
  }

  return ele;
}

/* collect elements around a node (ele->node [i] == node && *set == NULL initially assumed) */
void MESH_Elements_Around_Node (ELEMENT *ele, int node, SET **set)
{
  ELEMENT *nei;
  int i, j;

  SET_Insert (NULL, set, ele, NULL);

  for (i = 0; i < ele->neighs; i ++)
  {
    nei = ele->adj [i];
    for (j = 0; j < nei->type; j ++)
    {
      if (nei->nodes [j] == node)
      {
	if (!SET_Contains (*set, nei, NULL))
	{ 
	  MESH_Elements_Around_Node (nei, node, set);
	}
      }
    }
  }
}

/* find an element containing a referential point */
/* update mesh according to the given motion */
void MESH_Update (MESH *msh, void *body, void *shp, MOTION motion)
{
  ELEMENT *ele;
  FACE *fac;
  double (*ref) [3] = msh->ref_nodes,
	 (*cur) [3] = msh->cur_nodes;
  int m = msh->nodes_count, n;

  if (motion)
  {
    for (n = 0; n < m; n ++)
    {
      SGP sgp = {shp, NULL, GOBJ_DUMMY, NULL}; /* move current nodes (NULL for gobj implies nodal update) */
      motion (body, &sgp, ref [n], cur [n]);
    }
  }
  else /* restore reference configuration */
  {
    for (n = 0; n < m; n ++)
    {
      COPY (ref [n], cur [n]);
    }
  }

  for (ele = msh->surfeles; ele; ele = ele->next)
    for (fac = ele->faces; fac; fac = fac->next)
    {
      setup_normal (cur, fac); /* update normals */
    }
}

/* convert mesh into a list of convices;
 * ref > 0 => create referential mesh image;
 * if ele0ptr > 0 => cvx->ele [0] points to the source element;
 * otherwise => create current mesh image;
 * CONVEX->ele[0] == corresponding element */
CONVEX* MESH_Convex (MESH *msh, int ref, int ele0ptr)
{
  int vfac [30], surfaces [6], idx [6], nfac, *f, n, i;
  double nodes [8][3];
  ELEMENT *ele, *next;
  short bulk = 0;
  CONVEX *cvx;
  FACE *fac;

  cvx = NULL;
  for (ele = msh->surfeles; ele; ele = next)
  {
    load_nodes (ref ? msh->ref_nodes : msh->cur_nodes, ele->type, ele->nodes, nodes);
    nfac = neighs (ele->type); /* number of faces */

    for (i = 0; i < 6; i ++) idx [i] = 1; /* mark all faces */

    for (fac = ele->faces, f = vfac, n = 0; fac; fac = fac->next, n ++) /* surface faces first */
    {
      surfaces [n] = fac->surface; /* set surface identifiers */
      f = setup_face_vertices (ele, fac->index, f); /* write face vertex indices */
      idx [fac->index] = 0; /* unmark */
    }

    for (; n < nfac; n ++) /* internal faces next */
    {
      surfaces [n] = -INT_MAX; /* invalid pointer for internal faces */
      for (i = 0; i < 6; i ++)
      {
	if (idx [i]) /* mareked => spare internal face index */
	{
          f = setup_face_vertices (ele, i, f); /* write face vertex indices */
	  idx [i] = 0; /* unrmark */
	}
      }
    }

    cvx = CONVEX_Create (cvx, (double*)nodes, ele->type, vfac, nfac, surfaces, ele->volume); /* add new convex to the list */

    if (ele0ptr)
    {
      ERRMEM (cvx->ele = malloc (sizeof (ELEMENT*)));
      cvx->ele [0] = ele; /* store the corresponding element */ /* (&&&) */
      cvx->nele = 1;
    }

    if (ele->next) next = ele->next;
    else if (!bulk) { next = msh->bulkeles; bulk = 1; }
    else next = NULL;
  }

  return cvx;
}

/* compute extents of entire mesh */
void MESH_Extents (MESH *msh, double *extents)
{
  ELEMENT *ele;
  double e [6];

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;
    
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    element_extents (msh, ele, 0, e);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [3] > extents [3]) extents [3] = e [3];
    if (e [4] > extents [4]) extents [4] = e [4];
    if (e [5] > extents [5]) extents [5] = e [5];
  }

  extents [0] -= GEOMETRIC_EPSILON;
  extents [1] -= GEOMETRIC_EPSILON;
  extents [2] -= GEOMETRIC_EPSILON;
  extents [3] += GEOMETRIC_EPSILON;
  extents [4] += GEOMETRIC_EPSILON;
  extents [5] += GEOMETRIC_EPSILON;
}

/* compute oriented extents of entire mesh */
void MESH_Oriented_Extents (MESH *msh, double *vx, double *vy, double *vz, double *extents)
{
  double (*cur) [3], (*end) [3];
  double e [3];

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;
    
  for (cur = msh->cur_nodes, end = cur + msh->nodes_count; cur < end; cur ++)
  {
    e [0] = DOT (vx, cur[0]);
    e [1] = DOT (vy, cur[0]);
    e [2] = DOT (vz, cur[0]);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [0] > extents [3]) extents [3] = e [0];
    if (e [1] > extents [4]) extents [4] = e [1];
    if (e [2] > extents [5]) extents [5] = e [2];
  }
}

/* return first not NULL bulk material of an element */
void* MESH_First_Bulk_Material (MESH *msh)
{
  ELEMENT *ele;

  for (ele = msh->bulkeles; ele; ele = ele->next)
    if (ele->mat) return ele->mat;

  for (ele = msh->surfeles; ele; ele = ele->next)
    if (ele->mat) return ele->mat;

  return NULL;
}

/* Metis recursive partitioning (<= 8 partitions) */
void METIS_PartGraphRecursive (int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *); 

/* Metis k-way partitioning (> 8 partitions) */
void METIS_PartGraphKway (int *, int *, int *, int *, int *, int *, int *, int *, int *, int *, int *); 

/* partitioning related element copy */
static ELEMENT* copy_element (MEM *elemem, ELEMENT *ele, int part, MEM *facmem)
{
  FACE *fac, *gac;
  ELEMENT *out;
  int i;

  ERRMEM (out = MEM_Alloc (elemem));
  out->type = ele->type;
  for (i = 0; i < ele->type; i ++) out->nodes [i] = ele->nodes [i];
  out->volume = ele->volume;
  out->mat = ele->mat;
  for (i = 0; i < ele->neighs; i ++)
  {
    if (ele->adj [i]->domnum == part)
    {
      out->adj [out->neighs ++] = ele->adj [i]; /* to be later repaced by out->adj[]->dom (***) */
    }
  }

  ele->dom = (TRISURF*)out; /* map copy to original (***) */

  for (fac = ele->faces; fac; fac = fac->next)
  {
    ERRMEM (gac = MEM_Alloc (facmem));
    COPY (fac->normal, gac->normal);
    gac->type = fac->type;
    for (i = 0; i < fac->type; i ++) gac->nodes [i] = fac->nodes [i];
    gac->index = fac->index;
    gac->surface = fac->surface;
    gac->ele = out;
    gac->next = out->faces;
    out->faces = gac;
  }

  return out;
}

/* starting with ele gather elements with node and with ele->domnum != domnum into the out set */
static void elements_with_node_and_not_domnum (int node, int domnum, ELEMENT *ele, MEM *mem, SET **out)
{
  if (!SET_Contains (*out, ele, NULL))
  {
    SET_Insert (mem, out, ele, NULL);
    for (int i = 0; i < ele->neighs; i ++)
    {
      if (ele->adj[i]->domnum != domnum)
      {
	for (int j = 0; j < ele->adj [i]->type; j ++)
	{
	  if (ele->adj [i]->nodes [j] == node)
	    elements_with_node_and_not_domnum (node, domnum, ele->adj [i], mem, out);
	}
      }
    }
  }
}

/* partition mesh; return the resultant mesh parts; output a table of tuples (m1, m2, n1, n2) of gluing nodes,
 * where m1, m2 is a pair of the output meshes and n1, n2 are their corresponding coincident nodes;
 * additionally output tuples (m1, m2, e1, e2) of topologically adjacent surface element
 * pairs from the partitions boundaries (indexed as stored in lists and outputed by SGP_Create);
 * upon exit the 'domnum' element values of the input mesh indicate destination partitions of the elements */
MESH** MESH_Partition (MESH *msh, int nparts, int *numglue, int **gluenodes, int *numadj, int **adjeles)
{
  int i, j, n, m, *xadj, *adjncy, *vwgt, wgtflag, numflag, options, edgecut, *part;
  SET *aux, **n2p, *jtem, *ktem;
  ELEMENT *ele, *nel;
  MEM setmem, mapmem;
  MAP **nod, *item;
  double *x, *y;
  MESH **out;
  FACE *fac;

  for (n = 0, ele = msh->surfeles; ele; ele = ele->next) ele->domnum = n ++; /* number elements */
  for (ele = msh->bulkeles; ele; ele = ele->next) ele->domnum = n ++;

  ASSERT_DEBUG (n >= nparts, "Number of elements is smaller than the number of partitions");

  for (m = 0, ele = msh->surfeles; ele; ele = ele->next) m += ele->neighs; /* compute adjacency size */
  for (ele = msh->bulkeles; ele; ele = ele->next) m += ele->neighs;

  ERRMEM (xadj = MEM_CALLOC (sizeof (int [n+1])));
  ERRMEM (adjncy = MEM_CALLOC (sizeof (int [m])));
  ERRMEM (vwgt = MEM_CALLOC (sizeof (int [n])));
  ERRMEM (part = MEM_CALLOC (sizeof (int [n])));
  ERRMEM (out = MEM_CALLOC (nparts * sizeof (MESH*)));

  for (n = 0, ele = msh->surfeles; ele; ele = ele->next, n ++)
  {
    vwgt [n] = ele->type;
    xadj [n+1] = xadj [n] + ele->neighs;
    for (m = 0; m < ele->neighs; m ++) adjncy [xadj [n]+m] = ele->adj [m]->domnum;
  }
  for (ele = msh->bulkeles; ele; ele = ele->next, n ++)
  {
    vwgt [n] = ele->type;
    xadj [n+1] = xadj [n] + ele->neighs;
    for (m = 0; m < ele->neighs; m ++) adjncy [xadj [n]+m] = ele->adj [m]->domnum;
  }

  options = 0;
  wgtflag = 2;
  numflag = 0;

  if (nparts <= 8) /* partition */
    METIS_PartGraphRecursive (&n, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag, &nparts, &options, &edgecut, part);
  else METIS_PartGraphKway (&n, xadj, adjncy, vwgt, NULL, &wgtflag, &numflag, &nparts, &options, &edgecut, part);

  ELEMENT *head[] = {msh->surfeles, msh->bulkeles};

  for (m = 0; m < nparts; m ++) /* initialize output meshes */
  {
    ERRMEM (out [m] = MEM_CALLOC (sizeof (MESH)));
    MEM_Init (&out [m]->elemem, sizeof (ELEMENT), n);
    MEM_Init (&out [m]->facmem, sizeof (FACE), MEMCHUNK);
    MEM_Init (&out [m]->mapmem, sizeof (MAP), MIN (n, MEMCHUNK));
  }

  for (n = 0; n < 2; n ++) for (ele = head [n]; ele; ele = ele->next) /* set destination partitions */
  { 
    ele->domnum = part [ele->domnum];
  }

  for (n = 0; n < 2; n ++) for (ele = head [n]; ele; ele = ele->next) /* copy elements */
  { 
    m = ele->domnum; /* destination partition */

    nel = copy_element (&out [m]->elemem, ele, m, &out [m]->facmem);

    if (ele->faces)
    {
      out [m]->surfeles_count ++;
      nel->next = out [m]->surfeles;
      if (out [m]->surfeles) out [m]->surfeles->prev = nel;
      out [m]->surfeles = nel;
    }
    else
    {
      out [m]->bulkeles_count ++;
      nel->next = out [m]->bulkeles;
      if (out [m]->bulkeles) out [m]->bulkeles->prev = nel;
      out [m]->bulkeles = nel;
    }
  }

  for (m = 0; m < nparts; m ++) /* map adjacency */
  {
    ELEMENT *head[] = {out [m]->surfeles, out [m]->bulkeles};

    for (n = j = 0; n < 2; n ++) for (ele = head [n]; ele; ele = ele->next)
    {
      if (n == 0) ele->domnum = j ++; /* number surface elements (@@@) */

      for (i = 0; i < ele->neighs; i ++) ele->adj [i] = (ELEMENT*)ele->adj [i]->dom;
    }
  }

  j = 256;
  *numadj = 0;
  ERRMEM (*adjeles = malloc (j * sizeof (int [4])));
  MEM_Init (&setmem, sizeof (SET), MEMCHUNK);

  for (ele = msh->surfeles; ele; ele = ele->next) /* output partition boundary surface elements pairs */
  {
    for (i = 0; i < ele->neighs; i ++)
    {
      if (ele->domnum != ele->adj[i]->domnum && ele->adj[i]->faces)
      {
	for (int k = 0; k < ele->type; k ++)
	{
	  for (int l = 0; l < ele->adj[i]->type; l ++)
	  {
	    if (ele->nodes [k] == ele->adj [i]->nodes [l])
	    {
	      aux = NULL;
	      elements_with_node_and_not_domnum (ele->nodes [k], ele->domnum, ele->adj [i], &setmem, &aux);

	      for (jtem = SET_First (aux); jtem; jtem = SET_Next (jtem))
	      {
		(*numadj) ++;

		if (*numadj == j)
		{
		  j *= 2;
		  ERRMEM (*adjeles = realloc (*adjeles, j * sizeof (int [4])));
		}

		int *e = &(*adjeles)[4*((*numadj) - 1)];

		nel = jtem->data;
		e [0] = ele->domnum;
		e [1] = nel->domnum;
		e [2] = ((ELEMENT*)ele->dom)->domnum; /* (@@@) */
		e [3] = ((ELEMENT*)nel->dom)->domnum; /* (@@@) */
	      }

	      SET_Free (&setmem, &aux);
	    }
	  }
	}
      }
    }
  }

  for (ele = msh->surfeles; ele; ele = ele->next) ele->dom = NULL; /* clear pointers */
  for (ele = msh->bulkeles; ele; ele = ele->next) ele->dom = NULL;

  MEM_Init (&mapmem, sizeof (MAP), MEMCHUNK);
  ERRMEM (nod = MEM_CALLOC (nparts * sizeof (MAP*)));
  ERRMEM (n2p = MEM_CALLOC (msh->nodes_count * sizeof (SET*)));

  for (m = 0; m < nparts; m ++) /* map nodes */
  {
    ELEMENT *head[] = {out [m]->surfeles, out [m]->bulkeles};
    j = 0; /* local node numbering */

    for (n = 0; n < 2; n ++) for (ele = head [n]; ele; ele = ele->next)
    {
      if (ele->domnum) ele->domnum = 0; /* unmark (@@@) */

      for (i = 0; i < ele->type; i ++)
      {
        SET_Insert (&setmem, &n2p [ele->nodes [i]], (void*) (long) m, NULL); /* map old node number to partition */

	if ((item = MAP_Find_Node (nod [m], (void*) (long) ele->nodes [i], NULL)))
	{
	  ele->nodes [i] = (int) (long) item->data; /* map local number */
	}
	else
	{
	  MAP_Insert (&mapmem, &nod [m], (void*) (long) ele->nodes [i], (void*) (long) j, NULL);
	  ele->nodes [i] = j ++; /* increment */
	}
      }

      for (fac = ele->faces; fac; fac = fac->next)
      {
	for (i = 0; i < fac->type; i ++)
	{
	  if ((item = MAP_Find_Node (nod [m], (void*) (long) fac->nodes [i], NULL)))
	  {
	    fac->nodes [i] = (int) (long) item->data;
	  }
	  else
	  {
	    MAP_Insert (&mapmem, &nod [m], (void*) (long) fac->nodes [i], (void*) (long) j, NULL);
	    fac->nodes [i] = j ++ /* increment */;
	  }
	}
      }
    }

    ERRMEM (out [m]->ref_nodes = malloc (2 * j * sizeof (double [3])));
    out [m]->cur_nodes = out [m]->ref_nodes + j;
    out [m]->nodes_count = j;

    for (item = MAP_First (nod [m]); item; item = MAP_Next (item))
    {
      i = (int) (long) item->key;
      j = (int) (long) item->data;

      x = msh->ref_nodes [i],
      y = out [m]->ref_nodes [j];
      COPY (x, y);

      x = msh->cur_nodes [i],
      y = out [m]->cur_nodes [j];
      COPY (x, y);
    }
  }


  j = 256;
  *numglue = 0;
  ERRMEM (*gluenodes = malloc (j * sizeof (int [4])));

  for (i = 0; i < msh->nodes_count; i ++) /* generate gluing nodes pairs */
  {
    ASSERT_DEBUG (SET_Size (n2p [i]) >= 1, "Not all nodes has been mapped");

    for (jtem = SET_First (n2p [i]); jtem; jtem = SET_Next (jtem))
    {
      for (ktem = SET_First (n2p [i]); ktem; ktem = SET_Next (ktem))
      {
	if (jtem->data < ktem->data)
	{
	  (*numglue) ++;

	  if (*numglue == j)
	  {
	    j *= 2;
	    ERRMEM (*gluenodes = realloc (*gluenodes, j * sizeof (int [4])));
	  }

	  int *e = &(*gluenodes)[4*((*numglue) - 1)];

	  e [0] = (int) (long) jtem->data;
	  e [1] = (int) (long) ktem->data;
	  e [2] = (int) (long) MAP_Find (nod [e[0]], (void*) (long) i, NULL);
	  e [3] = (int) (long) MAP_Find (nod [e[1]], (void*) (long) i, NULL);
	}
      }
    }
  }

  for (i = 0; i < nparts; i ++) /* create mesh faces list */
  {
    FACE faces [6], *fac;
    int o;

    msh = out [i];

    ELEMENT *list [] = {msh->surfeles, msh->bulkeles};

    for (j = 0; j < 2; j ++)
    for (ele = list [j]; ele; ele = ele->next)
    {
      m = neighs (ele->type); 

      if (ele->neighs < m) /* discover new faces */
      {
	for (n = 0; n < m; n ++)
	  setup_face (ele, n, &faces [n], 0);

	for (n = 0; n < m; n ++) /* for each potential free face */
	{
	  for (fac = ele->faces; fac; fac = fac->next)
	    if (face_compare (&faces [n], fac) == 0) break; /* surface face found */

	  for (o = 0; o < ele->neighs; o ++)
	    if (element_has_nodes (ele->adj [o], faces [n].type, faces [n].nodes)) break; /* neighbor found */

	  if (fac == NULL && o == ele->neighs) /* neither a surface face nor having a neighbor => free surface */
	  {
	    ERRMEM (fac = MEM_Alloc (&msh->facmem));
	    fac->type = faces [n].type;
	    for (o = 0; o < fac->type; o ++) fac->nodes [o] = faces [n].nodes [o];
	    fac->index = n;
	    fac->ele = ele;
	    setup_normal (msh->ref_nodes, fac);
	    fac->next = ele->faces;
	    ele->faces = fac;
	  }
	}
      }

      for (fac = ele->faces; fac; fac = fac->next) /* append mesh faces list */
      {
	fac->n = msh->faces;
	msh->faces = fac;
      }
    }
  }

  free (xadj);
  free (adjncy);
  free (vwgt);
  free (part);
  free (nod);
  free (n2p);
  MEM_Release (&mapmem);
  MEM_Release (&setmem);

  return out;
}

/* delete element set */
void MESH_Delete_Elements (MESH *msh, SET *elements)
{
  /* FIXME / TODO */

#if 0
  int *node_map, m;

  ERRMEM (node_map = MEM_CALLOC (msh->nodes_count * sizeof (int)));

  for (ele = msh->surfeles; ele; ele = ele->next)
    for (n = 0; n < ele->type; n ++) node_map [ele->nodes [n]] ++;

  for (ele = msh->bulkeles; ele; ele = ele->next)
    for (n = 0; n < ele->type; n ++) node_map [ele->nodes [n]] ++;

  for (n = 0, m = 1; n < msh->nodes_count; n ++)
    if (node_map [n]) node_map [n] = m ++;

  if (m < (n+1)) /* there are not referenced nodes */
  {
    for (ele = msh->surfeles; ele; ele = ele->next)
    {
      for (n = 0; n < ele->type; n ++) ele->nodes [n] = node_map [ele->nodes [n]] - 1;
      for (fac = ele->faces; fac; fac = fac->next) for (n = 0; n < fac->type; n ++) fac->nodes [n] = node_map [fac->nodes [n]] - 1;
    }

    for (ele = msh->bulkeles; ele; ele = ele->next)
      for (n = 0; n < ele->type; n ++) ele->nodes [n] = node_map [ele->nodes [n]] - 1;

    double (*ref) [3], (*mref) [3] = msh->ref_nodes,
	   (*cur) [3], (*mcur) [3] = msh->cur_nodes;

    ERRMEM (ref = malloc (2 * (m-1) * sizeof (double [3])));
    cur = ref + (m-1);

    for (n = 0; n < msh->nodes_count; n ++)
    {
      if (node_map [n])
      { 
	COPY (mref [n], ref [node_map [n] - 1]);
	COPY (mcur [n], cur [node_map [n] - 1]);
      }
    }

    free (msh->ref_nodes);

    msh->nodes_count = m-1;
    msh->ref_nodes = ref;
    msh->cur_nodes = cur;
  }

  free (node_map);
#endif
}

/* free mesh memory */
void MESH_Destroy (MESH *msh)
{
  ELEMENT *ele;
  FACE *fac;
  int n;

  for (fac = msh->faces; fac; fac = fac->n)
  {
    if (fac->idata) free (fac->idata);
  }

  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    if (ele->dom)
    {
      for (n = 0; n < ele->domnum; n ++) free (ele->dom [n].tri);
      free (ele->dom);
    }
    if (ele->state) free (ele->state);
  }

  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    if (ele->dom)
    {
      for (n = 0; n < ele->domnum; n ++) free (ele->dom [n].tri);
      free (ele->dom);
    }
    if (ele->state) free (ele->state);
  }

  MEM_Release (&msh->facmem);
  MEM_Release (&msh->elemem);
  MEM_Release (&msh->mapmem);
  free (msh->ref_nodes);
  free (msh);
}

/* does the element contain a spatial or referential point? */
int ELEMENT_Contains_Point (MESH *msh, ELEMENT *ele, double *point, int ref)
{
  double nodes [8][3], q [3], d;

  load_nodes (ref ? msh->ref_nodes : msh->cur_nodes, ele->type, ele->nodes, nodes);

  d = gjk_convex_point ((double*)nodes, ele->type, point, q);

  return d <= GEOMETRIC_EPSILON;
}

/* does the element contain a spatial point? */
int ELEMENT_Contains_Spatial_Point (MESH *msh, ELEMENT *ele, double *point)
{
  double nodes [8][3], q [3], d;

  load_nodes (msh->cur_nodes, ele->type, ele->nodes, nodes);

  d = gjk_convex_point ((double*)nodes, ele->type, point, q);

  return d <= GEOMETRIC_EPSILON;
}

/* return >= node index if point == node[index] or -1 otherwise */
int ELEMENT_Ref_Point_To_Node (MESH *msh, ELEMENT *ele, double *point)
{
  double nodes [8][3], d [3], a;
  int  i;

  load_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  for (i = 0; i < ele->type; i ++)
  {
    SUB (nodes [i], point, d);
    MAXABS (d, a);
    if (a < GEOMETRIC_EPSILON) return ele->nodes [i];
  }

  return -1;
}

/* return distance of a spatial (ref == 0) or referential (ref == 1) point to the element */
double ELEMENT_Point_Distance (MESH *msh, ELEMENT *ele, double *point, int ref)
{
  double nodes [8][3], q [3], d;

  load_nodes (ref ? msh->ref_nodes : msh->cur_nodes, ele->type, ele->nodes, nodes);

  d = gjk_convex_point ((double*)nodes, ele->type, point, q);

  return d;
}

/* return distance of a spatial point to the element */
double ELEMENT_Spatial_Point_Distance (MESH *msh, ELEMENT *ele, double *point)
{
  double nodes [8][3], q [3], d;

  load_nodes (msh->cur_nodes, ele->type, ele->nodes, nodes);

  d = gjk_convex_point ((double*)nodes, ele->type, point, q);

  return d;
}

/* test wether two elements are adjacent
 * through a common face, edge or vertex */
int ELEMENT_Adjacent (ELEMENT *one, ELEMENT *two)
{
  int n, m;

  for (n = 0; n < one->type; n ++)
    for (m = 0; m < two->type; m ++)
      if (one->nodes [n] == two->nodes [m]) return 1; /* simly check whether a node is shared */

  /* one could copy node lists, sort them and compare in linear time,
   * which would give => n + m (copying) + n log n + m log m  (sorting) + n + m (comparing);
   * for the longest case, n = m = 8, this gives 90 steps, whereas the simple check above takes 64 steps */

  return 0;
}

/* update spatial extents of an individual element */
void ELEMENT_Extents (MESH *msh, ELEMENT *ele, double *extents)
{
  element_extents (msh, ele, 0, extents);
}


/* update referential extents of an individual element */
void ELEMENT_Ref_Extents (MESH *msh, ELEMENT *ele, double *extents)
{
  element_extents (msh, ele, 1, extents);
}

/* copy element vertices into 'ver' and return their count */
int ELEMENT_Vertices (MESH *msh, ELEMENT *ele, double *ver)
{
  typedef double (*node_type) [3];

  load_nodes (msh->cur_nodes, ele->type, ele->nodes, (node_type) ver);
  return ele->type;
}

/* return 6-vector (normal, point) planes of element faces,
 * where 'sur' are code of the surfaces of * first 'k' planes
 * correspond to the surface faces; return the total number of planes */
int ELEMENT_Planes (MESH *msh, ELEMENT *ele, double *pla, int *sur, int *k)
{
  FACE faces [6], *fac;
  int n, m, j, l;

  m = neighs (ele->type); 

  k = (k ? k : &l); /* in case of NULL */

  for (n = 0; n < m; n ++)
    setup_face (ele, n, &faces [n], 0);

  /* copy first 'k' surface planes */
  for ((*k) = 0, fac = ele->faces; fac; (*k) ++, fac = fac->next)
  {
    if (sur) sur [*k] = fac->surface;
    COPY (fac->normal, pla);
    COPY (msh->cur_nodes [fac->nodes[0]], pla + 3);
    pla += 6;
  }

  /* copy the remaining planes */
  for (n = 0, j = (*k); n < m; n ++)
  {
    for (fac = ele->faces; fac; fac = fac->next)
      if (face_compare (&faces [n], fac) == 0) break; /* surface face */

    if (fac == NULL) /* not a surface face */
    {
      setup_normal (msh->cur_nodes, &faces [n]);
      COPY (faces [n].normal, pla);
      COPY (msh->cur_nodes [faces [n].nodes[0]], pla + 3);
      pla += 6;
      j ++;
    }
  }

  return j;
}

/* copy element into a convex */
CONVEX* ELEMENT_Convex (MESH *msh, ELEMENT *ele, int ref)
{
  int fac [30], surfaces [6] = {INT_MAX, INT_MAX,
    INT_MAX, INT_MAX, INT_MAX, INT_MAX}, nfac, *f, n;
  double nodes [8][3];
  CONVEX *cvx;

  load_nodes (ref ? msh->ref_nodes : msh->cur_nodes, ele->type, ele->nodes, nodes);

  nfac = neighs (ele->type); /* number of faces */

  for (n = 0, f = fac; n < nfac; n ++) f = setup_face_vertices (ele, n, f); /* write face vertex indices */
  for (FACE *fac = ele->faces; fac; fac = fac->next) surfaces [fac->index] = fac->surface; /* set surface identifiers */
  cvx = CONVEX_Create (NULL, (double*)nodes, ele->type, fac, nfac, surfaces, ele->volume); /* create convex */

  return cvx;
}

/* compute element volume */
double ELEMENT_Volume (MESH *msh, ELEMENT *ele, int ref)
{
  int fac [30], surfaces [6] = {INT_MAX, INT_MAX,
    INT_MAX, INT_MAX, INT_MAX, INT_MAX}, nfac, *f, n;
  double nodes [8][3], volume;
  CONVEX *cvx;

  load_nodes (ref ? msh->ref_nodes : msh->cur_nodes, ele->type, ele->nodes, nodes);
  nfac = neighs (ele->type); /* number of faces */
  for (n = 0, f = fac; n < nfac; n ++) f = setup_face_vertices (ele, n, f); /* write face vertex indices */
  cvx = CONVEX_Create (NULL, (double*)nodes, ele->type, fac, nfac, surfaces, ele->volume); /* create convex */
  volume = CONVEX_Volume (cvx, ref);
  CONVEX_Destroy (cvx);
  return volume;
}

/* compute partial characteristic of an element: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void ELEMENT_Char_Partial (MESH *msh, ELEMENT *ele, int ref, double *vo, double *sx, double *sy, double *sz, double *eul)
{
  CONVEX *cvx;

  cvx = ELEMENT_Convex (msh, ele, ref);
  CONVEX_Char_Partial (cvx, 0, vo, sx, sy, sz, eul);
  CONVEX_Destroy (cvx);
}

/* pack face */
static void face_pack (FACE *fac, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_doubles (dsize, d, doubles, fac->normal, 3);
  pack_int (isize, i, ints, fac->type);
  pack_ints (isize, i, ints, fac->nodes, fac->type);
  pack_int (isize, i, ints, fac->index);
  pack_int (isize, i, ints, fac->surface);
}

/* unpack face */
static FACE* face_unpack (MESH *msh, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  FACE *fac;

  ERRMEM (fac = MEM_Alloc (&msh->facmem));

  unpack_doubles (dpos, d, doubles, fac->normal, 3);
  fac->type = unpack_int (ipos, i, ints);
  unpack_ints (ipos, i, ints, fac->nodes, fac->type);
  fac->index = unpack_int (ipos, i, ints);
  fac->surface = unpack_int (ipos, i, ints);

  return fac;
}

/* pack element */
static void element_pack (ELEMENT *ele, MAP *map, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  FACE *fac;
  int n;

  pack_int (isize, i, ints, ele->type);
  pack_int (isize, i, ints, ele->neighs);
  pack_int (isize, i, ints, ele->volume);

  pack_ints (isize, i, ints, ele->nodes, ele->type);

  /* rather than adjacency pack indices of neighbours in the output sequence */
  for (n = 0; n < ele->neighs; n ++) pack_int (isize, i, ints, (int) (long) MAP_Find (map, ele->adj [n], NULL));

  pack_int (isize, i, ints, ele->mat ? 1 : 0); /* pack material existence flag */
  if (ele->mat) pack_string (isize, i, ints, ele->mat->label);

  for (n = 0, fac = ele->faces; fac; fac = fac->next) n ++; /* faces count */

  pack_int (isize, i, ints, n);

  for (fac = ele->faces; fac; fac = fac->next)
    face_pack (fac, dsize, d, doubles, isize, i, ints);
}

/* unpack element */
static ELEMENT* element_unpack (void *solfec, MESH *msh, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  FACE *fac, *tail;
  ELEMENT *ele;
  int j, n;

  ERRMEM (ele = MEM_Alloc (&msh->elemem));

  ele->type = unpack_int (ipos, i, ints);
  ele->neighs = unpack_int (ipos, i, ints);
  ele->volume = unpack_int (ipos, i, ints);

  unpack_ints (ipos, i, ints, ele->nodes, ele->type);

  for (n = 0; n < ele->neighs; n ++) ele->adj [n] = (ELEMENT*) (long) unpack_int (ipos, i, ints);

  j = unpack_int (ipos, i, ints); /* unpack material existence flag */

  if (j)
  {
    ASSERT_TEXT (solfec, "Trying to unpack element material without the Solfec pointer");
    SOLFEC *sol = solfec;
    char *label = unpack_string (ipos, i, ints);
    ASSERT_DEBUG_EXT (ele->mat = MATSET_Find (sol->mat, label), "Failed to find material when unpacking an element");
    free (label);
  }

  n = unpack_int (ipos, i, ints); /* faces count */

  for (tail = NULL, j = 0; j < n; j ++)
  {
    fac = face_unpack (msh, dpos, d, doubles, ipos, i, ints);
    fac->ele = ele;
    if (tail) tail->next = fac;
    else ele->faces = fac;
    tail = fac;
    fac->n = msh->faces;
    msh->faces = fac;
  }

  return ele;
}

/* pack mesh */
void MESH_Pack (MESH *msh, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  ELEMENT *ele;
  MEM mem;
  MAP *map;
  int n;

  pack_int (isize, i, ints, msh->nodes_count);
  pack_int (isize, i, ints, msh->bulkeles_count);
  pack_int (isize, i, ints, msh->surfeles_count);

  pack_doubles (dsize, d, doubles, (double*)msh->cur_nodes, msh->nodes_count * 3);
  pack_doubles (dsize, d, doubles, (double*)msh->ref_nodes, msh->nodes_count * 3);

  MEM_Init (&mem, sizeof (MAP), msh->bulkeles_count + msh->surfeles_count);

  for (map = NULL, n = 0, ele = msh->bulkeles; ele; ele = ele->next, n ++)
    MAP_Insert (&mem, &map, ele, (void*) (long) n, NULL);

  for (ele = msh->surfeles; ele; ele = ele->next, n ++)
    MAP_Insert (&mem, &map, ele, (void*) (long) n, NULL);

  for (ele = msh->bulkeles; ele; ele = ele->next)
    element_pack (ele, map, dsize, d, doubles, isize, i, ints);

  for (ele = msh->surfeles; ele; ele = ele->next)
    element_pack (ele, map, dsize, d, doubles, isize, i, ints);

  MEM_Release (&mem);
}

/* unpack mesh */
MESH* MESH_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  ELEMENT *ele, **tab, *tail;
  int n, m, k;
  MESH *msh;

  ERRMEM (msh = MEM_CALLOC (sizeof (MESH)));

  msh->nodes_count = unpack_int (ipos, i, ints);
  msh->bulkeles_count = unpack_int (ipos, i, ints);
  msh->surfeles_count = unpack_int (ipos, i, ints);

  m = msh->bulkeles_count + msh->surfeles_count;

  MEM_Init (&msh->elemem, sizeof (ELEMENT), m);
  MEM_Init (&msh->facmem, sizeof (FACE), MEMCHUNK);
  MEM_Init (&msh->mapmem, sizeof (MAP), MIN (m, MEMCHUNK));

  ERRMEM (msh->ref_nodes = malloc (sizeof (double [3]) * (msh->nodes_count * 2)));
  msh->cur_nodes = msh->ref_nodes + msh->nodes_count;

  unpack_doubles (dpos, d, doubles, (double*)msh->cur_nodes, msh->nodes_count * 3);
  unpack_doubles (dpos, d, doubles, (double*)msh->ref_nodes, msh->nodes_count * 3);

  ERRMEM (tab = malloc (m * sizeof (ELEMENT*)));

  for (msh->bulkeles = tail = NULL, n = m = 0; n < msh->bulkeles_count; n ++, m ++)
  {
    ele = element_unpack (solfec, msh, dpos, d, doubles, ipos, i, ints);

    if (tail) ele->prev = tail, tail->next = ele;
    else msh->bulkeles = ele;
    tail = ele;

    tab [m] = tail = ele;
  }

  for (msh->surfeles = tail = NULL, n = 0; n < msh->surfeles_count; n ++, m ++)
  {
    ele = element_unpack (solfec, msh, dpos, d, doubles, ipos, i, ints);

    if (tail) ele->prev = tail, tail->next = ele;
    else  msh->surfeles = ele;
    tail = ele;

    tab [m] = ele;
  }

  for (n = 0; n < m; n ++) /* map adjacency */
  {
    ele = tab [n];

    for (k = 0; k < ele->neighs; k ++)
      ele->adj [k] = tab [(int) (long) ele->adj [k]];
  }

  free (tab);

  return msh;
}

/* export MBFCP definition */
void MESH_2_MBFCP (MESH *msh, FILE *out)
{
  ELEMENT *ele;
  FACE *fac;
  int n;

  fprintf (out, "NODES:\t%d\n", msh->nodes_count);

  for (n = 0; n < msh->nodes_count; n ++)
  {
    fprintf (out, "%g  %g  %g\n", msh->ref_nodes [n][0], msh->ref_nodes [n][1], msh->ref_nodes [n][2]);
  }

  fprintf (out, "ELEMENTS:\t%d\n", msh->surfeles_count + msh->bulkeles_count);

  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    fprintf (out, "%d", ele->type);
    for (n = 0; n < ele->type; n ++) fprintf (out, "  %d", ele->nodes [n]);
    fprintf (out, "\n");
  }

  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    fprintf (out, "%d", ele->type);
    for (n = 0; n < ele->type; n ++) fprintf (out, "  %d", ele->nodes [n]);
    fprintf (out, "\n");
  }

  for (fac = msh->faces, n = 0; fac; fac = fac->n, n ++);

  fprintf (out, "FACES:\t%d\n", n);

  for (fac = msh->faces; fac; fac = fac->n)
  {
    fprintf (out, "%d", fac->type);
    for (n = 0; n < fac->type; n ++) fprintf (out, "  %d", fac->nodes [n]);
    fprintf (out, "  %d\n", fac->surface);
  }
}

/* write mesh */
void MESH_Write (MESH *msh, char *path)
{
  int dsize, doubles, isize, ints, size;
  int *i, *data;
  double *d;
  FILE *f;
  XDR x;

  ASSERT (f = fopen (path, "w"), ERR_FILE_OPEN);
  xdrstdio_create (&x, f, XDR_ENCODE);

  dsize = doubles = isize = ints = 0;
  d = NULL;
  i = NULL;

  MESH_Pack (msh, &dsize, &d, &doubles, &isize, &i, &ints);
  data = compress (CMP_FASTLZ, d, doubles, i, ints, &size);

  ASSERT (xdr_int (&x, &size), ERR_PBF_WRITE);
  ASSERT (xdr_vector (&x, (char*)data, size, sizeof (int), (xdrproc_t)xdr_int), ERR_PBF_WRITE);

  free (data);
  free (d);
  free (i);
  xdr_destroy (&x);
  fclose (f);
}

/* read mesh */
MESH* MESH_Read (char *path)
{
  int dpos, doubles, ipos, ints, size;
  int *i, *data;
  double *d;
  MESH *msh;
  FILE *f;
  XDR x;

  if (!(f = fopen (path, "r"))) return NULL;
  xdrstdio_create (&x, f, XDR_DECODE);

  dpos = doubles = ipos = ints = 0;
  d = NULL;
  i = NULL;

  ASSERT (xdr_int (&x, &size), ERR_PBF_READ);
  ERRMEM (data = malloc (size * sizeof (int)));
  ASSERT (xdr_vector (&x, (char*)data, size, sizeof (int), (xdrproc_t)xdr_int), ERR_PBF_READ);

  decompress (data, size, &d, &doubles, &i, &ints);
  msh = MESH_Unpack (NULL, &dpos, d, doubles, &ipos, i, ints);

  free (data);
  free (d);
  free (i);
  xdr_destroy (&x);
  fclose (f);

  return msh;
}
