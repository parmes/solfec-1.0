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
#include "sol.h"
#include "spx.h"
#include "err.h"
#include "mem.h"
#include "map.h"
#include "alg.h"
#include "msh.h"
#include "pck.h"

/* used in some pools */
#define MEMCHUNK 1024

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
  n2 = fac->nodes [3];
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

      fac->ele = NULL; /* mark as the inner face */
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
int element_has_edge (ELEMENT *ele, int nod1, int nod2)
{
  int *nodes = ele->nodes,
      type = ele->type, n, j;

  for (n = j = 0; n < type; n ++)
    if (nodes [n] == nod1 ||
	nodes [n] == nod2) j ++;

  return (j == 2 ? 1 : 0);
}

/* compute planes - CCW oriented faces */
inline static void computeplanes (CONVEX *cvx)
{
  double *a, *b, *c, *nl;
  int n, m;

  for (n = m = 0; n < cvx->nfac; n ++, m += (cvx->fac [m] + 1))
  {
    a = &cvx->cur [cvx->fac [m + 1]];
    b = &cvx->cur [cvx->fac [m + 2]];
    c = &cvx->cur [cvx->fac [m + 3]];
    nl =  &cvx->pla [n * 4];
    NORMAL (a, b, c, nl);
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
  ERRMEM (msh = malloc (sizeof (MESH)));
  elemem = &msh->elemem;
 
  /* calculate elements */ 
  for (eleptr = elements; eleptr [0]; eleptr += (eleptr [0]+2)) elements_count ++;

  MEM_Init (elemem, sizeof (ELEMENT), elements_count);
  MEM_Init (&facmem, sizeof (FACE), MEMCHUNK);
  MEM_Init (&mapmem, sizeof (MAP), MEMCHUNK);

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
    if (fac->ele)
    {
      ele = fac->ele;

      cac = MEM_Alloc (&msh->facmem);
      setup_face (ele, fac->index, cac, 0); /* setup face nodes without sorting them */
      cac->index = fac->index;
      cac->ele = fac->ele;
      setup_normal (msh->cur_nodes, cac); /* calculate outer normal */
      cac->next = ele->faces; /* append element face list */
      ele->faces = cac;

      /* set the mapped surface kind if possible => otherwise the global one */
      gac = MAP_Find (smap, fac, (MAP_Compare) face_compare); 
      cac->surface = (gac ? gac->surface : surfaces [0]);
    }
  }

  /* clean up */
  MEM_Release (&facmem);
  MEM_Release (&mapmem);

  /* refine geometric epsilon */
  GEOMETRIC_EPSILON_ADAPT ((double*)msh->cur_nodes, msh->nodes_count);

  return msh;
}

/* create a meshed hexahedron by specifying its eight nodes and
 * division numbers along three edges adjacent to the 1st node */
MESH* MESH_Hex (double (*nodes) [3], int i, int j, int k, int *surfaces, int volume, double *dx, double *dy, double *dz)
{
  double (*nod) [3],
	 x, y, z,
	 ddx, ddy, ddz;
  int *ele, *ee, *ss,
      *sur, n, m,
      mx, my, mz,
      ii, jj, kk, nn;
  MESH *msh;


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
    if (ii == 0) /* face 2 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[1][m]];
      ss [5] = surfaces [1];
      ss += 6;
    }
    if (jj == (j-1)) /* face 3 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[2][m]];
      ss [5] = surfaces [2];
      ss += 6;
    }
    if (ii == (i-1)) /* face 4 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[3][m]];
      ss [5] = surfaces [3];
      ss += 6;
    }
    if (jj == 0) /* face 5 */
    {
      ss [0] = 4;
      for (m = 0; m < 4; m ++) ss [m+1] = ee[hex[4][m]];
      ss [5] = surfaces [4];
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

/* dummy adjacency update (needed in shp.c) */
void MESH_Update_Adjacency (MESH *msh)
{
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
  ERRMEM (ret = malloc (sizeof (MESH)));
  MEM_Init (&ret->elemem, sizeof (ELEMENT), msh->surfeles_count + msh->bulkeles_count);
  ret->surfeles_count = msh->surfeles_count;
  ret->bulkeles_count = msh->bulkeles_count;
  for (n = 0, ele = msh->surfeles; ele; ele = ele->next)
    for (fac = ele->faces; fac; fac = fac->next) n ++; /* count surface faces */
  MEM_Init (&ret->facmem, sizeof (FACE), n);
  ERRMEM (ret->ref_nodes = malloc (sizeof (double [3]) * (msh->nodes_count * 2)));
  ret->cur_nodes = ret->ref_nodes + msh->nodes_count;
  ret->nodes_count = msh->nodes_count;
  ret->surfeles = ret->bulkeles = NULL;
  
  /* copy vertices */
  memcpy (ret->ref_nodes, msh->ref_nodes, sizeof (double [3]) * (msh->nodes_count * 2));

  /* create, copy and map elements */
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    ERRMEM (cpy = MEM_Alloc (&ret->elemem));
    *cpy = *ele;

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

    /* maintain list */
    cpy->prev = NULL;
    cpy->next = ret->surfeles;
    if (ret->bulkeles) ret->bulkeles->prev = cpy;
    ret->bulkeles = cpy;
   
    /* map new and old elements */
    MAP_Insert (&mapmem, &map, ele, cpy, NULL);
  }

  /* maintain adjacency => use the element map */
  for (ele = ret->surfeles; ele; ele = ele->next)
  {
    for (n = 0; n < ele->neighs; n ++)
      ele->adj [n] = MAP_Find (map, ele->adj [n], NULL); /* find a new pointer corresponding to the old one */
  }
  for (ele = ret->bulkeles; ele; ele = ele->next)
  {
    for (n = 0; n < ele->neighs; n ++)
      ele->adj [n] = MAP_Find (map, ele->adj [n], NULL);
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
      setup_normal (cur, fac);
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
      setup_normal (cur, fac);
}

/* compute current partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void MESH_Char_Partial (MESH *msh, double *vo, double *sx, double *sy, double *sz, double *eul)
{
  double zero [3] = {0, 0, 0},
	 J, (*cur) [3] = msh->cur_nodes,
	 a [3], b [3], c [3];
  ELEMENT *ele;
  FACE *fac;

  /* loop over the surface faces and use simplex integration
   * in order to calculate the volume characteristics */
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    for (fac = ele->faces; fac; fac = fac->next)
    {
      COPY (cur [fac->nodes [0]], a);
      COPY (cur [fac->nodes [1]], b);
      COPY (cur [fac->nodes [2]], c);

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

      if (fac->type == 4)
      {
	COPY (c, b);
	COPY (cur [fac->nodes [3]], c);

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

/* get 'cur' characteristics of the meshed shape:
 * volume, mass center, and Euler tensor (centered) */
void MESH_Char (MESH *msh, double *volume, double *center, double *euler)
{
  double vo, sx, sy, sz,
	 cen [3], eul [9];

  vo = sx = sy = sz = 0.0;
  SET9 (eul, 0.0);

  MESH_Char_Partial (msh, &vo, &sx, &sy, &sz, eul);

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

/* find an element containing the point */
ELEMENT* MESH_Element_Containing_Point (MESH *msh, double *point, int ref)
{
  /* the linear search should be at some point optimised with a proper
     spatial search => move point to the reference configuration and
     query a search tree based on the reference configuration geometry */

  double (*nodes) [3] = ref ? msh->ref_nodes : msh->cur_nodes;
  double pla [24];
  ELEMENT *ele;

  /* first search surface elements */
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    create_element_planes (nodes, ele, pla);
    if (point_inside (neighs (ele->type), pla, point)) return ele;
  }

  /* then the bulk elements */
  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    create_element_planes (nodes, ele, pla);
    if (point_inside (neighs (ele->type), pla, point)) return ele;
  }

  return NULL;
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

  for (n = 0; n < m; n ++)
    motion (body, shp, NULL, ref [n], cur [n]); /* move current nodes (NULL for gobj implies nodal update) */

  for (ele = msh->surfeles; ele; ele = ele->next)
    for (fac = ele->faces; fac; fac = fac->next)
      setup_normal (cur, fac); /* update normals */
}

/* convert mesh into a list of convices;
 * surfonly > 0 => use only surface elements */
CONVEX* MESH_Convex (MESH *msh, int surfonly)
{
  int fac [30], surfaces [6] = {INT_MAX, INT_MAX,
    INT_MAX, INT_MAX, INT_MAX, INT_MAX}, nfac, *f, n;
  double nodes [8][3];
  ELEMENT *ele;
  CONVEX *cvx;

  cvx = NULL;
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    load_nodes (msh->cur_nodes, ele->type, ele->nodes, nodes);
    nfac = neighs (ele->type); /* number of faces */
    for (n = 0, f = fac; n < nfac; n ++) f = setup_face_vertices (ele, n, f); /* write face vertex indices */
    for (FACE *fac = ele->faces; fac; fac = fac->next) surfaces [fac->index] = fac->surface; /* set surface identifiers */
    cvx = CONVEX_Create (cvx, (double*)nodes, ele->type, fac, nfac, surfaces, ele->volume); /* add new convex to the list */
  }
  if (!surfonly) /* include bulk elements */
  {
    for (ele = msh->bulkeles; ele; ele = ele->next)
    {
      load_nodes (msh->cur_nodes, ele->type, ele->nodes, nodes);
      nfac = neighs (ele->type); /* number of faces */
      for (n = 0, f = fac; n < nfac; n ++) f = setup_face_vertices (ele, n, f); /* write face vertex indices */
      for (FACE *fac = ele->faces; fac; fac = fac->next) surfaces [fac->index] = fac->surface; /* set surface identifiers */
      cvx = CONVEX_Create (cvx, (double*)nodes, ele->type, fac, nfac, surfaces, ele->volume); /* add new convex to the list */
    }
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
    ELEMENT_Extents (msh, ele, e);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [3] > extents [3]) extents [3] = e [3];
    if (e [4] > extents [4]) extents [4] = e [4];
    if (e [5] > extents [5]) extents [5] = e [5];
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

/* free mesh memory */
void MESH_Destroy (MESH *msh)
{
  ELEMENT *ele;
  int n;

  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    for (n = 0; n < ele->domnum; n ++) free (ele->dom [n].tri);
    free (ele->dom);
  }

  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    for (n = 0; n < ele->domnum; n ++) free (ele->dom [n].tri);
    free (ele->dom);
  }

  MEM_Release (&msh->facmem);
  MEM_Release (&msh->elemem);
  free (msh->ref_nodes);
  free (msh);
}

/* does the element contain a spatial point? */
int ELEMENT_Contains_Point (MESH *msh, ELEMENT *ele, double *point)
{
  double pla [24];

  create_element_planes (msh->cur_nodes, ele, pla);
  return point_inside (neighs (ele->type), pla, point);
}

/* does the element contain a referential point? */
int ELEMENT_Contains_Ref_Point (MESH *msh, ELEMENT *ele, double *point)
{
  double pla [24];

  create_element_planes (msh->ref_nodes, ele, pla);
  return point_inside (neighs (ele->type), pla, point);
}


/* return distance of a spatial (ref == 0) or referential (ref == 1) point to the element */
double ELEMENT_Point_Distance (MESH *msh, ELEMENT *ele, double *point, int ref)
{
  double pla [24];

  create_element_planes (ref ? msh->ref_nodes : msh->cur_nodes, ele, pla);
  return point_distance (neighs (ele->type), pla, point);
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

/* update extents of an individual element */
void ELEMENT_Extents (MESH *msh, ELEMENT *ele, double *extents)
{
  double nodes [8][3];
  int n;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;
  
  load_nodes (msh->cur_nodes, ele->type, ele->nodes, nodes);

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
CONVEX* ELEMENT_Convex (MESH *msh, ELEMENT *ele)
{
  int fac [30], surfaces [6] = {INT_MAX, INT_MAX,
    INT_MAX, INT_MAX, INT_MAX, INT_MAX}, nfac, *f, n;
  double nodes [8][3];
  CONVEX *cvx;

  load_nodes (msh->cur_nodes, ele->type, ele->nodes, nodes);

  nfac = neighs (ele->type); /* number of faces */

  for (n = 0, f = fac; n < nfac; n ++) f = setup_face_vertices (ele, n, f); /* write face vertex indices */
  for (FACE *fac = ele->faces; fac; fac = fac->next) surfaces [fac->index] = fac->surface; /* set surface identifiers */
  cvx = CONVEX_Create (NULL, (double*)nodes, ele->type, fac, nfac, surfaces, ele->volume); /* add new convex to the list */

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
  cvx = CONVEX_Create (NULL, (double*)nodes, ele->type, fac, nfac, surfaces, ele->volume); /* add new convex to the list */
  volume = CONVEX_Volume (cvx, ref);
  CONVEX_Destroy (cvx);
  return volume;
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
  ELEMENT *ele;
  FACE *fac;
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
    SOLFEC *sol = solfec;
    char *label = unpack_string (ipos, i, ints);
    ASSERT_DEBUG_EXT (ele->mat = MATSET_Find (sol->mat, label), "Failed to find material when unpacking an element");
    free (label);
  }

  n = unpack_int (ipos, i, ints); /* faces count */

  for (ele->faces = NULL, j = 0; j < n; j ++)
  {
    fac = face_unpack (msh, dpos, d, doubles, ipos, i, ints);
    fac->ele = ele;
    fac->next = ele->faces;
    ele->faces = fac;
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
  ELEMENT *ele, **tab;
  int n, m, k;
  MESH *msh;

  ERRMEM (msh = malloc (sizeof (MESH)));

  msh->nodes_count = unpack_int (ipos, i, ints);
  msh->bulkeles_count = unpack_int (ipos, i, ints);
  msh->surfeles_count = unpack_int (ipos, i, ints);

  m = msh->bulkeles_count + msh->surfeles_count;

  MEM_Init (&msh->elemem, sizeof (ELEMENT), m);
  MEM_Init (&msh->facmem, sizeof (FACE), MEMCHUNK);

  ERRMEM (msh->ref_nodes = malloc (sizeof (double [3]) * (msh->nodes_count * 2)));
  msh->cur_nodes = msh->ref_nodes + msh->nodes_count;

  unpack_doubles (dpos, d, doubles, (double*)msh->cur_nodes, msh->nodes_count * 3);
  unpack_doubles (dpos, d, doubles, (double*)msh->ref_nodes, msh->nodes_count * 3);

  ERRMEM (tab = malloc (m * sizeof (ELEMENT*)));

  for (msh->bulkeles = NULL, n = m = 0; n < msh->bulkeles_count; n ++, m ++)
  {
    ele = element_unpack (solfec, msh, dpos, d, doubles, ipos, i, ints);
    ele->next = msh->bulkeles;
    if (msh->bulkeles)
      msh->bulkeles->prev = ele;
    msh->bulkeles = ele;
    tab [m] = ele;
  }

  for (msh->surfeles = NULL, n = 0; n < msh->surfeles_count; n ++, m ++)
  {
    ele = element_unpack (solfec, msh, dpos, d, doubles, ipos, i, ints);
    ele->next = msh->surfeles;
    if (msh->surfeles)
      msh->surfeles->prev = ele;
    msh->surfeles = ele;
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
