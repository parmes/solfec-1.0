/*
 * msh.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * mesh definition
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

#include "mat.h"
#include "mem.h"
#include "cvx.h"
#include "mot.h"
#include "tri.h"

#ifndef ELEMENT_TYPE
#define ELEMENT_TYPE
typedef struct element ELEMENT;
#endif

#ifndef __msh__
#define __msh__

typedef struct face FACE;
typedef struct mesh MESH;

/* triangular or quadrilateral face */
struct face
{
  double normal [3];

  int type, /* 3, 4 => triangle, quadrilateral */
      nodes [4],
      index, /* index within the element */
      surface; /* surface identifier */

  ELEMENT *ele;

  FACE *next;
};

/* finite element */
struct element
{
  int type, /* 4, 5, 6, 8 => tetrahedron, pyramid, wedge, hexahedron */
      nodes [8], /* node indices */
      neighs, /* number of neighbours */
      domnum, /* number of integration domains */
      volume; /* volume identifier */

  BULK_MATERIAL *mat;

  TRISURF *dom; /* integration domains */

  ELEMENT *prev,
	  *next,
	  *adj [6]; /* neighbouring elements */

  FACE *faces; /* corresponding surface faces */
};

/* general mesh */
struct mesh
{
  MEM facmem,
      elemem;

  double (*ref_nodes) [3],
	 (*cur_nodes) [3];

  ELEMENT *surfeles,
	  *bulkeles;
  
  int  surfeles_count,
       bulkeles_count,
       nodes_count;
};

/* create mesh from vector of nodes, element list in format =>
 * {nuber of nodes, node0, node1, ..., volume1}, {REPEAT}, ..., 0 (end of list); and surface kinds in format =>
 * global surface, {number of nodes, node0, node1, ..., surface}, {REPEAT}, ..., 0 (end of list); */
MESH* MESH_Create (double (*nodes) [3], int *elements, int *surfaces);

/* create a meshed hexahedron by specifying its eight nodes and
 * division numbers along three edges adjacent to the 1st node */
MESH* MESH_Hex (double (*nodes) [3], int i, int j, int k, int *surfaces, /* six surfaces are given */
                int volume, double *dx, double *dy, double *dz); /* spacing of divisions */

/* dummy adjacency update (needed in shp.c) */
void MESH_Update_Adjacency (MESH *msh);

/* create a copy of a mesh */
MESH* MESH_Copy (MESH *msh);

/* scaling of a mesh */
void MESH_Scale (MESH *msh, double *vector);

/* translation of a mesh */
void MESH_Translate (MESH *msh, double *vector);

/* rotation of a mesh */
void MESH_Rotate (MESH *msh, double *point, double *vector, double angle);

/* compute current partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void MESH_Char_Partial (MESH *msh, double *vo, double *sx, double *sy, double *sz, double *eul);

/* get 'cur' characteristics of the meshed shape:
 * volume, mass center, and Euler tensor (centered) */
void MESH_Char (MESH *msh, double *volume, double *center, double *euler);

/* find an element containing a spatial or referential point */
ELEMENT* MESH_Element_Containing_Point (MESH *msh, double *point, int ref);

/* update mesh according to the given motion */
void MESH_Update (MESH *msh, void *body, void *shp, MOTION motion);

/* convert mesh into a list of convices;
 * surfonly > 0 => use only surface elements */
CONVEX* MESH_Convex (MESH *msh, int surfonly);

/* compute extents of entire mesh */
void MESH_Extents (MESH *msh, double *extents);

/* return first not NULL bulk material of an element */
void* MESH_First_Bulk_Material (MESH *msh);

/* free mesh memory */
void MESH_Destroy (MESH *msh);
  
/* does the element contain the point? */
int ELEMENT_Contains_Point (MESH *msh, ELEMENT *ele, double *point, int ref);

/* test wether two elements are adjacent
 * through a common face, edge or vertex */
int ELEMENT_Adjacent (ELEMENT *one, ELEMENT *two);

/* update extents of an individual element */
void ELEMENT_Extents (MESH *msh, ELEMENT *ele, double *extents);

/* copy element vertices into 'ver' and return their count */
int ELEMENT_Vertices (MESH *msh, ELEMENT *ele, double *ver);

/* return 6-vector (normal, point) planes of element faces, where 'sur'
 * are code of the surfaces of * first 'k' planes correspond to the
 * surface faces (NULL accepted); return the total number of planes */
int ELEMENT_Planes (MESH *msh, ELEMENT *ele, double *pla, int *sur, int *k);

/* copy element into a convex */
CONVEX* ELEMENT_Convex (MESH *msh, ELEMENT *ele);

/* pack mesh into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void MESH_Pack (MESH *msh, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack mesh from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
MESH* MESH_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

#endif
