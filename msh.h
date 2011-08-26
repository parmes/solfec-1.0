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
#include "map.h"
#include "set.h"

#ifndef ELEMENT_TYPE
#define ELEMENT_TYPE
typedef struct element ELEMENT;
#endif

#ifndef __msh__
#define __msh__

typedef struct node NODE;
typedef struct face FACE;
typedef struct mesh MESH;

/* surface node */
struct node
{
  double *cur; /* points to the current node */

  FACE **fac; /* adjacent faces */

  int nfac; /* number of them */

  NODE *n; /* list */
};

/* triangular or quadrilateral face */
struct face
{
  double normal [6]; /* spatial, referential */

  int type, /* 3, 4 => triangle, quadrilateral */
      nodes [4], /* node numbers */
      index, /* index within the element */
      surface; /* surface identifier */

  NODE *nod [4];

  ELEMENT *ele;

  FACE *next, *n; /* element, mesh list */
};

/* finite element */
struct element
{
  int type, /* 4, 5, 6, 8 => tetrahedron, pyramid, wedge, hexahedron */
      nodes [8], /* node numbers */
      neighs, /* number of neighbours */
      domnum, /* number of integration domains; or destination partition after MESH_Partition */
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
  MEM nodmem,
      facmem,
      elemem,
      mapmem;

  double (*ref_nodes) [3],
	 (*cur_nodes) [3];

  NODE *surfnodes;

  ELEMENT *surfeles,
	  *bulkeles;
  
  FACE *faces;

  int  surfeles_count,
       bulkeles_count,
       nodes_count,
       surfnodes_count;

  MAP *map; /* MESH_Element_With_Node uses it */
};

/* create mesh from vector of nodes, element list in format =>
 * {nuber of nodes, node0, node1, ..., volume1}, {REPEAT}, ..., 0 (end of list); and surface kinds in format =>
 * global surface, {number of nodes, node0, node1, ..., surface}, {REPEAT}, ..., 0 (end of list); */
MESH* MESH_Create (double (*nodes) [3], int *elements, int *surfaces);

/* create a meshed hexahedron by specifying its eight nodes and
 * division numbers along three edges adjacent to the 1st node */
MESH* MESH_Hex (double (*nodes) [3], int i, int j, int k, int *surfaces, /* six surfaces are given */
                int volume, double *dx, double *dy, double *dz); /* spacing of divisions */

/* create pipe like mesh using a point, a direction, an inner radius,
 * a thickness and subdivison counts along the direction, radius and thickness */
MESH* MESH_Pipe (double *pnt, double *dir, double rin, double thi,
                 int ndir, int nrad, int nthi, int *surfaces, int volume); /* surfaces: bottom, top, inner, outer */

/* dummy adjacency update (needed in shp.c) */
void MESH_Update_Adjacency (MESH *msh);

/* dummy local plane adjacency breaking routine */
int MESH_Break_Adjacency (MESH *msh, double *point, double *normal);

/* create a copy of a mesh */
MESH* MESH_Copy (MESH *msh);

/* scaling of a mesh */
void MESH_Scale (MESH *msh, double *vector);

/* translation of a mesh */
void MESH_Translate (MESH *msh, double *vector);

/* rotation of a mesh */
void MESH_Rotate (MESH *msh, double *point, double *vector, double angle);

/* cut through mesh with a plane; return triangulated cross-section; vertices in the triangles
 * point to the memory allocated after the triangles memory; adjacency is not maintained;
 * TRI->ptr stores a pointer to the geometrical object that has been cut by the triangle */
TRI* MESH_Cut (MESH *msh, double *point, double *normal, int *m);

/* as above but this time the plane and the cut are in the reference configuration */
TRI* MESH_Ref_Cut (MESH *msh, double *point, double *normal, int *m);

/* split mesh in two with plane defined by (point, normal); output meshes are tetrahedral if some
 * elements are crossed; if only element boundaries are crossed then the original mesh is used;
 * topoadj != 0 implies cutting from the point and through the topological adjacency only */
void MESH_Split (MESH *msh, double *point, double *normal, short topoadj, int surfid, MESH **one, MESH **two);

/* compute partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void MESH_Char_Partial (MESH *msh, int ref, double *vo, double *sx, double *sy, double *sz, double *eul);

/* get characteristics of the meshed shape:
 * volume, mass center, and Euler tensor (centered) */
void MESH_Char (MESH *msh, int ref, double *volume, double *center, double *euler);

/* find an element containing a spatial or referential point */
ELEMENT* MESH_Element_Containing_Point (MESH *msh, double *point, int ref);

/* find an element containing a spatial point */
ELEMENT* MESH_Element_Containing_Spatial_Point (MESH *msh, double *point);

/* find an element with a given node */
ELEMENT* MESH_Element_With_Node (MESH *msh, int node);

/* collect elements around a node (ele->node [i] == node && *set == NULL initially assumed) */
void MESH_Elements_Around_Node (ELEMENT *ele, int node, SET **set);

/* update mesh according to the given motion */
void MESH_Update (MESH *msh, void *body, void *shp, MOTION motion);

/* convert mesh into a list of convices;
 * ref > 0 => create referential mesh image;
 * otherwise => create current mesh image;
 * CONVEX->ele[0] == corresponding element */
CONVEX* MESH_Convex (MESH *msh, int ref);

/* compute extents of entire mesh */
void MESH_Extents (MESH *msh, double *extents);

/* compute oriented extents of entire mesh */
void MESH_Oriented_Extents (MESH *msh, double *vx, double *vy, double *vz, double *extents);

/* return first not NULL bulk material of an element */
void* MESH_First_Bulk_Material (MESH *msh);

/* partition mesh; return the resultant mesh parts; output a table of tuples (m1, m2, n1, n2) of gluing nodes,
 * where m1, m2 is a pair of the output meshes and n1, n2 are their corresponding coincident nodes;
 * additionally output tuples (m1, m2, e1, e2) of topologically adjacent surface element
 * pairs from the partitions boundaries (indexed as stored in lists and outputed by SGP_Create);
 * upon exit the 'domnum' element values of the input mesh indicate destination partitions of the elements */
MESH** MESH_Partition (MESH *msh, int nparts, int *numglue, int **gluenodes, int *numadj, int **adjeles);

/* free mesh memory */
void MESH_Destroy (MESH *msh);
  
/* does the element contain a spatial/referential point? */
int ELEMENT_Contains_Point (MESH *msh, ELEMENT *ele, double *point, int ret);

/* does the element contain a spatial point? */
int ELEMENT_Contains_Spatial_Point (MESH *msh, ELEMENT *ele, double *point);

/* return >= node index if point == node[index] or -1 otherwise */
int ELEMENT_Ref_Point_To_Node (MESH *msh, ELEMENT *ele, double *point);

/* return distance of a spatial (ref == 0) or referential (ref == 1) point to the element */
double ELEMENT_Point_Distance (MESH *msh, ELEMENT *ele, double *point, int ref);

/* return distance of a spatial point to the element */
double ELEMENT_Spatial_Point_Distance (MESH *msh, ELEMENT *ele, double *point);

/* test wether two elements are adjacent
 * through a common face, edge or vertex */
int ELEMENT_Adjacent (ELEMENT *one, ELEMENT *two);

/* update spatial extents of an individual element */
void ELEMENT_Extents (MESH *msh, ELEMENT *ele, double *extents);

/* update referential extents of an individual element */
void ELEMENT_Ref_Extents (MESH *msh, ELEMENT *ele, double *extents);

/* copy element vertices into 'ver' and return their count */
int ELEMENT_Vertices (MESH *msh, ELEMENT *ele, double *ver);

/* return 6-vector (normal, point) planes of element faces, where 'sur'
 * are code of the surfaces of * first 'k' planes correspond to the
 * surface faces (NULL accepted); return the total number of planes */
int ELEMENT_Planes (MESH *msh, ELEMENT *ele, double *pla, int *sur, int *k);

/* copy element into a convex */
CONVEX* ELEMENT_Convex (MESH *msh, ELEMENT *ele, int ref);

/* compute element volume */
double ELEMENT_Volume (MESH *msh, ELEMENT *ele, int ref);

/* compute partial characteristic of an element: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void ELEMENT_Char_Partial (MESH *msh, ELEMENT *ele, int ref, double *vo, double *sx, double *sy, double *sz, double *eul);

/* update spatial extents of an individual node */
void NODE_Extents (MESH *msh, NODE *nod, double *extents);

/* pack mesh into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void MESH_Pack (MESH *msh, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack mesh from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
MESH* MESH_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

/* export MBFCP definition */
void MESH_2_MBFCP (MESH *msh, FILE *out);

/* write mesh */
void MESH_Write (MESH *msh, char *path);

/* read mesh */
MESH* MESH_Read (char *path);

#endif
