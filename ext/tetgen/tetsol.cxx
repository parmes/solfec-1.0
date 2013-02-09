/*
 * tetsol.cxx
 * Copyright (C) 2011, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * Tetgen to C interface
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
#include "tetgen.h"
#include "tetsol.h"

extern "C"
{
#include "../../mem.h"
#include "../../err.h"
#include "../../alg.h"
#include "../../kdt.h"
#include "../../spx.h"

/* generate tetrahedrons based on an input mesh object; pass -INT_MAX for (vol/surf)ids to inherit from the mesh */
MESH* tetrahedralize1 (MESH *shape, double volume, double quality, int volid, int surfid)
{
  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  double *a, *b, extents [6];
  int *elements, *surfaces;
  int i, j, n, *ele, *tet;
  char params [512];
  MESH *ret = NULL;
  MAP *map, *item;
  MEM mapmem;
  FACE *fac;

  /* map memory */
  MEM_Init (&mapmem, sizeof (MAP), 128);

  /* calculate faces and map face vertices */
  for (fac = shape->faces, n = j = 0, map = NULL; fac; fac = fac->n, n ++)
  {
    for (i = 0; i < fac->type; i ++)
    {
      if (!MAP_Find_Node (map, (void*) (long) fac->nodes [i], NULL))
      {
	MAP_Insert (&mapmem, &map, (void*) (long) fac->nodes [i], (void*) (long) j, NULL);
	j ++;
      }
    }
  }

  /* 0-based indexing */
  in.firstnumber = 0;

  /* input vertices */
  in.numberofpoints = MAP_Size (map);
  in.pointlist = new REAL[in.numberofpoints * 3];
  for (item = MAP_First (map); item; item = MAP_Next (item))
  {
    i = (int) (long) item->key;
    j = (int) (long) item->data;
    a = shape->cur_nodes [i];
    b = &in.pointlist [3*j];

    COPY (a, b);
  }

  /* input faces */
  in.numberoffacets = n;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  for (fac = shape->faces, n = 0; fac; fac = fac->n, n ++)
  {
    f = &in.facetlist [n];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = fac->type;
    p->vertexlist = new int[p->numberofvertices];
    for (i = 0; i < fac->type; i ++)
    {
      item = MAP_Find_Node (map, (void*) (long) fac->nodes [i], NULL);
      ASSERT_DEBUG (item, "Inconsistent face vertex mapping");
      j = (int) (long) item->data;
      p->vertexlist [i] = j;
    }
    in.facetmarkerlist [n] = fac->surface;
  }

  /* set up parameters */
  if (volume > 0.0 && quality > 0.0) sprintf (params, "Qpa%gq%g", volume, quality);
  else if (volume > 0.0) sprintf (params, "Qpa%g", volume);
  else if (quality > 0.0) sprintf (params, "Qpq%g", quality);
  else sprintf (params, "Qp");

  /* generate mesh */
  tetrahedralize (params, &in, &out);

  /* read triangles */
  ERRMEM (surfaces = (int*)malloc (out.numberoftrifaces * sizeof (int [5]) + sizeof (int [2])));
  if (surfid > -INT_MAX) surfaces [0] = surfid;
  else surfaces [0] = shape->faces->surface;

  for (i = 0, ele = surfaces+1, tet = out.trifacelist;
       i < out.numberoftrifaces; i ++, ele += 5, tet += 3)
  {
    ele [0] = 3;
    ele [1] = tet [0];
    ele [2] = tet [1];
    ele [3] = tet [2];
    ele [4] = out.trifacemarkerlist [i];
  }
  ele [0] = 0;

  /* creake input mesh based kd-tree for volume identifiers mapping */
  KDT *kd = KDT_Create (shape->nodes_count, (double*)shape->cur_nodes, 0.0);
  for (ELEMENT *ele = shape->surfeles; ele; ele = ele->next)
  { ELEMENT_Extents (shape, ele, extents); KDT_Drop (kd, extents, ele); }
  for (ELEMENT *ele = shape->bulkeles; ele; ele = ele->next)
  { ELEMENT_Extents (shape, ele, extents); KDT_Drop (kd, extents, ele); }

  /* read tetrahedrons */
  ERRMEM (elements =  (int*)malloc ((out.numberoftetrahedra+1) * sizeof (int [6])));
  if (out.numberofpoints == 0 || out.numberoftetrahedra == 0) goto err;
  for (i = 0, ele = elements, tet = out.tetrahedronlist;
       i < out.numberoftetrahedra; i ++, ele += 6, tet += 4)
  {
    ele [0] = 4;
    ele [1] = tet [0];
    ele [2] = tet [1];
    ele [3] = tet [2];
    ele [4] = tet [3];

    double *v [4] = {&out.pointlist [tet[0]*3],
                     &out.pointlist [tet[1]*3],
                     &out.pointlist [tet[2]*3],
                     &out.pointlist [tet[3]*3]};

    double p [3] = {.25*(v[0][0]+v[1][0]+v[2][0]+v[3][0]),
                    .25*(v[0][1]+v[1][1]+v[2][1]+v[3][1]),
                    .25*(v[0][2]+v[1][2]+v[2][2]+v[3][2])};

    KDT *leaf = KDT_Pick (kd, p);
    ASSERT_TEXT (leaf,  "Kd-tree based volid mapping has failed: please report this bug!");
    ELEMENT **ptr = (ELEMENT**) leaf->data, **end = ptr + leaf->n;
    for (; ptr != end; ptr ++)
      if (ELEMENT_Contains_Point (shape, *ptr, p, 0)) break;
    ASSERT_TEXT (ptr != end, "Element containing a spatial point has not been found: please report this bug!");
    ele [5] = (*ptr)->volume;
  }
  ele [0] = 0;

  /* create output mesh */
  ret = MESH_Create ((double (*)[3])out.pointlist, elements, surfaces);

  /* clean up */
err:
  free (elements);
  free (surfaces);
  KDT_Destroy (kd);
  MEM_Release (&mapmem);
  
  return ret;
}

/* generate tetrahedrons based on an input file; pass -INT_MAX for (vol/surf)ids to inherit from the input */
MESH* tetrahedralize2 (char *path, double volume, double quality, int volid, int surfid)
{
  int argc = 2;
  char *argv [2] = {"-p", path};
  tetgenbehavior b;
  tetgenio in, out;
  int *elements, *surfaces;
  int i, j, n, *ele, *tet;
  char params [512];
  MESH *ret = NULL;
  FACE *fac;

  /* command line */
  if (!b.parse_commandline(argc, argv)) return NULL;

  /* read from file */
  if (!in.load_plc(b.infilename, (int) b.object)) return NULL;

  /* set up parameters */
  if (volume > 0.0 && quality > 0.0) sprintf (params, "Qpa%gq%g", volume, quality);
  else if (volume > 0.0) sprintf (params, "Qpa%g", volume);
  else if (quality > 0.0) sprintf (params, "Qpq%g", quality);
  else sprintf (params, "Qp");

  /* generate mesh */
  tetrahedralize (params, &in, &out);

  /* read triangles */
  if (out.trifacemarkerlist)
  {
    ERRMEM (surfaces = (int*)malloc (out.numberoftrifaces * sizeof (int [5]) + sizeof (int [2])));

    for (i = 0, ele = surfaces+1, tet = out.trifacelist;
	 i < out.numberoftrifaces; i ++, ele += 5, tet += 3)
    {
      ele [0] = 3;
      ele [1] = tet [0] - out.firstnumber;
      ele [2] = tet [1] - out.firstnumber;
      ele [3] = tet [2] - out.firstnumber;
      ele [4] = out.trifacemarkerlist [i];
    }
    ele [0] = 0;
  }
  else ERRMEM (surfaces = (int*)MEM_CALLOC (sizeof (int [2])));
  if (surfid == -INT_MAX) surfaces [0] = 0;
  else surfaces [0] = surfid;

  /* read tetrahedrons */
  if (volid == -INT_MAX) volid = 0;
  ERRMEM (elements =  (int*)malloc ((out.numberoftetrahedra+1) * sizeof (int [6])));
  if (out.numberofpoints == 0 || out.numberoftetrahedra == 0) goto err;
  for (i = 0, ele = elements, tet = out.tetrahedronlist;
       i < out.numberoftetrahedra; i ++, ele += 6, tet += 4)
  {
    ele [0] = 4;
    ele [1] = tet [0] - out.firstnumber;
    ele [2] = tet [1] - out.firstnumber;
    ele [3] = tet [2] - out.firstnumber;
    ele [4] = tet [3] - out.firstnumber;
    if (out.tetrahedronattributelist) ele [5] = (int)
      out.tetrahedronattributelist [i*out.numberoftetrahedronattributes]; /* the first attribute as volid */
    else ele [5] = volid;
  }
  ele [0] = 0;

  /* create output mesh */
  ret = MESH_Create ((double (*)[3])out.pointlist, elements, surfaces);

err:
  /* clean up */
  free (elements);
  free (surfaces);
  
  return ret;
}

#if 0
/* recursively gather coplanar and adjacent triangles */
static void gathercoplanar (int *mark, TRI *tri, TRI *t, int fac, MEM *setmem, SET **copla)
{
  double *n0, *n1, prod [3];
  int i;

  SET_Insert (setmem, copla, t, NULL);
  mark [t-tri] = fac;
  n0 = t->out;

  for (i = 0; i < 3; i ++)
  {
    if (t->adj [i] && !mark [t->adj[i]-tri] && t->flg == t->adj [i]->flg) /* same face marker */
    {
      n1 = t->adj [i]->out;
      PRODUCT (n0, n1, prod);
      if (DOT (n0, n1) > 0.0 && LEN (prod) < GEOMETRIC_EPSILON) /* coplanar */
      {
	gathercoplanar (mark, tri, t->adj [i], fac, setmem, copla);
      }
    }
  }
}

/* returns set of sets of coplanar and adjacent triangles and a map of unique vertices */
static SET* coplanartriangles (MEM *setmem, TRI *tri, int m)
{
  SET *coplasets, *copla;
  int i, j, *mark, fac;
  TRI *t, *e;
  MAP *item;

  fac = 0;

  coplasets = NULL;

  TRI_Compadj (tri, m);

  ERRMEM (mark = (int*)MEM_CALLOC (sizeof (int [m])));

  for (t = tri, e = t + m; t != e; t ++) /* collect coplanar and adjacent triangle sets */
  {
    if (mark [t-tri] == 0)
    {
      fac ++;
      copla = NULL;
      gathercoplanar (mark, tri, t, fac, setmem, &copla);
      SET_Insert (setmem, &coplasets, copla, NULL);
    }
  }

  free (mark);

  return coplasets;
}
#endif

/* generate tetrahedrons bounded by triangular surfaces; TRI->flg store surfids */
MESH* tetrahedralize3 (TRI *tri, int m, double volume, double quality, int volid)
{
  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  int *elements, *surfaces;
  int i, j, n, *ele, *tet;
  MAP *ver, *pol, *item;
  char params [512];
  double *a, *b;
  MEM mapmem;
  TRI *t, *e;
  MESH *ret;

  ver = NULL;
  pol = NULL;
  ret = NULL;

  /* memory pools */
  MEM_Init (&mapmem, sizeof (MAP), 128);

  /* 0-based indexing */
  in.firstnumber = 0;

  /* calculate faces and map face vertices */
  for (t = tri, e = t + m, j = n = 0; t != e; t ++, n ++)
  {
    for (i = 0; i < 3; i ++)
    {
      if (!MAP_Find_Node (ver, t->ver [i], NULL))
      {
	MAP_Insert (&mapmem, &ver, t->ver [i], (void*) (long) j, NULL);
	j ++;
      }
    }
  }

  /* input vertices */
  in.numberofpoints = MAP_Size (ver);
  in.pointlist = new REAL[in.numberofpoints * 3];
  for (item = MAP_First (ver); item; item = MAP_Next (item))
  {
    j = (int) (long) item->data;
    a = (double*) item->key;
    b = &in.pointlist [3*j];

    COPY (a, b);
  }

  /* input faces */
  in.numberoffacets = n;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];
  for (t = tri, e = t + m, n = 0; t != e; t ++, n ++)
  {
    f = &in.facetlist [n];
    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;
    p = &f->polygonlist[0];
    p->numberofvertices = 3;
    p->vertexlist = new int [p->numberofvertices];
    for (i = 0; i < 3; i ++)
    {
      item = MAP_Find_Node (ver, t->ver [i], NULL);
      ASSERT_DEBUG (item, "Inconsistent face vertex mapping");
      j = (int) (long) item->data;
      p->vertexlist [i] = j;
    }
    in.facetmarkerlist [n] = t->flg;
  }

  /* set up parameters */
  if (volume > 0.0 && quality > 0.0) sprintf (params, "Qpa%gq%g", volume, quality);
  else if (volume > 0.0) sprintf (params, "Qpa%g", volume);
  else if (quality > 0.0) sprintf (params, "Qpq%g", quality);
  else sprintf (params, "Qp");

  /* generate mesh */
  tetrahedralize (params, &in, &out);

  /* read triangles */
  ERRMEM (surfaces = (int*)malloc (out.numberoftrifaces * sizeof (int [5]) + sizeof (int [2])));
  surfaces [0] = tri->flg;

  for (i = 0, ele = surfaces+1, tet = out.trifacelist;
       i < out.numberoftrifaces; i ++, ele += 5, tet += 3)
  {
    ele [0] = 3;
    ele [1] = tet [0];
    ele [2] = tet [1];
    ele [3] = tet [2];
    ele [4] = out.trifacemarkerlist [i];
  }
  ele [0] = 0;

  /* read tetrahedrons */
  ERRMEM (elements =  (int*)malloc ((out.numberoftetrahedra+1) * sizeof (int [6])));
  if (out.numberofpoints == 0 || out.numberoftetrahedra == 0) goto err;
  for (i = 0, ele = elements, tet = out.tetrahedronlist;
       i < out.numberoftetrahedra; i ++, ele += 6, tet += 4)
  {
#if DEBUG
    double *a = &out.pointlist [tet [0]*3],
           *b = &out.pointlist [tet [1]*3],
           *c = &out.pointlist [tet [2]*3],
           *d = &out.pointlist [tet [3]*3];

    ASSERT_TEXT (simplex_J (a, b, c, d) >= 1E-15, "Zero volume tetrahedron found: %d, %d, %d, %d", tet [0], tet [1], tet [2], tet [3]);
#endif
    ele [0] = 4;
    ele [1] = tet [0];
    ele [2] = tet [1];
    ele [3] = tet [2];
    ele [4] = tet [3];
    ele [5] = volid;
  }
  ele [0] = 0;

  /* create output mesh */
  ret = MESH_Create ((double (*)[3])out.pointlist, elements, surfaces);

  /* clean up */
err:
  free (elements);
  free (surfaces);
  MEM_Release (&mapmem);
  
  return ret;
}
}
