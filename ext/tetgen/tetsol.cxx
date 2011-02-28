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

/* generate tetrahedrons based on an input mesh object; pass -INT_MAX for (vol/surf)ids to inherit from the mesh */
MESH* tetrahedralize1 (MESH *shape, double volume, double quality, int volid, int surfid)
{
  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  int *elements, *surfaces;
  int i, j, n, *ele, *tet;
  char params [512];
  MESH *ret = NULL;
  MAP *map, *item;
  double *a, *b;
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
    if (surfid > -INT_MAX) in.facetmarkerlist [n] = surfid;
    else in.facetmarkerlist [n] = fac->surface;
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

  /* read tetrahedrons */
  if (volid == -INT_MAX) volid = shape->surfeles->volume;
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
    ele [5] = volid; /* FIXME: use kd-tree based mid-point to volume id mapping */
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
      ele [1] = tet [0];
      ele [2] = tet [1];
      ele [3] = tet [2];
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
    if (out.tetrahedronattributelist) ele [5] = (int) out.tetrahedronattributelist [i];
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
}
