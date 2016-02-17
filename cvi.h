/*
 * cvi.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * intersection of two convex polyhedrons
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

#include "tri.h"

#ifndef __cvi__
#define __cvi__

typedef enum {REGULARIZED,     /* filters out zero-volume intersections */ 
              NON_REGULARIZED} /* includes surface-to-surface zero-volume intersections */
	      CVIKIND;

/* compute intersection of two convex polyhedrons:
 * (va, nva) are vertices of polyhedron 'a' (3-vectors),
 * (pa, npa) are planes of polyhedron 'b' (6-vectors: normal, point),
 * similarly for polyhedron 'b'; the returned 'm' triangles bounding
 * the surface mesh of the intersection; the 'flg' member in TRI is set either
 * to a positive or to a negative index, depending on the origin of the triangle
 * in the plane of polyhedron 'a' (positive) or 'b' (negative); pointers in the
 * returned TRI table reference the memory placed in the same block;
 * the adjacency structure in the returned mesh is not set;
 * 'pv' if not NULL, points to the vertex memory (part of 'tri' memory block);
 * 'nv' if not NULL, is the number of vertices of the intersection convex */
TRI* cvi (double *va, int nva, double *pa, int npa,
          double *vb, int nvb, double *pb, int npb,
	  CVIKIND kind, int *m, double **pv, int *nv);

#endif
