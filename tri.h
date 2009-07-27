/*
 * tri.h
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

#ifndef __tri__
#define __tri__

typedef struct triangle TRI; /* surface triangle */
struct triangle
{
  double out [3], /* outward normal */
        *ver [3]; /* vertices (CCW) */
  TRI   *adj [3]; /* adjacent triangles (ordered: adj[0] <=> ver[0]-ver[1], ...) */
  int    flg;     /* flag(s) => used in external algorithms */
};

typedef struct polar_face_vertex PFV; /* vertex of a polar face */
struct polar_face_vertex
{
  double *coord, /* vertex coordinate => this is -pla[0-2]/pla[3] */
	 *nl;    /* outward normal => points to 'ver' in TRI, hence not normalised */
  PFV    *n;     /* next CCW vertex in a circular list */
};

/* sort adjacency as specified above */
void TRI_Sortadj (TRI *tri);

/* make 'p' and 'q' adjacent, assuming
 * there is enough NULL space in 'adj's;
 * return 0 on fault, 1 otherwise */
int TRI_Addadj (TRI *p, TRI *q);

/* copy triangles into a continuous memory block;
 * vertices are placed right after the returned table */
TRI* TRI_Copy (TRI *tri, int n);

/* compute adjacency structure */
void TRI_Compadj (TRI *tri, int n);

/* input => convex polyhedron containing zero (tri, n);
 * output => polar polyhedron defined by 'm' vertex lists;
 * a continuous block of memory is returned; 'nl's point to 'ver'
 * members in 'tri'; 'coord's point within the returned block */
PFV* TRI_Polarise (TRI *tri, int n, int *m);

/* extract vertices of triangulation (tri, n)
 * into a table of size (double [3]) x m */
double* TRI_Vertices (TRI *tri, int n, int *m);

/* extract planes of triangulation (tri, n)
 * into a table of size (double [6]) x m, where
 * each plane is represented by (normal, point) */
double* TRI_Planes (TRI *tri, int n, int *m);

#endif
