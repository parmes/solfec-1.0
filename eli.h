/*
 * eli.h
 * Copyright (C) 2011, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * ellipsoids
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
#include "mot.h"
#include "tri.h"
#include "cvx.h"

#ifndef __eli__
#define __eli__

typedef struct ellipsoid ELLIP;

/* analytical ellipsoid */
struct ellipsoid
{
  double ref_center [3],
	 ref_points [3][3];
  
  double cur_center [3],
	 cur_points [3][3];

  int surface, /* surface identifier */
      volume; /* volume identifier */

  BULK_MATERIAL *mat;
};

/* create an ellipsoid */
ELLIP* ELLIP_Create (double *center, double radius, int surface, int volume);

/*  dummy (needed in shp.c) */
void ELLIP_Update_Adjacency (ELLIP *eli);

/* dummy (needed in shp.c) */
int ELLIP_Break_Adjacency (ELLIP *eli, double *point, double *normal);

/* create a copy of a ellipsoid */
ELLIP* ELLIP_Copy (ELLIP *eli);

/* scaling of a ellipsoid; scale radius: r *= vector [0]; (set ref = cur) */
void ELLIP_Scale (ELLIP *eli, double *vector);

/* translation of a ellipsoid */
void ELLIP_Translate (ELLIP *eli, double *vector);

/* rotation of a ellipsoid */
void ELLIP_Rotate (ELLIP *eli, double *point, double *vector, double angle);

/* cut through ellipsoid with a plane; return triangulated cross-section; vertices in the triangles
 * point to the memory allocated after the triangles memory; adjacency is not maintained */
TRI* ELLIP_Cut (ELLIP *eli, double *point, double *normal, int *m);

/* split ellipsoid in two with plane defined by (point, normal); surfid corresponds to the new surface;
 * topoadj != 0 implies cutting from the point and through the topological adjacency only */
void ELLIP_Split (ELLIP *eli, double *point, double *normal, short topoadj, int surfid, CONVEX **one, CONVEX **two);

/* compute partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void ELLIP_Char_Partial (ELLIP *eli, int ref, double *vo, double *sx, double *sy, double *sz, double *eul);

/* get characteristics of the ellipsoid: volume, mass center, and Euler tensor (centered) */
void ELLIP_Char (ELLIP *eli, int ref, double *volume, double *center, double *euler);

/* update extents of an individual ellipsoid */
void ELLIP_Extents (void *data, ELLIP *eli, double *extents);

/* compute extents of a ellipsoid */
void ELLIP_Extents_2 (ELLIP *eli, double *extents);

/* compute oriented extents of a ellipsoid */
void ELLIP_Oriented_Extents (ELLIP *eli, double *vx, double *vy, double *vz, double *extents);

/* return first not NULL bulk material for a ellipsoid */
void* ELLIP_First_Bulk_Material (ELLIP *eli);

/* return ellipsoid containing a spatial point */
ELLIP* ELLIP_Containing_Point (ELLIP *eli, double *point);

/* does this ellipsoid contain the spatial point? */
int ELLIP_Contains_Point (void *dummy, ELLIP *eli, double *point);

/* return distance of a spatial point to the ellipsoid */
double ELLIP_Spatial_Point_Distance (void *dummy, ELLIP *eli, double *point);

/* update ellipsoid according to the given motion */
void ELLIP_Update (ELLIP *eli, void *body, void *shp, MOTION motion);

/* free ellipsoid */
void ELLIP_Destroy (ELLIP *eli);

/* pack ellipsoid into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void ELLIP_Pack (ELLIP *eli, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack ellipsoid from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
ELLIP* ELLIP_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

/* export MBFCP definition */
void ELLIP_2_MBFCP (ELLIP *eli, FILE *out);

#endif
