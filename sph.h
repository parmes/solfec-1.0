/*
 * sph.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * spheres
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

#ifndef __sph__
#define __sph__

typedef struct sphere SPHERE;

/* analytical sphere */
struct sphere
{
  double ref_center [3],
	 ref_points [3][3], /* three points for the sake of global deformation marking */
         ref_radius;
  
  double cur_center [3], /* the sphere is assumed to remain spherical (nearly) when deformed */
	 cur_points [3][3],
         cur_radius;

  int surface, /* surface identifier */
      volume; /* volume identifier */

  BULK_MATERIAL *mat;
};

/* create a sphere */
SPHERE* SPHERE_Create (double *center, double radius, int surface, int volume);

/*  dummy (needed in shp.c) */
void SPHERE_Update_Adjacency (SPHERE *sph);

/* dummy (needed in shp.c) */
int SPHERE_Break_Adjacency (SPHERE *sph, double *point, double *normal);

/* create a copy of a sphere */
SPHERE* SPHERE_Copy (SPHERE *sph);

/* scaling of a sphere; scale radius: r *= vector [0]; (set ref = cur) */
void SPHERE_Scale (SPHERE *sph, double *vector);

/* translation of a sphere */
void SPHERE_Translate (SPHERE *sph, double *vector);

/* rotation of a sphere */
void SPHERE_Rotate (SPHERE *sph, double *point, double *vector, double angle);

/* cut through sphere with a plane; return triangulated cross-section; vertices in the triangles
 * point to the memory allocated after the triangles memory; adjacency is not maintained */
TRI* SPHERE_Cut (SPHERE *sph, double *point, double *normal, int *m);

/* split sphere in two with plane defined by (point, normal); surfid corresponds to the new surface;
 * topoadj != 0 implies cutting from the point and through the topological adjacency only */
void SPHERE_Split (SPHERE *sph, double *point, double *normal, short topoadj, int surfid, SPHERE **one, SPHERE **two);

/* compute partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void SPHERE_Char_Partial (SPHERE *sph, int ref, double *vo, double *sx, double *sy, double *sz, double *eul);

/* get characteristics of the sphere: volume, mass center, and Euler tensor (centered) */
void SPHERE_Char (SPHERE *sph, int ref, double *volume, double *center, double *euler);

/* update extents of an individual sphere */
void SPHERE_Extents (void *data, SPHERE *sph, double *extents);

/* compute extents of a sphere */
void SPHERE_Extents_2 (SPHERE *sph, double *extents);

/* compute oriented extents of a sphere */
void SPHERE_Oriented_Extents (SPHERE *sph, double *vx, double *vy, double *vz, double *extents);

/* return first not NULL bulk material for a sphere */
void* SPHERE_First_Bulk_Material (SPHERE *sph);

/* return sphere containing a spatial point */
SPHERE* SPHERE_Containing_Point (SPHERE *sph, double *point);

/* does this sphere contain the spatial point? */
int SPHERE_Contains_Point (void *dummy, SPHERE *sph, double *point);

/* return distance of a spatial point to the sphere */
double SPHERE_Spatial_Point_Distance (void *dummy, SPHERE *sph, double *point);

/* update sphere according to the given motion */
void SPHERE_Update (SPHERE *sph, void *body, void *shp, MOTION motion);

/* free sphere */
void SPHERE_Destroy (SPHERE *sph);

/* pack sphere into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void SPHERE_Pack (SPHERE *sph, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack sphere from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
SPHERE* SPHERE_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

/* export MBFCP definition */
void SPHERE_2_MBFCP (SPHERE *sph, FILE *out);

#endif
