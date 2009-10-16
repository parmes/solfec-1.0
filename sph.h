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

  SPHERE **adj; /* adjacency */

  int surface, /* surface identifier */
      nadj, /* number of neighbours */
      volume; /* volume identifier */

  BULK_MATERIAL *mat;

  SPHERE *next;
};

/* create a sphere (sph == NULL) or append spheres list with another sphere */
SPHERE* SPHERE_Create (SPHERE *sph, double *center, double radius, int surface, int volume);

/* glue two sphere lists */
SPHERE* SPHERE_Glue (SPHERE *sph, SPHERE *spg);

/* update adjacency data of spheres;
 * no other function affects the adjacency */
void SPHERE_Update_Adjacency (SPHERE *sph);

/* create a copy of a list */
SPHERE* SPHERE_Copy (SPHERE *sph);

/* scaling of a list; scale radius: r *= vector [0]; (set ref = cur) */
void SPHERE_Scale (SPHERE *sph, double *vector);

/* translation of a list */
void SPHERE_Translate (SPHERE *sph, double *vector);

/* rotation of a list */
void SPHERE_Rotate (SPHERE *sph, double *point, double *vector, double angle);

/* split sphere lists in two lists with plane defined by (point, normal);
 * adjacencies between the split lists elements need to be recomputed */
void SPHERE_Split (SPHERE *sph, double *point, double *normal, SPHERE **one, SPHERE **two);

/* compute current partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void SPHERE_Char_Partial (SPHERE *sph, double *vo, double *sx, double *sy, double *sz, double *eul);

/* get 'cur' characteristics of the spheres in list:
 * volume, mass center, and Euler tensor (centered) */
void SPHERE_Char (SPHERE *sph, double *volume, double *center, double *euler);

/* update extents of an individual sphere */
void SPHERE_Extents (void *data, SPHERE *sph, double *extents);

/* compute extents of sphere list */
void SPHERE_List_Extents (SPHERE *sph, double *extents);

/* return first not NULL bulk material for a sphere list */
void* SPHERE_First_Bulk_Material (SPHERE *sph);

/* return sphere containing the point */
SPHERE* SPHERE_Containing_Point (SPHERE *sph, double *point);

/* does this sphere (not a list) contain the point? */
int SPHERE_Contains_Point (void *dummy, SPHERE *sph, double *point);

/* update sphere list according to the given motion */
void SPHERE_Update (SPHERE *sph, void *body, void *shp, MOTION motion);

/* test wether two spheres are adjacent */
int SPHERE_Adjacent (SPHERE *one, SPHERE *two);

/* free list of spheres */
void SPHERE_Destroy (SPHERE *sph);

/* pack sphere(s) into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void SPHERE_Pack (SPHERE *sph, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack sphere(s) from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
SPHERE* SPHERE_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

#endif
