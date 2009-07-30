/*
 * shp.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * shape definition
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

#include "mot.h"

#ifndef __shp__
#define __shp__

typedef struct shape SHAPE;
typedef struct shape_gobj_pair SGP;

/* general shape */
struct shape
{
  enum {SHAPE_MESH = 0,
        SHAPE_CONVEX,
        SHAPE_SPHERE} kind; /* kind of shape */

  void *data; /* representation */

  SHAPE *next;
};

/* shape and geometrical object pair */
struct shape_gobj_pair
{
  SHAPE *shp;
  void *gobj;
  void *box; /* points to the bounding box */
};

/* create a general shape */
SHAPE* SHAPE_Create (short kind, void *data);

/* create shape geometric object pairs */
SGP* SGP_Create (SHAPE *shp, int *nsgp);

/* glue two shape lists (gluing together basic shapes) */
SHAPE* SHAPE_Glue (SHAPE *shp, SHAPE *shq);

/* glue two shape lists (without gluing basic shapes) */
SHAPE* SHAPE_Glue_Simple (SHAPE *shp, SHAPE *shq);

/* update adjacency data of stored shapes;
 * no other function affects the adjacency */
void SHAPE_Update_Adjacency (SHAPE *shp);

/* scale cur shape => 
 * if MESH,  scale each: x *= vector [0], y *= vector [1], z *= vector [2];
 * if SPHERE, scale radius: r *= vector [0]; (set ref = cur) */
void SHAPE_Scale (SHAPE *shp, double *vector);

/* translate cur shape (set ref = cur) */
void SHAPE_Translate (SHAPE *shp, double *vector);

/* rotate cur shape (set ref = cur), around the line (point, vector) */
void SHAPE_Rotate (SHAPE *shp, double *point, double *vector, double angle);

/* get cur characteristics => volume, mass center, and Euler tensor (centered) */
void SHAPE_Char (SHAPE *shp, double *volume, double *center, double *euler);

/* return an object containing spatial point */
void* SHAPE_Gobj (SHAPE *shp, double *point, SHAPE **out);

/* return an index of object containing spatial point (or -1 on failure) */
int SHAPE_Sgp (SGP *sgp, int nsgp, double *point);

/* update current shape with given motion */
void SHAPE_Update (SHAPE *shp, void *body, MOTION motion);

/* copute shape extents */
void SHAPE_Extents (SHAPE *shp, double *extents);

/* release shape memory */
void SHAPE_Destroy (SHAPE *shp);

/* release shape wrapper memory (without data) */
void SHAPE_Destroy_Wrapper (SHAPE *shp);

/* pack shape into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void SHAPE_Pack (SHAPE *shp, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack shape from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
SHAPE* SHAPE_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

#endif
