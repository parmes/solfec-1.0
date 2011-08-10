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

#include <stdio.h>
#include "mot.h"
#include "tri.h"

#ifndef BOX_TYPE
#define BOX_TYPE
typedef struct box BOX;
#endif

#ifndef __shp__
#define __shp__

typedef struct shape SHAPE; /* general shape */
typedef struct shape_gobj_pair SGP; /* shape and geometric object pair */

/* geometrical objects */
enum gobj
{
  GOBJ_DUMMY   = 0x00,
  GOBJ_ELEMENT = 0x01, /* XXX: never edit these without looking into box.h: e.g. AABB_ELEMENT_ELEMENT, etc. */
  GOBJ_CONVEX  = 0x02,
  GOBJ_SPHERE  = 0x04,
  GOBJ_NODE    = 0x08
};

typedef enum gobj GOBJ; /* kind of geometrical object */

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
  GOBJ kind;
  BOX *box;
};

/* SGP creation flags */
enum sgp_flags
{
  SGP_MESH_NODES = 0x01
};

typedef enum sgp_flags SGP_FLAGS;

/* create a general shape */
SHAPE* SHAPE_Create (short kind, void *data);

/* create shape geometric object pairs */
SGP* SGP_Create (SHAPE *shp, int *nsgp, SGP_FLAGS flags);

/* SGP to GOBJ conversion for FEM purposes */
void* SGP_2_GOBJ (SGP *sgp);

/* get GOBJ type of given shape */
GOBJ SHAPE_2_GOBJ (SHAPE *shp);

/* copy shape */
SHAPE* SHAPE_Copy (SHAPE *shp);

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

/* cut through shape with a plane; return triangulated cross-section; all returned data
 * points to the memory allocated after the triangles memory; adjacency is not maintained;
 * TRI->ptr stores a pointer to the geometrical object that has been cut by the triangle;
 * (body, cur_to_ref, ref, cur, n) can be either all NULL or all valid pointers; if not NULL then
 * 'n' reference and current vertices are calculated (triagnle vertices are the current ones) */
TRI* SHAPE_Cut (SHAPE *shp, double *point, double *normal, int *m,
  void *body, MOTION cur_to_ref, SGP **sgp, double **ref, double **cur, int *n);

/* split shape by plane; output two parts of the split shape;
 * topoadj != 0 implies cutting from the point and through the topological adjacency only */
void SHAPE_Split (SHAPE *shp, double *point, double *normal, short topoadj, int surfid, SHAPE **one, SHAPE **two);

/* get spatial/referential characteristics => volume, mass center, and Euler tensor (centered) */
void SHAPE_Char (SHAPE *shp, int ref, double *volume, double *center, double *euler);

/* for the given shape (not a list) compute spatial/referential partial characteristic: 'vo'lume and
 * static momenta 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void SHAPE_Char_Partial (SHAPE *shp, int ref, double *vo, double *sx, double *sy, double *sz, double *eul);

/* return an object containing spatial point */
void* SHAPE_Gobj (SHAPE *shp, double *point, SHAPE **out);

/* return an index of object containing spatial point (or -1 on failure) */
int SHAPE_Sgp (SGP *sgp, int nsgp, double *point);

/* return an index of object closest to the spatial point (output distance in 'd' if not NULL) */
int SHAPE_Closest_Sgp (SGP *sgp, int nsgp, double *point, double *d);

/* update current shape with given motion;
 * body == motion == NULL => referential update */
void SHAPE_Update (SHAPE *shp, void *body, MOTION motion);

/* copute shape extents */
void SHAPE_Extents (SHAPE *shp, double *extents);

/* copute shape oriented extents in corrds given by three direction vectors */
void SHAPE_Oriented_Extents (SHAPE *shp, double *vx, double *vy, double *vz, double *extents);

/* return first bulk material recorded
 * in this individual shape (not a list) */
void* SHAPE_First_Bulk_Material (SHAPE *shp);

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

/* export MBFCP definition */
void SHAPE_2_MBFCP (SHAPE *shp, FILE *out);

#endif
