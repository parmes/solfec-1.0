/*
 * box.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * axis aligned bounding box overlap detection
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

#if MPI
#include <zoltan.h>
#endif

#include "shp.h"
#include "mem.h"
#include "map.h"
#include "set.h"

#ifndef __box__
#define __box__

/* geometrical
 * objects */
enum gobj
{
  GOBJ_DUMMY   = 0x00,
  GOBJ_ELEMENT = 0x01,
  GOBJ_CONVEX  = 0x02,
  GOBJ_SPHERE  = 0x04,
  GOBJ_NEW     = 0x08
};

/* detection
 * algorithms */
enum boxalg
{
  SWEEP_HASH2D_LIST = 0, /* same order as DRALG in dyr.h ... */
  SWEEP_HASH2D_XYTREE, 
  SWEEP_XYTREE,
  SWEEP_HASH1D_XYTREE, /* ... until here */
  HYBRID,
  HASH3D
};

typedef struct objpair OPR; /* pointer pair used for exclusion tests */
typedef struct box BOX; /* axis aligned box */
typedef struct aabb AABB; /* overlap detection driver data */
typedef enum gobj GOBJ; /* kind of geometrical object */
typedef enum boxalg BOXALG; /* type of overlap detection algorithm */
typedef void* (*BOX_Overlap_Create)  (void *data, BOX *one, BOX *two); /* created overlap callback => returns a user pointer */
typedef void  (*BOX_Overlap_Release) (void *data, BOX *one, BOX *two, void *user); /* released overlap callback => uses the user pointer */
typedef void  (*BOX_Extents_Update)  (void *data, void *gobj, double *extents); /* extents update callback */

/* object pair */
struct objpair
{
  unsigned int bod1,
	       bod2;
  int sgp1,
      sgp2;
};

/* bounding box */
struct box
{
  double extents [6]; /* min x, y, z, max x, y, z */

  void *data; /* extents update data, usually => sgp->shp->data */
  BOX_Extents_Update update; /* extents update callback => update (data, gobj, extents) */

  GOBJ kind; /* kind of a geometric object */

  SGP *sgp; /* shape and geometric object pair */

  void *body, /* owner of the shape */
       *mark; /* auxiliary marker used by hashing algorithms */

  MAP *adj; /* map of adjacent boxes to user pointers returned by the overlap create callback */

  BOX *prev, /* previous in list */
      *next; /* next in list */

#if MPI
  SET *children;
#endif
};

/* this routine returns a combined code of two
 * geomtric objects bounded by the respective boxes */
inline static short GOBJ_Pair_Code (BOX *one, BOX *two)
{ short a = one->kind, b = two->kind; return ((a << 8) | b); }

/* geometric object pair code based on object kinds */
#define GOBJ_Pair_Code_Ext(a, b) (((a) << 8) | (b))

/* and here are the codes */
/* and here are the codes */
#define AABB_ELEMENT_ELEMENT	0x0101
#define AABB_CONVEX_CONVEX	0x0202
#define AABB_SPHERE_SPHERE	0x0404

#define AABB_ELEMENT_CONVEX	0x0102
#define AABB_CONVEX_ELEMENT	0x0201

#define AABB_ELEMENT_SPHERE	0x0104
#define AABB_SPHERE_ELEMENT	0x0401

#define AABB_CONVEX_SPHERE	0x0204
#define AABB_SPHERE_CONVEX	0x0402

/* driver data */
struct aabb
{
  BOX *lst, /* current list of boxes */
      *in, /* list of boxes inserted before an update */
      *out, /* list of boxes deleted before an update */
     **tab; /* pointer table of boxes from the current list */

  MEM boxmem, /* box memory pool */
      mapmem, /* map memory pool */
      setmem, /* set memory pool */
      oprmem; /* object pair memory pool */

  SET *nobody, /* set of body pairs excluded from overlap tests */
      *nogobj; /* set of geometric object pairs excluded from overlap tests */

  int nin, /* length of the insertion list 'in' */
      nlst, /* length of the current box list 'lst' */
      ntab;  /* length of the pointer table 'tab' */

  char modified; /* modification flag => for time coherence */

  void *swp,  /* sweep plane data */
       *hsh;  /* hashing data */

  void *dom; /* the underlying domain */

#if MPI
  int nobody_modified, /* used to synchronise 'nobody' */
      nogobj_modified; /* used to synchronise 'nogobj' */

  BOX **aux; /* auxiliary table of boxes used during balancing */

  struct Zoltan_Struct *zol;
#endif
};

/* create box overlap driver data */
AABB* AABB_Create (int size);

/* insert geometrical object => return the associated box */
BOX* AABB_Insert (AABB *aabb, void *body, GOBJ kind, SGP *sgp, void *data, BOX_Extents_Update update);

/* delete an object (either 'gobj' or 'box' must be not NULL) */
void AABB_Delete (AABB *aabb, BOX *box);

/* delete all data related to this body */
void AABB_Delete_Body (AABB *aabb, void *body);

/* update state => detect created and released overlaps */
void AABB_Update (AABB *aabb, BOXALG alg, void *data, BOX_Overlap_Create create, BOX_Overlap_Release release);

/* never report overlaps betweem this pair of bodies (given by identifiers) */
void AABB_Exclude_Body_Pair (AABB *aabb, unsigned int id1, unsigned int id2);

/* never report overlaps betweem this pair of objects (bod1, sgp1), (bod1, sgp2) */
void AABB_Exclude_Gobj_Pair (AABB *aabb, unsigned int bod1, int sgp1, unsigned int bod2, int sgp2);

/* break box adjacency if boxes associated with the objects are adjacent
 * (this will cose re-detection of the overlap during the next update) */
void AABB_Break_Adjacency (AABB *aabb, BOX *one, BOX *two);

/* release memory */
void AABB_Destroy (AABB *aabb);

#endif
