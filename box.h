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

#if ZOLTAN
#include <zoltan.h>
#endif

#include "shp.h"
#include "mem.h"
#include "map.h"
#include "set.h"

#ifndef BODY_TYPE
#define BODY_TYPE
typedef struct general_body BODY;
#endif

#ifndef DOMAIN_TYPE
#define DOMAIN_TYPE
typedef struct domain DOM;
#endif

#ifndef __box__
#define __box__

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

#define BOXALG_COUNT (HASH3D+1) /* count of overlap algorithms */
typedef struct objpair OPR; /* pointer pair used for exclusion tests */
typedef struct aabb AABB; /* overlap detection driver data */
typedef enum boxalg BOXALG; /* type of overlap detection algorithm */
typedef void (*BOX_Overlap_Create)  (void *data, BOX *one, BOX *two); /* created overlap callback => returns a user pointer */
typedef void (*BOX_Extents_Update)  (void *data, void *gobj, double *extents); /* extents update callback */

#ifndef BOX_TYPE
#define BOX_TYPE
typedef struct box BOX; /* axis aligned box */
#endif

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

  BOX_Extents_Update update; /* extents update callback => update (data, gobj, extents) */

  GOBJ kind; /* kind of a geometric object */

  SGP *sgp; /* shape and geometric object pair */

  BODY *body; /* owner of the shape */

  void *mark; /* auxiliary marker used by hashing algorithms */

#if MPI
  SET *ranks; /* ranks where this box overalps */
#endif

  BOX *prev, *next; /* list pointers */
};

/* box pointer cast */
#define BOX(box) ((BOX*)(box))

/* this routine returns a combined code of two
 * geomtric objects bounded by the respective boxes */
inline static short GOBJ_Pair_Code (BOX *one, BOX *two)
{ short a = one->kind, b = two->kind; return ((a << 8) | b); }

/* geometric object pair code based on object kinds */
#define GOBJ_Pair_Code_Ext(a, b) (((a) << 8) | (b))

/* and here are the codes */
#define AABB_ELEMENT_ELEMENT	0x0101
#define AABB_CONVEX_CONVEX	0x0202
#define AABB_SPHERE_SPHERE	0x0404
#define AABB_ELLIP_ELLIP        0x1010

#define AABB_ELEMENT_CONVEX	0x0102
#define AABB_CONVEX_ELEMENT	0x0201

#define AABB_ELEMENT_SPHERE	0x0104
#define AABB_SPHERE_ELEMENT	0x0401

#define AABB_ELEMENT_ELLIP      0x0110
#define AABB_ELLIP_ELEMENT	0x1001

#define AABB_CONVEX_SPHERE	0x0204
#define AABB_SPHERE_CONVEX	0x0402

#define AABB_CONVEX_ELLIP       0x0210
#define AABB_ELLIP_CONVEX	0x1002

#define AABB_ELLIP_SPHERE	0x1004
#define AABB_SPHERE_ELLIP	0x0410

/* driver data */
struct aabb
{
  BOX *lst, /* current list of boxes */
     **tab; /* pointer table of boxes from the current list */

  MEM boxmem, /* box memory pool */
      mapmem, /* map memory pool */
      setmem, /* set memory pool */
      oprmem; /* object pair memory pool */

  SET *nobody, /* set of body pairs excluded from overlap tests */
      *nogobj; /* set of geometric object pairs excluded from overlap tests */

  int boxnum, /* number of boxes */
      tabsize; /* size of pointer table */

  char modified; /* modification flag => for time coherence */

  void *swp,  /* sweep plane data */
       *hsh;  /* hashing data */

  DOM *dom; /* the underlying domain */
};

/* get algorithm name string */
char* AABB_Algorithm_Name (BOXALG alg);

/* create box overlap driver data */
AABB* AABB_Create (int size);

/* insert geometrical object => return the associated box */
BOX* AABB_Insert (AABB *aabb, BODY *body, GOBJ kind, SGP *sgp, BOX_Extents_Update update);

/* delete geometrical gobject associated with the box */
void AABB_Delete (AABB *aabb, BOX *box);

/* insert a body */
void AABB_Insert_Body (AABB *aabb, BODY *body);

/* delete a body */
void AABB_Delete_Body (AABB *aabb, BODY *body);

/* update state and detect all created overlaps (aabb->dom != NULL) */
void AABB_Update (AABB *aabb, BOXALG alg, void *data, BOX_Overlap_Create create);

/* update state and detect all created overlaps (aabb->dom == NULL) */
void AABB_Simple_Detect (AABB *aabb, BOXALG alg, void *data, BOX_Overlap_Create create); 

/* never report overlaps betweem this pair of bodies (given by identifiers) */
void AABB_Exclude_Body_Pair (AABB *aabb, unsigned int id1, unsigned int id2);

/* undo the effect of AABB_Exclude_Body_Pair */
void AABB_Include_Body_Pair (AABB *aabb, unsigned int id1, unsigned int id2);

/* pack excluded body pairs ids */
void AABB_Pack_Body_Pairs (AABB *aabb, int *isize, int **i, int *ints);

/* unpack excluded body pairs ids (and include them into the set) */
void AABB_Unpack_Body_Pairs (AABB *aabb, int *ipos, int *i, int ints);

/* never report overlaps betweem this pair of objects (bod1, sgp1), (bod1, sgp2) */
void AABB_Exclude_Gobj_Pair (AABB *aabb, unsigned int bod1, int sgp1, unsigned int bod2, int sgp2);

/* release memory */
void AABB_Destroy (AABB *aabb);

/* get geometrical object extents update callback */
BOX_Extents_Update SGP_Extents_Update (SGP *sgp);

#endif
