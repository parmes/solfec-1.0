/*
 * ldy.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * the local dynamic problem
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

#include "mem.h"
#include "bod.h"
#include "sps.h"

#ifndef __ldy__
#define __ldy__

#define DOM_Z_SIZE          4      /* size of auxiliary storage in dom.h/constraint */

typedef struct offb OFFB;
typedef struct diab DIAB;
typedef struct locdyn LOCDYN;

enum upkind
{
  UPDIA, /* update diagonal blocks only */
  UPALL  /* update all blocks */
};

typedef enum upkind UPKIND;

enum locdyn_approach /* linearisation approaches */
{
  SEMISMOOTH_STRICT,
  SEMISMOOTH_HYBRID,
  VARIATIONAL_NONSMOOTH,
  VARIATIONAL_SMOOTHED,
};

typedef enum locdyn_approach LOCDYN_APPROACH;

#if MPI
/* eXternal Reaction */
typedef struct
{ 
  double R [3];
  int id;
  int rank; /* rank of parent block */
  short done;
  double *update; /* auxiliary update pointer */
} XR;

/* XR pointer cast */
#define XR(ptr) ((XR*)(ptr))

/* balancing algorithm */
typedef enum {LDB_OFF, LDB_GEOM, LDB_GRAPH} LDB;
#endif

/* off-diagonal
 * block entry */
struct offb
{
  double W [9], /* generalised inverse inertia block */
	 T [9], /* tangent operator block */
	 R [3]; /* prefetched dia->R */

  DIAB *dia; /* can be NULL for balanced boundary blocks */
  BODY *bod;
  OFFB *n;

#if MPI
  unsigned int id; /* adjacent constraint id */

  void *ext; /* external constraint (unbalanced W) */

  XR *x; /* external reaction for boundary blocks (dia == NULL) */
#endif
};

/* diagonal
 * block entry */
struct diab
{
  double    *R, /* average reaction => points to R[3] member of the underlying constraint */
	 U [3], /* final velocity */
	 V [3], /* initial velocity */
	 B [3], /* free velocity */
         W [9], /* generalised inverse inertia block */
	 T [9], /* tangent operator block */
	 rho;   /* scaling parameter */

  OFFB *adj;
  void *con;   /* the underlying constraint (and the owner od the reaction R[3]);
                  NULL for balanced constraints from aabb->diab */
  DIAB *p, *n;

#if MPI
  unsigned int id; /* this constraint id */

  OFFB *adjext; /* external adjacency */

  int rank; /* for a parent: rank of its child, and vice versa */

  int degree; /* the total number of blocks in this row */

  MAP *children; /* balanced boundary nodes map children ranks to local reaction REXT indices */

  SET *rext; /* balanced boundary receiving nodes store pointers to their XR here */

  /* local dynamic system entries can be migrated in parallel,
   * hence it is better to copy the necessary members of an
   * underlying constraint, in order to support independent migration */

  double REAC [3],
	 Z [DOM_Z_SIZE],
	 point [3],
	 base [9],
	 mpnt [3],
	 gap;

  short kind;

  SURFACE_MATERIAL mat;
#endif
};

/* local dynamics */
struct locdyn
{
  MEM offmem,
      diamem;

  void *dom; /* domain */
  DIAB *dia; /* list of diagonal blocks (unbalanced) */

  short modified; /* 1 if system structure has changed; otherwise 0 */

#if MPI
  MAP *insmap; /* maps newly inserted blocks to indices in 'ins' */
  DIAB **ins; /* newly inserted unbalanced blocks */
  int nins, sins; /* number of them and size of buffer */

  DIAB **del; /* recently deleted unbalanced blocks */
  int ndel, sdel; /* number and size */

  MEM mapmem, /* map items memory */
      setmem; /* set items memory */

  MAP *idbb; /* id-to-balanced diagonal block map */
  DIAB *diab; /* list of diagonal blocks (balanced) */
  int ndiab; /* number of balanced diagonal blocks */

  XR *REXT; /* table of reactions stored at other processors */
  int REXT_count; /* count of external reactions */

  LDB ldb, /* kind of balancing (default: LDB_OFF) */
      ldb_new; /* newly set kind of balancing */

  struct Zoltan_Struct *zol;

  double imbalance_tolerance;/* imbalance threshold */

  int nexpdia; /* number of exported rows */
#endif
};

/* create local dynamics for a domain */
LOCDYN* LOCDYN_Create (void *dom);

/* insert a 'con'straint between a pair of bodies =>
 * return the diagonal entry of the local dynamical system */
DIAB* LOCDYN_Insert (LOCDYN *ldy, void *con, BODY *one, BODY *two);

/* remove a diagonal entry from local dynamics */
void LOCDYN_Remove (LOCDYN *ldy, DIAB *dia);

/* update local dynamics => prepare for a solution */
void LOCDYN_Update_Begin (LOCDYN *ldy, UPKIND upkind);

/* update local dynamics => after the solution */
void LOCDYN_Update_End (LOCDYN *ldy);

#if MPI
/* change load balancing algorithm */
void LOCDYN_Balancing (LOCDYN *ldy, LDB ldb);

/* update mapping of balanced external reactions */
void LOCDYN_REXT_Update (LOCDYN *ldy);

/* return the union of 'inp' sets; return the communication 'pattern' used
 * to gather and scatter reactions in the union; if score < 0 the same union
 * set is created on all processors; if score >= 0 a single set is created
 * on the processor with the minimum score, while other processors get NULL */
SET* LOCDYN_Union_Create (LOCDYN *ldy, SET *inp, int score, void **pattern);

/* gather reactions */
void LOCDYN_Union_Gather (void *pattern);

/* scatter reactions */
void LOCDYN_Union_Scatter (void *pattern);

/* release memory used by the union set */
void LOCDYN_Union_Destroy (void *pattern);
#endif

/* set an approach to the linearisation of local dynamics */
void LOCDYN_Approach (LOCDYN *ldy, LOCDYN_APPROACH approach);

/* assemble tangent operator */
void LOCDYN_Tangent (LOCDYN *ldy);

/* compute merit function */
double LOCDYN_Merit (LOCDYN *ldy);

/* free memory */
void LOCDYN_Destroy (LOCDYN *ldy);

#endif
