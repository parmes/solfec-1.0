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

/* off-diagonal
 * block entry */
struct offb
{
  double W [9], /* generalised inverse inertia block */
	 T [9], /* tangent operator block */
	 R [3]; /* prefetched dia->R */

  DIAB *dia;
  BODY *bod;
  OFFB *n;
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
  void *con;   /* the underlying constraint an the owner od the reactopn R[3] */
  DIAB *p, *n;

#if MPI
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
  struct Zoltan_Struct *zol;

  DIAB *diab; /* list of diagonal blocks (balanced) */
#endif
};

/* create local dynamics for a domain */
LOCDYN* LOCDYN_Create (void *dom);

/* insert a 'con'straint between a pair of bodies =>
 * return the diagonal entry of the local dynamical system */
DIAB* LOCDYN_Insert (LOCDYN *ldy, void *con, BODY *one, BODY *two);

/* remove a diagonal entry from local dynamics */
void LOCDYN_Remove (LOCDYN *ldy, DIAB *dia);

#if MPI
/* insert an external constraint */
void LOCDYN_Insert_Ext (LOCDYN *ldy, void *con);

/* remove all external constraints */
void LOCDYN_Remove_Ext_All (LOCDYN *ldy);
#endif

/* update local dynamics => prepare for a solution */
void LOCDYN_Update_Begin (LOCDYN *ldy, UPKIND upkind);

/* update local dynamics => after the solution */
void LOCDYN_Update_End (LOCDYN *ldy);

/* free memory */
void LOCDYN_Destroy (LOCDYN *ldy);

#endif
