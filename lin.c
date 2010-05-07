/*
 * lin.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * constraints linearization routines
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

#include <complex.h>
#include <stdlib.h>
#include <time.h>
#include "dom.h"
#include "alg.h"
#include "bla.h"
#include "lap.h"
#include "lin.h"
#include "err.h"
#include "mtx.h"
#include "ext/csparse.h"
#include "ext/krylov/krylov.h"

#if MPI
#include "com.h"
#endif

#define USE_AT 0 /* disable for now */
#define USE_LU (!MPI && 0) /* user LU in serial case */

#define BLOCKS 256

struct vect
{
  double *x;

  int n;
};

#define vect(ptr) ((struct vect*)ptr)

struct linsys
{
  LINVAR variant;

  LOCDYN *ldy;

  SET *subset; /* used constraints */

  MEM setmem,
      Y_mem,
      T_mem;

  double fixed_point_tol; /* fixed point normal stress update error tolerance */

  double delta; /* Tikhonov regularaization */

  double epsilon; /* smoothing thickness */

#if USE_LU
  MEM mapmem; /* {DIAB, OFFB}->map */
  MX *A;
#endif

  struct vect *x, *b;

  hypre_FlexGMRESFunctions *gmres_functions;
  void *gmres_vdata;
  int iters;
  double resnorm;

#if MPI
  COMDATA *send, *recv;
  int nsend, nrecv;
  void *pattern;
  double **R; /* external con->R */
  int *basenum; /* initial constraint numbers on each processor */
#endif
};

#define DR(con) (con)->Y
#define DU(con) ((con)->Y+3)
#define RE(con) ((con)->Y+6)
#define RN(con) ((con)->Y+7) [0]

#if USE_LU
/* get compressed structure of system matrix */
static MX* matrix_create (LINSYS *sys)
{
  int i, j, k, l, n, ncon;
  short boundary;
  SET *jtem;
  MAP *item;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  MEM mem;
  MX *A;

  ncon = SET_Size (sys->subset);
  boundary = sys->variant & BOUNDARY_ONLY;
  MEM_Init (&sys->mapmem, sizeof (int [3]), BLOCKS);
  ERRMEM (A = MEM_CALLOC (sizeof (MX)));
  A->kind = MXCSC;

  /* matrix structure */
  MAP **pcol;

  /* allocate colum blocks mapping */
  ERRMEM (pcol = MEM_CALLOC (ncon * sizeof (MAP*)));
  MEM_Init (&mem, sizeof (MAP), BLOCKS);

  /* map colum blocks */
  for (jtem = SET_First (sys->subset); jtem; jtem = SET_Next (jtem))
  {
    int col, row, *map;

    con = jtem->data;
    dia = con->dia;

    col = row = con->num;
    ASSERT_DEBUG (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Diagonal block mapped twice");
    ERRMEM (map = MEM_Alloc (&sys->mapmem));
    ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL)); /* '+1' not to confuse with NULL */

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;

#if MPI
      if (boundary && con->dia->adjext == NULL) continue; /* skip internal adjacency */
#endif

      col = con->num;

      if (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL)) /* there are two W(i,j) blocks for two-body constraints */
      {
	ERRMEM (map = MEM_Alloc (&sys->mapmem));
	ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL));
      }
    }
  }

  for (j = 0, A->nzmax = 0; j < ncon; j ++)
  {
    A->nzmax += 9 * MAP_Size (pcol [j]); /* number of nonzeros */
  }
  A->n = A->m = ncon * 3; /* system dimension */
  A->nz = -1; /* compressed columns */

  /* eallocate compressed column storage */
  ERRMEM (A->x = MEM_CALLOC (sizeof (double) * A->nzmax));
  ERRMEM (A->i = malloc (sizeof (int) * A->nzmax));
  ERRMEM (A->p = MEM_CALLOC (sizeof (int) * (A->n+1))); /* '+ 1' as there is cols[0] == 0 always */

  int *rows = A->i,
      *cols = A->p,
      *aux;

  ERRMEM (aux = malloc (sizeof (int) * A->n));

  /* set up column sizes */
  for (j = 0; j < ncon; j ++)
  {
    n = 3 * MAP_Size (pcol [j]);

    cols [3*j+1] = n; /* '+1' so that cols[0] == 0 */
    cols [3*j+2] = n;
    cols [3*j+3] = n;
  }

  /* compute column pointers */
  for (n = 1; n <= A->n; n ++)
  {
    cols [n] += cols [n-1];
    aux [n-1] = cols [n-1]; /* initialize aux with cols  */
  }

  ASSERT_DEBUG (cols [A->n] == A->nzmax, "Inconsistent sparse storage");

  /* set up row pointers */
  for (j = 0; j < ncon; j ++)
  {
    for (item = MAP_First (pcol [j]); item; item = MAP_Next (item))
    {
      k = 3 * (((int) (long) item->key) - 1); /* base row block index (-1 as keys == row+1) */

      int *map = item->data;

      ASSERT_DEBUG (map [0] == 0 && map [1] == 0 && map [2] == 0, "Reusing the same map triplet: %d, %d, %d", map [0], map [1], map [2]);

      for (n = 0; n < 3; n ++) map [n] = aux [3*j+n]; /* set row blocks map */

      for (n = 0; n < 3; n ++) /* for each column block sub-column */
      {
	for (l = aux [3*j+n], i = 0; i < 3; i ++) /* for each row block sub-row */
	{
	  rows [l + i] = k + i; /* set up row index */
	}
	aux [3*j+n] += 3; /* increment relative row pointers */
      }
    }
  }

  /* assign sparse column maps to blocks */
  for (jtem = SET_First (sys->subset); jtem; jtem = SET_Next (jtem))
  {
    int col, row;

    con = jtem->data;
    dia = con->dia;

    col = row = con->num;

    ASSERT_DEBUG_EXT (dia->map = MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Inconsistent sparse storage");

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;

#if MPI
      if (boundary && con->dia->adjext == NULL) continue; /* skip internal adjacency */
#endif

      col = con->num;

      ASSERT_DEBUG_EXT (blk->map = MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Inconsistent sparse storage");
    }
  }

  /* free auxiliary space */
  MEM_Release (&mem);
  free (pcol);
  free (aux);

  return A; 
}

/* copy system matrix into A */
static void matrix_copy (LINSYS *sys, MX *A)
{
#if MPI
  short boundary = sys->variant & BOUNDARY_ONLY;
#endif
  double *a = A->x;
  SET *item;
  DIAB *dia;
  OFFB *blk;
  CON *con;

  for (double *b = a, *c = a + A->nzmax; b < c; b ++) *b = 0.0;

  for (item = SET_First (sys->subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;

    double *T = dia->T;
    int *map = dia->map;

    a [map[0]+0] += T[0];
    a [map[0]+1] += T[1];
    a [map[0]+2] += T[2];
    a [map[1]+0] += T[3];
    a [map[1]+1] += T[4];
    a [map[1]+2] += T[5];
    a [map[2]+0] += T[6];
    a [map[2]+1] += T[7];
    a [map[2]+2] += T[8];

    for (blk = dia->adj; blk; blk = blk->n)
    {
#if MPI
      if (boundary && blk->dia->adjext == NULL) continue;
#endif

      double *T = blk->T;
      int *map = blk->map;

      a [map[0]+0] += T[0];
      a [map[0]+1] += T[1];
      a [map[0]+2] += T[2];
      a [map[1]+0] += T[3];
      a [map[1]+1] += T[4];
      a [map[1]+2] += T[5];
      a [map[2]+0] += T[6];
      a [map[2]+1] += T[7];
      a [map[2]+2] += T[8];
    }
  }
}
#endif

/* gluing constraint stiffness */
static double glue_stiffness (CON *con)
{
  BULK_MATERIAL *a = con->master->mat,
		*b = con->slave->mat;

  return 2.0 / (1.0/a->young + 1.0/b->young);
}

/* update linearization of a non-contact constraint */
static void system_update_noncontact (CON *con, short dynamic, short boundary, double step, double *b)
{
  DIAB *dia;
  OFFB *blk;

  dia = con->dia;

  double *W = dia->W,
	 *U = dia->U,
	 *R = dia->R,
	 *T = dia->T,
	 *RE = RE(con); /* assumed updated */

  switch ((int) con->kind)
  {
  case FIXPNT:
  {
    double *V = dia->V;

    NNCOPY (W, T);

    if (dynamic)
    {
      b [0] = -U[0] - V[0] - RE[0];
      b [1] = -U[1] - V[1] - RE[1];
      b [2] = -U[2] - V[2] - RE[2];
    }
    else
    {
      b [0] = -U[0] - RE[0];
      b [1] = -U[1] - RE[1];
      b [2] = -U[2] - RE[2];
    }

    for (blk = dia->adj; blk; blk = blk->n)
    {
#if MPI
      if (boundary && blk->dia->adjext == NULL) continue;
#endif

      double *W = blk->W,
	     *T = blk->T;

      NNCOPY (W, T);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      double *W = blk->W,
	     *T = blk->T;

      NNCOPY (W, T);
    }
#endif
  }
  break;
  case FIXDIR:
  case VELODIR:
  case RIGLNK:
  {
    double *V = dia->V;

    T [0] = W[0]; /* set diagonal to diag (W) */
    T [4] = W[4];
    T [8] = W[8];
    T [1] = T [2] = T [3] =
    T [5] = T [6] = T [7] = 0;

    b [0] = -R[0] * W[0]; /* keep RT zero: W (-R) = W dR */
    b [1] = -R[1] * W[4];

    if (dynamic)
    {
      if (con->kind == VELODIR) b [2] = -U[2] + VELODIR(con->Z) - RE[2]; /* V(t+h) = W(R+dR) + B; U = WR + B; -U+V(t+h) = WdR */
      else b [2] = -U[2] - V[2] - RE[2]; /* FIXDIR, RIGLNK */
    }
    else
    {
      if (con->kind == VELODIR) b [2] = -U[2] + VELODIR(con->Z) - RE[2];
      else if (con->kind == FIXDIR) b [2] = -U[2] - RE[2]; 
      else /* RIGLNK: see doc/notes.lyx for explanation */
      {
	double C [3];

	ADD (U, RE, C);

	double d = RIGLNK_LEN (con->Z),
	       e = step * step * DOT2 (C, C);

	if (d*d > e)
	{
	  b[2] = -C[2] + (-d + sqrt (d*d - e)) / step;
	}
	else /* not possible to satisfy: do your best */
	{
	  b[2] = -C[2] - d / step;
	}
      }
    }

    for (blk = dia->adj; blk; blk = blk->n)
    {
#if MPI
      if (boundary && blk->dia->adjext == NULL) continue;
#endif

      double *W = blk->W,
	     *T = blk->T;

      T [2] = W[2];
      T [5] = W[5];
      T [8] = W[8];
      T [0] = T [1] = T [3] = T [4] = T [6] = T [7] = 0;
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      double *W = blk->W,
	     *T = blk->T;

      T [2] = W[2];
      T [5] = W[5];
      T [8] = W[8];
      T [0] = T [1] = T [3] = T [4] = T [6] = T [7] = 0;
    }
#endif
  }
  break;
  case GLUEPNT:
  {
    double C = 1.0 / (step * glue_stiffness (con));

    NNCOPY (W, T);
    T [0] += C;
    T [4] += C;
    T [8] += C;

    b [0] = -C*R[0] -U[0] - RE[0];
    b [1] = -C*R[1] -U[1] - RE[1];
    b [2] = -C*R[2] -U[2] - RE[2];

    for (blk = dia->adj; blk; blk = blk->n)
    {
#if MPI
      if (boundary && blk->dia->adjext == NULL) continue;
#endif
      double *W = blk->W,
	     *T = blk->T;

      NNCOPY (W, T);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      double *W = blk->W,
	     *T = blk->T;

      NNCOPY (W, T);
    }
#endif
  }
  break;
  }
}

/* update linear system for NONSMOOTH_HSW, NONSMOOTH_HYBRID, FIXED_POINT variants */
static void system_update_HSW_HYBRID_FIXED (LINVAR variant, DOM *dom, SET *subset, double *rhs)
{
  double d [3], norm, lim, udash, step, *b;
  short dynamic, pull, boundary;
  DIAB *dia;
  OFFB *blk;
  SET *item;
  CON *con;

  boundary = variant & BOUNDARY_ONLY;
  dynamic = dom->dynamic;
  step = dom->step;

#if MPI
  DOM_Update_External_Reactions (dom, 0); /* (###) */
#endif

  for (item = SET_First (subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;
    b = &rhs [3*con->num];

    double *W = dia->W,
	   *B = dia->B,
	   *U = dia->U,
	   *R = dia->R,
	   *T = dia->T,
	   *RE = RE(con);

    /* update residual */
    NVADDMUL (B, W, R, RE);
    SUB (RE, U, RE);
    for (blk = dia->adj; blk; blk = blk->n)
    {
      double *R = blk->dia->R,
             *W = blk->W;

      NVADDMUL (RE, W, R, RE);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      double *R = CON(blk->dia)->R, /* (###) */
             *W = blk->W;

      NVADDMUL (RE, W, R, RE);
    }
#endif

    switch (con->kind)
    {
    case CONTACT:
    {
      double *V = dia->V,
	     rho = dia->rho,
	     gap = con->gap,
	     fri = con->mat.base->friction,
	     res = con->mat.base->restitution;

      /* normal dashed velocity */
      if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
      else udash = ((MAX(gap, 0)/step) + U[2]);

      /* predict reaction */
      d [0] = R[0] - rho * U[0];
      d [1] = R[1] - rho * U[1];
      d [2] = R[2] - rho * udash;

      if (d [2] >= 0)  /* active normal set */
      {
	pull = 0;

	T [2] = W[2];
	T [5] = W[5];
	T [8] = W[8];

	b [2] = -udash - RE[2];

	for (blk = dia->adj; blk; blk = blk->n)
	{
#if MPI
	  if (boundary && blk->dia->adjext == NULL) continue;
#endif

	  double *W = blk->W,
		 *T = blk->T;

	  T [2] = W[2];
	  T [5] = W[5];
	  T [8] = W[8];
	}
#if MPI
	for (blk = dia->adjext; blk; blk = blk->n)
	{
	  double *W = blk->W,
		 *T = blk->T;

	  T [2] = W[2];
	  T [5] = W[5];
	  T [8] = W[8];
	}
#endif
      }
      else  /* inactive normal set */
      {
	pull = 1;

	T [2] = 0;
	T [5] = 0;
	T [8] = W[8];

	b [2] = -R[2] * W[8];

	for (blk = dia->adj; blk; blk = blk->n)
	{
#if MPI
	  if (boundary && blk->dia->adjext == NULL) continue;
#endif

	  double *T = blk->T;

	  T [2] = T [5] = T [8] = 0;
	}
#if MPI
	for (blk = dia->adjext; blk; blk = blk->n)
	{
	  double *T = blk->T;

	  T [2] = T [5] = T [8] = 0;
	}
#endif
      }

      /* tangential response */

      norm = sqrt (d[0]*d[0]+d[1]*d[1]); /* tangential force value */

      switch (LINEARIZATION_VARIANT (variant))
      {
	case NONSMOOTH_HSW:
	  lim = fri * MAX (0, d[2]);
	  break;
	case NONSMOOTH_HYBRID:
	  lim = fri * MAX (0, R[2]);
	  break;
	case FIXED_POINT:
	  lim = fri * MAX (0, RN(con));
	  break;
      }

      if ((variant & NONSMOOTH_HSW) && pull) goto ZERO_TANG; /* enforce AN = AT + IT */

      if (norm >= lim) /* frictional sliping */
      {
	double F [4], /* matrix associated with the derivative of an Euclidean norm in 2D */
	       M [4], H [4],  /* auxiliary metrices & vectors */
	       delta, alpha, beta, den, len, e; /* auxiliary scalars */

	/* active tangential set */

	if (lim > 0.0 && norm > 0.0) /* non-degenerate case */
	{
	  /* references below are linked to the paper "A primal-dual active set strategy
	   * for three-dimensional contact problems with Coulomb friction", S. Hueber, G. Stadler,
	   * B.I. Wohlmuth, SIAM Journal on Scientific Computing 30(2):572-596, 2007*/

	  len = sqrt (R[0]*R[0]+R[1]*R[1]);
	  den = MAX (lim, len) * norm;
	  e = lim / norm; /* (3.6) */
	  
#if 0
	  if (len == 0.0) beta = 1.0; /* commented in Section 3.3, although my conclusion is to set
					 'beta = 1' here rather than enforce homogenous tangential reaction */
	  else
	  {
	    /* formulae around (3.11) below */
	    alpha = (R[0]*d[0]+R[1]*d[1]) / (len*norm);
	    delta = MIN (len/lim, 1.0);
	    beta = (alpha < 0.0 ? 1.0 / (1.0 - alpha*delta) : 1.0); /* relaxation factor in case of direction change */
	  }
#else
	  alpha = delta = 0.0;
	  beta = 1.0;
#endif

	  /* (3.11) {or modified (3.6)} */
	  F [0] = (R[0]*d[0])/den;
	  F [1] = (R[1]*d[0])/den;
	  F [2] = (R[0]*d[1])/den;
	  F [3] = (R[1]*d[1])/den;

	  /* defined just above (3.10) */
	  M [0] = e * (1.0 - F[0]);
	  M [1] = - e * F[1];
	  M [2] = - e * F[2];
	  M [3] = e * (1.0 - F[3]);

	  /* inverse used in (3.10) */
	  H [0] = 1.0 - beta * M[0];
	  H [1] = - beta * M[1];
	  H [2] = - beta * M[2];
	  H [3] = 1.0 - beta * M[3];

	  /* scale */
	  M [0] *= W[0];
	  M [1] *= W[4];
	  M [2] *= W[0];
	  M [3] *= W[4];

	  H [0] *= W[0];
	  H [1] *= W[4];
	  H [2] *= W[0];
	  H [3] *= W[4];

	  T [0] = H[0] + rho*(M[0]*W[0] + M[2]*W[1]);
	  T [1] = H[1] + rho*(M[1]*W[0] + M[3]*W[1]);
	  T [3] = H[2] + rho*(M[0]*W[3] + M[2]*W[4]);
	  T [4] = H[3] + rho*(M[1]*W[3] + M[3]*W[4]);

	  if (variant & FIXED_POINT)
	  {
	    T [6] = rho*(M[0]*W[6] + M[2]*W[7]);
	    T [7] = rho*(M[1]*W[6] + M[3]*W[7]);

	    b [0] = (fri*(d[0]/norm)*RN(con) - R[0]) * W[0] - rho*(M[0]*RE[0] + M[2]*RE[1]); 
	    b [1] = (fri*(d[1]/norm)*RN(con) - R[1]) * W[4] - rho*(M[1]*RE[0] + M[3]*RE[1]); 
	  }
	  else /* HSW, HYBRID */
	  {
	    T [6] = rho*(M[0]*W[6] + M[2]*W[7]) - fri*(d[0]/norm) * W[0];
	    T [7] = rho*(M[1]*W[6] + M[3]*W[7]) - fri*(d[1]/norm) * W[4];

	    b [0] = (fri*(d[0]/norm)*R[2] - R[0]) * W[0] - rho*(M[0]*RE[0] + M[2]*RE[1]);
	    b [1] = (fri*(d[1]/norm)*R[2] - R[1]) * W[4] - rho*(M[1]*RE[0] + M[3]*RE[1]);
	  }

	  for (blk = dia->adj; blk; blk = blk->n)
	  {
#if MPI
	    if (boundary && blk->dia->adjext == NULL) continue;
#endif

	    double *W = blk->W,
		   *T = blk->T;

	    T[0] = rho*(M[0]*W[0] + M[2]*W[1]);
	    T[1] = rho*(M[1]*W[0] + M[3]*W[1]);
	    T[3] = rho*(M[0]*W[3] + M[2]*W[4]);
	    T[4] = rho*(M[1]*W[3] + M[3]*W[4]);
	    T[6] = rho*(M[0]*W[6] + M[2]*W[7]); 
	    T[7] = rho*(M[1]*W[6] + M[3]*W[7]);  
	  }
#if MPI
	  for (blk = dia->adjext; blk; blk = blk->n)
	  {
	    double *W = blk->W,
		   *T = blk->T;

	    T[0] = rho*(M[0]*W[0] + M[2]*W[1]);
	    T[1] = rho*(M[1]*W[0] + M[3]*W[1]);
	    T[3] = rho*(M[0]*W[3] + M[2]*W[4]);
	    T[4] = rho*(M[1]*W[3] + M[3]*W[4]);
	    T[6] = rho*(M[0]*W[6] + M[2]*W[7]); 
	    T[7] = rho*(M[1]*W[6] + M[3]*W[7]);  
	  }
#endif
	}
	else /* degenerate case => enforce homogenous tangential tractions.
		this is discussed in the Section 3.3 of the HSW paper */
	{
ZERO_TANG:
	  T [0] = W [0];
	  T [4] = W [4];
	  T [1] = T [2] = T [3] = 
	  T [5] = T [6] = T [7] = 0;

	  b [0] = -R[0] * W[0];
	  b [1] = -R[1] * W[4];

	  for (blk = dia->adj; blk; blk = blk->n)
	  {
#if MPI
	    if (boundary && blk->dia->adjext == NULL) continue;
#endif

	    double *T = blk->T;

	    T [0] = T [1] = T [3] = T [4] = T [6] = T [7] = 0;
	  }
#if MPI
	  for (blk = dia->adjext; blk; blk = blk->n)
	  {
	    double *T = blk->T;

	    T [0] = T [1] = T [3] = T [4] = T [6] = T [7] = 0;
	  }
#endif
	}
      }
      else /* frictional sticking: inactive tangential set */
      {

	T [0] = W[0];
	T [1] = W[1];
	T [3] = W[3];
	T [4] = W[4];

	if ((variant & NONSMOOTH_HYBRID) || (variant & FIXED_POINT))
	{
	  T [6] = W[6];
	  T [7] = W[7];

	  if (dynamic)
	  {
	    b [0] = -U[0] - V[0] - RE[0]; /* -V = W(R + dR) + B; U = WR + B; -U - V = W dR */
	    b [1] = -U[1] - V[1] - RE[1];
	  }
	  else
	  {
	    b [0] = -U[0] - RE[0];
	    b [1] = -U[1] - RE[1];
	  }
	}
	else /* HSW */
	{
	  /* FIXME: this is quasi-statics; work out dynamics */

	  T [6] = W[6]+U[0]/d[2];
	  T [7] = W[7]+U[1]/d[2];

	  b [0] = -(1.0 + rho*udash/d[2])*U[0] - RE[0];
	  b [1] = -(1.0 + rho*udash/d[2])*U[1] - RE[1];
	}

	for (blk = dia->adj; blk; blk = blk->n)
	{
#if MPI
	  if (boundary && blk->dia->adjext == NULL) continue;
#endif

	  double *W = blk->W,
		 *T = blk->T;

	  T [0] = W[0];
	  T [1] = W[1];
	  T [3] = W[3];
	  T [4] = W[4];
	  T [6] = W[6];
	  T [7] = W[7];
	}
#if MPI
	for (blk = dia->adjext; blk; blk = blk->n)
	{
	  double *W = blk->W,
		 *T = blk->T;

	  T [0] = W[0];
	  T [1] = W[1];
	  T [3] = W[3];
	  T [4] = W[4];
	  T [6] = W[6];
	  T [7] = W[7];
	}
#endif
      }
    }
    break;
    default:
    {
      system_update_noncontact (con, dynamic, boundary, step, b);
    }
    break;
    }
  }
}

/* real normal to friction cone */
inline static void real_n (double *S, double fri, double *n)
{
  double dot, len;

  dot = DOT2(S, S);
  len = sqrt(dot);

  if (len == 0 || len <= fri * S[2])
  {
    SET (n, 0.0);
  }
  else if (fri * len + S[2] < 0.0)
  {
    dot += S[2]*S[2];
    len = sqrt (dot);
    if (len == 0) { SET (n, 0.0); }
    else { DIV (S, len, n); }
  }
  else
  {
    dot = 1.0 / sqrt (1.0 + fri*fri);
    DIV2 (S, len, n);
    n [2] = -fri;
    SCALE (n, dot);
  }
}

/* imaginary i */
static double complex imaginary_i;

/* complex normal to friction cone */
inline static void complex_n (double complex *S, double complex fri, double complex *n)
{
  double complex dot, len;

  dot = DOT2(S, S);
  len = csqrt(dot);

  if (creal (len) == 0 || creal (len) <= creal (fri * S[2]))
  {
    SET (n, 0.0 + 0.0 * imaginary_i);
  }
  else if (creal (fri * len + S[2]) < 0.0)
  {
    dot += S[2]*S[2];
    len = csqrt (dot);
    if (creal (len) == 0) { SET(n, 0.0 + 0.0 * imaginary_i); }
    DIV (S, len, n);
  }
  else
  {
    dot = 1.0 / csqrt (1.0 + fri*fri);
    DIV2 (S, len, n);
    n [2] = -fri;
    SCALE (n, dot);
  }
}

/* real normal ray to friction cone */
inline static void real_m (double fri, short smooth, double *S, double eps, double *m)
{
  double n [3], fun;

  real_n (S, fri, n);
  fun = DOT (S, n);

  if (smooth == 1 && fun >= 0.0 && fun <= eps)
  {
    fun = ((2.0/eps) - (1.0/(eps*eps))*fun)*(fun*fun);
  }
  else if (smooth == 2)
  {
    if (fun >= 0.0 && fun <= eps)
    {
      fun = (fun*fun) / (2.0 * eps);
    }
    else if (fun > eps)
    {
      fun = fun - 0.5 * eps;
    }
  }

  MUL (n, fun, m)
}

/* complex normal ray to friction cone */
inline static void complex_m (double complex fri, short smooth, double complex *S, double complex eps, double complex *m)
{
  double complex n [3], fun;

  complex_n (S, fri, n);
  fun = DOT (S, n);

  if (smooth == 1 && creal (fun) >= 0.0 && creal (fun) <= creal (eps))
  {
    fun = ((2.0/eps) - (1.0/(eps*eps))*fun)*(fun*fun);
  }
  else if (smooth == 2)
  {
    if (creal (fun) >= 0.0 && creal (fun) <= creal (eps))
    {
      fun = (fun*fun) / (2.0 * eps);
    }
    else if (creal (fun) > creal (eps))
    {
      fun = fun - 0.5 * eps;
    }
  }

  MUL (n, fun, m)
}

/* real F = [UT, UN + fri |UT|]' */
inline static void real_F (double res, double fri, double gap, double step, short dynamic, double *V, double *U, double *F)
{
  double udash;

  if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
  else udash = ((MAX(gap, 0)/step) + U[2]);

  F [0] = U[0];
  F [1] = U[1];
  F [2] = (udash + fri * sqrt (DOT2(U, U)));
}
 
/* complex F = [UT, UN + fri |UT|]' */
inline static void complex_F (double res, double fri, double gap, double step, short dynamic, double *V, double complex *U, double complex *F)
{
  double complex udash;

  if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
  else udash = ((MAX(gap, 0)/step) + U[2]);

  F [0] = U[0];
  F [1] = U[1];
  F [2] = (udash + fri * csqrt (DOT2(U, U)));
}

/* update linear system for NONSMOOTH_VARIATIONAL, SMOOTHED_VARIATIONAL variants */
static void system_update_VARIATIONAL (LINVAR variant, double epsilon, DOM *dom, SET *subset, double *rhs)
{
  short smooth, dynamic, boundary;
  double step, h, *b;
  SET *item;
  OFFB *blk;
  DIAB *dia;
  CON *con;

  if (variant & SMOOTHED_VARIATIONAL) smooth = 1; /* (+++); TODO: decide on best smoothing */
  h = 1E-6 * epsilon; /* TODO: decide on best choice */
  boundary = variant & BOUNDARY_ONLY;
  dynamic = dom->dynamic;
  step = dom->step;

#if MPI
  DOM_Update_External_Reactions (dom, 0); /* (###) */
#endif

  for (item = SET_First (subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;
    b = &rhs [3*con->num];

    double *W = dia->W,
	   *B = dia->B,
	   *U = dia->U,
	   *R = dia->R,
	   *RE = RE(con);

    /* update residual */
    NVADDMUL (B, W, R, RE);
    SUB (RE, U, RE);
    for (blk = dia->adj; blk; blk = blk->n)
    {
      double *R = blk->dia->R,
             *W = blk->W;

      NVADDMUL (RE, W, R, RE);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      double *R = CON(blk->dia)->R, /* (###) */
             *W = blk->W;

      NVADDMUL (RE, W, R, RE);
    }
#endif

    switch (con->kind)
    {
    case CONTACT:
    {
      double *V = dia->V,
	     gap = con->gap,
	     fri = con->mat.base->friction,
	     res = con->mat.base->restitution,
	     *T  = dia->T,
	     TT [9],
	     dm [9],
	     dF [9],
             S [3],
	     F [3],
	     m [3],
	     H [3],
	     J [9];


      real_F (res, fri, gap, step, dynamic, V, U, F);
      SUB (R, F, S);
      real_m (fri, smooth, S, epsilon, m);
      ADD (F, m, H);

      double complex cU [3],
	             cS [3],
		     cF [3],
		     cm [3];

      for (int k = 0; k < 3; k ++)
      {
	cU [0] = U[0] + 0.0 * imaginary_i;
	cU [1] = U[1] + 0.0 * imaginary_i;
	cU [2] = U[2] + 0.0 * imaginary_i;
	cU [k] += h * imaginary_i;
        complex_F (res, fri, gap, step, dynamic, V, cU, cF);
	dF [3*k+0] = cimag (cF [0]) / h;
	dF [3*k+1] = cimag (cF [1]) / h;
	dF [3*k+2] = cimag (cF [2]) / h;

        cS [0] = S[0] + 0.0 * imaginary_i;
	cS [1] = S[1] + 0.0 * imaginary_i;
	cS [2] = S[2] + 0.0 * imaginary_i;
	cS [k] += h * imaginary_i;
        complex_m (fri, smooth, cS, epsilon, cm);
	dm [3*k+0] = cimag (cm [0]) / h;
	dm [3*k+1] = cimag (cm [1]) / h;
	dm [3*k+2] = cimag (cm [2]) / h;
      }

      IDENTITY (J);
      NNSUB (J, dm, J);
      NNMUL (dF, J, TT); /* dF/dU [I - dm/dS] */

      NVADDMUL (H, TT, RE, b);
      SCALE (b, -1.0); /* b = - [H(U, R) + dF/dU [I - dm/dS] RE] */

      NNMUL (TT, W, T);
      NNADD (T, dm, T);

      for (blk = dia->adj; blk; blk = blk->n)
      {
#if MPI
	if (boundary && blk->dia->adjext == NULL) continue;
#endif
	double *T = blk->T,
	       *W = blk->W;

        NNMUL (TT, W, T);
      }
#if MPI
      for (blk = dia->adjext; blk; blk = blk->n)
      {
	double *T = blk->T,
	       *W = blk->W;

        NNMUL (TT, W, T);
      }
#endif
    }
    break;
    default:
    {
      system_update_noncontact (con, dynamic, boundary, step, b);
    }
    break;
    }
  }
}

/* update noraml bounds and return relative error of the update */
static double update_normal_bounds (LINSYS *sys)
{
  double errup, errlo;
  SET *item;
  CON *con;

  errup = errlo = 0.0;

  for (item = SET_First (sys->subset); item; item = SET_Next (item))
  {
    con = item->data;
    if (con->kind == CONTACT)
    {
      double DRN,
	     RN;

      DRN = con->R[2] - RN(con);
      RN = RN(con) = con->R[2];

      errup += DRN*DRN;
      errlo += RN*RN;
    }
  }

#if MPI
  double errloc [2] = {errup, errlo}, errsum [2];
  MPI_Allreduce (errloc, errsum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  errup = errsum [0], errlo = errsum [1];
#endif

  return sqrt (errup / MAX (errlo, 1.0));
}

#if MPI
static void update_external_reactions (LINSYS *sys, double *x)
{
  SET *item, *jtem;
  double **r, *z;
  int i, *j, *k;
  COMDATA *d;
  CON *con;

  for (item = SET_First (sys->subset), d = sys->send; item; item = SET_Next (item))
  {
    con = item->data;
    for (jtem = SET_First (con->ext); jtem; jtem = SET_Next (jtem), d ++)
    {
      d->d = &x [3*con->num];
    }
  }

  COMALL_Repeat (sys->pattern);

  for (i = 0, d = sys->recv, r = sys->R; i < sys->nrecv; i ++, d ++)
  {
    for (j = d->i, k = d->i + d->ints, z = d->d; j < k; j ++, r ++, z += 3)
    {
      COPY (z, *r);
    }
  }
}
#endif

static void A_times_x_equals_y (LINSYS *sys, double *x, double *y)
{
  double *T, *z;
  short boundary;
  SET *item;
  CON *con;
  OFFB *blk;
  DIAB *dia;

#if MPI
  update_external_reactions (sys, x); /* (###) */
#endif

  boundary = sys->variant & BOUNDARY_ONLY;

  for (item = SET_First (sys->subset); item; item = SET_Next (item), y += 3)
  {
    con = item->data;
    dia = con->dia;
    T = dia->T;
    z = &x [3*con->num];

    NVMUL (T, z, y);

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
#if MPI
      if (boundary && blk->dia->adjext == NULL) continue;
#endif
      T = blk->T;
      z = &x [3*con->num];

      NVADDMUL (y, T, z, y);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = CON (blk->dia);
      T = blk->T;
      z = con->R; /* (###) */

      NVADDMUL (y, T, z, y);
    }
#endif
  }
}

#if USE_AT
static void AT_times_x_equals_y (LINSYS *sys, double *x, double *y)
{
#if MPI
  int *basenum, ncpu, rank;
  double *v;
#endif
  double *T, *z, *w;
  short boundary;
  SET *item;
  CON *con;
  OFFB *blk;
  DIAB *dia;
  DOM *dom;

  dom = sys->ldy->dom;
  boundary = sys->variant & BOUNDARY_ONLY;
#if MPI
  ncpu = dom->ncpu;
  rank = dom->rank;
  basenum = sys->basenum;
  ERRMEM (v = MEM_CALLOC (basenum [ncpu] * sizeof (double))); /* global result vector */
#else
  for (z = y, w = z + sys->b->n; z < w; z ++) (*z) = 0.0; /* zero y */
#endif

  for (item = SET_First (sys->subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;
    T = dia->T;
    z = &x [3*con->num];
#if MPI
    w = &v [3*(basenum [rank]+con->num)];
#else
    w = &y [3*con->num];
#endif

    TVADDMUL (w, T, z, w);

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
#if MPI
      if (boundary && blk->dia->adjext == NULL) continue;
#endif
      T = blk->T;
#if MPI
      w = &v [3*(basenum [rank]+con->num)];
#else
      w = &y [3*con->num];
#endif
      TVADDMUL (w, T, z, w);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = CON (blk->dia);
      T = blk->T;
      w = &v [3*(basenum [con->rank]+con->num)];

      TVADDMUL (w, T, z, w);
    }
#endif
  }

#if MPI
  MPI_Allreduce (MPI_IN_PLACE, v, 3*basenum [ncpu], MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (w = &v [3*basenum [rank]], z = &v [3*basenum [rank+1]]; w < z; w ++, y ++) (*y) = (*w);
  free (v);
#endif
}
#endif

static struct vect* newvect (int n)
{
  struct vect *v;

  ERRMEM (v = malloc (sizeof (struct vect)));
  ERRMEM (v->x = MEM_CALLOC (n * sizeof (double)));
  v->n = n;

  return v;
}

/* GMRES interface start */
static char* CAlloc (size_t count, size_t elt_size)
{
  char *ptr;

  ERRMEM (ptr = MEM_CALLOC (count * elt_size));
  return ptr;
}

static int Free (char *ptr)
{
  free (ptr);
  return 0;
}

static int CommInfo (void *A, int *my_id, int *num_procs)
{
#if MPI
  LINSYS *sys = (LINSYS*)A;
  DOM *dom = sys->ldy->dom;

  *my_id = dom->rank;
  *num_procs = dom->ncpu;
#else
  *my_id = 0;
  *num_procs = 1;
#endif
  return 0;
}

static void* CreateVector (void *vector)
{
  struct vect *a = vect (vector), *v;

  ERRMEM (v = malloc (sizeof (struct vect)));
  ERRMEM (v->x = MEM_CALLOC (a->n * sizeof (double)));
  v->n = a->n;

  return v;
}

static void* CreateVectorArray (int size, void *vectors)
{
  struct vect **v;
  int i;

  ERRMEM (v = malloc (size * sizeof (struct vect*)));
  for (i = 0; i < size; i ++)
  {
    v[i] = CreateVector (vectors);
  }

  return v;
}

static int DestroyVector (void *vector)
{
  struct vect *v = vect (vector);
  free (v->x);
  free (v);

  return 0;
}

static double InnerProd (void *vx, void *vy)
{
  double dot = 0.0, *x, *y, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n, y = vect (vy)->x; x < z; x ++, y ++)
  {
    dot += (*x) * (*y);
  }

#if MPI
  double val = dot;
  MPI_Allreduce (&val, &dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  return dot;
}

static int CopyVector (void *vx, void *vy)
{
  double *x, *y, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n, y = vect (vy)->x; x < z; x ++, y ++)
  {
    (*y) = (*x);
  }

  return 0;
}

static int ClearVector (void *vx)
{
  double *x, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n; x < z; x ++)
  {
    (*x) = 0.0;
  }

  return 0;
}

static int ScaleVector (double alpha, void *vx)
{
  double *x, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n; x < z; x ++)
  {
    (*x) *= alpha;
  }

  return 0;
}

static int  Axpy (double alpha, void *vx, void *vy )
{
  double *x, *y, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n, y = vect (vy)->x; x < z; x ++, y ++)
  {
    (*y) += alpha * (*x);
  }

  return 0;
}

static void *MatvecCreate (void *A, void *x)
{
  return NULL;
}

static int Matvec (void *matvec_data, double alpha, void *A, void *vx, double beta, void *vy)
{
  struct vect *vz = CreateVector (vx), *vu = CreateVector (vx);
  LINSYS *sys = (LINSYS*)A;
  double delta = sys->delta;

  ScaleVector (beta, vy);
  Axpy (delta, vx, vy);
  A_times_x_equals_y (sys, vect (vx)->x, vz->x);
#if USE_AT
  AT_times_x_equals_y (sys, vz->x, vu->x);
  Axpy (alpha, vu, vy);
#else
  Axpy (alpha, vz, vy);
#endif

  DestroyVector (vz);
  DestroyVector (vu);

  return 0;
}

static int MatvecDestroy (void *matvec_data)
{
  return 0;
}

static int PrecondSetup (void *vdata, void *A, void *vb, void *vx)
{
  return 0;
}

static int Precond (void *vdata, void *A, void *vb, void *vx)
{
#if 1
  LINSYS *sys = (LINSYS*)A;
  double *bb = vect (vb)->x, *xx = vect (vx)->x, *b, *x;
  int ipiv [3];
  DIAB *dia;
  SET *item;
  CON *con;

  for (item = SET_First (sys->subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;
    b = &bb [3*con->num];
    x = &xx [3*con->num];
    COPY (b, x);

    double *T = dia->T, S [9];

    NNCOPY (T, S);
    lapack_dgesv (3, 1, S, 3, ipiv, x, 3); /* TODO: inv (S) could be precomputed at the start */
  }

  return 0;
#endif
  return CopyVector (vb, vx);
}
/* GMRES interface end */

/* create linear system resulting from linearization of constraints */
LINSYS* LINSYS_Create (LINVAR variant, LOCDYN *ldy)
{
  DOM *dom = ldy->dom;
  short boundary;
  LINSYS *sys;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  int ncon;

  ncon = dom->ncon;
  imaginary_i = csqrt (-1);
  boundary = variant & BOUNDARY_ONLY;
  ERRMEM (sys = MEM_CALLOC (sizeof (LINSYS)));
  MEM_Init (&sys->setmem, sizeof (SET), MAX (ncon, BLOCKS));
  sys->variant = variant;
  sys->ldy = ldy;

  /* collect constraints */
  for (con = dom->con; con; con = con->next)
  {
#if MPI
    if (boundary)
    {
      if (con->dia->adjext)
      {
	SET_Insert (&sys->setmem, &sys->subset, con, NULL);
      }
    }
    else 
#endif
    SET_Insert (&sys->setmem, &sys->subset, con, NULL);
  }
  ncon = SET_Size (sys->subset);

  /* constraints local numbering */
  DOM_Number_Constraints (dom, 1, sys->subset);

  /* auxiliary memory pools */
  switch (LINEARIZATION_VARIANT (variant))
  {
    case NONSMOOTH_HSW:
    case NONSMOOTH_HYBRID:
    case NONSMOOTH_VARIATIONAL:
    case SMOOTHED_VARIATIONAL:
      MEM_Init (&sys->Y_mem, sizeof (double [9]), MAX (ncon, BLOCKS)); /* DR, DU, RE */
      break;
    case FIXED_POINT:
      MEM_Init (&sys->Y_mem, sizeof (double [10]), MAX (ncon, BLOCKS)); /* DR, DU, RE, RN */
      break;
      MEM_Init (&sys->Y_mem, sizeof (double [9]), MAX (ncon, BLOCKS)); /* DR, DU, RE */
      break;
  }
  MEM_Init (&sys->T_mem, sizeof (double [9]),  MAX (ncon, BLOCKS)); /* dHi/dRj */

  /* auxiliary constraint space */
  for (con = dom->con; con; con = con->next)
  {
#if MPI
    if (boundary && con->dia->adjext == NULL) continue;
#endif

    dia = con->dia;
    ERRMEM (con->Y = MEM_Alloc (&sys->Y_mem));
    ERRMEM (dia->T = MEM_Alloc (&sys->T_mem));

    for (blk = dia->adj; blk; blk = blk->n)
    {
#if MPI
      if (boundary && blk->dia->adjext == NULL) continue;
#endif
      ERRMEM (blk->T = MEM_Alloc (&sys->T_mem));
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      ERRMEM (blk->T = MEM_Alloc (&sys->T_mem));
    }
#endif
  }
#if MPI
  for (MAP *item = MAP_First (dom->conext); item; item = MAP_Next (item))
  {
    con = item->data;
    ERRMEM (con->Y = MEM_Alloc (&sys->Y_mem));
  }
#endif

  /* default fixed point tolerance */
  sys->fixed_point_tol = 1E-6;

  /* smoothing thickness */
  if (variant & SMOOTHED_VARIATIONAL)
  {
    double B [4] = {0, 0, 0, 0};

    for (con = dom->con; con; con = con->next)
    {
      if (con->kind == CONTACT)
      {
        dia = con->dia;
	ACC (dia->B, B);
	B [3] ++;
      }
    }

#if MPI
    double C [4];
    MPI_Allreduce (B, C, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    COPY4 (C, B);
#endif

    B [3] = 1.0 / MAX (B [3], 1.0);
    SCALE (B, B[3]); /* avergae free contact velocity */

    sys->epsilon = 1E-6 * LEN (B); /* TODO: work out a best choice */
  }
  else sys->epsilon = 1E-9;

  /* unknown and right hand size */
  sys->x = newvect (3 * ncon);
  sys->b = CreateVector (sys->x);

  /* create GMRES solver */
  sys->gmres_functions = hypre_FlexGMRESFunctionsCreate (CAlloc, Free, CommInfo, CreateVector, CreateVectorArray, DestroyVector,
    MatvecCreate, Matvec, MatvecDestroy, InnerProd, CopyVector, ClearVector, ScaleVector, Axpy, PrecondSetup, Precond);
  sys->gmres_vdata = hypre_FlexGMRESCreate (sys->gmres_functions);

#if MPI
  /* communication pattern */
  SET *item, *jtem;
  int i, *j, *k, n;
  COMDATA *d;
  double **r;

  for (item = SET_First (sys->subset); item; item = SET_Next (item))
  {
    con = item->data;
    for (jtem = SET_First (con->ext); jtem; jtem = SET_Next (jtem))
    {
      sys->nsend ++;
    }
  }

  ERRMEM (sys->send = MEM_CALLOC (sys->nsend * sizeof (COMDATA)));
  d = sys->send;

  for (item = SET_First (sys->subset); item; item = SET_Next (item))
  {
    con = item->data;
    for (jtem = SET_First (con->ext); jtem; jtem = SET_Next (jtem), d ++)
    {
      d->rank = (int) (long) jtem->data;
      d->doubles = 3;
      d->ints = 1;
      d->d = &vect (sys->x)->x [3*con->num];
      d->i = (int*) &con->id;
    }
  }

  sys->pattern = COMALL_Pattern (MPI_COMM_WORLD, sys->send, sys->nsend, &sys->recv, &sys->nrecv);

  COMALL_Repeat (sys->pattern);

  /* map received reaction pointers */
  for (i = n = 0, d = sys->recv; i < sys->nrecv; i ++, d ++)
  {
    for (j = d->i, k = j + d->ints; j < k; j ++)
    {
      n ++;
    }
  }

  ERRMEM (sys->R = MEM_CALLOC (n * sizeof (double*)));
  r = sys->R;

  for (i = 0, d = sys->recv; i < sys->nrecv; i ++, d ++)
  {
    for (j = d->i, k = j + d->ints; j < k; j ++, r ++)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*j), NULL), "Invalid constraint id");
      (*r) = con->R;
    }
  }

  /* initial constraint numbers */
  ERRMEM (sys->basenum = MEM_CALLOC ((dom->ncpu + 1) * sizeof (int)));

  MPI_Allgather  (&ncon, 1, MPI_INT, (sys->basenum + 1), 1, MPI_INT, MPI_COMM_WORLD);

  for (i = 1; i <= dom->ncpu; i ++) sys->basenum [i] += sys->basenum [i-1];
#endif

#if USE_LU
  sys->A = matrix_create (sys); /* used for preconditioning */
#endif

  return sys;
}

/* set fixed point approach normal stress update error tolerance */
void LINSYS_Fixed_Point_Tol (LINSYS *sys, double tol)
{
  sys->fixed_point_tol = tol;
}

/* update linear system at current reactions R */
void LINSYS_Update (LINSYS *sys)
{
  switch (LINEARIZATION_VARIANT (sys->variant))
  {
    case NONSMOOTH_HSW:
    case NONSMOOTH_HYBRID:
    case FIXED_POINT:
      system_update_HSW_HYBRID_FIXED (sys->variant, sys->ldy->dom, sys->subset, sys->b->x);
      break;
    case NONSMOOTH_VARIATIONAL:
    case SMOOTHED_VARIATIONAL:
      system_update_VARIATIONAL (sys->variant, sys->epsilon, sys->ldy->dom, sys->subset, sys->b->x);
      break;
  }

#if USE_AT
  /* b = (A^T) b */
  struct vect *c = CreateVector (sys->b);
  CopyVector (sys->b, c);
  AT_times_x_equals_y (sys, c, sys->b);
  DestroyVector (c);
#endif
}

/* solve for reaction increments DR */
void LINSYS_Solve (LINSYS *sys, double abstol, int maxiter)
{
  short boundary, variational;
  double *x, *z, *DR;
  DIAB *dia;
  SET *item;
  CON *con;

#if 0
  sys->delta = 0.1 * abstol; /* FIXME: decide on more inteligent choice */
  abstol -= sys->delta;
#endif

  boundary = sys->variant & BOUNDARY_ONLY;
  variational = sys->variant & (SMOOTHED_VARIATIONAL|NONSMOOTH_VARIATIONAL);

#if USE_LU
  double *xx = sys->x->x, *bb = sys->b->x;

  for (x = xx + sys->x->n; xx < x; xx ++, bb ++) *xx = *bb;

  matrix_copy (sys, sys->A);

  cs_lusol (1, sys->A, sys->x->x, 0);
#else
  hypre_FlexGMRESSetTol (sys->gmres_vdata, 0.0);
  hypre_FlexGMRESSetMaxIter (sys->gmres_vdata, maxiter);
  hypre_FlexGMRESSetAbsoluteTol (sys->gmres_vdata, abstol);
  hypre_FlexGMRESSetup (sys->gmres_vdata, sys, sys->b, sys->x);
  hypre_FlexGMRESSolve (sys->gmres_vdata, sys, sys->b, sys->x);
  hypre_FlexGMRESGetNumIterations (sys->gmres_vdata , &sys->iters);
  hypre_FlexGMRESGetFinalRelativeResidualNorm (sys->gmres_vdata, &sys->resnorm);
#endif

  for (item = SET_First (sys->subset), x = sys->x->x; item; item = SET_Next (item))
  {
    con = item->data;
    z = &x [3*con->num];
    DR = DR (con);
    COPY (z, DR);

    if (variational && con->kind == CONTACT) /* DR = proj (friction-cone, R+DR) - R */
    {
      double *R = con->R,
	      fri = con->mat.base->friction,
	      S [3],
	      m [3];

      ADD (R, DR, S);
      real_m (fri, 0, S, 0, m);
      SUB (S, m, S);
      SUB (S, R, DR);
    }
  }

#if MPI
  DOM_Update_External_Y (sys->ldy->dom, 3, boundary ? sys->subset : NULL); /* (&&&) */
#endif

  /* DU  */
  for (item = SET_First (sys->subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;

    double *DR = DR(con),
	   *DU = DU(con),
	   *W = dia->W;

    NVMUL (W, DR, DU);

    for (OFFB *blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;

#if MPI
      if (boundary && blk->dia->adjext == NULL) continue;
#endif

      double *DR = DR(con),
	     *W = blk->W;

      NVADDMUL (DU, W, DR, DU);
    }
#if MPI
    for (OFFB *blk = dia->adjext; blk; blk = blk->n)
    {
      con = CON (blk->dia);

      double *DR = DR (con), /* (&&&) */
	     *W = blk->W;

      NVADDMUL (DU, W, DR, DU);
    }
#endif
  }
}

/* compute merit function at (R + alpha * DR) */
double LINSYS_Merit (LINSYS *sys, double alpha)
{
  short variational, dynamic, smooth, variant;
  double H [3], value, step, epsilon;
  DIAB *dia;
  SET *item;
  CON *con;

  value = 0;
  step = sys->ldy->dom->step;
  dynamic = sys->ldy->dom->dynamic;
  variational = sys->variant & (SMOOTHED_VARIATIONAL|NONSMOOTH_VARIATIONAL);
  if (sys->variant & SMOOTHED_VARIATIONAL) smooth = 1; /* (+++); TODO: make consistent */
  variant = LINEARIZATION_VARIANT (sys->variant);
  epsilon = sys->epsilon;

  for (item = SET_First (sys->subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;

    if (con->kind == CONTACT)
    {
      double R [3], U [3],
	    *R0 = dia->R,
	    *U0 = dia->U;

      if (alpha > 0.0)
      {
	double *RE = RE(con),
	       *DR = DR(con),
	       *DU = DU(con);

	ADDMUL (R0, alpha, DR, R);
	ADDMUL (U0, alpha, DU, U);
	ADD (U, RE, U);
      }
      else
      {
	COPY (R0, R);
	COPY (U0, U);
      }

      double rho = dia->rho,
	     gap = con->gap,
	     fri = con->mat.base->friction,
	     res = con->mat.base->restitution,
	    *V = dia->V,
	     d [3],
	     bound,
	     udash,
	     norm;

      if (variational)
      {
        double S [3],
	       F [3],
	       m [3];

	real_F (res, fri, gap, step, dynamic, V, U, F);
	SUB (R, F, S);
	real_m (fri, smooth, S, epsilon, m);
	ADD (F, m, H);
      }
      else
      {
	if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
	else udash = ((MAX(gap, 0)/step) + U[2]);

	d [0] = R[0] - rho * U[0];
	d [1] = R[1] - rho * U[1];
	d [2] = R[2] - rho * udash;

	norm = sqrt (d[0]*d[0]+d[1]*d[1]);

	switch (variant)
	{
	  case NONSMOOTH_HSW:
	    bound = fri * MAX (0, d[2]);
	    break;
	  case NONSMOOTH_HYBRID:
	    bound = fri * MAX (0, R[2]);
	    break;
	  case FIXED_POINT:
	    bound = fri * MAX (0, RN(con));
	    break;
	}

	H [0] = MAX (bound, norm)*R[0] - bound*d[0];
	H [1] = MAX (bound, norm)*R[1] - bound*d[1];
	H [2] = R [2] - MAX (0, d[2]);
      }

      value += DOT (H, H);
    }
  }

#if MPI
  double val_i = value;
  MPI_Allreduce (&val_i, &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  return value;
}

/* advance solution R = R + alpha * DR; return |DR|/|R| */ 
double LINSYS_Advance (LINSYS *sys, double alpha)
{
  double errup, errlo, error;
  SET *item;
  CON *con;

  errup = errlo = 0.0;

  for (item = SET_First (sys->subset); item; item = SET_Next (item))
  {
    con = item->data;

    double *RE = RE(con),
	   *DR = DR(con),
	   *DU = DU(con),
	   *R = con->R,
	   *U = con->U;

    ADDMUL (R, alpha, DR, R);
    ADDMUL (U, alpha, DU, U);
    ADD (U, RE, U);

    errup += DOT(DR, DR); /* no alpha scaling here */
    errlo += DOT(R, R);
  }

#if MPI
  double errloc [2] = {errup, errlo}, errsum [2];
  MPI_Allreduce (errloc, errsum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  errup = errsum [0], errlo = errsum [1];
#endif

  error = sqrt (errup / MAX (errlo, 1.0));

  if ((sys->variant & FIXED_POINT) &&
      error < sys->fixed_point_tol) error = update_normal_bounds (sys);

  return error;
}

/* solve A x = b, where b = A [1, 1, ..., 1]' and return |x - [1, 1, ..., 1]| / |[1, 1, ..., 1]| */
double LINSYS_Test (LINSYS *sys, double abstol, int maxiter)
{
  double *x, *y, error;
  struct vect *vx, *vb;

  vx = CreateVector (sys->x);
  vb = CreateVector (sys->b);

  for (x = vx->x, y = x + vx->n; x < y; x ++) *x = 1.0;

  A_times_x_equals_y (sys, vx->x, vb->x);

  for (x = vx->x, y = x + vx->n; x < y; x ++) *x = 0.0;

  hypre_FlexGMRESSetTol (sys->gmres_vdata, 0.0);
  hypre_FlexGMRESSetMaxIter (sys->gmres_vdata, maxiter);
  hypre_FlexGMRESSetAbsoluteTol (sys->gmres_vdata, abstol);
  hypre_FlexGMRESSetup (sys->gmres_vdata, sys, vb, vx);
  hypre_FlexGMRESSolve (sys->gmres_vdata, sys, vb, vx);
  hypre_FlexGMRESGetNumIterations (sys->gmres_vdata , &sys->iters);
  hypre_FlexGMRESGetFinalRelativeResidualNorm (sys->gmres_vdata, &sys->resnorm);

  for (error = 0.0, x = vx->x, y = x + vx->n; x < y; x ++) 
  {
    error += (1.0 - (*x))*(1.0 - (*x));
  }

#if MPI
  double val_i [2] = {error, vx->n}, val_o [2];
  MPI_Allreduce (val_i, val_o, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  error = sqrt (val_o [0] / MAX (val_o [1], 1.0));
#else
  error = sqrt (error / (double) MAX (vx->n, 1));
#endif

  DestroyVector (vx);
  DestroyVector (vb);

  return error;
}

/* most recent iterations count */
int LINSYS_Iters (LINSYS *sys)
{
  return sys->iters;
}

/* most recent relative residual norm */
double LINSYS_Resnorm (LINSYS *sys)
{
  return sys->resnorm;
}

#if MPI
/* update computed external reactions (e.g. non-gluing) */
void LINSYS_Update_External_Reactions (LINSYS *sys)
{
  double *xx = sys->x->x, *x, *R;
  SET *item;
  CON *con;

  for (item = SET_First (item); item; item = SET_Next (item))
  {
    con = item->data;
    x = &xx [3*con->num];
    R = con->R;
    COPY (R, x);
  }

  update_external_reactions (sys, xx);
}
#endif

/* destroy linear system */
void LINSYS_Destroy (LINSYS *sys)
{
  MEM_Release (&sys->setmem);
  MEM_Release (&sys->Y_mem);
  MEM_Release (&sys->T_mem);
#if USE_LU
  MEM_Release (&sys->mapmem);
  if (sys->A) MX_Destroy (sys->A);
#endif
  DestroyVector (sys->x);
  DestroyVector (sys->b);
  hypre_FlexGMRESDestroy (sys->gmres_vdata);
#if MPI
  COMALL_Free (sys->pattern);
  free (sys->send);
  free (sys->recv);
  free (sys->R);
  free (sys->basenum);
#endif
  free (sys);
}