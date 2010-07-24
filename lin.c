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
#include "sol.h"
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

#define DUMP                    0      /* FIXME: debug dump */
#define DIFF_FACTOR             1E-10  /* TODO: test sensitivity */
#define EPSILON_BASE            1E-10  /* TODO: -||- */
#define EPSILON_FACTOR          1E-6   /* TODO: -||- */
#define ABSTOL_BASE             1E-15  /* TODO: -||- */ 
#define SMOOTHING               1      /* TODO: -||- */
#define DISABLE_NORM_SMOOTHING  1      /* TODO: -||- */
#define BLOCKS         256

/* timers */
#if TIMERS
#define INIT_TIMERS(sys) SOLFEC *solfec = sys->ldy->dom->solfec
#define S(LABEL) SOLFEC_Timer_Start (solfec, LABEL)
#define E(LABEL) SOLFEC_Timer_End (solfec, LABEL)
#else
#define INIT_TIMERS(sys)
#define S(LABEL)
#define E(LABEL)
#endif

typedef struct vect VECT;
typedef struct offt OFFT;
typedef struct offx OFFX;
typedef struct diat DIAT;

struct vect /* krylov solver's vector */
{
  LINSYS *sys;  

  double *x;

  int n;
};

#define vect(ptr) ((VECT*)ptr)

struct offt /* off-diagonal tangent block */
{
  double T [9],
	*W [2], /* pairs of same-body-pair constraints have two Wij entries in LOCDYN */
	*DR,    /* points to DIAT->DR for local constraints or to CON->R for external constraints */
	*R;     /* points to CON->R */

  int num,      /* local number: CON->num */
     *map;      /* maps columns of T in CSC storage */

#if MPI
  int rank; /* rank != dom->rank => external block => use R in matrix-vector product */
#endif

  OFFT *n;
};

struct offx /* off-diagonal inactive LOCDYN block */
{
  double *W,
	 *R;

  OFFX *n;
};

struct diat /* diagonal tangent block */
{
  double  T [9],
	 DR [3],
	 DU [3],
	 RE [3], /* residual */
	 RN,     /* normal stress: fixed point approach */
	 *R,
	 *U,
	 *V,
	 *W,
	  rho,
	  B [3]; /* assembled here, includes WijRj from outside of the constraints subset */

  CON *con;

  int num, /* local number: CON->num */
     *map;

  OFFT *adj;

  OFFX *adx; /* used to update free velocity only */

  DIAT *n;
};

struct linsys
{
  LINVAR variant;
  LOCDYN *ldy;
  DIAT *dia;  /* A: tangent matrix */
  int ndia;   /* blocks count: dim (A) / 3 */
  MX *A;      /* A: compressed */

  MEM diamem, /* DIAT */
      offmem, /* OFFT */
      ofxmem, /* OFFX */
      mapmem; /* {DIAT, OFFT}->map */

  double epsilon,         /* SMOOTHED_VARIATIONAL: smoothing thickness */
         delta;           /* Tikhonov regularaization */

  VECT *x, *b; /* A x = b, where x is composed of reaction increments DR */

  double resnorm, xnorm; /* |Ax - b|, |x| */
  int iters,            /* most recent iterations count */
      smooth;          /* SMOOTHED_VARIATIONAL: smoothing variant */

#if MPI
  COMDATA *send, *recv;
  int nsend, nrecv;
  void *pattern;  /* repetitice communication pattern in matrix-vector products */
  double **R;    /* external con->R */
  int *basenum; /* basenum [i] = SUM {0, i-1} ndia [i-1], where i is a processor number; basenum [ncpu] stores the total sum */
#endif
};

#if !MPI
/* compressed image of system matrix */
static MX* matrix_create (LINSYS *sys)
{
  int i, j, k, l, n, ncon, *rows, *cols, *aux;
  MAP **pcol, *item;
  DIAT *dia;
  OFFT *blk;
  MEM mem;
  MX *A;

  ncon = sys->ndia;
  ERRMEM (A = MEM_CALLOC (sizeof (MX)));
  A->kind = MXCSC;

  /* allocate colum blocks mapping */
  ERRMEM (pcol = MEM_CALLOC (ncon * sizeof (MAP*)));
  MEM_Init (&mem, sizeof (MAP), BLOCKS);

  /* map colum blocks */
  for (dia = sys->dia; dia; dia = dia->n)
  {
    int col, row, *map;

    col = row = dia->num;
    ASSERT_DEBUG (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Diagonal block mapped twice");
    ERRMEM (map = MEM_Alloc (&sys->mapmem));
    ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL)); /* '+1' not to confuse with NULL */

    for (blk = dia->adj; blk; blk = blk->n)
    {
      col = blk->num;
      ASSERT_DEBUG (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Diagonal block mapped twice");
      ERRMEM (map = MEM_Alloc (&sys->mapmem));
      ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL));
    }
  }

  for (j = 0, A->nzmax = 0; j < ncon; j ++)
  {
    A->nzmax += 9 * MAP_Size (pcol [j]); /* number of nonzeros */
  }
  A->n = A->m = ncon * 3; /* system dimension */
  A->nz = -1; /* compressed columns */

  /* allocate compressed column storage */
  ERRMEM (A->x = MEM_CALLOC (sizeof (double) * A->nzmax));
  ERRMEM (A->i = malloc (sizeof (int) * A->nzmax));
  ERRMEM (A->p = MEM_CALLOC (sizeof (int) * (A->n+1))); /* '+ 1' as there is cols[0] == 0 always */
  ERRMEM (aux = malloc (sizeof (int) * A->n));
  rows = A->i;
  cols = A->p;

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
  for (dia = sys->dia; dia; dia = dia->n)
  {
    int col, row;

    col = row = dia->num;
    ASSERT_DEBUG_EXT (dia->map = MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Inconsistent sparse storage");

    for (blk = dia->adj; blk; blk = blk->n)
    {
      col = blk->num;
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
  double *a = A->x;
  DIAT *dia;
  OFFT *blk;

  for (double *b = a, *c = a + A->nzmax; b < c; b ++) *b = 0.0;

  for (dia = sys->dia; dia; dia = dia->n)
  {
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

/* update linearization of a non-contact constraint */
static void system_update_noncontact (DIAT *dia, short dynamic, double step, double *b)
{
  OFFT *blk;

  CON *con = dia->con;

  double *W = dia->W,
	 *U = dia->U,
	 *R = dia->R,
	 *T = dia->T,
	 *RE = dia->RE; /* assumed updated */

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
      double *W = blk->W [0], *T = blk->T;

      NNCOPY (W, T);
    }
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
      double **W = blk->W, *T = blk->T;

      T [2] = W[0][2];
      T [5] = W[0][5];
      T [8] = W[0][8];

      if (W[1])
      {
	T [2] += W[1][2];
	T [5] += W[1][5];
	T [8] += W[1][8];
      }

      T [0] = T [1] = T [3] = T [4] = T [6] = T [7] = 0;
    }
  }
  break;
  }
}

/* update linear system for NONSMOOTH_HSW, NONSMOOTH_HYBRID, FIXED_POINT variants */
static void system_update_HSW_HYBRID_FIXED (LINSYS *sys, double *rhs)
{
  double d [3], norm, lim, udash, step, *b;
  short dynamic, pull;
  LINVAR variant;
  DIAT *dia;
  OFFT *blk;
  DOM *dom;
  CON *con;

  variant = sys->variant;
  dom = sys->ldy->dom;
  dynamic = dom->dynamic;
  step = dom->step;

  for (dia = sys->dia; dia; dia = dia->n)
  {
    con = dia->con;
    b = &rhs [3*dia->num];

    switch (con->kind)
    {
    case CONTACT:
    {
      double *V = dia->V,
	     *W = dia->W,
	     *U = dia->U,
	     *R = dia->R,
	     *T = dia->T,
	     *RE = dia->RE,
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
	  double **W = blk->W, *T = blk->T;

	  T [2] = W[0][2];
	  T [5] = W[0][5];
	  T [8] = W[0][8];

	  if (W[1])
	  {
	    T [2] += W[1][2];
	    T [5] += W[1][5];
	    T [8] += W[1][8];
	  }
	}
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
	  double *T = blk->T;

	  T [2] = T [5] = T [8] = 0;
	}
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
	  lim = fri * MAX (0, dia->RN);
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

	    b [0] = (fri*(d[0]/norm)*dia->RN - R[0]) * W[0] - rho*(M[0]*RE[0] + M[2]*RE[1]); 
	    b [1] = (fri*(d[1]/norm)*dia->RN - R[1]) * W[4] - rho*(M[1]*RE[0] + M[3]*RE[1]); 
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
	    double **W = blk->W, *T = blk->T;

	    T[0] = rho*(M[0]*W[0][0] + M[2]*W[0][1]);
	    T[1] = rho*(M[1]*W[0][0] + M[3]*W[0][1]);
	    T[3] = rho*(M[0]*W[0][3] + M[2]*W[0][4]);
	    T[4] = rho*(M[1]*W[0][3] + M[3]*W[0][4]);
	    T[6] = rho*(M[0]*W[0][6] + M[2]*W[0][7]); 
	    T[7] = rho*(M[1]*W[0][6] + M[3]*W[0][7]);  

	    if (W[1])
	    {
	      T[0] += rho*(M[0]*W[1][0] + M[2]*W[1][1]);
	      T[1] += rho*(M[1]*W[1][0] + M[3]*W[1][1]);
	      T[3] += rho*(M[0]*W[1][3] + M[2]*W[1][4]);
	      T[4] += rho*(M[1]*W[1][3] + M[3]*W[1][4]);
	      T[6] += rho*(M[0]*W[1][6] + M[2]*W[1][7]);
	      T[7] += rho*(M[1]*W[1][6] + M[3]*W[1][7]);
	    }
	  }
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
	    double *T = blk->T;

	    T [0] = T [1] = T [3] = T [4] = T [6] = T [7] = 0;
	  }
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
	  double **W = blk->W, *T = blk->T;

	  T [0] = W[0][0];
	  T [1] = W[0][1];
	  T [3] = W[0][3];
	  T [4] = W[0][4];
	  T [6] = W[0][6];
	  T [7] = W[0][7];

	  if (W[1])
	  {
	    T [0] += W[1][0];
	    T [1] += W[1][1];
	    T [3] += W[1][3];
	    T [4] += W[1][4];
	    T [6] += W[1][6];
	    T [7] += W[1][7];
	  }
	}
      }
    }
    break;
    default:
    {
      system_update_noncontact (dia, dynamic, step, b);
    }
    break;
    }
  }
}

/* imaginary i */
static double complex imaginary_i;

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
inline static void real_F (double res, double fri, double gap, double step, short dynamic, double epsilon, double *V, double *U, double *F)
{
  double udash;

  if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
  else udash = ((MAX(gap, 0)/step) + U[2]);

#if DISABLE_NORM_SMOOTHING
  epsilon = 0;
#endif

  F [0] = U[0];
  F [1] = U[1];
  F [2] = (udash + fri * sqrt (DOT2(U, U) + epsilon*epsilon));
}
 
/* complex F = [UT, UN + fri |UT|]' */
inline static void complex_F (double res, double fri, double gap, double step, short dynamic, double epsilon, double *V, double complex *U, double complex *F)
{
  double complex udash;

  if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
  else udash = ((MAX(gap, 0)/step) + U[2]);

#if DISABLE_NORM_SMOOTHING
  epsilon = 0;
#endif

  F [0] = U[0];
  F [1] = U[1];
  F [2] = (udash + fri * csqrt (DOT2(U, U) + epsilon*epsilon));
}

/* update linear system for NONSMOOTH_VARIATIONAL, SMOOTHED_VARIATIONAL variants */
static void system_update_VARIATIONAL (LINSYS *sys, double *rhs)
{
  double epsilon, step, h, *b;
  short smooth, dynamic;
  OFFT *blk;
  DIAT *dia;
  DOM *dom;
  CON *con;

  smooth = sys->smooth;
  epsilon = sys->epsilon;
  h = DIFF_FACTOR * (epsilon == 0.0 ? 1.0 : epsilon);
  dom = sys->ldy->dom;
  dynamic = dom->dynamic;
  step = dom->step;

  for (dia = sys->dia; dia; dia = dia->n)
  {
    con = dia->con;
    b = &rhs [3*dia->num];

    switch (con->kind)
    {
    case CONTACT:
    {
      double *V = dia->V,
	     *W = dia->W,
	     *U = dia->U,
	     *R = dia->R,
	     *T  = dia->T,
	     *RE = dia->RE,
	     gap = con->gap,
	     fri = con->mat.base->friction,
	     res = con->mat.base->restitution,
	     TT [9],
	     dm [9],
	     dF [9],
             S [3],
	     F [3],
	     m [3],
	     H [3],
	     J [9];


      real_F (res, fri, gap, step, dynamic, epsilon, V, U, F);
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
        complex_F (res, fri, gap, step, dynamic, epsilon, V, cU, cF);
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
	double **W = blk->W, *T = blk->T;

        NNMUL (TT, W[0], T);

	if (W[1])
	{
	  NNMUL (TT, W[1], J);
	  NNADD (T, J, T);
	}
      }
    }
    break;
    default:
    {
      system_update_noncontact (dia, dynamic, step, b);
    }
    break;
    }
  }
}

#if MPI
static void update_external_reactions (LINSYS *sys, double *x)
{
  double **r, *z;
  int i, *j, *k;
  COMDATA *d;
  SET *jtem;
  DIAT *dia;
  CON *con;

  for (dia = sys->dia, d = sys->send; dia; dia = dia->n)
  {
    con = dia->con;
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
  double *T, *z, *v;
  INIT_TIMERS (sys);
  OFFT *blk;
  DIAT *dia;

#if MPI
  S("LINCOM");
  update_external_reactions (sys, x); /* (###) */
  E("LINCOM");

  int rank = sys->ldy->dom->rank;
#endif

  S("LINMV");
  for (dia = sys->dia; dia; dia = dia->n)
  {
    v = &y [3*dia->num];
    z = &x [3*dia->num];
    T = dia->T;

    NVMUL (T, z, v);

    for (blk = dia->adj; blk; blk = blk->n)
    {
#if MPI
      if (blk->rank != rank) z = blk->R; /* (###) */
      else
#endif
      z = &x [3*blk->num];
      T = blk->T;

      NVADDMUL (v, T, z, v);
    }
  }
  E("LINMV");
}

static void AT_times_x_equals_y (LINSYS *sys, double *x, double *y)
{
#if MPI
  int *basenum, ncpu, rank, gdim;
  double *v;
#endif
  double *T, *z, *w;
  INIT_TIMERS (sys);
  OFFT *blk;
  DIAT *dia;
  DOM *dom;

  S("LINMV");
  dom = sys->ldy->dom;
#if MPI
  ncpu = dom->ncpu;
  rank = dom->rank;
  basenum = sys->basenum;
  gdim = 3 * basenum [ncpu];
  ERRMEM (v = MEM_CALLOC (2 * gdim * sizeof (double))); /* global result vector */
#else
  for (z = y, w = z + sys->b->n; z < w; z ++) (*z) = 0.0; /* zero y */
#endif

  for (dia = sys->dia; dia; dia = dia->n)
  {
#if MPI
    w = &v [3*(basenum [rank]+dia->num)];
#else
    w = &y [3*dia->num];
#endif
    z = &x [3*dia->num];
    T = dia->T;

    TVADDMUL (w, T, z, w);

    for (blk = dia->adj; blk; blk = blk->n)
    {
#if MPI
      w = &v [3*(basenum [blk->rank]+blk->num)];
#else
      w = &y [3*blk->num];
#endif
      T = blk->T;

      TVADDMUL (w, T, z, w);
    }
  }
  E("LINMV");

#if MPI
  S("LINCOM");
  blas_dcopy (gdim, v, 1, v + gdim, 1);
  MPI_Allreduce (v + gdim, v, gdim, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  for (w = &v [3*basenum [rank]], z = &v [3*basenum [rank+1]]; w < z; w ++, y ++) (*y) = (*w);
  free (v);
  E("LINCOM");
#endif
}

static VECT* newvect (int n, LINSYS *sys)
{
  VECT *v;

  ERRMEM (v = malloc (sizeof (VECT)));
  ERRMEM (v->x = MEM_CALLOC (n * sizeof (double)));
  v->sys = sys;
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

  *num_procs = dom->ncpu;
  *my_id = dom->rank;
#else
  *num_procs = 1;
  *my_id = 0;
#endif
  return 0;
}

static void* CreateVector (void *vector)
{
  VECT *a = vect (vector), *v;

  ERRMEM (v = malloc (sizeof (VECT)));
  ERRMEM (v->x = MEM_CALLOC (a->n * sizeof (double)));
  v->sys = a->sys;
  v->n = a->n;

  return v;
}

static void* CreateVectorArray (int size, void *vectors)
{
  VECT **v;
  int i;

  ERRMEM (v = malloc (size * sizeof (VECT*)));
  for (i = 0; i < size; i ++)
  {
    v[i] = CreateVector (vectors);
  }

  return v;
}

static int DestroyVector (void *vector)
{
  VECT *v = vect (vector);
  free (v->x);
  free (v);

  return 0;
}

static double InnerProd (void *vx, void *vy)
{
  INIT_TIMERS (vect (vx)->sys);
  double dot = 0.0, *x, *y, *z;

  S("LINRUN");
  for (x = vect (vx)->x, z = x + vect (vx)->n, y = vect (vy)->x; x < z; x ++, y ++)
  {
    dot += (*x) * (*y);
  }
  E("LINRUN");

#if MPI
  S("LINCOM");
  double val = dot;
  MPI_Allreduce (&val, &dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  E("LINCOM");
#endif

  return dot;
}

static int CopyVector (void *vx, void *vy)
{
  INIT_TIMERS (vect (vx)->sys);
  double *x, *y, *z;

  S("LINRUN");
  for (x = vect (vx)->x, z = x + vect (vx)->n, y = vect (vy)->x; x < z; x ++, y ++)
  {
    (*y) = (*x);
  }
  E("LINRUN");

  return 0;
}

static int ClearVector (void *vx)
{
  INIT_TIMERS (vect (vx)->sys);
  double *x, *z;

  S("LINRUN");
  for (x = vect (vx)->x, z = x + vect (vx)->n; x < z; x ++)
  {
    (*x) = 0.0;
  }
  E("LINRUN");

  return 0;
}

static int ScaleVector (double alpha, void *vx)
{
  INIT_TIMERS (vect (vx)->sys);
  double *x, *z;

  S("LINRUN");
  for (x = vect (vx)->x, z = x + vect (vx)->n; x < z; x ++)
  {
    (*x) *= alpha;
  }
  E("LINRUN");

  return 0;
}

static int  Axpy (double alpha, void *vx, void *vy )
{
  INIT_TIMERS (vect (vx)->sys);
  double *x, *y, *z;

  S("LINRUN");
  for (x = vect (vx)->x, z = x + vect (vx)->n, y = vect (vy)->x; x < z; x ++, y ++)
  {
    (*y) += alpha * (*x);
  }
  E("LINRUN");

  return 0;
}

static void *MatvecCreate (void *A, void *x)
{
  return NULL;
}

static int Matvec (void *matvec_data, double alpha, void *A, void *vx, double beta, void *vy)
{
  VECT *vz = CreateVector (vx);
  LINSYS *sys = (LINSYS*)A;
  double delta = sys->delta;

  ScaleVector (beta, vy);
  if (delta > 0.0) Axpy (alpha*delta, vx, vy);
  A_times_x_equals_y (sys, vect (vx)->x, vz->x);

  if (sys->variant & MULTIPLY_TRANSPOSED)
  {
    VECT *vu = CreateVector (vx);
    AT_times_x_equals_y (sys, vz->x, vu->x);
    Axpy (alpha, vu, vy);
    DestroyVector (vu);
  }
  else Axpy (alpha, vz, vy);

  DestroyVector (vz);

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
  double *bb = vect (vb)->x, *xx = vect (vx)->x, *b, *x;
  LINSYS *sys = (LINSYS*)A;
  double delta = sys->delta;
  short trans = sys->variant & MULTIPLY_TRANSPOSED;
  INIT_TIMERS (sys);
  int ipiv [3];
  DIAT *dia;

  S("LINPRE");
  for (dia = sys->dia; dia; dia = dia->n)
  {
    b = &bb [3*dia->num];
    x = &xx [3*dia->num];
    COPY (b, x);

    double *T = dia->T, S [9];

    if (trans) { TNMUL (T, T, S); }
    else { NNCOPY (T, S); }
    S [0] += delta;
    S [4] += delta;
    S [8] += delta;
    lapack_dgesv (3, 1, S, 3, ipiv, x, 3); /* TODO: inv (S) could be precomputed at the start */
  }
  E("LINPRE");

  return 0;
}
/* GMRES interface end */

/* compute residual and solution norms */
static void compute_resnorm_and_xnorm (LINSYS *sys)
{
  VECT *r = CreateVector (sys->b);

  A_times_x_equals_y (sys, sys->x->x, r->x);

  if (sys->variant & MULTIPLY_TRANSPOSED)
  {
    VECT *s = CreateVector (sys->b);
    CopyVector (r, s);
    AT_times_x_equals_y (sys, s->x, r->x);
    DestroyVector (s);
  }

  Axpy (-1.0, sys->b, r);

  sys->xnorm = sqrt (InnerProd (sys->x, sys->x));
  sys->resnorm = sqrt (InnerProd (r, r));

  DestroyVector (r);
}

/* project contact reactions back onto the friction cone */
static void project_contact_reactions (LINSYS *sys)
{
  double S [3], m [3], fri, *R, *a, *x = sys->x->x;
  DIAT *dia;
  CON *con;

  if (sys->variant & MULTIPLY_TRANSPOSED)
  {
    VECT *v = CreateVector (sys->x);
    double beta, dot0, dot1, *b, *y = v->x;

    dot1 = 1.0;
    beta = 1.0;
    do
    {
      for (dia = sys->dia; dia; dia = dia->n)
      {
	con = dia->con;
	if (con->kind == CONTACT)
	{
	  a = &x [3*dia->num];
	  b = &y [3*dia->num];
	  MUL (a, beta, b);
	  fri = con->mat.base->friction;
	  R = con->R;
	  ADD (R, b, S);
	  real_m (fri, 0, S, 0, m);
	  SUB (S, m, S);
	  SUB (S, R, b); /* v = proj (friction-cone, R+beta*DR) - R */
	}
	/* else it is zero */
      }

      dot0 = dot1;
      dot1 = -InnerProd (v, sys->b); /* b = - grad (f), hence <v, -b> = <v, grad (f) < 0 indicates descent */

      if (dot0 == 1.0 && dot1 < 0.0) break; /* skip if first projection is negative */
#if 0
      /* FIXME: remove */
      printf ("LINSYS: dot = %g, beta = %g\n", dot1, beta);
#endif
    }
    while (dot1 < dot0 && (beta *= 0.9) > 0.05);

    if (beta < 1.0) beta /= 0.9; /* corresponds to the minimum of the inner product */

    for (dia = sys->dia; dia; dia = dia->n)
    {
      con = dia->con;
      if (con->kind == CONTACT)
      {
        a = &x [3*dia->num];
        b = &y [3*dia->num];
	COPY (b, a);
      }
    }

    DestroyVector (v);
  }
  else
  {
    for (dia = sys->dia; dia; dia = dia->n)
    {
      con = dia->con;
      if (con->kind == CONTACT) /* DR = proj (friction-cone, R+DR) - R */
      {
	fri = con->mat.base->friction;
	R = con->R;
        a = &x [3*dia->num];
	ADD (R, a, S);
	real_m (fri, 0, S, 0, m);
	SUB (S, m, S);
	SUB (S, R, a);
      }
    }
  }
}

#if DUMP
static int dump_cmp (const void *a, const void *b)
{
  double *x = (double*)a, *y = (double*)b;
  int i;

  for (i = 0; i < 9; i ++)
  {
    if (x [i] < y [i]) return -1;
    else if (x [i] > y [i]) return 1;
  }

  return 0;
}

static void linsys_dump (LINSYS *sys, const char *path)
{
  char fullpath [512];
  double *x, *b;
  DIAT *dia;
  OFFT *blk;
  CON *con;
  FILE *f;
  int i, n;

#if MPI
  int rank;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  snprintf (fullpath, 512, "%s.%d", path, rank);
#else
  snprintf (fullpath, 512, "%s", path);
#endif
  ASSERT (f = fopen (fullpath, "w"), ERR_FILE_OPEN);

  for (dia = sys->dia; dia; dia = dia->n)
  {
    con = dia->con;

    b = &sys->b->x [3*dia->num];

    fprintf (f, "(%.6g, %.6g, %.6g) (%.6g, %.6g, %.6g)", con->point [0], con->point [1], con->point [2], b [0], b [1], b [2]);

    for (n = 1, blk = dia->adj; blk; blk = blk->n, n ++);

    ERRMEM (x = malloc (n * sizeof (double [9])));

    NNCOPY (dia->T, x);
    for (n = 9, blk = dia->adj; blk; blk = blk->n, n += 9)
    { 
      NNCOPY (blk->T, &x[n]);
    }

    qsort (x, n/9, sizeof (double [9]), dump_cmp);

    for (i = 0; i < n; i += 9)
    {
      fprintf (f, " (%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g )",
        x [i], x [i+1], x [i+2], x [i+3], x [i+4], x [i+5], x [i+6], x [i+7], x [i+8]);
    }

    free (x);

    fprintf (f, "\n");
  }

  fclose (f);
}
#endif

/* create linear system resulting from linearization of constraints */
LINSYS* LINSYS_Create (LINVAR variant, LOCDYN *ldy, SET *subset)
{
  DOM *dom = ldy->dom;
  MEM setmem, mapmem;
  int ncon, blocks;
  short dynamic;
  LINSYS *sys;
  SET *item;
  DIAB *dia;
  OFFB *blk;
  DIAT *dt;
  OFFT *bt;
  OFFX *bx;
  CON *con;
  MAP *map;

  dynamic = dom->dynamic;
  ERRMEM (sys = MEM_CALLOC (sizeof (LINSYS)));
  sys->variant = variant;
  sys->ldy = ldy;
  INIT_TIMERS (sys);

  S("LININIT");
  MEM_Init (&setmem, sizeof (SET), BLOCKS);
  MEM_Init (&mapmem, sizeof (MAP), BLOCKS);

  imaginary_i = csqrt (-1);

  if (!subset)
  {
    for (con = dom->con; con; con = con->next)
    {
      if (dynamic && con->kind == CONTACT && con->gap > 0) continue; /* skip open dynamic contact */

      SET_Insert (&setmem, &subset, con, NULL); /* insert all but the open contacts */
    }
  }

  /* number constraints */
  for (ncon = 0, item = SET_First (subset); item; item = SET_Next (item))
  {
    con = item->data;
    con->num = ncon ++;
  }
  sys->ndia = ncon;

  blocks = MAX (ncon, BLOCKS);
  MEM_Init (&sys->mapmem, sizeof (int [3]), blocks);
  MEM_Init (&sys->diamem, sizeof (DIAT), blocks);
  MEM_Init (&sys->offmem, sizeof (OFFT), blocks);
  MEM_Init (&sys->ofxmem, sizeof (OFFX), blocks);

#if MPI
  /* propagate information about the subset and numbering */
  COMDATA *send, *recv, *d;
  int nrecv, i, *j, *k;
  SET *subext, *jtem;

  ERRMEM (send = MEM_CALLOC (dom->ncpu * sizeof (COMDATA)));

  for (i = 0; i < dom->ncpu; i ++) send [i].rank = i;

  for (item = SET_First (subset); item; item = SET_Next (item))
  {
    con = item->data;

    for (jtem = SET_First (con->ext); jtem; jtem = SET_Next (jtem))
    {
      d = &send [(int) (long) jtem->data];

      d->ints += 2;
      ERRMEM (d->i = realloc (d->i, d->ints * sizeof (int)));
      d->i [d->ints-2] = con->id;
      d->i [d->ints-1] = con->num;
    }
  }

  COMALL (MPI_COMM_WORLD, send, dom->ncpu, &recv, &nrecv);

  for (i = 0, subext = NULL; i < nrecv; i ++)
  {
    d = &recv [i];
    for (j = d->i, k = j + d->ints; j < k; j += 2)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) j[0], NULL), "Invalid constraint id");
      con->num = j[1];
      SET_Insert (&setmem, &subext, con, NULL);
    }
  }

  for (i = 0; i < dom->ncpu; i ++) free (send [i].i);
  free (send);
  free (recv);
#endif

  /* create linearized system structure */
  for (map = NULL, item = SET_First (subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;

    ERRMEM (dt = MEM_Alloc (&sys->diamem));
    dt->R = dia->R;
    dt->U = dia->U;
    dt->V = dia->V;
    dt->W = dia->W;
    dt->rho = dia->rho;
    dt->con = con;
    dt->num = con->num;

    double *B = dt->B;
    COPY (dia->B, B);

    MAP_Free (&mapmem, &map);

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
      if (SET_Contains (subset, con, NULL))
      {
	if ((bt = MAP_Find (map, con, NULL))) /* second Wij adjacent to the same DIAB */
	{
	  ASSERT_DEBUG (bt->W [1] == NULL, "Inconsistent mappint of second Wij");
	  bt->W [1] = blk->W;
	}
	else
	{
	  ERRMEM (bt = MEM_Alloc (&sys->offmem));
          bt->W [0] = blk->W;
	  bt->num = con->num;
	  bt->R = con->R;
#if MPI
	  bt->rank = dom->rank;
#endif
	  MAP_Insert (&mapmem, &map, con, bt, NULL); /* map first Wij */

	  bt->n = dt->adj;
	  dt->adj = bt;
	}
      }
      else
      {
	double *R = blk->dia->R, *W = blk->W;

	NVADDMUL (B, W, R, B);

	ERRMEM (bx = MEM_Alloc (&sys->ofxmem));
	bx->W = W;
	bx->R = R;
	bx->n = dt->adx;
	dt->adx = bx;
      }
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = CON (blk->dia);
      if (SET_Contains (subext, con, NULL))
      {
	if ((bt = MAP_Find (map, con, NULL))) /* second Wij adjacent to the same DIAB */
	{
	  ASSERT_DEBUG (bt->W[1] == NULL, "Inconsistent mappint of second Wij");
	  bt->W [1] = blk->W;
	}
	else
	{
	  ERRMEM (bt = MEM_Alloc (&sys->offmem));
          bt->W [0] = blk->W;
	  bt->num = con->num;
	  bt->R = con->R;
	  bt->DR = con->R; /* same */
	  bt->rank = con->rank;

	  MAP_Insert (&mapmem, &map, blk->dia, bt, NULL); /* map first Wij */

	  bt->n = dt->adj;
	  dt->adj = bt;
	}
      }
      else
      {
	double *R = CON(blk->dia)->R, *W = blk->W;

	NVADDMUL (B, W, R, B);

	ERRMEM (bx = MEM_Alloc (&sys->ofxmem));
	bx->W = W;
	bx->R = R;
	bx->n = dt->adx;
	dt->adx = bx;
      }
    }
#endif

    dt->n = sys->dia;
    sys->dia = dt;
  }

  /* map off-diagonal local DR */
  MAP_Free (&mapmem, &map);
  for (dt = sys->dia; dt; dt = dt->n)
  {
    MAP_Insert (&mapmem, &map, (void*) (long) dt->num, dt->DR, NULL);
  }
  for (dt = sys->dia; dt; dt = dt->n)
  {
    for (bt = dt->adj; bt; bt = bt->n)
    {
      if (bt->DR == NULL)
      {
	ASSERT_DEBUG_EXT (bt->DR = MAP_Find (map, (void*) (long) bt->num, NULL), "Inconsistent OFFT numbering");
      }
    }
  }

  /* smoothing thickness */
  if (variant & SMOOTHED_VARIATIONAL)
  {
    double B [4] = {0, 0, 0, 0}, len;

    for (dt = sys->dia; dt; dt = dt->n)
    {
      con = dt->con;
      if (con->kind == CONTACT)
      {
        dia = con->dia;
	ACCABS (dia->B, B); /* use original free velocity */
	B [3] += 1.0;
      }
      else if (con->kind == VELODIR)
      {
	B [2] += fabs (VELODIR(con->Z));
	B [3] += 1.0;
      }
    }

#if MPI
    double C [4];
    MPI_Allreduce (B, C, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    COPY4 (C, B);
#endif

    B [3] = 1.0 / MAX (B [3], 1.0);
    SCALE (B, B[3]); /* avergae free contact velocity */
    len = LEN (B);

    sys->epsilon = EPSILON_FACTOR * (len == 0.0 ? EPSILON_BASE : len);
    sys->smooth = SMOOTHING;
  }
  else { sys->epsilon = 0; sys->smooth = 0; }

  /* unknown and right hand side */
  sys->x = newvect (3 * ncon, sys);
  sys->b = CreateVector (sys->x);
  sys->resnorm = 0.0; /* => delta = resnorm / xnorm == 0 at first */
  sys->xnorm = 1.0;

#if MPI
  /* communication pattern */
  double **r, *x = sys->x->x;
  int n;

  for (dt = sys->dia; dt; dt = dt->n)
  {
    con = dt->con;
    for (jtem = SET_First (con->ext); jtem; jtem = SET_Next (jtem))
    {
      sys->nsend ++;
    }
  }

  ERRMEM (sys->send = MEM_CALLOC (sys->nsend * sizeof (COMDATA)));
  d = sys->send;

  for (dt = sys->dia; dt; dt = dt->n)
  {
    con = dt->con;
    for (jtem = SET_First (con->ext); jtem; jtem = SET_Next (jtem), d ++)
    {
      d->rank = (int) (long) jtem->data;
      d->doubles = 3;
      d->ints = 1;
      d->d = &x [3*con->num];
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

  /* global constraint numbering bases */
  ERRMEM (sys->basenum = MEM_CALLOC ((dom->ncpu + 1) * sizeof (int)));
  MPI_Allgather  (&ncon, 1, MPI_INT, (sys->basenum + 1), 1, MPI_INT, MPI_COMM_WORLD);
  for (i = 1; i <= dom->ncpu; i ++) sys->basenum [i] += sys->basenum [i-1];
#endif

#if !MPI
  if (variant & DIRECT_SOLVE) sys->A = matrix_create (sys);
#endif

  MEM_Release (&setmem);
  MEM_Release (&mapmem);
  E("LININIT");

  return sys;
}

/* update normal reactions for the FIXED_POINT variant */
void LINSYS_Fixed_Point_Update (LINSYS *sys)
{
  INIT_TIMERS (sys);
  DIAT *dia;
  CON *con;

  S("LINRUN");
  for (dia = sys->dia; dia; dia = dia->n)
  {
    con = dia->con;
    if (con->kind == CONTACT)
    {
      dia->RN = con->R[2];
    }
  }
  E("LINRUN");
}

/* update free velocity in case of subset based system */
void LINSYS_Update_Free_Velocity (LINSYS *sys)
{
  INIT_TIMERS (sys);
  DIAT *dia;
  OFFX *blx;

  S("LINRUN");
  for (dia = sys->dia; dia; dia = dia->n)
  {
    double *B0 = dia->con->dia->B,
	   *B = dia->B;

    COPY (B0, B);

    for (blx = dia->adx; blx; blx = blx->n)
    {
      double *W = blx->W, *R = blx->R;

      NVADDMUL (B, W, R, B);
    }
  }
  E("LINRUN");
}

/* update linear system at current reactions R */
void LINSYS_Update (LINSYS *sys, short first)
{
  INIT_TIMERS (sys);
  DIAT *dia;
  OFFT *blk;

  if (first) /* update residual only for the first linear subproblem */
  {
#if MPI
    LINSYS_Update_External_Reactions (sys); /* (###) */ 
#endif

    S("LINUPD");
    for (dia = sys->dia; dia; dia = dia->n)
    {
      double *W = dia->W,
	     *B = dia->B,
	     *U = dia->U,
	     *R = dia->R,
	     *RE = dia->RE;

      NVADDMUL (B, W, R, RE);
      SUB (RE, U, RE);
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double **W = blk->W, *R = blk->R; /* (###) */

	NVADDMUL (RE, W[0], R, RE);
	if (W[1]) { NVADDMUL (RE, W[1], R, RE); }
      }
    }
    E("LINUPD");
  }
  else /* otherwise set to zero */
  {
    S("LINUPD");
    for (dia = sys->dia; dia; dia = dia->n)
    {
      double *RE = dia->RE;
      SET (RE, 0);
    }
    E("LINUPD");
  }

  S("LINUPD");
  switch (LINEARIZATION_VARIANT (sys->variant))
  {
    case NONSMOOTH_HSW:
    case NONSMOOTH_HYBRID:
    case FIXED_POINT:
      system_update_HSW_HYBRID_FIXED (sys, sys->b->x);
      break;
    case NONSMOOTH_VARIATIONAL:
    case SMOOTHED_VARIATIONAL:
      system_update_VARIATIONAL (sys, sys->b->x);
      break;
  }
  E("LINUPD");

  if (sys->variant & MULTIPLY_TRANSPOSED) /* b = (A^T) b */
  {
    S("LINMV");
    VECT *c = CreateVector (sys->b);
    CopyVector (sys->b, c);
    E("LINMV");
    AT_times_x_equals_y (sys, c->x, sys->b->x);
    DestroyVector (c);
  }

#if DUMP
  linsys_dump (sys, "LINSYS");
#endif
}

/* solve for reaction increments DR */
void LINSYS_Solve (LINSYS *sys, double beta, int maxiter)
{
  double *DR, *DU, *W, *x, *z, abstol;
  INIT_TIMERS (sys);
  DIAT *dia;
  OFFT *blk;
  CON *con;

  sys->delta = sys->resnorm / sys->xnorm; /* sqrt (L-curve) */

  abstol = beta * sys->resnorm; /* initially zero => maxiter is reached */

  if (abstol == 0)
  {
    abstol = ABSTOL_BASE * sqrt (InnerProd (sys->b, sys->b));
    if (abstol == 0) abstol = ABSTOL_BASE;
  }

#if !MPI
  if (sys->variant & DIRECT_SOLVE)
  {
    S("LINRUN");
    double *xx = sys->x->x, *bb = sys->b->x;

    for (x = xx + sys->x->n; xx < x; xx ++, bb ++) *xx = *bb;

    matrix_copy (sys, sys->A);

    cs_lusol (1, sys->A, sys->x->x, 0);
    E("LINRUN");
  }
  else
#endif
  if (sys->variant & MULTIPLY_TRANSPOSED)
  {
   hypre_PCGFunctions *pcg_functions;
    void *pcg_vdata;

    pcg_functions = hypre_PCGFunctionsCreate (CAlloc, Free, CommInfo, CreateVector, DestroyVector, MatvecCreate,
      Matvec, MatvecDestroy, InnerProd, CopyVector, ClearVector, ScaleVector, Axpy, PrecondSetup, Precond);
    pcg_vdata = hypre_PCGCreate (pcg_functions);

    hypre_PCGSetTol (pcg_vdata, 0.0);
    hypre_PCGSetMaxIter (pcg_vdata, maxiter);
    hypre_PCGSetAbsoluteTol (pcg_vdata, abstol);
    hypre_PCGSetup (pcg_vdata, sys, sys->b, sys->x);
    hypre_PCGSolve (pcg_vdata, sys, sys->b, sys->x);
    hypre_PCGGetNumIterations (pcg_vdata , &sys->iters);
    hypre_PCGDestroy (pcg_vdata);

    compute_resnorm_and_xnorm (sys);
  }
  else
  {
    hypre_FlexGMRESFunctions *gmres_functions;
    void *gmres_vdata;

    gmres_functions = hypre_FlexGMRESFunctionsCreate (CAlloc, Free, CommInfo, CreateVector, CreateVectorArray, DestroyVector,
      MatvecCreate, Matvec, MatvecDestroy, InnerProd, CopyVector, ClearVector, ScaleVector, Axpy, PrecondSetup, Precond);
    gmres_vdata = hypre_FlexGMRESCreate (gmres_functions);

    hypre_FlexGMRESSetTol (gmres_vdata, 0.0);
    hypre_FlexGMRESSetMinIter (gmres_vdata, 1);
    hypre_FlexGMRESSetMaxIter (gmres_vdata, maxiter);
    hypre_FlexGMRESSetAbsoluteTol (gmres_vdata, abstol);
    hypre_FlexGMRESSetup (gmres_vdata, sys, sys->b, sys->x);
    hypre_FlexGMRESSolve (gmres_vdata, sys, sys->b, sys->x);
    hypre_FlexGMRESGetNumIterations (gmres_vdata , &sys->iters);
    hypre_FlexGMRESDestroy (gmres_vdata);

    compute_resnorm_and_xnorm (sys);
  }

  S("LINRUN");
  /* project for variational formulations, but only in case of iterative solution */
  if ((sys->variant & (SMOOTHED_VARIATIONAL|NONSMOOTH_VARIATIONAL)) && (sys->variant & DIRECT_SOLVE) == 0)
  {
    project_contact_reactions (sys);
  }

  for (dia = sys->dia, x = sys->x->x; dia; dia = dia->n)
  {
    con = dia->con;
    z = &x [3*dia->num];
    DR = dia->DR;
    COPY (z, DR);
  }
  E("LINRUN");

#if MPI
  S("LINCOM");
  update_external_reactions (sys, sys->x->x); /* (&&&) */
  E("LINCOM");
#endif

  /* DU  */
  S("LINRUN");
  for (dia = sys->dia; dia; dia = dia->n)
  {
    DR = dia->DR,
    DU = dia->DU,
    W = dia->W;

    NVMUL (W, DR, DU);

    for (blk = dia->adj; blk; blk = blk->n)
    {
      double **W = blk->W;
      DR = blk->DR; /* (&&&) */

      NVADDMUL (DU, W[0], DR, DU);
      if (W[1]) { NVADDMUL (DU, W[1], DR, DU); }
    }
  }
  E("LINRUN");
}

/* compute merit function at (R + alpha * DR) */
double LINSYS_Merit (LINSYS *sys, double alpha)
{
  short variational, dynamic, smooth, variant;
  double H [3], value, step, epsilon;
  INIT_TIMERS (sys);
  DIAT *dia;
  CON *con;

  S("LINRUN");
  value = 0;
  step = sys->ldy->dom->step;
  dynamic = sys->ldy->dom->dynamic;
  variational = sys->variant & (SMOOTHED_VARIATIONAL|NONSMOOTH_VARIATIONAL);
  variant = LINEARIZATION_VARIANT (sys->variant);
  epsilon = sys->epsilon;
  smooth = sys->smooth;

  for (dia = sys->dia; dia; dia = dia->n)
  {
    con = dia->con;

    if (con->kind == CONTACT)
    {
      double R [3], U [3],
	    *R0 = dia->R,
	    *U0 = dia->U;

      if (alpha > 0.0)
      {
	double *RE = dia->RE,
	       *DR = dia->DR,
	       *DU = dia->DU;

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

	real_F (res, fri, gap, step, dynamic, epsilon, V, U, F);
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
	    bound = fri * MAX (0, dia->RN);
	    break;
	}

	H [0] = MAX (bound, norm)*R[0] - bound*d[0];
	H [1] = MAX (bound, norm)*R[1] - bound*d[1];
	H [2] = R [2] - MAX (0, d[2]);
      }

      value += DOT (H, H);
    }
  }
  E("LINRUN");

#if MPI
  S("LINCOM");
  double val_i = value;
  MPI_Allreduce (&val_i, &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  E("LINCOM");
#endif

  return value;
}

/* advance solution R = R + alpha * DR */
void LINSYS_Advance (LINSYS *sys, double alpha)
{
  INIT_TIMERS (sys);
  DIAT *dia;

  S("LINRUN");
  for (dia = sys->dia; dia; dia = dia->n)
  {
    double *RE = dia->RE,
	   *DR = dia->DR,
	   *DU = dia->DU,
	   *R = dia->R,
	   *U = dia->U;

    ADDMUL (R, alpha, DR, R);
    ADDMUL (U, alpha, DU, U);
    ADD (U, RE, U);
  }
  E("LINRUN");
}

/* solve A x = b, where b = A [1, 1, ..., 1]' and return |x - [1, 1, ..., 1]| / |[1, 1, ..., 1]| */
double LINSYS_Test (LINSYS *sys, double abstol, int maxiter)
{
  double *x, *y, error;
  VECT *vx, *vb;

  vx = CreateVector (sys->x);
  vb = CreateVector (sys->b);

  for (x = vx->x, y = x + vx->n; x < y; x ++) *x = 1.0;

  A_times_x_equals_y (sys, vx->x, vb->x);

  for (x = vx->x, y = x + vx->n; x < y; x ++) *x = 0.0;

  {
    hypre_FlexGMRESFunctions *gmres_functions;
    void *gmres_vdata;

    gmres_functions = hypre_FlexGMRESFunctionsCreate (CAlloc, Free, CommInfo, CreateVector, CreateVectorArray, DestroyVector,
      MatvecCreate, Matvec, MatvecDestroy, InnerProd, CopyVector, ClearVector, ScaleVector, Axpy, PrecondSetup, Precond);
    gmres_vdata = hypre_FlexGMRESCreate (gmres_functions);

    hypre_FlexGMRESSetTol (gmres_vdata, 0.0);
    hypre_FlexGMRESSetMaxIter (gmres_vdata, maxiter);
    hypre_FlexGMRESSetAbsoluteTol (gmres_vdata, abstol);
    hypre_FlexGMRESSetup (gmres_vdata, sys, vb, vx);
    hypre_FlexGMRESSolve (gmres_vdata, sys, vb, vx);
    hypre_FlexGMRESGetNumIterations (gmres_vdata , &sys->iters);
    hypre_FlexGMRESGetFinalRelativeResidualNorm (gmres_vdata, &sys->resnorm);
    hypre_FlexGMRESDestroy (gmres_vdata);
  }

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
/* update computed here external reactions */
void LINSYS_Update_External_Reactions (LINSYS *sys)
{
  double *R, *x, *y = sys->b->x;
  INIT_TIMERS (sys);
  DIAT *dia;

  S("LINCOM");
  for (dia = sys->dia; dia; dia = dia->n)
  {
    x = &y [3*dia->num];
    R = dia->R;
    COPY (R, x);
  }

  update_external_reactions (sys, y);
  E("LINCOM");
}
#endif

/* destroy linear system */
void LINSYS_Destroy (LINSYS *sys)
{
  MEM_Release (&sys->mapmem);
  MEM_Release (&sys->diamem);
  MEM_Release (&sys->offmem);
  MEM_Release (&sys->ofxmem);
#if !MPI
  if (sys->A) MX_Destroy (sys->A);
#endif
  DestroyVector (sys->x);
  DestroyVector (sys->b);
#if MPI
  COMALL_Free (sys->pattern);
  free (sys->send);
  free (sys->recv);
  free (sys->R);
  free (sys->basenum);
#endif
  free (sys);
}
