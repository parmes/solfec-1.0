/*
 * nts.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Newton constraints solver
 */

#include <stdlib.h>
#include <time.h>
#include "dom.h"
#include "bgs.h"
#include "nts.h"
#include "alg.h"
#include "err.h"

#define MEMBLK 256

#if 0
static double* create_system_matrix (MEM *mapmem, LOCDYN *ldy)
{
  DIAB *dia;
  OFFB *blk;
  DOM *dom;
  CON *con;

  dom->ldy;

  DOM_Number_Constraints (dom);

  for (con = dom->con; con; con->next)
  {
    dia = con->dia;
  }
}

static void destroy_system_matrix (MEM *mapmem, LOCDYN *ldy)
{
}
#endif

#if 0
/* compute linear system matrix and vector (the matrix 'a' is in the compressed row format) */
static void compute_matrix_and_vector (NTVARIANT variant, LOCDYN *ldy, double *a, double *b, double *r)
{
  double d [3], norm, lim, udash, step;
  short dynamic, pull;
  DIAB *dia;
  OFFB *blk;
  DOM *dom;

  dom = ldy->dom;
  dynamic = dom->dynamic;
  step = dom->step;

  for (dia = ldy->dia; dia; dia = dia->n, b += 3, r += 3)
  {
    CON *con = dia->con;

    double *R = dia->R,
	   *U = dia->U,
	   *V = dia->V,
	   *W = dia->W,
	   *B = dia->B,
	   rho = dia->rho,
	   gap = con->gap,
	   fri = con->mat.base->friction,
	   res = con->mat.base->restitution;

    int *map = dia->map;

    /* calculate residual */
    NVADDMUL (B, W, R, r);
    SUB (r, U, r);
    for (blk = dia->adj; blk; blk = blk->next)
    {
      double *W = blk->W,
	     *R = blk->dia->R;
      NVADDMUL (r, W, R, r);
    }

#if MPI
    /* TODO: remember to include adjext in all cases */
#endif

    /* normal dashed velocity */
    if (dynamic) udash = (U[2] + restitution * MIN (V[2], 0));
    else udash = ((MAX(gap, 0)/step) + U[2]);

    /* predict reaction */
    d [0] = R[0] - rho * U[0];
    d [1] = R[1] - rho * U[1];
    d [2] = R[2] - rho * udash;

    if (d [2] >= 0)
    {
      /* active normal set */
      pull = 0;

      a [map[0]+2] = W[2];
      a [map[1]+2] = W[5];
      a [map[2]+2] = W[8];

      b [2] = -udash - r[2];

      for (blk = dia->adj; blk; blk = blk->next)
      {
	double *W = blk->W;
	int *map = blk->map;

	a [map[0]+2] = W[2];
	a [map[1]+2] = W[5];
	a [map[2]+2] = W[8];
      }
    }
    else 
    {
      /* inactive normal set */
      pull = 1;

      a [map[0]+2] = 0.0;
      a [map[1]+2] = 0.0;
      a [map[2]+2] = 1.0 * W[8];

      b [2] = -R[2] * W[8];

      for (blk = dia->adj; blk; blk = blk->next)
      {
	int *map = blk->map;

	a [map[0]+2] = 0.0;
	a [map[1]+2] = 0.0;
	a [map[2]+2] = 0.0;
      }
    }

    /* tangential response */

    norm = sqrt (d[0]*d[0]+d[1]*d[1]); /* tangential force value */

    switch ((algo & (NEWT|HYB|FIX)))
    {
      case NEWT:
        lim = fri * MAX (0, d[2]);
	break;
      case HYB:
        lim = fri * MAX (0, R[2]);
	break;
      case FIX:
        lim = fri * MAX (0, dia->R2);
	break;
    }

    if ((algo & NEWT) && pull) goto ZERO_TANG; /* enforce AN = AT + IT */

    if (norm >= lim) /* frictional sliping */
    {
      double F [4], /* matrix associated with the derivative of an Euclidean norm in 2D */
	     M [4], H [4],  /* auxiliary metrices & vectors */
	     delta, alpha, beta, den, len, e; /* auxiliary scalars */

      /* active tangential set */

      if (lim > 0.0 && norm > 0.0) /* non-degenerate case */
      {
	/* references below are linked to the paper "A primal-dual active set
	 * strategy for three-dimensional contact problems with Coulomb friction",
	 * S. Hueber, G. Stadler, B.I. Wohlmuth, to appear in 2007 */

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

	if (algo & UNS) /* nonsymetric reduction */
	{
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

	  a [map[0]]   = H[0] + rho*(M[0]*W[0] + M[2]*W[1]);
	  a [map[0]+1] = H[1] + rho*(M[1]*W[0] + M[3]*W[1]);
	  a [map[1]]   = H[2] + rho*(M[0]*W[3] + M[2]*W[4]);
	  a [map[1]+1] = H[3] + rho*(M[1]*W[3] + M[3]*W[4]);

	  if (algo & FIX)
	  {
	    a [map[2]]   = rho*(M[0]*W[6] + M[2]*W[7]);
	    a [map[2]+1] = rho*(M[1]*W[6] + M[3]*W[7]);

	    b [0] = (fri*(d[0]/norm)*dia->R2 - R[0]) * W[0] - rho*(M[0]*r[0] + M[2]*r[1]); 
	    b [1] = (fri*(d[1]/norm)*dia->R2 - R[1]) * W[4] - rho*(M[1]*r[0] + M[3]*r[1]); 
	  }
	  else /* NEWT, HYB */
	  {
	    a [map[2]]   = rho*(M[0]*W[6] + M[2]*W[7]) - fri*(d[0]/norm) * W[0];
	    a [map[2]+1] = rho*(M[1]*W[6] + M[3]*W[7]) - fri*(d[1]/norm) * W[4];

	    b [0] = (fri*(d[0]/norm)*R[2] - R[0]) * W[0] - rho*(M[0]*r[0] + M[2]*r[1]);
	    b [1] = (fri*(d[1]/norm)*R[2] - R[1]) * W[4] - rho*(M[1]*r[0] + M[3]*r[1]);
	  }

	  for (blk = dia->adj; blk; blk = blk->next)
	  {
	    double *W = blk->W;
	    int *map = blk->map;

	    a [map[0]]   = rho*(M[0]*W[0] + M[2]*W[1]);
	    a [map[0]+1] = rho*(M[1]*W[0] + M[3]*W[1]);
	    a [map[1]]   = rho*(M[0]*W[3] + M[2]*W[4]);
	    a [map[1]+1] = rho*(M[1]*W[3] + M[3]*W[4]);
	    a [map[2]]   = rho*(M[0]*W[6] + M[2]*W[7]); 
	    a [map[2]+1] = rho*(M[1]*W[6] + M[3]*W[7]);  
	  }
	}
	else /* "symetric" reduction (symetric in the limit of iterates) */
	{
	  double T [2];
	  int perm [2];

	  ASSERT (lapack_dgetrf (2, 2, M, 2, perm) == 0, "Frictional tangent => inversion of M failed");
	  ASSERT (lapack_dgetri (2, M, 2, perm, T, 2) == 0, "Frictional tangent => inversion of M failed");

	  H [0] = fri*(d[0]/norm);
	  H [1] = fri*(d[1]/norm);
	  T [0] = (M[0]*H[0] + M[2]*H[1])/rho;
	  T [1] = (M[1]*H[0] + M[3]*H[1])/rho;

	  H [0] = (M[0] - 1.0)/rho;
	  H [1] = M[1]/rho;
	  H [2] = M[2]/rho;
	  H [3] = (M[3] - 1.0)/rho;

	  a [map[0]]   = W[0] + H[0];
	  a [map[0]+1] = W[1] + H[1];
	  a [map[1]]   = W[3] + H[2];
	  a [map[1]+1] = W[4] + H[3];

	  if (algo & FIX)
	  {
	    a [map[2]]   = W[6];
	    a [map[2]+1] = W[7];

	    b [0] = T[0]*dia->R2 - (M[0]*R[0] + M[2]*R[1])/rho - r[0];
	    b [1] = T[1]*dia->R2 - (M[1]*R[0] + M[3]*R[1])/rho - r[1];
	  }
	  else /* NEWT, HYB */
	  {
	    a [map[2]]   = W[6] - T[0];
	    a [map[2]+1] = W[7] - T[1];

	    b [0] = T[0]*R[2] - (M[0]*R[0] + M[2]*R[1])/rho - r[0];
	    b [1] = T[1]*R[2] - (M[1]*R[0] + M[3]*R[1])/rho - r[1];
	  }

	  for (blk = dia->adj; blk; blk = blk->next)
	  {
	    double *W = blk->W;
	    int *map = blk->map;

	    a [map[0]]   = W[0];
	    a [map[0]+1] = W[1];
	    a [map[1]]   = W[3];
	    a [map[1]+1] = W[4];
	    a [map[2]]   = W[6];
	    a [map[2]+1] = W[7];
	  }
	}
      }
      else /* degenerate case => enforce homogenous tangential tractions.
	      this is discussed in the Section 3.3 of the HSW paper */
      {
ZERO_TANG: /* this label is for the NEWT jump above (dirty) */
	
	a [map[0]]   = 1.0 * W[0];
	a [map[0]+1] = 0.0;
	a [map[1]]   = 0.0;
	a [map[1]+1] = 1.0 * W[4];
	a [map[2]]   = 0.0;
	a [map[2]+1] = 0.0;

	b [0] = -R[0] * W[0];
	b [1] = -R[1] * W[4];

	for (blk = dia->adj; blk; blk = blk->next)
	{
	  int *map = blk->map;

	  a [map[0]]   = 0.0;
	  a [map[0]+1] = 0.0;
	  a [map[1]]   = 0.0;
	  a [map[1]+1] = 0.0;
	  a [map[2]]   = 0.0;
	  a [map[2]+1] = 0.0;
	}
      }
    }
    else /* frictional sticking */
    {
      /* inactive tangential set */

      a [map[0]]   = W[0];
      a [map[0]+1] = W[1];
      a [map[1]]   = W[3];
      a [map[1]+1] = W[4];

      if ((algo & HYB) ||
	  (algo & FIX))
      {
	a [map[2]]   = W[6];
	a [map[2]+1] = W[7];

	b [0] = -U[0] - r[0];
	b [1] = -U[1] - r[1];
      }
      else /* NEWT */
      {
	a [map[2]]   = W[6]+U[0]/d[2];
	a [map[2]+1] = W[7]+U[1]/d[2];

	b [0] = -(1.0 + rho*udash/d[2])*U[0] - r[0];
	b [1] = -(1.0 + rho*udash/d[2])*U[1] - r[1];
      }

      for (blk = dia->adj; blk; blk = blk->next)
      {
	double *W = blk->W;
	int *map = blk->map;

	a [map[0]]   = W[0];
	a [map[0]+1] = W[1];
	a [map[1]]   = W[3];
	a [map[1]+1] = W[4];
	a [map[2]]   = W[6];
	a [map[2]+1] = W[7];
      }
    }
  }
}
#endif

/* create solver */
NEWTON* NEWTON_Create (NTVARIANT variant, double epsilon, int maxiter)
{
  NEWTON *nt;

  ERRMEM (nt = malloc (sizeof (NEWTON)));

  nt->variant = variant;
  nt->epsilon = epsilon;
  nt->maxiter = maxiter;

  MEM_Init (&nt->mapmem, sizeof (int [3]), MEMBLK);

  return NULL;
}

/* run solver */
void NEWTON_Solve (NEWTON *nt, LOCDYN *ldy)
{
}

/* write labeled satate values */
void NEWTON_Write_State (NEWTON *nt, PBF *bf)
{
}

/* destroy solver */
void NEWTON_Destroy (NEWTON *nt)
{
  MEM_Release (&nt->mapmem);
  free (nt);
}
