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
#include "ext/krylov/krylov.h"

#if MPI
#include "com.h"
#endif

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
      T_ij_mem,
     *T_ii_mem;

  double delta;

  struct vect *x, *b;

  hypre_FlexGMRESFunctions *gmres_functions;
  void *gmres_vdata;

#if MPI
  COMDATA *send, *recv;
  int nsend, nrecv;
  void *pattern;
  double **R; /* external con->R */
#endif
};

#define DR(con) (con)->Y
#define DU(con) ((con)->Y+3)
#define RE(con) ((con)->Y+6)
#define RN(con) ((con)->Y+7) [0]

/* gluing constraint stiffness */
static double glue_stiffness (CON *con)
{
  BULK_MATERIAL *a = con->master->mat,
		*b = con->slave->mat;

  return 2.0 / (1.0/a->young + 1.0/b->young);
}

/* update linearization of a non-contact constraint */
static void system_update_noncontact (CON *con, short dynamic, short nonglue, double step, double *b)
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
      if (nonglue && blk->dia->con->kind == GLUEPNT) continue;

      double *W = blk->W,
	     *T = blk->T;

      NNCOPY (W, T);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      if (nonglue && CON(blk->dia)->kind == GLUEPNT) continue;

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
      if (nonglue && blk->dia->con->kind == GLUEPNT) continue;

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
      if (nonglue && CON(blk->dia)->kind == GLUEPNT) continue;

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

    b [0] = -U[0] - RE[0];
    b [1] = -U[1] - RE[1];
    b [2] = -U[2] - RE[2];

    for (blk = dia->adj; blk; blk = blk->n)
    {
      if (nonglue && blk->dia->con->kind == GLUEPNT) continue;

      double *W = blk->W,
	     *T = blk->T;

      NNCOPY (W, T);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      if (nonglue && CON(blk->dia)->kind == GLUEPNT) continue;

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
  short dynamic, pull, nonglue;
  DIAB *dia;
  OFFB *blk;
  SET *item;
  CON *con;

  nonglue = variant & NON_GLUING;
  dynamic = dom->dynamic;
  step = dom->step;

#if MPI
  DOM_Update_External_Reactions (dom, 0); /* (###) */
#endif

  for (item = SET_First (subset); item; item = SET_Next (subset))
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
	  if (nonglue && blk->dia->con->kind == GLUEPNT) continue;

	  double *W = blk->W,
		 *T = blk->T;

	  T [2] = W[2];
	  T [5] = W[5];
	  T [8] = W[8];
	}
#if MPI
	for (blk = dia->adjext; blk; blk = blk->n)
	{
	  if (nonglue && CON(blk->dia)->kind == GLUEPNT) continue;

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
	    if (nonglue && blk->dia->con->kind == GLUEPNT) continue;

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
	    if (nonglue && CON(blk->dia)->kind == GLUEPNT) continue;

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
	    if (nonglue && blk->dia->con->kind == GLUEPNT) continue;

	    double *T = blk->T;

	    T [0] = T [1] = T [3] = T [4] = T [6] = T [7] = 0;
	  }
#if MPI
	  for (blk = dia->adjext; blk; blk = blk->n)
	  {
	    if (nonglue && CON(blk->dia)->kind == GLUEPNT) continue;

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
	  if (nonglue && blk->dia->con->kind == GLUEPNT) continue;

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
	  if (nonglue && CON(blk->dia)->kind == GLUEPNT) continue;

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
      system_update_noncontact (con, dynamic, nonglue, step, b);
    }
    break;
    }
  }
}

/* real normal to friction cone */
inline static void real_n (double *R, double fri, double *n)
{
  double dot, len;

  dot = DOT2(R, R);
  len = sqrt(dot);

  if (len == 0 || len <= fri * R[2])
  {
    SET (n, 0.0);
  }
  else if (fri * len + R[2] < 0.0)
  {
    dot += R[2]*R[2];
    len = sqrt (dot);
    if (len == 0) { SET (n, 0.0); }
    else { DIV (R, len, n); }
  }
  else
  {
    dot = 1.0 / sqrt (1.0 + fri*fri);
    DIV2 (R, len, n);
    n [2] = -fri;
    SCALE (n, dot);
  }
}

/* imaginary i */
static double imaginary_i;

/* complex normal to friction cone */
inline static void complex_n (double complex *R, double complex fri, double complex *n)
{
  double complex dot, len;

  dot = DOT2(R, R);
  len = csqrt(dot);

  if (creal (len) == 0 || creal (len) <= creal (fri * R[2]))
  {
    SET (n, 0.0 + 0.0 * imaginary_i);
  }
  else if (creal (fri * len + R[2]) < 0.0)
  {
    dot += R[2]*R[2];
    len = csqrt (dot);
    if (creal (len) == 0) { SET(n, 0.0 + 0.0 * imaginary_i); }
    DIV (R, len, n);
  }
  else
  {
    dot = 1.0 / csqrt (1.0 + fri*fri);
    DIV2 (R, len, n);
    n [2] = -fri;
    SCALE (n, dot);
  }
}

/* real normal ray to friction cone */
inline static void real_m (double fri, short smoothing, double *R, double eps, double *m)
{
  double n [3], fun;

  real_n (R, fri, n);
  fun = DOT (R, n);

  if (smoothing == 1 && fun >= 0.0 && fun <= eps)
  {
    fun = ((2.0/eps) - (1.0/(eps*eps))*fun)*(fun*fun);
  }
  else if (smoothing == 2)
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
inline static void complex_m (double complex fri, short smoothing, double complex *R, double complex eps, double complex *m)
{
  double complex n [3], fun;

  complex_n (R, fri, n);
  fun = DOT (R, n);

  if (smoothing == 1 && creal (fun) >= 0.0 && creal (fun) <= creal (eps))
  {
    fun = ((2.0/eps) - (1.0/(eps*eps))*fun)*(fun*fun);
  }
  else if (smoothing == 2)
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
inline static void real_F (CON *con, double step, short dynamic, double *U, double *R, double *F)
{
  DIAB *dia = con->dia;
  double res = con->mat.base->restitution,
	 fri = con->mat.base->friction,
	 gap = con->gap,
	  *V = dia->V,
	 udash;

  if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
  else udash = ((MAX(gap, 0)/step) + U[2]);

  F [0] = U[0];
  F [1] = U[1];
  F [2] = (udash + fri * sqrt (DOT2(U, U)));
}
 
/* complex F = [UT, UN + fri |UT|]' */
inline static void complex_F (CON *con, double step, short dynamic, double complex *U, double complex *R, double complex *F)
{
  DIAB *dia = con->dia;
  double res = con->mat.base->restitution,
	 fri = con->mat.base->friction,
	 gap = con->gap,
	  *V = dia->V;

  double complex udash;

  if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
  else udash = ((MAX(gap, 0)/step) + U[2]);

  F [0] = U[0];
  F [1] = U[1];
  F [2] = (udash + fri * csqrt (DOT2(U, U)));
}

/* update linear system for NONSMOOTH_VARIATIONAL, SMOOTHED_VARIATIONAL variants */
static void system_update_VARIATIONAL (LINVAR variant, DOM *dom, SET *subset, double *rhs)
{
#if 0
  short smooth, dynamic, nonglue;
  double epsilon, step, h;
  SET *item;
  OFFB *blk;
  DIAB *dia;
  CON *con;
#endif

  /* TODO */
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

static void A_times_x_equals_y_HSW_HYBRID_FIXED (LINSYS *sys, double *x, double *y)
{
  double *T, *z;
  short nonglue;
  SET *item;
  CON *con;
  OFFB *blk;
  DIAB *dia;

#if MPI
  update_external_reactions (sys, x); /* (###) */
#endif

  nonglue = sys->variant & NON_GLUING;

  for (item = SET_First (sys->subset); item; item = SET_Next (item), y += 3)
  {
    con = item->data;
    dia = con->dia;
    T = dia->T;
    z = &x [3*con->num];

    y [0] = T[0]*z[0] + T[3]*z[1] + T[6]*z[2];
    y [1] = T[1]*z[0] + T[4]*z[1] + T[7]*z[2];
    y [2] = T[2]*z[0] + T[5]*z[1] + T[8]*z[2];

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
      if (nonglue && con->kind == GLUEPNT) continue;
      T = blk->T;
      z = &x [3*con->num];
      y [0] += T[0]*z[0] + T[3]*z[1] + T[6]*z[2];
      y [1] += T[1]*z[0] + T[4]*z[1] + T[7]*z[2];
      y [2] += T[2]*z[0] + T[5]*z[1] + T[8]*z[2];
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = CON (blk->dia);
      if (nonglue && con->kind == GLUEPNT) continue;
      T = blk->T;
      z = con->R; /* (###) */
      y [0] += T[0]*z[0] + T[3]*z[1] + T[6]*z[2];
      y [1] += T[1]*z[0] + T[4]*z[1] + T[7]*z[2];
      y [2] += T[2]*z[0] + T[5]*z[1] + T[8]*z[2];
    }
#endif
  }
}

static void AT_times_x_equals_y_HSW_HYBRID_FIXED  (LINSYS *sys, double *x, double *y)
{
}

static void A_times_x_equals_y_VARIATIONAL (LINSYS *sys, double *x, double *y)
{
}

static void AT_times_x_equals_y_VARIATIONAL (LINSYS *sys, double *x, double *y)
{
}

static void A_times_x_equals_y (LINSYS *sys, struct vect *vx, struct vect *vy)
{
  switch (LINEARIZATION_VARIANT (sys->variant))
  {
    case NONSMOOTH_HSW:
    case NONSMOOTH_HYBRID:
    case FIXED_POINT:
      A_times_x_equals_y_HSW_HYBRID_FIXED (sys, vx->x, vy->x);
      break;
    case NONSMOOTH_VARIATIONAL:
    case SMOOTHED_VARIATIONAL:
      A_times_x_equals_y_VARIATIONAL (sys, vx->x, vy->x);
      break;
  }
}

static void AT_times_x_equals_y (LINSYS *sys, struct vect *vx, struct vect *vy)
{
  switch (LINEARIZATION_VARIANT (sys->variant))
  {
    case NONSMOOTH_HSW:
    case NONSMOOTH_HYBRID:
    case FIXED_POINT:
      AT_times_x_equals_y_HSW_HYBRID_FIXED (sys, vx->x, vy->x);
      break;
    case NONSMOOTH_VARIATIONAL:
    case SMOOTHED_VARIATIONAL:
      AT_times_x_equals_y_VARIATIONAL (sys, vx->x, vy->x);
      break;
  }
}

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

  Axpy (beta, vy, vy);
  Axpy (delta, vx, vy);
  A_times_x_equals_y (sys, vx, vz);
  AT_times_x_equals_y (sys, vz, vu);
  Axpy (alpha, vu, vy);

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
  return CopyVector (vb, vx);
}
/* GMRES interface end */

/* create linear system resulting from linearization of constraints */
LINSYS* LINSYS_Create (LINVAR variant, LOCDYN *ldy)
{
  DOM *dom = ldy->dom;
  short nonglue;
  LINSYS *sys;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  int ncon;

  ncon = dom->ncon;
  imaginary_i = csqrt (-1);
  nonglue = variant & NON_GLUING;
  ERRMEM (sys = MEM_CALLOC (sizeof (LINSYS)));
  MEM_Init (&sys->setmem, sizeof (SET), MAX (ncon, BLOCKS));
  sys->variant = variant;
  sys->ldy = ldy;

  /* collect constraints */
  for (con = dom->con; con; con = con->next)
  {
    if (nonglue)
    {
      if (con->kind != GLUEPNT)
      {
	SET_Insert (&sys->setmem, &sys->subset, con, NULL);
      }
    }
    else SET_Insert (&sys->setmem, &sys->subset, con, NULL);
  }
  ncon = SET_Size (sys->subset);

  /* constraints local numbering */
  DOM_Number_Constraints (dom, 1, sys->subset);

  /* auxiliary memory pools */
  switch (LINEARIZATION_VARIANT (variant))
  {
    case NONSMOOTH_HSW:
    case NONSMOOTH_HYBRID:
      MEM_Init (&sys->Y_mem, sizeof (double [9]), MAX (ncon, BLOCKS)); /* DR, DU, RE */
      MEM_Init (&sys->T_ij_mem, sizeof (double [9]),  MAX (ncon, BLOCKS)); /* dHi/dRj */
      sys->T_ii_mem = &sys->T_ij_mem;
      break;
    case FIXED_POINT:
      MEM_Init (&sys->Y_mem, sizeof (double [10]), MAX (ncon, BLOCKS)); /* DR, DU, RE, RN */
      MEM_Init (&sys->T_ij_mem, sizeof (double [9]),  MAX (ncon, BLOCKS)); /* dHi/dRj */
      sys->T_ii_mem = &sys->T_ij_mem;
      break;
    case NONSMOOTH_VARIATIONAL:
    case SMOOTHED_VARIATIONAL:
      MEM_Init (&sys->Y_mem, sizeof (double [9]), MAX (ncon, BLOCKS)); /* DR, DU, RE */
      MEM_Init (&sys->T_ij_mem, sizeof (double [9]),  MAX (ncon, BLOCKS)); /* dFi/dRj */
      ERRMEM (sys->T_ii_mem = MEM_CALLOC (sizeof (MEM)));
      MEM_Init (sys->T_ii_mem, sizeof (double [18]),  MAX (ncon, BLOCKS)); /* dFi/dRi, dmi/dRi */
      break;
  }

  /* auxiliary constraint space */
  for (con = dom->con; con; con = con->next)
  {
    if (nonglue && con->kind == GLUEPNT) continue;
    ERRMEM (con->Y = MEM_Alloc (&sys->Y_mem));
    dia = con->dia;
    ERRMEM (dia->T = MEM_Alloc (sys->T_ii_mem));
    for (blk = dia->adj; blk; blk = blk->n)
    {
      if (nonglue && blk->dia->con->kind == GLUEPNT) continue;
      ERRMEM (blk->T = MEM_Alloc (&sys->T_ij_mem));
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      if (nonglue && CON(blk->dia)->kind == GLUEPNT) continue;
      ERRMEM (blk->T = MEM_Alloc (&sys->T_ij_mem));
    }
#endif
  }
#if MPI
  for (MAP *item = MAP_First (dom->conext); item; item = MAP_Next (item))
  {
    con = item->data;
    if (nonglue && con->kind == GLUEPNT) continue;
    ERRMEM (con->Y = MEM_Alloc (&sys->Y_mem));
  }
#endif

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
#endif

  return sys;
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
      system_update_VARIATIONAL (sys->variant, sys->ldy->dom, sys->subset, sys->b->x);
      break;
  }
}

/* solve for reaction increments DR */
void LINSYS_Solve (LINSYS *sys, double abstol, int maxiter)
{
}

/* compute merit function at (R + alpha * DR) */
double LINSYS_Merit (LINSYS *sys, double alpha)
{
  return 0;
}

/* advance solution R = R + alpha * DR; return |DR|/|R| */ 
double LINSYS_Advance (LINSYS *sys, double alpha)
{
  return 0;
}

/* solve A x = b, where b = A [1, 1, ..., 1]' and return |x - [1, 1, ..., 1]| / |[1, 1, ..., 1]| */
double LINSYS_Test (LINSYS *sys, double abstol, int maxiter)
{
  return 0;
}

/* most recent iterations count */
int LINSYS_Iters (LINSYS *sys)
{
  return 0;
}

/* most recent residual norm */
double LINSYS_Resnorm (LINSYS *sys)
{
  return 0;
}

/* destroy linear system */
void LINSYS_Destroy (LINSYS *sys)
{
}
