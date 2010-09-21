/*
 * dbs.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * diagonal block solver
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

#include <stdlib.h>
#include <stdio.h>
#include "alg.h"
#include "pes.h"
#include "dom.h"
#include "lap.h"
#include "dbs.h"
#include "vic.h"
#include "err.h"

static int projected_gradient (short dynamic, double epsilon, int maxiter,
  double step, double friction, double restitution, double gap, double rho,
  double *W, double *B, double *V, double *U, double *R)
{
  double vector [3], scalar; /* auxiliary vector & scalar */
  int iter = 0; /* current iteration counter */
  double UN; /* dashed normal velocity */

  if (dynamic && gap > 0)
  {
    SET (R, 0);
    COPY (B, U);
    return 0;
  }

  do
  {
    /* store current
     * reaction */
    COPY (R, vector);

    /* update velocity */
    NVADDMUL (B, W, R, U);

    /* compute dashed normal velocity */
    if (dynamic) UN = (U[2] + restitution * MIN (V[2], 0));
    else UN = ((MAX(gap, 0)/step) + U[2]);

    /* predict new
     * reactions */
    R [0] -= rho * U[0];
    R [1] -= rho * U[1];
    R [2] -= rho * UN;
   
    /* project normal reaction
     * into its feasible domain */ 
    R [2] = MAX (0, R [2]);

    /* project tangential reaction
     * into the friction cone section */
    scalar = sqrt (R[0]*R[0]+R[1]*R[1]);
    if (scalar >= friction * R[2])
    {
      if (scalar > 0.0)
	scalar = friction * R[2] / scalar;

      R [0] *= scalar;
      R [1] *= scalar;
    }
    
    SUB (R, vector, vector); /* absolute difference */
    scalar = DOT (R, R); /* length of current solution */
    scalar = sqrt (DOT (vector, vector)/MAX (scalar, 1.0));
  }
  while (++ iter < maxiter && scalar > epsilon);

  return iter;
}

static int de_saxce_feng (short dynamic, double epsilon, int maxiter,
  double step, double friction, double restitution, double gap, double rho,
  double *W, double *B, double *V, double *U, double *R)
{
  double vector [3], scalar; /* auxiliary vector & scalar */
  int iter = 0; /* current iteration counter */
  double tau [3], UN;

  if (dynamic && gap > 0)
  {
    SET (R, 0);
    COPY (B, U);
    return 0;
  }

  do
  {
    /* store current
     * reaction */
    COPY (R, vector);

    /* update velocity */
    NVADDMUL (B, W, R, U);

    /* compute dashed normal velocity */
    if (dynamic) UN = (U[2] + restitution * MIN (V[2], 0));
    else UN = ((MAX(gap, 0)/step) + U[2]);

    /* predict new
     * reactions */
    tau [0] = R[0] - rho * U[0];
    tau [1] = R[1] - rho * U[1];
    tau [2] = R[2] - rho * (UN + friction * LEN2 (U));
 
    /* project onto friction cone */ 
    VIC_Project (friction, tau, R);

    SUB (R, vector, vector); /* absolute difference */
    scalar = DOT (R, R); /* length of current solution */
    scalar = sqrt (DOT (vector, vector)/MAX (scalar, 1.0));
  }
  while (++ iter < maxiter && scalar > epsilon);

  return iter;
}

static int semismooth_newton (short dynamic, double epsilon, int maxiter,
  double step, double friction, double restitution, double gap, double rho,
  double *W, double *B, double *V, double *U, double *R)
{
  double RES [3], UN, norm, lim, a [9], b [3], c [3], d [3], R0 [3], error;
  int divi, ipiv [3], iter;

  if (dynamic && gap > 0)
  {
    SET (R, 0);
    COPY (B, U);
    return 0;
  }

  divi = maxiter / 10;
  iter = 0;
  do
  {
    /* store current
     * reaction */
    COPY (R, R0);

    if (dynamic) UN = (U[2] + restitution * MIN (V[2], 0));
    else UN = ((MAX(gap, 0)/step) + U[2]);

    /* predict new
     * reactions */
    d [0] = R[0] - rho * U[0];
    d [1] = R[1] - rho * U[1];
    d [2] = R[2] - rho * UN;

    /* calculate residum RES = W*R + B - U */
    NVADDMUL (B, W, R, RES);
    SUB (RES, U, RES);

    if (d [2] >= 0)
    {
      norm = sqrt (d[0]*d[0]+d[1]*d[1]); /* tangential force value */
      lim = friction * MAX (0, d[2]); /* friction limit */

      if (norm >= lim) /* frictional sliping */
      {
	double F [4], /* matrix associated with the derivative of an Euclidean norm in 2D */
	       M [4], H [4],  /* auxiliary metrices & vectors */
	       delta, alfa, beta, den, len, e; /* auxiliary scalars */


	if (lim > 0.0) /* non-degenerate case */
	{
	  len = sqrt (R[0]*R[0]+R[1]*R[1]);
	  den = MAX (lim, len) * norm;
	  e = lim / norm;
	  if (len == 0.0) beta = 1.0;
	  else
	  {
	    alfa = (R[0]*d[0]+R[1]*d[1]) / (len*norm);
	    delta = MIN (len/lim, 1.0);
	    beta = (alfa < 0.0 ? 1.0 / (1.0 - alfa*delta) : 1.0); /* relaxation factor in case of direction change */
	  }

	  F [0] = (R[0]*d[0])/den;
	  F [1] = (R[1]*d[0])/den;
	  F [2] = (R[0]*d[1])/den;
	  F [3] = (R[1]*d[1])/den;

	  M [0] = e * (1.0 - F[0]);
	  M [1] = - e * F[1];
	  M [2] = - e * F[2];
	  M [3] = e * (1.0 - F[3]);

	  H [0] = 1.0 - beta * M[0];
	  H [1] = - beta * M[1];
	  H [2] = - beta * M[2];
	  H [3] = 1.0 - beta * M[3];

	  a [0] = H[0] + rho*(M[0]*W[0] + M[2]*W[1]);
	  a [1] = H[1] + rho*(M[1]*W[0] + M[3]*W[1]);
	  a [2] = W[2];
	  a [3] = H[2] + rho*(M[0]*W[3] + M[2]*W[4]);
	  a [4] = H[3] + rho*(M[1]*W[3] + M[3]*W[4]);
	  a [5] = W[5];
	  a [6] = rho*(M[0]*W[6] + M[2]*W[7]) - friction*(d[0]/norm);
	  a [7] = rho*(M[1]*W[6] + M[3]*W[7]) - friction*(d[1]/norm);
	  a [8] = W[8];

	  b [0] = friction*(d[0]/norm)*R[2] - R[0] - rho*(M[0]*RES[0] + M[2]*RES[1]);
	  b [1] = friction*(d[1]/norm)*R[2] - R[1] - rho*(M[1]*RES[0] + M[3]*RES[1]);
	  b [2] = -UN - RES[2];
	}
	else /* degenerate case => enforce homogenous tangential tractions */
	{
	  a [0] = 1.0;
	  a [1] = 0.0;
	  a [2] = W[2];
	  a [3] = 0.0;
	  a [4] = 1.0;
	  a [5] = W[5];
	  a [6] = 0.0;
	  a [7] = 0.0;
	  a [8] = W[8];

	  b [0] = -R[0] - RES[0];
	  b [1] = -R[1] - RES[1];
	  b [2] = -UN - RES[2];
	}
      }
      else /* frictional sticking */
      {
	a [0] = W[0];
	a [1] = W[1];
	a [2] = W[2];
	a [3] = W[3];
	a [4] = W[4];
	a [5] = W[5];
	a [6] = W[6]+U[0]/d[2];
	a [7] = W[7]+U[1]/d[2];
	a [8] = W[8];

	b [0] = -(1.0 + rho*U[2]/d[2])*U[0] - RES[0];
	b [1] = -(1.0 + rho*U[2]/d[2])*U[1] - RES[1];
	b [2] = -UN - RES[2];
      }
    }
    else 
    {
      a [0] = 1.0;
      a [1] = 0.0;
      a [2] = 0.0;
      a [3] = 0.0;
      a [4] = 1.0;
      a [5] = 0.0;
      a [6] = 0.0;
      a [7] = 0.0;
      a [8] = 1.0;

      b [0] = -R[0];
      b [1] = -R[1];
      b [2] = -R[2];
    }

    if (lapack_dgesv (3, 1, a, 3, ipiv, b, 3)) return -1;
    if (!isfinite (b[0]+b[1]+b[2])) return -1;

    NVADDMUL (RES, W, b, c);
    ADD (R, b, R);
    ADD (U, c, U);

    SUB (R, R0, R0); 
    error = DOT (R, R);
    error = sqrt (DOT (R0, R0) / MAX (error, 1.0));
    iter ++;

    if ((iter % divi) == 0)
    {
      rho *= 10.0; /* penalty scaling */
      if (isinf (rho)) return -1;
    }
  }
  while (iter < maxiter && error > epsilon);

  return iter;
}

static int fixpnt (short dynamic, double *W, double *B, double *V, double *U, double *R)
{
  double A [9];

  if (dynamic)
  {
    NNCOPY (W, A);
    COPY (V, U);
    SCALE (U, -1.0);
    SUB (U, B, R);
    if (lapack_dposv ('U', 3, 1, A, 3, R, 3)) return -1;
  }
  else
  {
    NNCOPY (W, A);
    SET (U, 0.0);
    SUB (U, B, R);
    if (lapack_dposv ('U', 3, 1, A, 3, R, 3)) return -1;
  }

  return 0;
}

static int fixdir (short dynamic, double *W, double *B, double *V, double *U, double *R)
{
  if (dynamic)
  {
    R [0] = 0.0;
    R [1] = 0.0;
    R [2] = -(V[2] + B[2]) / W[8];
    U [0] =  B [0]; 
    U [1] =  B [1]; 
    U [2] = -V [2];
  }
  else
  {
    R [0] = 0.0;
    R [1] = 0.0;
    R [2] = -B[2] / W[8];
    U [0] = B [0]; 
    U [1] = B [1]; 
    U [2] = 0.0;
  }

  return 0;
}

static int velodir (double *Z, double *W, double *B, double *U, double *R)
{
  R [0] = 0.0;
  R [1] = 0.0;
  R [2] = (VELODIR(Z) - B[2]) / W[8];
  U [0] = B [0]; 
  U [1] = B [1]; 
  U [2] = VELODIR(Z);

  return 0;
}

static int riglnk (short dynamic, double epsilon, int maxiter, double step,
  double *base, double *Z, double *W, double *B, double *V, double *U, double *R)
{
  double X [2],
	 Y [3],
         C [3],
	 D [4],
	 L [4],
	 l, len, error;
  int ipiv [4], iter;

  /* dynamic: q(n+1) = q(n) + (h/2) * (u(n) + u(n+1))  = q(n+1/2) + (h/2) * u(n+1);
   * static: q(n+1) = q(n) + h * u(n+1);
   * minimize 0.5 WNN RN^2 , subject to h(R) = |Z + h([WTN;WNN] RN + B)|^2 - d^2 = 0 */

  TVMUL (base, RIGLNK_VEC(Z), Y); /* local Z */
  if (dynamic) step *= 0.5;
  len = RIGLNK_LEN(Z);
  SET2 (R, 0);
  iter = 0;
  l = 0;
  do
  {
    /* gradient of h(R) = 2(Z + h([WTN; WNN] RN + B))' h [WTN; WNN] */

    ADDMUL (B, R[2], W+6, C);
    ADDMUL (Y, step, C, D);
    X [1] = DOT (D, D) - len*len;
    D [3] = DOT (W+6, D) * 2.0 * step;
    X[0] = W[8]*R[2] + l*D[3];
    L[0] = W[8]; L[2] = D[3];
    L[1] = D[3]; L[3] = 0.0;

    if (lapack_dgesv (2, 1, L, 2, ipiv, X, 2)) return -1;

    R[2] -= X[0];
    l -= X [1];

    error = sqrt (DOT2(X,X)/(R[2]*R[2] + l*l));

  } while (error > epsilon && ++iter < maxiter);

  if (dynamic && iter >= maxiter) /* fall back on the explicit velocity formula */
  {
    R [2] = (U[2] - B[2]) / W[8];
    iter = 0;
  }

  NVADDMUL (B, W, R, U);
 
  return iter;
}

/* diagsolver: diagonal solver kind
 * diagepsilon: relative accuracy on termination
 * diagmaxiter: maximal iterations count
 * dynamic: simulation kind (dom->dynamic)
 * step: time step (dom->step)
 * kind: constraint kind (con->kind)
 * mat: surface material (kind == CONTACT)
 * gap: constraint gap
 * Z: auxiliary Z storage (con->Z)
 * base: constraint local base (con->base)
 * dia: diagonal block of local dynamic (con->dia)
 * B: local free velocity (B = dia->B + sum [dia->adj] (W_i R_i));
 * diagonal block solver */
int DIAGONAL_BLOCK_Solver (DIAS diagsolver, double diagepsilon, int diagmaxiter,
  short dynamic, double step, short kind, SURFACE_MATERIAL *mat, double gap,
  double *Z, double *base, DIAB *dia, double *B)
{
  switch (kind)
  {
  case CONTACT:
    switch (mat->model)
    {
    case SIGNORINI_COULOMB:
      switch (diagsolver)
      {
      case DS_PROJECTED_GRADIENT:
	return projected_gradient (dynamic, diagepsilon, diagmaxiter, step, mat->friction,
			   mat->restitution, gap, dia->rho, dia->W, B, dia->V, dia->U, dia->R);
      case DS_DE_SAXCE_FENG:
	return de_saxce_feng (dynamic, diagepsilon, diagmaxiter, step, mat->friction,
			 mat->restitution, gap, dia->rho, dia->W, B, dia->V, dia->U, dia->R);
      case DS_SEMISMOOTH_NEWTON:
	return semismooth_newton (dynamic, diagepsilon, diagmaxiter, step, mat->friction,
			  mat->restitution, gap, dia->rho, dia->W, B, dia->V, dia->U, dia->R);
      }
      break;
    case SPRING_DASHPOT:
      {
	CON *con = dia->con;
        return PENALTY_Spring_Dashpot_Contact (con, 1, step, gap, mat->spring, mat->dashpot, mat->friction,
	                                       mat->cohesion, dia->W, dia->B, dia->V, dia->U, dia->R);
      }
      break;
    }
    break;
  case FIXPNT:
  case GLUE:
    return fixpnt (dynamic, dia->W, B, dia->V, dia->U, dia->R);
  case FIXDIR:
    return fixdir (dynamic, dia->W, B, dia->V, dia->U, dia->R);
  case VELODIR:
    return velodir (Z, dia->W, B, dia->U, dia->R);
  case RIGLNK:
    return riglnk (dynamic, diagepsilon, diagmaxiter,
	    step, base, Z, dia->W, B, dia->V, dia->U, dia->R);
  }

  return 0;
}
