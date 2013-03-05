/*
 * pes.c
 * Copyright (C) 2007-2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * penalty constraints solver
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

#include "alg.h"
#include "dom.h"
#include "bgs.h"
#include "pes.h"
#include "lap.h"
#include "err.h"

#define PENALTY_MAXITER 1000
#define PENALTY_EPSILON 1E-4 /* XXX */

/* spring and dashpot based explicit diagonal block contact solver */
int PENALTY_Spring_Dashpot_Contact (CON *con, short implicit, double step, double gap, double spring, double dashpot, double hpow,
                              double friction, double cohesion, double *W, double *B, double *V, double *U, double *R)
{
#if 1 
  double INV [4], WTT[4] = {W[0], W[1], W[3], W[4]}, BN, BT [2], det, len, g, s;
  short cohesive = con->state & CON_COHESIVE;

  BN = B[2] + W[2]*R[0] + W[5]*R[1];

  if (implicit)
  {
    g = (gap + 0.25 * step * (BN - V[2]));
    g = MIN (g, 0);
    s = spring * hpow * pow (-g, hpow - 1.0); /* current spring */
    if (dashpot < 0.0) dashpot = 2.0 * sqrt (step * s / W [8]); /* critical damping */
    R [2] = (spring * pow (-g, hpow) - 0.5 * dashpot * (BN + V[2]))
	  / (1.0 + (0.25  * step * s + 0.5 * dashpot) * W[8]);
  }
  else
  {
    g = MIN (gap, 0);
    s = spring * hpow * pow (-g, hpow - 1.0); /* current spring */
    if (dashpot < 0.0) dashpot = 2.0 * sqrt (step * s / W [8]); /* critical damping */
    R [2] = spring * pow (-g, hpow) - dashpot * V[2];
  }

  if (!cohesive && R[2] < 0.0)
  {
    SET (R, 0.0);
    SET (U, 0.0);
    return 0;
  }

  BT [0] = B[0] + W[6] * R[2];
  BT [1] = B[1] + W[7] * R[2];

  INVERT2(WTT, INV, det);

  if (det == 0.0) return -1;

  R [0] = - INV[0] * BT[0] - INV[2] * BT[1];
  R [1] = - INV[1] * BT[0] - INV[3] * BT[1];

  if (cohesive && R [2] < -cohesion * con->area)
  {
    cohesive = 0;
    con->state &= ~CON_COHESIVE;
    SURFACE_MATERIAL_Cohesion_Set (&con->mat, 0.0);
    R [2] = -cohesion * con->area;
  }

  len = LEN2 (R);
  det = friction * fabs (R[2]);

  if (len > det)
  {
    R[0] *= det / len;
    R[1] *= det / len;

    if (cohesive)
    {
      cohesive = 0;
      con->state &= ~CON_COHESIVE;
      SURFACE_MATERIAL_Cohesion_Set (&con->mat, 0.0);
    }

    /* TODO: it would be more "implicit" to solve local nonlinear problem (I + len/|UT| WTT) UT = BT,
     * TODO: but then the local iterations might fail; explicit stuff does not fail like this ... */
  }

  NVADDMUL (B, W, R, U); /* U = B + W R */

  return 0;
#else
  double T [9], Z [3], d [2], BN, rho, trc, det, norm, lim, error;
  int ipiv [3], iter;

  ASSERT_TEXT (cohesion == 0, "Cohesion is not yet handled in the PENALTY_SOLVER");

  BN = B[2] + W[2]*R[0] + W[5]*R[1];

  if (dashpot < 0.0) dashpot = sqrt (step * spring / W [8]); /* critical damping */

  if (implicit)
  {
    R [2] = (- spring * (gap + 0.25 * step * (BN - V[2])) - 0.5 * dashpot * (BN + V[2]))
	  / (1.0 + (0.25  * step * spring + 0.5 * dashpot) * W[8]);
  }
  else
  {
    R [2] = - spring * gap - dashpot * V[2];
  }

  trc = W[0]+W[4];
  det = W[0]*W[4] - W[1]*W[3];
  rho = 1.0 / (0.5*trc + sqrt (0.25*trc-det)); /* 1.0 / larger eigenvalue of W_TT */
  trc = 1.0 / (0.5*spring*step + dashpot);
  //trc = 1.0 / (0.25*spring*step + 0.5*dashpot);

  error = 1.0;
  iter = 0;

  do
  {

    NVADDMUL (B, W, R, U); /* U = B + W R */

    d [0] = R[0] - rho*U[0];
    d [1] = R[1] - rho*U[1];

    norm = LEN2 (d);
    lim = friction * (R[2]+cohesion);

    if (norm < lim)
    {
      NNCOPY (W, T);

      T [6] += U[0] / (R[2]+cohesion);
      T [7] += U[1] / (R[2]+cohesion);
      T [8] += trc;

      Z [0] = -U[0];
      Z [1] = -U[1];
    }
    else if (lim > 0)
    {
      double F [4], M [4], H [4], den, len, e;

      len = sqrt (R[0]*R[0]+R[1]*R[1]);
      den = lim * norm;
      e = lim / norm;

      F [0] = (R[0]*d[0])/den;
      F [1] = (R[1]*d[0])/den;
      F [2] = (R[0]*d[1])/den;
      F [3] = (R[1]*d[1])/den;

      M [0] = e * (1.0 - F[0]);
      M [1] = - e * F[1];
      M [2] = - e * F[2];
      M [3] = e * (1.0 - F[3]);

      H [0] = 1.0 - M[0];
      H [1] = - M[1];
      H [2] = - M[2];
      H [3] = 1.0 - M[3];

      T [0] = H[0] + rho*(M[0]*W[0] + M[2]*W[1]);
      T [1] = H[1] + rho*(M[1]*W[0] + M[3]*W[1]);
      T [2] = W[2];
      T [3] = H[2] + rho*(M[0]*W[3] + M[2]*W[4]);
      T [4] = H[3] + rho*(M[1]*W[3] + M[3]*W[4]);
      T [5] = W[5];
      T [6] = rho*(M[0]*W[6] + M[2]*W[7]) - friction*(d[0]/norm);
      T [7] = rho*(M[1]*W[6] + M[3]*W[7]) - friction*(d[1]/norm);
      T [8] = W[8] + trc;

      Z [0] = friction*(d[0]/norm)*(R[2]+cohesion) - R[0];
      Z [1] = friction*(d[1]/norm)*(R[2]+cohesion) - R[1];
    }
    else /* R_T = 0, R_N = spring * gap(t+h) + dashpot * U(t+h) */
    {
      IDENTITY (T);
      Z [0] = -R[0];
      Z [1] = -R[1];
      T [8] = W[8] + trc;
    }
 
    if (implicit)
    { 
      double r = spring*(gap + 0.5*step*U[2]) + dashpot*U[2];
      //double r = spring*(gap + 0.25*step*(U[2]-V[2])) + 0.5*dashpot*(U[2]+V[2]);
      Z [2] = -trc * ((R[2]+cohesion) -MAX(-r,0));
      ASSERT_TEXT (lapack_dgesv (3, 1, T, 3, ipiv, Z, 3) == 0, "Newton iterations failed in the PENALTY_SOLVER");
    }
    else
    {
      Z [2] = 0.0;
      ASSERT_TEXT (lapack_dgesv (2, 1, T, 3, ipiv, Z, 3) == 0, "Newton iterations failed in the PENALTY_SOLVER");
    }
    ADD (R, Z, R);

    error = DOT (R, R);
    error = sqrt (DOT (Z, Z) / MAX (error, 1.0));
    iter ++;

#if 0
    printf ("Iter = %d, Error = %g\n", iter, error);
#endif
  }
  while (error > 1E-10 && iter < 100);

  if (R[2] < 0.0)
  {
    SET (R, 0.0);
    COPY (B, U);
  }

  return 0;
#endif
}

/* create penalty solver */
PENALTY* PENALTY_Create (short implicit)
{
  PENALTY *ps;

  ERRMEM (ps = MEM_CALLOC (sizeof (PENALTY)));
  ps->implicit = implicit;

  return ps;
}

/* explcit constraint solver */
void PENALTY_Solve (PENALTY *ps, LOCDYN *ldy)
{
  GAUSS_SEIDEL *gs;
  short implicit;
  double step;
  LOCDYN *clo;
  DIAB *dia;
  CON *con;

  implicit = ps->implicit;
  step = ldy->dom->step;

#if MPI
  if (ldy->dom->rank == 0)
#endif
  if (ldy->dom->verbose) printf ("PENALTY_SOLVER: applying springs and dashpots...\n");

  /* first explicitly process contacts */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    con = dia->con;

    if (con->kind == CONTACT)
    {
      SURFACE_MATERIAL *bas = con->mat.base;
      PENALTY_Spring_Dashpot_Contact (con, implicit, step, con->gap, bas->spring, bas->dashpot, bas->hpow,
	                  bas->friction, bas->cohesion, dia->W, dia->B, dia->V, dia->U, dia->R);
    }
  }

#if MPI
  if (ldy->dom->rank == 0)
#endif
  if (ldy->dom->verbose) printf ("PENALTY_SOLVER: solving non-contact constraints...\n");

#if MPI
  DOM_Update_External_Reactions (ldy->dom, 0);
#endif
  clo = LOCDYN_Clone_Non_Contacts (ldy);
  gs = GAUSS_SEIDEL_Create (PENALTY_EPSILON, PENALTY_MAXITER, 1.0, GS_FAILURE_CONTINUE, 1E-9, 100, DS_SEMISMOOTH_NEWTON, NULL, NULL);
  gs->nomerit = 1;
  GAUSS_SEIDEL_Solve (gs, clo);
  GAUSS_SEIDEL_Destroy (gs);
  LOCDYN_Destroy (clo);
}

/* write labeled satate values */
void PENALTY_Write_State (PENALTY *ps, PBF *bf)
{
  /* nothing to write for the moment */
}

/* destroy penalty solver */
void PENALTY_Destroy (PENALTY *ps)
{
  free (ps);
}
