/*
 * scf.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * smoothed contact formulation
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
#include "scf.h"
#include "alg.h"
#include "err.h"

/* C(U,R) + X dU + Y dR, where (U,R) = F(U) + m(R - F(U)), F(U) = [UT, UN + mu|UT|] and m(S) = project-on-polar-cone (S) */
void SCF_Linearize (CON *con, double *U, double *R, double UT, double smoothing_omega, double *C, double *X, double *Y)
{
  DOM *dom = con->master->dom;
  short dynamic = dom->dynamic;
  double *V = con->V,
	  gap = con->gap,
	  step = dom->step,
	  fri = con->mat.base->friction,
	  res = con->mat.base->restitution,
	  coh = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area,
	  dF [9], S [3], F [3], m [3];
#define eps smoothing_omega 
  double udash, ulen, sdot, slen, l1, l2, u1[3], u2[3], eps2,
       fri2, onefri2, sq1, sq2, g1, g2, dg1, dg2, a, b, c, d;

  eps2 = eps*eps;
  fri2 = fri*fri;
  onefri2 = 1.0 + fri2;
  if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
  else udash = ((MAX(gap, 0)/step) + U[2]);
  ulen = sqrt (DOT2(U, U) + eps2);

  F [0] = U[0];
  F [1] = U[1];
  if (UT >= 0.0) F [2] = (udash + fri * UT);
  else F [2] = (udash + fri * (ulen - eps));

  SUB (R, F, S);
  S [2] += coh;

  sdot = DOT2 (S, S);
  slen = sqrt (sdot);
  l1 = -(S[2] + fri*slen) / onefri2;
  l2 =  (slen - fri*S[2]) / onefri2;
  if (slen != 0.0)
  {
    u2[0] = S[0]/slen;
    u2[1] = S[1]/slen;
    u2[2] = -fri;
    u1[0] = -fri*u2[0];
    u1[1] = -fri*u2[1];
    u1[2] = -1.0;
  }
  else
  {
    u2[0] =  1.0;
    u2[1] =  0.0;
    u2[2] = -fri;
    u1[0] = -fri;
    u1[1] =  0.0;
    u1[2] = -1.0;
  }
  sq1 = sqrt (l1*l1 + 4.0*eps2);
  sq2 = sqrt (l2*l2 + 4.0*eps2*fri2);
  g1 = 0.5 * (sq1 + l1);
  g2 = 0.5 * (sq2 + l2);

  m [0] = g1*u1[0] + g2*u2[0];
  m [1] = g1*u1[1] + g2*u2[1];
  m [2] = g1*u1[2] + g2*u2[2];

  ASSERT_DEBUG (C, "C needs to be not NULL in VIC_Linearize");
  ADD (F, m, C);

  if (X)
  {
    ASSERT_DEBUG (Y, "X needs to be accompanied by Y in VIC_Linearize");

    dF [1] = dF [3] = dF [6] = dF [7] = 0.0;
    dF [0] = dF [4] = dF [8] = 1.0;
    dF [2] = fri * U[0] / ulen;
    dF [5] = fri * U[1] / ulen;

    dg1 = 0.5*(1.0+l1/sq1);

    if (slen != 0.0)
    {
      dg2 = 0.5*(1.0+l2/sq2);

      /* a = 0.5*(1.0+(sq2 - fri*sq1)/(l2 - fri*l1)), and if
       * (...)/(...) is multiplied by (sq2 + fri*sq1) we get: */
      a = 0.5*(1.0+(l2 + fri*l1)/(sq2 + fri*sq1));
      b = (fri2 * dg1 + dg2) / onefri2;
      c = (fri * (dg1 - dg2)) / onefri2;
      d = (dg1 + fri2 * dg2) / onefri2;

      Y [0] = a + (b - a) * S[0]*S[0] / sdot;
      Y [1] = (b - a) * S[1]*S[0] / sdot;
      Y [2] = c * S[0] / slen;
      Y [3] = Y[1];
      Y [4] = a + (b - a) * S[1]*S[1] / sdot;
      Y [5] = c * S[1] / slen;
      Y [6] = Y[2];
      Y [7] = Y[5];
      Y [8] = d;
    }
    else
    {
      Y[1] = Y[2] = Y[3] = Y[5] = Y[6] = Y[7] = 0.0;
      Y[0] = Y[4] = Y[8] = dg1;
    }

    NNMUL (Y, dF, X);
    NNSUB (dF, X, X); /* X = [I - dm/dS] dF/dU */
  }
}

/* R = project-on-friction-cone (S) */
void SCF_Project (double friction, double cohesion, double *S, double *R)
{
  double m [3], slen, l1, l2, u1[3], u2[3], g1, g2, fri2;
#define fri friction

  S [2] += cohesion;

  fri2 = fri*fri;
  slen = LEN2 (S);
  l1 = -(S[2] + fri*slen) / (1.0 + fri2);
  l2 =  (slen - fri*S[2]) / (1.0 + fri2);
  if (slen != 0.0)
  {
    u2[0] = S[0]/slen;
    u2[1] = S[1]/slen;
    u2[2] = -fri;
    u1[0] = -fri*u2[0];
    u1[1] = -fri*u2[1];
    u1[2] = -1.0;
  }
  else
  {
    u2[0] =  1.0;
    u2[1] =  0.0;
    u2[2] = -fri;
    u1[0] = -fri;
    u1[1] =  0.0;
    u1[2] = -1.0;
  }
  g1 = MAX (l1, 0.0);
  g2 = MAX (l2, 0.0);
  m [0] = g1*u1[0] + g2*u2[0];
  m [1] = g1*u1[1] + g2*u2[1];
  m [2] = g1*u1[2] + g2*u2[2];

  S [2] -= cohesion;

  SUB (S, m, R);
}
