/*
 * mrf.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * constraints satisfaction merit function
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

#include "dom.h"
#include "mrf.h"
#include "alg.h"


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

/* real normal ray to friction cone */
inline static void real_m (double fri, double *S, double *m)
{
  double n [3], fun;

  real_n (S, fri, n);
  fun = DOT (S, n);
  MUL (n, fun, m)
}

/* constraint satisfaction merit function approximately indicates the
 * amount of spurious momentum due to constraint force inaccuracy;
 * (it is assumed that all (also external) reactions are updated) */
double MERIT_Function (LOCDYN *ldy)
{
  double step, upper, lower, A [9], det, Q [3], P [3];
  short dynamic;
  DIAB *dia;
  OFFB *blk;

  upper = lower = 0.0;
  dynamic = ldy->dom->dynamic;
  step = ldy->dom->step;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;

    double *W = dia->W,
	   *B = dia->B,
	   *V = dia->V,
	   *U = dia->U,
	   *R = dia->R;

    NVADDMUL (B, W, R, U);
    for (blk = dia->adj; blk; blk = blk->n)
    {
      double *W = blk->W, *R = blk->dia->R;
      NVADDMUL (U, W, R, U);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      double *W = blk->W, *R = CON(blk->dia)->R;
      NVADDMUL (U, W, R, U);
    }
#endif

    INVERT (W, A, det);
    NVMUL (A, B, Q);
    lower += DOT (Q, Q);

    switch (con->kind)
    {
    case CONTACT:
    {
      double res = con->mat.base->restitution,
	     fri = con->mat.base->friction,
	     gap = con->gap,
	     udash, m [3];

      if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
      else udash = ((MAX(gap, 0)/step) + U[2]);

      Q [0] = U[0];
      Q [1] = U[1];
      Q [2] = (udash + fri * LEN2(U));
      SUB (R, Q, P);
      real_m (fri, P, m);
      ADD (Q, m, P);
      NVMUL (A, P, Q);
      upper += DOT (Q, Q);
    }
    break;
    case FIXPNT:
    {
      if (dynamic) { ADD (U, V, P); }
      else { COPY (U, P); }
      NVMUL (A, P, Q);
      upper += DOT (Q, Q);
    }
    break;
    case FIXDIR:
    {
      if (dynamic) { P[2] = U[2] + V[2]; }
      else { P[2] = U[2]; }
      Q [2] = A[8] * P[2];
      upper += Q[2] * Q[2];
    }
    break;
    case VELODIR:
    {
      P [2] = VELODIR(con->Z) - U[2];
      Q [2] = A[8] * P[2];
      upper += Q[2] * Q[2];
    }
    break;
    case RIGLNK:
    {
      if (dynamic) { P[2] = RIGLNK_LEN(con->Z)/step + (U[2]+V[2]); }
      else { P [2] = RIGLNK_LEN(con->Z)/step + U[2]; }
      Q [2] = A[8] * P[2];
      upper += Q[2] * Q[2];
    }
    break;
    case GLUEPNT:
    {
      NVMUL (A, U, Q);
      upper += DOT (Q, Q);
    }
    break;
    }
  }

#if MPI
  double val_i [2] = {lower, upper}, val_o [2];
  MPI_Allreduce (val_i, val_o, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  lower = val_o [0], upper = val_o [1];
#endif

  return sqrt (upper / (lower == 0 ? 1 : lower));
}
