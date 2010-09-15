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
#include "vic.h"

/* constraint satisfaction merit function approximately indicates the
 * amount of spurious momentum due to constraint force inaccuracy;
 * update_U != 0 implies that U needs to be computed for current R;
 * (it is assumed that all (also external) reactions are updated) */
double MERIT_Function (LOCDYN *ldy, short update_U)
{
  double step, up, uplo [2], Q [3], P [3];
  short dynamic;
  DIAB *dia;
  OFFB *blk;
  CON *con;

  uplo [0] = 0.0;
  uplo [1] = ldy->free_energy;
  dynamic = ldy->dom->dynamic;
  step = ldy->dom->step;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    con = dia->con;

    double *W = dia->W,
	   *A = dia->A,
	   *B = dia->B,
	   *V = dia->V,
	   *U = dia->U,
	   *R = dia->R;

    if (update_U)
    {
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
    }

    switch (con->kind)
    {
    case CONTACT:
    {
      if (dynamic && con->gap > 0) /* open dynamic contact */
      {
	up = 0.0; /* has zero R regardless of U */
      }
      else
      {
	VIC_Linearize (con, 0, P, NULL, NULL);
	NVMUL (A, P, Q);
	up = DOT (Q, P);
      }
    }
    break;
    case FIXPNT:
    case GLUE:
    {
      if (dynamic) { ADD (U, V, P); }
      else { COPY (U, P); }
      NVMUL (A, P, Q);
      up = DOT (Q, P);
    }
    break;
    case FIXDIR:
    {
      if (dynamic) { P[2] = U[2] + V[2]; }
      else { P[2] = U[2]; }
      Q [2] = A[8] * P[2];
      up = Q[2] * P[2];
    }
    break;
    case VELODIR:
    {
      P [2] = VELODIR(con->Z) - U[2];
      Q [2] = A[8] * P[2];
      up = Q[2] * P[2];
    }
    break;
    case RIGLNK:
    {
      if (dynamic) { P[2] = 2.0*con->gap/step + U[2]; }
      else { P [2] = con->gap/step + U[2]; }
      Q [2] = A[8] * P[2];
      up = Q[2] * P[2];
    }
    break;
    }

    con->merit = up; /* per-constraint merit numerator */
    uplo [0] += up;
  }

#if MPI
  double inp [2] = {uplo [0], uplo [1]};
  MPI_Allreduce (inp, uplo, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); /* sum up */
#endif

  uplo [0] *= 0.5; /* was ommited above: E = 0.5 (AU, U) */
  uplo [1] = (uplo [1] == 0 ? 1 : uplo [1]);

  for (con = ldy->dom->con; con; con = con->next)
  {
    con->merit /= uplo [1]; /* per-constraint merit denominator */
  }

  return uplo [0] / uplo [1];
}
