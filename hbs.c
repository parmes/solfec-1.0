/*
 * hbs.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * hybrid constraint solver
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
#include <float.h>

#include "hbs.h"
#include "alg.h"
#include "bla.h"
#include "err.h"
#include "bgs.h"
#include "nts.h"
#include "dom.h"
#include "mrf.h"

#define BLOCKS 256


#if 0
static double gauss_seidel (LOCDYN *ldy, short dynamic, double step, int maxiter)
{
  double error;
  int iters;
  DIAB *dia;
  OFFB *blk;

  for (iters = 0; iters < maxiter; iters ++)
  {
    double errup = 0.0, errlo = 0.0;

    for (dia = ldy->dia; dia; dia = dia->n)
    {
      CON *con = dia->con;
#if MPI
      if (con->ext) continue; /* skip exporting constraints */
#endif
      double *R = con->R,
	     R0 [3],
	     B [3];

      /* compute local free velocity */
      COPY (dia->B, B);
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *W = blk->W,
	       *R = blk->dia->R;
	NVADDMUL (B, W, R, B);
      }
#if MPI
      for (blk = dia->adjext; blk; blk = blk->n)
      {
	double *W = blk->W,
	       *R = CON (blk->dia)->R;
	NVADDMUL (B, W, R, B);
      }
#endif
      
      COPY (R, R0); /* previous reaction */

      /* solve local diagonal block problem */
      DIAGONAL_BLOCK_Solver (GS_SEMISMOOTH_NEWTON, 1E-6, 100, dynamic, step, con->kind, con->mat.base, con->gap, con->Z, con->base, dia, B);

      /* accumulate relative
       * error components */
      SUB (R, R0, R0);
      errup += DOT (R0, R0);
      errlo += DOT (R, R);
    }

    /* calculate relative error */
    error = sqrt (errup / (errlo > 0 ? errlo : 1.0));
  }

  return error;
}
#endif

/* create solver */
HYBRID* HYBRID_Create (double epsilon, int maxiter, double meritval)
{
  HYBRID *hs;

  ERRMEM (hs = malloc (sizeof (HYBRID)));
  hs->epsilon = epsilon;
  hs->maxiter = maxiter;
  hs->meritval = meritval;

  return hs;
}

/* run solver */
void HYBRID_Solve (HYBRID *hs, LOCDYN *ldy)
{
  DOM *dom = ldy->dom;
  GAUSS_SEIDEL *gs;
  double maxrn;
  NEWTON *nt;
  MEM setmem;
  SET *subset;
  CON *con;

  MEM_Init (&setmem, sizeof (SET), BLOCKS);

  gs = GAUSS_SEIDEL_Create (1E-6, 3, 1E-6, GS_FAILURE_CONTINUE, 1E-6, 100, GS_PROJECTED_GRADIENT, NULL, NULL);

  GAUSS_SEIDEL_Solve (gs, ldy); /* few Gauss-Seidel iteratios identify strong forces */

  for (maxrn = 0.0, con = dom->con; con; con = con->next) /* find maximal normal reaction */
  {
    if (fabs (con->R [2]) > maxrn) maxrn = fabs (con->R [2]);
  }

  for (subset = NULL, con = dom->con; con; con = con->next) /* pick strong constraints */
  {
    if (fabs (con->R [2]) > 0.05 * maxrn)
    {
      SET_Insert (&setmem, &subset, con, NULL);
    }
  }
#if MPI
  if (dom->rank == 0)
#endif
  if (dom->verbose) printf ("HYBRID: picked %d%% of constraints for Newton\n", (100 * SET_Size (subset)) / dom->ncon);

  nt = NEWTON_Subset_Create (SMOOTHED_VARIATIONAL, ldy, subset, 1.0, 3, 1E-8);

  for (int iter = 0; iter < 1; iter ++)
  {
#if 1
    LINSYS_Update_Free_Velocity (nt->sys);
#endif

    NEWTON_Solve (nt, ldy);

#if MPI
    LINSYS_Update_External_Reactions (nt->sys);
#endif

    GAUSS_SEIDEL_Solve (gs, ldy); /* smooth out solution */
				  /* TODO: variant of GS that avoids initialization in parallel here */
  }

  NEWTON_Destroy (nt);
  MEM_Release (&setmem);
  GAUSS_SEIDEL_Destroy (gs);
}

/* write labeled satate values */
void HYBRID_Write_State (HYBRID *hs, PBF *bf)
{
}

/* destroy solver */
void HYBRID_Destroy (HYBRID *hs)
{
  free (hs);
}
