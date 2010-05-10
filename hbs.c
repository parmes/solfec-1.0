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
  double error, merit;
#if MPI
  NEWTON *nt;
  MEM setmem;
  SET *subset;
  CON *con;

  MEM_Init (&setmem, sizeof (SET), BLOCKS);

  subset = NULL;

  for (con = dom->con; con; con = con->next)
  {
    if (con->ext) SET_Insert (&setmem, &subset, con, NULL);
  }

  nt = NEWTON_Subset_Create (SMOOTHED_VARIATIONAL, ldy, subset, 1.0, 3, 1E-8);
#endif

  for (int i = 0; i < 5; i ++)
  {
    error = gauss_seidel (ldy, dom->dynamic, dom->step, 5);

#if MPI
    NEWTON_Solve (nt, ldy);

    LINSYS_Update_External_Reactions (nt->sys);
#endif

    merit =  MERIT_Function (ldy, 1);

#if MPI
    if (dom->rank == 0)
#endif
    if (dom->verbose) printf ("HYBRID: iter: %d, merit: %g\n", i, merit);
  }

#if MPI
  NEWTON_Destroy (nt);
  MEM_Release (&setmem);
#endif
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
