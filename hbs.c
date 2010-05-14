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

/* local subset relaxation solver */
static void gauss_seidel (SET *conset, DOM *dom, double epsilon, int maxiter)
{
  double error, step;
  short dynamic;
  int iters;
  SET *item;
  CON *con;

  dynamic = dom->dynamic;
  step = dom->step;
  iters = 0;
  do
  {
    double errup = 0.0, errlo = 0.0;
    OFFB *blk;
    DIAB *dia;
  
    for (item = SET_First (conset); item; item = SET_Next (item)) 
    {
      con = item->data;
      dia = con->dia;

      double R0 [3],
	     B [3],
	     *R = dia->R;

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
	       *R = CON(blk->dia)->R;
	NVADDMUL (B, W, R, B);
      }
#endif
      
      COPY (R, R0);

      DIAGONAL_BLOCK_Solver (DS_PROJECTED_GRADIENT, 1E-6, 1000,
        dynamic, step, con->kind, con->mat.base, con->gap, con->Z, con->base, dia, B);

      SUB (R, R0, R0);
      errup += DOT (R0, R0);
      errlo += DOT (R, R);
    }

    error = sqrt (errup) / sqrt (MAX (errlo, 1.0));
  }
  while (++ iters < maxiter && error > epsilon);
}

/* create solver */
HYBRID* HYBRID_Create (int presmooth, int refine, int postsmooth, double droptol, double meritval)
{
  HYBRID *hs;

  ERRMEM (hs = malloc (sizeof (HYBRID)));
  hs->presmooth = presmooth;
  hs->refine = refine;
  hs->postsmooth = postsmooth;
  hs->droptol = droptol;
  hs->meritval = meritval;

  return hs;
}

/* run solver */
void HYBRID_Solve (HYBRID *hs, LOCDYN *ldy)
{
#if 0
  SET *subset, *spc;
  GAUSS_SEIDEL *gs;
  short verbose;
  double RNMAX;
  MEM setmem;
  NEWTON *nt;
  DOM *dom;
  CON *con;
  int i;
#if MPI
  int rank;
#endif

  MEM_Init (&setmem, sizeof (SET), BLOCKS);

  nt = NULL;
  dom = ldy->dom;
  verbose = dom->verbose;

#if MPI
  rank = dom->rank;
#endif

  gs = GAUSS_SEIDEL_Subset_Create (ldy, NULL, 1E-3, 1, 1E-3); /* smoother */
  gs->variant = GS_MIDDLE_JACOBI;
  gs->verbose = 0;

  for (i = 0; i < hs->presmooth; i ++) /* presmoothing */
  {
    GAUSS_SEIDEL_Solve (gs, ldy); /* one step */

#if MPI
    if (rank == 0)
#endif
    if (verbose) printf ("HYBRID: presmoothing step %d, merit: %g\n", i, gs->merhist [0]);

    if (gs->merhist [0] < hs->meritval) goto out; /* done */
  }

  for (RNMAX = 0.0, con = dom->con; con; con = con->next) /* find maximal normal reaction */
  {
    if (fabs (con->R [2]) > RNMAX) RNMAX = fabs (con->R [2]);
  }

  for (RNMAX *= hs->droptol, subset = spc = NULL, con = dom->con; con; con = con->next) /* pick strong constraints */
  {
    if (fabs (con->R [2]) > RNMAX)
    {
      SET_Insert (&setmem, &subset, con, NULL);
    }

    if (!con->slave) SET_Insert (&setmem, &spc, con, NULL); /* single point constraints */
  }
#if MPI
  if (rank == 0)
#endif
  if (verbose) printf ("HYBRID: refining %d%% of constraints ...\n", (100 * SET_Size (subset)) / dom->ncon);

  nt = NEWTON_Subset_Create (SMOOTHED_VARIATIONAL, ldy, subset, 1E-10, hs->refine, 1E-10);  /* refiner */
  nt->linmaxiter = 5;
  nt->verbose = 1;

#if 0
  /* needed only if looping */
  LINSYS_Update_Free_Velocity (nt->sys);
#endif

  NEWTON_Solve (nt, ldy);

#if MPI
  LINSYS_Update_External_Reactions (nt->sys);
#endif

  gs->variant = GS_FULL;
  for (i = 0; i < hs->postsmooth; i ++) /* postsmoothing */
  {
    GAUSS_SEIDEL_Solve (gs, ldy); /* one step */

#if MPI
    if (rank == 0)
#endif
    if (verbose) printf ("HYBRID: postsmoothing step %d, merit: %g\n", i, gs->merhist [0]);

    if (gs->merhist [0] < hs->meritval) goto out; /* done */
  }

  gauss_seidel (spc, dom, 1E-4, 100); /* refine single point constraints solution */

out:
  if (nt) NEWTON_Destroy (nt);
  GAUSS_SEIDEL_Destroy (gs);
  MEM_Release (&setmem);
#else

#if MPI
  SET *inner, *outer;
  double merit;
  LINSYS *sys;
  MEM setmem;
  DOM *dom;
  CON *con;
  int i;

  MEM_Init (&setmem, sizeof (SET), BLOCKS);

  dom = ldy->dom;
  inner = outer = NULL;

  for (con = dom->con; con; con = con->next)
  {
    if (con->ext)
    {
      SET_Insert (&setmem, &outer, con, NULL);
    }
    else
    {
      SET_Insert (&setmem, &inner, con, NULL);
    }
  }

  sys = LINSYS_Create (SMOOTHED_VARIATIONAL, ldy, outer);

  for (i = 0; i < hs->refine; i ++)
  {
    LINSYS_Update_Free_Velocity (sys);

    LINSYS_Update (sys);

    LINSYS_Solve (sys, 1E-8, 2);

    LINSYS_Advance (sys, 1.0);

    LINSYS_Update_External_Reactions (sys);

    gauss_seidel (inner, dom, 0, 1);

    gauss_seidel (outer, dom, 0, 1);

    merit = MERIT_Function (ldy, 1);

    if (dom->rank == 0 && dom->verbose) printf ("HYBRID: step %d, merit: %g\n", i, merit);
  }

  MEM_Release (&setmem);
  LINSYS_Destroy (sys);
#endif

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
