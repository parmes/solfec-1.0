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
#endif

/* create solver */
HYBRID* HYBRID_Create (int refine, int smooth, double droptol, double meritval)
{
  HYBRID *hs;

  ERRMEM (hs = malloc (sizeof (HYBRID)));
  hs->refine = refine;
  hs->smooth = smooth;
  hs->droptol = droptol;
  hs->meritval = meritval;

  return hs;
}

/* run solver */
void HYBRID_Solve (HYBRID *hs, LOCDYN *ldy)
{
  GAUSS_SEIDEL *gs;
  short verbose;
  double RNMAX;
  SET *subset;
  MEM setmem;
  NEWTON *nt;
  DOM *dom;
  CON *con;
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

  for (RNMAX = 0.0, con = dom->con; con; con = con->next) /* find maximal normal reaction */
  {
    if (fabs (con->R [2]) > RNMAX) RNMAX = fabs (con->R [2]);
  }

  for (RNMAX *= hs->droptol, subset = NULL, con = dom->con; con; con = con->next) /* pick strong constraints */
  {
    if (fabs (con->R [2]) > RNMAX || con->kind != CONTACT)
    {
      SET_Insert (&setmem, &subset, con, NULL);
    }
  }

#if MPI
  if (rank == 0)
#endif
  if (verbose) printf ("HYBRID: refining %d%% of constraints ...\n", (100 * SET_Size (subset)) / dom->ncon);

  nt = NEWTON_Subset_Create (SMOOTHED_VARIATIONAL, ldy, subset, 1E-10, hs->refine, 1E-10);  /* refiner */
  nt->linmaxiter = 5;
  nt->verbose = 1;

#if 0
  LINSYS_Update_Free_Velocity (nt->sys);
#endif

  NEWTON_Solve (nt, ldy);

#if MPI
  LINSYS_Update_External_Reactions (nt->sys);
#endif

  gs = GAUSS_SEIDEL_Subset_Create (ldy, NULL, 1.0, hs->smooth, hs->meritval);
  gs->verbose = 1;
  GAUSS_SEIDEL_Solve (gs, ldy);

  if (nt) NEWTON_Destroy (nt);
  GAUSS_SEIDEL_Destroy (gs);
  MEM_Release (&setmem);
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
