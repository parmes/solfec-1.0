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
  GAUSS_SEIDEL *gs;
  short verbose;
  double RNMAX;
  SET *subset;
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

  for (RNMAX *= hs->droptol, subset = NULL, con = dom->con; con; con = con->next) /* pick strong constraints */
  {
    if (fabs (con->R [2]) > RNMAX)
    {
      SET_Insert (&setmem, &subset, con, NULL);
    }
  }
#if MPI
  if (rank == 0)
#endif
  if (verbose) printf ("HYBRID: refining %d%% of constraints ...\n", (100 * SET_Size (subset)) / dom->ncon);

  nt = NEWTON_Subset_Create (SMOOTHED_VARIATIONAL, ldy, subset, 1E-10, hs->refine, 1E-10);  /* refiner */
  nt->verbose = 0;

#if 0
  /* needed only if looping */
  LINSYS_Update_Free_Velocity (nt->sys);
#endif

  NEWTON_Solve (nt, ldy);

#if MPI
  LINSYS_Update_External_Reactions (nt->sys);
#endif

  for (i = 0; i < hs->postsmooth; i ++) /* postsmoothing */
  {
    GAUSS_SEIDEL_Solve (gs, ldy); /* one step */

#if MPI
    if (rank == 0)
#endif
    if (verbose) printf ("HYBRID: postsmoothing step %d, merit: %g\n", i, gs->merhist [0]);

    if (gs->merhist [0] < hs->meritval) goto out; /* done */
  }

out:
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
