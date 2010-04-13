/*
 * nts.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Newton constraints solver
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

#if MPI
#include <mpi.h>
#endif

#include <stdlib.h>
#include <time.h>
#include "dom.h"
#include "bgs.h"
#include "nts.h"
#include "alg.h"
#include "bla.h"
#include "err.h"

/* non-monotone globalization merit function value selection */
static double nonmonotone_merit_function (double *values, int length, double merit, int index)
{
  values [index % length] = merit;

  for (int n = 0; n < length; n ++)
    if (values [n] > merit) merit = values [n];

  return merit;
}

/* line search */
static double line_search (LINSYS *sys, double reference, double *merit)
{
  double auxiliary,
	 alpha,
	 gamma;

  int imax,
      iter;

  auxiliary = 2.0 * (*merit);
  alpha = 1.0;
  gamma = 0.1;
  imax = 32;

  for (iter = 0; iter < imax; iter ++)
  {
    (*merit) = LINSYS_Merit (sys, alpha);

    if ((*merit) <= (reference - gamma*alpha*auxiliary)) break;
    else alpha *= 0.9;
  }

  return alpha;
}

/* update constraint reactions */
static double reactions_update (NEWTON *nt, LINSYS *sys, LOCDYN *ldy, double *nonmonvalues, int iter, double *merit)
{
  double reference, alpha;

  if (iter && nt->variant != FIXED_POINT) /* skip line search during the first iteration */
  {
    (*merit) = LINSYS_Merit (sys, 0.0);

    reference = nonmonotone_merit_function (nonmonvalues, nt->nonmonlength, (*merit), iter);

    alpha = line_search (sys, reference, merit);
  }
  else alpha = 1.0;

  return LINSYS_Advance (sys, alpha); /* R = R + alpha * DR, U = U(R) */
}

/* update noraml bounds and return relative error of the update */
static double update_normal_bounds (LINSYS *sys, LOCDYN *ldy)
{
  double errup, errlo;
  CON *con;

  errup = errlo = 0.0;
  for (con = ldy->dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT)
    {
      double DRN,
	     RN;

      DRN = con->R[2] - RN(con);
      RN = RN(con) = con->R[2];

      errup += DRN*DRN;
      errlo += RN*RN;
    }
  }

#if MPI
  if (LINSYS_Global (sys)) /* sum up error */
  {
    double errloc [2] = {errup, errlo}, errsum [2];
    MPI_Allreduce (errloc, errsum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    errup = errsum [0], errlo = errsum [1];
  }
#endif

  return sqrt (errup) / MAX (sqrt (errlo), 1.0);
}

/* create solver */
NEWTON* NEWTON_Create (LINVAR variant, double epsilon, int maxiter, double meritval)
{
  NEWTON *nt;

  ERRMEM (nt = MEM_CALLOC (sizeof (NEWTON)));

  nt->variant = variant;
  nt->epsilon = epsilon;
  nt->maxiter = maxiter;
  nt->meritval = meritval;
  nt->nonmonlength = 10;
  nt->linmaxiter = 1000;
  nt->rerhist = NULL;
  nt->merhist = NULL;

  return nt;
}

/* run solver */
void NEWTON_Solve (NEWTON *nt, LOCDYN *ldy)
{
  double merit, error, *nonmonvalues;
  char fmt [512];
  LINSYS *sys;
  DOM *dom;

  dom = ldy->dom;

  if (dom->verbose)
  {
    sprintf (fmt, "NEWTON: (LIN its/res: %%%dd/%%.2e) iteration: %%%dd  error:  %%.2e  merit: %%.2e\n",
      (int)log10 (nt->linmaxiter) + 1, (int)log10 (nt->maxiter) + 1);
  }

  ERRMEM (nonmonvalues = MEM_CALLOC (nt->nonmonlength * sizeof (double)));
  nt->rerhist = realloc (nt->rerhist, nt->maxiter * sizeof (double));
  nt->merhist = realloc (nt->merhist, nt->maxiter * sizeof (double));

  sys = LINSYS_Create (nt->variant, 0, ldy);
  merit = error = 1.0;
  nt->iters = 0;

#if DEBUG
  error = LINSYS_Test (sys, nt->meritval, nt->linmaxiter);
#if MPI
  if (dom->rank == 0)
#endif
  printf ("NEWTON: W inversion relative accuract: %.2e\n", error);
#endif

  do
  {
    LINSYS_Update (sys); /* assemble A, b */

    LINSYS_Solve (sys, 0.01 * MIN (error, merit), nt->linmaxiter); /* x = A\b; TODO: develop into a rigorous inexact step */

    error = reactions_update (nt, sys, ldy, nonmonvalues, nt->iters, &merit); /* R(i+1) */

    if (nt->variant == FIXED_POINT && error < nt->epsilon) error = update_normal_bounds (sys, ldy);
      
    nt->rerhist [nt->iters] = error;
    nt->merhist [nt->iters] = merit;

#if MPI
    if (dom->rank == 0)
#endif
    if (dom->verbose) printf (fmt, LINSYS_Iters (sys), LINSYS_Resnorm (sys), nt->iters, error, merit);

  } while (++ nt->iters < nt->maxiter && (error > nt->epsilon || merit > nt->meritval));

  /* FIXME: remember about coordinate change at the end of *_VARIATIONAL case iterations */

  LINSYS_Destroy (sys);
  free (nonmonvalues);
}

/* write labeled satate values */
void NEWTON_Write_State (NEWTON *nt, PBF *bf)
{
}

/* return variant string */
char* NEWTON_Variant (NEWTON *nt)
{
  switch (nt->variant)
  {
    case NONSMOOTH_HSW: return "NONSMOOTH_HSW";
    case NONSMOOTH_HYBRID: return "NONMOOTH_HYBRID";
    case FIXED_POINT: return "FIXED_POINT";
    case NONSMOOTH_VARIATIONAL: return "NONSMOOTH_VARIATIONAL";
    case SMOOTHED_VARIATIONAL: return "SMOOTHED_VARIATIONAL";
  }

  return NULL;
}

/* destroy solver */
void NEWTON_Destroy (NEWTON *nt)
{
  free (nt->rerhist);
  free (nt->merhist);
  free (nt);
}
