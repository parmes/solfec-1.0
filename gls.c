/*
 * gls.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * gluing nonlinear constraint solver
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
#include "put.h"
#endif

#include "gls.h"
#include "glu.h"
#include "alg.h"
#include "bla.h"
#include "err.h"
#include "bgs.h"
#include "nts.h"
#include "dom.h"
#include "mrf.h"

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
static double reactions_update (int nonmonlength, LINSYS *sys, LOCDYN *ldy, double *nonmonvalues, int iter, double *merit)
{
  double reference, alpha;

  if (iter) /* skip line search during the first iteration */
  {
    (*merit) = LINSYS_Merit (sys, 0.0);

    reference = nonmonotone_merit_function (nonmonvalues, nonmonlength, (*merit), iter);

    alpha = line_search (sys, reference, merit);
  }
  else alpha = 1.0;

  return LINSYS_Advance (sys, alpha); /* R = R + alpha * DR, U = U(R) */
}

static void gluing_newton_solver (GLUING *gl, LOCDYN *ldy)
{
  double merit, error, *nonmonvalues;
  char fmt [512];
  LINSYS *sys;
  GLUE *glue;
  DOM *dom;
  int iters, maxiter, linmaxiter, nonmonlength;
  double epsilon, tol, meritval;

  dom = ldy->dom;


  sys = LINSYS_Create (NONSMOOTH_HSW|NON_GLUING, ldy);
  glue = GLUE_Create (ldy);
  merit = error = 1.0;
  iters = 0;
  maxiter = gl->maxiter;
  epsilon = gl->epsilon;
  linmaxiter = 1000;
  meritval = 1E-3 * gl->epsilon;
  nonmonlength = 1;

  ERRMEM (nonmonvalues = MEM_CALLOC (nonmonlength * sizeof (double)));

  if (dom->verbose)
  {
    sprintf (fmt, "GLUE: (LIN its/res: %%%dd/%%.2e) iteration: %%%dd  error:  %%.2e  merit: %%.2e\n",
      (int)log10 (linmaxiter) + 1, (int)log10 (maxiter) + 1);
  }

  do
  {
    LINSYS_Update (sys); /* assemble A, b */

    tol = 1E-3 * MIN (error, merit);
    tol = MAX (tol, 1E-15);
    LINSYS_Solve (sys, tol, linmaxiter);

    error = reactions_update (nonmonlength, sys, ldy, nonmonvalues, iters, &merit); /* R(i+1) */

    if (error < epsilon && merit < meritval)
    {
#if MPI
      DOM_Update_External_Reactions (dom, 0);
#endif
      error = merit = GLUE_Solve (glue, tol, linmaxiter);
#if MPI
      GLUE_Update_External_Reactions (glue); 
#endif
    }

#if MPI
    if (dom->rank == 0)
#endif
    if (dom->verbose) printf (fmt, LINSYS_Iters (sys), LINSYS_Resnorm (sys), iters, error, merit);

  } while (++ iters < maxiter && (error > epsilon || merit > meritval));

  LINSYS_Destroy (sys);
  GLUE_Destroy (glue);
  free (nonmonvalues);
}

/* create solver */
GLUING* GLUING_Create (double epsilon, int maxiter)
{
  GLUING *gl;

  ERRMEM (gl = malloc (sizeof (GLUING)));
  gl->epsilon = epsilon;
  gl->maxiter = maxiter;

  return gl;
}

/* run solver */
void GLUING_Solve (GLUING *gl, LOCDYN *ldy)
{
#if 0
  double error, merit, step;
  int verbose, diagiters;
  int div = 10, iters;
  short dynamic;
  char fmt [512];

  verbose = ldy->dom->verbose;

  if (verbose) sprintf (fmt, "GLUING: iteration: %%%dd  error:  %%.2e  merit:  %%.2e\n", (int)log10 (1000) + 1);

  int maxiter = 100;
  double epsilon = 1E-3,
	 meritval = 1.0E6;

  GLUE *glue = GLUE_Create (ldy);

  dynamic = ldy->dom->dynamic;
  step = ldy->dom->step;
  iters = 0;
  do
  {
    double errup = 0.0,
	   errlo = 0.0;
    OFFB *blk;
    DIAB *dia;
   
    for (dia = ldy->dia; dia; dia = dia->n)
    {
      double R0 [3],
	     B [3],
	     *R = dia->R;

      CON *con = dia->con;

      if (con->kind == GLUEPNT) continue;

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
	       *R = CON(blk->dia)->R;
	NVADDMUL (B, W, R, B);
      }
#endif
      
      COPY (R, R0); /* previous reaction */

      /* solve local diagonal block problem */
      diagiters = DIAGONAL_BLOCK_Solver (GS_PROJECTED_GRADIENT, 1E-6, 100,
	         dynamic, step, con->kind, con->mat.base, con->gap, con->Z, con->base, dia, B);

      /* accumulate relative
       * error components */
      SUB (R, R0, R0);
      errup += DOT (R0, R0);
      errlo += DOT (R, R);
    }

    /* merit function value */
    merit = MERIT_Function (ldy, 1);

    /* calculate relative error */
    error = sqrt (errup) / sqrt (MAX (errlo, 1.0));

#if MPI
    error = PUT_double_max (error);
#endif

    GLUE_Solve (glue, 1E-5, 100);

#if MPI
    DOM_Update_External_Reactions (ldy->dom, 0);
#endif

    if (iters % div == 0 && verbose) printf (fmt, iters, error, merit), div *= 2;
  }
  while (++ iters < maxiter && (error > epsilon || merit > meritval));

  if (verbose) printf (fmt, iters, error, merit);

   GLUE_Destroy (glue); 
#else
   gluing_newton_solver (gl, ldy);
#endif
}

/* write labeled satate values */
void GLUING_Write_State (GLUING *gl, PBF *bf)
{
}

/* destroy solver */
void GLUING_Destroy (GLUING *gl)
{
  free (gl);
}
