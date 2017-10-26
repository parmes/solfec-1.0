/*
 * bgs.h
 * Copyright (C) 2007, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * block gauss seidel constraints solver
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

#include "ldy.h"
#include "dbs.h"

#ifndef __bgs__
#define __bgs__

typedef struct gs GAUSS_SEIDEL;
typedef void (*GAUSS_SEIDEL_Callback) (void*);

enum gserror
{ 
  GS_OK,
  GS_DIVERGED,
  GS_DIAGONAL_DIVERGED,
  GS_DIAGONAL_FAILED
};

enum gsfail
{ 
  GS_FAILURE_CONTINUE,
  GS_FAILURE_EXIT,
  GS_FAILURE_CALLBACK
};

enum gsonoff
{
  GS_OFF = 0,
  GS_ON
};

enum gsvariant
{
  GS_FULL,
  GS_MIDDLE_JACOBI,
  GS_BOUNDARY_JACOBI
};

typedef enum gserror GSERROR;
typedef enum gsfail GSFAIL;
typedef enum gsonoff GSONOFF;
typedef enum gsvariant GSVARIANT;

struct gs
{
  double epsilon; /* relative accuracy */

  int maxiter; /* iterations bound */

  double meritval; /* merit function value sufficient for termination */

  GSFAIL failure; /* action on failure */

  void *data; /* failure callback data */

  GAUSS_SEIDEL_Callback callback; /* failure callback function */

  double diagepsilon; /* diagonal block solver relative accuracy */

  int diagmaxiter; /* diagonal block solver iterations bound */

  DIAS diagsolver; /* diagonal block problem solver type */

  GSERROR error; /* error code */

  int iters; /* most recent number of iterations */

  int *itershist; /* iterations history of all solver calls */

  int itershistcount; /* < 0: do not collect iterations history; >= 0 the size of the iterhistory */

  int itershistsize; /* iterations history buffer size */

#if MPI
  int colors, bot, mid, top, inn; /* processor colors, bottom, middle, top and inner set sizes */
#endif

  double *rerhist; /* relative error history */

  double *merhist; /* merit function history */

  GSONOFF reverse; /* iterate forward an backward alternately ? */

  GSVARIANT variant; /* parallel algorithm variant (ignored in serial mode) */

  int innerloops; /* number of inner GS loops per one global parallel step (ignored in serial mode) */

  short verbose; /* local verbosity flag */

  short nomerit; /* merit function evaluation flag */
};

/* create solver */
GAUSS_SEIDEL* GAUSS_SEIDEL_Create (double epsilon, int maxiter, double meritval, GSFAIL failure,
                                   double diagepsilon, int diagmaxiter, DIAS diagsolver,
				   void *data, GAUSS_SEIDEL_Callback callback);

/* run solver */
void GAUSS_SEIDEL_Solve (GAUSS_SEIDEL *gs, LOCDYN *ldy);

/* return faulure string */
char* GAUSS_SEIDEL_Failure (GAUSS_SEIDEL *gs);

/* return diagonal solver string */
char* GAUSS_SEIDEL_Diagsolver (GAUSS_SEIDEL *gs);

/* return error string */
char* GAUSS_SEIDEL_Error (GAUSS_SEIDEL *gs);

/* return reverse flag string */
char* GAUSS_SEIDEL_Reverse (GAUSS_SEIDEL *gs);

/* return variant string */
char* GAUSS_SEIDEL_Variant (GAUSS_SEIDEL *gs);

/* write labeled satate values */
void GAUSS_SEIDEL_Write_State (GAUSS_SEIDEL *gs, PBF *bf);

/* free solver */
void GAUSS_SEIDEL_Destroy (GAUSS_SEIDEL *gs);

#endif
