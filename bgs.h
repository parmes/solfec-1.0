/*
 * bgs.h
 * Copyright (C) 2007, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * block gauss seidel constraints solver
 */

#include "ldy.h"

#if MPI
#include <zoltan.h>
#endif

#ifndef __bgs__
#define __bgs__

typedef struct gs GAUSS_SEIDEL;
typedef void (*GAUSS_SEIDEL_Callback) (void*);

enum gserror
{ 
  GS_OK,
  GS_DIVERGED,
  GS_DIAGONAL_DIVERGED,
  GS_DIAGONAL_FAILED,
};

enum gsdias
{
  GS_PROJECTED_GRADIENT,
  GS_DE_SAXE_AND_FENG,
  GS_SEMISMOOTH_NEWTON
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
typedef enum gsdias GSDIAS;
typedef enum gsfail GSFAIL;
typedef enum gsonoff GSONOFF;
typedef enum gsvariant GSVARIANT;

struct gs
{
  double epsilon; /* relative accuracy */

  int maxiter; /* iterations bound */

  GSFAIL failure; /* action on failure */

  void *data; /* failure callback data */

  GAUSS_SEIDEL_Callback callback; /* failure callback function */

  double diagepsilon; /* diagonal block solver relative accuracy */

  int diagmaxiter; /* diagonal block solver iterations bound */

  GSDIAS diagsolver; /* diagonal block problem solver type */

  GSERROR error; /* error code */

  int iters; /* most recent number of iterations */

  GSONOFF history; /* error history recording flag */

  double *rerhist; /* relative error history */

  GSONOFF reverse; /* iterate forward an backward alternately ? */

  GSVARIANT variant; /* parallel algorithm variant (ignored in serial mode) */
};

/* create solver */
GAUSS_SEIDEL* GAUSS_SEIDEL_Create (double epsilon, int maxiter, GSFAIL failure,
                                   double diagepsilon, int diagmaxiter, GSDIAS diagsolver,
				   void *data, GAUSS_SEIDEL_Callback callback);

/* run solver */
void GAUSS_SEIDEL_Solve (GAUSS_SEIDEL *gs, LOCDYN *ldy);

/* return faulure string */
char* GAUSS_SEIDEL_Failure (GAUSS_SEIDEL *gs);

/* return diagonal solver string */
char* GAUSS_SEIDEL_Diagsolver (GAUSS_SEIDEL *gs);

/* return error string */
char* GAUSS_SEIDEL_Error (GAUSS_SEIDEL *gs);

/* return history flag string */
char* GAUSS_SEIDEL_History (GAUSS_SEIDEL *gs);

/* return reverse flag string */
char* GAUSS_SEIDEL_Reverse (GAUSS_SEIDEL *gs);

/* return variant string */
char* GAUSS_SEIDEL_Variant (GAUSS_SEIDEL *gs);

/* free solver */
void GAUSS_SEIDEL_Destroy (GAUSS_SEIDEL *gs);

/* diagonal block solver */
int DIAGONAL_BLOCK_Solver (GSDIAS diagsolver, double diagepsilon, int diagmaxiter,
  short dynamic, double step, short kind, SURFACE_MATERIAL *mat, double gap,
  double *Z, double *base, DIAB *dia, double *B);
#endif
