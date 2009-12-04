/*
 * tng.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * constraint root finding tangent operator routines
 */

#include "ldy.h"
#include "mtx.h"

#ifndef __tng__
#define __tng__

enum tangent_approach /* tangent linearisation approaches */
{
  SEMISMOOTH_STRICT,
  SEMISMOOTH_HYBRID,
  VARIATIONAL_NONSMOOTH,
  VARIATIONAL_SMOOTHED,
};

typedef enum tangent_approach TANGENT_APPROACH;

enum variational_smoothed_params /* 'params' indices for VARIATIONAL_SMOOTHED approach */
{
  DELTA = 0,
  EPSILON
};

#define TANGENT_PARAMS_SIZE 2  /* maximal size of 'params' argument below */

struct tangent_operator
{
  MX *tang; /* CSC tangent operator storage */

  MAP *idtoloc; /* map constraint identifiers to local indices */
};

typedef struct tangent_operator TANGENT_OPERATOR;

/* assemble a tangent operator for a specific linearisation approach */
TANGENT_OPERATOR* TANGENT_Assemble (LOCDYN *ldy, TANGENT_APPROACH approach, double *params);

/* update a tangent operator for a specific linearisation approach */
void TANGENT_Update (LOCDYN *ldy, TANGENT_APPROACH approach, double *params, TANGENT_OPERATOR *op);

/* compute a merit function for a specific linearisation approach */
double TANGENT_Merit (LOCDYN *ldy, TANGENT_APPROACH approach, double *params);

/* release tangent operator memory */
void TANGENT_Destroy (TANGENT_OPERATOR *op);

#endif
