/*
 * pbs.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * per-body constraints solver
 */

#include "ldy.h"

#ifndef __pbs__
#define __pbs__

typedef struct pb PER_BODY;

enum pbmethod
{
  PB_GAUSS_SEIDEL = 0,
  PB_NEWTON
};

typedef enum pbmethod PBMETHOD;

struct pb
{
  PBMETHOD method; /* solution method */

  double epsilon; /* relative accuracy */

  int maxiter; /* iterations bound */
};

/* create solver */
PER_BODY* PER_BODY_Create (PBMETHOD method, double epsilon, int maxiter);

/* run solver */
void PER_BODY_Solve (PER_BODY *pb, LOCDYN *ldy);

/* write labeled satate values */
void PER_BODY_Write_State (PER_BODY *pb, PBF *bf);

/* destroy solver */
void PER_BODY_Destroy (PER_BODY *pb);

#endif
