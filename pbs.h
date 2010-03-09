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

struct pb
{
  double epsilon; /* relative accuracy */

  int maxiter; /* iterations bound */
};

/* create solver */
PER_BODY* PER_BODY_Create (double epsilon, int maxiter);

/* run solver */
void PER_BODY_Solve (PER_BODY *pb, LOCDYN *ldy);

/* write labeled satate values */
void PER_BODY_Write_State (PER_BODY *pb, PBF *bf);

/* destroy solver */
void PER_BODY_Destroy (PER_BODY *pb);

#endif
