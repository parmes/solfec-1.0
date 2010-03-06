/*
 * pbs.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * per-body constraints solver
 */

#include "pbs.h"

/* create solver */
PER_BODY* PER_BODY_Create (PBMETHOD method, double epsilon, int maxiter)
{
  return NULL;
}

/* run solver */
void PER_BODY_Solve (PER_BODY *pb, LOCDYN *ldy)
{
}

/* write labeled satate values */
void PER_BODY_Write_State (PER_BODY *pb, PBF *bf)
{
}

/* destroy solver */
void PER_BODY_Destroy (PER_BODY *pb)
{
}
