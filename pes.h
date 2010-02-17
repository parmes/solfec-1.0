/*
 * pes.h
 * Copyright (C) 2007-2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * penalty constraints solver
 */

#include "ldy.h"

#ifndef __pes__
#define __pes__

typedef struct penalty PENALTY;

struct penalty
{
  short implicit;
};

/* create penalty solver */
PENALTY* PENALTY_Create (short implicit);

/* explcit constraint solver */
void PENALTY_Solve (PENALTY *ps, LOCDYN *ldy);

/* write labeled satate values */
void PENALTY_Write_State (PENALTY *ps, PBF *bf);

/* destroy penalty solver */
void PENALTY_Destroy (PENALTY *ps);

/* spring and dashpot based explicit diagonal block contact solver */
int PENALTY_Spring_Dashpot_Contact (CON *con, short implicit, double step, double gap, double spring, double dashpot,
                             double friction, double cohesion, double *W, double *B, double *V, double *U, double *R);
#endif
