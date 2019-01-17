/*
 * pes.h
 * Copyright (C) 2007-2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * penalty constraints solver
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

#ifndef __pes__
#define __pes__

#if __cplusplus
extern "C" { /* C */
#endif

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
int PENALTY_Spring_Dashpot_Contact (CON *con, short implicit, double step, double gap, double spring, double dashpot, double hpow,
                             double friction, double cohesion, double *W, double *B, double *V, double *U, double *R);
#if __cplusplus
} /* extern C */
#endif

#endif
