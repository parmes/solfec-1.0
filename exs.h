/*
 * exs.h
 * Copyright (C) 2007, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * explicit constraints solver
 */

#include "ldy.h"

#ifndef __exs__
#define __exs__

/* spring and dashpot based explicit diagonal block contact solver */
int EXPLICIT_Spring_Dashpot_Contact (CON *con, double step, double gap, double spring, double dashpot,
              double friction, double cohesion, double *W, double *B, double *V, double *U, double *R);

/* explcit constraint solver */
void EXPLICIT_Solve (LOCDYN *ldy);

#endif
