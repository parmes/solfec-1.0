/*
 * tts.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * test solver
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

#ifndef __tts__
#define __tts__

typedef struct test TEST;

struct test
{
  /* input */

  double meritval; /* merit function value sufficient for termination */

  int maxiter; /* iterations bound */

  int linmaxiter; /* linear solver iterations bound */

  /* output */

  double *merhist; /* merit function history */

  int iters; /* iterations count */
};

/* create solver */
TEST* TEST_Create (double meritval, int maxiter);

/* run solver */
void TEST_Solve (TEST *ts, LOCDYN *ldy);

/* write labeled state values */
void TEST_Write_State (TEST *ts, PBF *bf);

/* destroy solver */
void TEST_Destroy (TEST *ts);

#endif
