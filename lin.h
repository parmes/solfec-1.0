/*
 * lin.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * constraints linearization routines
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

#ifndef __lin__
#define __lin__

typedef struct linsys LINSYS;

enum linvar /* linearization variant */
{
  NONSMOOTH_HSW,
  NONSMOOTH_HYBRID,
  FIXED_POINT,
  NONSMOOTH_VARIATIONAL,
  SMOOTHED_VARIATIONAL
};

typedef enum linvar LINVAR;

enum linopt /* linear system options */
{
  LOCAL_SYSTEM = 0x01, /* in parallel create per-processor local systems */
  SYMMETRIZE = 0x02 /* symmetrize equations if possible */
};

typedef enum linopt LINOPT;

/* create linear system resulting from constraints linearization */
LINSYS* LINSYS_Create (LINVAR variant, LINOPT options, LOCDYN *ldy);

/* update linear system at current R and U */
void LINSYS_Update (LINSYS *sys);

/* solve linear system for reaction increments DR */
void LINSYS_Solve (LINSYS *sys, double accuracy, int maxiter);

/* compute merit function at (R + alpha * DR) */
double LINSYS_Merit (LINSYS *sys, double alpha);

/* advance solution R = R + alpha * DR; return |DR|/|R| */ 
double LINSYS_Advance (LINSYS *sys, double alpha);

/* test solution of W * x = b, where b = W * [1, 1, ..., 1];
 * return |x - [1, 1, ..., 1]| / |[1, 1, ..., 1]| */
double LINSYS_Test (LINSYS *sys, double accuracy, int maxiter);

/* returns 1 if a global system; 0 otherwise */
int LINSYS_Global (LINSYS *sys);

/* most recent iterations count */
int LINSYS_Iters (LINSYS *sys);

/* most recent residual norm */
double LINSYS_Resnorm (LINSYS *sys);

/* destroy linear system */
void LINSYS_Destroy (LINSYS *sys);

/* constraint satisfaction merit function;
 * (assumes that both dia->R and dia->U are valid) */
double MERIT_Function (LOCDYN *ldy);
#endif
