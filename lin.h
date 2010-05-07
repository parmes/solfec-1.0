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

enum linvar
{
  /* linearization variants */
  NONSMOOTH_HSW          = 0x0001,
  NONSMOOTH_HYBRID       = 0x0002,
  FIXED_POINT            = 0x0004,
  NONSMOOTH_VARIATIONAL  = 0x0008,
  SMOOTHED_VARIATIONAL   = 0x0010,
  /* additional flags */
  BOUNDARY_ONLY          = 0x0100
};

typedef enum linvar LINVAR;

/* extract linearization variant from linvar */
#define LINEARIZATION_VARIANT(var) ((var) & 0x00FF)

/* create linear system resulting from linearization of constraints */
LINSYS* LINSYS_Create (LINVAR variant, LOCDYN *ldy);

/* set fixed point approach normal stress update error tolerance */
void LINSYS_Fixed_Point_Tol (LINSYS *sys, double tol);

/* update linear system at current reactions R */
void LINSYS_Update (LINSYS *sys);

/* solve for reaction increments DR */
void LINSYS_Solve (LINSYS *sys, double abstol, int maxiter);

/* compute merit function at (R + alpha * DR) */
double LINSYS_Merit (LINSYS *sys, double alpha);

/* advance solution R = R + alpha * DR; return |DR|/|R| */ 
double LINSYS_Advance (LINSYS *sys, double alpha);

/* solve A x = b, where b = A [1, 1, ..., 1]' and return |x - [1, 1, ..., 1]| / |[1, 1, ..., 1]| */
double LINSYS_Test (LINSYS *sys, double abstol, int maxiter);

/* most recent iterations count */
int LINSYS_Iters (LINSYS *sys);

/* most recent relative residual norm */
double LINSYS_Resnorm (LINSYS *sys);

#if MPI
/* update computed external reactions (e.g. non-gluing) */
void LINSYS_Update_External_Reactions (LINSYS *sys);
#endif

/* destroy linear system */
void LINSYS_Destroy (LINSYS *sys);

#endif
