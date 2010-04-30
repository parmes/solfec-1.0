/*
 * gls.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * gluing nonlinear constraint solver
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

#ifndef __gls__
#define __gls__

typedef struct gluing GLUING;

struct gluing
{
  double epsilon; /* relative accuracy of velocity projection */

  int maxiter; /* iterations bound of velocity projection */
};

/* create solver */
GLUING* GLUING_Create (double epsilon, int maxiter);

/* run solver */
void GLUING_Solve (GLUING *gl, LOCDYN *ldy);

/* write labeled satate values */
void GLUING_Write_State (GLUING *gl, PBF *bf);

/* destroy solver */
void GLUING_Destroy (GLUING *gl);

#endif
