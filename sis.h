/*
 * sis.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Siconos 3D contact solvers interface
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

#ifndef __sis__
#define __sis__

typedef struct siconos SICONOS;

struct siconos
{
  int maxiter;
};

/* create solver */
SICONOS* SICONOS_Create (double meritval, int maxiter);

/* run solver */
void SICONOS_Solve (SICONOS *si, LOCDYN *ldy);

/* write labeled state values */
void SICONOS_Write_State (SICONOS *si, PBF *bf);

/* destroy solver */
void SICONOS_Destroy (SICONOS *si);

#endif
