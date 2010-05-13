/*
 * dbs.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * diagonal block constraints solver
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

#ifndef __dbs__
#define __dbs__

enum dias
{
  DS_PROJECTED_GRADIENT,
  DS_DE_SAXE_AND_FENG,
  DS_SEMISMOOTH_NEWTON
};

typedef enum dias DIAS;

/* diagsolver: diagonal solver kind
 * diagepsilon: relative accuracy on termination
 * diagmaxiter: maximal iterations count
 * dynamic: simulation kind (dom->dynamic)
 * step: time step (dom->step)
 * kind: constraint kind (con->kind)
 * mat: surface material (kind == CONTACT)
 * gap: constraint gap
 * Z: auxiliary Z storage (con->Z)
 * base: constraint local base (con->base)
 * dia: diagonal block of local dynamic (con->dia)
 * B: local free velocity (B = dia->B + sum [dia->adj] (W_i R_i));
 * diagonal block solver */
int DIAGONAL_BLOCK_Solver (DIAS diagsolver, double diagepsilon, int diagmaxiter,
  short dynamic, double step, short kind, SURFACE_MATERIAL *mat, double gap,
  double *Z, double *base, DIAB *dia, double *B);

#endif
