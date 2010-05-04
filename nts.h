/*
 * nts.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Newton constraints solver
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

#include "lin.h"

#ifndef __nts__
#define __nts__

typedef struct newton NEWTON;

struct newton
{
  short variant; /* linearization variant */

  double epsilon; /* relative accuracy sufficient for termination */

  int maxiter; /* iterations bound */

  double meritval; /* merit function value sufficient for termination */

  int nonmonlength; /* nonmonotone merit function buffer length */

  int linmaxiter; /* linear solver iterations bound */

  int iters; /* most recent number of iterations */

  double *rerhist; /* relative error history */

  double *merhist; /* merit function history */
};

/* create solver */
NEWTON* NEWTON_Create (LINVAR variant, double epsilon, int maxiter, double meritval);

/* run solver */
void NEWTON_Solve (NEWTON *nt, LOCDYN *ldy);

/* write labeled satate values */
void NEWTON_Write_State (NEWTON *nt, PBF *bf);

/* return variant string */
char* NEWTON_Variant (NEWTON *nt);

/* destroy solver */
void NEWTON_Destroy (NEWTON *nt);
#endif
