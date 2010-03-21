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

#include "ldy.h"

#ifndef __nts__
#define __nts__

typedef struct newton NEWTON;

enum ntvariant
{
  NT_NONSMOOTH_HSW          =  0x01,
  NT_NONSMOOTH_HYBRID       =  0x02,
  NT_FIXED_POINT            =  0x04,
  NT_NONSMOOTH_VARIATIONAL  =  0x08,
  NT_SMOOTHED_VARIATIONAL   =  0x10,
  NT_SYMMETRIZE             =  0x20   /* symmetrize linear equations if possible */
};

typedef enum ntvariant NTVARIANT;

struct newton
{
  NTVARIANT variant; /* method variant */

  double epsilon; /* relative accuracy sufficient for termination */

  int maxiter; /* iterations bound */

  double meritval; /* merit function value sufficient for termination */

  int nonmonlength; /* nonmonotone merit function buffer length */

  int linmaxiter; /* linear solver iterations bound */
};

/* create solver */
NEWTON* NEWTON_Create (NTVARIANT variant, double epsilon, int maxiter, double meritval);

/* run solver */
void NEWTON_Solve (NEWTON *nt, LOCDYN *ldy);

/* write labeled satate values */
void NEWTON_Write_State (NEWTON *nt, PBF *bf);

/* return variant string */
char* NEWTON_Variant (NEWTON *nt);

/* destroy solver */
void NEWTON_Destroy (NEWTON *nt);

/* constraint satisfaction merit function;
 * (assumes that both dia->R and dia->U are valid) */
double MERIT_Function (LOCDYN *ldy);
#endif
