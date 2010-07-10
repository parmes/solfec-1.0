/*
 * bss.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * body space solver
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

#include <stdlib.h>
#include <float.h>

#include "bss.h"
#include "alg.h"
#include "err.h"
#include "mrf.h"

/* create solver */
BSS* BSS_Create (int maxiter, double meritval)
{
  return NULL;
}

/* run solver */
void BSS_Solve (BSS *bs, LOCDYN *ldy)
{
}

/* write labeled satate values */
void BSS_Write_State (BSS *bs, PBF *bf)
{
}

/* destroy solver */
void BSS_Destroy (BSS *bs)
{
}
