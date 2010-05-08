/*
 * hbs.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * hybrid constraint solver
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

#include "hbs.h"
#include "alg.h"
#include "bla.h"
#include "err.h"
#include "bgs.h"
#include "nts.h"
#include "dom.h"
#include "mrf.h"

/* create solver */
HYBRID* HYBRID_Create (double epsilon, int maxiter)
{
  HYBRID *hs;

  ERRMEM (hs = malloc (sizeof (HYBRID)));
  hs->epsilon = epsilon;
  hs->maxiter = maxiter;

  return hs;
}

/* run solver */
void HYBRID_Solve (HYBRID *hs, LOCDYN *ldy)
{
}

/* write labeled satate values */
void HYBRID_Write_State (HYBRID *hs, PBF *bf)
{
}

/* destroy solver */
void HYBRID_Destroy (HYBRID *hs)
{
  free (hs);
}
