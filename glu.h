/*
 * glu.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * linear gluing solver
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

#ifndef __glu__
#define __glu__

typedef struct glue GLUE;

/* create gluing solver */
GLUE* GLUE_Create (LOCDYN *ldy);

/* compute gluing reactions; return their relative change */
double GLUE_Solve (GLUE *glu, double abstol, int maxiter);

#if MPI
/* update external gluing reactions */
void GLUE_Update_External_Reactions (GLUE *glu);
#endif

/* destroy gluing solver */
void GLUE_Destroy (GLUE *glu);

#endif