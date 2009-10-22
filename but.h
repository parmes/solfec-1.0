/*
 * but.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * body related utilities
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

#include "alg.h"

#ifndef __but__
#define __but__

/* convert Euler tensor to the inertia tensor */
inline static void euler2inertia (double *euler, double *inertia)
{
  double trace;

  trace = TRACE (euler);
  IDENTITY (inertia);
  SCALE9 (inertia, trace);
  NNSUB (inertia, euler, inertia); /* inertia = tr(euler)*one - euler */
}

/* Lame's lambda coefficient */
inline static double lambda (double young, double poisson)
{
  return young*poisson / ((1.0 + poisson)*(1.0 - 2.0*poisson));
}

/* Lame's mi coefficient */
inline static double mi (double young, double poisson)
{
  return young / (2.0*(1.0 + poisson));
}

#endif
