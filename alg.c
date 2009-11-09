/*
 * alg.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * basic operations on scalars, vectors, matrices ...
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

#define DELTA 1.0E-4
#define SAMPLE 8

double GEOMETRIC_EPSILON = DELTA;
static double DENOM = 1.0;

void GEOMETRIC_EPSILON_ADAPT (double *p, int n)
{
  double *a, *b, u [3], e, sum, add;
  int i, j, k;

  add = 0.0;
  k = MIN (SAMPLE, n);
  for (i = j =  0; i < k; i ++)
  {
    a = p + (rand () % n)*3;
    b = p + (rand () % n)*3;
    SUB (a, b, u);
    MAXABS (u, e); 
    add += e;
    if (e > 0.0) j ++;
  }

  if (j) 
  {
    sum = DENOM * (GEOMETRIC_EPSILON / DELTA);
    sum += add;
    DENOM += (double)j;
    GEOMETRIC_EPSILON = DELTA * (sum / DENOM);
  }
}
