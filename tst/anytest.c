/*
 * glvtest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of graphical viewer
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

#include <complex.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

int main (int argc, char **argv)
{
  double nan = 1.0 / 0.0;

  printf ("Result of 'isnan' on %g is %d\n", nan, isnan (nan));
  printf ("Result of 'isfinite' on %g is %d\n", nan, isfinite (nan));
  printf ("Size of 'long' is %ld\n", sizeof (long));
  printf ("Size of 'long long' is %ld\n", sizeof (long long));
#if !defined (__MINGW32__)
  printf ("Size of 'off_t' is %ld\n", sizeof (off_t));
#endif
  printf ("UINT_MAX is %u\n", UINT_MAX);

  /* 3x3 Gauss elimination */
  double x [3] = {1, 2, 3}, T [9] = {5, 2, -1, 4, 7, 2, -8, 1, 0.5};

  printf ("inv ([%g, %g, %g; %g, %g, %g; %g, %g, %g])*([%g, %g, %g]') =",
          T[0], T[3], T[6], T[1], T[4], T[7], T[2], T[5], T[8], x[0], x[1], x[2]);

  T [3] /= T[0]; T [6] /= T[0]; x [0] /= T[0];
  T [4] -= T[3]*T[1]; T [7] -= T[6]*T[1]; x [1] -= x[0]*T[1];
  T [5] -= T[3]*T[2]; T [8] -= T[6]*T[2]; x [2] -= x[0]*T[2];
  T [7] /= T [4]; x [1] /= T[4];
  T [8] -= T[7]*T[5]; x [2] -= x[1]*T[5];
  x [2] /= T [8];
  x [1] = x[1] - T[7]*x[2];
  x [0] = x[0] - T[3]*x[1] - T[6]*x[2];

  printf ("[%g, %g, %g]'\n", x[0], x[1], x[2]);

  return 0;
}
