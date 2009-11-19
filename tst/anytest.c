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
  printf ("UINT_MAX is %u\n", UINT_MAX);

  return 0;
}
