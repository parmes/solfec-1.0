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

double GEOMETRIC_EPSILON = 1E-6;

/* lexicographical comparison */
int POINTS_COMPARE (double *a, double *b)
{
  if (LT (a[0], b[0])) return -1;
  else if (EQ (a[0], b[0]))
  {
    if (LT (a[1], b[1])) return -1;
    else if (EQ (a[1], b[1]))
    {
      if (LT (a[2], b[2])) return -1;
      else if (EQ (a[2], b[2])) return 0;
    }
  }

  return 1;
}
