/*
 * gls.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Gajulapalli-Lasdon matrix scaling
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
#include <math.h>
#include "gls.h"

/* n: matrix dimension (input)
 * q: matrix coefficient (intput/output); scaled on exit
 * p: pointers to columns (input)
 * i: row indices (input)
 * x: row scaling (output)
 * y: column scaling (output)
 * return: 1 when done or 0 when out of memory;
 * try finding x and y such that x[i] |a[i, j]| y [j] is in [0.5, 1] */
int gls_csc (int n, double *q, int *p, int *i, double *x, double *y)
{
  double *mem, *mat, *xin, *yin, *a, *b, *r, *c, *e, t;
  int *col, *row, j, *k, m;

  mat = q;
  col = p;
  row = i;
  xin = x;
  yin = y;
  if (!(mem = calloc (4 * n, sizeof (double)))) return 0;
  a = mem;
  b = a + n;
  r = b + n;
  c = r + n;

  /* initialization */

  for (j = 0; j < n; j ++, p ++, b ++, c ++)
  {
    for (k = i + (*(p+1) - (*p)); i < k; i ++, q ++)
    {
      if ((*q) != 0.0)
      {
	t = - log2 (fabs (*q)) - 0.5;
	r [*i] += 1.0;
	(*c) += 1.0;
	a [*i] += t;
	(*b) += t;
      }
    }
  }

  a = mem;
  b = a + n;
  r = b + n;
  c = r + n;
  e = b;

  for (;a < e; a ++, b ++, r ++, c ++)
  {
    (*a) /= (*r);
    (*b) /= (*c);
  }

  /* computation */

  for (e = y + n; y < e; y ++) *y = 0.0;

  for (m = 0; m < 3; m ++)
  {
    a = mem;
    b = a + n;
    r = b + n;
    c = r + n;

    for (x = xin, e = x + n; x < e; x ++, a ++) (*x) = (*a);
   
    for (j = 0, p = col, i = row, x = xin, y = yin, q = mat; j < n; j ++, p ++, y ++)
    {
      for (k = i + (*(p+1) - (*p)); i < k; i ++, q ++)
      {
	if ((*q) != 0.0)
	{
	  x [*i] -= (*y) / r [*i];
	}
      }
    }

    for (x = xin, y = yin; x < e; x ++, y ++, b ++)
    {
      (*x) = (int) (*x);
      (*y) = (*b);
    }

    for (j = 0, p = col, i = row, x = xin, y = yin, q = mat; j < n; j ++, p ++, y ++, c ++)
    {
      for (k = i + (*(p+1) - (*p)); i < k; i ++, q ++)
      {
	if ((*q) != 0.0)
	{
	  (*y)  -= x[*i] / (*c);
	}
      }
    }

    for (y = yin, e = y + n; y < e; y ++) (*y) = (int) (*y);
  }

  for (x = xin, y = yin, e = x + n; x < e; x ++, y ++)
  {
    (*x) = exp2 (*x);
    (*y) = exp2 (*y);
  }

  for (j = 0, q = mat, p = col, i = row, x = xin, y = yin; j < n; j ++, p ++, y ++)
  {
    for (k = i + (*(p+1) - (*p)); i < k; i ++, q ++)
    {
      (*q) *= x[*i] * (*y);
    }
  }

  free (mem);

  return 1;
}