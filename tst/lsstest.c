/*
 * tst.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * -----------------------------------------------------------
 * linear system solver test
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
#include <ctype.h>
#include <stdio.h>
#include <math.h>
#include "ext/csparse.h"
#include "mtx.h"
#include "lss.h"
#include "err.h"

static void output (void *lss)
{
  int levels = LSS_Get (lss, LSS_LEVELS),
      iterations = LSS_Get (lss, LSS_ITERATIONS);

  double *ld = LSS_Getv (lss, LSS_LEVEL_DIMENSIONS),
	 *ln = LSS_Getv (lss, LSS_LEVEL_NONZEROS),
	 *ae = LSS_Getv (lss, LSS_ABSOLUTE_ERROR);

  FILE *aef = fopen ("lsstest-absolute-error.dat", "w");

  printf ("operator complexity = %f, grid complexity = %f\n",
    LSS_Get (lss, LSS_OPERATOR_COMPLEXITY), LSS_Get (lss, LSS_GRID_COMPLEXITY));

  for (int i = 0; i < levels; i ++)
  {
    printf ("level [%d]: dim = %d, nnz = %d\n", i, (int)ld[i], (int)ln[i]);
  }

  for (int i = 0; i < iterations; i ++)
  {
    fprintf (aef, "%e\n", ae [i]);
  }

  fclose (aef);
}

static void skipheader (FILE *f)
{
  int c;

  do
  {
    while (isspace (c = fgetc (f)));

    if (c != '%')
    {
      ungetc (c, f);
      return;
    }

    do c = fgetc (f); while (c != '\n' && c != EOF);

  } while (c != EOF);
}

static void sparsetest (char *path)
{
  int i, j, k, n, m, nnz;
  double a, *b, *x;
  void *lss;
  FILE *f;
  MX *AA, *A;

  ASSERT (f = fopen (path,"r"), ERR_FILE_OPEN);

  skipheader (f);

  printf ("Reading file ...\n");

  ASSERT (fscanf (f, "%d", &n) != EOF, ERR_FILE_FORMAT);
  ASSERT (fscanf (f, "%d", &m) != EOF, ERR_FILE_FORMAT);
  ASSERT (fscanf (f, "%d", &nnz) != EOF, ERR_FILE_FORMAT);
  ASSERT (n == m, ERR_FILE_FORMAT);

  ERRMEM (x = calloc (n, sizeof (double)));
  ERRMEM (b = calloc (n, sizeof (double)));
  AA = cs_spalloc (0, 0, 1, 1, 1) ;

  for (m = k = 0; k < nnz; k ++)
  {
    ASSERT (fscanf (f, "%d", &i) != EOF, ERR_FILE_FORMAT);
    ASSERT (fscanf (f, "%d", &j) != EOF, ERR_FILE_FORMAT);
    ASSERT (fscanf (f, "%lf", &a) != EOF, ERR_FILE_FORMAT);

    i --;
    j --;

    ERRMEM (cs_entry (AA, i, j, a));

    b [i] += a; /* all ones solution */

    if (i == j) m ++;
  }

  fclose (f);

  ASSERT (n == m, ERR_FILE_FORMAT);

  printf ("Setting up ...\n");

  A = cs_compress (AA);

  ASSERT (A->n == A->m && A->n == n, ERR_FILE_FORMAT);

  lss = LSS_Create (A->n, A->x, A->p, A->i);

  LSS_Set (lss, LSS_PRECONDITIONER, 3);
  LSS_Set (lss, LSS_SMOOTHING_STEPS, 3);
  LSS_Set (lss, LSS_DECIMATION, 8);
  LSS_Set (lss, LSS_RESTART, 20);
  LSS_Set (lss, LSS_CUTOFF, 16);
  LSS_Set (lss, LSS_ITERATIONS_BOUND, 200);
  LSS_Set (lss, LSS_RELATIVE_ACCURACY, 1E-6);
  LSS_Set (lss, LSS_ABSOLUTE_ACCURACY, 1E-6);

  printf ("Solving ...\n");

  if (LSS_Solve (lss, A->x, x, b) != LSSERR_NONE)
  {
    printf ("LSS ERROR: %s\n", LSS_Errmsg (lss));
  }

  for (a = 0.0, i = 0; i < n; i ++)
  {
    a += (x[i] - 1.0) * (x[i] - 1.0);
  }
  a = sqrt (a);

  printf ("Solution norm |x - 1| = %g\n", a);

  output (lss);

  LSS_Destroy (lss);
  MX_Destroy (AA);
  MX_Destroy (A);
  free (b);
  free (x);

  printf ("Done.\n");
}

int main (int argc, char **argv)
{
  if (argc != 2)
  {
    printf ("SYNOPSIS: lsstest path/to/matrix-market/square/matrix.mtx\n");
  }
  else
  {
    system ("rm -f *.dat *.eps");

    sparsetest (argv [1]);

    system ("gnuplot inp/lsstest.plt");
  }

  return 0;
}
