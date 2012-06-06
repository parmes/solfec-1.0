/*
 * cmptest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * compression test
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
#include "cmp.h"
#include "alg.h"
#include "err.h"

int main (int argc, char **argv)
{
  double *d;
  int doubles;

  int *i;
  int ints;

  int *data;
  int size;

  double *ud;
  int udl;

  int *ui;
  int uil;

  int n, length;

  if (argc > 1) length = atoi (argv [1]);
  
  doubles = ints = length = MAX (128, length);

  ERRMEM (d = malloc (sizeof (double [length])));
  ERRMEM (i = malloc (sizeof (int [length])));

  for (n = 0; n < length; n ++)
  {
    d [n] = (double) n;
    i [n] = n;
  }

  printf ("Compressing %d doubles and ints ... ", length);

  data = compress (CMP_FASTLZ, d, doubles, i, ints, &size);

  printf ("done.\n");

  printf ("Compression ratio: %g\n", (double) (sizeof (double [length])  + sizeof (int [ints])) / (double) sizeof (int [size]));

  printf ("Decompressing ... ");

  decompress (data, size, &ud, &udl, &ui, &uil);

  printf ("done.\n");

  printf ("Decompressed %d doubles and %d ints, ", udl, uil);
  if (udl != length || uil != length) printf ("ERROR\n");
  else
  {
    printf ("OK\nDecompressed values: ");

    for (n = 0; n < length; n ++)
    {
      if (ud [n] != d [n] || ui [n] != i [n]) break;
    }

    if (n == length) printf ("OK\n");
    else printf ("ERROR\n");
  }

  free (d);
  free (i);
  free (ud);
  free (ui);

  return 0;
}
