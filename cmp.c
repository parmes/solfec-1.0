/*
 * cmp.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * compression / decompression
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
#include <string.h>

#include "ext/fastlz.h"
#include "cmp.h"
#include "alg.h"
#include "err.h"

/* compress doubles and integers into an array of integers */
int* compress (CMP_ALG alg, double *d, int doubles, int *i, int ints, int *size)
{
  int length, outsize, remainder;
  char *input, *output;

  length = sizeof (double [doubles]) + sizeof (int [ints]);

  ERRMEM (input = malloc (length));
  memcpy (input , d, sizeof (double [doubles]));
  memcpy (input + sizeof (double [doubles]), i, sizeof (int [ints]));

  if (length >= 16)
  {
    outsize = (int) (1.1 * (double) length);
    outsize = MAX (66, outsize);
    ERRMEM (output = malloc (outsize));

    outsize = fastlz_compress (input, length, output);

    free (input);
  }
  else
  {
    outsize = length;
    output = input;
  }

  remainder = outsize % sizeof (int);

  outsize += sizeof (int [4]) + (remainder ? sizeof (int) - remainder : 0);

  ERRMEM (output = realloc (output, outsize));

  *size = outsize / sizeof (int);

  ((int*) output) [(*size) - 4] = alg;
  ((int*) output) [(*size) - 3] = doubles;
  ((int*) output) [(*size) - 2] = ints;
  ((int*) output) [(*size) - 1] = remainder;

  return (int*) output;
}

/* decompress doubles and integers from an array of integers */
void decompress (int *input, int size, double **d, int *doubles, int **i, int *ints)
{
  int length, outsize, remainder;
  char *output;
  CMP_ALG alg;

  alg       = input [size - 4];
  *doubles  = input [size - 3];
  *ints     = input [size - 2];
  remainder = input [size - 1];

  if (*doubles) { ERRMEM (*d = malloc (sizeof (double [*doubles]))); }
  else *d = NULL;

  if (*ints) { ERRMEM (*i = malloc (sizeof (int [*ints]))); }
  else *i = NULL;

  outsize = sizeof (int) * size;

  length = outsize - sizeof (int [4]) + (remainder ? sizeof (int) - remainder : 0);

  if (length >= 16)
  {
    outsize = sizeof (double [*doubles]) + sizeof (int [*ints]);
    ERRMEM (output = malloc (outsize));

    fastlz_decompress (input, length, output, outsize);
  }
  else output = (char*) input;

  if (*d) memcpy (*d, output, sizeof (double [*doubles]));
  if (*i) memcpy (*i, output + sizeof (double [*doubles]), sizeof (int [*ints]));

  if (output != (char*) input) free (output);
}
