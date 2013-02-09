/*
 * cmp.h
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

#ifndef __cmp__
#define __cmp__

typedef enum
{
  CMP_OFF,
  CMP_FASTLZ
} CMP_ALG;

/* compress doubles and integers into an array of integers */
int* compress (CMP_ALG alg, double *d, int doubles, int *i, int ints, int *size);

/* decompress doubles and integers from an array of integers */
void decompress (int *input, int size, double **d, int *doubles, int **i, int *ints);

#endif
