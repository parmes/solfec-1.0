/*
 * pck.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * pack and unpack data into doubles and ints
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

#ifndef __pck__
#define __pck__

/* pack vector of doubles */
void pack_doubles (int *dsize, double **d, int *doubles, double *input, int count);

/* pack vector of ints */
void pack_ints (int *isize, int **i, int *ints, int *input, int count);

/* pack string into ints */
void pack_string (int *isize, int **i, int *ints, char *input);

/* pack single double */
void pack_double (int *dsize, double **d, int *doubles, double input);

/* pack single int */
void pack_int (int *isize, int **i, int *ints, int input);

/* unpack vector of doubles */
void unpack_doubles (int *dpos, double *d, int doubles, double *output, int count);

/* unpack vector of ints */
void unpack_ints (int *ipos, int *i, int ints, int *output, int count);

/* unpack string */
char* unpack_string (int *ipos, int *i, int ints);

/* unpack single double */
double unpack_double (int *dpos, double *d, int doubles);

/* unpack single int */
int unpack_int (int *ipos, int *i, int ints);

#endif
