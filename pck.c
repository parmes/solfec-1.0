/*
 * pck.c
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

#include <string.h>
#include <stdlib.h>
#include "pck.h"
#include "err.h"

#define pack_many(size, array, stored, input, count, type)\
  do\
  {\
    type *slot, *iter;\
    int n = count;\
  \
    if ((*size) < (*stored) + count)\
    {\
      (*size) = 2 * ((*stored) + count);\
      ERRMEM ((*array) = realloc ((*array), (*size) * sizeof (type)));\
    }\
  \
    for (slot = &(*array)[*stored], iter = input, (*stored) += count;\
	 n > 0; n --, slot ++, iter ++) *slot = *iter;\
  }\
  while (0)

#define pack_one(size, array, stored, input, type)\
  do\
  {\
  \
    if ((*size) < (*stored) + 1)\
    {\
      (*size) = 2 * ((*stored) + 1);\
      ERRMEM ((*array) = realloc ((*array), (*size) * sizeof (type)));\
    }\
  \
    (*array)[*stored] = input;\
    (*stored) += 1;\
  }\
  while (0)

#define unpack_many(pos, array, stored, output, count, type)\
  do\
  {\
    type *slot, *iter;\
    int n = count;\
  \
    ASSERT (((*pos) + count) <= stored, ERR_PCK_UNPACK);\
  \
    for (iter = output, slot = &array[*pos], (*pos) += count;\
	 n > 0; n --, slot ++, iter ++) *iter = *slot;\
  }\
  while (0)

#define unpack_one(pos, array, stored, output, type)\
  do\
  {\
    ASSERT (((*pos) + 1) <= stored, ERR_PCK_UNPACK);\
  \
    output = array[*pos];\
    (*pos) += 1;\
  }\
  while (0)

/* pack vector of doubles */
void pack_doubles (int *dsize, double **d, int *doubles, double *input, int count)
{
  pack_many (dsize, d, doubles, input, count, double);
}

/* pack vector of ints */
void pack_ints (int *isize, int **i, int *ints, int *input, int count)
{
  pack_many (isize, i, ints, input, count, int);
}

/* pack string into ints */
void pack_string (int *isize, int **i, int *ints, char *input)
{
  if (input)
  {
    int len, n, *buf;

    len = strlen (input) + 1; /* + '\0' */

    ERRMEM (buf = malloc (len * sizeof (int)));

    for (n = 0; n < len; n ++) buf [n] = (int)input [n]; /* one-to-one for portability sake (heterogenous platforms) */

    pack_one (isize, i, ints, len, int); /* length of the integer packed string */
    pack_many (isize, i, ints, buf, len, int); /* the integer packed string itself */

    free (buf);
  }
  else pack_one (isize, i, ints, 0, int);
}

/* pack single double */
void pack_double (int *dsize, double **d, int *doubles, double input)
{
  pack_one (dsize, d, doubles, input, double);
}

/* pack single int */
void pack_int (int *isize, int **i, int *ints, int input)
{
  pack_one (isize, i, ints, input, int);
}

/* unpack vector of doubles */
void unpack_doubles (int *dpos, double *d, int doubles, double *output, int count)
{
  unpack_many (dpos, d, doubles, output, count, double);
}

/* unpack vector of ints */
void unpack_ints (int *ipos, int *i, int ints, int *output, int count)
{
  unpack_many (ipos, i, ints, output, count, int);
}

/* unpack string */
char* unpack_string (int *ipos, int *i, int ints)
{
  int len, n, *buf;
  char *str = NULL;

  unpack_one (ipos, i, ints, len, int);

  if (len > 0)
  {
    ERRMEM (buf = malloc (len * sizeof (int)));
    ERRMEM (str = malloc (len * sizeof (char)));

    unpack_many (ipos, i, ints, buf, len, int);

    for (n = 0; n < len; n ++) str [n] = (char) buf [n];

    free (buf);
  }

  return str;
}

/* unpack single double */
double unpack_double (int *dpos, double *d, int doubles)
{
  double output;

  unpack_one (dpos, d, doubles, output, double);

  return output;
}

/* unpack single int */
int unpack_int (int *ipos, int *i, int ints)
{
  int output;

  unpack_one (ipos, i, ints, output, int);

  return output;
}
