/*
 * pck.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * pack and unpack data into doubles and integers
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
    type *slot;\
  \
    if ((*size) < (*stored) + count)\
    {\
      (*size) = 2 * ((*stored) + count);\
      ERRMEM ((*array) = realloc ((*array), (*size) * sizeof (type)));\
    }\
  \
    for (slot = &(*array)[*stored], (*stored) += count;\
	 count > 0; count --, slot ++, input ++) *slot = *input;\
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
    type *slot;\
  \
    ASSERT (((*pos) + count) < stored, ERR_PCK_UNPACK);\
  \
    for (slot = &array[*pos], (*pos) += count;\
	 count > 0; count --, slot ++, output ++) *output = *slot;\
  }\
  while (0)

#define unpack_one(pos, array, stored, output, type)\
  do\
  {\
    ASSERT (((*pos) + 1) < stored, ERR_PCK_UNPACK);\
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
void pack_ints (int *isize, int **i, int *integers, int *input, int count)
{
  pack_many (isize, i, integers, input, count, int);
}

/* pack string into ints */
void pack_string (int *isize, int **i, int *integers, char *input)
{
  int len = strlen (input) + 1, /* + '\0' */
      bytlen = len * sizeof (char),
      remain = bytlen % sizeof (int),
      buflen = bytlen / sizeof (int) + (remain ? 1 : 0),
     *buf;

  if (input)
  {
    ERRMEM (buf = malloc (buflen * sizeof (int)));

    strcpy ((char*)buf, input);

    pack_one (isize, i, integers, buflen, int); /* length of the integer packed string */
    pack_many (isize, i, integers, buf, buflen, int); /* the integer packed string itself */

    free (buf);
  }
  else pack_one (isize, i, integers, 0, int);
}

/* pack single double */
void pack_double (int *dsize, double **d, int *doubles, double input)
{
  pack_one (dsize, d, doubles, input, double);
}

/* pack single int */
void pack_int (int *isize, int **i, int *integers, int input)
{
  pack_one (isize, i, integers, input, int);
}

/* unpack vector of doubles */
void unpack_doubles (int *dpos, double *d, int doubles, double *output, int count)
{
  unpack_many (dpos, d, doubles, output, count, double);
}

/* unpack vector of ints */
void unpack_ints (int *ipos, int *i, int integers, int *output, int count)
{
  unpack_many (ipos, i, integers, output, count, int);
}

/* unpack string */
char* unpack_string (int *ipos, int *i, int integers)
{
  int len,
     *buf = NULL;

  unpack_one (ipos, i, integers, len, int);

  if (len > 0)
  {
    ERRMEM (buf = malloc (len * sizeof (int)));
    unpack_many (ipos, i, integers, buf, len, int);
  }

  return (char*)buf;
}

/* unpack single double */
double unpack_double (int *dpos, double *d, int doubles)
{
  double output;

  unpack_one (dpos, d, doubles, output, double);

  return output;
}

/* unpack single int */
int unpack_int (int *ipos, int *i, int integers)
{
  int output;

  unpack_one (ipos, i, integers, output, int);

  return output;
}
