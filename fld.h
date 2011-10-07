/*
 * fld.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * scalar field (Python defined)
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

#include "mem.h"
#include "map.h"

#ifndef __fld__
#define __fld__

typedef struct field FIELD;
typedef struct fiset FISET;

struct field
{
  void *data; /* user data */
  void *call; /* user call */
  char *label;
};

struct fiset
{
  MEM fldmem,
      mapmem;

  MAP *map;   /* label based map */

  int size;   /* number of materials */
};

/* evaluate field */
double FIELD_Value (FIELD *fld, double x, double y, double z, double t);

/* create field set */
FISET* FISET_Create ();

/* insert new field */
FIELD* FISET_Insert (FISET *set, char *label, FIELD data);

/* find by label */
FIELD* FISET_Find (FISET *set, char *label);

/* release memory */
void FISET_Destroy (FISET *set);

#endif
