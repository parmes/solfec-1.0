/*
 * mat.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * bulk material
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

#include <stdio.h>
#include "mem.h"
#include "map.h"

#ifndef __mat__
#define __mat__

typedef struct bulkmat BULK_MATERIAL;
typedef struct matset MATSET;

struct bulkmat
{
  char *label;

  enum
  { KIRCHHOFF } model;

  double young,
         poisson,
         density;
};

struct matset
{
  MEM matmem,
      mapmem;

  MAP *map;   /* label based map */

  int size;   /* number of materials */
};

/* create bulk material set */
MATSET* MATSET_Create ();

/* insert new material */
BULK_MATERIAL* MATSET_Insert (MATSET *set, char *label, BULK_MATERIAL data);

/* find by label */
BULK_MATERIAL* MATSET_Find (MATSET *set, char *label);

/* release memory */
void MATSET_Destroy (MATSET *set);

/* export MBFCP definition */
void MATSET_2_MBFCP (MATSET *set, FILE *out);

#endif
