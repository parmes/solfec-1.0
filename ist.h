/*
 * ist.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara@gmai.com)
 * -------------------------------------------------------
 * rb-tree based set of integers
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

#ifndef __ist__
#define __ist__

typedef struct iset ISET;

struct iset
{
  ISET *p, *l, *r; /* parent, left, right */
  char colour; /* red, blue */
  int value;
};

/* insert an element into the set */
ISET* ISET_Insert (MEM *pool, ISET **root, int value);

/* postorder traverse and free set memory */
void ISET_Free (MEM *pool, ISET **root);

/* first element */
ISET* ISET_First (ISET *root);

/* next element */
ISET* ISET_Next (ISET *node);

#endif
