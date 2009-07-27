/*
 * set.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * rb-tree based set container
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

#ifndef __set__
#define __set__

typedef struct set SET; /* set type */
typedef int (*SET_Compare) (void*, void*); /* comparison callback */

struct set
{
  SET *p, *l, *r; /* parent, left, right */
  char colour; /* red, blue */
  void *data; /* user data */
};

/* insert an element into the set (return NULL if the data already is in the set) */
SET* SET_Insert (MEM *pool, SET **root,
  void *data, SET_Compare compare); /* 'compare' == NULL 
					=> direct pointer 'data' comparison 
					(this applies below as well) */

/* find an element */
void* SET_Find (SET *root, void *data , SET_Compare compare);

/* check if an element is already in the set */
int SET_Contains (SET *root, void *data , SET_Compare compare);

/* delete an element from the set */
void SET_Delete (MEM *pool, SET **root, void *data, SET_Compare compare);

/* postorder traverse and free set memory */
void SET_Free (MEM *pool, SET **root);

/* return number of items */
int SET_Size (SET *root);

/* first element */
SET* SET_First (SET *root);

/* last element */
SET* SET_Last (SET *root);

/* previous element */
SET* SET_Prev (SET *node);

/* next element */
SET* SET_Next (SET *node);

#endif
