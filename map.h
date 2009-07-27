/*
 * map.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * rb-tree based map container
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

#ifndef __map__
#define __map__

typedef struct map MAP; /* map type */
typedef int (*MAP_Compare) (void*, void*); /* comparison callback */

struct map
{
  MAP *p, *l, *r; /* parent, left, right */
  char colour; /* red, blue */
  void *key, *data; /* user data */
};

/* insert an element into the map (return NULL if the key is already mapped) */
MAP* MAP_Insert (MEM *pool, MAP **root,
  void *key, void *data, MAP_Compare compare); /* 'compare' == NULL 
						  => direct pointer 'key' comparison 
						  (this applies below as well) */
/* find data for a specific key */
void* MAP_Find (MAP *root, void *key, MAP_Compare compare);

/* find map node for a specific key */
MAP* MAP_Find_Node (MAP *root, void *key, MAP_Compare compare);

/* delete an element from the map */
void* MAP_Delete (MEM *pool, MAP **root, void *key, MAP_Compare compare);

/* delete a specific map node => return the next node by key */
MAP* MAP_Delete_Node (MEM *pool, MAP **root, MAP *node);

/* postorder traverse and free map memory */
void MAP_Free (MEM *pool, MAP **root);

/* return number of items */
int MAP_Size (MAP *root);

/* first element */
MAP* MAP_First (MAP *root);

/* last element */
MAP* MAP_Last (MAP *root);

/* previous element */
MAP* MAP_Prev (MAP *node);

/* next element */
MAP* MAP_Next (MAP *node);

#endif
