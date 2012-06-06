/*
 * skp.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * skip list based on [1]
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

#ifndef __skp__
#define __skp__

typedef struct skip_item ITEM;
typedef struct skip_list LIST;
typedef int (*LIST_Compare) (void*, void*);

struct skip_item
{
  void *key,
       *data;

  ITEM **forward;
};

struct skip_list
{
  MEM mem;

  int maxlevel,
      level,
      size;

  ITEM header;

  ITEM **update;

  LIST_Compare compare;
};

/* levels: number of skip levels
 * compare: item key comparison (NULL implies pointer comparison)
 * return: skip list or NULL when out of memory;
 * create skip list */
LIST* LIST_Create (int levels, LIST_Compare compare);

/* key: item key
 * data: item data
 * return: skip item or NULL when out of memory;
 * insert data item */
ITEM* LIST_Insert (LIST *list, void *key, void *data);

/* key: item key;
 * delete data item */
void  LIST_Delete (LIST *list, void *key);

/* key: item key
 * return: item data or NULL when not found;
 * find data item */
void* LIST_Find (LIST *list, void *key);

/* return list size */
int  LIST_Size (LIST *list);

/* destroy list */
void  LIST_Destroy (LIST *list);

/* get first item */
#define LIST_First(list) (list)->header.forward [0]

/* get next item */
#define LIST_Next(item) (item)->forward [0]

/* ========== REFERENCES ==========
 *
 * [1] William, Pugh. Skip lists: a probabilistic alternative to balanced
 *     trees.  Communications of the ACM 33 (6): 668–676, 1990.
 */

#endif
