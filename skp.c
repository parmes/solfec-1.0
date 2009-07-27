/*
 * skp.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ------------------------------------------------------------
 * skip list
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

#include "skp.h"

#define CHUNK 64

/* pseudorandom double from [0, 1) */
#define DRAND() ((double)rand()/((double)RAND_MAX + 1.0))

/* create new list item */
static ITEM* skip_item (LIST *list, void *key, void *data)
{
  ITEM *item;

  if (!(item = MEM_Alloc (&list->mem))) return NULL;

  item->key = key;
  item->data = data;
  item->forward = (ITEM**) (item + 1); /* TODO: rather than this, length of 'forward' should vary depending on the number of stored items;
                                          TODO: 'forward' could be 'malloc'ed here and 'realloc'ated later, when the number of items exceeds
                                          TODO: consecutive powers of two; the drawback would be in greater memory fragmentation */

  return item;
}

/* pick random list level */
static int random_level (int maxlevel)
{
  int n = 0;

  while (DRAND() < 0.5 && n < maxlevel) n ++;

  return n;
}

/* ========== INTERFACE ========== */

/* levels: number of skip levels
 * compare: item key comparison (NULL implies pointer comparison)
 * return: skip list or NULL when out of memory;
 * create skip list */
LIST* LIST_Create (int levels, LIST_Compare compare)
{
  LIST *list;

  if (!(list = malloc (sizeof (LIST)))) return NULL;
  MEM_Init (&list->mem, sizeof (ITEM) + levels * sizeof (ITEM*), CHUNK);
  list->maxlevel = levels - 1;
  list->level = 0;
  list->size = 0;

  if (!(list->header.forward = calloc (levels, sizeof (ITEM*)))) return NULL;
  if (!(list->update = malloc (levels * sizeof (ITEM*)))) return NULL;

  list->compare = compare;

  return list;
}

/* key: item key
 * data: item data
 * return: skip item or NULL when out of memory;
 * insert data item */
ITEM* LIST_Insert (LIST *list, void *key, void *data)
{
  ITEM *header, **update, *x;
  int maxlevel, level, i, n;
  LIST_Compare compare;

  maxlevel = list->maxlevel;
  compare = list->compare;
  header = &list->header;
  update = list->update;
  level = list->level;
  x = header;

  for (i = level; i >= 0; i --)
  {
    if (compare) while (x->forward [i] && compare (x->forward [i]->key, key) < 0) x = x->forward [i];
    else while (x->forward [i] && x->forward [i]->key < key) x = x->forward [i];

    update [i] = x;
  }

  x = x->forward [0];

  if (x && ((compare && compare (x->key, key) == 0) || ((!compare) && x->key == key))) x->data = data;
  else
  {
    n = random_level (maxlevel);

    if (n > level)
    {
      for (i = level + 1; i <= maxlevel; i ++) update [i] = header;

      list->level = n;
    }

    if(!(x = skip_item (list, key, data))) return NULL;

    for (i = 0; i <= n; i ++)
    {
      x->forward [i] = update [i]->forward [i];
      update [i]->forward [i] = x;
    }

    list->size ++;
  }

  return x;
}

/* key: item key;
 * delete data item */
void  LIST_Delete (LIST *list, void *key)
{
  ITEM *header, **update, *x;
  LIST_Compare compare;
  int level, i;

  compare = list->compare;
  header = &list->header;
  update = list->update;
  level = list->level;
  x = header;

  for (i = level; i >= 0; i --)
  {
    if (compare) while (x->forward [i] && compare (x->forward [i]->key, key) < 0) x = x->forward [i];
    else while (x->forward [i] && x->forward [i]->key < key) x = x->forward [i];

    update [i] = x;
  }

  x = x->forward [0];

  if (x && ((compare && compare (x->key, key) == 0) || ((!compare) && x->key == key)))
  {
    for (i = 0; i <= level; i ++)
    {
      if (update [i]->forward [i] != x) break;
      update [i]->forward [i] = x->forward [i];
    }

    MEM_Free (&list->mem, x);

    while (level > 0 && header->forward [level] == NULL) level --;
    list->level = level;


    list->size --;
  }
}

/* key: item key
 * return: item data or NULL when not found;
 * find data item */
void* LIST_Find (LIST *list, void *key)
{
  LIST_Compare compare;
  ITEM *x;
  int i;

  for (x = &list->header, i = list->level, compare = list->compare; i >= 0; i --)
  {
    if (compare) while (x->forward [i] && compare (x->forward [i]->key, key) < 0) x = x->forward [i];
    else while (x->forward [i] && x->forward [i]->key < key) x = x->forward [i];
  }

  x = x->forward [0];

  if (x && ((compare && compare (x->key, key) == 0) || ((!compare) && x->key == key))) return x->data;
  else return NULL;
}

/* return list size */
int  LIST_Size (LIST *list)
{
  return list->size;
}

/* destroy list */
void  LIST_Destroy (LIST *list)
{
  MEM_Release (&list->mem);
  free (list->header.forward);
  free (list->update);
  free (list);
}
