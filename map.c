/*
 * map.c
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

#include <stdlib.h>
#include <stdio.h>
#include "err.h"
#include "map.h"

enum {red, black};
#define NIL &map_sentinel
static MAP map_sentinel = { NULL, NIL, NIL, black, NULL, NULL };

inline static void map_rotate_l (MAP **root, MAP *x)
{
  MAP *y;

  y = x->r;
  x->r = y->l;
  if (y->l != NIL)
    y->l->p = x;
  if (y != NIL) y->p = x->p;
  if (x->p == NULL) *root = y;
  else if (x == x->p->l)
    x->p->l = y;
  else x->p->r = y;
  y->l = x;
  if (x != NIL) x->p = y;
}

inline static void map_rotate_r (MAP **root, MAP *x)
{
  MAP *y;

  y = x->l;
  x->l = y->r;
  if (y->r != NIL)
    y->r->p = x;
  if (y != NIL) y->p = x->p;
  if (x->p == NULL) *root = y;
  else if (x == x->p->r)
    x->p->r = y;
  else x->p->l = y;
  y->r = x;
  if (x != NIL) x->p = y;
}

inline static void map_delete_fixup (MAP **root, MAP *x)
{
  MAP *y;

  while (x != *root && x->colour == black)
  {
    if (x == x->p->l)
    {
      y = x->p->r;

      if (y->colour == red)
      {
        x->p->colour = red;
        y->colour = black;
        map_rotate_l (root, x->p);
	y = x->p->r;
      }
     
      if (y->r->colour == black && y->l->colour == black)
      {
         y->colour = red;
	 x = x->p;
      }
      else
      {
        if (y->r->colour == black)
	{
	  y->l->colour = black;
	  y->colour = red;
	  map_rotate_r (root, y);
	  y = x->p->r;
	}

	y->colour = x->p->colour;
	x->p->colour = black;
	y->r->colour = black;
	map_rotate_l (root, x->p);
	x = *root;
      }
    } 
    else
    {

      y = x->p->l;

      if (y->colour == red)
      {
        x->p->colour = red;
        y->colour = black;
        map_rotate_r (root, x->p);
	y = x->p->l;
      }
    
      if (y->l->colour == black && y->r->colour == black)
      {
         y->colour = red;
	 x = x->p;
      }
      else
      {
        if (y->l->colour == black)
	{
	  y->r->colour = black;
	  y->colour = red;
	  map_rotate_l (root, y);
	  y = x->p->l;
	}

	y->colour = x->p->colour;
	x->p->colour = black;
	y->l->colour = black;
	map_rotate_r (root, x->p);
	x = *root;
      }
    }
  }
  x->colour = black;
}

static void map_size (MAP *node, int *size)
{
  if (node != NULL && node != NIL)
  {
    map_size (node->l, size);
    (*size) ++;
    map_size (node->r, size);
  }
}

inline static MAP* item (MEM *pool, void *key, void *data)
{
  MAP *x;

  if (pool) { ERRMEM (x = MEM_Alloc (pool)); }
  else { ERRMEM (x = calloc (1, sizeof(MAP))); }

  x->l = x->r = NIL;
  x->p = NULL;
  x->key = key;
  x->data = data;

  return x;
}

MAP* MAP_Insert (MEM *pool, MAP **root, void *key, void *data, MAP_Compare compare)
{
  MAP *y, *x, *node;
  int cmp;

  if ((*root) == NULL || (*root) == NIL)
  {
    (*root) = x = item (pool, key, data);
  }
  else
  {
    node = *root;
    
    if (compare)
    {
      while (1)
      {
	cmp = compare (key, node->key);
	if (cmp < 0)
	{
	  if (node->l == NIL)
	  {
	    x = item (pool, key, data);
	    node->l = x;
	    x->p = node;
	    break;
	  }
	  else node = node->l;
	}
	else if (cmp > 0)
	{
	  if (node->r == NIL)
	  {
	    x = item (pool, key, data);
	    node->r = x;
	    x->p = node;
	    break;
	  }
	  else node = node->r;
	}
	else return NULL;
      }
    }
    else
    {
      while (1)
      {
	if (key < node->key)
	{
	  if (node->l == NIL)
	  {
	    x = item (pool, key, data);
	    node->l = x;
	    x->p = node;
	    break;
	  }
	  else node = node->l;
	}
	else if (key > node->key)
	{
	  if (node->r == NIL)
	  {
	    x = item (pool, key, data);
	    node->r = x;
	    x->p = node;
	    break;
	  }
	  else node = node->r;
	}
	else return NULL;
      }
    }
  }

  node = x;
  x->colour = red;
  while ((x != *root) && (x->p->colour == red))
  {
    if (x->p == x->p->p->l)
    {
      y = x->p->p->r;
      if (y->colour == red)
      {
        x->p->colour = black;
        y->colour = black;
        x->p->p->colour = red;
        x = x->p->p;
      }
      else
      {
        if (x == x->p->r)
	{
          x = x->p;
	  map_rotate_l (root, x);
        }

        x->p->colour = black;
        x->p->p->colour = red;
        map_rotate_r (root, x->p->p);
      }
    }
    else
    {
      y = x->p->p->l;
      if (y->colour == red)
      {
        x->p->colour = black;
        y->colour = black;
        x->p->p->colour = red;
        x = x->p->p;
      }
      else
      {
        if (x == x->p->l)
	{
          x = x->p;
          map_rotate_r (root, x);
        }
        x->p->colour = black;
        x->p->p->colour = red;
        map_rotate_l (root, x->p->p);
      }
    }
  }

  (*root)->colour = black;
  return node;
}

void* MAP_Find (MAP *root, void *key, MAP_Compare compare)
{
  MAP *node = root;
  int cmp;

  if (compare)
  {
    while (node != NIL && node != NULL)
    {
      cmp = compare (key, node->key);
      if (cmp == 0) return node->data;
      else if (cmp < 0) node = node->l;
      else node = node->r;
    }
  }
  else
  {
    while (node != NIL && node != NULL)
    {
      if (key == node->key) return node->data;
      else if (key < node->key) node = node->l;
      else node = node->r;
    }
  }

  return NULL;
}

MAP* MAP_Find_Node (MAP *root, void *key, MAP_Compare compare)
{
  MAP *node = root;
  int cmp;

  if (compare)
  {
    while (node != NIL && node != NULL)
    {
      cmp = compare (key, node->key);
      if (cmp == 0) return node;
      else if (cmp < 0) node = node->l;
      else node = node->r;
    }
  }
  else
  {
    while (node != NIL && node != NULL)
    {
      if (key == node->key) return node;
      else if (key < node->key) node = node->l;
      else node = node->r;
    }
  }

  return NULL;
}

void* MAP_Delete (MEM *pool, MAP **root, void *key, MAP_Compare compare)
{
  MAP *node = *root, *x, *y;
  void *data;
  int cmp;
 
  if (compare)
  {
    while (node != NIL && node != NULL)
    {
      cmp = compare (key, node->key);
      if (cmp == 0)
      {
        data = node->data;
        goto FOUND;
      }
      else if (cmp < 0) node = node->l;
      else node = node->r;
    }
  }
  else
  {
    while (node != NIL && node != NULL)
    {
      if (key == node->key)
      {
        data = node->data;
        goto FOUND;
      }
      else if (key < node->key) node = node->l;
      else node = node->r;
    }
  }
  return NULL;
  
FOUND:
  if (node->l == NIL || node->r == NIL)
  {
    y = node;
  }
  else
  {
    y = node->r;
    while (y->l != NIL) y = y->l;
  }
  
  if (y->l != NIL)
    x = y->l;
  else x = y->r;

  x->p = y->p;
  if (y->p)
    if (y == y->p->l)
      y->p->l = x;
    else
      y->p->r = x;
  else
    *root = x;

  if (y != node)
  {
    node->key = y->key;
    node->data = y->data;
  }

  if (y->colour == black)
    map_delete_fixup (root, x);
      
  if (pool) MEM_Free (pool, y);
  else free (y);

  return data;
}

MAP* MAP_Delete_Node (MEM *pool, MAP **root, MAP *node)
{
  MAP *x, *y;
 
  if (node->l == NIL || node->r == NIL)
  {
    y = node;
  }
  else
  {
    y = node->r;
    while (y->l != NIL) y = y->l;
  }
  
  if (y->l != NIL)
    x = y->l;
  else x = y->r;

  x->p = y->p;
  if (y->p)
    if (y == y->p->l)
      y->p->l = x;
    else
      y->p->r = x;
  else
    *root = x;

  if (y != node)
  {
    node->key = y->key;
    node->data = y->data;
  }
  else node = MAP_Next (node); /* as 'node == y' will be delted */

  if (y->colour == black)
    map_delete_fixup (root, x); /* this cannot change the next item after 'node' */
      
  if (pool) MEM_Free (pool, y);
  else free (y);

  return node; /* the next item after the input 'node', ordered by key */

  /* TODO: make sure that a correct node is returned */
}

void MAP_Free (MEM *pool, MAP **root)
{
  if (*root == NULL || *root == NIL) return;
  MAP_Free (pool, &(*root)->l);
  MAP_Free (pool, &(*root)->r);
  if (pool) MEM_Free (pool, *root);
  else free (*root);
  *root = NULL;
}

int MAP_Size (MAP *root)
{
  int size = 0;
  map_size (root, &size);
  return size;
}

MAP* MAP_First (MAP *root)
{
  if (root == NULL || root == NIL) return NULL;
  while (root->l != NIL) root = root->l;
  return root;
}

MAP* MAP_Last (MAP *root)
{
  if (root == NULL || root == NIL)
    return NULL;
  while (root->r != NIL) root = root->r;
  return root;
}

MAP* MAP_Prev (MAP *node)
{
  MAP *y;
  
  if (node == NULL)
    return NULL;

  if (node->l == NIL)  {

    if (node->p == NULL)
      return NULL; /* this is the only l node of root */
  
    if (node == node->p->r)
      return node->p; /* p is the predecessor */

    y = node->p;
    while (y->p && y == y->p->l)
      y = y->p; /* all the keys along this path are greater then 'node' key */

    if (y->p == NULL)
      return NULL; /* root is reached - node wast the very l tree element */

    return y->p; /* y was r descendant of y->p so the parent must be predecessor */
  }
  else {

    y = node->l;
    while (y->r != NIL) y = y->r;
    return y;
  }
}

/* symmetry of the above */
MAP* MAP_Next (MAP *node)
{
  MAP *y;
  
  if (node == NULL)
    return NULL;

  if (node->r == NIL)  {

    if (node->p == NULL)
      return NULL; /* this is the only r node of root */
  
    if (node == node->p->l)
      return node->p; /* p is the predecessor */

    y = node->p;
    while (y->p && y == y->p->r)
      y = y->p; /* all the keys along this path are smaller then 'node' key */

    if (y->p == NULL)
      return NULL; /* root is reached - node wast the very r tree element */

    return y->p; /* y was r descendant of y->p so the parent must be predecessor */
  }
  else {

    y = node->r;
    while (y->l != NIL) y = y->l;
    return y;
  }
}
