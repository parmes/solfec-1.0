/*
 * set.c
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

#include <stdlib.h>
#include <stdio.h>
#include "err.h"
#include "set.h"

enum {red, black};
#define NIL &set_sentinel
static SET set_sentinel = { NULL, NIL, NIL, black, NULL };

inline static void set_rotate_l (SET **root, SET *x)
{
  SET *y;

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

inline static void set_rotate_r (SET **root, SET *x)
{
  SET *y;

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

inline static void set_delete_fixup (SET **root, SET *x)
{
  SET *y;

  while (x != *root && x->colour == black)
  {
    if (x == x->p->l)
    {
      y = x->p->r;

      if (y->colour == red)
      {
        x->p->colour = red;
        y->colour = black;
        set_rotate_l (root, x->p);
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
	  set_rotate_r (root, y);
	  y = x->p->r;
	}

	y->colour = x->p->colour;
	x->p->colour = black;
	y->r->colour = black;
	set_rotate_l (root, x->p);
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
        set_rotate_r (root, x->p);
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
	  set_rotate_l (root, y);
	  y = x->p->l;
	}

	y->colour = x->p->colour;
	x->p->colour = black;
	y->l->colour = black;
	set_rotate_r (root, x->p);
	x = *root;
      }
    }
  }
  x->colour = black;
}

static void set_size (SET *node, int *size)
{
  if (node != NULL && node != NIL)
  {
    set_size (node->l, size);
    (*size) ++;
    set_size (node->r, size);
  }
}

inline static SET* item (MEM *pool, void *data)
{
  SET *x;

  if (pool) { ERRMEM (x = MEM_Alloc (pool)); }
  else { ERRMEM (x = MEM_CALLOC (sizeof (SET))); }

  x->l = x->r = NIL;
  x->p = NULL;
  x->data = data;

  return x;
}

SET* SET_Insert (MEM *pool, SET **root, void *data, SET_Compare compare)
{
  SET *y, *x, *node;
  int cmp;
  
  if ((*root) == NULL || (*root) == NIL)
  {
    (*root) = x = item (pool, data);
  }
  else
  {
    node = *root;
    
    if (compare)
    {
      while (1)
      {
	cmp = compare (data, node->data);
	if (cmp < 0)
	{
	  if (node->l == NIL)
	  {
	    x = item (pool, data);
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
	    x = item (pool, data);
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
	if (data < node->data)
	{
	  if (node->l == NIL)
	  {
	    x = item (pool, data);
	    node->l = x;
	    x->p = node;
	    break;
	  }
	  else node = node->l;
	}
	else if (data > node->data)
	{
	  if (node->r == NIL)
	  {
	    x = item (pool, data);
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
	  set_rotate_l (root, x);
        }

        x->p->colour = black;
        x->p->p->colour = red;
        set_rotate_r (root, x->p->p);
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
          set_rotate_r (root, x);
        }
        x->p->colour = black;
        x->p->p->colour = red;
        set_rotate_l (root, x->p->p);
      }
    }
  }

  (*root)->colour = black;
  return node;
}

void* SET_Find (SET *root, void *data , SET_Compare compare)
{
  SET *node = root;
  int cmp;

  if (compare)
  {
    while (node != NIL && node != NULL)
    {
      cmp = compare (data, node->data);
      if (cmp == 0) return node->data;
      else if (cmp < 0) node = node->l;
      else node = node->r;
    }
  }
  else
  {
    while (node != NIL && node != NULL)
    {
      if (data == node->data) return node->data;
      else if (data < node->data) node = node->l;
      else node = node->r;
    }
  }

  return 0;
}

int SET_Contains (SET *root, void *data, SET_Compare compare)
{
  SET *node = root;
  int cmp;

  if (compare)
  {
    while (node != NIL && node != NULL)
    {
      cmp = compare (data, node->data);
      if (cmp == 0) return 1;
      else if (cmp < 0) node = node->l;
      else node = node->r;
    }
  }
  else
  {
    while (node != NIL && node != NULL)
    {
      if (data == node->data) return 1;
      else if (data < node->data) node = node->l;
      else node = node->r;
    }
  }

  return 0;
}

void SET_Delete (MEM *pool, SET **root,
  void *data, SET_Compare compare)
{
  SET *node = *root, *x, *y;
  int cmp;
 
  if (compare)
  {
    while (node != NIL && node != NULL)
    {
      cmp = compare (data, node->data);
      if (cmp == 0)
      {
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
      if (data == node->data)
      {
        goto FOUND;
      }
      else if (data < node->data) node = node->l;
      else node = node->r;
    }
  }
  return;
  
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
    node->data = y->data;
  }

  if (y->colour == black)
    set_delete_fixup (root, x);
      
  if (pool) MEM_Free (pool, y);
  else free (y);
}

SET* SET_Delete_Node (MEM *pool, SET **root, SET *node)
{
  SET *x, *y, *next;
 
  next = SET_Next (node);

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
    node->data = y->data;
    next = node;
  }

  if (y->colour == black)
    set_delete_fixup (root, x);
      
  if (pool) MEM_Free (pool, y);
  else free (y);

  ASSERT_DEBUG (next != y, "Returning a freed item");

  return next; /* the next item after the input 'node', ordered by data */
}

void SET_Free (MEM *pool, SET **root)
{
  if (*root == NULL || *root == NIL) return;
  SET_Free (pool, &(*root)->l);
  SET_Free (pool, &(*root)->r);
  if (pool) MEM_Free (pool, *root);
  else free (*root);
  *root = NULL;
}

SET* SET_Copy (MEM *pool, SET *root)
{
  SET *x;

  if (root == NULL || root == NIL) return root;
  x = item (pool, root->data);
  x->colour = root->colour;
  x->l = SET_Copy (pool, root->l);
  if (x->l != NIL) x->l->p = x;
  x->r = SET_Copy (pool, root->r);
  if (x->r != NIL) x->r->p = x;

  return x;
}

int SET_Size (SET *root)
{
  int size = 0;
  set_size (root, &size);
  return size;
}

SET* SET_First (SET *root)
{
  if (root == NULL || root == NIL) return NULL;
  while (root->l != NIL) root = root->l;
  return root;
}

SET* SET_Last (SET *root)
{
  if (root == NULL || root == NIL)
    return NULL;
  while (root->r != NIL) root = root->r;
  return root;
}

SET* SET_Prev (SET *node)
{
  SET *y;
  
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
SET* SET_Next (SET *node)
{
  SET *y;
  
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
