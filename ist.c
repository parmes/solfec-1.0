/*
 * ist.c
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

#include <stdlib.h>
#include "ist.h"

enum {red, black};
#define NIL &set_sentinel
static ISET set_sentinel = {NULL, NIL, NIL, black, 0};

inline static void set_rotate_l (ISET **root, ISET *x)
{
  ISET *y;

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

inline static void set_rotate_r (ISET **root, ISET *x)
{
  ISET *y;

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

inline static ISET* item (MEM *pool, int value)
{
  ISET *x;

  if (!(x = MEM_Alloc (pool))) return NULL;

  x->l = x->r = NIL;
  x->value = value;
  x->p = NULL;

  return x;
}

/* insert an element into the set */
ISET* ISET_Insert (MEM *pool, ISET **root, int value)
{
  ISET *y, *x, *node;
  
  if ((*root) == NULL || (*root) == NIL)
  {
    (*root) = x = item (pool, value);
  }
  else
  {
    node = *root;
    
    while (1)
    {
      if (value < node->value)
      {
	if (node->l == NIL)
	{
	  if (!(x = item (pool, value))) return NULL;
	  node->l = x;
	  x->p = node;
	  break;
	}
	else node = node->l;
      }
      else if (value > node->value)
      {
	if (node->r == NIL)
	{
	  if (!(x = item (pool, value))) return NULL;
	  node->r = x;
	  x->p = node;
	  break;
	}
	else node = node->r;
      }
      else return node;
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

/* postorder traverse and free set memory */
void ISET_Free (MEM *pool, ISET **root)
{
  if (*root == NULL || *root == NIL) return;
  ISET_Free (pool, &(*root)->l);
  ISET_Free (pool, &(*root)->r);
  MEM_Free (pool, *root);
  *root = NULL;
}

/* first element */
ISET* ISET_First (ISET *root)
{
  if (root == NULL || root == NIL) return NULL;
  while (root->l != NIL) root = root->l;
  return root;
}

/* next element */
ISET* ISET_Next (ISET *node)
{
  ISET *y;
  
  if (node == NULL)
    return NULL;

  if (node->r == NIL)
  {

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
  else
  {
    y = node->r;
    while (y->l != NIL) y = y->l;
    return y;
  }
}
