/*
 * xyt.c
 * Copyright (C) 2005, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * Priority search tree
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
#include "mem.h"
#include "xyt.h"
#include "err.h"

enum { node, leaf };

enum { red, black };

/* lexicographical comparision for the binary tree structure */
inline static int compare (double x, void *key, XYTREE *nod)
{
  if (x < nod->x) return -1;
  else if (x > nod->x ) return 1;
  else if (key < nod->key) return -1;
  else if (key > nod->key) return 1;
  else return 0;
}

/* push a link down the heap structure */
inline static void push (XYTREE *tree, XYTREE *link)
{
   XYTREE *tmp;

   while (tree->link)
   {
     if (tree->link->y < link->y)
     {
       tmp = tree->link;
       tree->link = link;
       link = tmp;
     }

     if (compare (link->x, link->key, tree) <= 0)
       tree = tree->left;
     else tree = tree->right;
   }

   tree->link = link;
}

/* pull all links from a sobtree up */
inline static void pull (XYTREE *tree)
{
  XYTREE *tmp;

  do
  {
    if (tree->type == leaf)
    {
      tree->link = NULL;
      return;
    }
    else tmp = tree->left;
    
    if (tmp->link == NULL ||
       (tree->right->link &&
       (tree->right->link->y > tmp->link->y))) tmp = tree->right;
      
    tree->link = tmp->link;
    tree = tmp;
    
  } while (tree->link);
}

/* create a NIL node */
inline static XYTREE* nil (MEM *pool)
{
  XYTREE *x;

  ERRMEM (x = MEM_Alloc (pool));
  x->type = leaf;
  x->left = x->right = x;
  x->colour = black;
  return x;
}

/* swap links in two nodes */
inline static void swap_links (XYTREE *x, XYTREE *y)
{
  XYTREE *tmp;

  tmp = x->link;
  x->link = y->link;
  y->link = tmp;
}

/* standard tree rotation + heap maintenance */
inline static void rotate_left (XYTREE **root, XYTREE *x)
{
  XYTREE *y, *z;

  /* rotate */
  y = x->right;
  x->right = y->left;
  y->left->parent = x;
  y->parent = x->parent;
  if (x->parent == NULL) *root = y;
  else if (x == x->parent->left)
    x->parent->left = y;
  else x->parent->right = y;
  y->left = x;
  x->parent = y;

  /* update heap */
  swap_links (x, y);
  z = x->link;
  pull (x);
  if (z) push (y, z);
}

inline static void rotate_right (XYTREE **root, XYTREE *x)
{
  XYTREE *y, *z;

  y = x->left;
  x->left = y->right;
  y->right->parent = x;
  y->parent = x->parent;
  if (x->parent == NULL) *root = y;
  else if (x == x->parent->right)
    x->parent->right = y;
  else x->parent->left = y;
  y->right = x;
  x->parent = y;

  swap_links (x, y);
  z = x->link;
  pull (x);
  if (z) push (y, z);
}

/* standard delete fixup */
inline static void delete_fixup (XYTREE **root, XYTREE *x)
{
  XYTREE *y;

  while (x != *root && x->colour == black)
  {
    if (x == x->parent->left)
    {
      y = x->parent->right;

      if (y->colour == red)
      {
        x->parent->colour = red;
        y->colour = black;
	rotate_left (root, x->parent);
	y = x->parent->right;
      }

      if (y->right->colour == black && y->left->colour == black)
      {
         y->colour = red;
	 x = x->parent;
      }
      else
      {
        if (y->right->colour == black)
	{
	  y->left->colour = black;
	  y->colour = red;
	  rotate_right (root, y);
	  y = x->parent->right;
	}

	y->colour = x->parent->colour;
	x->parent->colour = black;
	y->right->colour = black;
	rotate_left (root, x->parent);
	x = *root;
      }
    }
    else
    {
      y = x->parent->left;

      if (y->colour == red)
      {
        x->parent->colour = red;
        y->colour = black;
	rotate_right (root, x->parent);
	y = x->parent->left;
      }

      if (y->left->colour == black && y->right->colour == black)
      {
         y->colour = red;
	 x = x->parent;
      }
      else
      {
        if (y->left->colour == black)
	{
	  y->right->colour = black;
	  y->colour = red;
	  rotate_left (root, y);
	  y = x->parent->left;
	}

	y->colour = x->parent->colour;
	x->parent->colour = black;
	y->left->colour = black;
	rotate_right (root, x->parent);
	x = *root;
      }
    }
  }
  x->colour = black;
}

/* after-insert red-black restructuring */
inline static void restore (XYTREE **root, XYTREE *x)
{
  XYTREE *y;

  x->colour = red;
  while (x != *root && x->parent->colour == red)
  {
    if (x->parent == x->parent->parent->left)
    {
      y = x->parent->parent->right;
      if (y->colour == red)
      {
        x->parent->colour = black;
        y->colour = black;
        x->parent->parent->colour = red;
        x = x->parent->parent;
      }
      else
      {
        if (x == x->parent->right)
	{
          x = x->parent;
	  rotate_left (root, x);
        }
        x->parent->colour = black;
        x->parent->parent->colour = red;
	rotate_right (root, x->parent->parent);
      }
    }
    else
    {
      y = x->parent->parent->left;
      if (y->colour == red)
      {
        x->parent->colour = black;
        y->colour = black;
        x->parent->parent->colour = red;
        x = x->parent->parent;
      }
      else
      {
        if (x == x->parent->left)
	{
          x = x->parent;
	  rotate_right (root, x);
        }
        x->parent->colour = black;
        x->parent->parent->colour = red;
	rotate_left (root, x->parent->parent);
      }
    }
  }
  (*root)->colour = black;
}

/* tree insert */
void XY_Insert (MEM *pool, XYTREE **tree, double x, double y, void *key, void *data)
{
  XYTREE *p, *q;
  int ret;

  if (*tree == NULL)
  {
    p = nil (pool);
    p->x = x;
    p->y = y;
    p->key = key;
    p->data = data;
    p->link = p;
    *tree = p;
    return;
  }

  q = *tree;

  while (1)
  {
    ret = compare (x, key, q);

    if (ret < 0)
    {
      if (q->left->type == leaf)
      {
        p = q->left;
	break;
      }
      else q = q->left;
    }
    else if (ret > 0)
    {
      if (q->right->type == leaf)
      {
        p = q->right;
	break;
      }
      else q = q->right;
    }
    else return;
  }

  ret = compare (x, key, p);

  if (ret == 0) return;

  p->parent = (p == q ? NULL : q);
  p->type = node;
  p->left = nil (pool);
  p->right = nil (pool);
  
  if (ret < 0)
  {
    p->left->x = x;
    p->left->y = y;
    p->left->key = key;
    p->left->data = data;
      
    p->right->x = p->x;
    p->right->y = p->y;
    p->right->key = p->key;
    p->right->data = p->data;

    for (q = p; q->link != p; q = q->parent);
    q->link = p->right;
    
    p->x = x;
    p->y = y;
    p->key = key;
    p->data = data;

    push (*tree, p->left);
  }
  else
  {
    p->right->x = x;
    p->right->y = y;
    p->right->key = key;
    p->right->data = data;

    p->left->x = p->x;
    p->left->y = p->y;
    p->left->key = p->key;
    p->left->data = p->data;

    for (q = p; q->link != p; q = q->parent);
    q->link = p->left;
    
    push (*tree, p->right);
  }

  restore (tree, p);
}

/* deletion */
void XY_Delete (MEM *pool, XYTREE **tree, double x, double y, void *key)
{
  XYTREE *p, *q, *r;

  r = *tree;
  if (r == NULL) return;
  p = NULL;
  q = NULL;

  while (r->type != leaf)
  {
    if (r->link && compare (x, key, r->link) == 0) p = r;

    q = r;

    if (compare (x, key, r) <= 0) r = r->left;
    else r = r->right;
  }

  if (compare (x, key, r) != 0) return; /* leaf not found */

  if (p) pull (p);

  if (q)
  {
    if (r == q->left) p = q->right;
    else p = q->left;

    p->parent = q->parent;
    if (q->parent == NULL) *tree = p;
    else if (q == q->parent->left) q->parent->left = p;
    else q->parent->right = p;

    if (q->link && q->link != q)
      push (p, q->link); /* care for the heap now */

    if (q->colour == black)
      delete_fixup (tree, p);

    MEM_Free (pool, q);
  }
  else *tree = NULL;

  MEM_Free (pool, r);
}

/* report all the cobtree, as we are sure that data lays in the
 * scope of our interest domain */
static void reportdownto (XYTREE *tree, double ymin, void *pointer, XY_Callback callback)
{
  if (tree->link == NULL) return;
  else if (tree->link->y < ymin) return;

  callback (pointer, tree->link->data);

  if (tree->type == leaf) return;

  reportdownto (tree->left, ymin, pointer, callback);
  reportdownto (tree->right, ymin, pointer, callback);
}

/* two-sided query */
void XY_Query (XYTREE *tree, double xmin, double ymin, void *pointer, XY_Callback callback)
{
  XYTREE *prev = NULL;

  while (tree != prev && tree->link)
  {
    prev = tree;

    if (tree->link->x > xmin &&
        tree->link->y > ymin)
	  callback (pointer, tree->link->data);

    if (xmin <= tree->x)
    {
      if (tree != tree->right) reportdownto (tree->right, ymin, pointer, callback);
      tree = tree->left;
    }
    else tree = tree->right;
  };
}
