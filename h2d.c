/*
 * h2d.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * 2D convex hull: Andrew's monotone chain algorithm
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
#include <float.h>
#include "mem.h"
#include "err.h"
#include "alg.h"
#include "lis.h"
#include "h2d.h"
#include "ext/predicates.h"

inline static int p2dcmp (P2D *a, P2D *b)
{
  if (a->coord [0] < b->coord [0]) return 1;
  else if (a->coord [0] == b->coord [0] && a->coord [1] <= b->coord [1]) return 1;
  else return 0;
}

/* quick circular point list sort */
IMPLEMENT_LIST_SORT (DOUBLY_LINKED, p2dsort, P2D, prev, next, p2dcmp)

/* input list containts valid 'point' values;
 * these are projected onto (p, normal) to give 'coord';
 * ouput list is made of elements of the input lists;
 * intput list structure gets destroyed upon return */
P2D* h2d (P2D *list, double *point, double *normal)
{
  double s [3], t [3], v [3], u [3], dot;
  P2D *minmin, *minmax, *maxmin, *maxmax;
  unsigned int depth = 0, save;
  P2D *p, *q, *r;

  /* project on (point, normal) plane */
  SET (s, 0.);
  for (q = list; q; q = q->next)
  {
    SUB (q->point, point, v);
    PRODUCT (v, normal, u); 
    dot = DOT (u, u);
    if (dot > GEOMETRIC_EPSILON*GEOMETRIC_EPSILON) break;
  }
  dot = 1./sqrt(dot);
  MUL (u, dot, s);
  PRODUCT (normal, s, t);
  NORMALIZE (t);
  for (q = list;  q; q = q->next)
  {
    SUB (q->point, point, v);
    q->coord [0] = DOT (v, s); 
    q->coord [1] = DOT (v, t); 
  }

  /* sort and remove duplicates */
  list = p2dsort (list);
  for (p = list, r = NULL; p; p = p->next)
  {
    if (r) r->next = p;
    else list = p;
    p->prev = r;
    r = p;
    while (p->next && (EQ (p->coord [0], p->next->coord [0]) && EQ (p->coord [1], p->next->coord [1]))) p = p->next;
  } r->next = NULL;

  /* get extremas */
  minmin = minmax = list;
  while (minmax->next &&
    minmax->next->coord [0] ==
    minmin->coord [0]) minmax = minmax->next;
  maxmin = maxmax = r;
  while (maxmin->prev &&
    maxmin->prev->coord [0] ==
    maxmax->coord [0]) maxmin = maxmin->prev;

  /* compute lower hull,
   * NULL-ify 'next' pointers
   * to mark not points in the lower hull,
   * therefore increase rubostness */
  p = minmax->next;
  minmin->next = NULL;
  list = minmin; depth ++;
  if (p) for (; p != maxmin; p = q)
  {
    if (orient2d (minmin->coord,
      maxmin->coord, p->coord) >= 0.) 
    { q = p->next; p->next = NULL; continue; }

    while (depth >= 2)
    {
      if (orient2d (list->next->coord,
        list->coord, p->coord) > 0.) break;
      q = list->next;
      list->next = NULL;
      list = q; depth --;
    }
    q = p->next;
    p->next = list;
    list = p; depth ++;
  }
  if (maxmin != minmin)
  { maxmin->next = list;
    list = maxmin; depth ++; }

  /* compute upper hull */
  save = depth;
  p = maxmin->prev;
  if (maxmax != maxmin)
  { maxmax->next = list;
    list = maxmax; depth ++; }
  if (p) for (; p != minmax; p = q)
  {
    if (orient2d (maxmax->coord, /* test also for not-NULL-ified 'next' */
      minmax->coord, p->coord) >= 0. || p->next)
    { q = p->prev; continue; }

    while ((depth - save) >= 2)
    {
      if (orient2d (list->next->coord,
        list->coord, p->coord) > 0.) break;
      list = list->next; depth --;
    }
    q = p->prev;
    p->next = list;
    list = p; depth ++;
  }
  if (minmax != minmin)
  { minmax->next = list;
    list = minmax; }

  /* again remove duplicates, XXX => do it once and well above */
  for (p = list, r = NULL; p; p = p->next)
  {
    if (r) r->next = p;
    else list = p;
    p->prev = r;
    r = p;
    while (p->next && (EQ (p->coord [0],p->next->coord [0]) && EQ (p->coord [1], p->next->coord [1]))) p = p->next;
  } r->next = NULL;

  return list;
}
