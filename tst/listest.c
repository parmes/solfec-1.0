/*
 * listest.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test list sorting
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
#include <time.h>
#include "lis.h"

typedef struct list LIST;

struct list
{
  int x;
  LIST *p, *n;
};

#define LE(i, j) ((i)->x <= (j)->x)
IMPLEMENT_LIST_SORT (DOUBLY_LINKED, merge_sort_list, LIST, p, n, LE)

static struct list* random_list (int length)
{
  struct list *l = NULL, *i;

  srand ((unsigned) time (NULL));

  for (; length > 0; length --)
  {
    i = malloc (sizeof (struct list));
    i->x = rand () % length;
    if (l) l->p = i;
    i->p = NULL;
    i->n = l;
    l = i;
  }

  return l;
}

int test_sorting (struct list *l)
{
  for (; l; l = l->n)
    if (l->p && l->p->x > l->x) return 0;
  return 1;
}

int main (int argc, char **argv)
{
  int length = 3;
  struct list *l, *i;

  if (argc > 1) 
  {
    if (atoi (argv [1]) > 0)
      length = atoi (argv [1]);
  }

  l = random_list (length);
  printf ("Random list of length %d:\n", length);
  for (i = l; i; i = i->n) printf ("%d\t", i->x);
  printf ("\n");
  l = merge_sort_list (l);
  printf ("Sorted list:\n");
  for (i = l; i; i = i->n) printf ("%d\t", i->x);
  printf ("\n");
  printf ("SORT => %s\n", test_sorting (l) ? "OK" : "FAILED");
}
