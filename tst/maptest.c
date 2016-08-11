/*
 * maptest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of MAP container
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
#include <math.h>
#include "map.h"
#include "alg.h"

static int map_test (int count)
{
  MAP *map, *item;
  MEM mem;
  int n, m, div;

  div = MAX (log (count), 2);

  MEM_Init (&mem, sizeof (MAP), 64);
  map = NULL;

  for (n = 0; n < count; n ++)
  {
    MAP_Insert (&mem, &map, (void*)(long)n, (void*)(long)n, NULL);
  }

  for (n = 0, item = MAP_First (map); item; n ++)
  {
    m = (int)item->key;

    if (rand () % div == 0 || n % div == 0)
    {
      item = MAP_Delete_Node (&mem, &map, item);
      if (item && (m+1) != (int)item->key) break;
      continue;
    }

    item = MAP_Next (item);
  }

  MEM_Release (&mem);

  return item ? 0 : 1;
}

int main (int argc, char **argv)
{
  int count = 128;

  if (argc > 1) count = atoi (argv [1]);

  if (count < 1) count = 1;

  if (map_test (count)) printf ("PASSED\n");
  else printf ("FAILED\n");

  return 0;
}
