/*
 * cmptest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * kd-tree test
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
#include "kdt.h"
#include "alg.h"
#include "err.h"

int main (int argc, char **argv)
{
  double *p, *q, *points, d [6], l;
  int i, n, k, m;
  void **data;
  KDT *kd;

  if (argc > 1) n = atoi (argv [1]);
  
  n = MAX (128, n);

  ERRMEM (points = malloc (n * sizeof (double [3])));

  printf ("Generating %d random points ...\n", n);
  for (i = 0, p = points; i < n; i ++, p += 3)
  {
    p [0] = DRAND ();
    p [1] = DRAND ();
    p [2] = DRAND ();
  }

  printf ("Building kd-tree ... "), fflush (stdout);
  kd = KDT_Create (n, points, 0.0);
  if (kd) printf ("OK\n");
  else printf ("FAILED\n");

  printf ("Finding nearest neighbours of the same points ... "), fflush (stdout);
  for (i = 0, p = points; i < n; i ++, p += 3)
  {
    q = KDT_Nearest (kd, p);
    SUB (p, q, d);
    l = LEN (d);
    if (l != 0.0) break;
  }
  if (i == n) printf ("OK\n");
  else printf ("FAILED: |q-p| = %g\n", l);

  printf ("Droping %d point centered random size boxes ... ", n), fflush (stdout);
  for (i = 0, p = points; i < n; i ++, p += 3)
  {
    d [0] = p [0];
    d [1] = p [1];
    d [2] = p [2];
    d [3] = p [0] + 1E-3 * DRAND();
    d [4] = p [1] + 1E-3 * DRAND();
    d [5] = p [2] + 1E-3 * DRAND();

    KDT_Drop (kd, d, p);
  }
  printf ("DONE\n");

  printf ("Picking boxes for %d random points ... ", n), fflush (stdout);
  for (i = m = 0; i < n; i ++)
  {
    d [0] = DRAND ();
    d [1] = DRAND ();
    d [2] = DRAND ();

    KDT_Pick (kd, d, &data, &k);
    if (k) free (data);
    m += k;
  }
  printf ("%d boxes picked\n", m);

  KDT_Destroy (kd);
  free (points);

  return 0;
}
