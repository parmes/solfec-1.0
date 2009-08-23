/*
 * put.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * parallel utilities
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
#include <limits.h>
#include <mpi.h>
#include "err.h"
#include "alg.h"

/* get statistics on single integer variable; return 1 for rank 0 and 0 for others */
int PUT_Int_Stats (int rank, int ncpu, int val, int *sum, int *min, int *avg, int *max)
{
  int *all;

  if (rank == 0) { ERRMEM (all = malloc (sizeof (int [ncpu]))); }
  else all = NULL;

  MPI_Gather (&val, 1, MPI_INT, all, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    int i, su, mi, av, ma;

    for (su = i = 0, ma = INT_MIN, mi = INT_MAX; i < ncpu; i ++)
    {
      su += all [i], ma = MAX (ma, all [i]), mi = MIN (mi, all [i]);
    }

    av = su / ncpu;

    if (sum) *sum = su;
    if (min) *min = mi;
    if (avg) *avg = av;
    if (max) *max = ma;

    return 1;
  }

  free (all);

  return 0;
}
