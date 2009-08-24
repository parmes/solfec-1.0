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
#include "put.h"

/* get statistics on single integer variable; return 1 for rank 0 and 0 for others */
int PUT_root_int_stats (int val, int *sum, int *min, int *avg, int *max)
{
  int rank, ncpu, *all;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &ncpu);

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

/* parallel timer end: get maximum of all calls; return 1 for rank 0 and 0 for others */
int PUT_root_timerend (TIMING *t, double *time)
{
  double local, timing;
  int rank;

  local = timerend (t);

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  if (!time) time = &timing;

  MPI_Reduce (&local, time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (rank == 0) return 1;
  else return 0;
}

/* parallel timer end: return maximum of all calls */
double PUT_timerend (TIMING *t)
{
  double local, time;

  local = timerend (t);

  MPI_Allreduce (&local, &time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return time;
}

/* return minimum of all calls */
int PUT_int_min (int val)
{
  int ret;

  MPI_Allreduce (&val, &ret, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

  return ret;
}
