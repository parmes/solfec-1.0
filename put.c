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

/* get statistics on a vector of integer variables */
void PUT_int_stats (int n, int *val, int *sum, int *min, int *avg, int *max)
{
  int ncpu, *all, i, j, va, su, mi, av, ma;

  MPI_Comm_size (MPI_COMM_WORLD, &ncpu);

  ERRMEM (all = malloc (sizeof (int [n * ncpu])));

  MPI_Allgather (val, n, MPI_INT, all, n, MPI_INT, MPI_COMM_WORLD);

  for (j = 0; j < n; j ++)
  {
    for (su = i = 0, ma = INT_MIN, mi = INT_MAX; i < ncpu; i ++)
    {
      va = all [i*n + j];
      su += va, ma = MAX (ma, va), mi = MIN (mi, va);
    }

    av = su / ncpu;

    if (sum) sum [j] = su;
    if (min) min [j] = mi;
    if (avg) avg [j] = av;
    if (max) max [j] = ma;
  }

  free (all);
}

/* get statistics on a vector of integer variables; return 1 for rank 0 and 0 for others */
int PUT_root_int_stats (int n, int *val, int *sum, int *min, int *avg, int *max)
{
  int rank, ncpu, *all;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &ncpu);

  if (rank == 0) { ERRMEM (all = malloc (sizeof (int [n * ncpu]))); }
  else all = NULL;

  MPI_Gather (val, n, MPI_INT, all, n, MPI_INT, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    int i, j, va, su, mi, av, ma;

    for (j = 0; j < n; j ++)
    {
      for (su = i = 0, ma = INT_MIN, mi = INT_MAX; i < ncpu; i ++)
      {
	va = all [i*n + j];
	su += va, ma = MAX (ma, va), mi = MIN (mi, va);
      }

      av = su / ncpu;

      if (sum) sum [j] = su;
      if (min) min [j] = mi;
      if (avg) avg [j] = av;
      if (max) max [j] = ma;
    }

    free (all);

    return 1;
  }

  return 0;
}

/* get statistics on a vector of double variables; return 1 for rank 0 and 0 for others */
int PUT_root_double_stats (int n, double *val, double *sum, double *min, double *avg, double *max)
{
  int rank, ncpu;
  double *all;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &ncpu);

  if (rank == 0) { ERRMEM (all = malloc (sizeof (double [n * ncpu]))); }
  else all = NULL;

  MPI_Gather (val, n, MPI_DOUBLE, all, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (rank == 0)
  {
    double va, su, mi, av, ma;
    int i, j;

    for (j = 0; j < n; j ++)
    {
      for (su = i = 0, ma = INT_MIN, mi = INT_MAX; i < ncpu; i ++)
      {
	va = all [i*n + j];
	su += va, ma = MAX (ma, va), mi = MIN (mi, va);
      }

      av = su / (double) ncpu;

      if (sum) sum [j] = su;
      if (min) min [j] = mi;
      if (avg) avg [j] = av;
      if (max) max [j] = ma;
    }

    free (all);

    return 1;
  }

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

/* MPI operator callback used to find minimal value rank */
static void int_min_rank (int *in, int *inout, int *len, MPI_Datatype *type)
{
  int i;

  for (i = 0; i < *len; i ++)
  {
    if (*in < *inout)
    {
      inout [0] = in [0];
      inout [1] = in [1];
    }

    inout += 2;
    in += 2;
  }
}

/* return minimum of all calls and its rank */
int PUT_int_min_rank (int val, int *rank)
{
  int in [2], out [2];
  MPI_Datatype type;
  MPI_Op op;

  in [0] = val;
  MPI_Comm_rank (MPI_COMM_WORLD, &in [1]);

  MPI_Type_contiguous (2, MPI_INT, &type);
  MPI_Type_commit (&type);
  MPI_Op_create ((MPI_User_function*)int_min_rank, 1, &op);
  MPI_Allreduce (in, out, 1, type, op, MPI_COMM_WORLD); /* compute (min(val), rank(min(val))) in out */
  MPI_Type_free (&type);
  MPI_Op_free (&op);

  if (rank) *rank = out [1];

  return out [0];
}
