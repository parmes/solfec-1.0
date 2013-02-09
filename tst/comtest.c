/*
 * comtest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of parallel comunication
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
#include "com.h"
#include "pck.h"
#include "err.h"

#define INTS 17
#define DOUBLES 56

struct object
{
  int i [INTS];
  double d [DOUBLES];
};

typedef struct object OBJECT;

static int RANK, SIZE;

static int test_simple (int comall)
{
  const int isize = 32,
	    dsize = 64,
	    nsend = 100;

  int i [isize];

  double d [dsize];

  COMDATA send [nsend], *recv;

  int k, n, nrecv;

  for (k = 0; k < isize; k ++) i[k] = RANK;

  for (k = 0; k < dsize; k ++) d[k] = RANK;

  for (k = 0; k < nsend; k ++)
  {
    send[k].rank = (RANK + 1) % SIZE;
    send[k].ints = isize;
    send[k].doubles = dsize;
    send[k].i = i;
    send[k].d = d;
  }

  if (comall) COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);
  else COM (MPI_COMM_WORLD, 0, send, nsend, &recv, &nrecv);

  for (k = 0; k < nrecv; k ++)
  {
    for (n = 0; n < recv [k].ints; n ++)
      if (recv [k].i [n] != recv [k].rank) return 0;

    for (n = 0; n < recv [k].doubles; n ++)
      if (recv [k].d [n] != (double)recv [k].rank) return 0;
  }

  free (recv);

  return 1;
}

static int test_pattern (int comall)
{
  const int isize = 32,
	    dsize = 64,
	    nsend = 100,
	    nrep = 16;

  int i [isize];

  double d [dsize];

  COMDATA send [nsend], *recv;

  int k, n, m, nrecv;

  void *pattern;

  for (k = 0; k < nsend; k ++)
  {
    send[k].rank = (RANK + 1) % SIZE;
    send[k].ints = isize;
    send[k].doubles = dsize;
    send[k].i = i;
    send[k].d = d;
  }

  if (comall) pattern = COMALL_Pattern (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);
  else pattern = COM_Pattern (MPI_COMM_WORLD, 0, send, nsend, &recv, &nrecv);

  for (m = 0; m < nrep; m ++)
  {
    for (k = 0; k < isize; k ++) i[k] = RANK + m;

    for (k = 0; k < dsize; k ++) d[k] = RANK + m;

    if (comall) COMALL_Repeat (pattern);
    else COM_Repeat (pattern);

    for (k = 0; k < nrecv; k ++)
    {
      for (n = 0; n < recv [k].ints; n ++)
	if (recv [k].i [n] != (recv [k].rank + m)) return 0;

      for (n = 0; n < recv [k].doubles; n ++)
	if (recv [k].d [n] != (double)(recv [k].rank + m)) return 0;
    }
  }

  if (comall) COMALL_Free (pattern);
  else COM_Free (pattern);

  free (recv);

  return 1;
}

static void pack (OBJECT *obj, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_ints (isize, i, ints, obj->i, INTS);
  pack_doubles (dsize, d, doubles, obj->d, DOUBLES);
}

static OBJECT* unpack (void *data, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  OBJECT *obj;

  ERRMEM (obj = malloc (sizeof (OBJECT)));

  unpack_ints (ipos, i, ints, obj->i, INTS);
  unpack_doubles (dpos, d, doubles, obj->d, DOUBLES);

  return obj;
}

static int test_object ()
{
  const int nsend = 100;

  OBJECT obj [nsend];

  COMOBJ send [nsend], *recv;

  int k, n, nrecv;

  for (k = 0; k < nsend; k ++)
  {
    for (n = 0; n < INTS; n ++) obj[k].i[n] = RANK;
    for (n = 0; n < DOUBLES; n ++) obj[k].d[n] = RANK;

    send[k].rank = (RANK + 1) % SIZE;
    send[k].o = &obj[k];
  }

  COMOBJS (MPI_COMM_WORLD, 0, (OBJ_Pack)pack, NULL, (OBJ_Unpack)unpack, send, nsend, &recv, &nrecv);

  for (k = 0; k < nrecv; k ++)
  {
    for (n = 0; n < INTS; n ++)
    {
      OBJECT *o = recv [k].o;
      if (o->i [n] != recv [k].rank) return 0;
    }

    for (n = 0; n < DOUBLES; n ++)
    {
      OBJECT *o = recv [k].o;
      if (o->d [n] != (double)recv [k].rank) return 0;
    }
  }

  for (k = 0; k < nrecv; k ++) free (recv [k].o);

  free (recv);

  return 1;
}

static int test_objall ()
{
  const int nsend = 1;

  OBJECT obj [nsend];

  COMOBJ send [nsend], *recv;

  int k, n, nrecv;

  for (k = 0; k < nsend; k ++)
  {
    for (n = 0; n < INTS; n ++) obj[k].i[n] = RANK;
    for (n = 0; n < DOUBLES; n ++) obj[k].d[n] = RANK;

    send[k].rank = (RANK + 1) % SIZE;
    send[k].o = &obj[k];
  }

  COMOBJALL (MPI_COMM_WORLD, (OBJ_Pack)pack, NULL, (OBJ_Unpack)unpack, obj, &recv, &nrecv);

  for (k = 0; k < nrecv; k ++)
  {
    for (n = 0; n < INTS; n ++)
    {
      OBJECT *o = recv [k].o;
      if (o->i [n] != recv [k].rank) return 0;
    }

    for (n = 0; n < DOUBLES; n ++)
    {
      OBJECT *o = recv [k].o;
      if (o->d [n] != (double)recv [k].rank) return 0;
    }
  }

  for (k = 0; k < nrecv; k ++) free (recv [k].o);

  free (recv);

  return 1;
}

int main (int argc, char **argv)
{
  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &RANK);
  MPI_Comm_size (MPI_COMM_WORLD, &SIZE);

  if (test_simple (0)) printf ("Rank %d of %d COM OK\n", RANK, SIZE);
  else printf ("Rank %d of %d COM FAILED\n", RANK, SIZE);

  if (test_pattern (0)) printf ("Rank %d of %d COM_Pattern OK\n", RANK, SIZE);
  else printf ("Rank %d of %d COM_Pattern FAILED\n", RANK, SIZE);

  if (test_object ()) printf ("Rank %d of %d COMOBJ OK\n", RANK, SIZE);
  else printf ("Rank %d of %d COMOBJ FAILED\n", RANK, SIZE);

  if (test_simple (1)) printf ("Rank %d of %d COMALL OK\n", RANK, SIZE);
  else printf ("Rank %d of %d COMALL FAILED\n", RANK, SIZE);

  if (test_pattern (1)) printf ("Rank %d of %d COMALL_Pattern OK\n", RANK, SIZE);
  else printf ("Rank %d of %d COMALL_Pattern FAILED\n", RANK, SIZE);

  if (test_objall ()) printf ("Rank %d of %d COMOBJALL OK\n", RANK, SIZE);
  else printf ("Rank %d of %d COMOBJALL FAILED\n", RANK, SIZE);

  MPI_Finalize ();

  return 0;
}
