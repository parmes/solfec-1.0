/*
The MIT License (MIT)

Copyright (c) 2017 EDF Energy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* Purpose: co-rotated FEM displacements sampling */

/* Contributors: Tomasz Koziara */

#include <string.h>
#include "sol.h"
#include "bcd.h"
#include "lng.h"
#include "mem.h"
#include "err.h"
#include "put.h"
#include "fem.h"

#if MPI
/* unify subsets across all ranks; return current rank  */
static int subset_unify (DOM *dom, SET** subset)
{
  int rank, size, *sendbuf, sendcount, *recvbuf, *recvcounts, *displs, i, recvsize, *ptr, j, max_dofs;
  SET *item;
  BODY *bod;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &size);

  sendcount = SET_Size (*subset);
  ERRMEM (recvcounts = malloc(size * sizeof(int)));
  ERRMEM (displs = malloc(size * sizeof(int)));

  MPI_Allgather (&sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

  for (displs[0] = recvsize = i = 0; i < size; i ++)
  {
    if (i) displs[i] = recvsize;
    recvsize += recvcounts[i];
  }

  ERRMEM (recvbuf = malloc(recvsize * sizeof(int)));
  ERRMEM (sendbuf = malloc(sendcount * sizeof(int)));

  if (*subset)
  {
    bod = MAP_Find (dom->idb, (*subset)->data, NULL);

    if (!bod)
    {
      WARNING(0, "Inconsistent body mapping");
      return -1;
    }
    if (bod->kind != FEM)
    {
      WARNING(0, "A non-FEM body have been passed");
      return -1;
    }

    j = FEM_Mesh_Dofs(bod);
  }
  else j = 0;

  max_dofs = PUT_int_max (j); /* establish globally a value of the number of DOFs */

  for (item = SET_First (*subset), ptr = sendbuf; item; item = SET_Next (item), ptr ++)
  {
    ptr[0] = (int) (long) item->data;

    bod = MAP_Find (dom->idb, item->data, NULL);
    if (!bod)
    {
      WARNING(0, "Inconsistent body mapping");
      return -1;
    }
    if (bod->kind != FEM)
    {
      WARNING(0, "A non-FEM body have been passed");
      return -1;
    }
    if (FEM_Mesh_Dofs(bod) != max_dofs) /* test consistency of all mesh sizes */
    {
      WARNING(0, "A mismatch between FEM mesh sizes");
      return -1;
    }
  }

  MPI_Allgatherv (sendbuf, sendcount, MPI_INT, recvbuf, recvcounts, displs, MPI_INT, MPI_COMM_WORLD);

  for (i = 0; i < size; i ++)
  {
    for (j = 0; j < recvcounts[i]; j ++)
    {
      long value = recvbuf[displs[i]+j];
      SET_Insert (NULL, subset, (void*) value, NULL);
    }
  }

  free (sendbuf);
  free (recvbuf);
  free (displs);
  free (recvcounts);

  return rank;
}
#endif

/* sort doubles */
static void quick_sort (double a[], int n)
{
  double p, t;
  int i, j;

  if (n < 2) return;

  p = a[n/2];

  for (i = 0, j = n - 1;; i++, j--)
  {
    while (a[i] < p) i++;
    while (p < a[j]) j--;
    if (i >= j) break;

    t = a[i];
    a[i] = a[j];
    a[j] = t;
  }

  quick_sort (a, i);
  quick_sort (a+i, n-i);
}

/* check if any sampling value is within a time step from the current time */
static int sampling_within_step (double *sampling, int length, double time, double step)
{
  int lo = 0, hi = length;

  while (hi > lo+1)
  {
    int mid = lo + (hi-lo)/2;
    if (time > sampling[mid]) lo = mid;
    else hi = mid;
  }

  if (fabs (sampling[lo]-time) < step) return 1;
  else return 0;
}

/* append 'SOLFEC->bcd' with new definition and return the appended list head; return 1 upon succes; 0 otherwise */
int BCD_Append (SOLFEC *sol, SET *subset, double *sampling, int length, void *output)
{
  BCD *bcd;

  if (length > 0) quick_sort (sampling, length);
   
  ERRMEM (bcd = MEM_CALLOC (sizeof (BCD)));

  bcd->subset = subset;
  bcd->sampling = sampling;
  bcd->length = length;
  bcd->pending = NULL;
  bcd->output = output;
  bcd->latest = 0.0;

#if MPI
  int rank = subset_unify (sol->dom, &bcd->subset);
  if (rank == 0)
  {
    lng_xincref (output);
    bcd->output = output;
  }
  else bcd->output = NULL;
#else
  lng_xincref (output);
  bcd->output = output;
#endif

  if (!bcd->subset)
  {
    WARNING(0, "Invalid (empty) body subset in corotated displacements definition");
    return 0;
  }

  BODY *bod = MAP_Find (sol->dom->idb, subset->data, NULL);

  if (!bod)
  {
    WARNING(0, "Inconsistent body mapping");
    return -1;
  }
  if (bod->kind != FEM)
  {
    WARNING(0, "A non-FEM body have been passed");
    return -1;
  }

  bcd->size = FEM_Mesh_Dofs(bod);
  bcd->next = sol->bcd;
  sol->bcd = bcd;

#if MPI
  if (rank >= 0) return 1;
  else return 0;
#else
  for (SET *item = SET_First (subset); item; item = SET_Next(item))
  {
    bod = MAP_Find (sol->dom->idb, item->data, NULL);
    if (!bod)
    {
      WARNING(0, "Inconsistent body mapping");
      return 0;
    }
    if (bod->kind != FEM)
    {
      WARNING(0, "A non-FEM body have been passed");
      return 0;
    }
    if (FEM_Mesh_Dofs(bod) != bcd->size)
    {
      WARNING(0, "A mismatch between FEM mesh sizes");
      return 0;
    }
  }
  return 1;
#endif
}

/* sample 'bcd' list at current time --> goes into pending samples */
void BCD_Sample (SOLFEC *sol, BCD *bcd)
{
  for (; bcd; bcd = bcd->next)
  {
    int sample_now = 0;

    if (bcd->length > 0 && sampling_within_step (bcd->sampling,
        bcd->length, sol->dom->time, sol->dom->step)) sample_now = 1;
    else
    {
      double interval = bcd->length == 0 ? bcd->sampling[0] :
        (sol->output_interval > 0.0 ? sol->output_interval : sol->dom->step);

      if (fabs(bcd->latest + interval - sol->dom->time) <= sol->dom->step)
      {
	bcd->latest += interval;
	sample_now = 1;
      }
    }

    if (sample_now)
    {
      for (SET *item = SET_First (bcd->subset); item; item = SET_Next (item))
      {
	BODY *bod = MAP_Find (sol->dom->idb, item->data, NULL);
	{
	  if (bod)
	  {
	    double *sample;
	    ERRMEM (sample = malloc (bcd->size * sizeof (double)));
	    FEM_Mesh_Corotated_Conf (bod, sample);
	    SET_Insert (NULL, &bcd->pending, sample, NULL);
	  }
	}
      }
    }
  }
}

/* output pending samples */
void BCD_Append_Output (BCD *bcd)
{
  for (; bcd; bcd = bcd->next)
  {
#if MPI
    int rank, size, sendcount, *recvcounts, *displs, i, recvsize, j;
    double *sendbuf, *recvbuf, *ptr, *buf;
    SET *item;

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    MPI_Comm_size (MPI_COMM_WORLD, &size);

    sendcount = SET_Size (bcd->pending) * bcd->size;

    if (rank == 0)
    {
      ERRMEM (recvcounts = malloc(size * sizeof(int)));
      ERRMEM (displs = malloc(size * sizeof(int)));
    }

    MPI_Gather (&sendcount, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    for (displs[0] = recvsize = i = 0; i < size; i ++)
    {
      if (i) displs[i] = recvsize;
      recvsize += recvcounts[i];
    }

    ERRMEM (sendbuf = malloc(sendcount * sizeof(double)));

    if (rank == 0)
    {
      ERRMEM (recvbuf = malloc(recvsize * sizeof(double)));
    }

    for (item = SET_First (bcd->pending), ptr = sendbuf; item; item = SET_Next (item), ptr += bcd->size)
    {
      memcpy (ptr, item->data, bcd->size * sizeof(double));
    }

    MPI_Gatherv (sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free (sendbuf);

    if (rank == 0)
    {
      for (i = 0; i < size; i ++)
      {
	for (j = 0; j < recvcounts[i]; j ++)
	{
	  ptr = &recvbuf[displs[i]+j];
	  ERRMEM (buf = malloc (bcd->size * sizeof(double)));
	  memcpy (buf, ptr, bcd->size * sizeof(double));
	  SET_Insert (NULL, &bcd->pending, (void*) buf, NULL);
	}
      }

      free (recvbuf);
      free (displs);
      free (recvcounts);
    }
#endif

    for (SET *item = SET_First (bcd->pending); item; item = SET_Next (item))
    {
#if MPI
      if (rank == 0) /* append on rank 0 */
#endif
      lng_append_list_with_list_of_doubles (bcd->output, item->data, bcd->size);

      free (item->data); /* free for all ranks */
    }

    SET_Free (NULL, &bcd->pending); /* free for all ranks */
  }
}

/* destroy sampling list */
void BCD_Destroy (BCD *bcd)
{
  BCD *next;

  for (; bcd; bcd = next)
  {
    next = bcd->next;
    SET_Free (NULL, &bcd->subset);
    free (bcd->sampling);
    lng_xdecref (bcd->output);
#if MPI
    for (SET *item = bcd->pending; item; item = SET_Next (item))
    {
      free (item->data);
    }
    SET_Free (NULL, &bcd->pending);
#endif
    free (bcd);
  }
}
