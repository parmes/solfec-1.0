/*
 * com.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * parallel communication
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
#include <string.h>
#include <mpi.h>
#include "mem.h"
#include "com.h"
#include "map.h"
#include "alg.h"
#include "put.h"
#include "err.h"

typedef struct compattern COMPATTERN;

/* integer and doubles point to point
 * communication pattern */
struct compattern
{
  MPI_Comm comm; 

  int tag,
     *rankmap, /* maps ranks to send buffers indices */
      send_size,
    (*send_sizes) [3],
     *send_position,
     *send_rank,
      send_count,
     *recv_rank,
    (*recv_sizes) [3],
      recv_count;

  char **send_data,
       **recv_data;

  MPI_Request *recv_req; /* for receiving */
  MPI_Status *recv_sta; /* for receiving */

  MPI_Request *send_req; /* for sending */
  MPI_Status *send_sta; /* for sending */

  COMDATA *send,
	  *recv;

  int nsend,
      nrecv;
};

typedef struct comallpattern COMALLPATTERN;

/* integer and doubles all to all
 * communication pattern */
struct comallpattern
{
  MPI_Comm comm; 

  int ncpu,
    (*send_sizes) [3],
     *send_counts,
     *send_disps,
     *send_position,
      send_size,
    (*recv_sizes) [3],
     *recv_counts,
     *recv_disps,
     *recv_position,
      recv_size;

  char *send_data,
       *recv_data;

  COMDATA *send,
	  *recv;

  int nsend,
      nrecv;
};

#if PARDEBUG
static int* next_int (COMDATA *data, int size, int *idx)
{
  int *out = NULL;

  while (idx [0] < size)
  {
    if (idx [1] < data [idx [0]].ints)
    {
      out = &data [idx [0]].i [idx [1]];
      idx [1] ++;
      break;
    }
    else idx [0] ++, idx [1] = 0;
  }

  return out;
}

static double* next_double (COMDATA *data, int size, int *idx)
{
  double *out = NULL;

  while (idx [0] < size)
  {
    if (idx [1] < data [idx [0]].doubles)
    {
      out = &data [idx [0]].d [idx [1]];
      idx [1] ++;
      break;
    }
    else idx [0] ++, idx [1] = 0;
  }

  return out;
}
#endif

/* communicate integers and doubles using point to point communication */
uint64_t COM (MPI_Comm comm, int tag,
         COMDATA *send, int nsend,
	 COMDATA **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMDATA *cd;
  int rank,
      ncpu,
    (*send_sizes) [3],
     *send_position,
     *send_rank,
      send_count,
     *send_rank_all,
     *send_count_all,
     *send_rank_disp,
     *recv_rank,
    (*recv_sizes) [3],
      recv_count,
      i, j, k, l;
  uint64_t send_size;
  char **send_data,
       **recv_data;
  MPI_Request *req;
  MPI_Status *sta;
  void *p;

  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &ncpu);

  ERRMEM (send_sizes = MEM_CALLOC (ncpu * sizeof (int [3])));
  ERRMEM (send_position = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (send_rank = malloc (ncpu * sizeof (int)));
  ERRMEM (send_data = malloc (ncpu * sizeof (char*)));

  /* compute send sizes */
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    send_sizes [cd->rank][0] += cd->ints;
    send_sizes [cd->rank][1] += cd->doubles;
    MPI_Pack_size (cd->ints, MPI_INT, comm, &j);
    MPI_Pack_size (cd->doubles, MPI_DOUBLE, comm, &k);
    send_sizes [cd->rank][2] += (j + k);
  }

  /* allocate send buffers */
  for (send_size = i = 0; i < ncpu; i ++)
  {
    if (send_sizes [i][2])
    {
      ERRMEM (send_data [i] = malloc (send_sizes [i][2]));
      send_position [i] = 0;
      send_size += send_sizes [i][2];
    }
  }

  /* pack ints */
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    if (cd->ints)
    {
      MPI_Pack (cd->i, cd->ints, MPI_INT, send_data [cd->rank], send_sizes [cd->rank][2], &send_position [cd->rank], comm);
    }
  }

  /* pack doubles */
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    if (cd->doubles)
    {
      MPI_Pack (cd->d, cd->doubles, MPI_DOUBLE, send_data [cd->rank], send_sizes [cd->rank][2], &send_position [cd->rank], comm);
    }
  }

#if DEBUG
  for (i = 0; i < ncpu; i ++)
  {
    ASSERT_DEBUG (send_position [i] <= send_sizes [i][2], "Incorrect packing");
  }
#endif

  /* compute send ranks and move data */
  for (send_count = i = 0; i < ncpu; i ++)
  {
    if (send_sizes [i][2])
    {
      send_rank [send_count] = i;
      send_data [send_count] = send_data [i];
      send_sizes [send_count][0] = send_sizes [i][0];
      send_sizes [send_count][1] = send_sizes [i][1];
      send_sizes [send_count][2] = send_sizes [i][2];
      send_count ++;
    }
  }

  ERRMEM (send_count_all = malloc (ncpu * sizeof (int)));
  ERRMEM (recv_rank = malloc (ncpu * sizeof (int)));

  /* gather all send ranks */
  MPI_Allgather (&send_count, 1, MPI_INT, send_count_all, 1, MPI_INT, comm);
  ERRMEM (send_rank_disp = malloc (ncpu * sizeof (int)));
  for (send_rank_disp [0] = l = i = 0; i < ncpu; i ++)
  { l += send_count_all [i]; if (i < ncpu-1) send_rank_disp [i+1] = l; }
  ERRMEM (send_rank_all = malloc (l * sizeof (int)));
  MPI_Allgatherv (send_rank, send_count, MPI_INT, send_rank_all, send_count_all, send_rank_disp, MPI_INT, comm);

  /* compute receive ranks */
  for (recv_count = k = i = 0; i < l; i += send_count_all [k], k ++)
  {
    for (j = 0; j < send_count_all [k]; j ++)
    {
      if (send_rank_all [i+j] == rank) /* 'k'th rank is sending here */
      {
	recv_rank [recv_count] = k;
	recv_count ++;
	break;
      }
    }
  }

  ERRMEM (recv_sizes = malloc (recv_count * sizeof (int [3])));
  ERRMEM (req = malloc (recv_count * sizeof (MPI_Request)));
  ERRMEM (sta = malloc (recv_count * sizeof (MPI_Status)));

  /* communicate receive sizes */
  for (i = 0; i < recv_count; i ++)
  {
    MPI_Irecv (recv_sizes [i], 3, MPI_INT, recv_rank [i], tag, comm, &req [i]);
  }
  MPI_Barrier (comm);
  for (i = 0; i < send_count; i ++)
  {
    MPI_Rsend (send_sizes [i], 3, MPI_INT, send_rank [i], tag, comm);
  }
  MPI_Waitall (recv_count, req, sta);

  /* contiguous receive size */
  uint64_t size = recv_count * sizeof (COMDATA);
  for (i = 0; i < recv_count; i ++)
  {
    size += recv_sizes [i][0] * sizeof (int) + 
            recv_sizes [i][1] * sizeof (double);
  }

  /* prepare receive buffers */
  ERRMEM (recv_data = malloc (recv_count * sizeof (char*)));
  ERRMEM ((*recv) = malloc (size));
  p = (*recv) + recv_count;
  *nrecv = recv_count;
  for (i = 0, cd = *recv; i < recv_count; i ++, cd ++)
  {
    cd->rank = recv_rank [i];
    cd->ints = recv_sizes [i][0];
    cd->doubles = recv_sizes [i][1];
    cd->i = p; p = (cd->i + cd->ints);
    cd->d = p; p = (cd->d + cd->doubles);
    ERRMEM (recv_data [i] = malloc (recv_sizes [i][2]));
  }

  /* communicate data */
  for (i = 0; i < recv_count; i ++)
  {
    MPI_Irecv (recv_data [i], recv_sizes [i][2], MPI_PACKED, recv_rank [i], tag, comm, &req [i]);
  }
  MPI_Barrier (comm);
  for (i = 0; i < send_count; i ++)
  {
    MPI_Rsend (send_data [i], send_sizes [i][2], MPI_PACKED, send_rank [i], tag, comm);
  }
  MPI_Waitall (recv_count, req, sta);

  /* unpack data */
  for (i = j = 0; i < recv_count; i ++, j = 0)
  {
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, (*recv) [i].i, (*recv) [i].ints, MPI_INT, comm);
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, (*recv) [i].d, (*recv) [i].doubles, MPI_DOUBLE, comm);
  }

  /* cleanup */
  free (send_rank_disp);
  free (send_sizes);
  free (send_position);
  free (send_rank);
  for (i = 0; i < send_count; i ++)
    free (send_data [i]);
  free (send_data);
  free (send_count_all);
  free (send_rank_all);
  free (recv_rank);
  free (recv_sizes);
  for (i = 0; i < recv_count; i ++)
    free (recv_data [i]);
  free (recv_data);
  free (req);
  free (sta);

  return send_size;
}

#if 0
/* simple version of the above */
uint64_t COM_simple (MPI_Comm comm, int tag,
         COMDATA *send, int nsend,
	 COMDATA **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMDATA *cd;
  int rank,
      ncpu,
      nsend_max,
     *ranks_ints_doubles,
      i, j, l, m, n;
  uint64_t send_size;
  void *p;

  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &ncpu);

  nsend_max = PUT_int_max (nsend);

  /* (ncpu + 1) since local send buffer is stored at the end */
  ERRMEM (ranks_ints_doubles = MEM_CALLOC (3 * nsend_max * (ncpu + 1) * sizeof (int)));

  /* compute send sizes */
  for (send_size = i = j = 0, l =  3 * nsend_max * ncpu, cd = send; i < nsend; i ++, cd ++)
  {
    if (cd->ints || cd->doubles)
    {
      send_size += cd->ints * sizeof(int) + cd->doubles * sizeof(double);
      ranks_ints_doubles [l+j+0] = cd->rank;
      ranks_ints_doubles [l+j+1] = cd->ints;
      ranks_ints_doubles [l+j+2] = cd->doubles;
      j += 3;
    }
  }

  /* gather send ranks, ints, doubles */
  MPI_Allgather (&ranks_ints_doubles[3 * nsend_max * ncpu], 3 * nsend_max, MPI_INT, ranks_ints_doubles, 3 * nsend_max, MPI_INT, comm);

  /* calculate contiguous receive size and prepare receive buffers */
  (*nrecv) = m = n = 0;
  for (i = 0; i < ncpu; i ++)
  {
    for (j = 0; j < 3 * nsend_max; j += 3)
    {
      l = 3 * nsend_max * i + j;

      if(ranks_ints_doubles[l+1] || ranks_ints_doubles[l+2])
      {
        if (ranks_ints_doubles[l+0] == rank)
	{
	  (*nrecv) ++;
	  m += ranks_ints_doubles[l+1];
	  n += ranks_ints_doubles[l+2];
	}
      }
      else break;
    }
  }
  uint64_t size = (*nrecv) * sizeof (COMDATA) + m * sizeof (int) + n * sizeof (double);
  ERRMEM ((*recv) = malloc (size));
  p = (*recv) + (*nrecv);
  cd = *recv;
  for (i = 0; i < ncpu; i ++)
  {
    for (j = 0; j < 3 * nsend_max; j += 3)
    {
      l = 3 * nsend_max * i + j;

      if (ranks_ints_doubles[l+1] || ranks_ints_doubles[l+2])
      {
        if (ranks_ints_doubles[l+0] == rank)
	{
	  cd->rank = i;
	  cd->ints = ranks_ints_doubles[l+1];
	  cd->doubles = ranks_ints_doubles[l+2];
	  cd->i = p; p = (cd->i + cd->ints);
	  cd->d = p; p = (cd->d + cd->doubles);
	  cd ++;
	}
      }
      else break;
    }
  }
  free (ranks_ints_doubles);

  /* communicate data */
  int recv_count = 0;
  MPI_Request *req;
  MPI_Status *sta;
  ERRMEM (req = malloc (2 * (*nrecv) * sizeof (MPI_Request)));
  ERRMEM (sta = malloc (2 * (*nrecv) * sizeof (MPI_Status)));
  for (i = 0, cd = *recv; i < *nrecv; i ++, cd ++)
  {
    if (cd->ints)
    {
      MPI_Irecv (cd->i, cd->ints, MPI_INT, cd->rank, tag, comm, &req[recv_count]);
      recv_count ++;
    }
    if (cd->doubles)
    {
      MPI_Irecv (cd->d, cd->doubles, MPI_DOUBLE, cd->rank, tag, comm, &req[recv_count]);
      recv_count ++;
    }
  }
  MPI_Barrier (comm);
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    if (cd->ints) 
    {
      MPI_Rsend (cd->i, cd->ints, MPI_INT, cd->rank, tag, comm);
    }
    if (cd->doubles) 
    {
      MPI_Rsend (cd->d, cd->doubles, MPI_DOUBLE, cd->rank, tag, comm);
    }
  }
  MPI_Waitall (recv_count, req, sta);
  free (req);
  free (sta);

  return send_size;
}
#endif

/* communicate integers and doubles using all to all communication */
uint64_t COMALL (MPI_Comm comm,
            COMDATA *send, int nsend,
	    COMDATA **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMDATA *cd;
  uint64_t send_size64,
          *send_disps64,
	   recv_size64,
	  *recv_disps64;
  int rank,
      ncpu,
    (*send_sizes) [3],
     *send_counts,
     *send_position,
    (*recv_sizes) [3],
     *recv_counts,
     *recv_position,
      i, j, k;
  char *send_data,
       *recv_data;
  void *p;

  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &ncpu);

  ERRMEM (send_sizes = MEM_CALLOC (ncpu * sizeof (int [3])));
  ERRMEM (send_counts = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (send_disps64 = MEM_CALLOC (ncpu * sizeof (uint64_t)));
  ERRMEM (send_position = MEM_CALLOC (ncpu * sizeof (int)));

  /* compute send sizes */
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    send_sizes [cd->rank][0] += cd->ints;
    send_sizes [cd->rank][1] += cd->doubles;
    MPI_Pack_size (cd->ints, MPI_INT, comm, &j);
    MPI_Pack_size (cd->doubles, MPI_DOUBLE, comm, &k);
    send_sizes [cd->rank][2] += (j + k);
  }

  /* compute send displacements */
  for (send_size64 = i = 0; i < ncpu; i ++)
  {
    send_counts [i] = send_sizes [i][2];
    send_size64 += send_counts [i]; /* note that this counter would overflow at 2.14GB if it was integer */
    if (i < (ncpu - 1)) send_disps64 [i+1] = send_size64;
  }
  send_disps64 [0] = 0;

  /* allocate send buffer */
  ERRMEM (send_data = malloc (send_size64));

  /* pack ints */
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    if (cd->ints)
    {
      MPI_Pack (cd->i, cd->ints, MPI_INT, &send_data [send_disps64 [cd->rank]], send_counts [cd->rank], &send_position [cd->rank], comm);
    }
  }

  /* pack doubles */
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    if (cd->doubles)
    {
      MPI_Pack (cd->d, cd->doubles, MPI_DOUBLE, &send_data [send_disps64 [cd->rank]], send_counts [cd->rank], &send_position [cd->rank], comm); 
    }
  }

#if DEBUG
  for (i = 0; i < ncpu; i ++)
  {
    ASSERT_DEBUG (send_position [i] <= send_counts [i], "Incorrect packing");
  }
#endif

  ERRMEM (recv_sizes = MEM_CALLOC (ncpu * sizeof (int [3])));
  ERRMEM (recv_counts = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (recv_disps64 = MEM_CALLOC (ncpu * sizeof (uint64_t)));
  ERRMEM (recv_position = MEM_CALLOC (ncpu * sizeof (int)));

  /* distribute send sizes into receive sizes */
  MPI_Alltoall (send_sizes, 3, MPI_INT, recv_sizes, 3, MPI_INT, comm);

  /* compute receive displacements */
  for (recv_size64 = i = 0; i < ncpu; i ++)
  {
    recv_counts [i] = recv_sizes [i][2];
    recv_size64 += recv_counts [i];
    if (i < (ncpu - 1)) recv_disps64 [i+1] = recv_size64;
  }
  recv_disps64 [0] = 0;

  /* allocate receive buffer */
  ERRMEM (recv_data = malloc (recv_size64));

  uint64_t sr64[2] = {send_size64, recv_size64}, sr64mx[2];

  MPI_Allreduce (sr64, sr64mx, 2, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);

  if (sr64mx[0] < INT_MAX && sr64mx[1] < INT_MAX)
  {
    int *send_disps,
        *recv_disps;

    ERRMEM (send_disps = MEM_CALLOC (ncpu * sizeof (int)));
    ERRMEM (recv_disps = MEM_CALLOC (ncpu * sizeof (int)));

    for (i = 0; i < ncpu; i ++)
    {
      send_disps[i] = send_disps64[i];
      recv_disps[i] = recv_disps64[i];
    }

    /* all to all send and receive */
    MPI_Alltoallv (send_data, send_counts, send_disps, MPI_PACKED, recv_data, recv_counts, recv_disps, MPI_PACKED, comm);

    free (send_disps);
    free (recv_disps);
  }
  else
  {
    WARNING (rank != 0, "COMALL MPI_Alltoallv displacement size > INT_MAX resulting in point-to-point workaround");
    int recv_count = 0;
    MPI_Request *req;
    MPI_Status *sta;
    ERRMEM (req = malloc (ncpu * sizeof (MPI_Request)));
    ERRMEM (sta = malloc (ncpu * sizeof (MPI_Status)));
    for (i = 0; i < ncpu; i ++)
    {
      if (recv_counts[i])
      {
	MPI_Irecv (&recv_data[recv_disps64[i]], recv_counts[i], MPI_PACKED, i, 0, comm, &req[recv_count]);
	recv_count ++;
      }
    }
    MPI_Barrier (comm);
    for (i = 0; i < ncpu; i ++)
    {
      if (send_counts[i])
      {
	MPI_Rsend (&send_data[send_disps64[i]], send_counts[i], MPI_PACKED, i, 0, comm);
      }
    }
    MPI_Waitall (recv_count, req, sta);
    free (req);
    free (sta);
  }

  if (recv_size64)
  {
    /* contiguous receive size */
    uint64_t size = ncpu * sizeof (COMDATA);
    for (i = 0; i < ncpu; i ++)
    {
      size += recv_sizes [i][0] * sizeof (int) + 
	      recv_sizes [i][1] * sizeof (double);
    }

    /* prepare output receive data */
    ERRMEM ((*recv) = malloc (size));
    p = (*recv) + ncpu;
    for (i = 0, cd = *recv; i < ncpu; i ++, cd ++)
    {
      cd->rank = i;
      cd->ints = recv_sizes [i][0];
      cd->doubles = recv_sizes [i][1];
      cd->i = p; p = (cd->i + cd->ints);
      cd->d = p; p = (cd->d + cd->doubles);
    }

    /* unpack data */
    for (i = 0; i < ncpu; i ++)
    {
      if (recv_counts [i])
      {
	MPI_Unpack (&recv_data [recv_disps64 [i]], recv_counts [i], &recv_position [i], (*recv) [i].i, (*recv) [i].ints, MPI_INT, comm);
	MPI_Unpack (&recv_data [recv_disps64 [i]], recv_counts [i], &recv_position [i], (*recv) [i].d, (*recv) [i].doubles, MPI_DOUBLE, comm);
      }
    }

    /* compress receive storage */
    for (*nrecv = i = 0; i < ncpu; i ++)
    {
      if (recv_counts [i])
      {
	(*recv) [*nrecv] = (*recv) [i];
	(*nrecv) ++;
      }
    }
  }
  else
  {
    *recv = NULL;
    *nrecv = 0;
  }

  /* cleanup */
  free (send_sizes);
  free (send_counts);
  free (send_disps64);
  free (send_position);
  free (recv_sizes);
  free (recv_counts);
  free (recv_disps64);
  free (recv_position);
  free (send_data);
  free (recv_data);

  return send_size64;
}

/* alternative approach to communicating integers and doubles using all to all communication;
   more appropriate for sending large data volumes */
uint64_t COMALL2 (MPI_Comm comm,
            COMDATA *send, int nsend,
	    COMDATA **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMDATA *cd;
  uint64_t *int_send_counts,
	   int_send_size,
	   *int_send_counts_max,
	   *int_block_position,
	   *double_send_counts,
	   double_send_size,
	   *double_send_counts_max,
	   *double_block_position;
  int rank, ncpu, i, j,
	   *int_block_send_counts,
	   *int_block_send_disps,
	   int_block_send_size,
	   *double_block_send_counts,
	   *double_block_send_disps,
	   double_block_send_size;
  int *int_send_data,
      *int_recv_data;
  double *double_send_data,
         *double_recv_data;
  void *p;

  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &ncpu);

  ERRMEM (int_send_counts = MEM_CALLOC (ncpu * sizeof (uint64_t)));
  ERRMEM (int_send_counts_max = MEM_CALLOC (ncpu * sizeof (uint64_t)));
  ERRMEM (int_block_position = MEM_CALLOC (ncpu * sizeof (uint64_t)));
  ERRMEM (double_send_counts = MEM_CALLOC (ncpu * sizeof (uint64_t)));
  ERRMEM (double_send_counts_max = MEM_CALLOC (ncpu * sizeof (uint64_t)));
  ERRMEM (double_block_position = MEM_CALLOC (ncpu * sizeof (uint64_t)));

  ERRMEM (int_block_send_counts = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (int_block_send_disps = MEM_CALLOC ((ncpu+1) * sizeof (int)));
  ERRMEM (double_block_send_counts = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (double_block_send_disps = MEM_CALLOC ((ncpu+1) * sizeof (int)));

  /* compute send counts */
  for (int_send_size = double_send_size = i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    int_send_counts [cd->rank] += cd->ints;
    double_send_counts [cd->rank] += cd->doubles;
    int_send_size += cd->ints;
    double_send_size += cd->doubles;
  }

  /* reduce total send sizes */
  uint64_t ids[2] = {int_send_size, double_send_size}, idsx[2];

  MPI_Allreduce (ids, idsx, 2, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);

  /* block sizes that allow to use integer displacements in MPI_Alltoallv */
  int int_block_size = (1 + idsx[0]/INT_MAX),
      double_block_size = (1 + idsx[1]/INT_MAX);

  /* contiguous block types */
  MPI_Datatype int_block_type, double_block_type;
  MPI_Type_contiguous(int_block_size, MPI_INT, &int_block_type);
  MPI_Type_contiguous(double_block_size, MPI_DOUBLE, &double_block_type);
  MPI_Type_commit(&int_block_type);
  MPI_Type_commit(&double_block_type);

  /* compute send displacements using these blocks */
  for (int_block_send_size = double_block_send_size = i = 0; i < ncpu; i ++)
  {
    int_block_send_counts [i] = 1 + int_send_counts [i] / int_block_size;
    int_block_send_size += int_block_send_counts [i];
    int_block_send_disps [i+1] = int_block_send_size;

    double_block_send_counts [i] = 1 + double_send_counts [i] / double_block_size;
    double_block_send_size += double_block_send_counts [i];
    double_block_send_disps [i+1] = double_block_send_size;
  }

  /* allocate send buffers */
  ERRMEM (int_send_data = malloc (int_block_send_size * sizeof(int [int_block_size])));
  ERRMEM (double_send_data = malloc (double_block_send_size * sizeof(double [double_block_size])));

  /* pack send data into send buffers */
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    j = cd->rank;
    memcpy (&int_send_data[int_block_size*int_block_send_disps[j] + int_block_position[j]], cd->i, sizeof(int [cd->ints]));
    int_block_position[j] += cd->ints;
    memcpy (&double_send_data[double_block_size*double_block_send_disps[j] + double_block_position[j]], cd->d, sizeof(double [cd->doubles]));
    double_block_position[j] += cd->doubles;
  }

  uint64_t *int_recv_counts,
	   *double_recv_counts;

  int      *int_block_recv_counts,
	   *int_block_recv_disps,
	   int_block_recv_size,
	   *double_block_recv_counts,
	   *double_block_recv_disps,
	   double_block_recv_size;

  /* allocate receive indexing buffers */
  ERRMEM (int_recv_counts = MEM_CALLOC (ncpu * sizeof (uint64_t)));
  ERRMEM (double_recv_counts = MEM_CALLOC (ncpu * sizeof (uint64_t)));

  ERRMEM (int_block_recv_counts = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (int_block_recv_disps = MEM_CALLOC ((ncpu+1) * sizeof (int)));
  ERRMEM (double_block_recv_counts = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (double_block_recv_disps = MEM_CALLOC ((ncpu+1) * sizeof (int)));

  /* distribute send sizes into receive sizes */
  MPI_Alltoall (int_send_counts, 1, MPI_UINT64_T, int_recv_counts, 1, MPI_UINT64_T, comm);
  MPI_Alltoall (double_send_counts, 1, MPI_UINT64_T, double_recv_counts, 1, MPI_UINT64_T, comm);

  MPI_Alltoall (int_block_send_counts, 1, MPI_INT, int_block_recv_counts, 1, MPI_INT, comm);
  MPI_Alltoall (double_block_send_counts, 1, MPI_INT, double_block_recv_counts, 1, MPI_INT, comm);

  /* compute receive block displacements */
  for (int_block_recv_size = double_block_recv_size = i = 0; i < ncpu; i ++)
  {
    int_block_recv_size += int_block_recv_counts [i];
    int_block_recv_disps [i+1] = int_block_recv_size;
    double_block_recv_size += double_block_recv_counts [i];
    double_block_recv_disps [i+1] = double_block_recv_size;
  }

  /* allocate receive buffers */
  ERRMEM (int_recv_data = malloc (int_block_recv_size * sizeof(int [int_block_size])));
  ERRMEM (double_recv_data = malloc (double_block_recv_size * sizeof(double [double_block_size])));

  /* all to all send and receive */
  MPI_Alltoallv (int_send_data, int_block_send_counts, int_block_send_disps, int_block_type,
                 int_recv_data, int_block_recv_counts, int_block_recv_disps, int_block_type, comm);
  MPI_Alltoallv (double_send_data, double_block_send_counts, double_block_send_disps, double_block_type,
                 double_recv_data, double_block_recv_counts, double_block_recv_disps, double_block_type, comm);

  if (int_block_recv_size || double_block_recv_size)
  {
    /* contiguous receive size */
    uint64_t size = ncpu * sizeof (COMDATA);
    for (i = 0; i < ncpu; i ++)
    {
      size += int_block_recv_counts [i] * sizeof (int [int_block_size]) + 
	      double_block_recv_counts [i] * sizeof (double [double_block_size]);
    }

    /* prepare output receive data */
    ERRMEM ((*recv) = malloc (size));
    p = (*recv) + ncpu;
    for (i = 0, cd = *recv; i < ncpu; i ++, cd ++)
    {
      cd->rank = i;
      cd->ints = int_recv_counts [i];
      cd->doubles = double_recv_counts [i];
      cd->i = p; p = (cd->i + cd->ints);
      cd->d = p; p = (cd->d + cd->doubles);
    }

    /* unpack data */
    for (i = 0; i < ncpu; i ++)
    {
      memcpy ((*recv)[i].i, &int_recv_data[int_block_size*int_block_recv_disps[i]], (*recv) [i].ints * sizeof(int));
      memcpy ((*recv)[i].d, &double_recv_data[double_block_size*double_block_recv_disps[i]], (*recv) [i].doubles * sizeof(double));
    }

    /* compress receive storage */
    for (*nrecv = i = 0; i < ncpu; i ++)
    {
      if (int_recv_counts [i] || double_recv_counts[i])
      {
	(*recv) [*nrecv] = (*recv) [i];
	(*nrecv) ++;
      }
    }
  }
  else
  {
    *recv = NULL;
    *nrecv = 0;
  }

  /* cleanup */
  free (int_send_counts);
  free (int_send_counts_max);
  free (int_block_position);
  free (double_send_counts);
  free (double_send_counts_max);
  free (double_block_position);
  free (int_block_send_counts);
  free (int_block_send_disps);
  free (double_block_send_counts);
  free (double_block_send_disps);
  free (int_send_data);
  free (double_send_data);
  free (int_recv_counts);
  free (double_recv_counts);
  free (int_block_recv_counts);
  free (int_block_recv_disps);
  free (double_block_recv_counts);
  free (double_block_recv_disps);
  free (int_recv_data);
  free (double_recv_data);

  return int_block_send_size * sizeof(int [int_block_size]) + double_block_send_size * sizeof(double [double_block_size]);
}

#define COMALLX COMALL2 /* COMALL or COMALL2 */

/* communicate one set of integers and doubles to all other processors */
uint64_t COMONEALL (MPI_Comm comm, COMDATA send,
	       COMDATA **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMDATA *send_data;
  int i, j, rank, ncpu, ret;

  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &ncpu);

  ERRMEM (send_data = MEM_CALLOC ((ncpu - 1) * sizeof (COMDATA)));

  for (i = 0; i < ncpu; i ++)
  {
    if (i != rank)
    {
      send_data [j] = send;
      send_data [j].rank = i;
      j ++;
    }
  }

  ret = COMALLX (comm, send_data, ncpu - 1, recv, nrecv);

  free (send_data);

  return ret;
}

/* communicate objects using point to point communication */
uint64_t COMOBJS (MPI_Comm comm, int tag,
	     OBJ_Pack pack,
	     void *data,
	     OBJ_Unpack unpack,
             COMOBJ *send, int nsend,
	     COMOBJ **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMDATA *send_data,
	  *recv_data,
	  *cd, *cc;
  int recv_count,
      i, n;
  uint64_t ret;
  COMOBJ *co;
  MAP *map;
  MEM mem;

  ERRMEM (send_data = malloc (nsend * sizeof (COMDATA)));
  MEM_Init (&mem, sizeof (MAP), MAX (nsend, 64));

  /* pack objects */
  for (i = 0, cd = send_data, co = send, map = NULL; i < nsend; i ++, cd ++, co ++)
  {
    int isize = 0,
	dsize = 0;

    cd->rank = co->rank;

    if ((cc = MAP_Find (map, co->o, NULL))) /* same object was already packed */
    {
      cd->ints = cc->ints;
      cd->doubles = cc->doubles;
      cd->i = cc->i;
      cd->d = cc->d;
    }
    else
    {
      cd->ints = 0;
      cd->doubles = 0;
      cd->i = NULL;
      cd->d = NULL;

      pack (co->o, &dsize, &cd->d, &cd->doubles, &isize, &cd->i, &cd->ints);

      MAP_Insert (&mem, &map, co->o, cd, NULL);
    }
  }

  /* send and receive packed data */
  if (tag == INT_MIN) ret = COMALLX (comm, send_data, nsend, &recv_data, &recv_count); /* all to all */
  else ret = COM (comm, tag, send_data, nsend, &recv_data, &recv_count); /* point to point */

#if PARDEBUG
  {
    int debug_send_count, *ip, *jp, ii [2], jj [2];
    COMDATA *debug_send_data;
    double *qq, *pp;

    /* send backwards */
    if (tag == INT_MIN) COMALLX (comm, recv_data, recv_count, &debug_send_data, &debug_send_count);
    else COM (comm, tag, recv_data, recv_count, &debug_send_data, &debug_send_count);

    ii[0] = 0, ii[1] = 0;
    jj[0] = 0, jj[1] = 0;
    do
    {
      ip = next_int (send_data, nsend, ii);
      jp = next_int (debug_send_data, debug_send_count, jj);
      if (ip && jp)
      {
	ASSERT_DEBUG (*ip == *jp, "Integer values mismatch");
      }
      else
      {
	ASSERT_DEBUG (!ip && !jp, "Integer count mismatch");
      }
    }
    while (ip && jp);

    ii[0] = 0, ii[1] = 0;
    jj[0] = 0, jj[1] = 0;
    do
    {
      qq = next_double (send_data, nsend, ii);
      pp = next_double (debug_send_data, debug_send_count, jj);
      if (qq && pp)
      {
	ASSERT_DEBUG (*qq == *pp, "Double values mismatch");
      }
      else
      {
	ASSERT_DEBUG (!qq && !pp, "Double count mismatch");
      }
    }
    while (qq && pp);

    free (debug_send_data);
  }
#endif

  if (recv_count)
  {
    *nrecv = recv_count;
    ERRMEM (*recv = malloc ((*nrecv) * sizeof (COMOBJ)));

    /* unpack received objects */
    for (n = i = 0, cd = recv_data, co = *recv; i < recv_count; i ++, cd ++)
    {
      int ipos = 0,
	  dpos = 0;

      do
      {
	if (n == *nrecv)
	{
	  *nrecv *= 2; /* resize the receive buffer */
	  ERRMEM (*recv = realloc (*recv, (*nrecv) * sizeof (COMOBJ)));
	  co = *recv + n; /* and reset the current pointer */
	}

	co->rank = cd->rank;
	co->o = unpack (data, &dpos, cd->d, cd->doubles, &ipos, cd->i, cd->ints);

	co ++;
	n ++;

      } while (ipos < cd->ints || dpos < cd->doubles); /* while something is left to unpack */
    }

    /* truncate output */
    if (n) ERRMEM (*recv = realloc (*recv, n * sizeof (COMOBJ)));
    *nrecv = n;
  }
  else
  {
    *recv = NULL;
    *nrecv = 0;
  }

  /* cleanup */
  for (MAP *item = MAP_First (map); item; item = MAP_Next (item))
  { cd = item->data; free (cd->i); free (cd->d); }
  MEM_Release (&mem);
  free (send_data);
  free (recv_data); /* contiguous */

  return ret;
}

/* communicate objects using all to all communication */
uint64_t COMOBJSALL (MPI_Comm comm,
	        OBJ_Pack pack,
	        void *data,
	        OBJ_Unpack unpack,
                COMOBJ *send, int nsend,
	        COMOBJ **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  return COMOBJS (comm, INT_MIN, pack, data, unpack, send, nsend, recv, nrecv); /* this is only a wrapper */
}

/* communicate an object to all other processors */
uint64_t COMOBJALL (MPI_Comm comm,
	       OBJ_Pack pack,
	       void *data,
	       OBJ_Unpack unpack,
	       void *object,
	       COMOBJ **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMOBJ *send_data;
  int i, j, rank, ncpu, ret;

  if (object)
  {
    MPI_Comm_rank (comm, &rank);
    MPI_Comm_size (comm, &ncpu);

    ERRMEM (send_data = MEM_CALLOC ((ncpu - 1) * sizeof (COMOBJ)));

    for (i = j = 0; i < ncpu; i ++)
    {
      if (i != rank)
      {
	send_data [j].rank = i;
	send_data [j].o = object;
	j ++;
      }
    }

    ret = COMOBJSALL (comm, pack, data, unpack, send_data, ncpu - 1, recv, nrecv);

    free (send_data);
  }
  else ret = COMOBJSALL (comm, pack, data, unpack, NULL, 0, recv, nrecv);

  return ret;
}

/* create a repetitive point to point communication pattern;
 * ranks and sizes must not change during the communication;
 * pointers to send and receive buffers data must not change */
void* COM_Pattern (MPI_Comm comm, int tag,
                   COMDATA *send, int nsend,
	           COMDATA **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMPATTERN *pattern;
  COMDATA *cd;
  int rank,
      ncpu,
     *send_rank_all,
     *send_count_all,
     *send_rank_disp,
      i, j, k, l;
  void *p;

  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &ncpu);

  ERRMEM (pattern = malloc (sizeof (COMPATTERN)));
  ERRMEM (pattern->rankmap = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (pattern->send_sizes = MEM_CALLOC (ncpu * sizeof (int [3])));
  ERRMEM (pattern->send_position = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (pattern->send_rank = malloc (ncpu * sizeof (int)));
  ERRMEM (pattern->send_data = malloc (ncpu * sizeof (char*)));
  pattern->nsend = nsend;
  pattern->send = send;
  pattern->comm = comm;
  pattern->tag = tag;

  /* compute send sizes */
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    pattern->send_sizes [cd->rank][0] += cd->ints;
    pattern->send_sizes [cd->rank][1] += cd->doubles;
    MPI_Pack_size (cd->ints, MPI_INT, comm, &j);
    MPI_Pack_size (cd->doubles, MPI_DOUBLE, comm, &k);
    pattern->send_sizes [cd->rank][2] += (j + k);
  }

  /* allocate send buffers and prepare rank map */
  for (pattern->send_size = i = j = 0; i < ncpu; i ++)
  {
    if (pattern->send_sizes [i][2])
    {
      ERRMEM (pattern->send_data [i] = malloc (pattern->send_sizes [i][2]));
      pattern->rankmap [i] = j;
      pattern->send_size += pattern->send_sizes [i][2];
      j ++;
    }
  }

  /* compute send ranks and move data */
  for (pattern->send_count = i = 0; i < ncpu; i ++)
  {
    if (pattern->send_sizes [i][2])
    {
      pattern->send_rank [pattern->send_count] = i;
      pattern->send_data [pattern->send_count] = pattern->send_data [i];
      pattern->send_sizes [pattern->send_count][0] = pattern->send_sizes [i][0];
      pattern->send_sizes [pattern->send_count][1] = pattern->send_sizes [i][1];
      pattern->send_sizes [pattern->send_count][2] = pattern->send_sizes [i][2];
      pattern->send_count ++;
    }
  }

  ERRMEM (send_count_all = malloc (ncpu * sizeof (int)));
  ERRMEM (pattern->recv_rank = malloc (ncpu * sizeof (int)));

  /* gather all send ranks */
  MPI_Allgather (&pattern->send_count, 1, MPI_INT, send_count_all, 1, MPI_INT, comm);
  ERRMEM (send_rank_disp = malloc (ncpu * sizeof (int)));
  for (send_rank_disp [0] = l = i = 0; i < ncpu; i ++)
  { l += send_count_all [i]; if (i < ncpu-1) send_rank_disp [i+1] = l; }
  ERRMEM (send_rank_all = malloc (l * sizeof (int)));
  MPI_Allgatherv (pattern->send_rank, pattern->send_count, MPI_INT, send_rank_all, send_count_all, send_rank_disp, MPI_INT, comm);

  /* compute receive ranks */
  for (pattern->recv_count = k = i = 0; i < l; i += send_count_all [k], k ++)
  {
    for (j = 0; j < send_count_all [k]; j ++)
    {
      if (send_rank_all [i+j] == rank) /* 'k'th rank is sending here */
      {
	pattern->recv_rank [pattern->recv_count] = k;
	pattern->recv_count ++;
	break;
      }
    }
  }

  ERRMEM (pattern->recv_sizes = malloc (pattern->recv_count * sizeof (int [3])));
  ERRMEM (pattern->recv_req = malloc (pattern->recv_count * sizeof (MPI_Request)));
  ERRMEM (pattern->recv_sta = malloc (pattern->recv_count * sizeof (MPI_Status)));
  ERRMEM (pattern->send_req = malloc (pattern->send_count * sizeof (MPI_Request)));
  ERRMEM (pattern->send_sta = malloc (pattern->send_count * sizeof (MPI_Status)));

  /* communicate receive sizes */
  for (i = 0; i < pattern->recv_count; i ++)
  {
    MPI_Irecv (pattern->recv_sizes [i], 3, MPI_INT, pattern->recv_rank [i], tag, comm, &pattern->recv_req [i]);
  }
  MPI_Barrier (comm);
  for (i = 0; i < pattern->send_count; i ++)
  {
    MPI_Rsend (pattern->send_sizes [i], 3, MPI_INT, pattern->send_rank [i], tag, comm);
  }
  MPI_Waitall (pattern->recv_count, pattern->recv_req, pattern->recv_sta);

  /* contiguous receive size */
  uint64_t size = pattern->recv_count * sizeof (COMDATA);
  for (i = 0; i < pattern->recv_count; i ++)
  {
    size += pattern->recv_sizes [i][0] * sizeof (int) + 
            pattern->recv_sizes [i][1] * sizeof (double);
  }

  /* prepare receive buffers */
  ERRMEM (pattern->recv_data = malloc (pattern->recv_count * sizeof (char*)));
  ERRMEM (pattern->recv = malloc (size));
  p = pattern->recv + pattern->recv_count;
  pattern->nrecv = pattern->recv_count;
  for (i = 0, cd = pattern->recv; i < pattern->recv_count; i ++, cd ++)
  {
    cd->rank = pattern->recv_rank [i];
    cd->ints = pattern->recv_sizes [i][0];
    cd->doubles = pattern->recv_sizes [i][1];
    cd->i = p; p = (cd->i + cd->ints);
    cd->d = p; p = (cd->d + cd->doubles);
    ERRMEM (pattern->recv_data [i] = malloc (pattern->recv_sizes [i][2]));
  }

  /* truncate */
  if (pattern->send_count)
  {
    ERRMEM (pattern->send_sizes = realloc (pattern->send_sizes, pattern->send_count * sizeof (int [3])));
    ERRMEM (pattern->send_position = realloc (pattern->send_position, pattern->send_count * sizeof (int)));
    ERRMEM (pattern->send_rank = realloc (pattern->send_rank, pattern->send_count * sizeof (int)));
    ERRMEM (pattern->send_data = realloc (pattern->send_data, pattern->send_count * sizeof (char*)));
  }
  if (pattern->recv_count) ERRMEM (pattern->recv_rank = realloc (pattern->recv_rank, pattern->recv_count * sizeof (int)));
 
  /* cleanup */
  free (send_rank_disp);
  free (send_count_all);
  free (send_rank_all);

  /* output */
  *nrecv = pattern->nrecv;
  *recv = pattern->recv;

  return pattern;
}

/* communicate integers and doubles accodring
 * to the pattern computed by COM_Pattern */
uint64_t COM_Repeat (void *pattern)
{
  COMPATTERN *cp = pattern;
  int *rankmap = cp->rankmap,
      *send_position = cp->send_position,
     (*send_sizes) [3] = cp->send_sizes,
      *send_rank = cp->send_rank,
       send_count = cp->send_count,
      *recv_rank = cp->recv_rank,
     (*recv_sizes) [3] = cp->recv_sizes,
       recv_count = cp->recv_count,
       tag = cp->tag,
       nsend = cp->nsend,
       i, j;
  char **send_data = cp->send_data,
       **recv_data = cp->recv_data;
  MPI_Request *recv_req = cp->recv_req;
  MPI_Status *recv_sta = cp->recv_sta;
  MPI_Comm comm = cp->comm;
  COMDATA *cd;

  for (i = 0; i < send_count; i ++)
    send_position [i] = 0;

  /* pack ints */
  for (i = 0, cd = cp->send; i < nsend; i ++, cd ++)
  {
    if (cd->ints)
    {
      j = rankmap [cd->rank];
      MPI_Pack (cd->i, cd->ints, MPI_INT, send_data [j], send_sizes [j][2], &send_position [j], comm);
    }
  }

  /* pack doubles */
  for (i = 0, cd = cp->send; i < nsend; i ++, cd ++)
  {
    if (cd->doubles)
    {
      j = rankmap [cd->rank];
      MPI_Pack (cd->d, cd->doubles, MPI_DOUBLE, send_data [j], send_sizes [j][2], &send_position [j], comm);
    }
  }

  /* communicate data */
  for (i = 0; i < recv_count; i ++)
  {
    MPI_Irecv (recv_data [i], recv_sizes [i][2], MPI_PACKED, recv_rank [i], tag, comm, &recv_req [i]);
  }
  MPI_Barrier (comm);
  for (i = 0; i < send_count; i ++)
  {
    MPI_Rsend (send_data [i], send_sizes [i][2], MPI_PACKED, send_rank [i], tag, comm);
  }
  MPI_Waitall (recv_count, recv_req, recv_sta);

  /* unpack data */
  for (i = j = 0, cd = cp->recv; i < recv_count; i ++, cd ++, j = 0)
  {
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, cd->i, cd->ints, MPI_INT, comm);
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, cd->d, cd->doubles, MPI_DOUBLE, comm);
  }

  return cp->send_size;
}

/* non-blocking send */
uint64_t COM_Send (void *pattern)
{
  COMPATTERN *cp = pattern;
  int *rankmap = cp->rankmap,
      *send_position = cp->send_position,
     (*send_sizes) [3] = cp->send_sizes,
      *send_rank = cp->send_rank,
      *recv_rank = cp->recv_rank,
     (*recv_sizes) [3] = cp->recv_sizes,
       send_count = cp->send_count,
       recv_count = cp->recv_count,
       tag = cp->tag,
       nsend = cp->nsend,
       i, j;
  char **send_data = cp->send_data,
       **recv_data = cp->recv_data;
  MPI_Request *send_req = cp->send_req,
              *recv_req = cp->recv_req;
  MPI_Comm comm = cp->comm;
  COMDATA *cd;

  for (i = 0; i < send_count; i ++)
    send_position [i] = 0;

  /* pack ints */
  for (i = 0, cd = cp->send; i < nsend; i ++, cd ++)
  {
    if (cd->ints)
    {
      j = rankmap [cd->rank];
      MPI_Pack (cd->i, cd->ints, MPI_INT, send_data [j], send_sizes [j][2], &send_position [j], comm);
    }
  }

  /* pack doubles */
  for (i = 0, cd = cp->send; i < nsend; i ++, cd ++)
  {
    if (cd->doubles)
    {
      j = rankmap [cd->rank];
      MPI_Pack (cd->d, cd->doubles, MPI_DOUBLE, send_data [j], send_sizes [j][2], &send_position [j], comm);
    }
  }

  /* send data */
  for (i = 0; i < send_count; i ++)
  {
    MPI_Isend (send_data [i], send_sizes [i][2], MPI_PACKED, send_rank [i], tag, comm, &send_req [i]);
  }

  /* receive data */
  for (i = 0; i < recv_count; i ++)
  {
    MPI_Irecv (recv_data [i], recv_sizes [i][2], MPI_PACKED, recv_rank [i], tag, comm, &recv_req [i]);
  }

  return cp->send_size;
}

/* blocking receive */
void COM_Recv (void *pattern)
{
  COMPATTERN *cp = pattern;
  int (*recv_sizes) [3] = cp->recv_sizes,
        recv_count = cp->recv_count,
        send_count = cp->send_count,
        i, j;
  char **recv_data = cp->recv_data;
  MPI_Request *send_req = cp->send_req,
	      *recv_req = cp->recv_req;
  MPI_Status *send_sta = cp->send_sta,
             *recv_sta = cp->recv_sta;
  MPI_Comm comm = cp->comm;
  COMDATA *cd;

  /* wait until until send is done */
  MPI_Waitall (send_count, send_req, send_sta);

  /* wait until until receive is done */
  MPI_Waitall (recv_count, recv_req, recv_sta);

  /* unpack data */
  for (i = j = 0, cd = cp->recv; i < recv_count; i ++, cd ++, j = 0)
  {
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, cd->i, cd->ints, MPI_INT, comm);
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, cd->d, cd->doubles, MPI_DOUBLE, comm);
  }
}

/* free point to point communication pattern */
void COM_Free (void *pattern)
{
  COMPATTERN *cp = pattern;
  int i;

  free (cp->rankmap);
  free (cp->send_sizes);
  free (cp->send_position);
  free (cp->send_rank);
  for (i = 0; i < cp->send_count; i ++)
    free (cp->send_data [i]);
  free (cp->send_data);
  free (cp->recv_rank);
  free (cp->recv_sizes);
  for (i = 0; i < cp->recv_count; i ++)
    free (cp->recv_data [i]);
  free (cp->recv_data);
  free (cp->recv_req);
  free (cp->recv_sta);
  free (cp->send_req);
  free (cp->send_sta);
  free (pattern);
}

/* create a repetitive all to all communication pattern;
 * ranks and sizes must not change during the communication;
 * pointers to send and receive buffers data must not change */
void* COMALL_Pattern (MPI_Comm comm,
                      COMDATA *send, int nsend,
	              COMDATA **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMALLPATTERN *pattern;
  COMDATA *cd;
  int rank,
      ncpu,
      i, j, k;
  void *p;

  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &ncpu);

  ERRMEM (pattern = MEM_CALLOC (sizeof (COMALLPATTERN)));
  ERRMEM (pattern->send_sizes = MEM_CALLOC (ncpu * sizeof (int [3])));
  ERRMEM (pattern->send_counts = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (pattern->send_disps = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (pattern->send_position = MEM_CALLOC (ncpu * sizeof (int)));

  pattern->ncpu = ncpu;
  pattern->comm = comm;
  pattern->send = send;
  pattern->nsend = nsend;

  /* compute send sizes */
  for (i = 0, cd = send; i < nsend; i ++, cd ++)
  {
    pattern->send_sizes [cd->rank][0] += cd->ints;
    pattern->send_sizes [cd->rank][1] += cd->doubles;
    MPI_Pack_size (cd->ints, MPI_INT, comm, &j);
    MPI_Pack_size (cd->doubles, MPI_DOUBLE, comm, &k);
    pattern->send_sizes [cd->rank][2] += (j + k);
  }

  /* compute send displacements */
  for (pattern->send_size = i = 0; i < ncpu; i ++)
  {
    pattern->send_counts [i] = pattern->send_sizes [i][2];
    pattern->send_size += pattern->send_counts [i];
    if (i < (ncpu - 1)) pattern->send_disps [i+1] = pattern->send_size;
  }
  pattern->send_disps [0] = 0;

  /* allocate send buffer */
  ERRMEM (pattern->send_data = malloc (pattern->send_size));

  ERRMEM (pattern->recv_sizes = MEM_CALLOC (ncpu * sizeof (int [3])));
  ERRMEM (pattern->recv_counts = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (pattern->recv_disps = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (pattern->recv_position = MEM_CALLOC (ncpu * sizeof (int)));

  /* distribute send sizes into receive sizes */
  MPI_Alltoall (pattern->send_sizes, 3, MPI_INT, pattern->recv_sizes, 3, MPI_INT, comm);

  /* compute receive displacements */
  for (pattern->recv_size = i = 0; i < ncpu; i ++)
  {
    pattern->recv_counts [i] = pattern->recv_sizes [i][2];
    pattern->recv_size += pattern->recv_counts [i];
    if (i < (ncpu - 1)) pattern->recv_disps [i+1] = pattern->recv_size;
  }
  pattern->recv_disps [0] = 0;

  /* allocate receive buffer */
  ERRMEM (pattern->recv_data = malloc (pattern->recv_size));

  if (pattern->recv_size)
  {
    /* contiguous receive size */
    uint64_t size = ncpu * sizeof (COMDATA);
    for (i = 0; i < ncpu; i ++)
    {
      size += pattern->recv_sizes [i][0] * sizeof (int) + 
	      pattern->recv_sizes [i][1] * sizeof (double);
    }

    /* prepare output receive data */
    ERRMEM ((*recv) = malloc (size));
    p = (*recv) + ncpu;
    for (i = 0, cd = *recv; i < ncpu; i ++, cd ++)
    {
      cd->rank = i;
      cd->ints = pattern->recv_sizes [i][0];
      cd->doubles = pattern->recv_sizes [i][1];
      cd->i = p; p = (cd->i + cd->ints);
      cd->d = p; p = (cd->d + cd->doubles);
    }

    /* compress receive storage */
    for (*nrecv = i = 0; i < ncpu; i ++)
    {
      if (pattern->recv_counts [i])
      {
	(*recv) [*nrecv] = (*recv) [i];
	(*nrecv) ++;
      }
    }
  }
  else
  {
    *recv = NULL;
    *nrecv = 0;
  }

  pattern->recv = *recv;
  pattern->nrecv = *nrecv;

  return pattern;
}

/* communicate integers and doubles accodring
 * to the pattern computed by COMALL_Pattern */
uint64_t COMALL_Repeat (void *pattern)
{
  COMALLPATTERN *pp = pattern;
  COMDATA *cd;
  int i;

  for (i = 0; i < pp->ncpu; i ++) pp->send_position [i] = pp->recv_position [i] = 0;

  /* pack ints */
  for (i = 0, cd = pp->send; i < pp->nsend; i ++, cd ++)
  {
    if (cd->ints)
    {
      MPI_Pack (cd->i, cd->ints, MPI_INT, &pp->send_data [pp->send_disps [cd->rank]], pp->send_counts [cd->rank], &pp->send_position [cd->rank], pp->comm);
    }
  }

  /* pack doubles */
  for (i = 0, cd = pp->send; i < pp->nsend; i ++, cd ++)
  {
    if (cd->doubles)
    {
      MPI_Pack (cd->d, cd->doubles, MPI_DOUBLE, &pp->send_data [pp->send_disps [cd->rank]], pp->send_counts [cd->rank], &pp->send_position [cd->rank], pp->comm); 
    }
  }

#if DEBUG
  for (i = 0; i < pp->ncpu; i ++)
  {
    ASSERT_DEBUG (pp->send_position [i] <= pp->send_counts [i], "Incorrect packing");
  }
#endif

  /* all to all send and receive */
  MPI_Alltoallv (pp->send_data, pp->send_counts, pp->send_disps, MPI_PACKED, pp->recv_data, pp->recv_counts, pp->recv_disps, MPI_PACKED, pp->comm);

  if (pp->recv_size)
  {
    /* unpack data */
    for (i = 0; i < pp->ncpu; i ++)
    {
      MPI_Unpack (&pp->recv_data [pp->recv_disps [i]], pp->recv_counts [i], &pp->recv_position [i], pp->recv [i].i, pp->recv [i].ints, MPI_INT, pp->comm);
      MPI_Unpack (&pp->recv_data [pp->recv_disps [i]], pp->recv_counts [i], &pp->recv_position [i], pp->recv [i].d, pp->recv [i].doubles, MPI_DOUBLE, pp->comm);
    }
  }

  return pp->send_size;
}

/* free all to all communication pattern */
void COMALL_Free (void *pattern)
{
  COMALLPATTERN *pp = pattern;

  free (pp->send_sizes);
  free (pp->send_counts);
  free (pp->send_disps);
  free (pp->send_position);
  free (pp->recv_sizes);
  free (pp->recv_counts);
  free (pp->recv_disps);
  free (pp->recv_position);
  free (pp->send_data);
  free (pp->recv_data);
  free (pp);
}
