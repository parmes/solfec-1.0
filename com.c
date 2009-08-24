/*
 * com.h
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

#if OWNASYNC
#include <pthread.h>
#endif

#include <stdlib.h>
#include <mpi.h>
#include "com.h"
#include "err.h"

typedef struct compattern COMPATTERN;

/* integer and doubles
 * communication pattern */
struct compattern
{
  MPI_Comm comm; 

  int tag,
     *rankmap, /* maps ranks to send buffers indices */
    (*send_sizes) [3],
     *send_position,
     *send_rank,
      send_count,
     *recv_rank,
    (*recv_sizes) [3],
      recv_count;

  char **send_data,
       **recv_data;

  MPI_Request *req; /* for receiving */
  MPI_Status *sta; /* for receiving */

  MPI_Request *reqs; /* for sending */
  MPI_Status *stas; /* for sending */

  COMDATA *send,
	  *recv;

  int nsend,
      nrecv;

#if OWNASYNC
  int active,
      wait,
      done;

  pthread_t thread;
#endif
};

/* communicate integers and doubles */
void COM (MPI_Comm comm, int tag,
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
  char **send_data,
       **recv_data;
  MPI_Request *req;
  MPI_Status *sta;
  void *p;

  MPI_Comm_rank (comm, &rank);
  MPI_Comm_size (comm, &ncpu);

  ERRMEM (send_sizes = calloc (ncpu, sizeof (int [3])));
  ERRMEM (send_position = calloc (ncpu, sizeof (int)));
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
  for (i = 0; i < ncpu; i ++)
  {
    if (send_sizes [i][2])
    {
      ERRMEM (send_data [i] = malloc (send_sizes [i][2]));
      send_position [i] = 0;
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
  j = recv_count * sizeof (COMDATA);
  for (i = 0; i < recv_count; i ++)
  {
    j += recv_sizes [i][0] * sizeof (int) + 
         recv_sizes [i][1] * sizeof (double);
  }

  /* prepare receive buffers */
  ERRMEM (recv_data = malloc (recv_count * sizeof (char*)));
  ERRMEM ((*recv) = malloc (j));
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
}

/* communicate objects */
void COMOBJS (MPI_Comm comm, int tag,
	      OBJ_Pack pack,
	      void *data,
	      OBJ_Unpack unpack,
              COMOBJ *send, int nsend,
	      COMOBJ **recv, int *nrecv) /* recv is contiguous => free (*recv) releases all memory */
{
  COMDATA *send_data,
	  *recv_data,
	  *cd;
  int recv_count,
      i, n;
  COMOBJ *co;

  ERRMEM (send_data = malloc (nsend * sizeof (COMDATA)));

  /* pack objects */
  for (i = 0, cd = send_data, co = send; i < nsend; i ++, cd ++, co ++)
  {
    int isize = 0,
	dsize = 0;

    cd->rank = co->rank;
    cd->ints = 0;
    cd->doubles = 0;
    cd->i = NULL;
    cd->d = NULL;

    pack (co->o, &dsize, &cd->d, &cd->doubles, &isize, &cd->i, &cd->ints);
  }

  /* send and receive packed data */
  COM (comm, tag, send_data, nsend, &recv_data, &recv_count);

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

  /* cleanup */
  for (i = 0, cd = send_data; i < nsend; i ++, cd ++)
  { free (cd->i);
    free (cd->d); }
  free (send_data);
  free (recv_data); /* contiguous */
}


#if OWNASYNC
static void* ownasync (void *pattern)
{
  COMPATTERN *cp = pattern;

  while (cp->active)
  {
    while (cp->wait);

    if (cp->active)
    {
      cp->wait = 1;

      COM_Repeat (cp);

      cp->done = 1;
    }
  }

  return NULL;
}
#endif

/* create a repetitive communication pattern;
 * ranks and sizes must not change during the
 * repetitive communication; pointers to send
 * and receive buffers data must not change */
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
  ERRMEM (pattern->rankmap = calloc (ncpu, sizeof (int)));
  ERRMEM (pattern->send_sizes = calloc (ncpu, sizeof (int [3])));
  ERRMEM (pattern->send_position = calloc (ncpu, sizeof (int)));
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
  for (i = j = 0; i < ncpu; i ++)
  {
    if (pattern->send_sizes [i][2])
    {
      ERRMEM (pattern->send_data [i] = malloc (pattern->send_sizes [i][2]));
      pattern->rankmap [i] = j;
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
  ERRMEM (pattern->req = malloc (pattern->recv_count * sizeof (MPI_Request)));
  ERRMEM (pattern->sta = malloc (pattern->recv_count * sizeof (MPI_Status)));
  ERRMEM (pattern->reqs = malloc (pattern->send_count * sizeof (MPI_Request)));
  ERRMEM (pattern->stas = malloc (pattern->send_count * sizeof (MPI_Status)));

  /* communicate receive sizes */
  for (i = 0; i < pattern->recv_count; i ++)
  {
    MPI_Irecv (pattern->recv_sizes [i], 3, MPI_INT, pattern->recv_rank [i], tag, comm, &pattern->req [i]);
  }
  MPI_Barrier (comm);
  for (i = 0; i < pattern->send_count; i ++)
  {
    MPI_Rsend (pattern->send_sizes [i], 3, MPI_INT, pattern->send_rank [i], tag, comm);
  }
  MPI_Waitall (pattern->recv_count, pattern->req, pattern->sta);

  /* contiguous receive size */
  j = pattern->recv_count * sizeof (COMDATA);
  for (i = 0; i < pattern->recv_count; i ++)
  {
    j += pattern->recv_sizes [i][0] * sizeof (int) + 
         pattern->recv_sizes [i][1] * sizeof (double);
  }

  /* prepare receive buffers */
  ERRMEM (pattern->recv_data = malloc (pattern->recv_count * sizeof (char*)));
  ERRMEM (pattern->recv = malloc (j));
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
  free (send_count_all);
  free (send_rank_all);

  /* output */
  *nrecv = pattern->nrecv;
  *recv = pattern->recv;

#if OWNASYNC
  pattern->active = 1;
  pattern->wait = 1;
  pattern->done = 0;

  pthread_create (&pattern->thread, NULL, ownasync, pattern);
#endif

  return pattern;
}

/* communicate integers and doubles accodring
 * to the pattern computed by COM_Pattern */
void COM_Repeat (void *pattern)
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
  MPI_Request *req = cp->req;
  MPI_Status *sta = cp->sta;
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
    MPI_Irecv (recv_data [i], recv_sizes [i][2], MPI_PACKED, recv_rank [i], tag, comm, &req [i]);
  }
  MPI_Barrier (comm);
  for (i = 0; i < send_count; i ++)
  {
    MPI_Rsend (send_data [i], send_sizes [i][2], MPI_PACKED, send_rank [i], tag, comm);
  }
  MPI_Waitall (recv_count, req, sta);

  /* unpack data */
  for (i = j = 0, cd = cp->recv; i < recv_count; i ++, cd ++, j = 0)
  {
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, cd->i, cd->ints, MPI_INT, comm);
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, cd->d, cd->doubles, MPI_DOUBLE, comm);
  }
}

/* non-blocking send */
void COM_Send (void *pattern)
{
  COMPATTERN *cp = pattern;
#if OWNASYNC
  cp->done = 0;
  cp->wait = 0;
#else
  int *rankmap = cp->rankmap,
      *send_position = cp->send_position,
     (*send_sizes) [3] = cp->send_sizes,
      *send_rank = cp->send_rank,
       send_count = cp->send_count,
       tag = cp->tag,
       nsend = cp->nsend,
       i, j;
  char **send_data = cp->send_data;
  MPI_Request *reqs = cp->reqs;
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
    MPI_Isend (send_data [i], send_sizes [i][2], MPI_PACKED, send_rank [i], tag, comm, &reqs [i]);
  }
#endif
}

/* blocking receive */
void COM_Recv (void *pattern)
{
  COMPATTERN *cp = pattern;
#if OWNASYNC
  while (!cp->done);
#else
  int *recv_rank = cp->recv_rank,
     (*recv_sizes) [3] = cp->recv_sizes,
       recv_count = cp->recv_count,
       send_count = cp->send_count,
       tag = cp->tag,
       i, j;
  char **recv_data = cp->recv_data;
  MPI_Request *reqs = cp->reqs;
  MPI_Status *stas = cp->stas;
  MPI_Status *sta = cp->sta;
  MPI_Comm comm = cp->comm;
  COMDATA *cd;

  /* wait until until send is done */
  MPI_Waitall (send_count, reqs, stas);

  /* receive data */
  for (i = 0; i < recv_count; i ++)
  {
    MPI_Recv (recv_data [i], recv_sizes [i][2], MPI_PACKED, recv_rank [i], tag, comm, &sta [i]);
  }

  /* unpack data */
  for (i = j = 0, cd = cp->recv; i < recv_count; i ++, cd ++, j = 0)
  {
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, cd->i, cd->ints, MPI_INT, comm);
    MPI_Unpack (recv_data [i], recv_sizes [i][2], &j, cd->d, cd->doubles, MPI_DOUBLE, comm);
  }
#endif
}

/* free communication pattern */
void COM_Free (void *pattern)
{
  COMPATTERN *cp = pattern;
  int i;

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
  free (cp->req);
  free (cp->sta);
  free (cp->reqs);
  free (cp->stas);

#if OWNASYNC
  void *r;

  cp->active = 0;
  cp->wait = 0;

  pthread_join (cp->thread, &r);
#endif

  free (pattern);
}
