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

#include <mpi.h>

#ifndef __com__
#define __com__

typedef struct comdata COMDATA; /* ints and doubles */
typedef struct comobj COMOBJ; /* objects */
typedef void  (*OBJ_Pack)   (void *obj, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints); /* pack an object */
typedef void* (*OBJ_Unpack) (void *data, int *dpos, double *d, int doubles, int *ipos, int *i, int ints); /* unpack an object */

struct comdata
{
  int rank, /* send or receive rank */
      ints, /* integers count */
      doubles; /* doubles count */

  int *i; /* integers */

  double *d; /* doubles */
};

struct comobj
{
  int rank; /* send or receive rank */

  void *o; /* object */
};

/* communicate integers and doubles using point to point communication */
void COM (MPI_Comm comm, int tag,
          COMDATA *send, int nsend,
	  COMDATA **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* communicate integers and doubles using all to all communication */
void COMALL (MPI_Comm comm,
             COMDATA *send, int nsend,
	     COMDATA **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* communicate one set of integers and doubles to all other processors */
void COMONEALL (MPI_Comm comm, COMDATA send,
	        COMDATA **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* communicate objects using point to point communication */
void COMOBJS (MPI_Comm comm, int tag,
	      OBJ_Pack pack,
	      void *data,
	      OBJ_Unpack unpack,
              COMOBJ *send, int nsend,
	      COMOBJ **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* communicate objects using all to all communication */
void COMOBJSALL (MPI_Comm comm,
	         OBJ_Pack pack,
	         void *data,
	         OBJ_Unpack unpack,
                 COMOBJ *send, int nsend,
	         COMOBJ **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* communicate an object to all other processors */
void COMOBJALL (MPI_Comm comm,
	        OBJ_Pack pack,
	        void *data,
	        OBJ_Unpack unpack,
		void *object,
	        COMOBJ **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* create a repetitive point to point communication pattern;
 * ranks and sizes must not change during the communication;
 * pointers to send and receive buffers data must not change */
void* COM_Pattern (MPI_Comm comm, int tag,
                   COMDATA *send, int nsend,
	           COMDATA **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* communicate integers and doubles accodring
 * to the pattern computed by COM_Pattern */
void COM_Repeat (void *pattern);

/* non-blocking send */
void COM_Send (void *pattern);

/* blocking receive */
void COM_Recv (void *pattern);

/* free communication pattern */
void COM_Free (void *pattern);

#endif
