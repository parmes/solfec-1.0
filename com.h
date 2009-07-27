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
typedef void  (*OBJ_Sizes)  (void *obj, int *ints, int *doubles); /* number of ints and doubles in an object */
typedef void  (*OBJ_Pack)   (void *obj, int *i, double *d); /* pack an object */
typedef void* (*OBJ_Unpack) (int *i, double *d); /* unpack an object */

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

/* communicate integers and doubles */
void COM (MPI_Comm comm, int tag,
          COMDATA *send, int nsend,
	  COMDATA **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* communicate objects */
void COMOBJS (MPI_Comm comm, int tag,
              OBJ_Sizes sizes,
	      OBJ_Pack pack,
	      OBJ_Unpack unpack,
              COMOBJ *send, int nsend,
	      COMOBJ **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* create a repetitive communication pattern;
 * ranks and sizes must not change during the
 * repetitive communication; pointers to send
 * and receive buffers data must not change */
void* COM_Pattern (MPI_Comm comm, int tag,
                   COMDATA *send, int nsend,
	           COMDATA **recv, int *nrecv); /* recv is contiguous => free (*recv) releases all memory */

/* communicate integers and doubles accodring
 * to the pattern computed by COM_Pattern */
void COM_Repeat (void *pattern);

/* free communication pattern */
void COM_Free (void *pattern);

#endif
