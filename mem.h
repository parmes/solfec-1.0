/*
 * mem.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * memory pool
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

#ifndef __mp__
#define __mp__

typedef struct memory_pool MEM;

struct memory_pool 
{
  void *blocks; /* list of allocated memory blocks */
  char *freechunk; /* next free chunk of memory */
  char *lastchunk; /* last chunk in current block */
  void *deadchunks; /* list of dealocated chunks of memory */
  size_t chunksize; /* size of a chunk */
  size_t chunksinblock; /* number of memory chunks in a block */
};

/* initialize memory pool */
void MEM_Init (MEM *pool, size_t chunksize, size_t chunksinblock);

/* allocate a chunk of memory from the pool */
void* MEM_Alloc (MEM *pool);

/* free a chunk of memory to the pool */
void MEM_Free (MEM *pool, void *chunk);

/* return amount of memory in the pool */
size_t MEM_Size (MEM *pool);

/* release memory pool memory back to system */
void MEM_Release (MEM *pool);

#endif
