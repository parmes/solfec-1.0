/*
 * mem.c
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

#include <string.h>
#include <stdlib.h>
#include "mem.h"
#include "err.h"

#if MEMDEBUG
#include "set.h"
#endif

typedef struct { void *p; size_t margin; } PTR; /* pointer with margin */

void* MEM_CALLOC (size_t size)
{
  void *chunk;

  if (!(chunk = malloc (size))) return NULL;

  memset (chunk, 0, size);

  return chunk;
}

void MEM_Init (MEM *pool, size_t chunksize, size_t chunksinblock)
{
  ASSERT_DEBUG (pool && chunksize > 0 && chunksinblock > 0, "A zero argument passed to MEM_Init");

  /* set chunksize not less than 'unsigned long' size
   * as wee plan to use chunks as items of 'deadchunks' list */
  pool->chunksize = (chunksize > sizeof(PTR) ? chunksize : sizeof(PTR));
  pool->chunksinblock = chunksinblock;
  pool->blocks = NULL;
  pool->freechunk = NULL;
  pool->lastchunk = NULL;
  pool->deadchunks = NULL;
}

void* MEM_Alloc (MEM *pool)
{
#if MEMDEBUG
  void *chunk;

  if (!(chunk = malloc (pool->chunksize))) return NULL;
  memset (chunk, 0, pool->chunksize);

  SET_Insert (NULL, (SET**) &pool->blocks, chunk, NULL);

  return chunk;
#else
  void *chunk, *block;

  if (pool->deadchunks)
  { /* if there are deallocated chunks, get one */
	  
    chunk = pool->deadchunks;
    pool->deadchunks = ((PTR*)pool->deadchunks)->p;
    memset (chunk, 0, pool->chunksize);
    return chunk;
  }
  else if (pool->freechunk == pool->lastchunk)
  { /* else if we need to allocate a new block ... */
   
    /* allocate a block of memory */
    block = malloc (pool->chunksize * pool->chunksinblock + sizeof(PTR));
    if (!block) return NULL; /* do not exit() here */
    memset (block, 0, pool->chunksize * pool->chunksinblock + sizeof(PTR));
   
    /* insert allocated block into the list */
    ((PTR*)block)->p = pool->blocks;
    pool->blocks = block;
    /* set free and last chunk */
    pool->freechunk = (char*)block + sizeof(PTR);
    pool->lastchunk = pool->freechunk + pool->chunksize * pool->chunksinblock;
  }

  chunk = pool->freechunk;
  pool->freechunk += pool->chunksize;
  return chunk;
#endif
}

void MEM_Free (MEM *pool, void *chunk)
{
#if MEMDEBUG
  if (!SET_Contains ((SET*)pool->blocks, chunk, NULL))
  {
    ASSERT_DEBUG (0, "Deletion from invalid memory pool");
  }

  SET_Delete (NULL, (SET**) &pool->blocks, chunk, NULL);

  free (chunk);
#else
  /* insert chunk into dead chunks list */
  ((PTR*)chunk)->p = pool->deadchunks;
  pool->deadchunks = chunk;
#endif
}

size_t MEM_Size (MEM *pool)
{
#ifdef MEMDEBUG
  return SET_Size (pool->blocks) * pool->chunksize;
#else
  size_t size = 0, chunk = pool->chunksize *
    pool->chunksinblock + sizeof(PTR);

  void *block = pool->blocks;

  /* loop over all blocks */  
  while (block)
  {
    block = ((PTR*)block)->p;
    size += chunk;
  }

  return size;
#endif
}

void MEM_Release (MEM *pool)
{
#if MEMDEBUG
  for (SET *item = SET_First (pool->blocks); item; item = SET_Next (item)) free (item->data);

  SET_Free (NULL, (SET**) &pool->blocks);
#else
  void *block = pool->blocks;
  size_t next;
  
  /* traverse and free all blocks of memory */
  while (block)
  {
    next = *((size_t*)block);
    free (block);
    block = (void*)next;
  }
#endif

  pool->blocks = NULL;
  pool->freechunk = NULL;
  pool->lastchunk = NULL;
  pool->deadchunks = NULL;
}
