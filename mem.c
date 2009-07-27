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

void MEM_Init (MEM *pool,
  size_t chunksize, size_t chunksinblock)
{
  /* set chunksize not less than 'unsigned long' size
   * as wee plan to use chunks as items of 'deadchunks' list */
  pool->chunksize = (chunksize > sizeof(size_t) ? chunksize : sizeof(size_t));
  pool->chunksinblock = chunksinblock;
  pool->blocks = NULL;
  pool->freechunk = NULL;
  pool->lastchunk = NULL;
  pool->deadchunks = NULL;
}

void* MEM_Alloc (MEM *pool)
{
#if MEMDEBUG
  return calloc (pool->chunksize, 1);
#else
  void *chunk,*block;

  if (pool->deadchunks)
  { /* if there are deallocated chunks, get one */
	  
    chunk = pool->deadchunks;
    pool->deadchunks = (void *)*((size_t*)pool->deadchunks);
    memset (chunk, 0, pool->chunksize);
    return chunk;
  }
  else if (pool->freechunk == pool->lastchunk)
  { /* else if we need to allocate a new block ... */
   
    /* allocate a block of memory */
    block = calloc (pool->chunksize * pool->chunksinblock + sizeof(size_t), 1);
    if (! block)
    {
      return NULL; /* do not exit() here */
    }
   
    /* insert allocated block into the list */
    *((size_t*)block) = (size_t)pool->blocks;
    pool->blocks = block;
    /* set free and last chunk */
    pool->freechunk = (char *)block + sizeof(size_t);
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
  free (chunk);
#else
  /* insert chunk into dead chunks list */
  *((size_t*)chunk) = (size_t)pool->deadchunks;
  pool->deadchunks = chunk;
#endif
}

size_t MEM_Size (MEM *pool)
{
#ifdef MEMDEBUG
  return 0;
#else
  size_t size = 0, chunk = pool->chunksize *
    pool->chunksinblock + sizeof(size_t);

  void *block = pool->blocks;
  size_t next;

  /* loop over all blocks */  
  while (block)
  {
    next = *((size_t*)block);
    block = (void*)next;
    size += chunk;
  }

  return size;
#endif
}

void MEM_Release (MEM *pool)
{
#ifndef MEMDEBUG
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
