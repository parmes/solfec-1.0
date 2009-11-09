/*
 * hsh.h
 * Copyright (C) 2005, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * ------------------------------------------------------------------------------
 * overlap detection  based on hashing box volumes and checking resulting list
 * of colliding hash table entries for overlaps.
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

#include "box.h"

#ifndef __hsh__
#define __hsh__

void* HASH_Create (int boxnum);

void HASH_Do (void *context, int boxnum, BOX **boxes, void *data, BOX_Overlap_Create report);

void HASH_Destroy (void *context);

#endif
