/*
 * xyt.h
 * Copyright (C) 2005, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * Priority search tree
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

#include "mem.h"

#ifndef __xyt__
#define __xyt__

typedef struct xytree XYTREE;
typedef void (*XY_Callback) (void *pointer, void *data);

struct xytree
{
  void *key, *data;
  double x, y;
  XYTREE *link; /* heap link */
  char type; /* node or leaf */
  char colour; /* red or black */
  XYTREE *parent, *left, *right; /* binary tree links */
};

/*
 * Insert an item into PST
 * Pass valid mempool pointer (adjusted to 'xytree' chunks),
 * tree root **pointer and data point. All the 'data' pointers
 * should be different if one wants to store points with repeated
 * x coordinates.
 * O (log N)
 */
void XY_Insert (MEM *pool, XYTREE **tree, double x, double y, void *key, void *data);

/*
 * Delete an item from PST
 * O (log N)
 */
void XY_Delete (MEM *pool, XYTREE **tree, double x, double y, void *key);

/*
 * Two sided query.
 * Report all the 'data' associated with points
 * above xmin and on right from ymin.
 * O (log N + K)
 */
void XY_Query (XYTREE *tree, double xmin, double ymin, void *pointer, XY_Callback callback);

/* This structure is also capable of enabling an efficient
 * three-sided query, which was is not needed here */

#endif
