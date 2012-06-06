/*
 * hul.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * convex hull in three dimensions according to the algorithm by
 * Barber et al. "The Quickhull Algorithm for Convex Hulls"
 * ACM Transactions on Mathematical Software, 1996
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

#include "tri.h"

#ifndef __hul__
#define __hul__

/* take n 'v'ertices and output m elements of the doubly connected edge list;
 * note that in the returend table of triangles of size 'm' all vertices point
 * to the memory in 'v'; return NULL if hull creation failed from geometrical
 * reasons; throw memory exception when out of memory */
TRI* hull (double *v, int n, int *m);

#endif
