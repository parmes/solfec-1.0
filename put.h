/*
 * put.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * parallel utilities
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

#include "tmr.h"

#ifndef __put__
#define __put__

/* used to extend body geometric extents for domain balancing */
#define PUT_GEOMEPS 0.001 /* XXX => quite arbitrary */

/* get statistics on a vector of integer variables */
void PUT_int_stats (int n, int *val, int *sum, int *min, int *avg, int *max);

/* get statistics on a vector of integer variables; return 1 for rank 0 and 0 for others */
int PUT_root_int_stats (int n, int *val, int *sum, int *min, int *avg, int *max);

/* get statistics on a vector of double variables; return 1 for rank 0 and 0 for others */
int PUT_root_double_stats (int n, double *val, double *sum, double *min, double *avg, double *max);

/* parallel timer end: get maximum of all calls; return 1 for rank 0 and 0 for others */
int PUT_root_timerend (TIMING *t, double *time);

/* parallel timer end: return maximum of all calls */
double PUT_timerend (TIMING *t);

/* return sum of all calls */
int PUT_int_sum (int val);

/* return maximum of all calls */
int PUT_int_max (int val);

/* return minimum of all calls */
int PUT_int_min (int val);

/* return minimum of all calls and its rank */
int PUT_int_min_rank (int val, int *rank);

/* return minimum of all calls */
double PUT_double_min (double val);

/* return maximum of all calls */
double PUT_double_max (double val);

#endif
