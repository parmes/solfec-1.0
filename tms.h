/*
 * tms.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * time series
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

#include <stdio.h>

#ifndef __tms__
#define __tms__

typedef struct time_series TMS;

struct time_series
{
  double value; /* constant value => used if size == 0 */
  double (*points) [2]; /* vector of (time, value) pairs */
  char *label; /* label for globally mapped time series */
  int marker; /* index of the last read interval */
  int size; /* number of pairs */
  char *path; /* file path for partially cached time series */
  long *offset; /* file offsets */
  double *time; /* time instants maching file offsets */
  int noffsets; /* number of offsets */
  int cache; /* partial cache size */
};

/* create a copy */
TMS* TMS_Copy (TMS *ts);

/* create time series */
TMS* TMS_Create (int size, double *times, double *values, char *label);

/* create time series from a text file */
TMS* TMS_File (char *path, char *label, int cache);

/* wrapper for a constant value */
TMS* TMS_Constant (double value, char *label);

/* create from another series through integration */
TMS* TMS_Integral (TMS *ts);

/* as above, but calculate derivative */
TMS* TMS_Derivative (TMS *ts);

/* obtain value at a specific time */
double TMS_Value (TMS *ts, double time);

/* output the series with a prescribed step => or
 * with the nominal acuracy if 'step' <= 0.0 */
void TMS_Output (TMS *ts, char *path, double step);

/* free time series memory */
void TMS_Destroy (TMS *ts);

/* pack time series into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void TMS_Pack (TMS *ts, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack time series from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
TMS* TMS_Unpack (int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

/* export MBFCP definition */
void TMS_2_MBFCP (TMS *tms, FILE *out);

/* released globally mapped time series */
void TMS_RELEASE_GLOBAL_MAP();

#endif
