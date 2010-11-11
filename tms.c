/*
 * tms.c
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

#include <stdlib.h>
#include <string.h>
#include "mem.h"
#include "tms.h"
#include "pck.h"
#include "err.h"

#define CHUNK 512

static int findmarker (double (*begin)[2], double (*end)[2], double time)
{
  double (*low)[2] = begin,
	 (*high)[2] = end,
	 (*mid)[2];

  while (low <= high)
  {
    mid = low + (high-low) / 2;

    if (time >= (*mid)[0] &&
        time < (*(mid+1))[0])
    {
      return (mid - begin);
    }
    else if (time < (*mid)[0])
    {
      high = mid - 1;
    }
    else
    {
      low = mid + 1;
    }
  }

  return 0;
}

static double linterp (double (*point) [2], double time)
{
  double dt = time - point[0][0];
  return point[0][1] + (point[1][1]-point[0][1]) * (dt / (point[1][0] - point[0][0]));
}

TMS* TMS_Copy (TMS *ts)
{
  TMS *out;

  ERRMEM (out = malloc (sizeof (TMS)));
  out->value = ts->value;
  ERRMEM (out->points = malloc (sizeof (double [2]) * ts->size));
  out->marker = ts->marker;
  out->size = ts->size;
  memcpy (out->points, ts->points, sizeof (double [2]) * ts->size);
  return out;
}

TMS* TMS_Create (int size, double *times, double *values)
{
  TMS *ts;
  int i;

  ASSERT_DEBUG (size > 1, "Time series size must be greater than one");
  ERRMEM (ts = malloc (sizeof (TMS)));
  ERRMEM (ts->points = malloc (sizeof (double [2]) * size));
  ts->marker = 0;
  ts->size = size;

  for (i = 0; i < size; i ++)
  {
    ts->points [i][0] = times [i];
    ts->points [i][1] = values [i];

    if (i)
    {
      ASSERT (times [i] > times [i-1], ERR_TMS_TIME_NOT_INCREASED);
    }
  }

  return ts;
}

TMS* TMS_File (char *path)
{
  int np;
  FILE *fp;
  TMS *ts;

  if (!(fp = fopen (path, "r"))) return NULL;
  ERRMEM (ts = malloc (sizeof (TMS)));
  ERRMEM (ts->points = MEM_CALLOC (sizeof (double [2]) * CHUNK));

  np = CHUNK;
  ts->size = 0;
  while (feof (fp) == 0)
  {
    if (fscanf (fp, "%lf", &ts->points [ts->size][0]) == EOF) break;
    if (fscanf (fp, "%lf", &ts->points [ts->size][1]) == EOF) break;

    if (ts->size)
    {
      ASSERT (ts->points [ts->size][0] > ts->points [ts->size-1][0], ERR_TMS_TIME_NOT_INCREASED);
    }

    if (++ ts->size >= np)
    {
      np += CHUNK;
      ERRMEM (ts->points = realloc (ts->points, sizeof (double [2]) * np));
    }
  }

  ts->size --;

  if (ts->size == 0)
  {
    free (ts->points);
    free (ts);
    ASSERT (0, ERR_FILE_EMPTY);
  }

  ts->marker = 0;

  return ts;
}

TMS* TMS_Constant (double value)
{
  TMS *ts;

  ERRMEM (ts = malloc (sizeof (TMS)));
  ts->value = value;
  ts->size = 0; /* indicates constant value */
  return ts;
}

TMS* TMS_Integral (TMS *ts)
{
  double (*pin) [2],
	 (*pout) [2];
  TMS *out;
  int n;

  ASSERT (ts->size > 0, ERR_TMS_INTEGRATE_CONSTANT);
  ERRMEM (out = malloc (sizeof (TMS)));
  out->size = ts->size;
  ERRMEM (out->points = MEM_CALLOC (sizeof (double [2]) * out->size));

  if (out->size == 0) return TMS_Constant ((ts->points[1][0] -
    ts->points[0][0]) * 0.5 *  (ts->points[0][1] + ts->points[1][1]));

  out->marker = 0;
  pin = ts->points;
  pout = out->points;

  pout [0][0] = pin[0][0];
  pout [0][1] = pin[0][1];

  for (n = 1; n < out->size; n ++)
  {
    pout [n][0] = pin [n][0];
    pout [n][1] = pout [n-1][1] + (pin[n][0] - pin[n-1][0]) * 0.5 * (pin[n-1][1] + pin[n][1]);
  }

  return out;
}

TMS* TMS_Derivative (TMS *ts)
{
  double (*pin) [2],
	 (*pout) [2];
  TMS *out;
  int n;

  if (ts->size == 0) return TMS_Constant (0.0);

  ERRMEM (out = malloc (sizeof (TMS)));
  out->size = ts->size - 1; /* one value per interval */
  ERRMEM (out->points = MEM_CALLOC (sizeof (double [2]) * out->size));

  if (out->size == 0) return TMS_Constant ((ts->points[1][0] -
    ts->points[0][0]) * 0.5 *  (ts->points[0][1] + ts->points[1][1]));

  out->marker = 0;
  pin = ts->points;
  pout = out->points;

  for (n = 0; n < out->size; n ++)
  {
    pout [n][0] = 0.5 * (pin[n][0] + pin[n+1][0]);
    pout [n][1] = (pin[n+1][1] - pin[n][1]) / (pin[n+1][0] - pin[n][0]);
  }

  return out;
}

double TMS_Value (TMS *ts, double time)
{
  double lo, hi;

  if (ts->size == 0) return ts->value;

  if (time < ts->points[0][0]) return ts->points[0][1];
  else if (time > ts->points[ts->size-1][0]) return ts->points[ts->size-1][1];

  lo = ts->points[ts->marker > 0 ? ts->marker - 1 : ts->marker][0];
  hi = ts->points[ts->marker < ts->size - 1 ? ts->marker + 1 : ts->marker][0];

  if (time < lo || time > hi)
  {
    ts->marker = findmarker (ts->points, ts->points + ts->size - 1, time);
  }
  else if (time >= lo && ts->marker &&
	   time < ts->points[ts->marker][0]) ts->marker --;
  else if (time >= ts->points[ts->marker+1][0] &&
	   time < hi && ts->marker < ts->size - 1) ts->marker ++;

  return linterp (&ts->points[ts->marker], time);
}

void TMS_Output (TMS *ts, char *path, double step)
{
  double (*p) [2], t, e;
  FILE *fp;
  int n;

  ASSERT (fp = fopen (path, "w"), ERR_FILE_OPEN);

  if (ts->size == 0)
  {
    fprintf (fp, "%e\n", TMS_Value (ts, 0.0));
  }
  else
  {
    if (step > 0.0)
    {
      for (t = ts->points[0][0],
	   e = ts->points[ts->size-1][0]; t < e; t += step)
	fprintf (fp, "%e\t%e\n", t, TMS_Value (ts, t));
    }
    else
    {
      p = ts->points;
      for (n = 0; n < ts->size; n ++)
	fprintf (fp, "%e\t%e\n", p[n][0], p[n][1]);
    }
  }

  fclose (fp);
}

void TMS_Destroy (TMS *ts)
{
  if (ts->size == 0) free (ts);
  else
  {
    free (ts->points);
    free (ts);
  }
}

void TMS_Pack (TMS *ts, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_double (dsize, d, doubles, ts->value);
  pack_int (isize, i, ints, ts->marker);
  pack_int (isize, i, ints, ts->size);
  if (ts->size) pack_doubles (dsize, d, doubles, (double*)ts->points, ts->size * 2);
}

TMS* TMS_Unpack (int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  TMS *ts;

  ERRMEM (ts = malloc (sizeof (TMS)));
  ts->value = unpack_double (dpos, d, doubles);
  ts->marker = unpack_int (ipos, i, ints);
  ts->size = unpack_int (ipos, i, ints);
  if (ts->size)
  {
    ERRMEM (ts->points = malloc (sizeof (double [2]) * ts->size));
    unpack_doubles (dpos, d, doubles, (double*)ts->points, ts->size * 2);
  }

  return ts;
}

/* export MBFCP definition */
void TMS_2_MBFCP (TMS *tms, FILE *out)
{
  int n;

  if (tms->size == 0)
  {
    fprintf (out, "CONSTANT:\t%g\n", tms->value);
  }
  else
  {
    fprintf (out, "TIME_SERIES:\t%d\n", tms->size);
    for (n = 0; n < tms->size; n ++) fprintf (out, "%g\t%g\n", tms->points [n][0], tms->points[n][1]);
  }
}
