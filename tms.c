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
#include <float.h>
#include "mem.h"
#include "tms.h"
#include "pck.h"
#include "err.h"
#include "map.h"

/* a map of globally labeled time series */
static MAP *labeled_time_series = NULL;

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

static int findoffset (double *begin, double *end, double time)
{
  double *low = begin,
	 *high = end,
	 *mid = begin + (end-begin) / 2;

  while (low < mid && mid < high)
  {
    if (time < (*mid))
    {
      high = mid;
    }
    else if (time >= (*mid))
    {
      low = mid;
    }

    mid = low + (high-low) / 2;
  }

  return mid - begin;
}

static void applyderivative (TMS *ts)
{
  double (*points) [2] = ts->points;
  double (*pin) [2], (*pout) [2];
  int size = ts->size, n;

#if 0
  ts->size = size - 1; /* central difference: one value per interval */

  ERRMEM (ts->points = MEM_CALLOC (sizeof (double [2]) * ts->size));

  ts->marker = 0;
  pin = points;
  pout = ts->points;

  for (n = 0; n < ts->size; n ++)
  {
    pout [n][0] = 0.5 * (pin[n][0] + pin[n+1][0]);
    pout [n][1] = (pin[n+1][1] - pin[n][1]) / (pin[n+1][0] - pin[n][0]);
  }
#else
  ts->size = 2*(size-1); /* constant value per interval with linear trend: two ends of the interval */

  ERRMEM (ts->points = MEM_CALLOC (sizeof (double [2]) * ts->size));

  ts->marker = 0;
  pin = points;
  pout = ts->points;

  for (n = 0; n < size-1; n ++)
  {
    double dt = pin[n+1][0] - pin[n][0], eps = 1E-10 * dt;

    pout [2*n][0] = pin[n][0];
    pout [2*n+1][0] = pin[n+1][0] - eps;
    pout [2*n][1] =
    pout [2*n+1][1] = (pin[n+1][1] - pin[n][1]) / dt;
  }
#endif

  free (points);
}

static double linterp (double (*point) [2], double time)
{
  double dt = time - point[0][0];
  return point[0][1] + (point[1][1]-point[0][1]) * (dt / (point[1][0] - point[0][0]));
}

#if 0
static TMS* bylabel (char *label)
{
  if (label)
  {
    return MAP_Find (labeled_time_series, label, (MAP_Compare)strcmp);
  }
  else return NULL;
}
#endif

static void addlabel (TMS *ts, char *label)
{
  if (label)
  {
    ERRMEM (ts->label = malloc(strlen(label)+1));
    strcpy (ts->label, label);
    ASSERT_TEXT (MAP_Find (labeled_time_series, label, (MAP_Compare)strcmp) == NULL,
                 "A time series object with label %s already exists", label);
    MAP_Insert (NULL, &labeled_time_series, ts->label, ts, (MAP_Compare)strcmp);
  }
  else ts->label = NULL;
}

static TMS* hardcopy (TMS *ts)
{
  TMS *out;

  ERRMEM (out = MEM_CALLOC (sizeof (TMS)));
  out->value = ts->value;
  ERRMEM (out->points = malloc (sizeof (double [2]) * ts->size));
  out->marker = ts->marker;
  out->size = ts->size;
  memcpy (out->points, ts->points, sizeof (double [2]) * ts->size);
  out->label = NULL;
  if (ts->cache)
  {
    out->cache = ts->cache;
    ERRMEM (out->path = malloc (strlen(ts->path)+1));
    strcpy (out->path, ts->path);
    ERRMEM (out->offset = malloc (ts->noffsets * sizeof(long)));
    memcpy (out->offset, ts->offset, ts->noffsets * sizeof(long));
    ERRMEM (out->time = malloc (ts->noffsets* sizeof(long)));
    memcpy (out->time, ts->time, ts->noffsets* sizeof(long));
    out->noffsets = ts->noffsets;
    out->op = ts->op;
  }

  return out;
}

TMS* TMS_Copy (TMS *ts)
{
  if (ts->label)
  {
    TMS *out = MAP_Find (labeled_time_series, ts->label, (MAP_Compare)strcmp);
    ASSERT_TEXT (out, "Labeled time series is not present in the global map");
    return NULL;
  }
  else return hardcopy(ts);
}

TMS* TMS_Create (int size, double *times, double *values, char *label)
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

  addlabel (ts, label);

  ts->cache = 0;
  ts->path = NULL; /* not partially cached */
  ts->offset = NULL;
  ts->time = NULL;
  ts->noffsets = 0;
  ts->op = 0;

  return ts;
}

TMS* TMS_File (char *path, char *label, int cache)
{
  int np, no, size;
  char buf [4096];
  double t0, t, v;
  long off0;
  FILE *fp;
  TMS *ts;

  if (!(fp = fopen (path, "r"))) return NULL;
  ERRMEM (ts = malloc (sizeof (TMS)));
  ERRMEM (ts->points = MEM_CALLOC (sizeof (double [2]) * CHUNK));

  if (cache)
  {
    ERRMEM (ts->path = malloc (strlen(path)+1));
    strcpy (ts->path, path);
    ERRMEM (ts->offset = MEM_CALLOC (sizeof (long) * CHUNK));
    ERRMEM (ts->time = MEM_CALLOC (sizeof (double) * CHUNK));
    ts->noffsets = 1;
  }

  ts->cache = cache;

  np = no = CHUNK;
  ts->size = 0;
  size = 0;
  t = -DBL_MAX;
  while ((off0 = ftell(fp)) >= 0 && fgets(buf, sizeof(buf), fp) != NULL) /* read line */
  {
    if (buf[0] == '#') continue; /* skip comments */

    if (cache == 0 || ts->size <= cache)
    {
      t0 = t;

      if (sscanf (buf, "%lf%lf", &t, &v) == EOF) break;

      ASSERT (t > t0, ERR_TMS_TIME_NOT_INCREASED);

      ts->points [ts->size][0] = t;
      ts->points [ts->size][1] = v;

      if (++ ts->size >= np)
      {
	np += CHUNK;
	ERRMEM (ts->points = realloc (ts->points, sizeof (double [2]) * np));
      }

      size ++;
    }
    else
    {
      t0 = t;

      if (sscanf (buf, "%lf%lf", &t, &v) == EOF) break;

      ASSERT (t > t0, ERR_TMS_TIME_NOT_INCREASED);

      size ++;
    }

    if (ts->size == 1)
    {
      ts->offset[0] = off0;
      ts->time[0] = t;
    }

    if (size == cache)
    {
      ts->offset[ts->noffsets] = off0;
      ts->time[ts->noffsets] = t;

      if (++ ts->noffsets >= no)
      {
	no += CHUNK;
	ERRMEM (ts->offset = realloc (ts->offset, sizeof (long) * no));
	ERRMEM (ts->time = realloc (ts->time, sizeof (double) * no));
      }

      size = 0;
    }
  }

  if (ts->size == 0)
  {
    free (ts->points);
    if (cache)
    {
      free (ts->path);
      free (ts->offset);
      free (ts->time);
    }
    free (ts);
    ASSERT (0, ERR_FILE_EMPTY);
  }
  else
  {
    ERRMEM (ts->points = realloc (ts->points, sizeof (double [2]) * ts->size)); /* trim */

    if (cache)
    {
      if (ts->time[ts->noffsets-1] < t)
      {
	ts->offset[ts->noffsets] = off0;
	ts->time[ts->noffsets] = t;
	ts->noffsets ++;
      }

      ERRMEM (ts->offset = realloc (ts->offset, sizeof (long) * ts->noffsets)); /* trim */
      ERRMEM (ts->time = realloc (ts->time, sizeof (double) * ts->noffsets));
    }
  }

  ts->marker = 0;

  addlabel (ts, label);

  if (cache == 0) /* not partially cached */
  {
    ts->cache = 0;
    ts->path = NULL;
    ts->offset = NULL;
    ts->time = NULL;
    ts->noffsets = 0;
  }

  ts->op = 0;

  fclose (fp);

  return ts;
}

TMS* TMS_Constant (double value, char *label)
{
  TMS *ts;

  ERRMEM (ts = malloc (sizeof (TMS)));
  ts->value = value;
  ts->size = 0; /* indicates constant value */

  addlabel (ts, label);

  ts->cache = 0;
  ts->path = NULL; /* not partially cached */
  ts->offset = NULL;
  ts->time = NULL;
  ts->noffsets = 0;
  ts->op = 0;

  return ts;
}

TMS* TMS_Integral (TMS *ts)
{
  double (*pin) [2],
	 (*pout) [2];
  TMS *out;
  int n;

  ASSERT_TEXT (ts->path == NULL, "Integral of partially cached time series objects is not supported");

  ASSERT (ts->size > 1, ERR_TMS_INTEGRATE_CONSTANT);
  ERRMEM (out = malloc (sizeof (TMS)));
  out->size = ts->size;

  if (out->size == 2)
  {
    free (out);
    return TMS_Constant ((ts->points[1][0] - ts->points[0][0]) * 0.5 *  (ts->points[0][1] + ts->points[1][1]), NULL);
  }

  ERRMEM (out->points = MEM_CALLOC (sizeof (double [2]) * out->size));

  out->marker = 0;
  pin = ts->points;
  pout = out->points;

  pout [0][0] = pin[0][0];
  pout [0][1] = 0.0;

  for (n = 1; n < out->size; n ++)
  {
    pout [n][0] = pin [n][0];
    pout [n][1] = pout [n-1][1] + (pin[n][0] - pin[n-1][0]) * 0.5 * (pin[n-1][1] + pin[n][1]);
  }

  out->label = NULL;
  out->cache = 0;
  out->path = NULL; /* not partially cached */
  out->offset = NULL;
  out->time = NULL;
  out->noffsets = 0;
  out->op = 0;

  return out;
}

TMS* TMS_Derivative (TMS *ts)
{
  TMS *out;

  if (ts->size == 0) return TMS_Constant (0.0, NULL);

  if (ts->size == 1)
  {
    return TMS_Constant ((ts->points[1][1] - ts->points[0][1])/(ts->points[1][0] - ts->points[0][0]), NULL);
  }

  out = hardcopy (ts);

  applyderivative (out);

  out->label = NULL;

  if (out->cache) out->op = 1; /* partially cached derivative */

  return out;
}

double TMS_Value (TMS *ts, double time)
{
  double lo, hi;

  if (ts->size == 0) return ts->value;

  /* if the time series is partially cached - when out of bounds - reload cache */
  if ((time < ts->points[0][0] || time > ts->points[ts->size-1][0]) && ts->path)
  {
    char buf [4096];
    int off = findoffset (ts->time, ts->time + ts->noffsets - 1, time);
    FILE *fp = fopen (ts->path, "r");
    ASSERT_TEXT (fp, "Opening data file [%s] for partially cached TIME_SERIES file has failed", ts->path);
    fseek (fp, ts->offset[off], SEEK_SET);
    ts->size = 0;
    while (ts->size <= ts->cache && fgets(buf, sizeof(buf), fp) != NULL) /* read line */
    {
      if (buf[0] == '#') continue; /* skip comments */

      if (sscanf (buf, "%lf%lf", &ts->points [ts->size][0], &ts->points [ts->size][1]) == EOF) break;

      ts->size ++;
    }
    fclose (fp);

    if (ts->op) applyderivative (ts);
  }

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
  size_t n;

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
      if (ts->path)
      {
	char buf [4096];
	FILE *f0 = fopen (ts->path, "r");
	ASSERT_TEXT (f0, "Opening data file [%s] for partially cached TIME_SERIES file has failed", ts->path);
        while ((n = fread(buf, sizeof(char), sizeof(buf), f0)) > 0)
        {
          ASSERT_TEXT (fwrite(buf, sizeof(char), n, fp) == n, "Writing data file [%s] has failed", path);
        }
	fclose (f0);
      }
      else
      {
	p = ts->points;
	for (n = 0; n < (size_t) ts->size; n ++)
	  fprintf (fp, "%e\t%e\n", p[n][0], p[n][1]);
      }
    }
  }

  fclose (fp);
}

void TMS_Destroy (TMS *ts)
{
  if (ts->label) return; /* labeled time series are released by TMS_RELEASE_GLOBAL_MAP */

  if (ts->size == 0) free (ts);
  else
  {
    free (ts->points);
    if (ts->path)
    {
      free (ts->path);
      free (ts->offset);
      free (ts->time);
    }
    free (ts);
  }
}

void TMS_Pack (TMS *ts, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  if (ts->label)
  {
    pack_int (isize, i, ints, 1); /* labeled */
    pack_string (isize, i, ints, ts->label);
  }
  else if (ts->path)
  {
    pack_int (isize, i, ints, 2); /* not labeled and partially cached */
    pack_int (isize, i, ints, ts->cache);
    pack_string (isize, i, ints, ts->path);
    pack_int (isize, i, ints, ts->op);
  }
  else
  {
    pack_int (isize, i, ints, 0); /* not labeled and not partially cached */
    pack_double (dsize, d, doubles, ts->value);
    pack_int (isize, i, ints, ts->marker);
    pack_int (isize, i, ints, ts->size);
    if (ts->size) pack_doubles (dsize, d, doubles, (double*)ts->points, ts->size * 2);
  }
}

TMS* TMS_Unpack (int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int type;
  TMS *ts;

  type = unpack_int (ipos, i, ints);

  if (type == 1) /* labeled */
  {
    char *lb = unpack_string (ipos, i, ints);
    ts = MAP_Find (labeled_time_series, lb, (MAP_Compare)strcmp);
    ASSERT_TEXT (ts, "Labeled time series is not present in the global map");
    free (lb);
  }
  else if (type == 2) /* not labeled and partially cached */
  {
    int cache = unpack_int (ipos, i, ints);
    char *path = unpack_string (ipos, i, ints);
    ts = TMS_File (path, NULL, cache);
    ts->op = unpack_int (ipos, i, ints);
    if (ts->op) applyderivative (ts);
    free (path);
  }
  else /* not labeled and not partially cached */
  {
    ERRMEM (ts = malloc (sizeof (TMS)));
    ts->value = unpack_double (dpos, d, doubles);
    ts->marker = unpack_int (ipos, i, ints);
    ts->size = unpack_int (ipos, i, ints);
    if (ts->size)
    {
      ERRMEM (ts->points = malloc (sizeof (double [2]) * ts->size));
      unpack_doubles (dpos, d, doubles, (double*)ts->points, ts->size * 2);
    }
    ts->cache = 0;
    ts->label = NULL;
    ts->path = NULL; /* not partially cached */
    ts->offset = NULL;
    ts->time = NULL;
    ts->noffsets = 0;
    ts->op = 0;
  }

  return ts;
}

/* export MBFCP definition */
void TMS_2_MBFCP (TMS *ts, FILE *out)
{
  int n;

  if (ts->size == 0)
  {
    fprintf (out, "CONSTANT:\t%g\n", ts->value);
  }
  else
  {
    fprintf (out, "TIME_SERIES:\t%d\n", ts->size);
    if (ts->path)
    {
      for (int i = 0; i < ts->noffsets; i ++)
      {
	TMS_Value (ts, ts->time[i]); /* load partial cache at this offset */

	/* output data points from the current cache */
	for (n = 0; n < ts->size; n ++) fprintf (out, "%g\t%g\n", ts->points [n][0], ts->points[n][1]);
      }
    }
    else for (n = 0; n < ts->size; n ++) fprintf (out, "%g\t%g\n", ts->points [n][0], ts->points[n][1]);
  }
}

/* released globally mapped time series */
void TMS_RELEASE_GLOBAL_MAP()
{
  MAP *item;
  TMS *ts;

  for (item = MAP_First (labeled_time_series); item; item = MAP_Next (item))
  {
    ts = item->data;
    ts->label = NULL;
    TMS_Destroy (ts); /* release TMS objects */
  }

  for (item = MAP_First (labeled_time_series); item; item = MAP_Next (item))
  {
    free (item->key); /* release labels */
  }

  MAP_Free (NULL, &labeled_time_series);  /* release map memory */
}
