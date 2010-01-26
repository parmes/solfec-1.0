/*
 * sol.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * solfec main module
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

#if MPI
#include <zoltan.h>
#include <time.h>
#include "put.h"
#endif

#if POSIX
#include <sys/stat.h>
#endif
#include <string.h>
#include <float.h>
#include "solfec.h"
#include "sol.h"
#include "err.h"
#include "tmr.h"

/* ============================= INPUT-OUTPUT VERSIOB ============================ */
/* === Version increments require precise records of causes and affected files === */
/* =============================================================================== */
/* 1                         an initial input-output version (dio.c, sol.c, bod.c) */
/* =============================================================================== */
#define IOVER 1                                                 /* current version */
/* =============================================================================== */

/* defulat initial amoung of boxes */
#define DEFSIZE 1024

/* clean labeled timers (total member of TIMING) */
static void clean_timers (SOLFEC *sol)
{
  for (MAP *item = MAP_First (sol->timers); item; item = MAP_Next (item))
  {
    TIMING *t = item->data;
    t->total = 0.0;
  }
}

/* turn on verbosity */
static int verbose_on (SOLFEC *sol, short kind, void *solver)
{
  if (sol->verbose)
  {
    sol->dom->verbose = 1;
  }

  return sol->verbose;
}

/* turn off verbosity */
static int verbose_off (SOLFEC *sol, short kind, void *solver)
{
  sol->dom->verbose = 0;

  return 0;
}

/* copy output path */
static char* copyoutpath (char *outpath)
{
  char *out = NULL;
  int l, k, m;

  if ((l = outpath ? strlen (outpath) : 0))
  {
    m = 0;
    k = OUTPUT_SUBDIR() ? strlen (OUTPUT_SUBDIR()) : 0;
    if (k && outpath [l-1] != '/') m = 1;
    ERRMEM (out = malloc (l + k + m + 1));
    strcpy (out, outpath);
    if (k)
    {
      if (outpath [l-1] != '/') out [l] = '/';
      strcpy (&out [l+m], OUTPUT_SUBDIR());
    }
  }

  return out;
}

/* from directory path get last name */
static char *lastname (char *path)
{
  int l = strlen (path);

  while (l > 0 && path [l-1] != '/') l --;

  return &path [l];
}

/* get file path from directory path */
static char *getpath (char *outpath)
{
  int l = strlen (outpath),
      n = l + strlen (lastname (outpath)) + 8;
  char *path;

  ERRMEM (path = malloc (n));
  strcpy (path, outpath);
  path [l] = '/';
  strcpy (path+l+1, lastname (outpath));

  return path;
}

/* attempt reading output path */
static PBF* readoutpath (char *outpath)
{
  char *path = getpath (outpath);
  PBF *bf = PBF_Read (path);

  free (path);
  return bf;
}

/* attempt writing output path */
static PBF* writeoutpath (char *outpath)
{
#if POSIX
  int i, l = strlen (outpath);

  for (i = 0; i < l; i ++) /* create all directories on the way */
  {
    if (outpath [i] == '/')
    {
       outpath [i] = '\0';
       mkdir (outpath, 0777); /* POSIX */
       outpath [i] = '/';
    }
  }
  mkdir (outpath, 0777); /* POSIX */
#endif

  char *path = getpath (outpath);
  PBF *bf = PBF_Write (path);

  free (path);
  return bf;
}

/* output state */
static void write_state (SOLFEC *sol)
{
  /* write time */

  PBF_Time (sol->bf, &sol->dom->time); /* the only domain member written outside of it */

  /* write initial flags */

  if (sol->iover < 0)
  {
    sol->iover = -sol->iover; /* make positive */
    PBF_Label (sol->bf, "IOVER");
    PBF_Int (sol->bf, &sol->iover, 1);
    PBF_Label (sol->bf, "IOPARALLEL");
    PBF_Int (sol->bf, &sol->ioparallel, 1);
  }

  /* write domain */

  DOM_Write_State (sol->dom, sol->bf, sol->output_compression);

  /* write timers */

  int numt = MAP_Size (sol->timers);
  PBF_Label (sol->bf, "TIMERS");
  PBF_Int (sol->bf, &numt, 1);
  
  for (MAP *item = MAP_First (sol->timers); item; item = MAP_Next (item))
  {
    TIMING *t = item->data;

    PBF_Label (sol->bf, item->key);
    PBF_Double (sol->bf, &t->total, 1);
    PBF_String (sol->bf, (char**) &item->key);
  }

  clean_timers (sol); /* restart total timing */
}

/* input state */
static void read_state (SOLFEC *sol)
{
  /* read time */

  PBF_Time (sol->bf, &sol->dom->time); /* the only domain member red outside of it */

  /* read initial flags */

  if (sol->iover < 0)
  {
    ASSERT (PBF_Label (sol->bf, "IOVER"), ERR_FILE_FORMAT);
    PBF_Int (sol->bf, &sol->iover, 1);
    ASSERT (PBF_Label (sol->bf, "IOPARALLEL"), ERR_FILE_FORMAT);
    PBF_Int (sol->bf, &sol->ioparallel, 1);
  }

  /* read domain */

  DOM_Read_State (sol->dom, sol->bf);

  /* read timers */

  clean_timers (sol); /* zero total timing */

  int n, numt, found = 0;
  for (PBF *bf = sol->bf; bf; bf = bf->next)
  {
    if (PBF_Label (bf, "TIMERS"))
    {
      found = 1;

      PBF_Int (bf, &numt, 1);
      
      for (n = 0; n < numt; n ++)
      {
	double total;
	char *label;
	TIMING *t;

	PBF_Double (bf, &total, 1); /* read timing */
	label = NULL; /* needs allocation */
	PBF_String (bf, &label); /* read label */

	if (!(t = MAP_Find (sol->timers, label, (MAP_Compare) strcmp))) /* add timer if missing */
	{
	  ERRMEM (t = MEM_Alloc (&sol->timemem));
	  MAP_Insert (&sol->mapmem, &sol->timers, label, t, (MAP_Compare) strcmp);
	}

	t->total = MAX (t->total, total); /* update to maximum across all processors */
      }
    }
  }
  ASSERT (found, ERR_FILE_FORMAT); /* the former root file should have this section */
}

/* read initial state if needed */
static int init (SOLFEC *sol)
{
  if (sol->iover < 0)
  {
    read_state (sol);
    return 0;
  }
  else return 1;
}

/* output statistics */
static void statsout (SOLFEC *sol)
{
  DOM *dom = sol->dom;
  time_t timer = time(NULL);
  double  elapsed = difftime (timer, sol->start), /* elapsed wall clock time for this run */
	  estimated = (elapsed / (sol->dom->time  - sol->t0)) * sol->duration - elapsed; /* estimated remaining time */
  int days = (int) (estimated / 86400.),
      hours = (int) ((estimated - days * 86400.) / 3600.),
      minutes = (int) ((estimated - days * 86400. - hours * 3600.) / 60.),
      seconds = (int) (estimated - days * 86400. - hours * 3600. - minutes * 60.);
  char string [32];

  ctime_r(&timer, string); 

#if MPI
  char *stapath;
  FILE *sta;
  int i;

  if (dom->rank == 0)
  {
    ERRMEM (stapath = malloc (strlen (sol->outpath) + 64));
    sprintf (stapath, "%s/STATE", sol->outpath);
    ASSERT (sta = fopen (stapath, "w"), ERR_FILE_OPEN);
    fprintf (sta, "----------------------------------------------------------------------------------------\n");
    fprintf (sta, "%sEstimated end in %d days, %d hours, %d minutes and %d seconds\n", string, days, hours, minutes, seconds);
    fprintf (sta, "----------------------------------------------------------------------------------------\n");
    fprintf (sta, "TIME: %g\n", sol->dom->time);
    fprintf (sta, "----------------------------------------------------------------------------------------\n");

    printf ("----------------------------------------------------------------------------------------\n");
    for (i = 0; i < dom->nstats; i ++) 
    {
      fprintf (sta, "%13s: SUM = %8d     MIN = %8d     AVG = %8d     MAX = %8d\n", dom->stats [i].name, dom->stats [i].sum, dom->stats [i].min, dom->stats [i].avg, dom->stats [i].max);
      printf ("%13s: SUM = %8d     MIN = %8d     AVG = %8d     MAX = %8d\n", dom->stats [i].name, dom->stats [i].sum, dom->stats [i].min, dom->stats [i].avg, dom->stats [i].max); 
    }
    fprintf (sta, "----------------------------------------------------------------------------------------\n");
    printf ("----------------------------------------------------------------------------------------\n");

    fclose (sta);
    free (stapath);
  }
#else
  const int N = 4;

  char *name [] = {"BODIES", "BOXES", "CONSTRAINTS", "SPARSIFIED"};

  int val [] = {dom->nbod, dom->aabb->boxnum, dom->ncon, dom->nspa}, i;

  for (i = 0; i < N; i ++) printf ("%11s: %8d\n", name [i], val [i]);
  printf ("%sEstimated end in %d days, %d hours, %d minutes and %d seconds\n", string, days, hours, minutes, seconds);
  printf ("----------------------------------------------------------------------------------------\n");
#endif
}

/* invoke constraint solver */
static void SOLVE (SOLFEC *sol, void *solver, SOLVER_KIND kind, LOCDYN *ldy, int verbose)
{
  SOLFEC_Timer_Start (sol, "CONSOL");

#if MPI
  if (ldy->dom->rank == 0)
#endif
  if (verbose) printf ("SOLVER ...\n");

  switch (kind)
  {
  case GAUSS_SEIDEL_SOLVER: GAUSS_SEIDEL_Solve (solver, ldy); break;
  case EXPLICIT_SOLVER: EXPLICIT_Solve (ldy); break;
  }

  SOLFEC_Timer_End (sol, "CONSOL");
}

/* create a solfec instance */
SOLFEC* SOLFEC_Create (short dynamic, double step, char *outpath)
{
  SOLFEC *sol;

  ERRMEM (sol = malloc (sizeof (SOLFEC)));
  sol->aabb = AABB_Create (DEFSIZE);
  sol->sps = SPSET_Create ();
  sol->mat = MATSET_Create ();
  sol->dom = DOM_Create (sol->aabb, sol->sps, dynamic, step);
  sol->dom->solfec = sol;

  sol->outpath = copyoutpath (outpath);
  sol->output_interval = 0;
  sol->output_time = 0;
  sol->output_compression = CMP_OFF;
  if ((sol->bf = readoutpath (sol->outpath))) sol->mode = SOLFEC_READ;
  else if ((sol->bf = writeoutpath (sol->outpath))) sol->mode = SOLFEC_WRITE;
  else THROW (ERR_FILE_OPEN);
  sol->iover = -IOVER; /* negative to indicate initial state */
#if MPI
  sol->ioparallel = 1; /* some data specific to parallel run will be outputed */
#else
  sol->ioparallel = 0;
#endif

  sol->callback_interval = DBL_MAX;
  sol->callback_time = DBL_MAX;
  sol->data = sol->call = NULL;
  sol->callback = NULL;

  MEM_Init (&sol->mapmem, sizeof (MAP), 128);
  MEM_Init (&sol->timemem, sizeof (TIMING), 128);
  sol->timers = NULL;
  sol->verbose = 1;

  return sol;
}

/* start a labeled timer */
void SOLFEC_Timer_Start (SOLFEC *sol, const char *label)
{
  TIMING *t;

  if (!(t = MAP_Find (sol->timers, (void*) label, (MAP_Compare) strcmp)))
  {
    ERRMEM (t = MEM_Alloc (&sol->timemem));
    MAP_Insert (&sol->mapmem, &sol->timers, (void*) label, t, (MAP_Compare) strcmp);
  }

  timerstart (t);
}

/* end a labeled timer (labeled timers are written to the output) */
void SOLFEC_Timer_End (SOLFEC *sol, const char *label)
{
  TIMING *t;

  if ((t = MAP_Find (sol->timers, (void*) label, (MAP_Compare) strcmp)))
  {
    timerend (t);
  }
}

/* get timing of a labeled timer */
double SOLFEC_Timing (SOLFEC *sol, const char *label)
{
  TIMING *t;

  if ((t = MAP_Find (sol->timers, (void*) label, (MAP_Compare) strcmp))) return t->total;

  return 0.0;
}

/* test whether a labeled timer exists */
int SOLFEC_Has_Timer (SOLFEC *sol, const char *label)
{
  if (MAP_Find (sol->timers, (void*) label, (MAP_Compare) strcmp)) return 1;
  else return 0;
}

/* solfec mode string */
char* SOLFEC_Mode (SOLFEC *sol)
{
  switch (sol->mode)
  {
  case SOLFEC_WRITE: return "WRITE";
  case SOLFEC_READ: return "READ";
  }

  return NULL;
}

/* run analysis with a specific constraint solver */
void SOLFEC_Run (SOLFEC *sol, SOLVER_KIND kind, void *solver, double duration)
{
  if (sol->mode == SOLFEC_WRITE)
  {
    UPKIND upkind;
    LOCDYN *ldy;
    int verbose;
    TIMING tim;
    double tt;

    verbose = verbose_on (sol, kind, solver);
    sol->duration = duration;
    sol->start = time (NULL);
    timerstart (&tim);

    for (sol->t0 = sol->dom->time; sol->dom->time < (sol->t0 + duration);)
    {
      upkind = (kind == EXPLICIT_SOLVER ? UPEXS : UPALL); /* here as user callback can change solver */

#if MPI
      if (sol->dom->rank == 0)
#endif
      if (verbose) printf ("TIME: %g ... ", sol->dom->time);

      /* begin update of domain */
      ldy = DOM_Update_Begin (sol->dom);

      /* begin update of local dynamics */
      LOCDYN_Update_Begin (ldy, upkind);

      /* solve constraints */
      SOLVE (sol, solver, kind, ldy, verbose);

      /* end update of local dynamics */
      LOCDYN_Update_End (ldy, upkind);

      /* end update of domain */
      DOM_Update_End (sol->dom);

      /* write output if needed */
      if (sol->dom->time >= sol->output_time)
      {
	sol->output_time += sol->output_interval;
	write_state (sol);
      }

      /* execute callback if needed */
      if (sol->dom->time >= sol->callback_time)
      {
	int ret;

	sol->callback_time += sol->callback_interval;
	ret = sol->callback (sol, sol->data, sol->call);
#if MPI
	ret = PUT_int_min (ret);
#endif
	if (!ret) break; /* interrupt run */
      }

      /* statistics are printed every
       * human perciveable period of time */
      tt = timerend (&tim);
      if (verbose) statsout (sol);
      if (tt < 1.0) verbose = verbose_off (sol, kind, solver);
      else if (tt >= 1.0) verbose = verbose_on (sol, kind, solver), timerstart (&tim);
    }
  }
  else /* READ */
  {
    if (init (sol))
    {
      PBF_Forward (sol->bf, 1);
      read_state (sol);
    }
  }
}

/* set results output interval */
void SOLFEC_Output (SOLFEC *sol, double interval, CMP_ALG compression)
{
  sol->output_interval = interval;
  sol->output_time = sol->dom->time + interval;
  sol->output_compression = compression;
}

/* the next time minus the current time */
double SOLFEC_Time_Skip (SOLFEC *sol)
{
  if (sol->mode == SOLFEC_READ)
  {
    double t;

    init (sol);
    PBF_Forward (sol->bf, 1);
    PBF_Time (sol->bf, &t);
    PBF_Backward (sol->bf, 1);
    return t - sol->dom->time;
  }
  else return sol->dom->step;
}

/* get analysis duration time limits */
void SOLFEC_Time_Limits (SOLFEC *sol, double *start, double *end)
{
  if (sol->mode == SOLFEC_READ)
  {
    double s, e;

    init (sol);
    PBF_Limits (sol->bf, &s, &e);
    if (start) *start = s;
    if (end) *end = e;
  }
  else 
  {
    if (start) *start = 0;
    if (end) *end = sol->dom->time;
  }
}

/* set up callback function */
void SOLFEC_Set_Callback (SOLFEC *sol, double interval, void *data, void *call, SOLFEC_Callback callback)
{
  sol->callback_interval = interval;
  sol->callback_time = sol->dom->time + interval;
  sol->data = data;
  sol->call = call;
  sol->callback = callback;
}

/* seek to specific time in READ mode */
void SOLFEC_Seek_To (SOLFEC *sol, double time)
{
  if (sol->mode == SOLFEC_READ)
  {
    init (sol);
    PBF_Seek (sol->bf, time);
    read_state (sol);
  }
}

/* step backward in READ modes */
void SOLFEC_Backward (SOLFEC *sol, int steps)
{
  if (sol->mode == SOLFEC_READ)
  {
    if (init (sol))
    {
      PBF_Backward (sol->bf, steps);
      read_state (sol);
    }
  }
}

/* step forward in READ modes */
void SOLFEC_Forward (SOLFEC *sol, int steps)
{
  if (sol->mode == SOLFEC_READ)
  {
    if (init (sol))
    {
      PBF_Forward (sol->bf, steps);
      read_state (sol);
    }
  }
}

/* read the history of an object (a labeled value, a body or
 * a constraint) and invoke the callback for every new state */
/* perform abort actions */
void SOLFEC_Abort (SOLFEC *sol)
{
  write_state (sol);

  if (sol->bf)
  {
    PBF_Flush (sol->bf);
    PBF_Close (sol->bf);
    sol->bf = NULL;
  }
}

/* free solfec memory */
void SOLFEC_Destroy (SOLFEC *sol)
{
  AABB_Destroy (sol->aabb); 
  SPSET_Destroy (sol->sps);
  MATSET_Destroy (sol->mat);
  DOM_Destroy (sol->dom);

  free (sol->outpath);

  if (sol->bf) PBF_Close (sol->bf);

  if (sol->mode == SOLFEC_READ)
  {
    for (MAP *item = MAP_First (sol->timers); item; item = MAP_Next (item))
    {
      free (item->key); /* labels were allocated in READ mode */
    }
  }

  MEM_Release (&sol->mapmem);
  MEM_Release (&sol->timemem);

  free (sol);
}

/* read histories of a set of requested items; allocate and fill 'history'  members
 * of those items; return table of times of the same 'size' as the 'history' members */
double* SOLFEC_History (SOLFEC *sol, SHI *shi, int nshi, double t0, double t1, int skip, int *size)
{
  if (sol->mode == SOLFEC_WRITE) return NULL;

  double save, *time;
  int cur, i;

  cur = 0;
  save = sol->dom->time;
  SOLFEC_Seek_To (sol, t0);
  *size = PBF_Span (sol->bf, t0, t1);
  ERRMEM (time = MEM_CALLOC (sizeof (double [(*size) + 4]))); /* safeguard */
  time [cur] = sol->dom->time;

  for (i = 0; i < nshi; i ++)
  {
    ERRMEM (shi[i].history = MEM_CALLOC (sizeof (double [(*size) + 4])));
  }

  do
  {
    for (i = 0; i < nshi; i ++)
    {
      switch (shi[i].item)
      {
      case BODY_ENTITY:
	{
	  double values [7];
          BODY_Point_Values (shi[i].bod, shi[i].point, shi[i].entity, values);
	  shi[i].history [cur] = values [shi[i].index];
	}
	break;
      case ENERGY_VALUE:
	{
	  for (SET *item = SET_First (shi[i].bodies); item; item = SET_Next (item))
	  {
	    BODY *bod = item->data;
	    shi[i].history [cur] += bod->energy [shi[i].index];
	  }
	}
	break;
      case TIMING_VALUE:
	{
          shi[i].history [cur] = SOLFEC_Timing (sol, shi[i].label);
	}
	break;
      }
    }

    SOLFEC_Forward (sol, skip);
    time [cur ++] = sol->dom->time;
  }
  while (sol->dom->time < t1);

  SOLFEC_Seek_To (sol, save); /* restore initial time frame */

  return time;
}
