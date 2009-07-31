/*
 * solfec.c
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
#endif

#include <sys/stat.h> /* POSIX */
#include <string.h>
#include <float.h>
#include "solfec.h"
#include "lng.h"
#include "err.h"
#include "glv.h"
#include "rnd.h"
#include "tmr.h"

/* ------------------------------------- */
#define IOVER 1 /* input-output version */
/* ----------------------------------- */

/* defulat initial amoung of boxes */
#define DEFSIZE 1024

/* turn on verbosity */
static int verbose_on (SOLFEC *sol, short kind, void *solver)
{
  sol->dom->verbose = 1;

  switch (kind)
  {
  case GAUSS_SEIDEL_SOLVER:
  {
    GAUSS_SEIDEL *gs = (GAUSS_SEIDEL*)solver;
    gs->verbose = 1;
  }
  break;
  case EXPLICIT_SOLVER:
  break;
  }

  return 1;
}

/* turn off verbosity */
static int verbose_off (SOLFEC *sol, short kind, void *solver)
{
  sol->dom->verbose = 0;

  switch (kind)
  {
  case GAUSS_SEIDEL_SOLVER:
  {
    GAUSS_SEIDEL *gs = (GAUSS_SEIDEL*)solver;
    gs->verbose = 0;
  }
  break;
  case EXPLICIT_SOLVER:
  break;
  }

  return 0;
}

/* copy output path */
static char* copyoutpath (char *outpath)
{
  char *out = NULL;
  int l;

  if ((l = outpath ? strlen (outpath) : 0))
  {
    ERRMEM (out = malloc (l + 1));
    strcpy (out, outpath);
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
  mkdir (outpath, 0777); /* POSIX */

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

  /* write version */

  if (sol->iover < 0)
  {
    sol->iover = -sol->iover; /* make positive */
    PBF_Label (sol->bf, "IOVER");
    PBF_Int (sol->bf, &sol->iover, 1);
  }

  /* write domain */

  DOM_Write_State (sol->dom, sol->bf);
}

/* input state */
static void read_state (SOLFEC *sol)
{
  /* read time */

  PBF_Time (sol->bf, &sol->dom->time); /* the only domain member red outside of it */

  /* read version */

  if (sol->iover < 0)
  {
    ASSERT (PBF_Label (sol->bf, "IOVER"), ERR_FILE_FORMAT);
    PBF_Int (sol->bf, &sol->iover, 1);
  }

  /* read domain */

  DOM_Read_State (sol->dom, sol->bf);
}

/* read initial state if needed */
static void init (SOLFEC *sol)
{
  if (sol->iover < 0) read_state (sol);
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
  sol->dom->owner = sol;

  sol->outpath = copyoutpath (outpath);
  sol->output_interval = step;
  sol->output_time = 0;
  if ((sol->bf = readoutpath (outpath))) sol->mode = SOLFEC_READ;
  else if ((sol->bf = writeoutpath (outpath))) sol->mode = SOLFEC_WRITE;
  else THROW (ERR_FILE_OPEN);
  sol->iover = -IOVER; /* negative to indicate initial state */

  sol->callback_interval = DBL_MAX;
  sol->callback_time = DBL_MAX;
  sol->data = sol->call = NULL;
  sol->callback = NULL;

  return sol;
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

#if MPI
void SOLFEC_Run (SOLFEC *sol, SOLVER_KIND kind, void *solver, double duration)
{
  int rank;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  verbose_on (sol, kind, solver);

  printf ("Temporary fake run of process %d\n", rank);

  if (sol->dom->time == 0.0) write_state (sol); /* write zero state */

  balance (sol->dom);

  verbose_off (sol, kind, solver);
}
#else
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

    if (sol->dom->time == 0.0) write_state (sol); /* write zero state */

    verbose = verbose_on (sol, kind, solver);
    timerstart (&tim);

    for (double t0 = sol->dom->time; sol->dom->time < (t0 + duration);)
    {
      upkind = (kind == EXPLICIT_SOLVER ? UPDIA : UPALL);

      /* output time */
      if (verbose) printf ("TIME: %.2e\n", sol->dom->time);

      /* begin update of domain */
      ldy = DOM_Update_Begin (sol->dom);

      /* begin update of local dynamics */
      LOCDYN_Update_Begin (ldy, upkind);

      /* solve constraints */
      switch (kind)
      {
      case GAUSS_SEIDEL_SOLVER: GAUSS_SEIDEL_Solve (solver, ldy); break;
      case EXPLICIT_SOLVER: EXPLICIT_Solve (ldy); break;
      }

      /* end update of local dynamics */
      LOCDYN_Update_End (ldy);

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
	sol->callback_time += sol->callback_interval;
	if (!sol->callback (sol, sol->data, sol->call)) break; /* interrupt run */
      }

      /* statistics are printed every human perciveable period of time */
      if ((tt = timerend (&tim)) < 0.5 && verbose) verbose = verbose_off (sol, kind, solver);
      else if (tt >= 0.5) { verbose = verbose_on (sol, kind, solver); timerstart (&tim); }
    }
  }
  else /* READ */
  {
    read_state (sol);
    PBF_Forward (sol->bf, 1);
  }
}
#endif

/* set results output interval */
void SOLFEC_Output (SOLFEC *sol, double interval)
{
  sol->output_interval = interval;
  sol->output_time = sol->dom->time + interval;
}

/* the next time minus the current time */
double SOLFEC_Time_Skip (SOLFEC *sol)
{
  if (sol->mode == SOLFEC_READ)
  {
    double t;

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
    init (sol);
    PBF_Backward (sol->bf, steps);
    read_state (sol);
  }
}

/* step forward in READ modes */
void SOLFEC_Forward (SOLFEC *sol, int steps)
{
  if (sol->mode == SOLFEC_READ)
  {
    init (sol);
    PBF_Forward (sol->bf, steps);
    read_state (sol);
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

  free (sol);
}

/*
 * main module
 */

/* return a file path from among the input arguments */
static char* getfile (int argc, char **argv)
{
  char *path;
  FILE *f;
  int n;

  for (n = 1, f = NULL; n < argc; n ++)
  {
    if ((f = fopen (argv [n], "r"))) break;
  }

  if (f)
  {
    fclose (f);
    path = argv [n];
  }
  else path = NULL;

  return path;
}

#if OPENGL
/* check whether the viewer option was set */
static int vieweron (int argc, char **argv)
{
  int n;

  for (n = 1; n < argc; n ++)
  {
    if (strcmp (argv [n], "-v") == 0) return 1;
  }

  return 0;
}
#endif

int main (int argc, char **argv)
{
  int error;

#if MPI
  float version;

  ASSERT (Zoltan_Initialize (argc, argv, &version) == ZOLTAN_OK, ERR_ZOLTAN_INIT);
#endif

  TRY ()
  {
    int lngerr = 1;

#if OPENGL
    if (vieweron (argc, argv)) RND_Viewer_On (); /* make renderer aware of viewer before calling interpreter */
    char *synopsis = "SYNOPSIS: solfec [-v] path\n";
#else
    char *synopsis = "SYNOPSIS: solfec path\n";
#endif

    if (!getfile (argc, argv)) printf (synopsis);
    else lngerr = lng (getfile (argc, argv)); /* call interpreter */

#if OPENGL
    if (vieweron (argc, argv) && !lngerr)
    {
      double extents [6] = {-1, -1, -1, 1, 1, 1};

      GLV (&argc, argv, "Solfec", 500, 500, extents, /* run viewer */
	   RND_Menu, RND_Init, RND_Idle, RND_Quit, RND_Render,
	   RND_Key, RND_Keyspec, RND_Mouse, RND_Motion, RND_Passive);
    }
#endif

    lngfinalize (); /* finalize interpreter */
  }
  CATCHANY (error)
  {
    fprintf (stderr, "Error: %s\n", errstring (error));
    return 1;
  }
  ENDTRY ()

  return 0;
}
