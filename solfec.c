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

#if ZOLTAN
#include <zoltan.h>
#endif

#include <signal.h>
#include <string.h>
#include <stdio.h>
#include "glv.h"
#include "rnd.h"
#include "lng.h"
#include "sol.h"
#include "err.h"
#include "ext/predicates.h"
#include "version.h"

/* global list of created SOLFEC objects */
static SOLFEC *solfec = NULL;

/* global output directory */
static char *OUTDIR = NULL;

/* global output sub-directory */
static char *SUBDIR = NULL;

/* global input file path */
static char *INPUTFILE = NULL;

/* global write mode flag */
static int WRITEMODEFLAG = 0;

/* global write mode flag */
static int CONTINUEWRITEFLAG = 0;

/* global wireframe rendering flag */
static int WIREFRAMEFLAG = 0;

/* global verbose output time interval */
static double VERBOSITYINTERVAL = 1.0;

/* global non-Solfec input arguments */
#define MAX_ARGC 128
static char *ARGV [MAX_ARGC];
static int ARGC = 0;

/* register new SOLFEC object */
void REGISTER_SOLFEC (SOLFEC *sol)
{
  sol->next = solfec;
  solfec = sol;
}

/* get output directory */
char* OUTPUT_DIR ()
{
  return OUTDIR;
}

/* get output sub-directory */
char* OUTPUT_SUBDIR ()
{
  return SUBDIR;
}

/* get input file path */
char* INPUT_FILE ()
{
  return INPUTFILE;
}

/* get write mode flag */
int WRITE_MODE_FLAG ()
{
  return WRITEMODEFLAG;
}

/* get continue write flag */
int CONTINUE_WRITE_FLAG ()
{
  return CONTINUEWRITEFLAG;
}

/* get wireframe flag */
int WIREFRAME_FLAG ()
{
  return WIREFRAMEFLAG;
}

/* verbose output time interval */
double VERBOSITY_INTERVAL ()
{
  return VERBOSITYINTERVAL;
}

/* get non-Solfec input arguments */
char** NON_SOLFEC_ARGV (int *argc)
{
  if (ARGC) 
  {
    *argc = ARGC;
    return ARGV;
  }
  else
  {
    *argc = 0;
    return NULL;
  }
}

#ifndef LIBSOLFEC /* executables */
#if OPENGL
static int WIDTH = 512, HEIGHT = 512; /* initial width and height of the viewer window */
#endif

#if MPI_VERSION >= 2
/* error handler callback */
static void MPI_error_handling (MPI_Comm *comm, int *err, ...)
{
  if (*err != MPI_SUCCESS)
  {
    if (*err == MPI_ERR_COMM)
    {
      fprintf (stderr, "MPI error: Invalid communicator\n"), fflush (stderr);
    }
    else if (*err == MPI_ERR_OTHER)
    {
      char str [MPI_MAX_ERROR_STRING];
      int len;

      MPI_Error_string (*err, str, &len);
      fprintf (stderr, "MPI error: %s\n", str), fflush (stderr);
    }

    for (; solfec; solfec = solfec->next) SOLFEC_Abort (solfec); /* abort SOLFEC (flush buffers) */

    TMS_RELEASE_GLOBAL_MAP(); /* release globally mapped time series */

    exit (1);
  }
}

/* set up MPI global error handler */
static void MPI_set_error_handling ()
{
  MPI_Errhandler handler;

  MPI_Comm_create_errhandler (MPI_error_handling, &handler);

  MPI_Comm_set_errhandler (MPI_COMM_WORLD, handler);
}
#endif

/* signal handler */
static void sighnd (int signal)
{
  for (; solfec; solfec = solfec->next) SOLFEC_Abort (solfec); /* abort SOLFEC (flush buffers) */

#if MPI
  MPI_Finalize ();
#endif

  TMS_RELEASE_GLOBAL_MAP(); /* release globally mapped time series */

  exit (1);
}

/* return a file path from among the input arguments */
static char* getfile (int argc, char **argv)
{
  char *path;
  FILE *f;
  int n;

  for (n = 1, path = NULL; n < argc; n ++)
  {
    if (strcmp (argv [n], "-s") == 0)
    {
      if (++ n < argc)
      {
	SUBDIR = argv [n]; /* set sub-directory */
      }
    }
    else if (strcmp (argv [n], "-i") == 0)
    {
      if (++ n < argc)
      {
	VERBOSITYINTERVAL = atof (argv[n]); /* verbose output time interval */
      }
    }
#if OPENGL
    else if (strcmp (argv [n], "-g") == 0)
    {
      if (++ n < argc)
      {
	sscanf (argv[n], "%dx%d", &WIDTH, &HEIGHT);
      }
    }
#endif
    else if (strcmp (argv [n], "-v") == 0) continue;
    else if (strcmp (argv [n], "-w") == 0) WRITEMODEFLAG = 1;
    else if (strcmp (argv [n], "-c") == 0) 
    {
      WRITEMODEFLAG = 1;
      CONTINUEWRITEFLAG = 1;
    }
    else if (strcmp (argv [n], "-f") == 0) WIREFRAMEFLAG = 1;
    else if (path == NULL && (f = fopen (argv [n], "r")))
    {
      path = argv [n];
      fclose (f);
    }
    else if (ARGC < MAX_ARGC) /* non-Solfec argument */
    {
      ARGV [ARGC ++] = argv [n];
    }
  }

  if (path) n = strlen (path);
  else n = 0;
  
  if (n > 3)
  {
    if (!(path[n-3] == '.' && path[n-2] == 'p' && path[n-1] == 'y'))
    {
      int m = 0;

      if (path[n-1] == '/') /* /a/directory/path/ */
      {
        path[n-1] = '\0';
	n--;
      }

      while (n > 0 && path[n-1] != '/')
      {
        n--;
	m++;
      }

      char *path1 = malloc (n+m+m+16);

      ERRMEM (path1);

      sprintf (path1, "%s/%s.py", path, &path[n]);

      OUTDIR = path; /* this will overwrite SOLFEC's outpath */

      path = path1; /* XXX: allow for memory leak */
    }
  }

  INPUTFILE = path;

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

  exactinit ();

#if !defined (__MINGW32__)
  signal (SIGHUP, sighnd);
  signal (SIGQUIT, sighnd);
  signal (SIGTSTP, sighnd);
#endif
  signal (SIGINT, sighnd);
  signal (SIGTERM, sighnd);
  signal (SIGSEGV, sighnd);

#if MPI
  MPI_Init (&argc, &argv);
#if ZOLTAN
  float version;
  ASSERT (Zoltan_Initialize (argc, argv, &version) == ZOLTAN_OK, ERR_ZOLTAN_INIT);
#endif

#if MPI_VERSION >= 2
  MPI_set_error_handling ();
#endif
#endif

  TRY ()
  {
    int lngerr = 0;

#if OPENGL
    if (vieweron (argc, argv)) RND_Switch_On (); /* make renderer aware of viewer before calling interpreter */
    #define synopsis "SYNOPSIS: solfec [-v] [-w] [-c] [-f] [-g WIDTHxHEIGHT] [-s sub-directory] [-i verbosity inteval] path\n"
#else
  #if MPI
    #define synopsis "SYNOPSIS: solfec-mpi [-c] [-s sub-directory] [-i verbosity inteval] path\n"
  #else
    #define synopsis "SYNOPSIS: solfec [-w] [-c] [-s sub-directory] [-i verbosity inteval] path\n"
  #endif
#endif

    char *path = getfile (argc, argv); /* parse input */

    if (!path) 
    {
      printf ("VERSION: Solfec-1.%s (%s)\n", VERSION_HASH, VERSION_DATE);
      printf (synopsis); /* print info */
    }
    else lngerr = lng (path); /* call interpreter */

    WARNING (lngerr == 0, "input parsing error");

#if OPENGL
    if (vieweron (argc, argv) && !lngerr)
    {
      double extents [6] = {-1, -1, -1, 1, 1, 1};

      GLV (&argc, argv, "Solfec", WIDTH, HEIGHT, extents, /* run viewer */
	   RND_Menu, RND_Init, RND_Idle, RND_Quit, RND_Render,
	   RND_Key, RND_Keyspec, RND_Mouse, RND_Motion, RND_Passive);

      /* FIXME: viewer should be called from lng_RUN in READ mode (for a SOLFEC object)
       * FIXME: so that any body creations that follow after the RUN command are prevented
       * FIXME: in Python and handled via reading of bodies from time > 0.0 */
    }
#endif

    lngfinalize (); /* finalize interpreter */
  }
  CATCHANY (error)
  {
    fprintf (stderr, "Error: %s\n", errstring (error));
    sighnd (0);
    return 1;
  }
  ENDTRY ()

#if MPI
  MPI_Finalize ();
#endif

  TMS_RELEASE_GLOBAL_MAP(); /* release globally mapped time series */

  return 0;
}
#endif /* ~LIBSOLFEC */
