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

#include <signal.h>
#include <string.h>
#include <stdio.h>
#include "glv.h"
#include "rnd.h"
#include "lng.h"
#include "sol.h"
#include "err.h"

/* global list of created SOLFEC objects */
static SOLFEC *solfec = NULL;

/* register new SOLFEC object */
void REGISTER_SOLFEC (SOLFEC *sol)
{
  sol->next = solfec;
  solfec = sol;
}

#ifndef LIBSOLFEC /* executables */
#if MPI_VERSION >= 2
/* error handler callback */
static void MPI_error_handling (MPI_Comm *comm, int *arg, ...)
{
  lngfinalize (); /* finalize interpreter: as a result
		     flush output buffers if possible */

  MPI_Abort (MPI_COMM_WORLD, 0);

  exit (1);
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

  lngfinalize (); /* finalize interpreter */
#if MPI
  MPI_Finalize ();
#endif

  exit (1);
}

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

  signal (SIGHUP, sighnd);
  signal (SIGINT, sighnd);
  signal (SIGQUIT, sighnd);
  signal (SIGTSTP, sighnd);
  signal (SIGTERM, sighnd);
  signal (SIGSEGV, sighnd);

#if MPI
  float version;

  MPI_Init (&argc, &argv);
  ASSERT (Zoltan_Initialize (argc, argv, &version) == ZOLTAN_OK, ERR_ZOLTAN_INIT);

#if MPI_VERSION >= 2
  MPI_set_error_handling ();
#endif
#endif

  TRY ()
  {
    int lngerr = 1;

#if OPENGL
    if (vieweron (argc, argv)) RND_Switch_On (); /* make renderer aware of viewer before calling interpreter */
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

      GLV (&argc, argv, "Solfec", 512, 512, extents, /* run viewer */
	   RND_Menu, RND_Init, RND_Idle, RND_Quit, RND_Render,
	   RND_Key, RND_Keyspec, RND_Mouse, RND_Motion, RND_Passive);
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

  return 0;
}
#endif /* ~LIBSOLFEC */
