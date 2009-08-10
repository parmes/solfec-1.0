/*
 * sol.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * solfec type
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

#include "alg.h"
#include "cvx.h"
#include "msh.h"
#include "sph.h"
#include "tms.h"
#include "bgs.h"
#include "exs.h"
#include "dom.h"
#include "mat.h"
#include "pbf.h"

#ifndef __sol__
#define __sol__

typedef struct solfec SOLFEC;
typedef int (*SOLFEC_Callback) (SOLFEC*, void*, void*);

enum solver_kind
{
  GAUSS_SEIDEL_SOLVER,
  EXPLICIT_SOLVER,
};

typedef enum solver_kind SOLVER_KIND;

enum solfec_mode
{
  SOLFEC_WRITE,
  SOLFEC_READ
};

typedef enum solfec_mode SOLFEC_MODE;

struct solfec
{
  SOLFEC_MODE mode;

  AABB *aabb; /* contact detection solver */

  SPSET *sps; /* surface pairs and materials */

  MATSET *mat; /* bulk materials */

  DOM *dom; /* bodies, constraints and time integration */

  int iover; /* input-output version */
  double output_interval,
	 output_time;
  char *outpath;
  PBF *bf;  

  /* callback data */
  double callback_interval,
	 callback_time;
  void *data, *call;
  SOLFEC_Callback callback;
};

/* create a solfec instance */
SOLFEC* SOLFEC_Create (short dynamic, double step, char *outpath);

/* solfec mode string */
char* SOLFEC_Mode (SOLFEC *sol);

/* run analysis with a specific constraint solver */
void SOLFEC_Run (SOLFEC *sol, SOLVER_KIND kind, void *solver, double duration);

/* set results output interval */
void SOLFEC_Output (SOLFEC *sol, double interval);

/* set up callback function */
void SOLFEC_Set_Callback (SOLFEC *sol, double interval, void *data, void *call, SOLFEC_Callback callback);

/* the next time minus the current time */
double SOLFEC_Time_Skip (SOLFEC *sol);

/* get analysis duration time limits */
void SOLFEC_Time_Limits (SOLFEC *sol, double *start, double *end);

/* seek to specific time in READ mode */
void SOLFEC_Seek_To (SOLFEC *sol, double time);

/* step backward in READ modes */
void SOLFEC_Backward (SOLFEC *sol, int steps);

/* step forward in READ modes */
void SOLFEC_Forward (SOLFEC *sol, int steps);

/* free solfec memory */
void SOLFEC_Destroy (SOLFEC *sol);

#endif
