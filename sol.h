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

#include <time.h>
#include "alg.h"
#include "cvx.h"
#include "msh.h"
#include "sph.h"
#include "tms.h"
#include "bgs.h"
#include "pes.h"
#include "nts.h"
#include "hbs.h"
#include "dom.h"
#include "mat.h"
#include "pbf.h"
#include "cmp.h"

#ifndef __sol__
#define __sol__

#ifndef SOLFEC_TYPE
#define SOLFEC_TYPE
typedef struct solfec SOLFEC;
#endif

typedef int (*SOLFEC_Callback) (SOLFEC*, void*, void*);

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
  int ioparallel; /* parallel output flag */
  CMP_ALG output_compression;
  double output_interval,
	 output_time;
  char *outpath;
  PBF *bf;  

  /* callback data */
  double callback_interval,
	 callback_time;
  void *data, *call;
  SOLFEC_Callback callback;

  /* labaled timers */
  MEM mapmem, timemem;
  MAP *timers;

  /* current run start and duration */
  double t0, duration;
  time_t start;

  /* global verbosity flag */
  short verbose;

  /* list structure */
  SOLFEC *next;
};

/* create a solfec instance */
SOLFEC* SOLFEC_Create (short dynamic, double step, char *outpath);

/* start a labeled timer */
void SOLFEC_Timer_Start (SOLFEC *sol, const char *label);

/* end a labeled timer (labeled timers are written to the output) */
void SOLFEC_Timer_End (SOLFEC *sol, const char *label);

/* get timing of a labeled timer */
double SOLFEC_Timing (SOLFEC *sol, const char *label);

/* test whether a labeled timer exists */
int SOLFEC_Has_Timer (SOLFEC *sol, const char *label);

/* solfec mode string */
char* SOLFEC_Mode (SOLFEC *sol);

/* run analysis with a specific constraint solver */
void SOLFEC_Run (SOLFEC *sol, SOLVER_KIND kind, void *solver, double duration);

/* set results output interval */
void SOLFEC_Output (SOLFEC *sol, double interval, CMP_ALG compression);

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

/* perform abort actions */
void SOLFEC_Abort (SOLFEC *sol);

/* free solfec memory */
void SOLFEC_Destroy (SOLFEC *sol);

/* SOLFEC time history item */
typedef struct solfec_history_item SHI;

/* history item data */
struct solfec_history_item
{
  BODY *bod;
  double point [3];
  short index;
  VALUE_KIND entity; /* (body, point, index, entity) <=> BODY_ENTITY */
  SET *bodies; /* (bodies, index) <=> ENERGY_VALUE */
  char *label; /* label <=> TIMING_VALUE or LABELED_ */
  enum {BODY_ENTITY, ENERGY_VALUE, TIMING_VALUE,
        LABELED_INT, LABELED_DOUBLE} item;
  enum {OP_SUM, OP_AVG, OP_MAX, OP_MIN} op; /* operation for parallel labeled data */
  unsigned char in_domain; /* flag > 0 when labeled data in inside domain data */
  double *history; /* output */
};

/* read histories of a set of requested items; allocate and fill 'history'  members
 * of those items; return table of times of the same 'size' as the 'history' members */
double* SOLFEC_History (SOLFEC *sol, SHI *shi, int nshi, double t0, double t1, int skip, int *size);
#endif
