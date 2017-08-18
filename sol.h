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
#include "sis.h"
#include "tts.h"
#include "dom.h"
#include "fld.h"
#include "mat.h"
#include "pbf.h"
#include "cmp.h"
#include "tmr.h"
#include "bcd.h"

#ifndef __sol__
#define __sol__

/* ============================= INPUT-OUTPUT VERSION ============================ */
/* === Version increments require precise records of causes and affected files === */
/* =============================================================================== */
/* 1 |                       an initial input-output version (dio.c, sol.c, bod.c) */
/* ------------------------------------------------------------------------------- */
/* 2 |            output of relative pre-impact velocities for contact constraints */
/* ------------------------------------------------------------------------------- */
/* 3 |  output of iover, body kind, conf and dofs size to enable RIGID_TO_FEM and  */
/*   | INITIALISE_STATE working when LOCAL_BODIES = yes is used in Config.mak      */
/* ------------------------------------------------------------------------------- */
/* 4 | con->Z written out also for the SPRING constraint in dio.c:write_constraint */
/* =============================================================================== */
#define IOVER 4                                                 /* current version */
/* =============================================================================== */

#ifndef FE_BASE_TYPE
#define FE_BASE_TYPE
typedef struct fe_base FE_BASE;
#endif

struct fe_base
{
  MX *evec;
  double *eval;
  char *label;
};

#ifndef SOLFEC_TYPE
#define SOLFEC_TYPE
typedef struct solfec SOLFEC;
#endif

typedef int (*SOLFEC_Callback) (SOLFEC*, void*, void*);

enum solver_kind
{
  NONE_SOLVER          = 0x00,
  GAUSS_SEIDEL_SOLVER  = 0x01,
  PENALTY_SOLVER       = 0x02,
  NEWTON_SOLVER        = 0x04,
  SICONOS_SOLVER       = 0x08,
  TEST_SOLVER          = 0x10,
  HYBRID_SOLVER_KIND   = 0x20
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

  FISET *fis; /* field set */

  SPSET *sps; /* surface pairs and materials */

  MATSET *mat; /* bulk materials */

  DOM *dom; /* bodies, constraints and time integration */

  int iover; /* input-output version */
  double output_interval,
	 output_time;
  char *outpath;
  PBF_FLG compression;
  PBF *bf;  

  /* body co-rotated FEM displacements sampling */
  BCD *bcd;

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
  TIMING verbose_timing;
  short verbose;

  /* output directory removal flag */
  short cleanup;

  /* current solver */
  SOLVER_KIND kind;
  void *solver;

  /* registered FE bases */
  MAP *registered_bases;

  /* list structure */
  SOLFEC *next;
};

/* create a solfec instance */
SOLFEC* SOLFEC_Create (short dynamic, double step, char *outpath);

/* allocate file name without extension */
char* SOLFEC_Alloc_File_Name (SOLFEC *sol, int extlen);

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
void SOLFEC_Output (SOLFEC *sol, double interval, PBF_FLG compression);

/* set up callback function */
void SOLFEC_Set_Callback (SOLFEC *sol, double interval, void *data, void *call, SOLFEC_Callback callback);

/* the next time minus the current time */
double SOLFEC_Time_Skip (SOLFEC *sol);

/* initialize solfec in READ state */
void SOLFEC_Read_Init (SOLFEC *sol);

/* get analysis duration time limits */
void SOLFEC_Time_Limits (SOLFEC *sol, double *start, double *end);

/* seek to specific time in READ mode */
void SOLFEC_Seek_To (SOLFEC *sol, double time);

/* step backward in READ modes */
int SOLFEC_Backward (SOLFEC *sol, int steps);

/* step forward in READ modes */
int SOLFEC_Forward (SOLFEC *sol, int steps, int corotated_displacements);

/* perform abort actions */
void SOLFEC_Abort (SOLFEC *sol);

/* free solfec memory */
void SOLFEC_Destroy (SOLFEC *sol);

/* SOLFEC time history item */
typedef struct solfec_history_item SHI;

/* auxiliary 'index' values */
enum {CONSTRAINT_GAP, CONSTRAINT_R, CONSTRAINT_U};

/* history item data */
struct solfec_history_item
{
  BODY *bod;
  short index;
  double point [3], vector [3];
  int surf1, surf2; /* surface pair */
  unsigned char contacts_only; /* contact filtering flag */
  VALUE_KIND entity; /* (body, point, index, entity) <=> BODY_ENTITY */
  SET *bodies; /* (bodies, index) <=> ENERGY_VALUE, CONSTRAINT_VALUE */
  char *label; /* label <=> TIMING_VALUE or LABELED_ */
  enum {BODY_ENTITY, ENERGY_VALUE, TIMING_VALUE,
        CONSTRAINT_VALUE, LABELED_INT, LABELED_DOUBLE} item;
  enum {OP_SUM, OP_AVG, OP_MAX, OP_MIN} op; /* operation for parallel labeled data */
  double *history; /* output */
};

/* read histories of a set of requested items; allocate and fill 'history'  members
 * of those items; return table of times of the same 'size' as the 'history' members;
 * skip every 'skip' steps; if 'skip' < 0 then print out a percentage based progress bar */
double* SOLFEC_History (SOLFEC *sol, SHI *shi, int nshi, double t0, double t1, int skip, int *size);

/* export MBFCP definition */
void SOLFEC_2_MBFCP (SOLFEC *sol, FILE *out);

/* initialize state from the ouput; return 1 on success, 0 otherwise;
 * optinal "subset" is a set of strings defining POSIX regular expressions to be matched
 * against body labels -- narrowing down the set of bodies whose state will  be initialized */
int SOLFEC_Initialize_State (SOLFEC *sol, char *path, double time, SET *subset);

/* map rigid motion onto FEM bodies; return 1 on success, 0 otherwise;
 * optinal "subset" is a set of strings defining POSIX regular expressions to be matched
 * against body labels -- narrowing down the set of bodies whose state will  be initialized */
int SOLFEC_Rigid_To_FEM (SOLFEC *sol, char *path, double time, SET *subset);

/* register FE base; (evec, eval, label) must be dynamically allocated */
void SOLFEC_Register_Base (SOLFEC *sol, MX *evec, double *eval, char *label);

#endif
