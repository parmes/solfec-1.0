/* HYBRID_SOLVER implementation */

/*
The MIT License (MIT)

Copyright (c) 2016 EDF Energy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* Contributors: Tomasz Koziara */

#include <complex.h>
#include <stdlib.h>
#include <float.h>
#include "solfec.h"
#include "hys.hpp" /* :) */
#include "hys.h"
#include "sol.h"
#include "dom.h"
#include "alg.h"
#include "err.h"

#define TWO_SIDED_FORCE 0

/* initialize boudary */
static void init_boundary (HYBRID_SOLVER *hs, DOM *dom)
{
  for (MAP *item = MAP_First(hs->solfec2parmec); item; item = MAP_Next (item))
  {
    BODY *bod = MAP_Find (dom->idb, item->key, NULL);
    ASSERT_TEXT (bod, "ERROR: Solfec-Parmec boundary body with Solfec id = %d has not been found", (int)(long) item->key);
#if TWO_SIDED_FORCE == 0
    ASSERT_TEXT (bod->kind == OBS, "ERROR: Solfec-Parmec boundary body with Solfec id = %d is not an obstacle", (int)(long) item->key);
#else
    ASSERT_TEXT (bod->kind == RIG, "ERROR: Solfec-Parmec boundary body with Solfec id = %d is not rigid", (int)(long) item->key);
    if (!bod->parmec) ERRMEM (bod->parmec = MEM_CALLOC (sizeof(PARMEC_FORCE)));
#endif
  }
}

/* apply constant boundary force from Solfec;
 * perform a number of Parmec time integration steps;
 * read back average velocity of boundary bodies into Solfec */
static void parmec_steps (HYBRID_SOLVER *hs, DOM *dom, double step, int nstep)
{
  MAP *item;

  for (item = MAP_First(hs->solfec2parmec); item; item = MAP_Next (item))
  {
    double force[3], torque[3];
    BODY *bod = MAP_Find (dom->idb, item->key, NULL);
    BODY_Rigid_Force (bod, dom->time, dom->step, force, torque);
    int num = (int) (long) item->data; /* parmec particle number */
    parmec_set_force_and_torque (num, force, torque);
#if TWO_SIDED_FORCE == 0
    SET6 (bod->velo, 0.0);
#endif
  }

  for (int i = 0; i < nstep; i ++)
  {
    parmec_one_step (step, hs->parmec_interval, hs->parmec_prefix);

#if TWO_SIDED_FORCE == 0
    for (item = MAP_First(hs->solfec2parmec); item; item = MAP_Next (item))
    {
      BODY *bod = MAP_Find (dom->idb, item->key, NULL);
      int num = (int) (long) item->data; /* parmec particle number */
      double angular[3], linear[3];
      parmec_get_angular_and_linear (num, angular, linear);
      ACC (angular, bod->velo);
      ACC (linear, bod->velo+3);
    }
#endif
  }

#if TWO_SIDED_FORCE
  for (item = MAP_First(hs->solfec2parmec); item; item = MAP_Next (item))
  {
    BODY *bod = MAP_Find (dom->idb, item->key, NULL);
    int num = (int) (long) item->data; /* parmec particle number */
    parmec_get_force_and_torque (num, nstep, bod->parmec->force, bod->parmec->torque);
  }
#else
  double inv = 1.0/(double)nstep;

  for (item = MAP_First(hs->solfec2parmec); item; item = MAP_Next (item))
  {
    BODY *bod = MAP_Find (dom->idb, item->key, NULL);
    SCALE6 (bod->velo, inv);
  }
#endif
}

/* create solver */
HYBRID_SOLVER* HYBRID_SOLVER_Create (char *parmec_file, double parmec_step, double parmec_interval[2],
                 char *parmec_prefix, MAP *parmec2solfec, void *solfec_solver, int solfec_solver_kind)
{
  HYBRID_SOLVER *hs;

  ERRMEM (hs = MEM_CALLOC (sizeof (HYBRID_SOLVER)));
  hs->parmec_file = parmec_file;
  hs->parmec_step = parmec_step;
  if (parmec_interval)
  {
    ERRMEM (hs->parmec_interval = malloc (2 * sizeof (double)));
    hs->parmec_interval[0] = parmec_interval[0];
    hs->parmec_interval[1] = parmec_interval[1];
  }
  else hs->parmec_interval = NULL;
  hs->parmec_prefix = parmec_prefix;
  hs->parmec2solfec = parmec2solfec;
  hs->solfec2parmec = NULL;
  for (MAP *item = MAP_First(parmec2solfec); item; item = MAP_Next (item))
  {
    MAP_Insert (NULL, &hs->solfec2parmec, item->data, item->key, NULL);
  }
  hs->solfec_solver = solfec_solver;
  hs->solfec_solver_kind = solfec_solver_kind;

  parmec_init(hs->parmec_file);

  /* TODO: perhaps invoke Solfec's python init so that Parmec's commands
           do not overwrite the Solfec's ones */

  return hs;
}

/* run solver */
void HYBRID_SOLVER_Run (HYBRID_SOLVER *hs, SOLFEC *sol, double duration)
{
#if HDF5
  if (CONTINUE_WRITE_FLAG())
  {
    fprintf (stderr, "ERROR: continued analysis [-c] is disabled with HYBRID_SOLVER\n");
    exit (1);

    /* TODO: see area in SOLFEC_Run for enabling this (use of 'duration' there and here);
             this also would need to be enabled on the parmec library level */
  }
#endif

  /* find Parmec time step */
  double actual_parmec_step = 0.5*sol->dom->step;
  int num_parmec_steps = 2;
  while (actual_parmec_step > hs->parmec_step)
  {
    actual_parmec_step *= 0.5;
    num_parmec_steps *= 2;
  }

  /* make sure that mapped Solfec bodies are skipped during time integration */
#if TWO_SIDED_FORCE == 0
  sol->dom->skipbodies = hs->solfec2parmec;
#endif

  /* exclude contact detection between mapped Solfec bodies (boundary);
   * we do want contact detection between boundary and interior bodies */
  for (MAP *item = MAP_First (hs->solfec2parmec); item; item = MAP_Next (item))
  {
    for (MAP *jtem = MAP_First (hs->solfec2parmec); jtem; jtem = MAP_Next (jtem))
    {
      if (item->key != jtem->key)
      {
	AABB_Exclude_Body_Pair (sol->aabb, (int)(long)item->key, (int)(long)jtem->key);
      }
    }
  }

  
  if (sol->dom->time == 0.0)
  {
    /* initialize boudary */
    init_boundary (hs, sol->dom);

    /* initial half-step */
    parmec_steps (hs, sol->dom, actual_parmec_step, num_parmec_steps/2);
  }

  double time0 = sol->dom->time;

  /* integrate over duration */
  sol->t0 = sol->dom->time; /* these can be overwritten when negative duration is used */
  sol->duration = duration; /* in SOLFEC_Run, allowing for a correct elapsed time estimate */
  while (sol->dom->time < time0 + duration)
  {
#if TWO_SIDED_FORCE == 0
    for (MAP *item = MAP_First(hs->solfec2parmec); item; item = MAP_Next (item))
    {
      BODY *bod = MAP_Find (sol->dom->idb, item->key, NULL);
      int num = (int) (long) item->data; /* parmec particle number */
      parmec_get_rotation_and_position (num, bod->conf, bod->conf+9);
      SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point);
    }
#endif

    SOLFEC_Run (sol, hs->solfec_solver_kind, hs->solfec_solver, -sol->dom->step); /* negative duration used */

    parmec_steps (hs, sol->dom, actual_parmec_step, num_parmec_steps);
  }

#if TWO_SIDED_FORCE == 0
  sol->dom->skipbodies = NULL; /* undo body skipping */
#endif

  /* include contact detection between mapped Solfec bodies;
   * this may be useful in case another solver is used following this one */
  for (MAP *item = MAP_First (hs->solfec2parmec); item; item = MAP_Next (item))
  {
    for (MAP *jtem = MAP_First (hs->solfec2parmec); jtem; jtem = MAP_Next (jtem))
    {
      if (item->key != jtem->key)
      {
	AABB_Include_Body_Pair (sol->aabb, (int)(long)item->key, (int)(long)jtem->key);
      }
    }
  }

  /* TODO --> MPI */
}

/* destroy solver */
void HYBRID_SOLVER_Destroy (HYBRID_SOLVER *hs)
{
  if (hs->parmec_interval) free (hs->parmec_interval);
  MAP_Free (NULL, &hs->parmec2solfec);
  MAP_Free (NULL, &hs->solfec2parmec);
}
