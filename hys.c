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
#include "hys.hpp" /* :) */
#include "hys.h"
#include "sol.h"
#include "dom.h"
#include "alg.h"
#include "err.h"

/* create solver */
HYBRID_SOLVER* HYBRID_SOLVER_Create (char *parmec_file, double parmec_step, double parmec_interval[2],
                 char *parmec_prefix, MAP *parmec2solfec, void *solfec_solver, int solfec_solver_kind)
{
  HYBRID_SOLVER *hs;

  ERRMEM (hs = MEM_CALLOC (sizeof (HYBRID_SOLVER)));
  hs->parmec_file = parmec_file;
  hs->parmec_step = parmec_step;
  hs->parmec_interval[0] = parmec_interval[0];
  hs->parmec_interval[1] = parmec_interval[1];
  hs->parmec_prefix = parmec_prefix;
  hs->parmec2solfec = parmec2solfec;
  hs->solfec2parmec = NULL;
  for (MAP *item = MAP_First(parmec2solfec); item; item = MAP_Next (item))
  {
    MAP_Insert (NULL, &hs->solfec2parmec, item->data, item->key, NULL);
  }
  hs->solfec_solver = solfec_solver;
  hs->solfec_solver_kind = solfec_solver_kind;

  parmec_init();

  return hs;
}

/* run solver */
void HYBRID_SOLVER_Run (HYBRID_SOLVER *hs, SOLFEC *sol, double duration)
{
  /* TODO */

  SOLFEC_Run (sol,
	      hs->solfec_solver_kind,
	      hs->solfec_solver,
	      duration);
}

/* destroy solver */
void HYBRID_SOLVER_Destroy (HYBRID_SOLVER *hs)
{
  MAP_Free (NULL, &hs->parmec2solfec);
  MAP_Free (NULL, &hs->solfec2parmec);
}
