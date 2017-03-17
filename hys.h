/* HYBRID_SOLVER interface */

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

#include "sol.h"
#include "map.h"

#ifndef __hys__
#define __hys__

typedef struct hybrid_solver HYBRID_SOLVER;

struct hybrid_solver
{
  char *parmec_file;
  double parmec_step;
  double *parmec_interval;
  void** parmec_interval_func;
  int* parmec_interval_tms;
  char *parmec_prefix;
  MAP *parmec2solfec;
  MAP *solfec2parmec;
  void *solfec_solver;
  int solfec_solver_kind;
};

/* create solver */
HYBRID_SOLVER* HYBRID_SOLVER_Create (char *parmec_file, double parmec_step,
           MAP *parmec2solfec, void *solfec_solver, int solfec_solver_kind);

/* run solver */
void HYBRID_SOLVER_Run (HYBRID_SOLVER *hs, SOLFEC *sol, double duration);

/* destroy solver */
void HYBRID_SOLVER_Destroy (HYBRID_SOLVER *hs);

#endif
