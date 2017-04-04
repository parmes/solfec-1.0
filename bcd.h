/*
The MIT License (MIT)

Copyright (c) 2017 EDF Energy

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

/* Purpose: co-rotated FEM displacements sampling */

#include "set.h"

#ifndef __bcd__
#define __bcd__

typedef struct corotated_displacements BCD;

struct corotated_displacements /* body co-rotated FEM displacements sampling list */
{
  SET *subset; /* set of body identifiers (MPI-ranks-global in parallel) */

  double *sampling; /* sampling of time instants */

  int length; /* if  0, sampling[0] stores output interval;
                 if -1, use SOLFEC-->output_interval;
		 otherise, sampling[0...length-1] stores a list of time instantcs at which to sample */

  int size; /* sample size == (BODY*) subset --> data --> dofs */

  SET *pending; /* pending samples */

  void *output; /* Python list of lists of output samples (!= NULL on MPI rank 0 in parallel) */

  double latest; /* latest sample time */

  BCD *next; /* next sampling definition/data in list */
};

/* append 'SOLFEC->bcd' with new definition and return the appended list head; return 1 upon succes; 0 otherwise */
int BCD_Append (SOLFEC *sol, SET *subset, double *sampling, int length, void *output);

/* sample 'bcd' list at current time --> goes into pending samples */
void BCD_Sample (SOLFEC *sol, BCD *bcd);

/* output pending samples */
void BCD_Append_Output (BCD *bcd);

/* destroy sampling list */
void BCD_Destroy (BCD *bcd);

#endif
