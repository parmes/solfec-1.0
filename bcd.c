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

/* Contributors: Tomasz Koziara */

#include "sol.h"
#include "bcd.h"

/* append sampling list 'bcd' with new definition and return the appended list head */
BCD* BCD_Append (BCD *bcd, void *definition, void *subset, double *sampling, int length, void *output)
{
  /* TODO */
  return NULL;
}

/* sample 'bcd' list at current time */
void BCD_Sample (SOLFEC *solfec, BCD *bcd)
{
  /* TODO */
}
#if MPI
/* gather output lists at rank 0 process */
void BCD_Rank0_Gather (BCD *bcd)
{
}
#endif

/* destroy sampling list */
void BCD_Destroy (BCD *bcd)
{
  /* TODO */
}
