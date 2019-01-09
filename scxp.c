/* Solfec's native export functionality */

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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "solfec.h"
#include "set.h"
#include "map.h"
#include "sol.h"
#include "dom.h"
#include "dio.h"
#include "iou.h"
#include "pbf.h"
#include "err.h"
#include "scxp.h"

/* mimicks sol.c:write_state */
static void write_state (SOLFEC *sol, PBF *bf, SET *subset)
{
  PBF_Time (bf, &sol->dom->time);

  dom_write_state (sol->dom, bf, subset);

  int numt = MAP_Size (sol->timers);
  PBF_Label (bf, "TIMERS");
  PBF_Int (bf, &numt, 1);
  
  for (MAP *item = MAP_First (sol->timers); item; item = MAP_Next (item))
  {
    TIMING *t = item->data;

    PBF_Label (bf, item->key);
    PBF_Double (bf, &t->total, 1);
    PBF_String (bf, (char**) &item->key);
  }
}

/* Export results in XDMF format;
 * ntimes > 0 --> number of individual time instances;
 * ntimes < 0 --> a time interval from times[0] to times[1];
 * ntimes = 0 --> export current geometry only without attributes;
 */
void solfec_export (SOLFEC *sol, double *times, int ntimes, char *path, SET *subset)
{
  PBF *bf;

  printf ("SOLFEC_EXPORT --> *** NOTE: parameter 'subset' is ignored for now. ***\n"); /* TODO */

  makedirspath (path);/* create all intermediate directories */
  char *outpath = getfilepath(path);
  bf = PBF_Write (outpath, PBF_OFF, PBF_OFF);
  char *inpout =  malloc (strlen(outpath) + 8);
  sprintf (inpout, "%s.py", outpath);
  copyfile (INPUT_FILE(), inpout);
  free (outpath);
  free (inpout);

  printf ("SOLFEC_EXPORT --> starting ...\n");

  if (ntimes < 0)
  {
    double start, end, t0, t1;

    SOLFEC_Time_Limits (sol, &start, &end);

    t0 = MAX (start, times[0]);

    t1 = MIN (end, times[1]);

    printf ("SOLFEC_EXPORT --> seeking to time %g ...\n", t0);

    SOLFEC_Seek_To (sol, t0);

    do
    {
      printf ("SOLFEC_EXPORT --> writing state at time %g ...\n", sol->dom->time);

      write_state (sol, bf, subset);

      printf ("SOLFEC_EXPORT --> moving to next time step ...\n");

      SOLFEC_Forward (sol, 1, 0);
    }
    while (sol->dom->time < t1);
  }
  else if (ntimes > 0)
  {
    for (int i = 0; i < ntimes; i ++)
    {
      printf ("SOLFEC_EXPORT --> seeking to time %g ...\n", times[i]);

      SOLFEC_Seek_To (sol, times[i]);

      printf ("SOLFEC_EXPORT --> writing state at time %g ...\n", times[i]);

      write_state (sol, bf, subset);
    }
  }
  else
  {
    printf ("SOLFEC_EXPORT --> writing state at time 0 ...\n");

    write_state (sol, bf, subset);
  }

  PBF_Close (bf);
}
