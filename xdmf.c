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

#if POSIX
#include <sys/stat.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sol.h"
#include "err.h"

/* returns xmf file path string; creatses intermediate directories */
static char *xmf_path_and_dirs (char *path)
{
  char *out;

  int l = strlen (path);

#if POSIX
  for (int i = 0; i < l; i ++) /* create all intermediate directories */
  {
    if (path [i] == '/')
    {
       path [i] = '\0';
       mkdir (path, 0777); /* POSIX */
       path [i] = '/';
    }
  }
  mkdir (path, 0777); /* POSIX */
#endif

  while (l > 0 && path [l-1] != '/') l --;
  char *lastname = &path [l];
  l = strlen (path);
  int n = l + strlen (lastname) + 8;
  ERRMEM (out = malloc (n));
  strcpy (out, path);
  out [l] = '/';
  strcpy (out+l+1, lastname);
  l = strlen (out);
  sprintf (&out[l], ".xmf");

  return out;
}

/* Export results in XMDF format;
 * ntimes > 0 --> number of individual time instances;
 * ntimes < 0 --> a time interval from times[0] to times[1];
 */
void xdmf_export (SOLFEC *sol, double *times, int ntimes, char *path)
{

  char *xmf_path = xmf_path_and_dirs (path);
  FILE *xmf_file = fopen (xmf_path, "w");
  ASSERT_TEXT (xmf_file, "Opening XDMF markup file %s has failed", xmf_path);

  fprintf (xmf_file, "<Xdmf>\n");
  fprintf (xmf_file, "<Domain>\n");
  fprintf (xmf_file, "<Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");

#if 0
  if (ntimes < 0)
  {
    double start, end, t0, t1;

    SOLFEC_Time_Limits (sol, &start, &end);

    t0 = MAX (start, times[0]);

    t1 = MIN (end, times[1]);

    SOLFEC_Seek_To (sol, t0);

    do
    {
      xdmf_step (sol->dom, xmf_file, path);

      SOLFEC_Forward (sol, 1);
    }
    while (sol->dom->time < t1);
  }
  else
  {
    for (int i = 0; i < ntimes; i ++)
    {
      SOLFEC_Seek_To (sol, times[i]);

      xdmf_step (sol->dom, xmf_file, path);
    }
  }
#endif

  fprintf (xmf_file, "</Grid>\n");
  fprintf (xmf_file, "</Domain>\n");
  fprintf (xmf_file, "</Xdmf>\n");

  fclose (xmf_file);
  free (xmf_path);
}
