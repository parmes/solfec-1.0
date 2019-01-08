/* file, directory and path utilites */

/*
The MIT License (MIT)

Copyright (c) 2008 Tomasz Koziara

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

#if POSIX
#include <sys/stat.h>
#endif

#include "err.h"

#ifndef __iou__
#define __iou__

/* copy a file */
static void copyfile (char *from, char *to)
{
  FILE *f, *t;
  char c;

  ASSERT (f = fopen (from, "rb"), ERR_FILE_OPEN);
  ASSERT (t = fopen (to, "wb"), ERR_FILE_OPEN);

  while (!feof(f))
  {
    c = fgetc(f);
    ASSERT (!ferror(f), ERR_FILE_READ);
    if(!feof(f)) fputc(c, t);
    ASSERT (!ferror(t), ERR_FILE_WRITE);
  }

  ASSERT (fclose (f) == 0, ERR_FILE_CLOSE);
  ASSERT (fclose (t) == 0, ERR_FILE_CLOSE);
}

static void makedirspath (char *path)
{
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
}

/* from directory path get last name */
static char *lastname (char *path)
{
  int l = strlen (path);

  while (l > 0 && path [l-1] != '/') l --;

  return &path [l];
}

/* get file path from directory path */
static char *getfilepath (char *outpath)
{
  int l = strlen (outpath),
      n = l + strlen (lastname (outpath)) + 8;
  char *path;

  ERRMEM (path = malloc (n));
  strcpy (path, outpath);
  path [l] = '/';
  strcpy (path+l+1, lastname (outpath));

  return path;
}

#endif
