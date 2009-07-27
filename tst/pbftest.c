/*
 * pbftest.c
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of the portable binary file interface
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

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "pbf.h"

static void write (PBF *bf, int n)
{
  char str [512], *pstr = str;
  double d = (double)n;
  float f = (float)n;
  unsigned long ul = (unsigned long)n;
  long l = (long)n;
  unsigned int ui = (unsigned int)n;
  int i = n;
  unsigned short us = (unsigned short)n;
  short s = (short)n;
  unsigned char uc = (unsigned char)n;
  char c = (char)n;

  sprintf (str, "CURRENT NUMBER IS %d", n);

  /* output time */
  PBF_Time (bf, &d);

  /* write unlabled data */
  PBF_String (bf, &pstr);
  PBF_Double (bf, &d, 1);
  PBF_Float (bf, &f, 1);
  PBF_Ulong (bf, &ul, 1);
  PBF_Long (bf, &l, 1);
  PBF_Uint (bf, &ui, 1);
  PBF_Int (bf, &i, 1);
  PBF_Ushort (bf, &us, 1);
  PBF_Short (bf, &s, 1);
  PBF_Uchar (bf, &uc, 1);
  PBF_Char (bf, &c, 1);

  /* write labeled data */
  PBF_Label (bf, "STRING");
  PBF_String (bf, &pstr);
  PBF_Label (bf, "DOUBLE");
  PBF_Double (bf, &d, 1);
  PBF_Label (bf, "FLOAT");
  PBF_Float (bf, &f, 1);
  PBF_Label (bf, "ULONG");
  PBF_Ulong (bf, &ul, 1);
  PBF_Label (bf, "LONG");
  PBF_Long (bf, &l, 1);
  PBF_Label (bf, "UINT");
  PBF_Uint (bf, &ui, 1);
  PBF_Label (bf, "INT");
  PBF_Int (bf, &i, 1);
  PBF_Label (bf, "USHORT");
  PBF_Ushort (bf, &us, 1);
  PBF_Label (bf, "SHORT");
  PBF_Short (bf, &s, 1);
  PBF_Label (bf, "UCHAR");
  PBF_Uchar (bf, &uc, 1);
  PBF_Label (bf, "CHAR");
  PBF_Char (bf, &c, 1);
}

static void read (PBF *bf)
{
  char *pstr = NULL;
  double d;
  float f;
  unsigned long ul;
  long l;
  unsigned int ui;
  int i;
  unsigned short us;
  short s;
  unsigned char uc;
  char c;

  printf ("-------------------------------------------\n");
  PBF_Time (bf, &d);
  printf ("TIME: %f\n", d);

  PBF_String (bf, &pstr);
  PBF_Double (bf, &d, 1);
  PBF_Float (bf, &f, 1);
  PBF_Ulong (bf, &ul, 1);
  PBF_Long (bf, &l, 1);
  PBF_Uint (bf, &ui, 1);
  PBF_Int (bf, &i, 1);
  PBF_Ushort (bf, &us, 1);
  PBF_Short (bf, &s, 1);
  PBF_Uchar (bf, &uc, 1);
  PBF_Char (bf, &c, 1);

  printf ("*** UNLABELED:\n");
  printf ("STRING: %s\n", pstr);
  printf ("DOUBLE: %lf\n", d);
  printf ("FLOAT: %f\n", f);
  printf ("ULONG: %lu\n", ul);
  printf ("LONG: %ld\n", l);
  printf ("UINT: %u\n", ui);
  printf ("INT: %d\n", i);
  printf ("USHORT: %u\n", us);
  printf ("SHORT: %d\n", s);
  printf ("UCHAR: %u\n", uc);
  printf ("CHAR: %d\n", c);

  free (pstr); pstr = NULL;
  d = 0.; f = 0.; ul = 0; l = 0;
  ui = 0; i = 0; us = 0; s = 0;
  uc = 0; c = 0;

  PBF_Label (bf, "STRING");
  PBF_String (bf, &pstr);
  PBF_Label (bf, "DOUBLE");
  PBF_Double (bf, &d, 1);
  PBF_Label (bf, "FLOAT");
  PBF_Float (bf, &f, 1);
  PBF_Label (bf, "ULONG");
  PBF_Ulong (bf, &ul, 1);
  PBF_Label (bf, "LONG");
  PBF_Long (bf, &l, 1);
  PBF_Label (bf, "UINT");
  PBF_Uint (bf, &ui, 1);
  PBF_Label (bf, "INT");
  PBF_Int (bf, &i, 1);
  PBF_Label (bf, "USHORT");
  PBF_Ushort (bf, &us, 1);
  PBF_Label (bf, "SHORT");
  PBF_Short (bf, &s, 1);
  PBF_Label (bf, "UCHAR");
  PBF_Uchar (bf, &uc, 1);
  PBF_Label (bf, "CHAR");
  PBF_Char (bf, &c, 1);

  printf ("*** LABELED:\n");
  printf ("STRING: %s\n", pstr);
  printf ("DOUBLE: %lf\n", d);
  printf ("FLOAT: %f\n", f);
  printf ("ULONG: %lu\n", ul);
  printf ("LONG: %ld\n", l);
  printf ("UINT: %u\n", ui);
  printf ("INT: %d\n", i);
  printf ("USHORT: %u\n", us);
  printf ("SHORT: %d\n", s);
  printf ("UCHAR: %u\n", uc);
  printf ("CHAR: %d\n", c);
  free (pstr);
}

int main (int argc, char **argv)
{
  PBF *bf;
  int n, nmax;
  double d;
  char *s;

  if (argc >= 2 && isdigit (argv [1][0]))
    nmax = atoi (argv [1]);
  else nmax = 10;

  bf = PBF_Write ("pbftest.state");
  for (n = 1; n <= nmax; n ++) write (bf, n);
  PBF_Close (bf);

  bf = PBF_Read ("pbftest.state");
  for (n = 1; n <= nmax; n ++) 
  {
    read (bf);
    PBF_Forward (bf, 1);
  }

  /* go every third and select one
   * labeled value */
  PBF_Seek (bf, 0);
  for (n = 1; n <= nmax; n += 3)
  {
    PBF_Label (bf, "DOUBLE");
    PBF_Double (bf, &d, 1);
    printf ("Selected every third DOUBLE = %f\n", d);
    PBF_Forward (bf, 3);
  }

  /* go every fifth and select one
   * labeled value */
  PBF_Seek (bf, 0);
  for (n = 1; n <= nmax; n += 5)
  {
    PBF_Label (bf, "STRING");
    s = NULL;
    PBF_String (bf, &s);
    printf ("Selected every fifth STRING = %s\n", s);
    PBF_Forward (bf, 5);
    free (s);
  }

  /* do the same backwards */
  PBF_Seek (bf, nmax);
  for (n = 1; n <= nmax; n += 3)
  {
    PBF_Label (bf, "DOUBLE");
    PBF_Double (bf, &d, 1);
    printf ("Selected every third DOUBLE = %f\n", d);
    PBF_Backward (bf, 3);
  }
  PBF_Seek (bf, nmax);
  for (n = 1; n <= nmax; n += 5)
  {
    PBF_Label (bf, "STRING");
    s = NULL;
    PBF_String (bf, &s);
    printf ("Selected every fifth STRING = %s\n", s);
    PBF_Backward (bf, 5);
    free (s);
  }
 
  PBF_Close (bf);
  return 0;
}
