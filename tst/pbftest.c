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

#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include "pbf.h"

static void write (PBF *bf, int n)
{
  char str [512], *pstr = str;
  double d = (double)n;
  short s = (short) (n % 32767);
  int i = n;
  unsigned int ui = (unsigned int)n;
#if !HDF5
  float f = (float)n;
  unsigned long ul = (unsigned long)n;
  long l = (long)n;
  unsigned short us = (unsigned short) (n % 32767);
  unsigned char uc = (unsigned char) (n % 127);
  char c = (char) (n % 127);
#endif

  sprintf (str, "CURRENT NUMBER IS %d", n);

  /* output time */
  PBF_Time (bf, &d);

  /* write unlabled data */
#if HDF5
  PBF_String (bf, &pstr);
  PBF_Short (bf, &s, 1);
  PBF_Int (bf, &i, 1);
  PBF_Uint (bf, &ui, 1);
  PBF_Double (bf, &d, 1);
#else
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
#endif

  /* write labeled data */
#if HDF5
  PBF_Label (bf, "STRING");
  PBF_String (bf, &pstr);
  PBF_Label (bf, "SHORT");
  PBF_Short (bf, &s, 1);
  PBF_Label (bf, "INT");
  PBF_Int (bf, &i, 1);
  PBF_Label (bf, "UINT");
  PBF_Uint (bf, &ui, 1);
  PBF_Label (bf, "DOUBLE");
  PBF_Double (bf, &d, 1);
#else
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
#endif
}

static int read (PBF *bf, int n)
{
  char *pstr = NULL;
  double d;
  short s;
  int i;
  unsigned int ui;
#if !HDF5
  float f;
  unsigned long ul;
  long l;
  unsigned short us;
  unsigned char uc;
  char c;
#endif

  char _str [512];
  double _d = (double)n;
  short _s = (short) (n % 32767);
  unsigned int _ui = (unsigned int)n;
  int _i = n;

#if !HDF5
  float _f = (float)n;
  unsigned long _ul = (unsigned long)n;
  long _l = (long)n;
  unsigned short _us = (unsigned short) (n % 32767);
  unsigned char _uc = (unsigned char) (n % 127);
  char _c = (char) (n % 127);
#endif

  sprintf (_str, "CURRENT NUMBER IS %d", n);

  PBF_Time (bf, &d);
  if (d != (double)n) return 0;

#if HDF5
  PBF_String (bf, &pstr);
  PBF_Short (bf, &s, 1);
  PBF_Int (bf, &i, 1);
  PBF_Uint (bf, &ui, 1);
  PBF_Double (bf, &d, 1);

  if (strcmp (pstr, _str) != 0 || d != _d
      || ui != _ui || i != _i || s != _s) return 0;
#else
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

  if (strcmp (pstr, _str) != 0 ||
      d != _d || f != _f || ul != _ul ||
      l != _l || ui != _ui || i != _i ||
      us != _us || s != _s || uc != _uc || c != _c) return 0;
#endif

  free (pstr); pstr = NULL;
  d = 0.; s = 0;  ui = 0; i = 0;
#if !HDF5
  f = 0.; ul = 0; l = 0; us = 0;
  uc = 0; c = 0;
#endif

#if HDF5
  PBF_Label (bf, "STRING");
  PBF_String (bf, &pstr);
  PBF_Label (bf, "SHORT");
  PBF_Short (bf, &s, 1);
  PBF_Label (bf, "INT");
  PBF_Int (bf, &i, 1);
  PBF_Label (bf, "DOUBLE");
  PBF_Double (bf, &d, 1);

  if (strcmp (pstr, _str) != 0 || d != _d
      || ui != _ui || i != _i || s != _s) return 0;
#else
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

  if (strcmp (pstr, _str) != 0 ||
      d != _d || f != _f || ul != _ul ||
      l != _l || ui != _ui || i != _i ||
      us != _us || s != _s || uc != _uc || c != _c) return 0;
#endif

  free (pstr);

  return 1;
}

int main (int argc, char **argv)
{
  char *s, str [512];
  int n, nmax;
  double d;
  PBF *bf;

  if (argc == 1)
  {
    printf ("pbftest [NUMBER OF SAMPLES] [COMPRESSION FLAG]\n");
    return 0;
  }

  if (argc >= 2 && isdigit (argv [1][0]))
    nmax = atoi (argv [1]);
  else nmax = 10;

  bf = PBF_Write ("pbftest.state", PBF_OFF, PBF_OFF);

  if (argc >= 3)
  {
    if (atoi (argv [2]) == 1)
    {
      bf->compression = PBF_ON; /* enable compression */
    }
  }

  for (n = 1; n <= nmax; n ++) write (bf, n);
  PBF_Close (bf);

  bf = PBF_Read ("pbftest.state");
  for (n = 1; n <= nmax; n ++) 
  {
    if (!read (bf, n))
    {
      fprintf (stderr, "FAILED (reading all)\n");
      return 1;
    }
    PBF_Forward (bf, 1);
  }

  /* read every third double */
  PBF_Seek (bf, 0);
  for (n = 1; n <= nmax; n += 3)
  {
    PBF_Label (bf, "DOUBLE");
    PBF_Double (bf, &d, 1);
    if (d != (double) n)
    {
      fprintf (stderr, "FAILED (forward reading every 3rd labeled double)\n");
      return 1;
    }
    PBF_Forward (bf, 3);
  }

  /* read every fifth string */
  PBF_Seek (bf, 0);
  for (n = 1; n <= nmax; n += 5)
  {
    PBF_Label (bf, "STRING");
    s = NULL;
    PBF_String (bf, &s);
    sprintf (str, "CURRENT NUMBER IS %d", n);
    if (strcmp (s, str) != 0)
    {
      fprintf (stderr, "FAILED (forward reading every 5th labeled string)\n");
      return 1;
    }
    PBF_Forward (bf, 5);
    free (s);
  }

  /* do the same backwards */
  PBF_Seek (bf, nmax);
  for (n = nmax; n >= 1; n -= 3)
  {
    PBF_Label (bf, "DOUBLE");
    PBF_Double (bf, &d, 1);
    if (d != (double) n)
    {
      fprintf (stderr, "FAILED (backward reading every 3rd labeled double)\n");
      return 1;
    }
    PBF_Backward (bf, 3);
  }
  PBF_Seek (bf, nmax);
  for (n = nmax; n >= 1; n -= 5)
  {
    PBF_Label (bf, "STRING");
    s = NULL;
    PBF_String (bf, &s);
    sprintf (str, "CURRENT NUMBER IS %d", n);
    if (strcmp (s, str) != 0)
    {
      fprintf (stderr, "FAILED (forward reading every 5th labeled string)\n");
      return 1;
    }
    PBF_Backward (bf, 5);
    free (s);
  }
 
  PBF_Close (bf);

  printf ("PASSED\n");
  return 0;
}
