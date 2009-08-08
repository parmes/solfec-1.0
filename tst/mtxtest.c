/*
 * mtxtest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of matrix routines
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

#include <stdio.h>
#include <stdlib.h>
#include "mtx.h"
#include "alg.h"

int main (int argc, char **argv)
{
  MX_DENSE (a, 3, 3);

  IDENTITY (a.x);
  SCALE9 (a.x, 2.0);

  MX_Inverse (&a, &a);

  printf ("[%.2f, %.2f, %.2f]\n[%.2f, %.2f, %.2f]\n[%.2f, %.2f, %.2f]\n",
	   a.x [0], a.x [3], a.x [6], a.x [1], a.x [4], a.x [7], a.x [2], a.x [5], a.x [8]);

  return 0;
}
