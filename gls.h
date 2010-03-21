/*
 * gls.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Gajulapalli-Lasdon matrix scaling
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

#ifndef __gls__
#define __gls__

/* n: matrix dimension (input)
 * q: matrix coefficient (intput/output); scaled on exit
 * p: pointers to columns (input)
 * i: row indices (input)
 * x: row scaling (output)
 * y: column scaling (output)
 * return: 1 when done or 0 when out of memory;
 * try finding x and y such that x[i] |a[i, j]| y [j] is in [0.5, 1] */
int gls_csc (int n, double *q, int *p, int *i, double *x, double *y);

#endif
