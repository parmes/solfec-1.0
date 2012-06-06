/*
 * tsi.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * triangle-sphere intersection
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

#ifndef __tsi__
#define __tsi__

enum
{
  TSI_A = 0,
  TSI_B,
  TSI_C,
  TSI_A_B,
  TSI_B_C,
  TSI_C_A,
  TSI_A_B_C,
  TSI_A_BC,
  TSI_B_CA,
  TSI_C_AB,
  TSI_AB,
  TSI_BC,
  TSI_CA,
  TSI_AB_BC,
  TSI_BC_CA,
  TSI_CA_AB,
  TSI_AB_BC_CA,
  TSI_BUBBLE,
  TSI_OUT
};

/* get status of intersection between triangle (a, b, c) and sphere (p, r) */
int TSI_Status (double *a, double *b, double *c, double *p, double r);

/* approximate intersection of triangle (a, b, c) and sphere (p, r) with second order triangular elements;
 * output vector 'cc' contains 'ncc' 2-dimensional points in the natural coordinates of the input triangle;
 * output triangulation 'tt' contains 'ntt' 6-vertex second order triangle vertex numbers, refering to 'cc';
 * output arrays are dynamically allocated; the status of the intersection is returned */
int TSI_Approx (double *a, double *b, double *c, double *p, double r, double **cc, int *ncc, int **tt, int *ntt);

#endif
