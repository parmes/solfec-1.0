/*
 * svk.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * Saint Venant - Kirchhoff material
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

#ifndef __svk__
#define __svk__

/* given Lame coefficients (lambda, mi) and deformation gradient 'F'
 * return internal hyperelastic energy (scaled by the 'volume') */
double SVK_Energy_R (double lambda, double mi, double volume, double *F); /* F is row-wise */
double SVK_Energy_C (double lambda, double mi, double volume, double *F); /* F is column-wise */

/* given Lame coefficients (lambda, mi) and deformation gradient 'F'
 * return det (F) end output first Piola-Kirchhoff stress 'P' (scaled by the 'volume') */
double SVK_Stress_R (double lambda, double mi, double volume, double *F, double *P); /* F, P are row-wise */
double SVK_Stress_C (double lambda, double mi, double volume, double *F, double *P); /* F, P are column-wise */
  
/* given Lame coefficients (lambda, mi) and deformation gradient 'F' output 9 x 9
 * tangent operator (scaled by the 'volume') into the matrix 'K' with leading column 'dim'ension (K is column-wise) */
void SVK_Tangent_R (double lambda, double mi, double volume, int dim, double *F, double *K); /* F is row-wise */
void SVK_Tangent_C (double lambda, double mi, double volume, int dim, double *F, double *K); /* F is column-wise */

/* as above, but returns only diagonal of the tangent matrix */
void SVK_Tangent_Diagonal_R (double lambda, double mi, double volume, int dim, double *F, double *K); /* F is row-wise */
void SVK_Tangent_Diagonal_C (double lambda, double mi, double volume, int dim, double *F, double *K); /* F is column-wise */

/* given Lame coefficients (lambda, mi) and deformation gradient 'F' output 9 x 9
 * tangent operator's partial derivative with respect to 'F[comp]' (scaled by the 'volume')
 * into the the matrix 'K' with leading column 'dim'ension (K is column-wise) */
void SVK_Tangent_Derivative_R (int comp, double lambda, double mi, double volume, int dim, double *F, double *K); /* F is row-wise */
void SVK_Tangent_Derivative_C (int comp, double lambda, double mi, double volume, int dim, double *F, double *K); /* F is column-wise */

#endif
