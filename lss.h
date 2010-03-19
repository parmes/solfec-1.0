/*
 * lss.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * -----------------------------------------------------------
 * linear system solver based on [1, 2, 3]
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

#ifndef __lss__
#define __lss__

/* ============ ERROR ============ | ============================ DESCRIPTION ============================ */
                                 /*|                                                                       */
enum lsserr                      /*|                                                                       */
{                                /*|                                                                       */
  LSSERR_NONE = 0,               /*| there was no error                                                    */ 
  LSSERR_OUT_OF_MEMORY,          /*| system memory has run out                                             */
  LSSERR_INVALID_ARGUMENT,       /*| invalid argument was used in function call                            */
  LSSERR_LACK_OF_CONVERGENCE,    /*| number of iterations has exceeded the prescribed bound                */
  LSSERR_EMPTY_COLUMN,           /*| system matrix has an empty column                                     */
  LSSERR_ZERO_ON_DIAGONAL        /*| system matrix has zero on the diagonal                                */
};                               /*|_______________________________________________________________________*/

typedef enum lsserr LSSERR;

/* ========== PARAMETER ========== | ====== VALUE ====== | ============================ DESCRIPTION ============================ */
                                 /*|  DEFAULT | ACCEPTED |                                                                       */
enum lsspar                      /*|----------|----------|-----------------------------------------------------------------------*/
{                                /*|          |          | Read-write parameters, affecting a next call to LSS_Run:              */
                                 /*|----------|----------|-----------------------------------------------------------------------*/
  LSS_ITERATIONS_BOUND,          /*| 1000     | any > 0  | k < VALUE for a sequence of approximate solutions x(k)                */
  LSS_RELATIVE_ACCURACY,         /*| 1E-6     | any > 0  | |A x(k) - b| < VALUE * |A x(0) - b| when LSS_Solve returns 0          */
  LSS_ABSOLUTE_ACCURACY,         /*| 1E-3     | any > 0  | |A x(k) - b| < VALUE when LSS_Solve returns 0                         */
  LSS_PRECONDITIONER,            /*| 1        | 0 ... 3  | switch between wavelet based preconditioners by Pereira et al. [5]    */
  LSS_DECIMATION,                /*| 4        | any >= 2 | preconditioner wavelet decimation (sub-sampling) value                */
  LSS_CUTOFF,                    /*| 16       | any >= 2 | coarsest level matrix dimension cutoff bound                          */
  LSS_RESTART,                   /*| 32       | any >= 1 | restart bound for top level matrix GMRES run                          */ 
  LSS_COARSE_RESTART,            /*| 16       | any >= 1 | restart bound for coarsest level matrix GMRES run                     */ 
  LSS_COARSE_ITERATIONS_BOUND,   /*| 16       | any >= 1 | iterations bound for coarsest level matrix GMRES run                  */ 
  LSS_COARSE_RELATIVE_ACCURACY,  /*| 1E-6     | any > 0  | relative accuracy for coarsest level GMRES run                        */
  LSS_COARSE_ABSOLUTE_ACCURACY,  /*| 1E-3     | any > 0  | absolute accuracy for coarsest level GMRES run                        */
  LSS_SMOOTHING_STEPS,           /*| 1        | any >= 0 | number of Gauss-Seidel sweeps on each level of preconditoning         */
                                 /*|----------|----------|-----------------------------------------------------------------------*/
                                 /*|          |          | Read-only parameters, corresponding to a recent call to LSS_Solve:    */
                                 /*|----------|----------|-----------------------------------------------------------------------*/
  LSS_ITERATIONS,                /*| none     | none     | number of performed iterations                                        */
  LSS_RELATIVE_ERROR,            /*| none     | none     | history (LSS_Ggetv) or last value (LSS_Get) of relative error         */
  LSS_ABSOLUTE_ERROR,            /*| none     | none     | history (LSS_Ggetv) or last value (LSS_Get) of absolute error         */ 
  LSS_LEVELS,                    /*| none     | none     | number of preconditioner levels                                       */
  LSS_LEVEL_DIMENSIONS,          /*| none     | none     | vector (LSS_Getv) of dimensions of reduced system matrices            */
  LSS_LEVEL_NONZEROS,            /*| none     | none     | vector (LSS_Getv) of numbers of nonzero entries in reduced matrices   */
  LSS_OPERATOR_COMPLEXITY,       /*| none     | none     | SUM (nnz(A_i)) / nnz (A_1), nnz(A) = number of nonzero elements in A  */
  LSS_GRID_COMPLEXITY            /*| none     | none     | SUM (dim(A_i)) / dim (A_1), dim(A) = dimension of matrix A            */
};                               /*|__________|__________|_______________________________________________________________________*/

typedef enum lsspar LSSPAR;

/* === PRECONDITIONER === | ========== DESCRIPTION =========== */
/*           0            | none                               */ 
/*           1            | Daubechies-2 (Haar)                */
/*           2            | Daubechies-4                       */
/*           3            | Daubechies-6                       */
/*________________________|____________________________________*/

/* n: dimension of n x n system A matrix
 * n < 0 => A in compressed row format (CSR)
 * n > 0 => A in compressed column formt (CSC)
 * a: system A matrix values
 * p: pointers to rows (CSR) or columns (CSC) of size (n+1)
 * i: indices of columns (CSR) or rows (CSC) corresponding to 'a';
 * create linear system solver */
void* LSS_Create (int n, double *a, int *p, int *i);

/* lss: returned by LSS_Create
 * p: parameter
 * v: value;
 * set parameter */
LSSERR LSS_Set (void *lss, LSSPAR p, double v);

/* lss: returned by LSS_Create
 * p: parameter
 * return: value;
 * get parameter */
double LSS_Get (void *lss, LSSPAR p);

/* lss: returned by LSS_Create 
 * p: parameter
 * return: vector;
 * get vector parameter */
double* LSS_Getv (void *lss, LSSPAR p);

/* lss: returned by LSS_Create  
 * a: system A matrix values (same (p, i) as in LSS_Create is assumed)
 * x: solution vector
 * b: right hand side vector
 * solve linear system A x = b, where A = (a, p, i) */
LSSERR LSS_Solve (void *lss, double *a, double *x, double *b);

/* lss: returned by LSS_Create;
 * return last error message (or NULL) */
char* LSS_Errmsg (void *lss);

/* lss: returned by LSS_Create;
 * release solver memory */
void LSS_Destroy (void *lss);

/* ========== USAGE ==========
 *
 * 1. Call LSS_Create.
 * 
 * 2. Call LSS_Set if needed.
 *
 * 3. Call LSS_Solve.
 *
 * 4. Call LSS_Get, LSS_Getv if needed.
 *
 * 5. Go to 3 if needed.
 *
 * 6. Call LSS_Destroy.
 *
 * ========== REFERENCES ==========
 *
 * [1] Y. Saad and M. Schultz. GMRES: a generalized minimal residual algorithm for solving
 *     nonsymmetric linear systems. SIAM J. Sci. Stat. Comput. Vol. 7, No. 3, pp. 856-869, 1986.
 *
 * [2] Y. Saad. A flexible inner-outer preconditioned GMRES algorithm.
 *     SIAM J. Sci. Comput. Vol. 14, No. 2, pp. 461-469, 1993.
 *
 * [3] F. Pereira, S. Verardi, S. Nabeta. A wavelet-based algebraic multigrid preconditioner
 *     for sparse linear systems. Applied Mathematics and Computation, Vol. 182, pp. 1098-1107, 2006.
 */

#endif
