/*
 * mtx.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * general sparse or dense matrix
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

#ifndef __mx__
#define __mx__

typedef struct general_matrix MX;

struct general_matrix
{
  enum {MXDENSE  = 0x01,        /* dense */
        MXBD     = 0x02,        /* block diagonal (square) */
	MXCSC    = 0x04} kind;  /* compressed columns */
  
  enum {MXTRANS  = 0x01,        /* transposed matrix */
	MXSTATIC = 0x02,        /* static matrix */
        MXDSUBLK = 0x04,        /* diagonal sub-block */
        MXIFAC   = 0x08} flags; /* factorised sparse inverse */

  int nzmax,   /* number of nonzero entries */
          m,   /* number of rows (DENSE, BD (and columns), CSC) */
	  n,   /* number of columns (DENSE, CSC), or number of blocks (BD) */
	 *p,   /* pointers to columns (CSC), or blocks (BD); p[n] == nzmax */
	 *i,   /* indices of column row entries (CSC, i[0...nzmax-1]), or indices
		  of the first block row/column (BD, i[n] = m) */
	 nz;   /* number of entries in triplet matrix, -1 for compressed-col */

  double *x;   /* values, x[0...nzmax-1] */

  void *sym,   /* symbolic factorisation for CSC inverse */
       *num;   /* numeric factorisation for CSC inverse */
};

/* static dense matrix */
#define MX_DENSE(name, m, n)\
  double __##name [m*n];\
  MX name = {MXDENSE, MXSTATIC, m*n, m, n, NULL, NULL, 0, __##name, NULL, NULL}

/* static dense matrix */
#define MX_DENSE_PTR(name, m, n, ptr)\
  MX name = {MXDENSE, MXSTATIC, m*n, m, n, NULL, NULL, 0, ptr, NULL, NULL}

/* static block diagonal matrix */
#define MX_BD(name, nzmax, m, n, p, i)\
  double __##name [nzmax];\
  MX name = {MXBD, MXSTATIC, nzmax, m, n, p, i, 0, __##name, NULL, NULL}

/* static sparse matrix */
#define MX_CSC(name, nzmax, m, n, p, i)\
  double __##name [nzmax];\
  MX name = {MXCSC, MXSTATIC, nzmax, m, n, p, i, 0, __##name, NULL, NULL}

/* create a matrix => structure tables (p, i) always have
 * to be provided; the tables 'p' and 'i' are coppied */
MX* MX_Create (short kind, int m, int n, int *p, int *i);

/* set to zero */
void MX_Zero (MX *a);

/* scale */
void MX_Scale (MX *a, double b);

/* set b = a; if 'b' == NULL return
 * new matrix; otherwise return 'b' */
MX* MX_Copy (MX *a, MX *b);

/* returned = transpose (a) => IMPORTANT: never use the transposition to set a pointer;
 * it can only be used "on the fly" in order to modify an input to other routines */
MX* MX_Tran (MX *a);

/* returned = diagonal sub-block (a) => IMPORTANT: never use this function to set a pointer;
 * it can only be used "on the fly" in order to modify an input to other routines; Applies
 * only to block diagonal (BD) matrices; (from, to) correspond to blocks range (inclusive) */
MX* MX_Diag (MX *a, int from, int to);

/* sum of two matrices => c = alpha * a + beta * b;
 * if 'c' == NULL return new matrix; otherwise return 'c' */
MX* MX_Add (double alpha, MX *a, double beta, MX *b, MX *c);

/* matrix matrix product => c = alpha * a * b + beta * c;
 * if 'c' == NULL return new matrix; otherwise return 'c' */
MX* MX_Matmat (double alpha, MX *a, MX *b, double beta, MX *c);

/* matrix vector product => c = alpha * a *b + beta * c */
void MX_Matvec (double alpha, MX *a, double *b, double beta, double *c);

/* triple matrix product => d = a * b * c;
 * if 'd' == NULL return new matrix; otherwise return 'd' */
MX* MX_Trimat (MX *a, MX *b, MX *c, MX *d);

/* inverse => b = inv (a); factorization is used for CSC;
 * if 'b' == NULL return new matrix; otherwise return 'b' */
MX* MX_Inverse (MX *a, MX *b);

/* compute |n| eigenvalues & eigenvectors (vec != NULL) in the upper or
 * lower range (n < 0 or n > 0) => symmetry of 'a' is assumed and the
 * results are outputed according to the ascending order of eigenvalues */
void MX_Eigen (MX *a, int n, double *val, MX *vec);

/* free matrix */
void MX_Destroy (MX *a);

#endif
