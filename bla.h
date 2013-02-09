/*
 * bla.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * BLAS interface
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

#ifndef __blas__
#define __blas__

double ddot_ (int *n, double *dx, int *incx, double *dy, int *incy);

void dscal_ (int *n, double *da, double *dx, int *incx);

void dcopy_ (int *n, double *dx, int *incx, double *dy, int *incy);

void daxpy_ (int *n, double *da, double *dx, int *incx, double *dy, int *incy);

void dgemm_ (char *transa, char *transb, int *m, int *n, int *k,
  double *alpha, double *a, int *lda, double *b, int *ldb,
  double *beta, double *c, int *ldc);

void dgemv_ (char *trans, int *m, int *n, double *alpha,
  double *a, int *lda, double *x, int *incx, double *beta,
  double *y, int *incy);

inline static double blas_ddot (int n, double *dx, int incx, double *dy, int incy)
{ return ddot_ (&n, dx, &incx, dy, &incy); }

inline static void blas_dscal (int n, double da, double *dx, int incx)
{ dscal_ (&n, &da, dx, &incx); }

inline static void blas_dcopy (int n, double *dx, int incx, double *dy, int incy)
{ dcopy_ (&n, dx, &incx, dy, &incy); }

inline static void blas_daxpy (int n, double da, double *dx, int incx, double *dy, int incy)
{ daxpy_ (&n, &da, dx, &incx, dy, &incy); }

inline static void blas_dgemm (char transa, char transb, int m, int n, int k,
  double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc)
{
  dgemm_ (&transa, &transb, &m, &n, &k,
    &alpha, a, &lda, b, &ldb,
    &beta, c, &ldc);
}

inline static void blas_dgemv (char trans, int m, int n, double alpha,
  double *a, int lda, double *x, int incx, double beta, double *y, int incy)
{
  dgemv_ (&trans, &m, &n, &alpha,
    a, &lda, x, &incx, &beta, y, &incy);
}

#endif
