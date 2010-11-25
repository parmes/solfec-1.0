/*
 * lap.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * LAPACK interface
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

#ifndef __lapack__
#define __lapack__

int dgetrf_ (int *m, int *n,
  double *a, int *lda, int *ipiv, int *info);

int dgetrs_ (char *trans, int *n, int *nrhs, double *a,
  int *lda, int *ipiv, double *b, int *ldb, int *info);

int dgesv_ (int *n, int *nrhs, double *a, int *lda,
  int *ipiv, double *b, int *ldb, int *info);

int dgels_ (char *trans, int *m, int *n, int *nrhs, double *a, int *lda,
            double *b, int *ldb, double *work, int *lwork, int *info);

int dgetri_ (int *n,double *a, int *lda,
  int *ipiv, double *work, int *lwork, int *info);

int dpotrf_ (char *uplo, int *n,
  double *a, int *lda, int *info);

int dpotrs_ (char *uplo, int *n, int *nrhs, double *a,
  int *lda, double *b, int *ldb, int *info);

int dposv_ (char *uplo, int *n, int *nrhs, double *a, int *lda,
  double *b, int *ldb, int *info);

int dpotri_ (char *uplo, int *n, double *a, int *lda, int *info);

int dgesvd_ (char *jobu, char *jobvt, int *m, int *n, 
  double *a, int *lda, double *s, double *u, int * ldu,
  double *vt, int *ldvt, double *work, int *lwork, int *info);

/* old simple eigenvalue driver */
int dsyev_ (char *jobz, char *uplo, int *n, double *a,
  int *lda, double *w, double *work, int *lwork, int *info);

/* the fastest current eigenvalue driver (Dhillon & Parlett) */
int dsyevr_ (char *jobz, char *range, char *uplo, int *n, double *a,
  int *lda, double *vl, double *vu, int *il, int *iu, double *abstol,
  int *m, double *w, double *z, int *ldz, int *isuppz, double *work,
  int *lwork, int *iwork, int *liwork, int *info);

inline static int lapack_dgetrf (int m,
  int n, double *a, int lda, int *ipiv)
{
  int info;
  dgetrf_ (&m, &n, a, &lda, ipiv, &info);
  return info;
}

inline static int lapack_dgetrs (char trans, int n,
  int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
  int info;
  dgetrs_ (&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  return info;
}

inline static int lapack_dgesv (int  n,
  int nrhs, double *a, int lda, int *ipiv, double *b, int ldb)
{
  int info;
  dgesv_ (&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
  return info;
}

inline static int lapack_dgels (char trans, int m, int n, int nrhs,
  double *a, int lda, double *b, int ldb, double *work, int lwork)
{
  int info;
  dgels_ (&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
  return info;
}

inline static int lapack_dgetri (int n,double *a,
  int lda, int *ipiv, double *work, int lwork)
{
  int info;
  dgetri_ (&n, a, &lda, ipiv, work, &lwork, &info);
  return info;
}

inline static int lapack_dpotrf (char uplo, int n, double *a, int lda)
{
  int info;
  dpotrf_ (&uplo, &n, a, &lda, &info);
  return info;
}

inline static int lapack_dpotrs (char uplo, int n, int nrhs, double *a, int lda, double *b, int ldb)
{
  int info;
  dpotrs_ (&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
  return info;
}

inline static int lapack_dposv (char uplo, int n, int nrhs, double *a, int lda, double *b, int ldb)
{
  int info;
  dposv_ (&uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
  return info;
}

inline static int lapack_dpotri (char uplo, int n, double *a, int lda)
{
  int info;
  dpotri_ (&uplo, &n, a, &lda, &info);
  return info;
}

inline static int lapack_dgesvd (char jobu, char jobvt,
  int m, int n, double *a, int lda, double *s, double *u,
  int ldu, double *vt, int ldvt, double *work, int lwork)
{
  int info;
  dgesvd_ (&jobu, &jobvt, &m, &n, a, &lda, s,
    u, &ldu, vt, &ldvt, work, &lwork, &info);
  return info;
}
 
inline static int lapack_dsyev (char jobz, char uplo, int n, double *a,
  int lda, double *w, double *work, int lwork)
{
  int info;
  dsyev_ (&jobz, &uplo, &n, a,
    &lda, w, work, &lwork, &info);
  return info;
}

inline static int lapack_dsyevr (char jobz, char range, char uplo, int n, double *a,
  int lda, double vl, double vu, int il, int iu, double abstol, int *m, double *w,
  double *z, int ldz, int *isuppz, double *work, int lwork, int *iwork, int liwork)
{
  int info;
  dsyevr_ (&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol,
    m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
  return info;
}

#endif
