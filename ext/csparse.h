/*
 * scparse.h
 * Copyright (c) Timothy A. Davis, 2006-2009
 */

/*
 * The original code repository:
 * http://www.cise.ufl.edu/research/sparse/CSparse/
 *
 * The original code is included into the monograph:
 * T. A. Davis. Direct Methods for Sparse Linear Systems. SIAM, 2006.
 */

/*
 * Modified: T. Koziara, 2009
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

#include "../mtx.h"

#ifndef __csparse__
#define __csparse__

#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

#define CS_VER 2                    /* CSparse Version 2.2.3 */
#define CS_SUBVER 2
#define CS_SUBSUB 3
#define CS_DATE "Jan 20, 2009"     /* CSparse release date */
#define CS_COPYRIGHT "Copyright (c) Timothy A. Davis, 2006-2009"

/* --- primary CSparse routines ------------------------- */
MX *cs_add (const MX *A, const MX *B, double alpha, double beta) ;
int cs_cholsol (int order, const MX *A, double *b) ;
MX *cs_compress (const MX *T) ;
int cs_dupl (MX *A) ;
int cs_entry (MX *T, int i, int j, double x) ;
int cs_gaxpy (const MX *A, const double *x, double *y) ;
MX *cs_load (FILE *f) ;
int cs_lusol (int order, const MX *A, double *b, double tol) ;
MX *cs_multiply (const MX *A, const MX *B) ;
double cs_norm (const MX *A) ;
int cs_print (const MX *A, int brief) ;
int cs_qrsol (int order, const MX *A, double *b) ;
MX *cs_transpose (const MX *A, int values) ;
MX* cs_transpose_ext (const MX *A, MX *C) ; /* return C = A' */

/* utilities */
void *cs_calloc (int n, size_t size) ;
void *cs_free (void *p) ;
void *cs_realloc (void *p, int n, size_t size, int *ok) ;
MX *cs_spalloc (int m, int n, int nzmax, int values, int triplet) ;
MX *cs_spfree (MX *A) ;
int cs_sprealloc (MX *A, int nzmax) ;
void *cs_malloc (int n, size_t size) ;

/* --- secondary CSparse routines and data structures ----------------------- */
typedef struct cs_symbolic  /* symbolic Cholesky, LU, or QR analysis */
{
    int *pinv ;     /* inverse row perm. for QR, fill red. perm for Chol */
    int *q ;        /* fill-reducing column permutation for LU and QR */
    int *parent ;   /* elimination tree for Cholesky and QR */
    int *cp ;       /* column pointers for Cholesky, row counts for QR */
    int *leftmost ; /* leftmost[i] = min(find(A(i,:))), for QR */
    int m2 ;        /* # of rows for QR, after adding fictitious rows */
    double lnz ;    /* # entries in L for LU or Cholesky; in V for QR */
    double unz ;    /* # entries in U for LU; in R for QR */
} css ;

typedef struct cs_numeric   /* numeric Cholesky, LU, or QR factorization */
{
    MX *L ;         /* L for LU and Cholesky, V for QR */
    MX *U ;         /* U for LU, R for QR, not used for Cholesky */
    int *pinv ;     /* partial pivoting for LU */
    double *B ;     /* beta [0..n-1] for QR */
} csn ;

typedef struct cs_dmperm_results    /* cs_dmperm or cs_scc output */
{
    int *p ;        /* size m, row permutation */
    int *q ;        /* size n, column permutation */
    int *r ;        /* size nb+1, block k is rows r[k] to r[k+1]-1 in A(p,q) */
    int *s ;        /* size nb+1, block k is cols s[k] to s[k+1]-1 in A(p,q) */
    int nb ;        /* # of blocks in fine dmperm decomposition */
    int rr [5] ;    /* coarse row decomposition */
    int cc [5] ;    /* coarse column decomposition */
} csd ;

int *cs_amd (int order, const MX *A) ;
csn *cs_chol (const MX *A, const css *S) ;
csd *cs_dmperm (const MX *A, int seed) ;
int cs_droptol (MX *A, double tol) ;
int cs_dropzeros (MX *A) ;
int cs_happly (const MX *V, int i, double beta, double *x) ;
int cs_ipvec (const int *p, const double *b, double *x, int n) ;
int cs_lsolve (const MX *L, double *x) ;
int cs_ltsolve (const MX *L, double *x) ;
csn *cs_lu (const MX *A, const css *S, double tol) ;
MX *cs_permute (const MX *A, const int *pinv, const int *q, int values) ;
int *cs_pinv (const int *p, int n) ;
int cs_pvec (const int *p, const double *b, double *x, int n) ;
csn *cs_qr (const MX *A, const css *S) ;
css *cs_schol (int order, const MX *A) ;
css *cs_sqr (int order, const MX *A, int qr) ;
MX *cs_symperm (const MX *A, const int *pinv, int values) ;
int cs_updown (MX *L, int sigma, const MX *C, const int *parent) ;
int cs_usolve (const MX *U, double *x) ;
int cs_utsolve (const MX *U, double *x) ;
/* utilities */
css *cs_sfree (css *S) ;
csn *cs_nfree (csn *N) ;
csd *cs_dfree (csd *D) ;

/* --- tertiary CSparse routines -------------------------------------------- */
int *cs_counts (const MX *A, const int *parent, const int *post, int ata) ;
double cs_cumsum (int *p, int *c, int n) ;
int cs_dfs (int j, MX *G, int top, int *xi, int *pstack, const int *pinv) ;
int cs_ereach (const MX *A, int k, const int *parent, int *s, int *w) ;
int *cs_etree (const MX *A, int ata) ;
int cs_fkeep (MX *A, int (*fkeep) (int, int, double, void *), void *other) ;
double cs_house (double *x, double *beta, int n) ;
int cs_leaf (int i, int j, const int *first, int *maxfirst, int *prevleaf,
    int *ancestor, int *jleaf) ;
int *cs_maxtrans (const MX *A, int seed) ;
int *cs_post (const int *parent, int n) ;
int *cs_randperm (int n, int seed) ;
int cs_reach (MX *G, const MX *B, int k, int *xi, const int *pinv) ;
int cs_scatter (const MX *A, int j, double beta, int *w, double *x, int mark,
    MX *C, int nz) ;
csd *cs_scc (MX *A) ;
int cs_spsolve (MX *G, const MX *B, int k, int *xi, double *x,
    const int *pinv, int lo) ;
int cs_tdfs (int j, int k, int *head, const int *next, int *post,
    int *stack) ;
/* utilities */
csd *cs_dalloc (int m, int n) ;
csd *cs_ddone (csd *D, MX *C, void *w, int ok) ;
MX *cs_done (MX *C, void *w, void *x, int ok) ;
int *cs_idone (int *p, MX *C, void *w, int ok) ;
csn *cs_ndone (csn *N, MX *C, void *w, void *x, int ok) ;

#define CS_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CS_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define CS_FLIP(i) (-(i)-2)
#define CS_UNFLIP(i) (((i) < 0) ? CS_FLIP(i) : (i))
#define CS_MARKED(w,j) (w [j] < 0)
#define CS_MARK(w,j) { w [j] = CS_FLIP (w [j]) ; }
#define CS_CSC(A) (A && (A->nz == -1))
#define CS_TRIPLET(A) (A && (A->nz >= 0))

#endif
