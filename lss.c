/*
 * lss.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * -----------------------------------------------------------
 * linear system solver
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

#include <limits.h>
#include <string.h>
#include <setjmp.h>
#include <float.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "mem.h"
#include "ist.h"
#include "lss.h"

/* ========== MACROS ========== */

/* private routine */
#define PRIVATE(__type__, __name__, ...)\
static __type__ __name__ (jmp_buf __env__, __VA_ARGS__)

/* private call */
#define CALL(__name__, ...)\
  __name__ (__env__, __VA_ARGS__)

/* assertion */
#define ASSERT(__test__, __code__)\
  if (! (__test__)) longjmp (__env__, __code__)

/* error handing */
#define ONERROR(__code__)\
  jmp_buf __env__; if (((__code__) = setjmp (__env__)))

/* auxiliary cast */
#define LSSOBJ(__ptr__) ((LSSOBJ*)(__ptr__))

/* minimium of two */
#define MIN(a, b) ((a) < (b) ? (a) : (b))

/* minimium of two */
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/* absolut value */
#define ABS(a) ((a) < 0 ? -(a) : (a))

/* parameters */
#define WAVELEN 8
#define CHUNKS 256
#define INF DBL_MAX
#define NIL INT_MIN

/* ========== TYPES ========== */

typedef struct matrix MATRIX;
typedef struct vector VECTOR;
typedef struct params PARAMS;
typedef struct worksp WORKSP;
typedef struct waveop WAVEOP;
typedef struct packet PACKET;
typedef struct jkpair JKPAIR;
typedef struct symbol SYMBOL;
typedef struct precnd PRECND;
typedef struct lssobj LSSOBJ;

struct matrix /* compressed format */
{
  enum {CSC, CSR} kind;

  int dimension;

  int adjust; /* number of appended diagonal ones */

  int *columns;

  int *rows;

  double *elements;

  double **updates;
};

struct vector
{
  int dimension;

  int adjust; /* number of appended ones or zeroes */

  double *elements;
};

struct params /* user parameters */
{
  int iterations_bound;

  double relative_accuracy;

  double absolute_accuracy;

  int preconditioner;

  int decimation;

  int cutoff;

  int restart;

  int coarse_restart;

  int coarse_iterations_bound;

  double coarse_relative_accuracy;

  double coarse_absolute_accuracy;

  int presmoother;

  int postsmoother;

  int presmoothing_steps;

  int postsmoothing_steps;

  int iterations;

  double *relative_error;

  double *absolute_error;

  int levels;

  double *level_dimensions;

  double *level_nonzeros;

  double operator_complexity;

  double grid_complexity;
};

struct worksp /* GMRES workspace */
{
  VECTOR *r;  /* residual */

  VECTOR *w;  /* work */

  VECTOR **v; /* raw base of size m */

  VECTOR **z; /* preconditioned base of size m */

  double *Q;  /* H = QR, where Q is orthogonal (here Q = {c0,s0,c1,s1,...} collects cosines and sines of consecutive rotations) */

  double *R;  /* H = QR, where R is upper triangular */

  double *H;  /* upper Hessenberg Gram-Schmidt matrix of size (m + 1) * m */

  double *y;  /* solution of a minimisation problem of size m */

  double *g;  /* right hand side in a minimisation problem of size (m + 1) */

  int m;      /* restart dimension */
};

struct waveop /* wavelet restriction/prolongation operator */
{
  int decimation;

  int length;

  double wavelet [WAVELEN];
};

struct packet /* represents columns of A to be multiplied by wavelet coefficients */
{
  double *column [WAVELEN]; /* columns of A */

  int *row [WAVELEN]; /* column rows of A */

  int length [WAVELEN]; /* lengths of compressed columns of A */
};

struct jkpair
{
  int j; /* row index in restriction operator P */

  int k; /* element index of wavelet coefficient in jth row */

  JKPAIR *n;
};

struct symbol /* symbolic representation of P * A * P^T, where P is restriction */
{
  int dimension; /* number of rows in P = A->dimension / decimation */

  double *wavelet; /* wavelet coefficients */

  int length; /* wavelet length */

  PACKET *packets; /* columns of packets multiplied by wavelet coefficients represent columns of A * P^T */

  JKPAIR **jk; /* map of [0, 1, ..., A->dimension-1] to (j, k) pairs of nonzero elements in columns of P */

  int *mapping; /* auxiliary mapping used when compressing columns of P * A * P^T */

  MEM mem;
};

struct precnd /* one level data of preconditioner */
{
  MATRIX *A;

  VECTOR *b;

  VECTOR *x;

  VECTOR *r;

  WAVEOP *waveop;

  SYMBOL *symbol; /* symbolic information on computation of lower level A */

  PARAMS *user;

  MATRIX *csr; /* compressed row image of A */

  PRECND *lower;

  WORKSP *worksp;

  PARAMS *params;
};

struct lssobj
{
  /* ==== implementation ==== */

  MATRIX *A;

  PRECND *M;

  VECTOR *b;

  VECTOR *x;

  WAVEOP *waveop;

  WORKSP *worksp;

  int modified;

  /* ==== interface ==== */

  int n,
     *p,
     *i; /* aguments of LSS_Create */

  PARAMS *params;

  LSSERR error;
};

/* ========== DECLARATIONS ========== */

static void residual (MATRIX *A, VECTOR *b, VECTOR *x, VECTOR *r);

static void gmres (MATRIX*, PRECND*, VECTOR*, VECTOR*, VECTOR*, VECTOR*, VECTOR**,
             VECTOR**, double*, double*, double*, double*, double*, int, PARAMS*);

/* ========== IMPLEMENTATION ========== */

/* create index pair */
PRIVATE (JKPAIR*, jkpair, MEM *mem, int j, int k)
{
  JKPAIR *q;

  ASSERT (q = MEM_Alloc (mem), LSSERR_OUT_OF_MEMORY);

  q->j = j;
  q->k = k;

  return q;
}

/* adjust matrix dimension so that A->dimension % decimation = 0 */
PRIVATE (void, adjust, int decimation, MATRIX *A)
{
  int i;

  i = A->dimension % decimation; 

  if (i)
  {
    A->adjust = decimation - i;

    ASSERT (A->columns = realloc (A->columns, (A->dimension + A->adjust + 1) * sizeof (int)), LSSERR_OUT_OF_MEMORY);
    
    for (i = 0; i < A->adjust; i ++)
    {
      A->columns [A->dimension + i + 1] = A->columns [A->dimension + i] + 1;
    }

    A->dimension += A->adjust;

    ASSERT (A->rows = realloc (A->rows, A->columns [A->dimension] * sizeof (int)), LSSERR_OUT_OF_MEMORY);
    ASSERT (A->elements = realloc (A->elements, A->columns [A->dimension] * sizeof (double)), LSSERR_OUT_OF_MEMORY);

    for (i = A->adjust; i > 0; i --)
    {
      A->rows [A->columns [A->dimension - i]] = A->dimension - i;
      A->elements [A->columns [A->dimension - i]] = 1.0;
    }
  }
}

/* create compressed column matrix from the LSS_Create input */
PRIVATE (MATRIX*, matrix, int n, double *a, int *p, int *i)
{
  MATRIX *A;

  ASSERT (n > 0, LSSERR_INVALID_ARGUMENT); /* FIXME: need to accept compressed row format */
  ASSERT (A = malloc (sizeof (MATRIX)), LSSERR_OUT_OF_MEMORY);
  ASSERT (A->columns = malloc ((n+1) * sizeof (int)), LSSERR_OUT_OF_MEMORY);
  ASSERT (A->rows = malloc (p [n] * sizeof (int)), LSSERR_OUT_OF_MEMORY);
  ASSERT (A->elements = malloc (p [n] * sizeof (double)), LSSERR_OUT_OF_MEMORY);

  memcpy (A->columns, p, (n+1) * sizeof (int));
  memcpy (A->rows, i, p [n] * sizeof (int));
  memcpy (A->elements, a, p [n] * sizeof (double));

  A->updates = NULL;
  A->dimension = n;
  A->adjust = 0;
  A->kind = CSC;

  return A;
}

/* update matrix from user storage */
PRIVATE (void, matrixupdate, MATRIX *A, double *a)
{
  double *x, *y, **z;
  int n, m;

  n = A->dimension - A->adjust;
  m = A->kind == CSC ? A->columns [n] : A->rows [n];

  x = A->elements;
  y = x + m;
  z = A->updates;

  if (a)
  {
    for (;x < y; x ++, a ++) *x = *a;
  }
  else
  {
    for (;x < y; x ++, z ++) *x = **z;
  }

  if (A->kind == CSR)
  {
    int *r = A->rows, *e = r + n;

    for (x = A->elements; r < e; r ++) ASSERT (x [*r] != 0.0, LSSERR_ZERO_ON_DIAGONAL);
  }
}

/* free matrix */
static void matrixfree (MATRIX *A)
{
  free (A->columns);
  free (A->rows);
  free (A->elements);
  free (A->updates);
  free (A);
}

/* create compressed row image of input matrix */
PRIVATE(MATRIX*, csr, MATRIX *A)
{
  MATRIX *B;
  int i, *n, *r, *e, *arows, *brows, *acolumns, *bcolumns, dimension, nnz;
  double *aelements, **bupdates, *a, **b;

  if (A->kind == CSR) return A;

  ASSERT (B = malloc (sizeof (MATRIX)), LSSERR_OUT_OF_MEMORY);
  ASSERT (B->columns = malloc (A->columns[A->dimension] * sizeof (int)), LSSERR_OUT_OF_MEMORY);
  ASSERT (B->rows = calloc (A->dimension + 1, sizeof (int)), LSSERR_OUT_OF_MEMORY);
  ASSERT (B->elements = malloc (A->columns[A->dimension] * sizeof (double)), LSSERR_OUT_OF_MEMORY); 
  ASSERT (B->updates = malloc (A->columns[A->dimension] * sizeof (double*)), LSSERR_OUT_OF_MEMORY); 
  ASSERT (n = calloc (A->dimension,  sizeof (int)), LSSERR_OUT_OF_MEMORY); 
  B->kind = CSR;
  B->dimension = A->dimension;
  B->adjust = 0;

  nnz = A->columns[A->dimension];
  dimension = A->dimension;
  arows = A->rows;
  acolumns = A->columns;
  aelements = A->elements;
  brows = B->rows;
  bcolumns = B->columns;
  bupdates = B->updates;

  for (i = 0; i < nnz; i ++) brows [arows [i] + 1] ++;

  for (i = 1; i <= dimension; i ++)  brows [i] += brows [i-1];

  for (i = 0; i < dimension; i ++)
  {
    for (r = &arows[acolumns [i]], e = &arows [acolumns[i+1]], a = &aelements [acolumns [i]]; r < e; r ++, a ++)
    {
      bcolumns [brows [*r] + n [*r]] = i;
      bupdates [brows [*r] + n [*r]] = a;
      n [*r] ++;
    }
  }

  /* swap diagonal row entry with first entry for each row */
  for (i = 0; i < dimension; i ++)
  {
    for (r = &bcolumns[brows [i]], e = &bcolumns [brows [i+1]], b = &bupdates[brows [i]]; r < e; r ++, b ++)
    {
      if (*r == i)
      {
	*r = bcolumns[brows[i]];
        bcolumns[brows[i]] = i;

	a = *b;
	*b = bupdates[brows[i]];
	bupdates[brows[i]] = a;

	break;
      }
    }
  }

  free (n);

  return B;
}

/* create vector of size n */
PRIVATE (VECTOR*, vectoralloc, int n)
{
  VECTOR *v;

  ASSERT (v = malloc (sizeof (VECTOR)), LSSERR_OUT_OF_MEMORY);
  ASSERT (v->elements = calloc (n, sizeof (double)), LSSERR_OUT_OF_MEMORY);
  v->dimension = n;
  v->adjust = 0;
  return v;
}

/* create vector wrapper */
PRIVATE (VECTOR*, vector, double *z, int n, int a)
{
  VECTOR *v;

  ASSERT (v = malloc (sizeof (VECTOR)), LSSERR_OUT_OF_MEMORY);
  ASSERT (v->elements = malloc ((n + a) * sizeof (double)), LSSERR_OUT_OF_MEMORY);
  memcpy (v->elements, z, n * sizeof (double));
  v->dimension = n + a;
  v->adjust = a;

  return v;
}

/* update vector form user storage */
PRIVATE (void, vectorupdate, VECTOR *v, double *z, double a)
{
  double *x, *y;

  x = v->elements;
  y = x + v->dimension - v->adjust;

  for (;x < y; x ++, z ++) *x = *z;

  for (int i = v->adjust; i > 0; i --, x ++) *x = a;
}

/* export vector to user storage */
PRIVATE (void, vectorexport, VECTOR *v, double *z)
{
  double *x, *y;

  x = v->elements;
  y = x + v->dimension - v->adjust;

  for (;x < y; x ++, z ++) *z = *x;
}

/* free vector */
static void vectorfree (VECTOR *v)
{
  free (v->elements);
  free (v);
}

/* create wavelet operator */
PRIVATE (WAVEOP*, waveop, int precon, int decimation)
{
  WAVEOP *op;

  ASSERT (op = malloc (sizeof (WAVEOP)), LSSERR_OUT_OF_MEMORY);
  op->decimation = decimation;

  switch (precon)
  {
    case 1: /* Daubechies-2 (Haar) */
      op->length = 2;
      op->wavelet [0] = 1.0 / sqrt (2.0);
      op->wavelet [1] = op->wavelet [0];
    break;
    case 2: /* Daubechies-4 */
    {
      double _s3 = sqrt (3.0),
	     _4_s2 = 4.0 * sqrt (2.0);

      op->length = 4;
      op->wavelet [0] = (1.0 + _s3)/_4_s2;
      op->wavelet [1] = (3.0 + _s3)/_4_s2;
      op->wavelet [2] = (3.0 - _s3)/_4_s2;
      op->wavelet [3] = (1.0 - _s3)/_4_s2;
    }
    break;
    case 3: /* Daubechies-6 */
    {
      double _16_s2 = 16.0 * sqrt (2.0),
	     _s10 = sqrt (10.0),
	     _s5_2_s10 = sqrt (5.0 + 2.0*_s10);

      op->length = 6;
      op->wavelet [0] = (1.0 + _s10 - _s5_2_s10) / _16_s2;
      op->wavelet [1] = (1.0 + _s10 - 3.0*_s5_2_s10) / _16_s2; 
      op->wavelet [2] = (1.0 - 2.0*_s10 - 2.0*_s5_2_s10) / _16_s2;
      op->wavelet [3] = (1.0 - 2.0*_s10 + 2.0*_s5_2_s10) / _16_s2; 
      op->wavelet [4] = (1.0 + _s10 + 3.0*_s5_2_s10) / _16_s2; 
      op->wavelet [5] = (1.0 + _s10 + _s5_2_s10) / _16_s2;
    }
    break;
  };

  return op;
}

/* free wavelet operator */
static void waveopfree (WAVEOP *op)
{
  free (op);
}

/* create default parameters */
static PARAMS* defparams ()
{
  PARAMS *params;

  if (!(params = malloc (sizeof (PARAMS)))) return NULL;

  params->iterations_bound = 1000;
  params->relative_accuracy = 1E-6;
  params->absolute_accuracy = 1E-3;
  params->preconditioner = 1;
  params->decimation = 4;
  params->cutoff = 16;
  params->restart = 32;
  params->coarse_restart = 16;
  params->coarse_iterations_bound = 16;
  params->coarse_relative_accuracy = 1E-6;
  params->coarse_absolute_accuracy = 1E-3;
  params->presmoother = 1;
  params->postsmoother = 1;
  params->presmoothing_steps = 1;
  params->postsmoothing_steps = 1;

  params->iterations = 0;
  if (!(params->relative_error = calloc (params->iterations_bound, sizeof (double)))) return NULL;
  if (!(params->absolute_error = calloc (params->iterations_bound, sizeof (double)))) return NULL;
  params->levels = 0;
  params->level_dimensions = NULL;
  params->level_nonzeros = NULL;
  params->operator_complexity = 0;
  params->grid_complexity = 0;

  return params;
}

/* create lowest level parameters */
static PARAMS* lowparams (PARAMS *user)
{
  PARAMS *params;

  if (!(params = calloc (1, sizeof (PARAMS)))) return NULL;

  params->iterations_bound = user->coarse_iterations_bound;
  params->relative_accuracy = user->coarse_relative_accuracy;
  params->absolute_accuracy = user->coarse_absolute_accuracy;

  return params;
}

/* udate lowest level parameters */
static void lowparamsupdate (PARAMS *params, PARAMS *user)
{
  params->iterations_bound = user->coarse_iterations_bound;
  params->relative_accuracy = user->coarse_relative_accuracy;
  params->absolute_accuracy = user->coarse_absolute_accuracy;
}

/* free parameters */
static void paramsfree (PARAMS *params)
{
  free (params->relative_error);
  free (params->absolute_error);
  free (params->level_dimensions);
  free (params->level_nonzeros);
  free (params);
}

/* create workspace */
PRIVATE (WORKSP*, workspace, int n, int m)
{
  WORKSP *worksp;

  ASSERT (worksp = malloc (sizeof (WORKSP)), LSSERR_OUT_OF_MEMORY);
  ASSERT (worksp->v = malloc (m * sizeof (VECTOR*)), LSSERR_OUT_OF_MEMORY);
  ASSERT (worksp->z = malloc (m * sizeof (VECTOR*)), LSSERR_OUT_OF_MEMORY);
  ASSERT (worksp->Q = malloc (m * sizeof (double [2])), LSSERR_OUT_OF_MEMORY);
  ASSERT (worksp->R = malloc ((m*m + m) * sizeof (double)), LSSERR_OUT_OF_MEMORY);
  ASSERT (worksp->H = malloc ((m*m + m) * sizeof (double)), LSSERR_OUT_OF_MEMORY);
  ASSERT (worksp->y = malloc (m * sizeof (double)), LSSERR_OUT_OF_MEMORY);
  ASSERT (worksp->g = malloc ((m + 1) * sizeof (double)), LSSERR_OUT_OF_MEMORY);
  worksp->r = CALL (vectoralloc, n);
  worksp->w = CALL (vectoralloc, n);
  worksp->m = m;

  for (m --; m >= 0; m --)
  {
    worksp->v [m] = CALL (vectoralloc, n);
    worksp->z [m] = CALL (vectoralloc, n);
  }

  return worksp;
}

/* free workspace */
static void workspfree (WORKSP *worksp)
{
  free (worksp->Q);
  free (worksp->R);
  free (worksp->H);
  free (worksp->y);
  free (worksp->g);

  vectorfree (worksp->r);
  vectorfree (worksp->w);

  for (worksp->m --; worksp->m >= 0; worksp->m --)
  {
    vectorfree (worksp->v [worksp->m]);
    vectorfree (worksp->z [worksp->m]);
  }

  free (worksp->v);
  free (worksp->z);
  free (worksp);
}

/* create symbolic representation of P * A * P^T, where P is restriction */
PRIVATE (MATRIX*, symbolic, WAVEOP *waveop, MATRIX *A, SYMBOL **symbol)
{
  int i, j, k, l, d, n, m, s;
  PACKET *pac, *pacend;
  JKPAIR *pair, **jk;
  int *column, *row;
  double *element;
  MEM *jkmem, mem;
  ISET *set, *it;
  MATRIX *B;

  ASSERT (*symbol = malloc (sizeof (SYMBOL)), LSSERR_OUT_OF_MEMORY); /* symbolic representaion of P * A * P^T */
  ASSERT (B = malloc (sizeof (MATRIX)), LSSERR_OUT_OF_MEMORY); /* numeric representation of P * A * P^T */

  (*symbol)->dimension = A->dimension / waveop->decimation;
  (*symbol)->wavelet = waveop->wavelet;
  (*symbol)->length = waveop->length;

  B->dimension = (*symbol)->dimension;
  B->adjust = 0;
  B->updates = NULL;
  B->kind = CSC;

  s = A->columns [A->dimension];
  ASSERT (B->rows = malloc (s * sizeof (int)), LSSERR_OUT_OF_MEMORY); /* initial guess (likely overallocated) */
  ASSERT (B->columns = malloc ((B->dimension+1) * sizeof (int)), LSSERR_OUT_OF_MEMORY);
  ASSERT ((*symbol)->packets = malloc ((*symbol)->dimension * sizeof (PACKET)), LSSERR_OUT_OF_MEMORY);
  ASSERT ((*symbol)->jk = calloc (A->dimension, sizeof (JKPAIR*)), LSSERR_OUT_OF_MEMORY);
  ASSERT ((*symbol)->mapping = calloc (B->dimension, sizeof (int)), LSSERR_OUT_OF_MEMORY);
  MEM_Init (&(*symbol)->mem, sizeof (JKPAIR), waveop->length * A->dimension);
  MEM_Init (&mem, sizeof (ISET), CHUNKS);

  m = (*symbol)->dimension;
  d = waveop->decimation;
  l = waveop->length;
  n = A->dimension;
  jkmem = &(*symbol)->mem;
  element = A->elements;
  column = A->columns;
  jk = (*symbol)->jk;
  row = A->rows;

  for (pac = (*symbol)->packets, j = 0; j < m; pac ++, j ++)
  {
    for (k = 0; k < l; k ++)
    {
      i = (j*d + k) % n;

      /* structure of X = A * P^T */
      pac->column [k] = &element [column [i]];
      pac->row [k] = &row [column [i]];
      pac->length [k] = column [i+1] - column [i];

      /* structure of P * X */
      pair = CALL (jkpair, jkmem, j, k);
      pair->n = jk [i];
      jk [i] = pair;
    }
  }

  row = B->rows;
  column = B->columns;
  pac = (*symbol)->packets;
  pacend = pac + m;

  for (*column = 0, column ++; pac < pacend; pac ++, column ++) /* j */
  {
    *column = *(column-1); /* column [j] = SUM (i < j) column [i] */

    set = NULL;

    for (k = 0; k < l; k ++)
    {
      for (i = 0; i < pac->length [k]; i ++)
      {
	for (pair = jk [pac->row [k][i]]; pair; pair = pair->n)
	{
	  ISET_Insert (&mem, &set, pair->j);  /* product of P with jth column of A * P^T has nonzero entry at pair->j */
	}
      }
    }

    for (it = ISET_First (set); it; it = ISET_Next (it))
    {
      if (*column > s) /* not enough memory for rows */
      { s *= 2; ASSERT (B->rows = realloc (B->rows, s * sizeof (int)), LSSERR_OUT_OF_MEMORY); row = B->rows; }

      row [*column] = it->value;
      (*column) ++;
    }

    ISET_Free (&mem, &set);
  }

  ASSERT (B->rows = realloc (B->rows, B->columns [B->dimension] * sizeof (int)), LSSERR_OUT_OF_MEMORY); /* trim */
  ASSERT (B->elements = malloc (B->columns [B->dimension] * sizeof (double)), LSSERR_OUT_OF_MEMORY);
  MEM_Release (&mem);

  return B;
}

/* compute numeric representation of B = P * A * P^T */
PRIVATE (void, numeric, SYMBOL *symbol, MATRIX *B)
{
  double *wavelet, *elements, *elemend, *ele, *e;
  int *mapping, *rows, *row, *rowend, *column;
  PACKET *pac, *pacend;
  JKPAIR *pair, **jk;
  int k, l;

  l = symbol->length;
  pac = symbol->packets;
  pacend = pac + symbol->dimension;
  wavelet = symbol->wavelet;
  mapping = symbol->mapping;
  elements = B->elements;
  column = B->columns;
  elemend = elements + column [B->dimension - B->adjust];
  rows = B->rows;
  jk = symbol->jk;

  for (e = elements; e < elemend; e ++) *e = 0.0;

  for (; pac < pacend; pac ++, column ++)
  {
    row = &rows [*column];
    rowend = &rows [*(column+1)];
    e = &elements [*column];

    for (k = 0; row < rowend; k ++, row ++) mapping [*row] = k; /* (0, 1, ..., B->dimension) to compressed column mapping */

    for (k = 0; k < l; k ++)
    {
      row = pac->row [k];
      rowend = row + pac->length [k];
      ele = pac->column [k];

      for (; row < rowend; row ++, ele ++)
      {
	for (pair = jk [*row]; pair; pair = pair->n)
	{
	  e [mapping [pair->j]] += wavelet [pair->k] * wavelet [k] * (*ele);
	}
      }
    }
  }
}

/* free symbolic data */
static void symbolfree (SYMBOL *symbol)
{
  free (symbol->packets);
  free (symbol->jk);
  free (symbol->mapping);
  MEM_Release (&symbol->mem);
  free (symbol);
}

/* create symbolic representation of preconditioner */
PRIVATE (PRECND*, precnd, WAVEOP *waveop, MATRIX *A, PARAMS *params)
{
  PRECND *M;

  ASSERT (M = malloc (sizeof (PRECND)), LSSERR_OUT_OF_MEMORY);

  M->user = params;

  if (A->dimension <= params->cutoff)
  {
    M->A = A;
    M->b = CALL (vectoralloc, A->dimension);
    M->x = CALL (vectoralloc, A->dimension);
    M->csr = NULL;
    M->waveop = NULL;
    M->symbol = NULL;
    M->lower = NULL;
    M->worksp = CALL (workspace, A->dimension, params->coarse_restart);
    ASSERT (M->params = lowparams (params), LSSERR_OUT_OF_MEMORY);
  }
  else
  {
    CALL (adjust, waveop->decimation, A);

    M->A = A;
    M->b = CALL (vectoralloc, A->dimension);
    M->x = CALL (vectoralloc, A->dimension);
    M->r = CALL (vectoralloc, A->dimension);
    M->b->adjust = A->adjust;
    M->x->adjust = A->adjust;
    M->csr = CALL (csr, A);
    M->waveop = waveop;
    M->lower = CALL (precnd, waveop, CALL (symbolic, waveop, A, &M->symbol), params);
    M->worksp = NULL;
    M->params = NULL;
  }

  return M;
}

/* update preconditioner */
PRIVATE (void, precndupdate, PRECND *M, PARAMS *params)
{
  if (M->lower)
  {
    CALL (matrixupdate, M->csr, NULL);
    CALL (numeric, M->symbol, M->lower->A);
    CALL (precndupdate, M->lower, params);
  }
  else 
  {
    lowparamsupdate (M->params, params);

    if (M->worksp->m != params->coarse_restart)
    {
      workspfree (M->worksp);
      M->worksp = CALL (workspace, M->A->dimension, params->coarse_restart);
    }
  }
}

/* preconditioner levels paramameters setup recursion */
PRIVATE (void, precndlevelsdo, PRECND *M, PARAMS *params)
{
  if (M)
  {
    params->levels ++;
    ASSERT (params->level_dimensions = realloc (params->level_dimensions, params->levels * sizeof (double)), LSSERR_OUT_OF_MEMORY);
    ASSERT (params->level_nonzeros = realloc (params->level_nonzeros, params->levels * sizeof (double)), LSSERR_OUT_OF_MEMORY);
    params->level_dimensions [params->levels - 1] = M->A->dimension;
    params->level_nonzeros [params->levels - 1] = M->A->columns [M->A->dimension];

    CALL (precndlevelsdo, M->lower, params);
  }
}

/* set up preconditioner levels paramameters */
PRIVATE (void, precndlevels, PRECND *M, PARAMS *params)
{
  params->levels = 0;
  CALL (precndlevelsdo, M, params);
}

/* sum of nonzero elements */
static double precndnnz (PRECND *M)
{
  if (M) return ((double)M->A->columns [M->A->dimension] + precndnnz (M->lower));
  else return 0.0;
}

/* ratio of total to top level nonzeros */
static double operator_complexity (PRECND *M)
{
  return (precndnnz (M) / (double)M->A->columns [M->A->dimension]);
}

/* sum of dimensions */
static double precndim (PRECND *M)
{
  if (M) return ((double)M->A->dimension + precndim (M->lower));
  else return 0.0;
}

/* ratio of total to top level dimensions */
static double grid_complexity (PRECND *M)
{
  return (precndim (M) / (double)M->A->dimension);
}

/* free preconditioner */
static void precndfree (PRECND *M)
{
  if (M->lower)
  {
    precndfree (M->lower);
    symbolfree (M->symbol);
    matrixfree (M->csr);
    vectorfree (M->r);
  }
  else
  {
    workspfree (M->worksp);
    paramsfree (M->params);
  }

  matrixfree (M->A);
  vectorfree (M->b);
  vectorfree (M->x);
  free (M);
}

/* y = x */
static void copy (VECTOR *x, VECTOR *y)
{
  double *a, *b, *c;

  a = x->elements;
  b = a + x->dimension;
  c = y->elements;

  for (;a < b; a ++, c ++) *c = *a;
}

/* restriction y = P * x */
static void restriction (WAVEOP *waveop, VECTOR *x, VECTOR *y)
{
  double *wavelet, *a, *b;
  int j, k, l, d, n, m;

  wavelet = waveop->wavelet;
  d = waveop->decimation;
  l = waveop->length;
  n = x->dimension;
  m = y->dimension - y->adjust;
  a = y->elements;
  for (b = a + m; a < b; a ++) *a = 0.0;
  a = x->elements;
  b = y->elements;

  for (j = 0; j < m; j ++, b ++)
  {
    for (k = 0; k < l; k ++) (*b) += wavelet [k] * a [(j*d + k) % n];
  }

  for (j = y->adjust; j > 0; j --, b ++) *b = 1.0;
}

/* prolongation y = y + P^T * x, where P is restriction */
static void prolongation (WAVEOP *waveop, VECTOR *x, VECTOR *y)
{
  double *wavelet, *a, *b;
  int j, k, l, d, n, m;

  wavelet = waveop->wavelet;
  d = waveop->decimation;
  l = waveop->length;
  m = x->dimension - x->adjust;
  n = y->dimension;
  a = x->elements;
  b = y->elements;

  for (j = 0; j < m; j ++, a ++)
  {
    for (k = 0; k < l; k ++) b [(j*d + k) % n] += wavelet [k] * (*a);
  }
}

/* Gauss-Seidel relaxation for A x = b */
static void gaussseidel (int steps, MATRIX *A, VECTOR *b, VECTOR *x)
{
  double *xelements, *belements, *aelements;
  double *xx, *bb, *dd, *aa, ss;
  int *cc, *ce, *rr, *re;
  int *arows, *acolumns;

  arows = A->rows;
  acolumns = A->columns;
  aelements = A->elements;
  xelements = x->elements;
  belements = b->elements;
  re = arows + A->dimension;

  for (; steps > 0; steps --)
  {
    for (xx = xelements, bb = belements, rr = arows; rr < re; xx ++, bb ++, rr ++)
    {
      for (dd = &aelements [*rr], aa = dd+1, cc = &acolumns [*rr+1],
	   ce = &acolumns [*(rr+1)], ss = 0.0; cc < ce; aa ++, cc ++) ss += (*aa) * xelements [*cc];

      *xx = (*bb - ss) / *dd;
    }
  }
}

/* Jacobi relaxation steps for A x = b */
static void jacobi (int steps, MATRIX *A, VECTOR *b, VECTOR *x) /* FIXME: optimize */
{
  double *xelements, *belements, *aelements, *selements, *delements;
  double *xx, *bb, *dd, *aa, *ss, *se;
  int *cc, *ce, *rr, *re;
  int *arows, *acolumns;
  VECTOR *s, *d;
  int err, j;

  ONERROR (err) return;

  s = CALL (vectoralloc, A->dimension);
  d = CALL (vectoralloc, A->dimension);

  arows = A->rows;
  acolumns = A->columns;
  aelements = A->elements;
  xelements = x->elements;
  belements = b->elements;
  selements = s->elements;
  delements = d->elements;
  ce = acolumns + A->dimension;
  se = selements + A->dimension;

  for (; steps > 0; steps --)
  {
    for (ss = selements; ss < se; ss ++) *ss = 0.0;

    for (cc = acolumns, ss = selements, dd = delements, xx = xelements, j = 0; cc < ce; cc ++, ss ++, dd ++, xx ++, j ++)
    {
      for (aa = &aelements [*cc], rr = &arows [*cc], re = &arows [*(cc+1)]; rr < re; aa ++, rr ++)
      {
	if (*rr != j) selements[*rr] += (*aa) * (*xx);
	else *dd = *aa;
      }
    }

    for (ss = selements, dd = delements, xx = xelements, bb = belements; ss < se; ss ++, dd ++, xx ++, bb ++) *xx = (*bb - *ss) / *dd;
  }
}

/* Kaczmarz relaxation steps for A x = b */
static void kaczmarz (int steps, MATRIX *A, VECTOR *b, VECTOR *x) /* FIXME: optimize, randomize */
{
  double *xelements, *belements, *aelements;
  double *bb, *aa, ss1, ss2;
  int *cc, *ce, *rr, *re;
  int *arows, *acolumns;
  int i, j, n;

  arows = A->rows;
  acolumns = A->columns;
  aelements = A->elements;
  xelements = x->elements;
  belements = b->elements;
  re = arows + A->dimension;
  n = A->dimension;

  srand ((unsigned) time (NULL));

  for (; steps > 0; steps --)
  {
    for (i = 0; i < n; i ++)
    {
      j = rand () % n;
      bb = &belements [j];
      rr = &arows [j];

      for (aa = &aelements [*rr], cc = &acolumns [*rr], ce = &acolumns [*(rr+1)], ss1 = 0.0, ss2 = 0.0; cc < ce; aa ++, cc ++)
      {
	ss1 += (*aa) * (*aa);
	ss2 += (*aa) * xelements [*cc];
      }

      ss1 = (*bb - ss2) / ss1;

      for (aa = &aelements [*rr], cc = &acolumns [*rr], ce = &acolumns [*(rr+1)]; cc < ce; aa ++, cc ++)
      {
        xelements [*cc] += ss1 * (*aa);
      }
    }
  }
}

/* relaxation steps for A x = b */
static void smoothing (int smoother, int steps, MATRIX *csc, MATRIX *csr, VECTOR *b, VECTOR *x)
{
  switch (smoother)
  {
    case 1:
      gaussseidel (steps, csr, b, x);
      break;
    case 2:
      jacobi (steps, csc, b, x);
      break;
    case 3:
      kaczmarz (steps, csr, b, x);
      break;
  }
}

/* approximately solve M x = b */
static void prevec (PRECND *M, VECTOR *b, VECTOR *x)
{
  if (M == NULL) { copy (b, x); return; }

  if (b) copy (b, M->b);

  if (M->lower)
  {
    smoothing (M->user->presmoother, M->user->presmoothing_steps, M->A, M->csr, M->b, M->x);

    residual (M->A, M->b, M->x, M->r);

    restriction (M->waveop, M->r, M->lower->b);

    prevec (M->lower, NULL, NULL);

    prolongation (M->waveop, M->lower->x, M->x);

    smoothing (M->user->postsmoother, M->user->postsmoothing_steps, M->A, M->csr, M->b, M->x);
  }
  else
  {
    WORKSP *s = M->worksp;

    gmres (M->A, NULL, M->b, M->x, s->r, s->w, s->v, s->z, s->Q, s->R, s->H, s->y, s->g, s->m, M->params);
  }

  if (x) copy (M->x, x);
}

/* y = A x */
static void matvec (MATRIX *A, VECTOR *x, VECTOR *y)
{
  double *ye, *ae, *xe, *yele, *yend, *xele, *aele;
  int *col, *row, *end, *acol, *aend, *arow;

  yele = y->elements;
  yend = yele + y->dimension;
  xele = x->elements;
  aele = A->elements;
  acol = A->columns;
  aend = acol + A->dimension;
  arow = A->rows;

  for (ye = yele; ye < yend; ye ++) *ye = 0.0;

  for (col = acol, ae = aele, xe = xele; col < aend; col ++, xe ++)
  {
    for (row = &arow [*col], end = &arow [*(col+1)]; row < end; row ++, ae ++) yele [*row] += (*ae) * (*xe);
  }
}

/* r = b - A x */
static void residual (MATRIX *A, VECTOR *b, VECTOR *x, VECTOR *r)
{
  double *bb = b->elements, *rr = r->elements;

  matvec (A, x, r);

  for (int n = r->dimension; n > 0; n --, rr ++, bb ++) *rr = *bb - *rr;
}

/* d = a + b * c */
static void vecaddmul (VECTOR *a, double b, VECTOR *c, VECTOR *d)
{
  double *x = a->elements, *y = c->elements, *z = d->elements;

  for (int n = a->dimension; n > 0; n --, x ++, y ++, z ++) *z = (*x) + b*(*y);
}

/* c = b * a */
static void vecmul (VECTOR *a, double b, VECTOR *c)
{
  double *x = a->elements, *y = c->elements;

  for (int n = a->dimension; n > 0; n --, x ++, y ++) *y = b*(*x);
}

/* dot product <a, b> */
static double vecvec (VECTOR *a, VECTOR *b)
{
  double *x = a->elements, *y = b->elements, dot = 0.0;

  for (int n = a->dimension; n > 0; n --, x ++, y ++) dot +=  (*y)*(*x);

  return dot;
}

/* |a| = sqrt (<a, a>) */
static double veclen (VECTOR *a)
{
  return sqrt (vecvec (a, a));
}

#define H(i, j) H [(j)*(m+1) + (i)]
#define R(i, j) R [(j)*(m+1) + (i)]

/* find min [y] |beta * e1 - H y|, return residual length |r| = |Ax - b| */
static double minimise (double beta, double *H, double *Q, double *R, double *y, double *g, int m, int j)
{
  double d, *q, *h, *r, *e;
  int k, l;

  if (j == 0) { g [0] = beta; for (q = g + 1, e = g + m + 1; q < e; q ++) *q = 0.0; } /* g = beta * e1 */

  for (h = &H (0, j), r = &R (0, j), e = r + j + 2; r < e; h ++, r ++) *r = *h; /* R(*, j) = H(*, j) */

  for (q = Q, r = &R (0, j), e = r + j; r < e; q += 2, r ++) /* R(*, j) = Q H (*, j) (apply previous rotations) */
  {
    d = r [0];

    r [0] = q [0] * r [0] - q [1] * r [1];
    r [1] = q [1] * d     + q [0] * r [1]; 
  }

  d = sqrt (r[0]*r[0] + r[1]*r[1]);

  q [0] =  r [0] / d; /* cos */
  q [1] = -r [1] / d; /* sin */

  r [0] = q [0] * r [0] - q [1] * r [1]; /* R (j+1, j) */

  d = g [j];

  g [j]   = q [0] * g [j] - q [1] * g [j+1]; /* g = Q beta * e1 */
  g [j+1] = q [1] * d     + q [0] * g [j+1]; 

  for (k = j; k >= 0; k --) /* min (y) |beta * e1 - H y| = min (y) |g - R y| (back substitution) */
  {
    for (d = g [k], l = j; l > k; l --) d -= R (k, l) * y [l];

    y [k] = d / R (k, k);
  }

  return fabs (g [j+1]); /* residual */

  /* TODO: use LAPACK's DLAIC1 for incremental condition number estimation of H,R */
  /* TODO: use LAPACK's DHSEQR to compute eigenvalues of Hessenberg matrix */
}

/* generalized minimal residual with extensions */
static void gmres (MATRIX *A, PRECND *M, VECTOR *b, VECTOR *x, VECTOR *r, VECTOR *w, VECTOR **v,
        VECTOR **z, double *Q, double *R, double *H, double *y, double *g, int m, PARAMS  *params)
{
  double alpha,
         beta,
	 gamma,
	 omega,
	 epsilon,
	 *relative_error,
	 *absolute_error;
  int iteration,
      maxiter,
      i, j;

  maxiter = params->iterations_bound;
  epsilon = params->relative_accuracy;
  gamma = params->absolute_accuracy;
  relative_error = params->relative_error;
  absolute_error = params->absolute_error;

  iteration = 0;

  do
  {
    residual (A, b, x, r); /* r = b - A x */
    beta = veclen (r);
    vecmul (r, 1.0/beta, v [0]); /* v[0] = r / |r| (first base vector) */

    if (!iteration) omega = beta;

    for (j = 0; j < m; j ++)
    {
      prevec (M, v [j], z [j]); /* z[j] = preconditioned jth base vector */
      matvec (A, z [j], w); /* w = ((A * M ** (-1)) ** (j+1)) * r */

      for (i = 0; i <= j; i ++) /* orthogonalize w so that <w, v[i]> = 0 */
      {
	H (i, j) = vecvec (w, v [i]);
	vecaddmul (w, -H (i, j), v [i], w); /* subtract components along previous base vectors (modified Gramm-Schmidt orthogonalisation) */
      }

      H (j+1, j) = veclen (w);
      if (j < (m-1)) vecmul (w, 1.0/H (j+1, j), v [j+1]); /* v[j+1] = w / |w| (next base vector) */

      alpha = minimise (beta, H, Q, R, y, g, m, j); /* solve min (y) [beta * e1 - H y] */

      if (relative_error) relative_error [iteration] = alpha / omega;
      if (absolute_error) absolute_error [iteration] = alpha;

      iteration ++;

      if ((alpha <= epsilon * omega && alpha <= gamma) || iteration >= maxiter) break;
    }

    for (i = 0; i <= j && i < m; i ++) vecaddmul (x, y [i], z [i], x); /* x = x + sum (i) [y [i] * z [i]] */
  }
  while ((alpha > epsilon * omega || alpha > gamma) && iteration < maxiter);

  params->iterations = iteration;
}

/* create internal data */
PRIVATE (void, create, LSSOBJ *lss, double *a, double *x, double *b)
{
  lss->A = CALL (matrix, lss->n, a, lss->p, lss->i);

  if (lss->params->preconditioner)
  {
    lss->waveop = CALL (waveop, lss->params->preconditioner, lss->params->decimation);
    lss->M = CALL (precnd, lss->waveop, lss->A, lss->params);
    lss->params->operator_complexity = operator_complexity (lss->M);
    lss->params->grid_complexity = grid_complexity (lss->M);
    CALL (precndlevels, lss->M, lss->params);
  }
  else
  {
    lss->waveop = NULL;
    lss->M = NULL;
  }

  lss->x = CALL (vector, x, ABS (lss->n), lss->A->adjust);
  lss->b = CALL (vector, b, ABS (lss->n), lss->A->adjust);
  lss->worksp = CALL (workspace, lss->A->dimension, lss->params->restart);
  lss->modified = 0;
}

/* update internal data */
PRIVATE (void, update, LSSOBJ *lss, double *a, double *x, double *b)
{
  CALL (matrixupdate, lss->A, a);
  CALL (vectorupdate, lss->x, x, 0.0);
  CALL (vectorupdate, lss->b, b, 1.0);
  if (lss->M) CALL (precndupdate, lss->M, lss->params);
}

/* destroy internal data */
static void destroy (LSSOBJ *lss)
{
  if (lss->A)
  {
    if (lss->M)
    {
      waveopfree (lss->waveop);
      precndfree (lss->M); /* also frees lss->A */
    }
    else matrixfree (lss->A);

    vectorfree (lss->b);
    vectorfree (lss->x);
    workspfree (lss->worksp);
  }
}

/* solve linear system */
PRIVATE (void, solve, LSSOBJ *lss, double *x)
{
  WORKSP *s = lss->worksp;

  /* solve linear system */
  gmres (lss->A, lss->M, lss->b, lss->x, s->r, s->w, s->v,
    s->z, s->Q, s->R, s->H, s->y, s->g, s->m, lss->params);

  /* copy solution to user storage */
  CALL (vectorexport, lss->x, x);

  /* test convergence */
  ASSERT (lss->params->iterations < lss->params->iterations_bound, LSSERR_LACK_OF_CONVERGENCE);
}

/* test for modification */
static int modified (LSSOBJ *lss)
{
  return lss->modified;
}

/* ========== INTERFACE ========== */

/* n: dimension of n x n system A matrix
 * n < 0 => A in compressed row format (CSR)
 * n > 0 => A in compressed column formt (CSC)
 * a: system A matrix values
 * p: pointers to rows (CSR) or columns (CSC) of size (n+1)
 * i: indices of columns (CSR) or rows (CSC) corresponding to 'a';
 * create linear system solver */
void* LSS_Create (int n, double *a, int *p, int *i)
{
  LSSOBJ *lss;

  if (!(lss = malloc (sizeof (LSSOBJ)))) return NULL;
  if (!(lss->params = defparams ())) return NULL; /* defaults */
  lss->error = LSSERR_NONE;
  lss->n = n;
  lss->p = p;
  lss->i = i;

  lss->A = NULL;
  lss->M = NULL;
  lss->b = NULL;
  lss->x = NULL;
  lss->waveop = NULL;
  lss->worksp = NULL;
  lss->modified = 1;

  return lss;
}

/* lss: returned by LSS_Create
 * p: parameter
 * v: value;
 * set parameter */
LSSERR LSS_Set (void *lss, LSSPAR p, double v)
{
  if (!lss) return LSSERR_OUT_OF_MEMORY;

  ONERROR (LSSOBJ(lss)->error) return LSSOBJ(lss)->error; /* an error has occured */

  switch (p)
  {
    case LSS_ITERATIONS_BOUND:
      ASSERT (v >= 1, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->iterations_bound = (int)v;
      free (LSSOBJ(lss)->params->relative_error);
      free (LSSOBJ(lss)->params->absolute_error);
      ASSERT (LSSOBJ(lss)->params->relative_error = calloc (LSSOBJ(lss)->params->iterations_bound, sizeof (double)), LSSERR_OUT_OF_MEMORY);
      ASSERT (LSSOBJ(lss)->params->absolute_error = calloc (LSSOBJ(lss)->params->iterations_bound, sizeof (double)), LSSERR_OUT_OF_MEMORY);
      break;
    case LSS_RELATIVE_ACCURACY:
      ASSERT (v > 0, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->relative_accuracy = v;
      break;
    case LSS_ABSOLUTE_ACCURACY:
      ASSERT (v > 0, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->absolute_accuracy = v;
      break;
    case LSS_PRECONDITIONER:
      ASSERT (v == 0 || v == 1 || v == 2 || v == 3, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->preconditioner = (int)v;
      LSSOBJ(lss)->modified = 1; /* need to rebuild preconditioner */
      break;
    case LSS_DECIMATION:
      ASSERT (v >= 2, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->decimation = (int)v;
      LSSOBJ(lss)->modified = 1; /* need to rebuild preconditioner */
      break;
    case LSS_CUTOFF:
      ASSERT (v >= 2, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->cutoff = (int)v;
      LSSOBJ(lss)->modified = 1; /* need to rebuild preconditioner */
      break;
    case LSS_RESTART:
      ASSERT (v >= 1, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->restart = (int)v;
      break;
    case LSS_COARSE_RESTART:
      ASSERT (v >= 1, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->coarse_restart = (int)v;
      break;
    case LSS_COARSE_ITERATIONS_BOUND:
      ASSERT (v >= 1, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->coarse_iterations_bound = (int)v;
      break;
    case LSS_COARSE_RELATIVE_ACCURACY:
      ASSERT (v >= 1, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->coarse_relative_accuracy = v;
      break;
    case LSS_COARSE_ABSOLUTE_ACCURACY:
      ASSERT (v >= 1, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->coarse_absolute_accuracy = v;
      break;
    case LSS_PRESMOOTHER:
      ASSERT (v == 1 || v == 2 || v == 3, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->presmoother = (int)v;
      break;
    case LSS_POSTSMOOTHER:
      ASSERT (v == 1 || v == 2 || v == 3, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->postsmoother = (int)v;
      break;
    case LSS_PRESMOOTHING_STEPS:
      ASSERT (v >= 0, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->presmoothing_steps = (int)v;
      break;
    case LSS_POSTSMOOTHING_STEPS:
      ASSERT (v >= 0, LSSERR_INVALID_ARGUMENT);
      LSSOBJ(lss)->params->postsmoothing_steps = (int)v;
      break;
    /* read-only */
    default:
      ASSERT (0, LSSERR_INVALID_ARGUMENT);
      break;
  }

  return LSSERR_NONE;
}

/* lss: returned by LSS_Create
 * p: parameter
 * return: value;
 * get parameter */
double LSS_Get (void *lss, LSSPAR p)
{
  if (!lss) return 0.0;

  switch (p)
  {
    case LSS_ITERATIONS_BOUND:
      return LSSOBJ(lss)->params->iterations_bound;
    case LSS_RELATIVE_ACCURACY:
      return LSSOBJ(lss)->params->relative_accuracy;
    case LSS_ABSOLUTE_ACCURACY:
      return LSSOBJ(lss)->params->absolute_accuracy;
    case LSS_PRECONDITIONER:
      return LSSOBJ(lss)->params->preconditioner;
    case LSS_DECIMATION:
      return LSSOBJ(lss)->params->decimation;
    case LSS_CUTOFF:
      return LSSOBJ(lss)->params->cutoff;
    case LSS_RESTART:
      return LSSOBJ(lss)->params->restart;
    case LSS_COARSE_RESTART:
      return LSSOBJ(lss)->params->coarse_restart;
    case LSS_COARSE_ITERATIONS_BOUND:
      return LSSOBJ(lss)->params->coarse_iterations_bound;
    case LSS_COARSE_RELATIVE_ACCURACY:
      return LSSOBJ(lss)->params->coarse_relative_accuracy;
    case LSS_COARSE_ABSOLUTE_ACCURACY:
      return LSSOBJ(lss)->params->coarse_absolute_accuracy;
    case LSS_PRESMOOTHER:
      return LSSOBJ(lss)->params->presmoother;
    case LSS_POSTSMOOTHER:
      return LSSOBJ(lss)->params->postsmoother;
    case LSS_PRESMOOTHING_STEPS:
      return LSSOBJ(lss)->params->presmoothing_steps;
    case LSS_POSTSMOOTHING_STEPS:
      return LSSOBJ(lss)->params->postsmoothing_steps;
    case LSS_ITERATIONS:
      return LSSOBJ(lss)->params->iterations;
    case LSS_RELATIVE_ERROR:
      return LSSOBJ(lss)->params->relative_error [LSSOBJ(lss)->params->iterations - 1];
    case LSS_ABSOLUTE_ERROR:
      return LSSOBJ(lss)->params->absolute_error [LSSOBJ(lss)->params->iterations - 1];
    case LSS_LEVELS:
      return LSSOBJ(lss)->params->levels;
    case LSS_OPERATOR_COMPLEXITY:
      return LSSOBJ(lss)->params->operator_complexity;
    case LSS_GRID_COMPLEXITY:
      return LSSOBJ(lss)->params->grid_complexity;
    default: break;
  }

  return 0.0;
}

/* lss: returned by LSS_Create 
 * p: parameter
 * return: vector;
 * get vector parameter */
double* LSS_Getv (void *lss, LSSPAR p)
{
  if (!lss) return NULL;

  switch (p)
  {
    case LSS_RELATIVE_ERROR:
      return LSSOBJ(lss)->params->relative_error;
    case LSS_ABSOLUTE_ERROR:
      return LSSOBJ(lss)->params->absolute_error;
    case LSS_LEVEL_DIMENSIONS:
      return LSSOBJ(lss)->params->level_dimensions;
    case LSS_LEVEL_NONZEROS:
      return LSSOBJ(lss)->params->level_nonzeros;
    default: return NULL;
  }

  return NULL;
}

/* lss: returned by LSS_Create  
 * a: system A matrix values (same (p, i) as in LSS_Create is assumed)
 * x: solution vector
 * b: right hand side vector
 * solve linear system A x = b, where A = (a, p, i) */
LSSERR LSS_Solve (void *lss, double *a, double *x, double *b)
{
  if (!lss) return LSSERR_OUT_OF_MEMORY;

  ONERROR (LSSOBJ(lss)->error) return LSSOBJ(lss)->error; /* an error has occured */

  if (modified (lss)) /* has system structure been modified? */
  {
    destroy (lss);  /* destroy current data */
    CALL (create, lss, a, x, b); /* create new data */
  }

  CALL (update, lss, a, x, b); /* update current data */
  CALL (solve, lss, x); /* solve and output */

  return LSSERR_NONE; /* successful termination */
}

/* lss: returned by LSS_Create;
 * release solver memory */
void LSS_Destroy (void *lss)
{
  if (!lss) return;

  destroy (lss);
  free (LSSOBJ(lss)->params);
  free (lss);
}
