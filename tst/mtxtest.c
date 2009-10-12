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

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "mtx.h"
#include "alg.h"
#include "err.h"

#define EPSILON 1E-9

static void skip_line (FILE *file)
{
  char c;

  while ((c = fgetc (file)) != '\n' && c != EOF);
}

static void skip_blanks (FILE *file)
{

  char c;

  while ((c = fgetc (file)) != EOF)
  {
    if (isblank (c) || c == '\n') continue;
    else if (c == '#') skip_line (file);
    else break;
  }
 
  if (c != EOF) ungetc (c, file);
}

static void dense_to_sparse (double *a, int m, int n, double **x, int **p, int **i)
{
  int j, k, l, o, nzmax;

  l = m * n;
  nzmax = 0;

  for (k = 0; k < l; k ++) if (a [k] != 0.0) nzmax ++;

  ERRMEM ((*x) = malloc (sizeof (double [nzmax])));
  ERRMEM ((*p) = malloc (sizeof (int [n+1])));
  ERRMEM ((*i) = malloc (sizeof (int [nzmax])));

  for (o = l = j = 0; j < n; j ++)
  {
    (*p) [j] = o;

    for (k = 0; k < m; k ++, l ++)
    {
      if (a [l] != 0.0)
      {
	(*x) [o] = a [l];
	(*i) [o] = k;
	o ++;
      }
    }
  }

  (*p) [j] = o;
}

typedef enum {INPUT, OUTPUT} OP;

static MX* read_matrix (OP op, FILE *file, int m, int n, short kind)
{
  double *a, *x, *y;
  int *p, *i, j, k;
  char name;
  MX *mat;

  skip_blanks (file);

  if (op == INPUT) fscanf (file, "%c(%d,%d)", &name, &m, &n);

  ERRMEM (a = malloc (sizeof (double [m*n])));

  for (k = 0; k < m; k ++)
  {
    for (j = 0; j < n; j ++)
    {
      fscanf (file, "%lf", &a [m*j+k]);
    }
  }

  switch (kind)
  {
  case MXDENSE:
  {
    mat = MX_Create (MXDENSE, m, n, NULL, NULL);

    for (k = 0, x = mat->x; k < mat->nzmax; k ++) x [k] = a [k];
  }
  break;
  case MXCSC:
  {
    dense_to_sparse (a, m, n, &y, &p, &i);

    mat = MX_Create (MXCSC, m, n, p, i);

    for (k = 0, x = mat->x; k < mat->nzmax; k ++) x [k] = y [k];

    free (y);
    free (p);
    free (i);
  }
  break;
  default:
    ASSERT (0, ERR_NOT_IMPLEMENTED);
  }

  if (op == INPUT) printf ("Red input %s matrix %c(%d,%d).\n", kind == MXDENSE ? "dense" : "sparse" , name, m, n);

  free (a);

  return mat;
}

static int test_run_0 (FILE *file, MX *A, MX *B, MX *C, MX *D, MX *E, double alpha, double beta, short kind)
{
  MX *X, *Y, *Z, *invA, *invB;

  printf ("TEST: alpha * A + beta * B ... ");
  X = MX_Add (alpha, A, beta, B, NULL);
  Y = read_matrix (OUTPUT, file, A->m, A->n, kind);
  Z = MX_Add (1.0, X, -1.0, Y, NULL);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: alpha * A' + beta * B ... ");
  MX_Add (alpha, MX_Tran (A), beta, B, X);
  Y = read_matrix (OUTPUT, file, A->m, A->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: alpha * A + beta * B' ... ");
  MX_Add (alpha, A, beta, MX_Tran (B), X);
  Y = read_matrix (OUTPUT, file, A->m, A->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: alpha * A' + beta * B' ... ");
  MX_Add (alpha, MX_Tran (A), beta, MX_Tran (B), X);
  Y = read_matrix (OUTPUT, file, A->m, A->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: A * B ... ");
  MX_Matmat (1.0, A, B, 0.0, X);
  Y = read_matrix (OUTPUT, file, A->m, B->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: A' * B ... ");
  MX_Matmat (1.0, MX_Tran (A), B, 0.0, X);
  Y = read_matrix (OUTPUT, file, A->n, B->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: A * B' ... ");
  MX_Matmat (1.0, A, MX_Tran (B), 0.0, X);
  Y = read_matrix (OUTPUT, file, A->m, B->m, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: A' * B' ... ");
  MX_Matmat (1.0, MX_Tran (A), MX_Tran (B), 0.0, X);
  Y = read_matrix (OUTPUT, file, A->n, B->m, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  invA = MX_Inverse (A, NULL);
  invB = MX_Inverse (B, NULL);

  MX_Destroy (X);
  MX_Destroy (Z);

  printf ("TEST: inv (A) * B ... ");
  X = MX_Matmat (1.0, invA, B, 0.0, NULL);
  Y = read_matrix (OUTPUT, file, invA->m, B->n, kind);
  Z = MX_Add (1.0, X, -1.0, Y, NULL);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: inv (A)' * B ... ");
  MX_Matmat (1.0, MX_Tran (invA), B, 0.0, X);
  Y = read_matrix (OUTPUT, file, invA->n, B->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: A * inv (B) ... ");
  MX_Matmat (1.0, A, invB, 0.0, X);
  Y = read_matrix (OUTPUT, file, A->m, invB->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: A * inv (B)' ... ");
  MX_Matmat (1.0, A, MX_Tran (invB), 0.0, X);
  Y = read_matrix (OUTPUT, file, A->m, invB->m, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: alpha * A * B + beta * C ... ");
  MX_Copy (C, X);
  MX_Matmat (alpha, A, B, beta, X);
  Y = read_matrix (OUTPUT, file, A->m, B->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: alpha * A * B + beta * C' ... ");
  MX_Copy (MX_Tran (C), X);
  MX_Matmat (alpha, A, B, beta, X);
  Y = read_matrix (OUTPUT, file, A->m, B->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: alpha * inv (A) * B + beta * C ... ");
  MX_Copy (C, X);
  MX_Matmat (alpha, invA, B, beta, X);
  Y = read_matrix (OUTPUT, file, invA->m, B->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: alpha * A * inv(B) + beta * C' ... ");
  MX_Copy (MX_Tran (C), X);
  MX_Matmat (alpha, A, invB, beta, X);
  Y = read_matrix (OUTPUT, file, A->m, invB->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: E * D ... ");
  MX_Matmat (1.0, E, D, 0.0, X);
  Y = read_matrix (OUTPUT, file, E->m, D->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: E * A * D ... ");
  MX_Trimat (E, A, D, X);
  Y = read_matrix (OUTPUT, file, E->m, D->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: D' * inv(B) * D ... ");
  MX_Trimat (MX_Tran (D), invB, D, X);
  Y = read_matrix (OUTPUT, file, D->n, D->n, kind);
  MX_Add (1.0, X, -1.0, Y, Z);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  printf ("TEST: E * inv(A) * E' ... ");
  X = MX_Trimat (E, invA, MX_Tran (E), NULL);
  Y = read_matrix (OUTPUT, file, E->m, E->m, kind);
  Z = MX_Add (1.0, X, -1.0, Y, NULL);
  MX_Destroy (Y);
  if (MX_Norm (Z) < EPSILON) printf ("OK\n");
  else { printf ("FAILED\n"); return 0; }

  MX_Destroy (X);
  MX_Destroy (Z);
  MX_Destroy (invA);
  MX_Destroy (invB);

  return 1;
}

static int test_set_0 (void)
{
  MX *A, *B, *C, *D, *E;
  double alpha, beta;
  FILE *file;

  ASSERT (file = fopen ("inp/mtx.0", "r"), ERR_FILE_OPEN);

  A = read_matrix (INPUT, file, 0, 0, MXDENSE);
  B = read_matrix (INPUT, file, 0, 0, MXDENSE);
  C = read_matrix (INPUT, file, 0, 0, MXDENSE);
  D = read_matrix (INPUT, file, 0, 0, MXDENSE);
  E = read_matrix (INPUT, file, 0, 0, MXDENSE);

  alpha = 9.81;
  beta = 0.123;

  if (!test_run_0 (file, A, B, C, D, E, alpha, beta, MXDENSE)) return 0;

  MX_Destroy (A);
  MX_Destroy (B);
  MX_Destroy (C);
  MX_Destroy (D);
  MX_Destroy (E);

  fseek (file, 0, SEEK_SET);

  A = read_matrix (INPUT, file, 0, 0, MXCSC);
  B = read_matrix (INPUT, file, 0, 0, MXCSC);
  C = read_matrix (INPUT, file, 0, 0, MXCSC);
  D = read_matrix (INPUT, file, 0, 0, MXCSC);
  E = read_matrix (INPUT, file, 0, 0, MXCSC);

  if (!test_run_0 (file, A, B, C, D, E, alpha, beta, MXCSC)) return 0;

  MX_Destroy (A);
  MX_Destroy (B);
  MX_Destroy (C);
  MX_Destroy (D);
  MX_Destroy (E);

  return 1;
}

static int test_all (void)
{
  if (!test_set_0 ()) return 0;

  return 1;
}

int main (int argc, char **argv)
{
  return test_all ();
}
