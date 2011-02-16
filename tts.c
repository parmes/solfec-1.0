/*
 * tts.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * test solver
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

#include <complex.h>
#include <stdlib.h>
#include <float.h>

#include "tts.h"
#include "dom.h"
#include "fem.h"
#include "alg.h"
#include "bla.h"
#include "lap.h"
#include "err.h"
#include "vic.h"
#include "mrf.h"
#include "ext/krylov/krylov.h"

typedef struct con_data CON_DATA;
typedef struct private PRIVATE;
typedef struct vector VECTOR;

struct con_data
{
  double X [9], /* U-linearization */
	 Y [9], /* R-linearization */
	 T [9]; /* diagonal preconditioner */

  CON *con; /* constraint */
};

struct private
{
  CON_DATA *start, *end; /* starte and end of data */

  VECTOR *b, /* right hand side */
	 *x; /* reactions incement */

  double epsilon, /* smoothing epsilon */
	 delta, /* egularization */
	 bnorm;

  int iters; /* linear solver iterations */

  int matvec; /* counter */

  DOM *dom; /* domain */
};

struct vector
{
  double *x;

  int n;
};

/* allocate vector */
static VECTOR* newvector (int n)
{
  VECTOR *v;

  ERRMEM (v = malloc (sizeof (VECTOR)));
  ERRMEM (v->x = MEM_CALLOC (n * sizeof (double)));
  v->n = n;

  return v;
}

/* GMRES interface start */
static char* CAlloc (size_t count, size_t elt_size)
{
  char *ptr;

  ERRMEM (ptr = MEM_CALLOC (count * elt_size));
  return ptr;
}

static int Free (char *ptr)
{
  free (ptr);
  return 0;
}

static int CommInfo (PRIVATE *A, int *my_id, int *num_procs)
{
#if MPI
  *num_procs = A->dom->ncpu;
  *my_id = A->dom->rank;
#else
  *num_procs = 1;
  *my_id = 0;
#endif
  return 0;
}

static void* CreateVector (VECTOR *a)
{
  VECTOR *v;

  ERRMEM (v = malloc (sizeof (VECTOR)));
  ERRMEM (v->x = MEM_CALLOC (a->n * sizeof (double)));
  v->n = a->n;

  return v;
}

static void* CreateVectorArray (int size, VECTOR *a)
{
  VECTOR **v;
  int i;

  ERRMEM (v = malloc (size * sizeof (VECTOR*)));
  for (i = 0; i < size; i ++)
  {
    v[i] = CreateVector (a);
  }

  return v;
}

static int DestroyVector (VECTOR *a)
{
  free (a->x);
  free (a);

  return 0;
}

static double InnerProd (VECTOR *a, VECTOR *b)
{
  double dot = 0.0, *x, *y, *z;

  for (x = a->x, z = x + a->n, y = b->x; x < z; x ++, y ++)
  {
    dot += (*x) * (*y);
  }

#if MPI
  double val = dot;
  MPI_Allreduce (&val, &dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  return dot;
}

static int CopyVector (VECTOR *a, VECTOR *b)
{
  double *x, *y, *z;

  for (x = a->x, z = x + a->n, y = b->x; x < z; x ++, y ++)
  {
    (*y) = (*x);
  }

  return 0;
}

static int ClearVector (VECTOR *a)
{
  double *x, *z;

  for (x = a->x, z = x + a->n; x < z; x ++)
  {
    (*x) = 0.0;
  }

  return 0;
}

static int ScaleVector (double alpha, VECTOR *a)
{
  double *x, *z;

  for (x = a->x, z = x + a->n; x < z; x ++)
  {
    (*x) *= alpha;
  }

  return 0;
}

static int  Axpy (double alpha, VECTOR *a, VECTOR *b)
{
  double *x, *y, *z;

  for (x = a->x, z = x + a->n, y = b->x; x < z; x ++, y ++)
  {
    (*y) += alpha * (*x);
  }

  return 0;
}

static void *MatvecCreate (void *A, void *x)
{
  return NULL;
}

static int Matvec (void *matvec_data, double alpha, PRIVATE *A, VECTOR *x, double beta, VECTOR *y)
{
  double Z [3], U [3], *X, *Y, *W, *R, *Q, *r, *v;
  CON *con, *adj;
  CON_DATA *dat;
  DIAB *dia;
  OFFB *blk;

#if MPI
  {
    COMDATA *send, *recv, *ptr;
    DOM *dom = A->dom;
    int nsend, nrecv;
    int i, *j, *k;
    SET *item;

    for (dat = A->start, nsend = 0; dat != A->end; dat ++) 
    {
      con = dat->con;
      nsend += SET_Size (con->ext);
    }
    ERRMEM (send = MEM_CALLOC (sizeof (COMDATA [nsend])));
    for (dat = A->start, ptr = send; dat != A->end; dat ++)
    {
      con = dat->con;
      for (item = SET_First (con->ext); item; item = SET_Next (item), ptr ++)
      {
	ptr->ints = 1;
	ptr->d = &x->x [con->num];
	ptr->doubles = 3;
	ptr->i = (int*) &con->id;
	ptr->rank = (int) (long) item->data;
      }
    }
    COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);
    for (i = 0; i < nrecv; i ++)
    {
      ptr = &recv [i];
      for (j = ptr->i, k = j + ptr->ints, R = ptr->d; j < k; j ++, R += 3)
      {
	ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*j), NULL), "Invalid con id: %d", *j);
	COPY (R, con->R);
      }
    }
    free (send);
    free (recv);
  }
#endif

  for (dat = A->start, R = r = x->x, Q = y->x; dat != A->end; dat ++, R += 3, Q += 3)
  {
    con = dat->con;
    dia = con->dia;
    v = &r [con->num];
    W = dia->W;
    NVMUL (W, v, U);
    for (blk = dia->adj; blk; blk = blk->n)
    {
      adj = blk->dia->con;
      if (adj->num >= 0)
      {
	v = &r [adj->num];
	W = blk->W;
	NVADDMUL (U, W, v, U);
      }
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      v = CON(blk->dia)->R;
      W = blk->W;
      NVADDMUL (U, W, v, U);
    }
#endif

    switch (con->kind)
    {
    case VELODIR:
    case FIXDIR:
    case RIGLNK:
    {
      U [0] = R [0];
      U [1] = R [1];
    }
    break;
    case CONTACT:
    {
      X = dat->X;
      Y = dat->Y;

      NVMUL (X, U, Z);
      NVADDMUL (Z, Y, R, U);  /* U = X dU + Y dR */
    }
    break;
    default:
    break;
    }

    SCALE (Q, beta);
    ADDMUL (Q, alpha, U, Q);
  }

  Axpy (alpha * A->delta, x, y);

  A->matvec ++;

  return 0;
}

static int MatvecDestroy (void *matvec_data)
{
  return 0;
}

static int PrecondSetup (void *vdata, void *A, void *b, void *x)
{
  return 0;
}

static int Precond0 (void *vdata, PRIVATE *A, VECTOR *b, VECTOR *x)
{
  double *T, *Q, *R;
  CON_DATA *dat;

  for (dat = A->start, R = x->x, Q = b->x; dat != A->end; dat ++, R += 3, Q += 3)
  {
    T = dat->T;
    NVMUL (T, Q, R);
  }

  return 0;
}

#if 0
static int Precond1 (void *vdata, PRIVATE *A, VECTOR *b, VECTOR *x)
{
  hypre_FlexGMRESFunctions *gmres_functions;
  void *gmres_vdata;
  int ret;

  gmres_functions = hypre_FlexGMRESFunctionsCreate (CAlloc, Free, (int (*) (void*,int*,int*)) CommInfo,
    (void* (*) (void*))CreateVector, (void* (*) (int, void*))CreateVectorArray, (int (*) (void*))DestroyVector,
    MatvecCreate, (int (*) (void*,double,void*,void*,double,void*))Matvec, MatvecDestroy,
    (double (*) (void*,void*))InnerProd, (int (*) (void*,void*))CopyVector, (int (*) (void*))ClearVector,
    (int (*) (double,void*))ScaleVector, (int (*) (double,void*,void*))Axpy,
    PrecondSetup, (int (*) (void*,void*,void*,void*))Precond0);
  gmres_vdata = hypre_FlexGMRESCreate (gmres_functions);

  hypre_FlexGMRESSetMinIter (gmres_vdata, 3);
  hypre_FlexGMRESSetMaxIter (gmres_vdata, 3);
  hypre_FlexGMRESSetup (gmres_vdata, A, b, x);
  ret = hypre_FlexGMRESSolve (gmres_vdata, A, b, x);
  hypre_FlexGMRESDestroy (gmres_vdata);

  return 0;
}
#endif

#if 0
static int Precond2 (void *vdata, PRIVATE *A, VECTOR *b, VECTOR *x)
{
  int n = A->end - A->start, i, j;
  MEM mapmem, dblmem;
  MAP **col;

  ERRMEM (col = MEM_CALLOC (3*n*sizeof (SET*)));
  MEM_Init (&mapmem, sizeof (MAP), 1204);
  MEM_Init (&dblmem, sizeof (double [9]), 256);

  for (CON_DATA *dat = A->start; dat != A->end; dat ++)
  {
    double *X = dat->X, *Y = dat->Y;
    CON *con = dat->con;
    DIAB *dia = con->dia;
    double *W = dia->W;
    i = con->num, j = i;
    double *T = MEM_Alloc (&dblmem);
    int ii, jj;
    NNMUL (X, W, T);
    NNADD (T, Y, T);
    T [0] += A->delta;
    T [4] += A->delta;
    T [8] += A->delta;
    for (ii = 0; ii < 3; ii ++)
      for (jj = 0; jj < 3; jj ++)
        MAP_Insert (&mapmem, &col [j+jj], (void*) (long) (i+ii), &T[3*jj+ii], NULL);

#if 0
    for (OFFB *blk = dia->adj; blk; blk = blk->n)
    {
      j = blk->dia->con->num;
      if (ABS (i-j) == 3)
      {
	double XW [9];
	W = blk->W;
	if (!(T = MAP_Find (col [j], (void*) (long) i, NULL))) T = MEM_Alloc (&dblmem);
	NNMUL (X, W, XW);
	NNADD (T, XW, T);
	for (ii = 0; ii < 3; ii ++)
	  for (jj = 0; jj < 3; jj ++)
	    MAP_Insert (&mapmem, &col [j+jj], (void*) (long) (i+ii), &T[3*jj+ii], NULL);
      }
    }
#endif
  }

  MX *P;
  ERRMEM (P = MEM_CALLOC (sizeof (MX)));
  P->nz = -1;
  P->kind = MXCSC;
  P->m = P->n = 3*n;
  ERRMEM (P->p = MEM_CALLOC (sizeof (int [P->m+1])));
  for (j = 1; j < P->m+1; j ++)
  {
    P->p [j] = P->p[j-1] + MAP_Size (col [j-1]);
  }
  P->nzmax = P->p [P->m];
  ERRMEM (P->i = MEM_CALLOC (sizeof (int [P->nzmax])));
  ERRMEM (P->x = MEM_CALLOC (sizeof (double [P->nzmax])));
  for (j = 0, i = 0; j < P->m; j ++)
  {
    MAP *item;
    for (item = MAP_First (col [j]); item; item = MAP_Next (item), i ++)
    {
      P->i [i] = (int) (long) item->key;
      P->x [i] = *((double*) item->data);
    }
  }

  MX_Inverse (P, P);
  MX_Matvec (1.0, P, b->x, 0.0, x->x);
  MX_Destroy (P);
  MEM_Release (&mapmem);
  MEM_Release (&dblmem);
  free (col);

  return 0;
}
#endif

/* GMRES interface end */

/* create private data */
static PRIVATE *create_data (TEST *ts, LOCDYN *ldy)
{
  CON_DATA *dat;
  PRIVATE *A;
  CON *con;
  DOM *dom;
  int n;

  dom = ldy->dom;
  ERRMEM (A = MEM_CALLOC (sizeof (PRIVATE)));
  ERRMEM (A->start = MEM_CALLOC (sizeof (CON_DATA [dom->ncon])));
  A->end = A->start;
  A->dom = dom;

  for (con = dom->con, dat = A->start, n = 0; con; con = con->next)
  {
    if (dom->dynamic && con->kind == CONTACT && con->gap > 0.0)
    {
      SET (con->R, 0.0);
      con->num = -1;
      continue;
    }
    dat->con = con;
    con->num = n;
    A->end ++;
    dat ++;
    n += 3;
  }

  A->b = newvector (n);
  A->x = newvector (n);

  /* U = W R + B */
  for (dat = A->start; dat != A->end; dat ++)
  {
    CON *con = dat->con;
    DIAB *dia = con->dia;
    double *W = dia->W, *U = con->U;
    NVADDMUL (dia->B, W, con->R, U);
    for (OFFB *blk = dia->adj; blk; blk = blk->n)
    {
      double *R = blk->dia->R;
      W = blk->W;
      NVADDMUL (U, W, R, U);
    }
#if MPI
    for (OFFB *blk = dia->adjext; blk; blk = blk->n)
    {
      double *R = CON(blk->dia)->R;
      W = blk->W;
      NVADDMUL (U, W, R, U);
    }
#endif
  }

  return A;
}

/* update linear system */
static void update_system (PRIVATE *A)
{
  double *b = A->b->x, epsilon = A->epsilon;
  DOM *dom = A->dom;
  short dynamic = dom->dynamic;
  CON_DATA *dat;

  for (dat = A->start; dat != A->end; dat ++, b += 3)
  {
    CON *con = dat->con;
    DIAB *dia = con->dia;
    double *U = con->U,
	   *V = con->V,
	   *R = con->R,
	   *W = dia->W,
	   *T = dat->T;

    switch (con->kind)
    {
    case FIXPNT:
    case GLUE:
    {
      if (dynamic)
      {
	b [0] = -V[0]-U[0];
	b [1] = -V[1]-U[1];
	b [2] = -V[2]-U[2];
      }
      else
      {
	b [0] = -U[0];
	b [1] = -U[1];
	b [2] = -U[2];
      }

      NNCOPY (W, T);
    }
    break;
    case FIXDIR:
    {
      b [0] = -R[0];
      b [1] = -R[1];
      if (dynamic) b [2] = -V[2]-U[2];
      else b [2] = -U[2];

      T [1] = T [3] = T [6] = T [7] = 0.0;
      T [0] = T [4] = 1.0;
      T [2] = W [2];
      T [5] = W [5];
      T [8] = W [8];
    }
    break;
    case VELODIR:
    {
      b [0] = -R[0];
      b [1] = -R[1];
      b [2] = VELODIR(dat->con->Z)-U[2];

      T [1] = T [3] = T [6] = T [7] = 0.0;
      T [0] = T [4] = 1.0;
      T [2] = W [2];
      T [5] = W [5];
      T [8] = W [8];
    }
    break;
    case RIGLNK:
    {
      double h = dom->step * (dynamic ? 0.5 : 1.0),
             d = RIGLNK_LEN (dat->con->Z),
	     delta;

      b [0] = -R[0];
      b [1] = -R[1];
      delta = d*d - h*h*DOT2(U,U);
      if (delta >= 0.0) b [2] = (sqrt (delta) - d)/h - U[2];
      else b[2] = -U[2];

      T [1] = T [3] = T [6] = T [7] = 0.0;
      T [0] = T [4] = 1.0;
      T [2] = W [2];
      T [5] = W [5];
      T [8] = W [8];
    }
    break;
    case CONTACT:
    {
      double *X = dat->X, *Y = dat->Y;

      VIC_Linearize (dat->con, U, R, -1, epsilon, b, X, Y);
      SCALE (b, -1.0);

      NNMUL (X, W, T);
      NNADD (T, Y, T);
    }
    break;
    }

    T [0] += A->delta;
    T [4] += A->delta;
    T [8] += A->delta;

    MX_DENSE_PTR (P, 3, 3, T);
    MX_Inverse (&P, &P);
  }

  A->bnorm = sqrt (InnerProd (A->b, A->b));
}

/* solve linear system */
static void linear_solve (PRIVATE *A, double abstol, int maxiter)
{
  hypre_FlexGMRESFunctions *gmres_functions;
  void *gmres_vdata;
  int ret;

  gmres_functions = hypre_FlexGMRESFunctionsCreate (CAlloc, Free, (int (*) (void*,int*,int*)) CommInfo,
    (void* (*) (void*))CreateVector, (void* (*) (int, void*))CreateVectorArray, (int (*) (void*))DestroyVector,
    MatvecCreate, (int (*) (void*,double,void*,void*,double,void*))Matvec, MatvecDestroy,
    (double (*) (void*,void*))InnerProd, (int (*) (void*,void*))CopyVector, (int (*) (void*))ClearVector,
    (int (*) (double,void*))ScaleVector, (int (*) (double,void*,void*))Axpy,
    PrecondSetup, (int (*) (void*,void*,void*,void*))Precond0);
  gmres_vdata = hypre_FlexGMRESCreate (gmres_functions);

  hypre_error_flag = 0;
  hypre_FlexGMRESSetTol (gmres_vdata, 0.0);
  hypre_FlexGMRESSetMinIter (gmres_vdata, 1);
  hypre_FlexGMRESSetMaxIter (gmres_vdata, maxiter);
  hypre_FlexGMRESSetAbsoluteTol (gmres_vdata, abstol);
  hypre_FlexGMRESSetup (gmres_vdata, A, A->b, A->x);
  ret = hypre_FlexGMRESSolve (gmres_vdata, A, A->b, A->x);
  hypre_FlexGMRESGetNumIterations (gmres_vdata , &A->iters);
  hypre_FlexGMRESDestroy (gmres_vdata);
}

/* update solution */
static void update_solution (PRIVATE *A)
{
  double *R, *dR;
  CON_DATA *dat;
  CON *con;

  for (dat = A->start, dR = A->x->x; dat != A->end; dat ++, dR += 3)
  {
    con = dat->con;
    R = con->R;
    ACC (dR, R);

    if (con->kind == CONTACT)
    {
      double c = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area;
      VIC_Project (con->mat.base->friction, c, R, R);
    }
  }

  /* U = W R + B */
#if MPI
  COMDATA *send, *recv, *ptr;
  DOM *dom = A->dom;
  int nsend, nrecv;
  int i, *j, *k;
  SET *item;

  for (dat = A->start, nsend = 0; dat != A->end; dat ++) 
  {
    con = dat->con;
    nsend += SET_Size (con->ext);
  }
  ERRMEM (send = MEM_CALLOC (sizeof (COMDATA [nsend])));
  for (dat = A->start, ptr = send; dat != A->end; dat ++)
  {
    con = dat->con;
    for (item = SET_First (con->ext); item; item = SET_Next (item), ptr ++)
    {
      ptr->ints = 1;
      ptr->d = con->R;
      ptr->doubles = 3;
      ptr->i = (int*) &con->id;
      ptr->rank = (int) (long) item->data;
    }
  }
  COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);
  for (i = 0; i < nrecv; i ++)
  {
    ptr = &recv [i];
    for (j = ptr->i, k = j + ptr->ints, R = ptr->d; j < k; j ++, R += 3)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*j), NULL), "Invalid con id: %d", *j);
      COPY (R, con->R);
    }
  }
  free (send);
  free (recv);
#endif

  for (dat = A->start; dat != A->end; dat ++)
  {
    CON *con = dat->con;
    DIAB *dia = con->dia;
    double *W = dia->W, *U = con->U;
    NVADDMUL (dia->B, W, con->R, U);
    for (OFFB *blk = dia->adj; blk; blk = blk->n)
    {
      double *R = blk->dia->R;
      W = blk->W;
      NVADDMUL (U, W, R, U);
    }
#if MPI
     for (OFFB *blk = dia->adjext; blk; blk = blk->n)
    {
      double *R = CON(blk->dia)->R;
      W = blk->W;
      NVADDMUL (U, W, R, U);
    }
#endif
  }
}

/* destroy private data */
static void destroy_data (PRIVATE *A)
{
  DestroyVector (A->b);
  DestroyVector (A->x);
  free (A->start);
  free (A);
}

/* return average abs (Wii) */
static double average_abs_Wii (LOCDYN *ldy)
{
  double Wii [2] = {0, 0};

  for (DIAB *dia = ldy->dia; dia; dia = dia->n)
  {
    double *W = dia->W;
    for (int i = 0; i < 9; i ++) Wii[0] += fabs (W[i]);
    Wii [1] += 9.0;
  }

#if MPI
  double Wjj [2] = {Wii [0], Wii [1]};
  MPI_Allreduce (Wjj, Wii, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  return Wii [0] / Wii [1];
}

/* create solver */
TEST* TEST_Create (double meritval, int maxiter)
{
  TEST *ts;

  ERRMEM (ts = MEM_CALLOC (sizeof (TEST)));
  ts->meritval = meritval;
  ts->maxiter = maxiter;
  ts->linmaxiter = maxiter * 10;

  return ts;
}

/* run solver */
void TEST_Solve (TEST *ts, LOCDYN *ldy)
{
  char fmt [512];
  double *merit;
  PRIVATE *A;
  DOM *dom;

  sprintf (fmt, "TEST_SOLVER: (LIN its: %%%dd) iteration: %%%dd merit: %%.2e\n",
           (int)log10 (ts->linmaxiter) + 1, (int)log10 (ts->maxiter) + 1);

  ERRMEM (ts->merhist = realloc (ts->merhist, ts->maxiter * sizeof (double)));
  dom = ldy->dom;
  A = create_data (ts, ldy);
  merit = &dom->merit;
  ts->iters = 0;
  A->delta = 0.1 * average_abs_Wii (ldy);
  A->epsilon = 1E-8;

  do
  {
    update_system (A);

    linear_solve (A, 0.25 * A->bnorm,  ts->linmaxiter);

    update_solution (A);

    *merit = MERIT_Function (ldy, 0);

    ts->merhist [ts->iters] = *merit;

#if MPI
    if (dom->rank == 0)
#endif
    if (dom->verbose) printf (fmt, A->iters, ts->iters, *merit);

  } while (++ ts->iters < ts->maxiter && *merit > ts->meritval);

#if MPI
    if (dom->rank == 0)
#endif
    if (dom->verbose) printf ("TEST_SOLVER: total MATVECs: %d\n", A->matvec);

  destroy_data (A);
}

/* write labeled state values */
void TEST_Write_State (TEST *ts, PBF *bf)
{
  PBF_Label (bf, "TSITERS");
  PBF_Int (bf, &ts->iters, 1);
}

/* destroy solver */
void TEST_Destroy (TEST *ts)
{
  free (ts->merhist);
  free (ts);
}
