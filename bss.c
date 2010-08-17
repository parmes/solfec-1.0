/*
 * bss.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * body space solver
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

#include "bss.h"
#include "dom.h"
#include "alg.h"
#include "bla.h"
#include "lap.h"
#include "err.h"
#include "ext/krylov/krylov.h"

#if MPI
#include "com.h"
#include "pck.h"
#endif

#define DIFF_FACTOR             1E-10  /* TODO: test sensitivity */
#define SMOOTHING		1      /* TODO: -||- */
#define SMOOTHING_EPSILON	1E-10  /* TODO: -||- */
#define DISABLE_NORM_SMOOTHING	1      /* TODO: -||- */
#define ABSTOL			1E-15  /* TODO: -||- */ 
#define CON_RHS_SET(con, value) (con)->dia->mprod = (MX*)(value)
#define CON_RHS_GET(con) ((double*)(con)->dia->mprod)
#define CON_X_SET(con, value) (con)->dia->sprod = (MX*)(value)
#define CON_X_GET(con) ((double*)(con)->dia->sprod)
#define CONTACT_X(con) (con)->dia->W
#define CONTACT_Y(con) (con)->dia->A

typedef struct bss_data BSS_DATA;
typedef struct vector VECTOR;

struct bss_data
{
  DOM *dom;

  int nprimal, /* SUM [bod in dom->dom] { bod->dofs } */
      ndual; /* SUM [con in dom->con] { size (con->R) } */

  VECTOR *b, /* right hand side */
         *x; /* nprimal, ndual unknowns */

  double *r, /* reaction workspace of size MAX [bod in dom->bod] { bod->dofs } */
          resnorm; /* linear (dual) residual norm */

  int iters; /* number of linear solver iterations */

#if MPI
  MEM mem; /* integer pairs memory */
  COMDATA *send, *recv; /* communication buffers */
  int ssend, nsend, nrecv; /* sizes */
  void *pattern; /* and pattern */
#endif
};

struct vector
{
  double *x;

  int n;
};

/* update needed before matrix vector product */
static void* update_constraints_data (BSS_DATA *A, double *x);

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

static int CommInfo (BSS_DATA *A, int *my_id, int *num_procs)
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

static int Matvec (void *matvec_data, double alpha, BSS_DATA *A, VECTOR *x, double beta, VECTOR *y)
{
  double c [3], *U, *R, *X, *Y, *a, *b, *r, step;
  DOM *dom = A->dom;
  BODY *bod;
  CON *con;
  int n;

  step = dom->step;

  update_constraints_data (A, x->x);

  ScaleVector (beta, y);

  for (bod = dom->bod, r = A->r, a = x->x, b = y->x; bod; a += n, b += n, bod = bod->next)
  {
    n = bod->dofs;
    BODY_Matvec (alpha, bod, a, 1.0, b);
    BODY_Reac (bod, r);
    blas_daxpy (n, -step*alpha, r, 1, b, 1);
  }

  for (con = dom->con, r = A->b->x; con; con = con->next)
  {
    n = CON_RHS_GET (con) - r;
    b = &y->x [n];
    U = con->U;
    R = con->R;

    switch (con->kind)
    {
    case FIXPNT:
    case GLUE:
      ADDMUL (b, alpha, U, b);
      break;
    case FIXDIR:
    case VELODIR:
      b [0] += alpha * U[2];
      break;
    case RIGLNK:
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
      break;
    case CONTACT:
      X = CONTACT_X (con);
      Y = CONTACT_Y (con);
      NVMUL (X, U, c);
      NVADDMUL (c, Y, R, c);
      ADDMUL (b, alpha, c, b);
      break;
    }
  }

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

static int Precond (void *vdata, BSS_DATA *A, VECTOR *b, VECTOR *x)
{
  double *p, *q;
  BODY *bod;
  int n;

  for  (bod = A->dom->bod, p = b->x, q = x->x; bod; p += n, q += n, bod = bod->next)
  {
    n = bod->dofs;
    BODY_Invvec (1.0, bod, p, 0.0, q);
  }

  blas_dcopy (A->ndual, &b->x[A->nprimal], 1, &x->x[A->nprimal], 1);

  return 0;
}
/* GMRES interface end */

#if MPI
/* allocate a send buffer item */
static COMDATA* send_item (COMDATA **send, int *size, int *count)
{
  if (++ (*count) >= (*size))
  {
    (*size) *= 2;
    ERRMEM ((*send) = realloc ((*send), (*size) * sizeof (COMDATA)));
  }

  return &(*send)[(*count) - 1];
}
#endif

/* update constraints velocities and reactions or
 * create the necessary for that communication pattern */
static void* update_constraints_data (BSS_DATA *A, double *x)
{
  double *R, *U, *B, *D, *u, *velo;
  DOM *dom = A->dom;
  BODY *bod;
  SET *item;
  CON *con;
  int i;
#if MPI
  COMDATA *ptr;
  int *j, *k;
#endif

#if MPI
  if (A->pattern)
  {
#endif
    /* update local velocities of local constraints */
    for (bod = dom->bod, u = x; bod; u += bod->dofs, bod = bod->next)
    {
      velo = bod->velo;
      bod->velo = u;
      for (item = SET_First (bod->con); item; item = SET_Next (item))
      {
	con = item->data;
#if MPI
	if (!(con->state & CON_EXTERNAL)) /* local */
	{
#endif
	  if (bod == con->master) BODY_Local_Velo (bod, mshp(con), mgobj(con), con->mpnt, con->base, NULL, con->U);
	  else BODY_Local_Velo (bod, sshp(con), sgobj(con), con->spnt, con->base, NULL, con->dia->B);
#if MPI
	}
#endif
      }
      bod->velo = velo;
    }

    /* update local constraints reactions */
    for (con = dom->con, u = A->x->x; con; con = con->next)
    {
      i = CON_X_GET (con) - u;
      D = &x [i];
      R = con->R;
      switch (con->kind)
      {
      case CONTACT:
      case FIXPNT:
      case GLUE:
        COPY (D, R);
	break;
      case VELODIR:
      case FIXDIR:
      case RIGLNK:
	R [2] = D [0];
	break;
      }
    }
#if MPI
  }
  else
  {
    A->nsend = 0;
    A->ssend = dom->ncon + 8;
    MEM_Init (&A->mem, sizeof (int [2]), 128);
    ERRMEM (A->send = malloc (sizeof (COMDATA [A->ssend])));

    for (bod = dom->bod; bod; bod = bod->next) /* for all parents */
    {
      for (item = SET_First (bod->con); item; item = SET_Next (item))
      {
	con = item->data;
	if (con->state & CON_EXTERNAL) /* needs local velocity update */
	{
	  ASSERT_DEBUG (MAP_Find (dom->conext, (void*) (long) con->id, NULL), "Invalid external constraint %s %d", CON_Kind (con), con->id);
	  ptr = send_item (&A->send, &A->ssend, &A->nsend);
	  ptr->rank = con->rank;
	  ERRMEM (ptr->i = MEM_Alloc (&A->mem));
	  ptr->i [0] = -con->id; /* negative indicates local velocity update */
	  ptr->ints = 2;
	  ptr->d = con->U;
	  ptr->doubles = 3;

	  if (bod == con->master)
	  {
	    ptr->i [1] = 1; /* master */
	    BODY_Local_Velo (bod, mshp(con), mgobj(con), con->mpnt, con->base, NULL, con->U);
	  }
	  else
	  {
	    ptr->i [1] = 0; /* slave */
	    BODY_Local_Velo (bod, sshp(con), sgobj(con), con->spnt, con->base, NULL, con->U);
	  }
	}
      }
    }

    for (con = dom->con; con; con = con->next) /* for all constraints */
    {
      for (item = SET_First (con->ext); item; item = SET_Next (item))
      {
	ptr = send_item (&A->send, &A->ssend, &A->nsend);
	ptr->rank = (int) (long) item->data; /* external rank */
	ERRMEM (ptr->i = MEM_Alloc (&A->mem));
	ptr->i [0] = con->id; /* positive indicates external reaction update */
	ptr->ints = 2;
	ptr->d = con->R;
	ptr->doubles = 3;
      }
    }

    A->pattern = COMALL_Pattern (MPI_COMM_WORLD, A->send, A->nsend, &A->recv, &A->nrecv);

    return A->pattern; /* initialized */
  }

  for (bod = dom->bod, u = x; bod; u += bod->dofs, bod = bod->next)
  {
    velo = bod->velo;
    bod->velo = u;
    for (item = SET_First (bod->con); item; item = SET_Next (item))
    {
      con = item->data;
      if (con->state & CON_EXTERNAL) /* update local velocities of external constraints attached to parents, before sending */
      {
	if (bod == con->master) BODY_Local_Velo (bod, mshp(con), mgobj(con), con->mpnt, con->base, NULL, con->U);
	else BODY_Local_Velo (bod, sshp(con), sgobj(con), con->spnt, con->base, NULL, con->U);
      }
    }
    bod->velo = velo;
  }

  COMALL_Repeat (A->pattern); /* send and receive */

  for (i = 0; i < A->nrecv; i ++)
  {
    ptr = &A->recv [i];
    for (j = ptr->i, k = j + ptr->ints, D = ptr->d; j < k; j += 2, D += 3)
    {
      if ((*j) < 0) /* local velocity */
      {
	ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*) (long) (-(*j)), NULL), "Invalid constraint id: %d", -(*j));
	if ((*(j+1))) /* master */
	{
	  U = con->U;
	  COPY (D, U);
	}
	else /* slave */
	{
	  B = con->dia->B;
	  COPY (D, B);
	}
      }
      else /* external reaction */
      {
	ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*j), NULL), "Invalid external constraint id: %d", (*j));
	R = con->R;
	COPY (D, R);
      }
    }
  }
#endif

  /* compute relative velocities */
  for (con = dom->con; con; con = con->next)
  {
    if (con->slave)
    {
      U = con->U;
      B = con->dia->B;

      SUB (U, B, U); /* relative = master - slave */
    }
  }

  return NULL;
}

/* update previous local velocities */
static void update_previous_local_velocities (DOM *dom)
{
  double *V, *B;
  CON *con;
#if MPI
  int nsend, nrecv, *isize, *dsize, i, *j, *k;
  COMDATA *send, *recv, *ptr;
  double X [3], *D;
  BODY *bod;
  SET *item;

  nsend = dom->ncpu;
  ERRMEM (send = MEM_CALLOC (sizeof (COMDATA [nsend])));
  ERRMEM (isize = MEM_CALLOC (sizeof (int [nsend])));
  ERRMEM (dsize = MEM_CALLOC (sizeof (int [nsend])));

  for (i = 0; i < nsend; i ++) send [i].rank = i;

  for (bod = dom->bod; bod; bod = bod->next) /* for all parents */
  {
    for (item = SET_First (bod->con); item; item = SET_Next (item))
    {
      con = item->data;
      if (con->state & CON_EXTERNAL) /* needs local velocity update */
      {
	ASSERT_DEBUG (MAP_Find (dom->conext, (void*) (long) con->id, NULL), "Invalid external constraint %s %d", CON_Kind (con), con->id);
	i = con->rank;
	ptr = &send [i];
	if (bod == con->master)
	{
	  pack_int (&isize [i], &ptr->i, &ptr->ints, con->id);
	  BODY_Local_Velo (bod, mshp(con), mgobj(con), con->mpnt, con->base, X, NULL);
	}
	else
	{
	  pack_int (&isize [i], &ptr->i, &ptr->ints, -con->id);
	  BODY_Local_Velo (bod, sshp(con), sgobj(con), con->spnt, con->base, X, NULL);
	}
        pack_doubles (&dsize [i], &ptr->d, &ptr->doubles, X, 3);
      }
    }
  }

  COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);

  for (i = 0; i < nrecv; i ++)
  {
    ptr = &recv [i];
    for (j = ptr->i, k = j + ptr->ints, D = ptr->d; j < k; j ++, D += 3)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*) (long) ABS (*j), NULL), "Invalid constraint id: %d", *j);
      if ((*j) < 0) /* slave */
      {
        B = con->dia->B;
        COPY (D, B);
      }
      else /* master */
      {
        V = con->dia->V;
        COPY (D, V);
      }
    }
  }

  for (i = 0; i < nsend; i ++)
  {
    ptr = &send [i];
    free (ptr->i);
    free (ptr->d);
  }
  free (send);
  free (isize);
  free (dsize);
  free (recv); /* includes recv[]->i and recv[]->d memory */
#endif

  /* remote master and slave previous local velocities are now stored in V and B members of con->dia;
   * compute the final values of con->dia->V, completing the local computations if needed */
  for (con = dom->con; con; con = con->next)
  {
    V = con->dia->V;

#if MPI
    if (con->master->flags & BODY_PARENT) /* local parent */
#endif
    {
      BODY_Local_Velo (con->master, mshp(con), mgobj(con), con->mpnt, con->base, V, NULL);
    }

    if (con->slave)
    {
      B = con->dia->B;

#if MPI
      if (con->slave->flags & BODY_PARENT) /* local slave */
#endif
      {
        BODY_Local_Velo (con->slave, sshp(con), sgobj(con), con->spnt, con->base, B, NULL);
      }

      SUB (V, B, V); /* relative = master - slave */
    }
  }
}

/* create BSS data */
static BSS_DATA *create_data (DOM *dom)
{
  short dynamic;
  double *b, *x;
  BSS_DATA *A;
  BODY *bod;
  CON *con;
  int n, m;

  ERRMEM (A = MEM_CALLOC (sizeof (BSS_DATA)));
  dynamic = dom->dynamic;
  A->dom = dom;
  for (bod = dom->bod, m = 0; bod; bod = bod->next)
  {
    n = bod->dofs;
    A->nprimal += n;
    if (n > m) m = n;
  }
  for (con = dom->con; con; con = con->next)
  {
    switch (con->kind)
    {
    case CONTACT:
    case FIXPNT:
    case GLUE: A->ndual += 3; break;
    case FIXDIR:
    case VELODIR:
    case RIGLNK: SET2 (con->R, 0.0); A->ndual += 1; break; /* only normal component */
    }
  }
  A->x = newvector (A->nprimal + A->ndual);
  A->b = newvector (A->nprimal + A->ndual);
  ERRMEM (A->r = MEM_CALLOC (sizeof (double [m])));

  /* initialize free momentum right hand side */
  for (bod = dom->bod, b = A->b->x; bod; b += bod->dofs, bod = bod->next)
  {
    BODY_Matvec (1.0, bod, bod->velo, 0.0, b);
  }

  /* update previous local velocities */
  if (dynamic) update_previous_local_velocities (dom);

  /* set up constraint right hand side and solution pointers and set constant values */
  for (con = dom->con, b = &A->b->x [A->nprimal], x = &A->x->x [A->nprimal]; con; con = con->next)
  {
    double *V = con->dia->V;

    CON_RHS_SET (con, b);
    CON_X_SET (con, x);

    switch ((int)con->kind)
    {
    case FIXPNT:
    case GLUE:
    {
      if (dynamic)
      {
	b [0] = -V[0];
	b [1] = -V[1];
	b [2] = -V[2];
      }
      else
      {
	b [0] = 0;
	b [1] = 0;
	b [2] = 0;
      }

      b += 3; /* increment */
      x += 3;
    }
    break;
    case FIXDIR:
    {
      if (dynamic) b [0] = -V[2];
      else b [0] = 0;

      b += 1; /* increment */
      x += 1;
    }
    break;
    case VELODIR:
    {
      b [0] = VELODIR(con->Z);

      b += 1; /* increment */
      x += 1;
    }
    break;
    case CONTACT:
    {
      b += 3; /* increment */
      x += 3;
    }
    break;
    case RIGLNK:
    {
      b += 2; /* increment */
      x += 2;
    }
    break;
    }
  }

#if MPI
  /* create communication pattern */
  A->pattern = update_constraints_data (A, NULL);
#endif

  return A;
}

/* imaginary i */
static double complex imaginary_i;

/* real normal to friction cone */
inline static void real_n (double *S, double fri, double *n)
{
  double dot, len;

  dot = DOT2(S, S);
  len = sqrt(dot);

  if (len == 0 || len <= fri * S[2])
  {
    SET (n, 0.0);
  }
  else if (fri * len + S[2] < 0.0)
  {
    dot += S[2]*S[2];
    len = sqrt (dot);
    if (len == 0) { SET (n, 0.0); }
    else { DIV (S, len, n); }
  }
  else
  {
    dot = 1.0 / sqrt (1.0 + fri*fri);
    DIV2 (S, len, n);
    n [2] = -fri;
    SCALE (n, dot);
  }
}

/* complex normal to friction cone */
inline static void complex_n (double complex *S, double complex fri, double complex *n)
{
  double complex dot, len;

  dot = DOT2(S, S);
  len = csqrt(dot);

  if (creal (len) == 0 || creal (len) <= creal (fri * S[2]))
  {
    SET (n, 0.0 + 0.0 * imaginary_i);
  }
  else if (creal (fri * len + S[2]) < 0.0)
  {
    dot += S[2]*S[2];
    len = csqrt (dot);
    if (creal (len) == 0) { SET(n, 0.0 + 0.0 * imaginary_i); }
    DIV (S, len, n);
  }
  else
  {
    dot = 1.0 / csqrt (1.0 + fri*fri);
    DIV2 (S, len, n);
    n [2] = -fri;
    SCALE (n, dot);
  }
}

/* real normal ray to friction cone */
inline static void real_m (double fri, short smooth, double *S, double eps, double *m)
{
  double n [3], fun;

  real_n (S, fri, n);
  fun = DOT (S, n);

  if (smooth == 1 && fun >= 0.0 && fun <= eps)
  {
    fun = ((2.0/eps) - (1.0/(eps*eps))*fun)*(fun*fun);
  }
  else if (smooth == 2)
  {
    if (fun >= 0.0 && fun <= eps)
    {
      fun = (fun*fun) / (2.0 * eps);
    }
    else if (fun > eps)
    {
      fun = fun - 0.5 * eps;
    }
  }

  MUL (n, fun, m)
}

/* complex normal ray to friction cone */
inline static void complex_m (double complex fri, short smooth, double complex *S, double complex eps, double complex *m)
{
  double complex n [3], fun;

  complex_n (S, fri, n);
  fun = DOT (S, n);

  if (smooth == 1 && creal (fun) >= 0.0 && creal (fun) <= creal (eps))
  {
    fun = ((2.0/eps) - (1.0/(eps*eps))*fun)*(fun*fun);
  }
  else if (smooth == 2)
  {
    if (creal (fun) >= 0.0 && creal (fun) <= creal (eps))
    {
      fun = (fun*fun) / (2.0 * eps);
    }
    else if (creal (fun) > creal (eps))
    {
      fun = fun - 0.5 * eps;
    }
  }

  MUL (n, fun, m)
}

/* real F = [UT, UN + fri |UT|]' */
inline static void real_F (double res, double fri, double gap, double step, short dynamic, double epsilon, double *V, double *U, double *F)
{
  double udash;

  if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
  else udash = ((MAX(gap, 0)/step) + U[2]);

#if DISABLE_NORM_SMOOTHING
  epsilon = 0;
#endif

  F [0] = U[0];
  F [1] = U[1];
  F [2] = (udash + fri * sqrt (DOT2(U, U) + epsilon*epsilon));
}
 
/* complex F = [UT, UN + fri |UT|]' */
inline static void complex_F (double res, double fri, double gap, double step, short dynamic, double epsilon, double *V, double complex *U, double complex *F)
{
  double complex udash;

  if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
  else udash = ((MAX(gap, 0)/step) + U[2]);

#if DISABLE_NORM_SMOOTHING
  epsilon = 0;
#endif

  F [0] = U[0];
  F [1] = U[1];
  F [2] = (udash + fri * csqrt (DOT2(U, U) + epsilon*epsilon));
}

/* update linear system */
static void update_system (BSS_DATA *A)
{
  double h, step;
  short dynamic;
  DOM *dom;
  CON *con;

  dom = A->dom;
  dynamic = dom->dynamic;
  step = dom->step;
  h = DIFF_FACTOR * SMOOTHING_EPSILON;
  imaginary_i = csqrt (-1);

  /* contact constraints linearisation */
  for (con = dom->con; con; con = con->next)
  {
    DIAB *dia = con->dia;
    double *U = con->U,
	   *R = con->R,
	   *b = CON_RHS_GET (con);

    switch ((int)con->kind)
    {
    case RIGLNK:
    {
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
    }
    break;
    case CONTACT:
    {
      double *V = dia->V,
	     *X = CONTACT_X (con),
	     *Y = CONTACT_Y (con),
	     gap = con->gap,
	     fri = con->mat.base->friction,
	     res = con->mat.base->restitution,
	     dF[9],
	     F [3],
             S [3],
	     m [3],
	     H [3],
	     J [9];

      double complex cU [3],
	             cS [3],
		     cF [3],
		     cm [3];


      real_F (res, fri, gap, step, dynamic, SMOOTHING_EPSILON, V, U, F);
      SUB (R, F, S);
      real_m (fri, SMOOTHING, S, SMOOTHING_EPSILON, m);
      ADD (F, m, H);

      for (int k = 0; k < 3; k ++)
      {
	cU [0] = U[0] + 0.0 * imaginary_i;
	cU [1] = U[1] + 0.0 * imaginary_i;
	cU [2] = U[2] + 0.0 * imaginary_i;
	cU [k] += h * imaginary_i;
        complex_F (res, fri, gap, step, dynamic, SMOOTHING_EPSILON, V, cU, cF);
	dF [3*k+0] = cimag (cF [0]) / h;
	dF [3*k+1] = cimag (cF [1]) / h;
	dF [3*k+2] = cimag (cF [2]) / h;

        cS [0] = S[0] + 0.0 * imaginary_i;
	cS [1] = S[1] + 0.0 * imaginary_i;
	cS [2] = S[2] + 0.0 * imaginary_i;
	cS [k] += h * imaginary_i;
        complex_m (fri, SMOOTHING, cS, SMOOTHING_EPSILON, cm);
	Y [3*k+0] = cimag (cm [0]) / h; /* Y = dm/dS */
	Y [3*k+1] = cimag (cm [1]) / h;
	Y [3*k+2] = cimag (cm [2]) / h;
      }

      IDENTITY (J);
      NNSUB (J, Y, J);
      NNMUL (dF, J, X); /* X = dF/dU [I - dm/dS] */

      NVMUL (X, U, b);
      NVADDMUL (b, Y, R, b);
      SUB (b, H, b); /* b = X U + Y R - H(U,R) */
    }
    break;
    }
  }
}

/* linear dual residual norm */
static double resnorm (BSS_DATA *A)
{
  double G [3], *U, *R, *X, *Y, *b, dot;
  CON *con;

  for (dot = 0, con = A->dom->con; con; con = con->next)
  {
    b = CON_RHS_GET (con);
    U = con->U;
    R = con->R;

    switch (con->kind)
    {
    case FIXPNT:
    case GLUE:
      SUB (U, b, G);
      dot += DOT (G, G);
      break;
    case FIXDIR:
    case VELODIR:
      G [0] = U[2] - b[0];
      dot += G[0]*G[0];
      break;
    case RIGLNK:
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
      break;
    case CONTACT:
      X = CONTACT_X (con);
      Y = CONTACT_Y (con);
      NVMUL (X, U, G);
      NVADDMUL (G, Y, R, G);
      SUB (G, b, G);
      dot += DOT (G, G);
      break;
    }
  }

  return sqrt (dot);
}

/* solve linear system */
static void linear_solve (BSS_DATA *A, double resdec, int maxiter)
{
  hypre_FlexGMRESFunctions *gmres_functions;
  void *gmres_vdata;
  double abstol;

  abstol = resdec * A->resnorm;

  if (abstol == 0.0) /* initially */
  {
    abstol = ABSTOL * sqrt (InnerProd (A->b, A->b));
    if (abstol == 0.0) abstol = ABSTOL;
  }

  gmres_functions = hypre_FlexGMRESFunctionsCreate (CAlloc, Free, (int (*) (void*,int*,int*)) CommInfo,
    (void* (*) (void*))CreateVector, (void* (*) (int, void*))CreateVectorArray, (int (*) (void*))DestroyVector,
    MatvecCreate, (int (*) (void*,double,void*,void*,double,void*))Matvec, MatvecDestroy,
    (double (*) (void*,void*))InnerProd, (int (*) (void*,void*))CopyVector, (int (*) (void*))ClearVector,
    (int (*) (double,void*))ScaleVector, (int (*) (double,void*,void*))Axpy,
    PrecondSetup, (int (*) (void*,void*,void*,void*))Precond);
  gmres_vdata = hypre_FlexGMRESCreate (gmres_functions);

  hypre_FlexGMRESSetTol (gmres_vdata, 0.0);
  hypre_FlexGMRESSetMinIter (gmres_vdata, 1);
  hypre_FlexGMRESSetMaxIter (gmres_vdata, maxiter);
  hypre_FlexGMRESSetAbsoluteTol (gmres_vdata, abstol);
  hypre_FlexGMRESSetup (gmres_vdata, A, A->b, A->x);
  hypre_FlexGMRESSolve (gmres_vdata, A, A->b, A->x);
  hypre_FlexGMRESGetNumIterations (gmres_vdata , &A->iters);
  hypre_FlexGMRESDestroy (gmres_vdata);

  A->resnorm = resnorm (A);
}

/* update solution */
static void update_solution (BSS_DATA *A)
{
  double m [3], *R, *x, fri;
  CON *con;

  for (con = A->dom->con; con; con = con->next)
  {
    x = CON_X_GET (con);
    R = con->R;

    switch (con->kind)
    {
    case CONTACT:
    {
      fri = con->mat.base->friction;
      real_m (fri, 0, x, 0, m);
      SUB (x, m, R); /* project onto the friction cone */
    }
    break;
    case FIXPNT:
    case GLUE:
      COPY (x, R);
      break;
    case FIXDIR:
    case VELODIR:
    case RIGLNK:
      R [2] = x[0];
      break;
    }
  }
}

/* compute merit function value */
static double merit_function (BSS_DATA *A)
{
  double G [3], *U, *b, step, dot;
  DOM *dom = A->dom;
  short dynamic;
  CON *con;

  dynamic = dom->dynamic;
  step = dom->step;

  for (dot = 0, con = dom->con; con; con = con->next)
  {
    b = CON_RHS_GET (con);
    U = con->U;

    switch (con->kind)
    {
    case FIXPNT:
    case GLUE:
      SUB (U, b, G);
      dot += DOT (G, G);
      break;
    case FIXDIR:
    case VELODIR:
      G [0] = U[2] - b[0];
      dot += G[0]*G[0];
      break;
    case RIGLNK:
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
      break;
    case CONTACT:
      {
	 double *V = con->dia->V,
		*R = con->R,
	         gap = con->gap,
	         fri = con->mat.base->friction,
	         res = con->mat.base->restitution,
	         F [3],
	         S [3],
	         m [3];

	real_F (res, fri, gap, step, dynamic, SMOOTHING_EPSILON, V, U, F);
	SUB (R, F, S);
	real_m (fri, SMOOTHING, S, SMOOTHING_EPSILON, m);
	ADD (F, m, G);
        dot += DOT (G, G);
      }
      break;
    }
  }

  return sqrt (dot);
}

/* destroy BSS data */
static void destroy_data (BSS_DATA *A)
{
  DestroyVector (A->b);
  DestroyVector (A->x);
  free (A->r);

#if MPI
  MEM_Release (&A->mem);
  free (A->send);
  free (A->recv);
#endif

  free (A);
}

/* create solver */
BSS* BSS_Create (double meritval, int maxiter)
{
  BSS *bs;

  ERRMEM (bs = MEM_CALLOC (sizeof (BSS)));
  bs->meritval = meritval;
  bs->maxiter = maxiter;
  bs->linminiter = 5;
  bs->resdec = 0.25;
  bs->verbose = 1;

  return bs;
}

/* run solver */
void BSS_Solve (BSS *bs, LOCDYN *ldy)
{
  char fmt [512];
  double *merit;
  BSS_DATA *A;
  DOM *dom;

  sprintf (fmt, "BODY_SPACE_SOLVER: (LIN its/res: %%%dd/%%.2e) iteration: %%%dd merit: %%.2e\n",
           (int)log10 (bs->linminiter * bs->maxiter) + 1, (int)log10 (bs->maxiter) + 1);

  ERRMEM (bs->merhist = realloc (bs->merhist, bs->maxiter * sizeof (double)));
  dom = ldy->dom;
  A = create_data (dom);
  merit = &dom->merit;
  bs->iters = 0;

  do
  {
    update_system (A);

    linear_solve (A, bs->resdec, bs->linminiter + bs->iters);

    update_solution (A);

    *merit = merit_function (A);

    bs->merhist [bs->iters] = *merit;

#if MPI
    if (dom->rank == 0)
#endif
    if (dom->verbose && bs->verbose) printf (fmt, A->iters, A->resnorm, bs->iters, *merit);

  } while (++ bs->iters < bs->maxiter && *merit > bs->meritval);

  destroy_data (A);
}

/* write labeled state values */
void BSS_Write_State (BSS *bs, PBF *bf)
{
  /* TODO */
}

/* destroy solver */
void BSS_Destroy (BSS *bs)
{
  free (bs->merhist);
  free (bs);
}
