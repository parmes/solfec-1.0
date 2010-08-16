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

#define ABSTOL_BASE 1E-15
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

  double *r, /* body space reaction of size MAX { bod->dofs } */
          delta, /* dual regularisation */
          resnorm, /* dual residual norm */
          xnorm; /* dual solution norm */

  int iters; /* number of linear solver iterations */

#if MPI
  COMDATA send, recv; /* communication buffers */
  int nsend, nrecv; /* sizes */
  void *pattern; /* and pattern */
#endif
};

struct vector
{
  double *x;

  int n;
};

#if MPI
/* update needed before matrix vector product */
static void* update_constraints_data (BSS_DATA *A, double *x);
#endif

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
  double *U, *R, *X, *Y, *a, *b, *r, step;
  BODY *bod;
  CON *con;
  DOM *dom;
  int n;

#if MPI
  update_constraints_data (A, x->x);
#endif

  dom = A->dom;
  step = dom->step;

  ScaleVector (beta, y);

  if (A->delta > 0.0)
  {
    blas_daxpy (A->ndual, A->delta, &x->x[A->nprimal], 1, &y->x[A->nprimal], 1);
  }

  for (bod = dom->bod, r = A->r, a = x->x, b = y->x; bod; a += n, b += n, bod = bod->next)
  {
    n = bod->dofs;
    BODY_Matvec (1.0, bod, a, 1.0, b);
    BODY_Reac (bod, r);
    blas_daxpy (n, -step, r, 1, b, 1);
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
      b [0] += U[0];
      b [1] += U[1];
      b [2] += U[2];
      break;
    case FIXDIR:
    case VELODIR:
      b [0] += U[2];
      break;
    case RIGLNK:
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
      break;
    case CONTACT:
      X = CONTACT_X (con);
      Y = CONTACT_Y (con);
      NVADDMUL (b, X, U, b);
      NVADDMUL (b, Y, R, b);
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
  double T [9], *Y, *p, *q, *r, delta;
  int ipiv [3], n;
  BODY *bod;
  CON *con;
  DOM *dom;

  delta = A->delta;
  dom = A->dom;


  for  (bod = dom->bod, p = b->x, q = x->x; bod; p += n, q += n, bod = bod->next)
  {
    n = bod->dofs;
    BODY_Invvec (1.0, bod, p, 0.0, q);
  }

  for (con = dom->con, r = &A->b->x [A->nprimal]; con; con = con->next)
  {
    n = CON_RHS_GET (con) - r;
    p = &b->x [n];
    q = &x->x [n];

    COPY (p, q); /* A == identity so far */

    switch ((int)con->kind)
    {
    case RIGLNK:
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
      break;
    case CONTACT:
      Y = CONTACT_Y (con);
      NNCOPY (Y, T);
      T [0] += delta;
      T [4] += delta;
      T [8] += delta;
      lapack_dgesv (3, 1, T, 3, ipiv, q, 3);
      break;
    }
  }

  return 0;
}
/* GMRES interface end */

#if MPI
/* update constraints velocities and reactions or
 * create the necessary for that communication pattern */
static void* update_constraints_data (BSS_DATA *A, double *x)
{
  if (A->pattern)
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  else
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }

  return A->pattern;
}
#endif

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
      x += 1;
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

/* update linear system */
static void update_system (BSS_DATA *A)
{
  short dynamic;
  DOM *dom;
  CON *con;

  dom = A->dom;
  dynamic = dom->dynamic;

  /* contact constraints linearisation */
  for (con = dom->con; con; con = con->next)
  {
    switch ((int)con->kind)
    {
    case RIGLNK:
    {
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
    }
    break;
    case CONTACT:
    {
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
    }
    break;
    }
  }
}

/* compute residual and solution norms */
static void compute_norms (BSS_DATA *A)
{
  double G [3], *U, *R, *X, *Y, *b, rdot, xdot;
  CON *con;

  for (rdot = xdot = 0, con = A->dom->con; con; con = con->next)
  {
    b = CON_RHS_GET (con);
    U = con->U;
    R = con->R;

    switch (con->kind)
    {
    case FIXPNT:
    case GLUE:
      SUB (U, b, G);
      rdot += DOT (G, G);
      break;
    case FIXDIR:
    case VELODIR:
      G [0] = U[2] - b[0];
      rdot += G[0]*G[0];
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
      rdot += DOT (G, G);
      break;
    }

    xdot += DOT (R, R);
  }

  A->resnorm = sqrt (rdot);
  A->xnorm = sqrt (xdot);
}

/* solve linear system */
static void linear_solve (BSS_DATA *A, double resdec, int maxiter)
{
  hypre_FlexGMRESFunctions *gmres_functions;
  void *gmres_vdata;
  double abstol;

  A->delta = A->resnorm / A->xnorm; /* sqrt (L-curve) */

  abstol = resdec * A->resnorm;

  if (abstol == 0.0) /* initially */
  {
    abstol = ABSTOL_BASE * sqrt (InnerProd (A->b, A->b));
    if (abstol == 0.0) abstol = ABSTOL_BASE;
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

  compute_norms (A);
}

/* update solution */
static void update_solution (BSS_DATA *A)
{
  double *R, *x;
  CON *con;

  for (con = A->dom->con; con; con = con->next)
  {
    x = CON_X_GET (con);
    R = con->R;

    switch (con->kind)
    {
    case CONTACT:
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
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  return 0.0;
}

/* destroy BSS data */
static void destroy_data (BSS_DATA *A)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
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
