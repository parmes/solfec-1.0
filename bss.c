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
#include "vic.h"
#include "ext/krylov/krylov.h"

#if MPI
#include "com.h"
#include "pck.h"
#endif

#define SMOOTHING_EPSILON	1E-10  /* TODO */
#define ABSTOL			1E-15  /* TODO */
#define CON_RHS_SET(con, value) (con)->dia->mprod = (MX*)(value)
#define CON_RHS(con) ((double*)(con)->dia->mprod)
#define CON_X_SET(con, value) (con)->dia->sprod = (MX*)(value)
#define CON_X(con) ((double*)(con)->dia->sprod)
#define CON_INDEX_SET(con, value) ((con)->dia->mH = (MX*)value)
#define CON_INDEX(con) (int)(long)(con)->dia->mH
#define CONTACT_X(con) (con)->dia->W
#define CONTACT_Y(con) (con)->dia->A

typedef struct bss_data BSS_DATA;
typedef struct vector VECTOR;

struct bss_data
{
  DOM *dom;

  int ndual; /* SUM [con in dom->con] { size (con->R) } */

  VECTOR *b, /* right hand side */
         *x, /* reactions */
	 *z; /* increments of reactions */

  double *r, /* reaction workspace of size MAX [bod in dom->bod] { bod->dofs } */
	 *u, /* global replacement velocity of size SUM [bod in dom->bod] { bod->dofs } */
	 *W, /* approximate diagonal W blocks */
          resnorm, /* linear residual norm */
	  xnorm, /* x norm */
	  delta; /* regularisation */

  int iters; /* number of linear solver iterations */

#if MPI
  MEM mem; /* communication memory */
  COMDATA *send, *recv; /* buffers */
  int ssend, nsend, nrecv; /* sizes */
  void *pattern; /* and pattern */
#endif
};

struct vector
{
  double *x;

  int n;
};

/* needed in matrix vector product */
static void* update_local_velocities (BSS_DATA *A);

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
  double c [3], *U, *R, *X, *Y, *b, *r, *u, step, delta;
  DOM *dom = A->dom;
  BODY *bod;
  CON *con;
  int n;

  step = dom->step;
  delta = A->delta;

  ScaleVector (beta, y);

  if (delta > 0.0) Axpy (alpha*delta, x, y);

  for (con = dom->con; con; con = con->next) /* copy x to con->R */
  {
    r = &x->x [CON_INDEX (con)];
    R = con->R;

    switch (con->kind)
    {
    case CONTACT:
    case FIXPNT:
    case GLUE:
      COPY (r, R);
      break;
    case VELODIR:
    case FIXDIR:
    case RIGLNK:
      R [2] = r [0];
      break;
    }
  }

#if MPI
  DOM_Update_External_Reactions  (dom, 0);
#endif

  for (bod = dom->bod, r = A->r, u = A->u; bod; u += n, bod = bod->next)
  {
    n = bod->dofs;
    BODY_Reac (bod, r);
    BODY_Invvec (step, bod, r, 0.0, u); /* would be velocity (***) */
  }

  update_local_velocities (A); /* compute con->U based on (***) */

  for (con = dom->con; con; con = con->next)
  {
    b = &y->x [CON_INDEX (con)];
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
#if 1
  double *W, *X, *Y, *bx, *xx, delta, T [9];
  int ipiv [3];
  CON *con;

  delta = A->delta;

  for (con = A->dom->con, W = A->W, bx = b->x, xx = x->x; con; con = con->next, W += 9)
  {
    switch (con->kind)
    {
    case FIXPNT:
    case GLUE:
      NNCOPY (W, T);
      T [0] += delta;
      T [4] += delta;
      T [8] += delta;
      COPY (bx, xx);
      lapack_dgesv (3, 1, T, 3, ipiv, xx, 3);
      bx += 3;
      xx += 3;
      break;
    case VELODIR:
    case FIXDIR:
      T[8] = W[8] + delta;
      xx [0] = bx [0] / T[8];
      bx += 1;
      xx += 1;
      break;
    case RIGLNK:
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
      break;
    case CONTACT:
      X = CONTACT_X (con);
      Y = CONTACT_Y (con);
      NNMUL (X, W, T);
      NNADD (T, Y, T);
      T [0] += delta;
      T [4] += delta;
      T [8] += delta;
      COPY (bx, xx);
      lapack_dgesv (3, 1, T, 3, ipiv, xx, 3);
      bx += 3;
      xx += 3;
      break;
    }
  }
#else
  CopyVector (b, x);
#endif

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

/* update local constraints velocities or create
 * the necessary for that communication pattern */
static void* update_local_velocities (BSS_DATA *A)
{
  double *U, *B, *u, *velo;
  DOM *dom = A->dom;
  BODY *bod;
  SET *item;
  CON *con;
#if MPI
  int i, *j, *k;
  COMDATA *ptr;
  double *D;
#endif

#if MPI
  if (A->pattern)
  {
#endif
    /* update local velocities of local constraints */
    for (bod = dom->bod, u = A->u; bod; u += bod->dofs, bod = bod->next)
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
	  else BODY_Local_Velo (bod, sshp(con), sgobj(con), con->spnt, con->base, NULL, con->dia->B+3);
#if MPI
	}
#endif
      }
      bod->velo = velo;
    }
#if MPI
  }
  else
  {
    A->nsend = 0;
    A->ssend = dom->ncon + 8;
    MEM_Init (&A->mem, sizeof (int), 128);
    ERRMEM (A->send = malloc (sizeof (COMDATA [A->ssend])));

    for (bod = dom->bod; bod; bod = bod->next) /* for all parents */
    {
      for (item = SET_First (bod->con); item; item = SET_Next (item))
      {
	con = item->data;
	if (con->state & CON_EXTERNAL) /* parent needs local velocity update */
	{
	  ASSERT_DEBUG (MAP_Find (dom->conext, (void*) (long) con->id, NULL), "Invalid external constraint %s %d", CON_Kind (con), con->id);
	  ptr = send_item (&A->send, &A->ssend, &A->nsend);
	  ptr->rank = con->rank;
	  ERRMEM (ptr->i = MEM_Alloc (&A->mem));
	  ptr->ints = 1;
	  ptr->d = con->U;
	  ptr->doubles = 3;

	  if (bod == con->master)
	  {
	    ptr->i [0] = -con->id; /* negative indicates master */
	  }
	  else
	  {
	    ptr->i [0] = con->id; /* positive indicates slave */
	  }
	}
      }
    }

    A->pattern = COMALL_Pattern (MPI_COMM_WORLD, A->send, A->nsend, &A->recv, &A->nrecv);

    return A->pattern; /* initialized */
  }

  for (bod = dom->bod, u = A->u; bod; u += bod->dofs, bod = bod->next)
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
    for (j = ptr->i, k = j + ptr->ints, D = ptr->d; j < k; j ++, D += 3)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*) (long) ABS(*j), NULL), "Invalid constraint id: %d", ABS(*j));

      if ((*j) < 0) /* master */
      {
	U = con->U;
	COPY (D, U);
      }
      else /* slave */
      {
	B = con->dia->B+3;
	COPY (D, B);
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
      B = con->dia->B+3;

      SUB (U, B, U); /* relative = master - slave */
    }
  }

  return NULL;
}

/* update previous and free local velocities */
static void update_const_local_velocities (DOM *dom)
{
  double *V, *B;
  CON *con;
#if MPI
  int nsend, nrecv, *isize, *dsize, i, *j, *k;
  COMDATA *send, *recv, *ptr;
  double X [6], *D;
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
	  pack_int (&isize [i], &ptr->i, &ptr->ints, -con->id);
	  BODY_Local_Velo (bod, mshp(con), mgobj(con), con->mpnt, con->base, X, X+3);
	}
	else
	{
	  pack_int (&isize [i], &ptr->i, &ptr->ints, con->id);
	  BODY_Local_Velo (bod, sshp(con), sgobj(con), con->spnt, con->base, X, X+3);
	}
        pack_doubles (&dsize [i], &ptr->d, &ptr->doubles, X, 6);
      }
    }
  }

  COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);

  for (i = 0; i < nrecv; i ++)
  {
    ptr = &recv [i];
    for (j = ptr->i, k = j + ptr->ints, D = ptr->d; j < k; j ++, D += 6)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*) (long) ABS(*j), NULL), "Invalid constraint id: %d", ABS(*j));
      V = con->dia->V;
      B = con->dia->B;
      if ((*j) < 0) /* master */
      {
        COPY (D, V);
        COPY (D+3, B);
      }
      else /* slave */
      {
        COPY (D, V+3);
        COPY (D+3, B+3);
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
    B = con->dia->B;

#if MPI
    if (con->master->flags & BODY_PARENT) /* local parent */
#endif
    {
      BODY_Local_Velo (con->master, mshp(con), mgobj(con), con->mpnt, con->base, V, B);
    }

    if (con->slave)
    {
#if MPI
      if (con->slave->flags & BODY_PARENT) /* local slave */
#endif
      {
        BODY_Local_Velo (con->slave, sshp(con), sgobj(con), con->spnt, con->base, V+3, B+3);
      }

      SUB (V, V+3, V); /* relative = master - slave */
      SUB (B, B+3, B);
    }
  }
}

/* approximate diagonal W blocks */
static void approximate_diagobal_blocks (CON *con, double *W)
{
  double step = con ? con->master->dom->step : 0;

  for (; con; con = con->next, W += 9)
  {
    BODY *m = con->master,
	 *s = con->slave;
    void *mgobj = mgobj(con),
	 *sgobj;
    SHAPE *mshp = mshp(con),
	  *sshp;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base;
    MX_DENSE_PTR (A, 3, 3, W);
    MX *H, *prod;

    if (s)
    {
      sgobj = sgobj(con);
      sshp = sshp(con);
    }

    H = BODY_Gen_To_Loc_Operator (m, mshp, mgobj, mpnt, base);
    prod = MX_Matmat (1.0, m->inverse, MX_Tran (H), 0.0, NULL);
    MX_Matmat (step, H, prod, 0.0, &A); /* H * inv (M) * H^T */
    MX_Destroy (prod);
    MX_Destroy (H);

    if (s)
    {
      H = BODY_Gen_To_Loc_Operator (s, sshp, sgobj, spnt, base);
      prod = MX_Matmat (1.0, s->inverse, MX_Tran (H), 0.0, NULL);
      MX_Matmat (step, H, prod, 1.0, &A); /* H * inv (M) * H^T */
      MX_Destroy (prod);
      MX_Destroy (H);
    }
  }
}

/* create BSS data */
static BSS_DATA *create_data (DOM *dom)
{
  double *R, *U, *B, *b, *x, *r, *u, step;
  int n, m, nprimal;
  short dynamic;
  BSS_DATA *A;
  BODY *bod;
  CON *con;

  ERRMEM (A = MEM_CALLOC (sizeof (BSS_DATA)));
  dynamic = dom->dynamic;
  step = dom->step;
  A->xnorm = 1.0;
  A->dom = dom;
  for (bod = dom->bod, nprimal = m = 0; bod; bod = bod->next)
  {
    n = bod->dofs;
    nprimal += n;
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
    case VELODIR: SET2 (con->R, 0.0); A->ndual += 1; break; /* only normal component */
    case RIGLNK: SET2 (con->R, 0.0); A->ndual += 2; break; /* normal component and multiplier */
    }
  }
  A->b = newvector (A->ndual);
  A->x = newvector (A->ndual);
  A->z = newvector (A->ndual);
  ERRMEM (A->r = MEM_CALLOC (sizeof (double [m])));
  ERRMEM (A->u = MEM_CALLOC (sizeof (double [nprimal])));
  ERRMEM (A->W = MEM_CALLOC (A->ndual * sizeof (double [9])));

  /* update previous and free local velocities V, B */
  if (dynamic) update_const_local_velocities (dom);

  /* set up constraint pointers and initialize x */
  for (con = dom->con, b = A->b->x, x = A->x->x; con; con = con->next)
  {
    n = x - A->x->x;
    R = con->R;

    CON_INDEX_SET (con, n);
    CON_RHS_SET (con, b);
    CON_X_SET (con, x);

    switch (con->kind)
    {
    case CONTACT:
    case FIXPNT:
    case GLUE: COPY (R, x); b += 3; x += 3; break;
    case FIXDIR:
    case VELODIR: x [0] = R [2]; b += 1; x += 1; break;
    case RIGLNK: x [0] = R [2]; b += 2; x += 2; break;
    }
  }

#if MPI
  A->pattern = update_local_velocities (A);
#endif

  /* since W and B are new, we update U with old reactions: U = W R + B */
  for (bod = dom->bod, r = A->r, u = A->u; bod; u += bod->dofs, bod = bod->next)
  {
    BODY_Reac (bod, r);
    BODY_Invvec (step, bod, r, 0.0, u); /* would be velocity (***) */
  }

  update_local_velocities (A); /* U = W R, based on (***) */

  for (con = dom->con; con; con = con->next)
  {
    B = con->dia->B;
    U = con->U;

    ACC (B, U); /* U = W R + B */
  }

  /* approximate diagonal W blocks */
  approximate_diagobal_blocks (dom->con, A->W);

  return A;
}

/* update linear system */
static void update_system (BSS_DATA *A)
{
  short dynamic = A->dom->dynamic;

  for (CON *con = A->dom->con; con; con = con->next)
  {
    double *U = con->U,
	   *V = con->dia->V,
	   *b = CON_RHS (con);

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
    }
    break;
    case FIXDIR:
    {
      if (dynamic) b [0] = -V[2]-U[2];
      else b [0] = -U[2];
    }
    break;
    case VELODIR:
    {
      b [0] = VELODIR(con->Z)-U[2];
    }
    break;
    case RIGLNK:
    {
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
    }
    break;
    case CONTACT:
    {
      double *X = CONTACT_X (con),
	     *Y = CONTACT_Y (con);

      VIC_Linearize (con, SMOOTHING_EPSILON, b, X, Y);

      SCALE (b, -1.0);
    }
    break;
    }
  }
}

/* solve linear system */
static void linear_solve (BSS_DATA *A, double resdec, int maxiter)
{
  hypre_FlexGMRESFunctions *gmres_functions;
  void *gmres_vdata;
  double abstol;

#if 1
  A->delta = A->resnorm / A->xnorm; /* sqrt (L-curve) */
#else
  A->delta = 0; /* FIXME */
#endif

  abstol = resdec * A->resnorm;

#if 1
  if (abstol == 0.0) /* initially */
  {
    abstol = ABSTOL * sqrt (InnerProd (A->b, A->b));
    if (abstol == 0.0) abstol = ABSTOL;
  }
#else
  abstol = 1E-10;
  maxiter = 100;
#endif

#if 0
  printf ("A->resnorm = %g, A->xnorm %g, A->delta = %g, abstol = %g\n", A->resnorm, A->xnorm, A->delta, abstol);
#endif

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
  hypre_FlexGMRESSetup (gmres_vdata, A, A->b, A->z);
  hypre_FlexGMRESSolve (gmres_vdata, A, A->b, A->z);
  hypre_FlexGMRESGetNumIterations (gmres_vdata , &A->iters);
  hypre_FlexGMRESGetFinalRelativeResidualNorm (gmres_vdata, &A->resnorm);
  hypre_FlexGMRESDestroy (gmres_vdata);

  A->resnorm *= sqrt (InnerProd (A->b, A->b)); /* relative = |Az-b|/|b| */ /* FIXME: compute yourself like in lin.c (seems to matter) */

#if 1
  {
    VECTOR *r = CreateVector (A->b);
    CopyVector (A->b, r);
    A->delta = 0;
    Matvec (NULL, -1, A, A->z, 1.0, r);
    double resnorm = sqrt (InnerProd (r, r));
    DestroyVector (r);
    //printf ("resnorm from GMRES = %g and from hand computation = %g\n", A->resnorm, resnorm);
    A->resnorm = resnorm;
  }
#endif

  A->xnorm = sqrt (InnerProd (A->z, A->z));
}

/* update solution */
static void update_solution (BSS_DATA *A)
{
  double *B, *U, *R, *x, *r, *u, step;
  DOM *dom = A->dom;
  BODY *bod;
  CON *con;

  step = dom->step;

  Axpy (1.0, A->z, A->x); /* x = x + dx */

  for (con = dom->con; con; con = con->next)
  {
    x = CON_X (con);
    R = con->R;

    switch (con->kind)
    {
    case CONTACT:
    {
      VIC_Project (con->mat.base->friction, x, R); /* project onto the friction cone */
    }
    break;
    case FIXPNT:
    case GLUE:
      COPY (x, R);
      break;
    case FIXDIR:
    case VELODIR:
    case RIGLNK:
      R [2] = x [0];
      break;
    }
  }

  /* FIXME: make sure that the last update from Matvec cannot be used instead */

  /* update U = W R + B */
  for (bod = dom->bod, r = A->r, u = A->u; bod; u += bod->dofs, bod = bod->next)
  {
    BODY_Reac (bod, r);
    BODY_Invvec (step, bod, r, 0.0, u); /* would be velocity (***) */
  }

  update_local_velocities (A); /* U = W R, based on (***) */

  for (con = dom->con; con; con = con->next)
  {
    B = con->dia->B;
    U = con->U;

    ACC (B, U); /* U = W R + B */
  }
}

/* compute merit function value */
static double merit_function (BSS_DATA *A)
{
  double G [3], *U, *V, dot;
  DOM *dom = A->dom;
  short dynamic = dom->dynamic;
  CON *con;

  for (dot = 0, con = dom->con; con; con = con->next)
  {
    V = con->dia->V;
    U = con->U;

    switch (con->kind)
    {
    case FIXPNT:
    case GLUE:
      if (dynamic) { ADD (U, V, G); }
      else { COPY (U, G); }
      dot += DOT (G, G);
      break;
    case FIXDIR:
      if (dynamic) G [0] = U[2] - V[2];
      else G [0] = U[2];
      dot += G[0]*G[0];
      break;
    case VELODIR:
      G [0] = VELODIR(con->Z) - U[2];
      dot += G[0]*G[0];
      break;
    case RIGLNK:
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
      break;
    case CONTACT:
      VIC_Linearize (con, 0, G, NULL, NULL);
      dot += DOT (G, G);
      break;
    }
  }

#if MPI
  /* TODO: sum up dots in parallel */
#endif

  return sqrt (dot);
}

/* destroy BSS data */
static void destroy_data (BSS_DATA *A)
{
  DestroyVector (A->b);
  DestroyVector (A->x);
  DestroyVector (A->z);
  free (A->r);
  free (A->u);
  free (A->W);

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
