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
#include "mrf.h"
#include "ext/krylov/krylov.h"

#if MPI
#include "com.h"
#include "pck.h"
#endif

#define SMOOTHING_EPSILON	1E-10  /* TODO */
#define ABSTOL			1E-15  /* TODO */

typedef struct bss_con_data BSS_CON_DATA;
typedef struct bss_data BSS_DATA;
typedef struct vector VECTOR;

struct bss_con_data
{
  MX *mH; /* master H */
  int mi, /* shift to mH in global velocity space */
     *mj; /* index mapping for mH converted by csc_to_dense () */

  MX *sH;
  int si,
     *sj;

  int n; /* constraint index */

  double *R, /* con->R */
	 *U, /* con->U */
         *V, /* con->dia->V */
	 *B, /* con->dia->B */
	 *W, /* con->dia->W */
	 *A, /* con->dia->A */
	  X [9], /* U-linearization */
	  Y [9], /* R-linearization */
	  T [9], /* diagonal preconditioner */
	  S [3]; /* S = W R + B - U */

  short kind; /* con->kind */

  CON *con; /* constraint */
};

struct bss_data
{
  int nprimal, /* SUM [bod in dom->bod] { bod->dofs } */ 
      ndual, /* 3 * ndat */
      ndat; /* number of active constraints */

  BSS_CON_DATA *dat; /* active constraints data */

  VECTOR *b, /* right hand side (of size ndual) */
         *x, /* reactions (-||-) */
	 *y, /* increments of velocities (-||-) */
	 *z; /* increments of reactions (-||-) */

  double *r, /* global reaction (of size nprimal) */
	 *u, /* global velocity (-||-) */
	 *a, /* auxiliary vector (of size MAX [bod in bod->dom] { bod->dofs }) */
          resnorm, /* |A z - b| */
	  znorm, /* |z| */
	  delta; /* regularisation = |A z - b| / |z| */

  int iters; /* linear solver iterations */

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

/* convert a sparse matrix into a dense one */
static MX* csc_to_dense (MX *a, int **map)
{
  int n, m, *i, *k, *p, *q, *f, *g;
  double *ax, *bx;
  MX *b;

  ERRMEM (f = MEM_CALLOC (sizeof (int [a->n]))); /* nonempty column flags */
  
  for (n = 0, m = a->m, p = a->p, q = p + a->n, ax = a->x, g = f; p < q; p ++, g ++)
  {
    if (*p < *(p+1)) /* if column is nonempty */
    {
      for (i = &a->i [*p], k = &a->i [*(p+1)]; i < k; i ++, ax ++)
      {
	if (*ax != 0.0) (*g) ++; /* and it really is */
      }

      if (*g > 0) n ++; /* account for it */
    }
  }

  b = MX_Create (MXDENSE, m, n, NULL, NULL);
  ERRMEM ((*map) = MEM_CALLOC (sizeof (int [n])));

  for (n = 0, p = a->p, bx = b->x, g = f; p < q; p ++, g ++)
  {
    if (*g > 0)  /* for each nonempty column */
    {
      for (i = &a->i [*p], k = &a->i [*(p+1)], ax = &a->x [*p]; i < k; i ++, ax ++)
      {
	bx [*i] = *ax; /* write it into dense storage */
      }
      
      (*map) [n] = p - a->p; /* map dense blocks to body-dofs indices */
      bx += m, n ++; /* next dense column */
    }
  }

  MX_Destroy (a);
  free (f);

  return b;
}

/* scatter and add to sparse vector */
inline static void scatter (double *a, int *j, int n, double *q)
{
  if (j)
  {
    for (int *k = j + n; j < k; j ++, a ++) q [*j] += *a;
  }
}

/* gather values of sparse vector */
inline static double* gather (double *q, int *j, int n, double *a)
{
  if (j)
  {
    for (int *k = j + n; j < k; j ++, a ++) *a = q [*j];

    return (a-n);
  }

  return q;
}

/* y = H' x */
static void H_trans_vector (double *a, BSS_CON_DATA *dat, int n, double *x, double *y)
{
  for (; n > 0; dat ++, n --)
  {
    double *p = &x [dat->n], *q, *r, c;

    if (dat->mH)
    {
      q = &y [dat->mi];

      if (dat->mj) c = 0.0, r = a;
      else c = 1.0, r = q;

      MX_Matvec (1.0, MX_Tran (dat->mH), p, c, r);

      scatter (r, dat->mj, dat->mH->n, q);
    }

    if (dat->sH)
    {
      q = &y [dat->si];

      if (dat->sj) c = 0.0, r = a;
      else c = 1.0, r = q;

      MX_Matvec (1.0, MX_Tran (dat->sH), p, c, r);

      scatter (r, dat->sj, dat->sH->n, q);
    }
  }
}

/* y = H' x */
static void H_times_vector (double *a, BSS_CON_DATA *dat, int n, double *x, double *y)
{
  for (; n > 0; dat ++, n --)
  {
    double *p, *q = &y [dat->n], *r;

    if (dat->mH)
    {
      p = &x [dat->mi];

      r = gather (p, dat->mj, dat->mH->n, a);

      MX_Matvec (1.0, dat->mH, r, 1.0, q);
    }

    if (dat->sH)
    {
      p = &x [dat->si];

      r = gather (p, dat->sj, dat->sH->n, a);

      MX_Matvec (1.0, dat->sH, r, 1.0, q);
    }
  }
}

/* y = alpha W x */
static void W_times_vector (BSS_DATA *A, double *x, double *y)
{
  double *r = A->r, *u = A->u, *a = A->a, step;
  BSS_CON_DATA *dat = A->dat;
  DOM *dom = A->dom;
  int n = A->ndat;
  BODY *bod;

  step = dom->step;
  SETN (r, A->nprimal, 0.0);
  H_trans_vector (a, dat, n, x, r);

#if MPI
  /* TODO: send x to external reactions
   * TODO: sum up external reaction contributions into r */
#endif

  for (bod = dom->bod; bod; r += bod->dofs, u += bod->dofs, bod = bod->next)
  {
    BODY_Invvec (step, bod, r, 0.0, u);
  }

  SETN (y, A->ndual, 0.0);
  H_times_vector (a, dat, n, A->u, y);

#if MPI
  /* TODO: compute local velocities needed by constraints established by child bodies
   * TODO: and send them from parent to local constraints => into y */
#endif

}

/* GMSS interface start */
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
  double Z [3], *Y, *X, *R, *U, *Q, *u, *q, *r;
  BSS_CON_DATA *dat = A->dat;
  int n = A->ndat;

  W_times_vector (A, x->x, A->y->x); /* dU = h W dR */

  for (u = A->y->x, r = x->x, q = y->x; n > 0; dat ++, n --)
  {
    U = &u [dat->n];
    R = &r [dat->n];
    Q = &q [dat->n];

    if (dat->kind == VELODIR || dat->kind == FIXDIR)
    {
      Z [0] = R [0];
      Z [1] = R [1];
      Z [2] = U [2];

      U = Z;
    }
    else if (dat->kind == CONTACT)
    {
      X = dat->X;
      Y = dat->Y;

      NVMUL (X, U, Z);
      NVADDMUL (Z, Y, R, Z);  /* Z = X dU + Y dR */

      U = Z;
    }
    else if (dat->kind == RIGLNK)
    {
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
    }

    SCALE (Q, beta);
    ADDMUL (Q, alpha, U, Q)
  }

  if (A->delta > 0.0) Axpy (alpha * A->delta, x, y);

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
  BSS_CON_DATA *dat, *end;
  double *T, *p, *q;

  for (dat = A->dat, end = dat + A->ndat; dat != end; dat ++)
  {
    p = &b->x[dat->n];
    q = &x->x[dat->n];
    T = dat->T;
    NVMUL (T, p, q);
  }

  return 0;
}
/* GMSS interface end */

/* update previous and free local velocities V, B */
static void update_V_and_B (DOM *dom)
{
  double X [6], *V, *B;
  CON *con;

  for (con = dom->con; con; con = con->next)
  {
    V = con->dia->V;
    B = con->dia->B;

    SET (V, 0);
    SET (B, 0);

#if MPI
    if (con->master->flags & BODY_PARENT) /* local parent */
#endif
    {
      BODY_Local_Velo (con->master, mshp(con), mgobj(con), con->mpnt, con->base, X, X+3);
      ADD (V, X, V);
      ADD (B, X+3, B);
    }

    if (con->slave)
    {
#if MPI
      if (con->slave->flags & BODY_PARENT) /* local slave */
#endif
      {
        BODY_Local_Velo (con->slave, sshp(con), sgobj(con), con->spnt, con->base, X, X+3);
	SUB (V, X, V); /* relative = master - slave */
	SUB (B, X+3, B);
      }
    }
  }

#if MPI
  /* in parallel, parent bodies calculate local velocities of their external constraints and send them to the constraint parents;
   * upon arrival these are added or subtracted from the constraint parents B and V members, depending on the master/slave relation */
  int nsend, nrecv, *isize, *dsize, i, *j, *k;
  COMDATA *send, *recv, *ptr;
  double *D;
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
	ADD (V, D, V);
	ADD (B, D+3, B);
      }
      else /* slave */
      {
	SUB (V, D, V);
	SUB (B, D+3, B);
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
}

/* create constraints data */
static BSS_CON_DATA *create_constraints_data (DOM *dom, int *ndat, VECTOR **v, double *free_energy)
{
  BSS_CON_DATA *out, *dat;
  double step, *x;
  short dynamic;
  BODY *bod;
  CON *con;
  MAP *map;
  int n;

  dynamic = dom->dynamic;
  (*free_energy) = 0.0;

  ERRMEM (out = MEM_CALLOC (dom->ncon * sizeof (BSS_CON_DATA)));
  (*v) = newvector (dom->ncon * 3);
  step = dom->step;
  x = (*v)->x;

  for (bod = dom->bod, map = NULL, n = 0; bod; n += bod->dofs, bod = bod->next)
  {
    MAP_Insert (NULL, &map, bod, (void*) (long) n, NULL);
  }

  for (con = dom->con, dat = out, n = 0; con; con = con->next)
  {
    if (dynamic && con->kind == CONTACT && con->gap > 0.0) continue; /* skip open dynamic contacts */

    DIAB *dia = con->dia;
    BODY *m = con->master,
	 *s = con->slave;
    void *mgobj = mgobj(con),
	 *sgobj;
    SHAPE *mshp = mshp(con),
	  *sshp;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base,
	   *B = dia->B,
	    X [3];
    MX_DENSE_PTR (W, 3, 3, dia->W);
    MX_DENSE_PTR (A, 3, 3, dia->A);
    short up = dia->W[8] == 0.0; /* needed only as a preconditioner => update once */
    MX *prod;

    if (s)
    {
      sgobj = sgobj(con);
      sshp = sshp(con);
    }

#if MPI
    /* TODO: if master is a child body then dat->mH is NULL */
#endif

    if (m->kind != OBS)
    {
      dat->mH = BODY_Gen_To_Loc_Operator (m, mshp, mgobj, mpnt, base);
      dat->mi = (int) (long) MAP_Find (map, m, NULL);

      if (up)
      {
	prod = MX_Matmat (1.0, m->inverse, MX_Tran (dat->mH), 0.0, NULL);
	MX_Matmat (step, dat->mH, prod, 0.0, &W); /* H * inv (M) * H^T */
	MX_Destroy (prod);
      }

      if (m->kind == FEM)
      {
	dat->mH = csc_to_dense (dat->mH, &dat->mj);
      }
    }

#if MPI
    /* TODO: if slave is a child body then dat->sH is NULL */
#endif

    if (s && s->kind != OBS)
    {
      dat->sH = BODY_Gen_To_Loc_Operator (s, sshp, sgobj, spnt, base);
      dat->si = (int) (long) MAP_Find (map, s, NULL);
      MX_Scale (dat->sH, -1.0);

      if (up)
      {
	prod = MX_Matmat (1.0, s->inverse, MX_Tran (dat->sH), 0.0, NULL);
	MX_Matmat (step, dat->sH, prod, 1.0, &W); /* H * inv (M) * H^T */
	MX_Destroy (prod);
      }

      if (s->kind == FEM)
      {
	dat->sH = csc_to_dense (dat->sH, &dat->sj);
      }
    }

    if (up)
    {
      MX_Copy (&W, &A);
      MX_Inverse (&A, &A);
    }

    NVMUL (A.x, B, X);
    (*free_energy) += DOT (X, B); /* sum up free energy */

    /* add up prescribed velocity contribution */
    if (con->kind == VELODIR) (*free_energy) += A.x[8] * VELODIR(con->Z) * VELODIR(con->Z);

    dat->n = n;
    dat->R = con->R;
    dat->U = con->U;
    dat->V = dia->V;
    dat->B = dia->B;
    dat->W = dia->W;
    dat->A = dia->A;
    dat->kind = con->kind;
    dat->con = con;

    COPY (dat->R, x); /* initialize with previous solution */

    dat ++;
    n += 3;
    x += 3;
  }

  (*free_energy) *= 0.5;
  (*ndat) = n / 3;
  (*v)->n = n;

  MAP_Free (NULL, &map);

  return out;
}

static void destroy_constraints_data (BSS_CON_DATA *dat, int n)
{
  BSS_CON_DATA *ptr = dat;

  for (; n > 0; dat ++, n --)
  {
    if (dat->mH) MX_Destroy (dat->mH);
    if (dat->sH) MX_Destroy (dat->sH);
    free (dat->mj);
    free (dat->sj);
  }

  free (ptr);
}

/* create BSS data */
static BSS_DATA *create_data (DOM *dom)
{
  BSS_CON_DATA *dat, *end;
  BSS_DATA *A;
  BODY *bod;
  int m, n;

  update_V_and_B (dom);

  ERRMEM (A = MEM_CALLOC (sizeof (BSS_DATA)));
  A->znorm = 1.0;
  A->dom = dom;
  for (bod = dom->bod, m = 0; bod; bod = bod->next)
  {
    n = bod->dofs;
    if (n > m) m = n;
    A->nprimal += n;
  }
  A->dat = create_constraints_data (dom, &A->ndat, &A->x, &dom->ldy->free_energy); /* A->x initialized with con->R */
  A->ndual = A->ndat * 3;
  A->b = newvector (A->ndual);
  A->y = newvector (A->ndual);
  A->z = newvector (A->ndual);
  ERRMEM (A->r = MEM_CALLOC (sizeof (double [A->nprimal])));
  ERRMEM (A->u = MEM_CALLOC (sizeof (double [A->nprimal])));
  ERRMEM (A->a = MEM_CALLOC (sizeof (double [m])));


  /* initial residual S = WR + B - U */
  double *Ayx = A->y->x;
  W_times_vector (A, A->x->x, Ayx);
  for (dat = A->dat, end = dat + A->ndat; dat != end; dat ++)
  {
    double *WR = &Ayx [dat->n],
	   *B  = dat->B,
	   *U  = dat->U,
	   *S  = dat->S;

    S [0] = WR[0] + B[0] - U[0];
    S [1] = WR[1] + B[1] - U[1];
    S [2] = WR[2] + B[2] - U[2];
  }

  return A;
}

/* update linear system */
static void update_system (BSS_DATA *A)
{
  double *Abx = A->b->x, delta = A->delta;
  DOM *dom = A->dom;
  short dynamic = dom->dynamic;
  BSS_CON_DATA *dat, *end;

  for (dat = A->dat, end = dat + A->ndat; dat != end; dat ++)
  {
    double *U = dat->U,
	   *R = dat->R,
	   *V = dat->V,
	   *W = dat->W,
	   *T = dat->T,
	   *S = dat->S,
	   *b = &Abx [dat->n];

    switch (dat->kind)
    {
    case FIXPNT:
    case GLUE:
    {
      if (dynamic)
      {
	b [0] = -V[0]-U[0]-S[0];
	b [1] = -V[1]-U[1]-S[1];
	b [2] = -V[2]-U[2]-S[2];
      }
      else
      {
	b [0] = -U[0]-S[0];
	b [1] = -U[1]-S[1];
	b [2] = -U[2]-S[2];
      }

      NNCOPY (W, T);
    }
    break;
    case FIXDIR:
    {
      b [0] = -R[0];
      b [1] = -R[1];
      if (dynamic) b [2] = -V[2]-U[2]-S[2];
      else b [2] = -U[2]-S[2];

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
      b [2] = VELODIR(dat->con->Z)-U[2]-S[2];

      T [1] = T [3] = T [6] = T [7] = 0.0;
      T [0] = T [4] = 1.0;
      T [2] = W [2];
      T [5] = W [5];
      T [8] = W [8];
    }
    break;
    case RIGLNK:
    {
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
    }
    break;
    case CONTACT:
    {
      double *X = dat->X, *Y = dat->Y, G [3];

      VIC_Linearize (dat->con, SMOOTHING_EPSILON, G, X, Y);

      NVADDMUL (G, X, S, b);
      SCALE (b, -1.0);

      NNMUL (X, W, T);
      NNADD (T, Y, T);
    }
    break;
    }

    T [0] += delta;
    T [4] += delta;
    T [8] += delta;

    MX_DENSE_PTR (P, 3, 3, T);
    MX_Inverse (&P, &P);
  }
}

/* solve linear system */
static void linear_solve (BSS_DATA *A, double resdec, int maxiter)
{
  hypre_FlexGMRESFunctions *gmres_functions;
  void *gmres_vdata;
  double abstol;
  VECTOR *r;
  int ret;

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

  hypre_error_flag = 0;
  hypre_FlexGMRESSetTol (gmres_vdata, 0.0);
  hypre_FlexGMRESSetMinIter (gmres_vdata, 1);
  hypre_FlexGMRESSetMaxIter (gmres_vdata, maxiter);
  hypre_FlexGMRESSetAbsoluteTol (gmres_vdata, abstol);
  hypre_FlexGMRESSetup (gmres_vdata, A, A->b, A->z);
  ret = hypre_FlexGMRESSolve (gmres_vdata, A, A->b, A->z);
  hypre_FlexGMRESGetNumIterations (gmres_vdata , &A->iters);
  hypre_FlexGMRESGetResidual (gmres_vdata, (void**) &r);

  if (ret & HYPRE_ERROR_CONV) /* failed to converge => invalid residual */
  {
    A->delta = 0;
    Matvec (NULL, -1.0, A, A->z, 1.0, A->b); /* residual without regularisation (delta = 0) */
    A->resnorm = sqrt (InnerProd (A->b, A->b));
  }
  else /* converged => valid residual */
  {
    Axpy (A->delta, A->z, r); /* residua without regularisation */
    A->resnorm = sqrt (InnerProd (r, r));
  }
  A->znorm = sqrt (InnerProd (A->z, A->z));
  A->delta = A->resnorm / A->znorm; /* sqrt (L-curve) */

  hypre_FlexGMRESDestroy (gmres_vdata);
}

/* update solution */
static void update_solution (BSS_DATA *A)
{
  double *R, *x, *z, *Azx = A->z->x, *Ayx = A->y->x, *Axx = A->x->x;
  BSS_CON_DATA *dat, *end;

  for (dat = A->dat, end = dat + A->ndat; dat != end; dat ++)
  {
    x = &Axx [dat->n];
    z = &Azx [dat->n];
    R = dat->R;

    if (dat->kind == CONTACT)
    {
      double S [3];
      ADD (x, z, S); /* S = R + DR */
      VIC_Project (dat->con->mat.base->friction, S, S); /* project onto the friction cone */
      SUB (S, x, z); /* DR = proj (R + DR) - R */
    }

    ACC (z, x); /* R = R + DR */
    COPY (x, R);
  }

  /* DU = W DR
   * U1 = U0 + DU + S
   * S  = W (R + DR) + B - U1
   *    = S0 + U0 + DU - U1 */
  W_times_vector (A, Azx, Ayx);
  for (dat = A->dat; dat != end; dat ++)
  {
    double *DU = &Ayx [dat->n],
	   *U  = dat->U,
	   *S  = dat->S,
	    U0 [3];

    COPY (U, U0);

    U [0] = (U [0] + DU [0]) + S [0];
    U [1] = (U [1] + DU [1]) + S [1];
    U [2] = (U [2] + DU [2]) + S [2];

    S [0] = ((S [0] + U0 [0]) + DU [0]) - U [0];
    S [1] = ((S [1] + U0 [1]) + DU [1]) - U [1];
    S [2] = ((S [2] + U0 [2]) + DU [2]) - U [2];
  }
}

/* destroy BSS data */
static void destroy_data (BSS_DATA *A)
{
  destroy_constraints_data (A->dat, A->ndat);
  DestroyVector (A->b);
  DestroyVector (A->x);
  DestroyVector (A->y);
  DestroyVector (A->z);
  free (A->r);
  free (A->u);
  free (A->a);
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

    *merit = MERIT_Function (ldy, 0);

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
