/*
 * nts.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * projected quasi-Newton solver
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

#include "bgs.h"
#include "nts.h"
#include "dom.h"
#include "fem.h"
#include "alg.h"
#include "bla.h"
#include "lap.h"
#include "err.h"
#include "scf.h"
#include "mrf.h"
#include "lis.h"
#include "ext/krylov/krylov.h"

#if MPI
#include "tag.h"
#endif

typedef struct con_data CON_DATA;
typedef struct private PRIVATE;
typedef struct vector VECTOR;

struct con_data
{
  MX *mH; /* master H */
  int mi, /* shift to mH in global velocity space */
     *mj; /* index mapping for mH converted by csc_to_dense () */

  MX *sH;
  int si,
     *sj;

  CON *con; /* constraint */

  double X [9], /* U-linearization */
	 Y [9], /* R-linearization */
	 T [9]; /* diagonal preconditioner */

  double R0 [3], /* initial reaction */
         RC [3]; /* current reaction */
};

struct private
{
  NEWTON *ns; /* solver data */

  DOM *dom; /* domain */

  CON_DATA *dat, *end; /* constraints data */

  MAP *bod; /* maps bodies to shifts in u and r */

  int ndofs; /* SUM { bod[]->dofs } */ 

  double *u, /* body space velocity */
         *r, /* body space reaction */
	 *a; /* auxiliary vector */

  int matvec; /* matrix vector products */

  VECTOR *dr, /* reactions increment */
	 *rhs; /* right hand side of linearization */
#if MPI
  SET *inner, *boundary;
  COMDATA *send, *recv;
  int nsend, nrecv;
  void *pattern; /* non-blocking communication pattern */
#endif
};

struct vector
{
  double *x;

  int n;
};

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
      
      (*map) [n] = p - a->p; /* map dense blocks to body->dofs indices */
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

/* y += H' R */
static void H_trans_R (double *a, CON_DATA *dat, CON_DATA *end, double *y)
{
  for (; dat != end; dat ++)
  {
    double *p = dat->con->R, *q, *r, c;

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

/* U += H u */
static void H_times_u (double *a, CON_DATA *dat, CON_DATA *end, double *u)
{
  for (; dat != end; dat ++)
  {
#if MPI
    if (!dat->con->dia) break; /* skip external */
#endif
    double *p, *q = dat->con->U, *r;

    if (dat->mH)
    {
      p = &u [dat->mi];

      r = gather (p, dat->mj, dat->mH->n, a);

      MX_Matvec (1.0, dat->mH, r, 1.0, q);
    }

    if (dat->sH)
    {
      p = &u [dat->si];

      r = gather (p, dat->sj, dat->sH->n, a);

      MX_Matvec (1.0, dat->sH, r, 1.0, q);
    }
  }
}

#if MPI
static void update_external_reactions (PRIVATE *A)
{
  DOM *dom = A->dom;
  int i, *j, *k;
  COMDATA *ptr;
  double *R;
  CON *con;

  for (i = 0; i < A->nrecv; i ++)
  {
    ptr = &A->recv [i];
    for (j = ptr->i, k = j + ptr->ints, R = ptr->d; j < k; j ++, R += 3)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*j), NULL), "Invalid con id: %d", *j);
      COPY (R, con->R);
    }
  }
}
#endif

/* U = W R + B */
static void U_WR_B (PRIVATE *A, short zero_B)
{
#if MPI
  if (A->pattern == NULL) /* initialize communication pattern */
  {
    CON_DATA *dat;
    COMDATA *ptr;
    SET *item;
    CON *con;

    for (dat = A->dat, A->nsend = 0, A->inner = A->boundary = NULL; dat != A->end; dat ++) 
    {
      con = dat->con;
      if (con->ext)
      {
        A->nsend += SET_Size (con->ext);
	SET_Insert (NULL, &A->boundary, con, NULL); /* boundary constraints */
      }
      else SET_Insert (NULL, &A->inner, con, NULL); /* inner constraints */
    }
    ERRMEM (A->send = MEM_CALLOC (sizeof (COMDATA [A->nsend])));

    for (dat = A->dat, ptr = A->send; dat != A->end; dat ++)
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

    A->pattern = COM_Pattern (MPI_COMM_WORLD, TAG_NEWTON, A->send, A->nsend, &A->recv, &A->nrecv);
  }
#endif

  if (A->ns->locdyn == LOCDYN_ON)
  {
#if MPI
    double *W, *R, *U, *B;
    SET *item;
    DIAB *dia;
    OFFB *blk;
    CON *con;

    COM_Send (A->pattern);

    for (item = SET_First (A->inner); item; item = SET_Next (item)) /* process inner constraints while sending */
    {
      con = item->data;
      dia = con->dia;
      W = dia->W;
      R = con->R;
      U = con->U;
      if (zero_B)
      {
	NVMUL (W, R, U);
      }
      else
      {
	B = dia->B;
	NVADDMUL (B, W, R, U);
      }
      for (blk = dia->adj; blk; blk = blk->n)
      {
	R = blk->dia->R;
	W = blk->W;
	NVADDMUL (U, W, R, U);
      }
      ASSERT_DEBUG (dia->adjext == NULL, "Inconsistent inner constraint");
    }

    COM_Recv (A->pattern);
    update_external_reactions (A);

    for (item = SET_First (A->boundary); item; item = SET_Next (item)) /* process boundary constraints */
    {
      con = item->data;
      dia = con->dia;
      W = dia->W;
      R = con->R;
      U = con->U;
      if (zero_B)
      {
	NVMUL (W, R, U);
      }
      else
      {
	B = dia->B;
	NVADDMUL (B, W, R, U);
      }
      for (blk = dia->adj; blk; blk = blk->n)
      {
	R = blk->dia->R;
	W = blk->W;
	NVADDMUL (U, W, R, U);
      }
      for (blk = dia->adjext; blk; blk = blk->n)
      {
	R = CON (blk->dia)->R;
	W = blk->W;
	NVADDMUL (U, W, R, U);
      }
    }
#else
    double *W, *R, *U, *B;
    CON_DATA *dat;
    DIAB *dia;
    OFFB *blk;
    CON *con;

    for (dat = A->dat; dat != A->end; dat ++)
    {
      con = dat->con;
      dia = con->dia;
      W = dia->W;
      R = con->R;
      U = con->U;
      if (zero_B)
      {
	NVMUL (W, R, U);
      }
      else
      {
	B = dia->B;
	NVADDMUL (B, W, R, U);
      }
      for (blk = dia->adj; blk; blk = blk->n)
      {
	R = blk->dia->R;
	W = blk->W;
	NVADDMUL (U, W, R, U);
      }
    }
#endif
  }
  else
  {
    double *B, *U, step;
    CON_DATA *dat;
    MAP *item;
    DIAB *dia;
    CON *con;
    int n;

#if MPI
    COM_Repeat (A->pattern);
    update_external_reactions (A);
#endif

    step = A->dom->step;
    SETN (A->r, A->ndofs, 0.0);
    H_trans_R (A->a, A->dat, A->end, A->r);

    for (item = MAP_First (A->bod); item; item = MAP_Next (item))
    {
      n = (int) (long) item->data;
      BODY_Invvec (step, item->key, &A->r [n], 0.0, &A->u [n]);
    }

    for (dat = A->dat; dat != A->end; dat ++)
    {
      con = dat->con;
      dia = con->dia;
#if MPI
      if (!dia) break; /* skip external */
#endif
      if (zero_B)
      {
	U = con->U;
	SET (U, 0.0);
      }
      else
      {
	B = dia->B;
	U = con->U;
	COPY (B, U); /* U = B */
      }
    }

    H_times_u (A->a, A->dat, A->end, A->u); /* U += H u */
  }
}

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
  double Z [3], *U, *X, *Y, *R, *Q;
  CON_DATA *dat;
  CON *con;

  for (dat = A->dat, R = x->x; dat != A->end; dat ++, R += 3)
  {
    con = dat->con;
#if MPI
    if (!con->dia) break; /* skip external */
#endif
    Q = con->R;
    COPY (R, Q);
  }

  U_WR_B (A, 1);

  for (dat = A->dat, R = x->x, Q = y->x; dat != A->end; dat ++, R += 3, Q += 3)
  {
    con = dat->con;
#if MPI
    if (!con->dia) break; /* skip external */
#endif
    U = con->U;

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

  Axpy (alpha * A->ns->delta, x, y);

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

static int Precond (void *vdata, PRIVATE *A, VECTOR *b, VECTOR *x)
{
  double *T, *Q, *R;
  CON_DATA *dat;

  for (dat = A->dat, R = x->x, Q = b->x; dat != A->end; dat ++, R += 3, Q += 3)
  {
#if MPI
    if (!dat->con->dia) break; /* skip external */
#endif
    T = dat->T;
    NVMUL (T, Q, R);
  }

  return 0;
}
/* GMRES interface end */

/* create constraints data for body-space mode */
static int body_space_constraints_data (DOM *dom, PRIVATE *A)
{
  LOCDYN *ldy = dom->ldy;
  MAP *item, *fem;
  int ndat, ret;
  CON_DATA *dat;
  double step;
  CON *con;
#if MPI
  MAP *jtem;
#endif

  ret = 0;
  fem = NULL;
  ldy->free_energy = 0.0;
  ndat = dom->ncon;
#if MPI
  ndat += MAP_Size (dom->conext);
#endif
  step = dom->step;
  ERRMEM (A->dat = MEM_CALLOC (ndat * sizeof (CON_DATA)));

  /* internal constraints */
  for (con = dom->con, dat = A->dat; con; con = con->next)
  {
    DIAB *dia = con->dia;
    BODY *m = con->master,
	 *s = con->slave;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base;
    MX_DENSE_PTR (W, 3, 3, dia->W);
    MX *prod, *inv;

    if (m->kind == FEM && MAP_Find (fem, m, NULL) == NULL)
    {
      MAP_Insert (NULL, &fem, m, FEM_Approx_Inverse (m), NULL); /* map approximate inverses */
    }

    if (MAP_Find (A->bod, m, NULL) == NULL)
    {
      MAP_Insert (NULL, &A->bod, m, (void*) (long) A->ndofs, NULL); /* map dofs shift */
      if (m->dofs > ret) ret = m->dofs;
      A->ndofs += m->dofs;
    }

    if (s)
    {
      if (s->kind == FEM && MAP_Find (fem, s, NULL) == NULL)
      {
	MAP_Insert (NULL, &fem, s, FEM_Approx_Inverse (s), NULL); /* map approximate inverses */
      }

      if (MAP_Find (A->bod, s, NULL) == NULL)
      {
	MAP_Insert (NULL, &A->bod, s, (void*) (long) A->ndofs, NULL); /* map dofs shift */
        if (s->dofs > ret) ret = s->dofs;
	A->ndofs += s->dofs;
      }
    }

    if (!(con->kind == CONTACT && m->kind == OBS))
    {
      dat->mH = BODY_Gen_To_Loc_Operator (m, con->kind, con->msgp, mpnt, base);
      dat->mi = (int) (long) MAP_Find (A->bod, m, NULL);

      if (m->kind == FEM) inv = MAP_Find (fem, m, NULL); else inv = m->inverse;
      prod = MX_Matmat (1.0, inv, MX_Tran (dat->mH), 0.0, NULL);
      MX_Matmat (step, dat->mH, prod, 0.0, &W); /* H * inv (M) * H^T */
      MX_Destroy (prod);

      if (m->kind == FEM && m->form != REDUCED_ORDER)
      {
	dat->mH = csc_to_dense (dat->mH, &dat->mj);
      }
    }

    if (s && !(con->kind == CONTACT && s->kind == OBS))
    {
      dat->sH = BODY_Gen_To_Loc_Operator (s, con->kind, con->ssgp, spnt, base);
      dat->si = (int) (long) MAP_Find (A->bod, s, NULL);
      MX_Scale (dat->sH, -1.0);

      if (s->kind == FEM) inv = MAP_Find (fem, s, NULL); else inv = s->inverse;
      prod = MX_Matmat (1.0, inv, MX_Tran (dat->sH), 0.0, NULL);
      MX_Matmat (step, dat->sH, prod, 1.0, &W); /* H * inv (M) * H^T */
      MX_Destroy (prod);

      if (s->kind == FEM && s->form != REDUCED_ORDER)
      {
	dat->sH = csc_to_dense (dat->sH, &dat->sj);
      }
    }

    COPY (con->R, dat->R0);
    dat->con = con;
    dat ++;
  }

#if MPI
  /* external constraints */
  for (item = MAP_First (dom->conext); item; item = MAP_Next (item))
  {
    con = item->data;
    BODY *m = con->master,
	 *s = con->slave;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base;

    if (!(con->kind == CONTACT && m->kind == OBS) && (jtem = MAP_Find_Node (A->bod, m, NULL)))
    {
      dat->mH = BODY_Gen_To_Loc_Operator (m, con->kind, con->msgp, mpnt, base);
      dat->mi = (int) (long) jtem->data;

      if (m->kind == FEM && m->form != REDUCED_ORDER)
      {
	dat->mH = csc_to_dense (dat->mH, &dat->mj);
      }
    }

    if (s && !(con->kind == CONTACT && s->kind == OBS) && (jtem = MAP_Find_Node (A->bod, s, NULL)))
    {
      dat->sH = BODY_Gen_To_Loc_Operator (s, con->kind, con->ssgp, spnt, base);
      dat->si = (int) (long) jtem->data;
      MX_Scale (dat->sH, -1.0);

      if (s->kind == FEM && s->form != REDUCED_ORDER)
      {
	dat->sH = csc_to_dense (dat->sH, &dat->sj);
      }
    }

    dat->con = con;
    dat ++;
  }
#endif

  /* mark end */
  A->end = dat;

  /* process diagonal blocks */
  for (dat = A->dat; dat != A->end; dat ++)
  {
    con = dat->con;
    DIAB *dia = con->dia;
    if (!dia) break; /* skip external */
    MX_DENSE_PTR (W, 3, 3, dia->W);
    MX_DENSE_PTR (A, 3, 3, dia->A);
    double *B = dia->B, X [3];

    MX_Copy (&W, &A);
    MX_Inverse (&A, &A);
    NVMUL (A.x, B, X);
    ldy->free_energy += DOT (X, B); /* sum up free energy */

    /* add up prescribed velocity contribution */
    if (con->kind == VELODIR)
    {
      ldy->free_energy += A.x[8] * VELODIR(con->Z) * VELODIR(con->Z);
    }
  }
  ldy->free_energy *= 0.5;

  for (item = MAP_First (fem); item; item = MAP_Next (item)) MX_Destroy (item->data);
  MAP_Free (NULL, &fem);

  return ret;
}

/* create constraints data for local dynamics mode */
static void locdyn_constraints_data (DOM *dom, PRIVATE *A)
{
  short dynamic;
  CON_DATA *dat;
  CON *con;

  dynamic = dom->dynamic;
  ERRMEM (A->dat = MEM_CALLOC (dom->ncon * sizeof (CON_DATA)));

  /* create internal constraints data */
  for (con = dom->con, dat = A->dat; con; con = con->next)
  {
    COPY (con->R, dat->R0);
    dat->con = con;
    dat ++;
  }
  A->end = dat;
}

/* destroy constraints data */
static void destroy_constraints_data (CON_DATA *dat, CON_DATA *end)
{
  CON_DATA *ptr = dat;

  for (; dat != end; dat ++)
  {
    if (dat->mH) MX_Destroy (dat->mH);
    if (dat->sH) MX_Destroy (dat->sH);
    free (dat->mj);
    free (dat->sj);
  }

  free (ptr);
}

/* create private data */
static PRIVATE *create_private_data (NEWTON *ns, LOCDYN *ldy)
{
  PRIVATE *A;
  int n;

  ERRMEM (A = MEM_CALLOC (sizeof (PRIVATE)));
  A->dom = ldy->dom;
  A->ns = ns;

  if (ns->locdyn == LOCDYN_OFF)
  {
    n = body_space_constraints_data (ldy->dom, A);
    ERRMEM (A->u = MEM_CALLOC (sizeof (double [A->ndofs])));
    ERRMEM (A->r = MEM_CALLOC (sizeof (double [A->ndofs])));
    ERRMEM (A->a = MEM_CALLOC (sizeof (double [n])));
  }
  else locdyn_constraints_data (ldy->dom, A);

  A->dr = newvector  (3 * (A->end - A->dat));
  A->rhs = CreateVector (A->dr);

  U_WR_B (A, 0); /* U = W R + B */

  return A;
}

/* destroy private data */
static void destroy_private_data (PRIVATE *A)
{
  destroy_constraints_data (A->dat, A->end);
  MAP_Free (NULL, &A->bod);
  free (A->u);
  free (A->r);
  free (A->a);
  DestroyVector (A->dr);
  DestroyVector (A->rhs);

#if MPI
  SET_Free (NULL, &A->boundary);
  SET_Free (NULL, &A->inner);
  free (A->send);
  free (A->recv);
  COM_Free (A->pattern);
#endif

  free (A);
}

/* single projected quasi-Newton step */
static int solve (PRIVATE *A, short linver, int linmaxiter, double epsilon, short dynamic,
           double step, double delta, double theta, double omega, VECTOR *dr, VECTOR *rhs)
{
  double *b  = rhs->x, *DR = dr->x, gamma = 1.0 - theta;
  int ipiv [3], iters = 0;
  CON_DATA *dat;

  for (dat = A->dat; dat != A->end; dat ++, b += 3, DR += 3)
  {
    CON *con = dat->con;
    DIAB *dia = con->dia;
#if MPI
    if (!dia) break; /* skip external */
#endif
    double *U = con->U,
	   *V = con->V,
	   *R = con->R,
	   *W = dia->W,
	   *T = dat->T;

    if (linver == PQN_GMRES)
    {
      double *RC = dat->RC;
      COPY (R, RC); /* save current reaction */
    }

    switch (con->kind)
    {
    case FIXPNT:
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
      double h = step * (dynamic ? 0.5 : 1.0),
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
    case SPRING:
    {
      WARNING (0, "SPRING not supported in NEWTON_SOLVER yet!");
      ASSERT (0, ERR_NOT_IMPLEMENTED); /* TODO */
    }
    break;
    case CONTACT:
    {
      double *X = dat->X, *Y = dat->Y;

      SCF_Linearize (dat->con, U, R, -1, omega, b, X, Y);
      SCALE (b, -1.0);

      NNMUL (X, W, T);
      NNADD (T, Y, T);
    }
    break;
    }

    if (linver == PQN_DIAG)
    {
      ASSERT (lapack_dgesv (3, 1, T, 3, ipiv, b, 3) == 0, ERR_MTX_LU_FACTOR); /* diagonalized solve */
      DR [0] = gamma * DR[0] + theta * b[0]; /* theta-averaging */
      DR [1] = gamma * DR[1] + theta * b[1];
      DR [2] = gamma * DR[2] + theta * b[2];
      ACC (DR, R);

      if (con->kind == CONTACT)
      {
	double c = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area;
	SCF_Project (con->mat.base->friction, c, R, R); /* projection */
      }
    }
    else /* PQN_GMRES */
    {
      T [0] += delta;
      T [4] += delta;
      T [8] += delta;

      MX_DENSE_PTR (P, 3, 3, T);
      MX_Inverse (&P, &P); /* preconditioner */
    }
  }

  if (linver == PQN_GMRES)
  {
    hypre_FlexGMRESFunctions *gmres_functions;
    void *gmres_vdata;
    int ret;

    gmres_functions = hypre_FlexGMRESFunctionsCreate (CAlloc, Free, (int (*) (void*,int*,int*)) CommInfo,
      (void* (*) (void*))CreateVector, (void* (*) (int, void*))CreateVectorArray, (int (*) (void*))DestroyVector,
      MatvecCreate, (int (*) (void*,double,void*,void*,double,void*))Matvec, MatvecDestroy,
      (double (*) (void*,void*))InnerProd, (int (*) (void*,void*))CopyVector, (int (*) (void*))ClearVector,
      (int (*) (double,void*))ScaleVector, (int (*) (double,void*,void*))Axpy,
      PrecondSetup, (int (*) (void*,void*,void*,void*))Precond);
    gmres_vdata = hypre_FlexGMRESCreate (gmres_functions);


    double bnorm = sqrt (InnerProd (rhs, rhs));

    hypre_error_flag = 0;
    hypre_FlexGMRESSetTol (gmres_vdata, 0.0);
    hypre_FlexGMRESSetMinIter (gmres_vdata, 1);
    hypre_FlexGMRESSetMaxIter (gmres_vdata, linmaxiter);
    hypre_FlexGMRESSetAbsoluteTol (gmres_vdata, epsilon * bnorm);
    hypre_FlexGMRESSetup (gmres_vdata, A, rhs, dr);
    ret = hypre_FlexGMRESSolve (gmres_vdata, A, rhs, dr); /* GMRES solve */
    hypre_FlexGMRESGetNumIterations (gmres_vdata , &iters);
    hypre_FlexGMRESDestroy (gmres_vdata);

    for (dat = A->dat, DR = dr->x; dat != A->end; dat ++, DR += 3)
    {
      CON *con = dat->con;
#if MPI
    if (!con->dia) break; /* skip external */
#endif
      double *RC = dat->RC,
	     *R = con->R;
      ADD (DR, RC, R);

      if (con->kind == CONTACT)
      {
	double c = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area;
	SCF_Project (con->mat.base->friction, c, R, R); /* project */
      }
    }
  }

  return iters;
}

/* reset reactions to the inital values */
static void restore_initial_R (PRIVATE *A)
{
  CON_DATA *dat;
  CON *con;

  for (dat = A->dat; dat != A->end; dat ++)
  {
#if MPI
    if (!dat->con->dia) break; /* skip external */
#endif
    con = dat->con;
    COPY (dat->R0, con->R);
  }
}

/* GMRES based solver */
static int gmres_based_solve (PRIVATE *A, NEWTON *ns, LOCDYN *ldy)
{
  double *merit, step;
  char fmt [512];
  short dynamic;
  int div;

  sprintf (fmt, "NEWTON_SOLVER: delta: %%6g iteration: %%%dd merit: %%.2e\n", (int)log10 (ns->maxiter) + 1);
  ERRMEM (ns->merhist = realloc (ns->merhist, ns->maxiter * sizeof (double)));
  ERRMEM (ns->mvhist = realloc (ns->mvhist, ns->maxiter * sizeof (double)));
  dynamic = ldy->dom->dynamic;
  step = ldy->dom->step;
  merit = &ldy->dom->merit;
  *merit = MERIT_Function (ldy, 0);
  ns->iters = 0;
  div = 1;

  while (ns->iters < ns->maxiter && A->matvec < ns->maxmatvec && *merit > ns->meritval)
  {
    solve (A, PQN_GMRES, ns->linmaxiter, ns->epsilon, dynamic, step, ns->delta, 0.0, ns->omega, A->dr, A->rhs);

    U_WR_B (A, 0);

    *merit = MERIT_Function (ldy, 0);

    if (isnan (*merit)) return 0; /* XXX */

    ns->merhist [ns->iters] = *merit;
    ns->mvhist [ns->iters] = A->matvec;

#if MPI
    if (ldy->dom->rank == 0)
#endif
    if (ldy->dom->verbose && ns->iters % div == 0) printf (fmt, ns->delta, ns->iters, *merit), div *= 2;

    ns->iters ++;
  }

#if MPI
  if (ldy->dom->rank == 0)
#endif
  if (ldy->dom->verbose) printf (fmt, ns->delta, ns->iters, *merit);

  if (*merit > ns->meritval) return 0;
  else return 1;


  return 1;
}

/* diagonalized solver */
static int diagonalized_solve (PRIVATE *A, NEWTON *ns, LOCDYN *ldy)
{
  double *merit, prevm, step;
  char fmt [512];
  short dynamic;
  int div, gt;

  sprintf (fmt, "NEWTON_SOLVER: theta: %%6g iteration: %%%dd merit: %%.2e\n", (int)log10 (ns->maxiter) + 1);
  ERRMEM (ns->merhist = realloc (ns->merhist, ns->maxiter * sizeof (double)));
  ERRMEM (ns->mvhist = realloc (ns->mvhist, ns->maxiter * sizeof (double)));
  dynamic = ldy->dom->dynamic;
  step = ldy->dom->step;
  merit = &ldy->dom->merit;
  *merit = MERIT_Function (ldy, 0);
  ns->iters = 0;
  div = 1;
  gt = 0;

  while (ns->iters < ns->maxiter && *merit > ns->meritval)
  {
    solve (A, PQN_DIAG, 0, 0.0, dynamic, step, 0.0, ns->theta, ns->omega, A->dr, A->rhs);

    U_WR_B (A, 0);

    prevm = *merit;

    *merit = MERIT_Function (ldy, 0);

    if (isnan (*merit)) return 0; /* XXX */

    ns->merhist [ns->iters] = *merit;
    ns->mvhist [ns->iters] = 0;

#if MPI
    if (ldy->dom->rank == 0)
#endif
    if (ldy->dom->verbose && ns->iters % div == 0) printf (fmt, ns->theta, ns->iters, *merit), div *= 2;

    ns->iters ++;
  }

#if MPI
  if (ldy->dom->rank == 0)
#endif
  if (ldy->dom->verbose) printf (fmt, ns->theta, ns->iters, *merit);

  if (*merit > ns->meritval) return 0;
  else return 1;
}

/* create solver */
NEWTON* NEWTON_Create (double meritval, int maxiter)
{
  NEWTON *ns;

  ERRMEM (ns = MEM_CALLOC (sizeof (NEWTON)));
  ns->meritval = meritval;
  ns->maxiter = maxiter;
  ns->locdyn = LOCDYN_ON;
  ns->linver = PQN_GMRES;
  ns->linmaxiter = 10;
  ns->maxmatvec = ns->linmaxiter * maxiter;
  ns->epsilon = 0.25;
  ns->delta = 0.0;
  ns->theta = 0.25;
  ns->omega = meritval * 1E-3;
  ns->merhist = NULL;
  ns->mvhist = NULL;

  return ns;
}

/* run solver */
void NEWTON_Solve (NEWTON *ns, LOCDYN *ldy)
{
  PRIVATE *A;
  int ret;

  A = create_private_data (ns, ldy);

  switch (ns->linver)
  {
  case PQN_GMRES: ret = gmres_based_solve (A, ns, ldy); break;
  case PQN_DIAG: ret = diagonalized_solve (A, ns, ldy); break;
  }

  if (ret == 0)
  {
    if (ns->locdyn == LOCDYN_ON)
    {
      GAUSS_SEIDEL *gs;

#if MPI
      if (ldy->dom->rank == 0)
#endif
      if (ldy->dom->verbose) printf ("NEWTON_SOLVER:  ************ FAILED => switching to GAUSS_SEIDEL! ************ \n");

      restore_initial_R (A);
      gs = GAUSS_SEIDEL_Create (1.0, ns->maxiter, ns->meritval, GS_FAILURE_CONTINUE, 1E-9, 100, DS_SEMISMOOTH_NEWTON, NULL, NULL);
      GAUSS_SEIDEL_Solve (gs, ldy);
      GAUSS_SEIDEL_Destroy (gs);
    }
    else 
    {
      fprintf (stderr, "NEWTON_SOLVER: ************ FAILED in body space mode => ABORDING! ************\n");
      EXIT (0);
    }
  }

  destroy_private_data (A);
}

/* write labeled state values */
void NEWTON_Write_State (NEWTON *ns, PBF *bf)
{
  PBF_Label (bf, "NTITERS");
  PBF_Int (bf, &ns->iters, 1);
}

/* destroy solver */
void NEWTON_Destroy (NEWTON *ns)
{
  free (ns->merhist);
  free (ns->mvhist);
  free (ns);
}
