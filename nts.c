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
#include "vic.h"
#include "mrf.h"
#include "ext/krylov/krylov.h"

typedef struct con_data CON_DATA;
typedef struct private PRIVATE;

struct con_data
{
  MX *mH; /* master H */
  int mi, /* shift to mH in global velocity space */
     *mj; /* index mapping for mH converted by csc_to_dense () */

  MX *sH;
  int si,
     *sj;

  CON *con; /* constraint */

  double DR [3], /* reaction increment */
	 R0 [3]; /* initial reaction */
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

/* U = W R + B */
static void U_WR_B (PRIVATE *A)
{
#if MPI
  DOM_Update_External_Reactions (A->dom, 0);
#endif

  if (A->ns->locdyn == LOCDYN_ON)
  {
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
      B = dia->B;
      NVADDMUL (B, W, R, U);
      for (blk = dia->adj; blk; blk = blk->n)
      {
	R = blk->dia->R;
	W = blk->W;
	NVADDMUL (U, W, R, U);
      }
#if MPI
      for (blk = dia->adjext; blk; blk = blk->n)
      {
	R = CON (blk->dia)->R;
	W = blk->W;
	NVADDMUL (U, W, R, U);
      }
#endif
    }
  }
  else
  {
    double *B, *U, step;
    CON_DATA *dat;
    MAP *item;
    DIAB *dia;
    CON *con;
    int n;

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
      B = dia->B;
      U = con->U;
      COPY (B, U);
    }

    H_times_u (A->a, A->dat, A->end, A->u);
  }
}

/* create constraints data for body-space mode */
static int body_space_constraints_data (DOM *dom, PRIVATE *A)
{
  LOCDYN *ldy = dom->ldy;
  MAP *item, *fem;
  int ndat, ret;
  short dynamic;
  CON_DATA *dat;
  double step;
  CON *con;
#if MPI
  MAP *jtem;
#endif

  ret = 0;
  fem = NULL;
  dynamic = dom->dynamic;
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
    if (dynamic && con->kind == CONTACT && con->gap > 0.0)
    {
      SET (con->R, 0.0); continue; /* skip open dynamic contacts */
    }

    DIAB *dia = con->dia;
    BODY *m = con->master,
	 *s = con->slave;
    void *mgobj = mgobj(con),
	 *sgobj;
    SHAPE *mshp = mshp(con),
	  *sshp;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base;
    MX_DENSE_PTR (W, 3, 3, dia->W);
    MX *prod, *inv;

#if MPI
    if (m->flags & BODY_CHILD)
    {
      if (dynamic) BODY_Dynamic_Init (m);
      else BODY_Static_Init (m);
    }
#endif

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
      sgobj = sgobj(con);
      sshp = sshp(con);

#if MPI
      if (s->flags & BODY_CHILD)
      {
	if (dynamic) BODY_Dynamic_Init (s);
	else BODY_Static_Init (s);
      }
#endif

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

    if (m->kind != OBS)
    {
      dat->mH = BODY_Gen_To_Loc_Operator (m, mshp, mgobj, mpnt, base);
      dat->mi = (int) (long) MAP_Find (A->bod, m, NULL);

      if (m->kind == FEM) inv = MAP_Find (fem, m, NULL); else inv = m->inverse;
      prod = MX_Matmat (1.0, inv, MX_Tran (dat->mH), 0.0, NULL);
      MX_Matmat (step, dat->mH, prod, 0.0, &W); /* H * inv (M) * H^T */
      MX_Destroy (prod);

      if (m->kind == FEM)
      {
	dat->mH = csc_to_dense (dat->mH, &dat->mj);
      }
    }

    if (s && s->kind != OBS)
    {
      dat->sH = BODY_Gen_To_Loc_Operator (s, sshp, sgobj, spnt, base);
      dat->si = (int) (long) MAP_Find (A->bod, s, NULL);
      MX_Scale (dat->sH, -1.0);

      if (s->kind == FEM) inv = MAP_Find (fem, s, NULL); else inv = s->inverse;
      prod = MX_Matmat (1.0, inv, MX_Tran (dat->sH), 0.0, NULL);
      MX_Matmat (step, dat->sH, prod, 1.0, &W); /* H * inv (M) * H^T */
      MX_Destroy (prod);

      if (s->kind == FEM)
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
    if (dynamic && con->kind == CONTACT && con->gap > 0.0) continue; /* skip open dynamic contacts */

    BODY *m = con->master,
	 *s = con->slave;
    void *mgobj = mgobj(con),
	 *sgobj;
    SHAPE *mshp = mshp(con),
	  *sshp;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base;

    if (s)
    {
      sgobj = sgobj(con);
      sshp = sshp(con);
    }

    if ((jtem = MAP_Find_Node (A->bod, m, NULL)))
    {
      dat->mH = BODY_Gen_To_Loc_Operator (m, mshp, mgobj, mpnt, base);
      dat->mi = (int) (long) jtem->data;

      if (m->kind == FEM)
      {
	dat->mH = csc_to_dense (dat->mH, &dat->mj);
      }
    }

    if ((jtem = MAP_Find_Node (A->bod, s, NULL)))
    {
      dat->sH = BODY_Gen_To_Loc_Operator (s, sshp, sgobj, spnt, base);
      dat->si = (int) (long) jtem->data;
      MX_Scale (dat->sH, -1.0);

      if (s->kind == FEM)
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
    if (dynamic && con->kind == CONTACT && con->gap > 0.0)
    {
      SET (con->R, 0.0);
      continue; /* skip open dynamic contacts */
    }

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

  U_WR_B (A); /* U = W R + B */

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
  free (A);
}

/* single projected semi-Newton step */
static void solve (CON_DATA *dat, CON_DATA *end, short dynamic, double step, double theta, double epsilon)
{
  double T [9], b [3], gamma = 1.0 - theta;
  int ipiv [3];

  for (; dat != end; dat ++)
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
	   *DR = dat->DR;

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
    case CONTACT:
    {
      double X [9], Y [9];

      VIC_Linearize (dat->con, U, R, -1, epsilon, b, X, Y);
      SCALE (b, -1.0);

      NNMUL (X, W, T);
      NNADD (T, Y, T);
    }
    break;
    }

    ASSERT (lapack_dgesv (3, 1, T, 3, ipiv, b, 3) == 0, ERR_MTX_LU_FACTOR);
    DR [0] = gamma * DR[0] + theta * b[0];
    DR [1] = gamma * DR[1] + theta * b[1];
    DR [2] = gamma * DR[2] + theta * b[2];
    ACC (DR, R);

    if (con->kind == CONTACT)
    {
      double c = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area;
      VIC_Project (con->mat.base->friction, c, R, R);
    }
  }
}

/* reset solution */
static void reset (PRIVATE *A)
{
  CON_DATA *dat;
  CON *con;

  for (dat = A->dat; dat != A->end; dat ++)
  {
#if MPI
    if (!dat->con->dia) break; /* skip external */
#endif
    con = dat->con;
    SET (dat->DR, 0.0);
    COPY (dat->R0, con->R);
  }

  U_WR_B (A);
}

/* create solver */
NEWTON* NEWTON_Create (double meritval, int maxiter)
{
  NEWTON *ns;

  ERRMEM (ns = MEM_CALLOC (sizeof (NEWTON)));
  ns->meritval = meritval;
  ns->maxiter = maxiter;
  ns->locdyn = LOCDYN_ON;
  ns->theta = 0.25;
  ns->epsilon = 1E-9;
  ns->presmooth = 10;

  return ns;
}

/* run solver */
void NEWTON_Solve (NEWTON *ns, LOCDYN *ldy)
{
  double *merit, prevm, step, theta0, merit0;
  GAUSS_SEIDEL *gs;
  char fmt [512];
  short dynamic;
  int div, gt;
  PRIVATE *A;

  if (ns->locdyn == LOCDYN_ON && ns->presmooth > 0)
  {
    gs = GAUSS_SEIDEL_Create (1E-10, ns->presmooth, 1.0, GS_FAILURE_CONTINUE, 1E-9, 100, DS_SEMISMOOTH_NEWTON, NULL, NULL);
    gs->verbose = 0;
    gs->nomerit = 1;
#if MPI
    if (ldy->dom->rank == 0)
#endif
    if (ldy->dom->verbose)
    {
      printf ("NEWTON_SOLVER: presmoothing ");
      for (gt = 0; gt < ns->presmooth; gt ++) printf (".");
      printf ("\n");
    }
    GAUSS_SEIDEL_Solve (gs, ldy);
    GAUSS_SEIDEL_Destroy (gs);
  }

  sprintf (fmt, "NEWTON_SOLVER: theta: %%6g iteration: %%%dd merit: %%.2e\n", (int)log10 (ns->maxiter) + 1);
  ERRMEM (ns->merhist = realloc (ns->merhist, ns->maxiter * sizeof (double)));
  A = create_private_data (ns, ldy);
  dynamic = ldy->dom->dynamic;
  merit = &ldy->dom->merit;
  step = ldy->dom->step;
  *merit = MERIT_Function (ldy, 0);
  theta0 = ns->theta;
  merit0 = *merit;
  ns->iters = 0;
  div = 1;
  gt = 0;

  while (ns->iters < ns->maxiter && *merit > ns->meritval)
  {
    solve (A->dat, A->end, dynamic, step, ns->theta, ns->epsilon);

    U_WR_B (A);

    prevm = *merit;

    *merit = MERIT_Function (ldy, 0);

    ns->merhist [ns->iters] = *merit;

    if (*merit > prevm && ++gt > 10 && *merit > 10)
    {
      if (ns->theta < 0.0009765625) ns->theta = 0.5; /* < 0.5**10 */
      else ns->theta *= 0.5;
      reset (A);
      gt = 0;
    }

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

  if (*merit > merit0)
  {
    reset (A);

    *merit = MERIT_Function (ldy, 0);

#if MPI
    if (ldy->dom->rank == 0)
#endif
    if (ldy->dom->verbose) printf ("NEWTON_SOLVER: DIVERGED => Reusing previous solution (merit: %.2e)\n", *merit);
  }

  destroy_private_data (A);

  ns->theta = theta0;
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
  free (ns);
}
