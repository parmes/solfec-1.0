/*
 * glu.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * linear gluing solver
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
#include "dom.h"
#include "alg.h"
#include "bla.h"
#include "lap.h"
#include "glu.h"
#include "err.h"
#include "ext/krylov/pcg.h"

#if MPI
#include "com.h"
#endif

#define BLOCKS 256

struct vect
{
  double *x;

  int n;
};

#define vect(ptr) ((struct vect*)ptr)

struct glue
{
  LOCDYN *ldy;

  SET *subset; /* GLUEPNT constraints */

  MEM setmem;

  struct vect *x, *b;

  hypre_PCGFunctions *pcg_functions;
  void *pcg_vdata;

#if MPI
  COMDATA *send, *recv;
  int nsend, nrecv;
  void *pattern;
  double **R; /* external con->R */
#endif
};

static double glue_stiffness (CON *con)
{
  BULK_MATERIAL *a = con->master->mat,
		*b = con->slave->mat;

  return 2.0 / (1.0/a->young + 1.0/b->young);
}

#if MPI
static void update_external_reactions (GLUE *glu, double *x)
{
  SET *item, *jtem;
  double **r, *z;
  int i, *j, *k;
  COMDATA *d;
  CON *con;

  for (item = SET_First (glu->subset), d = glu->send; item; item = SET_Next (item))
  {
    con = item->data;
    for (jtem = SET_First (con->ext); jtem; jtem = SET_Next (jtem), d ++)
    {
      d->d = &x [3*con->num];
    }
  }

  COMALL_Repeat (glu->pattern);

  for (i = 0, d = glu->recv, r = glu->R; i < glu->nrecv; i ++, d ++)
  {
    for (j = d->i, k = d->i + d->ints, z = d->d; j < k; j ++, r ++, z += 3)
    {
      COPY (z, *r);
    }
  }
}
#endif

static struct vect* newvect (int n)
{
  struct vect *v;

  ERRMEM (v = malloc (sizeof (struct vect)));
  ERRMEM (v->x = MEM_CALLOC (n * sizeof (double)));
  v->n = n;

  return v;
}

/* PCG interface start */
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

static int CommInfo (void *A, int *my_id, int *num_procs)
{
#if MPI
  GLUE *glu = (GLUE*)A;
  DOM *dom = glu->ldy->dom;

  *my_id = dom->rank;
  *num_procs = dom->ncpu;
#else
  *my_id = 0;
  *num_procs = 1;
#endif
  return 0;
}

static void* CreateVector (void *vector)
{
  struct vect *a = vect (vector), *v;

  ERRMEM (v = malloc (sizeof (struct vect)));
  ERRMEM (v->x = MEM_CALLOC (a->n * sizeof (double)));
  v->n = a->n;

  return v;
}

static int DestroyVector (void *vector)
{
  struct vect *v = vect (vector);
  free (v->x);
  free (v);

  return 0;
}

static void *MatvecCreate (void *A, void *x)
{
  return NULL;
}

static int Matvec (void *matvec_data, double alpha, void *A, void *vx, double beta, void *vy)
{
  GLUE *glu = (GLUE*)A;
  double *x = vect (vx)->x, *y = vect (vy)->x;
  double step = glu->ldy->dom->step;
  double *W, C, *z;
  SET *item;
  CON *con;
  OFFB *blk;
  DIAB *dia;

#if MPI
  update_external_reactions (glu, x); /* (###) */
#endif

  for (item = SET_First (glu->subset); item; item = SET_Next (item), y += 3)
  {
    con = item->data;
    dia = con->dia;
    W = dia->W;
    C = 1.0 / (step * glue_stiffness (dia->con));
    z = &x [3*con->num];

    y [0] = beta * y[0] + alpha * ((W[0]+C)*z[0] + W[3]*z[1] + W[6]*z[2]);
    y [1] = beta * y[1] + alpha * (W[1]*z[0] + (W[4]+C)*z[1] + W[7]*z[2]);
    y [2] = beta * y[2] + alpha * (W[2]*z[0] + W[5]*z[1] + (W[8]+C)*z[2]);

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
      if (con->kind == GLUEPNT)
      {
	W = blk->W;
	z = &x [3*con->num];
	y [0] += alpha * (W[0]*z[0] + W[3]*z[1] + W[6]*z[2]);
	y [1] += alpha * (W[1]*z[0] + W[4]*z[1] + W[7]*z[2]);
	y [2] += alpha * (W[2]*z[0] + W[5]*z[1] + W[8]*z[2]);
      }
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = CON (blk->dia);
      if (con->kind == GLUEPNT)
      {
	W = blk->W;
	z = con->R; /* (###) */
	y [0] += alpha * (W[0]*z[0] + W[3]*z[1] + W[6]*z[2]);
	y [1] += alpha * (W[1]*z[0] + W[4]*z[1] + W[7]*z[2]);
	y [2] += alpha * (W[2]*z[0] + W[5]*z[1] + W[8]*z[2]);
      }
    }
#endif
  }

  return 0;
}

static int MatvecDestroy (void *matvec_data)
{
  return 0;
}

static double InnerProd (void *vx, void *vy)
{
  double dot = 0.0, *x, *y, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n, y = vect (vy)->x; x < z; x ++, y ++)
  {
    dot += (*x) * (*y);
  }

#if MPI
  double val = dot;
  MPI_Allreduce (&val, &dot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  return dot;
}

static int CopyVector (void *vx, void *vy)
{
  double *x, *y, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n, y = vect (vy)->x; x < z; x ++, y ++)
  {
    (*y) = (*x);
  }

  return 0;
}

static int ClearVector (void *vx)
{
  double *x, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n; x < z; x ++)
  {
    (*x) = 0.0;
  }

  return 0;
}

static int ScaleVector (double alpha, void *vx)
{
  double *x, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n; x < z; x ++)
  {
    (*x) *= alpha;
  }

  return 0;
}

static int  Axpy (double alpha, void *vx, void *vy )
{
  double *x, *y, *z;

  for (x = vect (vx)->x, z = x + vect (vx)->n, y = vect (vy)->x; x < z; x ++, y ++)
  {
    (*y) += alpha * (*x);
  }

  return 0;
}

static int PrecondSetup (void *vdata, void *A, void *vb, void *vx)
{
  return 0;
}

static int Precond (void *vdata, void *A, void *vb, void *vx)
{
  return CopyVector (vb, vx);
}
/* PCG interface end */

static void compute_right_hand_side (SET *subset, double *b)
{
  double *y, *B, *W, *R;
  SET *item;
  CON *con;
  DIAB *dia;
  OFFB *blk;

  for (item = SET_First (subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;
    y = &b [3 * con->num];
    B = dia->B;
    COPY (B, y);

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
      if (con->kind != GLUEPNT)
      {
	W = blk->W;
	R = con->R;
	NVADDMUL (y, W, R, y);
      }
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = CON (blk->dia);
      if (con->kind != GLUEPNT)
      {
	W = blk->W;
	R = con->R;
	NVADDMUL (y, W, R, y);
      }
    }
#endif
    SCALE (y, -1.0);
  }
}

/* create gluing solver */
GLUE* GLUE_Create (LOCDYN *ldy)
{
  DOM *dom = ldy->dom;
  GLUE *glu;
  CON *con;
  int ncon;

  ncon = dom->ncon;
  ERRMEM (glu = MEM_CALLOC (sizeof (GLUE)));
  MEM_Init (&glu->setmem, sizeof (SET), MAX (ncon, BLOCKS));
  glu->ldy = ldy;

  /* collect gluing constraints */
  for (con = dom->con; con; con = con->next)
  {
    if (con->kind == GLUEPNT)
    {
      SET_Insert (&glu->setmem, &glu->subset, con, NULL);
    }
  }

  /* gluing constraints local numbering */
  DOM_Number_Constraints (dom, 1, glu->subset);

  /* unknown and right hand size */
  glu->x = newvect (3 * SET_Size (glu->subset));
  glu->b = CreateVector (glu->x);

  /* create PCG solver */
  glu->pcg_functions = hypre_PCGFunctionsCreate (CAlloc, Free, CommInfo, CreateVector, DestroyVector, MatvecCreate,
    Matvec, MatvecDestroy, InnerProd, CopyVector, ClearVector, ScaleVector, Axpy, PrecondSetup, Precond);
  glu->pcg_vdata = hypre_PCGCreate (glu->pcg_functions);

#if MPI
  /* communication pattern */
  SET *item, *jtem;
  int i, *j, *k, n;
  COMDATA *d;
  double **r;

  for (item = SET_First (glu->subset); item; item = SET_Next (item))
  {
    con = item->data;
    for (jtem = SET_First (con->ext); jtem; jtem = SET_Next (jtem))
    {
      glu->nsend ++;
    }
  }

  ERRMEM (glu->send = MEM_CALLOC (glu->nsend * sizeof (COMDATA)));
  d = glu->send;

  for (item = SET_First (glu->subset); item; item = SET_Next (item))
  {
    con = item->data;
    for (jtem = SET_First (con->ext); jtem; jtem = SET_Next (jtem), d ++)
    {
      d->rank = (int) (long) jtem->data;
      d->doubles = 3;
      d->ints = 1;
      d->d = &vect (glu->x)->x [3*con->num];
      d->i = (int*) &con->id;
    }
  }

  glu->pattern = COMALL_Pattern (MPI_COMM_WORLD, glu->send, glu->nsend, &glu->recv, &glu->nrecv);

  COMALL_Repeat (glu->pattern);

  /* map received reaction pointers */
  for (i = n = 0, d = glu->recv; i < glu->nrecv; i ++, d ++)
  {
    for (j = d->i, k = j + d->ints; j < k; j ++)
    {
      n ++;
    }
  }

  ERRMEM (glu->R = MEM_CALLOC (n * sizeof (double*)));
  r = glu->R;

  for (i = 0, d = glu->recv; i < glu->nrecv; i ++, d ++)
  {
    for (j = d->i, k = j + d->ints; j < k; j ++, r ++)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*j), NULL), "Invalid constraint id");
      (*r) = con->R;
    }
  }
#endif

  return glu;
}

/* compute gluing reactions */
void GLUE_Solve (GLUE *glu, double abstol, int maxiter)
{
  double *x, *z, *R;
  SET *item;
  CON *con;

  compute_right_hand_side (glu->subset, glu->b->x);

  hypre_PCGSetTol (glu->pcg_vdata, 0.0);
  hypre_PCGSetMaxIter (glu->pcg_vdata, maxiter);
  hypre_PCGSetAbsoluteTol (glu->pcg_vdata, abstol);
  hypre_PCGSetup (glu->pcg_vdata, glu, glu->b, glu->x);
  hypre_PCGSolve (glu->pcg_vdata, glu, glu->b, glu->x);

  for (item = SET_First (glu->subset), x = glu->x->x; item; item = SET_Next (item))
  {
    con = item->data;
    z = &x [3*con->num];
    R = con->R;
    COPY (z, R);
  }
}

#if MPI
/* update external gluing reactions */
void GLUE_Update_External_Reactions (GLUE *glu)
{
  update_external_reactions (glu, glu->x->x);
}
#endif

/* destroy gluing solver */
void GLUE_Destroy (GLUE *glu)
{
  MEM_Release (&glu->setmem);
  DestroyVector (glu->x);
  DestroyVector (glu->b);
  hypre_PCGDestroy (glu->pcg_vdata);
#if MPI
  COMALL_Free (glu->pattern);
  free (glu->send);
  free (glu->recv);
  free (glu->R);
#endif
  free (glu);
}
