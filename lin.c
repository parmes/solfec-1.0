/*
 * lin.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * constraints linearization routines
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

#if MPI
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_parcsr_ls.h>
#include <HYPRE_krylov.h>
#include <HYPRE_IJ_mv.h>
#include <HYPRE.h>
#endif

#if SPQR
#include <SuiteSparseQR_C.h>
#elif UMFPACK
#include <umfpack.h>
#else
#include "lss.h"
#endif

#include <complex.h>
#include <stdlib.h>
#include <time.h>
#include "dom.h"
#include "alg.h"
#include "bla.h"
#include "lap.h"
#include "lin.h"
#include "err.h"

#define MEMBLK 256

#if MPI
typedef struct hypre_data HYPRE_DATA;

struct hypre_data /* HYPRE IJ matrix data and the rest of linear system */
{
  MPI_Comm  comm; /* communicator (skips processors with no constrains) */

  int ilower, /* lower row index == dom->con->num * 3 */
      iupper; /* upper low index == ilower + dom->ncon * 3 - 1 */

  int jlower, /* lower column index == ilower */
      jupper; /* upper column index == iupper */

  int *ncols, /* number of columns in each row */
      *rows;  /* row indices: ilower, ilower + 1, ..., iupper */

  HYPRE_IJMatrix A;  /* HYPRE IJ matrix composed from the above data */
  HYPRE_ParCSRMatrix AP; /* reference to the parallel HYPRE CSR image of A */

  HYPRE_IJVector B; /* HYPRE IJ right hand side vector */
  HYPRE_ParVector BP; /* HYPRE parallel image of B */

  HYPRE_IJVector X; /* HYPRE IJ solution vector */
  HYPRE_ParVector XP; /* HYPRE parallel image of X */
};
#endif

struct linsys
{
  LINVAR variant; /* linearization variant */

  LINOPT options; /* linearization options */

  LOCDYN *ldy; /* local dynamics */

  MEM mapmem, /* memory of DIAB->map and OFFB->map int [3] vectors */
      datmem, /* constraint linearization data memory */
      ymem;   /* auxiliary con->Y storage memory */

  double *a; /* compressed column values */

  int *rows, /* row indices/pointers */
      *cols, /* column pointers/indices */
       nnz, /* number of nonzero entries */
       dim; /* dimension of the rectangular matrix */

  double *b; /* right hand side */

  double *x; /* solution */

#if MPI
  HYPRE_DATA *hyp;
#endif

#if SPQR
  cholmod_sparse *spqr_a; /* A */
  cholmod_dense *spqr_b; /* b */
  cholmod_common spqr_c; /* solver data */
#elif UMFPACK
  void *symbolic; /* symbolic factorization */
#else
  void *lss; /* LSS interface */
#endif

  int iters; /* iterations count */

  double resnorm; /* residual norm */
};

/* forward coordinates change */
static void variational_coord_change_forward (LINSYS *sys);

/* create linear system resulting from constraints linearization */
LINSYS* LINSYS_Create (LINVAR variant, LINOPT options, LOCDYN *ldy)
{
  int i, j, k, l, n;
  LINSYS *sys;
  MAP *item;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  DOM *dom;
  MEM mem;
#if MPI
  int *diacnt, /* this-processor coumn counts per row */
      *offcnt; /* off-processor column counts per row */
#endif

  dom = ldy->dom;
    
  ERRMEM (sys = MEM_CALLOC (sizeof (LINSYS)));

  sys->variant = variant;
  sys->options = options;
  sys->ldy = ldy;

  MEM_Init (&sys->mapmem, sizeof (int [3]), MEMBLK);
  MEM_Init (&sys->datmem, sizeof (LINDATA), MEMBLK);
  MEM_Init (&sys->ymem, sizeof (double [4]), MEMBLK);

  /* constraints numbering */
  DOM_Number_Constraints (dom, options & LOCAL_SYSTEM);

  /* auxiliary constraint space */
  for (con = dom->con; con; con = con->next)
  {
    ERRMEM (con->lin = MEM_Alloc (&sys->datmem));
    ERRMEM (con->Y = MEM_Alloc (&sys->ymem));
  }
#if MPI
  if (!(options & LOCAL_SYSTEM))
  {
    for (item = MAP_First (dom->conext); item; item = MAP_Next (item))
    {
      con = item->data;
      ERRMEM (con->Y = MEM_Alloc (&sys->ymem));
    }
  }
#endif

  /* system matrix */
#if MPI
  if (options & LOCAL_SYSTEM) /* create column mapping in DIAB->map and OFFB->map */
  {
#endif
    MAP **pcol;

    /* allocate colum blocks mapping */
    ERRMEM (pcol = MEM_CALLOC (dom->ncon * sizeof (MAP*)));
    MEM_Init (&mem, sizeof (MAP), MEMBLK);

    /* map colum blocks */
    for (dia = ldy->dia; dia; dia = dia->n)
    {
      int col, row, *map;

      col = row = dia->con->num;
      ASSERT_DEBUG (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Diagonal block mapped twice");
      ERRMEM (map = MEM_Alloc (&sys->mapmem));
      ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL)); /* '+1' not to confuse with NULL */

      for (blk = dia->adj; blk; blk = blk->n)
      {
	col = blk->dia->con->num;

	if (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL)) /* there are two W(i,j) blocks for two-body constraints */
	{
	  ERRMEM (map = MEM_Alloc (&sys->mapmem));
	  ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL));
	}
      }
    }

    for (j = 0, sys->nnz = 0; j < dom->ncon; j ++)
    {
      sys->nnz += 9 * MAP_Size (pcol [j]); /* number of nonzeros */
    }
    sys->dim = dom->ncon * 3; /* system dimension */

    /* eallocate compressed column storage */
    ERRMEM (sys->a = malloc (sizeof (double) * sys->nnz));
    ERRMEM (sys->rows = malloc (sizeof (int) * sys->nnz));
    ERRMEM (sys->cols = MEM_CALLOC (sizeof (int) * (sys->dim+1))); /* '+ 1' as there is cols[0] == 0 always */

    int *rows = sys->rows,
	*cols = sys->cols,
	*aux;

    ERRMEM (aux = malloc (sizeof (int) * sys->dim));

    /* set up column sizes */
    for (j = 0; j < dom->ncon; j ++)
    {
      n = 3 * MAP_Size (pcol [j]);

      cols [3*j+1] = n; /* '+1' so that cols[0] == 0 */
      cols [3*j+2] = n;
      cols [3*j+3] = n;
    }

    /* compute column pointers */
    for (n = 1; n <= sys->dim; n ++)
    {
      cols [n] += cols [n-1];
      aux [n-1] = cols [n-1]; /* initialize aux with cols  */
    }

    ASSERT_DEBUG (cols [sys->dim] == sys->nnz, "Inconsistent sparse storage");

    /* set up row pointers */
    for (j = 0; j < dom->ncon; j ++)
    {
      for (item = MAP_First (pcol [j]); item; item = MAP_Next (item))
      {
	k = 3 * (((int) (long) item->key) - 1); /* base row block index (-1 as keys == row+1) */

	int *map = item->data;

	ASSERT_DEBUG (map [0] == 0 && map [1] == 0 && map [2] == 0, "Reusing the same map triplet: %d, %d, %d", map [0], map [1], map [2]);

	for (n = 0; n < 3; n ++) map [n] = aux [3*j+n]; /* set row blocks map */

	for (n = 0; n < 3; n ++) /* for each column block sub-column */
	{
	  for (l = aux [3*j+n], i = 0; i < 3; i ++) /* for each row block sub-row */
	  {
	    rows [l + i] = k + i; /* set up row index */
	  }
	  aux [3*j+n] += 3; /* increment relative row pointers */
	}
      }
    }

    /* assign sparse column maps to blocks */
    for (dia = ldy->dia; dia; dia = dia->n)
    {
      int col, row;

      col = row = dia->con->num;

      ASSERT_DEBUG_EXT (dia->map = MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Inconsistent sparse storage");

      for (blk = dia->adj; blk; blk = blk->n)
      {
	col = blk->dia->con->num;

	ASSERT_DEBUG_EXT (blk->map = MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Inconsistent sparse storage");
      }
    }

    /* free auxiliary space */
    MEM_Release (&mem);
    free (pcol);
    free (aux);
#if MPI
  }
  else /* create row mapping in DIAB->map and OFFB->map */
  {
    int basenum = 0;
    MAP **prow;

    /* first constraint number */
    if (dom->con) basenum = dom->con->num;

    /* allocate colum blocks mapping */
    ERRMEM (diacnt = MEM_CALLOC (dom->ncon * sizeof (int [3])));
    ERRMEM (offcnt = MEM_CALLOC (dom->ncon * sizeof (int [3])));
    ERRMEM (prow = MEM_CALLOC (dom->ncon * sizeof (MAP*)));
    MEM_Init (&mem, sizeof (MAP), MEMBLK);

    /* map colum blocks */
    for (dia = ldy->dia; dia; dia = dia->n)
    {
      int col, row, *map;

      col = dia->con->num;
      row = col - basenum; /* use relative numbering for rows */
      ASSERT_DEBUG (!MAP_Find (prow [row], (void*) (long) (col+1), NULL), "Diagonal block mapped twice");
      ERRMEM (map = MEM_Alloc (&sys->mapmem));
      ERRMEM (MAP_Insert (&mem, &prow [row], (void*) (long) (col+1), map, NULL)); /* '+1' not to confuse with NULL */

      for (i = k = row * 3, j = k + 3; i < j; i ++) diacnt [i] = 3; /* diagonal blocks contribute 3 this-processor columns */

      for (blk = dia->adj; blk; blk = blk->n)
      {
	col = blk->dia->con->num;

	if (!MAP_Find (prow [row], (void*) (long) (col+1), NULL)) /* there are two W(i,j) blocks for two-body constraints */
	{
	  ERRMEM (map = MEM_Alloc (&sys->mapmem));
	  ERRMEM (MAP_Insert (&mem, &prow [row], (void*) (long) (col+1), map, NULL));

          for (i = k; i < j; i ++) diacnt [i] += 3; /* off-diagonal local blocks contribute 3 this-processor columns */
	}
      }
      for (blk = dia->adjext; blk; blk = blk->n)
      {
	col = ((CON*)blk->dia)->num; /* external constraint */

	if (!MAP_Find (prow [row], (void*) (long) (col+1), NULL)) /* there are two W(i,j) blocks for two-body constraints */
	{
	  ERRMEM (map = MEM_Alloc (&sys->mapmem));
	  ERRMEM (MAP_Insert (&mem, &prow [row], (void*) (long) (col+1), map, NULL));

          for (i = k; i < j; i ++) offcnt [i] += 3; /* off-diagonal external blocks contribute 3 off-processor columns */
	}
      }
    }

    for (j = 0, sys->nnz = 0; j < dom->ncon; j ++)
    {
      sys->nnz += 9 * MAP_Size (prow [j]); /* number of nonzeros */
    }
    sys->dim = dom->ncon * 3; /* system dimension */

    /* eallocate compressed column storage */
    ERRMEM (sys->a = malloc (sizeof (double) * sys->nnz));
    ERRMEM (sys->cols = malloc (sizeof (int) * sys->nnz));
    ERRMEM (sys->rows = MEM_CALLOC (sizeof (int) * (sys->dim+1))); /* '+ 1' as there is cols[0] == 0 always */

    int *rows = sys->rows,
	*cols = sys->cols,
	*aux;

    ERRMEM (aux = malloc (sizeof (int) * sys->dim));

    /* set up rows sizes */
    for (j = 0; j < dom->ncon; j ++)
    {
      n = 3 * MAP_Size (prow [j]);

      rows [3*j+1] = n; /* '+1' so that cols[0] == 0 */
      rows [3*j+2] = n;
      rows [3*j+3] = n;
    }

    /* compute rows pointers */
    for (n = 1; n <= sys->dim; n ++)
    {
      rows [n] += rows [n-1];
      aux [n-1] = rows [n-1]; /* initialize aux with cols  */
    }

    ASSERT_DEBUG (rows [sys->dim] == sys->nnz, "Inconsistent sparse storage");

    /* set up column pointers */
    for (j = 0; j < dom->ncon; j ++)
    {
      for (item = MAP_First (prow [j]); item; item = MAP_Next (item))
      {
	k = 3 * (((int) (long) item->key) - 1); /* base column block index (-1 as keys == col+1) */

	int *map = item->data;

	ASSERT_DEBUG (map [0] == 0 && map [1] == 0 && map [2] == 0, "Reusing the same map triplet: %d, %d, %d", map [0], map [1], map [2]);

	for (n = 0; n < 3; n ++) map [n] = aux [3*j+n]; /* set column blocks map */

	for (n = 0; n < 3; n ++) /* for each row block sub-row */
	{
	  for (l = aux [3*j+n], i = 0; i < 3; i ++) /* for each column block sub-column */
	  {
	    cols [l + i] = k + i; /* set up column index */
	  }
	  aux [3*j+n] += 3; /* increment relative column pointers */
	}
      }
    }

    /* assign sparse row maps to blocks */
    for (dia = ldy->dia; dia; dia = dia->n)
    {
      int col, row;

      col = dia->con->num;
      row = col - basenum; /* use relative numbering for rows */

      ASSERT_DEBUG_EXT (dia->map = MAP_Find (prow [row], (void*) (long) (col+1), NULL), "Inconsistent sparse storage");

      for (blk = dia->adj; blk; blk = blk->n)
      {
	col = blk->dia->con->num;

	ASSERT_DEBUG_EXT (blk->map = MAP_Find (prow [row], (void*) (long) (col+1), NULL), "Inconsistent sparse storage");
      }
      for (blk = dia->adjext; blk; blk = blk->n)
      {
	col = ((CON*)blk->dia)->num; /* external constraint */

	ASSERT_DEBUG_EXT (blk->map = MAP_Find (prow [row], (void*) (long) (col+1), NULL), "Inconsistent sparse storage");
      }
    }

    /* free auxiliary space */
    MEM_Release (&mem);
    free (prow);
    free (aux);
  }
#endif

  /* right hand side */
  ERRMEM (sys->b = MEM_CALLOC (sizeof (double) * sys->dim));

  /* solution vector */
  ERRMEM (sys->x = MEM_CALLOC (sizeof (double) * sys->dim));

  /* linear solver */
#if MPI
  if (options & LOCAL_SYSTEM) /* create column mapping in DIAB->map and OFFB->map */
  {
#endif

#if SPQR
    /* start CHOLMOD */
    cholmod_l_start (&sys->spqr_c);

    /* allocate A */
    ERRMEM (sys->spqr_a = MEM_CALLOC (sizeof (cholmod_sparse)));
    sys->spqr_a->nrow = sys->spqr_a->ncol = sys->dim;
    sys->spqr_a->nzmax = sys->nnz;
    sys->spqr_a->p = sys->cols;
    sys->spqr_a->i = sys->rows;
    sys->spqr_a->x = sys->a;
    sys->spqr_a->stype = 0;
    sys->spqr_a->itype = CHOLMOD_INT;
    sys->spqr_a->xtype = CHOLMOD_REAL;
    sys->spqr_a->dtype = CHOLMOD_DOUBLE;
    sys->spqr_a->sorted = TRUE;
    sys->spqr_a->packed = TRUE;

    /* allocate b */
    ERRMEM (sys->spqr_b = MEM_CALLOC (sizeof (cholmod_sparse)));
    sys->spqr_b->nrow = sys->spqr_b->nzmax = sys->spqr_b->d = sys->dim;
    sys->spqr_b->ncol = 1;
    sys->spqr_b->x = sys->b;
    sys->spqr_b->xtype = CHOLMOD_REAL;
    sys->spqr_b->dtype = CHOLMOD_DOUBLE;
#elif UMFPACK
    ASSERT (umfpack_di_symbolic (sys->dim, sys->dim, sys->cols, sys->rows,
            sys->a, &sys->symbolic, NULL, NULL) == UMFPACK_OK, ERR_UMFPACK_SOLVE);
#else
    ERRMEM (sys->lss = LSS_Create (sys->dim, sys->a, sys->cols, sys->rows));
    LSS_Set (sys->lss, LSS_RESTART, 100); /* FIXME: make user adjustable or automatic */
#endif

#if MPI
  }
  else /* initialize HYPRE */
  {
    MPI_Group allcpus, nonzero;
    int *ncon, *ranks;
    HYPRE_DATA *hyp;

    ERRMEM (ranks = MEM_CALLOC (sizeof (int [dom->ncpu])));
    ERRMEM (ncon = MEM_CALLOC (sizeof (int [dom->ncpu])));
    ERRMEM (hyp = MEM_CALLOC (sizeof (HYPRE_DATA)));
    sys->hyp = hyp;

    /* gather numbers of constraints from all subdomains */
    MPI_Allgather (&dom->ncon, 1, MPI_INT, ncon, 1, MPI_INT, MPI_COMM_WORLD);
    for (i = j = 0; i < dom->ncpu; i ++)
    {
      if (ncon [i]) ranks [j ++] = i;
    }

    /* create a group of processorss with nonzero constraint counts */
    MPI_Comm_group (MPI_COMM_WORLD, &allcpus);
    MPI_Group_incl (allcpus, j, ranks, &nonzero);

    /* create HYPRE communicator */
    MPI_Comm_create (MPI_COMM_WORLD, nonzero, &hyp->comm);

    /* clean */
    free (ranks);
    free (ncon);

    if (dom->con) /* there are some constraints here */
    {
      hyp->ilower = dom->con->num * 3;
      hyp->iupper = hyp->ilower + dom->ncon * 3 - 1;
      hyp->jlower = hyp->ilower;
      hyp->jupper = hyp->iupper;
      ERRMEM (hyp->ncols = malloc (sizeof (int [sys->dim])));
      ERRMEM (hyp->rows = malloc (sizeof (int [sys->dim])));
      for (j = 0, k = hyp->ilower; j < sys->dim; j ++, k ++)
      {
	hyp->ncols [j] = sys->rows [j+1] - sys->rows [j];
	hyp->rows [j] = k;
      }
    }

    if (hyp->comm != MPI_COMM_NULL) /* skip empty constraints set case */
    {
      /* create distributed matrix */
      HYPRE_IJMatrixCreate (hyp->comm, hyp->ilower, hyp->iupper, hyp->jlower, hyp->jupper, &hyp->A);
      HYPRE_IJMatrixSetObjectType (hyp->A, HYPRE_PARCSR);

      /* initialize row sizes and per/off-processor column counts to improve assembly efficiency */
      HYPRE_IJMatrixSetRowSizes (hyp->A, hyp->ncols);
      HYPRE_IJMatrixSetDiagOffdSizes (hyp->A, diacnt, offcnt);
      free (diacnt);
      free (offcnt);

      /* initialize matrix with zeros */
      HYPRE_IJMatrixInitialize (hyp->A);
      HYPRE_IJMatrixSetValues (hyp->A, sys->dim, hyp->ncols, hyp->rows, sys->cols, sys->a);
      HYPRE_IJMatrixAssemble (hyp->A);
      HYPRE_IJMatrixGetObject (hyp->A, (void **) &hyp->AP);

      /* b */
      HYPRE_IJVectorCreate (hyp->comm, hyp->jlower, hyp->jupper, &hyp->B);
      HYPRE_IJVectorSetObjectType (hyp->B, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize (hyp->B);
      HYPRE_IJVectorSetValues (hyp->B, sys->dim, hyp->rows, sys->b);
      HYPRE_IJVectorAssemble (hyp->B);
      HYPRE_IJVectorGetObject(hyp->B, (void **) &hyp->BP);

      /* x */
      HYPRE_IJVectorCreate (hyp->comm, hyp->jlower, hyp->jupper, &hyp->X);
      HYPRE_IJVectorSetObjectType (hyp->X, HYPRE_PARCSR);
      HYPRE_IJVectorInitialize (hyp->X);
      HYPRE_IJVectorSetValues (hyp->X, sys->dim, hyp->rows, sys->x);
      HYPRE_IJVectorAssemble (hyp->X);
      HYPRE_IJVectorGetObject(hyp->X, (void **) &hyp->XP);
    }
  }
#endif

  if (sys->variant == SMOOTHED_VARIATIONAL ||
      sys->variant == NONSMOOTH_VARIATIONAL)
  {
    variational_coord_change_forward (sys); /* R = R + m(R) */
  }

  return sys; 
}

/* update linear system for NT_NONSMOOTH_HSW, NT_NONSMOOTH_HYBRID, NT_FIXED_POINT variants */
static void system_update_HSW_HYBRID_FIXED (LINVAR variant, LINOPT options, LOCDYN *ldy, double *a, double *b)
{
  double d [3], norm, lim, udash, step;
  short dynamic, pull, global;
  DIAB *dia;
  OFFB *blk;
  DOM *dom;
  CON *con;

  dom = ldy->dom;
  dynamic = dom->dynamic;
  step = dom->step;
  global = !(options & LOCAL_SYSTEM);

#if MPI
  /* update external reactions (needed below) */
  if (global) DOM_Update_External_Reactions (ldy->dom, 0);
#endif

  for (con = dom->con; con; con = con->next, b += 3) /* constraints numbering corresponds to consecutive b[3]s */
  {
    dia = con->dia;

    double *W = dia->W,
	   *B = dia->B,
	   *U = dia->U,
	   *R = dia->R,
	   *RE = RE(con);

    /* update residual */
    NVADDMUL (B, W, R, RE);
    SUB (RE, U, RE);
    for (blk = dia->adj; blk; blk = blk->n)
    {
      double *R = blk->dia->R,
             *W = blk->W;

      NVADDMUL (RE, W, R, RE);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n) /* add external terms in both: local and global system modes */
    {
      double *R = ((CON*)blk->dia)->R, /* external reaction */
             *W = blk->W;

      NVADDMUL (RE, W, R, RE);
    }
#endif

    switch (con->kind)
    {
    case CONTACT:
    {
      double *V = dia->V,
	     rho = dia->rho,
	     gap = con->gap,
	     fri = con->mat.base->friction,
	     res = con->mat.base->restitution;

      int *map = dia->map;

      /* normal dashed velocity */
      if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
      else udash = ((MAX(gap, 0)/step) + U[2]);

      /* predict reaction */
      d [0] = R[0] - rho * U[0];
      d [1] = R[1] - rho * U[1];
      d [2] = R[2] - rho * udash;

      if (d [2] >= 0)  /* active normal set */
      {
	pull = 0;

#if MPI
	if (global)
	{
	  /* NOTE:
	   * in parallel global case 'map' references
	   * beginnings of compressed row blocks */

	  a [map[2]]   = W[2];
	  a [map[2]+1] = W[5];
	  a [map[2]+2] = W[8];
	}
	else
#endif
	{
	  /* NOTE:
	   * in sequential or local case 'map' references
	   * beginnings of compressed column blocks */

	  a [map[0]+2] = W[2];
	  a [map[1]+2] = W[5];
	  a [map[2]+2] = W[8];
	}

	b [2] = -udash - RE[2];

#if MPI
	if (global)
	{
	  OFFB *adj [2] = {dia->adj, dia->adjext};
	  for (int i = 0; i < 2; i ++)
	  for (blk = adj [i]; blk; blk = blk->n)
	  {
	    double *W = blk->W;
	    int *map = blk->map;

	    a [map[2]]   += W[2]; /* sum up off-diagonal blocks as pairs of W(i,j) are stored for two-body constraints */
	    a [map[2]+1] += W[5];
	    a [map[2]+2] += W[8];
	  }
	}
	else /* per-processor local system */
#endif
	for (blk = dia->adj; blk; blk = blk->n)
	{
	  double *W = blk->W;
	  int *map = blk->map;

	  a [map[0]+2] += W[2];
	  a [map[1]+2] += W[5];
	  a [map[2]+2] += W[8];
	}
      }
      else  /* inactive normal set */
      {
	pull = 1;

	a [map[2]+2] = W[8];

	b [2] = -R[2] * W[8];
      }

      /* tangential response */

      norm = sqrt (d[0]*d[0]+d[1]*d[1]); /* tangential force value */

      switch ((int) variant)
      {
	case NONSMOOTH_HSW:
	  lim = fri * MAX (0, d[2]);
	  break;
	case NONSMOOTH_HYBRID:
	  lim = fri * MAX (0, R[2]);
	  break;
	case FIXED_POINT:
	  lim = fri * MAX (0, RN(con));
	  break;
      }

      if (variant == NONSMOOTH_HSW && pull) goto ZERO_TANG; /* enforce AN = AT + IT */

      if (norm >= lim) /* frictional sliping */
      {
	double F [4], /* matrix associated with the derivative of an Euclidean norm in 2D */
	       M [4], H [4],  /* auxiliary metrices & vectors */
	       delta, alpha, beta, den, len, e; /* auxiliary scalars */

	/* active tangential set */

	if (lim > 0.0 && norm > 0.0) /* non-degenerate case */
	{
	  /* references below are linked to the paper "A primal-dual active set strategy
	   * for three-dimensional contact problems with Coulomb friction", S. Hueber, G. Stadler,
	   * B.I. Wohlmuth, SIAM Journal on Scientific Computing 30(2):572-596, 2007*/

	  len = sqrt (R[0]*R[0]+R[1]*R[1]);
	  den = MAX (lim, len) * norm;
	  e = lim / norm; /* (3.6) */
	  
#if 0
	  if (len == 0.0) beta = 1.0; /* commented in Section 3.3, although my conclusion is to set
					 'beta = 1' here rather than enforce homogenous tangential reaction */
	  else
	  {
	    /* formulae around (3.11) below */
	    alpha = (R[0]*d[0]+R[1]*d[1]) / (len*norm);
	    delta = MIN (len/lim, 1.0);
	    beta = (alpha < 0.0 ? 1.0 / (1.0 - alpha*delta) : 1.0); /* relaxation factor in case of direction change */
	  }
#else
	  alpha = delta = 0.0;
	  beta = 1.0;
#endif

	  /* (3.11) {or modified (3.6)} */
	  F [0] = (R[0]*d[0])/den;
	  F [1] = (R[1]*d[0])/den;
	  F [2] = (R[0]*d[1])/den;
	  F [3] = (R[1]*d[1])/den;

	  /* defined just above (3.10) */
	  M [0] = e * (1.0 - F[0]);
	  M [1] = - e * F[1];
	  M [2] = - e * F[2];
	  M [3] = e * (1.0 - F[3]);

	  if ((options & SYMMETRIZE) == 0) /* nonsymetric reduction */
	  {
	    /* inverse used in (3.10) */
	    H [0] = 1.0 - beta * M[0];
	    H [1] = - beta * M[1];
	    H [2] = - beta * M[2];
	    H [3] = 1.0 - beta * M[3];

	    /* scale */
	    M [0] *= W[0];
	    M [1] *= W[4];
	    M [2] *= W[0];
	    M [3] *= W[4];

	    H [0] *= W[0];
	    H [1] *= W[4];
	    H [2] *= W[0];
	    H [3] *= W[4];

#if MPI
	    if (global)
	    {
	      a [map[0]]   = H[0] + rho*(M[0]*W[0] + M[2]*W[1]);
	      a [map[1]]   = H[1] + rho*(M[1]*W[0] + M[3]*W[1]);
	      a [map[0]+1] = H[2] + rho*(M[0]*W[3] + M[2]*W[4]);
	      a [map[1]+1] = H[3] + rho*(M[1]*W[3] + M[3]*W[4]);
	    }
	    else
#endif
	    {
	      a [map[0]]   = H[0] + rho*(M[0]*W[0] + M[2]*W[1]);
	      a [map[0]+1] = H[1] + rho*(M[1]*W[0] + M[3]*W[1]);
	      a [map[1]]   = H[2] + rho*(M[0]*W[3] + M[2]*W[4]);
	      a [map[1]+1] = H[3] + rho*(M[1]*W[3] + M[3]*W[4]);
	    }

	    if (variant == FIXED_POINT)
	    {
#if MPI
	      if (global)
	      {
		a [map[0]+2] = rho*(M[0]*W[6] + M[2]*W[7]);
		a [map[1]+2] = rho*(M[1]*W[6] + M[3]*W[7]);
	      }
	      else
#endif
	      {
		a [map[2]]   = rho*(M[0]*W[6] + M[2]*W[7]);
		a [map[2]+1] = rho*(M[1]*W[6] + M[3]*W[7]);
	      }

	      b [0] = (fri*(d[0]/norm)*RN(con) - R[0]) * W[0] - rho*(M[0]*RE[0] + M[2]*RE[1]); 
	      b [1] = (fri*(d[1]/norm)*RN(con) - R[1]) * W[4] - rho*(M[1]*RE[0] + M[3]*RE[1]); 
	    }
	    else /* HSW, HYBRID */
	    {
#if MPI
	      if (global)
	      {
		a [map[0]+2] = rho*(M[0]*W[6] + M[2]*W[7]) - fri*(d[0]/norm) * W[0];
		a [map[1]+2] = rho*(M[1]*W[6] + M[3]*W[7]) - fri*(d[1]/norm) * W[4];
	      }
	      else
#endif
	      {
		a [map[2]]   = rho*(M[0]*W[6] + M[2]*W[7]) - fri*(d[0]/norm) * W[0];
		a [map[2]+1] = rho*(M[1]*W[6] + M[3]*W[7]) - fri*(d[1]/norm) * W[4];
	      }

	      b [0] = (fri*(d[0]/norm)*R[2] - R[0]) * W[0] - rho*(M[0]*RE[0] + M[2]*RE[1]);
	      b [1] = (fri*(d[1]/norm)*R[2] - R[1]) * W[4] - rho*(M[1]*RE[0] + M[3]*RE[1]);
	    }

#if MPI
	    if (global)
	    {
	      OFFB *adj [2] = {dia->adj, dia->adjext};
	      for (int i = 0; i < 2; i ++)
	      for (blk = adj [i]; blk; blk = blk->n)
	      {
		double *W = blk->W;
		int *map = blk->map;

		a [map[0]]   += rho*(M[0]*W[0] + M[2]*W[1]);
		a [map[1]]   += rho*(M[1]*W[0] + M[3]*W[1]);
		a [map[0]+1] += rho*(M[0]*W[3] + M[2]*W[4]);
		a [map[1]+1] += rho*(M[1]*W[3] + M[3]*W[4]);
		a [map[0]+2] += rho*(M[0]*W[6] + M[2]*W[7]); 
		a [map[1]+2] += rho*(M[1]*W[6] + M[3]*W[7]);  
	      }
	    }
	    else
#endif
	    for (blk = dia->adj; blk; blk = blk->n)
	    {
	      double *W = blk->W;
	      int *map = blk->map;

	      a [map[0]]   += rho*(M[0]*W[0] + M[2]*W[1]);
	      a [map[0]+1] += rho*(M[1]*W[0] + M[3]*W[1]);
	      a [map[1]]   += rho*(M[0]*W[3] + M[2]*W[4]);
	      a [map[1]+1] += rho*(M[1]*W[3] + M[3]*W[4]);
	      a [map[2]]   += rho*(M[0]*W[6] + M[2]*W[7]); 
	      a [map[2]+1] += rho*(M[1]*W[6] + M[3]*W[7]);  
	    }
	  }
	  else /* "symetric" reduction (symetric in the limit of iterates) */
	  {
	    double T [2];
	    int perm [2];

	    ASSERT (lapack_dgetrf (2, 2, M, 2, perm) == 0, ERR_MTX_MATRIX_INVERT);
	    ASSERT (lapack_dgetri (2, M, 2, perm, T, 2) == 0, ERR_MTX_MATRIX_INVERT);

	    H [0] = fri*(d[0]/norm);
	    H [1] = fri*(d[1]/norm);
	    T [0] = (M[0]*H[0] + M[2]*H[1])/rho;
	    T [1] = (M[1]*H[0] + M[3]*H[1])/rho;

	    H [0] = (M[0] - 1.0)/rho;
	    H [1] = M[1]/rho;
	    H [2] = M[2]/rho;
	    H [3] = (M[3] - 1.0)/rho;

#if MPI
	    if (global)
	    {
	      a [map[0]]   = W[0] + H[0];
	      a [map[1]]   = W[1] + H[1];
	      a [map[0]+1] = W[3] + H[2];
	      a [map[1]+1] = W[4] + H[3];
	    }
	    else
#endif
	    {
	      a [map[0]]   = W[0] + H[0];
	      a [map[0]+1] = W[1] + H[1];
	      a [map[1]]   = W[3] + H[2];
	      a [map[1]+1] = W[4] + H[3];
	    }

	    if (variant == FIXED_POINT)
	    {
#if MPI
	      if (global)
	      {
		a [map[0]+2] = W[6];
		a [map[1]+2] = W[7];
	      }
	      else
#endif
	      {
		a [map[2]]   = W[6];
		a [map[2]+1] = W[7];
	      }

	      b [0] = T[0]*RN(con) - (M[0]*R[0] + M[2]*R[1])/rho - RE[0];
	      b [1] = T[1]*RN(con) - (M[1]*R[0] + M[3]*R[1])/rho - RE[1];
	    }
	    else /* HSW, HYBRID */
	    {
#if MPI
	      if (global)
	      {
		a [map[0]+2] = W[6] - T[0];
		a [map[1]+2] = W[7] - T[1];
	      }
	      else
#endif
	      {
		a [map[2]]   = W[6] - T[0];
		a [map[2]+1] = W[7] - T[1];
	      }

	      b [0] = T[0]*R[2] - (M[0]*R[0] + M[2]*R[1])/rho - RE[0];
	      b [1] = T[1]*R[2] - (M[1]*R[0] + M[3]*R[1])/rho - RE[1];
	    }

#if MPI
	    if (global)
	    {
	      OFFB *adj [2] = {dia->adj, dia->adjext};
	      for (int i = 0; i < 2; i ++)
	      for (blk = adj [i]; blk; blk = blk->n)
	      {
		double *W = blk->W;
		int *map = blk->map;

		a [map[0]]   += W[0];
		a [map[1]]   += W[1];
		a [map[0]+1] += W[3];
		a [map[1]+1] += W[4];
		a [map[0]+2] += W[6];
		a [map[1]+2] += W[7];
	      }
	    }
	    else
#endif
	    for (blk = dia->adj; blk; blk = blk->n)
	    {
	      double *W = blk->W;
	      int *map = blk->map;

	      a [map[0]]   += W[0];
	      a [map[0]+1] += W[1];
	      a [map[1]]   += W[3];
	      a [map[1]+1] += W[4];
	      a [map[2]]   += W[6];
	      a [map[2]+1] += W[7];
	    }
	  }
	}
	else /* degenerate case => enforce homogenous tangential tractions.
		this is discussed in the Section 3.3 of the HSW paper */
	{
  ZERO_TANG:
	  a [map[0]]   = W[0];
	  a [map[1]+1] = W[4];

	  b [0] = -R[0] * W[0];
	  b [1] = -R[1] * W[4];
	}
      }
      else /* frictional sticking: inactive tangential set */
      {

#if MPI
	if (global)
	{
	  a [map[0]]   = W[0];
	  a [map[1]]   = W[1];
	  a [map[0]+1] = W[3];
	  a [map[1]+1] = W[4];
	}
	else
#endif
	{
	  a [map[0]]   = W[0];
	  a [map[0]+1] = W[1];
	  a [map[1]]   = W[3];
	  a [map[1]+1] = W[4];
	}

	if (variant == NONSMOOTH_HYBRID || variant == FIXED_POINT)
	{
#if MPI
	  if (global)
	  {
	    a [map[0]+2] = W[6];
	    a [map[1]+2] = W[7];
	  }
	  else
#endif
	  {
	    a [map[2]]   = W[6];
	    a [map[2]+1] = W[7];
	  }

	  if (dynamic)
	  {
	    b [0] = -U[0] - V[0] - RE[0]; /* -V = W(R + dR) + B; U = WR + B; -U - V = W dR */
	    b [1] = -U[1] - V[1] - RE[1];
	  }
	  else
	  {
	    b [0] = -U[0] - RE[0];
	    b [1] = -U[1] - RE[1];
	  }
	}
	else /* HSW */
	{
	  /* FIXME: this is quasi-statics; work out dynamics; it actually works pretty well in the dynamic case */
#if MPI
	  if (global)
	  {
	    a [map[0]+2] = W[6]+U[0]/d[2];
	    a [map[1]+2] = W[7]+U[1]/d[2];
	  }
	  else
#endif
	  {
	    a [map[2]]   = W[6]+U[0]/d[2];
	    a [map[2]+1] = W[7]+U[1]/d[2];
	  }

	  b [0] = -(1.0 + rho*udash/d[2])*U[0] - RE[0];
	  b [1] = -(1.0 + rho*udash/d[2])*U[1] - RE[1];
	}

#if MPI
	if (global)
	{
	  OFFB *adj [2] = {dia->adj, dia->adjext};
	  for (int i = 0; i < 2; i ++)
	  for (blk = adj [i]; blk; blk = blk->n)
	  {
	    double *W = blk->W;
	    int *map = blk->map;

	    a [map[0]]   += W[0];
	    a [map[1]]   += W[1];
	    a [map[0]+1] += W[3];
	    a [map[1]+1] += W[4];
	    a [map[0]+2] += W[6];
	    a [map[1]+2] += W[7];
	  }
	}
	else
#endif
	for (blk = dia->adj; blk; blk = blk->n)
	{
	  double *W = blk->W;
	  int *map = blk->map;

	  a [map[0]]   += W[0];
	  a [map[0]+1] += W[1];
	  a [map[1]]   += W[3];
	  a [map[1]+1] += W[4];
	  a [map[2]]   += W[6];
	  a [map[2]+1] += W[7];
	}
      }
    }
    break;
    case FIXPNT:
    {
      double *V = dia->V;

      int *map = dia->map;

#if MPI
      if (global)
      {
	a [map[0]]   = W[0];
	a [map[1]]   = W[1];
	a [map[2]]   = W[2];
	a [map[0]+1] = W[3];
	a [map[1]+1] = W[4];
	a [map[2]+1] = W[5];
	a [map[0]+2] = W[6];
	a [map[1]+2] = W[7];
	a [map[2]+2] = W[8];
      }
      else
#endif
      {
	a [map[0]]   = W[0];
	a [map[0]+1] = W[1];
	a [map[0]+2] = W[2];
	a [map[1]]   = W[3];
	a [map[1]+1] = W[4];
	a [map[1]+2] = W[5];
	a [map[2]]   = W[6];
	a [map[2]+1] = W[7];
	a [map[2]+2] = W[8];
      }

      if (dynamic)
      {
	b [0] = -U[0] - V[0] - RE[0];
	b [1] = -U[1] - V[1] - RE[1];
	b [2] = -U[2] - V[2] - RE[2];
      }
      else
      {
	b [0] = -U[0] - RE[0];
	b [1] = -U[1] - RE[1];
	b [2] = -U[2] - RE[2];
      }

#if MPI
      if (global)
      {
	OFFB *adj [2] = {dia->adj, dia->adjext};
	for (int i = 0; i < 2; i ++)
	for (blk = adj [i]; blk; blk = blk->n)
	{
	  double *W = blk->W;
	  int *map = blk->map;

	  a [map[0]]   += W[0];
	  a [map[1]]   += W[1];
	  a [map[2]]   += W[2];
	  a [map[0]+1] += W[3];
	  a [map[1]+1] += W[4];
	  a [map[2]+1] += W[5];
	  a [map[0]+2] += W[6];
	  a [map[1]+2] += W[7];
	  a [map[2]+2] += W[8];
	}
      }
      else
#endif
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *W = blk->W;
	int *map = blk->map;

	a [map[0]]   += W[0];
	a [map[0]+1] += W[1];
	a [map[0]+2] += W[2];
	a [map[1]]   += W[3];
	a [map[1]+1] += W[4];
	a [map[1]+2] += W[5];
	a [map[2]]   += W[6];
	a [map[2]+1] += W[7];
	a [map[2]+2] += W[8];
      }
    }
    break;
    case FIXDIR:
    case VELODIR:
    case RIGLNK:
    {
      double *V = dia->V;

      int *map = dia->map;

      a [map[0]]   = W[0]; /* set diagonal to diag (W) */
      a [map[1]+1] = W[4];
      a [map[2]+2] = W[8];

      b [0] = -R[0] * W[0]; /* keep RT zero: W (-R) = W dR */
      b [1] = -R[1] * W[4];

      if (dynamic)
      {
	if (con->kind == VELODIR) b [2] = -U[2] + VELODIR(con->Z) - RE[2]; /* V(t+h) = W(R+dR) + B; U = WR + B; -U+V(t+h) = WdR */
	else b [2] = -U[2] - V[2] - RE[2]; /* FIXDIR, RIGLNK */
      }
      else
      {
	if (con->kind == VELODIR) b [2] = -U[2] + VELODIR(con->Z) - RE[2];
	else if (con->kind == FIXDIR) b [2] = -U[2] - RE[2]; 
	else /* RIGLNK: see doc/notes.lyx for explanation */
	{
	  double C [3];

	  ADD (U, RE, C);

	  double d = RIGLNK_LEN (con->Z),
		 e = step * step * DOT2 (C, C);

	  if (d*d > e)
	  {
	    b[2] = -C[2] + (-d + sqrt (d*d - e)) / step;
	  }
	  else /* not possible to satisfy: do your best */
	  {
	    b[2] = -C[2] - d / step;
	  }
	}
      }

#if MPI
      if (global)
      {
	OFFB *adj [2] = {dia->adj, dia->adjext};
	for (int i = 0; i < 2; i ++)
	for (blk = adj [i]; blk; blk = blk->n)
	{
	  double *W = blk->W;
	  int *map = blk->map;

	  a [map[2]]   += W[2]; /* maintain tangentail coupling: W(3,:) */
	  a [map[2]+1] += W[5];
	  a [map[2]+2] += W[8];
	}
      }
      else
#endif
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *W = blk->W;
	int *map = blk->map;

	a [map[0]+2] += W[2];
	a [map[1]+2] += W[5];
	a [map[2]+2] += W[8];
      }
    }
    break;
    }
  }
}

/* real normal to friction cone */
inline static void real_n (double *R, double fri, double *n)
{
  double dot, len;

  dot = DOT2(R, R);
  len = sqrt(dot);

  if (len <= fri * R[2])
  {
    SET (n, 0.0);
  }
  else if (fri * len + R[2] < 0.0)
  {
    dot += R[2]*R[2];
    len = sqrt (dot);
    DIV (R, len, n);
  }
  else
  {
    dot = 1.0 / sqrt (1.0 + fri*fri);
    DIV2 (R, len, n);
    n [2] = -fri;
    SCALE (n, dot);
  }
}

/* complex normal to friction cone */
inline static void complex_n (double complex *R, double complex fri, double complex *n)
{
  double complex dot, len;

  dot = DOT2(R, R);
  len = csqrt(dot);

  if (creal (len) <= creal (fri * R[2]))
  {
    SET (n, 0.0);
  }
  else if (creal (fri * len + R[2]) < 0.0)
  {
    dot += R[2]*R[2];
    len = csqrt (dot);
    DIV (R, len, n);
  }
  else
  {
    dot = 1.0 / csqrt (1.0 + fri*fri);
    DIV2 (R, len, n);
    n [2] = -fri;
    SCALE (n, dot);
  }
}

#define SMOOTHING 1 /* FIXME */
#define SMOOTHING_EPSILON 1E-10 /* FIXME: get rid of */
#define DERIVATION_STEP 1E-15

/* real normal ray to friction cone */
inline static void real_m (CON *con, short smooth, double *R, double eps, double *m)
{
  double n [3], fun, fri;

  if (con->kind == CONTACT)
  {
    fri = con->mat.base->friction;

    real_n (R, fri, n);
    fun = DOT (R, n);

#if SMOOTHING == 1
    if (smooth && fun >= 0.0 && fun <= eps)
    {
      fun = ((2.0/eps) - (1.0/(eps*eps))*fun)*(fun*fun);
    }
#elif SMOOTHING == 2
    if (smooth)
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
#endif

    MUL (n, fun, m)
  }
  else
  {
    SET (m, 0.0);
  }
}

/* complex normal ray to friction cone */
inline static void complex_m (CON *con, short smooth, double complex img, double complex *R, double complex eps, double complex *m)
{
  double complex n [3], fun, fri;

  if (con->kind == CONTACT)
  {
#if MPI
    if (con->state & CON_EXTERNAL) fri = con->Y [3];
    else
#endif
    fri = con->mat.base->friction;

    complex_n (R, fri, n);
    fun = DOT (R, n);

#if SMOOTHING == 1
    if (smooth && creal (fun) >= 0.0 && creal (fun) <= creal (eps))
    {
      fun = ((2.0/eps) - (1.0/(eps*eps))*fun)*(fun*fun);
    }
#elif SMOOTHING == 2
    if (smooth)
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
#endif

    MUL (n, fun, m)
  }
  else
  {
    m [0] = 0.0 + 0.0 * img;
    m [1] = 0.0 + 0.0 * img;
    m [2] = 0.0 + 0.0 * img;
  }
}

/* real variational merit function (note that for contacts R can be outside of the friction cone here) */
inline static void real_H (CON *con, double step, short dynamic, short smooth, double *U, double *R, double epsilon, double *H)
{
  DIAB *dia = con->dia;
  double *V = dia->V;

  switch (con->kind)
  {
  case CONTACT:
  {
    double res = con->mat.base->restitution,
	   fri = con->mat.base->friction,
	   gap = con->gap,
           udash, m [3];

    if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
    else udash = ((MAX(gap, 0)/step) + U[2]);
 
    real_m (con, smooth, R, epsilon, m);

    H [0] = U[0] + m[0];
    H [1] = U[1] + m[1];
    H [2] = udash + fri * LEN2(U) + m[2];
  }
  break;
  case FIXPNT:
  {
    if (dynamic) { ADD (V, U, H); }
    else { COPY (U, H); }
  }
  break;
  case FIXDIR:
  {
    H [0] = R[0];
    H [1] = R[1];
    if (dynamic) H [2] = U[2] + V[2];
    else H [2] = U[2];
  }
  break;
  case VELODIR:
  {
    H [0] = R[0];
    H [1] = R[1];
    H [2] = U[2] - VELODIR(con->Z); /* U(t+h) = function of time stored at (t+h) in con->Z */
  }
  break;
  case RIGLNK:
  {
    H [0] = R[0];
    H [1] = R[1];
    if (dynamic) H[2] = U[2] + V[2];
    else /* RIGLNK: see doc/notes.lyx for explanation */
    {
      double d = RIGLNK_LEN (con->Z),
	     e = step * step * DOT2 (U, U);

      if (d*d > e)
      {
	H[2] = U[2] + (d - sqrt (d*d - e)) / step;
      }
      else /* not possible to satisfy: do your best */
      {
	H[2] = U[2] + d / step;
      }
    }
  }
  break;
  };
}

/* complex variational auxiliary function */
inline static void complex_F (CON *con, double step, short dynamic, double complex *U, double complex *R, double complex *H)
{
  DIAB *dia = con->dia;
  double *V = dia->V;

  switch (con->kind)
  {
  case CONTACT:
  {
    double res = con->mat.base->restitution,
	   fri = con->mat.base->friction,
	   gap = con->gap;

    double complex udash;

    if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
    else udash = ((MAX(gap, 0)/step) + U[2]);
 
    H [0] = U[0];
    H [1] = U[1];
    H [2] = udash + fri * csqrt (DOT2(U, U));
  }
  break;
  case FIXPNT:
  {
    if (dynamic) { ADD (V, U, H); }
    else { COPY (U, H); }
  }
  break;
  case FIXDIR:
  {
    H [0] = R ? R[0] : 0.0;
    H [1] = R ? R[1] : 0.0;
    if (dynamic) H [2] = U[2] + V[2];
    else H [2] = U[2];
  }
  break;
  case VELODIR:
  {
    H [0] = R ? R[0] : 0.0;
    H [1] = R ? R[1] : 0.0;
    H [2] = U[2] - VELODIR(con->Z); /* U(t+h) = function of time stored at (t+h) in con->Z */
  }
  break;
  case RIGLNK:
  {
    H [0] = R ? R[0] : 0.0;
    H [1] = R ? R[1] : 0.0;
    if (dynamic) H[2] = U[2] + V[2];
    else /* RIGLNK: see doc/notes.lyx for explanation */
    {
      double complex d = RIGLNK_LEN (con->Z),
	             e = step * step * DOT2 (U, U);

      if (creal (d*d) > creal (e))
      {
	H[2] = U[2] + (d - csqrt (d*d - e)) / step;
      }
      else /* not possible to satisfy: do your best */
      {
	H[2] = U[2] + d / step;
      }
    }
  }
  break;
  }
}

/* forward change of coordinates: R = R + m(R) */
static void variational_coord_change_forward (LINSYS *sys)
{
  CON *con;

  /* R_i = R_i - m(R_i) */
  for (con = sys->ldy->dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT)
    {
      double *R = con->R,
	     *Z = con->Z;

      COPY (Z, R);
    }
  }
}

/* backward change of coordinates: R = R - m(R) */
static void variational_coord_change_backward (LINSYS *sys)
{
  CON *con;
  short smooth = sys->variant & SMOOTHED_VARIATIONAL;
  double epsilon = SMOOTHING_EPSILON; /* FIXME */

  /* R_i = R_i - m(R_i) */
  for (con = sys->ldy->dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT)
    {
      double *R = con->R,
	     *Z = con->Z,
	      m [3];

      COPY (R, Z);
      real_m (con, smooth, R, epsilon, m); /* FIXME: epsilon-adhesion allowed */
      SUB (R, m, R);
    }
  }
}

/* update linear system NT_NONSMOOTH_VARIATIONAL, NT_SMOOTHED_VARIATIONAL variants */
static void system_update_VARIATIONAL (LINVAR variant, LINOPT options, LOCDYN *ldy, double *a, double *b)
{
  short smooth, global, dynamic;
  double epsilon, h, step;
  double complex img;
  OFFB *blk;
  DIAB *dia;
  CON *con;
  DOM *dom;

  smooth = variant & SMOOTHED_VARIATIONAL;
  global = !(options & LOCAL_SYSTEM);
  epsilon = SMOOTHING_EPSILON; /* FIXME: make into a user parameter */
  h = DERIVATION_STEP; /* FIXME: maker into a user parameter */ 
  dom = ldy->dom;
  img = csqrt (-1.0); /* imaginary i */
  step = dom->step;
  dynamic = dom->dynamic;

  /* S_i = R_i - m(R_i) */
  for (con = dom->con; con; con = con->next)
  {
    double *R = con->R,
	   *S = con->Y, /* use Y to store S */
	    m [3];

    real_m (con, smooth, R, epsilon, m);
    SUB (R, m, S);

#if MPI
     /* friction needs to be sent to external constraints //see complex_m()// */
    if (global && con->kind == CONTACT) S [3] = con->mat.base->friction;
#endif
  }

#if MPI
  if (global)
  {
    DOM_Update_External_Reactions (ldy->dom, 0); /* update external CON->R */

    DOM_Update_External_Y (dom, 4); /* update external CON->Y with (S, friction) */
  }
#endif

  /* U_i = SUM_j {W_ij * S_j} + B_i */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    double *W = dia->W,
	   *B = dia->B,
	   *S = dia->con->Y,
	   *U = dia->U;

    NVADDMUL (B, W, S, U);
    for (blk = dia->adj; blk; blk = blk->n)
    {
      double *W = blk->W,
	     *S = blk->dia->con->Y;

      NVADDMUL (U, W, S, U);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      double *W = blk->W,
	     *S = global ? CON(blk->dia)->Y : CON(blk->dia)->R;

      NVADDMUL (U, W, S, U);
    }
#endif
  }

  /* d H / d R */
  for (con = dom->con; con; con = con->next, b += 3)
  {
    dia = con->dia;

    double *W = dia->W,
	   *S = con->Y,
	   *R = dia->R,
	   *U = dia->U;

    int *map = dia->map, k;

    double complex cR [3],
	           cU [3],
	           cm [3],
		   cS [3],
		   cQ [3],
		   cF [3];

    for (k = 0; k < 3; k ++)
    {
      COPY (R, cR);
      cR [k] += h * img;

      complex_m (con, smooth, img, cR, epsilon, cm);
      SUB (cR, cm, cS);
      SUB (cS, S, cQ);
      NVADDMUL (U, W, cQ, cU);
      complex_F (con, step, dynamic, cU, cR, cF);
      if (con->kind == CONTACT) { ADD (cF, cm, cF); }

#if MPI
      if (global)
      {
	a  [map[0] + k] = cimag (cF [0]) / h;
	a  [map[1] + k] = cimag (cF [1]) / h;
	a  [map[2] + k] = cimag (cF [2]) / h;
      }
      else
#endif
      {
	a  [map[k] + 0] = cimag (cF [0]) / h;
	a  [map[k] + 1] = cimag (cF [1]) / h;
	a  [map[k] + 2] = cimag (cF [2]) / h;
      }
    }

    real_H (con, step, dynamic, smooth, U, R, epsilon, b);
    SCALE (b, -1.0);

#if MPI
    if (global)
    {
      OFFB *adj [2] = {dia->adj, dia->adjext};
      for (int i = 0; i < 2; i ++)
      for (blk = adj [i]; blk; blk = blk->n)
      {
	double *W = blk->W, *R, *S;

	int *map = blk->map;

	CON *adc; /* adjacent constraint */

	if (adj[i] == dia->adj) adc = blk->dia->con;
	else adc = CON(blk->dia);

	R = adc->R;
	S = adc->Y;

	for (k = 0; k < 3; k ++)
	{
	  COPY (R, cR);
	  cR [k] += h * img;

	  complex_m (adc, smooth, img, cR, epsilon, cm);
	  SUB (cR, cm, cS);
	  SUB (cS, S, cQ);
	  NVADDMUL (U, W, cQ, cU);
	  complex_F (con, step, dynamic, cU, NULL, cF);

	  a  [map[0] + k] += cimag (cF [0]) / h;
	  a  [map[1] + k] += cimag (cF [1]) / h;
	  a  [map[2] + k] += cimag (cF [2]) / h;
	}
      }
    }
    else
#endif
    for (blk = dia->adj; blk; blk = blk->n)
    {
      CON *adc = blk->dia->con;

      double *W = blk->W,
	     *R = adc->R,
	     *S = adc->Y;

      int *map = blk->map;

      for (k = 0; k < 3; k ++)
      {
	COPY (R, cR);
	cR [k] += h * img;

	complex_m (adc, smooth, img, cR, epsilon, cm);
	SUB (cR, cm, cS);
	SUB (cS, S, cQ);
	NVADDMUL (U, W, cQ, cU);
	complex_F (con, step, dynamic, cU, NULL, cF);

	a  [map[k] + 0] += cimag (cF [0]) / h;
	a  [map[k] + 1] += cimag (cF [1]) / h;
	a  [map[k] + 2] += cimag (cF [2]) / h;
      }
    }
  }
}

/* update linear system at current R and U */
void LINSYS_Update (LINSYS *sys)
{
  for (double *a = sys->a, *b = a + sys->nnz; a < b; a ++) *a = 0.0; /* zero system matrix */

  switch (sys->variant)
  {
    case NONSMOOTH_HSW:
    case NONSMOOTH_HYBRID:
    case FIXED_POINT:
      system_update_HSW_HYBRID_FIXED (sys->variant, sys->options, sys->ldy, sys->a, sys->b);
      break;
    case NONSMOOTH_VARIATIONAL:
    case SMOOTHED_VARIATIONAL:
      system_update_VARIATIONAL (sys->variant, sys->options, sys->ldy, sys->a, sys->b);
      break;
  }
}

/* solve linear system for reaction increments DR */
void LINSYS_Solve (LINSYS *sys, double accuracy, int maxiter)
{
  DOM *dom = sys->ldy->dom;

#if 0
  /* improve conditioning */
  {
    int j, k, q, *i, *p;
    double *a, delta;

    delta = 0.1 * accuracy;
    accuracy = accuracy - delta;

#if MPI
    if (!(sys->options & LOCAL_SYSTEM))
    {
      i = sys->cols;
      p = sys->rows;
      q = sys->hyp->ilower;
    }
    else
#endif
    {
      i = sys->rows;
      p = sys->cols;
      q = 0;
    }

    for (j = 0, a = sys->a; j < sys->dim; j ++)
    {
      for (k = p [j]; k < p [j+1]; k ++, a ++, i ++)
      {
	if ((q+j) == *i)
	{
	  *a += delta;
	}
      }
    }
  }
#endif

#if MPI
  if (sys->options & LOCAL_SYSTEM)
  {
#endif

#if SPQR
  {
    /* x = A\b */
    cholmod_dense *x = SuiteSparseQR_C_backslash_default (sys->spqr_a, sys->spqr_b, &sys->spqr_c);
    blas_dcopy (sys->dim, x->x, 1, sys->x, 1);
    cholmod_l_free_dense (&x, &sys->spqr_c) ;

    /* stats */
    sys->resnorm = 0.0;
    sys->iters = 0;
  }
#elif UMFPACK
  {
    /* x = A\b */
    void *numeric;
    ASSERT (umfpack_di_numeric (sys->cols, sys->rows, sys->a, sys->symbolic, &numeric, NULL, NULL) == UMFPACK_OK, ERR_UMFPACK_SOLVE);
    ASSERT (umfpack_di_solve (UMFPACK_A, sys->cols, sys->rows, sys->a, sys->x, sys->b, numeric, NULL, NULL) == UMFPACK_OK, ERR_UMFPACK_SOLVE);
    umfpack_di_free_numeric (&numeric);

    /* stats */
    sys->resnorm = 0.0;
    sys->iters = 0;
  }
#else
  {
    /* x = A\b */
    LSS_Set (sys->lss, LSS_ABSOLUTE_ACCURACY, accuracy);
    LSS_Set (sys->lss, LSS_ITERATIONS_BOUND, maxiter);
    LSS_Set (sys->lss, LSS_RELATIVE_ACCURACY, 1.0);
    LSS_Set (sys->lss, LSS_NORMALIZE_ROWS, 0);
#if 0
    LSS_Set (sys->lss, LSS_PRECONDITIONER, 3);
    LSS_Set (sys->lss, LSS_DECIMATION, 16);
    LSS_Set (sys->lss, LSS_SMOOTHING_STEPS, 5);
#endif
    if (LSS_Solve (sys->lss, sys->a, sys->x, sys->b) != LSSERR_NONE)
    {
      fprintf (stderr, "WARNING: LSS failed with message: %s\n", LSS_Errmsg (sys->lss));
    }

    /* stats */
    sys->resnorm = LSS_Get (sys->lss, LSS_ABSOLUTE_ERROR);
    sys->iters = LSS_Get (sys->lss, LSS_ITERATIONS);
  }
#endif

#if MPI
  }
  else /* parallel solve */
  {
    HYPRE_DATA *hyp = sys->hyp;

    if (hyp->comm != MPI_COMM_NULL) /* skip empty constraints case */
    {
      /* A */
      HYPRE_IJMatrixInitialize (hyp->A);
      HYPRE_IJMatrixSetValues (hyp->A, sys->dim, hyp->ncols, hyp->rows, sys->cols, sys->a);
      HYPRE_IJMatrixAssemble (hyp->A);

      /* b */
      HYPRE_IJVectorInitialize (hyp->B);
      HYPRE_IJVectorSetValues (hyp->B, sys->dim, hyp->rows, sys->b);
      HYPRE_IJVectorAssemble (hyp->B);

      /* Flexible GMRES with  AMG Preconditioner */
      {
	HYPRE_Solver solver;

	/* Create solver */
	HYPRE_ParCSRFlexGMRESCreate(hyp->comm, &solver);

	/* Set some parameters (See Reference Manual for more parameters) */
	HYPRE_FlexGMRESSetAbsoluteTol (solver, accuracy); /* conv. absolute tolerance */
	HYPRE_FlexGMRESSetMaxIter (solver, maxiter); /* max iterations */
	HYPRE_FlexGMRESSetTol (solver, 0.0); /* conv. tolerance */
	HYPRE_FlexGMRESSetKDim (solver, 32); /* restart */
	HYPRE_FlexGMRESSetPrintLevel (solver, 0); /* print solve info */
	HYPRE_FlexGMRESSetLogging (solver, 1); /* needed to get run info later */

#define PRECND 0
#define AMG 1
#define EUCLID 0
#define PARASAILS 0

#if PRECND
	HYPRE_Solver precond;

#if AMG
	/* Now set up the AMG preconditioner and specify any parameters */
	HYPRE_BoomerAMGCreate (&precond);
	HYPRE_BoomerAMGSetPrintLevel (precond, 0); /* print amg solution info */
	HYPRE_BoomerAMGSetCoarsenType (precond, 6);
	HYPRE_BoomerAMGSetRelaxType (precond, 6); /* Sym G.S./Jacobi hybrid */ 
	HYPRE_BoomerAMGSetNumSweeps (precond, 2);
	HYPRE_BoomerAMGSetTol (precond, 0.0); /* conv. tolerance zero */
	HYPRE_BoomerAMGSetMaxIter (precond, 1); /* do only one iteration! */

	/* Set the FlexGMRES preconditioner */
	HYPRE_FlexGMRESSetPrecond (solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
			    (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
#elif EUCLID
	HYPRE_EuclidCreate (hyp->comm, &precond);
        HYPRE_EuclidSetLevel(precond, 3);

        HYPRE_FlexGMRESSetPrecond (solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
	    (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, precond);
#elif PARASAILS 
      HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &precond);

      HYPRE_ParaSailsSetParams (precond, 0.1, 1);
      HYPRE_ParaSailsSetFilter (precond, 0.1);
      HYPRE_ParaSailsSetSym (precond, 0);
      HYPRE_ParaSailsSetLogging (precond, 0);

      HYPRE_FlexGMRESSetPrecond (solver, (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSolve,
                          (HYPRE_PtrToSolverFcn) HYPRE_ParaSailsSetup, precond);
#endif
#endif

	/* Now setup and solve! */
	HYPRE_ParCSRFlexGMRESSetup (solver, hyp->AP, hyp->BP, hyp->XP);
	HYPRE_ParCSRFlexGMRESSolve (solver, hyp->AP, hyp->BP, hyp->XP);

	/* Run info - needed logging turned on */
	HYPRE_FlexGMRESGetNumIterations (solver, &sys->iters);
	HYPRE_FlexGMRESGetFinalRelativeResidualNorm (solver, &sys->resnorm);

	/* Destory solver and preconditioner */
	HYPRE_ParCSRFlexGMRESDestroy (solver);
#if PRECND
#if AMG
	HYPRE_BoomerAMGDestroy (precond);
#elif EUCLID
	HYPRE_EuclidDestroy (precond);
#elif PARASAILS
      HYPRE_ParaSailsDestroy (precond);
#endif
#endif
      }

      /* x */
      HYPRE_IJVectorGetValues (hyp->X, sys->dim, hyp->rows, sys->x);
    }
  }
#endif

  double *x;
  CON *con;

  /* write DR */
  for (con = dom->con, x = sys->x; con; con = con->next, x += 3)
  {
    double *DR = DR(con);

    COPY (x, DR);
  }

#if MPI
  /* update external reactions increments in case of
   * NONSMOOTH_HSW, NONSMOOTH_HYBRID and FIXED_POINT */
  if (sys->variant <= FIXED_POINT && LINSYS_Global (sys))
  {
    for (con = dom->con; con; con = con->next)
    {
      double *DR = DR(con),
	     *Y  = con->Y;

      COPY (DR, Y);
    }

    /* send reaction increments to external
     * constraint CON->Y members: see (***) below */
    DOM_Update_External_Y (dom, 3);
  }
#endif

  /* compute DU in case of NONSMOOTH_HSW,
   * NONSMOOTH_HYBRID or FIXED_POINT variants */
  if (sys->variant <= FIXED_POINT)
  {
    for (con = dom->con; con; con = con->next)
    {
      DIAB *dia = con->dia;

      double *DR = DR(con),
	     *DU = DU(con),
	     *W = dia->W;

      NVMUL (W, DR, DU);

      for (OFFB *blk = dia->adj; blk; blk = blk->n)
      {
	double *DR = DR(blk->dia->con),
	       *W = blk->W;

	NVADDMUL (DU, W, DR, DU);
      }
#if MPI
      if (LINSYS_Global (sys))
      {
	for (OFFB *blk = dia->adjext; blk; blk = blk->n)
	{
	  double *DR = CON(blk->dia)->Y, /* (***) external reaction increment stored in CON->Y temporarily */
		 *W = blk->W;

	  NVADDMUL (DU, W, DR, DU);
	}
      }
#endif
    }
  }
}

/* compute merit function at (R + alpha * DR) */
double LINSYS_Merit (LINSYS *sys, double alpha)
{
  short dynamic, smooth, global;
  double value, step, epsilon;
  LINVAR variant;
  LOCDYN *ldy;
  DOM *dom;
  CON *con;

  value = 0.0;
  ldy = sys->ldy;
  dom = ldy->dom;
  step = dom->step;
  dynamic = dom->dynamic;
  variant = sys->variant;
  smooth = variant & SMOOTHED_VARIATIONAL;
  epsilon = SMOOTHING_EPSILON; /* FIXME: make into a user parameter */
  global = LINSYS_Global (sys);

  if (sys->variant == SMOOTHED_VARIATIONAL ||
      sys->variant == NONSMOOTH_VARIATIONAL)
  {
    for (con = dom->con; con; con = con->next)
    {
      double *R0 = con->R,
	     *DR = DR(con),
	     *S  = con->Y,
	      R [3], m [3];

      ADDMUL (R0, alpha, DR, R);
      real_m (con, smooth, R, epsilon, m);
      SUB (R, m, S);
    }

#if MPI
    if (global) DOM_Update_External_Y (dom, 3);
#endif
  }

  for (con = dom->con; con; con = con->next)
  {
    double G [3];
    DIAB *dia;
    OFFB *blk;

    dia = con->dia;

    if (variant <= FIXED_POINT)
    {
      double R [3], U [3],
	    *R0 = dia->R,
	    *U0 = dia->U;

      if (alpha > 0.0)
      {
	double *RE = RE(con),
	       *DR = DR(con),
	       *DU = DU(con);

	ADDMUL (R0, alpha, DR, R);
	ADDMUL (U0, alpha, DU, U);
	ADD (U, RE, U);
      }
      else
      {
	COPY (R0, R);
	COPY (U0, U);
      }

      if (con->kind == CONTACT)
      {
	double rho = dia->rho,
	       gap = con->gap,
	       fri = con->mat.base->friction,
	       res = con->mat.base->restitution,
	      *V = dia->V,
	       d [3],
	       bound,
	       udash,
	       norm;

	if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
	else udash = ((MAX(gap, 0)/step) + U[2]);

	d [0] = R[0] - rho * U[0];
	d [1] = R[1] - rho * U[1];
	d [2] = R[2] - rho * udash;

	norm = sqrt (d[0]*d[0]+d[1]*d[1]);

	switch ((int) variant)
	{
	  case NONSMOOTH_HSW:
	    bound = fri * MAX (0, d[2]);
	    break;
	  case NONSMOOTH_HYBRID:
	    bound = fri * MAX (0, R[2]);
	    break;
	  case FIXED_POINT:
	    bound = fri * MAX (0, RN(con));
	    break;
	}

	G [0] = MAX (bound, norm)*R[0] - bound*d[0];
	G [1] = MAX (bound, norm)*R[1] - bound*d[1];
	G [2] = R [2] - MAX (0, d[2]);
      }
      else /* non-contact */
      {
	real_H (con, step, dynamic, 0, U, R, 0.0, G); /* objective functions do not differ in this case */
      }
    }
    else /* SMOOTHED_VARIATIONAL, NONSMOOTH_VARIATIONAL */
    {
      double *S = con->Y,
	     *W = dia->W,
	     *B = dia->B,
	     *DR = DR(con),
	     *R0 = con->R,
	      U[3], R[3];

      NVADDMUL (B, W, S, U);
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *S = blk->dia->con->Y,
	       *W = blk->W;

	NVADDMUL (U, W, S, U);
      }
#if MPI
      for (blk = dia->adjext; blk; blk = blk->n)
      {
	double *W = blk->W,
	       *S = global ? CON(blk->dia)->Y : CON(blk->dia)->R;

	NVADDMUL (U, W, S, U);
      }
#endif
      ADDMUL (R0, alpha, DR, R);

      real_H (con, step, dynamic, smooth, U, R, epsilon, G);
    }

    value += DOT (G,G);
  }
  value *= 0.5;

#if MPI
  if (global)
  {
    double val_i = value;
    MPI_Allreduce (&val_i, &value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
#endif

  return value;
}

/* advance solution R = R + alpha * DR; return |DR|/|R| */ 
double LINSYS_Advance (LINSYS *sys, double alpha)
{
  double errup, errlo;
  CON *con;

  errup = errlo = 0.0;

  if (sys->variant <= FIXED_POINT)
  {
    for (con = sys->ldy->dom->con; con; con = con->next)
    {
      double *RE = RE(con),
	     *DR = DR(con),
	     *DU = DU(con),
	     *R = con->R,
	     *U = con->U;

      ADDMUL (R, alpha, DR, R);
      ADDMUL (U, alpha, DU, U);
      ADD (U, RE, U);

      errup += DOT(DR, DR); /* no alpha scaling here */
      errlo += DOT(R, R);
    }
  }
  else
  {
    for (con = sys->ldy->dom->con; con; con = con->next)
    {
      double *DR = DR(con),
	     *R = con->R;

      ADDMUL (R, alpha, DR, R);

      errup += DOT(DR, DR); /* no alpha scaling here */
      errlo += DOT(R, R);
    }
  }

#if MPI
  if (LINSYS_Global (sys)) /* sum up error */
  {
    double errloc [2] = {errup, errlo}, errsum [2];
    MPI_Allreduce (errloc, errsum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    errup = errsum [0], errlo = errsum [1];
  }
#endif

  return sqrt (errup) / sqrt (MAX (errlo, 1.0));
}

#if 0
/* qsort comparison */
static int dblcmp (const void *a, const void *b)
{
  if (*(double*)a < *(double*)b) return -1;
  else if (*(double*)a == *(double*)b) return 0;
  else return 1;
}

/* output system matrix in such a way that matrices from serial and parallel runs can be easily compared;
 * each row of matrix is precedded by the 5-tuple (body1, body2, sgp1, sgp2, row), where row = 1, 2 or 3;
 * the serial file can be line-sorted by the 5-tuples; the parallel files can be merged and similarly sorted */
static void output_system_matrix (LINSYS *sys, const char *path)
{
  FILE *file;
  DOM *dom;
  CON *con;

  ASSERT (file = fopen (path, "w"), ERR_FILE_OPEN);

  dom = sys->ldy->dom;

#if MPI
  if (sys->options & LOCAL_SYSTEM)
  {
#endif
  int i, j, *n, *r, *c, *e, *arows, *brows, *acolumns, *bcolumns, dimension, nnz;
  double *aelements, *belements, *a;

  ERRMEM (bcolumns = malloc (sys->cols [sys->dim] * sizeof (int)));
  ERRMEM (brows = calloc (sys->dim + 1, sizeof (int)));
  ERRMEM (belements = malloc (sys->cols [sys->dim] * sizeof (double)));
  ERRMEM (n = calloc (sys->dim,  sizeof (int)));

  /* (brows, bcolumns, belements) stores compressed row version of system matrix */

  nnz = sys->nnz;
  dimension = sys->dim;
  arows = sys->rows;
  acolumns = sys->cols;
  aelements = sys->a;

  for (i = 0; i < nnz; i ++) brows [arows [i] + 1] ++;

  for (i = 1; i <= dimension; i ++)  brows [i] += brows [i-1];

  for (i = 0; i < dimension; i ++)
  {
    for (r = &arows[acolumns [i]], e = &arows [acolumns[i+1]], a = &aelements [acolumns [i]]; r < e; r ++, a ++)
    {
      bcolumns [brows [*r] + n [*r]] = i;
      belements [brows [*r] + n [*r]] = *a;
      n [*r] ++;
    }
  }

  for (con = dom->con, i = 0; con; con = con->next)
  {
    for (j = i + 3; i < j; i ++)
    {
      c = &bcolumns[brows [i]];
      e = &bcolumns [brows[i+1]];
      a = &belements [brows [i]];

      qsort (a, e - c, sizeof (double), dblcmp);

      int body1,
	  body2,
	  sgp1,
	  sgp2;

      if (con->slave)
      {
	if (con->master->id < con->slave->id)
	{
	  body1 = con->master->id;
	  body2 = con->slave->id;
	  sgp1 = con->msgp - con->master->sgp;
	  sgp2 = con->ssgp - con->slave->sgp;
	}
	else
	{
	  body1 = con->slave->id;
	  body2 = con->master->id;
	  sgp1 = con->ssgp - con->slave->sgp;
	  sgp2 = con->msgp - con->master->sgp;
	}
      }
      else
      {
	body1 = con->master->id;
	body2 = -1;
	sgp1 = con->msgp - con->master->sgp;
	sgp2 = -1;
      }

      fprintf (file, "%d, %d, %d, %d, %d, %d, ", body1, body2, sgp1, sgp2, j-i, e - c);

      for (; c < e; c ++, a ++) fprintf (file, "%.1e, ", *a);

      fprintf (file, "\n");
    }
  }

  free (bcolumns);
  free (brows);
  free (belements);
  free (n);

#if MPI
  }
  else
  {
    int i, j, *c, *e, *rows, *columns;
    double *elements, *a, *onerow;

    ERRMEM (onerow = malloc (sys->dim * sizeof (double)));

    rows = sys->rows;
    columns = sys->cols;
    elements = sys->a;

    for (con = dom->con, i = 0; con; con = con->next)
    {
      for (j = i + 3; i < j; i ++)
      {
	c = &columns[rows [i]];
	e = &columns [rows[i+1]];
	a = &elements [rows [i]];

        blas_dcopy (e - c, a, 1, onerow, 1);
	qsort (onerow, e - c, sizeof (double), dblcmp);

	int body1,
	    body2,
	    sgp1,
	    sgp2;

	if (con->slave)
	{
	  if (con->master->id < con->slave->id)
	  {
	    body1 = con->master->id;
	    body2 = con->slave->id;
	    sgp1 = con->msgp - con->master->sgp;
	    sgp2 = con->ssgp - con->slave->sgp;
	  }
	  else
	  {
	    body1 = con->slave->id;
	    body2 = con->master->id;
	    sgp1 = con->ssgp - con->slave->sgp;
	    sgp2 = con->msgp - con->master->sgp;
	  }
	}
	else
	{
	  body1 = con->master->id;
	  body2 = -1;
	  sgp1 = con->msgp - con->master->sgp;
	  sgp2 = -1;
	}

	fprintf (file, "%d, %d, %d, %d, %d, %d, ", body1, body2, sgp1, sgp2, j-i, e - c);

	for (a = onerow; c < e; c ++, a ++) fprintf (file, "%.1e, ", *a);

	fprintf (file, "\n");
      }
    }
  }
#endif

  fclose (file);
}
#endif

/* test solution of W * x = b, where b = W * [1, 1, ..., 1];
 * return |x - [1, 1, ..., 1]| / |[1, 1, ..., 1]| */
double LINSYS_Test (LINSYS *sys, double accuracy, int maxiter)
{
  double *a, *b, *x, *e, error, one [3] = {1.0, 1.0, 1.0};
  LOCDYN *ldy;
  DIAB *dia;
  OFFB *blk;
  CON *con;

  ldy = sys->ldy;

  for (a = sys->a, b = a + sys->nnz; a < b; a ++) *a = 0.0; /* initialize with zeros */
  a = sys->a;
  b = sys->b;

  for (con = ldy->dom->con; con; con = con->next, b += 3)
  {
    dia = con->dia;

    double *W = dia->W;
    int *map = dia->map;

#if MPI
    a [map[0]]   = W[0];
    a [map[1]]   = W[1];
    a [map[2]]   = W[2];
    a [map[0]+1] = W[3];
    a [map[1]+1] = W[4];
    a [map[2]+1] = W[5];
    a [map[0]+2] = W[6];
    a [map[1]+2] = W[7];
    a [map[2]+2] = W[8];
#else
    a [map[0]]   = W[0];
    a [map[0]+1] = W[1];
    a [map[0]+2] = W[2];
    a [map[1]]   = W[3];
    a [map[1]+1] = W[4];
    a [map[1]+2] = W[5];
    a [map[2]]   = W[6];
    a [map[2]+1] = W[7];
    a [map[2]+2] = W[8];
#endif

    NVMUL (W, one, b);

#if MPI
    {
      OFFB *adj [2] = {dia->adj, dia->adjext};
      for (int i = 0; i < 2; i ++)
      for (blk = adj [i]; blk; blk = blk->n)
      {
	double *W = blk->W;
	int *map = blk->map;

	a [map[0]]   += W[0];
	a [map[1]]   += W[1];
	a [map[2]]   += W[2];
	a [map[0]+1] += W[3];
	a [map[1]+1] += W[4];
	a [map[2]+1] += W[5];
	a [map[0]+2] += W[6];
	a [map[1]+2] += W[7];
	a [map[2]+2] += W[8];

	NVADDMUL (b, W, one, b);
      }
    }
#else
    for (blk = dia->adj; blk; blk = blk->n)
    {
      double *W = blk->W;
      int *map = blk->map;

      a [map[0]]   += W[0];
      a [map[0]+1] += W[1];
      a [map[0]+2] += W[2];
      a [map[1]]   += W[3];
      a [map[1]+1] += W[4];
      a [map[1]+2] += W[5];
      a [map[2]]   += W[6];
      a [map[2]+1] += W[7];
      a [map[2]+2] += W[8];

      NVADDMUL (b, W, one, b);
    }
#endif
  }

#if DEBUG
  /* test for 0 on the diagonal */
  {
    int j, k, q, *i, *p;
    double *a;
    char *c;

#if MPI
    if (!(sys->options & LOCAL_SYSTEM))
    {
      i = sys->cols;
      p = sys->rows;
      c = "row";
      q = sys->hyp->ilower;
    }
    else
#endif
    {
      i = sys->rows;
      p = sys->cols;
      c = "column";
      q = 0;
    }

    for (j = 0, a = sys->a; j < sys->dim; j ++)
    {
      for (k = p [j]; k < p [j+1]; k ++, a ++, i ++)
      {
	if ((q+j) == *i)
	{
	  if (*a == 0.0)
	  {
#if MPI
	    printf ("LINSYS ERROR: rank %d\n", ldy->dom->rank);
#endif
	    printf ("LINSYS ERROR: diagonal zero for index %d\n", j);
	    printf ("LINSYS ERROR: system dimension is %d and number of nonzeros is %d\n", sys->dim, sys->nnz);
	    printf ("LINSYS ERROR: %s %d begins at %d and ends at %d\n", c, j, p [j], p [j+1]);
	  }
	}
      }
    }
  }
#endif

  LINSYS_Solve (sys, accuracy, maxiter);

  for (error = 0.0, x = sys->x, e = x + sys->dim; x < e; x ++)
  {
    error += (1.0 - (*x))*(1.0 - (*x));
  }

#if MPI
  double val_i [2] = {error, sys->dim}, val_o [2];
  MPI_Allreduce (val_i, val_o, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  error = sqrt (val_o [0] / MAX (val_o [1], 1.0));
#else
  error = sqrt (error / (double) MAX (sys->dim, 1));
#endif

#if 0
#if MPI
  char path [128];
  sprintf (path, "W_parallel_%d", ldy->dom->rank);
  output_system_matrix (sys, path);
#else
  output_system_matrix (sys, "W_serial");
#endif
  exit (0);
#endif

  return error;
}

/* returns 1 if a global system; 0 otherwise */
int LINSYS_Global (LINSYS *sys)
{
  return !(sys->options & LOCAL_SYSTEM);
}

/* most recent iterations count */
int LINSYS_Iters (LINSYS *sys)
{
  return sys->iters;
}

/* most recent residual norm */
double LINSYS_Resnorm (LINSYS *sys)
{
  return sys->resnorm;
}

/* destroy linear system */
void LINSYS_Destroy (LINSYS *sys)
{
  DIAB *dia;
  OFFB *blk;
  CON *con;

  if (sys->variant == SMOOTHED_VARIATIONAL ||
      sys->variant == NONSMOOTH_VARIATIONAL)
  {
    variational_coord_change_backward (sys); /* R = R - m(R) */
  }

  for (dia = sys->ldy->dia; dia; dia = dia->n)
  {
    con = dia->con;

    con->lin = NULL;
    con->Y = NULL;
    dia->map = NULL;

    for (blk = dia->adj; blk; blk = blk->n)
    {
      blk->map = NULL;
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      blk->map = NULL;
    }
#endif
  }
#if MPI
  if (!(sys->options & LOCAL_SYSTEM))
  {
    for (MAP *item = MAP_First (sys->ldy->dom->conext); item; item = MAP_Next (item))
    {
      CON *con = item->data;
      con->Y = NULL;
    }
  }
#endif

  MEM_Release (&sys->mapmem);
  MEM_Release (&sys->datmem);
  MEM_Release (&sys->ymem);

  free (sys->a);
  free (sys->x);
  free (sys->b);
  free (sys->rows);
  free (sys->cols);

#if SPQR
  cholmod_l_finish (&sys->spqr_c);
  free (sys->spqr_a);
  free (sys->spqr_b);
#elif UMFPACK
  umfpack_di_free_symbolic (&sys->symbolic);
#else
  LSS_Destroy (sys->lss);
#endif

#if  MPI
  if (sys->hyp)
  {
    HYPRE_DATA *hyp = sys->hyp;

    HYPRE_IJMatrixDestroy (hyp->A);
    HYPRE_IJVectorDestroy (hyp->B);
    HYPRE_IJVectorDestroy (hyp->X);
    MPI_Comm_free (&hyp->comm);
    free (hyp->ncols);
    free (hyp->rows);
    free (hyp);
  }
#endif

  free (sys);
}

/* constraint satisfaction merit function;
 * (assumes that both dia->R and dia->U are valid) */
double MERIT_Function (LOCDYN *ldy)
{
  LINSYS sys;

  sys.variant = NONSMOOTH_HSW;
  sys.ldy = ldy;

  return LINSYS_Merit (&sys, 0.0);
}
