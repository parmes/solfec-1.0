/*
 * glu.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * gluing solver
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
#include "glu.h"
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

struct glue
{
  LOCDYN *ldy; /* local dynamics */

  SET *subset; /* subset of gluing constraints */

  MEM mapmem, /* memory of DIAB->map and OFFB->map int [3] vectors */
      setmem; /* subset items memory */

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

  double resnorm;
  int iters;
};

static double glue_stiffness (CON *con, double step)
{
  BULK_MATERIAL *a = con->master->mat,
		*b = con->slave->mat;

  return step * 2.0 / (1.0/a->young + 1.0/b->young);
}

/* write system matrix */
static void write_system_matrix (DOM *dom, SET *subset, double *a)
{
  double *W, K, step;
  SET *item;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  int *map;

  step = dom->step;

  for (item = SET_First (subset); item; item = SET_Next (item))
  {
    con = item->data;
    dia = con->dia;
    W = dia->W;
    map = dia->map;
    K = glue_stiffness (con, step);

#if MPI
    /* NOTE:
     * in parallel global case 'map' references
     * beginnings of compressed row blocks */

    a [map[0]]   = W[0]*K + 1.0;
    a [map[1]]   = W[1]*K;
    a [map[2]]   = W[2]*K;
    a [map[0]+1] = W[3]*K;
    a [map[1]+1] = W[4]*K + 1.0;
    a [map[2]+1] = W[5]*K;
    a [map[0]+2] = W[6]*K;
    a [map[1]+2] = W[7]*K;
    a [map[2]+2] = W[8]*K + 1.0;

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;

      if (con->kind == GLUEPNT)
      {
	W = blk->W;
	map = blk->map;
	K = glue_stiffness (con, step);

	a [map[0]]   += W[0]*K; /* sum up off-diagonal blocks as pairs of W(i,j) are stored for two-body constraints */
	a [map[1]]   += W[1]*K;
	a [map[2]]   += W[2]*K;
	a [map[0]+1] += W[3]*K;
	a [map[1]+1] += W[4]*K;
	a [map[2]+1] += W[5]*K;
	a [map[0]+2] += W[6]*K;
	a [map[1]+2] += W[7]*K;
	a [map[2]+2] += W[8]*K;
      }
    }

    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = ((CON*)blk->dia);

      if (con->kind == GLUEPNT)
      {
	W = blk->W;
	map = blk->map;
	K = glue_stiffness (con, step);

	a [map[0]]   += W[0]*K; /* sum up off-diagonal blocks as pairs of W(i,j) are stored for two-body constraints */
	a [map[1]]   += W[1]*K;
	a [map[2]]   += W[2]*K;
	a [map[0]+1] += W[3]*K;
	a [map[1]+1] += W[4]*K;
	a [map[2]+1] += W[5]*K;
	a [map[0]+2] += W[6]*K;
	a [map[1]+2] += W[7]*K;
	a [map[2]+2] += W[8]*K;
      }
    }

#else
    /* NOTE:
     * in sequential or local case 'map' references
     * beginnings of compressed column blocks */

    a [map[0]]   = W[0]*K + 1.0;
    a [map[0]+1] = W[1]*K;
    a [map[0]+2] = W[2]*K;
    a [map[1]]   = W[3]*K;
    a [map[1]+1] = W[4]*K + 1.0;
    a [map[1]+2] = W[5]*K;
    a [map[2]]   = W[6]*K;
    a [map[2]+1] = W[7]*K;
    a [map[2]+2] = W[8]*K + 1.0;

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;

      if (con->kind == GLUEPNT)
      {
	W = blk->W;
	map = blk->map;
	K = glue_stiffness (con, step);

	a [map[0]]   += W[0]*K; /* sum up off-diagonal blocks as pairs of W(i,j) are stored for two-body constraints */
	a [map[0]+1] += W[1]*K;
	a [map[0]+2] += W[2]*K;
	a [map[1]]   += W[3]*K;
	a [map[1]+1] += W[4]*K;
	a [map[1]+2] += W[5]*K;
	a [map[2]]   += W[6]*K;
	a [map[2]+1] += W[7]*K;
	a [map[2]+2] += W[8]*K;
      }
    }
#endif
  }
}

/* update right hand side */
static void update_right_hand_side (DOM *dom, SET *subset, double *b)
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
    y = &b [con->num * 3];
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
      con = ((CON*)blk->dia);
      if (con->kind != GLUEPNT)
      {
	W = blk->W;
	R = con->R;
	NVADDMUL (y, W, R, y);
      }
    }
#endif
  }
}

/* create gluing solver */
GLUE* GLUE_Create (LOCDYN *ldy)
{
  int i, j, k, l, n, ncon;
  GLUE *glu;
  SET *jtem;
  MAP *item;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  DOM *dom;
  MEM mem;
#if MPI
  int *diacnt, /* this-processor column counts per row */
      *offcnt; /* off-processor column counts per row */
#endif

  dom = ldy->dom;
    
  ERRMEM (glu = MEM_CALLOC (sizeof (GLUE)));

  glu->ldy = ldy;

  MEM_Init (&glu->mapmem, sizeof (int [3]), MEMBLK);
  MEM_Init (&glu->setmem, sizeof (SET), MEMBLK);

  /* collect gluing constraints */
  for (con = dom->con; con; con = con->next)
  {
    if (con->kind == GLUEPNT)
    {
      SET_Insert (&glu->setmem, &glu->subset, con, NULL);
    }
  }
  ncon = SET_Size (glu->subset);

  /* gluing constraints numbering */
  DOM_Number_Constraints (dom, 0, glu->subset);

  /* matrix structure */
#if MPI
  int basenum = 0;
  MAP **prow; /* create row mapping in DIAB->map and OFFB->map */

  /* first constraint number */
  if (glu->subset)
  {
    jtem = SET_First (glu->subset);
    con = jtem->data;
    basenum = con->num;
  }

  /* allocate colum blocks mapping */
  ERRMEM (diacnt = MEM_CALLOC (ncon * sizeof (int [3])));
  ERRMEM (offcnt = MEM_CALLOC (ncon * sizeof (int [3])));
  ERRMEM (prow = MEM_CALLOC (ncon * sizeof (MAP*)));
  MEM_Init (&mem, sizeof (MAP), MEMBLK);

  /* map colum blocks */
  for (jtem = SET_First (glu->subset); jtem; jtem = SET_Next (jtem))
  {
    int col, row, *map;

    con = jtem->data;
    dia = con->dia;

    col = con->num;
    row = col - basenum; /* use relative numbering for rows */
    ASSERT_DEBUG (!MAP_Find (prow [row], (void*) (long) (col+1), NULL), "Diagonal block mapped twice");
    ERRMEM (map = MEM_Alloc (&glu->mapmem));
    ERRMEM (MAP_Insert (&mem, &prow [row], (void*) (long) (col+1), map, NULL)); /* '+1' not to confuse with NULL */

    for (i = k = row * 3, j = k + 3; i < j; i ++) diacnt [i] = 3; /* diagonal blocks contribute 3 this-processor columns */

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;

      if (con->kind == GLUEPNT) /* include only gluing adjacency */
      {
	col = con->num;

	if (!MAP_Find (prow [row], (void*) (long) (col+1), NULL)) /* there are two W(i,j) blocks for two-body constraints */
	{
	  ERRMEM (map = MEM_Alloc (&glu->mapmem));
	  ERRMEM (MAP_Insert (&mem, &prow [row], (void*) (long) (col+1), map, NULL));

	  for (i = k; i < j; i ++) diacnt [i] += 3; /* off-diagonal local blocks contribute 3 this-processor columns */
	}
      }
    }
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = ((CON*)blk->dia);

      if (con->kind == GLUEPNT) /* include only gluing adjacency */
      {
	col = con->num; /* external constraint */

	if (!MAP_Find (prow [row], (void*) (long) (col+1), NULL)) /* there are two W(i,j) blocks for two-body constraints */
	{
	  ERRMEM (map = MEM_Alloc (&glu->mapmem));
	  ERRMEM (MAP_Insert (&mem, &prow [row], (void*) (long) (col+1), map, NULL));

	  for (i = k; i < j; i ++) offcnt [i] += 3; /* off-diagonal external blocks contribute 3 off-processor columns */
	}
      }
    }
  }

  for (j = 0, glu->nnz = 0; j < ncon; j ++)
  {
    glu->nnz += 9 * MAP_Size (prow [j]); /* number of nonzeros */
  }
  glu->dim = ncon * 3; /* system dimension */

  /* eallocate compressed column storage */
  ERRMEM (glu->a = malloc (sizeof (double) * glu->nnz));
  ERRMEM (glu->cols = malloc (sizeof (int) * glu->nnz));
  ERRMEM (glu->rows = MEM_CALLOC (sizeof (int) * (glu->dim+1))); /* '+ 1' as there is cols[0] == 0 always */

  int *rows = glu->rows,
      *cols = glu->cols,
      *aux;

  ERRMEM (aux = malloc (sizeof (int) * glu->dim));

  /* set up rows sizes */
  for (j = 0; j < ncon; j ++)
  {
    n = 3 * MAP_Size (prow [j]);

    rows [3*j+1] = n; /* '+1' so that cols[0] == 0 */
    rows [3*j+2] = n;
    rows [3*j+3] = n;
  }

  /* compute rows pointers */
  for (n = 1; n <= glu->dim; n ++)
  {
    rows [n] += rows [n-1];
    aux [n-1] = rows [n-1]; /* initialize aux with cols  */
  }

  ASSERT_DEBUG (rows [glu->dim] == glu->nnz, "Inconsistent sparse storage");

  /* set up column pointers */
  for (j = 0; j < ncon; j ++)
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
  for (jtem = SET_First (glu->subset); jtem; jtem = SET_Next (jtem))
  {
    int col, row;

    con = jtem->data;
    dia = con->dia;

    col = dia->con->num;
    row = col - basenum; /* use relative numbering for rows */

    ASSERT_DEBUG_EXT (dia->map = MAP_Find (prow [row], (void*) (long) (col+1), NULL), "Inconsistent sparse storage");

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
      if (con->kind == GLUEPNT)
      {
	col = con->num;
	ASSERT_DEBUG_EXT (blk->map = MAP_Find (prow [row], (void*) (long) (col+1), NULL), "Inconsistent sparse storage");
      }
    }
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = ((CON*)blk->dia);
      if (con)
      {
	col = con->num;
	ASSERT_DEBUG_EXT (blk->map = MAP_Find (prow [row], (void*) (long) (col+1), NULL), "Inconsistent sparse storage");
      }
    }
  }

  /* free auxiliary space */
  MEM_Release (&mem);
  free (prow);
  free (aux);
#else
  MAP **pcol;

  /* allocate colum blocks mapping */
  ERRMEM (pcol = MEM_CALLOC (ncon * sizeof (MAP*)));
  MEM_Init (&mem, sizeof (MAP), MEMBLK);

  /* map colum blocks */
  for (jtem = SET_First (glu->subset); jtem; jtem = SET_Next (jtem))
  {
    int col, row, *map;

    con = jtem->data;
    dia = con->dia;

    col = row = con->num;
    ASSERT_DEBUG (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Diagonal block mapped twice");
    ERRMEM (map = MEM_Alloc (&glu->mapmem));
    ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL)); /* '+1' not to confuse with NULL */

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
      if (con->kind == GLUEPNT) /* include only gluing adjacency */
      {
	col = con->num;

	if (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL)) /* there are two W(i,j) blocks for two-body constraints */
	{
	  ERRMEM (map = MEM_Alloc (&glu->mapmem));
	  ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL));
	}
      }
    }
  }

  for (j = 0, glu->nnz = 0; j < ncon; j ++)
  {
    glu->nnz += 9 * MAP_Size (pcol [j]); /* number of nonzeros */
  }
  glu->dim = ncon * 3; /* system dimension */

  /* eallocate compressed column storage */
  ERRMEM (glu->a = MEM_CALLOC (sizeof (double) * glu->nnz));
  ERRMEM (glu->rows = malloc (sizeof (int) * glu->nnz));
  ERRMEM (glu->cols = MEM_CALLOC (sizeof (int) * (glu->dim+1))); /* '+ 1' as there is cols[0] == 0 always */

  int *rows = glu->rows,
      *cols = glu->cols,
      *aux;

  ERRMEM (aux = malloc (sizeof (int) * glu->dim));

  /* set up column sizes */
  for (j = 0; j < ncon; j ++)
  {
    n = 3 * MAP_Size (pcol [j]);

    cols [3*j+1] = n; /* '+1' so that cols[0] == 0 */
    cols [3*j+2] = n;
    cols [3*j+3] = n;
  }

  /* compute column pointers */
  for (n = 1; n <= glu->dim; n ++)
  {
    cols [n] += cols [n-1];
    aux [n-1] = cols [n-1]; /* initialize aux with cols  */
  }

  ASSERT_DEBUG (cols [glu->dim] == glu->nnz, "Inconsistent sparse storage");

  /* set up row pointers */
  for (j = 0; j < ncon; j ++)
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
  for (jtem = SET_First (glu->subset); jtem; jtem = SET_Next (jtem))
  {
    int col, row;

    con = jtem->data;
    dia = con->dia;

    col = row = con->num;

    ASSERT_DEBUG_EXT (dia->map = MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Inconsistent sparse storage");

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
      if (con->kind == GLUEPNT)
      {
	col = con->num;

	ASSERT_DEBUG_EXT (blk->map = MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Inconsistent sparse storage");
      }
    }
  }

  /* free auxiliary space */
  MEM_Release (&mem);
  free (pcol);
  free (aux);
#endif

  /* right hand side */
  ERRMEM (glu->b = MEM_CALLOC (sizeof (double) * glu->dim));

  /* solution vector */
  ERRMEM (glu->x = MEM_CALLOC (sizeof (double) * glu->dim));

  /* linear solver */
#if MPI
  MPI_Group allcpus, nonzero;
  int *nconall, *ranks;
  HYPRE_DATA *hyp;

  ERRMEM (ranks = MEM_CALLOC (sizeof (int [dom->ncpu])));
  ERRMEM (nconall = MEM_CALLOC (sizeof (int [dom->ncpu])));
  ERRMEM (hyp = MEM_CALLOC (sizeof (HYPRE_DATA)));
  glu->hyp = hyp;

  /* gather numbers of constraints from all subdomains */
  MPI_Allgather (&ncon, 1, MPI_INT, nconall, 1, MPI_INT, MPI_COMM_WORLD);
  for (i = j = 0; i < dom->ncpu; i ++)
  {
    if (nconall [i]) ranks [j ++] = i;
  }

  /* create a group of processorss with nonzero constraint counts */
  MPI_Comm_group (MPI_COMM_WORLD, &allcpus);
  MPI_Group_incl (allcpus, j, ranks, &nonzero);

  /* create HYPRE communicator */
  MPI_Comm_create (MPI_COMM_WORLD, nonzero, &hyp->comm);

  /* clean */
  free (ranks);
  free (nconall);

  if (glu->subset) /* there are some constraints here */
  {
    hyp->ilower = basenum * 3;
    hyp->iupper = hyp->ilower + ncon * 3 - 1;
    hyp->jlower = hyp->ilower;
    hyp->jupper = hyp->iupper;
    ERRMEM (hyp->ncols = malloc (sizeof (int [glu->dim])));
    ERRMEM (hyp->rows = malloc (sizeof (int [glu->dim])));
    for (j = 0, k = hyp->ilower; j < glu->dim; j ++, k ++)
    {
      hyp->ncols [j] = glu->rows [j+1] - glu->rows [j];
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
    HYPRE_IJMatrixSetValues (hyp->A, glu->dim, hyp->ncols, hyp->rows, glu->cols, glu->a);
    HYPRE_IJMatrixAssemble (hyp->A);
    HYPRE_IJMatrixGetObject (hyp->A, (void **) &hyp->AP);

    /* b */
    HYPRE_IJVectorCreate (hyp->comm, hyp->jlower, hyp->jupper, &hyp->B);
    HYPRE_IJVectorSetObjectType (hyp->B, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize (hyp->B);
    HYPRE_IJVectorSetValues (hyp->B, glu->dim, hyp->rows, glu->b);
    HYPRE_IJVectorAssemble (hyp->B);
    HYPRE_IJVectorGetObject(hyp->B, (void **) &hyp->BP);

    /* x */
    HYPRE_IJVectorCreate (hyp->comm, hyp->jlower, hyp->jupper, &hyp->X);
    HYPRE_IJVectorSetObjectType (hyp->X, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize (hyp->X);
    HYPRE_IJVectorSetValues (hyp->X, glu->dim, hyp->rows, glu->x);
    HYPRE_IJVectorAssemble (hyp->X);
    HYPRE_IJVectorGetObject(hyp->X, (void **) &hyp->XP);
  }
#else
#if SPQR
  /* start CHOLMOD */
  cholmod_l_start (&glu->spqr_c);

  /* allocate A */
  ERRMEM (glu->spqr_a = MEM_CALLOC (sizeof (cholmod_sparse)));
  glu->spqr_a->nrow = glu->spqr_a->ncol = glu->dim;
  glu->spqr_a->nzmax = glu->nnz;
  glu->spqr_a->p = glu->cols;
  glu->spqr_a->i = glu->rows;
  glu->spqr_a->x = glu->a;
  glu->spqr_a->stype = 0;
  glu->spqr_a->itype = CHOLMOD_INT;
  glu->spqr_a->xtype = CHOLMOD_REAL;
  glu->spqr_a->dtype = CHOLMOD_DOUBLE;
  glu->spqr_a->sorted = TRUE;
  glu->spqr_a->packed = TRUE;

  /* allocate b */
  ERRMEM (glu->spqr_b = MEM_CALLOC (sizeof (cholmod_sparse)));
  glu->spqr_b->nrow = glu->spqr_b->nzmax = glu->spqr_b->d = glu->dim;
  glu->spqr_b->ncol = 1;
  glu->spqr_b->x = glu->b;
  glu->spqr_b->xtype = CHOLMOD_REAL;
  glu->spqr_b->dtype = CHOLMOD_DOUBLE;
#elif UMFPACK
  ASSERT (umfpack_di_symbolic (glu->dim, glu->dim, glu->cols, glu->rows,
	  glu->a, &glu->symbolic, NULL, NULL) == UMFPACK_OK, ERR_UMFPACK_SOLVE);
#else
  ERRMEM (glu->lss = LSS_Create (glu->dim, glu->a, glu->cols, glu->rows));
  LSS_Set (glu->lss, LSS_RESTART, 100); /* FIXME: make user adjustable or automatic */
#endif
#endif

  /* prepare system matrix */
  write_system_matrix (dom, glu->subset, glu->a);

  return glu; 
}

/* compute gluing reactions */
void GLUE_Solve (GLUE *glu, double accuracy, int maxiter)
{
  DOM *dom = glu->ldy->dom;

  update_right_hand_side (dom, glu->subset, glu->b);

#if MPI
  {
    HYPRE_DATA *hyp = glu->hyp;

    if (hyp->comm != MPI_COMM_NULL) /* skip empty constraints case */
    {
      /* A */
      HYPRE_IJMatrixInitialize (hyp->A);
      HYPRE_IJMatrixSetValues (hyp->A, glu->dim, hyp->ncols, hyp->rows, glu->cols, glu->a);
      HYPRE_IJMatrixAssemble (hyp->A);

      /* b */
      HYPRE_IJVectorInitialize (hyp->B);
      HYPRE_IJVectorSetValues (hyp->B, glu->dim, hyp->rows, glu->b);
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
	HYPRE_FlexGMRESGetNumIterations (solver, &glu->iters);
	HYPRE_FlexGMRESGetFinalRelativeResidualNorm (solver, &glu->resnorm);

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
      HYPRE_IJVectorGetValues (hyp->X, glu->dim, hyp->rows, glu->x);
    }
  }
#else
#if SPQR
  {
    /* x = A\b */
    cholmod_dense *x = SuiteSparseQR_C_backslash_default (glu->spqr_a, glu->spqr_b, &glu->spqr_c);
    blas_dcopy (glu->dim, x->x, 1, glu->x, 1);
    cholmod_l_free_dense (&x, &glu->spqr_c) ;

    /* stats */
    glu->resnorm = 0.0;
    glu->iters = 0;
  }
#elif UMFPACK
  {
    /* x = A\b */
    void *numeric;
    ASSERT (umfpack_di_numeric (glu->cols, glu->rows, glu->a, glu->symbolic, &numeric, NULL, NULL) == UMFPACK_OK, ERR_UMFPACK_SOLVE);
    ASSERT (umfpack_di_solve (UMFPACK_A, glu->cols, glu->rows, glu->a, glu->x, glu->b, numeric, NULL, NULL) == UMFPACK_OK, ERR_UMFPACK_SOLVE);
    umfpack_di_free_numeric (&numeric);

    /* stats */
    glu->resnorm = 0.0;
    glu->iters = 0;
  }
#else
  {
    /* x = A\b */
    LSS_Set (glu->lss, LSS_ABSOLUTE_ACCURACY, accuracy);
    LSS_Set (glu->lss, LSS_ITERATIONS_BOUND, maxiter);
    LSS_Set (glu->lss, LSS_RELATIVE_ACCURACY, 1.0);
    LSS_Set (glu->lss, LSS_NORMALIZE_ROWS, 0);
#if 0
    LSS_Set (glu->lss, LSS_PRECONDITIONER, 3);
    LSS_Set (glu->lss, LSS_DECIMATION, 16);
    LSS_Set (glu->lss, LSS_SMOOTHING_STEPS, 5);
#endif
    if (LSS_Solve (glu->lss, glu->a, glu->x, glu->b) != LSSERR_NONE)
    {
      fprintf (stderr, "WARNING: LSS failed with message: %s\n", LSS_Errmsg (glu->lss));
    }

    /* stats */
    //FIXME: invalid write
    //glu->resnorm = LSS_Get (glu->lss, LSS_ABSOLUTE_ERROR);
    //glu->iters = LSS_Get (glu->lss, LSS_ITERATIONS);
  }
#endif
#endif

  double *x;
  SET *item;
  CON *con;

  /* write R, U */
  for (item = SET_First (glu->subset), x = glu->x; item; item = SET_Next (item), x += 3)
  {
    con = item->data;

    double *U = con->U,
	   *R = con->R,
	    K = -glue_stiffness (con, dom->step);

    COPY (x, U);
    MUL (U, K, R);
  }

#if MPI
  /* update external gluing R */
  DOM_Update_External_Reactions (dom, 0); /* FIXME: should be on glu->subset only */
#endif
}

/* destroy gluing solver */
void GLUE_Destroy (GLUE *glu)
{
  MEM_Release (&glu->mapmem);
  MEM_Release (&glu->setmem);

  free (glu->a);
  free (glu->x);
  free (glu->b);
  free (glu->rows);
  free (glu->cols);

#if SPQR
  cholmod_l_finish (&glu->spqr_c);
  free (glu->spqr_a);
  free (glu->spqr_b);
#elif UMFPACK
  umfpack_di_free_symbolic (&glu->symbolic);
#else
  LSS_Destroy (glu->lss);
#endif

#if  MPI
  if (glu->hyp)
  {
    HYPRE_DATA *hyp = glu->hyp;

    HYPRE_IJMatrixDestroy (hyp->A);
    HYPRE_IJVectorDestroy (hyp->B);
    HYPRE_IJVectorDestroy (hyp->X);
    MPI_Comm_free (&hyp->comm);
    free (hyp->ncols);
    free (hyp->rows);
    free (hyp);
  }
#endif

  free (glu);
}
