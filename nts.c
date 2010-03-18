/*
 * nts.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Newton constraints solver
 */

#if MPI
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_IJ_mv.h>
#elif UMFPACK
#include <umfpack.h>
#else
#include "lss.h"
#endif

#include <stdlib.h>
#include <time.h>
#include "dom.h"
#include "bgs.h"
#include "nts.h"
#include "alg.h"
#include "bla.h"
#include "lap.h"
#include "err.h"

#define MEMBLK 256

/* mask separating variants from variant options */
#define VARIANT_MASK (NT_NONSMOOTH_HSW|NT_NONSMOOTH_HYBRID|NT_FIXED_POINT|NT_NONSMOOTH_VARIATIONAL|NT_SMOOTHED_VARIATIONAL)

/* to be used in a switch */
#define VARIANT(variant) ((int)((variant) & VARIANT_MASK))

typedef struct linsys LINSYS; /* linear sysetm */

#if MPI
struct linsys /* HYPRE IJ matrix data and the rest of linear system */
{
  int ilower, /* lower row index == dom->con->num * 3 */
      iupper; /* upper low index == ilowr + dom->ncon * 3 - 1 */

  int jlower, /* lower column index == ilower */
      jupper; /* upper column index == iupper */

  int nrows, /* number of rows == dom->ncon * 3 */
     *ncols, /* number of columns in each row (depends on local dynamics) */
     *rows,  /* row indices: ilower, ilower + 1, ..., iupper */
     *cols;  /* column indices: depend on local dynamics */

  double *a; /* compressed row values */

  HYPRE_IJMatrix A;  /* HYPRE IJ matrix composed from the above data */
  HYPRE_ParCSRMatrix AP; /* reference to the parallel HYPRE CSR image of A */

  int nvalues, /* number of vector values */
     *indices; /* vector indices */

  double *b; /* right hand side */

  HYPRE_IJVector B; /* HYPRE IJ right hand side vector */
  HYPRE_ParVector BP; /* HYPRE parallel image of B */

  double *x; /* solution */

  HYPRE_IJVector X; /* HYPRE IJ solution vector */
  HYPRE_ParVector XP; /* HYPRE parallel image of X */

  int nnz, dim; /* FIXME */
};
#else
struct linsys /* SuperLU matrix data and the rest of linear system */
{
  double *a; /* compressed column values */

  int *cri, /* column row indices */
      *cbp, /* column beginning pointers */
       nnz, /* number of nonzero entries */
       dim; /* dimension of the rectangular matrix */

  double *b; /* right hand side */

  double *x; /* solution */

#if UMFPACK
  void *symbolic; /* UMFPACK symbolic factorization */
#else
  void *lss; /* LSS interface */
#endif
};
#endif

typedef struct conaux CONAUX; /* auxiliary constraint data */

struct conaux
{
  double RES [3], /* residual */
	 DR [3],  /* reaction increment */
	 DU [3],  /* velocity increment */
	 RN;      /* normal reaction bound during fixed point iterations */

  short sticks,   /* contact sticking history bits */
	hits;     /* nuber of alternating stick-slip patterns found */
};

#define CON_RES(con) ((CONAUX*)((con)->data))->RES
#define CON_DR(con) ((CONAUX*)((con)->data))->DR
#define CON_DU(con) ((CONAUX*)((con)->data))->DU
#define CON_RN(con) ((CONAUX*)((con)->data))->RN
#define CON_STICKS(con) ((CONAUX*)((con)->data))->sticks
#define CON_HITS(con) ((CONAUX*)((con)->data))->hits

/* create linear system */
#if MPI
static LINSYS* system_create (MEM *mapmem, MEM *auxmem, LOCDYN *ldy)
{
  /* TODO */

  return NULL;
}
#else
static LINSYS* system_create (MEM *mapmem, MEM *auxmem, LOCDYN *ldy)
{
  DOM *dom = ldy->dom;
  MAP **pcol, *item;
  int i, j, k, l, n;
  LINSYS *sys;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  MEM mem;
    
  ERRMEM (sys = MEM_CALLOC (sizeof (LINSYS)));

  /* number constraints */
  DOM_Number_Constraints (dom);

  /* allocate auxiliary space */
  for (con = dom->con; con; con = con->next)
  {
    ERRMEM (con->data = MEM_Alloc (auxmem));
  }

  /* allocate colum blocks mapping */
  MEM_Init (&mem, sizeof (MAP), MEMBLK);
  ERRMEM (pcol = MEM_CALLOC (dom->ncon * sizeof (MAP*)));

  /* map colum blocks */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    int col, row, *map;

    col = row = dia->con->num;
    ASSERT_DEBUG (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL), "Diagonal block mapped twice");
    ERRMEM (map = MEM_Alloc (mapmem));
    ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL)); /* '+1' not to confuse with NULL */

    for (blk = dia->adj; blk; blk = blk->n)
    {
      col = blk->dia->con->num;

      if (!MAP_Find (pcol [col], (void*) (long) (row+1), NULL)) /* there are two W(i,j) blocks for two-body constraints */
      {
        ERRMEM (map = MEM_Alloc (mapmem));
	ERRMEM (MAP_Insert (&mem, &pcol [col], (void*) (long) (row+1), map, NULL));
      }
    }
  }

#if DEBUG
  /* self consistency checks */
  for (j = 0; j < dom->ncon; j ++)
  {
    for (item = MAP_First (pcol [j]); item; item = MAP_Next (item))
    {
      int *map = item->data;

      ASSERT_DEBUG (map [0] == 0 && map [1] == 0 && map [2] == 0, "Nonzero triplet: %d, %d, %d", map [0], map [1], map [2]);
    }
  }

  for (j = 0; j < dom->ncon; j ++)
  {
    for (item = MAP_First (pcol [j]); item; item = MAP_Next (item))
    {
      int *map = item->data;
      int i = (((int) (long) item->key) - 1);

      ASSERT_DEBUG (map [0] == 0 && map [1] == 0 && map [2] == 0, "Duplicated map from: %d, %d at %d, %d", map [0], map [1], i, j);

      map [0] = i;
      map [1] = j;
    }
  }

  for (j = 0; j < dom->ncon; j ++)
  {
    for (item = MAP_First (pcol [j]); item; item = MAP_Next (item))
    {
      int *map = item->data;
      map [0] = map [1] = map [2] = 0;
    }
  }
#endif

  for (j = 0, sys->nnz = 0; j < dom->ncon; j ++)
  {
    sys->nnz += 9 * MAP_Size (pcol [j]); /* number of nonzeros */
  }
  sys->dim = dom->ncon * 3; /* system dimension */

  /* eallocate compressed column storage */
  ERRMEM (sys->a = malloc (sizeof (double) * sys->nnz));
  ERRMEM (sys->cri = malloc (sizeof (int) * sys->nnz));
  ERRMEM (sys->cbp = MEM_CALLOC (sizeof (int) * (sys->dim+1))); /* '+ 1' as there is cbp[0] == 0 always */

  int *cri = sys->cri,
      *cbp = sys->cbp,
      *aux;

  ERRMEM (aux = malloc (sizeof (int) * sys->dim));

  /* set up column sizes */
  for (j = 0; j < dom->ncon; j ++)
  {
    n = 3 * MAP_Size (pcol [j]);

    cbp [3*j+1] = n; /* '+1' so that cpb[0] == 0 */
    cbp [3*j+2] = n;
    cbp [3*j+3] = n;
  }

  /* compute column pointers */
  for (n = 1; n <= sys->dim; n ++)
  {
    cbp [n] += cbp [n-1];
    aux [n-1] = cbp [n-1]; /* initialize aux with cpb  */
  }

  ASSERT_DEBUG (cbp [sys->dim] == sys->nnz, "Inconsistent sparse storage");

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
          cri [l + i] = k + i; /* set up row index */
	}
	aux [3*j+n] += 3; /* increment relative row pointers */
      }
    }
  }

#if DEBUG
  /* self consistency check */
  for (j = 0; j < dom->ncon; j ++)
  {
    for (item = MAP_First (pcol [j]); item; item = MAP_Next (item))
    {
      int *map = item->data;

      for (n = 0; n < 3; n ++)
      {
	k = 3*j+n;

	ASSERT_DEBUG (map [n] >= cbp[k] && map [n] < cbp [k+1], "Inconsistent column blocks mapping");
      }
    }
  }
#endif

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

  /* clean */
  MEM_Release (&mem);
  free (pcol);
  free (aux);

  /* alloc right hand side and solution vectors */
  ERRMEM (sys->b = MEM_CALLOC (sizeof (double) * sys->dim));
  ERRMEM (sys->x = MEM_CALLOC (sizeof (double) * sys->dim));

  /* initialize linear solver */

#if UMFPACK
  ASSERT (umfpack_di_symbolic (sys->dim, sys->dim, sys->cbp, sys->cri, sys->a, &sys->symbolic, NULL, NULL) >= 0,  ERR_UMFPACK_SYMBOLIC);
#else
  ERRMEM (sys->lss = LSS_Create (sys->dim, sys->a, sys->cbp, sys->cri));

  LSS_Set (sys->lss, LSS_PRECONDITIONER, 3);
  LSS_Set (sys->lss, LSS_DECIMATION, 8);
  LSS_Set (sys->lss, LSS_RESTART, 100);
  LSS_Set (sys->lss, LSS_SMOOTHING_STEPS, 3);
#endif

  return sys; 
}
#endif

/* update linear system for NT_NONSMOOTH_HSW, NT_NONSMOOTH_HYBRID, NT_FIXED_POINT variants */
static void system_update_HSW_HYBRID_FIXED (LINSYS *sys, NTVARIANT variant, LOCDYN *ldy)
{
  double d [3], norm, lim, udash, step;
  short dynamic, pull;
  double *a, *b;
  DIAB *dia;
  OFFB *blk;
  DOM *dom;
  CON *con;

#if MPI
  /* update external reactions (needed below) */
  DOM_Update_External_Reactions (ldy->dom, 0);
#endif

  dom = ldy->dom;
  dynamic = dom->dynamic;
  step = dom->step;
  for (a = sys->a, b = a + sys->nnz; a < b; a ++) *a = 0.0; /* zero input matrix */
  a = sys->a;
  b = sys->b;

  for (con = dom->con; con; con = con->next, b += 3)
  {
    dia = con->dia;

    double *W = dia->W,
	   *B = dia->B,
	   *U = dia->U,
	   *R = dia->R,
	   *RES = CON_RES(con);

    /* update residual */
    NVADDMUL (B, W, R, RES);
    SUB (RES, U, RES);
    for (blk = dia->adj; blk; blk = blk->n)
    {
      double *R = blk->dia->R,
             *W = blk->W;

      NVADDMUL (RES, W, R, RES);
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      double *R = ((CON*)blk->dia)->R, /* external reaction */
             *W = blk->W;

      NVADDMUL (RES, W, R, RES);
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
	/* NOTE:
	 * in parallel case 'map' references
	 * beginnings of compressed row blocks */

	a [map[2]]   = W[2];
	a [map[2]+1] = W[5];
	a [map[2]+2] = W[8];
#else
	/* NOTE:
	 * in sequential case 'map' references
	 * beginnings of compressed column blocks */

	a [map[0]+2] = W[2];
	a [map[1]+2] = W[5];
	a [map[2]+2] = W[8];
#endif
	b [2] = -udash - RES[2];

#if MPI
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
#else
	for (blk = dia->adj; blk; blk = blk->n)
	{
	  double *W = blk->W;
	  int *map = blk->map;

	  a [map[0]+2] += W[2];
	  a [map[1]+2] += W[5];
	  a [map[2]+2] += W[8];
	}
#endif
      }
      else  /* inactive normal set */
      {
	pull = 1;

#if MPI
	a [map[2]]   = 0.0;
	a [map[2]+1] = 0.0;
	a [map[2]+2] = W[8];
#else
	a [map[0]+2] = 0.0;
	a [map[1]+2] = 0.0;
	a [map[2]+2] = W[8];
#endif
	b [2] = -R[2] * W[8];

	/* FIXME: mind zero initialization in the begining; hence below not needed */
#if 0
#if MPI
	{
	  OFFB *adj [2] = {dia->adj, dia->adjext};
	  for (int i = 0; i < 2; i ++)
	  for (blk = adj [i]; blk; blk = blk->n)
	  {
	    int *map = blk->map;

	    a [map[2]]   = 0.0;
	    a [map[2]+1] = 0.0;
	    a [map[2]+2] = 0.0;
	  }
	}
#else
	for (blk = dia->adj; blk; blk = blk->n)
	{
	  int *map = blk->map;

	  a [map[0]+2] = 0.0;
	  a [map[1]+2] = 0.0;
	  a [map[2]+2] = 0.0;
	}
#endif
#endif
      }

      /* tangential response */

      CON_STICKS(con) = CON_STICKS(con) << 1; /* rewind history */

      norm = sqrt (d[0]*d[0]+d[1]*d[1]); /* tangential force value */

      switch (VARIANT (variant))
      {
	case NT_NONSMOOTH_HSW:
	  lim = fri * MAX (0, d[2]);
	  break;
	case NT_NONSMOOTH_HYBRID:
	  lim = fri * MAX (0, R[2]);
	  break;
	case NT_FIXED_POINT:
	  lim = fri * MAX (0, CON_RN(con));
	  break;
      }

      if ((variant & NT_NONSMOOTH_HSW) && pull) goto ZERO_TANG; /* enforce AN = AT + IT */

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

	  if ((variant & NT_SYMMETRIZE) == 0) /* nonsymetric reduction */
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
	    a [map[0]]   = H[0] + rho*(M[0]*W[0] + M[2]*W[1]);
	    a [map[1]]   = H[1] + rho*(M[1]*W[0] + M[3]*W[1]);
	    a [map[0]+1] = H[2] + rho*(M[0]*W[3] + M[2]*W[4]);
	    a [map[1]+1] = H[3] + rho*(M[1]*W[3] + M[3]*W[4]);
#else
	    a [map[0]]   = H[0] + rho*(M[0]*W[0] + M[2]*W[1]);
	    a [map[0]+1] = H[1] + rho*(M[1]*W[0] + M[3]*W[1]);
	    a [map[1]]   = H[2] + rho*(M[0]*W[3] + M[2]*W[4]);
	    a [map[1]+1] = H[3] + rho*(M[1]*W[3] + M[3]*W[4]);
#endif

	    if (variant & NT_FIXED_POINT)
	    {
#if MPI
	      a [map[0]+2] = rho*(M[0]*W[6] + M[2]*W[7]);
	      a [map[1]+2] = rho*(M[1]*W[6] + M[3]*W[7]);
#else
	      a [map[2]]   = rho*(M[0]*W[6] + M[2]*W[7]);
	      a [map[2]+1] = rho*(M[1]*W[6] + M[3]*W[7]);
#endif
	      b [0] = (fri*(d[0]/norm)*CON_RN(con) - R[0]) * W[0] - rho*(M[0]*RES[0] + M[2]*RES[1]); 
	      b [1] = (fri*(d[1]/norm)*CON_RN(con) - R[1]) * W[4] - rho*(M[1]*RES[0] + M[3]*RES[1]); 
	    }
	    else /* HSW, HYBRID */
	    {
#if MPI
	      a [map[0]+2] = rho*(M[0]*W[6] + M[2]*W[7]) - fri*(d[0]/norm) * W[0];
	      a [map[1]+2] = rho*(M[1]*W[6] + M[3]*W[7]) - fri*(d[1]/norm) * W[4];
#else
	      a [map[2]]   = rho*(M[0]*W[6] + M[2]*W[7]) - fri*(d[0]/norm) * W[0];
	      a [map[2]+1] = rho*(M[1]*W[6] + M[3]*W[7]) - fri*(d[1]/norm) * W[4];
#endif
	      b [0] = (fri*(d[0]/norm)*R[2] - R[0]) * W[0] - rho*(M[0]*RES[0] + M[2]*RES[1]);
	      b [1] = (fri*(d[1]/norm)*R[2] - R[1]) * W[4] - rho*(M[1]*RES[0] + M[3]*RES[1]);
	    }

#if MPI
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
#else
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
#endif
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
	    a [map[0]]   = W[0] + H[0];
	    a [map[1]]   = W[1] + H[1];
	    a [map[0]+1] = W[3] + H[2];
	    a [map[1]+1] = W[4] + H[3];
#else
	    a [map[0]]   = W[0] + H[0];
	    a [map[0]+1] = W[1] + H[1];
	    a [map[1]]   = W[3] + H[2];
	    a [map[1]+1] = W[4] + H[3];
#endif

	    if (variant & NT_FIXED_POINT)
	    {
#if MPI
	      a [map[0]+2] = W[6];
	      a [map[1]+2] = W[7];
#else
	      a [map[2]]   = W[6];
	      a [map[2]+1] = W[7];
#endif
	      b [0] = T[0]*CON_RN(con) - (M[0]*R[0] + M[2]*R[1])/rho - RES[0];
	      b [1] = T[1]*CON_RN(con) - (M[1]*R[0] + M[3]*R[1])/rho - RES[1];
	    }
	    else /* HSW, HYBRID */
	    {
#if MPI
	      a [map[0]+2] = W[6] - T[0];
	      a [map[1]+2] = W[7] - T[1];
#else
	      a [map[2]]   = W[6] - T[0];
	      a [map[2]+1] = W[7] - T[1];
#endif
	      b [0] = T[0]*R[2] - (M[0]*R[0] + M[2]*R[1])/rho - RES[0];
	      b [1] = T[1]*R[2] - (M[1]*R[0] + M[3]*R[1])/rho - RES[1];
	    }

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
		a [map[0]+1] += W[3];
		a [map[1]+1] += W[4];
		a [map[0]+2] += W[6];
		a [map[1]+2] += W[7];
	      }
	    }
#else
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
#endif
	  }
	}
	else /* degenerate case => enforce homogenous tangential tractions.
		this is discussed in the Section 3.3 of the HSW paper */
	{
  ZERO_TANG:
#if MPI
	  a [map[0]]   = W[0];
	  a [map[1]]   = 0.0;
	  a [map[0]+1] = 0.0;
	  a [map[1]+1] = W[4];
	  a [map[0]+2] = 0.0;
	  a [map[1]+2] = 0.0;
#else
	  a [map[0]]   = W[0];
	  a [map[0]+1] = 0.0;
	  a [map[1]]   = 0.0;
	  a [map[1]+1] = W[4];
	  a [map[2]]   = 0.0;
	  a [map[2]+1] = 0.0;
#endif
	  b [0] = -R[0] * W[0];
	  b [1] = -R[1] * W[4];

	  /* FIXME: mind zero initialization at the beginning: hence below not needed */
#if 0
#if MPI
	  {
	    OFFB *adj [2] = {dia->adj, dia->adjext};
	    for (int i = 0; i < 2; i ++)
	    for (blk = adj [i]; blk; blk = blk->n)
	    {
	      int *map = blk->map;

	      a [map[0]]   = 0.0;
	      a [map[1]]   = 0.0;
	      a [map[0]+1] = 0.0;
	      a [map[1]+1] = 0.0;
	      a [map[0]+2] = 0.0;
	      a [map[1]+2] = 0.0;
	    }
	  }
#else
	  for (blk = dia->adj; blk; blk = blk->n)
	  {
	    int *map = blk->map;

	    a [map[0]]   = 0.0;
	    a [map[0]+1] = 0.0;
	    a [map[1]]   = 0.0;
	    a [map[1]+1] = 0.0;
	    a [map[2]]   = 0.0;
	    a [map[2]+1] = 0.0;
	  }
#endif
#endif
	}
      }
      else /* frictional sticking: inactive tangential set */
      {

        CON_STICKS(con)	 |= 0x1; /* feed in one bit */

#if MPI
	a [map[0]]   = W[0];
	a [map[1]]   = W[1];
	a [map[0]+1] = W[3];
	a [map[1]+1] = W[4];
#else
	a [map[0]]   = W[0];
	a [map[0]+1] = W[1];
	a [map[1]]   = W[3];
	a [map[1]+1] = W[4];
#endif

	if ((variant & NT_NONSMOOTH_HYBRID) ||
	    (variant & NT_FIXED_POINT))
	{
#if MPI
	  a [map[0]+2] = W[6];
	  a [map[1]+2] = W[7];
#else
	  a [map[2]]   = W[6];
	  a [map[2]+1] = W[7];
#endif
	  if (dynamic)
	  {
	    b [0] = -U[0] - V[0] - RES[0]; /* -V = W(R + dR) + B; U = WR + B; -U - V = W dR */
	    b [1] = -U[1] - V[1] - RES[1];
	  }
	  else
	  {
	    b [0] = -U[0] - RES[0];
	    b [1] = -U[1] - RES[1];
	  }
	}
	else /* HSW */
	{
	  /* TODO: this is quasi-statics; work out dynamics */
#if MPI
	  a [map[0]+2] = W[6]+U[0]/d[2];
	  a [map[1]+2] = W[7]+U[1]/d[2];
#else
	  a [map[2]]   = W[6]+U[0]/d[2];
	  a [map[2]+1] = W[7]+U[1]/d[2];
#endif
	  b [0] = -(1.0 + rho*udash/d[2])*U[0] - RES[0];
	  b [1] = -(1.0 + rho*udash/d[2])*U[1] - RES[1];
	}

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
	    a [map[0]+1] += W[3];
	    a [map[1]+1] += W[4];
	    a [map[0]+2] += W[6];
	    a [map[1]+2] += W[7];
	  }
	}
#else
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
#endif
      }

      if (isnan (b [2]))
      {
	ASSERT_DEBUG (0, "b [2] is NAN");
      }
    }
    break;
    case FIXPNT:
    {
      double *V = dia->V;

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
      if (dynamic)
      {
	b [0] = -U[0] - V[0] - RES[0];
	b [1] = -U[1] - V[1] - RES[1];
	b [2] = -U[2] - V[2] - RES[2];
      }
      else
      {
	b [0] = -U[0] - RES[0];
	b [1] = -U[1] - RES[1];
	b [2] = -U[2] - RES[2];
      }

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
      }
#endif
    }
    break;
    case FIXDIR:
    case VELODIR:
    case RIGLNK:
    {
      double *V = dia->V;

      int *map = dia->map;

#if MPI
      a [map[0]]   = W[0];
      a [map[1]]   = 0.0;
      a [map[2]]   = 0.0;
      a [map[0]+1] = 0.0;
      a [map[1]+1] = W[4];
      a [map[2]+1] = 0.0;
      a [map[0]+2] = 0.0;
      a [map[1]+2] = 0.0;
      a [map[2]+2] = W[8];
#else
      a [map[0]]   = W[0];
      a [map[0]+1] = 0.0;
      a [map[0]+2] = 0.0;
      a [map[1]]   = 0.0;
      a [map[1]+1] = W[4];
      a [map[1]+2] = 0.0;
      a [map[2]]   = 0.0;
      a [map[2]+1] = 0.0;
      a [map[2]+2] = W[8];
#endif
      b [0] = -R[0] * W[0]; /* keep RT zero: W (-R) = W dR */
      b [1] = -R[1] * W[4];

      if (dynamic)
      {
	if (con->kind == VELODIR) b [2] = -U[2] + VELODIR(con->Z) - RES[2]; /* V(t+h) = W(R+dR) + B; U = WR + B; -U+V(t+h) = WdR */
	else b [2] = -U[2] - V[2] - RES[2]; /* FIXDIR, RIGLNK */
      }
      else
      {
	if (con->kind == VELODIR) b [2] = -U[2] + VELODIR(con->Z) - RES[2];
	else if (con->kind == FIXDIR) b [2] = -U[2] - RES[2]; 
	else /* RIGLNK: see doc/notes.lyx for explanation */
	{
	  double C [3];

	  ADD (U, RES, C);

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
      {
	OFFB *adj [2] = {dia->adj, dia->adjext};
	for (int i = 0; i < 2; i ++)
	for (blk = adj [i]; blk; blk = blk->n)
	{
	  double *W = blk->W;
	  int *map = blk->map;

#if 0
	  a [map[0]]   = 0.0;
	  a [map[1]]   = 0.0;
#endif
	  a [map[2]]   += W[2]; /* maintain tangentail coupling: W(3,:) */
#if 0
	  a [map[0]+1] = 0.0;
	  a [map[1]+1] = 0.0;
#endif
	  a [map[2]+1] += W[5];
#if 0
	  a [map[0]+2] = 0.0;
	  a [map[1]+2] = 0.0;
#endif
	  a [map[2]+2] += W[8];
	}
      }
#else
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *W = blk->W;
	int *map = blk->map;

#if 0
	a [map[0]]   = 0.0;
	a [map[0]+1] = 0.0;
#endif
	a [map[0]+2] += W[2];
#if 0
	a [map[1]]   = 0.0;
	a [map[1]+1] = 0.0;
#endif
	a [map[1]+2] += W[5];
#if 0
	a [map[2]]   = 0.0;
	a [map[2]+1] = 0.0;
#endif
	a [map[2]+2] += W[8];
      }
#endif
    }
    break;
    }
  }
}

/* update linear system NT_NONSMOOTH_VARIATIONAL, NT_SMOOTHED_VARIATIONAL variants */
static void system_update_VARIATIONAL (LINSYS *sys, NTVARIANT variant, LOCDYN *ldy)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED); /* TODO */
}

/* update linear system */
static void system_update (LINSYS *sys, NTVARIANT variant, LOCDYN *ldy)
{
  switch (VARIANT (variant))
  {
    case NT_NONSMOOTH_HSW:
    case NT_NONSMOOTH_HYBRID:
    case NT_FIXED_POINT:
      system_update_HSW_HYBRID_FIXED (sys, variant, ldy);
      break;
    case NT_NONSMOOTH_VARIATIONAL:
    case NT_SMOOTHED_VARIATIONAL:
      system_update_VARIATIONAL (sys, variant, ldy);
      break;
  }
}

/* solve linear system */
#if MPI
static void system_solve (LINSYS *sys, LOCDYN *ldy, double accuracy)
{
  /* TODO: overwrite con->R with DR and use DOM_Update_External_Reactions in order to store in conext->R th DRs;
   *       then compute DU using them; note that next call to system_update will write con->R into conext->R again */
}
#else
static void system_solve (LINSYS *sys, LOCDYN *ldy, double accuracy)
{
  double *b;
  CON *con;

#if UMFPACK
  void *numeric;

  ASSERT (umfpack_di_numeric (sys->cbp, sys->cri, sys->a, sys->symbolic, &numeric, NULL, NULL) >= 0, ERR_UMFPACK_NUMERIC);
  ASSERT (umfpack_di_solve (UMFPACK_A, sys->cbp, sys->cri, sys->a, sys->x, sys->b, numeric, NULL, NULL) >= 0, ERR_UMFPACK_SOLVE);
  umfpack_di_free_numeric (&numeric);
#else

  LSS_Set (sys->lss, LSS_RESTART, 100);
  LSS_Set (sys->lss, LSS_SMOOTHING_STEPS, 3);
  LSS_Set (sys->lss, LSS_ITERATIONS_BOUND, 500);
  LSS_Set (sys->lss, LSS_RELATIVE_ACCURACY, 1E99);
  LSS_Set (sys->lss, LSS_ABSOLUTE_ACCURACY, accuracy);

  switch (LSS_Solve (sys->lss, sys->a, sys->x, sys->b))
  {
  case LSSERR_INVALID_ARGUMENT: printf ("LSS ERROR: invalid argument\n"); break;
  case LSSERR_OUT_OF_MEMORY: printf ("LSS ERROR: out of memory\n"); break;
  case LSSERR_LACK_OF_CONVERGENCE: printf ("LSS ERROR: lack of convergence\n"); break;
  case LSSERR_EMPTY_COLUMN: printf ("LSS ERROR: empty column\n"); break;
  case LSSERR_ZERO_ON_DIAGONAL: printf ("LSS ERROR: zero on diagonal\n"); break;
  case LSSERR_GMRES_BREAKDOWN: printf ("LSS ERROR: GMRES has broke down\n"); break;
  case LSSERR_NONE: break;
  }

  printf ("NEWTON: LSS error rel/abs: %g/%g and iters: %d\n", LSS_Get (sys->lss, LSS_RELATIVE_ERROR),
      LSS_Get (sys->lss, LSS_ABSOLUTE_ERROR), (int) LSS_Get (sys->lss, LSS_ITERATIONS));
#endif

  /* write DR */
  for (con = ldy->dom->con, b = sys->x; con; con = con->next, b += 3)
  {
    double *DR = CON_DR(con);

    COPY (b, DR);
  }

  /* compute DU */
  for (con = ldy->dom->con, b = sys->x; con; con = con->next, b += 3)
  {
    DIAB *dia = con->dia;

    double *DR = CON_DR(con),
	   *DU = CON_DU(con),
	   *W = dia->W;

    NVMUL (W, DR, DU);

    for (OFFB *blk = dia->adj; blk; blk = blk->n)
    {
      double *DR = CON_DR(blk->dia->con),
	     *W = blk->W;

      NVADDMUL (DU, W, DR, DU);
    }
  }
}
#endif

#if DEBUG
#if 0
static void write_system (LINSYS *sys, const char *matrix, const char *vector)
{
#if MPI
  /* TODO */
#else
  int *i, j, k;
  double *a;
  FILE *f;

  ASSERT (f = fopen (matrix, "w"), ERR_FILE_OPEN);
  fprintf (f, "%%%%MatrixMarket matrix coordinate real general\n");
  fprintf (f, "%d %d %d\n", sys->dim, sys->dim, sys->nnz);
  for (a = sys->a, i = sys->cri, j = 0; j < sys->dim; j ++)
  {
    for (k = sys->cbp [j]; k < sys->cbp [j+1]; k ++, a ++, i ++)
    {
      fprintf (f, "%d %d  %.15e\n", ((*i)+1),  (j+1),   *a);
    }
  }
  fclose (f);

  ASSERT (f = fopen (vector, "w"), ERR_FILE_OPEN);
  for (a = sys->b, j = 0; j < sys->dim; a ++, j ++)
  {
    fprintf (f, "%.15e\n", *a);
  }
  fclose (f);
#endif
}
#endif

static void test_linear_solver (LINSYS *sys, LOCDYN *ldy)
{
  double *a, *b, *x, *e, error, one [3] = {1.0, 1.0, 1.0};
  DIAB *dia;
  OFFB *blk;
  CON *con;

  for (a = sys->a, b = a + sys->nnz; a < b; a ++) *a = 0.0;
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

#if 0
    AUXDUMP ("B.vec", "%.15e\n", b [0]);
    AUXDUMP ("B.vec", "%.15e\n", b [1]);
    AUXDUMP ("B.vec", "%.15e\n", b [2]);
#endif
  }

#if !MPI && DEBUG
   {
      int j, k, *i;
      double *a;

      for (j = 0, a = sys->a, i = sys->cri; j < sys->dim; j ++)
      {
	for (k = sys->cbp [j]; k < sys->cbp [j+1]; k ++, a ++, i ++)
	{
	  if (j == *i)
	  {
	    if (*a == 0.0)
	    {
	      printf ("ERROR: diagonal zero for index %d\n", j);
	      printf ("ERROR: system dimension is %d and number of nonzeros is %d\n", sys->dim, sys->nnz);
	      printf ("ERROR: column %d begins at %d and ends at %d\n", j, sys->cbp [j], sys->cbp [j+1]);
	    }
	  }
	}
      }
    }
#endif

#if 0
  write_system (sys, "A.mtx", "B.vec");
#endif

  double accu = 1E-3;

  system_solve (sys, ldy, accu);

  for (error = 0.0, x = sys->x, e = x + sys->dim; x < e; x ++)
  {
    error += (1.0 - (*x))*(1.0 - (*x));
  }
  error = sqrt (error);

  printf ("NEWTON LINEAR SOLVER TEST: error = %g, %s\n", error, error < accu ? "PASSED" : "FAILED");
}
#endif

/* compute merit function at (R + alpha * dR) */
static double merit_function (NTVARIANT variant, LOCDYN *ldy, double alpha)
{
  double value, step;
  short dynamic;
  DIAB *dia;
  DOM *dom;
  CON *con;

  value = 0.0;
  dom = ldy->dom;
  step = dom->step;
  dynamic = dom->dynamic;
  for (con = dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT)
    {
      dia = con->dia;

      double R [3], U [3],
	     rho = dia->rho,
	     gap = con->gap,
	     fri = con->mat.base->friction,
	     res = con->mat.base->restitution,
	    *RES = CON_RES(con),
	    *DR = CON_DR(con),
	    *DU = CON_DU(con),
	    *R0 = dia->R,
	    *U0 = dia->U,
	    *V = dia->V,
	     G [3],
	     d [3],
	     bound,
	     udash,
	     norm;

      ADDMUL (R0, alpha, DR, R);
      ADDMUL (U0, alpha, DU, U);
      ADD (U, RES, U);

      if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
      else udash = ((MAX(gap, 0)/step) + U[2]);

      d [0] = R[0] - rho * U[0];
      d [1] = R[1] - rho * U[1];
      d [2] = R[2] - rho * udash;

      norm = sqrt (d[0]*d[0]+d[1]*d[1]);

      if (variant & (NT_NONSMOOTH_HSW|NT_NONSMOOTH_HYBRID|NT_FIXED_POINT))
      {
	switch (VARIANT(variant))
	{
	  case NT_NONSMOOTH_HSW:
	    bound = fri * MAX (0, d[2]);
	    break;
	  case NT_NONSMOOTH_HYBRID:
	    bound = fri * MAX (0, R[2]);
	    break;
	  case NT_FIXED_POINT:
	    bound = fri * MAX (0, CON_RN(con));
	    break;
	}

	G [0] = MAX (bound, norm)*R[0] - bound*d[0];
	G [1] = MAX (bound, norm)*R[1] - bound*d[1];
	G [2] = R [2] - MAX (0, d[2]);

	value += DOT (G,G);
      }
      else
      {
	ASSERT (0, ERR_NOT_IMPLEMENTED); /* TODO */
      }
    }
  }
  value *= 0.5;

  return value;
}

/* non-monotone globalization merit function value selection */
static double nonmonotone_merit_function (double *values, int length, double merit, int index)
{
  values [index % length] = merit;

  for (int n = 0; n < length; n ++)
    if (values [n] > merit) merit = values [n];

  return merit;
}

/* line search */
static double line_search (NTVARIANT variant, LOCDYN *ldy, double reference, double merit)
{
  double auxiliary,
	 alpha,
	 gamma;

  int imax,
      iter;

  auxiliary = 2.0 * merit;
  alpha = 1.0;
  gamma = 0.1;
  imax = 32;

  for (iter = 0; iter < imax; iter ++)
  {
    merit = merit_function (variant, ldy, alpha);

    if (merit <= (reference - gamma*alpha*auxiliary)) break;
    else alpha *= 0.9;
  }

  return alpha;
}

#if 0
/* scale regularization parameters */
static void scale_rhos (DOM *dom, double coef, char limit)
{
  CON *con;

  for (con = dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT)
    {
      if ((CON_STICKS(con) & 0x3) == 1 || (CON_STICKS(con) & 0x3) == 2)
      {
	if (CON_HITS(con) < limit)
	{
	  con->dia->rho *= coef;
	  CON_HITS(con) ++;
	}
      }
    }
  }
}
#endif

/* update constraint reactions */
static double reactions_update (LINSYS *sys, NTVARIANT variant, LOCDYN *ldy, int iter, int length, double *values)
{
  double merit, reference, alpha, errup, errlo;
  CON *con;
  DOM *dom;

  dom = ldy->dom;

  if (iter && VARIANT(variant) != NT_FIXED_POINT) /* TODO: the old code was skipping this for iter == 0 (test it) */
  {
    merit = merit_function (variant, ldy, 0.0);

    reference = nonmonotone_merit_function (values, length, merit, iter);

    alpha = line_search (variant, ldy, reference, merit);
  }
  else alpha = 1.0;

  /* R = R + DR, U = U + DU */
  errup = errlo = 0.0;
  for (con = dom->con; con; con = con->next)
  {
    double *RES = CON_RES(con),
           *DR = CON_DR(con),
	   *DU = CON_DU(con),
	   *R = con->R,
	   *U = con->U;

    ADDMUL (R, alpha, DR, R);
    ADDMUL (U, alpha, DU, U);
    ADD (U, RES, U);

    errup += alpha*alpha*DOT(DR, DR);
    errlo += DOT(R, R);
  }

#if 0
  /* scale dia->rho parameters */
  scale_rhos (dom, 10.0, 8);
#endif

  return sqrt (errup) / sqrt (MAX (errlo, 1.0));
}

/* update noraml bounds and return relative error of the update */
static double update_normal_bounds (LOCDYN *ldy)
{
  double errup, errlo;
  CON *con;

  errup = errlo = 0.0;
  for (con = ldy->dom->con; con; con = con->next)
  {
    if (con->kind == CONTACT)
    {
      double DRN,
	     RN;

      DRN = con->R[2] - CON_RN(con);
      RN = CON_RN(con) = con->R[2];

      errup += DRN*DRN;
      errlo += RN*RN;
    }
  }

  return sqrt (errup) / MAX (sqrt (errlo), 1.0);
}

/* destroy linear system */
#if MPI
static void system_destroy (LINSYS *sys, MEM *mapmem, MEM *auxmeme, LOCDYN *ldy)
{
  /* TODO */
}
#else
static void system_destroy (LINSYS *sys, MEM *mapmem, MEM *auxmem, LOCDYN *ldy)
{
  DIAB *dia;
  OFFB *blk;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    dia->con->data = NULL;

    dia->map = NULL;

    for (blk = dia->adj; blk; blk = blk->n)
    {
      blk->map = NULL;
    }
  }

  MEM_Release (mapmem);
  MEM_Release (auxmem);

  free (sys->a);
  free (sys->x);
  free (sys->b);
  free (sys->cri);
  free (sys->cbp);

#if UMFPACK
  umfpack_di_free_symbolic (&sys->symbolic);
#else
  LSS_Destroy (sys->lss);
#endif

  free (sys);
}
#endif

/* create solver */
NEWTON* NEWTON_Create (NTVARIANT variant, double epsilon, int maxiter)
{
  NEWTON *nt;

  ERRMEM (nt = MEM_CALLOC (sizeof (NEWTON)));

  nt->variant = variant;
  nt->epsilon = epsilon;
  nt->maxiter = maxiter;

  MEM_Init (&nt->mapmem, sizeof (int [3]), MEMBLK);
  MEM_Init (&nt->auxmem, sizeof (CONAUX), MEMBLK);

  nt->length = 10;

  return nt;
}

/* run solver */
void NEWTON_Solve (NEWTON *nt, LOCDYN *ldy)
{
  double merit, error, *values;
  char fmt [512];
  LINSYS *sys;
  int iter;
  DOM *dom;

  dom = ldy->dom;

  if (dom->verbose) sprintf (fmt, "NEWTON: iteration: %%%dd  error:  %%.2e  merit: %%.2e\n", (int)log10 (nt->maxiter) + 1);

  ERRMEM (values = MEM_CALLOC (sizeof (double [nt->length])));

  sys = system_create (&nt->mapmem, &nt->auxmem, ldy);

#if DEBUG
  test_linear_solver (sys, ldy);
#endif

  error = merit = 1.0;

  iter = 0;

  do
  {
    system_update (sys, nt->variant, ldy); /* assemble A, b */

#if 0
    write_system (sys, "A.mtx", "B.vec");
#endif

    system_solve (sys, ldy, MIN (error * 0.01, 1E-2)); /* solve x = inv (A) * b  */

    error = reactions_update (sys, nt->variant, ldy, iter, nt->length, values); /* R[iter+1] = R[iter] (x) */

    if ((nt->variant & NT_FIXED_POINT) && error < nt->epsilon)
    {
      error = update_normal_bounds (ldy);
    }

    merit = merit_function (nt->variant, ldy, 0.0);

    if (dom->verbose) printf (fmt, iter, error, merit);

  } while (error > nt->epsilon && ++ iter < nt->maxiter);

  system_destroy (sys, &nt->mapmem, &nt->auxmem, ldy);

  free (values);
}

/* write labeled satate values */
void NEWTON_Write_State (NEWTON *nt, PBF *bf)
{
}

/* destroy solver */
void NEWTON_Destroy (NEWTON *nt)
{
  MEM_Release (&nt->mapmem);
  MEM_Release (&nt->auxmem);
  free (nt);
}
