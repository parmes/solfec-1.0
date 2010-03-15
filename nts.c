/*
 * nts.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Newton constraints solver
 */

#if MPI
#include <HYPRE_parcsr_mv.h>
#include <HYPRE_IJ_mv.h>
#else
#include <slu_ddefs.h>
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
};
#else
struct linsys /* SuperLU matrix data and the rest of linear system */
{
  double *a; /* compressed column values */

  int *cri, /* column row indices */
      *cbp, /* column beginning pointers */
       nnz, /* number of nonzero entries */
       dim; /* dimension of the rectangular matrix */

  SuperMatrix A, /* SuperLU wrapper for 'a' */
	      AC; /* 'A' with permuted columns */

  double *b; /* right hand side */

  SuperMatrix B; /* SuperLU wrapper for 'b' */

  int *rowperm, /* permutation mappings */
      *colperm,
      *etree; /* elimination tree */

  superlu_options_t options;  /* SuperLU options */
  
  SuperMatrix L, U; /* numeric factors */
  void *work; /* work space */
  int lwork; /* size of work space */

  SuperLUStat_t stat; /* SuperLU statistics */
  int info; /* termination status */

  double *x; /* solution == b (overwrites) */
};
#endif

typedef struct conaux CONAUX; /* auxiliary constraint data */

struct conaux
{
  double RES [3], /* residual */
	 DR [3],  /* reaction increment */
	 DU [3],  /* velocity increment */
	 RN;      /* normal reaction bound during fixed point iterations */
};

#define CON_RES(con) ((CONAUX*)((con)->data))->RES
#define CON_DR(con) ((CONAUX*)((con)->data))->DR
#define CON_DU(con) ((CONAUX*)((con)->data))->DU
#define CON_RN(con) ((CONAUX*)((con)->data))->RN

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
  LINSYS *sys;
  DIAB *dia;
  OFFB *blk;
  CON *con;
    
  ERRMEM (sys = MEM_CALLOC (sizeof (LINSYS)));

  /* number constraints */
  DOM_Number_Constraints (dom);

  /* allocate auxiliary space */
  for (con = dom->con; con; con = con->next)
  {
    ERRMEM (con->data = MEM_Alloc (auxmem));
  }
  
  /* index diagonal entries and calculate number of nonzeros */
  for (sys->nnz = 0, dia = ldy->dia; dia; dia = dia->n)
  {
    sys->nnz += 9; /* 3x3 block space */
    for (blk = dia->adj; blk; blk = blk->n) sys->nnz += 9;
  }
  sys->dim = dom->ncon * 3; /* actual dimension */

  /* eallocate compressed column storage */
  ERRMEM (sys->a = malloc (sizeof (double) * sys->nnz));
  ERRMEM (sys->cri = malloc (sizeof (int) * sys->nnz));
  ERRMEM (sys->cbp = MEM_CALLOC (sizeof (int) * (sys->dim+1))); /* '+ 1' as there is cbp[0] == 0 always */

  /* eallocate permutation space and elimination tree */
  ERRMEM (sys->rowperm = MEM_CALLOC (sizeof (int) * sys->dim));
  ERRMEM (sys->colperm = malloc (sizeof (int) * sys->dim));
  ERRMEM (sys->etree = malloc (sizeof (int) * sys->dim));

  int *cri = sys->cri,
      *cbp = sys->cbp,
      *aux = sys->rowperm; /* used for temporary row entry indexing */

  /* count number of row entries in each column */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    int col;
   
    col = dia->con->num * 3;
    cbp [col+1] += 3; /* note the '+ 1' shift here => cbp[0] == 0 always */
    cbp [col+2] += 3;
    cbp [col+3] += 3;
    
    for (blk = dia->adj; blk; blk = blk->n)
    {
      col = blk->dia->con->num * 3;
      cbp [col+1] += 3;
      cbp [col+2] += 3;
      cbp [col+3] += 3;
    }
  }

  /* sum up column nonzeros in order to get the final structure  of 'cbp' */
  for (int n = 1; n <= sys->dim; n ++) cbp [n] += cbp [n-1];

  /* index row entries in each column */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    ERRMEM (dia->map = MEM_Alloc (mapmem));
    int col, row, *map = dia->map;

    /* diagonal entry */
    col = row = dia->con->num * 3;
    map [0] = cbp [col] + aux[col];
    map [1] = cbp [col+1] + aux[col+1];
    map [2] = cbp [col+2] + aux[col+2];
    cri [map[0]] = row;
    cri [map[0]+1] = row + 1;
    cri [map[0]+2] = row + 2;
    cri [map[1]] = row;
    cri [map[1]+1] = row + 1;
    cri [map[1]+2] = row + 2;
    cri [map[2]] = row;
    cri [map[2]+1] = row + 1;
    cri [map[2]+2] = row + 2;
    aux [col] += 3;
    aux [col+1] += 3;
    aux [col+2] += 3;
   
    /* off-diagonals */ 
    for (blk = dia->adj; blk; blk = blk->n)
    {
      ERRMEM (blk->map = MEM_Alloc (mapmem));
      int *map = blk->map;

      col = blk->dia->con->num * 3;
      map [0] = cbp [col] + aux[col];
      map [1] = cbp [col+1] + aux[col+1];
      map [2] = cbp [col+2] + aux[col+2];
      cri [map[0]] = row;
      cri [map[0]+1] = row + 1;
      cri [map[0]+2] = row + 2;
      cri [map[1]] = row;
      cri [map[1]+1] = row + 1;
      cri [map[1]+2] = row + 2;
      cri [map[2]] = row;
      cri [map[2]+1] = row + 1;
      cri [map[2]+2] = row + 2;
      aux [col] += 3;
      aux [col+1] += 3;
      aux [col+2] += 3;
    }
  }

  /* create SuperLU wrapper of tangent operator sys->A */
  dCreate_CompCol_Matrix (&sys->A, sys->dim, sys->dim, sys->nnz, sys->a, sys->cri, sys->cbp, SLU_NC, SLU_D, SLU_GE);

  /* realloc and prepare right hand side vector sys->B */
  ERRMEM (sys->b = realloc (sys->b, sizeof (double) * sys->dim));
  dCreate_Dense_Matrix (&sys->B, sys->dim, 1, sys->b, sys->dim, SLU_DN, SLU_D, SLU_GE); 

  /* set up factorization options */
  superlu_options_t *opt = &sys->options;

  /* may seem redundant, but is safer with respect
   * to option changes in different versions of SuperLU */
  set_default_options (opt);
 
  /* now do it more consciously
   * on a subset of options */ 
  opt->Fact = DOFACT; 
  opt->Equil = NO;
  opt->ColPerm = MMD_AT_PLUS_A;
  opt->Trans = NOTRANS;
  opt->IterRefine = NOREFINE;
  opt->PrintStat = YES;
  opt->DiagPivotThresh = 1.0;
  opt->PivotGrowth = NO;
  opt->ConditionNumber = NO;

  /* permute columns of matrix A and make it into AC */
  get_perm_c (opt->ColPerm, &sys->A, sys->colperm); /* fill-in minimising column ordering */
  sp_preorder (opt, &sys->A, sys->colperm, sys->etree, &sys->AC); /* elimination tree & AC */

  /* initialise statistics placeholder */
  StatInit (&sys->stat);

  /* solution vector */
  sys->x = sys->b;

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

	    a [map[2]]   = W[2];
	    a [map[2]+1] = W[5];
	    a [map[2]+2] = W[8];
	  }
	}
#else
	for (blk = dia->adj; blk; blk = blk->n)
	{
	  double *W = blk->W;
	  int *map = blk->map;

	  a [map[0]+2] = W[2];
	  a [map[1]+2] = W[5];
	  a [map[2]+2] = W[8];
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
      }

      /* tangential response */

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

		a [map[0]]   = rho*(M[0]*W[0] + M[2]*W[1]);
		a [map[1]]   = rho*(M[1]*W[0] + M[3]*W[1]);
		a [map[0]+1] = rho*(M[0]*W[3] + M[2]*W[4]);
		a [map[1]+1] = rho*(M[1]*W[3] + M[3]*W[4]);
		a [map[0]+2] = rho*(M[0]*W[6] + M[2]*W[7]); 
		a [map[1]+2] = rho*(M[1]*W[6] + M[3]*W[7]);  
	      }
	    }
#else
	    for (blk = dia->adj; blk; blk = blk->n)
	    {
	      double *W = blk->W;
	      int *map = blk->map;

	      a [map[0]]   = rho*(M[0]*W[0] + M[2]*W[1]);
	      a [map[0]+1] = rho*(M[1]*W[0] + M[3]*W[1]);
	      a [map[1]]   = rho*(M[0]*W[3] + M[2]*W[4]);
	      a [map[1]+1] = rho*(M[1]*W[3] + M[3]*W[4]);
	      a [map[2]]   = rho*(M[0]*W[6] + M[2]*W[7]); 
	      a [map[2]+1] = rho*(M[1]*W[6] + M[3]*W[7]);  
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

		a [map[0]]   = W[0];
		a [map[1]]   = W[1];
		a [map[0]+1] = W[3];
		a [map[1]+1] = W[4];
		a [map[0]+2] = W[6];
		a [map[1]+2] = W[7];
	      }
	    }
#else
	    for (blk = dia->adj; blk; blk = blk->n)
	    {
	      double *W = blk->W;
	      int *map = blk->map;

	      a [map[0]]   = W[0];
	      a [map[0]+1] = W[1];
	      a [map[1]]   = W[3];
	      a [map[1]+1] = W[4];
	      a [map[2]]   = W[6];
	      a [map[2]+1] = W[7];
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
	}
      }
      else /* frictional sticking */
      {
	/* inactive tangential set */
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

	    a [map[0]]   = W[0];
	    a [map[1]]   = W[1];
	    a [map[0]+1] = W[3];
	    a [map[1]+1] = W[4];
	    a [map[0]+2] = W[6];
	    a [map[1]+2] = W[7];
	  }
	}
#else
	for (blk = dia->adj; blk; blk = blk->n)
	{
	  double *W = blk->W;
	  int *map = blk->map;

	  a [map[0]]   = W[0];
	  a [map[0]+1] = W[1];
	  a [map[1]]   = W[3];
	  a [map[1]+1] = W[4];
	  a [map[2]]   = W[6];
	  a [map[2]+1] = W[7];
	}
#endif
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
      }
#else
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *W = blk->W;
	int *map = blk->map;

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

	  a [map[0]]   = 0.0;
	  a [map[1]]   = 0.0;
	  a [map[2]]   = W[2]; /* maintain tangentail coupling: W(3,:) */
	  a [map[0]+1] = 0.0;
	  a [map[1]+1] = 0.0;
	  a [map[2]+1] = W[5];
	  a [map[0]+2] = 0.0;
	  a [map[1]+2] = 0.0;
	  a [map[2]+2] = W[8];
	}
      }
#else
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *W = blk->W;
	int *map = blk->map;

	a [map[0]]   = 0.0;
	a [map[0]+1] = 0.0;
	a [map[0]+2] = W[2];
	a [map[1]]   = 0.0;
	a [map[1]+1] = 0.0;
	a [map[1]+2] = W[5];
	a [map[2]]   = 0.0;
	a [map[2]+1] = 0.0;
	a [map[2]+2] = W[8];
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
static void system_solve (LINSYS *sys, LOCDYN *ldy)
{
  /* TODO: overwrite con->R with DR and use DOM_Update_External_Reactions in order to store in conext->R th DRs;
   *       then compute DU using them; note that next call to system_update will write con->R into conext->R again */
}
#else
static void system_solve (LINSYS *sys, LOCDYN *ldy)
{
  SuperMatrix B;
  double *b;
  CON *con;

  ERRMEM (b = malloc (sizeof (double) * sys->dim));
  dCreate_Dense_Matrix (&B, sys->dim, 1, b, sys->dim, SLU_DN, SLU_D, SLU_GE); 

  if (sys->options.Fact == DOFACT) /* new */
  {
    /* symbolic factorization was done already: point it out
     * and go for the numeric factorization with row pivoting */

    sys->options.Fact = SamePattern; /* no column ordering */

    dgstrf (&sys->options, &sys->AC, sp_ienv (2), sp_ienv (1), sys->etree, NULL, -1,
      sys->colperm, sys->rowperm, &sys->L, &sys->U, &sys->stat, &sys->info); /* guess workspace size */
#if 0
    /* below after fixes internal SuperLU bug with the misassesment of space for small matrices */
    sys->lwork = sys->info + sp_ienv(1)*((sp_ienv(3)+sp_ienv(4)) + sys->dim) * sizeof (double);
    sys->lwork *= 2; /* just in case (from previous experience) */
#else
    sys->lwork = sys->info * sizeof (double); /* TODO: remove the #if after the current SuperLU proves rotbust in this respect */
#endif
    ERRMEM (sys->work = MEM_CALLOC (sys->lwork));
  }
  else /* consecutive */
  {
    /* to save time a factorization flag is changed here in order to
     * prevent row pivoting for consecutive solves (somewhat risky) */ 

    sys->options.Fact = SamePattern_SameRowPerm; /* no row pivoting */
  }

  /* update copy of b */
  blas_dcopy (sys->dim, sys->b, 1, b, 1);

  dgstrf (&sys->options, &sys->AC, sp_ienv (2), sp_ienv (1), sys->etree, sys->work, sys->lwork,
    sys->colperm, sys->rowperm, &sys->L, &sys->U, &sys->stat, &sys->info); /* numeric factorisation */
  ASSERT (sys->info == 0, ERR_SUPERLU_SOLVE);

  dgstrs (NOTRANS, &sys->L, &sys->U, sys->colperm, sys->rowperm, &sys->B, &sys->stat, &sys->info); /* system solve */
  ASSERT (sys->info == 0, ERR_SUPERLU_SOLVE);

  /* iterative refinement */
  {
    double ferr, berr;
    char equed = 'N';

    dgsrfs (NOTRANS, &sys->A, &sys->L, &sys->U, sys->colperm, sys->rowperm,
      &equed, NULL, NULL, &B, &sys->B, &ferr, &berr, &sys->stat, &sys->info);
    ASSERT (sys->info == 0, ERR_SUPERLU_SOLVE);
    free (b);
  }

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

/* update constraint reactions */
static double reactions_update (LINSYS *sys, NTVARIANT variant, LOCDYN *ldy, int iter, int length, double *values)
{
  double merit, reference, alpha, errup, errlo;
  CON *con;
  DOM *dom;

  dom = ldy->dom;

  if (VARIANT(variant) != NT_FIXED_POINT) /* TODO: the old code was skipping this for iter == 0 (test it) */
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
    MEM_Free (mapmem, dia->map);

    MEM_Free (auxmem, dia->con->data);

    for (blk = dia->adj; blk; blk = blk->n)
    {
      MEM_Free (mapmem, blk->map);
    }
  }

  free (sys->a);
  free (sys->cri);
  free (sys->cbp);
  free (sys->rowperm);
  free (sys->colperm);
  free (sys->etree);
  free (sys->b);
  free (sys->work);
  SUPERLU_FREE (sys->L.Store);
  SUPERLU_FREE (sys->U.Store);
  Destroy_CompCol_Permuted (&sys->AC);
  free (sys);
}
#endif

/* create solver */
NEWTON* NEWTON_Create (NTVARIANT variant, double epsilon, int maxiter)
{
  NEWTON *nt;

  ERRMEM (nt = malloc (sizeof (NEWTON)));

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
  double error, *values;
  char fmt [512];
  LINSYS *sys;
  int iter;
  DOM *dom;

  dom = ldy->dom;

  if (dom->verbose) sprintf (fmt, "NEWTON: iteration: %%%dd  error:  %%.2e\n", (int)log10 (nt->maxiter) + 1);

  ERRMEM (values = MEM_CALLOC (sizeof (double [nt->length])));

  sys = system_create (&nt->mapmem, &nt->auxmem, ldy);

  iter = 0;

  do
  {
    system_update (sys, nt->variant, ldy); /* assemble A, b */

    system_solve (sys, ldy); /* solve x = inv (A) * b  */

    error = reactions_update (sys, nt->variant, ldy, iter, nt->length, values); /* R[iter+1] = R[iter] (x) */

    if ((nt->variant & NT_FIXED_POINT) && error < nt->epsilon)
    {
      error = update_normal_bounds (ldy);
    }

    if (dom->verbose) printf (fmt, iter, error);

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
