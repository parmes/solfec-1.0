/*
 * pbs.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * per-body constraints solver
 */

#include <stdlib.h>
#include <time.h>
#include "dom.h"
#include "bgs.h"
#include "pbs.h"
#include "alg.h"
#include "mem.h"
#include "err.h"

#if MPI
#include "com.h"
#include "pck.h"
#endif

/* Gauss-Seidel per-body solver */
static void per_body_solver (SET *conset, double epsilon, int maxiter, short dynamic, double step)
{
  double R0 [3], B [3], *R, errup, errlo, error, diagepsilon;
  int iters, diagiters, diagmaxiter;
  GSDIAS diagsolver;
  SET *item;
  CON *con;
  OFFB *blk;
  DIAB *dia;

  diagsolver = GS_PROJECTED_GRADIENT;
  diagepsilon = 1E-6;
  diagmaxiter = 100;
  iters = 0;
  do
  {
    errup = errlo = 0.0;
  
    for (item = SET_First (conset); item; item = SET_Next (item)) 
    {
      con = item->data;
      dia = con->dia;
      R = dia->R;

      /* compute local free velocity */
      COPY (dia->B, B);
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *W = blk->W,
	       *R = blk->dia->R;

	NVADDMUL (B, W, R, B);
      }

#if MPI
      for (blk = dia->adjext; blk; blk = blk->n)
      {
	double *R = ((CON*)blk->dia)->R,
               *W = blk->W;

	NVADDMUL (B, W, R, B);
      }
#endif
      
      COPY (R, R0); /* previous reaction */

      /* solve local diagonal block problem */
      diagiters = DIAGONAL_BLOCK_Solver (diagsolver, diagepsilon, diagmaxiter, dynamic,
	           step, con->kind, con->mat.base, con->gap, con->Z, con->base, dia, B);

      if (diagiters >= diagmaxiter || diagiters < 0) /* failed */
      {
	if (con->kind == CONTACT)
	{
	  GSDIAS dias [3] = {GS_SEMISMOOTH_NEWTON, GS_PROJECTED_GRADIENT, GS_DE_SAXE_AND_FENG};

	  for (int i = 0; i < 3; i ++)
	  {
	    if (dias [i] != diagsolver) /* skip current diagonal solver */
	    {
	      COPY (R0, R); /* initialize with previous reaction */

	      diagiters = DIAGONAL_BLOCK_Solver (dias [i], diagepsilon, diagmaxiter, /* try another solver */
		dynamic, step, con->kind, con->mat.base, con->gap, con->Z, con->base, dia, B);

	      if (diagiters < diagmaxiter && diagiters >= 0) break; /* success */
	    }
	  }
	}

	if (diagiters >= diagmaxiter || diagiters < 0) /* failed */
	{
	  COPY (R0, R); /* use previous reaction */
	}
      }

      /* accumulate relative
       * error components */
      SUB (R, R0, R0);
      errup += DOT (R0, R0);
      errlo += DOT (R, R);
    }

    /* calculate relative error */
    error = sqrt (errup) / sqrt (MAX (errlo, 1.0));
  }
  while (++ iters < maxiter && error > epsilon);
}

#if MPI
/* pack data structure */
typedef struct pack_data { int rank; SET *set; } PD;

/* pack exported local dynamics rows */
static void local_dynamics_pack (PD *pd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  double B [3];
  SET *item;
  CON *con;
  DIAB *dia;
  OFFB *blk;
  int n;

  pack_int (isize, i, ints, SET_Size (pd->set));

  for (item = SET_First (pd->set); item; item = SET_Next (item))
  {
    con = item->data;

    pack_int (isize, i, ints, con->id);

    switch ((int) con->kind)
    {
      case CONTACT:
      SURFACE_MATERIAL_Pack_State (&con->mat, dsize, d, doubles, isize, i, ints);
      break;
      case VELODIR:
      case RIGLNK:
      pack_doubles (dsize, d, doubles, con->Z, DOM_Z_SIZE);
      break;
    }

    pack_double (dsize, d, doubles, con->gap);

    dia = con->dia;

    COPY (dia->B, B);

    for (n = 0, blk = dia->adjext; blk; blk = blk->n)
    {
      CON *con = (CON*) blk->dia; /* external blocks point directly to external contacts <= ldy.c:compute_adjext */

      if (con->rank == pd->rank) n ++; /* if the destination rank is the parent rank of this external constraint: export W */
      else /* otherwise contribute to B only */
      {
	double *W = blk->W,
	       *R = con->R;

	NVADDMUL (B, W, R, B);
      }
    }
    for (blk = dia->adj; blk; blk = blk->n)
    {
      CON *con = blk->dia->con;

      if (SET_Contains (con->ext, (void*) (long) pd->rank, NULL)) n ++; /* if the constraint was exported to the destination rank: export W */
      else /* otherwise contribute to B only */
      {
	double *W = blk->W,
	       *R = con->R;

	NVADDMUL (B, W, R, B);
      }
    }

    pack_int (isize, i, ints, n);

    for (blk = dia->adjext; blk; blk = blk->n)
    {
      CON *con = (CON*) blk->dia; /* external blocks point directly to external contacts <= ldy.c:compute_adjext */

      if (con->rank == pd->rank)
      {
	pack_int (isize, i, ints, con->id);
	pack_doubles (dsize, d, doubles, blk->W, 9);
	pack_doubles (dsize, d, doubles, con->R, 3);
      }
    }

    for (blk = dia->adj; blk; blk = blk->n)
    {
      CON *con = blk->dia->con;

      if (SET_Contains (con->ext, (void*) (long) pd->rank, NULL))
      {
	pack_int (isize, i, ints, con->id);
	pack_doubles (dsize, d, doubles, blk->W, 9);
	pack_doubles (dsize, d, doubles, con->R, 3);
      }
    }

    pack_doubles (dsize, d, doubles, dia->V, 3);
    pack_doubles (dsize, d, doubles, B, 3);
    pack_doubles (dsize, d, doubles, dia->W, 9);
    pack_double (dsize, d, doubles, dia->rho);
  }
}

/* unpack imported local dynamics rows */
static void* local_dynamics_unpack (LOCDYN *ldy, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  DOM *dom = ldy->dom;
  int k, n, j, m, id;
  double R [3];
  CON *con;
  DIAB *dia;
  OFFB *blk;

  n = unpack_int (ipos, i, ints);

  for (k = 0; k < n; k ++)
  {
    id = unpack_int (ipos, i, ints);

    ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid constraint id");

    switch ((int) con->kind)
    {
      case CONTACT:
      SURFACE_MATERIAL_Unpack_State (dom->sps, &con->mat, dpos, d, doubles, ipos, i, ints);
      break;
      case VELODIR:
      case RIGLNK:
      unpack_doubles (dpos, d, doubles, con->Z, DOM_Z_SIZE);
      break;
    }

    con->gap = unpack_double (dpos, d, doubles);

    ERRMEM (dia = MEM_Alloc (&ldy->diamem));
    dia->R = con->R;
    dia->U = con->U;
    con->dia = dia;

    m = unpack_int (ipos, i, ints);

    for (j = 0; j < m; j ++)
    {
      id = unpack_int (ipos, i, ints);

      CON *con = MAP_Find (dom->idc, (void*) (long) id, NULL);
      if (!con) { ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid constraint id"); }

      ERRMEM (blk = MEM_Alloc (&ldy->offmem));
      unpack_doubles (dpos, d, doubles, blk->W, 9);
      unpack_doubles (dpos, d, doubles, R, 3);

      blk->dia = (DIAB*) con; /* use the diagonal block pointer */

      blk->n = dia->adjext; /* append external adjacency list */
      dia->adjext = blk;
    }

    unpack_doubles (dpos, d, doubles, dia->V, 3);
    unpack_doubles (dpos, d, doubles, dia->B, 3);
    unpack_doubles (dpos, d, doubles, dia->W, 9);
    dia->rho = unpack_double (dpos, d, doubles);
  }

  return NULL;
}

/* pack boundary reactions */
static void boundary_reactions_pack (PD *pd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SET *item;
  CON *con;

  pack_int (isize, i, ints, SET_Size (pd->set));

  for (item = SET_First (pd->set); item; item = SET_Next (item))
  {
    con = item->data;

    pack_int (isize, i, ints, con->id);
    pack_doubles (dsize, d, doubles, con->R, 3);
  }
}

/* unpack and average boundary reactions */
static void* boundary_reactions_unpack (LOCDYN *ldy, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  DOM *dom = ldy->dom;
  double R [3];
  int j, n, id;
  CON *con;

  n = unpack_int (ipos, i, ints);

  for (j = 0; j < n; j ++)
  {
    id = unpack_int (ipos, i, ints);
    unpack_doubles (dpos, d, doubles, R, 3);

    con = MAP_Find (dom->idc, (void*) (long) id, NULL);
    if (!con) { ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) id, NULL), "Invalid constraint id"); }

    MID (con->R, R, con->R);
  }

  return NULL;
}
#endif

/* create solver */
PER_BODY* PER_BODY_Create (double epsilon, int maxiter)
{
  PER_BODY *pb;

  ERRMEM (pb = malloc (sizeof (PER_BODY)));

  pb->epsilon = epsilon;
  pb->maxiter = maxiter;

  return pb;
}

/* run solver */
void PER_BODY_Solve (PER_BODY *pb, LOCDYN *ldy)
{
  int maxiter = pb->maxiter;
  double epsilon = pb->epsilon;
  DOM *dom = ldy->dom;
  short dynamic = dom->dynamic;
  double step = dom->step;
  BODY *bod, **perm;
  int i, j, n;
#if MPI
  SET **conext, *item, *boundary, *fixed;
  MEM *setmem = &dom->setmem;
  COMDATA *send, *recv, *ptr;
  int nsend, nrecv, ncpu;
  COMOBJ *osend, *orecv;
  CON *con;
  PD *pd;

  ncpu = dom->ncpu;

  boundary = fixed = NULL;

  ERRMEM (conext = MEM_CALLOC (ncpu * sizeof (SET*))); /* rank sets of external constraints used by parent bodies */

  for (bod = dom->bod; bod; bod = bod->next) /* find external constraints attached to parent bodies */
  {
    if (bod->kind == OBS) continue;

    for (item = SET_First (bod->con); item; item = SET_Next (item))
    {
      con = item->data;

      if (!con->dia) /* external */
      {
	SET_Insert (setmem, &conext [con->rank], con, NULL); /* note that dom->conext, although contains all external constraints, might be
						                notabely larger as it includes constraints attached to children and parents */
      }

      if (!con->dia || con->ext) SET_Insert (setmem, &boundary, bod, NULL); /* mark as a boundary body */

      if (con->kind != CONTACT) SET_Insert (setmem, &fixed, con, NULL); /* mark as a fixed constraint */
    }
  }

  ERRMEM (send = MEM_CALLOC (ncpu * sizeof (COMDATA)));

  for (nsend = i = 0; i < ncpu; i ++)
  {
    if (conext [i])
    {
      send [nsend].ints = SET_Size (conext [i]);
      ERRMEM (send [nsend].i = malloc (send [nsend].ints * sizeof (int)));
      for (item = SET_First (conext [i]), j = 0; item; item = SET_Next (item), j ++)
      {
	con = item->data;
	send [nsend].i [j] = con->id;
      }
      send [nsend].rank = i;
      nsend ++;
    }
  }

  ERRMEM (osend = MEM_CALLOC (ncpu * sizeof (COMOBJ)));
  ERRMEM (pd = MEM_CALLOC (ncpu * sizeof (PD)));
  for (i = 0; i < ncpu; i ++) 
  {
    pd [i].rank = osend [i].rank = i;
    osend [i].o = &pd [i];
  }

  /* send ids of needed external constraints
   * (local dynamics rows) to parent ranks */
  COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);

  for (ptr = recv; ptr < recv + nrecv; ptr ++)
  {
    for (j = 0; j < ptr->ints; j ++)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*) (long) ptr->i [j], NULL), "Invalid constraint id");
      ASSERT_DEBUG (con->dia, "Inconsistent constraints: external constraint in place of parent one");
      SET_Insert (setmem, &pd [ptr->rank].set, con, NULL); /* schedule for sending (local dynamics row) */
    }
  }

  /* send local dynamics rows */
  COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack) local_dynamics_pack, ldy, (OBJ_Unpack) local_dynamics_unpack, osend, ncpu, &orecv, &nrecv);

  /* solve boundary bodies */
  for (item = SET_First (boundary); item; item = SET_Next (item))
  {
    bod = item->data;
    per_body_solver (bod->con, epsilon, maxiter, dynamic, step);
  }

  /* create complete send sets comprizing both the exporting and the importing constraints */
  for (i = 0; i < ncpu; i ++)
  {
    for (item = SET_First (conext [i]); item; item = SET_Next (item))
    {
      SET_Insert (setmem, &pd [i].set, item->data, NULL);
    }
  }
  free (orecv);

  /* negotiate boundary reactions */
  COMOBJSALL (MPI_COMM_WORLD, (OBJ_Pack) boundary_reactions_pack, ldy, (OBJ_Unpack) boundary_reactions_unpack, osend, ncpu, &orecv, &nrecv);
#endif

  srand (time (NULL));

  ERRMEM (perm = malloc (dom->nbod * sizeof (BODY*)));

  for (n = 0, bod = dom->bod; bod; bod = bod->next)
  {
#if MPI
    if (bod->kind != OBS && !SET_Contains (boundary, bod, NULL)) /* only inner bodies here */
#else
    if (bod->kind != OBS)
#endif
    {
      perm [n ++] = bod; /* initial permutation */
    }
  }

  for (i = 0; i < n; i ++)
  {
    j = rand () % n; /* permute bodies */
    bod = perm [j];
    perm [j] = perm [i];
    perm [i] = bod;
  }

  for (i = 0; i < n; i ++) /* run sequence of body solvers */
  {
    per_body_solver (perm [i]->con, epsilon, maxiter, dynamic, step);
  }

  free (perm);

#if MPI
  /* additionally iterate on bilateral constraints in order to ensure accuracy */
  per_body_solver (fixed, epsilon, maxiter, dynamic, step);

  /* clean up */
  for (i = 0; i < dom->ncpu; i ++)
  {
    if (conext [i])
    {
      for (item = SET_First (conext [i]); item; item = SET_Next (item))
      {
	OFFB *blk, *bne;

	con = item->data;

	/* free imported local dynamics rows */
	ASSERT_DEBUG (con->dia, "Inconsistent constraint data: missing imported row");
	for (blk = con->dia->adjext; blk; blk = bne)
	{
	  bne = blk->n;
	  MEM_Free (&ldy->offmem, blk);
	}
	MEM_Free (&ldy->diamem, con->dia);
	con->dia = NULL;
      }

      SET_Free (setmem, &conext [i]);
    }

    SET_Free (setmem, &pd [i].set);
    if (i < nsend) free (send [i].i);
  }
  SET_Free (setmem, &boundary);
  SET_Free (setmem, &fixed);
  free (conext);
  free (osend);
  free (orecv);
  free (send);
  free (recv);
  free (pd);
#endif
}

/* write labeled satate values */
void PER_BODY_Write_State (PER_BODY *pb, PBF *bf)
{
}

/* destroy solver */
void PER_BODY_Destroy (PER_BODY *pb)
{
  free (pb);
}
