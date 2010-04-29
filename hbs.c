/*
 * nts.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Hybrid constraints solver
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

#include "hbs.h"
#include "alg.h"
#include "bla.h"
#include "err.h"
#include "bgs.c"
#include "glu.h"

#define ENABLE 0 /* temporary flag before the code compiles */

#if ENABLE
typedef struct body_data BD;

struct body_data
{
  BODY *bod;

  CON **con; /* copy of bod->con */

  MX **H; /* corresponding H operators */

  int n; /* size of bod->con */

  MAP *contoi; /* maps CON* to an index in con[] */

  double *b; /* free velocity plus previous tangential contact response */

  double *u; /* free velocity projected on the constraints tangent set */

  BD *next;
};

typedef struct body_data_pair BDP;

struct body_data_pair
{
  BD *bd [2]; /* two body data sets */

  int **con [2]; /* indices in bd[i]->con of the common constraint subset */

  int n; /* size of the common constraint subset */

#if MPI
  int rank [2]; /* rank[i] >= 0 if bd[i] corresponds to a remote dummy body */
#endif

  BDP *next;
};

/* create a list of body data items */
static BD* body_data_create (DOM *dom)
{
  BD *data, *bd;
  BODY *bod;
  SET *item;
  CON *con;
  int i;
  
  for (data = NULL, bod = dom->bod; bod; bod = bod->next)
  {
    if (bod->kind == OBS) continue;

    ERRMEM (bd = MEM_CALLOC (sizeof (BD)));
    bd->n = SET_Size (bod->con);
    ERRMEM (bd->con = MEM_CALLOC (bd->n * sizeof (CON*)));
    ERRMEM (bd->H = MEM_CALLOC (bd->n * sizeof (MX*)));
    ERRMEM (bd->b = malloc (2 * bod->dofs * sizeof (double)));
    bd->u = bd->b + bod->dofs;
    bd->next = data;
    data = bd;

    for (item = SET_First (bod->con), i = 0; item; item = SET_Next (item), i ++)
    {
      con = item->data;
      bd->con [i] = con;
      MAP_Insert (NULL, &bd->contoi, con, (void*) (long) i, NULL);
      if (bod == con->master) bd->H [i] = BODY_Gen_To_Loc_Operator (bod, con->mshp, con->mgobj, con->mpnt, con->base);
      else bd->H [i] = BODY_Gen_To_Loc_Operator (bod, con->sshp, con->sgobj, con->spnt, con->base);
    }
  }

  return data;
}

/* destroy a list of body data items */
static void body_data_destroy (BD *data)
{
  BD *next;
  int i;

  for (; data; data = next)
  {
    next = data->next;
    for (i = 0; i < data->n; i ++) MX_Destroy (data->H [i]);
    MAP_Free (NULL, &data->contoi);
    free (data->con);
    free (data->H);
    free (data->b);
    free (data);
  }
}

/* create a list of body data pairs */
static BDP* body_data_pairs_create (BD *data)
{
  MAP *bodtobd, *map;
  BDP *pairs, *dp;
  CON *con;
  MEM mem;
  BD *bd;
  int i;

  MEM_Init (&mem, sizeof (MAP), 128);

  for (bodtobd = NULL, bd = data; bd; bd = bd->next) MAP_Insert (&mem, &bodtobd, bd->bod, bd, NULL);

#if MPI
  /* TODO: import needed dummy body data, create BD and map it to 'bodtobd' */
#endif

  for (pairs = NULL, map = NULL, bd = data; bd; bd = bd->next)
  {
    for (i = 0; i < bd->n; i ++)
    {
      con = bd->con [i];

      if (con->slave && bd->bod == con->master) /* pair owner */
      {
	if (!(dp = MAP_Find (map, con->slave, NULL)))
	{
	  ERRMEM (dp = MEM_CALLOC (sizeof (BDP)));
	  dp->bd [0] = bd->bod;
	  ASSERT_DEBUG_EXT (dp->bd [1] = MAP_Find (bodtobd, con->slave, NULL), "Inconsistent body mapping");
	  dp->next = pairs;
	  pairs = dp;
#if MPI
	  /* TODO: set up rank */
#endif
	}

	ERRMEM (dp->con [0] = realloc (dp->con[0], (dp->n+1) * sizeof (int)));
	ERRMEM (dp->con [1] = realloc (dp->con[1], (dp->n+1) * sizeof (int)));
	dp->con [0][dp->n] = i;
	ASSERT_DEBUG (MAP_Find_Item (dp->bd [1]->contoi, con, NULL), "Inconsistent constraint mapping");
	dp->con [1][dp->n] = (int) (long) MAP_Find (dp->bd [1]->contoi, con, NULL);
	dp->n ++;
      }
    }
  }

  MEM_Release (&mem);
  
  return pairs;
}

/* destroy a list of body data pairs */
static void body_data_pairs_destroy (BDP *pairs)
{
  BDP *next;

  for (; pairs; pairs = next)
  {
    next = pairs->next;
    free (pairs->con [0]);
    free (pairs->con [1]);
#if MPI
    /* TODO: free imported dummy BD */
#endif
    free (pairs);
  }
}

/* create and return H and d such that H * u <= d describes the constraints tangent space;
 * initialize bd->b with free velocity summed up with previous tangential contact response */
static MX* velocity_projection_setup (BD *bd, double **d)
{
  int i, j, *n, *k, m;
  double *x, *y, *b;
  MX **G, *H;
  CON *con;

  ERRMEM (*d = malloc (6 * bd->n * sizeof (double)));
  ERRMEM (G = malloc (3 * bd->n * sizeof (MX*)));
  ERRMEM (n = malloc (4 * bd->n * sizeof (int)));
  ERRMEM (y = malloc (bd->n * sizeof (double)));
  k = n + 2 * bd->n;

  /* initialize bd->b with bd->bod->u */
  m = bd->bod->dofs; b = bd->b;
  blas_dcopy (m, bd->bod->u, 1, b, 1);

  /* calculate H composition data and set up d */
  for (i = j = 0, x = *d; i < bd->n; i ++)
  {
    con = bd->con [i];
    H = bd->H [i];

    switch (con->kind)
    {
    case CONTACT:
    {
      {
	/* TODO */
      }

      G [j] = H;
      n [j] = 1;
      k [j] = 2;
      y [j] = 1.0;
      j += 1;
      x += 1;
    }
    break;
    case FIXPNT:
    {
      if (con->slave)
      {
	/* TODO */
      }
      else
      {
	SET6 (x, 0.0);
      }

      G [j] = H;
      n [j] = 3;
      k [j] = 0;
      y [j] = 1.0;
      G [j+1] = H;
      n [j+1] = 3;
      k [j+1] = 0;
      y [j+1] = -1.0;
      j += 2;
      x += 6;
    }
    break;
    case FIXDIR:
    {
      {
	x [0] = 0.0;
	x [1] = 0.0;
      }

      G [j] = H;
      n [j] = 1;
      k [j] = 2;
      y [j] = 1.0;
      G [j+1] = H;
      n [j+1] = 1;
      k [j+1] = 2;
      y [j+1] = -1.0;
      j += 2;
      x += 2;
    }
    break;
    case VELODIR:
    {
      {
	x [0] = VELODIR(con->Z);
	x [1] = -VELODIR(con->Z);
      }

      G [j] = H;
      n [j] = 1;
      k [j] = 2;
      y [j] = 1.0;
      G [j+1] = H;
      n [j+1] = 1;
      k [j+1] = 2;
      y [j+1] = -1.0;
      j += 2;
      x += 2;
    }
    break;
    case RIGLNK:
    {
      {
	x [0] = 0.0;
	x [1] = 0.0;
      }

      G [j] = H;
      n [j] = 1;
      k [j] = 2;
      y [j] = 1.0;
      G [j+1] = H;
      n [j+1] = 1;
      k [j+1] = 2;
      y [j+1] = -1.0;
      j += 2;
      x += 2;
    }
    break;
    }

    /* TODO: sum up bd->b */
  }

  H = MX_Compose_From_Rows (j, G, n, k, y);

  free (G);
  free (n);
  free (y);

  return H;
}

/* approximate bd->u = proj (bd->b on T), where T: H * u <= d */
static void velocity_projection (BD *bd)
{
}

/* create and return [H_i'; H_j'] and [M_i(u_i - b_i); M_j(u_j - b_j)]' */
static MX* normal_reactions_setup (BDP *dp, double **c)
{
}

/* compute R0 such that [M_i(u_i - b_i); M_j(u_j - b_j)]' = [H_i'; H_j'] R0;
 * for normal contact reactions do RN_i = proj (R0_i on R+); otherwise R_j = R0_j */
static void normal_reactions (BDP *dp)
{
}

/* compute tangential contact response */
static void tangential_reactions (BD *bd)
{
}

/* compute single-body equality constraints */
static void equality_constraints (BD *bd)
{
}
#endif

/* create solver */
HYBRID* HYBRID_Create (double epsilon, int maxiter)
{
  HYBRID *hb;

  ERRMEM (hb = malloc (sizeof (HYBRID)));
  hb->epsilon = epsilon;
  hb->maxiter = maxiter;

  return hb;
}

static void gluing_solve (LOCDYN *ldy)
{
  double error, merit, step;
  int verbose, diagiters;
  int div = 10, iters;
  short dynamic;
  char fmt [512];

  verbose = ldy->dom->verbose;

  if (verbose) sprintf (fmt, "GLUING: iteration: %%%dd  error:  %%.2e  merit:  %%.2e\n", (int)log10 (1000) + 1);

  int maxiter = 100;
  double epsilon = 1E-3,
	 meritval = 1.0E6;

  GLUE *glu = GLUE_Create (ldy);

  dynamic = ldy->dom->dynamic;
  step = ldy->dom->step;
  iters = 0;
  do
  {
    double errup = 0.0,
	   errlo = 0.0;
    OFFB *blk;
    DIAB *dia;
   
    for (dia = ldy->dia; dia; dia = dia->n)
    {
      double R0 [3],
	     B [3],
	     *R = dia->R;

      CON *con = dia->con;

      if (con->kind == GLUEPNT) continue;

      /* compute local free velocity */
      COPY (dia->B, B);
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *W = blk->W,
	       *R = blk->dia->R;
	NVADDMUL (B, W, R, B);
      }
      
      COPY (R, R0); /* previous reaction */

      /* solve local diagonal block problem */
      diagiters = DIAGONAL_BLOCK_Solver (GS_PROJECTED_GRADIENT, 1E-6, 100,
	         dynamic, step, con->kind, con->mat.base, con->gap, con->Z, con->base, dia, B);

      /* accumulate relative
       * error components */
      SUB (R, R0, R0);
      errup += DOT (R0, R0);
      errlo += DOT (R, R);
    }

    /* merit function value */
    merit = MERIT_Function (ldy);

    /* calculate relative error */
    error = sqrt (errup) / sqrt (MAX (errlo, 1.0));

    GLUE_Solve (glu, error*1E-2, 1000);

    if (iters % div == 0 && verbose) printf (fmt, iters, error, merit), div *= 2;
  }
  while (++ iters < maxiter && (error > epsilon || merit > meritval));

  if (verbose) printf (fmt, iters, error, merit);

  GLUE_Destroy (glu);
}

/* run solver */
void HYBRID_Solve (HYBRID *hb, DOM *dom)
{
#if ENABLE
  BDP *pairs, *dp;
  BD *data, *bd;

  data = body_data_create (dom);

  pairs = body_data_pairs_create (data);

  for (bd = data; bd; bd = bd->next)
  {
    velocity_projection (bd);
  }

#if MPI
  /* TODO: send BD->u to remote BDP->bd[]->u */
#endif

  for (dp = pairs; dp; dp = dp->next)
  {
    normal_reactions (dp);
  }

#if MPI
  /* TODO: send normal reactions to remote bodies */
#endif

  for (bd = data; bd; bd = bd->next)
  {
    tangential_reactions (bd);
  }

  for (bd = data; bd; bd = bd->next)
  {
    equality_constraints (bd);
  }

  body_data_pairs_destroy (pairs);

  body_data_destroy (data);

  /* TODO: scale all reactions by 1/step */
#else
  gluing_solve (dom->ldy);
#endif
}

/* write labeled satate values */
void HYBRID_Write_State (HYBRID *hb, PBF *bf)
{
}

/* destroy solver */
void HYBRID_Destroy (HYBRID *hb)
{
  free (hb);
}
