/*
 * ldy.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * the local dynamic problem
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

#include "sol.h"
#include "alg.h"
#include "dom.h"
#include "ldy.h"
#include "lap.h"
#include "msh.h"
#include "err.h"

/* memory block size */
#define BLKSIZE 128

enum update_kind /* update kind */
{
  UPPES, /* penalty solver specific update */
  UPNOTHING, /* skip update */
  UPALL /* update all data */
};

typedef enum update_kind UPKIND;

/* get update kind depending on a solver */
static UPKIND update_kind (SOLVER_KIND solver)
{
  switch (solver)
  {
    case PENALTY_SOLVER: return UPPES;
    default: return UPALL;
  }

  return UPNOTHING;
}

/* apply forward change of variables (nornal
 * contact forces) due to the cohesion, etc. */
static void variables_change_begin (LOCDYN *ldy)
{
  OFFB *blk;
  DIAB *dia;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;

    if (con->kind != CONTACT) continue; /* skip non-contacts */
    else if (con->mat.base->model == SPRING_DASHPOT) continue; /* skip spring-dashpots */

    double *B = dia->B; /* free velocity will
			   be eventually modified */

    if (con->state & CON_COHESIVE) /* cohesive state */
    {
      double c = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area,
	     *W = dia->W,
	     *R = dia->R;

      R [2] += c;       /* R_n_new = R_n + c <=> (R_n + c) >= 0 */
      B[0] -= (W[6]*c); /* in consequnce 'W_tn * c' gets subtracted */
      B[1] -= (W[7]*c); /* ... */
      B[2] -= (W[8]*c); /* and 'W_nn * c' here */
    }

    /* off-diagonal subtractions */
    for (blk = dia->adj; blk; blk = blk->n)
    {
      CON *con = blk->dia->con;
      if (con->state & CON_COHESIVE) /* cohesive state */
      {
	double c = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area,
	       *W = blk->W;

	B[0] -= (W[6]*c);
	B[1] -= (W[7]*c);
	B[2] -= (W[8]*c);
      }
    }
  }
}

/* apply back change of variables (nornal
 * contact forces) due to the cohesion, etc. */
static void variables_change_end (LOCDYN *ldy)
{
  DIAB *dia;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;

    if (con->kind != CONTACT) continue; /* skip non-contacts */
    else if (con->mat.base->model == SPRING_DASHPOT) continue; /* skip spring-dashpots */

    short state = con->state;

    if (state & CON_COHESIVE) /* cohesive state */
    {
      double c = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area,
	     f = con->mat.base->friction,
	     e = 1E-8,
	     *R = dia->R;

      if (R [2] < e || /* mode-I decohesion */
	  LEN2 (R) + e >= f * R[2]) /* mode-II decohesion */
      {
	con->state &= ~CON_COHESIVE;
	SURFACE_MATERIAL_Cohesion_Set (&con->mat, 0.0);
      }

      R [2] -= c; /* back change */
    }
  }
}

/* test whether two constraints are able to be adjacent */
static int adjacentable (BODY *bod, CON *one, CON *two)
{
  if (bod->kind == FEM)
  {
    if (bod->msh)
    {
      ELEMENT **e1, **f1, **e2, **f2;
      CONVEX *c1 = (bod == one->master ? mgobj(one) : sgobj(one)),
	     *c2 = (bod == two->master ? mgobj(two) : sgobj(two));

      for (e1 = c1->ele, f1 = e1 + c1->nele; e1 < f1; e1 ++)
      {
	for (e2 = c2->ele, f2 = e2 + c2->nele; e2 < f2; e2 ++)
	{
	  if (*e1 == *e2) return 1;
	  else if (ELEMENT_Adjacent (*e1, *e2)) return 1;
	}
      }

      return 0;
    }
    else
    {
      ELEMENT *e1 = (bod == one->master ? mgobj(one) : sgobj(one)),
	      *e2 = (bod == two->master ? mgobj(two) : sgobj(two));

      return ELEMENT_Adjacent (e1, e2); /* only in case of a common node W_one_two and W_two_one will be != 0 */
    }
  }

  return 1;
}

#if MPI
/* compute external adjacency */
static void compute_adjext (LOCDYN *ldy, UPKIND upkind)
{
  CON *con, *ext;
  OFFB *b, *n;
  BODY *bod;
  SET *item;
  MAP *jtem;
  DIAB *dia;
  int i;

  /* clear previous external adjacency */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    for (b = dia->adjext; b; b = n)
    {
      n = b->n;
      MEM_Free (&ldy->offmem, b);
    }

    dia->adjext = NULL;
  }

  /* walk over all external contacts and build new external adjacency */
  for (jtem = MAP_First (ldy->dom->conext); jtem; jtem = MAP_Next (jtem))
  {
    ext = jtem->data;

    BODY *bodies [] = {ext->master, ext->slave}; /* e.g. two remote child bodies of local parents might be in contact */

    for (i = 0, bod = bodies [i]; i < 2 && bod; i ++, bod = bodies [i]) /* (i < 2 && bod) skips NULL slaves of single-body constraints */
    {
      if (bod->kind == OBS) continue; /* obstacles do not trasnder adjacency */

      for (item = SET_First (bod->con); item; item = SET_Next (item)) 
      {
	con = item->data;

	if (con->state & CON_EXTERNAL) continue; /* for each regular constraint */

        if (upkind == UPPES && con->kind == CONTACT) continue; /* skip contacts during partial update (pes.c uses local dynamics only for non-contacts) */

	ASSERT_DEBUG (bod->flags & (BODY_PARENT|BODY_CHILD), "Regular constraint attached to a dummy"); /* we could skip dummies, but this reassures correctness */

	ASSERT_DEBUG (con->master == bod || con->slave == bod, "Incorrectly connected constraint in a body constraints list: "
	              "master->id = %d, slave->id = %d, bod->id = %d", con->master->id, con->slave->id, bod->id);

	dia = con->dia;
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = (DIAB*) ext; /* there is no diagonal block here, but we shall point directly to the external contact */
	b->bod = bod; /* adjacent through this body */
	b->n = dia->adjext;
	dia->adjext = b;
      }
    }
  }
}

#if PARDEBUG
/* return next pointer and realloc send memory if needed */
inline static COMDATA* sendnext (int nsend, int *size, COMDATA **send)
{
  if (nsend >= *size)
  {
    (*size) *= 2;
    ERRMEM (*send = realloc (*send, sizeof (COMDATA [*size])));
  }

  return &(*send)[nsend];
}

/* test consistency of external adjacency */
static int adjext_test (LOCDYN *ldy)
{
  int ssend, nsend, nrecv, i, j, *k, ret;
  COMDATA *ptr, *send, *recv;
  SET *ranks, *item;
  MEM setmem;
  DIAB *dia;
  OFFB *blk;
  CON *con;

  MEM_Init (&setmem, sizeof (SET), 128);

  ssend = 128;
  nsend = 0;
  ERRMEM (send = MEM_CALLOC (ssend * sizeof (COMDATA)));
  ptr = send;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    for (ranks = NULL, blk = dia->adjext; blk; blk = blk->n)
    { 
      con = (CON*) blk->dia;
      SET_Insert (&setmem, &ranks, (void*) (long) con->rank, NULL);
      con->state &= ~CON_DONE;
    }

    for (item = SET_First (ranks); item; item = SET_Next (item))
    {
      ptr->rank = (int) (long) item->data;
      ptr->ints = 1;
      ptr->doubles = 0;
      ptr->i = (int*) &dia->con->id;
      ptr->d = NULL;
      ptr= sendnext (++ nsend, &ssend, &send);
    }

    SET_Free (&setmem, &ranks);
  }

  COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);

  for (ret = 1, i = 0, ptr = recv; i < nrecv; i ++, ptr ++)
  {
    for (j = 0, k = ptr->i; j < ptr->ints; j ++, k ++)
    {
      if (!(con = MAP_Find (ldy->dom->conext, (void*) (long) (*k), NULL)))
      {
	ret = 0;
	WARNING_DEBUG (0, "External donstraint %d from rank %d not FOUND on rank %d", (*k), ptr->rank, ldy->dom->rank);
	goto out;
      }
      con->state |= CON_DONE;
    }
  }

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = (CON*) blk->dia;
      if ((con->state & CON_DONE) == 0)
      {
	ret = 0;
	WARNING_DEBUG (0, "External constraint %u not UPDATED on rank %d", con->id, ldy->dom->rank);
	goto out;
      }
    }
  }

out:
  MEM_Release (&setmem);
  free (send);
  free (recv);
  return ret;
}
#endif
#endif

/* dump comparison */
static int dumpcmp (CON *a, CON *b)
{
  for (int i = 0; i < 3; i ++)
  {
    if (a->point [i] < b->point [i]) return -1;
    else if (a->point [i] > b->point [i]) return 1;
  }

  return a < b ? -1 : a > b ? 1 : 0;
}

/* create local dynamics for a domain */
LOCDYN* LOCDYN_Create (DOM *dom)
{
  LOCDYN *ldy;

  ERRMEM (ldy = malloc (sizeof (LOCDYN)));
  MEM_Init (&ldy->offmem, sizeof (OFFB), BLKSIZE);
  MEM_Init (&ldy->diamem, sizeof (DIAB), BLKSIZE);
  ldy->dom = dom;
  ldy->dia = NULL;
  ldy->modified = 0;

  return ldy;
}

/* insert a 'con'straint between a pair of bodies =>
 * return the diagonal entry of the local dynamical system */
DIAB* LOCDYN_Insert (LOCDYN *ldy, CON *con, BODY *one, BODY *two)
{
  DIAB *dia, *nei;
  double *SYMW;
  SET *item;
  OFFB *b;
  CON *c;

  ERRMEM (dia = MEM_Alloc (&ldy->diamem));
  dia->R = con->R;
  dia->U = con->U;
  dia->con = con;

  /* insert into list */
  dia->n = ldy->dia;
  if (ldy->dia)
    ldy->dia->p = dia;
  ldy->dia = dia;

  if (one && one->kind != OBS) /* obstacles do not transfer adjacency */
  {
    for (item = SET_First (one->con); item; item = SET_Next (item))
    {
      c = item->data;

      if (c != con && c->dia && /* skip the coincident or unattached yet constraint */
	  adjacentable (one, con, c)) /* skip other cases where W_ij would be zero */
      {
	nei = c->dia;

        /* allocate block and put into 'nei->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = dia; /* adjacent with 'dia' */
	b->bod = one; /* adjacent trough body 'one' */
	SYMW = b->W, b->SYMW = NULL; /* compute when assembling */
	b->n = nei->adj; /* extend list ... */
	nei->adj = b; /* ... */

	/* allocate block and put into 'dia->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = nei; /* adjacent with 'nei' */
	b->bod = one; /* ... trough 'one' */
	b->SYMW = SYMW; /* compy when assembling */
	b->n = dia->adj;
	dia->adj = b;
      }
    }
  }

  if (two && two->kind != OBS) /* 'one' replaced with 'two' */
  {
    for (item = SET_First (two->con); item; item = SET_Next (item))
    {
      c = item->data;

      if (c != con && c->dia && /* skip the coincident or unattached yet constraint */
	  adjacentable (two, con, c)) /* skip other cases where W_ij would be zero */
      {
	nei = c->dia;

        /* allocate block and put into 'nei->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = dia; /* adjacent with 'dia' */
	b->bod = two; /* adjacent trough body 'two' */
	SYMW = b->W, b->SYMW = NULL; /* compute when assembling */
	b->n = nei->adj; /* extend list ... */
	nei->adj = b; /* ... */

	/* allocate block and put into 'dia->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = nei; /* adjacent with 'nei' */
	b->bod = two; /* ... trough 'two' */
	b->SYMW = SYMW; /* compy when assembling */
	b->n = dia->adj;
	dia->adj = b;
      }
    }
  }

  /* mark as modified */
  ldy->modified = 1;

  return dia;
}

/* remove a diagonal entry from local dynamics */
void LOCDYN_Remove (LOCDYN *ldy, DIAB *dia)
{
  OFFB *b, *c, *r;

  /* destroy blocks in
   * adjacent dia items */
  for (b = dia->adj; b; b = b->n)
  {
    c = b->dia->adj;
    if (c && c->dia == dia)
    {
      b->dia->adj = c->n;
      MEM_Free (&ldy->offmem, c); 
    }
    else for (; c; c = c->n)
    {
      if (c->n &&
	  c->n->dia == dia)
      {
	r = c->n;
        c->n = c->n->n;
        MEM_Free (&ldy->offmem, r); 
	break;
      }
    }
  }

  /* destroy directly
   * adjacent blocks */
  for (b = dia->adj; b; b = c)
  {
    c = b->n;
    MEM_Free (&ldy->offmem, b);
  }

#if MPI
  /* destroy externally
   * adjacent blocks */
  for (b = dia->adjext; b; b = c)
  {
    c = b->n;
    MEM_Free (&ldy->offmem, b);
  }
#endif

  /* remove from list */
  if (dia->p)
    dia->p->n = dia->n;
  else ldy->dia = dia->n;
  if (dia->n)
    dia->n->p = dia->p;

  /* destroy passed dia */
  MEM_Free (&ldy->diamem, dia);

  /* mark as modified */
  ldy->modified = 1;
}

/* updiae local dynamics => prepare for a solution */
void LOCDYN_Update_Begin (LOCDYN *ldy, SOLVER_KIND solver)
{
  UPKIND upkind = update_kind (solver);
  DOM *dom = ldy->dom;
  double step = dom->step;
  DIAB *dia;
  OFFB *blk;

  if (upkind == UPNOTHING) return; /* skip update */

#if MPI
  if (dom->rank == 0)
#endif
  if (dom->verbose) printf ("LOCDYN ... "), fflush (stdout);

  SOLFEC_Timer_Start (ldy->dom->solfec, "LOCDYN");

#if MPI
  compute_adjext (ldy, upkind);
#endif

  ldy->free_energy = 0.0;

  /* calculate local velocities and
   * assmeble the force-velocity 'W' operator */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;
    BODY *m = con->master,
	 *s = con->slave;
    void *mgobj = mgobj(con),
	 *sgobj;
    SHAPE *mshp = mshp(con),
	  *sshp;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base,
	   *V = dia->V,
	   *B = dia->B,
	   X0 [3], Y0 [3],
           X [3], Y [9];
    MX_DENSE_PTR (W, 3, 3, dia->W);
    MX_DENSE_PTR (A, 3, 3, dia->A);
    MX_DENSE (C, 3, 3);

    /* relative velocity = slave - master => outward master normal */
    BODY_Local_Velo (m, mshp, mgobj, mpnt, base, X0, X); /* master body pointer cannot be NULL */
    if (s)
    {
      sgobj = sgobj(con);
      sshp = sshp(con);
      BODY_Local_Velo (s, sshp, sgobj, spnt, base, Y0, Y); /* might be NULL for some constraints (one body) */
    }
    else { SET (Y0, 0.0); SET (Y, 0.0); }

    SUB (X0, Y0, V); /* previous time step velocity */
    SUB (X, Y, B); /* local free velocity */

    /* diagonal block */
    dia->mH = BODY_Gen_To_Loc_Operator (m, mshp, mgobj, mpnt, base);
#if MPI
    dia->mprod = MX_Matmat (1.0, dia->mH, m->inverse, 0.0, NULL);
    MX_Matmat (1.0, dia->mprod, MX_Tran (dia->mH), 0.0, &W); /* H * inv (M) * H^T */
#else
    dia->mprod = MX_Matmat (1.0, m->inverse, MX_Tran (dia->mH), 0.0, NULL);
    MX_Matmat (1.0, dia->mH, dia->mprod, 0.0, &W); /* H * inv (M) * H^T */
#endif
    if (s)
    { 
      dia->sH = BODY_Gen_To_Loc_Operator (s, sshp, sgobj, spnt, base);
#if MPI
      dia->sprod = MX_Matmat (1.0, dia->sH, s->inverse, 0.0, NULL);
      MX_Matmat (1.0, dia->sprod, MX_Tran (dia->sH), 0.0, &C); /* H * inv (M) * H^T */
#else
      dia->sprod = MX_Matmat (1.0, s->inverse, MX_Tran (dia->sH), 0.0, NULL);
      MX_Matmat (1.0, dia->sH, dia->sprod, 0.0, &C); /* H * inv (M) * H^T */
#endif
      NNADD (W.x, C.x, W.x);
    }
    SCALE9 (W.x, step); /* W = h * ( ... ) */

    if (upkind != UPPES) /* diagonal regularization (not needed by the explicit solver) */
    {
      NNCOPY (W.x, C.x); /* calculate regularisation parameter */
      ASSERT (lapack_dsyev ('N', 'U', 3, C.x, 3, X, Y, 9) == 0, ERR_LDY_EIGEN_DECOMP);
      dia->rho = 1.0 / X [2]; /* inverse of maximal eigenvalue */
    }

    NNCOPY (W.x, A.x);
    MX_Inverse (&A, &A); /* inverse of diagonal block */

    NVMUL (A.x, B, X);
    ldy->free_energy += DOT (X, B); /* sum up free energy */

    /* add up prescribed velocity contribution */
    if (con->kind == VELODIR) ldy->free_energy += A.x[8] * VELODIR(con->Z) * VELODIR(con->Z);
  }

  ldy->free_energy *= 0.5; /* 0.5 * DOT (AB, B) */

  for (dia = ldy->dia; dia; dia = dia->n) /* off-diagonal blocks update */
  {
    CON *con = dia->con;
    BODY *m = con->master,
	 *s = con->slave;

    if (upkind == UPPES && con->kind == CONTACT) continue; /* update only non-contact constraint blocks */

    /* off-diagonal local blocks */
    for (blk = dia->adj; blk; blk = blk->n)
    {
      if (upkind == UPALL && blk->SYMW) continue; /* skip blocks pointing to their symmetric copies (only during full updates) */

      MX *left, *right;
      DIAB *adj = blk->dia;
      CON *con = adj->con;
      BODY *bod = blk->bod;
      MX_DENSE_PTR (W, 3, 3, blk->W);
      double coef;

      ASSERT_DEBUG (bod == m || bod == s, "Off diagonal block is not connected!");
     
#if MPI
      left = (bod == m ? dia->mprod : dia->sprod);
#else
      left = (bod == m ? dia->mH : dia->sH);
#endif

      if (bod == con->master) /* master on the right */
      {
#if MPI
	right = adj->mH;
#else
	right =  adj->mprod;
#endif
	coef = (bod == s ? -step : step);
      }
      else /* blk->bod == con->slave (slave on the right) */
      {
#if MPI
	right = adj->sH;
#else
	right =  adj->sprod;
#endif
	coef = (bod == m ? -step : step);
      }

#if MPI
      MX_Matmat (1.0, left, MX_Tran (right), 0.0, &W);
#else
      MX_Matmat (1.0, left, right, 0.0, &W);
#endif
      SCALE9 (W.x, coef);
    }

#if MPI
    /* off-diagonal external blocks */
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      MX *left, *right;
      CON *ext = (CON*)blk->dia;
      BODY *bod = blk->bod;
      MX_DENSE_PTR (W, 3, 3, blk->W);
      double coef, *point;
      SGP *sgp;

      ASSERT_DEBUG (bod == m || bod == s, "Not connected external off-diagonal block");

      if (bod == ext->master)
      {
	sgp = ext->msgp, point = ext->mpnt;
	coef = (bod == s ? -step : step);
      }
      else
      {
	sgp = ext->ssgp, point = ext->spnt;
	coef = (bod == m ? -step : step);
      }
     
      left = (bod == m ? dia->mprod : dia->sprod);

      right =  BODY_Gen_To_Loc_Operator (bod, sgp->shp, sgp->gobj, point, ext->base);

      MX_Matmat (1.0, left, MX_Tran (right), 0.0, &W);
      SCALE9 (W.x, coef);
      MX_Destroy (right);
    }
#endif
  }

  /* use symmetry */
  if (upkind == UPALL)
  {
    for (dia = ldy->dia; dia; dia = dia->n)
    {  
      for (blk = dia->adj; blk; blk = blk->n)
      {
	if (blk->SYMW)
	{
	  TNCOPY (blk->SYMW, blk->W); /* transposed copy of a symmetric block */
	}
      }
    }
  }

  /* clean up */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    MX_Destroy (dia->mH);
    MX_Destroy (dia->mprod);
    if (dia->sH)
    {
      MX_Destroy (dia->sH);
      MX_Destroy (dia->sprod);
    }
  }

#if PARDEBUG
  if (upkind == UPALL)
  {
    if (adjext_test (ldy) == 0)
    {
      ASSERT_DEBUG (0, "Inconsistent adjext"); /* a debugger catchable assertion */
    }
  }
#endif

  /* forward variables change */
  if (upkind == UPALL) variables_change_begin (ldy);

  SOLFEC_Timer_End (ldy->dom->solfec, "LOCDYN");
}

/* updiae local dynamics => after the solution */
void LOCDYN_Update_End (LOCDYN *ldy, SOLVER_KIND solver)
{
  UPKIND upkind = update_kind (solver);

  if (upkind == UPNOTHING) return; /* skip update */

  SOLFEC_Timer_Start (ldy->dom->solfec, "LOCDYN");

  /* backward variables change */
  if (upkind == UPALL) variables_change_end (ldy);

  /* not modified */
  ldy->modified = 0;

  SOLFEC_Timer_End (ldy->dom->solfec, "LOCDYN");
}

/* dump local dynamics to file */
void LOCDYN_Dump (LOCDYN *ldy, const char *path)
{
  MEM mapmem, offmem;
  MAP *adj, *item;
  char *fullpath;
  OFFB *blk, *q;
  DIAB *dia;
  CON *con;
  FILE *f;

#if MPI
  ERRMEM (fullpath = malloc (strlen (path) + 64));
  snprintf (fullpath, strlen (path) + 64, "%s.%d", path, ldy->dom->rank);
#else
  fullpath = (char*) path;
#endif

  ASSERT (f = fopen (fullpath, "w"), ERR_FILE_OPEN);

  MEM_Init (&offmem, sizeof (OFFB), BLKSIZE);
  MEM_Init (&mapmem, sizeof (MAP), BLKSIZE);

  adj = NULL;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    con = dia->con;

    fprintf (f, "%s (%.6g, %.6g, %.6g) [%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g] (%u, %d, %u, %d) [%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g] => ",
      CON_Kind (con), con->point [0], con->point [1], con->point [2],
      con->base [0], con->base [1], con->base [2], con->base [3], con->base [4], con->base [5], con->base [6], con->base [7], con->base [8],
      con->master->id, con->msgp - con->master->sgp, con->slave ? con->slave->id : 0, con->slave ? con->ssgp - con->slave->sgp : 0,
      dia->W [0], dia->W [1], dia->W [2], dia->W [3], dia->W [4], dia->W [5], dia->W [6], dia->W [7], dia->W [8]);

    MAP_Free (&mapmem, &adj);
    MEM_Release (&offmem);

    for (blk = dia->adj; blk; blk = blk->n)
    {
      if (!(q = MAP_Find (adj, blk->dia->con, (MAP_Compare) dumpcmp)))
      {
	ERRMEM (q = MEM_Alloc (&offmem));
	ERRMEM (MAP_Insert (&mapmem, &adj, blk->dia->con, q, (MAP_Compare) dumpcmp));
      }
      NNADD (q->W, blk->W, q->W);
    }

#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      if (!(q = MAP_Find (adj, blk->dia, (MAP_Compare) dumpcmp)))
      {
	ERRMEM (q = MEM_Alloc (&offmem));
	ERRMEM (MAP_Insert (&mapmem, &adj, blk->dia, q, (MAP_Compare) dumpcmp));
      }
      NNADD (q->W, blk->W, q->W);
    }
#endif

    for (item = MAP_First (adj); item; item = MAP_Next (item))
    {
      con = item->key;
      q = item->data;

      fprintf (f, "%s (%.6g, %.6g, %.6g) [%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g] ",
        CON_Kind (con), con->point [0], con->point [1], con->point [2],
	q->W [0], q->W [1], q->W [2], q->W [3], q->W [4], q->W [5], q->W [6], q->W [7], q->W [8]);

    }

    fprintf (f, "\n");
  }

  MEM_Release (&mapmem);
  MEM_Release (&offmem);
  fclose (f);

#if MPI
  MPI_Barrier (MPI_COMM_WORLD);

  if (ldy->dom->rank == 0)
  {
    ASSERT (f = fopen (path, "w"), ERR_FILE_OPEN);
    for (int i = 0; i < ldy->dom->ncpu; i ++)
    {
      char *buf;
      long len;
      FILE *g;

      snprintf (fullpath, strlen (path) + 64, "%s.%d", path, i);
      ASSERT (g = fopen (fullpath, "r"), ERR_FILE_OPEN);
      fseek (g, 0, SEEK_END);
      len = ftell (g);
      ERRMEM (buf = malloc (len + 64));
      fseek (g, 0, SEEK_SET);
      fread (buf, 1, len, g);
      fwrite (buf, 1, len, f);
      fclose (g);
      free (buf);
      remove (fullpath);
    }
    fclose (f);
  }

  free (fullpath);
#endif
}

/* free memory */
void LOCDYN_Destroy (LOCDYN *ldy)
{
  MEM_Release (&ldy->diamem);
  MEM_Release (&ldy->offmem);

  free (ldy);
}
