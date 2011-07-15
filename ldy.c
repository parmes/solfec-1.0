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

#if MPI
#include "com.h"
#include "pck.h"
#include "put.h"
#endif

/* memory block size */
#define BLKSIZE 128

enum update_kind /* update kind */
{
  UPPES, /* penalty solver update */
  UPMIN, /* minimal update */
  UPALL /* update all data */
};

typedef enum update_kind UPKIND;

/* get update kind depending on a solver */
static UPKIND update_kind (SOLFEC *sol)
{
  switch (sol->kind)
  {
    case PENALTY_SOLVER: return UPPES;
    case NEWTON_SOLVER:
    {
      NEWTON *ns = sol->solver;
      if (ns->locdyn == LOCDYN_OFF) return UPMIN;
      else return UPALL;
    }
    default: return UPALL;
  }

  return UPALL;
}

/* sort out cohesion states */
static void update_cohesion (LOCDYN *ldy)
{
  DIAB *dia;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;

    if (con->kind != CONTACT) continue; /* skip non-contacts */
    else if (con->mat.base->model == SPRING_DASHPOT) continue; /* skip spring-dashpots */

    if (con->state & CON_COHESIVE) /* cohesive state */
    {
      double c = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area,
	     f = con->mat.base->friction,
	     e = COHESION_EPSILON * c,
	    *R = con->R;

      if ((R [2]+c) < e || /* mode-I decohesion */
	  LEN2 (R) + e >= f * (R[2]+c)) /* mode-II decohesion */
      {
	con->state &= ~CON_COHESIVE;
	SURFACE_MATERIAL_Cohesion_Set (&con->mat, 0.0);
      }
    }
  }
}

/* test whether two constraints are able to be adjacent */
static int adjacentable (BODY *bod, CON *one, CON *two)
{
  if (bod->kind == FEM && bod->scheme == SCH_DEF_EXP)
  {
    /* XXX: for diagonal inverse matrix (explicit integration) only
     * XXX: in case of a common node W_one_two and W_two_one will be != 0 */

    if (bod->msh) /* rough mesh */
    {
      ELEMENT **e1, **f1, **e2, **f2;
      CONVEX *c1 = (bod == one->master ? SGP_2_GOBJ (one->msgp) : SGP_2_GOBJ (one->ssgp)),
	     *c2 = (bod == two->master ? SGP_2_GOBJ (two->msgp) : SGP_2_GOBJ (two->ssgp));

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
    else /* regular mesh */
    {
      MESH *m1, *m2;
      ELEMENT *e1, *e2;
      double *p1, *p2;
      int n1, n2;

      if (bod == one->master)
      {
	m1 = mshp (one)->data;
	e1 = SGP_2_GOBJ (one->msgp);
	p1 = one->mpnt;
      }
      else
      {
	m1 = sshp (one)->data;
	e1 = SGP_2_GOBJ (one->ssgp);
	p1 = one->spnt;
      }

      if (bod == two->master)
      {
	m2 = mshp (two)->data;
	e2 = SGP_2_GOBJ (two->msgp);
	p2 = two->mpnt;
      }
      else
      {
	m2 = sshp (two)->data;
	e2 = SGP_2_GOBJ (two->ssgp);
	p2 = two->spnt;
      }

      n1 = ELEMENT_Ref_Point_To_Node (m1, e1, p1);
      n2 = ELEMENT_Ref_Point_To_Node (m2, e2, p2);

      if (n1 >= 0 && n2 >= 0 && n1 != n2) return 0; /* distinct mesh nodes */
      else return ELEMENT_Adjacent (e1, e2); /* (non)adjacent elements */
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

	if (adjacentable (bod, ext, con)) /* if constraints interact insert W block */
	{
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

  int _aid [2] = {(int)a->master->id, a->slave ? (int)a->slave->id : -1},
      _bid [2] = {(int)b->master->id, b->slave ? (int)b->slave->id : -1},
       aid [2] = {MIN (_aid[0], _aid[1]), MAX (_aid[0], _aid[1])},
       bid [2] = {MIN (_bid[0], _bid[1]), MAX (_bid[0], _bid[1])};

  for (int i = 0; i < 2; i ++)
  {
    if (aid [i] < bid [i]) return -1;
    else if (aid [i] > bid [i]) return 1;
  }

  if (a->kind < b->kind) return -1;
  else if (a->kind > b->kind) return 1;

  ASSERT_DEBUG (a == b, "Two different constraints between same pair of bodies have same spatial point");

  return 0;
}

/* update previous and free local velocities */
static void update_V_and_B (DOM *dom)
{
  double X [6], *V, *B;
  CON *con;

  for (con = dom->con; con; con = con->next)
  {
    V = con->V;
    B = con->dia->B;

    SET (V, 0);
    SET (B, 0);

#if MPI
    if (con->master->flags & BODY_PARENT) /* local parent */
#endif
    {
      BODY_Local_Velo (con->master, con->msgp, con->mpnt, con->base, X, X+3);
      ADD (V, X, V);
      ADD (B, X+3, B);
    }

    if (con->slave)
    {
#if MPI
      if (con->slave->flags & BODY_PARENT) /* local slave */
#endif
      {
        BODY_Local_Velo (con->slave, con->ssgp, con->spnt, con->base, X, X+3);
	SUB (V, X, V); /* relative = master - slave */
	SUB (B, X+3, B);
      }
    }
  }

#if MPI
  /* in parallel, parent bodies calculate local velocities of their external constraints and send them to the constraint parents;
   * upon arrival these are added or subtracted from the constraint parents B and V members, depending on the master/slave relation */
  int nsend, nrecv, *isize, *dsize, i, *j, *k;
  COMDATA *send, *recv, *ptr;
  double *D;
  BODY *bod;
  SET *item;

  nsend = dom->ncpu;
  ERRMEM (send = MEM_CALLOC (sizeof (COMDATA [nsend])));
  ERRMEM (isize = MEM_CALLOC (sizeof (int [nsend])));
  ERRMEM (dsize = MEM_CALLOC (sizeof (int [nsend])));

  for (i = 0; i < nsend; i ++) send [i].rank = i;

  for (bod = dom->bod; bod; bod = bod->next) /* for all parents */
  {
    for (item = SET_First (bod->con); item; item = SET_Next (item))
    {
      con = item->data;
      if (con->state & CON_EXTERNAL) /* needs local velocity update */
      {
	ASSERT_DEBUG (MAP_Find (dom->conext, (void*) (long) con->id, NULL), "Invalid external constraint %s %d", CON_Kind (con), con->id);
	i = con->rank;
	ptr = &send [i];
	if (bod == con->master)
	{
	  pack_int (&isize [i], &ptr->i, &ptr->ints, -con->id);
	  BODY_Local_Velo (bod, con->msgp, con->mpnt, con->base, X, X+3);
	}
	else
	{
	  pack_int (&isize [i], &ptr->i, &ptr->ints, con->id);
	  BODY_Local_Velo (bod, con->ssgp, con->spnt, con->base, X, X+3);
	}
        pack_doubles (&dsize [i], &ptr->d, &ptr->doubles, X, 6);
      }
    }
  }

  COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv); /* send V, B */

  for (i = 0; i < nrecv; i ++)
  {
    ptr = &recv [i];
    for (j = ptr->i, k = j + ptr->ints, D = ptr->d; j < k; j ++, D += 6)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*) (long) ABS(*j), NULL), "Invalid con id: %d", ABS(*j));
      V = con->V;
      B = con->dia->B;
      if ((*j) < 0) /* master */
      {
	ADD (V, D, V);
	ADD (B, D+3, B);
      }
      else /* slave */
      {
	SUB (V, D, V);
	SUB (B, D+3, B);
      }
    }
  }

  for (i = 0; i < nsend; i ++)
  {
    ptr = &send [i];
    free (ptr->i);
    free (ptr->d);
  }
  free (send);
  free (isize);
  free (dsize);
  free (recv); /* includes recv[]->i and recv[]->d memory */
#endif
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

  return ldy;
}

/* insert a 'con'straint between a pair of bodies =>
 * return the diagonal entry of the local dynamical system */
DIAB* LOCDYN_Insert (LOCDYN *ldy, CON *con, BODY *one, BODY *two)
{
  DIAB *dia, *nei;
  SET *item;
  OFFB *b;
  CON *c;

  ERRMEM (dia = MEM_Alloc (&ldy->diamem));
  dia->R = con->R;
  dia->U = con->U;
  dia->V = con->V;
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
	b->n = nei->adj; /* extend list ... */
	nei->adj = b; /* ... */

	/* allocate block and put into 'dia->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = nei; /* adjacent with 'nei' */
	b->bod = one; /* ... trough 'one' */
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
	b->n = nei->adj; /* extend list ... */
	nei->adj = b; /* ... */

	/* allocate block and put into 'dia->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = nei; /* adjacent with 'nei' */
	b->bod = two; /* ... trough 'two' */
	b->n = dia->adj;
	dia->adj = b;
      }
    }
  }

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
      if (c->n && c->n->dia == dia)
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
}

void LOCDYN_Update_Begin (LOCDYN *ldy)
{
  DOM *dom = ldy->dom;
  UPKIND upkind = update_kind (dom->solfec);
  double step = dom->step;
  OFFB *blk, *blj;
  DIAB *dia;

#if MPI
  if (dom->rank == 0)
#endif
  if (dom->verbose) printf ("LOCDYN ... "), fflush (stdout);

  SOLFEC_Timer_Start (ldy->dom->solfec, "LOCDYN");

  /* update previous and free velocites */
  update_V_and_B (dom);

  if (upkind == UPMIN) goto end; /* skip update */

#if MPI
  compute_adjext (ldy, upkind);
#endif

  ldy->free_energy = 0.0;

  /* calculate local velocities and assmeble
   * the diagonal force-velocity 'W' operator */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;
    BODY *m = con->master,
	 *s = con->slave;
    SGP *msgp = con->msgp,
	*ssgp = con->ssgp;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base,
	   *B = dia->B,
           X [3], Y [9];
    MX_DENSE_PTR (W, 3, 3, dia->W);
    MX_DENSE_PTR (A, 3, 3, dia->A);
    MX_DENSE (C, 3, 3);

    if (con->kind == GLUE && dia->adj == NULL && W.x [8] != 0.0)
    {
#if MPI
      if (dia->adjext == NULL)
#endif
      goto sumene; /* skip initialized explicit node-to-node gluing constraints */
    }

#if MPI
    if (m->flags & BODY_CHILD)
    {
      if (dom->dynamic) BODY_Dynamic_Init (m);
      else BODY_Static_Init (m);
    }
#endif

    if (s)
    {
#if MPI
      if (s->flags & BODY_CHILD)
      {
	if (dom->dynamic) BODY_Dynamic_Init (s);
	else BODY_Static_Init (s);
      }
#endif
    }

    /* diagonal block */
    if (m != s)
    {
      dia->mH = BODY_Gen_To_Loc_Operator (m, msgp, mpnt, base);
#if MPI
      dia->mprod = MX_Matmat (1.0, dia->mH, m->inverse, 0.0, NULL);
      MX_Matmat (1.0, dia->mprod, MX_Tran (dia->mH), 0.0, &W); /* H * inv (M) * H^T */
#else
      dia->mprod = MX_Matmat (1.0, m->inverse, MX_Tran (dia->mH), 0.0, NULL);
      MX_Matmat (1.0, dia->mH, dia->mprod, 0.0, &W); /* H * inv (M) * H^T */
#endif

      if (s)
      {
	dia->sH = BODY_Gen_To_Loc_Operator (s, ssgp, spnt, base);
	MX_Scale (dia->sH, -1.0);
#if MPI
	dia->sprod = MX_Matmat (1.0, dia->sH, s->inverse, 0.0, NULL);
	MX_Matmat (1.0, dia->sprod, MX_Tran (dia->sH), 0.0, &C); /* H * inv (M) * H^T */
#else
	dia->sprod = MX_Matmat (1.0, s->inverse, MX_Tran (dia->sH), 0.0, NULL);
	MX_Matmat (1.0, dia->sH, dia->sprod, 0.0, &C); /* H * inv (M) * H^T */
#endif
	NNADD (W.x, C.x, W.x);
      }
    }
    else /* eg. self-contact */
    {
      MX *mH = BODY_Gen_To_Loc_Operator (m, msgp, mpnt, base),
	 *sH = BODY_Gen_To_Loc_Operator (s, ssgp, spnt, base);

      dia->mH = MX_Add (1.0, mH, -1.0, sH, NULL);
      dia->sH = MX_Copy (dia->mH, NULL);

      MX_Destroy (mH);
      MX_Destroy (sH);
#if MPI
      dia->mprod = MX_Matmat (1.0, dia->mH, m->inverse, 0.0, NULL);
      dia->sprod = MX_Copy (dia->mprod, NULL);
      MX_Matmat (1.0, dia->mprod, MX_Tran (dia->mH), 0.0, &W); /* H * inv (M) * H^T */
#else
      dia->mprod = MX_Matmat (1.0, m->inverse, MX_Tran (dia->mH), 0.0, NULL);
      dia->sprod = MX_Copy (dia->mprod, NULL);
      MX_Matmat (1.0, dia->mH, dia->mprod, 0.0, &W); /* H * inv (M) * H^T */
#endif
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

sumene: 
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
      if (upkind == UPALL && blk->dia < dia) continue; /* skip lower triangle */

      MX *left, *right;
      DIAB *adj = blk->dia;
      BODY *bod = blk->bod;
      CON *con = adj->con;
      MX_DENSE_PTR (W, 3, 3, blk->W);

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
      }
      else /* blk->bod == con->slave (slave on the right) */
      {
#if MPI
	right = adj->sH;
#else
	right =  adj->sprod;
#endif
      }

#if MPI
      MX_Matmat (1.0, left, MX_Tran (right), 0.0, &W);
#else
      MX_Matmat (1.0, left, right, 0.0, &W);
#endif
      SCALE9 (W.x, step);
    }

#if MPI
    /* off-diagonal external blocks */
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      MX *left, *right;
      CON *ext = (CON*)blk->dia;
      BODY *bod = blk->bod;
      MX_DENSE_PTR (W, 3, 3, blk->W);

      ASSERT_DEBUG (bod == m || bod == s, "Not connected external off-diagonal block");

      if (bod == ext->master)
      {
	right = BODY_Gen_To_Loc_Operator (bod, ext->msgp, ext->mpnt, ext->base);

	if (bod == ext->slave) /* right self-contact */
	{
	  MX *a = right,
	     *b = BODY_Gen_To_Loc_Operator (bod, ext->ssgp, ext->spnt, ext->base);

	  right = MX_Add (1.0, a, -1.0, b, NULL);
	  MX_Destroy (a);
	}
      }
      else
      {
	right = BODY_Gen_To_Loc_Operator (bod, ext->ssgp, ext->spnt, ext->base);
	MX_Scale (right, -1.0);
      }
     
      left = (bod == m ? dia->mprod : dia->sprod);

      MX_Matmat (1.0, left, MX_Tran (right), 0.0, &W);
      SCALE9 (W.x, step);
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
	if (blk->dia < dia) /* lower triangle = transposed upper triangle */
	{
	  for (blj = blk->dia->adj; blj && (blj->dia != dia || blj->bod != blk->bod); blj = blj->n); /* find upper triangle symmetric block */
	  ASSERT_DEBUG (blj, "Inconsistent W adjacency");
	  TNCOPY (blj->W, blk->W); /* transposed copy of a symmetric block */
	}
      }
    }
  }

  /* clean up */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    if (dia->mH)
    {
      MX_Destroy (dia->mH);
      MX_Destroy (dia->mprod);
      dia->mH = NULL;
    }

    if (dia->sH)
    {
      MX_Destroy (dia->sH);
      MX_Destroy (dia->sprod);
      dia->sH = NULL;
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

end:
  SOLFEC_Timer_End (ldy->dom->solfec, "LOCDYN");
}

/* updiae local dynamics => after the solution */
void LOCDYN_Update_End (LOCDYN *ldy)
{
  UPKIND upkind = update_kind (ldy->dom->solfec);

  SOLFEC_Timer_Start (ldy->dom->solfec, "LOCDYN");

  /* update cohesion states */
  if (upkind != UPPES) update_cohesion (ldy);

  SOLFEC_Timer_End (ldy->dom->solfec, "LOCDYN");
}

/* dump local dynamics to file */
void LOCDYN_Dump (LOCDYN *ldy, const char *path)
{
#define WTOL 1E-15 /* XXX */
  MEM mapmem, offmem;
  MAP *adj, *item;
  double W [9], Z;
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

    NNCOPY (dia->W, W);
    MAXABSN (W, 9, Z);
    Z *= WTOL; /* drop tolerance */
    FILTERN (W, 9, Z); /* fill with zeros below the Z tolerance */

    fprintf (f, "%s (%.6g, %.6g, %.6g) (%d, %d) [%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g] [%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g] => ",
      CON_Kind (con), con->point [0], con->point [1], con->point [2],
      (int)con->master->id, con->slave ? (int)con->slave->id : -1,
      con->base [0], con->base [1], con->base [2], con->base [3], con->base [4], con->base [5], con->base [6], con->base [7], con->base [8],
      W [0], W [1], W [2], W [3], W [4], W [5], W [6], W [7], W [8]);

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

      NNCOPY (q->W, W);
      FILTERN (W, 9, Z); /* use diagonal Z */

      fprintf (f, "%s (%.6g, %.6g, %.6g) [%.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g, %.6g] ",
        CON_Kind (con), con->point [0], con->point [1], con->point [2],
	W [0], W [1], W [2], W [3], W [4], W [5], W [6], W [7], W [8]);

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

/* export W in MatrixMarket format */
void LOCDYN_W_MatrixMarket (LOCDYN *ldy, const char *path)
{
  int i, j, n, m, nnz;
  double *W;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  FILE *f;

#if MPI
  ASSERT (0, ERR_NOT_IMPLEMENTED);
#endif

  for (dia = ldy->dia, n = nnz = 0; dia; dia = dia->n)
  {
    con = dia->con;
    con->num = n+1;
    n += 3;
    nnz += 9;
    for (blk = dia->adj; blk; blk = blk->n) nnz += 9;
  }

  ASSERT (f = fopen (path, "w"), ERR_FILE_OPEN);
  fprintf (f, "%%%%MatrixMarket matrix coordinate real general\n");
  fprintf (f, "%d  %d  %d\n", n, n, nnz);

  for (dia = ldy->dia, n = nnz = 0; dia; dia = dia->n)
  {
    con = dia->con;
    n = m = con->num;
    W = dia->W;
    for (i = 0; i < 3; i ++)
      for (j = 0; j < 3; j ++)
	fprintf (f, "%d  %d  %.15g\n", n+i, m+j, W[3*j+i]);

    for (blk = dia->adj; blk; blk = blk->n)
    {
      con = blk->dia->con;
      m = con->num;
      W = blk->W;
      for (i = 0; i < 3; i ++)
	for (j = 0; j < 3; j ++)
	  fprintf (f, "%d  %d  %.15g\n", n+i, m+j, W[3*j+i]);
    }
  }

  fclose (f);
}

/* free memory */
void LOCDYN_Destroy (LOCDYN *ldy)
{
  MEM_Release (&ldy->diamem);
  MEM_Release (&ldy->offmem);
  free (ldy);
}
