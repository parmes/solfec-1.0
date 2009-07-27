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

#include "alg.h"
#include "dom.h"
#include "ldy.h"
#include "lap.h"
#include "err.h"

/* memory block size */
#define BLKSIZE 512

/* apply forward change of variables (nornal
 * contact forces) due to the cohesion, etc. */
static void variables_change_begin (LOCDYN *ldy)
{
  OFFB *blk;
  DIAB *dia;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;
    double *B = dia->B; /* free velocity will
			   be eventually modified */

    if (con->state & CON_COHESIVE) /* cohesive state */
    {
      double c = con->mat.cohesion * con->area,
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
	double c = con->mat.cohesion * con->area,
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
    short state = con->state;

    if (state & CON_COHESIVE) /* cohesive state */
    {
      double c = con->mat.cohesion * con->area,
	     *R = dia->R;

      R [2] -= c; /* back change */

      if ((state & CON_OPEN) || /* mode-I decohesion */
	(!(state & CON_STICK))) /* mode-II decohesion */
      {
	con->state &= ~CON_COHESIVE;
	con->mat.cohesion = 0.0;
      }
    }
  }
}

/* create local dynamics for a domain */
LOCDYN* LOCDYN_Create (void *dom)
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
DIAB* LOCDYN_Insert (LOCDYN *ldy, void *con, BODY *one, BODY *two)
{
  DIAB *dia, *nei;
  SET *item;
  OFFB *b;
  CON *c;

  ERRMEM (dia = MEM_Alloc (&ldy->diamem));
  dia->R = ((CON*)con)->R;
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
      if (item->data != con) /* skip the coincident constraint */
      {
	c = item->data;
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
      if (item->data != con) /* skip the coincident constraint */
      {
	c = item->data;
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

  /* mark as modified */
  ldy->modified = 1;

  return dia;
}

/* remove a diagonal entry from local dynamics */
void LOCDYN_Remove (LOCDYN *ldy, DIAB *dia)
{
  OFFB *b, *c, *r;

  /* destroy blocks in
   * adjacent diaa items */
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

  /* remove from list */
  if (dia->p)
    dia->p->n = dia->n;
  else ldy->dia = dia->n;
  if (dia->n)
    dia->n->p = dia->p;

  /* destroy passed diaa */
  MEM_Free (&ldy->diamem, dia);

  /* mark as modified */
  ldy->modified = 1;
}

/* updiae local dynamics => prepare for a solution */
void LOCDYN_Update_Begin (LOCDYN *ldy, UPKIND upkind)
{
  DOM *dom = ldy->dom;
  double step = dom->step;
  DIAB *dia;

  /* calculate local velocities and
   * assmeble the force-velocity 'W' operator */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;
    BODY *m = con->master,
	 *s = con->slave;
    void *mgobj = con->mgobj,
	 *sgobj = con->sgobj;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base,
	   *V = dia->V,
	   *B = dia->B,
           X [3], Y [9];
    MX_DENSE_PTR (W, 3, 3, dia->W);
    MX_DENSE (C, 3, 3);
    MX *mH, *sH;
    OFFB *blk;

    /* previous time step velocity */
    BODY_Local_Velo (m, PREVELO, mgobj, mpnt, base, X); /* master body pointer cannot be NULL */
    if (s) BODY_Local_Velo (s, PREVELO, sgobj, spnt, base, Y); /* might be NULL for some constraints (one body) */
    else { SET (Y, 0.0); }
    SUB (Y, X, V); /* relative = slave - master => outward master normal */

    /* local free velocity */
    BODY_Local_Velo (m, CURVELO, mgobj, mpnt, base, X);
    if (s) BODY_Local_Velo (s, CURVELO, sgobj, spnt, base, Y);
    else { SET (Y, 0.0); }
    SUB (Y, X, B);

    /* diagonal block */
    mH = BODY_Gen_To_Loc_Operator (m, mgobj, mpnt, base);
    MX_Trimat (mH, m->inverse, MX_Tran (mH), &W); /* H * inv (M) * H^T */
    if (s)
    { sH = BODY_Gen_To_Loc_Operator (s, sgobj, spnt, base);
      MX_Trimat (sH, s->inverse, MX_Tran (sH), &C); /* H * inv (M) * H^T */
      NNADD (W.x, C.x, W.x); }
    SCALE9 (W.x, step); /* W = h * ( ... ) */

    NNCOPY (W.x, C.x); /* calculate regularisation parameter */
    ASSERT (lapack_dsyev ('N', 'U', 3, C.x, 3, X, Y, 9) == 0, ERR_LDY_EIGEN_DECOMP);
    dia->rho = 1.0 / X [2]; /* inverse of maximal eigenvalue */

    /* off-diagonal blocks if requested */
    for (blk = dia->adj; upkind == UPALL && blk; blk = blk->n)
    {
      DIAB *dia = blk->dia;
      CON *con = dia->con;
      BODY *bod = blk->bod;
      MX *lH, *rH, *inv;
      MX_DENSE_PTR (W, 3, 3, blk->W);
      double coef;

      ASSERT_DEBUG (bod == m || bod == s, "Off diagonal block is not connected!");
     
      lH = (bod == m ? mH : sH); /* dia->bod is a valid body (not an obstacle)
                                   as it was inserted into the dual graph */
      inv = bod->inverse;

      if (bod == con->master)
      {
	rH =  BODY_Gen_To_Loc_Operator (bod, con->mgobj, con->mpnt, con->base);
	coef = (bod == s ? -step : step);
      }
      else /* blk->bod == dia->slave */
      {
	rH =  BODY_Gen_To_Loc_Operator (bod, con->sgobj, con->spnt, con->base);
	coef = (bod == m ? -step : step);
      }

      MX_Trimat (lH, inv, MX_Tran (rH), &W);
      SCALE9 (W.x, coef);
      MX_Destroy (rH);
    }

    MX_Destroy (mH);
    if (s) MX_Destroy (sH);
  }

  /* forward variables change */
  variables_change_begin (ldy);
}

/* updiae local dynamics => after the solution */
void LOCDYN_Update_End (LOCDYN *ldy)
{
  /* backward variables change */
  variables_change_end (ldy);

  /* not modified */
  ldy->modified = 0;
}

/* free memory */
void LOCDYN_Destroy (LOCDYN *ldy)
{
  MEM_Release (&ldy->diamem);
  MEM_Release (&ldy->offmem);
  free (ldy);
}
