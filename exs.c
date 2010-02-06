/*
 * exs.c
 * Copyright (C) 2007, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * explicit constraints solver
 */

#include "alg.h"
#include "dom.h"
#include "bgs.h"
#include "exs.h"
#include "err.h"

/* spring and dashpot based explicit diagonal block contact solver */
int EXPLICIT_Spring_Dashpot_Contact (CON *con, double step, double gap, double spring, double dashpot, double friction,
                                     double cohesion, double *W, double *B, double *V, double *U, double *R)
{
  double INV [4], WTT[4] = {W[0], W[1], W[3], W[4]}, BN, BT [2], det, len;
  short cohesive = con->state & CON_COHESIVE;

  BN = B[2] + W[2]*R[0] + W[5]*R[1];

  R [2] = (- spring * (gap + 0.25 * step * (BN - V[2])) - 0.5 * dashpot * (BN + V[2]))
        / (1.0 + (0.25  * step * spring + 0.5 * dashpot) * W[8]);

  if (!cohesive && R[2] < 0.0)
  {
    SET (R, 0.0);
    SET (U, 0.0);
    return 0;
  }

  BT [0] = B[0] + W[6] * R[2];
  BT [1] = B[1] + W[7] * R[2];

  INVERT2(WTT, INV, det);

  if (det == 0.0) return -1;

  R [0] = - INV[0] * BT[0] - INV[2] * BT[1];
  R [1] = - INV[1] * BT[0] - INV[3] * BT[1];

  if (cohesive && R [2] < -cohesion * con->area)
  {
    cohesive = 0;
    con->state &= ~CON_COHESIVE;
    SURFACE_MATERIAL_Cohesion_Set (&con->mat, 0.0);
    R [2] = -cohesion * con->area;
  }

  len = LEN2 (R);
  det = friction * fabs (R[2]);

  if (len > det)
  {
    R[0] *= det / len;
    R[1] *= det / len;

    if (cohesive)
    {
      cohesive = 0;
      con->state &= ~CON_COHESIVE;
      SURFACE_MATERIAL_Cohesion_Set (&con->mat, 0.0);
    }

    /* TODO: it would be more "implicit" to solve local nonlinear problem (I + len/|UT| WTT) UT = BT,
     * TODO: but then the local iterations might fail; explicit stuff does not fail like this ... */
  }

  NVADDMUL (B, W, R, U); /* U = B + W R */

  return 0;
}

#if MPI
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

/* receive non-contact reactions */
static void receive_noncontact_reactions (DOM *dom, COMDATA *recv, int nrecv)
{
  COMDATA *ptr;
  int i, j, *k;
  double *R;
  CON *con;

  for (i = 0, ptr = recv; i < nrecv; i ++, ptr ++)
  {
    for (j = 0, k = ptr->i, R = ptr->d; j < ptr->ints; j ++, k ++, R += 3)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*k), NULL), "Invalid constraint id");
      COPY (R, con->R);
    }
  }
}
#endif

/* explcit constraint solver */
void EXPLICIT_Solve (LOCDYN *ldy)
{
  double error, step;
  short dynamic;
  int iters;
  DIAB *dia;
  CON *con;

  step = ldy->dom->step;

  /* first explicitly process contacts */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    con = dia->con;

    if (con->kind == CONTACT)
    {
      EXPLICIT_Spring_Dashpot_Contact (con, step, con->gap, con->mat.base->spring, con->mat.base->dashpot,
	         con->mat.base->friction, con->mat.base->cohesion, dia->W, dia->B, dia->V, dia->U, dia->R);
    }
  }

#if MPI
  /* send contact reactions (well also other among them, minority) */
  DOM_Update_External_Reactions (ldy->dom, 0);
#endif

  /* Gauss-Seidel sweep for non-contacts now */
  dynamic = ldy->dom->dynamic;
  iters = 0;
  do
  {
    double errup = 0.0,
	   errlo = 0.0;
    OFFB *blk;
    DIAB *dia;
   
    for (dia = ldy->dia; dia; dia = dia->n)
    {
      CON *con = dia->con;

      if (con->kind == CONTACT) continue; /* skip contacts */

      double R0 [3],
	     B [3],
	     *R = dia->R;

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
	CON *con = (CON*) blk->dia;
	double *W = blk->W,
	       *R = con->R;
	NVADDMUL (B, W, R, B);
      }
#endif
      
      COPY (R, R0); /* previous reaction */

      /* solve local diagonal block problem */
      DIAGONAL_BLOCK_Solver (GS_PROJECTED_GRADIENT, 1E-6, 1000, dynamic, step,
	         con->kind, con->mat.base, con->gap, con->Z, con->base, dia, B);

      /* accumulate relative
       * error components */
      SUB (R, R0, R0);
      errup += DOT (R0, R0);
      errlo += DOT (R, R);
    }

    /* calculate relative error */
    error = sqrt (errup) / sqrt (MAX (errlo, 1.0));
  }
  while (++ iters < 1000 && error > 1E-3);

  ASSERT_DEBUG (iters < 1000 && error < 1E-3, "Gauss-Seidel part of EXPLICIT_SOLVER not convergent");

#if MPI
  COMDATA *send, *recv, *ptr;
  int nsend, nrecv, size;
  SET *ranks, *item;
  MEM setmem;
  OFFB *blk;

  size = 128;
  nsend = 0;
  ERRMEM (send = MEM_CALLOC (size * sizeof (COMDATA)));
  MEM_Init (&setmem, sizeof (SET), 128);
  ptr = send;

  /* create send sets for non-contact reactions */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;

    if (con->kind == CONTACT) continue;

    for (ranks = NULL, blk = dia->adjext; blk; blk = blk->n)
    { 
      CON *con = (CON*) blk->dia;
      SET_Insert (&setmem, &ranks, (void*) (long) con->rank, NULL);
    }

    for (item = SET_First (ranks); item; item = SET_Next (item))
    {
      ptr->rank = (int) (long) item->data;
      ptr->ints = 1;
      ptr->doubles = 3;
      ptr->i = (int*) &con->id;
      ptr->d = con->R;
      ptr= sendnext (++ nsend, &size, &send);
    }

    SET_Free (&setmem, &ranks);
  }

  /* send non-contact reactions */
  COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);

  /* receive non-contact reations */
  receive_noncontact_reactions (ldy->dom, recv, nrecv);

  MEM_Release (&setmem);
  free (send);
  free (recv);
#endif
}

/* write labeled satate values */
void EXPLICIT_Write_State (GAUSS_SEIDEL *gs, PBF *bf)
{
  /* nothing to write for the moment */
}
