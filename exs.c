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
int EXPLICIT_Spring_Dashpot_Contact (CON *con, double gap, double spring, double dashpot, double friction,
                                     double cohesion, double *W, double *B, double *V, double *U, double *R)
{
  double INV [4], WTT[4] = {W[0], W[1], W[3], W[4]}, BT [2], det, len;
  short cohesive = con->state & CON_COHESIVE;

  if (!cohesive && gap >= 0)
  {
    SET (R, 0);
    COPY (B, U);
    return 0;
  }

  R [2] = - spring * gap - 0.5 * dashpot * (B[2] + W[2]*R[0] + W[5]*R[1] + V[2]) / (1.0 + 0.5 * dashpot * W[8]);

  if (!cohesive && R[2] < 0.0) R [2] = 0.0;

  BT [0] = B[0] + W[6] * R[2];
  BT [1] = B[1] + W[7] * R[2];

  INVERT2(WTT, INV, det);

  if (det == 0.0) return -1;

  R [0] = - INV[0] * BT[0] - INV[2] * BT[1];
  R [1] = - INV[1] * BT[0] - INV[3] * BT[1];

  if (cohesive && fabs (R[2]) > cohesion)
  {
    cohesive = 0;
    con->state &= ~CON_COHESIVE;
    SURFACE_MATERIAL_Cohesion_Set (&con->mat, 0.0);
  }

  if (!cohesive)
  {
    len = LEN2 (R);
    det = friction * R[2];

    if (len > det)
    {
      R[0] *= det / len;
      R[1] *= det / len;

      /* TODO: it would be more "implicit" to solve local nonlinear problem (I + len/|UT| WTT) UT = BT,
       * TODO: but then the local iterations might fail; explicit stuff does not fail like this ... */
    }
  }

  return 0;
}

/* explcit constraint solver */
void EXPLICIT_Solve (LOCDYN *ldy)
{
  double error, step;
  short dynamic;
  int iters;
  DIAB *dia;
  CON *con;

  /* first explicitly process contacts */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    con = dia->con;

    if (con->kind == CONTACT)
    {
      EXPLICIT_Spring_Dashpot_Contact (con, con->gap, con->mat.base->spring, con->mat.base->dashpot,
	  con->mat.base->friction, con->mat.base->cohesion, dia->W, dia->B, dia->V, dia->U, dia->R);
    }
  }

  /* Gauss-Seidel sweep for non-contacts now */
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
  DOM_Update_External_Reactions (ldy->dom, 0);
#endif
}
