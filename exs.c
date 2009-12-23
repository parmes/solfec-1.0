/*
 * exs.c
 * Copyright (C) 2007, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * explicit constraints solver
 */

#include "alg.h"
#include "dom.h"
#include "exs.h"

static int fixpnt (double spring, double *X, double *x, double *W, double *B, double *U, double *R)
{
  SUB (x, X, R);
  SCALE (R, spring);
  NVADDMUL (B, W, R, U);
  return 0;
}

static int fixdir (double spring, double *X, double *x, double *n, double *W, double *B, double *U, double *R)
{
  double len;

  SUB (x, X, R);
  len = DOT (n, R);
  SUBMUL (R, len, n, R);
  SCALE (R, spring);
  NVADDMUL (B, W, R, U);
  return 0;
}

static int velodir (double *Z, double *W, double *B, double *U, double *R)
{
  R [0] = 0.0;
  R [1] = 0.0;
  R [2] = (VELODIR(Z) - B[2]) / W[8];
  U [0] = B [0]; 
  U [1] = B [1]; 
  U [2] = VELODIR(Z);
  return 0;
}

static int riglnk (double gap, double spring, double *n, double *W, double *B, double *U, double *R)
{
  double force;

  force = gap * spring;
  COPY (n, R);
  SCALE (R, force);
  NVADDMUL (B, W, R, U);
  return 0;
}

static int solver (short kind, SURFACE_MATERIAL *mat, double gap, double *Z,
                   double *mpnt, double *point, double *base, DIAB *dia)
{
  switch (kind)
  {
  case CONTACT:
    return EXPLICIT_Spring_Dashpot_Contact (gap, mat->spring, mat->dashpot,
	                       mat->friction, dia->W, dia->B, dia->V, dia->U, dia->R);
  case FIXPNT:
    return fixpnt (mat->spring, mpnt, point, dia->W, dia->B, dia->U, dia->R);
  case FIXDIR:
    return fixdir (mat->spring, mpnt, point, base+6, dia->W, dia->B, dia->U, dia->R);
  case VELODIR:
    return velodir (Z, dia->W, dia->B, dia->U, dia->R);
  case RIGLNK:
    return riglnk (gap, mat->spring, base+6, dia->W, dia->B, dia->U, dia->R);
  }

  return 0;
}

/* spring and dashpot based explicit diagonal block contact solver */
int EXPLICIT_Spring_Dashpot_Contact (double gap, double spring, double dashpot,
         double friction, double *W, double *B, double *V, double *U, double *R)
{
  double INV [4], WTT[4] = {W[0], W[1], W[3], W[4]}, BT [2], det, len;

  if (gap >= 0)
  {
    SET (R, 0);
    COPY (B, U);
    return 0;
  }

  R [2] = - spring * gap - dashpot * V[2];

  BT [0] = B[0] + W[6] * R[2];
  BT [1] = B[1] + W[7] * R[2];

  INVERT2(WTT, INV, det);

  if (det == 0.0) return -1;

  R [0] = - INV[0] * BT[0] - INV[2] * BT[1];
  R [1] = - INV[1] * BT[0] - INV[3] * BT[1];

  len = LEN2 (R);
  det = friction * R[2];

  if (len > det)
  {
    /* TODO: solve local nonlinear problem (I + len/|UT| WTT) UT = BT */
    R[0] *= det / len;
    R[1] *= det / len;
  }

  /* TODO: make it aware of cohesion */

  return 0;
}

/* explcit constraint solver */
void EXPLICIT_Solve (LOCDYN *ldy)
{
  DIAB *dia;
  CON *con;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    con = dia->con;

    /* FIXME: non-contacts will have NULL material */

    solver (con->kind, con->mat.base, con->gap, con->Z, con->mpnt, con->point, con->base, dia);
  }

  /* TODO: use Guss-Seidel for non-contacts (not RIGLNK) and penalty for rest */
}
