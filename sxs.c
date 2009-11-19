/*
 * sxs.c
 * Copyright (C) 2007, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * semi-explicit constraints solver
 */

#include "alg.h"
#include "dom.h"
#include "exs.h"
#include "bgs.h"

#define EPSILON 1E-9

/* diagonal block solver */
static void solver (short dynamic, double step, short kind, SURFACE_MATERIAL *mat, double gap, double *Z, double *base, DIAB *dia, double *B)
{
  if (DIAGONAL_BLOCK_Solver (GS_PROJECTED_GRADIENT, EPSILON, 32, dynamic, step, kind, mat, gap, Z, base, dia, B) >= 32)
  {
    if (DIAGONAL_BLOCK_Solver (GS_SEMISMOOTH_NEWTON, EPSILON, 16, dynamic, step, kind, mat, gap, Z, base, dia, B) >= 16)
    {
      if (DIAGONAL_BLOCK_Solver (GS_DE_SAXE_AND_FENG, EPSILON, 1024, dynamic, step, kind, mat, gap, Z, base, dia, B) >= 1024)
      {
	SET (dia->R, 0.0); /* failed to solve - leave it for now */
      }
    }
  }
}

/* semi-explcit constraint solver */
void SEMI_EXPLICIT_Solve (LOCDYN *ldy)
{
  short dynamic;
  double step;
  DIAB *dia;
  CON *con;

  dynamic = DOM(ldy->dom)->dynamic;
  step = DOM(ldy->dom)->step;

#if MPI 
  for (dia = ldy->diab; dia; dia = dia->n) /* use balanced blocks */
#else
  for (dia = ldy->dia; dia; dia = dia->n)
#endif
  {
    con = dia->con;

#if MPI
    if (con) solver (dynamic, step, con->kind, &con->mat, con->gap, con->Z, con->base, dia, dia->B); /* LDB_OFF */
    else solver (dynamic, step, dia->kind, &dia->mat, dia->gap, dia->Z, dia->base, dia, dia->B);
#else
    solver (dynamic, step, con->kind, &con->mat, con->gap, con->Z, con->base, dia, dia->B);
#endif
  }
}
