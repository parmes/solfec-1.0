/*
 * pes.c
 * Copyright (C) 2007-2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * penalty constraints solver
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
#include "bgs.h"
#include "pes.h"
#include "err.h"

#define PENALTY_MAXITER 1000
#define PENALTY_EPSILON 1E-4 /* XXX */

/* spring and dashpot based explicit diagonal block contact solver */
int PENALTY_Spring_Dashpot_Contact (CON *con, short implicit, double step, double gap, double spring, double dashpot,
                              double friction, double cohesion, double *W, double *B, double *V, double *U, double *R)
{
  double INV [4], WTT[4] = {W[0], W[1], W[3], W[4]}, BN, BT [2], det, len;
  short cohesive = con->state & CON_COHESIVE;

  BN = B[2] + W[2]*R[0] + W[5]*R[1];

  if (implicit)
  {
    R [2] = (- spring * (gap + 0.25 * step * (BN - V[2])) - 0.5 * dashpot * (BN + V[2]))
	  / (1.0 + (0.25  * step * spring + 0.5 * dashpot) * W[8]);
  }
  else
  {
    R [2] = - spring * gap - dashpot * V[2];
  }

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

/* create penalty solver */
PENALTY* PENALTY_Create (short implicit)
{
  PENALTY *ps;

  ERRMEM (ps = MEM_CALLOC (sizeof (PENALTY)));
  ps->implicit = implicit;

  return ps;
}

/* explcit constraint solver */
void PENALTY_Solve (PENALTY *ps, LOCDYN *ldy)
{
  GAUSS_SEIDEL *gs;
  short implicit;
  double step;
  LOCDYN *clo;
  DIAB *dia;
  CON *con;

  implicit = ps->implicit;
  step = ldy->dom->step;

#if MPI
  if (ldy->dom->rank == 0)
#endif
  if (ldy->dom->verbose) printf ("PENALTY_SOLVER: applying springs and dashpots...\n");

  /* first explicitly process contacts */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    con = dia->con;

    if (con->kind == CONTACT)
    {
      SURFACE_MATERIAL *bas = con->mat.base;
      PENALTY_Spring_Dashpot_Contact (con, implicit, step, con->gap, bas->spring, bas->dashpot,
	                  bas->friction, bas->cohesion, dia->W, dia->B, dia->V, dia->U, dia->R);
    }
  }

#if MPI
  if (ldy->dom->rank == 0)
#endif
  if (ldy->dom->verbose) printf ("PENALTY_SOLVER: solving non-contact constraints...\n");

#if MPI
  DOM_Update_External_Reactions (ldy->dom, 0);
#endif
  clo = LOCDYN_Clone_Non_Contacts (ldy);
  gs = GAUSS_SEIDEL_Create (PENALTY_EPSILON, PENALTY_MAXITER, 1.0, GS_FAILURE_CONTINUE, 1E-9, 100, DS_SEMISMOOTH_NEWTON, NULL, NULL);
  gs->nomerit = 1;
  GAUSS_SEIDEL_Solve (gs, clo);
  GAUSS_SEIDEL_Destroy (gs);
  LOCDYN_Destroy (clo);
}

/* write labeled satate values */
void PENALTY_Write_State (PENALTY *ps, PBF *bf)
{
  /* nothing to write for the moment */
}

/* destroy penalty solver */
void PENALTY_Destroy (PENALTY *ps)
{
  free (ps);
}
