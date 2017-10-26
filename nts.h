/*
 * nts.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * projected quasi-Newton solver
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

#include "ldy.h"

#ifndef __nts__
#define __nts__

typedef struct newton NEWTON;

struct newton
{
  /* input */

  double meritval; /* value of merit function sufficient for termination */

  int maxiter; /* iterations bound */

  enum {LOCDYN_ON, LOCDYN_OFF} locdyn; /* local dynamics assembling */

  enum {PQN_GMRES, PQN_DIAG} linver; /* linear solver version */

  int linmaxiter; /* linear solver iterations bound */

  int maxmatvec; /* matrix-vector products bound */

  double epsilon; /* linear solver epsilon */

  double delta; /* diagonal regularization */

  double theta; /* relaxation parameter for PQN_DIAG version of the solver */

  double omega; /* smoothing omega */

  short gsflag;

  /* output */

  double *merhist; /* merit function history */

  int *mvhist; /* matrix-vector products history */

  int iters; /* iterations count */

  int *itershist; /* iterations history of all solver calls */

  int itershistcount; /* < 0: do not collect iterations history; >= 0 the size of the iterhistory */

  int itershistsize; /* iterations history buffer size */
};

/* create solver */
NEWTON* NEWTON_Create (double meritval, int maxiter);

/* run solver */
void NEWTON_Solve (NEWTON *ns, LOCDYN *ldy);

/* write labeled state values */
void NEWTON_Write_State (NEWTON *ns, PBF *bf);

/* destroy solver */
void NEWTON_Destroy (NEWTON *bs);

#endif
