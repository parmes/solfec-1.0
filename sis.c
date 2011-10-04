/*
 * sis.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Siconos 3D contact solvers interface
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

#include <FrictionContact3D_Solvers.h>
#include "sis.h"
#include "dom.h"
#include "alg.h"
#include "err.h"

#define COPY_REV(a, b)\
  (b) [0] = (a) [2];\
  (b) [1] = (a) [1];\
  (b) [2] = (a) [0]

/* translate local dynamics into Siconos FrictionContactProblem */
static FrictionContactProblem* LOCDYN2FCP (LOCDYN *ldy)
{
  SparseBlockStructuredMatrix *W;
  FrictionContactProblem *fc;
  double *q, *mu, **b, *b0;
  NumericsMatrix *M;
  size_t *pp, *ii;
  int i, j, k;
  DIAB *dia;
  OFFB *blk;
  CON *con;

  for (i = j = 0, dia = ldy->dia; dia; dia = dia->n)
  {
    con = dia->con;
    if (con->kind == CONTACT)
    {
      con->num = i ++; /* i = number of contact points */
      j ++; /* j = number of blocks; diagonal block */
      for (blk = dia->adj; blk; blk = blk->n)
      {
	CON *adj = blk->dia->con;
	if (adj->kind == CONTACT) j ++; /* off diagonal block */
      }
    }
  }

  if (i == 0) return NULL;

  ERRMEM (fc = malloc (sizeof (FrictionContactProblem)));
  fc->dimension = 3;
  fc->numberOfContacts = i;
  ERRMEM (fc->q = malloc (sizeof (double [3 * i])));
  ERRMEM (fc->mu = malloc (sizeof (double [i])));

  for (dia = ldy->dia, q = fc->q, mu = fc->mu; dia; dia = dia->n)
  {
    CON *con = dia->con;
    if (con->kind == CONTACT)
    {
      double *B = dia->B;
      COPY_REV (B, q);
      *mu = con->mat.base->friction;
      q += 3;
      mu ++;
    }
  }

  ERRMEM (M = malloc (sizeof (NumericsMatrix)));
  ERRMEM (W = malloc (sizeof (SparseBlockStructuredMatrix)));
  M->storageType = 1;
  M->size0 = 3 * i;
  M->size1 = M->size0;
  M->matrix0 = NULL;
  M->matrix1 = W;
  fc->M = M;

  ERRMEM (W->block = malloc (j * sizeof (double*)));
  ERRMEM (b0 = malloc (j * sizeof (double [9])));
  for (dia = ldy->dia, b = W->block; dia; dia = dia->n) /* set block pointers */
  {
    con = dia->con;
    if (con->kind == CONTACT)
    {
      double *W = dia->W;
      for (k = 0; k < 9; k ++) b0 [k] = W[8-k];
      *b = b0;
      b ++;
      b0 += 9;
      for (blk = dia->adj; blk; blk = blk->n)
      {
	if (blk->dia->con->kind == CONTACT)
	{
	  double *W = blk->W;
          for (k = 0; k < 9; k ++) b0 [k] = W[8-k];
	  *b = b0;
	  b ++;
	  b0 += 9;
	}
      }
    }
  }

  W->nbblocks = j;
  W->blocknumber0 = i;
  W->blocknumber1 = i;
  ERRMEM (W->blocksize0 = malloc (i * sizeof (int)));
  ERRMEM (W->blocksize1 = malloc (i * sizeof (int)));

  for (k = 0; k < i; k ++)
  {
    W->blocksize0 [k] = (k+1) * 3;
    W->blocksize1 [k] = (k+1) * 3;
  }

  W->filled1 = i + 1;
  W->filled2 = j;

  ERRMEM (W->index1_data = malloc (W->filled1 * sizeof (int)));
  ERRMEM (W->index2_data = malloc (W->filled2 * sizeof (int)));

  pp = W->index1_data;
  *pp = 0;
  for (dia = ldy->dia, pp ++, ii = W->index2_data; dia; dia = dia->n) /* set compressed block row indexing */
  {
    con = dia->con;
    if (con->kind == CONTACT)
    {
      *ii = con->num;
      ii ++;
      k = 1;
      for (blk = dia->adj; blk; blk = blk->n)
      {
	CON *adj = blk->dia->con;
	if (adj->kind == CONTACT)
	{
	  *ii = adj->num; /* block column indices */
	  ii ++;
	  k ++;
	}
      }
      *pp = (*(pp-1)) + k; /* block row pointers */
      pp ++;
    }
  }
  
  return fc;
}

/* free the problem */
static void FCP_free (FrictionContactProblem *fc)
{
  free (fc->q);
  free (fc->mu);
  free (fc->M->matrix1->block [0]);
  free (fc->M->matrix1->block);
  free (fc->M->matrix1->blocksize0);
  free (fc->M->matrix1->blocksize1);
  free (fc->M->matrix1->index1_data);
  free (fc->M->matrix1->index2_data);
  free (fc->M->matrix1);
  free (fc->M);
  free (fc);
}

/* create solver */
SICONOS* SICONOS_Create (double epsilon, int maxiter)
{
  SICONOS *si;

  ERRMEM (si = MEM_CALLOC (sizeof (SICONOS)));
  si->epsilon = epsilon;
  si->maxiter = maxiter;

  return si;
}

/* run solver */
void SICONOS_Solve (SICONOS *si, LOCDYN *ldy)
{
  FrictionContactProblem *fc;
  double *R, *U, *R1, *U1;
  NumericsOptions no;
  SolverOptions so;
  DIAB *dia;
  CON *con;

  fc = LOCDYN2FCP (ldy);

  if (fc)
  {
    ERRMEM (R = malloc (sizeof (double [fc->M->size0])));
    ERRMEM (U = malloc (sizeof (double [fc->M->size0])));

    for (dia = ldy->dia, R1 = R, U1 = U; dia; dia = dia->n)
    {
      con = dia->con;
      if (con->kind == CONTACT)
      {
	double *R0 = con->R,
	       *U0 = con->U;
	COPY_REV (R0, R1);
	COPY_REV (U0, U1);
	R1 += 3;
	U1 += 3;
      }
    }

    no.verboseMode = si->verbose;
    frictionContact3D_setDefaultSolverOptions(&so, SICONOS_FRICTION_3D_NSGS);
    so.iparam [0] = si->maxiter;
    so.dparam [0] = si->epsilon;

    frictionContact3D_driver (fc, R, U, &so, &no);

    deleteSolverOptions (&so);

    for (dia = ldy->dia, R1 = R, U1 = U; dia; dia = dia->n)
    {
      con = dia->con;
      if (con->kind == CONTACT)
      {
	double *R0 = con->R,
	       *U0 = con->U;
	COPY_REV (R1, R0);
	COPY_REV (U1, U0);
	R1 += 3;
	U1 += 3;
      }
    }

    FCP_free (fc);
    free (R);
    free (U);
  }
}

/* write labeled state values */
void SICONOS_Write_State (SICONOS *si, PBF *bf)
{
}

/* destroy solver */
void SICONOS_Destroy (SICONOS *si)
{
  free (si);
}
