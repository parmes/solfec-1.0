/*
 * pqns.cu
 * Copyright (C) 2010, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
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

#include <iostream>
#include "alg.h"
using std::cerr;
using std::endl;

/* compressed row block matrix - vector product */
__global__ void U_WR_B (const int num, const int *ptr, const int *adj,
  const float *W0, const float *R0, const float *B0, float *U0)
{
  int con = blockDim.x * blockIdx.x + threadIdx.x;

  if (con < num)
  {
    const float *W = W0;
    const float *B = &B0[con*3];
    float WR[3] = {B[0], B[1], B[2]};
    int j, j0 = ptr [con], j1 = ptr [con+1];
    for (j = j0; j < j1; j ++, W += 9)
    { 
      const float *R = &R0[adj[j]*3];
      NVADDMUL (WR, W, R, WR);
    }
    float *U = &U0[con*3];
    COPY (WR, U);
  }
}

/* constraint kinds */
enum {_CONTACT_ = 0, _FIXPNT_, _FIXDIR_, _VELODIR_, _RIGLNK_, _GLUE_};

/* diagonal solver */
__global__ void solve (const int dynamic, const float step, const float theta, const float eps,
  const int num, const int *kind, const float *mat0, const int *ptr, const float *W0,
  const float *V0, const float *U0, float *DR0, float *R0)
{
  int con = blockDim.x * blockIdx.x + threadIdx.x;

  const float *U = &U0[con*3],
              *V = &V0[con*3],
              *mat = &mat0[con*4],
              *W = &W0[ptr[con]*9];

  float *R = &R0[con*3],
        *DR = &DR0[con*3];

  float x [3], T [9], gamma = 1.0 - theta;

  switch (kind [con])
  {
  case _FIXPNT_:
  case _GLUE_:
  {
    if (dynamic)
    {
      x [0] = -V[0]-U[0];
      x [1] = -V[1]-U[1];
      x [2] = -V[2]-U[2];
    }
    else
    {
      x [0] = -U[0];
      x [1] = -U[1];
      x [2] = -U[2];
    }

    NNCOPY (W, T);
  }
  break;
  case _FIXDIR_:
  {
    x [0] = -R[0];
    x [1] = -R[1];
    if (dynamic) x [2] = -V[2]-U[2];
    else x [2] = -U[2];

    T [1] = T [3] = T [6] = T [7] = 0.0;
    T [0] = T [4] = 1.0;
    T [2] = W [2];
    T [5] = W [5];
    T [8] = W [8];
  }
  break;
  case _VELODIR_:
  {
    x [0] = -R[0];
    x [1] = -R[1];
    x [2] = mat[0]-U[2];

    T [1] = T [3] = T [6] = T [7] = 0.0;
    T [0] = T [4] = 1.0;
    T [2] = W [2];
    T [5] = W [5];
    T [8] = W [8];
  }
  break;
  case _RIGLNK_:
  {
    float h = step * (dynamic ? 0.5 : 1.0),
	  d = mat [0],
	  delta;

    x [0] = -R[0];
    x [1] = -R[1];
    delta = d*d - h*h*DOT2(U,U);
    if (delta >= 0.0) x [2] = (sqrt (delta) - d)/h - U[2];
    else x [2] = -U[2];

    T [1] = T [3] = T [6] = T [7] = 0.0;
    T [0] = T [4] = 1.0;
    T [2] = W [2];
    T [5] = W [5];
    T [8] = W [8];
  }
  break;
  case _CONTACT_:
  {
    float X [9], Y [9], dF [9], S [3], F [3], m [3],
	  fri = mat [0],
	  res = mat [1],
	  coh = mat [2],
	  gap = mat [3];

    float udash, ulen, sdot, slen, l1, l2, u1[3], u2[3], eps2,
	 fri2, onefri2, sq1, sq2, g1, g2, dg1, dg2, a, b, c, d;

    eps2 = eps*eps;
    fri2 = fri*fri;
    onefri2 = 1.0 + fri2;
    if (dynamic) udash = (U[2] + res * MIN (V[2], 0));
    else udash = ((MAX(gap, 0)/step) + U[2]);
    ulen = sqrt (DOT2(U, U) + eps2);

    F [0] = U[0];
    F [1] = U[1];
    F [2] = (udash + fri * (ulen - eps));

    SUB (R, F, S);
    S [2] += coh;

    sdot = DOT2 (S, S);
    slen = sqrt (sdot);
    l1 = -(S[2] + fri*slen) / onefri2;
    l2 =  (slen - fri*S[2]) / onefri2;
    if (slen != 0.0)
    {
      u2[0] = S[0]/slen;
      u2[1] = S[1]/slen;
      u2[2] = -fri;
      u1[0] = -fri*u2[0];
      u1[1] = -fri*u2[1];
      u1[2] = -1.0;
    }
    else
    {
      u2[0] =  1.0;
      u2[1] =  0.0;
      u2[2] = -fri;
      u1[0] = -fri;
      u1[1] =  0.0;
      u1[2] = -1.0;
    }
    sq1 = sqrt (l1*l1 + 4.0*eps2);
    sq2 = sqrt (l2*l2 + 4.0*eps2*fri2);
    g1 = 0.5 * (sq1 + l1);
    g2 = 0.5 * (sq2 + l2);

    m [0] = g1*u1[0] + g2*u2[0];
    m [1] = g1*u1[1] + g2*u2[1];
    m [2] = g1*u1[2] + g2*u2[2];

    ADD (F, m, x);
    SCALE (x, -1.0);

    dF [1] = dF [3] = dF [6] = dF [7] = 0.0;
    dF [0] = dF [4] = dF [8] = 1.0;
    dF [2] = fri * U[0] / ulen;
    dF [5] = fri * U[1] / ulen;

    dg1 = 0.5*(1.0+l1/sq1);
    dg2 = 0.5*(1.0+l2/sq2);
    a = 0.5*(1.0+(l2 + fri*l1)/(sq2 + fri*sq1));
    b = (fri2 * dg1 + dg2) / onefri2;
    c = (fri * (dg1 - dg2)) / onefri2;
    d = (dg1 + fri2 * dg2) / onefri2;

    if (slen != 0.0)
    {
      Y [0] = a + (b - a) * S[0]*S[0] / sdot;
      Y [1] = (b - a) * S[1]*S[0] / sdot;
      Y [2] = c * S[0] / slen;
      Y [3] = Y[1];
      Y [4] = a + (b - a) * S[1]*S[1] / sdot;
      Y [5] = c * S[1] / slen;
      Y [6] = Y[2];
      Y [7] = Y[5];
      Y [8] = d;
    }
    else
    {
      Y[1] = Y[2] = Y[3] = Y[5] = Y[6] = Y[7] = 0.0;
      Y[0] = Y[4] = Y[8] = dg1;
    }

    NNMUL (Y, dF, X);
    NNSUB (dF, X, X); /* X = [I - dm/dS] dF/dU */

    NNMUL (X, W, T);
    NNADD (T, Y, T);
  }
  break;
  }

  /* 3x3 Gauss elimination */
  T [3] /= T[0]; T [6] /= T[0]; x [0] /= T[0];
  T [4] -= T[3]*T[1]; T [7] -= T[6]*T[1]; x [1] -= x[0]*T[1];
  T [5] -= T[3]*T[2]; T [8] -= T[6]*T[2]; x [2] -= x[0]*T[2];
  T [7] /= T [4]; x [1] /= T[4];
  T [8] -= T[7]*T[5]; x [2] -= x[1]*T[5];
  x [2] /= T [8];
  x [1] = x[1] - T[7]*x[2];
  x [0] = x[0] - T[3]*x[1] - T[6]*x[2];

  DR [0] = gamma * DR[0] + theta * x[0];
  DR [1] = gamma * DR[1] + theta * x[1];
  DR [2] = gamma * DR[2] + theta * x[2];
  ACC (DR, R);

  if (kind [con] == _CONTACT_) /* projection */
  {
    float fri = mat [0],
	  coh = mat [2],
          slen, l1, l2, u1[3], u2[3], g1, g2, fri2, m [3], S [3];

    COPY (R, S);
    S [2] += coh;
    fri2 = fri*fri;
    slen = LEN2 (S);
    l1 = -(S[2] + fri*slen) / (1.0 + fri2);
    l2 =  (slen - fri*S[2]) / (1.0 + fri2);
    if (slen != 0.0)
    {
      u2[0] = S[0]/slen;
      u2[1] = S[1]/slen;
      u2[2] = -fri;
      u1[0] = -fri*u2[0];
      u1[1] = -fri*u2[1];
      u1[2] = -1.0;
    }
    else
    {
      u2[0] =  1.0;
      u2[1] =  0.0;
      u2[2] = -fri;
      u1[0] = -fri;
      u1[1] =  0.0;
      u1[2] = -1.0;
    }
    g1 = MAX (l1, 0.0);
    g2 = MAX (l2, 0.0);
    m [0] = g1*u1[0] + g2*u2[0];
    m [1] = g1*u1[1] + g2*u2[1];
    m [2] = g1*u1[2] + g2*u2[2];
    S [2] -= coh;
    SUB (S, m, R);
  }
}

/* merit function */
__global__ void merit (const int dynamic, const float step, const int num, const int *kind,
  const float *mat0, const float *A0, const float *V0, const float *U0, const float *R0, float *value)
{
  int con = blockDim.x * blockIdx.x + threadIdx.x;

  /* TODO */
}

extern "C"
{
#include "../err.h"
#include "../mem.h"
#include "../dom.h"
#include "../ldy.h"
#include "../map.h"
#include "../tmr.h"
#include "../lis.h"

#define ASSERT_CUDA(call) \
if((call) != cudaSuccess)\
{\
  cudaError_t err = cudaGetLastError(); \
  cerr << "CUDA error calling \""#call"\", code is " << err << endl; \
  THROW (ERR_CUDA);\
}

/* LOCDYN diagobal block list sorting */
#define DIABLE(i, j) ((i)->con->kind <= (j)->con->kind)
IMPLEMENT_LIST_SORT (DOUBLY_LINKED, locdyn_sort, DIAB, p, n, DIABLE)

/* PQN solver; returns the number of iterations and writes the merit function history */
int CUDA_PQN_Solve (LOCDYN *ldy, float meritval, int maxiter, float theta, float epsilon, double *merhist)
{
  float *d_R,
        *d_U,
        *d_V,
        *d_B,
        *d_A,
        *d_W,
        *d_DR,
        *d_mat;

  int *d_kind,
      *d_ptr,
      *d_adj;

  int size, num, n, *imem;
  float *fmem;
  void *mem;
  DIAB *dia;
  OFFB *blk;
  CON *con;

  /* sort constraints by kind */
  ldy->dia = locdyn_sort (ldy->dia);

  /* number constraints */
  for (dia = ldy->dia, size = num = 0; dia; dia = dia->n, num ++)
  {
    con = dia->con;
    con->num = num;
    size += dia->nadj + 1;
  }

  /* allocate device memory */
  ASSERT_CUDA (cudaMalloc((void**)&d_R, num * sizeof(float [3])));
  ASSERT_CUDA (cudaMalloc((void**)&d_U, num * sizeof(float [3])));
  ASSERT_CUDA (cudaMalloc((void**)&d_V, num * sizeof(float [3])));
  ASSERT_CUDA (cudaMalloc((void**)&d_B, num * sizeof(float [3])));
  ASSERT_CUDA (cudaMalloc((void**)&d_A, num * sizeof(float [9])));
  ASSERT_CUDA (cudaMalloc((void**)&d_W, size * sizeof(float [9])));
  ASSERT_CUDA (cudaMalloc((void**)&d_DR, num * sizeof(float [3])));
  ASSERT_CUDA (cudaMalloc((void**)&d_mat, num * sizeof(float [4])));
  ASSERT_CUDA (cudaMalloc((void**)&d_kind, num * sizeof(int)));
  ASSERT_CUDA (cudaMalloc((void**)&d_ptr, (num+1) * sizeof(int)));
  ASSERT_CUDA (cudaMalloc((void**)&d_adj, (size) * sizeof(int)));

  /* allocate host memory */
  ERRMEM (mem = malloc (size * sizeof (float [9])));

  /* copy R */
  for (dia = ldy->dia, fmem = (float*) mem; dia; dia = dia->n, fmem += 3) { double *R = dia->R; COPY (R, fmem); }
  ASSERT_CUDA (cudaMemcpy(d_R, mem, num * sizeof(float [3]), cudaMemcpyHostToDevice));
  /* copy U */
  for (dia = ldy->dia, fmem = (float*) mem; dia; dia = dia->n, fmem += 3) { double *U = dia->U; COPY (U, fmem); }
  ASSERT_CUDA (cudaMemcpy(d_U, mem, num * sizeof(float [3]), cudaMemcpyHostToDevice));
  /* copy V */
  for (dia = ldy->dia, fmem = (float*) mem; dia; dia = dia->n, fmem += 3) { double *V = dia->V; COPY (V, fmem); }
  ASSERT_CUDA (cudaMemcpy(d_V, mem, num * sizeof(float [3]), cudaMemcpyHostToDevice));
  /* copy B */
  for (dia = ldy->dia, fmem = (float*) mem; dia; dia = dia->n, fmem += 3) { double *B = dia->B; COPY (B, fmem); }
  ASSERT_CUDA (cudaMemcpy(d_B, mem, num * sizeof(float [3]), cudaMemcpyHostToDevice));
  /* copy A */
  for (dia = ldy->dia, fmem = (float*) mem; dia; dia = dia->n, fmem += 9) { double *A = dia->A; NNCOPY (A, fmem); }
  ASSERT_CUDA (cudaMemcpy(d_A, mem, num * sizeof(float [9]), cudaMemcpyHostToDevice));
  /* copy W */
  for (dia = ldy->dia, fmem = (float*) mem; dia; dia = dia->n)
  {
    double *W = dia->W;
    NNCOPY (W, fmem);
    for (blk = dia->adj, fmem += 9; blk; blk = blk->n, fmem += 9)
    {
      W = blk->W;
      NNCOPY (W, fmem); 
    }
  }
  ASSERT_CUDA (cudaMemcpy(d_W, mem, size * sizeof(float [9]), cudaMemcpyHostToDevice));
  /* zero DR */
  ASSERT_CUDA (cudaMemset(d_DR, 0, num * sizeof(float [3])));
  /* copy mat */
  for (dia = ldy->dia, fmem = (float*) mem; dia; dia = dia->n, fmem += 4)
  {
    con = dia->con;
    switch (con->kind)
    {
    case _VELODIR_:
      fmem [0] = VELODIR(con->Z);
      break;
    case _RIGLNK_:
      fmem [0] = RIGLNK_LEN(con->Z);
      break;
    case _CONTACT_:
      fmem [0] = con->mat.base->friction;
      fmem [1] = con->mat.base->restitution;
      fmem [2] = SURFACE_MATERIAL_Cohesion_Get (&con->mat) * con->area;
      fmem [3] = con->gap;
      break;
    default:
      break;
    }
  }
  ASSERT_CUDA (cudaMemcpy(d_mat, mem, num * sizeof(float [4]), cudaMemcpyHostToDevice));
  /* copy kind */
  for (dia = ldy->dia, imem = (int*) mem; dia; dia = dia->n, imem ++) { imem [0] = dia->con->kind; }
  ASSERT_CUDA (cudaMemcpy(d_kind, mem, num * sizeof(int), cudaMemcpyHostToDevice));
  /* copy ptr */
  for (dia = ldy->dia, imem = (int*) mem, imem [0] = 0; dia; dia = dia->n, imem ++) { imem [1] = dia->nadj + 1; }
  for (n = 0, imem = (int*) mem; n < num; n ++) imem [n+1] += imem [n];
  ASSERT_CUDA (cudaMemcpy(d_ptr, mem, (num+1) * sizeof(int), cudaMemcpyHostToDevice));
  /* copy adj */
  for (dia = ldy->dia, imem = (int*) mem; dia; dia = dia->n)
  {
    imem [0] = dia->con->num;
    for (blk = dia->adj, imem ++; blk; blk = blk->n, imem ++)
      imem [0] = blk->dia->con->num;
  }
  ASSERT_CUDA (cudaMemcpy(d_adj, mem, size * sizeof(int), cudaMemcpyHostToDevice));


  /* TODO */

  return 0;
}
}
