/*
 * spmv.cu
 * Copyright (C) 2010, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * sparse matrix-vector product
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
using std::cerr;
using std::endl;

#if 1
__global__ void csrspmv (
    const int rows,
    const int *ptr,
    const float *a,
    const float *x,
    float *y)
{
  int row = blockDim.x * blockIdx.x + threadIdx.x;

  if(row < rows)
  {
    float dot = 0;
    int row_start = ptr [row];
    int row_end = ptr [row+1];
    for (int j = row_start; j < row_end; j ++) dot += a [j] * x [j];
    y [row] += dot;
  }
}
#else
__global__ void csrspmv (
    const int rows,
    const int *ptr,
    const float *a,
    const float *x,
    float *y)
{
  __shared__ float vals [256];
  int thread_id = blockDim.x * blockIdx.x + threadIdx.x ; /* global thread index */
  int warp_id = thread_id / 32; /* global warp index */
  int lane = thread_id & (32 - 1); /* thread index within the warp */

  /* one warp per row */
  int row = warp_id ;

  if (row < rows)
  {
    int row_start = ptr [row];
    int row_end = ptr [row+1];

    /* compute running sum per thread  */
    vals [threadIdx.x] = 0;
    for (int j = row_start + lane ; j < row_end ; j += 32) vals [threadIdx.x] += a [j] * x [j];

    if (lane < 16) vals [threadIdx.x] += vals [threadIdx.x + 16];
    if (lane < 8) vals [threadIdx.x] += vals [threadIdx.x + 8];
    if (lane < 4) vals [threadIdx.x] += vals [threadIdx.x + 4];
    if (lane < 2) vals [threadIdx.x] += vals [threadIdx.x + 2];
    if (lane < 1) vals [threadIdx.x] += vals [threadIdx.x + 1];

    if (lane == 0) y [row] += vals [threadIdx.x];
  }
}
#endif

extern "C"
{
#include "../err.h"
#include "../mem.h"
#include "../ldy.h"
#include "../tmr.h"

#define ASSERT_CUDA(call) \
if((call) != cudaSuccess)\
{\
  cudaError_t err = cudaGetLastError(); \
  cerr << "CUDA error calling \""#call"\", code is " << err << endl; \
  THROW (ERR_CUDA);\
}

typedef struct
{
  int rows;
  int nnz;
  float *d_a;
  int *d_ptr, *h_ptr;
  float *d_x, *h_x;
  float *d_y, *h_y;
  LOCDYN *ldy;
} U_WR_B_DATA;

void* CUDA_U_WR_B_Create (LOCDYN *ldy)
{
  int i, n, m, *h_ptr;
  U_WR_B_DATA *u;
  float *h_a;
  DIAB *dia;
  OFFB *blk;

  ERRMEM (u = (U_WR_B_DATA*) MEM_CALLOC (sizeof (U_WR_B_DATA)));
  u->ldy = ldy;

  for (dia = ldy->dia, i = 0; dia; dia = dia->n, i += 3);
  u->rows = i;
 
  ERRMEM (u->h_ptr = (int*) malloc ((u->rows+1) * sizeof (int)));
  ERRMEM (u->h_y = (float*) malloc (u->rows * sizeof (float)));
  h_ptr = u->h_ptr;
  h_ptr [0] = 0;

  for (dia = ldy->dia, n = 0, i = 1; dia; dia = dia->n, i += 3)
  {
    m = (1+dia->nadj);
#if MPI
    m += dia->nadjext;
#endif
    h_ptr [i] = h_ptr [i+1] = h_ptr [i+2] = 3*m; /* nonzeros per row */
    n += m;
  }
  u->nnz = 9*n;

  ERRMEM (h_a = (float*) malloc (u->nnz * sizeof (float)));
  ERRMEM (u->h_x = (float*) malloc (u->nnz * sizeof (float)));

  for (i = 1; i <= u->rows; i ++) h_ptr [i] += h_ptr [i-1]; /* host pointers done */

  for (dia = ldy->dia, i = 0; dia; dia = dia->n, i += 3)
  {
    double *W = dia->W;
    float *a0 = &h_a  [h_ptr [i]],
          *a1 = &h_a  [h_ptr [i+1]],
          *a2 = &h_a  [h_ptr [i+2]];

    a0[0] = W[0]; a0[1] = W[3]; a0[2] = W[6];
    a1[0] = W[1]; a1[1] = W[4]; a1[2] = W[7];
    a2[0] = W[2]; a2[1] = W[5]; a2[2] = W[8];
    a0 += 3; a1 += 3; a2 += 3;
    for (blk = dia->adj; blk; blk = blk->n)
    {
      W = blk->W;
      a0[0] = W[0]; a0[1] = W[3]; a0[2] = W[6];
      a1[0] = W[1]; a1[1] = W[4]; a1[2] = W[7];
      a2[0] = W[2]; a2[1] = W[5]; a2[2] = W[8];
      a0 += 3; a1 += 3; a2 += 3;
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      W = blk->W;
      a0[0] = W[0]; a0[1] = W[3]; a0[2] = W[6];
      a1[0] = W[1]; a1[1] = W[4]; a1[2] = W[7];
      a2[0] = W[2]; a2[1] = W[5]; a2[2] = W[8];
      a0 += 3; a1 += 3; a2 += 3;
    }
#endif
  }
  
  ASSERT_CUDA (cudaMalloc((void**)&u->d_a, u->nnz * sizeof(float)));
  ASSERT_CUDA (cudaMalloc((void**)&u->d_x, u->nnz * sizeof(float)));
  ASSERT_CUDA (cudaMalloc((void**)&u->d_y, u->rows * sizeof(float)));
  ASSERT_CUDA (cudaMalloc((void**)&u->d_ptr, (u->rows+1) * sizeof(int)));

  /* copy h_a  to u->d_a */
  ASSERT_CUDA (cudaMemcpy(u->d_a, h_a, u->nnz * sizeof(float), cudaMemcpyHostToDevice));

  /* copy h_ptr  to u->d_ptr */
  ASSERT_CUDA (cudaMemcpy(u->d_ptr, h_ptr, (u->rows+1) * sizeof(int), cudaMemcpyHostToDevice));

  free (h_a);

  return u;
}

void CUDA_U_WR_B (void *U_WR_B)
{
  U_WR_B_DATA *u = (U_WR_B_DATA*) U_WR_B;
  float *h_x, *h_y;
  int *h_ptr;
  DIAB *dia;
  OFFB *blk;
  TIMING tt;
  double t[2];

  timerstart (&tt);

  /* copy reactions */
  for (dia = u->ldy->dia, h_ptr = u->h_ptr, h_x = u->h_x; dia; dia = dia->n, h_ptr += 3)
  {
    double *R = dia->R;
    float *a0 = &h_x  [h_ptr [0]],
          *a1 = &h_x  [h_ptr [1]],
          *a2 = &h_x  [h_ptr [2]];

    a0[0] = a1[0] = a2 [0] = R[0];
    a0[1] = a1[1] = a2 [1] = R[1];
    a0[2] = a1[2] = a2 [2] = R[2];
    a0 += 3; a1 += 3; a2 += 3;
    for (blk = dia->adj; blk; blk = blk->n)
    {
      R = blk->dia->R;
      a0[0] = a1[0] = a2 [0] = R[0];
      a0[1] = a1[1] = a2 [1] = R[1];
      a0[2] = a1[2] = a2 [2] = R[2];
      a0 += 3; a1 += 3; a2 += 3;
    }
#if MPI
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      R = CON(blk->dia)->R;
      a0[0] = a1[0] = a2 [0] = R[0];
      a0[1] = a1[1] = a2 [1] = R[1];
      a0[2] = a1[2] = a2 [2] = R[2];
      a0 += 3; a1 += 3; a2 += 3;
    }
#endif
  }

  /* copy u->h_x  to u->d_x */
  ASSERT_CUDA (cudaMemcpy(u->d_x, u->h_x, u->nnz * sizeof(float), cudaMemcpyHostToDevice));

  /* zero u->d_y */
  ASSERT_CUDA (cudaMemset(u->d_y, 0, u->rows * sizeof(float)));

  t[0] = timerend (&tt);
  timerstart (&tt);

  /* GPU matrix vector product */  
  int tpb = 256;
  int bpg = (u->rows + tpb - 1) / tpb;
  csrspmv <<<bpg, tpb>>> (u->rows, u->d_ptr, u->d_a, u->d_x, u->d_y);

  t[1] = timerend (&tt);
  timerstart (&tt);

  /* copy u->d_y to u->h_y */
  ASSERT_CUDA (cudaMemcpy(u->h_y, u->d_y, u->rows * sizeof(float), cudaMemcpyDeviceToHost));

  /* update U */
  for (dia = u->ldy->dia, h_y = u->h_y; dia; dia = dia->n, h_y += 3)
  {
    double *B = dia->B,
           *U = dia->U;

    U [0] = B[0] + h_y[0];
    U [1] = B[1] + h_y[1];
    U [2] = B[2] + h_y[2];
  }

  t[0] += timerend (&tt);
  printf ("CUDA TIMING: host %g, devi %g\n", t[0], t[1]);
}

void CUDA_U_WR_B_Destroy (void *U_WR_B)
{
  U_WR_B_DATA *u = (U_WR_B_DATA*) U_WR_B;

  ASSERT_CUDA (cudaFree(u->d_a));
  ASSERT_CUDA (cudaFree(u->d_x));
  ASSERT_CUDA (cudaFree(u->d_y));
  ASSERT_CUDA (cudaFree(u->d_ptr));
  free (u->h_ptr);
  free (u->h_x);
  free (u->h_y);
  free (u);
}
}
