/*
 * fem.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * finite element method
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

#include <float.h>
#include "mem.h"
#include "dom.h"
#include "fem.h"
#include "but.h"
#include "alg.h"
#include "bla.h"
#include "hyb.h"
#include "cvi.h"
#include "svk.h"
#include "gjk.h"
#include "err.h"

typedef double (*node_t) [3]; /* mesh node */

#define DOM_TOL 0.150
#define CUT_TOL 0.015
#define MAX_NODES_COUNT 64
#define FEM_VEL0(bod) ((bod)->velo + (bod)->dofs)
#define FEM_FORCE(bod) ((bod)->velo + (bod)->dofs * 2)
#define FEM_FEXT(bod) ((bod)->velo + (bod)->dofs * 3)
#define FEM_FINT(bod) ((bod)->velo + (bod)->dofs * 4)
#define FEM_MESH(bod) ((bod)->msh ? (bod)->msh : (bod)->shape->data)
#define FEM_MATERIAL(bod, ele) ((ele)->mat ? (ele)->mat : (bod)->mat)

static const double I_TET1_X[] = {0.25};
static const double I_TET1_Y[] = {0.25};
static const double I_TET1_Z[] = {0.25};
static const double I_TET1_W[] = {0.16666666666666666};
static const int    I_TET1_N   =  1;

static const double I_TET2_X [] = {0.13819660112501052, 0.13819660112501052, 0.13819660112501052, 0.58541019662496840};
static const double I_TET2_Y [] = {0.13819660112501052, 0.13819660112501052, 0.58541019662496840, 0.13819660112501052};
static const double I_TET2_Z [] = {0.13819660112501052, 0.58541019662496840, 0.13819660112501052, 0.13819660112501052};
static const double I_TET2_W [] = {0.04166666666666666, 0.04166666666666666, 0.04166666666666666, 0.04166666666666666};
static const int    I_TET2_N    =  4;

#define ISQR3 0.57735026918962584 
static const double I_HEX2_X [] = {-ISQR3, ISQR3, ISQR3, -ISQR3, -ISQR3, ISQR3, ISQR3, -ISQR3};
static const double I_HEX2_Y [] = {-ISQR3, -ISQR3, ISQR3, ISQR3, -ISQR3, -ISQR3, ISQR3, ISQR3};
static const double I_HEX2_Z [] = {-ISQR3, -ISQR3, -ISQR3, -ISQR3, ISQR3, ISQR3, ISQR3, ISQR3};
static const double I_HEX2_W [] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
static const int    I_HEX2_N    = 8;

/* linear tetrahedron shape functions */
inline static void tet_o1_shapes (double *point, double *shapes)
{
  shapes [0] = point [0];
  shapes [1] = point [1];
  shapes [2] = point [2];
  shapes [3] = 1.0 - (point [0] + point [1] + point [2]);
}

/* linear hexahedron shape functions */
inline static void hex_o1_shapes (double *point, double *shapes)
{
  shapes [0] = 0.125 * (1.0 - point [0]) * (1.0 - point [1]) * (1.0 - point [2]);
  shapes [1] = 0.125 * (1.0 + point [0]) * (1.0 - point [1]) * (1.0 - point [2]);
  shapes [2] = 0.125 * (1.0 + point [0]) * (1.0 + point [1]) * (1.0 - point [2]);
  shapes [3] = 0.125 * (1.0 - point [0]) * (1.0 + point [1]) * (1.0 - point [2]);
  shapes [4] = 0.125 * (1.0 - point [0]) * (1.0 - point [1]) * (1.0 + point [2]);
  shapes [5] = 0.125 * (1.0 + point [0]) * (1.0 - point [1]) * (1.0 + point [2]);
  shapes [6] = 0.125 * (1.0 + point [0]) * (1.0 + point [1]) * (1.0 + point [2]);
  shapes [7] = 0.125 * (1.0 - point [0]) * (1.0 + point [1]) * (1.0 + point [2]);
}

/* linear tetrahedron shape functions */
inline static void tet_o1_derivs (double *point, double *derivs)
{
  derivs [0] = derivs [4] = derivs [8] = 1.0;
  derivs [1] = derivs [2] = derivs [3] = derivs [5] = derivs [6] = derivs [7] = 0.0;
  derivs [9] = derivs [10] = derivs [11] = -1.0;
}

/* linear hexahedron shape functions */
inline static void hex_o1_derivs (double *point, double *derivs)
{
  derivs[0] = -0.125 * (1 - point[1]) * (1 - point[2]);
  derivs[1] = -0.125 * (1 - point[0]) * (1 - point[2]);
  derivs[2] = -0.125 * (1 - point[0]) * (1 - point[1]);

  derivs[3] =  0.125 * (1 - point[1]) * (1 - point[2]);
  derivs[4] = -0.125 * (1 + point[0]) * (1 - point[2]);
  derivs[5] = -0.125 * (1 + point[0]) * (1 - point[1]);

  derivs[6] =  0.125 * (1 + point[1]) * (1 - point[2]);
  derivs[7] =  0.125 * (1 + point[0]) * (1 - point[2]);
  derivs[8] = -0.125 * (1 + point[0]) * (1 + point[1]);

  derivs[9] = -0.125 * (1 + point[1]) * (1 - point[2]);
  derivs[10] =  0.125 * (1 - point[0]) * (1 - point[2]);
  derivs[11] = -0.125 * (1 - point[0]) * (1 + point[1]);

  derivs[12] = -0.125 * (1 - point[1]) * (1 + point[2]);
  derivs[13] = -0.125 * (1 - point[0]) * (1 + point[2]);
  derivs[14] =  0.125 * (1 - point[0]) * (1 - point[1]);

  derivs[15] =  0.125 * (1 - point[1]) * (1 + point[2]);
  derivs[16] = -0.125 * (1 + point[0]) * (1 + point[2]);
  derivs[17] =  0.125 * (1 + point[0]) * (1 - point[1]);

  derivs[18] =  0.125 * (1 + point[1]) * (1 + point[2]);
  derivs[19] =  0.125 * (1 + point[0]) * (1 + point[2]);
  derivs[20] =  0.125 * (1 + point[0]) * (1 + point[1]);

  derivs[21] = -0.125 * (1 + point[1]) * (1 + point[2]);
  derivs[22] =  0.125 * (1 - point[0]) * (1 + point[2]);
  derivs[23] =  0.125 * (1 - point[0]) * (1 + point[1]);  
}

/* linear tetrahedron local to global point transformation */
inline static void tet_local_to_global (node_t nodes, double *local, double *global)
{
  double shapes [4];
  int i, j;

  tet_o1_shapes (local, shapes);

  SET (global, 0.0);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 4; j ++)
      global [i] += nodes [j][i] * shapes [j];
}

/* linear tetrahedron global to local point transformation */
inline static void tet_global_to_local (node_t nodes, double *global, double *local)
{
  double A [9], B [3], I [9], det;

  SUB (nodes [0], nodes [3], A);
  SUB (nodes [1], nodes [3], A+3);
  SUB (nodes [2], nodes [3], A+6);
  SUB (global, nodes [3], B);
  INVERT (A, I, det);
  ASSERT (det > 0.0, ERR_FEM_COORDS_INVERT);
  NVMUL (I, B, local);

#if 0
  ASSERT_DEBUG (local [0] >= 0.0 && local [1] >= 0.0 local [2] >= 0.0 &&
  (local [0] + local [1] + local [2]) <= 1.0, "Local coords out of bounds");
#endif
}

/* linear hexahedron local to global point transformation */
inline static void hex_local_to_global (node_t nodes, double *local, double *global)
{
  double shapes [8];
  int i, j;

  hex_o1_shapes (local, shapes);

  SET (global, 0.0);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 8; j ++)
      global [i] += nodes [j][i] * shapes [j];
}

/* linear hexahedron global to local point transformation */
inline static void hex_global_to_local (node_t nodes, double *global, double *local)
{
  double A [9], B [3], I [9], det, error;
  double shapes [8], derivs [24];
  int i, j, k, l;

  SET (local, 0.0);
  l = 0;

  do
  {
    hex_o1_shapes (local, shapes);

    COPY (global, B);

    for (i = 0; i < 3; i ++)
      for (j = 0; j < 8; j ++)
	B [i] -= nodes [j][i] * shapes [j];

    hex_o1_derivs (local, derivs);

    SET9 (A, 0.0);

    for (i = 0; i < 3; i ++)
      for (j = 0; j < 3; j ++)
	for (k = 0; k < 8; k ++) A [3*j+i] += nodes[k][i] * derivs [3*k+j];

    INVERT (A, I, det);
    ASSERT (det > 0.0, ERR_FEM_COORDS_INVERT);
    NVMUL (I, B, A);
    ADD (local, A, local);
    error = sqrt (DOT (A,A) / (1.0 + DOT (local, local)));

  } while (++l < 64 && error > 1E-9);

  ASSERT (l < 64, ERR_FEM_COORDS_INVERT);

#if 0
  ASSERT_DEBUG (local [0] >= 0.0 && local [0] <= 1.0 &&
                local [1] >= 0.0 && local [1] <= 1.0 &&
		local [2] >= 0.0 && local [2] <= 1.0, "Local coords out of bounds");
#endif
}

/* linear tetrahedron transformation determinant at local point */
inline static double tet_o1_det (node_t nodes, double *point, double *F)
{
  double derivs [12], G [9];
  int i, j, k;

  tet_o1_derivs (point, derivs);

  if (!F) F = G;

  SET9 (F, 0.0);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 3; j ++)
      for (k = 0; k < 4; k ++) F [3*j+i] += nodes[k][i] * derivs [3*k+j];

  return DET (F);
}

/* linear hexahedron transformation determinant at local point */
inline static double hex_o1_det (node_t nodes, double *point, double *F)
{
  double derivs [24], G [9];
  int i, j, k;

  hex_o1_derivs (point, derivs);

  if (!F) F = G;

  SET9 (F, 0.0);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 3; j ++)
      for (k = 0; k < 8; k ++) F [3*j+i] += nodes[k][i] * derivs [3*k+j];

  return DET (F);
}

/* linear tetrahedron deformation determinant at local point */
inline static void tet_o1_gradient (node_t q, double *point, double *F0, double *derivs, double *F)
{
  double local_derivs [12], IF0 [9], det, *l, *d;
  int i, j, k;

  tet_o1_derivs (point, local_derivs);

  INVERT (F0, IF0, det);
  ASSERT (det > 0.0, ERR_FEM_COORDS_INVERT);
  for (k = 0, l = local_derivs, d = derivs; k < 4; k ++, l += 3, d += 3) { TVMUL (IF0, l, d); }

  IDENTITY (F);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 3; j ++)
      for (k = 0; k < 4; k ++) F [3*j+i] += q[k][i] * derivs [3*k+j];
}

/* linear hexahedron deformation determinant at local point */
inline static void hex_o1_gradient (node_t q, double *point, double *F0, double *derivs, double *F)
{
  double local_derivs [24], IF0 [9], det, *l, *d;
  int i, j, k;

  hex_o1_derivs (point, local_derivs);

  INVERT (F0, IF0, det);
  ASSERT (det > 0.0, ERR_FEM_COORDS_INVERT);
  for (k = 0, l = local_derivs, d = derivs; k < 8; k ++, l += 3, d += 3) { TVMUL (IF0, l, d); }

  IDENTITY (F);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 3; j ++)
      for (k = 0; k < 8; k ++) F [3*j+i] += q[k][i] * derivs [3*k+j];
}

/* integration macro */
#define INTEGRATE(SCHEME, SUBORD, dom, domnum, nodes, ...)\
{\
  double point [3], weight;\
  int __k__, __l__;\
\
  if (dom)\
  {\
    double __subnodes__ [4][3], __subJ__;\
    double __subpoint__ [3];\
    TRI *__t__, *__e__;\
    double __volume__;\
\
    for (__volume__ = __k__ = 0; __k__ < domnum; __k__ ++) __volume__ += dom [__k__].volume;\
\
    for (__l__ = 0; __l__ < domnum; dom ++, __l__ ++)\
    {\
      if (dom->volume > DOM_TOL * __volume__)\
      {\
	COPY (dom->center, __subnodes__ [3]);\
\
	for (__t__ = dom->tri, __e__ = __t__ + dom->m; __t__ < __e__; __t__ ++)\
	{\
	  COPY (__t__->ver [0], __subnodes__ [0]);\
	  COPY (__t__->ver [1], __subnodes__ [1]);\
	  COPY (__t__->ver [2], __subnodes__ [2]);\
\
	  for (__k__ = 0; __k__ < I_TET##SUBORD##_N; __k__ ++)\
	  {\
	    __subpoint__ [0] = I_TET##SUBORD##_X [__k__];\
	    __subpoint__ [1] = I_TET##SUBORD##_Y [__k__];\
	    __subpoint__ [2] = I_TET##SUBORD##_Z [__k__];\
	    __subJ__ = tet_o1_det (__subnodes__, __subpoint__, NULL);\
	    tet_local_to_global (__subnodes__, __subpoint__, point);\
	    weight = __subJ__ * I_TET##SUBORD##_W [__k__];\
\
	    __VA_ARGS__\
	  }\
	}\
      }\
      else\
      {\
	COPY (dom->center, point);\
	weight = dom->volume;\
\
	__VA_ARGS__\
      }\
    }\
  }\
  else\
  {\
    for (__k__ = 0; __k__ < SCHEME##_N; __k__ ++)\
    {\
      point [0] = SCHEME##_X [__k__];\
      point [1] = SCHEME##_Y [__k__];\
      point [2] = SCHEME##_Z [__k__];\
      weight = SCHEME##_W [__k__];\
\
      __VA_ARGS__\
    }\
  }\
}

/* lump linear tetrahedron mass */
inline static void tet_o1_lump (TRISURF *dom, int domnum, node_t nodes, double density, double **out)
{
  double J, coef, integral;
  double shapes [4];
  int i, j;

  INTEGRATE (I_TET2, 2, dom, domnum, nodes,

    tet_o1_shapes (point, shapes);
    J = tet_o1_det (nodes, point, NULL);
    coef = density * J * weight;

    for (i = 0; i < 4; i ++)
    {
      for (j = 0; j < 4; j ++)
      {
	integral = coef * shapes [i] * shapes [j];

	out [i][0] += integral;
	out [i][1] += integral;
	out [i][2] += integral;
      }
    }
  )
}

/* lump linear hexahedron mass */
inline static void hex_o1_lump (TRISURF *dom, int domnum, node_t nodes, double density, double **out)
{
  double J, coef, integral;
  double shapes [8];
  int i, j;

  INTEGRATE (I_HEX2, 2, dom, domnum, nodes,

    hex_o1_shapes (point, shapes);
    J = hex_o1_det (nodes, point, NULL);
    coef = density * J * weight;

    for (i = 0; i < 8; i ++)
    {
      for (j = 0; j < 8; j ++)
      {
	integral = coef * shapes [i] * shapes [j];

	out [i][0] += integral;
	out [i][1] += integral;
	out [i][2] += integral;
      }
    }
  )
}

/* compute linear tetrahedron body force contribution */
inline static void tet_o1_body_force (TRISURF *dom, int domnum, node_t nodes, double density, double *f, double *g)
{
  double J, coef, integral;
  double shapes [4];
  int i;

  for (i = 0; i < 12; i ++) g [i] = 0.0;

  INTEGRATE (I_TET1, 1, dom, domnum, nodes,

    tet_o1_shapes (point, shapes);
    J = tet_o1_det (nodes, point, NULL);
    coef = density * J * weight;

    for (i = 0; i < 4; i ++)
    {
      integral = coef * shapes [i];

      g [3*i+0] += f [0] * integral;
      g [3*i+1] += f [1] * integral;
      g [3*i+2] += f [2] * integral;
    }
  )
}

/* compute linear hexahedron body force contribution */
inline static void hex_o1_body_force (TRISURF *dom, int domnum, node_t nodes, double density, double *f, double *g)
{
  double J, coef, integral;
  double shapes [8];
  int i;

  for (i = 0; i < 24; i ++) g [i] = 0.0;

  INTEGRATE (I_HEX2, 1, dom, domnum, nodes,

    hex_o1_shapes (point, shapes);
    J = hex_o1_det (nodes, point, NULL);
    coef = density * J * weight;

    for (i = 0; i < 8; i ++)
    {
      integral = coef * shapes [i];

      g [3*i+0] += f [0] * integral;
      g [3*i+1] += f [1] * integral;
      g [3*i+2] += f [2] * integral;
    }
  )
}

/* compute linear tetrahedron internal force or force derivative contribution */
inline static void tet_o1_internal_force (short derivative, TRISURF *dom, int domnum, node_t nodes, BULK_MATERIAL *mat, double (*q) [3], double *g)
{
  double derivs [12], F0 [9], F [9], P [9], K [81], KB [9], *B, *p;
  double J, integral, mat_lambda, mat_mi;
  int i, j;

  for (i = 0, j = 12 * (derivative ? 12 : 1); i < j; i ++) g [i] = 0.0;
  mat_lambda = lambda (mat->young, mat->poisson);
  mat_mi  = mi (mat->young, mat->poisson);

  INTEGRATE (I_TET1, 1, dom, domnum, nodes,

    J = tet_o1_det (nodes, point, F0);
    tet_o1_gradient (q, point, F0, derivs, F);
    integral = J * weight;

    if (derivative)
    {
      SVK_Tangent_C (mat_lambda, mat_mi, integral, 9, F, K);

      for (i = 0; i < 12; i ++) /* see doc/notes.lyx for details */
      {
	SET9 (KB, 0);
	for (j = 0; j < 3; j ++)
	{
	  p = &K [9*((i%3) + (3*j))];
	  integral = derivs [3*(i/3)+j];
	  NNADDMUL (KB, integral, p, KB);
	}

	for (j = 0, B = derivs, p = &g[12*i]; j < 4; j ++, B += 3, p += 3) { NVADDMUL (p, KB, B, p); }
      }
    }
    else
    {
      SVK_Stress_C (mat_lambda, mat_mi, integral, F, P); /* column-wise */

      for (i = 0, B = derivs, p = g; i < 4; i ++, B += 3, p += 3) { NVADDMUL (p, P, B, p); }
    }
  )
}

/* compute linear hexahedron internal force or force derivative contribution */
inline static void hex_o1_internal_force (short derivative, TRISURF *dom, int domnum, node_t nodes, BULK_MATERIAL *mat, double (*q) [3], double *g)
{
  double derivs [24], F0 [9], F [9], P [9], K [81], KB [9], *B, *p;
  double J, integral, mat_lambda, mat_mi;
  int i, j;

  for (i = 0, j = 24 * (derivative ? 24 : 1); i < j; i ++) g [i] = 0.0;
  mat_lambda = lambda (mat->young, mat->poisson);
  mat_mi  = mi (mat->young, mat->poisson);

  INTEGRATE (I_HEX2, 1, dom, domnum, nodes,

    J = hex_o1_det (nodes, point, F0);
    hex_o1_gradient (q, point, F0, derivs, F);
    integral = J * weight;

    if (derivative)
    {
      SVK_Tangent_C (mat_lambda, mat_mi, integral, 9, F, K);

      for (i = 0; i < 24; i ++) /* see doc/notes.lyx for details */
      {
	SET9 (KB, 0);
	for (j = 0; j < 3; j ++)
	{
	  p = &K [9*((i%3) + (3*j))];
	  integral = derivs [3*(i/3)+j];
	  NNADDMUL (KB, integral, p, KB);
	}

	for (j = 0, B = derivs, p = &g[24*i]; j < 8; j ++, B += 3, p += 3) { NVADDMUL (p, KB, B, p); }
      }
    }
    else
    {
      SVK_Stress_C (mat_lambda, mat_mi, integral, F, P); /* column-wise */

      for (i = 0, B = derivs, p = g; i < 8; i ++, B += 3, p += 3) { NVADDMUL (p, P, B, p); }
    }
  )
}

/* copy node coordinates into a local table */
static void load_nodes (node_t heap, int type, int *nodes, node_t stack)
{
  int n;

  for (n = 0; n < type; n ++)
  { COPY (heap [nodes [n]], stack [n]); }
}

/* lump contribution of the element mass into the diagonal */
static void lump_mass_matrix (BODY *bod, MESH *msh, ELEMENT *ele, double *x)
{
  double nodes [8][3],
	 density,
	*out [8];
  int i;

  density = ele->mat ? ele->mat->density : bod->mat->density;

  load_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  for (i = 0; i < ele->type; i ++) out [i] = &x [3 * ele->nodes [i]];
 
  switch (bod->form)
  {
  case FEM_O1:
  {
    switch (ele->type)
    {
    case 4: tet_o1_lump (ele->dom, ele->domnum, nodes, density, out); break;
    case 8: hex_o1_lump (ele->dom, ele->domnum, nodes, density, out); break;
    case 5:
      COPY (nodes [4], nodes [5]); out [5] = out [4];
      COPY (nodes [4], nodes [6]); out [6] = out [4];
      COPY (nodes [4], nodes [7]); out [7] = out [4];
      hex_o1_lump (ele->dom, ele->domnum, nodes, density, out);
      break;
    case 6:
      COPY (nodes [5], nodes [7]); out [7] = out [5];
      COPY (nodes [5], nodes [6]); out [6] = out [5];
      COPY (nodes [4], nodes [5]); out [5] = out [4];
      COPY (nodes [3], nodes [4]); out [4] = out [3];
      COPY (nodes [2], nodes [3]); out [3] = out [2];
      hex_o1_lump (ele->dom, ele->domnum, nodes, density, out);
      break;
    }
  }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }
}

/* compute deformation gradient at a local point */
static void deformation_gradient (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *F)
{
  double derivs [24], nodes [8][3], q [8][3], F0 [9], *p;
  int i;

  load_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  for (i = 0; i < ele->type; i ++)
  {
    p = &bod->conf [3 * ele->nodes [i]];
    COPY (p, q[i]);
  }

  switch (bod->form)
  {
  case FEM_O1:
  {
    switch (ele->type)
    {
    case 4: 
      tet_o1_det (nodes, point, F0);
      tet_o1_gradient (q, point, F0, derivs, F);
      break;
    case 8:
      hex_o1_det (nodes, point, F0);
      hex_o1_gradient (q, point, F0, derivs, F);
      break;
    case 5:
      COPY (nodes [4], nodes [5]);
      COPY (nodes [4], nodes [6]);
      COPY (nodes [4], nodes [7]);
      COPY (q [4], q [5]);
      COPY (q [4], q [6]);
      COPY (q [4], q [7]);
      hex_o1_det (nodes, point, F0);
      hex_o1_gradient (q, point, F0, derivs, F);
      break;
    case 6:
      COPY (nodes [5], nodes [7]);
      COPY (nodes [5], nodes [6]);
      COPY (nodes [4], nodes [5]);
      COPY (nodes [3], nodes [4]);
      COPY (nodes [2], nodes [3]);
      COPY (q [5], q [7]);
      COPY (q [5], q [6]);
      COPY (q [4], q [5]);
      COPY (q [3], q [4]);
      COPY (q [2], q [3]);
      hex_o1_det (nodes, point, F0);
      hex_o1_gradient (q, point, F0, derivs, F);
      break;
    }
  }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }
}

/* compute element shape functions at a local point and return global matrix */
static MX* shape_functions (BODY *bod, MESH *msh, ELEMENT *ele, double *point)
{
  MX *N = NULL;

  switch (bod->form)
  {
  case FEM_O1:
  {
    static int *p, *i, *q, *u, k, n, m;
    double shapes [8], sum, *x, *y;

    n = bod->dofs;

    ERRMEM (p = MEM_CALLOC ((2*n+1) * sizeof (int)));
    i = p + n + 1;

    for (k = 0, u = i; k < ele->type; k ++, u += 3)
    {
      m = ele->nodes [k] * 3;
      q = &p [m + 1];
      q [0] = 1; q [1] = 1; q [2] = 1;
      u [0] = 0; u [1] = 1; u [2] = 2;
    }

    for (q = p + 1; q < i; q ++) *q += *(q-1);

    N = MX_Create (MXCSC, 3, n, p, i);
    x = N->x;
    free (p);
    p = N->p;

    switch (ele->type)
    {
      case 4:
	tet_o1_shapes (point, shapes);
	for (k = 0; k < 4; k ++)
	{
	  m = ele->nodes [k] * 3;
	  y = &x [p [m]];
	  SET (y, shapes [k]);
	}
	break;
      case 8:
	hex_o1_shapes (point, shapes);
	for (k = 0; k < 8; k ++)
	{
	  m = ele->nodes [k] * 3;
	  y = &x [p [m]];
	  SET (y, shapes [k]);
	}
	break;
      case 5:
	hex_o1_shapes (point, shapes);
	for (k = 0; k < 4; k ++)
	{
	  m = ele->nodes [k] * 3;
	  y = &x [p [m]];
	  SET (y, shapes [k]);
	}
	m = ele->nodes [k] * 3;
	y = &x [p [m]];
	sum = (shapes [4] + shapes [5] + shapes [6] + shapes [7]);
	SET (y, sum);
	break;
      case 6:
	hex_o1_shapes (point, shapes);
	for (k = 0; k < 2; k ++)
	{
	  m = ele->nodes [k] * 3;
	  y = &x [p [m]];
	  SET (y, shapes [k]);
	}
	m = ele->nodes [2] * 3;
	y = &x [p [m]];
	sum = (shapes [2] + shapes [3]);
	SET (y, sum);
	for (k = 4; k < 6; k ++)
	{
	  m = ele->nodes [k-1] * 3;
	  y = &x [p [m]];
	  SET (y, shapes [k]);
	}
	m = ele->nodes [5] * 3;
	y = &x [p [m]];
	sum = (shapes [6] + shapes [7]);
	SET (y, sum);
	break;
    }
  }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }

  return N;
}

/* compute element shape functions at a local point and output them into a local matrix */
static int shape_functions_ext (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *shapes)
{

  switch (bod->form)
  {
  case FEM_O1:
  {
    switch (ele->type)
    {
      case 4: tet_o1_shapes (point, shapes); return 4;
      case 8: hex_o1_shapes (point, shapes); return 8;
      case 5: hex_o1_shapes (point, shapes); shapes [4] += shapes [5] + shapes [6] + shapes [7]; return 5;
      case 6: hex_o1_shapes (point, shapes); shapes [2] += shapes [3];  shapes [3] = shapes [4];
	                          shapes [4] = shapes [5]; shapes [5] = shapes [6] + shapes [7]; return 6;
    }
  }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }

  return 0;
}

/* load displacemnts into a local buffer */
static void load_displacements (BODY *bod, MESH *msh, ELEMENT *ele, double (*q) [3])
{
  double *conf = bod->conf, *p;
  int i;

  switch (bod->form)
  {
  case FEM_O1:
  {
    for (i = 0; i < ele->type; i ++)
    {
      p = &conf [3 * ele->nodes [i]];
      COPY (p, q[i]);
    }
  }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }
}

/* load velocities into a local buffers */
static void load_velocities (BODY *bod, MESH *msh, ELEMENT *ele, double (*u0) [3], double (*u) [3])
{
  double *vel0 = FEM_VEL0 (bod), *velo = bod->velo, *p;
  int i, j;

  switch (bod->form)
  {
  case FEM_O1:
  {
    for (i = 0; i < ele->type; i ++)
    {
      j = 3 * ele->nodes [i];
      p = &vel0 [j];
      COPY (p, u0[i]);
      p = &velo [j];
      COPY (p, u[i]);
    }
  }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }
}

/* load node number into local buffer */
static int load_node_numbers (BODY *bod, MESH *msh, ELEMENT *ele, int *numbers)
{
  int i;

  switch (bod->form)
  {
  case FEM_O1:
  {
    for (i = 0; i < ele->type; i ++) numbers [i] = ele->nodes [i];

    return i;
  }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }

  return 0;
}
 
/* copute point force contribution */
static void point_force (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *f, double *force)
{
  MX *N = shape_functions (bod, msh, ele, point);
  MX_Matvec (1.0, MX_Tran (N), f, 1.0, force);
  MX_Destroy (N);
}

/* copute body force contribution */
static void body_force (BODY *bod, MESH *msh, ELEMENT *ele, double *f, double *g)
{
  double nodes [8][3],
	 density;

  density = ele->mat ? ele->mat->density : bod->mat->density;

  load_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  switch (bod->form)
  {
  case FEM_O1:
    switch (ele->type)
    {
    case 4: tet_o1_body_force (ele->dom, ele->domnum, nodes, density, f, g); break;
    case 8: hex_o1_body_force (ele->dom, ele->domnum, nodes, density, f, g); break;
    case 5:
      COPY (nodes [4], nodes [5]);
      COPY (nodes [4], nodes [6]);
      COPY (nodes [4], nodes [7]);
      hex_o1_body_force (ele->dom, ele->domnum, nodes, density, f, g);
      break;
    case 6:
      COPY (nodes [5], nodes [7]);
      COPY (nodes [5], nodes [6]);
      COPY (nodes [4], nodes [5]);
      COPY (nodes [3], nodes [4]);
      COPY (nodes [2], nodes [3]);
      hex_o1_body_force (ele->dom, ele->domnum, nodes, density, f, g);
      break;
    }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }
}

/* copute internal force or force derivative contribution */
static void internal_force (int derivative, BODY *bod, MESH *msh, ELEMENT *ele, double *g)
{
  BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
  double nodes [8][3], q [8][3], *p;
  int i;

  load_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  for (i = 0; i < ele->type; i ++)
  {
    p = &bod->conf [3 * ele->nodes [i]];
    COPY (p, q[i]);
  }

  switch (bod->form)
  {
  case FEM_O1:
    switch (ele->type)
    {
    case 4: tet_o1_internal_force (derivative, ele->dom, ele->domnum, nodes, mat, q, g); break;
    case 8: hex_o1_internal_force (derivative, ele->dom, ele->domnum, nodes, mat, q, g); break;
    case 5:
      COPY (nodes [4], nodes [5]);
      COPY (nodes [4], nodes [6]);
      COPY (nodes [4], nodes [7]);
      COPY (q [4], q [5]);
      COPY (q [4], q [6]);
      COPY (q [4], q [7]);
      hex_o1_internal_force (derivative, ele->dom, ele->domnum, nodes, mat, q, g);
      break;
    case 6:
      COPY (nodes [5], nodes [7]);
      COPY (nodes [5], nodes [6]);
      COPY (nodes [4], nodes [5]);
      COPY (nodes [3], nodes [4]);
      COPY (nodes [2], nodes [3]);
      COPY (q [5], q [7]);
      COPY (q [5], q [6]);
      COPY (q [4], q [5]);
      COPY (q [3], q [4]);
      COPY (q [2], q [3]);
      hex_o1_internal_force (derivative, ele->dom, ele->domnum, nodes, mat, q, g);
      break;
    }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }
}

/* compute local coordinates of a global point */
static void global_to_local (double (*mesh_nodes) [3], ELEMENT *ele, double *point, double *local)
{
  double nodes [8][3];

  load_nodes (mesh_nodes, ele->type, ele->nodes, nodes);

  switch (ele->type)
  {
    case 4: tet_global_to_local (nodes, point, local); break;
    case 8: hex_global_to_local (nodes, point, local); break;
    case 5:
    {
      COPY (nodes [4], nodes [5]);
      COPY (nodes [4], nodes [6]);
      COPY (nodes [4], nodes [7]);
      hex_global_to_local (nodes, point, local);
    }
    break;
    case 6:
    {
      COPY (nodes [5], nodes [7]);
      COPY (nodes [5], nodes [6]);
      COPY (nodes [4], nodes [5]);
      COPY (nodes [3], nodes [4]);
      COPY (nodes [2], nodes [3]);
      hex_global_to_local (nodes, point, local);
    }
    break;
  }
}

/* compute local coordinates of a spatial point */
#define spatial_to_local(msh, ele, x, point) global_to_local ((msh)->cur_nodes, ele, x, point)

/* compute local coordinates of a referential point */
#define referential_to_local(msh, ele, X, point) global_to_local ((msh)->ref_nodes, ele, X, point)

/* copute global coordinates of a local point */
static void local_to_global (double (*mesh_nodes) [3], ELEMENT *ele, double *point, double *X)
{
  double nodes [8][3];

  load_nodes (mesh_nodes, ele->type, ele->nodes, nodes);

  switch (ele->type)
  {
    case 4: tet_local_to_global (nodes, point, X); break;
    case 8: hex_local_to_global (nodes, point, X); break;
    case 5:
    {
      COPY (nodes [4], nodes [5]);
      COPY (nodes [4], nodes [6]);
      COPY (nodes [4], nodes [7]);
      hex_local_to_global (nodes, point, X); break;
    }
    break;
    case 6:
    {
      COPY (nodes [5], nodes [7]);
      COPY (nodes [5], nodes [6]);
      COPY (nodes [4], nodes [5]);
      COPY (nodes [3], nodes [4]);
      COPY (nodes [2], nodes [3]);
      hex_local_to_global (nodes, point, X); break;
    }
    break;
  }
}

/* compute spatial coordinates of a local point */
#define local_to_spatial(msh, ele, x, point) local_to_global ((msh)->cur_nodes, ele, x, point)

/* compute referential coordinates of a local point */
#define local_to_referential(msh, ele, X, point) local_to_global ((msh)->ref_nodes, ele, X, point)

/* compute Cauchy stress at local point */
static void fem_element_cauchy (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *values)
{
  BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
  double P [9], J, F [9];

  deformation_gradient (bod, msh, ele, point, F);

  J = SVK_Stress_C (lambda (mat->young, mat->poisson), mi (mat->young, mat->poisson), 1.0, F, P); /* column-wise, per unit volume */

  values [0] = (F[0]*P[0]+F[3]*P[1]+F[6]*P[2])/J; /* sx  */
  values [1] = (F[1]*P[3]+F[4]*P[4]+F[7]*P[5])/J; /* sy  */
  values [2] = (F[2]*P[6]+F[5]*P[7]+F[8]*P[8])/J; /* sz  */
  values [3] = (F[0]*P[3]+F[3]*P[4]+F[6]*P[5])/J; /* sxy */
  values [4] = (F[0]*P[6]+F[3]*P[7]+F[6]*P[8])/J; /* sxz */
  values [5] = (F[2]*P[0]+F[5]*P[1]+F[8]*P[2])/J; /* syz */
}

/* return element stabbed by a spatial point */
static ELEMENT* stabbed_element_cur (MESH *msh, ELEMENT **ele, int nele, double *x)
{
  ELEMENT **cur = ele, *ret;
  double dist, d;

  for (; nele > 0; ele ++, nele --)
    if (ELEMENT_Contains_Point (msh, *ele, x)) return *ele;

  for (dist = DBL_MAX, ret = NULL; cur < ele; cur ++)
  {
    d = ELEMENT_Point_Distance (msh, *cur, x, 0);
    if (d < dist) dist = d, ret = *cur;
  }

  return ret;
}

/* return element stabbed by a referential point */
static ELEMENT* stabbed_element (MESH *msh, ELEMENT **ele, int nele, double *X)
{
  ELEMENT **cur = ele, *ret;
  double dist, d;

  for (; nele > 0; ele ++, nele --)
    if (ELEMENT_Contains_Ref_Point (msh, *ele, X)) return *ele;

  for (dist = DBL_MAX, ret = NULL; cur < ele; cur ++)
  {
    d = ELEMENT_Point_Distance (msh, *cur, X, 1);
    if (d < dist) dist = d, ret = *cur;
  }

  return ret;
}

/* return transformation operator from the generalised to the local velocity space */
static MX* gen_to_loc_operator (BODY *bod, MESH *msh, ELEMENT *ele, double *X, double *base)
{
  int i [] = {0, 1, 2, 0, 1, 2, 0, 1, 2}, p [] = {0, 3, 6, 9};
  MX_CSC (base_trans, 9, 3, 3, p, i);
  double point [3];
  MX *N, *H;

  TNCOPY (base, base_trans.x);
  referential_to_local (msh, ele, X, point);
  N = shape_functions (bod, msh, ele, point);
  H = MX_Matmat (1.0, &base_trans, N, 0.0, NULL);
  MX_Destroy (N);
  return H;
}

/* accumulate constraints reaction */
inline static void fem_constraints_force_accum (BODY *bod, MESH *msh, ELEMENT *ele, double *X, double *base, double *R, short isma, double *force)
{
  double shapes [MAX_NODES_COUNT], P [3], point [3], *f;
  int numbers [MAX_NODES_COUNT], i, n;

  referential_to_local (msh, ele, X, point);
  n = shape_functions_ext (bod, msh, ele, point, shapes);
  load_node_numbers (bod, msh, ele, numbers);

  NVMUL (base, R, P);

  if (isma) for (i = 0; i < n; i ++) { f = &force [3 * numbers [i]]; ADDMUL (f, shapes [i], P, f); }
  else for (i = 0; i < n; i ++) { f = &force [3 * numbers [i]]; SUBMUL (f, shapes [i], P, f); }
}

/* compute constraints force */
static void fem_constraints_force (BODY *bod, double *force)
{
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  SET *node;
  int i;

  msh = FEM_MESH (bod);

  for (i = 0; i < bod->dofs; i ++) force [i] = 0.0;

  for (node = SET_First (bod->con); node; node = SET_Next (node))
  {
    CON *con = node->data;
    short isma = (bod == con->master);
    double *X = (isma ? con->mpnt : con->spnt);

    if (bod->msh)
    {
      cvx = (isma ? mgobj(con) : sgobj(con));
      ele = stabbed_element (msh, cvx->ele, cvx->nele, X); /* TODO: optimize */
    }
    else ele = (isma ? mgobj(con) : sgobj(con));

    fem_constraints_force_accum (bod, msh, ele, X, con->base, con->R, isma, force);
  }
}

/* compute inernal force */
static void fem_internal_force (BODY *bod, double *fint)
{
  MESH *msh = FEM_MESH (bod);
  double g [24], *v, *w;
  ELEMENT *ele;
  int bulk, i;

  for (i = 0; i < bod->dofs; i ++) fint [i] = 0.0;

  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    internal_force (0, bod, msh, ele, g);

    for (i = 0, v = g; i < ele->type; i ++, v += 3)
    {
      w = &fint [ele->nodes [i] * 3];
      ADD (w, v, w);
    }

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }
}
 
/* compute out of balance force = fext - fint */
static void fem_dynamic_force (BODY *bod, double time, double step, double *fext, double *fint, double *force)
{
  MESH *msh = FEM_MESH (bod);
  ELEMENT *ele;
  double g [24],
	 f [3],
	 point [3],
	 value,
	*v,
	*w;
  int bulk,
      i;

  /* zero forces */
  for (i = 0; i < bod->dofs; i ++) fext [i] = fint [i] = 0.0;

  /* point forces */
  for (FORCE *frc = bod->forces; frc; frc = frc->next)
  {
    if (frc->func)
    {
      ERRMEM (v = MEM_CALLOC (bod->dofs * sizeof (double)));
      frc->func (frc->data, frc->call, bod->dofs, bod->conf, bod->dofs, bod->velo, time, step, v);
      blas_daxpy (bod->dofs, 1.0, v, 1, fext, 1);
      free (v);
    }
    else
    {
      value = TMS_Value (frc->data, time);
      COPY (frc->direction, f);
      SCALE (f, value);

      ele = MESH_Element_Containing_Point (msh, frc->ref_point, 1);

      if (ele)
      {
	referential_to_local (msh, ele, frc->ref_point, point);

	if (frc->kind & CONVECTED)
	{ 
	  deformation_gradient (bod, msh, ele, point, g);
	  NVMUL (g, f, g+9);
	  COPY (g+9, f);
	}

	point_force (bod, msh, ele, point, f, fext);
      }
    }
  }

  /* gravitation */
  if (bod->dom->gravity [0])
  {
    f [0] = TMS_Value (bod->dom->gravity [0], time);
    f [1] = TMS_Value (bod->dom->gravity [1], time);
    f [2] = TMS_Value (bod->dom->gravity [2], time);

    for (ele = msh->surfeles, bulk = 0; ele; )
    {
      body_force (bod, msh, ele, f, g);

      for (i = 0, v = g; i < ele->type; i ++, v += 3)
      {
	w = &fext [ele->nodes [i] * 3];
	ADD (w, v, w);
      }

      if (bulk) ele = ele->next;
      else if (ele->next) ele = ele->next;
      else ele = msh->bulkeles, bulk = 1;
    }
  }

  /* internal forces */
  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    internal_force (0, bod, msh, ele, g);

    for (i = 0, v = g; i < ele->type; i ++, v += 3)
    {
      w = &fint [ele->nodes [i] * 3];
      ADD (w, v, w);
    }

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }

  /* force = fext - fint */
  for (double *x = fext, *y = fint, *z = force, *u = z + bod->dofs; z < u; x ++, y ++, z ++) *z = (*x) - (*y);
}

/* the smame computation for the static case */
#define fem_static_force(bod, time, step, fext, fint, force) fem_dynamic_force (bod,time,step,fext,fint,force)

/* compute global tangent stiffness */
static MX* tangent_stiffness (BODY *bod)
{
  struct colblock { int row; double val [3]; } *cb; /* column block */
  int i, j, k, l, *pp, *ii, *kk;
  double K [576], *A;
  MAP **col, *item;
  ELEMENT *ele;
  short bulk;
  MESH *msh;
  MEM blkmem,
      mapmem;
  MX *tang;

  MEM_Init  (&blkmem, sizeof (struct colblock), bod->dofs * 64);
  ERRMEM (col = MEM_CALLOC (sizeof (MAP*) * bod->dofs)); /* sparse columns */
  MEM_Init  (&mapmem, sizeof (MAP), bod->dofs * 64);
  msh = FEM_MESH (bod);

  for (ele = msh->surfeles, bulk = 0; ele;
       ele = (ele->next ? ele->next : bulk ? NULL : msh->bulkeles),
       bulk = (ele == msh->bulkeles ? 1 : bulk)) /* for each element in mesh */
  {
    internal_force (1, bod, msh, ele, K); /* compute internal force derivartive: K */

#if 0
  {
    double max = K [0];

    for (i = 1, l = ele->type * 3; i < l*l; i ++) max = MAX (max, K [i]);

    for (i = 0; i < l; i ++)
    {
      for (j = i+1; j < l; j ++)
      {
	ASSERT_DEBUG (fabs (K[l*j+i] - K[l*i+j]) < 1E-10 * max, "Unsymmetry of K: %e, %e", K[l*j+i], K[l*i+j]);
      }
    }
  }
#endif

    for (k = 0, A = K; k < ele->type; k ++) /* initialize K column block pointer; for element each node */
    {
      for (l = 0; l < 3; l ++) /* for each nodal degree of freedom */
      {
	j = 3 * ele->nodes [k] + l; /* for each global column index */

	for (i = 0; i < ele->type; i ++, A += 3) /* for each column row-block; shift column block pointer A */
	{
	  if (!(cb = MAP_Find (col [j], (void*) (long) ele->nodes [i], NULL))) /* if this row-block was not mapped */
	  {
	    ERRMEM (cb = MEM_Alloc (&blkmem));
	    cb->row = 3 * ele->nodes [i]; /* row-block initial index */
	    MAP_Insert (&mapmem, &col [j], (void*) (long) ele->nodes [i], cb, NULL); /* map it */
	  }
	  ACC (A, cb->val); /* accumulate values */
	}
      }
    }
  }

  ERRMEM (pp = malloc (sizeof (int [bod->dofs + 1]))); /* column pointers */

  for (pp [0] = 0, j = 0; j < bod->dofs; j ++) pp [j+1] = pp [j] + 3 * MAP_Size (col [j]); /* computed */

  ERRMEM (ii = malloc (sizeof (int [pp [bod->dofs]]))); /* row indices */

  for (j = 0, kk = ii; j < bod->dofs; j ++) /* initialize row index pointer; for each column */
  {
    for (item = MAP_First (col [j]); item; item = MAP_Next (item), kk += 3) /* for each row-block; increment column pointer by 3 */
    {
      cb = item->data;
      kk [0] = cb->row; /* set row indices */
      kk [1] = kk[0] + 1;
      kk [2] = kk[1] + 1;
    }
  }

  tang = MX_Create (MXCSC, bod->dofs, bod->dofs, pp, ii); /* create tangent matrix structure */

  for (j = 0, A = tang->x; j < bod->dofs; j ++) /* initialize column values pointer A; for each column */
  {
    for (item = MAP_First (col [j]); item; item = MAP_Next (item), A += 3) /* for each row-block, increment A */
    {
      cb = item->data;
      COPY (cb->val, A); /* copy values */
    }
  }

  free (ii);
  free (pp);
  free (col);
  MEM_Release (&mapmem);
  MEM_Release (&blkmem);

#if 0
  {
    double max = tang->x [0];

    for (i = 0; i < tang->p [tang->n]; i ++) max = MAX(max, tang-> x[i]);

    for (j = 0; j < tang->n; j ++)
    {
      for (k = tang->p [j]; k < tang->p [j+1]; k ++)
      {
	i = tang->i [k];
	for (l = tang->p [i]; l < tang->p [i+1]; l ++)
	{
	  if (tang->i [l] == j)
	  {
	    ASSERT_DEBUG (fabs (tang->x [k] - tang->x [l]) < 1E-10 * max, "Unsymmetric K: %e, %e\n", tang->x [k], tang->x [l]);
	    break;
	  }
	}
	ASSERT_DEBUG (l < tang->p [i+1], "Unsymmetric K: matching entry not found");
      }
    }
  }
#endif

  return tang;
}

/* compute diagonalized inertia operator */
static MX* diagonal_inertia (BODY *bod)
{
  MESH *msh = FEM_MESH (bod);
  int bulk,
     *p,
     *i,
      n,
      k;
  ELEMENT *ele;
  double *x;
  MX *M;

  n = bod->dofs;

  ERRMEM (p = malloc (sizeof (int [n+1])));
  ERRMEM (i = malloc (sizeof (int [n])));

  for (k = 0, p [n] = n; k < n; k ++) p [k] = i [k] = k; /* diagonal pattern */

  M = MX_Create (MXCSC, n, n, p, i);
  x = M->x;
  free (p);
  free (i);

  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    lump_mass_matrix (bod, msh, ele, x);

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }

  return M; 
}

/* set up inverse operator for the explicit dynamic time stepping */
static void fem_dynamic_explicit_inverse (BODY *bod)
{
  double *x, *y;

  ASSERT_DEBUG (!bod->M && !bod->inverse, "Body matrices not nil");

  bod->M = diagonal_inertia (bod);

  bod->inverse = MX_Copy (bod->M, NULL);

  for (x = bod->inverse->x, y = x + bod->dofs; x < y; x ++)
  {
    ASSERT (*x > 0.0, ERR_FEM_MASS_NOT_SPD);
    (*x) = 1.0 / (*x); /* invert diagonal */
  }
}

/* set up inverse operator for the implicit dynamic time stepping */
static void fem_dynamic_implicit_inverse (BODY *bod, double step, double *force)
{
  MX *M, *K, *A;

  if (bod->M) M = bod->M; else bod->M = M = diagonal_inertia (bod);

  if (bod->inverse) MX_Destroy (bod->inverse);

  K = tangent_stiffness (bod);

  if (force)
  {
    /* account for the previous velocity */
    MX_Matvec (1.0 / step, M, bod->velo, 1.0, force);

    /* account for the internal force increment */
    MX_Matvec (-0.25 * step, K, bod->velo, 1.0, force);
  }

  /* calculate tangent operator A = M + h*h/4 K */
  bod->inverse = A = MX_Add (1.0, M, 0.25*step*step, K, NULL);

  /* invert A */
  MX_Inverse (A, A);

  /* clean up */
  MX_Destroy (K);
}

/* solve implicit ingetration nonlineare problem */
static void fem_dynamic_implicit_solve (BODY *bod, double time, double step, double *fext, double *f, short begin)
{
  int n = bod->dofs,
      imax = 16,
      iter,
      i;

  double half = 0.5 * step,
	*fint = FEM_FINT (bod),
	*u0 = FEM_VEL0 (bod),
	*u = bod->velo,
	*q = bod->conf,
        *qorig, *aux, *res,
	errup, errlo, error;

  ERRMEM (qorig = malloc (sizeof (double [3 * n])));
  aux = qorig + n; res = aux + n;

  if (begin)
  {
    blas_dcopy (n, q, 1, qorig, 1); /* qorig = q (t) */
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    fem_dynamic_force (bod, time+half, step, fext, fint, f);  /* f(t+h/2,q(t+h/2)) = fext (t+h/2) - fint (q(t+h/2)) */
    fem_dynamic_implicit_inverse (bod, step, NULL); /* A = M + (h*h/4) * K */
    MX_Matvec (step, bod->inverse, f, 1.0, u); /* u(t+h) = u(t) + inv (A) * h * f */
  }
  else
  {
    blas_dcopy (n, q, 1, qorig, 1); /* qorig = q (t+h/2) */
    blas_daxpy (n, -half, u0, 1, qorig, 1); /* qorig = q(t) = q(t+h/2) - (h/2) * u(t) */
    MX_Matvec (step, bod->inverse, f, 1.0, u); /* u(t+h) += h * inv (M) * force */
  }

  iter = 0;
  do
  {
    /* in Simo and Tarnow paper, they show that PK2 should be taken as S=C[E(q(t+h)) + E(q(t))]/2,
     * nevertheless the simpler approach below where we use C[E((q(t+h) + q(t))/2)] did not lead
     * to energy inconsistent results in our test, but to the contratey. Hence we stick with it */
    for (i = 0; i < n; i ++) q [i] = qorig [i] + 0.25 * step * (u[i] + u0[i]); /* overwrite bod->conf ... */
    fem_internal_force (bod, fint); /* ... as it is used in there */
    fem_dynamic_implicit_inverse (bod, step, NULL);
    for (i = 0; i < n; i ++) { res [i] = step * (fext [i] - fint [i]); aux [i] = u [i] - u0[i]; }
    MX_Matvec (-1.0, bod->M, aux, 1.0, res);
    MX_Matvec (1.0, bod->inverse, res, 0.0, aux);
    for (i = 0; i < n; i ++) u [i] += aux [i];
    errlo = blas_ddot (n, u, 1, u, 1);
    errup = blas_ddot (n, aux, 1, aux, 1);
    error = sqrt (errup / MAX (errlo, 1.0));
  }
  while (error > 1E-8 && ++ iter < imax);

#if 0
  printf ("DEF_IMP: iter = %d, error = %e\n", iter, error);
#endif

  ASSERT (iter < imax, ERR_BOD_SCHEME_NOT_CONVERGED);

  if (begin)
  {
    for (i = 0; i < n; i ++) q [i] = qorig [i] + 0.5 * step * u0[i]; /* q(t+h/2) = q(t) + (h/2) * u(t) */
  }
  else
  {
    for (i = 0; i < n; i ++) q [i] = qorig [i] + 0.5 * step * (u0[i] + u [i]); /* q(t+h/2) = q(t) + (h/2) * (u(t+h) + u(t)) */
  }

  free (qorig);
}

/* compute inv (M) * K for an element */
static MX* inverse_mass_times_stiffencess (int explicit, BODY *bod, MESH *msh, ELEMENT *ele)
{
  double mass [24];
  double *x, *y;
  MX *IMK, *IM;
  int i, j, n;

  n = ele->type * 3;
  IM = explicit ? bod->inverse : bod->M;
  IMK = MX_Create (MXDENSE, n, n, NULL, NULL);
  internal_force (1, bod, msh, ele, IMK->x);

  for (i = 0; i < ele->type; i ++)
  {
    for (j = 0; j < 3; j ++)
    {
      mass [3*i+j] = IM->x [ele->nodes [i] * 3 + j];
      if (explicit == 0) mass [3*i+j] = 1.0 / mass [3*i+j];
    }
  }

  for (j = 0, x = IMK->x; j < n; j ++) /* compute IMK = IM * K */
  {
    for (i = 0, y = mass; i < n; i ++, x ++, y ++) (*x) *= (*y); /* scale each column by diagonal IM entries */
  }

  return IMK;
}

/* static time-stepping inverse */
static void fem_static_inverse (BODY *bod, double step)
{
  MX *M, *K, *A;

  if (bod->M) M = bod->M; else bod->M = M = diagonal_inertia (bod);

  if (bod->inverse) MX_Destroy (bod->inverse);

  K = tangent_stiffness (bod);

#if 0
  MESH *msh = FEM_MESH (bod);
  double eigmax;
  ELEMENT *ele;
  MX *IMK;

  /* estimate maximal eigenvalue of inv (M) * K based on just one element */
  ele = msh->surfeles;
  IMK = inverse_mass_times_stiffencess (0, bod, msh, ele); /* element inv (M) * K */
  MX_Scale (IMK, step * step * step);
  MX_Eigen (IMK, 1, &eigmax, NULL); /* maximal eigenvalue */
  ASSERT (eigmax > 0.0, ERR_BOD_MAX_FREQ_LE0);
  MX_Destroy (IMK);
#endif

  /* calculate tangent operator A = coef * M + h*h/4 K, where picking coef
   * seems tricky (eigmax/4.0 implies good damping, but we need to take care
   * for allowing a fair amount of the rigid motion as well: TODO: figure out) */
  bod->inverse = A = MX_Add (step*step*1E+6, M, step*step, K, NULL);

  /* invert A */
  MX_Inverse (A, A);

  /* clean up */
  MX_Destroy (K);
}

/* attach (element, local point) pairs to cvx->epn placeholder so that
 * current vertices can be updated fast by using FEM_Cur_Point_Ext */
static void post_proces_convices (SHAPE *shp, MESH *msh)
{
  CONVEX *cvx;
  ELEPNT *epn;
  double *ref;
  int n;

  for (; shp; shp = shp->next)
  {
    for (cvx = shp->data; cvx; cvx = cvx->next)
    {
      ERRMEM (cvx->epn = malloc (cvx->nver * sizeof (ELEPNT)));

      for (epn = cvx->epn, ref = cvx->ref, n = 0; n < cvx->nver; epn ++, ref += 3, n ++)
      {
	epn->ele = stabbed_element (msh, cvx->ele, cvx->nele, ref);
	ASSERT_DEBUG (epn->ele, "Invalid referential point stabbing an element");
        referential_to_local (msh, epn->ele, ref, epn->pnt);
      }
    }
  }
}

/* there are few things to be done here:
 * 1. Test whether the shape_volume == volume cut of the mesh
 * 2. Delete elements whose volume (ele->dom) < TOL * volume (ele)
 *    or delte ele->dom when volume (ele->dom) == volume (ele)
 * 3. Delete nodes unattached to elements
 * 4. Transform global ele->dom points into local element coordinates */
static void post_process_intersections (double shape_volume, MESH *msh)
{
  double cut_volume, ele_volume, point [3];
  ELEMENT *ele, *next;
  int bulk, n, k;
  TRISURF *surf;
  FACE *fac;
  TRI *tri;

  /* 1. */
 
  for (cut_volume = 0.0, ele = msh->surfeles; ele; ele = ele->next)
    for (surf = ele->dom, n = 0; n < ele->domnum; surf ++, n ++)
      cut_volume += TRI_Char (surf->tri, surf->m, surf->center); /* compute subdomain center as well */

  for (ele = msh->bulkeles; ele; ele = ele->next)
    for (surf = ele->dom, n = 0; n < ele->domnum; surf ++, n ++)
      cut_volume += TRI_Char (surf->tri, surf->m, surf->center);

  /* make sure that cut error is not too large, so that the mesh really contains the shape */
#if DEBUG
  if (fabs (shape_volume - cut_volume) >= CUT_TOL * shape_volume)
    printf ("shape_volume = %g, cut_volume = %g, error = %g\n",
      shape_volume, cut_volume, fabs (shape_volume - cut_volume) / shape_volume);
#endif
  ASSERT (fabs (shape_volume - cut_volume) < CUT_TOL * shape_volume, ERR_FEM_CUT_VOLUME);

  /* 2. */

  for (ele = msh->surfeles, bulk = 0; ele; ele = next)
  {
    next = ele->next;

    if (!ele->dom) /* delete element */
    {
      if (ele->prev) ele->prev->next = next;
      else if (bulk) msh->bulkeles = next;
      else msh->surfeles = next;
      if (next) next->prev = ele->prev;

      MEM_Free (&msh->elemem, ele);
    }
    else
    {
      ele_volume = ELEMENT_Volume (msh, ele, 1);
      cut_volume = 0.0;

      for (surf = ele->dom, n = 0; n < ele->domnum; surf ++, n ++)
	cut_volume += TRI_Char (surf->tri, surf->m, NULL);

      if (fabs (ele_volume - cut_volume) < GEOMETRIC_EPSILON * ele_volume) /* clear ele->dom */
      {
	for (surf = ele->dom, n = 0; n < ele->domnum; surf ++, n ++) free (surf->tri);
	free (ele->dom);
	ele->dom = NULL;
	ele->domnum = 0;
      }
    }

    if (!next && !bulk) next = msh->bulkeles, bulk = 1;
  }

  /* 3. */

  int *node_map, m;

  ERRMEM (node_map = MEM_CALLOC (msh->nodes_count * sizeof (int)));

  for (ele = msh->surfeles; ele; ele = ele->next)
    for (n = 0; n < ele->type; n ++) node_map [ele->nodes [n]] ++;

  for (ele = msh->bulkeles; ele; ele = ele->next)
    for (n = 0; n < ele->type; n ++) node_map [ele->nodes [n]] ++;

  for (n = 0, m = 1; n < msh->nodes_count; n ++)
    if (node_map [n]) node_map [n] = m ++;

  if (m < (n+1)) /* there are not referenced nodes */
  {
    for (ele = msh->surfeles; ele; ele = ele->next)
    {
      for (n = 0; n < ele->type; n ++) ele->nodes [n] = node_map [ele->nodes [n]] - 1;
      for (fac = ele->faces; fac; fac = fac->next) for (n = 0; n < fac->type; n ++) fac->nodes [n] = node_map [fac->nodes [n]] - 1;
    }

    for (ele = msh->bulkeles; ele; ele = ele->next)
      for (n = 0; n < ele->type; n ++) ele->nodes [n] = node_map [ele->nodes [n]] - 1;

    double (*ref) [3], (*mref) [3] = msh->ref_nodes,
	   (*cur) [3], (*mcur) [3] = msh->cur_nodes;

    ERRMEM (ref = malloc (2 * (m-1) * sizeof (double [3])));
    cur = ref + (m-1);

    for (n = 0; n < msh->nodes_count; n ++)
    {
      if (node_map [n])
      { 
	COPY (mref [n], ref [node_map [n] - 1]);
	COPY (mcur [n], cur [node_map [n] - 1]);
      }
    }

    free (msh->ref_nodes);

    msh->nodes_count = m-1;
    msh->ref_nodes = ref;
    msh->cur_nodes = cur;
  }

  free (node_map);

 /* 4. */

  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    for (surf = ele->dom, n = 0; n < ele->domnum; surf ++, n ++)
    {
      SET *done = NULL;

      COPY (surf->center, point); global_to_local (msh->ref_nodes, ele, point, surf->center); 

      for (tri = surf->tri, m = 0; m < surf->m; tri ++, m ++)
      {
	for (k = 0; k < 3; k ++)
	{ 
	  if (!SET_Find (done, tri->ver [k], NULL)) /* note that multiple triangles reference same nodes */
	  {
	    COPY (tri->ver [k], point);
	    global_to_local (msh->ref_nodes, ele, point, tri->ver [k]);
	    SET_Insert (NULL, &done, tri->ver [k], NULL);
	  }
	}
      }

      surf->volume = TRI_Char (surf->tri, surf->m, NULL); /* local element sub-volume for simplified integration */

      SET_Free (NULL, &done);
    }

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }
}

#if 0
/* dump intersection result for manual inverstigation for with tst/cvitest */
static void dump_intersection (CONVEX *cvx, double *vertices, double *planes, int n, int k, int m, TRI *tri, double *pla)
{
  double d, p[3], q[3];

  d = gjk (cvx->cur, cvx->nver, vertices, n, p, q);

  if ((!tri && d < GEOMETRIC_EPSILON) || (tri && TRI_Char (tri, m, NULL) < 0.0))
  {
    int i;
    printf ("%.24e\n", GEOMETRIC_EPSILON);
    printf ("%d   %d\n", cvx->nver, cvx->nfac);
    for (i = 0; i < cvx->nver; i ++) printf ("%.24e   %.24e   %.24e\n", cvx->cur[3*i], cvx->cur[3*i+1], cvx->cur[3*i+2]);
    for (i = 0; i < cvx->nfac; i ++) printf ("%.24e   %.24e   %.24e   %.24e   %.24e   %.24e\n", pla[6*i], pla[6*i+1], pla[6*i+2], pla[6*i+3], pla[6*i+4], pla[6*i+5]);
    printf ("%d   %d\n", n, k);
    for (i = 0; i < n; i ++) printf ("%.24e   %.24e   %.24e\n", vertices[3*i], vertices[3*i+1], vertices[3*i+2]);
    for (i = 0; i < k; i ++) printf ("%.24e   %.24e   %.24e   %.24e   %.24e   %.24e\n", planes[6*i], planes[6*i+1], planes[6*i+2], planes[6*i+3], planes[6*i+4], planes[6*i+5]);
  }
}
#endif

/* element and convex bounding boxe intersection callback; compute
 * their volumetric intersection to be later used for integration */
static void overlap (void *data, BOX *one, BOX *two)
{
  double vertices [24], planes [36], *pla;
  ELEMENT *ele;
  CONVEX *cvx;
  int m, n, k;
  TRI *tri;

  switch (GOBJ_Pair_Code (one, two))
  {
  case AABB_ELEMENT_CONVEX: ele = one->sgp->gobj; cvx = two->sgp->gobj; break;
  case AABB_CONVEX_ELEMENT: ele = two->sgp->gobj; cvx = one->sgp->gobj; break;
  default: ASSERT_DEBUG (0, "Unexpected pair code"); break;
  }

  n = ELEMENT_Vertices (data, ele, vertices);
  k = ELEMENT_Planes (data, ele, planes, NULL, NULL);
  pla = CONVEX_Planes (cvx);
  tri = cvi (cvx->cur, cvx->nver, pla, cvx->nfac, vertices, n, planes, k, REGULARIZED, &m);
#if 0
  dump_intersection (cvx, vertices, planes, n, k, m, tri, pla);
#endif
  free (pla);

  if (tri)
  {
#if DEBUG
    for (n = 0; n < cvx->nele; n ++) { ASSERT_DEBUG (cvx->ele [n] != ele, "CONVEX-ELEMENT intersection detected twice: this should not happen"); }
#endif

    ERRMEM (cvx->ele = realloc (cvx->ele, (++cvx->nele) * sizeof (ELEMENT*)));
    cvx->ele [cvx->nele-1] = ele;
    ERRMEM (ele->dom = realloc (ele->dom, (++ele->domnum) * sizeof (TRISURF)));
    ele->dom [ele->domnum-1].tri = tri;
    ele->dom [ele->domnum-1].m = m;
  }
}

/* create FEM internals for a body */
void FEM_Create (FEMFORM form, MESH *msh, SHAPE *shp, BULK_MATERIAL *mat, BODY *bod)
{
  /* compute shape characteristics */
  SHAPE_Char (shp, &bod->ref_volume, bod->ref_center, bod->ref_tensor);
  bod->ref_mass = bod->ref_volume * mat->density;
  SCALE9 (bod->ref_tensor, mat->density);

  if (msh) /* the given mesh is assumed to properly contain the shape */
  {
    SHAPE msh_shp = {SHAPE_MESH, msh, NULL, NULL};
    BOX **msh_boxes, **shp_boxes, **box;
    SGP *msh_sgp, *shp_sgp, *sgp, *sge;
    BOX_Extents_Update update;
    int msh_nsgp, shp_nsgp;
    ELEMENT *ele;
    MEM boxmem;

    msh_nsgp = msh->surfeles_count + msh->bulkeles_count;
    ERRMEM (msh_sgp = sgp = MEM_CALLOC (msh_nsgp * sizeof (SGP)));
    for (ele = msh->surfeles; ele; ele = ele->next, sgp ++) sgp->shp = &msh_shp, sgp->gobj = ele;
    for (ele = msh->bulkeles; ele; ele = ele->next, sgp ++) sgp->shp = &msh_shp, sgp->gobj = ele;
    shp_sgp = SGP_Create (shp, &shp_nsgp);
    MEM_Init (&boxmem, sizeof (BOX), msh_nsgp + shp_nsgp);
    ERRMEM (msh_boxes = malloc (msh_nsgp * sizeof (AABB*)));
    ERRMEM (shp_boxes = malloc (shp_nsgp * sizeof (AABB*)));

    for (sgp = msh_sgp, box = msh_boxes, sge = sgp + msh_nsgp; sgp < sge; sgp ++, box ++)
    {
      ERRMEM ((*box) = MEM_Alloc (&boxmem));
      update = SGP_Extents_Update (sgp);
      update (sgp->shp->data, sgp->gobj, (*box)->extents);
      (*box)->sgp = sgp;
      (*box)->kind = GOBJ_ELEMENT;
    }

    for (sgp = shp_sgp, box = shp_boxes, sge = sgp + shp_nsgp; sgp < sge; sgp ++, box ++)
    {
      ERRMEM ((*box) = MEM_Alloc (&boxmem));
      update = SGP_Extents_Update (sgp);
      update (sgp->shp->data, sgp->gobj, (*box)->extents);
      (*box)->sgp = sgp;
      (*box)->kind = GOBJ_CONVEX;
    }

    /* find overlaps between bounding boxes of mesh elements and convices */
    hybrid_ext (msh_boxes, msh_nsgp, shp_boxes, shp_nsgp, msh, overlap);

    free (shp_boxes);
    free (msh_boxes);
    MEM_Release (&boxmem);
    free (shp_sgp);
    free (msh_sgp);

    /* post-process intersections */
    post_process_intersections (bod->ref_volume, msh);

    /* post-proces convices */
    post_proces_convices (shp, msh);
  }
  else msh = shp->data; /* retrive the mesh pointer from the shape */

  switch (form)
  {
    case FEM_O1:
    {
      bod->dofs = msh->nodes_count * 3;
      ERRMEM (bod->conf = MEM_CALLOC (6 * bod->dofs * sizeof (double))); /* configuration, velocity, previous velocity, force, fext, fint */
      bod->velo = bod->conf + bod->dofs;
    }
    break;
    case FEM_O2:
    {
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
    }
    break;
  }
}

/* overwrite state */
void FEM_Overwrite_State (BODY *bod, double *q, double *u)
{
  blas_dcopy (bod->dofs, q, 1, bod->conf, 1);
  blas_dcopy (bod->dofs, u, 1, bod->velo, 1);
}

/* set initial rigid motion velocity */
void FEM_Initial_Velocity (BODY *bod, double *linear, double *angular)
{
  MESH *msh = FEM_MESH (bod);
  double (*nod) [3] = msh->ref_nodes,
	 (*end) [3] = nod + msh->nodes_count,
	 *velo = bod->velo,
	 *X0 = bod->ref_center,
	 A [3];

  for (; nod < end; nod ++, velo += 3)
  {
    if (linear) { COPY (linear, velo); }

    if (angular)
    {
      SUB (nod [0], X0, A);
      PRODUCTADD (angular, A, velo);
    }
  }
}

/* initialise dynamic time stepping */
void FEM_Dynamic_Init (BODY *bod)
{
  if (bod->scheme == SCH_DEF_EXP)
  {
    if (!bod->inverse) fem_dynamic_explicit_inverse (bod); /* initialize once */
  }
  else fem_dynamic_implicit_inverse (bod, bod->dom->step, NULL); /* update every time */
}

/* estimate critical step for the dynamic scheme */
double FEM_Dynamic_Critical_Step (BODY *bod)
{
  if (bod->scheme == SCH_DEF_EXP)
  {
    MESH *msh = FEM_MESH (bod);
    double step, tcrit, eigmax;
    ELEMENT *ele;
    int bulk;
    MX *IMK;

    for (ele = msh->surfeles, bulk = 0, step = DBL_MAX; ele; )
    {
      IMK = inverse_mass_times_stiffencess (1, bod, msh, ele); /* element inv (M) * K */
      MX_Eigen (IMK, 1, &eigmax, NULL); /* maximal eigenvalue */
      MX_Destroy (IMK);
      ASSERT (eigmax > 0.0, ERR_BOD_MAX_FREQ_LE0);
      tcrit = 2.0 / sqrt (eigmax); /* limit of stability => t_crit <= 2.0 / omega_max */
      if (tcrit < step) step = tcrit;

      if (bulk) ele = ele->next;
      else if (ele->next) ele = ele->next;
      else ele = msh->bulkeles, bulk = 1;
    }

    return step;
  }
  else return DBL_MAX;
}

/* perform the initial half-step of the dynamic scheme */
void FEM_Dynamic_Step_Begin (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  double half = 0.5 * step,
	 c = bod->damping,
	*x = bod->inverse->x,
	*u0 = FEM_VEL0 (bod),
	*f = FEM_FORCE (bod),
	*fext = FEM_FEXT (bod),
	*fint = FEM_FINT (bod),
	*q = bod->conf,
	*u = bod->velo,
	*e = u + n;

  blas_dcopy (n, u, 1, u0, 1); /* save u (t) */

  switch (bod->scheme)
  {
  case SCH_DEF_EXP:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    fem_dynamic_force (bod, time+half, step, fext, fint, f);  /* f = fext (t+h/2) - fint (q(t+h/2)) */
    for (; u < e; u ++, x ++, f ++) (*u) += step * (*x) * (*f); /* u(t+h) = u(t) + inv (M) * h * f */
  }
  break;
  case SCH_DEF_LIM:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    fem_dynamic_force (bod, time+half, step, fext, fint, f);  /* f = fext (t+h/2) - fint (q(t+h/2)) */
    fem_dynamic_implicit_inverse (bod, step, NULL); /* A = M + (h*h/4) * K */
    MX_Matvec (step, bod->inverse, f, 1.0, u); /* u(t+h) = u(t) + inv (A) * h * f */
  }
  break;
  case SCH_DEF_LIM2:
  {
    fem_dynamic_force (bod, time+half, step, fext, fint, f);  /* f = fext (t+h/2) - fint (q(t)) */
    fem_dynamic_implicit_inverse (bod, step, f); /* f += (1/h) M u(t) - (h/4) K u (t) */
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    MX_Matvec (step, bod->inverse, f, 0.0, u); /* u(t+h) = inv (A) * h * force */
  }
  break;
  case SCH_DEF_IMP:
  {
    /* q(t+h/2) = q(t) + (h/2) * u(t)
     * f = fext (t+h/2) - fint ([q(t) + q(t+h)]/2) 
     * A = M + (h*h/4) * K ([q(t) + q(t+h)]/2) 
     * u (t+h) = u (t) + inv (A) * h * f */
    fem_dynamic_implicit_solve (bod, time, step, fext, f, 1);
  }
  break;
  default:
  break;
  }

  if (c > 0.0) for (u = bod->velo; u < e; u ++, u0++) (*u) -= c * (*u0); /* u(t+h) -= c * u (t) */
}

/* perform the final half-step of the dynamic scheme */
void FEM_Dynamic_Step_End (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  double half = 0.5 * step,
	*energy = bod->energy,
	*x = bod->inverse->x,
	*r = FEM_FORCE (bod),
	*ir = r,
	*fext = FEM_FEXT (bod),
	*fint = FEM_FINT (bod),
	*u0 = FEM_VEL0 (bod),
	*iu0 = u0,
	*u = bod->velo,
	*iu = u,
        *q = bod->conf,
	*e = u + n;

  fem_constraints_force (bod, r); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  blas_daxpy (n, 1.0, r, 1, fext, 1);  /* fext += r */

  switch (bod->scheme)
  {
  case SCH_DEF_EXP:
  {
    for (; iu < e; iu ++, x ++, ir ++) (*iu) += step * (*x) * (*ir); /* u(t+h) += inv (M) * h * r */
    blas_daxpy (n, half, u, 1, q, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
  }
  break;
  case SCH_DEF_LIM:
  case SCH_DEF_LIM2:
  {
    MX_Matvec (step, bod->inverse, r, 1.0, u); /* u(t+h) += h * inv (M) * force */
    blas_daxpy (n, half, u, 1, q, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
  }
  break;
  case SCH_DEF_IMP:
  {
    /* f = fext (t+h/2) - fint ([q(t) + q(t+h)]/2) 
     * A = M + (h*h/4) * K ([q(t) + q(t+h)]/2) 
     * u (t+h) = u (t) + inv (A) * h * f
     * q(t+h) = q(t+h/2) + (h/2) * u(t+h) */
    fem_dynamic_implicit_solve (bod, time, step, fext, r, 0);
  }
  break;
  default:
  break;
  }

  /* energy */
  for (ir = r, iu = u; iu < e; ir ++, iu ++, iu0 ++) *ir = half * ((*iu) + (*iu0)); /* dq = (h/2) * {u(t) + u(t+h)} */
  energy [EXTERNAL] += blas_ddot (n, r, 1, fext, 1);
  energy [INTERNAL] += blas_ddot (n, r, 1, fint, 1);

  if (bod->msh) /* in such case SHAPE_Update will not update "rough" mesh */
  {
    MESH *msh = bod->msh;
    double (*cur) [3] = msh->cur_nodes,
	   (*ref) [3] = msh->ref_nodes,
	   (*end) [3] = ref + msh->nodes_count;

    for (; ref < end; cur ++, ref ++, q += 3) { ADD (ref[0], q, cur[0]); }
  }
}

/* initialise static time stepping */
void FEM_Static_Init (BODY *bod)
{
  fem_static_inverse (bod, bod->dom->step);
}

/* perform the initial half-step of the static scheme */
void FEM_Static_Step_Begin (BODY *bod, double time, double step)
{
  double *force = FEM_FORCE (bod);

  fem_static_inverse (bod, step); /* compute inverse of static tangent operator */
  fem_static_force (bod, time+step, step, FEM_FEXT(bod), FEM_FINT(bod), force);  /* f(t+h) = fext (t+h) - fint (q(t+h)) */
  MX_Matvec (step, bod->inverse, force, 0.0, bod->velo); /* u(t+h) = inv (A) * h * f(t+h) */
}

/* perform the final half-step of the static scheme */
void FEM_Static_Step_End (BODY *bod, double time, double step)
{
  double *force = FEM_FORCE (bod);

  fem_constraints_force (bod, force); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) += inv (A) * h * r */
  blas_daxpy (bod->dofs, step, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t) + h * u(t+h) */

  if (bod->msh) /* in such case SHAPE_Update will not update "rough" mesh */
  {
    MESH *msh = bod->msh;
    double (*cur) [3] = msh->cur_nodes,
	   (*ref) [3] = msh->ref_nodes,
	   (*end) [3] = ref + msh->nodes_count,
	    *q = bod->conf;

    for (; ref < end; cur ++, ref ++, q += 3) { ADD (ref[0], q, cur[0]); }
  }
}

/* motion x = x (X, t) */
void FEM_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x)
{
  double shapes [MAX_NODES_COUNT], q [MAX_NODES_COUNT][3], point [3];
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  int n, i;

  if (bod->msh)
  {
    msh = bod->msh;
    ASSERT_DEBUG_EXT ((cvx = gobj), "NULL geometric object for the separate mesh FEM scenario");
    ele = stabbed_element (msh, cvx->ele, cvx->nele, X);
    ASSERT_DEBUG (ele, "Invalid referential point stabbing an element");
  }
  else
  {
    msh = shp->data;
    ele = gobj;
  }

  if (ele)
  {
    referential_to_local (msh, ele, X, point);
    n = shape_functions_ext (bod, msh, ele, point, shapes);
    load_displacements (bod, msh, ele, q);

    COPY (X, x); for (i = 0; i < n; i ++) { ADDMUL (x, shapes [i], q [i], x); } /* x = X + N q */
  }
  else /* NULL implies nodal update (X is within msh->ref_nodes) */
  {
    int n = (node_t) X - msh->ref_nodes;
    double *q = &bod->conf [n * 3];

    ADD (msh->ref_nodes [n], q, x);
  }
}

/* motion x = x (element, ref point, local point) */
void FEM_Cur_Point_Ext (BODY *bod, ELEMENT *ele, double *X, double *point, double *x)
{
  double shapes [MAX_NODES_COUNT], q [MAX_NODES_COUNT][3];
  MESH *msh = FEM_MESH (bod);
  int n, i;

  n = shape_functions_ext (bod, msh, ele, point, shapes);
  load_displacements (bod, msh, ele, q);

  COPY (X, x); for (i = 0; i < n; i ++) { ADDMUL (x, shapes [i], q [i], x); } /* x = X + N q */
}

/* inverse motion X = X (x, state) */
void FEM_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X)
{
  double point [3];
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;

  if (bod->msh)
  {
    msh = bod->msh;
    ASSERT_DEBUG_EXT ((cvx = gobj), "NULL geometric object for the separate mesh FEM scenario");
    ele = stabbed_element_cur (msh, cvx->ele, cvx->nele, x);
    ASSERT_DEBUG (ele, "Invalid spatial point stabbing an element");
  }
  else
  {
    msh = shp->data;
    ele = gobj;
  }

  spatial_to_local (msh, ele, x, point);

  local_to_referential (msh, ele, point, X);
}

/* obtain spatial velocity at (gobj, referential point), expressed in the local spatial 'base' */
void FEM_Local_Velo (BODY *bod, SHAPE *shp, void *gobj, double *X, double *base, double *prevel, double *curvel)
{
  double shapes [MAX_NODES_COUNT], u0 [MAX_NODES_COUNT][3], u [MAX_NODES_COUNT][3], point [3], vglo [3];
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  int i, n;

  if (bod->msh)
  {
    msh = bod->msh;
    ASSERT_DEBUG_EXT ((cvx = gobj), "NULL geometric object for the separate mesh FEM scenario");
    ele = stabbed_element (msh, cvx->ele, cvx->nele, X);
    ASSERT_DEBUG (ele, "Invalid referential point stabbing an element");
  }
  else
  {
    msh = shp->data;
    ele = gobj;
  }

  referential_to_local (msh, ele, X, point);
  n = shape_functions_ext (bod, msh, ele, point, shapes);
  load_velocities (bod, bod->msh, ele, u0, u);

  if (prevel)
  {
    SET (vglo, 0); for (i = 0; i < n; i ++) { ADDMUL (vglo, shapes [i], u0 [i], vglo); } /* vglo = N u0 */
    TVMUL (base, vglo, prevel); /* prevel = base' vglo */
  }

  if (curvel)
  {
    SET (vglo, 0); for (i = 0; i < n; i ++) { ADDMUL (vglo, shapes [i], u [i], vglo); } /* vglo = N u */
    TVMUL (base, vglo, curvel); /* prevel = base' vglo */
  }
}

/* return transformation operator from the generalised to the local velocity space at (element, ref. point, base) */
MX* FEM_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *X, double *base)
{
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;

  if (bod->msh)
  {
    msh = bod->msh;
    ASSERT_DEBUG_EXT ((cvx = gobj), "NULL geometric object for the separate mesh FEM scenario");
    ele = stabbed_element (msh, cvx->ele, cvx->nele, X);
    ASSERT_DEBUG (ele, "Invalid referential point stabbing an element");
  }
  else
  {
    msh = shp->data;
    ele = gobj;
  }

  return gen_to_loc_operator (bod, msh, ele, X, base);
}

/* compute current kinetic energy */
double FEM_Kinetic_Energy (BODY *bod)
{
  if (bod->M)
  {
    double *x = bod->M->x,
	   *y = x + bod->dofs,
	   *u = bod->velo,
	   sum;

    for (sum = 0.0; x < y; x ++, u ++) sum += (*u)*(*u) * (*x);

    return 0.5 * sum;
  }

  return 0.0;
}

/* get some values at a referential point */
void FEM_Point_Values (BODY *bod, double *X, VALUE_KIND kind, double *values)
{
  MESH *msh = FEM_MESH (bod);
  double point [3];
  ELEMENT *ele;

  ele = MESH_Element_Containing_Point (msh, X, 1);

  if (ele)
  {
    referential_to_local (msh, ele, X, point);

    switch (kind)
    {
    case VALUE_DISPLACEMENT:
    {
      MX *N = shape_functions (bod, msh, ele, point);
      MX_Matvec (1.0, N, bod->conf, 0.0, values);
      MX_Destroy (N);
    }
    break;
    case VALUE_VELOCITY:
    {
      MX *N = shape_functions (bod, msh, ele, point);
      MX_Matvec (1.0, N, bod->velo, 0.0, values);
      MX_Destroy (N);
    }
    break;
    case VALUE_STRESS:
    {
      fem_element_cauchy (bod, msh, ele, point, values);
    }
    break;
    case VALUE_MISES:
    {
      double stress [6];

      fem_element_cauchy (bod, msh, ele, point, stress);
      MISES (stress, values [0]);
    }
    break;
    case VALUE_STRESS_AND_MISES:
    {
      fem_element_cauchy (bod, msh, ele, point, values);
      MISES (values, values [6]);
    }
    break;
    }
  }
}

/* get some values at a local point of an element */
void FEM_Element_Point_Values (BODY *bod, ELEMENT *ele, double *point, VALUE_KIND kind, double *values)
{
  MESH *msh = FEM_MESH (bod);

  switch (kind)
  {
  case VALUE_DISPLACEMENT:
  {
    MX *N = shape_functions (bod, msh, ele, point);
    MX_Matvec (1.0, N, bod->conf, 0.0, values);
    MX_Destroy (N);
  }
  break;
  case VALUE_VELOCITY:
  {
    MX *N = shape_functions (bod, msh, ele, point);
    MX_Matvec (1.0, N, bod->velo, 0.0, values);
    MX_Destroy (N);
  }
  break;
  case VALUE_STRESS:
  {
    fem_element_cauchy (bod, msh, ele, point, values);
  }
  break;
  case VALUE_MISES:
  {
    double stress [6];

    fem_element_cauchy (bod, msh, ele, point, stress);
    MISES (stress, values [0]);
  }
  break;
  case VALUE_STRESS_AND_MISES:
  {
    fem_element_cauchy (bod, msh, ele, point, values);
    MISES (values, values [6]);
  }
  break;
  }
}

/* get some values at a curent mesh node */
void FEM_Cur_Node_Values (BODY *bod, double *node, VALUE_KIND kind, double *values)
{
  MESH *msh = FEM_MESH (bod);
  int n = (node_t) node - msh->cur_nodes;
  ELEMENT *ele = MAP_Find (msh->map, node, NULL),
	  *start [2] = {msh->surfeles, msh->bulkeles};
  double point [3] = {0, 0, 0};
  int i, j;

  if (ele == NULL) /* if this node has not been yet mapped to an element */
  {
    for(j = 0; j < 2; j ++)
    {
      for (ele = start [j]; ele; ele = ele->next) /* make a costly linear search for an element */
      {
	for (i = 0; i < ele->type; i ++)
	{
	  if (ele->nodes [i] == n)
	  {
	    switch (ele->type)
	    {
	    case 4: if (i < 3) point [i] = 1.0; break;
	    case 8:
	      switch (i)
	      {
	      case 0: point [0] = -1; point [1] = -1; point [2] = -1; break;
	      case 1: point [0] =  1; point [1] = -1; point [2] = -1; break;
	      case 2: point [0] =  1; point [1] =  1; point [2] = -1; break;
	      case 3: point [0] = -1; point [1] =  1; point [2] = -1; break;
	      case 4: point [0] = -1; point [1] = -1; point [2] =  1; break;
	      case 5: point [0] =  1; point [1] = -1; point [2] =  1; break;
	      case 6: point [0] =  1; point [1] =  1; point [2] =  1; break;
	      case 7: point [0] = -1; point [1] =  1; point [2] =  1; break;
	      }
	      break;
	    case 6:
	      switch (i)
	      {
	      case 0: point [0] = -1; point [1] = -1; point [2] = -1; break;
	      case 1: point [0] =  1; point [1] = -1; point [2] = -1; break;
	      case 2: point [0] =  1; point [1] =  1; point [2] = -1; break;
	      case 3: point [0] = -1; point [1] = -1; point [2] =  1; break;
	      case 4: point [0] =  1; point [1] = -1; point [2] =  1; break;
	      case 5: point [0] =  1; point [1] =  1; point [2] =  1; break;
	      }
	      break;
	    case 5:
	      switch (i)
	      {
	      case 0: point [0] = -1; point [1] = -1; point [2] = -1; break;
	      case 1: point [0] =  1; point [1] = -1; point [2] = -1; break;
	      case 2: point [0] =  1; point [1] =  1; point [2] = -1; break;
	      case 3: point [0] = -1; point [1] =  1; point [2] = -1; break;
	      case 4: point [0] = -1; point [1] = -1; point [2] =  1; break;
	      }
	      break;
	    }

	    MAP_Insert (&msh->mapmem, &msh->map, node, ele, NULL); /* and map it */
	    goto ok;
	  }
	}
      }
    }
  }

  ASSERT_DEBUG (ele, "Element with specified node number was not found");

ok:
  switch (kind)
  {
  case VALUE_DISPLACEMENT:
  {
    double *q = &bod->conf [3 * n];
    COPY (q, values);
  }
  break;
  case VALUE_VELOCITY:
  {
    double *u = &bod->velo [3 * n];
    COPY (u, values);
  }
  break;
  case VALUE_STRESS:
  {
    fem_element_cauchy (bod, msh, ele, point, values); /* TODO: average from neighouring elements */
  }
  break;
  case VALUE_MISES:
  {
    double stress [6];

    fem_element_cauchy (bod, msh, ele, point, stress); /* TODO: average from neighouring elements */
    MISES (stress, values [0]);
  }
  break;
  case VALUE_STRESS_AND_MISES:
  {
    fem_element_cauchy (bod, msh, ele, point, values); /* TODO: average from neighouring elements */
    MISES (values, values [6]);
  }
  break;
  }
}

/* issued by state reading routines of body interface */
void FEM_Update_Rough_Mesh (BODY *bod)
{
  MESH *msh = bod->msh;
  double *q = bod->conf,
	(*cur) [3] = msh->cur_nodes,
	(*ref) [3] = msh->ref_nodes,
	(*end) [3] = ref + msh->nodes_count;

  for (; ref < end; ref ++, cur ++, q += 3) { ADD (ref[0], q, cur[0]); }
}

/* release FEM memory */
void FEM_Destroy (BODY *bod)
{
  free (bod->conf);
}

#if MPI
/* get configuration packing size */
int FEM_Conf_Pack_Size (BODY *bod)
{
  return bod->dofs;
}

/* get velocity packing size */
int FEM_Velo_Pack_Size (BODY *bod)
{
  return 5 * bod->dofs;
}
#endif
