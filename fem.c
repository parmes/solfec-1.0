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

#include <string.h>
#include <float.h>
#include "lap.h"
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

#define IMP_EPS 1E-8
#define MAX_ITERS 16
#define MAX_NODES 20
#define DOM_TOL 0.150
#define CUT_TOL 0.015
#define FEM_VEL0(bod) ((bod)->velo + (bod)->dofs)     /* previous velocity */
#define FEM_FEXT(bod) ((bod)->velo + (bod)->dofs * 2) /* external force */
#define FEM_FINT(bod) ((bod)->velo + (bod)->dofs * 3) /* internal force */
#define FEM_FBOD(bod) ((bod)->velo + (bod)->dofs * 4) /* unit body force */
#define FEM_ROT(bod) ((bod)->conf + (bod)->dofs)      /* rotation */
#define FEM_MESH(bod) ((bod)->msh ? (bod)->msh : (bod)->shape->data)
#define FEM_MATERIAL(bod, ele) ((ele)->mat ? (ele)->mat : (bod)->mat)

/* ==================== INTEGRATION ======================= */

/* order 1 */
static const double I_TET1_X[] = {0.25};
static const double I_TET1_Y[] = {0.25};
static const double I_TET1_Z[] = {0.25};
static const double I_TET1_W[] = {0.16666666666666666};
#define             I_TET1_N      1

static const double I_PYR1_X[] = {0}; /* TODO */
static const double I_PYR1_Y[] = {0};
static const double I_PYR1_Z[] = {0};
static const double I_PYR1_W[] = {0};
#define             I_PYR1_N      0

static const double I_WED1_X[] = {0}; /* TODO */
static const double I_WED1_Y[] = {0};
static const double I_WED1_Z[] = {0};
static const double I_WED1_W[] = {0};
#define             I_WED1_N      0

static const double I_HEX1_X[] = {0}; /* TODO */
static const double I_HEX1_Y[] = {0};
static const double I_HEX1_Z[] = {0};
static const double I_HEX1_W[] = {0};
#define             I_HEX1_N      0

static const double I_TRI1_X[] = {0.333333333333333333};
static const double I_TRI1_Y[] = {0.333333333333333333};
static const double I_TRI1_W[] = {0.5};
#define             I_TRI1_N      1

static const double I_QUA1_X [] = {0.0};
static const double I_QUA1_Y [] = {0.0};
static const double I_QUA1_W [] = {4.0};
#define             I_QUA1_N       1

/* order 2 */
static const double I_TET2_X [] = {0.13819660112501052, 0.13819660112501052, 0.13819660112501052, 0.58541019662496840};
static const double I_TET2_Y [] = {0.13819660112501052, 0.13819660112501052, 0.58541019662496840, 0.13819660112501052};
static const double I_TET2_Z [] = {0.13819660112501052, 0.58541019662496840, 0.13819660112501052, 0.13819660112501052};
static const double I_TET2_W [] = {0.04166666666666666, 0.04166666666666666, 0.04166666666666666, 0.04166666666666666};
#define             I_TET2_N       4

static const double I_PYR2_X[] = {0}; /* TODO */
static const double I_PYR2_Y[] = {0};
static const double I_PYR2_Z[] = {0};
static const double I_PYR2_W[] = {0};
#define             I_PYR2_N      0

static const double I_WED2_X[] = {0}; /* TODO */
static const double I_WED2_Y[] = {0};
static const double I_WED2_Z[] = {0};
static const double I_WED2_W[] = {0};
#define             I_WED2_N      0

#define ISQR3 0.57735026918962584 
static const double I_HEX2_X [] = {-ISQR3,  ISQR3,  ISQR3, -ISQR3, -ISQR3,  ISQR3, ISQR3, -ISQR3};
static const double I_HEX2_Y [] = {-ISQR3, -ISQR3,  ISQR3,  ISQR3, -ISQR3, -ISQR3, ISQR3,  ISQR3};
static const double I_HEX2_Z [] = {-ISQR3, -ISQR3, -ISQR3, -ISQR3,  ISQR3,  ISQR3, ISQR3,  ISQR3};
static const double I_HEX2_W [] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
#define             I_HEX2_N       8

static const double I_TRI2_X[] = {0.66666666666666667, 0.16666666666666667, 0.16666666666666667};
static const double I_TRI2_Y[] = {0.16666666666666667, 0.66666666666666667, 0.16666666666666667};
static const double I_TRI2_W[] = {0.16666666666666667, 0.16666666666666667, 0.16666666666666667};
#define             I_TRI2_N      3

static const double I_QUA2_X [] = {-ISQR3,  ISQR3, ISQR3, -ISQR3};
static const double I_QUA2_Y [] = {-ISQR3, -ISQR3, ISQR3,  ISQR3};
static const double I_QUA2_W [] = {1.0, 1.0, 1.0, 1.0};
#define             I_QUA2_N       4

/* order 3 */
static const double I_TET3_X[] = {0}; /* TODO */
static const double I_TET3_Y[] = {0};
static const double I_TET3_Z[] = {0};
static const double I_TET3_W[] = {0};
#define             I_TET3_N      0

static const double I_PYR3_X[] = {0}; /* TODO */
static const double I_PYR3_Y[] = {0};
static const double I_PYR3_Z[] = {0};
static const double I_PYR3_W[] = {0};
#define             I_PYR3_N      0

static const double I_WED3_X[] = {0}; /* TODO */
static const double I_WED3_Y[] = {0};
static const double I_WED3_Z[] = {0};
static const double I_WED3_W[] = {0};
#define             I_WED3_N      0

static const double I_HEX3_X[] = {0}; /* TODO */
static const double I_HEX3_Y[] = {0};
static const double I_HEX3_Z[] = {0};
static const double I_HEX3_W[] = {0};
#define             I_HEX3_N      0

static const double I_TRI3_X[] = {0}; /* TODO */
static const double I_TRI3_Y[] = {0};
static const double I_TRI3_W[] = {0};
#define             I_TRI3_N      0

static const double I_QUA3_X [] = {0}; /* TODO */
static const double I_QUA3_Y [] = {0};
static const double I_QUA3_W [] = {0};
#define             I_QUA3_N       0

/* order 4 */
static const double I_TET4_X[] = {0}; /* TODO */
static const double I_TET4_Y[] = {0};
static const double I_TET4_Z[] = {0};
static const double I_TET4_W[] = {0};
#define             I_TET4_N      0

static const double I_PYR4_X[] = {0}; /* TODO */
static const double I_PYR4_Y[] = {0};
static const double I_PYR4_Z[] = {0};
static const double I_PYR4_W[] = {0};
#define             I_PYR4_N      0

static const double I_WED4_X[] = {0}; /* TODO */
static const double I_WED4_Y[] = {0};
static const double I_WED4_Z[] = {0};
static const double I_WED4_W[] = {0};
#define             I_WED4_N      0

static const double I_HEX4_X[] = {0}; /* TODO */
static const double I_HEX4_Y[] = {0};
static const double I_HEX4_Z[] = {0};
static const double I_HEX4_W[] = {0};
#define             I_HEX4_N      0

static const double I_TRI4_X[] = {0}; /* TODO */
static const double I_TRI4_Y[] = {0};
static const double I_TRI4_W[] = {0};
#define             I_TRI4_N      0

static const double I_QUA4_X [] = {0}; /* TODO */
static const double I_QUA4_Y [] = {0};
static const double I_QUA4_W [] = {0};
#define             I_QUA4_N       0

#define MAX_ORDER 4

static const double *I_TET_X [] = {NULL, I_TET1_X, I_TET2_X, I_TET3_X, I_TET4_X};
static const double *I_TET_Y [] = {NULL, I_TET1_Y, I_TET2_Y, I_TET3_Y, I_TET4_Y};
static const double *I_TET_Z [] = {NULL, I_TET1_Z, I_TET2_Z, I_TET3_Z, I_TET4_Z};
static const double *I_TET_W [] = {NULL, I_TET1_W, I_TET2_W, I_TET3_W, I_TET4_W};
static const int     I_TET_N [] = {   0, I_TET1_N, I_TET2_N, I_TET3_N, I_TET4_N};

static const double *I_PYR_X [] = {NULL, I_PYR1_X, I_PYR2_X, I_PYR3_X, I_PYR4_X};
static const double *I_PYR_Y [] = {NULL, I_PYR1_Y, I_PYR2_Y, I_PYR3_Y, I_PYR4_Y};
static const double *I_PYR_Z [] = {NULL, I_PYR1_Z, I_PYR2_Z, I_PYR3_Z, I_PYR4_Z};
static const double *I_PYR_W [] = {NULL, I_PYR1_W, I_PYR2_W, I_PYR3_W, I_PYR4_W};
static const int     I_PYR_N [] = {   0, I_PYR1_N, I_PYR2_N, I_PYR3_N, I_PYR4_N};

static const double *I_WED_X [] = {NULL, I_WED1_X, I_WED2_X, I_WED3_X, I_WED4_X};
static const double *I_WED_Y [] = {NULL, I_WED1_Y, I_WED2_Y, I_WED3_Y, I_WED4_Y};
static const double *I_WED_Z [] = {NULL, I_WED1_Z, I_WED2_Z, I_WED3_Z, I_WED4_Z};
static const double *I_WED_W [] = {NULL, I_WED1_W, I_WED2_W, I_WED3_W, I_WED4_W};
static const int     I_WED_N [] = {   0, I_WED1_N, I_WED2_N, I_WED3_N, I_WED4_N};

static const double *I_HEX_X [] = {NULL, I_HEX1_X, I_HEX2_X, I_HEX3_X, I_HEX4_X};
static const double *I_HEX_Y [] = {NULL, I_HEX1_Y, I_HEX2_Y, I_HEX3_Y, I_HEX4_Y};
static const double *I_HEX_Z [] = {NULL, I_HEX1_Z, I_HEX2_Z, I_HEX3_Z, I_HEX4_Z};
static const double *I_HEX_W [] = {NULL, I_HEX1_W, I_HEX2_W, I_HEX3_W, I_HEX4_W};
static const int     I_HEX_N [] = {   0, I_HEX1_N, I_HEX2_N, I_HEX3_N, I_HEX4_N};

static const double *I_TRI_X [] = {NULL, I_TRI1_X, I_TRI2_X, I_TRI3_X, I_TRI4_X};
static const double *I_TRI_Y [] = {NULL, I_TRI1_Y, I_TRI2_Y, I_TRI3_Y, I_TRI4_Y};
static const double *I_TRI_W [] = {NULL, I_TRI1_W, I_TRI2_W, I_TRI3_W, I_TRI4_W};
static const int     I_TRI_N [] = {   0, I_TRI1_N, I_TRI2_N, I_TRI3_N, I_TRI4_N};

static const double *I_QUA_X [] = {NULL, I_QUA1_X, I_QUA2_X, I_QUA3_X, I_QUA4_X};
static const double *I_QUA_Y [] = {NULL, I_QUA1_Y, I_QUA2_Y, I_QUA3_Y, I_QUA4_Y};
static const double *I_QUA_W [] = {NULL, I_QUA1_W, I_QUA2_W, I_QUA3_W, I_QUA4_W};
static const int     I_QUA_N [] = {   0, I_QUA1_N, I_QUA2_N, I_QUA3_N, I_QUA4_N};

/* load 3D integrator data */
inline static int integrator3d_load (int type, int order, const double **X, const double **Y, const double **Z, const double **W)
{
  int N;

  ASSERT_DEBUG (order >= 1 && order < MAX_ORDER, "Integration order out of bounds");

  switch (type)
  {
  case 4:
  case 10:
  {
    *X = I_TET_X [order];
    *Y = I_TET_Y [order];
    *Z = I_TET_Z [order];
    *W = I_TET_W [order];
     N = I_TET_N [order];
  }
  break;
  case 5:
  case 13:
  {
    *X = I_PYR_X [order];
    *Y = I_PYR_Y [order];
    *Z = I_PYR_Z [order];
    *W = I_PYR_W [order];
     N = I_PYR_N [order];
  }
  break;
  case 6:
  case 15:
  {
    *X = I_WED_X [order];
    *Y = I_WED_Y [order];
    *Z = I_WED_Z [order];
    *W = I_WED_W [order];
     N = I_WED_N [order];
  }
  break;
  case 8:
  case 20:
  {
    *X = I_HEX_X [order];
    *Y = I_HEX_Y [order];
    *Z = I_HEX_Z [order];
    *W = I_HEX_W [order];
     N = I_HEX_N [order];
  }
  break;
  }

  return N;
}

/* predict 3D integration order */
typedef enum {MASS, BODF, INTF} ENTITY3D;
static int integrator3d_order (int type, ENTITY3D entity)
{
  /* FIXME */

  return 2;
}

/* load 2D integrator data */
inline static int integrator2d_load (int type, int order, const double **X, const double **Y, const double **W)
{
  int N;

  ASSERT_DEBUG (order >= 1 && order < MAX_ORDER, "Integration order out of bounds");

  switch (type)
  {
  case 3:
  case 6:
  {
    *X = I_TRI_X [order];
    *Y = I_TRI_Y [order];
    *W = I_TRI_W [order];
     N = I_TRI_N [order];
  }
  break;
  case 4:
  case 8:
  {
    *X = I_QUA_X [order];
    *Y = I_QUA_Y [order];
    *W = I_QUA_W [order];
     N = I_QUA_N [order];
  }
  break;
  }

  return N;
}

/* predict 2D integration order */
typedef enum {SURF, DROT} ENTITY2D;
static int integrator2d_order (int type, ENTITY2D entity)
{
  /* FIXME */

  return 1;
}

/* element integration */
#define INTEGRATE3D(TYPE, ENTITY, DOM, DOMNUM, ...)\
{\
  const double *__X__, *__Y__, *__Z__, *__W__;\
  double point [3], weight;\
  int __N__, __k__, __l__;\
  TRISURF *__d__ = DOM;\
\
  if (__d__)\
  {\
    double __subnodes__ [4][3], __subJ__;\
    double __subpoint__ [3];\
    int __domnum__ = DOMNUM;\
    TRI *__t__, *__e__;\
    double __volume__;\
\
    __N__ = integrator3d_load (4, integrator3d_order (TYPE, ENTITY), &__X__, &__Y__, &__Z__, &__W__);\
\
    for (__volume__ = __k__ = 0; __k__ < __domnum__; __k__ ++) __volume__ += __d__ [__k__].volume;\
\
    for (__l__ = 0; __l__ < __domnum__; __d__ ++, __l__ ++)\
    {\
      if (__d__->volume > DOM_TOL * __volume__)\
      {\
	COPY (__d__->center, __subnodes__ [3]);\
\
	for (__t__ = __d__->tri, __e__ = __t__ + __d__->m; __t__ < __e__; __t__ ++)\
	{\
	  COPY (__t__->ver [0], __subnodes__ [0]);\
	  COPY (__t__->ver [1], __subnodes__ [1]);\
	  COPY (__t__->ver [2], __subnodes__ [2]);\
\
	  for (__k__ = 0; __k__ < __N__; __k__ ++)\
	  {\
	    __subpoint__ [0] = __X__ [__k__];\
	    __subpoint__ [1] = __Y__ [__k__];\
	    __subpoint__ [2] = __Z__ [__k__];\
	    __subJ__ = element_det (4, __subnodes__, __subpoint__, NULL);\
	    tet_o1_local_to_global (__subnodes__, __subpoint__, point);\
	    weight = __subJ__ * __W__ [__k__];\
\
	    __VA_ARGS__\
	  }\
	}\
      }\
      else\
      {\
	COPY (__d__->center, point);\
	weight = __d__->volume;\
\
	__VA_ARGS__\
      }\
    }\
  }\
  else\
  {\
    __N__ = integrator3d_load (TYPE, integrator3d_order (TYPE, ENTITY), &__X__, &__Y__, &__Z__, &__W__);\
    for (__k__ = 0; __k__ < __N__; __k__ ++)\
    {\
      point [0] = __X__ [__k__];\
      point [1] = __Y__ [__k__];\
      point [2] = __Z__ [__k__];\
      weight = __W__ [__k__];\
\
      __VA_ARGS__\
    }\
  }\
}

/* face integral begins */
#define INTEGRAL2D_BEGIN(TYPE, ENTITY)\
{\
  const double *__X__, *__Y__, *__W__;\
  double point [2], weight;\
  int __N__, __k__;\
\
  __N__ = integrator2d_load (TYPE, integrator2d_order (TYPE, ENTITY), &__X__, &__Y__, &__W__);\
  for (__k__ = 0; __k__ < __N__; __k__ ++)\
  {\
    point [0] = __X__ [__k__];\
    point [1] = __Y__ [__k__];\
    weight = __W__ [__k__];\

/* integration code in between
 * uses point and weight values */

/* face integral ends */
#define INTEGRAL2D_END()\
  }\
}

/* ==================== FACE ======================= */

/* face shape functions at a local 2-point */
inline static int face_shapes (FACE *fac, double *point, double *shapes)
{
  switch (fac->type)
  {
  case 3: /* 1st order triangle */
    shapes [0] = point [0];
    shapes [1] = point [1];
    shapes [2] = 1.0 - point [0] - point [1];
    return 3;
  case 4: /* 1st order quadrilateral */
    shapes [0] = 0.25 * (1.0 - point [0]) * (1.0 - point [1]);
    shapes [1] = 0.25 * (1.0 + point [0]) * (1.0 - point [1]);
    shapes [2] = 0.25 * (1.0 + point [0]) * (1.0 + point [1]);
    shapes [3] = 0.25 * (1.0 - point [0]) * (1.0 + point [1]);
    return 4;
  }

  return 0;
}

/* face shape derivatices at a local 2-point */
inline static int face_derivs (FACE *fac, double *point, double *derivs)
{
  switch (fac->type)
  {
  case 3:
    derivs [0] =  1.0;
    derivs [1] =  0.0;
    derivs [2] =  0.0;
    derivs [3] =  1.0;
    derivs [4] = -1.0;
    derivs [5] = -1.0;
    return 3;
  case 4:
    derivs [0] =  0.25 * (point [1] - 1.0);
    derivs [1] =  0.25 * (point [0] - 1.0);
    derivs [2] =  0.25 * (1.0 - point [1]);
    derivs [3] = -0.25 * (1.0 + point [0]);
    derivs [4] =  0.25 * (1.0 + point [1]);
    derivs [5] =  0.25 * (1.0 + point [0]);
    derivs [6] = -0.25 * (1.0 + point [1]);
    derivs [7] =  0.25 * (1.0 - point [0]);
    return 4;
  }

  return 0;
}

/* compute face coordinates transformation determinant at a local 2-point */
inline static double face_det (FACE *fac, node_t nodes, double *point)
{
  double derivs [16], d0 [3], d1[3], normal [3], *d;
  int i, n;

  n = face_derivs (fac, point, derivs);

  SET (d0, 0);
  SET (d1, 0);

  for (i = 0, d = derivs; i < n; i ++, d += 2)
  {
    ADDMUL (d0, d[0], nodes[i], d0);
    ADDMUL (d1, d[1], nodes[i], d1);
  }

  PRODUCT (d0, d1, normal);

  return LEN (normal);
}

/* load face node coordinates into a local table */
#define face_nodes(heap, type, nodes, stack) element_nodes (heap, type, nodes, stack)

/* load face displacemnts into a local table */
inline static void face_displacements (double *heap, FACE *fac, double (*q) [3])
{
  double *p;
  int i;

  for (i = 0; i < fac->type; i ++)
  {
    p = &heap [3 * fac->nodes [i]];
    COPY (p, q[i]);
  }
}

/* =============================== ELEMENT ================================ */

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

/* linear pyramid shape functions */
inline static void pyr_o1_shapes (double *point, double *shapes)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* linear wedge shape functions */
inline static void wed_o1_shapes (double *point, double *shapes)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* linear tetrahedron shape derivatives */
inline static void tet_o1_derivs (double *point, double *derivs)
{
  derivs [0] = derivs [4] = derivs [8] = 1.0;
  derivs [1] = derivs [2] = derivs [3] = derivs [5] = derivs [6] = derivs [7] = 0.0;
  derivs [9] = derivs [10] = derivs [11] = -1.0;
}

/* linear hexahedron shape derivatives */
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

  derivs[9]  = -0.125 * (1 + point[1]) * (1 - point[2]);
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

/* linear pyramid shape derivatives  */
inline static void pyr_o1_derivs (double *point, double *derivs)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* linear wedge shape derivatives  */
inline static void wed_o1_derivs (double *point, double *derivs)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* linear tetrahedron local to global point transformation */
static void tet_o1_local_to_global (node_t nodes, double *local, double *global)
{
  double shapes [4];
  int i, j;

  tet_o1_shapes (local, shapes);

  SET (global, 0.0);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 4; j ++)
      global [i] += nodes [j][i] * shapes [j];
}

/* finite element shape functions; returns the number of nodes */
inline static int element_shapes (int type, double *point, double *shapes)
{
  switch (type)
  {
  case 4: tet_o1_shapes (point, shapes); return 4;
  case 5: pyr_o1_shapes (point, shapes); return 5;
  case 6: wed_o1_shapes (point, shapes); return 6;
  case 8: hex_o1_shapes (point, shapes); return 8;
  }

  return 0;
}

/* finite element shape derivatives; returns the number of nodes */
inline static int element_derivs (int type, double *point, double *derivs)
{
  switch (type)
  {
  case 4: tet_o1_derivs (point, derivs); return 4;
  case 5: pyr_o1_derivs (point, derivs); return 5;
  case 6: wed_o1_derivs (point, derivs); return 6;
  case 8: hex_o1_derivs (point, derivs); return 8;
  }

  return 0;
}

/* copy element node coordinates into a local table */
inline static int element_nodes (node_t heap, int type, int *nodes, node_t stack)
{
  int i;

  for (i = 0; i < type; i ++)
  { 
    COPY (heap [nodes [i]], stack [i]);
  }

  return i;
}

/* copy element node displacemnts into a local table */
static int element_displacements (double *heap, ELEMENT *ele, double (*q) [3])
{
  double *p;
  int i;

  for (i = 0; i < ele->type; i ++)
  {
    p = &heap [3 * ele->nodes [i]];
    COPY (p, q[i]);
  }

  return i;
}

/* copy element node velocities into a local table */
static int element_velocities (BODY *bod, ELEMENT *ele, double *velo, double (*u) [3])
{
  double *p;
  int i, j;

  for (i = 0; i < ele->type; i ++)
  {
    j = 3 * ele->nodes [i];
    p = &velo [j];
    COPY (p, u[i]);
  }

  return i;
}

/* copy element node numbers into a local table */
static int element_node_numbers (ELEMENT *ele, int *numbers)
{
  int i;

  for (i = 0; i < ele->type; i ++) numbers [i] = ele->nodes [i];

  return i;
}
 
/* local to global point transformation */
static void local_to_global (node_t heap, ELEMENT *ele, double *local, double *global)
{
  double nodes [MAX_NODES][3], shapes [MAX_NODES];
  int i, j, n;

  n = element_nodes (heap, ele->type, ele->nodes, nodes);

  element_shapes (ele->type, local, shapes);

  SET (global, 0.0);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < n; j ++)
      global [i] += nodes [j][i] * shapes [j];
}

/* global to local point transformation */
static void global_to_local (node_t heap, ELEMENT *ele, double *global, double *local)
{
  double nodes [MAX_NODES][3], derivs [3*MAX_NODES], shapes [MAX_NODES];
  double A [9], B [3], I [9], det, error;
  int i, j, k, l, n;

  element_nodes (heap, ele->type, ele->nodes, nodes);

  SET (local, 0.0);
  l = 0;

  do
  {
    n = element_shapes (ele->type, local, shapes);

    COPY (global, B);

    for (i = 0; i < 3; i ++)
      for (j = 0; j < n; j ++)
	B [i] -= nodes [j][i] * shapes [j];

    element_derivs (ele->type, local, derivs);

    SET9 (A, 0.0);

    for (i = 0; i < 3; i ++)
      for (j = 0; j < 3; j ++)
	for (k = 0; k < n; k ++) A [3*j+i] += nodes[k][i] * derivs [3*k+j];

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

/* compute local coordinates of a spatial point */
#define spatial_to_local(msh, ele, x, point) global_to_local ((msh)->cur_nodes, ele, x, point)

/* compute local coordinates of a referential point */
#define referential_to_local(msh, ele, X, point) global_to_local ((msh)->ref_nodes, ele, X, point)

/* compute spatial coordinates of a local point */
#define local_to_spatial(msh, ele, x, point) local_to_global ((msh)->cur_nodes, ele, x, point)

/* compute referential coordinates of a local point */
#define local_to_referential(msh, ele, X, point) local_to_global ((msh)->ref_nodes, ele, X, point)

/* coordinates transformation determinant at local point */
inline static double element_det (int type, node_t nodes, double *point, double *F)
{
  double derivs [3*MAX_NODES], G [9];
  int i, j, k, n;

  n = element_derivs (type, point, derivs);

  if (!F) F = G;

  SET9 (F, 0.0);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 3; j ++)
      for (k = 0; k < n; k ++) F [3*j+i] += nodes[k][i] * derivs [3*k+j];

  return DET (F);
}

/* element deformation determinant at local point */
inline static void element_gradient (int type, node_t q, double *point, double *F0, double *derivs, double *F)
{
  double local_derivs [3*MAX_NODES], IF0 [9], det, *l, *d;
  int i, j, k, n;

  n = element_derivs (type, point, local_derivs);

  INVERT (F0, IF0, det);
  ASSERT (det > 0.0, ERR_FEM_COORDS_INVERT);
  for (k = 0, l = local_derivs, d = derivs; k < n; k ++, l += 3, d += 3) { TVMUL (IF0, l, d); }

  IDENTITY (F);

  for (i = 0; i < 3; i ++)
    for (j = 0; j < 3; j ++)
      for (k = 0; k < n; k ++) F [3*j+i] += q[k][i] * derivs [3*k+j];
}

/* compute element shape functions at a local point and return global matrix */
static MX* element_shapes_matrix (BODY *bod, MESH *msh, ELEMENT *ele, double *point)
{
  static int *p, *i, *q, *u, k, n, m, o;
  double shapes [MAX_NODES], *x, *y;
  MX *N;

  o = element_shapes (ele->type, point, shapes);

  n = bod->dofs;

  ERRMEM (p = MEM_CALLOC ((2*n+1) * sizeof (int)));
  i = p + n + 1;

  for (k = 0, u = i; k < o; k ++, u += 3)
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

  for (k = 0; k < o; k ++)
  {
    m = ele->nodes [k] * 3;
    y = &x [p [m]];
    SET (y, shapes [k]);
  }

  return N;
}

/* lump contribution of the element mass into the global diagonal matrix */
static void element_lump_mass (BODY *bod, MESH *msh, ELEMENT *ele, double *x)
{
  double J, coef, integral, density,
         nodes [MAX_NODES][3],
	 shapes [MAX_NODES],
	*out [MAX_NODES];
  int i, j, n;

  density = ele->mat ? ele->mat->density : bod->mat->density;

  n = element_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  for (i = 0; i < n; i ++) out [i] = &x [3 * ele->nodes [i]];

  INTEGRATE3D (ele->type, MASS, ele->dom, ele->domnum,

    element_shapes (ele->type, point, shapes);
    J = element_det (ele->type, nodes, point, NULL);
    coef = density * J * weight;

    for (i = 0; i < n; i ++)
    {
      for (j = 0; j < n; j ++)
      {
	integral = coef * shapes [i] * shapes [j];

	out [i][0] += integral;
	out [i][1] += integral;
	out [i][2] += integral;
      }
    }
  )
}

/* copute element body force contribution */
static void element_body_force (BODY *bod, MESH *msh, ELEMENT *ele, double *f, double *g)
{
  double J, coef, integral, density,
         nodes [MAX_NODES][3],
         shapes [MAX_NODES];
  int i, n;

  density = ele->mat ? ele->mat->density : bod->mat->density;

  n = element_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  for (i = 0; i < 3*n; i ++) g [i] = 0.0;

  INTEGRATE3D (ele->type, BODF, ele->dom, ele->domnum, 

    element_shapes (ele->type, point, shapes);
    J = element_det (ele->type, nodes, point, NULL);
    coef = density * J * weight;

    for (i = 0; i < n; i ++)
    {
      integral = coef * shapes [i];

      g [3*i+0] += f [0] * integral;
      g [3*i+1] += f [1] * integral;
      g [3*i+2] += f [2] * integral;
    }
  )
}

/* copute element internal force or force derivative contribution */
static void element_internal_force (int derivative, BODY *bod, MESH *msh, ELEMENT *ele, double *g)
{
  double nodes [MAX_NODES][3], q [MAX_NODES][3], derivs [3*MAX_NODES],
	 F0 [9], F [9], P [9], K [81], KB [9], J, integral, *B, *p;
  BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
  double mat_lambda, mat_mi;
  int i, j, n, m;

  n = element_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  m = 3 * n;

  for (i = 0; i < n; i ++)
  {
    p = &bod->conf [3 * ele->nodes [i]];
    COPY (p, q[i]);
  }

  for (i = 0, j = m * (derivative ? m : 1); i < j; i ++) g [i] = 0.0;
  mat_lambda = lambda (mat->young, mat->poisson);
  mat_mi  = mi (mat->young, mat->poisson);

  INTEGRATE3D (ele->type, INTF, ele->dom, ele->domnum,

    J = element_det (ele->type, nodes, point, F0);
    element_gradient (ele->type, q, point, F0, derivs, F);
    integral = J * weight;

    if (derivative)
    {
      SVK_Tangent_C (mat_lambda, mat_mi, integral, 9, F, K); /* TODO: generalize in BULK_MATERIAL interface */

      for (i = 0; i < m; i ++) /* see doc/notes.lyx for details */
      {
	SET9 (KB, 0);
	for (j = 0; j < 3; j ++)
	{
	  p = &K [9*((i%3) + (3*j))];
	  integral = derivs [3*(i/3)+j];
	  NNADDMUL (KB, integral, p, KB);
	}

	for (j = 0, B = derivs, p = &g[m*i]; j < n; j ++, B += 3, p += 3) { NVADDMUL (p, KB, B, p); }
      }
    }
    else
    {
      SVK_Stress_C (mat_lambda, mat_mi, integral, F, P); /* TODO: generalize in BULK_MATERIAL interface */

      for (i = 0, B = derivs, p = g; i < n; i ++, B += 3, p += 3) { NVADDMUL (p, P, B, p); }
    }
  )
}

/* compute inv (M) * K for an element */
static MX* element_inv_M_K (BODY *bod, MESH *msh, ELEMENT *ele)
{
  double mass [24];
  double *x, *y;
  MX *IMK, *IM;
  int i, j, n;

  n = ele->type * 3;
  IM = bod->M;
  IMK = MX_Create (MXDENSE, n, n, NULL, NULL);
  element_internal_force (1, bod, msh, ele, IMK->x);

  for (i = 0; i < ele->type; i ++)
  {
    for (j = 0; j < 3; j ++)
    {
      mass [3*i+j] = 1.0 / IM->x [ele->nodes [i] * 3 + j];
    }
  }

  for (j = 0, x = IMK->x; j < n; j ++) /* compute IMK = IM * K */
  {
    for (i = 0, y = mass; i < n; i ++, x ++, y ++) (*x) *= (*y); /* scale each column by diagonal IM entries */
  }

  return IMK;
}

/* =========================== GENERAL ================================== */

/* compute deformation gradient at a local point */
static void deformation_gradient (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *F)
{
  double derivs [3*MAX_NODES], nodes [MAX_NODES][3], q [MAX_NODES][3], F0 [9], *p;
  int i, n;

  n = element_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  for (i = 0; i < n; i ++)
  {
    p = &bod->conf [3 * ele->nodes [i]];
    COPY (p, q[i]);
  }

  element_det (ele->type, nodes, point, F0);
  element_gradient (ele->type, q, point, F0, derivs, F);
}

/* compute Cauchy stress at a local point */
static void cauchy_stress (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *values)
{
  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
    {
      BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
      double P [9], J, F [9];

      deformation_gradient (bod, msh, ele, point, F);

      J = SVK_Stress_C (lambda (mat->young, mat->poisson), mi (mat->young, mat->poisson), 1.0, F, P);  /* TODO: generalize in BULK_MATERIAL interface */

      values [0] = (F[0]*P[0]+F[3]*P[1]+F[6]*P[2])/J; /* sx  */
      values [1] = (F[1]*P[3]+F[4]*P[4]+F[7]*P[5])/J; /* sy  */
      values [2] = (F[2]*P[6]+F[5]*P[7]+F[8]*P[8])/J; /* sz  */
      values [3] = (F[0]*P[3]+F[3]*P[4]+F[6]*P[5])/J; /* sxy */
      values [4] = (F[0]*P[6]+F[3]*P[7]+F[6]*P[8])/J; /* sxz */
      values [5] = (F[2]*P[0]+F[5]*P[1]+F[8]*P[2])/J; /* syz */
    }
    break;
    case BODY_COROTATIONAL:
    {
      /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
    }
    break;
  }
}

/* accumulate point force contribution into the body force vector */
static void accumulate_point_force (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *f, double *force)
{
  MX *N = element_shapes_matrix (bod, msh, ele, point);
  MX_Matvec (1.0, MX_Tran (N), f, 1.0, force);
  MX_Destroy (N);
}

/* return element stabbed by a spatial point */
static ELEMENT* stabbed_spatial_element (MESH *msh, ELEMENT **ele, int nele, double *x)
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
static ELEMENT* stabbed_referential_element (MESH *msh, ELEMENT **ele, int nele, double *X)
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

/* accumulate constraints reaction */
inline static void accumulate_reac (BODY *bod, MESH *msh, ELEMENT *ele, double *X, double *base, double *R, short isma, double *force)
{
  double shapes [MAX_NODES], P [3], point [3], *f;
  int numbers [MAX_NODES], i, n;

  referential_to_local (msh, ele, X, point);
  n = element_shapes (ele->type, point, shapes);
  element_node_numbers (ele, numbers);

  NVMUL (base, R, P);

  if (isma) for (i = 0; i < n; i ++) { f = &force [3 * numbers [i]]; ADDMUL (f, shapes [i], P, f); }
  else for (i = 0; i < n; i ++) { f = &force [3 * numbers [i]]; SUBMUL (f, shapes [i], P, f); }
}

/* compute constraints force r = SUM H' R */
static void fem_constraints_force (BODY *bod, double *r)
{
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  SET *node;
  int i;

  msh = FEM_MESH (bod);

  for (i = 0; i < bod->dofs; i ++) r [i] = 0.0;

  for (node = SET_First (bod->con); node; node = SET_Next (node))
  {
    CON *con = node->data;
    short isma = (bod == con->master);
    double *X = (isma ? con->mpnt : con->spnt);

    if (bod->msh)
    {
      cvx = (isma ? mgobj(con) : sgobj(con));
      ele = stabbed_referential_element (msh, cvx->ele, cvx->nele, X); /* TODO: optimize */
    }
    else ele = (isma ? mgobj(con) : sgobj(con));

    accumulate_reac (bod, msh, ele, X, con->base, con->R, isma, r);

    if (isma && bod == con->slave) /* self-contact */
    {
      X = con->spnt;

      if (bod->msh)
      {
	cvx = sgobj(con);
	ele = stabbed_referential_element (msh, cvx->ele, cvx->nele, X); /* TODO: optimize */
      }
      else ele = sgobj(con);

      accumulate_reac (bod, msh, ele, X, con->base, con->R, 0, r);
    }
  }
}

/* compute inernal force */
static void internal_force (BODY *bod, double *fint)
{
  MESH *msh = FEM_MESH (bod);
  double g [24], *v, *w;
  ELEMENT *ele;
  int bulk, i;

  for (i = 0; i < bod->dofs; i ++) fint [i] = 0.0;

  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    element_internal_force (0, bod, msh, ele, g);

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

/* initialize unit body force */
static void unit_body_force (BODY *bod)
{
  double g [24], f [3], *v, *w;
  double *ubf = FEM_FBOD (bod);
  MESH *msh = FEM_MESH (bod);
  ELEMENT *ele;
  int bulk, i;

  SET (f, 1.0);
  blas_dscal (bod->dofs, 0.0, ubf, 1);
  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    element_body_force (bod, msh, ele, f, g);

    for (i = 0, v = g; i < ele->type; i ++, v += 3)
    {
      w = &ubf [ele->nodes [i] * 3];
      ADD (w, v, w);
    }

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }
}

/* compute external force */
static void external_force (BODY *bod, double time, double step, double *fext)
{
  double g [24], f [3], point [3], value, *v;
  double *ubf = FEM_FBOD (bod);
  MESH *msh = FEM_MESH (bod);
  ELEMENT *ele;
  int i;

  /* zero forces */
  for (i = 0; i < bod->dofs; i ++) fext [i] = 0.0;

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

      ele = MESH_Element_Containing_Point (msh, frc->ref_point, 1); /* TODO: optimize */

      if (ele)
      {
	referential_to_local (msh, ele, frc->ref_point, point);

	if (frc->kind & CONVECTED)
	{ 
	  deformation_gradient (bod, msh, ele, point, g);
	  NVMUL (g, f, g+9);
	  COPY (g+9, f);
	}

	accumulate_point_force (bod, msh, ele, point, f, fext);
      }
#if DEBUG
      else WARNING_DEBUG (0, "Point force outside of FE mesh");
#endif
    }
  }

  /* gravitation */
  if (bod->dom->gravity [0])
  {
    f [0] = TMS_Value (bod->dom->gravity [0], time);
    f [1] = TMS_Value (bod->dom->gravity [1], time);
    f [2] = TMS_Value (bod->dom->gravity [2], time);

    i = bod->dofs / 3;
    blas_daxpy (i, f[0], ubf, 3, fext, 3);
    blas_daxpy (i, f[1], ubf+1, 3, fext+1, 3);
    blas_daxpy (i, f[2], ubf+2, 3, fext+2, 3);
  }
}
 
/* compute out of balance force = fext - fint */
static void dynamic_force (BODY *bod, double time, double step, double *fext, double *fint, double *force)
{

  external_force (bod, time, step, fext);

  internal_force (bod, fint);

  for (double *x = fext, *y = fint, *z = force, *u = z + bod->dofs; z < u; x ++, y ++, z ++) *z = (*x) - (*y);
}

/* the smame computation for the static case */
#define static_force(bod, time, step, fext, fint, force) dynamic_force (bod,time,step,fext,fint,force)

/* compute global tangent stiffness */
static MX* tangent_stiffness (BODY *bod, short spd)
{
  struct colblock { int row; double val [3]; } *cb; /* column block */
  int i, j, k, l, n, *pp, *ii, *kk;
  double K [576], *A;
  MAP **col, *item;
  ELEMENT *ele;
  short bulk;
  MESH *msh;
  MEM blkmem,
      mapmem;
  MX *tang;

  if (spd) spd = 1;
  MEM_Init  (&blkmem, sizeof (struct colblock), bod->dofs);
  ERRMEM (col = MEM_CALLOC (sizeof (MAP*) * bod->dofs)); /* sparse columns */
  MEM_Init  (&mapmem, sizeof (MAP), bod->dofs);
  msh = FEM_MESH (bod);

  for (ele = msh->surfeles, bulk = 0; ele;
       ele = (ele->next ? ele->next : bulk ? NULL : msh->bulkeles),
       bulk = (ele == msh->bulkeles ? 1 : bulk)) /* for each element in mesh */
  {
    element_internal_force (1, bod, msh, ele, K); /* compute internal force derivartive: K */

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
	  if (spd && ele->nodes [k] > ele->nodes [i]) continue; /* lower triangle */

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

  for (pp [0] = 0, j = 0; j < bod->dofs; j ++) pp [j+1] = pp [j] + 3 * MAP_Size (col [j]) - spd * (j % 3); /* subtract upper triangular j % 3 sticking out bits */

  ERRMEM (ii = malloc (sizeof (int [pp [bod->dofs]]))); /* row indices */

  for (j = 0, kk = ii; j < bod->dofs; j ++) /* initialize row index pointer; for each column */
  {
    for (item = MAP_First (col [j]); item; item = MAP_Next (item)) /* for each row-block */
    {
      cb = item->data;
      if (spd && cb->row < j) /* diagonal block with sticking out upper triangle */
      {
	for (n = 0; n < 3 - j % 3; n ++, kk ++) kk [0] = j + n;
      }
      else /* lower triangle block */
      {
	kk [0] = cb->row;
	kk [1] = kk[0] + 1;
	kk [2] = kk[1] + 1;
	kk += 3;
      }
    }
  }

  tang = MX_Create (MXCSC, bod->dofs, bod->dofs, pp, ii); /* create tangent matrix structure */
  if (spd) tang->flags |= MXSPD;

  for (j = 0, A = tang->x; j < bod->dofs; j ++) /* initialize column values pointer A; for each column */
  {
    for (item = MAP_First (col [j]); item; item = MAP_Next (item)) /* for each row-block */
    {
      cb = item->data;
      if (spd && cb->row < j) /* diagonal block with sticking out upper triangle */
      {
	for (n = 0; n < 3 - j % 3; n ++, A ++) A [0] = cb->val [j % 3 + n];
      }
      else /* lower triangle block */
      {
        COPY (cb->val, A); /* copy values */
	A += 3;
      }
    }
  }

  free (ii);
  free (pp);
  free (col);
  MEM_Release (&mapmem);
  MEM_Release (&blkmem);

  return tang;
}

/* compute diagonalized inertia operator */
static MX* diagonal_inertia (BODY *bod, short spd)
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
  if (spd) M->flags |= MXSPD;
  x = M->x;
  free (p);
  free (i);

  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    element_lump_mass (bod, msh, ele, x);

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }

  return M; 
}

/* =================== TOAL LAGRANGIAN =================== */

/* compute inverse operator for the implicit dynamic time stepping */
static void TL_dynamic_inverse (BODY *bod, double step, double *force)
{
  MX *M, *K;

  if (bod->inverse) MX_Destroy (bod->inverse);

  K = tangent_stiffness (bod, 0);

  M = bod->M;

  if (force)
  {
    /* account for the previous velocity */
    MX_Matvec (1.0 / step, M, bod->velo, 1.0, force);

    /* account for the internal force increment */
    MX_Matvec (-0.25 * step, K, bod->velo, 1.0, force);
  }

  /* calculate tangent operator A = M + h*h/4 K */
  bod->inverse = MX_Add (1.0, M, 0.25*step*step, K, NULL);

  /* invert A */
  MX_Inverse (bod->inverse, bod->inverse);

  /* clean up */
  MX_Destroy (K);
}

/* solve implicit ingetration nonlinear equations */
static void TL_dynamic_solve (BODY *bod, double time, double step, double *fext, double *f, short begin)
{
  int n = bod->dofs,
      iter,
      i;

  double half = 0.5 * step,
	 quad = 0.25 * step,
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
    dynamic_force (bod, time+half, step, fext, fint, f);  /* f(t+h/2,q(t+h/2)) = fext (t+h/2) - fint (q(t+h/2)) */
    TL_dynamic_inverse (bod, step, NULL); /* A = M + (h*h/4) * K */
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
    for (i = 0; i < n; i ++) q [i] = qorig [i] + quad * (u0[i] + u[i]); /* overwrite bod->conf ... */
    internal_force (bod, fint); /* ... as it is used in there ... */
    TL_dynamic_inverse (bod, step, NULL); /* ... and there */
    for (i = 0; i < n; i ++) { res [i] = step * (fext [i] - fint [i]); aux [i] = u [i] - u0[i]; }
    MX_Matvec (-1.0, bod->M, aux, 1.0, res);
    MX_Matvec (1.0, bod->inverse, res, 0.0, aux);
    for (i = 0; i < n; i ++) u [i] += aux [i];
    errlo = blas_ddot (n, u, 1, u, 1);
    errup = blas_ddot (n, aux, 1, aux, 1);
    error = sqrt (errup / MAX (errlo, 1.0));
  }
  while (error > IMP_EPS && ++ iter < MAX_ITERS);

#if 0
  printf ("DEF_IMP: iter = %d, error = %e\n", iter, error);
#endif

  ASSERT (iter < MAX_ITERS, ERR_BOD_SCHEME_NOT_CONVERGED);

  if (begin)
  {
    for (i = 0; i < n; i ++) q [i] = qorig [i] + half * u0[i]; /* q(t+h/2) = q(t) + (h/2) * u(t) */
  }
  else
  {
    for (i = 0; i < n; i ++) q [i] = qorig [i] + half * (u0[i] + u [i]); /* q(t+h/2) = q(t) + (h/2) * (u(t+h) + u(t)) */
  }

  free (qorig);
}

/* static time-stepping inverse */
static void TL_static_inverse (BODY *bod, double step)
{
  MX *M, *K, *A;

  if (bod->M) M = bod->M; else bod->M = M = diagonal_inertia (bod, 0);

  if (bod->inverse) MX_Destroy (bod->inverse);

  K = tangent_stiffness (bod, 0);

#if 0
  MESH *msh = FEM_MESH (bod);
  double eigmax;
  ELEMENT *ele;
  MX *IMK;

  /* estimate maximal eigenvalue of inv (M) * K based on just one element */
  ele = msh->surfeles;
  IMK = element_inv_M_K (bod, msh, ele); /* element inv (M) * K */
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

/* total lagrangian initialise dynamic time stepping */
static void TL_dynamic_init (BODY *bod)
{
  if (!bod->M) bod->M = diagonal_inertia (bod, 0);

  if (bod->scheme == SCH_DEF_EXP)
  {
    if (!bod->inverse) /* initialize once */
    {
      double *x, *y;

      bod->inverse = MX_Copy (bod->M, NULL);

      for (x = bod->inverse->x, y = x + bod->dofs; x < y; x ++)
      {
	ASSERT (*x > 0.0, ERR_FEM_MASS_NOT_SPD);
	(*x) = 1.0 / (*x); /* invert diagonal */
      }
    }
  }
  else TL_dynamic_inverse (bod, bod->dom->step, NULL); /* update every time */
}

/* total lagrangian estimate critical step for the dynamic scheme */
static double TL_dynamic_critical_step (BODY *bod)
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
      IMK = element_inv_M_K (bod, msh, ele); /* element inv (M) * K */
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

/* total lagrangian perform the initial half-step of the dynamic scheme */
static void TL_dynamic_step_begin (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  double half = 0.5 * step,
	 c = bod->damping,
	*x = bod->inverse->x,
	*u0 = FEM_VEL0 (bod),
	*fext = FEM_FEXT (bod),
	*fint = FEM_FINT (bod),
	*q = bod->conf,
	*u = bod->velo,
	*e = u + n,
	*f, *g;

  ERRMEM (f = malloc (sizeof (double [n])));

  blas_dcopy (n, u, 1, u0, 1); /* save u (t) */

  switch (bod->scheme)
  {
  case SCH_DEF_EXP:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    dynamic_force (bod, time+half, step, fext, fint, f);  /* f = fext (t+h/2) - fint (q(t+h/2)) */
    for (g = f; u < e; u ++, x ++, g ++) (*u) += step * (*x) * (*g); /* u(t+h) = u(t) + inv (M) * h * f */
  }
  break;
  case SCH_DEF_LIM:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    dynamic_force (bod, time+half, step, fext, fint, f);  /* f = fext (t+h/2) - fint (q(t+h/2)) */
    TL_dynamic_inverse (bod, step, NULL); /* A = M + (h*h/4) * K */
    MX_Matvec (step, bod->inverse, f, 1.0, u); /* u(t+h) = u(t) + inv (A) * h * f */
  }
  break;
  case SCH_DEF_LIM2:
  {
    dynamic_force (bod, time+half, step, fext, fint, f);  /* f = fext (t+h/2) - fint (q(t)) */
    TL_dynamic_inverse (bod, step, f); /* f += (1/h) M u(t) - (h/4) K u (t) */
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
    TL_dynamic_solve (bod, time, step, fext, f, 1);
  }
  break;
  default:
  break;
  }

  if (c > 0.0) for (u = bod->velo; u < e; u ++, u0++) (*u) -= c * (*u0); /* u(t+h) -= c * u (t) */

  free (f);
}

/* total lagrangian perform the final half-step of the dynamic scheme */
static void TL_dynamic_step_end (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  double half = 0.5 * step,
	*x = bod->inverse->x,
	*fext = FEM_FEXT (bod),
	*u = bod->velo,
	*iu = u,
	*q = bod->conf,
	*e = u + n,
	*r, *ir;

  ERRMEM (ir = r = malloc (sizeof (double [n])));

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
    TL_dynamic_solve (bod, time, step, fext, r, 0);
  }
  break;
  default:
  break;
  }

  free (r);
}

/* total lagrangian initialise static time stepping */
static void TL_static_init (BODY *bod)
{
  TL_static_inverse (bod, bod->dom->step);
}

/* total lagrangian perform the initial half-step of the static scheme */
static void TL_static_step_begin (BODY *bod, double time, double step)
{
  double *f;

  ERRMEM (f = malloc (sizeof (double [bod->dofs])));
  TL_static_inverse (bod, step); /* compute inverse of static tangent operator */
  static_force (bod, time+step, step, FEM_FEXT(bod), FEM_FINT(bod), f);  /* f(t+h) = fext (t+h) - fint (q(t+h)) */
  MX_Matvec (step, bod->inverse, f, 0.0, bod->velo); /* u(t+h) = inv (A) * h * f(t+h) */
  free (f);
}

/* total lagrangian perform the final half-step of the static scheme */
static void TL_static_step_end (BODY *bod, double time, double step)
{
  double *r;

  ERRMEM (r = malloc (sizeof (double [bod->dofs])));
  fem_constraints_force (bod, r); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  MX_Matvec (step, bod->inverse, r, 1.0, bod->velo); /* u(t+h) += inv (A) * h * r */
  blas_daxpy (bod->dofs, step, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t) + h * u(t+h) */
  free (r);
}

/* =================== BODY COROTATIONAL =================== */

/* compute surface integral = INT { skew [N] A [X + Shapes (X) q] } */
static void BC_surface_integral (BODY *bod, MESH *msh, double *conf, int num, double *A, double *integral)
{
  double q [8][3], nodes [8][3], B [3], C [3], shapes [8], J, coef, *N, *Y, *i, *e;
  int j, n;
  FACE *fac;

  for (i = integral, e = i + 3 * num; i < e; i += 3) SET (i, 0);

  for (fac = msh->faces; fac; fac = fac->n)
  {
    face_nodes (msh->ref_nodes, fac->type, fac->nodes, nodes);
    face_displacements (conf, fac, q);
    for (j = 0; j < fac->type; j ++) ADD (nodes [j], q [j], nodes [j]);
    N = fac->normal + 3;

    INTEGRAL2D_BEGIN (fac->type, DROT) /* defines point and weight */
    {
      n = face_shapes (fac, point, shapes);
      J = face_det (fac, nodes, point);
      coef = J * weight;

      SET (B, 0);
      for (j = 0; j < n; j ++)
	ADDMUL (B, shapes [j], nodes [j], B);
      SCALE (B, coef);

      for (i = integral, Y = A; i < e; i += 3, Y += 9)
      {
	NVMUL (Y, B, C);
	PRODUCTADD (N, C, i);
      }
    }
    INTEGRAL2D_END ()
  }
}

/* update rotation for given configuration;
 * bod - input body
 * msh - input mesh
 * q   - input configuration
 * R   - input/output rotation */
static void BC_update_rotation (BODY *bod, MESH *msh, double *q, double *R)
{
  double h [12], *dh, dJ [3], ddJ [9], O [3], A [36], dO [3][9], error;
  int iter;

  SET (O, 0);
  iter = 0;
  dh = h+3;

  do
  {
    SCALE (O, -1);
    EXPMAP (O, A+9);
    TNMUL (R, A+9, A); /* A = [exp (O) R]' */

    EXPMAP123 (O, dO[0], dO[1], dO[2]);
    SCALE9 (dO[0], -1);
    SCALE9 (dO[1], -1);
    SCALE9 (dO[2], -1);
    TNMUL (R, dO[0], A+9); /* A = R' dexp(-O) / dO */
    TNMUL (R, dO[1], A+18);
    TNMUL (R, dO[2], A+27);
    SCALE (O, -1);

    BC_surface_integral (bod, msh, q, 4, A, h); /* h = INT { skew [N] R' exp (-O) x }; dh = INT { skew [N] R' dexp (-O) / dO x } */

    TVMUL (dh, h, dJ); /* dJ = h' dh */
    TNMUL (dh, dh, ddJ); /* dJ = dh' dh */

    if (lapack_dposv ('U', 3, 1, ddJ, 3, dJ, 3) != 0)
    {
      ASSERT_DEBUG (0, "SINGULAR ddJ at iter %d", iter); /* FIXME */
    }

    SUB (O, dJ, O);
    error = sqrt (DOT (dJ, dJ) / (1.0 + DOT (O, O)));
  }
  while (++ iter < 64 && error > 1E-10);

  EXPMAP (O, A);
  NNCOPY (R, A+9);
  NNMUL (A, A+9, R); /* R = exp (O) R */

#if 0
  printf ("O = %g, %g, %g after %d iterations\n", O[0], O[1], O[2], iter);
#endif

  ASSERT_DEBUG (iter < 64, "DIVERGED rotation update"); /* FIXME */
}

/* fint = R K R' [(I-R)Z + q + (h/4) u] */
static void BC_internal_force (BODY *bod, double *R, double *q, double *fint)
{
  double *a, *b, *x, *y, *z, *d, (*Z) [3], (*e) [3], Y [3];
  MESH *msh = FEM_MESH (bod);
  MX *K = bod->K;
  int n = K->n;

  ERRMEM (a = MEM_CALLOC (2 * sizeof (double [n])));
  b = a + n;

  for (Z = msh->ref_nodes, e = Z + msh->nodes_count, d = b; Z < e; Z ++, d += 3, q += 3)
  {
    NVMUL (R, Z[0], Y);
    SUB (Z[0], Y, Y);
    ADD (Y, q, d); /* d = (I-R)Z + q */
  }

  for (x = b, y = b + n, z = a; x < y; x += 3, z += 3)
  {
    TVMUL (R, x, z); /* a = R' [(I-R)Z + q] */
  }

  MX_Matvec (1.0, K, a, 0.0, b); /* b = K a */

  for (x = b, y = b + n, z = fint; x < y; x += 3, z += 3)
  {
    NVMUL (R, x, z); /* fint = R b  */
  }

  free (a);
}

/* compute inverse of (M + (h*h/4) K) */
static void BC_inverse (BODY *bod, double step, MX *M, MX *K)
{
  if (!bod->inverse) 
  {
    if (bod->scheme == SCH_DEF_EXP)
    {
      double *x, *y;

      bod->inverse = MX_Copy (M, NULL);

      for (x = bod->inverse->x, y = x + bod->dofs; x < y; x ++)
      {
	ASSERT (*x > 0.0, ERR_FEM_MASS_NOT_SPD);
	(*x) = 1.0 / (*x); /* invert diagonal */
      }
    }
    else
    {
      bod->inverse = MX_Add (1.0, M, 0.25*step*step, K, NULL);

      MX_Inverse (bod->inverse, bod->inverse);
    }
  }
}

/* compute u = alpha * R A R' b + beta * u  */
static void BC_matvec (double alpha, MX *A, double *R, double *b, double beta, double *u)
{
  double *x, *y, *z, *w;
  int n = A->n;

  ERRMEM (x = MEM_CALLOC (2 * sizeof (double [n])));

  for (y = x, z = x + n, w = b; y < z; y += 3, w += 3)
  {
    TVMUL (R, w, y);
  }

  y = x + n;
  MX_Matvec (alpha, A, x, 0.0, y);

  for (z = y + n, w = u; y < z; y += 3, w += 3)
  {
    SCALE (w, beta);
    NVADDMUL (w, R, y, w);
  }

  free (x);
}

/* solve implicit ingetration nonlinear equations */
static void BC_dynamic_solve (BODY *bod, double time, double step, double *fext, double *f, short begin)
{
  MESH *msh = FEM_MESH (bod);

  int n = bod->dofs,
      iter,
      i;

  double half = 0.5 * step,
	 quad = 0.25 * step,
	*fint = FEM_FINT (bod),
	*u0 = FEM_VEL0 (bod),
	*R = FEM_ROT (bod),
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
    BC_update_rotation (bod, msh, q, R); /* R1 = R(q(t+h/2)) */
    external_force (bod, time+half, step, fext); /* fext = fext (t+h/2) */
    BC_internal_force (bod, R, q, fint); /* fint = R1 K R1' [(I-R1)Z + q(t+h/2)] */
    memset (f, 0, sizeof (double [n])); 
    blas_daxpy (n, step, fext, 1, f, 1);
    blas_daxpy (n, -step, fint, 1, f, 1); /* f = h (fext - fint) */
    BC_matvec (1.0, bod->inverse, R, f, 1.0, u); /* u(t+h) = u(t) + inv (M + (h*h/4) R1 K R1') f */
  }
  else
  {
    blas_dcopy (n, q, 1, qorig, 1); /* qorig = q (t+h/2) */
    blas_daxpy (n, -half, u0, 1, qorig, 1); /* qorig = q(t) = q(t+h/2) - (h/2) * u(t) */
    BC_matvec (step, bod->inverse, R, f, 1.0, u); /* u(t+h) += inv (M + (h*h/4) R1 K R1') f */
  }

  iter = 0;
  do
  {
    for (i = 0; i < n; i ++) q [i] = qorig [i] + quad * (u0[i] + u[i]); /* overwrite bod->conf */
    BC_update_rotation (bod, msh, q, R); /* R1 = R(q(t+h/2)) */
    BC_internal_force (bod, R, q, fint); /* fint = R1 K R1' [(I-R1)Z + q(t+h/2)] */
    for (i = 0; i < n; i ++) { res [i] = step * (fext [i] - fint [i]); aux [i] = u [i] - u0[i]; }
    MX_Matvec (-1.0, bod->M, aux, 1.0, res);
    BC_matvec (1.0, bod->inverse, R, res, 0.0, aux);
    for (i = 0; i < n; i ++) u [i] += aux [i];
    errlo = blas_ddot (n, u, 1, u, 1);
    errup = blas_ddot (n, aux, 1, aux, 1);
    error = sqrt (errup / MAX (errlo, 1.0));
  }
  while (error > IMP_EPS && ++ iter < MAX_ITERS);

#if 0
  printf ("DEF_IMP: iter = %d, error = %e\n", iter, error);
#endif

  ASSERT (iter < MAX_ITERS, ERR_BOD_SCHEME_NOT_CONVERGED);

  if (begin)
  {
    for (i = 0; i < n; i ++) q [i] = qorig [i] + half * u0[i]; /* q(t+h/2) = q(t) + (h/2) * u(t) */
  }
  else
  {
    for (i = 0; i < n; i ++) q [i] = qorig [i] + half * (u0[i] + u [i]); /* q(t+h/2) = q(t) + (h/2) * (u(t+h) + u(t)) */
  }

  free (qorig);
}

/* body co-rotational initialise dynamic time stepping */
static void BC_dynamic_init (BODY *bod)
{
  if (!bod->M && !bod->K)
  {
    bod->M = diagonal_inertia (bod, 1);

    bod->K = tangent_stiffness (bod, 1);

    BC_inverse (bod, bod->dom->step, bod->M, bod->K); /* initialize once */
  }
}

/* body co-rotational estimate critical step for the dynamic scheme */
static double BC_dynamic_critical_step (BODY *bod)
{
  return TL_dynamic_critical_step (bod);
}

/* body co-rotational perform the initial half-step of the dynamic scheme */
static void BC_dynamic_step_begin (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  MESH *msh = FEM_MESH (bod);
  double half = 0.5 * step,
	*u0 = FEM_VEL0 (bod),
	*R = FEM_ROT (bod),
	*fext = FEM_FEXT (bod),
	*fint = FEM_FINT (bod),
	*q = bod->conf,
	*u = bod->velo,
	*b;

  blas_dcopy (n, u, 1, u0, 1); /* save u (t) */

  ERRMEM (b = MEM_CALLOC (sizeof (double [n])));

  switch (bod->scheme)
  {
  case SCH_DEF_EXP:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    BC_update_rotation (bod, msh, q, R); /* R1 = R(q(t+h/2)) */

    external_force (bod, time+half, step, fext); /* fext = fext (t+h/2) */
    BC_internal_force (bod, R, q, fint); /* fint = R1 K R1' [(I-R1)Z + q(t+h/2)] */
    blas_daxpy (n, step, fext, 1, b, 1);
    blas_daxpy (n, -step, fint, 1, b, 1); /* b = h (fext - fint) */

    MX_Matvec (1.0, bod->inverse, b, 1.0, u); /* u(t+h) = u(t) + inv (M) * b */
  }
  break;
  case SCH_DEF_LIM:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    BC_update_rotation (bod, msh, q, R); /* R1 = R(q(t+h/2)) */

    external_force (bod, time+half, step, fext); /* fext = fext (t+h/2) */
    BC_internal_force (bod, R, q, fint); /* fint = R1 K R1' [(I-R1)Z + q(t+h/2)] */
    blas_daxpy (n, step, fext, 1, b, 1);
    blas_daxpy (n, -step, fint, 1, b, 1); /* b = h (fext - fint) */

    BC_matvec (1.0, bod->inverse, R, b, 1.0, u); /* u(t+h) = u(t) + inv (M + (h*h/4) R1 K R1') b */
  }
  break;
  case SCH_DEF_LIM2:
  {
    external_force (bod, time+half, step, fext); /* fext = fext (t+h/2) */
    BC_internal_force (bod, R, q, fint); /* fint = R K R' [(I-R)Z + q(t)] */
    blas_daxpy (n, 1.0, fext, 1, b, 1);
    blas_daxpy (n, -1.0, fint, 1, b, 1); /* b = fext - fint */
    MX_Matvec (1.0 / step, bod->M, u, 1.0, b); /* b += (1/h) M u (t) */
    BC_matvec (-0.25 * step, bod->K, R, u, 1.0, b); /* b -= (h/4) K u (t) */
    /* FIXME: dampnig effect is absent here - unlike in TL */

    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    BC_matvec (step, bod->inverse, R, b, 0.0, u); /* u(t+h) = inv (A) * h * b */
  }
  break;
  case SCH_DEF_IMP:
  { 
    /* q(t+h/2) = q(t) + (h/2) * u(t)
     * f = fext (t+h/2) - fint ([q(t) + q(t+h)]/2) 
     * A = M + (h*h/4) * K ([q(t) + q(t+h)]/2) 
     * u (t+h) = u (t) + inv (A) * h * f */
    BC_dynamic_solve (bod, time, step, fext, b, 1);
  }
  break;
  default:
  break;
  }

  free (b);
}

/* body co-rotational perform the final half-step of the dynamic scheme */
static void BC_dynamic_step_end (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  MESH *msh = FEM_MESH (bod);
  double half = 0.5 * step,
	*fext = FEM_FEXT (bod),
	*R = FEM_ROT (bod),
	*q = bod->conf,
	*u = bod->velo,
	*r;

  ERRMEM (r = malloc (sizeof (double [n])));

  fem_constraints_force (bod, r); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */

  blas_daxpy (n, 1.0, r, 1, fext, 1);  /* fext += r */

  switch (bod->scheme)
  {
  case SCH_DEF_EXP:
  {
    MX_Matvec (step, bod->inverse, r, 1.0, u); /* u(t+h) += inv (M) * h * r */

    blas_daxpy (n, half, u, 1, q, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
    BC_update_rotation (bod, msh, q, R); /* R(t+h) = R(q(t+h)) */
  }
  break;
  case SCH_DEF_LIM:
  case SCH_DEF_LIM2:
  {
    BC_matvec (step, bod->inverse, R, r, 1.0, u); /* u(t+h) += inv (A) * h * r */

    blas_daxpy (n, half, u, 1, q, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
    BC_update_rotation (bod, msh, q, R); /* R(t+h) = R(q(t+h)) */
  }
  break;
  case SCH_DEF_IMP:
  {
    /* f = fext (t+h/2) - fint ([q(t) + q(t+h)]/2) 
     * A = M + (h*h/4) * K ([q(t) + q(t+h)]/2) 
     * u (t+h) = u (t) + inv (A) * h * f
     * q(t+h) = q(t+h/2) + (h/2) * u(t+h) */
    BC_dynamic_solve (bod, time, step, fext, r, 0);
  }
  break;
  default:
  break;
  }

  free (r);
}

/* body co-rotational initialise static time stepping */
static void BC_static_init (BODY *bod)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* body co-rotational perform the initial half-step of the static scheme */
static void BC_static_step_begin (BODY *bod, double time, double step)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* body co-rotational perform the final half-step of the static scheme */
static void BC_static_step_end (BODY *bod, double time, double step)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* ================== INTERSECTIONS ==================== */

/* attach (element, local point) pairs to cvx->epn placeholder so that
 * current vertices can be updated fast by using FEM_Cur_Point_Ext */
static void post_proces_convices (SHAPE *shp, MESH *msh, FEMFORM form)
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
	epn->ele = stabbed_referential_element (msh, cvx->ele, cvx->nele, ref);
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

  /* 1. Test whether the shape_volume == volume cut of the mesh */
 
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

  /* 2. Delete elements whose volume (ele->dom) < TOL * volume (ele)
   *    or delte ele->dom when volume (ele->dom) == volume (ele) */

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

  /* 3. Delete nodes unattached to elements */

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

  /* 4. Transform global ele->dom points into local element coordinates */

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

/* ================== INTERFACE ==================== */

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
    post_proces_convices (shp, msh, form);
  }
  else msh = shp->data; /* retrive the mesh pointer from the shape */

  /* allocate dofs */
  bod->dofs = msh->nodes_count * 3;
  switch (form)
  {
    case TOTAL_LAGRANGIAN:
    {
      ERRMEM (bod->conf = MEM_CALLOC (6 * bod->dofs * sizeof (double))); /* configuration, velocity, previous velocity, fext, fint, fbod */
      bod->velo = bod->conf + bod->dofs;
    }
    break;
    case BODY_COROTATIONAL:
    {
      ERRMEM (bod->conf = MEM_CALLOC ((6 * bod->dofs + 9) * sizeof (double))); /* configuration, rotation, velocity, previous velocity, fext, fint, fbod */
      bod->velo = bod->conf + bod->dofs + 9;
      double *R = FEM_ROT (bod);
      IDENTITY (R);
    }
    break;
  }

  /* save formulation */
  bod->form = form;

  /* save rought mesh if needed */
  if (msh != shp->data) bod->msh = msh;
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
  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
      TL_dynamic_init (bod);
      break;
    case BODY_COROTATIONAL:
      BC_dynamic_init (bod);
      break;
  }

  unit_body_force (bod);
}

/* estimate critical step for the dynamic scheme */
double FEM_Dynamic_Critical_Step (BODY *bod)
{
  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
      return TL_dynamic_critical_step (bod);
    case BODY_COROTATIONAL:
      return BC_dynamic_critical_step (bod);
  }

  return DBL_MAX;
}

/* perform the initial half-step of the dynamic scheme */
void FEM_Dynamic_Step_Begin (BODY *bod, double time, double step)
{
  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
      TL_dynamic_step_begin (bod, time, step);
      break;
    case BODY_COROTATIONAL:
      BC_dynamic_step_begin (bod, time, step);
      break;
  }

}

/* perform the final half-step of the dynamic scheme */
void FEM_Dynamic_Step_End (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  double half = 0.5 * step,
	*energy = bod->energy,
	*fext = FEM_FEXT (bod),
	*fint = FEM_FINT (bod),
	*u0 = FEM_VEL0 (bod),
	*u = bod->velo,
	*ue = u + n,
	*iu0 = u0,
	*iu = u,
	*dq, *idq;

  ERRMEM (idq = dq = malloc (sizeof (double [n])));

  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
      TL_dynamic_step_end (bod, time, step);
      break;
    case BODY_COROTATIONAL:
      BC_dynamic_step_end (bod, time, step);
      break;
  }

  for (; iu < ue; idq ++, iu ++, iu0 ++) *idq = half * ((*iu) + (*iu0)); /* dq = (h/2) * {u(t) + u(t+h)} */
  energy [EXTERNAL] += blas_ddot (n, dq, 1, fext, 1);
  energy [INTERNAL] += blas_ddot (n, dq, 1, fint, 1);

  if (bod->msh) /* in such case SHAPE_Update will not update "rough" mesh */
  {
    MESH *msh = bod->msh;
    double (*cur) [3] = msh->cur_nodes,
	   (*ref) [3] = msh->ref_nodes,
	   (*end) [3] = ref + msh->nodes_count,
	    *q = bod->conf;

    for (; ref < end; cur ++, ref ++, q += 3) { ADD (ref[0], q, cur[0]); }
  }

  free (dq);
}

/* initialise static time stepping */
void FEM_Static_Init (BODY *bod)
{
  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
      TL_static_init (bod);
      break;
    case BODY_COROTATIONAL:
      BC_static_init (bod);
      break;
  }

  unit_body_force (bod);
}

/* perform the initial half-step of the static scheme */
void FEM_Static_Step_Begin (BODY *bod, double time, double step)
{
  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
      TL_static_step_begin (bod, time, step);
      break;
    case BODY_COROTATIONAL:
      BC_static_step_begin (bod, time, step);
      break;
  }
}

/* perform the final half-step of the static scheme */
void FEM_Static_Step_End (BODY *bod, double time, double step)
{
  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
      TL_static_step_end (bod, time, step);
      break;
    case BODY_COROTATIONAL:
      BC_static_step_end (bod, time, step);
      break;
  }

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
  double shapes [MAX_NODES], q [MAX_NODES][3], point [3];
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  int n, i;

  if (bod->msh)
  {
    msh = bod->msh;
    ASSERT_DEBUG_EXT ((cvx = gobj), "NULL geometric object for the separate mesh FEM scenario");
    ele = stabbed_referential_element (msh, cvx->ele, cvx->nele, X);
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
    n = element_shapes (ele->type, point, shapes);
    element_displacements (bod->conf, ele, q);

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
  double shapes [MAX_NODES], q [MAX_NODES][3];
  int i, n;

  n = element_shapes (ele->type, point, shapes);
  element_displacements (bod->conf, ele, q);

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
    ele = stabbed_spatial_element (msh, cvx->ele, cvx->nele, x);
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
  double shapes [MAX_NODES], u [MAX_NODES][3], point [3], vglo [3];
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  int i, n;

  if (bod->msh)
  {
    msh = bod->msh;
    ASSERT_DEBUG_EXT ((cvx = gobj), "NULL geometric object for the separate mesh FEM scenario");
    ele = stabbed_referential_element (msh, cvx->ele, cvx->nele, X);
    ASSERT_DEBUG (ele, "Invalid referential point stabbing an element");
  }
  else
  {
    msh = shp->data;
    ele = gobj;
  }

  referential_to_local (msh, ele, X, point);
  n = element_shapes (ele->type, point, shapes);

  if (prevel)
  {
    element_velocities (bod, ele, FEM_VEL0 (bod), u);
    SET (vglo, 0); for (i = 0; i < n; i ++) { ADDMUL (vglo, shapes [i], u [i], vglo); } /* vglo = N u0 */
    TVMUL (base, vglo, prevel); /* prevel = base' vglo */
  }

  if (curvel)
  {
    element_velocities (bod, ele, bod->velo, u);
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
    ele = stabbed_referential_element (msh, cvx->ele, cvx->nele, X);
    ASSERT_DEBUG (ele, "Invalid referential point stabbing an element");
  }
  else
  {
    msh = shp->data;
    ele = gobj;
  }

  int i [] = {0, 1, 2, 0, 1, 2, 0, 1, 2}, p [] = {0, 3, 6, 9};
  MX_CSC (base_trans, 9, 3, 3, p, i);
  double point [3];
  MX *N, *H;

  TNCOPY (base, base_trans.x);
  referential_to_local (msh, ele, X, point);
  N = element_shapes_matrix (bod, msh, ele, point);
  H = MX_Matmat (1.0, &base_trans, N, 0.0, NULL);
  MX_Destroy (N);

  if (bod->form == BODY_COROTATIONAL && bod->scheme != SCH_DEF_EXP
      && bod->dom->solver != BODY_SPACE_SOLVER) /* BODY_SPACE_SOLVER must see the regular H = E' N , rather than H R */
  {
    double *x = H->x, *y = x + H->nzmax,
           *R = FEM_ROT (bod), T [9];

    ASSERT_DEBUG ((y-x) % 9 == 0, "Number of nonzeros in H not divisble by 9");

    for (; x < y; x += 9)
    {
      NNMUL (x, R, T);
      NNCOPY (T, x); /* H = E' N R <=> rotaions gets shifted to H */
    }
  }

  return H;
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

/* get some values at a local point of an element */
void FEM_Element_Point_Values (BODY *bod, ELEMENT *ele, double *point, VALUE_KIND kind, double *values)
{
  MESH *msh = FEM_MESH (bod);

  switch (kind)
  {
  case VALUE_DISPLACEMENT:
  {
    MX *N = element_shapes_matrix (bod, msh, ele, point);
    MX_Matvec (1.0, N, bod->conf, 0.0, values);
    MX_Destroy (N);
  }
  break;
  case VALUE_VELOCITY:
  {
    MX *N = element_shapes_matrix (bod, msh, ele, point);
    MX_Matvec (1.0, N, bod->velo, 0.0, values);
    MX_Destroy (N);
  }
  break;
  case VALUE_STRESS:
  {
    cauchy_stress (bod, msh, ele, point, values);
  }
  break;
  case VALUE_MISES:
  {
    double stress [6];

    cauchy_stress (bod, msh, ele, point, stress);
    MISES (stress, values [0]);
  }
  break;
  case VALUE_STRESS_AND_MISES:
  {
    cauchy_stress (bod, msh, ele, point, values);
    MISES (values, values [6]);
  }
  break;
  }
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

    FEM_Element_Point_Values (bod, ele, point, kind, values);
  }
}

/* get some values at a curent mesh node (node points inside MESH->cur_nodes) */
void FEM_Cur_Node_Values (BODY *bod, double *node, VALUE_KIND kind, double *values)
{
  int i, j;
  double point [3];
  MESH *msh = FEM_MESH (bod);
  int n = (node_t) node - msh->cur_nodes;
  ELEMENT *ele = MESH_Element_With_Node (msh, n, point);
  ASSERT (ele, ERR_MSH_ELEMENT_WITH_NODE);

  if (kind >= VALUE_STRESS) /* average from neigbouring elements */
  {
    double X [3], v [7];
    SET *set, *item;
    int m;

    local_to_referential (msh, ele, point, X);

    set = NULL;

    MESH_Elements_Around_Node (ele, n, &set);

    switch ((int) kind)
    {
    case VALUE_STRESS: m = 6; break;
    case VALUE_MISES: m = 1; break;
    case VALUE_STRESS_AND_MISES: m = 7; break;
    }

    for (j = 0; j < m; j ++) values [j] = 0;

    for (item = SET_First (set); item; item = SET_Next (item))
    {
      ele = item->data;

      referential_to_local (msh, ele, X, point);

      FEM_Element_Point_Values (bod, ele, point, kind, v);

      for (j = 0; j < m; j ++) values [j] += v [j];
    }

    i = SET_Size (set);

    for (j = 0; j < m; j ++) values [j] /= (double)i;

    SET_Free (NULL, &set);
  }
  else FEM_Element_Point_Values (bod, ele, point, kind, values);
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
  return bod->dofs + (bod->form == BODY_COROTATIONAL ? 9 : 0);
}

/* get velocity packing size */
int FEM_Velo_Pack_Size (BODY *bod)
{
  return 4 * bod->dofs; /* velo, vel0, fext, fint */
}
#endif

/* compute c = alpha * INVERSE (bod) * b + beta * c */
void FEM_Invvec (double alpha, BODY *bod, double *b, double beta, double *c)
{
  switch (bod->form)
  {
  case TOTAL_LAGRANGIAN:
    MX_Matvec (alpha, bod->inverse, b, beta, c);
    break;
  case BODY_COROTATIONAL:
    BC_matvec (alpha, bod->inverse, FEM_ROT (bod), b, beta, c);
    break;
  }
}
