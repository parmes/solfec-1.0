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
#include "sol.h"
#include "fem.h"
#include "alg.h"
#include "bla.h"
#include "hyb.h"
#include "cvi.h"
#include "gjk.h"
#include "kdt.h"
#include "svk.h"
#include "but.h"
#include "err.h"
#if OMP
#include <omp.h>
#include "ompu.h"
#endif

typedef double (*node_t) [3]; /* mesh node */

#define IMP_EPS 1E-9
#define MAX_ITERS 64
#define MAX_NODES 20
#define DOM_TOL 0.1
#define CUT_TOL 0.001
#define MESH_DOFS(msh) ((msh)->nodes_count * 3)
#define FEM_MESH_CONF(bod) ((bod)->form < BODY_COROTATIONAL_MODAL ? (bod)->conf : (bod)->conf + (bod)->dofs + 9) /* mesh space configuration */
#define FEM_MESH_VELO(bod) ((bod)->form < BODY_COROTATIONAL_MODAL ? (bod)->velo : (bod)->velo + (bod)->dofs * 4) /* mesh space velocity */
#define FEM_MESH_VEL0(bod) (FEM_MESH_VELO(bod) + MESH_DOFS(FEM_MESH(bod))) /* mesh space previous velocity */
#define FEM_MESH_MASS(bod) ((bod)->form < BODY_COROTATIONAL_MODAL ? (bod)->M->x : ((bod)->conf + (bod)->dofs + 9 + MESH_DOFS(FEM_MESH(bod)))) /* mesh space velocity */
#define FEM_VEL0(bod) ((bod)->velo + (bod)->dofs)     /* previous velocity */
#define FEM_FEXT(bod) ((bod)->velo + (bod)->dofs * 2) /* external force */
#define FEM_FINT(bod) ((bod)->velo + (bod)->dofs * 3) /* internal force */
#define FEM_FBOD(bod) ((bod)->form < BODY_COROTATIONAL_MODAL ? (bod)->velo + (bod)->dofs * 4 : FEM_MESH_VEL0(bod) + MESH_DOFS(FEM_MESH(bod))) /* unit body force */
#define FEM_ROT(bod) ((bod)->conf + (bod)->dofs) /* rotation */

/* ==================== INTEGRATION ======================= */

/* order 1 */
static const double I_TET1_X[] = {0.25};
static const double I_TET1_Y[] = {0.25};
static const double I_TET1_Z[] = {0.25};
static const double I_TET1_W[] = {0.16666666666666666};
#define             I_TET1_N      1

static const double I_PYR1_X[] = {0.00};
static const double I_PYR1_Y[] = {0.00};
static const double I_PYR1_Z[] = {0.25};
static const double I_PYR1_W[] = {1.333333333333333333};
#define             I_PYR1_N      1

static const double I_WED1_X[] = {0.25};
static const double I_WED1_Y[] = {0.25};
static const double I_WED1_Z[] = {0.00};
static const double I_WED1_W[] = {0.333333333333333333};
#define             I_WED1_N      1

static const double I_HEX1_X[] = {0.0};
static const double I_HEX1_Y[] = {0.0};
static const double I_HEX1_Z[] = {0.0};
static const double I_HEX1_W[] = {8.0};
#define             I_HEX1_N      1

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

#define ISQR3_PYR1 0.263184055569714
#define ISQR3_PYR2 0.506616303349788
#define _1SUB_PYR1 0.544151844011225
#define _1SUB_PYR2 0.122514822655441
#define _1_WPYR1 0.100785882079825
#define _1_WPYR2 0.232547451253508
static const double I_PYR2_X[] = {-ISQR3_PYR1,  ISQR3_PYR1, ISQR3_PYR1, -ISQR3_PYR1, -ISQR3_PYR2,  ISQR3_PYR2, ISQR3_PYR2, -ISQR3_PYR2};
static const double I_PYR2_Y[] = {-ISQR3_PYR1, -ISQR3_PYR1, ISQR3_PYR1,  ISQR3_PYR1, -ISQR3_PYR2, -ISQR3_PYR2, ISQR3_PYR2,  ISQR3_PYR2};
static const double I_PYR2_Z[] = {_1SUB_PYR1,  _1SUB_PYR1,  _1SUB_PYR1,  _1SUB_PYR1,  _1SUB_PYR2,  _1SUB_PYR2, _1SUB_PYR2,  _1SUB_PYR2};
static const double I_PYR2_W[] = {_1_WPYR1, _1_WPYR1, _1_WPYR1, _1_WPYR1, _1_WPYR2, _1_WPYR2, _1_WPYR2,  _1_WPYR2};
#define             I_PYR2_N      8

#define ISQR3 0.57735026918962584 
static const double I_WED2_X[] = {0.66666666666666667, 0.16666666666666667, 0.16666666666666667, 0.66666666666666667, 0.16666666666666667, 0.16666666666666667};
static const double I_WED2_Y[] = {0.16666666666666667, 0.66666666666666667, 0.16666666666666667, 0.16666666666666667, 0.66666666666666667, 0.16666666666666667};
static const double I_WED2_Z[] = {-ISQR3, -ISQR3, -ISQR3, ISQR3, ISQR3, ISQR3};
static const double I_WED2_W[] = {0.16666666666666667, 0.16666666666666667, 0.16666666666666667, 0.16666666666666667, 0.16666666666666667, 0.16666666666666667};
#define             I_WED2_N      6

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
static const double I_HEX3_X [] = {0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834};
static const double I_HEX3_Y [] = {-0.7745966692414834, -0.7745966692414834, -0.7745966692414834, 0.0, 0.0, 0.0, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834, 0.0, 0.0, 0.0, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834, -0.7745966692414834, -0.7745966692414834, -0.7745966692414834, 0.0, 0.0, 0.0, 0.7745966692414834, 0.7745966692414834, 0.7745966692414834};
static const double I_HEX3_Z [] = {-0.7745966692414834, 0.0, 0.7745966692414834, -0.7745966692414834, 0.0, 0.7745966692414834, -0.7745966692414834, 0.0, 0.7745966692414834, -0.7745966692414834, 0.0, 0.7745966692414834, -0.7745966692414834, 0.0, 0.7745966692414834, -0.7745966692414834, 0.0, 0.7745966692414834, -0.7745966692414834, 0.0, 0.7745966692414834, -0.7745966692414834, 0.0, 0.7745966692414834, -0.7745966692414834, 0.0, 0.7745966692414834};
static const double I_HEX3_W [] = {0.1714677640603567, 0.27434842249657065, 0.1714677640603567, 0.27434842249657065, 0.43895747599451296, 0.27434842249657065, 0.1714677640603567, 0.27434842249657065, 0.1714677640603567, 0.27434842249657065, 0.43895747599451296, 0.27434842249657065, 0.438957475994513, 0.7023319615912208, 0.438957475994513, 0.27434842249657065, 0.43895747599451296, 0.27434842249657065, 0.1714677640603567, 0.27434842249657065, 0.1714677640603567, 0.27434842249657065, 0.43895747599451296, 0.27434842249657065, 0.1714677640603567, 0.27434842249657065, 0.1714677640603567};
#define             I_HEX3_N    27 

/* order 4 */

static const double I_HEX4_X [] = {-0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526};
static const double I_HEX4_Y [] = {-0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.8611363115940526, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, -0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.3399810435848563, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526, 0.8611363115940526};
static const double I_HEX4_Z [] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526, -0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
static const double I_HEX4_W [] = {0.04209147749053145, 0.07891151579507055, 0.07891151579507055, 0.04209147749053145, 0.07891151579507055, 0.1479403360567813, 0.1479403360567813, 0.07891151579507055, 0.07891151579507055, 0.1479403360567813, 0.1479403360567813, 0.07891151579507055, 0.04209147749053145, 0.07891151579507055, 0.07891151579507055, 0.04209147749053145, 0.07891151579507055, 0.1479403360567813, 0.1479403360567813, 0.07891151579507055, 0.1479403360567813, 0.27735296695391304, 0.27735296695391304, 0.1479403360567813, 0.1479403360567813, 0.27735296695391304, 0.27735296695391304, 0.1479403360567813, 0.07891151579507055, 0.1479403360567813, 0.1479403360567813, 0.07891151579507055, 0.07891151579507055, 0.1479403360567813, 0.1479403360567813, 0.07891151579507055, 0.1479403360567813, 0.27735296695391304, 0.27735296695391304, 0.1479403360567813, 0.1479403360567813, 0.27735296695391304, 0.27735296695391304, 0.1479403360567813, 0.07891151579507055, 0.1479403360567813, 0.1479403360567813, 0.07891151579507055, 0.04209147749053145, 0.07891151579507055, 0.07891151579507055, 0.04209147749053145, 0.07891151579507055, 0.1479403360567813, 0.1479403360567813, 0.07891151579507055, 0.07891151579507055, 0.1479403360567813, 0.1479403360567813, 0.07891151579507055, 0.04209147749053145, 0.07891151579507055, 0.07891151579507055, 0.04209147749053145};
#define             I_HEX4_N    64

/* collective rules */
static const double *I_TET_X [] = {NULL, I_TET1_X, I_TET2_X};
static const double *I_TET_Y [] = {NULL, I_TET1_Y, I_TET2_Y};
static const double *I_TET_Z [] = {NULL, I_TET1_Z, I_TET2_Z};
static const double *I_TET_W [] = {NULL, I_TET1_W, I_TET2_W};
static const int     I_TET_N [] = {   0, I_TET1_N, I_TET2_N};

static const double *I_PYR_X [] = {NULL, I_PYR1_X, I_PYR2_X};
static const double *I_PYR_Y [] = {NULL, I_PYR1_Y, I_PYR2_Y};
static const double *I_PYR_Z [] = {NULL, I_PYR1_Z, I_PYR2_Z};
static const double *I_PYR_W [] = {NULL, I_PYR1_W, I_PYR2_W};
static const int     I_PYR_N [] = {   0, I_PYR1_N, I_PYR2_N};

static const double *I_WED_X [] = {NULL, I_WED1_X, I_WED2_X};
static const double *I_WED_Y [] = {NULL, I_WED1_Y, I_WED2_Y};
static const double *I_WED_Z [] = {NULL, I_WED1_Z, I_WED2_Z};
static const double *I_WED_W [] = {NULL, I_WED1_W, I_WED2_W};
static const int     I_WED_N [] = {   0, I_WED1_N, I_WED2_N};

static const double *I_HEX_X [] = {NULL, I_HEX1_X, I_HEX2_X, I_HEX3_X, I_HEX4_X};
static const double *I_HEX_Y [] = {NULL, I_HEX1_Y, I_HEX2_Y, I_HEX3_Y, I_HEX4_Y};
static const double *I_HEX_Z [] = {NULL, I_HEX1_Z, I_HEX2_Z, I_HEX3_Z, I_HEX4_Z};
static const double *I_HEX_W [] = {NULL, I_HEX1_W, I_HEX2_W, I_HEX3_W, I_HEX4_W};
static const int     I_HEX_N [] = {   0, I_HEX1_N, I_HEX2_N, I_HEX3_N, I_HEX4_N};

static const double *I_TRI_X [] = {NULL, I_TRI1_X, I_TRI2_X};
static const double *I_TRI_Y [] = {NULL, I_TRI1_Y, I_TRI2_Y};
static const double *I_TRI_W [] = {NULL, I_TRI1_W, I_TRI2_W};
static const int     I_TRI_N [] = {   0, I_TRI1_N, I_TRI2_N};

static const double *I_QUA_X [] = {NULL, I_QUA1_X, I_QUA2_X};
static const double *I_QUA_Y [] = {NULL, I_QUA1_Y, I_QUA2_Y};
static const double *I_QUA_W [] = {NULL, I_QUA1_W, I_QUA2_W};
static const int     I_QUA_N [] = {   0, I_QUA1_N, I_QUA2_N};

#define MAX_ORDER 4 /* integration order bound */

/* minimal local coords */
static double mincoord [9][3] =
{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0},
 {0, 0, 0}, /* tet */
 {-1, -1, 0}, /* pyr */
 {0, 0, -1}, /* wed */
 {0, 0, 0},
 {-1, -1, -1}}; /* hex */

/* maximal local coords */
static double maxcoord [9][3] =
{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0},
 {1, 1, 1}, /* tet */
 {1, 1, 1}, /* pyr */
 {1, 1, 1}, /* wed */
 {0, 0, 0},
 {1, 1, 1}}; /* hex */

/* load 3D integrator data */
inline static int integrator3d_load (int type, int order, const double **X, const double **Y, const double **Z, const double **W)
{
  int N;

  ASSERT_DEBUG (order >= 1 && order <= MAX_ORDER, "Integration order out of bounds");

  switch (type)
  {
  case 4:
  {
    *X = I_TET_X [order];
    *Y = I_TET_Y [order];
    *Z = I_TET_Z [order];
    *W = I_TET_W [order];
     N = I_TET_N [order];
  }
  break;
  case 5:
  {
    *X = I_PYR_X [order];
    *Y = I_PYR_Y [order];
    *Z = I_PYR_Z [order];
    *W = I_PYR_W [order];
     N = I_PYR_N [order];
  }
  break;
  case 6:
  {
    *X = I_WED_X [order];
    *Y = I_WED_Y [order];
    *Z = I_WED_Z [order];
    *W = I_WED_W [order];
     N = I_WED_N [order];
  }
  break;
  case 8:
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
  if (type == 8 && entity == MASS) return 4;
  return 2; /* XXX */
}

/* load 2D integrator data */
inline static int integrator2d_load (int type, int order, const double **X, const double **Y, const double **W)
{
  int N;

  ASSERT_DEBUG (order >= 1 && order <= MAX_ORDER, "Integration order out of bounds");

  switch (type)
  {
  case 3:
  {
    *X = I_TRI_X [order];
    *Y = I_TRI_Y [order];
    *W = I_TRI_W [order];
     N = I_TRI_N [order];
  }
  break;
  case 4:
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
static int integrator2d_order (int type)
{
  if (type == 3) return 1;
  else return 2; /* quads */
}

/* avoid warnings related to unsed point[], weight variables */
#if defined(__clang__)
#elif defined(__GNUC__) /* GCC */
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

/* element integration; note that below __t__->ver, __t__->center are in local element coordinates already */
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
    __N__ = integrator3d_load (4, 1, &__X__, &__Y__, &__Z__, &__W__);\
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
	  COPY (__t__->ver [0], __subnodes__ [2]);\
	  COPY (__t__->ver [1], __subnodes__ [1]);\
	  COPY (__t__->ver [2], __subnodes__ [0]);\
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
#define INTEGRAL2D_BEGIN(TYPE)\
{\
  const double *__X__, *__Y__, *__W__;\
  double point [2], weight;\
  int __N__, __k__;\
\
  __N__ = integrator2d_load (TYPE, integrator2d_order (TYPE), &__X__, &__Y__, &__W__);\
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
inline static double face_det (FACE *fac, node_t nodes, double *point, double *normal)
{
  double derivs [16], d0 [3], d1[3], tmp [3], *d;
  int i, n;

  n = face_derivs (fac, point, derivs);

  if (!normal) normal = tmp;

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
inline static int element_nodes (node_t heap, int type, int *nodes, node_t stack);
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

/* create face integration data;
 * FACE->idata = {shapes0, normal0, J0*wgt0, shapes1, normal1, J1*wgt1, shapes2, normal2, J2*wgt2, shapes3, normal3, J3*wgt3} */
#define FACE_SHAPES(fac, i) ((fac)->idata+(i)*((fac)->type+4))
#define FACE_NORMAL(fac, i) (FACE_SHAPES(fac,i)+(fac)->type)
static void create_face_integration_data (MESH *msh)
{
  double refn [4][3], J, *N, *shapes;
  FACE *fac;
  int n;

  for (fac = msh->faces; fac; fac = fac->n)
  {
    n = 0; INTEGRAL2D_BEGIN (fac->type) { n ++; } INTEGRAL2D_END ()

    ERRMEM (fac->idata = malloc (sizeof (double [n * (fac->type + 4)])));

    face_nodes (msh->ref_nodes, fac->type, fac->nodes, refn);

    INTEGRAL2D_BEGIN (fac->type) /* defines point and weight */
    {
      shapes = FACE_SHAPES (fac, __k__);
      N = FACE_NORMAL (fac, __k__);
      face_shapes (fac, point, shapes);
      J = face_det (fac, refn, point, N);
      N [3] = J * weight;
      NORMALIZE (N);
    }
    INTEGRAL2D_END ()
  }
}

/* =============================== ELEMENT ================================ */

/* linear tetrahedron shape functions */
inline static void tet_o1_shapes (double *point, double *shapes)
{
  shapes [0] = 1.0 - (point [0] + point [1] + point [2]);
  shapes [1] = point [0];
  shapes [2] = point [1];
  shapes [3] = point [2];
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
  double ratio;

  if (point [2] != 1.0) ratio = (point[0]*point[1]*point[2]) / (1 - point[2]);
  else ratio = 0.0; /* XXX */

  shapes [0] = 0.25 * ((1.0 + point[0])*(1.0 + point[1]) - point[2] + ratio); /* rabbit functions */
  shapes [1] = 0.25 * ((1.0 - point[0])*(1.0 + point[1]) - point[2] - ratio);
  shapes [2] = 0.25 * ((1.0 - point[0])*(1.0 - point[1]) - point[2] + ratio);
  shapes [3] = 0.25 * ((1.0 + point[0])*(1.0 - point[1]) - point[2] - ratio);
  shapes [4] = point [2];
}

/* linear wedge shape functions */
inline static void wed_o1_shapes (double *point, double *shapes)
{
  shapes [0] = 0.5 * (1.0 - point[0] - point[1])*(1.0 - point[2]);
  shapes [1] = 0.5 * point[0]*(1.0 - point[2]);
  shapes [2] = 0.5 * point[1]*(1.0 - point[2]);
  shapes [3] = 0.5 * (1.0 - point[0] - point[1])*(1.0 + point[2]);
  shapes [4] = 0.5 * point[0]*(1.0 + point[2]);
  shapes [5] = 0.5 * point[1]*(1.0 + point[2]);
}

/* linear tetrahedron shape derivatives */
inline static void tet_o1_derivs (double *point, double *derivs)
{
  derivs [0] = -1.0;
  derivs [1] = -1.0;
  derivs [2] = -1.0;

  derivs [3] = 1.0;
  derivs [4] = 0.0;
  derivs [5] = 0.0;

  derivs [6] = 0.0;
  derivs [7] = 1.0;
  derivs [8] = 0.0;

  derivs [9]  = 0.0;
  derivs [10] = 0.0;
  derivs [11] = 1.0;
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
  double drat0, drat1, drat2;

  if (point [2] != 1.0)
  {
    drat0 = (point[1]*point[2]) / (1.0 - point[2]),
    drat1 = (point[0]*point[2]) / (1.0 - point[2]),
    drat2 = (point[0]*point[1]) / (1.0 - point[2])*(1.0 - point[2]);
  }
  else drat0 = drat1 = drat2 = 0.0; /* XXX */

  derivs [0] = 0.25 * ((1.0 + point[1]) + drat0);
  derivs [1] = 0.25 * ((1.0 + point[0]) + drat1);
  derivs [2] = 0.25 * (-1.0 + drat2);

  derivs [3] = 0.25 * (-(1.0 + point[1]) - drat0);
  derivs [4] = 0.25 * ((1.0 - point[0]) - drat1);
  derivs [5] = 0.25 * (-1.0 - drat2);

  derivs [6] = 0.25 * (-(1.0 - point[1]) + drat0);
  derivs [7] = 0.25 * (-(1.0 - point[0]) + drat1);
  derivs [8] = 0.25 * (-1.0 + drat2);

  derivs [ 9] = 0.25 * ((1.0 - point[1]) - drat0);
  derivs [10] = 0.25 * (-(1.0 + point[0]) - drat1);
  derivs [11] = 0.25 * (-1.0 - drat2);

  derivs [12] = 0.0;
  derivs [13] = 0.0;
  derivs [14] = 1.0;
}

/* linear wedge shape derivatives  */
inline static void wed_o1_derivs (double *point, double *derivs)
{
  derivs [0] = -0.5 * (1.0 - point[2]);
  derivs [1] = -0.5 * (1.0 - point[2]);
  derivs [2] = -0.5 * (1.0 - point[0] - point[1]);

  derivs [3] =  0.5 * (1.0 - point[2]);
  derivs [4] =  0.0;
  derivs [5] = -0.5 * point[0];

  derivs [6] =  0.0;
  derivs [7] =  0.5 * (1.0 - point[2]);
  derivs [8] = -0.5 * point[1];

  derivs [ 9] = -0.5 * (1.0 + point[2]);
  derivs [10] = -0.5 * (1.0 + point[2]);
  derivs [11] =  0.5 * (1.0 - point[0] - point[1]);

  derivs [12] = 0.5 * (1.0 + point[2]);
  derivs [13] = 0.0;
  derivs [14] = 0.5 * point[0];

  derivs [15] = 0.0;
  derivs [16] = 0.5 * (1.0 + point[2]);
  derivs [17] = 0.5 * point[1];
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

/* copy element nodal values into a local table */
static int element_nodal_values (double *heap, ELEMENT *ele, double (*q) [3])
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
  double A [9], B [3], error;
  int i, j, k, l, n, ipiv [3];

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

    ASSERT (lapack_dgesv (3, 1, A, 3, ipiv, B, 3) == 0, ERR_FEM_COORDS_INVERT);
    ADD (local, B, local);
    error = LEN (B);

  } while (++l < 64 && error > 1E-9);

  ASSERT (l < 64, ERR_FEM_COORDS_INVERT);

  for (i = 0; i < 3; i ++) /* project back onto the element if the point is off due to roundoff */
  {
    if (local [i] < mincoord [ele->type][i]) local [i] = mincoord [ele->type][i];
    else if (local [i] > maxcoord [ele->type][i]) local [i] = maxcoord [ele->type][i];
  }
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
  int dofs = MESH_DOFS (msh);
  MX *N;

  o = element_shapes (ele->type, point, shapes);

  n = dofs;

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

/* return number of integration points for internal force computation */
static int number_of_integration_points (ELEMENT *ele)
{
  int i = 0;

  INTEGRATE3D (ele->type, INTF, ele->dom, ele->domnum, i ++;)

  return i;
}

/* allocate bulk material states at integration points */
static void allocate_element_states (MESH *msh, BULK_MATERIAL *mat)
{
  ELEMENT *ele;
  int nip;

  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    if (mat->nstate && ele->state == NULL)
    {
      nip = number_of_integration_points (ele);
      ERRMEM (ele->state = MEM_CALLOC (nip * sizeof (double [mat->nstate])));
    }
  }
}

/* simplex integrated volume ==
 * shape functions integrated volume test */
static void test_volume_integral (MESH *msh, double ref_volume, int body_id)
{
  double J, volume, nodes [MAX_NODES][3];
  ELEMENT *ele;
  int bulk;

  for (ele = msh->surfeles, bulk = 0, volume = 0.0; ele; )
  {
    element_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

    INTEGRATE3D (ele->type, MASS, ele->dom, ele->domnum,

      J = element_det (ele->type, nodes, point, NULL);
      volume += J * weight;
    )

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }

  WARNING (fabs (volume - ref_volume) < CUT_TOL * ref_volume,
    "FEM BODY %d:\nShape volume is %g.\nIntegrated volume is %g.\n"
    "Error %g is beyond the tolerance of %g.\n"
    "This issue occurs when your background hexahedrons are not rectilinear.\n"
    "Refine the background mesh or make it rectilinear. Alternately, use a tetrahedral background mesh.\n",
    body_id, ref_volume, volume, fabs (ref_volume - volume) / ref_volume, CUT_TOL);
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
  double nodes [MAX_NODES][3], q [MAX_NODES][3], derivs [3*MAX_NODES], field [MAX_NFIELD],
	 shapes [MAX_NODES], F0 [9], F [9], P [9], K [81], KB [9], J, integral, *B, *p;
  BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
  double *bfld = bod->field, *conf = FEM_MESH_CONF (bod);
  int i, j, n, m, ip = 0,
      nbfld = mat->nfield,
      *nod = ele->nodes;

  ASSERT_TEXT (mat->nfield < MAX_NFIELD, "The maximum of %d field variables has been exceeded.\n", MAX_NFIELD);

  n = element_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  m = 3 * n;

  if (bod->form == BODY_COROTATIONAL ||
      bod->form >= BODY_COROTATIONAL_MODAL)
  {
    for (i = 0; i < n; i ++)
    {
      SET (q[i], 0); /* initial displacement */
    }
  }
  else
  {
    for (i = 0; i < n; i ++)
    {
      p = &conf [3 * nod [i]];
      COPY (p, q[i]); /* current displacement */
    }
  }

  for (i = 0, j = m * (derivative ? m : 1); i < j; i ++) g [i] = 0.0;

  INTEGRATE3D (ele->type, INTF, ele->dom, ele->domnum,

    J = element_det (ele->type, nodes, point, F0);
    element_gradient (ele->type, q, point, F0, derivs, F);
    integral = J * weight;

    if (bfld)
    {
      element_shapes (ele->type, point, shapes);

      for (i = 0; i < nbfld; i ++)
      {
        for (field [i] = 0.0, j = 0; j < n; j ++)
	{
	  field [i] += shapes [j] * bfld [nbfld * nod [j] + i]; /* interpolate fields from nodes */
	}
      }
    }

    if (derivative)
    {
      BULK_MATERIAL_ROUTINE (mat, ele->state + ip * mat->nstate, field, F, integral, NULL, K);

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
      BULK_MATERIAL_ROUTINE (mat, ele->state + ip * mat->nstate, field, F, integral, P, NULL);

      for (i = 0, B = derivs, p = g; i < n; i ++, B += 3, p += 3) { NVADDMUL (p, P, B, p); }
    }

    ip ++;
  )
}

/* compute elastic energy of individual element (and its volume if pvol != NULL) */
double FEM_Element_Internal_Energy (BODY *bod, MESH *msh, ELEMENT *ele, double *pvol)
{
  double nodes [MAX_NODES][3], q [MAX_NODES][3], derivs [3*MAX_NODES], F0 [9], F [9], J, integral;
  BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
  double *conf = FEM_MESH_CONF (bod), *p;
  int i, n, *nod = ele->nodes;

  n = element_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  for (i = 0; i < n; i ++)
  {
    p = &conf [3 * nod [i]];
    COPY (p, q[i]); /* current mesh space displacement */
  }

  integral = 0.0;

  if (pvol) *pvol = 0.0;

  INTEGRATE3D (ele->type, INTF, ele->dom, ele->domnum,

    J = element_det (ele->type, nodes, point, F0);
    element_gradient (ele->type, q, point, F0, derivs, F);
    integral += J * weight * SVK_Energy_C (lambda (mat->young, mat->poisson), mi (mat->young, mat->poisson), 1.0, F);

    if (pvol) *pvol += J * weight;
  )

  /* XXX/TODO => SVK material fixed above */

  return integral;
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

/* return closest integration point state */
static double* element_closest_ipoint_state (BODY *bod, ELEMENT *ele, double *p)
{
  if (ele->state == NULL) return NULL;

  BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
  double d [3], dot, dmin = DBL_MAX;
  int i = 0, j = 0;

  INTEGRATE3D (ele->type, INTF, ele->dom, ele->domnum,

      SUB (p, point, d);
      dot = DOT (d, d);
      if (dot < dmin)
      {
       dmin = dot;
       j = i;
      }

      i ++;
  )

  return &ele->state [j * mat->nstate];
}

/* =========================== GENERAL ================================== */

/* compute deformation gradient at a local point */
static void deformation_gradient (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *F)
{
  double derivs [3*MAX_NODES], nodes [MAX_NODES][3], q [MAX_NODES][3], F0 [9], *p;
  double *conf = FEM_MESH_CONF (bod);
  int i, n;

  n = element_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  for (i = 0; i < n; i ++)
  {
    p = &conf [3 * ele->nodes [i]];
    COPY (p, q[i]);
  }

  element_det (ele->type, nodes, point, F0);
  element_gradient (ele->type, q, point, F0, derivs, F);
}

/* compute Cauchy stress at a local point */
static void cauchy_stress (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *values)
{
  BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
  double P [9], F [9], J, shapes [MAX_NODES], field [MAX_NFIELD], *bfld = bod->field;
  double *state = element_closest_ipoint_state (bod, ele, point);
  int i, j, n, nbfld = mat->nfield, *nod = ele->nodes;

  if (bfld)
  {
    n = element_shapes (ele->type, point, shapes);

    for (i = 0; i < nbfld; i ++)
    {
      for (field [i] = 0.0, j = 0; j < n; j ++)
      {
	field [i] += shapes [j] * bfld [nbfld * nod [j] + i]; /* interpolate fields from nodes */
      }
    }
  }

  deformation_gradient (bod, msh, ele, point, F);

  J = BULK_MATERIAL_ROUTINE (mat, state, field, F, 1.0, P, NULL);

  values [0] = (F[0]*P[0]+F[3]*P[1]+F[6]*P[2])/J; /* sx  */
  values [1] = (F[1]*P[3]+F[4]*P[4]+F[7]*P[5])/J; /* sy  */
  values [2] = (F[2]*P[6]+F[5]*P[7]+F[8]*P[8])/J; /* sz  */
  values [3] = (F[0]*P[3]+F[3]*P[4]+F[6]*P[5])/J; /* sxy */
  values [4] = (F[0]*P[6]+F[3]*P[7]+F[6]*P[8])/J; /* sxz */
  values [5] = (F[2]*P[0]+F[5]*P[1]+F[8]*P[2])/J; /* syz */

  /* XXX: BODY_COROTATIONAL case is handled the same ways as TOTAL_LAGRANGIAN (simplification) */
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
    if (ELEMENT_Contains_Point (msh, *ele, x, 0)) return *ele;

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
    if (ELEMENT_Contains_Point (msh, *ele, X, 1)) return *ele;

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

/* compute constraints force r = SUM H' R, where r is in the modeling space (could be reduced) */
static void fem_constraints_force (BODY *bod, double *r)
{
  double *rmsh; /* reaction in the mesh space */
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  SET *node;
  int i, dofs;

  msh = FEM_MESH (bod);
  dofs = MESH_DOFS (msh);

  if (bod->form >= BODY_COROTATIONAL_MODAL)
  {
    ERRMEM (rmsh = malloc (dofs * sizeof (double)));
  }
  else rmsh = r;

  for (i = 0; i < dofs; i ++) rmsh [i] = 0.0;

  for (node = SET_First (bod->con); node; node = SET_Next (node))
  {
    CON *con = node->data;
    short isma = (bod == con->master);
    double *X = (isma ? con->mpnt : con->spnt);

    if (bod->msh)
    {
      cvx = (isma ? con->msgp->gobj : con->ssgp->gobj);
      ele = stabbed_referential_element (msh, cvx->ele, cvx->nele, X); /* TODO: optimize */
    }
    else ele = (isma ? con->msgp->gobj : con->ssgp->gobj);

    accumulate_reac (bod, msh, ele, X, con->base, con->R, isma, rmsh);

    if (isma && bod == con->slave) /* self-contact */
    {
      X = con->spnt;

      if (bod->msh)
      {
	cvx = con->ssgp->gobj;
	ele = stabbed_referential_element (msh, cvx->ele, cvx->nele, X); /* TODO: optimize */
      }
      else ele = con->ssgp->gobj;

      accumulate_reac (bod, msh, ele, X, con->base, con->R, 0, rmsh);
    }
  }

  if (bod->form >= BODY_COROTATIONAL_MODAL) /* H' = bod->evec' * R' * N' * E */
  {
    double *R = FEM_ROT (bod),
	   *x = rmsh,
	   *y = x + dofs,
	    z [3];

    for (;x < y; x += 3)
    {
      COPY (x, z);
      TVMUL (R, z, x);
    }

    MX_Matvec (1.0, MX_Tran (bod->evec), rmsh, 0.0, r); /* r = bod->evec' * R' * SUM { N' * E * R } */

    free (rmsh); /* clean */
  }
}

/* compute inernal force */
static void internal_force (BODY *bod, double *fint)
{
  MESH *msh = FEM_MESH (bod);
  int dofs = MESH_DOFS (msh);
  double g [24], *v, *w;
  int i;

  for (i = 0; i < dofs; i ++) fint [i] = 0.0;

#if OMP
  int j, n;
  ELEMENT **pele = ompu_elements (msh, &n);
  omp_lock_t *locks = ompu_locks (msh->nodes_count);
  #pragma omp parallel for shared (pele, fint, locks) private (g, i, v, w)
  for (j = 0; j < n; j ++)
  {
    element_internal_force (0, bod, msh, pele[j], g);

    for (i = 0, v = g; i < pele[j]->type; i ++, v += 3)
    {
      int k = pele[j]->nodes[i];
      w = &fint [k * 3];
      omp_set_lock(&locks[k]);
      ACC (v, w);
      omp_unset_lock (&locks[k]);
    }
  }
  ompu_locks_free (locks, msh->nodes_count);
  free (pele);
#else
  ELEMENT *ele;
  int bulk;
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
#endif
}

/* compute inernal energy */
static double internal_energy (BODY *bod)
{
  MESH *msh = FEM_MESH (bod);
  double energy = 0.0;

#if OMP
  int j, n;
  ELEMENT **pele = ompu_elements (msh, &n);
  #pragma omp parallel for shared (pele) reduction(+:energy)
  for (j = 0; j < n; j ++)
  {
    energy += FEM_Element_Internal_Energy (bod, msh, pele[j], NULL);
  }
  free (pele);
#else
  ELEMENT *ele;
  int bulk;
  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    energy += FEM_Element_Internal_Energy (bod, msh, ele, NULL);

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }
#endif

  return energy;
}

/* initialize unit body force */
static void unit_body_force (BODY *bod)
{
  double g [24], f [3], *v, *w;
  double *ubf = FEM_FBOD (bod);
  MESH *msh = FEM_MESH (bod);
  int dofs = MESH_DOFS (msh);
  ELEMENT *ele;
  int bulk, i;

  SET (f, 1.0);
  memset (ubf, 0, dofs * sizeof (double));
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

/* accumulate surface pressure */
static void accumulate_pressure (BODY *bod, MESH *msh, double *conf, int surfid, double value, double *fext)
{
  double q [4][3], nodes [4][3], *shapes, *N, p;
  int i, j, k, n;
  FACE *fac;

  for (fac = msh->faces; fac; fac = fac->n)
  {
    if (fac->surface == surfid)
    {
      face_nodes (msh->ref_nodes, fac->type, fac->nodes, nodes);
      face_displacements (conf, fac, q);
      for (j = 0; j < fac->type; j ++) ADD (nodes [j], q [j], nodes [j]);

      INTEGRAL2D_BEGIN (fac->type) /* defines point and weight */
      {
	shapes = FACE_SHAPES (fac, __k__);
	N = FACE_NORMAL (fac, __k__);
	p = value * N [3];

        for (i = 0; i < fac->type; i ++)
	{
	  k = fac->nodes [i];
	  for (n = 0; n < 3; n ++)
	    fext [3*k+n] += N [n] * p * shapes [i]; /* TODO/XXX: verify */
	}
      }
      INTEGRAL2D_END ()
    }
  }
}

/* compute external force */
static void external_force (BODY *bod, double time, double step, double *fext)
{
  double g [24], f [3], point [3], value, *v;
  double *ubf = FEM_FBOD (bod);
  MESH *msh = FEM_MESH (bod);
  int dofs = MESH_DOFS (msh);
  ELEMENT *ele;
  int i;

  /* zero forces */
  for (i = 0; i < dofs; i ++) fext [i] = 0.0;

  /* point forces */
  for (FORCE *frc = bod->forces; frc; frc = frc->next)
  {
    if (frc->func)
    {
      ERRMEM (v = MEM_CALLOC (bod->dofs * sizeof (double)));
      frc->func (frc->data, frc->call, FEM_Conf_Size (bod), bod->conf, bod->dofs, bod->velo, time, step, v);
      if (bod->form >= BODY_COROTATIONAL_MODAL) MX_Matvec (1.0, bod->evec, v, 1.0, fext); /* accumulate from reduced to global space */
      else blas_daxpy (bod->dofs, 1.0, v, 1, fext, 1);
      free (v);
    }
    else
    {
      value = TMS_Value (frc->data, time);

      if (frc->kind & PRESSURE)
      {
        accumulate_pressure (bod, bod->shape->data, FEM_MESH_CONF (bod), frc->surfid, value, fext);
      }
      else
      {
	COPY (frc->direction, f);
	SCALE (f, value);

	ele = MESH_Element_Containing_Point (msh, frc->ref_point, 1); /* TODO: optimize */

	ASSERT_TEXT (ele, "Point force seems to be applied outside of FEM mesh.\n"
			  "This should have been detected when interpreting input.\n"
			  "Please report this bug!\n");

	referential_to_local (msh, ele, frc->ref_point, point);

	if (frc->kind & CONVECTED)
	{ 
	  deformation_gradient (bod, msh, ele, point, g);
	  NVMUL (g, f, g+9);
	  COPY (g+9, f);
	}

	accumulate_point_force (bod, msh, ele, point, f, fext);
      }
    }
  }

  /* gravitation */
  if (bod->dom->gravity [0])
  {
    f [0] = TMS_Value (bod->dom->gravity [0], time);
    f [1] = TMS_Value (bod->dom->gravity [1], time);
    f [2] = TMS_Value (bod->dom->gravity [2], time);

    i = dofs / 3;
    blas_daxpy (i, f[0], ubf, 3, fext, 3);
    blas_daxpy (i, f[1], ubf+1, 3, fext+1, 3);
    blas_daxpy (i, f[2], ubf+2, 3, fext+2, 3);
  }
}
 
/* compute global tangent stiffness */
static MX* tangent_stiffness (BODY *bod, short spd)
{
  int i, j, k, l, n, dofs, *pp, *ii, *kk;
  double K [576], *A, *rowblk;
  MAP **col, *item;
  ELEMENT *ele;
  MESH *msh;
#if OMP
  MEM *blkmem,
      *mapmem;
#else
  MEM blkmem,
      mapmem;
#endif
  MX *tang;

  if (spd) spd = 1;
  msh = FEM_MESH (bod);
  dofs = MESH_DOFS (msh);
  ERRMEM (col = MEM_CALLOC (sizeof (MAP*) * dofs)); /* sparse columns */
#if OMP
  int threads;
  #pragma omp parallel
  #pragma omp master
  threads = omp_get_num_threads();
  ERRMEM (blkmem = MEM_CALLOC (sizeof(MEM) * threads));
  ERRMEM (mapmem = MEM_CALLOC (sizeof(MEM) * threads));
  for (i = 0; i < threads; i ++)
  {
    MEM_Init  (&blkmem[i], sizeof (double [3]), dofs);
    MEM_Init  (&mapmem[i], sizeof (MAP), dofs);
  }
#else
  MEM_Init  (&blkmem, sizeof (double [3]), dofs);
  MEM_Init  (&mapmem, sizeof (MAP), dofs);
#endif

#if OMP
  int ei, en;
  ELEMENT **pele = ompu_elements (msh, &en);
  omp_lock_t *locks = ompu_locks (msh->nodes_count);
  #pragma omp parallel for shared (bod, spd, msh, pele, locks, col, blkmem, mapmem) private (ele, K, k, A, l, j, n, i, rowblk)
  for (ei = 0; ei < en; ei ++)
#else
  short bulk;
  for (ele = msh->surfeles, bulk = 0; ele;
       ele = (ele->next ? ele->next : bulk ? NULL : msh->bulkeles),
       bulk = (ele == msh->bulkeles ? 1 : bulk)) /* for each element in mesh */
#endif
  {
#if OMP
    ele = pele[ei];
#endif
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
#if OMP
      omp_set_lock (&locks[ele->nodes[k]]);
#endif
      for (l = 0; l < 3; l ++) /* for each nodal degree of freedom */
      {
	j = 3 * ele->nodes [k] + l; /* for each global column index */

	for (n = 0; n < ele->type; n ++, A += 3) /* for each column row-block; shift column block pointer A */
	{
	  i = 3 * ele->nodes [n];

	  if (spd && i+2 < j) continue; /* skip upper triangle (leave diagonal overlaping blocks) */

	  if (!(rowblk = MAP_Find (col [j], (void*) (long) i, NULL))) /* if this row-block was not mapped */
	  {
#if OMP
	    ERRMEM (rowblk = MEM_Alloc (&blkmem[omp_get_thread_num()]));
	    MAP_Insert (&mapmem[omp_get_thread_num()], &col [j], (void*) (long) i, rowblk, NULL); /* map it */
#else
	    ERRMEM (rowblk = MEM_Alloc (&blkmem));
	    MAP_Insert (&mapmem, &col [j], (void*) (long) i, rowblk, NULL); /* map it */
#endif
	  }
	  ACC (A, rowblk); /* accumulate values */
	}
      }
#if OMP
      omp_unset_lock (&locks[ele->nodes[k]]);
#endif
    }
  }
#if OMP
  ompu_locks_free (locks, msh->nodes_count);
  free (pele);
#endif

  ERRMEM (pp = malloc (sizeof (int [dofs + 1]))); /* column pointers */

  for (pp [0] = 0, j = 0; j < dofs; j ++) pp [j+1] = pp [j] + 3 * MAP_Size (col [j]) - spd * (j % 3); /* subtract upper triangular j % 3 sticking out bits */

  ERRMEM (ii = malloc (sizeof (int [pp [dofs]]))); /* row indices */

  for (j = 0, kk = ii; j < dofs; j ++) /* initialize row index pointer; for each column */
  {
    for (item = MAP_First (col [j]); item; item = MAP_Next (item)) /* for each row-block */
    {
      i = (int) (long) item->key;
      if (spd && i < j) /* diagonal block with sticking out upper triangle */
      {
	for (n = 0; n < 3 - j % 3; n ++, kk ++) kk [0] = j + n;
      }
      else /* lower triangle block */
      {
	kk [0] = i;
	kk [1] = kk[0] + 1;
	kk [2] = kk[1] + 1;
	kk += 3;
      }
    }
  }

  tang = MX_Create (MXCSC, dofs, dofs, pp, ii); /* create tangent matrix structure */
  if (spd) tang->flags |= MXSPD;

  for (j = 0, A = tang->x; j < dofs; j ++) /* initialize column values pointer A; for each column */
  {
    for (item = MAP_First (col [j]); item; item = MAP_Next (item)) /* for each row-block */
    {
      i = (int) (long) item->key;
      rowblk = item->data;
      if (spd && i < j) /* diagonal block with sticking out upper triangle */
      {
	for (n = 0; n < 3 - j % 3; n ++, A ++) A [0] = rowblk [j % 3 + n];
      }
      else /* lower triangle block */
      {
        COPY (rowblk, A); /* copy values */
	A += 3;
      }
    }
  }

  free (ii);
  free (pp);
  free (col);
#if OMP
  for (i = 0; i < threads; i ++)
  {
    MEM_Release (&mapmem[i]);
    MEM_Release (&blkmem[i]);
  }
  free (mapmem);
  free (blkmem);
#else
  MEM_Release (&mapmem);
  MEM_Release (&blkmem);
#endif

  return tang;
}

/* compute diagonalized inertia operator */
static MX* diagonal_inertia (BODY *bod, short spd)
{
  MESH *msh = FEM_MESH (bod);
  int n = MESH_DOFS (msh),
      bulk,
     *p,
     *i,
      k;
  ELEMENT *ele;
  double *x;
  MX *M;

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

/* update nodal field variables */
static void update_fields (BODY *bod, double t)
{
  MESH *msh = FEM_MESH (bod);
  double (*cur) [3] = msh->cur_nodes,
	 (*end) [3] = cur + msh->nodes_count,
	 *bfld = bod->field;
  int i, nfld = bod->mat->nfield;
  FIELD **pfld = bod->mat->fld;

  if (bfld)
  {
    for (; cur < end; cur ++, bfld += nfld)
    {
      for (i = 0; i < nfld; i ++)
      {
	bfld [i] = FIELD_Value (pfld [i], cur[0][0], cur[0][1], cur[0][2], t);
      }
    }
  }
}

/* =================== TOAL LAGRANGIAN =================== */

/* compute out of balance force = fext - fint */
static void TL_dynamic_force (BODY *bod, double time, double step, double *fext, double *fint, double *force)
{

  external_force (bod, time, step, fext);

  internal_force (bod, fint);

  for (double *x = fext, *y = fint, *z = force, *u = z + bod->dofs; z < u; x ++, y ++, z ++) *z = (*x) - (*y);
}

/* the smame computation for the static case */
#define TL_static_force(bod, time, step, fext, fint, force) TL_dynamic_force (bod,time,step,fext,fint,force)

/* compute inverse operator for the implicit dynamic time stepping */
static void TL_dynamic_inverse (BODY *bod, double step, double *force)
{
  if (bod->inverse) MX_Destroy (bod->inverse);

  if (bod->K) MX_Destroy (bod->K);

  bod->K = tangent_stiffness (bod, 1);

  if (force)
  {
    /* account for the previous velocity */
    MX_Matvec (1.0 / step, bod->M, bod->velo, 1.0, force);

    /* account for the internal force increment */
    MX_Matvec (-0.25 * step, bod->K, bod->velo, 1.0, force);
  }

  /* calculate tangent operator A = M + (damping*h/2 + h*h/4) K */
  bod->inverse = MX_Add (1.0, bod->M, 0.5*bod->damping*step + 0.25*step*step, bod->K, NULL);

  /* invert A */
  MX_Inverse (bod->inverse, bod->inverse);
}

/* static time-stepping inverse */
static void TL_static_inverse (BODY *bod, double step)
{
  MX *M, *K;

  if (bod->M) M = bod->M; else bod->M = M = diagonal_inertia (bod, 1);

  if (bod->inverse) MX_Destroy (bod->inverse);

  K = tangent_stiffness (bod, 1);

  bod->inverse = MX_Add (1.0, M, step*step, K, NULL); /* TODO: figure out alpha and beta scaling */

  MX_Inverse (bod->inverse, bod->inverse);

  MX_Destroy (K);
}

/* total lagrangian initialise dynamic time stepping */
static void TL_dynamic_init (BODY *bod)
{
  if (!bod->M) /* once */
  {
    bod->M = diagonal_inertia (bod, bod->scheme != SCH_DEF_EXP);
  }

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
    double vmi = DBL_MAX, vol, tcrit, eigmax;
    MESH *msh = FEM_MESH (bod);
    ELEMENT *emi = NULL;
    MX *IMK;

#if OMP
    int j, n;
    ELEMENT **pele = ompu_surfeles (msh, &n);
    omp_lock_t *lock = ompu_locks (1);
    #pragma omp parallel for private (vol) shared (pele, lock, msh, emi, vmi)
    for (j = 0; j < n; j ++)
    {
      vol = ELEMENT_Volume (msh, pele[j], 0);

      if (vol < vmi)
      {
	omp_set_lock (lock);
	emi = pele[j];
	vmi = vol;
	omp_unset_lock (lock);
      }
    }
    ompu_locks_free (lock, 1);
    free (pele);
#else
    ELEMENT *ele;
    int bulk;
    /* find element with smallest volume */
    for (ele = msh->surfeles, emi = NULL, bulk = 0, vmi = DBL_MAX; ele; )
    {
      vol = ELEMENT_Volume (msh, ele, 0);
      if (vol < vmi)
      {
	emi = ele;
	vmi = vol;

        /* TODO: actually find an element with the smallest edge, rather than volume */
      }

      if (bulk) ele = ele->next;
      else if (ele->next) ele = ele->next;
      else ele = msh->bulkeles, bulk = 1;
    }
#endif

    IMK = element_inv_M_K (bod, msh, emi); /* element inv (M) * K */
    MX_Eigen (IMK, 1, &eigmax, NULL); /* maximal eigenvalue */
    MX_Destroy (IMK);
    ASSERT (eigmax > 0.0, ERR_BOD_MAX_FREQ_LE0);
    tcrit = 2.0 / sqrt (eigmax); /* limit of stability => t_crit <= 2.0 / omega_max */

    return tcrit;
  }
  else return DBL_MAX;
}

/* total lagrangian perform the initial half-step of the dynamic scheme */
static void TL_dynamic_step_begin (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  double half = 0.5 * step,
	*x = bod->inverse->x,
	*u0 = FEM_VEL0 (bod),
	*fext = FEM_FEXT (bod),
	*fint = FEM_FINT (bod),
	*q = bod->conf,
	*u = bod->velo,
	*e = u + n,
	*f, *g;

  update_fields (bod, time + half); /* update nodal field variables */

  ERRMEM (f = malloc (sizeof (double [n])));

  blas_dcopy (n, u, 1, u0, 1); /* save u (t) */

  switch (bod->scheme)
  {
  case SCH_DEF_EXP:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    TL_dynamic_force (bod, time+half, step, fext, fint, f);  /* f = fext (t+h/2) - fint (q(t+h/2)) */
    if (bod->damping > 0.0)
    {
      if (bod->K) MX_Destroy (bod->K);
      bod->K = tangent_stiffness (bod, 1);
      MX_Matvec (-bod->damping, bod->K, u, 1.0, f); /* f -= damping K u (t) */
    }
    for (g = f; u < e; u ++, x ++, g ++) (*u) += step * (*x) * (*g); /* u(t+h) = u(t) + inv (M) * h * f */
  }
  break;
  case SCH_DEF_LIM:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    TL_dynamic_force (bod, time+half, step, fext, fint, f);  /* f = fext (t+h/2) - fint (q(t+h/2)) */
    TL_dynamic_inverse (bod, step, NULL); /* A = M + (h*h/4) * K */
    if (bod->damping > 0.0) MX_Matvec (-bod->damping, bod->K, u, 1.0, f); /* f -= damping K u (t) */
    MX_Matvec (step, bod->inverse, f, 1.0, u); /* u(t+h) = u(t) + inv (A) * h * f */
  }
  break;
  default:
  break;
  }

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

  ERRMEM (r = malloc (sizeof (double [n])));

  fem_constraints_force (bod, r); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  blas_daxpy (n, 1.0, r, 1, fext, 1);  /* fext += r */

  switch (bod->scheme)
  {
  case SCH_DEF_EXP:
  {
    for (ir = r; iu < e; iu ++, x ++, ir ++) (*iu) += step * (*x) * (*ir); /* u(t+h) += inv (M) * h * r */
    blas_daxpy (n, half, u, 1, q, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
  }
  break;
  case SCH_DEF_LIM:
  {
    MX_Matvec (step, bod->inverse, r, 1.0, u); /* u(t+h) += h * inv (M) * force */
    blas_daxpy (n, half, u, 1, q, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
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
  if (!bod->M && !bod->K)
  {
    TL_static_inverse (bod, bod->dom->step);
  }
}

/* total lagrangian perform the initial half-step of the static scheme */
static void TL_static_step_begin (BODY *bod, double time, double step)
{
  double *f;

  update_fields (bod, time+step); /* update nodal field variables */

  ERRMEM (f = malloc (sizeof (double [bod->dofs])));
  TL_static_inverse (bod, step); /* compute inverse of static tangent operator */
  TL_static_force (bod, time+step, step, FEM_FEXT(bod), FEM_FINT(bod), f);  /* f(t+h) = fext (t+h) - fint (q(t+h)) */
  MX_Matvec (step, bod->inverse, f, 0.0, bod->velo); /* u(t+h) = inv (A) * h * f(t+h) */
  free (f);
}

/* total lagrangian perform the final half-step of the static scheme */
static void TL_static_step_end (BODY *bod, double time, double step)
{
  double *r;

  ERRMEM (r = malloc (sizeof (double [bod->dofs])));
  fem_constraints_force (bod, r); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  blas_daxpy (bod->dofs, 1.0, r, 1, FEM_FEXT (bod), 1);  /* fext += r */
  MX_Matvec (step, bod->inverse, r, 1.0, bod->velo); /* u(t+h) += inv (A) * h * r */
  blas_daxpy (bod->dofs, step, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t) + h * u(t+h) */
  free (r);
}

/* =================== BODY COROTATIONAL =================== */

/* compute surface integral = INT { skew [N] A [X + Shapes (X) q] } */
static void BC_surface_integral (BODY *bod, MESH *msh, double *conf, int num, double *A, double *integral)
{
  double q [4][3], refn [4][3], curn [4][3], B [3], C [3], *N, *shapes, *Y, *i, *e;
  FACE *fac;
  int j;

  for (i = integral, e = i + 3 * num; i < e; i += 3) SET (i, 0);

#if OMP
  int integral_size = 3 * num;
  double *local_integral;
  int threads;
  int k, m;
  #pragma omp parallel
  #pragma omp master
  threads = omp_get_num_threads();
  FACE **pfac = ompu_faces (msh, &m);
  ERRMEM (local_integral = MEM_CALLOC (threads * sizeof (double [integral_size])));
  #pragma omp parallel
  {
  double *integral = &local_integral[omp_get_thread_num() * integral_size];
  double *e = integral + integral_size;
  #pragma omp for private (q,refn,curn,B,C,N,shapes,Y,fac,i,j)
  for (k = 0; k < m; k ++)
#else
  for (fac = msh->faces; fac; fac = fac->n)
#endif
  {
#if OMP
    fac = pfac[k];
#endif
    face_nodes (msh->ref_nodes, fac->type, fac->nodes, refn);
    face_displacements (conf, fac, q);
    for (j = 0; j < fac->type; j ++) ADD (refn [j], q [j], curn [j]);

    INTEGRAL2D_BEGIN (fac->type) /* defines point and weight */
    {
      shapes = FACE_SHAPES (fac, __k__);
      N = FACE_NORMAL (fac, __k__);

      SET (B, 0);
      for (j = 0; j < fac->type; j ++)
	ADDMUL (B, shapes [j], curn [j], B);
      SCALE (B, N[3]);

      for (i = integral, Y = A; i < e; i += 3, Y += 9)
      {
	NVMUL (Y, B, C);
	PRODUCTADD (N, C, i);
      }
    }
    INTEGRAL2D_END ()
  }
#if OMP
  }
  double *li = local_integral;
  for (int i = 0; i < threads; i ++)
    for (j = 0; j < integral_size; j ++, li ++)
      integral [j] += *li;
  free (local_integral);
  free (pfac);
#endif
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

    ASSERT (lapack_dposv ('U', 3, 1, ddJ, 3, dJ, 3) == 0, ERR_FEM_ROT_SINGULAR_JACOBIAN);

    SUB (O, dJ, O);
    error = sqrt (DOT (dJ, dJ) / (1.0 + DOT (O, O)));
  }
  while (error > IMP_EPS && ++ iter < MAX_ITERS);

  EXPMAP (O, A);
  NNCOPY (R, A+9);
  NNMUL (A, A+9, R); /* R = exp (O) R */

#if 0
  printf ("O = %g, %g, %g after %d iterations\n", O[0], O[1], O[2], iter);
#endif

  ASSERT (iter < MAX_ITERS, ERR_FEM_ROT_NEWTON_DIVERGENCE);
}

/* fint = R K R' [(I-R)Z + q] */
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

/* energy = 0.5 * [(I-R)Z + q] R K R' [(I-R)Z + q] */
static double BC_internal_energy (BODY *bod)
{
  double *a, *b, *x, *y, *d, *z, (*Z) [3], (*e) [3], Y [3], intene;
  double *R = FEM_ROT (bod), *q = bod->conf;
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

  intene = 0.5 * blas_ddot (n, a, 1, b, 1);

  free (a);

  return intene;
}

/* compute u = alpha * R A R' b + beta * u  */
static void BC_matvec (double alpha, MX *A, double *R, double *b, double beta, double *u)
{
  double *x, *y, *z;
  int n = A->n;

  ERRMEM (x = MEM_CALLOC (2 * sizeof (double [n])));

  for (y = x, z = x + n; y < z; y += 3, b += 3)
  {
    TVMUL (R, b, y);
  }

  MX_Matvec (alpha, A, x, 0.0, y);

  if (beta != 1.0) blas_dscal (n, beta, u, 1);

  for (z = y + n; y < z; y += 3, u += 3)
  {
    NVADDMUL (u, R, y, u);
  }

  free (x);
}

/* body co-rotational initialise dynamic time stepping */
static void BC_dynamic_init (BODY *bod)
{
  if (!bod->M && !bod->K)
  {
    bod->M = diagonal_inertia (bod, bod->scheme != SCH_DEF_EXP);

    bod->K = tangent_stiffness (bod, 1);

    if (bod->scheme == SCH_DEF_EXP)
    {
      double *x, *y;

      bod->inverse = MX_Copy (bod->M, NULL);

      for (x = bod->inverse->x, y = x + bod->dofs; x < y; x ++)
      {
	ASSERT (*x > 0.0, ERR_FEM_MASS_NOT_SPD);
	(*x) = 1.0 / (*x); /* invert diagonal */
      }
    }
    else
    {
      double step = bod->dom->step;

      /* calculate initial tangent operator A(0) = M + (damping*h/2 + h*h/4) K(q(0)) */
      bod->inverse = MX_Add (1.0, bod->M, 0.5*bod->damping*step + 0.25*step*step, bod->K, NULL);

      MX_Inverse (bod->inverse, bod->inverse);
    }
  }
}

/* body co-rotational estimate critical step for the dynamic scheme */
static double BC_dynamic_critical_step (BODY *bod)
{
  if (bod->cristep0 == 0.0)
  {
    bod->cristep0 = TL_dynamic_critical_step (bod);
  }

  return bod->cristep0;
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

  ERRMEM (b = malloc (sizeof (double [n])));

  switch (bod->scheme)
  {
  case SCH_DEF_EXP:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    BC_update_rotation (bod, msh, q, R); /* R1 = R(q(t+h/2)) */

    external_force (bod, time+half, step, fext); /* fext = fext (t+h/2) */
    BC_internal_force (bod, R, q, fint); /* fint = R1 K R1' [(I-R1)Z + q(t+h/2)] */
    blas_dcopy (n, fext, 1, b, 1);
    blas_daxpy (n, -1.0, fint, 1, b, 1); /* b = h (fext - fint) */
    if (bod->damping > 0.0) BC_matvec (-bod->damping, bod->K, R, u, 1.0, b); /* b -= damping K u (t) */

    MX_Matvec (step, bod->inverse, b, 1.0, u); /* u(t+h) = u(t) + inv (M) * b */
  }
  break;
  case SCH_DEF_LIM:
  {
    blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    BC_update_rotation (bod, msh, q, R); /* R1 = R(q(t+h/2)) */

    external_force (bod, time+half, step, fext); /* fext = fext (t+h/2) */
    BC_internal_force (bod, R, q, fint); /* fint = R1 K R1' [(I-R1)Z + q(t+h/2)] */
    blas_dcopy (n, fext, 1, b, 1);
    blas_daxpy (n, -1.0, fint, 1, b, 1); /* b = h (fext - fint) */
    if (bod->damping > 0.0) BC_matvec (-bod->damping, bod->K, R, u, 1.0, b); /* b -= damping K u (t) */

    BC_matvec (step, bod->inverse, R, b, 1.0, u); /* u(t+h) = u(t) + inv (M + (h*h/4) R1 K R1') b */
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
  {
    BC_matvec (step, bod->inverse, R, r, 1.0, u); /* u(t+h) += inv (A) * h * r */

    blas_daxpy (n, half, u, 1, q, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
    BC_update_rotation (bod, msh, q, R); /* R(t+h) = R(q(t+h)) */
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
  if (!bod->M && !bod->K)
  {
    double step = bod->dom->step;

    bod->M = diagonal_inertia (bod, 1);

    bod->K = tangent_stiffness (bod, 1);

    bod->inverse = MX_Add (1.0, bod->M, step*step, bod->K, NULL); /* TODO: figure out alpha and beta scaling */

    MX_Inverse (bod->inverse, bod->inverse);
  }
}

/* body co-rotational perform the initial half-step of the static scheme */
static void BC_static_step_begin (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  double *fext = FEM_FEXT (bod),
	 *fint = FEM_FINT (bod),
	 *R = FEM_ROT (bod),
	 *q = bod->conf,
	 *u = bod->velo,
         *b;

  ERRMEM (b = MEM_CALLOC (sizeof (double [n])));

  external_force (bod, time+step, step, fext); /* fext = fext (t+h) */
  BC_internal_force (bod, R, q, fint); /* fint = R K R' [(I-R)Z + q(t)] */
  blas_daxpy (n, step, fext, 1, b, 1);
  blas_daxpy (n, -step, fint, 1, b, 1); /* b = h (fext - fint) */

  BC_matvec (1.0, bod->inverse, R, b, 0.0, u); /* u(t+h) = inv (A) b */

  free (b);
}

/* body co-rotational perform the final half-step of the static scheme */
static void BC_static_step_end (BODY *bod, double time, double step)
{
  MESH *msh = FEM_MESH (bod);
  int n = bod->dofs;
  double *fext = FEM_FEXT (bod),
	 *R = FEM_ROT (bod),
	 *q = bod->conf,
	 *u = bod->velo,
	 *r;

  ERRMEM (r = malloc (sizeof (double [n])));

  fem_constraints_force (bod, r); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  blas_daxpy (n, 1.0, r, 1, fext, 1);  /* fext += r */

  BC_matvec (step, bod->inverse, R, r, 1.0, u); /* u(t+h) += inv (A) * h * r */

  blas_daxpy (n, step, u, 1, q, 1); /* q (t+h) = q(t) + h * u(t+h) */
  BC_update_rotation (bod, msh, q, R); /* R(t+h) = R(q(t+h)) */

  free (r);
}

/* =================== REDUCED ORDER =================== */

/* x = R x */
static void RO_rotate_forward (double *R, double *x, int n)
{
  double *y = x + n, z [3];

  for (;x < y; x += 3)
  {
    COPY (x, z);
    NVMUL (R, z, x);
  }
}

/* x = R' x */
static void RO_rotate_backward (double *R, double *x, int n)
{
  double *y = x + n, z [3];

  for (;x < y; x += 3)
  {
    COPY (x, z);
    TVMUL (R, z, x);
  }
}

/* qm = REq - (I-R)Z */
static void RO_lift_conf (BODY *bod, MX *E, MESH *msh, double *R, double *q, double *qm)
{
  double (*Z) [3], Y [3], *x, *y;

  MX_Matvec (1.0, E, q, 0.0, qm); /* d0 = E q */

  for (x = qm, y = qm+E->m, Z = msh->ref_nodes; x < y; x += 3, Z ++)
  {
    COPY (x, Y);
    NVMUL (R, Y, x); /* d = R d0 */

    NVMUL (R, Z[0], Y);
    SUB (Z[0], Y, Y);
    SCC (Y, x); /* qm = d - (I-R)Z */
  }

  /* FIXME --> 'BC-RO' case needs to be taken into account here;
               otherwise in 'READ' mode the rigid motion is lost;
	       this has to do with mesh mass scaling in projections; */
}

/* q = project_onto_E_in_the_M_norm (d = R'[(I-R)Z+qm]) */
static void RO_project_conf (BODY *bod, MX *E, MESH *msh, double *tmp, double *R, double *qm, double *q)
{
  double (*Z) [3], Y [3], *x, *y, *z;

  blas_dcopy (E->m, qm, 1, tmp, 1);

  for (x = tmp, y = tmp+E->m, Z = msh->ref_nodes, z = FEM_MESH_MASS (bod); x < y; x += 3, z += 3, Z ++)
  {
    NVMUL (R, Z[0], Y);
    SUB (Z[0], Y, Y);
    ACC (Y, x); /* d = (I-R)Z + qm */

    COPY (x, Y);
    TVMUL (R, Y, x); /* R' d */

    if (bod->form == BODY_COROTATIONAL_MODAL)
    {
      HADAMARD (x, z, x); /* tmp = M R' d */
    }
  }

  MX_Matvec (1.0, MX_Tran (E), tmp, 0.0, q); /* q = E' M R' d */
}

/* u = E'M R' um */
static void RO_project_velo (BODY *bod, MX *E, double *tmp, double *R, double *um, double *u)
{
  double Y [3], *x, *y, *z;

  blas_dcopy (E->m, um, 1, tmp, 1);

  for (x = tmp, y = tmp+E->m, z = FEM_MESH_MASS (bod); x < y; x += 3, z += 3)
  {
    COPY (x, Y);
    TVMUL (R, Y, x); /* R' um */

    if (bod->form == BODY_COROTATIONAL_MODAL)
    {
      HADAMARD (x, z, x); /* tmp = M R' um */
    }
  }

  MX_Matvec (1.0, MX_Tran (E), tmp, 0.0, u); /* u = E' M R' um */
}

/* fint = diagonal (K) q {in the reduced space} */
static void RO_internal_force (BODY *bod, double *q, double *fint)
{
#if 0
  double *eval = bod->K->x;
  int i, n = bod->dofs;

  for (i = 0; i < n; i ++)
  {
    fint [i] = eval [i] * q [i];
  }
#else
  MX_Matvec (1.0, bod->K, q, 0.0, fint);
#endif
}

/* initialize reduced velocity at time == 0.0 */
static void RO_initialize_velocity (BODY *bod)
{
  /* project initial velocity onto the reduced basis */
  double *velo = FEM_MESH_VELO (bod), *M = FEM_MESH_MASS (bod), *p;
  MX *E = bod->evec;
  int i;

  if (bod->form == BODY_COROTATIONAL_MODAL)
  {
    ERRMEM (p = malloc (sizeof (double [E->m])));
    for (i = 0; i < E->m; i ++) p [i] = M[i]*velo[i]; /* p = full initial momentum */
    MX_Matvec (1.0, MX_Tran (E), p, 0.0, bod->velo); /* p_reduced(0) = E' p_full(0) */
    /* since reduced mass is identity, reduced velocity is by magnitude the same as momentum */
    free (p);
  }
  else if (bod->form == BODY_COROTATIONAL_REDUCED_ORDER)
  {
    MX_Matvec (1.0, MX_Tran (E), velo, 0.0, bod->velo); /* u_reduced(0) = E' u_full(0) */
  }
}

/* reduced order initialise dynamic time stepping */
static void RO_dynamic_init (BODY *bod)
{
  if (!bod->M && !bod->K)
  {
    DOM *dom = bod->dom;
    double step = dom->step,
	   time = dom->time;
    MX *E = bod->evec;
    MX *M = NULL;
    int i;

    ASSERT_TEXT (bod->evec, "Reduced base does not exist. Call MODAL_ANALYSIS after creating a reduced order body!");

    if (bod->form == BODY_COROTATIONAL_MODAL)
    {
      /* bod->M = E' M E  = I */
      bod->M = MX_Identity (MXCSC, E->n);

      /* bod->K = E' K E = diag (lambdas) */
      bod->K = MX_Identity (MXCSC, E->n);
      for (i = 0; i < E->n; i ++) bod->K->x [i] = i < 6 ? 0.0 : bod->eval [i];

      bod->inverse = MX_Identity (MXCSC, E->n);
      /* calculate initial tangent operator A(0) = E' (M + (damping*h/2 + h*h/4) K(q(0))) E */
      for (i = 0; i < E->n; i ++) bod->inverse->x [i] = 1.0 / (1.0 + (0.5*bod->damping*step + 0.25*step*step)*bod->K->x[i]);
    }
    else if (bod->form == BODY_COROTATIONAL_REDUCED_ORDER)
    {
      /* bod->M = E' M E  = I */
      M = diagonal_inertia (bod, 1);
      MX *EtM = MX_Matmat (1.0, MX_Tran(E), M, 0.0, NULL);
      bod->M = MX_Matmat (1.0, EtM, E, 0.0, NULL);
      MX_Destroy (EtM);

      /* bod->K = E' K E = diag (lambdas) */
      MX *K = tangent_stiffness (bod, 1);
      MX *EtK = MX_Matmat (1.0, MX_Tran(E), K, 0.0, NULL);
      bod->K = MX_Matmat (1.0, EtK, E, 0.0, NULL);
      MX_Destroy (EtK);
      MX_Destroy (K);

      double step = bod->dom->step;

      /* calculate initial tangent operator A(0) = M + (damping*h/2 + h*h/4) K(q(0)) */
      bod->inverse = MX_Add (1.0, bod->M, 0.5*bod->damping*step + 0.25*step*step, bod->K, NULL);

      MX_Inverse (bod->inverse, bod->inverse);
    }
    else
    {
      ASSERT_TEXT (0, "This should not happen --> please report a bug");
    }

    /* store diagonal mesh space inertia matrix */
    if (!M) M = diagonal_inertia (bod, 1);
    double *x = FEM_MESH_MASS (bod);
    blas_dcopy (E->m, M->x, 1, x, 1);

    /* clean up */
    MX_Destroy (M); 
    
    /* initialize reduced velocity */
    if (time == 0.0) RO_initialize_velocity (bod);
  }
}

/* reduced order estimate critical step for the dynamic scheme */
static double RO_dynamic_critical_step (BODY *bod)
{
  return DBL_MAX;
}

/* reduced order perform the initial half-step of the dynamic scheme */
static void RO_dynamic_step_begin (BODY *bod, double time, double step)
{
  MX *E = bod->evec;
  MESH *msh = FEM_MESH (bod);
  int n = bod->dofs,
      nm = MESH_DOFS (msh);
  double half = 0.5 * step,
	*u0 = FEM_VEL0 (bod),
	*u0m = FEM_MESH_VEL0 (bod),
	*R = FEM_ROT (bod),
	*fext = FEM_FEXT (bod),
	*fint = FEM_FINT (bod),
	*qm = FEM_MESH_CONF (bod),
	*um = FEM_MESH_VELO (bod),
	*u = bod->velo,
	*q = bod->conf,
	*b, *v, *tmp;

  ERRMEM (b = MEM_CALLOC (sizeof (double [2*n+nm])));
  v = b + n;
  tmp = v + n;

  blas_dcopy (n, u, 1, u0, 1); /* save u (t) */
  blas_dcopy (nm, um, 1, u0m, 1);

  blas_daxpy (nm, half, um, 1, qm, 1); /* qm(t+h/2) = qm(t) + (h/2)um(t) */
  BC_update_rotation (bod, msh, qm, R); /* R1 = R(qm(t+h/2)) */
  RO_project_conf (bod, E, msh, tmp, R, qm, q); /* q(t+h/2) = proj_M_E (qm(t+h/2)) */

  external_force (bod, time+half, step, tmp);
  RO_rotate_backward (R, tmp, nm); /* f */
  MX_Matvec (1.0, MX_Tran (E), tmp, 0.0, fext); /* fext = E'R'f */
  RO_internal_force (bod, q, fint); /* fint = K q */
  blas_dcopy (n, fext, 1, b, 1);
  blas_daxpy (n, -1.0, fint, 1, b, 1); /* b = fext - fint */
  if (bod->damping > 0.0)
  {
    RO_project_velo (bod, E, tmp, R, um, v); /* v = E'MR1'um <= half-step rotation is used */
    MX_Matvec (-bod->damping, bod->K, v, 1.0, b); /* b -= damping K u (t) */
  }

  MX_Matvec (step, bod->inverse, b, 1.0, u); /* u(t+h) = u(t) + h inv (A) b */

  blas_dcopy (n, u, 1, b, 1);
  blas_daxpy (n, -1.0, u0, 1, b, 1); /* du */
  MX_Matvec (1.0, E, b, 0.0, tmp);
  RO_rotate_forward (R, tmp, nm); /* dum */
  blas_daxpy (nm, 1.0, tmp, 1, um, 1); /* um = u0m + REdu */

  free (b);
}

/* reduced order perform the final half-step of the dynamic scheme */
static void RO_dynamic_step_end (BODY *bod, double time, double step)
{
  MX *E = bod->evec;
  MESH *msh = FEM_MESH (bod);
  int n = bod->dofs,
      nm = MESH_DOFS (msh);
  double half = 0.5 * step,
	*R = FEM_ROT (bod),
	*fext = FEM_FEXT (bod),
	*qm = FEM_MESH_CONF (bod),
	*um = FEM_MESH_VELO (bod),
	*u0m = FEM_MESH_VEL0 (bod),
	*u0 = FEM_VEL0 (bod),
	*u = bod->velo,
	*q = bod->conf,
	*r, *tmp;

  ERRMEM (r = malloc (sizeof (double [n+nm])));
  tmp = r + n;

  fem_constraints_force (bod, r); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  blas_daxpy (n, 1.0, r, 1, fext, 1);  /* fext += r */

  MX_Matvec (step, bod->inverse, r, 1.0, u); /* u(t+h) += inv (A) * h * r */

  blas_dcopy (n, u, 1, r, 1);
  blas_daxpy (n, -1.0, u0, 1, r, 1); /* du */
  MX_Matvec (1.0, E, r, 0.0, tmp);
  RO_rotate_forward (R, tmp, nm); /* dum */
  blas_dcopy (nm, u0m, 1, um, 1);
  blas_daxpy (nm, 1.0, tmp, 1, um, 1); /* um = u0m + REdu => this produces "good" deformation */

  blas_dcopy (nm, qm, 1, tmp, 1); /* qm(t+h/2) */
  blas_daxpy (nm, half, um, 1, tmp, 1); /* qm(t+h) = qm(t+h/2) + (h/2)um(t+h) */
  BC_update_rotation (bod, msh, tmp, R); /* R(t+h) = R(qm(t+h)) */

  if (bod->form == BODY_COROTATIONAL_MODAL)
  {
    RO_project_velo (bod, E, tmp, R, um, r); /* r[i>=6] store the "good" deformation components */
    for (int i = 6; i < n; i ++) u[i] = r[i]; /* while u[i<6] store the "good" rotation components */
  }
  else
  {
    RO_project_velo (bod, E, tmp, R, um, u); /* project um onto reduced u */
  }

  MX_Matvec (1.0, E, u, 0.0, um); /* resolve back um from u, using non-incremental formula */
  RO_rotate_forward (R, um, nm); /* now um(t+h) combines "good" rotation and deformation */

  blas_daxpy (nm, half, um, 1, qm, 1); /* qm(t+h) = qm(t+h/2) + (h/2)um(t+h) */
  RO_project_conf (bod, E, msh, tmp, R, qm, q); /* qm(t+h) = resolve_back (proj(qm(t+h))) */

  free (r);
}

/* reduced order initialise static time stepping */
static void RO_static_init (BODY *bod)
{
  RO_dynamic_init (bod);
}

/* reduced order perform the initial half-step of the static scheme */
static void RO_static_step_begin (BODY *bod, double time, double step)
{
  MESH *msh = FEM_MESH (bod);
  int nm = MESH_DOFS (msh),
      n = bod->dofs,
      i;
  double *u0m = FEM_MESH_VEL0 (bod),
	 *um = FEM_MESH_VELO (bod),
	 *u0 = FEM_VEL0 (bod),
	 *u = bod->velo;

  /* zero velocity */
  for (i = 0; i < nm; i ++) { u0m[i] = 0.0; um[i] = 0.0; }
  for (i = 0; i < n; i ++) { u0[i] = 0.0; u[i] = 0.0; }

  RO_dynamic_step_begin (bod, time, step);
}

/* reduced order perform the final half-step of the static scheme */
static void RO_static_step_end (BODY *bod, double time, double step)
{
  RO_dynamic_step_end (bod, time, step);
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
  int bulk, n, k, m;
  TRISURF *surf;
  ELEMENT *ele;
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

  SET *eledel = NULL;

  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    if (!ele->dom) /* schedule element deletion */
    {
      SET_Insert (NULL, &eledel, ele, NULL);
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

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }

  MESH_Delete_Elements (msh, eledel);
  SET_Free (NULL, &eledel);

  /* 3. Transform global ele->dom points into local element coordinates */

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
  double vertices [24], planes [72], *pla;
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
  tri = cvi (cvx->cur, cvx->nver, pla, cvx->nfac, vertices, n, planes, k, REGULARIZED, &m, NULL, NULL);
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

/* map m1 values onto m2 values */
void FEM_Map_State (MESH *m1, double *q1, double *u1, MESH *m2, double *q2, double *u2)
{
  double shapes [MAX_NODES], val [MAX_NODES][3], point [3];
  double extents [6], (* nod) [3], (*end) [3];
  ELEMENT *ele;
  int i, n;
  KDT *kd;

  /* creake input mesh based kd-tree */
  kd = KDT_Create (m1->nodes_count, (double*)m1->ref_nodes, 0.0);
  for (ele = m1->surfeles; ele; ele = ele->next)
  { ELEMENT_Ref_Extents (m1, ele, extents); KDT_Drop (kd, extents, ele); }
  for (ele = m1->bulkeles; ele; ele = ele->next)
  { ELEMENT_Ref_Extents (m1, ele, extents); KDT_Drop (kd, extents, ele); }

  /* map m2 nodal values */
  for (nod = m2->ref_nodes, end = nod + m2->nodes_count; nod != end; nod ++, q2 +=3, u2 += 3)
  {
    KDT *q = KDT_Pick (kd, nod [0]);
    ASSERT_DEBUG (q, "Inconsistent kd-tree query");
    ELEMENT **ptr = (ELEMENT**) q->data, **qtr = ptr + q->n;
    for (; ptr != qtr; ptr ++)
      if (ELEMENT_Contains_Point (m1, *ptr, nod [0], 1)) break;
    ASSERT_DEBUG (ptr != qtr, "Element containing a referential point has not been found");

    referential_to_local (m1, *ptr, nod [0], point);
    n = element_shapes ((*ptr)->type, point, shapes);

    element_nodal_values (q1, *ptr, val);
    SET (q2, 0); for (i = 0; i < n; i ++) { ADDMUL (q2, shapes [i], val [i], q2); }

    element_nodal_values (u1, *ptr, val);
    SET (u2, 0); for (i = 0; i < n; i ++) { ADDMUL (u2, shapes [i], val [i], u2); }
  }

  /* clean up */
  KDT_Destroy (kd);
}

/* directly map m1 values onto m2 values */
void FEM_Map_State_Direct (MESH *m1, double *q1, double *u1, MESH *m2, double *q2, double *u2, int *mapping)
{
  for (int i = 0; i < m2->nodes_count; i ++)
  {
    double *vq2 = &q2[3*i],
           *vu2 = &u2[3*i],
	   *vq1 = &q1[3*mapping[i]],
	   *vu1 = &u1[3*mapping[i]];

    COPY (vq1, vq2);
    COPY (vu1, vu2);
  }
}

/* ================== Utilities ==================== */

/* modal motion update callback */
static void modal_motion (BODY *bod, ELEMENT *ele, double *point, double *x)
{
  if (bod->msh)
  {
    local_to_spatial (bod->msh, ele, point, x);
  }
}

/* ================== INTERFACE ==================== */

/* create FEM internals for a body */
void FEM_Create (FEMFORM form, MESH *msh, SHAPE *shp, BULK_MATERIAL *mat, BODY *bod)
{
  /* compute shape characteristics */
  SHAPE_Char (shp, 1, &bod->ref_volume, bod->ref_center, bod->ref_tensor);
  bod->ref_mass = bod->ref_volume * mat->density;
  SCALE9 (bod->ref_tensor, mat->density);

  if (msh) /* the given mesh is assumed to properly contain the shape */
  {
    SHAPE msh_shp = {SHAPE_MESH, msh, NULL};
    BOX **msh_boxes, **shp_boxes, **box;
    SGP *msh_sgp, *shp_sgp, *sgp, *sge;
    int msh_nsgp, shp_nsgp, shp_nsgpall;
    BOX_Extents_Update update;
    ELEMENT *ele;
    MEM boxmem;
  
    SHAPE_Update (shp, NULL, NULL); /* restore reference shape */
    MESH_Update (msh, NULL, NULL, NULL); /* same for the background mesh */

    msh_nsgp = msh->surfeles_count + msh->bulkeles_count;
    ERRMEM (msh_sgp = sgp = MEM_CALLOC (msh_nsgp * sizeof (SGP)));
    for (ele = msh->surfeles; ele; ele = ele->next, sgp ++) sgp->shp = &msh_shp, sgp->gobj = ele, sgp->kind = GOBJ_ELEMENT;
    for (ele = msh->bulkeles; ele; ele = ele->next, sgp ++) sgp->shp = &msh_shp, sgp->gobj = ele, sgp->kind = GOBJ_ELEMENT;
    shp_sgp = SGP_Create (shp, &shp_nsgp, &shp_nsgpall); /* bacause shape is CONVEX based SGPs will cover the complete volume */
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

  /* surface faces integration data */
  create_face_integration_data (msh);

  /* allocate dofs */
  switch (form)
  {
    case TOTAL_LAGRANGIAN:
    {
      bod->dofs = MESH_DOFS (msh);
      ERRMEM (bod->conf = MEM_CALLOC (6 * bod->dofs * sizeof (double))); /* configuration, velocity, previous velocity, fext, fint, fbod */
      bod->velo = bod->conf + bod->dofs;
      if (mat->nfield) ERRMEM (bod->field = MEM_CALLOC (mat->nfield * (bod->dofs / 3) * sizeof (double))); /* field variables */
      else bod->field = NULL;
    }
    break;
    case BODY_COROTATIONAL:
    {
      bod->dofs = MESH_DOFS (msh);
      ERRMEM (bod->conf = MEM_CALLOC ((6 * bod->dofs + 9) * sizeof (double))); /* configuration, rotation, velocity, previous velocity, fext, fint, fbod */
      bod->velo = bod->conf + bod->dofs + 9;
      double *R = bod->conf + bod->dofs;
      IDENTITY (R);
      bod->field = NULL; /* linear material only */
    }
    break;
    case BODY_COROTATIONAL_MODAL:
    case BODY_COROTATIONAL_REDUCED_ORDER:
    {
      ASSERT_DEBUG (bod->evec && bod->eval, "Modal analysis results must exist!");

      /* ------------------------------------------------------------------------------------------------------------
       * bod->dofs = bod->evec->n
       * bod->con = {reduced configuration} {rotation} {full configuration} {mesh space diagonal mass}
       * bod->velo = {reduced => velocity, previous velocity, fext, fint} {full => velocity, previous velocity, fbod}
       * ------------------------------------------------------------------------------------------------------------ */

      MX *E = bod->evec;
      bod->dofs = E->n;
      ERRMEM (bod->conf = MEM_CALLOC (((E->n+9+2*E->m) + (4*E->n+3*E->m)) * sizeof (double)));
      bod->velo = bod->conf + (E->n+9+2*E->m);
      double *R = bod->conf + E->n;
      IDENTITY (R); /* initialise rotation */
      bod->field = NULL; /* linear material only */
    }
    break;
  }

  /* save formulation */
  bod->form = form;

  /* save rought mesh if needed */
  if (msh != shp->data)
  {
    bod->msh = msh;

    /* simplex integrated volume ==
     * shape functions integrated volume test */
    test_volume_integral (msh, bod->ref_volume, bod->id);
  }

  /* allocate bulk material
   * states at integration points */
  allocate_element_states (msh, mat);
}

/* overwrite state */
void FEM_Overwrite_State (BODY *bod, double *q, double *u)
{
  blas_dcopy (MESH_DOFS(FEM_MESH(bod)), q, 1, bod->conf, 1); /* always in the mesh space */
  blas_dcopy (bod->dofs, u, 1, bod->velo, 1); /* can be mesh or reduced space */

  if (bod->form >= BODY_COROTATIONAL_MODAL) /* overwrite the mesh space */
  {
    FEM_Post_Read (bod);
  }
}

/* set initial rigid motion velocity */
void FEM_Initial_Velocity (BODY *bod, double *linear, double *angular)
{
  MESH *msh = FEM_MESH (bod);
  double (*nod) [3] = msh->ref_nodes,
	 (*end) [3] = nod + msh->nodes_count,
	 *velo = FEM_MESH_VELO (bod),
	 *X0 = bod->ref_center,
	 A [3];

  for (; nod < end; nod ++, velo += 3)
  {
    if (linear) { COPY (linear, velo); }
    else { SET (velo, 0.0); }

    if (angular)
    {
      SUB (nod [0], X0, A);
      PRODUCTADD (angular, A, velo);
    }
  }
}

/* set rigid motion */
void FEM_From_Rigid (BODY *bod, double *rotation, double *position, double *angular, double *linear)
{
  ASSERT_TEXT (bod->kind == FEM, "This is not an FEM body!");

  MESH *msh = FEM_MESH (bod);
  double (*ref) [3] = msh->ref_nodes,
	 (*cur) [3] = msh->cur_nodes;
  int m = msh->nodes_count, n;
  double *center = bod->ref_center,
         *conf = bod->conf;
  double A[3];

  for (n = 0; n < m; n ++, conf += 3)
  {
    SUB (ref[n], center, A);
    NVADDMUL (position, rotation, A, cur[n]);
    SUB (cur[n], ref[n], conf);
  }

  FEM_Initial_Velocity (bod, linear, angular);
}

/* initialise dynamic time stepping */
void FEM_Dynamic_Init (BODY *bod)
{
  short noubf = 1;
  if (bod->M) noubf = 0; /* unit body force already computed */

  bod->cristep0 = 0.0; /* used by BODY_COROTATIONAL */

  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
      TL_dynamic_init (bod);
      break;
    case BODY_COROTATIONAL:
      BC_dynamic_init (bod);
      break;
    case BODY_COROTATIONAL_MODAL:
    case BODY_COROTATIONAL_REDUCED_ORDER:
      RO_dynamic_init (bod);
      break;
  }

  if (noubf) unit_body_force (bod); /* compute once (after initialization so that buffers are right for the RO model) */
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
    case BODY_COROTATIONAL_MODAL:
    case BODY_COROTATIONAL_REDUCED_ORDER:
      return RO_dynamic_critical_step (bod);
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
    case BODY_COROTATIONAL_MODAL:
    case BODY_COROTATIONAL_REDUCED_ORDER:
      RO_dynamic_step_begin (bod, time, step);
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
    case BODY_COROTATIONAL_MODAL:
    case BODY_COROTATIONAL_REDUCED_ORDER:
      RO_dynamic_step_end (bod, time, step);
      break;
  }

  for (; iu < ue; idq ++, iu ++, iu0 ++) *idq = half * ((*iu) + (*iu0)); /* dq = (h/2) * {u(t) + u(t+h)} */
  energy [EXTERNAL] += blas_ddot (n, dq, 1, fext, 1); /* XXX: may not be too good for the reduced order model (save q0 and copute dq = q1-q0) */
  if (bod->form == BODY_COROTATIONAL_MODAL)
  {
    energy [INTERNAL] = 0.0;
    for (double *q = bod->conf, *x = bod->K->x, *y = x + n; x < y; q ++, x ++)
      energy [INTERNAL] += 0.5 * (*x) * (*q) * (*q);
  }
  else if (bod->form == BODY_COROTATIONAL_REDUCED_ORDER)
  {
    double *q = bod->conf, *p;
    int n = bod->dofs;
    ERRMEM (p = malloc (n * sizeof (double)));
    MX_Matvec (1.0, bod->K, q, 0.0, p);
    energy[INTERNAL] = 0.5 * blas_ddot (n, q, 1, p, 1);
    free (p);
  }
  else if (bod->form == BODY_COROTATIONAL)
  {
    energy[INTERNAL] = BC_internal_energy (bod);
  }
  else
  {
    /* energy [INTERNAL] += blas_ddot (n, dq, 1, fint, 1);
     * computing internal energy like above may produce negative energy increments during
     * impacts since fint is computed at q(t+h/2) whereas dq includes impact correction; 
     * this is effect is present when the time integration step is excessively large */

    energy [INTERNAL] = internal_energy (bod);
  }

  free (dq);
}

/* initialise static time stepping */
void FEM_Static_Init (BODY *bod)
{
  if (!bod->M) unit_body_force (bod); /* compute once */

  switch (bod->form)
  {
    case TOTAL_LAGRANGIAN:
      TL_static_init (bod);
      break;
    case BODY_COROTATIONAL:
      BC_static_init (bod);
      break;
    case BODY_COROTATIONAL_MODAL:
    case BODY_COROTATIONAL_REDUCED_ORDER:
      RO_static_init (bod);
      break;
  }
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
    case BODY_COROTATIONAL_MODAL:
    case BODY_COROTATIONAL_REDUCED_ORDER:
      RO_static_step_begin (bod, time, step);
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
    case BODY_COROTATIONAL_MODAL:
    case BODY_COROTATIONAL_REDUCED_ORDER:
      RO_static_step_end (bod, time, step);
      break;
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
    element_nodal_values (FEM_MESH_CONF (bod), ele, q);

    COPY (X, x); for (i = 0; i < n; i ++) { ADDMUL (x, shapes [i], q [i], x); } /* x = X + N q */
  }
  else /* NULL implies nodal update (X is within msh->ref_nodes) */
  {
    int n = (node_t) X - msh->ref_nodes;
    double *conf = FEM_MESH_CONF (bod);
    double *q = &conf [n * 3];

    ADD (msh->ref_nodes [n], q, x);
  }
}

/* motion x = x (element, ref point, local point) */
void FEM_Cur_Point_Ext (BODY *bod, ELEMENT *ele, double *X, double *point, double *x)
{
  double shapes [MAX_NODES], q [MAX_NODES][3];
  int i, n;

  n = element_shapes (ele->type, point, shapes);
  element_nodal_values (FEM_MESH_CONF (bod), ele, q);

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

/* pull-forward v = {dx/dX} V (X, state) */
void FEM_Cur_Vector (BODY *bod, ELEMENT *ele, double *X, double *V, double *v)
{
  MESH *msh = FEM_MESH (bod);
  double point [3], F [9];

  if (ele == NULL)
  {
    ele = MESH_Element_Containing_Point (msh, X, 1);
  }

  ASSERT_TEXT (ele, "Element containing referential point was not found.\n"
                    "Please report this bug!\n");

  referential_to_local (msh, ele, X, point);
  deformation_gradient (bod, msh, ele, point, F);
  NVMUL (F, V, v);
}

/* push-back V = {dX/dx} v (x, state) */
void FEM_Ref_Vector (BODY *bod, ELEMENT *ele, double *x, double *v, double *V)
{
  MESH *msh = FEM_MESH (bod);
  double point [3], F [9];

  if (ele == NULL)
  {
    ele = MESH_Element_Containing_Point (msh, x, 0);
  }

  ASSERT_TEXT (ele, "Element containing referential point was not found.\n"
                    "Please report this bug!\n");

  spatial_to_local (msh, ele, x, point);
  deformation_gradient (bod, msh, ele, point, F);
  TVMUL (F, v, V);
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
    element_nodal_values (FEM_MESH_VEL0 (bod), ele, u);
    SET (vglo, 0); for (i = 0; i < n; i ++) { ADDMUL (vglo, shapes [i], u [i], vglo); } /* vglo = N u0 */
    TVMUL (base, vglo, prevel); /* prevel = base' vglo */
  }

  if (curvel)
  {
    element_nodal_values (FEM_MESH_VELO (bod), ele, u);
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
  short body_space = 0;
  DOM *dom = bod->dom;
  SOLFEC *sol = dom->solfec;
  double point [3];
  MX *N, *H;

  if (sol->kind == NEWTON_SOLVER)
  {
    NEWTON *ns = sol->solver;
    body_space = (ns->locdyn == LOCDYN_OFF);
  }

  TNCOPY (base, base_trans.x);
  referential_to_local (msh, ele, X, point);
  N = element_shapes_matrix (bod, msh, ele, point);
  H = MX_Matmat (1.0, &base_trans, N, 0.0, NULL);
  MX_Destroy (N);

  if ((bod->form == BODY_COROTATIONAL && bod->scheme != SCH_DEF_EXP && body_space == 0)  /* XXX => NEWTON_SOLVER sees the regular H = E' N ,
											    rather than H R when using the body-space mode,
											    since FEM_Invvec already incorportes R */
      || bod->form >= BODY_COROTATIONAL_MODAL) /* in this case H = E' N R * bod->evec, hence it is a projection of pulled-back E' N onto the reduced base */
  {
    double *x = H->x, *y = x + H->nzmax, *R = FEM_ROT (bod), T [9];

    ASSERT_DEBUG ((y-x) % 9 == 0, "Number of nonzeros in H not divisble by 9");

    for (; x < y; x += 9)
    {
      NNMUL (x, R, T);
      NNCOPY (T, x); /* H = E' N R <=> rotaions gets shifted to H */
    }

    if (bod->form >= BODY_COROTATIONAL_MODAL)
    {
      N = H; /* tentative */
      H = MX_Matmat (1.0, H, bod->evec, 0.0, NULL); /* H = E' N R * bod->evec */
      MX_Destroy (N);
    }
  }

  return H;
}

/* compute current kinetic energy */
double FEM_Kinetic_Energy (BODY *bod)
{
  if (bod->M)
  {
    switch (bod->form)
    {
    case BODY_COROTATIONAL_MODAL:
    case BODY_COROTATIONAL_REDUCED_ORDER:
    {
      ASSERT_DEBUG (bod->evec, "Reduced base must exist");

      double *u = bod->velo, *v, ene;
      int n = bod->dofs;

      ERRMEM (v = malloc (n * sizeof (double)));
      MX_Matvec (1.0, bod->M, u, 0.0, v);
      ene = 0.5 * blas_ddot (n, u, 1, v, 1);
      free (v);

      return ene;
    }
    break;
    default: /* TL, BC */
    {
      double *x = bod->M->x,
	     *y = x + bod->dofs,
	     *u = bod->velo,
	     sum;

      for (sum = 0.0; x < y; x ++, u ++) sum += (*u)*(*u) * (*x);

      return 0.5 * sum;
    }
    break;
    }
  }

  return 0.0;
}

/* get some values at a local point of an element */
void FEM_Element_Point_Values (BODY *bod, ELEMENT *ele, double *point, VALUE_KIND kind, double *values)
{
  MESH *msh = FEM_MESH (bod);

  switch (kind)
  {
  case VALUE_COORD:
  {
    MX *N = element_shapes_matrix (bod, msh, ele, point);
    MX_Matvec (1.0, N, FEM_MESH_CONF (bod), 1.0, values); /* XXX => (@@@) when called from within here, at input: values = X */
    MX_Destroy (N);                                      /* on the other hand, when called from rendering VALUE_COORD is not used */
  }
  break;
  case VALUE_DISPLACEMENT:
  {
    MX *N = element_shapes_matrix (bod, msh, ele, point);
    MX_Matvec (1.0, N, FEM_MESH_CONF (bod), 0.0, values);
    MX_Destroy (N);
  }
  break;
  case VALUE_DISP_NORM:
  {
    double disp[3];
    MX *N = element_shapes_matrix (bod, msh, ele, point);
    MX_Matvec (1.0, N, FEM_MESH_CONF (bod), 0.0, disp);
    MX_Destroy (N);
    values[0] = LEN(disp);
  }
  break;
  case VALUE_VELOCITY:
  {
    MX *N = element_shapes_matrix (bod, msh, ele, point);
    MX_Matvec (1.0, N, FEM_MESH_VELO (bod), 0.0, values);
    MX_Destroy (N);
  }
  break;
  case VALUE_VELO_NORM:
  {
    double velo[3];
    MX *N = element_shapes_matrix (bod, msh, ele, point);
    MX_Matvec (1.0, N, FEM_MESH_VELO (bod), 0.0, velo);
    MX_Destroy (N);
    values[0] = LEN(velo);
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

/* get some values at a referential point (ele can be NULL if not known) */
void FEM_Point_Values (BODY *bod, ELEMENT *ele, double *X, VALUE_KIND kind, double *values)
{
  MESH *msh = FEM_MESH (bod);
  double point [3];

  if (ele == NULL)
  {
    ele = MESH_Element_Containing_Point (msh, X, 1);
  }

  ASSERT_TEXT (ele, "point outside of referential mesh");

  int nod = ELEMENT_Ref_Point_To_Node (msh, ele, X);

  if (nod == -1)
  {
    ASSERT (ele, ERR_FEM_POINT_OUTSIDE);

    referential_to_local (msh, ele, X, point);

    if (kind == VALUE_COORD) { COPY (X, values); } /* see (@@@) */

    FEM_Element_Point_Values (bod, ele, point, kind, values);
  }
  else
  {
    FEM_Cur_Node_Values (bod, msh->cur_nodes[nod], kind, values);
  }
}

/* get some values at a curent mesh node (node points inside MESH->cur_nodes) */
void FEM_Cur_Node_Values (BODY *bod, double *node, VALUE_KIND kind, double *values)
{
  int i, j;
  double point [3];
  MESH *msh = FEM_MESH (bod);
  int n = (node_t) node - msh->cur_nodes;

  if (kind >= VALUE_STRESS) /* average from neigbouring elements */
  {
    ELEMENT *ele = MESH_Element_With_Node (msh, n);
    ASSERT (ele, ERR_MSH_ELEMENT_WITH_NODE);
    double X [3], v [7];
    SET *set, *item;
    int m;

    COPY (msh->ref_nodes [n], X);

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
  else
  {
    switch ((int)kind)
    {
    case VALUE_COORD:
    {
      double *conf = FEM_MESH_CONF (bod),
	     *q = &conf [3*n],
	     *p = msh->ref_nodes [n];
      ADD (p, q, values);
    }
    break;
    case VALUE_DISPLACEMENT:
    {
      double *conf = FEM_MESH_CONF (bod),
	     *q = &conf [3*n];
      COPY (q, values);
    }
    break;
    case VALUE_DISP_NORM:
    {
      double *conf = FEM_MESH_CONF (bod),
	     *q = &conf [3*n];
      values[0] = LEN(q);
    }
    break;
    case VALUE_VELOCITY:
    {
      double *velo = FEM_MESH_VELO (bod),
	     *u = &velo [3*n];
      COPY (u, values);
    }
    break;
    case VALUE_VELO_NORM:
    {
      double *velo = FEM_MESH_VELO (bod),
	     *u = &velo [3*n];
      values[0] = LEN(u);
    }
    break;
    }
  }
}

/* issued by state reading routines of body interface */
void FEM_Update_Rough_Mesh (BODY *bod)
{
  MESH *msh = bod->msh;
  double *q = FEM_MESH_CONF (bod),
	(*cur) [3] = msh->cur_nodes,
	(*ref) [3] = msh->ref_nodes,
	(*end) [3] = ref + msh->nodes_count;

  for (; ref < end; ref ++, cur ++, q += 3) { ADD (ref[0], q, cur[0]); }
}

/* split body by a referential plane; output one body with new boundary or two bodies if fragmentation occurs */
void FEM_Split (BODY *bod, double *point, double *normal, short topoadj, int surfid[2], BODY **one, BODY **two)
{
  SHAPE *copy, *sone, *stwo;
  MESH *mone, *mtwo;
  char *label;

  *one = *two = NULL;
  copy = SHAPE_Copy (bod->shape);
  SHAPE_Update (copy, NULL, NULL); /* restore reference configuration */
  SHAPE_Split (copy, point, normal, topoadj, surfid, &sone, &stwo); /* split in reference configuration */
  SHAPE_Destroy (copy);

  if (bod->msh)
  {
    MESH *copy = MESH_Copy (bod->msh);
    MESH_Update (copy, NULL, NULL, NULL);
    MESH_Split (copy, point, normal, topoadj, surfid, 1, &mone, &mtwo);
    MESH_Destroy (copy);
  }
  else
  {
    mone = NULL;
    mtwo = NULL;
  }

  if (bod->label) ERRMEM (label = malloc (strlen (bod->label) + 8));
  else label = NULL;

  if (sone)
  {
    ASSERT_DEBUG (!bod->msh || (bod->msh && mone), "Cut shape but not rought mesh");
    if (bod->label) sprintf (label, "%s/1", bod->label);
    (*one) = BODY_Create (bod->kind, sone, bod->mat, label, bod->flags & BODY_PERMANENT_FLAGS, bod->form, mone, NULL, NULL, NULL);
    FEM_Map_State (FEM_MESH (bod), FEM_MESH_CONF (bod), FEM_MESH_VELO (bod), FEM_MESH (*one), (*one)->conf, (*one)->velo);
    if (bod->form >= BODY_COROTATIONAL_MODAL)
    {
      /* TODO: compute and map reduced state of 'sone' */
      ASSERT (0, ERR_NOT_IMPLEMENTED); /* FIXME */
    }
    SHAPE_Update ((*one)->shape, (*one), (MOTION)BODY_Cur_Point); 
    if (mone) FEM_Update_Rough_Mesh (*one);
  }

  if (stwo)
  {
    ASSERT_DEBUG (!bod->msh || (bod->msh && mtwo), "Cut shape but not rought mesh");
    if (bod->label) sprintf (label, "%s/2", bod->label);
    (*two) = BODY_Create (bod->kind, stwo, bod->mat, label, bod->flags & BODY_PERMANENT_FLAGS, bod->form, mtwo, NULL, NULL, NULL);
    FEM_Map_State (FEM_MESH (bod), FEM_MESH_CONF (bod), FEM_MESH_VELO (bod), FEM_MESH (*two), (*two)->conf, (*two)->velo);
    if (bod->form >= BODY_COROTATIONAL_MODAL)
    {
      /* TODO: compute and map reduced state of 'stwo' */
      ASSERT (0, ERR_NOT_IMPLEMENTED); /* FIXME */
    }
    SHAPE_Update ((*two)->shape, (*two), (MOTION)BODY_Cur_Point); 
    if (mtwo) FEM_Update_Rough_Mesh (*two);
  }

  if (bod->label) free (label);
}

/* separate body whose shape is separable into sub-bodies */
BODY** FEM_Separate (BODY *bod, int *m)
{
  BODY **out = NULL;
  SHAPE **shp;
  MESH **msh;
  char *label;
  int i;

  shp = SHAPE_Separate (bod->shape, m);

  if (!shp) return NULL;

  ERRMEM (out = malloc ((*m) * sizeof (BODY*)));

  if (bod->msh)
  {
    msh = MESH_Separate (bod->msh, &i, 0);
    ASSERT_TEXT ((*m) == i, "Shape and background mesh fragmented into different numbers of fragments: %d, %d", *m, i);

    /* how do we know that the sequence of meshes maches the sequence of separated shapes ?
     * we need to do this mapping here => find meshes that contain inner convex points */

    for (i = 0; i < (*m); i ++) /* for each shape */
    {
      CONVEX *x = shp[i]->data;
      double *p = x->cur,
	     *e = p + 3 * x->nver,
	     q [3] = {0, 0, 0};

      for (; p < e; p += 3) { ACC (p, q); }
      DIV (q, (double) x->nver, q); /* centroid of the first convex */

      for (int j = i; j < (*m); j ++) /* for all unmatched meshes */
      {
	if (MESH_Element_Containing_Spatial_Point (msh [j], q)) /* if it contains the centroid */
	{
	  MESH *ptr = msh [i]; /* save this mesh pointer */
	  msh [i] = msh [j]; /* make ith mesh matching the ith shape */
	  msh [j] = ptr; /* put the saved pointer in the spare place */
	}
      }
    }
  }
  else msh = NULL;

  if (bod->label) ERRMEM (label = malloc (strlen (bod->label) + 8));
  else label = NULL;

  for (i = 0; i < (*m); i ++)
  {
    if (bod->label) sprintf (label, "%s/%d", bod->label, i);
    MESH *backmesh = (msh ? msh [i] : NULL);
    out [i] = BODY_Create (bod->kind, shp [i], bod->mat, label, bod->flags & BODY_PERMANENT_FLAGS, bod->form, backmesh, NULL, NULL, NULL);
    FEM_Map_State (FEM_MESH (bod), FEM_MESH_CONF (bod), FEM_MESH_VELO(bod), FEM_MESH (out [i]), out [i]->conf, out [i]->velo);
    if (bod->form >= BODY_COROTATIONAL_MODAL)
    {
      /* TODO: compute and map reduced state of 'out [i]' */
      ASSERT (0, ERR_NOT_IMPLEMENTED); /* FIXME */
    }
    SHAPE_Update (out [i]->shape, out [i], (MOTION)BODY_Cur_Point); 
    if (backmesh) FEM_Update_Rough_Mesh (out [i]);
  }

  if (bod->label) free (label);
  if (msh) free (msh);
  free (shp);

  return out;
}

/* release FEM memory */
void FEM_Destroy (BODY *bod)
{
  free (bod->conf);
  if (bod->field) free (bod->field);
}

/* get configuration write/read size */
int FEM_Conf_Size (BODY *bod)
{
  switch (bod->form)
  {
  case BODY_COROTATIONAL_MODAL:
  case BODY_COROTATIONAL_REDUCED_ORDER:
    return bod->dofs + 9; /* reduced conf, rotation */
  default: return bod->dofs; /* conf */
  }
}

#if MPI
/* get configuration packing size */
int FEM_Conf_Pack_Size (BODY *bod)
{
  switch (bod->form)
  {
  case BODY_COROTATIONAL_MODAL:
  case BODY_COROTATIONAL_REDUCED_ORDER:
    ASSERT_DEBUG (bod->evec, "Reduced base must be present at this point");
    return bod->dofs + 9 + bod->evec->m; /* reduced conf, R, full conf */
  case BODY_COROTATIONAL: return bod->dofs + 9; /* conf, R */
  default: return bod->dofs; /* conf */
  }
}

/* get velocity packing size */
int FEM_Velo_Pack_Size (BODY *bod)
{
  switch (bod->form)
  {
  case BODY_COROTATIONAL_MODAL:
  case BODY_COROTATIONAL_REDUCED_ORDER:
    ASSERT_DEBUG (bod->evec, "Reduced base must be present at this point");
    return 4 * bod->dofs + 2 * bod->evec->m; /* reduced: velo, vel0, fext, fint; full: velo, vel0 */
  default: return 4 * bod->dofs; /* velo, vel0, fext, fint */
  }
}
#endif

/* compute c = alpha * INVERSE (bod) * b + beta * c */
void FEM_Invvec (double alpha, BODY *bod, double *b, double beta, double *c)
{
  switch (bod->form)
  {
  case TOTAL_LAGRANGIAN:
  case BODY_COROTATIONAL_MODAL:
  case BODY_COROTATIONAL_REDUCED_ORDER:
    MX_Matvec (alpha, bod->inverse, b, beta, c);
    break;
  case BODY_COROTATIONAL:
    BC_matvec (alpha, bod->inverse, FEM_ROT (bod), b, beta, c);
    break;
  }
}

/* create approximate inverse operator */
MX* FEM_Approx_Inverse (BODY *bod)
{
  if (bod->form >= BODY_COROTATIONAL_MODAL || /* dense */
      bod->scheme == SCH_DEF_EXP) return MX_Copy (bod->inverse, NULL); /* dense or diagonal */
  else 
  {
    int *p, *i, n, k;
    double *x, *y, *z;
    MX *I;

    n = bod->dofs;
    ERRMEM (p = malloc (sizeof (int [n+1])));
    ERRMEM (i = malloc (sizeof (int [n])));
    ERRMEM (x = malloc (sizeof (double [n])));
    for (k = 0, p [n] = n; k < n; k ++) p [k] = i [k] = k; /* diagonal pattern */
    ERRMEM (I = malloc (sizeof (MX)));
    I->kind = MXCSC;
    I->flags = 0;
    I->nzmax = n;
    I->m = n;
    I->n = n;
    I->p = p;
    I->i = i;
    I->nz = -1;
    I->x = x;
    I->sym = I->num = NULL;

    for (y = x + n, z = bod->M->x; x != y; x ++, z ++) *x = 1.0 / (*z);

    return I;
  }
}

/* compute n lowest modal eigenvalues, given an absolute tolerance and iterations bound;
 * returns the corresponding eigenvectors in columns of a dense matrix, or NULL if the eigenvalue solver has failed */
MX* FEM_Modal_Analysis (BODY *bod, int n, double abstol, int maxiter, int verbose, double *val)
{
  MX *K = tangent_stiffness (bod, 1),
     *M = diagonal_inertia (bod, 1),
     *V = MX_Create (MXDENSE, K->n, ABS (n), NULL, NULL);
  int iters;

  iters = MX_CSC_Geneigen (K, M, n, abstol, maxiter, verbose, val, V);

  MX_Destroy (K);
  MX_Destroy (M);

  if (iters >= 0 && iters < maxiter) return V;
  else
  {
    MX_Destroy (V);
    return NULL;
  }
}

/* load an eigen mode as the current shape */
void FEM_Load_Mode (BODY *bod, int mode, double scale)
{
  if (bod->evec && mode >= 0 && mode < bod->evec->m)
  {
    double *vec0 = bod->evec->x + (bod->evec->m * mode), *vec, vmax;
    MESH *msh = FEM_MESH (bod);
    double (*ref) [3] = msh->ref_nodes, (*cur) [3] = msh->cur_nodes;
    int i, n = msh->nodes_count;

    for (i = 0, vec = vec0, vmax = 0.0; i < n; i ++, vec += 3)
    {
      if (vec [0] > vmax) vmax = vec [0];
      if (vec [1] > vmax) vmax = vec [1];
      if (vec [2] > vmax) vmax = vec [2];
    }

    if (vmax != 0.0) scale *= 1.0 / vmax; /* max displacement == scale */

    for (i = 0, vec = vec0; i < n; i ++, ref ++, cur ++, vec += 3)
    {
      ADDMUL (ref [0], scale, vec, cur [0]);
    }

    SHAPE_Update (bod->shape, bod, (MOTION)modal_motion); 
  }
}

/* export M and K in MatrixMarket formats; in 'spd' mode only lower tirangle is used */
void FEM_MatrixMarket_M_K (BODY *bod, short spdM, char *pathM, short spdK, char *pathK)
{
  MX *K = tangent_stiffness (bod, spdK),
     *M = diagonal_inertia (bod, spdM);

  MX_MatrixMarket (K, pathK);
  MX_MatrixMarket (M, pathM);

  MX_Destroy (K);
  MX_Destroy (M);
}

/* called after reading to post-process internal data */
void FEM_Post_Read (BODY *bod)
{
  if (bod->form >= BODY_COROTATIONAL_MODAL)
  {
    RO_lift_conf (bod, bod->evec, FEM_MESH (bod), FEM_ROT (bod), bod->conf, FEM_MESH_CONF (bod));
    MX_Matvec (1.0, bod->evec, bod->velo, 0.0, FEM_MESH_VELO (bod));
    RO_rotate_forward (FEM_ROT (bod), FEM_MESH_VELO (bod), bod->evec->m);
  }
}

/* return mesh dofs count */
int FEM_Mesh_Dofs (BODY *bod)
{
  return MESH_DOFS(FEM_MESH(bod));
}

/* return mesh displacements */
double* FEM_Mesh_Conf (BODY *bod)
{
  return FEM_MESH_CONF (bod);
}

/* output mesh co-rotated displacements */
void FEM_Mesh_Corotated_Conf (BODY *bod, double *disp)
{
  double *qm = FEM_MESH_CONF(bod), *R, TMP[9], (*Z) [3], Y [3], *x, *y;
  MESH *msh = FEM_MESH(bod);
  int qmsize = MESH_DOFS(msh);

  if (bod->form >= BODY_COROTATIONAL_MODAL)
  {
    R = FEM_ROT(bod);
  }
  else
  {
    IDENTITY (TMP); /* FIXME --> in time-sequencial use previous rotation is preffered */
    R = TMP;
    BC_update_rotation (bod, msh, qm, R);
  }

  blas_dcopy (qmsize, qm, 1, disp, 1);

  for (x = disp, y = disp+qmsize, Z = msh->ref_nodes; x < y; x += 3, Z ++)
  {
    NVMUL (R, Z[0], Y);
    SUB (Z[0], Y, Y);
    ACC (Y, x); /* d = (I-R)Z + qm */

    COPY (x, Y);
    TVMUL (R, Y, x); /* R' d */
  }
}

/* output six rigid body displacements */
void FEM_Mesh_Rigid_Displacements (BODY *bod, double *disp)
{
  double (*Z) [3], *C, Y[3], *dx, *dy, *dz, *rx, *ry, *rz, inv;
  MESH *msh = FEM_MESH(bod);
  int i, n = msh->nodes_count, dofs = 3*n;

  Z = msh->ref_nodes;
  C = bod->ref_center;
  dx = disp;
  dy = dx+dofs;
  dz = dy+dofs;
  rx = dz+dofs;
  ry = rx+dofs;
  rz = ry+dofs;

  for (i = 0; i < n; i ++, dx += 3, dy += 3, dz += 3, rx += 3, ry += 3, rz += 3, Z ++)
  {
    dx[0] = 1.0;
    dx[1] = 0.0;
    dx[2] = 0.0;

    dy[0] = 0.0;
    dy[1] = 1.0;
    dy[2] = 0.0;

    dz[0] = 0.0;
    dz[1] = 0.0;
    dz[2] = 1.0;

    SUB (Z[0], C, Y);

    PRODUCT (dx, Y, rx);
    PRODUCT (dy, Y, ry);
    PRODUCT (dz, Y, rz);
  }

  dx = disp;
  dy = dx+dofs;
  dz = dy+dofs;
  rx = dz+dofs;
  ry = rx+dofs;
  rz = ry+dofs;

  inv = blas_ddot (dofs, dx, 1, dx, 1);
  inv = 1.0/sqrt(inv);
  blas_dscal (dofs, inv, dx, 1);

  inv = blas_ddot (dofs, dy, 1, dy, 1);
  inv = 1.0/sqrt(inv);
  blas_dscal (dofs, inv, dy, 1);

  inv = blas_ddot (dofs, dz, 1, dz, 1);
  inv = 1.0/sqrt(inv);
  blas_dscal (dofs, inv, dz, 1);

  inv = blas_ddot (dofs, rx, 1, rx, 1);
  inv = 1.0/sqrt(inv);
  blas_dscal (dofs, inv, rx, 1);

  inv = blas_ddot (dofs, ry, 1, ry, 1);
  inv = 1.0/sqrt(inv);
  blas_dscal (dofs, inv, ry, 1);

  inv = blas_ddot (dofs, rz, 1, rz, 1);
  inv = 1.0/sqrt(inv);
  blas_dscal (dofs, inv, rz, 1);
}

/* count mesh surface integration points;
 * surface == INT_MAX --> entire surface */
int FEM_Mesh_Surface_Integration_Points_Count (MESH *msh, int surface)
{
  FACE *fac;
  int n = 0;

  for (fac = msh->faces; fac; fac = fac->n)
  {
    if (surface == INT_MAX || fac->surface == surface)
    {
      INTEGRAL2D_BEGIN (fac->type) { n ++; } INTEGRAL2D_END ()
    }
  }

  return n;
}

/* get mesh surface integration points (referential coordinates);
 * surface == INT_MAX --> entire surface */
void FEM_Mesh_Surface_Integration_Points_Get (MESH *msh, int surface, double *refpnt)
{
  double refn [4][3], shapes [4];
  FACE *fac;

  for (fac = msh->faces; fac; fac = fac->n)
  {
    if (surface == INT_MAX || fac->surface == surface)
    {
      face_nodes (msh->ref_nodes, fac->type, fac->nodes, refn);

      INTEGRAL2D_BEGIN (fac->type) /* defines point and weight */
      {
	face_shapes (fac, point, shapes);

	refpnt[0] = shapes[0]*refn[0][0]+shapes[1]*refn[1][0]+shapes[2]*refn[2][0];
	refpnt[1] = shapes[0]*refn[0][1]+shapes[1]*refn[1][1]+shapes[2]*refn[2][1];
	refpnt[2] = shapes[0]*refn[0][2]+shapes[1]*refn[1][2]+shapes[2]*refn[2][2];

	if (fac->type == 4)
	{
	  refpnt[0] += shapes[3]*refn[3][0];
	  refpnt[1] += shapes[3]*refn[3][1];
	  refpnt[2] += shapes[3]*refn[3][2];
	}

	refpnt += 3;
      }
      INTEGRAL2D_END ()
    }
  }
}

/* count mesh volume integration points;
 * volume == INT_MAX --> entire volume */
int FEM_Mesh_Volume_Integration_Points_Count (MESH *msh, int volume)
{
  ELEMENT *ele;
  int n = 0;

  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    if (volume == INT_MAX || ele->volume == volume)
    {
      n += number_of_integration_points (ele);
    }
  }

  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    if (volume == INT_MAX || ele->volume == volume)
    {
      n += number_of_integration_points (ele);
    }
  }
 
  return n;
}

/* get mesh volume integration points (referential coordinates);
 * volume == INT_MAX --> entire volume */
void FEM_Mesh_Volume_Integration_Points_Get (MESH *msh, int volume, double *refpnt)
{
  double refn [MAX_NODES][3], shapes [MAX_NODES];
  ELEMENT *ele;
  int bulk, i;

  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    if (volume == INT_MAX || ele->volume == volume)
    {
      element_nodes (msh->ref_nodes, ele->type, ele->nodes, refn);

      INTEGRATE3D (ele->type, INTF, ele->dom, ele->domnum,

        element_shapes (ele->type, point, shapes);

	refpnt[0] = shapes[0]*refn[0][0]+shapes[1]*refn[1][0]+shapes[2]*refn[2][0]+shapes[3]*refn[3][0];
	refpnt[1] = shapes[0]*refn[0][1]+shapes[1]*refn[1][1]+shapes[2]*refn[2][1]+shapes[3]*refn[3][1];
	refpnt[2] = shapes[0]*refn[0][2]+shapes[1]*refn[1][2]+shapes[2]*refn[2][2]+shapes[3]*refn[3][2];

	for (i = 4; i < ele->type; i ++)
	{
	  refpnt[0] += shapes[i]*refn[i][0];
	  refpnt[1] += shapes[i]*refn[i][1];
	  refpnt[2] += shapes[i]*refn[i][2];
	}

	refpnt += 3;
      )
    }

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }
}
