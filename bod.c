/*
 * bod.c
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * general body
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
#include "sol.h"
#include "lap.h"
#include "bla.h"
#include "mem.h"
#include "alg.h"
#include "err.h"
#include "bod.h"
#include "svk.h"
#include "dom.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "eli.h"
#include "pck.h"
#include "err.h"
#include "lng.h"
#include "fem.h"
#include "but.h"
#include "rnd.h"
#include "put.h"

/* implicit PRB integration */
#define IMP_EPS 1E-8
#define MAX_ITERS 16

/* energy conservation tolerance */
#define ENE_TOL 0.01

/* sizes */
#define RIG_CONF_SIZE	15	/* rotation matrix, mass center position, auxiliary vector of size 3 */
#define RIG_VELO_SIZE   18      /* referential angular velocity, spatial linear velocity, copies of both at (t-h), fext */

#define PRB_CONF_SIZE	12	/* deformation gradient (row-wise), mass cener position */
#define PRB_VELO_SIZE	45	/* deformation velocity, spatial linear velocity, copies of both at (t-h), fext, fint */

/* parameters */
#define RIG_SOLVER_EPSILON (DBL_EPSILON * 1000.0)
#define RIG_SOLVER_MAXITER 64

/* member access macros */
#define BOD_M0(bod) (bod)->ref_mass
#define BOD_V0(bod) (bod)->ref_volume
#define BOD_X0(bod) (bod)->ref_center
#define BOD_E0(bod) (bod)->ref_tensor
#define BOD_J0(bod) (bod)->ref_tensor
#define RIG_ROTATION(bod) (bod)->conf
#define RIG_CENTER(bod) ((bod)->conf+9)
#define RIG_CENTER1(conf) ((conf)+9)
#define RIG_AUXILIARY(bod) ((bod)->conf+12)
#define RIG_ANGVEL(bod) (bod)->velo
#define RIG_LINVEL(bod) ((bod)->velo+3)
#define RIG_ANGVEL0(bod) ((bod)->velo+6)
#define RIG_LINVEL0(bod) ((bod)->velo+9)
#define RIG_FEXT(bod) ((bod)->velo+12)
#define PRB_GRADIENT(bod) (bod)->conf
#define PRB_CENTER(bod) ((bod)->conf+9)
#define PRB_CENTER1(conf) ((conf)+9)
#define PRB_GRADVEL(bod) (bod)->velo
#define PRB_LINVEL(bod) ((bod)->velo+9)
#define PRB_GRADVEL0(bod) ((bod)->velo+12)
#define PRB_LINVEL0(bod) ((bod)->velo+21)
#define PRB_FEXT(bod) ((bod)->velo+24)
#define PRB_FINT(bod) ((bod)->velo+36)

/* allocate memory for a body */
static void* alloc_body (short kind)
{
  switch (kind)
  {
    case RIG:
    /* alloc rigid body structure  */
    return MEM_CALLOC (sizeof (BODY) + sizeof (double [RIG_CONF_SIZE + RIG_VELO_SIZE]));

    case PRB:
    /* alloc pseudo-rigid body structure */
    return MEM_CALLOC (sizeof (BODY) + sizeof (double [PRB_CONF_SIZE + PRB_VELO_SIZE]));

    case FEM:
    /* alloc finite element discretised body =>
     * configuration and velocity are allocated individually */
    return MEM_CALLOC (sizeof (BODY));
  }

  return NULL;
}

/* find resultant pressure vector */
static void resultant_pressure (MESH *msh, FORCE *frc)
{
  int surfid = frc->surfid;
  double *point = frc->ref_point,
	 *dir = frc->direction,
	 (*cur) [3] = msh->cur_nodes,
	 a0, a1, dv;
  FACE *fac;

  SET (point, 0.0);
  SET (dir, 0.0);
  dv = 0.0;

  for (fac = msh->faces; fac; fac = fac->n)
  {
    if (fac->surface == surfid)
    {
      double *a = cur [fac->nodes [0]],
	     *b = cur [fac->nodes [1]],
	     *c = cur [fac->nodes [2]],
	     *d = fac->type == 4 ? cur [fac->nodes [3]] : NULL;

      TRIANGLE_AREA (a, b, c, a0);
      if (fac->type == 4)
      {
        TRIANGLE_AREA (c, d, a, a1);
        a0 += a1;
	point [0] += 0.25 * (a[0]+b[0]+c[0]+d[0]);
	point [1] += 0.25 * (a[1]+b[1]+c[1]+d[1]);
	point [2] += 0.25 * (a[2]+b[2]+c[2]+d[2]);
      }
      else
      {
	point [0] += (a[0]+b[0]+c[0])/3.0;
	point [1] += (a[1]+b[1]+c[1])/3.0;
	point [2] += (a[2]+b[2]+c[2])/3.0;
      }

      ADDMUL (dir, a0, fac->normal, dir); /* resultant direction = normal * area */
      dv += 1.0;
    }
  }

  if (dv > 0.0)
  {
    DIV (point, dv, point); /* average point */
  }
}

/* ------------------- RIG --------------------- */

/* implicit solver of => exp[hW]JW = G, outputing W and A = exp[hW] */
static double SOLVE (double h, double *J, double *W, double *G, double *A)
{
  int ipiv [3], i = 0;
  double O [3],
         Z0 [9],
	 Z1 [9],
	 Z2 [9],
	 JW [3],
	 Z [9],
	 B [9],
	 R [3],
	 error,
	 level;

  COPY   (W, O);
  SCALE  (O, h);
  EXPMAP (O, A);
  NNMUL  (A, J, B);
  NVMUL  (B, W, R);
  SUB    (R, G, R) /* initial residual => R = exp [hW]JW - G */

  level = 1E-12 * h;
  MAXABS (R, error);
  while (error > level)
  {
    NVMUL     (J, W, JW);
    EXPMAP123 (O, Z0, Z1, Z2);
    NVMUL     (Z0, JW, Z);
    NVMUL     (Z1, JW, Z+3);
    NVMUL     (Z2, JW, Z+6);
    SCALE9    (Z, h);
    NNADD     (Z, B, Z); /* tangent => dexp[hW]*h*JW + exp[hW]*J */

    ASSERT (lapack_dgesv (3, 1, Z, 3, ipiv, R, 3) == 0, ERR_BOD_NEW3_SINGULAR_JACOBIAN);

    SUB (W, R, W); 

    COPY   (W, O);
    SCALE  (O, h);
    EXPMAP (O, A);
    NNMUL  (A, J, B);
    NVMUL  (B, W, R);
    SUB    (R, G, R) /* residual => R = exp [hW]JW - G */

    MAXABS (R, error);
    ASSERT (++i < 50, ERR_BOD_NEW3_NEWTON_DIVERGENCE);
  }

  return LEN (O);
}

/* calculates transformation operator from the generalized velocity
 * space to the local Cartesian space at 'X' in given the 'base' */
static void rig_operator_H (BODY *bod, double *X, double *base, double *H)
{
  double A [3], S [9], T [9],
	 *R = RIG_ROTATION (bod),
	 *X0 = BOD_X0 (bod);
  

  SUB (X, X0, A);
  VECSKEW (A, S);
  NTMUL (R, S, T);
  TNMUL (base, T, H);
  TNCOPY (base, H + 9);
}

/* set up the inverse of the inertia
 * for the dynamic time stepping */
static void rig_dynamic_inverse (BODY *bod)
{
  double *J, *x, m;
  MX *M, *I;

  if (!bod->inverse)
  {
    int p [] = {0, 9, 18},
	i [] = {0, 3, 6};

    M = MX_Create (MXBD, 6, 2, p, i);
    bod->M = M;

    I = MX_Create (MXBD, 6, 2, p, i);
    bod->inverse = I;
  }
  else M = bod->M, I = bod->inverse;

  J = BOD_J0 (bod);
  m = BOD_M0 (bod);
  x = M->x;
  NNCOPY (J, x);
  x += 9;
  IDENTITY (x);
  SCALEDIAG (x, m);
  MX_Inverse (M, I);
}
 
/* no difference for the static inverse routine */
#define rig_static_inverse(bod) rig_dynamic_inverse(bod)

/* calculate linear resultant, spatial torque, referential torque at time t */
static void rig_force (BODY *bod, double *q, double *u, double t, double h,
                       double *linforc, double *spatorq, double *reftorq)
{
  double *X0 = BOD_X0 (bod),
	 A [3], a [3], f [3], b [3], v;
  short kind;

  SET (linforc, 0);
  SET (spatorq, 0);
  SET (reftorq, 0);

  for (FORCE *frc = bod->forces; frc; frc = frc->next)
  {
    if (frc->func)
    {
      int nq = BODY_Conf_Size (bod),
	  nu = bod->dofs;
      double f [9] = {0,0,0,0,0,0,0,0,0};
      frc->func (frc->data, frc->call, nq, q, nu, u, t, h, f);
      ADD (f, linforc, linforc);
      ADD (f+3, spatorq, spatorq);
      ADD (f+6, reftorq, reftorq);
    }
    else
    {
      if (frc->kind & PRESSURE)
      {
	resultant_pressure (bod->shape->data, frc);
	kind = 0;
      }
      else kind = frc->kind;
      v = TMS_Value (frc->data, t);
      COPY (frc->direction, f);
      SCALE (f, v);
      
      if (kind & CONVECTED)
      { 
	if (kind & TORQUE)
	{ 
	  ADD (reftorq, f, reftorq);
	}
	else
	{
	  NVMUL (q, f, b);             /* b => spatial force */
	  SUB (frc->ref_point, X0, A); /* (X - X0) => referential force arm */
	  NVMUL (q, A, a);             /* a => spatial force arm */
	  PRODUCTADD (a, b, spatorq);  /* (x-x0) x b => spatial torque */
	  ADD (linforc, b, linforc);
	}
      }
      else
      {
	if (kind & TORQUE)
	{ 
	  ADD (spatorq, f, spatorq);
	}
	else
	{
	  SUB (frc->ref_point, X0, A); /* (X - X0) => referential force arm */
	  NVMUL (q, A, a);             /* a => spatial force arm */
	  PRODUCTADD (a, f, spatorq);  /* (x-x0) x b => spatial torque */
	  ADD (linforc, f, linforc);
	}
      }
    }
  }

  if (bod->dom->gravity [0])
  {
    f [0] = TMS_Value (bod->dom->gravity [0], t);
    f [1] = TMS_Value (bod->dom->gravity [1], t);
    f [2] = TMS_Value (bod->dom->gravity [2], t);

    ADDMUL (linforc, BOD_M0 (bod), f, linforc);
  }
}

/* calculate force (time + step) = external (time + step) - internal (time + step),
 * provided that the current state of the body is given as R(time), W(time); */
static void rig_static_force (BODY *bod, double time, double step, double *force)
{
  double *J  = BOD_J0(bod),
	 *I  = bod->inverse->x,
	 *W  = RIG_ANGVEL(bod),
	 *R  = RIG_ROTATION(bod),
	 *v  = RIG_LINVEL(bod),
	 *x  = RIG_CENTER(bod),
	 O [9], DR [9], W1 [3], q [12],
	 spatorq [3], reftorq [3];
 
  COPY (W, O);
  SCALE (O, step);
  EXPMAP (O, DR); 
  NNCOPY (R, O);
  NNMUL (R, DR, q);          /* R(t+h) = R(t) exp [h * W(t)] */
  ADDMUL (x, step, v, q+9);  /* x(t+h) = x(t) + h * v(t) */

  rig_force (bod, q, bod->velo,
    time+step, step, force + 3, spatorq, reftorq);
  TVADDMUL (reftorq, R, spatorq, force);

  COPY (force, O);
  SCALE (O, step); 
  NVMUL (J, W, W1);
  TVADDMUL (O, DR, W1, O);
  NVMUL (I, O, W1);            /* W1 = I * [exp [- h * W(t)]*J*W(t) + h*T(t+h)] */

  NVMUL (J, W1, O);
  PRODUCTSUB (W1, O, force);   /* force [0..2] = T - W1 x J * W1 */
}

/* accumulate constraints reaction */
inline static void rig_constraints_force_accum (BODY *bod, double *point, double *base, double *R, short isma, double *force)
{
  double H [18], r [6];

  rig_operator_H (bod, point, base, H);
  blas_dgemv ('T', 3, 6, 1.0, H, 3, R, 1, 0.0, r, 1);

  if (isma)
  {
    force [0] += r[0];
    force [1] += r[1];
    force [2] += r[2];
    force [3] += r[3];
    force [4] += r[4];
    force [5] += r[5];
  }
  else
  {
    force [0] -= r[0];
    force [1] -= r[1];
    force [2] -= r[2];
    force [3] -= r[3];
    force [4] -= r[4];
    force [5] -= r[5];
  }
}

/* calculate constraints reaction */
static void rig_constraints_force (BODY *bod, double *force)
{
  SET *node;

  force [0] = force [1] = force [2] =
  force [3] = force [4] = force [5] = 0.0;

  for (node = SET_First (bod->con); node; node = SET_Next (node))
  {
    CON *con = node->data;
    short isma = (bod == con->master);
    double *point = (isma ? con->mpnt : con->spnt);

    rig_constraints_force_accum (bod, point, con->base, con->R, isma, force);
  }
}

/* ------------------- PRB --------------------- */

/* calculates transformation operator from generalized velocity
 * space to the local cartesian space at 'X' in given 'base' */
static void prb_operator_H (BODY *bod, double *X, double *base, double *H)
{
  double A [3];

  memset (H, 0, sizeof (double [36]));

  SUB (X, bod->ref_center, A); 

  H [0] = base [0] * A [0];
  H [3] = base [0] * A [1];
  H [6] = base [0] * A [2];
  H [9] = base [1] * A [0];
  H [12] = base [1] * A [1];
  H [15] = base [1] * A [2];
  H [18] = base [2] * A [0];
  H [21] = base [2] * A [1];
  H [24] = base [2] * A [2];
  H [27] = base [0];
  H [30] = base [1];
  H [33] = base [2];

  H [1] = base [3] * A [0];
  H [4] = base [3] * A [1];
  H [7] = base [3] * A [2];
  H [10] = base [4] * A [0];
  H [13] = base [4] * A [1];
  H [16] = base [4] * A [2];
  H [19] = base [5] * A [0];
  H [22] = base [5] * A [1];
  H [25] = base [5] * A [2];
  H [28] = base [3];
  H [31] = base [4];
  H [34] = base [5];

  H [2] = base [6] * A [0];
  H [5] = base [6] * A [1];
  H [8] = base [6] * A [2];
  H [11] = base [7] * A [0];
  H [14] = base [7] * A [1];
  H [17] = base [7] * A [2];
  H [20] = base [8] * A [0];
  H [23] = base [8] * A [1];
  H [26] = base [8] * A [2];
  H [29] = base [6];
  H [32] = base [7];
  H [35] = base [8];
}

/* set up inverse of inertia for
 * the explicit dynamic time stepping */
static void prb_dynamic_explicit_inverse (BODY *bod)
{
  double *E0, m, *x;
  MX *M, *I;

  if (!bod->inverse)
  {
    int p [] = {0, 9, 18, 27, 36},
        i [] = {0, 3, 6, 9, 12};

    M = MX_Create (MXBD, 12, 4, p, i);
    bod->M = M;

    I = MX_Create (MXBD, 12, 4, p, i);
    bod->inverse = I;
  }
  else M = bod->M, I = bod->inverse;

  E0 = bod->ref_tensor; /* referential Euler tensor */
  m = bod->ref_mass;
  x = M->x;
  NNCOPY (E0, x); x += 9;
  NNCOPY (E0, x); x += 9;
  NNCOPY (E0, x); x += 9;
  IDENTITY (x);
  SCALEDIAG (x, m);
  MX_Inverse (M, I);
}

/* set up inverse operator for the implicit dynamic time stepping */
static void prb_dynamic_implicit_inverse (BODY *bod, double step, double *conf, double *force)
{
  MX *M, *A, *K;

  /* place-holder for the
   * tangent operator */
  if (!bod->inverse)
  {
    int pm [] = {0, 9, 18, 27, 36},
	im [] = {0, 3, 6, 9, 12},
        pa [] = {0, 81, 90},
        ia [] = {0, 9, 12};
    double *E0, m, *x;

    M = MX_Create (MXBD, 12, 4, pm, im);
    bod->M = M;

    K = MX_Create (MXDENSE, 9, 9, NULL, NULL);
    bod->K = K;

    /* set up mass matrix */
    E0 = bod->ref_tensor;
    m = bod->ref_mass;
    x = M->x;
    NNCOPY (E0, x); x += 9;
    NNCOPY (E0, x); x += 9;
    NNCOPY (E0, x); x += 9;
    IDENTITY (x);
    SCALEDIAG (x, m);

    A = MX_Create (MXBD, 12, 2, pa, ia);
    bod->inverse = A;
  }
  else M = bod->M, K = bod->K, A = bod->inverse;

  /* calculate stiffness matrix */
  SVK_Tangent_R (lambda (bod->mat->young, bod->mat->poisson),
    mi (bod->mat->young, bod->mat->poisson),
    bod->ref_volume, 9, conf ? conf : bod->conf, K->x);

  if (force)
  {
    /* account for the previous velocity */
    MX_Matvec (1.0 / step, M, bod->velo, 1.0, force);

    /* account for the internal force increment */
    MX_Matvec (-0.25 * step, K, bod->velo, 1.0, force);
  }

  /* calculate tangent operator A = M + (damping*h + h*h/4) K */
  MX_Add (1.0, MX_Diag(M, 0, 2), bod->damping*step + 0.25*step*step, K, MX_Diag(A, 0, 0));
  MX_Copy (MX_Diag (M, 3, 3), MX_Diag (A, 1, 1));

  /* invert A */
  MX_Inverse (A, A);
}

/* update tangent stiffness */
static void prb_tangent_stiffness (BODY *bod)
{
  if (!bod->K) bod->K = MX_Create (MXDENSE, 9, 9, NULL, NULL);

  SVK_Tangent_R (lambda (bod->mat->young, bod->mat->poisson),
    mi (bod->mat->young, bod->mat->poisson),
    bod->ref_volume, 9, bod->conf, bod->K->x);
}

/* set up inverse of inertia for
 * the implicit static time stepping */
static void prb_static_inverse (BODY *bod, double step)
{
  int p [] = {0, 9, 18, 27, 36},
      i [] = {0, 3, 6, 9, 12};
  double *E0, m, *x;
  MX_BD (M, 36, 12, 4, p, i);
  MX_DENSE (K, 9, 9);
  double eigmax;
  MX *IM, *IMK, *A;

  /* place-holder for the
   * tangent operator */
  if (!bod->inverse)
  {
    int p [] = {0, 81, 90},
        i [] = {0, 9, 12};

    A = MX_Create (MXBD, 12, 2, p, i);
    bod->inverse = A;
  }
  else A = bod->inverse;

  /* set up mass matrix */
  E0 = bod->ref_tensor;
  m = bod->ref_mass;
  x = M.x;
  NNCOPY (E0, x); x += 9;
  NNCOPY (E0, x); x += 9;
  NNCOPY (E0, x); x += 9;
  IDENTITY (x);
  SCALEDIAG (x, m);

  /* calculate stiffness matrix */
  SVK_Tangent_R (lambda (bod->mat->young, bod->mat->poisson),
    mi (bod->mat->young, bod->mat->poisson),
    bod->ref_volume, 9, bod->conf, K.x);

  /* inverse "deformable" sub-block of mass matrix */
  IM = MX_Inverse (MX_Diag (&M, 0, 2), NULL);

  /* get maximal eigenvalue of inv (M) * K */
  IMK = MX_Matmat (1.0, IM, &K, 0.0, NULL);
  MX_Eigen (IMK, 1, &eigmax, NULL); /* compute maximal eigenvalue */
  ASSERT (eigmax > 0.0, ERR_BOD_MAX_FREQ_LE0);

  /* calculate scaled tangent operator */
  MX_Add (eigmax / 4.0, MX_Diag(&M, 0, 2), step*step, &K, MX_Diag(A, 0, 0));
  MX_Copy (MX_Diag (&M, 3, 3), MX_Diag (A, 1, 1));

  /* invert A */
  MX_Inverse (A, A);

  /* clean up */
  MX_Destroy (IMK);
  MX_Destroy (IM);
}

/* compute internal force for the pseudo-rigid model */
inline static void prb_internal_force (BODY *bod, double *F, double *P)
{
  SVK_Stress_R (lambda (bod->mat->young, bod->mat->poisson),
                mi (bod->mat->young, bod->mat->poisson),
		bod->ref_volume, F, P);
}
 
/* compute force(time) = external(time) - internal(time), provided
 * that the configuration(time) has already been set elsewhere */
static void prb_dynamic_force (BODY *bod, double time, double step, double *fext, double *fint, double *force)
{
  double *F = PRB_GRADIENT (bod),
	 *X0 = BOD_X0 (bod),
	 P [9], A [3], f [12],
	 value;

  SET12 (force, 0);

  for (FORCE *frc = bod->forces; frc; frc = frc->next)
  {
    if (frc->func)
    {
      int nq = BODY_Conf_Size (bod);
      frc->func (frc->data, frc->call, nq, bod->conf, bod->dofs, bod->velo, time, step, f);
      blas_daxpy (12, 1.0, f, 1, force, 1);
    }
    else
    {
      if (frc->kind & PRESSURE) resultant_pressure (bod->shape->data, frc);
      value = TMS_Value (frc->data, time);
      COPY (frc->direction, f);
      SUB (frc->ref_point, X0, A);
      SCALE (f, value);

      if (frc->kind & CONVECTED) /* obtain spatial force */
      { 
	NVMUL (F, f, P);
	COPY (P, f);
      }

      force [0] += A [0]*f [0];
      force [1] += A [1]*f [0];
      force [2] += A [2]*f [0];
      force [3] += A [0]*f [1];
      force [4] += A [1]*f [1];
      force [5] += A [2]*f [1];
      force [6] += A [0]*f [2];
      force [7] += A [1]*f [2];
      force [8] += A [2]*f [2];
      force [9] += f [0];
      force [10] += f [1];
      force [11] += f [2];
    }
  }

  if (bod->dom->gravity [0])
  {
    f [0] = TMS_Value (bod->dom->gravity [0], time);
    f [1] = TMS_Value (bod->dom->gravity [1], time);
    f [2] = TMS_Value (bod->dom->gravity [2], time);

    ADDMUL (force+9, BOD_M0 (bod), f, force+9);
  }

  COPY12 (force, fext);

  SVK_Stress_R (lambda (bod->mat->young, bod->mat->poisson),
    mi (bod->mat->young, bod->mat->poisson), bod->ref_volume, F, P); /* internal force */

  NNCOPY (P, fint);

  NNSUB (force, P, force); /* force = external - internal */
}

/* the smame computation for the static case */
#define prb_static_force(bod, time, step, fext, fint, force) prb_dynamic_force (bod,time,step,fext,fint,force)

/* accumulate constraints reaction */
inline static void prb_constraints_force_accum (BODY *bod, double *point, double *base, double *R, short isma, double *force)
{
  double H [36], r [12];

  prb_operator_H (bod, point, base, H);
  blas_dgemv ('T', 3, 12, 1.0, H, 3, R, 1, 0.0, r, 1);

  if (isma)
  {
    force [0]  += r[0];
    force [1]  += r[1];
    force [2]  += r[2];
    force [3]  += r[3];
    force [4]  += r[4];
    force [5]  += r[5];
    force [6]  += r[6];
    force [7]  += r[7];
    force [8]  += r[8];
    force [9]  += r[9];
    force [10] += r[10];
    force [11] += r[11];
  }
  else
  {
    force [0]  -= r[0];
    force [1]  -= r[1];
    force [2]  -= r[2];
    force [3]  -= r[3];
    force [4]  -= r[4];
    force [5]  -= r[5];
    force [6]  -= r[6];
    force [7]  -= r[7];
    force [8]  -= r[8];
    force [9]  -= r[9];
    force [10] -= r[10];
    force [11] -= r[11];
  }
}

/* calculate constraints reaction */
static void prb_constraints_force (BODY *bod, double *force)
{
  SET *node;

  force [0] = force [1]  = force [2]  =
  force [3] = force [4]  = force [5]  =
  force [6] = force [7]  = force [8]  =
  force [9] = force [10] = force [11] = 0.0;

  for (node = SET_First (bod->con); node; node = SET_Next (node))
  {
    CON *con = node->data;
    short isma = (bod == con->master);
    double *point = (isma ? con->mpnt : con->spnt);

    prb_constraints_force_accum (bod, point, con->base, con->R, isma, force);
  }
}

/* calculate cauche stress for a pseudo-rigid body */
void prb_cauchy (BODY *bod, double *stress)
{
  double P [9], J, *F;

  F = bod->conf;
  J = SVK_Stress_R (lambda (bod->mat->young, bod->mat->poisson),
    mi (bod->mat->young, bod->mat->poisson), 1.0, F, P); /* per unit volume */

  stress [0] = (F [0]*P [0] + F [1]*P [3] + F [2]*P [6]) / J;
  stress [1] = (F [3]*P [1] + F [4]*P [4] + F [5]*P [7]) / J;
  stress [2] = (F [6]*P [2] + F [7]*P [5] + F [8]*P [8]) / J;
  stress [3] = (F [0]*P [1] + F [1]*P [4] + F [2]*P [7]) / J;
  stress [4] = (F [0]*P [2] + F [1]*P [5] + F [2]*P [8]) / J;
  stress [5] = (F [3]*P [2] + F [4]*P [5] + F [5]*P [8]) / J;
}

/* -------------------------------------- */

/* copy user label */
static char* copylabel (char *label)
{
  char *out = NULL;
  int l;

  if ((l = label ? strlen (label) : 0))
  {
    ERRMEM (out = malloc (l + 1));
    strcpy (out, label);
  }

  return out;
}

#if MPI
static int conf_pack_size (BODY *bod)
{
  switch (bod->kind)
  {
  case OBS: return 0;
  case RIG: return RIG_CONF_SIZE;
  case PRB: return PRB_CONF_SIZE;
  case FEM: return FEM_Conf_Pack_Size (bod);
  }

  return 0;
}

static int velo_pack_size (BODY *bod)
{
  switch (bod->kind)
  {
  case OBS: return 0;
  case RIG: return RIG_VELO_SIZE;
  case PRB: return PRB_VELO_SIZE;
  case FEM: return FEM_Velo_Pack_Size (bod);
  }

  return 0;
}
#endif

/* compute work of contact constraints */
static void compute_contacts_work (BODY *bod, double step)
{
  double DU [3], *R, *energy = bod->energy, coef;
  short dynamic = bod->dom->dynamic;
  SET *item;
  CON *con;

  coef = dynamic ? 0.5 * step : step;
  for (item = SET_First (bod->con); item; item = SET_Next (item))
  {
    con = item->data;
    if (con->kind == CONTACT)
    {
      R = con->R;
      if (dynamic) { ADD (con->U, con->V, DU); }
      else { COPY (con->U, DU); }
      energy [FRICWORK] += coef * DOT2 (DU, R);
      energy [CONTWORK] += coef * DU [2] * R [2];
    }
  }
}

void overwrite_state (BODY *src, BODY *dst)
{
  switch (src->kind)
  {
    case OBS: break;
    case RIG:
    {
      double *cs = RIG_CENTER(src),
             *cd = RIG_CENTER(dst),
	     *Cs = BOD_X0(src),
	     *Cd = BOD_X0(dst),
	      A [3];

      SUB (cs, Cs, A);
      ADD (Cd, A, cd);
      memcpy (dst->conf, src->conf, sizeof (double [9]));
      memcpy (dst->velo, src->velo, sizeof (double [RIG_VELO_SIZE]));
    }
    break;
    case PRB:
    {
      double *cs = PRB_CENTER(src),
             *cd = PRB_CENTER(dst),
	     *Cs = BOD_X0(src),
	     *Cd = BOD_X0(dst),
	      A [3];

      SUB (cs, Cs, A);
      ADD (Cd, A, cd);
      memcpy (dst->conf, src->conf, sizeof (double [9]));
      memcpy (dst->velo, src->velo, sizeof (double [PRB_VELO_SIZE]));
    }
    break;
    case FEM: break;
  }
}

/* -------------- interface ------------- */

BODY* BODY_Create (short kind, SHAPE *shp, BULK_MATERIAL *mat, char *label, BODY_FLAGS flags, short form, MESH *msh)
{
  BODY *bod;

  switch (kind)
  {
    case OBS:
    {
      ERRMEM (bod = MEM_CALLOC (sizeof (BODY)));
    }
    break;
    case RIG:
    {
      double euler [9];

      ERRMEM (bod = alloc_body (RIG));
      bod->conf = (double*) (bod + 1);
      bod->velo = bod->conf + RIG_CONF_SIZE;
      SHAPE_Char (shp, 1,
	&bod->ref_volume,
	 bod->ref_center, euler);
      SCALE9 (euler, mat->density);
      euler2inertia (euler, bod->ref_tensor);
      bod->ref_mass = bod->ref_volume * mat->density;
      bod->id = 0;
      bod->dofs = 6;
      bod->inverse = NULL; 
      bod->forces = NULL;
      
      /* configuration & velocity */
      IDENTITY (RIG_ROTATION(bod));
      COPY (BOD_X0(bod), RIG_CENTER(bod));
      SET (RIG_ANGVEL(bod), 0);
      SET (RIG_LINVEL(bod), 0);
    }
    break;
    case PRB:
    {
      ERRMEM (bod = alloc_body (PRB));
      bod->conf = (double*) (bod + 1);
      bod->velo = bod->conf + PRB_CONF_SIZE;
      SHAPE_Char (shp, 1,
	&bod->ref_volume,
	 bod->ref_center,
	 bod->ref_tensor);
      SCALE9 (bod->ref_tensor, mat->density);
      bod->ref_mass = bod->ref_volume * mat->density;
      bod->id = 0;
      bod->dofs = 12;
      bod->inverse = NULL; 
      bod->forces = NULL;

      /* configuration & velocity */
      IDENTITY (PRB_GRADIENT(bod));
      COPY (BOD_X0(bod), PRB_CENTER(bod));
      SET (PRB_GRADVEL(bod), 0);
      SET (PRB_LINVEL(bod), 0);
    }
    break;
    case FEM:
      ERRMEM (bod = alloc_body (FEM));
      FEM_Create (form, msh, shp, mat, bod);
    break;
    default:
      ASSERT (0, ERR_BOD_KIND);
    break;
  }

  /* set kind */
  bod->kind = kind;

  /* set material */
  bod->mat = mat;

  /* set shape */ 
  bod->shape = shp;

  /* set flags */
  bod->flags = flags;

  /* create piars table */
  SGP_FLAGS sgp_flags = 0;
  if (bod->flags & BODY_DETECT_NODE_CONTACT) sgp_flags |= SGP_MESH_NODES;
  bod->sgp = SGP_Create (shp, &bod->nsgp, sgp_flags);

  /* update body extents */
  SHAPE_Extents (shp, bod->extents);

  /* not in a list */
  bod->next =
  bod->prev = NULL;

  /* set label */
  bod->label = copylabel (label);

  /* default integration scheme */
  bod->scheme = kind == RIG ? SCH_RIG_NEG : SCH_DEF_EXP;

  /* initial damping */
  bod->damping = 0.0;

#if MPI
  bod->children = NULL;

  MPI_Comm_rank (MPI_COMM_WORLD, &bod->rank);
#endif

  return bod;
}

char* BODY_Kind (BODY *bod)
{
  switch (bod->kind)
  {
  case OBS: return "OBSTACLE";
  case RIG: return "RIGID";
  case PRB: return "PSEUDO_RIGID";
  case FEM: return "FINITE_ELEMENT";
  }

  return NULL;
}

int BODY_Conf_Size (BODY *bod)
{
  switch (bod->kind)
  {
  case OBS: return 0;
  case RIG: return 12;
  case PRB: return 12;
  case FEM: return FEM_Conf_Size (bod);
  }

  return 0;
}

void BODY_Overwrite_Chars (BODY *bod, double mass, double volume, double *center, double *tensor)
{
  bod->ref_mass = mass;
  bod->ref_volume = volume;
  COPY (center, bod->ref_center);
  NNCOPY (tensor, bod->ref_tensor);
}

void BODY_Overwrite_State (BODY *bod, double *q, double *u)
{
  switch (bod->kind)
  {
    case OBS: break;
    case RIG:
      memcpy (bod->conf, q, sizeof (double [12]));
      memcpy (bod->velo, u, sizeof (double [6]));
    break;
    case PRB:
      memcpy (bod->conf, q, sizeof (double [12]));
      memcpy (bod->velo, u, sizeof (double [12]));
    break;
    case FEM:
      FEM_Overwrite_State (bod, q, u);
    break;
  }
}

void BODY_Initial_Velocity (BODY *bod, double *linear, double *angular)
{
  switch (bod->kind)
  {
    case OBS: break;
    case RIG:
      if (angular) {COPY (angular, RIG_ANGVEL(bod));}
      if (linear) {COPY (linear, RIG_LINVEL(bod));}
    break;
    case PRB:
      if (angular)
      {
	double copy [3];
	COPY (angular, copy);
	SCALE (copy, -1.0); /* otherwise reversed directions => TODO: why? */
	VECSKEW (copy, PRB_GRADVEL(bod));
      }
      if (linear) {COPY (linear, PRB_LINVEL(bod));}
    break;
    case FEM:
      FEM_Initial_Velocity (bod, linear, angular);
    break;
  }
}

void BODY_Apply_Force (BODY *bod, short kind, double *point, double *direction, TMS *data, void *call, FORCE_FUNC func, int surfid)
{
  FORCE *frc;

  if (kind & PRESSURE)
  {
    ASSERT_DEBUG (data, "NULL pointer passed incorectly");
  }
  else
  {
    ASSERT_DEBUG ((kind & SPATIAL) || (kind & CONVECTED), "Invalid force kind");

    if (kind & TORQUE)
    {
      ASSERT_DEBUG (bod->kind == RIG, "Torque can be only applied to rigid bodies");
      ASSERT_DEBUG ((direction && data) || func, "NULL pointer passed incorectly");
    }
    else
    {
      ASSERT_DEBUG ((point && direction && data) || func, "NULL pointer passed incorectly");
    }
  }

  ERRMEM (frc = MEM_CALLOC (sizeof (FORCE)));

  /* set up force */
  frc->kind = kind;
  if (direction) { NORMALIZE (direction); COPY (direction, frc->direction); }
  if (point) { COPY (point, frc->ref_point); }
  frc->data = data;
  frc->call = call;
  frc->func = func;
  frc->surfid = surfid;
  
  /* append body forces list */
  frc->next = bod->forces;
  bod->forces = frc;
}

void BODY_Clear_Forces (BODY *bod)
{
  FORCE *frc, *nxt;

  /* delete list items */
  for (frc = bod->forces; frc; frc = nxt)
  {
    nxt = frc->next;
    if (frc->data) TMS_Destroy (frc->data);
    free (frc);
  }

  bod->forces = NULL;
}

void BODY_Material (BODY *bod, int volume, BULK_MATERIAL *mat)
{
  SHAPE *shp;
  MESH *msh;
  CONVEX *cvx;
  SPHERE *sph;
  ELLIP *eli;
  ELEMENT *ele;

  for (shp = bod->shape; shp; shp = shp->next)
  {
    switch (shp->kind)
    {
    case SHAPE_MESH:

      msh = shp->data;

      for (ele = msh->surfeles; ele; ele = ele->next)
	if (ele->volume == volume)
	  ele->mat = mat;

      for (ele = msh->bulkeles; ele; ele = ele->next)
	if (ele->volume == volume)
	  ele->mat = mat;

      break;
    case SHAPE_CONVEX:

      cvx = shp->data;

      for (; cvx; cvx = cvx->next)
	if (cvx->volume == volume)
	  cvx->mat = mat;

      break;
    case SHAPE_SPHERE:

      sph = shp->data;
      if (sph->volume == volume)
	sph->mat = mat;

      break;
    case SHAPE_ELLIP:

      eli = shp->data;
      if (eli->volume == volume)
	eli->mat = mat;

      break;
    }
  }
}

void BODY_Dynamic_Init (BODY *bod)
{
  switch (bod->kind)
  {
    case OBS:
      if (!bod->inverse)  /* initialize once */
      {
	bod->inverse = MX_Create (MXDENSE, 3, 3, NULL, NULL);
        MX_Zero (bod->inverse); /* fake zero matrix */
      }
      break;
    case RIG: 
      if (!bod->inverse) rig_dynamic_inverse (bod); /* initialize once */
      break;
    case PRB: 
      if (bod->scheme == SCH_DEF_EXP) 
      {
	if (!bod->inverse) prb_dynamic_explicit_inverse (bod); /* initialize once */
      }
      else prb_dynamic_implicit_inverse (bod, bod->dom->step, NULL, NULL); /* update every time */
      break;
    case FEM:
      FEM_Dynamic_Init (bod);
      break;
  }
}

double BODY_Dynamic_Critical_Step (BODY *bod)
{
  double step = 0.0;

  switch (bod->kind)
  {
    case OBS:
    case RIG:
      step = DBL_MAX;
    break;
    case PRB:
    {
      if (bod->scheme == SCH_DEF_EXP)
      {
	MX_DENSE (K, 9, 9);
	double eigmax;
	MX *IMK;

	SVK_Tangent_R (
	  lambda (bod->mat->young, bod->mat->poisson), mi (bod->mat->young, bod->mat->poisson),
	  bod->ref_volume, 9, bod->conf, K.x); /* calculate stiffness matrix */
	IMK = MX_Matmat (1.0, MX_Diag (bod->inverse, 0, 2), &K, 0.0, NULL); /* inv(M) * K => done on a sub-block of M */
	MX_Eigen (IMK, 1, &eigmax, NULL); /* compute maximal eigenvalue */
	MX_Destroy (IMK);
	ASSERT (eigmax > 0.0, ERR_BOD_MAX_FREQ_LE0);
	step = 2.0 / sqrt (eigmax); /* limit of stability => t_crit <= 2.0 / omega_max */
      }
      else step = DBL_MAX;
    }
    break;
    case FEM:
      step = FEM_Dynamic_Critical_Step (bod);
    break;
  }

  return 0.5 * step; /* XXX: coefficient */
}

void BODY_Dynamic_Step_Begin (BODY *bod, double time, double step)
{
  switch (bod->kind)
  {
    case OBS: break;
    case RIG:
    {
      double half = 0.5 * step,
	     force [6],
	     *J  = BOD_J0(bod),
	     *I  = bod->inverse->x,
	     *W  = RIG_ANGVEL(bod),
	     *W0 = RIG_ANGVEL0(bod),
	     *v  = RIG_LINVEL(bod),
	     *v0 = RIG_LINVEL0(bod),
	     *R  = RIG_ROTATION(bod),
	     *x  = RIG_CENTER(bod),
	     *fext = RIG_FEXT(bod),
	     O [9], DR [9], W05 [3],
	     spatorq [3],
	     reftorq [3];
     
      COPY (W, W0);
      COPY (v, v0);
      COPY (W, O);
      SCALE (O, half);
      EXPMAP (O, DR); 
      NNCOPY (R, O);
      NNMUL (O, DR, R);       /* R(t+h/2) = R(t) exp [(h/2) * W(t)] */
      ADDMUL (x, half, v, x); /* x(t+h/2) = x(t) + (h/2) * v(t) */

      rig_force (bod, bod->conf, bod->velo,
	time+half, step, force + 3, spatorq, reftorq);
      TVADDMUL (reftorq, R, spatorq, force);

      if (bod->scheme > SCH_RIG_POS)
      {
	double *A = RIG_AUXILIARY(bod);

	NVMUL (J, W, O);
	TVMUL (DR, O, A);
	ADDMUL (A, step, force, A);  /* exp [-(h/2)W(t)]JW(t) + hT(t+h/2) */
      }

      COPY (force, O);
      SCALE (O, half); 
      NVMUL (J, W, W05);
      TVADDMUL (O, DR, W05, O);
      NVMUL (I, O, W05);            /* W05 = I * [exp [-(h/2) * W(t)]*J*W(t) + (h/2)*T(t+h/2)] */

      NVMUL (J, W05, O);
      PRODUCTSUB (W05, O, force); /* force [0..2] = T - W05 x J * W05 */
      COPY6 (force, fext); /* fext = force */
      
      MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) = u(t) + inv (M) * h * f(t+h/2) */
    }
    break;
    case PRB:
    {
      double half = 0.5 * step,
	     *L  = PRB_GRADVEL(bod),
	     *L0 = PRB_GRADVEL0(bod),
	     *v  = PRB_LINVEL(bod),
	     *v0 = PRB_LINVEL0(bod),
	     force [12];

      NNCOPY (L, L0);
      COPY (v, v0);

      switch (bod->scheme)
      {
      case SCH_DEF_EXP:
      {
	blas_daxpy (12, half, bod->velo, 1, bod->conf, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
	prb_dynamic_force (bod, time+half, step, PRB_FEXT(bod), PRB_FINT(bod), force);  /* f(t+h/2) = fext (t+h/2) - fint (q(t+h/2)) */
	if (bod->damping > 0.0)
	{
	  prb_tangent_stiffness (bod);
	  MX_Matvec (-bod->damping, bod->K, bod->velo, 1.0, force); /* f -= damping K u (t) */
	}
	MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) = u(t) + inv (M) * h * f(t+h/2) */
      }
      break;
      case SCH_DEF_LIM:
      {
	blas_daxpy (12, half, bod->velo, 1, bod->conf, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
	prb_dynamic_force (bod, time+half, step, PRB_FEXT(bod), PRB_FINT(bod), force);  /* f = fext (t+h/2) - fint (q(t+h/2)) */
	prb_dynamic_implicit_inverse (bod, step, NULL, NULL); /* A = M + (h*h/4) * K */
        if (bod->damping > 0.0) MX_Matvec (-bod->damping, bod->K, bod->velo, 1.0, force); /* f -= damping K u (t) */
	MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) = u (t) + inv (A) * h * f */
      }
      break;
      default:
      break;
      }
    }
    break;
    case FEM:
      FEM_Dynamic_Step_Begin (bod, time, step);
    break;
  }

  /* update shape */
  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point);
  if (bod->msh) FEM_Update_Rough_Mesh (bod);
}

void BODY_Dynamic_Step_End (BODY *bod, double time, double step)
{
  double *energy = bod->energy;

  switch (bod->kind)
  {
    case OBS: break;
    case RIG:
    {
      SCHEME sch = bod->scheme;
      double half = 0.5 * step,
	     force [6],
	     *x = RIG_CENTER(bod),
	     *v = RIG_LINVEL(bod),
	     *R = RIG_ROTATION(bod),
	     *W = RIG_ANGVEL(bod),
	     *velo = bod->velo,
	     *vel0 = RIG_ANGVEL0(bod),
	     *fext = RIG_FEXT(bod),
	     O [9], DR [9];

      rig_constraints_force (bod, force); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
      MX_Matvec (step, bod->inverse, force, 1.0, velo); /* u(t+h) += inv (M) * h * r */
      ADDMUL (x, half, v, x); /* x(t+h) = x(t+h/2) + (h/2) * v(t+h) */

      if (sch <= SCH_RIG_NEG)
      {
	COPY (W, O);
	SCALE (O, half);
	EXPMAP (O, DR); 
	NNCOPY (R, O);
	NNMUL (O, DR, R); /* R(t) = R(t+h/2) exp [(h/2) * W(t+h)] */
      }

      if (sch > SCH_RIG_POS)
      {
	double *A = RIG_AUXILIARY(bod);
	
	ADDMUL (A, step, force, A);

	if (sch == SCH_RIG_NEG)
	{
          double *I = bod->inverse->x;

	  TVMUL (DR, A, O);
	  NVMUL (I, O, W); /* W2(t+h) = I * exp [-(h/2)W1(t+h)][exp [-(h/2)W(t)]JW(t) + hT(t+h/2)] */
	}
	else
	{
          double *J  = BOD_J0(bod);

	  SOLVE (half, J, W, A, DR);
	  NNCOPY (R, O);
	  NNMUL (O, DR, R); /* R(t) = R(t+h/2) exp [(h/2) * W(t+h)] */
	}
      }

      /* energy */
      ACC6 (force, fext);
      ADD6 (velo, vel0, force);
      SCALE6 (force, half); /* dq = (h/2) * {u(t) + u(t+h)} */
      energy [EXTERNAL] += DOT6 (force, fext);
    }
    break;
    case PRB:
    {
      double half = 0.5 * step,
	    *fext = PRB_FEXT(bod),
	    *fint = PRB_FINT(bod),
	    *velo = PRB_GRADVEL(bod),
	    *vel0 = PRB_GRADVEL0(bod),
             force [12], dq [12];

      prb_constraints_force (bod, force); /* force = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
      ACC12 (force, fext); /* fext += H^T R */

      switch (bod->scheme)
      {
      default: /* SCH_DEF_EXP, SCH_DEF_LIM */
      {
        MX_Matvec (step, bod->inverse, force, 1.0, velo); /* u(t+h) += h * inv (M) * force */
        blas_daxpy (12, half, velo, 1, bod->conf, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
      }
      break;
      }

      /* energy */
      ADD12 (velo, vel0, dq);
      SCALE12 (dq, half); /* dq = (h/2) * {u(t) + u(t+h)} */
      energy [EXTERNAL] += DOT12 (dq, fext);
      energy [INTERNAL] += DOT9 (dq, fint);
      /* computing internal energy like above may produce negative energy increments during
       * impacts since fint is computed at q(t+h/2) whereas dq includes impact correction; 
       * this is effect is present when the time integration step is excessively large */
    }
    break;
    case FEM:
      FEM_Dynamic_Step_End (bod, time, step);
    break;
  }

  if (bod->kind != OBS)
  {
    energy [KINETIC] = BODY_Kinetic_Energy (bod);

    compute_contacts_work (bod, step); /* CONTWORK and FRICWORK */

#if 0 /* TODO: remove condition when proved robust */
    double emax, etot;

    MAXABS (energy, emax);

    if (emax > DBL_EPSILON && bod->damping == 0.0) /* discard tiny energy balance and damping */
    {
      etot = energy[KINETIC] + energy[INTERNAL] - energy[EXTERNAL];

      if (!(etot < ENE_TOL * emax || emax < 0))
	fprintf (stderr, "KIN = %g, INT = %g, EXT = %g, TOT = %g\n", energy [KINETIC], energy [INTERNAL], energy [EXTERNAL], etot);

      ASSERT (etot < ENE_TOL * emax || emax < 0, ERR_BOD_ENERGY_CONSERVATION);
    }
#endif

    /* update shape */
    SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point);
    if (bod->msh) FEM_Update_Rough_Mesh (bod);
  }
}

void BODY_Static_Init (BODY *bod)
{
  switch (bod->kind)
  {
    case OBS:
      if (!bod->inverse) 
      {
	bod->inverse =
	MX_Create (MXDENSE, 3, 3, NULL, NULL);
        MX_Zero (bod->inverse); /* fake zero matrix */
      }
      break;
    case RIG:
    {
      double *W = RIG_ANGVEL(bod),
	     *v = RIG_LINVEL(bod);
     
      /* cancel initial
       * velocities */ 
      SET (W, 0);
      SET (v, 0);
     
      /* initialise inverse */ 
      rig_static_inverse (bod);
    }
    break;
    case PRB:
    {
      double *L = PRB_GRADVEL(bod),
	     *v = PRB_LINVEL(bod);
     
      /* cancel initial
       * velocities */ 
      SET9 (L, 0);
      SET (v, 0);

      /* initialise inverse */ 
      prb_static_inverse (bod, bod->dom->step);
    }
    break;
    case FEM:
      FEM_Static_Init (bod);
    break;
  }
}

void BODY_Static_Step_Begin (BODY *bod, double time, double step)
{
  switch (bod->kind)
  {
    case OBS: break;
    case RIG:
    {
      double force [6];
     
      rig_static_force (bod, time+step, step, force);  /* f(t+h) = [R(t+h)^T * torque (t+h) - W(t+h) x J * W(t+h); linforc (t+h)] */
      MX_Matvec (step, bod->inverse, force, 0.0, bod->velo); /* u(t+h) = inv (M) * h * f(t+h) */
    }
    break;
    case PRB:
    {
      double force [12];

      prb_static_inverse (bod, step); /* compute inverse of static tangent operator */
      prb_static_force (bod, time+step, step, PRB_FEXT(bod), PRB_FINT(bod), force);  /* f(t+h) = fext (t+h) - fint (q(t+h)) */
      MX_Matvec (step, bod->inverse, force, 0.0, bod->velo); /* u(t+h) = inv (A) * h * f(t+h) */
    }
    break;
    case FEM:
      FEM_Static_Step_Begin (bod, time, step);
    break;
  }
}

void BODY_Static_Step_End (BODY *bod, double time, double step)
{
  switch (bod->kind)
  {
    case OBS: break;
    case RIG:
    {
      double force [6],
	     *x = RIG_CENTER(bod),
	     *v = RIG_LINVEL(bod),
	     *R = RIG_ROTATION(bod),
	     *W = RIG_ANGVEL(bod),
	     O [9], DR [9];

      rig_constraints_force (bod, force); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
      MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) += inv (M) * h * r */
      ADDMUL (x, step, v, x); /* x(t+h) = x(t) + h * v(t+h) */
      COPY (W, O);
      SCALE (O, step);
      EXPMAP (O, DR); 
      NNCOPY (R, O);
      NNMUL (O, DR, R); /* R(t) = R(t) exp [h * W(t+h)] */
    }
    break;
    case PRB:
    {
      double force [12];
      
      prb_constraints_force (bod, force); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
      MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) += inv (A) * h * r */
      blas_daxpy (12, step, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t) + h * u(t+h) */
    }
    break;
    case FEM:
      FEM_Static_Step_End (bod, time, step);
    break;
  }

  /* update shape */
  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point);
  if (bod->msh) FEM_Update_Rough_Mesh (bod);
}

void BODY_Update_Extents (BODY *bod)
{
  double *e = bod->extents;

  SHAPE_Extents (bod->shape, e);

#if MPI
  double d [3], *p;
  SET *item;
  CON *con;

  /* now make sure that all attached constraint points are within the extents;
   * due to roundof this could be compromised sometimes, so that constraints might
   * migrate to partitions where bodies would have no representation (child/parent) */

  for (item = SET_First (bod->con); item; item = SET_Next (item))
  {
    con = item->data;
    p = con->point;
    if (p [0] < e [0]) e [0] = p [0];
    if (p [1] < e [1]) e [1] = p [1];
    if (p [2] < e [2]) e [2] = p [2];
    if (p [0] > e [3]) e [3] = p [0];
    if (p [1] > e [4]) e [4] = p [1];
    if (p [2] > e [5]) e [5] = p [2];
  }

 /* extend extents of the body as it might have been altered above;
  * note that we are not using GEOMETRIC_EPSILON here as this can be set
  * by the user which in turn may cause migration consitency problems */
  SUB (e+3, e, d);
  SCALE (d, PUT_GEOMEPS);
  SUB (e, d, e);
  ADD (e+3, d, e+3);
#endif
}

void BODY_Cur_Point (BODY *bod, SGP *sgp, double *X, double *x)
{
  switch (bod->kind)
  {
    case OBS:
      COPY (X, x);
    break;
    case RIG:
    {
      double *R = RIG_ROTATION(bod),
	     *c = RIG_CENTER(bod),
	     *C = BOD_X0(bod),
	     A [3];
      SUB (X, C, A);
      NVADDMUL (c, R, A, x);
    }
    break;
    case PRB:
    {
      double *F = PRB_GRADIENT(bod),
	     *c = PRB_CENTER(bod),
	     *C = BOD_X0(bod),
	     A [3];

      SUB (X, C, A);
      TVADDMUL (c, F, A, x); /* transpose, as F is stored row-wise */
    }
    break;
    case FEM:
    {
      FEM_Cur_Point (bod, sgp->shp, SGP_2_GOBJ (sgp), X, x);
    }
    break;
  }
}

void BODY_Ref_Point (BODY *bod, SGP *sgp, double *x, double *X)
{
  switch (bod->kind)
  {
    case OBS:
      COPY (x, X);
    break;
    case RIG:
    {
      double *R = RIG_ROTATION(bod),
	     *c = RIG_CENTER(bod),
	     *C = BOD_X0(bod),
	     a [3];
      SUB (x, c, a);
      TVADDMUL (C, R, a, X);
    }
    break;
    case PRB:
    {
      double *F = PRB_GRADIENT(bod),
	     *c = PRB_CENTER(bod),
	     *C = BOD_X0(bod),
	     FT [9], IF [9], a [3], det;

      SUB (x, c, a);
      NTCOPY (F, FT);
      INVERT (FT, IF, det);
      ASSERT (det > 0.0, ERR_BOD_MOTION_INVERT);
      NVADDMUL (C, IF, a, X);
    }
    break;
    case FEM:
      FEM_Ref_Point (bod, sgp->shp, SGP_2_GOBJ (sgp), x, X);
    break;
  }
}

void BODY_Cur_Vector (BODY *bod, void *ele, double *X, double *V, double *v)
{
  switch (bod->kind)
  {
    case OBS:
      COPY (V, v);
    break;
    case RIG:
    {
      double *R = RIG_ROTATION(bod);
      NVMUL (R, V, v);
    }
    break;
    case PRB:
    {
      double *F = PRB_GRADIENT(bod);
      NVMUL (F, V, v);
    }
    break;
    case FEM:
      FEM_Cur_Vector (bod, ele, X, V, v);
    break;
  }
}

void BODY_Ref_Vector (BODY *bod, void *ele, double *x, double *v, double *V)
{
  switch (bod->kind)
  {
    case OBS:
      COPY (v, V);
    break;
    case RIG:
    {
      double *R = RIG_ROTATION(bod);
      TVMUL (R, v, V);
    }
    break;
    case PRB:
    {
      double *F = PRB_GRADIENT(bod);
      TVMUL (F, v, V);
    }
    break;
    case FEM:
      FEM_Ref_Vector (bod, ele, x, v, V);
    break;
  }
}

void BODY_Local_Velo (BODY *bod, SGP *sgp, double *point, double *base, double *prevel, double *curvel)
{
  switch (bod->kind)
  {
    case OBS:
      if (prevel) SET (prevel, 0.0);
      if (curvel) SET (curvel, 0.0);
    break;
    case RIG:
    {
      double H [18];

      rig_operator_H (bod, point, base, H);
      if (prevel) blas_dgemv ('N', 3, 6, 1.0, H, 3, bod->velo + 6, 1, 0.0, prevel, 1);
      if (curvel) blas_dgemv ('N', 3, 6, 1.0, H, 3, bod->velo, 1, 0.0, curvel, 1);
    }
    break;
    case PRB:
    {
      double H [36];

      prb_operator_H (bod, point, base, H);
      if (prevel) blas_dgemv ('N', 3, 12, 1.0, H, 3, bod->velo + 12, 1, 0.0, prevel, 1);
      if (curvel) blas_dgemv ('N', 3, 12, 1.0, H, 3, bod->velo, 1, 0.0, curvel, 1);
    }
    break;
    case FEM:
      FEM_Local_Velo (bod, sgp->shp, SGP_2_GOBJ (sgp), point, base, prevel, curvel);
    break;
  }
}

MX* BODY_Gen_To_Loc_Operator (BODY *bod, SGP *sgp, double *point, double *base)
{
  MX *H = NULL;

  switch (bod->kind)
  {
    case OBS:
      H = MX_Create (MXDENSE, 3, 3, NULL, NULL); /* fake zero matrix */
      MX_Zero (H);
      break;
    case RIG:
      H = MX_Create (MXDENSE, 3, 6, NULL, NULL);
      rig_operator_H (bod, point, base, H->x);
    break;
    case PRB:
      H = MX_Create (MXDENSE, 3, 12, NULL, NULL);
      prb_operator_H (bod, point, base, H->x);
    break;
    case FEM:
      H = FEM_Gen_To_Loc_Operator (bod, sgp->shp, SGP_2_GOBJ (sgp), point, base);
    break;
  }

  return H;
}

double BODY_Kinetic_Energy (BODY *bod)
{
  double energy = 0.0;

  switch (bod->kind)
  {
    case OBS: break;
    case RIG:
    {
      double *J = BOD_J0(bod),
	      m = BOD_M0(bod),
	     *W = RIG_ANGVEL(bod),
	     *v = RIG_LINVEL(bod),
	     JW [3];

      NVMUL (J, W, JW);
      energy = 0.5 * (DOT(W, JW) + m*DOT(v, v));
    }
    break;
    case PRB:
    {
      double *E = BOD_E0(bod),
	      m = BOD_M0(bod),
	     *L = PRB_GRADVEL(bod),
	     *v = PRB_LINVEL(bod),
	     EL [9];

      NNMUL (E, L, EL);
      energy = 0.5 * (DOT9(L, EL) + m*DOT(v, v));
    }
    break;
    case FEM:
      energy = FEM_Kinetic_Energy (bod);
    break;
  }

  return energy;
}

void BODY_Point_Values (BODY *bod, double *point, VALUE_KIND kind, double *values)
{
  switch (bod->kind)
  {
  case OBS: break;
  case RIG:
  case PRB:
  {
    switch (kind)
    {
    case VALUE_COORD:
    {
      BODY_Cur_Point (bod, NULL, point, values);
    }
    break;
    case VALUE_DISPLACEMENT:
    {
      double cur_point [3];

      BODY_Cur_Point (bod, NULL, point, cur_point);
      SUB (cur_point, point, values);
    }
    break;
    case VALUE_VELOCITY:
    {
      double base [9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

      BODY_Local_Velo (bod, NULL, point, base, NULL, values);
    }
    break;
    case VALUE_STRESS:
    {
      if (bod->kind == PRB)
	prb_cauchy (bod, values);
    }
    break;
    case VALUE_MISES:
    {
      double stress [6];

      if (bod->kind == PRB)
      {
	prb_cauchy (bod, stress);
	MISES (stress, values [0]);
      }
    }
    break;
    case VALUE_STRESS_AND_MISES:
    {
      if (bod->kind == PRB)
      {
	prb_cauchy (bod, values);
	MISES (values, values [6]);
      }
    }
    break;
    }
  }
  break;
  case FEM:
    FEM_Point_Values (bod, NULL, point, kind, values);
  break;
  }
}

void BODY_Split (BODY *bod, double *point, double *normal, short topoadj, int surfid, BODY **one, BODY **two)
{
  switch (bod->kind)
  {
  case OBS:
  case RIG:
  case PRB:
    {
      SHAPE *copy, *sone, *stwo;
      char *label;

      *one = *two = NULL;
      copy = SHAPE_Copy (bod->shape);
      SHAPE_Update (copy, NULL, NULL); /* restore reference configuration */
      SHAPE_Split (copy, point, normal, topoadj, surfid, &sone, &stwo); /* split in reference configuration */
      SHAPE_Destroy (copy);

      if (bod->label) ERRMEM (label = malloc (strlen (bod->label) + 8));
      else label = NULL;

      if (sone)
      {
	if (bod->label) sprintf (label, "%s/1", bod->label);
	(*one) = BODY_Create (bod->kind, sone, bod->mat, label, bod->flags & BODY_PERMANENT_FLAGS, 0, NULL);
	overwrite_state (bod, *one);
        SHAPE_Update ((*one)->shape, (*one), (MOTION)BODY_Cur_Point); 
      }

      if (stwo)
      {
	if (bod->label) sprintf (label, "%s/2", bod->label);
	(*two) = BODY_Create (bod->kind, stwo, bod->mat, label, bod->flags & BODY_PERMANENT_FLAGS, 0, NULL);
	overwrite_state (bod, *two);
        SHAPE_Update ((*two)->shape, (*two), (MOTION)BODY_Cur_Point); 
      }

      if (bod->label) free (label);
    }
    break;
  case FEM:
    FEM_Split (bod, point, normal, topoadj, surfid, one, two);
    break;
  }

  BODY *out [] = {*one, *two};
  for (int i = 0; i < 2; i ++)
  {
    if (out [i])
    {
      COPY6 (bod->extents, out[i]->extents);
      out [i]->scheme = bod->scheme;
      out [i]->damping = bod->damping;
    }
  }

  if (out [0] && out [1])
  {
    double votot = out [0]->ref_volume + out [1]->ref_volume,
           coef [] = {out [0]->ref_volume / votot, out [1]->ref_volume / votot};

    for (int i = 0; i < 2; i ++)
      for (int j = 0; j < BODY_ENERGY_SPACE; j ++)
        out [i]->energy [j] = coef [i] * bod->energy [j]; /* XXX: volume proportional split is not best in all cases */
  }

  /* TODO: transfer forces to parts */
}

BODY** BODY_Separate (BODY *bod, int *m)
{
  BODY **out = NULL;
  *m = 0;

  switch (bod->kind)
  {
  case OBS:
  case RIG:
  case PRB:
    {
      SHAPE **shp;
      char *label;

      shp = SHAPE_Separate (bod->shape, m);

      if (shp)
      {
	ERRMEM (out = malloc ((*m) * sizeof (BODY*)));

	if (bod->label) ERRMEM (label = malloc (strlen (bod->label) + 8));
	else label = NULL;

	for (int i = 0; i < (*m); i ++)
	{
	  if (bod->label) sprintf (label, "%s/%d", bod->label, i);
	  out [i] = BODY_Create (bod->kind, shp [i], bod->mat, label, bod->flags & BODY_PERMANENT_FLAGS, 0, NULL);
	  overwrite_state (bod, out [i]);
	  SHAPE_Update (out [i]->shape, out [i], (MOTION)BODY_Cur_Point); 
	}

	if (bod->label) free (label);
	free (shp);
      }
    }
    break;
  case FEM:
    out = FEM_Separate (bod, m);
    break;
  }

  double vtot = 0;

  for (int i = 0; i < (*m); i ++)
  {
    COPY6 (bod->extents, out[i]->extents);
    out [i]->scheme = bod->scheme;
    out [i]->damping = bod->damping;

    vtot += out [i]->ref_volume;
  }

  for (int i = 0; i < (*m); i ++)
  {
    double coef = out [i]->ref_volume / vtot;

    for (int j = 0; j < BODY_ENERGY_SPACE; j ++)
      out [i]->energy [j] = coef * bod->energy [j]; /* XXX: volume proportional split is not best in all cases */
  }

  /* TODO: transfer forces to parts */

  return out;
}

void BODY_Write_State (BODY *bod, PBF *bf)
{
  PBF_Double (bf, bod->conf, BODY_Conf_Size (bod));
  PBF_Double (bf, bod->velo, bod->dofs);
  PBF_Double (bf, bod->energy, BODY_ENERGY_SIZE(bod));

#if MPI
  PBF_Int (bf, &bod->dom->rank, 1);
#endif
}

void BODY_Read_State (BODY *bod, PBF *bf)
{
  PBF_Double (bf, bod->conf, BODY_Conf_Size (bod));
  PBF_Double (bf, bod->velo, bod->dofs);
  PBF_Double (bf, bod->energy, BODY_ENERGY_SIZE(bod));

  if (bf->parallel == PBF_ON)
  {
    if (bod->dom->solfec->mode == SOLFEC_READ) /* XXX: longish ->->-> */
    {
      PBF_Int (bf, &bod->rank, 1);
    }
    else /* fake it => ranks are actually used in WRITE mode */
    {
      int rank;
      PBF_Int (bf, &rank, 1); /* branching down here from Python's INITIALIZE_STATE */
    }
  }

  /* update shape */
  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point); 
  if (bod->msh) FEM_Update_Rough_Mesh (bod);
}

void BODY_Destroy (BODY *bod)
{
  FORCE *forc, *next;

  for (forc = bod->forces; forc; forc = next)
  { 
    next = forc->next;

    if (forc->data && !forc->func) TMS_Destroy (forc->data);

    free (forc);
  }

  SHAPE_Destroy (bod->shape);

  free (bod->sgp);

  CRACK_Destroy_List (bod->cra);

  if (bod->inverse) MX_Destroy (bod->inverse);

  if (bod->M) MX_Destroy (bod->M);

  if (bod->K) MX_Destroy (bod->K);

  if (bod->kind == FEM) FEM_Destroy (bod);

  if (bod->msh) MESH_Destroy (bod->msh);

  if (bod->eval) free (bod->eval);

  if (bod->evec) MX_Destroy (bod->evec);

#if OPENGL
  if (bod->rendering) RND_Free_Rendering_Data (bod->rendering);
#endif

  free (bod);
}

/* pack force list */
static void pack_forces (FORCE *forces, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  FORCE *frc;
  int count, id;

  for (count = 0, frc = forces; frc; frc = frc->next) count ++;

  pack_int (isize, i, ints, count);

  for (frc = forces; frc; frc = frc->next)
  {
    pack_int (isize, i, ints, frc->kind);
    pack_doubles (dsize, d, doubles, frc->ref_point, 3);
    pack_doubles (dsize, d, doubles, frc->direction, 3);

    if (frc->func)
    {
      pack_int (isize, i, ints, 1);
      ASSERT_DEBUG_EXT (id = lngcallback_id (frc->data, frc->call), "Invalid callback pair");
      pack_int (isize, i, ints, id);
    }
    else
    {
      pack_int (isize, i, ints, 0);
      TMS_Pack (frc->data, dsize, d, doubles, isize, i, ints);
    }

    if (frc->kind == PRESSURE) pack_int (isize, i, ints, frc->surfid);
  }
}

/* unpack force list */
static FORCE* unpack_forces (int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  FORCE *forces, *frc;
  int count, id, n, func;

  forces = NULL;

  count = unpack_int (ipos, i, ints);

  for (n = 0; n < count; n ++)
  {
    ERRMEM (frc = MEM_CALLOC (sizeof (FORCE)));

    frc->kind = unpack_int (ipos, i, ints);
    unpack_doubles (dpos, d, doubles, frc->ref_point, 3);
    unpack_doubles (dpos, d, doubles, frc->direction, 3);

    func = unpack_int (ipos, i, ints);

    if (func)
    {
      id = unpack_int (ipos, i, ints);
      ASSERT_DEBUG_EXT (lngcallback_set (id, (void**) &frc->data, &frc->call), "Invalid callback id");
    }
    else
    {
      frc->data = TMS_Unpack (dpos, d, doubles, ipos, i, ints);
    }

    if (frc->kind == PRESSURE) frc->surfid = unpack_int (ipos, i, ints);

    frc->next = forces;
    forces = frc;
  }

  return forces;
}

/* pack body */
void BODY_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  /* these are arguments of BODY_Create */
  pack_int (isize, i, ints, bod->kind);
  pack_int (isize, i, ints, bod->flags & BODY_PERMANENT_FLAGS);
  if (bod->kind == FEM) pack_int (isize, i, ints, bod->msh ? -bod->form : bod->form);
  if (bod->msh) MESH_Pack (bod->msh, dsize, d, doubles, isize, i, ints);
  SHAPE_Pack (bod->shape, dsize, d, doubles, isize, i, ints);
  pack_string (isize, i, ints, bod->mat->label);
  pack_string (isize, i, ints, bod->label);

  /* characteristics will be overwritten when unpacking */
  pack_double (dsize, d, doubles, bod->ref_mass);
  pack_double (dsize, d, doubles, bod->ref_volume);
  pack_doubles (dsize, d, doubles, bod->ref_center, 3);
  pack_doubles (dsize, d, doubles, bod->ref_tensor, 9);

  /* body id */
  pack_int (isize, i, ints, bod->id);

  /* configuration and velocity */
  pack_doubles (dsize, d, doubles, bod->conf, BODY_Conf_Size (bod));
  pack_doubles (dsize, d, doubles, bod->velo, bod->dofs);

  /* pack the list of forces */
  pack_forces (bod->forces, dsize, d, doubles, isize, i, ints);

  /* pack scheme */
  pack_int (isize, i, ints, bod->scheme);

  /* damping */
  pack_double (dsize, d, doubles, bod->damping);

  /* pack cracks */
  CRACKS_Pack (bod->cra, dsize, d, doubles, isize, i, ints); 

  /* modal analysis */
  if (bod->eval && bod->evec)
  {
    pack_int (isize, i, ints, bod->evec->n);
    pack_doubles (dsize, d, doubles, bod->eval, bod->evec->n);
    MX_Pack (bod->evec, dsize, d, doubles, isize, i, ints);
  }
  else pack_int (isize, i, ints, 0);
}

/* unpack body */
BODY* BODY_Unpack (SOLFEC *sol, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  BODY *bod;
  int kind;
  MESH *msh;
  SHAPE *shp;
  char *label;
  BULK_MATERIAL *mat;
  BODY_FLAGS flags;
  short form;

  /* unpack BODY_Create arguments and create body */
  kind = unpack_int (ipos, i, ints);
  flags = unpack_int (ipos, i, ints);
  if (kind == FEM) form = unpack_int (ipos, i, ints); else form = 0;
  if (form < 0) msh = MESH_Unpack (sol, dpos, d, doubles, ipos, i, ints), form = -form; else msh = NULL;
  shp = SHAPE_Unpack (sol, dpos, d, doubles, ipos, i, ints);
  label = unpack_string (ipos, i, ints);
  ASSERT_DEBUG_EXT (mat = MATSET_Find (sol->mat, label), "Invalid bulk material label");
  free (label);
  label = unpack_string (ipos, i, ints);
  bod = BODY_Create (kind, shp, mat, label, flags, form, msh);
  free (label);

  /* overwritte characteristics */
  bod->ref_mass = unpack_double (dpos, d, doubles);
  bod->ref_volume = unpack_double (dpos, d, doubles);
  unpack_doubles (dpos, d, doubles, bod->ref_center, 3);
  unpack_doubles (dpos, d, doubles, bod->ref_tensor, 9);

  /* body id */
  bod->id = unpack_int (ipos, i, ints);

  /* configuration and velocity */
  unpack_doubles (dpos, d, doubles, bod->conf, BODY_Conf_Size (bod));
  unpack_doubles (dpos, d, doubles, bod->velo, bod->dofs);

  /* unpack the list of forces */
  bod->forces = unpack_forces (dpos, d, doubles, ipos, i, ints);

  /* unpack scheme */
  bod->scheme = unpack_int (ipos, i, ints);

  /* damping */
  bod->damping = unpack_double (dpos, d, doubles);

  /* unpack cracks */
  bod->cra = CRACKS_Unpack (dpos, d, doubles, ipos, i, ints);

  /* modal analysis */
  int modal = unpack_int (ipos, i, ints);
  if (modal)
  {
    ERRMEM (bod->eval = malloc (sizeof (double [modal])));
    unpack_doubles (dpos, d, doubles, bod->eval, modal);
    bod->evec = MX_Unpack (dpos, d, doubles, ipos, i, ints);
  }

  return bod;
}

#if MPI
/* pack parent body */
void BODY_Parent_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  /* configuration and velocity */
  pack_doubles (dsize, d, doubles, bod->conf, conf_pack_size (bod));
  pack_doubles (dsize, d, doubles, bod->velo, velo_pack_size (bod));

  /* pack children ranks */
  pack_int (isize, i, ints, SET_Size (bod->children));
  for (SET *item = SET_First (bod->children); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int) (long) item->data);

  /* pack scheme and flags */
  pack_int (isize, i, ints, bod->scheme);
  pack_int (isize, i, ints, bod->flags & BODY_PERMANENT_FLAGS);

  /* damping */
  pack_double (dsize, d, doubles, bod->damping);

  /* pack energy */
  pack_doubles (dsize, d, doubles, bod->energy, BODY_ENERGY_SIZE(bod));
}

/* unpack parent body */
void BODY_Parent_Unpack (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  short dynamic = bod->dom->dynamic;
  int j, k;

  /* configuration and velocity */
  unpack_doubles (dpos, d, doubles, bod->conf, conf_pack_size (bod));
  unpack_doubles (dpos, d, doubles, bod->velo, velo_pack_size (bod));

  /* unpack children ranks */
  SET_Free (&bod->dom->setmem, &bod->children);
  k = unpack_int (ipos, i, ints);
  for (j = 0; j < k; j ++)
    SET_Insert (&bod->dom->setmem, &bod->children, (void*) (long) unpack_int (ipos, i, ints), NULL);

  /* unpack scheme and flags */
  bod->scheme = unpack_int (ipos, i, ints);
  bod->flags |= unpack_int (ipos, i, ints);

  /* damping */
  bod->damping = unpack_double (dpos, d, doubles);

  /* unpack energy */
  unpack_doubles (dpos, d, doubles, bod->energy, BODY_ENERGY_SIZE(bod));

  /* init inverse */
  if (dynamic) BODY_Dynamic_Init (bod);
  else BODY_Static_Init (bod);

  /* update shape */
  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point); 
  if (bod->msh) FEM_Update_Rough_Mesh (bod);
}

/* pack child body */
void BODY_Child_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, bod->rank); /* pack parent rank */

  pack_int (isize, i, ints, SET_Size (bod->children)); /* pack children ranks */
  for (SET *item = SET_First (bod->children); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int) (long) item->data);

  pack_doubles (dsize, d, doubles, bod->conf, conf_pack_size (bod)); /* configuration */

  pack_int (isize, i, ints, bod->scheme); /* pack integration scheme */
  
  pack_double (dsize, d, doubles, bod->damping); /* damping */
}

/* unpack child body */
void BODY_Child_Unpack (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  short dynamic = bod->dom->dynamic;
  int j, k, l;

  bod->rank = unpack_int (ipos, i, ints); /* unpack parent rank */
 
  SET_Free (&bod->dom->setmem, &bod->children); /* unpack other children ranks */
  k = unpack_int (ipos, i, ints);
  for (j = 0; j < k; j ++)
  {
    l = unpack_int (ipos, i, ints);
    if (l != bod->dom->rank) /* ommit own rank */
      SET_Insert (&bod->dom->setmem, &bod->children, (void*) (long) l, NULL);
  }

  unpack_doubles (dpos, d, doubles, bod->conf, conf_pack_size (bod)); /* configuration */

  bod->scheme = unpack_int (ipos, i, ints);  /* unpack integration scheme */

  bod->damping = unpack_double (dpos, d, doubles); /* damping */

  /* init inverse */
  if (dynamic) BODY_Dynamic_Init (bod);
  else BODY_Static_Init (bod);

  /* update shape */
  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point); /* update shape */
  if (bod->msh) FEM_Update_Rough_Mesh (bod);
}

/* pack child update */
void BODY_Child_Update_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_doubles (dsize, d, doubles, bod->conf, conf_pack_size (bod));
}

/* unpack child update */
void BODY_Child_Update_Unpack (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  short dynamic = bod->dom->dynamic;

  unpack_doubles (dpos, d, doubles, bod->conf, conf_pack_size (bod));

  /* init inverse */
  if (dynamic) BODY_Dynamic_Init (bod);
  else BODY_Static_Init (bod);

  /* update shape */
  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point); 
  if (bod->msh) FEM_Update_Rough_Mesh (bod);
}
#endif

/* compute c = alpha * INVERSE (bod) * b + beta * c */
void BODY_Invvec (double alpha, BODY *bod, double *b, double beta, double *c)
{
  switch (bod->kind)
  {
    case RIG:
    case PRB:
      MX_Matvec (alpha, bod->inverse, b, beta, c);
      break;
    case FEM:
      FEM_Invvec (alpha, bod, b, beta, c);
      break;
    case OBS:
      break;
  }
}

/* clone body by first rotating (point, vector, angle) it and then translating */
BODY* BODY_Clone (BODY *bod, double *translate, double *point, double *vector, double angle, char *label)
{
  SHAPE *shp;
  MESH *msh;
  BODY *out;

  shp = SHAPE_Copy (bod->shape);
  if (bod->msh) msh = MESH_Copy (bod->msh);
  else msh = NULL;

  if (angle != 0.0)
  {
    SHAPE_Rotate (shp, point, vector, angle);
    if (msh) MESH_Rotate (msh, point, vector, angle);
  }

  SHAPE_Translate (shp, translate);
  if (msh) MESH_Translate (msh, translate);

  out = BODY_Create (bod->kind, shp, bod->mat, label, bod->flags, bod->form, msh);

  out->scheme = bod->scheme;
  out->damping = bod->damping;

  if (bod->eval && bod->evec) /* modal analysis results */
  {
    ERRMEM (out->eval = malloc (sizeof (double [bod->evec->n])));
    blas_dcopy (bod->evec->n, bod->eval, 1, out->eval, 1);
    out->evec = MX_Copy (bod->evec, NULL);

    if (angle != 0.0)
    {
      double R [9], q [3], *x, *y;
      int i;

      ROTATION_MATRIX (vector, angle, R);

      for (i = 0; i < out->evec->n; i ++)
      {
	for (x = &out->evec->x[i*out->evec->m], y = x + out->evec->m; x < y; x += 3)
	{
	  COPY (x, q);
	  NVMUL (R, q, x); /* rotate eigenshapes */
	}
      }
    }

    FEM_Init_Reduced_Order (bod);
  }

  return out;
}

/* export MBFCP definition */
void BODY_2_MBFCP (BODY *bod, FILE *out)
{
  char *kind [] = {"OBSTACLE", "RIGID", "PSEUDO_RIGID", "FINITE_ELEMENT"};
  FORCE *f;
  int n;

  fprintf (out, "ID:\t%d\n", bod->id);
  fprintf (out, "LABEL:\t%s\n", bod->label);
  fprintf (out, "KINEMATICS:\t%s\n", kind [bod->kind]);
  fprintf (out, "BULK_MATERIAL:\t%s\n", bod->mat->label);

  SHAPE_2_MBFCP (bod->shape, out);

  fprintf (out, "VELOCITY:");
  for (n = 0; n < bod->dofs; n ++) fprintf (out, "  %g", bod->velo [n]);
  fprintf (out, "\n");

  for (f = bod->forces, n = 0; f; f = f->next, n ++);
  fprintf (out, "FORCES:\t%d\n", n);
  for (f = bod->forces; f; f = f->next)
  {
    fprintf (out, "KIND:\t");
    if (f->kind & SPATIAL) fprintf (out, "SPATIAL");
    else fprintf (out, "CONVECTED");
    if (f->kind & TORQUE) fprintf (out, "_TORQUE\n");
    else fprintf (out, "\n");

    fprintf (out, "POINT:\t%g  %g  %g\n", f->ref_point [0], f->ref_point [1], f->ref_point [2]);
    fprintf (out, "DIRECTION:\t%g  %g  %g\n", f->direction [0], f->direction [1], f->direction [2]);

    TMS_2_MBFCP (f->data, out);
  }

  fprintf (out, "\n");
}

