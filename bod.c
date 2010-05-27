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
#include "pck.h"
#include "err.h"
#include "lng.h"
#include "fem.h"
#include "but.h"
#include "rnd.h"

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
#define RIG_AUXILIARY(bod) ((bod)->conf+12)
#define RIG_ANGVEL(bod) (bod)->velo
#define RIG_LINVEL(bod) ((bod)->velo+3)
#define RIG_ANGVEL0(bod) ((bod)->velo+6)
#define RIG_LINVEL0(bod) ((bod)->velo+9)
#define RIG_FEXT(bod) ((bod)->velo+12)
#define PRB_GRADIENT(bod) (bod)->conf
#define PRB_CENTER(bod) ((bod)->conf+9)
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
  MX *M;

  if (!bod->inverse)
  {
    int p [] = {0, 9, 18},
	i [] = {0, 3, 6};

    M = MX_Create (MXBD, 6, 2, p, i);
    bod->inverse = M;
  }
  else M = bod->inverse;

  J = BOD_J0 (bod);
  m = BOD_M0 (bod);
  x = M->x;
  NNCOPY (J, x);
  x += 9;
  IDENTITY (x);
  SCALEDIAG (x, m);
  MX_Inverse (M, M);
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
      v = TMS_Value (frc->data, t);
      COPY (frc->direction, f);
      kind = frc->kind;
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
inline static void rig_constraints_force_accum (BODY *bod, double *point, double *base, double *R, short subr, double *force)
{
  double H [18], r [6];

  rig_operator_H (bod, point, base, H);
  blas_dgemv ('T', 3, 6, 1.0, H, 3, R, 1, 0.0, r, 1);

  if (subr)
  {
    force [0] -= r[0];
    force [1] -= r[1];
    force [2] -= r[2];
    force [3] -= r[3];
    force [4] -= r[4];
    force [5] -= r[5];
  }
  else
  {
    force [0] += r[0];
    force [1] += r[1];
    force [2] += r[2];
    force [3] += r[3];
    force [4] += r[4];
    force [5] += r[5];
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
    short tpc = TWO_POINT_CONSTRAINT (con);
    double *point = (isma ? con->mpnt : con->spnt);

    rig_constraints_force_accum (bod, point, con->base, con->R, isma && tpc, force);
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
  MX *M;

  if (!bod->inverse)
  {
    int p [] = {0, 9, 18, 27, 36},
        i [] = {0, 3, 6, 9, 12};

    M = MX_Create (MXBD, 12, 4, p, i);
    bod->inverse = M;
  }
  else M = bod->inverse;

  E0 = bod->ref_tensor; /* referential Euler tensor */
  m = bod->ref_mass;
  x = M->x;
  NNCOPY (E0, x); x += 9;
  NNCOPY (E0, x); x += 9;
  NNCOPY (E0, x); x += 9;
  IDENTITY (x);
  SCALEDIAG (x, m);
  MX_Inverse (M, M);
}

/* set up inverse operator for the implicit dynamic time stepping */
static void prb_dynamic_implicit_inverse (BODY *bod, double step, double *conf, double *force)
{
  MX_DENSE (K, 9, 9);
  MX *M, *A;

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
  else M = bod->M, A = bod->inverse;

    /* calculate stiffness matrix */
  SVK_Tangent_R (lambda (bod->mat->young, bod->mat->poisson),
    mi (bod->mat->young, bod->mat->poisson),
    bod->ref_volume, 9, conf ? conf : bod->conf, K.x);

  if (force)
  {
    /* account for the previous velocity */
    MX_Matvec (1.0 / step, M, bod->velo, 1.0, force);

    /* account for the internal force increment */
    MX_Matvec (-0.25 * step, &K, bod->velo, 1.0, force);
  }

  /* calculate tangent operator A = M + h*h/4 K */
  MX_Add (1.0, MX_Diag(M, 0, 2), 0.25*step*step, &K, MX_Diag(A, 0, 0));
  MX_Copy (MX_Diag (M, 3, 3), MX_Diag (A, 1, 1));

  /* invert A */
  MX_Inverse (A, A);
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
      value = TMS_Value (frc->data, time);
      COPY (frc->direction, f);
      SUB (frc->ref_point, X0, A);
      SCALE (f, value);

      if (frc->kind & CONVECTED) /* obtain spatial force */
      { NVMUL (F, f, P);
	COPY (P, f); }

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

/* solve implicit ingetration nonlineare problem */
static void prb_dynamic_implicit_solve (BODY *bod, double time, double step, double *fext, double *f, short begin)
{
  int iter, imax = 16;

  double half = 0.5 * step,
	*fint = PRB_FINT(bod),
	*velo = PRB_GRADVEL(bod),
	*vel0 = PRB_GRADVEL0(bod),
	 qorig [12],
	 aux [12],
	 res [12],
	 error;

  if (begin)
  {
    blas_dcopy (12, bod->conf, 1, qorig, 1);
    blas_daxpy (12, half, velo, 1, bod->conf, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
    prb_dynamic_force (bod, time+half, step, fext, fint, f);  /* f = fext (t+h/2) - fint (q(t+h/2)) */
    prb_dynamic_implicit_inverse (bod, step, NULL, NULL); /* A = M + (h*h/4) * K */
    MX_Matvec (step, bod->inverse, f, 1.0, velo); /* u(t+h) = u(t) + inv (A) * h * f */
  }
  else
  {
    blas_dcopy (12, bod->conf, 1, qorig, 1);
    blas_daxpy (12, -half, vel0, 1, qorig, 1); /* q(t) = q(t+h/2) - (h/2) * u(t) */
    MX_Matvec (step, bod->inverse, f, 1.0, velo); /* u(t+h) += h * inv (M) * f */
  }
  
  iter = 0;
  do
  {
    NNADD (velo, vel0, aux);
    SCALE9 (aux, 0.25 * step);
    NNADD (qorig, aux, aux);
    prb_internal_force (bod, aux, fint);
    prb_dynamic_implicit_inverse (bod, step, aux, NULL);
    COPY12 (fext, res);
    NNSUB (res, fint, res);
    SCALE12 (res, step);
    SUB12 (velo, vel0, aux);
    MX_Matvec (-1.0, bod->M, aux, 1.0, res);
    MX_Matvec (1.0, bod->inverse, res, 0.0, aux);
    ADD12 (velo, aux, velo);
    error = DOT12 (velo, velo);
    error = sqrt (DOT (aux, aux) / MAX (error, 1.0));
  }
  while (error > 1E-8 && ++ iter < imax);

#if 0
  printf ("DEF_IMP %c: iter = %d, error = %e\n", begin ? 'B' : 'E', iter, error);
#endif

  ASSERT (iter < imax, ERR_BOD_SCHEME_NOT_CONVERGED);

  if (!begin) blas_daxpy (12, half, velo, 1, bod->conf, 1); /* q(t+h) = q(t+h/2) + (h/2) * u(t+h) */
}

/* the smame computation for the static case */
#define prb_static_force(bod, time, step, fext, fint, force) prb_dynamic_force (bod,time,step,fext,fint,force)

/* accumulate constraints reaction */
inline static void prb_constraints_force_accum (BODY *bod, double *point, double *base, double *R, short subr, double *force)
{
  double H [36], r [12];

  prb_operator_H (bod, point, base, H);
  blas_dgemv ('T', 3, 12, 1.0, H, 3, R, 1, 0.0, r, 1);

  if (subr)
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
  else
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
    short tpc = TWO_POINT_CONSTRAINT (con);
    double *point = (isma ? con->mpnt : con->spnt);

    prb_constraints_force_accum (bod, point, con->base, con->R, isma && tpc, force);
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
  double DU [3], *R, *energy = bod->energy;
  DIAB *dia;
  SET *item;
  CON *con;

  for (item = SET_First (bod->con); item; item = SET_Next (item))
  {
    con = item->data;
    if (con->kind == CONTACT)
    {
      dia = con->dia;
#if MPI
      if (dia) /* NULL for external constraints; FIXME: figure out a way of not skipping them */
      {
#endif
      R = con->R;
      ADD (dia->U, dia->V, DU);
      energy [FRICWORK] += 0.5 * step * DOT2 (DU, R);
      energy [CONTWORK] += 0.5 * step * DU [2] * R [2];
#if MPI
      }
#endif
    }
  }
}

/* -------------- interface ------------- */

BODY* BODY_Create (short kind, SHAPE *shp, BULK_MATERIAL *mat, char *label, short form, MESH *msh)
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
      SHAPE_Char (shp,
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
      SHAPE_Char (shp,
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
      bod->form = form;
      bod->msh = msh;
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

  /* update shape adjacency */
  SHAPE_Update_Adjacency (shp);

  /* create piars table */
  bod->sgp = SGP_Create (shp, &bod->nsgp);

  /* update body extents */
  SHAPE_Extents (shp, bod->extents);

  /* not in a list */
  bod->next =
  bod->prev = NULL;

  /* set label */
  bod->label = copylabel (label);

  bod->flags = 0; /* no flags here */

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
  case FEM: return bod->dofs;
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

void BODY_Apply_Force (BODY *bod, short kind, double *point, double *direction, TMS *data, void *call, FORCE_FUNC func)
{
  FORCE *frc;

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

  ERRMEM (frc = MEM_CALLOC (sizeof (FORCE)));

  /* set up force */
  frc->kind = kind;
  if (direction) { NORMALIZE (direction); COPY (direction, frc->direction); }
  if (point) { COPY (point, frc->ref_point); }
  frc->data = data;
  frc->call = call;
  frc->func = func;
  
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

      for (; sph; sph = sph->next)
	if (sph->volume == volume)
	  sph->mat = mat;

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

  return step;
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
	     c = bod->damping,
	     force [12],
	     *L  = PRB_GRADVEL(bod),
	     *L0 = PRB_GRADVEL0(bod),
	     *v  = PRB_LINVEL(bod),
	     *v0 = PRB_LINVEL0(bod);

      NNCOPY (L, L0);
      COPY (v, v0);

      switch (bod->scheme)
      {
      case SCH_DEF_EXP:
      {
	blas_daxpy (12, half, bod->velo, 1, bod->conf, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
	prb_dynamic_force (bod, time+half, step, PRB_FEXT(bod), PRB_FINT(bod), force);  /* f(t+h/2) = fext (t+h/2) - fint (q(t+h/2)) */
	MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) = u(t) + inv (M) * h * f(t+h/2) */
	if (c > 0.0) for (v = L+9; L < v; L ++, L0++) (*L) -= c * (*L0); /* u(t+h) -= c * u (t) (deforomable part) */
      }
      break;
      case SCH_DEF_LIM:
      {
	blas_daxpy (12, half, bod->velo, 1, bod->conf, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
	prb_dynamic_force (bod, time+half, step, PRB_FEXT(bod), PRB_FINT(bod), force);  /* f = fext (t+h/2) - fint (q(t+h/2)) */
	prb_dynamic_implicit_inverse (bod, step, NULL, NULL); /* A = M + (h*h/4) * K */
	MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) = u (t) + inv (A) * h * f */
      }
      break;
      case SCH_DEF_LIM2:
      {
	prb_dynamic_force (bod, time+half, step, PRB_FEXT(bod), PRB_FINT(bod), force);  /* f = fext (t+h/2) - fint (q(t)) */
	prb_dynamic_implicit_inverse (bod, step, NULL, force); /* f += (1/h) M u(t) - (h/4) K u (t) */
	blas_daxpy (12, half, bod->velo, 1, bod->conf, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
	MX_Matvec (step, bod->inverse, force, 0.0, bod->velo); /* u(t+h) = inv (A) * h * force */
      }
      break;
      case SCH_DEF_IMP:
      {
	/* q(t+h/2) = q(t) + (h/2) * u(t)
	 * f = fext (t+h/2) - fint ([q(t) + q(t+h)]/2) 
	 * A = M + (h*h/4) * K ([q(t) + q(t+h)]/2) 
	 * u (t+h) = u (t) + inv (A) * h * f */
        prb_dynamic_implicit_solve (bod, time, step, PRB_FEXT(bod), force, 1);
      }
      default:
      break;
      }
    }
    break;
    case FEM:
      FEM_Dynamic_Step_Begin (bod, time, step);
    break;
  }

  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point);
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
      default: /* SCH_DEF_EXP, SCH_DEF_LIM, SCH_DEF_LIM2 */
      {
        MX_Matvec (step, bod->inverse, force, 1.0, velo); /* u(t+h) += h * inv (M) * force */
        blas_daxpy (12, half, velo, 1, bod->conf, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
      }
      break;
      case SCH_DEF_IMP:
      {
	/* f = fext (t+h/2) - fint ([q(t) + q(t+h)]/2) 
	 * A = M + (h*h/4) * K ([q(t) + q(t+h)]/2) 
	 * u (t+h) = u (t) + inv (A) * h * f
	 * q(t+h) = q(t+h/2) + (h/2) * u(t+h) */
        prb_dynamic_implicit_solve (bod, time, step, fext, force, 0);
      }
      break;
      }

      /* energy */
      ADD12 (velo, vel0, dq);
      SCALE12 (dq, half); /* dq = (h/2) * {u(t) + u(t+h)} */
      energy [EXTERNAL] += DOT12 (dq, fext);
      energy [INTERNAL] += DOT9 (dq, fint);
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

    SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point);
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

  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point);
}

void BODY_Update_Extents (BODY *bod)
{
  double *e = bod->extents;

  SHAPE_Extents (bod->shape, e);

#if MPI
  double *p;
  SET *item;
  CON *con;

  /* now make sure that all attached constraint points are within the extents;
   * due to roundof this could be compromised sometimes, so that constraints might
   * migrate to partitions where bodies would have no representation (child/parent) */

  for (item = SET_First (bod->con); item; item = SET_Next (item))
  {
    con = item->data;
    p = con->point;
    if (p [0] < e [0]) e [0] = p [0] - GEOMETRIC_EPSILON;
    if (p [1] < e [1]) e [1] = p [1] - GEOMETRIC_EPSILON;
    if (p [2] < e [2]) e [2] = p [2] - GEOMETRIC_EPSILON;
    if (p [0] > e [3]) e [3] = p [0] + GEOMETRIC_EPSILON;
    if (p [1] > e [4]) e [4] = p [1] + GEOMETRIC_EPSILON;
    if (p [2] > e [5]) e [5] = p [2] + GEOMETRIC_EPSILON;
  }
#endif
}

void BODY_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x)
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
      FEM_Cur_Point (bod, shp, gobj, X, x);
    break;
  }
}

void BODY_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X)
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
      FEM_Ref_Point (bod, shp, gobj, x, X);
    break;
  }
}

void BODY_Local_Velo (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base, double *prevel, double *curvel)
{
  switch (bod->kind)
  {
    case OBS:
      SET (prevel, 0.0);
      SET (curvel, 0.0);
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
      FEM_Local_Velo (bod, shp, gobj, point, base, prevel, curvel);
    break;
  }
}

MX* BODY_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base)
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
      H = FEM_Gen_To_Loc_Operator (bod, shp, gobj, point, base);
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
    case VALUE_DISPLACEMENT:
    {
      double cur_point [3];

      BODY_Cur_Point (bod, NULL, NULL, point, cur_point);
      SUB (cur_point, point, values);
    }
    break;
    case VALUE_VELOCITY:
    {
      double base [9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

      BODY_Local_Velo (bod, NULL, NULL, point, base, NULL, values);
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
    FEM_Point_Values (bod, point, kind, values);
  break;
  }
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

  if (bod->dom->solfec->ioparallel)
  {
    PBF_Int (bf, &bod->rank, 1);
  }

  if (bod->shape) SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point); 

  if (bod->msh) FEM_Update_Rough_Mesh (bod);
}

void BODY_Pack_State (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_doubles (dsize, d, doubles, bod->conf, BODY_Conf_Size (bod));
  pack_doubles (dsize, d, doubles, bod->velo, bod->dofs);
  pack_doubles (dsize, d, doubles, bod->energy, BODY_ENERGY_SIZE(bod));

#if MPI
  pack_int (isize, i, ints, bod->dom->rank);
#endif
}

void BODY_Unpack_State (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  unpack_doubles (dpos, d, doubles, bod->conf, BODY_Conf_Size (bod));
  unpack_doubles (dpos, d, doubles, bod->velo, bod->dofs);
  unpack_doubles (dpos, d, doubles, bod->energy, BODY_ENERGY_SIZE(bod));

  if (bod->dom->solfec->ioparallel)
  {
    bod->rank = unpack_int (ipos, i, ints);
  }

  if (bod->shape) SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point); 

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

  if (bod->inverse) MX_Destroy (bod->inverse);

  if (bod->M) MX_Destroy (bod->M);

  if (bod->kind == FEM) FEM_Destroy (bod);

  if (bod->msh) MESH_Destroy (bod->msh);

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

  /* constraints: pack their integer ids */
  pack_int (isize, i, ints, SET_Size (bod->con));
  for (SET *item = SET_First (bod->con); item; item = SET_Next (item))
    pack_int (isize, i, ints, CON(item->data)->id);

  /* pack the list of forces */
  pack_forces (bod->forces, dsize, d, doubles, isize, i, ints);

  /* pack scheme and flags */
  pack_int (isize, i, ints, bod->scheme);
  pack_int (isize, i, ints, bod->flags & BODY_PERMANENT_FLAGS);

  /* damping */
  pack_double (dsize, d, doubles, bod->damping);
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
  int ncon, n, id;
  short form;
  DOM *dom;

  /* unpack BODY_Create arguments and create body */
  kind = unpack_int (ipos, i, ints);
  if (kind == FEM) form = unpack_int (ipos, i, ints); else form = 0;
  if (form < 0) msh = MESH_Unpack (sol, dpos, d, doubles, ipos, i, ints), form = -form; else msh = NULL;
  shp = SHAPE_Unpack (sol, dpos, d, doubles, ipos, i, ints);
  label = unpack_string (ipos, i, ints);
  ASSERT_DEBUG_EXT (mat = MATSET_Find (sol->mat, label), "Invalid bulk material label");
  free (label);
  label = unpack_string (ipos, i, ints);
  bod = BODY_Create (kind, shp, mat, label, form, msh);
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

  /* unpack constraints */
  dom = sol->dom;
  ncon = unpack_int (ipos, i, ints);
  for (n = 0; n < ncon; n ++)
  {
    CON *con;

    id = unpack_int (ipos, i, ints);
    ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*) (long) id, NULL), "Invalid constraint id");
    SET_Insert (&dom->setmem, &bod->con, con, NULL);
  }

  /* unpack the list of forces */
  bod->forces = unpack_forces (dpos, d, doubles, ipos, i, ints);

  /* unpack scheme and flags */
  bod->scheme = unpack_int (ipos, i, ints);
  bod->flags = unpack_int (ipos, i, ints);

  /* damping */
  bod->damping = unpack_double (dpos, d, doubles);

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
  short dynamic;
  int j, k;

  dynamic = bod->dom->dynamic;

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
}

/* pack child body */
void BODY_Child_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_doubles (dsize, d, doubles, bod->conf, BODY_Conf_Size (bod));
  pack_doubles (dsize, d, doubles, bod->velo, 2 * bod->dofs); /* current and previous velocity */
  pack_int (isize, i, ints, bod->rank); /* pack parent rank */

  /* pack children ranks */
  pack_int (isize, i, ints, SET_Size (bod->children));
  for (SET *item = SET_First (bod->children); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int) (long) item->data);

  /* pack integration scheme */
  pack_int (isize, i, ints, bod->scheme);
}

/* unpack child body */
void BODY_Child_Unpack (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  short dynamic;
  int j, k, l;

  dynamic = bod->dom->dynamic;

  unpack_doubles (dpos, d, doubles, bod->conf, BODY_Conf_Size (bod));
  unpack_doubles (dpos, d, doubles, bod->velo, 2 * bod->dofs); /* current and previous velocity */
  bod->rank = unpack_int (ipos, i, ints); /* unpack parent rank */

  /* unpack other children ranks */
  SET_Free (&bod->dom->setmem, &bod->children);
  k = unpack_int (ipos, i, ints);
  for (j = 0; j < k; j ++)
  {
    l = unpack_int (ipos, i, ints);
    if (l != bod->dom->rank) /* ommit own rank */
      SET_Insert (&bod->dom->setmem, &bod->children, (void*) (long) l, NULL);
  }

  /* unpack integration scheme */
  bod->scheme = unpack_int (ipos, i, ints);

  /* init inverse */
  if (dynamic) BODY_Dynamic_Init (bod);
  else BODY_Static_Init (bod);

  /* update shape */
  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point);
}

/* pack child update */
void BODY_Child_Update_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_doubles (dsize, d, doubles, bod->conf, BODY_Conf_Size (bod));
}

/* unpack child update */
void BODY_Child_Update_Unpack (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  unpack_doubles (dpos, d, doubles, bod->conf, BODY_Conf_Size (bod));

  /* update shape */
  SHAPE_Update (bod->shape, bod, (MOTION)BODY_Cur_Point); 
}
#endif
