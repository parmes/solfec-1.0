/*
 * epr.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * extended pseudo-rigid model
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

#include "epr.h"
#include "dom.h"
#include "hyb.h"
#include "goc.h"
#include "alg.h"
#include "bla.h"
#include "svk.h"
#include "msh.h"
#include "sph.h"
#include "cvx.h"
#include "err.h"
#include "pck.h"
#include "but.h"

#define EPR_M0(bod) (bod)->ref_mass
#define EPR_V0(bod) (bod)->ref_volume
#define EPR_A0(bod) (bod)->ref_center
#define EPR_J0(bod) (bod)->ref_tensor
#define EPR_CONF_SIZE(bod) ((bod)->dofs + 6)
#define EPR_U(bod) (bod)->conf
#define EPR_CENTER(bod) ((bod)->conf + (bod)->dofs - 6)
#define EPR_ROTATION(bod) ((bod)->conf + (bod)->dofs - 3) 
#define EPR_UVEL(bod) (bod)->velo
#define EPR_LINVEL(bod) ((bod)->velo + (bod)->dofs - 6)
#define EPR_ANGVEL(bod) ((bod)->velo + (bod)->dofs - 3) 
#define EPR_VEL0(bod) ((bod)->velo + (bod)->dofs)
#define EPR_FORCE(bod) ((bod)->velo + (bod)->dofs * 2)
#define EPR_REFTORQ(bod) ((bod)->velo + (bod)->dofs * 3 - 3)

typedef struct epr_element EPR_ELEMENT;
typedef struct epr_mesh EPR_MESH;

/* extended pseudo-rigid element */
struct epr_element
{
  SHAPE *head, *tail; /* beginning and end of shapes within this element */

  EPR_ELEMENT **adj; /* self-inclusive adjacency */

  BULK_MATERIAL *mat;

  int nadj; /* >= 1 */

  double A [3]; /* referential auxiliary point */

  double vol; /* element volume */
};

/* extended pseudo-rigid mesh */
struct epr_mesh
{
  EPR_ELEMENT *ele;

  int nele;

  double INVJ [9], /* global rigid inverse inertia operator (for slow rotations integration) */
	 RHS [3]; /* exp [-(h/2)W(t)]JW(t) + hT(t+h/2) used by rigid integrator */
};

/* overlap callback for element adjacency */
static void* overlap (void *data, BOX *one, BOX *two)
{
  EPR_ELEMENT *epr = one->sgp->shp->epr, *epq = two->sgp->shp->epr, **x, **y;
  double p [3], q [3];

  if (epr == epq) return NULL; /* exclude self-contact */

  if (gobjdistance (GOBJ_Pair_Code (one, two), one->sgp, two->sgp, p, q) < GEOMETRIC_EPSILON) /* if they touch */
  {
    for (x = epr->adj, y = x + epr->nadj; x < y; x ++)
    {
      if (*x == epq) return NULL; /* already adjacent */
    }

    ERRMEM (epr->adj = realloc (epr->adj, (++epr->nadj) * sizeof (EPR_ELEMENT*)));  /* extend adjacency */
    epr->adj [epr->nadj-1] = epq;
    ERRMEM (epq->adj = realloc (epq->adj, (++epq->nadj) * sizeof (EPR_ELEMENT*))); 
    epq->adj [epq->nadj-1] = epr;
  }

  return NULL;
}

/* set up adjacency pointers */
void set_up_adjacency (SHAPE *shp, EPR_MESH *msh)
{
  BOX_Extents_Update update;
  EPR_ELEMENT *ele, *end;
  SGP *sgp, *sge, *sgx;
  BOX **boxes;
  MEM mem;
  int num;

  sgp = SGP_Create (shp, &num);
  MEM_Init (&mem, sizeof (BOX), num);
  ERRMEM (boxes = malloc (sizeof (AABB*) * num));
  for (sgx = sgp, sge = sgp + num, num = 0; sgx < sge; sgx ++, num ++)
  {
    ERRMEM (boxes [num] = MEM_Alloc (&mem));
    update = SGP_Extents_Update (sgx);
    update (sgx->shp->data, sgx->gobj, boxes [num]->extents); /* set up extents */
    boxes [num]->sgp = sgx;
  }

  hybrid (boxes, num, NULL, overlap); /* detect boxoverlaps => set adjacency inside the callback */

  /* append self-reference to the adjacency list */
  for (ele = msh->ele, end = ele + msh->nele; ele < end; ele ++)
  {
    ERRMEM (ele->adj = realloc (ele->adj, (++ele->nadj) * sizeof (EPR_ELEMENT)));
    ele->adj [ele->nadj-1] = ele;
  }

  MEM_Release (&mem); /* done */
  free (boxes);
}

/* dive into shape primitives and try setting up ele[]->mat */
void set_up_bulk_materials (EPR_ELEMENT *ele, int nele)
{
  EPR_ELEMENT *end;
  SHAPE *shp;

  for (end = ele + nele; ele < end; ele ++)
  {
    for (shp = ele->head; shp != ele->tail; shp = shp->next)
    {
      if ((ele->mat = SHAPE_First_Bulk_Material (shp))) break; /* set element material to first encountered primitive shape material */
    }
  }
}

/* compute ele[]->A, ele[]->vol and return a table of static
 * momenta and Euler tensors: {[Sx,Sy,Sz,E(3,3)], [-||-], ...} */
static double* compute_epr_chars (EPR_ELEMENT *ele, int nele)
{
  double *A, *chars, *sta, *eul, coef, sum;
  EPR_ELEMENT *cur, *end, **nei, **nee;
  SHAPE *shp;

  ERRMEM (chars = calloc (nele, sizeof (double [12])));

  /* compute characteristics and volumes */
  for (cur = ele, end = cur + nele, sta = chars, eul = sta + 3; cur < end; cur ++, sta += 12, eul += 12)
  {
    for (shp = cur->head; shp != cur->tail; shp = shp->next)
    {
      SHAPE_Char_Partial (shp, &cur->vol, &sta[0], &sta[1], &sta[2], eul);
    }
  }

  /* compute ele[]->A */
  for (cur = ele, end = cur + nele, A = cur->A; cur < end; cur ++, A = cur->A)
  {
    SET (A, 0.0);
    sum = 0.0;

    for (nei = cur->adj, nee = nei + cur->nadj; nei < nee; nei ++)
    {
      sta = &chars [12 * ((*nei) - ele)];
      coef = 1.0 / (double) (*nei)->nadj;
      ADDMUL (A, coef, sta, A);
      sum += coef * (*nei)->vol;
    }
   
    DIV (A, sum, A); 
  }

  return chars;
}

/* sum up a single value into a column map of values */
static void sum_up (int i, int j, double x, MEM *map, MEM *val, MAP **col)
{
  double *y;
  MAP *item;

  if (!(item = MAP_Find_Node (col [j], (void*) (long) i, NULL)))
  {
    ERRMEM (y = MEM_Alloc (val));
    item = MAP_Insert (map, &col [j], (void*) (long) i, y, NULL);
  }

  y = item->data, *y += x;
}

/* sum up element inertia into a column map of values */
static void sum_up_inertia (int i, int j, double *Ai, double *Aj, double *eul, double *sta, int nadj, double density, MEM *map, MEM *val, MAP **col)
{
  double XAj [9], AiX[9], AiAj [9], MijT [9], coef;
  double I, J, k;

  DIADIC (sta, Aj, XAj);
  DIADIC (Ai, sta, AiX);
  DIADIC (Ai, Aj, AiAj);

  NTSUB (eul, XAj, MijT);
  NTSUB (MijT, AiX, MijT);
  NTADD (MijT, AiAj, MijT);

  coef = density / (double) nadj * nadj;

  SCALE9 (MijT, coef);

  i *= 9; /* from element numbering to matrix block numbering */
  j *= 9;

  for (k = 0; k < 9; k += 3) /* three diagonal blocks */
  {
    I = i + k;
    J = j + k;

    sum_up (I+0, J+0, MijT [0], map, val, col);
    sum_up (I+1, J+0, MijT [1], map, val, col);
    sum_up (I+2, J+0, MijT [2], map, val, col);
    sum_up (I+0, J+1, MijT [3], map, val, col);
    sum_up (I+1, J+1, MijT [4], map, val, col);
    sum_up (I+2, J+1, MijT [5], map, val, col);
    sum_up (I+0, J+2, MijT [6], map, val, col);
    sum_up (I+1, J+2, MijT [7], map, val, col);
    sum_up (I+2, J+2, MijT [8], map, val, col);
  }
}

/* get shape stabbed by a referential point */
static SHAPE* epr_shape_stab (BODY *bod, double *point)
{
  SHAPE *shp;

  if (SHAPE_Gobj (bod->shape, point, &shp)) return shp; /* TODO: optimize by an EPR-specific spatial tree based search */
  else return NULL;
}

/* get (self-inclusive) average deformation gradient over adjacency */
static void adj_average_gradient (EPR_ELEMENT *base, EPR_ELEMENT *ele, double *R, double *conf, double F[9])
{
  EPR_ELEMENT **nei, **nee;
  double U[9], *Ui, coef;
  int shift;

  coef = 1.0 / (double) ele->nadj;

  IDENTITY (U);

  for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
  {
    shift =  9 * ((*nei) - base);
    Ui = &conf [shift];
    NNADD (U, Ui, U);
  }

  SCALE9 (U, coef);

  NNMUL (R, U, F);
}

/* get lame coefficients */
static void lame_coefs (BODY *bod, EPR_ELEMENT *ele, double *lame_lambda, double *lame_mi)
{
  if (ele->mat)
  {
    *lame_lambda = lambda (ele->mat->young, ele->mat->poisson);
    *lame_mi = mi (ele->mat->young, ele->mat->poisson);
  }
  else
  {
    *lame_lambda = lambda (bod->mat->young, bod->mat->poisson);
    *lame_mi = mi (bod->mat->young, bod->mat->poisson);
  }
}

/* compute force(time) = external(time) - internal(time), provided
 * that the configuration(time) has already been set elsewhere */
static void epr_dynamic_force (BODY *bod, double time, double step, double *force)
{
  EPR_ELEMENT *ele, *end, **nei, **nee;
  EPR_MESH *msh;
  SHAPE *shp;
  double f [3],
	 g [3],
	*point,
	*cur,
	*last,
	*lin,
	 F [9],
	 P [9],
	 value,
	 B [3],
	 *R,
	 *A0,
	 *spatorq,
	 b [3],
	 coef;
  int shift,
      kind;

  msh = bod->priv;

  blas_dscal (bod->dofs, 0.0, force, 1); /* zero force */

  lin = force + bod->dofs - 6;

  spatorq = lin + 3;

  R = EPR_ROTATION (bod); /* current rigid rotation */
  A0 = EPR_A0 (bod); /* referential mass center */

  /* point forces */
  for (FORCE *frc = bod->forces; frc; frc = frc->next)
  {
    if (frc->func)
    {
      ERRMEM (cur = calloc (bod->dofs, sizeof (double)));
      frc->func (frc->data, frc->call, bod->dofs, bod->conf, bod->dofs, bod->velo, time, step, cur);
      blas_daxpy (bod->dofs, 1.0, cur, 1, force, 1);
      free (cur);
    }
    else
    {
      value = TMS_Value (frc->data, time);
      COPY (frc->direction, f);
      SCALE (f, value);
      point = frc->ref_point;
      kind = frc->kind;

      if (kind & CONVECTED) /* obtain spatial force */
      {
	NVMUL (R, f, g);
	COPY (g, f);
      }

      shp = epr_shape_stab (bod, point); /* get element stabbed by the point */

      if (shp && !(kind & TORQUE)) /* deformable term */
      {
	ele = shp->epr;

	for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
	{
	  shift = 9 * ((*nei) - msh->ele);
	  cur = &force [shift];
	  coef = 1.0 / (double) (*nei)->nadj;

	  SUB (point, (*nei)->A, B);

	  cur [0] += B [0] * f [0] * coef;
	  cur [1] += B [1] * f [0] * coef;
	  cur [2] += B [2] * f [0] * coef;
	  cur [3] += B [0] * f [1] * coef;
	  cur [4] += B [1] * f [1] * coef;
	  cur [5] += B [2] * f [1] * coef;
	  cur [6] += B [0] * f [2] * coef;
	  cur [7] += B [1] * f [2] * coef;
	  cur [8] += B [2] * f [2] * coef;
	}
      }

      if (kind & CONVECTED) /* rigid terms */
      { 
	if (kind & TORQUE)
	{ 
          NVADDMUL (spatorq, R, f, spatorq); /* add spatial image of it */
	}
	else
	{
	  SUB (point, A0, B);          /* (X - A0) => referential force arm */
	  NVMUL (R, B, b);             /* a => spatial force arm */
	  PRODUCTADD (b, f, spatorq);  /* (x-x0) x b => spatial torque */
	  ADD (lin, f, lin);           /* linear force */
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
	  SUB (point, A0, B);          /* (X - X0) => referential force arm */
	  NVMUL (R, B, b);             /* a => spatial force arm */
	  PRODUCTADD (b, f, spatorq);  /* (x-x0) x b => spatial torque */
	  ADD (lin, f, lin);           /* linear force */
	}
      }
    }
  }

  /* gravity */
  if (DOM(bod->dom)->gravval)
  {
    COPY (DOM(bod->dom)->gravdir, f);
    value = TMS_Value (DOM(bod->dom)->gravval, time);
    SCALE (f, value);
    ADDMUL (lin, bod->ref_mass, f, lin);
  }

  /* internal forces */
  for (ele = msh->ele, end = ele + msh->nele; ele < end; ele ++)
  {
    double lambda, mi;

    for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
    {
      lame_coefs (bod, *nei, &lambda, &mi);
      adj_average_gradient (msh->ele, *nei, R, bod->conf, F);
      SVK_Stress_R (lambda, mi, (*nei)->vol, F, P);
      coef = 1.0 / (double) (*nei)->nadj;
      NNSUBMUL (cur, coef, P, cur); /* force = external - internal */
    }
  }

  /* multiply everything but the linear term by R^T
   * (mind that U-corresponding force blocks are row-wise) */

  for (cur = force, last = force + 9 * msh->nele; cur < last; cur += 9)
  {
    TNCOPY (cur, P);
    TNMUL (R, P, F);
    NTCOPY (F, cur);
  }

  COPY (spatorq, f);
  TVMUL (R, f, spatorq);
}

/* static unbalanced force computation */
#define epr_static_force(bod, time, step, force) epr_dynamic_force (bod, time, step, force)

/* compute static inverse operator */
static void epr_static_inverse (BODY *bod, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* accumulate constraints reaction */
inline static void epr_constraints_force_accum (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base, double *R, short isma, double *force)
{
  MX *H;

  H = EPR_Gen_To_Loc_Operator (bod, shp, gobj, point, base);

  if (isma) MX_Matvec (-1.0, MX_Tran (H), R, 1.0, force);
  else MX_Matvec (1.0, MX_Tran (H), R, 1.0, force);

  MX_Destroy (H);
}

/* calculate constraints reaction */
static void epr_constraints_force (BODY *bod, double *force)
{
  SET *node;

  blas_dscal (bod->dofs, 0.0, force, 1); /* zero force */

  for (node = SET_First (bod->con); node; node = SET_Next (node))
  {
    CON *con = node->data;
    short isma = (bod == con->master);
    double *point = (isma ? con->mpnt : con->spnt);
    SHAPE *shp = (isma ? con->mshp : con->sshp);
    void *gobj = (isma ? con->mgobj : con->sgobj);

    epr_constraints_force_accum (bod, shp, gobj, point, con->base, con->R, isma, force);
  }

#if MPI
  for (MAP *node = MAP_First (bod->conext); node; node = MAP_Next (node))
  {
    CONEXT *con = node->data;

    epr_constraints_force_accum (bod, con->sgp->shp, con->sgp->gobj, con->point, con->base, con->R, con->isma, force);
  }
#endif
}

/* calculates a block of transformation operator from generalized
 * velocity space to the local cartesian space at 'X' in given 'base^T * R' */
static void epr_operator_block_H (double *X, double *A, double *baseTR, double *H) /* H is assumed all zero on entry */
{
  double B [3];

  SUB (X, A, B); 

  H [0] = baseTR [0] * B [0];
  H [3] = baseTR [0] * B [1];
  H [6] = baseTR [0] * B [2];
  H [9] = baseTR [3] * B [0];
  H [12] = baseTR [3] * B [1];
  H [15] = baseTR [3] * B [2];
  H [18] = baseTR [6] * B [0];
  H [21] = baseTR [6] * B [1];
  H [24] = baseTR [6] * B [2];

  H [1] = baseTR [1] * B [0];
  H [4] = baseTR [1] * B [1];
  H [7] = baseTR [1] * B [2];
  H [10] = baseTR [4] * B [0];
  H [13] = baseTR [4] * B [1];
  H [16] = baseTR [4] * B [2];
  H [19] = baseTR [7] * B [0];
  H [22] = baseTR [7] * B [1];
  H [25] = baseTR [7] * B [2];

  H [2] = baseTR [2] * B [0];
  H [5] = baseTR [2] * B [1];
  H [8] = baseTR [2] * B [2];
  H [11] = baseTR [5] * B [0];
  H [14] = baseTR [5] * B [1];
  H [17] = baseTR [5] * B [2];
  H [20] = baseTR [8] * B [0];
  H [23] = baseTR [8] * B [1];
  H [26] = baseTR [8] * B [2];
}

/* calculate cauche stress for a pseudo-rigid body */
void epr_cauchy (BODY *bod, EPR_MESH *msh, EPR_ELEMENT *ele, double *stress)
{
  double P [9], J, F [9], lambda, mi;

  lame_coefs (bod, ele, &lambda, &mi);
  adj_average_gradient (msh->ele, ele, EPR_ROTATION (bod), bod->conf, F);

  J = SVK_Stress_R (lambda, mi , 1.0, F, P); /* per unit volume */

  stress [0] = (F [0]*P [0] + F [1]*P [3] + F [2]*P [6]) / J;
  stress [1] = (F [3]*P [1] + F [4]*P [4] + F [5]*P [7]) / J;
  stress [2] = (F [6]*P [2] + F [7]*P [5] + F [8]*P [8]) / J;
  stress [3] = (F [0]*P [1] + F [1]*P [4] + F [2]*P [7]) / J;
  stress [4] = (F [0]*P [2] + F [1]*P [5] + F [2]*P [8]) / J;
  stress [5] = (F [3]*P [2] + F [4]*P [5] + F [5]*P [8]) / J;
}

#if DEBUG
/* test symmetry of a sparse matrix */
static int test_symmetry (MAP **col, int n)
{
  MAP *item, *jtem;
  double xij, xji;
  int i, j;

  for (j = 0; j < n; j ++)
  {
    for (item = MAP_First (col [j]); item; item = MAP_Next (item))
    {
      i = (int) (long) item->key;

      xji = * (double*) item->data;

      if (!(jtem = MAP_Find_Node (col [i], (void*) (long) j, NULL))) return 0;

      xij = * (double*) jtem->data;

      if (xij != xji) return 0;
    }
  }

  return 1;
}
#endif

/* create EPR internals for a body */
void EPR_Create (SHAPE *shp, BULK_MATERIAL *mat, BODY *bod)
{
  EPR_ELEMENT *ele;
  double euler [9];
  EPR_MESH *msh;
  SHAPE *shq;

  /* get some global mass characteristics */
  SHAPE_Char (shp, &bod->ref_volume, bod->ref_center, euler);
  bod->ref_mass = bod->ref_volume * mat->density;
  SCALE9 (bod->ref_tensor, mat->density);
  euler2inertia (euler, bod->ref_tensor);

  /* create mesh */ 
  ERRMEM (msh = malloc (sizeof (EPR_MESH)));

  for (shq = shp, msh->nele = 1; shq; shq = shq->next) /* initialize with 1 for the sake of the last NULL-terimed element */
  {
    if (shq->epr) msh->nele ++; /* extended element list end marker */
  }

  ERRMEM (msh->ele = calloc (msh->nele, sizeof (EPR_ELEMENT)));

  for (shq = shp, ele = msh->ele, ele->head = shq; shq; shq = shq->next)
  {
    shq->epr = ele; /* current shape points to current element */

    if (shq->epr) /* end of element list marker (inclusive) */
    {
      ele->tail = shq->next; /* record list tail (exclusive) */

      ele ++; /* iterate to next element */

      if (shq->next) ele->head = shq->next; /* record list head */
    }
  }

  set_up_adjacency (shp, msh);

  bod->dofs = 9 * msh->nele + 6; /* a number of U_i (fast), angular velocity (slow), linear velocity (fast) */

  ERRMEM (bod->conf = calloc (4 * bod->dofs + 6, sizeof (double))); /* conf, velo, vel0, force */
  bod->velo = bod->conf + EPR_CONF_SIZE (bod);
 
  IDENTITY (EPR_ROTATION (bod)); /* R(0) is 1, U_i(0) are 0 */
  COPY (bod->ref_center, EPR_CENTER (bod)); /* referential mass center */

  set_up_bulk_materials (msh->ele, msh->nele);

  bod->priv = msh;
}

/* return configuration size */
int EPR_Conf_Size (BODY *bod)
{
  return EPR_CONF_SIZE (bod);
}

/* overwrite state */
void EPR_Overwrite_State (BODY *bod, double *q, double *u)
{
  double *x, *y;

  for (x = bod->conf, y = x + EPR_CONF_SIZE (bod); x < y; x ++, q ++) *x = *q;

  for (x = bod->velo, y = x + bod->dofs; x < y; x ++, u ++) *x = *u;
}

/* set initial rigid motion velocity */
void EPR_Initial_Velocity (BODY *bod, double *linear, double *angular)
{
  double *l = EPR_LINVEL (bod),
	 *a = EPR_ANGVEL (bod);

  if (linear) { COPY (linear, l); }

  if (angular) { COPY (angular, a); }
}

/* initialise dynamic time stepping */
void EPR_Dynamic_Init (BODY *bod, SCHEME scheme)
{
  EPR_ELEMENT *ele, *end, **nei, **nee, **nej, **nff;
  double *xx, *chars, *sta, *eul, *j0, den;
  int m, n, *p, *i, ii, jj, *pp, *pi;
  MEM mapmem, valmem;
  MAP **col, *item;
  EPR_MESH *msh;

  msh = bod->priv;

  if (bod->inverse) MX_Destroy (bod->inverse);

  chars = compute_epr_chars (msh->ele, msh->nele);

  den = bod->mat->density;

  m = n = bod->dofs;

  ERRMEM (col = calloc (n, sizeof (MAP*)));
  MEM_Init (&mapmem, sizeof (MAP), MAX (256, 8 * n));
  MEM_Init (&valmem, sizeof (double), MAX (256, 8 * n));

  /* sum up element inertia into the column values map */
  for (ele = msh->ele, end = ele + msh->nele; ele < end; ele ++)
  {
    sta = &chars [12 * (ele - msh->ele)];
    eul = sta + 3;

    for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
    {
      ii = (*nei) - msh->ele;
      for (nej = ele->adj, nff = nej + ele->nadj; nej < nff; nej ++)
      {
	jj = (*nej) - msh->ele;
	sum_up_inertia (ii, jj, (*nei)->A, (*nej)->A, eul, sta, ele->nadj, den, &mapmem, &valmem, col);
      }
    }
  }

  /* sum up scalar mass at the end */
  sum_up (m-6, n-6, EPR_M0(bod), &mapmem, &valmem, col);
  sum_up (m-5, n-5, EPR_M0(bod), &mapmem, &valmem, col);
  sum_up (m-4, n-4, EPR_M0(bod), &mapmem, &valmem, col);

  /* sum up global inertia */
  j0 = EPR_J0 (bod);
  sum_up (m-3, n-3, j0[0], &mapmem, &valmem, col);
  sum_up (m-2, n-3, j0[1], &mapmem, &valmem, col);
  sum_up (m-1, n-3, j0[2], &mapmem, &valmem, col);
  sum_up (m-3, n-2, j0[3], &mapmem, &valmem, col);
  sum_up (m-2, n-2, j0[4], &mapmem, &valmem, col);
  sum_up (m-1, n-2, j0[5], &mapmem, &valmem, col);
  sum_up (m-3, n-1, j0[6], &mapmem, &valmem, col);
  sum_up (m-2, n-1, j0[7], &mapmem, &valmem, col);
  sum_up (m-1, n-1, j0[8], &mapmem, &valmem, col);

  /* test symmetry */
  ASSERT_DEBUG (test_symmetry (col, n), "Symmetry test failed");

  /* compute sparse storage */
  for (ii = jj = 0; jj < n; jj ++)
  {
    ii += MAP_Size (col [jj]);
  }

  ERRMEM (p = malloc (n * sizeof (int)));
  ERRMEM (i = malloc (ii * sizeof (int)));

  p [0] = 0;
  pp = &p[1];

  for (jj = 0, pi = i; jj < n; jj ++, pp ++)
  {
    for (item = MAP_First (col[jj]); item; item = MAP_Next (item), pi ++)
    {
      *pi = (int) (long) item->key;
    }

    *pp = (*(pp-1)) + MAP_Size (col [jj]);
  }

  /* create matrix */
  bod->inverse = MX_Create (MXCSC, m, n, p, i);

  /* set values and scale by volume density */
  for (jj = 0, xx = bod->inverse->x; jj < n; jj ++)
  {
    for (item = MAP_First (col[jj]); item; item = MAP_Next (item), xx ++)
    {
      *xx = *((double*) item->data);
    }
  }

  /* compute inverse factorization */
  MX_Inverse (bod->inverse, bod->inverse);

  /* compute inverse inertia */
  INVERT (EPR_J0 (bod), msh->INVJ, den);
  ASSERT_DEBUG (den != 0.0, "Inverting inertia tensor failed");

  /* set integration scheme */
  ASSERT (scheme == SCH_DEFAULT, ERR_BOD_SCHEME);
  bod->scheme = scheme;

  /* clean up */
  MEM_Release (&mapmem);
  MEM_Release (&valmem);
  free (chars);
  free (col);
  free (p);
  free (i);
}

/* estimate critical step for the dynamic scheme */
double EPR_Dynamic_Critical_Step (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
  return 0.0;
}

/* perform the initial half-step of the dynamic scheme */
void EPR_Dynamic_Step_Begin (BODY *bod, double time, double step)
{
  EPR_MESH *msh = bod->priv;
  double half = 0.5 * step,
        *force = EPR_FORCE (bod),
	*J = EPR_J0(bod),
	*I = msh->INVJ,
	*S = msh->RHS,
	*T = EPR_REFTORQ (bod), /* some shited storage within 'force' */
	*W = EPR_ANGVEL(bod),
	*R = EPR_ROTATION(bod),
	O [9], DR [9], W05 [3];

  blas_dcopy (bod->dofs, bod->velo, 1, EPR_VEL0 (bod), 1); /* save u(t) */
  blas_daxpy (bod->dofs - 3, half, bod->velo, 1, bod->conf, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */

  COPY (W, O);
  SCALE (O, half);
  EXPMAP (O, DR); 
  NNCOPY (R, O);
  NNMUL (O, DR, R);       /* R(t+h/2) = R(t) exp [(h/2) * W(t)] */

  epr_dynamic_force (bod, time+half, step, force);  /* f(t+h/2) = fext (t+h/2) - fint (q(t+h/2)) */
  
  NVMUL (J, W, O);
  TVMUL (DR, O, S);
  ADDMUL (S, step, T, S);  /* exp [-(h/2)W(t)]JW(t) + hT(t+h/2) */

  COPY (T, O);
  SCALE (O, half); 
  NVMUL (J, W, W05);
  TVADDMUL (O, DR, W05, O);
  NVMUL (I, O, W05);            /* W05 = I * [exp [-(h/2) * W(t)]*J*W(t) + (h/2)*T(t+h/2)] */

  NVMUL (J, W05, O);
  PRODUCTSUB (W05, O, T);   /* force [0..2] = T - W05 x J * W05 */

  MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) = u(t) + inv (M) * h * f(t+h/2) */
}

/* perform the final half-step of the dynamic scheme */
void EPR_Dynamic_Step_End (BODY *bod, double time, double step)
{
  EPR_MESH *msh = bod->priv;
  double half = 0.5 * step,
        *force = EPR_FORCE (bod),
        *R = EPR_ROTATION(bod),
        *W = EPR_ANGVEL(bod),
	*I = msh->INVJ,
	*S = msh->RHS,
	*T = EPR_REFTORQ (bod),
         O [9], DR [9];

  epr_constraints_force (bod, force); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) += inv (M) * h * r */
  blas_daxpy (bod->dofs - 3, half, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */

  COPY (W, O);
  SCALE (O, half);
  EXPMAP (O, DR); 
  NNCOPY (R, O);
  NNMUL (O, DR, R); /* R(t) = R(t+h/2) exp [(h/2) * W(t+h)] */

  ADDMUL (S, step, T, S);
  TVMUL (DR, S, O);
  NVMUL (I, O, W); /* W2(t+h) = I * exp [-(h/2)W1(t+h)][exp [-(h/2)W(t)]JW(t) + hT(t+h/2)] */

  /* symmetrize Ui */
  blas_dcopy (bod->dofs - 6, bod->velo, 1, force, 1); /* copy {Ui} int 'force' then copy back and symmetrize */
  for (double *IN = force, *OUT = bod->velo, *END = IN + bod->dofs - 6; IN < END; IN += 9, OUT += 9) { SYMM (IN, OUT); }
}

/* initialise static time stepping */
void EPR_Static_Init (BODY *bod)
{
  /* cancel initial velocity */
  blas_dscal (bod->dofs, 0.0, bod->velo, 1);
}

/* perform the initial half-step of the static scheme */
void EPR_Static_Step_Begin (BODY *bod, double time, double step)
{
  double *force = EPR_FORCE (bod);

  epr_static_inverse (bod, step); /* compute inverse of static tangent operator */
  epr_static_force (bod, time+step, step, force);  /* f(t+h) = fext (t+h) - fint (q(t+h)) */
  MX_Matvec (step, bod->inverse, force, 0.0, bod->velo); /* u(t+h) = inv (A) * h * f(t+h) */
}

/* perform the final half-step of the static scheme */
void EPR_Static_Step_End (BODY *bod, double time, double step)
{
  double *force = EPR_FORCE (bod),
	 *R = EPR_ROTATION(bod),
	 *W = EPR_ANGVEL(bod),
	 O [9], DR [9];
      
  epr_constraints_force (bod, force); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) += inv (A) * h * r */
  blas_daxpy (bod->dofs - 3, step, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t) + h * u(t+h) */

  COPY (W, O);
  SCALE (O, step);
  EXPMAP (O, DR); 
  NNCOPY (R, O);
  NNMUL (O, DR, R); /* R(t) = R(t) exp [h * W(t+h)] */

  /* symmetrize Ui */
  blas_dcopy (bod->dofs - 6, bod->velo, 1, force, 1); /* copy {Ui} int 'force' then copy back and symmetrize */
  for (double *IN = force, *OUT = bod->velo, *END = IN + bod->dofs - 6; IN < END; IN += 9, OUT += 9) { SYMM (IN, OUT); }
}

/* motion x = x (X, state) */
void EPR_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x)
{
  double *R = EPR_ROTATION(bod),
	 *c = EPR_CENTER(bod),
	 *C = EPR_A0(bod),
	 A [3];

  SUB (X, C, A);
  NVADDMUL (c, R, A, x);
}

/* inverse motion X = X (x, state) */
void EPR_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X)
{
  double *R = EPR_ROTATION(bod),
	 *c = EPR_CENTER(bod),
	 *C = EPR_A0(bod),
	 a [3];

  SUB (x, c, a);
  TVADDMUL (C, R, a, X);
}

/* obtain velocity at (element, point), expressed in the local 'base' => all entities are spatial */
void EPR_Local_Velo (BODY *bod, VELOTIME time, SHAPE *shp, void *gobj, double *point, double *base, double *velo)
{
  double *genvel, *W, *U, *a, *A, *R, coef, B [3], v [3], w [3];
  EPR_ELEMENT *ele, **nei, **nee;
  EPR_MESH *msh;
  int shift;

  msh = bod->priv;
  ele = shp->epr;
  genvel = time == CURVELO ? bod->velo : EPR_VEL0 (bod);
  a = genvel + bod->dofs - 6;
  W = a + 3;
  coef = 1.0  / (double) ele->nadj;
  R = EPR_ROTATION (bod);

  SET (w, 0.0);
  for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
  {
    shift = 9 * ((*nei) - msh->ele);
    U = &genvel [shift];
    A = (*nei)->A;

    SUB (point, A, B);
    NVADDMUL (w, U, B, w); /* SUM Ni dUi/dt (point - Ai) */
  }
  SCALE (w, coef); /* scale by shape function value */

  NVADDMUL (a, R, w, v); /* v = R * SUM Ni dUi/dt (point - Ai) + da/dt */

  A = EPR_A0 (bod);
  SUB (point, A, B);
  PRODUCT (W, B, w);
  ADD (v, w, v); /* v = dR/dT (point - A0) + R * SUM Ni dUi/dt (point - Ai) + da/dt */

  TVMUL (base, v, velo); /* express in local base: velo = base^T * v */
}

/* return transformation operator from the generalised to the local velocity space at (element, point, base) */
MX* EPR_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base)
{
  double *tail, coef, *R, baseTR [9], *A0, A [3], S [9];
  EPR_ELEMENT *ele, **nei, **nee;
  EPR_MESH *msh;
  int shift;
  MX *H;

  msh = bod->priv;
  ele = shp->epr;
  coef = 1.0 / (double) ele->nadj;
  R = EPR_ROTATION (bod);
  A0 = EPR_A0 (bod);

  TNMUL (base, R, baseTR);

  H = MX_Create (MXDENSE, 3, bod->dofs, NULL, NULL); /* initially a zero matrix */

  /* deformable terms */
  for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
  {
    shift = 27 * ((*nei) - msh->ele);
    epr_operator_block_H (point, (*nei)->A, baseTR, &H->x[shift]);
  }

  MX_Scale (H, coef); /* scale by shape function value */

  /* linear term */
  tail = &H->x[3 * bod->dofs - 18];
  TNCOPY (base, tail); /* copy last bit from the linear motion */

  /* angular term */
  tail += 9;
  SUB (point, A0, A);
  VECSKEW (A, S);
  NTMUL (baseTR, S, tail);

  return H;
}

/* compute current kinetic energy */
double EPR_Kinetic_Energy (BODY *bod)
{
  double *tmp = EPR_FORCE (bod);

  ASSERT_DEBUG (bod->inverse, "Invalid inertia operator");

  MX_Matvec (1.0, MX_Uninv (bod->inverse), bod->velo, 0.0, tmp);

  return 0.5 * blas_ddot (bod->dofs, tmp, 1, bod->velo, 1);
}

/* get some values at a node of a geometrical object */
void EPR_Nodal_Values (BODY *bod, SHAPE *shp, void *gobj, int node, VALUE_KIND kind, double *values)
{
  double ref_point [3];

  switch (shp->kind)
  {
  case SHAPE_MESH:
  {
    MESH *msh = shp->data;
    ELEMENT *ele = gobj;
    int n = ele->nodes [node];

    COPY (msh->ref_nodes [n], ref_point);
  }
  break;
  case SHAPE_CONVEX:
  {
    CONVEX *cvx = gobj;
    double *ref = cvx->ref + 3 * node;

    COPY (ref, ref_point);
  }
  break;
  case SHAPE_SPHERE:
  {
    SPHERE *sph = gobj;

    COPY (sph->ref_center, ref_point);
  }
  break;
  }

  switch (kind)
  {
  case VALUE_DISPLACEMENT:
  {
    double cur_point [3];

    EPR_Cur_Point (bod, shp, gobj, ref_point, cur_point);
    SUB (cur_point, ref_point, values);
  }
  break;
  case VALUE_VELOCITY:
  {
    double base [9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    EPR_Local_Velo (bod, CURVELO, shp, gobj, ref_point, base, values);
  }
  break;
  case VALUE_STRESS:
  {
    epr_cauchy (bod, bod->priv, shp->epr, values);
  }
  break;
  case VALUE_MISES:
  {
    double stress [6];

    epr_cauchy (bod, bod->priv, shp->epr, stress);
    MISES (stress, values [0]);
  }
  break;
  case VALUE_STRESS_AND_MISES:
  {
    epr_cauchy (bod, bod->priv, shp->epr, values);
    MISES (values, values [6]);
  }
  break;
  }
}

/* get some values at a referential point */
void EPR_Point_Values (BODY *bod, double *point, VALUE_KIND kind, double *values)
{
  SHAPE *shp;

  if ((shp = epr_shape_stab (bod, point)))
  {
    switch (kind)
    {
    case VALUE_DISPLACEMENT:
    {
      double cur_point [3];

      EPR_Cur_Point (bod, shp, NULL, point, cur_point);
      SUB (cur_point, point, values);
    }
    break;
    case VALUE_VELOCITY:
    {
      double base [9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

      EPR_Local_Velo (bod, CURVELO, shp, NULL, point, base, values);
    }
    break;
    case VALUE_STRESS:
    {
      epr_cauchy (bod, bod->priv, shp->epr, values);
    }
    break;
    case VALUE_MISES:
    {
      double stress [6];

      epr_cauchy (bod, bod->priv, shp->epr, stress);
      MISES (stress, values [0]);
    }
    break;
    case VALUE_STRESS_AND_MISES:
    {
      epr_cauchy (bod, bod->priv, shp->epr, values);
      MISES (values, values [6]);
    }
    break;
    }
  }
}

/* release EPR memory */
void EPR_Destroy (BODY *bod)
{
  EPR_ELEMENT *ele, *end;
  EPR_MESH *msh;

  msh = bod->priv;

  for (ele = msh->ele, end = ele + msh->nele; ele < end; ele ++) free (ele->adj);

  free (msh->ele);

  free (msh);

  bod->priv = NULL;

  if (bod->inverse) MX_Destroy (bod->inverse);

  bod->inverse = NULL;
}

/* write state */
void EPR_Write_State (BODY *bod, PBF *bf)
{
  double *U = bod->conf,
	 *V = bod->velo,
	 *W = U + bod->dofs - 6,
	  Z [12];

  for (; U < W; U += 9, V += 9)
  {
    SYMMLOTR (U, Z);
    SYMMLOTR (V, Z+6);

    PBF_Double (bf, Z, 12);
  }

  PBF_Double (bf, U, 12);
  PBF_Double (bf, V, 6);
}

/* read state */
void EPR_Read_State (BODY *bod, PBF *bf)
{
  double *U = bod->conf,
	 *V = bod->velo,
	 *W = U + bod->dofs - 6,
	  Z [12];

  for (; U < W; U += 9, V += 9)
  {
    PBF_Double (bf, Z, 12);

    LOTRSYMM (Z, U);
    LOTRSYMM (Z+6, V);
  }

  PBF_Double (bf, U, 12);
  PBF_Double (bf, V, 6);
}

/* pack state */
void EPR_Pack_State (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  double *U = bod->conf,
	 *V = bod->velo,
	 *W = U + bod->dofs - 6,
	  Z [12];

  for (; U < W; U += 9, V += 9)
  {
    SYMMLOTR (U, Z);
    SYMMLOTR (V, Z+6);

    pack_doubles (dsize, d, doubles, Z, 12);
  }

  pack_doubles (dsize, d, doubles, U, 12);
  pack_doubles (dsize, d, doubles, V, 6);
}

/* unpack state */
void EPR_Unpack_State (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  double *U = bod->conf,
	 *V = bod->velo,
	 *W = U + bod->dofs - 6,
	  Z [12];

  for (; U < W; U += 9, V += 9)
  {
    unpack_doubles (dpos, d, doubles, Z, 12);

    LOTRSYMM (Z, U);
    LOTRSYMM (Z+6, V);
  }

  unpack_doubles (dpos, d, doubles, U, 12);
  unpack_doubles (dpos, d, doubles, V, 6);
}
