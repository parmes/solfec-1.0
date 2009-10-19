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

#define EPR_VEL0(bod) ((bod)->conf + 2 * (bod)->dofs)
#define EPR_FORCE(bod) ((bod)->conf + 3 * (bod)->dofs)

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
static void sum_up_inertia (int i, int j, double *Ai, double *Aj, double *eul, double *sta, int nadj, MEM *map, MEM *val, MAP **col)
{
  double XAj [9], AiX[9], AiAj [9], MijT [9], coef;
  double I, J, k;

  DIADIC (sta, Aj, XAj);
  DIADIC (Ai, sta, AiX);
  DIADIC (Ai, Aj, AiAj);

  NTSUB (eul, XAj, MijT);
  NTSUB (MijT, AiX, MijT);
  NTADD (MijT, AiAj, MijT);

  coef = 1.0 / (double) nadj * nadj;

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
static void adj_average_gradient (EPR_ELEMENT *base, EPR_ELEMENT *ele, double *conf, double F[9])
{
  EPR_ELEMENT **nei, **nee;
  double *grad, coef;
  int shift;

  coef = 1.0 / (double) ele->nadj;

  SET9 (F, 0.0);

  for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
  {
    shift =  9 * ((*nei) - base);
    grad = &conf [shift];
    NNADD (grad, F, F);
  }

  SCALE9 (F, coef);
}

/* Lame's lambda coefficient */
static double lambda (double young, double poisson)
{
  return young*poisson / ((1.0 + poisson)*(1.0 - 2.0*poisson));
}

/* Lame's mi coefficient */
static double mi (double young, double poisson)
{
  return young / (2.0*(1.0 + poisson));
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
	*lin,
	 F [9],
	 P [9],
	 value,
	 B [3],
	 coef;
  int shift;

  msh = bod->priv;

  blas_dscal (bod->dofs, 0.0, force, 1); /* zero force */

  lin = force + bod->dofs - 3;

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

      shp = epr_shape_stab (bod, point); /* get element stabbed by the point */

      if (shp)
      {
	ele = shp->epr;

	if (frc->kind & CONVECTED) /* obtain spatial force */
	{
	  adj_average_gradient (msh->ele, ele, bod->conf, F);
	  NVMUL (F, f, g);
	  COPY (g, f);
	}

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

	lin [0] += f [0];
	lin [1] += f [1];
	lin [2] += f [2];
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
      adj_average_gradient (msh->ele, *nei, bod->conf, F);
      SVK_Stress_R (lambda, mi, (*nei)->vol, F, P);
      coef = 1.0 / (double) (*nei)->nadj;
      NNSUBMUL (cur, coef, P, cur); /* force = external - internal */
    }
  }
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
 * velocity space to the local cartesian space at 'X' in given 'base' */
static void epr_operator_block_H (double *X, double *A, double *base, double *H) /* H is assumed all zero on entry */
{
  double B [3];

  SUB (X, A, B); 

  H [0] = base [0] * B [0];
  H [3] = base [0] * B [1];
  H [6] = base [0] * B [2];
  H [9] = base [1] * B [0];
  H [12] = base [1] * B [1];
  H [15] = base [1] * B [2];
  H [18] = base [2] * B [0];
  H [21] = base [2] * B [1];
  H [24] = base [2] * B [2];

  H [1] = base [3] * B [0];
  H [4] = base [3] * B [1];
  H [7] = base [3] * B [2];
  H [10] = base [4] * B [0];
  H [13] = base [4] * B [1];
  H [16] = base [4] * B [2];
  H [19] = base [5] * B [0];
  H [22] = base [5] * B [1];
  H [25] = base [5] * B [2];

  H [2] = base [6] * B [0];
  H [5] = base [6] * B [1];
  H [8] = base [6] * B [2];
  H [11] = base [7] * B [0];
  H [14] = base [7] * B [1];
  H [17] = base [7] * B [2];
  H [20] = base [8] * B [0];
  H [23] = base [8] * B [1];
  H [26] = base [8] * B [2];
}

/* calculate cauche stress for a pseudo-rigid body */
void epr_cauchy (BODY *bod, EPR_MESH *msh, EPR_ELEMENT *ele, double *stress)
{
  double P [9], J, F [9], lambda, mi;

  lame_coefs (bod, ele, &lambda, &mi);
  adj_average_gradient (msh->ele, ele, bod->conf, F);

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
  EPR_MESH *msh;
  SHAPE *shq;

  /* get some global mass characteristics */
  SHAPE_Char (shp, &bod->ref_volume, bod->ref_center, bod->ref_tensor);
  bod->ref_mass = bod->ref_volume * mat->density;

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

  bod->dofs = 9 * msh->nele + 3;

  ERRMEM (bod->conf = calloc (4 * bod->dofs, sizeof (double))); /* conf, velo, vel0, force */
  bod->velo = bod->conf + bod->dofs;

  for (double *F = bod->conf, *E = F + 9 * msh->nele; F < E; F += 9) { IDENTITY (F); }

  set_up_bulk_materials (msh->ele, msh->nele);

  bod->priv = msh;
}

/* return configuration size */
int EPR_Conf_Size (BODY *bod)
{
  return bod->dofs;
}

/* overwrite state */
void EPR_Overwrite_State (BODY *bod, double *q, double *u)
{
  double *x, *y;

  for (x = bod->conf, y = x + bod->dofs; x < y; x ++, q ++) *x = *q;
  for (x = bod->velo, y = x + bod->dofs; x < y; x ++, u ++) *x = *u;
}

/* set initial rigid motion velocity */
void EPR_Initial_Velocity (BODY *bod, double *linear, double *angular)
{
  double *u = bod->velo, *e = u + bod->dofs - 3;

  if (angular)
  {
    for (; u < e; u += 9) { VECSKEW (angular, u); }
  }

  if (linear) { COPY (linear, e); }
}

/* initialise dynamic time stepping */
void EPR_Dynamic_Init (BODY *bod, SCHEME scheme)
{
  EPR_ELEMENT *ele, *end, **nei, **nee, **nej, **nff;
  double *xx, *chars, *sta, *eul, den;
  int m, n, *p, *i, ii, jj, *pp, *pi;
  MEM mapmem, valmem;
  MAP **col, *item;
  EPR_MESH *msh;

  if (bod->inverse) MX_Destroy (bod->inverse);

  msh = bod->priv;

  chars = compute_epr_chars (msh->ele, msh->nele);

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
	sum_up_inertia (ii, jj, (*nei)->A, (*nej)->A, eul, sta, ele->nadj, &mapmem, &valmem, col);
      }
    }
  }

  /* sum up scalar mass at the end */
  sum_up (m-3, n-3, bod->ref_volume, &mapmem, &valmem, col);
  sum_up (m-2, n-2, bod->ref_volume, &mapmem, &valmem, col);
  sum_up (m-1, n-1, bod->ref_volume, &mapmem, &valmem, col);

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
  for (jj = 0, xx = bod->inverse->x, den = bod->mat->density; jj < n; jj ++)
  {
    for (item = MAP_First (col[jj]); item; item = MAP_Next (item), xx ++)
    {
      *xx = (*((double*) item->data)) * den;
    }
  }

  /* compute inverse factorization */
  MX_Inverse (bod->inverse, bod->inverse);

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
  double half = 0.5 * step,
        *force = EPR_FORCE (bod);

  blas_dcopy (bod->dofs, bod->velo, 1, EPR_VEL0 (bod), 1); /* save u(t) */
  blas_daxpy (bod->dofs, half, bod->velo, 1, bod->conf, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
  epr_dynamic_force (bod, time+half, step, force);  /* f(t+h/2) = fext (t+h/2) - fint (q(t+h/2)) */
  MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) = u(t) + inv (M) * h * f(t+h/2) */
}

/* perform the final half-step of the dynamic scheme */
void EPR_Dynamic_Step_End (BODY *bod, double time, double step)
{
  double half = 0.5 * step,
        *force = EPR_FORCE (bod);
      
  epr_constraints_force (bod, force); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) += inv (M) * h * r */
  blas_daxpy (bod->dofs, half, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
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
  double *force = EPR_FORCE (bod);
      
  epr_constraints_force (bod, force); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  MX_Matvec (step, bod->inverse, force, 1.0, bod->velo); /* u(t+h) += inv (A) * h * r */
  blas_daxpy (bod->dofs, step, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t) + h * u(t+h) */
}

/* motion x = x (X, state) */
void EPR_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x)
{
  double *conf, *F, *a, *A, coef, B [3];
  EPR_ELEMENT *ele, **nei, **nee;
  EPR_MESH *msh;
  int shift;

  msh = bod->priv;
  ele = shp->epr;
  conf = bod->conf;
  a = conf + bod->dofs - 3;
  coef = 1.0  / (double) ele->nadj;
  SET (x, 0.0);

  for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
  {
    shift = 9 * ((*nei) - msh->ele);
    F = &conf [shift];
    A = (*nei)->A;

    SUB (X, A, B);
    TVADDMUL (x, F, B, x); /* transpose, as F is stored row-wise */
  }

  SCALE (x, coef); /* scale by shape function value */

  ADD (x, a, x); /* add motion of the mass center */
}

/* inverse motion X = X (x, state) */
void EPR_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X)
{
  double *conf, *F, *a, *A, coef, B [3], C [3], FAVG [9], IF [9], det;
  EPR_ELEMENT *ele, **nei, **nee;
  EPR_MESH *msh;
  int shift;

  msh = bod->priv;
  ele = shp->epr;
  conf = bod->conf;
  a = conf + bod->dofs - 3;
  coef = 1.0  / (double) ele->nadj;

  SET9 (FAVG, 0.0);
  SET (B, 0.0);

  for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
  {
    shift = 9 * ((*nei) - msh->ele);
    F = &conf [shift];
    A = (*nei)->A;

    TVADDMUL (B, F, A, B); /* B = (x-a) + SUM Ni Fi Ai; transpose, as F is stored row-wise */
    NTADD (FAVG, F, FAVG); /* FAVG = SUM Ni Fi */
  }

  SCALE (B, coef);  /* scale by shape function value */
  SUB (x, a, C);
  ADD (B, C, B); /* B = (x-a) + SUM Ni Fi Ai */

  SCALE9 (FAVG, coef); /* FAVG = SUM Ni Fi */

  INVERT (FAVG, IF, det);
  ASSERT (det > 0.0, ERR_BOD_MOTION_INVERT);
  NVMUL (IF, B, X);
  NVADDMUL (C, IF, a, X); /* X = inv (SUM Ni Fi) * {(x - a) + SUM Ni Fi Ai} */
}

/* obtain velocity at (element, point), expressed in the local 'base' => all entities are spatial */
void EPR_Local_Velo (BODY *bod, VELOTIME time, SHAPE *shp, void *gobj, double *point, double *base, double *velo)
{
  double *genvel, *F, *a, *A, coef, B [3], v [3];
  EPR_ELEMENT *ele, **nei, **nee;
  EPR_MESH *msh;
  int shift;

  msh = bod->priv;
  ele = shp->epr;
  genvel = time == CURVELO ? bod->velo : EPR_VEL0 (bod);
  a = genvel + bod->dofs - 3;
  coef = 1.0  / (double) ele->nadj;
  SET (v, 0.0);

  for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
  {
    shift = 9 * ((*nei) - msh->ele);
    F = &genvel [shift];
    A = (*nei)->A;

    SUB (point, A, B);
    TVADDMUL (v, F, B, v); /* v = SUM Ni dF/dti (point - Ai) + da/dt; transpose, as F is stored row-wise */
  }

  SCALE (v, coef); /* scale by shape function value */

  ADD (v, a, v); /* v = SUM Ni dF/dti (point - Ai) + da/dt */

  TVMUL (base, v, velo); /* express in local base */
}

/* return transformation operator from the generalised to the local velocity space at (element, point, base) */
MX* EPR_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base)
{
  EPR_ELEMENT *ele, **nei, **nee;
  double *tail, coef;
  EPR_MESH *msh;
  int shift;
  MX *H;

  msh = bod->priv;
  ele = shp->epr;
  coef = 1.0 / (double) ele->nadj;

  H = MX_Create (MXDENSE, 3, bod->dofs, NULL, NULL); /* initially a zero matrix */

  for (nei = ele->adj, nee = nei + ele->nadj; nei < nee; nei ++)
  {
    shift = 27 * ((*nei) - msh->ele);
    epr_operator_block_H (point, (*nei)->A, base, &H->x[shift]);
  }

  MX_Scale (H, coef); /* scale by shape function value */

  tail = &H->x[3 * bod->dofs - 9];
  TNCOPY (base, tail); /* copy last bit from the linear motion */

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
