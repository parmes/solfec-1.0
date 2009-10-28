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
#include "dom.h"
#include "fem.h"
#include "but.h"
#include "alg.h"
#include "bla.h"
#include "hyb.h"
#include "cvi.h"
#include "svk.h"
#include "err.h"

typedef double (*node_t) [3]; /* mesh node */

#define FEM_VEL0(bod) ((bod)->velo + (bod)->dofs)
#define FEM_FORCE(bod) ((bod)->velo + (bod)->dofs * 2)
#define FEM_MESH(bod) ((bod)->msh ? (bod)->msh : (bod)->shape->data)
#define FEM_MATERIAL(bod, ele) ((ele)->mat ? (ele)->mat : (bod)->mat)

/* copy node coordinates into a local table */
static void load_nodes (node_t heap, int type, int *nodes, node_t stack)
{
  int n;

  for (n = 0; n < type; n ++)
  { COPY (heap [nodes [n]], stack [n]); }
}

/* linear tetrahedron local to global point transformation */
inline static void tet_local_to_global (node_t nodes, double *local, double *point)
{
  /* TODO */
}

/* linear tetrahedron global to local point transformation */
inline static void tet_global_to_local (node_t nodes, double *point, double *local)
{
  /* TODO */
}

/* linear hexahedron local to global point transformation */
inline static void hex_local_to_global (node_t nodes, double *local, double *point)
{
  /* TODO */
}

/* linear hexahedron global to local point transformation */
inline static void hex_global_to_local (node_t nodes, double *point, double *local)
{
  /* TODO */
}

/* linear tetrahedron shape functions */
inline static void tet_o1_shapes (double *point, double *shapes)
{
  /* TODO */
}

/* linear hexahedron shape functions */
inline static void hex_o1_shapes (double *point, double *shapes)
{
  /* TODO */
}

/* linear tetrahedron shape functions */
inline static void tet_o1_derivs (double *point, double *derivs)
{
  /* TODO */
}

/* linear hexahedron shape functions */
inline static void hex_o1_derivs (double *point, double *derivs)
{
  /* TODO */
}

/* lump linear tetrahedron mass */
inline static void tet_o1_lump (node_t nodes, double density, double **out)
{
  /* TODO */
}

/* lump linear hexahedron mass */
inline static void hex_o1_lump (node_t nodes, double density, double **out)
{
  /* TODO */
}

/* compute linear tetrahedron body force contribution */
inline static void tet_o1_body_force (node_t nodes, double density, double *f, double *g)
{
  /* TODO */
}

/* compute linear hexahedron body force contribution */
inline static void hex_o1_body_force (node_t nodes, double density, double *f, double *g)
{
  /* TODO */
}

/* compute linear tetrahedron internal force contribution */
inline static void tet_o1_internal_force (node_t nodes, BULK_MATERIAL *mat, double *g)
{
  /* TODO */
}

/* compute linear hexahedron internal force contribution */
inline static void hex_o1_internal_force (node_t nodes, BULK_MATERIAL *mat, double *g)
{
  /* TODO */
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
    case 4: tet_o1_lump (nodes, density, out); break;
    case 8: hex_o1_lump (nodes, density, out); break;
    case 5:
      COPY (nodes [4], nodes [5]); out [5] = out [4];
      COPY (nodes [4], nodes [6]); out [6] = out [4];
      COPY (nodes [4], nodes [7]); out [7] = out [4];
      hex_o1_lump (nodes, density, out);
      break;
    case 6:
      COPY (nodes [5], nodes [7]); out [7] = out [5];
      COPY (nodes [5], nodes [6]); out [6] = out [5];
      COPY (nodes [4], nodes [5]); out [5] = out [4];
      COPY (nodes [3], nodes [4]); out [4] = out [3];
      COPY (nodes [2], nodes [3]); out [3] = out [2];
      hex_o1_lump (nodes, density, out);
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
  double *q = bod->conf;
  int i, j, k, n;

  IDENTITY (F);

  switch (bod->form)
  {
  case FEM_O1:
  {
    double derivs [24],
	   *p [8];

    n = ele->type;

    for (i = 0; i < n; i ++) p [i] = &q [ele->nodes [i] * 3];

    switch (n)
    {
    case 4: tet_o1_derivs (point, derivs); break; /* N1,1 N1,2 N1,3 N2,1 N2,3 ... */
    case 8: hex_o1_derivs (point, derivs); break;
    case 5:
      hex_o1_derivs (point, derivs);
      ADD (derivs + 15, derivs + 12, derivs + 12);
      ADD (derivs + 18, derivs + 12, derivs + 12);
      ADD (derivs + 21, derivs + 12, derivs + 12);
      break;
    case 6:
      hex_o1_derivs (point, derivs);
      ADD (derivs + 9, derivs + 6, derivs + 6);
      COPY (derivs + 12, derivs + 9);
      COPY (derivs + 15, derivs + 12);
      ADD (derivs + 21, derivs + 18, derivs + 15);
      break;
    }

    for (i = 0; i < 3; i ++)
      for (j = 0; j < 3; j ++)
	for (k = 0; k < n; k ++) F [3*j+i] += p[k][i] * derivs [3*k+j];
  }
  break;
  case FEM_O2:
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
  }
  break;
  }
}

/* compute element shape functions at a local point */
static MX* shape_functions (FEMFORM form, MESH *msh, ELEMENT *ele, double *point)
{
  MX *N = NULL;

  switch (form)
  {
  case FEM_O1:
  {
    double shapes [8], *x;
    static int p [] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24},
	       i [] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2};

    switch (ele->type)
    {
      case 4:
	tet_o1_shapes (point, shapes);
	N = MX_Create (MXCSC, 3, 12, p, i);
	x = N->x;
	x [0] = x [1]  = x [2]  = shapes [0];
	x [3] = x [4]  = x [5]  = shapes [1];
	x [6] = x [7]  = x [8]  = shapes [2];
	x [9] = x [10] = x [11] = shapes [3];
	break;
      case 8:
	hex_o1_shapes (point, shapes);
	N = MX_Create (MXCSC, 3, 24, p, i);
	x = N->x;
	x [0]  = x [1]  = x [2]  = shapes [0];
	x [3]  = x [4]  = x [5]  = shapes [1];
	x [6]  = x [7]  = x [8]  = shapes [2];
	x [9]  = x [10] = x [11] = shapes [3];
	x [12] = x [13] = x [14] = shapes [4];
	x [15] = x [16] = x [17] = shapes [5];
	x [18] = x [19] = x [20] = shapes [6];
	x [21] = x [22] = x [23] = shapes [7];
	break;
      case 5:
	hex_o1_shapes (point, shapes);
	N = MX_Create (MXCSC, 3, 15, p, i);
	x = N->x;
	x [0]  = x [1]  = x [2]  = shapes [0];
	x [3]  = x [4]  = x [5]  = shapes [1];
	x [6]  = x [7]  = x [8]  = shapes [2];
	x [9]  = x [10] = x [11] = shapes [3];
	x [12] = x [13] = x [14] = (shapes [4] + shapes [5] + shapes [6] + shapes [7]);
	break;
      case 6:
	hex_o1_shapes (point, shapes);
	N = MX_Create (MXCSC, 3, 18, p, i);
	x = N->x;
	x [0]  = x [1]  = x [2]  = shapes [0];
	x [3]  = x [4]  = x [5]  = shapes [1];
	x [6]  = x [7]  = x [8]  = (shapes [2] + shapes [3]);
	x [9]  = x [10] = x [11] = shapes [4];
	x [12] = x [13] = x [14] = shapes [5];
	x [15] = x [16] = x [17] = (shapes [6] + shapes [7]);
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

/* copute point force contribution */
static void point_force (FEMFORM form, MESH *msh, ELEMENT *ele, double *point, double *f, double *g)
{
  MX *N = shape_functions (form, msh, ele, point);
  MX_Matvec (1.0, MX_Tran (N), f, 0.0, g);
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
    case 4: tet_o1_body_force (nodes, density, f, g); break;
    case 8: hex_o1_body_force (nodes, density, f, g); break;
    case 5:
      COPY (nodes [4], nodes [5]);
      COPY (nodes [4], nodes [6]);
      COPY (nodes [4], nodes [7]);
      hex_o1_body_force (nodes, density, f, g);
      break;
    case 6:
      COPY (nodes [5], nodes [7]);
      COPY (nodes [5], nodes [6]);
      COPY (nodes [4], nodes [5]);
      COPY (nodes [3], nodes [4]);
      COPY (nodes [2], nodes [3]);
      hex_o1_body_force (nodes, density, f, g);
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

/* copute internal force contribution */
static void internal_force (BODY *bod, MESH *msh, ELEMENT *ele, double *g)
{
  BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
  double nodes [8][3];

  load_nodes (msh->ref_nodes, ele->type, ele->nodes, nodes);

  switch (bod->form)
  {
  case FEM_O1:
    switch (ele->type)
    {
    case 4: tet_o1_internal_force (nodes, mat, g); break;
    case 8: hex_o1_internal_force (nodes, mat, g); break;
    case 5:
      COPY (nodes [4], nodes [5]);
      COPY (nodes [4], nodes [6]);
      COPY (nodes [4], nodes [7]);
      hex_o1_internal_force (nodes, mat, g);
      break;
    case 6:
      COPY (nodes [5], nodes [7]);
      COPY (nodes [5], nodes [6]);
      COPY (nodes [4], nodes [5]);
      COPY (nodes [3], nodes [4]);
      COPY (nodes [2], nodes [3]);
      hex_o1_internal_force (nodes, mat, g);
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

/* gather elements with given node */
static void elements_around_node (int node, ELEMENT *ele, MEM *mem, SET **set)
{
  int n;

  for (n = 0; n < ele->type; n ++)
  {
    if (ele->nodes [n] == node)
    {
      ele->type = -ele->type; /* mark as visited */
      SET_Insert (mem, set, ele, NULL); /* add to set */
      for (n = 0; n < ele->neighs; n ++) elements_around_node (node, ele->adj [n], mem, set); /* recurse */
      break; /* done */
    }
  }
}

/* compute Cauchy stress at node */
static void fem_nodal_cauchy (BODY *bod, MESH *msh, ELEMENT *ele, int node, double *values)
{
  double F [9], J, P [9], X [3], point [3];
  SET *set, *item;
  MEM mem;

  set = NULL;
  MEM_Init (&mem, sizeof (SET), 64);
  elements_around_node (node, ele, &mem, &set);
  COPY (msh->ref_nodes [node], X);

  SET6 (values, 0.0);
  for (item = SET_First (set); item; item = SET_Next (item)) /* sum up stresses around the node */
  {
    ele = item->data;
    ASSERT_DEBUG (ele->type < 0, "Inconsistent elements around a node");
    ele->type = -ele->type; /* unmark */
    BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);

    referential_to_local (msh, ele, X, point);
    deformation_gradient (bod, msh, ele, point, F);

    J = SVK_Stress_C (lambda (mat->young, mat->poisson), mi (mat->young, mat->poisson), 1.0, F, P); /* column-wise, per unit volume */

    values [0] += (F[0]*P[0]+F[3]*P[1]+F[6]*P[2])/J; /* sx  */
    values [1] += (F[1]*P[3]+F[4]*P[4]+F[7]*P[5])/J; /* sy  */
    values [2] += (F[2]*P[6]+F[5]*P[7]+F[8]*P[8])/J; /* sz  */
    values [3] += (F[0]*P[3]+F[3]*P[4]+F[6]*P[5])/J; /* sxy */
    values [4] += (F[0]*P[6]+F[3]*P[7]+F[6]*P[8])/J; /* sxz */
    values [5] += (F[2]*P[0]+F[5]*P[1]+F[8]*P[2])/J; /* syz */
  }

  if (set)
  {
    X [0] = 1.0 / (double) SET_Size (set);
    SCALE6 (values, X[0]); /* averaged stress */
  }

  MEM_Release (&mem);
}

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

/* return element stabbed by a referential point */
static ELEMENT* stabbed_element (MESH *msh, ELEMENT **ele, int nele, double *X)
{
  for (; nele > 0; ele ++, nele --)
    if (ELEMENT_Contains_Point (msh, *ele, X, 1)) return *ele;

  return NULL;
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
  N = shape_functions (bod->form, msh, ele, point);
  H = MX_Matmat (1.0, &base_trans, N, 0.0, NULL);
  MX_Destroy (N);
  return H;
}

/* accumulate constraints reaction */
inline static void fem_constraints_force_accum (BODY *bod, MESH *msh, ELEMENT *ele, double *X, double *base, double *R, short isma, double *forc)
{
  MX *H = gen_to_loc_operator (bod, msh, ele, X, base);

  if (isma) MX_Matvec (1.0, MX_Tran (H), R, -1.0, forc);
  else MX_Matvec (1.0, MX_Tran (H), R, 1.0, forc);

  MX_Destroy (H);
}

/* compute constraints force */
static void fem_constraints_force (BODY *bod, double *forc)
{
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  SET *node;

  msh = FEM_MESH (bod);

  blas_dscal (bod->dofs, 0.0, forc, 1);

  for (node = SET_First (bod->con); node; node = SET_Next (node))
  {
    CON *con = node->data;
    short isma = (bod == con->master);
    double *X = (isma ? con->mpnt : con->spnt);

    if (bod->msh)
    {
      cvx = (isma ? con->mgobj : con->sgobj);
      ele = stabbed_element (msh, cvx->ele, cvx->nele, X);
    }
    else ele = (isma ? con->mgobj : con->sgobj);

    fem_constraints_force_accum (bod, msh, ele, X, con->base, con->R, isma, forc);
  }

#if MPI
  for (MAP *node = MAP_First (bod->conext); node; node = MAP_Next (node))
  {
    CONEXT *con = node->data;

    if (bod->msh)
    {
      cvx = con->sgp->gobj;
      ele = stabbed_element (msh, cvx->ele, cvx->nele, con->point);
    }
    else ele = con->sgp->gobj;

    fem_constraints_force_accum (bod, msh, ele, con->point, con->base, con->R, con->isma, forc);
  }
#endif
}

/* compute out of balance force = fext - fint */
static void fem_dynamic_force (BODY *bod, double time, double step, double *force)
{
  MESH *msh = FEM_MESH (bod);
  FEMFORM form = bod->form;
  ELEMENT *ele;
  double g [24],
	 f [3],
	 point [3],
	 value,
	*v,
	*w;
  int bulk,
      i;

  /* add point forces */
  for (FORCE *frc = bod->forces; frc; frc = frc->next)
  {
    if (frc->func)
    {
      ERRMEM (v = calloc (bod->dofs, sizeof (double)));
      frc->func (frc->data, frc->call, bod->dofs, bod->conf, bod->dofs, bod->velo, time, step, v);
      blas_daxpy (bod->dofs, 1.0, v, 1, force, 1);
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

	point_force (form, msh, ele, point, f, g);

	for (i = 0, v = g; i < ele->type; i ++, v += 3)
	{
	  w = &force [ele->nodes [i] * 3];
	  ADD (w, v, w);
	}
      }
    }
  }

  /* add gravitation */
  if (DOM(bod->dom)->gravval)
  {
    COPY (DOM(bod->dom)->gravdir, f);
    value = TMS_Value (DOM(bod->dom)->gravval, time);
    SCALE (f, value);

    for (ele = msh->surfeles, bulk = 0; ele; )
    {
      body_force (bod, msh, ele, f, g);

      for (i = 0, v = g; i < ele->type; i ++, v += 3)
      {
	w = &force [ele->nodes [i] * 3];
	ADD (w, v, w);
      }

      if (bulk) ele = ele->next;
      else if (ele->next) ele = ele->next;
      else ele = msh->bulkeles, bulk = 1;
    }
  }

  /* subtract internal forces */
  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    internal_force (bod, msh, ele, g);

    for (i = 0, v = g; i < ele->type; i ++, v += 3)
    {
      w = &force [ele->nodes [i] * 3];
      SUB (w, v, w);
    }

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }
}

/* compute inv (M) * K for an element */
static MX* inverse_mass_times_stiffencess (BODY *bod, MESH *msh, ELEMENT *ele)
{
  /* TODO */
  return NULL;
}

/* element and convex bounding boxe intersection callback; compute
 * their volumetric intersection to be later used for integration */
static void* overlap (void *data, BOX *one, BOX *two)
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
  tri = cvi (cvx->cur, cvx->nver, pla, cvx->nfac, vertices, n, planes, k, &m);
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

  return NULL;
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
    MEM boxmem;

    msh_sgp = SGP_Create (&msh_shp, &msh_nsgp);
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
  }
  else msh = shp->data; /* retrive the mesh pointer from the shape */

  if (form == FEM_O1)
  {
    bod->dofs = msh->nodes_count * 3;
    ERRMEM (bod->conf = calloc (4, bod->dofs * sizeof (double))); /* configuration, velocity, previous velocity, force */
    bod->velo = bod->conf + bod->dofs;
  }
  else
  {
    /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
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
void FEM_Dynamic_Init (BODY *bod, SCHEME scheme)
{
  MESH *msh = FEM_MESH (bod);
  int bulk,
     *p,
     *i,
      n,
      k;
  ELEMENT *ele;
  double *x, *y;

  if (bod->inverse) MX_Destroy (bod->inverse);

  n = bod->dofs;

  ERRMEM (p = malloc (sizeof (int [n+1])));
  ERRMEM (i = malloc (sizeof (int [n])));

  for (k = 0, p [n] = n; k < n; k ++) p [k] = i [k] = k; /* diagonal pattern */

  bod->inverse = MX_Create (MXCSC, n, n, p, i);
  x = bod->inverse->x;

  for (ele = msh->surfeles, bulk = 0; ele; )
  {
    lump_mass_matrix (bod, msh, ele, x);

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }

  /* invert diagonal */
  for (y = x + bod->dofs; x < y; x ++)
  {
    ASSERT (*x > 0.0, ERR_FEM_MASS_NOT_SPD);
    (*x) = 1.0 / (*x);
  }
}

/* estimate critical step for the dynamic scheme */
double FEM_Dynamic_Critical_Step (BODY *bod)
{
  MESH *msh = FEM_MESH (bod);
  double step, tcrit, eigmax;
  ELEMENT *ele;
  int bulk;
  MX *IMK;

  for (ele = msh->surfeles, bulk = 0, step = DBL_MAX; ele; )
  {
    IMK = inverse_mass_times_stiffencess (bod, msh, ele); /* element inv (M) * K */
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

/* perform the initial half-step of the dynamic scheme */
void FEM_Dynamic_Step_Begin (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  double half = 0.5 * step,
	*x = bod->inverse->x,
	*u0 = FEM_VEL0 (bod),
	*f = FEM_FORCE (bod),
	*q = bod->conf,
	*u = bod->velo,
	*e = u + n;

  blas_dcopy (n, u, 1, u0, 1); /* save u (t) */
  blas_daxpy (n, half, u, 1, q, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
  fem_dynamic_force (bod, time+half, step, f);  /* f(t+h/2) = fext (t+h/2) - fint (q(t+h/2)) */
  for (; u < e; u ++, x ++, f ++) (*u) += step * (*x) * (*f); /* u(t+h) = u(t) + inv (M) * h * f(t+h/2) */
}

/* perform the final half-step of the dynamic scheme */
void FEM_Dynamic_Step_End (BODY *bod, double time, double step)
{
  int n = bod->dofs;
  double half = 0.5 * step,
	*x = bod->inverse->x,
	*r = FEM_FORCE (bod),
	*u = bod->velo,
	*e = u + n;
  
  fem_constraints_force (bod, r); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  for (; u < e; u ++, x ++, r ++) (*u) += step * (*x) * (*r); /* u(t+h) += inv (M) * h * r */
  blas_daxpy (n, half, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */

  if (bod->msh) /* update current mesh nodes in case of a separate mesh */
  {
    MESH *msh = bod->msh;
    double *q = bod->conf,
	  (*ref) [3] = msh->ref_nodes,
	  (*cur) [3] = msh->cur_nodes,
	  (*end) [3] = cur + msh->nodes_count;

    for (;cur < end; q += 3, ref ++, cur ++) { ADD (ref [0], q, cur [0]); }
  }
}

/* initialise static time stepping */
void FEM_Static_Init (BODY *bod)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* perform the initial half-step of the static scheme */
void FEM_Static_Step_Begin (BODY *bod, double time, double step)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* perform the final half-step of the static scheme */
void FEM_Static_Step_End (BODY *bod, double time, double step)
{
  /* TODO */ ASSERT (0, ERR_NOT_IMPLEMENTED);
}

/* motion x = x (X, t) */
void FEM_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x)
{
  double point [3];
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  MX *N;

  if (bod->msh)
  {
    msh = bod->msh;
    ASSERT_DEBUG_EXT ((cvx = gobj), "NULL geometric object for the separate mesh FEM scenario");
    ele = stabbed_element (msh, cvx->ele, cvx->nele, X);
  }
  else
  {
    msh = shp->data;
    ele = gobj;
  }

  if (ele)
  {
    COPY (X, x);
    referential_to_local (msh, ele, X, point);
    N = shape_functions (bod->form, msh, ele, point);
    MX_Matvec (1.0, N, bod->conf, 1.0, x); /* x = X + N q */
    MX_Destroy (N);
  }
  else /* NULL implies nodal update (X is within msh->ref_nodes) */
  {
    int n = (node_t) X - msh->ref_nodes;

    COPY (msh->cur_nodes [n], x);
  }
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
    ele = stabbed_element (msh, cvx->ele, cvx->nele, X);
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
void FEM_Local_Velo (BODY *bod, VELOTIME time, SHAPE *shp, void *gobj, double *X, double *base, double *velo)
{
  double *u = (time == CURVELO ? bod->velo : FEM_VEL0 (bod));
  double point [3], vglo [3];
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  MX *N;

  if (bod->msh)
  {
    msh = bod->msh;
    ASSERT_DEBUG_EXT ((cvx = gobj), "NULL geometric object for the separate mesh FEM scenario");
    ele = stabbed_element (msh, cvx->ele, cvx->nele, X);
  }
  else
  {
    msh = shp->data;
    ele = gobj;
  }

  referential_to_local (msh, ele, X, point);
  N = shape_functions (bod->form, msh, ele, point);
  MX_Matvec (1.0, N, u, 0.0, vglo); /* vglo = N u */
  TVMUL (base, vglo, velo); /* velo = base' vglo */
  MX_Destroy (N);
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
  if (bod->inverse)
  {
    double *x = bod->inverse->x,
	   *y = x + bod->dofs,
	   *u = bod->velo,
	   sum;

    for (sum = 0.0; x < y; x ++, u ++) sum += (*u)*(*u) / (*x);

    return 0.5 * sum;
  }

  return 0.0;
}

/* get some values at a node of a geometrical object */
void FEM_Nodal_Values (BODY *bod, SHAPE *shp, void *gobj, int node, VALUE_KIND kind, double *values)
{
  double *X, point [3];
  ELEMENT *ele;
  CONVEX *cvx;
  MESH *msh;
  int n;

  if (bod->msh)
  {
    msh = bod->msh;
    ASSERT_DEBUG_EXT ((cvx = gobj), "NULL geometric object for the separate mesh FEM scenario");
    ASSERT_DEBUG (node >= 0 && node < cvx->nver, "CONVEX node number out of bounds");
    X = cvx->ref + node * 3;
    ele = stabbed_element (msh, cvx->ele, cvx->nele, X);
    referential_to_local (msh, ele, X, point);
    n = 0;
  }
  else
  {
    msh = shp->data;
    ele = gobj;
    ASSERT_DEBUG (ele && node >= 0 && node < ele->type, "Invalid input when obtaining nodal FEM value");
    n = ele->nodes [node];
  }

  switch (kind)
  {
  case VALUE_DISPLACEMENT:
  {
    if (bod->msh)
    {
      MX *N = shape_functions (bod->form, msh, ele, point);
      MX_Matvec (1.0, N, bod->conf, 0.0, values);
      MX_Destroy (N);
    }
    else
    {
      double *q = &bod->conf [3 * n];
      COPY (q, values);
    }
  }
  break;
  case VALUE_VELOCITY:
  {
    if (bod->msh)
    {
      MX *N = shape_functions (bod->form, msh, ele, point);
      MX_Matvec (1.0, N, bod->velo, 0.0, values);
      MX_Destroy (N);
    }
    else
    {
      double *u = &bod->velo [3 * n];
      COPY (u, values);
    }
  }
  break;
  case VALUE_STRESS:
  {
    if (bod->msh) fem_element_cauchy (bod, msh, ele, point, values);
    else fem_nodal_cauchy (bod, msh, ele, n, values);
  }
  break;
  case VALUE_MISES:
  {
    double stress [6];

    if (bod->msh) fem_element_cauchy (bod, msh, ele, point, stress);
    else fem_nodal_cauchy (bod, msh, ele, n, stress);
    MISES (stress, values [0]);
  }
  break;
  case VALUE_STRESS_AND_MISES:
  {
    if (bod->msh) fem_element_cauchy (bod, msh, ele, point, values);
    else fem_nodal_cauchy (bod, msh, ele, n, values);
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

    switch (kind)
    {
    case VALUE_DISPLACEMENT:
    {
      MX *N = shape_functions (bod->form, msh, ele, point);
      MX_Matvec (1.0, N, bod->conf, 0.0, values);
      MX_Destroy (N);
    }
    break;
    case VALUE_VELOCITY:
    {
      MX *N = shape_functions (bod->form, msh, ele, point);
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

/* release FEM memory */
void FEM_Destroy (BODY *bod)
{
  free (bod->conf);
}
