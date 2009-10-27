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

#include "dom.h"
#include "fem.h"
#include "but.h"
#include "alg.h"
#include "bla.h"
#include "hyb.h"
#include "cvi.h"
#include "err.h"

typedef double (*node_t) [3]; /* mesh node */

#define FEM_VEL0(bod) ((bod)->velo + (bod)->dofs)
#define FEM_FORCE(bod) ((bod)->velo + (bod)->dofs * 2)

/* lump linear tetrahedron mass */
static void tet_o1_lump (node_t nodes, double density, double **out)
{
}

/* lump linear hexahedron mass */
static void hex_o1_lump (node_t nodes, double density, double **out)
{
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
	 copy [4][3],
	 density,
	*out [8];
  int i;

  density = ele->mat ? ele->mat->density : bod->mat->density;

  load_nodes (msh->ref_nodes, ele->type, ele->nodes, (node_t) nodes);

  for (i = 0; i < ele->type; i ++) out [i] = &x [3 * ele->nodes [i]];
 
  switch (bod->form)
  {
  case FEM_O1:
  {
    switch (ele->type)
    {
    case 4: tet_o1_lump ((node_t) nodes, density, out); break;
    case 8: hex_o1_lump ((node_t) nodes, density, out); break;
    case 5:
      COPY (nodes [0], copy [0]);
      COPY (nodes [1], copy [1]);
      COPY (nodes [2], copy [2]);
      COPY (nodes [4], copy [3]);
      out [3] = &x [3 * ele->nodes [4]];
      tet_o1_lump ((node_t) copy, density, out);

      COPY (nodes [2], copy [1]);
      COPY (nodes [3], copy [2]);
      out [1] = &x [3 * ele->nodes [2]];
      out [2] = &x [3 * ele->nodes [3]];
      tet_o1_lump ((node_t) copy, density, out);
      break;
    case 6:
      tet_o1_lump ((node_t) nodes, density, out);

      COPY (nodes [1], copy [0]);
      COPY (nodes [4], copy [1]);
      COPY (nodes [2], copy [2]);
      COPY (nodes [3], copy [3]);
      out [1] = &x [3 * ele->nodes [1]];
      out [2] = &x [3 * ele->nodes [4]];
      tet_o1_lump ((node_t) copy, density, out);

      COPY (nodes [4], copy [0]);
      COPY (nodes [5], copy [1]);
      out [1] = &x [3 * ele->nodes [4]];
      out [2] = &x [3 * ele->nodes [5]];
      tet_o1_lump ((node_t) copy, density, out);
      break;
    }
  }
  break;
  case FEM_O2:
  {
    ASSERT (0, ERR_NOT_IMPLEMENTED);
    /* TODO */
  }
  break;
  }
}

/* compute deformation gradient */
static void deformation_gradient (MESH *msh, ELEMENT *ele, double *point, double *F)
{
}

/* copute point force contribution */
static void point_force (MESH *msh, ELEMENT *ele, double *point, double *f, double *g)
{
}

/* copute body force contribution */
static void body_force (MESH *msh, ELEMENT *ele, double *f, double *g)
{
}

/* copute internal force contribution */
static void internal_force (BODY *bod, MESH *msh, ELEMENT *ele, double *g)
{
}

/* out of balance force = fext - fint */
static void fem_dynamic_force (BODY *bod, double time, double step, double *forc)
{
  MESH *msh = bod->shape->data;
  ELEMENT *ele;
  double g [24],
	 f [3],
	 value,
	*v,
	*w;
  int bulk,
      i;

  /* add external forces */
  for (FORCE *frc = bod->forces; frc; frc = frc->next)
  {
    if (frc->func)
    {
      ERRMEM (v = calloc (bod->dofs, sizeof (double)));
      frc->func (frc->data, frc->call, bod->dofs, bod->conf, bod->dofs, bod->velo, time, step, v);
      blas_daxpy (bod->dofs, 1.0, v, 1, forc, 1);
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
	if (frc->kind & CONVECTED)
	{ 
	  deformation_gradient (msh, ele, frc->ref_point, g);
	  NVMUL (g, f, g+9);
	  COPY (g+9, f);
	}

	point_force (msh, ele, frc->ref_point, f, g);

	for (i = 0, v = g; i < ele->type; i ++, v += 3)
	{
	  w = &forc [ele->nodes [i] * 3];
	  ADD (w, v, w);
	}
      }
    }
  }

  /* gravitation */
  if (DOM(bod->dom)->gravval)
  {
    COPY (DOM(bod->dom)->gravdir, f);
    value = TMS_Value (DOM(bod->dom)->gravval, time);
    SCALE (f, value);

    for (ele = msh->surfeles, bulk = 0; ele; )
    {
      body_force (msh, ele, f, g);

      for (i = 0, v = g; i < ele->type; i ++, v += 3)
      {
	w = &forc [ele->nodes [i] * 3];
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
      w = &forc [ele->nodes [i] * 3];
      SUB (w, v, w);
    }

    if (bulk) ele = ele->next;
    else if (ele->next) ele = ele->next;
    else ele = msh->bulkeles, bulk = 1;
  }
}

/* accumulate constraints reaction */
inline static void fem_constraints_force_accum (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *base, double *R, short isma, double *forc)
{
#if 0
  MX *H = FEM_Gen_To_Loc_Operator (bod, msh, ele, point, base);

  if (isma) MX_Matvec (1.0, MX_Tran (H), R, -1.0, forc);
  else MX_Matvec (1.0, MX_Tran (H), R, 1.0, forc);

  MX_Destroy (H);
#endif
}

/* calculate constraints force */
static void fem_constraints_force (BODY *bod, double *forc)
{
  SET *node;

  blas_dscal (bod->dofs, 0.0, forc, 1);

  for (node = SET_First (bod->con); node; node = SET_Next (node))
  {
    CON *con = node->data;
    short isma = (bod == con->master);
    double *point = (isma ? con->mpnt : con->spnt);
    MESH *msh = (isma ? con->mshp->data : con->sshp->data);
    ELEMENT *ele = (isma ? con->mgobj : con->sgobj);

    fem_constraints_force_accum (bod, msh, ele, point, con->base, con->R, isma, forc);
  }

#if MPI
  for (MAP *node = MAP_First (bod->conext); node; node = MAP_Next (node))
  {
    CONEXT *con = node->data;

    fem_constraints_force_accum (bod, con->sgp->shp->data, con->sgp->gobj, con->point, con->base, con->R, con->isma, forc);
  }
#endif
}

/* compute local coordinates of a spatial point */
static void spatial_to_local (MESH *msh, ELEMENT *ele, double *x, double *point)
{
}

/* copute referential coordinates of a local point */
static void local_to_referential (MESH *msh, ELEMENT *ele, double *point, double *X)
{
}

/* get shape functions of an element at given referential point */
static MX* shape_functions (MESH *msh, ELEMENT *ele, double *point)
{
  return NULL;
}

/* compute Cauchy stress at node */
static void fem_nodal_cauchy (BODY *bod, MESH *msh, int n, double *values)
{
}

/* compute Cauchy stress at point */
static void fem_element_cauchy (BODY *bod, MESH *msh, ELEMENT *ele, double *point, double *values)
{
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

  bod->dofs = msh->nodes_count * 3;
  ERRMEM (bod->conf = calloc (4, bod->dofs * sizeof (double))); /* configuration, velocity, previous velocity, force */
  bod->velo = bod->conf + bod->dofs;
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
  MESH *msh = bod->shape->data;
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
  MESH *msh = bod->shape->data;
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
    ASSERT_DEBUG (*x > 0.0, "Singular lumped mass matrix");
    (*x) = 1.0 / (*x);
  }
}

/* estimate critical step for the dynamic scheme */
double FEM_Dynamic_Critical_Step (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
  return 0.0;
}

/* perform the initial half-step of the dynamic scheme */
void FEM_Dynamic_Step_Begin (BODY *bod, double time, double step)
{
  double half = 0.5 * step,
	*velp = FEM_VEL0 (bod),
	*forc = FEM_FORCE (bod);

  blas_dcopy (bod->dofs, bod->velo, 1, velp, 1); /* save u (t) */
  blas_daxpy (bod->dofs, half, bod->velo, 1, bod->conf, 1); /* q(t+h/2) = q(t) + (h/2) * u(t) */
  fem_dynamic_force (bod, time+half, step, forc);  /* f(t+h/2) = fext (t+h/2) - fint (q(t+h/2)) */
  MX_Matvec (step, bod->inverse, forc, 1.0, bod->velo); /* u(t+h) = u(t) + inv (M) * h * f(t+h/2) */
}

/* perform the final half-step of the dynamic scheme */
void FEM_Dynamic_Step_End (BODY *bod, double time, double step)
{
  double half = 0.5 * step,
	*forc = FEM_FORCE (bod);
  
  fem_constraints_force (bod, forc); /* r = SUM (over constraints) { H^T * R (average, [t, t+h]) } */
  MX_Matvec (step, bod->inverse, forc, 1.0, bod->velo); /* u(t+h) += inv (M) * h * r */
  blas_daxpy (12, half, bod->velo, 1, bod->conf, 1); /* q (t+h) = q(t+h/2) + (h/2) * u(t+h) */
}

/* initialise static time stepping */
void FEM_Static_Init (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* perform the initial half-step of the static scheme */
void FEM_Static_Step_Begin (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* perform the final half-step of the static scheme */
void FEM_Static_Step_End (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* motion x = x (X, state) */
void FEM_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x)
{
#if 0
  if (ele)
  {
    double base [9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    MX *H = FEM_Gen_To_Loc_Operator (bod, msh, ele, X, base);

    COPY (X, x);
    MX_Matvec (1.0, H, bod->conf, 1.0, x); /* x = X + N q */
    MX_Destroy (H);
  }
  else /* ele == NULL implies nodal update (X is within the mesh->ref_nodes) */
  {
    int n = (node_t) X - msh->ref_nodes;

    COPY (msh->cur_nodes [n], x);
  }
#endif
}

/* inverse motion X = X (x, state) */
void FEM_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X)
{
#if 0
  double point [3];

  spatial_to_local (msh, ele, x, point);

  local_to_referential (msh, ele, point, X);
#endif
}

/* obtain velocity at (element, point), expressed in the local 'base' => all entities are spatial */
void FEM_Local_Velo (BODY *bod, VELOTIME time, SHAPE *shp, void *gobj, double *point, double *base, double *velo)
{
#if 0
  double *u = (time == CURVELO ? bod->velo : FEM_VEL0 (bod));
  MX *H = FEM_Gen_To_Loc_Operator (bod, msh, ele, point, base);

  MX_Matvec (1.0, H, u, 0.0, velo);
  MX_Destroy (H);
#endif
}

/* return transformation operator from the generalised to the local velocity space at (element, point, base) */
MX* FEM_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base)
{
#if 0
  double i [] = {0, 1, 2, 0, 1, 2, 0, 1, 2}, p [] = {0, 3, 6, 9};
  MX_CSC (base_trans, 9, 3, 3, p, i);
  MX *N, *H;

  TNCOPY (base, base_trans.x);

  N = shape_functions (msh, ele, point);

  H = MX_Matmat (1.0, &base_trans, N, 0.0, NULL);

  MX_Destroy (N);

  return H;
#endif
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
#if 0
  int n = ele->nodes [node];

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
    fem_nodal_cauchy (bod, msh, n, values);
  }
  break;
  case VALUE_MISES:
  {
    double stress [6];

    fem_nodal_cauchy (bod, msh, n, stress);
    MISES (stress, values [0]);
  }
  break;
  case VALUE_STRESS_AND_MISES:
  {
    fem_nodal_cauchy (bod, msh, n, values);
    MISES (values, values [6]);
  }
  break;
  }
#endif
}

/* get some values at a referential point */
void FEM_Point_Values (BODY *bod, double *point, VALUE_KIND kind, double *values)
{
  ELEMENT *ele;
  MESH *msh;

  ele = MESH_Element_Containing_Point (msh, point, 1);
  msh = bod->shape->data;

  if (ele)
  switch (kind)
  {
  case VALUE_DISPLACEMENT:
  {
    double x [3];

    FEM_Cur_Point (bod, msh, ele, point, x);
    SUB (x, point, values);
  }
  break;
  case VALUE_VELOCITY:
  {
    double base [9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    FEM_Local_Velo (bod, CURVELO, msh, ele, point, base, values);
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

/* release FEM memory */
void FEM_Destroy (BODY *bod)
{
  free (bod->conf);
}
