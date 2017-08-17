/*
 * eli.c
 * Copyright (C) 2011, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * ellipsoids
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

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "sol.h"
#include "shp.h"
#include "eli.h"
#include "alg.h"
#include "lap.h"
#include "gjk.h"
#include "pck.h"
#include "err.h"

/* compute the linear operator that transforms a unit sphere into an ellipsoid;
 * extract from it the rotation oprator and scaling coefficients for unit sphere */
static void sca_rot (double *c, double (*p) [3], double *sca, double *rot)
{
  double V [9] = {p[0][0] - c[0], p[0][1] - c[1], p[0][2] - c[2],
                  p[1][0] - c[0], p[1][1] - c[1], p[1][2] - c[2],
                  p[2][0] - c[0], p[2][1] - c[1], p[2][2] - c[2]}, U [9];
  int i;

  POLAR (V, 1E-10, rot, U, i); /* XXX: find out a better way */
  sca [0] = U [0];
  sca [1] = U [4];
  sca [2] = U [8];
}

/* create an ellipsoid */
ELLIP* ELLIP_Create (double *center, double *radii, int surface, int volume)
{
  ELLIP *eli;

  ERRMEM (eli = malloc (sizeof (ELLIP)));

  COPY (center, eli->ref_center);
  COPY (center, eli->ref_point [0]); eli->ref_point[0][0] += radii [0];
  COPY (center, eli->ref_point [1]); eli->ref_point[1][1] += radii [1];
  COPY (center, eli->ref_point [2]); eli->ref_point[2][2] += radii [2];

  COPY (center, eli->cur_center);
  COPY (center, eli->cur_point [0]); eli->cur_point[0][0] += radii [0];
  COPY (center, eli->cur_point [1]); eli->cur_point[1][1] += radii [1];
  COPY (center, eli->cur_point [2]); eli->cur_point[2][2] += radii [2];

  COPY (radii, eli->ref_sca);
  COPY (radii, eli->cur_sca);
  IDENTITY (eli->ref_rot);
  IDENTITY (eli->cur_rot);

  eli->surface = surface;
  eli->volume = volume;
  eli->mat = NULL;

  return eli;
}

/* create a copy of a ellipsoid */
ELLIP* ELLIP_Copy (ELLIP *eli)
{
  ELLIP *out;

  ERRMEM (out = malloc (sizeof (ELLIP)));
  memcpy (out, eli, sizeof (ELLIP));

  return out;
}

/* scaling of a ellipsoid; scale radius: r *= vector [0]; (set ref = cur) */
void ELLIP_Scale (ELLIP *eli, double *vector)
{
  double *ref [4] = {eli->ref_center,
                     eli->ref_point [0],
                     eli->ref_point [1],
                     eli->ref_point [2]};

  for (int i = 0; i < 4; i ++)
  {
    ref [i][0] *= vector [0];
    ref [i][1] *= vector [1];
    ref [i][2] *= vector [2];
  }

  ELLIP_Update (eli, NULL, NULL, NULL);
}

/* translation of a ellipsoid */
void ELLIP_Translate (ELLIP *eli, double *vector)
{
  double *ref [4] = {eli->ref_center,
                     eli->ref_point [0],
                     eli->ref_point [1],
                     eli->ref_point [2]};

  for (int i = 0; i < 4; i ++)
  {
    ref [i][0] += vector [0];
    ref [i][1] += vector [1];
    ref [i][2] += vector [2];
  }

  ELLIP_Update (eli, NULL, NULL, NULL);
}

/* rotation of a ellipsoid */
void ELLIP_Rotate (ELLIP *eli, double *point, double *vector, double angle)
{
  double R [9], omega [3];
  double *ref [4] = {eli->ref_center,
                     eli->ref_point [0],
                     eli->ref_point [1],
                     eli->ref_point [2]};

  angle *=  ALG_PI / 180.0;
  COPY (vector, omega); 
  NORMALIZE (omega); 
  SCALE (omega, angle);
  EXPMAP (omega, R);

  for (int i = 0; i < 4; i ++)
  {
    SUB (ref [i], point, omega);
    NVADDMUL (point, R, omega, ref [i]);
  }

  ELLIP_Update (eli, NULL, NULL, NULL);
}

/* cut through ellipsoid with a plane; return triangulated cross-section; vertices in the triangles
 * point to the memory allocated after the triangles memory; adjacency is not maintained */
TRI* ELLIP_Cut (ELLIP *eli, double *point, double *normal, int *m)
{
  /* TODO */
  WARNING_DEBUG (0, "Ellipsoid cutting has not been implemented yet!");
  *m = 0;
  return NULL;
}

/* split ellipsoid in two with plane defined by (point, normal); surfid corresponds to the new surface;
 * topoadj != 0 implies cutting from the point and through the topological adjacency only */
void ELLIP_Split (ELLIP *eli, double *point, double *normal, short topoadj, int surfid[2], CONVEX **one, CONVEX **two)
{
  /* TODO */
  WARNING_DEBUG (0, "Ellipsoid splitting has not been implemented yet!");
  *one = *two = NULL;
}

/* compute partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void ELLIP_Char_Partial (ELLIP *eli, int ref, double *vo, double *sx, double *sy, double *sz, double *eul)
{
  double v, a, b, c, e, *x, *R, tmp [9], tmq [9], eye [9];

  if (ref)
  {
    ELLIP *elj = ELLIP_Copy (eli);
    ELLIP_Update (elj, NULL, NULL, NULL);
    ELLIP_Char_Partial (elj, 0, vo, sx, sy, sz, eul);
    ELLIP_Destroy (elj);
    return;
  }

  R = eli->cur_rot;
  x = eli->cur_center;
  a = eli->cur_sca [0];
  b = eli->cur_sca [1];
  c = eli->cur_sca [2];
  v = (4.0/3.0)*ALG_PI*a*b*c;

  /* Stainer's theorem =>
   * J(i, j) = I(i, j) + M*(<x,x> * delta (i, j) - diadic (x, x))
   * hence it is possible to produce inertia tensor with respect
   * to the origin from one with respect to the center point */

  e = DOT (x, x);
  IDENTITY (eye);
  SCALEDIAG (eye, e);
  DIADIC (x, x, tmp);
  NNSUB (eye, tmp, tmp);
  SCALE9 (tmp, v);
  IDENTITY (eye); /* prepare I0 (i, j) */
  eye [0] = v * (b*b + c*c) / 5.0;
  eye [4] = v * (a*a + c*c) / 5.0;
  eye [8] = v * (b*b + a*a) / 5.0;
  NTMUL (eye, R, tmq);
  NNMUL (R, tmq, eye); /* I1 (i, j) = R I0 R' */
  NNADD (eye, tmp, tmp); /* tmp = inertia of the sphere with respect to x, y, z passing 0 (see above) */

  /* note that Inertia = Trace(Euler)*Identity - Euler,
   * hence Euler = 0.5 * Trace (Inertia)*Identity - Inertia,
   * as Trace(Inertia) = 3*Trace(Euler) - Trace(Euler) */

  e = 0.5 * TRACE (tmp);
  IDENTITY (eye);
  SCALEDIAG (eye, e);
  NNSUB (eye, tmp, tmp); /* tmp = euler tensor with repsect to x, y, z passing 0 */

  /* sum up */
  *vo += v;
  *sx += v * x [0];
  *sy += v * x [1];
  *sz += v * x [2];
  NNADD (tmp, eul, eul);

  /* TODO: make sure the above is correct */
}

/* get characteristics of the ellipsoid: volume, mass center, and Euler tensor (centered) */
void ELLIP_Char (ELLIP *eli, int ref, double *volume, double *center, double *euler)
{
  double vo, sx, sy, sz, cen [3], eul [9];

  vo = sx = sy = sz = 0.0;
  SET9 (eul, 0.0);

  ELLIP_Char_Partial (eli, ref, &vo, &sx, &sy, &sz, eul);

  cen [0] = sx / vo;
  cen [1] = sy / vo;
  cen [2] = sz / vo;

  eul [0] -= (2*sx - cen[0]*vo)*cen[0];
  eul [4] -= (2*sy - cen[1]*vo)*cen[1];
  eul [8] -= (2*sz - cen[2]*vo)*cen[2];
  eul [3] -= cen[0]*sy + cen[1]*sx - cen[0]*cen[1]*vo;
  eul [6] -= cen[0]*sz + cen[2]*sx - cen[0]*cen[2]*vo;
  eul [7] -= cen[1]*sz + cen[2]*sy - cen[1]*cen[2]*vo;
  eul [1] = eul[3];
  eul [2] = eul[6];
  eul [5] = eul[7];

  if (volume) *volume = vo;
  if (center) COPY (cen, center);
  if (euler) NNCOPY (eul, euler);
}

/* update extents of an individual ellipsoid */
void ELLIP_Extents (void *data, ELLIP *eli, double *extents)
{
  double *c = eli->cur_center,
	 *r = eli->cur_sca,
	 *R = eli->cur_rot,
          p [8][3] = {{- r[0], - r[1], - r[2]},
                      {+ r[0], - r[1], - r[2]}, 
                      {+ r[0], + r[1], - r[2]}, 
                      {- r[0], + r[1], - r[2]}, 
                      {- r[0], - r[1], + r[2]},
                      {+ r[0], - r[1], + r[2]}, 
                      {+ r[0], + r[1], + r[2]}, 
                      {- r[0], + r[1], + r[2]}},
	  q [3];

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;

  for (int i = 0; i < 8; i ++)
  {
    NVADDMUL (c, R, p[i], q);

    if (q [0] < extents [0]) extents [0] = q [0];
    if (q [1] < extents [1]) extents [1] = q [1];
    if (q [2] < extents [2]) extents [2] = q [2];
    if (q [0] > extents [3]) extents [3] = q [0];
    if (q [1] > extents [4]) extents [4] = q [1];
    if (q [2] > extents [5]) extents [5] = q [2];
  }

  extents [0] -= GEOMETRIC_EPSILON;
  extents [1] -= GEOMETRIC_EPSILON;
  extents [2] -= GEOMETRIC_EPSILON;
  extents [3] += GEOMETRIC_EPSILON;
  extents [4] += GEOMETRIC_EPSILON;
  extents [5] += GEOMETRIC_EPSILON;
}

/* compute extents of a ellipsoid */
void ELLIP_Extents_2 (ELLIP *eli, double *extents)
{
  ELLIP_Extents (NULL, eli, extents);
}

/* compute oriented extents of a ellipsoid */
void ELLIP_Oriented_Extents (ELLIP *eli, double *vx, double *vy, double *vz, double *extents)
{
  double *c = eli->cur_center,
	 *r = eli->cur_sca,
	 *R = eli->cur_rot,
          p [8][3] = {{- r[0], - r[1], - r[2]},
                      {+ r[0], - r[1], - r[2]}, 
                      {+ r[0], + r[1], - r[2]}, 
                      {- r[0], + r[1], - r[2]}, 
                      {- r[0], - r[1], + r[2]},
                      {+ r[0], - r[1], + r[2]}, 
                      {+ r[0], + r[1], + r[2]}, 
                      {- r[0], + r[1], + r[2]}},
	  q [3],
	  e [3];

  for (int i = 0; i < 8; i ++)
  {
    NVADDMUL (c, R, p[i], q);

    e [0] = DOT (vx, q);
    e [1] = DOT (vy, q);
    e [2] = DOT (vz, q);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [0] > extents [3]) extents [3] = e [0];
    if (e [1] > extents [4]) extents [4] = e [1];
    if (e [2] > extents [5]) extents [5] = e [2];
  }
}

/* return first not NULL bulk material for a ellipsoid */
void* ELLIP_First_Bulk_Material (ELLIP *eli)
{
  return eli->mat;
}

/* return ellipsoid containing a spatial point */
ELLIP* ELLIP_Containing_Point (ELLIP *eli, double *point)
{
  double *S = eli->cur_sca,
	 *R = eli->cur_rot,
	 *c = eli->cur_center,
	 oneps = 1.0 + GEOMETRIC_EPSILON,
	 v0 [3],
	 v1 [3];

  SUB (point, c, v0);
  TVMUL (R, v0, v1);
  v1 [0] /=  S [0];
  v1 [1] /=  S [1];
  v1 [2] /=  S [2];
  if (DOT (v1, v1) <= oneps*oneps) return eli;
  else return NULL;
}

/* does this ellipsoid contain the spatial point? */
int ELLIP_Contains_Point (void *dummy, ELLIP *eli, double *point)
{
  if (ELLIP_Containing_Point (eli, point)) return 1;
  else return 0;
}

/* return distance of a spatial point to the ellipsoid */
double ELLIP_Spatial_Point_Distance (void *dummy, ELLIP *eli, double *point)
{
  double q [3];

  return gjk_ellip_point (eli->cur_center, eli->cur_sca, eli->cur_rot, point, q);
}

/* update ellipsoid according to the given motion */
void ELLIP_Update (ELLIP *eli, void *body, void *shp, MOTION motion)
{
  SGP sgp = {shp, eli, GOBJ_ELLIP, NULL};
  double *ref = eli->ref_center,
	 (*ref_pnt) [3] = eli->ref_point,
	 *cur = eli->cur_center,
	 (*cur_pnt) [3] = eli->cur_point;

  if (motion)
  { 
    motion (body, &sgp, ref, cur);
    motion (body, &sgp, ref_pnt [0], cur_pnt [0]);
    motion (body, &sgp, ref_pnt [1], cur_pnt [1]);
    motion (body, &sgp, ref_pnt [2], cur_pnt [2]);

    BODY *bod = body;

    switch (bod->kind)
    {
    case OBS:
    case RIG:
    {
      double *R1 = bod->conf, *R0 = eli->ref_rot, *rot = eli->cur_rot;

      NNMUL (R1, R0, rot);
    }
    break;
    case PRB:
    {
      double *F = bod->conf, *sca0 = eli->ref_sca, *rot0 = eli->ref_rot, *sca1 = eli->cur_sca, *rot1 = eli->cur_rot;
      double U[9] = {1.0/(sca0[0]*sca0[0]), 0.0, 0.0, 0.0, 1.0/(sca0[1]*sca0[1]), 0.0, 0.0, 0.0, 1.0/(sca0[2]*sca0[2])};
      double A0[9], iF[9], det, X[3], Y[9], A[9];

      NTMUL (U, rot0, Y);
      NNMUL (rot0, Y, A0);

      TNCOPY (F, Y); /* T --> since deformation gradient is stored row-wise */

      INVERT (Y, iF, det);
      ASSERT_TEXT (det > 0.0, "det(F) <= 0.0 during ellipsoid update");

      NNMUL (A0, iF, Y);
      TNMUL (iF, Y, A);

      ASSERT_TEXT (lapack_dsyev ('V', 'U', 3, A, 3, X, Y, 9) == 0, "Eigen decomposition failed during ellipsoid update");

      if (DET(A) < 0.0) /* det(A) is 1.0 or -1.0 */
      {
	SCALE9 (A, -1.0); /* keep positive space orientation */
      }

      NNCOPY (A, rot1);

      sca1[0] = 1.0/sqrt(X[0]);
      sca1[1] = 1.0/sqrt(X[1]);
      sca1[2] = 1.0/sqrt(X[2]);
    }
    break;
    default:
    {
      ASSERT_TEXT (0, "Invalid body kind during ellipsoid update");
    }
    break;
    }
  }
  else
  {
    COPY (ref, cur);
    COPY (ref_pnt [0], cur_pnt [0]);
    COPY (ref_pnt [1], cur_pnt [1]);
    COPY (ref_pnt [2], cur_pnt [2]);

    sca_rot (eli->ref_center, eli->ref_point, eli->ref_sca, eli->ref_rot);

    COPY (eli->ref_sca, eli->cur_sca);
    NNCOPY (eli->ref_rot, eli->cur_rot);
  }
}

/* free ellipsoid */
void ELLIP_Destroy (ELLIP *eli)
{
  free (eli);
}

/* pack ellipsoid into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void ELLIP_Pack (ELLIP *eli, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, eli->surface);
  pack_int (isize, i, ints, eli->volume);

  pack_doubles (dsize, d, doubles, eli->cur_center, 3);
  pack_doubles (dsize, d, doubles, (double*)eli->cur_point, 9);

  pack_doubles (dsize, d, doubles, eli->ref_center, 3);
  pack_doubles (dsize, d, doubles, (double*)eli->ref_point, 9);

  pack_doubles (dsize, d, doubles, eli->ref_sca, 3);
  pack_doubles (dsize, d, doubles, eli->ref_rot, 9);

  pack_doubles (dsize, d, doubles, eli->cur_sca, 3);
  pack_doubles (dsize, d, doubles, eli->cur_rot, 9);

  pack_int (isize, i, ints, eli->mat ? 1 : 0); /* pack material existence flag */
  if (eli->mat) pack_string (isize, i, ints, eli->mat->label);
}

/* unpack ellipsoid from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
ELLIP* ELLIP_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  ELLIP *eli;
  int j;

  ERRMEM (eli = MEM_CALLOC (sizeof (ELLIP)));

  eli->surface = unpack_int (ipos, i, ints);
  eli->volume = unpack_int (ipos, i, ints);

  unpack_doubles (dpos, d, doubles, eli->cur_center, 3);
  unpack_doubles (dpos, d, doubles, (double*)eli->cur_point, 9);

  unpack_doubles (dpos, d, doubles, eli->ref_center, 3);
  unpack_doubles (dpos, d, doubles, (double*)eli->ref_point, 9);

  unpack_doubles (dpos, d, doubles, eli->ref_sca, 3);
  unpack_doubles (dpos, d, doubles, eli->ref_rot, 9);

  unpack_doubles (dpos, d, doubles, eli->cur_sca, 3);
  unpack_doubles (dpos, d, doubles, eli->cur_rot, 9);

  j = unpack_int (ipos, i, ints); /* unpack material existence flag */

  if (j)
  {
    SOLFEC *sol = solfec;
    char *label = unpack_string (ipos, i, ints);
    ASSERT_DEBUG_EXT (eli->mat = MATSET_Find (sol->mat, label), "Failed to find material when unpacking a eliere");
    free (label);
  }

  return eli;
}

/* export MBFCP definition */
void ELLIP_2_MBFCP (ELLIP *eli, FILE *out)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED); /* TODO => ELLIP */
}
