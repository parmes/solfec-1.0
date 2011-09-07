/*
 * sph.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * spheres
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
#include "sph.h"
#include "alg.h"
#include "pck.h"
#include "err.h"

/* point in sphere test */
inline static int point_inside (double *center, double radius, double *point)
{
  double v [3];

  SUB (point, center, v);
  if (LEN (v) <= radius + GEOMETRIC_EPSILON) return 1;
  else return 0;
}

/* create a sphere */
SPHERE* SPHERE_Create (double *center, double radius, int surface, int volume)
{
  SPHERE *out;
  double (*rp) [3],
	 (*cp) [3];

  ERRMEM (out = malloc (sizeof (SPHERE)));
  COPY (center, out->ref_center);
  out->ref_radius = radius;
  COPY (center, out->cur_center);
  out->cur_radius = radius;
  out->surface = surface;
  out->volume = volume;
  rp = out->ref_point;
  cp = out->cur_point;
  rp [0][0] = radius;
  rp [0][1] = 0.0;
  rp [0][2] = 0.0;
  ADD (rp[0], center, rp[0]);
  COPY (rp[0], cp[0]);
  rp [1][0] = 0.0;
  rp [1][1] = radius;
  rp [1][2] = 0.0;
  ADD (rp[1], center, rp[1]);
  COPY (rp[1], cp[1]);
  rp [2][0] = 0.0;
  rp [2][1] = 0.0;
  rp [2][2] = radius;
  ADD (rp[2], center, rp[2]);
  COPY (rp[2], cp[2]);
  out->mat = NULL;

  return out;
}

/* dummy (needed in shp.c) */
void SPHERE_Update_Adjacency (SPHERE *sph)
{
}


/* dummy (needed in shp.c) */
int SPHERE_Break_Adjacency (SPHERE *sph, double *point, double *normal)
{
  return 0;
}

/* create a copy of a sphere */
SPHERE* SPHERE_Copy (SPHERE *sph)
{
  SPHERE *twin;

  ERRMEM (twin = malloc (sizeof (SPHERE)));
  memcpy (twin, sph, sizeof (SPHERE));

  return twin;
}

/* scaling of a sphere */
void SPHERE_Scale (SPHERE *sph, double *vector)
{
  double (*ref_pnt) [3] = sph->ref_point,
	 (*cur_pnt) [3] = sph->cur_point,
	 omega [3];

  sph->cur_radius *= vector [0];
  sph->ref_radius = sph->cur_radius;

  for (int i = 0; i < 3; i ++)
  {
    SUB (cur_pnt [i], sph->cur_center, omega);
    SCALE (omega, vector [0]);
    ADD (sph->cur_center, omega, cur_pnt [i]);
    COPY (cur_pnt [i], ref_pnt [i]);
  }
}

/* translation of a sphere */
void SPHERE_Translate (SPHERE *sph, double *vector)
{
  double (*ref_pnt) [3] = sph->ref_point,
	 (*cur_pnt) [3] = sph->cur_point;

  ADD (sph->cur_center, vector, sph->cur_center);
  COPY (sph->cur_center, sph->ref_center);

  for (int i = 0; i < 3; i ++)
  {
    ADD (cur_pnt [i], vector, cur_pnt [i]);
    COPY (cur_pnt [i], ref_pnt [i]);
  }
}

/* rotation of a sphere */
void SPHERE_Rotate (SPHERE *sph, double *point, double *vector, double angle)
{
  double R [9], omega [3];
  double (*ref_pnt) [3] = sph->ref_point,
	 (*cur_pnt) [3] = sph->cur_point;

  angle *=  ALG_PI / 180.0;
  COPY (vector, omega); 
  NORMALIZE (omega); 
  SCALE (omega, angle);
  EXPMAP (omega, R);
  SUB (sph->cur_center, point, omega);
  NVADDMUL (point, R, omega, sph->cur_center);
  COPY (sph->cur_center, sph->ref_center);

  for (int i = 0; i < 3; i ++)
  {
    SUB (cur_pnt [i], point, omega);
    NVADDMUL (point, R, omega, cur_pnt [i]);
    COPY (cur_pnt [i], ref_pnt [i]);
  }
}

/* cut through sphere with a plane; return triangulated cross-section; vertices in the triangles
 * point to the memory allocated after the triangles memory; adjacency is not maintained */
TRI* SPHERE_Cut (SPHERE *sph, double *point, double *normal, int *m)
{
  /* TODO */
  WARNING_DEBUG (0, "Sphere cutting has not been implemented yet!");
  *m = 0;
  return NULL;
}

/* split sphere in two with plane defined by (point, normal); surfid corresponds to the new surface;
 * topoadj != 0 implies cutting from the point and through the topological adjacency only */
void SPHERE_Split (SPHERE *sph, double *point, double *normal, short topoadj, int surfid, CONVEX **one, CONVEX **two)
{
  /* TODO */
  WARNING_DEBUG (0, "Sphere splitting has not been implemented yet!");
  *one = *two = NULL;
}

/* compute partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void SPHERE_Char_Partial (SPHERE *sph, int ref, double *vo, double *sx, double *sy, double *sz, double *eul)
{
  double v, e, r, *a, tmp [9], eye [9];

  r = ref ? sph->ref_radius : sph->cur_radius;
  a = ref ? sph->ref_center : sph->cur_center;
  v = (4.0/3.0)*ALG_PI*r*r*r;

  /* Stainer's theorem =>
   * J(i, j) = I(i, j) + M*(a^2 * delta (i, j) - diadic (a, a))
   * hence it is possible to produce inertia tensor with respect
   * to the origin from one with respect to the center point */

  e = DOT (a, a);
  IDENTITY (eye);
  SCALEDIAG (eye, e);
  DIADIC (a, a, tmp);
  NNSUB (eye, tmp, tmp);
  SCALE9 (tmp, v);
  e = (2.0/5.0)*v*r*r; /* diagonal entry of sphere inertia tensor */
  IDENTITY (eye);
  SCALEDIAG (eye, e);
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
  *sx += v * a [0];
  *sy += v * a [1];
  *sz += v * a [2];
  NNADD (tmp, eul, eul);

  /* TODO: make sure the above is correct */
}

/* get characteristics of a sphere: volume, mass center, and Euler tensor (centered) */
void SPHERE_Char (SPHERE *sph, int ref, double *volume, double *center, double *euler)
{
  double vo, sx, sy, sz, cen [3], eul [9];

  vo = sx = sy = sz = 0.0;
  SET9 (eul, 0.0);

  SPHERE_Char_Partial (sph, ref, &vo, &sx, &sy, &sz, eul);

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

/* update extents of an individual sphere */
void SPHERE_Extents (void *data, SPHERE *sph, double *extents)
{
  double *center, radius;

  radius = sph->cur_radius;
  center = sph->cur_center;
  extents [0] = center [0] - radius - GEOMETRIC_EPSILON;
  extents [1] = center [1] - radius - GEOMETRIC_EPSILON;
  extents [2] = center [2] - radius - GEOMETRIC_EPSILON;
  extents [3] = center [0] + radius + GEOMETRIC_EPSILON;
  extents [4] = center [1] + radius + GEOMETRIC_EPSILON;
  extents [5] = center [2] + radius + GEOMETRIC_EPSILON;
}

/* compute extents of a sphere */
void SPHERE_Extents_2 (SPHERE *sph, double *extents)
{
  SPHERE_Extents (NULL, sph, extents);
}

/* compute oriented extents of a sphere */
void SPHERE_Oriented_Extents (SPHERE *sph, double *vx, double *vy, double *vz, double *extents)
{
  double e [6], len [3], r;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;

  len [0] = LEN (vx);
  len [1] = LEN (vy);
  len [2] = LEN (vz);
    
  e [0] = DOT (sph->cur_center, vx);
  e [1] = DOT (sph->cur_center, vy);
  e [2] = DOT (sph->cur_center, vz);
  COPY (e, e + 3);
  r = sph->cur_radius;
  e [0] -= r / len [0];
  e [1] -= r / len [1];
  e [2] -= r / len [2];
  e [4] += r / len [0];
  e [5] += r / len [1];
  e [6] += r / len [2];

  if (e [0] < extents [0]) extents [0] = e [0];
  if (e [1] < extents [1]) extents [1] = e [1];
  if (e [2] < extents [2]) extents [2] = e [2];
  if (e [3] > extents [3]) extents [3] = e [3];
  if (e [4] > extents [4]) extents [4] = e [4];
  if (e [5] > extents [5]) extents [5] = e [5];
}

/* return first not NULL bulk material for a sphere */
void* SPHERE_First_Bulk_Material (SPHERE *sph)
{
  return sph->mat;
}

/* return sphere containing a spatial point */
SPHERE* SPHERE_Containing_Point (SPHERE *sph, double *point)
{
  if (point_inside (sph->cur_center, sph->cur_radius, point)) return sph;
  else return NULL;
}

/* does this sphere contain the point? */
int SPHERE_Contains_Point (void *dummy, SPHERE *sph, double *point)
{
  return point_inside (sph->cur_center, sph->cur_radius, point);
}

/* return distance of a spatial point to the sphere */
double SPHERE_Spatial_Point_Distance (void *dummy, SPHERE *sph, double *point)
{
  double v [3], d;

  SUB (point, sph->cur_center, v);
  d = LEN (v) - sph->cur_radius;
  return MIN (0.0, d);
}

/* update sphere according to the given motion */
void SPHERE_Update (SPHERE *sph, void *body, void *shp, MOTION motion)
{
  SGP sgp = {shp, sph, GOBJ_SPHERE, NULL};
  double *ref = sph->ref_center,
	 (*ref_pnt) [3] = sph->ref_point,
	 *cur = sph->cur_center,
	 (*cur_pnt) [3] = sph->cur_point;

  if (motion)
  { 
    motion (body, &sgp, ref, cur); /* move center */
    motion (body, &sgp, ref_pnt [0], cur_pnt [0]); /* move marker points */
    motion (body, &sgp, ref_pnt [1], cur_pnt [1]);
    motion (body, &sgp, ref_pnt [2], cur_pnt [2]);
  }
  else
  {
    COPY (ref, cur);
    COPY (ref_pnt [0], cur_pnt [0]);
    COPY (ref_pnt [1], cur_pnt [1]);
    COPY (ref_pnt [2], cur_pnt [2]);
  }
}

/* free sphere */
void SPHERE_Destroy (SPHERE *sph)
{
  free (sph);
}

/* pack sphere into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void SPHERE_Pack (SPHERE *sph, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, sph->surface);
  pack_int (isize, i, ints, sph->volume);

  pack_doubles (dsize, d, doubles, sph->cur_center, 3);
  pack_doubles (dsize, d, doubles, (double*)sph->cur_point, 9);
  pack_double  (dsize, d, doubles, sph->cur_radius);

  pack_doubles (dsize, d, doubles, sph->ref_center, 3);
  pack_doubles (dsize, d, doubles, (double*)sph->ref_point, 9);
  pack_double  (dsize, d, doubles, sph->ref_radius);

  pack_int (isize, i, ints, sph->mat ? 1 : 0); /* pack material existence flag */
  if (sph->mat) pack_string (isize, i, ints, sph->mat->label);
}

/* unpack sphere from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
SPHERE* SPHERE_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  SPHERE *sph;
  int j;

  ERRMEM (sph = MEM_CALLOC (sizeof (SPHERE)));

  sph->surface = unpack_int (ipos, i, ints);
  sph->volume = unpack_int (ipos, i, ints);

  unpack_doubles (dpos, d, doubles, sph->cur_center, 3);
  unpack_doubles (dpos, d, doubles, (double*)sph->cur_point, 9);
  sph->cur_radius = unpack_double  (dpos, d, doubles);

  unpack_doubles (dpos, d, doubles, sph->ref_center, 3);
  unpack_doubles (dpos, d, doubles, (double*)sph->ref_point, 9);
  sph->ref_radius = unpack_double  (dpos, d, doubles);

  j = unpack_int (ipos, i, ints); /* unpack material existence flag */

  if (j)
  {
    SOLFEC *sol = solfec;
    char *label = unpack_string (ipos, i, ints);
    ASSERT_DEBUG_EXT (sph->mat = MATSET_Find (sol->mat, label), "Failed to find material when unpacking a sphere");
    free (label);
  }

  return sph;
}

/* export MBFCP definition */
void SPHERE_2_MBFCP (SPHERE *sph, FILE *out)
{
  fprintf (out, "CENTER:\t%g  %g  %g\n", sph->ref_center [0], sph->ref_center [1], sph->ref_center [2]);
  fprintf (out, "RADIUS:\t%g\n", sph->ref_radius);
  fprintf (out, "SURFID:\t%d\n", sph->surface);
}
