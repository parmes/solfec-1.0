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
#include "spx.h"
#include "alg.h"
#include "mem.h"
#include "map.h"
#include "hyb.h"
#include "gjk.h"
#include "err.h"
#include "pck.h"

/* point in sphere test */
inline static int point_inside (double *center, double radius, double *point)
{
  double v [3];

  SUB (point, center, v);
  if (LEN (v) <= radius + GEOMETRIC_EPSILON) return 1;
  else return 0;
}

/* overlap callback for convex adjacency */
static void overlap (void *data, BOX *one, BOX *two)
{
  double p [3], q [3];
  SPHERE *sph = (SPHERE*)one->sgp,
	 *spg = (SPHERE*)two->sgp;

  if (gjk_sphere_sphere (sph->cur_center, sph->cur_radius, spg->cur_center, spg->cur_radius, p, q) < GEOMETRIC_EPSILON) /* if they touch */
  {
    ERRMEM (sph->adj = realloc (sph->adj, (++sph->nadj) * sizeof (BOX*)));  /* extend adjacency */
    sph->adj [sph->nadj-1] = spg;
    ERRMEM (spg->adj = realloc (spg->adj, (++spg->nadj) * sizeof (BOX*))); 
    spg->adj [spg->nadj-1] = sph;
  }
}

/* create a sphere (sph == NULL) or append spheres list with another sphere */
SPHERE* SPHERE_Create (SPHERE *sph, double *center, double radius, int surface, int volume)
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
  out->adj = NULL;
  out->next = sph;
  rp = out->ref_points;
  cp = out->cur_points;
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

/* glue two sphere lists */
SPHERE* SPHERE_Glue (SPHERE *sph, SPHERE *spg)
{
  SPHERE *out = sph;

  for (; sph->next; sph = sph->next);
  sph->next = spg;
  return out;
}

/* update adjacency data of spheres;
 * no other function affects the adjacency */
void SPHERE_Update_Adjacency (SPHERE *sph)
{
  MEM mem;
  BOX **boxes;
  SPHERE *spg;
  int num;
  
  for (spg = sph, num = 0; spg; spg = spg->next)
  {
    spg->nadj = 0;
    num ++;
  }

  if (num < 2) return;

  MEM_Init (&mem, sizeof (BOX), num);
  ERRMEM (boxes = malloc (sizeof (AABB*) * num));
  for (spg = sph, num = 0; spg; spg = spg->next, num ++)
  {
    ERRMEM (boxes [num] = MEM_Alloc (&mem));
    SPHERE_Extents (NULL, spg, boxes [num]->extents); /* set up extents */
    boxes [num]->sgp = (SGP*)spg;
  }

  hybrid (boxes, num, NULL, overlap); /* detect boxoverlaps => set adjacency inside the callback */

  MEM_Release (&mem); /* done */
  free (boxes);
}

/* create a copy of a list */
SPHERE* SPHERE_Copy (SPHERE *sph)
{
  SPHERE *twin, *tail;

  for (tail = NULL; sph; sph = sph->next)
  {
    ERRMEM (twin = malloc (sizeof (SPHERE)));
    memcpy (twin, sph, sizeof (SPHERE));
    twin->adj = NULL;
    twin->nadj = 0;
    twin->next = tail;
    tail = twin;
  }

  return twin;
}

/* scaling of a list  */
void SPHERE_Scale (SPHERE *sph, double *vector)
{
  double (*ref_pnt) [3] = sph->ref_points,
	 (*cur_pnt) [3] = sph->cur_points,
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

/* translation of a list */
void SPHERE_Translate (SPHERE *sph, double *vector)
{
  double (*ref_pnt) [3] = sph->ref_points,
	 (*cur_pnt) [3] = sph->cur_points;

  ADD (sph->cur_center, vector, sph->cur_center);
  COPY (sph->cur_center, sph->ref_center);

  for (int i = 0; i < 3; i ++)
  {
    ADD (cur_pnt [i], vector, cur_pnt [i]);
    COPY (cur_pnt [i], ref_pnt [i]);
  }
}

/* rotation of a list */
void SPHERE_Rotate (SPHERE *sph, double *point, double *vector, double angle)
{
  double R [9], omega [3];
  double (*ref_pnt) [3] = sph->ref_points,
	 (*cur_pnt) [3] = sph->cur_points;

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

/* cut through spheres with a plane; return triangulated cross-section; vertices in the triangles
 * point to the memory allocated after the triangles memory; adjacency is not maintained;
 * TRI->adj[0] stores a pointer to the geometrical object that has been cut by the triangle */
TRI* SPHERE_Cut (SPHERE *sph, double *point, double *normal, int *m)
{
  WARNING_DEBUG (0, "Sphere cutting has not been implemented yet");
  /* TODO */
  *m = 0;
  return NULL;
}

/* split sphere lists in two lists with plane defined by (point, normal); adjacencies between
 * the split lists elements need to be recomputed; surfid corresponds to the new surface */
void SPHERE_Split (SPHERE *sph, double *point, double *normal, int surfid, SPHERE **one, SPHERE **two)
{
  double v [3];
  SPHERE *o;

  *one = *two = NULL;

  for (; sph; sph = sph->next)
  {
    o = SPHERE_Copy (sph);

    SUB (o->cur_center, point, v);

    if (DOT (v, normal) >= 0)
    {
      o->next = *two;
      *two = o;
    }
    else
    {
      o->next = *one;
      *one = o;
    }
  }
}

/* compute current partial characteristic: 'vo'lume and static momenta
 * 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void SPHERE_Char_Partial (SPHERE *sph, double *vo, double *sx, double *sy, double *sz, double *eul)
{
  double v, e, r, *a, tmp [9], eye [9];

  for (; sph; sph = sph->next)
  {
    r = sph->cur_radius;
    a = sph->cur_center;
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
  }

  /* TODO: make sure the above is correct */
}

/* get 'cur' characteristics of the spheres in list:
 * volume, mass center, and Euler tensor (centered) */
void SPHERE_Char (SPHERE *sph, double *volume, double *center, double *euler)
{
  double vo, sx, sy, sz,
	 cen [3], eul [9];

  vo = sx = sy = sz = 0.0;
  SET9 (eul, 0.0);

  SPHERE_Char_Partial (sph, &vo, &sx, &sy, &sz, eul);

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

/* compute extents of sphere list */
void SPHERE_List_Extents (SPHERE *sph, double *extents)
{
  double e [6];

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;
    
  for (; sph; sph = sph->next)
  {
    SPHERE_Extents (NULL, sph, e);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [3] > extents [3]) extents [3] = e [3];
    if (e [4] > extents [4]) extents [4] = e [4];
    if (e [5] > extents [5]) extents [5] = e [5];
  }
}

/* compute oriented extents of sphere list */
void SPHERE_List_Oriented_Extents (SPHERE *sph, double *vx, double *vy, double *vz, double *extents)
{
  double e [6], len [3], r;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;

  len [0] = LEN (vx);
  len [1] = LEN (vy);
  len [2] = LEN (vz);
    
  for (; sph; sph = sph->next)
  {
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
}

/* return first not NULL bulk material for a sphere list */
void* SPHERE_First_Bulk_Material (SPHERE *sph)
{
  for (; sph; sph = sph->next)
    if (sph->mat) return sph->mat;

  return NULL;
}

/* return sphere containing a spatial point */
SPHERE* SPHERE_Containing_Point (SPHERE *sph, double *point)
{
  for (; sph; sph = sph->next)
    if (point_inside (sph->cur_center, sph->cur_radius, point)) return sph;

  return NULL;
}

/* does this sphere (not a list) contain the point? */
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

/* update sphere list according to the given motion */
void SPHERE_Update (SPHERE *sph, void *body, void *shp, MOTION motion)
{
  for (; sph; sph = sph->next)
  {
    double *ref = sph->ref_center,
	   (*ref_pnt) [3] = sph->ref_points,
	   *cur = sph->cur_center,
	   (*cur_pnt) [3] = sph->cur_points;

    if (motion) motion (body, shp, sph, ref, cur); /* move center */
    else { COPY (ref, cur); }

    if (motion)
    {
      motion (body, shp, sph, ref_pnt [0], cur_pnt [0]); /* move marker points */
      motion (body, shp, sph, ref_pnt [1], cur_pnt [1]);
      motion (body, shp, sph, ref_pnt [2], cur_pnt [2]);
    }
    else
    {
      COPY (ref_pnt [0], cur_pnt [0]);
      COPY (ref_pnt [1], cur_pnt [1]);
      COPY (ref_pnt [2], cur_pnt [2]);
    }
  }
}

/* test wether two spheres are adjacent */
int SPHERE_Adjacent (SPHERE *one, SPHERE *two)
{
  int n;

  for (n = 0; n < two->nadj; n ++)
    if (two->adj [n] == one) return 1; /* enough to compare one way (adjacency lists are symmetric) */

  return 0;
}

/* free list of spheres */
void SPHERE_Destroy (SPHERE *sph)
{
  SPHERE *nxt;

  for (;sph; sph = nxt)
  {
    nxt = sph->next;
    free (sph->adj);
    free (sph);
  }
}

void SPHERE_Pack (SPHERE *sph, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SPHERE *ptr;
  int count, n;
  MAP *map;

  for (count = 0, ptr = sph; ptr; ptr = ptr->next) count ++; /* number of spheres */

  for (map = NULL, n = 0, ptr = sph; ptr; ptr = ptr->next, n ++)
    MAP_Insert (NULL, &map, ptr, (void*) (long) n, NULL); /* map pointers to table indices */

  pack_int (isize, i, ints, count); /* number of spheres */

  for (; sph; sph = sph->next)
  {
    pack_int (isize, i, ints, sph->surface);
    pack_int (isize, i, ints, sph->nadj);
    pack_int (isize, i, ints, sph->volume);

    pack_doubles (dsize, d, doubles, sph->cur_center, 3);
    pack_doubles (dsize, d, doubles, (double*)sph->cur_points, 9);
    pack_double  (dsize, d, doubles, sph->cur_radius);

    pack_doubles (dsize, d, doubles, sph->ref_center, 3);
    pack_doubles (dsize, d, doubles, (double*)sph->ref_points, 9);
    pack_double  (dsize, d, doubles, sph->ref_radius);

    /* rather than adjacency pack indices of neighbours in the output sequence */
    for (n = 0; n < sph->nadj; n ++) pack_int (isize, i, ints, (int) (long) MAP_Find (map, sph->adj [n], NULL));

    pack_int (isize, i, ints, sph->mat ? 1 : 0); /* pack material existence flag */
    if (sph->mat) pack_string (isize, i, ints, sph->mat->label);
  }

  MAP_Free (NULL, &map);
}

SPHERE* SPHERE_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int count, k, surface, nadj, volume, n, j;
  SPHERE *ptr, **tab, *tail = NULL, *head;
  
  count = unpack_int (ipos, i, ints); /* number of spheres */

  ERRMEM (tab = malloc (count * sizeof (SPHERE*)));
  
  for (k = 0; k < count; k ++) /* unpack spheres */
  {
    surface = unpack_int (ipos, i, ints);
    nadj = unpack_int (ipos, i, ints);
    volume = unpack_int (ipos, i, ints);

    ERRMEM (ptr = MEM_CALLOC (sizeof (SPHERE)));
    ERRMEM (ptr->adj = malloc (nadj * sizeof (SPHERE*)));
    ptr->surface = surface;
    ptr->nadj = 0;
    ptr->volume = volume;
    tab [k] = ptr;

    if (tail)
    {
      tail->next = ptr;
      tail = ptr;
    }
    else head = tail = ptr;

    unpack_doubles (dpos, d, doubles, ptr->cur_center, 3);
    unpack_doubles (dpos, d, doubles, (double*)ptr->cur_points, 9);
    ptr->cur_radius = unpack_double  (dpos, d, doubles);

    unpack_doubles (dpos, d, doubles, ptr->ref_center, 3);
    unpack_doubles (dpos, d, doubles, (double*)ptr->ref_points, 9);
    ptr->ref_radius = unpack_double  (dpos, d, doubles);

    for (n = 0; n < nadj; n ++)
    {
      j = unpack_int (ipos, i, ints);
      ptr->adj [ptr->nadj ++] = (SPHERE*) (long) j; /* store index in 'tab' for the moment */
    }

    j = unpack_int (ipos, i, ints); /* unpack material existence flag */

    if (j)
    {
      SOLFEC *sol = solfec;
      char *label = unpack_string (ipos, i, ints);
      ASSERT_DEBUG_EXT (ptr->mat = MATSET_Find (sol->mat, label), "Failed to find material when unpacking a sphere");
      free (label);
    }
  }

  /* now map adjacency */
  for (k = 0; k < count; k ++)
  {
    ptr = tab [k];

    for (n = 0; n < ptr->nadj; n ++)
    {
      ASSERT_DEBUG (0 <= (int) (long) ptr->adj [n] && (int) (long) ptr->adj [n] < count, "Adjacent sphere index out of bounds");
      ptr->adj [n] = tab [(int) (long) ptr->adj [n]];
    }
  }

  free (tab);

  return head;
}

/* export MBFCP definition */
void SPHERE_2_MBFCP (SPHERE *sph, FILE *out)
{
  SPHERE *ptr;
  int n;

  for (ptr = sph, n = 0; ptr; ptr = ptr->next, n ++);

  fprintf (out, "SPHERES:\t%d\n", n);

  for (ptr = sph; ptr; ptr = ptr->next)
  {
    fprintf (out, "CENTER:\t%g  %g  %g\n", ptr->ref_center [0], ptr->ref_center [1], ptr->ref_center [2]);
    fprintf (out, "RADIUS:\t%g\n", ptr->ref_radius);
    fprintf (out, "SURFID:\t%d\n", ptr->surface);
  }
}
