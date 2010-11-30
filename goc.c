/*
 * goc.c
 * Copyright (C) 2008, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * geometric object contact detection
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

#include <limits.h>
#include <float.h>
#include "alg.h"
#include "box.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "cvi.h"
#include "gjk.h"
#include "goc.h"
#include "err.h"

/* line-plane intersection => intersection = point + direction * coef if 1 is returned */
inline static int lineplane (double *plane, double *point, double *direction, double *coef)
{
  double d;
 
  d = DOT (plane, direction); 

  if (fabs (d) < GEOMETRIC_EPSILON) return 0; /* the penalty solver is sensitice to ill-computed gaps, hece our caution */

  *coef = - PLANE (plane, point) / d;

  return 1;
}

/* compute semi-negative convex-sphere gap */
inline static double convex_sphere_gap (double *pla, int npla, double *center, double radius, double *normal)
{
  double plane [4],
         max = -DBL_MAX,
         coef;

  for (; npla > 0; pla += 6, npla --)
  {
    COPY (pla, plane);
    plane [3] = - DOT (plane, pla+3);

    if (lineplane (plane, center, normal, &coef))
    {
      if (coef > max) max = coef;
    }
  }

  max = fabs (max);

  return (radius > max ? max - radius : 0);
}

/* compute semi-negative sphere-sphere gap */
inline static double sphere_sphere_gap (double *ca, double ra, double *cb, double rb, double *normal)
{
  double x [3], d, e;

  SUB (cb, ca, x);
  d = DOT (x, normal);
  e = ra + rb;

  return (e > d ? d - e : 0);
}

/* return a surface code with the input point closest to an input plane */
inline static int nearest_surface (double *pnt, double *pla, int *sur, int n)
{
  double v [3], d, min = DBL_MAX;
  int ret = sur [0];

  for (; n > 0; n --, pla += 6, sur ++)
  {
    SUB (pnt, pla+3, v);
    d = DOT (pla, v);

    if (fabs (d) < min)
    {
      min = d;
      ret = *sur;
    }
  }

  return ret;
}

/* if m > 0 return surface outward normal index (1 or 2);
 * if m < 0 assume outward a-normal and return 1 if spair did not change or 2 otherwise;
 * return 0 if a geometrical error occured */
inline static int point_normal_spair_area_gap
 (TRI *tri, int m, double *pa, int npa, double *pb, int npb, int *sa, int nsa, int *sb, int nsb, /* input */
  double *point, double *normal, int *spair, double *area, double *gap) /* output */
{
  double plane [4], v [3], p [3], pos = DBL_MAX, neg = -pos, a, b;
  TRI *t, *e;
  int j, k;

  SET (normal, 0.0);
  SET (point, 0.0);
  SET (p, 0.0);
  *area = 0.0;
  k = 1;

  for (t = tri, e = t + ABS (m); t != e; t ++)
  {
    TRIANGLE_AREA (t->ver[0], t->ver[1], t->ver[2], a); b = a*a; /* eliminate influence of manute areas */
    if (t->flg > 0 &&  t->flg < nsa + 1) ADDMUL (normal, b, t->out, normal); /* add from one surface */
    if (t->flg < 0 && -t->flg < nsb + 1) SUBMUL (normal, b, t->out, normal); /* subtract from another */
    MID3 (t->ver[0], t->ver[1], t->ver[2], v);
    ADDMUL (point, a, v, point);
    ACC (v, p);
    *area += a;
  }
  DIV (point, *area, point); /* surface mass center */
  NORMALIZE (normal); /* resultant normal */
  *area *= 0.5; /* half surface */

  for (t = tri; t != e; t ++) /* test whether point falls outside of intersection volume */
  {
    SUB (t->ver [0], point, v);
    if (DOT (t->out, v) < 0.0) /* it does */
    {
      a = (double) ABS (m);
      DIV (p, a, point); /* use average point instead */
      break;
    }
  }

  if (m > 0)
  {
    spair [0] = nearest_surface (point, pa, sa, nsa);
    spair [1] = nearest_surface (point, pb, sb, nsb);
  }
  else
  {
    if ((j = nearest_surface (point, pa, sa, nsa)) != spair [0]) spair [0] = j, k = 2;
    if ((j = nearest_surface (point, pb, sb, nsb)) != spair [1]) spair [1] = j, k = 2;
  }

  for (t = tri; t != e; t ++) /* compute semi-negative penetration gap from surface triangles */
  {
    if ((t->flg > 0 &&  t->flg < nsa + 1) ||
	(t->flg < 0 && -t->flg < nsb + 1)) /* belonging to outer surfaces */
    {
      COPY (t->out, plane);
      plane [3] = - DOT (plane, t->ver [0]);

      if (lineplane (plane, point, normal, &a))
      {
	if (a >= 0.0 && a < pos) pos = a;
	else if (a <= 0.0 && a > neg) neg = a;
      }
    }
  }
  if (pos == DBL_MAX || neg == -DBL_MAX) *gap = 0.0;
  else *gap = neg - pos;

  return k;
}

/* detect contact between two convex polyhedrons 'a' and 'b', where
 * va, nva, vb, nbv: vertices and vertex counts
 * pa, npa, pb, npb: 6-coord planes (normal, point) and plane counts
 * sa, nsa, sb, nsb: surface identifiers and surface planes counts
 *                   (the first nsa, nsb planes are on the surface)
 * ----------------------------------------------------------------
 * ramaining paramteres behave like in gobjcontact routine 
 */
static int detect_convex_convex (
  double *va, int nva, double *pa, int npa, int *sa, int nsa,
  double *vb, int nvb, double *pb, int npb, int *sb, int nsb,
  double onepnt [3],
  double twopnt [3],
  double normal [3],
  double *gap,
  double *area,
  int spair [2])
{
  double sanity;
  int k = 0, m;
  TRI *tri;

  sanity = gjk (va, nva, vb, nvb, onepnt, twopnt);

  if (sanity < GEOMETRIC_EPSILON)
  {
    if (!(tri = cvi (va, nva, pa, npa, vb, nvb, pb, npb, NON_REGULARIZED, &m))) return 0;

    k = point_normal_spair_area_gap (tri, m, pa, npa, pb, npb, sa, nsa, sb, nsb, onepnt, normal, spair, area, gap);
    sanity = (onepnt[0]+onepnt[1]+onepnt[2]+normal[0]+normal[1]+normal[2]+(*area)+(*gap));
    COPY (onepnt, twopnt);
    free (tri);
  }

  if (!isfinite (sanity)) return 0;
  else return k;
}

/* detect contact between a convex and a sphere */
static int detect_convex_sphere (
  double *vc, int nvc, double *pc, int npc, int *sc, int nsc, /* same as above */
  double *c, double r, int s, /* center, radius, surface */
  double onepnt [3],
  double twopnt [3],
  double normal [3],
  double *gap,
  double *area,
  int spair [2])
{
  double g, h;

  if ((g = gjk_convex_sphere (vc, nvc, c, r, onepnt, twopnt)) < GEOMETRIC_EPSILON)
  {
    for (h = r; g < GEOMETRIC_EPSILON && h > GEOMETRIC_EPSILON; h *= 0.5)
      g = gjk_convex_sphere (vc, nvc, c, h, onepnt, twopnt); /* shrink sphere and redo => find projection
								of the sphere center on the convex surface */
    SUB (c, onepnt, normal);
    NORMALIZE (normal);
    ADDMUL (c, -r, normal, twopnt);

    spair [0] = nearest_surface (onepnt, pc, sc, nsc);
    spair [1] = s;
    *area = 1.0;
    *gap = convex_sphere_gap (pc, npc, c, r, normal);
    return 1;
  }

  return 0;
}

/* compare two points */
inline static int pntcmp (double *a, double *b)
{
  for (int i = 0; i < 3 ; i ++)
  {
    if (a [i] < b [i]) return -1;
    else if (a [i] > b [i]) return 1;
  }

  return 0;
}

/* detect contact between spheres 'a' and 'b' */
static int detect_sphere_sphere (
  double *ca, double ra, int sa, /* center, radius, surface */
  double *cb, double rb, int sb,
  double onepnt [3],
  double twopnt [3],
  double normal [3],
  double *gap,
  double *area,
  int spair [2])
{
  if (((*gap) = gjk_sphere_sphere (ca, ra, cb, rb, onepnt, twopnt)) < GEOMETRIC_EPSILON)
  {
    spair [0] = sa;
    spair [1] = sb;

    *area = 1.0;

    if (pntcmp (ca, cb) <= 0) /* same normal orientation regardless of sphere processing order */
    {
      SUB (onepnt, ca, normal);
      NORMALIZE (normal);
      *gap = sphere_sphere_gap (ca, ra, cb, rb, normal);
      return 1;
    }
    else
    {
      SUB (twopnt, cb, normal);
      NORMALIZE (normal);
      *gap = sphere_sphere_gap (cb, rb, ca, ra, normal);
      return 2;
    }
  }

  return 0;
}

/* update contact between two convex polyhedrons 'a' and 'b', where
 * va, nva, vb, nbv: vertices and vertex counts
 * pa, npa, pb, npb: 6-coord planes (normal, point) and plane counts
 * sa, nsa, sb, nsb: surface identifiers and surface planes counts
 *                   (the first nsa, nsb planes are on the surface)
 * ----------------------------------------------------------------
 * ramaining paramteres behave like in gobjcontact routine 
 */
static int update_convex_convex (
  double *va, int nva, double *pa, int npa, int *sa, int nsa,
  double *vb, int nvb, double *pb, int npb, int *sb, int nsb,
  double onepnt [3],
  double twopnt [3],
  double normal [3], /* outward with restpect to the 'a' body (master) */
  double *gap,
  double *area,
  int spair [2])
{
  double sanity;
  int k = 0, m;
  TRI *tri;

  sanity = gjk (va, nva, vb, nvb, onepnt, twopnt);

  if (sanity < GEOMETRIC_EPSILON)
  {
    if (!(tri = cvi (va, nva, pa, npa, vb, nvb, pb, npb, NON_REGULARIZED, &m))) return 0;

    k = point_normal_spair_area_gap (tri, -m, pa, npa, pb, npb, sa, nsa, sb, nsb, onepnt, normal, spair, area, gap);
    sanity = (onepnt[0]+onepnt[1]+onepnt[2]+normal[0]+normal[1]+normal[2]+(*area)+(*gap));
    COPY (onepnt, twopnt);
    free (tri);
  }

  if (!isfinite (sanity)) return 0;
  else return k;
}

/* update contact between a convex and a sphere */
static int update_convex_sphere (
  double *vc, int nvc, double *pc, int npc, int *sc, int nsc, /* same as above */
  double *c, double r, int s, /* center, radius, surface */
  double onepnt [3],
  double twopnt [3],
  double normal [3],
  double *gap,
  double *area,
  int spair [2])
{
  double g, h;
  int s0;

  if (((*gap) = gjk_convex_sphere (vc, nvc, c, r, onepnt, twopnt)) < GEOMETRIC_EPSILON)
  {
    for (h = r, g = *gap; g < GEOMETRIC_EPSILON && h > GEOMETRIC_EPSILON; h *= 0.5)
      g = gjk_convex_sphere (vc, nvc, c, h, onepnt, twopnt); /* shrink sphere and redo => find projection
								of the sphere center on the convex surface */
    SUB (c, onepnt, normal);
    NORMALIZE (normal);
    ADDMUL (c, -r, normal, twopnt);

    s0 = spair [0];
    spair [0] = nearest_surface (onepnt, pc, sc, nsc);
    *gap = convex_sphere_gap (pc, npc, c, r, normal);

    if (s0 == spair [0]) return 1;
    else return 2;
  }

  return 0;
}

/* update contact between spheres 'a' and 'b' */
static int update_sphere_sphere (
  double *ca, double ra, int sa, /* center, radius, surface */
  double *cb, double rb, int sb,
  double onepnt [3],
  double twopnt [3],
  double normal [3],
  double *gap,
  double *area,
  int spair [2])
{
  if (((*gap) = gjk_sphere_sphere (ca, ra, cb, rb, onepnt, twopnt)) < GEOMETRIC_EPSILON)
  {
    SUB (onepnt, ca, normal);
    NORMALIZE (normal);
    *gap = sphere_sphere_gap (ca, ra, cb, rb, normal);
    return 1;
  }

  return 0;
}

/* initialize convex a convex representation */
inline static void convex_init (CONVEX *cvx, double **v, int *nv, double **p, int *np, int **s, int *ns)
{
  *v = cvx->cur;
  *nv = cvx->nver;
  *p = CONVEX_Planes (cvx);
  *np = cvx->nfac;
  *s = cvx->surface;
  *ns = cvx->nfac;
}

/* finalize a convex representation */
inline static void convex_done (CONVEX *cvx, double **v, int *nv, double **p, int *np, int **s, int *ns)
{
  free (*p);
}

/* swap surface pairs */
inline static void swap (int spair [2])
{
  int tmp = spair [0];

  spair [0] = spair [1];
  spair [1] = tmp;
}

/* swap back surface pair and swap renturn value */
inline static int detect_swap (int ret, int spair [2])
{
  int tmp = spair [0];

  spair [0] = spair [1];
  spair [1] = tmp;

  return (ret == 1 ? 2 : (ret == 2 ? 1 : 0));
}

/* swap back surface pair and copy return */
inline static int update_swap (int ret, int spair [2])
{
  int tmp = spair [0];

  spair [0] = spair [1];
  spair [1] = tmp;

  return ret;
}

/* detect contact */
static int detect (
    short paircode,
    SHAPE *oneshp, void *onegobj,
    SHAPE *twoshp, void *twogobj,
    double onepnt [3],
    double twopnt [3],
    double normal [3],
    double *gap,
    double *area,
    int spair [2])
{
  switch (paircode)
  {
    case AABB_ELEMENT_ELEMENT:
    {
      double va [24], pa [36],
             vb [24], pb [36];
      int nva, npa, sa [6], nsa,
	  nvb, npb, sb [6], nsb;

      nva = ELEMENT_Vertices (oneshp->data, onegobj, va);
      npa = ELEMENT_Planes (oneshp->data, onegobj, pa, sa, &nsa);
      nvb = ELEMENT_Vertices (twoshp->data, twogobj, vb);
      npb = ELEMENT_Planes (twoshp->data, twogobj, pb, sb, &nsb);

      return detect_convex_convex (va, nva, pa, npa, sa, nsa,
                                   vb, nvb, pb, npb, sb, nsb,
                                   onepnt, twopnt, normal,
				   gap, area, spair);
    }
    break;
    case AABB_CONVEX_CONVEX:
    {
      double *va, *pa,
             *vb, *pb;
      int nva, npa, *sa, nsa,
	  nvb, npb, *sb, nsb,
	  ret;

      convex_init (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);
      convex_init (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      ret = detect_convex_convex (va, nva, pa, npa, sa, nsa,
                                  vb, nvb, pb, npb, sb, nsb,
                                  onepnt, twopnt, normal,
				  gap, area, spair);

      convex_done (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);
      convex_done (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      return ret;
    }
    break;
    case AABB_SPHERE_SPHERE:
    {
      SPHERE *a = onegobj,
	     *b = twogobj;

      return detect_sphere_sphere (a->cur_center, a->cur_radius, a->surface,
                                   b->cur_center, b->cur_radius, b->surface,
                                   onepnt, twopnt, normal,
				   gap, area, spair);
    }
    break;
    case AABB_ELEMENT_CONVEX:
    {
      double va [24], pa [36],
             *vb, *pb;
      int nva, npa, sa [6], nsa,
	  nvb, npb, *sb, nsb,
	  ret;

      nva = ELEMENT_Vertices (oneshp->data, onegobj, va);
      npa = ELEMENT_Planes (oneshp->data, onegobj, pa, sa, &nsa);

      convex_init (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      ret = detect_convex_convex (va, nva, pa, npa, sa, nsa,
                                  vb, nvb, pb, npb, sb, nsb,
                                  onepnt, twopnt, normal,
				  gap, area, spair);

      convex_done (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      return ret;
    }
    break;
    case AABB_CONVEX_ELEMENT:
    {
      double *va, *pa,
             vb [24], pb [36];
      int nva, npa, *sa, nsa,
	  nvb, npb, sb [6], nsb,
	  ret;

      convex_init (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);

      nvb = ELEMENT_Vertices (twoshp->data, twogobj, vb);
      npb = ELEMENT_Planes (twoshp->data, twogobj, pb, sb, &nsb);

      ret = detect_convex_convex (va, nva, pa, npa, sa, nsa,
                                  vb, nvb, pb, npb, sb, nsb,
                                  onepnt, twopnt, normal,
				  gap, area, spair);

      convex_done (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);

      return ret;
    }
    break;
    case AABB_ELEMENT_SPHERE:
    {
      double va [24], pa [36];
      int nva, npa, sa [6], nsa;
      SPHERE *b = twogobj;

      nva = ELEMENT_Vertices (oneshp->data, onegobj, va);
      npa = ELEMENT_Planes (oneshp->data, onegobj, pa, sa, &nsa);

      return detect_convex_sphere (va, nva, pa, npa, sa, nsa,
                                   b->cur_center, b->cur_radius, b->surface,
                                   onepnt, twopnt, normal,
				   gap, area, spair);
    }
    break;
    case AABB_SPHERE_ELEMENT:
    {
      SPHERE *a = onegobj;
      double vb [24], pb [36];
      int nvb, npb, sb [6], nsb,
	  ret;

      nvb = ELEMENT_Vertices (twoshp->data, twogobj, vb);
      npb = ELEMENT_Planes (twoshp->data, twogobj, pb, sb, &nsb);

      swap (spair);

      ret = detect_convex_sphere (vb, nvb, pb, npb, sb, nsb,
                                  a->cur_center, a->cur_radius, a->surface,
                                  twopnt, onepnt, normal,
				  gap, area, spair);

      return detect_swap (ret, spair);

    }
    break;
    case AABB_CONVEX_SPHERE:
    {
      double *va, *pa;
      SPHERE *b = twogobj;
      int nva, npa, *sa, nsa,
	  ret;

      convex_init (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);

      ret = detect_convex_sphere (va, nva, pa, npa, sa, nsa,
                                  b->cur_center, b->cur_radius, b->surface,
                                  onepnt, twopnt, normal,
				  gap, area, spair);

      convex_done (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);

      return ret;
    }
    break;
    case AABB_SPHERE_CONVEX:
    {
      SPHERE *a = onegobj;
      double *vb, *pb;
      int nvb, npb, *sb, nsb,
	  ret;

      convex_init (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      swap (spair);

      ret = detect_convex_sphere (vb, nvb, pb, npb, sb, nsb,
                                  a->cur_center, a->cur_radius, a->surface,
                                  twopnt, onepnt, normal,
				  gap, area, spair);

      convex_done (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      return detect_swap (ret, spair);
    }
    break;
  }

  return 0;
}

/* update contact */
static int update (
    short paircode,
    SHAPE *oneshp, void *onegobj,
    SHAPE *twoshp, void *twogobj,
    double onepnt [3],
    double twopnt [3],
    double normal [3],
    double *gap,
    double *area,
    int spair [2])
{
  switch (paircode)
  {
    case AABB_ELEMENT_ELEMENT:
    {
      double va [24], pa [36],
             vb [24], pb [36];
      int nva, npa, sa [6], nsa,
	  nvb, npb, sb [6], nsb;

      nva = ELEMENT_Vertices (oneshp->data, onegobj, va);
      npa = ELEMENT_Planes (oneshp->data, onegobj, pa, sa, &nsa);
      nvb = ELEMENT_Vertices (twoshp->data, twogobj, vb);
      npb = ELEMENT_Planes (twoshp->data, twogobj, pb, sb, &nsb);

      return update_convex_convex (va, nva, pa, npa, sa, nsa,
                                   vb, nvb, pb, npb, sb, nsb,
                                   onepnt, twopnt, normal,
				   gap, area, spair);
    }
    break;
    case AABB_CONVEX_CONVEX:
    {
      double *va, *pa,
             *vb, *pb;
      int nva, npa, *sa, nsa,
	  nvb, npb, *sb, nsb,
	  ret;

      convex_init (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);
      convex_init (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      ret = update_convex_convex (va, nva, pa, npa, sa, nsa,
                                  vb, nvb, pb, npb, sb, nsb,
                                  onepnt, twopnt, normal,
				  gap, area, spair);

      convex_done (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);
      convex_done (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      return ret;
    }
    break;
    case AABB_SPHERE_SPHERE:
    {
      SPHERE *a = onegobj,
	     *b = twogobj;

      return update_sphere_sphere (a->cur_center, a->cur_radius, a->surface,
                                   b->cur_center, b->cur_radius, b->surface,
                                   onepnt, twopnt, normal,
				   gap, area, spair);
    }
    break;
    case AABB_ELEMENT_CONVEX:
    {
      double va [24], pa [36],
             *vb, *pb;
      int nva, npa, sa [6], nsa,
	  nvb, npb, *sb, nsb,
	  ret;

      nva = ELEMENT_Vertices (oneshp->data, onegobj, va);
      npa = ELEMENT_Planes (oneshp->data, onegobj, pa, sa, &nsa);

      convex_init (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      ret = update_convex_convex (va, nva, pa, npa, sa, nsa,
                                  vb, nvb, pb, npb, sb, nsb,
                                  onepnt, twopnt, normal,
				  gap, area, spair);

      convex_done (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      return ret;
    }
    break;
    case AABB_CONVEX_ELEMENT:
    {
      double *va, *pa,
             vb [24], pb [36];
      int nva, npa, *sa, nsa,
	  nvb, npb, sb [6], nsb,
	  ret;

      convex_init (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);

      nvb = ELEMENT_Vertices (twoshp->data, twogobj, vb);
      npb = ELEMENT_Planes (twoshp->data, twogobj, pb, sb, &nsb);

      ret = update_convex_convex (va, nva, pa, npa, sa, nsa,
                                  vb, nvb, pb, npb, sb, nsb,
                                  onepnt, twopnt, normal,
				  gap, area, spair);

      convex_done (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);

      return ret;
    }
    break;
    case AABB_ELEMENT_SPHERE:
    {
      double va [24], pa [36];
      int nva, npa, sa [6], nsa;
      SPHERE *b = twogobj;

      nva = ELEMENT_Vertices (oneshp->data, onegobj, va);
      npa = ELEMENT_Planes (oneshp->data, onegobj, pa, sa, &nsa);

      return update_convex_sphere (va, nva, pa, npa, sa, nsa,
                                   b->cur_center, b->cur_radius, b->surface,
                                   onepnt, twopnt, normal,
				   gap, area, spair);
    }
    break;
    case AABB_SPHERE_ELEMENT:
    {
      SPHERE *a = onegobj;
      double vb [24], pb [36];
      int nvb, npb, sb [6], nsb,
	  ret;

      nvb = ELEMENT_Vertices (twoshp->data, twogobj, vb);
      npb = ELEMENT_Planes (twoshp->data, twogobj, pb, sb, &nsb);

      swap (spair);

      ret = update_convex_sphere (vb, nvb, pb, npb, sb, nsb,
                                  a->cur_center, a->cur_radius, a->surface,
                                  twopnt, onepnt, normal,
				  gap, area, spair);

      return update_swap (ret, spair);

    }
    break;
    case AABB_CONVEX_SPHERE:
    {
      double *va, *pa;
      SPHERE *b = twogobj;
      int nva, npa, *sa, nsa,
	  ret;

      convex_init (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);

      ret = update_convex_sphere (va, nva, pa, npa, sa, nsa,
                                  b->cur_center, b->cur_radius, b->surface,
                                  onepnt, twopnt, normal,
				  gap, area, spair);

      convex_done (onegobj, &va, &nva, &pa, &npa, &sa, &nsa);

      return ret;
    }
    break;
    case AABB_SPHERE_CONVEX:
    {
      SPHERE *a = onegobj;
      double *vb, *pb;
      int nvb, npb, *sb, nsb,
	  ret;

      convex_init (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      swap (spair);

      ret = update_convex_sphere (vb, nvb, pb, npb, sb, nsb,
                                  a->cur_center, a->cur_radius, a->surface,
                                  twopnt, onepnt, normal,
				  gap, area, spair);

      convex_done (twogobj, &vb, &nvb, &pb, &npb, &sb, &nsb);

      return update_swap (ret, spair);
    }
    break;
  }

  return 0;
}

/* detect or update contact data
 * between two geometric objects */
int gobjcontact (
    GOCDO action, short paircode,
    SHAPE *oneshp, void *onegobj,
    SHAPE *twoshp, void *twogobj,
    double onepnt [3],
    double twopnt [3],
    double normal [3],
    double *gap,
    double *area,
    int spair [2])
{
  if (action == CONTACT_DETECT)
    return detect (paircode, oneshp, onegobj, twoshp,
      twogobj, onepnt, twopnt, normal, gap, area, spair);
  else return update (paircode, oneshp, onegobj, twoshp,
    twogobj, onepnt, twopnt, normal, gap, area, spair);
}


/* get distance between two objects (output closest point pair in p, q) */
double gobjdistance (short paircode, SGP *one, SGP *two, double *p, double *q)
{
  switch (paircode)
  {
    case AABB_ELEMENT_ELEMENT:
    {
      double va [24],
             vb [24];
      int nva,
	  nvb;

      nva = ELEMENT_Vertices (one->shp->data, one->gobj, va);
      nvb = ELEMENT_Vertices (two->shp->data, two->gobj, vb);

      return gjk (va, nva, vb, nvb, p, q);
    }
    break;
    case AABB_CONVEX_CONVEX:
    {
      CONVEX *a = one->gobj,
	     *b = two->gobj;

      return gjk (a->cur, a->nver, b->cur, b->nver, p, q);
    }
    break;
    case AABB_SPHERE_SPHERE:
    {
      SPHERE *a = one->gobj,
	     *b = two->gobj;

      return gjk_sphere_sphere (a->cur_center, a->cur_radius, b->cur_center, b->cur_radius, p, q);
    }
    break;
    case AABB_ELEMENT_CONVEX:
    {
      double va [24];
      int nva;
      CONVEX *b = two->gobj;

      nva = ELEMENT_Vertices (one->shp->data, one->gobj, va);

      return gjk (va, nva, b->cur, b->nver, p, q);
    }
    break;
    case AABB_CONVEX_ELEMENT:
    {
      CONVEX *a = one->gobj;
      double vb [24];
      int nvb;

      nvb = ELEMENT_Vertices (two->shp->data, two->gobj, vb);

      return gjk (a->cur, a->nver, vb, nvb, p, q);
    }
    break;
    case AABB_ELEMENT_SPHERE:
    {
      double va [24];
      int nva;
      SPHERE *b = two->gobj;

      nva = ELEMENT_Vertices (one->shp->data, one->gobj, va);

      return gjk_convex_sphere (va, nva, b->cur_center, b->cur_radius, p, q);
    }
    break;
    case AABB_SPHERE_ELEMENT:
    {
      SPHERE *a = one->gobj;
      double vb [24];
      int nvb;

      nvb = ELEMENT_Vertices (two->shp->data, two->gobj, vb);

      return gjk_convex_sphere (vb, nvb, a->cur_center, a->cur_radius, q, p);
    }
    break;
    case AABB_CONVEX_SPHERE:
    {
      CONVEX *a = one->gobj;
      SPHERE *b = two->gobj;

      return gjk_convex_sphere (a->cur, a->nver, b->cur_center, b->cur_radius, p, q);
    }
    break;
    case AABB_SPHERE_CONVEX:
    {
      SPHERE *a = one->gobj;
      CONVEX *b = two->gobj;

      return gjk_convex_sphere (b->cur, b->nver, a->cur_center, a->cur_radius, q, p);
    }
    break;
  }

  return 0;
}
