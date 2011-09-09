/*
 * gjk.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------------------
 * distance between a pair of convex polyhedra according to the classical algortihm
 * by Gilbert et al. IEEE J. of Robotics and Automation, 4/2, 1988, pp. 193-203
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
#include "gjk.h"
#include "alg.h"
#include "err.h"

/* auxiliary point structure used to describe
 * points in the Minkowski difference set (A - B) */
typedef struct { double w[3], *a, *b; } point;

/* simple bit lighting macro for vertex sets */
#define set(i,j,k,l) ((i*8)|(j*4)|(k*2)|(l*1))

/* number of bits lightened in numbers from 0 to 15;
 * they correspond to the cardinality of a vertex set
 * describing the inner simplex approximation of (A - B) */
static const int card [] =
{0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};

/* vertex schemes; the first three sets use only indices from {0, 1},
 * the first seven schemes use only indices {0, 1, 2}, and the complete
 * fifteen of them employs the complete set of indices {0, 1, 2, 3} */
static const int base [] =
{set(0,0,0,1),
 set(0,0,1,0),
 set(0,0,1,1),
 set(0,1,0,0),
 set(0,1,0,1),
 set(0,1,1,0),
 set(0,1,1,1),
 set(1,0,0,0),
 set(1,0,0,1),
 set(1,0,1,0),
 set(1,1,0,0),
 set(1,0,1,1),
 set(1,1,0,1),
 set(1,1,1,0),
 set(1,1,1,1)};

/* last [n] is the number of last base sets
 * to be checked during the projection on
 * a simplex with n vertices */
static const int last [] = 
{0, 0, 3, 7, 15};

/* find minimal point in set (c, n) along the direction of 'v' */
inline static double* minimal_support_point (double *c, int n, double *v)
{
  double dot, dotmin = DBL_MAX, *out = NULL;

  for (; n > 0; n --, c += 3)
  {
    dot = DOT (c, v);
    if (dot < dotmin) {dotmin = dot; out = c;}
  }

  return out;
}

/* find maximal point in set (c, n) along the direction of 'v' */
inline static double* maximal_support_point (double *c, int n, double *v)
{
  double dot, dotmax = -DBL_MAX, *out = NULL;

  for (; n > 0; n --, c += 3)
  {
    dot = DOT (c, v);
    if (dot > dotmax) {dotmax = dot; out = c;}
  }

  return out;
}

/* allocate output point for curved primitives */
inline static double* output_point (point *w, int n, double x [4][3], short maximal)
{
  double *out = NULL;
  int i, j;

  /* identify a spare return
   * pointer inside of 'x' */
  for (j = 0; j < 4; j ++)
  {
    for (i = 0; i < n; i ++)
    {
      if (maximal && w [i].b == x [j]) break; /* nope => this one is being used */
      else if (!maximal && w [i].a == x [j]) break;
    }

    if (i == n) /* found */
    {
      out = x [j];
      break;
    }
  }

  ASSERT_DEBUG (out, "Failed to identify a spare return pointer in maximal_sphere_support_point of gjk");

  return out;
}

/* find minimal point in the sphere (c, r) along the direction of 'v' */
inline static double* minimal_sphere_support_point (point *w, int n, double (*x) [3], double *c, double r, double *v)
{
  double *out = output_point (w, n, x, 0);

  COPY (v, out);
  NORMALIZE (out);
  SCALE (out, -r);
  ADD (c, out, out);

  return out;
}

/* find maximal point in the sphere (c, r) along the direction of 'v' */
inline static double* maximal_sphere_support_point (point *w, int n, double (*x) [3], double *c, double r, double *v)
{
  double *out = output_point (w, n, x, 1);

  COPY (v, out);
  NORMALIZE (out);
  SCALE (out, r);
  ADD (c, out, out);

  return out;
}

/* find minimal point in the ellipsoid (c, sca, rot) along the direction of 'v' */
inline static double* minimal_ellip_support_point (point *w, int n, double (*x) [3], double *c, double *sca, double *rot, double *v)
{
  /* (paraphrase of a comment from OpenTissue engine)
   * An ellipsoid E is a scaled, rotated and translated unit, zero-centered ball B.
   * Scaling and rotation can be respresented as a linear transformation T.
   * We can write it in a general way as
   *
   * E = T(B) + c, where T = {rot} * diag {sca}
   *
   * That means we can create a support function S of an ellipsoid from that of a unit sphere ball by
   *
   * S_E(v) = S_{T(B) + c}(v)
   *
   * Further we know that for any affine transformation, A(v) = T v + c, we have
   *
   * S_{A(B)}(v)  =   A(S_B(T'v))
   */

  double *out = output_point (w, n, x, 0);
  double T [9], q [3];

  NNCOPY (rot, T);
  SCALE (T, sca[0]);
  SCALE (T+3, sca[1]);
  SCALE (T+6, sca[2]);

  TVMUL (T, v, q);
  NORMALIZE (q);
  SCALE (q, -1.0); /* minimal support point on the unit sphere */
  NVMUL (T, q, out);
  ADD (c, out, out);

  return out;
}

/* find maximal point in the ellipsoid (c, sca, rot) along the direction of 'v' */
inline static double* maximal_ellip_support_point (point *w, int n, double (*x) [3], double *c, double *sca, double *rot, double *v)
{
  /* (paraphrase of a comment from OpenTissue engine)
   * An ellipsoid E is a scaled, rotated and translated unit, zero-centered ball B.
   * Scaling and rotation can be respresented as a linear transformation T.
   * We can write it in a general way as
   *
   * E = T(B) + c, where T = {rot} * diag {sca}
   *
   * That means we can create a support function S of an ellipsoid from that of a unit sphere ball by
   *
   * S_E(v) = S_{T(B) + c}(v)
   *
   * Further we know that for any affine transformation, A(v) = T v + c, we have
   *
   * S_{A(B)}(v)  =   A(S_B(T'v))
   */

  double *out = output_point (w, n, x, 1);
  double T [9], q [3];

  NNCOPY (rot, T);
  SCALE (T, sca[0]);
  SCALE (T+3, sca[1]);
  SCALE (T+6, sca[2]);

  TVMUL (T, v, q);
  NORMALIZE (q);
  NVMUL (T, q, out);
  ADD (c, out, out);

  return out;
}

/* projection of the zero point (0, 0, 0) onto the convex
 * hull spanned by points (w, n); on exit set 'w' is overwritten
 * by the smallest simplex in (w, n) such that projection of (0, 0, 0) 
 * onto it is the same as the projection onto the original hull; the
 * new dimension of 'w' is returend; 'l' contains barycentric coordinates
 * of the projection point; 'v' contains the point itslef */
static int project (point *w, int n, double *l, double *v)
{
  double dot [4][4], delta [16][4], sum;
  int i, j, k, m, s, c, o, f;

  if (n == 1)
  {
    l[0] = 1.0;
    COPY (w[0].w, v);
    return 1;
  }

  /* calculate all necessary dot products;
   * this could be optimised by calculating
   * dots related only to the new point */
  for (i = 0; i < n; i ++)
  {
    for (j = i; j < n; j ++)
    {
      dot [i][j] = DOT (w[i].w,w[j].w);
      dot [j][i] = dot [i][j];
    }
  }

  /* loop and calculate components
   * of the determinant expansion
   * according to the recursive formula
   * from the reference paper */
  for (m = 0; m < last [n]; m ++)
  {
    s = base [m];
    c = ~s;
    f = 1;
    if (card [s] == 1)
    {
      for (i = 0; i < n; i ++)
      {
	if (s & (1<<i))
	{
	  delta [s][i] = 1.0;
	  k = i;
	  break;
	}
      }
      for (j = 0; j < n; j ++)
      {
	if (c & (1<<j))
	{
	  delta [s][j] = dot[i][k] - dot[i][j];
	  if (delta [s][j] > 0.0) f = 0; /* cancel the flag => this is not yet the right sub-simplex */
	}
      }
    }
    else
    {
      for (i = 0, k = -1; i < n; i ++)
      {
	if (s & (1<<i))
	{
	  o = s;
	  o &= ~(1<<i);
	  if (k < 0) k = i;
	  delta [s][i] =  delta[o][i];
	  if (delta [s][i] <= 0.0) f = 0; 
	}
      }
      for (j = 0; j < n; j ++)
      {
	if (c & (1<<j))
	{
	  delta [s][j] = 0.0;
	  for (i = 0; i < n; i ++)
	  {
	    if (s & (1<<i))
	    {
	      delta [s][j] += delta [s][i] * (dot [i][k] - dot [i][j]); 
	    }
	  }
	  if (delta [s][j] > 0.0) f = 0; 
	}
      }
    }

    /* the flag was not canceled, that is
     * the projection on the current feature
     * was succesful => refresh vertex table
     * and compute barycentric coordinates
     * of the projection point */
    if (f)
    {
      sum = 0.0;
      for (i = o = 0; i < n; i ++)
      {
	if (s & (1<<i))
	{
	  w[o] = w[i];
	  l[o] = delta [s][i];
          sum += l[o];
	  o ++;
	}
      }
      SET (v, 0.0);
      for(i = 0; i < o; i ++)
      {
	l[i] /= sum; /* current coordinate */
	ADDMUL (v, l[i], w[i].w, v); /* convex combination == projection point */
      }
      return o;
    }
  }

  return 0;
}

/* public driver routine => input two polytopes A = (a, na) and B = (b, nb); outputs
 * p in A and q in B such that d = |p - q| is minimal; the distance d is returned */
double gjk (double *a, int na, double *b, int nb, double *p, double *q)
{
  point w [4];
  double v [3],
	 vlen,
	 delta,
	 mi = 0.0,
	 l [4];

  int toofar = 1,
      n = 0,
      j = 0,
      k = (na+nb)*(na+nb);

  SUB (a, b, v);
  vlen = LEN (v);

  while (toofar && vlen > GEOMETRIC_EPSILON && n < 4 && j ++ < k) /* (#) see below */
  {
    w[n].a = minimal_support_point (a, na, v);
    w[n].b = maximal_support_point (b, nb, v);
    SUB (w[n].a, w[n].b, w[n].w);
    delta = DOT (v, w[n].w) / vlen;
    mi = MAX (mi, delta);
    toofar = (vlen - mi) > GEOMETRIC_EPSILON;
    if (toofar)
    {
      n = project (w, ++n, l, v); /* (#) n = 4 can happen due to roundoff */
      vlen = LEN (v);
    }
  }

  if (n)
  {
    SET (p, 0);
    SET (q, 0);
    for (n --; n >= 0; n --)
    {
      ADDMUL (p, l[n], w[n].a, p);
      ADDMUL (q, l[n], w[n].b, q);
    }
  }
  else /* the while loop was never entered */
  {
    COPY (a, p);
    COPY (b, q);
  }

  return vlen;
}

/* public driver routine => input polytope A = (a, na) and sphere B = (c, r); outputs
 * p in A and q in B such that d = |p - q| is minimal; the distance d is returned */
double gjk_convex_sphere (double *a, int na, double *c, double r, double *p, double *q)
{
  point w [4];
  double v [3],
	 b [4][3], /* support points for the sphere */
	 vlen,
	 delta,
	 mi = 0.0,
	 l [4];

  int toofar = 1,
      n = 0,
      j = 0,
      k = 4*na*na;

  COPY (c, b [0]);
  b [0][rand () % 3] += r; /* be is now a point on the sphere */
  SUB (a, b [0], v); /* an initial point in the set A-B */
  vlen = LEN (v);

  while (toofar && vlen > GEOMETRIC_EPSILON && n < 4 && j ++ < k) /* (#) see below */
  {
    w[n].a = minimal_support_point (a, na, v);
    w[n].b = maximal_sphere_support_point (w, n, b, c, r, v);
    SUB (w[n].a, w[n].b, w[n].w);
    delta = DOT (v, w[n].w) / vlen;
    mi = MAX (mi, delta);
    toofar = (vlen - mi) > GEOMETRIC_EPSILON;
    if (toofar)
    {
      n = project (w, ++n, l, v); /* (#) n = 4 can happen due to roundoff */
      vlen = LEN (v);
    }
  }

  if (n)
  {
    SET (p, 0);
    SET (q, 0);
    for (n --; n >= 0; n --)
    {
      ADDMUL (p, l[n], w[n].a, p);
      ADDMUL (q, l[n], w[n].b, q);
    }
  }
  else /* while loop never entered */
  {
    COPY (a, p);
    COPY (b [0], q); 
  }

  return vlen;
}

/* (a,na) and (b,bsca, brot) are the input polyhedron and ellipsoid; 'p' and 'q' are the two outputed
 * closest points, respectively in polyhedron (a,na) and ellipsoid (b, bsca, brot); the distance is returned */
double gjk_convex_ellip (double *a, int na, double *b, double *bsca, double *brot, double *p, double *q)
{
  point w [4];
  double v [3],
	 z [4][3], /* support points for the ellipsoid */
	 vlen,
	 delta,
	 mi = 0.0,
	 l [4];

  int toofar = 1,
      n = 0,
      j = 0,
      k = 4*na*na,
      i;

  SET (v, 0);
  i = rand () % 3;
  v[i] = bsca [i]; /* a point on a scaled unit sphere */
  NVMUL (brot, v, z [0]); /* a point on a scaled, rotated unit sphere */
  ADD (b, z [0], z [0]); /* a point on a scaled, rotated and translated unit sphere => ellipsoid */

  SUB (a, z [0], v); /* an initial point in the set A-B */
  vlen = LEN (v);

  while (toofar && vlen > GEOMETRIC_EPSILON && n < 4 && j ++ < k) /* (#) see below */
  {
    w[n].a = minimal_support_point (a, na, v);
    w[n].b = maximal_ellip_support_point (w, n, z, b, bsca, brot, v);
    SUB (w[n].a, w[n].b, w[n].w);
    delta = DOT (v, w[n].w) / vlen;
    mi = MAX (mi, delta);
    toofar = (vlen - mi) > GEOMETRIC_EPSILON;
    if (toofar)
    {
      n = project (w, ++n, l, v); /* (#) n = 4 can happen due to roundoff */
      vlen = LEN (v);
    }
  }

  if (n)
  {
    SET (p, 0);
    SET (q, 0);
    for (n --; n >= 0; n --)
    {
      ADDMUL (p, l[n], w[n].a, p);
      ADDMUL (q, l[n], w[n].b, q);
    }
  }
  else /* while loop never entered */
  {
    COPY (a, p);
    COPY (z [0], q); 
  }

  return vlen;
}

/* (a,na) and p are the input polyhedron and point; 'q' is the outputed
 * closest point on the polyhedron; the distance is returned */
double gjk_convex_point (double *a, int na, double *p, double *q)
{
  point w [4];
  double v [3],
	 vlen,
	 delta,
	 mi = 0.0,
	 l [4];

  int toofar = 1,
      n = 0,
      j = 0,
      k = 4*na*na;

  SUB (a, p, v); /* an initial point in the set A-B */
  vlen = LEN (v);

  while (toofar && vlen > GEOMETRIC_EPSILON && n < 4 && j ++ < k) /* (#) see below */
  {
    w[n].a = minimal_support_point (a, na, v);
    w[n].b = p;
    SUB (w[n].a, w[n].b, w[n].w);
    delta = DOT (v, w[n].w) / vlen;
    mi = MAX (mi, delta);
    toofar = (vlen - mi) > GEOMETRIC_EPSILON;
    if (toofar)
    {
      n = project (w, ++n, l, v); /* (#) n = 4 can happen due to roundoff */
      vlen = LEN (v);
    }
  }

  if (n)
  {
    SET (q, 0);
    for (n --; n >= 0; n --)
    {
      ADDMUL (q, l[n], w[n].a, q);
    }
  }
  else /* while loop never entered */
  {
    COPY (a, q);
  }

  return vlen;
}

/* public driver routine => input sphere A = (a, ra) and sphere B = (b, rb); outputs
 * p in A and q in B such that d = |p - q| is minimal; the distance d is returned */
double gjk_sphere_sphere (double *a, double ra, double *b, double rb, double *p, double *q)
{
  double d [3], dlen, deno;

  SUB (a, b, d);
  dlen = LEN (d);
  if (dlen == 0.0)
  { COPY (a, p);
    COPY (b, q);
    return 0.0; }
  deno = ra / dlen;
  SUBMUL (a, deno, d, p);
  deno = rb / dlen;
  ADDMUL (b, deno, d, q);
  deno = ra + rb;

  if (dlen > deno)
    return (dlen - deno);
  else
  {
    MID (p, q, p);
    COPY (p, q);
    return 0.0;
  }
}

/* (a,ra) and (b,bsca,brot)) are the input sphere and ellipsoid; 'p' and 'q' are the two outputed
 * closest points, respectively in sphere (a,ra) and ellipsoid (b,bsca,brot); the distance is returned */
double gjk_sphere_ellip (double *a, double ra, double *b, double *bsca, double *brot, double *p, double *q)
{
  point w [4];
  double v [3],
	 y [4][3], /* support points for the sphere */
	 z [4][3], /* support points for the ellipsoid */
	 vlen,
	 delta,
	 mi = 0.0,
	 l [4];

  int toofar = 1,
      n = 0,
      j = 0,
      k = 128,
      i;

  i = rand () % 3;
  COPY (a, y [0]);
  y [0][i] += ra; /* y[0] is now a point on the sphere */

  SET (v, 0);
  v[i] = bsca [i]; /* a point on a scaled unit sphere */
  NVMUL (brot, v, z [0]); /* a point on a scaled, rotated unit sphere */
  ADD (b, z [0], z [0]); /* a point on a scaled, rotated and translated unit sphere => ellipsoid */

  SUB (y [0], z [0], v); /* an initial point in the set A-B */
  vlen = LEN (v);

  while (toofar && vlen > GEOMETRIC_EPSILON && n < 4 && j ++ < k) /* (#) see below */
  {
    w[n].a = minimal_sphere_support_point (w, n, y, a, ra, v);
    w[n].b = maximal_ellip_support_point (w, n, z, b, bsca, brot, v);
    SUB (w[n].a, w[n].b, w[n].w);
    delta = DOT (v, w[n].w) / vlen;
    mi = MAX (mi, delta);
    toofar = (vlen - mi) > GEOMETRIC_EPSILON;
    if (toofar)
    {
      n = project (w, ++n, l, v); /* (#) n = 4 can happen due to roundoff */
      vlen = LEN (v);
    }
  }

  if (n)
  {
    SET (p, 0);
    SET (q, 0);
    for (n --; n >= 0; n --)
    {
      ADDMUL (p, l[n], w[n].a, p);
      ADDMUL (q, l[n], w[n].b, q);
    }
  }
  else /* while loop never entered */
  {
    COPY (y [0], p);
    COPY (z [0], q); 
  }

  return vlen;
}

/* (a,asca,arot) and (b,bsca,brot) are the two input ellipsoids; 'p' and 'q' are the two outputed
 * closest points, respectively in (a,asca,arot) and (b,bsca,brot); the distance is returned */
double gjk_ellip_ellip (double *a, double *asca, double *arot, double *b, double *bsca, double *brot, double *p, double *q)
{
  point w [4];
  double v [3],
	 y [4][3], /* support points for the first ellipsoid */
	 z [4][3], /* support points for the second ellipsoid */
	 vlen,
	 delta,
	 mi = 0.0,
	 l [4];

  int toofar = 1,
      n = 0,
      j = 0,
      k = 128,
      i;

  i = rand () % 3;

  SET (v, 0);
  v[i] = asca [i]; /* a point on a scaled unit sphere */
  NVMUL (arot, v, y [0]); /* a point on a scaled, rotated unit sphere */
  ADD (a, y [0], y [0]); /* a point on a scaled, rotated and translated unit sphere => ellipsoid */

  SET (v, 0);
  v[i] = bsca [i]; /* a point on a scaled unit sphere */
  NVMUL (brot, v, z [0]); /* a point on a scaled, rotated unit sphere */
  ADD (b, z [0], z [0]); /* a point on a scaled, rotated and translated unit sphere => ellipsoid */

  SUB (y [0], z [0], v); /* an initial point in the set A-B */
  vlen = LEN (v);

  while (toofar && vlen > GEOMETRIC_EPSILON && n < 4 && j ++ < k) /* (#) see below */
  {
    w[n].a = minimal_ellip_support_point (w, n, y, a, asca, arot, v);
    w[n].b = maximal_ellip_support_point (w, n, z, b, bsca, brot, v);
    SUB (w[n].a, w[n].b, w[n].w);
    delta = DOT (v, w[n].w) / vlen;
    mi = MAX (mi, delta);
    toofar = (vlen - mi) > GEOMETRIC_EPSILON;
    if (toofar)
    {
      n = project (w, ++n, l, v); /* (#) n = 4 can happen due to roundoff */
      vlen = LEN (v);
    }
  }

  if (n)
  {
    SET (p, 0);
    SET (q, 0);
    for (n --; n >= 0; n --)
    {
      ADDMUL (p, l[n], w[n].a, p);
      ADDMUL (q, l[n], w[n].b, q);
    }
  }
  else /* while loop never entered */
  {
    COPY (y [0], p);
    COPY (z [0], q); 
  }

  return vlen;
}

/* (a,asca,arot) and p are the input ellipsoid and point; 'q' is the outputed
 * closest point on the ellipsoid; the distance is returned */
double gjk_ellip_point (double *a, double *asca, double *arot, double *p, double *q)
{
  point w [4];
  double v [3],
	 y [4][3], /* support points for the ellipsoid */
	 vlen,
	 delta,
	 mi = 0.0,
	 l [4];

  int toofar = 1,
      n = 0,
      j = 0,
      k = 128,
      i;

  i = rand () % 3;

  SET (v, 0);
  v[i] = asca [i]; /* a point on a scaled unit sphere */
  NVMUL (arot, v, y [0]); /* a point on a scaled, rotated unit sphere */
  ADD (a, y [0], y [0]); /* a point on a scaled, rotated and translated unit sphere => ellipsoid */

  SUB (y [0], p, v); /* an initial point in the set A-B */
  vlen = LEN (v);

  while (toofar && vlen > GEOMETRIC_EPSILON && n < 4 && j ++ < k) /* (#) see below */
  {
    w[n].a = minimal_ellip_support_point (w, n, y, a, asca, arot, v);
    w[n].b = p;
    SUB (w[n].a, w[n].b, w[n].w);
    delta = DOT (v, w[n].w) / vlen;
    mi = MAX (mi, delta);
    toofar = (vlen - mi) > GEOMETRIC_EPSILON;
    if (toofar)
    {
      n = project (w, ++n, l, v); /* (#) n = 4 can happen due to roundoff */
      vlen = LEN (v);
    }
  }

  if (n)
  {
    SET (q, 0);
    for (n --; n >= 0; n --)
    {
      ADDMUL (q, l[n], w[n].a, q);
    }
  }
  else /* while loop never entered */
  {
    COPY (y [0], q);
  }

  return vlen;
}

/* compute gap function betwen two primitives along the given unit normal;
 * the normal direction is assumed to be outward to the first primitive */
double gjk_convex_convex_gap (double *a, int na, double *b, int nb, double *normal)
{
  double *p, *q, d [3];

  p = maximal_support_point (a, na, normal);
  q = minimal_support_point (b, nb, normal);
  SUB (q, p, d);

  return DOT (normal, d);
}

double gjk_convex_sphere_gap (double *a, int na, double *b, double rb, double *normal)
{
  double *p, *q, d [3], y [1][3];

  p = maximal_support_point (a, na, normal);
  q = minimal_sphere_support_point (NULL, 0, y, b, rb, normal);
  SUB (q, p, d);

  return DOT (normal, d);
}

double gjk_convex_ellip_gap (double *a, int na, double *b, double *bsca, double *brot, double *normal)
{
  double *p, *q, d [3], y [1][3];

  p = maximal_support_point (a, na, normal);
  q = minimal_ellip_support_point (NULL, 0, y, b, bsca, brot, normal);
  SUB (q, p, d);

  return DOT (normal, d);
}

double gjk_sphere_sphere_gap (double *a, double ra, double *b, double rb, double *normal)
{
  double *p, *q, d [3], y [1][3], z [1][3];

  p = maximal_sphere_support_point (NULL, 0, y, a, ra, normal);
  q = minimal_sphere_support_point (NULL, 0, z, b, rb, normal);
  SUB (q, p, d);

  return DOT (normal, d);
}

double gjk_sphere_ellip_gap (double *a, double ra, double *b, double *bsca, double *brot, double *normal)
{
  double *p, *q, d [3], y [1][3], z [1][3];

  p = maximal_sphere_support_point (NULL, 0, y, a, ra, normal);
  q = minimal_ellip_support_point (NULL, 0, z, b, bsca, brot, normal);
  SUB (q, p, d);

  return DOT (normal, d);
}

double gjk_ellip_ellip_gap (double *a, double *asca, double *arot, double *b, double *bsca, double *brot, double *normal)
{
  double *p, *q, d [3], y [1][3], z [1][3];

  p = maximal_ellip_support_point (NULL, 0, y, a, asca, arot, normal);
  q = minimal_ellip_support_point (NULL, 0, z, b, bsca, brot, normal);
  SUB (q, p, d);

  return DOT (normal, d);
}

/* compute furthest or closest (near == 0 or 1) point 'p' of a primitive along given normal direction */
void gjk_ellip_support_point (double *a, double *sca, double *rot, double *normal, short near, double *p)
{
  if (near) minimal_ellip_support_point (NULL, 0, (double (*) [3])p, a, sca, rot, normal);
  else maximal_ellip_support_point (NULL, 0, (double (*) [3])p, a, sca, rot, normal);
}
