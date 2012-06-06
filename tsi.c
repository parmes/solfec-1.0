/*
 * tsi.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * triangle-sphere intersection
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
#include <math.h>
#include "tsi.h"
#include "alg.h"
#include "err.h"

#define UV(pnt, u, v)\
  SUB (pnt, a, l);\
  uv = - DOT (ba, ac);\
  lv = - DOT (l, ac);\
  lu =  DOT (l, ba);\
  tmp = uv*uv - ba2*ac2;\
  u = (uv*lv - ac2*lu) / tmp;\
  v = (uv*lu - ba2*lv) / tmp;\
 
#define UVW(pnt, u, v, w)\
  UV (pnt, u, v);\
  w = u + v

#define SPHCENTER(pnt)\
  PRODUCT (ba, ac, n);\
  ADD (p, n, s);\
  SUB (a, s, d);\
  SUB (p, s, l);\
  tmp = DOT (n, d) / DOT (n, l);\
  ADDMUL (s, tmp, l, pnt)

/* get status of intersection between triangle (a, b, c) and sphere (p, r) */
int TSI_Status (double *a, double *b, double *c, double *p, double r)
{
  /* macro internals */
  double n [3], s [3], d [3], l [3];
  double uv, lv, lu, tmp;

  /* local variables */
  double u, v, w, ba2, ac2, rsq, rsqeps;
  double pa [3], pb [3], pc [3];
  double ba [3], cb [3], ac [3];
  double x [3], px [3];
  unsigned char abc = 0x00;
  
  SUB (p, a, pa);
  SUB (p, b, pb);
  SUB (p, c, pc);

  rsq = r * r;
  rsqeps = rsq - GEOMETRIC_EPSILON * rsq;

  /* vertices inside */
  if (DOT (pa, pa) < rsqeps) abc |= 0x01;
  if (DOT (pb, pb) < rsqeps) abc |= 0x02;
  if (DOT (pc, pc) < rsqeps) abc |= 0x04;

  SUB (b, a, ba);
  SUB (c, b, cb);
  SUB (a, c, ac);

  ba2 = DOT (ba, ba);
  ac2 = DOT (ac, ac);

  u = DOT (pa, ba) / ba2;
  v = DOT (pb, cb) / DOT (cb, cb);
  w = DOT (pc, ac) / ac2;

  /* edges inside */
  ADDMUL (a, u, ba, x); SUB (p, x, px);
  if (0. <= u && u <= 1. && DOT (px, px) < rsqeps) abc |= 0x10;
  ADDMUL (b, v, cb, x); SUB (p, x, px);
  if (0. <= v && v <= 1. && DOT (px, px) < rsqeps) abc |= 0x20;
  ADDMUL (c, w, ac, x); SUB (p, x, px);
  if (0. <= w && w <= 1. && DOT (px, px) < rsqeps) abc |= 0x40;

  switch (abc)
  {
  case 0x0:
  {
    SPHCENTER (x);
    SUB (p, x, px);
    UVW (x, u, v, w);

    /* area inside */
    if (u > 0. && v > 0. && w < 1. && DOT (px, px) < rsqeps) return TSI_BUBBLE;
  }
  break;
  case 0x21:
  case 0x31:
  case 0x61:
  case 0x71:
	    return TSI_A_BC;
  case 0x42:
  case 0x62:
  case 0x52:
  case 0x72:
	     return TSI_B_CA;
  case 0x14:
  case 0x34:
  case 0x54:
  case 0x74:
	     return TSI_C_AB;

  case 0x10: return TSI_AB;
  case 0x20: return TSI_BC;
  case 0x40: return TSI_CA;
  case 0x30: return TSI_AB_BC;
  case 0x60: return TSI_BC_CA;
  case 0x50: return TSI_CA_AB;
  case 0x70: return TSI_AB_BC_CA;
  }

  switch (abc & 0x0f)
  {
  case 0x01: return TSI_A;
  case 0x02: return TSI_B;
  case 0x04: return TSI_C;
  case 0x03: return TSI_A_B;
  case 0x06: return TSI_B_C;
  case 0x05: return TSI_C_A;
  case 0x07: return TSI_A_B_C;
  }

  return TSI_OUT;
}

#define CIRCLE(I)\
  SPHCENTER (pp);\
  SUB (I, pp, l);\
  rp = LEN (l)

#define ONCIRCLE(A, I)\
  SUB (A, pp, l);\
  NORMALIZE (l);\
  ADDMUL (pp, rp, l, I)

#define SPHERE_SEGMENT_1(A, B, pnt)\
  SUB (B, A, n);\
  SUB (A, p, s);\
  d [0] = DOT (n, n);\
  d [1] = 2. * DOT (n, s);\
  d [2] = DOT (p, p) + DOT (A, A) -\
    2. * DOT (p, A) - rsq;\
  tmp = sqrt (d [1]*d [1] - 4.*d [0] *d [2]);\
  uv = (- d [1] + tmp) / (2. * d [0]);\
  ADDMUL (A, uv, n, pnt)

#define SPHERE_SEGMENT_2(A, B, pnt1, pnt2)\
  SPHERE_SEGMENT_1(A, B, pnt2);\
  lv = (- d [1] - tmp) / (2. * d [0]);\
  ADDMUL (A, lv, n, pnt1)

#define ALLOC(NV, NT)\
  ERRMEM ((*cc) = malloc (sizeof (double [2]) * (NV)));\
  (*ncc) = NV;\
  ERRMEM ((*tt) = malloc (sizeof (int [6]) * (NT)));\
  (*ntt) = NT

#define TRI(N, V0, V1, V2, V3, V4, V5)\
  (*tt)[N*6] = V0;\
  (*tt)[N*6 + 1] = V1;\
  (*tt)[N*6 + 2] = V2;\
  (*tt)[N*6 + 3] = V3;\
  (*tt)[N*6 + 4] = V4;\
  (*tt)[N*6 + 5] = V5

#define TSI_V(A, B, C)\
  SPHERE_SEGMENT_1(A, B, i1);\
  SPHERE_SEGMENT_1(A, C, i2);\
  MID (i1, i2, i3);\
  CIRCLE (i1);\
  ONCIRCLE (i3, i3);\
  ALLOC (6, 1);\
  UV (A, (*cc)[0], (*cc)[1]);\
  UV (i1, (*cc)[2], (*cc)[3]);\
  UV (i2, (*cc)[4], (*cc)[5]);\
  MID (A, i1, i1);\
  MID (A, i2, i2);\
  UV (i1, (*cc)[6], (*cc)[7]);\
  UV (i3, (*cc)[8], (*cc)[9]);\
  UV (i2, (*cc)[10], (*cc)[11]);\
  TRI (0, 0, 1, 2, 3, 4, 5)

#define TSI_V_V(A, B, C)\
  SPHERE_SEGMENT_1(A, C, i1);\
  SPHERE_SEGMENT_1(B, C, i2);\
  MID (i1, i2, i3);\
  CIRCLE (i1);\
  ONCIRCLE (i3, i3);\
  ALLOC (9, 2);\
  UV (A, (*cc)[0], (*cc)[1]);\
  UV (B, (*cc)[2], (*cc)[3]);\
  UV (i2, (*cc)[4], (*cc)[5]);\
  UV (i1, (*cc)[6], (*cc)[7]);\
  MID (A, B, n);\
  UV (n, (*cc)[8], (*cc)[9]);\
  MID (B, i2, n);\
  UV (n, (*cc)[10], (*cc)[11]);\
  UV (i3, (*cc)[12], (*cc)[13]);\
  MID (A, i1, n);\
  UV (n, (*cc)[14], (*cc)[15]);\
  MID (i1, B, n);\
  UV (n, (*cc)[16], (*cc)[17]);\
  TRI (0, 0, 1, 3, 4, 8, 7);\
  TRI (1, 1, 2, 3, 5, 6, 8)

#define TSI_V_V_V(A, B, C)\
  ALLOC (6, 1);\
  UV (A, (*cc)[0], (*cc)[1]);\
  UV (B, (*cc)[2], (*cc)[3]);\
  UV (C, (*cc)[4], (*cc)[5]);\
  MID (A, B, n);\
  UV (n, (*cc)[6], (*cc)[7]);\
  MID (B, C, n);\
  UV (n, (*cc)[8], (*cc)[9]);\
  MID (C, A, n);\
  UV (n, (*cc)[10], (*cc)[11]);\
  TRI (0, 0, 1, 2, 3, 4, 5)

#define TSI_V_VV(A, B, C)\
  SPHERE_SEGMENT_1(A, B, i1);\
  SPHERE_SEGMENT_2(B, C, i2, i3);\
  SPHERE_SEGMENT_1(A, C, i4);\
  CIRCLE (i1);\
  ALLOC (12, 3);\
  UV (A, (*cc)[0], (*cc)[1]);\
  UV (i1, (*cc)[2], (*cc)[3]);\
  UV (i2, (*cc)[4], (*cc)[5]);\
  UV (i3, (*cc)[6], (*cc)[7]);\
  UV (i4, (*cc)[8], (*cc)[9]);\
  MID (A, i1, n);\
  UV (n, (*cc)[10], (*cc)[11]);\
  MID (i1, i2, n);\
  ONCIRCLE (n, n);\
  UV (n, (*cc)[12], (*cc)[13]);\
  MID (i2, i3, n);\
  UV (n, (*cc)[14], (*cc)[15]);\
  MID (i3, i4, n);\
  ONCIRCLE (n, n);\
  UV (n, (*cc)[16], (*cc)[17]);\
  MID (i4, A, n);\
  UV (n, (*cc)[18], (*cc)[19]);\
  MID (A, i2, n);\
  UV (n, (*cc)[20], (*cc)[21]);\
  MID (i3, A, n);\
  UV (n, (*cc)[22], (*cc)[23]);\
  TRI (0, 0, 1, 2, 5, 6, 10);\
  TRI (1, 0, 2, 3, 10, 7, 11);\
  TRI (2, 0, 3, 4, 11, 8, 9)

#define TSI_VV(A, B, C)\
  SPHERE_SEGMENT_2(A, B, i1, i2);\
  MID (i1, i2, i4);\
  CIRCLE (i1);\
  UVW (pp, u, v, w);\
  if (u > 0. && v > 0. && w < 1.)\
  {\
    ONCIRCLE (C, i3);\
    ALLOC (10, 3);\
    UV (i1, (*cc)[0], (*cc)[1]);\
    UV (i2, (*cc)[2], (*cc)[3]);\
    UV (i3, (*cc)[4], (*cc)[5]);\
    UV (i4, (*cc)[6], (*cc)[7]);\
    MID (i2, i3, n);\
    ONCIRCLE (n, n);\
    UV (n, (*cc)[8], (*cc)[9]);\
    MID (i3, i1, n);\
    ONCIRCLE (n, n);\
    UV (n, (*cc)[10], (*cc)[11]);\
    MID (i1, pp, n);\
    UV (n, (*cc)[12], (*cc)[13]);\
    MID (i2, pp, n);\
    UV (n, (*cc)[14], (*cc)[15]);\
    MID (i3, pp, n);\
    UV (n, (*cc)[16], (*cc)[17]);\
    UV (pp, (*cc)[18], (*cc)[19]);\
    TRI (0, 0, 9, 2, 6, 8, 5);\
    TRI (1, 0, 1, 9, 3, 7, 6);\
    TRI (2, 1, 2, 9, 4, 8, 7);\
  }\
  else\
  {\
    SUB (A, B, i5);\
    NORMAL (A, B, C, n);\
    PRODUCT (i5, n, i6);\
    SUB (C, i4, i5);\
    if (DOT (i5, i6) < 0.)\
    { SCALE (i6, -1.); }\
    ADD (i4, i6, i6);\
    ONCIRCLE (i6, i3);\
    ALLOC (6, 1);\
    UV (i1, (*cc)[0], (*cc)[1]);\
    UV (i2, (*cc)[2], (*cc)[3]);\
    UV (i3, (*cc)[4], (*cc)[5]);\
    MID (i2, i3, i2);\
    ONCIRCLE (i2, i2);\
    MID (i3, i1, i1);\
    ONCIRCLE (i1, i1);\
    UV (i4, (*cc)[6], (*cc)[7]);\
    UV (i2, (*cc)[8], (*cc)[9]);\
    UV (i1, (*cc)[10], (*cc)[11]);\
    TRI (0, 0, 1, 2, 3, 4, 5);\
  }
  
#define TSI_VV_VV(A, B, C)\
  SPHERE_SEGMENT_2(A, B, i1, i2);\
  SPHERE_SEGMENT_2(B, C, i3, i4);\
  CIRCLE (i1);\
  UVW (pp, u, v, w);\
  if (u > 0. && v > 0. && w < 1.)\
  {\
    ALLOC (16, 5);\
    MID (A, C, d);\
    ONCIRCLE (d, d);\
    UV (i1, (*cc)[0], (*cc)[1]);\
    UV (i2, (*cc)[2], (*cc)[3]);\
    UV (i3, (*cc)[4], (*cc)[5]);\
    UV (i4, (*cc)[6], (*cc)[7]);\
    MID (i1, i2, n);\
    UV (n, (*cc)[8], (*cc)[9]);\
    MID (i2, i3, n);\
    ONCIRCLE (n, n);\
    UV (n, (*cc)[10], (*cc)[11]);\
    MID (i3, i4, n);\
    UV (n, (*cc)[12], (*cc)[13]);\
    MID (i4, d, n);\
    ONCIRCLE (n, n);\
    UV (n, (*cc)[14], (*cc)[15]);\
    UV (d, (*cc)[16], (*cc)[17]);\
    MID (i1, d, n);\
    ONCIRCLE (n, n);\
    UV (n, (*cc)[18], (*cc)[19]);\
    MID (i1, pp, n);\
    UV (n, (*cc)[20], (*cc)[21]);\
    MID (i2, pp, n);\
    UV (n, (*cc)[22], (*cc)[23]);\
    MID (i3, pp, n);\
    UV (n, (*cc)[24], (*cc)[25]);\
    MID (i4, pp, n);\
    UV (n, (*cc)[26], (*cc)[27]);\
    MID (d, pp, n);\
    UV (n, (*cc)[28], (*cc)[29]);\
    UV (pp, (*cc)[30], (*cc)[31]);\
    TRI (0, 0, 1, 15, 4, 11, 10);\
    TRI (1, 1, 2, 15, 5, 12, 11);\
    TRI (2, 2, 3, 15, 6, 13, 12);\
    TRI (3, 3, 8, 15, 7, 14, 13);\
    TRI (4, 8, 0, 15, 9, 10, 14);\
  }\
  else\
  {\
    ALLOC (9, 2);\
    UV (i1, (*cc)[0], (*cc)[1]);\
    UV (i2, (*cc)[2], (*cc)[3]);\
    UV (i3, (*cc)[4], (*cc)[5]);\
    UV (i4, (*cc)[6], (*cc)[7]);\
    MID (i1, i2, n);\
    UV (n, (*cc)[8], (*cc)[9]);\
    MID (i2, i3, n);\
    ONCIRCLE (n, n);\
    UV (n, (*cc)[10], (*cc)[11]);\
    MID (i3, i4, n);\
    UV (n, (*cc)[12], (*cc)[13]);\
    MID (i1, i4, n);\
    ONCIRCLE (n, n);\
    UV (n, (*cc)[14], (*cc)[15]);\
    MID (i2, i4, n);\
    UV (n, (*cc)[16], (*cc)[17]);\
    TRI (0, 0, 1, 3, 4, 8, 7);\
    TRI (1, 1, 2, 3, 5, 6, 8);\
  }

#define TSI_VV_VV_VV(A, B, C)\
  SPHERE_SEGMENT_2(A, B, i1, i2);\
  SPHERE_SEGMENT_2(B, C, i3, i4);\
  SPHERE_SEGMENT_2(C, A, i5, i6);\
  CIRCLE (i1);\
  ALLOC (19, 6);\
  UV (i1, (*cc)[0], (*cc)[1]);\
  UV (i2, (*cc)[2], (*cc)[3]);\
  UV (i3, (*cc)[4], (*cc)[5]);\
  UV (i4, (*cc)[6], (*cc)[7]);\
  UV (i5, (*cc)[8], (*cc)[9]);\
  UV (i6, (*cc)[10], (*cc)[11]);\
  MID (i1, i2, n);\
  UV (n, (*cc)[12], (*cc)[13]);\
  MID (i2, i3, n);\
  ONCIRCLE (n, n);\
  UV (n, (*cc)[14], (*cc)[15]);\
  MID (i3, i4, n);\
  UV (n, (*cc)[16], (*cc)[17]);\
  MID (i4, i5, n);\
  ONCIRCLE (n, n);\
  UV (n, (*cc)[18], (*cc)[19]);\
  MID (i5, i6, n);\
  UV (n, (*cc)[20], (*cc)[21]);\
  MID (i6, i1, n);\
  ONCIRCLE (n, n);\
  UV (n, (*cc)[22], (*cc)[23]);\
  MID (i1, pp, n);\
  UV (n, (*cc)[24], (*cc)[25]);\
  MID (i2, pp, n);\
  UV (n, (*cc)[26], (*cc)[27]);\
  MID (i3, pp, n);\
  UV (n, (*cc)[28], (*cc)[29]);\
  MID (i4, pp, n);\
  UV (n, (*cc)[30], (*cc)[31]);\
  MID (i5, pp, n);\
  UV (n, (*cc)[32], (*cc)[33]);\
  MID (i6, pp, n);\
  UV (n, (*cc)[34], (*cc)[35]);\
  UV (pp, (*cc)[36], (*cc)[37]);\
  TRI (0, 0, 1, 18, 6, 13, 12);\
  TRI (1, 1, 2, 18, 7, 14, 13);\
  TRI (2, 2, 3, 18, 8, 15, 14);\
  TRI (3, 3, 4, 18, 9, 16, 15);\
  TRI (4, 4, 5, 18, 10, 17, 16);\
  TRI (5, 5, 0, 18, 11, 12, 17)

#define TSI_DOBUBBLE(A, B, C)\
  SPHCENTER (pp);\
  SPHERE_SEGMENT_1(pp, A, i1);\
  SPHERE_SEGMENT_1(pp, B, i2);\
  SPHERE_SEGMENT_1(pp, C, i3);\
  SUB (i1, pp, l);\
  rp = LEN (l);\
  MID (i1, i2, i4);\
  ONCIRCLE (i4, i4);\
  MID (i2, i3, i5);\
  ONCIRCLE (i5, i5);\
  MID (i3, i1, i6);\
  ONCIRCLE (i6, i6);\
  ALLOC (10, 3);\
  UV (i1, (*cc)[0], (*cc)[1]);\
  UV (i2, (*cc)[2], (*cc)[3]);\
  UV (i3, (*cc)[4], (*cc)[5]);\
  UV (i4, (*cc)[6], (*cc)[7]);\
  UV (i5, (*cc)[8], (*cc)[9]);\
  UV (i6, (*cc)[10], (*cc)[11]);\
  MID (i1, pp, n);\
  UV (n, (*cc)[12], (*cc)[13]);\
  MID (i2, pp, n);\
  UV (n, (*cc)[14], (*cc)[15]);\
  MID (i3, pp, n);\
  UV (n, (*cc)[16], (*cc)[17]);\
  UV (pp, (*cc)[18], (*cc)[19]);\
  TRI (0, 0, 1, 9, 3, 7, 6);\
  TRI (1, 1, 2, 9, 4, 8, 7);\
  TRI (2, 2, 0, 9, 5, 6, 8)

/* approximate intersection of triangle (a, b, c) and sphere (p, r) with second order triangular elements */
int TSI_Approx (double *a, double *b, double *c, double *p, double r, double **cc, int *ncc, int **tt, int *ntt)
{
  /* macro internals */
  double i1 [3], i2 [3], i3 [3], i4 [3];
  double uv, lv, lu, tmp, pp [3], rp;
  double n [3], s [3], d [3], l [3];
  double i5 [3], i6 [3];

  /* local variables */
  double u, v, w, ba2, ac2, rsq, rsqeps;
  double pa [3], pb [3], pc [3];
  double ba [3], cb [3], ac [3];
  double x [3], px [3];
  unsigned char abc = 0x00;

  /* nullify output */
  *cc = NULL; *ncc = 0;
  *tt = NULL; *ntt = 0;
  
  SUB (p, a, pa);
  SUB (p, b, pb);
  SUB (p, c, pc);

  rsq = r * r;
  rsqeps = rsq - GEOMETRIC_EPSILON * rsq;

  /* vertices inside */
  if (DOT (pa, pa) < rsqeps) abc |= 0x01;
  if (DOT (pb, pb) < rsqeps) abc |= 0x02;
  if (DOT (pc, pc) < rsqeps) abc |= 0x04;

  SUB (b, a, ba);
  SUB (c, b, cb);
  SUB (a, c, ac);

  ba2 = DOT (ba, ba);
  ac2 = DOT (ac, ac);

  /* project center on triangle edges */
  u = DOT (pa, ba) / ba2;
  v = DOT (pb, cb) / DOT (cb, cb);
  w = DOT (pc, ac) / ac2;

  /* edges inside */
  ADDMUL (a, u, ba, x); SUB (p, x, px);
  if (0. <= u && u <= 1. && DOT (px, px) < rsqeps) abc |= 0x10;
  ADDMUL (b, v, cb, x); SUB (p, x, px);
  if (0. <= v && v <= 1. && DOT (px, px) < rsqeps) abc |= 0x20;
  ADDMUL (c, w, ac, x); SUB (p, x, px);
  if (0. <= w && w <= 1. && DOT (px, px) < rsqeps) abc |= 0x40;

  switch (abc)
  {
  case 0x0:
  {
    SPHCENTER (x);
    SUB (p, x, px);
    UVW (x, u, v, w);

    /* area inside */
    if (u > 0. && v > 0. && w < 1. && DOT (px, px) < rsqeps)
    {
      TSI_DOBUBBLE (a, b, c);
      return TSI_BUBBLE;
    }
  }
  break;
  case 0x21:
  case 0x31:
  case 0x61:
  case 0x71:
    TSI_V_VV (a, b, c);
    return TSI_A_BC;
  case 0x42:
  case 0x62:
  case 0x52:
  case 0x72:
    TSI_V_VV (b, c, a);
    return TSI_B_CA;
  case 0x14:
  case 0x34:
  case 0x54:
  case 0x74:
    TSI_V_VV (c, a, b);
    return TSI_C_AB;
  case 0x10:
    TSI_VV (a, b, c);
    return TSI_AB;
  case 0x20:
    TSI_VV (b, c, a);
    return TSI_BC;
  case 0x40:
    TSI_VV (c, a, b);
    return TSI_CA;
  case 0x30:
    TSI_VV_VV (a, b, c);
    return TSI_AB_BC;
  case 0x60:
    TSI_VV_VV (b, c, a);
    return TSI_BC_CA;
  case 0x50:
    TSI_VV_VV (c, a, b);
    return TSI_CA_AB;
  case 0x70:
    TSI_VV_VV_VV (a, b, c);
    return TSI_AB_BC_CA;
  }

  switch (abc & 0x0f)
  {
  case 0x01:
    TSI_V (a, b, c);
    return TSI_A;
  case 0x02:
    TSI_V (b, c, a);
    return TSI_B;
  case 0x04:
    TSI_V (c, a, b);
    return TSI_C;
  case 0x03:
    TSI_V_V (a, b, c);
    return TSI_A_B;
  case 0x06:
    TSI_V_V (b, c, a);
    return TSI_B_C;
  case 0x05:
    TSI_V_V (c, a, b);
    return TSI_C_A;
  case 0x07:
    TSI_V_V_V (a, b, c);
    return TSI_A_B_C;
  }

  return TSI_OUT;
}
