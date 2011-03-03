/*
 * alg.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * basic operations on scalars, vectors, matrices ...
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

#ifndef __alg__
#define __alg__

/* using geometric tollerance is
 * an Ahilles heal of Solfec */

extern double GEOMETRIC_EPSILON; /* 1.0E-4 by default */

/* global functuins */

/* lexicographical point comparison accounting for the GEOMETRIC_EPSILON */
int POINTS_COMPARE (double *a, double *b);

/* some constants and small,
 * scalar macros follow */

#define ALG_PI 3.14159265358979323846 

#ifndef MIN
  #define MIN(v, w) ((v) < (w) ? (v) : (w))
#endif

#ifndef MAX
  #define MAX(v, w) ((v) > (w) ? (v) : (w))
#endif

#define SGN(v) ((v) > 0 ? 1 : ((v) < 0 ? -1 : 0))
#define ABS(v) ((v) > 0 ? (v) : -(v))

/* comparisons accounting for the GEOMETRIC_EPSILON */
#define LT(x, y) ((x)+GEOMETRIC_EPSILON < (y)-GEOMETRIC_EPSILON)
#define GT(x, y) ((x)-GEOMETRIC_EPSILON > (y)+GEOMETRIC_EPSILON)
#define GE(x, y) ((y) - (x) <= 2.0*GEOMETRIC_EPSILON)
#define LE(x, y) ((x) - (y) <= 2.0*GEOMETRIC_EPSILON)
#define EQ(x, y) (fabs ((x) - (y)) <= 2.0*GEOMETRIC_EPSILON)
#define NE(x, y) (!EQ(X, y))

#define DRAND() ((double) rand () / (double) RAND_MAX)
#define DRANDEXT(x, y) ((x) + ((y) - (x))*DRAND())

/*
 * hashing function macros based on =>
 * M. Teschner, B. Heidelberger, M. Muller, D. Pomeranets, M. Gross,
 * Optimized Spatial Hashing for Collision Detection of Deformable Objects,
 * 2003, Proc. Vision, Modelling, Visualization VMV03, pp. 47-54.
 */

#define C1 73856093
#define C2 19349663
#define C3 83492791
#define INTEGER(x, avgsize) ((int) ((x) / (avgsize)))
#define HASH1(x, hashsize) ABS (((x) * C1) % (hashsize))
#define HASH2(x, y, hashsize) ABS ((((x) * C1) ^ ((y) * C2)) % (hashsize))
#define HASH3(x, y, z, hashsize) ABS ((((x) * C1) ^ ((y) * C2) ^ ((z) * C3)) % (hashsize))

/* various algebraic macros for small
 * vectors and matrices follow below */

#define VECTOR(a, b, c, d)\
{\
  (a) [0] = b;\
  (a) [1] = c;\
  (a) [2] = d;\
}
				      
#define ADD(a, b, c)\
{\
  (c) [0] = (a) [0] + (b) [0];\
  (c) [1] = (a) [1] + (b) [1];\
  (c) [2] = (a) [2] + (b) [2];\
}

#define ADD4(a, b, c)\
{\
  (c) [0] = (a) [0] + (b) [0];\
  (c) [1] = (a) [1] + (b) [1];\
  (c) [2] = (a) [2] + (b) [2];\
  (c) [3] = (a) [3] + (b) [3];\
}

#define ADDMUL(a, eps, b, c)\
{\
  (c) [0] = (a) [0] + (eps)*(b) [0];\
  (c) [1] = (a) [1] + (eps)*(b) [1];\
  (c) [2] = (a) [2] + (eps)*(b) [2];\
}

#define SUB(a, b, c)\
{\
  (c) [0] = (a) [0] - (b) [0];\
  (c) [1] = (a) [1] - (b) [1];\
  (c) [2] = (a) [2] - (b) [2];\
}

#define SUB2(a, b, c)\
{\
  (c) [0] = (a) [0] - (b) [0];\
  (c) [1] = (a) [1] - (b) [1];\
}

#define SUB4(a, b, c)\
{\
  (c) [0] = (a) [0] - (b) [0];\
  (c) [1] = (a) [1] - (b) [1];\
  (c) [2] = (a) [2] - (b) [2];\
  (c) [3] = (a) [3] - (b) [3];\
}

#define SUBMUL(a, eps, b, c)\
{\
  (c) [0] = (a) [0] - (eps)*(b) [0];\
  (c) [1] = (a) [1] - (eps)*(b) [1];\
  (c) [2] = (a) [2] - (eps)*(b) [2];\
}

#define MUL(a, eps, b)\
{\
  (b) [0] = (a) [0] * (eps);\
  (b) [1] = (a) [1] * (eps);\
  (b) [2] = (a) [2] * (eps);\
}

#define AXPY(a, eps, b)\
{\
  (b) [0] += (a) [0] * (eps);\
  (b) [1] += (a) [1] * (eps);\
  (b) [2] += (a) [2] * (eps);\
}

#define MUL2(a, eps, b)\
{\
  (b) [0] = (a) [0] * (eps);\
  (b) [1] = (a) [1] * (eps);\
}

#define DIV(a, eps, b)\
{\
  (b) [0] = (a) [0] / (eps);\
  (b) [1] = (a) [1] / (eps);\
  (b) [2] = (a) [2] / (eps);\
}

#define DIV2(a, eps, b)\
{\
  (b) [0] = (a) [0] / (eps);\
  (b) [1] = (a) [1] / (eps);\
}

#define HADAMARD(a, b, c)\
{\
  (c) [0] = (a) [0] * (b) [0];\
  (c) [1] = (a) [1] * (b) [1];\
  (c) [2] = (a) [2] * (b) [2];\
}

#define DOT(a, b)\
  ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2])

#define DOT2(a, b)\
  ((a)[0]*(b)[0] + (a)[1]*(b)[1])

#define DOT4(a, b) (DOT(a,b) + (a)[3]*(b)[3])

#define DOT5(a, b) (DOT4(a,b) + (a)[4]*(b)[4])

#define LEN(a) (sqrt (DOT (a, a)))

#define LEN2(a) (sqrt (DOT2 (a, a)))

#define LEN4(a) (sqrt (DOT4 (a, a)))

#define PRODUCT(a, b, c)\
{\
  (c) [0] = (a) [1]*(b) [2] - (a) [2]*(b) [1];\
  (c) [1] = (a) [2]*(b) [0] - (a) [0]*(b) [2];\
  (c) [2] = (a) [0]*(b) [1] - (a) [1]*(b) [0];\
}

#define PRODUCTADD(a, b, c)\
{\
  (c) [0] += (a) [1]*(b) [2] - (a) [2]*(b) [1];\
  (c) [1] += (a) [2]*(b) [0] - (a) [0]*(b) [2];\
  (c) [2] += (a) [0]*(b) [1] - (a) [1]*(b) [0];\
}

#define PRODUCTSUB(a, b, c)\
{\
  (c) [0] -= (a) [1]*(b) [2] - (a) [2]*(b) [1];\
  (c) [1] -= (a) [2]*(b) [0] - (a) [0]*(b) [2];\
  (c) [2] -= (a) [0]*(b) [1] - (a) [1]*(b) [0];\
}

#define SCALE(a, b)\
{\
  (a) [0] *= b;\
  (a) [1] *= b;\
  (a) [2] *= b;\
}

#define SCALE2(a, b)\
{\
  (a) [0] *= b;\
  (a) [1] *= b;\
}

#define SCALE4(a, b)\
{\
  (a) [0] *= b;\
  (a) [1] *= b;\
  (a) [2] *= b;\
  (a) [3] *= b;\
}

#define COPY(a, b)\
{\
  (b) [0] = (a) [0];\
  (b) [1] = (a) [1];\
  (b) [2] = (a) [2];\
}

#define COPY2(a, b)\
{\
  (b) [0] = (a) [0];\
  (b) [1] = (a) [1];\
}

#define COPY4(a, b)\
{\
  (b) [0] = (a) [0];\
  (b) [1] = (a) [1];\
  (b) [2] = (a) [2];\
  (b) [3] = (a) [3];\
}

#define COPY6(a, b)\
{\
  (b) [0] = (a) [0];\
  (b) [1] = (a) [1];\
  (b) [2] = (a) [2];\
  (b) [3] = (a) [3];\
  (b) [4] = (a) [4];\
  (b) [5] = (a) [5];\
}

#define ACC(a, b)\
{\
  (b) [0] += (a) [0];\
  (b) [1] += (a) [1];\
  (b) [2] += (a) [2];\
}

#define SCC(a, b)\
{\
  (b) [0] -= (a) [0];\
  (b) [1] -= (a) [1];\
  (b) [2] -= (a) [2];\
}

#define ACCABS(a, b)\
{\
  (b) [0] += fabs ((a) [0]);\
  (b) [1] += fabs ((a) [1]);\
  (b) [2] += fabs ((a) [2]);\
}

#define ACC6(a, b)\
{\
  (b) [0] += (a) [0];\
  (b) [1] += (a) [1];\
  (b) [2] += (a) [2];\
  (b) [3] += (a) [3];\
  (b) [4] += (a) [4];\
  (b) [5] += (a) [5];\
}

#define COPY12(a, b)\
{\
  (b) [0] = (a) [0];\
  (b) [1] = (a) [1];\
  (b) [2] = (a) [2];\
  (b) [3] = (a) [3];\
  (b) [4] = (a) [4];\
  (b) [5] = (a) [5];\
  (b) [6] = (a) [6];\
  (b) [7] = (a) [7];\
  (b) [8] = (a) [8];\
  (b) [9] = (a) [9];\
  (b) [10] = (a) [10];\
  (b) [11] = (a) [11];\
}

#define ACC12(a, b)\
{\
  (b) [0] += (a) [0];\
  (b) [1] += (a) [1];\
  (b) [2] += (a) [2];\
  (b) [3] += (a) [3];\
  (b) [4] += (a) [4];\
  (b) [5] += (a) [5];\
  (b) [6] += (a) [6];\
  (b) [7] += (a) [7];\
  (b) [8] += (a) [8];\
  (b) [9] += (a) [9];\
  (b) [10] += (a) [10];\
  (b) [11] += (a) [11];\
}

#define SET(a, b)\
{\
  (a) [0] = b;\
  (a) [1] = b;\
  (a) [2] = b;\
}

#define SETN(a, n, b) for (double *x = (a), *y = x + (n); x != y; x ++) (*x) = (b)

#define SETRAND(a, b)\
{\
  (a) [0] = DRANDEXT(-b, b);\
  (a) [1] = DRANDEXT(-b, b);\
  (a) [2] = DRANDEXT(-b, b);\
}

#define SET2(a, b)\
{\
  (a) [0] = b;\
  (a) [1] = b;\
}

#define MID(a, b, c)\
{\
  (c) [0] = .5*((a) [0] + (b) [0]);\
  (c) [1] = .5*((a) [1] + (b) [1]);\
  (c) [2] = .5*((a) [2] + (b) [2]);\
}

#define MID3(a, b, c, d)\
{\
  (d) [0] = .333333333333333*((a) [0] + (b) [0] + (c) [0]);\
  (d) [1] = .333333333333333*((a) [1] + (b) [1] + (c) [1]);\
  (d) [2] = .333333333333333*((a) [2] + (b) [2] + (c) [2]);\
}

#define DAMP(beta, a, b, c)\
{\
  (c) [0] = (beta)*(a) [0] + (1.0 - (beta))*(b) [0];\
  (c) [1] = (beta)*(a) [1] + (1.0 - (beta))*(b) [1];\
  (c) [2] = (beta)*(a) [2] + (1.0 - (beta))*(b) [2];\
}

/* triangle (a, b, c) area
 * reads |((a-b)x(c-b))|/2 */
#define TRIANGLE_AREA(a, b, c, area)\
{\
  double ab [3], cb [3], nl [3];\
  SUB (a, b, ab);\
  SUB (c, b, cb);\
  PRODUCT (ab, cb, nl);\
  (area) = 0.5 * LEN (nl);\
}

/* counter clock wise (a, b, c) normal */
#define NORMAL(a, b, c, normal)\
{\
  double ba [3], cb [3];\
  SUB (b, a, ba);\
  SUB (c, b, cb);\
  PRODUCT (ba, cb, normal);\
}

#define NORMALIZE(a)\
{\
  double len = LEN (a);\
  (a) [0] /= len;\
  (a) [1] /= len;\
  (a) [2] /= len;\
}

#define PLANE(pln, pnt) (DOT (pln, pnt) + (pln) [3])

/* p = 4-plane and line segment (a, b) intersection; call only for intersecting entities */
inline static void PLANESEG (double *plane, double *a, double *b, double *p)
{
  double l [3], s;
  
  l [0] = b [0] - a [0];
  l [1] = b [1] - a [1];
  l [2] = b [2] - a [2];
  s = PLANE (plane, a) / DOT (plane, l);
  p [0] = a [0] - l [0] * s;
  p [1] = a [1] - l [1] * s;
  p [2] = a [2] - l [2] * s;
}

#define PLANECOPY(a, b)\
{\
  (b) [0] = (a) [0];\
  (b) [1] = (a) [1];\
  (b) [2] = (a) [2];\
  (b) [3] = (a) [3];\
}

#define MAXABS(a, maximum)\
{\
  double a0 = ABS((a)[0]),\
	 a1 = ABS((a)[1]),\
	 a2 = ABS((a)[2]);\
  (maximum) = a0;\
  if ((maximum) < a1) (maximum) = a1;\
  if ((maximum) < a2) (maximum) = a2;\
}

#define MAXABSIDX(a, imax)\
{\
  double a0 = ABS((a)[0]),\
	 a1 = ABS((a)[1]),\
	 a2 = ABS((a)[2]),\
         maximum = a0;\
  (imax) = 0;\
  if (maximum < a1) {maximum = a1; (imax) = 1;}\
  if (maximum < a2) {maximum = a2; (imax) = 2;}\
}

#define MAXABS2(a, maximum)\
{\
  double a0 = ABS((a)[0]),\
	 a1 = ABS((a)[1]);\
  (maximum) = a0;\
  if ((maximum) < a1) (maximum) = a1;\
}

#define MAXABSN(a, n, maximum)\
{\
  (maximum) = ABS((a)[0]);\
  for (int i = 1; i < (n); i ++)\
  {\
    double ai = ABS((a)[i]);\
    if ((maximum) < ai) (maximum) = ai;\
  }\
}

#define FILTERN(a, n, tol)\
{\
  for (int i = 0; i < (n); i ++)\
  {\
    double ai = ABS((a)[i]);\
    if (ai <= tol) (a)[i] = 0;\
  }\
}

#define NORMALIZE4(a)\
{\
  double len = LEN4 (a);\
  (a) [0] /= len;\
  (a) [1] /= len;\
  (a) [2] /= len;\
  (a) [3] /= len;\
}

#define SCALEDIAG(A, EPS)\
{\
  (A) [0] *= (EPS);\
  (A) [4] *= (EPS);\
  (A) [8] *= (EPS);\
}

#define SET6(A, EPS)\
{\
  (A) [0] = (EPS);\
  (A) [1] = (EPS);\
  (A) [2] = (EPS);\
  (A) [3] = (EPS);\
  (A) [4] = (EPS);\
  (A) [5] = (EPS);\
}

#define SET9(A, EPS)\
{\
  (A) [0] = (EPS);\
  (A) [1] = (EPS);\
  (A) [2] = (EPS);\
  (A) [3] = (EPS);\
  (A) [4] = (EPS);\
  (A) [5] = (EPS);\
  (A) [6] = (EPS);\
  (A) [7] = (EPS);\
  (A) [8] = (EPS);\
}
				  
#define SET12(A, EPS)\
{\
  (A) [0] = (EPS);\
  (A) [1] = (EPS);\
  (A) [2] = (EPS);\
  (A) [3] = (EPS);\
  (A) [4] = (EPS);\
  (A) [5] = (EPS);\
  (A) [6] = (EPS);\
  (A) [7] = (EPS);\
  (A) [8] = (EPS);\
  (A) [9] = (EPS);\
  (A) [10] = (EPS);\
  (A) [11] = (EPS);\
}

#define SCALE6(A, EPS)\
{\
  (A) [0] *= (EPS);\
  (A) [1] *= (EPS);\
  (A) [2] *= (EPS);\
  (A) [3] *= (EPS);\
  (A) [4] *= (EPS);\
  (A) [5] *= (EPS);\
}

#define SCALE9(A, EPS)\
{\
  (A) [0] *= (EPS);\
  (A) [1] *= (EPS);\
  (A) [2] *= (EPS);\
  (A) [3] *= (EPS);\
  (A) [4] *= (EPS);\
  (A) [5] *= (EPS);\
  (A) [6] *= (EPS);\
  (A) [7] *= (EPS);\
  (A) [8] *= (EPS);\
}

#define SCALE12(A, EPS)\
{\
  (A) [0] *= (EPS);\
  (A) [1] *= (EPS);\
  (A) [2] *= (EPS);\
  (A) [3] *= (EPS);\
  (A) [4] *= (EPS);\
  (A) [5] *= (EPS);\
  (A) [6] *= (EPS);\
  (A) [7] *= (EPS);\
  (A) [8] *= (EPS);\
  (A) [9] *= (EPS);\
  (A) [10] *= (EPS);\
  (A) [11] *= (EPS);\
}

#define DOT6(A, B)\
((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2]+\
 (A)[3]*(B)[3]+(A)[4]*(B)[4]+(A)[5]*(B)[5])\

#define LEN6(A) (sqrt (DOT6 (A, A)))

#define DOT9(A, B)\
((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2]+\
 (A)[3]*(B)[3]+(A)[4]*(B)[4]+(A)[5]*(B)[5]+\
 (A)[6]*(B)[6]+(A)[7]*(B)[7]+(A)[8]*(B)[8])\

#define LEN9(A) (sqrt (DOT9 (A, A)))

#define DOT12(A, B)\
((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2]+\
 (A)[3]*(B)[3]+(A)[4]*(B)[4]+(A)[5]*(B)[5]+\
 (A)[6]*(B)[6]+(A)[7]*(B)[7]+(A)[8]*(B)[8]+\
 (A)[9]*(B)[9]+(A)[10]*(B)[10]+(A)[11]*(B)[11])\

#define LEN12(A) (sqrt (DOT12 (A, A)))

#define ADD6(A, B, C)\
{\
  (C) [0] = (A)[0] + (B)[0];\
  (C) [1] = (A)[1] + (B)[1];\
  (C) [2] = (A)[2] + (B)[2];\
  (C) [3] = (A)[3] + (B)[3];\
  (C) [4] = (A)[4] + (B)[4];\
  (C) [5] = (A)[5] + (B)[5];\
}

#define SUB6(A, B, C)\
{\
  (C) [0] = (A)[0] - (B)[0];\
  (C) [1] = (A)[1] - (B)[1];\
  (C) [2] = (A)[2] - (B)[2];\
  (C) [3] = (A)[3] - (B)[3];\
  (C) [4] = (A)[4] - (B)[4];\
  (C) [5] = (A)[5] - (B)[5];\
}

#define ADD12(A, B, C)\
{\
  (C) [0] = (A)[0] + (B)[0];\
  (C) [1] = (A)[1] + (B)[1];\
  (C) [2] = (A)[2] + (B)[2];\
  (C) [3] = (A)[3] + (B)[3];\
  (C) [4] = (A)[4] + (B)[4];\
  (C) [5] = (A)[5] + (B)[5];\
  (C) [6] = (A)[6] + (B)[6];\
  (C) [7] = (A)[7] + (B)[7];\
  (C) [8] = (A)[8] + (B)[8];\
  (C) [9] = (A)[9] + (B)[9];\
  (C) [10] = (A)[10] + (B)[10];\
  (C) [11] = (A)[11] + (B)[11];\
}

#define SUB12(A, B, C)\
{\
  (C) [0] = (A)[0] - (B)[0];\
  (C) [1] = (A)[1] - (B)[1];\
  (C) [2] = (A)[2] - (B)[2];\
  (C) [3] = (A)[3] - (B)[3];\
  (C) [4] = (A)[4] - (B)[4];\
  (C) [5] = (A)[5] - (B)[5];\
  (C) [6] = (A)[6] - (B)[6];\
  (C) [7] = (A)[7] - (B)[7];\
  (C) [8] = (A)[8] - (B)[8];\
  (C) [9] = (A)[9] - (B)[9];\
  (C) [10] = (A)[10] - (B)[10];\
  (C) [11] = (A)[11] - (B)[11];\
}

#define MAX9(A, NORM)\
{\
  double __AUX__;\
  (NORM) = 0.0;\
  __AUX__ = fabs ((A) [0]);\
  if (__AUX__ > (NORM)) (NORM) = __AUX__;\
  __AUX__ = fabs ((A) [1]);\
  if (__AUX__ > (NORM)) (NORM) = __AUX__;\
  __AUX__ = fabs ((A) [2]);\
  if (__AUX__ > (NORM)) (NORM) = __AUX__;\
  __AUX__ = fabs ((A) [3]);\
  if (__AUX__ > (NORM)) (NORM) = __AUX__;\
  __AUX__ = fabs ((A) [4]);\
  if (__AUX__ > (NORM)) (NORM) = __AUX__;\
  __AUX__ = fabs ((A) [5]);\
  if (__AUX__ > (NORM)) (NORM) = __AUX__;\
  __AUX__ = fabs ((A) [6]);\
  if (__AUX__ > (NORM)) (NORM) = __AUX__;\
  __AUX__ = fabs ((A) [7]);\
  if (__AUX__ > (NORM)) (NORM) = __AUX__;\
  __AUX__ = fabs ((A) [8]);\
  if (__AUX__ > (NORM)) (NORM) = __AUX__;\
}

#define HADAMARD9(A, B, C)\
{\
  (C) [0] = (A)[0]*(B)[0];\
  (C) [1] = (A)[1]*(B)[1];\
  (C) [2] = (A)[2]*(B)[2];\
  (C) [3] = (A)[3]*(B)[3];\
  (C) [4] = (A)[4]*(B)[4];\
  (C) [5] = (A)[5]*(B)[5];\
  (C) [6] = (A)[6]*(B)[6];\
  (C) [7] = (A)[7]*(B)[7];\
  (C) [8] = (A)[8]*(B)[8];\
}

#define NNCOPY(A, B)\
{\
  (B) [0] = (A) [0];\
  (B) [1] = (A) [1];\
  (B) [2] = (A) [2];\
  (B) [3] = (A) [3];\
  (B) [4] = (A) [4];\
  (B) [5] = (A) [5];\
  (B) [6] = (A) [6];\
  (B) [7] = (A) [7];\
  (B) [8] = (A) [8];\
}

#define TNCOPY(A, B)\
{\
  (B) [0] = (A) [0];\
  (B) [1] = (A) [3];\
  (B) [2] = (A) [6];\
  (B) [3] = (A) [1];\
  (B) [4] = (A) [4];\
  (B) [5] = (A) [7];\
  (B) [6] = (A) [2];\
  (B) [7] = (A) [5];\
  (B) [8] = (A) [8];\
}

#define NTCOPY(A, B)\
{\
  (B) [0] = (A) [0];\
  (B) [3] = (A) [1];\
  (B) [6] = (A) [2];\
  (B) [1] = (A) [3];\
  (B) [4] = (A) [4];\
  (B) [7] = (A) [5];\
  (B) [2] = (A) [6];\
  (B) [5] = (A) [7];\
  (B) [8] = (A) [8];\
}

#define NNSUB(A, B, C)\
{\
  (C) [0] = (A)[0] - (B)[0];\
  (C) [1] = (A)[1] - (B)[1];\
  (C) [2] = (A)[2] - (B)[2];\
  (C) [3] = (A)[3] - (B)[3];\
  (C) [4] = (A)[4] - (B)[4];\
  (C) [5] = (A)[5] - (B)[5];\
  (C) [6] = (A)[6] - (B)[6];\
  (C) [7] = (A)[7] - (B)[7];\
  (C) [8] = (A)[8] - (B)[8];\
}

#define TNSUB(A, B, C)\
{\
  (C) [0] = (A)[0] - (B)[0];\
  (C) [1] = (A)[3] - (B)[1];\
  (C) [2] = (A)[6] - (B)[2];\
  (C) [3] = (A)[1] - (B)[3];\
  (C) [4] = (A)[4] - (B)[4];\
  (C) [5] = (A)[7] - (B)[5];\
  (C) [6] = (A)[2] - (B)[6];\
  (C) [7] = (A)[5] - (B)[7];\
  (C) [8] = (A)[8] - (B)[8];\
}

#define NTSUB(A, B, C)\
{\
  (C) [0] = (A)[0] - (B)[0];\
  (C) [1] = (A)[1] - (B)[3];\
  (C) [2] = (A)[2] - (B)[6];\
  (C) [3] = (A)[3] - (B)[1];\
  (C) [4] = (A)[4] - (B)[4];\
  (C) [5] = (A)[5] - (B)[7];\
  (C) [6] = (A)[6] - (B)[2];\
  (C) [7] = (A)[7] - (B)[5];\
  (C) [8] = (A)[8] - (B)[8];\
}

#define NNADD(A, B, C)\
{\
  (C) [0] = (A)[0] + (B)[0];\
  (C) [1] = (A)[1] + (B)[1];\
  (C) [2] = (A)[2] + (B)[2];\
  (C) [3] = (A)[3] + (B)[3];\
  (C) [4] = (A)[4] + (B)[4];\
  (C) [5] = (A)[5] + (B)[5];\
  (C) [6] = (A)[6] + (B)[6];\
  (C) [7] = (A)[7] + (B)[7];\
  (C) [8] = (A)[8] + (B)[8];\
}

#define TNADD(A, B, C)\
{\
  (C) [0] = (A)[0] + (B)[0];\
  (C) [1] = (A)[3] + (B)[1];\
  (C) [2] = (A)[6] + (B)[2];\
  (C) [3] = (A)[1] + (B)[3];\
  (C) [4] = (A)[4] + (B)[4];\
  (C) [5] = (A)[7] + (B)[5];\
  (C) [6] = (A)[2] + (B)[6];\
  (C) [7] = (A)[5] + (B)[7];\
  (C) [8] = (A)[8] + (B)[8];\
}

#define NTADD(A, B, C)\
{\
  (C) [0] = (A)[0] + (B)[0];\
  (C) [1] = (A)[1] + (B)[3];\
  (C) [2] = (A)[2] + (B)[6];\
  (C) [3] = (A)[3] + (B)[1];\
  (C) [4] = (A)[4] + (B)[4];\
  (C) [5] = (A)[5] + (B)[7];\
  (C) [6] = (A)[6] + (B)[2];\
  (C) [7] = (A)[7] + (B)[5];\
  (C) [8] = (A)[8] + (B)[8];\
}

#define NNSUBMUL(A, EPS, B, C)\
{\
  (C) [0] = (A)[0] - (EPS)*(B)[0];\
  (C) [1] = (A)[1] - (EPS)*(B)[1];\
  (C) [2] = (A)[2] - (EPS)*(B)[2];\
  (C) [3] = (A)[3] - (EPS)*(B)[3];\
  (C) [4] = (A)[4] - (EPS)*(B)[4];\
  (C) [5] = (A)[5] - (EPS)*(B)[5];\
  (C) [6] = (A)[6] - (EPS)*(B)[6];\
  (C) [7] = (A)[7] - (EPS)*(B)[7];\
  (C) [8] = (A)[8] - (EPS)*(B)[8];\
}

#define TNSUBMUL(A, EPS, B, C)\
{\
  (C) [0] = (A)[0] - (EPS)*(B)[0];\
  (C) [1] = (A)[3] - (EPS)*(B)[1];\
  (C) [2] = (A)[6] - (EPS)*(B)[2];\
  (C) [3] = (A)[1] - (EPS)*(B)[3];\
  (C) [4] = (A)[4] - (EPS)*(B)[4];\
  (C) [5] = (A)[7] - (EPS)*(B)[5];\
  (C) [6] = (A)[2] - (EPS)*(B)[6];\
  (C) [7] = (A)[5] - (EPS)*(B)[7];\
  (C) [8] = (A)[8] - (EPS)*(B)[8];\
}

#define NTSUBMUL(A, EPS, B, C)\
{\
  (C) [0] = (A)[0] - (EPS)*(B)[0];\
  (C) [1] = (A)[1] - (EPS)*(B)[3];\
  (C) [2] = (A)[2] - (EPS)*(B)[6];\
  (C) [3] = (A)[3] - (EPS)*(B)[1];\
  (C) [4] = (A)[4] - (EPS)*(B)[4];\
  (C) [5] = (A)[5] - (EPS)*(B)[7];\
  (C) [6] = (A)[6] - (EPS)*(B)[2];\
  (C) [7] = (A)[7] - (EPS)*(B)[5];\
  (C) [8] = (A)[8] - (EPS)*(B)[8];\
}

#define NNADDMUL(A, EPS, B, C)\
{\
  (C) [0] = (A)[0] + (EPS)*(B)[0];\
  (C) [1] = (A)[1] + (EPS)*(B)[1];\
  (C) [2] = (A)[2] + (EPS)*(B)[2];\
  (C) [3] = (A)[3] + (EPS)*(B)[3];\
  (C) [4] = (A)[4] + (EPS)*(B)[4];\
  (C) [5] = (A)[5] + (EPS)*(B)[5];\
  (C) [6] = (A)[6] + (EPS)*(B)[6];\
  (C) [7] = (A)[7] + (EPS)*(B)[7];\
  (C) [8] = (A)[8] + (EPS)*(B)[8];\
}

#define TNADDMUL(A, EPS, B, C)\
{\
  (C) [0] = (A)[0] + (EPS)*(B)[0];\
  (C) [1] = (A)[3] + (EPS)*(B)[1];\
  (C) [2] = (A)[6] + (EPS)*(B)[2];\
  (C) [3] = (A)[1] + (EPS)*(B)[3];\
  (C) [4] = (A)[4] + (EPS)*(B)[4];\
  (C) [5] = (A)[7] + (EPS)*(B)[5];\
  (C) [6] = (A)[2] + (EPS)*(B)[6];\
  (C) [7] = (A)[5] + (EPS)*(B)[7];\
  (C) [8] = (A)[8] + (EPS)*(B)[8];\
}

#define NTADDMUL(A, EPS, B, C)\
{\
  (C) [0] = (A)[0] + (EPS)*(B)[0];\
  (C) [1] = (A)[1] + (EPS)*(B)[3];\
  (C) [2] = (A)[2] + (EPS)*(B)[6];\
  (C) [3] = (A)[3] + (EPS)*(B)[1];\
  (C) [4] = (A)[4] + (EPS)*(B)[4];\
  (C) [5] = (A)[5] + (EPS)*(B)[7];\
  (C) [6] = (A)[6] + (EPS)*(B)[2];\
  (C) [7] = (A)[7] + (EPS)*(B)[5];\
  (C) [8] = (A)[8] + (EPS)*(B)[8];\
}

#define NNMUL(A, B, C)\
{\
 (C) [0] = (A)[0]*(B)[0]+(A)[3]*(B)[1]+(A)[6]*(B)[2];\
 (C) [1] = (A)[1]*(B)[0]+(A)[4]*(B)[1]+(A)[7]*(B)[2];\
 (C) [2] = (A)[2]*(B)[0]+(A)[5]*(B)[1]+(A)[8]*(B)[2];\
 (C) [3] = (A)[0]*(B)[3]+(A)[3]*(B)[4]+(A)[6]*(B)[5];\
 (C) [4] = (A)[1]*(B)[3]+(A)[4]*(B)[4]+(A)[7]*(B)[5];\
 (C) [5] = (A)[2]*(B)[3]+(A)[5]*(B)[4]+(A)[8]*(B)[5];\
 (C) [6] = (A)[0]*(B)[6]+(A)[3]*(B)[7]+(A)[6]*(B)[8];\
 (C) [7] = (A)[1]*(B)[6]+(A)[4]*(B)[7]+(A)[7]*(B)[8];\
 (C) [8] = (A)[2]*(B)[6]+(A)[5]*(B)[7]+(A)[8]*(B)[8];\
}

#define TNMUL(A, B, C)\
{\
 (C) [0] = (A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2];\
 (C) [1] = (A)[3]*(B)[0]+(A)[4]*(B)[1]+(A)[5]*(B)[2];\
 (C) [2] = (A)[6]*(B)[0]+(A)[7]*(B)[1]+(A)[8]*(B)[2];\
 (C) [3] = (A)[0]*(B)[3]+(A)[1]*(B)[4]+(A)[2]*(B)[5];\
 (C) [4] = (A)[3]*(B)[3]+(A)[4]*(B)[4]+(A)[5]*(B)[5];\
 (C) [5] = (A)[6]*(B)[3]+(A)[7]*(B)[4]+(A)[8]*(B)[5];\
 (C) [6] = (A)[0]*(B)[6]+(A)[1]*(B)[7]+(A)[2]*(B)[8];\
 (C) [7] = (A)[3]*(B)[6]+(A)[4]*(B)[7]+(A)[5]*(B)[8];\
 (C) [8] = (A)[6]*(B)[6]+(A)[7]*(B)[7]+(A)[8]*(B)[8];\
}

#define NTMUL(A, B, C)\
{\
 (C) [0] = (A)[0]*(B)[0]+(A)[3]*(B)[3]+(A)[6]*(B)[6];\
 (C) [1] = (A)[1]*(B)[0]+(A)[4]*(B)[3]+(A)[7]*(B)[6];\
 (C) [2] = (A)[2]*(B)[0]+(A)[5]*(B)[3]+(A)[8]*(B)[6];\
 (C) [3] = (A)[0]*(B)[1]+(A)[3]*(B)[4]+(A)[6]*(B)[7];\
 (C) [4] = (A)[1]*(B)[1]+(A)[4]*(B)[4]+(A)[7]*(B)[7];\
 (C) [5] = (A)[2]*(B)[1]+(A)[5]*(B)[4]+(A)[8]*(B)[7];\
 (C) [6] = (A)[0]*(B)[2]+(A)[3]*(B)[5]+(A)[6]*(B)[8];\
 (C) [7] = (A)[1]*(B)[2]+(A)[4]*(B)[5]+(A)[7]*(B)[8];\
 (C) [8] = (A)[2]*(B)[2]+(A)[5]*(B)[5]+(A)[8]*(B)[8];\
}

#define NVMUL(A, B, C)\
{\
 (C) [0] = (A)[0]*(B)[0]+(A)[3]*(B)[1]+(A)[6]*(B)[2];\
 (C) [1] = (A)[1]*(B)[0]+(A)[4]*(B)[1]+(A)[7]*(B)[2];\
 (C) [2] = (A)[2]*(B)[0]+(A)[5]*(B)[1]+(A)[8]*(B)[2];\
}

#define NVMUL2(A, B, C)\
{\
 (C) [0] = (A)[0]*(B)[0]+(A)[2]*(B)[1];\
 (C) [1] = (A)[1]*(B)[0]+(A)[3]*(B)[1];\
}

#define NVMUL4(A, B, C)\
{\
 (C) [0] = (A)[0]*(B)[0]+(A)[4]*(B)[1]+(A)[8]*(B)[2]+(A)[12]*(B)[3];\
 (C) [1] = (A)[1]*(B)[0]+(A)[5]*(B)[1]+(A)[9]*(B)[2]+(A)[13]*(B)[3];\
 (C) [2] = (A)[2]*(B)[0]+(A)[6]*(B)[1]+(A)[10]*(B)[2]+(A)[14]*(B)[3];\
 (C) [3] = (A)[3]*(B)[0]+(A)[7]*(B)[1]+(A)[11]*(B)[2]+(A)[15]*(B)[3];\
}

#define TVMUL(A, B, C)\
{\
 (C) [0] = (A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2];\
 (C) [1] = (A)[3]*(B)[0]+(A)[4]*(B)[1]+(A)[5]*(B)[2];\
 (C) [2] = (A)[6]*(B)[0]+(A)[7]*(B)[1]+(A)[8]*(B)[2];\
}

#define NVADDMUL(C, A, B, D)\
{\
 (D) [0] = (C)[0] + ((A)[0]*(B)[0]+(A)[3]*(B)[1]+(A)[6]*(B)[2]);\
 (D) [1] = (C)[1] + ((A)[1]*(B)[0]+(A)[4]*(B)[1]+(A)[7]*(B)[2]);\
 (D) [2] = (C)[2] + ((A)[2]*(B)[0]+(A)[5]*(B)[1]+(A)[8]*(B)[2]);\
}

#define TVADDMUL(C, A, B, D)\
{\
 (D) [0] = (C)[0] + ((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2]);\
 (D) [1] = (C)[1] + ((A)[3]*(B)[0]+(A)[4]*(B)[1]+(A)[5]*(B)[2]);\
 (D) [2] = (C)[2] + ((A)[6]*(B)[0]+(A)[7]*(B)[1]+(A)[8]*(B)[2]);\
}

#define NVSUBMUL(C, A, B, D)\
{\
 (D) [0] = (C)[0] - ((A)[0]*(B)[0]+(A)[3]*(B)[1]+(A)[6]*(B)[2]);\
 (D) [1] = (C)[1] - ((A)[1]*(B)[0]+(A)[4]*(B)[1]+(A)[7]*(B)[2]);\
 (D) [2] = (C)[2] - ((A)[2]*(B)[0]+(A)[5]*(B)[1]+(A)[8]*(B)[2]);\
}

#define TVSUBMUL(C, A, B, D)\
{\
 (D) [0] = (C)[0] - ((A)[0]*(B)[0]+(A)[1]*(B)[1]+(A)[2]*(B)[2]);\
 (D) [1] = (C)[1] - ((A)[3]*(B)[0]+(A)[4]*(B)[1]+(A)[5]*(B)[2]);\
 (D) [2] = (C)[2] - ((A)[6]*(B)[0]+(A)[7]*(B)[1]+(A)[8]*(B)[2]);\
}

#define TRACE2(A) ((A)[0] + (A)[3])

#define TRACE(A) ((A)[0] + (A)[4] + (A)[8])

#define IDENTITY(A)\
{\
  (A) [0] = 1;\
  (A) [1] = 0;\
  (A) [2] = 0;\
  (A) [3] = 0;\
  (A) [4] = 1;\
  (A) [5] = 0;\
  (A) [6] = 0;\
  (A) [7] = 0;\
  (A) [8] = 1;\
}

#define DET2(F) (((F)[0]*(F)[3]-(F)[1]*(F)[2]))

#define INVERT2(F, INV, DET)\
if (((DET)=((F)[0]*(F)[3]-(F)[1]*(F)[2])) != 0.0)\
{\
  (INV) [0] =  (F)[3] / (DET);\
  (INV) [1] = -(F)[1] / (DET);\
  (INV) [2] = -(F)[2] / (DET);\
  (INV) [3] =  (F)[0] / (DET);\
}

#define DET(F) ((F)[0]*((F)[4]*(F)[8]-(F)[5]*(F)[7])+\
                (F)[3]*((F)[2]*(F)[7]-(F)[1]*(F)[8])+\
                (F)[6]*((F)[1]*(F)[5]-(F)[2]*(F)[4]))

#define INVERT(F, INV, DET)\
if (((DET) =\
 ((F)[0]*((F)[4]*(F)[8]-(F)[5]*(F)[7])+\
  (F)[3]*((F)[2]*(F)[7]-(F)[1]*(F)[8])+\
 ((F)[1]*(F)[5]-(F)[2]*(F)[4])*(F)[6])) != 0.0)\
{\
  (DET) = 1.0 / (DET);\
  (INV) [0] = ((F)[4]*(F)[8]-(F)[5]*(F)[7])*(DET);\
  (INV) [1] = ((F)[2]*(F)[7]-(F)[1]*(F)[8])*(DET);\
  (INV) [2] = ((F)[1]*(F)[5]-(F)[2]*(F)[4])*(DET);\
  (INV) [3] = ((F)[5]*(F)[6]-(F)[3]*(F)[8])*(DET);\
  (INV) [4] = ((F)[0]*(F)[8]-(F)[2]*(F)[6])*(DET);\
  (INV) [5] = ((F)[2]*(F)[3]-(F)[0]*(F)[5])*(DET);\
  (INV) [6] = ((F)[3]*(F)[7]-(F)[4]*(F)[6])*(DET);\
  (INV) [7] = ((F)[1]*(F)[6]-(F)[0]*(F)[7])*(DET);\
  (INV) [8] = ((F)[0]*(F)[4]-(F)[1]*(F)[3])*(DET);\
  (DET) = 1.0 / (DET);\
}

#define SYMM(A, SYMM)\
{\
  (SYMM) [0] = (A) [0];\
  (SYMM) [1] = 0.5 * ((A) [1] + (A) [3]);\
  (SYMM) [2] = 0.5 * ((A) [2] + (A) [6]);\
  (SYMM) [3] = (SYMM) [1];\
  (SYMM) [4] = (A) [4];\
  (SYMM) [5] = 0.5 * ((A) [5] + (A) [7]);\
  (SYMM) [6] = (SYMM) [2];\
  (SYMM) [7] = (SYMM) [5];\
  (SYMM) [8] = (A) [8];\
}

#define SKEW(A, SKEW)\
{\
  (SKEW) [0] = 0;\
  (SKEW) [1] = 0.5 * ((A) [1] - (A) [3]);\
  (SKEW) [2] = -0.5 * ((A) [2] - (A) [6]);\
  (SKEW) [3] = -(SKEW) [1];\
  (SKEW) [4] = 0;\
  (SKEW) [5] = 0.5 * ((A) [5] - (A) [7]);\
  (SKEW) [6] = -(SKEW) [2];\
  (SKEW) [7] = -(SKEW) [5];\
  (SKEW) [8] = 0;\
}

#define SKEWVEC(A, VSKEW)\
{\
  (VSKEW) [0] = 0.5 * ((A) [5] - (A) [7]);\
  (VSKEW) [1] = 0.5 * ((A) [6] - (A) [2]);\
  (VSKEW) [2] = 0.5 * ((A) [1] - (A) [3]);\
}

#define VECSKEW(VSKEW, A)\
{\
  (A) [0] = 0;\
  (A) [1] = (VSKEW) [2];\
  (A) [2] = -(VSKEW) [1];\
  (A) [3] = -(VSKEW) [2];\
  (A) [4] = 0;\
  (A) [5] = (VSKEW) [0];\
  (A) [6] = (VSKEW) [1];\
  (A) [7] = -(VSKEW) [0];\
  (A) [8] = 0;\
}

#define SYMMLOTR(SYMM,LO)\
  (LO)[0] = (SYMM)[0], (LO)[1] = (SYMM)[1], (LO)[2] = (SYMM)[2],\
  (LO)[3] = (SYMM)[4], (LO)[4] = (SYMM)[5], (LO)[5] = (SYMM)[8]

#define LOTRSYMM(LO,SYMM)\
  (SYMM)[0] = (LO)[0], (SYMM)[1] = (LO)[1], (SYMM)[2] = (LO)[2],\
  (SYMM)[3] = (LO)[1], (SYMM)[4] = (LO)[3], (SYMM)[5] = (LO)[4],\
  (SYMM)[6] = (LO)[2], (SYMM)[7] = (LO)[4], (SYMM)[8] = (LO)[5]

#define DIADIC(a, b, C)\
{\
  (C) [0] = (a) [0]*(b) [0];\
  (C) [1] = (a) [1]*(b) [0];\
  (C) [2] = (a) [2]*(b) [0];\
  (C) [3] = (a) [0]*(b) [1];\
  (C) [4] = (a) [1]*(b) [1];\
  (C) [5] = (a) [2]*(b) [1];\
  (C) [6] = (a) [0]*(b) [2];\
  (C) [7] = (a) [1]*(b) [2];\
  (C) [8] = (a) [2]*(b) [2];\
}

/* polar decomposition of 'F' with 'EPS' tollerance,
 * outputing 'R' and 'U' computed after 'ITERS' iterations */
#define POLAR(F, EPS, R, U, ITERS)\
{\
  double __TMP__ [9],\
	 __L2__, __MX__,\
	 __IL2__, __IMX__,\
	 __GAMMA__, __DET__,\
	 __EPS2__ = (EPS)*(EPS);\
  NNCOPY (F, R);\
  (ITERS) = 0;\
  do\
  {\
    (ITERS) ++;\
    NNCOPY (R, U);\
    INVERT (R, __TMP__, __DET__);\
    if (__DET__ == 0.0)\
    { (ITERS) = -1; break; }\
    __L2__ = LEN9 (R);\
    NNMAX (R, __MX__);\
    __IL2__ = LEN9 (__TMP__);\
    NNMAX (__TMP__, __IMX__);\
    __GAMMA__ = sqrt(sqrt((__IL2__*__IMX__)/\
	(__L2__*__MX__))) * 0.5;\
    NNSCALE (R, __GAMMA__);\
    __GAMMA__ = 0.25 / __GAMMA__;\
    NNSCALE (__TMP__, __GAMMA__);\
    NTADD (R, __TMP__, R);\
    NNSUB (R, U, __TMP__);\
  } while (NNDOT (__TMP__, __TMP__) >\
    __EPS2__ * NNDOT (U, U));\
  TNMUL (R, F, __TMP__);\
  SYMM (__TMP__, U);\
}

/* ten degrees squared in radians =>
 * below ten degrees the exponential map
 * below will use Taylor expansions */
#define DEG_10_SQ 3.0461741978671E-02

/* compute an exponential map; given
 * the rotationa vector 'VSKEW' output
 * the orthogonal rotation matrix 'R' */
#define EXPMAP(VSKEW, R)\
{\
  double __ANG_2__, __SIN_X__,\
    __1_COS_XX__, __0__, __1__, __2__,\
    __01__, __02__, __12__, __S0__,\
    __S1__, __S2__;\
  __ANG_2__ = DOT(VSKEW, VSKEW);\
  if (__ANG_2__ < DEG_10_SQ)\
  {\
  __SIN_X__ = 1.0 +\
    (-1.666666666666667E-1 +\
    (8.333333333333333E-3 +\
    (-1.984126984126984E-4 +\
    (2.755731922398589E-6 +\
    (-2.505210838544172E-8 +\
     1.605904383682161E-10 * __ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__;\
  __1_COS_XX__ = 0.5 +\
    (-4.166666666666667E-2 +\
    (1.388888888888889E-3 +\
    (-2.480158730158730E-5 +\
    (2.755731922398589E-7 +\
    (-2.087675698786810E-9 +\
     1.147074559772972E-11 * __ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__;\
  }\
  else\
  {\
    __0__ = __ANG_2__;\
    __ANG_2__ = sqrt (__ANG_2__);\
    __SIN_X__ = sin (__ANG_2__) / __ANG_2__;\
    __1_COS_XX__ = (1.0 - cos (__ANG_2__)) / __0__;\
  }\
  __0__ = (VSKEW) [0] * (VSKEW) [0];\
  __1__ = (VSKEW) [1] * (VSKEW) [1];\
  __2__ = (VSKEW) [2] * (VSKEW) [2];\
  __01__ = (VSKEW) [0] * (VSKEW) [1];\
  __02__ = (VSKEW) [0] * (VSKEW) [2];\
  __12__ = (VSKEW) [1] * (VSKEW) [2];\
  __S0__ = __SIN_X__ * (VSKEW) [0];\
  __S1__ = __SIN_X__ * (VSKEW) [1];\
  __S2__ = __SIN_X__ * (VSKEW) [2];\
  (R) [0] = - __1_COS_XX__*(__2__+__1__);\
  (R) [1] = __1_COS_XX__*__01__;\
  (R) [2] = __1_COS_XX__*__02__;\
  (R) [3] = (R) [1];\
  (R) [4] = - __1_COS_XX__*(__2__+__0__);\
  (R) [5] = __1_COS_XX__*__12__;\
  (R) [6] = (R) [2];\
  (R) [7] = (R) [5];\
  (R) [8] = - __1_COS_XX__*(__1__+__0__);\
  (R) [0] += 1.0;\
  (R) [1] += __S2__;\
  (R) [2] -= __S1__;\
  (R) [3] -= __S2__;\
  (R) [4] += 1.0;\
  (R) [5] += __S0__;\
  (R) [6] += __S1__;\
  (R) [7] -= __S0__;\
  (R) [8] += 1.0;\
}

/* compute differential of the exponential map =>
 * this is useful for moving vectors between tangent spaces
 * spanned at identity and 'VECSKEW' */
#define DEXPMAP(VSKEW, R)\
{\
  double __ANG_2__, __1_XX_SIN_XXX__,\
    __1_COS_XX__, __0__, __1__, __2__,\
    __S0__, __S1__, __S2__;\
  __ANG_2__ = DOT(VSKEW, VSKEW);\
  if (__ANG_2__ < DEG_10_SQ)\
  {\
  __1_XX_SIN_XXX__ = \
      1.666666666666667E-1 +\
    (-8.333333333333333E-3 +\
    (1.984126984126984E-4 +\
    (-2.755731922398589E-6 +\
    (2.505210838544172E-8 +\
    (-1.605904383682161E-10+\
      7.647163731819816E-13 * __ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__;\
  __1_COS_XX__ = 0.5 +\
    (-4.166666666666667E-2 +\
    (1.388888888888889E-3 +\
    (-2.480158730158730E-5 +\
    (2.755731922398589E-7 +\
    (-2.087675698786810E-9 +\
     1.147074559772972E-11 * __ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__\
    )*__ANG_2__;\
  }\
  else\
  {\
    __0__ = __ANG_2__;\
    __ANG_2__ = sqrt (__ANG_2__);\
    __1_XX_SIN_XXX__ = (1.0 / __0__) - (sin (__ANG_2__) / (__0__*__ANG_2__));\
    __1_COS_XX__ = (1.0 - cos (__ANG_2__)) / __0__;\
  }\
  __0__ = (VSKEW)[0]*(VSKEW)[0];\
  __1__ = (VSKEW)[1]*(VSKEW)[1];\
  __2__ = (VSKEW)[2]*(VSKEW)[2];\
  __S0__ = __1_COS_XX__ * (VSKEW) [0];\
  __S1__ = __1_COS_XX__ * (VSKEW) [1];\
  __S2__ = __1_COS_XX__ * (VSKEW) [2];\
  (R) [0] = __1_XX_SIN_XXX__ * (-__2__-__1__);\
  (R) [1] = __1_XX_SIN_XXX__ * (VSKEW)[0]*(VSKEW)[1];\
  (R) [2] = __1_XX_SIN_XXX__ * (VSKEW)[0]*(VSKEW)[2];\
  (R) [3] = (R) [1];\
  (R) [4] = __1_XX_SIN_XXX__ * (-__2__-__0__); \
  (R) [5] =  __1_XX_SIN_XXX__ * (VSKEW)[1]*(VSKEW)[2];\
  (R) [6] = (R) [2];\
  (R) [7] = (R) [5];\
  (R) [8] = __1_XX_SIN_XXX__ * (-__1__-__0__);\
  (R) [0] += 1.0;\
  (R) [1] += __S2__;\
  (R) [2] -= __S1__;\
  (R) [3] -= __S2__;\
  (R) [4] += 1.0;\
  (R) [5] += __S0__;\
  (R) [6] += __S1__;\
  (R) [7] -= __S0__;\
  (R) [8] += 1.0;\
}

/* partial derivatices of the exponential
 * map => useful for Newton iterations */
#define EXPMAP123(OMEGA, R0, R1, R2)\
{\
  double x,\
         a,\
         b,\
         c,\
         b_ac,\
         ac_2cc_2bc,\
         c_b,\
	 aux [9],\
         tmp [3];\
  x = DOT(OMEGA, OMEGA);\
  if (x < 1E-16 * DEG_10_SQ)\
  {\
    R0[0]=R0[1]=R0[2]=R0[3]=R0[4]=R0[6]=R0[8]=0;\
    R1[0]=R1[1]=R1[3]=R1[4]=R1[5]=R1[7]=R1[8]=0;\
    R2[0]=R2[2]=R2[4]=R2[5]=R2[6]=R2[7]=R2[8]=0;\
    R0[5]=1.0; R0[7]=-1.0;\
    R1[2]=-1.0; R1[6]=1.0;\
    R2[1]=1.0; R2[3]=-1.0;\
  }\
  else\
  {\
    c = 1.0 / x;\
    if (x < DEG_10_SQ)\
    {\
      a = 1.0 +\
      (-1.666666666666667E-1 +\
      (8.333333333333333E-3 +\
      (-1.984126984126984E-4 +\
      (2.755731922398589E-6 +\
      (-2.505210838544172E-8 +\
       1.605904383682161E-10 * x\
      )* x\
      )* x\
      )* x\
      )* x\
      )* x;\
      b = c - 0.5 +\
      (4.166666666666667E-2 +\
      (-1.388888888888889E-3 +\
      (2.480158730158730E-5 +\
      (-2.755731922398589E-7 +\
      (2.087675698786810E-9 -\
       1.147074559772972E-11 * x\
      )* x\
      )* x\
      )* x\
      )* x\
      )* x;\
    }\
    else\
    {\
      x = sqrt (x);\
      a = sin (x) / x;\
      b = c * cos (x);\
    }\
    b_ac = b - a*c;\
    ac_2cc_2bc = (a - 2.0 * (c - b)) * c;\
    c_b = c - b;\
    DIADIC (OMEGA, OMEGA, aux);\
    SCALE9 (aux, ac_2cc_2bc);\
    COPY (OMEGA, tmp);\
    SCALE (tmp, b_ac);\
    aux [0] -= a;\
    aux [1] += tmp[2];\
    aux [2] -= tmp[1];\
    aux [3] -= tmp[2];\
    aux [4] -=  a;\
    aux [5] += tmp[0];\
    aux [6] += tmp[1];\
    aux [7] -= tmp[0];\
    aux [8] -= a;\
    COPY (OMEGA, tmp);\
    SCALE (tmp, c_b);\
    (R0) [0] = aux[0] * (OMEGA)[0] + 2.0 * tmp[0];\
    (R0) [1] = aux[1] * (OMEGA)[0] + tmp[1];\
    (R0) [2] = aux[2] * (OMEGA)[0] + tmp[2];\
    (R0) [3] = aux[3] * (OMEGA)[0] + tmp[1];\
    (R0) [4] = aux[4] * (OMEGA)[0];\
    (R0) [5] = aux[5] * (OMEGA)[0] + a;\
    (R0) [6] = aux[6] * (OMEGA)[0] + tmp[2];\
    (R0) [7] = aux[7] * (OMEGA)[0] - a;\
    (R0) [8] = aux[8] * (OMEGA)[0];\
    (R1) [0] = aux[0] * (OMEGA)[1];\
    (R1) [1] = aux[1] * (OMEGA)[1] + tmp[0];\
    (R1) [2] = aux[2] * (OMEGA)[1] - a;\
    (R1) [3] = aux[3] * (OMEGA)[1] + tmp[0];\
    (R1) [4] = aux[4] * (OMEGA)[1] + 2.0 * tmp[1];\
    (R1) [5] = aux[5] * (OMEGA)[1] + tmp[2];\
    (R1) [6] = aux[6] * (OMEGA)[1] + a;\
    (R1) [7] = aux[7] * (OMEGA)[1] + tmp[2];\
    (R1) [8] = aux[8] * (OMEGA)[1];\
    (R2) [0] = aux[0] * (OMEGA)[2];\
    (R2) [1] = aux[1] * (OMEGA)[2] + a;\
    (R2) [2] = aux[2] * (OMEGA)[2] + tmp[0];\
    (R2) [3] = aux[3] * (OMEGA)[2] - a;\
    (R2) [4] = aux[4] * (OMEGA)[2];\
    (R2) [5] = aux[5] * (OMEGA)[2] + tmp[1];\
    (R2) [6] = aux[6] * (OMEGA)[2] + tmp[0];\
    (R2) [7] = aux[7] * (OMEGA)[2] + tmp[1];\
    (R2) [8] = aux[8] * (OMEGA)[2] + 2.0 * tmp[2];\
  }\
}

/* compute Mises norm of a Cauchy stress */
#define MISES(s, v)\
  do {\
  double a = (s [0] - s [1])*(s [0] - s [1]);\
  double b = (s [0] - s [2])*(s [0] - s [2]);\
  double c = (s [2] - s [1])*(s [2] - s [1]);\
  double d = 6. * (s [3]*s [3] + s [4]*s [4] + s [5]*s [5]);\
  v = .707106781186548 * sqrt (a + b + c + d);\
  } while (0)

#define PROJECT_POINT_ON_LINE(point, line_point, line_direction, projection)\
{\
  double dif [3], dot1, dot2, eps;\
  \
  dot1 = DOT (line_direction, line_direction);\
  if (dot1 == 0.0) COPY (line_point, projection);\
  SUB (point, line_point, dif);\
  dot2 = DOT (dif, line_direction);\
  eps = dot2 / dot1;\
  ADDMUL (line_point, eps, line_direction, projection);\
}

#endif
