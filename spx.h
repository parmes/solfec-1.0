/*
 * spx.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * simplex integration related routines
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

#ifndef __spx__
#define __spx__

/* simplex integration of few
 * usefull global coordinate functions;
 * self-explanatory
 */

double simplex_J (double *a, double *b, double *c, double *d);

#define simplex_1(J, a, b, c, d) ((J)/6.)

#define simplex_x(J, a, b, c, d) (((J)/24.)*((a) [0] + (b) [0] + (c) [0] + (d) [0]))

#define simplex_y(J, a, b, c, d) (((J)/24.)*((a) [1] + (b) [1] + (c) [1] + (d) [1]))

#define simplex_z(J, a, b, c, d) (((J)/24.)*((a) [2] + (b) [2] + (c) [2] + (d) [2]))

#define simplex_xx(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [0]*(a) [0] + (a) [0]*(b) [0] + (a) [0]*(c) [0] + (a) [0]*(d) [0] +\
(b) [0]*(a) [0] + 2.*(b) [0]*(b) [0] + (b) [0]*(c) [0] + (b) [0]*(d) [0] +\
(c) [0]*(a) [0] + (c) [0]*(b) [0] + 2.*(c) [0]*(c) [0] + (c) [0]*(d) [0] +\
(d) [0]*(a) [0] + (d) [0]*(b) [0] + (d) [0]*(c) [0] + 2.*(d) [0]*(d) [0] \
))

#define simplex_xy(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [0]*(a) [1] + (a) [0]*(b) [1] + (a) [0]*(c) [1] + (a) [0]*(d) [1] +\
(b) [0]*(a) [1] + 2.*(b) [0]*(b) [1] + (b) [0]*(c) [1] + (b) [0]*(d) [1] +\
(c) [0]*(a) [1] + (c) [0]*(b) [1] + 2.*(c) [0]*(c) [1] + (c) [0]*(d) [1] +\
(d) [0]*(a) [1] + (d) [0]*(b) [1] + (d) [0]*(c) [1] + 2.*(d) [0]*(d) [1] \
))

#define simplex_xz(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [0]*(a) [2] + (a) [0]*(b) [2] + (a) [0]*(c) [2] + (a) [0]*(d) [2] +\
(b) [0]*(a) [2] + 2.*(b) [0]*(b) [2] + (b) [0]*(c) [2] + (b) [0]*(d) [2] +\
(c) [0]*(a) [2] + (c) [0]*(b) [2] + 2.*(c) [0]*(c) [2] + (c) [0]*(d) [2] +\
(d) [0]*(a) [2] + (d) [0]*(b) [2] + (d) [0]*(c) [2] + 2.*(d) [0]*(d) [2] \
))

#define simplex_yy(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [1]*(a) [1] + (a) [1]*(b) [1] + (a) [1]*(c) [1] + (a) [1]*(d) [1] +\
(b) [1]*(a) [1] + 2.*(b) [1]*(b) [1] + (b) [1]*(c) [1] + (b) [1]*(d) [1] +\
(c) [1]*(a) [1] + (c) [1]*(b) [1] + 2.*(c) [1]*(c) [1] + (c) [1]*(d) [1] +\
(d) [1]*(a) [1] + (d) [1]*(b) [1] + (d) [1]*(c) [1] + 2.*(d) [1]*(d) [1] \
))

#define simplex_yz(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [1]*(a) [2] + (a) [1]*(b) [2] + (a) [1]*(c) [2] + (a) [1]*(d) [2] +\
(b) [1]*(a) [2] + 2.*(b) [1]*(b) [2] + (b) [1]*(c) [2] + (b) [1]*(d) [2] +\
(c) [1]*(a) [2] + (c) [1]*(b) [2] + 2.*(c) [1]*(c) [2] + (c) [1]*(d) [2] +\
(d) [1]*(a) [2] + (d) [1]*(b) [2] + (d) [1]*(c) [2] + 2.*(d) [1]*(d) [2] \
))

#define simplex_zz(J, a, b, c, d)\
(((J)/120.)*(\
2.*(a) [2]*(a) [2] + (a) [2]*(b) [2] + (a) [2]*(c) [2] + (a) [2]*(d) [2] +\
(b) [2]*(a) [2] + 2.*(b) [2]*(b) [2] + (b) [2]*(c) [2] + (b) [2]*(d) [2] +\
(c) [2]*(a) [2] + (c) [2]*(b) [2] + 2.*(c) [2]*(c) [2] + (c) [2]*(d) [2] +\
(d) [2]*(a) [2] + (d) [2]*(b) [2] + (d) [2]*(c) [2] + 2.*(d) [2]*(d) [2] \
))

#endif
