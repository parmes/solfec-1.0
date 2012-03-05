/*
 * h2d.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * 2D convex hull: Andrew's monotone chain algorithm
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

#ifndef __h2d__
#define __h2d__

typedef struct p2d P2D;

/* auxiliary point list */
struct p2d
{
  double point [3];
  double coord [2];
  P2D *prev, *next;
};


/* input list containts valid 'point' values;
 * these are projected onto (point, normal) to give 'coord';
 * ouput list is made of elements of the input lists;
 * intput list structure gets destroyed upon return */
P2D* h2d (P2D *list, double *point, double *normal);

#endif
