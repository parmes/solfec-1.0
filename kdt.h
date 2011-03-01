/*
 * kdt.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * kd-tree based map
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

#ifndef __kdt__
#define __kdt__

typedef struct kdt KDT;

struct kdt
{
  double p [3]; /* splitting point */

  int d; /* splitting dimension */

  KDT *l, *r; /* left and right branch */

  int n; /* leaf data items */

  void **data; /* leaf data */
};

/* create kd-tree for n points separated by epsilon
 * margin; returns NULL for coincident points */
KDT* KDT_Create (int n, double *p, double epsilon);

/* drop data down the kd-tree */
void KDT_Drop (KDT *kd, double *extents, void *data);

/* pick data for a point; free buffer after use */
void KDT_Pick (KDT *kd, double *p, void ***data, int *n);

/* return nearest point in kd-tree */
double* KDT_Nearest (KDT *kd, double *p);

/* destroy kd-tree */
void KDT_Destroy (KDT *kd);

#endif
