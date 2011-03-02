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
  double p [3]; /* point or invalid for a leaf */

  int d; /* splitting dimension or -1 for a leaf */

  KDT *u, *l, *r; /* upper, left and right branch */

  int n; /* tree node index or leaf data items count */

  void **data; /* leaf data items */
};

/* create kd-tree for n points; epsilon separation is ensured
 * between the input points and the remaining points are filtered our */
KDT* KDT_Create (int n, double *p, double epsilon);

/* drop data down the kd-tree */
void KDT_Drop (KDT *kd, double *extents, void *data);

/* pick data for a point; free buffer after use */
void KDT_Pick (KDT *kd, double *p, void ***data, int *n);

/* return nearest node in kd-tree within epsilon radius */
KDT* KDT_Nearest (KDT *kd, double *p, double epsilon);

/* return the number kd-tree nodes; note that kd->n indices
 * become valid for tree nodes only after KDT_Size was called */
int KDT_Size (KDT *kd);

/* first node */
KDT* KDT_First (KDT *kd);

/* next node */
KDT* KDT_Next (KDT *kd);

/* destroy kd-tree */
void KDT_Destroy (KDT *kd);

#endif
