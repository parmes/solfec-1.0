/*
 * fem.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * finite element method
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

#include "bod.h"
#include "msh.h"

#ifndef __fem__
#define __fem__

/* formulation kind */
typedef enum {FEM_O1 = 1, FEM_O2} FEMFORM; /* must start with > 1 (see BODY_Pack) */

/* create FEM internals for a body (note that 'msh' might be NULL so that shp->data is a mesh) */
void FEM_Create (FEMFORM form, MESH *msh, SHAPE *shp, BULK_MATERIAL *mat, BODY *bod);

/* overwrite state */
void FEM_Overwrite_State (BODY *bod, double *q, double *u);

/* set initial rigid motion velocity */
void FEM_Initial_Velocity (BODY *bod, double *linear, double *angular);

/* initialise dynamic time stepping */
void FEM_Dynamic_Init (BODY *bod, SCHEME scheme);

/* estimate critical step for the dynamic scheme */
double FEM_Dynamic_Critical_Step (BODY *bod);

/* perform the initial half-step of the dynamic scheme */
void FEM_Dynamic_Step_Begin (BODY *bod, double time, double step);

/* perform the final half-step of the dynamic scheme */
void FEM_Dynamic_Step_End (BODY *bod, double time, double step);

/* initialise static time stepping */
void FEM_Static_Init (BODY *bod);

/* perform the initial half-step of the static scheme */
void FEM_Static_Step_Begin (BODY *bod, double time, double step);

/* perform the final half-step of the static scheme */
void FEM_Static_Step_End (BODY *bod, double time, double step);

/* motion x = x (X, state) */
void FEM_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x);

/* inverse motion X = X (x, state) */
void FEM_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X);

/* obtain spatial velocity at (gobj, referential point), expressed in the local spatial 'base' */
void FEM_Local_Velo (BODY *bod, VELOTIME time, SHAPE *shp, void *gobj, double *X, double *base, double *velo);

/* return transformation operator from the generalised to the local velocity space at (element, ref. point, base) */
MX* FEM_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *X, double *base);

/* compute current kinetic energy */
double FEM_Kinetic_Energy (BODY *bod);

/* get some values at a node of a geometrical object */
void FEM_Nodal_Values (BODY *bod, SHAPE *shp, void *gobj, int node, VALUE_KIND kind, double *values);

/* get some values at a referential point */
void FEM_Point_Values (BODY *bod, double *X, VALUE_KIND kind, double *values);

/* release FEM memory */
void FEM_Destroy (BODY *bod);

#endif
