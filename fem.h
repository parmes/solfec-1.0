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
#include "cvx.h"

#ifndef __fem__
#define __fem__

/* formulation */
typedef enum {TOTAL_LAGRANGIAN = 1, BODY_COROTATIONAL} FEMFORM; /* must be > 1 (see BODY_Pack in bod.c) */ 

/* create FEM internals for a body (note that 'msh' might be NULL so that shp->data is a mesh) */
void FEM_Create (FEMFORM form, MESH *msh, SHAPE *shp, BULK_MATERIAL *mat, BODY *bod);

/* overwrite state */
void FEM_Overwrite_State (BODY *bod, double *q, double *u);

/* set initial rigid motion velocity */
void FEM_Initial_Velocity (BODY *bod, double *linear, double *angular);

/* initialise dynamic time stepping */
void FEM_Dynamic_Init (BODY *bod);

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

/* motion x = x (element, local point) */
void FEM_Cur_Point_Ext (BODY *bod, ELEMENT *ele, double *X, double *point, double *x);

/* inverse motion X = X (x, state) */
void FEM_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X);

/* pull-forward v = {dx/dX} V (X, state) */
void FEM_Cur_Vector (BODY *bod, ELEMENT *ele, double *X, double *V, double *v);

/* push-back V = {dX/dx} v (x, state) */
void FEM_Ref_Vector (BODY *bod, ELEMENT *ele, double *x, double *v, double *V);

/* obtain spatial velocity at (gobj, referential point), expressed in the local spatial 'base' */
void FEM_Local_Velo (BODY *bod, SHAPE *shp, void *gobj, double *X, double *base, double *prevel, double *curvel);

/* return transformation operator from the generalised to the local velocity space at (element, ref. point, base) */
MX* FEM_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *X, double *base);

/* compute current kinetic energy */
double FEM_Kinetic_Energy (BODY *bod);

/* get some values at a local point of an element */
void FEM_Element_Point_Values (BODY *bod, ELEMENT *ele, double *point, VALUE_KIND kind, double *values);

/* get some values at a referential point (ele cam be NULL if not known) */
void FEM_Point_Values (BODY *bod, ELEMENT *ele, double *X, VALUE_KIND kind, double *values);

/* get some values at a curent mesh node (node points inside MESH->cur_nodes) */
void FEM_Cur_Node_Values (BODY *bod, double *node, VALUE_KIND kind, double *values);

/* issued by state reading routines of body interface */
void FEM_Update_Rough_Mesh (BODY *bod);

/* split body by referential plane; output two bodies with inherited state of the input body */
void FEM_Split (BODY *bod, double *point, double *normal, int surfid, BODY **one, BODY **two);

/* release FEM memory */
void FEM_Destroy (BODY *bod);

#if MPI
/* get configuration packing size */
int FEM_Conf_Pack_Size (BODY *bod);

/* get velocity packing size */
int FEM_Velo_Pack_Size (BODY *bod);
#endif

/* compute c = alpha * INVERSE (bod) * b + beta * c */
void FEM_Invvec (double alpha, BODY *bod, double *b, double beta, double *c);

/* create approximate inverse operator */
MX* FEM_Approx_Inverse (BODY *bod);

#endif
