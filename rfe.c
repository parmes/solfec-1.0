/*
 * rfe.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * rough finite element method
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

#include "rfe.h"
#include "err.h"

/* create RFE internals for a body */
void RFE_Create (FEMFORM form, MESH *msh, SHAPE *shp, BULK_MATERIAL *mat, BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* get rough mesh */
MESH* RFE_Mesh (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
  return NULL;
}

/* overwrite state */
void RFE_Overwrite_State (BODY *bod, double *q, double *u)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* set initial rigid motion velocity */
void RFE_Initial_Velocity (BODY *bod, double *linear, double *angular)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* initialise dynamic time stepping */
void RFE_Dynamic_Init (BODY *bod, SCHEME scheme)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* estimate critical step for the dynamic scheme */
double RFE_Dynamic_Critical_Step (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
  return 0.0;
}

/* perform the initial half-step of the dynamic scheme */
void RFE_Dynamic_Step_Begin (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* perform the final half-step of the dynamic scheme */
void RFE_Dynamic_Step_End (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* initialise static time stepping */
void RFE_Static_Init (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* perform the initial half-step of the static scheme */
void RFE_Static_Step_Begin (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* perform the final half-step of the static scheme */
void RFE_Static_Step_End (BODY *bod, double time, double step)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* motion x = x (X, state) */
void RFE_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO: ele == NULL implies nodal update (X is within the mesh->ref_nodes) */
}

/* inverse motion X = X (x, state) */
void RFE_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* obtain velocity at (element, point), expressed in the local 'base' => all entities are spatial */
void RFE_Local_Velo (BODY *bod, VELOTIME time, SHAPE *shp, void *gobj, double *point, double *base, double *velo)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* return transformation operator from the generalised to the local velocity space at (element, point, base) */
MX* RFE_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
  return NULL;
}

/* compute current kinetic energy */
double RFE_Kinetic_Energy (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
  return 0.0;
}

/* get some values at a node of a geometrical object */
void RFE_Nodal_Values (BODY *bod, SHAPE *shp, void *gobj, int node, VALUE_KIND kind, double *values)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* get some values at a referential point */
void RFE_Point_Values (BODY *bod, double *point, VALUE_KIND kind, double *values)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}

/* release RFE memory */
void RFE_Destroy (BODY *bod)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED);
  /* TODO */
}
