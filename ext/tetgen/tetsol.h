/*
 * tetsol.h
 * Copyright (C) 2011, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * Tetgen to C interface
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

#ifdef __cplusplus
extern "C"
{
#endif
#include "../../msh.h"
#ifdef __cplusplus
}
#endif

#ifndef __tetsol__
#define __tetsol__
#ifdef __cplusplus
extern "C"
{
#endif

/* generate tetrahedrons based on an input mesh object; pass -INT_MAX for (vol/surf)ids to inherit from the mesh */
MESH* tetrahedralize1 (MESH *shape, double volume, double quality, int volid, int surfid, double min_angle, double max_angle, double ref_length);

/* generate tetrahedrons based on an input file; pass -INT_MAX for (vol/surf)ids to inherit from the input */
MESH* tetrahedralize2 (char *path, double volume, double quality, int volid, int surfid);

/* generate tetrahedrons bounded by triangular surfaces; TRI->flg store surfids */
MESH* tetrahedralize3 (TRI *tri, int m, double volume, double quality, int volid);

#ifdef __cplusplus
}
#endif
#endif
