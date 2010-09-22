/*
 * sps.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * surface pair set (surface material)
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

#include "mem.h"
#include "map.h"
#include "set.h"
#include "pbf.h"

#ifndef __sps__
#define __sps__

typedef struct surfmat_state SURFACE_MATERIAL_STATE;
typedef struct surfmat SURFACE_MATERIAL;
typedef struct spset SPSET;

#define COHESION_EPSILON 1E-8

struct surfmat_state
{
  SURFACE_MATERIAL *base; /* base material */

  double *state; /* material state */
};

struct surfmat
{
  short index; /* tab index */

  int surf1,
      surf2; /* surface indices */

  char *label;

  enum
  { SIGNORINI_COULOMB,
    SPRING_DASHPOT } model;

  double friction,
         cohesion,
         restitution,
	 spring,
	 dashpot;
};

struct spset /* surface pair set */
{
  MEM matmem,
      setmem,
      mapmem;

  SURFACE_MATERIAL def; /* default data */

  SURFACE_MATERIAL **tab; /* table of surface materials */

  SET *set;   /* surface pair index based set */

  MAP *map;   /* label based map */

  int size;   /* number of materials */
};

/* The idea is to map surface pairs (identified by two integer numbers)
 * with data records defining material parameters. As the data record (SURFACE_MATERIAL)
 * might extend in the future it is more convenient to pass the whole record
 * by to the member routines. It is passed by value and if necessary - copied.
 * An obvious interface follows below */

/* create surface pair set */
SPSET* SPSET_Create ();

/* set up default material */
void SPSET_Default (SPSET *set, SURFACE_MATERIAL data);

/* insert new material */
SURFACE_MATERIAL* SPSET_Insert (SPSET *set, int surf1, int surf2, char *label, SURFACE_MATERIAL data);

/* find by surface pair */
SURFACE_MATERIAL* SPSET_Find (SPSET *set, int surf1, int surf2);

/* find by label */
SURFACE_MATERIAL* SPSET_Find_Label (SPSET *set, char *label);

/* is there a pair like this? */
int SPSET_Has (SPSET *set, int surf1, int surf2);

/* release memory */
void SPSET_Destroy (SPSET *set);

/* transfer data from the source 'src' to a destiny 'dst'
 * at the given time => e.g. cohesion will be >0 only for 'time' == 0 */
short SURFACE_MATERIAL_Transfer (double time, SURFACE_MATERIAL *src, SURFACE_MATERIAL_STATE *dst);

/* write surface material state */
void SURFACE_MATERIAL_Write_State (SURFACE_MATERIAL_STATE *mat, PBF *bf);

/* read surface material state; return contact state flags */
short SURFACE_MATERIAL_Read_State (SPSET *set, SURFACE_MATERIAL_STATE *mat, PBF *bf);

/* pack surface material state */
void SURFACE_MATERIAL_Pack_State (SURFACE_MATERIAL_STATE *mat, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack surface material state; return contact state flags */
short SURFACE_MATERIAL_Unpack_State (SPSET *set, SURFACE_MATERIAL_STATE *mat, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

/* free material state memory */
void SURFACE_MATERIAL_Destroy_State (SURFACE_MATERIAL_STATE *mat);

/* cohesion state get */
double SURFACE_MATERIAL_Cohesion_Get (SURFACE_MATERIAL_STATE *mat);

/* cohesion state set */
void SURFACE_MATERIAL_Cohesion_Set (SURFACE_MATERIAL_STATE *mat, double value);
#endif
