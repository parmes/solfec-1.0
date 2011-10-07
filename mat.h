/*
 * mat.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * bulk material
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

#include <stdio.h>
#include "mem.h"
#include "map.h"
#include "fld.h"

#ifndef __mat__
#define __mat__

/* UMAT: user material subroutine mimicking the ABAQUS counterpart
 * ----------------------------------------------------------------
 * stress (ntens): Cauchy stress at the beginning (IN) or end (OUT) of the increment;
 * statev (nstatv): solution dependent state variables at the beginning (IN) or end (OUT) of the increment;
 * ddsdde (ntens, ntens): Jacobian matrix of the constitutive model: d (stress increment) / d (strain increment);
 * sse, spd, scd: specific elastic strain energy, plastic dissipation and creep dissipation, respectively;
 * rpl: mechanical volumetric heat generation per unit time at the end of the increment;
 * ddsddt (ntens): variation of stress increments with respect to the temperature;
 * drplde (ntens): variation of 'prl' with respect to the strain increments;
 * drpldt: variation of 'prl' with respect to the temperature;
 * stran (ntens): the total strains minus the thermal strains; rotated according to rigid motion approximations to logarithmic strain;
 * dstran (ntens): the total strain increments minus the thermal strain increments;
 * time (2): value of step time and total time, respectively,  at the beginning of the current increment;
 * dtime: time increment;
 * temp: temperature at the start of the increment;
 * dtemp: increment of temperature;
 * predef: array of interpolated values of predefined field variables at this point at the start of the increment, based on the values read in at the nodes;
 * dpred: array of increments of predefined field variables;
 * cmname: user-defined material name;
 * ndi: number of direct stress components at this point;
 * nshr: number of engineering shear stress components at this point;
 * ntens: size of the stress or strain component array (ndi + nshr);
 * nstatv: number of solution-dependent state variables that are associated with this material type;
 * props (nprops): user-specified array of material constants associated with this user material;
 * nprops: user-defined number of material constants associated with this user material;
 * coords: an array containing the current coordinates of this point;
 * drot (3, 3): rotation increment matrix; provided so that vector- or tensor-valued state variables can be rotated in this subroutine;
 * pnewdt: ratio of suggested new time increment to the time increment being used;
 * celent: characteristic element length; 
 * dfgrd0 (3, 3): array containing the deformation gradient at the beginning of the increment;
 * dfgrd1 (3, 3): array containing the deformation gradient at the end of the increment;
 * noel: element number;
 * npt: integration point number;
 * layer: layer number (for composite shells and layered solids);
 * kspt: section point number within the current layer;
 * kstep: step number;
 * kinc: increment number;
 */

typedef void (*UMAT) (double *stress, double *statev, double *ddsdde, double *sse, double *spd, double *scd,
                      double *rpl, double *ddsddt, double *drplde, double *drpldt, double *stran, double *dstran,
		      double *time, double *dtime, double *temp, double *dtemp, double *predef, double *dpred,
		      char *cmname, int *ndi, int *nshr, int *ntens, int *nstatv, double *props, int *nprops,
		      double *coords, double *drot, double *pnewdt, double *celent, double *dfgrd0, double *dfgrd1,
		      int *noel, int *npt, int *layer, int *kspt, int *kstep, int *kinc);

typedef struct bulkmat BULK_MATERIAL;
typedef struct matset MATSET;

struct bulkmat
{
  char *label;

  enum
  {
    KIRCHHOFF,      /* elastic material */
    TSANG_MARSDEN   /* nuclear graphite */
  } model;

  double young,
         poisson,
         density;

  UMAT umat;

  int nfield, /* number of fields (stored at mesh nodes) */
      nstate; /* number of state variables (stored at integration points) */

  FIELD **fld; /* fields */
};

struct matset
{
  MEM matmem,
      mapmem;

  MAP *map;   /* label based map */

  int size;   /* number of materials */
};

/* bulk material routine
 * ---------------------
 *  mat (IN) - material
 *  F (IN) - deformation gradient (column-wise)
 *  a (IN) - coefficient that will scale P and K
 *  P (OUT) - first Piola tensor (column-wise); NULL allowed
 *  K (OUT) - tangent dP/dF; NULL allowed
 * --------------------------------------
 *  return det (F)
 */
double BULK_MATERIAL_ROUTINE (BULK_MATERIAL *mat, double *F, double a, double *P, double *K);

/* create bulk material set */
MATSET* MATSET_Create ();

/* insert new material */
BULK_MATERIAL* MATSET_Insert (MATSET *set, char *label, BULK_MATERIAL data);

/* find by label */
BULK_MATERIAL* MATSET_Find (MATSET *set, char *label);

/* release memory */
void MATSET_Destroy (MATSET *set);

/* export MBFCP definition */
void MATSET_2_MBFCP (MATSET *set, FILE *out);

#endif
