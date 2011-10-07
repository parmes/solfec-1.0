/*
 * mat.c
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

#include <string.h>
#include "mat.h"
#include "but.h"
#include "svk.h"
#include "alg.h"
#include "err.h"

#define CHUNK 128

/* create new label */
static char* newlabel (int size, char *label)
{
  char *out;
  int l;

  if ((l = label ? strlen (label) : 0))
  {
    ERRMEM (out = malloc (l + 1));
    strcpy (out, label);
  }
  else
  {
    ERRMEM (out = malloc (256));
    snprintf (out, 256, "BULK_MATERIAL_%d", size);
    ERRMEM (out = realloc (out, strlen (out) + 1));
  }

  return out;
}

/* create bulk material set */
MATSET* MATSET_Create ()
{
  MATSET *set;

  ERRMEM (set = malloc (sizeof (MATSET)));
  MEM_Init (&set->matmem, sizeof (BULK_MATERIAL), CHUNK);
  MEM_Init (&set->mapmem, sizeof (MAP), CHUNK);
  set->map = NULL;
  set->size = 0;
  return set;
}

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
double BULK_MATERIAL_ROUTINE (BULK_MATERIAL *mat, double *F, double a, double *P, double *K)
{
  double J;

  switch (mat->model)
  {
  case KIRCHHOFF:
  {
    double mat_lambda, mat_mi;

    mat_lambda = lambda (mat->young, mat->poisson);
    mat_mi = mi (mat->young, mat->poisson);
    if (K) SVK_Tangent_C (mat_lambda, mat_mi, a, 9, F, K);
    if (P) J = SVK_Stress_C (mat_lambda, mat_mi, a, F, P);
    else J = DET (F);
  }
  break;
  case TSANG_MARSDEN:
  {
    J = DET (F);

    /* TODO */
  }
  break;
  }

  return J;
}

/* insert new material */
BULK_MATERIAL* MATSET_Insert (MATSET *set, char *label, BULK_MATERIAL data)
{
  BULK_MATERIAL *out;

  if (!label || !(out = MAP_Find (set->map, label, (MAP_Compare) strcmp)))
  {
    ERRMEM (out = MEM_Alloc (&set->matmem)); 
    out->label = newlabel (set->size, label);
    MAP_Insert (&set->mapmem, &set->map, out->label, out, (MAP_Compare) strcmp);
    set->size ++;
  }

  /* copy all data from the 'model' onwards */
  memcpy (&out->model, &data.model, sizeof (BULK_MATERIAL) - sizeof (char*));

  return out;
}

/* find by label */
BULK_MATERIAL* MATSET_Find (MATSET *set, char *label)
{
  return MAP_Find (set->map, label, (MAP_Compare) strcmp);
}

/* release memory */
void MATSET_Destroy (MATSET *set)
{
  BULK_MATERIAL *mat;
  MAP *item;

  for (item = MAP_First (set->map); item; item = MAP_Next (item))
  {
    mat = item->data;
    free (mat->label);
  }

  MEM_Release (&set->matmem);
  MEM_Release (&set->mapmem);
  free (set);
}

/* export MBFCP definition */
void MATSET_2_MBFCP (MATSET *set, FILE *out)
{
  BULK_MATERIAL *mat;
  MAP *item;

  fprintf (out, "BULK_MATERIALS:\t%d\n", set->size);

  for (item = MAP_First (set->map); item; item = MAP_Next (item))
  {
    mat = item->data;
    switch (mat->model)
    {
    case KIRCHHOFF:
    case TSANG_MARSDEN:
      fprintf (out, "LABEL:\t%s\n", mat->label);
      fprintf (out, "MODEL:\t%s\n", "KIRCHHOFF");
      fprintf (out, "YOUNG:\t%g\n", mat->young);
      fprintf (out, "POISSON:\t%g\n", mat->poisson);
      fprintf (out, "DENSITY:\t%g\n", mat->density);
      break;
    }
  }

  fprintf (out, "\n");
}
