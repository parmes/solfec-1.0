/*
 * sps.c
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

#include <string.h>
#include "alg.h"
#include "dom.h"
#include "sps.h"
#include "err.h"
#include "pck.h"

/* pools memory chunk size */
#define SPCHUNK 64

/* state indices */
#define COHESION 0

/* set item comparison */
static int setcmp (SURFACE_MATERIAL *a, SURFACE_MATERIAL *b)
{
  if (a->surf1 < b->surf1) return -1;
  else if (a->surf1 == b->surf1 && a->surf2 < b->surf2) return -1;
  else if (a->surf1 == b->surf1 && a->surf2 == b->surf2) return 0;
  else return 1;
}

/* create new label */
static char* newlabel (int size, char *label)
{
  int l = label ? strlen (label) : 0;
  char *out;

  if (l)
  {
    ERRMEM (out = malloc (l + 1));
    strcpy (out, label);
  }
  else
  {
    ERRMEM (out = malloc (256));
    snprintf (out, 256, "SURFACE_MATERIAL_%d", size);
    ERRMEM (out = realloc (out, strlen (out) + 1));
  }

  return out;
}

/* get state size */
int state_size (SURFACE_MATERIAL *mat)
{
  if (mat->cohesion > 0.0) return 1;
  else return 0;
}

/* get state flags */
short state_flags (SURFACE_MATERIAL_STATE *mat)
{
  short flags = 0;

  if (mat->base->cohesion > 0.0 && mat->state [COHESION] > 0.0) flags |= CON_COHESIVE;

  return flags;
}

/* ------------- interface -------------- */

/* create surface pair set */
SPSET* SPSET_Create ()
{
  SPSET *set;

  ERRMEM (set = malloc (sizeof (SPSET)));
  MEM_Init (&set->matmem, sizeof (SURFACE_MATERIAL), SPCHUNK);
  MEM_Init (&set->setmem, sizeof (SET), SPCHUNK);
  MEM_Init (&set->mapmem, sizeof (MAP), SPCHUNK);
  set->tab = NULL;
  set->set = NULL;
  set->map = NULL;
  set->size = 0;

  SURFACE_MATERIAL data = {0, INT_MAX, INT_MAX, NULL, SIGNORINI_COULOMB, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  SPSET_Insert (set, INT_MAX, INT_MAX, "DEFAULT", data);

  return set;
}

/* insert new material; surf1, surf2 can be set to INT_MAX to indicate default pairing;
 * e.g. (INT_MAX, INT_MAX) is same as SPSET_Default; (INT_MAX, 1) is all surfaces with surface 1, etc. */
SURFACE_MATERIAL* SPSET_Insert (SPSET *set, int surf1, int surf2, char *label, SURFACE_MATERIAL data)
{
  data.surf1 = MIN (surf1, surf2); /* ascending sort */
  data.surf2 = MAX (surf1, surf2); /* ascending sort */
 
  SURFACE_MATERIAL *out = SET_Find (set->set, &data, (SET_Compare) setcmp);

  if (out == NULL)
  {
    ERRMEM (out = MEM_Alloc (&set->matmem)); 

    out->surf1 = data.surf1;
    out->surf2 = data.surf2;

    out->label = newlabel (set->size, label);

    SET_Insert (&set->setmem, &set->set, out, (SET_Compare) setcmp);
    MAP_Insert (&set->mapmem, &set->map, out->label, out, (SET_Compare) strcmp);
  }
  else if (label && strcmp (out->label, label) != 0)
  {
    MAP_Delete (&set->mapmem, &set->map, out->label, (SET_Compare) strcmp);

    free (out->label);

    out->label = newlabel (set->size, label);

    MAP_Insert (&set->mapmem, &set->map, out->label, out, (SET_Compare) strcmp);
  }

  /* copy all data from the 'model' onwards */
  memcpy (&out->model, &data.model, sizeof (SURFACE_MATERIAL) - ((char*)&data.model - (char*)&data.index));

  ERRMEM (set->tab = realloc (set->tab, (set->size + 1) * sizeof (SURFACE_MATERIAL*)));
  set->tab [set->size] = out;
  out->index = set->size ++;


  return out;
}

/* find by surface pair */
SURFACE_MATERIAL* SPSET_Find (SPSET *set, int surf1, int surf2)
{
  SURFACE_MATERIAL *out, tmp;
  
  tmp.surf1 = MIN (surf1, surf2);
  tmp.surf2 = MAX (surf1, surf2);

  out = SET_Find (set->set, &tmp, (SET_Compare) setcmp);

  if (!out)
  {
    tmp.surf1 = surf1;
    tmp.surf2 = INT_MAX;
    out = SET_Find (set->set, &tmp, (SET_Compare) setcmp);
  }

  if (!out)
  {
    tmp.surf1 = surf2;
    tmp.surf2 = INT_MAX;
    out = SET_Find (set->set, &tmp, (SET_Compare) setcmp);
  }

  if (!out)
  {
    tmp.surf1 = tmp.surf2 = INT_MAX;
    out = SET_Find (set->set, &tmp, (SET_Compare) setcmp);
    ASSERT_DEBUG (out, "Default surface material not found!");
  }

  return out;
}

/* find by label */
SURFACE_MATERIAL* SPSET_Find_Label (SPSET *set, char *label)
{
  return MAP_Find (set->map, label, (MAP_Compare) strcmp);
}

/* is there a pair like this? */
int SPSET_Has (SPSET *set, int surf1, int surf2)
{
  SURFACE_MATERIAL tmp;
  
  tmp.surf1 = MIN (surf1, surf2);
  tmp.surf2 = MAX (surf1, surf2);

  if (SET_Contains (set->set, &tmp, (SET_Compare) setcmp)) return 1;
  return 0;
}

/* release memory */
void SPSET_Destroy (SPSET *set)
{
  SURFACE_MATERIAL *mat;
  SET *item;

  for (item = SET_First (set->set); item; item = SET_Next (item))
  {
    mat = item->data;
    free (mat->label);
  }

  MEM_Release (&set->matmem);
  MEM_Release (&set->setmem);
  MEM_Release (&set->mapmem);
  free (set->tab);
  free (set);
}

/* export MBFCP definition */
void SPSET_2_MBFCP (SPSET *set, FILE *out)
{
  SURFACE_MATERIAL **tab = set->tab;

  fprintf (out, "SURFACE_MATERIALS:\t%d\n", set->size);

  for (int i = 0; i < set->size + 1; i ++)
  {
    if (i)
    {
      fprintf (out, "SURF1:\t%d\n", tab [i]->surf1);
      fprintf (out, "SURF2:\t%d\n", tab [i]->surf2);
    }
    else
    {
      fprintf (out, "SURF1:\t%s\n", "ANY");
      fprintf (out, "SURF2:\t%s\n", "ANY");
    }

    switch (tab [i]->model)
    {
    case SIGNORINI_COULOMB:
      fprintf (out, "MODEL:\t%s\n", "SIGNORINI_COULOMB");
      fprintf (out, "FRICTION:\t%g\n", tab [i]->friction);
      fprintf (out, "COHESION:\t%g\n", tab [i]->cohesion);
      fprintf (out, "RESTITUTION:\t%g\n", tab [i]->restitution);
      break;
    case SPRING_DASHPOT:
      fprintf (out, "MODEL:\t%s\n", "SPRING_DASHPOT");
      fprintf (out, "FRICTION:\t%g\n", tab [i]->friction);
      fprintf (out, "COHESION:\t%g\n", tab [i]->cohesion);
      fprintf (out, "SPRING:\t%g\n", tab [i]->spring);
      fprintf (out, "DASHPOT:\t%g\n", tab [i]->dashpot);
      break;
    }
  }

  fprintf (out, "\n");
  free (tab);
}

/* transfer data from the source 'src' to a destiny 'dst' at the given 'time' */
short SURFACE_MATERIAL_Transfer (double time, SURFACE_MATERIAL *src, SURFACE_MATERIAL_STATE *dst)
{
  dst->base = src;

  if (dst->state) free (dst->state), dst->state = NULL;

  int size = state_size (src);

  if (size)
  {
    ERRMEM (dst->state = MEM_CALLOC (sizeof (double [size])));
    
    if (time == 0.0 && src->cohesion > 0.0) /* contacts created after the initial moment cannot get cohesive */
    {
      dst->state [COHESION] = src->cohesion;
    }
  }

  return state_flags (dst);
}

/* write surface material state */
void SURFACE_MATERIAL_Write_State (SURFACE_MATERIAL_STATE *mat, PBF *bf)
{
  int size = state_size (mat->base);

  PBF_Short (bf, &mat->base->index, 1);
  PBF_Int (bf, &size, 1);

  if (size) PBF_Double (bf, mat->state, size);
}

/* read surface material state; return contact state flags */
short SURFACE_MATERIAL_Read_State (SPSET *set, SURFACE_MATERIAL_STATE *mat, PBF *bf)
{
  short index;
  int size;

  if (mat->state) free (mat->state), mat->state = NULL;

  PBF_Short (bf, &index, 1);
  PBF_Int (bf, &size, 1);
  ASSERT_DEBUG (index <= set->size, "Invalid material index");
  mat->base = set->tab [index];
  if (size)
  {
    ERRMEM (mat->state = malloc (sizeof (double [size])));
    PBF_Double (bf, mat->state, size);
  }

  return state_flags (mat);
}

/* pack surface material state */
void SURFACE_MATERIAL_Pack_State (SURFACE_MATERIAL_STATE *mat, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  int size = state_size (mat->base);

  pack_int (isize, i, ints, mat->base->index);
  pack_int (isize, i, ints, size);

  if (size) pack_doubles (dsize, d, doubles, mat->state, size);
}

/* unpack surface material state; return contact state flags */
short SURFACE_MATERIAL_Unpack_State (SPSET *set, SURFACE_MATERIAL_STATE *mat, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int index;
  int size;

  if (mat->state) free (mat->state), mat->state = NULL;

  index = unpack_int (ipos, i, ints);
  size = unpack_int (ipos, i, ints);
  ASSERT_DEBUG (index <= set->size, "Invalid material index");
  mat->base = set->tab [index];
  if (size)
  {
    ERRMEM (mat->state = malloc (sizeof (double [size])));
    unpack_doubles (dpos, d, doubles, mat->state, size);
  }

  return state_flags (mat);
}

/* free material state memory */
void SURFACE_MATERIAL_Destroy_State (SURFACE_MATERIAL_STATE *mat)
{
  if (mat->state) free (mat->state);
}

/* cohesion state get */
double SURFACE_MATERIAL_Cohesion_Get (SURFACE_MATERIAL_STATE *mat)
{
  if (mat->state) return mat->state [COHESION];
  else return 0.0;
}

/* cohesion state set */
void SURFACE_MATERIAL_Cohesion_Set (SURFACE_MATERIAL_STATE *mat, double value)
{
  if (mat->state) mat->state [COHESION] = value;
}
