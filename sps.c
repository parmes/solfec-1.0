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

#define SPCHUNK 128

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
    snprintf (out, 256, "SURFACE_MATERIAL_%d", size);
    ERRMEM (out = realloc (out, strlen (out) + 1));
  }

  return out;
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
  set->def.friction = 
  set->def.cohesion = 
  set->def.restitution = 0.0;
  set->set = NULL;
  set->map = NULL;
  set->size = 0;
  return set;
}

/* set up default material */
void SPSET_Default (SPSET *set, SURFACE_MATERIAL data)
{ 
  set->def = data;
  set->def.label = "DEFAULT";
}

/* insert new material */
SURFACE_MATERIAL* SPSET_Insert (SPSET *set, int surf1, int surf2, char *label, SURFACE_MATERIAL data)
{
  SURFACE_MATERIAL *out;

  data.surf1 = MIN (surf1, surf2); /* this is not expected from the user */
  data.surf2 = MAX (surf1, surf2); /* ... */
  
  if (! (out = SET_Find (set->set, &data, (SET_Compare) setcmp)))
  {
    ERRMEM (out = MEM_Alloc (&set->matmem)); 
    out->surf1 = data.surf1;
    out->surf2 = data.surf2;
    out->label = newlabel (set->size, label);
    SET_Insert (&set->setmem, &set->set, out, (SET_Compare) setcmp);
    MAP_Insert (&set->mapmem, &set->map, out->label, out, (SET_Compare) strcmp);
    set->size ++;
  }

  /* copy all data from the 'model' onwards */
  memcpy (&out->model, &data.model, sizeof (SURFACE_MATERIAL) - 2 * sizeof (int) - sizeof (char*));

  return out;
}

/* find by surface pair */
SURFACE_MATERIAL* SPSET_Find (SPSET *set, int surf1, int surf2)
{
  SURFACE_MATERIAL *out, tmp;
  
  tmp.surf1 = MIN (surf1, surf2);
  tmp.surf2 = MAX (surf1, surf2);

  if ((out = SET_Find (set->set, &tmp, (SET_Compare) setcmp))) return out;
  else return &set->def;
}

/* find by label */
SURFACE_MATERIAL* SPSET_Find_Label (SPSET *set, char *label)
{
  SURFACE_MATERIAL *mat;
  
  mat = MAP_Find (set->map, label, (MAP_Compare) strcmp);

  if (!mat && strcmp (set->def.label, label) == 0) return &set->def;
  else return mat;
}

/* is there a pair like this? */
int SPSET_Has (SPSET *set, int surf1, int surf2)
{
  SURFACE_MATERIAL tmp;
  
  tmp.surf1 = MIN (surf1, surf2);
  tmp.surf2 = MAX (surf1, surf2);

  if (SET_Find (set->set, &tmp, (SET_Compare) setcmp)) return 1;
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
  free (set);
}

/* transfer data from the source 'src' to a destiny 'dst' at the given 'time' */
short SURFACE_MATERIAL_Transfer (double time, SURFACE_MATERIAL *src, SURFACE_MATERIAL *dst)
{
  short state = 0;

  *dst = *src;

  if (time > 0.0)
  {
    dst->cohesion = 0.0; /* contacts created after the initial moment cannot get cohesive */
  }

  if (dst->cohesion > 0.0)
  {
    state |= CON_COHESIVE;
  }

  return state;
}
