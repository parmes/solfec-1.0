/*
 * cra.c
 * Copyright (C) 2011, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * body cracking
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

#include "dom.h"
#include "cra.h"
#include "mem.h"
#include "alg.h"
#include "err.h"

/* pseudo-rigid body cracking */
static void prb_crack (BODY *bod, BODY **one, BODY **two)
{
  double values [6], cauchy [9], pnt [3], vec [3], tension;
  CRACK *cra, *crb, *c;

  *one = *two = NULL;

  for (cra = bod->cra; cra; cra = cra->next)
  {
    BODY_Point_Values (bod, cra->point, VALUE_STRESS, values);
    cauchy [0] = values [0];
    cauchy [1] = values [3];
    cauchy [2] = values [4];
    cauchy [3] = cauchy [1];
    cauchy [4] = values [1];
    cauchy [5] = values [5];
    cauchy [6] = cauchy [2];
    cauchy [7] = cauchy [5];
    cauchy [8] = values [2];
    NVMUL (cauchy, cra->normal, vec);
    tension = DOT (cra->normal, vec);

    if (tension > cra->ft)
    {
      BODY_Cur_Vector (bod, NULL, cra->point, cra->normal, vec);
      BODY_Cur_Point (bod, NULL, NULL, cra->point, pnt);
      BODY_Split (bod, pnt, vec, cra->surfid, one, two);
      break;
    }
  }

  crb = cra;

  if (crb && *one && *two)
  {
    for (cra = bod->cra; cra; cra = cra->next)
    {
      if (cra == crb) continue;

      c = CRACK_Create();
      *c = *cra; /* XXX: maintain correctnes when CRACK stores more data */
      c->next = (*one)->cra;
      (*one)->cra = c;

      c = CRACK_Create();
      *c = *cra; /* XXX: maintain correctnes when CRACK stores more data */
      c->next = (*two)->cra;
      (*two)->cra = c;
    }
  }
}

/* finite element body cracking */
static void fem_crack (BODY *bod, BODY **one, BODY **two)
{
  /* TODO */
}

/* map contraints to new bodies */
static void map_constraints  (DOM *dom, BODY *bod, BODY *one, BODY *two)
{
  /* TODO */
}

/* create crack object */
CRACK* CRACK_Create ()
{
  return MEM_CALLOC (sizeof (CRACK));
}

/* delete crack object */
void CRACK_Destroy (CRACK *cra)
{
  free (cra);
}

/* delete list of cracks */
void CRACK_Destroy_List (CRACK *cra)
{
  CRACK *next;

  for (; cra; cra = next)
  {
    next = cra->next;
    CRACK_Destroy (cra);
  }
}

/* propagate cracks and adjust the domain */
void Propagate_Cracks (DOM *dom)
{
  BODY *bod, *one, *two, *next;

  for (bod = dom->bod; bod; bod = next)
  {
    next = bod->next;
    one = two = NULL;

    if (bod->cra)
    {
      switch (bod->kind)
      {
      case PRB: prb_crack (bod, &one, &two); break;
      case FEM: fem_crack (bod, &one, &two); break;
      default: break;
      }
    }

    if (one && two)
    {
      map_constraints (dom, bod, one, two);
      DOM_Remove_Body (dom, bod);
      DOM_Insert_Body (dom, one);
      DOM_Insert_Body (dom, two);
      BODY_Destroy (bod);
    }
    else
    {
      ASSERT_TEXT (one == NULL && two == NULL, "A body cracked, but only one part was created");
    }
  }
}
