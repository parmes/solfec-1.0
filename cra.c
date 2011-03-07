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
#include "err.h"

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
  /* TODO */
}
