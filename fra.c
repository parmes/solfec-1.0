/*
 * fra.c
 * Copyright (C) 2013, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * fracture code
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
#include "fem.h"
#include "fra.h"
#include "mem.h"
#include "alg.h"
#include "err.h"
#include "fem.h"
#include "pck.h"

/* check fracture criterion */
void Fracture_Check (DOM *dom)
{
  BODY *bod;

  for (bod = dom->bod; bod; bod = bod->next)
  {
    if (bod->flags & BODY_CHECK_FRACTURE)
    {
      MESH *msh = FEM_MESH (bod);
      double energy, volume;
      ELEMENT *ele;
      int bulk;

      for (ele = msh->surfeles, bulk = 0; ele; )
      {
        BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);

	energy = FEM_Element_Internal_Energy (bod, msh, ele, &volume);

	if (energy > mat->criten * volume)
	{
	  FRACTURE_TIME *ft = MEM_Alloc (&dom->ftlmem);
          ft->time = dom->time;
	  ft->next = bod->ftl;
	  bod->ftl = ft;
	}

	if (bulk) ele = ele->next;
	else if (ele->next) ele = ele->next;
	else ele = msh->bulkeles, bulk = 1;
      }
    }
  }
}

/* export data for fracture analysis in Yaffems (return number of exported analysis instances) */
int Fracute_Export_Yaffems (BODY *bod, double voume, double quality, FILE *output)
{
  ASSERT (0, ERR_NOT_IMPLEMENTED); /* FIXME/TODO */

  return 0;
}
