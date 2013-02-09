/*
 * psc.c
 * Copyright (C) 2013, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * parallel self-consistency tests
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

/* write body data to file before sending the body via MPI calles;
 * the file is SOLFEC->outpath/bodyID.data; */

#include <stdio.h>
#include "sol.h"
#include "dom.h"
#include "bod.h"
#include "psc.h"
#include "err.h"

static void psc_write_shape (SHAPE *shp)
{
  /* TODO */
}

static void psc_write_matrix (MX *mx)
{
  /* TODO */
}

/* write body data to file before sending the body via MPI calles;
 * the file is SOLFEC->outpath/bodyID.data; */
void PSC_Write_Body (BODY *bod)
{
  char path [1024];
  FILE *f;

  snprintf (path, 1024, "%s/body%d.data", bod->dom->solfec->outpath, bod->id);
  ASSERT (f = fopen (path, "w"), ERR_FILE_OPEN);

  fwrite (&bod->kind, sizeof (bod->kind), 1, f);

  fwrite (bod->mat->label, 1, strlen (bod->mat->label), f);

  fwrite (&bod->ref_mass, sizeof (double), 1, f);
  fwrite (&bod->ref_volume, sizeof (double), 1, f);
  fwrite (bod->ref_center, sizeof (double), 3, f);
  fwrite (bod->ref_tensor, sizeof (double), 9, f);

  fwrite (&bod->dofs, sizeof (int), 1, f);

  fwrite (&bod->form, sizeof (bod->form), 1, f);

  int confsize = bod->kind != FEM ? 12 : bod->form == REDUCED_ORDER ? bod->dofs + 9 : bod->dofs;

  fwrite (bod->conf, sizeof (double), confsize, f);
  fwrite (bod->velo, sizeof (double), bod->dofs, f);

  /* TODO => bod->forces */

  /* TODO => bod->cra */

  psc_write_shape (bod->shape);

  fwrite (bod->extents, sizeof (double), 6, f);

  fwrite (&bod->scheme, sizeof (bod->scheme), 1, f);

  psc_write_matrix (bod->inverse);

  psc_write_matrix (bod->M);

  psc_write_matrix (bod->K);

  fwrite (&bod->damping, sizeof (double), 1, f);

  if (bod->evec)
  {
    psc_write_matrix (bod->evec);

    fwrite (bod->eval, sizeof (double), bod->evec->m, f);
  }

  /* XXX: skip bod->label as non-essential */

  fwrite (&bod->form, sizeof (bod->form), 1, f);

  /* XXX: skip bod->mesh for the moment */

  /* XXX: skip bod->energy as non-essential */

  /* XXX: skip bod->fracture as non-essential */

  fclose (f);
}

/* after the body has been imported via MPI calls,
 * read body data from file and compare */
void PSC_Test_Body (BODY *bod)
{
}
