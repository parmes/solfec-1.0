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

#include <limits.h>
#include "ext/tetgen/tetsol.h"
#include "sol.h"
#include "dom.h"
#include "dio.h"
#include "fem.h"
#include "fra.h"
#include "mem.h"
#include "alg.h"
#include "err.h"
#include "fem.h"
#include "pck.h"

typedef struct fracture_state FS;

struct fracture_state
{
  double area;
  double point [3];
  double R [3];

  FS *next;

  FS *inext;
};

/* free list */
static void fracture_state_free (FS *list)
{
  FS *next;

  for (;list; list = next)
  {
    next = list->next;
    free (list);
  }
}

/* write fracture state */
static void fracture_state_write (DOM *dom)
{
  char path [1024];
  SET *item;
  BODY *bod;
  CON *con;
  FILE *f;
  int n;
  XDR x;

#if MPI
  snprintf (path, 1024, "%s/fracture%d.dat", dom->solfec->outpath, dom->rank);
#else
  snprintf (path, 1024, "%s/fracture.dat", dom->solfec->outpath);
#endif
  ASSERT (f = fopen (path, "a"), ERR_FILE_OPEN);
  xdrstdio_create (&x, f, XDR_ENCODE);

  for (bod = dom->bod; bod; bod = bod->next)
  {
    if (bod->fracture)
    {
      xdr_u_int (&x, &bod->id);
      n = SET_Size (bod->con);
      xdr_int (&x, &n);
      for (item = SET_First (bod->con); item; item = SET_Next (item))
      {
	con = item->data;
	xdr_double (&x, &con->area);
	xdr_double (&x, &con->point[0]);
	xdr_double (&x, &con->point[1]);
	xdr_double (&x, &con->point[2]);
	xdr_double (&x, &con->R[0]);
	xdr_double (&x, &con->R[1]);
	xdr_double (&x, &con->R[2]);
      }

      bod->fracture = 0;
    }
  }

  xdr_destroy (&x);
  fclose (f);
}

/* read fracture state */
static FS* fracture_state_read (BODY *bod)
{
  FS *out = NULL, *item, *instance;
  char path [1024];
  unsigned int id;
  int i, n;
  FILE *f;
  XDR x;

  snprintf (path, 1024, "%s/fracture.dat", bod->dom->solfec->outpath);
  f = fopen (path, "r");
  /* TODO: read MPI mode data in case f == NULL but fractureRANK.dat exit */
  if (f)
  {
    xdrstdio_create (&x, f, XDR_DECODE);

    while (! feof (f))
    {
      xdr_u_int (&x, &id);
      xdr_int (&x, &n);
      for (i = 0, instance = NULL; i < n; i ++)
      {
        ERRMEM (item = malloc (sizeof (FS)));

	xdr_double (&x, &item->area);
	xdr_double (&x, &item->point[0]);
	xdr_double (&x, &item->point[1]);
	xdr_double (&x, &item->point[2]);
	xdr_double (&x, &item->R[0]);
	xdr_double (&x, &item->R[1]);
	xdr_double (&x, &item->R[2]);

	if (id == bod->id)
	{
	  item->inext = instance;
	  instance = item;

	  if (i == n-1)
	  {
	    item->next = out;
	    out = item;
	  }
	}
	else free (item);
      }
    }

    xdr_destroy (&x);
    fclose (f);
  }

  return out;
}


/* check fracture criterion */
void Fracture_Check (DOM *dom)
{
  short on = 0;
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
	  bod->fracture = 1;
	  on = 1;
	}

	if (bulk) ele = ele->next;
	else if (ele->next) ele = ele->next;
	else ele = msh->bulkeles, bulk = 1;
      }
    }
  }

  if (on) fracture_state_write (dom);
}

/* export data for fracture analysis in Yaffems (return number of exported analysis instances) */
int Fracture_Export_Yaffems (BODY *bod, double volume, double quality, FILE *output)
{
  SOLFEC *sol = bod->dom->solfec;
  FS *list, *it, *jt;
  MESH *msh;
  int n, m;

  if (!(bod->flags & BODY_CHECK_FRACTURE) || sol->mode == SOLFEC_WRITE) return 0;

  list = fracture_state_read (bod);

  if (list)
  {
    msh = tetrahedralize1 (bod->shape->data, volume, quality, -INT_MAX, -INT_MAX);

    MESH_2_MBFCP (msh, output);

    for (it = list, m = 0; it; it = it->next) m ++;

    fprintf (output, "INSTANCES:\t%d\n", m);

    for (it = list; it; it= it->next)
    {
      for (jt = it, n = 0; jt; jt = jt->inext) n ++;

      fprintf (output, "FORCES: %d\n", n);

      for (jt = it; jt; jt = jt->inext)
      {
	fprintf (output, "AREA %g POINT %g %g %g FORCE %g %g %g\n", jt->area,
	    jt->point[0], jt->point[1], jt->point[2], jt->R[0], jt->R[1], jt->R[2]);
      }
    }

    fracture_state_free (list);

    return m;
  }

  return 0;
}
