/*
 * fld.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * scalar field (Python defined)
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

#include <Python.h>
#include <structmember.h>
#include <string.h>
#include "fld.h"
#include "lng.h"
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
    snprintf (out, 256, "FIELD_%d", size);
    ERRMEM (out = realloc (out, strlen (out) + 1));
  }

  return out;
}

/* evaluate field */
double FIELD_Value (FIELD *fld, double x, double y, double z, double t)
{
  PyObject *result, *args;
  double ret = 0.0;
  int n;

  if (fld->data)
  {
    if (PyTuple_Check (fld->data))
    {
      if (!(args = PyTuple_New (PyTuple_Size (fld->data) + 4))) goto err;

      for (n = 0; n < PyTuple_Size (fld->data); n ++)
	PyTuple_SetItem (args, n, PyTuple_GetItem (fld->data, n));

      PyTuple_SetItem (args, n, PyFloat_FromDouble (x));
      PyTuple_SetItem (args, n+1, PyFloat_FromDouble (y));
      PyTuple_SetItem (args, n+2, PyFloat_FromDouble (z));
      PyTuple_SetItem (args, n+3, PyFloat_FromDouble (t));
    }
    else args = Py_BuildValue ("(d, d, d, d)", x, y, z, t);
  }
  else args = Py_BuildValue ("(d, d, d, d)", x, y, z, t);

  result = PyObject_CallObject (fld->call, args); /* call user callback */

  Py_DECREF (args);

  if (result)
  {
    ret = PyFloat_AsDouble (result);

    Py_DECREF (result);
  }
  else /* error during the Python callback run */
  {
err:
    PyErr_Print (); /* print traceback */
#if MPI
    MPI_Abort (MPI_COMM_WORLD, 3000);
#else
    lngfinalize ();
#endif
    exit (1);
  }

  return ret;
}

/* create field set */
FISET* FISET_Create ()
{
  FISET *set;

  ERRMEM (set = malloc (sizeof (FISET)));
  MEM_Init (&set->fldmem, sizeof (FIELD), CHUNK);
  MEM_Init (&set->mapmem, sizeof (MAP), CHUNK);
  set->map = NULL;
  set->size = 0;
  return set;
}

/* insert new field */
FIELD* FISET_Insert (FISET *set, char *label, FIELD data)
{
  FIELD *out;

  if (!label || !(out = MAP_Find (set->map, label, (MAP_Compare) strcmp)))
  {
    ERRMEM (out = MEM_Alloc (&set->fldmem)); 
    out->label = newlabel (set->size, label);
    MAP_Insert (&set->mapmem, &set->map, out->label, out, (MAP_Compare) strcmp);
    set->size ++;
  }

  out->data = data.data;
  out->call = data.call;

  return out;
}

/* find by label */
FIELD* FISET_Find (FISET *set, char *label)
{
  return MAP_Find (set->map, label, (MAP_Compare) strcmp);
}

/* release memory */
void FISET_Destroy (FISET *set)
{
  FIELD *fld;
  MAP *item;

  for (item = MAP_First (set->map); item; item = MAP_Next (item))
  {
    fld = item->data;
    free (fld->label);
  }

  MEM_Release (&set->fldmem);
  MEM_Release (&set->mapmem);
  free (set);
}
