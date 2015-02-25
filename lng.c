/*
 * lng.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * Python based input language parser
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

#if POSIX
#include <sys/stat.h>
#endif
#include <Python.h>
#include <structmember.h>
#include <limits.h>
#include <float.h>
#include "ext/tetgen/tetsol.h"
#include "solfec.h"
#include "alg.h"
#include "sol.h"
#include "rnd.h"
#include "lng.h"
#include "fem.h"
#include "box.h"
#include "goc.h"
#include "err.h"
#include "eli.h"
#include "fra.h"

#if MPI
#include <mpi.h>
#endif

#ifndef Py_RETURN_FALSE
#define Py_RETURN_FALSE return Py_INCREF(Py_False), Py_False
#endif

#ifndef Py_RETURN_TRUE
#define Py_RETURN_TRUE return Py_INCREF(Py_True), Py_True
#endif

#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE return Py_INCREF(Py_None), Py_None
#endif

typedef struct callback_pair CALLBACK_PAIR;

struct callback_pair
{
  int id;
  PyObject *data;
  PyObject *call;
  CALLBACK_PAIR *next;
};

static int callback_id = 1;

static CALLBACK_PAIR *callback_pairs = NULL; /* global stack of callback pairs */

/* string buffer length */
#define BUFLEN 512

/* minimal type initialization */
#define TYPEINIT(typedesc, type, name, flags, dealloc, new, methods, members, getset)\
memset (&(typedesc), 0, sizeof (PyTypeObject));\
(typedesc).tp_basicsize = sizeof (type);\
(typedesc).tp_name = name;\
(typedesc).tp_flags = flags;\
(typedesc).tp_dealloc = (destructor)dealloc;\
(typedesc).tp_new = new;\
(typedesc).tp_methods = methods;\
(typedesc).tp_members = members;\
(typedesc).tp_getset = getset

/*
 * callback stack
 */

/* push callback pair data on global stack */
static void callback_pair_push (PyObject *data, PyObject *call)
{
  CALLBACK_PAIR *pair;

  ERRMEM (pair = malloc (sizeof (CALLBACK_PAIR)));

  Py_XINCREF (data);
  Py_XINCREF (call);

  pair->id = callback_id ++;
  pair->data = data;
  pair->call = call;

  pair->next = callback_pairs;
  callback_pairs = pair;
}

/* pop all data and decrement reference counts */
static void callback_pairs_destroy ()
{
  CALLBACK_PAIR *pair, *next;

  for (pair = callback_pairs; pair; pair = next)
  {
    next = pair->next;

    Py_XDECREF (pair->data);
    Py_XDECREF (pair->call);

    free (pair);
  }
}

/*
 * utilities
 */

/* test whether an object is a string */
static int is_string (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyString_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a string object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* return NULL-termined string (obj was verified with is_string) */
static char* as_string (PyObject *obj)
{
  if (obj) return PyString_AsString (obj);

  return NULL;
}

/* a positive number test */
static int is_positive (double num, char *var)
{
  if (num <= 0)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be positive", var);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

/* a non-negative number test */
static int is_non_negative (double num, char *var)
{
  if (num < 0)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be non-negative", var);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

/* a in range number test */
static int is_in_range (double num, char *var, double min, double max)
{
  if (num < min || num > max)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be in range [%.2f, %.2f]", var, min, max);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

#if !MPI
/* a bigger test */
static int is_gt (double num, char *var, double val)
{
  if (num <= val)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be > %g", var, val);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}
#endif

/* test whether a number <= val */
static int is_le (double num, char *var, double val)
{
  char buf [BUFLEN];

  if (num > val)
  {
    sprintf (buf, "'%s' must be a number <= %g", var, val);
    PyErr_SetString (PyExc_TypeError, buf);
    return 0;
  }

  return 1;
}

/* test whether a number is >= val */
static int is_ge (double num, char *var, double val)
{
  char buf [BUFLEN];

  if (num < val)
  {
    sprintf (buf, "'%s' must be >= %g", var, val);
    PyErr_SetString (PyExc_TypeError, buf);
    return 0;
  }

  return 1;
}

/* test whether a num is in (lo, hi] */
static int is_gt_le (double num, char *var, double lo, double hi)
{
  char buf [BUFLEN];

  if (num <= lo || num > hi)
  {
    sprintf (buf, "'%s' must be a number > %g and <= %g", var, lo, hi);
    PyErr_SetString (PyExc_TypeError, buf);
    return 0;
  }

  return 1;
}

/* test whether a num is in [lo, hi] */
static int is_ge_le (double num, char *var, double lo, double hi)
{
  char buf [BUFLEN];

  if (num < lo || num > hi)
  {
    sprintf (buf, "'%s' must be a number >= %g and <= %g", var, lo, hi);
    PyErr_SetString (PyExc_TypeError, buf);
    return 0;
  }

  return 1;
}

/* test whether an object is a list and if so check
 * if list length is divisible by div and >= len */
static int is_list (PyObject *obj, char *var, int div, int len)
{
  if (obj)
  {
    if (!PyList_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (div > 0 && len > 0 &&
	!(PyList_Size (obj) % div == 0 && PyList_Size (obj) >= len))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must have N * %d elements, where N >= %d", var, div, len / div);
      PyErr_SetString (PyExc_ValueError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether an object is a tuple of length len */
static int is_tuple (PyObject *obj, char *var, int len)
{
  if (obj)
  {
    if (!PyTuple_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a tuple object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (len > 0 && PyTuple_Size (obj) != len)
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must have %d elements", var, len);
      PyErr_SetString (PyExc_ValueError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether an object is a number */
static int is_number (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyNumber_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a number", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether an object is a number > val */
static int is_number_gt (PyObject *obj, char *var, double val)
{
  if (obj)
  {
    char buf [BUFLEN];

    if (!PyNumber_Check (obj))
    {
      sprintf (buf, "'%s' must be a number", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
    else
    {
      double num = PyFloat_AsDouble (obj);

      if (num <= val)
      {
	sprintf (buf, "'%s' must be a number > %g", var, val);
	PyErr_SetString (PyExc_TypeError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* test whether an object is a number >= val */
static int is_number_ge (PyObject *obj, char *var, double val)
{
  if (obj)
  {
    char buf [BUFLEN];

    if (!PyNumber_Check (obj))
    {
      sprintf (buf, "'%s' must be a number", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
    else
    {
      double num = PyFloat_AsDouble (obj);

      if (num < val)
      {
	sprintf (buf, "'%s' must be a number >= %g", var, val);
	PyErr_SetString (PyExc_TypeError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* test whether an object is a number >= val0 and < val1 */
static int is_number_gt_le (PyObject *obj, char *var, double lo, double hi)
{
  if (obj)
  {
    char buf [BUFLEN];

    if (!PyNumber_Check (obj))
    {
      sprintf (buf, "'%s' must be a number", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
    else
    {
      double num = PyFloat_AsDouble (obj);

      if (num <= lo || num > hi)
      {
	sprintf (buf, "'%s' must be a number > %g and <= %g", var, lo, hi);
	PyErr_SetString (PyExc_TypeError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* test whether an object is a list (details as above) or a number */
static int is_list_or_number (PyObject *obj, char *var, int div, int len)
{
  if (obj)
  {
    if (!(PyList_Check (obj) || PyNumber_Check (obj)))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list or a number object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (PyList_Check (obj))
    {
      if (!(PyList_Size (obj) % div == 0 && PyList_Size (obj) >= len))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must have N * %d elements, where N >= %d", var, div, len / div);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* test whether an object is a list (details as above) or a string */
static int is_list_or_string (PyObject *obj, char *var, int div, int len)
{
  if (obj)
  {
    if (!(PyString_Check (obj) || PyList_Check (obj)))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list or a string object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (PyList_Check (obj))
    {
      if (!(PyList_Size (obj) % div == 0 && PyList_Size (obj) >= len))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must have N * %d elements, where N >= %d", var, div, len / div);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* test whether obj is a Python function pointer */
static int is_callable (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyCallable_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be callable", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* define keywords */
#define KEYWORDS(...) char *kwl [] = {__VA_ARGS__, NULL}

/* parse arguments with keywords */
#define PARSEKEYS(fmt, ...) if (!PyArg_ParseTupleAndKeywords (args, kwds, fmt, kwl, __VA_ARGS__)) return NULL

/* parse arguments without keywords */
#define PARSE(fmt, ...) if (!PyArg_ParseTuple (args, fmt, __VA_ARGS__)) return NULL

/* object types assertion */
#define TYPETEST(test) if(!(test)) return NULL

/* string argument if block comparison */
#define IFIS(obj, val) if (strcmp (PyString_AsString (obj), val) == 0)
#define ELIF(obj, val) else if (strcmp (PyString_AsString (obj), val) == 0)
#define ELSE else

/*
 * CONVEX => object
 */

typedef struct lng_CONVEX lng_CONVEX;

static PyTypeObject lng_CONVEX_TYPE; /* type descriptor */

struct lng_CONVEX
{
  PyObject_HEAD
  CONVEX *cvx;
};

/* test whether an object is of CONVEX type */
static int is_convex_test (PyObject *obj)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_CONVEX_TYPE)) return 0;
  }

  return 1;
}

/* test whether an object is of CONVEX type */
static int is_convex (lng_CONVEX *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_CONVEX_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a CONVEX object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* constructor */
static PyObject* lng_CONVEX_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  int volid, nver, *fac, nfac, *surfaces, i, j, k, l, n, m, o;
  KEYWORDS ("vertices", "faces", "volid", "convex");
  PyObject *vertices, *faces;
  lng_CONVEX *self, *convex;
  double *ver;
  CONVEX *cvx; /* convex->cvx */

  self = (lng_CONVEX*)type->tp_alloc (type, 0);

  if (self)
  {
    convex = NULL;
    cvx = NULL;

    PARSEKEYS ("OOi|O", &vertices, &faces, &volid, &convex);

    TYPETEST (is_list (vertices, kwl[0], 3, 12) && is_list (faces, kwl[1], 1, 0) && is_convex (convex, kwl[3]));

    /* test face definitions */
    l = PyList_Size (faces);
    for (nfac = i = 0; i < l; i ++)
    {
      k = PyInt_AsLong (PyList_GetItem (faces, i));

      if (k < 3) /* face size */
      {
	PyErr_SetString (PyExc_ValueError, "A face must have more than 2 vertices");
	return NULL;
      }

      /* add one more for the
       * surface definition */
      i += (k + 1);

      if (i >= l) /* incomplete */
      {
	PyErr_SetString (PyExc_ValueError, "The last face definition is incomplete");
	return NULL;
      }

      nfac ++; /* only when not incomplete */
    }

    if (nfac < 4) /* finished before defining 4 faces */
    {
      PyErr_SetString (PyExc_ValueError, "At least 4 faces must be defined");
      return NULL;
    }

    /* try to allocate and read some data now */
    ERRMEM (fac = malloc (sizeof (int) * (l - nfac)));
    ERRMEM (surfaces = malloc (sizeof (int) * nfac));
    nver = PyList_Size (vertices) / 3; /* vertices count */

    for (i = n = m = 0; i < l;)
    {
      fac [n ++] = k = PyInt_AsLong (PyList_GetItem (faces, i ++));

      for (j = 0; j < k; j ++)
      {
        fac [n ++] = PyInt_AsLong (PyList_GetItem (faces, i ++));
	if (fac [n-1] < 0 || fac [n-1] >= nver) /* must be within the right range */
	{
	  char buf [BUFLEN];
	  sprintf (buf, "Vertex %d in face %d is outside of range [0, %d]",j , m, nver-1);
	  PyErr_SetString (PyExc_ValueError, buf);
	  return NULL;
	}
      }

      /* test for repeated vertices in face definition */
      for (j = 1; j <= k; j ++)
      {
	for (o = j + 1; o <= k; o ++)
	{
	  if (fac [n-j] == fac [n-o])
	  {
	    char buf [BUFLEN];
	    sprintf (buf, "Vertices %d and %d in face %d are the same", k-j, k-o, m);
	    PyErr_SetString (PyExc_ValueError, buf);
	    return NULL;
	  }
	}
      }

      /* read surface identifier for the current face */
      surfaces [m ++] = PyInt_AsLong (PyList_GetItem (faces, i ++));
    }

    l = nver * 3;
    ERRMEM (ver = malloc (sizeof (double) * l));
    for (i = 0; i < l; i ++) ver [i] = PyFloat_AsDouble (PyList_GetItem (vertices, i));

    if (convex) cvx = convex->cvx; /* we are appending a list */

    self->cvx = CONVEX_Create (cvx, ver, nver, fac, nfac, surfaces, volid);

    if (cvx)
    {
      convex->cvx = NULL; /* empty */
    }

    /* clean up */
    free (surfaces);
    free (ver);
    free (fac);
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_CONVEX_dealloc (lng_CONVEX *self)
{
  if (self->cvx) CONVEX_Destroy (self->cvx);

  self->ob_type->tp_free ((PyObject*)self);
}

/* return vertex */
static PyObject* lng_CONVEX_vertex (lng_CONVEX *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("n");
  int n;

  PARSEKEYS ("i", &n);

  if (!self->cvx)
  {
    PyErr_SetString (PyExc_ValueError, "The CONVEX object is empty");
    return NULL;
  }
  else if (n < 0 || n >= self->cvx->nver)
  {
    PyErr_SetString (PyExc_ValueError, "Vertex index out of bounds");
    return NULL;
  }

  return Py_BuildValue ("(d, d, d)", self->cvx->cur [3*n], self->cvx->cur [3*n+1], self->cvx->cur [3*n+2]);
}

/* number of vertices */
static PyObject* lng_CONVEX_get_nver (lng_CONVEX *self, void *closure)
{
  if (!self->cvx)
  {
    PyErr_SetString (PyExc_ValueError, "The CONVEX object is empty");
    return NULL;
  }

  return PyInt_FromLong (self->cvx->nver);
}

static int lng_CONVEX_set_nver (lng_CONVEX *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* CONVEX methods */
static PyMethodDef lng_CONVEX_methods [] =
{
  {"vertex", (PyCFunction)lng_CONVEX_vertex, METH_VARARGS|METH_KEYWORDS, "Return a vertex point of a convex"},
  {NULL, NULL, 0, NULL}
};

/* CONVEX members */
static PyMemberDef lng_CONVEX_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* CONVEX getset */
static PyGetSetDef lng_CONVEX_getset [] =
{
  {"nver", (getter)lng_CONVEX_get_nver, (setter)lng_CONVEX_set_nver, "vertices count", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * MESH => object
 */

typedef struct lng_MESH lng_MESH;

static PyTypeObject lng_MESH_TYPE; /* type descriptor */

struct lng_MESH
{
  PyObject_HEAD
  MESH *msh;
};

/* test whether an object is of MESH type */
static int is_mesh (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance (obj, (PyObject*)&lng_MESH_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a MESH object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* constructor */
static PyObject* lng_MESH_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("nodes", "elements", "surfids");
  int *ele, *sur, i, j, k, l, n, m, nn, o;
  PyObject *nodes, *elements, *surfids;
  double (*nod) [3];
  lng_MESH *self;

  self = (lng_MESH*)type->tp_alloc (type, 0);

  if (self)
  {
    PARSEKEYS ("OOO", &nodes, &elements, &surfids);

    TYPETEST (is_list (nodes, kwl[0], 3, 12) && is_list (elements, kwl[1], 1, 0) && is_list_or_number (surfids, kwl[2], 1, 1));

    /* test element definitions */
    l = PyList_Size (elements);
    for (i = 0; i < l; i ++)
    {
      k = PyInt_AsLong (PyList_GetItem (elements, i));

      if (!(k == 4 || k == 5 || k == 6 || k == 8))
      {
	PyErr_SetString (PyExc_ValueError, "An element must have 4, 5, 6, or 8 nodes");
	return NULL;
      }

      /* add one more for the
       * volume definition */
      i += (k + 1);

      if (i >= l) /* incomplete */
      {
	PyErr_SetString (PyExc_ValueError, "The last element definition is incomplete");
	return NULL;
      }
    }

    /* read elements */
    ERRMEM (ele = malloc ((l + 1) * sizeof (int)));
    nn = PyList_Size (nodes) / 3; /* nodes count */

    for (m = n = i = 0; i < l; m ++)
    {
      ele [n ++] = k = PyInt_AsLong (PyList_GetItem (elements, i ++));

      for (j = 0; j < k; j ++)
      {
        ele [n ++] = PyInt_AsLong (PyList_GetItem (elements, i ++));

	if (ele  [n-1] < 0 || ele [n-1] >= nn) /* must be within the right range */
	{
	  char buf [BUFLEN];
	  sprintf (buf, "Node %d in element %d is outside of range [0, %d]",j , m, nn-1);
	  PyErr_SetString (PyExc_ValueError, buf);
	  return NULL;
	}
      }

      /* test for repeated nodes in element definition */
      for (j = 1; j <= k; j ++)
      {
	for (o = j + 1; o <= k; o ++)
	{
	  if (ele [n-j] == ele [n-o])
	  {
	    char buf [BUFLEN];
	    sprintf (buf, "Nodes %d and %d in element %d are the same", k-j, k-o, m);
	    PyErr_SetString (PyExc_ValueError, buf);
	    return NULL;
	  }
	}
      }

      ele [n ++] = PyInt_AsLong (PyList_GetItem (elements, i ++)); /* volid */
    }
    ele [n] = 0; /* end of list */

   if (PyList_Check (surfids))
   {
      /* test surface definitions */
      l = PyList_Size (surfids);
      for (i = 1; i < l; i ++)
      {
	k = PyInt_AsLong (PyList_GetItem (surfids, i));

	if (!(k == 3 || k == 4))
	{
	  PyErr_SetString (PyExc_ValueError, "A face must have 3 or 4 nodes");
	  return NULL;
	}

	/* add one more for the
	 * surface definition */
	i += (k + 1);

	if (i >= l) /* incomplete */
	{
	  PyErr_SetString (PyExc_ValueError, "The last face definition is incomplete");
	  return NULL;
	}
      }

      /* read surfaces */
      ERRMEM (sur = malloc ((l + 1) * sizeof (int)));
      sur [0] = PyInt_AsLong (PyList_GetItem (surfids, 0)); /* gid */

      for (m = 0, n = i = 1; i < l; m ++)
      {
	sur [n ++] = k = PyInt_AsLong (PyList_GetItem (surfids, i ++));

	for (j = 0; j < k; j ++)
	{
	  sur [n ++] = PyInt_AsLong (PyList_GetItem (surfids, i ++));

	  if (sur [n-1] < 0 || sur [n-1] >= nn) /* must be within the right range */
	  {
	    char buf [BUFLEN];
	    sprintf (buf, "Node %d in face %d is outside of range [0, %d]", j, m, nn-1);
	    PyErr_SetString (PyExc_ValueError, buf);
	    return NULL;
	  }
	}

	/* test for repeated nodes in face definition */
	for (j = 1; j <= k; j ++)
	{
	  for (o = j + 1; o <= k; o ++)
	  {
	    if (sur [n-j] == sur [n-o])
	    {
	      char buf [BUFLEN];
	      sprintf (buf, "Nodes %d and %d in face %d are the same", k-j, k-o, m);
	      PyErr_SetString (PyExc_ValueError, buf);
	      return NULL;
	    }
	  }
	}

	sur [n ++] = PyInt_AsLong (PyList_GetItem (surfids, i ++)); /* surfid */
      }
      sur [n] = 0; /* end of list */
   }
   else
   {
      ERRMEM (sur = malloc (2 * sizeof (int)));
      sur [0] = PyInt_AsLong (surfids);
      sur [1] = 0; /* end of list */
   }

    /* nodes */
    ERRMEM (nod = malloc (nn * sizeof (double [3])));
    for (i = 0; i < nn; i ++) 
      for (j = 0; j < 3; j ++)
	nod [i][j] = PyFloat_AsDouble (PyList_GetItem (nodes, 3*i + j));

    self->msh = MESH_Create (nod, ele, sur);

    /* clean up */
    free (nod);
    free (ele);
    free (sur);
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_MESH_dealloc (lng_MESH *self)
{
  if (self->msh) MESH_Destroy (self->msh);

  self->ob_type->tp_free ((PyObject*)self);
}

/* set volid for all faces. returns the volid set */
static PyObject* lng_MESH_setvolid (lng_MESH *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("volid");
  int v;
  ELEMENT *ele;
  PARSEKEYS ("i", &v);
  if (!self->msh)
  {
    PyErr_SetString (PyExc_ValueError, "The MESH object is empty");
    return NULL;
  }
  
  // modify bulk elements
  for (ele = self->msh->bulkeles; ele; ele = ele->next)
  {
    ele->volume = v;
  }
  
  // modify surface elements
  for (ele = self->msh->surfeles; ele; ele = ele->next)
  {
    ele->volume = v;
  }
  
  return Py_BuildValue ("i", v);
  
}

/* return node */
static PyObject* lng_MESH_node (lng_MESH *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("n", "x", "y", "z");
  double x, y, z;
  int n;

  x = y = z = DBL_MAX;

  PARSEKEYS ("i|ddd", &n, &x, &y, &z);

  if (!self->msh)
  {
    PyErr_SetString (PyExc_ValueError, "The MESH object is empty");
    return NULL;
  }
  else if (n < 0 || n >= self->msh->nodes_count)
  {
    PyErr_SetString (PyExc_ValueError, "Node index out of bounds");
    return NULL;
  }

  if (x != DBL_MAX) 
  {
    self->msh->ref_nodes [n][0] =
    self->msh->cur_nodes [n][0] = x;
  }
  if (y != DBL_MAX) 
  {
    self->msh->ref_nodes [n][1] =
    self->msh->cur_nodes [n][1] = y;
  }
  if (z != DBL_MAX) 
  {
    self->msh->ref_nodes [n][2] =
    self->msh->cur_nodes [n][2] = z;
  }

  return Py_BuildValue ("(d, d, d)", self->msh->cur_nodes [n][0], self->msh->cur_nodes[n][1], self->msh->cur_nodes[n][2]);
}

/* return nodes numbers on surface */
static PyObject* lng_MESH_nodes_on_surface (lng_MESH *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("surfid");
  PyObject *list, *obj;
  SET *set, *item;
  int surfid;
  FACE *fac;
  int n;

  PARSEKEYS ("i", &surfid);

  if (!self->msh)
  {
    PyErr_SetString (PyExc_ValueError, "The MESH object is empty");
    return NULL;
  }

  for (set = NULL, fac = self->msh->faces; fac; fac = fac->n)
  {
    if (fac->surface == surfid)
    {
      for (n = 0; n < fac->type; n ++)
      {
        SET_Insert (NULL, &set, (void*) (long) fac->nodes[n], NULL);
      }
    }
  }

  if (set == NULL) Py_RETURN_NONE;
  else
  {
    if (!(list = PyList_New (SET_Size (set)))) return NULL;

    for (item = SET_First (set), n = 0; item; item = SET_Next (item), n ++)
    {
      if (!(obj = PyInt_FromLong((long)item->data))) return NULL;

      PyList_SetItem (list, n, obj);
    }

    SET_Free (NULL, &set);
    return list;
  }
}

/* number of nodes */
static PyObject* lng_MESH_get_nnod (lng_MESH *self, void *closure)
{
  if (!self->msh)
  {
    PyErr_SetString (PyExc_ValueError, "The MESH object is empty");
    return NULL;
  }

  return PyInt_FromLong (self->msh->nodes_count);
}

static int lng_MESH_set_nnod (lng_MESH *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* return the nodes, elements and surfid data in the form that MESH() requires */
static PyObject* lng_MESH_get_meshdata (lng_MESH *self, PyObject *args, PyObject *kwds)
{
  PyObject *ellist, *nodelist, *surfidlist;
  ELEMENT *ele;
  FACE *fac;
  int n, d;
  double v;
  
  if (!self->msh)
  {
    PyErr_SetString (PyExc_ValueError, "The MESH object is empty");
    return NULL;
  }
  
  // extract nodes
  nodelist = PyList_New(0);
  for (n=0; n < self->msh->nodes_count; n++)
  {
    for (d=0; d < 3; d++)
    {
      v = self->msh->cur_nodes[n][d]; // use the CURRENT config, not the starting one
      PyList_Append(nodelist, PyFloat_FromDouble(v));
    }
  }

  // extract bulk elements
  ellist = PyList_New(0);
  for (ele = self->msh->bulkeles; ele; ele = ele->next)
  {
    PyList_Append(ellist, PyInt_FromLong(ele->type)); // number of nodes
    for (n = 0; n < ele->type; n ++)
    {
      PyList_Append(ellist, PyInt_FromLong(ele->nodes[n])); // node numbers
    }
    PyList_Append(ellist, PyInt_FromLong(ele->volume)); // volid
  }
  
  // extract suface elements 
  for (ele = self->msh->surfeles; ele; ele = ele->next)
  {
    PyList_Append(ellist, PyInt_FromLong(ele->type)); // number of nodes
    for (n = 0; n < ele->type; n ++)
    {
      PyList_Append(ellist, PyInt_FromLong(ele->nodes[n])); // node numbers
    }
    PyList_Append(ellist, PyInt_FromLong(ele->volume)); // volid
  }
  
  // extract surfids - assume GID is zero, if this is the only sid it is returned as [GID]
  surfidlist = PyList_New(1);
  PyList_SetItem(surfidlist, 0, PyInt_FromLong(0)); // set GID
  for (fac = self->msh->faces; fac; fac = fac->n)
  {
    if (fac->surface != 0){
      PyList_Append(surfidlist, PyInt_FromLong(fac->type)); // f1 etc
      for (n = 0; n < fac->type; n++)
      {
        PyList_Append(surfidlist, PyInt_FromLong(fac->nodes[n])); // n1 etc
      }
      PyList_Append(surfidlist, PyInt_FromLong(fac->surface)); // sid
    }
  }
  
  return PyTuple_Pack(3, nodelist, ellist, surfidlist);
}

/* MESH methods */
static PyMethodDef lng_MESH_methods [] =
{ 
  {"node", (PyCFunction)lng_MESH_node, METH_VARARGS|METH_KEYWORDS, "Return a node point of a mesh"},
  {"nodes_on_surface", (PyCFunction)lng_MESH_nodes_on_surface, METH_VARARGS|METH_KEYWORDS, "Return list of node numbers belonging to a given surface"},
  {"get_data", (PyCFunction)lng_MESH_get_meshdata, METH_NOARGS, "Return mesh data - very developmental!"},
  {"set_volid", (PyCFunction)lng_MESH_setvolid, METH_VARARGS|METH_KEYWORDS, "Set volid for all faces"},
  {NULL, NULL, 0, NULL}
};

/* MESH members */
static PyMemberDef lng_MESH_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* MESH getset */
static PyGetSetDef lng_MESH_getset [] =
{ 
  {"nnod", (getter)lng_MESH_get_nnod, (setter)lng_MESH_set_nnod, "nodes count", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * SPHERE => object
 */

typedef struct lng_SPHERE lng_SPHERE;

static PyTypeObject lng_SPHERE_TYPE;

struct lng_SPHERE
{
  PyObject_HEAD
  SPHERE *sph;
};

#if 0
static int is_sphere (lng_SPHERE *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_SPHERE_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a SPHERE object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}
#endif

/* constructor */
static PyObject* lng_SPHERE_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("center", "radius", "volid", "surfid");
  double radius, c [3];
  int volid, surfid;
  PyObject *center;
  lng_SPHERE *self;

  self = (lng_SPHERE*)type->tp_alloc (type, 0);

  if (self)
  {
    PARSEKEYS ("Odii", &center, &radius, &volid, &surfid);

    TYPETEST (is_tuple (center, kwl [0], 3));

    c [0] = PyFloat_AsDouble (PyTuple_GetItem (center, 0));
    c [1] = PyFloat_AsDouble (PyTuple_GetItem (center, 1));
    c [2] = PyFloat_AsDouble (PyTuple_GetItem (center, 2));

    self->sph = SPHERE_Create (c, radius, surfid, volid);
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_SPHERE_dealloc (lng_SPHERE *self)
{
  if (self->sph) SPHERE_Destroy (self->sph);

  self->ob_type->tp_free ((PyObject*)self);
}

/* center */
static PyObject* lng_SPHERE_get_center (lng_SPHERE *self, void *closure)
{
  if (!self->sph)
  {
    PyErr_SetString (PyExc_ValueError, "The SPHERE object is empty");
    return NULL;
  }

  return Py_BuildValue ("(d, d, d)", self->sph->cur_center [0], self->sph->cur_center [1], self->sph->cur_center [2]);
}

static int lng_SPHERE_set_center (lng_SPHERE *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* radius */
static PyObject* lng_SPHERE_get_radius (lng_SPHERE *self, void *closure)
{
  if (!self->sph)
  {
    PyErr_SetString (PyExc_ValueError, "The SPHERE object is empty");
    return NULL;
  }

  return PyFloat_FromDouble (self->sph->cur_radius);
}

static int lng_SPHERE_set_radius (lng_SPHERE *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* SPHERE methods */
static PyMethodDef lng_SPHERE_methods [] =
{ {NULL, NULL, 0, NULL} };

/* SPHERE members */
static PyMemberDef lng_SPHERE_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* SPHERE getset */
static PyGetSetDef lng_SPHERE_getset [] =
{ 
  {"center", (getter)lng_SPHERE_get_center, (setter)lng_SPHERE_set_center, "sphere center", NULL},
  {"radius", (getter)lng_SPHERE_get_radius, (setter)lng_SPHERE_set_radius, "sphere radius", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * ELLIP => object
 */

typedef struct lng_ELLIP lng_ELLIP;

static PyTypeObject lng_ELLIP_TYPE;

struct lng_ELLIP
{
  PyObject_HEAD
  ELLIP *eli;
};

#if 0
static int is_ellip (lng_ELLIP *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_ELLIP_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a ELLIP object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}
#endif

/* constructor */
static PyObject* lng_ELLIP_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("center", "radii", "volid", "surfid");
  double c [3], r [3];
  int volid, surfid;
  PyObject *center, *radii;
  lng_ELLIP *self;

  self = (lng_ELLIP*)type->tp_alloc (type, 0);

  if (self)
  {
    PARSEKEYS ("OOii", &center, &radii, &volid, &surfid);

    TYPETEST (is_tuple (center, kwl [0], 3) && is_tuple (radii, kwl [1], 3));

    c [0] = PyFloat_AsDouble (PyTuple_GetItem (center, 0));
    c [1] = PyFloat_AsDouble (PyTuple_GetItem (center, 1));
    c [2] = PyFloat_AsDouble (PyTuple_GetItem (center, 2));

    r [0] = PyFloat_AsDouble (PyTuple_GetItem (radii, 0));
    r [1] = PyFloat_AsDouble (PyTuple_GetItem (radii, 1));
    r [2] = PyFloat_AsDouble (PyTuple_GetItem (radii, 2));

    self->eli = ELLIP_Create (c, r, surfid, volid);
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_ELLIP_dealloc (lng_ELLIP *self)
{
  if (self->eli) ELLIP_Destroy (self->eli);

  self->ob_type->tp_free ((PyObject*)self);
}

/* center */
static PyObject* lng_ELLIP_get_center (lng_ELLIP *self, void *closure)
{
  if (!self->eli)
  {
    PyErr_SetString (PyExc_ValueError, "The ELLIP object is empty");
    return NULL;
  }

  return Py_BuildValue ("(d, d, d)", self->eli->cur_center [0], self->eli->cur_center [1], self->eli->cur_center [2]);
}

static int lng_ELLIP_set_center (lng_ELLIP *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* radii */
static PyObject* lng_ELLIP_get_radii (lng_ELLIP *self, void *closure)
{
  if (!self->eli)
  {
    PyErr_SetString (PyExc_ValueError, "The ELLIP object is empty");
    return NULL;
  }

  return Py_BuildValue ("(d, d, d)", self->eli->cur_sca [0], self->eli->cur_sca [1], self->eli->cur_sca [2]);
}

static int lng_ELLIP_set_radii (lng_ELLIP *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* rotation */
static PyObject* lng_ELLIP_get_rot (lng_ELLIP *self, void *closure)
{
  if (!self->eli)
  {
    PyErr_SetString (PyExc_ValueError, "The ELLIP object is empty");
    return NULL;
  }

  return Py_BuildValue ("(d, d, d, d, d, d, d, d, d)",
    self->eli->cur_rot [0], self->eli->cur_rot [1], self->eli->cur_rot [2],
    self->eli->cur_rot [3], self->eli->cur_rot [4], self->eli->cur_rot [5],
    self->eli->cur_rot [6], self->eli->cur_rot [7], self->eli->cur_rot [8]);
}

static int lng_ELLIP_set_rot (lng_ELLIP *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* ELLIP methods */
static PyMethodDef lng_ELLIP_methods [] =
{ {NULL, NULL, 0, NULL} };

/* ELLIP members */
static PyMemberDef lng_ELLIP_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* ELLIP getset */
static PyGetSetDef lng_ELLIP_getset [] =
{ 
  {"center", (getter)lng_ELLIP_get_center, (setter)lng_ELLIP_set_center, "ellipsoid center", NULL},
  {"radii", (getter)lng_ELLIP_get_radii, (setter)lng_ELLIP_set_radii, "ellipsoid radii", NULL},
  {"rot", (getter)lng_ELLIP_get_rot, (setter)lng_ELLIP_set_rot, "ellipsoid rotation", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * SOLFEC => object
 */

static PyObject* lng_CONSTRAINT_WRAPPER (CON *con); /* constraint wrapper constructor */

static PyObject* lng_BODY_WRAPPER (BODY *bod); /* body wrapper constructor */

typedef struct lng_SOLFEC lng_SOLFEC; /* solfec type */

static PyTypeObject lng_SOLFEC_TYPE; /* type descriptor */

struct lng_SOLFEC /* solfec object */
{
  PyObject_HEAD
  SOLFEC *sol;
};

/* test whether an object is of SOLFEC type */
static int is_solfec (lng_SOLFEC *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_SOLFEC_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a SOLFEC object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* constructor */
static PyObject* lng_SOLFEC_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("analysis", "step", "output");
  PyObject *analysis, *output;
  lng_SOLFEC *self;
  char *outpath;
  double step;

  self = (lng_SOLFEC*)type->tp_alloc (type, 0);

  if (self)
  {
    step = 1E-3;
    outpath = "solfec.out";
    output = NULL;

    PARSEKEYS ("OdO|O", &analysis, &step, &output);

    TYPETEST (is_string (analysis, kwl [0]) && is_positive (step, kwl[1]) && is_string (output, kwl [2]));

    outpath = PyString_AsString (output);

    IFIS (analysis, "QUASI_STATIC")
    {
      self->sol = SOLFEC_Create (0, step, outpath);
      REGISTER_SOLFEC (self->sol);
#if OPENGL
      if (RND_Is_On ()) RND_Domain (self->sol->dom); /* just in case a viewer is enabled (last created SOLFEC object) */
#endif
    }
    ELIF (analysis, "DYNAMIC")
    {
      self->sol = SOLFEC_Create (1, step, outpath);
      REGISTER_SOLFEC (self->sol);
#if OPENGL
      if (RND_Is_On ()) RND_Domain (self->sol->dom); /* pass last created SOLFEC object to the rendering */
#endif
    }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid analysis kind");
      return NULL;
    }
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_SOLFEC_dealloc (lng_SOLFEC *self)
{
#if OPENGL
  if (RND_Is_On ()) return; /* do not delete in viewer mode */
  else
#endif
  SOLFEC_Destroy (self->sol);

  self->ob_type->tp_free ((PyObject*)self);
}

/* setgets */

static PyObject* lng_SOLFEC_get_analysis (lng_SOLFEC *self, void *closure)
{
  return PyString_FromString (self->sol->dom->dynamic ? "DYNAMIC" : "QUASI_STATIC");
}

static int lng_SOLFEC_set_analysis (lng_SOLFEC *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SOLFEC_get_time (lng_SOLFEC *self, void *closure)
{
  return PyFloat_FromDouble (self->sol->dom->time);
}

static int lng_SOLFEC_set_time (lng_SOLFEC *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SOLFEC_get_mode (lng_SOLFEC *self, void *closure)
{
  return PyString_FromString (SOLFEC_Mode (self->sol));
}

static int lng_SOLFEC_set_mode (lng_SOLFEC *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SOLFEC_get_constraints (lng_SOLFEC *self, void *closure)
{
  PyObject *list, *obj;
  CON *con;
  int n;

  if (!(list = PyList_New (self->sol->dom->ncon))) return NULL;

  for (n = 0, con = self->sol->dom->con; con; n ++, con = con->next)
  {
    if (!(obj = lng_CONSTRAINT_WRAPPER (con))) return NULL;

    PyList_SetItem (list, n, obj);
  }

  return list;
}

static int lng_SOLFEC_set_constraints (lng_SOLFEC *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SOLFEC_get_ncon (lng_SOLFEC *self, void *closure)
{
  return PyInt_FromLong (self->sol->dom->ncon);
}

static int lng_SOLFEC_set_ncon (lng_SOLFEC *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SOLFEC_get_bodies (lng_SOLFEC *self, void *closure)
{
  PyObject *list, *obj;
  BODY *bod;
  int n;

  if (!(list = PyList_New (self->sol->dom->nbod))) return NULL;

  for (n = 0, bod = self->sol->dom->bod; bod; n ++, bod = bod->next)
  {
    if (!(obj = lng_BODY_WRAPPER (bod))) return NULL;

    PyList_SetItem (list, n, obj);
  }

  return list;
}

static int lng_SOLFEC_set_bodies (lng_SOLFEC *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SOLFEC_get_nbod (lng_SOLFEC *self, void *closure)
{
  return PyInt_FromLong (self->sol->dom->nbod);
}

static int lng_SOLFEC_set_nbod (lng_SOLFEC *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SOLFEC_get_step (lng_SOLFEC *self, void *closure)
{
  return PyFloat_FromDouble (self->sol->dom->step);
}

static int lng_SOLFEC_set_step (lng_SOLFEC *self, PyObject *value, void *closure)
{
  double step;

  if (!is_number_gt (value, "step", 0)) return -1;
  step = PyFloat_AsDouble (value);
  self->sol->dom->step = step;

  return 0;
}

static PyObject* lng_SOLFEC_get_verbose (lng_SOLFEC *self, void *closure)
{
  if (self->sol->verbose) return PyString_FromString ("ON");
  else return PyString_FromString ("OFF");
}

static int lng_SOLFEC_set_verbose (lng_SOLFEC *self, PyObject *value, void *closure)
{
  if (!is_string (value, "verbose")) return -1;

  IFIS (value, "ON") self->sol->verbose = 1;
  ELIF (value, "OFF") self->sol->verbose = 0;
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid verbose value (ON/OFF accepted)");
    return -1;
  }

  return 0;
}

static PyObject* lng_SOLFEC_get_outpath (lng_SOLFEC *self, void *closure)
{
  return PyString_FromString (self->sol->outpath);
}

static int lng_SOLFEC_set_outpath (lng_SOLFEC *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* SOLFEC methods */
static PyMethodDef lng_SOLFEC_methods [] =
{ {NULL, NULL, 0, NULL} };

/* SOLFEC members */
static PyMemberDef lng_SOLFEC_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* SOLFEC getset */
static PyGetSetDef lng_SOLFEC_getset [] =
{ 
  {"analysis", (getter)lng_SOLFEC_get_analysis, (setter)lng_SOLFEC_set_analysis, "analysis kind", NULL},
  {"time", (getter)lng_SOLFEC_get_time, (setter)lng_SOLFEC_set_time, "current time", NULL},
  {"mode", (getter)lng_SOLFEC_get_mode, (setter)lng_SOLFEC_set_mode, "analysis mode", NULL},
  {"constraints", (getter)lng_SOLFEC_get_constraints, (setter)lng_SOLFEC_set_constraints, "list of constraints", NULL},
  {"ncon", (getter)lng_SOLFEC_get_ncon, (setter)lng_SOLFEC_set_ncon, "constraints count", NULL},
  {"bodies", (getter)lng_SOLFEC_get_bodies, (setter)lng_SOLFEC_set_bodies, "list of bodies", NULL},
  {"nbod", (getter)lng_SOLFEC_get_nbod, (setter)lng_SOLFEC_set_nbod, "bodies count", NULL},
  {"step", (getter)lng_SOLFEC_get_step, (setter)lng_SOLFEC_set_step, "time step", NULL},
  {"verbose", (getter)lng_SOLFEC_get_verbose, (setter)lng_SOLFEC_set_verbose, "verbosity", NULL},
  {"outpath", (getter)lng_SOLFEC_get_outpath, (setter)lng_SOLFEC_set_outpath, "verbosity", NULL},
  {NULL, 0, 0, NULL, NULL} 
};

/*
 * FIELD => object
 */

typedef struct lng_FIELD lng_FIELD; /* surface material type */

static PyTypeObject lng_FIELD_TYPE; /* type descriptor */

struct lng_FIELD /* surface material object */
{
  PyObject_HEAD
  FIELD *fld;
};

/* constructor */
static PyObject* lng_FIELD_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS (
    "solfec",
    "filed_callback",
    "label",
    "data");

  PyObject *callback, *label, *data;
  lng_SOLFEC *solfec;
  lng_FIELD *self;
  SOLFEC *sol;
  FIELD fld;

  self = (lng_FIELD*)type->tp_alloc (type, 0);

  if (self)
  {
    label = NULL;
    data = NULL;

    PARSEKEYS ("OO|OO", &solfec, &callback, &label, &data);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_callable (callback, kwl[1]) && is_string (label, kwl[2]));

    sol = solfec->sol;
    fld.call = callback;
    fld.data = data;

    self->fld = FISET_Insert (sol->fis, as_string (label), fld); /* insert data into field set */
  }

  return  (PyObject*)self;
}

/* create material wrapper */
static PyObject* lng_FIELD_WRAPPER (FIELD *fld)
{
  lng_FIELD *self;

  self = (lng_FIELD*)lng_FIELD_TYPE.tp_alloc (&lng_FIELD_TYPE, 0);

  self->fld = fld;

  return (PyObject*)self;
}

/* destructor */
static void lng_FIELD_dealloc (lng_FIELD *self)
{
  self->ob_type->tp_free ((PyObject*)self);
}

/* setgets */

static PyObject* lng_FIELD_get_label (lng_FIELD *self, void *closure)
{
  return PyString_FromString (self->fld->label);
}

static int lng_FIELD_set_label (lng_FIELD *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* FIELD methods */
static PyMethodDef lng_FIELD_methods [] =
{ {NULL, NULL, 0, NULL} };

/* FIELD members */
static PyMemberDef lng_FIELD_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* FIELD getset */
static PyGetSetDef lng_FIELD_getset [] =
{ 
  {"label", (getter)lng_FIELD_get_label, (setter)lng_FIELD_set_label, "label", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/* test whether an object is a bulk material or a bulk material label */
static int is_field (SOLFEC *sol, PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_FIELD_TYPE))
    {
      if (PyString_Check (obj) && FISET_Find (sol->fis, PyString_AsString (obj))) return 1;

      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a FIELD object or a valid field label", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* get field related to an object or a label */
static FIELD* get_field (SOLFEC *sol, PyObject *obj)
{
  if (obj)
  {
    if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_FIELD_TYPE))
    {
      lng_FIELD *f = (lng_FIELD*)obj;
      return f->fld;
    }
    else return FISET_Find (sol->fis, PyString_AsString (obj));
  }

  return NULL;
}

/*
 * SURFACE_MATERIAL => object
 */

typedef struct lng_SURFACE_MATERIAL lng_SURFACE_MATERIAL; /* surface material type */

static PyTypeObject lng_SURFACE_MATERIAL_TYPE; /* type descriptor */

struct lng_SURFACE_MATERIAL /* surface material object */
{
  PyObject_HEAD
  SURFACE_MATERIAL *mat;
};

/* constructor */
static PyObject* lng_SURFACE_MATERIAL_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS (
    "solfec",
    "surf1",
    "surf2",
    "model",
    "label",
    "friction",
    "cohesion",
    "restitution",
    "spring",
    "dashpot",
    "hpow");

  SURFACE_MATERIAL mat = 
  {
   0,
   INT_MAX,    /* surf1 */
   INT_MAX,    /* surf2 */
   NULL, /* label */
   SIGNORINI_COULOMB, /* model */
   0.0,  /* friction */
   0.0,  /* cohesion */
   0.0,  /* restitution */
   0.0,  /* spring */
   0.0,  /* dashpot */
   1.0,  /* Hertz power */
  };

  PyObject *model, *label;
  lng_SURFACE_MATERIAL *self;
  lng_SOLFEC *solfec;
  SOLFEC *sol;

  self = (lng_SURFACE_MATERIAL*)type->tp_alloc (type, 0);

  if (self)
  {
    model = NULL;
    label = NULL;

    PARSEKEYS ("O|iiOOdddddd", &solfec, &mat.surf1, &mat.surf2, &model, &label,
      &mat.friction, &mat.cohesion, &mat.restitution, &mat.spring, &mat.dashpot, &mat.hpow);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_string (model, kwl[3]) && is_string (label, kwl[4]) &&
	      is_non_negative (mat.friction, kwl[5]) && is_non_negative (mat.cohesion, kwl[6]) &&
	      is_in_range (mat.restitution, kwl[7], 0, 1) && is_non_negative (mat.spring, kwl[8]) &&
              is_non_negative (mat.hpow, kwl[10]));

    sol = solfec->sol;

    if (model)
    {
      IFIS (model, "SIGNORINI_COULOMB") mat.model = SIGNORINI_COULOMB;
      ELIF (model, "SPRING_DASHPOT") mat.model = SPRING_DASHPOT;
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Unknown SURFACE_MATERIAL model");
	return NULL;
      }
    }

    if (mat.surf1 != INT_MAX && mat.surf2 != INT_MAX) /* not a default material */
      self->mat = SPSET_Insert (sol->sps, mat.surf1, mat.surf2, as_string (label), mat); /* insert data into the surface pair set */
    else SPSET_Default (sol->sps, mat);
  }

  return  (PyObject*)self;
}

/* create material wrapper */
static PyObject* lng_SURFACE_MATERIAL_WRAPPER (SURFACE_MATERIAL *mat)
{
  lng_SURFACE_MATERIAL *self;

  self = (lng_SURFACE_MATERIAL*)lng_SURFACE_MATERIAL_TYPE.tp_alloc (&lng_SURFACE_MATERIAL_TYPE, 0);

  self->mat = mat;

  return (PyObject*)self;
}

/* destructor */
static void lng_SURFACE_MATERIAL_dealloc (lng_SURFACE_MATERIAL *self)
{
  self->ob_type->tp_free ((PyObject*)self);
}

/* setgets */

static PyObject* lng_SURFACE_MATERIAL_get_surf1 (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyInt_FromLong (self->mat->surf1);
}

static int lng_SURFACE_MATERIAL_set_surf1 (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SURFACE_MATERIAL_get_surf2 (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyInt_FromLong (self->mat->surf2);
}

static int lng_SURFACE_MATERIAL_set_surf2 (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SURFACE_MATERIAL_get_label (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyString_FromString (self->mat->label);
}

static int lng_SURFACE_MATERIAL_set_label (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_SURFACE_MATERIAL_get_model (lng_SURFACE_MATERIAL *self, void *closure)
{
  switch (self->mat->model)
  {
  case SIGNORINI_COULOMB:
    return PyString_FromString ("SIGNORINI_COULOMB");
  case SPRING_DASHPOT:
    return PyString_FromString ("SPRING_DASHPOT");
  }

  return NULL;
}

static int lng_SURFACE_MATERIAL_set_model (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_string (value, "model")) return -1;

  IFIS (value, "SIGNORINI_COULOMB") self->mat->model = SIGNORINI_COULOMB;
  ELIF (value, "SPRING_DASHPOT") self->mat->model = SPRING_DASHPOT;
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Unknown SURFACE_MATERIAL model");
    return -1;
  }

  return 0;
}

static PyObject* lng_SURFACE_MATERIAL_get_friction (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->friction);
}

static int lng_SURFACE_MATERIAL_set_friction (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "friction", 0)) return -1;
  self->mat->friction = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_SURFACE_MATERIAL_get_cohesion (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->cohesion);
}

static int lng_SURFACE_MATERIAL_set_cohesion (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "cohesion", 0)) return -1;
  self->mat->cohesion = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_SURFACE_MATERIAL_get_restitution (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->restitution);
}

static int lng_SURFACE_MATERIAL_set_restitution (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "restitution", 0)) return -1;
  self->mat->restitution = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_SURFACE_MATERIAL_get_spring (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->spring);
}

static int lng_SURFACE_MATERIAL_set_spring (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "spring", 0)) return -1;
  self->mat->spring = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_SURFACE_MATERIAL_get_dashpot (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->dashpot);
}

static int lng_SURFACE_MATERIAL_set_dashpot (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "dashpot", 0)) return -1;
  self->mat->dashpot = PyFloat_AsDouble (value);
  return 0;
}

/* SURFACE_MATERIAL methods */
static PyMethodDef lng_SURFACE_MATERIAL_methods [] =
{ {NULL, NULL, 0, NULL} };

/* SURFACE_MATERIAL members */
static PyMemberDef lng_SURFACE_MATERIAL_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* SURFACE_MATERIAL getset */
static PyGetSetDef lng_SURFACE_MATERIAL_getset [] =
{ 
  {"surf1", (getter)lng_SURFACE_MATERIAL_get_surf1, (setter)lng_SURFACE_MATERIAL_set_surf1, "surf1", NULL},
  {"surf2", (getter)lng_SURFACE_MATERIAL_get_surf2, (setter)lng_SURFACE_MATERIAL_set_surf2, "surf2", NULL},
  {"label", (getter)lng_SURFACE_MATERIAL_get_label, (setter)lng_SURFACE_MATERIAL_set_label, "label", NULL},
  {"model", (getter)lng_SURFACE_MATERIAL_get_model, (setter)lng_SURFACE_MATERIAL_set_model, "model", NULL},
  {"friction", (getter)lng_SURFACE_MATERIAL_get_friction, (setter)lng_SURFACE_MATERIAL_set_friction, "friction", NULL},
  {"cohesion", (getter)lng_SURFACE_MATERIAL_get_cohesion, (setter)lng_SURFACE_MATERIAL_set_cohesion, "cohesion", NULL},
  {"restitution", (getter)lng_SURFACE_MATERIAL_get_restitution, (setter)lng_SURFACE_MATERIAL_set_restitution, "restitution", NULL},
  {"spring", (getter)lng_SURFACE_MATERIAL_get_spring, (setter)lng_SURFACE_MATERIAL_set_spring, "spring", NULL},
  {"dashpot", (getter)lng_SURFACE_MATERIAL_get_dashpot, (setter)lng_SURFACE_MATERIAL_set_dashpot, "dashpot", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * BULK_MATERIAL => object
 */

typedef struct lng_BULK_MATERIAL lng_BULK_MATERIAL; /* surface material type */

static PyTypeObject lng_BULK_MATERIAL_TYPE; /* type descriptor */

struct lng_BULK_MATERIAL /* surface material object */
{
  PyObject_HEAD
  BULK_MATERIAL *mat;
};

/* constructor */
static PyObject* lng_BULK_MATERIAL_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS (
    "solfec",
    "model",
    "label",
    "young",
    "poisson",
    "density",
    "tensile",
    "fields",
    "fracene");

  BULK_MATERIAL mat;

  mat.label = NULL;
  mat.model = KIRCHHOFF;
  mat.young = 1E9;
  mat.poisson = 0.25;
  mat.density = 1E3;
  mat.tensile = DBL_MAX;
  mat.fracene = DBL_MAX;
  mat.umat = NULL;
  mat.nstate = 0;
  mat.nfield = 0;

  PyObject *model, *label, *fields;
  lng_BULK_MATERIAL *self;
  lng_SOLFEC *solfec;
  SOLFEC *sol;

  self = (lng_BULK_MATERIAL*)type->tp_alloc (type, 0);

  if (self)
  {
    model = NULL;
    label = NULL;
    fields = NULL;

    PARSEKEYS ("O|OOddddOd", &solfec, &model, &label, &mat.young, &mat.poisson, &mat.density, &mat.tensile, &fields, &mat.fracene);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_string (model, kwl[1]) && is_string (label, kwl[2]) &&
	      is_non_negative (mat.young, kwl[3]) && is_non_negative (mat.poisson, kwl[4]) &&
	      is_non_negative (mat.density, kwl[5]) && is_non_negative (mat.tensile, kwl[6]) &&
	      is_list (fields, kwl [7], 1, 1) && is_non_negative (mat.fracene, kwl[8]));

    sol = solfec->sol;

    if (model)
    {
      IFIS (model, "KIRCHHOFF") mat.model = KIRCHHOFF;
      ELIF (model, "TSANG_MARSDEN")
      {
	mat.model = TSANG_MARSDEN;
	mat.nstate = 19;
	mat.nfield = 4;

	if (!fields || PyList_Size (fields) != 4)
	{
	  PyErr_SetString (PyExc_ValueError, "4 fields must be given");
	  return NULL;
	}

	for (int i = 0; i < 4; i ++)
	{
	  PyObject *o = PyList_GetItem (fields, i);
	  if (!is_field (sol, o, "fields []")) return NULL; /* sets Err string internally */
	  mat.fld [i] = get_field (sol, o);
	}
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Unknown BULK_MATERIAL model");
	return NULL;
      }
    }

    self->mat = MATSET_Insert (sol->mat, as_string (label), mat); /* insert data into bulk material set */
  }

  return  (PyObject*)self;
}

/* create material wrapper */
static PyObject* lng_BULK_MATERIAL_WRAPPER (BULK_MATERIAL *mat)
{
  lng_BULK_MATERIAL *self;

  self = (lng_BULK_MATERIAL*)lng_BULK_MATERIAL_TYPE.tp_alloc (&lng_BULK_MATERIAL_TYPE, 0);

  self->mat = mat;

  return (PyObject*)self;
}

/* destructor */
static void lng_BULK_MATERIAL_dealloc (lng_BULK_MATERIAL *self)
{
  self->ob_type->tp_free ((PyObject*)self);
}

/* setgets */

static PyObject* lng_BULK_MATERIAL_get_label (lng_BULK_MATERIAL *self, void *closure)
{
  return PyString_FromString (self->mat->label);
}

static int lng_BULK_MATERIAL_set_label (lng_BULK_MATERIAL *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BULK_MATERIAL_get_model (lng_BULK_MATERIAL *self, void *closure)
{
  switch (self->mat->model)
  {
  case KIRCHHOFF:
    return PyString_FromString ("KIRCHHOFF");
  case TSANG_MARSDEN:
    return PyString_FromString ("TSANG_MARSDEN");
  }

  return NULL;
}

static int lng_BULK_MATERIAL_set_model (lng_BULK_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_string (value, "model")) return -1;

  IFIS (value, "KIRCHHOFF") self->mat->model = KIRCHHOFF;
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Unknown BULK_MATERIAL model");
    return -1;
  }

  return 0;
}

static PyObject* lng_BULK_MATERIAL_get_young (lng_BULK_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->young);
}

static int lng_BULK_MATERIAL_set_young (lng_BULK_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "young", 0)) return -1;
  self->mat->young = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_BULK_MATERIAL_get_poisson (lng_BULK_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->poisson);
}

static int lng_BULK_MATERIAL_set_poisson (lng_BULK_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "poisson", 0)) return -1;
  self->mat->poisson = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_BULK_MATERIAL_get_density (lng_BULK_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->density);
}

static int lng_BULK_MATERIAL_set_density (lng_BULK_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "density", 0)) return -1;
  self->mat->density = PyFloat_AsDouble (value);
  return 0;
}

/* BULK_MATERIAL methods */
static PyMethodDef lng_BULK_MATERIAL_methods [] =
{ {NULL, NULL, 0, NULL} };

/* BULK_MATERIAL members */
static PyMemberDef lng_BULK_MATERIAL_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* BULK_MATERIAL getset */
static PyGetSetDef lng_BULK_MATERIAL_getset [] =
{ 
  {"label", (getter)lng_BULK_MATERIAL_get_label, (setter)lng_BULK_MATERIAL_set_label, "label", NULL},
  {"model", (getter)lng_BULK_MATERIAL_get_model, (setter)lng_BULK_MATERIAL_set_model, "model", NULL},
  {"young", (getter)lng_BULK_MATERIAL_get_young, (setter)lng_BULK_MATERIAL_set_young, "young", NULL},
  {"poisson", (getter)lng_BULK_MATERIAL_get_poisson, (setter)lng_BULK_MATERIAL_set_poisson, "poisson", NULL},
  {"density", (getter)lng_BULK_MATERIAL_get_density, (setter)lng_BULK_MATERIAL_set_density, "density", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * BODY => object
 */

typedef struct lng_BODY lng_BODY; /* body type */

static PyTypeObject lng_BODY_TYPE;

struct lng_BODY 
{
  PyObject_HEAD

#if MPI
  unsigned int id;

  DOM *dom;
#endif

  BODY *bod;
};

#if MPI
/* is a body on this processor? */
static int IS_HERE (lng_BODY *body)
{
  if (body->id)
  {
    BODY *bod = MAP_Find (body->dom->idb, (void*) (long) body->id, NULL);

    if (bod)
    {
      body->bod = bod;
      return 1;
    }
  }

  return 0;
}
#endif

/* test whether an object is a bulk material or a bulk material label */
static int is_bulk_material (SOLFEC *sol, PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_BULK_MATERIAL_TYPE))
    {
      if (PyString_Check (obj) && MATSET_Find (sol->mat, PyString_AsString (obj))) return 1;

      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a BULK_MATERAIL object or a valid material label", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether an object is a basic shape */
static int is_basic_shape (PyObject *obj)
{
  if (obj)
  {
    if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_CONVEX_TYPE))
    {
      lng_CONVEX *convex = (lng_CONVEX*) obj;
      if (convex->cvx == NULL) return 0; /* empty */
    }
    else if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_MESH_TYPE))
    {
      lng_MESH *mesh = (lng_MESH*) obj;
      if (mesh->msh == NULL) return 0; /* empty */
    }
    else if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_SPHERE_TYPE))
    {
      lng_SPHERE *sphere = (lng_SPHERE*) obj;
      if (sphere->sph == NULL) return 0; /* empty */
    }
    else if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_ELLIP_TYPE))
    {
      lng_ELLIP *ellip = (lng_ELLIP*) obj;
      if (ellip->eli == NULL) return 0; /* empty */
    }
    else return 0;
  }

  return 1;
}

/* test whether an object is a valid shape representation */
static int is_shape (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!is_basic_shape (obj))
    {
      if (PyList_Check (obj))
      {
	int i, n = PyList_Size (obj);

	for (i = 0; i < n; i ++)
	  if (!is_basic_shape (PyList_GetItem (obj, i))) break;

	if (i == n) return 1;
      }

      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a non-empty CONVEX/MESH/SPHERE/ELLIP object or a list of those", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether an object is a valid shape composed of convices */
static int is_shape_convex (PyObject *obj, char *var)
{
  if (obj)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be a non-empty CONVEX object or a list of those", var);

    if (!is_basic_shape (obj)) /* may be a list */
    {
      if (PyList_Check (obj))
      {
	int i, n = PyList_Size (obj);

	for (i = 0; i < n; i ++)
	  if (!is_convex_test (PyList_GetItem (obj, i))) break;

	if (i == n) return 1;
      }

      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
    else if (!is_convex_test (obj)) /* must be a convex */
    {
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether an object is a valid shape or a 3-vector */
static int is_shape_or_vector (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!is_basic_shape (obj))
    {
      if (PyList_Check (obj))
      {
	int i, n = PyList_Size (obj);

	for (i = 0; i < n; i ++)
	  if (!is_basic_shape (PyList_GetItem (obj, i))) break;

	if (i == n) return 1;
      }

      if (!PyTuple_Check (obj))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must be a non-empty CONVEX/MESH/SPHERE/ELLIP object, a list of those or a (x, y, z) tuple", var);
	PyErr_SetString (PyExc_TypeError, buf);
	return 0;
      }

      if (PyTuple_Size (obj) != 3)
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must have %d elements", var, 3);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* is this a BODY object? */
static int is_body_check (PyObject *obj)
{
  if (PyObject_IsInstance (obj, (PyObject*)&lng_BODY_TYPE))
  {
    lng_BODY *body = (lng_BODY*)obj;
    return body->bod != NULL;
  }
  else return 0;
}

/* test whether an object is a valid shape or a body object */
static int is_shape_or_body (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!is_basic_shape (obj))
    {
      if (PyList_Check (obj))
      {
	int i, n = PyList_Size (obj);

	for (i = 0; i < n; i ++)
	  if (!is_basic_shape (PyList_GetItem (obj, i))) break;

	if (i == n) return 1;
      }

      if (!is_body_check (obj))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must be a non-empty CONVEX/MESH/SPHERE/ELLIP object, a list of those or a BODY object", var);
	PyErr_SetString (PyExc_TypeError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* get bulk material related to an object or a label */
static BULK_MATERIAL* get_bulk_material (SOLFEC *sol, PyObject *obj)
{
  if (obj)
  {
    if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_BULK_MATERIAL_TYPE))
    {
      lng_BULK_MATERIAL *mat = (lng_BULK_MATERIAL*)obj;
      return mat->mat;
    }
    else return MATSET_Find (sol->mat, PyString_AsString (obj));
  }

  return NULL;
}

/* return basic shape kind */
static int shape_kind (PyObject *obj)
{
  if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_CONVEX_TYPE)) return SHAPE_CONVEX;
  else if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_MESH_TYPE)) return SHAPE_MESH;
  else if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_SPHERE_TYPE)) return SHAPE_SPHERE;
  else if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_ELLIP_TYPE)) return SHAPE_ELLIP;

  return -1;
}

/* return basic shape and empty the container */
static void* get_shape (PyObject *obj, short empty)
{
  void *out;

  out = NULL;

  if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_CONVEX_TYPE))
  {
    lng_CONVEX *convex = (lng_CONVEX*)obj;
    out = convex->cvx;
    if (empty) convex->cvx = NULL; /* empty */
  }
  else if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_MESH_TYPE))
  {
    lng_MESH *mesh = (lng_MESH*)obj;
    out = mesh->msh;
    if (empty) mesh->msh = NULL; /* empty */
  }
  else if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_SPHERE_TYPE))
  {
    lng_SPHERE *sphere = (lng_SPHERE*)obj;
    out = sphere->sph;
    if (empty) sphere->sph = NULL; /* empty */
  }
  else if (PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_ELLIP_TYPE))
  {
    lng_ELLIP *ellip = (lng_ELLIP*)obj;
    out = ellip->eli;
    if (empty) ellip->eli = NULL; /* empty */
  }

  return out;
}

/* create a shape out of basic shapes */
static SHAPE* create_shape (PyObject *obj, short empty)
{
  if (obj)
  {
    if (PyList_Check (obj))
    {
      int i, n = PyList_Size (obj);
      SHAPE *out = NULL;
      PyObject *it;

      for (i = 0; i < n; i ++)
      {
	it = PyList_GetItem (obj, i);
	if (empty > 0) out = SHAPE_Glue (SHAPE_Create (shape_kind (it), get_shape (it, 1)), out); /* empty and glue simple shapes (destructive for simple shape lists) */
	else if (empty < 0) out = SHAPE_Glue_Simple (SHAPE_Create (shape_kind (it), get_shape (it, 1)), out); /* empty and do not glue simple shape (non destructive) */
	else out = SHAPE_Glue_Simple (SHAPE_Create (shape_kind (it), get_shape (it, 0)), out); /* do not empty and do not glue simple shape (non destructive) */
      }

      return out;
    }
    else return SHAPE_Create (shape_kind (obj), get_shape (obj, empty));
  }

  return NULL;
}

/* test whether an object is of BODY type */
static int is_body (lng_BODY *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_BODY_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a BODY object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
    else
    {
      lng_BODY *body = (lng_BODY*)obj;
      if (body->bod == NULL)
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' is a recently deleted BODY object", var);
	PyErr_SetString (PyExc_TypeError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* body object constructor */
static PyObject* lng_BODY_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "kind", "shape", "material", "label", "form", "mesh", "modal");
  PyObject *kind, *shape, *material, *label, *formulation, *modal;
  lng_SOLFEC *solfec;
  lng_BODY *self;
  lng_MESH *mesh;
  MESH *msh;
  short form;
  char *lab;
  MX *E; /* modal base */
  double *val; /* modal eigenvectors */

  self = (lng_BODY*)type->tp_alloc (type, 0);

#if MPI && LOCAL_BODIES
  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  if (self && rank != 0) /* all bodies created on rank 0 */
  {
    PARSEKEYS ("OOOO|OOOO", &solfec, &kind, &shape, &material, &label, &formulation, &mesh, &modal);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_string (kind, kwl[1]) &&
	      is_shape (shape, kwl[2]) && is_bulk_material (solfec->sol, material, kwl[3]));

    self->dom = solfec->sol->dom;
    self->id = self->dom->bid ++;
    self->bod = (BODY*)1; /* XXX is_body will not complain */
  }
  else
#endif
  if (self)
  {
    label = NULL;
    formulation = NULL;
    form = TOTAL_LAGRANGIAN;
    mesh = NULL;
    msh = NULL;
    modal = NULL;
    E = NULL;
    val = NULL;

    PARSEKEYS ("OOOO|OOOO", &solfec, &kind, &shape, &material, &label, &formulation, &mesh, &modal);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_string (kind, kwl[1]) && is_shape (shape, kwl[2]) &&
	      is_bulk_material (solfec->sol, material, kwl[3]) && is_string (label, kwl[4]) &&
	      is_string (formulation, kwl[5]) && is_mesh ((PyObject*)mesh, kwl[6]) && is_tuple (modal, kwl[7], 2));

    if (label) lab = PyString_AsString (label);
    else lab = NULL;

    IFIS (kind, "FINITE_ELEMENT")
    {
      if (mesh) TYPETEST (is_shape_convex (shape, kwl[2]));
    }

    SHAPE *shp = create_shape (shape, 1);

    for (SHAPE *x = shp; x; x = x->next)
      if (x->kind == SHAPE_CONVEX) CONVEX_Compute_Adjacency (x->data); /* deformable convex juxtaposition needs adjacency information;
								          it is also needed for sparsification in this and other cases */
    IFIS (kind, "RIGID")
    {
      self->bod = BODY_Create (RIG, shp, get_bulk_material (solfec->sol, material), lab, 0, form, NULL, NULL, NULL);
    }
    ELIF (kind, "PSEUDO_RIGID")
    {
      self->bod = BODY_Create (PRB, shp, get_bulk_material (solfec->sol, material), lab, 0, form, NULL, NULL, NULL);
    }
    ELIF (kind, "FINITE_ELEMENT")
    {
      if (mesh)
      {
        msh = mesh->msh;
	mesh->msh = NULL; /* empty */
      }
      else
      {
        TYPETEST (is_mesh (shape, kwl[2]));
      }

      if (modal)
      {
	PyObject *vlist = PyTuple_GetItem (modal, 0),
		 *Elist = PyTuple_GetItem (modal, 1);

        TYPETEST (is_list (vlist, "modal[0]", 0, 0) && is_list (Elist, "modal[1]", 0, 0));

	int n = PyList_Size (vlist),
	    m = PyList_Size (Elist) / n,
	    OK, i;

	if (msh)
	{
	  OK = (m == msh->nodes_count * 3);
	}
	else
	{
	  MESH *tmp = (MESH*)shp->data;
	  OK = (m == tmp->nodes_count * 3);
	}

	if (!OK)
	{
	  PyErr_SetString (PyExc_ValueError, "Modal analysis data size and mesh size do not match");
	  return NULL;
	}

	ERRMEM (val = malloc (sizeof (double [n])));
	E = MX_Create (MXDENSE, m, n, NULL, NULL);

	for (i = 0; i < n; i ++) val [i] = PyFloat_AsDouble (PyList_GetItem (vlist, i));
	for (i = 0; i < n*m; i ++) E->x [i] = PyFloat_AsDouble (PyList_GetItem (Elist, i));
      }

      if (formulation)
      {
        IFIS (formulation, "TL")
	{
	  form = TOTAL_LAGRANGIAN;
	}
	ELIF (formulation, "BC")
	{
	  form = BODY_COROTATIONAL;
	}
	ELIF (formulation, "RO")
	{
	  form = REDUCED_ORDER;

	  if (!(E && val))
	  {
	    PyErr_SetString (PyExc_ValueError, "Modal data must be passed for RO formulation");
	    return NULL;
	  }
	}
	ELSE
	{
	  PyErr_SetString (PyExc_ValueError, "Invalid FEM formulation");
	  return NULL;
	}
      }

      self->bod = BODY_Create (FEM, shp, get_bulk_material (solfec->sol, material), lab, 0, form, msh, E, val);
    }
    ELIF (kind, "OBSTACLE")
    {
      self->bod = BODY_Create (OBS, shp, get_bulk_material (solfec->sol, material), lab, 0, form, NULL, NULL, NULL);
    }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Unknown BODY kind");
      return NULL;
    }

#if MPI
    self->dom = solfec->sol->dom;
    self->id = self->dom->bid; /* before inserting, the body may be deleted in LOCAL_BODIES mode */
#endif

    if ((solfec->sol->dom->dynamic == 0 && self->bod->kind != RIG) || /* XXX => LIM is closest to the quasi-static time stepping; some code parts test body->scheme without checking for quasi-statics */
	 self->bod->form == REDUCED_ORDER) self->bod->scheme = SCH_DEF_LIM; /* reduced order model uses only the 'DEF_LIM' scheme (no advantage in using 'DEF_EXP' since bod->M is dense anyway) */

    DOM_Insert_Body (solfec->sol->dom, self->bod); /* insert body into the domain */
  }

  return (PyObject*)self;
}

/* body wrapper constructor */
static PyObject* lng_BODY_WRAPPER (BODY *bod)
{
  lng_BODY *self;

  self = (lng_BODY*)lng_BODY_TYPE.tp_alloc (&lng_BODY_TYPE, 0);

#if MPI
  self->id = bod->id;
  self->dom = bod->dom;
#endif

  self->bod = bod;

  return (PyObject*)self;
}

/* destructor */
static void lng_BODY_dealloc (lng_BODY *self)
{
  self->ob_type->tp_free ((PyObject*)self);
}

/* getsets */

static PyObject* lng_BODY_get_id (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  return PyInt_FromLong (self->bod->id);

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_id (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_kind (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  return PyString_FromString (BODY_Kind (self->bod));

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_kind (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_label (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  if (self->bod->label)
    return PyString_FromString (self->bod->label);
  else return PyString_FromString (""); /* an empty label */

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_label (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_conf (lng_BODY *self, void *closure)
{
  PyObject *conf;
  int i, size;
  double *q;

#if MPI
  if (IS_HERE (self))
  {
#endif

  size = BODY_Conf_Size (self->bod);
  q = self->bod->conf;

  if ((conf = PyTuple_New (size)))
  {
    for (i = 0; i < size; i ++)
      PyTuple_SetItem (conf, i, PyFloat_FromDouble (q [i]));
  }

  return conf;

#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_conf (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_velo (lng_BODY *self, void *closure)
{
  PyObject *velo;
  int i, size;
  double *u;

#if MPI
  if (IS_HERE (self))
  {
#endif

  size = self->bod->dofs;
  u = self->bod->velo;

  if ((velo = PyTuple_New (size)))
  {
    for (i = 0; i < size; i ++)
      PyTuple_SetItem (velo, i, PyFloat_FromDouble (u [i]));
  }

  return velo;

#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_velo (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_mass (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  return PyFloat_FromDouble (self->bod->ref_mass);

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_mass (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_volume (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  return PyFloat_FromDouble (self->bod->ref_volume);

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_volume (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_center (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  double *c = self->bod->ref_center;
  return Py_BuildValue ("(d, d, d)", c[0], c[1], c[2]);

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_center (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_tensor (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  double *t = self->bod->ref_tensor;
  return Py_BuildValue ("(d, d, d, d, d, d, d, d, d)", t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8]);

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_tensor (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_selfcontact (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  if (self->bod->flags & BODY_DETECT_SELF_CONTACT)
    return PyString_FromString ("ON");
  else return PyString_FromString ("OFF");

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_selfcontact (lng_BODY *self, PyObject *value, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  if (!is_string (value, "selfcontact")) return -1;
  else if (self->bod->kind != FEM)
  {
    PyErr_SetString (PyExc_ValueError, "Self-contact is only valid for FINITE_ELEMENT bodies");
    return -1;
  }

  IFIS (value, "ON") self->bod->flags |= BODY_DETECT_SELF_CONTACT;
  ELIF (value, "OFF") self->bod->flags &= ~BODY_DETECT_SELF_CONTACT;
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid self-contact flag");
    return -1;
  }

#if MPI && LOCAL_BODIES
  }
#endif

  return 0;
}

static PyObject* lng_BODY_get_scheme (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  switch (self->bod->scheme)
  {
  case SCH_RIG_NEG: return PyString_FromString ("RIG_NEG");
  case SCH_RIG_POS: return PyString_FromString ("RIG_POS");
  case SCH_RIG_IMP: return PyString_FromString ("RIG_IMP");
  case SCH_DEF_EXP: return PyString_FromString ("DEF_EXP");
  case SCH_DEF_LIM: return PyString_FromString ("DEF_LIM");
  }

  return NULL;

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_scheme (lng_BODY *self, PyObject *value, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  if (!is_string (value, "scheme")) return -1;

  if (self->bod->dom->dynamic == 0)
  {
    PyErr_Warn (NULL, "Unable to set integration scheme for quasi-statics");
    return 0;
  }

  IFIS (value, "DEFAULT") self->bod->scheme = self->bod->kind == RIG ? SCH_RIG_NEG : SCH_DEF_EXP;
  ELIF (value, "RIG_POS")
  {
    if (self->bod->kind != RIG)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid integration scheme");
      return -1;
    }

    self->bod->scheme = SCH_RIG_POS;
  }
  ELIF (value, "RIG_NEG")
  {
    if (self->bod->kind != RIG)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid integration scheme");
      return -1;
    }

    self->bod->scheme = SCH_RIG_NEG;
  }
  ELIF (value, "RIG_IMP")
  {
    if (self->bod->kind != RIG)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid integration scheme");
      return -1;
    }

    self->bod->scheme = SCH_RIG_IMP;
  }
  ELIF (value, "DEF_EXP")
  {
    if (self->bod->kind == RIG || self->bod->form == REDUCED_ORDER)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid integration scheme");
      return -1;
    }

    self->bod->scheme = SCH_DEF_EXP;
  }
  ELIF (value, "DEF_LIM")
  {
    if (self->bod->kind == RIG)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid integration scheme");
      return -1;
    }

    self->bod->scheme = SCH_DEF_LIM;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid integration scheme");
    return -1;
  }

#if MPI && LOCAL_BODIES
  }
#endif

  return 0;
}

static PyObject* lng_BODY_get_damping (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  return PyFloat_FromDouble (self->bod->damping);

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_damping (lng_BODY *self, PyObject *value, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  if (!is_number_ge (value, "damping", 0.0)) return -1;

  self->bod->damping = PyFloat_AsDouble (value);  

#if MPI && LOCAL_BODIES
  }
#endif

  return 0;
}

static PyObject* lng_BODY_get_constraints (lng_BODY *self, void *closure)
{
#if MPI
  if (IS_HERE (self))
  {
#endif
  PyObject *list, *obj;
  SET *item;
  CON *con;
  int n;

  if (!(list = PyList_New (SET_Size (self->bod->con)))) return NULL;

  for (n = 0, item = SET_First (self->bod->con), con = item->data;
       item; n ++, item = SET_Next (item), con = item->data)
  {
    if (!(obj = lng_CONSTRAINT_WRAPPER (con))) return NULL;

    PyList_SetItem (list, n, obj);
  }

  return list;
#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_constraints (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_ncon (lng_BODY *self, void *closure)
{
#if MPI
  if (IS_HERE (self))
  {
#endif
  return PyInt_FromLong (SET_Size (self->bod->con));
#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_ncon (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_material (lng_BODY *self, void *closure)
{
#if MPI
  if (IS_HERE (self))
  {
#endif
  return lng_BULK_MATERIAL_WRAPPER (self->bod->mat);
#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_material (lng_BODY *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_BODY_get_fracturecheck (lng_BODY *self, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  if (self->bod->flags & BODY_CHECK_FRACTURE)
    return PyString_FromString ("ON");
  else return PyString_FromString ("OFF");

#if MPI && LOCAL_BODIES
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_fracturecheck (lng_BODY *self, PyObject *value, void *closure)
{
#if MPI && LOCAL_BODIES
  if (IS_HERE (self))
  {
#endif

  if (!is_string (value, "scheme")) return -1;

  if (self->bod->kind != FEM)
  {
    PyErr_Warn (NULL, "Unable to check fracture criterion for non FEM bodies!");
    return 0;
  }

  IFIS (value, "ON") self->bod->flags |= BODY_CHECK_FRACTURE;
  ELIF (value, "OFF") self->bod->flags &= ~BODY_CHECK_FRACTURE;
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Neither 'ON' nor 'OFF'!");
    return -1;
  }

#if MPI && LOCAL_BODIES
  }
#endif

  return 0;
}

/* BODY methods */
static PyMethodDef lng_BODY_methods [] =
{ {NULL, NULL, 0, NULL} };

/* BODY members */
static PyMemberDef lng_BODY_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* BODY getset */
static PyGetSetDef lng_BODY_getset [] =
{
  {"id", (getter)lng_BODY_get_id, (setter)lng_BODY_set_id, "id", NULL},
  {"kind", (getter)lng_BODY_get_kind, (setter)lng_BODY_set_kind, "kind", NULL},
  {"label", (getter)lng_BODY_get_label, (setter)lng_BODY_set_label, "label", NULL},
  {"conf", (getter)lng_BODY_get_conf, (setter)lng_BODY_set_conf, "configuration", NULL},
  {"velo", (getter)lng_BODY_get_velo, (setter)lng_BODY_set_velo, "velocity", NULL},
  {"mass", (getter)lng_BODY_get_mass, (setter)lng_BODY_set_mass, "referential mass", NULL},
  {"volume", (getter)lng_BODY_get_volume, (setter)lng_BODY_set_volume, "referential volume", NULL},
  {"center", (getter)lng_BODY_get_center, (setter)lng_BODY_set_center, "referential mass center", NULL},
  {"tensor", (getter)lng_BODY_get_tensor, (setter)lng_BODY_set_tensor, "referential Euler/inertia tensor", NULL},
  {"selfcontact", (getter)lng_BODY_get_selfcontact, (setter)lng_BODY_set_selfcontact, "selfcontact", NULL},
  {"scheme", (getter)lng_BODY_get_scheme, (setter)lng_BODY_set_scheme, "scheme", NULL},
  {"damping", (getter)lng_BODY_get_damping, (setter)lng_BODY_set_damping, "damping", NULL},
  {"constraints", (getter)lng_BODY_get_constraints, (setter)lng_BODY_set_constraints, "constraints list", NULL},
  {"ncon", (getter)lng_BODY_get_ncon, (setter)lng_BODY_set_ncon, "constraints count", NULL},
  {"material", (getter)lng_BODY_get_material, (setter)lng_BODY_set_material, "global body material", NULL},
  {"fracturecheck", (getter)lng_BODY_get_fracturecheck, (setter)lng_BODY_set_fracturecheck, "fracture check", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * TIME_SERIES => object
 */

typedef struct lng_TIME_SERIES lng_TIME_SERIES;

static PyTypeObject lng_TIME_SERIES_TYPE;

struct lng_TIME_SERIES
{
  PyObject_HEAD
  TMS *ts;
};

/* test whether an object is of TIME_SERIES type */
static int is_time_series (lng_TIME_SERIES *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance ((PyObject*)obj, (PyObject*)&lng_TIME_SERIES_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a TIME_SERIES object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether an object is a number or a TIME_SERIES or a callable object */
static int is_number_or_time_series_or_callable (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyNumber_Check (obj) &&
        !PyObject_IsInstance (obj, (PyObject*)&lng_TIME_SERIES_TYPE) &&
	!PyCallable_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a number or a TIME_SERIES object or a callback routine", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}


/* test whether an object is a number or a TIME_SERIES */
static int is_number_or_time_series (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyNumber_Check (obj) &&
        !PyObject_IsInstance (obj, (PyObject*)&lng_TIME_SERIES_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a number or a TIME_SERIES object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* constructor */
static PyObject* lng_TIME_SERIES_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("points");
  lng_TIME_SERIES *self;
  PyObject *points;

  self = (lng_TIME_SERIES*)type->tp_alloc (type, 0);

  if (self)
  {
    PARSEKEYS ("O", &points);

    TYPETEST (is_list_or_string (points, kwl [0], 2, 4));

    if (PyString_Check (points))
    {
      if (!(self->ts = TMS_File (PyString_AsString (points))))
      {
	PyErr_SetString (PyExc_ValueError, "Could not open file");
	return NULL;
      }
    }
    else
    {
      double *times,
	     *values;
      int i, n;

      n = PyList_Size (points) / 2;

      ERRMEM (times = malloc (sizeof (double [n])));
      ERRMEM (values = malloc (sizeof (double [n])));

      for (i = 0; i < n; i ++)
      {
	times [i] = PyFloat_AsDouble (PyList_GetItem (points, 2*i));
	values [i] = PyFloat_AsDouble (PyList_GetItem (points, 2*i + 1));
      }

      self->ts = TMS_Create (n, times, values);
      free (times);
      free (values);
    }
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_TIME_SERIES_dealloc (lng_TIME_SERIES *self)
{
  TMS_Destroy (self->ts);

  self->ob_type->tp_free ((PyObject*)self);
}

/* TIME_SERIES methods */
static PyMethodDef lng_TIME_SERIES_methods [] =
{ {NULL, NULL, 0, NULL} };

/* TIME_SERIES members */
static PyMemberDef lng_TIME_SERIES_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* TIME_SERIES getset */
static PyGetSetDef lng_TIME_SERIES_getset [] =
{ {NULL, 0, 0, NULL, NULL} };

/*
 * GAUSS_SEIDEL_SOLVER => object
 */

typedef struct lng_GAUSS_SEIDEL_SOLVER lng_GAUSS_SEIDEL_SOLVER;

static PyTypeObject lng_GAUSS_SEIDEL_SOLVER_TYPE;

struct lng_GAUSS_SEIDEL_SOLVER
{
  PyObject_HEAD
  GAUSS_SEIDEL *gs;
  PyObject *data;
  PyObject *callback;
};

/* failure callback */
static void lng_GAUSS_SEIDEL_callback (lng_GAUSS_SEIDEL_SOLVER *sol)
{
  PyObject *result;
  PyObject *args;

  if (sol->data)
  {
    if (PyTuple_Check (sol->data))
    {
      if (!(args = PyTuple_New (PyTuple_Size (sol->data) + 1))) goto err;

      PyTuple_SetItem (args, 0, (PyObject*)sol);

      for (int n = 0; n < PyTuple_Size (sol->data); n ++)
	PyTuple_SetItem (args, n+1, PyTuple_GetItem (sol->data, n));
    }
    else args = Py_BuildValue ("(O, O)", (PyObject*)sol, sol->data);
  }
  else args = Py_BuildValue ("(O)", (PyObject*)sol);

  result = PyObject_CallObject (sol->callback, args); /* call user callback */

  Py_DECREF (args);

  if (result)
  {
    if (PyInt_AsLong (result) == 0) /* exit */
    {
      fprintf (stderr, "GAUSS_SEIDEL_SOLVER failed with error code %s\n", GAUSS_SEIDEL_Error (sol->gs));
#if MPI
      MPI_Abort (MPI_COMM_WORLD, 1000);
#else
      lngfinalize ();
#endif
      exit (1);
    }

    Py_DECREF (result);
  }
  else /* error during the Python callback run */
  {
err:
    PyErr_Print (); /* print traceback */
#if MPI
    MPI_Abort (MPI_COMM_WORLD, 1001);
#else
    lngfinalize ();
#endif
    exit (1);
  }
}

/* constructor */
static PyObject* lng_GAUSS_SEIDEL_SOLVER_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("epsilon", "maxiter", "meritval", "failure", "diagepsilon", "diagmaxiter", "diagsolver", "data", "callback");
  double epsilon, meritval, diagepsilon;
  PyObject *failure, *diagsolver;
  lng_GAUSS_SEIDEL_SOLVER *self;
  int maxiter, diagmaxiter;
  GSFAIL gsfail;
  DIAS dias;

  self = (lng_GAUSS_SEIDEL_SOLVER*)type->tp_alloc (type, 0);

  if (self)
  {
    failure = NULL;
    diagepsilon = DBL_MAX;
    diagmaxiter = INT_MAX;
    diagsolver = NULL;
    gsfail = GS_FAILURE_CONTINUE;
    dias = DS_SEMISMOOTH_NEWTON;
    self->data = NULL;
    self->callback = NULL;
    meritval = 1.0;

    PARSEKEYS ("di|dOdiOOO", &epsilon, &maxiter, &meritval, &failure,
      &diagepsilon, &diagmaxiter, &diagsolver, &self->data, &self->callback);

    TYPETEST (is_positive (epsilon, kwl[0]) && is_positive (maxiter, kwl[1]) && is_positive (epsilon, kwl[2]) &&
      is_string (failure, kwl[3]) && is_positive (diagepsilon, kwl[4]) && is_positive (diagmaxiter, kwl[5]) &&
      is_string (diagsolver, kwl[6]) && is_callable (self->callback, kwl[8]));

    if (failure)
    {
      IFIS (failure, "CONTINUE")
      {
	gsfail = GS_FAILURE_CONTINUE;
      }
      ELIF (failure, "EXIT")
      {
	gsfail = GS_FAILURE_EXIT;
      }
      ELIF (failure, "CALLBACK")
      {
	if (!self->callback)
	{
	  PyErr_SetString (PyExc_ValueError, "A callback routine must be given with CALLBACK failure action");
	  return NULL;
	}

	callback_pair_push (self->data, self->callback);
	gsfail = GS_FAILURE_CALLBACK;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid failure action");
	return NULL;
      }
    }

    if (diagsolver)
    {
      IFIS (diagsolver, "SEMISMOOTH_NEWTON")
      {
	dias = DS_SEMISMOOTH_NEWTON;
      }
      ELIF (diagsolver, "PROJECTED_GRADIENT")
      {
	dias = DS_PROJECTED_GRADIENT;
      }
      ELIF (diagsolver, "DE_SAXCE_FENG")
      {
	dias = DS_DE_SAXCE_FENG;
      }
      ELIF (diagsolver, "PROJECTED_NEWTON")
      {
	dias = DS_PROJECTED_NEWTON;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid diagonal solver");
	return NULL;
      }
    }

    if (diagepsilon == DBL_MAX)
    {
      diagepsilon = 0.01 * MIN (epsilon, MIN (meritval, 1E-4));
      if (diagepsilon < DBL_EPSILON) diagepsilon = DBL_EPSILON;
    }

    if (diagmaxiter == INT_MAX)
      diagmaxiter = MAX (100, maxiter / 100);

    self->gs = GAUSS_SEIDEL_Create (epsilon, maxiter, meritval, gsfail, diagepsilon,
      diagmaxiter, dias, self, (GAUSS_SEIDEL_Callback)lng_GAUSS_SEIDEL_callback);
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_GAUSS_SEIDEL_SOLVER_dealloc (lng_GAUSS_SEIDEL_SOLVER *self)
{
#if OPENGL
  if (RND_Is_On ()) return; /* do not delete in viewer mode */
  else
#endif
  GAUSS_SEIDEL_Destroy (self->gs);

  Py_XDECREF (self->callback);
  Py_XDECREF (self->data);

  self->ob_type->tp_free ((PyObject*)self);
}

/* getsets */

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_failure (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyString_FromString (GAUSS_SEIDEL_Failure (self->gs));
}

static int lng_GAUSS_SEIDEL_SOLVER_set_failure (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_error (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyString_FromString (GAUSS_SEIDEL_Error (self->gs));
}

static int lng_GAUSS_SEIDEL_SOLVER_set_error (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_iters (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->gs->iters);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_iters (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_rerhist (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  PyObject *list;
  int i;

  ERRMEM (list = PyList_New (self->gs->iters));

  for (i = 0; i < self->gs->iters; i ++)
    PyList_SetItem (list, i, PyFloat_FromDouble (self->gs->rerhist [i]));

  return list;
}

static int lng_GAUSS_SEIDEL_SOLVER_set_rerhist (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_merhist (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  PyObject *list;
  int i;

  ERRMEM (list = PyList_New (self->gs->iters));

  for (i = 0; i < self->gs->iters; i ++)
    PyList_SetItem (list, i, PyFloat_FromDouble (self->gs->merhist [i]));

  return list;
}

static int lng_GAUSS_SEIDEL_SOLVER_set_merhist (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_epsilon (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->gs->epsilon);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_epsilon (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "epsilon", 0)) return -1;
  self->gs->epsilon = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_maxiter (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->gs->maxiter);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_maxiter (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "maxiter", 0)) return -1;
  self->gs->maxiter = PyInt_AsLong (value);
  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_diagepsilon (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->gs->diagepsilon);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_diagepsilon (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "diagepsilon", 0)) return -1;
  self->gs->diagepsilon = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_diagmaxiter (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->gs->diagmaxiter);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_diagmaxiter (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "diagmaxiter", 0)) return -1;
  self->gs->diagmaxiter = PyInt_AsLong (value);
  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_diagsolver (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyString_FromString (GAUSS_SEIDEL_Diagsolver (self->gs));
}

static int lng_GAUSS_SEIDEL_SOLVER_set_diagsolver (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_string (value, "diagsolver")) return -1;

  IFIS (value, "SEMISMOOTH_NEWTON")
  {
    self->gs->diagsolver = DS_SEMISMOOTH_NEWTON;
  }
  ELIF (value, "PROJECTED_GRADIENT")
  {
    self->gs->diagsolver = DS_PROJECTED_GRADIENT;
  }
  ELIF (value, "DE_SAXCE_FENG")
  {
    self->gs->diagsolver = DS_DE_SAXCE_FENG;
  }
  ELIF (value, "PROJECTED_NEWTON")
  {
    self->gs->diagsolver = DS_PROJECTED_NEWTON;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid diagonal solver");
    return -1;
  }

  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_reverse (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyString_FromString (GAUSS_SEIDEL_Reverse (self->gs));
}

static int lng_GAUSS_SEIDEL_SOLVER_set_reverse (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_string (value, "reverse")) return -1;

  IFIS (value, "ON")
  {
    self->gs->reverse = GS_ON;
  }
  ELIF (value, "OFF")
  {
    self->gs->reverse = GS_OFF;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid reverse switch (ON/OFF accepted)");
    return -1;
  }

  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_variant (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyString_FromString (GAUSS_SEIDEL_Variant (self->gs));
}

static int lng_GAUSS_SEIDEL_SOLVER_set_variant (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_string (value, "variant")) return -1;

  IFIS (value, "FULL")
  {
    self->gs->variant = GS_FULL;
  }
  ELIF (value, "MIDDLE_JACOBI")
  {
    self->gs->variant = GS_MIDDLE_JACOBI;
  }
  ELIF (value, "BOUNDARY_JACOBI")
  {
    self->gs->variant = GS_BOUNDARY_JACOBI;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid variant");
    return -1;
  }

  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_innerloops (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->gs->innerloops);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_innerloops (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "innerloops", 1)) return -1;

  self->gs->innerloops = (int) PyInt_AsLong (value);

  return 0;
}

/* GAUSS_SEIDEL_SOLVER methods */
static PyMethodDef lng_GAUSS_SEIDEL_SOLVER_methods [] =
{ {NULL, NULL, 0, NULL} };

/* GAUSS_SEIDEL_SOLVER members */
static PyMemberDef lng_GAUSS_SEIDEL_SOLVER_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* GAUSS_SEIDEL_SOLVER getset */
static PyGetSetDef lng_GAUSS_SEIDEL_SOLVER_getset [] =
{
  {"failure", (getter)lng_GAUSS_SEIDEL_SOLVER_get_failure, (setter)lng_GAUSS_SEIDEL_SOLVER_set_failure, "failure action", NULL},
  {"error", (getter)lng_GAUSS_SEIDEL_SOLVER_get_error, (setter)lng_GAUSS_SEIDEL_SOLVER_set_error, "error code", NULL},
  {"iters", (getter)lng_GAUSS_SEIDEL_SOLVER_get_iters, (setter)lng_GAUSS_SEIDEL_SOLVER_set_iters, "number of iterations", NULL},
  {"rerhist", (getter)lng_GAUSS_SEIDEL_SOLVER_get_rerhist, (setter)lng_GAUSS_SEIDEL_SOLVER_set_rerhist, "relative error history", NULL},
  {"merhist", (getter)lng_GAUSS_SEIDEL_SOLVER_get_merhist, (setter)lng_GAUSS_SEIDEL_SOLVER_set_merhist, "merit function history", NULL},
  {"epsilon", (getter)lng_GAUSS_SEIDEL_SOLVER_get_epsilon, (setter)lng_GAUSS_SEIDEL_SOLVER_set_epsilon, "relative accuracy", NULL},
  {"maxiter", (getter)lng_GAUSS_SEIDEL_SOLVER_get_maxiter, (setter)lng_GAUSS_SEIDEL_SOLVER_set_maxiter, "iterations bound", NULL},
  {"diagepsilon", (getter)lng_GAUSS_SEIDEL_SOLVER_get_diagepsilon, (setter)lng_GAUSS_SEIDEL_SOLVER_set_diagepsilon, "diagonal solver relative accuracy", NULL},
  {"diagmaxiter", (getter)lng_GAUSS_SEIDEL_SOLVER_get_diagmaxiter, (setter)lng_GAUSS_SEIDEL_SOLVER_set_diagmaxiter, "diagonal solver iterations bound", NULL},
  {"diagsolver", (getter)lng_GAUSS_SEIDEL_SOLVER_get_diagsolver, (setter)lng_GAUSS_SEIDEL_SOLVER_set_diagsolver, "diagonal solver kind", NULL},
  {"reverse", (getter)lng_GAUSS_SEIDEL_SOLVER_get_reverse, (setter)lng_GAUSS_SEIDEL_SOLVER_set_reverse, "iteration reversion flag", NULL},
  {"variant", (getter)lng_GAUSS_SEIDEL_SOLVER_get_variant, (setter)lng_GAUSS_SEIDEL_SOLVER_set_variant, "parallel update variant", NULL},
  {"innerloops", (getter)lng_GAUSS_SEIDEL_SOLVER_get_innerloops, (setter)lng_GAUSS_SEIDEL_SOLVER_set_innerloops, "number of inner loops per one parallel step", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * PENALTY_SOLVER => object
 */

typedef struct lng_PENALTY_SOLVER lng_PENALTY_SOLVER;

static PyTypeObject lng_PENALTY_SOLVER_TYPE;

struct lng_PENALTY_SOLVER
{
  PyObject_HEAD

  PENALTY *ps;
};

/* constructor */
static PyObject* lng_PENALTY_SOLVER_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("variant");
  lng_PENALTY_SOLVER *self;
  PyObject *variant;
  short implicit;

  self = (lng_PENALTY_SOLVER*)type->tp_alloc (type, 0);

  if (self)
  {
    variant = NULL;
    implicit = 1;

    PARSEKEYS ("|O", &variant);

    TYPETEST (is_string (variant, kwl[0]));

    if (variant)
    {
      IFIS (variant, "IMPLICIT")
      {
	implicit = 1;
      }
      ELIF (variant, "EXPLICIT")
      {
	implicit = 0;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid variant");
	return NULL;
      }
    }

    self->ps = PENALTY_Create (implicit);
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_PENALTY_SOLVER_dealloc (lng_PENALTY_SOLVER *self)
{
#if OPENGL
  if (RND_Is_On ()) return; /* do not delete in viewer mode */
  else
#endif
  PENALTY_Destroy (self->ps);

  self->ob_type->tp_free ((PyObject*)self);
}

/* PENALTY_SOLVER methods */
static PyMethodDef lng_PENALTY_SOLVER_methods [] =
{ {NULL, NULL, 0, NULL} };

/* PENALTY_SOLVER members */
static PyMemberDef lng_PENALTY_SOLVER_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* PENALTY_SOLVER getset */
static PyGetSetDef lng_PENALTY_SOLVER_getset [] =
{ {NULL, 0, 0, NULL, NULL} };

/*
 * NEWTON_SOLVER => object
 */

typedef struct lng_NEWTON_SOLVER lng_NEWTON_SOLVER;

static PyTypeObject lng_NEWTON_SOLVER_TYPE;

struct lng_NEWTON_SOLVER
{
  PyObject_HEAD

  NEWTON *ns;
};

/* constructor */
static PyObject* lng_NEWTON_SOLVER_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("meritval", "maxiter", "locdyn", "linver", "linmaxiter", "maxmatvec", "epsilon", "delta", "theta", "omega", "gsflag");
  double meritval, epsilon, delta, theta, omega;
  int maxiter, linmaxiter, maxmatvec;
  PyObject *locdyn, *linver, *gsflag;
  lng_NEWTON_SOLVER *self;

  self = (lng_NEWTON_SOLVER*)type->tp_alloc (type, 0);

  if (self)
  {
    meritval = 1E-8;
    maxiter = 1000;
    locdyn = NULL;
    linver = NULL;
    linmaxiter = 10;
    maxmatvec = 10000;
    epsilon = 0.25;
    delta = 0.0;
    theta = 0.25;
    omega = 1E-10;
    gsflag = NULL;

    PARSEKEYS ("|diOOiiddddO", &meritval, &maxiter, &locdyn, &linver, &linmaxiter, &maxmatvec, &epsilon, &delta, &theta, &omega, &gsflag);

    TYPETEST (is_positive (meritval, kwl[0]) && is_positive (maxiter, kwl[1]) && is_string (locdyn, kwl[2]) && is_string (locdyn, kwl[3]) &&
      is_positive (linmaxiter, kwl[4]) && is_positive (maxmatvec, kwl[5]) && is_positive (epsilon, kwl[6]) && is_non_negative (delta, kwl[7]) &&
      is_gt_le (theta, kwl[8], 0, 1.0) && is_positive (omega, kwl[9]));

    self->ns = NEWTON_Create (meritval, maxiter);
    self->ns->linmaxiter = linmaxiter;
    self->ns->maxmatvec = maxmatvec;
    self->ns->epsilon = epsilon;
    self->ns->delta = delta;
    self->ns->theta = theta;
    self->ns->omega = omega;

    if (locdyn)
    {
      IFIS (locdyn, "ON")
      {
	self->ns->locdyn = LOCDYN_ON;
      }
      ELIF (locdyn, "OFF")
      {
	self->ns->locdyn = LOCDYN_OFF;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid locdyn value: neither ON nor OFF");
	return NULL;
      }
    }

    if (linver)
    {
      IFIS (linver, "GMRES")
      {
	self->ns->linver = PQN_GMRES;
      }
      ELIF (linver, "DIAG")
      {
	self->ns->linver = PQN_DIAG;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid linver value: neither GMRES nor DIAG");
	return NULL;
      }
    }

    if (gsflag)
    {
      IFIS (gsflag, "ON")
      {
	self->ns->gsflag = GS_ON;
      }
      ELIF (gsflag, "OFF")
      {
	self->ns->gsflag = GS_OFF;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid gsflag value: neither ON nor OFF");
	return NULL;
      }
    }
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_NEWTON_SOLVER_dealloc (lng_NEWTON_SOLVER *self)
{
#if OPENGL
  if (RND_Is_On ()) return; /* do not delete in viewer mode */
  else
#endif
  NEWTON_Destroy (self->ns);

  self->ob_type->tp_free ((PyObject*)self);
}

static PyObject* lng_NEWTON_SOLVER_get_meritval (lng_NEWTON_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->ns->meritval);
}

static int lng_NEWTON_SOLVER_set_meritval (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "meritval", 0)) return -1;
  self->ns->meritval = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_maxiter (lng_NEWTON_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->ns->maxiter);
}

static int lng_NEWTON_SOLVER_set_maxiter (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "maxiter", 0)) return -1;
  self->ns->maxiter = PyInt_AsLong (value);
  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_locdyn (lng_NEWTON_SOLVER *self, void *closure)
{
  if (self->ns->locdyn == LOCDYN_ON) return PyString_FromString ("ON");
  else return PyString_FromString ("OFF");
}

static int lng_NEWTON_SOLVER_set_locdyn (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_string (value, "locdyn")) return -1;

  IFIS (value, "ON")
  {
    self->ns->locdyn = LOCDYN_ON;
  }
  ELIF (value, "OFF")
  {
    self->ns->locdyn = LOCDYN_OFF;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid locdyn value: neither ON nor OFF");
    return -1;
  }

  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_linver (lng_NEWTON_SOLVER *self, void *closure)
{
  if (self->ns->linver == PQN_GMRES) return PyString_FromString ("GMRES");
  else return PyString_FromString ("DIAG");
}

static int lng_NEWTON_SOLVER_set_linver (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_string (value, "linver")) return -1;

  IFIS (value, "GMRES")
  {
    self->ns->linver = PQN_GMRES;
  }
  ELIF (value, "DIAG")
  {
    self->ns->linver = PQN_DIAG;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid linver value: neither GMRES nor DIAG");
    return -1;
  }

  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_linmaxiter (lng_NEWTON_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->ns->linmaxiter);
}

static int lng_NEWTON_SOLVER_set_linmaxiter (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "linmaxiter", 0)) return -1;
  self->ns->linmaxiter = PyInt_AsLong (value);
  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_maxmatvec (lng_NEWTON_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->ns->maxmatvec);
}

static int lng_NEWTON_SOLVER_set_maxmatvec (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "maxmatvec", 0)) return -1;
  self->ns->maxmatvec = PyInt_AsLong (value);
  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_epsilon (lng_NEWTON_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->ns->epsilon);
}

static int lng_NEWTON_SOLVER_set_epsilon (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "epsilon", 0)) return -1;
  self->ns->epsilon = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_delta (lng_NEWTON_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->ns->delta);
}

static int lng_NEWTON_SOLVER_set_delta (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "delta", 0)) return -1;
  self->ns->delta  = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_theta (lng_NEWTON_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->ns->theta);
}

static int lng_NEWTON_SOLVER_set_theta (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt_le (value, "theta", 0, 0.5)) return -1;
  self->ns->theta = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_omega (lng_NEWTON_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->ns->omega);
}

static int lng_NEWTON_SOLVER_set_omega (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "omega", 0)) return -1;
  self->ns->omega = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_NEWTON_SOLVER_get_merhist (lng_NEWTON_SOLVER *self, void *closure)
{
  PyObject *list;
  int i;

  ERRMEM (list = PyList_New (self->ns->iters));

  for (i = 0; i < self->ns->iters; i ++)
    PyList_SetItem (list, i, PyFloat_FromDouble (self->ns->merhist [i]));

  return list;
}

static int lng_NEWTON_SOLVER_set_merhist (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_NEWTON_SOLVER_get_mvhist (lng_NEWTON_SOLVER *self, void *closure)
{
  PyObject *list;
  int i;

  ERRMEM (list = PyList_New (self->ns->iters));

  for (i = 0; i < self->ns->iters; i ++)
    PyList_SetItem (list, i, PyFloat_FromDouble (self->ns->mvhist [i]));

  return list;
}

static int lng_NEWTON_SOLVER_set_mvhist (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_NEWTON_SOLVER_get_iters (lng_NEWTON_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->ns->iters);
}

static int lng_NEWTON_SOLVER_set_iters (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_NEWTON_SOLVER_get_gsflag (lng_NEWTON_SOLVER *self, void *closure)
{
  if (self->ns->gsflag == GS_ON) return PyString_FromString ("ON");
  else return PyString_FromString ("OFF");
}

static int lng_NEWTON_SOLVER_set_gsflag (lng_NEWTON_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_string (value, "locdyn")) return -1;

  IFIS (value, "ON")
  {
    self->ns->gsflag = GS_ON;
  }
  ELIF (value, "OFF")
  {
    self->ns->gsflag = GS_OFF;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid locdyn value: neither ON nor OFF");
    return -1;
  }

  return 0;
}

/* NEWTON_SOLVER methods */
static PyMethodDef lng_NEWTON_SOLVER_methods [] =
{ {NULL, NULL, 0, NULL} };

/* NEWTON_SOLVER members */
static PyMemberDef lng_NEWTON_SOLVER_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* NEWTON_SOLVER getset */
static PyGetSetDef lng_NEWTON_SOLVER_getset [] =
{ 
  {"meritval", (getter)lng_NEWTON_SOLVER_get_meritval, (setter)lng_NEWTON_SOLVER_set_meritval, "merit function accuracy", NULL},
  {"maxiter", (getter)lng_NEWTON_SOLVER_get_maxiter, (setter)lng_NEWTON_SOLVER_set_maxiter, "iterations bound", NULL},
  {"locdyn", (getter)lng_NEWTON_SOLVER_get_locdyn, (setter)lng_NEWTON_SOLVER_set_locdyn, "local dynamics assembling", NULL},
  {"linver", (getter)lng_NEWTON_SOLVER_get_linver, (setter)lng_NEWTON_SOLVER_set_linver, "linearization version", NULL},
  {"linmaxiter", (getter)lng_NEWTON_SOLVER_get_linmaxiter, (setter)lng_NEWTON_SOLVER_set_linmaxiter, "GMRES iterations bound", NULL},
  {"maxmatvec", (getter)lng_NEWTON_SOLVER_get_maxmatvec, (setter)lng_NEWTON_SOLVER_set_maxmatvec, "GMRES matrix-vector products bound", NULL},
  {"epsilon", (getter)lng_NEWTON_SOLVER_get_epsilon, (setter)lng_NEWTON_SOLVER_set_epsilon, "GMRES relative accuracy", NULL},
  {"delta", (getter)lng_NEWTON_SOLVER_get_delta, (setter)lng_NEWTON_SOLVER_set_delta, "diagonal regularization", NULL},
  {"theta", (getter)lng_NEWTON_SOLVER_get_theta, (setter)lng_NEWTON_SOLVER_set_theta, "relaxation parameter", NULL},
  {"omega", (getter)lng_NEWTON_SOLVER_get_omega, (setter)lng_NEWTON_SOLVER_set_omega, "equation smoothing parameter", NULL},
  {"merhist", (getter)lng_NEWTON_SOLVER_get_merhist, (setter)lng_NEWTON_SOLVER_set_merhist, "merit function history", NULL},
  {"mvhist", (getter)lng_NEWTON_SOLVER_get_mvhist, (setter)lng_NEWTON_SOLVER_set_mvhist, "matrix-vector products history", NULL},
  {"iters", (getter)lng_NEWTON_SOLVER_get_iters, (setter)lng_NEWTON_SOLVER_set_iters, "iterations count", NULL},
  {"gsflag", (getter)lng_NEWTON_SOLVER_get_gsflag, (setter)lng_NEWTON_SOLVER_set_gsflag, "Gauss-Seidel failure iterations flag", NULL},
  {NULL, 0, 0, NULL, NULL}
};

#if WITHSICONOS

/*
 * SICONOS_SOLVER => object
 */

typedef struct lng_SICONOS_SOLVER lng_SICONOS_SOLVER;

static PyTypeObject lng_SICONOS_SOLVER_TYPE;

struct lng_SICONOS_SOLVER
{
  PyObject_HEAD

  SICONOS *si;
};

/* constructor */
static PyObject* lng_SICONOS_SOLVER_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("epsilon", "maxiter", "verbose");
  lng_SICONOS_SOLVER *self;
  PyObject *verbose;
  double epsilon;
  int maxiter;

  self = (lng_SICONOS_SOLVER*)type->tp_alloc (type, 0);

  if (self)
  {
    epsilon = 1E-4;
    maxiter = 1000;
    verbose = NULL;

    PARSEKEYS ("|diO", &epsilon, &maxiter, &verbose);

    TYPETEST (is_positive (epsilon, kwl[0]) && is_positive (maxiter, kwl[1]) && is_string (verbose, kwl[2]));

    self->si = SICONOS_Create (epsilon, maxiter);

    if (verbose)
    {
      IFIS (verbose, "ON")
      {
	self->si->verbose = 1;
      }
      ELIF (verbose, "OFF")
      {
	self->si->verbose = 0;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid verbose value ('ON' or 'OFF' are valid only)");
	return NULL;
      }
    }
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_SICONOS_SOLVER_dealloc (lng_SICONOS_SOLVER *self)
{
#if OPENGL
  if (RND_Is_On ()) return; /* do not delete in viewer mode */
  else
#endif
  SICONOS_Destroy (self->si);

  self->ob_type->tp_free ((PyObject*)self);
}

static PyObject* lng_SICONOS_SOLVER_get_epsilon (lng_SICONOS_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->si->epsilon);
}

static int lng_SICONOS_SOLVER_set_epsilon (lng_SICONOS_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "epsilon", 0)) return -1;
  self->si->epsilon = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_SICONOS_SOLVER_get_maxiter (lng_SICONOS_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->si->maxiter);
}

static int lng_SICONOS_SOLVER_set_maxiter (lng_SICONOS_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "maxiter", 0)) return -1;
  self->si->maxiter = PyInt_AsLong (value);
  return 0;
}

/* SICONOS_SOLVER methods */
static PyMethodDef lng_SICONOS_SOLVER_methods [] =
{ {NULL, NULL, 0, NULL} };

/* SICONOS_SOLVER members */
static PyMemberDef lng_SICONOS_SOLVER_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* SICONOS_SOLVER getset */
static PyGetSetDef lng_SICONOS_SOLVER_getset [] =
{ 
  {"epsilon", (getter)lng_SICONOS_SOLVER_get_epsilon, (setter)lng_SICONOS_SOLVER_set_epsilon, "relative reaction change bound", NULL},
  {"maxiter", (getter)lng_SICONOS_SOLVER_get_maxiter, (setter)lng_SICONOS_SOLVER_set_maxiter, "iterations bound", NULL},
  {NULL, 0, 0, NULL, NULL}
};

#endif

/*
 * TEST_SOLVER => object
 */

typedef struct lng_TEST_SOLVER lng_TEST_SOLVER;

static PyTypeObject lng_TEST_SOLVER_TYPE;

struct lng_TEST_SOLVER
{
  PyObject_HEAD

  TEST *ts;
};

/* constructor */
static PyObject* lng_TEST_SOLVER_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("meritval", "maxiter", "maxmatvec", "linmaxiter", "epsilon", "delta", "omega");
  lng_TEST_SOLVER *self;
  int maxiter, linmaxiter, maxmatvec;
  double meritval, epsilon, delta, omega;

  self = (lng_TEST_SOLVER*)type->tp_alloc (type, 0);

  if (self)
  {
    meritval = 1E-8;
    maxiter = 1000;
    maxmatvec = 1000;
    linmaxiter = 10;
    epsilon = 0.25;
    delta = 0.51E-7;
    omega = 1E-9;

    PARSEKEYS ("|diiiddd", &meritval, &maxiter, &maxmatvec, &linmaxiter, &epsilon, &delta, &omega);

    TYPETEST (is_positive (meritval, kwl[0]) && is_positive (maxiter, kwl[1]) && is_positive (maxiter, kwl[2]) &&
	      is_positive (maxiter, kwl[3]) && is_positive (epsilon, kwl[4]) && is_positive (delta, kwl[5]) &&
	      is_positive (omega, kwl[6]));

    self->ts = TEST_Create (meritval, maxiter);
    self->ts->maxmatvec = maxmatvec;
    self->ts->linmaxiter = linmaxiter;
    self->ts->epsilon = epsilon;
    self->ts->delta = delta;
    self->ts->omega = omega;
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_TEST_SOLVER_dealloc (lng_TEST_SOLVER *self)
{
#if OPENGL
  if (RND_Is_On ()) return; /* do not delete in viewer mode */
  else
#endif
  TEST_Destroy (self->ts);

  self->ob_type->tp_free ((PyObject*)self);
}

static PyObject* lng_TEST_SOLVER_get_meritval (lng_TEST_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->ts->meritval);
}

static int lng_TEST_SOLVER_set_meritval (lng_TEST_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "meritval", 0)) return -1;
  self->ts->meritval = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_TEST_SOLVER_get_maxiter (lng_TEST_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->ts->maxiter);
}

static int lng_TEST_SOLVER_set_maxiter (lng_TEST_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_gt (value, "maxiter", 0)) return -1;
  self->ts->maxiter = PyInt_AsLong (value);
  return 0;
}

static PyObject* lng_TEST_SOLVER_get_linmaxiter (lng_TEST_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->ts->linmaxiter);
}

static int lng_TEST_SOLVER_set_linmaxiter (lng_TEST_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number_ge (value, "linmaxiter", 1)) return -1;
  self->ts->linmaxiter = PyInt_AsLong (value);
  return 0;
}

static PyObject* lng_TEST_SOLVER_get_merhist (lng_TEST_SOLVER *self, void *closure)
{
  PyObject *list;
  int i;

  ERRMEM (list = PyList_New (self->ts->iters));

  for (i = 0; i < self->ts->iters; i ++)
    PyList_SetItem (list, i, PyFloat_FromDouble (self->ts->merhist [i]));

  return list;
}

static int lng_TEST_SOLVER_set_merhist (lng_TEST_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_TEST_SOLVER_get_mvhist (lng_TEST_SOLVER *self, void *closure)
{
  PyObject *list;
  int i;

  ERRMEM (list = PyList_New (self->ts->iters));

  for (i = 0; i < self->ts->iters; i ++)
    PyList_SetItem (list, i, PyInt_FromLong (self->ts->mvhist [i]));

  return list;
}

static int lng_TEST_SOLVER_set_mvhist (lng_TEST_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_TEST_SOLVER_get_iters (lng_TEST_SOLVER *self, void *closure)
{
  return PyInt_FromLong (self->ts->iters);
}

static int lng_TEST_SOLVER_set_iters (lng_TEST_SOLVER *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* TEST_SOLVER methods */
static PyMethodDef lng_TEST_SOLVER_methods [] =
{ {NULL, NULL, 0, NULL} };

/* TEST_SOLVER members */
static PyMemberDef lng_TEST_SOLVER_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* TEST_SOLVER getset */
static PyGetSetDef lng_TEST_SOLVER_getset [] =
{ 
  {"meritval", (getter)lng_TEST_SOLVER_get_meritval, (setter)lng_TEST_SOLVER_set_meritval, "merit function accuracy", NULL},
  {"maxiter", (getter)lng_TEST_SOLVER_get_maxiter, (setter)lng_TEST_SOLVER_set_maxiter, "iterations bound", NULL},
  {"linmaxiter", (getter)lng_TEST_SOLVER_get_linmaxiter, (setter)lng_TEST_SOLVER_set_linmaxiter, "linear solver maximal iterations count", NULL},
  {"merhist", (getter)lng_TEST_SOLVER_get_merhist, (setter)lng_TEST_SOLVER_set_merhist, "merit function history", NULL},
  {"mvhist", (getter)lng_TEST_SOLVER_get_mvhist, (setter)lng_TEST_SOLVER_set_mvhist, "merit function history", NULL},
  {"iters", (getter)lng_TEST_SOLVER_get_iters, (setter)lng_TEST_SOLVER_set_iters, "iterations count", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * CONSTRAINT => object
 */

typedef struct lng_CONSTRAINT lng_CONSTRAINT; /* constraint type */

static PyTypeObject lng_CONSTRAINT_TYPE;

struct lng_CONSTRAINT
{
  PyObject_HEAD

  unsigned int id;

  DOM *dom;

  CON *con;
};

/* try assigning a constraint pointer to the id (constraints get
 * deleted and recreated in READ mode and they migrate in parallel) */
static int ID_TO_CONSTRAINT (DOM *dom, lng_CONSTRAINT *constraint)
{
  if (dom && constraint->id)
  {
    CON *con = MAP_Find (dom->idc, (void*) (long) constraint->id, NULL);

    if (con)
    {
      constraint->con = con;
      return 1;
    }
  }

  return 0;
}

/* test whether an object is of BODY or CONSTRAINT type */
static int is_body_or_constraint (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance (obj, (PyObject*)&lng_BODY_TYPE) &&
        !PyObject_IsInstance (obj, (PyObject*)&lng_CONSTRAINT_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a BODY or a CONSTRAINT object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
    else if (PyObject_IsInstance (obj, (PyObject*)&lng_BODY_TYPE))
    {
      lng_BODY *body = (lng_BODY*)obj;
      if (body->bod == NULL)
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' is a recently deleted BODY object", var);
	PyErr_SetString (PyExc_TypeError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* constructor */
static PyObject* lng_CONSTRAINT_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  PyErr_SetString (PyExc_RuntimeError, "Only specific CONSTRAINT objects can be created (e.g. FIX_POINT)");
  return NULL;
}

/* constraint wrapper constructor */
static PyObject* lng_CONSTRAINT_WRAPPER (CON *con)
{
  lng_CONSTRAINT *self;

  self = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  self->dom = con->master->dom;
  self->id = con->id;
  self->con = con;

  return (PyObject*)self;
}

/* destructor */
static void lng_CONSTRAINT_dealloc (lng_CONSTRAINT *self)
{
  self->ob_type->tp_free ((PyObject*)self);
}

/* getsets */

static PyObject* lng_CONSTRAINT_get_kind (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    return PyString_FromString (CON_Kind (self->con));
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_kind (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_R (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    double *R = self->con->R;
    return Py_BuildValue ("(d, d, d)", R[0], R[1], R[2]);
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_R (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_U (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    double *U = self->con->U;
    return Py_BuildValue ("(d, d, d)", U[0], U[1], U[2]);
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_U (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_V (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    double *V = self->con->V;
    return Py_BuildValue ("(d, d, d)", V[0], V[1], V[2]);
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_V (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_base (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    double *base = self->con->base;
    return Py_BuildValue ("(d, d, d, d, d, d, d, d, d)", base [0], base [1],
       base [2], base [3], base [4], base [5], base [6], base [7], base [8]);
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_base (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
    return -1;
  }
  else return -1;
}

static PyObject* lng_CONSTRAINT_get_point (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    double *point = self->con->point;

    if (self->con->kind == RIGLNK)
    {
      double other [3];

      ADD (point, RIGLNK_VEC (self->con->Z), other);
      return Py_BuildValue ("(d, d, d, d, d, d)", point[0],
	  point[1], point[2], other[0], other[1], other[2]);
    }

    return Py_BuildValue ("(d, d, d)", point[0], point[1], point[2]);
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_point (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_area (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    return PyFloat_FromDouble (self->con->area);
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_area (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_gap (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    return PyFloat_FromDouble (self->con->gap);
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_gap (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_merit (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    return PyFloat_FromDouble (self->con->merit);
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_merit (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_adjbod (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    BODY *master = self->con->master,
	 *slave = self->con->slave;
    PyObject *m, *s;


    if (master && slave)
    {
      m = lng_BODY_WRAPPER (master);
      s = lng_BODY_WRAPPER (slave);

      if (!m || !s) return NULL;
      else return Py_BuildValue ("(O, O)", m, s);
    }
    else return lng_BODY_WRAPPER (master);
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_adjbod (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_matlab (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    if (self->con->kind == CONTACT)
      return PyString_FromString (self->con->mat.base->label);
    else Py_RETURN_NONE;
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_matlab (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_spair (lng_CONSTRAINT *self, void *closure)
{
  if (ID_TO_CONSTRAINT (self->dom, self))
  {
    CON *con = self->con;

    if (con->kind == CONTACT)
    {
      return Py_BuildValue ("(i, i)", con->spair[0], con->spair[1]);
    }
    else Py_RETURN_NONE; 
  }
  else Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_spair (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

/* CONSTRAINT methods */
static PyMethodDef lng_CONSTRAINT_methods [] =
{ {NULL, NULL, 0, NULL} };

/* CONSTRAINT members */
static PyMemberDef lng_CONSTRAINT_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* CONSTRAINT getset */
static PyGetSetDef lng_CONSTRAINT_getset [] =
{ 
  {"kind", (getter)lng_CONSTRAINT_get_kind, (setter)lng_CONSTRAINT_set_kind, "constraint kind", NULL},
  {"R", (getter)lng_CONSTRAINT_get_R, (setter)lng_CONSTRAINT_set_R, "constraint reaction", NULL},
  {"U", (getter)lng_CONSTRAINT_get_U, (setter)lng_CONSTRAINT_set_U, "constraint output velocity", NULL},
  {"V", (getter)lng_CONSTRAINT_get_V, (setter)lng_CONSTRAINT_set_V, "contact input velocity", NULL},
  {"base", (getter)lng_CONSTRAINT_get_base, (setter)lng_CONSTRAINT_set_base, "constraint local base", NULL},
  {"point", (getter)lng_CONSTRAINT_get_point, (setter)lng_CONSTRAINT_set_point, "constraint spatial point", NULL},
  {"area", (getter)lng_CONSTRAINT_get_area, (setter)lng_CONSTRAINT_set_area, "constraint area", NULL},
  {"gap", (getter)lng_CONSTRAINT_get_gap, (setter)lng_CONSTRAINT_set_gap, "constraint gap", NULL},
  {"merit", (getter)lng_CONSTRAINT_get_merit, (setter)lng_CONSTRAINT_set_merit, "constraint merit function", NULL},
  {"adjbod", (getter)lng_CONSTRAINT_get_adjbod, (setter)lng_CONSTRAINT_set_adjbod, "constraint adjacent bodies", NULL},
  {"matlab", (getter)lng_CONSTRAINT_get_matlab, (setter)lng_CONSTRAINT_set_matlab, "contact constraint material", NULL},
  {"spair", (getter)lng_CONSTRAINT_get_spair, (setter)lng_CONSTRAINT_set_spair, "contact surface pairing", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/* 
 * Solfec => module
 */

/* create convex hull of a point set */
static PyObject* lng_HULL (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("points", "volid", "surfid", "convex");
  int volid, surfid, npnt, i, l;
  lng_CONVEX *out, *convex;
  PyObject *points;
  CONVEX *cvx;
  double *pnt;
  int error;

  out = (lng_CONVEX*)lng_CONVEX_TYPE.tp_alloc (&lng_CONVEX_TYPE, 0);

  if (out)
  {
    convex = NULL;
    cvx = NULL;

    PARSEKEYS ("Oii|O", &points, &volid, &surfid, &convex);

    TYPETEST (is_list (points, kwl[0], 3, 12) && is_convex (convex, kwl[3]));

    l = PyList_Size (points);
    ERRMEM (pnt = malloc (l * sizeof (double)));
    npnt = l / 3;

    for (i = 0; i < l; i ++) pnt [i] = PyFloat_AsDouble (PyList_GetItem (points, i));

    if (convex) cvx = convex->cvx; /* we are appending a list */

    TRY ()
      out->cvx = CONVEX_Hull (cvx, pnt, npnt, surfid, volid);
    CATCHANY (error)
    {
      PyErr_SetString (PyExc_RuntimeError, errstring (error));
#if MPI
      PyErr_Print ();
      MPI_Abort (MPI_COMM_WORLD, 2000+error);
#endif
      return NULL;
    }
    ENDTRY ()

    if (out->cvx == NULL)
    {
      PyErr_SetString (PyExc_RuntimeError, "Failed to create convex hull");
      return NULL;
    }

    if (cvx) convex->cvx = NULL; /* empty */

    free (pnt);
  }

  return (PyObject*)out;
}

/* create convex list out of a MESH object */
static PyObject* lng_MESH2CONVEX (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("mesh");
  lng_MESH *mesh;
  lng_CONVEX *out;

  out = (lng_CONVEX*)lng_CONVEX_TYPE.tp_alloc (&lng_CONVEX_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("O", &mesh);

    TYPETEST (is_mesh ((PyObject*)mesh, kwl[0]));

    if (mesh->msh == NULL)
    {
      PyErr_SetString (PyExc_RuntimeError, "The mesh object is empty");
      return NULL;
    }

    out->cvx = MESH_Convex (mesh->msh, 0, 0);
  }

  return (PyObject*)out;
}

/* create a mesh object */
static PyObject* lng_HEX (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("nodes", "i", "j", "k", "volid", "surids", "dx", "dy", "dz");
  PyObject *nodes, *surfids, *dx, *dy, *dz;
  double (*nod) [3], *vdx, *vdy, *vdz;
  int i, j, k, n, m, surfaces [6], volid;
  lng_MESH *out;

  out = (lng_MESH*)lng_MESH_TYPE.tp_alloc (&lng_MESH_TYPE, 0);

  if (out)
  {
    dx = dy = dz = NULL;

    PARSEKEYS ("OiiiiO|OOO", &nodes, &i, &j, &k, &volid, &surfids, &dx, &dy, &dz);

    TYPETEST (is_list (nodes, kwl[0], 3, 24) && is_list (surfids, kwl[5], 1, 6) &&
	      is_positive (i, kwl[1]) && is_positive (j, kwl[2]) && is_positive (k, kwl[3]) &&
	      is_list (dx, kwl[6], 1, i) && is_list (dy, kwl[7], 1, j) && is_list (dz, kwl[8], 1, k));

    ERRMEM (nod = malloc (8 * sizeof (double [3])));
    for (n = 0; n < 8; n ++)
      for (m = 0; m < 3; m ++)
	nod [n][m] = PyFloat_AsDouble (PyList_GetItem (nodes, 3*n + m));

    for (n = 0; n < 6; n ++)
      surfaces [n] = PyInt_AsLong (PyList_GetItem (surfids, n));

    if (dx)
    {
      ERRMEM (vdx = malloc (i * sizeof (double)));
      for (n = 0; n < i; n ++)
      {
        vdx [n] = PyFloat_AsDouble (PyList_GetItem (dx, n));
	TYPETEST (is_positive (vdx [n], "dx[i] for all i"));
      }
    }
    else vdx = NULL;

    if (dy)
    {
      ERRMEM (vdy = malloc (j * sizeof (double)));
      for (n = 0; n < j; n ++)
      {
        vdy [n] = PyFloat_AsDouble (PyList_GetItem (dy, n));
	TYPETEST (is_positive (vdy [n], "dy[j] for all j"));
      }
    }
    else vdy = NULL;

    if (dz)
    {
      ERRMEM (vdz = malloc (k * sizeof (double)));
      for (n = 0; n < k; n ++)
      {
        vdz [n] = PyFloat_AsDouble (PyList_GetItem (dz, n));
	TYPETEST (is_positive (vdz [n], "dz[k] for all k"));
      }
    }
    else vdz = NULL;

    out->msh = MESH_Hex (nod, i, j, k, surfaces, volid, vdx, vdy, vdz);

    free (nod);
    free (vdx);
    free (vdy);
    free (vdz);
  }

  return (PyObject*)out;
}

/* create tetrahedral mesh */
static PyObject* lng_TETRAHEDRALIZE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "path", "volume", "quality", "volid", "surid");
  double volume, quality;
  PyObject *shape, *path;
  int volid, surfid;
  lng_MESH *out;

  out = (lng_MESH*)lng_MESH_TYPE.tp_alloc (&lng_MESH_TYPE, 0);

  if (out)
  {
    volume = -DBL_MAX;
    quality = -DBL_MAX; 
    volid = -INT_MAX;
    surfid = -INT_MAX;

    PARSEKEYS ("OO|ddii", &shape, &path, &volume, &quality, &volid, &surfid);

    TYPETEST (is_string (path, kwl [1]));

#if MPI /* avoid writing to the same file from several processes  -- begin */
    int counter, rank;

    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    for (counter = 0; counter < 2; counter ++)
    {
    if ((rank == 0 && counter == 0) || /* let the rank 0 process go first, while others skip at first */
	(rank > 0 && counter > 0)) /* at second run, let the rank 0 process skip, while others enter */
    {
#endif

    out->msh = MESH_Read (PyString_AsString (path));

    if (out->msh == NULL)
    {
      if (volume != -DBL_MAX && volume <= 0.0)
      {
	PyErr_SetString (PyExc_ValueError, "Maximal volume must be positive");
	return NULL;
      }
      if (volume == -DBL_MAX) volume = 0.0;

      if (quality != -DBL_MAX && quality <= 1.0)
      {
	PyErr_SetString (PyExc_ValueError, "Quality must be > 1.0");
	return NULL;
      }
      if (quality == -DBL_MAX) quality = 0.0;

      if (PyObject_IsInstance (shape, (PyObject*)&lng_MESH_TYPE))
      {
	out->msh = tetrahedralize1 (((lng_MESH*)shape)->msh, volume, quality, volid, surfid);
	if (!out->msh)
	{
	  PyErr_SetString (PyExc_ValueError, "Mesh generation has failed");
	  return NULL;
	}
      }
      else if (PyString_Check (shape))
      {
	out->msh = tetrahedralize2 (PyString_AsString (shape), volume, quality, volid, surfid);
	if (!out->msh)
	{
	  PyErr_SetString (PyExc_ValueError, "Mesh generation has failed");
	  return NULL;
	}
      }
      else
      {
	PyErr_SetString (PyExc_ValueError, "Shape must be aither a MESH object or string");
	return NULL;
      }

      MESH_Write (out->msh, PyString_AsString (path));
    }

#if MPI /* avoid writing to the same file from several processes  -- end */
    }
    MPI_Barrier (MPI_COMM_WORLD); /* all processes meet here twice */
    }
#endif

  }

  return (PyObject*)out;
}

/* create a mesh object */
static PyObject* lng_PIPE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("pnt", "dir", "rin", "thi", "ndri", "nrad", "nthi", "volid", "surids");
  int ndir, nrad, nthi, volid, surfaces [4];
  double point [3], direct [3], rin, thi;
  PyObject *pnt, *dir, *surfids;
  lng_MESH *out;

  out = (lng_MESH*)lng_MESH_TYPE.tp_alloc (&lng_MESH_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOddiiiiO", &pnt, &dir, &rin, &thi, &ndir, &nrad, &nthi, &volid, &surfids);

    TYPETEST (is_tuple (pnt, kwl[0], 3) && is_tuple (dir, kwl[1], 3) &&
	      is_positive (rin, kwl[2]) && is_positive (thi, kwl[3]) && is_positive (ndir, kwl[4]) &&
	      is_positive (nrad, kwl[5]) && is_positive (nthi, kwl[6]) && is_list (surfids, kwl[8], 1, 4));

    point [0] = PyFloat_AsDouble (PyTuple_GetItem (pnt, 0));
    point [1] = PyFloat_AsDouble (PyTuple_GetItem (pnt, 1));
    point [2] = PyFloat_AsDouble (PyTuple_GetItem (pnt, 2));

    direct [0] = PyFloat_AsDouble (PyTuple_GetItem (dir, 0));
    direct [1] = PyFloat_AsDouble (PyTuple_GetItem (dir, 1));
    direct [2] = PyFloat_AsDouble (PyTuple_GetItem (dir, 2));

    surfaces [0] = PyInt_AsLong (PyList_GetItem (surfids, 0));
    surfaces [1] = PyInt_AsLong (PyList_GetItem (surfids, 1));
    surfaces [2] = PyInt_AsLong (PyList_GetItem (surfids, 2));
    surfaces [3] = PyInt_AsLong (PyList_GetItem (surfids, 3));

    if (!(out->msh = MESH_Pipe (point, direct, rin, thi, ndir, nrad, nthi, surfaces, volid))) return NULL;
  }

  return (PyObject*)out;
}

/* create a rough mesh object */
static PyObject* lng_ROUGH_HEX (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "i", "j", "k", "dx", "dy", "dz");
  PyObject *shape, *dx, *dy, *dz;
  double (*nod) [3], *vdx, *vdy, *vdz;
  int i, j, k, n, surfaces [6], volid;
  lng_MESH *out;
  SHAPE *shp;

  out = (lng_MESH*)lng_MESH_TYPE.tp_alloc (&lng_MESH_TYPE, 0);

  if (out)
  {
    dx = dy = dz = NULL;

    PARSEKEYS ("Oiii|OOO", &shape, &i, &j, &k, &dx, &dy, &dz);

    TYPETEST (is_shape_convex (shape, kwl[0]) && is_positive (i, kwl[1]) && is_positive (j, kwl[2]) && is_positive (k, kwl[3]) &&
	      is_list (dx, kwl[4], 1, i) && is_list (dy, kwl[5], 1, j) && is_list (dz, kwl[6], 1, k));

    MX_DENSE (euler, 3, 3);
    MX_DENSE (eigv, 3, 3);
    double eval [3],
	   extents [6],
	   point [3],
	   *vx = eigv.x,
	   *vy = vx + 3,
	   *vz = vy + 3;
    int o [4] = {0, 1, 2, 3};

    shp = create_shape (shape, 0); /* do not empty */
    SHAPE_Char (shp, 0, NULL, NULL, euler.x);
    MX_Eigen (&euler, 3, eval, &eigv);
    SHAPE_Oriented_Extents (shp, vx, vy, vz, extents);
    SHAPE_Destroy_Wrapper (shp);

    PRODUCT (vx, vy, point);
    if (DOT (point, vz) < 0.0) o[3] = 0, o[2] = 1, o[1] = 2, o[0] = 3; /* keep the orientation right */

    ERRMEM (nod = malloc (8 * sizeof (double [3])));
    COPY (extents, point);
    NVMUL (eigv.x, point, nod [o[0]]);
    point [0] = extents [3];
    NVMUL (eigv.x, point, nod [o[1]]);
    point [1] = extents [4];
    NVMUL (eigv.x, point, nod [o[2]]);
    point [0] = extents [0];
    NVMUL (eigv.x, point, nod [o[3]]);
    point [1] = extents [1];
    point [2] = extents [5];
    NVMUL (eigv.x, point, nod [4+o[0]]);
    point [0] = extents [3];
    NVMUL (eigv.x, point, nod [4+o[1]]);
    point [1] = extents [4];
    NVMUL (eigv.x, point, nod [4+o[2]]);
    point [0] = extents [0];
    NVMUL (eigv.x, point, nod [4+o[3]]);

    for (volid = INT_MAX, n = 0; n < 6; n ++) surfaces [n] = INT_MAX; /* fake ids */

    if (dx)
    {
      ERRMEM (vdx = malloc (i * sizeof (double)));
      for (n = 0; n < i; n ++)
      {
        vdx [n] = PyFloat_AsDouble (PyList_GetItem (dx, n));
	TYPETEST (is_positive (vdx [n], "dx[i] for all i"));
      }
    }
    else vdx = NULL;

    if (dy)
    {
      ERRMEM (vdy = malloc (j * sizeof (double)));
      for (n = 0; n < j; n ++)
      {
        vdy [n] = PyFloat_AsDouble (PyList_GetItem (dy, n));
	TYPETEST (is_positive (vdy [n], "dy[j] for all j"));
      }
    }
    else vdy = NULL;

    if (dz)
    {
      ERRMEM (vdz = malloc (k * sizeof (double)));
      for (n = 0; n < k; n ++)
      {
        vdz [n] = PyFloat_AsDouble (PyList_GetItem (dz, n));
	TYPETEST (is_positive (vdz [n], "dz[k] for all k"));
      }
    }
    else vdz = NULL;

    out->msh = MESH_Hex (nod, i, j, k, surfaces, volid, vdx, vdy, vdz);

    free (nod);
    free (vdx);
    free (vdy);
    free (vdz);
  }

  return (PyObject*)out;
}

/* create fixed point constraint */
static PyObject* lng_FIX_POINT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point", "strength");
  double p [3], strength;
  lng_CONSTRAINT *out;
  PyObject *point;
  lng_BODY *body;

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    strength = DBL_MAX;

    PARSEKEYS ("OO|d", &body, &point, &strength);

    TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3));

#if MPI
    if (IS_HERE (body))
    {
#endif

    out->dom = body->bod->dom;

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    if (!(out->con = DOM_Fix_Point (body->bod->dom, body->bod, p, strength)))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }
    else out->id = out->con->id;

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create fixed direction constraint */
static PyObject* lng_FIX_DIRECTION (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point", "direction", "body2", "point2");
  lng_CONSTRAINT *out;
  lng_BODY *body, *body2;
  PyObject *point, *direction, *point2;
  double p [3], d [3], p2 [3];

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    body2 = NULL;
    point2 = NULL;

    PARSEKEYS ("OOO|OO", &body, &point, &direction, &body2, &point2);

    TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3) && is_tuple (direction, kwl[2], 3) &&
              is_body (body2, kwl[3]) && is_tuple (point2, kwl[4], 3));

#if MPI
    if (IS_HERE (body))
    {
#endif

    out->dom = body->bod->dom;

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
    d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
    d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

    if (body2 && point2)
    {
      p2 [0] = PyFloat_AsDouble (PyTuple_GetItem (point2, 0));
      p2 [1] = PyFloat_AsDouble (PyTuple_GetItem (point2, 1));
      p2 [2] = PyFloat_AsDouble (PyTuple_GetItem (point2, 2));
    }

    if (!(out->con = DOM_Fix_Direction (body->bod->dom, body->bod, p, d, body2?body2->bod:NULL, p2)))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }
    else out->id = out->con->id;

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create prescribed displacement constraint */
static PyObject* lng_SET_DISPLACEMENT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point", "direction", "tms");
  lng_CONSTRAINT *out;
  lng_BODY *body;
  PyObject *point, *direction;
  lng_TIME_SERIES *tms;
  double p [3], d [3];

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOOO", &body, &point, &direction, &tms);

    TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3) &&
	      is_tuple (direction, kwl[2], 3) && is_time_series (tms, kwl[3]));

#if MPI
    if (IS_HERE (body))
    {
#endif

    out->dom = body->bod->dom;

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
    d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
    d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

    if (!(out->con = DOM_Set_Velocity (body->bod->dom, body->bod, p, d, TMS_Derivative (tms->ts))))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }
    else out->id = out->con->id;

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create prescribed velocity constraint */
static PyObject* lng_SET_VELOCITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point", "direction", "value");
  lng_CONSTRAINT *out;
  lng_BODY *body;
  PyObject *point, *direction, *value;
  double p [3], d [3];
  TMS *ts;

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOOO", &body, &point, &direction, &value);

    TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3) &&
	      is_tuple (direction, kwl[2], 3) && is_number_or_time_series (value, kwl[3]));

#if MPI
    if (IS_HERE (body))
    {
#endif

    if (body->bod->dom->dynamic == 0 && body->bod->kind == OBS)
    {
      PyErr_SetString (PyExc_ValueError, "Moving obstacles do not work in the quasi-static case");
      return NULL;
    }

    out->dom = body->bod->dom;

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
    d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
    d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

    if (PyNumber_Check (value)) ts = TMS_Constant (PyFloat_AsDouble (value));
    else ts = TMS_Copy (((lng_TIME_SERIES*)value)->ts);

    if (!(out->con = DOM_Set_Velocity (body->bod->dom, body->bod, p, d, ts)))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }
    else out->id = out->con->id;

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create prescribed acceleration constraint */
static PyObject* lng_SET_ACCELERATION (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point", "direction", "tms");
  lng_CONSTRAINT *out;
  lng_BODY *body;
  PyObject *point, *direction;
  lng_TIME_SERIES *tms;
  double p [3], d [3];

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOOO", &body, &point, &direction, &tms);

    TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3) &&
	      is_tuple (direction, kwl[2], 3) && is_time_series (tms, kwl[3]));

#if MPI
    if (IS_HERE (body))
    {
#endif

    if (body->bod->dom->dynamic == 0 && body->bod->kind == OBS)
    {
      PyErr_SetString (PyExc_ValueError, "Moving obstacles do not work in the quasi-static case");
      return NULL;
    }

    out->dom = body->bod->dom;

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
    d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
    d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

    if (!(out->con = DOM_Set_Velocity (body->bod->dom, body->bod, p, d, TMS_Integral (tms->ts))))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }
    else out->id = out->con->id;

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create rigid link constraint */
static PyObject* lng_PUT_RIGID_LINK (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body1", "body2", "point1", "point2", "strength");
  double p1 [3], p2 [3], strength;
  PyObject *point1, *point2;
  lng_BODY *body1, *body2;
  lng_CONSTRAINT *out;

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    strength = DBL_MAX;

    PARSEKEYS ("OOOO|d", &body1, &body2, &point1, &point2, &strength);

    if ((PyObject*)body1 == Py_None)
    {
      TYPETEST (is_body (body2, kwl[1]) && is_tuple (point1, kwl[2], 3) && is_tuple (point2, kwl[3], 3));
    }
    else if ((PyObject*)body2 == Py_None)
    {
      TYPETEST (is_body (body1, kwl[0]) && is_tuple (point1, kwl[2], 3) && is_tuple (point2, kwl[3], 3));
    }
    else
    {
      TYPETEST (is_body (body1, kwl[0]) && is_body (body2, kwl[1]) &&
		is_tuple (point1, kwl[2], 3) && is_tuple (point2, kwl[3], 3));

      if (body1->bod->dom != body2->bod->dom)
      {
	PyErr_SetString (PyExc_ValueError, "Cannot link bodies from different domains");
	return NULL;
      }

      if (body1->bod->kind == OBS && body2->bod->kind == OBS)
      {
	PyErr_SetString (PyExc_ValueError, "Cannot constrain an obstacle with the rigid link");
	return NULL;
      }
    }

    p1 [0] = PyFloat_AsDouble (PyTuple_GetItem (point1, 0));
    p1 [1] = PyFloat_AsDouble (PyTuple_GetItem (point1, 1));
    p1 [2] = PyFloat_AsDouble (PyTuple_GetItem (point1, 2));

    p2 [0] = PyFloat_AsDouble (PyTuple_GetItem (point2, 0));
    p2 [1] = PyFloat_AsDouble (PyTuple_GetItem (point2, 1));
    p2 [2] = PyFloat_AsDouble (PyTuple_GetItem (point2, 2));

#if MPI
    if (((PyObject*)body1 == Py_None && IS_HERE (body2)) ||
	((PyObject*)body2 == Py_None && IS_HERE (body1)))
    {
#endif

    if ((PyObject*)body1 == Py_None) out->con = DOM_Put_Rigid_Link (body2->bod->dom, NULL, body2->bod, p1, p2, strength);
    else if ((PyObject*)body2 == Py_None) out->con = DOM_Put_Rigid_Link (body1->bod->dom, body1->bod, NULL, p1, p2, strength);
    else out->con = DOM_Put_Rigid_Link (body1->bod->dom, body1->bod, body2->bod, p1, p2, strength);

    out->dom = (PyObject*)body1 == Py_None ? body2->bod->dom : body1->bod->dom;

    if (!out->con)
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }
    else out->id = out->con->id;

#if MPI
    }
    else /* both bodies passed */
    {
      if (body1->bod->dom->time != 0.0)
      {
	PyErr_SetString (PyExc_ValueError, "Rigid links can be inserted only at time zero");
	return NULL;
      }

      if (!DOM_Pending_Constraint (body1->bod->dom, RIGLNK, body1->bod, body2->bod, p1, p2, NULL, NULL, -1, -1, strength))
      {
	PyErr_SetString (PyExc_ValueError, "Point outside of domain");
	return NULL;
      }
    }
#endif
  }

  return (PyObject*)out;
}

/* create slider constraint */
static PyObject* lng_PUT_SPRING (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body1", "point1", "body2", "point2", "function", "limits");
  PyObject *point1, *point2, *function, *limits;
  double p1 [3], p2 [3], lim [2];
  lng_BODY *body1, *body2;
  lng_CONSTRAINT *out;

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOOOOO", &body1, &point1, &body2, &point2, &function, &limits);

    TYPETEST (is_body (body1, kwl[0]) && is_body (body2, kwl[1]) &&
	      is_tuple (point1, kwl[2], 3) && is_tuple (point2, kwl[3], 3) &&
              is_callable (function, kwl[4]) && is_tuple (limits, kwl[5], 2));

    if (body1->bod->dom != body2->bod->dom)
    {
      PyErr_SetString (PyExc_ValueError, "Cannot link bodies from different domains");
      return NULL;
    }

    if (!(body1->bod->kind = RIG || body1->bod->kind == PRB) ||
        !(body2->bod->kind = RIG || body2->bod->kind == PRB))
    {
      PyErr_SetString (PyExc_ValueError, "Slider only works with rigid and pseudo-rigid bodies");
      return NULL;
    }

    p1 [0] = PyFloat_AsDouble (PyTuple_GetItem (point1, 0));
    p1 [1] = PyFloat_AsDouble (PyTuple_GetItem (point1, 1));
    p1 [2] = PyFloat_AsDouble (PyTuple_GetItem (point1, 2));

    p2 [0] = PyFloat_AsDouble (PyTuple_GetItem (point2, 0));
    p2 [1] = PyFloat_AsDouble (PyTuple_GetItem (point2, 1));
    p2 [2] = PyFloat_AsDouble (PyTuple_GetItem (point2, 2));

    lim [0] = PyFloat_AsDouble (PyTuple_GetItem (limits, 0));
    lim [1] = PyFloat_AsDouble (PyTuple_GetItem (limits, 1));

#if MPI
    if (((PyObject*)body1 == Py_None && IS_HERE (body2)) ||
	((PyObject*)body2 == Py_None && IS_HERE (body1)))
    {
#endif

    out->con = DOM_Put_Spring (body1->bod->dom, body1->bod, p1, body2->bod, p2, function, lim);

    out->dom = body1->bod->dom;

    if (!out->con)
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }
    else out->id = out->con->id;

#if MPI
    }
    else /* both bodies passed */
    {
      if (body1->bod->dom->time != 0.0)
      {
	PyErr_SetString (PyExc_ValueError, "Springs can be inserted only at time zero");
	return NULL;
      }

      if (!DOM_Pending_Constraint (body1->bod->dom, SPRING, body1->bod, body2->bod, p1, p2, lim, (TMS*)function, -1, -1, DBL_MAX))
      {
	PyErr_SetString (PyExc_ValueError, "Point outside of domain");
	return NULL;
      }
    }
#endif
  }

  return (PyObject*)out;
}

/* set gravity */
static PyObject* lng_GRAVITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "vector");
  PyObject *vector, *val [3];
  lng_SOLFEC *solfec;
  DOM *dom;
  TMS *ts;
  int i;

  PARSEKEYS ("OO", &solfec, &vector);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_tuple (vector, kwl[1], 3));

  
  val [0] = PyTuple_GetItem (vector, 0);
  val [1] = PyTuple_GetItem (vector, 1);
  val [2] = PyTuple_GetItem (vector, 2);

  TYPETEST (is_number_or_time_series (val[0], kwl[1]) &&
            is_number_or_time_series (val[1], kwl[1]) &&
            is_number_or_time_series (val[2], kwl[1]));

  dom = solfec->sol->dom;

  for (i = 0; i < 3; i ++)
  {
    if (PyNumber_Check (val [i])) ts = TMS_Constant (PyFloat_AsDouble(val [i]));
    else ts = TMS_Copy (((lng_TIME_SERIES*)val [i])->ts);

    if (dom->gravity [i]) TMS_Destroy (dom->gravity [i]);

    dom->gravity [i] = ts;
  }

  Py_RETURN_NONE;
}

/* callback for forces defined by a callable object */
static void lng_FORCE_callback (PyObject *data, PyObject *call, int nq, double *q, int nu, double *u, double t, double h, double *f)
{
  PyObject *qtup, *utup;
  PyObject *result;
  PyObject *args;
  int n;

  if (!(qtup = PyTuple_New (nq))) goto err;
  for (n = 0; n < nq; n ++) PyTuple_SetItem (qtup, n, PyFloat_FromDouble (q [n]));

  if (!(utup = PyTuple_New (nq))) goto err;
  for (n = 0; n < nu; n ++) PyTuple_SetItem (utup, n, PyFloat_FromDouble (u [n]));

  if (data)
  {
    if (PyTuple_Check (data))
    {
      if (!(args = PyTuple_New (PyTuple_Size (data) + 4))) goto err;

      for (n = 0; n < PyTuple_Size (data); n ++)
	PyTuple_SetItem (args, n, PyTuple_GetItem (data, n));

      PyTuple_SetItem (args, n, qtup);
      PyTuple_SetItem (args, n+1, utup);
      PyTuple_SetItem (args, n+2, PyFloat_FromDouble (t));
      PyTuple_SetItem (args, n+3, PyFloat_FromDouble (h));
    }
    else args = Py_BuildValue ("(O, O, O, d, d)", data, qtup, utup, t, h);
  }
  else args = Py_BuildValue ("(O, O, d, d)", qtup, utup, t, h);

  result = PyObject_CallObject (call, args); /* call user callback */

  Py_DECREF (args);
  Py_DECREF (qtup);
  Py_DECREF (utup);

  if (result)
  {
    if (!PyTuple_Check (result))
    {
      PyErr_SetString (PyExc_ValueError, "FORCE callback routine must reurn a tuple");
      goto err;
    }

    if (PyTuple_Size (result) != (nu == 6 ? 9 : nu))
    {
      PyErr_SetString (PyExc_ValueError, "Tuple of invalid size was returned from FORCE callback");
      goto err;
    }

    for (n = 0; n < PyTuple_Size (result); n ++)
      f [n] = PyFloat_AsDouble (PyTuple_GetItem (result, n));

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
}

/* apply point force */
static PyObject* lng_FORCE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "kind", "point", "direction", "value", "data");
  PyObject *kind, *point, *direction, *value, *data;
  double p [3], d [3];
  lng_BODY *body;
  short k;
  TMS *ts;
  FORCE_FUNC func;
  PyObject *call;

  data = NULL;
  func = NULL;
  call = NULL;

  PARSEKEYS ("OOOOO|O", &body, &kind, &point, &direction, &value, &data);

  TYPETEST (is_body (body, kwl[0]) && is_string (kind, kwl[1]) && is_tuple (point, kwl[2], 3) &&
            is_tuple (direction, kwl[3], 3) && is_number_or_time_series_or_callable (value, kwl[4]));

#if MPI && LOCAL_BODIES
  if (IS_HERE (body))
  {
#endif
 
  if (body->bod->kind == OBS)
  {
    PyErr_SetString (PyExc_ValueError, "Cannot load an obstacle");
    return NULL;
  }

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
  d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
  d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

  if (PyNumber_Check (value)) ts = TMS_Constant (PyFloat_AsDouble(value));
  else if (PyCallable_Check (value))
  {
    ts = (TMS*) data;
    call = value;
    func = (FORCE_FUNC) lng_FORCE_callback;
    callback_pair_push (data, value);
  }
  else ts = TMS_Copy (((lng_TIME_SERIES*)value)->ts);

  IFIS (kind, "SPATIAL")
  {
    k = SPATIAL;
  }
  ELIF (kind, "CONVECTED")
  {
    k = CONVECTED;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid force kind");
    return NULL;
  }

  if (body->bod->kind == FEM)
  {
    if (SHAPE_Sgp (body->bod->sgp, body->bod->nsgp, p) < 0)
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of a finite-element body shape");
      return NULL;
    }
  }

  BODY_Apply_Force (body->bod, k, p, d, ts, call, func, 0);

#if MPI && LOCAL_BODIES
  }
#endif

  Py_RETURN_NONE;
}

/* apply point torque */
static PyObject* lng_TORQUE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "kind", "direction", "value");
  PyObject *kind, *direction, *value;
  lng_BODY *body;
  double d [3];
  short k;
  TMS *ts;

  PARSEKEYS ("OOOO", &body, &kind, &direction, &value);

  TYPETEST (is_body (body, kwl[0]) && is_string (kind, kwl[1]) &&
            is_tuple (direction, kwl[2], 3) &&
	    is_number_or_time_series (value, kwl[3]));

#if MPI && LOCAL_BODIES
  if (IS_HERE (body))
  {
#endif

  if (body->bod->kind != RIG)
  {
    PyErr_SetString (PyExc_ValueError, "Torque can be appliedy to RIGID bodies only");
    return NULL;
  }

  d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
  d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
  d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

  if (PyNumber_Check (value)) ts = TMS_Constant (PyFloat_AsDouble(value));
  else ts = TMS_Copy (((lng_TIME_SERIES*)value)->ts);

  IFIS (kind, "SPATIAL")
  {
    k = SPATIAL;
  }
  ELIF (kind, "CONVECTED")
  {
    k = CONVECTED;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid force kind");
    return NULL;
  }

  BODY_Apply_Force (body->bod, k | TORQUE, NULL, d, ts, NULL, NULL, 0);

#if MPI && LOCAL_BODIES
  }
#endif

  Py_RETURN_NONE;
}

/* apply pressure */
static PyObject* lng_PRESSURE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "surfid", "value");
  PyObject *value;
  lng_BODY *body;
  int surfid;
  BODY *bod;
  TMS *ts;

  PARSEKEYS ("OiO", &body, &surfid, &value);

  TYPETEST (is_body (body, kwl[0]) && is_number_or_time_series (value, kwl[2]));

#if MPI && LOCAL_BODIES
  if (IS_HERE (body))
  {
#endif

  bod = body->bod;

  if (bod->shape->kind != SHAPE_MESH || bod->shape->next != NULL)
  {
    PyErr_SetString (PyExc_ValueError, "Pressure can be applied to a body with a single MESH based shape");
    return NULL;
  }

  if (PyNumber_Check (value)) ts = TMS_Constant (PyFloat_AsDouble(value));
  else ts = TMS_Copy (((lng_TIME_SERIES*)value)->ts);

  BODY_Apply_Force (body->bod, PRESSURE, NULL, NULL, ts, NULL, NULL, surfid);

#if MPI && LOCAL_BODIES
  }
#endif

  Py_RETURN_NONE;
}

/* simplified crack */
static PyObject* lng_SIMPLIFIED_CRACK (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point", "normal", "surfid", "criterion", "topoadj", "ft", "Gf");
  PyObject *point, *normal, *surfid, *criterion, *topoadj;
  double p [3], n [3], ft, Gf;
  lng_BODY *body;
  CRACK *cra;

  topoadj = NULL;
  ft = 0.0;
  Gf = 0.0;

  PARSEKEYS ("OOOOO|Odd", &body, &point, &normal,
             &surfid, &criterion, &topoadj, &ft, &Gf);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3) && is_tuple (normal, kwl[2], 3) &&
            is_tuple (surfid, kwl[3], 2) && is_string (criterion, kwl[4]) && is_string (topoadj, kwl[5]));

#if MPI && LOCAL_BODIES
  if (IS_HERE (body))
  {
#endif

  if (!(body->bod->kind == PRB || body->bod->kind == FEM))
  {
    PyErr_SetString (PyExc_ValueError, "Invalid body kind: only PSEUDO_RIGID and FINITE_ELEMENT bodies can crack");
    return NULL;
  }

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  n [0] = PyFloat_AsDouble (PyTuple_GetItem (normal, 0));
  n [1] = PyFloat_AsDouble (PyTuple_GetItem (normal, 1));
  n [2] = PyFloat_AsDouble (PyTuple_GetItem (normal, 2));

  NORMALIZE (n);

  cra = CRACK_Create ();
  COPY (p, cra->point);
  COPY (n, cra->normal);
  cra->surfid[0] = PyInt_AsLong (PyTuple_GetItem (surfid, 0));
  cra->surfid[1] = PyInt_AsLong (PyTuple_GetItem (surfid, 1));

  if (topoadj)
  {
    IFIS (topoadj, "ON") cra->topoadj = 1;
    ELIF (topoadj, "OFF") cra->topoadj = 0;
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid topoadj value: neither 'ON' nor 'OFF'");
      return NULL;
    }
  }
  else cra->topoadj = 0;

  IFIS (criterion, "TENSILE")
  {
    cra->crit = TENSILE;
    cra->ft = ft;
    cra->Gf = Gf;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid criterion");
    return NULL;
  }

  if (SHAPE_Cut_Possible (body->bod->shape, 1, cra->point, cra->normal, cra->topoadj))
  {
    cra->next = body->bod->cra;
    body->bod->cra = cra;
  }
  else
  {
    if (cra->topoadj && SHAPE_Cut_Possible (body->bod->shape, 1, cra->point, cra->normal, 0))
      PyErr_SetString (PyExc_ValueError, "The input point does not belong the body shape.");
    else PyErr_SetString (PyExc_ValueError, "The crack plane does not intersect the body shape.");
    CRACK_Destroy (cra);
    return NULL;
  }

#if MPI && LOCAL_BODIES
  }
#endif

  Py_RETURN_NONE;
}

/* set imbalance tolerances */
static PyObject* lng_IMBALANCE_TOLERANCE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "tolerance", "weightfactor", "updatefreq");
  lng_SOLFEC *solfec;
  double tolerance,
	 weightfactor;
  int updatefreq;

  weightfactor = 1.0;
  updatefreq = 10;

  PARSEKEYS ("Od|di", &solfec, &tolerance, &weightfactor, &updatefreq);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_positive (tolerance, kwl[1]) &&
            is_ge_le (weightfactor, kwl [2], 0.0, 1.0) && is_ge (updatefreq, kwl [3], 1));
#if MPI
  solfec->sol->dom->weight_factor = weightfactor;
  solfec->sol->dom->updatefreq = updatefreq;
#endif

  Py_RETURN_TRUE;
}

/* return current rank */
static PyObject* lng_RANK (PyObject *self, PyObject *args, PyObject *kwds)
{
  int rank;

#if MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
#else
  rank = 0;
#endif

  return PyInt_FromLong (rank);
}

/* set up parallel barrier */
static PyObject* lng_BARRIER (PyObject *self, PyObject *args, PyObject *kwds)
{
#if MPI
  MPI_Barrier (MPI_COMM_WORLD);
#endif

  Py_RETURN_TRUE;
}

/* return number of CPUs */
static PyObject* lng_NCPU (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec");
  lng_SOLFEC *solfec;
  int ncpu;

  PARSEKEYS ("O", &solfec);

  TYPETEST (is_solfec (solfec, kwl[0]));

#if MPI
  MPI_Comm_size (MPI_COMM_WORLD, &ncpu);
#else
  PBF *bf;

  for (bf = solfec->sol->bf, ncpu = 0; bf; bf = bf->next) ncpu ++;
#endif

  return PyInt_FromLong (ncpu);
}

/* test whether an object is on this processor */
static PyObject* lng_HERE (PyObject *self, PyObject *args, PyObject *kwds)
{
#if MPI
  KEYWORDS ("solfec", "object");
  lng_CONSTRAINT *constraint;
  lng_SOLFEC *solfec;
  PyObject *object;
  lng_BODY *body;

  PARSEKEYS ("OO", &solfec, &object);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_body_or_constraint (object, kwl[1]));

  if (PyObject_IsInstance (object, (PyObject*)&lng_BODY_TYPE))
  {
    body = (lng_BODY*)object;

    if (IS_HERE (body)) Py_RETURN_TRUE;
    else Py_RETURN_FALSE;
  }
  else
  {
    constraint = (lng_CONSTRAINT*)object;

    if (ID_TO_CONSTRAINT (solfec->sol->dom, constraint)) Py_RETURN_TRUE;
    else Py_RETURN_FALSE;
  }
#else
  Py_RETURN_TRUE;
#endif
}

/* test whether the viewer is enabled */
static PyObject* lng_VIEWER (PyObject *self, PyObject *args, PyObject *kwds)
{
#if OPENGL
  if (RND_Is_On ()) Py_RETURN_TRUE; /* return true and maintain reference count of Py_True */
  else
#endif
  Py_RETURN_FALSE; /* return false and maintain reference count of Py_False */
}

/* overwrite body characteristics */
static PyObject* lng_BODY_CHARS (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "mass", "volume", "center", "tensor");
  lng_BODY *body;
  double mass, volume;
  PyObject *center, *tensor;
  double c [3], t [9];
  int i;

  PARSEKEYS ("OddOO", &body, &mass, &volume, &center, &tensor);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (center, kwl[3], 3) && is_tuple (tensor, kwl[4], 9));

#if MPI && LOCAL_BODIES
  if (IS_HERE (body))
  {
#endif

  c [0] = PyFloat_AsDouble (PyTuple_GetItem (center, 0));
  c [1] = PyFloat_AsDouble (PyTuple_GetItem (center, 1));
  c [2] = PyFloat_AsDouble (PyTuple_GetItem (center, 2));

  for (i = 0; i < 9; i ++) t [i] = PyFloat_AsDouble (PyTuple_GetItem (tensor, i));

  BODY_Overwrite_Chars (body->bod, mass, volume, c, t);

#if MPI && LOCAL_BODIES
  }
#endif

  Py_RETURN_NONE;
}

/* apply initial velocity */
static PyObject* lng_INITIAL_VELOCITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "linear", "angular");
  lng_BODY *body;
  PyObject *linear, *angular;
  double l [3], a [3];

  PARSEKEYS ("OOO", &body, &linear, &angular);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (linear, kwl[1], 3) && is_tuple (angular, kwl[2], 3));

#if MPI && LOCAL_BODIES
  if (IS_HERE (body))
  {
#endif

  l [0] = PyFloat_AsDouble (PyTuple_GetItem (linear, 0));
  l [1] = PyFloat_AsDouble (PyTuple_GetItem (linear, 1));
  l [2] = PyFloat_AsDouble (PyTuple_GetItem (linear, 2));

  a [0] = PyFloat_AsDouble (PyTuple_GetItem (angular, 0));
  a [1] = PyFloat_AsDouble (PyTuple_GetItem (angular, 1));
  a [2] = PyFloat_AsDouble (PyTuple_GetItem (angular, 2));

  BODY_Initial_Velocity (body->bod, l, a);

#if MPI && LOCAL_BODIES
  }
#endif

  Py_RETURN_NONE;
}

/* apply bulk material */
static PyObject* lng_MATERIAL (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "body", "volid", "material");
  PyObject *material;
  lng_SOLFEC *solfec;
  lng_BODY *body;
  int volid;

  PARSEKEYS ("OOdO", &solfec, &body, &volid, &material);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body, kwl[1]) &&
            is_bulk_material (solfec->sol, material, kwl[2]));

  BULK_MATERIAL *mat0 = body->bod->mat,
                *mat1 = get_bulk_material (solfec->sol, material);

  if (mat0->nfield != mat1->nfield || mat0->nstate != mat1->nstate)
  {
    PyErr_SetString (PyExc_ValueError, "Number of field variables and state variables must be the same as for the body global material!");
    return NULL;
  }

  for (int i = 0; i < mat0->nfield; i ++)
  {
    if (mat0->fld [i] != mat1->fld [i])
    {
      PyErr_SetString (PyExc_ValueError, "Fields must match those of the body global material!");
      return NULL;
    }
  }

  BODY_Material (body->bod, volid, mat1);

  Py_RETURN_NONE;
}

/* delete body or constraint */
static PyObject* lng_DELETE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "object");
  lng_CONSTRAINT *constraint;
  lng_SOLFEC *solfec;
  PyObject *object;
  lng_BODY *body;

  PARSEKEYS ("OO", &solfec, &object);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_body_or_constraint (object, kwl[1]));

  if (PyObject_IsInstance (object, (PyObject*)&lng_BODY_TYPE))
  {
    body = (lng_BODY*)object;

#if MPI
    if (IS_HERE (body))
    {
#endif

    DOM_Remove_Body (solfec->sol->dom, body->bod);
    BODY_Destroy (body->bod); /* used only when body was removed from the domain */
    body->bod = NULL;

#if MPI
    }
#endif
  }
  else
  {
    constraint = (lng_CONSTRAINT*)object;

    if (ID_TO_CONSTRAINT (solfec->sol->dom, constraint))
    {
      DOM_Remove_Constraint (solfec->sol->dom, constraint->con);
    }
  }

  Py_RETURN_NONE;
}

/* scale shape */
static PyObject* lng_SCALE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "coefs");
  PyObject *shape, *coefs;
  double c [3], s [3];
  SHAPE *shp;

  PARSEKEYS ("OO", &shape, &coefs);

  TYPETEST (is_shape_or_vector (shape, kwl[0]) && is_tuple (coefs, kwl[1], 3));

  c [0] = PyFloat_AsDouble (PyTuple_GetItem (coefs, 0));
  c [1] = PyFloat_AsDouble (PyTuple_GetItem (coefs, 1));
  c [2] = PyFloat_AsDouble (PyTuple_GetItem (coefs, 2));

  if (PyTuple_Check (shape)) /* vector */
  {
    s [0] = PyFloat_AsDouble (PyTuple_GetItem (shape, 0));
    s [1] = PyFloat_AsDouble (PyTuple_GetItem (shape, 1));
    s [2] = PyFloat_AsDouble (PyTuple_GetItem (shape, 2));

    /* doing PyTuple_SetItem (shape, ..., PyFLoat_FromDouble (...)) is not
     * possible as tuple is immutable (cannot be modified after creation) */

    return Py_BuildValue ("(d, d, d)", s[0] * c[0], s[1] * c[1], s[2] * c[2]);
  }
  else /* shape */
  {
    shp = create_shape (shape, 0);

    SHAPE_Scale (shp, c);

    SHAPE_Destroy_Wrapper (shp);
  }

  Py_INCREF (shape);
  return shape;
}

/* translate shape */
static PyObject* lng_TRANSLATE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "vector");
  PyObject *shape, *vector;
  double v [3], s [3];
  SHAPE *shp;

  PARSEKEYS ("OO", &shape, &vector);

  TYPETEST (is_shape_or_vector (shape, kwl[0]) && is_tuple (vector, kwl[1], 3));

  v [0] = PyFloat_AsDouble (PyTuple_GetItem (vector, 0));
  v [1] = PyFloat_AsDouble (PyTuple_GetItem (vector, 1));
  v [2] = PyFloat_AsDouble (PyTuple_GetItem (vector, 2));

  if (PyTuple_Check (shape)) /* vector */
  {
    s [0] = PyFloat_AsDouble (PyTuple_GetItem (shape, 0));
    s [1] = PyFloat_AsDouble (PyTuple_GetItem (shape, 1));
    s [2] = PyFloat_AsDouble (PyTuple_GetItem (shape, 2));

    /* doing PyTuple_SetItem (shape, ..., PyFLoat_FromDouble (...)) is not
     * possible as tuple is immutable (cannot be modified after creation) */

    return Py_BuildValue ("(d, d, d)", s[0] + v[0], s[1] + v[1], s[2] + v[2]);
  }
  else /* shape */
  {
    shp = create_shape (shape, 0);

    SHAPE_Translate (shp, v);

    SHAPE_Destroy_Wrapper (shp);
  }

  Py_INCREF (shape);
  return shape;
}

/* rotate shape */
static PyObject* lng_ROTATE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "point", "vector", "angle");
  PyObject *shape, *point, *vector;
  double r [9], t [3], s [3], p [3], v [3], angle;
  SHAPE *shp;

  PARSEKEYS ("OOOd", &shape, &point, &vector, &angle);

  TYPETEST (is_shape_or_vector (shape, kwl[0]) && is_tuple (point, kwl[1], 3) && is_tuple (vector, kwl[2], 3));

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  v [0] = PyFloat_AsDouble (PyTuple_GetItem (vector, 0));
  v [1] = PyFloat_AsDouble (PyTuple_GetItem (vector, 1));
  v [2] = PyFloat_AsDouble (PyTuple_GetItem (vector, 2));

  if (LEN (v) == 0)
  {
    PyErr_SetString (PyExc_ValueError, "Rotation direction is zero");
    return NULL;
  }

  if (PyTuple_Check (shape)) /* vector */
  {
    s [0] = PyFloat_AsDouble (PyTuple_GetItem (shape, 0)) - p[0];
    s [1] = PyFloat_AsDouble (PyTuple_GetItem (shape, 1)) - p[1];
    s [2] = PyFloat_AsDouble (PyTuple_GetItem (shape, 2)) - p[2];

    COPY (v, t); 
    NORMALIZE (t);
    SCALE (t, (ALG_PI * angle / 180.0));
    EXPMAP (t, r);
    NVADDMUL (p, r, s, t);

    /* doing PyTuple_SetItem (shape, ..., PyFLoat_FromDouble (...)) is not
     * possible as tuple is immutable (cannot be modified after creation) */

    return Py_BuildValue ("(d, d, d)", t[0], t[1], t[2]);
  }
  else  /* shape */
  {
    shp = create_shape (shape, 0);

    SHAPE_Rotate (shp, p, v, angle);

    SHAPE_Destroy_Wrapper (shp);
  }

  Py_INCREF (shape);
  return shape;
}

/* convert shape into PyList of basic shapes */
static PyObject* shape_to_list (SHAPE *shp)
{
  if (shp)
  {
    lng_MESH *mesh;
    lng_CONVEX *convex;
    lng_SPHERE *sphere;
    lng_ELLIP *ellip;
    PyObject *list;
    SHAPE *shq;
    int n;

    for (n = 0, shq = shp; shq; shq = shq->next) n ++;

    if (n > 1) /* create list of object */
    {
      if (!(list = PyList_New (n))) return NULL;

      for (n = 0, shq = shp; shq; shq = shq->next, n ++)
      {
	switch (shq->kind)
	{
	case SHAPE_MESH:
	  mesh = (lng_MESH*)lng_MESH_TYPE.tp_alloc (&lng_MESH_TYPE, 0);
	  if (mesh)
	  {
	    mesh->msh = shq->data;
	    PyList_SetItem (list, n, (PyObject*)mesh);
	  }
	  else return NULL;
	  break;
	case SHAPE_CONVEX:
	  convex = (lng_CONVEX*)lng_CONVEX_TYPE.tp_alloc (&lng_CONVEX_TYPE, 0);
	  if (convex)
	  {
	    convex->cvx = shq->data;
	    PyList_SetItem (list, n, (PyObject*)convex);
	  }
	  else return NULL;
	  break;
	case SHAPE_SPHERE:
	  sphere = (lng_SPHERE*)lng_SPHERE_TYPE.tp_alloc (&lng_SPHERE_TYPE, 0);
	  if (sphere)
	  {
	    sphere->sph = shq->data;
	    PyList_SetItem (list, n, (PyObject*)sphere);
	  }
	  else return NULL;
	  break;
        case SHAPE_ELLIP:
	  ellip = (lng_ELLIP*)lng_ELLIP_TYPE.tp_alloc (&lng_ELLIP_TYPE, 0);
	  if (ellip)
	  {
	    ellip->eli = shq->data;
	    PyList_SetItem (list, n, (PyObject*)ellip);
	  }
	  else return NULL;
	  break;

	}
      }

      return list;
    }
    else /* create single object */
    {
      switch (shp->kind)
      {
      case SHAPE_MESH:
	mesh = (lng_MESH*)lng_MESH_TYPE.tp_alloc (&lng_MESH_TYPE, 0);
	if (mesh)
	{
	  mesh->msh = shp->data;
	  return (PyObject*)mesh;
	}
	else return NULL;
	break;
      case SHAPE_CONVEX:
	convex = (lng_CONVEX*)lng_CONVEX_TYPE.tp_alloc (&lng_CONVEX_TYPE, 0);
	if (convex)
	{
	  convex->cvx = shp->data;
	  return (PyObject*)convex;
	}
	else return NULL;
	break;
      case SHAPE_SPHERE:
	sphere = (lng_SPHERE*)lng_SPHERE_TYPE.tp_alloc (&lng_SPHERE_TYPE, 0);
	if (sphere)
	{
	  sphere->sph = shp->data;
	  return (PyObject*)sphere;
	}
	else return NULL;
	break;
      case SHAPE_ELLIP:
	ellip = (lng_ELLIP*)lng_ELLIP_TYPE.tp_alloc (&lng_ELLIP_TYPE, 0);
	if (ellip)
	{
	  ellip->eli = shp->data;
	  return (PyObject*)ellip;
	}
	else return NULL;
	break;
      }
    }
  }

  Py_RETURN_NONE;
}

/* split shape by plane */
static PyObject* lng_SPLIT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "point", "normal", "surfid", "topoadj", "remesh");
  PyObject *shape, *point, *normal, *surfid, *topoadj, *remesh, *out, *b, *f;
  SHAPE *shp, *shq, *back, *front, *next;
  double p [3], n [3];
  int error, surf[2];
  short tadj, rmsh;

  surf [0] = surf [1] = 0;
  surfid = NULL;
  topoadj = NULL;
  remesh = NULL;

  PARSEKEYS ("OOO|OOO", &shape, &point, &normal, &surfid, &topoadj, &remesh);

  TYPETEST (is_shape (shape, kwl[0]) && is_tuple (point, kwl[1], 3) && is_tuple (normal, kwl[2], 3)
         && is_tuple (surfid, kwl[3], 2) && is_string (topoadj, kwl [4]) && is_string (remesh, kwl[5]));

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  n [0] = PyFloat_AsDouble (PyTuple_GetItem (normal, 0));
  n [1] = PyFloat_AsDouble (PyTuple_GetItem (normal, 1));
  n [2] = PyFloat_AsDouble (PyTuple_GetItem (normal, 2));

  if (surfid)
  {
    surf [0] = PyInt_AsLong (PyTuple_GetItem (surfid, 0));
    surf [1] = PyInt_AsLong (PyTuple_GetItem (surfid, 1));
  }

  if (topoadj)
  {
    IFIS (topoadj, "ON") tadj = 1;
    ELIF (topoadj, "OFF") tadj = 0;
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid topoadj value: neither 'ON' nor 'OFF'");
      return NULL;
    }
  }
  else tadj = 0;

  if (remesh)
  {
    IFIS (remesh, "ON") rmsh = 1;
    ELIF (remesh, "OFF") rmsh = 0;
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid remesh value: neither 'ON' nor 'OFF'");
      return NULL;
    }
  }
  else rmsh = 1;

  shp = create_shape (shape, 1); /* empty */

  back = front = NULL;

  TRY ()
  {
    for (shq = shp; shq; shq = shq->next)
    {
      switch (shq->kind)
      {
      case SHAPE_CONVEX:
      {
	CONVEX *one = NULL, *two = NULL;

	CONVEX_Split (shq->data, p, n, tadj, surf, &one, &two);
	if (one) back = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, one), back);
	if (two) front = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, two), front);
      }
      break;
      case SHAPE_SPHERE:
      {
	CONVEX *one = NULL, *two = NULL;

	SPHERE_Split (shq->data, p, n, tadj, surf, &one, &two);
	if (one && two)
	{
	  back = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, one), back);
	  front = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, two), front);
	}
	else if (one) back = SHAPE_Glue (SHAPE_Create (SHAPE_SPHERE, one), back);
	else if (two) front = SHAPE_Glue (SHAPE_Create (SHAPE_SPHERE, two), front);
      }
      break;
      case SHAPE_ELLIP:
      {
	CONVEX *one = NULL, *two = NULL;

	ELLIP_Split (shq->data, p, n, tadj, surf, &one, &two);
	if (one && two)
	{
	  back = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, one), back);
	  front = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, two), front);
	}
	else if (one) back = SHAPE_Glue (SHAPE_Create (SHAPE_ELLIP, one), back);
	else if (two) front = SHAPE_Glue (SHAPE_Create (SHAPE_ELLIP, two), front);
      }
      break;
      case SHAPE_MESH:
      {
	MESH *one = NULL, *two = NULL;

	if (MESH_Split (shq->data, p, n, tadj, surf, rmsh, &one, &two) == 1)
	{
	  PyErr_SetString (PyExc_ValueError, "Inter-element boundary mesh splitting failed while remeshing is disabled");
	  return NULL;
	}

	if (one) back = SHAPE_Glue (SHAPE_Create (SHAPE_MESH, one), back);
	if (two) front = SHAPE_Glue (SHAPE_Create (SHAPE_MESH, two), front);
      }
      break;
      }
    }
  }
  CATCHANY (error)
  {
    PyErr_SetString (PyExc_RuntimeError, errstring (error));
#if MPI
      PyErr_Print ();
      MPI_Abort (MPI_COMM_WORLD, 4000+error);
#endif
    return NULL;
  }
  ENDTRY ()

  SHAPE_Destroy (shp);

  shp = SHAPE_Glue (back, front); /* glue front and back into one list */
  back = front = NULL; /* and empty these lists for the moment */

  for (shq = shp; shq; shq = next) /* for each item */
  {
    double limits [2];

    next = shq->next;

    if (SHAPE_Separable (shq)) /* check if separable */
    {
      SHAPE **s;
      int m;

      s = SHAPE_Separate (shq, &m); /* separate */
      SHAPE_Destroy (shq);

      for (m --; m >= 0; m --) /* and classify each part */
      {
        SHAPE_Limits_Along_Line (s [m], p, n, limits);
	if (fabs (limits [0]) < fabs (limits [1])) { s [m]->next = back; back = s [m]; }
	else { s [m]->next = front; front = s [m]; }
      }
    }
    else /* or classify the current sub-shape */
    {
      SHAPE_Limits_Along_Line (shq, p, n, limits);
      if (fabs (limits [0]) < fabs (limits [1])) { shq->next = back; back = shq; }
      else { shq->next = front; front = shq; }
    }
  }

  b =  shape_to_list (back);
  f =  shape_to_list (front);

  if (!(b||f)) return NULL;

  out = Py_BuildValue ("(O, O)", b, f);

  if (back) SHAPE_Destroy_Wrapper (back);
  if (front) SHAPE_Destroy_Wrapper (front);

  return out;
}

/* split mesh by node set */
static PyObject* lng_MESH_SPLIT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("mesh", "nodeset", "surfid");
  PyObject *nodeset, *surfid, *list;
  int i, nin, nout;
  lng_MESH *mesh;
  SET *nodes;
  MESH **out;

  surfid = NULL;

  PARSEKEYS ("OO|O", &mesh, &nodeset, &surfid);

  TYPETEST (is_mesh ((PyObject*)mesh, kwl[0]) && is_list (nodeset, kwl[1], 0, 0));

  if (mesh->msh == NULL)
  {
    PyErr_SetString (PyExc_RuntimeError, "The mesh object is empty");
    return NULL;
  }

  nodes = NULL;
  nin = PyList_Size (nodeset);
  for (i = 0; i < nin; i ++)
  {
    SET_Insert (NULL, &nodes, (void*)PyInt_AsLong (PyList_GetItem (nodeset,i)), NULL);
  }

  if (surfid) i = PyInt_AsLong (surfid);
  else i = 0;

  out = MESH_Split_By_Nodes (mesh->msh, nodes, i, &nout);

  SET_Free (NULL, &nodes);

  if (out)
  {
    list = PyList_New (nout);
    for (i = 0; i < nout; i ++)
    {
      lng_MESH *item = (lng_MESH*)lng_MESH_TYPE.tp_alloc (&lng_MESH_TYPE, 0);
      if (item)
      {
	item->msh = out[i];
	PyList_SetItem (list, i, (PyObject*)item);
      }
      else return NULL;
    }

    free (out);

    return list;
  }
  else Py_RETURN_NONE;
}

/* copy shape */
static PyObject* lng_COPY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape");
  PyObject *shape, *out;
  SHAPE *shp, *shq;

  PARSEKEYS ("O", &shape);

  TYPETEST (is_shape (shape, kwl[0]));

  shp = create_shape (shape, 0);

  for (shq = shp; shq; shq = shq->next)
  {
    switch (shq->kind)
    {
    case SHAPE_MESH:
      shq->data = MESH_Copy (shq->data);
      break;
    case SHAPE_CONVEX:
      shq->data = CONVEX_Copy (shq->data);
      break;
    case SHAPE_SPHERE:
      shq->data = SPHERE_Copy (shq->data);
      break;
    case SHAPE_ELLIP:
      shq->data = ELLIP_Copy (shq->data);
      break;
    }
  }

  out = shape_to_list (shp);

  SHAPE_Destroy_Wrapper (shp);

  return out;
}

/* bend shape */
static PyObject* lng_BEND (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "point", "direction", "angle");
  PyObject *point, *direction;
  double p [3], d [3], angle;
  lng_MESH *shape;
  MESH *msh;

  PARSEKEYS ("OOOd", &shape, &point, &direction, &angle);

  TYPETEST (is_mesh ((PyObject*)shape, kwl[0]) && is_tuple (point, kwl[1], 3)
         && is_tuple (direction, kwl[2], 3) && is_positive (angle, kwl[3]));

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
  d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
  d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

  msh = shape->msh;

  angle = ALG_PI * angle / 180.0;

  double (*nod) [3], (*end) [3], (*sel) [3], proj [3],
	 dif [3], q [3], v [3], dot, dmin, rmin, r, alpha;
  ELEMENT *ele;

  /* test for element-axis intersection */
  for (ele = msh->surfeles, nod = msh->ref_nodes; ele; ele = ele->next)
  {
    SET (q, 0.0);
    for (int i = 0; i < ele->type; i ++)
    {
      sel = &nod [ele->nodes [i]];
      ACC (sel[0], q);
    }
    DIV (q, (double) ele->type, q);
    PROJECT_POINT_ON_LINE (q, p, d, proj); /* project element center on the bending axis */
    if (ELEMENT_Contains_Point (msh, ele, proj, 1))
    {
      PyErr_SetString (PyExc_ValueError, "Bending axis is stabbing the mesh");
      return NULL;
    }
  }

  dmin = DBL_MAX;
  for (nod = msh->ref_nodes, end = nod + msh->nodes_count; nod != end; nod ++)
  {
    PROJECT_POINT_ON_LINE (nod[0], p, d, proj);
    SUB (nod[0], proj, dif);
    dot = DOT (dif, dif);
    if (dot < dmin) { dmin = dot; sel = nod; } /* find closest mesh node */
  }
  PROJECT_POINT_ON_LINE (sel[0], p, d, proj);
  SUB (sel[0], proj, dif);
  PRODUCT (d, dif, v);
  rmin = sqrt (dmin); /* shortest bending radius */
  NORMALIZE (v); /* v is the direction of angle growth */

  for (nod = msh->ref_nodes, sel = msh->cur_nodes, end = nod + msh->nodes_count; nod != end; nod ++, sel ++)
  {
    PROJECT_POINT_ON_LINE (nod[0], p, d, proj);
    SUB (nod[0], proj, dif);
    dot = DOT (v, dif); /* arc length */
    if (dot > 0)
    {
      SUBMUL (nod[0], dot, v, q); /* project node on (axis, v) plane */
      PROJECT_POINT_ON_LINE (q, p, d, proj);
      SUB (q, proj, dif);
      r = LEN (dif); /* bending radius  = |proj (line, proj (plane, node)) - proj (plane, node)| */
      alpha = dot / rmin;

      if (alpha <= angle)
      {
        nod[0][0] = proj[0] + cos (alpha) * dif [0] + r * sin (alpha) * v [0];
        nod[0][1] = proj[1] + cos (alpha) * dif [1] + r * sin (alpha) * v [1];
        nod[0][2] = proj[2] + cos (alpha) * dif [2] + r * sin (alpha) * v [2];
      }
      else
      {
	double u [3], len;

	/* bend up to the angle */
        nod[0][0] = proj[0] + cos (angle) * dif [0] + r * sin (angle) * v [0]; 
        nod[0][1] = proj[1] + cos (angle) * dif [1] + r * sin (angle) * v [1];
        nod[0][2] = proj[2] + cos (angle) * dif [2] + r * sin (angle) * v [2];

	/* and extend along the normalize curve derivative */
	u [0] = - sin (angle) * dif [0] + r * cos (angle) * v [0];
	u [1] = - sin (angle) * dif [1] + r * cos (angle) * v [1];
	u [2] = - sin (angle) * dif [2] + r * cos (angle) * v [2];
	NORMALIZE (u);

	/* by the remaining length */
	len = dot - angle * rmin;
	ADDMUL (nod[0], len, u, nod[0]);
      }

      COPY (nod[0], sel[0]);
    }
  }

  /* update face normals */
  MESH_Update (msh, NULL, NULL, NULL);

  Py_INCREF (shape);
  return (PyObject*)shape;
}


/* get object by label */
static PyObject* lng_BYLABEL (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "kind", "label");
  PyObject *kind, *label, *out;
  lng_SOLFEC *solfec;

  PARSEKEYS ("OOO", &solfec, &kind, &label);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_string (kind, kwl[1]) && is_string (label, kwl[2]));

  IFIS (kind, "SURFACE_MATERIAL")
  {
    SURFACE_MATERIAL *mat;

    if ((mat = SPSET_Find_Label (solfec->sol->sps, PyString_AsString (label))))
      out = lng_SURFACE_MATERIAL_WRAPPER (mat);
    else Py_RETURN_NONE;
  }
  ELIF (kind, "BULK_MATERIAL")
  {
    BULK_MATERIAL *mat;

    if ((mat = MATSET_Find (solfec->sol->mat, PyString_AsString (label))))
      out = lng_BULK_MATERIAL_WRAPPER (mat);
    else Py_RETURN_NONE;
  }
  ELIF (kind, "FIELD")
  {
    FIELD *fld;

    if ((fld = FISET_Find (solfec->sol->fis, PyString_AsString (label))))
      out = lng_FIELD_WRAPPER (fld);
    else Py_RETURN_NONE;
  }
  ELIF (kind, "BODY")
  {
    BODY *bod;

    if ((bod = DOM_Find_Body (solfec->sol->dom, PyString_AsString (label))))
      out = lng_BODY_WRAPPER (bod);
    else Py_RETURN_NONE;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid object kind");
    return NULL;
  }

  return out;
}

/* calculate mass center */
static PyObject* lng_MASS_CENTER (PyObject *self, PyObject *args, PyObject *kwds)
{
  double v, c [3] = {0, 0, 0}, e [9];
  KEYWORDS ("shape");
  PyObject *shape;
  SHAPE *shp;

  PARSEKEYS ("O", &shape);

  TYPETEST (is_shape_or_body (shape, kwl[0]));

  if (is_body_check (shape)) /* body */
  {
    BODY *bod = ((lng_BODY*)shape)->bod;
    SHAPE_Char (bod->shape, 0, &v, c, e);
  }
  else /* shape */
  {
    shp = create_shape (shape, 0);
    SHAPE_Char (shp, 0, &v, c, e);
    SHAPE_Destroy_Wrapper (shp);
  }

  return Py_BuildValue ("(d, d, d)", c[0], c[1], c[2]);
}

/* exclude a pair of bodies from contact detection */
static PyObject* lng_CONTACT_EXCLUDE_BODIES (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body1", "body2");
  lng_BODY *body1, *body2;
  SOLFEC *sol;

  PARSEKEYS ("OO", &body1, &body2);

  TYPETEST (is_body (body1, kwl[0]) && is_body (body2, kwl[1]));

#if MPI
  if (body1->dom != body2->dom)
  {
    PyErr_SetString (PyExc_ValueError, "Bodies from different domains");
    return NULL;
  }

  sol = body1->dom->solfec;

  AABB_Exclude_Body_Pair (sol->aabb, body1->id, body2->id);
#else
  if (body1->bod->dom != body2->bod->dom)
  {
    PyErr_SetString (PyExc_ValueError, "Bodies from different domains");
    return NULL;
  }

  sol = body1->bod->dom->solfec;

  AABB_Exclude_Body_Pair (sol->aabb, body1->bod->id, body2->bod->id);
#endif

  Py_RETURN_NONE;
}

/* exclude a pair of surfaces from contact detection */
static PyObject* lng_CONTACT_EXCLUDE_SURFACES (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "surf1", "surf2");
  lng_SOLFEC *solfec;
  int surf1, surf2;

  PARSEKEYS ("Oii", &solfec, &surf1, &surf2);

  TYPETEST (is_solfec (solfec, kwl[0]));

  DOM_Exclude_Contact (solfec->sol->dom, surf1, surf2);

  Py_RETURN_NONE;
}

/* set contact sparsification threshold */
static PyObject* lng_CONTACT_SPARSIFY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "threshold", "minarea", "mindist");
  double threshold, minarea, mindist;
  lng_SOLFEC *solfec;

  mindist = GEOMETRIC_EPSILON;
  minarea = 0.0;

  PARSEKEYS ("Od|dd", &solfec, &threshold, &minarea, &mindist);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_non_negative (minarea, kwl[2]) && is_non_negative (mindist, kwl[3]));;

  if (threshold < 0 || threshold > 1.0)
  {
    PyErr_SetString (PyExc_ValueError, "threshold outside of [0, 1] bounds");
    return NULL;
  }

  solfec->sol->dom->threshold = threshold;
  solfec->sol->dom->minarea = minarea;
  solfec->sol->dom->mindist = mindist;

  Py_RETURN_NONE;
}

/* test whether an object is a constraint solver */
static int is_solver (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance (obj, (PyObject*)&lng_GAUSS_SEIDEL_SOLVER_TYPE) &&
        !PyObject_IsInstance (obj, (PyObject*)&lng_PENALTY_SOLVER_TYPE) &&
        !PyObject_IsInstance (obj, (PyObject*)&lng_NEWTON_SOLVER_TYPE) &&
#if WITHSICONOS
        !PyObject_IsInstance (obj, (PyObject*)&lng_SICONOS_SOLVER_TYPE) &&
#endif
        !PyObject_IsInstance (obj, (PyObject*)&lng_TEST_SOLVER_TYPE))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a constraint solver object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* return solver kind */
static int get_solver_kind (PyObject *obj)
{
  if (PyObject_IsInstance (obj, (PyObject*)&lng_GAUSS_SEIDEL_SOLVER_TYPE))
    return GAUSS_SEIDEL_SOLVER;
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_PENALTY_SOLVER_TYPE))
    return PENALTY_SOLVER;
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_NEWTON_SOLVER_TYPE))
    return NEWTON_SOLVER;
#if WITHSICONOS
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_SICONOS_SOLVER_TYPE))
    return SICONOS_SOLVER;
#endif
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_TEST_SOLVER_TYPE))
    return TEST_SOLVER;
  else return -1;
}

/* return solver interface */
static void* get_solver (PyObject *obj)
{
  if (PyObject_IsInstance (obj, (PyObject*)&lng_GAUSS_SEIDEL_SOLVER_TYPE))
    return ((lng_GAUSS_SEIDEL_SOLVER*)obj)->gs;
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_PENALTY_SOLVER_TYPE))
    return ((lng_PENALTY_SOLVER*)obj)->ps;
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_NEWTON_SOLVER_TYPE))
    return ((lng_NEWTON_SOLVER*)obj)->ns;
#if WITHSICONOS
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_SICONOS_SOLVER_TYPE))
    return ((lng_SICONOS_SOLVER*)obj)->si;
#endif
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_TEST_SOLVER_TYPE))
    return ((lng_TEST_SOLVER*)obj)->ts;
  else return NULL;
}

/* run analysis */
static PyObject* lng_RUN (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "solver", "duration");
  PyObject *solver;
  lng_SOLFEC *solfec;
  double duration;
  int error;

  PARSEKEYS ("OOd", &solfec, &solver, &duration);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_solver (solver, kwl[1]) && is_positive (duration, kwl[2]));

  if (solfec->sol->mode == SOLFEC_READ) Py_RETURN_NONE; /* skip READ mode */

#if OPENGL 
  if (!RND_Is_On ()) /* otherwise interactive run is controlled by the viewer */
#endif
  {
    TRY ()
    {
      SOLFEC_Run (solfec->sol,
		  get_solver_kind (solver),
		  get_solver (solver),
		  duration);
    }
    CATCHANY (error)
    {
      PyErr_SetString (PyExc_RuntimeError, errstring (error));
#if MPI
      PyErr_Print ();
      MPI_Abort (MPI_COMM_WORLD, 5000+error);
#endif
      return NULL;
    }
    ENDTRY ()
  }
#if OPENGL
  else /* only map a domain to a specific solver */
  {
    RND_Solver (solfec->sol->dom,
		get_solver_kind (solver),
		get_solver (solver));
  }
#endif

  Py_RETURN_NONE;
}

/* set output frequency */
static PyObject* lng_OUTPUT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "interval", "compression");
  PyObject *compression;
  lng_SOLFEC *solfec;
  double interval;
  PBF_FLG cmp;

  compression = NULL;
  cmp = PBF_OFF;

  PARSEKEYS ("Od|O", &solfec, &interval, &compression);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_non_negative (interval, kwl[1]) && is_string (compression, kwl [2]));

  if (solfec->sol->mode == SOLFEC_READ) Py_RETURN_NONE; /* skip READ mode */

  if (compression)
  {
    IFIS (compression, "OFF")
    {
      cmp = PBF_OFF;
    }
    ELIF (compression, "ON")
    {
      cmp = PBF_ON;
    }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid compression mode");
      return NULL;
    }
  }

  SOLFEC_Output (solfec->sol, interval, cmp);

  Py_RETURN_NONE;
}

/* set scene extents */
static PyObject* lng_EXTENTS (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "extents");
  lng_SOLFEC *solfec;
  PyObject *extents;
  double e [6];

  PARSEKEYS ("OO", &solfec, &extents);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_tuple (extents, kwl[1], 6));

  if (solfec->sol->mode == SOLFEC_READ) Py_RETURN_NONE; /* skip READ mode */

  e [0] = PyFloat_AsDouble (PyTuple_GetItem (extents, 0));
  e [1] = PyFloat_AsDouble (PyTuple_GetItem (extents, 1));
  e [2] = PyFloat_AsDouble (PyTuple_GetItem (extents, 2));
  e [3] = PyFloat_AsDouble (PyTuple_GetItem (extents, 3));
  e [4] = PyFloat_AsDouble (PyTuple_GetItem (extents, 4));
  e [5] = PyFloat_AsDouble (PyTuple_GetItem (extents, 5));

  DOM_Extents (solfec->sol->dom, e);

  Py_RETURN_NONE;
}

/* analysis callback */
static int lng_CALLBACK_callback (SOLFEC *sol, PyObject *data, PyObject *callback)
{
  PyObject *result;
  PyObject *args;
  int ret;

  ret = 1; /* do not stop time stepping run */

  if (PyTuple_Check (data)) args = data;
  else args = Py_BuildValue ("(O)", data);

  result = PyObject_CallObject (callback, args); /* call user callback */

  if (!PyTuple_Check (data)) { Py_DECREF (args); }

  if (result)
  {
    ret = PyInt_AsLong (result);
    Py_DECREF (result);
  }
  else /* error during the Python callback run */
  {
    PyErr_Print (); /* print traceback */
#if MPI
    MPI_Abort (MPI_COMM_WORLD, 6000);
#else
    lngfinalize ();
#endif
    exit (1);
  }

  return ret;
}

/* set analysis callback */
static PyObject* lng_CALLBACK (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "interval", "data", "callback");
  PyObject *data, *callback;
  lng_SOLFEC *solfec;
  double interval;

  data = NULL;

  PARSEKEYS ("OdOO", &solfec, &interval, &data, &callback);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_positive (interval, kwl[1]) &&
            is_callable (callback, kwl[3]));

  if (solfec->sol->mode == SOLFEC_READ) Py_RETURN_NONE; /* skip READ mode */

  SOLFEC_Set_Callback (solfec->sol, interval, data, callback, (SOLFEC_Callback) lng_CALLBACK_callback);
  callback_pair_push (data, callback);

  Py_RETURN_NONE;
}

/* set unphysical penetration depth */
static PyObject* lng_UNPHYSICAL_PENETRATION (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "depth");
  lng_SOLFEC *solfec;
  double depth;

  PARSEKEYS ("Od", &solfec, &depth);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_positive (depth, kwl[1]));

  if (solfec->sol->mode == SOLFEC_READ) Py_RETURN_NONE; /* skip READ mode */

  solfec->sol->dom->depth = -depth; /* make negative */

  Py_RETURN_NONE;
}

/* set geometric epsilon */
static PyObject* lng_GEOMETRIC_EPSILON (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("epsilon");
  double epsilon;

  PARSEKEYS ("d", &epsilon);

  TYPETEST (is_positive (epsilon, kwl[0]));

  GEOMETRIC_EPSILON = epsilon;

  Py_RETURN_NONE;
}

/* enable/dsiable warnings */
static PyObject* lng_WARNINGS (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("state");
  PyObject *state;

  PARSEKEYS ("O", &state);

  TYPETEST (is_string (state, kwl[0]));

  IFIS (state, "ON")
  {
    WARNINGS_ENABLED = 1;
  }
  ELIF (state, "OFF")
  {
    WARNINGS_ENABLED = 0;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Only 'ON' or 'OFF' are valid states");
    return NULL;
  }

  Py_RETURN_NONE;
}

/* initialize state */
static PyObject* lng_INITIALIZE_STATE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "path", "time");
  lng_SOLFEC *solfec;
  PyObject *path;
  double time;

  PARSEKEYS ("OOd", &solfec, &path, &time);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_string (path, kwl[1]));

  if (solfec->sol->mode == SOLFEC_WRITE)
  {
    if (SOLFEC_Initialize_State (solfec->sol, PyString_AsString (path), time) == 0)
    {
      PyErr_SetString (PyExc_RuntimeError, "State initialization has failed");
      return NULL;
    }
  }
  else
  {
    WARNING (0, "INITIALIZE_STATE has been ingnored in the 'READ' mode");
  }

  Py_RETURN_NONE;
}

/* dump local dynamics */
static PyObject* lng_LOCDYN_DUMP (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "path");
  lng_SOLFEC *solfec;
  PyObject *path;
  char *pstr;

  PARSEKEYS ("OO", &solfec, &path);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_string (path, kwl[1]));

  pstr = PyString_AsString (path);

  LOCDYN_Dump (solfec->sol->dom->ldy, pstr);

  Py_RETURN_NONE;
}

/* overlap callback data */
typedef struct 
{
  SET *set;
  BODY *bod;
  double gap;
} OCD;

/* overlap callback */
static void overlap_create (OCD *ocd, BOX *one, BOX *two)
{
  double onepnt [3], twopnt [3], normal [3], gap, area;
  int state, spair [2];

  if (one->body == two->body) return;

  state = gobjcontact (
    CONTACT_DETECT, GOBJ_Pair_Code (one, two),
    one->sgp->shp, one->sgp->gobj,
    two->sgp->shp, two->sgp->gobj,
    onepnt, twopnt, normal, &gap, &area, spair, NULL, NULL);

  if (state && gap <= ocd->gap)
  {
    if (ocd->bod == one->body)
    {
      SET_Insert (NULL, &ocd->set, one->sgp->shp, NULL);
    }
    else
    {
      SET_Insert (NULL, &ocd->set, two->sgp->shp, NULL);
    }
  }
}

/* calculate overlaps */
static PyObject* lng_OVERLAPPING (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("obstacles", "shapes", "not", "gap");
  PyObject *obstacles, *shapes, *notobj;
  MAP *imap, *omap, *jtem;
  SET *item, *add, *sub;
  SHAPE *outshp, *ptr;
  BODY *obs, *shp;
  SGP *sgp, *sgpe;
  AABB *aabb;
  int not, i;
  OCD ocd;

  notobj = NULL;
  not = 0;
  ocd.gap = 0;

  PARSEKEYS ("OO|Od", &obstacles, &shapes, &notobj, &ocd.gap);

  TYPETEST (is_shape (obstacles, kwl[0]) && is_shape (shapes, kwl[1]) &&
            is_string (notobj, kwl [2]) && is_le (ocd.gap, kwl[3], 0));

  if (notobj)
  {
    IFIS (notobj, "NOT")
    {
      not = 1;
    }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Only 'NOT' is valid");
      return NULL;
    }
  }

  ERRMEM (obs = MEM_CALLOC (sizeof (BODY)));
  ERRMEM (shp = MEM_CALLOC (sizeof (BODY)));
  obs->shape = create_shape (obstacles, 0);
  shp->shape = create_shape (shapes, -1); /* empty and simple glue */
  obs->sgp = SGP_Create (obs->shape, &obs->nsgp);
  shp->sgp = SGP_Create (shp->shape, &shp->nsgp);
  aabb = AABB_Create (obs->nsgp + shp->nsgp);
  obs->kind = shp->kind = RIG;
  ocd.bod = shp;
  ocd.set = NULL;
  add = NULL;
  sub = NULL;

  /* sequentially number input shapes */
  for (ptr = shp->shape, i = 0, imap = NULL; ptr; ptr = ptr->next, i ++)
  {
    MAP_Insert (NULL, &imap, ptr, (void*) (long) i, NULL);
  }

  for (sgp = obs->sgp, sgpe = sgp + obs->nsgp; sgp < sgpe; sgp ++)
  {
    AABB_Insert (aabb, obs, sgp->kind, sgp, SGP_Extents_Update (sgp));
  }
  for (sgp = shp->sgp, sgpe = sgp + shp->nsgp; sgp < sgpe; sgp ++)
  {
    AABB_Insert (aabb, shp, sgp->kind, sgp, SGP_Extents_Update (sgp));
  }

  /* detect overlaps */
  AABB_Simple_Detect (aabb, HYBRID, &ocd, (BOX_Overlap_Create) overlap_create);

  if (not)
  {
    for (sgp = shp->sgp, sgpe = sgp + shp->nsgp; sgp < sgpe; sgp ++)
    {
      if (SET_Contains (ocd.set, sgp->shp, NULL)) SET_Insert (NULL, &sub, sgp->shp, NULL);
      else SET_Insert (NULL, &add, sgp->shp, NULL);
      sgp->shp->next = NULL;
    }
  }
  else
  {
    for (sgp = shp->sgp, sgpe = sgp + shp->nsgp; sgp < sgpe; sgp ++)
    {
      if (SET_Contains (ocd.set, sgp->shp, NULL)) SET_Insert (NULL, &add, sgp->shp, NULL);
      else SET_Insert (NULL, &sub, sgp->shp, NULL);
      sgp->shp->next = NULL;
    }
  }

  /* create output map of sequentially numbered shapes */
  for (item = SET_First (add), omap = NULL; item; item = SET_Next (item))
  {
    jtem = MAP_Find_Node (imap, item->data, NULL);
    ASSERT_DEBUG (jtem, "Invalid shape mapping");
    MAP_Insert (NULL, &omap, jtem->data, jtem->key, NULL); /* use initial numbering */
  }


  /* create output shape list using a subseqnence of the initial sequence of shapes; this assures repeteability
   * of the output in parallel, regardless of the randomized nature of the output of the AABB_Simple_Detect */
  for (outshp = NULL, jtem = MAP_First (omap); jtem; jtem = MAP_Next (jtem))
  {
    ptr = jtem->data;
    ptr->next = outshp;
    outshp = ptr;
  }

  /* destroy remaining shapes */
  for (item = SET_First (sub); item; item = SET_Next (item))
  {
    SHAPE_Destroy (item->data);
  }

  free (obs->sgp);
  free (shp->sgp);
  AABB_Destroy (aabb);
  SET_Free (NULL, &add);
  SET_Free (NULL, &sub);
  SET_Free (NULL, &ocd.set);
  MAP_Free (NULL, &imap);
  MAP_Free (NULL, &omap);
  SHAPE_Destroy_Wrapper (obs->shape);
  free (obs);
  free (shp);

  return shape_to_list (outshp);
}

/* MBFCP export */
static PyObject* lng_MBFCP_EXPORT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "path");
  lng_SOLFEC *solfec;
  PyObject *path;
  char *outpath;
  FILE *out;
#if MPI
  int rank;
#endif

  PARSEKEYS ("OO", &solfec, &path);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_string (path, kwl[1]));

#if MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  ERRMEM (outpath = malloc (strlen (PyString_AsString (path)) + 64));
  sprintf (outpath, "%s.%d", PyString_AsString (path), rank);
#else
  outpath = PyString_AsString (path);
#endif

  if (!(out = fopen (outpath, "w")))
  {
    PyErr_SetString (PyExc_RuntimeError, "File open failed");
    return NULL;
  }

  SOLFEC_2_MBFCP (solfec->sol, out);

#if MPI
  free (outpath);
#endif

  Py_RETURN_NONE;
}

/* non-Solfec arguments */
static PyObject* lng_NON_SOLFEC_ARGV (PyObject *self, PyObject *args, PyObject *kwds)
{
  char **argv;
  int argc;

  if ((argv = NON_SOLFEC_ARGV (&argc)))
  {
    PyObject *list, *obj;
    int n;

    if (!(list = PyList_New (argc))) return NULL;

    for (n = 0; n < argc; n ++)
    {
      obj = PyString_FromString (argv [n]);
      PyList_SetItem (list, n, obj);
    }

    return list;
  }
  else Py_RETURN_NONE;
}

/* read modal data */
static MX* read_modal_analysis (const char *path, double **v)
{
  int m, n;
  MX *E;

#if HDF5
  PBF *f;

  f = PBF_Read (path);
  if (!f) return NULL;
  PBF_Int (f, &m, 1);
  PBF_Int (f, &n, 1);
  E = MX_Create (MXDENSE, m, n, NULL, NULL);
  ERRMEM (*v = malloc (sizeof (double [n])));
  PBF_Double (f, *v, n);
  PBF_Double (f, E->x, n*m);
  PBF_Close (f);
#else
  FILE *f;
  XDR xdr;

  f = fopen (path, "r");
  if (!f) return NULL;
  xdrstdio_create (&xdr, f, XDR_DECODE);
  ASSERT (xdr_int (&xdr, &m), ERR_FILE_READ);
  ASSERT (xdr_int (&xdr, &n), ERR_FILE_READ);
  E = MX_Create (MXDENSE, m, n, NULL, NULL);
  ERRMEM (*v = malloc (sizeof (double [n])));
  ASSERT (xdr_vector (&xdr, (char*)(*v), n, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_READ);
  ASSERT (xdr_vector (&xdr, (char*)E->x, n*m, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_READ);
  xdr_destroy (&xdr);
  fclose (f);
#endif

  return E;
}

/* write modal data */
static void write_modal_analysis (MX *E, double *v, char *path)
{
#if POSIX
  int i, l = strlen (path);

  for (i = 0; i < l; i ++) /* create all directories on the way */
  {
    if (path [i] == '/')
    {
       path [i] = '\0';
       mkdir (path, 0777); /* POSIX */
       path [i] = '/';
    }
  }
#endif

#if HDF5
  PBF *f;

  f = PBF_Write (path, PBF_OFF, PBF_OFF);
  PBF_Int (f, &E->m, 1);
  PBF_Int (f, &E->n, 1);
  PBF_Double (f, v, E->n);
  PBF_Double (f, E->x, E->n*E->m);
  PBF_Close (f);
#else
  FILE *f;
  XDR xdr;

  ASSERT (f = fopen (path, "w"), ERR_FILE_OPEN);
  xdrstdio_create (&xdr, f, XDR_ENCODE);
  ASSERT (xdr_int (&xdr, &E->m), ERR_FILE_WRITE);
  ASSERT (xdr_int (&xdr, &E->n), ERR_FILE_WRITE);
  ASSERT (xdr_vector (&xdr, (char*)v, E->n, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_WRITE);
  ASSERT (xdr_vector (&xdr, (char*)E->x, E->n*E->m, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_WRITE);
  xdr_destroy (&xdr);
  fclose (f);
#endif
}

/* model analysis of FEM bodies */
static PyObject* lng_MODAL_ANALYSIS (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "num", "path", "abstol", "maxiter", "verbose");
  PyObject *path, *val, *vec, *verbose;
  int num, i, maxiter, vrb;
  double *v, abstol;
  lng_BODY *body;
  MX *V;

  vrb = 0;
  abstol = 1E-11;
  maxiter = 100;
  verbose = NULL;
  V = NULL;
  v = NULL;

  PARSEKEYS ("OiO|diO", &body, &num, &path, &abstol, &maxiter, &verbose);

  TYPETEST (is_body (body, kwl[0]) && is_positive (num, kwl [1]) && is_string (path, kwl [2]) &&
      is_positive (abstol, kwl [3]) && is_positive (maxiter, kwl [4]) && is_string (verbose, kwl [5]));

#if MPI && !LOCAL_BODIES /* avoid writing to the same file from several processes */
  int counter, rank;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  for (counter = 0; counter < 2; counter ++)
  {
  if ((rank == 0 && counter == 0) || /* let the rank 0 process go first, while others skip at first */
      (rank > 0 && counter > 0)) /* at second run, let the rank 0 process skip, while others enter */
  {
#endif

  V = read_modal_analysis (PyString_AsString (path), &v);
  if (V && v) /* read without refering to a body et all (precomputed results, parallel use, etc.) */
  {
    if (V->n != num) /* recompute anew */
    {
      MX_Destroy (V);
      free (v);
      V = NULL;
      v = NULL;
    }
  }

#if MPI && LOCAL_BODIES
  if (IS_HERE (body))
  {
#endif

  if (body->bod->kind != FEM)
  {
    PyErr_SetString (PyExc_RuntimeError, "Input body must be of Finite Element kind");
    return NULL;
  }

  if (verbose)
  {
    IFIS (verbose, "ON")
    {
      vrb = 1;
    }
    ELIF (verbose, "OFF")
    {
      vrb = 0;
    }
    ELSE
    {
      PyErr_SetString (PyExc_RuntimeError, "The verbose value can be either 'ON' or 'OFF'");
      return NULL;
    }
  }

  if (!(V && v))
  {
    ERRMEM (v = malloc (sizeof (double [num])));

    V = FEM_Modal_Analysis (body->bod, num, abstol, maxiter, vrb, v);

    if (V)
    {
      write_modal_analysis (V, v, PyString_AsString (path));
    }
    else
    {
      free (v);
      PyErr_SetString (PyExc_RuntimeError, "Eigenvalue solver has failed!");
      return NULL;
    }
  }

  body->bod->evec = V;
  body->bod->eval = v;

#if MPI && LOCAL_BODIES
  }
#endif

#if MPI && !LOCAL_BODIES /* avoid writing to the same file from several processes  -- end */
  }
  MPI_Barrier (MPI_COMM_WORLD); /* all processes meet here twice */
  }
#endif

  if (V && v)
  {
    ERRMEM (val = PyList_New (V->n));
    ERRMEM (vec = PyList_New (V->nzmax));

    for (i = 0; i < V->n; i ++) PyList_SetItem (val, i, PyFloat_FromDouble (v [i]));
    for (i = 0; i < V->nzmax; i ++) PyList_SetItem (vec, i, PyFloat_FromDouble (V->x [i]));

    return Py_BuildValue ("(O, O)", val, vec);
  }
  else Py_RETURN_NONE;
}

/* export body matrices in MatrixMarket format */
static PyObject* lng_BODY_MM_EXPORT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "pathM", "pathK", "spdM", "spdK");
  PyObject *pathM, *pathK, *spdM, *spdK;
  short spd_M, spd_K;
  lng_BODY *body;

  spd_M = 1;
  spd_K = 1;
  spdM = NULL;
  spdK = NULL;

  PARSEKEYS ("OOO|OO", &body, &pathM, &pathK, &spdM, &spdK);

  TYPETEST (is_body (body, kwl[0]) && is_string (pathM, kwl [1]) &&
    is_string (pathK, kwl [2]) && is_string (spdM, kwl [3]) && is_string (spdK, kwl [4]));

#if MPI && LOCAL_BODIES
  if (IS_HERE (body))
  {
#endif

  if (body->bod->kind != FEM)
  {
    PyErr_SetString (PyExc_RuntimeError, "Input body must be of Finite Element kind");
    return NULL;
  }

  if (spdM)
  {
    IFIS (spdM, "ON")
    {
      spd_M = 1;
    }
    ELIF (spdM, "OFF")
    {
      spd_M = 0;
    }
    ELSE
    {
      PyErr_SetString (PyExc_RuntimeError, "The spdM value can be either 'ON' or 'OFF'");
      return NULL;
    }
  }

if (spdK)
  {
    IFIS (spdK, "ON")
    {
      spd_K = 1;
    }
    ELIF (spdK, "OFF")
    {
      spd_K = 0;
    }
    ELSE
    {
      PyErr_SetString (PyExc_RuntimeError, "The spdK value can be either 'ON' or 'OFF'");
      return NULL;
    }
  }

  FEM_MatrixMarket_M_K (body->bod, spd_M, PyString_AsString (pathM), spd_K, PyString_AsString (pathK));

#if MPI && LOCAL_BODIES
  }
#endif

  Py_RETURN_NONE;
}

/* add display point */
static PyObject* lng_DISPLAY_POINT (PyObject *self, PyObject *args, PyObject *kwds)
{
#if !MPI
  KEYWORDS ("body", "point", "label");
  PyObject *point, *label;
  DISPLAY_POINT *dp;
  lng_BODY *body;
  int n;

  label = NULL;

  PARSEKEYS ("OO|O", &body, &point, &label);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3) && is_string (label, kwl[2]));

  ERRMEM (dp = MEM_CALLOC (sizeof (DISPLAY_POINT)));

  dp->X [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  dp->X [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  dp->X [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  if ((n = SHAPE_Sgp (body->bod->sgp, body->bod->nsgp, dp->X)) < 0)
  {
    PyErr_SetString (PyExc_ValueError, "Point outside of the body");
    return NULL;
  }

  COPY (dp->X, dp->x);

  dp->sgp = &body->bod->sgp [n];

  if (label)
  {
    char *string = PyString_AsString (label);
    ERRMEM (dp->label = malloc (strlen (string) + 1));
    strcpy (dp->label, string);
  }

  SET_Insert (NULL, &body->bod->displaypoints, dp, NULL);
#endif

  Py_RETURN_NONE;
}

/* export data for fracture analysis in Yaffems */
static PyObject* lng_FRACTURE_EXPORT_YAFFEMS (PyObject *self, PyObject *args, PyObject *kwds)
{
#if !MPI
  KEYWORDS ("body", "path", "volume", "quality");
  double volume, quality;
  char text [1024];
  PyObject *path;
  lng_BODY *body;
  FILE *output;
  int num;

  volume = DBL_MAX;
  quality = 1.3;

  PARSEKEYS ("OO|dd", &body, &path, &volume, &quality);

  TYPETEST (is_body (body, kwl[0]) && is_string (path, kwl[1]) && 
     is_non_negative (volume, kwl[2]) && is_gt (quality, kwl[3], 1.0));

  if ((body->bod->flags & BODY_CHECK_FRACTURE) == 0)
  {
    PyErr_SetString (PyExc_ValueError, "Not a FEM body with fracturecheck enabled!");
    return NULL;
  }

  snprintf (text, 1024, "%s.vtk", PyString_AsString (path));

  output = fopen (text, "w");

  if (output == NULL)
  {
    PyErr_SetString (PyExc_ValueError, "Failed to create the output file!");
    return NULL;
  }

  num = Fracture_Export_Yaffems (body->bod, volume, quality, output);

  fclose (output);

  return PyInt_FromLong (num);
#endif

  Py_RETURN_NONE;
}

/* simulation duration */
static PyObject* lng_DURATION (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec");
  lng_SOLFEC *solfec;

  PARSEKEYS ("O", &solfec);

  TYPETEST (is_solfec (solfec, kwl[0]));

  if (solfec->sol->mode == SOLFEC_WRITE) return PyFloat_FromDouble (solfec->sol->dom->time);
  else
  {
    double start, end;
    SOLFEC_Time_Limits (solfec->sol, &start, &end);
    return Py_BuildValue ("(d, d)", start, end);
  }
}

/* skip forward */
static PyObject* lng_FORWARD (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "steps");
  lng_SOLFEC *solfec;
  int steps;

  PARSEKEYS ("Oi", &solfec, &steps);

  TYPETEST (is_solfec (solfec, kwl[0]));

  if (solfec->sol->mode == SOLFEC_READ) SOLFEC_Forward (solfec->sol, steps); 

  Py_RETURN_NONE;
}

/* skip backward */
static PyObject* lng_BACKWARD (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "steps");
  lng_SOLFEC *solfec;
  int steps;

  PARSEKEYS ("Oi", &solfec, &steps);

  TYPETEST (is_solfec (solfec, kwl[0]));

  if (solfec->sol->mode == SOLFEC_READ) SOLFEC_Backward (solfec->sol, steps); 

  Py_RETURN_NONE;
}

/* seek to time */
static PyObject* lng_SEEK (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "time");
  lng_SOLFEC *solfec;
  double time;

  PARSEKEYS ("Od", &solfec, &time);

  TYPETEST (is_solfec (solfec, kwl[0]));

  if (solfec->sol->mode == SOLFEC_READ) SOLFEC_Seek_To (solfec->sol, time); 

  Py_RETURN_NONE;
}

/* get displacement */
static PyObject* lng_DISPLACEMENT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point");
  double p [3], x [3];
  lng_BODY *body;
  PyObject *point;
  int error;

  PARSEKEYS ("OO", &body, &point);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3));

#if MPI
  if (IS_HERE (body))
  {
#endif

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  TRY ()
    BODY_Point_Values (body->bod, p, VALUE_DISPLACEMENT, x);
  CATCHANY (error)
  {
    PyErr_SetString (PyExc_RuntimeError, errstring (error));
#if MPI
    PyErr_Print ();
    MPI_Abort (MPI_COMM_WORLD, 2000+error);
#endif
    return NULL;
  }
  ENDTRY ()

  return Py_BuildValue ("(d, d, d)", x[0], x[1], x[2]);

#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

/* get velocity */
static PyObject* lng_VELOCITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point");
  double p [3], x [3];
  lng_BODY *body;
  PyObject *point;
  int error;

  PARSEKEYS ("OO", &body, &point);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3));

#if MPI
  if (IS_HERE (body))
  {
#endif

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  TRY ()
    BODY_Point_Values (body->bod, p, VALUE_VELOCITY, x);
  CATCHANY (error)
  {
    PyErr_SetString (PyExc_RuntimeError, errstring (error));
#if MPI
    PyErr_Print ();
    MPI_Abort (MPI_COMM_WORLD, 2000+error);
#endif
    return NULL;
  }
  ENDTRY ()



  return Py_BuildValue ("(d, d, d)", x[0], x[1], x[2]);

#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

/* get stress */
static PyObject* lng_STRESS (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point");
  double p [3], x [7];
  lng_BODY *body;
  PyObject *point;
  int error;

  PARSEKEYS ("OO", &body, &point);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3));

#if MPI
  if (IS_HERE (body))
  {
#endif

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  TRY ()
    BODY_Point_Values (body->bod, p, VALUE_STRESS_AND_MISES, x);
  CATCHANY (error)
  {
    PyErr_SetString (PyExc_RuntimeError, errstring (error));
#if MPI
    PyErr_Print ();
    MPI_Abort (MPI_COMM_WORLD, 2000+error);
#endif
    return NULL;
  }
  ENDTRY ()



  return Py_BuildValue ("(d, d, d, d, d, d, d)", x[0], x[1], x[2], x[3], x[4], x[5], x[6]);

#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

/* check whether an object is a SOLFEC, a BODY or a list of bodies */
static int is_solfec_or_body_or_list_of_bodies (PyObject *obj, char *var)
{
  int i, l;
  char buf [BUFLEN];
  PyObject *item;

  if (!obj) return 1;
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_SOLFEC_TYPE)) return 1;
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_BODY_TYPE)) return 1;
  else if (PyList_Check (obj))
  {
    l = PyList_Size (obj);
    for (i = 0; i < l; i ++)
    {
      item = PyList_GetItem (obj, i);
      if (!PyObject_IsInstance (item, (PyObject*)&lng_BODY_TYPE)) goto err;
#if MPI
      if (!IS_HERE ((lng_BODY*)item))
      {
	sprintf (buf, "One of the bodies in '%s' is not on the current processor", var);
	PyErr_SetString (PyExc_TypeError, buf);
	return 0;
      }
#endif
    }
    return 1;
  }
 
err: 
  sprintf (buf, "'%s' must be a SOLFEC a BODY object or a list of BODY objects", var);
  PyErr_SetString (PyExc_TypeError, buf);
  return 0;
}

/* conver object to a body set */
static SET* object_to_body_set (PyObject *obj, MEM *setmem, SOLFEC *sol)
{
  SET *ret;
  int i, l;

  ret = NULL;

  ASSERT_DEBUG (obj, "Invalid NULL object argument to object_to_body_set");

  if (PyObject_IsInstance (obj, (PyObject*)&lng_BODY_TYPE))
  {
    lng_BODY *body = (lng_BODY*)obj;
    SET_Insert (setmem, &ret, body->bod, NULL);
  }
  else if (PyList_Check (obj))
  {
    l = PyList_Size (obj);
    for (i = 0; i < l; i ++)
    {
      lng_BODY *body = (lng_BODY*)PyList_GetItem (obj, i);
      SET_Insert (setmem, &ret, body->bod, NULL);
    }
  }

  return ret;
}

/* render a list of bodies */
static PyObject* lng_RENDER (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "object"); // a body or a list of bodies
  SET *bodies;
  lng_SOLFEC *solfec;
  PyObject *object;
  MEM setmem;
  
  object = NULL;
  
  PARSEKEYS ("O|O", &solfec, &object);
  
  TYPETEST (is_solfec (solfec, kwl[0]) && is_solfec_or_body_or_list_of_bodies (object, kwl[1]));
  
  if (object)
  {
    MEM_Init (&setmem, sizeof (SET), 128);
    bodies = object_to_body_set (object, &setmem, solfec->sol);
  }
  #if OPENGL
  if (RND_Is_On ()) select_id (bodies);
  #endif
  Py_RETURN_NONE;
}


/* energy */
static PyObject* lng_ENERGY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "object");
  double energy [BODY_ENERGY_SPACE];
  SET *bodies, *item;
  lng_SOLFEC *solfec;
  PyObject *object;
  MEM setmem;
  BODY *bod;
  int i;

  object = NULL;

  memset (energy, 0, sizeof (energy));

  PARSEKEYS ("O|O", &solfec, &object);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_solfec_or_body_or_list_of_bodies (object, kwl[1]));

  if (object)
  {
    MEM_Init (&setmem, sizeof (SET), 128);
    bodies = object_to_body_set (object, &setmem, solfec->sol);
    for (item = SET_First (bodies); item; item = SET_Next (item))
    {
      bod = item->data;
      for (i = 0; i < BODY_ENERGY_SIZE(bod); i ++) energy [i] += bod->energy [i];
    }
    MEM_Release (&setmem);
  }
  else
  {
    for (bod = solfec->sol->dom->bod; bod; bod = bod->next)
    {
      for (i = 0; i < BODY_ENERGY_SIZE(bod); i ++) energy [i] += bod->energy [i];
    }
  }

  return Py_BuildValue ("(d, d, d, d, d)", energy [0], energy [1], energy [2], energy [3], energy [4]);
}

/* timing */
static PyObject* lng_TIMING (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "kind");
  lng_SOLFEC *solfec;
  PyObject *kind;
  char *label;

  PARSEKEYS ("OO", &solfec, &kind);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_string (kind, kwl[1]));

  label = PyString_AsString (kind);

  if (!SOLFEC_Has_Timer (solfec->sol, label))
  {
    PyErr_SetString (PyExc_ValueError, "Invalid timing kind");
    return NULL;
  }

  return PyFloat_FromDouble (SOLFEC_Timing (solfec->sol, label));
}

/* parse constraint related history item */
static int parse_constraint_history_item (SHI *shi, PyObject *entity, double *direction, int surf1, int surf2)
{
  if (direction)
  {
    COPY (direction, shi->vector);
  }
  else
  {
    SET (shi->vector, 0.0);
  }

  shi->surf1 = surf1;
  shi->surf2 = surf2;

  IFIS (entity, "GAP")
  {
    shi->item = CONSTRAINT_VALUE;
    shi->index = CONSTRAINT_GAP;
    shi->op = OP_MIN;
  }
  ELIF (entity, "R")
  {
    shi->item = CONSTRAINT_VALUE;
    shi->index = CONSTRAINT_R;
    shi->op = OP_SUM;
  }
  ELIF (entity, "CR")
  {
    shi->item = CONSTRAINT_VALUE;
    shi->index = CONSTRAINT_R;
    shi->contacts_only = 1;
    shi->op = OP_SUM;
  }
  ELIF (entity, "U")
  {
    shi->item = CONSTRAINT_VALUE;
    shi->index = CONSTRAINT_U;
    shi->op = OP_AVG;
  }
  ELIF (entity, "CU")
  {
    shi->item = CONSTRAINT_VALUE;
    shi->index = CONSTRAINT_U;
    shi->contacts_only = 1;
    shi->op = OP_AVG;
  }
  ELSE return 0;

  return 1;
}

/* parse single item of hitory items list */
static int parse_history_item (PyObject *obj, MEM *setmem, SOLFEC *sol, SHI *shi)
{
  if (PyTuple_Check (obj))
  {
    if (PyTuple_Size (obj) == 3)
    {
      PyObject *point, *entity;
      lng_BODY *body;

      body = (lng_BODY*)PyTuple_GetItem (obj, 0);
      point = PyTuple_GetItem (obj, 1);
      entity = PyTuple_GetItem (obj, 2);

      if (!(is_body (body, "body") && is_tuple (point, "point", 3) && is_string (entity, "entity"))) return 0;

      shi->bod = body->bod;
      shi->point [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
      shi->point [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
      shi->point [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

      IFIS (entity, "CX")
      {
	shi->entity = VALUE_COORD;
	shi->index = 0;
      }
      ELIF (entity, "CY")
      {
	shi->entity = VALUE_COORD;
	shi->index = 1;
      }
      ELIF (entity, "CZ")
      {
	shi->entity = VALUE_COORD;
	shi->index = 2;
      }
      ELIF (entity, "DX")
      {
	shi->entity = VALUE_DISPLACEMENT;
	shi->index = 0;
      }
      ELIF (entity, "DY")
      {
	shi->entity = VALUE_DISPLACEMENT;
	shi->index = 1;
      }
      ELIF (entity, "DZ")
      {
	shi->entity = VALUE_DISPLACEMENT;
	shi->index = 2;
      }
      ELIF (entity, "VX")
      {
	shi->entity = VALUE_VELOCITY;
	shi->index = 0;
      }
      ELIF (entity, "VY")
      {
	shi->entity = VALUE_VELOCITY;
	shi->index = 1;
      }
      ELIF (entity, "VZ")
      {
	shi->entity = VALUE_VELOCITY;
	shi->index = 2;
      }
      ELIF (entity, "SX")
      {
	shi->entity = VALUE_STRESS;
	shi->index = 0;
      }
      ELIF (entity, "SY")
      {
	shi->entity = VALUE_STRESS;
	shi->index = 1;
      }
      ELIF (entity, "SZ")
      {
	shi->entity = VALUE_STRESS;
	shi->index = 2;
      }
      ELIF (entity, "SXY")
      {
	shi->entity = VALUE_STRESS;
	shi->index = 3;
      }
      ELIF (entity, "SXZ")
      {
	shi->entity = VALUE_STRESS;
	shi->index = 4;
      }
      ELIF (entity, "SYZ")
      {
	shi->entity = VALUE_STRESS;
	shi->index = 5;
      }
      ELIF (entity, "MISES")
      {
	shi->entity = VALUE_MISES;
	shi->index = 0;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid entity kind");
	return 0;
      }

      shi->item = BODY_ENTITY;
    }
    else if (PyTuple_Size (obj) == 2)
    {
      PyObject *object, *kind;

      object = PyTuple_GetItem (obj, 0);
      kind = PyTuple_GetItem (obj, 1);

      if (!(is_solfec_or_body_or_list_of_bodies (object, "object") && is_string (kind, "kind"))) return 0;

      IFIS (kind, "KINETIC")
      {
	shi->index = KINETIC;
        shi->item = ENERGY_VALUE;
      }
      ELIF (kind, "INTERNAL")
      {
	shi->index = INTERNAL;
        shi->item = ENERGY_VALUE;
      }
      ELIF (kind, "EXTERNAL")
      {
	shi->index = EXTERNAL;
        shi->item = ENERGY_VALUE;
      }
      ELIF (kind, "CONTACT")
      {
	shi->index = CONTWORK;
        shi->item = ENERGY_VALUE;
      }
      ELIF (kind, "FRICTION")
      {
	shi->index = FRICWORK;
        shi->item = ENERGY_VALUE;
      }
      ELSE
      {
	if (!parse_constraint_history_item (shi, kind, NULL, INT_MAX, INT_MAX))
	{
	  char text [BUFLEN];
	  snprintf (text, BUFLEN, "Invalid entity kind: %s", PyString_AsString (kind));
	  PyErr_SetString (PyExc_ValueError, text);
	  return 0;
	}
      }

      shi->bodies = object_to_body_set (object, setmem, sol);
    }
    else if (PyTuple_Size (obj) == 4)
    {
      PyObject *object, *direction, *pair, *entity, *surf1, *surf2;
      double d [3];
      int s [2];

      object = PyTuple_GetItem (obj, 0);
      direction = PyTuple_GetItem (obj, 1);
      pair = PyTuple_GetItem (obj, 2);
      entity = PyTuple_GetItem (obj, 3);

      if (!(is_solfec_or_body_or_list_of_bodies (object, "object") && is_string (entity, "entity"))) return 0;
      if (direction != Py_None && !is_tuple (direction, "direction", 3)) return 0;
      if (pair != Py_None && !is_tuple (pair, "pair", 2)) return 0;

      if (pair != Py_None)
      {
	surf1 = PyTuple_GetItem (pair, 0); surf2 = PyTuple_GetItem (pair, 1);
	if (!(is_number (surf1, "surf1") && is_number (surf2, "surf2"))) return 0;
	s [0] = PyInt_AsLong (surf1);
	s [1] = PyInt_AsLong (surf2);
      }
      else  s [0] = s [1] = INT_MAX;

      if (direction != Py_None)
      {
	for (int i = 0; i < 3; i ++)
	  d [i] = PyFloat_AsDouble (PyTuple_GetItem (direction, i));
      }
      else SET (d, 0);

      if (!parse_constraint_history_item (shi, entity, d, s[0], s[1]))
      {
	char text [BUFLEN];
	snprintf (text, BUFLEN, "Invalid entity kind: %s", PyString_AsString (entity));
	PyErr_SetString (PyExc_ValueError, text);
	return 0;
      }

      shi->bodies = object_to_body_set (object, setmem, sol);
    }
    else
    {
       PyErr_SetString (PyExc_ValueError, "Invalid tuple size");
       return 0;
    }
  }
  else if (PyString_Check (obj))
  {
    shi->label = PyString_AsString (obj);

    IFIS (obj, "TIMINT") { shi->item = TIMING_VALUE; }
    ELIF (obj, "CONUPD") { shi->item = TIMING_VALUE; }
    ELIF (obj, "CONDET") { shi->item = TIMING_VALUE; }
    ELIF (obj, "LOCDYN") { shi->item = TIMING_VALUE; }
    ELIF (obj, "CONSOL") { shi->item = TIMING_VALUE; }
    ELIF (obj, "PARBAL") { shi->item = TIMING_VALUE; }
    ELIF (obj, "GSINIT") { shi->item = TIMING_VALUE; }
    ELIF (obj, "GSRUN") { shi->item = TIMING_VALUE; }
    ELIF (obj, "GSCOM") { shi->item = TIMING_VALUE; }
    ELIF (obj, "GSMCOM") { shi->item = TIMING_VALUE; }
    ELIF (obj, "LININIT") { shi->item = TIMING_VALUE; }
    ELIF (obj, "LINUPD") { shi->item = TIMING_VALUE; }
    ELIF (obj, "LINMV") { shi->item = TIMING_VALUE; }
    ELIF (obj, "LINPRE") { shi->item = TIMING_VALUE; }
    ELIF (obj, "LINRUN") { shi->item = TIMING_VALUE; }
    ELIF (obj, "LINCOM") { shi->item = TIMING_VALUE; }
    ELIF (obj, "STEP") { shi->item = LABELED_DOUBLE; shi->op = OP_MIN; }
    ELIF (obj, "CONS") { shi->item = LABELED_INT; shi->op = OP_SUM; }
    ELIF (obj, "BODS") { shi->item = LABELED_INT; shi->op = OP_SUM; }
    ELIF (obj, "NEWCONS") { shi->item = LABELED_INT; shi->op = OP_SUM; }
    ELIF (obj, "NEWBODS") { shi->item = LABELED_INT; shi->op = OP_SUM; }
    ELIF (obj, "GSITERS") { shi->item = LABELED_INT; shi->op = OP_MAX; }
    ELIF (obj, "GSCOLORS") { shi->item = LABELED_INT; shi->op = OP_MAX; }
    ELIF (obj, "GSBOT") { shi->item = LABELED_INT; shi->op = OP_SUM; }
    ELIF (obj, "GSMID") { shi->item = LABELED_INT; shi->op = OP_SUM; }
    ELIF (obj, "GSTOP") { shi->item = LABELED_INT; shi->op = OP_SUM; }
    ELIF (obj, "GSINN") { shi->item = LABELED_INT; shi->op = OP_SUM; }
    ELIF (obj, "MERIT") { shi->item = LABELED_DOUBLE; shi->op = OP_MAX; }
    ELIF (obj, "NTITERS") { shi->item = LABELED_INT; shi->op = OP_MAX; }
    ELIF (obj, "BSITERS") { shi->item = LABELED_INT; shi->op = OP_MAX; }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid string value");
      return 0;
    }
  }
  else
  {
    PyErr_SetString (PyExc_ValueError, "Invalid history item");
    return 0;
  }

  return 1;
}

/* parse history items list */
static SHI* parse_history_items (PyObject *list, MEM *setmem, SOLFEC *sol, int *nshi)
{
  SHI *shi;
  int i;

  if (PyList_Check (list))
  {
    *nshi = PyList_Size (list);
    ERRMEM (shi = MEM_CALLOC (sizeof (SHI [*nshi])));
    for (i = 0; i < *nshi; i ++)
    {
      if (!parse_history_item (PyList_GetItem (list, i), setmem, sol, &shi [i]))
      {
        free (shi);
	return NULL;
      }
    }
  }
  else
  {
    ERRMEM (shi = MEM_CALLOC (sizeof (SHI)));
    *nshi = 1;
    if (!parse_history_item (list, setmem, sol, shi))
    {
      free (shi);
      return NULL;
    }
  }

  return shi;
}

/* history of an entity */
static PyObject* lng_HISTORY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "list", "t0", "t1", "skip", "progress");
  double start, end, t0, t1, *time;
  PyObject *list, *tuple, *vals, *progress;
  int skip, nshi, i, j, size;
  lng_SOLFEC *solfec;
  SOLFEC *sol;
  MEM setmem;
  SHI *shi;
  int error;

  progress = NULL;
  skip = 1;

  PARSEKEYS ("OOdd|iO", &solfec, &list, &t0, &t1, &skip, &progress);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_positive (skip, kwl [4]) && is_string (progress, kwl[5]));

  sol = solfec->sol;

  if (sol->mode == SOLFEC_WRITE) Py_RETURN_NONE;
  else
  {
    MEM_Init (&setmem, sizeof (SET), 128);

    if (!(shi = parse_history_items (list, &setmem, sol, &nshi)))
    {
      MEM_Release (&setmem);
      return NULL;
    }

    SOLFEC_Time_Limits (sol, &start, &end);

    if (t0 < start) t0 = start;

    if (t1 > end) t1 = end;

    if (progress)
    {
      IFIS (progress, "ON")
      {
	skip = -skip;
      }
      ELIF (progress, "OFF")
      {
	/* do nothing */
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid progress value (neither 'ON' nor 'OFF')");
	return NULL;
      }
    }

    TRY ()
      time = SOLFEC_History (sol, shi, nshi, t0, t1, skip, &size);
    CATCHANY (error)
    {
      PyErr_SetString (PyExc_RuntimeError, errstring (error));
#if MPI
      PyErr_Print ();
      MPI_Abort (MPI_COMM_WORLD, 2000+error);
#endif
      return NULL;
    }
    ENDTRY ()

    ERRMEM (tuple = PyTuple_New (nshi + 1));
    ERRMEM (vals = PyList_New (size));
    for (j = 0; j < size; j ++) PyList_SetItem (vals, j, PyFloat_FromDouble (time [j]));
    PyTuple_SetItem (tuple, 0, vals);
    free (time);

    for (i = 0; i < nshi; i ++)
    {
      ERRMEM (vals = PyList_New (size));
      for (j = 0; j < size; j ++) PyList_SetItem (vals, j, PyFloat_FromDouble (shi [i].history [j]));
      PyTuple_SetItem (tuple, i + 1, vals);
      free (shi [i].history);
    }

    MEM_Release (&setmem);
    free (shi);

    return tuple;
  }
}

static PyMethodDef lng_methods [] =
{
  {"HULL", (PyCFunction)lng_HULL, METH_VARARGS|METH_KEYWORDS, "Create convex hull from a point set"},
  {"MESH2CONVEX", (PyCFunction)lng_MESH2CONVEX, METH_VARARGS|METH_KEYWORDS, "Create list of CONVEX objects from a MESH object"},
  {"HEX", (PyCFunction)lng_HEX, METH_VARARGS|METH_KEYWORDS, "Create mesh of a hexahedral shape"},
  {"TETRAHEDRALIZE", (PyCFunction)lng_TETRAHEDRALIZE, METH_VARARGS|METH_KEYWORDS, "Generate tetrahedral mesh"},
  {"PIPE", (PyCFunction)lng_PIPE, METH_VARARGS|METH_KEYWORDS, "Create mesh of a pipe shape"},
  {"ROUGH_HEX", (PyCFunction)lng_ROUGH_HEX, METH_VARARGS|METH_KEYWORDS, "Create a rough hexahedral mesh containing a given shape"},
  {"FIX_POINT", (PyCFunction)lng_FIX_POINT, METH_VARARGS|METH_KEYWORDS, "Create a fixed point constraint"},
  {"FIX_DIRECTION", (PyCFunction)lng_FIX_DIRECTION, METH_VARARGS|METH_KEYWORDS, "Create a fixed direction constraint"},
  {"SET_DISPLACEMENT", (PyCFunction)lng_SET_DISPLACEMENT, METH_VARARGS|METH_KEYWORDS, "Create a prescribed displacement constraint"},
  {"SET_VELOCITY", (PyCFunction)lng_SET_VELOCITY, METH_VARARGS|METH_KEYWORDS, "Create a prescribed velocity constraint"},
  {"SET_ACCELERATION", (PyCFunction)lng_SET_ACCELERATION, METH_VARARGS|METH_KEYWORDS, "Create a prescribed acceleration constraint"},
  {"PUT_RIGID_LINK", (PyCFunction)lng_PUT_RIGID_LINK, METH_VARARGS|METH_KEYWORDS, "Create a rigid linke constraint"},
  {"PUT_SPRING", (PyCFunction)lng_PUT_SPRING, METH_VARARGS|METH_KEYWORDS, "Create a spring constraint"},
  {"GRAVITY", (PyCFunction)lng_GRAVITY, METH_VARARGS|METH_KEYWORDS, "Set gravity acceleration"},
  {"FORCE", (PyCFunction)lng_FORCE, METH_VARARGS|METH_KEYWORDS, "Apply point force"},
  {"TORQUE", (PyCFunction)lng_TORQUE, METH_VARARGS|METH_KEYWORDS, "Apply point torque"},
  {"PRESSURE", (PyCFunction)lng_PRESSURE, METH_VARARGS|METH_KEYWORDS, "Apply pressure"},
  {"SIMPLIFIED_CRACK", (PyCFunction)lng_SIMPLIFIED_CRACK, METH_VARARGS|METH_KEYWORDS, "Prescribe crack"},
  {"IMBALANCE_TOLERANCE", (PyCFunction)lng_IMBALANCE_TOLERANCE, METH_VARARGS|METH_KEYWORDS, "Adjust parallel imbalance tolerance"},
  {"RANK", (PyCFunction)lng_RANK, METH_NOARGS, "Get current processor rank"},
  {"BARRIER", (PyCFunction)lng_BARRIER, METH_NOARGS, "Set up parallel barrier"},
  {"NCPU", (PyCFunction)lng_NCPU, METH_VARARGS|METH_KEYWORDS, "Get the number of processors"},
  {"HERE", (PyCFunction)lng_HERE, METH_VARARGS|METH_KEYWORDS, "Test whether an object is located on the current processor"},
  {"VIEWER", (PyCFunction)lng_VIEWER, METH_NOARGS, "Test whether the viewer is enabled"},
  {"BODY_CHARS", (PyCFunction)lng_BODY_CHARS, METH_VARARGS|METH_KEYWORDS, "Overwrite body characteristics"},
  {"INITIAL_VELOCITY", (PyCFunction)lng_INITIAL_VELOCITY, METH_VARARGS|METH_KEYWORDS, "Apply initial velocity"},
  {"MATERIAL", (PyCFunction)lng_MATERIAL, METH_VARARGS|METH_KEYWORDS, "Apply bulk material"},
  {"DELETE", (PyCFunction)lng_DELETE, METH_VARARGS|METH_KEYWORDS, "Delete a body or a constraint"},
  {"SCALE", (PyCFunction)lng_SCALE, METH_VARARGS|METH_KEYWORDS, "Scale shape"},
  {"TRANSLATE", (PyCFunction)lng_TRANSLATE, METH_VARARGS|METH_KEYWORDS, "Translate shape"},
  {"ROTATE", (PyCFunction)lng_ROTATE, METH_VARARGS|METH_KEYWORDS, "Rotate shape"},
  {"SPLIT", (PyCFunction)lng_SPLIT, METH_VARARGS|METH_KEYWORDS, "Split shape by plane"},
  {"MESH_SPLIT", (PyCFunction)lng_MESH_SPLIT, METH_VARARGS|METH_KEYWORDS, "Split mesh by node set"},
  {"COPY", (PyCFunction)lng_COPY, METH_VARARGS|METH_KEYWORDS, "Copy shape"},
  {"BEND", (PyCFunction)lng_BEND, METH_VARARGS|METH_KEYWORDS, "Bend shape"},
  {"BYLABEL", (PyCFunction)lng_BYLABEL, METH_VARARGS|METH_KEYWORDS, "Get object by label"},
  {"MASS_CENTER", (PyCFunction)lng_MASS_CENTER, METH_VARARGS|METH_KEYWORDS, "Get mass center"},
  {"CONTACT_EXCLUDE_BODIES", (PyCFunction)lng_CONTACT_EXCLUDE_BODIES, METH_VARARGS|METH_KEYWORDS, "Exclude body pair from contact detection"},
  {"CONTACT_EXCLUDE_SURFACES", (PyCFunction)lng_CONTACT_EXCLUDE_SURFACES, METH_VARARGS|METH_KEYWORDS, "Exclude surface pair from contact detection"},
  {"CONTACT_SPARSIFY", (PyCFunction)lng_CONTACT_SPARSIFY, METH_VARARGS|METH_KEYWORDS, "Adjust contact sparsification"},
  {"RUN", (PyCFunction)lng_RUN, METH_VARARGS|METH_KEYWORDS, "Run analysis"},
  {"OUTPUT", (PyCFunction)lng_OUTPUT, METH_VARARGS|METH_KEYWORDS, "Set data output interval"},
  {"EXTENTS", (PyCFunction)lng_EXTENTS, METH_VARARGS|METH_KEYWORDS, "Set scene extents"},
  {"CALLBACK", (PyCFunction)lng_CALLBACK, METH_VARARGS|METH_KEYWORDS, "Set analysis callback"},
  {"UNPHYSICAL_PENETRATION", (PyCFunction)lng_UNPHYSICAL_PENETRATION, METH_VARARGS|METH_KEYWORDS, "Set unphysical penetration bound"},
  {"GEOMETRIC_EPSILON", (PyCFunction)lng_GEOMETRIC_EPSILON, METH_VARARGS|METH_KEYWORDS, "Set geometric epsilon"},
  {"WARNINGS", (PyCFunction)lng_WARNINGS, METH_VARARGS|METH_KEYWORDS, "Enable or disable warnings"},
  {"INITIALIZE_STATE", (PyCFunction)lng_INITIALIZE_STATE, METH_VARARGS|METH_KEYWORDS, "Initialize Solfec state"},
  {"LOCDYN_DUMP", (PyCFunction)lng_LOCDYN_DUMP, METH_VARARGS|METH_KEYWORDS, "Dump local dynamics"},
  {"OVERLAPPING", (PyCFunction)lng_OVERLAPPING, METH_VARARGS|METH_KEYWORDS, "Detect shapes (not) overlapping obstacles"},
  {"MBFCP_EXPORT", (PyCFunction)lng_MBFCP_EXPORT, METH_VARARGS|METH_KEYWORDS, "Export MBFCP definition"},
  {"NON_SOLFEC_ARGV", (PyCFunction)lng_NON_SOLFEC_ARGV, METH_NOARGS, "Return non-Solfec input arguments"},
  {"MODAL_ANALYSIS", (PyCFunction)lng_MODAL_ANALYSIS, METH_VARARGS|METH_KEYWORDS, "Perform modal analysis of a FEM body"},
  {"BODY_MM_EXPORT", (PyCFunction)lng_BODY_MM_EXPORT, METH_VARARGS|METH_KEYWORDS, "Export MatrixMarket M and K matrices of a FEM body"},
  {"DISPLAY_POINT", (PyCFunction)lng_DISPLAY_POINT, METH_VARARGS|METH_KEYWORDS, "Add display point"},
  {"FRACTURE_EXPORT_YAFFEMS", (PyCFunction)lng_FRACTURE_EXPORT_YAFFEMS, METH_VARARGS|METH_KEYWORDS, "Export fracture data to Yaffems"},
  {"DURATION", (PyCFunction)lng_DURATION, METH_VARARGS|METH_KEYWORDS, "Get analysis duration"},
  {"RENDER", (PyCFunction)lng_RENDER, METH_VARARGS|METH_KEYWORDS, "Render bodies"},
  {"FORWARD", (PyCFunction)lng_FORWARD, METH_VARARGS|METH_KEYWORDS, "Set forward in READ mode"},
  {"BACKWARD", (PyCFunction)lng_BACKWARD, METH_VARARGS|METH_KEYWORDS, "Set backward in READ mode"},
  {"SEEK", (PyCFunction)lng_SEEK, METH_VARARGS|METH_KEYWORDS, "Seek to time in READ mode"},
  {"DISPLACEMENT", (PyCFunction)lng_DISPLACEMENT, METH_VARARGS|METH_KEYWORDS, "Get displacement of a referential point"},
  {"VELOCITY", (PyCFunction)lng_VELOCITY, METH_VARARGS|METH_KEYWORDS, "Get velocity of a referential point"},
  {"STRESS", (PyCFunction)lng_STRESS, METH_VARARGS|METH_KEYWORDS, "Get stress of a referential point"},
  {"ENERGY", (PyCFunction)lng_ENERGY, METH_VARARGS|METH_KEYWORDS, "Get energy"},
  {"TIMING", (PyCFunction)lng_TIMING, METH_VARARGS|METH_KEYWORDS, "Get timing"},
  {"HISTORY", (PyCFunction)lng_HISTORY, METH_VARARGS|METH_KEYWORDS, "Get history of entites"},
  {NULL, 0, 0, NULL}
};

/* 
 * initialization
 */

static void initlng (const char *path)
{
  PyObject *m;
  PyObject *ppath = PyString_FromString(path);

  TYPEINIT (lng_CONVEX_TYPE, lng_CONVEX, "solfec.CONVEX",
    Py_TPFLAGS_DEFAULT, lng_CONVEX_dealloc, lng_CONVEX_new,
    lng_CONVEX_methods, lng_CONVEX_members, lng_CONVEX_getset);

  TYPEINIT (lng_MESH_TYPE, lng_MESH, "solfec.MESH",
    Py_TPFLAGS_DEFAULT, lng_MESH_dealloc, lng_MESH_new,
    lng_MESH_methods, lng_MESH_members, lng_MESH_getset);

  TYPEINIT (lng_SPHERE_TYPE, lng_SPHERE, "solfec.SPHERE",
    Py_TPFLAGS_DEFAULT, lng_SPHERE_dealloc, lng_SPHERE_new,
    lng_SPHERE_methods, lng_SPHERE_members, lng_SPHERE_getset);

  TYPEINIT (lng_ELLIP_TYPE, lng_ELLIP, "solfec.ELLIP",
    Py_TPFLAGS_DEFAULT, lng_ELLIP_dealloc, lng_ELLIP_new,
    lng_ELLIP_methods, lng_ELLIP_members, lng_ELLIP_getset);

  TYPEINIT (lng_SOLFEC_TYPE, lng_SOLFEC, "solfec.SOLFEC",
    Py_TPFLAGS_DEFAULT, lng_SOLFEC_dealloc, lng_SOLFEC_new,
    lng_SOLFEC_methods, lng_SOLFEC_members, lng_SOLFEC_getset);

  TYPEINIT (lng_FIELD_TYPE, lng_FIELD, "solfec.FIELD",
    Py_TPFLAGS_DEFAULT, lng_FIELD_dealloc, lng_FIELD_new,
    lng_FIELD_methods, lng_FIELD_members, lng_FIELD_getset);

  TYPEINIT (lng_SURFACE_MATERIAL_TYPE, lng_SURFACE_MATERIAL, "solfec.SURFACE_MATERIAL",
    Py_TPFLAGS_DEFAULT, lng_SURFACE_MATERIAL_dealloc, lng_SURFACE_MATERIAL_new,
    lng_SURFACE_MATERIAL_methods, lng_SURFACE_MATERIAL_members, lng_SURFACE_MATERIAL_getset);

  TYPEINIT (lng_BULK_MATERIAL_TYPE, lng_BULK_MATERIAL, "solfec.BULK_MATERIAL",
    Py_TPFLAGS_DEFAULT, lng_BULK_MATERIAL_dealloc, lng_BULK_MATERIAL_new,
    lng_BULK_MATERIAL_methods, lng_BULK_MATERIAL_members, lng_BULK_MATERIAL_getset);

  TYPEINIT (lng_BODY_TYPE, lng_BODY, "solfec.BODY",
    Py_TPFLAGS_DEFAULT, lng_BODY_dealloc, lng_BODY_new,
    lng_BODY_methods, lng_BODY_members, lng_BODY_getset);

  TYPEINIT (lng_TIME_SERIES_TYPE, lng_TIME_SERIES, "solfec.TIME_SERIES",
    Py_TPFLAGS_DEFAULT, lng_TIME_SERIES_dealloc, lng_TIME_SERIES_new,
    lng_TIME_SERIES_methods, lng_TIME_SERIES_members, lng_TIME_SERIES_getset);

  TYPEINIT (lng_GAUSS_SEIDEL_SOLVER_TYPE, lng_GAUSS_SEIDEL_SOLVER, "solfec.GAUSS_SEIDEL_SOLVER",
    Py_TPFLAGS_DEFAULT, lng_GAUSS_SEIDEL_SOLVER_dealloc, lng_GAUSS_SEIDEL_SOLVER_new,
    lng_GAUSS_SEIDEL_SOLVER_methods, lng_GAUSS_SEIDEL_SOLVER_members, lng_GAUSS_SEIDEL_SOLVER_getset);

  TYPEINIT (lng_PENALTY_SOLVER_TYPE, lng_PENALTY_SOLVER, "solfec.PENALTY_SOLVER",
    Py_TPFLAGS_DEFAULT, lng_PENALTY_SOLVER_dealloc, lng_PENALTY_SOLVER_new,
    lng_PENALTY_SOLVER_methods, lng_PENALTY_SOLVER_members, lng_PENALTY_SOLVER_getset);

  TYPEINIT (lng_NEWTON_SOLVER_TYPE, lng_NEWTON_SOLVER, "solfec.NEWTON_SOLVER",
    Py_TPFLAGS_DEFAULT, lng_NEWTON_SOLVER_dealloc, lng_NEWTON_SOLVER_new,
    lng_NEWTON_SOLVER_methods, lng_NEWTON_SOLVER_members, lng_NEWTON_SOLVER_getset);

#if WITHSICONOS
  TYPEINIT (lng_SICONOS_SOLVER_TYPE, lng_SICONOS_SOLVER, "solfec.SICONOS_SOLVER",
    Py_TPFLAGS_DEFAULT, lng_SICONOS_SOLVER_dealloc, lng_SICONOS_SOLVER_new,
    lng_SICONOS_SOLVER_methods, lng_SICONOS_SOLVER_members, lng_SICONOS_SOLVER_getset);
#endif

  TYPEINIT (lng_TEST_SOLVER_TYPE, lng_TEST_SOLVER, "solfec.TEST_SOLVER",
    Py_TPFLAGS_DEFAULT, lng_TEST_SOLVER_dealloc, lng_TEST_SOLVER_new,
    lng_TEST_SOLVER_methods, lng_TEST_SOLVER_members, lng_TEST_SOLVER_getset);

  TYPEINIT (lng_CONSTRAINT_TYPE, lng_CONSTRAINT, "solfec.CONSTRAINT",
    Py_TPFLAGS_DEFAULT, lng_CONSTRAINT_dealloc, lng_CONSTRAINT_new,
    lng_CONSTRAINT_methods, lng_CONSTRAINT_members, lng_CONSTRAINT_getset);

  if (PyType_Ready (&lng_CONVEX_TYPE) < 0) return;
  if (PyType_Ready (&lng_MESH_TYPE) < 0) return;
  if (PyType_Ready (&lng_SPHERE_TYPE) < 0) return;
  if (PyType_Ready (&lng_ELLIP_TYPE) < 0) return;
  if (PyType_Ready (&lng_SOLFEC_TYPE) < 0) return;
  if (PyType_Ready (&lng_FIELD_TYPE) < 0) return;
  if (PyType_Ready (&lng_SURFACE_MATERIAL_TYPE) < 0) return;
  if (PyType_Ready (&lng_BULK_MATERIAL_TYPE) < 0) return;
  if (PyType_Ready (&lng_BODY_TYPE) < 0) return;
  if (PyType_Ready (&lng_TIME_SERIES_TYPE) < 0) return;
  if (PyType_Ready (&lng_GAUSS_SEIDEL_SOLVER_TYPE) < 0) return;
  if (PyType_Ready (&lng_PENALTY_SOLVER_TYPE) < 0) return;
  if (PyType_Ready (&lng_NEWTON_SOLVER_TYPE) < 0) return;
#if WITHSICONOS
  if (PyType_Ready (&lng_SICONOS_SOLVER_TYPE) < 0) return;
#endif
  if (PyType_Ready (&lng_TEST_SOLVER_TYPE) < 0) return;
  if (PyType_Ready (&lng_CONSTRAINT_TYPE) < 0) return;

  if (!(m =  Py_InitModule3 ("solfec", lng_methods, "Solfec module"))) return;

  Py_INCREF (&lng_CONVEX_TYPE);
  Py_INCREF (&lng_MESH_TYPE);
  Py_INCREF (&lng_SPHERE_TYPE);
  Py_INCREF (&lng_ELLIP_TYPE);
  Py_INCREF (&lng_SOLFEC_TYPE);
  Py_INCREF (&lng_FIELD_TYPE);
  Py_INCREF (&lng_SURFACE_MATERIAL_TYPE);
  Py_INCREF (&lng_BULK_MATERIAL_TYPE);
  Py_INCREF (&lng_BODY_TYPE);
  Py_INCREF (&lng_TIME_SERIES_TYPE);
  Py_INCREF (&lng_GAUSS_SEIDEL_SOLVER_TYPE);
  Py_INCREF (&lng_PENALTY_SOLVER_TYPE);
  Py_INCREF (&lng_NEWTON_SOLVER_TYPE);
#if WITHSICONOS
  Py_INCREF (&lng_SICONOS_SOLVER_TYPE);
#endif
  Py_INCREF (&lng_TEST_SOLVER_TYPE);
  Py_INCREF (&lng_CONSTRAINT_TYPE);

  PyModule_AddObject (m, "CONVEX", (PyObject*)&lng_CONVEX_TYPE);
  PyModule_AddObject (m, "MESH", (PyObject*)&lng_MESH_TYPE);
  PyModule_AddObject (m, "SPHERE", (PyObject*)&lng_SPHERE_TYPE);
  PyModule_AddObject (m, "ELLIP", (PyObject*)&lng_ELLIP_TYPE);
  PyModule_AddObject (m, "SOLFEC", (PyObject*)&lng_SOLFEC_TYPE);
  PyModule_AddObject (m, "FIELD", (PyObject*)&lng_FIELD_TYPE);
  PyModule_AddObject (m, "SURFACE_MATERIAL", (PyObject*)&lng_SURFACE_MATERIAL_TYPE);
  PyModule_AddObject (m, "BULK_MATERIAL", (PyObject*)&lng_BULK_MATERIAL_TYPE);
  PyModule_AddObject (m, "BODY", (PyObject*)&lng_BODY_TYPE);
  PyModule_AddObject (m, "TIME_SERIES", (PyObject*)&lng_TIME_SERIES_TYPE);
  PyModule_AddObject (m, "GAUSS_SEIDEL_SOLVER", (PyObject*)&lng_GAUSS_SEIDEL_SOLVER_TYPE);
  PyModule_AddObject (m, "PENALTY_SOLVER", (PyObject*)&lng_PENALTY_SOLVER_TYPE);
  PyModule_AddObject (m, "NEWTON_SOLVER", (PyObject*)&lng_NEWTON_SOLVER_TYPE);
#if WITHSICONOS
  PyModule_AddObject (m, "SICONOS_SOLVER", (PyObject*)&lng_SICONOS_SOLVER_TYPE);
#endif
  PyModule_AddObject (m, "TEST_SOLVER", (PyObject*)&lng_TEST_SOLVER_TYPE);
  PyModule_AddObject (m, "CONSTRAINT", (PyObject*)&lng_CONSTRAINT_TYPE);
  PyModule_AddObject (m, "__file__", ppath);
}

/* 
 * input file parsing
 */

/* interpret an input file (return 0 on success)*/
int lng (const char *path)
{
  FILE *file;
  char *line;
  int error;

  ASSERT (file = fopen (path, "r"), ERR_FILE_OPEN);

  Py_Initialize();

  initlng (path);

  PyRun_SimpleString("from solfec import CONVEX\n"
                     "from solfec import HULL\n"
                     "from solfec import MESH2CONVEX\n"
                     "from solfec import MESH\n"
                     "from solfec import HEX\n"
                     "from solfec import TETRAHEDRALIZE\n"
                     "from solfec import PIPE\n"
                     "from solfec import ROUGH_HEX\n"
                     "from solfec import SPHERE\n"
                     "from solfec import ELLIP\n"
                     "from solfec import SOLFEC\n"
                     "from solfec import FIELD\n"
                     "from solfec import SURFACE_MATERIAL\n"
                     "from solfec import BULK_MATERIAL\n"
                     "from solfec import BODY\n"
                     "from solfec import TIME_SERIES\n"
                     "from solfec import GAUSS_SEIDEL_SOLVER\n"
                     "from solfec import PENALTY_SOLVER\n"
                     "from solfec import NEWTON_SOLVER\n"
                     "from solfec import TEST_SOLVER\n"
                     "from solfec import FIX_POINT\n"
                     "from solfec import FIX_DIRECTION\n"
                     "from solfec import SET_DISPLACEMENT\n"
                     "from solfec import SET_VELOCITY\n"
                     "from solfec import SET_ACCELERATION\n"
                     "from solfec import PUT_RIGID_LINK\n"
                     "from solfec import PUT_SPRING\n"
                     "from solfec import GRAVITY\n"
                     "from solfec import FORCE\n"
                     "from solfec import TORQUE\n"
                     "from solfec import PRESSURE\n"
                     "from solfec import SIMPLIFIED_CRACK\n"
                     "from solfec import IMBALANCE_TOLERANCE\n"
                     "from solfec import RANK\n"
                     "from solfec import BARRIER\n"
                     "from solfec import NCPU\n"
                     "from solfec import HERE\n"
                     "from solfec import VIEWER\n"
                     "from solfec import BODY_CHARS\n"
                     "from solfec import INITIAL_VELOCITY\n"
                     "from solfec import MATERIAL\n"
                     "from solfec import DELETE\n"
                     "from solfec import SCALE\n"
                     "from solfec import TRANSLATE\n"
                     "from solfec import ROTATE\n"
                     "from solfec import SPLIT\n"
                     "from solfec import MESH_SPLIT\n"
                     "from solfec import COPY\n"
                     "from solfec import BEND\n"
                     "from solfec import BYLABEL\n"
                     "from solfec import MASS_CENTER\n"
                     "from solfec import CONTACT_EXCLUDE_BODIES\n"
                     "from solfec import CONTACT_EXCLUDE_SURFACES\n"
                     "from solfec import CONTACT_SPARSIFY\n"
                     "from solfec import RUN\n"
                     "from solfec import OUTPUT\n"
                     "from solfec import EXTENTS\n"
                     "from solfec import CALLBACK\n"
                     "from solfec import UNPHYSICAL_PENETRATION\n"
                     "from solfec import GEOMETRIC_EPSILON\n"
                     "from solfec import WARNINGS\n"
                     "from solfec import INITIALIZE_STATE\n"
                     "from solfec import LOCDYN_DUMP\n"
                     "from solfec import OVERLAPPING\n"
                     "from solfec import MBFCP_EXPORT\n"
                     "from solfec import NON_SOLFEC_ARGV\n"
                     "from solfec import MODAL_ANALYSIS\n"
                     "from solfec import BODY_MM_EXPORT\n"
                     "from solfec import DISPLAY_POINT\n"
                     "from solfec import FRACTURE_EXPORT_YAFFEMS\n"
                     "from solfec import DURATION\n"
                     "from solfec import RENDER\n"
                     "from solfec import FORWARD\n"
                     "from solfec import BACKWARD\n"
                     "from solfec import SEEK\n"
                     "from solfec import DISPLACEMENT\n"
                     "from solfec import VELOCITY\n"
                     "from solfec import STRESS\n"
                     "from solfec import ENERGY\n"
                     "from solfec import TIMING\n"
                     "from solfec import HISTORY\n"
                     "from solfec import __file__\n");

#if WITHSICONOS
  PyRun_SimpleString ("from solfec import SICONOS_SOLVER\n");
#endif

  // add a python function:
  // Split a string into a script name + args, and try to run this script
  // from the pwd passing it args.
  // NB: adding this function here is hacky, but it is much easier to define in Python
  PyRun_SimpleString("def runscript(cmd_str):\n"
                     "  import sys\n"
                     "  sys.argv = cmd_str.split()\n"
                     "  if not sys.argv[0].endswith('py'): sys.argv[0] += '.py'\n"
                     "  execfile ('%s' % sys.argv[0])\n");

  ERRMEM (line = MEM_CALLOC (128 + strlen (path)));
  sprintf (line, "execfile ('%s')", path);

  error = PyRun_SimpleString (line); /* we do not run a file directly because FILE destriptors differe
					between WIN32 and UNIX while Python is often provided in binary form */
  fclose (file);
  free (line);

  return error;
}

/* finalize interpreter */
void lngfinalize ()
{
  callback_pairs_destroy (); /* destroy all callback pairs */

  Py_Finalize(); /* delay this until the last moment => so that objects used in the viewer are not deallocated */
}

/* get positive id of a callback pointer pair;
 * return 0 if the pair was not found  */
int  lngcallback_id (void *data, void *call)
{
  CALLBACK_PAIR *pair;

  for (pair = callback_pairs; pair; pair = pair->next)
  {
    if (((void*)pair->data == data) &&
	((void*)pair->call == call)) return pair->id;
  }

  return 0;
}

/* set callback pointers for a given id; return 1 when
 * the id was found, or return 0 otherwise */
int  lngcallback_set (int id, void **data, void **call)
{
  CALLBACK_PAIR *pair;

  for (pair = callback_pairs; pair; pair = pair->next)
  {
    if (pair->id == id)
    {
      *data = pair->data;
      *call = pair->call;
      return 1;
    }
  }

  return 0;
}

/* handle PUT_SPRING spring Python callback */
double springcallback (void *call, double stroke, double velocity)
{
  double force = 0.0;
  PyObject *result;
  PyObject *args;

  args = Py_BuildValue ("(d, d)", stroke, velocity);

  result = PyObject_CallObject (call, args); /* call user function */

  Py_DECREF (args);

  if (result)
  {
    force  = PyFloat_AsDouble (result);

    Py_DECREF (result);
  }
  else /* Python call failed */
  {
    PyErr_Print (); /* print traceback */
#if MPI
    MPI_Abort (MPI_COMM_WORLD, 3000);
#endif
    exit (1);
  }

  return force;
}
