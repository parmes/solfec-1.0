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

#include <Python.h>
#include <structmember.h>
#include <limits.h>
#include <float.h>
#include "alg.h"
#include "sol.h"
#include "rnd.h"
#include "lng.h"
#include "fem.h"
#include "err.h"

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

    if (!(PyList_Size (obj) % div == 0 && PyList_Size (obj) >= len))
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

    if (PyTuple_Size (obj) != len)
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
	  sprintf (buf, "Vertex %d in face %d is outside of range [0, %d)",j , m, nver);
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
	  sprintf (buf, "Node %d in element %d is outside of range [0, %d)",j , m, nn);
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
	    sprintf (buf, "Node %d in face %d is outside of range [0, %d)", j, m, nn);
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

/* return node */
static PyObject* lng_MESH_node (lng_MESH *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("n");
  int n;

  PARSEKEYS ("i", &n);

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

  return Py_BuildValue ("(d, d, d)", self->msh->cur_nodes [n][0], self->msh->cur_nodes[n][1], self->msh->cur_nodes[n][2]);
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

/* MESH methods */
static PyMethodDef lng_MESH_methods [] =
{ 
  {"node", (PyCFunction)lng_MESH_node, METH_VARARGS|METH_KEYWORDS, "Return a node point of a mesh"},
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

/* constructor */
static PyObject* lng_SPHERE_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("center", "radius", "volid", "surfid", "sphere");
  lng_SPHERE *self, *sphere;
  double radius, c [3];
  int volid, surfid;
  PyObject *center;
  SPHERE *sph;

  self = (lng_SPHERE*)type->tp_alloc (type, 0);

  if (self)
  {
    sphere = NULL;

    PARSEKEYS ("Odii|O", &center, &radius, &volid, &surfid, &sphere);

    TYPETEST (is_tuple (center, kwl [0], 3) && is_sphere (sphere, kwl [4]));

    if (sphere) sph = sphere->sph;
    else sph = NULL;

    c [0] = PyFloat_AsDouble (PyTuple_GetItem (center, 0));
    c [1] = PyFloat_AsDouble (PyTuple_GetItem (center, 1));
    c [2] = PyFloat_AsDouble (PyTuple_GetItem (center, 2));

    self->sph = SPHERE_Create (sph, c, radius, surfid, volid);

    if (sph)
    {
      sphere->sph = NULL; /* empty */
    }
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

    PARSEKEYS ("OdO", &analysis, &step, &output);

    TYPETEST (is_string (analysis, kwl [0]) && is_string (output, kwl [2]));

    outpath = PyString_AsString (output);

    IFIS (analysis, "QUASI_STATIC")
    {
      self->sol = SOLFEC_Create (0, step, outpath);
#if OPENGL
      RND_Domain (self->sol->dom); /* just in case a viewer is enabled (last created SOLFEC object) */
#endif
    }
    ELIF (analysis, "DYNAMIC")
    {
      self->sol = SOLFEC_Create (1, step, outpath);
#if OPENGL
      RND_Domain (self->sol->dom); /* pass last created SOLFEC object to the rendering */
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

  if (!is_number (value, "step")) return -1;
  
  step = PyFloat_AsDouble (value);

  if (!is_positive (step, "step")) return -1;

  self->sol->dom->step = step;

  return 0;
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
  {NULL, 0, 0, NULL, NULL} 
};

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
    "dashpot");

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

    PARSEKEYS ("O|iiOOddddd", &solfec, &mat.surf1, &mat.surf2, &model, &label,
      &mat.friction, &mat.cohesion, &mat.restitution, &mat.spring, &mat.dashpot);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_string (model, kwl[3]) && is_string (label, kwl[4]) &&
	      is_non_negative (mat.friction, kwl[5]) && is_non_negative (mat.cohesion, kwl[6]) &&
	      is_in_range (mat.restitution, kwl[7], 0, 1) && is_non_negative (mat.spring, kwl[8]) &&
	      is_non_negative (mat.dashpot, kwl[9]));

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
  if (!is_number (value, "friction")) return -1;
  self->mat->friction = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_SURFACE_MATERIAL_get_cohesion (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->cohesion);
}

static int lng_SURFACE_MATERIAL_set_cohesion (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number (value, "cohesion")) return -1;
  self->mat->cohesion = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_SURFACE_MATERIAL_get_restitution (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->restitution);
}

static int lng_SURFACE_MATERIAL_set_restitution (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number (value, "restitution")) return -1;
  self->mat->restitution = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_SURFACE_MATERIAL_get_spring (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->spring);
}

static int lng_SURFACE_MATERIAL_set_spring (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number (value, "spring")) return -1;
  self->mat->spring = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_SURFACE_MATERIAL_get_dashpot (lng_SURFACE_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->dashpot);
}

static int lng_SURFACE_MATERIAL_set_dashpot (lng_SURFACE_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number (value, "dashpot")) return -1;
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
    "density");

  BULK_MATERIAL mat = 
  {
   NULL, /* label */
   KIRCHHOFF, /* model */
   1E6,  /* young */
   0.25, /* poisson */
   1E3   /* density */
  };

  PyObject *model, *label;
  lng_BULK_MATERIAL *self;
  lng_SOLFEC *solfec;
  SOLFEC *sol;

  self = (lng_BULK_MATERIAL*)type->tp_alloc (type, 0);

  if (self)
  {
    model = NULL;
    label = NULL;

    PARSEKEYS ("O|OOddd", &solfec, &model, &label, &mat.young, &mat.poisson, &mat.density);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_string (model, kwl[1]) && is_string (label, kwl[2]) &&
	      is_non_negative (mat.young, kwl[3]) && is_non_negative (mat.poisson, kwl[4]) &&
	      is_non_negative (mat.density, kwl[5]));

    sol = solfec->sol;

    if (model)
    {
      IFIS (model, "KIRCHHOFF") mat.model = KIRCHHOFF;
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
  if (!is_number (value, "young")) return -1;
  self->mat->young = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_BULK_MATERIAL_get_poisson (lng_BULK_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->poisson);
}

static int lng_BULK_MATERIAL_set_poisson (lng_BULK_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number (value, "poisson")) return -1;
  self->mat->poisson = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_BULK_MATERIAL_get_density (lng_BULK_MATERIAL *self, void *closure)
{
  return PyFloat_FromDouble (self->mat->density);
}

static int lng_BULK_MATERIAL_set_density (lng_BULK_MATERIAL *self, PyObject *value, void *closure)
{
  if (!is_number (value, "density")) return -1;
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

  short dodestroy;

  BODY *bod;
};

#if MPI
/* try assigning a parent body pointer to the id */
static int ID_TO_BODY (lng_BODY *body)
{
  if (body->id)
  {
    BODY *bod = MAP_Find (body->dom->idb, (void*) (long) body->id, NULL);

    if (bod && !(bod->flags & BODY_CHILD))
    {
      body->bod = bod;
      return 1;
    }
  }

  return 0;
}

/* get rank */
static int RANK ()
{
  int rank;

  MPI_Comm_rank (MPI_COMM_WORLD, &rank);

  return rank;
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
      sprintf (buf, "'%s' must be a non-empty CONVEX/MESH/SPHERE object or a list of those", var);
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
    if (!is_basic_shape (obj))
    {
      if (PyList_Check (obj))
      {
	int i, n = PyList_Size (obj);

	for (i = 0; i < n; i ++)
	  if (!is_convex_test (PyList_GetItem (obj, i))) break;

	if (i == n) return 1;
      }

      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a non-empty CONVEX object or a list of those", var);
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
	sprintf (buf, "'%s' must be a non-empty CONVEX/MESH/SPHERE object, a list of those or a (x, y, z) tuple", var);
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
  return PyObject_IsInstance (obj, (PyObject*)&lng_BODY_TYPE);
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
	sprintf (buf, "'%s' must be a non-empty CONVEX/MESH/SPHERE object, a list of those or a BODY object", var);
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
	if (empty) out = SHAPE_Glue (SHAPE_Create (shape_kind (it), get_shape (it, empty)), out); /* glue simple shapes (destructive for simple shape lists) */
	else out = SHAPE_Glue_Simple (SHAPE_Create (shape_kind (it), get_shape (it, empty)), out); /* do not glue simple shape (non destructive) */
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
  }

  return 1;
}

/* body object constructor */
static PyObject* lng_BODY_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "kind", "shape", "material", "label", "formulation", "mesh");
  PyObject *kind, *shape, *material, *label, *formulation;
  lng_SOLFEC *solfec;
  lng_BODY *self;
  lng_MESH *mesh;
  short form;
  char *lab;

  self = (lng_BODY*)type->tp_alloc (type, 0);

  if (self)
  {
    self->dodestroy = 0; /* never destroy by default as the body will be imediately owned by a domain */

    label = NULL;
    formulation = NULL;
    form = FEM_O1;
    mesh = NULL;

    PARSEKEYS ("OOOO|OOO", &solfec, &kind, &shape, &material, &label, &formulation, &mesh);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_string (kind, kwl[1]) && is_shape (shape, kwl[2]) &&
	      is_bulk_material (solfec->sol, material, kwl[3]) && is_string (label, kwl[4]) &&
	      is_string (formulation, kwl[5]) && is_mesh ((PyObject*)mesh, kwl[6]));

#if MPI
    if (RANK() > 0) /* bodies can only be created on process zero */
    {
      self->dom = solfec->sol->dom;
      self->id = 0; /* identifiers start from 1 => always invalid (ID_TO_BODY returns 0) */
      return (PyObject*)self; /* return a zero initialized objecy (tp_alloc) */
    }
#endif

    if (label) lab = PyString_AsString (label);
    else lab = NULL;

    IFIS (kind, "RIGID")
    {
      self->bod = BODY_Create (RIG, create_shape (shape, 1), get_bulk_material (solfec->sol, material), lab, form, NULL);
    }
    ELIF (kind, "PSEUDO_RIGID")
    {
      self->bod = BODY_Create (PRB, create_shape (shape, 1), get_bulk_material (solfec->sol, material), lab, form, NULL);
    }
    ELIF (kind, "ROUGH_FINITE_ELEMENT")
    {
      TYPETEST (is_shape_convex (shape, kwl[2]) && is_mesh ((PyObject*)mesh, kwl[6]));

      self->bod = BODY_Create (RFE, create_shape (shape, 1), get_bulk_material (solfec->sol, material), lab, form, mesh->msh);
      mesh->msh = NULL; /* empty */
    }
    ELIF (kind, "FINITE_ELEMENT")
    {
      TYPETEST (is_mesh (shape, kwl[2]));

      if (formulation)
      {
	IFIS (formulation, "FEM_O1")
	{
	  form = FEM_O1;
	}
	ELIF (formulation, "FEM_O2")
	{
	  form = FEM_O2;
	}
	ELSE
	{
	  PyErr_SetString (PyExc_ValueError, "Invalid FEM formulation");
	  return NULL;
	}
      }

      self->bod = BODY_Create (FEM, create_shape (shape, 1), get_bulk_material (solfec->sol, material), lab, form, NULL);
    }
    ELIF (kind, "OBSTACLE")
    {
      self->bod = BODY_Create (OBS, create_shape (shape, 1), get_bulk_material (solfec->sol, material), lab, form, NULL);
    }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Unknown BODY kind");
      return NULL;
    }

    DOM_Insert_Body (solfec->sol->dom, self->bod); /* insert body into the domain */

#if MPI
    self->id = self->bod->id;
    self->dom = solfec->sol->dom;
#endif
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
  if (self->dodestroy) BODY_Destroy (self->bod); /* used only when body was removed from the domain */

  self->ob_type->tp_free ((PyObject*)self);
}

/* getsets */

static PyObject* lng_BODY_get_kind (lng_BODY *self, void *closure)
{
#if MPI
  if (ID_TO_BODY (self))
  {
#endif

  return PyString_FromString (BODY_Kind (self->bod));

#if MPI
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
#if MPI
  if (ID_TO_BODY (self))
  {
#endif

  if (self->bod->label)
    return PyString_FromString (self->bod->label);
  else return PyString_FromString (""); /* an empty label */

#if MPI
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
  if (ID_TO_BODY (self))
  {
#endif

  if (self->bod->kind == OBS) Py_RETURN_NONE;

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
  if (ID_TO_BODY (self))
  {
#endif

  if (self->bod->kind == OBS) Py_RETURN_NONE;

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

static PyObject* lng_BODY_get_selfcontact (lng_BODY *self, void *closure)
{
#if MPI
  if (ID_TO_BODY (self))
  {
#endif

  if (self->bod->flags & BODY_DETECT_SELF_CONTACT)
    return PyString_FromString ("ON");
  else return PyString_FromString ("OFF");

#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

static int lng_BODY_set_selfcontact (lng_BODY *self, PyObject *value, void *closure)
{
  if (!is_string (value, "selfcontact")) return -1;

#if MPI
  if (ID_TO_BODY (self))
  {
#endif

  IFIS (value, "ON") self->bod->flags |= BODY_DETECT_SELF_CONTACT;
  ELIF (value, "OFF") self->bod->flags &= ~BODY_DETECT_SELF_CONTACT;
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid self-contact flag");
    return -1;
  }

#if MPI
  }
#endif

  return 0;
}

static PyObject* lng_BODY_get_scheme (lng_BODY *self, void *closure)
{
#if MPI
  if (ID_TO_BODY (self))
  {
#endif

  switch (self->bod->scheme)
  {
  case SCH_DEFAULT: return PyString_FromString ("DEFAULT");
  case SCH_RIG_NEG: return PyString_FromString ("RIG_NEG");
  case SCH_RIG_POS: return PyString_FromString ("RIG_POS");
  case SCH_RIG_IMP: return PyString_FromString ("RIG_IMP");
  }

#if MPI
  }
  else Py_RETURN_NONE;
#endif

  return NULL;
}

static int lng_BODY_set_scheme (lng_BODY *self, PyObject *value, void *closure)
{
  if (!is_string (value, "scheme")) return -1;

#if MPI
  if (ID_TO_BODY (self))
  {
#endif

  IFIS (value, "DEFAULT") self->bod->scheme = SCH_DEFAULT;
  ELIF (value, "RIG_POS")
  {
    if (self->bod->kind != RIG)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid integration scheme for a rigid body");
      return -1;
    }

    self->bod->scheme = SCH_RIG_POS;
  }
  ELIF (value, "RIG_NEG")
  {
    if (self->bod->kind != RIG)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid integration scheme for a rigid body");
      return -1;
    }

    self->bod->scheme = SCH_RIG_NEG;
  }
  ELIF (value, "RIG_IMP")
  {
    if (self->bod->kind != RIG)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid integration scheme for a rigid body");
      return -1;
    }

    self->bod->scheme = SCH_RIG_IMP;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid integration scheme");
    return -1;
  }

#if MPI
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
  {"kind", (getter)lng_BODY_get_kind, (setter)lng_BODY_set_kind, "kind", NULL},
  {"label", (getter)lng_BODY_get_label, (setter)lng_BODY_set_label, "label", NULL},
  {"conf", (getter)lng_BODY_get_conf, (setter)lng_BODY_set_conf, "configuration", NULL},
  {"velo", (getter)lng_BODY_get_velo, (setter)lng_BODY_set_velo, "velocity", NULL},
  {"selfcontact", (getter)lng_BODY_get_selfcontact, (setter)lng_BODY_set_selfcontact, "selfcontact", NULL},
  {"scheme", (getter)lng_BODY_get_scheme, (setter)lng_BODY_set_scheme, "scheme", NULL},
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
      exit (1);
    }

    Py_DECREF (result);
  }
  else /* error during the Python callback run */
  {
err:
    PyErr_Print (); /* print traceback */
    exit (1);
  }
}

/* constructor */
static PyObject* lng_GAUSS_SEIDEL_SOLVER_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("epsilon", "maxiter", "failure", "diagepsilon", "diagmaxiter", "diagsolver", "data", "callback");
  PyObject *failure, *diagsolver;
  lng_GAUSS_SEIDEL_SOLVER *self;
  double epsilon, diagepsilon;
  int maxiter, diagmaxiter;
  GSFAIL gsfail;
  GSDIAS gsdias;

  self = (lng_GAUSS_SEIDEL_SOLVER*)type->tp_alloc (type, 0);

  if (self)
  {
    failure = NULL;
    diagepsilon = DBL_MAX;
    diagmaxiter = INT_MAX;
    diagsolver = NULL;
    gsfail = GS_FAILURE_CONTINUE;
    gsdias = GS_SEMISMOOTH_NEWTON;
    self->data = NULL;
    self->callback = NULL;

    PARSEKEYS ("di|OdiOOO", &epsilon, &maxiter, &failure,
      &diagepsilon, &diagmaxiter, &diagsolver, &self->data, &self->callback);

    TYPETEST (is_positive (epsilon, kwl[0]) && is_positive (maxiter, kwl[1]) &&
      is_string (failure, kwl[2]) && is_positive (diagepsilon, kwl[3]) &&
      is_positive (diagmaxiter, kwl[4]) && is_string (diagsolver, kwl[5]) &&
      is_callable (self->callback, kwl[7]));

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
	gsdias = GS_SEMISMOOTH_NEWTON;
      }
      ELIF (diagsolver, "PROJECTED_GRADIENT")
      {
	gsdias = GS_PROJECTED_GRADIENT;
      }
      ELIF (diagsolver, "DE_SAXE_AND_FENG")
      {
	gsdias = GS_DE_SAXE_AND_FENG;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid diagonal solver");
	return NULL;
      }
    }

    if (diagepsilon == DBL_MAX)
      diagepsilon = 0.01 * epsilon;

    if (diagmaxiter == INT_MAX)
      diagmaxiter = MAX (100, maxiter / 100);

    self->gs = GAUSS_SEIDEL_Create (epsilon, maxiter, gsfail, diagepsilon,
      diagmaxiter, gsdias, self, (GAUSS_SEIDEL_Callback)lng_GAUSS_SEIDEL_callback);
  }

  return (PyObject*)self;
}

/* destructor */
static void lng_GAUSS_SEIDEL_SOLVER_dealloc (lng_GAUSS_SEIDEL_SOLVER *self)
{
  Py_XDECREF (self->callback);
  Py_XDECREF (self->data);

  GAUSS_SEIDEL_Destroy (self->gs);

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

  if (self->gs->history == GS_OFF)
  {
    PyErr_SetString (PyExc_ValueError, "Hisoty recording is switched off");
    return NULL;
  }

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

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_epsilon (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->gs->epsilon);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_epsilon (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number (value, "epsilon")) return -1;
  self->gs->epsilon = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_maxiter (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->gs->maxiter);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_maxiter (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number (value, "maxiter")) return -1;
  self->gs->maxiter = PyInt_AsLong (value);
  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_diagepsilon (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->gs->diagepsilon);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_diagepsilon (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number (value, "diagepsilon")) return -1;
  self->gs->diagepsilon = PyFloat_AsDouble (value);
  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_diagmaxiter (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyFloat_FromDouble (self->gs->diagmaxiter);
}

static int lng_GAUSS_SEIDEL_SOLVER_set_diagmaxiter (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_number (value, "diagmaxiter")) return -1;
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
    self->gs->diagsolver = GS_SEMISMOOTH_NEWTON;
  }
  ELIF (value, "PROJECTED_GRADIENT")
  {
    self->gs->diagsolver = GS_PROJECTED_GRADIENT;
  }
  ELIF (value, "DE_SAXE_AND_FENG")
  {
    self->gs->diagsolver = GS_DE_SAXE_AND_FENG;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid diagonal solver");
    return -1;
  }

  return 0;
}

static PyObject* lng_GAUSS_SEIDEL_SOLVER_get_history (lng_GAUSS_SEIDEL_SOLVER *self, void *closure)
{
  return PyString_FromString (GAUSS_SEIDEL_History (self->gs));
}

static int lng_GAUSS_SEIDEL_SOLVER_set_history (lng_GAUSS_SEIDEL_SOLVER *self, PyObject *value, void *closure)
{
  if (!is_string (value, "history")) return -1;

  IFIS (value, "ON")
  {
    self->gs->history = GS_ON;
  }
  ELIF (value, "OFF")
  {
    self->gs->history = GS_OFF;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid history switch (ON/OFF accepted)");
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

  IFIS (value, "MID_LOOP")
  {
    self->gs->variant = GS_MID_LOOP;
  }
  ELIF (value, "MID_THREAD")
  {
    self->gs->variant = GS_MID_THREAD;
  }
  ELIF (value, "MID_TO_ALL")
  {
    self->gs->variant = GS_MID_TO_ALL;
  }
  ELIF (value, "MID_TO_ONE")
  {
    self->gs->variant = GS_MID_TO_ONE;
  }
  ELIF (value, "NOB_MID_LOOP")
  {
    self->gs->variant = GS_NOB_MID_LOOP;
  }
  ELIF (value, "NOB_MID_THREAD")
  {
    self->gs->variant = GS_NOB_MID_THREAD;
  }
  ELIF (value, "NOB_MID_TO_ALL")
  {
    self->gs->variant = GS_NOB_MID_TO_ALL;
  }
  ELIF (value, "NOB_MID_TO_ONE")
  {
    self->gs->variant = GS_NOB_MID_TO_ONE;
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
  {"epsilon", (getter)lng_GAUSS_SEIDEL_SOLVER_get_epsilon, (setter)lng_GAUSS_SEIDEL_SOLVER_set_epsilon, "relative accuracy", NULL},
  {"maxiter", (getter)lng_GAUSS_SEIDEL_SOLVER_get_maxiter, (setter)lng_GAUSS_SEIDEL_SOLVER_set_maxiter, "iterations bound", NULL},
  {"diagepsilon", (getter)lng_GAUSS_SEIDEL_SOLVER_get_diagepsilon, (setter)lng_GAUSS_SEIDEL_SOLVER_set_diagepsilon, "diagonal solver relative accuracy", NULL},
  {"diagmaxiter", (getter)lng_GAUSS_SEIDEL_SOLVER_get_diagmaxiter, (setter)lng_GAUSS_SEIDEL_SOLVER_set_diagmaxiter, "diagonal solver iterations bound", NULL},
  {"diagsolver", (getter)lng_GAUSS_SEIDEL_SOLVER_get_diagsolver, (setter)lng_GAUSS_SEIDEL_SOLVER_set_diagsolver, "diagonal solver kind", NULL},
  {"history", (getter)lng_GAUSS_SEIDEL_SOLVER_get_history, (setter)lng_GAUSS_SEIDEL_SOLVER_set_history, "solution history recording", NULL},
  {"variant", (getter)lng_GAUSS_SEIDEL_SOLVER_get_variant, (setter)lng_GAUSS_SEIDEL_SOLVER_set_variant, "parallel method variant", NULL},
  {"reverse", (getter)lng_GAUSS_SEIDEL_SOLVER_get_reverse, (setter)lng_GAUSS_SEIDEL_SOLVER_set_reverse, "iteration reversion flag", NULL},
  {NULL, 0, 0, NULL, NULL}
};

/*
 * EXPLICIT_SOLVER => object
 */

typedef struct lng_EXPLICIT_SOLVER lng_EXPLICIT_SOLVER;

static PyTypeObject lng_EXPLICIT_SOLVER_TYPE;

struct lng_EXPLICIT_SOLVER
{
  PyObject_HEAD
};

/* constructor */
static PyObject* lng_EXPLICIT_SOLVER_new (PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  lng_EXPLICIT_SOLVER *self;

  self = (lng_EXPLICIT_SOLVER*)type->tp_alloc (type, 0);

  return (PyObject*)self;
}

/* destructor */
static void lng_EXPLICIT_SOLVER_dealloc (lng_EXPLICIT_SOLVER *self)
{
  self->ob_type->tp_free ((PyObject*)self);
}

/* EXPLICIT_SOLVER methods */
static PyMethodDef lng_EXPLICIT_SOLVER_methods [] =
{ {NULL, NULL, 0, NULL} };

/* EXPLICIT_SOLVER members */
static PyMemberDef lng_EXPLICIT_SOLVER_members [] =
{ {NULL, 0, 0, 0, NULL} };

/* EXPLICIT_SOLVER getset */
static PyGetSetDef lng_EXPLICIT_SOLVER_getset [] =
{ {NULL, 0, 0, NULL, NULL} };

/*
 * CONSTRAINT => object
 */

typedef struct lng_CONSTRAINT lng_CONSTRAINT; /* constraint type */

static PyTypeObject lng_CONSTRAINT_TYPE;

struct lng_CONSTRAINT
{
  PyObject_HEAD

#if MPI
  unsigned int id;
#endif

  CON *con;
};

#if MPI
/* try assigning a constraint pointer to the id */
static int ID_TO_CONSTRAINT (lng_SOLFEC *solfec, lng_CONSTRAINT *constraint)
{
  if (constraint->id)
  {
    CON *con = MAP_Find (solfec->sol->dom->idc, (void*) (long) constraint->id, NULL);

    if (con)
    {
      constraint->con = con;
      return 1;
    }
  }

  return 0;
}
#endif

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

#if MPI
  self->id = con->id;
#endif

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
  return PyString_FromString (CON_Kind (self->con));
}

static int lng_CONSTRAINT_set_kind (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_R (lng_CONSTRAINT *self, void *closure)
{
  double *R = self->con->dia->R;
  return Py_BuildValue ("(d, d, d)", R[0], R[1], R[2]);
}

static int lng_CONSTRAINT_set_R (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_base (lng_CONSTRAINT *self, void *closure)
{
  double *base = self->con->base;
  return Py_BuildValue ("(d, d, d, d, d, d, d, d, d)", base [0], base [1],
     base [2], base [3], base [4], base [5], base [6], base [7], base [8]);
}

static int lng_CONSTRAINT_set_base (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_point (lng_CONSTRAINT *self, void *closure)
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

static int lng_CONSTRAINT_set_point (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_adjbod (lng_CONSTRAINT *self, void *closure)
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

static int lng_CONSTRAINT_set_adjbod (lng_CONSTRAINT *self, PyObject *value, void *closure)
{
  PyErr_SetString (PyExc_ValueError, "Writing to a read-only member");
  return -1;
}

static PyObject* lng_CONSTRAINT_get_matlab (lng_CONSTRAINT *self, void *closure)
{
  if (self->con->kind == CONTACT)
    return PyString_FromString (self->con->mat.label);

  Py_RETURN_NONE;
}

static int lng_CONSTRAINT_set_matlab (lng_CONSTRAINT *self, PyObject *value, void *closure)
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
  {"base", (getter)lng_CONSTRAINT_get_base, (setter)lng_CONSTRAINT_set_base, "constraint local base", NULL},
  {"point", (getter)lng_CONSTRAINT_get_point, (setter)lng_CONSTRAINT_set_point, "constraint spatial point", NULL},
  {"adjbod", (getter)lng_CONSTRAINT_get_adjbod, (setter)lng_CONSTRAINT_set_adjbod, "constraint adjacent bodies", NULL},
  {"matlab", (getter)lng_CONSTRAINT_get_matlab, (setter)lng_CONSTRAINT_set_matlab, "contact constraint material", NULL},
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
      return NULL;
    }
    ENDTRY ()

    if (cvx)
    {
      convex->cvx = NULL; /* empty */
    }

    free (pnt);
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

/* create fixed point constraint */
static PyObject* lng_FIX_POINT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "body", "point");
  lng_CONSTRAINT *out;
  lng_SOLFEC *solfec;
  lng_BODY *body;
  PyObject *point;
  double p [3];

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOO", &solfec, &body, &point);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body, kwl[1]) && is_tuple (point, kwl[2], 3));

#if MPI
    if (ID_TO_BODY (body))
    {
#endif

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    if (!(out->con = DOM_Fix_Point (solfec->sol->dom, body->bod, p)))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create fixed direction constraint */
static PyObject* lng_FIX_DIRECTION (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "body", "point", "direction");
  lng_CONSTRAINT *out;
  lng_SOLFEC *solfec;
  lng_BODY *body;
  PyObject *point, *direction;
  double p [3], d [3];

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOOO", &solfec, &body, &point, &direction);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body, kwl[1]) &&
      is_tuple (point, kwl[2], 3) && is_tuple (direction, kwl[3], 3));

#if MPI
    if (ID_TO_BODY (body))
    {
#endif

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
    d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
    d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

    if (!(out->con = DOM_Fix_Direction (solfec->sol->dom, body->bod, p, d)))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create prescribed displacement constraint */
static PyObject* lng_SET_DISPLACEMENT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "body", "point", "direction", "tms");
  lng_CONSTRAINT *out;
  lng_SOLFEC *solfec;
  lng_BODY *body;
  PyObject *point, *direction;
  lng_TIME_SERIES *tms;
  double p [3], d [3];

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOOO", &solfec, &body, &point, &direction, &tms);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body, kwl[1]) &&
      is_tuple (point, kwl[2], 3) && is_tuple (direction, kwl[3], 3) &&
      is_time_series (tms, kwl[4]));

#if MPI
    if (ID_TO_BODY (body))
    {
#endif

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
    d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
    d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

    if (!(out->con = DOM_Set_Velocity (solfec->sol->dom, body->bod, p, d, TMS_Derivative (tms->ts))))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create prescribed velocity constraint */
static PyObject* lng_SET_VELOCITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "body", "point", "direction", "tms");
  lng_CONSTRAINT *out;
  lng_SOLFEC *solfec;
  lng_BODY *body;
  PyObject *point, *direction;
  lng_TIME_SERIES *tms;
  double p [3], d [3];

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOOO", &solfec, &body, &point, &direction, &tms);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body, kwl[1]) &&
      is_tuple (point, kwl[2], 3) && is_tuple (direction, kwl[3], 3) &&
      is_time_series (tms, kwl[4]));

#if MPI
    if (ID_TO_BODY (body))
    {
#endif

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
    d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
    d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

    if (!(out->con = DOM_Set_Velocity (solfec->sol->dom, body->bod, p, d, TMS_Copy (tms->ts))))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create prescribed acceleration constraint */
static PyObject* lng_SET_ACCELERATION (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "body", "point", "direction", "tms");
  lng_CONSTRAINT *out;
  lng_SOLFEC *solfec;
  lng_BODY *body;
  PyObject *point, *direction;
  lng_TIME_SERIES *tms;
  double p [3], d [3];

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOOO", &solfec, &body, &point, &direction, &tms);

    TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body, kwl[1]) &&
      is_tuple (point, kwl[2], 3) && is_tuple (direction, kwl[3], 3) &&
      is_time_series (tms, kwl[4]));

#if MPI
    if (ID_TO_BODY (body))
    {
#endif

    p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
    d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
    d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

    if (!(out->con = DOM_Set_Velocity (solfec->sol->dom, body->bod, p, d, TMS_Integral (tms->ts))))
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* create rigid link constraint */
static PyObject* lng_PUT_RIGID_LINK (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "body1", "body2", "point1", "point2");
  lng_CONSTRAINT *out;
  lng_SOLFEC *solfec;
  lng_BODY *body1, *body2;
  PyObject *point1, *point2;
  double p1 [3], p2 [3];

  out = (lng_CONSTRAINT*)lng_CONSTRAINT_TYPE.tp_alloc (&lng_CONSTRAINT_TYPE, 0);

  if (out)
  {
    PARSEKEYS ("OOOOO", &solfec, &body1, &body2, &point1, &point2);

    if ((PyObject*)body1 == Py_None)
    {
      TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body2, kwl[2]) &&
	        is_tuple (point1, kwl[3], 3) && is_tuple (point2, kwl[4], 3));
    }
    else if ((PyObject*)body2 == Py_None)
    {
      TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body1, kwl[1]) &&
	        is_tuple (point1, kwl[3], 3) && is_tuple (point2, kwl[4], 3));
    }
    else
    {
      TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body1, kwl[1]) && is_body (body2, kwl[2]) &&
		is_tuple (point1, kwl[3], 3) && is_tuple (point2, kwl[4], 3));
    }

#if MPI
    if (((PyObject*)body1 == Py_None && ID_TO_BODY (body2)) ||
	((PyObject*)body2 == Py_None && ID_TO_BODY (body1)) ||
	(ID_TO_BODY (body1) && ID_TO_BODY (body2)))
    {
#endif

    p1 [0] = PyFloat_AsDouble (PyTuple_GetItem (point1, 0));
    p1 [1] = PyFloat_AsDouble (PyTuple_GetItem (point1, 1));
    p1 [2] = PyFloat_AsDouble (PyTuple_GetItem (point1, 2));

    p2 [0] = PyFloat_AsDouble (PyTuple_GetItem (point2, 0));
    p2 [1] = PyFloat_AsDouble (PyTuple_GetItem (point2, 1));
    p2 [2] = PyFloat_AsDouble (PyTuple_GetItem (point2, 2));

    if ((PyObject*)body1 == Py_None) out->con = DOM_Put_Rigid_Link (solfec->sol->dom, NULL, body2->bod, p1, p2);
    else if ((PyObject*)body2 == Py_None) out->con = DOM_Put_Rigid_Link (solfec->sol->dom, body1->bod, NULL, p1, p2);
    else out->con = DOM_Put_Rigid_Link (solfec->sol->dom, body1->bod, body2->bod, p1, p2);

    if (!out->con)
    {
      PyErr_SetString (PyExc_ValueError, "Point outside of domain");
      return NULL;
    }

#if MPI
    }
#endif
  }

  return (PyObject*)out;
}

/* set gravity */
static PyObject* lng_GRAVITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "direction", "value");
  PyObject *direction, *value;
  lng_SOLFEC *solfec;
  double d [3];
  TMS *ts;

  PARSEKEYS ("OOO", &solfec, &direction, &value);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_tuple (direction, kwl[1], 3) &&
	    is_number_or_time_series (value, kwl[2]));

  d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
  d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
  d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

  if (PyNumber_Check (value))
    ts = TMS_Constant (PyFloat_AsDouble(value));
  else ts = TMS_Copy (((lng_TIME_SERIES*)value)->ts);

  COPY (d, solfec->sol->dom->gravdir);
  if  (solfec->sol->dom->gravval)
    TMS_Destroy (solfec->sol->dom->gravval);
  solfec->sol->dom->gravval = ts;

  Py_RETURN_NONE;
}

/* callback for forces defined by a callable object */
static void lng_FORCE_callback (PyObject *data, PyObject *call, int nq, int nu, double *q, double *u, double t, double h, double *f)
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

      PyTuple_SetItem (args, n, PyTuple_GetItem (data, n));
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

    if (!PyTuple_Size (result) != (nu == 6 ? 9 : nu))
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

#if MPI
  if (ID_TO_BODY (body))
  {
#endif 

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  d [0] = PyFloat_AsDouble (PyTuple_GetItem (direction, 0));
  d [1] = PyFloat_AsDouble (PyTuple_GetItem (direction, 1));
  d [2] = PyFloat_AsDouble (PyTuple_GetItem (direction, 2));

  if (PyNumber_Check (value))
    ts = TMS_Constant (PyFloat_AsDouble(value));
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

  BODY_Apply_Force (body->bod, k, p, d, ts, call, func);

#if MPI
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

#if MPI
  if (ID_TO_BODY (body))
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

  if (PyNumber_Check (value))
    ts = TMS_Constant (PyFloat_AsDouble(value));
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

  BODY_Apply_Force (body->bod, k | TORQUE, NULL, d, ts, NULL, NULL);

#if MPI
  }
#endif

  Py_RETURN_NONE;
}

/* set imbalance tolerances */
static PyObject* lng_IMBALANCE_TOLERANCES (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "timint", "condet", "locdyn");
  double timint, condet, locdyn;
  lng_SOLFEC *solfec;

  PARSEKEYS ("Oddd", &solfec, &timint, &condet, &locdyn);

  TYPETEST (is_solfec (solfec, kwl[0]));

#if MPI
  solfec->sol->dom->imbalance_tolerance = timint;
  solfec->sol->aabb->imbalance_tolerance = condet;
  solfec->sol->dom->ldy->imbalance_tolerance = locdyn;
#endif

  Py_RETURN_TRUE;
}

/* change LOCDYN balancing approach */
static PyObject* lng_LOCDYN_BALANCING (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "mode");
  lng_SOLFEC *solfec;
  PyObject *mode;

  PARSEKEYS ("OO", &solfec, &mode);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_string (mode, kwl[1]));

  IFIS (mode, "OFF")
  {
#if MPI
    LOCDYN_Balancing (solfec->sol->dom->ldy, LDB_OFF);
#endif
  }
  ELIF (mode, "GEOM")
  {
#if MPI
    LOCDYN_Balancing (solfec->sol->dom->ldy, LDB_GEOM);
#endif
  }
  ELIF (mode, "GRAPH")
  {
#if MPI
    LOCDYN_Balancing (solfec->sol->dom->ldy, LDB_GRAPH);
#endif
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid balancing mode");
    return NULL;
  }

  Py_RETURN_TRUE;
}

/* return number of CPUs */
static PyObject* lng_NCPU (PyObject *self, PyObject *args, PyObject *kwds)
{
#if MPI
  int ncpu;

  MPI_Comm_size (MPI_COMM_WORLD, &ncpu);
  return PyInt_FromLong (ncpu);
#else
  return PyInt_FromLong (1);
#endif
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

  PARSEKEYS ("O", &object);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_body_or_constraint (object, kwl[1]));

  if (PyObject_IsInstance (object, (PyObject*)&lng_BODY_TYPE))
  {
    body = (lng_BODY*)object;

    if (ID_TO_BODY (body)) Py_RETURN_TRUE;
    else Py_RETURN_FALSE;
  }
  else
  {
    constraint = (lng_CONSTRAINT*)object;

    if (ID_TO_CONSTRAINT (solfec, constraint)) Py_RETURN_TRUE;
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
  if (RND_On ()) Py_RETURN_TRUE; /* return true and maintain reference count of Py_True */
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

#if MPI
  if (ID_TO_BODY (body))
  {
#endif

  c [0] = PyFloat_AsDouble (PyTuple_GetItem (center, 0));
  c [1] = PyFloat_AsDouble (PyTuple_GetItem (center, 1));
  c [2] = PyFloat_AsDouble (PyTuple_GetItem (center, 2));

  for (i = 0; i < 9; i ++) t [i] = PyFloat_AsDouble (PyTuple_GetItem (tensor, i));

  BODY_Overwrite_Chars (body->bod, mass, volume, c, t);

#if MPI
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

#if MPI
  if (ID_TO_BODY (body))
  {
#endif

  l [0] = PyFloat_AsDouble (PyTuple_GetItem (linear, 0));
  l [1] = PyFloat_AsDouble (PyTuple_GetItem (linear, 1));
  l [2] = PyFloat_AsDouble (PyTuple_GetItem (linear, 2));

  a [0] = PyFloat_AsDouble (PyTuple_GetItem (angular, 0));
  a [1] = PyFloat_AsDouble (PyTuple_GetItem (angular, 1));
  a [2] = PyFloat_AsDouble (PyTuple_GetItem (angular, 2));

  BODY_Initial_Velocity (body->bod, l, a);

#if MPI
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

#if MPI
  if (ID_TO_BODY (body))
  {
#endif

  BODY_Material (body->bod, volid, get_bulk_material (solfec->sol, material));

#if MPI
  }
#endif

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
    if (ID_TO_BODY (body))
    {
#endif

    DOM_Remove_Body (solfec->sol->dom, body->bod);
    body->dodestroy = 1; /* set destruction flag as this body will no longer be destroyed by the domain */

#if MPI
    }
#endif
  }
  else
  {
    constraint = (lng_CONSTRAINT*)object;

#if MPI
    if (ID_TO_CONSTRAINT (solfec, constraint))
    {
#endif

    DOM_Remove_Constraint (solfec->sol->dom, constraint->con);

#if MPI
    }
#endif
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
    s [0] = PyFloat_AsDouble (PyTuple_GetItem (shape, 0));
    s [1] = PyFloat_AsDouble (PyTuple_GetItem (shape, 1));
    s [2] = PyFloat_AsDouble (PyTuple_GetItem (shape, 2));

    COPY (v, t); 
    NORMALIZE (t);
    SCALE (t, (ALG_PI * angle / 180.0));
    EXPMAP (t, r);
    NVMUL (r, s, t);

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
      }
    }
  }

  Py_RETURN_NONE;
}

/* split shape by plane */
static PyObject* lng_SPLIT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("shape", "point", "normal");
  PyObject *shape, *point, *normal, *out, *b, *f;
  SHAPE *shp, *shq, *back, *front;
  double p [3], n [3];
  int error;

  PARSEKEYS ("OOO", &shape, &point, &normal);

  TYPETEST (is_shape (shape, kwl[0]) && is_tuple (point, kwl[1], 3) && is_tuple (normal, kwl[2], 3));

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  n [0] = PyFloat_AsDouble (PyTuple_GetItem (normal, 0));
  n [1] = PyFloat_AsDouble (PyTuple_GetItem (normal, 1));
  n [2] = PyFloat_AsDouble (PyTuple_GetItem (normal, 2));

  shp = create_shape (shape, 1); /* empty */

  for (shq = shp; shq; shq = shq->next)
  {
    if (shq->kind == SHAPE_MESH)
    {
      PyErr_SetString (PyExc_ValueError, "MESH object cannot be SPLIT");
      return NULL;
    }
  }

  back = front = NULL;

  TRY ()
  {
    for (shq = shp; shq; shq = shq->next)
    {
      if (shq->kind == SHAPE_CONVEX)
      {
	CONVEX *one = NULL, *two = NULL;

	CONVEX_Split (shq->data, p, n, &one, &two);
	if (one) back = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, one), back);
	if (two) front = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, two), front);
      }
      else
      {
	SPHERE *one = NULL, *two = NULL;

	SPHERE_Split (shq->data, p, n, &one, &two);
	if (one) back = SHAPE_Glue (SHAPE_Create (SHAPE_SPHERE, one), back);
	if (two) front = SHAPE_Glue (SHAPE_Create (SHAPE_SPHERE, two), front);
      }
    }
  }
  CATCHANY (error)
  {
    PyErr_SetString (PyExc_RuntimeError, errstring (error));
    return NULL;
  }
  ENDTRY ()

  SHAPE_Destroy (shp);

  b =  shape_to_list (back);
  f =  shape_to_list (front);

  if (!(b||f)) return NULL;

  out = Py_BuildValue ("(O, O)", b, f);

  if (back) SHAPE_Destroy_Wrapper (back);
  if (front) SHAPE_Destroy_Wrapper (front);

  return out;
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
    }
  }

  out = shape_to_list (shp);

  SHAPE_Destroy_Wrapper (shp);

  return out;
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
#if MPI
    if (ID_TO_BODY ((lng_BODY*)shape))
    {
#endif

    BODY *bod = ((lng_BODY*)shape)->bod;
    SHAPE_Char (bod->shape, &v, c, e);

#if MPI
    }
#endif
  }
  else /* shape */
  {
    shp = create_shape (shape, 0);
    SHAPE_Char (shp, &v, c, e);
    SHAPE_Destroy_Wrapper (shp);
  }

  return Py_BuildValue ("(d, d, d)", c[0], c[1], c[2]);
}

/* exclude a pair of bodies from contact detection */
static PyObject* lng_CONTACT_EXCLUDE_BODIES (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "body1", "body2");
  lng_BODY *body1, *body2;
  lng_SOLFEC *solfec;

  PARSEKEYS ("OOO", &solfec, &body1, &body2);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body1, kwl[1]) && is_body (body2, kwl[2]));

#if MPI
  AABB_Exclude_Body_Pair (solfec->sol->aabb, body1->id, body2->id);
#else
  AABB_Exclude_Body_Pair (solfec->sol->aabb, body1->bod->id, body2->bod->id);
#endif

  Py_RETURN_NONE;
}

/* exclude a pair of geometric objects from contact detection */
static PyObject* lng_CONTACT_EXCLUDE_OBJECTS (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "body1", "point1", "body2", "point2");
  PyObject *point1, *point2;
  lng_BODY *body1, *body2;
  double p1 [3], p2 [3];
  lng_SOLFEC *solfec;
  int sgp1, sgp2;

  PARSEKEYS ("OOOOO", &solfec, &body1, &point1, &body2, &point2);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_body (body1, kwl[1]) &&
    is_tuple (point1, kwl[2], 3) && is_body (body2, kwl[3]) && is_tuple (point2, kwl[4], 3));

#if MPI
  if (ID_TO_BODY (body1) && ID_TO_BODY (body2))
  {
#endif

  p1 [0] = PyFloat_AsDouble (PyTuple_GetItem (point1, 0));
  p1 [1] = PyFloat_AsDouble (PyTuple_GetItem (point1, 1));
  p1 [2] = PyFloat_AsDouble (PyTuple_GetItem (point1, 2));

  p2 [0] = PyFloat_AsDouble (PyTuple_GetItem (point2, 0));
  p2 [1] = PyFloat_AsDouble (PyTuple_GetItem (point2, 1));
  p2 [2] = PyFloat_AsDouble (PyTuple_GetItem (point2, 2));

  sgp1 = SHAPE_Sgp (body1->bod->sgp, body1->bod->nsgp, p1);

  if (sgp1 < 0)
  {
    PyErr_SetString (PyExc_ValueError, "First point outside of body one");
    return NULL;
  }

  sgp2 = SHAPE_Sgp (body2->bod->sgp, body2->bod->nsgp, p2);

  if (sgp2 < 0)
  {
    PyErr_SetString (PyExc_ValueError, "Second point outside of body two");
    return NULL;
  }

  AABB_Exclude_Gobj_Pair (solfec->sol->aabb, body1->bod->id, sgp1, body2->bod->id, sgp2);

#if MPI
  }
#endif

  Py_RETURN_NONE;
}

/* set contact sparsification threshold */
static PyObject* lng_CONTACT_SPARSIFY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "threshold");
  lng_SOLFEC *solfec;
  double threshold;

  PARSEKEYS ("Od", &solfec, &threshold);

  TYPETEST (is_solfec (solfec, kwl[0]));

  if (threshold < 0 || threshold > 1.0)
  {
    PyErr_SetString (PyExc_ValueError, "threshold outside of [0, 1] bounds");
    return NULL;
  }

  solfec->sol->dom->threshold = threshold;

  Py_RETURN_NONE;
}

/* test whether an object is a constraint solver */
static int is_solver (PyObject *obj, char *var)
{
  if (obj)
  {
    if (!PyObject_IsInstance (obj, (PyObject*)&lng_GAUSS_SEIDEL_SOLVER_TYPE) &&
        !PyObject_IsInstance (obj, (PyObject*)&lng_EXPLICIT_SOLVER_TYPE))
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
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_EXPLICIT_SOLVER_TYPE))
    return EXPLICIT_SOLVER;
  else return -1;
}

/* return solver interface */
static void* get_solver (PyObject *obj)
{
  if (PyObject_IsInstance (obj, (PyObject*)&lng_GAUSS_SEIDEL_SOLVER_TYPE))
    return ((lng_GAUSS_SEIDEL_SOLVER*)obj)->gs;
  else if (PyObject_IsInstance (obj, (PyObject*)&lng_EXPLICIT_SOLVER_TYPE))
    return NULL;
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
  if (!RND_On ()) /* otherwise interactive run is controlled by the viewer */
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
  CMP_ALG cmp;

  compression = NULL;
  cmp = CMP_OFF;

  PARSEKEYS ("Od|O", &solfec, &interval, &compression);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_positive (interval, kwl[1]) && is_string (compression, kwl [2]));

  if (solfec->sol->mode == SOLFEC_READ) Py_RETURN_NONE; /* skip READ mode */

  if (compression)
  {
    IFIS (compression, "OFF")
    {
      cmp = CMP_OFF;
    }
    ELIF (compression, "FASTLZ")
    {
      cmp = CMP_FASTLZ;
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

  PARSEKEYS ("OO", &body, &point);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3));

#if MPI
  if (ID_TO_BODY (body))
  {
#endif

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  BODY_Point_Values (body->bod, p, VALUE_DISPLACEMENT, x);

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

  PARSEKEYS ("OO", &body, &point);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3));

#if MPI
  if (ID_TO_BODY (body))
  {
#endif

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  BODY_Point_Values (body->bod, p, VALUE_VELOCITY, x);

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

  PARSEKEYS ("OO", &body, &point);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3));

#if MPI
  if (ID_TO_BODY (body))
  {
#endif

  p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
  p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
  p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

  BODY_Point_Values (body->bod, p, VALUE_STRESS_AND_MISES, x);

  return Py_BuildValue ("(d, d, d, d, d, d, d)", x[0], x[1], x[2], x[3], x[4], x[5], x[6]);

#if MPI
  }
  else Py_RETURN_NONE;
#endif
}

/* callback data */
typedef struct
{
  BODY *bod;
  double p [3];
  VALUE_KIND kind;
  int index;
  PyObject *times;
  PyObject *values;
}
lng_HISTORY_DATA;

/* callback for SOLFEC_History */
static void lng_HISTORY_callback (lng_HISTORY_DATA *data, double time)
{
  double x [6];

  BODY_Point_Values (data->bod, data->p, data->kind, x);

  PyList_Append (data->times, PyFloat_FromDouble (time));
  PyList_Append (data->values, PyFloat_FromDouble (x [data->index]));
}

/* history of an entity */
static PyObject* lng_HISTORY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("body", "point", "entity", "t0", "t1");
  double start, end, t0, t1;
  PyObject *point, *entity;
  lng_HISTORY_DATA data;
  lng_BODY *body;

  PARSEKEYS ("OOOdd", &body, &point, &entity, &t0, &t1);

  TYPETEST (is_body (body, kwl[0]) && is_tuple (point, kwl[1], 3) && is_string (entity, kwl[2]));

#if MPI
  if (ID_TO_BODY (body))
  {
#endif

  SOLFEC *sol = DOM(body->bod->dom)->owner;

  if (sol->mode == SOLFEC_WRITE) Py_RETURN_NONE;
  else
  {
    SOLFEC_Time_Limits (sol, &start, &end);

    if (t0 < start || t1 > end)
    {
      PyErr_SetString (PyExc_ValueError, "t0 or t1 outside of simulation duration limits");
      return NULL;
    }

    IFIS (entity, "DX")
    {
      data.kind = VALUE_DISPLACEMENT;
      data.index = 0;
    }
    ELIF (entity, "DY")
    {
      data.kind = VALUE_DISPLACEMENT;
      data.index = 1;
    }
    ELIF (entity, "DZ")
    {
      data.kind = VALUE_DISPLACEMENT;
      data.index = 2;
    }
    ELIF (entity, "VX")
    {
      data.kind = VALUE_VELOCITY;
      data.index = 0;
    }
    ELIF (entity, "VY")
    {
      data.kind = VALUE_VELOCITY;
      data.index = 1;
    }
    ELIF (entity, "VZ")
    {
      data.kind = VALUE_VELOCITY;
      data.index = 2;
    }
    ELIF (entity, "SX")
    {
      data.kind = VALUE_STRESS;
      data.index = 0;
    }
    ELIF (entity, "SY")
    {
      data.kind = VALUE_STRESS;
      data.index = 1;
    }
    ELIF (entity, "SZ")
    {
      data.kind = VALUE_STRESS;
      data.index = 2;
    }
    ELIF (entity, "SXY")
    {
      data.kind = VALUE_STRESS;
      data.index = 3;
    }
    ELIF (entity, "SXZ")
    {
      data.kind = VALUE_STRESS;
      data.index = 4;
    }
    ELIF (entity, "SYZ")
    {
      data.kind = VALUE_STRESS;
      data.index = 5;
    }
    ELIF (entity, "MISES")
    {
      data.kind = VALUE_MISES;
      data.index = 0;
    }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid entity kind");
      return NULL;
    }

    data.bod = body->bod;

    data.p [0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    data.p [1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    data.p [2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    ERRMEM (data.times = PyList_New (0));
    ERRMEM (data.values = PyList_New (0));

    SOLFEC_History (sol, NULL, NULL, NULL, 0, body->bod, NULL,
      t0, t1, &data, (void (*) (void*, double)) lng_HISTORY_callback);

    return Py_BuildValue ("(O,O)", data.times, data.values);
  }

#if MPI
  }
  else Py_RETURN_NONE;
#endif
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

/* callback data */
typedef struct
{
  double val;
  PyObject *times;
  PyObject *values;
}
lng_TIMING_HISTORY_DATA;

/* callback for SOLFEC_History */
static void lng_TIMING_HISTORY_callback (lng_TIMING_HISTORY_DATA *data, double time)
{
  PyList_Append (data->times, PyFloat_FromDouble (time));
  PyList_Append (data->values, PyFloat_FromDouble (data->val));
}

/* timing history */
static PyObject* lng_TIMING_HISTORY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("solfec", "kind", "t0", "t1");
  lng_TIMING_HISTORY_DATA data;
  double start, end, t0, t1;
  lng_SOLFEC *solfec;
  PyObject *kind;
  char *label;

  PARSEKEYS ("OOdd", &solfec, &kind, &t0, &t1);

  TYPETEST (is_solfec (solfec, kwl[0]) && is_string (kind, kwl[1]));

  if (solfec->sol->mode == SOLFEC_WRITE) Py_RETURN_NONE;
  else
  {
    SOLFEC_Time_Limits (solfec->sol, &start, &end);

    if (t0 < start || t1 > end)
    {
      PyErr_SetString (PyExc_ValueError, "t0 or t1 outside of simulation duration limits");
      return NULL;
    }

    label = PyString_AsString (kind);

    if (!SOLFEC_Has_Timer (solfec->sol, label))
    {
      PyErr_SetString (PyExc_ValueError, "Invalid timing kind");
      return NULL;
    }

    ERRMEM (data.times = PyList_New (0));
    ERRMEM (data.values = PyList_New (0));

    SOLFEC_History (solfec->sol, label, &data.val, NULL, 1, NULL, NULL,
      t0, t1, &data, (void (*) (void*, double)) lng_TIMING_HISTORY_callback);

    return Py_BuildValue ("(O,O)", data.times, data.values);
  }
}

static PyMethodDef lng_methods [] =
{
  {"HULL", (PyCFunction)lng_HULL, METH_VARARGS|METH_KEYWORDS, "Create convex hull from a point set"},
  {"HEX", (PyCFunction)lng_HEX, METH_VARARGS|METH_KEYWORDS, "Create mesh of a hexahedral shape"},
  {"FIX_POINT", (PyCFunction)lng_FIX_POINT, METH_VARARGS|METH_KEYWORDS, "Create a fixed point constraint"},
  {"FIX_DIRECTION", (PyCFunction)lng_FIX_DIRECTION, METH_VARARGS|METH_KEYWORDS, "Create a fixed direction constraint"},
  {"SET_DISPLACEMENT", (PyCFunction)lng_SET_DISPLACEMENT, METH_VARARGS|METH_KEYWORDS, "Create a prescribed displacement constraint"},
  {"SET_VELOCITY", (PyCFunction)lng_SET_VELOCITY, METH_VARARGS|METH_KEYWORDS, "Create a prescribed velocity constraint"},
  {"SET_ACCELERATION", (PyCFunction)lng_SET_ACCELERATION, METH_VARARGS|METH_KEYWORDS, "Create a prescribed acceleration constraint"},
  {"PUT_RIGID_LINK", (PyCFunction)lng_PUT_RIGID_LINK, METH_VARARGS|METH_KEYWORDS, "Create a rigid linke constraint"},
  {"GRAVITY", (PyCFunction)lng_GRAVITY, METH_VARARGS|METH_KEYWORDS, "Set gravity acceleration"},
  {"FORCE", (PyCFunction)lng_FORCE, METH_VARARGS|METH_KEYWORDS, "Apply point force"},
  {"TORQUE", (PyCFunction)lng_TORQUE, METH_VARARGS|METH_KEYWORDS, "Apply point torque"},
  {"IMBALANCE_TOLERANCES", (PyCFunction)lng_IMBALANCE_TOLERANCES, METH_VARARGS|METH_KEYWORDS, "Adjust imbalance tolerances"},
  {"LOCDYN_BALANCING", (PyCFunction)lng_LOCDYN_BALANCING, METH_VARARGS|METH_KEYWORDS, "Set local dynamics balancing mode"},
  {"NCPU", (PyCFunction)lng_NCPU, METH_NOARGS, "Get the number of processors"},
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
  {"COPY", (PyCFunction)lng_COPY, METH_VARARGS|METH_KEYWORDS, "Copy shape"},
  {"BYLABEL", (PyCFunction)lng_BYLABEL, METH_VARARGS|METH_KEYWORDS, "Get object by label"},
  {"MASS_CENTER", (PyCFunction)lng_MASS_CENTER, METH_VARARGS|METH_KEYWORDS, "Get mass center"},
  {"CONTACT_EXCLUDE_BODIES", (PyCFunction)lng_CONTACT_EXCLUDE_BODIES, METH_VARARGS|METH_KEYWORDS, "Exclude body pair from contact detection"},
  {"CONTACT_EXCLUDE_OBJECTS", (PyCFunction)lng_CONTACT_EXCLUDE_OBJECTS, METH_VARARGS|METH_KEYWORDS, "Exclude geometric object pair from contact detection"},
  {"CONTACT_SPARSIFY", (PyCFunction)lng_CONTACT_SPARSIFY, METH_VARARGS|METH_KEYWORDS, "Adjust contact sparsification"},
  {"RUN", (PyCFunction)lng_RUN, METH_VARARGS|METH_KEYWORDS, "Run analysis"},
  {"OUTPUT", (PyCFunction)lng_OUTPUT, METH_VARARGS|METH_KEYWORDS, "Set data output interval"},
  {"EXTENTS", (PyCFunction)lng_EXTENTS, METH_VARARGS|METH_KEYWORDS, "Set scene extents"},
  {"CALLBACK", (PyCFunction)lng_CALLBACK, METH_VARARGS|METH_KEYWORDS, "Set analysis callback"},
  {"UNPHYSICAL_PENETRATION", (PyCFunction)lng_UNPHYSICAL_PENETRATION, METH_VARARGS|METH_KEYWORDS, "Set analysis callback"},
  {"DURATION", (PyCFunction)lng_DURATION, METH_VARARGS|METH_KEYWORDS, "Get analysis duration"},
  {"FORWARD", (PyCFunction)lng_FORWARD, METH_VARARGS|METH_KEYWORDS, "Set forward in READ mode"},
  {"BACKWARD", (PyCFunction)lng_BACKWARD, METH_VARARGS|METH_KEYWORDS, "Set backward in READ mode"},
  {"SEEK", (PyCFunction)lng_SEEK, METH_VARARGS|METH_KEYWORDS, "Seek to time in READ mode"},
  {"DISPLACEMENT", (PyCFunction)lng_DISPLACEMENT, METH_VARARGS|METH_KEYWORDS, "Get displacement of a referential point"},
  {"VELOCITY", (PyCFunction)lng_VELOCITY, METH_VARARGS|METH_KEYWORDS, "Get velocity of a referential point"},
  {"STRESS", (PyCFunction)lng_STRESS, METH_VARARGS|METH_KEYWORDS, "Get stress of a referential point"},
  {"HISTORY", (PyCFunction)lng_HISTORY, METH_VARARGS|METH_KEYWORDS, "Get history of an entity for a body"},
  {"TIMING", (PyCFunction)lng_TIMING, METH_VARARGS|METH_KEYWORDS, "Get timing"},
  {"TIMING_HISTORY", (PyCFunction)lng_TIMING_HISTORY, METH_VARARGS|METH_KEYWORDS, "Get timing history"},
  {NULL, 0, 0, NULL}
};

/* 
 * initialization
 */

static void initlng (void)
{
  PyObject *m;

  TYPEINIT (lng_CONVEX_TYPE, lng_CONVEX, "solfec.CONVEX",
    Py_TPFLAGS_DEFAULT, lng_CONVEX_dealloc, lng_CONVEX_new,
    lng_CONVEX_methods, lng_CONVEX_members, lng_CONVEX_getset);

  TYPEINIT (lng_MESH_TYPE, lng_MESH, "solfec.MESH",
    Py_TPFLAGS_DEFAULT, lng_MESH_dealloc, lng_MESH_new,
    lng_MESH_methods, lng_MESH_members, lng_MESH_getset);

  TYPEINIT (lng_SPHERE_TYPE, lng_SPHERE, "solfec.SPHERE",
    Py_TPFLAGS_DEFAULT, lng_SPHERE_dealloc, lng_SPHERE_new,
    lng_SPHERE_methods, lng_SPHERE_members, lng_SPHERE_getset);

  TYPEINIT (lng_SOLFEC_TYPE, lng_SOLFEC, "solfec.SOLFEC",
    Py_TPFLAGS_DEFAULT, lng_SOLFEC_dealloc, lng_SOLFEC_new,
    lng_SOLFEC_methods, lng_SOLFEC_members, lng_SOLFEC_getset);

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

  TYPEINIT (lng_EXPLICIT_SOLVER_TYPE, lng_EXPLICIT_SOLVER, "solfec.EXPLICIT_SOLVER",
    Py_TPFLAGS_DEFAULT, lng_EXPLICIT_SOLVER_dealloc, lng_EXPLICIT_SOLVER_new,
    lng_EXPLICIT_SOLVER_methods, lng_EXPLICIT_SOLVER_members, lng_EXPLICIT_SOLVER_getset);

  TYPEINIT (lng_CONSTRAINT_TYPE, lng_CONSTRAINT, "solfec.CONSTRAINT",
    Py_TPFLAGS_DEFAULT, lng_CONSTRAINT_dealloc, lng_CONSTRAINT_new,
    lng_CONSTRAINT_methods, lng_CONSTRAINT_members, lng_CONSTRAINT_getset);

  if (PyType_Ready (&lng_CONVEX_TYPE) < 0) return;
  if (PyType_Ready (&lng_MESH_TYPE) < 0) return;
  if (PyType_Ready (&lng_SPHERE_TYPE) < 0) return;
  if (PyType_Ready (&lng_SOLFEC_TYPE) < 0) return;
  if (PyType_Ready (&lng_SURFACE_MATERIAL_TYPE) < 0) return;
  if (PyType_Ready (&lng_BULK_MATERIAL_TYPE) < 0) return;
  if (PyType_Ready (&lng_BODY_TYPE) < 0) return;
  if (PyType_Ready (&lng_TIME_SERIES_TYPE) < 0) return;
  if (PyType_Ready (&lng_GAUSS_SEIDEL_SOLVER_TYPE) < 0) return;
  if (PyType_Ready (&lng_EXPLICIT_SOLVER_TYPE) < 0) return;
  if (PyType_Ready (&lng_CONSTRAINT_TYPE) < 0) return;

  if (!(m =  Py_InitModule3 ("solfec", lng_methods, "Solfec module"))) return;

  Py_INCREF (&lng_CONVEX_TYPE);
  Py_INCREF (&lng_MESH_TYPE);
  Py_INCREF (&lng_SPHERE_TYPE);
  Py_INCREF (&lng_SOLFEC_TYPE);
  Py_INCREF (&lng_SURFACE_MATERIAL_TYPE);
  Py_INCREF (&lng_BULK_MATERIAL_TYPE);
  Py_INCREF (&lng_BODY_TYPE);
  Py_INCREF (&lng_TIME_SERIES_TYPE);
  Py_INCREF (&lng_GAUSS_SEIDEL_SOLVER_TYPE);
  Py_INCREF (&lng_EXPLICIT_SOLVER_TYPE);
  Py_INCREF (&lng_CONSTRAINT_TYPE);

  PyModule_AddObject (m, "CONVEX", (PyObject*)&lng_CONVEX_TYPE);
  PyModule_AddObject (m, "MESH", (PyObject*)&lng_MESH_TYPE);
  PyModule_AddObject (m, "SPHERE", (PyObject*)&lng_SPHERE_TYPE);
  PyModule_AddObject (m, "SOLFEC", (PyObject*)&lng_SOLFEC_TYPE);
  PyModule_AddObject (m, "SURFACE_MATERIAL", (PyObject*)&lng_SURFACE_MATERIAL_TYPE);
  PyModule_AddObject (m, "BULK_MATERIAL", (PyObject*)&lng_BULK_MATERIAL_TYPE);
  PyModule_AddObject (m, "BODY", (PyObject*)&lng_BODY_TYPE);
  PyModule_AddObject (m, "TIME_SERIES", (PyObject*)&lng_TIME_SERIES_TYPE);
  PyModule_AddObject (m, "GAUSS_SEIDEL_SOLVER", (PyObject*)&lng_GAUSS_SEIDEL_SOLVER_TYPE);
  PyModule_AddObject (m, "EXPLICIT_SOLVER", (PyObject*)&lng_EXPLICIT_SOLVER_TYPE);
  PyModule_AddObject (m, "CONSTRAINT", (PyObject*)&lng_CONSTRAINT_TYPE);
}

/* 
 * input file parsing
 */

/* interpret an input file (return 0 on success)*/
int lng (const char *path)
{
  FILE *file;
  int error;

  ASSERT (file = fopen (path, "r"), ERR_FILE_OPEN);

  Py_Initialize();

  initlng ();

  PyRun_SimpleString("from solfec import CONVEX\n"
                     "from solfec import HULL\n"
                     "from solfec import MESH\n"
                     "from solfec import HEX\n"
                     "from solfec import SPHERE\n"
                     "from solfec import SOLFEC\n"
                     "from solfec import SURFACE_MATERIAL\n"
                     "from solfec import BULK_MATERIAL\n"
                     "from solfec import BODY\n"
                     "from solfec import TIME_SERIES\n"
                     "from solfec import GAUSS_SEIDEL_SOLVER\n"
                     "from solfec import EXPLICIT_SOLVER\n"
                     "from solfec import FIX_POINT\n"
                     "from solfec import FIX_DIRECTION\n"
                     "from solfec import SET_DISPLACEMENT\n"
                     "from solfec import SET_VELOCITY\n"
                     "from solfec import SET_ACCELERATION\n"
                     "from solfec import PUT_RIGID_LINK\n"
                     "from solfec import GRAVITY\n"
                     "from solfec import FORCE\n"
                     "from solfec import TORQUE\n"
                     "from solfec import IMBALANCE_TOLERANCES\n"
                     "from solfec import LOCDYN_BALANCING\n"
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
                     "from solfec import COPY\n"
                     "from solfec import BYLABEL\n"
                     "from solfec import MASS_CENTER\n"
                     "from solfec import CONTACT_EXCLUDE_BODIES\n"
                     "from solfec import CONTACT_EXCLUDE_OBJECTS\n"
                     "from solfec import CONTACT_SPARSIFY\n"
                     "from solfec import RUN\n"
                     "from solfec import OUTPUT\n"
                     "from solfec import EXTENTS\n"
                     "from solfec import CALLBACK\n"
                     "from solfec import UNPHYSICAL_PENETRATION\n"
                     "from solfec import DURATION\n"
                     "from solfec import FORWARD\n"
                     "from solfec import BACKWARD\n"
                     "from solfec import SEEK\n"
                     "from solfec import DISPLACEMENT\n"
                     "from solfec import VELOCITY\n"
                     "from solfec import STRESS\n"
                     "from solfec import HISTORY\n"
                     "from solfec import TIMING\n"
                     "from solfec import TIMING_HISTORY\n");

  error = PyRun_SimpleFile (file, path);

  fclose (file);

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
