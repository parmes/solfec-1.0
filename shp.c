/*
 * shp.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * shape implementation
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

#include <float.h>
#include "sol.h"
#include "alg.h"
#include "shp.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "err.h"
#include "pck.h"

/* specific shape interface hooks */
typedef void* (*copy_func) (void*);
static copy_func copy [] = {(copy_func)MESH_Copy, (copy_func)CONVEX_Copy, (copy_func)SPHERE_Copy};
typedef void (*adjup_func) (void*);
static adjup_func adjup [] = {(adjup_func)MESH_Update_Adjacency, (adjup_func)CONVEX_Update_Adjacency, (adjup_func)SPHERE_Update_Adjacency};
typedef void (*scale_func) (void*, double*);
static scale_func scale [] = {(scale_func)MESH_Scale, (scale_func)CONVEX_Scale, (scale_func)SPHERE_Scale};
typedef void (*translate_func) (void*, double*);
static translate_func translate [] = {(translate_func)MESH_Translate, (translate_func)CONVEX_Translate, (translate_func)SPHERE_Translate};
typedef void (*rotate_func) (void*, double*, double*, double);
static rotate_func rotate [] = {(rotate_func)MESH_Rotate, (rotate_func)CONVEX_Rotate, (rotate_func)SPHERE_Rotate};
typedef TRI* (*cut_func) (void*, double*, double*, int*);
static cut_func cut [] = {(cut_func)MESH_Cut, (cut_func)CONVEX_Cut, (cut_func)SPHERE_Cut};
typedef void (*gcha_func) (void*, int, double*, double*, double*, double*, double*);
static gcha_func gcha [] = {(gcha_func)MESH_Char_Partial, (gcha_func)CONVEX_Char_Partial, (gcha_func)SPHERE_Char_Partial};
typedef void* (*gobj_func) (void*, double*);
static gobj_func gobj [] = {(gobj_func)MESH_Element_Containing_Spatial_Point, (gobj_func)CONVEX_Containing_Point, (gobj_func)SPHERE_Containing_Point};
typedef int (*gobjs_func) (void*, void*, double*);
static gobjs_func gobjs [] = {(gobjs_func)ELEMENT_Contains_Spatial_Point, (gobjs_func)CONVEX_Contains_Point, (gobjs_func)SPHERE_Contains_Point};
typedef double (*gobjdst_func) (void*, void*, double*);
static gobjdst_func gobjdst [] = {(gobjdst_func)ELEMENT_Spatial_Point_Distance, (gobjdst_func)CONVEX_Spatial_Point_Distance, (gobjdst_func)SPHERE_Spatial_Point_Distance};
typedef void (*update_func) (void*, void*, void*, MOTION);
static update_func update [] = {(update_func)MESH_Update, (update_func)CONVEX_Update, (update_func)SPHERE_Update};
typedef void (*extents_func) (void*, double*);
static extents_func objextents [] = {(extents_func)MESH_Extents, (extents_func)CONVEX_List_Extents, (extents_func)SPHERE_List_Extents};
typedef void (*oextents_func) (void*, double*, double*, double*, double*);
static oextents_func objorientedextents [] = {(oextents_func)MESH_Oriented_Extents, (oextents_func)CONVEX_List_Oriented_Extents, (oextents_func)SPHERE_List_Oriented_Extents};
typedef void* (*first_bulk_func) (void*);
static first_bulk_func firstbulk [] = {(first_bulk_func)MESH_First_Bulk_Material, (first_bulk_func)CONVEX_First_Bulk_Material, (first_bulk_func)SPHERE_First_Bulk_Material};
typedef void (*destroy_func) (void*);
static destroy_func destroy [] = {(destroy_func)MESH_Destroy, (destroy_func)CONVEX_Destroy, (destroy_func)SPHERE_Destroy};
typedef void (*pack_func) (void*, int*, double**, int*, int*, int**, int*);
static pack_func pack [] = {(pack_func)MESH_Pack, (pack_func)CONVEX_Pack, (pack_func)SPHERE_Pack};
typedef void* (*unpack_func) (void*, int*, double*, int, int*, int*, int);
static unpack_func unpack [] = {(unpack_func)MESH_Unpack, (unpack_func)CONVEX_Unpack, (unpack_func)SPHERE_Unpack};

/* append shape */
static SHAPE* append (SHAPE *list, short kind, void *data)
{
  SHAPE *shp;

  shp = SHAPE_Create (kind, data);
  shp->next = list;

  return shp;
}

/* create a general shape */
SHAPE* SHAPE_Create (short kind, void *data)
{
  SHAPE *shp;

  ERRMEM (shp = malloc (sizeof (SHAPE)));
  shp->kind = kind;
  shp->data = data;
  shp->next = NULL;

  return shp;
}

/* create shape geometric object pairs */
SGP* SGP_Create (SHAPE *shp, int *nsgp, SGP_FLAGS flags)
{
  SGP *sgp, *ptr;
  int n = 0;

  /* compute geomerical objects */
  for (SHAPE *shq = shp; shq; shq = shq->next)
  {
    switch (shq->kind)
    {
      case SHAPE_MESH:
      {
	MESH *msh = shq->data;
	for (ELEMENT *ele = msh->surfeles; ele; ele = ele->next) n ++;
	if (flags & SGP_MESH_NODES) n += msh->surfnodes_count;
      }
      break;
      case SHAPE_CONVEX:
      {
	CONVEX *cvx = shq->data;
	for (; cvx; cvx = cvx->next) n ++;
      }
      break;
      case SHAPE_SPHERE:
      {
	SPHERE *sph = shq->data;
	for (; sph; sph = sph->next) n ++;
      }
      break;
    }
  }

  /* allocate */
  ERRMEM (ptr = sgp = MEM_CALLOC (n * sizeof (SGP)));
  *nsgp = n;

  /* set pointers */
  for (SHAPE *shq = shp; shq; shq = shq->next)
  {
    switch (shq->kind)
    {
      case SHAPE_MESH:
      {
	MESH *msh = shq->data;
	for (ELEMENT *ele = msh->surfeles; ele; ele = ele->next, ptr ++)
	  ptr->shp = shq, ptr->gobj = ele, ptr->kind = GOBJ_ELEMENT;

	if (flags & SGP_MESH_NODES)
	{
	  for (NODE *nod = msh->surfnodes; nod; nod = nod->n, ptr ++)
	    ptr->shp = shq, ptr->gobj = nod, ptr->kind = GOBJ_NODE;
	}
      }
      break;
      case SHAPE_CONVEX:
      {
	CONVEX *cvx = shq->data;
	for (; cvx; cvx = cvx->next, ptr ++) ptr->shp = shq, ptr->gobj = cvx, ptr->kind = GOBJ_CONVEX;
      }
      break;
      case SHAPE_SPHERE:
      {
	SPHERE *sph = shq->data;
	for (; sph; sph = sph->next, ptr ++) ptr->shp = shq, ptr->gobj = sph, ptr->kind = GOBJ_SPHERE;
      }
      break;
    }
  }

  return sgp;
}

/* SGP to GOBJ conversion for FEM purposes */
void* SGP_2_GOBJ (SGP *sgp)
{
  switch (sgp->kind)
  {
  case GOBJ_NODE: return ((NODE*)sgp->gobj)->fac [0]->ele;
  default: return sgp->gobj;
  }
}

/* get GOBJ type of given shape */
GOBJ SHAPE_2_GOBJ (SHAPE *shp)
{
  switch (shp->kind)
  {
  case SHAPE_MESH: return GOBJ_ELEMENT;
  case SHAPE_CONVEX: return GOBJ_CONVEX;
  case SHAPE_SPHERE: return GOBJ_SPHERE;
  }

  ASSERT_TEXT (0, "Unknown shape kind");

  return 0;
}

/* copy shape */
SHAPE* SHAPE_Copy (SHAPE *shp)
{
  SHAPE *out, *q;
  void *data;

  for (out = NULL; shp; shp = shp->next)
  {
    data = copy [shp->kind] (shp->data);
    q = SHAPE_Create (shp->kind, data);
    q->next = out;
    out = q;
  }

  return out;
}

/* glue two shape lists (gluing together basic shapes) */
SHAPE* SHAPE_Glue (SHAPE *shp, SHAPE *shq)
{
  SHAPE *out, *next;
  CONVEX *cvx;
  SPHERE *sph;

  for (out = NULL, cvx = NULL, sph = NULL; shp; shp = next)
  {
    next = shp->next;

    switch (shp->kind)
    {
    case SHAPE_MESH:
      shp->next = out; /* meshes are copied */
      out = shp;
      break;
    case SHAPE_CONVEX:
      cvx = CONVEX_Glue (shp->data, cvx); /* convices are lumped together */
      free (shp);
      break;
    case SHAPE_SPHERE:
      sph = SPHERE_Glue (shp->data, sph); /* spheres are also lumped */
      free (shp);
      break;
    }
  }

  for (; shq; shq = next)
  {
    next = shq->next;

    switch (shq->kind)
    {
    case SHAPE_MESH:
      shq->next = out; /* meshes are copied */
      out = shq;
      break;
    case SHAPE_CONVEX:
      cvx = CONVEX_Glue (shq->data, cvx); /* convices are lumped together */
      free (shq);
      break;
    case SHAPE_SPHERE:
      sph = SPHERE_Glue (shq->data, sph); /* spheres are also lumped */
      free (shq);
      break;
    }
  }

  if (cvx) out = append (out, SHAPE_CONVEX, cvx); /* append with convices */

  if (sph) out = append (out, SHAPE_SPHERE, sph); /* append with spheres */

  return out;
}

/* glue two shape lists (without gluing basic shapes) */
SHAPE* SHAPE_Glue_Simple (SHAPE *shp, SHAPE *shq)
{
  SHAPE *shr = shp;

  for (; shp->next; shp = shp->next);
  shp->next = shq;
  return shr;
}

/* update adjacency data of stored shapes;
 * no other function affects the adjacency */
void SHAPE_Update_Adjacency (SHAPE *shp)
{
  for (; shp; shp = shp->next)
    adjup [shp->kind] (shp->data);
}

/* scale cur shape => 
 * if MESH,  scale each: x *= vector [0], y *= vector [1], z *= vector [2];
 * if SPHERE, scale radius: r *= vector [0]; (set ref = cur) */
void SHAPE_Scale (SHAPE *shp, double *vector)
{
  for (; shp; shp = shp->next)
    scale [shp->kind] (shp->data, vector);
}

/* translate cur shape (set ref = cur) */
void SHAPE_Translate (SHAPE *shp, double *vector)
{
  for (; shp; shp = shp->next)
    translate [shp->kind] (shp->data, vector);
}

/* rotate cur shape (set ref = cur), around the line (point, vector) */
void SHAPE_Rotate (SHAPE *shp, double *point, double *vector, double angle)
{
  for (; shp; shp = shp->next)
    rotate [shp->kind] (shp->data, point, vector, angle);
}

/* cut through shape with a plane; return triangulated cross-section; all returned data
 * points to the memory allocated after the triangles memory; adjacency is not maintained;
 * TRI->adj[0] stores a pointer to the geometrical object that has been cut by the triangle;
 * (cur_to_ref, sgp, ref, cur, n) can be either all NULL or all valid pointers; if not NULL then
 * 'n' reference and current vertices are calculated (triagnle vertices are the current ones) */
TRI* SHAPE_Cut (SHAPE *shp, double *point, double *normal, int *m,
  void *body, MOTION cur_to_ref, SGP **sgp, double **ref, double **cur, int *n)
{
  TRI **tri, *out, *t, *e, *q;
  int *ntr, l, s, i, j, k;
  MAP *vertices, *item;
  MEM mapmem, setmem;
  double *v;
  SET *set;

  ASSERT_DEBUG ((body && cur_to_ref && sgp && ref && cur && n) ||
                (!body && !cur_to_ref && !sgp && !ref && !cur && !n),
		"cur_to_ref, sgp, ref, cur, n must all be either NULL or valid pointers");

  if (shp->next == NULL && cur_to_ref == NULL) /* one subshape case */
  {
    return cut [shp->kind] (shp->data, point, normal, m);
  }

  j = 0;
  l = 0;
  s = 128;
  (*m) = 0;
  set = NULL;
  out = NULL;
  vertices = NULL;
  MEM_Init (&mapmem, sizeof (MAP), 128);
  MEM_Init (&setmem, sizeof (SET), 128);
  ERRMEM (tri = malloc (s * sizeof (TRI*)));
  ERRMEM (ntr = malloc (s * sizeof (int)));

  for (; shp; shp = shp->next)
  {
    tri [l] = cut [shp->kind] (shp->data, point, normal, &ntr [l]); /* find intersection */

    if (tri [l]) /* found */
    {
      for (t = tri [l], e = t + ntr [l]; t != e; t ++)
      {
	for (i = 0; i < 3; i++)
	{
	  if (!MAP_Find_Node (vertices, t->ver [i], NULL)) /* not mapped yet */
	  {
	    item = MAP_Insert (&mapmem, &vertices, t->ver [i], (void*) (long) j, NULL); /* map it */
	    ASSERT_DEBUG (item, "Failed to map vertex during shape cutting");
	    if (item) j ++;
	  }
	}
	t->adj [1] = (TRI*) shp; /* record the corresponding shape */
      }
      (*m) += ntr [l];
      l ++;
    }

    if (l == s)
    {
      s *= 2;
      ERRMEM (tri = realloc (tri, s * sizeof (TRI*)));
      ERRMEM (ntr = realloc (ntr, s * sizeof (int)));
    }
  }

  if (l == 0) goto out;

  i = MAP_Size (vertices);

  if (cur_to_ref)
  {
    ERRMEM (out = malloc ((*m) * sizeof (TRI) + 2 * i * sizeof (double [3]) + i * sizeof (SGP)));
    v = (double*) (out + (*m));
    (*sgp) = (SGP*) (v + 6 * i);
    (*ref) = v + 3 * i;
    (*cur) = v;
    (*n) = i;
  }
  else
  {
    ERRMEM (out = malloc ((*m) * sizeof (TRI) + i * sizeof (double [3])));
    v = (double*) (out + (*m));
  }

  /* copy vertices  */
  for (item = MAP_First (vertices); item; item = MAP_Next (item))
  {
    i = (int)(long)item->data;
    double *a = item->key, *b = &v [3 * i];
    COPY (a, b);
  }

  /* compile output */
  for (i = 0, q = out; i < l; i ++)
  {
    for (t = tri [i], e = t + ntr [i]; t != e; t ++, q ++)
    {
      q->adj [0] = t->adj [0];
      COPY (t->out, q->out);

      for (k = 0; k < 3; k ++)
      {
	item = MAP_Find_Node (vertices, t->ver [k], NULL);
	ASSERT_DEBUG (item, "Invalid vertex mapping during shape cutting");
	j = 3*((int) (long) item->data);
	q->ver [k] = &v [j]; /* map to new storage */

        /* map referential to current if needed,
	 * but without repititions; set up SGP pair */
	if (cur_to_ref && !SET_Contains (set, item, NULL))
	{
	  double *Ver = &(*ref) [j];
	  SGP *sg = &(*sgp) [j/3];
	  sg->shp = (SHAPE*)t->adj [1];
	  sg->gobj = t->adj [0];
	  sg->kind = SHAPE_2_GOBJ (sg->shp);
	  cur_to_ref (body, sg, q->ver [k], Ver);
	  SET_Insert (&setmem, &set, item, NULL);
	}
      }
    }
  }
out:
  for (i = 0; i < l; i ++) free (tri [i]);
  MEM_Release (&mapmem);
  MEM_Release (&setmem);
  free (tri);
  free (ntr);

  return out;
}

/* split shape by plane; output two parts of the split shape;
 * topoadj != 0 implies cutting from the point and through the topological adjacency only */
void SHAPE_Split (SHAPE *shp, double *point, double *normal, short topoadj, int surfid, SHAPE **one, SHAPE **two)
{
  SHAPE *shq, *back, *front;

  back = front = NULL;

  for (shq = shp; shq; shq = shq->next)
  {
    if (shq->kind == SHAPE_CONVEX)
    {
      CONVEX *one = NULL, *two = NULL;

      CONVEX_Split (shq->data, point, normal, topoadj, surfid, &one, &two);
      if (one) back = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, one), back);
      if (two) front = SHAPE_Glue (SHAPE_Create (SHAPE_CONVEX, two), front);
    }
    else if (shq->kind == SHAPE_SPHERE)
    {
      SPHERE *one = NULL, *two = NULL;

      SPHERE_Split (shq->data, point, normal, topoadj, surfid, &one, &two);
      if (one) back = SHAPE_Glue (SHAPE_Create (SHAPE_SPHERE, one), back);
      if (two) front = SHAPE_Glue (SHAPE_Create (SHAPE_SPHERE, two), front);
    }
    else if (shq->kind == SHAPE_MESH)
    {
      MESH *one = NULL, *two = NULL;

      MESH_Split (shq->data, point, normal, topoadj, surfid, &one, &two);
      if (one) back = SHAPE_Glue (SHAPE_Create (SHAPE_MESH, one), back);
      if (two) front = SHAPE_Glue (SHAPE_Create (SHAPE_MESH, two), front);
    }
  }

  *one = back;
  *two = front;
}

/* get spatial/referential characteristics => volume, mass center, and Euler tensor (centered) */
void SHAPE_Char (SHAPE *shp, int ref, double *volume, double *center, double *euler)
{
  double vo, sx, sy, sz,
	 cen [3], eul [9];

  vo = sx = sy = sz = 0.0;
  SET9 (eul, 0.0);

  for (; shp; shp = shp->next)
    gcha [shp->kind] (shp->data, ref, &vo, &sx, &sy, &sz, eul);

  cen [0] = sx / vo;
  cen [1] = sy / vo;
  cen [2] = sz / vo;

  eul [0] -= (2*sx - cen[0]*vo)*cen[0];
  eul [4] -= (2*sy - cen[1]*vo)*cen[1];
  eul [8] -= (2*sz - cen[2]*vo)*cen[2];
  eul [3] -= cen[0]*sy + cen[1]*sx - cen[0]*cen[1]*vo;
  eul [6] -= cen[0]*sz + cen[2]*sx - cen[0]*cen[2]*vo;
  eul [7] -= cen[1]*sz + cen[2]*sy - cen[1]*cen[2]*vo;
  eul [1] = eul[3];
  eul [2] = eul[6];
  eul [5] = eul[7];

  if (volume) *volume = vo;
  if (center) COPY (cen, center);
  if (euler) NNCOPY (eul, euler);
}

/* for the given shape (not a list) compute spatial/referential partial characteristic: 'vo'lume and
 * static momenta 'sx', 'sy, 'sz' and 'eul'er tensor; assume that all input data is initially zero; */
void SHAPE_Char_Partial (SHAPE *shp, int ref, double *vo, double *sx, double *sy, double *sz, double *eul)
{
  double v = 0, s[3] = {0, 0, 0}, e [9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  gcha [shp->kind] (shp->data, ref, &v, &s[0], &s[1], &s[2], e);

  if (vo) *vo += v;
  if (sx) *sx += s [0];
  if (sy) *sy += s [1];
  if (sz) *sz += s [2];
  if (eul)
  {
    eul [0] += e[0];
    eul [4] += e[4];
    eul [8] += e[8];
    eul [3] += e[3];
    eul [6] += e[6];
    eul [7] += e[7];
    eul [1] += e[3];
    eul [2] += e[6];
    eul [5] += e[7];
  }
}

/* return an object containing spatial point */
void* SHAPE_Gobj (SHAPE *shp, double *point, SHAPE **out)
{
  void *obj;

  /* TODO: optimize by using a deformable spatial tree */

  for (obj = NULL; shp; shp = shp->next)
  {
    obj = gobj [shp->kind] (shp->data, point);

    if (obj)
    {
      if (out) *out = shp;
      break;
    }
  }

  return obj;
}

/* return an index of object containing spatial point (or -1 on failure) */
int SHAPE_Sgp (SGP *sgp, int nsgp, double *point)
{
  int i;

  /* TODO: optimize by using a deformable spatial tree */

  for (i = 0; i < nsgp; i ++, sgp ++)
  {
    if (sgp->kind & GOBJ_NODE) continue; /* XXX: skip nodes */

    if (gobjs [sgp->shp->kind] (sgp->shp->data, sgp->gobj, point)) return i;
  }

  return -1;
}

/* return an index of object closest to the spatial point (output distance in 'd' if not NULL) */
int SHAPE_Closest_Sgp (SGP *sgp, int nsgp, double *point, double *d)
{
  double dst, dstmin;
  int i, j;

  /* TODO: optimize by using a deformable spatial tree */

  for (i = j = 0, dstmin = DBL_MAX; i < nsgp; i ++, sgp ++)
  {
    if (sgp->kind == GOBJ_NODE) /* XXX: future cases (e.g. edges) need to be handled */
    {
      NODE *nod = sgp->gobj;
      double a [3];
      SUB (point, nod->cur, a);
      dst = LEN (a);
    }
    else dst = gobjdst [sgp->shp->kind] (sgp->shp->data, sgp->gobj, point);

    if (dst < dstmin)
    {
      dstmin = dst;
      j = i;
    }
  }

  if (d) *d = dstmin;

  return j;
}

/* update current shape with given motion */
void SHAPE_Update (SHAPE *shp, void *body, MOTION motion)
{
  for (; shp; shp = shp->next)
    update [shp->kind] (shp->data, body, shp, motion);
}

/* copute shape extents */
void SHAPE_Extents (SHAPE *shp, double *extents)
{
  double e [6], margin;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;

  for (; shp; shp = shp->next)
  {
    objextents [shp->kind] (shp->data, e);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [3] > extents [3]) extents [3] = e [3];
    if (e [4] > extents [4]) extents [4] = e [4];
    if (e [5] > extents [5]) extents [5] = e [5];
  }

  margin = 10.0 * GEOMETRIC_EPSILON;
  extents [0] -= margin;
  extents [1] -= margin;
  extents [2] -= margin;
  extents [3] += margin;
  extents [4] += margin;
  extents [5] += margin;
}

/* copute shape oriented extents in corrds given by three direction vectors */
void SHAPE_Oriented_Extents (SHAPE *shp, double *vx, double *vy, double *vz, double *extents)
{
  double e [6], margin;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;

  for (; shp; shp = shp->next)
  {
    objorientedextents [shp->kind] (shp->data, vx, vy, vz, e);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [3] > extents [3]) extents [3] = e [3];
    if (e [4] > extents [4]) extents [4] = e [4];
    if (e [5] > extents [5]) extents [5] = e [5];
  }

  margin = 10.0 * GEOMETRIC_EPSILON;
  extents [0] -= margin;
  extents [1] -= margin;
  extents [2] -= margin;
  extents [3] += margin;
  extents [4] += margin;
  extents [5] += margin;
}

/* return first bulk material recorded
 * in this individual shape (not a list) */
void* SHAPE_First_Bulk_Material (SHAPE *shp)
{
  return firstbulk [shp->kind] (shp->data);
}

/* release shape memory */
void SHAPE_Destroy (SHAPE *shp)
{
  SHAPE *next;

  for (; shp; shp = next)
  {
    next = shp->next;
    destroy [shp->kind] (shp->data);
    free (shp);
  }
}

/* release shape wrapper memory (without data) */
void SHAPE_Destroy_Wrapper (SHAPE *shp)
{
  SHAPE *next;

  for (; shp; shp = next)
  {
    next = shp->next;
    free (shp);
  }
}

/* pack shape */
void SHAPE_Pack (SHAPE *shp, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  SHAPE *ptr;
  int count;

  for (count = 0, ptr = shp; ptr; ptr = ptr->next) count ++;

  pack_int (isize, i, ints, count);

  for (ptr = shp; ptr; ptr = ptr->next)
  {
    pack_int (isize, i, ints, ptr->kind);
    pack [ptr->kind] (ptr->data, dsize, d, doubles, isize, i, ints);
  }
}

/* unpack shape */
SHAPE* SHAPE_Unpack (void *solfec, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int count;
  SHAPE *tail = NULL,
	*head,
	*ptr;

  count = unpack_int (ipos, i, ints);

  for (; count > 0; count --)
  {
    int kind = unpack_int (ipos, i, ints);

    ptr = SHAPE_Create (kind, unpack [kind] (solfec, dpos, d, doubles, ipos, i, ints));

    if (tail)
    {
      tail->next = ptr;
      tail = ptr;
    }
    else head = tail = ptr;
  }

  return head;
}

/* export MBFCP definition */
void SHAPE_2_MBFCP (SHAPE *shp, FILE *out)
{
  SET *msh = NULL, *item;
  SHAPE *ptr;
  int n;

  for (ptr = shp, n = 0; ptr; ptr = ptr->next, n ++);

  fprintf (out, "SHAPES:\t%d\n", n);

  for (ptr = shp; ptr; ptr = ptr->next)
  {
    switch (ptr->kind)
    {
    case SHAPE_MESH:
      SET_Insert (NULL, &msh, ptr->data, NULL);
      break;
    case SHAPE_CONVEX:
      CONVEX_2_MBFCP (ptr->data, out);
      break;
    case SHAPE_SPHERE:
      SPHERE_2_MBFCP (ptr->data, out);
      break;
    }
  }

  if (SET_Size (msh))
  {
    fprintf (out, "MESHES:\t%d\n", SET_Size (msh));

    for (item = SET_First (msh); item; item = SET_Next (item))
    {
      MESH_2_MBFCP (item->data, out);
    }

    SET_Free (NULL, &msh);
  }
}
