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
#include "fem.h"
#include "pck.h"

/* pseudo-rigid body cracking */
static CRACK* prb_crack (BODY *bod, BODY **one, BODY **two)
{
  double values [6], cauchy [9], vec [3], tension;
  CRACK *cra;

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
      BODY_Split (bod, cra->point, cra->normal, cra->topoadj, cra->surfid, one, two);

      /* TODO: energy decrease in the parts */

      return cra;
    }
  }

  return NULL;
}

/* cut through the mesh and create element and referential point pairs */
static ELEPNT* element_point_pairs (MESH *msh, double *point, double *normal, short topoadj, int *nepn)
{
  ELEPNT *epn;
  int i, n;
  TRI *tri;

  tri = MESH_Ref_Cut (msh, point, normal, nepn); /* referential cut */

  if (topoadj)
  {
    TRI_Compadj (tri, *nepn);
    ASSERT_TEXT (tri = TRI_Topoadj (tri, *nepn, point, nepn), "Failed to cut through the mesh wiht the crack plane ");
  }
  
  ERRMEM (epn = MEM_CALLOC (sizeof (ELEPNT [*nepn])));
  for (n = 0; n < *nepn; n ++)
  {
    epn [n].ele = tri [n].ptr;
    for (i = 0; i < 3; i ++)
    {
      epn [n].pnt [i] = (tri [n].ver [0][i] + tri [n].ver [1][i] + tri [n].ver [2][i]) / 3.0;
    }
  }

  free (tri);
  return epn;
}

/* finite element body cracking */
static CRACK* fem_crack (BODY *bod, BODY **one, BODY **two)
{
  double values [6], cauchy [9], vec [3], tension;
  ELEPNT *epn;
  CRACK *cra;
  int i;

  *one = *two = NULL;

  for (cra = bod->cra; cra; cra = cra->next)
  {
    if (!cra->epn) cra->epn = element_point_pairs (bod->msh ?  bod->msh : bod->shape->data,
	                                 cra->point, cra->normal, cra->topoadj, &cra->nepn);

    for (i = 0, epn = cra->epn; i < cra->nepn; i ++, epn ++)
    {
      FEM_Point_Values (bod, epn->ele, epn->pnt, VALUE_STRESS, values);
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
	BODY_Split (bod, cra->point, cra->normal, cra->topoadj, cra->surfid, one, two);

	/* TODO: energy decrease in the parts */

	return cra;
      }
    }
  }

  return NULL;
}

/* remap contraints to new bodies */
static void remap_constraints (DOM *dom, BODY *bod, CRACK *cra, BODY *one, BODY *two)
{
  SET *item, *con1, *con2;
  double a [3], d;
  CON *con;

  for (con1 = con2 = NULL, item = SET_First (bod->con); item; item = SET_Next (item))
  {
    con = item->data;
    if (con->kind == RIGLNK && bod == con->slave)
    {
      double s [3], *z = RIGLNK_VEC (con->Z);
      SUB (con->point, z, s);
      SUB (s, cra->point, a);
    }
    else
    {
      SUB (con->point, cra->point, a);
    }

    d = DOT (cra->normal, a);
    if (d <= 0.0) SET_Insert (NULL, &con1, con, NULL);
    else SET_Insert (NULL, &con2, con, NULL);
  }

  for (item = SET_First (con1); item; item = SET_Next (item))
  {
    DOM_Transfer_Constraint (dom, item->data, bod, one);
  }

  for (item = SET_First (con2); item; item = SET_Next (item))
  {
    DOM_Transfer_Constraint (dom, item->data, bod, two ? two : one); /* two is NULL for half-cracks */
  }

  SET_Free (NULL, &con1);
  SET_Free (NULL, &con2);
}

/* copy crack data */
static void copy_crack (CRACK *src, CRACK *dst)
{
  *dst = *src; /* XXX: maintain correctnes when CRACK stores more data */
  dst->epn = NULL;
  dst->nepn = 0;
}

/* create crack object */
CRACK* CRACK_Create ()
{
  return MEM_CALLOC (sizeof (CRACK));
}

/* delete crack object */
void CRACK_Destroy (CRACK *cra)
{
  free (cra->epn);
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

/* pack cracks */
void CRACKS_Pack (CRACK *list, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  CRACK *cra;
  int n;

  for (n = 0, cra = list; cra; cra = cra->next) n ++;

  pack_int (isize, i, ints, n);

  for (cra = list; cra; cra = cra->next)
  {
    pack_doubles (dsize, d, doubles, cra->point, 6);
    pack_int (isize, i, ints, cra->crit);
    pack_int (isize, i, ints, cra->surfid);
    pack_double (dsize, d, doubles, cra->ft);
    pack_double (dsize, d, doubles, cra->Gf);
  }
}

/* unpack cracks */
CRACK* CRACKS_Unpack (int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  CRACK *out, *cra;
  int j, n;

  n = unpack_int (ipos, i, ints);

  for (j = 0, out = NULL; j < n; j ++)
  {
    cra = CRACK_Create ();
    unpack_doubles (dpos, d, doubles, cra->point, 6);
    cra->crit = unpack_int (ipos, i, ints);
    cra->surfid = unpack_int (ipos, i, ints);
    cra->ft = unpack_double (dpos, d, doubles);
    cra->Gf = unpack_double (dpos, d, doubles);
    cra->next = out;
    out = cra;
  }

  return out;
}

/* propagate cracks and adjust the domain */
void Propagate_Cracks (DOM *dom)
{
  BODY *bod, *one, *two, *next;
  SET *item, *newbod = NULL;
  CRACK *cra, *crb, *c;

  for (bod = dom->bod; bod; bod = next)
  {
    next = bod->next;
    one = two = NULL;
    cra = NULL;

    if (bod->cra) /* try to crack bodies */
    {
      switch (bod->kind)
      {
      case PRB: cra = prb_crack (bod, &one, &two); break;
      case FEM: cra = fem_crack (bod, &one, &two); break;
      default: break;
      }
    }

    if (cra && one && two) /* cracked and fragmented in two */
    {
      for (crb = bod->cra; crb; crb = crb->next) /* copy remaining cracks */
      {
	if (crb == cra) continue; /* skip current */

	if (SHAPE_Cut_Possible (one->shape, 1, crb->point, crb->normal, crb->topoadj)) /* crack overlaps first fragment */
	{
	  c = CRACK_Create();
	  copy_crack (crb, c);
	  c->next = one->cra;
	  one->cra = c;
	}

	if (SHAPE_Cut_Possible (two->shape, 1, crb->point, crb->normal, crb->topoadj)) /* crack overlaps second fragment */
	{
	  c = CRACK_Create();
	  copy_crack (crb, c);
	  c->next = two->cra;
	  two->cra = c;
	}
      }

      remap_constraints (dom, bod, cra, one, two); /* remap or delete constraints */

      DOM_Remove_Body (dom, bod); /* remove from domain */
      BODY_Destroy (bod); /* destroys cracks store at this body */

      if (SHAPE_Separable (one->shape))
      {
         /* FIXME => TODO */
      }
      else
      {
        if (one->shape->nbradj) one->flags |= BODY_DETECT_SELF_CONTACT; /* FIXME => remove when the above is made work */

	SET_Insert (NULL, &newbod, one, NULL);
      }

      if (SHAPE_Separable (two->shape))
      {
         /* FIXME => TODO */
      }
      else 
      {
        if (two->shape->nbradj) two->flags |= BODY_DETECT_SELF_CONTACT; /* FIXME => remove when the above is made work */

	SET_Insert (NULL, &newbod, two, NULL);
      }
    }
    else if (cra && cra->topoadj && one) /* half-crack with new boundary */
    {
      for (crb = bod->cra; crb; crb = crb->next) /* copy remaining cracks */
      {
	if (crb == cra) continue; /* skip current */

	if (SHAPE_Cut_Possible (one->shape, 1, crb->point, crb->normal, crb->topoadj)) /* crack overlaps new body */
	{
	  c = CRACK_Create();
	  copy_crack (crb, c);
	  c->next = one->cra;
	  one->cra = c;
	}
      }

      remap_constraints (dom, bod, cra, one, NULL); /* remap or delete constraints */

      DOM_Remove_Body (dom, bod); /* remove from domain */
      BODY_Destroy (bod); /* destroys cracks store at this body */

      if (SHAPE_Separable (one->shape))
      {
	/* FIXME => TODO */
      }
      else
      {
	one->flags |= BODY_DETECT_SELF_CONTACT; /* enable self-contact checking */

	SET_Insert (NULL, &newbod, one, NULL);
      }
    }
    else
    {
      ASSERT_TEXT (cra == NULL && one == NULL && two == NULL,
        "A body cracked, but body splitting has failed.\n"
	"Adjust GEOMETRIC_EPSILON or slightly shift the crack plane.");
    }
  }

  for (item = SET_First (newbod); item; item = SET_Next (item)) /* insert new bodies into domain */
  {
#if MPI
    DOM_Pending_Body (dom, item->data);
#else
    DOM_Insert_Body (dom, item->data);
#endif
  }

  SET_Free (NULL, &newbod);
}
