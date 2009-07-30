/*
 * box.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * axis aligned bounding box overlap detection
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

#include "box.h"
#include "hyb.h"
#include "alg.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "bod.h"
#include "swp.h"
#include "hsh.h"
#include "err.h"

#define SIZE 512

/* auxiliary data */
struct auxdata
{
  SET *nobody;
  SET *nogobj;
  MEM *mapmem;
  void *data;
  BOX_Overlap_Create create;
};

/* body pair comparison */
static int bodcmp (OPR *a, OPR *b)
{
  if (a->bod1 < b->bod1) return -1;
  else if (a->bod1 == b->bod1 && a->bod2 < b->bod2) return -1;
  else if (a->bod1 == b->bod1 && a->bod2 == b->bod2) return 0;
  else return 1;
}

/* geomtric object pair comparison */
static int gobjcmp (OPR *a, OPR *b)
{
  int cmp = bodcmp (a, b);

  if (cmp == 0)
  {
    if (a->sgp1 < b->sgp1) return -1;
    else if (a->sgp1 == b->sgp1 && a->sgp2 < b->sgp2) return -1;
    else if (a->sgp1 == b->sgp1 && a->sgp2 == b->sgp2) return 0;
    else return 1;
  }
  else return cmp;
}

/* check for a lack of overlap */
static int no_overlap (double *a, double *b)
{
  if (a[0] > b[3] || b[0] > a[3]) return 1;
  if (a[1] > b[4] || b[1] > a[4]) return 1;
  if (a[2] > b[5] || b[2] > a[5]) return 1;
  return 0;
}

/* local overlap creation callback => filters our unwnated adjacency */
static void* local_create (struct auxdata *aux, BOX *one, BOX *two)
{
  void *user;

  /* check if these are two obstacles => no need to report overlap */
  if (BODY(one->body)->kind == OBS && BODY(two->body)->kind == OBS) return NULL;

  /* boxes already adjacent ? */
  if (MAP_Find (one->adj, two, NULL)) return NULL;

  /* self-contact ? <=> one->body == two->body */
  if (one->body == two->body)
  {
    if (!(((BODY*)one->body)->flags & BODY_DETECT_SELF_CONTACT)) return NULL;
    /* it's much more efficient to handle this case separately with an on-body
     * flag rathert than use the exclusion set within the AABB structure */
  }

  int id1 = BODY (one->body)->id,
      id2 = BODY (two->body)->id,
      no1 = one->sgp - BODY (one->body)->sgp,
      no2 = two->sgp - BODY (two->body)->sgp;

  OPR pair = {MIN (id1, id2), MAX (id1, id2), MIN (no1, no2), MAX (no1, no2)};

  /* an excluded body pair ? */
  if (SET_Find (aux->nobody, &pair, (SET_Compare) bodcmp)) return NULL;

  /* an excluded object pair ? */
  if (SET_Find (aux->nogobj, &pair, (SET_Compare) gobjcmp)) return NULL;

  /* test topological adjacency */
  switch (GOBJ_Pair_Code (one, two))
  {
    case AABB_ELEMENT_ELEMENT:
    {
      if (one->body == two->body) /* assuming only one mesh per body is possible,
				     the two elements are within the same mesh;
				     exclude topologically adjacent elements */
	if (ELEMENT_Adjacent (one->sgp->gobj, two->sgp->gobj)) return NULL;
    }
    break;
    case AABB_CONVEX_CONVEX:
    {
      if (one->body == two->body) /* exclude topologically adjacent convices */
	if (CONVEX_Adjacent (one->sgp->gobj, two->sgp->gobj)) return NULL;
    }
    break;
    case AABB_SPHERE_SPHERE:
    {
      if (one->body == two->body) /* exclude topologically adjacent spheres */
	if (SPHERE_Adjacent (one->sgp->gobj, two->sgp->gobj)) return NULL;
    }
  }
 
  /* report overlap creation and get the user pointer */
  user = aux->create (aux->data, one, two);

  if (user) /* if the user pointer is NULL keep re-detecting the overlap */
  {
    MAP_Insert (aux->mapmem, &one->adj, two, user, NULL); /* map overlap symmetrically */
    MAP_Insert (aux->mapmem, &two->adj, one, user, NULL);
    /* one might think that only one adjacency (say box with bigger pointer
     * stores box with smaller pointer) could be used. nevertheless, for box
     * removal one needs to release all overlaps, and then the complete adjacency
     * informaltion is necessary; thus the symmetrical mapping */
  }

  return user;
}

/* create box overlap driver data */
AABB* AABB_Create (int size)
{
  AABB *aabb;

  ERRMEM (aabb = malloc (sizeof (AABB)));
  aabb->lst = aabb->in = aabb->out= NULL;
  aabb->tab = NULL;

  MEM_Init (&aabb->boxmem, sizeof (BOX), MIN (size, SIZE));
  MEM_Init (&aabb->mapmem, sizeof (MAP), MIN (size, SIZE));
  MEM_Init (&aabb->setmem, sizeof (SET), MIN (size, SIZE));
  MEM_Init (&aabb->oprmem, sizeof (OPR), MIN (size, SIZE));

  aabb->nobody = aabb->nogobj= NULL;

  aabb->nlst = 0;
  aabb->ntab = 0;
  aabb->modified = 0;
  aabb->swp = aabb->hsh = NULL;

  return aabb;
}

/* insert geometrical object */
BOX* AABB_Insert (AABB *aabb, void *body, GOBJ kind, SGP *sgp, void *data, BOX_Extents_Update update)
{
  BOX *box;

  ERRMEM (box = MEM_Alloc (&aabb->boxmem));
  box->data = data;
  box->update = update;
  box->kind = kind;
  box->body = body;
  box->sgp = sgp;

  /* include into the
   * insertion list */
  box->next = aabb->in;
  if (aabb->in) aabb->in->prev= box;
  aabb->in = box;

  /* set modified */
  aabb->modified = 1;

  return box;
}

/* delete an object */
void AABB_Delete (AABB *aabb, BOX *box)
{
  /* remove box from the current list */
  if (box->prev) box->prev->next = box->next;
  else aabb->lst = box->next;
  if (box->next) box->next->prev = box->prev;
  aabb->nlst --; /* decrease count */

  /* insert it into
   * the deletion list */
  box->next = aabb->out;
  aabb->out = box;

  /* set modified */
  aabb->modified = 1;
}

/* update state => detect created and released overlaps */
void AABB_Update (AABB *aabb, BOXALG alg, void *data, BOX_Overlap_Create create, BOX_Overlap_Release release)
{
  struct auxdata aux = {aabb->nobody, aabb->nogobj, &aabb->mapmem, data, create};
  BOX *box, *next, *adj;
  MAP *item, *jtem;
  int nin;

  /* update extents */
  for (box = aabb->lst; box; box = box->next) /* for each current box */
    box->update (box->data, box->sgp->gobj, box->extents);

  for (box = aabb->in, nin = 0; box; box = box->next, nin ++) /* for each inserted box */
    box->update (box->data, box->sgp->gobj, box->extents);

  /* check for released overlaps */
  for (box = aabb->out; box; box = next) /* for each deleted box */
  {
    next = box->next;

    for (item = MAP_First (box->adj); item; item = MAP_Next (item)) /* for each adjacent box */
    {
      adj = item->key;
      release (data, box, adj, item->data); /* release overlap */
      MAP_Delete (&aabb->mapmem, &adj->adj, box, NULL); /* remove 'box' from 'adj's adjacency */
    }

    MAP_Free (&aabb->mapmem, &box->adj); /* free adjacency */
    MEM_Free (&aabb->boxmem, box); /* free box */
  }
  aabb->out = NULL; /* list emptied */

  for (box = aabb->lst; box; box = box->next) /* for each current box */
  {
    for (item = MAP_First (box->adj); item; item = jtem) /* for each adjacent box */
    {
      adj = item->key;
      if (no_overlap (box->extents, adj->extents))
      {
	release (data, box, adj, item->data); /* release overlap */
	MAP_Delete (&aabb->mapmem, &adj->adj, box, NULL); /* remove 'box' from 'adj's adjacency */
	jtem = MAP_Delete_Node (&aabb->mapmem, &box->adj, item); /* remove 'adj' from 'box's adjacency */
      }
      else jtem = MAP_Next (item);
    }
  }

  if (aabb->modified) /* merge insertion and curent lists, update pointer table */
  {
    BOX **b;

    if (aabb->in)
    {
      for (box = aabb->in; box->next; box = box->next); /* rewind till the last inserted one */
      box->next = aabb->lst; /* link it with the current list */
      if (aabb->lst) aabb->lst->prev = box; /* set previous link */
      aabb->lst = aabb->in; /* list appended */
      aabb->in = NULL; /* list emptied */
      aabb->nlst += nin; /* increase number of current boxes */
      aabb->ntab = aabb->nlst; /* resize the pointer table */
      ERRMEM (aabb->tab = realloc (aabb->tab, sizeof (BOX*) * aabb->ntab));
    }

    for (box = aabb->lst, b = aabb->tab; box; box = box->next, b ++) *b = box; /* overwrite box pointers */
  }

  /* the algorithm
   * specific part */
  switch (alg)
  {
    case HYBRID:
    {
      hybrid (aabb->tab, aabb->ntab, &aux, (BOX_Overlap_Create)local_create);
    }
    break;
    case HASH3D:
    {
      if (!aabb->hsh) aabb->hsh = HASH_Create (aabb->ntab);

      HASH_Do (aabb->hsh, aabb->ntab, aabb->tab,
	       &aux, (BOX_Overlap_Create)local_create); 
    }
    break;
    case SWEEP_HASH2D_LIST:
    case SWEEP_HASH1D_XYTREE:
    case SWEEP_HASH2D_XYTREE:
    case SWEEP_XYTREE:
    {
      if (!aabb->swp) aabb->swp = SWEEP_Create (aabb->ntab, (DRALG)alg);

      SWEEP_Do (aabb->swp, (DRALG)alg, aabb->modified, aabb->ntab,
	        aabb->tab, &aux, (BOX_Overlap_Create)local_create); 
    }
    break;
  }

  /* set unmodified */
  aabb->modified = 0;
}

/* never report overlaps betweem this pair of bodies (given by identifiers) */
void AABB_Exclude_Body_Pair (AABB *aabb, unsigned int id1, unsigned int id2)
{
  OPR *opr;
  
  ERRMEM (opr = MEM_Alloc (&aabb->oprmem));
  opr->bod1 = MIN (id1, id2);
  opr->bod2 = MAX (id1, id2);
  SET_Insert (&aabb->setmem, &aabb->nobody, opr, (SET_Compare) bodcmp);
}

/* never report overlaps betweem this pair of objects (bod1, sgp1), (bod1, sgp2) */
void AABB_Exclude_Gobj_Pair (AABB *aabb, unsigned int bod1, int sgp1, unsigned int bod2, int sgp2)
{
  OPR *opr;
  
  ERRMEM (opr = MEM_Alloc (&aabb->oprmem));
  opr->bod1 = MIN (bod1, bod2);
  opr->bod2 = MAX (bod1, bod2);
  opr->sgp1 = MIN (sgp1, sgp2);
  opr->sgp2 = MAX (sgp1, sgp2);
  SET_Insert (&aabb->setmem, &aabb->nogobj, opr, (SET_Compare) gobjcmp);
}

/* break box adjacency if boxes associated with the objects are adjacent */
void AABB_Break_Adjacency (AABB *aabb, BOX *one, BOX *two)
{
  /* mutually delete from adjacency lists */
  MAP_Delete (&aabb->mapmem, &one->adj, two, NULL);
  MAP_Delete (&aabb->mapmem, &two->adj, one, NULL);
}

/* release memory */
void AABB_Destroy (AABB *aabb)
{
  free (aabb->tab);
  MEM_Release (&aabb->boxmem);
  MEM_Release (&aabb->mapmem);
  MEM_Release (&aabb->setmem);
  MEM_Release (&aabb->oprmem);
  if (aabb->swp) SWEEP_Destroy (aabb->swp);
  if (aabb->hsh) HASH_Destroy (aabb->hsh);
  free (aabb);
}
