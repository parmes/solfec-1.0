/*
 * ldy.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * the local dynamic problem
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

#include "alg.h"
#include "dom.h"
#include "ldy.h"
#include "lap.h"
#include "err.h"

#if MPI
#include <limits.h>
#include "pck.h"
#include "com.h"
#include "tag.h"
#endif

/* memory block size */
#define BLKSIZE 512

/* apply forward change of variables (nornal
 * contact forces) due to the cohesion, etc. */
static void variables_change_begin (LOCDYN *ldy)
{
  OFFB *blk;
  DIAB *dia;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;
    double *B = dia->B; /* free velocity will
			   be eventually modified */

    if (con->state & CON_COHESIVE) /* cohesive state */
    {
      double c = con->mat.cohesion * con->area,
	     *W = dia->W,
	     *R = dia->R;

      R [2] += c;       /* R_n_new = R_n + c <=> (R_n + c) >= 0 */
      B[0] -= (W[6]*c); /* in consequnce 'W_tn * c' gets subtracted */
      B[1] -= (W[7]*c); /* ... */
      B[2] -= (W[8]*c); /* and 'W_nn * c' here */
    }

    /* off-diagonal subtractions */
    for (blk = dia->adj; blk; blk = blk->n)
    {
      CON *con = blk->dia->con;
      if (con->state & CON_COHESIVE) /* cohesive state */
      {
	double c = con->mat.cohesion * con->area,
	       *W = blk->W;

	B[0] -= (W[6]*c);
	B[1] -= (W[7]*c);
	B[2] -= (W[8]*c);
      }
    }
  }
}

/* apply back change of variables (nornal
 * contact forces) due to the cohesion, etc. */
static void variables_change_end (LOCDYN *ldy)
{
  DIAB *dia;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;
    short state = con->state;

    if (state & CON_COHESIVE) /* cohesive state */
    {
      double c = con->mat.cohesion * con->area,
	     *R = dia->R;

      R [2] -= c; /* back change */

      if ((state & CON_OPEN) || /* mode-I decohesion */
	(!(state & CON_STICK))) /* mode-II decohesion */
      {
	con->state &= ~CON_COHESIVE;
	con->mat.cohesion = 0.0;
      }
    }
  }
}

#if MPI
/* swap two ids */
inline static void ID_swap (ZOLTAN_ID_PTR a, ZOLTAN_ID_PTR b)
{
  ZOLTAN_ID_TYPE c = *a; *a = *b; *b = c;
}

/* sort a list of identifiers */
static void ID_sort (ZOLTAN_ID_PTR begin, ZOLTAN_ID_PTR end)
{
  ZOLTAN_ID_TYPE *lower = begin,
	         *upper = end,
                  bound = *(begin+(end-begin)/2);
  
  while (lower <= upper)
  {
    while (*lower < bound) lower++;

    while (bound < *upper) upper--;

    if (lower < upper) ID_swap (lower ++, upper--);
    else lower ++;
  }

  if (begin < upper) ID_sort (begin, upper);
  if (upper < end) ID_sort (upper+1, end);
}

/* number of vertices in the local W graph */
static int vertex_count (LOCDYN *ldy, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return ldy->ndiab + ldy->nins; /* balanced + newly inserted */
}

/* list of vertices in the local W graph */
static void vertex_list (LOCDYN *ldy, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr)
{
  DIAB *dia;
  int i, j;

  for (dia = ldy->diab, i = 0; dia; i ++, dia = dia->n) /* balanced vertices */
  {
    global_ids [i * num_gid_entries] = dia->id;
    local_ids [i * num_lid_entries] = UINT_MAX; /* mapped */
  }

  for (j = 0; j < ldy->nins; i ++, j ++) /* newly inserted vertices */
  {
    global_ids [i * num_gid_entries] = ldy->ins [j]->id;
    local_ids [i * num_lid_entries] = j; /* indexed */
  }

  *ierr = ZOLTAN_OK;
}

/* sizes of edges (compressed rows) of the local W graph */
static void edge_sizes (LOCDYN *ldy, int *num_lists, int *num_pins, int *format, int *ierr)
{
  DIAB *dia;
  OFFB *b;
  int j, n = 0;

  for (dia = ldy->diab; dia; dia = dia->n, n ++)
    for (b = dia->adj; b; b = b->n) n ++;

  for (j = 0; j < ldy->nins; j ++, n ++)
  {
    dia = ldy->ins [j];
    for (b = dia->adj; b; b = b->n) n ++;
    for (b = dia->adjext; b; b = b->n) n ++;
  }

  *num_lists = ldy->ndiab + ldy->nins;
  *num_pins =  n;
  *format = ZOLTAN_COMPRESSED_EDGE;
  *ierr = ZOLTAN_OK;
}

/* graph edges (compressed rows) of W */
static void edge_list (LOCDYN *ldy, int num_gid_entries, int num_vtx_edge, int num_pins,
  int format, ZOLTAN_ID_PTR vtxedge_GID, int *vtxedge_ptr, ZOLTAN_ID_PTR pin_GID, int *ierr)
{
  ZOLTAN_ID_PTR GID = pin_GID;
  DIAB *dia;
  OFFB *b;
  int j, n = 0;

  for (dia = ldy->diab; dia; dia = dia->n, vtxedge_GID ++, vtxedge_ptr ++)
  {
    *vtxedge_GID = dia->id;
    *vtxedge_ptr = n;
    *GID = dia->id; n ++, GID ++;

    for (b = dia->adj; b; b = b->n, n ++, GID ++) *GID = b->id;

    ID_sort (&pin_GID [*vtxedge_ptr], GID-1);
  }

  for (j = 0; j < ldy->nins; j ++, vtxedge_GID ++, vtxedge_ptr ++)
  {
    dia = ldy->ins [j];
    *vtxedge_GID = dia->id;
    *vtxedge_ptr = n;
    *GID = dia->id; n ++, GID ++;

    for (b = dia->adj; b; b = b->n, n ++, GID ++) *GID = b->id;
    for (b = dia->adjext; b; b = b->n, n ++, GID ++) *GID = b->id;

    ID_sort (&pin_GID [*vtxedge_ptr], GID-1);
  }

  *ierr = ZOLTAN_OK;
}

/* pack off-diagonal block ids */
static void pack_block_offids (DIAB *dia, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  OFFB *b;
  int n;

  for (n = 0, b = dia->adj; b; b = b->n) n ++; /* number of internal blocks */
  for (b = dia->adjext; b; b = b->n) n ++; /* number of external blocks */

  pack_int (isize, i, ints, n); /* total number of off-diagonal blocks */

  for (b = dia->adj; b; b = b->n) /* pack internal blocks */
    pack_int (isize, i, ints, b->id);

  for (b = dia->adjext; b; b = b->n) /* pack external blocks */
    pack_int (isize, i, ints, b->id);
}

/* unpack off-diagonal block ids */
static void unpack_block_offids (DIAB *dia, MEM *offmem, MAP *idbb, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  OFFB *b, *n;
  int m;

  /* free current off-diagonal blocks */
  for (b = dia->adj; b; b = n)
  {
    n = b->n;
    MEM_Free (offmem, b);
  }
  dia->adj = NULL; /* clear list */

  m = unpack_int (ipos, i, ints); /* number of off-diagonal blocks */
  for (; m > 0; m --) /* unpack off-diagonal block ids */
  {
    ERRMEM (b = MEM_Alloc (offmem));

    b->id = unpack_int (ipos, i, ints);
    b->n = dia->adj;
    dia->adj = b;
  }

  /* TODO: remove duplicates */
}

/* pack diagonal block */
static void pack_block (DIAB *dia, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  OFFB *b;
  int n;

  pack_doubles (dsize, d, doubles, dia->R, 3);
  pack_doubles (dsize, d, doubles, dia->U, 3);
  pack_doubles (dsize, d, doubles, dia->V, 3);
  pack_doubles (dsize, d, doubles, dia->B, 3);
  pack_doubles (dsize, d, doubles, dia->W, 9);
  pack_double  (dsize, d, doubles, dia->rho);

  pack_doubles (dsize, d, doubles, dia->Z, DOM_Z_SIZE);
  pack_doubles (dsize, d, doubles, dia->point, 3);
  pack_doubles (dsize, d, doubles, dia->base, 9);
  pack_doubles (dsize, d, doubles, dia->mpnt, 3);
  pack_double  (dsize, d, doubles, dia->gap);
  pack_int     (isize, i, ints, dia->kind);
  SURFACE_MATERIAL_Pack_Data (&dia->mat, dsize, d, doubles, isize, i, ints);

  for (n = 0, b = dia->adj; b; b = b->n) n ++; /* number of internal blocks */
  for (b = dia->adjext; b; b = b->n) n ++; /* number of external blocks */

  pack_int (isize, i, ints, n); /* total number of off-diagonal blocks */

  for (b = dia->adj; b; b = b->n) /* pack internal blocks */
  {
    pack_doubles (dsize, d, doubles, b->W, 9);
    pack_int (isize, i, ints, b->id);
  }

  for (b = dia->adjext; b; b = b->n) /* pack external  blocks */
  {
    pack_doubles (dsize, d, doubles, b->W, 9);
    pack_int (isize, i, ints, b->id);
  }
}

/* unpack diagonal block */
static void unpack_block (DIAB *dia, MEM *offmem, MAP *idbb, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  OFFB *b, *n;
  int m;

  unpack_doubles (dpos, d, doubles, dia->R, 3);
  unpack_doubles (dpos, d, doubles, dia->U, 3);
  unpack_doubles (dpos, d, doubles, dia->V, 3);
  unpack_doubles (dpos, d, doubles, dia->B, 3);
  unpack_doubles (dpos, d, doubles, dia->W, 9);
  dia->rho = unpack_double  (dpos, d, doubles);

  unpack_doubles (dpos, d, doubles, dia->Z, DOM_Z_SIZE);
  unpack_doubles (dpos, d, doubles, dia->point, 3);
  unpack_doubles (dpos, d, doubles, dia->base, 9);
  unpack_doubles (dpos, d, doubles, dia->mpnt, 3);
  dia->gap = unpack_double  (dpos, d, doubles);
  dia->kind = unpack_int (ipos, i, ints);
  SURFACE_MATERIAL_Unpack_Data (&dia->mat, dpos, d, doubles, ipos, i, ints);

  /* free current off-diagonal blocks */
  for (b = dia->adj; b; b = n)
  {
    n = b->n;
    MEM_Free (offmem, b);
  }
  dia->adj = NULL; /* clear list */

  m = unpack_int (ipos, i, ints); /* number of off-diagonal blocks */
  for (; m > 0; m --) /* unpack off-diagonal blocks */
  {
    ERRMEM (b = MEM_Alloc (offmem));

    unpack_doubles (dpos, d, doubles, b->W, 9);
    b->id = unpack_int (ipos, i, ints);
    b->dia = MAP_Find (idbb, (void*)b->id, NULL); /* might be NULL for boundary nodes of the local graph portion */
    b->n = dia->adj;
    dia->adj = b;
  }

  /* TODO: sum up and remove duplicates (happen due to multiple contacts between same two bodies)*/
}

/* delete balanced block */
static void delete_balanced_block (LOCDYN *ldy, DIAB *dia)
{
  OFFB *b, *n;

  MAP_Delete (&ldy->mapmem, &ldy->idbb, (void*)dia->id, NULL);

  if (dia->p) dia->p->n = dia->n;
  else ldy->diab = dia->n;
  if (dia->n) dia->n->p = dia->p;

  for (b = dia->adj; b; b = n) /* delete off-diagonal blocks */
  {
    n = b->n;
    MEM_Free (&ldy->offmem, b);
  }

  MEM_Free (&ldy->diamem, dia);

  ldy->ndiab --;
}

/* pack deletion data */
static void pack_delete (SET *del, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  /* ids of blocks to be deleted */
  pack_int (isize, i, ints, SET_Size (del));
  for (SET *item = SET_First (del); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int)item->data);
}

/* unpack deletion data */
static void* unpack_delete  (LOCDYN *ldy, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int ndel;

  /* delete blocks */
  ndel = unpack_int (ipos, i, ints);
  for (; ndel > 0; ndel --)
  {
    DIAB *dia;
    int id;

    id = unpack_int (ipos, i, ints);
    dia = MAP_Find (ldy->idbb, (void*)id, NULL);
    if (dia) delete_balanced_block (ldy, dia); /* could be NULL if deletion followed insertion before balancing */
  }

  return NULL;
}

/* pack migration data */
static void pack_migrate (DIAB *dia, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  pack_int (isize, i, ints, dia->id);
  pack_int (isize, i, ints, dia->rank);
}

/* unpack migration data */
static void* unpack_migrate (LOCDYN *ldy, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  DIAB *dia;

  ERRMEM (dia = MEM_Alloc (&ldy->diamem));
  dia->id = unpack_int (ipos, i, ints);
  dia->rank = unpack_int (ipos, i, ints);
  dia->R = dia->REAC;

  dia->n = ldy->diab;
  if (ldy->diab) ldy->diab->p = dia;
  ldy->diab = dia;

  MAP_Insert (&ldy->mapmem, &ldy->idbb, (void*)dia->id, dia, NULL);

  ldy->ndiab ++;

  return dia;
}

/* pack off-diagonal ids data */
static void pack_offids (SET *upd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  /* pack updated blocks */
  pack_int (isize, i, ints, SET_Size (upd));
  for (SET *item = SET_First (upd); item; item = SET_Next (item))
  {
    DIAB *dia = item->data;

    pack_int (isize, i, ints, dia->id);
    pack_block_offids (dia, dsize, d, doubles, isize, i, ints);
  }
}

/* unpack off-diagonal ids data */
static void* unpack_offids (LOCDYN *ldy, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int nupd;

  nupd = unpack_int (ipos, i, ints);
  for (; nupd > 0; nupd --)
  {
    DIAB *dia;
    int id;

    id = unpack_int (ipos, i, ints);
    ASSERT_DEBUG_EXT (dia = MAP_Find (ldy->idbb, (void*)id, NULL), "Invalid block id");
    unpack_block_offids (dia, &ldy->offmem, ldy->idbb, dpos, d, doubles, ipos, i, ints);
  }

  return NULL;
}

/* pack update data */
static void pack_update (SET *upd, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  /* pack updated blocks */
  pack_int (isize, i, ints, SET_Size (upd));
  for (SET *item = SET_First (upd); item; item = SET_Next (item))
  {
    DIAB *dia = item->data;

    pack_int (isize, i, ints, dia->id);
    pack_block (dia, dsize, d, doubles, isize, i, ints);
  }
}

/* unpack update data */
static void* unpack_update (LOCDYN *ldy, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
{
  int nupd;

  nupd = unpack_int (ipos, i, ints);
  for (; nupd > 0; nupd --)
  {
    DIAB *dia;
    int id;

    id = unpack_int (ipos, i, ints);
    ASSERT_DEBUG_EXT (dia = MAP_Find (ldy->idbb, (void*)id, NULL), "Invalid block id");
    unpack_block (dia, &ldy->offmem, ldy->idbb, dpos, d, doubles, ipos, i, ints);
  }

  return NULL;
}

/* clear external adjacency of a block */
static void clear_adjext (LOCDYN *ldy, DIAB *dia)
{
  OFFB *b, *n;

  for (b = dia->adjext; b; b = n)
  {
    n = b->n;
    MEM_Free (&ldy->offmem, b);
  }

  dia->adjext = NULL;
}

/* build external adjacency */
static void locdyn_adjext (LOCDYN *ldy)
{
  CONEXT *ext;
  SET *item;
  BODY *bod;
  DIAB *dia;
  CON *con;
  OFFB *b;

  for (dia = ldy->dia; dia; dia = dia->n) clear_adjext (ldy, dia);

  for (ext = DOM(ldy->dom)->conext; ext; ext = ext->next) /* for each external constraint in the domain */
  {
    bod = ext->bod; /* for all involved bodies (parents and children) */

    if (bod->kind == OBS) continue; /* obstacles do not trasnder adjacency */

    for (item = SET_First (bod->con); item; item = SET_Next (item))  /* for each regular constraint */
    {
      con = item->data;
      dia = con->dia;

      ERRMEM (b = MEM_Alloc (&ldy->offmem));
      b->dia = NULL; /* there is no diagonal block here */
      b->bod = bod; /* adjacent through this body */
      b->id = ext->id; /* useful when migrated */
      b->ext = ext;
      b->n = dia->adjext;
      dia->adjext = b;
    }
  }
}

/* balance local dynamics */
static int locdyn_balance (LOCDYN *ldy)
{
  MEM setmem;
  DIAB *dia;
  DOM *dom;
  MAP *map, /* maps ranks to sets */
      *set;
  int i, j;

  dom = ldy->dom;
  MEM_Init (&setmem, sizeof (SET), BLKSIZE);
  map = NULL;

  /* map deleted blocks to ranks */
  for (i = 0; i < ldy->ndel; i ++)
  {
    int rank = ldy->del [i]->rank;

    if (!(set = MAP_Find_Node (map, (void*)rank, NULL)))
      set = MAP_Insert (&ldy->mapmem, &map, (void*)rank, NULL, NULL);

    SET_Insert (&setmem, (SET**)&set->data, (void*)ldy->del [i]->id, NULL);
  }

  COMOBJ *send, *recv, *ptr;
  int nsend, nrecv;
  MAP *item;

  nsend = MAP_Size (map);
  ERRMEM (send = malloc (sizeof (COMOBJ [nsend])));

  for (ptr = send, item = MAP_First (map); item; item = MAP_Next (item), ptr ++)
  {
    ptr->rank = (int)item->key;
    ptr->o = item->data;
  }

  /* communicate deleted blocks */
  COMOBJS (MPI_COMM_WORLD, TAG_LOCDYN_DELETE, (OBJ_Pack)pack_delete, ldy, (OBJ_Unpack)unpack_delete, send, nsend, &recv, &nrecv);

  MAP_Free (&ldy->mapmem, &map);
  MEM_Release (&setmem);
  free (send);
  free (recv);

  /* prepare update of off-diagonal block ids of balanced blocks */
  for (map = NULL, dia = ldy->dia; dia; dia = dia->n)
  {
    if (!MAP_Find_Node (ldy->insmap, dia, NULL)) /* not newly inserted */
    {
      if (!(set = MAP_Find_Node (map, (void*)dia->rank, NULL)))
	set = MAP_Insert (&ldy->mapmem, &map, (void*)dia->rank, NULL, NULL);

      SET_Insert (&setmem, (SET**)&set->data, dia, NULL);
    }
  }

  nsend = MAP_Size (map);
  ERRMEM (send = malloc (sizeof (COMOBJ [nsend])));

  for (ptr = send, item = MAP_First (map); item; item = MAP_Next (item), ptr ++)
  {
    ptr->rank = (int)item->key;
    ptr->o = item->data;
  }

  /* communicate updated off-diagonal block ids to balanced blocks */
  COMOBJS (MPI_COMM_WORLD, TAG_LOCDYN_OFFIDS, (OBJ_Pack)pack_offids, ldy, (OBJ_Unpack)unpack_offids, send, nsend, &recv, &nrecv);

  MAP_Free (&ldy->mapmem, &map);
  MEM_Release (&setmem);
  free (send);
  free (recv);

  /* graph balancing data */
  int changes,
      num_gid_entries,
      num_lid_entries,
      num_import,
      *import_procs,
      num_export,
      *export_procs;

  ZOLTAN_ID_PTR import_global_ids,
		import_local_ids,
		export_global_ids,
		export_local_ids;

  /* balance graph comprising the old balanced blocks and the newly inserted unbalanced blocks */
  ASSERT (Zoltan_LB_Balance (ldy->zol, &changes, &num_gid_entries, &num_lid_entries,
	  &num_import, &import_global_ids, &import_local_ids, &import_procs,
	  &num_export, &export_global_ids, &export_local_ids, &export_procs) == ZOLTAN_OK, ERR_ZOLTAN);

  nsend = num_export;
  ERRMEM (send = malloc (sizeof (COMOBJ [nsend])));
  map = NULL;

  for (ptr = send, i = 0; i < num_export; i ++, ptr ++)
  {
    ZOLTAN_ID_TYPE m = export_local_ids [i];

    if (m == UINT_MAX) /* mapped migrated block */
    {
      ASSERT_DEBUG_EXT (dia = MAP_Find (ldy->idbb, (void*)export_global_ids [i], NULL), "Invalid block id");

      /* parent block will need to have updated rank of its migrated copy */
      if (!(set = MAP_Find_Node (map, (void*)dia->rank, NULL)))
	set = MAP_Insert (&ldy->mapmem, &map, (void*)dia->rank, NULL, NULL);

      SET_Insert (&setmem, (SET**)&set->data, ptr, NULL); /* map (newrank, dia) the rank set of its parent */
    }
    else /* use local index */
    {
      ASSERT_DEBUG (m < (unsigned)ldy->nins, "Invalid local index");
      dia = ldy->ins [m];
    }

    ptr->rank = export_procs [i];
    ptr->o = dia;
  }

  /* communicate migration data (unpacking inserts the imported blocks) */
  COMOBJS (MPI_COMM_WORLD, TAG_LOCDYN_BALANCE, (OBJ_Pack)pack_migrate, ldy, (OBJ_Unpack)unpack_migrate, send, nsend, &recv, &nrecv);

  /* communicate rank updates to parents */
  COMDATA *dsend, *drecv, *dtr;
  int dnsend, dnrecv;

  dnsend = MAP_Size (map);
  ERRMEM (dsend = malloc (sizeof (COMDATA [dnsend])));

  for (dtr = dsend, item = MAP_First (map); item; item = MAP_Next (item), dtr ++)
  {
    SET *set = item->data,
	*jtem;

    dtr->rank = (int)item->key;
    dtr->ints = 2 * SET_Size (set);
    ERRMEM (dtr->i = malloc (sizeof (int [dtr->ints])));
    dtr->doubles = 0;

    for (i = 0, jtem = SET_First (set); jtem; i += 2, jtem = SET_Next (jtem))
    {
      ptr = jtem->data;
      dia = ptr->o;
      dtr->i [i] = dia->id;
      dtr->i [i + 1] = ptr->rank;
    }
  }

  /* communicate rank updates */
  COM (MPI_COMM_WORLD, TAG_LOCDYN_RANKS, dsend, dnsend, &drecv, &dnrecv);

  /* update ranks of blocks of received ids */
  for (i = 0, dtr = drecv; i < dnrecv; i ++, dtr ++)
  {
    CON *con;

    for (j = 0; j < dtr->ints; j += 2)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*)dtr->i [j], NULL), "Invalid constraint id");
      con->dia->rank = dtr->i [j + 1]; /* the current rank of the child copy of this block */
    }
  }

  /* update ranks of exported local blocks */
  for (i = 0; i < num_export; i ++)
  {
    ZOLTAN_ID_TYPE m = export_local_ids [i];

    if (m < UINT_MAX) /* local index */
    {
      dia = ldy->ins [m];
      dia->rank = export_procs [i];
    }
  }

  MEM_Release (&setmem);
  MAP_Free (&ldy->mapmem, &map);
  for (i = 0; i < dnsend; i ++) free (dsend [i].i);
  free (dsend);
  free (drecv);

  /* delete migrated balanced blocks */
  for (i = 0, ptr = send; i < nsend; i ++, ptr ++)
  {
    dia = ptr->o;
    if (dia->con == NULL) /* if balanced */
      delete_balanced_block (ldy, dia);
  }

  Zoltan_LB_Free_Data (&import_global_ids, &import_local_ids, &import_procs,
                       &export_global_ids, &export_local_ids, &export_procs);

  free (send);
  free (recv);

  /* internally 'migrate' blocks from the insertion list to the balanced
   * block list (those that have not been exported away from here) */
  int dsize = 0, isize = 0;
  double *dd = NULL;
  int *ii = NULL;
  for (i = 0; i < ldy->nins; i ++)
  {
    DIAB *dia = ldy->ins [i];

    if (dia->rank == dom->rank) /* same rank => not exported */
    {
      int doubles = 0, ints = 0, dpos = 0, ipos = 0;

      /* use migration routines as a wrapper here */
      pack_migrate (dia, &dsize, &dd, &doubles, &isize, &ii, &ints);
      unpack_migrate (ldy, &dpos, dd, doubles, &ipos, ii, ints);
    }
  }
  free (dd);
  free (ii);

  /* update constraint coppied data of unbalanced blocks and create update
   * sets at the same time; update all blocks without regard for the rank
   * (local or exported) => use communication to update all */
  for (map = NULL, dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;

    COPY (dia->R, dia->REAC);
    for (i = 0; i < DOM_Z_SIZE; i ++)
      dia->Z [i] = con->Z [i];
    COPY (con->point, dia->point);
    NNCOPY (con->base, dia->base);
    COPY (con->mpnt, dia->mpnt);
    dia->gap = con->gap;
    dia->kind = con->kind;
    dia->mat = con->mat;

    if (!(set = MAP_Find_Node (map, (void*)dia->rank, NULL)))
      set = MAP_Insert (&ldy->mapmem, &map, (void*)dia->rank, NULL, NULL);

    SET_Insert (&setmem, (SET**)&set->data, dia, NULL);
  }

  nsend = MAP_Size (map);
  ERRMEM (send = malloc (sizeof (COMOBJ [nsend])));

  for (ptr = send, item = MAP_First (map); item; item = MAP_Next (item), ptr ++)
  {
    ptr->rank = (int)item->key;
    ptr->o = item->data;
  }

  /* communicate updated blocks */
  COMOBJS (MPI_COMM_WORLD, TAG_LOCDYN_UPDATE, (OBJ_Pack)pack_update, ldy, (OBJ_Unpack)unpack_update, send, nsend, &recv, &nrecv);

  MAP_Free (&ldy->mapmem, &map);
  MEM_Release (&setmem);
  free (send);
  free (recv);

  /* empty the deleted parents table and free parents */
  for (i = 0; i < ldy->ndel; i ++) MEM_Free (&ldy->diamem, ldy->del [i]);
  ldy->ndel = 0;

  /* empty the inserted parents table */
  MAP_Free (&ldy->mapmem, &ldy->insmap);
  ldy->nins = 0;

  return changes;
}

/* cummunicate reactions from balanced 
 * (migrated) to unbalanced (local) systems */
static void locdyn_gossip (LOCDYN *ldy)
{
  double *R;
  MEM setmem;
  DIAB *dia;
  MAP *map, /* map ranks to sets of balanced blocks */
      *set;
  DOM *dom; 
  CON *con;
  int i, j;

  dom = DOM (ldy->dom);

  MEM_Init (&setmem, sizeof (SET), BLKSIZE);

  for (map = NULL, dia = ldy->diab; dia; dia = dia->n)
  {
    if (dia->rank != dom->rank) /* send */
    {
      if (!(set = MAP_Find_Node (map, (void*)dia->rank, NULL)))
	set = MAP_Insert (&ldy->mapmem, &map, (void*)dia->rank, NULL, NULL);

      SET_Insert (&setmem, (SET**)&set->data, dia, NULL); /* map balanced blocks to their parent ranks */
    }
    else /* do not send */
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*)dia->id, NULL), "Invalid constraint id"); /* find constraint */
      R = con->dia->R;
      COPY (dia->R, R); /* update reaction */
    }
  }

  COMDATA *send, *recv, *ptr;
  int nsend, nrecv;
  SET *item;

  nsend = MAP_Size (map); /* the number of distinct parent ranks */
  ERRMEM (send = malloc (sizeof (COMDATA [nsend])));

  for (set = MAP_First (map), ptr = send; set; set = MAP_Next (set), ptr ++)
  {
    ptr->rank = (int)set->key;
    ptr->ints = SET_Size ((SET*)set->data);
    ptr->doubles = 3 * ptr->ints;
    ERRMEM (ptr->i = malloc (sizeof (int [ptr->ints]))); /* ids */
    ERRMEM (ptr->d = malloc (sizeof (double [ptr->doubles]))); /* reactions */

    for (i = 0, item = SET_First ((SET*)set->data); item; i ++, item = SET_Next (item))
    {
      dia = item->data;

      ptr->i [i] = dia->id;
      COPY (dia->R, &ptr->d [3*i]);
    }
  }

  /* communicate reactions */
  COM (MPI_COMM_WORLD, TAG_LOCDYN_REAC, send, nsend, &recv, &nrecv);

  for (i = 0, ptr = recv; i < nrecv; i ++, ptr ++)
  {
    for (j = 0; j < ptr->ints; j ++)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*)ptr->i[j], NULL), "Invalid constraint id"); /* find constraint */
      dia = con->dia;
      R = &ptr->d [3*j];
      COPY (R, dia->R); /* update reaction */
    }
  }

  MEM_Release (&setmem);
  MAP_Free (&ldy->mapmem, &map);
  for (i = 0, ptr = send; i < nsend; i ++, ptr ++)
  { free (ptr->i); free (ptr->d); }
  free (send);
  free (recv);
}

/* append buffer with a new item */
static void append (DIAB ***buf, int *n, int *s, DIAB *dia)
{
  int i = *n;

  (*n) ++;

  if ((*n) >= (*s))
  {
    (*s) *= 2;
    ERRMEM ((*buf) = realloc (*buf, (*s) * sizeof (DIAB*)));
  }

  (*buf) [i] = dia;
}

/* create MPI related data */
static void create_mpi (LOCDYN *ldy)
{
  ERRMEM (ldy->ins = malloc (BLKSIZE * sizeof (DIAB*)));
  ldy->insmap = NULL;
  ldy->sins = BLKSIZE;
  ldy->nins = 0;

  ERRMEM (ldy->del = malloc (BLKSIZE * sizeof (DIAB*)));
  ldy->sdel = BLKSIZE;
  ldy->ndel = 0;

  MEM_Init (&ldy->mapmem, sizeof (MAP), BLKSIZE);
  ldy->idbb = NULL;
  ldy->diab = NULL;
  ldy->ndiab = 0;

  /* create Zoltan context */
  ASSERT (ldy->zol = Zoltan_Create (MPI_COMM_WORLD), ERR_ZOLTAN);

  /* general parameters */
  Zoltan_Set_Param (ldy->zol, "DEBUG_LEVEL", "0");
  Zoltan_Set_Param (ldy->zol, "DEBUG_MEMORY", "0");
  Zoltan_Set_Param (ldy->zol, "NUM_GID_ENTRIES", "1");
  Zoltan_Set_Param (ldy->zol, "NUM_LID_ENTRIES", "1");

  /* load balaninc parameters */
  Zoltan_Set_Param (ldy->zol, "LB_METHOD", "HYPERGRAPH");
  Zoltan_Set_Param (ldy->zol, "HYPERGRAPH_PACKAGE", "PHG");
  Zoltan_Set_Param (ldy->zol, "AUTO_MIGRATE", "FALSE");
  Zoltan_Set_Param (ldy->zol, "RETURN_LISTS", "EXPORT");

  /* PHG parameters */
  Zoltan_Set_Param (ldy->zol, "PHG_OUTPUT_LEVEL", "0");

  /* callbacks */
  Zoltan_Set_Fn (ldy->zol, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) vertex_count, ldy);
  Zoltan_Set_Fn (ldy->zol, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)()) vertex_list, ldy);
  Zoltan_Set_Fn (ldy->zol, ZOLTAN_HG_SIZE_CS_FN_TYPE, (void (*)()) edge_sizes, ldy);
  Zoltan_Set_Fn (ldy->zol, ZOLTAN_HG_CS_FN_TYPE, (void (*)()) edge_list, ldy);
}

/* destroy MPI related data */
static void destroy_mpi (LOCDYN *ldy)
{
  free (ldy->ins);
  free (ldy->del);
  Zoltan_Destroy (&ldy->zol);
}
#endif

/* create local dynamics for a domain */
LOCDYN* LOCDYN_Create (void *dom)
{
  LOCDYN *ldy;

  ERRMEM (ldy = malloc (sizeof (LOCDYN)));
  MEM_Init (&ldy->offmem, sizeof (OFFB), BLKSIZE);
  MEM_Init (&ldy->diamem, sizeof (DIAB), BLKSIZE);
  ldy->dom = dom;
  ldy->dia = NULL;
  ldy->modified = 0;

#if MPI
  create_mpi (ldy);
#endif

  return ldy;
}

/* insert a 'con'straint between a pair of bodies =>
 * return the diagonal entry of the local dynamical system */
DIAB* LOCDYN_Insert (LOCDYN *ldy, void *con, BODY *one, BODY *two)
{
  DIAB *dia, *nei;
  SET *item;
  OFFB *b;
  CON *c;

  ERRMEM (dia = MEM_Alloc (&ldy->diamem));
  dia->R = CON(con)->R;
  dia->con = con;

#if MPI
  dia->id = CON(con)->id;
  dia->rank = DOM(ldy->dom)->rank;

  /* schedule for balancing use */
  append (&ldy->ins, &ldy->nins, &ldy->sins, dia);
  MAP_Insert (&ldy->mapmem, &ldy->insmap, dia, (void*)(ldy->nins-1), NULL);
#endif

  /* insert into list */
  dia->n = ldy->dia;
  if (ldy->dia)
    ldy->dia->p = dia;
  ldy->dia = dia;

  if (one && one->kind != OBS) /* obstacles do not transfer adjacency */
  {
    for (item = SET_First (one->con); item; item = SET_Next (item))
    {
      c = item->data;

      if (c != con && c->dia) /* skip the coincident or unattached yet constraint */
      {
	nei = c->dia;

        /* allocate block and put into 'nei->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = dia; /* adjacent with 'dia' */
	b->bod = one; /* adjacent trough body 'one' */
	b->n = nei->adj; /* extend list ... */
	nei->adj = b; /* ... */
#if MPI
	b->id = dia->id;
#endif

	/* allocate block and put into 'dia->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = nei; /* adjacent with 'nei' */
	b->bod = one; /* ... trough 'one' */
	b->n = dia->adj;
	dia->adj = b;
#if MPI
	b->id = nei->id;
#endif
      }
    }
  }

  if (two && two->kind != OBS) /* 'one' replaced with 'two' */
  {
    for (item = SET_First (two->con); item; item = SET_Next (item))
    {
      c = item->data;

      if (c != con && c->dia) /* skip the coincident or unattached yet constraint */
      {
	nei = c->dia;

        /* allocate block and put into 'nei->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = dia; /* adjacent with 'dia' */
	b->bod = two; /* adjacent trough body 'two' */
	b->n = nei->adj; /* extend list ... */
	nei->adj = b; /* ... */
#if MPI
	b->id = dia->id;
#endif

	/* allocate block and put into 'dia->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = nei; /* adjacent with 'nei' */
	b->bod = two; /* ... trough 'two' */
	b->n = dia->adj;
	dia->adj = b;
#if MPI
	b->id = nei->id;
#endif
      }
    }
  }

  /* mark as modified */
  ldy->modified = 1;

  return dia;
}

/* remove a diagonal entry from local dynamics */
void LOCDYN_Remove (LOCDYN *ldy, DIAB *dia)
{
  OFFB *b, *c, *r;

  /* destroy blocks in
   * adjacent dia items */
  for (b = dia->adj; b; b = b->n)
  {
    c = b->dia->adj;
    if (c && c->dia == dia)
    {
      b->dia->adj = c->n;
      MEM_Free (&ldy->offmem, c); 
    }
    else for (; c; c = c->n)
    {
      if (c->n &&
	  c->n->dia == dia)
      {
	r = c->n;
        c->n = c->n->n;
        MEM_Free (&ldy->offmem, r); 
	break;
      }
    }
  }

  /* destroy directly
   * adjacent blocks */
  for (b = dia->adj; b; b = c)
  {
    c = b->n;
    MEM_Free (&ldy->offmem, b);
  }

  /* remove from list */
  if (dia->p)
    dia->p->n = dia->n;
  else ldy->dia = dia->n;
  if (dia->n)
    dia->n->p = dia->p;

#if MPI
  MAP *item, *jtem;

  if ((item = MAP_Find_Node (ldy->insmap, dia, NULL))) /* was inserted and now is deleted before balancing */
  {
    ldy->ins [(int)item->data] = ldy->ins [-- ldy->nins]; /* replace this item withe the last one */
    ASSERT_DEBUG (jtem = MAP_Find_Node (ldy->insmap, ldy->ins [(int)item->data], NULL), "Failed to find an inserted block");
    jtem->data = item->data; /* update block to index mapping */
    MAP_Delete_Node (&ldy->mapmem, &ldy->insmap, item); /* remove from map */
    MEM_Free (&ldy->diamem, dia); /* and free */
  }
  else
  {
    clear_adjext (ldy, dia); /* clear external adjacency */
    append (&ldy->del, &ldy->ndel, &ldy->sdel, dia); /* schedule for balancing use and a later deletion */
  }
#else
  /* destroy passed dia */
  MEM_Free (&ldy->diamem, dia);
#endif

  /* mark as modified */
  ldy->modified = 1;
}

/* updiae local dynamics => prepare for a solution */
void LOCDYN_Update_Begin (LOCDYN *ldy, UPKIND upkind)
{
  DOM *dom = ldy->dom;
  double step = dom->step;
  DIAB *dia;

#if MPI
  locdyn_adjext (ldy);
#endif

  /* calculate local velocities and
   * assmeble the force-velocity 'W' operator */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    CON *con = dia->con;
    BODY *m = con->master,
	 *s = con->slave;
    void *mgobj = con->mgobj,
	 *sgobj = con->sgobj;
    SHAPE *mshp = con->mshp,
	  *sshp = con->sshp;
    double *mpnt = con->mpnt,
	   *spnt = con->spnt,
	   *base = con->base,
	   *V = dia->V,
	   *B = dia->B,
           X [3], Y [9];
    MX_DENSE_PTR (W, 3, 3, dia->W);
    MX_DENSE (C, 3, 3);
    MX *mH, *sH;
    OFFB *blk;

    /* previous time step velocity */
    BODY_Local_Velo (m, PREVELO, mshp, mgobj, mpnt, base, X); /* master body pointer cannot be NULL */
    if (s) BODY_Local_Velo (s, PREVELO, sshp, sgobj, spnt, base, Y); /* might be NULL for some constraints (one body) */
    else { SET (Y, 0.0); }
    SUB (Y, X, V); /* relative = slave - master => outward master normal */

    /* local free velocity */
    BODY_Local_Velo (m, CURVELO, mshp, mgobj, mpnt, base, X);
    if (s) BODY_Local_Velo (s, CURVELO, sshp, sgobj, spnt, base, Y);
    else { SET (Y, 0.0); }
    SUB (Y, X, B);

    /* diagonal block */
    mH = BODY_Gen_To_Loc_Operator (m, mshp, mgobj, mpnt, base);
    MX_Trimat (mH, m->inverse, MX_Tran (mH), &W); /* H * inv (M) * H^T */
    if (s)
    { sH = BODY_Gen_To_Loc_Operator (s, sshp, sgobj, spnt, base);
      MX_Trimat (sH, s->inverse, MX_Tran (sH), &C); /* H * inv (M) * H^T */
      NNADD (W.x, C.x, W.x); }
    SCALE9 (W.x, step); /* W = h * ( ... ) */

    NNCOPY (W.x, C.x); /* calculate regularisation parameter */
    ASSERT (lapack_dsyev ('N', 'U', 3, C.x, 3, X, Y, 9) == 0, ERR_LDY_EIGEN_DECOMP);
    dia->rho = 1.0 / X [2]; /* inverse of maximal eigenvalue */

    /* off-diagonal blocks if requested */
    for (blk = dia->adj; upkind == UPALL && blk; blk = blk->n)
    {
      DIAB *dia = blk->dia;
      CON *con = dia->con;
      BODY *bod = blk->bod;
      MX *lH, *rH, *inv;
      MX_DENSE_PTR (W, 3, 3, blk->W);
      double coef;

      ASSERT_DEBUG (bod == m || bod == s, "Off diagonal block is not connected!");
     
      lH = (bod == m ? mH : sH); /* dia->bod is a valid body (not an obstacle)
                                   as it was inserted into the dual graph */
      inv = bod->inverse;

      if (bod == con->master)
      {
	rH =  BODY_Gen_To_Loc_Operator (bod, con->mshp, con->mgobj, con->mpnt, con->base);
	coef = (bod == s ? -step : step);
      }
      else /* blk->bod == dia->slave */
      {
	rH =  BODY_Gen_To_Loc_Operator (bod, con->sshp, con->sgobj, con->spnt, con->base);
	coef = (bod == m ? -step : step);
      }

      MX_Trimat (lH, inv, MX_Tran (rH), &W);
      SCALE9 (W.x, coef);
      MX_Destroy (rH);
    }

#if MPI
    /* off-diagonal external blocks if requested */
    for (blk = dia->adjext; upkind == UPALL && blk; blk = blk->n)
    {
      CONEXT *ext = blk->ext;
      BODY *bod = blk->bod;
      MX *lH, *rH, *inv;
      MX_DENSE_PTR (W, 3, 3, blk->W);
      double coef;

      ASSERT_DEBUG (bod == m || bod == s, "Off diagonal block is not connected!");
     
      lH = (bod == m ? mH : sH);
                               
      inv = bod->inverse;

      rH =  BODY_Gen_To_Loc_Operator (bod, ext->sgp->shp, ext->sgp->gobj, ext->point, ext->base);

      if (ext->isma) coef = (bod == s ? -step : step);
      else coef = (bod == m ? -step : step);

      MX_Trimat (lH, inv, MX_Tran (rH), &W);
      SCALE9 (W.x, coef);
      MX_Destroy (rH);
    }
#endif

    MX_Destroy (mH);
    if (s) MX_Destroy (sH);
  }

  /* forward variables change */
  variables_change_begin (ldy);

#if MPI
  locdyn_balance (ldy);
#endif
}

/* updiae local dynamics => after the solution */
void LOCDYN_Update_End (LOCDYN *ldy)
{
  /* backward variables change */
  variables_change_end (ldy);

  /* not modified */
  ldy->modified = 0;

#if MPI
  locdyn_gossip (ldy);
#endif
}

/* free memory */
void LOCDYN_Destroy (LOCDYN *ldy)
{
  MEM_Release (&ldy->diamem);
  MEM_Release (&ldy->offmem);

#if MPI
  destroy_mpi (ldy);
#endif

  free (ldy);
}
