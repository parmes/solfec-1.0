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
    local_ids [i * num_lid_entries] = -1; /* mapped */
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
    for (b = dia->adjext; b; b = b->n) n ++;
    for (b = dia->adj; b; b = b->n) n ++;
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
  DIAB *dia;
  OFFB *b;
  int j, n = 0;

  for (dia = ldy->diab; dia; dia = dia->n, vtxedge_GID ++, vtxedge_ptr ++)
  {
    *vtxedge_GID = dia->id;
    *vtxedge_ptr = n;
    *pin_GID = dia->id; pin_GID ++;

    for (b = dia->adj; b; b = b->n, n ++, pin_GID ++) *pin_GID = b->id;
  }

  for (j = 0; j < ldy->nins; j ++, vtxedge_GID ++, vtxedge_ptr ++)
  {
    dia = ldy->ins [j];
    *vtxedge_GID = dia->id;
    *vtxedge_ptr = n;
    *pin_GID = dia->id; pin_GID ++;

    for (b = dia->adjext; b; b = b->n, n ++, pin_GID ++) *pin_GID = b->id;
    for (b = dia->adj; b; b = b->n, n ++, pin_GID ++) *pin_GID = b->id;
  }

  *ierr = ZOLTAN_OK;
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

  for (n = 0, b = dia->adjext; b; b = b->n) n ++; /* number of external blocks */
  for (b = dia->adj; b; b = b->n) n ++; /* number of internal blocks */

  pack_int (isize, i, ints, n); /* total number of off-diagonal blocks */

  for (b = dia->adjext; b; b = b->n) /* pack external blocks */
  {
    pack_doubles (dsize, d, doubles, b->W, 9);
    pack_int (isize, i, ints, b->id);
  }

  for (b = dia->adj; b; b = b->n) /* pack internal blocks */
  {
    pack_doubles (dsize, d, doubles, b->W, 9);
    pack_int (isize, i, ints, b->id);
  }
}

/* unpack diagonal block */
static void unpack_block (DIAB *dia, MEM *offmem, int *dpos, double *d, int doubles, int *ipos, int *i, int ints)
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
    b->n = dia->adj;
    dia->adj = b;
  }
}

/* pack deletion data */
static void pack_delete (SET *del, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints)
{
  /* ids of blocks to be deleted */
  pack_int (isize, i, ints, SET_Size (del));
  for (SET *item = SET_First (del); item; item = SET_Next (item))
    pack_int (isize, i, ints, (int)item->data);
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
    unpack_block (dia, &ldy->offmem, dpos, d, doubles, ipos, i, ints);
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

  dia->n = ldy->diab;
  if (ldy->diab) ldy->diab->p = dia;
  ldy->diab = dia;

  MAP_Insert (&ldy->mapmem, &ldy->idbb, (void*)dia->id, dia, NULL);

  ldy->ndiab ++;

  return dia;
}

/* balance local dynamics */
static int balance (LOCDYN *ldy)
{
  MEM setmem;
  DIAB *dia;
  DOM *dom;
  MAP *map, /* maps ranks to sets */
      *set;
  int i, j;

  MEM_Init (&setmem, sizeof (SET), BLKSIZE);
  map = NULL;

  /* map deleted blocks to ranks */
  for (i = 0; i < ldy->ndel; i ++)
  {
    int rank = ldy->del [i]->rank;

    if (!(set = MAP_Find_Node (map, (void*)rank, NULL)))
      set = MAP_Insert (&ldy->mapmem, &map, (void*)rank, NULL, NULL);

    SET_Insert (&setmem, (SET**)&set, (void*)ldy->del [i]->id, NULL);
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
    int m = export_local_ids [i];

    if (m < 0) /* mapped migrated block */
    {
      ASSERT_DEBUG_EXT (dia = MAP_Find (ldy->idbb, (void*)export_global_ids [i], NULL), "Invalid block id");

      /* parent block will need to have updated rank of its migrated copy */
      if (!(set = MAP_Find_Node (map, (void*)dia->rank, NULL)))
	set = MAP_Insert (&ldy->mapmem, &map, (void*)dia->rank, NULL, NULL);
      SET_Insert (&setmem, (SET**)&set, (void*)dia->id, NULL); /* map the id of this block to the rank set of its parent */
    }
    else /* use local index */
    {
      ASSERT_DEBUG (m < ldy->nins, "Invalid local index");
      dia = ldy->ins [m];
    }

    ptr->rank = export_procs [i];
    ptr->o = dia;
  }

  /* communicate migration data */
  COMOBJS (MPI_COMM_WORLD, TAG_LOCDYN_BALANCE, (OBJ_Pack)pack_migrate, ldy, (OBJ_Unpack)unpack_migrate, send, nsend, &recv, &nrecv);

  /* communicate rank updates to parents */
  COMDATA *dsend, *drecv, *dtr;
  int dnsend, dnrecv;

  dnsend = MAP_Size (map);
  ERRMEM (dsend = malloc (sizeof (COMOBJ [dnsend])));

  for (dtr = dsend, item = MAP_First (map); item; item = MAP_Next (item), dtr ++)
  {
    SET *set = item->data,
	*jtem;

    dtr->rank = (int)item->key;
    dtr->ints = SET_Size (set);
    ERRMEM (dtr->i = malloc (sizeof (int [dtr->ints])));

    for (i = 0, jtem = SET_First (set); jtem; i ++, jtem = SET_Next (jtem))
      dtr->i [i] = (int)jtem->data;
  }

  /* communicate rank updates */
  COM (MPI_COMM_WORLD, TAG_LOCDYN_RANKS, dsend, dnsend, &drecv, &dnrecv);

  /* update ranks */
  dom = ldy->dom;
  for (i = 0, dtr = drecv; i < dnrecv; i ++, dtr ++)
  {
    CON *con;

    for (j = 0; i < dtr->ints; j ++)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*)dtr->i [j], NULL), "Invalid constraint id");
      con->dia->rank = dtr->rank; /* the incoming rank is the current rank of the child copy of this block */
    }
  }

  MAP_Free (&ldy->mapmem, &map);
  MEM_Release (&setmem);
  for (i = 0; i < dnsend; i ++) free (dsend [i].i);
  free (dsend);
  free (drecv);

  /* delete migrated blocks */
  for (i = 0, ptr = send; i < nsend; i ++, ptr ++) delete_balanced_block (ldy, ptr->o);

  Zoltan_LB_Free_Data (&import_global_ids, &import_local_ids, &import_procs,
                       &export_global_ids, &export_local_ids, &export_procs);

  free (send);
  free (recv);

  /* update constraint coppied data of unbalanced blocks 
   * and create update sets at the same time */
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
      MAP_Insert (&ldy->mapmem, &map, (void*)dia->rank, NULL, NULL);

    SET_Insert (&setmem, (SET**)&set, dia, NULL);
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
  ldy->nins = 0;

  return changes;
}

/* cummunicate reactions from balanced 
 * (migrated) to unbalanced (local) systems */
static void gossip (LOCDYN *ldy)
{
  /* TODO */
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

/* create MPI context */
static void create_mpi (LOCDYN *ldy)
{
  ERRMEM (ldy->ins = malloc (BLKSIZE * sizeof (DIAB*)));
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

/* destroy MPI context */
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

  /* insert into list */
  dia->n = ldy->dia;
  if (ldy->dia)
    ldy->dia->p = dia;
  ldy->dia = dia;

  if (one && one->kind != OBS) /* obstacles do not transfer adjacency */
  {
    for (item = SET_First (one->con); item; item = SET_Next (item))
    {
      if (item->data != con) /* skip the coincident constraint */
      {
	c = item->data;
	nei = c->dia;

        /* allocate block and put into 'nei->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = dia; /* adjacent with 'dia' */
	b->bod = one; /* adjacent trough body 'one' */
	b->n = nei->adj; /* extend list ... */
	nei->adj = b; /* ... */

	/* allocate block and put into 'dia->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = nei; /* adjacent with 'nei' */
	b->bod = one; /* ... trough 'one' */
	b->n = dia->adj;
	dia->adj = b;
      }
    }
  }

  if (two && two->kind != OBS) /* 'one' replaced with 'two' */
  {
    for (item = SET_First (two->con); item; item = SET_Next (item))
    {
      if (item->data != con) /* skip the coincident constraint */
      {
	c = item->data;
	nei = c->dia;

        /* allocate block and put into 'nei->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = dia; /* adjacent with 'dia' */
	b->bod = two; /* adjacent trough body 'two' */
	b->n = nei->adj; /* extend list ... */
	nei->adj = b; /* ... */

	/* allocate block and put into 'dia->adj' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = nei; /* adjacent with 'nei' */
	b->bod = two; /* ... trough 'two' */
	b->n = dia->adj;
	dia->adj = b;
      }
    }
  }

  /* mark as modified */
  ldy->modified = 1;

#if MPI
  dia->id = CON(con)->id;
  dia->rank = DOM(ldy->dom)->rank;

  /* schedule for balancing use */
  append (&ldy->ins, &ldy->nins, &ldy->sins, dia);
#endif

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
  /* schedule for balancing use and a later deletion */
  append (&ldy->del, &ldy->ndel, &ldy->sdel, dia);
#else
  /* destroy passed dia */
  MEM_Free (&ldy->diamem, dia);
#endif

  /* mark as modified */
  ldy->modified = 1;
}

#if MPI
/* insert an external constraint */
void LOCDYN_Insert_Ext (LOCDYN *ldy, void *con)
{
  BODY *bod [2];
  CONEXT *ext;
  SET *item;
  DIAB *dia;
  OFFB *b;
  CON *c;
  int i;

  ext = con;
  bod [0] = ext->master;
  bod [1] = ext->slave;

  for (i = 0; i < 2; i ++)
  {
    if (bod [i] && bod [i]->kind != OBS) /* obstacles do not transfer adjacency */
    {
      for (item = SET_First (bod [i]->con); item; item = SET_Next (item))
      {
	c = item->data;
	dia = c->dia;

	/* allocate block and put into 'dia->adjext' list */ 
	ERRMEM (b = MEM_Alloc (&ldy->offmem));
	b->dia = NULL; /* there is no diagonal block here */
	b->ext = ext; /* but there is this external constraint instead */
	b->bod = bod [i]; /* adjacent through this body */
	b->id = ext->id; /* useful when migrated */
	b->n = dia->adjext;
	dia->adjext = b;
      }
    }
  }
}

/* remove all external constraints */
void LOCDYN_Remove_Ext_All (LOCDYN *ldy)
{
  OFFB *b, *n;
  DIAB *dia;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    for (b = dia->adjext; b; b = n)
    {
      n = b->n;
      MEM_Free (&ldy->offmem, b); /* free external off-diagonal blocks */
    }

    dia->adjext = NULL;
  }
}
#endif

/* updiae local dynamics => prepare for a solution */
void LOCDYN_Update_Begin (LOCDYN *ldy, UPKIND upkind)
{
  DOM *dom = ldy->dom;
  double step = dom->step;
  DIAB *dia;

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
    /* external off-diagonal blocks if requested */
    for (blk = dia->adjext; upkind == UPALL && blk; blk = blk->n)
    {
      CONEXT *ext = blk->ext;
      BODY *bod = blk->bod;
      MX *lH, *rH, *inv;
      MX_DENSE_PTR (W, 3, 3, blk->W);
      double coef;

      ASSERT_DEBUG (bod == m || bod == s, "Off diagonal block is not connected!");
     
      lH = (bod == m ? mH : sH); /* dia->bod is a valid body (not an obstacle)
                                   as it was inserted into the dual graph */
      inv = bod->inverse;

      if (bod == ext->master)
      {
	rH =  BODY_Gen_To_Loc_Operator (bod, ext->msgp->shp, ext->msgp->gobj, ext->mpnt, ext->base);
	coef = (bod == s ? -step : step);
      }
      else /* blk->bod == ext->slave */
      {
	rH =  BODY_Gen_To_Loc_Operator (bod, ext->ssgp->shp, ext->ssgp->gobj, ext->spnt, ext->base);
	coef = (bod == m ? -step : step);
      }

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
  balance (ldy);
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
  gossip (ldy);
#endif
}

/* free memory */
void LOCDYN_Destroy (LOCDYN *ldy)
{
#if MPI
  destroy_mpi (ldy);
#endif

  MEM_Release (&ldy->diamem);
  MEM_Release (&ldy->offmem);

  free (ldy);
}
