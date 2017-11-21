/*
 * bgs.c
 * Copyright (C) 2007, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * block gauss seidel solver
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

#include <stdlib.h>
#include <stdio.h>
#include "alg.h"
#include "dom.h"
#include "lap.h"
#include "bgs.h"
#include "pes.h"
#include "err.h"
#include "mrf.h"
#include "sol.h"

#if MPI
#include "tag.h"
#include "com.h"
#include "lis.h"
#endif

/* timers */
#if TIMERS
#define S(LABEL) SOLFEC_Timer_Start (ldy->dom->solfec, LABEL)
#define E(LABEL) SOLFEC_Timer_End (ldy->dom->solfec, LABEL)
#else
#define S(LABEL)
#define E(LABEL)
#endif

#if MPI
/* create rank coloring using adjacency graph between processors derived from the W graph */
static int* processor_coloring (GAUSS_SEIDEL *gs, LOCDYN *ldy)
{
  int i, n, m, ncpu, rank, *color, *size, *disp, *adj;
  SET *adjcpu, *item;
  MEM setmem;
  DIAB *dia;
  OFFB *blk;
  CON *con;

  adjcpu = NULL;
  rank = ldy->dom->rank;
  ncpu = ldy->dom->ncpu;
  MEM_Init (&setmem, sizeof (SET), 128);
  ERRMEM (color = MEM_CALLOC (ncpu * sizeof (int)));
  ERRMEM (disp = malloc (sizeof (int [ncpu + 1])));
  ERRMEM (size = malloc (sizeof (int [ncpu])));

  /* collaps W adjacency into processor adjacency */
  for (dia = ldy->dia; dia; dia = dia->n)
  {
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = (CON*) blk->dia;
      SET_Insert (&setmem, &adjcpu, (void*) (long) con->rank, NULL);
    }
  }

  n = SET_Size (adjcpu);
  MPI_Allgather (&n, 1, MPI_INT, size, 1, MPI_INT, MPI_COMM_WORLD);

  for (i = disp [0] = 0; i < ncpu - 1; i ++) disp [i+1] = disp [i] + size [i];
  for (i = 0, item = SET_First (adjcpu); item; i ++, item = SET_Next (item)) color [i] = (int) (long) item->data;

  m = disp [ncpu] = (disp [ncpu-1] + size [ncpu-1]);
  ERRMEM (adj = malloc (sizeof (int [m])));

  MPI_Allgatherv (color, n, MPI_INT, adj, size, disp, MPI_INT, MPI_COMM_WORLD); /* gather graph adjacency */

  for (i = 0; i < ncpu; i ++) color [i] = 0; /* zero colors */

  for (i = 0; i < ncpu; i ++) /* simple BFS coloring */
  {
    int *j, *k;

    do
    {
      color [i] ++; /* start from first color */

      for (j = &adj[disp[i]], k = &adj[disp[i+1]]; j < k; j ++) /* for each adjacent vertex */
      {
	if (color [*j] == color [i]) break; /* see whether the trial color exists in the adjacency */
      }
    }
    while (j < k); /* if so try next color */
  }

  for (m = i = 0; i < ncpu; i ++) m = MAX (m, color [i]); /* compute number of colors */
  gs->colors = m; /* record number of colors */

  if (rank == 0 && ldy->dom->verbose && gs->verbose)
  {
#if DEBUG
  for (i = 0; i < ncpu; i ++)
  {
    int *j, *k;

    printf ("GAUSS_SEIDEL: RANK %d [%d] ADJCPU:", i, color [i]);
    for (j = &adj[disp[i]], k = &adj[disp[i+1]]; j < k; j ++) printf (" %d [%d]", *j, color [*j]);
    printf ("\n");
  }
#endif
    printf ("GAUSS_SEIDEL: PROCESSOR COLORS = %d\n", m);
  }

  MEM_Release (&setmem);
  free (size);
  free (disp);
  free (adj);

  return color;
}

/* return next pointer and realloc send memory if needed */
inline static COMDATA* sendnext (int nsend, int *size, COMDATA **send)
{
  if (nsend >= *size)
  {
    (*size) *= 2;
    ERRMEM (*send = realloc (*send, sizeof (COMDATA [*size])));
  }

  return &(*send)[nsend];
}

/* receive external reactions */
static void receive_reactions (DOM *dom, COMDATA *recv, int nrecv)
{
  COMDATA *ptr;
  int i, j, *k;
  double *R;
  CON *con;

  for (i = 0, ptr = recv; i < nrecv; i ++, ptr ++)
  {
    for (j = 0, k = ptr->i, R = ptr->d; j < ptr->ints; j ++, k ++, R += 3)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*k), NULL), "Invalid constraint id");
      COPY (R, con->R);
      con->state |= CON_DONE;
    }
  }
}

/* receive reactions updated by middle nodes */
static void receive_middle_reactions (DOM *dom, COMDATA *recv, int nrecv, MEM *setmem, SET **midupd)
{
  COMDATA *ptr;
  int i, j, *k;
  CON *con;

  for (i = 0, ptr = recv; i < nrecv; i ++, ptr ++)
  {
    for (j = 0, k = ptr->i; j < ptr->ints; j ++, k ++)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->conext, (void*) (long) (*k), NULL), "Invalid constraint id");
      SET_Insert (setmem, midupd, con, NULL);
    }
  }
}

/* a single row Gauss-Seidel step */
static int gauss_seidel (GAUSS_SEIDEL *gs, short dynamic, double step, DIAB *dia, double *errup, double *errlo)
{
  double R0 [3], B [3], *R, *W;
  int diagiters;
  OFFB *blk;
  CON *con;

  /* compute local velocity */
  COPY (dia->B, B);
  for (blk = dia->adj; blk; blk = blk->n)
  {
    W = blk->W;
    R = blk->dia->R;
    NVADDMUL (B, W, R, B);
  }
  for (blk = dia->adjext; blk; blk = blk->n)
  {
    con = (CON*) blk->dia;
    W = blk->W;
    R = con->R;
    NVADDMUL (B, W, R, B);
  }
 
  R = dia->R; 
  COPY (R, R0); /* previous reaction */

  /* solve local diagonal block problem */
  con = dia->con;
  diagiters = DIAGONAL_BLOCK_Solver (gs->diagsolver, gs->diagepsilon, gs->diagmaxiter, dynamic,
                    step, con->kind, &con->mat, con->gap, con->area, con->Z, con->base, dia, B);

  if (diagiters >= gs->diagmaxiter || diagiters < 0) /* failed */
  {
    if (con->kind == CONTACT)
    {
      DIAS dias [4] = {DS_SEMISMOOTH_NEWTON, DS_PROJECTED_GRADIENT, DS_DE_SAXCE_FENG, DS_PROJECTED_NEWTON};

      for (int i = 0; i < 4; i ++)
      {
	if (dias [i] != gs->diagsolver) /* skip current diagonal solver */
	{
	  COPY (R0, R); /* initialize with previous reaction */

	  diagiters = DIAGONAL_BLOCK_Solver (dias [i], gs->diagepsilon, gs->diagmaxiter, /* try another solver */
	    dynamic, step, con->kind, &con->mat, con->gap, con->area, con->Z, con->base, dia, B);

	  if (diagiters < gs->diagmaxiter && diagiters >= 0) break; /* success */
	}
      }
    }

    if (diagiters >= gs->diagmaxiter || diagiters < 0) /* failed */
    {
      COPY (R0, R); /* use previous reaction */
    }
  }

  /* accumulate relative
   * error components */
  SUB (R, R0, R0);
  *errup += DOT (R0, R0);
  *errlo += DOT (R, R);

  return diagiters;
}

/* a Guss-Seidel sweep over a set of blocks */
static int gauss_seidel_sweep (SET *set, int reverse, GAUSS_SEIDEL *gs, short dynamic, double step, int loops, double *errup, double *errlo)
{
  SET* (*first) (SET*);
  SET* (*next) (SET*);
  int di, dimax, n;
  double up, lo;

  if (reverse) first = SET_Last, next = SET_Prev;
  else first = SET_First, next = SET_Next;

  dimax = 0;

  for (SET *item = first (set); item; item = next (item)) /* first loop contributes to the outputed error components */
  {
    di = gauss_seidel (gs, dynamic, step, item->data, errup, errlo);
    dimax = MAX (dimax, di);
  }

  for (n = 0, up = lo = 0.0; n < loops-1; n ++)  /* remaining inner loops do not contribute to the outputed error components */
  {
    for (SET *item = first (set); item; item = next (item))
    {
      di = gauss_seidel (gs, dynamic, step, item->data, &up, &lo);
      dimax = MAX (dimax, di);
    }
  }

  return dimax;
}

/* middle node list needs score-based sorting */
typedef struct middle_list MIDDLE_NODE;

struct middle_list
{
  DIAB *dia;
  int score;
  MIDDLE_NODE *next;
};

/* middle node list sorting */
#define MLLE(i, j) ((i)->score <= (j)->score)
IMPLEMENT_LIST_SORT (SINGLE_LINKED, middle_list_sort, MIDDLE_NODE, prev, next, MLLE)

/* perform a Guss-Seidel loop over a set of blocks */
static int gauss_seidel_loop (SET *middle, SET *midupd, int reverse, MEM *setmem, int mycolor, int *color,
                   GAUSS_SEIDEL *gs, LOCDYN *ldy, short dynamic, double step, double *errup, double *errlo)
{
  SET *requs, *ranks, *item, *jtem;
  MIDDLE_NODE *list, *cur;
  int di, dimax;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  MEM lstmem, reqmem;

  MEM_Init (&lstmem, sizeof (MIDDLE_NODE), 128);
  MEM_Init (&reqmem, sizeof (MPI_Request), 128);

  dimax = 0;
  list = NULL;
  requs = NULL;

  /* post receives first */
  for (item = SET_First (midupd); item; item = SET_Next (item))
  {
    MPI_Request *req;
    con = item->data;
    ASSERT_DEBUG ((con->state & CON_DONE) == 0 && con->dia == NULL, "Invalid external constraint");
    ERRMEM (req = MEM_Alloc (&reqmem));
    MPI_Irecv (con->R, 3, MPI_DOUBLE, con->rank, TAG_LAST+con->id, MPI_COMM_WORLD, req);
    con->dia = (DIAB*) req; /* use spare (NULL) DIAB pointer for the request */
  }

  /* create middle node list */
  for (item = SET_First (middle); item; item = SET_Next (item))
  {
    ERRMEM (cur = MEM_Alloc (&lstmem));
    cur->dia = item->data;
    cur->score = 0;
    for (blk = cur->dia->adjext; blk; blk = blk->n)
    {
      con = (CON*) blk->dia;
      if (reverse) cur->score = MIN (cur->score, color [con->rank]); /* smallest color first in sorted list */
      else cur->score = MIN (cur->score, -color [con->rank]); /* largest color first in sorted list */
    }
    cur->next = list;
    list = cur;
  }

  /* sort middle node list */
  list = middle_list_sort (list);

  /* process middle nodes */
  for (cur = list; cur; cur = cur->next)
  {
    dia = cur->dia;

    for (ranks = NULL, blk = dia->adjext; blk; blk = blk->n)
    {
      con = (CON*) blk->dia;
      if ((reverse && mycolor > color [con->rank] && (con->state & CON_DONE) == 0) || /* if reversed iterations receive from lower colors */
          (reverse == 0 && mycolor < color [con->rank] && (con->state & CON_DONE) == 0)) /* else receive from higher colors */
      {
	MPI_Status sta;
	S("GSMCOM"); MPI_Wait ((MPI_Request*)con->dia, &sta); E("GSMCOM");
	con->state |= CON_DONE;
      }

      SET_Insert (setmem, &ranks, (void*) (long) con->rank, NULL); /* schedule for sending to this rank after the reaction is coputed */
    }

    S("GSRUN"); di = gauss_seidel (gs, dynamic, step, dia, errup, errlo); E("GSRUN"); /* compute reaction */
    dimax = MAX (dimax, di);

    con = dia->con;
    for (jtem = SET_First (ranks); jtem; jtem = SET_Next (jtem)) /* update remote external reactions */
    {
      MPI_Request *req;
      ERRMEM (req = MEM_Alloc (&reqmem));
      MPI_Isend (con->R, 3, MPI_DOUBLE, (int) (long) jtem->data, TAG_LAST+con->id, MPI_COMM_WORLD, req); /* send to remote ranks */
      SET_Insert (setmem, &requs, req, NULL);
    }

    SET_Free (setmem, &ranks);
  }

  for (item = SET_First (requs); item; item = SET_Next (item))
  {
    MPI_Status sta;
    S("GSMCOM"); MPI_Wait (item->data, &sta); E("GSMCOM"); /* wait until all sends complete */
  }

  /* process set of blocks updated by middle nodes and look for undone external reactions */
  for (item = SET_First (midupd); item; item = SET_Next (item))
  {
    con = item->data;
    if ((con->state & CON_DONE) == 0) /* undone external reaction found */
    {
      MPI_Status sta;
      S("GSMCOM"); MPI_Wait ((MPI_Request*)con->dia, &sta); E("GSMCOM"); /* receive update */
      con->state |= CON_DONE; /* mark done */
    }
    con->dia = NULL; /* release request pointer */
  }

  MEM_Release (&lstmem);
  MEM_Release (&reqmem);

  return dimax;
}

/* undo all external reactions */
static void undo_all (LOCDYN *ldy)
{
  DIAB *dia;
  OFFB *blk;
  CON *con;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = (CON*) blk->dia;
      con->state &= ~CON_DONE;
    }
  }
}

#if DEBUG
/* test whether all external reactions are done */
static int all_done  (LOCDYN *ldy)
{
  DIAB *dia;
  OFFB *blk;
  CON *con;

  for (dia = ldy->dia; dia; dia = dia->n)
  {
    for (blk = dia->adjext; blk; blk = blk->n)
    {
      con = (CON*) blk->dia;
      if ((con->state & CON_DONE) == 0)
      {
	return 0;
      }
    }
  }

  return 1;
}
#endif
#endif /* MPI */

/* create solver */
GAUSS_SEIDEL* GAUSS_SEIDEL_Create (double epsilon, int maxiter, double meritval, GSFAIL failure,
                                   double diagepsilon, int diagmaxiter, DIAS diagsolver,
				   void *data, GAUSS_SEIDEL_Callback callback)
{
  GAUSS_SEIDEL *gs;

  ERRMEM (gs = malloc (sizeof (GAUSS_SEIDEL)));
  gs->epsilon = epsilon;
  gs->maxiter = maxiter;
  gs->meritval = meritval;
  gs->failure = failure;
  gs->diagepsilon = diagepsilon;
  gs->diagmaxiter = diagmaxiter;
  gs->diagsolver = diagsolver;
  gs->data = data;
  gs->callback = callback;
  gs->rerhist = NULL;
  gs->merhist = NULL;
  gs->reverse = GS_OFF;
  gs->error = GS_OK;
  gs->variant = GS_FULL;
  gs->innerloops = 1;
  gs->verbose = 1;
  gs->nomerit = 0;
  gs->itershist = NULL;
  gs->itershistcount = -1;
  gs->itershistsize = 0;

  return gs;
}

#if MPI
/* run parallel solver */
void GAUSS_SEIDEL_Solve (GAUSS_SEIDEL *gs, LOCDYN *ldy)
{
  int div = 10, di, dimax, diagiters, mycolor, rank, *color;
  short dynamic, verbose, nomerit;
  double error, *merit, step;
  char fmt [512];
  SET *ranks;
  MEM setmem;
  DIAB *dia;
  OFFB *blk;
  CON *con;
  DOM *dom;

  SET *bottom    = NULL,
      *top       = NULL,
      *middle    = NULL,
      *internal  = NULL,
      *midupd    = NULL,
      *int1      = NULL,
      *int2      = NULL,
      *all       = NULL;

  int size1 = 0,
      size2 = 0,
      size3 = 0,
      size4,
      size5;

  void *bot_pattern, /* communication pattern when sending from lower to higher processors */
       *top_pattern, /* the reverse communication pattern */
       *mid_pattern; /* used in case of GS_MIDDLE_JACOBI variant */

  COMDATA *send_bot, *recv_bot, *ptr_bot,
	  *send_top, *recv_top, *ptr_top,
	  *send_mid, *recv_mid, *ptr_mid;

  int size_bot, nsend_bot, nrecv_bot,
      size_top, nsend_top, nrecv_top,
      size_mid, nsend_mid, nrecv_mid;

  COMDATA *send, *recv, *ptr;

  int nsend, nrecv, size;

  S("GSINIT");

  dom = ldy->dom;
  rank = dom->rank;
  merit = &dom->merit;
  verbose = dom->verbose && gs->verbose;

  nomerit = gs->nomerit ? 1 : gs->meritval >= 1.0 ? 1 : 0;

  if (nomerit) *merit = 0.0;

  if (rank == 0 && verbose) sprintf (fmt, "GAUSS_SEIDEL: iteration: %%%dd  error:  %%.2e  merit:  %%.2e\n", (int)log10 (gs->maxiter) + 1);

  gs->rerhist = realloc (gs->rerhist, gs->maxiter * sizeof (double));
  gs->merhist = realloc (gs->merhist, gs->maxiter * sizeof (double));

  MEM_Init (&setmem, sizeof (SET), 256);

  if (gs->variant < GS_BOUNDARY_JACOBI)
  {
    color = processor_coloring (gs, ldy); /* color processors */
    mycolor = color [rank];

    /* create block sets */
    for (dia = ldy->dia; dia; dia = dia->n)
    {
      int lo = 0, hi = 0;

      for (blk = dia->adjext; blk; blk = blk->n)
      {
	con = (CON*) blk->dia;
	int adjcolor = color [con->rank];

	if (adjcolor < mycolor) hi ++;
	else lo ++;
      }

      if (lo && hi) SET_Insert (&setmem, &middle, dia, NULL);
      else if (lo) SET_Insert (&setmem, &bottom, dia, NULL), size1 ++;
      else if (hi) SET_Insert (&setmem, &top, dia, NULL), size2 ++;
      else SET_Insert (&setmem, &internal, dia, NULL), size3 ++;
    }

    /* size1 + |int2| = size2 + |int1|
     * |int1| + |int2| = size3
     * -------------------------------
     * |int2| = (size3 + size2 - size1) / 2
     */

    size4 = (size3 + size2 - size1) / 2;
    size4 = MAX (0, size4);
    size5 = 0;

    /* create int1 and int2 such that: |bot| + |int2| = |top| + |int1| */
    for (SET *item = SET_First (internal); item; item = SET_Next (item))
    {
      dia = item->data;

      if (size5 < size4) SET_Insert (&setmem, &int2, dia, NULL), size5 ++; /* TODO: += |adj| rather than += 1 */
      else SET_Insert (&setmem, &int1, dia, NULL);
    }

    int sizes [5] = {SET_Size (bottom), SET_Size (middle), SET_Size (top), SET_Size (int1), SET_Size (int2)}, result [5];

    gs->bot = sizes [0];
    gs->mid = sizes [1];
    gs->top = sizes [2];
    gs->inn = sizes [3] + sizes [4];

    MPI_Reduce (sizes, result, 5, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0 && verbose) printf ("GAUSS_SEIDEL: |BOTTOM| = %d, |MIDDLE| = %d, |TOP| = %d, |INT1| = %d, |INT2| = %d\n",
				       result [0], result [1], result [2], result [3], result [4]);

    size_bot = size_top = 512;
    nsend_bot = nsend_top = 0;

    ERRMEM (send_bot = MEM_CALLOC (size_bot * sizeof (COMDATA)));
    ERRMEM (send_top = MEM_CALLOC (size_top * sizeof (COMDATA)));

    ptr_bot = send_bot;
    ptr_top = send_top;

    /* prepare bottom send buffer */
    for (SET *item = SET_First (bottom); item; item = SET_Next (item))
    {
      dia = item->data;

      for (ranks = NULL, blk = dia->adjext; blk; blk = blk->n)
      { 
	con = (CON*) blk->dia;
	SET_Insert (&setmem, &ranks, (void*) (long) con->rank, NULL);
      }

      for (SET *jtem = SET_First (ranks); jtem; jtem = SET_Next (jtem))
      {
	ptr_bot->rank = (int) (long) jtem->data;
	ptr_bot->ints = 1;
	ptr_bot->doubles = 3;
	ptr_bot->i = (int*) &dia->con->id;
	ptr_bot->d = dia->R;
	ptr_bot = sendnext (++ nsend_bot, &size_bot, &send_bot);
      }

      SET_Free (&setmem, &ranks);
    }

    /* prepare top send buffer */
    for (SET *item = SET_First (top); item; item = SET_Next (item))
    {
      dia = item->data;

      for (ranks = NULL, blk = dia->adjext; blk; blk = blk->n)
      { 
	con = (CON*) blk->dia;
	SET_Insert (&setmem, &ranks, (void*) (long) con->rank, NULL);
      }

      for (SET *jtem = SET_First (ranks); jtem; jtem = SET_Next (jtem))
      {
	ptr_top->rank = (int) (long) jtem->data;
	ptr_top->ints = 1;
	ptr_top->doubles = 3;
	ptr_top->i = (int*) &dia->con->id;
	ptr_top->d = dia->R;
	ptr_top = sendnext (++ nsend_top, &size_top, &send_top);
      }

      SET_Free (&setmem, &ranks);
    }

    bot_pattern = COM_Pattern (MPI_COMM_WORLD, TAG_GAUSS_SEIDEL_BOTTOM, send_bot, nsend_bot, &recv_bot, &nrecv_bot);
    top_pattern = COM_Pattern (MPI_COMM_WORLD, TAG_GAUSS_SEIDEL_TOP, send_top, nsend_top, &recv_top, &nrecv_top);

    if (gs->variant == GS_FULL)
    {
      size = 128;
      nsend = 0;
      ERRMEM (send = MEM_CALLOC (size * sizeof (COMDATA)));
      ptr = send;

      /* create send sets of external reactions updated by middle nodes */
      for (SET *item = SET_First (middle); item; item = SET_Next (item))
      {
	dia = item->data;

	for (ranks = NULL, blk = dia->adjext; blk; blk = blk->n)
	{ 
	  con = (CON*) blk->dia;
	  SET_Insert (&setmem, &ranks, (void*) (long) con->rank, NULL);
	}

	for (SET *jtem = SET_First (ranks); jtem; jtem = SET_Next (jtem))
	{
	  ptr->rank = (int) (long) jtem->data;
	  ptr->ints = 1;
	  ptr->doubles = 0;
	  ptr->i = (int*) &dia->con->id;
	  ptr->d = NULL;
	  ptr= sendnext (++ nsend, &size, &send);
	}

	SET_Free (&setmem, &ranks);
      }

      /* send ranks of middle nodes */
      COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);

      /* discover which external constraints are updated by middle nodes */
      receive_middle_reactions (dom, recv, nrecv, &setmem, &midupd);

      free (send);
      free (recv);
    }
    else /* GS_MIDDLE_JACOBI */
    {
      size_mid = 512;
      nsend_mid = 0;

      ERRMEM (send_mid = MEM_CALLOC (size_mid * sizeof (COMDATA)));

      ptr_mid = send_mid;

      /* prepare middle send buffer */
      for (SET *item = SET_First (middle); item; item = SET_Next (item))
      {
	dia = item->data;

	for (ranks = NULL, blk = dia->adjext; blk; blk = blk->n)
	{ 
	  con = (CON*) blk->dia;
	  SET_Insert (&setmem, &ranks, (void*) (long) con->rank, NULL);
	}

	for (SET *jtem = SET_First (ranks); jtem; jtem = SET_Next (jtem))
	{
	  ptr_mid->rank = (int) (long) jtem->data;
	  ptr_mid->ints = 1;
	  ptr_mid->doubles = 3;
	  ptr_mid->i = (int*) &dia->con->id;
	  ptr_mid->d = dia->R;
	  ptr_mid = sendnext (++ nsend_mid, &size_mid, &send_mid);
	}

	SET_Free (&setmem, &ranks);
      }

      mid_pattern = COM_Pattern (MPI_COMM_WORLD, TAG_GAUSS_SEIDEL_BOTTOM, send_mid, nsend_mid, &recv_mid, &nrecv_mid);
    }
  }
  else if (gs->variant == GS_BOUNDARY_JACOBI)
  {
    for (dia = ldy->dia; dia; dia = dia->n)
    {
      SET_Insert (&setmem, &all, dia, NULL);
    }
  }

  dynamic = dom->dynamic;
  step = dom->step;
  gs->error = GS_OK;
  gs->iters = 0;
  dimax = 0;

  E("GSINIT");

  do
  {
    double errup = 0.0,
	   errlo = 0.0,
	   errloc [2],
	   errsum [2];

    S("GSRUN"); undo_all (ldy); E("GSRUN"); 

    if (gs->reverse && gs->iters % 2)
    {
      if (gs->variant != GS_BOUNDARY_JACOBI)
      {
	S("GSRUN"); di = gauss_seidel_sweep (bottom, 1, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
	S("GSCOM"); COM_Send (bot_pattern); E("GSCOM");
	S("GSRUN"); di = gauss_seidel_sweep (int1, 1, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
	S("GSCOM"); COM_Recv (bot_pattern); receive_reactions (dom, recv_bot, nrecv_bot); E("GSCOM");

	if (gs->variant == GS_FULL)
	{
	  di = gauss_seidel_loop (middle, midupd, 1, &setmem, mycolor, color, gs, ldy, dynamic, step, &errup, &errlo); dimax = MAX (dimax, di);
	}
	else /* GS_MIDDLE_JACOBI */
	{
	  S("GSRUN"); di = gauss_seidel_sweep (middle, 1, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
	  S("GSCOM"); COM_Repeat (mid_pattern); receive_reactions (dom, recv_mid, nrecv_mid); E("GSCOM");
	}

	S("GSRUN"); di = gauss_seidel_sweep (top, 1, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
	S("GSCOM"); COM_Send (top_pattern); E("GSCOM");
	S("GSRUN"); di = gauss_seidel_sweep (int2, 1, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
	S("GSCOM"); COM_Recv (top_pattern); receive_reactions (dom, recv_top, nrecv_top); E("GSCOM");
      }
      else
      {
	S("GSRUN"); di = gauss_seidel_sweep (all, 1, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
      }
    }
    else
    {
      if (gs->variant != GS_BOUNDARY_JACOBI)
      {
	S("GSRUN"); di = gauss_seidel_sweep (top, 0, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
	S("GSCOM"); COM_Send (top_pattern); E("GSCOM");
	S("GSRUN"); di = gauss_seidel_sweep (int2, 0, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN"); /* large |top| => large |int2| */
	S("GSCOM"); COM_Recv (top_pattern); receive_reactions (dom, recv_top, nrecv_top); E("GSCOM");

	if (gs->variant == GS_FULL)
	{
	  di = gauss_seidel_loop (middle, midupd, 0, &setmem, mycolor, color, gs, ldy, dynamic, step, &errup, &errlo); dimax = MAX (dimax, di);
	}
	else /* GS_MIDDLE_JACOBI */
	{
	  S("GSRUN"); di = gauss_seidel_sweep (middle, 0, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
	  S("GSCOM"); COM_Repeat (mid_pattern); receive_reactions (dom, recv_mid, nrecv_mid); E("GSCOM");
	}

	S("GSRUN"); di = gauss_seidel_sweep (bottom, 0, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
	S("GSCOM"); COM_Send (bot_pattern); E("GSCOM");
	S("GSRUN"); di = gauss_seidel_sweep (int1, 0, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
	S("GSCOM"); COM_Recv (bot_pattern); receive_reactions (dom, recv_bot, nrecv_bot); E("GSCOM");
      }
      else
      {
	S("GSRUN"); di = gauss_seidel_sweep (all, 0, gs, dynamic, step, gs->innerloops, &errup, &errlo); dimax = MAX (dimax, di); E("GSRUN");
      }
    }

    if (gs->variant == GS_BOUNDARY_JACOBI) DOM_Update_External_Reactions (dom, 0);
#if DEBUG
    else ASSERT_DEBUG (all_done (ldy), "Not all external reactions were updated");
#endif

    /* merit function */
    if (!nomerit)
    {
      S("GSRUN"); *merit = MERIT_Function (ldy, 1); E("GSRUN"); 
    }

    /* sum up error */
    S("GSCOM"); 
    errloc [0] = errup, errloc [1] = errlo;
    MPI_Allreduce (errloc, errsum, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    errup = errsum [0], errlo = errsum [1];
    E("GSCOM"); 

    /* calculate relative error */
    error = sqrt (errup) / sqrt (errlo == 0.0 ? 1.0 : errlo);

    /* record values */
    gs->rerhist [gs->iters] = error;
    gs->merhist [gs->iters] = *merit;

    if (gs->iters % div == 0 && rank == 0 && verbose) printf (fmt, gs->iters, error, *merit), div *= 2;
  }
  while (++ gs->iters < gs->maxiter && (error > gs->epsilon || *merit > gs->meritval));

  if (gs->itershistcount >= 0)
  {
    if (gs->itershistcount >= gs->itershistsize)
    {
      gs->itershistsize += 1024;
      ERRMEM (gs->itershist = realloc (gs->itershist, gs->itershistsize * sizeof(int)));
    }
    gs->itershist[gs->itershistcount] = gs->iters;
    gs->itershistcount ++;
  }

  if (rank == 0 && verbose) printf (fmt, gs->iters, error, *merit);

  if (gs->variant < GS_BOUNDARY_JACOBI)
  {
    COM_Free (bot_pattern);
    COM_Free (top_pattern);
    free (send_bot);
    free (recv_bot);
    free (send_top);
    free (recv_top);
    free (color);

    if (gs->variant == GS_MIDDLE_JACOBI)
    {
      COM_Free (mid_pattern);
      free (send_mid);
      free (recv_mid);
    }
  }
  MEM_Release (&setmem);

  /* get maximal iterations count of a diagonal block solver (this has been
   * delayed until here to minimize small communication within the loop) */
  MPI_Allreduce (&dimax, &diagiters, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  if (diagiters >= gs->diagmaxiter || diagiters < 0)
  {
    if (diagiters < 0) gs->error = GS_DIAGONAL_FAILED;
    else gs->error = GS_DIAGONAL_DIVERGED;

    switch ((int) gs->failure)
    {
    case GS_FAILURE_EXIT:
      THROW (ERR_GAUSS_SEIDEL_DIAGONAL_DIVERGED);
      break;
    case GS_FAILURE_CALLBACK:
      gs->callback (gs->data);
      break;
    }
  }
  else if (gs->iters >= gs->maxiter)
  {
    gs->error = GS_DIVERGED;

    switch (gs->failure)
    {
    case GS_FAILURE_CONTINUE:
      break;
    case GS_FAILURE_EXIT:
      THROW (ERR_GAUSS_SEIDEL_DIVERGED);
      break;
    case GS_FAILURE_CALLBACK:
      gs->callback (gs->data);
      break;
    }
  }
}
#else
/* run serial solver */
void GAUSS_SEIDEL_Solve (GAUSS_SEIDEL *gs, LOCDYN *ldy)
{
  double error, *merit, step;
  int verbose, diagiters;
  short dynamic, nomerit;
  char fmt [512];
  int div = 10;
  DIAB *end;

  S("GSRUN");

  verbose = ldy->dom->verbose && gs->verbose;
  merit = &ldy->dom->merit;

  nomerit = gs->nomerit ? 1 : gs->meritval >= 1.0 ? 1 : 0;

  if (nomerit) *merit = 0.0;

  if (verbose) sprintf (fmt, "GAUSS_SEIDEL: iteration: %%%dd  error:  %%.2e  merit:  %%.2e\n", (int)log10 (gs->maxiter) + 1);

  gs->rerhist = realloc (gs->rerhist, gs->maxiter * sizeof (double));
  gs->merhist = realloc (gs->merhist, gs->maxiter * sizeof (double));

  if (gs->reverse && ldy->dia) for (end = ldy->dia; end->n; end = end->n); /* find last block for the backward run */
  else end = NULL;

  dynamic = ldy->dom->dynamic;
  step = ldy->dom->step;
  gs->error = GS_OK;
  gs->iters = 0;
  do
  {
    double errup = 0.0,
	   errlo = 0.0;
    OFFB *blk;
    DIAB *dia;
   
    for (dia = end && gs->iters % 2 ? end : ldy->dia; dia; dia = end && gs->iters % 2 ? dia->p : dia->n) /* run forward and backward alternately */
    {
      double R0 [3],
	     B [3],
	     *R = dia->R;

      /* compute local free velocity */
      COPY (dia->B, B);
      for (blk = dia->adj; blk; blk = blk->n)
      {
	double *W = blk->W,
	       *R = blk->dia->R;
	NVADDMUL (B, W, R, B);
      }
      
      COPY (R, R0); /* previous reaction */

      /* solve local diagonal block problem */
      CON *con = dia->con;
      diagiters = DIAGONAL_BLOCK_Solver (gs->diagsolver, gs->diagepsilon, gs->diagmaxiter, dynamic,
	                step, con->kind, &con->mat, con->gap, con->area, con->Z, con->base, dia, B);

      if (diagiters >= gs->diagmaxiter || diagiters < 0)
      {
	if (diagiters < 0) gs->error = GS_DIAGONAL_FAILED;
	else gs->error = GS_DIAGONAL_DIVERGED;

	switch (gs->failure)
	{
	case GS_FAILURE_CONTINUE:

	  if (con->kind == CONTACT)
	  {
	    DIAS dias [4] = {DS_SEMISMOOTH_NEWTON, DS_PROJECTED_GRADIENT, DS_DE_SAXCE_FENG, DS_PROJECTED_NEWTON};

	    for (int i = 0; i < 4; i ++)
	    {
	      if (dias [i] != gs->diagsolver) /* skip current diagonal solver */
	      {
		COPY (R0, R); /* initialize with previous reaction */

		diagiters = DIAGONAL_BLOCK_Solver (dias [i], gs->diagepsilon, gs->diagmaxiter, /* try another solver */
		  dynamic, step, con->kind, &con->mat, con->gap, con->area, con->Z, con->base, dia, B);

		if (diagiters < gs->diagmaxiter && diagiters >= 0) break; /* success */
	      }
	    }
	  }

	  if (diagiters >= gs->diagmaxiter || diagiters < 0) /* failed */
	  {
	    COPY (R0, R); /* use previous reaction */
	  }

	  break;
	case GS_FAILURE_EXIT:
	  THROW (ERR_GAUSS_SEIDEL_DIAGONAL_DIVERGED);
	  break;
	case GS_FAILURE_CALLBACK:
	  gs->callback (gs->data);
	  break;
	}
      }

      /* accumulate relative
       * error components */
      SUB (R, R0, R0);
      errup += DOT (R0, R0);
      errlo += DOT (R, R);
    }

    /* merit function value */
    if (!nomerit)
    {
      *merit = MERIT_Function (ldy, 1);
    }

    /* calculate relative error */
    error = sqrt (errup) / sqrt (errlo == 0.0 ? 1.0 : errlo);

    /* record values */
    gs->rerhist [gs->iters] = error;
    gs->merhist [gs->iters] = *merit;

    if (gs->iters % div == 0 && verbose) printf (fmt, gs->iters, error, *merit), div *= 2;
  }
  while (++ gs->iters < gs->maxiter && (error > gs->epsilon || *merit > gs->meritval));

  if (gs->itershistcount >= 0)
  {
    if (gs->itershistcount >= gs->itershistsize)
    {
      gs->itershistsize += 1024;
      ERRMEM (gs->itershist = realloc (gs->itershist, gs->itershistsize * sizeof(int)));
    }
    gs->itershist[gs->itershistcount] = gs->iters;
    gs->itershistcount ++;
  }

  if (verbose) printf (fmt, gs->iters, error, *merit);

  E("GSRUN");

  if (gs->iters >= gs->maxiter)
  {
    gs->error = GS_DIVERGED;

    switch (gs->failure)
    {
    case GS_FAILURE_CONTINUE:
      break;
    case GS_FAILURE_EXIT:
      THROW (ERR_GAUSS_SEIDEL_DIVERGED);
      break;
    case GS_FAILURE_CALLBACK:
      gs->callback (gs->data);
      break;
    }
  }
}
#endif

/* return faulure string */
char* GAUSS_SEIDEL_Failure (GAUSS_SEIDEL *gs)
{
  switch (gs->failure)
  {
  case GS_FAILURE_CONTINUE: return "FAILURE_CONTINUE";
  case GS_FAILURE_EXIT: return "FAILURE_EXIT";
  case GS_FAILURE_CALLBACK: return "FAILURE_CALLBACK";
  }

  return NULL;
}

/* return diagonal solver string */
char* GAUSS_SEIDEL_Diagsolver (GAUSS_SEIDEL *gs)
{
  switch (gs->diagsolver)
  {
  case DS_PROJECTED_GRADIENT: return "PROJECTED_GRADIENT";
  case DS_DE_SAXCE_FENG: return "DE_SAXCE_FENG";
  case DS_SEMISMOOTH_NEWTON: return "SEMISMOOTH_NEWTON";
  case DS_PROJECTED_NEWTON: return "PROJECTED_NEWTON";
  }

  return NULL;
}

/* return error string */
char* GAUSS_SEIDEL_Error (GAUSS_SEIDEL *gs)
{
  switch (gs->error)
  {
  case GS_OK: return "OK";
  case GS_DIVERGED: return "DIVERGED";
  case GS_DIAGONAL_DIVERGED: return "DIAGONAL_DIVERGED";
  case GS_DIAGONAL_FAILED: return "DIAGONAL_FAILED";
  }

  return NULL;
}

/* return reverse flag string */
char* GAUSS_SEIDEL_Reverse (GAUSS_SEIDEL *gs)
{
  switch (gs->reverse)
  {
  case GS_ON: return "ON";
  case GS_OFF: return "OFF";
  }

  return NULL;
}

/* return variant string */
char* GAUSS_SEIDEL_Variant (GAUSS_SEIDEL *gs)
{
  switch (gs->variant)
  {
  case GS_FULL: return "FULL";
  case GS_MIDDLE_JACOBI: return "MIDDLE_JACOBI";
  case GS_BOUNDARY_JACOBI: return "BOUNDARY_JACOBI";
  }

  return NULL;
}

/* write labeled satate values */
void GAUSS_SEIDEL_Write_State (GAUSS_SEIDEL *gs, PBF *bf)
{
  PBF_Label (bf, "GSITERS");
  PBF_Int (bf, &gs->iters, 1);
#if MPI
  PBF_Label (bf, "GSCOLORS");
  PBF_Int (bf, &gs->colors, 1);
  PBF_Label (bf, "GSBOT");
  PBF_Int (bf, &gs->bot, 1);
  PBF_Label (bf, "GSMID");
  PBF_Int (bf, &gs->mid, 1);
  PBF_Label (bf, "GSTOP");
  PBF_Int (bf, &gs->top, 1);
  PBF_Label (bf, "GSINN");
  PBF_Int (bf, &gs->inn, 1);
#endif
}

/* free solver */
void GAUSS_SEIDEL_Destroy (GAUSS_SEIDEL *gs)
{
  free (gs->rerhist);
  free (gs->merhist);
  free (gs->itershist);
  free (gs);
}
