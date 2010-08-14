/*
 * bss.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * body space solver
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
#include <float.h>

#include "bss.h"
#include "alg.h"
#include "err.h"
#include "mrf.h"

#if 0
/* update local velocities of remote parent constraints, attached to local parent bodies */
static void update_remote_local_velocities (DOM *dom)
{
  int nsend, nrecv, *isize, *dsize, i, *j, *k;
  COMDATA *send, *recv, *ptr;
  double X [6], *D, *V, *B;
  BODY *bod;
  SET *item;
  CON *con;

  nsend = dom->ncpu;
  ERRMEM (send = MEM_CALLOC (sizeof (COMDATA [nsend])));
  ERRMEM (isize = MEM_CALLOC (sizeof (int [nsend])));
  ERRMEM (dsize = MEM_CALLOC (sizeof (int [nsend])));

  for (i = 0; i < nsend; i ++) send [i].rank = i;

  for (bod = dom->bod; bod; bod = bod->next) /* for all parents */
  {
    for (item = SET_First (bod->con); item; item = SET_Next (item))
    {
      con = item->data;
      if (con->state & CON_EXTERNAL) /* needs local velocity update */
      {
	if (!MAP_Find (dom->conext, (void*) (long) con->id, NULL))
	{
	  CON *tst;
	  if ((tst = MAP_Find (dom->idc, (void*) (long) con->id, NULL)))
	  {
	    printf ("Time %g, tank %d, ", dom->time, dom->rank);
	    if (tst != con) printf ("External found as a different regular\n");
	    else  printf ("External found as same regular\n");
	  }
	}
	ASSERT_DEBUG (MAP_Find (dom->conext, (void*) (long) con->id, NULL), "Invalid external constraint %s %d", CON_Kind (con), con->id);
	i = con->rank;
	ptr = &send [i];
	if (bod == con->master)
	{
	  pack_int (&isize [i], &ptr->i, &ptr->ints, con->id);
	  BODY_Local_Velo (bod, mshp(con), mgobj(con), con->mpnt, con->base, X, X+3);
	}
	else
	{
	  pack_int (&isize [i], &ptr->i, &ptr->ints, -con->id);
	  BODY_Local_Velo (bod, sshp(con), sgobj(con), con->spnt, con->base, X, X+3);
	}
        pack_doubles (&dsize [i], &ptr->d, &ptr->doubles, X, 6);
      }
    }
  }

  dom->bytes += COMALL (MPI_COMM_WORLD, send, nsend, &recv, &nrecv);

  for (i = 0; i < nrecv; i ++)
  {
    ptr = &recv [i];
    for (j = ptr->i, k = j + ptr->ints, D = ptr->d; j < k; j ++, D += 6)
    {
      ASSERT_DEBUG_EXT (con = MAP_Find (dom->idc, (void*) (long) ABS (*j), NULL), "Invalid constraint id: %d", *j);
      V = con->dia->V;
      B = con->dia->B;
      if ((*j) < 0) {V += 3; B += 3;} /* slave */
      COPY (D, V);
      COPY (D+3, B);
    }
  }

  for (i = 0; i < nsend; i ++)
  {
    ptr = &send [i];
    free (ptr->i);
    free (ptr->d);
  }
  free (send);
  free (isize);
  free (dsize);
  free (recv); /* includes recv[]->i and recv[]->d memory */
}
#endif

/* create solver */
BSS* BSS_Create (int maxiter, double meritval)
{
  BSS *bs;

  ERRMEM (bs = MEM_CALLOC (sizeof (BSS)));
  bs->meritval = meritval;
  bs->maxiter = maxiter;

  return bs;
}

/* run solver */
void BSS_Solve (BSS *bs, LOCDYN *ldy)
{
}

/* write labeled state values */
void BSS_Write_State (BSS *bs, PBF *bf)
{
}

/* destroy solver */
void BSS_Destroy (BSS *bs)
{
  free (bs);
}
