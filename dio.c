/*
 * dio.c
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * domain input-output
 */

#include <string.h>
#include "sol.h"
#include "dio.h"
#include "pck.h"
#include "err.h"

/* write constraint state */
static void write_constraint (CON *con, PBF *bf)
{
  unsigned int zero = 0;
  int kind = con->kind;

  PBF_Uint (bf, &con->id, 1);
  PBF_Int (bf, &kind, 1);

  PBF_Double (bf, con->R, 3);
  PBF_Double (bf, con->U, 3);
  PBF_Double (bf, con->point, 3);
  PBF_Double (bf, con->base, 9);
  PBF_Double (bf, &con->merit, 1);

  PBF_Uint (bf, &con->master->id, 1);
  if (con->slave) PBF_Uint (bf, &con->slave->id, 1);
  else PBF_Uint (bf, &zero, 1);

  if (kind == CONTACT)
  {
    SURFACE_MATERIAL_Write_State (&con->mat, bf);
    PBF_Double (bf, &con->area, 1);
    PBF_Double (bf, &con->gap, 1);
  }

  if (kind == RIGLNK || kind == VELODIR) PBF_Double (bf, con->Z, DOM_Z_SIZE);

#if MPI
  PBF_Int (bf, &con->master->dom->rank, 1);
#endif
}

/* read constraint state */
static CON* read_constraint (DOM *dom, int iover, PBF *bf)
{
  unsigned int id;
  CON *con;
  int kind;

  ERRMEM (con = MEM_Alloc (&dom->conmem));

  PBF_Uint (bf, &con->id, 1);
  PBF_Int (bf, &kind, 1);
  con->kind = kind;

  PBF_Double (bf, con->R, 3);
  PBF_Double (bf, con->U, 3);
  PBF_Double (bf, con->point, 3);
  PBF_Double (bf, con->base, 9);
  PBF_Double (bf, &con->merit, 1);

  PBF_Uint (bf, &id, 1);
  ASSERT_DEBUG_EXT (con->master = MAP_Find (dom->idb, (void*) (long) id, NULL), "Invalid master id");
  PBF_Uint (bf, &id, 1);
  if (id) ASSERT_DEBUG_EXT (con->slave = MAP_Find (dom->idb, (void*) (long) id, NULL), "Invalid slave id");

  if (kind == CONTACT)
  {
    con->state |= SURFACE_MATERIAL_Read_State (dom->sps, &con->mat, bf);
    PBF_Double (bf, &con->area, 1);
    PBF_Double (bf, &con->gap, 1);
  }

  if (kind == RIGLNK || kind == VELODIR) PBF_Double (bf, con->Z, DOM_Z_SIZE);

  if (dom->solfec->ioparallel)
  {
    PBF_Int (bf, &con->rank, 1);
  }

  return con;
}

/* attach constraints to bodies after reading */
static void dom_attach_constraints (DOM *dom)
{
  BODY *bod;
  CON *con;

  for (bod = dom->bod; bod; bod = bod->next) SET_Free (&dom->setmem, &bod->con);

  for (con = dom->con; con; con = con->next)
  {
    if (con->master) SET_Insert (&dom->setmem, &con->master->con, con, NULL);

    if (con->slave) SET_Insert (&dom->setmem, &con->slave->con, con, NULL);
  }
}

/* write uncompressed domain state */
void dom_write_state (DOM *dom, PBF *bf)
{
  /* mark domain output */

  PBF_Label (bf, "DOM");

  /* write time step */

  PBF_Label (bf, "STEP");

  PBF_Double (bf, &dom->step, 1);

  /* write contacts */

  PBF_Label (bf, "CONS");
 
  PBF_Int (bf, &dom->ncon, 1);

  for (CON *con = dom->con; con; con = con->next)
  {
    write_constraint (con, bf);
  }

  /* write ids of bodies that have been deleted and empty the deleted bodies ids set */

  unsigned int id;
  SET *item;
  int size;

  PBF_Label (bf, "DELBODS");

  size = SET_Size (dom->delb);

  PBF_Int (bf, &size, 1);

  for (item = SET_First (dom->delb); item; item = SET_Next (item))
  {
    id = (int) (long) item->data;
    PBF_Uint (bf, &id, 1);
  }

  SET_Free (&dom->setmem, &dom->delb);

  /* write complete data of newly created bodies and empty the newly created bodies set */

  int dsize = 0, doubles;
  double *d = NULL;

  int isize = 0, ints, *i = NULL;

  PBF_Label (bf, "NEWBODS");

  size = SET_Size (dom->newb);

  PBF_Int (bf, &size, 1);

  for (item = SET_First (dom->newb); item; item = SET_Next (item))
  {
    doubles = ints = 0;

    BODY_Pack (item->data, &dsize, &d, &doubles, &isize, &i, &ints);

    PBF_Int (bf, &doubles, 1);
    PBF_Double (bf, d, doubles);

    PBF_Int (bf, &ints, 1);
    PBF_Int (bf, i, ints);
  }

  free (d);
  free (i);

  SET_Free (&dom->setmem, &dom->newb);

  /* write regular bodies (this also includes states of newly created ones) */

  PBF_Label (bf, "BODS");

  PBF_Int (bf, &dom->nbod, 1);

  for (BODY *bod = dom->bod; bod; bod = bod->next)
  {
    PBF_Uint (bf, &bod->id, 1);

    if (bod->label) PBF_Label (bf, bod->label); /* label body record for fast access */

    BODY_Write_State (bod, bf);
  }
}

/* read uncompressed domain state */
void dom_read_state (DOM *dom, PBF *bf)
{
  int iover, ncon;

  /* input version */
  iover = dom->solfec->iover;

  /* clear contacts */
  MAP_Free (&dom->mapmem, &dom->idc);
  MEM_Release (&dom->conmem);
  dom->con = NULL;
  dom->ncon = 0;

  for (; bf; bf = bf->next)
  {
    if (PBF_Label (bf, "DOM"))
    {
      /* read time step */

      ASSERT (PBF_Label (bf, "STEP"), ERR_FILE_FORMAT);

      PBF_Double (bf, &dom->step, 1);

      /* read contacts */

      ASSERT (PBF_Label (bf, "CONS"), ERR_FILE_FORMAT);
    
      PBF_Int (bf, &ncon, 1);

      for (int n = 0; n < ncon; n ++)
      {
	CON *con;
	
	con = read_constraint (dom, iover, bf);
	MAP_Insert (&dom->mapmem, &dom->idc, (void*) (long) con->id, con, NULL);
	con->next = dom->con;
	if (dom->con) dom->con->prev = con;
	dom->con = con;
      }

      dom->ncon += ncon;

      /* read ids of bodies that need to be deleted and remove them from all containers */

      unsigned int id;
      int size, n;

      ASSERT (PBF_Label (bf, "DELBODS"), ERR_FILE_FORMAT);

      PBF_Int (bf, &size, 1);

      for (n = 0; n < size; n ++)
      {
	BODY *bod;

	PBF_Uint (bf, &id, 1);

	ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) id, NULL), "Invalid body id");

	if (bod->label) MAP_Delete (&dom->mapmem, &dom->lab, bod->label, (MAP_Compare) strcmp);
	MAP_Delete (&dom->mapmem, &dom->idb, (void*) (long) id, NULL);
	if (bod->next) bod->next->prev = bod->prev;
	if (bod->prev) bod->prev->next = bod->next;
	else dom->bod = bod->next;
	dom->nbod --;
	BODY_Destroy (bod);
      }

      /* read complete data of newly created bodies and insert them into all containers */

      int dpos, doubles;
      double *d = NULL;

      int ipos, ints, *i = NULL;

      ASSERT (PBF_Label (bf, "NEWBODS"), ERR_FILE_FORMAT);

      PBF_Int (bf, &size, 1);

      for (n = 0; n < size; n ++)
      {
	BODY *bod;

	PBF_Int (bf, &doubles, 1);
	ERRMEM (d  = realloc (d, doubles * sizeof (double)));
	PBF_Double (bf, d, doubles);

	PBF_Int (bf, &ints, 1);
	ERRMEM (i  = realloc (i, ints * sizeof (int)));
	PBF_Int (bf, i, ints);

	dpos = ipos = 0;

	bod = BODY_Unpack (dom->solfec, &dpos, d, doubles, &ipos, i, ints);

	if (bod->label) MAP_Insert (&dom->mapmem, &dom->lab, bod->label, bod, (MAP_Compare) strcmp);
	MAP_Insert (&dom->mapmem, &dom->idb, (void*) (long) bod->id, bod, NULL);
	bod->next = dom->bod;
	if (dom->bod) dom->bod->prev = bod;
	dom->bod = bod;
	dom->nbod ++;
      }

      free (d);
      free (i);

      /* read regular bodies */

      ASSERT (PBF_Label (bf, "BODS"), ERR_FILE_FORMAT);

      int nbod;

      PBF_Int (bf, &nbod, 1);

      for (int n = 0; n < nbod; n ++)
      {
	unsigned int id;
	BODY *bod;

	PBF_Uint (bf, &id, 1);
	ASSERT_DEBUG_EXT (bod = MAP_Find (dom->idb, (void*) (long) id, NULL), "Body id invalid");
	BODY_Read_State (bod, bf);
      }
    }
  }

  dom_attach_constraints (dom); /* attach constraints to bodies */
}

/* read uncompressed state of an individual body */
int dom_read_body (DOM *dom, PBF *bf, BODY *bod)
{
  if (bod->label)
  {
    for (; bf; bf = bf->next)
    {
      if (PBF_Label (bf, bod->label))
      {
	BODY_Read_State (bod, bf);
	return 1;
      }
    }
  }
  else
  {
    for (; bf; bf = bf->next)
    {
      if (PBF_Label (bf, "BODS"))
      {
	int nbod;

	PBF_Int (bf, &nbod, 1);

	for (int n = 0; n < nbod; n ++)
	{
	  unsigned int id;
	  BODY *obj;

	  PBF_Uint (bf, &id, 1);
	  ASSERT_DEBUG_EXT (obj = MAP_Find (dom->idb, (void*) (long) id, NULL), "Body id invalid");
	  if (bod->id == obj->id) 
	  {
	    BODY_Read_State (bod, bf);
	    return 1;
	  }
	  else /* skip body and continue */
	  {
	    BODY fake;

	    ERRMEM (fake.conf = malloc (sizeof (double [BODY_Conf_Size (obj)])));
	    ERRMEM (fake.velo = malloc (sizeof (double [obj->dofs])));
	    fake.shape = NULL;

	    BODY_Read_State (&fake, bf);

	    free (fake.conf);
	    free (fake.velo);
	  }
	}
      }
    }
  }

  return 0;
}

/* read uncompressed state of an individual constraint */
int dom_read_constraint (DOM *dom, PBF *bf, CON *con)
{
  int iover = dom->solfec->iover;

  for (; bf; bf = bf->next)
  {
    int ncon;

    if (PBF_Label (bf, "CONS"))
    {
      PBF_Int (bf, &ncon, 1);

      for (int n = 0; n < ncon; n ++)
      {
	CON *obj = read_constraint (dom, iover, bf);

	if (con->id == obj->id)
	{
	  *con = *obj;
          MEM_Free (&dom->conmem, obj); /* not needed */
	  return 1;
	}
	else MEM_Free (&dom->conmem, obj); /* skip and continue */
      }
    }
  }

  return 0;
}
