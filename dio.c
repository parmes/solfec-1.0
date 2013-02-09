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
#if IOVER > 1
  if (kind == CONTACT)
  {
    PBF_Double (bf, con->V, 3);
  }
#endif
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
    PBF_Int (bf, con->spair, 2);
  }

  if (kind == RIGLNK || kind == VELODIR) PBF_Double (bf, con->Z, DOM_Z_SIZE);

#if MPI
  PBF_Int (bf, &con->master->dom->rank, 1);
#endif
}

/* read constraint state */
static CON* read_constraint (DOM *dom, int iover, PBF *bf)
{
  SOLFEC_MODE mode = dom->solfec->mode;
  unsigned int id;
  CON *con;
  int kind;

  ERRMEM (con = MEM_Alloc (&dom->conmem));

  PBF_Uint (bf, &con->id, 1);
  PBF_Int (bf, &kind, 1);
  con->kind = kind;

  PBF_Double (bf, con->R, 3);
  PBF_Double (bf, con->U, 3);
  if (iover > 1 && kind == CONTACT)
  {
    PBF_Double (bf, con->V, 3);
  }
  PBF_Double (bf, con->point, 3);
  PBF_Double (bf, con->base, 9);
  PBF_Double (bf, &con->merit, 1);

  PBF_Uint (bf, &id, 1);
  ASSERT_DEBUG_EXT (con->master = MAP_Find (dom->allbodies, (void*) (long) id, NULL), "Invalid master id");
  PBF_Uint (bf, &id, 1);
  if (id) ASSERT_DEBUG_EXT (con->slave = MAP_Find (dom->allbodies, (void*) (long) id, NULL), "Invalid slave id");

  if (kind == CONTACT)
  {
    con->state |= SURFACE_MATERIAL_Read_State (dom->sps, &con->mat, bf);
    PBF_Double (bf, &con->area, 1);
    PBF_Double (bf, &con->gap, 1);
    PBF_Int (bf, con->spair, 2);
  }

  if (kind == RIGLNK || kind == VELODIR) PBF_Double (bf, con->Z, DOM_Z_SIZE);

  if (bf->parallel == PBF_ON)
  {
    if (mode == SOLFEC_READ)
    {
      PBF_Int (bf, &con->rank, 1);
    }
    else /* fake it => ranks are actually used in WRITE mode */
    {
      int rank;
      PBF_Int (bf, &rank, 1);
    }
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
    if (con->master)
    {
      ASSERT_DEBUG (MAP_Find_Node (dom->idb, (void*) (long) con->master->id, NULL), "Invalid master id");
      SET_Insert (&dom->setmem, &con->master->con, con, NULL);
    }

    if (con->slave)
    {
      ASSERT_DEBUG (MAP_Find_Node (dom->idb, (void*) (long) con->slave->id, NULL), "Invalid slave id");
      SET_Insert (&dom->setmem, &con->slave->con, con, NULL);
    }
  }
}

/* write new bodies data */
static void write_new_bodies (DOM *dom)
{
  char *path, *ext;
  FILE *file;
  XDR xdr;

  if (dom->newb == NULL) return; /* nothing to write */

  path = SOLFEC_Alloc_File_Name (dom->solfec, 16);
  ext = path + strlen (path);

#if MPI
  sprintf (ext, ".bod.%d", dom->rank);
#else
  sprintf (ext, ".bod");
#endif

  ASSERT (file = fopen (path, "a"), ERR_FILE_OPEN);
  xdrstdio_create (&xdr, file, XDR_ENCODE);

  int isize = 0, ints, *i = NULL;
  int dsize = 0, doubles;
  double *d = NULL;
  SET *item;

  for (item = SET_First (dom->newb); item; item = SET_Next (item))
  {
    doubles = ints = 0;

    BODY_Pack (item->data, &dsize, &d, &doubles, &isize, &i, &ints);

    ASSERT (xdr_int (&xdr, &doubles), ERR_PBF_WRITE);
    ASSERT (xdr_vector (&xdr, (char*)d, doubles, sizeof (double), (xdrproc_t)xdr_double), ERR_PBF_WRITE);

    ASSERT (xdr_int (&xdr, &ints), ERR_PBF_WRITE);
    ASSERT (xdr_vector (&xdr, (char*)i, ints, sizeof (int), (xdrproc_t)xdr_int), ERR_PBF_WRITE);
  }

  free (d);
  free (i);

  xdr_destroy (&xdr);
  fclose (file);
  free (path);
}

/* read new bodies data */
static void read_new_bodies (DOM *dom, PBF *bf)
{
  char *path, *ext;
  FILE *file;
  int m, n;
  XDR xdr;

  if (dom->solfec->verbose)
    printf ("Reading all bodies ...\n");

  dom->allbodiesread = 1; /* mark as read */

  path = SOLFEC_Alloc_File_Name (dom->solfec, 16);
  ext = path + strlen (path);
  
  for (m = 0; bf; bf = bf->next) m ++; /* count input files */

  for (n = 0; n < m; n ++)
  {
    if (n || m > 1)
    {
      sprintf (ext, ".bod.%d", n);
      if (!(file = fopen (path, "r"))) continue; /* no new bodies for this rank */
    }
    else /* n == 0 && m == 1 */
    {
      sprintf (ext, ".bod.%d", n);
      if (!(file = fopen (path, "r"))) /* either prallel with "mpirun -np 1" */
      {
	sprintf (ext, ".bod");
	if (!(file = fopen (path, "r"))) continue; /* or serial */
      }
    }

    xdrstdio_create (&xdr, file, XDR_DECODE);

    int ipos, ints, *i, dpos, doubles;
    double *d;
    BODY *bod;

    for (;;)
    {
      if (xdr_int (&xdr, &doubles))
      {
	ERRMEM (d = malloc (sizeof (double [doubles])));
	if (xdr_vector (&xdr, (char*)d, doubles, sizeof (double), (xdrproc_t)xdr_double))
	{
	  if (xdr_int (&xdr, &ints))
	  {
	    ERRMEM (i = malloc (sizeof (int [ints])));
	    if (xdr_vector (&xdr, (char*)i, ints, sizeof (int), (xdrproc_t)xdr_int))
	    {
	      ipos = dpos = 0;

	      bod = BODY_Unpack (dom->solfec, &dpos, d, doubles, &ipos, i, ints);

	      if (!MAP_Find (dom->allbodies, (void*) (long) bod->id, NULL))
	      {
	        MAP_Insert (&dom->mapmem, &dom->allbodies, (void*) (long) bod->id, bod, NULL);
	      }
	      else BODY_Destroy (bod); /* FIXME: bodies created in input files at time > 0;
					  FIXME: perhaps there is no need of moving GLV to the fist lng_RUN call,
					  FIXME: but rather bodies created in Python should not be put into the 'dom->newb' set;
					  FIXME: this way, as now, all Python created bodies will be anyway read at time 0 */
	      free (d);
	      free (i);
	    }
	    else
	    {
	      free (d);
	      free (i);
	      break;
	    }
	  }
	  else
	  {
	    free (d);
	    break;
	  }
	}
	else
	{
	  free (d);
	  break;
	}
      }
      else break;
    }

    xdr_destroy (&xdr);
    fclose (file);
  }

  free (path);
}

/* write domain state */
void dom_write_state (DOM *dom, PBF *bf)
{
  /* mark domain output */

  PBF_Label (bf, "DOM");

  /* write time step */

  PBF_Label (bf, "STEP");

  PBF_Double (bf, &dom->step, 1);

  /* write constraints merit */

  PBF_Label (bf, "MERIT");

  PBF_Double (bf, &dom->merit, 1);

  /* write complete data of newly created bodies and empty the newly created bodies set */

  write_new_bodies (dom); /* writing is done to a separate file */

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

  /* write constraints */

  PBF_Label (bf, "CONS");
 
  PBF_Int (bf, &dom->ncon, 1);

  for (CON *con = dom->con; con; con = con->next)
  {
    write_constraint (con, bf);
  }
}

/* read domain state */
void dom_read_state (DOM *dom, PBF *bf)
{
  BODY *bod, *next;
  int iover, ncon;

  /* input version */
  iover = dom->solfec->iover;

  /* clear contacts */
  MAP_Free (&dom->mapmem, &dom->idc);
  MEM_Release (&dom->conmem);
  dom->con = NULL;
  dom->ncon = 0;

  /* read all bodies if needed */
  if (!dom->allbodiesread) read_new_bodies (dom, bf);

  /* mark all bodies as absent */
  for (bod = dom->bod; bod; bod = bod->next) bod->flags |= BODY_ABSENT;

  for (; bf; bf = bf->next)
  {
    if (PBF_Label (bf, "DOM"))
    {
      /* read time step */

      ASSERT (PBF_Label (bf, "STEP"), ERR_FILE_FORMAT);

      PBF_Double (bf, &dom->step, 1);

      /* read constraints merit */

      ASSERT (PBF_Label (bf, "MERIT"), ERR_FILE_FORMAT);

      PBF_Double (bf, &dom->merit, 1);

      /* read body states */

      ASSERT (PBF_Label (bf, "BODS"), ERR_FILE_FORMAT);

      int nbod;

      PBF_Int (bf, &nbod, 1);

      for (int n = 0; n < nbod; n ++)
      {
	unsigned int id;

	PBF_Uint (bf, &id, 1);
	bod = MAP_Find (dom->idb, (void*) (long) id, NULL);

	if (bod == NULL) /* pick from all bodies set */
	{
	  ASSERT_DEBUG_EXT (bod = MAP_Find (dom->allbodies, (void*) (long) id, NULL), "Body id invalid");

	  if (bod->label) MAP_Insert (&dom->mapmem, &dom->lab, bod->label, bod, (MAP_Compare) strcmp);
	  MAP_Insert (&dom->mapmem, &dom->idb, (void*) (long) bod->id, bod, NULL);
	  bod->next = dom->bod;
	  if (dom->bod) dom->bod->prev = bod;
	  dom->bod = bod;
	  bod->dom = dom;
	  dom->nbod ++;
	}

	BODY_Read_State (bod, bf);
	bod->flags &= ~BODY_ABSENT;
      }

      /* read constraints */

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
    }
  }

  /* remove absent bodies */
  for (bod = dom->bod; bod; bod = next)
  {
    next = bod->next;

    if (bod->flags & BODY_ABSENT)
    {
      if (bod->label) MAP_Delete (&dom->mapmem, &dom->lab, bod->label, (MAP_Compare) strcmp);
      MAP_Delete (&dom->mapmem, &dom->idb, (void*) (long) bod->id, NULL);
      if (bod->next) bod->next->prev = bod->prev;
      if (bod->prev) bod->prev->next = bod->next;
      else dom->bod = bod->next;
      dom->nbod --;
    }
  }

  /* attach constraints to bodies */
  dom_attach_constraints (dom);
}

/* read state of an individual body */
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

/* read state of an individual constraint */
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


/* initialize domain state */
int dom_init_state (DOM *dom, PBF *bf)
{
  for (; bf; bf = bf->next)
  {
    if (PBF_Label (bf, "DOM"))
    {
      /* read body states */

      ASSERT (PBF_Label (bf, "BODS"), ERR_FILE_FORMAT);

      int nbod;

      PBF_Int (bf, &nbod, 1);

      for (int n = 0; n < nbod; n ++)
      {
	unsigned int id;
	BODY *bod;

	PBF_Uint (bf, &id, 1);
	ASSERT_TEXT (bod = MAP_Find (dom->allbodies, (void*) (long) id, NULL),
	             "Invalid body identifier => most likely due to a mismatched output file.");
	BODY_Read_State (bod, bf); /* XXX: we need to read all bodies since this can also be called
				           in parallel and only some bodies may be present; yet in order
					   to maintain the read consitency we need to read everything */
      }
    }
  }

  return 1;
}
