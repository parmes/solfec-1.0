/*
 * fra.c
 * Copyright (C) 2013, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * fracture code
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

#include <limits.h>
#include "ext/tetgen/tetsol.h"
#include "sol.h"
#include "dom.h"
#include "dio.h"
#include "fem.h"
#include "fra.h"
#include "mem.h"
#include "alg.h"
#include "err.h"
#include "fem.h"
#include "pck.h"
#include "lap.h"
#include "kdt.h"

/* free list */
void fracture_state_free (FS *list)
{
  FS *next;

  for (;list; list = next)
  {
    next = list->next;
    free (list->disp);
    free (list);
  }
}

/* write fracture state */
static void fracture_state_write (DOM *dom)
{
  char path [1024];
  double R[3], r, (*disp) [3];
  int i, n, dofs;
  MESH *msh;
  SET *item;
  BODY *bod;
  CON *con;

#if HDF5
  int numbod;
  PBF *f;

  snprintf (path, 1024, "%s/fracture", dom->solfec->outpath);
  ASSERT (f = PBF_Write (path, PBF_ON, PBF_ON), ERR_FILE_OPEN);

  PBF_Time (f, &dom->time);

  for (numbod = 0, bod = dom->bod; bod; bod = bod->next)
  {
    if (bod->fracture)
    {
      msh = bod->shape->data;
      dofs = 3 * msh->nodes_count;
      ERRMEM (disp = malloc (msh->nodes_count * sizeof (double [3])));
      for (i = 0; i < msh->nodes_count; i ++)
      {
        SUB (msh->cur_nodes [i], msh->ref_nodes [i], disp [i]);
      }

      PBF_Uint (f, &bod->id, 1);
      PBF_Int (f, &dofs, 1);
      PBF_Double (f, (double*)disp, dofs);

      n = SET_Size (bod->con);
      PBF_Int (f, &n, 1);
      for (item = SET_First (bod->con); item; item = SET_Next (item))
      {
	con = item->data;
	r = sqrt (con->area/ALG_PI);
	PBF_Double (f, &r, 1);

	if (bod == con->master)
	{
	  PBF_Double (f, con->mpnt, 3);
	}
	else
	{
	  PBF_Double (f, con->spnt, 3);
	}

        NVMUL (con->base, con->R, R);
	PBF_Double (f, R, 3);
      }

      bod->fracture = 0;

      free (disp);
      numbod ++;
    }

    PBF_Int2 (f, "numbod", &numbod, 1);
  }

  PBF_Close (f);
#else
  FILE *f;
  XDR x;

#if MPI
  snprintf (path, 1024, "%s/fracture%d.dat", dom->solfec->outpath, dom->rank);
#else
  snprintf (path, 1024, "%s/fracture.dat", dom->solfec->outpath);
#endif
  ASSERT (f = fopen (path, "a"), ERR_FILE_OPEN);
  xdrstdio_create (&x, f, XDR_ENCODE);

  for (bod = dom->bod; bod; bod = bod->next)
  {
    if (bod->fracture)
    {
      msh = bod->shape->data;
      dofs = 3 * msh->nodes_count;
      ERRMEM (disp = malloc (msh->nodes_count * sizeof (double [3])));
      for (i = 0; i < msh->nodes_count; i ++)
      {
        SUB (msh->cur_nodes [i], msh->ref_nodes [i], disp [i]);
      }

      ASSERT (xdr_u_int (&x, &bod->id), ERR_FILE_WRITE);
      ASSERT (xdr_int (&x, &dofs), ERR_FILE_WRITE);
      ASSERT (xdr_vector (&x, (char*)disp, dofs, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_WRITE);
      n = SET_Size (bod->con);
      ASSERT (xdr_int (&x, &n), ERR_FILE_WRITE);
      for (item = SET_First (bod->con); item; item = SET_Next (item))
      {
	con = item->data;
	r = sqrt (con->area/ALG_PI);
	ASSERT (xdr_double (&x, &r), ERR_FILE_WRITE);

	if (bod == con->master)
	{
          ASSERT (xdr_vector (&x, (char*)con->mpnt, 3, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_WRITE);
	}
	else
	{
          ASSERT (xdr_vector (&x, (char*)con->spnt, 3, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_WRITE);
	}

        NVMUL (con->base, con->R, R);
        ASSERT (xdr_vector (&x, (char*)R, 3, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_WRITE);
      }

      bod->fracture = 0;

      free (disp);
    }
  }

  xdr_destroy (&x);
  fclose (f);
#endif
}

/* read fracture state */
FS* fracture_state_read (BODY *bod)
{
  FS *out = NULL, *item, *instance;
  char path [1024];
  unsigned int id;
  int i, n, dofs;
  double *disp;

#if HDF5
  PBF *f, *g;

  snprintf (path, 1024, "%s/fracture", bod->dom->solfec->outpath);
  g = PBF_Read (path);

  do
  {
    double time;

    PBF_Time (g, &time); /* unused, but could be useful at some point */

    for (f = g; f; f = f->next)
    {
      int numbod;

      PBF_Int2 (f, "numbod", &numbod, 1);

      while (numbod > 0)
      {
	PBF_Uint (f, &id, 1);
	PBF_Int (f, &dofs, 1);
	ERRMEM (disp = malloc (dofs * sizeof (double)));
	PBF_Double (f, disp, dofs);
	PBF_Int (f, &n, 1);
	for (i = 0, instance = NULL; i < n; i ++)
	{
	  ERRMEM (item = MEM_CALLOC (sizeof (FS)));

	  PBF_Double (f, &item->radius, 1);
	  PBF_Double (f, item->point, 3);
	  PBF_Double (f, item->force, 3);

	  if (id == bod->id)
	  {
	    item->inext = instance;
	    instance = item;

	    if (i == (n-1))
	    {
	      item->disp = disp; /* put displacements into first element of instance list */
	      item->next = out;
	      out = item;
	    }
	  }
	  else free (item);
	}

	if (!out || out->disp != disp) free (disp);  /* not used */

	numbod --;
      }
    }
  } while (PBF_Forward (g, 1));

  PBF_Close (g);
#else
  FILE *f;
  XDR x;

  snprintf (path, 1024, "%s/fracture.dat", bod->dom->solfec->outpath);
  f = fopen (path, "r");
  /* TODO: read MPI mode data in case f == NULL but fractureRANK.dat exit */
  if (f)
  {
    xdrstdio_create (&x, f, XDR_DECODE);

    while (! feof (f))
    {
      if (xdr_u_int (&x, &id) == 0) break;
      ASSERT (xdr_int (&x, &dofs), ERR_FILE_READ);
      ERRMEM (disp = malloc (dofs * sizeof (double)));
      ASSERT (xdr_vector (&x, (char*)disp, dofs, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_READ);
      ASSERT (xdr_int (&x, &n), ERR_FILE_READ);
      for (i = 0, instance = NULL; i < n; i ++)
      {
        ERRMEM (item = MEM_CALLOC (sizeof (FS)));

	ASSERT (xdr_double (&x, &item->radius), ERR_FILE_READ);
        ASSERT (xdr_vector (&x, (char*)item->point, 3, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_READ);
        ASSERT (xdr_vector (&x, (char*)item->force, 3, sizeof (double), (xdrproc_t)xdr_double), ERR_FILE_READ);

	if (id == bod->id)
	{
	  item->inext = instance;
	  instance = item;

	  if (i == (n-1))
	  {
	    item->disp = disp; /* put displacements into first element of instance list */
	    item->next = out;
	    out = item;
	  }
	}
	else free (item);
      }

      if (!out || out->disp != disp) free (disp);  /* not used */
    }

    xdr_destroy (&x);
    fclose (f);
  }
#endif

  return out;
}

/* check fracture criterion */
void Fracture_Check (DOM *dom)
{
  short on = 0;
  BODY *bod;

#if 0
  for (bod = dom->bod; bod; bod = bod->next)
  {
    if (bod->flags & BODY_CHECK_FRACTURE)
    {
      double p [3], v [6], s [9], w [9];
      MESH *msh = FEM_MESH (bod);
      ELEMENT *ele;
      int bulk;

      VECTOR (p, 0.5, 0.5, 0.5);

      for (ele = msh->surfeles, bulk = 0; ele; )
      {
        BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);

        FEM_Element_Point_Values (bod, ele, p, VALUE_STRESS, v);

	s [0] = v [0];
	s [1] = v [3];
	s [2] = v [4];
	s [3] = s [1];
	s [4] = v [1];
	s [5] = v [5];
	s [6] = s [2];
	s [7] = s [5];
	s [8] = v [2];
        
	ASSERT (lapack_dsyev ('N', 'U', 3, s, 3, v, w, 9) == 0, ERR_LDY_EIGEN_DECOMP);

	if (v [2] > mat->tensile) /* maximal eigenvalue larger than tensile strength */
	{
	  bod->fracture = 1;
	  on = 1;
	}

	if (bulk) ele = ele->next;
	else if (ele->next) ele = ele->next;
	else ele = msh->bulkeles, bulk = 1;
      }
    }
  }
#else
  for (bod = dom->bod; bod; bod = bod->next)
  {
    if (bod->flags & BODY_CHECK_FRACTURE)
    {
      MESH *msh = FEM_MESH (bod);
      ELEMENT *ele;
      int bulk;
      double body_energy = 0.0;
      double frac_energy = 0.0;

      for (ele = msh->surfeles, bulk = 0; ele; )
      {
        BULK_MATERIAL *mat = FEM_MATERIAL (bod, ele);
	double volume;
        body_energy += FEM_Element_Internal_Energy (bod, msh, ele, &volume);
	frac_energy += mat->fracene * volume;

	  if (bulk) ele = ele->next;
	  else if (ele->next) ele = ele->next;
	  else ele = msh->bulkeles, bulk = 1;
      }

      // fracture condition
      if (body_energy >= frac_energy)
      {
        bod->fracture = 1;
        on = 1; 
      }
    }
  }
#endif

  if (on) fracture_state_write (dom);
}

/* export data for fracture analysis in Yaffems (return number of exported analysis instances) */
int Fracture_Export_Yaffems (BODY *bod, double volume, double quality, FILE *output)
{
  double extents [6], *q, *u, (*p) [3];
  SOLFEC *sol = bod->dom->solfec;
  int n, m, elno, fano;
  FS *list, *it, *jt;
  ELEMENT *ele;
  FACE *fac;
  MESH *msh;
  KDT *kd;

  if (!(bod->flags & BODY_CHECK_FRACTURE) || sol->mode == SOLFEC_WRITE) return 0;

  list = fracture_state_read (bod);

  if (list)
  {
    MESH *copy = MESH_Copy (bod->shape->data);
    MESH_Update (copy, NULL, NULL, NULL); /* reference configuration */
    msh = tetrahedralize1 (copy, volume, quality, -INT_MAX, -INT_MAX); /* generate tet mesh in reference configuration */
    MESH_Destroy (copy);

    /* allocate displacements on the tet mesh */
    ERRMEM (q = malloc (6 * msh->nodes_count * sizeof (double)));
    u = q + 3 * msh->nodes_count;

    /* map faces to a kd-tree for quick point queries */
    kd = KDT_Create (msh->nodes_count, (double*)msh->ref_nodes, 0.0);
    for (ele = msh->surfeles; ele; ele = ele->next)
    { 
      ELEMENT_Ref_Extents (msh, ele, extents);
      for (fac = ele->faces; fac; fac = fac->next) KDT_Drop (kd, extents, fac);
    }

    fprintf (output, "%s\n", "# vtk DataFile Version 2.0");
    fprintf (output, "%s\n", "Test Title");
    fprintf (output, "ASCII\n");
    fprintf (output, "\n");
    fprintf (output, "DATASET UNSTRUCTURED_GRID\n");
    fprintf (output, "POINTS %d float\n", msh->nodes_count);

    for (n = 0; n < msh->nodes_count; n ++)
    {
      fprintf (output, "%f %f %f\n", msh->ref_nodes [n][0], msh->ref_nodes [n][1], msh->ref_nodes [n][2]);
    }

    for (fano = 0, fac = msh->faces; fac; fano ++, fac = fac->n) fac->index = fano; /* count and index faces */

    ERRMEM (p = malloc (fano * sizeof (double [3]))); /* allocate face pressures */

    elno = msh->surfeles_count + msh->bulkeles_count;

    fprintf (output, "\n");
    fprintf (output, "CELLS %d %d\n", elno + fano, elno*5 + fano*4);

    for (ele = msh->surfeles; ele; ele = ele->next)
    {
      fprintf (output, "4 %d %d %d %d\n", ele->nodes[0], ele->nodes[1], ele->nodes[2], ele->nodes[3]);
    }

    for (ele = msh->bulkeles; ele; ele = ele->next)
    {
      fprintf (output, "4 %d %d %d %d\n", ele->nodes[0], ele->nodes[1], ele->nodes[2], ele->nodes[3]);
    }

    for (fac = msh->faces; fac; fac = fac->n)
    {
      fprintf (output, "3 %d  %d  %d\n", fac->nodes[0], fac->nodes[1], fac->nodes[2]);
    }

    fprintf (output, "\n");
    fprintf (output, "CELL_TYPES %d\n", elno + fano);

    for (n = 0; n < elno; n++)
    {
      fprintf (output, "10\n");
    }
    
    for (n = 0; n < fano; n++)
    {
      fprintf (output, "5\n");
    }

    fprintf (output, "\n");
    fprintf (output, "POINT_DATA %d\n", msh->nodes_count);
    for (it = list, m = 0; it; it = it->next, m ++)
    {
      /* map displacements from the hex to the tet mesh */
      FEM_Map_State (bod->shape->data, it->disp, bod->velo, msh, q, u); /* only it->disp to q mapping is used */

      fprintf (output, "\n");
      fprintf (output, "VECTORS disp%d float\n", m+1);

      for (n = 0; n < msh->nodes_count; n ++)
      {
	fprintf (output, "%f %f %f\n", q[3*n], q[3*n+1], q[3*n+2]);
      }
    }

    fprintf (output, "\n");
    fprintf (output, "CELL_DATA %d\n", elno + fano);
    for (it = list, m = 0; it; it = it->next, m ++)
    {
      fprintf (output, "\n");
      fprintf (output, "VECTORS pres%d float\n", m);

      for (n = 0; n < elno; n ++) /* skip elements */
      {
        fprintf (output, "0 0 0\n");
      }

      for (n = 0; n < fano; n ++)
      {
        SET (p [n], 0.0); /* zero face pressures */
      }

      for (jt = it; jt; jt = jt->inext) /* for each point force in this instance */
      {
        double (*ref) [3] = msh->ref_nodes;
	double a [3], b [3], c [3], area;
        SET *set = NULL, *item;
	double *qa, *qb, *qc;

        extents [0] = jt->point[0] - jt->radius - GEOMETRIC_EPSILON; /* set up search extents */
        extents [1] = jt->point[1] - jt->radius - GEOMETRIC_EPSILON;
        extents [2] = jt->point[2] - jt->radius - GEOMETRIC_EPSILON;
        extents [3] = jt->point[0] + jt->radius + GEOMETRIC_EPSILON;
        extents [4] = jt->point[1] + jt->radius + GEOMETRIC_EPSILON;
        extents [5] = jt->point[2] + jt->radius + GEOMETRIC_EPSILON;

	KDT_Pick_Extents (kd, extents, &set); /* pick kd-tree leaves within the extents */

	for (item = SET_First (set); item; item = SET_Next (item))
	{
	  KDT *leaf = item->data;
	  for (n = 0; n < leaf->n; n ++)
	  {
	    fac = leaf->data [n]; /* face dropped into this leaf */
            qa = &q[3*fac->nodes[0]];
            qb = &q[3*fac->nodes[1]];
            qc = &q[3*fac->nodes[2]];
	    ADD (ref[fac->nodes[0]], qa, a); /* current face nodes */
	    ADD (ref[fac->nodes[1]], qb, b);
	    ADD (ref[fac->nodes[2]], qc, c);
	    TRIANGLE_AREA (a, b, c, area); /* current face area */

	    if (area > 0.0) /* XXX */
	    {
	      p [fac->index][0] += jt->force [0] / area; /* add up pressure */
	      p [fac->index][1] += jt->force [1] / area;
	      p [fac->index][2] += jt->force [2] / area;
	    }
	  }
	}

	SET_Free (NULL, &set);
      }

      for (n = 0; n < fano; n ++)
      {
        fprintf (output, "%g %g %g\n", p[n][0], p[n][1], p[n][2]);
      }
    }

    fracture_state_free (list);

    free (q);

    free (p);

    return m;
  }

  return 0;
}

/* export data for fracture analysis in MoFEM (return number of exported analysis instances) */
int Fracture_Export_MoFEM (BODY *bod, double volume, double quality, FILE *output)
{
  double extents [6], *q, *u, (*p) [3];
  SOLFEC *sol = bod->dom->solfec;
  int n, elno, fano;
  FS *list, *it, *jt;
  ELEMENT *ele;
  FACE *fac;
  MESH *msh;
  KDT *kd;

  if (!(bod->flags & BODY_CHECK_FRACTURE) || sol->mode == SOLFEC_WRITE) return 0;

  list = fracture_state_read (bod);

  if (list)
  {
    MESH *copy = MESH_Copy (bod->shape->data);
    MESH_Update (copy, NULL, NULL, NULL); /* reference configuration */
    msh = tetrahedralize1 (copy, volume, quality, -INT_MAX, -INT_MAX); /* generate tet mesh in reference configuration */
    MESH_Destroy (copy);

    /* allocate displacements on the tet mesh */
    ERRMEM (q = malloc (6 * msh->nodes_count * sizeof (double)));
    u = q + 3 * msh->nodes_count;

    /* map faces to a kd-tree for quick point queries */
    kd = KDT_Create (msh->nodes_count, (double*)msh->ref_nodes, 0.0);
    for (ele = msh->surfeles; ele; ele = ele->next)
    { 
      ELEMENT_Ref_Extents (msh, ele, extents);
      for (fac = ele->faces; fac; fac = fac->next) KDT_Drop (kd, extents, fac);
    }

    //______________________________________________________
    // output file start
    for (fano = 0, fac = msh->faces; fac; fano ++, fac = fac->n) fac->index = fano; /* count and index faces */
    ERRMEM (p = malloc (fano * sizeof (double [3]))); /* allocate face pressures */

    elno = msh->surfeles_count + msh->bulkeles_count;

    fprintf (output, "mOFF %d %d %d\n", msh->nodes_count, fano, elno); // file header

    /* map displacements from the hex to the tet mesh */
    FEM_Map_State (bod->shape->data, bod->conf, bod->velo, msh, q, u); /* only bod->disp to q mapping is used */

    for (n = 0; n < msh->nodes_count; n ++)
    {
      fprintf (output, "%f %f %f %f %f %f\n", msh->ref_nodes [n][0], msh->ref_nodes [n][1], msh->ref_nodes [n][2], q[3*n], q[3*n+1], q[3*n+2]);
    }

    //______________________________________________________
    /* rewind the list to the end to find the last element,
       which corresponds to the earliest in time fracture instance */
    for (it = list; it->next; it = it->next) continue;

    /* for (it = list; it; it = it->next) */
    /* FIXME -- FIXME -- FIXME -- FIXME */
    {
      for (n = 0; n < fano; n ++)
      {
        SET (p [n], 0.0); /* zero face pressures */
      }

      for (jt = it; jt; jt = jt->inext) /* for each point force in this instance */
      {
        double (*ref) [3] = msh->ref_nodes;
        double a [3], b [3], c [3], area;
        SET *set = NULL, *item;
        double *qa, *qb, *qc;

        extents [0] = jt->point[0] - jt->radius - GEOMETRIC_EPSILON; /* set up search extents */
        extents [1] = jt->point[1] - jt->radius - GEOMETRIC_EPSILON;
        extents [2] = jt->point[2] - jt->radius - GEOMETRIC_EPSILON;
        extents [3] = jt->point[0] + jt->radius + GEOMETRIC_EPSILON;
        extents [4] = jt->point[1] + jt->radius + GEOMETRIC_EPSILON;
        extents [5] = jt->point[2] + jt->radius + GEOMETRIC_EPSILON;

        KDT_Pick_Extents (kd, extents, &set); /* pick kd-tree leaves within the extents */

	for (item = SET_First (set); item; item = SET_Next (item))
	{
	  KDT *leaf = item->data;
	  for (n = 0; n < leaf->n; n ++)
	  {
	    fac = leaf->data [n]; /* face dropped into this leaf */

	    qa = &q[3*fac->nodes[0]];
	    qb = &q[3*fac->nodes[1]];
	    qc = &q[3*fac->nodes[2]];

	    ADD (ref[fac->nodes[0]], qa, a); /* current face nodes */
	    ADD (ref[fac->nodes[1]], qb, b);
	    ADD (ref[fac->nodes[2]], qc, c);

	    TRIANGLE_AREA (a, b, c, area); /* current face area */

	    if (area > 0.0) /* XXX */
	    {
	      p [fac->index][0] += jt->force [0] / area; /* add up pressure */
	      p [fac->index][1] += jt->force [1] / area; /* FIXME: seems to be adding up to much pressure -> divided by area */
	      p [fac->index][2] += jt->force [2] / area;
	    }
	  }
	}

        SET_Free (NULL, &set);
      }

      for (fac = msh->faces, n=0; fac; fac = fac->n, n ++)
      {
        fprintf (output, "3 %d %d %d %g %g %g\n", fac->nodes[0], fac->nodes[1], fac->nodes[2], p[n][0], p[n][1], p[n][2]);
      }
    }

    //______________________________________________________
    for (ele = msh->surfeles; ele; ele = ele->next)
    {
      fprintf (output, "4 %d %d %d %d\n", ele->nodes[0], ele->nodes[1], ele->nodes[2], ele->nodes[3]);
    }

    for (ele = msh->bulkeles; ele; ele = ele->next)
    {
      fprintf (output, "4 %d %d %d %d\n", ele->nodes[0], ele->nodes[1], ele->nodes[2], ele->nodes[3]);
    }
    // output file complete
    //______________________________________________________

    fracture_state_free (list);
    free (q);
    free (p);
    MESH_Destroy (msh);
  }

  return 0;
}

/* Export body reference mesh into VTK */
int EXPORT_SOLFEC_MESH (BODY *bod, FILE *output)
{
  int n = 0, el_size = 0;
  ELEMENT *ele;

  MESH *msh = MESH_Copy (bod->shape->data);
  MESH_Update (msh, NULL, NULL, NULL); /* reference configuration */

  // file output start
  fprintf (output, "%s\n", "# vtk DataFile Version 2.0");
  fprintf (output, "%s\n", "SOLFEC ORIGINAL BODY MESH");
  fprintf (output, "ASCII\n");
  fprintf (output, "\n");
  fprintf (output, "DATASET UNSTRUCTURED_GRID\n");
  fprintf (output, "POINTS %d float\n", msh->nodes_count);

  for (; n < msh->nodes_count; n ++) 
  {
    fprintf (output, "%f %f %f\n", msh->ref_nodes [n][0], msh->ref_nodes [n][1], msh->ref_nodes [n][2]);
  }
  
  //______________________________________________________
  int elno = msh->surfeles_count + msh->bulkeles_count; /* element no */

  // calculate size
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    if (ele->type == 4) el_size += 5;
    else if (ele->type == 8) el_size += 9;
  }

  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    if (ele->type == 4) el_size += 5;
    else if (ele->type == 8) el_size += 9;
  }

  fprintf (output, "\n");
  fprintf (output, "CELLS %d %d\n", elno, el_size);

  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    if (ele->type == 4) fprintf (output, "4 %d %d %d %d\n", ele->nodes[0], ele->nodes[1], ele->nodes[2], ele->nodes[3]);
    else if (ele->type == 8) fprintf (output, "8 %d %d %d %d %d %d %d %d\n", ele->nodes[0], ele->nodes[1], ele->nodes[2], ele->nodes[3], ele->nodes[4], ele->nodes[5], ele->nodes[6], ele->nodes[7]);
  }

  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    if (ele->type == 4) fprintf (output, "4 %d %d %d %d\n", ele->nodes[0], ele->nodes[1], ele->nodes[2], ele->nodes[3]);
    else if (ele->type == 8) fprintf (output, "8 %d %d %d %d %d %d %d %d\n", ele->nodes[0], ele->nodes[1], ele->nodes[2], ele->nodes[3], ele->nodes[4], ele->nodes[5], ele->nodes[6], ele->nodes[7]);
  }

  //______________________________________________________
  fprintf (output, "\n");
  fprintf (output, "CELL_TYPES %d\n", elno);

  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    if (ele->type == 4) fprintf (output, "10\n");
    else if (ele->type == 8) fprintf (output, "12\n");
  }

  for (ele = msh->bulkeles; ele; ele = ele->next)
  {
    if (ele->type == 4) fprintf (output, "10\n");
    else if (ele->type == 8) fprintf (output, "12\n");
  }
  // file output complete

  MESH_Destroy (msh);

  return 0;
}
