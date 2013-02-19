/*
 * psc.c
 * Copyright (C) 2013, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * parallel self-consistency tests
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

/* write body data to file before sending the body via MPI calles;
 * the file is SOLFEC->outpath/bodyID.data; */

#include <stdio.h>
#include <float.h>
#include "sol.h"
#include "dom.h"
#include "bod.h"
#include "psc.h"
#include "err.h"

#define DEQ(x,y) (fabs((x)-(y))<=DBL_EPSILON)

static void psc_write_face (FACE *fac, FILE *f)
{
  fwrite (fac->normal, sizeof (double), 3, f);
  fwrite (&fac->type, sizeof (int), 1, f);
  fwrite (fac->nodes, sizeof (int), fac->type, f);
  fwrite (&fac->index, sizeof (int), 1, f);
  fwrite (&fac->surface, sizeof (int), 1, f);

  /* XXX: skip fac->idata */
}

static void psc_write_element (ELEMENT *ele, FILE *f)
{
  FACE *fac;
  short i;

  fwrite (&ele->type, sizeof (short), 1, f);
  fwrite (&ele->neighs, sizeof (short), 1, f);
  fwrite (&ele->volume, sizeof (int), 1, f);
  fwrite (ele->nodes, sizeof (int), ele->type, f);

  /* XXX: skip ele->{mat, state} as they are not used in practice */
  /* XXX: skip ele->{domnum, dom} due to the same reasons */

  for (i = 0; i < ele->neighs; i ++) fwrite (&ele->adj[i]->flag, sizeof (short), 1, f);

  for (fac = ele->faces, i = 0; fac; fac = fac->next) i ++;

  fwrite (&i, sizeof (short), 1, f);

  for (fac = ele->faces; fac; fac = fac->next) psc_write_face (fac, f);
}

static void psc_write_mesh (MESH *msh, FILE *f)
{
  ELEMENT *ele;
  int i;

  fwrite (&msh->nodes_count, sizeof (int), 1, f);
  fwrite (msh->ref_nodes, sizeof (double [3]), msh->nodes_count, f);
  fwrite (msh->cur_nodes, sizeof (double [3]), msh->nodes_count, f);

  fwrite (&msh->surfeles_count, sizeof (int), 1, f);
  fwrite (&msh->bulkeles_count, sizeof (int), 1, f);

  /* number elements and store their index in ele->flag */
  for (ele = msh->surfeles, i = 0; ele; ele = ele->next, i ++) ele->flag = i;
  for (ele = msh->bulkeles; ele; ele = ele->next, i ++) ele->flag = i;

  for (ele = msh->surfeles; ele; ele = ele->next) psc_write_element (ele, f);
  for (ele = msh->bulkeles; ele; ele = ele->next) psc_write_element (ele, f);
}

static void psc_write_shape (SHAPE *shp, FILE *f)
{
  SHAPE *ptr;
  int n;

  for (ptr = shp, n = 0; ptr; ptr = ptr->next) n ++;

  ASSERT_TEXT (n == 1, "Compund shapes are not supported in PSC mode");
  ASSERT_TEXT (shp->kind == SHAPE_MESH, "Only MESH shapes are supported in PSC mode");
  /* TODO: support all shapes */

  psc_write_mesh (shp->data, f);
}

static void psc_write_matrix (MX *a, FILE *f)
{
  fwrite (&a->kind, sizeof (a->kind), 1, f);
  fwrite (&a->flags, sizeof (a->flags), 1, f);

  fwrite (&a->nzmax, sizeof (int), 1, f);
  fwrite (&a->m, sizeof (int), 1, f);
  fwrite (&a->n, sizeof (int), 1, f);
  fwrite (&a->nz, sizeof (int), 1, f);

  switch (a->kind)
  {
  case MXDENSE:
  {
    fwrite (a->x, sizeof (double), a->nzmax, f);
  }
  break;
  case MXBD:
  {
    fwrite (a->p, sizeof (int), a->n + 1, f);
    fwrite (a->i, sizeof (int), a->n + 1, f);
    fwrite (a->x, sizeof (double), a->nzmax, f);
  }
  break;
  case MXCSC:
  {
    fwrite (a->p, sizeof (int), a->n + 1, f);
    fwrite (a->i, sizeof (int), a->nzmax, f);
    fwrite (a->x, sizeof (double), a->nzmax, f);
  }
  break;
  }
}

static FACE* psc_read_face (MEM *facmem, FILE *f)
{
  FACE *fac;

  ERRMEM (fac = MEM_Alloc (facmem));

  fread (fac->normal, sizeof (double), 3, f);
  fread (&fac->type, sizeof (int), 1, f);
  fread (fac->nodes, sizeof (int), fac->type, f);
  fread (&fac->index, sizeof (int), 1, f);
  fread (&fac->surface, sizeof (int), 1, f);

  return fac;
}

static ELEMENT* psc_read_element (MEM *elemem, MEM *facmem, FILE *f)
{
  FACE *fac, *tail;
  ELEMENT *ele;
  short j;
  int i;

  ERRMEM (ele = MEM_Alloc (elemem));

  fread (&ele->type, sizeof (short), 1, f);
  fread (&ele->neighs, sizeof (short), 1, f);
  fread (&ele->volume, sizeof (int), 1, f);
  fread (ele->nodes, sizeof (int), ele->type, f);

  for (i = 0; i < ele->neighs; i ++)
  {
    fread (&j, sizeof (short), 1, f);
    ele->adj [i] = (void*) (long) j;
  }

  fread (&j, sizeof (short), 1, f);

  for (i = 0, tail = NULL; i < j; i ++)
  {
    fac = psc_read_face (facmem, f);
    fac->ele = ele;

    if (tail) tail->next = fac;
    else ele->faces = fac;
    tail = fac;
  }

  return ele;
}

static MESH* psc_read_mesh (FILE *f)
{
  ELEMENT *ele, **tab, *tail;
  MESH *msh;
  int i, j;

  ERRMEM (msh = MEM_CALLOC (sizeof (MESH)));
  MEM_Init (&msh->elemem, sizeof (ELEMENT), 128);
  MEM_Init (&msh->facmem, sizeof (FACE), 128);
  MEM_Init (&msh->mapmem, sizeof (MAP), 128);

  fread (&msh->nodes_count, sizeof (int), 1, f);

  ERRMEM (msh->ref_nodes = malloc (2 * msh->nodes_count * sizeof (double [3])));
  msh->cur_nodes = msh->ref_nodes + msh->nodes_count;

  fread (msh->ref_nodes, sizeof (double [3]), msh->nodes_count, f);
  fread (msh->cur_nodes, sizeof (double [3]), msh->nodes_count, f);

  fread (&msh->surfeles_count, sizeof (int), 1, f);
  fread (&msh->bulkeles_count, sizeof (int), 1, f);

  ERRMEM (tab = malloc ((msh->surfeles_count + msh->bulkeles_count) * sizeof (ELEMENT*)));

  for (i = 0, tail = NULL; i < msh->surfeles_count; i ++)
  {
    ele = psc_read_element (&msh->elemem, &msh->facmem, f);
    ele->flag = i;

    if (tail) ele->prev = tail, tail->next = ele;
    else msh->surfeles = ele;
    tail = ele;
    tab [i] = ele;
  }

  for (tail = NULL; i < msh->surfeles_count + msh->bulkeles_count; i ++)
  {
    ele = psc_read_element (&msh->elemem, &msh->facmem, f);
    ele->flag = i;

    if (tail) ele->prev = tail, tail->next = ele;
    else msh->bulkeles = ele;
    tail = ele;
    tab [i] = ele;
  }

  for (i = 0; i < msh->surfeles_count + msh->bulkeles_count; i ++)
  {
    ele = tab [i];

    for (j = 0; j < ele->neighs; j ++)
    {
      int idx = (long) (void*) ele->adj[j];
      ASSERT_TEXT (idx >= 0 && idx < msh->surfeles_count + msh->bulkeles_count, "INCONSITENT ELEMENT INDEXING");
      ele->adj[j] = tab [idx];
    }
  }

  free (tab);

  return msh;
}

static SHAPE* psc_read_shape (FILE *f)
{
  SHAPE *shp;

  ERRMEM (shp = MEM_CALLOC (sizeof (SHAPE)));

  shp->kind = SHAPE_MESH;
  shp->data = psc_read_mesh (f);

  return shp;
}

static int psc_compare_faces (FACE *a, FACE *b)
{
  int i;

  if (a->normal [0] != b->normal [0] ||
      a->normal [1] != b->normal [1] ||
      a->normal [2] != b->normal [2]) return 0;

  if (a->type != b->type) return 0;

  for (i = 0; i < a->type; i ++)
  {
    if (a->nodes [i] != b->nodes [i]) return 0;
  }

  if (a->index != b->index) return 0;

  if (a->surface != b->surface) return 0;

  if (a->ele->flag != b->ele->flag) return 0;

  return 1;
}

static int psc_compare_elements (ELEMENT *a, ELEMENT *b)
{
  FACE *fa, *fb;
  int i;

  if (a->type != b->type)
  {
    fprintf (stderr, "\nPSC: ELEMENT => type");
    return 0;
  }

  if (a->neighs != b->neighs)
  {
    fprintf (stderr, "\nPSC: ELEMENT => neighbours count");
    return 0;
  }

  if (a->volume != b->volume)
  {
    fprintf (stderr, "\nPSC: ELEMENT => volume id");
    return 0;
  }

  for (i = 0; i < a->type; i ++)
  {
    if (a->nodes [i] != b->nodes [i])
    {
      fprintf (stderr, "\nPSC: ELEMENT => node index");
      return 0;
    }
  }

  for (i = 0; i < a->neighs; i ++)
  {
    if (a->adj[i]->flag != b->adj[i]->flag)
    {
      fprintf (stderr, "\nPSC: ELEMENT => neighbour index");
      return 0;
    }
  }

  for (fa = a->faces, fb = b->faces; fa || fb; fa = fa->next, fb = fb->next)
  {
    if (fa && !fb)
    {
      fprintf (stderr, "\nPSC: ELEMENT => number of faces");
      return 0;
    }
    if (!fa && fb)
    {
      fprintf (stderr, "\nPSC: ELEMENT => number of faces");
      return 0;
    }

    if (psc_compare_faces (fa, fb) == 0)
    {
      fprintf (stderr, "\nPSC: ELEMENT => face");
      return 0;
    }
  }

  return 1;
}

static int psc_compare_meshes (MESH *a, MESH *b)
{
  ELEMENT *ele, *y;
  int i, j;

  if (a->nodes_count != b->nodes_count)
  {
    fprintf (stderr, "\nPSC: MESH => nodes_count\n");
    return 0;
  }

  for (i = 0; i < a->nodes_count; i ++)
  {
    for (j = 0; j < 3; j ++)
    {
      if (a->ref_nodes[i][j] != b->ref_nodes[i][j])
      {
        fprintf (stderr, "\nPSC: MESH => ref_nodes\n");
        return 0;
      }
    }
  }

  for (i = 0; i < a->nodes_count; i ++)
  {
    for (j = 0; j < 3; j ++)
    {
      if (a->cur_nodes[i][j] != b->cur_nodes[i][j])
      {
        fprintf (stderr, "\nPSC: MESH => cur_nodes\n");
	return 0;
      }
    }
  }

  if (a->surfeles_count != b->surfeles_count) return 0;

  if (a->bulkeles_count != b->bulkeles_count) return 0;

  /* number elements and store their index in ele->flag */
  for (ele = a->surfeles, i = 0; ele; ele = ele->next, i ++) ele->flag = i;
  for (ele = a->bulkeles; ele; ele = ele->next, i ++) ele->flag = i;
  for (ele = b->surfeles, i = 0; ele; ele = ele->next, i ++) ele->flag = i;
  for (ele = b->bulkeles; ele; ele = ele->next, i ++) ele->flag = i;

  for (ele = a->surfeles, y = b->surfeles; ele; ele = ele->next, y = y->next)
  {
    if (psc_compare_elements (ele, y) == 0)
    {
      fprintf (stderr, "\nPSC: MESH => surface element\n");
      return 0;
    }
  }

  for (ele = a->bulkeles, y = b->bulkeles; ele; ele = ele->next, y = y->next)
  {
    if (psc_compare_elements (ele, y) == 0)
    {
      fprintf (stderr, "\nPSC: MESH => bulk element\n");
      return 0;
    }
  }

  return 1;
}

static int psc_compare_shapes (SHAPE *a, SHAPE *b)
{
  if (a->kind != b->kind) return 0;

  return psc_compare_meshes (a->data, b->data);
}

static MX* psc_read_matrix (FILE *f)
{
  MX *a;

  ERRMEM (a = MEM_CALLOC (sizeof (MX)));

  fread (&a->kind, sizeof (a->kind), 1, f);
  fread (&a->flags, sizeof (a->flags), 1, f);

  fread (&a->nzmax, sizeof (int), 1, f);
  fread (&a->m, sizeof (int), 1, f);
  fread (&a->n, sizeof (int), 1, f);
  fread (&a->nz, sizeof (int), 1, f);

  switch (a->kind)
  {
  case MXDENSE:
  {
    ERRMEM (a->x = malloc (a->nzmax * sizeof (double)));
    fread (a->x, sizeof (double), a->nzmax, f);
  }
  break;
  case MXBD:
  {
    ERRMEM (a->p = malloc ((a->n+1) * sizeof (int)));
    ERRMEM (a->i = malloc ((a->n+1) * sizeof (int)));
    ERRMEM (a->x = malloc (a->nzmax * sizeof (double)));
    fread (a->p, sizeof (int), a->n + 1, f);
    fread (a->i, sizeof (int), a->n + 1, f);
    fread (a->x, sizeof (double), a->nzmax, f);
  }
  break;
  case MXCSC:
  {
    ERRMEM (a->p = malloc ((a->n+1) * sizeof (int)));
    ERRMEM (a->i = malloc (a->nzmax * sizeof (int)));
    ERRMEM (a->x = malloc (a->nzmax * sizeof (double)));
    fread (a->p, sizeof (int), a->n + 1, f);
    fread (a->i, sizeof (int), a->nzmax, f);
    fread (a->x, sizeof (double), a->nzmax, f);
  }
  break;
  }

  return a;
}

static int psc_compare_matrices (MX *a, MX *b)
{
  int j;

  if (a->kind != b->kind) return 0;
  if (a->flags != b->flags) return 0;

  if (a->nzmax != b->nzmax) return 0;
  if (a->m != b->m) return 0;
  if (a->n != b->n) return 0;
  if (a->nz != b->nz) return 0;

  switch (a->kind)
  {
  case MXDENSE:
  {
    for (j = 0; j < a->nzmax; j ++)
    {
      if (a->x[j] != b->x[j]) return 0;
    }
  }
  break;
  case MXBD:
  {
    for (j = 0; j < a->n + 1; j ++)
    {
      if (a->p[j] != b->p[j]) return 0;
    }
    for (j = 0; j < a->n + 1; j ++)
    {
      if (a->i[j] != b->i[j]) return 0;
    }
    for (j = 0; j < a->nzmax; j ++)
    {
      if (a->x[j] != b->x[j]) return 0;
    }
  }
  break;
  case MXCSC:
  {
    for (j = 0; j < a->n + 1; j ++)
    {
      if (a->p[j] != b->p[j]) return 0;
    }
    for (j = 0; j < a->nzmax; j ++)
    {
      if (a->i[j] != b->i[j]) return 0;
    }
    for (j = 0; j < a->nzmax; j ++)
    {
      if (a->x[j] != b->x[j]) return 0;
    }
  }
  break;
  }

  return 0;
}

static void psc_matrix_free (MX *a)
{
  switch (a->kind)
  {
  case MXDENSE:
  {
    free (a->x);
  }
  break;
  case MXBD:
  case MXCSC:
  {
    free (a->p);
    free (a->i);
    free (a->x);

  }
  break;
  }

  free (a);
}

/* write body data to file before sending the body via MPI calles;
 * the file is SOLFEC->outpath/bodyID.data; */
void PSC_Write_Body (BODY *bod)
{
  char txt [1024];
  FILE *f;
  int i;

  snprintf (txt, 1024, "%s/body%d.data", bod->dom->solfec->outpath, bod->id);
  ASSERT (f = fopen (txt, "w"), ERR_FILE_OPEN);

  fwrite (&bod->kind, sizeof (bod->kind), 1, f);

#if 0
  i = strlen(bod->mat->label);
  ASSERT_TEXT (i < 1024, "Material label is too long!");
  fwrite (&i, sizeof (int), i, f);
  fwrite (bod->mat->label, sizeof (char), i, f);
#endif

  fwrite (&bod->ref_mass, sizeof (double), 1, f);
  fwrite (&bod->ref_volume, sizeof (double), 1, f);
  fwrite (bod->ref_center, sizeof (double), 3, f);
  fwrite (bod->ref_tensor, sizeof (double), 9, f);

  fwrite (&bod->dofs, sizeof (int), 1, f);

  fwrite (&bod->form, sizeof (bod->form), 1, f);

  int confsize = bod->kind != FEM ? 12 : bod->form == REDUCED_ORDER ? bod->dofs + 9 : bod->dofs;

  fwrite (bod->conf, sizeof (double), confsize, f);
  fwrite (bod->velo, sizeof (double), bod->dofs, f);

  /* TODO => bod->forces */

  /* TODO => bod->cra */

  psc_write_shape (bod->shape, f);

  fwrite (bod->extents, sizeof (double), 6, f);

  fwrite (&bod->scheme, sizeof (bod->scheme), 1, f);

  psc_write_matrix (bod->inverse, f);

  psc_write_matrix (bod->M, f);

  psc_write_matrix (bod->K, f);

  fwrite (&bod->damping, sizeof (double), 1, f);

  if (bod->evec)
  {
    i = 1;

    fwrite (&i, sizeof (int), 1, f);

    psc_write_matrix (bod->evec, f);

    fwrite (bod->eval, sizeof (double), bod->evec->n, f);
  }
  else
  {
    i = 0;

    fwrite (&i, sizeof (int), 1, f);
  }

  /* XXX: skip bod->label as non-essential */

  /* XXX: skip bod->mesh for the moment */

  /* XXX: skip bod->energy as non-essential */

  /* XXX: skip bod->fracture as non-essential */

  fclose (f);
}

/* after the body has been imported via MPI calls,
 * read body data from file and compare */
void PSC_Test_Body (BODY *bod)
{
  char txt [1024];
  FILE *f;
  BODY a;
  int i;

  snprintf (txt, 1024, "%s/body%d.data", bod->dom->solfec->outpath, bod->id);
  ASSERT (f = fopen (txt, "r"), ERR_FILE_OPEN);

  fread (&a.kind, sizeof (a.kind), 1, f);

  if (a.kind != bod->kind)
  {
    ASSERT_TEXT (0, "PSC ERROR: kind");
  }

#if 0
  fread (&i, sizeof (int), 1, f);
  fread (txt, sizeof (char), i, f);
  txt[i] = '\0';

  if (strcmp (txt, bod->mat->label) != 0)
  {
    ASSERT_TEXT (0, "PSC ERROR: material");
  }
#endif

  fread (&a.ref_mass, sizeof (double), 1, f);

  if (a.ref_mass != bod->ref_mass)
  {
    ASSERT_TEXT (0, "PSC ERROR: ref_mass");
  }

  fread (&a.ref_volume, sizeof (double), 1, f);

  if (a.ref_volume != bod->ref_volume)
  {
    ASSERT_TEXT (0, "PSC ERROR: ref_volume");
  }

  fread (&a.ref_center, sizeof (double), 3, f);

  for (i = 0; i < 3; i ++)
  {
    if (a.ref_center [i]!= bod->ref_center [i])
    {
      ASSERT_TEXT (0, "PSC ERROR: ref_center");
    }
  }

  fread (&a.ref_tensor, sizeof (double), 9, f);

  for (i = 0; i < 9; i ++)
  {
    if (a.ref_tensor [i]!= bod->ref_tensor [i])
    {
      ASSERT_TEXT (0, "PSC ERROR: ref_tensor");
    }
  }

  fread (&a.dofs, sizeof (int), 1, f);

  if (a.dofs != bod->dofs)
  {
    ASSERT_TEXT (0, "PSC ERROR: dofs");
  }

  fread (&a.form, sizeof (a.form), 1, f);

  if (a.form != bod->form)
  {
    ASSERT_TEXT (0, "PSC ERROR: form");
  }

  int confsize = bod->kind != FEM ? 12 : bod->form == REDUCED_ORDER ? bod->dofs + 9 : bod->dofs;

  ERRMEM (a.conf = malloc (sizeof (double [confsize])));

  fread (a.conf, sizeof (double), confsize, f);

  for (i = 0; i < confsize; i ++)
  {
    if (a.conf [i] != bod->conf [i])
    {
      ASSERT_TEXT (0, "PSC ERROR: conf");
    }
  }

  free (a.conf);

  ERRMEM (a.velo = malloc (sizeof (double [a.dofs])));

  fread (a.velo, sizeof (double), a.dofs, f);

  for (i = 0; i < bod->dofs; i ++)
  {
    if (!DEQ(a.velo[i], bod->velo[i])) /* XXX: differs after 15th decimal place => why is velocity giving this kind of trouble? */
    {
      double x = fabs (a.velo[i]-bod->velo[i]),
             y = DBL_EPSILON;
      printf ("%.17f > %.17f\n", x, y);
      ASSERT_TEXT (0, "PSC ERROR: velo => %d => %.15f != %.15f", i, a.velo[i], bod->velo[i]);
    }
  }

  free (a.velo);

  /* TODO => bod->forces */

  /* TODO => bod->cra */

  a.shape = psc_read_shape (f);

  if (psc_compare_shapes (a.shape, bod->shape) == 0)
  {
    ASSERT_TEXT (0, "PSC ERROR: shape");
  }

  SHAPE_Destroy (a.shape);

  fread (a.extents, sizeof (double), 6, f);

  for (i = 0; i < 6; i ++)
  {
    if (a.extents [i] != bod->extents [i])
    {
      ASSERT_TEXT (0, "PSC ERROR: extents => %d => %.15f != %.15f", i, a.extents[i], bod->extents[i]);
    }
  }

  fread (&a.scheme, sizeof (a.scheme), 1, f);

  if (a.scheme != bod->scheme)
  {
    ASSERT_TEXT (0, "PSC ERROR: scheme");
  }

  a.inverse = psc_read_matrix (f);

  if (psc_compare_matrices (a.inverse, bod->inverse) == 0)
  {
    ASSERT_TEXT (0, "PSC ERROR: inverse");
  }

  psc_matrix_free (a.inverse);

  a.M = psc_read_matrix (f);

  if (psc_compare_matrices (a.M, bod->M) == 0)
  {
    ASSERT_TEXT (0, "PSC ERROR: M");
  }

  psc_matrix_free (a.M);

  a.K = psc_read_matrix (f);

  if (psc_compare_matrices (a.K, bod->K) == 0)
  {
    ASSERT_TEXT (0, "PSC ERROR: K");
  }

  psc_matrix_free (a.K);

  fread (&a.damping, sizeof (double), 1, f);

  if (a.damping != bod->damping)
  {
    ASSERT_TEXT (0, "PSC ERROR: damping");
  }

  fwrite (&i, sizeof (int), 1, f);

  if (i && bod->evec == NULL)
  {
    ASSERT_TEXT (0, "PSC ERROR: evec existence");
  }

  if (i)
  {
    a.evec = psc_read_matrix (f);

    if (psc_compare_matrices (a.evec, bod->evec) == 0)
    {
      ASSERT_TEXT (0, "PSC ERROR: evec");
    }

    ERRMEM (a.eval = malloc (sizeof (double [a.evec->n])));

    fread (a.eval, sizeof (double), a.evec->n, f);

    for (i = 0; i < a.evec->m; i ++)
    {
      if (a.eval [i]!= bod->eval [i])
      {
	ASSERT_TEXT (0, "PSC ERROR: eval");
      }
    }

    psc_matrix_free (a.evec);

    free (a.eval);
  }

  /* XXX: skip bod->label as non-essential */

  /* XXX: skip bod->mesh for the moment */

  /* XXX: skip bod->energy as non-essential */

  /* XXX: skip bod->fracture as non-essential */

  fclose (f);
}
