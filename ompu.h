/* OpenMP utilities */

/*
The MIT License (MIT)

Copyright (c) 2017 EDF Energy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* Contributors: Tomasz Koziara */

#include "msh.h"
#include "dom.h"
#include "box.h"
#include "ldy.h"

inline static BODY** ompu_bodies (DOM *dom, int *n)
{
  int j;
  BODY **pbod, *bod;
  ERRMEM (pbod = malloc ((dom->nbod) * sizeof(BODY*)));
  for (bod = dom->bod, j = 0; bod; bod = bod->next, j++) pbod[j] = bod;
  *n = dom->nbod;
  return pbod;
}

inline static DIAB** ompu_diab (LOCDYN *ldy, int *n)
{
  int j = 0;
  DIAB *dia, **pdia;
  for (dia = ldy->dia; dia; dia = dia->n) j ++;
  *n = j;
  ERRMEM (pdia = malloc ((*n) * sizeof(DIAB*)));
  for (dia = ldy->dia, j = 0; dia; dia = dia->n, j++) pdia[j] = dia;
  return pdia;
}

inline static BOX** ompu_boxes (AABB *aabb, int *n)
{
  int j = 0;
  BOX *box, **pbox;
  for (box = aabb->lst; box; box = box->next) j ++;
  *n = j;
  ERRMEM (pbox = malloc ((*n) * sizeof(BOX*)));
  for (box = aabb->lst, j = 0; box; box= box->next, j++) pbox[j] = box;
  return pbox;
}

inline static CON** ompu_constraints (DOM *dom, int *n)
{
  int j;
  CON **pcon, *con;
  ERRMEM (pcon = malloc ((dom->ncon) * sizeof(CON*)));
  for (con = dom->con, j = 0; con; con = con->next, j++) pcon[j] = con;
  *n = dom->ncon;
  return pcon;
}

inline static FACE** ompu_faces (MESH *msh, int *n)
{
  int j = 0;
  FACE *fac, **pfac;
  for (fac = msh->faces; fac; fac = fac->n) j ++;
  *n = j;
  ERRMEM (pfac = malloc ((*n) * sizeof(FACE*)));
  for (fac = msh->faces, j = 0; fac; fac = fac->n, j++) pfac[j] = fac;
  return pfac;
}

inline static ELEMENT** ompu_elements (MESH *msh, int *n)
{
  int j;
  ELEMENT *ele, **pele;
  *n = msh->surfeles_count + msh->bulkeles_count;
  ERRMEM (pele = malloc ((*n) * sizeof(ELEMENT*)));
  for (ele = msh->surfeles, j = 0; ele; ele = ele->next, j++) pele[j] = ele;
  for (ele = msh->bulkeles; ele; ele = ele->next, j++) pele[j] = ele;
  return pele;
}

inline static ELEMENT** ompu_surfeles (MESH *msh, int *n)
{
  int j;
  ELEMENT *ele, **pele;
  *n = msh->surfeles_count;
  ERRMEM (pele = malloc ((*n) * sizeof(ELEMENT*)));
  for (ele = msh->surfeles, j = 0; ele; ele = ele->next, j++) pele[j] = ele;
  return pele;
}

inline static omp_lock_t* ompu_locks (int n)
{
  omp_lock_t *locks;
  ERRMEM (locks = malloc (n * sizeof(omp_lock_t)));
  for (int j = 0; j < n; j ++) omp_init_lock (&locks[j]);
  return locks;
}

inline static void ompu_locks_free (omp_lock_t *locks, int n)
{
  for (int j = 0; j < n; j ++) omp_destroy_lock (&locks[j]);
  free (locks);
}
