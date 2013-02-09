/*
 * cvitest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of convex intersection
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

#if OPENGL
#if __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif
#include "glv.h"

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdio.h>
#include <time.h>
#include "alg.h"
#include "hul.h"
#include "cvi.h"
#include "err.h"

#define minim 4
#define limit 128
int counter = 0;
FILE *input = NULL;
short GLVON = 0;
enum  {GEN, HUL, CVI} mode = HUL;
double apoint [limit][3],
       bpoint [limit][3],
       aplane [limit][6],
       bplane [limit][6];
int DOGEN = 1; /* set 0 in case of a file input */
CVIKIND kind = NON_REGULARIZED;
int asize, bsize,
    anpla, bnpla;
TRI *a = NULL,
    *b = NULL,
    *c = NULL;
int alength,
    blength,
    clength;

static void cube (double *centre, double edge, double (*point) [3])
{
  double half [3] = {0.5, 0.5, 0.5};

  point [0][0] = 0;
  point [0][1] = 0;
  point [0][2] = 0;

  point [1][0] = 0;
  point [1][1] = 1;
  point [1][2] = 0;
  
  point [2][0] = 1;
  point [2][1] = 1;
  point [2][2] = 0;
  
  point [3][0] = 1;
  point [3][1] = 0;
  point [3][2] = 0;
  
  point [4][0] = 0;
  point [4][1] = 0;
  point [4][2] = 1;
  
  point [5][0] = 0;
  point [5][1] = 1;
  point [5][2] = 1;
  
  point [6][0] = 1;
  point [6][1] = 1;
  point [6][2] = 1;
  
  point [7][0] = 1;
  point [7][1] = 0;
  point [7][2] = 1;

  for (int n = 0; n < 8; n ++)
  {
    SCALE (point [n], edge);
    SUB (point [n], half, point [n]);
    ADD (point [n], centre, point [n]);
  }
}

static void gen ()
{
  double move [3];

  if (rand () % 2 == 0) /* cube */
  {
    double a [3], b [3], r [9];
    int i;
   
    asize = 8;
    bsize = 8;
    
    SETRAND (a, 0.5);
    SETRAND (b, 0.5);

    cube (a, DRANDEXT (0.5, 1.5), apoint);
    cube (b, DRANDEXT (0.5, 1.5), bpoint);

    SETRAND (a, 1.0);
    SETRAND (b, 1.0);

    EXPMAP (a, r);
    for (i = 0; i < 8; i ++)
    { COPY (apoint [i], a);
      NVMUL (r, a, apoint [i]); }

    EXPMAP (b, r);
    for (i = 0; i < 8; i ++)
    { COPY (bpoint [i], b);
      NVMUL (r, b, bpoint [i]); }
  }
  else /* points */
  {
    asize = bsize = 0;

    while (asize < minim) asize = rand () % limit;

    while (bsize < minim) bsize = rand () % limit;

    SETRAND (move, 0.05);

    for (int n = 0; n < asize; n ++)
    { SETRAND (apoint [n], 1.0);
      ADD (apoint [n], move, apoint [n]); }

    for (int n = 0; n < bsize; n ++)
    { SETRAND (bpoint [n], 1.0);
      SUB (bpoint [n], move, bpoint [n]); }
  }
}

static void export (void)
{
  int i;

  printf ("%.24e\n", GEOMETRIC_EPSILON);
  printf ("%d   %d\n", asize, anpla);
  for (i = 0; i < asize; i ++) printf ("%.24e   %.24e   %.24e\n", apoint[i][0], apoint[i][1], apoint[i][2]);
  for (i = 0; i < anpla; i ++) printf ("%.24e   %.24e   %.24e   %.24e   %.24e   %.24e\n", aplane[i][0], aplane[i][1], aplane[i][2], aplane[i][3], aplane[i][4], aplane[i][5]);
  printf ("%d   %d\n", bsize, bnpla);
  for (i = 0; i < bsize; i ++) printf ("%.24e   %.24e   %.24e\n", bpoint[i][0], bpoint[i][1], bpoint[i][2]);
  for (i = 0; i < bnpla; i ++) printf ("%.24e   %.24e   %.24e   %.24e   %.24e   %.24e\n", bplane[i][0], bplane[i][1], bplane[i][2], bplane[i][3], bplane[i][4], bplane[i][5]);
}

static void inp (void)
{
  double extents [6], d [3];
  int i, j;

  printf ("COUNTER: %d\n", counter ++);

  if (input && !feof (input))
  {
    if (fscanf (input, "%lf", &GEOMETRIC_EPSILON) == EOF) goto LEOF;

    if (fscanf (input, "%d", &asize) == EOF) goto LEOF;
    if (fscanf (input, "%d", &anpla) == EOF) goto LEOF;
    if (asize > limit || anpla > limit) { fprintf (stderr, "Increase 'limit' constant and recompile 'cvitest'\n"); exit (0); }
    for (i = 0; i < asize; i ++) fscanf (input, "%lf", &apoint [i][0]), fscanf (input, "%lf", &apoint [i][1]), fscanf (input, "%lf", &apoint [i][2]);
    for (i = 0; i < anpla; i ++) fscanf (input, "%lf", &aplane [i][0]), fscanf (input, "%lf", &aplane [i][1]), fscanf (input, "%lf", &aplane [i][2]),
                                 fscanf (input, "%lf", &aplane [i][3]), fscanf (input, "%lf", &aplane [i][4]), fscanf (input, "%lf", &aplane [i][5]);

    if (fscanf (input, "%d", &bsize) == EOF) goto LEOF;
    if (fscanf (input, "%d", &bnpla) == EOF) goto LEOF;
    if (bsize > limit || bnpla > limit) { fprintf (stderr, "Increase 'limit' constant and recompile 'cvitest'\n"); exit (0); }
    for (i = 0; i < bsize; i ++) fscanf (input, "%lf", &bpoint [i][0]), fscanf (input, "%lf", &bpoint [i][1]), fscanf (input, "%lf", &bpoint [i][2]);
    for (i = 0; i < bnpla; i ++) fscanf (input, "%lf", &bplane [i][0]), fscanf (input, "%lf", &bplane [i][1]), fscanf (input, "%lf", &bplane [i][2]),
                                 fscanf (input, "%lf", &bplane [i][3]), fscanf (input, "%lf", &bplane [i][4]), fscanf (input, "%lf", &bplane [i][5]);
  }
  else
  {
LEOF:
    printf ("END OF FILE\n");
  }

  if (GLVON)
  {
    extents [0] = extents [1] = extents [2] =  DBL_MAX;
    extents [3] = extents [4] = extents [5] = -DBL_MAX;

    for (i = 0; i < asize; i ++)
    {
      for (j = 0; j < 3; j ++)
      {
	if (apoint [i][j] < extents [j]) extents [j] = apoint [i][j];
	if (apoint [i][j] > extents [3+j]) extents [3+j] = apoint [i][j];
      }
    }

    for (i = 0; i < bsize; i ++)
    {
      for (j = 0; j < 3; j ++)
      {
	if (bpoint [i][j] < extents [j]) extents [j] = bpoint [i][j];
	if (bpoint [i][j] > extents [3+j]) extents [3+j] = bpoint [i][j];
      }
    }

    for (i = 0; i < 3; i ++) d [i] = extents [3+i] - extents [i];

    for (i = 0; i < 3; i ++)
    {
      extents [3+i] += 0.25 * d [i];
      extents [i] -= 0.25 * d [i];
    }

    GLV_Reset_Extents (extents);
  }
}

static void render (void)
{
  glPointSize (3.0);

  if (mode == CVI)
  {
    glEnable (GL_NORMALIZE);

    glColor3d (1.0, 0.0, 0.0);
    glBegin (GL_TRIANGLES);
    for (TRI *t = a, *e = t + alength; t < e; t ++)
    {
      glNormal3dv (t->out);
      glVertex3dv (t->ver [0]);
      glVertex3dv (t->ver [1]);
      glVertex3dv (t->ver [2]);
    }
    glEnd ();

    glColor3d (0.0, 1.0, 0.0);
    glBegin (GL_TRIANGLES);
    for (TRI *t = b, *e = t + blength; t < e; t ++)
    {
      glNormal3dv (t->out);
      glVertex3dv (t->ver [0]);
      glVertex3dv (t->ver [1]);
      glVertex3dv (t->ver [2]);
    }
    glEnd ();
  }
  else if (mode == GEN)
  {
    glEnable (GL_NORMALIZE);

    glBegin (GL_TRIANGLES);
    for (TRI *t = c, *e = t + clength; t < e; t ++)
    {
      if (t->flg > 0) glColor3d (1.0, 0.0, 0.0);
      else glColor3d (0.0, 1.0, 0.0);
      glNormal3dv (t->out);
      glVertex3dv (t->ver [0]);
      glVertex3dv (t->ver [1]);
      glVertex3dv (t->ver [2]);
    }
    glEnd ();
  }
  else if (mode == HUL)
  {
    glDisable (GL_LIGHTING);
    glColor3d (1.0, 0.5, 0.5);
    glBegin (GL_POINTS);
    for (int n = 0; n < asize; n ++)
      glVertex3dv (apoint [n]);
    glEnd ();

    glColor3d (0.5, 1.0, 0.5);
    glBegin (GL_POINTS);
    for (int n = 0; n < bsize; n ++)
      glVertex3dv (bpoint [n]);
    glEnd ();
    glEnable (GL_LIGHTING);

    if (!DOGEN)
    {
      double p [3], eps;

      eps = GLV_Minimal_Extent () * 0.1;
      glLineWidth (2.0);
      glDisable (GL_LIGHTING);
      glColor3d (1.0, 0.5, 0.5);
      glBegin (GL_LINES);
      for (int n = 0; n < anpla; n ++)
      {
	ADDMUL (aplane[n]+3, eps, aplane[n], p);
	glVertex3dv (aplane[n]+3);
	glVertex3dv (p);
      }
      glEnd ();

      glColor3d (0.5, 1.0, 0.5);
      glBegin (GL_LINES);
      for (int n = 0; n < bnpla; n ++)
      {
	ADDMUL (bplane[n]+3, eps, bplane[n], p);
	glVertex3dv (bplane[n]+3);
	glVertex3dv (p);
      }
      glEnd ();
      glEnable (GL_LIGHTING);
      glLineWidth (1.0);
    }
  }

  GLVON = 1;
}

static void key (int key, int x, int y)
{
  if (key == ' ')
  {

    switch (mode)
    {
      case GEN:
	if (DOGEN) gen ();
	else inp ();
	mode = HUL;
      break;
      case HUL:
        free (a);
        free (b);
        a = hull ((double*)apoint, asize, &alength);
        b = hull ((double*)bpoint, bsize, &blength);
	mode = CVI;
      break;
      case CVI:
      {
	double *va, *pa,
	       *vb, *pb;
	int nva, npa,
	    nvb, npb;

	if (DOGEN)
	{
	  va = TRI_Vertices (a, alength, &nva);
	  pa = TRI_Planes (a, alength, &npa);
	  vb = TRI_Vertices (b, blength, &nvb);
	  pb = TRI_Planes (b, blength, &npb);

	  if (nva < limit && nvb < limit &&
	      npa < limit && npb < limit)
	  {
	    asize = nva;
	    anpla = npa;
	    bsize = nvb;
	    bnpla = npb;
	    memcpy (apoint, va, nva * sizeof (double [3]));
	    memcpy (aplane, pa, npa * sizeof (double [6]));
	    memcpy (bpoint, vb, nvb * sizeof (double [3]));
	    memcpy (bplane, pb, npb * sizeof (double [6]));
	  }
	  else
	  {
	    asize = anpla =
	    bsize = bnpla = 0;
	  }
	}
	else
	{
	  va = (double*)apoint;
	  pa = (double*)aplane;
	  vb = (double*)bpoint;
	  pb = (double*)bplane;
	  nva = asize;
	  npa = anpla;
	  nvb = bsize;
	  npb = bnpla;
	}

        free (c);
	c = cvi (va, nva, pa, npa,
	         vb, nvb, pb, npb,
		 kind, &clength);
	mode = GEN;

	if (DOGEN)
	{
	  free (va);
	  free (pa);
	  free (vb);
	  free (pb);
	}
      }
      break;
    }

    GLV_Redraw_All ();
  }
  else if (key == 'e') export ();
  else if (key == 'k')
  {
    if (kind == REGULARIZED) kind = NON_REGULARIZED, printf ("NON_REGULARIZED\n");
    else kind = REGULARIZED, printf ("REGULARIZED\n");
  }
}

int main (int argc, char **argv)
{
  double extents [6], d [3];
  int i, j;
  
  printf ("SPACE - iterate over the test stages\n");
  printf ("e - export current case\n");
  printf ("k - switch between regularized and non-regularized intersection kind\n");

  srand ((unsigned) time (NULL));

  if (argc == 2)
  {
    ASSERT (input = fopen (argv [1], "r"), ERR_FILE_OPEN);
    DOGEN = 0;
    inp ();
  }
  else gen ();

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;

  for (i = 0; i < asize; i ++)
  {
    for (j = 0; j < 3; j ++)
    {
      if (apoint [i][j] < extents [j]) extents [j] = apoint [i][j];
      if (apoint [i][j] > extents [3+j]) extents [3+j] = apoint [i][j];
    }
  }

  for (i = 0; i < bsize; i ++)
  {
    for (j = 0; j < 3; j ++)
    {
      if (bpoint [i][j] < extents [j]) extents [j] = bpoint [i][j];
      if (bpoint [i][j] > extents [3+j]) extents [3+j] = bpoint [i][j];
    }
  }

  for (i = 0; i < 3; i ++) d [i] = extents [3+i] - extents [i];

  for (i = 0; i < 3; i ++)
  {
    extents [3+i] += 0.25 * d [i];
    extents [i] -= 0.25 * d [i];
  }

  GLV (&argc, argv, "cvitest", 400, 400, extents, NULL,
    NULL, NULL, NULL, render, key, NULL, NULL, NULL, NULL);

  if (input) fclose (input);

  return 0;
}
#else
int main (int argc, char **atgv) { return 0; }
#endif
