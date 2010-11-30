/*
 * gjktest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of convex proximity
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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "alg.h"
#include "hul.h"
#include "gjk.h"

#define minim 4
#define limit 128
enum  {GEN, HUL, GJK} mode = HUL;
double apoint [limit][3],
       bpoint [limit][3];
int asize, bsize;
TRI *a = NULL,
    *b = NULL;
double center [3], radius; /* sphere if 'b' is NULL */
int alength,
    blength;
double p [3],
       q [3]; /* two closest points */

static void gen ()
{
  double move [3];

  if (rand () % 2 == 0) /* convex and sphere */
  {
    asize = bsize = 0;

    while (asize < minim) asize = rand () % limit;

    SETRAND (move, 1.0);

    for (int n = 0; n < asize; n ++)
    { SETRAND (apoint [n], 0.75);
      ADD (apoint [n], move, apoint [n]); }


    SETRAND (center, 1.0);
    SUB (center, move, center);
    radius = 0.75 * DRAND ();
  }
  else /* convex and convex */
  {
    asize = bsize = 0;

    while (asize < minim) asize = rand () % limit;

    while (bsize < minim) bsize = rand () % limit;

    SETRAND (move, 1.0);

    for (int n = 0; n < asize; n ++)
    { SETRAND (apoint [n], 0.75);
      ADD (apoint [n], move, apoint [n]); }

    for (int n = 0; n < bsize; n ++)
    { SETRAND (bpoint [n], 0.75);
      SUB (bpoint [n], move, bpoint [n]); }
  }
}

static void render (void)
{
  glPointSize (3.0);

  if (mode == GJK || mode == GEN)
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
    if (b)
    {
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
    else
    {
      glPushMatrix ();
      glTranslated (center [0], center [1], center [2]);
      glutSolidSphere (radius, 32, 32);
      glPopMatrix ();
    }

    if (mode == GEN)
    {
      glPointSize (3.0);
      glColor3d (0, 0, 1);
      glDisable (GL_LIGHTING);
      glBegin (GL_POINTS);
      glVertex3dv (p);
      glVertex3dv (q);
      glEnd ();
      glBegin (GL_LINES);
      glVertex3dv (p);
      glVertex3dv (q);
      glEnd ();
      glEnable (GL_LIGHTING);
      glPointSize (1.0);
    }
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
    if (bsize)
    {
      glBegin (GL_POINTS);
      for (int n = 0; n < bsize; n ++)
	glVertex3dv (bpoint [n]);
      glEnd ();
      glEnable (GL_LIGHTING);
    }
    else
    {
      glPushMatrix ();
      glTranslated (center [0], center [1], center [2]);
      glutSolidSphere (radius, 12, 12);
      glPopMatrix ();
    }
  }
}

static void key (int key, int x, int y)
{
  if (key == ' ')
  {

    switch (mode)
    {
      case GEN:
	gen ();
	mode = HUL;
      break;
      case HUL:
        free (a);
        free (b);
        a = hull ((double*)apoint, asize, &alength);
	if (bsize > 0)
          b = hull ((double*)bpoint, bsize, &blength);
	else { b = NULL; blength = 0; }
	mode = GJK;
      break;
      case GJK:
      {
	double *va, *vb, d [3];
	int nva, nvb;

	va = TRI_Vertices (a, alength, &nva);

	if (blength)
	{
	  vb = TRI_Vertices (b, blength, &nvb);
	  gjk (va, nva, vb, nvb, p, q);
	  free (vb);
	}
	else gjk_convex_sphere (va, nva, center, radius, p, q);

	SUB (p, q, d);
	printf ("|p-q|=%g\n", LEN (d));

	free (va);
	mode = GEN;
      }
      break;
    }

    GLV_Redraw_All ();
  }
}

int main (int argc, char **argv)
{
  double extents [6] =
  { -2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  
  printf ("SPACE - iterate over the test stages\n");

  srand ((unsigned) time (NULL));

  gen ();

  GLV (&argc, argv, "gjktest", 400, 400, extents, NULL,
    NULL, NULL, NULL, render, key, NULL, NULL, NULL, NULL);

  return 0;
}
#else
int main (int argc, char **argv) { return 0; }
#endif
