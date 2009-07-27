/*
 * hultest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of convex hull solver
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
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "alg.h"
#include "hul.h"

#define minim 4
#define limit 128
enum  {GEN, HUL} mode = HUL;
double point [limit][3];
int size;
TRI *tri = NULL;
int length;


static void gen ()
{
  while (size < minim) size = rand () % limit;

  for (int n = 0; n < size; n ++)
  { SETRAND (point [n], 1.0); }
}

static void render (void)
{
  if (mode == GEN)
  {
    glColor3d (0.5, 0.5, 0.5);
    glEnable (GL_NORMALIZE);
    glBegin (GL_TRIANGLES);
    for (TRI *t = tri, *e = t + length; t < e; t ++)
    {
      glNormal3dv (t->out);
      glVertex3dv (t->ver [0]);
      glVertex3dv (t->ver [1]);
      glVertex3dv (t->ver [2]);
    }
    glEnd ();

    glDisable (GL_LIGHTING);
    glColor3d (1.0, 0.0, 0.0);
    glBegin (GL_LINES);
    for (TRI *t = tri, *e = t + length; t < e; t ++)
    {
      double p [3], q [3];
      ADD (t->ver [0], t->ver [1], p);
      ADD (t->ver [2], p, p);
      DIV (p, 3.0, p);
      ADDMUL (p, 0.5, t->out, q);
      glVertex3dv (p);
      glVertex3dv (q);
    }
    glEnd ();
    glDisable (GL_LIGHTING);
  }

  glDisable (GL_LIGHTING);
  glPointSize (3.0);
  glColor3d (0.3, 1.0, 0.3);
  glBegin (GL_POINTS);
  for (int n = 0; n < size; n ++)
    glVertex3dv (point [n]);
  glEnd ();
  glPointSize (1.0);
  glEnable (GL_LIGHTING);
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
        free (tri);
        tri = hull ((double*)point, size, &length);
	mode = GEN;
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

  GLV (&argc, argv, "hultest", 400, 400, extents, NULL,
    NULL, NULL, NULL, render, key, NULL, NULL, NULL, NULL);

  return 0;
}
