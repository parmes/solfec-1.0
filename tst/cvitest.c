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
#endif

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "alg.h"
#include "hul.h"
#include "cvi.h"

#define minim 4
#define limit 128
enum  {GEN, HUL, CVI} mode = HUL;
double apoint [limit][3],
       bpoint [limit][3];
int asize, bsize;
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
        b = hull ((double*)bpoint, bsize, &blength);
	mode = CVI;
      break;
      case CVI:
      {
	double *va, *pa,
	       *vb, *pb;
	int nva, npa,
	    nvb, npb;

	va = TRI_Vertices (a, alength, &nva);
	pa = TRI_Planes (a, alength, &npa);
	vb = TRI_Vertices (b, blength, &nvb);
	pb = TRI_Planes (b, blength, &npb);
        free (c);
	c = cvi (va, nva, pa, npa,
	         vb, nvb, pb, npb,
		 &clength);
	mode = GEN;
	free (va);
	free (pa);
	free (vb);
	free (pb);
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

  GLV (&argc, argv, "cvitest", 400, 400, extents, NULL,
    NULL, NULL, NULL, render, key, NULL, NULL, NULL, NULL);

  return 0;
}
