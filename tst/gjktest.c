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

enum {GJK_CONVEX_CONVEX = 1,
      GJK_CONVEX_SPHERE,
      GJK_CONVEX_POINT,
      GJK_CONVEX_ELLIP,
      GJK_SPHERE_ELLIP,
      GJK_ELLIP_ELLIP,
      GJK_ELLIP_POINT} mode = GJK_CONVEX_CONVEX; /* test mode */

/* convex data */
#define minim 4
#define limit 128
double apoint [limit][3], bpoint [limit][3];
int asize, bsize;
TRI *a = NULL, *b = NULL;
int alength, blength;

/* sphere data */
double center [3], radius; /* sphere if 'b' is NULL */

/* ellipsoid data */
double el1_center [3], el1_sca [3], el1_rot [9];
double el2_center [3], el2_sca [3], el2_rot [9];

/* two closest points */
double p [3], q [3];

/* generate test */
static void gen ()
{
  double d [3], move [3];

  switch (mode)
  {
  case GJK_CONVEX_CONVEX:
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

    free (a);
    free (b);
    a = hull ((double*)apoint, asize, &alength);
    b = hull ((double*)bpoint, bsize, &blength);

    double *va, *vb;
    int nva, nvb;

    va = TRI_Vertices (a, alength, &nva);
    vb = TRI_Vertices (b, blength, &nvb);

    gjk (va, nva, vb, nvb, p, q);

    free (va);
    free (vb);
  }
  break;
  case GJK_CONVEX_SPHERE:
  {
    asize = bsize = 0;

    while (asize < minim) asize = rand () % limit;

    SETRAND (move, 1.0);

    for (int n = 0; n < asize; n ++)
    { SETRAND (apoint [n], 0.75);
      ADD (apoint [n], move, apoint [n]); }

    free (a);
    a = hull ((double*)apoint, asize, &alength);

    SETRAND (center, 1.0);
    SUB (center, move, center);
    radius = 0.75 * DRAND ();

    double *va;
    int nva;

    va = TRI_Vertices (a, alength, &nva);

    gjk_convex_sphere (va, nva, center, radius, p, q);

    free (va);
  }
  break;
  case GJK_CONVEX_POINT:
  {
    asize = bsize = 0;

    while (asize < minim) asize = rand () % limit;

    SETRAND (move, 1.0);

    for (int n = 0; n < asize; n ++)
    { SETRAND (apoint [n], 0.75);
      ADD (apoint [n], move, apoint [n]); }

    free (a);
    a = hull ((double*)apoint, asize, &alength);

    SETRAND (center, 1.0);
    SUB (center, move, center);
    radius = 0.0;

    double *va;
    int nva;

    va = TRI_Vertices (a, alength, &nva);

    COPY (center, p);
    gjk_convex_point (va, nva, p, q);

    free (va);
  }
  break;
  case GJK_CONVEX_ELLIP:
  {
    asize = bsize = 0;

    while (asize < minim) asize = rand () % limit;

    SETRAND (move, 1.0);

    for (int n = 0; n < asize; n ++)
    { SETRAND (apoint [n], 0.75);
      ADD (apoint [n], move, apoint [n]); }

    free (a);
    a = hull ((double*)apoint, asize, &alength);

    SETRAND (el1_center, 1.0);
    SUB (el1_center, move, el1_center);
    el1_sca [0] = 0.75 * DRAND ();
    el1_sca [1] = 0.75 * DRAND ();
    el1_sca [2] = 0.75 * DRAND ();
    EXPMAP (el1_sca, el1_rot);

    double *va;
    int nva;

    va = TRI_Vertices (a, alength, &nva);

    gjk_convex_ellip (va, nva, el1_center, el1_sca, el1_rot, p, q);

    free (va);
  }
  break;
  case GJK_SPHERE_ELLIP:
  {
    SETRAND (move, 1.0);

    SETRAND (center, 1.0);
    radius = 0.75 * DRAND ();

    SETRAND (el1_center, 1.0);
    SUB (el1_center, move, el1_center);
    el1_sca [0] = 0.75 * DRAND ();
    el1_sca [1] = 0.75 * DRAND ();
    el1_sca [2] = 0.75 * DRAND ();
    EXPMAP (el1_sca, el1_rot);

    gjk_sphere_ellip (center, radius, el1_center, el1_sca, el1_rot, p, q);
  }
  break;
  case GJK_ELLIP_ELLIP:
  {
    SETRAND (move, 1.0);

    SETRAND (el1_center, 1.0);
    el1_sca [0] = 0.75 * DRAND ();
    el1_sca [1] = 0.75 * DRAND ();
    el1_sca [2] = 0.75 * DRAND ();
    EXPMAP (el1_sca, el1_rot);

    SETRAND (el2_center, 1.0);
    SUB (el2_center, move, el2_center);
    el2_sca [0] = 0.75 * DRAND ();
    el2_sca [1] = 0.75 * DRAND ();
    el2_sca [2] = 0.75 * DRAND ();
    EXPMAP (el2_sca, el2_rot);

    gjk_ellip_ellip (el1_center, el1_sca, el1_rot, el2_center, el2_sca, el2_rot, p, q);
  }
  break;
  case GJK_ELLIP_POINT:
  {
    SETRAND (move, 1.0);

    SETRAND (el1_center, 1.0);
    el1_sca [0] = 0.75 * DRAND ();
    el1_sca [1] = 0.75 * DRAND ();
    el1_sca [2] = 0.75 * DRAND ();
    EXPMAP (el1_sca, el1_rot);

    SETRAND (center, 1.0);
    SUB (center, move, center);
    radius = 0.0;

    COPY (center, p);
    gjk_ellip_point (el1_center, el1_sca, el1_rot, p, q);
  }
  break;
  }

  SUB (p, q, d);
  printf ("|p-q|=%g\n", LEN (d));
}

static void render (void)
{
  glPointSize (3.0);

  glEnable (GL_NORMALIZE);

  if (mode >= GJK_CONVEX_CONVEX &&
      mode <= GJK_CONVEX_ELLIP)
  {
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
  }

  if (mode == GJK_CONVEX_CONVEX)
  {
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

  if (mode == GJK_CONVEX_SPHERE || mode == GJK_SPHERE_ELLIP)
  {
    if (mode == GJK_CONVEX_SPHERE) glColor3d (0.0, 1.0, 0.0);
    else glColor3d (1.0, 0.0, 0.0);
    glPushMatrix ();
    glTranslated (center [0], center [1], center [2]);
    glutSolidSphere (radius, 32, 32);
    glPopMatrix ();
  }

  if (mode == GJK_CONVEX_POINT || mode == GJK_ELLIP_POINT)
  { 
    glColor3d (0.0, 1.0, 0.0);
    glPointSize (3.0);
    glBegin (GL_POINTS);
    glVertex3dv (center);
    glEnd ();
    glPointSize (1.0);
  }

  if (mode >= GJK_CONVEX_ELLIP && mode <= GJK_ELLIP_POINT)
  {
    if (mode <= GJK_SPHERE_ELLIP) glColor3d (0.0, 1.0, 0.0);
    else glColor3d (1.0, 0.0, 0.0);
    glPushMatrix ();
    glTranslated (el1_center [0], el1_center [1], el1_center [2]);
    GLfloat m [16] = {el1_rot [0], el1_rot [1], el1_rot [2], 0,
                      el1_rot [3], el1_rot [4], el1_rot [5], 0,
		      el1_rot [6], el1_rot [7], el1_rot [8], 0,
		      0          , 0          , 0          , 1};
    glMultMatrixf (m);
    glScaled (el1_sca [0], el1_sca [1], el1_sca [2]);
    glutSolidSphere (1.0, 32, 32);
    glPopMatrix ();
  }

  if (mode == GJK_ELLIP_ELLIP)
  {
    glColor3d (0.0, 1.0, 0.0);
    glPushMatrix ();
    glTranslated (el2_center [0], el2_center [1], el2_center [2]);
    GLfloat m [16] = {el2_rot [0], el2_rot [1], el2_rot [2], 0,
                      el2_rot [3], el2_rot [4], el2_rot [5], 0,
		      el2_rot [6], el2_rot [7], el2_rot [8], 0,
		      0          , 0          , 0          , 1};
    glMultMatrixf (m);
    glScaled (el2_sca [0], el2_sca [1], el2_sca [2]);
    glutSolidSphere (1.0, 32, 32);
    glPopMatrix ();
  }

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

static void key (int key, int x, int y)
{
  switch (key)
  {
  case '1': mode = GJK_CONVEX_CONVEX; break;
  case '2': mode = GJK_CONVEX_SPHERE; break;
  case '3': mode = GJK_CONVEX_POINT; break;
  case '4': mode = GJK_CONVEX_ELLIP; break;
  case '5': mode = GJK_SPHERE_ELLIP; break;
  case '6': mode = GJK_ELLIP_ELLIP; break;
  case '7': mode = GJK_ELLIP_POINT; break;
  }

  gen ();

  GLV_Redraw_All ();
}

int main (int argc, char **argv)
{
  double extents [6] =
  { -2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  
  printf ("1 - convex to convex\n");
  printf ("2 - convex to sphere\n");
  printf ("3 - convex to point\n");
  printf ("4 - convex to ellipsoid\n");
  printf ("5 - sphere to ellipsoid\n");
  printf ("6 - ellipsoid to ellipsoid\n");
  printf ("7 - ellipsoid to point\n");
  printf ("OTHER - regeneratge current test case\n");

  srand ((unsigned) time (NULL));

  gen ();

  GLV (&argc, argv, "gjktest", 400, 400, extents, NULL,
    NULL, NULL, NULL, render, key, NULL, NULL, NULL, NULL);

  return 0;
}
#else
int main (int argc, char **argv) { return 0; }
#endif
