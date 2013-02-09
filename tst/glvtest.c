/*
 * glvtest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of graphical viewer
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

#include <stdlib.h>
#include <time.h>
#include "glv.h"

static void render (void)
{
  glColor3d (0.5, 0.5, 0.5);
  glutSolidCube (1.0);
}

int main (int argc, char **argv)
{
  double extents [6] =
  { -2.0, -2.0, -2.0, 2.0, 2.0, 2.0};
  
  srand ((unsigned) time (NULL));

  GLV (&argc, argv, "glvtest", 400, 400, extents, NULL,
    NULL, NULL, NULL, render, NULL, NULL, NULL, NULL, NULL);

  return 0;
}
#else
int main (int argc, char **argv) { return 0; }
#endif
