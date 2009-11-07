/*
 * rnd.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * OpenGL rendering
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

#if __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "set.h"
#include "alg.h"
#include "glv.h"
#include "sol.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "shp.h"
#include "lng.h"
#include "rnd.h"
#include "err.h"

typedef struct body_data BODY_DATA; /* body rendering data */

struct body_data
{
  enum {HIDDEN     = 0x01,         /* hidden body */
        OFF        = 0x02,         /* disabled body */
	ROUGH_MESH = 0x04} flags;  /* rough mesh rendering */
};

int RND_Menu (char ***names, int **codes)
{
  return 0;
}

void RND_Init ()
{
}

int  RND_Idle ()
{
  return 0;
}

void RND_Quit ()
{
}

void RND_Render ()
{
}

void RND_Key (int key, int x, int y)
{
}

void RND_Keyspec (int key, int x, int y)
{
}

void RND_Mouse (int button, int state, int x, int y)
{
}

void RND_Motion (int x, int y)
{
}

void RND_Passive (int x, int y)
{
}

void RND_Viewer_On ()
{
}

void RND_Domain (DOM *dom)
{
}

void RND_Solver (DOM *dom, int kind, void *solver)
{
}

int  RND_On ()
{
  return 1;
}

void RND_Free_Rendering_Data (void *data)
{
}
