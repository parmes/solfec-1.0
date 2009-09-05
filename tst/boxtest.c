/*
 * boxtest.c
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * test of BOX overlap detection
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

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#include "mem.h"
#include "err.h"
#include "alg.h"
#include "box.h"
#include "bod.h"

#if OPENGL
#if __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif
#include "glv.h"
#endif

#define INSIDE(point, extents)\
  ((point) [0] >= (extents) [0] &&\
   (point) [0] <= (extents) [3] &&\
   (point) [1] >= (extents) [1] &&\
   (point) [1] <= (extents) [4] &&\
   (point) [2] >= (extents) [2] &&\
   (point) [2] <= (extents) [5])

typedef struct testbox TESTBOX;

struct testbox
{
  double coord [6]; /* box position */

  double velo [3]; /* box velocity */

  BOX *box; /* related box object */

  SGP *sgp;

  TESTBOX *p, *n; /* list links */
};

/* global bounding box */
static double extents [6] = {-1000., -1000., -1000., 1000., 1000., 1000.};

/* extents volume */
static double volume = 8e9;

/* box extents extremas */
static double box_wx = 100.;
static double box_wy = 100.;
static double box_wz = 100.;

/* global dummy body */
static BODY bod;

/* shape object pairs */
SGP *sgp = NULL;

/* vector of boxes */
static TESTBOX *box = NULL;

/* number of particles */
static int boxsize = 4096;

/* list of deleted boxes */
static TESTBOX *deleted = NULL;

/* global AABB context */
static AABB *aabb = NULL;

/* insertion/deletion probability */
static double insdelprob = 0.01;

/* gravity acceleration */
static double gravity [3] = {0., 0., -.25};

/* gravity flag */
static int gravity_exists = 0;

/* time step */
static double timestep = 1.;

/* collision restitution coefficient */
static double restitution = .9;

/* current frame number */
static int frame = 0;

/* algorithm kind */
BOXALG algorithm = HYBRID;

/* number of overlaps */
int noverlaps = 0;

/* initial box arrangements */
enum {BRAND, BADJ};

#if OPENGL
/* view extents */
static double vextents [6] = {-2000., -2000., -2000., 2000., 2000., 2000.};

/* FPS rate */
static int fps = 0;

/* boxes on/off */
static int boxeson = 1;

/* overlaps graph on/off */
static int overlapson = 1;

/* time stepping on/off */
static int dostep = 1;

/* current arrangement */
static short curarr = BRAND;
#endif

/* BOX related callbacks */
static void box_extents_update (void *data, TESTBOX *box, double *extents)
{
  double epsilon [3];

  SET (epsilon, GEOMETRIC_EPSILON);
  COPY6 (box->coord, extents);
  SUB (extents, epsilon, extents);
  ADD (extents+3, epsilon, extents+3);
}

static void* box_overlap_create (void *data, BOX *one, BOX *two)
{
  noverlaps ++;
  return one; /* return a valid pointer */
}

static void box_overlap_release (void *data, BOX *one, BOX *two) 
{
  noverlaps --;
}

/* assign a random coordinate within the
 * globally specified extents */
static void random_coord (double *coord)
{
  coord [0] = extents [0] + (extents [3] - extents [0] - box_wx) * DRAND ();  
  coord [1] = extents [1] + (extents [4] - extents [1] - box_wy) * DRAND ();  
  coord [2] = extents [2] + (extents [5] - extents [2] - box_wz) * DRAND ();  
  coord [3] = coord [0] + box_wx * DRAND ();
  coord [4] = coord [1] + box_wy * DRAND ();
  coord [5] = coord [2] + box_wz * DRAND ();
}

/* assign a random velocity */
static void random_velo (double *velo)
{
  velo [0] = -.5 + DRAND (); 
  velo [1] = -.5 + DRAND ();   
  velo [2] = -.5 + DRAND ();  
}

/* initialise boxes and BOX set */
static void generate_box_set (int howmany, short arrange)
{
  int i;

  /* enable self-contact (we use one dummy body) */
  bod.flags |= BODY_DETECT_SELF_CONTACT;

  /* adjust edge widths to the number of boxes */
  box_wx = box_wy = box_wz = pow (volume/(double)howmany, 0.33),

  /* initialise random generator */
  srand ((unsigned)time (NULL));

  /* initialise particles memory */
  howmany = MAX (howmany, 1);
  ERRMEM (box = realloc (box, sizeof (TESTBOX) * howmany));
  ERRMEM (sgp = realloc (sgp, sizeof (SGP) * howmany));
  boxsize = howmany;
  deleted = NULL; /* empty deleted boxes list */

  /* create initial BOX set */
  if (aabb) AABB_Destroy (aabb);
  aabb = AABB_Create (howmany);

  /* let the dummy body be rigid */
  bod.kind = RIG;
 
  switch (arrange)
  {
  case BRAND:
  {
    /* initialise particles */
    for (i = 0; i < howmany; i ++)
    {
      random_coord (box [i].coord);
      random_velo (box [i].velo);

      sgp [i].gobj = &box [i];
      box [i].sgp = &sgp [i];
      box [i].box = AABB_Insert (aabb, &bod, GOBJ_DUMMY, &sgp [i], NULL, (BOX_Extents_Update)box_extents_update);
    }
  }
  break;
  case BADJ:
  {
    double step = 0.99 * pow (volume/(double)howmany, 1.0/3.0),
	   x, y, z;

    i = 0;
    for (x = extents [0]; x < extents [3] - step; x += step)
    for (y = extents [1]; y < extents [4] - step; y += step)
    for (z = extents [2]; z < extents [5] - step; z += step)
    {
      if (i < howmany)
      {
	box [i].coord [0] = x;
	box [i].coord [1] = y;
	box [i].coord [2] = z;
	box [i].coord [3] = x + step;
	box [i].coord [4] = y + step;
	box [i].coord [5] = z + step;
        random_velo (box [i].velo);

        sgp [i].gobj = &box [i];
        box [i].sgp = &sgp [i];
	box [i].box = AABB_Insert (aabb, &bod, GOBJ_DUMMY, &sgp [i], NULL, (BOX_Extents_Update)box_extents_update);
	i ++;
      }
    }
    boxsize = i;
  }
  break;
  }
}

/* change velocity after the wall impact */
static void velocity_change (double *velo, double nx, double ny, double nz)
{
  double proj [3], len, dif [3];

  len = velo [0] * nx + velo [1] * ny + velo [2] * nz;
  proj [0] = - nx * len;
  proj [1] = - ny * len;
  proj [2] = - nz * len;
  ADD (velo, proj, dif);
  ADD (proj, dif, velo);
  SCALE (velo, restitution);
}

/* advance box motion by one time step */
static void box_motion_step ()
{
  double pmid [3], rmid [3], dir [3], rel [3];
  MAP *item;
  BOX *obi, *obj;
  TESTBOX *p, *r;
  double *q, *u;

  /* update velocities in order to avoid collisions */
  for (obj = aabb->lst; obj; obj = obj->next)
  {
    for (item = MAP_First (obj->adj); item; item = MAP_Next (item))
    {
      obi = item->key;
      if (obi < obj)
      {
	p = obi->sgp->gobj;
	r = obj->sgp->gobj;
	MID (p->coord, p->coord + 3, pmid);
	MID (r->coord, r->coord + 3, rmid);
	SUB (rmid, pmid, dir); NORMALIZE (dir);
	SUB (r->velo, p->velo, rel);
	if (DOT (rel, dir) < 0.)
	{
	  velocity_change (p->velo, dir [0], dir [1], dir [2]);
	  velocity_change (r->velo, -dir [0], -dir [1], -dir [2]);
	}
      }
    }
  }

  /* update positions */
  for (obj = aabb->lst; obj; obj = obj->next)
  {
    p = obj->sgp->gobj;
    q = p->coord;
    u = p->velo;
    
    if (gravity_exists)
    {
      ADDMUL (u, timestep, gravity, u);
    }

    ADDMUL (q, timestep, u, q);
    ADDMUL (q+3, timestep, u, q+3);

    if (q [0] < extents [0])
    {
      ADDMUL (q, -timestep, u, q);
      ADDMUL (q+3, -timestep, u, q+3);
      velocity_change (u, -1, 0, 0);
    }
    else if (q [1] < extents [1])
    {
      ADDMUL (q, -timestep, u, q);
      ADDMUL (q+3, -timestep, u, q+3);
      velocity_change (u, 0, -1, 0);
    }
    else if (q [2] < extents [2])
    {
      ADDMUL (q, -timestep, u, q);
      ADDMUL (q+3, -timestep, u, q+3);
      velocity_change (u, 0, 0, -1);
    }
    else if (q [3] > extents [3])
    {
      ADDMUL (q, -timestep, u, q);
      ADDMUL (q+3, -timestep, u, q+3);
      velocity_change (u, 1, 0, 0);
    }
    else if (q [4] > extents [4])
    {
      ADDMUL (q, -timestep, u, q);
      ADDMUL (q+3, -timestep, u, q+3);
      velocity_change (u, 0, 1, 0);
    }
    else if (q [5] > extents [5])
    {
      ADDMUL (q, -timestep, u, q);
      ADDMUL (q+3, -timestep, u, q+3);
      velocity_change (u, 0, 0, 1);
    }
  }
}

/* single computational step */
static void single_computational_step ()
{
  BOX *obj;
  TESTBOX *p, *n;
  double prob;
  int i;
  
  /* advance motion */
  box_motion_step ();

  /* do some insertion deletions */
  for (i = 0, p = deleted; p; p = p->n, i ++);
  for (p = deleted, prob = insdelprob / (double)i; p; p = n)
  {
    n = p->n;

    if (DRAND () < prob)
    {
      p = deleted;
      random_coord (p->coord);
      random_velo (p->velo);
      deleted = p->n; 

      p->box = AABB_Insert (aabb, &bod, GOBJ_DUMMY, p->sgp, NULL, (BOX_Extents_Update)box_extents_update);
      boxsize ++;
    } 
  }
 
  prob = insdelprob / (double)boxsize;
  for (obj = aabb->lst; obj; obj = obj->next)
  {
    if (DRAND () < prob)
    {
      p = obj->sgp->gobj;
      AABB_Delete (aabb, p->box);
      p->n = deleted;
      deleted = p;
      boxsize --;
    }	
  }

  /* update box overlaps */
  AABB_Update (aabb, algorithm, NULL,
    (BOX_Overlap_Create)box_overlap_create,
    (BOX_Overlap_Release)box_overlap_release);

  /* iterate frame */
  frame ++;
}

/* free everything */
static void free_all_data ()
{
  free (box);
  free (sgp);
  AABB_Destroy (aabb);
}

#if OPENGL
static void box_color (int boxnum, GLfloat *color)
{
  color [0] = (double)boxnum / (double)boxsize;
  color [1] = sin (color [0] * ALG_PI);
  color [2] = 1. - color [0];
}

/* FPS timer function */
static void view_fps_callback (int value)
{
  fps = frame;
  frame = 0;

  /* set up next trigger */
  glutTimerFunc (1000, view_fps_callback, 0);
}

/* render 2d statistics */
static void view_render2d (void)
{
  glClearColor (1, 1, 1, 1);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glColor3d (0., 0., 0.);
  GLV_Print (8, 8, 0, GLV_FONT_12,
    "FPS: %d, NCN: %d, NBX: %d, ALG: %s",
    fps, noverlaps, boxsize, AABB_Algorithm_Name (algorithm));
  glutSwapBuffers();
}

/* view initialisation */
static void view_init (void)
{
  /* open 2d window */
  GLV_Open_Window (0, 0, 1600, 24, view_render2d);

  /* set up first FPS trigger */
  glutTimerFunc (1000, view_fps_callback, 0);
}

#define GLBOX(extents)\
  glNormal3d (0, 0, -1);\
  glVertex3d (extents [0], extents [1], extents [2]);\
  glVertex3d (extents [0], extents [4], extents [2]);\
  glVertex3d (extents [3], extents [4], extents [2]);\
  glVertex3d (extents [3], extents [1], extents [2]);\
  glNormal3d (0, 0, 1);\
  glVertex3d (extents [0], extents [1],extents [5]);\
  glVertex3d (extents [3], extents [1],extents [5]);\
  glVertex3d (extents [3], extents [4],extents [5]);\
  glVertex3d (extents [0], extents [4],extents [5]);\
  glNormal3d (-1, 0, 0);\
  glVertex3d (extents [0], extents [1], extents [2]);\
  glVertex3d (extents [0], extents [1], extents [5]);\
  glVertex3d (extents [0], extents [4], extents [5]);\
  glVertex3d (extents [0], extents [4], extents [2]);\
  glNormal3d (1, 0, 0);\
  glVertex3d (extents [3], extents [1], extents [2]);\
  glVertex3d (extents [3], extents [4], extents [2]);\
  glVertex3d (extents [3], extents [4], extents [5]);\
  glVertex3d (extents [3], extents [1], extents [5]);\
  glNormal3d (0, -1, 0);\
  glVertex3d (extents [0], extents [1], extents [2]);\
  glVertex3d (extents [3], extents [1], extents [2]);\
  glVertex3d (extents [3], extents [1], extents [5]);\
  glVertex3d (extents [0], extents [1], extents [5]);\
  glNormal3d (0, 1, 0);\
  glVertex3d (extents [0], extents [4], extents [2]);\
  glVertex3d (extents [0], extents [4], extents [5]);\
  glVertex3d (extents [3], extents [4], extents [5]);\
  glVertex3d (extents [3], extents [4], extents [2])

/* render 3d scene */
static void view_render3d (void)
{
  double m1 [3], m2 [3], *ext;
  GLfloat color [3];
  TESTBOX *p;
  BOX *u, *v;
  MAP *item;
  int n;

  /* draw outline of the scene
   * bounding extents */
  glColor3d (.6, .6, .6);
  glBegin (GL_LINE_STRIP);
    glVertex3d (extents [0], extents [1], extents [2]);
    glVertex3d (extents [3], extents [1], extents [2]);
    glVertex3d (extents [3], extents [4], extents [2]);
    glVertex3d (extents [0], extents [4], extents [2]);
    glVertex3d (extents [0], extents [1], extents [2]);
    glVertex3d (extents [0], extents [1], extents [5]);
    glVertex3d (extents [3], extents [1], extents [5]);
    glVertex3d (extents [3], extents [4], extents [5]);
    glVertex3d (extents [0], extents [4], extents [5]);
    glVertex3d (extents [0], extents [1], extents [5]);
  glEnd ();
  glBegin (GL_LINES);
    glVertex3d (extents [3], extents [1], extents [2]);
    glVertex3d (extents [3], extents [1], extents [5]);
    glVertex3d (extents [3], extents [4], extents [2]);
    glVertex3d (extents [3], extents [4], extents [5]);
    glVertex3d (extents [0], extents [4], extents [2]);
    glVertex3d (extents [0], extents [4], extents [5]);
  glEnd ();

  /* draw boxes */
  if (boxeson)
  {
    glBegin (GL_QUADS);
    for (n = 0, u = aabb->lst; u; n ++, u = u->next)
    {
      box_color (n, color);
      glColor3fv (color);
      p = u->sgp->gobj;
      GLBOX (p->coord);
    }
    glEnd ();
  }

  if (overlapson)
  {
    /* render overlap graph */
    glBegin (GL_LINES);
    glColor3d (1., 0., 0.);
    for (u = aabb->lst; u; u = u->next)
    {
      for (item = MAP_First (u->adj); item; item = MAP_Next (item))
      {
	v = item->key;
	if (u < v)
	{
	  ext = u->extents;
	  MID (ext, ext + 3, m1);
	  ext = v->extents;
	  MID (ext, ext + 3, m2);
	  glVertex3dv (m1);
	  glVertex3dv (m2);
	}
      }
    }
    glEnd ();
  }
}

static void view_key (int key, int x, int y)
{
  switch (key)
  {
  case '1':
    {
      algorithm = HYBRID;
    }
    break;
  case '2':
    {
      algorithm = HASH3D;
    }
    break;
  case '3':
    {
      algorithm = SWEEP_HASH2D_LIST;
    }
    break;
  case '4':
    {
      algorithm = SWEEP_HASH1D_XYTREE;
    }
    break;
  case '5':
    {
      algorithm = SWEEP_HASH2D_XYTREE;
    }
    break;
  case '6':
    {
      algorithm = SWEEP_XYTREE;
    }
    break;
  case 'g':
    {
      gravity_exists = !gravity_exists;
    }
    break;
  case 'b':
    {
      boxeson = !boxeson;
    }
    break;
  case 'o':
    {
      overlapson = !overlapson;
    }
    break;
  case ' ':
    {
      dostep = !dostep;
    }
    break;
  };
}

static void view_key_spec (int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_F1:
      curarr = BRAND;
      generate_box_set (boxsize, BRAND);
    break;
    case GLUT_KEY_F2:
      curarr = BADJ;
      generate_box_set (boxsize, BADJ);
    break;
    case GLUT_KEY_UP:
      boxsize *= 2;
      generate_box_set (boxsize, curarr);
    break;
    case GLUT_KEY_DOWN:
      if (boxsize > 8) boxsize /= 2;
      generate_box_set (boxsize, curarr);
    break;
  }
}

static int view_idle (void)
{
  if (dostep) single_computational_step ();
  else
  {
    /* update box overlaps only */
    AABB_Update (aabb, algorithm, NULL,
      (BOX_Overlap_Create)box_overlap_create,
      (BOX_Overlap_Release)box_overlap_release);
  }

  return 1;
}

static void view_quit (void)
{
  free_all_data ();
}
#endif

int main (int argc, char **argv)
{
  if (argc == 2)
    generate_box_set (MAX (atoi (argv [1]), 8), BRAND);
  else generate_box_set (boxsize, BRAND);

#if OPENGL
  printf ("boxtest key bindings:\n");
  printf ("1 - HYBRID algorithm\n");
  printf ("2 - HASH3D algorithm\n");
  printf ("3 - SWEEP_HASH2D_LIST algorithm\n");
  printf ("4 - SWEEP_HASH1D_XYTREE algorithm\n");
  printf ("5 - SWEEP_HASH2D_XYTREE algorithm\n");
  printf ("6 - SWEEP_XYTREE algorithm\n");
  printf ("g - gravity on/off\n");
  printf ("b - boxes drawing on/off\n");
  printf ("o - overlaps graph drawing on/off\n");
  printf ("SPACE - time stepping on/off\n");
  printf ("UP or DOWN - number of boxes *2 or /2\n");
  printf ("F1 - random set generation\n");
  printf ("F2 - adjacent set generation\n");

  /* start viewer */
  GLV (&argc, argv, "BOXES test", 400, 400, vextents, NULL,
    view_init, view_idle, view_quit, view_render3d,
    view_key, view_key_spec, NULL, NULL, NULL);

#else
  for (argc = 0; argc < 512; argc ++)
    single_computational_step ();
  free_all_data ();
#endif

  return 0;
}
