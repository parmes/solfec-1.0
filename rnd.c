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

typedef union value_source VALUE_SOURCE; /* vertex value source */

union value_source
{
  ELEPNT *epn;

  double *pnt;
};

typedef struct body_data BODY_DATA; /* body rendering data */

struct body_data
{
  enum {TRANSPARENT = 0x01,         /* transparency flag */
        HIDDEN      = 0x02,         /* hidden state */
	ROUGH_MESH  = 0x04} flags;  /* rough mesh rendering */

  GLuint triangles, /* VBO of triangle vertices, normals and colors */
	 lines; /* VBO of line vertices */

  GLsizei triangles_count,
	  lines_count;

  double **vertex_sources,
	 **normal_sources,
	 **vertex_value_sources,
	 **line_sources;

  double *values; /* vertex scalar field values */

  VALUE_SOURCE *value_sources; /* value sources */

  int values_count;

  ELEPNT **color_sources; /* only for FEM */

  double **spheres; /* center and radius and three points pointer tuples */

  GLfloat *sphere_colors; /* sphere surface colors */

  int spheres_count;
};

typedef struct solver_data SOLVER_DATA; /* solver interface */

struct solver_data
{
  int kind;

  void *solver;
};

typedef struct pointer_pair POINTER_PAIR; /* pointer pair */

struct pointer_pair
{
  double *one, *two;
};

typedef struct selection SELECTION; /* body selection set */

struct selection
{
  SET *set;

  SELECTION *prev;
};

enum {MENU_DOMAIN = 0, MENU_RENDER, MENU_TOOLS, MENU_ANALYSIS, MENU_KINDS, MENU_RESULTS, MENU_LAST}; /* menu identifiers */

enum /* menu items */
{
  DOMAIN_NEXT,
  DOMAIN_PREVIOUS,
  RENDER_SELECTION_2D,
  RENDER_SELECTION_3D,
  RENDER_PREVIOUS_SELECTION,
  TOOLS_TRANSPARENT,
  TOOLS_HIDE,
  TOOLS_SHOW_ALL,
  TOOLS_CLIPPING_PLANE,
  TOOLS_ROUGH_MESH,
  ANALYSIS_RUN,
  ANALYSIS_STOP,
  ANALYSIS_STEP,
  ANALYSIS_SEEKTO,
  ANALYSIS_FORWARD,
  ANALYSIS_BACKWARD,
  ANALYSIS_SKIP,
  KINDS_OF_CONSTRAINTS,
  KINDS_OF_FORCES,
  KINDS_OF_BODIES,
  KINDS_OF_SURFACES,
  KINDS_OF_VOLUMES,
  RESULTS_DX,
  RESULTS_DY,
  RESULTS_DZ,
  RESULTS_VX,
  RESULTS_VY,
  RESULTS_VZ,
  RESULTS_SX,
  RESULTS_SY,
  RESULTS_SZ,
  RESULTS_SXY,
  RESULTS_SXZ,
  RESULTS_SYZ,
  RESULTS_MISES,
  RESULTS_RT,
  RESULTS_RN,
  RESULTS_R
};

/* declarations */

static void menu_analysis (int);

/* global data */

#define CHUNK 256 /* memory pool items chunk */

#define BIGCHUNK 1024 /* memory pool items big chunk */

static short enabled = 0; /* renderin on/off */

static DOM *domain = NULL; /* current domain */

#define solfec ((SOLFEC*)domain->owner) /* current solfec */

static MAP *solvers = NULL; /* maps domains to solvers */

static char* menu_name [MENU_LAST];  /* menu names */

static int menu_code [MENU_LAST]; /* menu codes */

static SELECTION *selection = NULL; /* current selection */

static MEM selsetmem; /* selection sets memory */

static int skip_steps = 1; /* number of steps to skip when rewinding analysis */

/* push new selection on stack */
static void selection_push (SET *set)
{
  SELECTION *s;

  ERRMEM (s = malloc (sizeof (SELECTION)));
  s->prev = selection;
  s->set = set;

  selection = s;
}

/* translate scalar value into color */
inline static void value_to_color (double value, GLfloat *color)
{
  SET (color, 0.8);
}

/* obtain scalar point value */
static double point_value (BODY *bod, double *point)
{
  return 0.0;
}

/* compare pointer pair */
static int paircompare (POINTER_PAIR *a, POINTER_PAIR *b)
{
  if (a->one < b->one) return -1;
  else if (a->one == b->one)
  {
    if (a->two < b->two) return -1;
    else if (a->two == b->two) return 0;
  }

  return 1;
}

/* insert a pointer pair into the set */
static void register_line (MEM *pairmem, MEM *setmem, SET **lset, double *a, double *b)
{
  POINTER_PAIR *pair;

  ERRMEM (pair = MEM_Alloc (pairmem));

  pair->one = (a < b ? a : b);
  pair->two = (a > b ? a : b);

  if (!SET_Insert (setmem, lset, pair, (SET_Compare) paircompare)) MEM_Free (pairmem, pair);
}

/* create body rendering data */
static BODY_DATA* create_body_data (BODY *bod)
{
  double **vsr, **nsr, **lsr, **vvs, **end, *pla;
  GLfloat *ver, *v, *nor, *n, *col, *c, *lin, *l;
  MEM pairmem, mapmem, setmem;
  VALUE_SOURCE *source;
  POINTER_PAIR *pair;
  MAP *vmap, *jtem;
  SET *lset, *item;
  BODY_DATA *data;
  int i, j, *f;
  ELEMENT *ele;
  SPHERE *sph;
  CONVEX *cvx;
  SHAPE *shp;
  MESH *msh;
  FACE *fac;

  ERRMEM (data = calloc (sizeof (BODY_DATA), 1));

  MEM_Init (&pairmem, sizeof (POINTER_PAIR), CHUNK);
  MEM_Init (&mapmem, sizeof (MAP), CHUNK);
  MEM_Init (&setmem, sizeof (SET), CHUNK);
  vmap = NULL;
  lset = NULL;

  for (shp = bod->shape; shp; shp = shp->next)
  {
    switch (shp->kind)
    {
    case SHAPE_MESH:
      msh = shp->data;
      for (ele = msh->surfeles; ele; ele = ele->next)
      {
	for (fac = ele->faces; fac; fac = fac->next)
	{
	  data->triangles_count += (fac->type - 2);

	  for (i = 0; i < fac->type; i ++) MAP_Insert (&mapmem, &vmap, &msh->cur_nodes [fac->nodes [i]][0], NULL, NULL);

	  for (i = 0; i < fac->type - 1; i ++)
	    register_line (&pairmem, &setmem, &lset, &msh->cur_nodes [fac->nodes [i]][0], &msh->cur_nodes[fac->nodes [i+1]][0]);
	  register_line (&pairmem, &setmem, &lset, &msh->cur_nodes [fac->nodes [i]][0], &msh->cur_nodes [fac->nodes [0]][0]);
	}
      }
      break;
    case SHAPE_CONVEX:
      for (cvx = shp->data; cvx; cvx = cvx->next)
      {
	for (f = cvx->fac, j = 0; j < cvx->nfac; f += f[0]+1, j ++)
	{
	  data->triangles_count += (f[0] - 2);

	  for (i = 1; i <= f[0]; i ++) MAP_Insert (&mapmem, &vmap, &cvx->cur [f[i]], cvx, NULL);

	  for (i = 1; i <= f[0]-1; i ++)
	    register_line (&pairmem, &setmem, &lset, &cvx->cur [f[i]], &cvx->cur [f[i+1]]);
	  register_line (&pairmem, &setmem, &lset, &cvx->cur [f[i]], &cvx->cur [f[1]]);
	}
      }
      break;
    case SHAPE_SPHERE:
      for (sph = shp->data; sph; sph = sph->next)
      {
	data->spheres_count ++;

	ERRMEM (data->spheres = realloc (data->spheres, data->spheres_count * sizeof (double*) * 5));
	j = (data->spheres_count - 1) * 5;
	data->spheres [j] = sph->cur_center;
	data->spheres [j+1] = &sph->cur_radius;
	data->spheres [j+2] = &sph->cur_points [0][0];
	data->spheres [j+3] = &sph->cur_points [1][0];
	data->spheres [j+4] = &sph->cur_points [2][0];
      }
      break;
    }
  }

  data->lines_count = SET_Size (lset);
  ERRMEM (lin = malloc (data->lines_count * sizeof (GLfloat) * 6));
  ERRMEM (ver = malloc (data->triangles_count * sizeof (GLfloat) * 27));
  nor = ver + data->triangles_count * 9;
  col = nor + data->triangles_count * 9;
  ERRMEM (data->vertex_sources = malloc (data->triangles_count * sizeof (double*) * 9 + data->lines_count * sizeof (double*) * 2));
  data->normal_sources = data->vertex_sources + data->triangles_count * 3;
  data->vertex_value_sources = data->normal_sources + data->triangles_count * 3;
  data->line_sources = data->vertex_value_sources + data->triangles_count * 3;

  for (item = SET_First (lset), lsr = data->line_sources; item; item = SET_Next (item), lsr += 2)
  {
    pair = item->data;
    lsr [0] = pair->one;
    lsr [1] = pair->two;
  }

  vsr = data->vertex_sources;
  nsr = data->normal_sources;

  for (shp = bod->shape; shp; shp = shp->next)
  {
    switch (shp->kind)
    {
    case SHAPE_MESH:
      msh = shp->data;
      for (ele = msh->surfeles; ele; ele = ele->next)
      {
	for (fac = ele->faces; fac; fac = fac->next)
	{
	  for (i = 1; i < fac->type - 1; i ++, vsr += 3, nsr += 3)
	  {
	     vsr [0] = &msh->cur_nodes [fac->nodes [0]][0];
	     vsr [1] = &msh->cur_nodes [fac->nodes [i]][0];
	     vsr [2] = &msh->cur_nodes [fac->nodes [i+1]][0];
	     nsr [0] = nsr [1] = nsr [2] = fac->normal;
	  }
	}
      }
      break;
    case SHAPE_CONVEX:
      for (cvx = shp->data; cvx; cvx = cvx->next)
      {
	for (f = cvx->fac, j = 0, pla = cvx->pla; j < cvx->nfac; f += f[0]+1, j ++, pla += 4)
	{
	  for (i = 2; i <= f[0]-1; i ++, vsr += 3, nsr += 3)
	  {
	    vsr [0] = &cvx->cur [f[1]];
	    vsr [1] = &cvx->cur [f[i]];
	    vsr [2] = &cvx->cur [f[i+1]];
	    nsr [0] = nsr [1] = nsr [2] = pla;
	  }
	}
      }
      break;
    case SHAPE_SPHERE: break;
    }
  }

  data->values_count = MAP_Size (vmap); /* number of unique vertices */
  ERRMEM (data->values = calloc (data->values_count, sizeof (double)));
  ERRMEM (data->value_sources = malloc (data->values_count * sizeof (VALUE_SOURCE)));

  for (source = data->value_sources, jtem = MAP_First (vmap); jtem; source ++, jtem = MAP_Next (jtem))
  {
    if (bod->msh)
    {
      cvx = jtem->data; /* must have been mapped to a convex */
      source->epn = &cvx->epn [((double*)jtem->key - cvx->cur) / 3]; /* extract ELEPNT */
    }
    else source->pnt = jtem->key;

    jtem->data = &data->values [source - data->value_sources]; /* map to source */
  }

  for (vsr = data->vertex_sources,
       nsr = data->normal_sources,
       vvs = data->vertex_value_sources,
       end = vsr + data->triangles_count * 3,
       v = ver, n = nor, c = col; vsr < end;
       vsr ++, nsr ++, vvs ++, v += 3, n += 3, c += 3)
  {
    COPY (*vsr, v);
    COPY (*nsr, n);
    ASSERT_DEBUG_EXT (*vvs = MAP_Find (vmap, *vsr, NULL), "Inconsistent vertex mapping");
    value_to_color (**vvs, c);
  }

  for (lsr = data->line_sources,
       end = lsr + data->lines_count * 2,
       l = lin; lsr < end; lsr ++, l += 3)
  {
    COPY (*lsr, l);
  }

  if (data->spheres_count)
  {
    ERRMEM (data->sphere_colors = malloc (data->spheres_count * sizeof (GLfloat) * 3));

    for (c = data->sphere_colors, vsr = data->spheres, v = c + data->spheres_count * 3; c < v; c += 3, vsr += 5)
      value_to_color (point_value (bod, *vsr), c);
  }

  glGenBuffersARB (1, &data->lines);
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->lines);
  glBufferDataARB (GL_ARRAY_BUFFER_ARB, data->lines_count * sizeof (GLfloat) * 6, lin, GL_DYNAMIC_DRAW_ARB);

  glGenBuffersARB (1, &data->triangles);
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);
  glBufferDataARB (GL_ARRAY_BUFFER_ARB, data->triangles_count * sizeof (GLfloat) * 27, ver, GL_DYNAMIC_DRAW_ARB);

  free (lin);
  free (ver);
  MEM_Release (&setmem);
  MEM_Release (&mapmem);
  MEM_Release (&pairmem);

  return data;
}

/* update body rendering data */
static void update_body_data (BODY *bod, BODY_DATA *data)
{
  GLfloat *ver, *v, *nor, *n, *col, *c, *lin, *l;
  double **vsr, **nsr, **lsr, **vvs, **end;

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->lines);
  lin = glMapBufferARB (GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB);

  for (lsr = data->line_sources,
       end = lsr + data->lines_count * 2,
       l = lin; lsr < end; lsr ++, l += 3)
  {
    COPY (*lsr, l);
  }

  glUnmapBufferARB (GL_ARRAY_BUFFER_ARB);

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);
  ver = glMapBufferARB (GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB);
  nor = ver + data->triangles_count * 9;
  col = nor + data->triangles_count * 9;

  for (vsr = data->vertex_sources,
       nsr = data->normal_sources,
       vvs = data->vertex_value_sources,
       end = vsr + data->triangles_count * 3,
       v = ver, n = nor, c = col; vsr < end;
       vsr ++, nsr ++, vvs ++, v += 3, n += 3, c += 3)
  {
    COPY (*vsr, v);
    COPY (*nsr, n);
    value_to_color (**vvs, c);
  }

  glUnmapBufferARB (GL_ARRAY_BUFFER_ARB);
}

/* render sphere triangles */
static void render_sphere_triangles (double *center, double radius, GLfloat color [3])
{
  glMatrixMode (GL_MODELVIEW_MATRIX);
  glPushMatrix ();
    glTranslated (center[0], center[1], center[2]);
    glColor3fv (color);
    glutSolidSphere (radius, 12, 12);
  glPopMatrix ();
}

/* render sphere points */
static void render_sphere_points (double *a, double *b, double *c)
{
  glBegin (GL_POINTS);
    glVertex3dv (a);
    glVertex3dv (b);
    glVertex3dv (c);
  glEnd ();
}

/* render body triangles */
static void render_body_triangles (BODY *bod)
{
  BODY_DATA *data;

  if (bod->rendering == NULL) bod->rendering = create_body_data (bod);

  data = bod->rendering;

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);

  glEnableClientState (GL_VERTEX_ARRAY);
  glEnableClientState (GL_NORMAL_ARRAY);
  glEnableClientState (GL_COLOR_ARRAY);

  glVertexPointer (3, GL_FLOAT, 0, 0);
  glNormalPointer (GL_FLOAT, 0, (void*) (data->triangles_count * sizeof (GLfloat) * 9));
  glColorPointer (3, GL_FLOAT, 0, (void*) (data->triangles_count * sizeof (GLfloat) * 18));

  glDrawArrays (GL_TRIANGLES, 0, data->triangles_count * 3);

  glDisableClientState (GL_VERTEX_ARRAY);
  glDisableClientState (GL_NORMAL_ARRAY);
  glDisableClientState (GL_COLOR_ARRAY);

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);

  double **sph, **end;
  GLfloat *col;

  for (sph = data->spheres, end = sph + data->spheres_count * 5, col = data->sphere_colors; sph < end; sph += 5, col += 3)
    render_sphere_triangles (sph[0], *sph[1], col);
}

/* render body lines */
static void render_body_lines (BODY *bod)
{
  BODY_DATA *data;

  if (bod->rendering == NULL) bod->rendering = create_body_data (bod);

  data = bod->rendering;

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->lines);

  glEnableClientState (GL_VERTEX_ARRAY);

  glVertexPointer (3, GL_FLOAT, 0, 0);

  glDrawArrays (GL_LINES, 0, data->lines_count * 2);

  glDisableClientState (GL_VERTEX_ARRAY);

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);

  double **sph, **end;

  for (sph = data->spheres, end = sph + data->spheres_count * 5; sph < end; sph += 5)
    render_sphere_points (sph[2], sph[3], sph[4]);
}

/* render body set */
static void render_body_set (SET *set)
{
  GLfloat color [4] = {0.0, 0.0, 0.0};
  SET *item;

  glDisable (GL_LIGHTING);
  glColor3fv (color);

  for (item = SET_First (set); item; item = SET_Next (item))
  {
    render_body_lines (item->data);
  }

  glEnable (GL_LIGHTING);
  glEnable (GL_POLYGON_OFFSET_FILL);
  glPolygonOffset (1.0, 1.0);

  for (item = SET_First (set); item; item = SET_Next (item))
  {
    render_body_triangles (item->data);
  }
}

/* update scene extents */
static void update_extents ()
{
  double e [6], extents [6];
  SET *item;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;

  for (item = SET_First (selection->set); item; item = SET_Next (item))
  {
    BODY *bod = item->data;
    BODY_DATA *data = bod->rendering;

    if (data && data->flags & HIDDEN) continue;

    SHAPE_Extents (bod->shape, e);

    if (e [0] < extents [0]) extents [0] = e [0];
    if (e [1] < extents [1]) extents [1] = e [1];
    if (e [2] < extents [2]) extents [2] = e [2];
    if (e [3] > extents [3]) extents [3] = e [3];
    if (e [4] > extents [4]) extents [4] = e [4];
    if (e [5] > extents [5]) extents [5] = e [5];
  }

  GLV_Update_Extents (extents);
}

/* update bodies */
static void update ()
{
  for (BODY *bod = domain->bod; bod; bod = bod->next) update_body_data (bod, bod->rendering);

  GLV_Redraw_All ();
}

/* one simulation step */
static void step ()
{
  SOLVER_DATA *s = MAP_Find (solvers, domain, NULL);

  if (s) SOLFEC_Run (solfec, s->kind, s->solver, domain->step);
  else SOLFEC_Run (solfec, EXPLICIT_SOLVER, NULL, domain->step); /* default and in read mode */

  update ();
}

/* run simulation */
static void run (int dummy)
{
  step ();

  if (domain->flags & DOM_RUN_ANALYSIS) glutTimerFunc (1000 * SOLFEC_Time_Skip (solfec), run, 0);
}

/* seek to specific time frame */
static void seek_to_time (char *text)
{
  double t, s, e;

  if (text)
  {
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);

    SOLFEC_Time_Limits (solfec, &s, &e);

    t = atof (text);

    if (t < s) t = s;
    if (t > e) t = e;

    SOLFEC_Seek_To (solfec, t);

    update ();
  }
}

/* set skip steps value */
static void set_skip_steps (char *text)
{
  if (text)
  {
    skip_steps = ABS (atoi (text));
    skip_steps = MAX (skip_steps, 1);
  }
}

/* menus */

static void menu_domain (int item)
{
}

static void menu_render (int item)
{
}

static void menu_tools (int item)
{
}

static void menu_kinds (int item)
{
}

static void menu_analysis (int item)
{
  switch (item)
  {
  case ANALYSIS_RUN:
    domain->flags |= DOM_RUN_ANALYSIS;
    run (0);
    glutSetMenu (menu_code [MENU_ANALYSIS]);
    glutChangeToMenuEntry (1, "stop /RETURN/", ANALYSIS_STOP);
    break;
  case ANALYSIS_STOP:
    domain->flags &= ~DOM_RUN_ANALYSIS;
    glutSetMenu (menu_code [MENU_ANALYSIS]);
    glutChangeToMenuEntry (1, "run /RETURN/", ANALYSIS_RUN);
    break;
  case ANALYSIS_STEP:
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    step ();
    break;
  case ANALYSIS_SEEKTO:
    {
      static char caption [256];
      double s, e;

      SOLFEC_Time_Limits (solfec, &s, &e);
      snprintf (caption, 256, "Seek to time from [%g, %g]", s, e);
      GLV_Read_Text (caption, seek_to_time);
    }
    break;
  case ANALYSIS_FORWARD:
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    SOLFEC_Forward (solfec, skip_steps);
    update ();
    break;
  case ANALYSIS_BACKWARD:
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    SOLFEC_Backward (solfec, skip_steps);
    update ();
    break;
  case ANALYSIS_SKIP:
    GLV_Read_Text ("FORWARD and BACKWARD skip", set_skip_steps);
    break;
  }
}

static void menu_results (int item)
{
}

/* callabacks */

int RND_Menu (char ***names, int **codes)
{
  int local [4];

  ASSERT (domain, ERR_RND_NO_DOMAIN);

  menu_name [MENU_DOMAIN] = "domain";
  menu_code [MENU_DOMAIN] = glutCreateMenu (menu_domain);
  glutAddMenuEntry ("previous /</", DOMAIN_PREVIOUS);
  glutAddMenuEntry ("next />/", DOMAIN_NEXT);

  menu_name [MENU_RENDER] = "render";
  menu_code [MENU_RENDER] = glutCreateMenu (menu_render);
  glutAddMenuEntry ("2D selection /2/", RENDER_SELECTION_2D);
  glutAddMenuEntry ("3D selection /3/", RENDER_SELECTION_3D);
  glutAddMenuEntry ("previous selection /p/", RENDER_PREVIOUS_SELECTION);

  menu_name [MENU_TOOLS] = "tools";
  menu_code [MENU_TOOLS] = glutCreateMenu (menu_tools);
  glutAddMenuEntry ("toggle transparent /t/", TOOLS_TRANSPARENT);
  glutAddMenuEntry ("hide /h/", TOOLS_HIDE);
  glutAddMenuEntry ("show all /a/", TOOLS_SHOW_ALL);
  glutAddMenuEntry ("clipping plane /c/", TOOLS_CLIPPING_PLANE);
  glutAddMenuEntry ("toggle rough mesh /r/", TOOLS_ROUGH_MESH);

  menu_name [MENU_KINDS] = "kinds of";
  menu_code [MENU_KINDS] = glutCreateMenu (menu_kinds);
  glutAddMenuEntry ("constraints", KINDS_OF_CONSTRAINTS);
  glutAddMenuEntry ("forces", KINDS_OF_FORCES);
  glutAddMenuEntry ("bodies", KINDS_OF_BODIES);
  glutAddMenuEntry ("surfaces", KINDS_OF_SURFACES);
  glutAddMenuEntry ("volumes", KINDS_OF_VOLUMES);

  menu_name [MENU_ANALYSIS] = "analysis";
  menu_code [MENU_ANALYSIS] = glutCreateMenu (menu_analysis);
  glutAddMenuEntry ("run /RETURN/", ANALYSIS_RUN);
  glutAddMenuEntry ("step /SPACE/", ANALYSIS_STEP);
  if (solfec->mode == SOLFEC_READ)
  {
    glutAddMenuEntry ("seek to /UP/", ANALYSIS_SEEKTO);
    glutAddMenuEntry ("forward /LEFT/", ANALYSIS_FORWARD);
    glutAddMenuEntry ("backward /RIGHT/", ANALYSIS_BACKWARD);
    glutAddMenuEntry ("skip /DOWN/", ANALYSIS_SKIP);
  }

  local [0] = glutCreateMenu (menu_results);
  glutAddMenuEntry ("dx", RESULTS_DX);
  glutAddMenuEntry ("dy", RESULTS_DY);
  glutAddMenuEntry ("dz", RESULTS_DZ);
  local [1] = glutCreateMenu (menu_results);
  glutAddMenuEntry ("vx", RESULTS_VX);
  glutAddMenuEntry ("vy", RESULTS_VY);
  glutAddMenuEntry ("vz", RESULTS_VZ);
  local [2] = glutCreateMenu (menu_results);
  glutAddMenuEntry ("sx", RESULTS_SX);
  glutAddMenuEntry ("sy", RESULTS_SY);
  glutAddMenuEntry ("sz", RESULTS_SZ);
  glutAddMenuEntry ("sxy", RESULTS_SXY);
  glutAddMenuEntry ("sxz", RESULTS_SXZ);
  glutAddMenuEntry ("syz", RESULTS_SYZ);
  glutAddMenuEntry ("mises", RESULTS_MISES);
  local [3] = glutCreateMenu (menu_results);
  glutAddMenuEntry ("normal", RESULTS_RN);
  glutAddMenuEntry ("tangent", RESULTS_RT);
  glutAddMenuEntry ("resultant", RESULTS_R);

  menu_name [MENU_RESULTS] = "results";
  menu_code [MENU_RESULTS] = glutCreateMenu (menu_results);
  glutAddSubMenu ("displacements", local [0]);
  glutAddSubMenu ("velocities", local [1]);
  glutAddSubMenu ("stresses", local [2]);
  glutAddSubMenu ("reactions", local [3]);

  *names = menu_name;
  *codes = menu_code;
  return MENU_LAST;
}

void RND_Init ()
{
  BODY *bod;
  SET *set;

  ASSERT (domain, ERR_RND_NO_DOMAIN);

  MEM_Init (&selsetmem, sizeof (SET), BIGCHUNK);
  set = NULL;

  for (bod = domain->bod; bod; bod = bod->next) SET_Insert (&selsetmem, &set, bod, NULL);

  selection_push (set);

  update_extents ();
}

int  RND_Idle ()
{
  return 0;
}

void RND_Quit ()
{
  for (MAP *item = MAP_First (solvers); item; item = MAP_Next (item)) free (item->data); /* free solver interfaces */
  MAP_Free (NULL, &solvers);

  lngfinalize (); /* finalize Python */
}

void RND_Render ()
{
  render_body_set (selection->set);
}

void RND_Key (int key, int x, int y)
{
  switch (key)
  {
  case 27:
    break;
  case '\r':
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    else menu_analysis (ANALYSIS_RUN);
    break;
  case ' ':
    menu_analysis (ANALYSIS_STEP);
    break;
  case '<':
    menu_domain (DOMAIN_PREVIOUS);
    break;
  case '>':
    menu_domain (DOMAIN_NEXT);
    break;
  case '2':
    menu_render (RENDER_SELECTION_2D);
    break;
  case '3':
    menu_render (RENDER_SELECTION_3D);
    break;
  case 'p':
    menu_render (RENDER_PREVIOUS_SELECTION);
    break;
  case 't':
    menu_tools (TOOLS_TRANSPARENT);
    break;
  case 'h':
    menu_tools (TOOLS_HIDE);
    break;
  case 'a':
    menu_tools (TOOLS_SHOW_ALL);
    break;
  case 'c':
    menu_tools (TOOLS_CLIPPING_PLANE);
    break;
  case 'r':
    menu_tools (TOOLS_ROUGH_MESH);
    break;
  case '=':
    break;
  case '-':
    break;
  }

}

void RND_Keyspec (int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_UP:
      menu_analysis (ANALYSIS_SEEKTO);
      break;
    case GLUT_KEY_LEFT:
      menu_analysis (ANALYSIS_BACKWARD);
      break;
    case GLUT_KEY_RIGHT:
      menu_analysis (ANALYSIS_FORWARD);
      break;
    case GLUT_KEY_DOWN:
      menu_analysis (ANALYSIS_SKIP);
      break;
  }
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

/* user callable routines */

void RND_Switch_On ()
{
  enabled = 1;
}

int  RND_Is_On ()
{
  return enabled;
}

void RND_Domain (DOM *dom)
{
  dom->next = domain;
  if (domain) domain->prev = dom;
  domain = dom;
}

void RND_Solver (DOM *dom, int kind, void *solver)
{
  SOLVER_DATA *data;
  MAP *item;

  ERRMEM (data = malloc (sizeof (SOLVER_DATA)));
  data->solver = solver;
  data->kind = kind;

  if ((item = MAP_Find (solvers, dom, NULL))) free (item->data); /* if mapped, free it */

  MAP_Insert (NULL, &solvers, dom, data, NULL); /* map to domain */
}

void RND_Free_Rendering_Data (void *ptr)
{
  BODY_DATA *data;

  data = ptr;

  free (data->vertex_sources);
  free (data->values);
  free (data->value_sources);
  free (data->spheres);
  free (data->sphere_colors);

  glDeleteBuffersARB (1, &data->triangles);
  glDeleteBuffersARB (1, &data->lines);

  free (data);
}
