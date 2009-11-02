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
#include "sol.h"
#include "set.h"
#include "alg.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "shp.h"
#include "glv.h"
#include "err.h"
#include "rnd.h"
#include "lng.h"

/* ------------------------------------------------- */
static DOM *domain = NULL; /* global domain pointer */
/* ----------------------------------------------- */

#define solfec ((SOLFEC*)domain->owner)

typedef struct solver_interface SOLVER_INTERFACE;

struct solver_interface
{
  int kind;
  short *solver;
};

static MAP *solver_interfaces = NULL; /* maps domains to solver interfaces */

#define SHAPE_ELEMENT 1978 /* an arbirary number higher then shapes defined in shp.h */

static short vieweron = 0; /* viewer flag  */

static int domain_menu,
           render_menu,
	   tools_menu,
	   kindsof_menu,
	   analysis_menu,
	   analysis_menu_items,
	   results_menu; /* menu handles */

static int analysis_skip_steps = 1; /* skip steps when rewinding forward or backward */

static MAP *selection = NULL; /* selected bodies */

enum /* rendering flags */
{
  OUTLINE = 0x01,
  SELECTION = 0x02
};

enum /* menu items */
{
  DOMAIN_NEXT,
  DOMAIN_PREV,
  RENDER_ALL,
  RENDER_SELECTION_2D,
  RENDER_SELECTION_3D,
  TOOLS_TOGGLE_VISIBLE,
  TOOLS_CLIPPING_PLANE,
  TOOLS_SWITCH_OFF,
  TOOLS_SWITCH_ON_ALL,
  TOOLS_ROUGH_MESH_ON,
  ANALYSIS_RUN,
  ANALYSIS_STOP,
  ANALYSIS_STEP,
  ANALYSIS_SEEKTO,
  ANALYSIS_FORWARD,
  ANALYSIS_BACKWARD,
  ANALYSIS_SKIP,
  KINDSOF_CONSTRAINTS,
  KINDSOF_FORCES,
  KINDSOF_BODIES,
  KINDSOF_SURFACES,
  KINDSOF_VOLUMES,
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

static enum /* mouse mode */
{
  MOUSE_MODE_NONE,
  MOUSE_SELECTION_2D_BEGIN,
  MOUSE_SELECTION_2D_END,
  MOUSE_SELECTION_3D_BEGIN,
  MOUSE_SELECTION_3D_END,
} mousemode;

static int mouse_x,
	   mouse_y; /* selection coordinates */

static int current_tool = 0;

static enum /* tool modifier flags */
{
  TOOL_FLAG_SWITCH_OFF = 1, /* mutate toggle visible into switch off */
  TOOL_FLAG_ROUGH_MESH_ON,  /* rough mesh on tool */
} tool_flag = 0;

BODY *toggle_visible = NULL;

static double clipping_factor;
static double clipping_plane [4];  /* clipping plane */
static double clipping_quad [8][3]; /* plane clipping quadrangle */

static int volumetric_map = 0; /* volumetric map flag */
static int volumetric_map_range = 1; /* number of items in the volumetric map */
static int volumetric_map_kind = 0; /* kind of item drawn */
static int volumetric_map_window = 0; /* map legend window */
static int volumetric_map_is_discrete = 0; /* dicrete map flag */
static SET *volumetric_map_discrete; /* set of discrete data */
static double volumetric_map_min = 0; /* minimum of the volumetric map values */
static double volumetric_map_max = 1; /* maximum */
#define MAP_LEGEND_ROWS 8 /* number of rows in the map legend */
#define WIDTH_DISC 50 /* legend width for discrete data */
#define WIDTH_CONT 100 /* legend width for continuous data */
#define LEGEND_FONT GLV_FONT_8_BY_13

static int time_window; /* time window handler */
#define TIME_HEIGHT 16 /* time window height */
#define TIME_FONT GLV_FONT_8_BY_13

typedef void (*FUNC_OFF) (); /* turning off function callback */

#define FUNC_OFF_MAX 512

static FUNC_OFF off_stack [FUNC_OFF_MAX]; /* stack of off functions */
static int last_off = 0;

/* push on off stack */
static void pushoff (FUNC_OFF func)
{
  ASSERT (last_off < FUNC_OFF_MAX, ERR_RND_STACK_OVERFLOW);
  off_stack [last_off ++] = func;
}

/* pop from off stack */
static void popoff ()
{
  if (last_off > 0) off_stack [-- last_off] ();
}

/* find function on stack, call it and remove */
static void deloff (FUNC_OFF func)
{
  int n;

  for (n = 0; n < last_off; n ++)
  {
    if (off_stack [n] == func)
    {
      func ();

      for (; n + 1 < last_off; n ++)
	off_stack [n] = off_stack [n + 1];

      last_off --;
    }
  }
}

/* Hue to RBG color mapping */
static GLfloat Hue_2_RGB (double v1, double  v2, double vH)
{
  if (vH < 0.0) vH += 1.0;
  if (vH > 1.0) vH -= 1.0;
  if ((6.0 * vH) < 1.0) return (v1 + (v2 - v1) * 6.0 * vH);
  if ((2.0 * vH) < 1.0) return (v2);
  if ((3.0 * vH) < 2.0) return (v1 + (v2 - v1) * ((2.0/3.0) - vH) * 6.0);
  return v1;
}

/* HSL to RBG color mapping */
static void HSL_2_RGB (double H, double S, double L, GLfloat *RGB)
{
  if (S == 0)
  {
    RGB [0] = L;
    RGB [1] = L;
    RGB [2] = L;
  }
  else
  {
    double var_1, var_2;

    if (L < 0.5) var_2 = L * (1.0 + S);
    else var_2 = (L + S) - (S * L);

    var_1 = 2.0 * L - var_2;

    RGB [0] = Hue_2_RGB (var_1, var_2, H + (1.0/3.0));
    RGB [1] = Hue_2_RGB (var_1, var_2, H);
    RGB [2] = Hue_2_RGB (var_1, var_2, H - (1.0/3.0));
  }
}

/* value to color mapping */
inline static void mapcolor (double min, double val, double max, GLfloat *color)
{
  HSL_2_RGB (0.69 * (1.0-(val-min)/(max-min)), 1.0, 0.45, color);
}

/* switch off tools */
static void tools_off ()
{
  current_tool = 0;
  tool_flag = 0;
  toggle_visible = NULL;
  GLV_Redraw_All ();
}

/* time window width */
static int time_width ()
{
  return GLV_Print_Width (TIME_FONT, "t=%g", domain->time) + 6;
}

/* render current time */
static void render_time ()
{
  glDisable (GL_LIGHTING);
  glDisable (GL_DEPTH_TEST);

  glColor3f (1, 1, 1);
  glRecti (0, 0, time_width (), TIME_HEIGHT);
  glColor3f (0, 0, 0);
  GLV_Print (3, 3, 0, TIME_FONT, "t=%g", domain->time);

  glEnable (GL_LIGHTING);
  glEnable (GL_DEPTH_TEST);
}

/* get caption text */
static char* map_legend_caption ()
{
  switch (volumetric_map_kind)
  {
  case KINDSOF_CONSTRAINTS: return "CONSTRAINT KINDS";
  case KINDSOF_FORCES: return "FORCE KINDS";
  case KINDSOF_BODIES: return "BODY KINDS";
  case KINDSOF_SURFACES: return "SURFACE KINDS";
  case KINDSOF_VOLUMES: return "VOLUME KINDS";
  case RESULTS_DX: return "DX";
  case RESULTS_DY: return "DY";
  case RESULTS_DZ: return "DZ";
  case RESULTS_VX: return "VX";
  case RESULTS_VY: return "VY";
  case RESULTS_VZ: return "VZ";
  case RESULTS_SX: return "SX";
  case RESULTS_SY: return "SY";
  case RESULTS_SZ: return "SZ";
  case RESULTS_SXY: return "SXY";
  case RESULTS_SXZ: return "SXZ";
  case RESULTS_SYZ: return "SYZ";
  case RESULTS_MISES: return "MISES";
  case RESULTS_RN: return "RN";
  case RESULTS_RT: return "RT";
  case RESULTS_R: return "R";
  }

  return NULL;
}

/* get caption widtth */
static int map_caption_width ()
{
  return GLV_Print_Width (LEGEND_FONT, map_legend_caption ()) + 9;
}

/* get legend width */
static int map_legend_width ()
{
  int i, j;

  i = volumetric_map_range / MAP_LEGEND_ROWS;
  if (volumetric_map_range % MAP_LEGEND_ROWS) i ++;

  if (volumetric_map_is_discrete) i *= WIDTH_DISC;
  else i *= WIDTH_CONT;

  j = map_caption_width ();

  return MAX (i, j);
}

/* get legend height */
static int map_legend_height ()
{
  int i;

  i = MIN (volumetric_map_range, MAP_LEGEND_ROWS) + 1;

  return MAX (i, 1) * 16 + 6;
}

/* get map legend value string (if any) */
char *get_map_legend_value_string (void *data)
{
  if (volumetric_map_kind == KINDSOF_CONSTRAINTS)
  {
    switch ((int)data)
    {
      case CONTACT: return "CNT";
      case FIXPNT: return "PNT";
      case FIXDIR: return "DIR";
      case VELODIR: return "VEL";
      case RIGLNK: return "LNK";
      default: return "???";
    }
  }
  else if (volumetric_map_kind == KINDSOF_BODIES)
  {
    switch ((int)data)
    {
      case OBS: return "OBS";
      case RIG: return "RIG";
      case PRB: return "PRB";
      case FEM: return "FEM";
      default: return "???";
    }
  }

  return NULL;
}

/* render map legend */
static void render_map_legend ()
{
  GLfloat color [4] = {1, 1, 1, 1};
  double value, step;
  int i, j, k, l;
  GLint v [4];
  char *str;
  SET *item;

  glDisable (GL_LIGHTING);
  glDisable (GL_DEPTH_TEST);

  glGetIntegerv (GL_VIEWPORT, v);
  glColor3f (1, 1, 1);
  glRecti (v[0] + 3, v[1] + 3, 3 + map_caption_width (), 19);
  glColor3f (0, 0, 0);
  GLV_Print (v[0] + 6, v[1] + 6, 0, LEGEND_FONT, "%s", map_legend_caption ());

  if (volumetric_map_is_discrete)
  {
    glPushMatrix ();
    glTranslated (3, 3, 0);
    for (i = 1, j = 0, item = SET_First (volumetric_map_discrete); item; item = SET_Next (item))
    {
      mapcolor (volumetric_map_min, (double)(int)item->data, volumetric_map_max, color);
      glColor4fv (color);
      glRecti (v[0] + j * WIDTH_DISC, v[1] + i * 16, v[0] + j * WIDTH_DISC + 16, v[1] + i * 16 + 16);
      if ((str = get_map_legend_value_string (item->data)))
      {
	glColor3f (1, 1, 1);
	l = GLV_Print_Width (LEGEND_FONT, str) + 5;
	glRecti (v[0] + j * WIDTH_DISC + 16, v[1] + i * 16, v[0] + j * WIDTH_DISC + 16 + l, v[1] + (i+1) * 16);
        glColor3f (0, 0, 0);
        GLV_Print (v[0] + j * WIDTH_DISC + 18, v[1] + i * 16 + 3, 0, LEGEND_FONT, str);
      }
      else
      {
	glColor3f (1, 1, 1);
	l = GLV_Print_Width (LEGEND_FONT, "%d", (int)item->data) + 5;
	glRecti (v[0] + j * WIDTH_DISC + 16, v[1] + i * 16, v[0] + j * WIDTH_DISC + 16 + l, v[1] + (i+1) * 16);
        glColor3f (0, 0, 0);
	GLV_Print (v[0] + j * WIDTH_DISC + 18, v[1] + i * 16 + 3, 0, LEGEND_FONT, "%d", (int)item->data);
      }
      if (i ++ == MAP_LEGEND_ROWS) { i = 1; j ++; }
    }
    glPopMatrix ();
  }
  else
  {
    if (volumetric_map_min < DBL_MAX)
    {
      step = (volumetric_map_max - volumetric_map_min) / (double) volumetric_map_range;
      value = volumetric_map_min + 0.5 * step;

      glPushMatrix ();
      glTranslated (3, 3, 0);
      for (i = 1, j = k = 0; k < volumetric_map_range; k ++, value += step)
      {
	mapcolor (volumetric_map_min, value, volumetric_map_max, color);
	glColor4fv (color);
	glRecti (v[0] + j * WIDTH_CONT, v[1] + i * 16, v[0] + j * WIDTH_CONT + 16, v[1] + i * 16 + 16);
	glColor3f (1, 1, 1);
	l = GLV_Print_Width (LEGEND_FONT, "%.2e", value) + 5;
	glRecti (v[0] + j * WIDTH_CONT + 16, v[1] + i * 16, v[0] + j * WIDTH_CONT + 16 + l, v[1] + (i+1) * 16);
	glColor3f (0, 0, 0);
	GLV_Print (v[0] + j * WIDTH_CONT + 18, v[1] + i * 16 + 3, 0, LEGEND_FONT, "%.2e", value);
	if (i ++ == MAP_LEGEND_ROWS) { i = 1; j ++; }
      }
      glPopMatrix ();
    }
  }

  glEnable (GL_LIGHTING);
  glEnable (GL_DEPTH_TEST);
}

/* get nodal value from displacements onwards */
inline static double get_bod_shp_gobj_node (BODY *bod, SHAPE *shp, void *gobj, int node)
{
  switch (volumetric_map_kind)
  {
  case RESULTS_DX:
  case RESULTS_DY:
  case RESULTS_DZ:
    {
      double d [3];

      SET (d, DBL_MAX);

      BODY_Nodal_Values (bod, shp, gobj, node, VALUE_DISPLACEMENT, d);

      return d [volumetric_map_kind - RESULTS_DX];
    }
    break;
  case RESULTS_VX:
  case RESULTS_VY:
  case RESULTS_VZ:
    {
      double v [3];

      SET (v, DBL_MAX);

      BODY_Nodal_Values (bod, shp, gobj, node, VALUE_VELOCITY, v);

      return v [volumetric_map_kind - RESULTS_VX];
    }
    break;
  case RESULTS_SX:
  case RESULTS_SY:
  case RESULTS_SZ:
  case RESULTS_SXY:
  case RESULTS_SXZ:
  case RESULTS_SYZ:
    {
      double s [6];

      SET6 (s, DBL_MAX);

      BODY_Nodal_Values (bod, shp, gobj, node, VALUE_STRESS, s);

      return s [volumetric_map_kind - RESULTS_SX];
    }
    break;
  case RESULTS_MISES:
    {
      double mises = DBL_MAX;

      BODY_Nodal_Values (bod, shp, gobj, node, VALUE_MISES, &mises);

      return mises;
    }
    break;
  }

  return DBL_MAX;
}

/* value of element node */
static double get_element_node_value (BODY *bod, SHAPE *shp, ELEMENT *ele, FACE *fac, int node) /* node indexing is element/face local (depending if fac == NULL) */
{
  if (fac) /* translate face node number to element node number */
  {
    int n = fac->nodes [node], i;

    for (i = 0; i < ele->type; i ++)
      if (ele->nodes [i] == n) break;

    ASSERT_DEBUG (i < ele->type, "Inconsitency in face and element node numbering");

    node = i;
  }

  switch (volumetric_map_kind)
  {
  case KINDSOF_CONSTRAINTS: return DBL_MAX;
  case KINDSOF_FORCES: return DBL_MAX;
  case KINDSOF_BODIES: return bod->kind;
  case KINDSOF_SURFACES:
    if (fac) return fac->surface;
    else return DBL_MAX;
    break;
  case KINDSOF_VOLUMES: return ele->volume;
  default: return get_bod_shp_gobj_node (bod, shp, ele, node);
  }

  return DBL_MAX;
}

/* value of convex face vertex */
static double get_convex_face_vertex_value (BODY *bod, SHAPE *shp, CONVEX *cvx, int *fac, int f, int v)
{
  int node = fac [v] / 3;

  switch (volumetric_map_kind)
  {
  case KINDSOF_CONSTRAINTS: return DBL_MAX;
  case KINDSOF_FORCES: return DBL_MAX;
  case KINDSOF_BODIES: return bod->kind;
  case KINDSOF_SURFACES: return cvx->surface [f];
  case KINDSOF_VOLUMES: return cvx->volume;
  default: return get_bod_shp_gobj_node (bod, shp, cvx, node);
  }

  return DBL_MAX;
}

/* value of sphere surface */
static double get_sphere_surface_value (BODY *bod, SHAPE *shp, SPHERE *sph)
{
  switch (volumetric_map_kind)
  {
  case KINDSOF_CONSTRAINTS: return DBL_MAX;
  case KINDSOF_FORCES: return DBL_MAX;
  case KINDSOF_BODIES: return bod->kind;
  case KINDSOF_SURFACES: return sph->surface;
  case KINDSOF_VOLUMES: return sph->volume;
  default: return get_bod_shp_gobj_node (bod, shp, sph, 0);
  }

  return DBL_MAX;
}

/* register new value */
inline static void volumetric_map_register_value (double value)
{
  if (value < DBL_MAX)
  {
    if (volumetric_map_is_discrete)
    {
      SET_Insert (&domain->setmem, &volumetric_map_discrete, (void*)(int)value, NULL);
    }

    if (value < volumetric_map_min) volumetric_map_min = value;
    if (value > volumetric_map_max) volumetric_map_max = value;
  }
}

/* process body to get extrema of volumetric map */
static void volumetric_map_extrema_process_body (BODY *bod)
{
  double value;
  ELEMENT *ele;
  CONVEX *cvx;
  SPHERE *sph;
  SHAPE *shp;
  FACE *fac;
  MESH *msh;
  int n, m, l;

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
	  for (n = 0; n < fac->type; n ++)
	  {
	    value = get_element_node_value (bod, shp, ele, fac, n);
	    volumetric_map_register_value (value);
	  }
	}
      }

      for (ele = msh->bulkeles; ele; ele = ele->next)
      {
	for (n = 0; n < ele->type; n ++)
	{
	  value = get_element_node_value (bod, shp, ele, NULL, n);
	  volumetric_map_register_value (value);
	}
      }

      break;
    case SHAPE_CONVEX:

      for (cvx = shp->data; cvx; cvx = cvx->next)
      {
	for (n = m = 0; n < cvx->nfac; n ++, m += (cvx->fac [m] + 1))
	{
	  for (l = 1; l <= cvx->fac [m]; l ++)
	  {
	    value = get_convex_face_vertex_value (bod, shp, cvx, &cvx->fac [m], n, l);
	    volumetric_map_register_value (value);
	  }
	}
      }

      break;
    case SHAPE_SPHERE:

      for (sph = shp->data; sph; sph = sph->next)
      {
	value = get_sphere_surface_value (bod, shp, sph);
	volumetric_map_register_value (value);
      }

      break;
    }
  }
}

/* skip this constraint when calculating
 * extrema or rendering constraints ? */
static int skip_constraint (CON *con)
{
  if (selection)
  {
    void *one = MAP_Find (selection, (void*)con->master->id, NULL),
	 *two = con->slave ? MAP_Find (selection, (void*)con->slave->id, NULL) : NULL;

    if ((!(one || two)) ||
	((!one) && (con->slave && (con->slave->flags & (BODY_HIDDEN|BODY_OFF)))) ||
	((!two) && con->master->flags & (BODY_HIDDEN|BODY_OFF))) return 1;

  }

  if ((con->master->flags & (BODY_HIDDEN|BODY_OFF)) &&
      (con->slave == NULL || (con->slave->flags & (BODY_HIDDEN|BODY_OFF)))) return 1;

  return 0;
}


/* shall we draw on constraints ? */
static int volumetric_map_constraint_based ()
{
  if (!volumetric_map) return 0;

  switch (volumetric_map_kind)
  {
  case KINDSOF_CONSTRAINTS:
  case RESULTS_RT:
  case RESULTS_RN:
  case RESULTS_R:
    return 1;
  }

  return 0;
}

/* render mesh surface */
/* set map extrema */
static void volumetric_map_extrema ()
{
  if (volumetric_map_is_discrete)
  {
    volumetric_map_min = 0;
    volumetric_map_max = 1;
  }
  else
  {
    volumetric_map_min = DBL_MAX;
    volumetric_map_max = -DBL_MAX;
  }


  if (volumetric_map_constraint_based ())
  {
    CON *con;

    for (con = domain->con; con; con = con->next)
    {
      if (skip_constraint (con)) continue;

      switch (volumetric_map_kind)
      {
      case KINDSOF_CONSTRAINTS:
	volumetric_map_register_value (con->kind);
	break;
      case RESULTS_RT:
	volumetric_map_register_value (LEN2 (con->R));
	break;
      case RESULTS_RN:
	volumetric_map_register_value (con->R [2]);
	break;
      case RESULTS_R:
	volumetric_map_register_value (LEN (con->R));
	break;
      }
    }
  }
  else
  {
    BODY *bod;

    if (selection)
    {
      for (MAP *item = MAP_First (selection); item; item = MAP_Next (item))
      {
	bod = item->data;
	if (bod->flags & (BODY_HIDDEN|BODY_OFF)) continue;
	volumetric_map_extrema_process_body (bod);
      }
    }
    else
    {
      for (bod = domain->bod; bod; bod = bod->next)
      {
	if (bod->flags & (BODY_HIDDEN|BODY_OFF)) continue;
	volumetric_map_extrema_process_body (bod);
      }
    }
  }

 if (volumetric_map_min == volumetric_map_max)
 {
   double d = 1E-9 * fabs (volumetric_map_min);
   volumetric_map_min -= d;
   volumetric_map_max += d;
 }
}

/* turn off volumetric map */
static void volumetric_map_off ()
{
  volumetric_map = 0;
  GLV_Close_Viewport (volumetric_map_window);
}

/* turn on volumetric map */
static void volumetric_map_on (int kind)
{
  if (!volumetric_map) pushoff (volumetric_map_off); /* push switch off on stack (ESCAPE pops those functions) */

  volumetric_map = 1;
  volumetric_map_kind = kind;

  switch (kind)
  {
  case KINDSOF_CONSTRAINTS:
  case KINDSOF_FORCES:
  case KINDSOF_BODIES:
  case KINDSOF_SURFACES:
  case KINDSOF_VOLUMES:
    SET_Free (&domain->setmem, &volumetric_map_discrete);
    volumetric_map_is_discrete = 1;
    break;
  default:
    volumetric_map_is_discrete = 0;
    break;
  }

  volumetric_map_extrema ();

  if (volumetric_map_is_discrete)
    volumetric_map_range = SET_Size (volumetric_map_discrete);
  else volumetric_map_range = 8;

  if (volumetric_map_window) GLV_Close_Viewport (volumetric_map_window);
  volumetric_map_window = GLV_Open_Viewport (0, 0, map_legend_width (), map_legend_height (), 0, render_map_legend);
}

/* color of element node */
static int get_element_node_color (BODY *bod, SHAPE *shp, ELEMENT *ele, FACE *fac, int node, GLfloat *color)
{
  double val;

  if (!volumetric_map) return 0;
  val = get_element_node_value (bod, shp, ele, fac, node);
  if (val < DBL_MAX)
  {
    mapcolor (volumetric_map_min, val, volumetric_map_max, color);
    return 1;
  }
  return 0;
}

/* color of convex face vertex */
static int get_convex_face_vertex_color (BODY *bod, SHAPE *shp, CONVEX *cvx, int *fac, int f, int v, GLfloat *color)
{
  double val;

  if (!volumetric_map) return 0;
  val = get_convex_face_vertex_value (bod, shp, cvx, fac, f, v);

  if (val < DBL_MAX)
  {
    mapcolor (volumetric_map_min, val, volumetric_map_max, color);
    return 1;
  }
  return 0;
}

/* color of sphere surface */
static void get_sphere_surface_color (BODY *bod, SHAPE *shp, SPHERE *sph, GLfloat *color)
{
  double val;

  if (!volumetric_map) return;
  val = get_sphere_surface_value (bod, shp, sph);
  if (val < DBL_MAX) mapcolor (volumetric_map_min, val, volumetric_map_max, color);
}

/* color of constraint based rendering */
static void get_constraint_based_color (CON *con, GLfloat *color)
{
  double value = volumetric_map_min;

  switch (volumetric_map_kind)
  {
  case KINDSOF_CONSTRAINTS:
    value = con->kind;
    break;
  case RESULTS_RT:
    value = LEN2 (con->R);
    break;
  case RESULTS_RN:
    value = con->R [2];
    break;
  case RESULTS_R:
    value = LEN (con->R);
    break;
  }

  if (value < DBL_MAX) mapcolor (volumetric_map_min, value, volumetric_map_max, color);
}

/* render mesh surface */
static void render_mesh (BODY *bod, SHAPE *shp, MESH *mesh, int flags, GLfloat color [4])
{
  GLfloat outcol [4] = {0.0, 0.0, 0.0, 1.0};
  double (*cur) [3] = mesh->cur_nodes;
  ELEMENT *ele;
  FACE *fac;
  GLboolean blend;

  glGetBooleanv (GL_BLEND, &blend);
  if (blend)
  { SET (outcol, 0.65); }

  for (ele = mesh->surfeles; ele; ele = ele->next)
  {
    for (fac = ele->faces; fac; fac = fac->next)
    {
      /* outline */
      if (flags & OUTLINE)
      {
	glDisable (GL_LIGHTING);
	glColor4fv (outcol);
	if (fac->type == 3)
	{
	  glBegin (GL_LINE_LOOP);
	  glVertex3dv (cur[fac->nodes[0]]);
	  glVertex3dv (cur[fac->nodes[1]]);
	  glVertex3dv (cur[fac->nodes[2]]);
	  glEnd ();
	}
	else
	{
	  glBegin (GL_LINE_LOOP);
	  glVertex3dv (cur[fac->nodes[0]]);
	  glVertex3dv (cur[fac->nodes[1]]);
	  glVertex3dv (cur[fac->nodes[2]]);
	  glVertex3dv (cur[fac->nodes[3]]);
	  glEnd ();
	}
	glEnable (GL_LIGHTING);

	glEnable (GL_POLYGON_OFFSET_FILL);
	glPolygonOffset (1.0, 1.0);
      }

      /* fill */
      glColor4fv (color);
      if (fac->type == 3)
      {
	if (flags & SELECTION)
	{
          glBegin (GL_TRIANGLES);
	  glNormal3dv (fac->normal);
	  glVertex3dv (cur[fac->nodes[0]]);
	  glVertex3dv (cur[fac->nodes[1]]);
	  glVertex3dv (cur[fac->nodes[2]]);
	  glEnd ();
	}
	else
	{
	  glBegin (GL_TRIANGLES);
	  glNormal3dv (fac->normal);
	  if (shp && get_element_node_color (bod, shp, ele, fac, 0, color)) glColor4fv (color);
	  glVertex3dv (cur[fac->nodes[0]]);
	  if (shp && get_element_node_color (bod, shp, ele, fac, 1, color)) glColor4fv (color);
	  glVertex3dv (cur[fac->nodes[1]]);
	  if (shp && get_element_node_color (bod, shp, ele, fac, 2, color)) glColor4fv (color);
	  glVertex3dv (cur[fac->nodes[2]]);
	  glEnd ();
	}
      }
      else
      {
	if (flags & SELECTION)
	{
	  glBegin (GL_QUADS);
	  glNormal3dv (fac->normal);
	  glVertex3dv (cur[fac->nodes[0]]);
	  glVertex3dv (cur[fac->nodes[1]]);
	  glVertex3dv (cur[fac->nodes[2]]);
	  glVertex3dv (cur[fac->nodes[3]]);
	  glEnd ();
	}
	else
	{
	  glBegin (GL_QUADS);
	  glNormal3dv (fac->normal);
	  if (shp && get_element_node_color (bod, shp, ele, fac, 0, color)) glColor4fv (color);
	  glVertex3dv (cur[fac->nodes[0]]);
	  if (shp && get_element_node_color (bod, shp, ele, fac, 1, color)) glColor4fv (color);
	  glVertex3dv (cur[fac->nodes[1]]);
	  if (shp && get_element_node_color (bod, shp, ele, fac, 2, color)) glColor4fv (color);
	  glVertex3dv (cur[fac->nodes[2]]);
	  if (shp && get_element_node_color (bod, shp, ele, fac, 3, color)) glColor4fv (color);
	  glVertex3dv (cur[fac->nodes[3]]);
	  glEnd ();
	}
      }

      if (flags & OUTLINE) glDisable (GL_POLYGON_OFFSET_FILL);
    }
  }
}

/* render single convex */
static void render_convex (BODY *bod, SHAPE *shp, CONVEX *cvx, int flags, GLfloat color [4])
{
  GLfloat outcol [4] = {0.0, 0.0, 0.0, 1.0};
  int n, m, l;
  GLboolean blend;

  glGetBooleanv (GL_BLEND, &blend);
  if (blend)
  { SET (outcol, 0.65); }

  /* outline */
  if (flags & OUTLINE)
  {
    glDisable (GL_LIGHTING);
    glColor4fv (outcol);
    for (n = m = 0; n < cvx->nfac;
      n ++, m += (cvx->fac [m] + 1))
    {
      glBegin (GL_LINE_LOOP);
      for (l = 1; l <= cvx->fac [m]; l ++)
	glVertex3dv (&cvx->cur [cvx->fac [m + l]]);
      glEnd ();
    }
    glEnable (GL_LIGHTING);

    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1.0, 1.0);
  }

  /* fill */
  glColor4fv (color);
  for (n = m = 0; n < cvx->nfac;
    n ++, m += (cvx->fac [m] + 1))
  {
    glBegin (GL_POLYGON);
    glNormal3dv (&cvx->pla [n * 4]);
    for (l = 1; l <= cvx->fac [m]; l ++)
    {
      if (!(flags & SELECTION) && get_convex_face_vertex_color (bod, shp, cvx, &cvx->fac [m], n, l, color)) glColor4fv (color);
      glVertex3dv (&cvx->cur [cvx->fac [m + l]]);
    }
    glEnd ();
  }

  if (flags & OUTLINE) glDisable (GL_POLYGON_OFFSET_FILL);
}

/* render single sphere */
static void render_sphere (BODY *bod, SHAPE *shp, SPHERE *sphere, int flags, GLfloat color [4])
{
  GLfloat outcol [4] = {0.0, 0.0, 0.0, 1.0};
  double *c = sphere->cur_center;
  GLboolean blend;

  glGetBooleanv (GL_BLEND, &blend);
  if (blend)
  { SET (outcol, 0.65); }

  if (flags & OUTLINE)
  {
    glPointSize (2.0);
    glColor4fv (outcol);
    glBegin (GL_POINTS);
      glVertex3dv (sphere->cur_points[0]);
      glVertex3dv (sphere->cur_points[1]);
      glVertex3dv (sphere->cur_points[2]);
    glEnd ();

    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1.0, 1.0);
    glPointSize (1.0);
  }

  glMatrixMode (GL_MODELVIEW_MATRIX);
  glPushMatrix ();
    glTranslated (c[0], c[1], c[2]);
    if (!(flags & SELECTION)) get_sphere_surface_color (bod, shp, sphere, color);
    glColor4fv (color);
    glutSolidSphere (sphere->cur_radius, 12, 12);
  glPopMatrix ();

  if (flags & OUTLINE) glDisable (GL_POLYGON_OFFSET_FILL);
}

/* render single element */
static void render_element (BODY *bod, SHAPE *shp, ELEMENT *ele, int flags, GLfloat color [4])
{
  GLfloat outcol [4] = {0.0, 0.0, 0.0, 1.0};
  CONVEX *cvx;
  int n, m, l;

  /* get convex of element */
  cvx = ELEMENT_Convex (shp->data, ele);

  /* outline */
  if (flags & OUTLINE)
  {
    glDisable (GL_LIGHTING);
    glColor4fv (outcol);
    for (n = m = 0; n < cvx->nfac;
      n ++, m += (cvx->fac [m] + 1))
    {
      glBegin (GL_LINE_LOOP);
      for (l = 1; l <= cvx->fac [m]; l ++)
	glVertex3dv (&cvx->cur [cvx->fac [m + l]]);
      glEnd ();
    }
    glEnable (GL_LIGHTING);

    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1.0, 1.0);
  }

  /* fill */
  glColor4fv (color);
  for (n = m = 0; n < cvx->nfac;
    n ++, m += (cvx->fac [m] + 1))
  {
    glBegin (GL_POLYGON);
    glNormal3dv (&cvx->pla [n * 4]);
    for (l = 1; l <= cvx->fac [m]; l ++)
    {
      if (!(flags & SELECTION) && get_element_node_color (bod, shp, ele, NULL, cvx->fac [m+l] / 3, color)) glColor4fv (color);
      glVertex3dv (&cvx->cur [cvx->fac [m+l]]);
    }
    glEnd ();
  }

  if (flags & OUTLINE) glDisable (GL_POLYGON_OFFSET_FILL);

  CONVEX_Destroy (cvx);
}

/* render contact */
static void render_contact (CON *con, GLfloat color [4])
{
  double other [3],
	 scal = GLV_Minimal_Extent() * 0.03;

  glColor4fv (color);

  glPointSize (3.0);
  glBegin (GL_POINTS);
    glVertex3dv (con->point);
  glEnd ();
  glPointSize (1.0);

  ADDMUL (con->point, 0.5*scal, con->base, other);

  glLineWidth (2.0);
  glBegin (GL_LINES);
    glVertex3dv (con->point);
    glVertex3dv (other);
  glEnd ();
  glLineWidth (1.0);

  ADDMUL (con->point, 0.5*scal, con->base+3, other);

  glLineWidth (2.0);
  glBegin (GL_LINES);
    glVertex3dv (con->point);
    glVertex3dv (other);
  glEnd ();
  glLineWidth (1.0);

  ADDMUL (con->point, scal, con->base+6, other);

  glLineWidth (2.0);
  glBegin (GL_LINES);
    glVertex3dv (con->point);
    glVertex3dv (other);
  glEnd ();
  glLineWidth (1.0);
}

/* render fixed point */
static void render_fixpnt (CON *con, GLfloat color [4])
{
  glColor4fv (color);
  glPointSize (4.0);
  glBegin (GL_POINTS);
    glVertex3dv (con->point);
  glEnd ();
  glPointSize (1.0);
}

/* render fixed direction */
static void render_fixdir (CON *con, GLfloat color [4])
{
  double other [3],
	 scal = GLV_Minimal_Extent() * 0.03;

  glColor4fv (color);

  glPointSize (3.0);
  glBegin (GL_POINTS);
    glVertex3dv (con->point);
  glEnd ();
  glPointSize (1.0);

  ADDMUL (con->point, scal, con->base+6, other);

  glLineWidth (2.0);
  glBegin (GL_LINES);
    glVertex3dv (con->point);
    glVertex3dv (other);
  glEnd ();
  glLineWidth (1.0);
}

/* render velocity constraint */
static void render_velodir (CON *con, GLfloat color [4])
{
  double other [3],
	 scal = GLV_Minimal_Extent() * 0.03;

  glColor4fv (color);

  glPointSize (3.0);
  glBegin (GL_POINTS);
    glVertex3dv (con->point);
  glEnd ();
  glPointSize (1.0);

  ADDMUL (con->point, scal, con->base+6, other);

  glLineWidth (2.0);
  glBegin (GL_LINES);
    glVertex3dv (con->point);
    glVertex3dv (other);
  glEnd ();
  glLineWidth (1.0);
}

/* render rigid link constraint */
static void render_riglnk (CON *con, GLfloat color [4])
{
  double other [3];

  ADD (con->point, RIGLNK_VEC (con->Z), other);

  glColor4fv (color);

  glPointSize (3.0);
  glBegin (GL_POINTS);
    glVertex3dv (con->point);
    glVertex3dv (other);
  glEnd ();
  glPointSize (1.0);
  
  glLineWidth (2.0);
  glBegin (GL_LINES);
    glVertex3dv (con->point);
    glVertex3dv (other);
  glEnd ();
  glLineWidth (1.0);
}

/* draw arrow from p to q */
static void arrow3d (double *p, double *q)
{
  double R [9],
	 a [3],
	 b [3],
	 c [3],
	 d [3],
	 x [3],
	 y [3],
	 t [3],
	 r [3],
	 n [3],
	 l,
	 o,
	 s,
	 angle;
  int i;

  SUB (q, p, d);
  MAXABSIDX (d, i);
  SET (x, 1.0);
  ADD (x, d, x);
  x [i] = 0.0;
  PRODUCT (d, x, t);
  l = LEN (d);
  DIV (d, l, d);
  COPY (d, n);
  SCALE (n, -1);
  NORMALIZE (t);
  o = 0.65 * l;
  ADDMUL (p, o, d, r);
  angle = ALG_PI / 8.0;
  SCALE (d, angle);
  EXPMAP (d, R);
  o = 0.075 * l;
  s = 0.150 * l;

  glDisable (GL_LIGHTING);

  ADDMUL (p, o, t, a);
  ADDMUL (r, o, t, b);

  for (i = 0; i < 16; i ++)
  {
    SUB (b, r, x);
    NVMUL (R, x, y);
    ADD (r, y, c);

    SUB (a, r, x);
    NVMUL (R, x, y);
    ADD (r, y, d);

    glNormal3dv (y);
    glBegin (GL_QUADS);
    glVertex3dv (d);
    glVertex3dv (c);
    glVertex3dv (b);
    glVertex3dv (a);
    glEnd ();

    COPY (c, b);
    COPY (d, a);
  }

  ADDMUL (r, o, t, a);
  ADDMUL (r, s, t, b);
  for (i = 0; i < 16; i ++)
  {
    SUB (b, r, x);
    NVMUL (R, x, y);
    ADD (r, y, c);

    SUB (a, r, x);
    NVMUL (R, x, y);
    ADD (r, y, d);

    glNormal3dv (n);
    glBegin (GL_QUADS);
    glVertex3dv (d);
    glVertex3dv (c);
    glVertex3dv (b);
    glVertex3dv (a);
    glEnd ();

    NORMAL (b, c, a, x);
    glNormal3dv (x);
    glBegin (GL_TRIANGLES);
    glVertex3dv (b);
    glVertex3dv (c);
    glVertex3dv (q);
    glEnd ();

    COPY (c, b);
    COPY (d, a);
  }

  ADDMUL (p, o, t, a);
  glNormal3dv (n);
  glBegin (GL_POLYGON);
  glVertex3dv (a);
  for (i = 0; i < 16; i ++)
  {
    SUB (a, p, x);
    TVMUL (R, x, y);
    ADD (p, y, b);
    glVertex3dv (b);
    COPY (b, a);
  }
  glEnd ();


  glEnable (GL_LIGHTING);
}

/* render tangential reactions */
static void render_rt (CON *con, GLfloat color [4])
{
  double r [3],
	 other [3],
	 ext = GLV_Minimal_Extent() * 0.1,
	 eps, len;

  COPY (con->base, r);
  SCALE (r, con->R[0]);
  ADDMUL (r, con->R[1], con->base+3, r);
  len = LEN (r);
  eps = (ext  / len) * (1.0 + (len - volumetric_map_min) / (volumetric_map_max - volumetric_map_min + 1.0));
  ADDMUL (con->point, -eps, r, other);

  glColor4fv (color);
  arrow3d (other, con->point);
}

/* render normal reactions */
static void render_rn (CON *con, GLfloat color [4])
{
  double r [3],
	 other [3],
	 ext = GLV_Minimal_Extent() * 0.1,
	 eps, len;

  COPY (con->base + 6, r);
  SCALE (r, con->R[2]);
  len = LEN (r);
  eps = (ext  / len) * (1.0 + (len - volumetric_map_min) / (volumetric_map_max - volumetric_map_min + 1.0));
  ADDMUL (con->point, -eps, r, other);

  glColor4fv (color);
  arrow3d (other, con->point);
}

/* render rigid link constraint */
static void render_r (CON *con, GLfloat color [4])
{
  double r [3],
	 other [3],
	 ext = GLV_Minimal_Extent() * 0.1,
	 eps,
	 len;

  COPY (con->base, r);
  SCALE (r, con->R[0]);
  ADDMUL (r, con->R[1], con->base+3, r);
  ADDMUL (r, con->R[2], con->base+6, r);
  len = LEN (r);
  eps = (ext  / len) * (1.0 + (len - volumetric_map_min) / (volumetric_map_max - volumetric_map_min + 1.0));
  ADDMUL (con->point, -eps, r, other);

  glColor4fv (color);
  arrow3d (other, con->point);
}

/* compute scene extents */
static void get_scene_extents (double *extents)
{
  double e [6];
  BODY *bod;
  CON *con;

  extents [0] = extents [1] = extents [2] =  DBL_MAX;
  extents [3] = extents [4] = extents [5] = -DBL_MAX;

  if (selection)
  {
    for (MAP *item = MAP_First (selection); item; item = MAP_Next (item))
    {
      BODY *bod = item->data;

      if (bod->flags & (BODY_HIDDEN|BODY_OFF)) continue;

      SHAPE_Extents (bod->shape, e);

      if (e [0] < extents [0]) extents [0] = e [0];
      if (e [1] < extents [1]) extents [1] = e [1];
      if (e [2] < extents [2]) extents [2] = e [2];
      if (e [3] > extents [3]) extents [3] = e [3];
      if (e [4] > extents [4]) extents [4] = e [4];
      if (e [5] > extents [5]) extents [5] = e [5];
    }
  }
  else
  {
    for (bod = domain->bod; bod; bod = bod->next)
    {
      if (bod->flags & (BODY_HIDDEN|BODY_OFF)) continue;

      SHAPE_Extents (bod->shape, e);

      if (e [0] < extents [0]) extents [0] = e [0];
      if (e [1] < extents [1]) extents [1] = e [1];
      if (e [2] < extents [2]) extents [2] = e [2];
      if (e [3] > extents [3]) extents [3] = e [3];
      if (e [4] > extents [4]) extents [4] = e [4];
      if (e [5] > extents [5]) extents [5] = e [5];
    }
  }

  for (con = domain->con; con; con = con->next)
  {
    if (skip_constraint (con)) continue;

    if (con->kind == RIGLNK)
    {
      if (con->slave == NULL)
      {
	ADD (con->point, RIGLNK_VEC (con->Z), e);

	if (e [0] < extents [0]) extents [0] = e [0];
	if (e [1] < extents [1]) extents [1] = e [1];
	if (e [2] < extents [2]) extents [2] = e [2];
	if (e [0] > extents [3]) extents [3] = e [0];
	if (e [1] > extents [4]) extents [4] = e [1];
	if (e [2] > extents [5]) extents [5] = e [2];
      }
    }
  }
}

/* update viwer extents */
static void update_extents ()
{
  double extents [6];

  get_scene_extents (extents);

  GLV_Update_Extents (extents);
}

/* clipping reaction to mouse motion */
static void move_clipping_plane (int x, int y)
{
  GLint viewport [4];
  double u;
  int xmid, ymid;
  double *a = clipping_quad [0],
	 *b = clipping_quad [1],
	 *c = clipping_quad [2],
	 *d = clipping_quad [3],
         *a0 = clipping_quad [4],
	 *b0 = clipping_quad [5],
	 *c0 = clipping_quad [6],
	 *d0 = clipping_quad [7];


  glGetIntegerv (GL_VIEWPORT, viewport);

  xmid = viewport [2] / 2;
  ymid = viewport [3] / 2;

  u = -clipping_factor * (double)(y - ymid) / (double) MIN (viewport[2], viewport[3]);

  ADDMUL (a0, u, clipping_plane, a);
  ADDMUL (b0, u, clipping_plane, b);
  ADDMUL (c0, u, clipping_plane, c);
  ADDMUL (d0, u, clipping_plane, d);

  clipping_plane [3] = -DOT (clipping_plane, a);
}

/* calculate clipping plane */
static void update_clipping_plane (double *normal)
{
  double extents [6],
	 *a = clipping_quad [0],
	 *b = clipping_quad [1],
	 *c = clipping_quad [2],
	 *d = clipping_quad [3],
	 tang [2][3],
	 mid [3],
	 width [3], w;

  get_scene_extents (extents);

  w = DOT (normal, extents + 3);
  clipping_factor = DOT (normal, extents) - w;
  clipping_factor = fabs (clipping_factor); /* extents along the normal */

  COPY (normal, tang [0]);
  tang [0][0] += DRANDEXT (-1, 1);
  tang [0][1] += DRANDEXT (-1, 1);
  tang [0][2] += DRANDEXT (-1, 1);
  PRODUCT (normal, tang [0], tang [1]);
  PRODUCT (normal, tang [1], tang [0]);
  NORMALIZE (tang [0]);
  NORMALIZE (tang [1]);

  MID (extents, extents + 3, mid);
  SUB (extents + 3, extents, width);
  w = LEN (width);

  COPY (mid, a);
  ADDMUL (a, -w, tang[0], a);
  ADDMUL (a, -w, tang[1], a);

  COPY (mid, b);
  ADDMUL (b, -w, tang[0], b);
  ADDMUL (b,  w, tang[1], b);

  COPY (mid, c);
  ADDMUL (c, w, tang[0], c);
  ADDMUL (c, w, tang[1], c);

  COPY (mid, d);
  ADDMUL (d,  w, tang[0], d);
  ADDMUL (d, -w, tang[1], d);

  COPY (a, clipping_quad [4]);
  COPY (b, clipping_quad [5]);
  COPY (c, clipping_quad [6]);
  COPY (d, clipping_quad [7]);
}

/* read clipping normal */
static void read_clipping_normal (char *text)
{
  double normal [3];
  int i, j, k, n;

  if (text)
  {
    n = strlen (text);

    for (j = k = 0; k < 3; k ++, j = i)
    {
      for (i = j; i < n; i ++)
      {
	if (text [i] == ',')
	{
	  text [i] = '\0';
	  i ++;
	  break;
	}
      }
      normal [k] = atof (&text [j]);
    }
  }

  if (LEN (normal) > 0.0)
  {
    NORMALIZE (normal);
    COPY (normal, clipping_plane);
    update_clipping_plane (normal);
    if (!current_tool) pushoff (tools_off);
    current_tool = TOOLS_CLIPPING_PLANE;
    GLV_Redraw_All ();
  }
}

/* declare menu callback */
static void menu_analysis (int);

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

    if (volumetric_map) volumetric_map_on (volumetric_map_kind);
    else GLV_Redraw_All ();

    /* resize time viewport so to fit the text output */
    GLV_Resize_Viewport (time_window, time_width (), TIME_HEIGHT);
  }
}

/* set skip steps value */
static void set_skip_steps (char *text)
{
  if (text)
  {
    analysis_skip_steps = ABS (atoi (text));
    analysis_skip_steps = MAX (analysis_skip_steps, 1);
  }
}

/* set up body rendering => return 1 if this body should be skipped or 0 otherwise */
inline static int render_body_begin (BODY *bod, GLfloat color [4])
{
  if (bod->flags & BODY_OFF) return 1;
  else if (bod == toggle_visible)
  {
    if (!bod->flags & BODY_HIDDEN)
    {
      glEnable (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }

    if (tool_flag) color [1] = 1.0;
    else color [0] = 1.0;
  }
  else if (current_tool == TOOLS_TOGGLE_VISIBLE &&
	   bod->flags & BODY_HIDDEN)
  {
      glEnable (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  else if (bod->flags & BODY_HIDDEN) return 1;

  return 0;
}


/* finalize body rendering */
inline static void render_body_end (BODY *bod, GLfloat color [4])
{
  if (bod == toggle_visible)
  {
    if (! bod->flags & BODY_HIDDEN)
      glDisable (GL_BLEND);

    if (tool_flag) color [1] = 0.7;
    else color [0] = 0.7;
  }
  else if (current_tool == TOOLS_TOGGLE_VISIBLE &&
	   bod->flags & BODY_HIDDEN) glDisable (GL_BLEND);
}

/* render volumetric object - basic */
inline static void render_object_basic (BODY *bod, SHAPE *shp, short kind, void *obj, int flags, GLfloat color [4])
{
  switch (kind)
  {
    case SHAPE_MESH: render_mesh (bod, shp, obj, flags, color); break;
    case SHAPE_CONVEX: render_convex (bod, shp, obj, flags, color); break;
    case SHAPE_SPHERE: render_sphere (bod, shp, obj, flags, color); break;
    case SHAPE_ELEMENT: render_element (bod, shp, obj, flags, color); break;
  }
}

/* render volumetric object - complete (with clipping, etc.) */
static void render_object (BODY *bod, SHAPE *shp, short kind, void *obj, int flags, GLfloat color [4])
{
  if (current_tool == TOOLS_CLIPPING_PLANE)
  {
    double normal [3];

    glEnable (GL_CLIP_PLANE0);
    glClipPlane (GL_CLIP_PLANE0, clipping_plane);

    glDisable (GL_DEPTH_TEST);
    glColorMask (GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

    glEnable (GL_STENCIL_TEST);
    glClear (GL_STENCIL_BUFFER_BIT);

    glStencilFunc (GL_ALWAYS, 0, 0);

    glStencilOp (GL_KEEP, GL_KEEP, GL_INCR);
    glCullFace (GL_FRONT);
    render_object_basic (bod, shp, kind, obj, flags, color);

    glStencilOp (GL_KEEP, GL_KEEP, GL_DECR);
    glCullFace (GL_BACK); 
    render_object_basic (bod, shp, kind, obj, flags, color);

    glEnable (GL_DEPTH_TEST);
    glDisable (GL_CLIP_PLANE0);
    glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

    glStencilFunc (GL_NOTEQUAL, 0, ~0); 

    COPY (clipping_plane, normal);
    SCALE (normal, -1);

    glBegin (GL_QUADS);
    glNormal3dv (normal);
    for (int j = 3; j >= 0; j--) glVertex3dv (clipping_quad [j]);
    glEnd ();

    glDisable (GL_STENCIL_TEST);
    glEnable (GL_CLIP_PLANE0);
  }

  render_object_basic (bod, shp, kind, obj, flags, color);

  if (current_tool == TOOLS_CLIPPING_PLANE)
  {
    glDisable (GL_CLIP_PLANE0);
  }
}

/* render shape */
static void render_shape (BODY *bod, SHAPE *shp, int flags, GLfloat color [4])
{
  switch (shp->kind)
  {
  case SHAPE_MESH:
    if (current_tool == TOOLS_CLIPPING_PLANE && volumetric_map) /* draw bulk elements if clipping and volumetric maps are on */
    {
      ELEMENT *ele;
      MESH *msh;

      msh = shp->data;

      for (ele = msh->surfeles; ele; ele = ele->next)
      {
	render_object (bod, shp, SHAPE_ELEMENT, ele, flags, color);
      }

      for (ele = msh->bulkeles; ele; ele = ele->next)
      {
	render_object (bod, shp, SHAPE_ELEMENT, ele, flags, color);
      }
    }
    else render_object (bod, shp, SHAPE_MESH, shp->data, flags, color);
    break;
  case SHAPE_CONVEX:
    for (CONVEX *cvx = shp->data; cvx; cvx = cvx->next)
      render_object (bod, shp, SHAPE_CONVEX, cvx, flags, color);
    break;
  case SHAPE_SPHERE:
    for (SPHERE *sph = shp->data; sph; sph = sph->next)
      render_object (bod, shp, SHAPE_SPHERE, sph, flags, color);
    break;
  }
}

/* convert body id into an RGBA code */
static void idtorgba (int id, GLfloat color [4])
{
  unsigned char rgba [4];

  rgba [0] = ((unsigned char*)&id) [0];
  rgba [1] = ((unsigned char*)&id) [1];
  rgba [2] = ((unsigned char*)&id) [2];
  rgba [3] = ((unsigned char*)&id) [3];

  color [0] = (GLfloat) rgba [0] / 255.0f;
  color [1] = (GLfloat) rgba [1] / 255.0f;
  color [2] = (GLfloat) rgba [2] / 255.0f;
  color [3] = (GLfloat) rgba [3] / 255.0f;
}

/* convert RGBA to body id */
static int rgbatoid (unsigned char *rgba)
{
  int id;

  ((unsigned char*)&id) [0] = rgba [0];
  ((unsigned char*)&id) [1] = rgba [1];
  ((unsigned char*)&id) [2] = rgba [2];
  ((unsigned char*)&id) [3] = rgba [3];

  return id;
}

/* rendering in selection mode */
static void render_selection (int skip_flags)
{
  GLfloat color [4];

  if (selection)
  {
    for (MAP *item = MAP_First (selection); item; item = MAP_Next (item))
    {
      BODY *bod = item->data;

      if (bod->flags & skip_flags) continue;

      glPushName (bod->id);
      idtorgba (bod->id, color);
      for (SHAPE *shp = bod->shape; shp; shp = shp->next)
      {
	switch (shp->kind)
	{
	  case SHAPE_MESH:
	    render_mesh (NULL, shp, shp->data, SELECTION, color);
	    break;
	  case SHAPE_CONVEX:
            for (CONVEX *cvx = shp->data; cvx; cvx = cvx->next)
              render_convex (NULL, shp, cvx, SELECTION, color);
	    break;
	  case SHAPE_SPHERE:
            for (SPHERE *sph = shp->data; sph; sph = sph->next)
	      render_sphere (NULL, shp, sph, SELECTION, color);
	    break;
	}
      }
      glPopName ();
    }
  }
  else
  {
    for (BODY *bod = domain->bod; bod; bod = bod->next)
    {
      if (bod->flags & skip_flags) continue;

      glPushName (bod->id);
      idtorgba (bod->id, color);
      for (SHAPE *shp = bod->shape; shp; shp = shp->next)
      {
	switch (shp->kind)
	{
	  case SHAPE_MESH:
	    render_mesh (NULL, shp, shp->data, SELECTION, color);
	    break;
	  case SHAPE_CONVEX:
            for (CONVEX *cvx = shp->data; cvx; cvx = cvx->next)
              render_convex (NULL, shp, cvx, SELECTION, color);
	    break;
	  case SHAPE_SPHERE:
            for (SPHERE *sph = shp->data; sph; sph = sph->next)
	      render_sphere (NULL, shp, sph, SELECTION, color);
	    break;
	}
      }
      glPopName ();
    }
  }
}

/* selection using selection buffer */
static void select_bodies_3d (int x1, int y1, int x2, int y2)
{
  GLuint *sel, ssel, nsel, n, m, k;
  GLint viewport[4];
  int x, y, w, h; 

  w = ABS (x1 - x2);
  h = ABS (y1 - y2);
  x = (x1 + x2) / 2;
  y = (y1 + y2) / 2;
  w = MAX (w, 2);
  h = MAX (h, 2);
  ssel = domain->nbod * 4; /* assumes each body id stored in a separate hit, which is excessive and must be enough */
  ERRMEM (sel = malloc (sizeof (GLuint [ssel])));

  glSelectBuffer (ssel, sel);
  glRenderMode (GL_SELECT);
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glLoadIdentity ();
  glGetIntegerv (GL_VIEWPORT, viewport);
  gluPickMatrix (x, viewport [3] - y, w, h, viewport);
  GLV_SetProjectionMatrix (viewport [2], viewport [3]);
  glMatrixMode (GL_MODELVIEW);
  glInitNames ();
  render_selection (BODY_HIDDEN|BODY_OFF);
  glMatrixMode (GL_PROJECTION);
  glPopMatrix ();
  glMatrixMode (GL_MODELVIEW);
  glFlush ();
  nsel = glRenderMode (GL_RENDER);

  if (nsel > 0 && nsel < INT_MAX)
  {
    if (selection)
    {
      MAP *newsel = NULL;
      BODY *bod;

      for (n = m = 0; n < nsel; n ++, m += sel [m]+3)
      {
	for (k = 0; k < sel [m]; k ++)
	{
	  bod = MAP_Find (selection, (void*) sel [m+3+k], NULL);
	  WARNING_DEBUG (bod, "Invalid selection name: %d", sel [m+3+k]);
	  if (bod) MAP_Insert (&domain->mapmem, &newsel, (void*) bod->id, bod, NULL);
	}
      }

      MAP_Free (&domain->mapmem, &selection);
      selection = newsel;
    }
    else
    {
      BODY *bod;

      for (n = m = 0; n < nsel; n ++, m += sel [m]+3)
      {
	for (k = 0; k < sel [m]; k ++)
	{
	  bod = MAP_Find (domain->idb, (void*) sel [m+3+k], NULL);
	  WARNING_DEBUG (bod, "Invalid selection name: %d", sel [m+3+k]);
	  if (bod) MAP_Insert (&domain->mapmem, &selection, (void*) bod->id, bod, NULL);
	}
      }
    }
  }

  free (sel);
}

/* selection using coloring */
static void select_bodies_2d (int x1, int y1, int x2, int y2)
{
  unsigned char (*pix) [4];
  int x, y, w, h, n, m;
  GLint viewport[4];
  SET *item, *set;
  BODY *bod;

  glGetIntegerv (GL_VIEWPORT, viewport);

  x = MIN (x1, x2);
  y = MIN (viewport[3]-y1, viewport[3]-y2);
  w = MAX (x1, x2) - x;
  h = MAX (viewport[3]-y1, viewport[3]-y2) - y;
  w = MAX (w, 1);
  h = MAX (h, 1);
  m = w * h;

  glDisable (GL_LIGHTING);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  render_selection (BODY_HIDDEN|BODY_OFF);
  glEnable (GL_LIGHTING);

  ERRMEM (pix = malloc (m * sizeof (unsigned char [4])));
  glReadPixels (x, y, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pix);

  for (set = NULL, n = 0; n < m; n ++)
    SET_Insert (&domain->setmem, &set, (void*) rgbatoid (pix [n]), NULL);

  free (pix);

  if (set)
  {
    if (selection)
    {
      MAP *newsel = NULL;

      for (item = SET_First (set); item; item = SET_Next (item))
      {
	bod = MAP_Find (selection, item->data, NULL);
	if (bod) MAP_Insert (&domain->mapmem, &newsel, (void*) bod->id, bod, NULL);
      }

      MAP_Free (&domain->mapmem, &selection);
      selection = newsel;
    }
    else
    {
      for (item = SET_First (set); item; item = SET_Next (item))
      {
	bod = MAP_Find (domain->idb, item->data, NULL);
	if (bod) MAP_Insert (&domain->mapmem, &selection, (void*) bod->id, bod, NULL);
      }
    }
  }

  SET_Free (&domain->setmem, &set);
}

/* selection one body using coloring */
static BODY* select_body (int x, int y)
{
  unsigned char pix [4];
  GLint viewport[4];
  BODY *bod;

  glDisable (GL_LIGHTING);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  render_selection (BODY_OFF);
  glEnable (GL_LIGHTING);

  glGetIntegerv (GL_VIEWPORT, viewport);
  glReadPixels (x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pix);

  if (selection) bod = MAP_Find (selection, (void*) rgbatoid (pix), NULL);
  else bod = MAP_Find (domain->idb, (void*) rgbatoid (pix), NULL);

  return bod;
}

/* switch off selection */
static void selection_off ()
{

  mousemode = MOUSE_MODE_NONE;
  GLV_Rectangle_Off ();
  GLV_Release_Mouse ();

  /* update volumetric map (and its extrema) after selection */
  if (volumetric_map) volumetric_map_on (volumetric_map_kind);
  else GLV_Redraw_All ();
}

/* one solfec step in write mode */
static void solfec_step ()
{
  SOLVER_INTERFACE *si = MAP_Find (solver_interfaces, domain, NULL);

  if (si) SOLFEC_Run (solfec, si->kind, si->solver, domain->step);
  else SOLFEC_Run (solfec, EXPLICIT_SOLVER, NULL, domain->step); /* use explicit solver by default */

  /* update volumetric map (and its extrema) after updating state */
  if (volumetric_map) volumetric_map_on (volumetric_map_kind);
  else GLV_Redraw_All ();

  /* resize time viewport so to fit the text output */
  GLV_Resize_Viewport (time_window, time_width (), TIME_HEIGHT);
}

/* solfe run timer callback */
static void solfec_run (int value)
{
  solfec_step ();

  if (domain->flags & DOM_RUN_ANALYSIS)
    glutTimerFunc (1000 * SOLFEC_Time_Skip (solfec), solfec_run, 0);
}

/* menu callbacks */

static void menu_render (int value)
{
  switch (value)
  {
  case RENDER_ALL:

    MAP_Free (&domain->mapmem, &selection);
    mousemode = MOUSE_MODE_NONE;
    update_extents ();

    if (current_tool == TOOLS_CLIPPING_PLANE)
      update_clipping_plane (clipping_plane);

    /* update volumetric map (and its extrema) after selection */
    if (volumetric_map) volumetric_map_on (volumetric_map_kind);

    break;
  case RENDER_SELECTION_2D:
    deloff (tools_off);
    mousemode = MOUSE_SELECTION_2D_BEGIN;
    GLV_Hold_Mouse ();
    pushoff (selection_off);
    break;
  case RENDER_SELECTION_3D:
    deloff (tools_off);
    mousemode = MOUSE_SELECTION_3D_BEGIN;
    GLV_Hold_Mouse ();
    pushoff (selection_off);
    break;
  }
}

static void menu_domain (int value)
{
  switch (value)
  {
  case DOMAIN_NEXT:
    if (domain->next)
    {
      domain = domain->next;
      GLV_Redraw_All ();
    }
    break;
  case DOMAIN_PREV:
    if (domain->prev)
    {
      domain = domain->prev;
      GLV_Redraw_All ();
    }
    break;
  }

  /* update analysis menu */

  glutSetMenu (analysis_menu);

  if (domain->flags & DOM_RUN_ANALYSIS)
    glutChangeToMenuEntry (1, "stop /RETURN/", ANALYSIS_STOP);
  else glutChangeToMenuEntry (1, "run /RETURN/", ANALYSIS_RUN);

  if (solfec->mode == SOLFEC_READ && analysis_menu_items == 2)
  {
    glutAddMenuEntry ("seek to /UP/", ANALYSIS_SEEKTO);
    glutAddMenuEntry ("forward /LEFT/", ANALYSIS_FORWARD);
    glutAddMenuEntry ("backward /RIGHT/", ANALYSIS_BACKWARD);
    glutAddMenuEntry ("skip /DOWN/", ANALYSIS_SKIP);
    analysis_menu_items = 6;
  }
  else if (solfec->mode == SOLFEC_WRITE && analysis_menu_items == 6)
  {
    glutRemoveMenuItem (6);
    glutRemoveMenuItem (5);
    glutRemoveMenuItem (4);
    glutRemoveMenuItem (3);
    analysis_menu_items = 2;
  }
}

static void menu_tools (int value)
{
  deloff (selection_off);

  switch (value)
  {
  case TOOLS_TOGGLE_VISIBLE:
    if (!current_tool) pushoff (tools_off);
    current_tool = value;
    tool_flag = 0;
    break;
  case TOOLS_CLIPPING_PLANE:
    GLV_Read_Text ("Plane normal nx, ny, nz", read_clipping_normal);
    break;
  case TOOLS_SWITCH_OFF:
    if (!current_tool) pushoff (tools_off);
    current_tool = TOOLS_TOGGLE_VISIBLE;
    tool_flag = TOOL_FLAG_SWITCH_OFF;
    break;
  case TOOLS_SWITCH_ON_ALL:
    for (BODY *bod = domain->bod; bod; bod = bod->next)
      bod->flags &= ~(BODY_HIDDEN|BODY_OFF);
    GLV_Redraw_All ();
    break;
  case TOOLS_ROUGH_MESH_ON:
    if (!current_tool) pushoff (tools_off);
    current_tool = TOOLS_TOGGLE_VISIBLE;
    tool_flag = TOOL_FLAG_ROUGH_MESH_ON;
    break;
  }
}

static void menu_kinds (int value)
{
  volumetric_map_on (value);
  GLV_Redraw_All ();
}

static void menu_analysis (int value)
{
  switch (value)
  {
  case ANALYSIS_RUN:
    domain->flags |= DOM_RUN_ANALYSIS;
    solfec_run (0);
    glutSetMenu (analysis_menu);
    glutChangeToMenuEntry (1, "stop /RETURN/", ANALYSIS_STOP);
    break;
  case ANALYSIS_STOP:
    domain->flags &= ~DOM_RUN_ANALYSIS;
    glutSetMenu (analysis_menu);
    glutChangeToMenuEntry (1, "run /RETURN/", ANALYSIS_RUN);
    break;
  case ANALYSIS_STEP:
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    solfec_step ();
    break;
  case ANALYSIS_SEEKTO:
    {
      static char caption [128];
      double s, e;

      SOLFEC_Time_Limits (solfec, &s, &e);
      snprintf (caption, 128, "Seek to time from [%.2e, %.2e]", s, e);
      GLV_Read_Text (caption, seek_to_time);
    }
    break;
  case ANALYSIS_FORWARD:
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    SOLFEC_Forward (solfec, analysis_skip_steps);
    if (volumetric_map) volumetric_map_on (volumetric_map_kind);
    else GLV_Redraw_All ();
    /* resize time viewport so to fit the text output */
    GLV_Resize_Viewport (time_window, time_width (), TIME_HEIGHT);
    break;
  case ANALYSIS_BACKWARD:
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    SOLFEC_Backward (solfec, analysis_skip_steps);
    if (volumetric_map) volumetric_map_on (volumetric_map_kind);
    else GLV_Redraw_All ();
    /* resize time viewport so to fit the text output */
    GLV_Resize_Viewport (time_window, time_width (), TIME_HEIGHT);
    break;
  case ANALYSIS_SKIP:
    GLV_Read_Text ("FORWARD and BACKWARD skip", set_skip_steps);
    break;
  }
}

static void menu_results (int value)
{
  volumetric_map_on (value);
  GLV_Redraw_All ();
}

int RND_Menu (char ***names, int **codes)
{
  static int code [9];
  static char *name [9];
  int local [9];

  ASSERT (domain, ERR_RND_NO_DOMAIN);

  name [0] = "domain";
  code [0] = domain_menu = glutCreateMenu (menu_domain);
  glutAddMenuEntry ("previous /</", DOMAIN_PREV);
  glutAddMenuEntry ("next />/", DOMAIN_NEXT);

  name [1] = "render";
  code [1] = render_menu = glutCreateMenu (menu_render);
  glutAddMenuEntry ("all bodies /a/", RENDER_ALL);
  glutAddMenuEntry ("2D selection /2/", RENDER_SELECTION_2D);
  glutAddMenuEntry ("3D selection /3/", RENDER_SELECTION_3D);

  name [2] = "tools";
  code [2] = tools_menu = glutCreateMenu (menu_tools);
  glutAddMenuEntry ("toggle visible /v/", TOOLS_TOGGLE_VISIBLE);
  glutAddMenuEntry ("clipping plane /c/", TOOLS_CLIPPING_PLANE);
  glutAddMenuEntry ("switch off /f/", TOOLS_SWITCH_OFF);
  glutAddMenuEntry ("switch on all /o/", TOOLS_SWITCH_ON_ALL);
  glutAddMenuEntry ("rough mesh on /r/", TOOLS_ROUGH_MESH_ON);

  name [3] = "kinds of";
  code [3] = kindsof_menu = glutCreateMenu (menu_kinds);
  glutAddMenuEntry ("constraints", KINDSOF_CONSTRAINTS);
  glutAddMenuEntry ("forces", KINDSOF_FORCES);
  glutAddMenuEntry ("bodies", KINDSOF_BODIES);
  glutAddMenuEntry ("surfaces", KINDSOF_SURFACES);
  glutAddMenuEntry ("volumes", KINDSOF_VOLUMES);

  name [4] = "analysis";
  code [4] = analysis_menu = glutCreateMenu (menu_analysis);
  glutAddMenuEntry ("run /RETURN/", ANALYSIS_RUN);
  glutAddMenuEntry ("step /SPACE/", ANALYSIS_STEP);
  analysis_menu_items = 2;
  if (solfec->mode == SOLFEC_READ)
  {
    glutAddMenuEntry ("seek to /UP/", ANALYSIS_SEEKTO);
    glutAddMenuEntry ("forward /LEFT/", ANALYSIS_FORWARD);
    glutAddMenuEntry ("backward /RIGHT/", ANALYSIS_BACKWARD);
    glutAddMenuEntry ("skip /DOWN/", ANALYSIS_SKIP);
    analysis_menu_items = 6;
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

  name [5] = "results";
  code [5] = results_menu = glutCreateMenu (menu_results);
  glutAddSubMenu ("displacements", local [0]);
  glutAddSubMenu ("velocities", local [1]);
  glutAddSubMenu ("stresses", local [2]);
  glutAddSubMenu ("reactions", local [3]);

  *codes = code;
  *names = name;
  return 6;
}

/* initialize rendering */
void RND_Init ()
{
  double extents [6];
  int w, h;

  get_scene_extents (extents);

  GLV_Reset_Extents (extents);

  GLV_Sizes (&w, &h);
  time_window = GLV_Open_Viewport (0, -(h - TIME_HEIGHT), time_width (), TIME_HEIGHT, 0, render_time);
}

/* idle actions */
int  RND_Idle ()
{
  return 0;
}

/* finalize rendering */
void RND_Quit ()
{
  MAP *item;

  for (item = MAP_First (solver_interfaces); item; item = MAP_Next (item)) free (item->data); /* free solver interfaces */
  MAP_Free (NULL, &solver_interfaces); /* and the map itself */

  lngfinalize (); /* finalize Python context */
}

/* render scene */
void RND_Render ()
{
  int flags = OUTLINE;

  if (volumetric_map_constraint_based ())
  {
    GLfloat color [4] = {1.0, 0.0, 0.0, 1.0};

    for (CON *con = domain->con; con; con = con->next)
    {
      if (skip_constraint (con)) continue;

      get_constraint_based_color (con, color);

      switch (volumetric_map_kind)
      {
      case KINDSOF_CONSTRAINTS:
	switch (con->kind)
	{
	  case CONTACT: render_contact (con, color); break;
	  case FIXPNT: render_fixpnt (con, color); break;
	  case FIXDIR: render_fixdir (con, color); break;
	  case VELODIR: render_velodir (con, color); break;
	  case RIGLNK: render_riglnk (con, color); break;
	}
	break;
      case RESULTS_RT: render_rt (con, color); break;
      case RESULTS_RN: render_rn (con, color); break;
      case RESULTS_R:  render_r  (con, color); break;
      }
    }

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }

  if (selection)
  {
    for (MAP *item = MAP_First (selection); item; item = MAP_Next (item))
    {
      BODY *bod = item->data;

      GLfloat color [4] = {0.7, 0.7, 0.7, 0.2};

      if (render_body_begin (bod, color)) continue;

      for (SHAPE *shp = bod->shape; shp; shp = shp->next)
	render_shape (bod, shp, flags, color);

      render_body_end (bod, color);

      if (bod->flags & BODY_SHOW_ROUGH_MESH)
      {
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	render_object (bod, NULL, SHAPE_MESH, bod->msh, flags, color);
	glDisable (GL_BLEND);
      }
    }
  }
  else
  {
    for (BODY *bod = domain->bod; bod; bod = bod->next)
    {
      GLfloat color [4] = {0.7, 0.7, 0.7, 0.2};

      if (render_body_begin (bod, color)) continue;

      for (SHAPE *shp = bod->shape; shp; shp = shp->next)
	render_shape (bod, shp, flags, color);

      render_body_end (bod, color);

      if (bod->flags & BODY_SHOW_ROUGH_MESH)
      {
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	render_object (bod, NULL, SHAPE_MESH, bod->msh, flags, color);
	glDisable (GL_BLEND);
      }
    }
  }

  if (volumetric_map_constraint_based ())
  {
    glDisable (GL_BLEND);
  }
  else
  {
    for (CON *con = domain->con; con; con = con->next)
    {
      if (skip_constraint (con)) continue;
      if (con->kind == RIGLNK) 
      {
        GLfloat color [4] = {0.7, 0.7, 0.7, 0.2};

	render_riglnk (con, color); /* render rigid links */
      }
    }
  }
}

/* react to key pressed */
void RND_Key (int key, int x, int y)
{
  switch (key)
  {
  case 27:
    popoff (); /* switch off stuff */
    break;
  case '\r':
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    else menu_analysis (ANALYSIS_RUN);
    break;
  case ' ':
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    solfec_step ();
    break;
  case '<':
    menu_domain (DOMAIN_PREV);
    break;
  case '>':
    menu_domain (DOMAIN_NEXT);
    break;
  case 'a':
    menu_render (RENDER_ALL);
    break;
  case '2':
    menu_render (RENDER_SELECTION_2D);
    break;
  case '3':
    menu_render (RENDER_SELECTION_3D);
    break;
  case 'v':
    menu_tools (TOOLS_TOGGLE_VISIBLE);
    break;
  case 'c':
    menu_tools (TOOLS_CLIPPING_PLANE);
    break;
  case 'f':
    menu_tools (TOOLS_SWITCH_OFF);
    break;
  case 'o':
    menu_tools (TOOLS_SWITCH_ON_ALL);
    break;
  case 'r':
    menu_tools (TOOLS_ROUGH_MESH_ON);
    break;
  case '=':
    if (volumetric_map)
    {
      if (volumetric_map_kind < RESULTS_R)
	volumetric_map_on (volumetric_map_kind + 1);
    }
    else volumetric_map_on (KINDSOF_CONSTRAINTS);
    GLV_Redraw_All ();
    break;
  case '-':
    if (volumetric_map)
    {
      if (volumetric_map_kind > KINDSOF_CONSTRAINTS)
	volumetric_map_on (volumetric_map_kind - 1);
    }
    else volumetric_map_on (RESULTS_R);
    GLV_Redraw_All ();
    break;
  }
}

/* react to special key pressed */
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

/* react to mouse click */
void RND_Mouse (int button, int state, int x, int y)
{
  switch (state)
  {
    case GLUT_UP:
    break;

    case GLUT_DOWN:
    break;
  }

  switch (button)
  {
    case GLUT_LEFT_BUTTON:
    switch (state)
    {
    case GLUT_UP:
      if (mousemode == MOUSE_SELECTION_2D_END||
          mousemode == MOUSE_SELECTION_3D_END)
      {
	if (mousemode == MOUSE_SELECTION_2D_END)
	  select_bodies_2d (mouse_x, mouse_y, x, y);
	else select_bodies_3d (mouse_x, mouse_y, x, y);

        update_extents ();

	deloff (selection_off);
      }
      break;
    case GLUT_DOWN:
      if (mousemode == MOUSE_SELECTION_2D_BEGIN)
      {
	mouse_x = x;
	mouse_y = y;
	mousemode = MOUSE_SELECTION_2D_END;
      }
      else if (mousemode == MOUSE_SELECTION_3D_BEGIN)
      {
	mouse_x = x;
	mouse_y = y;
	mousemode = MOUSE_SELECTION_3D_END;
      }
      else if (current_tool == TOOLS_TOGGLE_VISIBLE)
      {
        if (toggle_visible)
	{
	  if ((glutGetModifiers () & GLUT_ACTIVE_CTRL)
	      || tool_flag == TOOL_FLAG_SWITCH_OFF)
	  {
	    toggle_visible->flags |= BODY_OFF; /* switch off */
	  }
	  else if (tool_flag == TOOL_FLAG_ROUGH_MESH_ON)
	  {
	    if (toggle_visible->flags & BODY_SHOW_ROUGH_MESH)
	      toggle_visible->flags &= ~BODY_SHOW_ROUGH_MESH;
	    else if (toggle_visible->msh) toggle_visible->flags |= BODY_SHOW_ROUGH_MESH;
	  }
	  else /* hide */
	  {
	    if (toggle_visible->flags & BODY_HIDDEN)
	      toggle_visible->flags &= ~BODY_HIDDEN; /* switch to visible */
	    else toggle_visible->flags |= BODY_HIDDEN; /* hide */

	    /* update volumetric map (and its extrema) after toggling state */
	    if (volumetric_map) volumetric_map_on (volumetric_map_kind);
	  }
	}
      }
      break;
    }
    break;

    case GLUT_MIDDLE_BUTTON:
    switch (state)
    {
    case GLUT_UP:
      break;
    case GLUT_DOWN:
      break;
    }
    break;

    case GLUT_RIGHT_BUTTON:
    switch (state)
    {
    case GLUT_UP:
      break;
    case GLUT_DOWN:
      break;
    }
    break;
  }
}

/* react to mouse motion */
void RND_Motion (int x, int y)
{
  if (mousemode == MOUSE_SELECTION_2D_END||
      mousemode == MOUSE_SELECTION_3D_END)
  {
    GLV_Rectangle_On (mouse_x, mouse_y, x, y);
  }
}

/* react on passive mouse motion */
void RND_Passive (int x, int y)
{
  switch (current_tool)
  {
  case TOOLS_TOGGLE_VISIBLE:
    toggle_visible = select_body (x, y);
    GLV_Redraw_All ();
    break;
  case TOOLS_CLIPPING_PLANE:
    move_clipping_plane (x, y);
    GLV_Redraw_All ();
    break;
  }
}

/* enable renedring before opening viewer */
void RND_Viewer_On ()
{
  vieweron = 1; /* [-v] option was specified by the user => viewer is enabled */
}

/* set domain to be rendered */
void RND_Domain (DOM *dom)
{
  dom->next = domain;
  if (domain) domain->prev = dom;
  domain = dom;
}

/* map solver to a domain */
void RND_Solver (DOM *dom, int kind, void *solver)
{
  SOLVER_INTERFACE *si;
  MAP *item;

  ERRMEM (si = malloc (sizeof (SOLVER_INTERFACE))); /* create solver interfaces */
  si->kind = kind;
  si->solver = solver;

  if ((item = MAP_Find (solver_interfaces, dom, NULL))) free (item->data); /* if previously mapped: free it */

  MAP_Insert (NULL, &solver_interfaces, dom, si, NULL); /* map solver interface */
}

/* check whether rendering is on */
int  RND_On ()
{
  return vieweron;
}
