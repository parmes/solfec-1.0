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
  #include <GL/glext.h>
#endif
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "solfec.h"
#include "set.h"
#include "alg.h"
#include "glv.h"
#include "sol.h"
#include "msh.h"
#include "cvx.h"
#include "sph.h"
#include "eli.h"
#include "shp.h"
#include "lng.h"
#include "rnd.h"
#include "fem.h"
#include "err.h"

typedef struct value_source VALUE_SOURCE; /* vertex value source */

struct value_source
{
  ELEPNT *epn;

  double *pnt;
};

typedef struct cut_data CUT_DATA;

struct cut_data
{
  enum {EULER, LAGRANGE} kind;

  double point [3],
	 normal [3];

  BODY *bod;

  TRI *tri;

  SGP *sgp;

  double *ref, /* referential counterpartes of triangle vertices */
	 *cur, /* pointed in 'tri' */
	 *val; /* vertex value */

  int m, /* triangles */
      n; /* vertices */

  short values_updated; /* update flag */

  CUT_DATA *next;
};

typedef struct body_data BODY_DATA; /* body rendering data */

struct body_data
{
  enum {SEETHROUGH  = 0x01,         /* transparency flag */
        HIDDEN      = 0x02,         /* hidden state */
	ROUGH_MESH  = 0x04,         /* rough mesh rendering */
        WIREFRAME   = 0x08} flags;  /* wireframe mode */

#if VBO
  GLuint triangles, /* VBO of triangle vertices, normals and colors */
	 lines; /* VBO of line vertices */
#else
  GLfloat *ver, /* vertex array data */
	  *lin;
#endif

  GLsizei triangles_count,
	  lines_count;

  double **vertex_sources,
	 **normal_sources,
	  *vertex_values, /* as many values as vertices (vertex_vlues_count) */
	 **vertex_value_sources, /* maps vertices to 'values' */
	 **line_sources;

  double *values; /* unique vertex scalar field values */

  VALUE_SOURCE *value_sources; /* value sources */

  int values_count,
      vertex_values_count; /* = triangles_count * 3 or lines_count * 2 in WIREFRAME mode */

  short values_updated; /* update flag */

  ELEPNT **color_sources; /* only for FEM */

  MAP *surfaces; /* map surface ids to sets of &vertex_values [i] pointers */

  MAP *volumes; /* map volume ids to sets of &vertex_values [i] pointers */

  SPHERE **spheres; /* spheres */

  GLfloat *sphere_colors; /* sphere surface colors */

  int spheres_count;

  ELLIP **ellips; /* ellipsoids */

  GLfloat *ellip_colors; /* ellipsoid surface colors */

  int ellips_count;

  BODY_DATA *rough; /* rough mesh rendering */
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
  RENDER_WIREFRAME_2D,
  RENDER_WIREFRAME_3D,
  RENDER_PREVIOUS_SELECTION,
  RENDER_BODIES,
  TOOLS_TRANSPARENT,
  TOOLS_TRANSPARENT_ALL,
  TOOLS_TRANSPARENT_NONE,
  TOOLS_EULER_CUT,
  TOOLS_LAGRANGE_CUT,
  TOOLS_CLEAR_CUTS,
  TOOLS_HIDE,
  TOOLS_SHOW_ALL,
  TOOLS_ROUGH_MESH,
  TOOLS_PREVIOUS_RESULT,
  TOOLS_NEXT_RESULT,
  TOOLS_SMALLER_SCALING,
  TOOLS_BIGGER_SCALING,
  TOOLS_OUTPATH,
  TOOLS_POINTS_COORDS,
  TOOLS_POINTS_DISTANCE,
  TOOLS_POINTS_ANGLE,
  TOOLS_TRACKBALL_CENTER,
  TOOLS_WIREFRAME_ALL,
  TOOLS_WIREFRAME_NONE,
  TOOLS_NEXT_MODE,
  TOOLS_PREVIOUS_MODE,
  ANALYSIS_RUN,
  ANALYSIS_STOP,
  ANALYSIS_STEP,
  ANALYSIS_SEEKTO,
  ANALYSIS_FORWARD,
  ANALYSIS_BACKWARD,
  ANALYSIS_SKIP,
  KINDS_OF_CONSTRAINTS,
  KINDS_OF_CONSTRAINT_RANKS,
  KINDS_OF_FORCES,
  KINDS_OF_BODIES,
  KINDS_OF_BODY_RANKS,
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
  RESULTS_R,
  RESULTS_UT,
  RESULTS_UN,
  RESULTS_U,
  RESULTS_GAP,
  RESULTS_MERIT,
  RESULTS_MODAL_ANALYSIS
};

enum mouse_mode
{
  MOUSE_NONE,
  MOUSE_SELECTION_BEGIN,
  MOUSE_SELECTION_END,
  MOUSE_PICK_BODY,
  MOUSE_PICK_POINT,
  MOUSE_CUTTING_PLANE
};

typedef enum mouse_mode MOUSE_MODE;

enum selection_mode
{
  SELECTION_NONE,
  SELECTION_2D,
  SELECTION_3D
};

typedef enum selection_mode SELECTION_MODE;

typedef struct legend_data LEGEND_DATA;

struct legend_data
{
  short preventity;

  short entity;

  double extents [2];

  int range; /* number of items in the legend */

  int window;

  SET *discrete;

  double constant; /* value_to_color mapping regularisation constant */
};

/* declarations */

static void menu_analysis (int);
static void menu_kinds (int);
static void menu_results (int);
static void update ();
static int  time_width ();

/* global data */

#define CHUNK 256 /* memory pool items chunk */

#define BIGCHUNK 1024 /* memory pool items big chunk */

static short enabled = 0; /* renderin on/off */

static DOM *domain = NULL; /* current domain */

#define solfec domain->solfec /* current solfec */

static MAP *solvers = NULL; /* maps domains to solvers */

static char* menu_name [MENU_LAST];  /* menu names */

static int menu_code [MENU_LAST]; /* menu codes */

static SELECTION *selection = NULL; /* current selection */

static MEM rndsetmem; /* sets memory */

static MEM rndmapmem;  /* maps memory */

static int skip_steps = 1; /* number of steps to skip when rewinding analysis */

static MOUSE_MODE mouse_mode = MOUSE_NONE;

static SELECTION_MODE selection_mode = SELECTION_NONE;

static GLfloat neutral_color [3] = {0.8, 0.8, 0.8};

static int mouse_start [2] = {0, 0};

static short tool_mode = 0; /* current tool */

#define PICKED_BODY_TIP_LEN 1024
static BODY *picked_body = NULL; /* currently picked body */
static char picked_body_tip [PICKED_BODY_TIP_LEN];

static double *picked_point = NULL, /* currently picked point */
	      *picked_point_hist [2] = {NULL, NULL}; /* history of picked points */

static double arrow_factor = 0.05; /* arrow drawing constant */

static int time_window = 0; /* time window handle */
#define TIME_HEIGHT 16 /* time window height */
#define TIME_FONT GLV_FONT_8_BY_13

static int coord_window = 0; /* coordinates window handle */
#define COORD_WIDTH 64
#define COORD_HEIGHT 36

static LEGEND_DATA legend; /* legend data */
#define LEGEND_ROWS 8 /* number of rows in the legend */
#define LEGEND_WIDTH_DISC 50 /* legend width for discrete data */
#define LEGEND_WIDTH_CONT 100 /* legend width for continuous data */
#define LEGEND_FONT GLV_FONT_8_BY_13

static char *tip; /* tip upper */

/* body transparency test */
#define SEETHROUGH(body) (((BODY_DATA*)((BODY*)(body))->rendering)->flags & SEETHROUGH)

/* body rough mesh flag test */
#define ROUGH_MESH(body) (((BODY_DATA*)((BODY*)(body))->rendering)->flags & ROUGH_MESH)

static short show_outpath = 0; /* output path printing flag */

static double global_extents [6] = {0, 0, 0, 1, 1, 1}; /* global scene extents */

static short cut_sketch = 0; /* sketch cut plane */
static double cut_point [3] = {0, 0, 0}; /* current cut point */
static double cut_normal [3] = {0, 0, 1}; /* current cut normal */
static short cut_kind = EULER; /* cut kind */
static CUT_DATA *cuts = NULL; /* cuts through bodies */
static SET *euler_cuts = NULL; /* points and normal of current Euler cuts */

static short render_bodies = 1; /* body rendering flag */

static short wireframeon = 0; /* wireframe selection flag */

static short modal_analysis_menu = 0; /* modal analysis menu item flag */
static double eigenshape_factor = 0.15; /* eigenshapes scaling factor */
static int current_eigenmode = -1; /* current eigenmode */
static char modetip [512]; /* mode number and eigenvalue tip */


/* menu modal analysis callback */
static void menu_modal_analysis (int mode)
{
  current_eigenmode = mode;

  if (SET_Prev (selection->set) == NULL &&
      SET_Next (selection->set) == NULL) /* just one body */
  {
    BODY *bod = selection->set->data;
    double *e = bod->extents, d [3], scale;

    SUB (e+3, e, d);
    scale = eigenshape_factor * (d[0] + d[1] + d [2]) / 3.0;

    FEM_Load_Mode (bod, mode, scale);

    snprintf (modetip, 512, "mode %d, %g", mode+1, bod->eval [mode]); 
    tip = modetip;

    update ();
  }
}

/* modal analaysis results menu item turn on; this will depend
 * on selection => it can only be turned on when just one body
 * is selected and this body has the modal analysis results attached to it */
static void modal_analysis_results ()
{
  /* computing set is much more time consuming than
   * checking for previous and next item pointers */

  if (SET_Prev (selection->set) == NULL &&
      SET_Next (selection->set) == NULL) /* just one body */
  {
    BODY *bod = selection->set->data;

    if (bod->eval && bod->evec && modal_analysis_menu == 0)
    {
      int local = glutCreateMenu (menu_modal_analysis);

      for (int i = 0; i < bod->evec->n; i ++)
      {
        char mode [256];
	snprintf (mode, 256, "mode %d   (%g)", i+1, bod->eval [i]);
        glutAddMenuEntry (mode, i);
      }

      glutSetMenu (menu_code [MENU_RESULTS]);
      glutAddSubMenu ("modal analysis", local);

      glutSetMenu (menu_code [MENU_TOOLS]);
      glutChangeToMenuEntry (12, "bigger eigenshapes /]/", TOOLS_BIGGER_SCALING);
      glutChangeToMenuEntry (13, "smaller eigenshapes /[/", TOOLS_SMALLER_SCALING);
      glutAddMenuEntry ("next eigenshape /}/", TOOLS_NEXT_MODE);
      glutAddMenuEntry ("previous eigenshape /{/", TOOLS_PREVIOUS_MODE);

      modal_analysis_menu = 1;
    }
  }
  else
  {
    if (modal_analysis_menu)
    {
      glutSetMenu (menu_code [MENU_RESULTS]);
      glutRemoveMenuItem (5);

      glutSetMenu (menu_code [MENU_TOOLS]);
      glutChangeToMenuEntry (12, "bigger arrows /]/", TOOLS_BIGGER_SCALING);
      glutChangeToMenuEntry (13, "smaller arrows /[/", TOOLS_SMALLER_SCALING);
      glutRemoveMenuItem (22);
      glutRemoveMenuItem (21);

      if (tip == modetip) tip = NULL;

      modal_analysis_menu = 0;
    }
  }
}

/* initialize selection */
static void selection_init ()
{
  SELECTION *prev;
  BODY *bod;

  while (selection)
  {
    SET_Free (&rndsetmem, &selection->set);
    prev = selection->prev;
    free (selection);
    selection = prev;
  }
  
  ERRMEM (selection = malloc (sizeof (SELECTION)));

  for (selection->set = NULL, bod = domain->bod; bod; bod = bod->next) SET_Insert (&rndsetmem, &selection->set, bod, NULL);

  selection->prev = NULL;

  modal_analysis_results ();
}

/* push new selection on stack */
static void selection_push (SET *set)
{
  SELECTION *s;

  ERRMEM (s = malloc (sizeof (SELECTION)));
  s->prev = selection;
  s->set = set;

  selection = s;

  modal_analysis_results ();
}

/* pop most recent selection */
static void selection_pop ()
{
  if (selection->prev)
  {
    SELECTION *top = selection;
    SET_Free (&rndsetmem, &selection->set);
    selection = selection->prev;
    free (top);
  }

  modal_analysis_results ();
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
#if __MINGW32__ || OSTYPE_WIN32 || OSTYPE_LINUX
  ((unsigned char*)&id) [3] = 0; /* FIXME: alpha always one */
#else
  ((unsigned char*)&id) [3] = rgba [3];
#endif

  return id;
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

/* translate scalar value into color */
inline static void value_to_color (double value, GLfloat *color)
{
  HSL_2_RGB (0.69 * (1.0-(value-legend.extents[0])/(legend.extents[1]-legend.extents[0]+legend.constant)), 1.0, 0.45, color);
}

/* is legend constraint based */
static short legend_constraint_based (void)
{
  if ((legend.entity >= KINDS_OF_CONSTRAINTS &&
       legend.entity <= KINDS_OF_FORCES) ||
      (legend.entity >= RESULTS_RT &&
       legend.entity <= RESULTS_MERIT)) return 1;

  return 0;
}

/* obtain (referential) point value */
static double point_value (BODY *bod, SHAPE *shp, void *gobj, double *X)
{
  double values [7];
  VALUE_KIND kind;
  short index;

  switch (legend.entity)
  {
  case RESULTS_DX:
  case RESULTS_DY:
  case RESULTS_DZ:
    kind = VALUE_DISPLACEMENT;
    index = legend.entity - RESULTS_DX;
    break;
  case RESULTS_VX:
  case RESULTS_VY:
  case RESULTS_VZ:
    kind = VALUE_VELOCITY;
    index = legend.entity - RESULTS_VX;
    break;
  case RESULTS_SX:
  case RESULTS_SY:
  case RESULTS_SZ:
  case RESULTS_SXY:
  case RESULTS_SXZ:
  case RESULTS_SYZ:
    kind = VALUE_STRESS;
    index = legend.entity - RESULTS_SX;
    break;
  case RESULTS_MISES:
    kind = VALUE_MISES;
    index = 0;
    break;
  default:
    return 0.0;
  }

  if (bod->kind == FEM)
  {
    ELEMENT *ele = NULL;

    if (shp->kind == SHAPE_CONVEX) /* convices with background mesh */
    {
      MESH *msh = bod->msh;
      CONVEX *cvx = gobj;
      double dist, d;
      int i;

      ASSERT_DEBUG (msh, "Background mesh is missing!");

      for (i = 0; i < cvx->nele; i ++)
      {
	if (ELEMENT_Contains_Point (msh, cvx->ele [i], X, 1)) /* look through overlapped elements */
	{
	  ele = cvx->ele [i];
	  break;
	}
      }

      if (!ele) /* none was contained X => find the closest one */
      {
        for (dist = DBL_MAX, i = 0; i < cvx->nele; i ++)
	{
	  d = ELEMENT_Point_Distance (msh, cvx->ele [i], X, 1);
	  if (d < dist) dist = d, ele = cvx->ele [i];
	}
      }
    }
    else ele = gobj;

    FEM_Point_Values (bod, ele, X, kind, values);
  }
  else BODY_Point_Values (bod, X, kind, values);

  return values [index];
}

/* obtain scalar sphere point value */
static double sphere_value (BODY *bod, SPHERE *sph)
{
  switch (legend.entity)
  {
  case KINDS_OF_BODIES: return bod->kind;
  case KINDS_OF_SURFACES:  return sph->surface;
  case KINDS_OF_VOLUMES: return sph->volume;
  case KINDS_OF_BODY_RANKS: return bod->rank;
  default: return point_value (bod, NULL, NULL, sph->ref_center);
  }

  return 0.0;
}

/* obtain scalar ellipsoid point value */
static double ellip_value (BODY *bod, ELLIP *eli)
{
  switch (legend.entity)
  {
  case KINDS_OF_BODIES: return bod->kind;
  case KINDS_OF_SURFACES:  return eli->surface;
  case KINDS_OF_VOLUMES: return eli->volume;
  case KINDS_OF_BODY_RANKS: return bod->rank;
  default: return point_value (bod, NULL, NULL, eli->ref_center);
  }

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

/* register an identifier in a map of sets */
static void register_identifier (MAP **map, int identifier, double *val)
{
  MAP *item;

  if (!(item = MAP_Find_Node (*map, (void*) (long) identifier, NULL)))
  {
    ERRMEM (item = MAP_Insert (&rndmapmem, map, (void*) (long) identifier, NULL, NULL));
  }

  SET_Insert (&rndsetmem, (SET**) &item->data, val, NULL);
}

/* create body rendering data */
static BODY_DATA* create_body_data (BODY *bod, int wireframe)
{
  double **vsr, **nsr, **lsr, **vvs, **end, *pla, *val;
  GLfloat *ver, *v, *nor, *n, *col, *c, *lin, *l;
  MAP *fmap, *emap, *vmap, *jtem;
  MEM pairmem, mapmem, setmem;
  VALUE_SOURCE *source;
  POINTER_PAIR *pair;
  SET *lset, *item;
  BODY_DATA *data;
  int i, j, *f;
  SPHERE *sph;
  CONVEX *cvx;
  ELLIP *eli;
  SHAPE *shp;
  MESH *msh;
  FACE *fac;

  ERRMEM (data = MEM_CALLOC (sizeof (BODY_DATA)));
  if (wireframe) data->flags |= WIREFRAME;

  MEM_Init (&pairmem, sizeof (POINTER_PAIR), CHUNK);
  MEM_Init (&mapmem, sizeof (MAP), CHUNK);
  MEM_Init (&setmem, sizeof (SET), CHUNK);
  fmap = NULL;
  emap = NULL;
  vmap = NULL;
  lset = NULL;

  for (shp = bod->shape; shp; shp = shp->next)
  {
    switch (shp->kind)
    {
    case SHAPE_MESH:
      msh = shp->data;
      for (fac = msh->faces; fac; fac = fac->n)
      {
	data->triangles_count += (fac->type - 2);

	for (i = 0; i < fac->type; i ++) MAP_Insert (&mapmem, &vmap, &msh->cur_nodes [fac->nodes [i]][0], NULL, NULL);

	for (i = 0; i < fac->type - 1; i ++)
	  register_line (&pairmem, &setmem, &lset, &msh->cur_nodes [fac->nodes [i]][0], &msh->cur_nodes[fac->nodes [i+1]][0]);
	register_line (&pairmem, &setmem, &lset, &msh->cur_nodes [fac->nodes [i]][0], &msh->cur_nodes [fac->nodes [0]][0]);
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
      {
        sph = shp->data;
	data->spheres_count ++;

	ERRMEM (data->spheres = realloc (data->spheres, data->spheres_count * sizeof (SPHERE*)));
	j = (data->spheres_count - 1);
	data->spheres [j] = sph;
      }
      break;
    case SHAPE_ELLIP:
      {
        eli = shp->data;
	data->ellips_count ++;

	ERRMEM (data->ellips = realloc (data->ellips, data->ellips_count * sizeof (ELLIP*)));
	j = (data->ellips_count - 1);
	data->ellips [j] = eli;
      }
      break;
    }
  }

  if (wireframe)
  {
    data->lines_count = SET_Size (lset);
    ERRMEM (lin = malloc (data->lines_count * sizeof (GLfloat) * 12));
    col = lin + data->lines_count * 6;
    ERRMEM (data->vertex_sources = malloc (data->lines_count * (sizeof (double) * 2) + data->lines_count * sizeof (double*) * 4));
    data->vertex_values = (double*) data->vertex_sources; /* here used for line vertex values */
    data->vertex_value_sources = (double**) (data->vertex_values + data->lines_count * 2);
    data->line_sources = data->vertex_value_sources + data->lines_count * 2;
    data->vertex_values_count = data->lines_count * 2;
    data->surfaces = NULL;
    data->volumes = NULL;
    data->triangles_count = 0;

    for (item = SET_First (lset), lsr = data->line_sources; item; item = SET_Next (item), lsr += 2)
    {
      pair = item->data;
      lsr [0] = pair->one;
      lsr [1] = pair->two;
    }

    data->values_count = MAP_Size (vmap); /* number of unique vertices */
    ERRMEM (data->values = MEM_CALLOC (data->values_count * sizeof (double)));
    ERRMEM (data->value_sources = malloc (data->values_count * sizeof (VALUE_SOURCE)));

    for (source = data->value_sources, jtem = MAP_First (vmap); jtem; source ++, jtem = MAP_Next (jtem))
    {
      if (bod->msh)
      {
	cvx = jtem->data; /* must have been mapped to a convex */
	source->epn = &cvx->epn [((double*)jtem->key - cvx->cur) / 3]; /* extract ELEPNT */
      }

      source->pnt = jtem->key; /* needed for both scalar field rendering and point picking */

      jtem->data = &data->values [source - data->value_sources]; /* map to source */
    }

    for (shp = bod->shape; shp; shp = shp->next)
    {
      switch (shp->kind)
      {
      case SHAPE_MESH:
	msh = shp->data;
	for (fac = msh->faces; fac; fac = fac->n)
	{
	  for (i = 1; i < fac->type - 1; i ++)
	  {
	    val = &msh->cur_nodes [fac->nodes [0]][0];
	    MAP_Insert (&mapmem, &fmap, val, (void*) (long) fac->surface, NULL); /* map vertices to surfaces */
	    MAP_Insert (&mapmem, &emap, val, (void*) (long) fac->ele->volume, NULL); /* and volumes */
	    val = &msh->cur_nodes [fac->nodes [i]][0];
	    MAP_Insert (&mapmem, &fmap, val, (void*) (long) fac->surface, NULL);
	    MAP_Insert (&mapmem, &emap, val, (void*) (long) fac->ele->volume, NULL);
	    val = &msh->cur_nodes [fac->nodes [i+1]][0];
	    MAP_Insert (&mapmem, &fmap, val, (void*) (long) fac->surface, NULL);
	    MAP_Insert (&mapmem, &emap, val, (void*) (long) fac->ele->volume, NULL);
	  }
	}
	break;
      case SHAPE_CONVEX:
	for (cvx = shp->data; cvx; cvx = cvx->next)
	{
	  for (f = cvx->fac, j = 0, pla = cvx->pla; j < cvx->nfac; f += f[0]+1, j ++, pla += 4)
	  {
	    for (i = 2; i <= f[0]-1; i ++)
	    {
	      val = &cvx->cur [f[1]];
	      MAP_Insert (&mapmem, &fmap, val, (void*) (long) cvx->surface [j], NULL);
	      MAP_Insert (&mapmem, &emap, val, (void*) (long) cvx->volume, NULL);
	      val = &cvx->cur [f[i]];
	      MAP_Insert (&mapmem, &fmap, val, (void*) (long) cvx->surface [j], NULL);
	      MAP_Insert (&mapmem, &emap, val, (void*) (long) cvx->volume, NULL);
	      val = &cvx->cur [f[i+1]];
	      MAP_Insert (&mapmem, &fmap, val, (void*) (long) cvx->surface [j], NULL);
	      MAP_Insert (&mapmem, &emap, val, (void*) (long) cvx->volume, NULL);
	    }
	  }
	}
	break;
      case SHAPE_SPHERE: break;
      case SHAPE_ELLIP: break;
      }
    }

    for (lsr = data->line_sources,
	 val = data->vertex_values,
	 vvs = data->vertex_value_sources,
	 end = lsr + data->lines_count * 2,
	 l = lin, c = col; lsr < end;
	 lsr ++, vvs ++, val ++, l += 3, c += 3)
    {
      ASSERT_DEBUG_EXT (*vvs = MAP_Find (vmap, *lsr, NULL), "Inconsistent vertex mapping");
      register_identifier (&data->surfaces, (int) (long) MAP_Find (fmap, *lsr, NULL), val); /* map surfaces to vertex values */
      register_identifier (&data->volumes, (int) (long) MAP_Find (emap, *lsr, NULL), val); /* map volumes to vertex values */
      COPY (*lsr, l);
      COPY (neutral_color, c);
    }

    if (data->spheres_count)
    {
      ERRMEM (data->sphere_colors = malloc (data->spheres_count * sizeof (GLfloat) * 3));

      for (c = data->sphere_colors, v = c + data->spheres_count * 3; c < v; c += 3)
      {
	COPY (neutral_color, c);
      }
    }

    if (data->ellips_count)
    {
      ERRMEM (data->ellip_colors = malloc (data->ellips_count * sizeof (GLfloat) * 3));

      for (c = data->ellip_colors, v = c + data->ellips_count * 3; c < v; c += 3)
      {
	COPY (neutral_color, c);
      }
    }

#if VBO
    glGenBuffersARB (1, &data->lines);
    glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->lines);
    glBufferDataARB (GL_ARRAY_BUFFER_ARB, data->lines_count * sizeof (GLfloat) * 12, lin, GL_DYNAMIC_DRAW_ARB);

    glGenBuffersARB (1, &data->triangles);
    glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);
    glBufferDataARB (GL_ARRAY_BUFFER_ARB, 0, NULL, GL_DYNAMIC_DRAW_ARB);

    free (lin);
#else
    data->lin = lin;
#endif
  }
  else
  {
    data->lines_count = SET_Size (lset);
    ERRMEM (lin = malloc (data->lines_count * sizeof (GLfloat) * 6));
    ERRMEM (ver = malloc (data->triangles_count * sizeof (GLfloat) * 27));
    nor = ver + data->triangles_count * 9;
    col = nor + data->triangles_count * 9;
    ERRMEM (data->vertex_sources = malloc (data->triangles_count * (sizeof (double*) * 9 + sizeof (double) * 3) + data->lines_count * sizeof (double*) * 2));
    data->normal_sources = data->vertex_sources + data->triangles_count * 3;
    data->vertex_values = (double*) (data->normal_sources + data->triangles_count * 3);
    data->vertex_value_sources = (double**) (data->vertex_values + data->triangles_count * 3);
    data->line_sources = data->vertex_value_sources + data->triangles_count * 3;
    data->vertex_values_count = data->triangles_count * 3;
    data->surfaces = NULL;
    data->volumes = NULL;

    for (item = SET_First (lset), lsr = data->line_sources; item; item = SET_Next (item), lsr += 2)
    {
      pair = item->data;
      lsr [0] = pair->one;
      lsr [1] = pair->two;
    }

    vsr = data->vertex_sources;
    nsr = data->normal_sources;
    val = data->vertex_values;

    for (shp = bod->shape; shp; shp = shp->next)
    {
      switch (shp->kind)
      {
      case SHAPE_MESH:
	msh = shp->data;
	for (fac = msh->faces; fac; fac = fac->n)
	{
	  for (i = 1; i < fac->type - 1; i ++, vsr += 3, nsr += 3, val += 3)
	  {
	    vsr [0] = &msh->cur_nodes [fac->nodes [0]][0];
	    vsr [1] = &msh->cur_nodes [fac->nodes [i]][0];
	    vsr [2] = &msh->cur_nodes [fac->nodes [i+1]][0];
	    nsr [0] = nsr [1] = nsr [2] = fac->normal;
	    register_identifier (&data->surfaces, fac->surface, &val [0]);
	    register_identifier (&data->surfaces, fac->surface, &val [1]);
	    register_identifier (&data->surfaces, fac->surface, &val [2]);
	    register_identifier (&data->volumes, fac->ele->volume, &val [0]);
	    register_identifier (&data->volumes, fac->ele->volume, &val [1]);
	    register_identifier (&data->volumes, fac->ele->volume, &val [2]);
	  }
	}
	break;
      case SHAPE_CONVEX:
	for (cvx = shp->data; cvx; cvx = cvx->next)
	{
	  for (f = cvx->fac, j = 0, pla = cvx->pla; j < cvx->nfac; f += f[0]+1, j ++, pla += 4)
	  {
	    for (i = 2; i <= f[0]-1; i ++, vsr += 3, nsr += 3, val += 3)
	    {
	      vsr [0] = &cvx->cur [f[1]];
	      vsr [1] = &cvx->cur [f[i]];
	      vsr [2] = &cvx->cur [f[i+1]];
	      nsr [0] = nsr [1] = nsr [2] = pla;
	      register_identifier (&data->surfaces, cvx->surface [j], &val [0]);
	      register_identifier (&data->surfaces, cvx->surface [j], &val [1]);
	      register_identifier (&data->surfaces, cvx->surface [j], &val [2]);
	      register_identifier (&data->volumes, cvx->volume, &val [0]);
	      register_identifier (&data->volumes, cvx->volume, &val [1]);
	      register_identifier (&data->volumes, cvx->volume, &val [2]);
	    }
	  }
	}
	break;
      case SHAPE_SPHERE: break;
      case SHAPE_ELLIP: break;
      }
    }

    data->values_count = MAP_Size (vmap); /* number of unique vertices */
    ERRMEM (data->values = MEM_CALLOC (data->values_count * sizeof (double)));
    ERRMEM (data->value_sources = malloc (data->values_count * sizeof (VALUE_SOURCE)));

    for (source = data->value_sources, jtem = MAP_First (vmap); jtem; source ++, jtem = MAP_Next (jtem))
    {
      if (bod->msh)
      {
	cvx = jtem->data; /* must have been mapped to a convex */
	source->epn = &cvx->epn [((double*)jtem->key - cvx->cur) / 3]; /* extract ELEPNT */
      }

      source->pnt = jtem->key; /* needed for both scalar field rendering and point picking */

      jtem->data = &data->values [source - data->value_sources]; /* map to source */
    }

    for (vsr = data->vertex_sources,
	 nsr = data->normal_sources,
	 vvs = data->vertex_value_sources,
	 end = vsr + data->vertex_values_count,
	 v = ver, n = nor, c = col; vsr < end;
	 vsr ++, nsr ++, vvs ++, v += 3, n += 3, c += 3)
    {
      COPY (*vsr, v);
      COPY (*nsr, n);
      ASSERT_DEBUG_EXT (*vvs = MAP_Find (vmap, *vsr, NULL), "Inconsistent vertex mapping");
      COPY (neutral_color, c);
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

      for (c = data->sphere_colors, v = c + data->spheres_count * 3; c < v; c += 3)
      {
	COPY (neutral_color, c);
      }
    }

    if (data->ellips_count)
    {
      ERRMEM (data->ellip_colors = malloc (data->ellips_count * sizeof (GLfloat) * 3));

      for (c = data->ellip_colors, v = c + data->ellips_count * 3; c < v; c += 3)
      {
	COPY (neutral_color, c);
      }
    }

#if VBO
    glGenBuffersARB (1, &data->lines);
    glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->lines);
    glBufferDataARB (GL_ARRAY_BUFFER_ARB, data->lines_count * sizeof (GLfloat) * 6, lin, GL_DYNAMIC_DRAW_ARB);

    glGenBuffersARB (1, &data->triangles);
    glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);
    glBufferDataARB (GL_ARRAY_BUFFER_ARB, data->triangles_count * sizeof (GLfloat) * 27, ver, GL_DYNAMIC_DRAW_ARB);

    free (lin);
    free (ver);
#else
    data->lin = lin;
    data->ver = ver;
#endif
  }

  MEM_Release (&setmem);
  MEM_Release (&mapmem);
  MEM_Release (&pairmem);

  return data;
}

/* update body set constraint or force legend values */
static void update_body_constraint_or_force_values (BODY *bod)
{
  BODY_DATA *data = bod->rendering;
  double value;
  FORCE *force;
  SET *item;
  CON *con;

  if (data->flags & HIDDEN) return;

  if (legend.entity != KINDS_OF_FORCES)
  {
    for (item = SET_First (bod->con); item; item = SET_Next (item))
    {
      con = item->data;

      if (con->state & CON_DONE) continue;

      switch (legend.entity)
      {
      case KINDS_OF_CONSTRAINTS: 
        SET_Insert (&rndsetmem, &legend.discrete, (void*) (long) con->kind, NULL);
	value = con->kind;
	break;
      case KINDS_OF_CONSTRAINT_RANKS: 
        SET_Insert (&rndsetmem, &legend.discrete, (void*) (long) con->rank, NULL);
	value = con->rank;
	break;
      case RESULTS_RT: value = LEN2 (con->R); break;
      case RESULTS_RN: value = con->R[2]; break;
      case RESULTS_R: value = SGN (con->R[2]) * LEN (con->R); break;
      case RESULTS_UT: value = LEN2 (con->U); break;
      case RESULTS_UN: value = con->U[2]; break;
      case RESULTS_U: value = SGN (con->U[2]) * LEN (con->U); break;
      case RESULTS_GAP: value = con->gap; break;
      case RESULTS_MERIT: value = con->merit; break;
      }

      if (value < legend.extents [0]) legend.extents [0] = value;
      if (value > legend.extents [1]) legend.extents [1] = value;

      con->state |= CON_DONE;
    }
  }
  else
  {
    for (force = bod->forces; force; force = force->next)
    {
      if (force->data == NULL) continue;

      value = TMS_Value (force->data, domain->time);

      if (value < legend.extents [0]) legend.extents [0] = value;
      if (value > legend.extents [1]) legend.extents [1] = value;
    }
  }
}

/* update body values */
static void update_body_values (BODY *bod, BODY_DATA *data)
{
  double **vvs, *val, *end, values [7] = {0, 0, 0, 0, 0, 0, 0};
  VALUE_SOURCE *src, *last;
  VALUE_KIND kind;
  short index;
  MAP *item;
  SET *jtem;

  if (legend.entity >= KINDS_OF_BODIES && legend.entity < RESULTS_RT && render_bodies)
  {
    switch (legend.entity)
    {
    case KINDS_OF_BODIES:
      for (val = data->vertex_values, end = val + data->vertex_values_count; val < end; val ++) *val = (double) bod->kind;
      if (bod->kind < legend.extents [0]) legend.extents [0] = bod->kind;
      if (bod->kind > legend.extents [1]) legend.extents [1] = bod->kind;
      SET_Insert (&rndsetmem, &legend.discrete, (void*) (long) bod->kind, NULL);
      break;
    case KINDS_OF_SURFACES:
      for (item = MAP_First (data->surfaces); item; item = MAP_Next (item))
      {
	for (jtem = SET_First (item->data); jtem; jtem = SET_Next (jtem)) val = jtem->data, *val = (double) (long) item->key;
	if ((double) (long) item->key < legend.extents [0]) legend.extents [0] = (double) (long) item->key;
	if ((double) (long) item->key > legend.extents [1]) legend.extents [1] = (double) (long) item->key;
        SET_Insert (&rndsetmem, &legend.discrete, item->key, NULL);
      }
      break;
    case KINDS_OF_VOLUMES:
      for (item = MAP_First (data->volumes); item; item = MAP_Next (item))
      {
	for (jtem = SET_First (item->data); jtem; jtem = SET_Next (jtem)) val = jtem->data, *val = (double) (long) item->key;
	if ((double) (long) item->key < legend.extents [0]) legend.extents [0] = (double) (long) item->key;
	if ((double) (long) item->key > legend.extents [1]) legend.extents [1] = (double) (long) item->key;
        SET_Insert (&rndsetmem, &legend.discrete, item->key, NULL);
      }
      break;
    case KINDS_OF_BODY_RANKS:
      for (val = data->vertex_values, end = val + data->vertex_values_count; val < end; val ++) *val = (double) bod->rank;
      if (bod->rank < legend.extents [0]) legend.extents [0] = bod->rank;
      if (bod->rank > legend.extents [1]) legend.extents [1] = bod->rank;
      SET_Insert (&rndsetmem, &legend.discrete, (void*) (long) bod->rank, NULL);
      break;
    case RESULTS_DX:
    case RESULTS_DY:
    case RESULTS_DZ:
      kind = VALUE_DISPLACEMENT;
      index = legend.entity - RESULTS_DX;
      break;
    case RESULTS_VX:
    case RESULTS_VY:
    case RESULTS_VZ:
      kind = VALUE_VELOCITY;
      index = legend.entity - RESULTS_VX;
      break;
    case RESULTS_SX:
    case RESULTS_SY:
    case RESULTS_SZ:
    case RESULTS_SXY:
    case RESULTS_SXZ:
    case RESULTS_SYZ:
      kind = VALUE_STRESS;
      index = legend.entity - RESULTS_SX;
      break;
    case RESULTS_MISES:
      kind = VALUE_MISES;
      index = 0;
      break;
    }

    if (legend.entity >= RESULTS_DX)
    {
      src = data->value_sources;
      last = src + data->values_count;
      val = data->values;

      if (bod->kind == FEM)
      {
	if (bod->msh)
	{
	  for (; src < last; src ++, val ++)
	  {
	    FEM_Element_Point_Values (bod, src->epn->ele, src->epn->pnt, kind, values);
	    *val = values [index];
	    if (*val < legend.extents [0])  legend.extents [0] = *val;
	    if (*val > legend.extents [1])  legend.extents [1] = *val;
	  }
	}
	else
	{
	  for (; src < last; src ++, val ++)
	  {
	    FEM_Cur_Node_Values (bod, src->pnt, kind, values);
	    *val = values [index];
	    if (*val < legend.extents [0])  legend.extents [0] = *val;
	    if (*val > legend.extents [1])  legend.extents [1] = *val;
	  }
	}
      }
      else if (bod->kind != OBS && !(bod->kind == RIG && legend.entity >= RESULTS_SX))
      {
	for (; src < last; src ++, val ++)
	{
	  BODY_Point_Values (bod, src->pnt, kind, values);
	  *val = values [index];
	  if (*val < legend.extents [0])  legend.extents [0] = *val;
	  if (*val > legend.extents [1])  legend.extents [1] = *val;
	}
      }

      /* update vertex values */
      for (val = data->vertex_values, vvs = data->vertex_value_sources,
	   end = val + data->vertex_values_count; val < end; vvs ++, val ++) *val = **vvs;
    }

    for (SPHERE **sph = data->spheres, **end = sph + data->spheres_count; sph < end; sph ++)
    {
      values [0] = sphere_value (bod, *sph);
      if (values[0] < legend.extents [0])  legend.extents [0] = values[0];
      if (values[0] > legend.extents [1])  legend.extents [1] = values[0];

      switch (legend.entity)
      {
      case KINDS_OF_SURFACES:
        SET_Insert (&rndsetmem, &legend.discrete, (void*) (long) (*sph)->surface, NULL);
	break;
      case KINDS_OF_VOLUMES:
        SET_Insert (&rndsetmem, &legend.discrete, (void*) (long) (*sph)->volume, NULL);
	break;
      }
    }

    for (ELLIP **eli = data->ellips, **end = eli + data->ellips_count; eli < end; eli ++)
    {
      values [0] = ellip_value (bod, *eli);
      if (values[0] < legend.extents [0])  legend.extents [0] = values[0];
      if (values[0] > legend.extents [1])  legend.extents [1] = values[0];

      switch (legend.entity)
      {
      case KINDS_OF_SURFACES:
        SET_Insert (&rndsetmem, &legend.discrete, (void*) (long) (*eli)->surface, NULL);
	break;
      case KINDS_OF_VOLUMES:
        SET_Insert (&rndsetmem, &legend.discrete, (void*) (long) (*eli)->volume, NULL);
	break;
      }
    }


    switch (bod->kind)
    {
    case OBS:
      data->values_updated = legend.entity < RESULTS_DX;
      break;
    case RIG:
      data->values_updated = legend.entity < RESULTS_SX;
      break;
    default:
      data->values_updated = 1;
      break;
    }
  }
  else if (legend.entity) update_body_constraint_or_force_values (bod);
}

/* update body rendering data */
static void update_body_data (BODY *bod, BODY_DATA *data)
{
  double **vsr, **nsr, **lsr, **end, *val, *tail;
  GLfloat *ver, *v, *nor, *n, *col, *c, *lin, *l;

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->lines);
  lin = glMapBufferARB (GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB);
#else
  lin = data->lin;
#endif

  for (lsr = data->line_sources,
       end = lsr + data->lines_count * 2,
       l = lin; lsr < end; lsr ++, l += 3)
  {
    COPY (*lsr, l);
  }

  if (data->flags & WIREFRAME)
  {
    col = lin + data->lines_count * 6;
  }
  else
  {
#if VBO
    glUnmapBufferARB (GL_ARRAY_BUFFER_ARB);

    glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);
    ver = glMapBufferARB (GL_ARRAY_BUFFER_ARB, GL_WRITE_ONLY_ARB);
#else
    ver = data->ver;
#endif

    nor = ver + data->triangles_count * 9;
    col = nor + data->triangles_count * 9;

    for (vsr = data->vertex_sources,
	 nsr = data->normal_sources,
	 end = vsr + data->triangles_count * 3,
	 v = ver, n = nor; vsr < end;
	 vsr ++, nsr ++, val ++, v += 3, n += 3)
    {
      COPY (*vsr, v);
      COPY (*nsr, n);
    }
  }

  if (data->values_updated)
  {
    for (val = data->vertex_values, tail = val + data->vertex_values_count, c = col; val < tail;  val ++, c += 3) value_to_color (*val, c);
  }
  else for (c = col, l = c + data->vertex_values_count * 3; c < l;  c += 3) { COPY (neutral_color, c); }

#if VBO
  glUnmapBufferARB (GL_ARRAY_BUFFER_ARB);
#endif

  if (data->spheres_count)
  {
    SPHERE **sph;

    if (data->values_updated)
    {
      for (c = data->sphere_colors, sph = data->spheres, v = c + data->spheres_count * 3; c < v; c += 3, sph ++)
      {
	value_to_color (sphere_value (bod, *sph), c);
      }
    }
    else
    {
      for (c = data->sphere_colors, v = c + data->spheres_count * 3; c < v; c += 3)
      {
	COPY (neutral_color, c);
      }
    }
  }

  if (data->ellips_count)
  {
    ELLIP **eli;

    if (data->values_updated)
    {
      for (c = data->ellip_colors, eli = data->ellips, v = c + data->ellips_count * 3; c < v; c += 3, eli ++)
      {
	value_to_color (ellip_value (bod, *eli), c);
      }
    }
    else
    {
      for (c = data->ellip_colors, v = c + data->ellips_count * 3; c < v; c += 3)
      {
	COPY (neutral_color, c);
      }
    }
  }

  if (data->rough) update_body_data (bod, data->rough);
}

/* update cut */
static void update_cuts (void)
{
  CUT_DATA *cut, *next, *list;

  /* delete current Euler cuts */
  for (list = NULL, cut = cuts; cut; cut = next)
  {
    next = cut->next;

    if (cut->kind == EULER)
    {
      BODY_DATA *data = cut->bod->rendering;
      data->flags &= ~SEETHROUGH;
      free (cut->tri);
      free (cut);
    }
    else cut->next = list, list = cut; /* put aside Lagrange cuts */
  }

  /* update Lagrange cuts */
  for (cut = list; cut; cut = cut->next)
  {
    for (int i = 0; i < cut->n; i ++)
    {
      double *X = &cut->ref[3*i], *x = &cut->cur[3*i];
      SGP *sgp = &cut->sgp [i];
      BODY_Cur_Point (cut->bod, sgp, X, x);
    }
  }

  /* create new Euler cuts */
  for (SET *jtem = SET_First (euler_cuts); jtem; jtem = SET_Next (jtem))
  {
    double *point = jtem->data, *normal = point + 3;
    for (SET *item = SET_First (selection->set); item; item = SET_Next (item))
    {
      BODY *bod = item->data;
      BODY_DATA *data = bod->rendering;
      double *ref, *cur;
      CUT_DATA *cut;
      SGP *sgp;
      TRI *tri;
      int m, n;

      if (bod->kind == OBS) continue;
      if (data->flags & HIDDEN) continue;

      tri = SHAPE_Cut (bod->shape, point, normal, &m, bod, (MOTION)BODY_Ref_Point, &sgp, &ref, &cur, &n);

      if (tri)
      {
	data->flags |= SEETHROUGH; /* make the body transparent if it was cut */
	ERRMEM (cut = MEM_CALLOC (sizeof (CUT_DATA) + sizeof (double [n])));
	cut->val = (double*) (cut+1);
	cut->kind = cut_kind;
	COPY (point, cut->point);
	COPY (normal, cut->normal);
	cut->bod = bod;
	cut->tri = tri;
	cut->sgp = sgp;
	cut->ref = ref;
	cut->cur = cur;
	cut->m = m;
	cut->n = n;
	cut->next = list;
	list = cut;
      }
    }
  }

  /* reset list */
  cuts = list;
}

/* update cut values */
static void update_cut_values (CUT_DATA *cut)
{
  if (legend.entity >= RESULTS_DX && legend.entity < RESULTS_RT)
  {
    for (int i = 0; i < cut->n; i ++)
    {
      double val = point_value (cut->bod, cut->sgp[i].shp, cut->sgp[i].gobj, &cut->ref [3*i]);
      if (val < legend.extents [0])  legend.extents [0] = val;
      if (val > legend.extents [1])  legend.extents [1] = val;
      cut->val [i] = val;
    }
  }
  else
  {
    for (int i = 0; i < cut->n; i ++) cut->val [i] = 0.0;
  }

  switch (cut->bod->kind)
  {
  case OBS:
    cut->values_updated = legend.entity < RESULTS_DX;
    break;
  case RIG:
    cut->values_updated = legend.entity < RESULTS_SX;
    break;
  default:
    cut->values_updated = 1;
    break;
  }
}

/* update cuts values */
static void update_cuts_values (void)
{
  /* update cut point values */
  for (CUT_DATA *cut = cuts; cut; cut = cut->next)
  {
    update_cut_values (cut);
  }
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

/* render ellipsoid triangles */
static void render_ellip_triangles (double *center, double *sca, double *rot, GLfloat color [3])
{
  glMatrixMode (GL_MODELVIEW_MATRIX);
  glPushMatrix ();
    glTranslated (center[0], center[1], center[2]);
    GLfloat m [16] = {rot [0], rot [1], rot [2], 0,
                      rot [3], rot [4], rot [5], 0,
		      rot [6], rot [7], rot [8], 0,
		      0      , 0      , 0      , 1};
    glMultMatrixf (m);
    glScaled (sca [0], sca [1], sca [2]);
    glColor3fv (color);
    glutSolidSphere (1.0, 12, 12);
  glPopMatrix ();
}

/* render sphere triangles for selection */
static void selection_render_sphere_triangles (double *center, double radius)
{
  glMatrixMode (GL_MODELVIEW_MATRIX);
  glPushMatrix ();
    glTranslated (center[0], center[1], center[2]);
    glutSolidSphere (radius, 12, 12);
  glPopMatrix ();
}

/* render ellipsoid triangles for selection */
static void selection_render_ellip_triangles (double *center, double *sca, double *rot)
{
  glMatrixMode (GL_MODELVIEW_MATRIX);
  glPushMatrix ();
    glTranslated (center[0], center[1], center[2]);
    GLfloat m [16] = {rot [0], rot [1], rot [2], 0,
                      rot [3], rot [4], rot [5], 0,
		      rot [6], rot [7], rot [8], 0,
		      0      , 0      , 0      , 1};
    glMultMatrixf (m);
    glScaled (sca [0], sca [1], sca [2]);
    glutSolidSphere (1.0, 12, 12);
  glPopMatrix ();
}

/* get legend caption */
static char* legend_caption ()
{
  switch (legend.entity)
  {
  case KINDS_OF_CONSTRAINTS: return "CONSTRAINT KINDS";
  case KINDS_OF_CONSTRAINT_RANKS: return "CONSTRAINT RANKS";
  case KINDS_OF_FORCES: return "FORCE KINDS";
  case KINDS_OF_BODIES: return "BODY KINDS";
  case KINDS_OF_BODY_RANKS: return "BODY RANKS";
  case KINDS_OF_SURFACES: return "SURFACE KINDS";
  case KINDS_OF_VOLUMES: return "VOLUME KINDS";
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
  case RESULTS_RT: return "RT";
  case RESULTS_RN: return "RN";
  case RESULTS_R: return "R";
  case RESULTS_UT: return "UT";
  case RESULTS_UN: return "UN";
  case RESULTS_U: return "U";
  case RESULTS_GAP: return "GAP";
  case RESULTS_MERIT: return "MERIT";
  }

  return NULL;
}

/* get legend caption widtth */
static int legend_caption_width ()
{
  return GLV_Print_Width (LEGEND_FONT, legend_caption ()) + 9;
}

/* get legend width */
static int legend_width ()
{
  int i, j;

  i = legend.range / LEGEND_ROWS;
  if (legend.range % LEGEND_ROWS) i ++;

  if (legend.discrete) i *= LEGEND_WIDTH_DISC;
  else i *= LEGEND_WIDTH_CONT;

  j = legend_caption_width ();

  return MAX (i, j);
}

/* get legend height */
static int legend_height ()
{
  int i;

  i = MIN (legend.range, LEGEND_ROWS) + 1;

  return MAX (i, 1) * 16 + 6;
}

/* get legend value string */
char *legend_value_string (void *data)
{
  if (legend.entity == KINDS_OF_CONSTRAINTS)
  {
    switch ((int)data)
    {
      case CONTACT: return "CNT";
      case FIXPNT: return "PNT";
      case FIXDIR: return "DIR";
      case VELODIR: return "VEL";
      case RIGLNK: return "LNK";
      case GLUE: return "GLU";
      default: return "???";
    }
  }
  else if (legend.entity == KINDS_OF_BODIES)
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

/* render legend */
static void legend_render ()
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
  glRecti (v[0] + 3, v[1] + 3, 3 + legend_caption_width (), 19);
  glColor3f (0, 0, 0);
  GLV_Print (v[0] + 6, v[1] + 6, 1, LEGEND_FONT, "%s", legend_caption ());

  if (legend.discrete)
  {
    glPushMatrix ();
    glTranslated (3, 3, 0);
    for (i = 1, j = 0, item = SET_First (legend.discrete); item; item = SET_Next (item))
    {
      value_to_color ((double)(int)item->data, color);
      glColor3fv (color);
      glRecti (v[0] + j * LEGEND_WIDTH_DISC, v[1] + i * 16, v[0] + j * LEGEND_WIDTH_DISC + 16, v[1] + i * 16 + 16);
      if ((str = legend_value_string (item->data)))
      {
	glColor3f (1, 1, 1);
	l = GLV_Print_Width (LEGEND_FONT, str) + 5;
	glRecti (v[0] + j * LEGEND_WIDTH_DISC + 16, v[1] + i * 16, v[0] + j * LEGEND_WIDTH_DISC + 16 + l, v[1] + (i+1) * 16);
        glColor3f (0, 0, 0);
        GLV_Print (v[0] + j * LEGEND_WIDTH_DISC + 18, v[1] + i * 16 + 3, 1, LEGEND_FONT, str);
      }
      else
      {
	glColor3f (1, 1, 1);
	l = GLV_Print_Width (LEGEND_FONT, "%d", (int)item->data) + 5;
	glRecti (v[0] + j * LEGEND_WIDTH_DISC + 16, v[1] + i * 16, v[0] + j * LEGEND_WIDTH_DISC + 16 + l, v[1] + (i+1) * 16);
        glColor3f (0, 0, 0);
	GLV_Print (v[0] + j * LEGEND_WIDTH_DISC + 18, v[1] + i * 16 + 3, 1, LEGEND_FONT, "%d", (int)item->data);
      }
      if (i ++ == LEGEND_ROWS) { i = 1; j ++; }
    }
    glPopMatrix ();
  }
  else
  {
    if (legend.extents [0] < DBL_MAX)
    {
      step = (legend.extents [1] - legend.extents [0]) / (double) (legend.range - 1);
      value = legend.extents [0];

      glPushMatrix ();
      glTranslated (3, 3, 0);
      for (i = 1, j = k = 0; k < legend.range; k ++, value += step)
      {
	value_to_color (value, color);
	glColor3fv (color);
	glRecti (v[0] + j * LEGEND_WIDTH_CONT, v[1] + i * 16, v[0] + j * LEGEND_WIDTH_CONT + 16, v[1] + i * 16 + 16);
	glColor3f (1, 1, 1);
	l = GLV_Print_Width (LEGEND_FONT, "%.2e", value) + 5;
	glRecti (v[0] + j * LEGEND_WIDTH_CONT + 16, v[1] + i * 16, v[0] + j * LEGEND_WIDTH_CONT + 16 + l, v[1] + (i+1) * 16);
	glColor3f (0, 0, 0);
	GLV_Print (v[0] + j * LEGEND_WIDTH_CONT + 18, v[1] + i * 16 + 3, 1, LEGEND_FONT, "%.2e", value);
	if (i ++ == LEGEND_ROWS) { i = 1; j ++; }
      }
      glPopMatrix ();
    }
  }

  glEnable (GL_LIGHTING);
  glEnable (GL_DEPTH_TEST);
}

/* disable_legend */
static void legend_disable ()
{
  if (legend.window)
  {
    GLV_Close_Viewport (legend.window);

    legend.window = 0;
  }

  if (legend.entity)
  {
    for (BODY *bod = domain->bod; bod; bod = bod->next)
    {
      BODY_DATA *data = bod->rendering;
      data->values_updated = 0;
    }

    legend.preventity = legend.entity;

    legend.entity = 0;
  }
}

/* enable legend */
static void legend_enable ()
{
  if (legend.discrete) legend.range = SET_Size (legend.discrete);
  else
  {
    if (legend.extents [0] == legend.extents [1]) legend.range = 1;
    else legend.range = 8;
  }

  legend.constant = legend.range <= 1 ? 1.0 : 0.0;

  if (legend.window) GLV_Close_Viewport (legend.window);
  legend.window = GLV_Open_Viewport (0, 0, legend_width (), legend_height (), 0, legend_render);
}

/* pop previous legend entity */
static void legend_pop ()
{
  if ((legend.entity = legend.preventity)) update ();
}

/* render 3 points */
static void render_3_points (double *a, double *b, double *c)
{
  glBegin (GL_POINTS);
    glVertex3dv (a);
    glVertex3dv (b);
    glVertex3dv (c);
  glEnd ();
}

/* render body triangles */
static void render_body_triangles (BODY *bod, short skip)
{
  BODY_DATA *data;

  if (bod->rendering == NULL) bod->rendering = create_body_data (bod, WIREFRAME_FLAG());

  data = bod->rendering;

  if (bod == picked_body ||           /* do not render a picked body */
      data->flags & skip) return;

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);
#endif

  glEnableClientState (GL_VERTEX_ARRAY);
  glEnableClientState (GL_NORMAL_ARRAY);
  glEnableClientState (GL_COLOR_ARRAY);

#if VBO
  glVertexPointer (3, GL_FLOAT, 0, 0);
  glNormalPointer (GL_FLOAT, 0, (void*) (data->triangles_count * sizeof (GLfloat) * 9));
  glColorPointer (3, GL_FLOAT, 0, (void*) (data->triangles_count * sizeof (GLfloat) * 18));
#else
  glVertexPointer (3, GL_FLOAT, 0, data->ver);
  glNormalPointer (GL_FLOAT, 0, data->ver +  (data->triangles_count * 9));
  glColorPointer (3, GL_FLOAT, 0, data->ver + (data->triangles_count * 18));
#endif

  glDrawArrays (GL_TRIANGLES, 0, data->triangles_count * 3);

  glDisableClientState (GL_VERTEX_ARRAY);
  glDisableClientState (GL_NORMAL_ARRAY);
  glDisableClientState (GL_COLOR_ARRAY);

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);
#endif

  SPHERE **sph, **end0;
  ELLIP **eli, **end1;
  GLfloat *col;

  for (sph = data->spheres, end0 = sph + data->spheres_count, col = data->sphere_colors; sph < end0; sph ++, col += 3)
    render_sphere_triangles ((*sph)->cur_center, (*sph)->cur_radius, col);

  for (eli = data->ellips, end1 = eli + data->ellips_count, col = data->ellip_colors; eli < end1; eli ++, col += 3)
    render_ellip_triangles ((*eli)->cur_center, (*eli)->cur_sca, (*eli)->cur_rot, col);
}

/* render body triangles without colors */
static void render_body_triangles_plain (BODY *bod, short skip)
{
  BODY_DATA *data = bod->rendering;

  if (bod == picked_body ||           /* do not render a picked body */
      data->flags & skip) return;

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);
#endif

  glEnableClientState (GL_VERTEX_ARRAY);
  glEnableClientState (GL_NORMAL_ARRAY);

#if VBO
  glVertexPointer (3, GL_FLOAT, 0, 0);
  glNormalPointer (GL_FLOAT, 0, (void*) (data->triangles_count * sizeof (GLfloat) * 9));
#else
  glVertexPointer (3, GL_FLOAT, 0, data->ver);
  glNormalPointer (GL_FLOAT, 0, data->ver + (data->triangles_count * 9));
#endif

  glDrawArrays (GL_TRIANGLES, 0, data->triangles_count * 3);

  glDisableClientState (GL_VERTEX_ARRAY);
  glDisableClientState (GL_NORMAL_ARRAY);

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);
#endif

  for (SPHERE **sph = data->spheres, **end = sph + data->spheres_count; sph < end; sph ++)
    selection_render_sphere_triangles ((*sph)->cur_center, (*sph)->cur_radius);

  for (ELLIP **eli = data->ellips, **end = eli + data->ellips_count; eli < end; eli ++)
    selection_render_ellip_triangles ((*eli)->cur_center, (*eli)->cur_sca, (*eli)->cur_rot);
}

/* render body lines */
static void render_body_lines (BODY *bod, short skip)
{
  BODY_DATA *data;

  if (bod->rendering == NULL) bod->rendering = create_body_data (bod, WIREFRAME_FLAG());

  data = bod->rendering;

  if (data->flags & skip) return;

  int wire = data->flags & WIREFRAME;

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->lines);
#endif

  glEnableClientState (GL_VERTEX_ARRAY);
  if (wire) glEnableClientState (GL_COLOR_ARRAY);

#if VBO
  glVertexPointer (3, GL_FLOAT, 0, 0);
  if (wire) glColorPointer (3, GL_FLOAT, 0, (void*) (data->lines_count * sizeof (GLfloat) * 6));
#else
  glVertexPointer (3, GL_FLOAT, 0, data->lin);
  if (wire) glColorPointer (3, GL_FLOAT, 0, data->lin + (data->lines_count * 6));
#endif

  glDrawArrays (GL_LINES, 0, data->lines_count * 2);

  glDisableClientState (GL_VERTEX_ARRAY);
  if (wire) glDisableClientState (GL_COLOR_ARRAY);

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);
#endif

  for (SPHERE **sph = data->spheres, **end = sph + data->spheres_count; sph < end; sph ++)
    render_3_points ((*sph)->cur_point[0], (*sph)->cur_point[1], (*sph)->cur_point[2]);

  for (ELLIP **eli = data->ellips, **end = eli + data->ellips_count; eli < end; eli ++)
    render_3_points ((*eli)->cur_point[0], (*eli)->cur_point[1], (*eli)->cur_point[2]);
}

/* render body for selection */
static void selection_render_body (BODY *bod)
{
  BODY_DATA *data = bod->rendering;

  if (data->flags & HIDDEN) return;

  if (data->flags & WIREFRAME)
  {
    glLineWidth (4.0);

#if VBO
    glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->lines);
#endif

    glEnableClientState (GL_VERTEX_ARRAY);

#if VBO
    glVertexPointer (3, GL_FLOAT, 0, 0);
#else
    glVertexPointer (3, GL_FLOAT, 0, data->lin);
#endif

    glDrawArrays (GL_LINES, 0, data->lines_count * 2);

    glDisableClientState (GL_VERTEX_ARRAY);

#if VBO
    glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);
#endif

    glLineWidth (1.0);

    glPointSize (4.0);

    for (SPHERE **sph = data->spheres, **end = sph + data->spheres_count; sph < end; sph ++)
      render_3_points ((*sph)->cur_point[0], (*sph)->cur_point[1], (*sph)->cur_point[2]);

    for (ELLIP **eli = data->ellips, **end = eli + data->ellips_count; eli < end; eli ++)
      render_3_points ((*eli)->cur_point[0], (*eli)->cur_point[1], (*eli)->cur_point[2]);

    glPointSize (1.0);
  }
  else
  {
#if VBO
    glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);
#endif

    glEnableClientState (GL_VERTEX_ARRAY);
    glEnableClientState (GL_NORMAL_ARRAY);

#if VBO
    glVertexPointer (3, GL_FLOAT, 0, 0);
    glNormalPointer (GL_FLOAT, 0, (void*) (data->triangles_count * sizeof (GLfloat) * 9));
#else
    glVertexPointer (3, GL_FLOAT, 0, data->ver);
    glNormalPointer (GL_FLOAT, 0, data->ver + (data->triangles_count * 9));
#endif

    glDrawArrays (GL_TRIANGLES, 0, data->triangles_count * 3);

    glDisableClientState (GL_VERTEX_ARRAY);
    glDisableClientState (GL_NORMAL_ARRAY);

#if VBO
    glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);
#endif

    for (SPHERE **sph = data->spheres, **end = sph + data->spheres_count; sph < end; sph ++)
      selection_render_sphere_triangles ((*sph)->cur_center, (*sph)->cur_radius);

    for (ELLIP **eli = data->ellips, **end = eli + data->ellips_count; eli < end; eli ++)
      selection_render_ellip_triangles ((*eli)->cur_center, (*eli)->cur_sca, (*eli)->cur_rot);
  }
}

/* render rough mesh */
static void render_rough_mesh (BODY *bod)
{
  BODY_DATA *data = bod->rendering,
	    *rough = data->rough;

  if (data->flags & HIDDEN) return;

  if (!rough)
  {
    SHAPE shape = {SHAPE_MESH, bod->msh, NULL};
    BODY body = bod [0];
    body.shape = &shape;
    body.msh = NULL;
    data->rough = create_body_data (&body, 0);
    rough = data->rough;
  }

  glColor4f (0.0, 0.0, 0.0, 0.5);

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, rough->lines);
#endif

  glEnableClientState (GL_VERTEX_ARRAY);

#if VBO
  glVertexPointer (3, GL_FLOAT, 0, 0);
#else
  glVertexPointer (3, GL_FLOAT, 0, rough->lin);
#endif

  glDrawArrays (GL_LINES, 0, rough->lines_count * 2);

  glDisableClientState (GL_VERTEX_ARRAY);

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);
#endif

  glColor4f (0.9, 0.9, 0.9, 0.5);

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, rough->triangles);
#endif

  glEnableClientState (GL_VERTEX_ARRAY);
  glEnableClientState (GL_NORMAL_ARRAY);

#if VBO 
  glVertexPointer (3, GL_FLOAT, 0, 0);
  glNormalPointer (GL_FLOAT, 0, (void*) (rough->triangles_count * sizeof (GLfloat) * 9));
#else
  glVertexPointer (3, GL_FLOAT, 0, rough->ver);
  glNormalPointer (GL_FLOAT, 0, rough->ver + (rough->triangles_count * 9));
#endif

  glDrawArrays (GL_TRIANGLES, 0, rough->triangles_count * 3);

  glDisableClientState (GL_VERTEX_ARRAY);
  glDisableClientState (GL_NORMAL_ARRAY);

#if VBO
  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);
#endif
}

/* render contact */
static void render_contact (CON *con, GLfloat color [3])
{
  double other [3],
	 scal = GLV_Minimal_Extent() * 0.03;

  glColor3fv (color);

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
static void render_fixpnt (CON *con, GLfloat color [3])
{
  glColor3fv (color);
  glPointSize (4.0);
  glBegin (GL_POINTS);
    glVertex3dv (con->point);
  glEnd ();
  glPointSize (1.0);
}

/* render fixed direction */
static void render_fixdir (CON *con, GLfloat color [3])
{
  double other [3],
	 scal = GLV_Minimal_Extent() * 0.03;

  glColor3fv (color);

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
static void render_velodir (CON *con, GLfloat color [3])
{
  double other [3],
	 scal = GLV_Minimal_Extent() * 0.03;

  glColor3fv (color);

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
static void render_riglnk (CON *con, GLfloat width, GLfloat color [3])
{
  double other [3];

  SUB (con->point, RIGLNK_VEC (con->Z), other);

  glColor3fv (color);

  glPointSize (3.0);
  glBegin (GL_POINTS);
    glVertex3dv (con->point);
    glVertex3dv (other);
  glEnd ();
  glPointSize (1.0);
  
  glLineWidth (width);
  glBegin (GL_LINES);
    glVertex3dv (con->point);
    glVertex3dv (other);
  glEnd ();
  glLineWidth (1.0);
}

/* render glue constraint */
static void render_glue (CON *con, GLfloat color [3])
{
  glColor3fv (color);
  glPointSize (4.0);
  glBegin (GL_POINTS);
    glVertex3dv (con->point);
  glEnd ();
  glPointSize (1.0);
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
  int i, div = 2;

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
  angle = ALG_PI / (double) div;
  SCALE (d, angle);
  EXPMAP (d, R);
  o = 0.075 * l;
  s = 0.150 * l;

  ADDMUL (p, o, t, a);
  ADDMUL (r, o, t, b);

  div *= 2;

  for (i = 0; i < div; i ++)
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
  for (i = 0; i < div; i ++)
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
  for (i = 0; i < div; i ++)
  {
    SUB (a, p, x);
    TVMUL (R, x, y);
    ADD (p, y, b);
    glVertex3dv (b);
    COPY (b, a);
  }
  glEnd ();
}

/* render tangential vector */
static void render_vt (CON *con, double *V, double value)
{
  GLfloat color [3];
  double r [3],
	 other [3],
	 ext = GLV_Minimal_Extent() * arrow_factor,
	 eps, len;

  value_to_color (value, color);
  COPY (con->base, r);
  SCALE (r, V[0]);
  ADDMUL (r, V[1], con->base+3, r);
  len = LEN (r);
  if (len == 0.0) len = 1.0;
  eps = (ext  / len) * (0.5 + (value - legend.extents[0]) / (legend.extents[1] - legend.extents[0] + 1.0));
  ADDMUL (con->point, -eps, r, other);

  glColor3fv (color);
  arrow3d (other, con->point);
}

/* render normal vector */
static void render_vn (CON *con, double *V, double value)
{
  GLfloat color [3];
  double r [3],
	 other [3],
	 ext = GLV_Minimal_Extent() * arrow_factor,
	 eps, len;

  value_to_color (value, color);
  COPY (con->base + 6, r);
  SCALE (r, V[2]);
  len = LEN (r);
  if (len == 0.0) len = 1.0;
  eps = (ext  / len) * (0.5 + (value - legend.extents[0]) / (legend.extents[1] - legend.extents[0] + 1.0));
  ADDMUL (con->point, -eps, r, other);

  glColor3fv (color);
  arrow3d (other, con->point);
}

/* render vector */
static void render_v (CON *con, double *V, double value)
{
  GLfloat color [3];
  double r [3],
	 other [3],
	 *base = con->base,
	 ext = GLV_Minimal_Extent() * arrow_factor,
	 eps,
	 len;

  value_to_color (value, color);
  NVMUL (base, V, r);
  len = LEN (r);
  if (len == 0.0) len = 1.0;
  eps = (ext  / len) * (0.5 + (value - legend.extents[0]) / (legend.extents[1] - legend.extents[0] + 1.0));
  ADDMUL (con->point, -eps, r, other);

  glColor3fv (color);
  arrow3d (other, con->point);
}

/* draw tube from p to q */
static void tube3d (double *p, double *q)
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
  int i, div = 2;

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
  ADDMUL (p, l, d, r);
  angle = ALG_PI / (double) div;
  SCALE (d, angle);
  EXPMAP (d, R);
  o = 0.075 * l;
  s = 0.150 * l;

  ADDMUL (p, o, t, a);
  ADDMUL (r, o, t, b);

  div *= 2;

  for (i = 0; i < div; i ++)
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

  ADDMUL (q, o, t, a);
  glNormal3dv (n);
  glBegin (GL_POLYGON);
  glVertex3dv (a);
  for (i = 0; i < div; i ++)
  {
    SUB (a, q, x);
    NVMUL (R, x, y);
    ADD (q, y, b);
    glVertex3dv (b);
    COPY (b, a);
  }
  glEnd ();

  ADDMUL (p, o, t, a);
  glNormal3dv (n);
  glBegin (GL_POLYGON);
  glVertex3dv (a);
  for (i = 0; i < div; i ++)
  {
    SUB (a, p, x);
    TVMUL (R, x, y);
    ADD (p, y, b);
    glVertex3dv (b);
    COPY (b, a);
  }
  glEnd ();
}

/* render constraint scalar */
static void render_scalar (CON *con, double scalar)
{
  GLfloat color [3];
  double r [3],
	 first [3],
	 other [3],
	 ext = GLV_Minimal_Extent() * arrow_factor,
	 eps,
	 len;

  value_to_color (scalar, color);
  COPY (con->base+6, r);
  SCALE (r, scalar);
  len = LEN (r);
  if (len == 0.0) len = 1.0;
  eps = (ext  / len) * (1.0 + (len - legend.extents[0]) / (legend.extents[1] - legend.extents[0] + 1.0));
  ADDMUL (con->point, -0.5 * eps, r, first);
  ADDMUL (first,  eps, r, other);

  glColor3fv (color);
  tube3d (first, other);
}

/* render force */
static void render_force (BODY *bod, FORCE *force, GLfloat color [3])
{
  double r [3],
	 other [3],
	 point [3],
	 ext = GLV_Minimal_Extent() * arrow_factor,
	 eps,
	 len;

  if (force->kind & PRESSURE)
  {
    MESH *msh = bod->shape->data;
    FACE *fac = msh->faces;
    int surfid = force->surfid;
    double (*cur) [3] = msh->cur_nodes;

    for (; fac; fac = fac->n)
    {
      if (fac->surface == surfid)
      {
	double *a = cur [fac->nodes [0]],
	       *b = cur [fac->nodes [1]],
	       *c = cur [fac->nodes [2]],
	       *d = fac->type == 4 ? cur [fac->nodes [3]] : NULL;

	if (fac->type == 4)
	{
	  glBegin (GL_QUADS);
          glColor3fv (color);
	  glNormal3dv (fac->normal);
	  glVertex3dv (a);
	  glVertex3dv (b);
	  glVertex3dv (c);
	  glVertex3dv (d);
	  glEnd ();
	}
	else
	{
	  glBegin (GL_TRIANGLES);
          glColor3fv (color);
	  glNormal3dv (fac->normal);
	  glVertex3dv (a);
	  glVertex3dv (b);
	  glVertex3dv (c);
	  glEnd ();
	}
      }
    }
  }
  else
  {
    if (bod->kind == FEM)
    {
      SGP *sgp;
      int n;

      if ((n = SHAPE_Sgp (bod->sgp, bod->nsgp, force->ref_point)) < 0) return; /* TODO: optimize */
      sgp = &bod->sgp [n];
      BODY_Cur_Point (bod, sgp, force->ref_point, point); /* TODO: optimize */
    }
    else if (force->kind & TORQUE)
    {
      BODY_Cur_Point (bod, NULL, bod->ref_center, point);
    }
    else BODY_Cur_Point (bod, NULL, force->ref_point, point);

    COPY (force->direction, r);
    len = TMS_Value (force->data, domain->time);
    if (len == 0.0) len = 1.0;
    SCALE (r, len);
    len = LEN (r);
    eps = (ext  / len) * (1.0 + (len - legend.extents[0]) / (legend.extents[1] - legend.extents[0] + 1.0));
    ADDMUL (point, -eps, r, other);

    glColor3fv (color);
    arrow3d (other, point);
  }
}

/* render body set constraints or forces */
static void render_body_set_constraints_or_forces (SET *set)
{
  GLfloat color [3];
  SET *item, *jtem;
  BODY_DATA *data;
  FORCE *force;
  BODY *bod;
  CON *con;

  for (con = domain->con; con; con = con->next) con->state &= ~CON_DONE; /* all undone */

  for (item = SET_First (set); item; item = SET_Next (item))
  {
    bod = item->data;
    data = bod->rendering;

    if (data->flags & HIDDEN) continue;

    if (legend.entity != KINDS_OF_FORCES) /* redner constraint */
    {
      for (jtem = SET_First (bod->con); jtem; jtem = SET_Next (jtem))
      {
	con = jtem->data;

	if (con->state & CON_DONE) continue;

	switch (legend.entity)
	{
	case KINDS_OF_CONSTRAINTS:
	case KINDS_OF_CONSTRAINT_RANKS:

	  if (legend.entity == (short) KINDS_OF_CONSTRAINTS) value_to_color (con->kind, color);
	  else value_to_color (con->rank, color);

	  switch (con->kind)
	  {
	    case CONTACT: render_contact (con, color); break;
	    case FIXPNT: render_fixpnt (con, color); break;
	    case FIXDIR: render_fixdir (con, color); break;
	    case VELODIR: render_velodir (con, color); break;
	    case RIGLNK: render_riglnk (con, 2.0, color); break;
	    case GLUE: render_glue (con, color); break;
	  }

	  break;
	case RESULTS_RT: render_vt (con, con->R, LEN2 (con->R)); break;
	case RESULTS_RN: render_vn (con, con->R, con->R[2]); break;
	case RESULTS_R:  render_v (con, con->R, SGN (con->R[2])*LEN (con->R)); break;
	case RESULTS_UT: render_vt (con, con->U, LEN2 (con->U)); break;
	case RESULTS_UN: render_vn (con, con->U, con->U[2]); break;
	case RESULTS_U:  render_v (con, con->U, SGN (con->U[2])*LEN (con->U)); break;
	case RESULTS_GAP:  render_scalar (con, con->gap); break;
	case RESULTS_MERIT:  render_scalar (con, con->merit); break;
	}

	con->state |= CON_DONE;
      }
    }
    else /* render force */
    {
      for (force = bod->forces; force; force = force->next)
      {
	if (force->data == NULL) continue;

	value_to_color (TMS_Value (force->data, domain->time), color);
	render_force (bod, force, color);
      }
    }
  }
}

/* regular rendering of rigid links */
static void render_rigid_links (SET *set, GLfloat *color)
{
  SET *item, *jtem;
  BODY_DATA *data;
  BODY *bod;
  CON *con;

  /* TODO: optimize iteration through all bodies by pre-selecting 
   * TODO: a set of rigid links when initializing a time step */

  for (item = SET_First (set); item; item = SET_Next (item))
  {
    bod = item->data;
    data = bod->rendering;

    if (data->flags & HIDDEN) continue;

    for (jtem = SET_First (bod->con); jtem; jtem = SET_Next (jtem))
    {
      con = jtem->data;

      if (con->kind == RIGLNK) render_riglnk (con, 1.0, color);
    }
  }
}

/* render body set */
static void render_body_set (SET *set)
{
  GLfloat color [4] = {0.0, 0.0, 0.0, 0.4};
  SET *item;

  if (legend_constraint_based ()) /* render constraints or forces over transparent volumes */
  {
    glDisable (GL_LIGHTING);
    render_body_set_constraints_or_forces (set);
    glEnable (GL_LIGHTING);

    if (!render_bodies) return;

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    COPY (neutral_color, color);

    for (item = SET_First (set); item; item = SET_Next (item))
    {
      if (item->data == picked_body) continue;

      glColor4fv (color);
      render_body_lines (item->data, HIDDEN);
      render_body_triangles_plain (item->data, HIDDEN|WIREFRAME);
    }

    glDisable (GL_BLEND);
  }
  else /* regular rendering */
  {
    if (!render_bodies) return;

    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1.0, 1.0);

    for (item = SET_First (set); item; item = SET_Next (item))
    {
      glDisable (GL_LIGHTING);
      glColor3fv (color);
      render_body_lines (item->data, SEETHROUGH|HIDDEN);
      glEnable (GL_LIGHTING);
      render_body_triangles (item->data, SEETHROUGH|HIDDEN|WIREFRAME);
    }

    render_rigid_links (set, color);

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    COPY (neutral_color, color);

    for (item = SET_First (set); item; item = SET_Next (item))
    {
      if (item->data == picked_body) continue;

      if (SEETHROUGH (item->data))
      {
	glColor4fv (color);
	render_body_lines (item->data, HIDDEN|WIREFRAME);
        render_body_triangles_plain (item->data, HIDDEN|WIREFRAME);
      }

      if (ROUGH_MESH (item->data)) render_rough_mesh (item->data);
    }

    glDisable (GL_BLEND);
  }
}

/* render body set for 2D selection */
static void selection_2D_render_body_set (SET *set)
{
  GLfloat color [4];
  BODY *bod;
  SET *item;

  for (item = SET_First (set); item; item = SET_Next (item))
  {
    bod = item->data;
    idtorgba (bod->id, color);
    glColor4fv (color);
    selection_render_body (item->data);
  }
}

/* render point set for 2D selection */
static void selection_2D_render_point_set (BODY_DATA *rendering)
{
  GLfloat color [4];
  VALUE_SOURCE *src;
  int i;

  glPointSize (16.0);
  glBegin (GL_POINTS);
  for (i = 1, src = rendering->value_sources; i <= rendering->values_count; i ++, src ++)
  {
    idtorgba (i, color);
    glColor4fv (color);
    glVertex3dv (src->pnt);
  }
  glEnd ();
  glPointSize (1.0);
}

/* render body set for 3D selection */
static void selection_3D_render_body_set (SET *set)
{
  BODY *bod;
  SET *item;

  for (item = SET_First (set); item; item = SET_Next (item))
  {
    bod = item->data;
    glPushName (bod->id);
    selection_render_body (item->data);
    glPopName ();
  }
}

/* render picked body if any */
static void render_picked_body (void)
{
  GLfloat color [3] = {0.5, 0.5, 0.5};

  if (picked_body)
  {
    switch (tool_mode)
    {
    case TOOLS_TRANSPARENT: color [0] = 1.0; break;
    case TOOLS_HIDE: color [1] = 1.0; break;
    case TOOLS_ROUGH_MESH: color [2] = 1.0; break;
    }

    glDisable (GL_LIGHTING);
    if (SEETHROUGH (picked_body))
    {
      glColor3f (0., 0., 0.);
      render_body_lines (picked_body, 0);
    }
    glColor3fv (color);
    selection_render_body (picked_body);
    glEnable (GL_LIGHTING);
  }
}

/* render picked point if any */
static void render_picked_point (void)
{
  GLfloat color [3] = {0.5, 0.5, 0.5};

  if (picked_point)
  {
    switch (tool_mode)
    {
    case TOOLS_TRACKBALL_CENTER:
    case TOOLS_POINTS_COORDS: color [0] = 1.0; break;
    case TOOLS_POINTS_DISTANCE: color [2] = 1.0; break;
    case TOOLS_POINTS_ANGLE: color [1] = 1.0; break;
    }

    glDisable (GL_LIGHTING);
    glDisable (GL_DEPTH_TEST);
    glPointSize (8.0);
    glBegin (GL_POINTS);
    glColor3fv (color);
    glVertex3dv (picked_point);
    glEnd ();
    glPointSize (1.0);
    glEnable (GL_DEPTH_TEST);
    glEnable (GL_LIGHTING);
  }
}

/* render the sketch of the cut plane */
static void render_cut_sketch (void)
{
  if (!cut_sketch) return;

  GLfloat color [4] = {1.0, 0.75, 0.75, 0.4};
  double v [3], w [3], u [3], q [4][3], d;
  int i;


  MAXABS (cut_normal, d);
  i = 0; v[0] = fabs (cut_normal [0]);
  if (fabs (cut_normal [1]) < v[0]) i = 1, v[0] = cut_normal [1];
  if (fabs (cut_normal [2]) < v[0]) i = 2, v[0] = cut_normal [2];
  COPY (cut_normal, w);
  w [i] += d;
  PRODUCT (cut_normal, w, v);
  PRODUCT (cut_normal, v, w);
  SUB (global_extents+3, global_extents, u);
  MAXABS (u, d);
  d *= 10.0;
  ADDMUL (cut_point, -d, v, q[0]);
  ADDMUL (cut_point, -d, w, q[1]);
  ADDMUL (cut_point,  d, v, q[2]);
  ADDMUL (cut_point,  d, w, q[3]);

  glDisable (GL_LIGHTING);
  glDisable (GL_CULL_FACE);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glNormal3dv (cut_normal);
  glColor4fv (color);
  glBegin (GL_QUADS);
  glVertex3dv (q[0]);
  glVertex3dv (q[1]);
  glVertex3dv (q[2]);
  glVertex3dv (q[3]);
  glEnd ();
  glDisable (GL_BLEND);
  glEnable (GL_CULL_FACE);
  glEnable (GL_LIGHTING);
}

/* render cuts */
static void render_cuts (void)
{
  GLfloat color [3] = {0, 0.75, 0};
  double *cur, *val;
  CUT_DATA *cut;
  TRI *t, *e;
  int i, j;

  glDisable (GL_LIGHTING);
  glDisable (GL_CULL_FACE);
  for (cut = cuts; cut; cut = cut->next)
  {
    if (!SET_Contains (selection->set, cut->bod, NULL)) continue;
    BODY_DATA *data = cut->bod->rendering;
    if (data->flags & HIDDEN) continue;

    cur = cut->cur;
    val = cut->val;

    for (t = cut->tri, e = t + cut->m; t != e; t ++)
    {
      glBegin (GL_TRIANGLES);
      glNormal3dv (t->out);
      if (cut->values_updated)
      {
	for (i = 0; i < 3; i ++)
	{
	  j = (t->ver [i] - cur)/3;
	  value_to_color (val [j], color);
	  glColor3fv (color);
	  glVertex3dv (t->ver[i]);
	}
      }
      else
      {
	COPY (neutral_color, color);
	glColor3fv (color);
	for (i = 0; i < 3; i ++) glVertex3dv (t->ver[i]);
      }
      glEnd ();
    }
  }
  glEnable (GL_CULL_FACE);
  glEnable (GL_LIGHTING);
}

/* update scene extents */
static void update_extents ()
{
  double e [6], *extents, margin;
  SET *item;

  extents = global_extents;

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

  SUB (extents+3, extents, e);
  margin = 1.0 / exp (MAX (e[0], MAX (e[1], e[2])) /  MIN (e[0], MIN (e[1], e[2]))); /* 0.36 for aspect ratio 1.0 extents */
  SUBMUL (extents, margin, e, extents);
  ADDMUL (extents+3, margin, e, extents+3);

  GLV_Update_Extents (extents);
}

/* update bodies */
static void update ()
{
  update_cuts ();

  for (BODY *bod = domain->bod; bod; bod = bod->next)
  {
    if (bod->rendering == NULL) bod->rendering = create_body_data (bod, WIREFRAME_FLAG ());
  }

  if (legend.entity)
  {
    legend.extents [0] =  DBL_MAX;
    legend.extents [1] = -DBL_MAX;

    SET_Free (&rndsetmem, &legend.discrete);

    if (legend_constraint_based ()) 
    {
      for (CON *con = domain->con; con; con = con->next) con->state &= ~CON_DONE; /* all undone */
    }

    for (SET *item = SET_First (selection->set); item; item = SET_Next (item))
    {
      BODY *bod = item->data;
      update_body_values (bod, bod->rendering);
    }

    update_cuts_values ();

    legend_enable ();
  }

  if (render_bodies)
  {
    for (BODY *bod = domain->bod; bod; bod = bod->next) update_body_data (bod, bod->rendering);
  }

  GLV_Resize_Viewport (time_window, time_width (), TIME_HEIGHT); /* stretch time window to fit text */

  GLV_Redraw_All (); /* redraw all widgets */
}

/* XXX => this may be vournable <= XXX */
#define SELECTION_REINIT_BEGIN()\
  double idsum0 = domain->nbod;\
  for (BODY *bod = domain->bod; bod; bod = bod->next) idsum0 += bod->id\

/* XXX => this may be vournable <= XXX */
#define SELECTION_REINIT_END()\
  double idsum1 = domain->nbod;\
  for (BODY *bod = domain->bod; bod; bod = bod->next) idsum1 += bod->id;\
  if (idsum1 != idsum0) selection_init ()

/* one simulation step */
static void step ()
{
  SOLVER_DATA *s = MAP_Find (solvers, domain, NULL);

  double epsilon = DBL_EPSILON;

  SELECTION_REINIT_BEGIN ();

  /* find a small number such that added to current time it makes a difference */
  while (domain->time + epsilon == domain->time) epsilon += DBL_EPSILON;

  /* (***) note that domain->step might be decreased due to stability issues */
  if (s) SOLFEC_Run (solfec, s->kind, s->solver, epsilon); /* use epsilon as the duration in order to make just one step; see (***) */
  else 
  {
    PENALTY *ps = PENALTY_Create (1);
    SOLFEC_Run (solfec, PENALTY_SOLVER, ps, epsilon); /* default and in read mode */
    PENALTY_Destroy (ps);
  }

  SELECTION_REINIT_END ();

  update ();
}

static int forward ()
{
  SELECTION_REINIT_BEGIN ();

  int ret = SOLFEC_Forward (solfec, skip_steps);

  SELECTION_REINIT_END ();

  return ret;
}

static int backward ()
{
  SELECTION_REINIT_BEGIN ();

  int ret = SOLFEC_Backward (solfec, skip_steps);

  SELECTION_REINIT_END ();

  return ret;
}

/* run simulation */
static void run (int dummy)
{
  if (solfec->mode == SOLFEC_WRITE) step ();
  else 
  {
    if (forward ()) update ();
    else GLV_AVI_Stop ();
  }

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

    SELECTION_REINIT_BEGIN ();

    SOLFEC_Seek_To (solfec, t);

    SELECTION_REINIT_END ();

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

/* read cut normal */
static void read_cut_normal (char *text)
{
  if (text)
  {
    SET (cut_normal, 0.0)
    sscanf (text, "%lf%lf%lf", cut_normal, cut_normal+1, cut_normal+2);
    if (LEN (cut_normal) > 0.0)
    {
      NORMALIZE (cut_normal);
      MID (global_extents, global_extents + 3, cut_point);
      cut_sketch = 1;
      mouse_mode = MOUSE_CUTTING_PLANE;
      tip = "Press 'c' in order create the cut";
      GLV_Redraw_All ();
    }
  }
}

/* 2D selection using indexed coloring */
static void select_2D (int x1, int y1, int x2, int y2)
{
  unsigned char (*pix) [4];
  SET *ids, *set, *item;
  int x, y, w, h, n, m;
  GLint viewport [4];
  BODY *bod;

  glDisable (GL_LIGHTING);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  selection_2D_render_body_set (selection->set);
  glEnable (GL_LIGHTING);

  glGetIntegerv (GL_VIEWPORT, viewport);

  x = MIN (x1, x2);
  y = MIN (viewport[3]-y1, viewport[3]-y2);
  w = MAX (x1, x2) - x;
  h = MAX (viewport[3]-y1, viewport[3]-y2) - y;
  w = MAX (w, 1);
  h = MAX (h, 1);
  m = w * h;

  ERRMEM (pix = malloc (m * sizeof (unsigned char [4])));
  glReadPixels (x, y, w, h, GL_RGBA, GL_UNSIGNED_BYTE, pix);
  for (ids = NULL, n = 0; n < m; n ++) SET_Insert (&rndsetmem, &ids, (void*) (long) rgbatoid (pix [n]), NULL);
  free (pix);

  if (ids)
  {
    for (set = NULL, item = SET_First (ids); item; item = SET_Next (item))
    {
      bod = MAP_Find (domain->idb, item->data, NULL);
      if (bod) SET_Insert (&rndsetmem, &set, bod, NULL);
    }

    if (set) selection_push (set);

    SET_Free (&rndsetmem, &ids);
  }
}

/* 3D selection using indexed coloring */
static void select_3D (int x1, int y1, int x2, int y2)
{
  GLuint *sel, selsize, selcount, n, m, k;
  GLint viewport [4];
  int x, y, w, h; 
  BODY *bod;
  SET *set;

  w = ABS (x1 - x2);
  h = ABS (y1 - y2);
  x = (x1 + x2) / 2;
  y = (y1 + y2) / 2;
  w = MAX (w, 2);
  h = MAX (h, 2);
  selsize = SET_Size (selection->set) * 4; /* assumes that each body id stored in a separate hit, which is excessive and must be enough */
  ERRMEM (sel = malloc (sizeof (GLuint [selsize])));

  glSelectBuffer (selsize, sel);
  glRenderMode (GL_SELECT);
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glLoadIdentity ();
  glGetIntegerv (GL_VIEWPORT, viewport);
  gluPickMatrix (x, viewport [3] - y, w, h, viewport);
  GLV_SetProjectionMatrix (viewport [2], viewport [3]);
  glMatrixMode (GL_MODELVIEW);
  glInitNames ();
  selection_3D_render_body_set (selection->set);
  glMatrixMode (GL_PROJECTION);
  glPopMatrix ();
  glMatrixMode (GL_MODELVIEW);
  glFlush ();
  selcount = glRenderMode (GL_RENDER);

  if (selcount > 0 && selcount < INT_MAX)
  {
    for (n = m = 0, set = NULL; n < selcount; n ++, m += sel [m]+3)
    {
      for (k = 0; k < sel [m]; k ++)
      {
	bod = MAP_Find (domain->idb, (void*) sel [m+3+k], NULL);
	if (bod) SET_Insert (&rndsetmem, &set, bod, NULL);
      }
    }

    if (set) selection_push (set);
  }

  free (sel);
}

/* pick one body using 2D selection */
static BODY* pick_body (int x, int y)
{
  unsigned char pix [4];
  GLint viewport [4];

  glDisable (GL_LIGHTING);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  selection_2D_render_body_set (selection->set);
  glEnable (GL_LIGHTING);

  glGetIntegerv (GL_VIEWPORT, viewport);
  glReadPixels (x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pix);

  return MAP_Find (domain->idb, (void*) rgbatoid (pix), NULL);
}

/* pick one point using 2D selection */
static double* pick_point (int x, int y)
{
  unsigned char pix [4];
  GLint viewport [4];
  double *point;
  BODY *bod;

  glDisable (GL_LIGHTING);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  selection_2D_render_body_set (selection->set);
  glEnable (GL_LIGHTING);

  glGetIntegerv (GL_VIEWPORT, viewport);
  glReadPixels (x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pix);

  bod = MAP_Find (domain->idb, (void*) rgbatoid (pix), NULL);
  point = NULL;

  if (bod)
  {
    glDisable (GL_LIGHTING);
    glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    selection_2D_render_point_set (bod->rendering);
    glEnable (GL_LIGHTING);

    glReadPixels (x, viewport[3] - y, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, pix);
    BODY_DATA *data = bod->rendering;
    int i = rgbatoid (pix);
    if (i > 0 && i <= data->values_count) point = data->value_sources [i-1].pnt;
  }

  return point;
}

/* switch wireframe modes */
static void switchwireframe (SET *bodies)
{
  BODY_DATA *data;
  SET *item;
  BODY *bod;

  for (item = SET_First (bodies); item; item = SET_Next (item))
  {
    bod = item->data;
    data = bod->rendering;
    if (data)
    {
      int wire = data->flags & WIREFRAME;
      RND_Free_Rendering_Data (data);
      if (wire) bod->rendering = create_body_data (bod, 0);
      else bod->rendering = create_body_data (bod, 1);
    }
  }
}

/* disable all modes */
static int modes_off ()
{
  int ret = tool_mode;
  mouse_mode = MOUSE_NONE;
  selection_mode = SELECTION_NONE;
  tool_mode = 0;
  tip = NULL;
  cut_sketch = 0;
  picked_body = NULL;
  picked_point = NULL;
  picked_point_hist [0] = NULL;
  picked_point_hist [1] = NULL;
  GLV_Rectangle_Off ();
  GLV_Release_Mouse ();
  update ();
  return ret;
}

/* time window width */
static int time_width ()
{
  if (show_outpath) return GLV_Print_Width (TIME_FONT, "t=%g [%s]", domain->time, solfec->outpath) + 6;
  else if (tip) return GLV_Print_Width (TIME_FONT, "t=%g (%s)", domain->time, tip) + 6;
  else return GLV_Print_Width (TIME_FONT, "t=%g", domain->time) + 6;
}

/* render current time */
static void time_render ()
{
  glDisable (GL_LIGHTING);
  glDisable (GL_DEPTH_TEST);

  glColor3f (1, 1, 1);
  glRecti (0, 0, time_width (), TIME_HEIGHT);
  glColor3f (0, 0, 0);
  if (show_outpath) GLV_Print (3, 3, 1, TIME_FONT, "t=%g [%s]", domain->time, solfec->outpath);
  else if (tip) GLV_Print (3, 3, 1, TIME_FONT, "t=%g (%s)", domain->time, tip);
  else GLV_Print (3, 3, 1, TIME_FONT, "t=%g", domain->time);

  glEnable (GL_LIGHTING);
  glEnable (GL_DEPTH_TEST);
}

/* render coordinates */
static void coord_render ()
{
  glDisable (GL_LIGHTING);

  glColor3d (0., 0., 0.);
  glBegin (GL_LINES);
  glVertex3d (0., 0., 0.);
  glVertex3d (.7, 0., 0.);
  glVertex3d (0., 0., 0.);
  glVertex3d (0., .7, 0.);
  glVertex3d (0., 0., 0.);
  glVertex3d (0., 0., .7);
  glEnd ();

  GLV_Print (.7, 0, .2, TIME_FONT, "x");
  GLV_Print (0, .7, .2, TIME_FONT, "y");
  GLV_Print (0,  0, .8, TIME_FONT, "z");

  glEnable (GL_LIGHTING);
}

/* menus */
static void menu_tools (int item);

static void menu_domain (int item)
{
  int prevmode = solfec->mode;
  DOM *prevdom = domain;
  
  switch (item)
  {
  case DOMAIN_NEXT:
    if (domain->next) domain = domain->next;
    break;
  case DOMAIN_PREVIOUS:
    if (domain->prev) domain = domain->prev;
    break;
  }

  /* update analysis menu */

  if (domain != prevdom)
  {
    menu_tools (TOOLS_CLEAR_CUTS);

    selection_init ();

    update_extents ();

    update ();

    glutSetMenu (menu_code [MENU_ANALYSIS]);

    if (domain->flags & DOM_RUN_ANALYSIS)
      glutChangeToMenuEntry (1, "stop /RETURN/", ANALYSIS_STOP);
    else glutChangeToMenuEntry (1, "run /RETURN/", ANALYSIS_RUN);

    if (solfec->mode == SOLFEC_READ && prevmode == SOLFEC_WRITE)
    {
      glutAddMenuEntry ("seek to /UP/", ANALYSIS_SEEKTO);
      glutAddMenuEntry ("forward /LEFT/", ANALYSIS_FORWARD);
      glutAddMenuEntry ("backward /RIGHT/", ANALYSIS_BACKWARD);
      glutAddMenuEntry ("skip /DOWN/", ANALYSIS_SKIP);
    }
    else if (solfec->mode == SOLFEC_WRITE && prevmode == SOLFEC_READ)
    {
      glutRemoveMenuItem (6);
      glutRemoveMenuItem (5);
      glutRemoveMenuItem (4);
      glutRemoveMenuItem (3);
    }
  }
}

static void menu_render (int item)
{
  switch (item)
  {
  case RENDER_SELECTION_2D:
  case RENDER_WIREFRAME_2D:
    legend_disable ();
    modes_off ();
    mouse_mode = MOUSE_SELECTION_BEGIN;
    selection_mode = SELECTION_2D;
    if (item == RENDER_WIREFRAME_2D) wireframeon = 1;
    else wireframeon = 0;
    GLV_Hold_Mouse ();
    break;
  case RENDER_SELECTION_3D:
  case RENDER_WIREFRAME_3D:
    legend_disable ();
    modes_off ();
    mouse_mode = MOUSE_SELECTION_BEGIN;
    selection_mode = SELECTION_3D;
    if (item == RENDER_WIREFRAME_3D) wireframeon = 1;
    else wireframeon = 0;
    GLV_Hold_Mouse ();
    break;
  case RENDER_PREVIOUS_SELECTION:
    selection_pop ();
    update_extents ();
    update ();
    break;
  case RENDER_BODIES:
    render_bodies = !render_bodies;
    update ();
    break;
  }
}

static void menu_tools (int item)
{
  switch (item)
  {
  case TOOLS_EULER_CUT:
  case TOOLS_LAGRANGE_CUT:
    if (mouse_mode == MOUSE_CUTTING_PLANE)
    {
      for (SET *item = SET_First (selection->set); item; item = SET_Next (item))
      {
	BODY *bod = item->data;
	BODY_DATA *data = bod->rendering;
	double *ref, *cur;
	CUT_DATA *cut;
	SGP *sgp;
	TRI *tri;
	int m, n;

	if (bod->kind == OBS) continue;
	if (data->flags & HIDDEN) continue;

	tri = SHAPE_Cut (bod->shape, cut_point, cut_normal, &m, bod, (MOTION)BODY_Ref_Point, &sgp, &ref, &cur, &n);

	if (tri)
	{
	  data->flags |= SEETHROUGH; /* make the body transparent if it was cut */
	  ERRMEM (cut = MEM_CALLOC (sizeof (CUT_DATA) + sizeof (double [n])));
	  cut->val = (double*) (cut+1);
	  cut->kind = cut_kind;
	  COPY (cut_point, cut->point);
	  COPY (cut_normal, cut->normal);
	  cut->bod = bod;
	  cut->tri = tri;
	  cut->sgp = sgp;
	  cut->ref = ref;
          cut->cur = cur;
	  cut->m = m;
	  cut->n = n;
	  cut->next = cuts;
	  cuts = cut;
	  update_cut_values (cut);
	}
      }

      if (cut_kind == EULER)
      {
	double *point_normal;
	ERRMEM (point_normal = malloc (sizeof (double [6])));
	COPY (cut_point, point_normal);
	COPY (cut_normal, point_normal + 3);
	SET_Insert (NULL, &euler_cuts, point_normal, NULL); /* record Euler cut */
      }

      tip = NULL;
      cut_sketch = 0;
      mouse_mode = MOUSE_NONE;
      GLV_Redraw_All ();
    }
    else 
    {
      if (item == TOOLS_EULER_CUT)
      {
	if (!GLV_Reading_Text ())
	{
	  cut_kind = EULER;
	  GLV_Read_Text ("Euler cut normal (format: nx ny nz)", read_cut_normal);
	}
      }
      else
      {
	if (!GLV_Reading_Text ())
	{
	  cut_kind = LAGRANGE;
	  GLV_Read_Text ("Lagrange cut normal (format: nx ny nz)", read_cut_normal);
	}
      }
    }
    break;
  case TOOLS_CLEAR_CUTS:
    {
      CUT_DATA *cut, *next;
      SET *item;

      for (cut = cuts; cut; cut = next)
      {
        BODY_DATA *data = cut->bod->rendering;	
	data->flags &= ~SEETHROUGH;
	next = cut->next;
	free (cut->tri);
	free (cut);
      }

      for (item = SET_First (euler_cuts); item; item = SET_Next (item))
      {
	free (item->data);
      }

      SET_Free (NULL, &euler_cuts);
      cuts = NULL; /* clear cuts list */
      update ();
    }
    break;
  case TOOLS_TRANSPARENT:
  case TOOLS_ROUGH_MESH:
  case TOOLS_HIDE:
    modes_off ();
    mouse_mode = MOUSE_PICK_BODY;
    tool_mode = item;
    GLV_Hold_Mouse ();
    switch (item)
    {
      case TOOLS_TRANSPARENT: GLV_Window_Title ("Solfec: transparent"); break;
      case TOOLS_ROUGH_MESH: GLV_Window_Title ("Solfec: rough mesh"); break;
      case TOOLS_HIDE: GLV_Window_Title ("Solfec: hide"); break;
    }
    break;
  case TOOLS_TRANSPARENT_ALL:
    for (BODY *bod = domain->bod; bod; bod = bod->next)
    { BODY_DATA *data = bod->rendering; data->flags |= SEETHROUGH; }
    GLV_Redraw_All ();
    break;
  case TOOLS_TRANSPARENT_NONE:
    for (BODY *bod = domain->bod; bod; bod = bod->next)
    { BODY_DATA *data = bod->rendering; data->flags &= ~SEETHROUGH; }
    for (CUT_DATA *cut = cuts; cut; cut = cut->next)
    { BODY_DATA *data = cut->bod->rendering;	data->flags |= SEETHROUGH; }
    GLV_Redraw_All ();
    break;
  case TOOLS_SHOW_ALL:
    for (BODY *bod = domain->bod; bod; bod = bod->next)
    { BODY_DATA *data = bod->rendering; data->flags &= ~HIDDEN; }
    update_extents ();
    break;
  case TOOLS_NEXT_RESULT:
    if (legend.entity)
    {
      if (legend.entity < RESULTS_MERIT)
      {
	legend.entity ++;
	update ();
      }
    }
    else menu_kinds (KINDS_OF_CONSTRAINTS);
    break;
  case TOOLS_PREVIOUS_RESULT:
    if (legend.entity)
    {
      if (legend.entity > KINDS_OF_CONSTRAINTS)
      {
	legend.entity --;
	update ();
      }
    }
    else menu_results (RESULTS_MERIT);
    break;
  case TOOLS_SMALLER_SCALING:
    if (modal_analysis_menu)
    {
      if (current_eigenmode >= 0)
      {
	if (eigenshape_factor > 0.02) eigenshape_factor -= 0.01;
	menu_modal_analysis (current_eigenmode);
      }
    }
    else
    {
      if (arrow_factor > 0.02) arrow_factor -= 0.01;
    }

    GLV_Redraw_All ();
    break;
  case TOOLS_BIGGER_SCALING:
    if (modal_analysis_menu)
    {
      if (current_eigenmode >= 0)
      {
	if (eigenshape_factor < 0.5) eigenshape_factor += 0.01;
	menu_modal_analysis (current_eigenmode);
      }
    }
    else
    {
      if (arrow_factor < 0.10) arrow_factor += 0.01;
    }
    GLV_Redraw_All ();
    break;
  case TOOLS_OUTPATH:
    show_outpath = !show_outpath;
    GLV_Redraw_All ();
    break;
  case TOOLS_POINTS_COORDS:
  case TOOLS_POINTS_DISTANCE:
  case TOOLS_POINTS_ANGLE:
  case TOOLS_TRACKBALL_CENTER:
    modes_off ();
    mouse_mode = MOUSE_PICK_POINT;
    tool_mode = item;
    GLV_Hold_Mouse ();
    switch (item)
    {
      case TOOLS_POINTS_COORDS: GLV_Window_Title ("Solfec: points coords"); break;
      case TOOLS_POINTS_DISTANCE: GLV_Window_Title ("Solfec: points distance"); break;
      case TOOLS_POINTS_ANGLE: GLV_Window_Title ("Solfec: points angle"); break;
      case TOOLS_TRACKBALL_CENTER: GLV_Window_Title ("Solfec: trackball center"); break;
    }
    break;
  case TOOLS_WIREFRAME_ALL:
    for (SET *item = SET_First (selection->set); item; item = SET_Next (item))
    {
      BODY *bod = item->data;
      BODY_DATA *data = bod->rendering;
      if (data)
      {
	if (!(data->flags & WIREFRAME))
	{
	  RND_Free_Rendering_Data (data);
	  bod->rendering = create_body_data (bod, 1);
	}
      }
    }
    update ();
    break;
  case TOOLS_WIREFRAME_NONE:
    for (SET *item = SET_First (selection->set); item; item = SET_Next (item))
    {
      BODY *bod = item->data;
      BODY_DATA *data = bod->rendering;
      if (data)
      {
	if (data->flags & WIREFRAME)
	{
	  RND_Free_Rendering_Data (data);
	  bod->rendering = create_body_data (bod, 0);
	}
      }
    }
    update ();
    break;
  case TOOLS_NEXT_MODE:
    if (modal_analysis_menu)
    {
      BODY *bod = selection->set->data;
      if (current_eigenmode < (bod->evec->n-1)) current_eigenmode ++;
      else current_eigenmode = 0;
      menu_modal_analysis (current_eigenmode);
    }
    GLV_Redraw_All ();
    break;
  case TOOLS_PREVIOUS_MODE:
    if (modal_analysis_menu)
    {
      BODY *bod = selection->set->data;
      if (current_eigenmode > 0) current_eigenmode --;
      else current_eigenmode = bod->evec->n-1;
      menu_modal_analysis (current_eigenmode);
    }
    GLV_Redraw_All ();
    break;
  }
}

static void menu_kinds (int item)
{
  legend.entity = item;
  update ();
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
      if (!GLV_Reading_Text ())
      {
	static char caption [256];
	double s, e;

	SOLFEC_Time_Limits (solfec, &s, &e);
	snprintf (caption, 256, "Seek to time from [%g, %g]", s, e);
	GLV_Read_Text (caption, seek_to_time);
      }
    }
    break;
  case ANALYSIS_FORWARD:
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    forward ();
    update ();
    break;
  case ANALYSIS_BACKWARD:
    if (domain->flags & DOM_RUN_ANALYSIS) menu_analysis (ANALYSIS_STOP);
    backward ();
    update ();
    break;
  case ANALYSIS_SKIP:
    if (!GLV_Reading_Text ())
    {
      GLV_Read_Text ("FORWARD and BACKWARD skip", set_skip_steps);
    }
    break;
  }
}

static void menu_results (int item)
{
  legend.entity = item;
  update ();
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
  glutAddMenuEntry ("2D wireframe /4/", RENDER_WIREFRAME_2D);
  glutAddMenuEntry ("3D wireframe /5/", RENDER_WIREFRAME_3D);
  glutAddMenuEntry ("previous selection /p/", RENDER_PREVIOUS_SELECTION);
  glutAddMenuEntry ("toggle bodies /b/", RENDER_BODIES);

  menu_name [MENU_TOOLS] = "tools";
  menu_code [MENU_TOOLS] = glutCreateMenu (menu_tools);
  glutAddMenuEntry ("toggle transparent /t/", TOOLS_TRANSPARENT);
  glutAddMenuEntry ("transparent all /T/", TOOLS_TRANSPARENT_ALL);
  glutAddMenuEntry ("transparent none /n/", TOOLS_TRANSPARENT_NONE);
  glutAddMenuEntry ("Euler cut /c/", TOOLS_EULER_CUT);
  glutAddMenuEntry ("Lagrange cut /C/", TOOLS_LAGRANGE_CUT);
  glutAddMenuEntry ("clear cuts /s/", TOOLS_CLEAR_CUTS);
  glutAddMenuEntry ("hide /h/", TOOLS_HIDE);
  glutAddMenuEntry ("show all /a/", TOOLS_SHOW_ALL);
  glutAddMenuEntry ("toggle rough mesh /r/", TOOLS_ROUGH_MESH);
  glutAddMenuEntry ("next result /+/", TOOLS_NEXT_RESULT);
  glutAddMenuEntry ("previous result /-/", TOOLS_PREVIOUS_RESULT);
  glutAddMenuEntry ("bigger arrows /]/", TOOLS_BIGGER_SCALING);
  glutAddMenuEntry ("smaller arrows /[/", TOOLS_SMALLER_SCALING);
  glutAddMenuEntry ("toggle output path /o/", TOOLS_OUTPATH);
  glutAddMenuEntry ("point coordiantes /x/", TOOLS_POINTS_COORDS);
  glutAddMenuEntry ("points distance /d/", TOOLS_POINTS_DISTANCE);
  glutAddMenuEntry ("points angle /g/", TOOLS_POINTS_ANGLE);
  glutAddMenuEntry ("trackball center /L/", TOOLS_TRACKBALL_CENTER);
  glutAddMenuEntry ("wireframe all /w/", TOOLS_WIREFRAME_ALL);
  glutAddMenuEntry ("wireframe none /W/", TOOLS_WIREFRAME_NONE);

  menu_name [MENU_KINDS] = "kinds of";
  menu_code [MENU_KINDS] = glutCreateMenu (menu_kinds);
  glutAddMenuEntry ("constraints", KINDS_OF_CONSTRAINTS);
  glutAddMenuEntry ("constraint ranks", KINDS_OF_CONSTRAINT_RANKS);
  glutAddMenuEntry ("forces", KINDS_OF_FORCES);
  glutAddMenuEntry ("bodies", KINDS_OF_BODIES);
  glutAddMenuEntry ("body ranks", KINDS_OF_BODY_RANKS);
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
  glutAddMenuEntry ("RT", RESULTS_RT);
  glutAddMenuEntry ("RN", RESULTS_RN);
  glutAddMenuEntry ("R", RESULTS_R);
  glutAddMenuEntry ("UT", RESULTS_UT);
  glutAddMenuEntry ("UN", RESULTS_UN);
  glutAddMenuEntry ("U", RESULTS_U);
  glutAddMenuEntry ("gap", RESULTS_GAP);
  glutAddMenuEntry ("merit", RESULTS_MERIT);

  menu_name [MENU_RESULTS] = "results";
  menu_code [MENU_RESULTS] = glutCreateMenu (menu_results);
  glutAddSubMenu ("displacements", local [0]);
  glutAddSubMenu ("velocities", local [1]);
  glutAddSubMenu ("stresses", local [2]);
  glutAddSubMenu ("constraints", local [3]);

  *names = menu_name;
  *codes = menu_code;
  return MENU_LAST;
}

void RND_Init ()
{
  int w, h;

  ASSERT (domain, ERR_RND_NO_DOMAIN);

  MEM_Init (&rndsetmem, sizeof (SET), BIGCHUNK);

  MEM_Init (&rndmapmem, sizeof (MAP), BIGCHUNK);

  selection_init ();

  update_extents ();

  legend.preventity = 0;
  legend.entity = 0;
  legend.window = 0;
  legend.discrete = NULL;

  GLV_Sizes (&w, &h);

  time_window = GLV_Open_Viewport (0, -(h - TIME_HEIGHT), time_width (), TIME_HEIGHT, 0, time_render);

  coord_window = GLV_Open_Viewport (-(w - COORD_WIDTH), 0, COORD_WIDTH, COORD_HEIGHT, 1, coord_render);
}

int  RND_Idle ()
{
  return 0;
}

void RND_Quit ()
{
  DOM *next;

  /* clear cuts memory */
  menu_tools (TOOLS_CLEAR_CUTS);

  for (MAP *item = MAP_First (solvers); item; item = MAP_Next (item)) free (item->data); /* free solver interfaces */
  MAP_Free (NULL, &solvers);

  /* delete all SOLFEC objects */
  while (domain->prev) domain = domain->prev; /* rewind back */
  for (; domain; domain = next)
  {
    next = domain->next;
    SOLFEC_Destroy (solfec); /* see macro at the top */
  }

  lngfinalize (); /* finalize Python */
}

void RND_Render ()
{
  render_cuts ();

  render_body_set (selection->set);

  render_picked_body ();

  render_picked_point ();

  render_cut_sketch ();
}

void RND_Key (int key, int x, int y)
{
  switch (key)
  {
  case 27:
    legend_disable ();
    if (modes_off ()) legend_pop ();
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
  case '4':
    menu_render (RENDER_WIREFRAME_2D);
    break;
  case '5':
    menu_render (RENDER_WIREFRAME_3D);
    break;
  case 'p':
    menu_render (RENDER_PREVIOUS_SELECTION);
    break;
  case 'b':
    menu_render (RENDER_BODIES);
    break;
  case 't':
    menu_tools (TOOLS_TRANSPARENT);
    break;
  case 'T':
    menu_tools (TOOLS_TRANSPARENT_ALL);
    break;
  case 'n':
    menu_tools (TOOLS_TRANSPARENT_NONE);
    break;
  case 'c':
    menu_tools (TOOLS_EULER_CUT);
    break;
  case 'C':
    menu_tools (TOOLS_LAGRANGE_CUT);
    break;
  case 's':
    menu_tools (TOOLS_CLEAR_CUTS);
    break;
  case 'h':
    menu_tools (TOOLS_HIDE);
    break;
  case 'a':
    menu_tools (TOOLS_SHOW_ALL);
    break;
  case 'r':
    menu_tools (TOOLS_ROUGH_MESH);
    break;
  case '+':
    menu_tools (TOOLS_NEXT_RESULT);
    break;
  case '-':
    menu_tools (TOOLS_PREVIOUS_RESULT);
    break;
  case '[':
    menu_tools (TOOLS_SMALLER_SCALING);
    break;
  case ']':
    menu_tools (TOOLS_BIGGER_SCALING);
    break;
  case '{':
    menu_tools (TOOLS_PREVIOUS_MODE);
    break;
  case '}':
    menu_tools (TOOLS_NEXT_MODE);
    break;
  case 'o':
    menu_tools (TOOLS_OUTPATH);
    break;
  case 'x':
    menu_tools (TOOLS_POINTS_COORDS);
    break;
  case 'd':
    menu_tools (TOOLS_POINTS_DISTANCE);
    break;
  case 'g':
    menu_tools (TOOLS_POINTS_ANGLE);
    break;
  case 'L':
    menu_tools (TOOLS_TRACKBALL_CENTER);
    break;
  case 'w':
    menu_tools (TOOLS_WIREFRAME_ALL);
    break;
  case 'W':
    menu_tools (TOOLS_WIREFRAME_NONE);
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
  switch (button)
  {
    case GLUT_LEFT_BUTTON:
    switch (state)
    {
    case GLUT_UP:
      if (mouse_mode == MOUSE_SELECTION_END)
      {
	switch (selection_mode)
	{
	case SELECTION_2D: 
	  select_2D (mouse_start [0], mouse_start [1], x, y); 
	  if (wireframeon && selection->prev)
	  {
	    switchwireframe (selection->set);
            selection_pop ();
	  }
	  break;
	case SELECTION_3D:
	  select_3D (mouse_start [0], mouse_start [1], x, y); 
	  if (wireframeon && selection->prev)
	  {
	    switchwireframe (selection->set);
            selection_pop ();
	  }
	  break;
	default: break;
	}

	modes_off ();
        update_extents ();
	legend_pop ();
      }
      break;
    case GLUT_DOWN:
      if (mouse_mode == MOUSE_SELECTION_BEGIN)
      {
	mouse_mode = MOUSE_SELECTION_END;
	mouse_start [0] = x;
	mouse_start [1] = y;
      }
      else if (mouse_mode == MOUSE_PICK_BODY && picked_body)
      {
	BODY_DATA *data = picked_body->rendering;

	switch (tool_mode)
	{
	case TOOLS_TRANSPARENT:
	  if (data->flags & SEETHROUGH) data->flags &= ~SEETHROUGH;
	  else data->flags |= SEETHROUGH;
	  break;
	case TOOLS_HIDE: data->flags |= HIDDEN; break;
	case TOOLS_ROUGH_MESH:
	  if (data->flags & ROUGH_MESH) data->flags &= ~ROUGH_MESH;
	  else if (picked_body->msh) data->flags |= ROUGH_MESH; /* only if it has rough mesh */
	  break;
	}

	GLV_Redraw_All ();
      }
      else if (mouse_mode == MOUSE_PICK_POINT && picked_point)
      {
	switch (tool_mode)
	{
	case TOOLS_POINTS_COORDS:
	  printf ("POINT: %g, %g, %g\n", picked_point [0], picked_point [1], picked_point [2]);
	  break;
	case TOOLS_POINTS_DISTANCE:
	  if (picked_point_hist [0])
	  {
	    double d [3], len;

	    SUB (picked_point_hist [0], picked_point, d);
	    len = LEN (d);

	    printf ("DISTANCE from [%g, %g, %g] to [%g, %g, %g]: %g\n",
	      picked_point_hist [0][0], picked_point_hist [0][1], picked_point_hist [0][2], 
	      picked_point [0], picked_point [1], picked_point [2], len);

	    picked_point_hist [0] = NULL;

	    GLV_Redraw_All ();
	  }
	  else picked_point_hist [0] = picked_point;
	  break;
	case TOOLS_POINTS_ANGLE:
	  if (picked_point_hist [0] && picked_point_hist [1])
	  {
	    double d0 [3], d1 [3], l0, l1, angle;

	    SUB (picked_point_hist [1], picked_point_hist [0], d0);
	    SUB (picked_point, picked_point_hist [0], d1);
	    l0 = LEN (d0);
	    l1 = LEN (d1);
	    if (l0 > 0.0 && l1 > 0.0)
	    {
	      angle = 180.0 * acos (DOT (d0, d1) / (l0*l1)) / ALG_PI;

	      printf ("ANGLE between [%g, %g, %g] and [%g, %g, %g] and [%g, %g, %g]: %g\n",
		picked_point_hist [0][0], picked_point_hist [0][1], picked_point_hist [0][2], 
		picked_point_hist [1][0], picked_point_hist [1][1], picked_point_hist [1][2], 
		picked_point [0], picked_point [1], picked_point [2], angle);
	    }
	    else printf ("Coincident points has been picked. Try again.\n");

	    picked_point_hist [0] = NULL;
	    picked_point_hist [1] = NULL;

	    GLV_Redraw_All ();
	  }
	  else 
	  {
	    picked_point_hist [1] = picked_point_hist [0];
	    picked_point_hist [0] = picked_point;
	  }
	  break;
	case TOOLS_TRACKBALL_CENTER:
	  GLV_Trackball_Center (picked_point);
          RND_Key (27, 0, 0); /* emulate ESC */
	  break;
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

void RND_Motion (int x, int y)
{
  if (mouse_mode == MOUSE_SELECTION_END)
  {
    GLV_Rectangle_On (mouse_start [0], mouse_start [1], x, y);
  }
}

void RND_Passive (int x, int y)
{
  if (mouse_mode == MOUSE_PICK_BODY)
  {
    picked_body = pick_body (x, y);

    if (picked_body)
    {
      if (picked_body->label)
	snprintf (picked_body_tip, PICKED_BODY_TIP_LEN,
	"BODY LABEL: %s, ID: %d", picked_body->label, picked_body->id);
      else snprintf (picked_body_tip, PICKED_BODY_TIP_LEN, "BODY ID: %d", picked_body->id);

      tip = picked_body_tip; /* tip in the upper left part of the window (after time) */
      GLV_Resize_Viewport (time_window, time_width (), TIME_HEIGHT); /* stretch time window to fit text */
    }

    GLV_Redraw_All ();
  }
  else if (mouse_mode == MOUSE_PICK_POINT)
  {
    picked_point = pick_point (x, y);

    GLV_Redraw_All ();
  }
  else if (mouse_mode == MOUSE_CUTTING_PLANE)
  {
    int w, h, ystart;
    double v [3], u, d;

    GLV_Sizes (&w, &h);
    ystart = h / 2;
    SUB (global_extents+3, global_extents, v);
    MAXABS (v, d);
    u = -(double)(y - ystart) / (double) MIN (w, h); 
    MID (global_extents+3, global_extents, v);
    ADDMUL (v, u*d, cut_normal, cut_point);
    cut_sketch = 1;
    GLV_Redraw_All ();
  }
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
  static DOM *tail = NULL;

  if (tail)
  {
    tail->next = dom;
    dom->prev = tail;
    tail = dom;
  }
  else domain = tail = dom;
}

void RND_Solver (DOM *dom, int kind, void *solver)
{
  SOLVER_DATA *data;
  MAP *item;

  if ((item = MAP_Find (solvers, dom, NULL))) data = item->data; /* already mapped */
  else
  {
    ERRMEM (data = malloc (sizeof (SOLVER_DATA)));
    MAP_Insert (NULL, &solvers, dom, data, NULL); /* map to domain */
  }

  data->solver = solver;
  data->kind = kind;
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
  free (data->ellips);
  free (data->ellip_colors);

#if VBO
  glDeleteBuffersARB (1, &data->triangles);
  glDeleteBuffersARB (1, &data->lines);
#else
  free (data->ver);
  free (data->lin);
#endif

  if (data->rough) RND_Free_Rendering_Data (data->rough);

  free (data);
}
