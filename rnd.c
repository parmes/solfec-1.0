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
#include "fem.h"
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
	  *vertex_values, /* as many values as vertices (triangles_count * 3) */
	 **vertex_value_sources, /* maps vertices to 'values' */
	 **line_sources;

  double *values; /* unique vertex scalar field values */

  VALUE_SOURCE *value_sources; /* value sources */

  int values_count;

  short values_updated; /* update flag */

  ELEPNT **color_sources; /* only for FEM */

  MAP *surfaces; /* map surface ids to sets of &vertex_values [i] pointers */

  MAP *volumes; /* map volume ids to sets of &vertex_values [i] pointers */

  SPHERE **spheres; /* sphere spheres */

  GLfloat *sphere_colors; /* sphere surface colors */

  int spheres_count;

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
  RENDER_PREVIOUS_SELECTION,
  TOOLS_TRANSPARENT,
  TOOLS_HIDE,
  TOOLS_SHOW_ALL,
  TOOLS_ROUGH_MESH,
  TOOLS_PREVIOUS_RESULT,
  TOOLS_NEXT_RESULT,
  TOOLS_SMALLER_ARROWS,
  TOOLS_BIGGER_ARROWS,
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

enum mouse_mode
{
  MOUSE_NONE,
  MOUSE_SELECTION_BEGIN,
  MOUSE_SELECTION_END,
  MOUSE_PICK_BODY
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

static MEM rndsetmem; /* sets memory */

static MEM rndmapmem;  /* maps memory */

static int skip_steps = 1; /* number of steps to skip when rewinding analysis */

static MOUSE_MODE mouse_mode = MOUSE_NONE;

static SELECTION_MODE selection_mode = SELECTION_NONE;

static GLfloat neutral_color [3] = {0.8, 0.8, 0.8};

static int mouse_start [2] = {0, 0};

static short tool_mode = 0; /* current tool */

static BODY *picked_body = NULL; /* currently picked body */

static double arrow_factor = 0.05; /* arrow drawing constant */

static int time_window = 0; /* time window handler */
#define TIME_HEIGHT 16 /* time window height */
#define TIME_FONT GLV_FONT_8_BY_13

static LEGEND_DATA legend; /* legend data */
#define LEGEND_ROWS 8 /* number of rows in the legend */
#define LEGEND_WIDTH_DISC 50 /* legend width for discrete data */
#define LEGEND_WIDTH_CONT 100 /* legend width for continuous data */
#define LEGEND_FONT GLV_FONT_8_BY_13

/* body transparency test */
#define TRANSPARENT(body) (((BODY_DATA*)((BODY*)(body))->rendering)->flags & TRANSPARENT)

/* body rough mesh flag test */
#define ROUGH_MESH(body) (((BODY_DATA*)((BODY*)(body))->rendering)->flags & ROUGH_MESH)

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
}

/* push new selection on stack */
static void selection_push (SET *set)
{
  SELECTION *s;

  ERRMEM (s = malloc (sizeof (SELECTION)));
  s->prev = selection;
  s->set = set;

  selection = s;
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
       legend.entity <= RESULTS_R)) return 1;

  return 0;
}

/* obtain scalar sphere point value */
static double sphere_value (BODY *bod, SPHERE *sph)
{
  double values [7];
  VALUE_KIND kind;
  short index;

  switch (legend.entity)
  {
  case KINDS_OF_BODIES: return bod->kind;
  case KINDS_OF_SURFACES:  return sph->surface;
  case KINDS_OF_VOLUMES: return sph->volume;
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

  BODY_Point_Values (bod, sph->ref_center, kind, values);

  return values [index];
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
static BODY_DATA* create_body_data (BODY *bod)
{
  double **vsr, **nsr, **lsr, **vvs, **end, *pla, *val;
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

  ERRMEM (data = MEM_CALLOC (sizeof (BODY_DATA)));

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

	ERRMEM (data->spheres = realloc (data->spheres, data->spheres_count * sizeof (SPHERE*)));
	j = (data->spheres_count - 1);
	data->spheres [j] = sph;
      }
      break;
    }
  }

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
      for (ele = msh->surfeles; ele; ele = ele->next)
      {
	for (fac = ele->faces; fac; fac = fac->next)
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
	    register_identifier (&data->volumes, ele->volume, &val [0]);
	    register_identifier (&data->volumes, ele->volume, &val [1]);
	    register_identifier (&data->volumes, ele->volume, &val [2]);
	  }
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

      if (con->state & CON_DONERND) continue;

      switch (legend.entity)
      {
      case KINDS_OF_CONSTRAINTS: 
        SET_Insert (&rndsetmem, &legend.discrete, (void*) (long) con->kind, NULL);
	value = con->kind;
	break;
      case RESULTS_RT: value = LEN2 (con->R); break;
      case RESULTS_RN: value = fabs (con->R[2]); break;
      case RESULTS_R: value = LEN (con->R); break;
      }

      if (value < legend.extents [0]) legend.extents [0] = value;
      if (value > legend.extents [1]) legend.extents [1] = value;

      con->state |= CON_DONERND;
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

  if (legend.entity >= KINDS_OF_BODIES && legend.entity < RESULTS_RT)
  {
    switch (legend.entity)
    {
    case KINDS_OF_BODIES:
      for (val = data->vertex_values, end = val + data->triangles_count * 3; val < end; val ++) *val = (double) bod->kind;
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
      else if (bod->kind != OBS)
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
	   end = val + data->triangles_count * 3; val < end; vvs ++, val ++) *val = **vvs;
    }

    SPHERE **sph, **end;

    for (sph = data->spheres, end = sph + data->spheres_count; sph < end; sph ++)
    {
      values [0] = sphere_value (bod, *sph);
      if (values[0] < legend.extents [0])  legend.extents [0] = values[0];
      if (values[0] > legend.extents [1])  legend.extents [1] = values[0];
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
       end = vsr + data->triangles_count * 3,
       v = ver, n = nor; vsr < end;
       vsr ++, nsr ++, val ++, v += 3, n += 3)
  {
    COPY (*vsr, v);
    COPY (*nsr, n);
  }

  if (data->values_updated)
  {
    for (val = data->vertex_values, tail = val + data->triangles_count * 3, c = col; val < tail;  val ++, c += 3) value_to_color (*val, c);
  }
  else for (c = col, l = c + data->triangles_count * 9; c < l;  c += 3) { COPY (neutral_color, c); }

  glUnmapBufferARB (GL_ARRAY_BUFFER_ARB);

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

  if (data->rough) update_body_data (bod, data->rough);
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

/* render sphere triangles for selection */
static void selection_render_sphere_triangles (double *center, double radius)
{
  glMatrixMode (GL_MODELVIEW_MATRIX);
  glPushMatrix ();
    glTranslated (center[0], center[1], center[2]);
    glutSolidSphere (radius, 12, 12);
  glPopMatrix ();
}

/* get legend caption */
static char* legend_caption ()
{
  switch (legend.entity)
  {
  case KINDS_OF_CONSTRAINTS: return "CONSTRAINT KINDS";
  case KINDS_OF_FORCES: return "FORCE KINDS";
  case KINDS_OF_BODIES: return "BODY KINDS";
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
  case RESULTS_RN: return "RN";
  case RESULTS_RT: return "RT";
  case RESULTS_R: return "R";
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
  GLV_Print (v[0] + 6, v[1] + 6, 0, LEGEND_FONT, "%s", legend_caption ());

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
        GLV_Print (v[0] + j * LEGEND_WIDTH_DISC + 18, v[1] + i * 16 + 3, 0, LEGEND_FONT, str);
      }
      else
      {
	glColor3f (1, 1, 1);
	l = GLV_Print_Width (LEGEND_FONT, "%d", (int)item->data) + 5;
	glRecti (v[0] + j * LEGEND_WIDTH_DISC + 16, v[1] + i * 16, v[0] + j * LEGEND_WIDTH_DISC + 16 + l, v[1] + (i+1) * 16);
        glColor3f (0, 0, 0);
	GLV_Print (v[0] + j * LEGEND_WIDTH_DISC + 18, v[1] + i * 16 + 3, 0, LEGEND_FONT, "%d", (int)item->data);
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
	GLV_Print (v[0] + j * LEGEND_WIDTH_CONT + 18, v[1] + i * 16 + 3, 0, LEGEND_FONT, "%.2e", value);
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
  if (legend.discrete)legend.range = SET_Size (legend.discrete);
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
static void render_body_triangles (BODY *bod, short skip)
{
  BODY_DATA *data;

  if (bod->rendering == NULL) bod->rendering = create_body_data (bod);

  data = bod->rendering;

  if (bod == picked_body ||           /* do not render a picked body */
      data->flags & skip) return;

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

  SPHERE **sph, **end;
  GLfloat *col;

  for (sph = data->spheres, end = sph + data->spheres_count, col = data->sphere_colors; sph < end; sph ++, col += 3)
    render_sphere_triangles ((*sph)->cur_center, (*sph)->cur_radius, col);
}

/* render body lines */
static void render_body_lines (BODY *bod, short skip)
{
  BODY_DATA *data;

  if (bod->rendering == NULL) bod->rendering = create_body_data (bod);

  data = bod->rendering;

  if (data->flags & skip) return;

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->lines);

  glEnableClientState (GL_VERTEX_ARRAY);

  glVertexPointer (3, GL_FLOAT, 0, 0);

  glDrawArrays (GL_LINES, 0, data->lines_count * 2);

  glDisableClientState (GL_VERTEX_ARRAY);

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);

  SPHERE **sph, **end;

  for (sph = data->spheres, end = sph + data->spheres_count; sph < end; sph ++)
    render_sphere_points ((*sph)->cur_points[0], (*sph)->cur_points[1], (*sph)->cur_points[2]);
}

/* render body triangles for selection */
static void selection_render_body_triangles (BODY *bod)
{
  BODY_DATA *data = bod->rendering;

  if (data->flags & HIDDEN) return;

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, data->triangles);

  glEnableClientState (GL_VERTEX_ARRAY);
  glEnableClientState (GL_NORMAL_ARRAY);

  glVertexPointer (3, GL_FLOAT, 0, 0);
  glNormalPointer (GL_FLOAT, 0, (void*) (data->triangles_count * sizeof (GLfloat) * 9));

  glDrawArrays (GL_TRIANGLES, 0, data->triangles_count * 3);

  glDisableClientState (GL_VERTEX_ARRAY);
  glDisableClientState (GL_NORMAL_ARRAY);

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);

  SPHERE **sph, **end;

  for (sph = data->spheres, end = sph + data->spheres_count; sph < end; sph ++)
    selection_render_sphere_triangles ((*sph)->cur_center, (*sph)->cur_radius);
}

/* render rough mesh */
static void render_rough_mesh (BODY *bod)
{
  BODY_DATA *data = bod->rendering,
	    *rough = data->rough;

  if (data->flags & HIDDEN) return;

  if (!rough)
  {
    SHAPE shape = {SHAPE_MESH, bod->msh, NULL, NULL};
    BODY body = bod [0];
    body.shape = &shape;
    body.msh = NULL;
    data->rough = create_body_data (&body);
    rough = data->rough;
  }

  glColor4f (0.0, 0.0, 0.0, 0.2);

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, rough->lines);

  glEnableClientState (GL_VERTEX_ARRAY);

  glVertexPointer (3, GL_FLOAT, 0, 0);

  glDrawArrays (GL_LINES, 0, rough->lines_count * 2);

  glDisableClientState (GL_VERTEX_ARRAY);

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);

  glColor4f (0.9, 0.9, 0.9, 0.3);

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, rough->triangles);

  glEnableClientState (GL_VERTEX_ARRAY);
  glEnableClientState (GL_NORMAL_ARRAY);

  glVertexPointer (3, GL_FLOAT, 0, 0);
  glNormalPointer (GL_FLOAT, 0, (void*) (rough->triangles_count * sizeof (GLfloat) * 9));

  glDrawArrays (GL_TRIANGLES, 0, rough->triangles_count * 3);

  glDisableClientState (GL_VERTEX_ARRAY);
  glDisableClientState (GL_NORMAL_ARRAY);

  glBindBufferARB (GL_ARRAY_BUFFER_ARB, 0);
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
static void render_riglnk (CON *con, GLfloat color [3])
{
  double other [3];

  ADD (con->point, RIGLNK_VEC (con->Z), other);

  glColor3fv (color);

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

/* render tangential reactions */
static void render_rt (CON *con, GLfloat color [3])
{
  double r [3],
	 other [3],
	 ext = GLV_Minimal_Extent() * arrow_factor,
	 eps, len;

  COPY (con->base, r);
  SCALE (r, con->R[0]);
  ADDMUL (r, con->R[1], con->base+3, r);
  len = LEN (r);
  if (len == 0.0) len = 1.0;
  eps = (ext  / len) * (1.0 + (len - legend.extents[0]) / (legend.extents[1] - legend.extents[0] + 1.0));
  ADDMUL (con->point, -eps, r, other);

  glColor3fv (color);
  arrow3d (other, con->point);
}

/* render normal reactions */
static void render_rn (CON *con, GLfloat color [3])
{
  double r [3],
	 other [3],
	 ext = GLV_Minimal_Extent() * arrow_factor,
	 eps, len;

  COPY (con->base + 6, r);
  SCALE (r, con->R[2]);
  len = LEN (r);
  if (len == 0.0) len = 1.0;
  eps = (ext  / len) * (1.0 + (len - legend.extents[0]) / (legend.extents[1] - legend.extents[0] + 1.0));
  ADDMUL (con->point, -eps, r, other);

  glColor3fv (color);
  arrow3d (other, con->point);
}

/* render reaction */
static void render_r (CON *con, GLfloat color [3])
{
  double r [3],
	 other [3],
	 ext = GLV_Minimal_Extent() * arrow_factor,
	 eps,
	 len;

  COPY (con->base, r);
  SCALE (r, con->R[0]);
  ADDMUL (r, con->R[1], con->base+3, r);
  ADDMUL (r, con->R[2], con->base+6, r);
  len = LEN (r);
  if (len == 0.0) len = 1.0;
  eps = (ext  / len) * (1.0 + (len - legend.extents[0]) / (legend.extents[1] - legend.extents[0] + 1.0));
  ADDMUL (con->point, -eps, r, other);

  glColor3fv (color);
  arrow3d (other, con->point);
}

/* render force constraint */
static void render_force (BODY *bod, FORCE *force, GLfloat color [3])
{
  double r [3],
	 other [3],
	 point [3],
	 ext = GLV_Minimal_Extent() * arrow_factor,
	 eps,
	 len;

  if (bod->kind == FEM)
  {
    SGP *sgp;
    int n;

    if ((n = SHAPE_Sgp (bod->sgp, bod->nsgp, force->ref_point)) < 0) return; /* TODO: optimize */
    BODY_Cur_Point (bod, sgp->shp, sgp->gobj, force->ref_point, point); /* TODO: optimize */
  }
  else BODY_Cur_Point (bod, NULL, NULL, force->ref_point, point);

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

/* render body set constraints or forces */
static void render_body_set_constraints_or_forces (SET *set)
{
  GLfloat color [3];
  SET *item, *jtem;
  BODY_DATA *data;
  FORCE *force;
  BODY *bod;
  CON *con;

  for (con = domain->con; con; con = con->next) con->state &= ~CON_DONERND; /* all undone */

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

	if (con->state & CON_DONERND) continue;

	switch (legend.entity)
	{
	case KINDS_OF_CONSTRAINTS:

	  value_to_color (con->kind, color);

	  switch (con->kind)
	  {
	    case CONTACT: render_contact (con, color); break;
	    case FIXPNT: render_fixpnt (con, color); break;
	    case FIXDIR: render_fixdir (con, color); break;
	    case VELODIR: render_velodir (con, color); break;
	    case RIGLNK: render_riglnk (con, color); break;
	  }

	  break;
	case RESULTS_RT: value_to_color (LEN2 (con->R), color); render_rt (con, color); break;
	case RESULTS_RN: value_to_color (fabs (con->R[2]), color); render_rn (con, color); break;
	case RESULTS_R:  value_to_color (LEN (con->R), color); render_r (con, color); break;
	}

	con->state |= CON_DONERND;
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

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    COPY (neutral_color, color);

    for (item = SET_First (set); item; item = SET_Next (item))
    {
      if (item->data == picked_body) continue;

      glColor4fv (color);
      render_body_lines (item->data, HIDDEN);
      selection_render_body_triangles (item->data);
    }

    glDisable (GL_BLEND);
  }
  else /* regular rendering */
  {
    glDisable (GL_LIGHTING);
    glColor3fv (color);

    for (item = SET_First (set); item; item = SET_Next (item))
    {
      render_body_lines (item->data, TRANSPARENT|HIDDEN);
    }

    glEnable (GL_LIGHTING);
    glEnable (GL_POLYGON_OFFSET_FILL);
    glPolygonOffset (1.0, 1.0);

    for (item = SET_First (set); item; item = SET_Next (item))
    {
      render_body_triangles (item->data, TRANSPARENT|HIDDEN);
    }

    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    COPY (neutral_color, color);

    for (item = SET_First (set); item; item = SET_Next (item))
    {
      if (item->data == picked_body) continue;

      if (TRANSPARENT (item->data))
      {
	glColor4fv (color);
	render_body_lines (item->data, HIDDEN);
	selection_render_body_triangles (item->data);
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
    selection_render_body_triangles (item->data);
  }
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
    selection_render_body_triangles (item->data);
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
    if (TRANSPARENT (picked_body))
    {
      glColor3f (0., 0., 0.);
      render_body_lines (picked_body, 0);
    }
    glColor3fv (color);
    selection_render_body_triangles (picked_body);
    glEnable (GL_LIGHTING);
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
  if (legend.entity)
  {
    legend.extents [0] =  DBL_MAX;
    legend.extents [1] = -DBL_MAX;

    SET_Free (&rndsetmem, &legend.discrete);

    if (legend_constraint_based ()) 
    {
      for (CON *con = domain->con; con; con = con->next) con->state &= ~CON_DONERND; /* all undone */
    }

    for (SET *item = SET_First (selection->set); item; item = SET_Next (item))
    {
      BODY *bod = item->data;
      update_body_values (bod, bod->rendering);
    }

    legend_enable ();
  }

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
  if (solfec->mode == SOLFEC_WRITE) step ();
  else 
  {
    SOLFEC_Forward (solfec, skip_steps);
    update ();
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

/* disable all modes */
static int modes_off ()
{
  int ret = tool_mode;
  mouse_mode = MOUSE_NONE;
  selection_mode = SELECTION_NONE;
  tool_mode = 0;
  picked_body = NULL;
  GLV_Rectangle_Off ();
  GLV_Release_Mouse ();
  update ();
  return ret;
}

/* time window width */
static int time_width ()
{
  return GLV_Print_Width (TIME_FONT, "t=%g", domain->time) + 6;
}

/* render current time */
static void time_render ()
{
  GLV_Resize_Viewport (time_window, time_width (), TIME_HEIGHT);

  glDisable (GL_LIGHTING);
  glDisable (GL_DEPTH_TEST);

  glColor3f (1, 1, 1);
  glRecti (0, 0, time_width (), TIME_HEIGHT);
  glColor3f (0, 0, 0);
  GLV_Print (3, 3, 0, TIME_FONT, "t=%g", domain->time);

  glEnable (GL_LIGHTING);
  glEnable (GL_DEPTH_TEST);
}

/* menus */

static void menu_domain (int item)
{
  int prevmode = solfec->mode;
  DOM *prevdom = domain;
  
  switch (item)
  {
  case DOMAIN_NEXT:
    if (domain->next) domain = domain->next; break;
  case DOMAIN_PREVIOUS:
    if (domain->prev) domain = domain->prev; break;
  }

  /* update analysis menu */

  if (domain != prevdom)
  {
    selection_init ();

    update_extents ();

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
    legend_disable ();
    modes_off ();
    mouse_mode = MOUSE_SELECTION_BEGIN;
    selection_mode = SELECTION_2D;
    GLV_Hold_Mouse ();
    break;
  case RENDER_SELECTION_3D:
    legend_disable ();
    modes_off ();
    mouse_mode = MOUSE_SELECTION_BEGIN;
    selection_mode = SELECTION_3D;
    GLV_Hold_Mouse ();
    break;
  case RENDER_PREVIOUS_SELECTION:
    selection_pop ();
    update_extents ();
    update ();
    break;
  }
}

static void menu_tools (int item)
{
  switch (item)
  {
  case TOOLS_TRANSPARENT:
  case TOOLS_ROUGH_MESH:
  case TOOLS_HIDE:
    modes_off ();
    mouse_mode = MOUSE_PICK_BODY;
    tool_mode = item;
    GLV_Hold_Mouse ();
    break;
  case TOOLS_SHOW_ALL:
    for (BODY *bod = domain->bod; bod; bod = bod->next)
    { BODY_DATA *data = bod->rendering; data->flags &= ~HIDDEN; }
    update_extents ();
    break;
  case TOOLS_NEXT_RESULT:
    if (legend.entity)
    {
      if (legend.entity < RESULTS_R)
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
    else menu_results (RESULTS_R);
    break;
  case TOOLS_SMALLER_ARROWS:
    if (arrow_factor > 0.02) arrow_factor -= 0.01;
    GLV_Redraw_All ();
    break;
  case TOOLS_BIGGER_ARROWS:
    if (arrow_factor < 0.10) arrow_factor += 0.01;
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
  glutAddMenuEntry ("previous selection /p/", RENDER_PREVIOUS_SELECTION);

  menu_name [MENU_TOOLS] = "tools";
  menu_code [MENU_TOOLS] = glutCreateMenu (menu_tools);
  glutAddMenuEntry ("toggle transparent /t/", TOOLS_TRANSPARENT);
  glutAddMenuEntry ("hide /h/", TOOLS_HIDE);
  glutAddMenuEntry ("show all /a/", TOOLS_SHOW_ALL);
  glutAddMenuEntry ("toggle rough mesh /r/", TOOLS_ROUGH_MESH);
  glutAddMenuEntry ("next result /+/", TOOLS_NEXT_RESULT);
  glutAddMenuEntry ("previous result /-/", TOOLS_PREVIOUS_RESULT);
  glutAddMenuEntry ("bigger arrows /]/", TOOLS_BIGGER_ARROWS);
  glutAddMenuEntry ("smaller arrows /[/", TOOLS_SMALLER_ARROWS);

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

  render_picked_body ();
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
  case 'r':
    menu_tools (TOOLS_ROUGH_MESH);
    break;
  case '=':
    menu_tools (TOOLS_NEXT_RESULT);
    break;
  case '-':
    menu_tools (TOOLS_PREVIOUS_RESULT);
    break;
  case '[':
    menu_tools (TOOLS_SMALLER_ARROWS);
    break;
  case ']':
    menu_tools (TOOLS_BIGGER_ARROWS);
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
	case SELECTION_2D: select_2D (mouse_start [0], mouse_start [1], x, y); break;
	case SELECTION_3D: select_3D (mouse_start [0], mouse_start [1], x, y); break;
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
	  if (data->flags & TRANSPARENT) data->flags &= ~TRANSPARENT;
	  else data->flags |= TRANSPARENT;
	  break;
	case TOOLS_HIDE: data->flags |= HIDDEN; break;
	case TOOLS_ROUGH_MESH:
	  if (data->flags & ROUGH_MESH) data->flags &= ~ROUGH_MESH;
	  else if (picked_body->msh) data->flags |= ROUGH_MESH;
	  break;
	}

	GLV_Redraw_All ();
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

  if (data->rough) RND_Free_Rendering_Data (data->rough);

  free (data);
}
