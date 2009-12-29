/*
 * bod.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * general body
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
#include "pbf.h"
#include "shp.h"
#include "tms.h"
#include "set.h"
#include "mtx.h"
#include "mat.h"
#include "msh.h"

#ifndef DOMAIN_TYPE
#define DOMAIN_TYPE
typedef struct domain DOM;
#endif

#ifndef SOLFEC_TYPE
#define SOLFEC_TYPE
typedef struct solfec SOLFEC;
#endif

#ifndef __bod__
#define __bod__

typedef struct general_force FORCE;
typedef void (*FORCE_FUNC) (void *data, void *call, /* user data and user callback pointers */
                            int nq, double *q, int nu, double *u,   /* user defined data, configuration, velocity, time, time step */
                            double t, double h, double *f);  /* for rigid bodies 'f' comprises [spatial force; spatial torque; referential torque];
                                                                for other types of bodies 'f' is the generalised force */
#ifndef BODY_TYPE
#define BODY_TYPE
typedef struct general_body BODY;
#endif

/* results
 * value kinds */
typedef enum
{
  VALUE_DISPLACEMENT,
  VALUE_VELOCITY,
  VALUE_STRESS,
  VALUE_MISES,
  VALUE_STRESS_AND_MISES
} VALUE_KIND;

/* time integration
 * schemes */
typedef enum 
{
  SCH_RIG_POS,   /* NEW1 with positive energy drift (high accuracy, approximate momentum conservation) */
  SCH_RIG_NEG,   /* NEW2 with with negative energy drift (exact momentum conservation) (DEFAULT) */
  SCH_RIG_IMP,   /* NEW3 semi-simplict and stable (no energy drift, extact momentum conservation) */
  SCH_DEF_EXP,   /* explicit scheme for deformable kinematics (DEFAULT) */
  SCH_DEF_IMP,   /* linearly implicit scheme for deformable kinematics */
} SCHEME;

struct general_force
{
  enum {SPATIAL   = 0x01,
        CONVECTED = 0x02,
        TORQUE    = 0x04} kind; /* force kind (torque applies only to rigid bodies) */

  double ref_point [3],         /* referential point */
         direction [3];         /* spatial or referential */

  TMS *data;

  void *call;

  FORCE_FUNC func;

  FORCE *next;
};

/* energy kinds */
#define EXTERNAL 0
#define KINETIC 1
#define INTERNAL 2
#define BODY_ENERGY_SIZE(bod) ((bod)->kind == OBS ? 0 : (bod)->kind == RIG ? 2 : 3)

/* body flags */
typedef enum
{
  BODY_DETECT_SELF_CONTACT = 0x01, /* enable self contact detection */
  BODY_PARENT              = 0x02, /* a parent body */
  BODY_CHILD               = 0x04, /* a child body */
  BODY_CHILD_UPDATED       = 0x08, /* an updated child */
} BODY_FLAGS;

/* flags that are migrated with bodies (the rest is filtered out) */
#define BODY_PERMANENT_FLAGS (BODY_DETECT_SELF_CONTACT)

struct general_body
{
  enum {OBS, RIG, PRB, FEM} kind; /* obstacle, rigid, pseudo-rigid, finite element */

  unsigned int id;  /* unique identifier (for serialization & parallel processing) */

  BULK_MATERIAL *mat; /* default material */

  double ref_mass,
         ref_volume,
         ref_center [3],
         ref_tensor [9]; /* RIG => Inertia tensor
			    PRB => Euler tensor */
  double *conf,
         *velo;

  int dofs;        /* velocity degrees of freedom count */

  SET *con;        /* adjacent constraints */

  FORCE *forces;   /* applied external forces */
  
  SHAPE *shape;    /* shape of the body */

  SGP *sgp;        /* shape and geometric object pairs */

  int nsgp;        /* number of them (above) */

  double extents [6];  /* shape extents */

  SCHEME scheme;    /* integration scheme */

  MX *inverse;      /* a suitable inverse oprator (e.g. inertia for dynamics) */

  double damping;   /* mass proportional damping */

  DOM *dom;        /* domain storing the body */

  BODY *prev,       /* list */
       *next;

  char *label;      /* user specified label */

  BODY_FLAGS flags;      /* flags */

  short form; /* FEM formulation */

  MESH *msh; /* FEM mesh when 'shape' is made of CONVEX objects ("rough mesh") */

  double energy [3]; /* external, kinetic, potential */

#if MPI
  SET *children; /* set of children ids for a parent; set of other children for a child */

  int rank; /* parent => new/current rank; child => parent's rank */
#else
  void *rendering; /* rendering data */
#endif
};

/* body pointer cast */
#define BODY(bod) ((BODY*)(bod))

/* create a body */
BODY* BODY_Create (short kind, SHAPE *shp, BULK_MATERIAL *mat, char *label, short form, MESH *msh);

/* get body kind string */
char* BODY_Kind (BODY *bod);

/* get configuration size */
int BODY_Conf_Size (BODY *bod);

/* overwrite mass and volume characteristics */
void BODY_Overwrite_Chars (BODY *bod, double mass, double volume, double *center, double *tensor);

/* overwrite body state */
void BODY_Overwrite_State (BODY *bod, double *q, double *u);

/* apply an initial rigid motion velocity */
void BODY_Initial_Velocity (BODY *bod, double *linear, double *angular);

/* apply a force (if 'func' is given, 'data' is regarded as the user data pointer to the callback 'func') */
void BODY_Apply_Force (BODY *bod, short kind, double *point, double *direction, TMS *data, void *call, FORCE_FUNC func);

/* remove all forces */
void BODY_Clear_Forces (BODY *bod);

/* set new mapterial */
void BODY_Material (BODY *bod, int volume, BULK_MATERIAL *mat);

/* initialise dynamic time stepping */
void BODY_Dynamic_Init (BODY *bod);

/* estimate critical step for the dynamic scheme */
double BODY_Dynamic_Critical_Step (BODY *bod);

/* perform the initial half-step of the dynamic scheme */
void BODY_Dynamic_Step_Begin (BODY *bod, double time, double step);

/* perform the final half-step of the dynamic scheme */
void BODY_Dynamic_Step_End (BODY *bod, double time, double step);

/* initialise static time stepping */
void BODY_Static_Init (BODY *bod);

/* perform the initial half-step of the static scheme */
void BODY_Static_Step_Begin (BODY *bod, double time, double step);

/* perform the final half-step of the static scheme */
void BODY_Static_Step_End (BODY *bod, double time, double step);

/* motion x = x (X, state) */
void BODY_Cur_Point (BODY *bod, SHAPE *shp, void *gobj, double *X, double *x);

/* inverse motion X = X (x, state) */
void BODY_Ref_Point (BODY *bod, SHAPE *shp, void *gobj, double *x, double *X);

/* obtain spatial velocity at (gobj, referential point), expressed in the local spatial 'base' */
void BODY_Local_Velo (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base, double *prevel, double *curvel);

/* return transformation operator from the generalised to the local velocity space at (element, point, base) */
MX* BODY_Gen_To_Loc_Operator (BODY *bod, SHAPE *shp, void *gobj, double *point, double *base);

/* compute current kinetic energy */
double BODY_Kinetic_Energy (BODY *bod);

/* get some values at a referential point */
void BODY_Point_Values (BODY *bod, double *point, VALUE_KIND kind, double *values);

/* write body state */
void BODY_Write_State (BODY *bod, PBF *bf);

/* read body state */
void BODY_Read_State (BODY *bod, PBF *bf);

/* pack body state */
void BODY_Pack_State (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack body state */
void BODY_Unpack_State (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

/* release body memory */
void BODY_Destroy (BODY *bod);

/* pack body into double and integer buffers (d and i buffers are of initial
 * dsize and isize, while the final numberof of doubles and ints is packed) */
void BODY_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);

/* unpack body from double and integer buffers (unpacking starts at dpos and ipos in
 * d and i and no more than a specific number of doubles and ints can be red) */
BODY* BODY_Unpack (SOLFEC *sol, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

#if MPI
/* parent bodies store all body data and serve for time stepping */
void BODY_Parent_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);
void BODY_Parent_Unpack (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

/* child bodies store a minimal subset of needed data and serve for constraint solution */
void BODY_Child_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);
void BODY_Child_Unpack (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);

/* child body updates pack and unpack configurations and update shapes */
void BODY_Child_Update_Pack (BODY *bod, int *dsize, double **d, int *doubles, int *isize, int **i, int *ints);
void BODY_Child_Update_Unpack (BODY *bod, int *dpos, double *d, int doubles, int *ipos, int *i, int ints);
#endif
#endif
