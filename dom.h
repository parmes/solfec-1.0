/*
 * dom.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * a domain gathers bodies and constraints
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

#if MPI
#include <zoltan.h>
#include "com.h"
#endif

#include "mem.h"
#include "map.h"
#include "box.h"
#include "sps.h"
#include "bod.h"
#include "ldy.h"
#include "pbf.h"
#include "cmp.h"

#ifndef SOLFEC_TYPE
#define SOLFEC_TYPE
typedef struct solfec SOLFEC;
#endif

#ifndef __dom__
#define __dom__

#ifndef CONSTRAINT_TYPE
#define CONSTRAINT_TYPE
typedef struct constraint CON;
#endif

#ifndef DOMAIN_TYPE
#define DOMAIN_TYPE
typedef struct domain DOM;
#endif

#define DOM_Z_SIZE          4      /* size of auxiliary storage */
#define RIGLNK_VEC(Z)   (Z)        /* rigid link vector */
#define RIGLNK_LEN(Z)   ((Z)[3])   /* rigid link length */
#define VELODIR(Z)      ((Z)[0])   /* prescribed velocity at (t+h) */

struct constraint
{
  double R [3], /* average constraint reaction */
	 U [3]; /* relative velocity */

  DIAB *dia; /* diagonal entry in the local dynamical system */

  double point [3], /* spatial point */
         base [9], /* local orthogonal base */
	 area, /* contact area */
	 gap; /* contact gap */

  double Z [DOM_Z_SIZE]; /* auxiliary storage */

  void *data; /* auxiliary solver data */

  unsigned int id; /* identifier */

  int num; /* one of: 0, 1, ..., total constraints count - 1 */

  enum {CONTACT, FIXPNT, FIXDIR, VELODIR, RIGLNK} kind; /* constraint kind */

  enum {CON_COHESIVE = 0x01,
        CON_NEW      = 0x02,
	CON_IDLOCK   = 0x04, /* locked ID cannot be freed to the pool */
	CON_EXTERNAL = 0x08, /* a boundary constraint migrated in from another processor */
        CON_DONE     = 0x10} state; /* constraint state */

  short paircode; /* geometric object pair code for a contact */

  SURFACE_MATERIAL_STATE mat; /* surface pair material data */

  TMS *tms; /* time series data (if any) */

  BODY *master, /* master body */
       *slave; /* slave body */

  double mpnt [3], /* master referential point */
	 spnt [3]; /* slave referential point */

  SGP *msgp, /* master (shape, gobj, box) triplet */
      *ssgp; /* slave triplet (if any) */

  int rank; /* parallel: origin rank for an external constraint;
               serial read: rank of residence during parallel run */
#if MPI
  SET *ext; /* ranks of remote external images of this constraint */
#endif

  CON *prev, *next; /* list */
};

#define mshp(con) ((con)->msgp->shp)
#define sshp(con) ((con)->ssgp->shp)
#define mgobj(con) ((con)->msgp->gobj)
#define sgobj(con) ((con)->ssgp->gobj)
#define mkind(con) GOBJ_Kind ((con)->msgp)
#define skind(con) GOBJ_Kind ((con)->ssgp)

#if MPI
typedef struct domain_statistics DOMSTATS;

/* parallel statistics */
struct domain_statistics
{
  char *name;
  int sum;
  int min;
  int avg;
  int max;
};

typedef struct domain_balancing_data DBD;

/* load balancing send sets */
struct domain_balancing_data
{
  int rank;
  DOM *dom;
  SET *bodies;
  SET *children;
  SET *constraints;
  SET *remove;
  SET *update;
  SET *glue;
  SET *ext;
};
#endif

/* box overlap algorithm selection data */
typedef struct aabb_data AABB_DATA;

struct aabb_data
{
  double aabb_timings [BOXALG_COUNT],
	 aabb_limits [BOXALG_COUNT+1];

  int aabb_counter;

  BOXALG aabb_algo;
};

/* domain flags */
typedef enum
{
  DOM_RUN_ANALYSIS   = 0x01, /* on when the viewer runs analysis for this domain */
  DOM_DEPTH_VIOLATED = 0x02  /* on when unphysical penetration has occured */
} DOM_FLAGS;

/* domain data */
struct domain
{
  MEM conmem, /* constraints memory pool */
      mapmem, /* map items memory pool */
      setmem; /* set items memory pool */

  AABB *aabb; /* box overlap engine */
  SPSET *sps; /* surface pairs */

  short dynamic; /* 1 for dynamics, 0 for quas-statics */
  double step; /* time step size */
  double time; /* current time */

  unsigned int bid; /* last free body identifier */
  SET *sparebid; /* spare constraint ids */
  MAP *lab; /* bodies mapped by labels (a subset) */
  MAP *idb; /* bodies mapped by identifiers (all bodies) */
  BODY *bod; /* list of bodies */
  int nbod; /* number of bodies */
  SET *delb; /* set of deleted body ids for time > 0 and before state write */
  SET *newb; /* set of newly created bodies for time > 0 and before state write */

  unsigned int cid;  /* last free constraint identifier */
  SET *sparecid; /* spare constraint ids */
  MAP *idc; /* constraints mapped by identifiers */
  CON *con; /* list of constraints */
  int ncon; /* number of constraints */
  int nspa; /* number of sparsified contacts */

  LOCDYN *ldy; /* local dynamics */
  SOLFEC *solfec; /* SOLFEC context */
  AABB_DATA *aabb_data; /* box ovrlap algorithm selection data */

  TMS *gravity [3]; /* global gravity value */

  double extents [6]; /* scene extents */
  double threshold; /* sparsification threshold */
  double minarea; /* minimal contact point area */
  double depth; /* unphisical interpenetration depth bound (negative) */

  DOM_FLAGS flags; /* some flags */
  short verbose; /* verbosity flag */

#if MPI
  int rank; /* communicator rank */
  int ncpu; /* cummunicator size */
  MAP *allbodies; /* all created bodies mapped by ids (minimises load balancing communication) */
  SET *children; /* current children */
  struct Zoltan_Struct *zol; /* load balancing */
  double imbalance_tolerance; /* imbalance threshold */
  int lock_directions; /* locked direactions flag */
  double degenerate_ratio; /* degeneration ratio for domain elongation */
  double weight_factor; /* local dynamics weight factor */
  unsigned int noid; /* constraint id generation ommition flag */
  MAP *conext; /* id based map of external constraints */
  int bytes; /* bytes sent during load balancing */
  DOMSTATS *stats; /* domain statistics */
  int nstats; /* statistics count */
  DBD *dbd; /* load balancing send sets */
#endif

  DOM *prev, *next; /* list */
};

/* constraint kind string */
char* CON_Kind (CON *con);

/* constraint pointer cast */
#define CON(con) ((CON*)con)

/* create a domain => 'aabb' is the box overlap solver, 'sps' is the
 * surface pair set, 'dynamic' == 1  indicates a dynamic simulation
 * (quasi-static otherwise), finally 'step' is the time step */
DOM* DOM_Create (AABB *aabb, SPSET *sps, short dynamic, double step);

/* a body needs to be inserted into the domain, before any other routine
 * involving the body pointer can be called; inserting bodies into the
 * domain enables automatic maintenance of contacts constraints */

/* insert a body into the domain */
void DOM_Insert_Body (DOM *dom, BODY *bod);

/* remove a body from the domain (do not destroy it) */
void DOM_Remove_Body (DOM *dom, BODY *bod);

/* find labeled body */
BODY* DOM_Find_Body (DOM *dom, char *label);

/* fix a referential point of the body along all directions */
CON* DOM_Fix_Point (DOM *dom, BODY *bod, double *pnt);

/* fix a referential point of the body along the spatial direction */
CON* DOM_Fix_Direction (DOM *dom, BODY *bod, double *pnt, double *dir);

/* prescribe a velocity of the referential point along the spatial direction */
CON* DOM_Set_Velocity (DOM *dom, BODY *bod, double *pnt, double *dir, TMS *vel);

/* insert rigid link constraint between two (referential) points of bodies; if one of the body
 * pointers is NULL then the link acts between the other body and the fixed (spatial) point */
CON* DOM_Put_Rigid_Link (DOM *dom, BODY *master, BODY *slave, double *mpnt, double *spnt);

/* remove a constraint from the domain (destroy it) */
void DOM_Remove_Constraint (DOM *dom, CON *con);

/* set simulation scene extents */
void DOM_Extents (DOM *dom, double *extents);

/* go over contact points and remove those whose corresponding
 * areas are much smaller than those of other points related to
 * objects directly topologically adjacent in their shape definitions */
void DOM_Sparsify_Contacts (DOM *dom);

/* domain update initial half-step => bodies and constraints are
 * updated and the unupdated local dynamic problem is returned */
LOCDYN* DOM_Update_Begin (DOM *dom);

/* (end update of local dynamics before calling this routine)
 * domain update final half-step => once the local dynamic
 * problem has been solved (externally), motion of bodies
 * is updated with the help of new constraint reactions */
void DOM_Update_End (DOM *dom);

#if MPI
/* send boundary reactions to their external receivers;
 * if 'normal' is > 0 only normal components are sent */
void DOM_Update_External_Reactions (DOM *dom, short normal);
#endif

/* assign con->num values */
void DOM_Number_Constraints (DOM *dom);

/* write domain state */
void DOM_Write_State (DOM *dom, PBF *bf, CMP_ALG alg);

/* read domain state */
void DOM_Read_State (DOM *dom, PBF *bf, CMP_ALG alg);

/* read state of an individual body */
int  DOM_Read_Body (DOM *dom, PBF *bf, BODY *bod);

/* read state of an individual constraint */
int  DOM_Read_Constraint (DOM *dom, PBF *bf, CON *con);

/* release memory */
void DOM_Destroy (DOM *dom);
#endif
