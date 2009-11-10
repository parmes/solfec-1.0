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
#endif

#include "mem.h"
#include "map.h"
#include "box.h"
#include "sps.h"
#include "bod.h"
#include "ldy.h"
#include "pbf.h"
#include "cmp.h"

#ifndef __dom__
#define __dom__

typedef struct constraint CON;
typedef struct domain DOM;

#define RIGLNK_VEC(Z)   (Z)        /* rigid link vector */
#define RIGLNK_LEN(Z)   ((Z)[3])   /* rigid link length */
#define VELODIR(Z)      ((Z)[0])   /* prescribed velocity at (t+h) */

struct constraint
{
  double R [3]; /* average constraint reaction */

  DIAB *dia; /* diagonal entry in the local dynamical system */

  double point [3], /* spatial point */
         base [9], /* local orthogonal base */
	 area, /* contact area */
	 gap; /* contact gap */

  double Z [DOM_Z_SIZE]; /* auxiliary storage */

  unsigned int id; /* identifier */

  enum {CONTACT, FIXPNT, FIXDIR, VELODIR, RIGLNK} kind; /* constraint kind */

  enum {CON_STICK    = 0x01,
        CON_OPEN     = 0x02,
        CON_COHESIVE = 0x04,
        CON_NEW      = 0x08,
        CON_DONERND  = 0x10} state; /* constraint state */

  short paircode; /* geometric object pair code for a contact */

  SURFACE_MATERIAL mat; /* surface pair material data */

  TMS *tms; /* time series data */

  BODY *master, /* master body */
       *slave; /* slave body */

  SHAPE *mshp, /* shape of mgobj */
        *sshp; /* shape of sgobj */

  void *mgobj, /* master geometrical object */
       *sgobj; /* slave geometrical object */

  GOBJ mkind, /* kind of mobj */
       skind; /* kind of sobj */

  double mpnt [3], /* master referential point */
	 spnt [3]; /* slave referential point */

  SGP *msgp, /* master (shape, gobj, box) triplet */
      *ssgp; /* slave triplet (if any) */

  CON *prev, /* list */
      *next;
};

#if MPI
typedef struct conext CONEXT;

/* external constraints are those created
 * by child bodies and attached to parents */
struct conext
{
  double R [3], /* reaction */
         point [3],  /* referential point */
         base [9]; /* spatial base */

  unsigned int id; /* constraint id */

  BODY *bod; /* parent body */

  int isma; /* 1 if the paranet body is a master; 0 otherwise */

  SGP *sgp; /* shape geometric object pair */

  CONEXT *next; /* in-domain list for fast deletion */

  int rank; /* child rank */
};
#endif

/* domain flags */
typedef enum
{
  DOM_RUN_ANALYSIS   = 0x01, /* on when the viewer runs analysis for this domain */
  DOM_DEPTH_VIOLATED = 0x02  /* on when unphysical penetration has occured */
} DOM_FLAGS;

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
  MAP *lab; /* bodies mapped by labels (a subset) */
  MAP *idb; /* bodies mapped by identifiers (all bodies) */
  BODY *bod; /* list of bodies */
  int nbod; /* number of bodies */
  SET *delb; /* set of deleted body ids for time > 0 and before state write */
  SET *newb; /* set of newly created bodies for time > 0 and before state write */

  unsigned int cid;  /* last free constraint identifier */
  MAP *idc; /* constraints mapped by identifiers */
  CON *con; /* list of constraints */
  int ncon; /* number of constraints */
  int nspa; /* number of sparsified contacts */

  LOCDYN *ldy; /* local dynamics */

  double gravdir [3]; /* global gravity direction */
  TMS *gravval; /* global gravity value */

  double extents [6]; /* scene extents */

  void *owner; /* SOLFEC */

  void *data; /* private data */

  DOM *prev, *next; /* list */

  DOM_FLAGS flags;

  double threshold; /* sparsification threshold */

  double depth; /* unphisical interpenetration depth bound (negative) */

  short verbose; /* verbosity flag */

#if MPI
  int rank; /* communicator rank */

  int ncpu; /* cummunicator size */

  struct Zoltan_Struct *zol; /* body partitioning */

  double imbalance_tolerance; /* imbalance threshold */

  SET *sparecid; /* spare constraint ids */

  SET *sparebid; /* spare body ids */

  unsigned int noid; /* constraint id generation ommition flag */

  MAP *children; /* id-to-child map */

  SET **delch; /* id sets of children to be deleted for each rank (after removed parent) */

  MEM extmem; /* memory pool of external constraints */

  CONEXT *conext; /* external constraints (parents and children) */

  int nexpbod, /* number of exported parent bodies */
      nexpchild; /* number of exported child bodies */
#endif
};

/* constraint kind string */
char* CON_Kind (CON *con);

/* constraint pointer cast */
#define CON(con) ((CON*)(con))

/* domain pointer cast */
#define DOM(dom) ((DOM*)(dom))

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

/* domain update initial half-step => bodies and constraints are
 * updated and the unupdated local dynamic problem is returned */
LOCDYN* DOM_Update_Begin (DOM *dom);

/* (end update of local dynamics before calling this routine)
 * domain update final half-step => once the local dynamic
 * problem has been solved (externally), motion of bodies
 * is updated with the help of new constraint reactions */
void DOM_Update_End (DOM *dom);

#if MPI
/* balance children according to a given geometric partitioning */
void DOM_Balance_Children (DOM *dom, struct Zoltan_Struct *zol);

/* update children shapes */
void DOM_Update_Children (DOM *dom);
#endif

/* write domain state */
void DOM_Write_State (DOM *dom, PBF *bf, CMP_ALG alg);

/* read domain state */
void DOM_Read_State (DOM *dom, PBF *bf);

/* read state of an individual body */
int  DOM_Read_Body (DOM *dom, PBF *bf, BODY *bod);

/* read state of an individual constraint */
int  DOM_Read_Constraint (DOM *dom, PBF *bf, CON *con);

/* release memory */
void DOM_Destroy (DOM *dom);

#endif
