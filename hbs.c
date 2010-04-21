/*
 * nts.c
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * Hybrid constraints solver
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

#include "hbs.h"
#include "err.h"

#define ENABLE 0 /* temporary flag before the code compiles */

#if ENABLE
typedef struct body_data BD;

struct body_data
{
  BODY *bod;

  CON **con; /* copy of bod->con */

  MX **H; /* corresponding H operators */

  int n; /* size of bod->con */

  double *b; /* free velocity plus previous tangential contact response */

  double *u; /* free velocity projected on the constraints tangent set */

  BD *next;
};

typedef struct body_data_pair BDP;

struct body_data_pair
{
  BD *bd [2]; /* two body data sets */

  CON **con [2]; /* two subsets of constraints that the corresponding bodies share */

  int n [2]; /* sizes of the constraint subsets */

#if MPI
  int rank [2]; /* rank[i] >= 0 if db[i] corresponds to a remote dummy body */
#endif

  BDP *next;
};

/* create a list of body data items */
static BD* body_data_create (DOM *dom)
{
  /* TODO: remember to scale H by 1/h */
}

/* destroy a list of body data items */
static void body_data_destroy (BD *data)
{
}

/* create a list of body data pairs */
static DBP* body_data_pairs_create (BD *data)
{
#if MPI
  /* TODO: bring here needed dummy body data */
#endif
}

/* destroy a list of body data pairs */
static void body_data_pairs_destroy (DBP *pairs)
{
}

/* create and return H and d such that H * u <= d describes the constraints tangent space;
 * initialize db->b with free velocity summed up with previous tangential contact response */
static MX* velocity_projection_setup (BD *bd, double **d)
{
} 

/* approximate db->u = proj (db->b on T), where T: H * u <= d */
static void velocity_projection (DB *db)
{
}

/* create and return [H_i'; H_j'] and [M_i(u_i - b_i); M_j(u_j - b_j)]' */
static MX* normal_reactions_setup (DBP *dp, double **c)
{
}

/* compute R0 such that [M_i(u_i - b_i); M_j(u_j - b_j)]' = [H_i'; H_j'] R0;
 * for normal contact reactions do RN_i = proj (R0_i on R+); otherwise R_j = R0_j */
static void normal_reactions (DBP *dp)
{
}

/* compute tangential contact response */
static void tangential_reactions (DB *db)
{
}

/* compute single-body equality constraints */
static void equality_constraints (DB *db)
{
}
#endif

/* create solver */
HYBRID* HYBRID_Create (double epsilon, int maxiter)
{
  HYBRID *hb;

  ERRMEM (hb = malloc (sizeof (HYBRID)));
  hb->epsilon = epsilon;
  hb->maxiter = maxiter;

  return hb;
}

/* run solver */
void HYBRID_Solve (HYBRID *hb, DOM *dom)
{
#if ENABLE
  DBP *pairs, *dp;
  DB *data, *db;

  data = body_data_create (dom);

  pairs = body_data_pairs_create (data);

  for (db = data; db; db = db->next)
  {
    velocity_projection (db);
  }

#if MPI
  /* TODO: send DB->u to remote DBP->db[]->u */
#endif

  for (dp = pairs; dp; dp = dp->next)
  {
    normal_reactions (dp);
  }

#if MPI
  /* TODO: send normal reactions to remote bodies */
#endif

  for (db = data; db; db = db->next)
  {
    tangential_reactions (db);
  }

  for (db = data; db; db = db->next)
  {
    equality_constraints (db);
  }

  body_data_pairs_destroy (pairs);

  body_data_destroy (data);
#endif
}

/* write labeled satate values */
void HYBRID_Write_State (HYBRID *hb, PBF *bf)
{
}

/* destroy solver */
void HYBRID_Destroy (HYBRID *hb)
{
  free (hb);
}
