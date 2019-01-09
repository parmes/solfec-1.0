/*
 * dio.h
 * Copyright (C) 2009, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * domain input-output
 */

#include "dom.h"
#include "cmp.h"

#ifndef __dio__
#define __dio__

/* write domain state */
void dom_write_state (DOM *dom, PBF *bf, SET *subset);

/* read domain state */
void dom_read_state (DOM *dom, PBF *bf);

/* read state of an individual body */
int dom_read_body (DOM *dom, PBF *bf, BODY *bod);

/* read state of an individual constraint */
int dom_read_constraint (DOM *dom, PBF *bf, CON *con);

/* initialize domain state */
int dom_init_state (DOM *dom, PBF *bf, SET *subset);

/* map rigid onto FEM state */
int dom_rigid_to_fem (DOM *dom, PBF *bf, SET *subset);

#endif
