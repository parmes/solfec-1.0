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

/* write uncompressed domain state */
void dom_write_state (DOM *dom, PBF *bf);

/* read uncompressed domain state */
void dom_read_state (DOM *dom, PBF *bf);

/* read uncompressed state of an individual body */
int dom_read_body (DOM *dom, PBF *bf, BODY *bod);

/* read uncompressed state of an individual constraint */
int dom_read_constraint (DOM *dom, PBF *bf, CON *con);

#endif
