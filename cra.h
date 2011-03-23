/*
 * cra.h
 * Copyright (C) 2011, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * body cracking
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

#ifndef DOMAIN_TYPE
#define DOMAIN_TYPE
typedef struct domain DOM;
#endif

#include "cvx.h" /* ELEPNT */

#ifndef __cra__
#define __cra__

typedef struct crack CRACK;

struct crack
{
  double point [3],
	 normal [3];

  enum {TENSILE} crit;

  int surfid,
      nepn;

  double ft,
	 Gf;

  ELEPNT *epn; /* for cuts through FEM meshes */

  CRACK *next;
};

/* create crack object */
CRACK* CRACK_Create ();

/* delete crack object */
void CRACK_Destroy (CRACK *cra);

/* delete list of cracks */
void CRACK_Destroy_List (CRACK *cra);

/* propagate cracks and adjust the domain */
void Propagate_Cracks (DOM *dom);

#endif
