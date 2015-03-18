/*
 * fra.h
 * Copyright (C) 2013, Tomasz Koziara (t.koziara AT gmail.com)
 * ---------------------------------------------------------------
 * fracture code
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

#ifndef __fra__
#define __fra__

typedef struct fracture_state FS;

struct fracture_state
{
  /* instance data */
  double *disp;
  FS *inext; /* instances list */

  /* contact points data */
  double radius;
  double point [3];
  double force [3];
  FS *next; /* contact forces list within instance */
};

/* free list */
void fracture_state_free (FS *list);

/* read fracture state */
FS* fracture_state_read (BODY *bod);

/* check fracture criterion */
void Fracture_Check (DOM *dom);

/* export data for fracture analysis in Yaffems (return number of exported analysis instances) */
int Fracture_Export_Yaffems (BODY *bod, double volume, double quality, FILE *output);

#endif
