/*
 * dyr.h
 * Copyright (C) 2005, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------------
 * dynamic rectangle structure
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


#include "box.h"

#ifndef __dyr__
#define __dyr__

enum dralg /* an approach for the approximation of dynamic rectangle structure */
{
  HASH2D_LIST = 0,
  HASH2D_XYTREE,
  XYTREE_ALONE,
  HASH1D_XYTREE
};

typedef enum dralg DRALG;

/*
 * Create dynamic rectangle (DR) structure return 'dr' pointer.
 */
void* DR_Create (int boxnum, DRALG algo);

/*
 * Set up avarage edge lengths for hashing approaches.
 */
void DR_Params (void *dr, double avglen [2]);

/*
 * Insert a rectnagle into DR, performing the overlap query before.
 */
void DR_Insert (void *dr, double xmin, double ymin, double xmax, double ymax,
                BOX *box, void *data, BOX_Overlap_Create report);
/*
 * Delete a rectangle from the DR.
 */
void DR_Delete (void *dr, double xmin, double ymin, double xmax, double ymax, BOX *box);

/*
 * Clean up.
 */
void DR_Destroy (void *dr);

#endif
