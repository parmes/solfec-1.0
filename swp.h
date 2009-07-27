/*
 * swp.h
 * Copyright (C) 2005, 2009 Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------------
 * sweep plane box intersection
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
#include "dyr.h"

#ifndef __swp__
#define __swp__

/*
 * Initialize seep-plane alogorithm, using a DR algo.
 * Return the algorithm context.
 */
void* SWEEP_Create (int boxnum, DRALG algo);

/*
 * Perform overlap detection test.
 */
void SWEEP_Do (void *context, DRALG algo, int changed, int boxnum, BOX **boxes, void *data, BOX_Overlap_Create report);

/*
 * Clean up.
 * Release context data (context no more valid).
 */
void SWEEP_Destroy (void *context);

#endif
