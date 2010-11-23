/*
 * mrf.h
 * Copyright (C) 2010 Tomasz Koziara (t.koziara AT gmail.com)
 * -------------------------------------------------------------------
 * constraints satisfaction merit function
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

#include "ldy.h"

#ifndef __mrf__
#define __mrf__

/* constraint satisfaction merit function approximately indicates the
 * amount of spurious momentum due to constraint force inaccuracy;
 * update_U != 0 implies that U needs to be computed for current R;
 * (it is assumed that all (also external) reactions are updated) */
double MERIT_Function (LOCDYN *ldy, short update_U);

#endif
