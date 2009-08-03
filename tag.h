/*
 * tag.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * communication tags
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

#ifndef __tag__
#define __tag__

enum
{
  TAG_CONAUX,
  TAG_PARENTS,
  TAG_ORPHANS,
  TAG_CHILDREN_DELETE,
  TAG_CHILDREN_INSERT
};

#endif
