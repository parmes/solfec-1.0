/*
 * rnd.h
 * Copyright (C) 2008, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * OpenGL rendering
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

#include "dom.h"

#ifndef __rnd__
#define __rnd__

#if __cplusplus
extern "C" { /* C */
#endif

/* functions below match
 * GLUT hooks from glv.h */
int  RND_Menu (char ***names, int **codes);
void RND_Init ();
int  RND_Idle ();
void RND_Quit ();
void RND_Render ();
void RND_Key (int key, int x, int y);
void RND_Keyspec (int key, int x, int y);
void RND_Mouse (int button, int state, int x, int y);
void RND_Motion (int x, int y);
void RND_Passive (int x, int y);

/* enable renedring before opening a viewer */
void RND_Switch_On ();

/* check whether rendering is on */
int  RND_Is_On ();

/* add domain to be rendered */
void RND_Domain (DOM *dom);

/* map solver to a domain */
void RND_Solver (DOM *dom, int kind, void *solver);

/* free body associated rendering data */
void RND_Free_Rendering_Data (void *ptr);

void select_id (SET *bodies);

#if __cplusplus
} /* extern C */
#endif

#endif
