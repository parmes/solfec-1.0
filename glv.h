/*
 * glv.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * simple opengl viewer
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

#ifndef __glv__
#define __glv__

typedef int  (*View_Menu) (char***, int**); /* menu set up (CALLED FIRST) */
typedef void (*View_Init) (void); /* initialisation (CALLED SECOND) */
typedef int  (*View_Idle) (void); /* idle callback => return != 0 to cose view update */
typedef void (*View_Quit) (void); /* on viewer exit */
typedef void (*View_Render) (void); /* window rendering */
typedef void (*View_Key) (int key, int x, int y); /* keyboard */
typedef void (*View_Keyspec) (int key, int x, int y); /* special keyboard */
typedef void (*View_Mouse) (int button, int state, int x, int y); /* mouse */
typedef void (*View_Motion) (int x, int y); /* mouse motion */
typedef void (*View_Passive_Motion) (int x, int y); /* passive mouse motion */

/* create main window */
void GLV (
  int *argc,
  char **argv,
  char *title,
  int width,
  int height,
  double *extents,
  View_Menu menu,
  View_Init init,
  View_Idle idle,
  View_Quit quit,
  View_Render render,
  View_Key key,
  View_Keyspec keyspec,
  View_Mouse mouse,
  View_Motion motion,
  View_Passive_Motion passive);


/* redraw all */
void GLV_Redraw_All (void);

/* reset scene extents */
void GLV_Reset_Extents (double *extents);

/* update scene extents */
void GLV_Update_Extents (double *extents);

/* get minimal view extent */
double GLV_Minimal_Extent ();

/* open a 2D window with 1-to-1
 * coordinate-to-pixel mapping;
 * window's content is not exported
 * in BMP and AVI format */
int GLV_Open_Window (
  int x,
  int y,
  int w,
  int h,
  View_Render render);

/* close window */
void GLV_Close_Window (int window);

/* open a viewport whose content
 * is exported in BMP and AVI format */
int GLV_Open_Viewport (
  int x,
  int y,
  int w,
  int h,
  int is3D,
  View_Render render);

/* move viewport */
void GLV_Move_Viewport (int viewport, int x, int y, int w, int h);

/* close window */
void GLV_Close_Viewport (int viewport);

/* show tiled text intput window 
 * and return read text by callback */
void GLV_Read_Text (char *title, void (*done) (char *text));

/* output text at specified coordinates */
enum {GLV_FONT_10 = 10, GLV_FONT_12 = 12, GLV_FONT_18 = 18};
void GLV_Print (double x, double y, double z, int font, char *fmt, ...);

/* output screen shot bitmap */
void GLV_Screen_Bitmap (char *path);

/* take over mouse */
void GLV_Hold_Mouse ();

/* release mouse takeover */
void GLV_Release_Mouse ();

/* viewer specific projection matrix (useful for picking) */
void GLV_SetProjectionMatrix (int w, int h);

/* enable drawing rectangle in screen coordinates */
void GLV_Rectangle_On (int x1, int y1, int x2, int y2);

/* disable drawing rectangle */
void GLV_Rectangle_Off ();

#endif
