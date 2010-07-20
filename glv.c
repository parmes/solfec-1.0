/*
 * glv.c
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * simple general purpose opengl viewer
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

#if __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <time.h>
#include <math.h>
#include "ext/gl2ps.h"
#include "err.h"
#include "bmp.h"
#include "glv.h"
#include "alg.h"

#define MAXPRINTLEN 2048

/* window identifiers */
#define MAXWINDOWS 128
static int windows [MAXWINDOWS],
	   windowscount = 0;

static int main_menu = 0,
	   view_menu = 0,
	   export_menu = 0;

#define MAXVIEWPORTS 128
static struct
{
  int id;

  int x, y, w, h;

  double rx, ry, rw, rh; /* original ratios */

  View_Render render;

  int is3D;
}
viewports [MAXVIEWPORTS];
static int viewportscount = 0;

/* user callbacks */
static struct
{
  View_Init init;
  View_Idle idle;
  View_Quit quit;
  View_Render render;
  View_Key key;
  View_Keyspec keyspec;
  View_Mouse mouse;
  View_Motion motion;
  View_Passive_Motion passive;
  short holdmouse;
} user;

/* point of view */
static struct
{
  GLdouble from [3],
           to [3],
           up [3],
	   rhs [3],
	   left,
	   right,
	   bottom,
	   top,
	   neardst,
	   fardst;
} oldlook,
  look;

/* mouse modes */
static enum
{
  MOUSE_MODE_NONE,
  MOUSE_MODE_TRACKBALL,
  MOUSE_MODE_SPIN,
  MOUSE_MODE_ZOOM,
  MOUSE_MODE_MOVE
} activemode = MOUSE_MODE_TRACKBALL, mousemode;

/* reshape modes */
static enum 
{
  RESHAPE_MODE_ORTHO,
  RESHAPE_MODE_PERSPECTIVE
} reshapemode = RESHAPE_MODE_ORTHO;

/* menu items */
enum
{
  MENU_QUIT,
  MENU_TRACKBALL,
  MENU_SPIN,
  MENU_ZOOM,
  MENU_MOVE,
  MENU_VIEW_PERSPECTIVE,
  MENU_VIEW_ORTHO,
  MENU_VIEW_OUTLINE,
  MENU_VIEW_FILL,
  MENU_VIEW_FRONT,
  MENU_VIEW_BACK,
  MENU_VIEW_LEFT,
  MENU_VIEW_RIGHT,
  MENU_VIEW_TOP,
  MENU_VIEW_BOTTOM,
  MENU_EXPORT_AVI_START,
  MENU_EXPORT_AVI_STOP,
  MENU_EXPORT_PDF,
  MENU_EXPORT_EPS,
  MENU_EXPORT_BMP
};

/* menu shift caused by user menu */
static int menu_shift = 0;

/* mouse coordinates */
static int xstart, ystart;

/* main window dimensions */
static int width, height;

/* text input window */
#define TEXTLEN 512
static struct
{
  int id;
  char *title;
  void (*done) (char*); /* callback called when the text has been read */
  int length;
  char text [TEXTLEN+1];
  short cursor;
  short visible;
} input;

/* AVI movie context */
void *AVI = NULL;

/* rectangle drawing */
static struct
{
  int enabled,
      x1, y1,
      x2, y2;

} rectangle = {0, 0, 0, 0, 0};

/* render rectangle */
static void render_rectangle ()
{
  glDisable (GL_DEPTH_TEST);
  glMatrixMode (GL_PROJECTION);
  glPushMatrix ();
  glLoadIdentity ();
  gluOrtho2D (0, width, 0, height);
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix ();
  glLoadIdentity ();
  glColor4d (0.0, 1.0, 0.0, 0.2);
  glEnable (GL_BLEND);
  glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glRecti (rectangle.x1, rectangle.y1,
	   rectangle.x2, rectangle.y2);
  glDisable (GL_BLEND);
  glBegin (GL_LINE_LOOP);
  glVertex2d (rectangle.x1+1, rectangle.y1+1);
  glVertex2d (rectangle.x1+1, rectangle.y2-1);
  glVertex2d (rectangle.x2-1, rectangle.y2-1);
  glVertex2d (rectangle.x2-1, rectangle.y1+1);
  glEnd ();
  glPopMatrix ();
  glMatrixMode (GL_PROJECTION);
  glPopMatrix ();
  glEnable (GL_DEPTH_TEST);
}

/* global setup */
static void globals3D ()
{
  GLfloat pos [4] =
   {look.from [0],
    look.from [1],
    look.from [2], 1.0},
    specular [4] =
    {0.0, 0.0, 0.0, 1.0},
    emission [4] =
    {0.0, 0.0, 0.0, 1.0};

  glLightfv (GL_LIGHT0, GL_POSITION, pos);
  glEnable (GL_LIGHTING);
  glEnable (GL_LIGHT0);
  glEnable (GL_COLOR_MATERIAL);
  glShadeModel (GL_SMOOTH);
  glFrontFace (GL_CCW);
  glEnable (GL_DEPTH_TEST);
  glEnable (GL_NORMALIZE);
  glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, specular);
  glMaterialfv (GL_FRONT_AND_BACK, GL_EMISSION, emission);
  glEnable (GL_CULL_FACE);
  glCullFace (GL_BACK);
}

/* basic viewport reshape */
static void basic_reshape2D (int w, int h)
{
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluOrtho2D (0., (double) w, 0., (double) h);
  glMatrixMode (GL_MODELVIEW);
}

/* reshape 2D subwindow */
static void reshape2D (int w, int h)
{
  glViewport (0, 0, w, h);
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  gluOrtho2D (0., (double) w, 0., (double) h);
  glMatrixMode (GL_MODELVIEW);
}

/* basic reshape 3D */
static void basic_reshape3D (int w, int h)
{
  double from [3], to [3], up [3],
	 left, right, bottom, top,
	 neardst, fardst, d;

  SUB (look.from, look.to, from);
  NORMALIZE (from);
  SCALE (from, 2);
  SET (to, 0);
  COPY (look.up, up);
  d = MAX (MAX (look.left, look.right),
           MAX (look.top, look.bottom));
  left = look.left / d;
  right = look.right / d;
  bottom = look.bottom / d;
  top = look.top / d;
  neardst = look.neardst;
  fardst = look.fardst;

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();

  switch (reshapemode)
  {
    case RESHAPE_MODE_ORTHO:
    {
      if (w <= h)
	glOrtho (left, right,
	  bottom * (GLdouble)h/(GLdouble)w,
	  top * (GLdouble)h/(GLdouble)w,
	  neardst, fardst);
      else
	glOrtho (left * (GLdouble)w/(GLdouble)h,
	  right * (GLdouble)w/(GLdouble)h,
	  bottom, top, neardst, fardst);
    }
    break;
    case RESHAPE_MODE_PERSPECTIVE:
    {
      gluPerspective (70.,
	(double)w/(double)h,
	neardst,
	fardst);
    }
    break;
  }

  gluLookAt (from [0], from [1], from [2],
             to [0], to [1], to [2],
	     up [0], up [1], up [2]);
 
  globals3D ();
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();
}

/* reshape the main window */
static void reshape3D (int w, int h)
{
  if (h == 0) h = 1;
  if (w == 0) w = 1;
  width = w;
  height = h;
  glViewport (0, 0, w, h);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();

  switch (reshapemode)
  {
    case RESHAPE_MODE_ORTHO:
    {
      if (w <= h)
	glOrtho (look.left, look.right,
	  look.bottom * (GLdouble)h/(GLdouble)w,
	  look.top * (GLdouble)h/(GLdouble)w,
	  look.neardst, look.fardst);
      else
	glOrtho (look.left * (GLdouble)w/(GLdouble)h,
	  look.right * (GLdouble)w/(GLdouble)h,
	  look.bottom, look.top,
	  look.neardst, look.fardst);
    }
    break;
    case RESHAPE_MODE_PERSPECTIVE:
    {
      gluPerspective (70.,
	(double)w/(double)h,
	look.neardst,
	look.fardst);
    }
    break;
  }

  gluLookAt (
    look.from [0], look.from [1], look.from [2],
    look.to [0], look.to [1], look.to [2],
    look.up [0], look.up [1], look.up [2]);
 
  globals3D ();
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();

  for (int i = 0; i < viewportscount; i ++)
  {
    if (viewports [i].w < 0) viewports [i].w = - (int)(viewports [i].rw * (double) w);

    if (viewports [i].h < 0) viewports [i].h = - (int)(viewports [i].rh * (double) h);

    if (viewports [i].x < 0)
    {
      if (viewports [i].rx < 0.8) viewports [i].x = - (int)(viewports [i].rx * (double) w);
      else viewports [i].x = - (w - ABS (viewports [i].w));
    }

    if (viewports [i].y < 0)
    {
      if (viewports [i].ry < 0.8) viewports [i].y = - (int)(viewports [i].ry * (double) h);
      else viewports [i].y = - (h - ABS (viewports [i].h));
    }
  }
}

/* update all windows */
static void updateall ()
{
  for (int w = windowscount - 1; w >= 0; w --) 
  {
    glutSetWindow (windows [w]);
    glutPostRedisplay();
  }

  reshape3D (width, height);
}

/* draw input window */
static void render2D ()
{
  GLint mode [2];
  glGetIntegerv (GL_POLYGON_MODE, mode);
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  glDisable (GL_LIGHTING);
  glDisable (GL_DEPTH_TEST);

  glColor3d (1, 1, 1);
  glRecti (0, 0, width, 24);
  glColor3d (0, 0, 0);
  glBegin (GL_LINE_LOOP);
  glVertex2d (2, 2);
  glVertex2d (width - 2, 2);
  glVertex2d (width - 2, 22);
  glVertex2d (2, 22);
  glEnd ();

  if (input.cursor) GLV_Print (8, 8, 0, GLV_FONT_10, "%s: %s|", input.title, input.text);
  else GLV_Print (8, 8, 0, GLV_FONT_10, "%s: %s", input.title, input.text);

  glEnable (GL_LIGHTING);
  glEnable (GL_DEPTH_TEST);
  glPolygonMode (GL_FRONT_AND_BACK, mode [0]);
}

/* draw main window */
static void render3D ()
{
  GLint viewport [4];

  glClearColor (1, 1, 1, 1);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glGetIntegerv (GL_VIEWPORT, viewport);
  gl2psBeginViewport (viewport);

  if (user.render) user.render ();

  gl2psEndViewport ();

  if (rectangle.enabled) render_rectangle ();

  for (int i = 0; i < viewportscount; i ++)
  {
    glViewport (ABS (viewports [i].x), ABS (viewports [i].y),
	        ABS (viewports [i].w), ABS (viewports [i].h));

    if (viewports [i].is3D) basic_reshape3D (ABS (viewports [i].w), ABS (viewports [i].h));
    else basic_reshape2D (ABS (viewports [i].w), ABS (viewports [i].h));

    glGetIntegerv (GL_VIEWPORT, viewport);
    gl2psBeginViewport (viewport);

    viewports [i].render ();

    gl2psEndViewport ();
  }

  glutSwapBuffers();

  if (AVI)
  {
    char *rgb = RGB_Alloc (width, height);
    glReadPixels (0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, rgb);
    AVI_Frame (AVI, rgb);
    RGB_Free (rgb);
  }

  if (viewportscount)
    reshape3D (width, height); /* restore main viewport projection */
}

/* mouse button events */
static void mouse3D (int button, int state, int x, int y)
{
  if (user.holdmouse) goto callback;

  switch (state)
  {
    case GLUT_UP:
    break;

    case GLUT_DOWN:
    xstart = x;
    ystart = y;
    break;
  }

  switch (button)
  {
    case GLUT_LEFT_BUTTON:
    switch (state)
    {
    case GLUT_UP:
      mousemode = MOUSE_MODE_NONE;
      break;
    case GLUT_DOWN:
      mousemode = activemode;
      memcpy (&oldlook, &look, sizeof (look));
      break;
    }
    break;

    case GLUT_MIDDLE_BUTTON:
    switch (state)
    {
    case GLUT_UP:
      break;
    case GLUT_DOWN:
      break;
    }
    break;

    case GLUT_RIGHT_BUTTON:
    switch (state)
    {
    case GLUT_UP:
      break;
    case GLUT_DOWN:
      break;
    }
    break;
  }

callback:
  if (user.mouse)
    user.mouse (button, state, x, y);
}

/* mouse motion */
static void motion3D (int x, int y)
{
  if (user.holdmouse) goto callback;

  switch (mousemode)
  {
    case MOUSE_MODE_TRACKBALL:
    {
      double R [9], a [3], b [3], e [3], r, u, s;

      r = (double)(x - xstart) / (double) MIN (width, height);
      u = -(double)(y - ystart) / (double) MIN (width, height); 
      COPY (oldlook.rhs, b);
      SCALE (b, r);
      ADDMUL (b, u, oldlook.up, b);
      ADD (b, oldlook.from, b);
      SUB (b, oldlook.to, b);
      SUB (oldlook.from, oldlook.to, a);
      PRODUCT (b, a, e);
      s = ALG_PI * sqrt (r*r+u*u) / LEN (e);
      SCALE (e, s);
      EXPMAP (e, R);
      NVMUL (R, a, b);
      ADD (look.to, b, look.from);
      NVMUL (R, oldlook.up, look.up);
      NVMUL (R, oldlook.rhs, look.rhs);
      updateall ();
    }
    break;
    case MOUSE_MODE_SPIN:
    {
      double R [9], a [3], b [3], e [3], u, s;

      u = -(double)(y - ystart) / (double) MIN (width, height); 
      SUB (oldlook.from, oldlook.to, a);
      COPY (a, e);
      s = ALG_PI * u / LEN (e);
      SCALE (e, s);
      EXPMAP (e, R);
      NVMUL (R, a, b);
      ADD (look.to, b, look.from);
      COPY (oldlook.up, a);
      NVMUL (R, oldlook.up, look.up);
      NVMUL (R, oldlook.rhs, look.rhs);
      updateall ();
    }
    break;
    case MOUSE_MODE_ZOOM:
    {
      double a [3], c [3], d [3], u, l, r, b, t;

      u = -(double)(y - ystart) / (double) MIN (width, height); 

      switch (reshapemode)
      {
        case RESHAPE_MODE_ORTHO:
	  l = oldlook.left - (oldlook.right - oldlook.left) * u;
	  r = oldlook.right + (oldlook.right - oldlook.left) * u;
	  b = oldlook.bottom - (oldlook.top - oldlook.bottom) * u;
	  t = oldlook.top + (oldlook.top - oldlook.bottom) * u;
	  if (l < r) { look.left = l; look.right = r; }
	  if (b < t) { look.bottom = b; look.top = t; }
        break;
        case RESHAPE_MODE_PERSPECTIVE:
	  SUB (oldlook.from, oldlook.to, a);
	  ADDMUL (oldlook.from, u, a, c);
	  SUB (c, look.to, d);
	  if (DOT (a, d) > 0.0)
	  {
	    COPY (c, look.from);
	  }
        break;
      }
      
      updateall ();
    }
    break;
    case MOUSE_MODE_MOVE:
    {
      double b [3], r, u;

      r = (oldlook.right - oldlook.left) * (double)(x - xstart) / (double) MIN (width, height);
      u = -(oldlook.top - oldlook.bottom) * (double)(y - ystart) / (double) MIN (width, height); 
      COPY (oldlook.rhs, b);
      SCALE (b, r);
      ADDMUL (b, u, oldlook.up, b);
      SUB (oldlook.from, b, look.from);
      SUB (oldlook.to, b, look.to);
      updateall ();
    }
    break;
    case MOUSE_MODE_NONE:
    break;
  }

callback:
  if (user.motion)
    user.motion (x, y);
}

/* passive mouse motion */
static void passive3D (int x, int y)
{
  if (user.passive)
    user.passive (x, y);
}
 
/* handle special keys */
static void keyspec3D (int key, int x, int y)
{
  if (user.keyspec)
    user.keyspec (key, x, y);
}

/* handle regular keys */
static void key3D (unsigned char key, int x, int y)
{
  /* input gets all focus */
  if (input.visible)
  {
    if (key == '\r')
    {
      input.visible = 0;
      GLV_Close_Viewport (input.id);
      return;
    }
    else if (key == 27)
    {
      input.visible = 0;
      GLV_Close_Viewport (input.id);
      input.length = 0;
      input.text [0] = '\0';
      return;
    }
    else if (key == '\t')
    {
      if ((input.length + 1) < TEXTLEN)
      {
	input.text [input.length ++] = ' ';
	input.text [input.length ++] = ' ';
	input.text [input.length] = '\0';
      }
    }
    else if (key == '\b' || key == 127) /* 127 is DEL (Mac) */
    {
      if (input.length)
        input.text [-- input.length] = '\0';
    }
    else if (input.length < TEXTLEN)
    {
      input.text [input.length] = key;
      input.text [++ input.length] = '\0';
    }

    updateall ();
  }
  
  else /* focus back to
	  the main window */
 
  if (user.key)
    user.key (key, x, y);
}

/* idle processing */
static void idle3D ()
{
  if (user.idle)
    if (user.idle ())
      updateall ();
}

/* point of view initialisation */
static void init3D (double *extents)
{
  double a [3], b [3], rl, tb, nf, mx;

  COPY (extents + 3, look.from);
  MID (extents, extents + 3, look.to);
  SUB (look.to, look. from, a);
  VECTOR (b, 0.0, 0.0, 1.0);
  PRODUCT (a, b, look.rhs);
  PRODUCT (look.rhs, a, look.up);
  NORMALIZE (look.rhs);
  NORMALIZE (look.up);

  rl = extents [3] - extents [0];
  tb = extents [4] - extents [1];
  nf = extents [5] - extents [2];
  mx = MAX (rl, MAX (tb, nf));
  look.left = -0.5*mx;
  look.right = 0.5*mx;
  look.bottom = -0.5*mx;
  look.top = 0.5*mx;
  look.fardst = 1000.0 * mx;
  look.neardst = 0.001 * MIN (rl, MIN (tb, nf));

  if (user.init)
    user.init ();
}

/* export PDF */
static void pdf (char *path)
{
  int state = GL2PS_OVERFLOW, buffsize = 0, len;
  GLint viewport [4] = {0, 0, width, height};
  FILE *fp;

  if (path)
  {
    len = strlen (path);
    if (strcmp (&path[len-4], ".pdf"))
      sprintf (&path [len], ".pdf");

    ASSERT (fp = fopen (path, "w"), ERR_FILE_OPEN);
   
    while (state == GL2PS_OVERFLOW)
    {
      buffsize += 16777216;
      gl2psBeginPage("Solfec picture", "Solfec", viewport, GL2PS_PDF, GL2PS_SIMPLE_SORT, 0, GL_RGBA, 0, NULL, 8, 8, 8, buffsize, fp, path);
      render3D ();
      state = gl2psEndPage();
    }

    fclose(fp);
  }
}

/* export EPS */
static void eps (char *path)
{
  int state = GL2PS_OVERFLOW, buffsize = 0, len;
  GLint viewport [4] = {0, 0, width, height};
  FILE *fp;

  if (path)
  {
    len = strlen (path);
    if (strcmp (&path[len-4], ".eps"))
      sprintf (&path [len], ".eps");

    ASSERT (fp = fopen (path, "w"), ERR_FILE_OPEN);
   
    while (state == GL2PS_OVERFLOW)
    {
      buffsize += 16777216;
      gl2psBeginPage("Solfec picture", "Solfec", viewport, GL2PS_EPS, GL2PS_SIMPLE_SORT, 0, GL_RGBA, 0, NULL, 8, 8, 8, buffsize, fp, path);
      render3D ();
      state = gl2psEndPage();
    }

    fclose(fp);
  }
}

/* export BMP */
static void bmp (char *path)
{
  int len;
  void *rgb;

  if (path)
  {
    len = strlen (path);
    if (strcmp (&path[len-4], ".bmp"))
      sprintf (&path [len], ".bmp");
   
    rgb = RGB_Alloc (width, height);
    updateall ();
    glReadPixels (0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, rgb);
    BMP_Output (width, height, rgb, path);
    RGB_Free (rgb);
  }
}

/* export AVI */
static void avi (char *path)
{
  int len;

  len = strlen (path);
  if (strcmp (&path[len-4], ".avi"))
  sprintf (&path [len], ".avi");

  if (AVI) AVI_Close (AVI);
  AVI = AVI_Open (width, height, 16, path);
}

/* view menu */
static void menu_view3D (int value)
{
  switch (value)
  {
    case MENU_VIEW_PERSPECTIVE:
      glutSetMenu (view_menu);
      glutChangeToMenuEntry (1, "ortho", MENU_VIEW_ORTHO);
      reshapemode = RESHAPE_MODE_PERSPECTIVE;
      updateall ();
      break;
    case MENU_VIEW_ORTHO:
      glutSetMenu (view_menu);
      glutChangeToMenuEntry (1, "perspective", MENU_VIEW_PERSPECTIVE);
      reshapemode = RESHAPE_MODE_ORTHO;
      updateall ();
      break;
    case MENU_VIEW_OUTLINE:
      glutSetMenu (view_menu);
      glutChangeToMenuEntry (2, "fill", MENU_VIEW_FILL);
      glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
      updateall ();
      break;
    case MENU_VIEW_FILL:
      glutSetMenu (view_menu);
      glutChangeToMenuEntry (2, "outline", MENU_VIEW_OUTLINE);
      glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
      updateall ();
      break;
    case MENU_VIEW_FRONT:
      {
	double vec [3], len;
	SUB (look.from, look. to, vec);
	len = LEN (vec);
	COPY (look.to, look.from);
	look.from [1] += len;
	SET (look.up, 0);
	look.up [2] = 1;
	SUB (look.from, look.to, vec);
	PRODUCT (look.up, vec, look.rhs);
	NORMALIZE (look.rhs);
	NORMALIZE (look.up);
      }
      updateall ();
      break;
    case MENU_VIEW_BACK:
      {
	double vec [3], len;
	SUB (look.from, look. to, vec);
	len = LEN (vec);
	COPY (look.to, look.from);
	look.from [1] -= len;
	SET (look.up, 0);
	look.up [2] = 1;
	SUB (look.from, look.to, vec);
	PRODUCT (look.up, vec, look.rhs);
	NORMALIZE (look.rhs);
	NORMALIZE (look.up);
      }
      updateall ();
      break;
    case MENU_VIEW_LEFT:
      {
	double vec [3], len;
	SUB (look.from, look. to, vec);
	len = LEN (vec);
	COPY (look.to, look.from);
	look.from [0] -= len;
	SET (look.up, 0);
	look.up [2] = 1;
	SUB (look.from, look.to, vec);
	PRODUCT (look.up, vec, look.rhs);
	NORMALIZE (look.rhs);
	NORMALIZE (look.up);
      }
      updateall ();
      break;
    case MENU_VIEW_RIGHT:
      {
	double vec [3], len;
	SUB (look.from, look. to, vec);
	len = LEN (vec);
	COPY (look.to, look.from);
	look.from [0] += len;
	SET (look.up, 0);
	look.up [2] = 1;
	SUB (look.from, look.to, vec);
	PRODUCT (look.up, vec, look.rhs);
	NORMALIZE (look.rhs);
	NORMALIZE (look.up);
      }
      updateall ();
      break;
    case MENU_VIEW_TOP:
      {
	double vec [3], len;
	SUB (look.from, look. to, vec);
	len = LEN (vec);
	COPY (look.to, look.from);
	look.from [2] += len;
	SET (look.up, 0);
	look.up [1] = 1;
	SUB (look.from, look.to, vec);
	PRODUCT (look.up, vec, look.rhs);
	NORMALIZE (look.rhs);
	NORMALIZE (look.up);
      }
      updateall ();
      break;
    case MENU_VIEW_BOTTOM:
      {
	double vec [3], len;
	SUB (look.from, look. to, vec);
	len = LEN (vec);
	COPY (look.to, look.from);
	look.from [2] -= len;
	SET (look.up, 0);
	look.up [1] = 1;
	SUB (look.from, look.to, vec);
	PRODUCT (look.up, vec, look.rhs);
	NORMALIZE (look.rhs);
	NORMALIZE (look.up);
      }
      updateall ();
      break;
  }
}

/* export menu */
static void menu_export3D (int value)
{
  switch (value)
  {
    case MENU_EXPORT_AVI_START:
      glutSetMenu (export_menu);
      glutChangeToMenuEntry (1, "AVI STOP", MENU_EXPORT_AVI_STOP);
      GLV_Read_Text ("AVI FILE NAME", avi);
      break;
    case MENU_EXPORT_AVI_STOP:
      glutSetMenu (export_menu);
      glutChangeToMenuEntry (1, "AVI START", MENU_EXPORT_AVI_START);
      if (AVI) { AVI_Close (AVI); AVI = NULL; }
      break;
    case MENU_EXPORT_PDF:
      GLV_Read_Text ("PDF FILE NAME", pdf);
      break;
    case MENU_EXPORT_EPS:
      GLV_Read_Text ("EPS FILE NAME", eps);
      break;
    case MENU_EXPORT_BMP:
      GLV_Read_Text ("BMP FILE NAME", bmp);
      break;
  }
}

/* main menu */
static void menu3D (int value)
{
  switch (value)
  {
    case  MENU_TRACKBALL:
      activemode = MOUSE_MODE_TRACKBALL;
      break;
    case  MENU_SPIN:
      activemode = MOUSE_MODE_SPIN;
      break;
    case  MENU_ZOOM:
      activemode = MOUSE_MODE_ZOOM;
      break;
    case  MENU_MOVE:
      activemode = MOUSE_MODE_MOVE;
      break;
    case MENU_QUIT:
      if (AVI) AVI_Close (AVI);
      if (user.quit)
	user.quit ();
      exit (0);
      break;
  }
}

static void timer (int value)
{
  if (input.visible)
  {
    /* make the cursor blink */
    input.cursor = !input.cursor;
    updateall ();

    /* set up next trigger */
    glutTimerFunc (750, timer, 0);
  }
  else /* last timer run executes the 'done' callback => ... */
  {   /* ... this prevents capturing the input window when reading buffers */
    if (input.done)
    {
      input.done (input.text);
      input.done = NULL;
    }
  }
}

/* ------------ interface ------------ */

void GLV (
  int *argc,
  char **argv,
  char *title,
  int wdt,
  int hgh,
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
  View_Passive_Motion passive)
{
  int i, *codes;
  char **names;

  user.init = init;
  user.idle = idle;
  user.quit = quit;
  user.render = render;
  user.key = key;
  user.keyspec = keyspec;
  user.mouse = mouse;
  user.motion = motion;
  user.passive = passive;

  height = hgh;
  width = wdt;
  glutInit (argc, argv);
  glutInitWindowPosition (0, 0); 
  glutInitWindowSize (width, height);
  glutInitDisplayMode (GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE | GLUT_STENCIL);

  /* set up the main window */
  ASSERT (windows [windowscount ++] =
    glutCreateWindow (title), ERR_GLV_WINDOW);
  glutReshapeFunc (reshape3D);
  glutDisplayFunc (render3D);
  glutSpecialFunc (keyspec3D);
  glutKeyboardFunc (key3D);
  glutMouseFunc (mouse3D);
  glutMotionFunc (motion3D);
  glutPassiveMotionFunc (passive3D);
  glutIdleFunc (idle3D);

  view_menu = glutCreateMenu (menu_view3D);
  glutAddMenuEntry ("perspective", MENU_VIEW_PERSPECTIVE);
  glutAddMenuEntry ("outline", MENU_VIEW_OUTLINE);
  glutAddMenuEntry ("front", MENU_VIEW_FRONT);
  glutAddMenuEntry ("back", MENU_VIEW_BACK);
  glutAddMenuEntry ("left", MENU_VIEW_LEFT);
  glutAddMenuEntry ("right", MENU_VIEW_RIGHT);
  glutAddMenuEntry ("top", MENU_VIEW_TOP);
  glutAddMenuEntry ("bottom", MENU_VIEW_BOTTOM);

  export_menu = glutCreateMenu (menu_export3D);
  glutAddMenuEntry ("AVI START", MENU_EXPORT_AVI_START);
  glutAddMenuEntry ("PDF", MENU_EXPORT_PDF);
  glutAddMenuEntry ("EPS", MENU_EXPORT_EPS);
  glutAddMenuEntry ("BMP", MENU_EXPORT_BMP);

  if (menu) menu_shift = menu (&names, &codes) + 2;
  else menu_shift = 2;

  main_menu = glutCreateMenu (menu3D);
  for (i = 0; i < menu_shift - 2; i ++) glutAddSubMenu (names [i], codes [i]);

  glutAddSubMenu ("view", view_menu);
  glutAddSubMenu ("export", export_menu);

  glutAddMenuEntry ("trackball", MENU_TRACKBALL);
  glutAddMenuEntry ("spin", MENU_SPIN);
  glutAddMenuEntry ("zoom", MENU_ZOOM);
  glutAddMenuEntry ("move", MENU_MOVE);
  glutAddMenuEntry ("quit", MENU_QUIT);
  glutAttachMenu (GLUT_RIGHT_BUTTON);

  init3D (extents);
  glutMainLoop ();
}

void GLV_Redraw_All (void)
{ 
  updateall ();
}

void GLV_Reset_Extents (double *extents)
{
  View_Init init;

  init = user.init;
  user.init = NULL;
  init3D (extents); /* let it not call user init */
  user.init = init;
  updateall ();
}

void GLV_Update_Extents (double *extents)
{
  double a [3], b [3], rl, tb, nf, mx;

  /* align to the middle point of extents
   * but along the same direction */
  SUB (look.from, look.to, a);
  NORMALIZE (a);
  SUB (extents + 3, extents, b);
  mx = LEN (b);
  look.to [0] = 0.5 * (extents [0] + extents [3]);
  look.to [1] = 0.5 * (extents [1] + extents [4]);
  look.to [2] = 0.5 * (extents [2] + extents [5]);
  ADDMUL (look.to, mx, a, look.from);

  rl = extents [3] - extents [0];
  tb = extents [4] - extents [1];
  nf = extents [5] - extents [2];
  mx = MAX (rl, MAX (tb, nf));
  look.left = -0.5*mx;
  look.right = 0.5*mx;
  look.bottom = -0.5*mx;
  look.top = 0.5*mx;
  look.fardst = 1000.0 * mx;
  look.neardst = 0.001 * MIN (rl, MIN (tb, nf));

  updateall ();
}

double GLV_Minimal_Extent ()
{
  double a = look.right - look.left,
	 b = look.top - look.bottom;

  return MIN (a, b);
}

void GLV_Sizes (int *w, int *h)
{
  *w = width;
  *h = height;
}

int GLV_Open_Window (
  int x,
  int y,
  int w,
  int h,
  View_Render render)
{
  if (windowscount < MAXWINDOWS)
  {
    ASSERT (windows [windowscount ++] =
      glutCreateSubWindow (windows [0], x, y, w, h),
      ERR_GLV_WINDOW);

    glutDisplayFunc (render);
    glutReshapeFunc (reshape2D);
    return windows [windowscount -1];
  }

  return -1;
}

void GLV_Close_Window (int window)
{
  for (int n = 0; n < windowscount; n ++)
  {
    if (windows [n] == window)
    {
      for (; n + 1 < windowscount; n ++)
	windows [n] = windows [n + 1];
      glutDestroyWindow (window);
      windowscount --;
      updateall ();
      return;
    }
  }
}

int GLV_Open_Viewport (
  int x,
  int y,
  int w,
  int h,
  int is3D,
  View_Render render)
{
  if (viewportscount < MAXVIEWPORTS)
  {
    int i, j;

    srand (time (NULL));

    while (1) /* generate unique viewport id */
    {
      j = rand ();

      for (i = 0; i < viewportscount; i ++)
	if (viewports [i].id == j) break;

      if (i == viewportscount) break;
    }

    viewports [i].id = j;
    viewports [i].x = x;
    viewports [i].y = y;
    viewports [i].w = w;
    viewports [i].h = h;
    viewports [i].rx = (double) ABS(x) / (double) width;
    viewports [i].ry = (double) ABS(y) / (double) height;
    viewports [i].rw = (double) ABS(w) / (double) width;
    viewports [i].rh = (double) ABS(h) / (double) height;
    viewports [i].is3D = is3D;
    viewports [i].render = render;
    viewportscount ++;

    updateall ();

    return viewports [viewportscount-1].id;
  }

  return -1;
}

void GLV_Move_Viewport (int viewport, int x, int y, int w, int h)
{
  for (int n = 0; n < viewportscount; n ++)
  {
    if (viewports [n].id == viewport)
    {
      viewports [n].x = x;
      viewports [n].y = y;
      viewports [n].w = w;
      viewports [n].h = h;
      return;
    }
  }
}

void GLV_Resize_Viewport (int viewport, int w, int h)
{
  for (int n = 0; n < viewportscount; n ++)
  {
    if (viewports [n].id == viewport)
    {
      viewports [n].w = w;
      viewports [n].h = h;
      return;
    }
  }
}

void GLV_Close_Viewport (int viewport)
{
  for (int n = 0; n < viewportscount; n ++)
  {
    if (viewports [n].id == viewport)
    {
      for (; n + 1 < viewportscount; n ++)
	viewports [n] = viewports [n + 1];
      viewportscount --;
      updateall ();
      return;
    }
  }
}

void GLV_Read_Text (char *title, void (*done) (char *text))
{
  input.id = GLV_Open_Viewport (0, -(height / 2 - 6), -width, 24, 0, render2D);
  input.title = title;
  input.done = done;
  input.cursor = 0;
  input.length = 0;
  input.text [0] = '\0';
  input.visible = 1;
  glutTimerFunc (750, timer, 0);
}

void GLV_Print (double x, double y, double z, int font, char *fmt, ...)
{
  va_list arg;
  char buff [MAXPRINTLEN];
  int i;

  va_start (arg, fmt);
  vsnprintf (buff, MAXPRINTLEN, fmt, arg); /* read formated string */
  va_end (arg);

  glRasterPos3d (x, y, z); /* string position */

  switch (font)
  {
    case GLV_FONT_8_BY_13:
      gl2psText (buff, "Courier", 13);
      for (i = 0; buff[i]; i++)
	glutBitmapCharacter (GLUT_BITMAP_8_BY_13, buff[i]);
      break;
    case GLV_FONT_9_BY_15:
      gl2psText (buff, "Courier", 15);
      for (i = 0; buff[i]; i++)
	glutBitmapCharacter (GLUT_BITMAP_9_BY_15, buff[i]);
      break;
    case GLV_FONT_10:
      gl2psText (buff, "Courier", 10);
      for (i = 0; buff[i]; i++)
	glutBitmapCharacter (GLUT_BITMAP_HELVETICA_10, buff[i]); /* output characters */
      break;
    case GLV_FONT_12:
      gl2psText (buff, "Courier", 12);
      for (i = 0; buff[i]; i++)
	glutBitmapCharacter (GLUT_BITMAP_HELVETICA_12, buff[i]);
      break;
    case GLV_FONT_18:
      gl2psText (buff, "Courier", 18);
      for (i = 0; buff[i]; i++)
	glutBitmapCharacter (GLUT_BITMAP_HELVETICA_18, buff[i]);
      break;
  }
}

int GLV_Print_Width (int font, char *fmt, ...)
{
  va_list arg;
  char buff [MAXPRINTLEN];

  va_start (arg, fmt);
  vsnprintf (buff, MAXPRINTLEN, fmt, arg); /* read formated string */
  va_end (arg);

  return font * strlen (buff);
}

void GLV_Screen_Bitmap (char *path)
{ bmp (path); }

void GLV_Hold_Mouse ()
{ user.holdmouse = 1; }

void GLV_Release_Mouse ()
{ user.holdmouse = 0; }

void GLV_SetProjectionMatrix (int w, int h)
{
  switch (reshapemode)
  {
    case RESHAPE_MODE_ORTHO:
    {
      if (w <= h)
	glOrtho (look.left, look.right,
	  look.bottom * (GLdouble)h/(GLdouble)w,
	  look.top * (GLdouble)h/(GLdouble)w,
	  look.neardst, look.fardst);
      else
	glOrtho (look.left * (GLdouble)w/(GLdouble)h,
	  look.right * (GLdouble)w/(GLdouble)h,
	  look.bottom, look.top,
	  look.neardst, look.fardst);
    }
    break;
    case RESHAPE_MODE_PERSPECTIVE:
    {
      gluPerspective (70.,
	(double)w/(double)h,
	look.neardst,
	look.fardst);
    }
    break;
  }

  gluLookAt (
    look.from [0], look.from [1], look.from [2],
    look.to [0], look.to [1], look.to [2],
    look.up [0], look.up [1], look.up [2]);
}

void GLV_Rectangle_On (int x1, int y1, int x2, int y2)
{
  rectangle.enabled = 1;
  rectangle.x1 = MIN (x1, x2);
  rectangle.y1 = MIN (height-y1, height-y2);
  rectangle.x2 = MAX (x1, x2);
  rectangle.y2 = MAX (height-y1, height-y2);
  updateall ();
}

void GLV_Rectangle_Off ()
{
  rectangle.enabled = 0;
}
