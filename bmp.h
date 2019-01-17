/*
 * bmp.h
 * Copyright (C) 2007, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * minimalistic bitmap and avi output
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

#ifndef __bmp__
#define __bmp__

#if __cplusplus
extern "C" { /* C */
#endif

/* allocate RGB buffer */
void* RGB_Alloc (int width, int height);

/* free the buffer */
void RGB_Free (void *rgb);

/* output an RGB bitmap file (24 bits per pixel) */
void BMP_Output (int width, int height, void *rgb, const char *path);

/* return expected duration */
double AVI_Duration (int frames, int fps);

/* open an RGB avi file (24 bits per pixel, 'fps' frames per second) */
void* AVI_Open (int width, int height, int fps, const char *path);

/* output one frame  */
void AVI_Frame (void *avi, void *rgb);

/* close the avi file */
void AVI_Close (void *avi);

#if __cplusplus
} /* extern C */
#endif

#endif
