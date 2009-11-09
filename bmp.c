/*
 * bmp.c
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

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "err.h"
#include "bmp.h"

/* offsets to AVI header entries updated
 * at the end of the write process */
#define FILE_SIZE_OFFSET 4
#define FRAMES_1_OFFSET 48
#define FRAMES_2_OFFSET 140
#define MOVIE_SIZE_OFFSET 216

/* movie
 * handle */
typedef
struct
{
  FILE *file;
  int width;
  int height;  
  int frames;
} AVI;

/* void pointer casting */
#define HANDLE(ptr) ((AVI*)ptr)

/* round up to the next 4-divisible width */
int ROUNDED_WIDTH (int width)
{
  width *= 3;
  while (width % 4) width ++;
  return width;
}

/* rounded up image size */
#define ROUNDED_SIZE(width, height)\
  (ROUNDED_WIDTH (width) * (height))

/* little-endian DWORD output */
void DWORD (int32_t v, FILE *f)
{
  putc (v, f);
  putc (v >> 8, f);
  putc (v >> 16, f);
  putc (v >> 24, f);
}

/* little-endian WORD output */
void WORD (int16_t v, FILE *f)
{
  putc (v, f);
  putc (v >> 8, f);
}

/* string output (no byte reversal) */
#define STRING(s, f) fwrite (s, 1, strlen (s), f)

/* little-endian RGB data output */
void RGB (int w, int h, char *rgb, FILE *f)
{
  w = ROUNDED_WIDTH (w);
  for (; h; h --)
  {
    for(char *b = rgb; b < rgb + w; b += 3)
    {
      char tmp = b [2];
      b [2] = b [0];
      b [0] = tmp;
    }

    fwrite (rgb, 1, w, f);
    rgb += w;
  }
}

/*----- interface -----*/

void* RGB_Alloc (int width, int height)
{
  void *rgb;
  ERRMEM (rgb = malloc (ROUNDED_SIZE (width, height)));
  return rgb;
}

void RGB_Free (void *rgb)
{ free (rgb); }

void BMP_Output (int width, int height, void *rgb, const char *path)
{
  FILE *f;

  ASSERT (f = fopen (path, "wb"), ERR_FILE_OPEN);

  WORD (0x4D42, f); /* BM magic */
  DWORD (ROUNDED_SIZE (width, height) + 54, f); /* header size */
  WORD (0, f); /* reserved */
  WORD (0, f); /* reserved */
  DWORD (54, f); /* offset to data */
  DWORD (40, f); /* DIB header size */
  DWORD (width, f);
  DWORD (height, f);
  WORD (1, f); /* number of planes */
  WORD (24, f); /* bits per pixel */
  DWORD (0, f); /* compression */
  DWORD (0, f); /* image size (will be calculated) */
  DWORD (0, f); /* horisontal pixels per meter */
  DWORD (0, f); /* vertical ... */
  DWORD (0, f); /* number of colors */
  DWORD (0, f); /* number of important colors */
  RGB (width, height, rgb, f);
  fclose (f);
}

double AVI_Duration (int frames, int fps)
{  return (double) frames / (double) fps; }

void* AVI_Open (int width, int height, int fps, const char *path)
{
  AVI *avi;
  FILE *f;

  ERRMEM (avi = calloc (sizeof (AVI), 1));
  ASSERT (avi->file = fopen (path, "wb"), ERR_FILE_OPEN);

  avi->width = width;
  avi->height = height;
  avi->frames = 0;
  f = avi->file;
  
  STRING ("RIFF", f);
  DWORD (0, f); /* file size => updated at the end */
  STRING ("AVI ", f);
  STRING ("LIST", f);
  DWORD (192, f); /* length of this list */
  STRING ("hdrl", f); /* headers list */
  STRING ("avih", f); /* avi header */
  DWORD (56, f); /* size of this header */
  DWORD (1000000.0 / (double)fps, f);
  DWORD (0, f); /* maximal bytes per second (leave default) */
  DWORD (0, f); /* reserved value */
  DWORD (16, f); /* file will have an index */
  DWORD (0, f); /* number of frames => updated at the end */
  DWORD (0, f); /* initial number of frames */
  DWORD (1, f); /* only one video stram */
  DWORD (ROUNDED_SIZE (width, height), f);
  DWORD (width, f);
  DWORD (height, f);
  DWORD (0, f); /* scale */
  DWORD (0, f); /* rate */
  DWORD (0, f); /* start */
  DWORD (0, f); /* length */
  STRING ("LIST", f);
  DWORD (116, f); /* size of this list */
  STRING ("strl", f); /* video stream list */
  STRING ("strh", f);  /* video stream header */
  DWORD (56, f); /* size of this header */
  STRING ("vids", f);
  STRING ("DIB ", f); /* Device Independent Bitmap */
  DWORD (0, f); /* flags */
  DWORD (0, f); /* priority */
  DWORD (0, f); /* initial frames */
  DWORD (1000000.0 / (double)fps, f); /* micro seconds per frame */
  DWORD (1000000, f); /* rate */
  DWORD (0, f); /* start */
  DWORD (0, f); /* number of frames => updated at the end */
  DWORD (ROUNDED_SIZE (width, height), f);
  DWORD (0, f); /* quality */
  DWORD (0, f); /* sample size */
  DWORD (0, f); /* reserved */
  WORD (width, f);
  WORD (height, f);
  STRING ("strf", f); /* video stream format */
  DWORD (40, f); /* size of header */
  DWORD (40, f); /* size of DIB header */
  DWORD (width, f);
  DWORD (height, f);
  WORD (1, f); /* number of planes */
  WORD (24, f); /* bits per pixel */
  DWORD (0, f); /* compression */
  DWORD (ROUNDED_SIZE (width, height), f); /* image size */
  DWORD (0, f); /* pixels per meter along x */
  DWORD (0, f); /* ... along y */
  DWORD (0, f); /* number of colors */
  DWORD (0, f); /* number of important colors */
  STRING ("LIST", f);
  DWORD (0, f); /* movie size => updated at the end */
  STRING ("movi", f); /* the movie frames begin here */
  return avi;
}

void AVI_Frame (void *avi, void *rgb)
{
  STRING ("00db", HANDLE (avi)->file);
  DWORD (ROUNDED_SIZE (HANDLE (avi)->width,
    HANDLE (avi)->height), HANDLE (avi)->file);
  RGB (HANDLE (avi)->width, HANDLE (avi)->height,
    rgb, HANDLE (avi)->file);
  
  HANDLE (avi)->frames ++;
}

void AVI_Close (void *avi)
{
  int offset = 4, frame = HANDLE (avi)->frames, index_size = frame * 16,
    frame_size = ROUNDED_SIZE (HANDLE (avi)->width, HANDLE (avi)->height),
    movie_size = 4 + (frame * (frame_size + 8)),
    file_size = index_size + movie_size + 212;

  STRING ("idx1", HANDLE (avi)->file);
  DWORD (index_size, HANDLE (avi)->file);

  for (frame = 0; frame < HANDLE (avi)->frames; frame ++)
  {
    STRING ("00db", HANDLE (avi)->file);
    DWORD (16, HANDLE (avi)->file);
    DWORD (offset, HANDLE (avi)->file);
    DWORD (frame_size, HANDLE (avi)->file);
    offset += frame_size + 8;
  }
  
  fseek (HANDLE (avi)->file, FILE_SIZE_OFFSET, SEEK_SET);
  DWORD (file_size, HANDLE (avi)->file);
  fseek (HANDLE (avi)->file, FRAMES_1_OFFSET, SEEK_SET);
  DWORD (HANDLE (avi)->frames, HANDLE (avi)->file);
  fseek (HANDLE (avi)->file, FRAMES_2_OFFSET, SEEK_SET);
  DWORD (HANDLE (avi)->frames, HANDLE (avi)->file);
  fseek (HANDLE (avi)->file, MOVIE_SIZE_OFFSET, SEEK_SET);
  DWORD (movie_size, HANDLE (avi)->file);
  fclose (HANDLE (avi)->file);
  free (avi);
}
