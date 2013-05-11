/*
 * pbf.h
 * Copyright (C) 2006, Tomasz Koziara (t.koziara AT gmail.com)
 * --------------------------------------------------------------
 * portable binary format (PBF)
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

#ifndef __pbf__
#define __pbf__

#if HDF5 /* FIXME new implementation FIXME */

#include <hdf5.h>
#include <hdf5_hl.h>

#define PBF_MAXSTACK 128 /* maximal group stack */

typedef struct pbf PBF; /* file type */

/* access mode */
typedef enum {PBF_READ, PBF_WRITE} PBF_ACC;

/* compression flag */
typedef enum {PBF_ON, PBF_OFF} PBF_FLG;

/* file */
struct pbf
{
  PBF_ACC mode; /* access mode */
  PBF_FLG compression; /* compression flag */
  PBF_FLG parallel; /* parallel flag */

  double time, *times; /* current time, all times (READ) */

  int frame, count; /* current time frame number, total frames count */

  int *i; /* raw ints space */
  int ipos, ints; /* raw ints position and size */
  double *d; /* raw doubles space */
  int dpos, doubles; /* raw doubles position and size */

  hid_t stack [PBF_MAXSTACK]; /* file id followed by groups stack */
  short top; /* index of the stack top item */

  PBF *next; /* list of parallel files (READ) */
};

/* open for writing */
PBF* PBF_Write (const char *path, PBF_FLG append, PBF_FLG parallel);

/* open for reading */
PBF* PBF_Read (const char *path);

/* close file */
void PBF_Close (PBF *bf);

/* initialize new time frame for writing/reading */
void PBF_Time (PBF *bf, double *time);

/* set current label; in read mode return
 * positive value if the label was found */
int PBF_Label (PBF *bf, const char *label);

/* write/read raw data <=> order of writing and reading must be preserved */
void PBF_Short (PBF *bf, short *value, int length);
void PBF_Int (PBF *bf, int *value, int length);
void PBF_Uint (PBF *bf, unsigned int *value, int length);
void PBF_Double (PBF *bf, double *value, int length);

/* read/write NULL-termined string */
void PBF_String (PBF *bf, char **value);

/* has group check */
int PBF_Has_Group (PBF *bf, const char *name);

/* push/pop group on stak */
void PBF_Push (PBF *bf, const char *name);
void PBF_Pop (PBF *bf);

/* write/read named datasets (length > 1) or attributes (length == 1) */
void PBF_Int2 (PBF *bf, const char *name, int *value, hsize_t length);
void PBF_Double2 (PBF *bf, const char *name, double *value, hsize_t length);
void PBF_String2 (PBF *bf, const char *name, char **value);

/* get time limits in read mode */
void PBF_Limits (PBF *bf, double *start, double *end);

/* seek to time in read mode */
void PBF_Seek (PBF *bf, double time);

/* make 'steps' backward in read mode */
int PBF_Backward (PBF *bf, int steps);

/* make 'steps' forward in read mode */
int PBF_Forward (PBF *bf, int steps);

/* get number of time instants spanned by [t0, t1] */
unsigned int PBF_Span (PBF *bf, double t0, double t1);

#else /* old XDR based implementation */

#include <stdint.h>
#include <stdio.h>
#include <rpc/types.h>
#include <rpc/xdr.h>
#include "map.h"
#include "mem.h"

#if __MINGW32__
  #define FSEEK fseeko64
  #define FTELL ftello64
  #define OFF_T off64_t
#elif OSTYPE_LINUX
  #define FSEEK fseeko64
  #define FTELL ftello64
  #define OFF_T __off64_t
#else
  #define FSEEK fseeko
  #define FTELL ftello
  #define OFF_T off_t
#endif

#if __APPLE__
  #define xdr_uint64_t xdr_u_int64_t
#endif

#define PBF_MAXSTRING 4096 /* maximal string length used */

typedef struct pbf_marker PBF_MARKER; /* file marker */
typedef struct pbf_label PBF_LABEL; /* label type */
typedef struct pbf PBF; /* file type */

/* marker */
struct pbf_marker
{
  double time; /* time moment */
  u_int ipos; /* index position */
  uint64_t doff; /* data offset */
};

/* label */
struct pbf_label
{
  char *name; /* label name */
  int index; /* unique index */
  u_int dpos; /* data position */
};

/* access mode */
typedef enum {PBF_READ, PBF_WRITE} PBF_ACC;

/* compression flag */
typedef enum {PBF_ON, PBF_OFF} PBF_FLG;

/* file */
struct pbf
{
  PBF_ACC mode; /* access mode */
  char *dph, /* data path */
       *iph, /* index path */
       *lph; /* label path */
  FILE *dat, /* data file */
       *idx, /* index file */
       *lab; /* label file */
  XDR x_dat, /* data stream */
      x_idx, /* index stream */
      x_lab; /* labels stream */
  char *mem; /* read/write memory */
  u_int membase, /* memory base */
	memsize; /* memory size */
  MEM mappool, /* map items pool */
      labpool; /* labels pool */
  PBF_LABEL *ltab; /* table of labels */
  MAP *labels; /* name mapped labels */
  PBF_MARKER *mtab; /* markers */
  int lsize; /* free index (WRITE) or ltab size (READ) */
  u_int msize, /* mtab size (READ) */
        cur; /* index of current time frame */
  double time; /* current time */
  PBF_FLG compression; /* compression flag */
  PBF_FLG parallel; /* parallel flag */
  PBF *next; /* list of parallel files (READ) */
};

/* open for writing */
PBF* PBF_Write (const char *path, PBF_FLG append, PBF_FLG parallel);

/* open for reading */
PBF* PBF_Read (const char *path);

/* close file */
void PBF_Close (PBF *bf);

/* read/write current time */
void PBF_Time (PBF *bf, double *time);

/* set current label; in read mode return
 * positive value if the label was found */
int PBF_Label (PBF *bf, const char *label);

/* read/write characters */
void PBF_Char (PBF *bf, char *value, unsigned int length);

/* read/write unsigned characters */
void PBF_Uchar (PBF *bf, unsigned char *value, unsigned int length);

/* read/write other simple types ... */
void PBF_Short (PBF *bf, short *value, unsigned int length);
void PBF_Ushort (PBF *bf, unsigned short *value, unsigned int length);
void PBF_Int (PBF *bf, int *value, unsigned int length);
void PBF_Uint (PBF *bf, unsigned int *value, unsigned int length);
void PBF_Long (PBF *bf, long *value, unsigned int length);
void PBF_Ulong (PBF *bf, unsigned long *value, unsigned int length);
void PBF_Float (PBF *bf, float *value, unsigned int length);
void PBF_Double (PBF *bf, double *value, unsigned int length);

/* read/write NULL-termined string */
void PBF_String (PBF *bf, char **value);

/* get time limits in read mode */
void PBF_Limits (PBF *bf, double *start, double *end);

/* seek to time in read mode */
void PBF_Seek (PBF *bf, double time);

/* make 'steps' backward in read mode */
int PBF_Backward (PBF *bf, unsigned int steps);

/* make 'steps' forward in read mode */
int PBF_Forward (PBF *bf, unsigned int steps);

/* get number of time instants spanned by [t0, t1] */
unsigned int PBF_Span (PBF *bf, double t0, double t1);
#endif
#endif
