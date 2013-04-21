/*
 * pbf.c
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

#if MPI
#include <mpi.h>
#endif

#include <string.h>
#include <limits.h>
#include <float.h>
#include "ext/fastlz.h"
#include "pbf.h"
#include "err.h"
#include "alg.h"

/* memory increment */
#define CHUNK 1024

/* memory margin */
#define MARGIN 64

/* DAT file format:
 * ----------------
 *  [FRAME_0]
 *  [FRAME_1]
 *  ...
 *  [FRAME_N]
 * ----------
 *  FRAME_i:
 * -------------------
 *    [UNLABELED_DATA]
 *    [LABELED_0]
 *    [LABELED_1]
 *    ...
 *    [LABELED_K]
 * --------------
 */

/* IDX file format:
 * ----------------
 *  [FRAME_0]
 *  [FRAME_1]
 *  ...
 *  [FRAME_N]
 *  [FRAME_INF]
 * ----------
 *  FRAME_i:
 * -------------------------------
 *  [TIME] (double) {current time}
 *  [DOFF] (uint64_t) {offest of FRAME_i in DAT file}
 *  [IDX_0] (int) {index of first label}
 *  [POS_0] (u_int) {relative position of labeled data in XDR FRAME_i stream}
 *  [IDX_1]
 *  [POS_1]
 *  ...
 *  [IDX_K]
 *  [POS_K]
 *  [-1] (int) {end of labels marker}
 * ----------------------------------
 *  FRAME_INF:
 * -------------------------------
 *  [TIME] (double) {DBL_MAX time}
 *  [DOFF] (uint64_t) {offest to last data in DAT file}
 *  [-2] (int) {INF frame marker}
 * ----------------------------------
 */

/* LAB file format:
 * ----------------
 *  [LABEL_0] (string) {first label name}
 *  [LABEL_1]
 *  ...
 *  [LABEL_M]
 * ----------
 */

/* write to data file */
static u_int fileread (char **mem, u_int size, FILE *f)
{
  char cmp;

  fread (&cmp, 1, 1, f); /* read 1 byte */

  size --; /* subtract 1 byte */

  if (cmp)
  {
    int maxout, outsize;
    char *inp;

    ERRMEM (inp = malloc (size));
    ASSERT (fread (inp, 1, size, f) == size, ERR_PBF_READ);
    maxout = 2 * size;
    do
    {
      maxout *= 2;
      free (*mem);
      ERRMEM (*mem = malloc (maxout));
      outsize = fastlz_decompress (inp, size, *mem, maxout);
    } while (outsize == 0);
    size = outsize;
    free (inp);
  }
  else
  {
    free (*mem);
    ERRMEM (*mem = malloc (size));
    ASSERT (fread (*mem, 1, size, f) == size, ERR_PBF_READ);
  }

  return size;
}

/* read from data file */
static void filewrite (char *mem, u_int size, FILE *f, char cmp)
{
  if (cmp && size < 16) cmp = 0;  /* see ext/fastlz.h */

  fwrite (&cmp, 1, 1, f); /* write compresion flag (adds 1 byte per frame) */

  if (cmp)
  {
    int nbuf, num;
    char *out;

    nbuf = 2 * MAX (size, 66); /* see ext/fastlz.h */
    ERRMEM (out = malloc (nbuf));
    num = fastlz_compress (mem, size, out);
    WARNING (num < (int)size, "Compression increased the buffer size => Consider disabling it.");
    ASSERT_TEXT (num < nbuf, "Compression increased the buffer size by more than 100%%.");
    ASSERT (fwrite (out, 1, num, f) == (unsigned)num, ERR_PBF_WRITE);
    free (out);
  }
  else ASSERT (fwrite (mem, 1, size, f) == size, ERR_PBF_WRITE);
}

/* grow memory buffer in WRITE mode */
static void growmem (PBF *bf, u_int nb)
{
  /* increment memory base */
  bf->membase += xdr_getpos (&bf->x_dat);

  /* resize buffer */
  bf->memsize += 2 * bf->memsize + nb;
  ERRMEM (bf->mem = realloc (bf->mem, bf->memsize));

  /* create new bigger stream starting at base */
  xdr_destroy (&bf->x_dat);
  xdrmem_create (&bf->x_dat, bf->mem + bf->membase, bf->memsize - bf->membase, XDR_ENCODE);
}

/* initialize a frame to be red */
static void initialise_frame (PBF *bf, int frm)
{
  PBF_LABEL *l;
  int index;

  bf->cur = frm; /* set current frame */
  bf->time = bf->mtab [frm].time; /* and time */

  /* create new memory XDR stream for DATA chunk */
  bf->memsize = bf->mtab [frm+1].doff - bf->mtab [frm].doff;
  FSEEK (bf->dat, (OFF_T) bf->mtab [frm].doff, SEEK_SET);
  bf->memsize = fileread (&bf->mem, bf->memsize, bf->dat);
  xdr_destroy (&bf->x_dat);
  xdrmem_create (&bf->x_dat, bf->mem, bf->memsize, XDR_DECODE);

  /* empty current labels set */
  MAP_Free (&bf->mappool, &bf->labels);
 
  /* seek to the frame in IDX file */ 
  ASSERT (xdr_setpos (&bf->x_idx, bf->mtab [frm].ipos), ERR_PBF_INDEX_FILE_CORRUPTED);

  /* read labels */
  ASSERT (xdr_int (&bf->x_idx, &index), ERR_PBF_INDEX_FILE_CORRUPTED);
  while (index >= 0)
  {
    ASSERT (index < bf->lsize, ERR_PBF_INDEX_FILE_CORRUPTED);
    l = &bf->ltab [index];

    /* map label as current one */
    MAP_Insert (&bf->mappool, &bf->labels, l->name, l, (MAP_Compare) strcmp);

    /* read position */
    ASSERT (xdr_u_int (&bf->x_idx, &l->dpos), ERR_PBF_INDEX_FILE_CORRUPTED);

    /* get next label */
    ASSERT (xdr_int (&bf->x_idx, &index), ERR_PBF_INDEX_FILE_CORRUPTED);
  }
}

/* initialise labels and time index */
static void initialise_reading (PBF *bf)
{
  int index, num, siz;
  u_int dpos;
  
  /* create labels table */
  num = 0; siz = CHUNK;
  ERRMEM (bf->ltab = malloc (sizeof (PBF_LABEL) * siz));
  while (! feof (bf->lab))
  {
    bf->ltab [num].name = NULL;
    if (xdr_string (&bf->x_lab, &bf->ltab [num].name, PBF_MAXSTRING))
    {
      bf->ltab [num].index = num;
      if (++ num >= siz) { siz += CHUNK; ERRMEM (bf->ltab = realloc (bf->ltab, sizeof (PBF_LABEL) * siz)); }
    }
    else break;
  }
  bf->ltab = realloc (bf->ltab, sizeof (PBF_LABEL) * num); /* shrink */
  bf->lsize = num;
  
  /* create markers table */
  num = 0; siz = CHUNK;
  ERRMEM (bf->mtab = malloc (sizeof (PBF_MARKER) * siz));
  while (! feof (bf->idx))
  {
    /* time and unlabeled data position */
    ASSERT (xdr_double (&bf->x_idx, &bf->mtab [num].time), ERR_PBF_INDEX_FILE_CORRUPTED);
    ASSERT (xdr_uint64_t (&bf->x_idx, &bf->mtab [num].doff), ERR_PBF_INDEX_FILE_CORRUPTED);
    bf->mtab [num].ipos = xdr_getpos (&bf->x_idx);

    /* skip labels */
    ASSERT (xdr_int (&bf->x_idx, &index), ERR_PBF_INDEX_FILE_CORRUPTED);
    while (index >= 0)
    {
      ASSERT (index < bf->lsize, ERR_PBF_INDEX_FILE_CORRUPTED);
      ASSERT (xdr_u_int (&bf->x_idx, &dpos), ERR_PBF_INDEX_FILE_CORRUPTED);
      ASSERT (xdr_int (&bf->x_idx, &index), ERR_PBF_INDEX_FILE_CORRUPTED); 
    }

    if (index == -2) break; /* infinite frame */
    
    if (++ num >= siz) { siz += CHUNK; ERRMEM (bf->mtab = realloc (bf->mtab, sizeof (PBF_MARKER) * siz)); }
  }
  ASSERT (index == -2, ERR_PBF_INDEX_FILE_CORRUPTED);
  bf->mtab = realloc (bf->mtab, sizeof (PBF_MARKER) * (num + 1)); /* shrink (add INF frame) */
  bf->msize = num;

  /* read first frame */
  initialise_frame (bf, 0);
}

/* write last frame data */
static void write_frame (PBF *bf)
{
  if (xdr_getpos (&bf->x_idx) > 0)
  {
    /* mark end of frame labels */
    int index = -1;
    ASSERT (xdr_int (&bf->x_idx, &index), ERR_PBF_WRITE);

    /* write from XDR stream to file */
    filewrite (bf->mem, bf->membase + xdr_getpos (&bf->x_dat),
	       bf->dat, bf->compression == PBF_ON);

    /* rewind XDR */
    bf->membase = 0;
    xdr_destroy (&bf->x_dat);
    xdrmem_create (&bf->x_dat, bf->mem, bf->memsize, XDR_ENCODE);
  }
}

/* finalize frames after last write */
static void finalize_frames (PBF *bf)
{
  if (bf->mode == PBF_WRITE)
  {
    uint64_t doff;
    double time = DBL_MAX;
    int index = -2;

    /* write last frame */
    write_frame (bf);

    /* write infinite frame marker */
    ASSERT (xdr_double (&bf->x_idx, &time), ERR_PBF_WRITE);
    doff = (uint64_t) FTELL (bf->dat);
    ASSERT (xdr_uint64_t (&bf->x_idx, &doff), ERR_PBF_WRITE);
    ASSERT (xdr_int (&bf->x_idx, &index), ERR_PBF_WRITE);
  }
}
 
/* test before closing */
static int is_empty (FILE *f)
{
  fseek (f, 0, SEEK_END);
  return ftell (f) == 0;
}

/* copy path */
static char* copypath (const char *path)
{
  char *out = NULL;
  int l;

  if ((l = path ? strlen (path) : 0))
  {
    ERRMEM (out = malloc (l + 1));
    strcpy (out, path);
  }

  return out;
}

PBF* PBF_Write (const char *path)
{
  char *txt;
  PBF *bf;

#if MPI
  int rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
#endif

  ERRMEM (bf = malloc (sizeof (PBF)));
  ERRMEM (txt = malloc (strlen (path) + 16));
  bf->compression = PBF_OFF;
  bf->membase = 0;
  bf->memsize = CHUNK;
  ERRMEM (bf->mem = malloc (bf->memsize));

  /* openin files */
#if MPI
  sprintf (txt, "%s.dat.%d", path, rank);
#else
  sprintf (txt, "%s.dat", path);
#endif
  if (! (bf->dat = fopen (txt, "w"))) goto failure;
  xdrmem_create (&bf->x_dat, bf->mem, bf->memsize, XDR_ENCODE);
  bf->dph = copypath (txt);

#if MPI
  sprintf (txt, "%s.idx.%d", path, rank);
#else
  sprintf (txt, "%s.idx", path);
#endif
  if (! (bf->idx = fopen (txt, "w"))) goto failure;
  xdrstdio_create (&bf->x_idx, bf->idx, XDR_ENCODE);
  bf->iph = copypath (txt);

#if MPI
  sprintf (txt, "%s.lab.%d", path, rank);
#else
  sprintf (txt, "%s.lab", path);
#endif
  if (! (bf->lab = fopen (txt, "w"))) goto failure;
  xdrstdio_create (&bf->x_lab, bf->lab, XDR_ENCODE);
  bf->lph = copypath (txt);

  /* initialise the rest */
  MEM_Init (&bf->mappool, sizeof (MAP), CHUNK);
  MEM_Init (&bf->labpool, sizeof (PBF_LABEL), CHUNK);
  bf->ltab = NULL;
  bf->labels = NULL;
  bf->mtab = NULL;
  bf->mode = PBF_WRITE;
  bf->time = 0.;
  bf->lsize = 0;
  bf->msize = 0;
  bf->cur = 0;
#if MPI
  bf->parallel = PBF_ON;
#else
  bf->parallel = PBF_OFF;
#endif
  bf->next = NULL;
#if HDF5
#if MPI
  sprintf (txt, "%s.h5.%d", path, rank);
#else
  sprintf (txt, "%s.h5", path);
#endif
  if ((bf->stack[0] = H5Fcreate(txt, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) goto failure;
  bf->name[0] = copypath (txt);
  bf->frame = 0;
  bf->top = 0;
#endif
  free (txt);
  return bf;
  
failure: 
  free (bf);
  free (txt);
  return NULL;
}

PBF* PBF_Read (const char *path)
{
  PBF *bf, *out;
  FILE *dat;
  char *txt;
  int n, m;

  ERRMEM (txt = malloc (strlen (path) + 16));

  /* count input files */
  m = 0;
  do
  {
    sprintf (txt, "%s.dat.%d", path, m);
    dat = fopen (txt, "r");
  } while (dat && fclose (dat) == 0 && ++ m); /* m incremented as last */

  /* open input files */
  out = NULL;
  n = m-1;
  do
  {
    ERRMEM (bf = malloc (sizeof (PBF)));
    ERRMEM (bf->mem = malloc (CHUNK));
    bf->compression = PBF_OFF;
    bf->memsize = CHUNK;
    bf->membase = 0;

    /* openin files */
    if (m) sprintf (txt, "%s.dat.%d", path, n);
    else sprintf (txt, "%s.dat", path);
    if (! (bf->dat = fopen (txt, "r"))) goto failure;
    xdrmem_create (&bf->x_dat, bf->mem, bf->memsize, XDR_DECODE);
    bf->dph = copypath (txt);
    if (m) sprintf (txt, "%s.idx.%d", path, n);
    else sprintf (txt, "%s.idx", path);
    if (! (bf->idx = fopen (txt, "r"))) goto failure;
    xdrstdio_create (&bf->x_idx, bf->idx, XDR_DECODE);
    bf->iph = copypath (txt);
    if (m) sprintf (txt, "%s.lab.%d", path, n);
    else sprintf (txt, "%s.lab", path);
    if (! (bf->lab = fopen (txt, "r"))) goto failure;
    xdrstdio_create (&bf->x_lab, bf->lab, XDR_DECODE);
    bf->lph = copypath (txt);

#if HDF5
    if (m) sprintf (txt, "%s.h5.%d", path, n);
    else sprintf (txt, "%s.h5", path);
    if ((bf->stack[0] = H5Fopen(txt, H5F_ACC_RDONLY, H5P_DEFAULT)) < 0) goto failure;
    bf->name[0] = copypath (txt);
    bf->frame = 0;
    bf->top = 0;
#endif

    /* initialise the rest */
    MEM_Init (&bf->mappool, sizeof (MAP), CHUNK);
    MEM_Init (&bf->labpool, sizeof (PBF_LABEL), CHUNK);
    bf->mem = NULL;
    bf->ltab = NULL;
    bf->labels = NULL;
    bf->mtab = NULL;
    bf->mode = PBF_READ;
    bf->time = 0.;
    bf->lsize = 0;
    bf->msize = 0;
    bf->cur = 0;
    initialise_reading (bf);
    if (m) bf->parallel = PBF_ON;
    else bf->parallel = PBF_OFF;
    bf->next = out;
    out = bf;

  } while (-- n >= 0); /* the first item in the returned list corresponds to rank 0 */

  free (txt);
  return out;
  
failure: 
  free (bf->mem);
  free (bf);
  free (txt);
  return NULL;
}

void PBF_Close (PBF *bf)
{
  PBF *next;

  for (; bf; bf = next)
  {
    int empty;

    /* finalize writing */
    finalize_frames (bf);

    /* close streams */
    xdr_destroy (&bf->x_dat);
    xdr_destroy (&bf->x_idx);
    xdr_destroy (&bf->x_lab);

    empty = is_empty (bf->dat);

    fclose (bf->dat);
    fclose (bf->idx);
    fclose (bf->lab);

#if HDF5
    while (bf->top > 0) PBF_Pop_h5 (bf);
    H5Fclose (bf->stack[0]);
    free (bf->name[0]);
#endif

    if (empty) /* remove empty files */
    {
      remove (bf->dph);
      remove (bf->iph);
      remove (bf->lph);
    }
  
    free (bf->dph);
    free (bf->iph);
    free (bf->lph);
    free (bf->mem);

    /* free labels & markers */ 
    if (bf->mode == PBF_READ)
    {
      int k;

      for (k = 0; k < bf->lsize; k ++)
      {
	free (bf->ltab [k].name);
      }

      free (bf->ltab);
      free (bf->mtab);
    }
    else
    {
      MAP *k;

      for (k = MAP_First (bf->labels); k; k = MAP_Next (k))
      {
	PBF_LABEL *l = (PBF_LABEL*)k->data;
	free (l->name);
      }
      MEM_Release (&bf->labpool);
    }

    MEM_Release (&bf->mappool);
    next = bf->next;
    free (bf);
  }
}

/* read/write current time */

void PBF_Time (PBF *bf, double *time)
{
  uint64_t doff;

  if (bf->mode == PBF_WRITE)
  {
    ASSERT ((*time) >= bf->time, ERR_PBF_OUTPUT_TIME_DECREASED);

    /* write last frame */
    write_frame (bf);

    /* get current data offset */
    doff = (uint64_t) FTELL (bf->dat);

    /* output current data position and time */
    ASSERT (xdr_double (&bf->x_idx, time), ERR_PBF_WRITE);
    ASSERT (xdr_uint64_t (&bf->x_idx, &doff), ERR_PBF_WRITE);

    /* set time */
    bf->time = *time;

#if HDF5
    char name [128];
    snprintf (name, 128, "/%llu", bf->frame);
    while (bf->top > 0) PBF_Pop_h5 (bf);
    PBF_Push_h5 (bf, name);
    PBF_Double_h5 (bf, "time", time, 1);
    bf->frame ++;
#endif
  }
  else
  { 
    *time = bf->time;
  }
}

int PBF_Label (PBF *bf, const char *label)
{
  PBF_LABEL *l;
  u_int dpos;

  if (bf->mode == PBF_WRITE)
  {
    if (!(l = MAP_Find (bf->labels, (void*)label, (MAP_Compare) strcmp)))
    {
      /* create new label */
      l = MEM_Alloc (&bf->labpool);
      ERRMEM (l->name = malloc (strlen (label) + 1));
      strcpy (l->name, label);
      l->index = bf->lsize ++;
      MAP_Insert (&bf->mappool, &bf->labels, l->name, l, (MAP_Compare) strcmp);

      /* output definition */
      ASSERT (xdr_string (&bf->x_lab, (char**)&label, PBF_MAXSTRING), ERR_PBF_WRITE);
    }

    /* record label and position in the index file */
    ASSERT (xdr_int (&bf->x_idx, &l->index), ERR_PBF_WRITE);
    dpos = bf->membase + xdr_getpos (&bf->x_dat);
    ASSERT (xdr_u_int (&bf->x_idx, &dpos), ERR_PBF_WRITE);
  }
  else
  {
    if ((l = MAP_Find (bf->labels, (void*)label, (MAP_Compare) strcmp)))
    {
      /* seek to labeled data begining (relative displacement) */
      ASSERT (xdr_setpos (&bf->x_dat, l->dpos), ERR_PBF_READ);
    }
    else return 0;
  }

  return 1;
}

#define IO(type, call)\
  if (bf->mode == PBF_WRITE)\
  {\
    u_int nb = sizeof (type) * length + MARGIN;\
    if (bf->membase + xdr_getpos (&bf->x_dat) + nb >= bf->memsize) growmem (bf, nb);\
    ASSERT (xdr_vector (&bf->x_dat, (char*)value, length, sizeof (type), (xdrproc_t)call), ERR_PBF_WRITE);\
  }\
  else ASSERT (xdr_vector (&bf->x_dat, (char*)value, length, sizeof (type), (xdrproc_t)call), ERR_PBF_READ)

void PBF_Char (PBF *bf, char *value, unsigned int length)
{ 
  IO (char, xdr_char);
}

void PBF_Uchar (PBF *bf, unsigned char *value, unsigned int length)
{ 
  IO (unsigned char, xdr_u_char);
}

void PBF_Short (PBF *bf, short *value, unsigned int length)
{ 
  IO (short, xdr_short);
}

void PBF_Ushort (PBF *bf, unsigned short *value, unsigned int length)
{ 
  IO (unsigned short, xdr_u_short);
}

void PBF_Int (PBF *bf, int *value, unsigned int length)
{
  IO (int, xdr_int);
}

void PBF_Uint (PBF *bf, unsigned int *value, unsigned int length)
{ 
  IO (unsigned int, xdr_u_int);
}

void PBF_Long (PBF *bf, long *value, unsigned int length)
{ 
  IO (long, xdr_long);
}

void PBF_Ulong (PBF *bf, unsigned long *value, unsigned int length)
{ 
  IO (unsigned long, xdr_u_long);
}

void PBF_Float (PBF *bf, float *value, unsigned int length)
{ 
  IO (float, xdr_float);
}

void PBF_Double (PBF *bf, double *value, unsigned int length)
{ 
  IO (double, xdr_double);
}

void PBF_String (PBF *bf, char **value)
{ 
  if (bf->mode == PBF_WRITE)
  {
    u_int nb = strlen (*value) + MARGIN;
    if (bf->membase + xdr_getpos (&bf->x_dat) + nb >= bf->memsize) growmem (bf, nb);
    ASSERT (xdr_string (&bf->x_dat, value, PBF_MAXSTRING), ERR_PBF_WRITE);
  }
  else ASSERT (xdr_string (&bf->x_dat, value, PBF_MAXSTRING), ERR_PBF_READ);
}

#if HDF5
/* push group on stak */
void PBF_Push_h5 (PBF *bf, const char *name)
{
  bf->top ++;

  if (bf->mode == PBF_WRITE)
  {
    ASSERT (bf->top < PBF_MAXSTACK, ERR_PBF_WRITE);
    ASSERT ((bf->stack[bf->top] = H5Gcreate (bf->stack[bf->top-1], name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0, ERR_PBF_WRITE);
  }
  else
  {
    ASSERT (bf->top < PBF_MAXSTACK, ERR_PBF_READ);
    ASSERT ((bf->stack [bf->top] = H5Gopen (bf->stack[bf->top-1], name, H5P_DEFAULT)) >= 0,  ERR_PBF_WRITE);
  }

  bf->name[bf->top] = copypath (name);
}

/* then write datasets or attributes */
void PBF_Char_h5 (PBF *bf, const char *name, char *value, hsize_t length)
{
}

void PBF_Short_h5 (PBF *bf, const char *name, short *value, hsize_t length)
{
}

void PBF_Int_h5 (PBF *bf, const char *name, int *value, hsize_t length)
{
}

void PBF_Long_h5 (PBF *bf, const char *name, long *value, hsize_t length)
{
}

void PBF_Float_h5 (PBF *bf, const char *name, float *value, hsize_t length)
{
}

void PBF_Double_h5 (PBF *bf, const char *name, double *value, hsize_t length)
{
  if (bf->mode == PBF_WRITE)
  {
    if (length == 1)
    {
      ASSERT (H5LTset_attribute_double (bf->stack [bf->top], bf->name [bf->top], name, value, length) >= 0, ERR_PBF_WRITE);
    }
    else
    {
      ASSERT (H5LTmake_dataset_double (bf->stack [bf->top], name, 1, &length, value) >= 0, ERR_PBF_WRITE);
    }
  }
  else
  {
#if 0
    if (length == 1)
    {
      ASSERT (H5LTget_attribute_double (bf->stack [bf->top], bf->name [bf->top], name, value) >= 0, ERR_PBF_WRITE);
    }
    else
    {
      ASSERT (H5LTread_dataset_double (bf->stack [bf->top], name, value) >= 0, ERR_PBF_WRITE);
    }
#endif
  }
}

void PBF_String_h5 (PBF *bf, const char *name, char **value)
{
}

/* pop group from stack */
void PBF_Pop_h5 (PBF *bf)
{
  ASSERT_DEBUG (bf->top > 0, "PBF ERROR: too many pops!\n");

  H5Gclose (bf->stack [bf->top]);
  free (bf->name[bf->top]);

  bf->top --;
}
#endif

void PBF_Limits (PBF *bf, double *start, double *end)
{
  if (bf->mode == PBF_READ)
  {
    *start = bf->mtab [0].time;
    *end = bf->mtab [bf->msize - 1].time;
  }
}

void PBF_Seek (PBF *bf, double time)
{
  if (bf->mode == PBF_READ)
  {
    for (; bf; bf = bf->next)
    {
      PBF_MARKER *l, *h, *m;

      /* binary search
       * of marker */
      l = bf->mtab;
      h = l + bf->msize - 1;
      while (l <= h)
      {
	m = l + (h - l) / 2;
	if (time == m->time) break;
	else if (time < m->time) h = m - 1;
	else l = m + 1;
      }

      /* handle limit cases */
      if (h < bf->mtab) m = l;
      else if (l > (bf->mtab + bf->msize - 1)) m = h;
      initialise_frame (bf, m - bf->mtab);
    }
  }
}

int PBF_Backward (PBF *bf, unsigned int steps)
{
  if (bf->mode == PBF_READ)
  {
    int pos, ret;

    for (; bf; bf = bf->next)
    {
      if (bf->cur < steps) pos = 0, ret = 0;
      else pos = bf->cur - steps, ret = 1;
      initialise_frame (bf, pos);
    }

    return ret;
  }

  return 0;
}

int PBF_Forward (PBF *bf, unsigned int steps)
{
  if (bf->mode == PBF_READ)
  {
    int pos, ret;

    for (; bf; bf = bf->next)
    {
      if (bf->cur + steps >=  bf->msize) pos = bf->msize - 1, ret = 0;
      else pos = bf->cur + steps, ret = 1;
      initialise_frame (bf, pos);
    }

    return ret;
  }

  return 0;
}

unsigned int PBF_Span (PBF *bf, double t0, double t1)
{
  if (bf->mode == PBF_READ)
  {
    PBF_MARKER *l, *h, *m0, *m1;

    ASSERT_DEBUG (t0 <= t1, "t0 > t1");

    /* binary search
     * of marker */
    l = bf->mtab;
    h = l + bf->msize - 1;
    while (l <= h)
    {
      m0 = l + (h - l) / 2;
      if (t0 == m0->time) break;
      else if (t0 < m0->time) h = m0 - 1;
      else l = m0 + 1;
    }

    /* handle limit cases */
    if (h < bf->mtab) m0 = l;
    else if (l > (bf->mtab + bf->msize - 1)) m0 = h;

    /* binary search
     * of marker */
    l = bf->mtab;
    h = l + bf->msize - 1;
    while (l <= h)
    {
      m1 = l + (h - l) / 2;
      if (t1 == m1->time) break;
      else if (t1 < m1->time) h = m1 - 1;
      else l = m1 + 1;
    }

    /* handle limit cases */
    if (h < bf->mtab) m1 = l;
    else if (l > (bf->mtab + bf->msize - 1)) m1 = h;

    return m1 - m0;
  }

  return 0;
}
