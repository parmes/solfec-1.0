/*
The MIT License (MIT)

Copyright (c) 2016 EDF Energy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* Contributors: Tomasz Koziara */

#if HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif
#if POSIX
#include <sys/stat.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "sol.h"
#include "pck.h"
#include "err.h"

#if HDF5
/* returns xmf file path string; creatses intermediate directories */
static char *xmf_path_and_dirs (char *path)
{
  char *out;

  int l = strlen (path);

#if POSIX
  for (int i = 0; i < l; i ++) /* create all intermediate directories */
  {
    if (path [i] == '/')
    {
       path [i] = '\0';
       mkdir (path, 0777); /* POSIX */
       path [i] = '/';
    }
  }
  mkdir (path, 0777); /* POSIX */
#endif

  while (l > 0 && path [l-1] != '/') l --;
  char *lastname = &path [l];
  l = strlen (path);
  int n = l + strlen (lastname) + 8;
  ERRMEM (out = malloc (n));
  strcpy (out, path);
  out [l] = '/';
  strcpy (out+l+1, lastname);
  l = strlen (out);
  sprintf (&out[l], ".xmf");

  return out;
}

/* export mesh */
static void xdmf_mesh_export (MESH *msh, char *name, int bid, double time, FILE *xmf_file, hid_t h5_path, hid_t h5_file, short h5_first)
{
  int topo_count = 0;
  ELEMENT *ele;

  if (h5_first)
  {
    int *topo = NULL, topo_size = 0;
    for (ele = msh->surfeles; ele; ele = ele->next)
    {
      switch (ele->type)
      {
      case 4:
	pack_int (&topo_size, &topo, &topo_count, 6);
      break;
      case 5:
	pack_int (&topo_size, &topo, &topo_count, 7);
      break;
      case 6:
	pack_int (&topo_size, &topo, &topo_count, 8);
      break;
      case 8:
	pack_int (&topo_size, &topo, &topo_count, 9);
      break;
      }
      pack_ints (&topo_size, &topo, &topo_count, ele->nodes, ele->type);
    }
    for (ele = msh->bulkeles; ele; ele = ele->next)
    {
      switch (ele->type)
      {
      case 4:
	pack_int (&topo_size, &topo, &topo_count, 6);
      break;
      case 5:
	pack_int (&topo_size, &topo, &topo_count, 7);
      break;
      case 6:
	pack_int (&topo_size, &topo, &topo_count, 8);
      break;
      case 8:
	pack_int (&topo_size, &topo, &topo_count, 9);
      break;
      }
      pack_ints (&topo_size, &topo, &topo_count, ele->nodes, ele->type);
    }
    hsize_t length = topo_count;
    ASSERT_TEXT (H5LTmake_dataset_int (h5_file, "TOPO", 1, &length, topo) >= 0, "HDF5 write error");
    free (topo);
  }
  else
  {
    for (ele = msh->surfeles; ele; ele = ele->next)
    {
      topo_count += ele->type + 1;
    }
    for (ele = msh->bulkeles; ele; ele = ele->next)
    {
      topo_count += ele->type + 1;
    }
  }

  hsize_t dims[2] = {msh->nodes_count, 3};
  ASSERT_TEXT (H5LTmake_dataset_double (h5_path, "GEOM", 2, dims, (double*)msh->cur_nodes) >= 0, "HDF5 write error");

  if (name) fprintf (xmf_file, "<Grid Name=\"%s\" Type=\"Uniform\">\n", name);
  else fprintf (xmf_file, "<Grid Name=\"%d\" Type=\"Uniform\">\n", bid);
  fprintf (xmf_file, "<Time Type=\"Single\" Value=\"%f\" />\n", time);
  fprintf (xmf_file, "<Topology Type=\"Mixed\" NumberOfElements=\"%d\">\n", msh->surfeles_count + msh->bulkeles_count);
  fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\">\n", topo_count);
  fprintf (xmf_file, "%d.h5:/TOPO\n", bid);
  fprintf (xmf_file, "</DataStructure>\n");
  fprintf (xmf_file, "</Topology>\n");
  fprintf (xmf_file, "<Geometry GeometryType=\"XYZ\">\n");
  fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", msh->nodes_count);
  fprintf (xmf_file, "%d.h5:/%f/GEOM\n", bid, time);
  fprintf (xmf_file, "</DataStructure>\n");
  fprintf (xmf_file, "</Geometry>\n");
  fprintf (xmf_file, "</Grid>\n");
}

/* export domain state at current time */
static void xdmf_step (DOM *dom, FILE *xmf_file, char *path)
{
  BODY *bod;

  fprintf (xmf_file, "<Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n");

  for (bod = dom->bod; bod; bod = bod->next)
  {
    int h5_text_len = strlen(path) + 128;
    SHAPE *shp = bod->shape;
    char *h5_text;
    hid_t h5_file;
    hid_t h5_path;
    FILE *h5_test;
    short h5_first = 0;

    if (shp->kind != SHAPE_MESH) continue; /* FIXME */

    ERRMEM (h5_text = malloc(h5_text_len));

    snprintf (h5_text, h5_text_len, "%s/%d.h5", path, bod->id);
    if ((h5_test = fopen (h5_text, "r")) != NULL)
    {
      fclose (h5_test);
      ASSERT_TEXT ((h5_file = H5Fopen(h5_text, H5F_ACC_RDWR, H5P_DEFAULT)) >= 0, "Opening HDF5 file failed");
    }
    else 
    {
      h5_first = 1;
      ASSERT_TEXT ((h5_file = H5Fcreate(h5_text, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "Opening HDF5 file failed");
    }

    snprintf (h5_text, h5_text_len, "/%f", dom->time);
    ASSERT_TEXT ((h5_path = H5Gcreate (h5_file, h5_text, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 write error");

    switch (shp->kind)
    {
    case SHAPE_MESH:
      xdmf_mesh_export (shp->data, bod->label, bod->id, dom->time, xmf_file, h5_path, h5_file, h5_first);
    break;
    case SHAPE_CONVEX:
      /* TODO */
    break;
    case SHAPE_SPHERE:
      /* TODO */
    break;
    case SHAPE_ELLIP:
      /* TODO */
    break;
    }

    H5Gclose (h5_path);
    H5Fclose (h5_file);
    free (h5_text);
  }

  fprintf (xmf_file, "</Grid>\n");
}

/* Export results in XMDF format;
 * ntimes > 0 --> number of individual time instances;
 * ntimes < 0 --> a time interval from times[0] to times[1];
 */
void xdmf_export (SOLFEC *sol, double *times, int ntimes, char *path)
{
  /* TODO: overwrite behavior following Solfec input flags */
  char *xmf_path = xmf_path_and_dirs (path);
  FILE *xmf_file = fopen (xmf_path, "w");
  ASSERT_TEXT (xmf_file, "Opening XDMF markup file %s has failed", xmf_path);

  fprintf (xmf_file, "<Xdmf>\n");
  fprintf (xmf_file, "<Domain>\n");
  fprintf (xmf_file, "<Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");

#if 1
  if (ntimes < 0)
  {
    double start, end, t0, t1;

    SOLFEC_Time_Limits (sol, &start, &end);

    t0 = MAX (start, times[0]);

    t1 = MIN (end, times[1]);

    SOLFEC_Seek_To (sol, t0);

    do
    {
      xdmf_step (sol->dom, xmf_file, path);

      SOLFEC_Forward (sol, 1);
    }
    while (sol->dom->time < t1);
  }
  else
  {
    for (int i = 0; i < ntimes; i ++)
    {
      SOLFEC_Seek_To (sol, times[i]);

      xdmf_step (sol->dom, xmf_file, path);
    }
  }
#else
  SOLFEC_Seek_To (sol, 0.0);

  xdmf_step (sol->dom, xmf_file, path);
#endif

  fprintf (xmf_file, "</Grid>\n");
  fprintf (xmf_file, "</Domain>\n");
  fprintf (xmf_file, "</Xdmf>\n");

  fclose (xmf_file);
  free (xmf_path);
}
#else
void xdmf_export (SOLFEC *sol, double *times, int ntimes, char *path)
{
  fprintf (stderr, "Error: XDMF export is not supported without HDF5 --> Re-compile Soflec with HDF5 support.\n");
}
#endif
