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
static char *path_and_dirs (char *path, char *ext)
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
  int n = l + strlen (lastname) + strlen(ext) + 8;
  ERRMEM (out = malloc (n));
  strcpy (out, path);
  out [l] = '/';
  strcpy (out+l+1, lastname);
  l = strlen (out);
  sprintf (&out[l], "%s", ext);

  return out;
}

/* read int attribute from HDF5 file path */
static int read_int (hid_t file, const char *path, const char *name)
{
  int value;

  ASSERT_TEXT (H5LTget_attribute_int (file, path, name, &value) >= 0, "HDF5 file write error");

  return value;
}

/* read double attribute from HDF5 file path */
static double read_double (hid_t file, const char *path, const char *name)
{
  double value;

  ASSERT_TEXT (H5LTget_attribute_double (file, path, name, &value) >= 0, "HDF5 file write error");

  return value;
}

/* read string attribute from HDF5 file path */
static char* read_string (hid_t file, const char *path, const char *name)
{
  char *value;
  H5T_class_t c;
  hsize_t d;
  size_t s;

  ASSERT_TEXT (H5LTget_attribute_info (file, path, name,  &d, &c, &s) >= 0, "HDF5 file write error"); 
  ERRMEM (value = malloc (s));
  ASSERT_TEXT (H5LTget_attribute_string (file, path, name, value) >= 0, "HDF5 file write error");

  return value;
}

/* write current mesh state */
static void write_mesh (MESH *msh, char *name, double time, int step, hid_t h5_body)
{
  if (step == 0)
  {
    int *topo = NULL, topo_size = 0, topo_count = 0;
    ELEMENT *ele;

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
    ASSERT_TEXT (H5LTmake_dataset_int (h5_body, "TOPO", 1, &length, topo) >= 0, "HDF5 write error");
    free (topo);
    ASSERT_TEXT (H5LTset_attribute_int (h5_body, ".", "TOPO_SIZE", &topo_count, 1) >= 0, "HDF5 file write error");
    int elements = msh->surfeles_count + msh->bulkeles_count;
    ASSERT_TEXT (H5LTset_attribute_int (h5_body, ".", "ELEMENTS", &elements, 1) >= 0, "HDF5 file write error");
    ASSERT_TEXT (H5LTset_attribute_int (h5_body, ".", "NODES", &msh->nodes_count, 1) >= 0, "HDF5 file write error");
  }

  char h5_text[1024];
  hid_t h5_step;
  
  snprintf (h5_text, 1024, "%d", step);
  ASSERT_TEXT ((h5_step = H5Gcreate (h5_body, h5_text, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file write error");
  ASSERT_TEXT (H5LTset_attribute_double (h5_step, ".", "TIME", &time, 1) >= 0, "HDF5 file write error");
  hsize_t dims[2] = {msh->nodes_count, 3};
  ASSERT_TEXT (H5LTmake_dataset_double (h5_step, "GEOM", 2, dims, (double*)msh->cur_nodes) >= 0, "HDF5 write error");
  H5Gclose (h5_step);
}

/* write bodies state at current time step */
static void write_bodies (DOM *dom, hid_t h5_file, int **bid, int *bid_count, int *bid_size)
{
  hid_t h5_bodies;
  BODY *bod;

  if (H5Lexists (h5_file, "BODIES", H5P_DEFAULT))
  {
    ASSERT_TEXT (h5_bodies = H5Gopen (h5_file, "BODIES", H5P_DEFAULT), "HDF5 file write error");
  }
  else
  {
    ASSERT_TEXT ((h5_bodies = H5Gcreate (h5_file, "BODIES", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file write error");
  }

  for (bod = dom->bod; bod; bod = bod->next)
  {
    SHAPE *shp = bod->shape;
    char h5_text[1024];
    hid_t h5_body;
    int steps;

    if (shp->kind != SHAPE_MESH) continue; /* FIXME */

    snprintf (h5_text, 1024, "%d", bod->id);

    if (H5Lexists (h5_bodies, h5_text, H5P_DEFAULT))
    {
      ASSERT_TEXT (h5_body = H5Gopen (h5_bodies, h5_text, H5P_DEFAULT), "HDF5 file write error");
      ASSERT_TEXT (H5LTget_attribute_int (h5_body, ".", "STEPS", &steps) >= 0, "HDF5 file write error");
      steps ++;
      ASSERT_TEXT (H5LTset_attribute_int (h5_body, ".", "STEPS", &steps, 1) >= 0, "HDF5 file write error");
    }
    else 
    {
      steps = 1;
      pack_int (bid_size, bid, bid_count, bod->id);
      ASSERT_TEXT ((h5_body = H5Gcreate (h5_bodies, h5_text, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file write error");
      ASSERT_TEXT (H5LTset_attribute_int (h5_body, ".", "STEPS", &steps, 1) >= 0, "HDF5 file write error");
      if (bod->label) ASSERT_TEXT (H5LTset_attribute_string (h5_body, ".", "LABEL", bod->label) >= 0, "HDF5 file write error");
      else ASSERT_TEXT (H5LTset_attribute_string (h5_body, ".", "LABEL", h5_text) >= 0, "HDF5 file write error");
    }

    switch (shp->kind)
    {
    case SHAPE_MESH:
      write_mesh (shp->data, bod->label, dom->time, steps-1, h5_body);
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

    H5Gclose (h5_body);
  }

  H5Gclose (h5_bodies);
}

/* Export results in XMDF format;
 * ntimes > 0 --> number of individual time instances;
 * ntimes < 0 --> a time interval from times[0] to times[1];
 */
void xdmf_export (SOLFEC *sol, double *times, int ntimes, char *path)
{
  /* First --> write heavu data into a HDF5 file */
  int *bid = NULL, bid_count = 0, bid_size = 0;
  char *h5_path = path_and_dirs (path, ".h5");
  hsize_t h5_size;
  hid_t h5_file;

  ASSERT_TEXT ((h5_file = H5Fcreate(h5_path, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file open error");

  if (ntimes < 0)
  {
    double start, end, t0, t1;

    SOLFEC_Time_Limits (sol, &start, &end);

    t0 = MAX (start, times[0]);

    t1 = MIN (end, times[1]);

    SOLFEC_Seek_To (sol, t0);

    do
    {
      write_bodies (sol->dom, h5_file, &bid, &bid_count, &bid_size);

      SOLFEC_Forward (sol, 1);
    }
    while (sol->dom->time < t1);
  }
  else
  {
    for (int i = 0; i < ntimes; i ++)
    {
      SOLFEC_Seek_To (sol, times[i]);

      write_bodies (sol->dom, h5_file, &bid, &bid_count, &bid_size);
    }
  }

  ASSERT_TEXT (H5LTset_attribute_int (h5_file, ".", "BODY_COUNT", &bid_count, 1) >= 0, "HDF5 file write error");
  h5_size = bid_count;
  ASSERT_TEXT (H5LTmake_dataset_int (h5_file, "BODY_IDS", 1, &h5_size, bid) >= 0, "HDF5 file write error");

  /* Second --> using the information from the HDF5 file write XDMF file */
  char *xmf_path = path_and_dirs (path, ".xmf");
  FILE *xmf_file = fopen (xmf_path, "w");
  ASSERT_TEXT (xmf_file, "Opening XDMF markup file %s has failed", xmf_path);

  fprintf (xmf_file, "<Xdmf>\n");
  fprintf (xmf_file, "<Domain>\n");
  fprintf (xmf_file, "<Grid GridType=\"Collection\" CollectionType=\"Spatial\">\n\n");

  int k = strlen (h5_path) - 1;
  while (k >= 0 && h5_path[k] != '/') k --;

  for (int i = 0; i < bid_count; i ++)
  {
    char h5_text [1024];

    snprintf (h5_text, 1024, "/BODIES/%d", bid[i]);

    fprintf (xmf_file, "<Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");

    int steps = read_int  (h5_file, h5_text, "STEPS");
    int elements = read_int (h5_file, h5_text, "ELEMENTS");
    int nodes = read_int (h5_file, h5_text, "NODES");
    int topo_size = read_int (h5_file, h5_text, "TOPO_SIZE");
    char *label = read_string (h5_file, h5_text, "LABEL");

    for (int j = 0; j < steps; j ++)
    {
      fprintf (xmf_file, "<Grid Name=\"%s\" Type=\"Uniform\">\n", label);
      snprintf (h5_text, 1024, "/BODIES/%d/%d", bid[i], j);
      double time = read_double (h5_file, h5_text, "TIME");
      fprintf (xmf_file, "<Time Type=\"Single\" Value=\"%f\" />\n", time);
      fprintf (xmf_file, "<Topology Type=\"Mixed\" NumberOfElements=\"%d\">\n", elements);
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\">\n", topo_size);
      fprintf (xmf_file, "%s:/BODIES/%d/TOPO\n", &h5_path[k+1], bid[i]);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Topology>\n");
      fprintf (xmf_file, "<Geometry GeometryType=\"XYZ\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/BODIES/%d/%d/GEOM\n", &h5_path[k+1], bid[i], j);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Geometry>\n");
      fprintf (xmf_file, "</Grid>\n");
    }

    free (label);

    fprintf (xmf_file, "</Grid>\n\n");
  }

  fprintf (xmf_file, "</Grid>\n");
  fprintf (xmf_file, "</Domain>\n");
  fprintf (xmf_file, "</Xdmf>\n");

  /* Clean up */
  free (bid);
  free (h5_path);
  free (xmf_path);
  fclose (xmf_file);
  H5Fclose (h5_file);
}
#else
void xdmf_export (SOLFEC *sol, double *times, int ntimes, char *path)
{
  fprintf (stderr, "Error: XDMF export is not supported without HDF5 --> Re-compile Soflec with HDF5 support.\n");
}
#endif
