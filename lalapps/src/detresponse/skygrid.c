/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#include <stdarg.h>

#include <ctype.h>
#include "config.h"
#include "cmdline.h"
#include "skygrid.h"

#include <lal/LALConfig.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/Units.h>

#include <lal/AVFactories.h>
#include <lal/VectorOps.h>


static UINT4 num_ra;
static UINT4 num_dec;
static UINT4 grid_lim;

static const double rad2deg = 180./LAL_PI;
static UINT4Vector *skygrid_dims = NULL;

void init_skygrid(LALStatus *status)
{
  num_ra = args_info.n_ra_arg;
  num_dec = args_info.n_dec_arg;
  grid_lim = num_ra * num_dec;

  LALU4CreateVector(status, &skygrid_dims, 2);
  skygrid_dims->data[0] = num_ra;
  skygrid_dims->data[1] = num_dec;
}

skygrid_t * alloc_skygrid(LALStatus *status, skygrid_t *g)
{
  UINT4 i = 0;

  if (g != NULL)
    LALDDestroyArray(status, g);

  LALDCreateArray(status, g, skygrid_dims);
    
  return g;
}

void free_skygrid(LALStatus *status, skygrid_t *skygrid)
{
  LALDDestroyArray(status, skygrid);
}

void cleanup_skygrid(LALStatus *status)
{
  LALU4DestroyVector(status, &skygrid_dims);
}

REAL4 skygrid_avg(LALStatus *status, const skygrid_t response)
{
  UINT4 i;
  REAL4 retval = 0.;
  
  for (i = 0; i < grid_lim; ++i)
    retval += response->data[i];
  
  retval /= grid_lim;
  
  return retval;
}




void skygrid_square(LALStatus *status, skygrid_t square, const skygrid_t input)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    square->data[i] = (input->data[i]) * (input->data[i]);
  
}




void skygrid_sqrt(LALStatus *status, skygrid_t result, const skygrid_t input)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    result->data[i] = (REAL4)sqrt((double)(input->data[i]));
}



REAL4 skygrid_rms(LALStatus *status, const skygrid_t input)
{
  skygrid_t tmpgrid;
  REAL4 retval;
  
  if (alloc_skygrid(status, &tmpgrid) == NULL)
  {
    fprintf(stderr, "skygrid_rms: failed to allocate\n");
    exit(1);
  }
  
  skygrid_square(status, tmpgrid, input);
  retval = (REAL4)(sqrt(skygrid_avg(status, tmpgrid)));
  
  free_skygrid(status, &tmpgrid);
  
  return retval;
}



INT4 skygrid_copy(LALStatus *status, skygrid_t dest, const skygrid_t src)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    dest->data[i] = src->data[i];
  
  return i;
}


void skygrid_print(LALStatus *status, 
                   const LIGOTimeGPS * gps,
                   const skygrid_t input, 
                   const char * filename)
{
  FILE * outfile = NULL;
  char write_mode[2];
  UINT4 i, j;
  
  write_mode[0] = '\0';
  write_mode[1] = '\0';
  
  /* we only have to be careful for time series output;
   * snapshot and average are just written to a single file */
  if (args_info.timeseries_given && 
      (tolower(args_info.format_arg[0]) == 's'))
    write_mode[0] = 'a';
  else
    write_mode[0] = 'w';
    
  outfile = xfopen(filename, write_mode);
  
  if (args_info.format_arg[1] == 'a')
  {
    if (args_info.format_arg[2] == 'l')
    {
      if (gps != (LIGOTimeGPS *)NULL)
        fprintf(outfile, "%.9d", gps->gpsSeconds);
      
      for (i = 0; i < grid_lim; ++i)
        fprintf(outfile, " %20.14e", input->data[i]); 
      fprintf(outfile, "\n");
      
    }
    else if (args_info.format_arg[2] == 'f')
    {
      if (gps != (LIGOTimeGPS *)NULL)
        fprintf(outfile, "%c %.9d\n", '%', gps->gpsSeconds);
      
      for (i = 0; i < grid_lim; ++i)
        fprintf(outfile, "% 20.14e\t", input->data[i]); 
      fprintf(outfile, "\n");      
    }
    else if (args_info.format_arg[2] == 'm')
    {
      if (gps != (LIGOTimeGPS *)NULL)
        fprintf(outfile, "%c %.9d\n", '%', gps->gpsSeconds);
      
      for (i = 0; i < num_ra; ++i)
      {
        for (j = 0; j < num_dec; ++j)
          fprintf(outfile, "% 20.14e\t", input->data[i*num_dec + j]);
        fprintf(outfile, "\n");
      }      
    }
    else
    {
      fprintf(stderr, "unknown format requested\n");
      exit(56);
    }
  }
  else if (args_info.format_arg[1] == 'b')
  {
    /* no comments for binary output */
    if (gps == (LIGOTimeGPS *)NULL)
    {
      fprintf(stderr, "must specify GPS timestamp for binary format\n");
      exit(22);
    }
    fwrite(&(gps->gpsSeconds), sizeof(INT4), 1, outfile);
    fwrite(&grid_lim, sizeof(UINT4), 1, outfile);
    fwrite(input->data, sizeof(REAL4), grid_lim, outfile);    
  }
  else
  {
    fprintf(stderr, "unknown format requested\n");
    exit(57);
  }

  xfclose(outfile);
} /* skygrid_print() */


void skygrid_fabs(LALStatus *status, skygrid_t absgrid, const skygrid_t input)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    absgrid->data[i] = fabs(input->data[i]);
}


void skygrid_add(LALStatus *status, 
                 skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    sum->data[i] = a->data[i] + b->data[i];
}

void skygrid_subtract(LALStatus *status, 
                      skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    sum->data[i] = a->data[i] - b->data[i];
}

void skygrid_scalar_mult(LALStatus *status, 
                         skygrid_t result, const skygrid_t a, REAL8 b)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    result->data[i] = b * a->data[i];
}

void skygrid_zero(LALStatus *status, skygrid_t a)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    a->data[i] = 0.;
}
