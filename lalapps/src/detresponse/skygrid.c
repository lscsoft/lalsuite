/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#include <ctype.h>
#include "config.h"
#include "skygrid.h"

static const INT4 grid_lim = NUM_RA * NUM_DEC;
static const INT4 dec_lim  = (NUM_DEC - 1)/2;

static const double rad2deg = 180./LAL_PI;

REAL4 skygrid_avg(const skygrid_t response)
{
  INT4 i, j, cnt;
  REAL4 retval = 0.;
  
  for (j = 0; j < NUM_RA; ++j)
  {
    for (i = -dec_lim+1; i <= dec_lim-1; ++i)
    {
      cnt = j*NUM_DEC + i + dec_lim;
      retval += response[cnt];
    }
  }
  
  retval /= grid_lim;
  
  return retval;
}




void skygrid_square(skygrid_t square, const skygrid_t input)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    square[i] = (input[i]) * (input[i]);
  
}




void skygrid_sqrt(skygrid_t result, const skygrid_t input)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    result[i] = (REAL4)sqrt((double)(input[i]));
}



REAL4 skygrid_rms(const skygrid_t input)
{
  skygrid_t tmpgrid;
  
  skygrid_square(tmpgrid, input);
  return (REAL4)(sqrt(skygrid_avg(tmpgrid)));
}



INT4 skygrid_copy(skygrid_t dest, const skygrid_t src)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    dest[i] = src[i];
  
  return i;
}


void skygrid_print(const LIGOTimeGPS * gps,
                   const skygrid_t input, 
                   const char * filename)
{
  FILE * outfile = NULL;
  char write_mode[2];
  int i, j;
  
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
      const int tot = NUM_RA * NUM_DEC;
      if (gps != (LIGOTimeGPS *)NULL)
        fprintf(outfile, "%.9d", gps->gpsSeconds);
      
      for (i = 0; i < tot; ++i)
        fprintf(outfile, " %14.8e", input[i]); 
      fprintf(outfile, "\n");
      
    }
    else if (args_info.format_arg[2] == 'f')
    {
      const int tot = NUM_RA * NUM_DEC;
      if (gps != (LIGOTimeGPS *)NULL)
        fprintf(outfile, "%c %.9d\n", '%', gps->gpsSeconds);
      
      for (i = 0; i < tot; ++i)
        fprintf(outfile, "% 14.8e\t", input[i]); 
      fprintf(outfile, "\n");      
    }
    else if (args_info.format_arg[2] == 'm')
    {
      if (gps != (LIGOTimeGPS *)NULL)
        fprintf(outfile, "%c %.9d\n", '%', gps->gpsSeconds);
      
      for (i = 0; i < NUM_RA; ++i)
      {
        for (j = 0; j < NUM_DEC; ++j)
          fprintf(outfile, "% 14.8e\t", input[i*NUM_DEC + j]);
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
    const size_t num_el = NUM_DEC*NUM_RA;
    if (gps == (LIGOTimeGPS *)NULL)
    {
      fprintf(stderr, "must specify GPS timestamp for binary format\n");
      exit(22);
    }
    fwrite(&(gps->gpsSeconds), sizeof(INT4), 1, outfile);
    fwrite(&num_el, sizeof(INT4), 1, outfile);
    fwrite(input, sizeof(REAL4), num_el, outfile);    
  }
  else
  {
    fprintf(stderr, "unknown format requested\n");
    exit(57);
  }

  xfclose(outfile);
} /* skygrid_print() */


void skygrid_fabs(skygrid_t absgrid, const skygrid_t input)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    absgrid[i] = fabs(input[i]);
}


void skygrid_add(skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    sum[i] = a[i] + b[i];
}

void skygrid_subtract(skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    sum[i] = a[i] - b[i];
}

void skygrid_scalar_mult(skygrid_t result, const skygrid_t a, REAL4 b)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    result[i] = b * a[i];
}

void skygrid_zero(skygrid_t a)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    a[i] = 0.;
}
