/*
 * Author: David Chin <dwchin@umich.edu> +1-734-709-9119
 * $Id$
 */
#include <stdarg.h>

#include <ctype.h>
#include <assert.h>
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
static UINT4 start_ra;
static UINT4 count_ra;
static UINT4 num_dec;
static UINT4 start_dec;
static UINT4 count_dec;
static UINT4 grid_lim;

static const double rad2deg = 180./LAL_PI;
static UINT4Vector *skygrid_dims = NULL;

static const size_t PathNameLength = 256;  /* this should be long enough for path names */


void init_ephemeris(LALStatus *status, EphemerisData *p_ephemeris_data)
{
  CHAR *lal_prefix = NULL;
  CHAR  earthdat_path[PathNameLength];
  CHAR  sundat_path[PathNameLength];
  CHAR *share_path = "/share/lal/";
  CHAR *fn_earthdat = "/earth03-06.dat";
  CHAR *fn_sundat = "/sun03-06.dat";
  int i;
  
  for (i = 0; i < LALNameLength; ++i)
  {
    earthdat_path[i] = '\0';
    sundat_path[i] = '\0';
  }
  
  lal_prefix = getenv("LAL");
  
  if (lal_prefix != NULL)
  {
    if (lalDebugLevel > 3)
      printf("lal_prefix = %s\n", lal_prefix);
    
    if (strlen(earthdat_path) >= mystrlcat(earthdat_path, lal_prefix, LALNameLength))
    {
      fprintf(stderr, "oops!\n");
      exit(7);
    }
  
    if (strlen(sundat_path) >= mystrlcat(sundat_path, lal_prefix, LALNameLength))
    {
      fprintf(stderr, "oops!\n");
      exit(7);
    }
  }
  else
  {
    if (strlen(earthdat_path) >= mystrlcat(earthdat_path, ".", LALNameLength))
    {
      fprintf(stderr, "oops!\n");
      exit(7);
    }
    
    if (strlen(sundat_path) >= mystrlcat(sundat_path, ".", LALNameLength))
    {
      fprintf(stderr, "oops!\n");
      exit(7);
    }
  }
  
  if (strlen(earthdat_path) >= mystrlcat(earthdat_path, share_path, LALNameLength))
  {
    fprintf(stderr, "bah!\n");
    exit(7);
  }
    
  if (strlen(earthdat_path) >= mystrlcat(earthdat_path, fn_earthdat, LALNameLength))
  {
    fprintf(stderr, "bah!\n");
    exit(7);   
  }
    
  if (strlen(sundat_path) >= mystrlcat(sundat_path, share_path, LALNameLength))
  {
    fprintf(stderr, "bah!\n");
    exit(7);
  }
    
  if (strlen(sundat_path) >= mystrlcat(sundat_path, fn_sundat, LALNameLength))
  {
    fprintf(stderr, "bah!\n");
    exit(7);
  }
  
  if (lalDebugLevel > 3)
  {
    printf("LALNameLength = %d\n", LALNameLength);
    printf("earthdat_path = %s\n", earthdat_path);
    printf("strlen(earthdat_path) = %zd\n", strlen(earthdat_path));
    printf("sundat_path = %s\n", sundat_path);
    printf("strlen(sundat_path) = %zd\n", strlen(sundat_path));
  }
  
  p_ephemeris_data->ephiles.earthEphemeris = (char *)xcalloc(PathNameLength, sizeof(char));
  p_ephemeris_data->ephiles.sunEphemeris = (char *)xcalloc(PathNameLength, sizeof(char));
  
  if (args_info.earth_ephemeris_given)
  {
    (void)mystrlcpy(p_ephemeris_data->ephiles.earthEphemeris, 
                    args_info.earth_ephemeris_arg, PathNameLength);
  }
  else 
  {
    (void)mystrlcpy(p_ephemeris_data->ephiles.earthEphemeris, earthdat_path, PathNameLength);
  }
  
  if(args_info.sun_ephemeris_given)
  {
    (void)mystrlcpy(p_ephemeris_data->ephiles.sunEphemeris, 
                    args_info.sun_ephemeris_arg, PathNameLength);
  } 
  else 
  {
    (void)mystrlcpy(p_ephemeris_data->ephiles.sunEphemeris, sundat_path, PathNameLength);
  }
  
  LALInitBarycenter(status, p_ephemeris_data);
} /* END: init_ephemeris() */

void cleanup_ephemeris(LALStatus *status, EphemerisData *p_ephemeris_data)
{
  free(p_ephemeris_data->ephiles.earthEphemeris);
  free(p_ephemeris_data->ephiles.sunEphemeris); 
  LALFree(p_ephemeris_data->ephemE);
  LALFree(p_ephemeris_data->ephemS);
} /* END: cleanup_ephemeris() */

void init_skygrid(LALStatus *status)
{
  num_ra    = args_info.n_ra_arg;
  start_ra  = args_info.start_ra_arg;
  count_ra  = args_info.count_ra_arg;
  
  num_dec   = args_info.n_dec_arg;
  start_dec = args_info.start_dec_arg;
  count_dec = args_info.count_dec_arg;
  
  grid_lim  = count_ra * count_dec;
  
  if (lalDebugLevel)
  {
    printf("start_ra  = %d\n", start_ra);
    printf("count_ra  = %d\n", count_ra);
    printf("start_dec = %d\n", start_dec);
    printf("count_dec = %d\n", count_dec);
    printf("grid_lim  = %d\n", grid_lim);
  }
  
  /*
   * error checking
   */
  assert(count_ra <= num_ra - start_ra);
  assert(count_dec <= num_dec - start_dec);

  LALU4CreateVector(status, &skygrid_dims, 2);
  skygrid_dims->data[0] = count_ra;
  skygrid_dims->data[1] = count_dec;
} /* END: init_skygrid() */


skygrid_t * alloc_skygrid(LALStatus *status, skygrid_t *g)
{
  if (*g != NULL)
    LALDDestroyArray(status, g);

  LALDCreateArray(status, g, skygrid_dims);
    
  return g;
} /* END: alloc_skygrid() */


void free_skygrid(LALStatus *status, skygrid_t *g)
{
  if (lalDebugLevel)
  {
    printf("(*g)->data = %x\n", (*g)->data);
    printf("(*g)->dimLength->length = %x\n", (*g)->dimLength->length);
    printf("(*g)->dimLength->data[] = %d x %d\n", 
           (*g)->dimLength->data[0], (*g)->dimLength->data[1]);
  }
  
  LALDDestroyArray(status, g);
} /* END: free_skygrid() */


void cleanup_skygrid(LALStatus *status)
{
  LALU4DestroyVector(status, &skygrid_dims);
} /* END: cleanup_skygrid() */


REAL8 skygrid_avg(LALStatus *status, const skygrid_t response)
{
  UINT4 i;
  REAL8 retval = 0.;
  
  for (i = 0; i < grid_lim; ++i)
    retval += response->data[i];
  
  retval /= grid_lim;
  
  return retval;
} /* END: skygrid_avg() */




void skygrid_square(LALStatus *status, skygrid_t square, const skygrid_t input)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    square->data[i] = (input->data[i]) * (input->data[i]);
  
} /* END: skygrid_square() */




void skygrid_sqrt(LALStatus *status, skygrid_t result, const skygrid_t input)
{
  INT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    result->data[i] = (REAL4)sqrt((double)(input->data[i]));
} /* END: skygrid_sqrt() */



REAL8 skygrid_rms(LALStatus *status, const skygrid_t input)
{
  skygrid_t tmpgrid;
  REAL8 retval;
  
  if (alloc_skygrid(status, &tmpgrid) == NULL)
  {
    fprintf(stderr, "skygrid_rms: failed to allocate\n");
    exit(1);
  }
  
  skygrid_square(status, tmpgrid, input);
  retval = sqrt(skygrid_avg(status, tmpgrid));
  
  free_skygrid(status, &tmpgrid);
  
  return retval;
} /* END: skygrid_rms() */



INT4 skygrid_copy(LALStatus *status, skygrid_t dest, const skygrid_t src)
{
  UINT4 i;
  
  /* FIXME: replace with memcpy(3) */
  for (i = 0; i < grid_lim; ++i)
    dest->data[i] = src->data[i];
  
  return i;
} /* END: skygrid_copy() */


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
      
      for (i = 0; i < count_ra; ++i)
      {
        for (j = 0; j < count_dec; ++j)
          fprintf(outfile, "% 20.14e\t", input->data[i*count_dec + j]);
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
    fwrite(input->data, sizeof(REAL8), grid_lim, outfile);    
  }
  else
  {
    fprintf(stderr, "unknown format requested\n");
    exit(57);
  }

  xfclose(outfile);
} /* END: skygrid_print() */


void skygrid_fabs(LALStatus *status, skygrid_t absgrid, const skygrid_t input)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    absgrid->data[i] = fabs(input->data[i]);
} /* END: skygrid_fabs() */


void skygrid_add(LALStatus *status, 
                 skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    sum->data[i] = a->data[i] + b->data[i];
} /* END: skygrid_add() */

void skygrid_subtract(LALStatus *status, 
                      skygrid_t sum, const skygrid_t a, const skygrid_t b)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    sum->data[i] = a->data[i] - b->data[i];
} /* END: skygrid_subtract() */

void skygrid_scalar_mult(LALStatus *status, 
                         skygrid_t result, const skygrid_t a, REAL8 b)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    result->data[i] = b * a->data[i];
} /* END: skygrid_scalar() */

void skygrid_zero(LALStatus *status, skygrid_t a)
{
  UINT4 i;
  
  for (i = 0; i < grid_lim; ++i)
    a->data[i] = 0.;
} /* END: skygrid_zero() */
