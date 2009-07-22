/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <math.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/PrintFTSeries.h>

#include <lal/LALCalibration.h>
#include <lal/FrameCache.h>
#include <lal/LALFrameIO.h>
#include <lal/LALString.h>

static int parseargs( int argc, char **argv );
int epoch;
const char *cacheFile;
const char *calibrationFile;
const char *readoutChannel;
char defaultReadoutChannel[] = "Xn:LSC-DARM_ERR";

static int generate_file_name( char *fname, size_t size, const char *sname, int t, int dt );
static int write_REAL4TimeSeries( REAL4TimeSeries *series );
static int write_COMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series );

const char *program;
void err_exit( const char *fmt, ... )
{
  va_list ap;
  fprintf( stderr, "%s error: ", program );
  va_start( ap, fmt );
  vfprintf( stderr, fmt, ap );
  va_end( ap );
  fprintf( stderr, "\n" );
  exit( 1 );
}


int main( int argc, char *argv[] )
{
  LIGOTimeGPS tstart;
  /* COMPLEX8FrequencySeries *response  = NULL; */
  LALCalData *caldata;
  int t0, dt;
  char ifo[3];
  char site;
  const char *basename;

  /*
  if ( argc != 2 )
  {
    fprintf( stderr, "usage: %s calibration_frame_file\n", argv[0] );
    return 1;
  }
  */

  program = argv[0];
  lalDebugLevel = 7;
  XLALSetErrorHandler( XLALAbortErrorHandler );

  parseargs( argc, argv );
  if ( ! calibrationFile && ! cacheFile )
    err_exit( "no calibration frame file specified: try %s --help\n", program );

  if ( calibrationFile )
  {
    basename = strrchr( calibrationFile, '/' );
    basename = basename ? basename : calibrationFile;

    ifo[2] = 0;
    sscanf( basename, "%c-%c%c_%*s-%d-%d.gwf", &site, &ifo[0], &ifo[1], &t0, &dt );

    tstart.gpsSeconds     = epoch ? epoch : t0;
    tstart.gpsNanoSeconds = 0;
    memcpy( defaultReadoutChannel, ifo, 2 );
    readoutChannel = readoutChannel ? readoutChannel : defaultReadoutChannel;

    caldata = XLALFrGetCalData( &tstart, readoutChannel, calibrationFile );
    write_REAL4TimeSeries( caldata->cavityFactors );
    write_REAL4TimeSeries( caldata->openLoopFactors );
    write_COMPLEX8FrequencySeries( caldata->responseReference );
    write_COMPLEX8FrequencySeries( caldata->openLoopGainReference );
    write_COMPLEX8FrequencySeries( caldata->cavityGainReference );
    if ( caldata->actuationReference )
      write_COMPLEX8FrequencySeries( caldata->actuationReference );
    if ( caldata->digitalFilterReference )
      write_COMPLEX8FrequencySeries( caldata->digitalFilterReference );
  }
  else /* cache file */
  {
    FrCache *cache;
    basename = strrchr( cacheFile, '/' );
    basename = basename ? basename : cacheFile;

    ifo[2] = 0;
    sscanf( basename, "%c-%c%c_%*s-%d-%d.cache", &site, &ifo[0], &ifo[1], &t0, &dt );

    tstart.gpsSeconds     = epoch ? epoch : t0;
    tstart.gpsNanoSeconds = 0;
    memcpy( defaultReadoutChannel, ifo, 2 );
    readoutChannel = readoutChannel ? readoutChannel : defaultReadoutChannel;

    cache = XLALFrImportCache( cacheFile );
    caldata = XLALFrCacheGetCalData( &tstart, readoutChannel, cache );
    XLALFrDestroyCache( cache );

    write_REAL4TimeSeries( caldata->cavityFactors );
    write_REAL4TimeSeries( caldata->openLoopFactors );
    write_COMPLEX8FrequencySeries( caldata->responseReference );
    write_COMPLEX8FrequencySeries( caldata->openLoopGainReference );
    write_COMPLEX8FrequencySeries( caldata->cavityGainReference );
    if ( caldata->actuationReference )
      write_COMPLEX8FrequencySeries( caldata->actuationReference );
    if ( caldata->digitalFilterReference )
      write_COMPLEX8FrequencySeries( caldata->digitalFilterReference );
  }

  /*
  tstart.gpsSeconds += 1001;
  response = XLALCreateCOMPLEX8FrequencySeries( "response", &tstart, 0, 0.5 * 16384.0 / (1024*1024), &caldata->responseReference->sampleUnits, 1024*1024+1 );
  XLALUpdateResponse( response, 0.0, caldata );
  write_COMPLEX8FrequencySeries( response );
  */


  XLALDestroyCalData( caldata );

  return 0;
}


static int parseargs( int argc, char **argv )
{
  struct option long_options[] = 
  {
    { "help", no_argument, 0, 'h' },
    { "cache", required_argument, 0, 'c' },
    { "file", required_argument, 0, 'f' },
    { "epoch", required_argument, 0, 'e' },
    { "readout-channel", required_argument, 0, 'r' },
    { 0, 0, 0, 0 }
  };
  char args[] = "hc:f:e:r:";
  while ( 1 )
  {
    int option_index = 0;
    int c;

    c = getopt_long_only( argc, argv, args, long_options, &option_index );
    if ( c == -1 ) /* end of options */
      break;

    switch ( c )
    {
      case 0: /* if option set a flag, nothing else to do */
        if ( long_options[option_index].flag )
          break;
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
      case 'h': /* help */
        fprintf( stderr, "usage: %s [options]\n", argv[0] );
        fprintf( stderr, "options:\n" );
        fprintf( stderr, "\t-h, --help        \tprint this message and exit\n" );
        fprintf( stderr, "\t--cache CACHEFILE \tread calibration frame file from cache file\n" );
        fprintf( stderr, "\t--file FILENAME   \tread specified calibration frame file\n" );
        fprintf( stderr, "\t--epoch EPOCH     \tget calibration data for specified epoch\n" );
        fprintf( stderr, "\t--readout CHANNEL \tget calibration data for specified readout channel\n" );
        exit( 0 );
      case 'c': /* cache */
        cacheFile = optarg;
        break;
      case 'f': /* file */
        calibrationFile = optarg;
        break;
      case 'e': /* epoch */
        epoch = atoi( optarg );
        break;
      case 'r': /* readout-channel */
        readoutChannel = optarg;
        break;
      case '?':
      default:
        fprintf( stderr, "unknown error while parsing options\n" );
        exit( 1 );
    }
  }
  return 0;
}


static int generate_file_name( char *fname, size_t size, const char *sname, int t, int dt )
{
  char *tmp_name;
  char *c;

  tmp_name = XLALStringDuplicate( sname );

  /* slashes are not allowed */
  if ( strchr( tmp_name, '/' ) )
    err_exit( "slashes are not allowed in output file name %s\n", tmp_name );

  /* convert hyphens to underscores */
  while ( ( c = strchr( tmp_name, '-' ) ) )
    *c = '_';

  /* convert colons to underscores */
  while ( ( c = strchr( tmp_name, ':' ) ) )
    *c = '_';

  /* convert spaces to underscores */
  while ( ( c = strchr( tmp_name, ' ' ) ) )
    *c = '_';

  snprintf( fname, size, "%c-%s-%d-%d.dat", *tmp_name, tmp_name, t, dt );

  LALFree( tmp_name );

  return 0;
}



/* routine to write a time series */
static int write_REAL4TimeSeries( REAL4TimeSeries *series )
{
  char fname[FILENAME_MAX];
  int t0, dt;
  REAL8 t;
  size_t i;
  FILE *fp;
  t0 = series->epoch.gpsSeconds;
  dt = ceil(1e-9*series->epoch.gpsNanoSeconds + series->data->length*series->deltaT);
  t = XLALGPSGetREAL8( &series->epoch );
  generate_file_name( fname, sizeof( fname ), series->name, t0, dt );
  XLALPrintInfo( "writing series %s to file %s\n", series->name, fname );
  fp = fopen( fname, "w" );
  for ( i = 0; i < series->data->length; ++i )
    fprintf( fp, "%9f\t%e\n", t + i * series->deltaT, series->data->data[i] );
  fclose( fp );
  return 0;
}


/* routine to write a complex frequency series */
static int write_COMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series )
{
  char fname[FILENAME_MAX];
  int t0, dt;
  size_t i;
  FILE *fp;
  t0 = series->epoch.gpsSeconds;
  dt = ceil(1e-9*series->epoch.gpsNanoSeconds + 1.0/series->deltaF);
  generate_file_name( fname, sizeof( fname ), series->name, t0, dt );
  XLALPrintInfo( "writing series %s to file %s\n", series->name, fname );
  fp = fopen( fname, "w" );
  /* NOTE: omit first point! */
  for ( i = 1; i < series->data->length; ++i )
    fprintf( fp, "%9f\t%e\t%e\n", series->f0 + i * series->deltaF,
        hypot( series->data->data[i].im, series->data->data[i].re ),
        atan2( series->data->data[i].im, series->data->data[i].re ) );
  fclose( fp );
  return 0;
}

