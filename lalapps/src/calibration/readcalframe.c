#include <math.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/PrintFTSeries.h>

#include <lal/LALCalibration.h>
#include <lal/LALFrameIO.h>
#include <lal/LALString.h>

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
  char fname[FILENAME_MAX] = "H-H1_CAL_S4_V4-793152133-2834520.gwf";
  char readoutChannel[] = "Xn:LSC-DARM_ERR";
  LIGOTimeGPS tstart;
  /* COMPLEX8FrequencySeries *response  = NULL; */
  LALCalData *caldata;
  int t0, dt;
  char ifo[3];
  char site;
  char *basename;

  if ( argc != 2 )
  {
    fprintf( stderr, "usage: %s calibration_frame_file\n", argv[0] );
    return 1;
  }

  program = argv[0];
  lalDebugLevel = 7;
  XLALSetErrorHandler( XLALAbortErrorHandler );

  basename = strrchr( argv[1], '/' );
  basename = basename ? basename : argv[1];

  ifo[2] = 0;
  sscanf( basename, "%c-%c%c_%*s-%d-%d.gwf", &site, &ifo[0], &ifo[1], &t0, &dt );

  tstart.gpsSeconds     = t0;
  tstart.gpsNanoSeconds = 0;
  memcpy( readoutChannel, ifo, 2 );

  caldata = XLALFrGetCalData( &tstart, readoutChannel, fname );
  write_REAL4TimeSeries( caldata->cavityFactors );
  write_REAL4TimeSeries( caldata->openLoopFactors );
  write_COMPLEX8FrequencySeries( caldata->responseReference );
  write_COMPLEX8FrequencySeries( caldata->openLoopGainReference );
  write_COMPLEX8FrequencySeries( caldata->cavityGainReference );

  /*
  tstart.gpsSeconds += 1001;
  response = XLALCreateCOMPLEX8FrequencySeries( "response", &tstart, 0, 0.5 * 16384.0 / (1024*1024), &caldata->responseReference->sampleUnits, 1024*1024+1 );
  XLALUpdateResponse( response, 0.0, caldata );
  write_COMPLEX8FrequencySeries( response );
  */


  XLALDestroyCalData( caldata );

  return 0;
}


static int generate_file_name( char *fname, size_t size, const char *sname, int t, int dt )
{
  char *c;
  char *tmp_name;

  tmp_name = XLALStringDuplicate( sname );

  /* slashes are not allowed */
  if ( strchr( tmp_name, '/' ) )
    err_exit( "slashes are not allowed in output file name %s\n", tmp_name );

  /* convert hyphens to underscores */
  while ( ( c = strchr( tmp_name, '-' ) ) )
    *c = '_';

  /* convert colons to hypens */
  while ( ( c = strchr( tmp_name, ':' ) ) )
    *c = '-';

  /* convert spaces to underscores */
  while ( ( c = strchr( tmp_name, ' ' ) ) )
    *c = '_';

  LALSnprintf( fname, size, "%s-%d-%d.dat", tmp_name, t, dt );

  LALFree( tmp_name );

  return 0;
}



/* routine to write a time series */
static int write_REAL4TimeSeries( REAL4TimeSeries *series )
{
  char fname[FILENAME_MAX];
  int t, dt;
  t  = series->epoch.gpsSeconds;
  dt = ceil(1e-9*series->epoch.gpsNanoSeconds + series->data->length*series->deltaT);
  generate_file_name( fname, sizeof( fname ), series->name, t, dt );
  XLALPrintInfo( "writing series %s to file %s\n", series->name, fname );
  LALSPrintTimeSeries( series, fname );
  return 0;
}


/* routine to write a complex frequency series */
static int write_COMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series )
{
  char fname[FILENAME_MAX];
  int t, dt;
  t  = series->epoch.gpsSeconds;
  dt = ceil(1e-9*series->epoch.gpsNanoSeconds + 1.0/series->deltaF);
  generate_file_name( fname, sizeof( fname ), series->name, t, dt );
  XLALPrintInfo( "writing series %s to file %s\n", series->name, fname );
  LALCPrintFrequencySeries( series, fname );
  return 0;
}

