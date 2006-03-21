int isinf( double );
int isnan( double );
#include <math.h>
#include <limits.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

#include "LALASCIIFileRead.h"
#include "LALFrameIO.h"
#include "LALString.h"

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
  char rcsid[] = "$Id$";
  char *history = NULL;
  LALDataFileNameFields *fields = NULL;
  char fname[FILENAME_MAX] = "";
  char run[FILENAME_MAX] = "";
  char ifo[3] = "";
  int version = 9999;
  int tstart = 0;
  int tend = 0;
  int have_fac = 0;
  int have_oloop_gain = 0;
  int have_cav_gain_as_q = 0;
  int have_response_as_q = 0;
  int have_cav_gain_darm_err = 0;
  int have_response_darm_err = 0;
  int detectorFlags;
  LIGOTimeGPS epoch;
  FrFile *frfile;
  FrameH *frame;
  int arg;

  program = argv[0];
  lalDebugLevel = 7;
  XLALSetErrorHandler( XLALAbortErrorHandler );

  history = XLALStringAppend( history, argv[0] );
  for ( arg = 1; arg < argc; ++arg )
  {
    history = XLALStringAppend( history, " " );
    history = XLALStringAppend( history, argv[arg] );
  }

  fprintf( stderr, "Input files:\n\n" );

  /* first pass: just parse filenames for sanity */
  fields = XLALCalloc( argc, sizeof( *fields ) );
  for ( arg = 1; arg < argc; ++arg )
  {
    const char *filename = argv[arg];

    fprintf( stderr, "\t%s\n", filename );

    XLALDataFileNameParse( fields + arg, filename );

    if ( strstr( fields[arg].description, "CAL_REF" ) )
    {
      LALCalRefFileNameDescriptionFields dscfields;
      XLALCalRefFileNameDescriptionParse( &dscfields, fields[arg].description );
      if ( arg == 1 )
      {
        version = dscfields.version;
        XLALStringCopy( ifo, dscfields.ifo, sizeof( ifo ) );
        XLALStringCopy( run, dscfields.run, sizeof( run ) );
      }
      else
      {
        if ( dscfields.version != version )
          err_exit( "incompatible version numbers" );
        if ( strcmp( ifo, dscfields.ifo ) )
          err_exit( "incompatible interferometers" );
        if ( strcmp( run, dscfields.run ) )
          err_exit( "incompatible runs" );
      }
      if ( ! strcmp( dscfields.channelPostfix, "OLOOP_GAIN" ) )
        ++have_oloop_gain;
      else if ( ! strcmp( dscfields.channelPostfix, "CAV_GAIN_AS_Q" ) )
        ++have_cav_gain_as_q;
      else if ( ! strcmp( dscfields.channelPostfix, "RESPONSE_AS_Q" ) )
        ++have_response_as_q;
      else if ( ! strcmp( dscfields.channelPostfix, "CAV_GAIN_DARM_ERR" ) )
        ++have_cav_gain_darm_err;
      else if ( ! strcmp( dscfields.channelPostfix, "RESPONSE_DARM_ERR" ) )
        ++have_response_darm_err;
      else
        fprintf( stderr, "%s warning: unrecognized reference function for file %s\n", program, filename );
    }
    else if ( strstr( fields[arg].description, "CAL_FAC" ) )
    {
      LALCalFacFileNameDescriptionFields dscfields;
      XLALCalFacFileNameDescriptionParse( &dscfields, fields[arg].description );
      if ( arg == 1 )
      {
        version = dscfields.version;
        XLALStringCopy( ifo, dscfields.ifo, sizeof( ifo ) );
        XLALStringCopy( run, dscfields.run, sizeof( run ) );
      }
      else
      {
        if ( dscfields.version != version )
          err_exit( "incompatible version numbers for file %s", filename );
        if ( strcmp( ifo, dscfields.ifo ) )
          err_exit( "incompatible interferometers for file %s", filename );
        if ( strcmp( run, dscfields.run ) )
          err_exit( "incompatible runs for file %s", filename );
      }
      ++have_fac;
    }
    else
      err_exit( "invalid file name %s", filename );

    /* other sanity checks */
    if ( fields[arg].site[0] != ifo[0] )
      err_exit( "site field not consitent with ifo for file %s", filename );
    if ( strcmp( fields[arg].extension, "txt" ) )
      err_exit( "incorrect extention for file %s", filename );

    /* set start and end time */
    if ( arg == 1 )
    {
      tstart = fields[arg].tstart;
      tend   = fields[arg].tstart + fields[arg].duration;
    }
    else
    {
      if ( fields[arg].tstart < tstart )
        tstart = fields[arg].tstart;
      if ( fields[arg].tstart + fields[arg].duration > tend )
        tend = fields[arg].tstart + fields[arg].duration;
    }
  }


  if ( ! have_fac )
    err_exit( "missing a calibration factors file" );
  if ( ! have_oloop_gain )
    err_exit( "missing an open loop gain reference function file" );
  if ( ! ( have_cav_gain_as_q && have_response_as_q ) &&
      ! ( have_cav_gain_darm_err && have_response_darm_err ) )
    err_exit( "missing a complete set of reference function files" );

  fprintf( stderr, "\nCalibration:\n\n" );
  fprintf( stderr, "\tInterferometer:\t%s\n", ifo );
  fprintf( stderr, "\tData Run:      \t%s\n", run );
  if ( have_cav_gain_as_q && have_response_as_q )
    fprintf( stderr, "\tChannel:       \tAS_Q\n" );
  if ( have_cav_gain_darm_err && have_response_darm_err )
    fprintf( stderr, "\tChannel:       \tDARM_ERR\n" );
  fprintf( stderr, "\tVersion:       \t%d\n", version );
  fprintf( stderr, "\n" );


  /* based on IFO name, choose the correct detector */
  if ( strstr( ifo, "H2" ) )
    detectorFlags = LAL_LHO_2K_DETECTOR_BIT;
  else if ( strstr( ifo, "H1" ) )
    detectorFlags = LAL_LHO_4K_DETECTOR_BIT;
  else if ( strstr( ifo, "L1" ) )
    detectorFlags = LAL_LLO_4K_DETECTOR_BIT;
  else
  {
    fprintf( stderr, "%s: error: %s is not a LIGO detector (don't use this program)\n", argv[0], ifo );
    return 1;
  }


  XLALGPSSetREAL8( &epoch, tstart );
  frame = XLALFrameNew( &epoch, tend - tstart, "LIGO", 0, 0, detectorFlags );
  FrHistoryAdd( frame, rcsid );
  FrHistoryAdd( frame, history );
  XLALFree( history );


  /* now read the data files */
  for ( arg = 1; arg < argc; ++arg )
  {
    const char *filename = argv[arg];
    if ( strstr( fields[arg].description, "CAL_FAC" ) ) /* factors file */
    {
      REAL4TimeSeries *alpha = NULL;
      REAL4TimeSeries *gamma = NULL;
      XLALASCIIFileReadCalFac( &alpha, &gamma, filename );
      XLALFrameAddCalFac( frame, alpha );
      XLALFrameAddCalFac( frame, gamma );
      XLALDestroyREAL4TimeSeries( gamma );
      XLALDestroyREAL4TimeSeries( alpha );
    }
    else if ( strstr( fields[arg].description, "CAL_REF" ) ) /* reference function file */
    {
      COMPLEX8FrequencySeries *series = NULL;
      REAL8 duration;
      XLALASCIIFileReadCalRef( &series, &duration, filename );
      XLALFrameAddCalRef( frame, series, version, duration );
      XLALDestroyCOMPLEX8FrequencySeries( series );
    }
    else
      err_exit( "unknown calibration file type for file %s", filename );
  }

  /* now write the frame to a framefile */
  if ( version == -1 ) /* unity factors */
    snprintf( fname, sizeof( fname ), "%c-%s_CAL_%s_VU-%d-%d.gwf", ifo[0], ifo, run, (int)floor(tstart), (int)ceil(tend) - (int)floor(tstart) );
  else
    snprintf( fname, sizeof( fname ), "%c-%s_CAL_%s_V%d-%d-%d.gwf", ifo[0], ifo,run, version, (int)floor(tstart), (int)ceil(tend) - (int)floor(tstart) );
  fprintf( stderr, "\nOutput: %s\n\n", fname );
  frfile = FrFileONew( fname, 0 );
  FrameWrite( frame, frfile );
  FrFileOEnd( frfile );

  return 0;
}
