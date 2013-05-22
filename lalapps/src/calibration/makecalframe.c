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
#include <lal/Units.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

#include <lal/LALFrameIO.h>
#include <lal/LALString.h>

#include "LALASCIIFileRead.h"

#define OPEN_ENDED_DURATION (999999999)
#define NUM_OPEN_ENDED_POINTS (3)
#define OPEN_ENDED_DELTA_T (OPEN_ENDED_DURATION/NUM_OPEN_ENDED_POINTS)

const char *program;
static void err_exit( const char *fmt, ... )
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
  /** \deprecated FIXME: the following code uses obsolete CVS ID tags.
   *  It should be modified to use git version information. */
  char rcsid[] = "$Id$";
  char rcsname[] = "$Name$";
  char *arguments = NULL;
  LALDataFileNameFields *fields = NULL;
  char fname[FILENAME_MAX] = "";
  char run[FILENAME_MAX] = "";
  char ifo[3] = "";
  int version = 9999;
  int tstart = 0;
  int tend = 0;
  int tstartunity = 0;
  int open_ended_reference = 0;
  int open_ended_factors = 0;
  int have_fac = 0;
  int have_oloop_gain = 0;
  int have_actuation  = 0;
  int have_cav_gain_as_q = 0;
  int have_response_as_q = 0;
  int have_digflt_as_q   = 0;
  int have_cav_gain_darm_err = 0;
  int have_response_darm_err = 0;
  int have_digflt_darm_err   = 0;
  int detectorFlags;
  LIGOTimeGPS epoch;
  FrFile *frfile;
  FrameH *frame;
  int arg;

  program = argv[0];
  XLALSetErrorHandler( XLALAbortErrorHandler );

  arguments = XLALStringAppend( arguments, argv[0] );
  for ( arg = 1; arg < argc; ++arg )
  {
    arguments = XLALStringAppend( arguments, " " );
    arguments = XLALStringAppend( arguments, argv[arg] );
  }

  fprintf( stderr, "Input files:\n\n" );

  /* first pass: just parse filenames for sanity */
  fields = XLALCalloc( argc, sizeof( *fields ) );
  for ( arg = 1; arg < argc; ++arg )
  {
    const char *filename = argv[arg];
    /* skip all files that are not _CAL_REF_ or _CAL_FAC_ */
    if ( ! ( strstr(filename, "_CAL_REF_") || strstr(filename, "_CAL_FAC_") ) )
    {
      fprintf( stderr, "\tignoring: %s\n", filename );
      /* ignore this by setting argv[] to NULL */
      argv[arg] = NULL;
      continue;
    }

    fprintf( stderr, "\tusing: %s\n", filename );

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
      else if ( ! strcmp( dscfields.channelPostfix, "ACTUATION" ) )
        ++have_actuation;
      else if ( ! strcmp( dscfields.channelPostfix, "CAV_GAIN_AS_Q" ) )
        ++have_cav_gain_as_q;
      else if ( ! strcmp( dscfields.channelPostfix, "RESPONSE_AS_Q" ) )
        ++have_response_as_q;
      else if ( ! strcmp( dscfields.channelPostfix, "DIGFLT_AS_Q" ) )
        ++have_digflt_as_q;
      else if ( ! strcmp( dscfields.channelPostfix, "CAV_GAIN_DARM_ERR" ) )
        ++have_cav_gain_darm_err;
      else if ( ! strcmp( dscfields.channelPostfix, "RESPONSE_DARM_ERR" ) )
        ++have_response_darm_err;
      else if ( ! strcmp( dscfields.channelPostfix, "DIGFLT_DARM_ERR" ) )
        ++have_digflt_darm_err;
      else
      {
        fprintf( stderr, "%s warning: unrecognized reference function for file %s\n", program, filename );
        /* ignore this by setting argv[] to NULL */
        argv[arg] = NULL;
        continue;
      }
      if ( fields[arg].duration == OPEN_ENDED_DURATION )
        open_ended_reference = 1;
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

      /* set start and end time */
      if ( ! have_fac++ )
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
      if ( fields[arg].duration == OPEN_ENDED_DURATION )
        open_ended_factors = 1;
    }
    else
      err_exit( "invalid file name %s", filename );

    /* other sanity checks */
    if ( fields[arg].site[0] != ifo[0] )
      err_exit( "site field not consitent with ifo for file %s", filename );
    if ( strcmp( fields[arg].extension, "txt" ) )
      err_exit( "incorrect extention for file %s", filename );
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


  tstartunity = tend;
  if ( open_ended_reference || open_ended_factors )
    tend = tstart + OPEN_ENDED_DURATION;

  XLALGPSSetREAL8( &epoch, tstart );
  frame = XLALFrameNew( &epoch, tend - tstart, "LIGO", 0, 0, detectorFlags );
  XLALFrHistoryAdd( frame, "Program Id", rcsid );
  XLALFrHistoryAdd( frame, "Program Tag", rcsname );
  XLALFrHistoryAdd( frame, "Program Args", arguments );
  XLALFree( arguments );


  /* now read the data files */
  for ( arg = 1; arg < argc; ++arg )
  {
    const char *filename = argv[arg];
    if ( ! filename ) /* skip this one */
      continue;
    if ( strstr( fields[arg].description, "CAL_FAC" ) ) /* factors file */
    {
#if 0
      CHAR factorStatusName[] = "Xn:CAL-FAC_STATUS";
      CHAR seconds[] = "s";
      CHAR dimensionless[] = "";
#endif
      REAL4TimeSeries *alpha = NULL;
      REAL4TimeSeries *lal_gamma = NULL;
#if 0
      INT4TimeSeries *factorStatus = NULL;
      FrProcData *proc;
      FrVect *vect;
      UINT4 i;
#endif
      XLALASCIIFileReadCalFac( &alpha, &lal_gamma, filename );
#if 0
      memcpy( factorStatusName, ifo, 2 );
      factorStatus = XLALCreateINT4TimeSeries( factorStatusName, &alpha->epoch, alpha->f0, alpha->deltaT, &lalDimensionlessUnit, alpha->data->length );
      for ( i = 0; i < alpha->data->length; ++i )
        if ( alpha->data->data[i] == 0.0 && lal_gamma->data->data[i] == 0.0 )
          factorStatus->data->data[i] = -1; /* indicates invalid calibration */
        else
          factorStatus->data->data[i] = 0;
      if ( open_ended ) /* append unity */
      {
        UINT4 newLength;
        UINT4 origLength;
        origLength = alpha->data->length;
        newLength = floor( (REAL8)OPEN_ENDED_DURATION / alpha->deltaT );
        XLALResizeREAL4TimeSeries( alpha, 0, newLength );
        XLALResizeREAL4TimeSeries( lal_gamma, 0, newLength );
        XLALResizeINT4TimeSeries( factorStatus, 0, newLength );
        for ( i = origLength; i < newLength; ++i )
        {
          alpha->data->data[i] = 1.0;
          lal_gamma->data->data[i] = 1.0;
          factorStatus->data->data[i] = 1; /* indicates unity factor */
        }
      }
#endif
      /*
      XLALFrameAddCalFac( frame, alpha );
      XLALFrameAddCalFac( frame, lal_gamma );
      */
      XLALFrameAddCalFac( frame, alpha, version );
      XLALFrameAddCalFac( frame, lal_gamma, version );
      if ( open_ended_reference || ! open_ended_factors ) /* add stat data with unity factors */
      {
        LIGOTimeGPS unityEpoch;
        REAL4TimeSeries *unityAlpha;
        REAL4TimeSeries *unityGamma;
        UINT4 i;
        open_ended_factors = 1; /* only do this once! */
        unityEpoch.gpsSeconds = tstartunity;
        unityEpoch.gpsNanoSeconds = 0;
        unityAlpha = XLALCreateREAL4TimeSeries( alpha->name, &unityEpoch, alpha->f0, OPEN_ENDED_DELTA_T, &lalDimensionlessUnit, NUM_OPEN_ENDED_POINTS );
        unityGamma = XLALCreateREAL4TimeSeries( lal_gamma->name, &unityEpoch, lal_gamma->f0, OPEN_ENDED_DELTA_T, &lalDimensionlessUnit, NUM_OPEN_ENDED_POINTS );
        for ( i = 0; i < NUM_OPEN_ENDED_POINTS; ++i )
          unityAlpha->data->data[i] = unityGamma->data->data[i] = 1.0;
        /*
        XLALFrameAddCalFac( frame, unityAlpha, version );
        XLALFrameAddCalFac( frame, unityGamma, version );
        */
        /* REVISION: version for unity factors is 0 */
        XLALFrameAddCalFac( frame, unityAlpha, 0 );
        XLALFrameAddCalFac( frame, unityGamma, 0 );
        XLALDestroyREAL4TimeSeries( unityGamma );
        XLALDestroyREAL4TimeSeries( unityAlpha );
      }
      /* add alphaStatus channel */
#if 0
      vect = FrVectNew1D( factorStatusName, FR_VECT_4S, factorStatus->data->length, factorStatus->deltaT, seconds, dimensionless );
      vect->startX[0] = 0.0;
      memcpy( vect->data, factorStatus->data->data, factorStatus->data->length * sizeof( *factorStatus->data->data ) );
      proc = FrProcDataNewV( frame, vect );
      proc->timeOffset = 0;
      proc->type       = 1;
      proc->subType    = 0;
      proc->tRange     = factorStatus->data->length * factorStatus->deltaT;
      proc->fShift     = 0.0;
      proc->phase      = 0.0;
      proc->BW         = 0.0;
      XLALDestroyINT4TimeSeries( factorStatus );
#endif
      XLALDestroyREAL4TimeSeries( lal_gamma );
      XLALDestroyREAL4TimeSeries( alpha );
    }
    else if ( strstr( fields[arg].description, "CAL_REF" ) ) /* reference function file */
    {
      COMPLEX8FrequencySeries *series = NULL;
      REAL8 duration;
      XLALASCIIFileReadCalRef( &series, &duration, filename );
      /* modify duration so that the series fits into existing frame */
      /*
      if ( XLALGPSGetREAL8( &series->epoch ) + duration > tend )
        duration = floor( tend - XLALGPSGetREAL8( &series->epoch ) );
        */
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
  frfile = FrFileONew( fname, 1 ); /* 1 = GZIP */
  FrameWrite( frame, frfile );
  FrFileOEnd( frfile );

  return 0;
}
