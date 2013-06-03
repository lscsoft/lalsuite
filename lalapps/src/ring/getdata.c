/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin, Matt Pitkin
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
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/Date.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/FrameStream.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LIGOMetadataRingdownUtils.h>

#include "lalapps.h"
#include "getdata.h"
#include "errutil.h"

/* create simulated data */
REAL4TimeSeries * get_simulated_data(
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          dataType,
    REAL8        sampleRate,
    UINT4        simSeed,
    REAL4        simScale
    )
{
  RandomParams    *ranpar = NULL;
  REAL4TimeSeries *series;
  UINT4 npoints;
  UINT4 j;

  verbose( "creating simulated white gaussian noise with random seed %u\n",
      simSeed );

  series = LALCalloc( 1, sizeof( *series ) );
  if ( ! series )
    return NULL;

  npoints = duration * sampleRate;

  /* populate data with gaussian random numbers */
  series->data = XLALCreateREAL4Vector( npoints );
  ranpar = XLALCreateRandomParams( simSeed );
  XLALNormalDeviates( series->data, ranpar );
  XLALDestroyRandomParams( ranpar );

  /* set metadata */
  snprintf( series->name, sizeof( series->name ), "%s_SIM", channelName );
  series->epoch  = *epoch;
  series->deltaT = 1.0/sampleRate;
  if ( dataType == LALRINGDOWN_DATATYPE_SIM )
    series->sampleUnits = lalStrainUnit;
  else
    series->sampleUnits = lalADCCountUnit;

  /* scale series data */
  for ( j = 0; j < series->data->length; ++j )
    series->data->data[j] *= simScale;

  return series;
}

/* create simulated data */
REAL4TimeSeries * get_simulated_data_new(
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    REAL8        sampleRate,
    UINT4        simSeed,
    REAL8FrequencySeries  *psd
    )
{
  REAL8TimeSeries *series;
  REAL4TimeSeries *output;
  UINT4 npoints;
  UINT4 j;

  npoints = duration * sampleRate;

  gsl_rng *rng;
  gsl_rng_env_setup();
  rng = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set( rng, simSeed );

  verbose( "creating simulated gaussian noise with random seed %u\n",
      simSeed );

  series = XLALCreateREAL8TimeSeries("TEMP",epoch,0.,1.0/sampleRate,&lalStrainUnit,npoints);
  output = XLALCreateREAL4TimeSeries("TEMP",epoch,0.,1.0/sampleRate,&lalStrainUnit,npoints);
  snprintf( output->name, sizeof( output->name ), "%s_SIM", channelName );
  
  XLALSimNoise(series, 0 , psd, rng);

  for ( j = 0; j < series->data->length; ++j )
    output->data->data[j] = series->data->data[j];

  XLALDestroyREAL8Vector( series->data );
  LALFree( series);

  return output;
}



/* create no noise data */
REAL4TimeSeries * get_zero_data(
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          dataType,
    REAL8        sampleRate
    )
{
  REAL4TimeSeries *series;
  UINT4 npoints;
  UINT4 j;

  verbose( "creating zero data\n");

  series = LALCalloc( 1, sizeof( *series ) );
  if ( ! series )
    return NULL;

  npoints = duration * sampleRate;

  /* populate data with gaussian random numbers */
  series->data = XLALCreateREAL4Vector( npoints );

  /* set metadata */
  snprintf( series->name, sizeof( series->name ), "%s_ZERO", channelName );
  series->epoch  = *epoch;
  series->deltaT = 1.0/sampleRate;
  if ( dataType == LALRINGDOWN_DATATYPE_ZERO )
    series->sampleUnits = lalStrainUnit;
  else
    series->sampleUnits = lalADCCountUnit;
  
  for ( j = 0; j < series->data->length; ++j )
    series->data->data[j] = 0;

  return series;
}



/* read frame data */
REAL4TimeSeries * ring_get_frame_data(
    const char  *cacheName,
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          dataType
    )
{
  LALCache        *cache  = NULL;
  FrStream        *stream = NULL;
  int              mode   = LAL_FR_VERBOSE_MODE;
  REAL4TimeSeries *series;

  verbose( "get data from cache file %s\n", cacheName );

  /* open the data cache and use it to get a frame stream */
  cache  = XLALCacheImport( cacheName );
  stream = XLALFrCacheOpen( cache );
  XLALDestroyCache( cache );

  /* set the mode of the frame stream */
  XLALFrSetMode( stream, mode );

  /* read the series */
  series = fr_get_REAL4TimeSeries( channelName, epoch, duration, stream );

  /* close the frame stream */
  XLALFrClose( stream );

  /* if this is strain data, correct the units */
  if ( dataType == LALRINGDOWN_DATATYPE_HT_REAL4 )
    series->sampleUnits = lalStrainUnit;

  return series;
}


/* read double-precision frame data and convert to single-precision data */
REAL4TimeSeries * get_frame_data_dbl_convert(
    const char  *cacheName,
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          dataType,
    REAL8        dblHighPassFreq
    )
{
  REAL4TimeSeries *series;
  REAL8TimeSeries *dblser;
  LALCache        *cache  = NULL;
  FrStream        *stream = NULL;
  int              mode   = LAL_FR_VERBOSE_MODE;
  UINT4 j;

  verbose( "get data from cache file %s\n", cacheName );

  /* open the data cache and use it to get a frame stream */
  cache  = XLALCacheImport( cacheName );
  stream = XLALFrCacheOpen( cache );
  XLALDestroyCache( cache );

  /* set the mode of the frame stream */
  XLALFrSetMode( stream, mode );

  /* read double-precision series */
  dblser = fr_get_REAL8TimeSeries( channelName, epoch, duration, stream );

  /* close the frame stream */
  XLALFrClose( stream );

  /* highpass the double-precision series */
  highpass_REAL8TimeSeries( dblser, dblHighPassFreq );

  /* now create single-precision series */
  series = LALCalloc( 1, sizeof( *series ) );

  /* copy metadata */
  snprintf( series->name, sizeof( series->name ), "%s_CNVRT", dblser->name );
  series->epoch       = dblser->epoch;
  series->deltaT      = dblser->deltaT;
  series->f0          = dblser->f0;
  series->sampleUnits = dblser->sampleUnits;

  /* scale and copy data */
  series->data = XLALCreateREAL4Vector( dblser->data->length );
  for ( j = 0; j < series->data->length; ++j )
    series->data->data[j] = (REAL4)( 1.0 * dblser->data->data[j] );

  /* if this is strain data, correct the units */
  if ( dataType == LALRINGDOWN_DATATYPE_HT_REAL8 )
    series->sampleUnits = lalStrainUnit;
  
  /* destroy REAL8 time series */
  XLALDestroyREAL8Vector( dblser->data );
  LALFree(dblser);

  return series;
}


/* low-level routine to read single-precision frame data */
REAL4TimeSeries * fr_get_REAL4TimeSeries( const char *channelName,
    LIGOTimeGPS *epoch, REAL8 duration, FrStream *stream )
{
  REAL4TimeSeries *series;
  REAL8 srate;
  UINT4 npoints;

  series = LALCalloc( 1, sizeof( *series ) );
  if ( ! series )
    return NULL;

  /* if epoch is not NULL then seek to the appropriate time */
  if ( epoch )
  {
    verbose( "seek to gps time %d.%09d\n", epoch->gpsSeconds,
        epoch->gpsNanoSeconds );
    XLALFrSeek( stream, epoch );
    series->epoch = *epoch;
  }

  strncpy( series->name, channelName, sizeof( series->name ) - 1 );

  /* call first time to get sample rate */
  XLALFrGetREAL4TimeSeriesMetadata( series, stream );

  /* compute sample rate and number of points to request */
  srate   = 1.0/series->deltaT;
  npoints = floor( duration*srate + 0.5 ); /* round to nearest integer */

  /* create the data */
  series->data = XLALCreateREAL4Vector( npoints );

  /* now get the data */
  verbose( "read %g seconds of data (%u points at sample rate %g Hz)\n",
      duration, npoints, srate );
  XLALFrGetREAL4TimeSeries( series, stream );

  return series;
}


/* low-level routine to read double-precision frame data */
REAL8TimeSeries * fr_get_REAL8TimeSeries( const char *channelName,
    LIGOTimeGPS *epoch, REAL8 duration, FrStream *stream )
{
  REAL8TimeSeries *series;
  REAL8 srate;
  UINT4 npoints;

  series = LALCalloc( 1, sizeof( *series ) );
  if ( ! series )
    return NULL;

  /* if epoch is not NULL then seek to the appropriate time */
  if ( epoch )
  {
    verbose( "seek to gps time %d.%09d\n", epoch->gpsSeconds,
        epoch->gpsNanoSeconds );
    XLALFrSeek( stream, epoch );
    series->epoch = *epoch;
  }

  strncpy( series->name, channelName, sizeof( series->name ) - 1 );

  /* call first time to get sample rate */
  XLALFrGetREAL8TimeSeriesMetadata( series, stream );

  /* compute sample rate and number of points to request */
  srate   = 1.0/series->deltaT;
  npoints = floor( duration*srate + 0.5 ); /* round to nearest integer */

  /* create the data */
  series->data = XLALCreateREAL8Vector( npoints );

  /* now get the data */
  verbose( "read %g seconds of data (%u points at sample rate %g Hz)\n",
      duration, npoints, srate );
  XLALFrGetREAL8TimeSeries( series, stream );

  return series;
}


/* resample time series */
int resample_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 sampleRate )
{
  LALStatus status = blank_status;
  char name[LALNameLength];
  if ( sampleRate > 0.0 && sampleRate * series->deltaT < 1.0 )
  {
    ResampleTSParams resamplepar;
    memset( &resamplepar, 0, sizeof( resamplepar ) );
    resamplepar.deltaT     = 1.0/sampleRate;
    resamplepar.filterType = defaultButterworth;
    verbose( "resampling data from %g Hz to %g Hz\n", 1.0/series->deltaT,
        sampleRate );
    LAL_CALL( LALResampleREAL4TimeSeries( &status, series, &resamplepar ),
        &status );
    strncpy( name, series->name, LALNameLength * sizeof(char) );
    snprintf( series->name, sizeof( series->name ),
        "%s_RSMPL", name );
  }
  return 0;
}


/* highpass filter time series */
int highpass_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 frequency )
{
  LALStatus status = blank_status;
  char name[LALNameLength];
  PassBandParamStruc highpasspar;

  if ( frequency > 0.0 )
  {
    highpasspar.nMax = 8;
    highpasspar.f1   = -1;
    highpasspar.a1   = -1;
    highpasspar.f2   = frequency;
    highpasspar.a2   = 0.9; /* this means 10% attenuation at f2 */
    verbose( "highpass filtering data at %g Hz\n", frequency );
    LAL_CALL( LALDButterworthREAL4TimeSeries( &status, series, &highpasspar ),
        &status );
    strncpy( name, series->name, LALNameLength * sizeof(char) );
    snprintf( series->name, sizeof( series->name ),
        "%s_HPFLTR", name );
  }
  return 0;
}


/* highpass filter double-precision time series */
int highpass_REAL8TimeSeries( REAL8TimeSeries *series, REAL8 frequency )
{
  LALStatus status = blank_status;
  char name[LALNameLength];
  PassBandParamStruc highpasspar;
  if ( frequency > 0.0 )
  {
    highpasspar.nMax = 8;
    highpasspar.f1   = -1;
    highpasspar.a1   = -1;
    highpasspar.f2   = frequency;
    highpasspar.a2   = 0.9; /* this means 10% attenuation at f2 */
    verbose( "highpass filtering data at %g Hz\n", frequency );
    LAL_CALL( LALButterworthREAL8TimeSeries( &status, series, &highpasspar ),
        &status );
    strncpy( name, series->name, LALNameLength * sizeof(char) );
    snprintf( series->name, sizeof( series->name ),
        "%s_HPFLTR", name );
  }
  return 0;
}


/* trim padding from data */
int trimpad_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 padData )
{
  char name[LALNameLength];
  UINT4 padSamples = floor( padData / series->deltaT + 0.5 );
  UINT4 blockSamples = series->data->length - 2 * padSamples;

  if ( padData > 0 )
  {
    series = XLALResizeREAL4TimeSeries(series, padSamples, blockSamples);
    strncpy( name, series->name, LALNameLength * sizeof(char) );
    snprintf( series->name, sizeof( series->name ),
        "%s_STRIPPAD", name );
  }

  return 0;
}
