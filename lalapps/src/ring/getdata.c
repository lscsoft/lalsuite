#include <math.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/FrameStream.h>
#include <lal/Units.h>

#include "lalapps.h"
#include "getdata.h"
#include "errutil.h"

RCSID( "$Id$" );

REAL4TimeSeries * get_simulated_data(
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData,
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

  series->data = XLALCreateREAL4Vector( npoints );
  ranpar = XLALCreateRandomParams( simSeed );
  XLALNormalDeviates( series->data, ranpar );
  XLALDestroyRandomParams( ranpar );

  /* set metadata */
  LALSnprintf( series->name, sizeof( series->name ), "%s_SIM", channelName );
  series->epoch  = *epoch;
  series->deltaT = 1.0/sampleRate;
  if ( strainData )
    series->sampleUnits = lalStrainUnit;
  else
    series->sampleUnits = lalADCCountUnit;

  /* scale series data */
  for ( j = 0; j < series->data->length; ++j )
    series->data->data[j] *= simScale;

  return series;
}

REAL4TimeSeries * get_frame_data(
    const char  *cacheName,
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData
    )
{
  FrCache         *cache  = NULL;
  FrStream        *stream = NULL;
  int              mode   = LAL_FR_VERBOSE_MODE;
  REAL4TimeSeries *series;

  verbose( "get data from cache file %s\n", cacheName );

  /* open the data cache and use it to get a frame stream */
  cache  = XLALFrImportCache( cacheName );
  stream = XLALFrCacheOpen( cache );
  XLALFrDestroyCache( cache );

  /* set the mode of the frame stream */
  XLALFrSetMode( stream, mode );

  /* read the series */
  series = fr_get_REAL4TimeSeries( channelName, epoch, duration, stream );

  /* close the frame stream */
  XLALFrClose( stream );

  /* if this is strain data, correct the units */
  if ( strainData )
    series->sampleUnits = lalStrainUnit;

  return series;
}

REAL4TimeSeries * get_frame_data_dbl_convert(
    const char  *cacheName,
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData,
    REAL8        dblHighPassFreq,
    REAL8        dblScale
    )
{
  REAL4TimeSeries *series;
  REAL8TimeSeries *dblser;
  FrCache         *cache  = NULL;
  FrStream        *stream = NULL;
  int              mode   = LAL_FR_VERBOSE_MODE;
  UINT4 j;

  verbose( "get data from cache file %s\n", cacheName );

  /* open the data cache and use it to get a frame stream */
  cache  = XLALFrImportCache( cacheName );
  stream = XLALFrCacheOpen( cache );
  XLALFrDestroyCache( cache );

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
  LALSnprintf( series->name, sizeof( series->name ), "%s_CNVRT", dblser->name );
  series->epoch       = dblser->epoch;
  series->deltaT      = dblser->deltaT;
  series->f0          = dblser->f0;
  series->sampleUnits = dblser->sampleUnits;

  /* scale and copy data */
  series->data = XLALCreateREAL4Vector( dblser->data->length );
  for ( j = 0; j < series->data->length; ++j )
    series->data->data[j] = (REAL4)( dblScale * dblser->data->data[j] );

  /* if this is strain data, correct the units */
  if ( strainData )
    series->sampleUnits = lalStrainUnit;

  return series;
}


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


int resample_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 sampleRate )
{
  LALStatus status = blank_status;
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
    LALSnprintf( series->name, sizeof( series->name ),
        "%s_RSMPL", series->name );
  }
  return 0;
}


int highpass_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 frequency )
{
  LALStatus status = blank_status;
  PassBandParamStruc highpasspar;
  if ( frequency > 0.0 )
  {
    highpasspar.nMax = 8;
    highpasspar.f1   = -1;
    highpasspar.a1   = -1;
    highpasspar.f2   = frequency;
    highpasspar.a1   = 0.9; /* this means 10% attenuation at f2 */
    verbose( "highpass filtering data at %g Hz\n", frequency );
    LAL_CALL( LALDButterworthREAL4TimeSeries( &status, series, &highpasspar ),
        &status );
    LALSnprintf( series->name, sizeof( series->name ),
        "%s_HPFLTR", series->name );
  }
  return 0;
}

int highpass_REAL8TimeSeries( REAL8TimeSeries *series, REAL8 frequency )
{
  LALStatus status = blank_status;
  PassBandParamStruc highpasspar;
  if ( frequency > 0.0 )
  {
    highpasspar.nMax = 8;
    highpasspar.f1   = -1;
    highpasspar.a1   = -1;
    highpasspar.f2   = frequency;
    highpasspar.a1   = 0.9; /* this means 10% attenuation at f2 */
    verbose( "highpass filtering data at %g Hz\n", frequency );
    LAL_CALL( LALButterworthREAL8TimeSeries( &status, series, &highpasspar ),
        &status );
    LALSnprintf( series->name, sizeof( series->name ),
        "%s_HPFLTR", series->name );
  }
  return 0;
}
