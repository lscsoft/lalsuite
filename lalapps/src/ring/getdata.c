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

static FrChanIn blank_chanin;

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
  LALStatus        status = blank_status;
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

  LAL_CALL( LALSCreateVector( &status, &series->data, npoints ), &status );
  LAL_CALL( LALCreateRandomParams( &status, &ranpar, simSeed ),
      &status );
  LAL_CALL( LALNormalDeviates( &status, series->data, ranpar ), &status );
  LAL_CALL( LALDestroyRandomParams( &status, &ranpar ), &status );

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
    int          channelType,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData
    )
{
  LALStatus        status = blank_status;
  FrCache         *cache  = NULL;
  FrStream        *stream = NULL;
  int              mode   = LAL_FR_VERBOSE_MODE;
  REAL4TimeSeries *series;

  verbose( "get data from cache file %s\n", cacheName );

  /* open the data cache and use it to get a frame stream */
  LAL_CALL( LALFrCacheImport( &status, &cache, cacheName ), &status );
  LAL_CALL( LALFrCacheOpen( &status, &stream, cache ), &status );
  LAL_CALL( LALDestroyFrCache( &status, &cache ), &status );

  /* set the mode of the frame stream */
  LAL_CALL( LALFrSetMode( &status, mode, stream ), &status );

  /* read the series */
  series = fr_get_REAL4TimeSeries( channelName, channelType, epoch, duration,
      stream );

  /* close the frame stream */
  LAL_CALL( LALFrClose( &status, &stream ), &status );

  /* if this is strain data, correct the units */
  if ( strainData )
    series->sampleUnits = lalStrainUnit;

  return series;
}

REAL4TimeSeries * get_frame_data_dbl_convert(
    const char  *cacheName,
    const char  *channelName,
    int          channelType,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData,
    REAL8        dblHighPassFreq,
    REAL8        dblScale
    )
{
  REAL4TimeSeries *series;
  REAL8TimeSeries *dblser;
  LALStatus        status = blank_status;
  FrCache         *cache  = NULL;
  FrStream        *stream = NULL;
  int              mode   = LAL_FR_VERBOSE_MODE;
  UINT4 j;

  verbose( "get data from cache file %s\n", cacheName );

  /* open the data cache and use it to get a frame stream */
  LAL_CALL( LALFrCacheImport( &status, &cache, cacheName ), &status );
  LAL_CALL( LALFrCacheOpen( &status, &stream, cache ), &status );
  LAL_CALL( LALDestroyFrCache( &status, &cache ), &status );

  /* set the mode of the frame stream */
  LAL_CALL( LALFrSetMode( &status, mode, stream ), &status );

  /* read double-precision series */
  dblser = fr_get_REAL8TimeSeries( channelName, channelType, epoch, duration,
      stream );

  /* close the frame stream */
  LAL_CALL( LALFrClose( &status, &stream ), &status );

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
  LAL_CALL( LALSCreateVector( &status, &series->data, dblser->data->length ),
      &status );
  for ( j = 0; j < series->data->length; ++j )
    series->data->data[j] = (REAL4)( dblScale * dblser->data->data[j] );

  /* if this is strain data, correct the units */
  if ( strainData )
    series->sampleUnits = lalStrainUnit;

  return series;
}


REAL4TimeSeries * fr_get_REAL4TimeSeries( const char *channelName,
    int channelType, LIGOTimeGPS *epoch, REAL8 duration, FrStream *stream )
{
  LALStatus        status = blank_status;
  FrChanIn         chanin = blank_chanin; 
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
    LAL_CALL( LALFrSeek( &status, epoch, stream ), &status );
    series->epoch = *epoch;
  }

  strncpy( series->name, channelName, sizeof( series->name ) - 1 );
  chanin.name = channelName;
  chanin.type = channelType;

  /* call first time to get sample rate */
  LAL_CALL( LALFrGetREAL4TimeSeries( &status, series, &chanin, stream ),
      &status );

  /* compute sample rate and number of points to request */
  srate   = 1.0/series->deltaT;
  npoints = floor( duration*srate + 0.5 ); /* round to nearest integer */

  /* create the data */
  LAL_CALL( LALSCreateVector( &status, &series->data, npoints ), &status );

  /* now get the data */
  verbose( "read %g seconds of data (%u points at sample rate %g Hz)\n",
      duration, npoints, srate );
  LAL_CALL( LALFrGetREAL4TimeSeries( &status, series, &chanin, stream ),
      &status );

  return series;
}


REAL8TimeSeries * fr_get_REAL8TimeSeries( const char *channelName,
    int channelType, LIGOTimeGPS *epoch, REAL8 duration, FrStream *stream )
{
  LALStatus        status = blank_status;
  FrChanIn         chanin = blank_chanin; 
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
    LAL_CALL( LALFrSeek( &status, epoch, stream ), &status );
    series->epoch = *epoch;
  }

  strncpy( series->name, channelName, sizeof( series->name ) - 1 );
  chanin.name = channelName;
  chanin.type = channelType;

  /* call first time to get sample rate */
  LAL_CALL( LALFrGetREAL8TimeSeries( &status, series, &chanin, stream ),
      &status );

  /* compute sample rate and number of points to request */
  srate   = 1.0/series->deltaT;
  npoints = floor( duration*srate + 0.5 ); /* round to nearest integer */

  /* create the data */
  LAL_CALL( LALDCreateVector( &status, &series->data, npoints ), &status );

  /* now get the data */
  verbose( "read %g seconds of data (%u points at sample rate %g Hz)\n",
      duration, npoints, srate );
  LAL_CALL( LALFrGetREAL8TimeSeries( &status, series, &chanin, stream ),
      &status );

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
