#include <math.h>
#include <string.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/AVFactories.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>
#include <lal/Calibration.h>
#include <lal/Units.h>

#include "lalapps.h"
#include "getresp.h"
#include "errutil.h"
#include "gpstime.h"

RCSID( "$Id$" );

static LALUnit      strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
static FrCacheSieve blank_sieve;


COMPLEX8FrequencySeries * get_response(
    const char  *cacheName,
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale,
    int          impulseResponse
    )
{
  COMPLEX8FrequencySeries *response;

  if ( impulseResponse )
    response = get_impulse_response( ifoName, epoch, dataDuration,
        dataSampleRate, responseScale );
  else
    response = get_frame_response( cacheName, ifoName, epoch, dataDuration,
        dataSampleRate, responseScale );

  return response;
}


COMPLEX8FrequencySeries * get_impulse_response(
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale
    )
{
  COMPLEX8FrequencySeries *response;
  UINT4 npoints;
  UINT4 k;

  response = LALCalloc( 1, sizeof( *response ) );
  if ( ! response )
    return NULL;

  npoints = floor( dataDuration * dataSampleRate + 0.5 ); /* round */
  LALSnprintf( response->name, sizeof( response->name ),
      "%s:CAL-RESPONSE", ifoName );
  response->deltaF      = 1.0/dataDuration;
  response->epoch       = *epoch;
  response->sampleUnits = strainPerCount;
  response->data        = XLALCreateCOMPLEX8Vector( npoints/2 + 1 );

  for ( k = 0; k < response->data->length; ++k )
  {
    response->data->data[k].re = responseScale;
    response->data->data[k].im = 0.0;
  }

  return response;
}

COMPLEX8FrequencySeries * get_frame_response(
    const char  *cacheName,
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale
    )
{
  LALStatus                status   = blank_status;
  FrCache                 *cache    = NULL;
  COMPLEX8FrequencySeries *response;
  COMPLEX8FrequencySeries *refResponse;
  COMPLEX8FrequencySeries *refSensing;
  COMPLEX8TimeSeries      *alpha;
  COMPLEX8TimeSeries      *alphabeta;
  CalibrationFunctions     calfuncs;
  CalibrationUpdateParams  calfacts;
  char ifo[3];
  UINT4 npoints;
  UINT4 k;

  response = LALCalloc( 1, sizeof( *response ) );
  if ( ! response )
    return NULL;

  npoints = floor( dataDuration * dataSampleRate + 0.5 ); /* round */
  LALSnprintf( response->name, sizeof( response->name ),
      "%s:CAL-RESPONSE", ifoName );
  response->deltaF      = 1.0/dataDuration;
  response->epoch       = *epoch;
  response->sampleUnits = strainPerCount;
  response->data        = XLALCreateCOMPLEX8Vector( npoints/2 + 1 );

  verbose("obtaining calibration information from cache file %s\n", cacheName);

  cache = XLALFrImportCache( cacheName );

  refResponse = get_reference_response_function( cache, ifoName );
  refSensing  = get_reference_sensing_function( cache, ifoName );
  alpha       = get_cavity_gain_factor( cache, epoch, ifoName );
  alphabeta   = get_open_loop_gain_factor( cache, epoch, ifoName );

  XLALFrDestroyCache( cache );

  strncpy( ifo, ifoName, 2 );
  ifo[2] = 0;

  calfuncs.responseFunction = refResponse;
  calfuncs.sensingFunction  = refSensing;

  calfacts.openLoopFactor   = alphabeta;
  calfacts.sensingFactor    = alpha;
  calfacts.epoch            = alpha->epoch;
  calfacts.ifo              = ifo;
  ns_to_epoch( &calfacts.duration, sec_to_ns( alpha->deltaT ) );

  verbose( "updating reference calibration to epoch %d.%09d\n",
      calfacts.epoch.gpsSeconds, calfacts.epoch.gpsNanoSeconds );
  LAL_CALL( LALUpdateCalibration( &status, &calfuncs, &calfuncs, &calfacts ),
      &status );

  verbose( "calibrating with alpha=%g and alphabeta=%g\n", calfacts.alpha.re,
      calfacts.alphabeta.re );
  LAL_CALL( LALResponseConvert( &status, response, refResponse ), &status );

  /* scale response function */
  for ( k = 0; k < response->data->length; ++k )
  {
    response->data->data[k].re *= responseScale;
    response->data->data[k].im *= responseScale;
  }

  /* cleanup memory in reference functions */
  XLALDestroyCOMPLEX8Vector( alphabeta->data );
  LALFree( alphabeta );
  XLALDestroyCOMPLEX8Vector( alpha->data );
  LALFree( alpha );
  XLALDestroyCOMPLEX8Vector( refSensing->data );
  LALFree( refSensing );
  XLALDestroyCOMPLEX8Vector( refResponse->data );
  LALFree( refResponse );

  return response;
}


COMPLEX8FrequencySeries * get_reference_response_function( FrCache *calCache,
    const char *ifoName )
{
  COMPLEX8FrequencySeries *refResponse;
  FrCache                 *refCache = NULL;
  FrCacheSieve             sieve    = blank_sieve;
  FrStream                *stream   = NULL;

  refResponse = LALCalloc( 1, sizeof( *refResponse ) );
  if ( ! refResponse )
    return NULL;

  LALSnprintf( refResponse->name, sizeof( refResponse->name ),
      "%s:CAL-RESPONSE", ifoName );

  sieve.dscRegEx = "CAL_REF";
  refCache = XLALFrSieveCache( calCache, &sieve );
  stream = XLALFrCacheOpen( refCache );
  XLALFrGetCOMPLEX8FrequencySeries( refResponse, stream );
  refResponse->sampleUnits = strainPerCount;
  XLALFrClose( stream );
  XLALFrDestroyCache( refCache );

  return refResponse;
}

COMPLEX8FrequencySeries * get_reference_sensing_function( FrCache *calCache,
    const char *ifoName )
{
  COMPLEX8FrequencySeries *refSensing;
  FrCache                 *refCache = NULL;
  FrCacheSieve             sieve    = blank_sieve;
  FrStream                *stream   = NULL;

  refSensing = LALCalloc( 1, sizeof( *refSensing ) );
  if ( ! refSensing )
    return NULL;

  LALSnprintf( refSensing->name, sizeof( refSensing->name ),
      "%s:CAL-CAV_GAIN", ifoName );

  sieve.dscRegEx = "CAL_REF";
  refCache = XLALFrSieveCache( calCache, &sieve );
  stream = XLALFrCacheOpen( refCache );
  XLALFrGetCOMPLEX8FrequencySeries( refSensing, stream );
  refSensing->sampleUnits = lalDimensionlessUnit;
  XLALFrClose( stream ); 
  XLALFrDestroyCache( refCache );

  return refSensing;
}

COMPLEX8TimeSeries * get_cavity_gain_factor( FrCache *calCache,
    LIGOTimeGPS *epoch, const char *ifoName )
{
  COMPLEX8TimeSeries *alpha;
  FrCache        *facCache = NULL;
  FrCacheSieve    sieve    = blank_sieve;
  FrStream       *stream   = NULL;

  sieve.dscRegEx = "CAL_FAC";

  alpha       = LALCalloc( 1, sizeof( *alpha ) );
  alpha->data = XLALCreateCOMPLEX8Vector( 1 );  /* only 1 point */

  LALSnprintf( alpha->name, sizeof( alpha->name ), "%s:CAL-CAV_FAC", ifoName );

  facCache = XLALFrSieveCache( calCache, &sieve );
  stream = XLALFrCacheOpen( facCache );
  XLALFrSeek( stream, epoch );
  XLALFrGetCOMPLEX8TimeSeries( alpha, stream );
  alpha->sampleUnits = lalDimensionlessUnit;
  XLALFrClose( stream );
  XLALFrDestroyCache( facCache );

  return alpha;
}

COMPLEX8TimeSeries * get_open_loop_gain_factor( FrCache *calCache,
    LIGOTimeGPS *epoch, const char *ifoName )
{
  COMPLEX8TimeSeries *alphabeta;
  FrCache        *facCache = NULL;
  FrCacheSieve    sieve    = blank_sieve;
  FrStream       *stream   = NULL;

  sieve.dscRegEx = "CAL_FAC";
  alphabeta = LALCalloc( 1, sizeof( *alphabeta ) );
  alphabeta->data = XLALCreateCOMPLEX8Vector( 1 ); /* only 1 point */

  LALSnprintf( alphabeta->name, sizeof( alphabeta->name ),
      "%s:CAL-OLOOP_FAC", ifoName );

  facCache = XLALFrSieveCache( calCache, &sieve );
  stream = XLALFrCacheOpen( facCache );
  XLALFrSeek( stream, epoch );
  XLALFrGetCOMPLEX8TimeSeries( alphabeta, stream );
  alphabeta->sampleUnits = lalDimensionlessUnit;
  XLALFrClose( stream );
  XLALFrDestroyCache( facCache );

  return alphabeta;
}

