/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSimulation.c
 *
 * Author: Brown, D. A., and Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <lal/SeqFactories.h>
#include <lal/DetectorSite.h>
#include <lal/LALNoiseModels.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpEngine.h>

NRCSID( FINDCHIRPSIMULATIONC, "$Id$" );

void
LALFindChirpInjectSignals (
    LALStatus                  *status,
    REAL4TimeSeries            *chan,
    SimInspiralTable           *events,
    COMPLEX8FrequencySeries    *resp
    )
{
  UINT4                 k;
  DetectorResponse      detector;
  SimInspiralTable     *thisEvent = NULL;
  PPNParamStruc         ppnParams;
  CoherentGW            waveform;
  INT8                  waveformStartTime;
  INT8                  chanStartTime;
  REAL4TimeSeries       signal;
  COMPLEX8Vector       *unity = NULL;


  INITSTATUS( status, "LALFindChirpInjectSignals", FINDCHIRPSIMULATIONC );
  ATTATCHSTATUSPTR( status );

  ASSERT( chan, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( chan->data, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( chan->data->data, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );

  ASSERT( events, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );

  ASSERT( resp, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( resp->data, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( resp->data->data, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );


  /*
   *
   * set up structures and parameters needed
   *
   */


  LALGPStoINT8( status->statusPtr, &chanStartTime, &(chan->epoch) );
  CHECKSTATUSPTR( status );

  /* fixed waveform injection parameters */
  ppnParams.deltaT   = chan->deltaT;
  ppnParams.lengthIn = 0;
  ppnParams.ppn      = NULL;


  /*
   *
   * compute the transfer function from the given response function
   *
   */


  /* allocate memory and copy the parameters describing the freq series */
  memset( &detector, 0, sizeof( DetectorResponse ) );
  detector.transfer = (COMPLEX8FrequencySeries *)
    LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
  if ( ! detector.transfer ) 
  {
    ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
  }
  memcpy( &(detector.transfer->epoch), &(resp->epoch),
      sizeof(LIGOTimeGPS) );
  detector.transfer->f0 = resp->f0;
  detector.transfer->deltaF = resp->deltaF;

  detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );
  /* set the detector site */
  switch ( chan->name[0] )
  {
    case 'H':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexLHODIFF];
      LALWarning( status, "computing waveform for Hanford." );
      break;
    case 'L':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
      LALWarning( status, "computing waveform for Livingston." );
      break;
    default:
      LALFree( detector.site );
      detector.site = NULL;
      LALWarning( status, "Unknown detector site, computing plus mode "
          "waveform with no time delay" );
      break;
  }

  /* set up units for the transfer function */
  {
    RAT4 negOne = { -1, 0 };
    LALUnit unit;
    LALUnitPair pair;
    pair.unitOne = &lalADCCountUnit;
    pair.unitTwo = &lalStrainUnit;
    LALUnitRaise( status->statusPtr, &unit, pair.unitTwo, &negOne );
    CHECKSTATUSPTR( status );
    pair.unitTwo = &unit;
    LALUnitMultiply( status->statusPtr, &(detector.transfer->sampleUnits),
        &pair );
    CHECKSTATUSPTR( status );
  }

  /* invert the response function to get the transfer function */
  LALCCreateVector( status->statusPtr, &( detector.transfer->data ),
      resp->data->length );
  CHECKSTATUSPTR( status );

  LALCCreateVector( status->statusPtr, &unity, resp->data->length );
  CHECKSTATUSPTR( status );
  for ( k = 0; k < resp->data->length; ++k ) 
  {
    unity->data[k].re = 1.0;
    unity->data[k].im = 0.0;
  }

  LALCCVectorDivide( status->statusPtr, detector.transfer->data, unity,
      resp->data );
  CHECKSTATUSPTR( status );

  LALCDestroyVector( status->statusPtr, &unity );
  CHECKSTATUSPTR( status );


  /*
   *
   * loop over the signals and inject them into the time series
   *
   */


  for ( thisEvent = events; thisEvent; thisEvent = thisEvent->next )
  {


    /* 
     *
     * populate ppn parameter structure from injection event 
     *
     */


    /* input fields */
    ppnParams.mTot = thisEvent->mtotal;
    ppnParams.eta  = thisEvent->eta;
    ppnParams.d    = thisEvent->distance;
    ppnParams.inc  = thisEvent->inclination;
    ppnParams.phi  = thisEvent->coa_phase;

    /* frequency cutoffs */
    ppnParams.fStartIn = 40.0;
    ppnParams.fStopIn  = -1.0 / 
      (6.0 * sqrt(6.0) * LAL_PI * ppnParams.mTot * LAL_MTSUN_SI);

    /* passed fields */
    ppnParams.position.longitude   = thisEvent->longitude;
    ppnParams.position.latitude    = thisEvent->latitude;
    ppnParams.position.system      = COORDINATESYSTEM_EQUATORIAL;
    ppnParams.psi                  = thisEvent->polarization;
    ppnParams.epoch.gpsSeconds     = 0;
    ppnParams.epoch.gpsNanoSeconds = 0;


    /* 
     *
     * generate waveform and inject it into the data
     *
     */


    /* clear the waveform structure */
    memset( &waveform, 0, sizeof(CoherentGW) );

    LALGeneratePPNInspiral( status->statusPtr, &waveform, &ppnParams );
    CHECKSTATUSPTR( status );

    if ( ppnParams.dfdt > 2.0 ) 
    {
      fprintf( stderr, "Waveform sampling interval is too large:\n"
          "\tmaximum df*dt = %f", ppnParams.dfdt );
      fflush( stderr );
    }

    /* get the gps start time of the signal to inject */
    LALGPStoINT8( status->statusPtr, &waveformStartTime, 
        &(thisEvent->geocent_end_time) );
    CHECKSTATUSPTR( status );
    waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );

    /* clear the signal structure */
    memset( &signal, 0, sizeof(REAL4TimeSeries) );

    /* set the start times for injection */
    LALINT8toGPS( status->statusPtr, &(signal.epoch), &waveformStartTime );
    CHECKSTATUSPTR( status );
    memcpy( &(waveform.a->epoch), &(signal.epoch), 
        sizeof(LIGOTimeGPS) );
    memcpy( &(waveform.f->epoch), &(signal.epoch), 
        sizeof(LIGOTimeGPS) );
    memcpy( &(waveform.phi->epoch), &(signal.epoch), 
        sizeof(LIGOTimeGPS) );

    /* set the parameters for the signal time series */
    signal.deltaT = chan->deltaT;
    if ( ( signal.f0 = chan->f0 ) != 0 )
    {
      ABORT( status, FINDCHIRPENGINEH_EHETR, FINDCHIRPENGINEH_MSGEHETR );
    }
    signal.sampleUnits = lalADCCountUnit;

    /* simulate the detectors response to the inspiral */
    LALSCreateVector( status->statusPtr, &(signal.data), 
        (UINT4) ppnParams.length );
    CHECKSTATUSPTR( status );

    LALSimulateCoherentGW( status->statusPtr, 
        &signal, &waveform, &detector );
    CHECKSTATUSPTR( status );

    /* inject the signal into the data channel */
    LALSSInjectTimeSeries( status->statusPtr, chan, &signal );
    CHECKSTATUSPTR( status );

    /* destroy the signal */
    LALSDestroyVector( status->statusPtr, &(signal.data) );
    CHECKSTATUSPTR( status );

    LALSDestroyVectorSequence( status->statusPtr, &(waveform.a->data) );
    CHECKSTATUSPTR( status );

    LALSDestroyVector( status->statusPtr, &(waveform.f->data) );
    CHECKSTATUSPTR( status );

    LALDDestroyVector( status->statusPtr, &(waveform.phi->data) );
    CHECKSTATUSPTR( status );

    LALFree( waveform.a );
    LALFree( waveform.f );
    LALFree( waveform.phi );
  }

  LALCDestroyVector( status->statusPtr, &( detector.transfer->data ) );
  CHECKSTATUSPTR( status );

  if ( detector.site ) LALFree( detector.site );
  LALFree( detector.transfer );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
