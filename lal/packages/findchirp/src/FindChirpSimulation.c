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

#if 0 
<lalVerbatim file="FindChirpSimulationCV">
Author: Brown, D. A. and Creighton, T. D
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpSimulation.c}}
\label{ss:FindChirpSimulation.c}

\noindent Provides an interface between code build from \texttt{findchirp} and
various simulation packages for injecting chirps into data.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpSimulationCP}
\idx{LALFindChirpInjectSignals()}
\idx{LALRandomPPNParamStruc()}

\begin{description}
\item[\texttt{LALFindChirpInjectSignals()}] injects the signals described
in the linked list of \texttt{SimInspiralTable} structures \texttt{events}
into the data \texttt{chan}. The response function \texttt{resp} should
contain the response function to use when injecting the signals into the data.
\end{description}

\subsubsection*{Algorithm}

\noindent None.

\subsubsection*{Notes}
\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpSimulationCV}}
</lalLaTeX> 
#endif

#include <lal/Units.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/VectorOps.h>
#include <lal/SeqFactories.h>
#include <lal/DetectorSite.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/Inject.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>

NRCSID( FINDCHIRPSIMULATIONC, "$Id$" );

/* <lalVerbatim file="FindChirpSimulationCP"> */
void
LALFindChirpInjectSignals (
    LALStatus                  *status,
    REAL4TimeSeries            *chan,
    SimInspiralTable           *events,
    COMPLEX8FrequencySeries    *resp
    )
/* </lalVerbatim> */
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
  CHAR                  warnMsg[512];

  INITSTATUS( status, "LALFindChirpInjectSignals", FINDCHIRPSIMULATIONC );
  ATTATCHSTATUSPTR( status );

  ASSERT( chan, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( chan->data, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( chan->data->data, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  ASSERT( events, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  ASSERT( resp, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( resp->data, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( resp->data->data, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );


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
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
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
    case 'G':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      LALWarning( status, "computing waveform for GEO600." );
      break;
    case 'T':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      LALWarning( status, "computing waveform for TAMA300." );
      break;
    case 'V':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
      LALWarning( status, "computing waveform for Virgo." );
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
    ppnParams.mTot = thisEvent->mass1 + thisEvent->mass2;
    ppnParams.eta  = thisEvent->eta;
    ppnParams.d    = thisEvent->distance * 1.0e6 * LAL_PC_SI;
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

    LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
        "Injected waveform parameters:\n"
        "ppnParams.mTot = %e\n"
        "ppnParams.eta = %e\n"
        "ppnParams.d = %e\n"
        "ppnParams.inc = %e\n"
        "ppnParams.phi = %e\n"
        "ppnParams.psi = %e\n"
        "ppnParams.position.longitude = %e\n"
        "ppnParams.position.latitude = %e\n", 
        ppnParams.mTot, 
        ppnParams.eta, 
        ppnParams.d,
        ppnParams.inc,
        ppnParams.phi,
        ppnParams.psi, 
        ppnParams.position.longitude, 
        ppnParams.position.latitude );
    LALInfo( status, warnMsg );


    /* 
     *
     * generate waveform and inject it into the data
     *
     */


    /* clear the waveform structure */
    memset( &waveform, 0, sizeof(CoherentGW) );

    if ( !strncmp("GeneratePPN",thisEvent->waveform,11))
    {
      /* generate the waveform for injection */
      LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
          "Injecting waveform using LALGeneratePPNInspiral\n");
      LALInfo( status, warnMsg );
      LALGeneratePPNInspiral( status->statusPtr, &waveform, &ppnParams );
      CHECKSTATUSPTR( status );
    }
    else if ( !strncmp("EOB",thisEvent->waveform,3) ||
              !strncmp("PadeT1",thisEvent->waveform,6) ||
              !strncmp("TaylorT1",thisEvent->waveform,8) ||
              !strncmp("TaylorT2",thisEvent->waveform,8) ||
              !strncmp("TaylorT3",thisEvent->waveform,8) ||
              !strncmp("SpinTaylor",thisEvent->waveform,10) )
    {
      /* create the waveform using the new function from inject */
      LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
          "Injecting waveform using LALGenerateInspiral\n");
      LALInfo( status, warnMsg );
      LALGenerateInspiral(status->statusPtr, &waveform, thisEvent, &ppnParams );
      CHECKSTATUSPTR( status );
    }
    else
    {
      LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
          "Unknown injection waveform\n");
      LALInfo( status, warnMsg );
      ABORT( status, FINDCHIRPH_EWVFM, FINDCHIRPH_MSGEWVFM );
    }

    LALInfo( status, ppnParams.termDescription );

    if ( ppnParams.dfdt > 2.0 ) 
    {
      LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
          "Waveform sampling interval is too large:\n"
          "\tmaximum df*dt = %f", ppnParams.dfdt );
      LALInfo( status, warnMsg );
      ABORT( status, FINDCHIRPH_EDFDT, FINDCHIRPH_MSGEDFDT );
    }

    if ( thisEvent->geocent_end_time.gpsSeconds )
    {
      /* get the gps start time of the signal to inject */
      LALGPStoINT8( status->statusPtr, &waveformStartTime, 
          &(thisEvent->geocent_end_time) );
      CHECKSTATUSPTR( status );
      waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );
    }
    else
    {
      /* center the waveform in the data segment */
      LALGPStoINT8( status->statusPtr, &waveformStartTime, 
          &(chan->epoch) );
      CHECKSTATUSPTR( status );

      waveformStartTime += (INT8) ( 1000000000.0 * 
          ((REAL8) (chan->data->length - ppnParams.length) / 2) * chan->deltaT
          );
    }

    LALSnprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg), 
        "Injected waveform timing:\n"
        "thisEvent->geocent_end_time.gpsSeconds = %d\n"
        "thisEvent->geocent_end_time.gpsNanoSeconds = %d\n"
        "ppnParams.tc = %e\n"
        "waveformStartTime = %lld\n",
        thisEvent->geocent_end_time.gpsSeconds,
        thisEvent->geocent_end_time.gpsNanoSeconds,
        ppnParams.tc,
        waveformStartTime );
    LALInfo( status, warnMsg );

    /* clear the signal structure */
    memset( &signal, 0, sizeof(REAL4TimeSeries) );

    /* set the start times for injection */
    LALINT8toGPS( status->statusPtr, &(waveform.a->epoch), &waveformStartTime );
    CHECKSTATUSPTR( status );
    memcpy( &(waveform.f->epoch), &(waveform.a->epoch), 
        sizeof(LIGOTimeGPS) );
    memcpy( &(waveform.phi->epoch), &(waveform.a->epoch), 
        sizeof(LIGOTimeGPS) );

    /* set the start time of the signal vector to the start time of the chan */
    signal.epoch = chan->epoch;

    /* set the parameters for the signal time series */
    signal.deltaT = chan->deltaT;
    if ( ( signal.f0 = chan->f0 ) != 0 )
    {
      ABORT( status, FINDCHIRPH_EHETR, FINDCHIRPH_MSGEHETR );
    }
    signal.sampleUnits = lalADCCountUnit;

    /* simulate the detectors response to the inspiral */
    LALSCreateVector( status->statusPtr, &(signal.data), chan->data->length );
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
