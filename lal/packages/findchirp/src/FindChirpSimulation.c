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
#include <lal/GenerateInspiral.h>
#include <lal/GeneratePPNInspiral.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/Inject.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/GenerateInspiral.h>


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
     * generate waveform and inject it into the data
     *
     */


    /* clear the waveform structure */
    memset( &waveform, 0, sizeof(CoherentGW) );


    LALGenerateInspiral(status->statusPtr, &waveform, thisEvent, &ppnParams );
    CHECKSTATUSPTR( status );


    LALInfo( status, ppnParams.termDescription );


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

    LALSDestroyVector( status->statusPtr, &(waveform.shift->data) );
    CHECKSTATUSPTR( status );
    
    LALFree( waveform.a );
    LALFree( waveform.f );
    LALFree( waveform.phi );
    LALFree( waveform.shift );
  }

  LALCDestroyVector( status->statusPtr, &( detector.transfer->data ) );
  CHECKSTATUSPTR( status );

  if ( detector.site ) LALFree( detector.site );
  LALFree( detector.transfer );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpSimulationCP"> */
INT4
XLALFindChirpSetAnalyzeSegment (
    DataSegmentVector          *dataSegVec,
    SimInspiralTable           *injections
    )
/* </lalVerbatim> */
{  
  DataSegment      *currentSegment;
  SimInspiralTable *thisInjection;
  INT8              chanStartTime;
  INT8              chanEndTime;
  UINT4             i;
  UINT4             k;

  /* set all segments not to be analyzed by default */
  for ( i = 0; i < dataSegVec->length; ++i )
  {
    /* point to current segment */
    currentSegment = dataSegVec->data + i;
    currentSegment->analyzeSegment = 0;
  }

  /* make sure the sim inspirals are time ordered */
  XLALSortSimInspiral(&injections, XLALCompareSimInspiralByGeocentEndTime);

  /* loop over segments checking for injections into each */
  for ( i = 0; i < dataSegVec->length; ++i )
  {
    /* point to current segment */
    currentSegment = dataSegVec->data + i;

    /* compute the start and end of segment */
    chanStartTime = XLALGPStoINT8( &currentSegment->chan->epoch );
    chanEndTime = chanStartTime + 
      (INT8) (1e9 * currentSegment->chan->data->length * 
              currentSegment->chan->deltaT);

    /* look for injection into segment */
    k = 0;
    thisInjection=injections;
    while (thisInjection)
    {
      INT8 ta = XLALGPStoINT8( &thisInjection->geocent_end_time );

      if ( ta > chanStartTime && ta <= chanEndTime )
      {
        currentSegment->analyzeSegment += (UINT4)(pow(2.0,(double)(k)));
      }

      if (ta > chanEndTime)
        break;

      thisInjection=thisInjection->next;
      k = k + 1;
    }
  }

  return 0;
}


/* <lalVerbatim file="FindChirpSimulationCP"> */
void
LALFindChirpSetAnalyseTemplate (
    LALStatus                  *status,
    UINT4                      *analyseThisTmplt,
    REAL4                      mmFast,
    REAL8                      deltaF,
    INT4                       sampleRate,
    FindChirpDataParams        *fcDataParams,
    int                        numTmplts,
    InspiralTemplateNode       *tmpltHead,
    int                        numInjections,
    SimInspiralTable           *injections
    )
/* </lalVerbatim> */
{
  InspiralTemplateNode  *tmpltCurrent = NULL;
  REAL8FrequencySeries  *mmFshf = NULL;
  InspiralTemplate      *mmFTemplate = NULL;
  InspiralMetric        mmFmetric;
  InspiralMomentsEtc    mmFmoments;
  SimInspiralTable      *mmFInjection=NULL;
  UINT4                 mmF_i, ki, kj;
  REAL4                 dt0, dt3, metricDist, match;
  CHAR                  myMsg[8192];
  UINT4                 approximant;

  INITSTATUS( status, "LALFindChirpSetAnalyseTemplate", FINDCHIRPSIMULATIONC );
  ATTATCHSTATUSPTR( status );

  ASSERT( analyseThisTmplt, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcDataParams->wtildeVec->data, status, 
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  /*ASSERT( injections, status,                  */ 
  /*    FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL ); */

  /* Get the approximant. If the approximant is not BCV or BCVSpin, then */
  /* we try to tag the templates ... the BCV waveforms are not included  */
  /* yet in the scheme of things                                         */
  LALGetApproximantFromString(status->statusPtr, injections->waveform, 
      &approximant);
  CHECKSTATUSPTR (status);


  if (mmFast >= 0.0                  && 
      (approximant != (UINT4)BCV     &&
       approximant != (UINT4)BCVSpin &&
       approximant != (UINT4)BCVC )    
     )
  {
    /* If mmFast option is used, assume NONE of the templates need to be */
    /* analysed to begin with. Thus re-init analyseThisTmplt elements to */
    /* zero                                                              */
    for ( tmpltCurrent = tmpltHead, kj = 0; tmpltCurrent; 
        tmpltCurrent = tmpltCurrent->next, kj++ )
    {
      analyseThisTmplt[kj] = 0;
    }

    /* Go through the next steps only if injections is non-NULL */
    if (injections)
    {
        /* Calculate the noise moments here. This is a once and for all     */
        /* thing as the psd won't change                                    */
        mmFshf = (REAL8FrequencySeries *) 
                LALCalloc (1, sizeof(REAL8FrequencySeries));
        mmFshf->f0     = 0.0;
        mmFshf->deltaF = deltaF;
        mmFshf->data   = NULL;
        LALDCreateVector (status->statusPtr, &mmFshf->data, 
                fcDataParams->wtildeVec->length);
        CHECKSTATUSPTR (status);

        /* Populate the shf vector from the wtilde vector */
        for (ki=0; ki<fcDataParams->wtildeVec->length ; ki++) 
        {
            if (fcDataParams->wtildeVec->data[ki].re) 
            {
                mmFshf->data->data[ki] = 1./fcDataParams->wtildeVec->data[ki].re;
            }
            else 
            {
                /* Note that we can safely set shf to be zero as this is correctly */
                /* handled in the LALGetInspiralMoments function                   */
                mmFshf->data->data[ki] = 0.0;
            }
        }
        /* Init the template */
        mmFTemplate = (InspiralTemplate *) LALCalloc(1, sizeof(InspiralTemplate));
        mmFTemplate->fLower     = fcDataParams->fLow;
        mmFTemplate->fCutoff    = (REAL4)(sampleRate/2) - mmFshf->deltaF;
        mmFTemplate->tSampling  = (REAL4)(sampleRate);
        mmFTemplate->massChoice = m1Andm2;
        mmFTemplate->ieta       = 1.L;
        LALGetApproximantFromString(status->statusPtr, injections->waveform, 
                &(mmFTemplate->approximant));
        CHECKSTATUSPTR (status);
        LALGetOrderFromString(status->statusPtr, injections->waveform, 
                &(mmFTemplate->order));
        CHECKSTATUSPTR (status);

        LALSnprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
                "%d Injections, Order = %d, Approx = %d\n\n",
                numInjections, mmFTemplate->order, mmFTemplate->approximant);
        LALInfo (status, myMsg);

        mmFTemplate->mass1      = injections->mass1;
        mmFTemplate->mass2      = injections->mass2;
        LALInspiralParameterCalc( status->statusPtr, mmFTemplate );
        CHECKSTATUSPTR (status);

        LALSnprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
                "%d Injections, t0 = %e, t3 = %e\n",
                numInjections, mmFTemplate->t0, mmFTemplate->t3);
        LALInfo (status, myMsg);

        /* Now we are ready to calculate the noise moments */
        LALGetInspiralMoments( status->statusPtr, &mmFmoments, mmFshf, mmFTemplate);
        CHECKSTATUSPTR (status);

        /* We already have the noise moments so the shf series is no longer */
        /* required - free it                                               */
        LALDDestroyVector (status->statusPtr, &mmFshf->data);
        CHECKSTATUSPTR (status);
        LALFree (mmFshf);

        /* Starting with the head node of injections */
        mmFInjection = injections;

        /* Loop over all the injections */
        for (mmF_i = 0; (int)mmF_i < numInjections; mmF_i++)
        {
            /* For this injection we must                                     */
            /* (a) Calculate the metric at the point of injection             */
            /* (b) Loop over all the templates and tag theie level to 0       */
            /* or 1 with respect to this injection.                           */
            mmFTemplate->mass1      = mmFInjection->mass1;
            mmFTemplate->mass2      = mmFInjection->mass2;
            LALInspiralParameterCalc( status->statusPtr, mmFTemplate );
            CHECKSTATUSPTR (status);

            LALSnprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
                    "%d Injections, m1 = %e, m2 = %e eta = %e\n",
                    numInjections, mmFTemplate->mass1, mmFTemplate->mass2, 
                    mmFTemplate->eta);
            LALInfo (status, myMsg);
            LALSnprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
                    "%d Injections, t0 = %e, t3 = %e\n",
                    numInjections, mmFTemplate->t0, mmFTemplate->t3);
            LALInfo (status, myMsg);

            LALInspiralComputeMetric( status->statusPtr, &mmFmetric, mmFTemplate, 
                    &mmFmoments );
            CHECKSTATUSPTR (status);

            LALSnprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
                    "%d Injections, G00 = %e, G01 = %e, G11 = %e\n\n",
                    numInjections, mmFmetric.G00, mmFmetric.G01, mmFmetric.G11);
            LALInfo (status, myMsg);


            /* Now that the metric has been calculated we loop over the       */
            /* templates to mark their level.                                 */

            /* kj is the index on the templates. Note that kj always starts   */
            /* from zero even if startTemplate and stopTemplate are specified */
            /* by user while running lalapps_inspiral program                 */
            kj = 0;

            for ( tmpltCurrent = tmpltHead; tmpltCurrent; 
                    tmpltCurrent = tmpltCurrent->next)
            {
                dt0 = tmpltCurrent->tmpltPtr->t0 - mmFTemplate->t0;
                dt3 = tmpltCurrent->tmpltPtr->t3 - mmFTemplate->t3;

                metricDist  = (mmFmetric.G00 * (dt0*dt0));
                metricDist += (2.0* mmFmetric.G01 * (dt0*dt3));
                metricDist += (mmFmetric.G11 * (dt3*dt3));
                match       = 1.0 - metricDist;

                if (match >= mmFast && match <= 1.0)
                {
                    analyseThisTmplt[kj] += (UINT4)(pow(2.0, (double)(mmF_i)));
                }

                /* Advance kj for the next template */
                kj = kj + 1;

                LALSnprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
                        "%-5d %d %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
                        kj-1,
                        analyseThisTmplt[kj-1],
                        mmFTemplate->t0,
                        mmFTemplate->t3,
                        tmpltCurrent->tmpltPtr->t0,
                        tmpltCurrent->tmpltPtr->t3,
                        match,
                        dt0,
                        dt3,
                        mmFmetric.G00,
                        mmFmetric.G01,
                        mmFmetric.G11,
                        mmFmetric.theta,
                        mmFmetric.g00,
                        mmFmetric.g11
                        );
                LALInfo (status, myMsg);
            } /* End of loop over templates */

            /* Point to the next injection */
            mmFInjection = mmFInjection->next;
        }
    } /* End of if (injections) */
  }
  else 
  {
      for ( tmpltCurrent = tmpltHead, kj = 0; tmpltCurrent; 
              tmpltCurrent = tmpltCurrent->next, kj++)
      {
          analyseThisTmplt[kj] = pow(2.0, (double)(numInjections)) - 1 ;
      }
  }
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpSimulationCP"> */
UINT4
XLALCmprSgmntTmpltFlags (
        UINT4 numInjections, 
        UINT4 TmpltFlag, 
        UINT4 SgmntFlag
        )
/* </lalVerbatim> */
{
    UINT4 k1, bitTmplt, bitSgmnt, analyseTag;

    /* To begin with assume that the comparison is false. */
    analyseTag = 0;

    /* Loop over all the injections */
    for (k1=0; k1<numInjections; k1++) 
    {
        bitTmplt = (TmpltFlag>>k1)&01;
        bitSgmnt = (SgmntFlag>>k1)&01;
        if (bitSgmnt && (bitTmplt == bitSgmnt)) 
        {
            analyseTag = 1;
            break;
        }
    }

    return analyseTag;
}

