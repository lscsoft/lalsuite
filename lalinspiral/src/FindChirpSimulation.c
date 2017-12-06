/*
*  Copyright (C) 2007 Duncan Brown, Eirini Messaritaki, Gareth Jones, Jolien Creighton, Patrick Brady, Reinhard Prix, Anand Sengupta, Stephen Fairhurst, Craig Robinson , Thomas Cokelaer
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

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpSimulation.c
 *
 * Author: Brown, D. A., and Creighton, T. D.
 *
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A. and Creighton, T. D
 * \file
 * \ingroup FindChirp_h
 *
 * \brief Provides an interface between code build from \c findchirp and
 * various simulation packages for injecting chirps into data.
 *
 * ### Prototypes ###
 *
 * <dl>
 * <dt><tt>LALFindChirpInjectSignals()</tt></dt><dd> injects the signals described
 * in the linked list of \c SimInspiralTable structures \c events
 * into the data \c chan. The response function \c resp should
 * contain the response function to use when injecting the signals into the data.</dd>
 * </dl>
 *
 * ### Algorithm ###
 *
 * None.
 *
 * ### Notes ###
 *
 *
 * ### Uses ###
 *
 * \code
 * LALCalloc()
 * LALFree()
 * \endcode
 *
 * ### Notes ###
 *
 */

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
#include <lal/LIGOMetadataInspiralUtils.h>
#include <lal/LIGOMetadataRingdownUtils.h>
#include <lal/LALInspiralBank.h>
#include <lal/FindChirp.h>
#include <lal/LALStdlib.h>
#include <lal/LALInspiralBank.h>
#include <lal/GenerateInspiral.h>
#include <lal/NRWaveInject.h>
#include <lal/GenerateInspRing.h>
#include <lal/LALErrno.h>
#include <math.h>
#include <lal/LALInspiral.h>
#include <lal/LALError.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimulation.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static int FindTimeSeriesStartAndEnd (
              REAL4Vector *signalvec,
              UINT4 *start,
              UINT4 *end
             );


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
  REAL4TimeSeries       signalvec;
  COMPLEX8Vector       *unity = NULL;
  CHAR                  warnMsg[512];
  CHAR                  ifo[LIGOMETA_IFO_MAX];
  REAL8                 timeDelay;

  INITSTATUS(status);
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


  /* fixed waveform injection parameters */
  memset( &ppnParams, 0, sizeof(PPNParamStruc) );
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
      strcpy(ifo, "H1");
      break;
    case 'L':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
      LALWarning( status, "computing waveform for Livingston." );
      strcpy(ifo, "L1");
      break;
    case 'G':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
      LALWarning( status, "computing waveform for GEO600." );
      strcpy(ifo, "G1");
      break;
    case 'T':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexTAMA300DIFF];
      LALWarning( status, "computing waveform for TAMA300." );
      strcpy(ifo, "T1");
      break;
    case 'V':
      *(detector.site) = lalCachedDetectors[LALDetectorIndexVIRGODIFF];
      LALWarning( status, "computing waveform for Virgo." );
      strcpy(ifo, "V1");
      break;
    default:
      LALFree( detector.site );
      detector.site = NULL;
      LALWarning( status, "Unknown detector site, computing plus mode "
          "waveform with no time delay" );
      break;
  }

  /* set up units for the transfer function */
  if (XLALUnitDivide( &(detector.transfer->sampleUnits),
                      &lalADCCountUnit, &lalStrainUnit ) == NULL) {
    ABORTXLAL(status);
  }

  /* invert the response function to get the transfer function */
  LALCCreateVector( status->statusPtr, &( detector.transfer->data ),
      resp->data->length );
  CHECKSTATUSPTR( status );

  LALCCreateVector( status->statusPtr, &unity, resp->data->length );
  CHECKSTATUSPTR( status );
  for ( k = 0; k < resp->data->length; ++k )
  {
    unity->data[k] = 1.0;
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

    if ( strstr( thisEvent->waveform, "KludgeIMR") ||
         strstr( thisEvent->waveform, "KludgeRingOnly") )
     {
       CoherentGW *wfm;
       SimRingdownTable *ringEvent;
       int injectSignalType = LALRINGDOWN_IMR_INJECT;


       ringEvent = (SimRingdownTable *)
         LALCalloc( 1, sizeof(SimRingdownTable) );
       wfm = XLALGenerateInspRing( &waveform, thisEvent, ringEvent,
           injectSignalType);
       LALFree(ringEvent);

       if ( !wfm )
       {
         LALInfo( status, "Unable to generate merger/ringdown, "
             "injecting inspiral only");
         ABORT( status, FINDCHIRPH_EIMRW, FINDCHIRPH_MSGEIMRW );
       }
       waveform = *wfm;
     }


    if ( thisEvent->geocent_end_time.gpsSeconds )
    {
      /* get the gps start time of the signal to inject */
      waveformStartTime = XLALGPSToINT8NS( &(thisEvent->geocent_end_time) );
      waveformStartTime -= (INT8) ( 1000000000.0 * ppnParams.tc );
    }
    else
    {
      LALInfo( status, "Waveform start time is zero: injecting waveform "
          "into center of data segment" );

      /* center the waveform in the data segment */
      waveformStartTime = XLALGPSToINT8NS( &(chan->epoch) );

      waveformStartTime += (INT8) ( 1000000000.0 *
          ((REAL8) (chan->data->length - ppnParams.length) / 2.0) * chan->deltaT
          );
    }

    snprintf( warnMsg, sizeof(warnMsg)/sizeof(*warnMsg),
        "Injected waveform timing:\n"
        "thisEvent->geocent_end_time.gpsSeconds = %d\n"
        "thisEvent->geocent_end_time.gpsNanoSeconds = %d\n"
        "ppnParams.tc = %e\n"
        "waveformStartTime = %" LAL_INT8_FORMAT "\n",
        thisEvent->geocent_end_time.gpsSeconds,
        thisEvent->geocent_end_time.gpsNanoSeconds,
        ppnParams.tc,
        waveformStartTime );
    LALInfo( status, warnMsg );

    /* if we generated waveform.h then use LALInjectStrainGW
       otherwise use SimulateCoherentGW */
    if( waveform.h == NULL)
    {
      /* clear the signal structure */
      memset( &signalvec, 0, sizeof(REAL4TimeSeries) );

      /* set the start time of the signal vector to the appropriate start time of the injection */
      if ( detector.site )
      {
        timeDelay = XLALTimeDelayFromEarthCenter( detector.site->location, thisEvent->longitude,
          thisEvent->latitude, &(thisEvent->geocent_end_time) );
        if ( XLAL_IS_REAL8_FAIL_NAN( timeDelay ) )
        {
          ABORTXLAL( status );
        }
      }
      else
      {
        timeDelay = 0.0;
      }
      /* Give a little more breathing space to aid band-passing */
      XLALGPSSetREAL8( &(signalvec.epoch), (waveformStartTime * 1.0e-9) - 0.25 + timeDelay );
      UINT4 signalvecLength=waveform.phi->data->length + (UINT4)ceil((0.5+timeDelay)/waveform.phi->deltaT);
      

      /* set the parameters for the signal time series */
      signalvec.deltaT = chan->deltaT;
      if ( ( signalvec.f0 = chan->f0 ) != 0 )
      {
        ABORT( status, FINDCHIRPH_EHETR, FINDCHIRPH_MSGEHETR );
      }
      signalvec.sampleUnits = lalADCCountUnit;

      /* set the start times for injection */
      XLALINT8NSToGPS( &(waveform.a->epoch), waveformStartTime );
      /* put a rug on a polished floor? */
      waveform.f->epoch = waveform.a->epoch;
      waveform.phi->epoch = waveform.a->epoch;
      /* you might as well set a man trap */
      if ( waveform.shift )
      {
        waveform.shift->epoch = waveform.a->epoch;
      }
      /* and to think he'd just come from the hospital */

      /* simulate the detectors response to the inspiral */
      LALSCreateVector( status->statusPtr, &(signalvec.data), signalvecLength );
      CHECKSTATUSPTR( status );

      LALSimulateCoherentGW( status->statusPtr,
          &signalvec, &waveform, &detector );
      CHECKSTATUSPTR( status );
      
      /* Taper the signal */
      {

          if ( ! strcmp( "TAPER_START", thisEvent->taper ) )
          {
              XLALSimInspiralREAL4WaveTaper( signalvec.data, LAL_SIM_INSPIRAL_TAPER_START );
          }
          else if (  ! strcmp( "TAPER_END", thisEvent->taper ) )
          {
              XLALSimInspiralREAL4WaveTaper( signalvec.data, LAL_SIM_INSPIRAL_TAPER_END );
          }
          else if (  ! strcmp( "TAPER_STARTEND", thisEvent->taper ) )
          {
              XLALSimInspiralREAL4WaveTaper( signalvec.data, LAL_SIM_INSPIRAL_TAPER_STARTEND );
          }
          else if ( strcmp( "TAPER_NONE", thisEvent->taper ) )
          {
              XLALPrintError( "Invalid injection tapering option specified: %s\n",
                 thisEvent->taper );
              ABORT( status, LAL_BADPARM_ERR, LAL_BADPARM_MSG );
          }
      }

      /* Band pass the signal */
      if ( thisEvent->bandpass )
      {
          UINT4 safeToBandPass = 0;
          UINT4 start=0, end=0;
          REAL4Vector *bandpassVec = NULL;

          safeToBandPass = FindTimeSeriesStartAndEnd (
                  signalvec.data, &start, &end );

          if ( safeToBandPass )
          {
              /* Check if we can grab some padding at the extremeties.
               * This will make the bandpassing better
               */

              if (((INT4)start - (int)(0.25/chan->deltaT)) > 0 )
                    start -= (int)(0.25/chan->deltaT);
              else
                    start = 0;

              if ((end + (int)(0.25/chan->deltaT)) < signalvec.data->length )
                    end += (int)(0.25/chan->deltaT);
              else
                    end = signalvec.data->length - 1;

              bandpassVec = (REAL4Vector *)
                      LALCalloc(1, sizeof(REAL4Vector) );

              bandpassVec->length = (end - start + 1);
              bandpassVec->data = signalvec.data->data + start;

              if ( XLALBandPassInspiralTemplate( bandpassVec,
                          1.1*thisEvent->f_lower,
                          1.05*thisEvent->f_final,
                          1./chan->deltaT) != XLAL_SUCCESS )
              {
                  LALError( status, "Failed to Bandpass signal" );
                  ABORT (status, LALINSPIRALH_EBPERR, LALINSPIRALH_MSGEBPERR);
              };

              LALFree( bandpassVec );
          }
      }

      /* inject the signal into the data channel */
      int retcode=XLALSimAddInjectionREAL4TimeSeries(chan, &signalvec, NULL);

      if(retcode!=XLAL_SUCCESS){
	ABORTXLAL(status);
      }

    }
    else
    {
      const INT4 wfmLength = waveform.h->data->length;
      INT4 i = wfmLength;
      REAL4 *dataPtr = waveform.h->data->data;
      REAL4 *tmpdata;

      /* XXX This code will BREAK if the first element of the frequency XXX *
       * XXX series does not contain dynRange. This is the case for     XXX *
       * XXX calibrated strain data, but will not be the case when      XXX *
       * XXX filtering uncalibrated data.                               XXX */
      REAL8 dynRange;
      LALWarning (status, "Attempting to calculate dynRange: Will break if un-calibrated strain-data is used.");
      dynRange = 1.0/(crealf(resp->data->data[0]));

      /* set the start times for injection */
      XLALINT8NSToGPS( &(waveform.h->epoch), waveformStartTime );

      /*
       * We are using functions from NRWaveInject which do not take a vector
       * with alternating h+ & hx but rather one that stores h+ in the first
       * half and hx in the second half. That is, we must transpose h.
       *
       * We also must multiply the strain by the distance to be compatible
       * with NRWaveInject.
       */
      if (waveform.h->data->vectorLength != 2)
          LALAbort("expected alternating h+ and hx");
      tmpdata = XLALCalloc(2 * wfmLength, sizeof(*tmpdata));
      for (i=0;i<wfmLength;i++) {
            tmpdata[i] = dataPtr[2*i] * thisEvent->distance;
            tmpdata[wfmLength+i] = dataPtr[2*i+1] * thisEvent->distance;
      }
      memcpy(dataPtr, tmpdata, 2 * wfmLength * sizeof(*tmpdata));
      XLALFree(tmpdata);
      waveform.h->data->vectorLength = wfmLength;
      waveform.h->data->length = 2;

      LALInjectStrainGW( status->statusPtr ,
                                      chan ,
                                waveform.h ,
                                 thisEvent ,
                                       ifo ,
                                  dynRange  );
      CHECKSTATUSPTR( status );

    }




    if ( waveform.shift )
    {
      LALSDestroyVector( status->statusPtr, &(waveform.shift->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.shift );
    }

    if( waveform.h )
    {
      LALSDestroyVectorSequence( status->statusPtr, &(waveform.h->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.h );
    }
    if( waveform.a )
    {
      LALSDestroyVectorSequence( status->statusPtr, &(waveform.a->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.a );
      /*
       * destroy the signal only if waveform.h is NULL as otherwise it won't
       * be created
       * */
      if ( waveform.h == NULL )
      {
	LALSDestroyVector( status->statusPtr, &(signalvec.data) );
        CHECKSTATUSPTR( status );
      }
    }
    if( waveform.f )
    {
      LALSDestroyVector( status->statusPtr, &(waveform.f->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.f );
    }
    if( waveform.phi )
    {
      LALDDestroyVector( status->statusPtr, &(waveform.phi->data) );
      CHECKSTATUSPTR( status );
      LALFree( waveform.phi );
    }
  }

  LALCDestroyVector( status->statusPtr, &( detector.transfer->data ) );
  CHECKSTATUSPTR( status );

  if ( detector.site ) LALFree( detector.site );
  LALFree( detector.transfer );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}



INT4
XLALFindChirpSetAnalyzeSegment (
    DataSegmentVector          *dataSegVec,
    SimInspiralTable           *injections
    )

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
    chanStartTime = XLALGPSToINT8NS( &currentSegment->chan->epoch );
    chanEndTime = chanStartTime +
      (INT8) (1e9 * currentSegment->chan->data->length *
              currentSegment->chan->deltaT);

    /* look for injection into segment */
    k = 0;
    thisInjection=injections;
    while (thisInjection)
    {
      INT8 ta = XLALGPSToINT8NS( &thisInjection->geocent_end_time );

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


INT4
XLALFindChirpTagTemplateAndSegment (
        DataSegmentVector       *dataSegVec,
        InspiralTemplate        *tmpltHead,
        SnglInspiralTable       **events,
        CHAR                    *ifo,
        REAL4                   tdFast,
        UINT4                   *analyseThisTmplt
        )

{
    UINT4                s, t;  /* s over segments and t over templates */
    SnglInspiralTable    *thisEvent = NULL;
    InspiralTemplate     *thisTmplt = NULL;
    DataSegment          *currentSegment = NULL;
    UINT4                flag = 0;
    INT8                 tc, chanStartTime, chanEndTime=0;

#ifndef LAL_NDEBUG
    /* Sanity checks on input arguments for debugging */
    if (!dataSegVec || !tmpltHead || !events || !(*events) ||
          !ifo || !analyseThisTmplt )
       XLAL_ERROR( XLAL_EFAULT );

    if ( tdFast < 0.0 || tdFast > 1.0 )
       XLAL_ERROR( XLAL_EINVAL );
#endif

    /* Do a IFO cut on the coinc list */
    (*events) =  XLALIfoCutSingleInspiral( events, ifo);

    if ( XLALClearErrno() ) {
       XLAL_ERROR( XLAL_EFUNC );
    }

    /* make sure the sngl inspirals are time ordered */
    *events = XLALSortSnglInspiral(*events, LALCompareSnglInspiralByTime);

    /* TODO: Make sure triggers are unique in masses or time */

    /* Loop over the coinc events */
    for (thisEvent=*events; thisEvent; thisEvent=thisEvent->next)
    {
        REAL8    g11, g12, g22;

        tc = XLALGPSToINT8NS( &thisEvent->end_time );
        flag = 0;

        /* Loop over segments */
        for (s = 0; s < dataSegVec->length; ++s )
        {

           /* point to current segment */
           currentSegment = dataSegVec->data + s;

           /* compute the start and end of segment */
           chanStartTime = XLALGPSToINT8NS( &currentSegment->chan->epoch );
           chanEndTime = chanStartTime +
                (INT8) (1e9 * currentSegment->chan->data->length *
                currentSegment->chan->deltaT);

           if ( tc > chanStartTime && tc <= chanEndTime )
           {
              flag += (1 << s);
           }


        } /* loop over segs */

        /* Check we haven't gone beyond the time of interest */
        if ( !flag )
        {
           if ( tc > chanEndTime )
              break;
           else
              continue;
        }


        /* Projected metric g_ij for this event */
        g11 = thisEvent->Gamma[3] -
           thisEvent->Gamma[1] * thisEvent->Gamma[1]/thisEvent->Gamma[0];

        g12 = thisEvent->Gamma[4] -
           thisEvent->Gamma[1]*thisEvent->Gamma[2]/thisEvent->Gamma[0];

        g22 = thisEvent->Gamma[5] -
           thisEvent->Gamma[2]*thisEvent->Gamma[2]/thisEvent->Gamma[0];

        /* Loop over templates */
        for ( thisTmplt = tmpltHead, t = 0; thisTmplt;
                  thisTmplt = thisTmplt->next, t++ )
        {

            REAL8    dt0, dt3;
            REAL8    match;

            dt0   = thisTmplt->t0 - thisEvent->tau0;
            dt3   = thisTmplt->t3 - thisEvent->tau3;
            match = g11*dt0*dt0 + 2.0*g12*dt0*dt3 + g22*dt3*dt3;
            match = 1.0 - match;


            /* If the match between the template and the event is high enough,
             * we mark the template to analyse the segments. This is achieved
             * using bitwise or. */
            if ( match >= tdFast )
            {
                analyseThisTmplt[t] = analyseThisTmplt[t] | flag;
            }

        } /* loop over templts */


    } /* loop over events */

    return XLAL_SUCCESS;

}




INT4
XLALFindChirpSetFollowUpSegment (
    DataSegmentVector          *dataSegVec,
    SnglInspiralTable          **events
    )

{
  DataSegment       *currentSegment;
  SnglInspiralTable *thisEvent;
  INT8               chanStartTime;
  INT8               chanEndTime;
  UINT4              i;
  UINT4              k;

  /* set all segments not to be analyzed by default */
  for ( i = 0; i < dataSegVec->length; ++i )
  {
    /* point to current segment */
    currentSegment = dataSegVec->data + i;
    currentSegment->analyzeSegment = 0;
  }

  /* make sure the sngl inspirals are time ordered */
  *events = XLALSortSnglInspiral(*events, LALCompareSnglInspiralByTime);

  /* loop over segments checking for injections into each */
  for ( i = 0; i < dataSegVec->length; ++i )
  {
    /* point to current segment */
    currentSegment = dataSegVec->data + i;

    /* compute the start and end of segment */
    chanStartTime = XLALGPSToINT8NS( &currentSegment->chan->epoch );
    chanEndTime = chanStartTime +
      (INT8) (1e9 * currentSegment->chan->data->length *
              currentSegment->chan->deltaT);

    /* look for event into segment */
    k = 0;
    thisEvent = *events;
    while (thisEvent)
    {
      INT8 ta = XLALGPSToINT8NS( &thisEvent->end_time );

      if ( ta > chanStartTime && ta <= chanEndTime )
      {
        currentSegment->analyzeSegment += 1;
      }

      if (ta > chanEndTime)
        break;

      thisEvent=thisEvent->next;
      k = k + 1;
    }
  }

  return 0;
}


void
LALFindChirpSetAnalyseTemplate (
    LALStatus                  *status,
    UINT4                      *analyseThisTmplt,
    REAL4                      mmFast,
    REAL8                      deltaF,
    INT4                       sampleRate,
    FindChirpDataParams        *fcDataParams,
    int                        UNUSED numTmplts,
    InspiralTemplate           *tmpltHead,
    int                        numInjections,
    SimInspiralTable           *injections
    )

{
  InspiralTemplate      *tmpltCurrent = NULL;
  REAL8FrequencySeries  *mmFshf = NULL;
  InspiralTemplate      *mmFTemplate = NULL;
  InspiralMetric        mmFmetric;
  InspiralMomentsEtc    mmFmoments;
  SimInspiralTable      *mmFInjection=NULL;
  UINT4                 mmF_i, ki, kj;
  REAL4                 dt0, dt3, metricDist, match;
  CHAR                  myMsg[8192];
  UINT4                 approximant;
  int                   oldxlalErrno;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  ASSERT( analyseThisTmplt, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( fcDataParams->wtildeVec->data, status,
      FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* Go through the next steps only if mmFast is used */
  /* AND injections is non-NULL                       */

  if ( (mmFast >= 0.0) && (injections) )
  {
    /* If mmFast option is used, assume NONE of the templates need to be */
    /* analysed to begin with. Thus re-init analyseThisTmplt elements to */
    /* zero                                                              */
    for ( tmpltCurrent = tmpltHead, kj = 0; tmpltCurrent;
        tmpltCurrent = tmpltCurrent->next, kj++ )
    {
      analyseThisTmplt[kj] = 0;
    }

    /* Get the approximant. If the approximant is not BCV or BCVSpin, then */
    /* we try to tag the templates ... the BCV waveforms are not included  */
    /* yet in the scheme of things                                         */
    oldxlalErrno = xlalErrno;
    xlalErrno = 0;
    approximant = XLALGetApproximantFromString(injections->waveform);
    if ( (int) approximant == XLAL_FAILURE)
      ABORTXLAL(status);
    xlalErrno = oldxlalErrno;

    if (approximant != (UINT4)BCV &&
        approximant != (UINT4)BCVSpin &&
        approximant != (UINT4)BCVC )
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
        if (crealf(fcDataParams->wtildeVec->data[ki]))
        {
          mmFshf->data->data[ki] =
            1./crealf(fcDataParams->wtildeVec->data[ki]);
        }
        else
        {
          /* Note that we can safely set shf to be zero as this is    */
          /* correctly handled in the LALGetInspiralMoments function  */
          mmFshf->data->data[ki] = 0.0;
        }
      }
      /* Init the template */
      mmFTemplate = (InspiralTemplate *)
        LALCalloc(1, sizeof(InspiralTemplate));
      mmFTemplate->fLower     = fcDataParams->fLow;
      mmFTemplate->fCutoff    = (REAL4)(sampleRate/2) - mmFshf->deltaF;
      mmFTemplate->tSampling  = (REAL4)(sampleRate);
      mmFTemplate->massChoice = m1Andm2;
      mmFTemplate->ieta       = 1.L;
      oldxlalErrno = xlalErrno;
      xlalErrno = 0;
      mmFTemplate->approximant = XLALGetApproximantFromString(
          injections->waveform);
      if ( (int) mmFTemplate->approximant == XLAL_FAILURE)
        ABORTXLAL(status);

      mmFTemplate->order = XLALGetOrderFromString(injections->waveform);
      if ( (int) mmFTemplate->order == XLAL_FAILURE)
        ABORTXLAL(status);
      xlalErrno = oldxlalErrno;

      snprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
          "%d Injections, Order = %d, Approx = %d\n\n",
          numInjections, mmFTemplate->order, mmFTemplate->approximant);
      LALInfo (status, myMsg);

      mmFTemplate->mass1      = injections->mass1;
      mmFTemplate->mass2      = injections->mass2;
      LALInspiralParameterCalc( status->statusPtr, mmFTemplate );
      CHECKSTATUSPTR (status);

      snprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
          "%d Injections, t0 = %e, t3 = %e\n",
          numInjections, mmFTemplate->t0, mmFTemplate->t3);
      LALInfo (status, myMsg);

      /* Now we are ready to calculate the noise moments */
      LALGetInspiralMoments( status->statusPtr, &mmFmoments,
          mmFshf, mmFTemplate);
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
        /* For this injection we must                                   */
        /* (a) Calculate the metric at the point of injection           */
        /* (b) Loop over all the templates and tag theie level to 0     */
        /* or 1 with respect to this injection.                         */
        mmFTemplate->mass1      = mmFInjection->mass1;
        mmFTemplate->mass2      = mmFInjection->mass2;
        LALInspiralParameterCalc( status->statusPtr, mmFTemplate );
        CHECKSTATUSPTR (status);

        snprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
            "%d Injections, m1 = %e, m2 = %e eta = %e\n",
            numInjections, mmFTemplate->mass1, mmFTemplate->mass2,
            mmFTemplate->eta);
        LALInfo (status, myMsg);
        snprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
            "%d Injections, t0 = %e, t3 = %e\n",
            numInjections, mmFTemplate->t0, mmFTemplate->t3);
        LALInfo (status, myMsg);

        LALInspiralComputeMetric( status->statusPtr, &mmFmetric,
            mmFTemplate, &mmFmoments );
        CHECKSTATUSPTR (status);

        snprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
            "%d Injections, G00 = %e, G01 = %e, G11 = %e\n\n",
            numInjections, mmFmetric.G00, mmFmetric.G01, mmFmetric.G11);
        LALInfo (status, myMsg);


        /* Now that the metric has been calculated we loop over the      */
        /* templates to mark their level.                                */

        /* kj is the index on the templates. Note that kj always starts  */
        /* from zero even if startTemplate and stopTemplate are specified*/
        /* by user while running lalapps_inspiral program                */
        kj = 0;

        for ( tmpltCurrent = tmpltHead; tmpltCurrent;
            tmpltCurrent = tmpltCurrent->next)
        {
          dt0 = tmpltCurrent->t0 - mmFTemplate->t0;
          dt3 = tmpltCurrent->t3 - mmFTemplate->t3;

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

          snprintf (myMsg, sizeof(myMsg)/sizeof(*myMsg),
              "%-5d %d %e %e %e %e %e %e %e %e %e %e %e %e %e\n",
              kj-1,
              analyseThisTmplt[kj-1],
              mmFTemplate->t0,
              mmFTemplate->t3,
              tmpltCurrent->t0,
              tmpltCurrent->t3,
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

      } /* End of loop over injections */

    } /* End of if (approximant is in the OK list) */
    else
    {
        LALWarning (status, "Requested --fast option for unsupported family of waveforms");
        LALWarning (status, "Tagging all templates to be analysed");

        for ( tmpltCurrent = tmpltHead, kj = 0; tmpltCurrent;
                tmpltCurrent = tmpltCurrent->next, kj++)
        {
            analyseThisTmplt[kj] = pow(2.0, (double)(numInjections)) - 1 ;
        }
    }


  } /* End of if (mmFast >= 0.0 && injections) */

  else if ( (mmFast >= 0.0) && (!injections) )
  {
    /* If mmFast option is used,  BUT the injections pointer is NULL, then */
    /* NONE of the templates need to be analysed.                          */
    for ( tmpltCurrent = tmpltHead, kj = 0; tmpltCurrent;
        tmpltCurrent = tmpltCurrent->next, kj++ )
    {
      analyseThisTmplt[kj] = 0;
    }
  }

  else /* if the --fast option is not used */
  {
    for ( tmpltCurrent = tmpltHead, kj = 0; tmpltCurrent;
        tmpltCurrent = tmpltCurrent->next, kj++)
    {
      if ( numInjections )
      {
        analyseThisTmplt[kj] = pow(2.0, (double)(numInjections)) - 1 ;
      }
      else
      {
        analyseThisTmplt[kj] = 1;
      }
    }
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}



UINT4
XLALCmprSgmntTmpltFlags (
    UINT4 numInjections,
    UINT4 TmpltFlag,
    UINT4 SgmntFlag
    )

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



UINT4
XLALFindChirpBankSimInitialize (
    REAL4FrequencySeries       *spec,
    COMPLEX8FrequencySeries    *resp,
    REAL8                       fLow
    )

{
  UINT4 k, cut;
  REAL4 psdMin = 0;
  const REAL8 psdScaleFac = 1.0e-40;

  /* set low frequency cutoff index of the psd and    */
  /* the value the psd at cut                         */
  cut = fLow / spec->deltaF > 1 ? fLow / spec->deltaF : 1;

  psdMin = spec->data->data[cut] *
    ( ( crealf(resp->data->data[cut]) * crealf(resp->data->data[cut]) +
        cimagf(resp->data->data[cut]) * cimagf(resp->data->data[cut]) ) / psdScaleFac );

  /* calibrate the input power spectrum, scale to the */
  /* range of REAL4 and store the as S_v(f)           */
  for ( k = 0; k < cut; ++k )
  {
    spec->data->data[k] = psdMin;
  }
  for ( k = cut; k < spec->data->length; ++k )
  {
    REAL4 respRe = crealf(resp->data->data[k]);
    REAL4 respIm = cimagf(resp->data->data[k]);
    spec->data->data[k] = spec->data->data[k] *
      ( ( respRe * respRe + respIm * respIm ) / psdScaleFac );
  }

  /* set the response function to the sqrt of the inverse */
  /* of the psd scale factor since S_h = |R|^2 S_v        */
  for ( k = 0; k < resp->data->length; ++k )
  {
    resp->data->data[k] = crectf( sqrt( psdScaleFac ), 0 );
  }

  return cut;
}

SimInspiralTable *
XLALFindChirpBankSimInjectSignal (
    DataSegmentVector          *dataSegVec,
    COMPLEX8FrequencySeries    *resp,
    SimInspiralTable           *injParams,
    FindChirpBankSimParams     *simParams
    )
{
  LALStatus             status;
  SimInspiralTable     *bankInjection;
  CHAR                  tmpChName[LALNameLength];
  REAL4                 M, mu;
#if 0
  LALFrStream             *frStream = NULL;
  REAL4TimeSeries       frameData;
  INT8                  waveformStartTime = 0;
  UINT4                 waveformLengthCtr = 0;
  UINT4                 i;
#endif

  memset( &status, 0, sizeof(LALStatus) );

  if ( injParams && simParams )
  {
    XLALPrintError(
        "XLAL Error: specify either a sim_inspiral table or sim params\n" );
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }

  /* create memory for the bank simulation */
  bankInjection = (SimInspiralTable *)
    LALCalloc( 1, sizeof(SimInspiralTable) );

  if ( simParams && simParams->approx == FrameFile )
  {
#ifdef LAL_FRAME_ENABLED

#if 0
    /* add the waveform from a frame file */
    if ( ! simParams->frameName || ! simParams->frameChan )
    {
      XLALPrintError(
          "XLAL Error: frame name and channel name must be specified\n" );
      XLAL_ERROR_NULL( XLAL_EINVAL );
    }

    fprintf( stderr, "reading data from %s %s\n",
        simParams->frameName, simParams->frameChan );
    snprintf( bankInjection->waveform, LIGOMETA_WAVEFORM_MAX * sizeof(CHAR),
        "%s", simParams->frameName );
    snprintf( bankInjection->source, LIGOMETA_SOURCE_MAX * sizeof(CHAR),
        "%s", simParams->frameChan );

    frStream = XLALFrStreamOpen( NULL, simParams->frameName );
    XLALFrStreamSetMode( frStream, LAL_FR_STREAM_VERBOSE_MODE );

    memset( &frameData, 0, sizeof(REAL4TimeSeries) );
    snprintf( frameData.name,
        sizeof(frameData.name) / sizeof(*(frameData.name)), "%s",
        simParams->frameChan );

    frameData.data =
      XLALCreateREAL4Vector( dataSegVec->data->chan->data->length );

    XLALFrStreamGetREAL4TimeSeries( &frameData, frStream );
    if ( xlalErrno == XLAL_EIO )
    {
      /* ignore end of data: we need the channel */
      /* to be longer than the waveform          */
      xlalErrno = 0;
    }

    fprintf( stderr, "template frame epoch is %d.%d\n",
        frameData.epoch.gpsSeconds, frameData.epoch.gpsNanoSeconds );
    fprintf( stderr, "template frame sampling rate = %le\n",
        frameData.deltaT );
    fprintf( stderr, "expected template sampling rate = %le\n",
        dataSegVec->data->chan->deltaT );

    XLALFrStreamClose( frStream );

    /* center the waveform in the data segment */
    waveformStartTime = XLALGPSToINT8NS( &(dataSegVec->data->chan->epoch) );

    for ( i = 0; i < frameData.data->length; ++i )
    {
      if ( frameData.data->data[i] ) ++waveformLengthCtr;
    }

    waveformStartTime += (INT8) ( 1000000000.0 *
        ((REAL8) ( dataSegVec->data->chan->data->length - waveformLengthCtr ) /
         2.0) * dataSegVec->data->chan->deltaT );

    XLALINT8NSToGPS( &(frameData.epoch), waveformStartTime );

    LALSSInjectTimeSeries( &status, dataSegVec->data->chan, &frameData );
    if ( status.statusCode )
    {
      REPORTSTATUS( &status );
      XLAL_ERROR_NULL( XLAL_EFAILED );
    }

    XLALDestroyREAL4Vector( frameData.data );
#endif
    XLALPrintError( "XLAL Error: frame reading not implemented\n" );
    XLAL_ERROR_NULL( XLAL_EINVAL );
#else
    XLALPrintError( "XLAL Error: LAL not compiled with frame support\n" );
    XLAL_ERROR_NULL( XLAL_EINVAL );
#endif
  }
  else
  {
    /* we read from the injParams or generate our own params */

    if ( injParams )
    {
      /* use injParams so copy the parameters from the input table */
      memcpy( bankInjection, injParams, sizeof(SimInspiralTable) );
      bankInjection->next = NULL;
    }
    else
    {
      /* use simParams so generate our own injection parameters */
      bankInjection->f_lower = simParams->f_lower;
      /* set up the injection masses */
      if ( simParams->maxMass == simParams->minMass )
      {
        bankInjection->mass1 = simParams->maxMass;
        bankInjection->mass2 = simParams->maxMass;
      }
      else
      {
        /* generate random parameters for the injection */
        bankInjection->mass1 = XLALUniformDeviate( simParams->randParams );
        bankInjection->mass1 *= (simParams->maxMass - simParams->minMass);
        bankInjection->mass1 += simParams->minMass;

        bankInjection->mass2 = XLALUniformDeviate( simParams->randParams );
        bankInjection->mass2 *= (simParams->maxMass - simParams->minMass);
        bankInjection->mass2 += simParams->minMass;
      }

      /* set up derived mass quantities */
      M = bankInjection->mass1 + bankInjection->mass2;
      mu = bankInjection->mass1 * bankInjection->mass2 / M;
      bankInjection->eta =  mu / M;
      bankInjection->mchirp = pow( mu, 3.0/5.0) * pow( M, 2.0/5.0 );

      /* set the correct waveform approximant */
      if ( simParams->approx == TaylorT1 )
      {
        snprintf( bankInjection->waveform, LIGOMETA_WAVEFORM_MAX,
            "TaylorT1twoPN" );
      }
      else if ( simParams->approx == TaylorT2 )
      {
        snprintf( bankInjection->waveform, LIGOMETA_WAVEFORM_MAX,
            "TaylorT2twoPN" );
      }
      else if ( simParams->approx == TaylorT3 )
      {
        snprintf( bankInjection->waveform, LIGOMETA_WAVEFORM_MAX,
            "TaylorT3twoPN" );
      }
      else if ( simParams->approx == PadeT1 )
      {
        snprintf( bankInjection->waveform, LIGOMETA_WAVEFORM_MAX,
            "PadeT1twoPN" );
      }
      else if ( simParams->approx == EOB )
      {
        snprintf( bankInjection->waveform, LIGOMETA_WAVEFORM_MAX,
            "EOBtwoPN" );
      }
      else if ( simParams->approx == GeneratePPN )
      {
        snprintf( bankInjection->waveform, LIGOMETA_WAVEFORM_MAX,
            "GeneratePPNtwoPN" );
      }
      else
      {
        XLALPrintError(
            "error: unknown waveform for bank simulation injection\n" );
        XLAL_ERROR_NULL( XLAL_EINVAL );
      }

      /* set the injection distance to 1 Mpc */
      bankInjection->distance = 1.0;
    }

    /* inject the signals, preserving the channel name (Tev mangles it) */
    snprintf( tmpChName, LALNameLength * sizeof(CHAR), "%s",
        dataSegVec->data->chan->name );

    /* make sure the injection is hplus with no time delays */
    dataSegVec->data->chan->name[0] = 'P';
    LALFindChirpInjectSignals( &status,
        dataSegVec->data->chan, bankInjection, resp );
    if ( status.statusCode )
    {
      REPORTSTATUS( &status );
      XLAL_ERROR_NULL( XLAL_EFAILED );
    }

    /* restore the saved channel name */
    snprintf( dataSegVec->data->chan->name,
        LALNameLength * sizeof(CHAR), "%s", tmpChName );
  }

  /* return a pointer to the created sim_inspiral table */
  return bankInjection;
}


REAL4
XLALFindChirpBankSimSignalNorm(
    FindChirpDataParams         *fcDataParams,
    FindChirpSegmentVector      *fcSegVec,
    UINT4                        cut
    )
{
  /* compute the minimal match normalization */
  UINT4 k;
  REAL4 matchNorm = 0;
  REAL4 *tmpltPower = fcDataParams->tmpltPowerVec->data;
  COMPLEX8 *wtilde = fcDataParams->wtildeVec->data;
  COMPLEX8 *fcData = fcSegVec->data->data->data->data;

  /* compute sigmasq for the injected waveform and store in matchNorm */
  switch ( fcDataParams->approximant )
  {
    case FindChirpSP:
      /* up to Nyquist */
      for ( k = cut; k < fcDataParams->tmpltPowerVec->length; ++k )
      {
        if ( tmpltPower[k] ) matchNorm += ( crealf(fcData[k]) * crealf(fcData[k]) +
            cimagf(fcData[k]) * cimagf(fcData[k]) ) / tmpltPower[k];
      }
      break;

    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case PadeT1:
    case EOB:
    case GeneratePPN:
    case BCVSpin:
    case FindChirpPTF:
      /* integrated up to Nyquist*/
      for ( k = cut; k < fcDataParams->wtildeVec->length; ++k )
      {
        if ( crealf(wtilde[k]) ) matchNorm += ( crealf(fcData[k]) * crealf(fcData[k]) +
            cimagf(fcData[k]) * cimagf(fcData[k]) ) / crealf(wtilde[k]);
      }
      break;

    default:
      XLALPrintError( "Error: unknown approximant\n" );
      XLAL_ERROR_REAL4( XLAL_EINVAL );
  }

  matchNorm *= ( 4.0 * (REAL4) fcSegVec->data->deltaT ) /
    (REAL4) ( 2 * (fcSegVec->data->data->data->length - 1) );

  /* square root to get minimal match normalization \sqrt{(s|s)} */
  matchNorm = sqrt( matchNorm );
  return matchNorm;
}


SimInstParamsTable *
XLALFindChirpBankSimMaxMatch (
    SnglInspiralTable         **bestTmplt,
    REAL4                       matchNorm
    )
{
  SnglInspiralTable     *loudestEvent = NULL;

  /* allocate memory for the loudest event over the template bank */
  loudestEvent = (SnglInspiralTable *) LALCalloc( 1, sizeof(SnglInspiralTable) );

  /* find the loudest snr over the template bank */
  while ( *bestTmplt )
  {
    SnglInspiralTable *tmpEvent = *bestTmplt;

    if ( (*bestTmplt)->snr > loudestEvent->snr )
      memcpy( loudestEvent, *bestTmplt, sizeof(SnglInspiralTable) );

    *bestTmplt = (*bestTmplt)->next;
    LALFree( tmpEvent );
  }

  /* make sure we return a null terminated list with only the best tmplt */
  loudestEvent->next = NULL;
  *bestTmplt = loudestEvent;

  /* return the match for the loudest event */
  return XLALFindChirpBankSimComputeMatch( *bestTmplt, matchNorm );
}


SimInstParamsTable *
XLALFindChirpBankSimComputeMatch (
    SnglInspiralTable   *tmplt,
    REAL4                matchNorm
    )
{
  SimInstParamsTable    *maxMatch = NULL;

  /* create the match output table */
  maxMatch = (SimInstParamsTable *) LALCalloc( 1, sizeof(SimInstParamsTable) );
  snprintf( maxMatch->name, LIGOMETA_SIMINSTPARAMS_NAME_MAX, "match" );

  /* store the match in the sim_inst_params structure */
  maxMatch->value = tmplt->snr / matchNorm;

  return maxMatch;

}

static int FindTimeSeriesStartAndEnd (
        REAL4Vector *signalvec,
        UINT4 *start,
        UINT4 *end
        )
{
  UINT4 i; /* mid, n; indices */
  UINT4 flag, safe = 1;
  UINT4 length;

#ifndef LAL_NDEBUG
  if ( !signalvec )
    XLAL_ERROR( XLAL_EFAULT );

  if ( !signalvec->data )
    XLAL_ERROR( XLAL_EFAULT );
#endif

  length = signalvec->length;

  /* Search for start and end of signal */
  flag = 0;
  i = 0;
  while(flag == 0 && i < length )
  {
      if( signalvec->data[i] != 0.)
      {
          *start = i;
          flag = 1;
      }
      i++;
  }
  if ( flag == 0 )
  {
      return flag;
  }

  flag = 0;
  i = length - 1;
  while(flag == 0)
  {
      if( signalvec->data[i] != 0.)
      {
          *end = i;
          flag = 1;
      }
      i--;
  }

  /* Check we have more than 2 data points */
  if(((*end) - (*start)) <= 1)
  {
      XLALPrintWarning( "Data less than 3 points in this signal!\n" );
      safe = 0;
  }

  return safe;

}
