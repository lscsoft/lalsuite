/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSlave.c
 *
 * Author: Brown, D. A., and Creighton, T. D.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/RealFFT.h>
#include <lal/BandPassTimeSeries.h>
#include <lal/TFTransform.h>
#include <lal/RingSearch.h>
#include <lal/FindChirpEngine.h>

NRCSID( FINDCHIRPSLAVEC, "$Id$" );

/*-----------------------------------------------------------------------*/
static void
CreateRandomPPNParamStruc (
    LALStatus          *status,
    PPNParamStruc      *params,
    REAL4              *m1,
    REAL4              *m2,
    REAL4               mMin,
    REAL4               mMax,
    REAL4               fmin,
    REAL8               deltaT,
    RandomParams       *rpar
    )
{
  REAL4         mDiff      = mMax - mMin;
  REAL4         cannonDist = 1.0e6;       /* cannonical distance in pc */

  INITSTATUS( status, "CreateRandomPPNParamStruc", FINDCHIRPSLAVEC );
  ATTATCHSTATUSPTR( status );

  /* fixed parameters. */
  params->position.latitude = params->position.longitude = 0.0;
  params->position.system = COORDINATESYSTEM_EQUATORIAL;
  params->psi = 0.0;
  params->lengthIn = 0;
  params->epoch.gpsSeconds = 0;
  params->epoch.gpsNanoSeconds = 0;
  params->deltaT = deltaT;

  /* set up the masses */
  TRY( LALUniformDeviate( status->statusPtr, m1, rpar ), status );
  TRY( LALUniformDeviate( status->statusPtr, m2, rpar ), status );
  *m1 = mMin + mDiff * *m1;
  *m2 = mMin + mDiff * *m2;
  params->mTot = *m1 + *m2;
  params->eta = *m1 * *m2 / ( params->mTot * params->mTot );

  /* other params */
  params->inc = 0.0;
  params->phi = 0.0;
  params->d = cannonDist * LAL_PC_SI;
  params->fStartIn = fmin;

  /* fStopIn is negative to allow bypass of the bad pn waveform stop */
  params->fStopIn = -1.0 / 
    (6.0 * sqrt(6.0) * LAL_PI * params->mTot * LAL_MTSUN_SI);

  /* ppn parameter */
  params->ppn = NULL;
  
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*-----------------------------------------------------------------------*/
static void 
GetTemplateBankMPI ( 
    LALStatus          *status,
    InspiralTemplate **tmpltBankHead,
    UINT4             *numTmplts,
    InitExchParams     initExchParams
    )
{
  ExchParams    exchInspiralTemplates;
  ExchParams   *thisExchPtr = NULL;

  INITSTATUS( status, "GetTemplateBankMPI", FINDCHIRPSLAVEC );
  ATTATCHSTATUSPTR( status );

  if ( *tmpltBankHead )
  {
    ABORT( status, FINDCHIRPENGINEH_ENNUL, FINDCHIRPENGINEH_MSGENNUL );
  }

  /* set up mpi exchange inspiral template type */
  exchInspiralTemplates.exchObjectType = ExchInspiralTemplate;
  exchInspiralTemplates.send           = 0; /* master sends */
  exchInspiralTemplates.numObjects     = 1;
  exchInspiralTemplates.partnerProcNum = 0; /* master */

  LALInitializeExchange( status->statusPtr, &thisExchPtr, 
      &exchInspiralTemplates, &initExchParams );
  CHECKSTATUSPTR( status );

  /* see how many templates we are going to be sent */
  LALExchangeUINT4( status->statusPtr, numTmplts, thisExchPtr );
  CHECKSTATUSPTR( status );

  if ( *numTmplts )
  {
    /* get the template bank */
    LALExchangeTemplateBank( status->statusPtr, 
        tmpltBankHead, thisExchPtr );
    CHECKSTATUSPTR( status );
  }

  LALFinalizeExchange( status->statusPtr, &thisExchPtr );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*-----------------------------------------------------------------------*/
static void 
OutputEventsMPI ( 
    LALStatus          *status,
    InspiralEvent      *eventList,
    InitExchParams     initExchParams
    )
{
  ExchParams     exchInspiralEvents;
  ExchParams    *thisExchPtr = NULL;

  INITSTATUS( status, "OutputEventsMPI", FINDCHIRPSLAVEC );
  ATTATCHSTATUSPTR( status );

  if ( ! eventList )
  {
    ABORT( status, FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  }

  /* set up mpi exchange inspiral event type */
  exchInspiralEvents.exchObjectType    = ExchInspiralEvent;
  exchInspiralEvents.send              = 1; /* I send */
  exchInspiralEvents.numObjects        = 1; /* I send */
  exchInspiralEvents.partnerProcNum    = 0; /* master */

  LALInitializeExchange( status->statusPtr, &thisExchPtr, 
      &exchInspiralEvents, &initExchParams );
  CHECKSTATUSPTR( status );

  LALExchangeInspiralEventList (status->statusPtr, &eventList, thisExchPtr);

  LALFinalizeExchange( status->statusPtr, &thisExchPtr );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*-----------------------------------------------------------------------*/
static void 
ExchNumTmpltsFilteredMPI ( 
    LALStatus          *status,
    UINT4               numTmplts,
    UINT4               numEvents,
    InitExchParams      initExchParams
    )
{
  ExchParams    exchNumTmpltsFiltered;
  ExchParams   *thisExchPtr = NULL;

  INITSTATUS( status, "ExchNumTmpltsFiltered", FINDCHIRPSLAVEC );
  ATTATCHSTATUSPTR( status );

  exchNumTmpltsFiltered.exchObjectType = ExchNumTmpltsFiltered;
  exchNumTmpltsFiltered.send           = 1; /* I send */
  exchNumTmpltsFiltered.numObjects     = 1;
  exchNumTmpltsFiltered.partnerProcNum = 0; /* master */

  LALInitializeExchange( status->statusPtr, &thisExchPtr, 
      &exchNumTmpltsFiltered, &initExchParams );
  CHECKSTATUSPTR( status );

  LALExchangeUINT4( status->statusPtr, &numTmplts, thisExchPtr );
  CHECKSTATUSPTR( status );

  LALExchangeUINT4( status->statusPtr, &numEvents, thisExchPtr );
  CHECKSTATUSPTR( status );

  LALFinalizeExchange( status->statusPtr, &thisExchPtr );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*-----------------------------------------------------------------------*/
static REAL4 
MedianSpec(
    LALStatus  *status,
    REAL4      *p, 
    UINT4       j, 
    UINT4       flength, 
    UINT4       numSegs
    )
{
  /* p points to array of power spectra data over time slices */
  /* j is desired frequency offset into power spectra array   */
  /* flength is size of frequency series obtained from DFT    */
  /* numSegs is the number of time slices to be evaluated     */
  /* status points to LAL status struct passed into main      */
  /* returns the median value, over time slice at given freq. */

  INT4  outer  = 0;       /* local loop counter */
  INT4  middle = 0;       /* local loop counter */
  INT4  inner  = 0;       /* local loop counter */
  REAL4 returnVal = 0.0;  /* holder for return value */

  ASSERT( p, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( flength, status, 
      FINDCHIRPENGINEH_EZERO, FINDCHIRPENGINEH_MSGEZERO );
  ASSERT( j < flength, status, 
      FINDCHIRPENGINEH_EAOVR, FINDCHIRPENGINEH_MSGEAOVR );
  ASSERT( flength, status, 
      FINDCHIRPENGINEH_EZERO, FINDCHIRPENGINEH_MSGEZERO );

  /* allocate memory array for insert sort, test for success */
  REAL4 *s = (REAL4 *) LALMalloc( numSegs * sizeof(REAL4) );
  if ( ! s )
  {
    ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
  }

  /* zero out the sort array */
  memset( s, 0, numSegs * sizeof(REAL4) );

  /* scan time slices for a given frequency */
  for ( outer = 0; outer < numSegs; ++outer )
  {
    /* insert power value into sort array */
    REAL4 tmp = p[outer * flength + j]; /* obtain value to insert */
    for ( middle = 0; middle < numSegs; ++middle )
    {
      if ( tmp > s[middle] )
      {
        /* insert taking place of s[middle] */
        for ( inner = numSegs - 1; inner > middle; --inner )
        {
          s[inner] = s [inner - 1];  /* move old values */
        }
        s[middle] = tmp;   /* insert new value */
        break;  /* terminate loop */
      }
    }
  }  /* done inserting into sort array */

  /* check for odd or even number of segments */
  if ( numSegs % 2 )
  {
    /* if odd number of segments, return median */
    returnVal = s[numSegs / 2];
  }
  else
  {
    /* if even number of segments, return average of two medians */
    returnVal = 0.5 * (s[numSegs/2] + s[(numSegs/2) - 1]);
  }

  /* free memory used for sort array */
  LALFree( s );

  return returnVal;
}
/*-----------------------------------------------------------------------*/
void
LALFindChirpSlave (
    LALStatus                  *status, 
    InspiralEvent             **outputEventHandle,
    DataSegmentVector          *dataSegVec,
    FindChirpSlaveParams       *params 
    )
{
  UINT4                         i, j, k;
  UINT4                         simCount = 0;
  UINT4                         numEvents = 0;
  INT4                         *filterSegment    = NULL;
  InspiralTemplate             *tmpltBankHead    = NULL;
  InspiralTemplate             *currentTmplt     = NULL;
  InspiralEvent                *eventList        = NULL;
  InspiralEvent               **eventListHandle  = NULL;
  UINT4                         numberOfTemplates;
  InitExchParams                initExchParams;
  MPI_Comm                     *mpiComm;


  INITSTATUS( status, "FindChirpSlave", FINDCHIRPSLAVEC );
  ATTATCHSTATUSPTR( status );

  ASSERT( outputEventHandle, status, 
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( !*outputEventHandle, status, 
      FINDCHIRPENGINEH_ENNUL, FINDCHIRPENGINEH_MSGENNUL );

  ASSERT( dataSegVec, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( dataSegVec->data, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( dataSegVec->data->chan, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( dataSegVec->data->spec, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( dataSegVec->data->resp, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  
  ASSERT( params, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params->rhosqThreshVec, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params->chisqThreshVec, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params->fcSegVec, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params->fcSegVec->data, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params->dataParams, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params->tmpltParams, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params->filterParams, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params->filterInput, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  ASSERT( params->notFinished, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );
  

  /*
   * 
   * set up the exchange params
   *
   */


  ASSERT( params->mpiComm, status,
      FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );

  mpiComm = (MPI_Comm *) params->mpiComm;
  MPI_Comm_rank( initExchParams.mpiComm = *mpiComm, 
      &(initExchParams.myProcNum) );


  /*
   *
   * get template parameter bank
   *
   */


  /* get template bank from master using MPI */
  GetTemplateBankMPI( status->statusPtr, &tmpltBankHead, 
      &numberOfTemplates, initExchParams );
  CHECKSTATUSPTR( status );

  /* if no template bank is returned, tell the master that the slave */
  /* is shutting down                                                */
  if ( ! tmpltBankHead )
  {
    ExchParams       *thisExchPtr = NULL;
    ExchParams        exchFinished;

    /* exchange finished message */
    exchFinished.exchObjectType         = ExchFinished;
    exchFinished.send                   = 0; /* master sends */
    exchFinished.numObjects             = 0; /* irrelevant */
    exchFinished.partnerProcNum         = 0; /* master */

    /* tell the master the slave is shutting down */
    LALInitializeExchange( status->statusPtr, &thisExchPtr, 
        &exchFinished, &initExchParams );
    CHECKSTATUSPTR( status );

    LALFinalizeExchange( status->statusPtr, &thisExchPtr );
    CHECKSTATUSPTR( status );
  }

  /* if we have no template bank, we are finished and should exit */
  /* setting notFinished to zero (i.e. finished)                  */
  if ( ! tmpltBankHead )
  {
    *(params->notFinished) = 0;
    goto exit;
  }


  /*
   *
   * main loop: execute until we run out of data (simulated or otherwise)
   *
   */


  /* set the eventlist handle to contain the address of the eventList   */
  /* passed in by the calling function. if no events are ever found     */
  /* this remains NULL. If this is a bankSim, we free the events and    */
  /* reset the value to be NULL                                         */
  eventListHandle = outputEventHandle;

  /* if this is a simulation, initialize the simulation counter         */
  if ( params->simParams )
  {
    simCount = params->simParams->simCount;
  }

  while ( 1 )
  {

   /*
    *
    * simulated data generation
    *
    */


    /* if we are doing a simulation, regenerate the  data internally    */
    /* according to the rules below                                     */
    if ( params->simParams )
    {
      DataSegment    *currentDataSeg = dataSegVec->data;

      REAL8   deltaT = currentDataSeg->chan->deltaT;
      REAL8   deltaF = currentDataSeg->resp->deltaF;

      UINT4   tdLength = currentDataSeg->chan->data->length;
      UINT4   fdLength = currentDataSeg->resp->data->length;

      /* strain PSD should be scaled by this factor so it fits in a REAL4 */
      REAL8   psdFactor = 9.0e-46;

      REAL8   minPsd;
      REAL4   response;
      REAL4   transfer;
      REAL4   spectrum;

      if ( params->simParams->simType == fcBankMinimalMatch )
      {
        InspiralEvent  *loudestEvent = NULL;

        /* check that we have been passed a strain PSD by the datacondAPI */
        if ( ! params->haveSpec )
        {
          ABORT( status, FINDCHIRPENGINEH_ESPEC, FINDCHIRPENGINEH_MSGESPEC );
        }

        /* create storeage for loudest event and (s|s) for each segment */
        loudestEvent = params->simParams->loudestEvent = (InspiralEvent *)
          LALCalloc( dataSegVec->length, sizeof(InspiralEvent) );
        params->simParams->signalNorm = (REAL4 *)
          LALCalloc( dataSegVec->length, sizeof(REAL4) );
        if ( ! params->simParams->loudestEvent || 
            ! params->simParams->signalNorm )
        {
          ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
        }

        for ( i = 0; i < dataSegVec->length; ++i, ++currentDataSeg )
        {
          UINT4         waveStart;
          CoherentGW    waveform;
          PPNParamStruc ppnParams;
          REAL4         mass1, mass2;
          REAL4        *data = currentDataSeg->chan->data->data;
          COMPLEX8     *resp = currentDataSeg->resp->data->data;

          /* set the next pointers in the array so ExchInspiralEvent    */
          /* can treat it as a linked list later                        */
          if ( i ) (loudestEvent + i - 1)->next = loudestEvent + 1;

          /* data to zero, spectrum to LIGO I psd, response to constant */
          response = (REAL4) sqrt( psdFactor );
          transfer = 1.0 / response;
          LALLIGOIPsd( NULL, &minPsd, (REAL8) params->dataParams->fLow );
          memset( data, 0, tdLength * sizeof(REAL4) );
          for( k = 0; k < fdLength ; ++k )
          {
            resp[k].re = response;
            resp[k].im = 0.0;
          }

          /* generate a random pair of masses and the ppn params struct */
          memset( &ppnParams, 0, sizeof(PPNParamStruc) );
          CreateRandomPPNParamStruc( status->statusPtr, 
              &ppnParams, &mass1, &mass2, 
              params->simParams->mMin, params->simParams->mMax, 
              params->simParams->fLowInject, deltaT,
              params->simParams->randomParams );
          CHECKSTATUSPTR( status );

          /* set the parameters in the loudest event structure          */
          loudestEvent[i].id = loudestEvent[i].segmentNumber = 
            loudestEvent[i].tmplt.number  = i;
          loudestEvent[i].tmplt.mass1     = (REAL8) mass1;
          loudestEvent[i].tmplt.mass2     = (REAL8) mass2;
          loudestEvent[i].tmplt.eta       = (REAL8) ppnParams.eta;
          loudestEvent[i].tmplt.totalMass = (REAL8) ppnParams.mTot;

          /* generate the inspiral waveform                             */
          memset( &waveform, 0, sizeof(CoherentGW) );
          LALGeneratePPNInspiral( status->statusPtr, &waveform, &ppnParams );
          CHECKSTATUSPTR( status );
#ifdef INSPIRAL_SO_DEBUG
          fprintf( stdout, "%s\n", ppnParams.termDescription );
          fprintf( stdout, "fStartIn = %e\n", ppnParams.fStartIn );
          fprintf( stdout, "fStopIn  = %e\n", ppnParams.fStopIn );
          fprintf( stdout, "fStart   = %e\n", ppnParams.fStart );
          fprintf( stdout, "fStop    = %e\n", ppnParams.fStop );
          fprintf( stdout, "length   = %d\n", ppnParams.length );
          fprintf( stdout, "duration = %e\n", 
              (REAL4) ppnParams.length * deltaT );
          fflush( stdout );
#endif

          /* place the waveform in the middle of the data segment       */
          if ( ppnParams.length > tdLength )
          {
            ABORT( status, FINDCHIRPENGINEH_EWAVL, FINDCHIRPENGINEH_MSGEWAVL );
          }
          waveStart = (tdLength - ppnParams.length) / 2;

          /* now generate v(t): first we generate h_plus(t) with code   */
          /* stolen from Tev, then we multiply by 1/R (it's ok to do    */
          /* this in the time domain as R is frequency independent;     */
          /* this saves a call to LALSimulateCoherentGW                 */
          {
            REAL8 t, x;
            REAL8 dx = 1.0;
            REAL8 xMax = waveform.a->data->length - 1;
            REAL8 *phiData = waveform.phi->data->data;
            REAL4 *aData = waveform.a->data->data;
            for ( x = 0.0, t = 0.0 ; x < xMax; x += dx, t += ppnParams.deltaT ) 
            {
              UINT4 jj = floor( x );
              REAL8 frac = x - jj;
              REAL8 p = frac*phiData[jj+1] + ( 1.0 - frac )*phiData[jj];
              REAL8 ap = frac*aData[2*jj+2] + ( 1.0 - frac )*aData[2*jj];
              REAL8 hp = ap * cos( p );
              data[jj + waveStart] = transfer * (REAL4) hp;
            }
          }

          /* free waveform memory generated by LALGeneratePPNInspiral() */
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

        /* condition the data segments: this will compute S_h( |f_k| )  */
        LALFindChirpSPData( status->statusPtr, params->fcSegVec, dataSegVec, 
            params->dataParams );
        CHECKSTATUSPTR( status );
        params->dataConditioned = 1;

        /* compute the normalisation of the injected signal (s|s)       */
        for ( i = 0; i < dataSegVec->length; ++i )
        {
          REAL4     sigNorm = 0;
          REAL4    *tmpltPower = params->dataParams->tmpltPowerVec->data;
          COMPLEX8 *fcData = (params->fcSegVec->data + i)->data->data->data;

          for ( k = 0; k < fdLength; ++k )
          {
            if ( tmpltPower[k] )
              sigNorm += ( fcData[k].re * fcData[k].re +
                  fcData[k].im * fcData[k].im ) / tmpltPower[k];
          }
          sigNorm *= ( 4.0 * (REAL4) deltaT ) / (REAL4) tdLength;

          params->simParams->signalNorm[i] = sigNorm;
        }
      } /* end if ( fcBankMinimalMatch ) */
      else if ( params->simParams->simType == fcGaussianNoise ||
          params->simParams->simType == fcGaussianNoiseInject )
      {
        /* replace the data with gaussian noise */
        REAL4Vector  *chanVec = params->simParams->chan->data;

        REAL4        *chan = chanVec->data;
        REAL4        *spec = dataSegVec->data->spec->data->data;
        COMPLEX8     *resp = dataSegVec->data->resp->data->data;

        REAL4         gaussianVarsq = params->simParams->gaussianVarsq;
        REAL4         gaussianVar = sqrt( gaussianVarsq );

        REAL4         respFactor = 1.0 / 
          ( gaussianVar * (REAL4) sqrt( 2.0 * deltaT ) );

        /* fill data channel with gaussian noise */
        LALNormalDeviates( status->statusPtr, 
            chanVec, params->simParams->randomParams );
        CHECKSTATUSPTR( status );

        for ( j = 0; j < chanVec->length; ++j )
        {
          chan[j] *= gaussianVar;
        }

        spectrum = 2.0 * gaussianVarsq * deltaT;

        LALLIGOIPsd( NULL, &minPsd, (REAL8) params->dataParams->fLow );

        for( k = 0; k < fdLength; ++k )
        {
          REAL8 psd;
          REAL8 freq = (REAL8) k * deltaF;
          if ( freq < params->dataParams->fLow )
          {
            resp[k].re = 1.0 / (REAL4) sqrt( psdFactor * minPsd );
            resp[k].re *= respFactor;
          }
          else
          {
            LALLIGOIPsd( NULL, &psd, freq );
            resp[k].re = 1.0 / (REAL4) sqrt( psdFactor * psd );
            resp[k].re *= respFactor;
          }
          resp[k].im = 0.0;
          spec[k] = spectrum;
        }

        /* we have created a psd. don't recreate it later */
        params->haveSpec = 1;
      }

      /* if ConditionData() has created any events to inject, inject them */
      if ( params->simParams->injectEvent )
      {
        /* inject signals into the input data stream */
        LALFindChirpInjectSignals( status->statusPtr, 
            params->simParams->chan, params->simParams->injectEvent, 
            dataSegVec->data->resp );
        CHECKSTATUSPTR( status );
      }

      /* unless this is a min match sim we need to (re)condition the data */
      if ( params->simParams->simType != fcBankMinimalMatch )
      {
        params->dataConditioned = 0;
      }

      /* unless this is a minimal match or a gaussian noise simulation, */
      /* (re)bandpass the data and (re)compule the psd                  */
      if ( params->simParams->simType != fcBankMinimalMatch &&
          params->simParams->simType != fcGaussianNoise &&
          params->simParams->simType != fcGaussianNoiseInject )
      {
        params->haveSpec = 0;
      }
    }


   /*
    *
    * band pass the data to get rid of low frequency crap
    *
    */

    
    if ( ! params->bandPassed && ( params->specType == fcSpecMeanBandPassed || 
          params->specType == fcSpecMedianBandPassed ) )
    {
      REAL4                     fsafety = 0;
      PassBandParamStruc        highpassParam;
      REAL4TimeSeries          *rawChannel;

      if ( params->simParams && params->simParams->injectEvent )
      {
        rawChannel = params->simParams->chan;
      }
      else
      {
        rawChannel = (REAL4TimeSeries *) 
          LALCalloc( 1, sizeof(REAL4TimeSeries) );
        memcpy( rawChannel, dataSegVec->data->chan, sizeof(REAL4TimeSeries) );
        rawChannel->data = &(params->dataChannel);
      }

      /* Set up for a highpass filter */
      highpassParam.nMax = 4;
      fsafety = params->dataParams->fLow - 10.0;
      highpassParam.f1 = fsafety > 150.0 ? 150.0 : fsafety;
      highpassParam.f2 = -1.0;
      highpassParam.a1 = 0.1;
      highpassParam.a2 = -1.0;
      
      LALButterworthREAL4TimeSeries( status->statusPtr, 
          rawChannel, &highpassParam );
      CHECKSTATUSPTR (status);

      if ( ! (params->simParams && params->simParams->injectEvent) )
      {
        LALFree( rawChannel );
      }

      params->bandPassed = 1;
    }
    else
    {
      params->bandPassed = 1;
    }


   /*
    *
    * compute the PSD of the input data
    *
    */


    if ( ! params->haveSpec )
    {
      REAL4TimeSeries          *chanPtr = dataSegVec->data->chan;
      REAL4FrequencySeries     *specPtr = dataSegVec->data->spec;
      REAL4                     deltaT = chanPtr->deltaT;
      UINT4                     tdLength = chanPtr->data->length;
      UINT4                     fdLength = tdLength / 2 + 1;
      RAT4                      negRootTwo = { -1, 1 };
      LALUnit                   unit;
      LALUnitPair               pair;

      /* set the spectrum frequency series parameters */
      specPtr->epoch.gpsSeconds = chanPtr->epoch.gpsSeconds;
      specPtr->epoch.gpsNanoSeconds = chanPtr->epoch.gpsNanoSeconds;
      specPtr->f0 = 0;
      specPtr->deltaF = 1.0 / ((REAL8) chanPtr->data->length * deltaT);
      specPtr->sampleUnits = lalADCCountUnit;
      pair.unitOne = &lalADCCountUnit;
      pair.unitTwo = &lalHertzUnit;
      TRY( LALUnitRaise( status->statusPtr, &unit, 
            pair.unitTwo, &negRootTwo ), status );
      pair.unitTwo = &unit;
      TRY( LALUnitMultiply( status->statusPtr, &(specPtr->sampleUnits),
            &pair ), status );

      if ( params->specType == fcSpecMean || 
          params->specType == fcSpecMeanBandPassed )
      {
        /* compute a mean power spectrum estimate from the data */
        COMPLEX8FrequencySeries   segmentFFT;
        LALWindowParams           winParams;
        RealDFTParams            *dftParams = NULL;

        memset( &segmentFFT, 0, sizeof(COMPLEX8FrequencySeries) );
        memset( &winParams, 0, sizeof(LALWindowParams) );

        winParams.type   = Hann;
        winParams.length = tdLength;

        LALCCreateVector( status->statusPtr, &(segmentFFT.data), fdLength );
        CHECKSTATUSPTR( status );
        
        memset( specPtr->data->data, 0, fdLength * sizeof(REAL4) );

        LALCreateRealDFTParams( status->statusPtr , &dftParams, &winParams, 1 );
        CHECKSTATUSPTR( status );

        for ( i = 0; i < dataSegVec->length; ++i )
        {
          /* compute the fft of each data segment */
          LALComputeFrequencySeries( status->statusPtr, &segmentFFT, 
              dataSegVec->data[i].chan, dftParams );
          CHECKSTATUSPTR (status);

          for ( k = 0; k < fdLength; ++k )
          {
            REAL4 fftRe = segmentFFT.data->data[k].re;
            REAL4 fftIm = segmentFFT.data->data[k].im;
            specPtr->data->data[k] += fftRe * fftRe + fftIm * fftIm;
          }
        }

        /* renormalize the psd to conventions document standard */
        {
          REAL4 psdnorm = ( 2.0 * deltaT ) / (REAL4) dataSegVec->length;
          for ( k = 0; k < fdLength; ++k )
          {
            specPtr->data->data[k] *= psdnorm;
          }
        }

        LALDestroyRealDFTParams( status->statusPtr , &dftParams );
        CHECKSTATUSPTR( status );

        LALCDestroyVector( status->statusPtr, &(segmentFFT.data) );
        CHECKSTATUSPTR( status );
      }
      else if ( params->specType == fcSpecMedian ||
          params->specType == fcSpecMedianBandPassed )
      {
#if 0
        /* compute a median power spectrum estimate from the data */
        AvgSpecParams           avgParams;
        REAL4TimeSeries         dataChanF;
        RealFFTPlan            *fftPlan = NULL;

        LALCreateForwardRealFFTPlan( status->statusPtr, &fftPlan, tdLength, 0);
        CHECKSTATUSPTR( status );

        memcpy( &dataChanF, chanPtr, sizeof(REAL4TimeSeries) );
        dataChanF.data = &(params->dataChannel);

        avgParams.segsize = tdLength;
        avgParams.fwdplan = fftPlan;
        avgParams.wintype = Hann;

        LALMedianSpectrum( status->statusPtr, specPtr, &dataChanF, &avgParams );
        CHECKSTATUSPTR( status );

        LALDestroyRealFFTPlan (status->statusPtr, &fftPlan );
        CHECKSTATUSPTR (status);
#endif
        /* compute a median power spectrum estimate using power code */
        COMPLEX8FrequencySeries   segmentFFT;
        REAL4                    *psdSeg;
        LALWindowParams           winParams;
        RealDFTParams            *dftParams = NULL;

        memset( &segmentFFT, 0, sizeof(COMPLEX8FrequencySeries) );
        memset( &winParams, 0, sizeof(LALWindowParams) );

        LALCCreateVector( status->statusPtr, &(segmentFFT.data), fdLength );
        CHECKSTATUSPTR( status );
        
        memset( specPtr->data->data, 0, fdLength * sizeof(REAL4) );

        winParams.type   = Hann;
        winParams.length = tdLength;

        psdSeg = (REAL4 *) 
          LALCalloc( dataSegVec->length, fdLength * sizeof(REAL4) );

        LALCreateRealDFTParams( status->statusPtr , &dftParams, &winParams, 1 );
        CHECKSTATUSPTR( status );

        for ( i = 0; i < dataSegVec->length; ++i )
        {
          /* compute the fft of each data segment */
          LALComputeFrequencySeries( status->statusPtr, &segmentFFT, 
              dataSegVec->data[i].chan, dftParams );
          CHECKSTATUSPTR (status);

          for ( k = 0; k < fdLength; ++k )
          {
            REAL4 fftRe = segmentFFT.data->data[k].re;
            REAL4 fftIm = segmentFFT.data->data[k].im;
            psdSeg[i * fdLength + k] = 2.0 * (REAL4) deltaT * 
              (fftRe * fftRe + fftIm * fftIm);
          }
        }

        for ( k = 0; k < fdLength; ++k )
        {
          specPtr->data->data[k] = MedianSpec( status->statusPtr, 
                psdSeg, k, fdLength, dataSegVec->length );
          CHECKSTATUSPTR( status );

          specPtr->data->data[k] /= LAL_LN2;
        }
        
        LALFree( psdSeg );

        LALDestroyRealDFTParams( status->statusPtr , &dftParams );
        CHECKSTATUSPTR( status );

        LALCDestroyVector( status->statusPtr, &(segmentFFT.data) );
        CHECKSTATUSPTR( status );
      }
      else
      {
        ABORT( status, FINDCHIRPENGINEH_ESTPE, FINDCHIRPENGINEH_MSGESTPE );
      }

      params->haveSpec = 1;
    }


   /*
    *
    * if the data has not been conditioned yet, condition it
    *
    */


    /* check that we generated a PSD by now as we can't proceed without one */
    if ( ! params->haveSpec )
    {
      ABORT( status, FINDCHIRPENGINEH_ESPEC, FINDCHIRPENGINEH_MSGESPEC );
    }

    /* eventually there will be a choice of ways to condition the data based */
    /* on the type of template used (stationary phase, time domain, etc.)    */
    if ( ! params->dataConditioned )
    {
      LALFindChirpSPData( status->statusPtr, params->fcSegVec, dataSegVec, 
          params->dataParams );
      CHECKSTATUSPTR( status );

      params->dataConditioned = 1;
    }    


    /*
     *
     * if we have a place to store the standard candle, create it
     *
     */


    if ( params->candlePtr )
    {
      FindChirpStandardCandle   *candle = params->candlePtr;

      REAL4 cannonDist = 1.0; /* Mpc */
      REAL4 m  = (REAL4) candle->tmplt.totalMass;
      REAL4 mu = (REAL4) candle->tmplt.mu;
      REAL4 distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1e6 * LAL_PC_SI);
      REAL4 candleTmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
        pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) *
        pow( LAL_MTSUN_SI / (REAL4) params->tmpltParams->deltaT, -1.0/6.0 );

      distNorm *= params->tmpltParams->dynRange;
      candleTmpltNorm *= candleTmpltNorm;
      candleTmpltNorm *= distNorm * distNorm;

      candle->sigmasq = 4.0 * ( (REAL4) params->filterParams->deltaT / 
          (REAL4) params->filterParams->qVec->length );
      candle->sigmasq *= candleTmpltNorm * params->fcSegVec->data->segNorm;

      candle->effDistance = sqrt( candle->sigmasq / candle->rhosq );
    }


    /* 
     *
     * loop over linked list of template nodes 
     *
     */


    for ( currentTmplt = tmpltBankHead; currentTmplt; 
        currentTmplt = currentTmplt->next )
    {
      /* set the template pointer the address of the current template */
      params->filterInput->tmplt = currentTmplt;

      /* pointer to the array of segments to filter against the */
      /* current template: this is what we base the decision    */
      /* rule on                                                */
      filterSegment = currentTmplt->segmentIdVec->data;

      /* set the filter thresholds depending on the level of the template */
      params->filterParams->rhosqThresh = 
        params->rhosqThreshVec[currentTmplt->level];
      params->filterParams->chisqThresh = 
        params->chisqThreshVec[currentTmplt->level];


      /*
       *
       * generate the template to filter against
       *
       */
      

      /* enventually this will be replaced so that multiple forms of */
      /* template generation are allowed (not just stationary phase) */
      LALFindChirpSPTemplate( status->statusPtr, params->filterInput->fcTmplt, 
          currentTmplt, params->tmpltParams );
      CHECKSTATUSPTR( status );


      /* 
       *
       * loop over data segments 
       *
       */


      for ( i = 0; i < params->fcSegVec->length; i++ )
      {


        /*
         *
         * filter data segment
         *
         */

        
        /* is this segment marked for filtering against this template? */
        if ( filterSegment[i] )
        {
          /* point the filter input to this data segment */
          params->filterInput->segment = params->fcSegVec->data + i;

          /* filter the data segment against the template */
          LALFindChirpFilterSegment( status->statusPtr, eventListHandle, 
              params->filterInput, params->filterParams );
          CHECKSTATUSPTR( status );

          /* dereference the eventListHandle to get a pointer the events */
          eventList = *eventListHandle;

          /* process list of returned events */
          if ( eventList )
          {
            if ( params->simParams && 
             params->simParams->simType == fcBankMinimalMatch )
            {
              /* if we just want the loudest event, save only the loudest */
              /* in the list for the data segment                         */
              InspiralEvent *loudestEvent = params->simParams->loudestEvent;

              while ( eventList )
              {
                InspiralEvent *thisEvent = eventList;

                if ( thisEvent->snrsq > loudestEvent[i].snrsq )
                {
                  memcpy( &(loudestEvent[i].time), &(thisEvent->time),
                      sizeof(LIGOTimeGPS) );
                  memcpy( loudestEvent[i].ifoName, thisEvent->ifoName,
                      2 * sizeof(CHAR) );
                  memcpy( loudestEvent[i].channel, thisEvent->channel,
                      LALNameLength * sizeof(CHAR) );
                  loudestEvent[i].snrsq   = thisEvent->snrsq;
                  loudestEvent[i].chisq   = thisEvent->chisq;
                  loudestEvent[i].sigma   = thisEvent->sigma;
                  loudestEvent[i].effDist = thisEvent->effDist;

                  /* XXX remove these two lines XXX */
                  loudestEvent[i].tmplt.t0 = thisEvent->tmplt.mass1;
                  loudestEvent[i].tmplt.t2 = thisEvent->tmplt.mass2;
                }

                eventList = eventList->next;
                LALFree( thisEvent );
                thisEvent = NULL;
              }

              /* reset the output event list handle */
              *outputEventHandle = NULL;

              ASSERT( ! *eventListHandle, status, 
                  FINDCHIRPENGINEH_ENNUL, FINDCHIRPENGINEH_MSGENNUL );
            }
            else
            {
              /* send a UINT4Vector of seg and tmplt numbers to the master */
            }
            
            if ( eventList )
            {
              ++numEvents;
              /* set the event list handle to the lastEvent->next pointer */
              while ( eventList->next ) 
              {
                ++numEvents;
                eventList = eventList->next;
              }

              eventListHandle = &(eventList->next);

              ASSERT( ! *eventListHandle, status, 
                  FINDCHIRPENGINEH_ENNUL, FINDCHIRPENGINEH_MSGENNUL );
            }

          } /* if ( eventList ) */

        } /* if ( filterSegment[i] ) */

      } /* end loop over data segments */

    } /* end loop over templates */


    if ( params->simParams && 
        params->simParams->simType == fcBankMinimalMatch )
    {
      InspiralEvent    *loudestEvent = params->simParams->loudestEvent;
      REAL4            *signalNorm   = params->simParams->signalNorm;
      for ( i = 0; i < dataSegVec->length; ++i )
      {
        /* compute match amd store in snrsq field */
        loudestEvent[i].snrsq /= signalNorm[i];
      }

      /* send loudest events to the master */
      OutputEventsMPI ( status->statusPtr, params->simParams->loudestEvent, 
          initExchParams );
      CHECKSTATUSPTR( status );

      /* free the loudest event array and the normalisation */
      LALFree( params->simParams->signalNorm );
      LALFree( params->simParams->loudestEvent );
    }

    /* break out of the main loop if we have no more data */
    if ( ! params->simParams )
    {
      break;
    }
    else if ( ! --simCount ) 
    {
      break;
    }


  } /* end main loop */


  /*
   *
   * free memory and exit with appropiate finished status
   *
   */


  /* free the template bank */
  LALFindChirpDestroyInspiralBank( status->statusPtr, &tmpltBankHead );
  CHECKSTATUSPTR( status );

  /* we are not finished unless the master returns no templates.      */
  *(params->notFinished) = 1;

  /* return the number of templates that we have just filtered to the */
  /* master so that it can caluclate progress information for the     */
  /* wrapperAPI                                                       */
  ExchNumTmpltsFilteredMPI( status->statusPtr, numberOfTemplates, 
      numEvents, initExchParams );
  CHECKSTATUSPTR( status );

  /* normal exit */
exit:
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


