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

#include <lal/FindChirpEngine.h>

NRCSID( FINDCHIRPSLAVEC, "$Id$" );

/*----------------------------------------------------------------------*/
static void
I8ToLIGOTimeGPS( LIGOTimeGPS *output, INT8 input )
{
  /* convert INT8 nanoseconds to LIGOTimeGPS */
  INT8 s = input / 1000000000LL;
  output->gpsSeconds = (INT4)( s );
  output->gpsNanoSeconds = (INT4)( input - 1000000000LL*s );
  return;
}

/*----------------------------------------------------------------------*/
static void
LIGOTimeGPSToI8( INT8 *output, LIGOTimeGPS input )
{
  /* convert LIGOTimeGPS to INT8 nanoseconds */
  *output = (INT8) input.gpsNanoSeconds 
    + 1000000000LL * (INT8) input.gpsSeconds;
  return;
}

/*-----------------------------------------------------------------------*/
static void
CreateRandomPPNParamStruc (
    LALStatus          *status,
    PPNParamStruc      *params,
    REAL4              *m1,
    REAL4              *m2,
    REAL4               mMin,
    REAL4               mMax,
    REAL8               deltaT,
    RandomParams       *rand
    )
{
  REAL4         mDiff      = mMax - mMin;
  REAL4         cannonDist = 1.0e6;       /* cannonical distance in pc */
  REAL4         fmin       = 25.0;
  REAL4         fmax       = 500.0;

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
  TRY( LALUniformDeviate( status->statusPtr, m1, rand ), status );
  TRY( LALUniformDeviate( status->statusPtr, m2, rand ), status );
  *m1 = mMin + mDiff * *m1;
  *m2 = mMin + mDiff * *m2;
  params->mTot = *m1 + *m2;
  params->eta = *m1 * *m2 / ( params->mTot * params->mTot );

  /* other params */
  params->inc = 0.0;
  params->phi = 0.0;
  params->d = cannonDist * LAL_PC_SI;
  params->fStartIn = fmin;
  params->fStopIn = fmax;

  /* ppn parameter */
  params->ppn = NULL;
  
  DETATCHSTATUSPTR( status );
  RETURN( status );
}

/*-----------------------------------------------------------------------*/
#ifdef LAL_MPI_ENABLED
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
    ABORT( status, FINDCHIRPENGINEH_ENNUL, FINDCHIRPENGINEH_MSGENNUL );

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
#endif /* LAL_MPI_ENABLED */


/*-----------------------------------------------------------------------*/
static void
GetTemplateBankFile ( 
    LALStatus                  *status,
    InspiralTemplate          **tmpltBankHead,
    FILE                       *tmpltBankFilePtr
    )
{
  INITSTATUS( status, "GetTemplateBankFile", FINDCHIRPSLAVEC );

  /* read template bank from a specified file */

  RETURN( status );
}

/*-----------------------------------------------------------------------*/
#ifdef LAL_MPI_ENABLED
static void 
OutputEventsMPI ( 
    LALStatus          *status,
    InspiralEvent      *eventList,
    InitExchParams     initExchParams
    )
{
  ExchParams     exchInspiralEvents;
  ExchParams    *thisExchPtr = NULL;
  InspiralEvent *event;

  INITSTATUS( status, "OutputEventsMPI", FINDCHIRPSLAVEC );
  ATTATCHSTATUSPTR( status );

  if ( ! eventList )
    ABORT( status, FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL );

  /* set up mpi exchange inspiral event type */
  exchInspiralEvents.exchObjectType    = ExchInspiralEvent;
  exchInspiralEvents.send              = 1; /* I send */
  exchInspiralEvents.numObjects        = 1; /* I send */
  exchInspiralEvents.partnerProcNum    = 0; /* master */

  LALInitializeExchange( status->statusPtr, &thisExchPtr, 
      &exchInspiralEvents, &initExchParams );
  CHECKSTATUSPTR( status );

  for ( event = eventList; event; event = event->next )
  {
    LALExchangeInspiralEvent( status->statusPtr, event, thisExchPtr );
    CHECKSTATUSPTR( status );
  }

  LALFinalizeExchange( status->statusPtr, &thisExchPtr );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
#endif /* LAL_MPI_ENABLED */

/*-----------------------------------------------------------------------*/
static void
OutputEventsFile ( 
    LALStatus                  *status,
    InspiralEvent              *eventList,
    FILE                       *eventFilePtr
    )
{
  INITSTATUS( status, "OutputEventsFile", FINDCHIRPSLAVEC );

  /* write inspiral events to a file */

  RETURN( status );
}

/*-----------------------------------------------------------------------*/
#ifdef LAL_MPI_ENABLED
static void 
ExchNumTmpltsFilteredMPI ( 
    LALStatus          *status,
    UINT4               numTmplts,
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

  LALFinalizeExchange( status->statusPtr, &thisExchPtr );
  CHECKSTATUSPTR( status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
#endif /* LAL_MPI_ENABLED */

/*-----------------------------------------------------------------------*/
void
LALFindChirpSlave (
    LALStatus                  *status, 
    BOOLEAN                    *notFinished,
    DataSegmentVector          *dataSegVec,
    FindChirpSlaveParams       *params 
    )
{
  UINT4                         i, j, k;
  UINT4                        *filterSegment    = NULL;
  InspiralTemplate             *tmpltBankHead    = NULL;
  InspiralTemplate             *currentTmplt     = NULL;
  InspiralTemplateNode         *tmpltNodeHead    = NULL;
  InspiralTemplateNode         *currentTmpltNode = NULL;
  InspiralEvent                *eventList        = NULL;
  InspiralEvent                *event            = NULL;
#ifdef LAL_MPI_ENABLED
  UINT4                         numberOfTemplates;
  InitExchParams                initExchParams;
  MPI_Comm                     *mpiComm;
#endif /* LAL_MPI_ENABLED */


  INITSTATUS( status, "FindChirpSlave", FINDCHIRPSLAVEC );
  ATTATCHSTATUSPTR( status );

  umask( 0 );

  /*
   * 
   * if MPI is defined and we are using it, set up the exchange params
   *
   */


#ifdef LAL_MPI_ENABLED
  if ( params->useMPI )
  {
    mpiComm = (MPI_Comm *) params->mpiComm;
    MPI_Comm_rank( initExchParams.mpiComm = *mpiComm, 
        &(initExchParams.myProcNum) );
  }
#endif /* HAVE_MPI */


  /*
   *
   * get template parameter bank
   *
   */


  /* If we have mpi enabled, there is a choice between getting the      */
  /* template bank from a file or getting it from the master via mpi.   */
  /* If we do not have mpi enabled, then the only choice is to get it   */
  /* from the file.                                                     */
#ifdef LAL_MPI_ENABLED
  if ( params->useMPI )
  {
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
  }
  else
#endif /* LAL_MPI_ENABLED */
  {
    /* get template bank from a specified file */
    GetTemplateBankFile( status->statusPtr, 
        &tmpltBankHead, params->tmpltBankFilePtr );
    CHECKSTATUSPTR( status );
  }

  /* if we have no template bank, we are finished and should exit */
  /* setting notFinished to zero (i.e. finished)                  */
  if ( ! tmpltBankHead )
  {
    *notFinished = 0;
    goto exit;
  }


 /*
  *
  * create a linked list of template nodes for the coarse templates
  *
  */


  LALFindChirpCreateTmpltNode( status->statusPtr, tmpltBankHead, 
      &tmpltNodeHead );
  CHECKSTATUSPTR( status );

  currentTmpltNode = tmpltNodeHead;

  for ( currentTmplt = tmpltBankHead->next; currentTmplt; 
      currentTmplt = currentTmplt->next )
  {
    LALFindChirpCreateTmpltNode( status->statusPtr, currentTmplt, 
        &currentTmpltNode );
    CHECKSTATUSPTR( status );
  }


  /*
   *
   * main loop: execute until we run out of data (simulated or otherwise)
   *
   */


  while ( 1 )
  {
#if 0
    fprintf( stdout, "---- top of while(1) -------------------\n" );
    if ( params->simParams ) 
      fprintf( stdout, "simCount = %u\n", params->simParams->simCount );
    fflush( stdout );
#endif


    /*
     *
     * simulated data generation
     *
     */


    /* if we are doing a simulation, regenerate the  data internally    */
    /* according to the rules below                                     */
    if ( params->simParams )
    {
      if ( params->simParams->simType == bankMinimalMatch )
      {
        DataSegment    *currentDataSeg = dataSegVec->data;
        InspiralEvent  *loudestEvent = NULL;

        REAL4   dynRange   = params->dataParams->dynRange;
        REAL4   dynRangeSq = dynRange * dynRange;

        REAL8   deltaT = currentDataSeg->real4Data->deltaT;
        REAL8   deltaF = currentDataSeg->spec->deltaF;

        UINT4   tdLength = currentDataSeg->real4Data->data->length;
        UINT4   fdLength = currentDataSeg->spec->data->length;

        REAL4 fourDeltaToverN  = ( 4.0 * (REAL4) deltaT ) / (REAL4) tdLength;

        /* create storeage for loudest event and (s|s) for each segment */
#if 1
        loudestEvent = params->simParams->loudestEvent = (InspiralEvent *)
          LALCalloc( dataSegVec->length, sizeof(InspiralEvent) );
        params->simParams->signalNorm = (REAL4 *)
          LALCalloc( dataSegVec->length, sizeof(REAL4) );
        if ( ! params->simParams->loudestEvent || 
            ! params->simParams->signalNorm )
        {
          ABORT( status, FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
        }
#endif

        for ( i = 0; i < dataSegVec->length; ++i )
        {
          REAL8         psdFactor = 9.0e-46;
          REAL8         minPsd;
          REAL4         response;
          REAL4         transfer;
          UINT4         waveStart;
          CoherentGW    waveform;
          PPNParamStruc ppnParams;
          REAL4         mass1, mass2;
          REAL4        *data = currentDataSeg->real4Data->data->data;
          REAL4        *spec = currentDataSeg->spec->data->data;
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
            REAL8 psd;
            REAL8 freq = (REAL8) k * deltaF;
            if ( freq < params->dataParams->fLow )
            {
              spec[k] = minPsd;
            }
            else
            {
              LALLIGOIPsd( NULL, &psd, freq );
              spec[k] = (REAL4) psd;
            }
            resp[k].re = response;
            resp[k].im = 0.0;
          }

          /* generate a random pair of masses and the ppn params struct */
          memset( &ppnParams, 0, sizeof(PPNParamStruc) );
          CreateRandomPPNParamStruc( status->statusPtr, 
              &ppnParams, &mass1, &mass2, 
              params->simParams->mMin, params->simParams->mMax, deltaT,
             params->simParams->randomParams );
          CHECKSTATUSPTR( status );

          /* XXX set params so they match Bruce's 1.4,1.4 chirp */
          mass1 = mass2 = 1.4;
          ppnParams.mTot = mass1 + mass2;
          ppnParams.eta = mass1*mass2/( ppnParams.mTot*ppnParams.mTot );

          /* set the parameters in the loudest event structure          */
          loudestEvent[i].id = loudestEvent[i].segmentNumber = 
            loudestEvent[i].tmplt.number  = i;
          loudestEvent[i].tmplt.mass1     = (REAL8) mass1;
          loudestEvent[i].tmplt.mass2     = (REAL8) mass2;
          loudestEvent[i].tmplt.eta       = (REAL8) ppnParams.eta;
          loudestEvent[i].tmplt.totalMass = (REAL8) ppnParams.mTot;
#if 0
          fprintf( stdout, "random chirp generated: m1 = %e, m2 = %e\n", 
              mass1, mass2 );
          fflush( stdout );
#endif

          /* generate the inspiral waveform                             */
          memset( &waveform, 0, sizeof(CoherentGW) );
          LALGeneratePPNInspiral( status->statusPtr, &waveform, &ppnParams );
#if 0
          fprintf( stdout, "%d: %s\n", ppnParams.termCode,
	       ppnParams.termDescription );
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
            REAL4 *fData = waveform.f->data->data;
            REAL4 *aData = waveform.a->data->data;
            for ( x = 0.0, t = 0.0 ; x < xMax; x += dx, t += ppnParams.deltaT ) 
            {
              UINT4 jj = floor( x );
              REAL8 frac = x - jj;
              REAL8 p = frac*phiData[jj+1] + ( 1.0 - frac )*phiData[jj];
              REAL8 f = frac*fData[jj+1] + ( 1.0 - frac )*fData[jj];
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
          LALFree( waveform.a );
          LALFree( waveform.f );
          LALFree( waveform.phi );

          currentDataSeg++;
        }

        /* condition the data segments: this will compute S_h( |f_k| )  */
        LALFindChirpSPData( status->statusPtr, params->fcSegVec, dataSegVec, 
            params->dataParams );
        CHECKSTATUSPTR( status );
        params->dataConditioned = 1;

#if 0
        {
          FILE  *tdfp, *fdfp;

          if ( ! (tdfp = fopen( "/home/duncan/incomming/td.dat", "w" ))
              || ! (fdfp = fopen( "/home/duncan/incomming/fd.dat", "w" )) )
          {
            ABORT( status, 999, "an egregious error has occoured" );
          }

          for ( j = 0; j < tdLength; ++j )
          {
            fprintf( tdfp, "%u\t%e\n", j, 
                dataSegVec->data->real4Data->data->data[j] );
          }
          fclose( tdfp );

          for ( k = 0; k < fdLength; ++k )
          {
            fprintf( fdfp, "%u\t%e\t%e\t%e\t%e\n", k, 
                dataSegVec->data->spec->data->data[k],
                dataSegVec->data->resp->data->data[k].re,
                dataSegVec->data->resp->data->data[k].im,
                params->dataParams->wtildeVec->data[k].re,
                params->dataParams->wtildeVec->data[k].re );
          }
          fclose( tdfp );
          fclose( fdfp );
        }
#endif

        /* compute the normalisation of the injected signal (s|s)       */
        for ( i = 0; i < dataSegVec->length; ++i )
        {
          REAL4     sigNorm = params->simParams->signalNorm[i];
          REAL4    *tmpltPower = params->dataParams->tmpltPowerVec->data;
          COMPLEX8 *fcData = (params->fcSegVec->data + i)->data->data->data;

          for ( k = 0; k < fdLength; ++k )
          {
            if ( tmpltPower[k] )
              sigNorm += (fcData[k].re * fcData[k].re) / tmpltPower[k];
          }
          sigNorm *= fourDeltaToverN;
#if 0
          fprintf( stdout, "sigNorm[%u] = %e\n", i, sigNorm );
          fflush( stdout );
#endif
        }
      }
      else if ( params->simParams->simType == gaussianNoise ||
          params->simParams->simType == gaussianNoiseInject )
      {
        /* fill input data stream with gaussian noise */
      }

      if ( params->simParams->simType == gaussianNoiseInject ||
          params->simParams->simType == realDataInject )
      {
        /* inject a signal into the input data stream */
      }
    }


   /*
    *
    * if the data has not been conditioned yet, condition it
    *
    */


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
     * loop over linked list of template nodes 
     *
     */


    for ( currentTmpltNode = tmpltNodeHead; currentTmpltNode; 
        currentTmpltNode = currentTmpltNode->next )
    {
      UINT4 insertedTemplates = 0;
      
      /* set the template pointer the address of the current template */
      /* in the list of template nodes                                */
      currentTmplt = params->filterInput->tmplt = currentTmpltNode->tmpltPtr;

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
#if 0
          fprintf( stdout, "filtering segment %d against template %d\n",
              params->filterInput->segment->number,
              currentTmplt->number );
          fflush( stdout );
#endif
          /* filter the data segment against the template */
          LALFindChirpFilterSegment( status->statusPtr, &eventList, 
              params->filterInput, params->filterParams );
          CHECKSTATUSPTR( status );
#if 0
          fprintf( stdout, "events found at %p\n", eventList );
          fflush( stdout );
#endif
          /* process list of returned events */
          if ( eventList )
          {
            if ( currentTmplt->fine && ! insertedTemplates )
            {
              /* we need to insert the fine list of templates and tag */
              /* the inserted template to that this data segment      */
              /* is filtered against them                             */
              InspiralTemplate       *fineTmplt;
              InspiralTemplateNode   *insertNode = currentTmpltNode;

              for ( fineTmplt = currentTmplt; fineTmplt;
                  fineTmplt = fineTmplt->next )
              {
                LALFindChirpCreateTmpltNode( status->statusPtr, 
                    fineTmplt, &insertNode );
                CHECKSTATUSPTR( status );

                fineTmplt->segmentIdVec->data[i] = 1;
              }

              insertedTemplates = 1;
            }
            else if ( currentTmplt->fine && insertedTemplates )
            {
              /* we do not need to insert fine templates, but we need     */
              /* to tag the inserted fine templates so that this data     */
              /* segment is filtered against them                         */
              InspiralTemplateNode     *fineTmpltNode;
              /* record the level of this template                        */
              UINT4 fineTmpltLevel = currentTmpltNode->next->tmpltPtr->level;

              /* start at the fist fine template that we have inserted    */
              /* and continue until we reach the next coarse template     */
              for ( fineTmpltNode = currentTmpltNode->next;
                  fineTmpltNode->tmpltPtr->level == fineTmpltLevel;
                  fineTmpltNode = fineTmpltNode->next )
              {
                fineTmpltNode->tmpltPtr->segmentIdVec->data[i] = 1;
                fineTmpltNode = fineTmpltNode->next;
              }
            }
            else if ( params->simParams && 
             params->simParams->simType == bankMinimalMatch )
            {
              /* if we just want the loudest event, save only the loudest */
              /* in the list for the data segment                         */
              InspiralEvent *loudestEvent = params->simParams->loudestEvent;
              InspiralEvent *thisEvent;
#if 0
              fprintf( stdout, "searching for loudest event\n" );
              fflush( stdout );
#endif
              for ( thisEvent = eventList; thisEvent; 
                  thisEvent = thisEvent->next )
              {
                if ( thisEvent->snrsq > loudestEvent[i].snrsq )
                {
#if 0
                  fprintf( stdout, "found new loudest event: %e\n",
                      thisEvent->snrsq );
                  fflush( stdout );
#endif
                  memcpy( &(loudestEvent[i].time), &(thisEvent->time),
                      sizeof(LIGOTimeGPS) );
                  loudestEvent[i].snrsq   = thisEvent->snrsq;
                  loudestEvent[i].chisq   = thisEvent->chisq;
                  loudestEvent[i].sigma   = thisEvent->sigma;
                  loudestEvent[i].effDist = thisEvent->effDist;
                }
              }
            }
            else
            {
#ifdef LAL_MPI_ENABLED
              if ( params->useMPI )
              {
                OutputEventsMPI ( status->statusPtr, eventList, 
                    initExchParams );
                CHECKSTATUSPTR( status );
              }
              else
#endif /* LAL_MPI_ENABLED */
              {
                OutputEventsFile ( status->statusPtr, eventList, 
                    params->eventFilePtr );
                CHECKSTATUSPTR( status );
              }                
            }

            /* free the events */
            while ( eventList )
            {
              InspiralEvent *current = eventList;
              eventList = eventList->next;
              LALFree( current );
              current = NULL;
            }

          } /* if ( eventList ) */

        } /* if ( filterSegment[i] ) */

      } /* end loop over data segments */


      /* 
       *
       * if the next template up a level, remove all the inserted templates 
       *
       */


      if ( currentTmpltNode->next )
      {
        InspiralTemplate       *nextTmplt = currentTmpltNode->next->tmpltPtr;
        UINT4                   currentLevel = currentTmplt->level;

        if ( nextTmplt->level < currentLevel )
        {
          while ( currentTmpltNode->tmpltPtr->level == currentLevel )
          {
            memset( currentTmpltNode->tmpltPtr->segmentIdVec->data, 0, 
                currentTmpltNode->tmpltPtr->segmentIdVec->length 
                * sizeof( UINT4 ) );
            LALFindChirpDestroyTmpltNode( status->statusPtr, 
                &currentTmpltNode );
            CHECKSTATUSPTR( status );
          }
        }
      }


    } /* end loop over templates */


    if ( params->simParams && 
        params->simParams->simType == bankMinimalMatch )
    {
      InspiralEvent    *loudestEvent = params->simParams->loudestEvent;
      REAL4            *signalNorm   = params->simParams->signalNorm;
      for ( i = 0; i < dataSegVec->length; ++i )
      {
        /* compute match amd store in snrsq field */
        loudestEvent[i].snrsq /= signalNorm[i];
#if 0
        fprintf( stdout, "snrsq = %e, snr = %e\n", loudestEvent[i].snrsq,
            sqrt( loudestEvent[i].snrsq ) );
#endif
      }

      /* send loudest events to the master */
#ifdef LAL_MPI_ENABLED
      if ( params->useMPI )
      {
        OutputEventsMPI ( status->statusPtr, params->simParams->loudestEvent, 
            initExchParams );
        CHECKSTATUSPTR( status );
      }
      else
#endif /* LAL_MPI_ENABLED */
      {
        OutputEventsFile ( status->statusPtr, params->simParams->loudestEvent, 
            params->eventFilePtr );
        CHECKSTATUSPTR( status );
      }                

      /* free the loudest event array and the normalisation */
      LALFree( params->simParams->signalNorm );
      LALFree( params->simParams->loudestEvent );
    }

    /* break out of the main loop if we have no more data */
    if ( ! params->simParams )
    {
      break;
    }
    else if ( ! --(params->simParams->simCount) ) 
    {
      break;
    }


  } /* end main loop */


  /*
   *
   * free memory and exit with appropiate finished status
   *
   */


  /* free linked list of template nodes */
  while ( tmpltNodeHead )
  {
    currentTmpltNode = tmpltNodeHead;
    tmpltNodeHead = tmpltNodeHead->next;
    LALFree( currentTmpltNode );
    currentTmpltNode = NULL;
  }

  /* free the template bank */
  LALFindChirpDestroyInspiralBank( status->statusPtr, &tmpltBankHead );
  CHECKSTATUSPTR( status );

  /* if we are using mpi then we are not finished unless the master     */
  /* returns no templates. if we are not using mpi, then we are done.   */
#ifdef LAL_MPI_ENABLED
  if ( params->useMPI )
  {
    *notFinished = 1;
    /* return the number of templates that we have just filtered to the */
    /* master so that it can caluclate progress information for the     */
    /* wrapperAPI                                                       */
    ExchNumTmpltsFilteredMPI( status->statusPtr, numberOfTemplates, 
        initExchParams );
    CHECKSTATUSPTR( status );
  }
  else
#endif /* LAL_MPI_ENABLED */
  {
    *notFinished = 0;
  }

  /* normal exit */
exit:
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
