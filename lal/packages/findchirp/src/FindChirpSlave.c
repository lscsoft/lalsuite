/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSlave.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/FindChirpEngine.h>

NRCSID( FINDCHIRPSLAVEC, "$Id$" );

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
  UINT4                         i;
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
   * if the data has not been conditioned yet, condition it
   *
   */


  /* eventually there will be a choice of ways to condition the data based */
  /* on the type of template used (stationary phase, time domain, etc.)    */
  if ( ! params->dataConditioned )
  {
    LALFindChirpSPData( status->statusPtr, params->fcSegVec, 
       dataSegVec, params->dataParams );
    CHECKSTATUSPTR( status );

    params->dataConditioned = 1;
  }    


  /*
   *
   * main loop: execute until we run out of data (simulated or otherwise)
   *
   */


  while ( 1 )
  {


    /*
     *
     * simulation data generation
     *
     */


    /* if we are doing a simulation, regenerate the */
    /* data internally according to the rules below */
    if ( params->simParams )
    {
      /* decide what to do with the input data */
      if ( params->simParams->simType == bankMinimalMatch )
      {
        /* zero the input data */

        /* create a vector to store the loudest event for each segment */
      }
      else if ( params->simParams->simType == gaussianNoise ||
          params->simParams->simType == gaussianNoiseInject )
      {
        /* fill the data with gaussian noise */
      }

      /* decide what to inject into the input data */
      if ( params->simParams->simType == bankMinimalMatch ||
          params->simParams->simType == gaussianNoiseInject ||
          params->simParams->simType == realDataInject )
      {
        /* inject a signal into the data */
      }
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

          /* filter the data segment against the template */
          LALFindChirpFilterSegment( status->statusPtr, &eventList, 
              params->filterInput, params->filterParams );
          CHECKSTATUSPTR( status );

#if 0
          fprintf( stdout, "returned from filter segments with events at %p\n",
              eventList );
          fflush( stdout );
          fprintf( stdout, "insertedTemplates = %d\n", insertedTemplates );
          fflush( stdout );
          fprintf( stdout, "currentTmplt at %p\n", currentTmplt );
          fflush( stdout );
          fprintf( stdout, "fine template at %p\n", currentTmplt->fine );
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
      /* compute (s|s) */

      /* send events to the master */

      /* free loudest event array  */
    }

    /* break out of the main loop if we have no more data */
    if ( ! params->simParams )
    {
      break;
    }
    else if ( params->simParams->simCount-- ) 
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
