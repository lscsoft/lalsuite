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

NRCSID( INSPIRALSLAVEC, "$Id$" );

void
LALFindChirpSlave (
    LALStatus                  *status, 
    BOOLEAN                    *notFinished,
    FindChirpSegmentVector     *fcSegVec,
    FindChirpSlaveParams       *params 
    )
{
  InitExchParams                initExchParams;
  ExchParams                    exchInspiralTemplates;
  ExchParams                    exchInspiralEvents;
  ExchParams                    exchFinished;
  ExchParams                   *thisExch = NULL;

  InspiralTemplate             *bankHead     = NULL;
  InspiralTemplate             *bankCurrent  = NULL;

  InspiralEvent                *eventList = NULL;

  INT4                          myRank;
  INT4                          numTmplts;
  UINT4                         i;

  
  INITSTATUS( status, "FindChirpSlave", INSPIRALSLAVEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * set up mpi and exchange types
   *
   */


  MPI_Comm_rank( *(params->mpiComm), &myRank );

  /* init exchange params */
  initExchParams.mpiComm   = *(params->mpiComm);
  initExchParams.myProcNum = myRank;

  /* exchange inspiral templates */
  exchInspiralTemplates.exchObjectType = ExchInspiralTemplate;
  exchInspiralTemplates.send           = 0; /* master sends */
  exchInspiralTemplates.numObjects     = 1;
  exchInspiralTemplates.partnerProcNum = 0; /* master */

  /* exchange inspiral events */
  exchInspiralEvents.exchObjectType    = ExchInspiralEvent;
  exchInspiralEvents.send              = 1; /* I send */
  exchInspiralEvents.numObjects        = 1; /* I send */
  exchInspiralEvents.partnerProcNum    = 0; /* master */
  
  /* exchange finished message */
  exchFinished.exchObjectType         = ExchFinished;
  exchFinished.send                   = 0; /* master sends */
  exchFinished.numObjects             = 0; /* irrelevant */
  exchFinished.partnerProcNum         = 0; /* master */


  /*
   *
   * get template parameter bank from master
   *
   */


  /* initialize the exchange of templates */
  LALInitializeExchange( status->statusPtr, &thisExch, 
      &exchInspiralTemplates, &initExchParams );
  CHECKSTATUSPTR( status );

  /* see how many templates we are going to be sent */
  LALExchangeUINT4( status->statusPtr, &numTmplts, thisExch );
  CHECKSTATUSPTR( status );

  /* if the master did not return any templates shut down */
  if ( !numTmplts )
  {
    /* close the template bank exchange */
    LALFinalizeExchange( status->statusPtr, &thisExch );
    CHECKSTATUSPTR( status );

    /* tell the master the slave is shutting down */
    LALInitializeExchange( status->statusPtr, &thisExch, 
        &exchFinished, &initExchParams );
    CHECKSTATUSPTR( status );

    LALFinalizeExchange( status->statusPtr, &thisExch );
    CHECKSTATUSPTR( status );

    /* and tell the wrapper the slave is shut down */
    *notFinished = 0;

    /* return control to ApplySearch() */
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }

  /* get the template bank */
  LALExchangeTemplateBank( status->statusPtr, &bankHead, thisExch );
  CHECKSTATUSPTR( status );

  /* finalize the template bank exchange */
  LALFinalizeExchange( status->statusPtr, &thisExch );
  CHECKSTATUSPTR( status );


  /*
   *
   * filter engine
   *
   */


  for ( bankCurrent = bankHead; bankCurrent; bankCurrent = bankCurrent->next)
  {

    /* generate a stationary phase template */
    LALFindChirpSPTemplate( status->statusPtr, params->filterInput->fcTmplt, 
        bankCurrent, params->tmpltParams );
    CHECKSTATUSPTR( status );

    params->filterInput->tmplt = bankCurrent;

    /* set the thresholds depending on the level of the template */
    params->filterParams->rhosqThresh = 
      params->rhosqThreshVec[bankCurrent->level];
    params->filterParams->chisqThresh = 
      params->chisqThreshVec[bankCurrent->level];


    /* 
     *
     * loop over data segments 
     *
     */


    for ( i = 0; i < fcSegVec->length; i++ )
    {
      eventList = NULL;

      /* should we filter this segment against this template? */
      if ( bankCurrent->segmentIdVec->data[fcSegVec->data[i].number] )
      {
        params->filterInput->segment = fcSegVec->data + i;

        /* filter the data segment against the template */
        LALFindChirpFilterSegment( status->statusPtr, &eventList, 
            params->filterInput, params->filterParams );
        CHECKSTATUSPTR( status );

        /* if the filter returned any events and we don't want the just the */
        /* loudest event, send them to the serach master                    */
        if ( eventList && ! (params->inspiralDebugFlag == 2) )
        {
          InspiralEvent *event;

          LALInitializeExchange( status->statusPtr, &thisExch, 
              &exchInspiralEvents, &initExchParams );
          CHECKSTATUSPTR( status );

          /* exchange all the events returned by the filter */
          for ( event = eventList; event; event = event->next )
          {
            LALExchangeInspiralEvent( status->statusPtr, event, thisExch );
            CHECKSTATUSPTR( status );
          }

          LALFinalizeExchange( status->statusPtr, &thisExch );
          CHECKSTATUSPTR( status );

          /* free the events */
          while ( eventList )
          {
            InspiralEvent *current;
            current = eventList;
            eventList = eventList->next;
            LALFree( current );
            current = NULL;
          }
        } 
        else if ( eventList && (params->inspiralDebugFlag == 2) )
        {
          InspiralEvent *thisEvent = eventList;
          InspiralEvent *lastEvent;

          while ( thisEvent )
          {
#if 0
            fprintf( stderr, "slave: searching for loudest event...\n" );
#endif

            /* save a copy of the loudest event */
            if ( thisEvent->snrsq > params->loudestEvent->snrsq )
            {
#if 0
              fprintf( stderr, "slave: found a new loudest event!\n" );
#endif
              memcpy( params->loudestEvent, thisEvent, sizeof(InspiralEvent) );
            }

            /* free the events */
            lastEvent = thisEvent;
            thisEvent = thisEvent->next;
            LALFree( lastEvent );
            lastEvent = NULL;
          }

          /* make sure no events are returned to the wrapperAPI */
          eventList = NULL;

        } /* end if... eventList */

      } /* end if... filter */

    } /* end loop over data segments */

  } /* end loop over linked list */


  if ( params->inspiralDebugFlag == 2 )
  {
    /* send the loudest event to that search master */
#if 0
    fprintf( stderr, "slave: sending loudest event to master\n" );
#endif

    LALInitializeExchange( status->statusPtr, &thisExch, 
        &exchInspiralEvents, &initExchParams );
    CHECKSTATUSPTR( status );

    LALExchangeInspiralEvent( status->statusPtr, params->loudestEvent, thisExch );
    CHECKSTATUSPTR( status );

    LALFinalizeExchange( status->statusPtr, &thisExch );
    CHECKSTATUSPTR( status );
  }
    

  LALFindChirpDestroyInspiralBank( status->statusPtr, &bankHead );
  CHECKSTATUSPTR( status );

  *notFinished = 1;

  DETATCHSTATUSPTR( status );
  RETURN( status );
}
