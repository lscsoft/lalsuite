/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpMaster.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/FindChirpEngine.h>

NRCSID( INSPIRALMASTERC, "$Id$" );

void
LALFindChirpMaster (
    LALStatus                  *status, 
    InspiralEvent             **eventList,
    FindChirpMasterParams      *params 
               )
{
  InitExchParams                initExchParams;
  

  INT4                          myRank;
  UINT4                         i;
  UINT4                         numTmpltExch = params->numTmpltExch;
  InspiralTemplate             *currentTmplt      = NULL;
  InspiralTemplate             *fineTmplt         = NULL;
  InspiralTemplateNode         *thisTmpltNode     = NULL;
  InspiralTemplateNode         *insertTmpltNode   = NULL;


  INITSTATUS( status, "LALFindChirpMaster", INSPIRALMASTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * start up the inspiral master
   *
   */


  /* get my rank */
  MPI_Comm_rank( *(params->mpiComm), &myRank );

  /* set the communicators for the exchanges */
  initExchParams.mpiComm = *(params->mpiComm);
  initExchParams.myProcNum = myRank;


  /*
   *
   * if there are any slaves left, respond to their requests
   *
   */


  /* keep looping waiting for exchange requests from the slaves         */
  /* we break out of this loop only if we recieve events of finished    */
  /* messages from the slaves                                           */

  while ( *(params->numSlaves) )
  {
    ExchParams *thisExch;

    thisExch = NULL;

    /* break out of the loop if there are no slaves */
    if ( ! *(params->numSlaves) )
    {
      break;
    }
    
    /* wait for a request from a slave */
    LALInitializeExchange( status->statusPtr, &thisExch, 
        NULL, &initExchParams );
    CHECKSTATUSPTR (status);

    /* case the recieved exchange type */
    switch ( thisExch->exchObjectType )
    {


      /*
       *
       * exchange templates
       *
       */


      case ExchInspiralTemplate:

        if ( *(params->inspiralDebugFlagPtr) == fcBankMinimalMatch )
        {
          UINT4         tmpNumTmplts = 1;

#if 0
          fprintf( stdout, "master: slave %d requested bank simulation\n", 
              thisExch->partnerProcNum );
          fflush( stdout );
#endif

          /* check that the bank simulation vector exists       */
          ASSERT( params->bankSentVec, status, 
              FINDCHIRPENGINEH_ENULL, FINDCHIRPENGINEH_MSGENULL ); 

          /* for a bank simulation, we want to send the whole   */
          /* bank once only to every slave                      */
          if ( ! params->bankSentVec->data[thisExch->partnerProcNum] )
          {
            /* tell the slave how many templates we have for it */
            LALExchangeUINT4( status->statusPtr, &(params->numTmpltsTotal), 
                thisExch );
            CHECKSTATUSPTR( status );

#if 0
            fprintf( stdout, "master: sending whole bank to slave %d\n", 
                thisExch->partnerProcNum );
            fflush( stdout );
#endif

            /* exchange the whole template bank...              */
            LALExchangeTemplateBank( status->statusPtr, 
                &(params->tmpltBankHead), thisExch );
            CHECKSTATUSPTR( status );

            /* update the progress info                         */
            params->numTmpltsToFilter += params->numTmpltsTotal;
            
            /* ...and remember that we have sent it             */
            ++params->bankSentVec->data[thisExch->partnerProcNum];
          }
          else
          {
            /* tell the slave that there are no templates       */
#if 0
            fprintf( stdout, "master: no more templates for slave %d\n", 
                thisExch->partnerProcNum );
            fflush( stdout );
#endif
            numTmpltExch = 0;
            LALExchangeUINT4( status->statusPtr, &numTmpltExch, thisExch );
            CHECKSTATUSPTR( status );
          }
        }
        else
        {
          /* if there are any templates left send them to the slave */
          if ( params->currentTmpltNode )
          {
            UINT4                        tmpNumTmplts = 0;
            InspiralTemplate            *tmpBankHead = NULL;
            InspiralTemplateNode        *tempCurrentTmpltNode;

            /* count the number of templates to send to the slave */
            for ( tempCurrentTmpltNode = params->currentTmpltNode; 
                tempCurrentTmpltNode && tmpNumTmplts < numTmpltExch;
                tempCurrentTmpltNode = tempCurrentTmpltNode->next )
            { 
              ++tmpNumTmplts;
            }

            /* tell the slave how many templates we have for it   */
            LALExchangeUINT4( status->statusPtr, &tmpNumTmplts, thisExch );
            CHECKSTATUSPTR( status );

            /* allocate some memory for a temporary inspiral bank */
            tmpBankHead = (InspiralTemplate *) 
              LALCalloc( tmpNumTmplts, sizeof(InspiralTemplate) );
            if ( ! tmpBankHead )
            {
              ABORT( status, 
                  FINDCHIRPENGINEH_EALOC, FINDCHIRPENGINEH_MSGEALOC );
            }

            /* copy the templates from the bank to the temporary  */
            /* template storage                                   */
            for ( i = 0; i < tmpNumTmplts; ++i )
            {
              memcpy( tmpBankHead + i, params->currentTmpltNode->tmpltPtr,
                  sizeof(InspiralTemplate) );

              /* point the segmentIdVec at the address of the     */
              /* id vector in the real bank                       */
              (tmpBankHead + i)->segmentIdVec = 
                params->currentTmpltNode->tmpltPtr->segmentIdVec;

              /* set the next pointer so the array is also a linked list */
              (tmpBankHead + i)->next = NULL;
              (tmpBankHead + i)->fine = NULL;
              if ( i ) (tmpBankHead + i - 1)->next = (tmpBankHead + i);

              /* increment the current template node pointer in the list */
              params->currentTmpltNode = params->currentTmpltNode->next;
            }

            /* exchange the temporary template bank... */
            LALExchangeTemplateBank( status->statusPtr, 
                &tmpBankHead, thisExch );
            CHECKSTATUSPTR( status );

            /* ...and destroy it */
            LALFree( tmpBankHead );
          }
          else /* no templates */
          {
            /* tell the slave that there are no templates */
            numTmpltExch = 0;
            LALExchangeUINT4( status->statusPtr, &numTmpltExch, thisExch );
            CHECKSTATUSPTR( status );
          }
        }

        break;


      /*
       *
       * exchange events
       *
       */


      case ExchInspiralEvent:

        /* recieve a linked list of inspiral events from a slave */
        LALExchangeInspiralEventList( status->statusPtr, eventList, thisExch );
        CHECKSTATUSPTR( status );


#if 0
        /*
         *
         * this is the master part of the heirarchical search
         *
         */


        /* we are guaranteed that a the linked list of events will      */
        /* come from only _one_ template and _one_ data segment         */
        /* so we can just look at the first event in the list for       */
        /* the information that we need (template id and segment id)    */

        /* template node to insert additional nodes after */
        insertTmpltNode = params->currentTmpltNode;

        /* loop over the template bank looking for a template that      */
        /* with the same id the id of the template in the retuened      */
        /* event: this is ineffienct for large template banks and       */
        /* should be optimised                                          */
        for ( thisTmpltNode = params->tmpltNodeHead; thisTmpltNode; 
            thisTmpltNode = thisTmpltNode->next )
        {
          if ( (*eventList)->tmplt.number == thisTmpltNode->tmpltPtr->number )
          {
            /* insert the fine bank into the list to filter */
            for ( fineTmplt = thisTmpltNode->tmpltPtr->fine; fineTmplt;
                fineTmplt = fineTmplt->next )
            {
              INT4 *filterSegment = fineTmplt->segmentIdVec->data;

              /* create a template node for this fine template */
              LALFindChirpCreateTmpltNode( status->statusPtr, 
                  fineTmplt, &insertTmpltNode );
              CHECKSTATUSPTR( status );

              /* set the id of the segment to filter against this template */
              filterSegment[(*eventList)->segmentNumber] = 1;

              /* increase the count of the total number of templates */
              params->numTmpltsToFilter++;

              /* clear the event list since this is not the lowest template */
              while ( *eventList )
              {
                InspiralEvent *tmpEvent;
                tmpEvent = *eventList;
                *eventList = (*eventList)->next;
                LALFree (tmpEvent);
                tmpEvent = NULL;
              }

            }

            /* break out of the for loop over the template bank */
            break;
          }
        }
#endif


        break;


        /*
         *
         * slave is returning number of templates filtered
         *
         */

      case ExchNumTmpltsFiltered:


        /* add the number of templates returned by the slave to the     */
        /* total number that have been filtered and update the progress */
        {
          UINT4 numTmpltsReturned = 0;
          
          LALExchangeUINT4( status->statusPtr, &numTmpltsReturned, 
              thisExch );
          CHECKSTATUSPTR( status );

          params->numTmpltsFiltered += numTmpltsReturned;
        }

        break;


        /*
         *
         * slave is finished
         *
         */


      case ExchFinished:

        /* decrement the number of active slaves */
        --*(params->numSlaves);

        break;


        /*
         *
         * unrecognized exchange type
         *
         */


      default:
        ABORT( status, FINDCHIRPENGINEH_EUEXT, FINDCHIRPENGINEH_MSGEUEXT );

    }


    /*
     *
     * break out of loop if a slave returns events or filtered tmplts
     *
     */


    if ( thisExch->exchObjectType == ExchNumTmpltsFiltered ||
        thisExch->exchObjectType == ExchInspiralEvent )
    {
      /* finish the exchange and return control to the wrapperAPI */
      LALFinalizeExchange (status->statusPtr, &thisExch);
      CHECKSTATUSPTR(status);

      break;
    }
    else
    {
      /* just finish the exchange and loop again */
      LALFinalizeExchange (status->statusPtr, &thisExch);
      CHECKSTATUSPTR(status);
    }

  } /* end while( *(params->numSlave) ) */


  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
