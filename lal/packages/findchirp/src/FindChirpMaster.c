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
  InspiralEvent                *thisEvent = NULL;
  InspiralTemplate             *fineBank  = NULL;
  InspiralTemplateNode         *thisTmplt = NULL;


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


  if ( *(params->numSlaves) )
  {
    ExchParams *thisExch = NULL;

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

        /* if there are any templates left send them to the slave */
        if ( params->tmpltCurrent )
        {
          UINT4                        tmpNumTmplts;
          InspiralTemplate            *tmpBankHead = NULL;
          InspiralTemplateNode        *tmpTmpltCurrent;

          /* count the number of templates to send */
          for ( tmpTmpltCurrent = params->tmpltCurrent, tmpNumTmplts = 0; 
              tmpTmpltCurrent && tmpNumTmplts < numTmpltExch;
              tmpTmpltCurrent = tmpTmpltCurrent->next )
          { 
            ++tmpNumTmplts;
          }

          /* tell the slave how many templates we have for it. actually */
          /* it doesn't care about the number, as long as it is > 0     */
          LALExchangeUINT4( status->statusPtr, &tmpNumTmplts, thisExch );
          CHECKSTATUSPTR( status );

          /* allocate some memory for a temporary inspiral bank */
          tmpBankHead = (InspiralTemplate *) 
            LALCalloc( tmpNumTmplts, sizeof(InspiralTemplate) );

          /* copy the templates from the bank to the temporary storage */
          /* and set the linked list pointers in the array             */
          for ( i = 0; i < tmpNumTmplts; ++i )
          {
            memcpy( tmpBankHead + i, params->tmpltCurrent->tmpltPtr,
                sizeof(InspiralTemplate) );
            (tmpBankHead + i)->next = NULL;
            (tmpBankHead + i)->fine = NULL;
            if ( i ) (tmpBankHead + i - 1)->next = (tmpBankHead + i);
            params->tmpltCurrent = params->tmpltCurrent->next;
          }

          /* exchange the temporary template bank... */
          LALExchangeTemplateBank( status->statusPtr, &tmpBankHead, thisExch );
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


        /*
         *
         * this is the master part of the heirarchical search
         *
         */


        /* fprintf( stderr, "event handler started\n" ); */

        /* this is a dog: i should think of a better way to do this */

        /* loop through event list */
        for ( thisEvent = *eventList; thisEvent; 
            thisEvent = thisEvent->next )
        { 
          /* look for a template that matches the template id of the event */
          /* fprintf( stderr, "searching event at %p\n", thisEvent ); */
          for ( thisTmplt = params->tmpltHead; thisTmplt; 
              thisTmplt = thisTmplt->next )
          {
            /* for each tmplt that we find, have we inserted that template? */
            /* fprintf( stderr, "searching tmplt at %p\n", thisTmplt ); */
            if ( thisEvent->tmplt.number == thisTmplt->tmpltPtr->number && 
                ! thisTmplt->inserted )
            {
              /* insert the fine bank into the list to filter */
              /* fprintf( stderr, "inserting fine bank %p \n", */
              /*    thisTmplt->tmpltPtr->fine ); */
              for ( fineBank = thisTmplt->tmpltPtr->fine; fineBank;
                  fineBank = fineBank->next )
              {
                /* fprintf( stderr, "fine tmplt %p\n", fineBank ); */
                LALFindChirpCreateTmpltNode( status->statusPtr, 
                    fineBank, &(params->tmpltCurrent) );
                CHECKSTATUSPTR( status );

                params->numTmplts++;
              }

              thisTmplt->inserted = 1;
            }
          }
        }

        /* fprintf( stderr, "event handler done\n" ); */

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

    /* finish the exchange */
    LALFinalizeExchange (status->statusPtr, &thisExch);
    CHECKSTATUSPTR(status);


    /*
     *
     * calaculate the progress information for wrapperAPI
     *
     */


    /* if there are any remaining templates, calculate the fraction */
    /* remaining by the progress through the template bank          */
    if ( params->tmpltCurrent )
    {
      /* this is arse */
      *(params->fracRemaining) = 1.0;
    }
    else /* we are all the way through the bank */
    {
      *(params->fracRemaining) = 0.0;
    }

    /* we're not finished, as we need to wait for the slaves to report finished */
    *(params->notFinished) = 1;

  }
  else /* there are no active slaves */
  {
    /* tell the wrapperAPI that we are finished */
    *(params->fracRemaining) = 0.0;
    *(params->notFinished)   = 0;
  }


  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
