/*----------------------------------------------------------------------- 
 * 
 * File Name: EPExch.c
 *
 * Author: 
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <lal/AVFactories.h>
#include <lal/Comm.h>
#include <lal/ExcessPower.h>
#include <lal/BurstSearch.h>
#include "EPExch.h"

NRCSID (EPEXCHC, "$Id$");

void
LALExchangeEPEventList (
    LALStatus     *status,
    BurstEvent   **eventHead,
    INT4           numEvents,
    ExchParams    *exchParams
    )
{
  INT4      i, code;
  BurstEvent   *eventCurrent = NULL;
  BurstEvent   *prev         = NULL;
  BurstEvent   *events        = NULL;
  MPI_Status stat;

  INITSTATUS (status, "LALExchangeEPEventList", EPEXCHC);
  ATTATCHSTATUSPTR (status);

  if (exchParams->send)  /* I am sending */
  {

    /* check that we have events to send */
    ASSERT (*eventHead, status, EPEXCH_ENULL, EPEXCH_MSGENULL);

    /* allocate memory for array of pointers to events */
    events = NULL;
    events = (BurstEvent *) LALMalloc (numEvents * sizeof(BurstEvent));

    /*  Make sure that the allocation was succesful */
    if ( !(events) ){
      ABORT (status, EPEXCH_ENULL, EPEXCH_MSGENULL);
    }

    /* copy out events into array */
    i=0;
    eventCurrent = *eventHead;
    while (eventCurrent != NULL && i < numEvents)
    {
      i++;
      memcpy( (events + i-1), eventCurrent, sizeof(BurstEvent));
      eventCurrent = eventCurrent->nextEvent;
    }

    code = MPI_Send( events, numEvents * sizeof(BurstEvent), MPI_BYTE, 
        exchParams->partnerProcNum, MPIZ, exchParams->mpiComm );
    if ( code != MPI_SUCCESS )
    {
      ABORT( status, EPEXCH_EMPIE, EPEXCH_MSGEMPIE);
    }

    LALFree(events);

  }
  else /* I am receiving */
  {

    /* check that this is a new list */
    ASSERT (!*eventHead, status, EPEXCH_ENNUL, EPEXCH_MSGENNUL);

    /* allocate memory for array of pointers to events */
    events = NULL;
    events = (BurstEvent *) LALMalloc (numEvents * sizeof(BurstEvent));

    /*  Make sure that the allocation was succesful */
    if ( !(events) ){
      ABORT (status, EPEXCH_ENULL, EPEXCH_MSGENULL);
    }

    code = MPI_Recv( events, numEvents * sizeof(BurstEvent), MPI_BYTE, 
        exchParams->partnerProcNum, MPIZ, exchParams->mpiComm, &stat );
    if ( code != MPI_SUCCESS || stat.MPI_ERROR != MPI_SUCCESS )
    {
      ABORT( status, EPEXCH_EMPIE, EPEXCH_MSGEMPIE);
    }

    i=0;
    /* Parse the events into a linked list*/
    do
    {
      /* create memory for the first event */
      eventCurrent = (BurstEvent *) LALMalloc (sizeof(BurstEvent));
      ASSERT (eventCurrent, status, EPEXCH_ENULL, EPEXCH_MSGENULL);

      /* make a note of the first node in the list to return */
      if ( *eventHead == NULL ) *eventHead = eventCurrent;

      /* copy the event from the array */
      memcpy( eventCurrent, (events + i), sizeof(BurstEvent));
      
      /* point the previous node to this node */
      if ( prev != NULL ) prev->nextEvent = eventCurrent;
      prev = eventCurrent;

      /* increment pointer */
      ++i;
    } while ( eventCurrent->nextEvent != NULL && i<numEvents);

    /* Terminate the linked list */
    eventCurrent->nextEvent = NULL;

    LALFree(events);

  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


