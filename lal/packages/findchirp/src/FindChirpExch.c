/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpExch.c
 *
 * Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpExchCV">
Author: Allen, B., Brown, D. A. and Creighton, J. D. E.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpExch.c}}
\label{ss:FindChirpExch.c}

This module functions to exchange \texttt{findchirp} specific data types
using MPI.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpExchCP}
\idx{LALExchangeDataSegment()}
\idx{LALExchangeInspiralTemplate()}
\idx{LALExchangeInspiralEventList()}
\idx{LALExchangeTemplateBank()}

\subsubsection*{Description}


\subsubsection*{Algorithm}

None.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpExchCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/DataBuffer.h>
#include <lal/FindChirpExch.h>

NRCSID (FINDCHIRPEXCHC, "$Id$");

#pragma <lalVerbatim file="FindChirpExchCP">
void
LALExchangeDataSegment (
    LALStatus          *status,
    DataSegment        *segment,
    ExchParams         *exchParams
    )
#pragma </lalVerbatim>
{
  CHARVector box; /* a box to hold some bytes of data */

  INITSTATUS (status, "LALExchangeDataSegment", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  /* only do a minimal check to see if arguments are somewhat reasonable */
  ASSERT (exchParams, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);
  ASSERT (segment, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);
  ASSERT (segment->chan, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);
  ASSERT (segment->spec, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);
  ASSERT (segment->resp, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);

  if (exchParams->send)  /* I am sending */
  {
    INT4 dest = exchParams->partnerProcNum;

    /* stuff the data segment into a box */
    box.length = sizeof (DataSegment);
    box.data   = (CHAR *) segment;

    /* send box (this sends too much, but it is simple) */
    LALMPISendCHARVector (status->statusPtr, &box, dest, exchParams->mpiComm);
    CHECKSTATUSPTR (status);

    /* send pointer fields of data segment */

    /* data */
    LALMPISendREAL4TimeSeries (status->statusPtr, segment->chan, dest,
        exchParams->mpiComm);
    CHECKSTATUSPTR (status);

    /* spec */
    LALMPISendREAL4FrequencySeries (status->statusPtr, segment->spec, dest,
        exchParams->mpiComm);
    CHECKSTATUSPTR (status);

    /* resp */
    LALMPISendCOMPLEX8FrequencySeries (status->statusPtr, segment->resp, dest,
        exchParams->mpiComm);
    CHECKSTATUSPTR (status);
  }
  else /* I am receiving */
  {
    DataSegment tmpSegment;
    INT4        source = exchParams->partnerProcNum;

    box.length = sizeof (DataSegment);
    box.data   = (CHAR *) &tmpSegment;

    /* receive box */
    LALMPIRecvCHARVector (status->statusPtr, &box, source,
        exchParams->mpiComm);
    CHECKSTATUSPTR (status);

    /* receive remaining fields */

    /* data */
    LALMPIRecvREAL4TimeSeries (status->statusPtr, segment->chan, source,
        exchParams->mpiComm);
    CHECKSTATUSPTR (status);

    /* spec */
    LALMPIRecvREAL4FrequencySeries (status->statusPtr, segment->spec, source,
        exchParams->mpiComm);
    CHECKSTATUSPTR (status);

    /* resp */
    LALMPIRecvCOMPLEX8FrequencySeries (status->statusPtr, segment->resp, source,
        exchParams->mpiComm);
    CHECKSTATUSPTR (status);
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


#pragma <lalVerbatim file="FindChirpExchCP">
void
LALExchangeInspiralTemplate (
    LALStatus          *status,
    InspiralTemplate   *tmplt,
    ExchParams         *exchParams
    )
#pragma </lalVerbatim>
{
  CHARVector box; /* a box to hold some bytes of data */

  INITSTATUS (status, "LALExchangeInspiralTemplate", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  /* only do a minimal check to see if arguments are somewhat reasonable */
  ASSERT (exchParams, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);
  ASSERT (tmplt, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);

  /* stuff the template into a box */
  box.length = sizeof (InspiralTemplate);
  box.data   = (CHAR *) tmplt;

  if (exchParams->send)  /* I am sending */
  {
    INT4 dest = exchParams->partnerProcNum;

    /* send box */
    LALMPISendCHARVector (status->statusPtr, &box, dest, exchParams->mpiComm);
    CHECKSTATUSPTR (status);

    /* exchange the length of the segmentId vector */
    LALExchangeUINT4( status->statusPtr, &(tmplt->segmentIdVec->length),
        exchParams );
    CHECKSTATUSPTR (status);

    /* and then send the vector */
    LALMPISendINT4Vector( status->statusPtr, tmplt->segmentIdVec,
        dest, exchParams->mpiComm);
    CHECKSTATUSPTR (status);
  }
  else /* I am receiving */
  {
    INT4 source = exchParams->partnerProcNum;
    UINT4      length;

    /* receive box */
    LALMPIRecvCHARVector (status->statusPtr, &box, source, exchParams->mpiComm);
    CHECKSTATUSPTR (status);

    /* exchange the length of the segmentId vector */
    LALExchangeUINT4( status->statusPtr, &length, exchParams );
    CHECKSTATUSPTR (status);
    
    /* allocate the memory for it */
    tmplt->segmentIdVec = NULL;
    LALI4CreateVector( status->statusPtr, &(tmplt->segmentIdVec), length);
    CHECKSTATUSPTR( status );

    /* and then recv the vector */
    LALMPIRecvINT4Vector( status->statusPtr, tmplt->segmentIdVec,
        source, exchParams->mpiComm);
    CHECKSTATUSPTR (status);
    
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


#pragma <lalVerbatim file="FindChirpExchCP">
void
LALExchangeInspiralEvent (
    LALStatus          *status,
    InspiralEvent      *event,
    ExchParams         *exchParams
    )
#pragma </lalVerbatim>
{
  CHARVector box; /* a box to hold some bytes of data */

  INITSTATUS (status, "LALExchangeInspiralEvent", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  /* only do a minimal check to see if arguments are somewhat reasonable */
  ASSERT (exchParams, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);
  ASSERT (event, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);

  /* stuff the event into a box */
  box.length = sizeof (InspiralEvent);
  box.data   = (CHAR *) event;

  if (exchParams->send)  /* I am sending */
  {
    INT4 dest = exchParams->partnerProcNum;

    /* send box */
    LALMPISendCHARVector (status->statusPtr, &box, dest, exchParams->mpiComm);
    CHECKSTATUSPTR (status);
  }
  else /* I am receiving */
  {
    INT4 source = exchParams->partnerProcNum;

    /* receive box */
    LALMPIRecvCHARVector (status->statusPtr, &box, source, exchParams->mpiComm);
    CHECKSTATUSPTR (status);
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


#pragma <lalVerbatim file="FindChirpExchCP">
void
LALExchangeInspiralEventList (
    LALStatus          *status,
    InspiralEvent     **eventHead,
    ExchParams         *exchParams
    )
#pragma </lalVerbatim>
{
  InspiralEvent   *eventCurrent = NULL;
  InspiralEvent   *prev         = NULL;
  InspiralEvent   *events       = NULL;
  UINT4            numEvents    = 0;
  INT4             i, code;
  MPI_Status       stat;

  INITSTATUS (status, "LALExchangeInspiralEvent", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  if (exchParams->send)  /* I am sending */
  {

    /* check that we have a bank to send */
    ASSERT (*eventHead, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);

    /* count the events */
    for (eventCurrent = *eventHead;
        eventCurrent != NULL;
        eventCurrent = eventCurrent->next
        )
    {
      numEvents++;
    }

    /* exchange the number of events */
    LALExchangeUINT4( status->statusPtr, &numEvents, exchParams );
    CHECKSTATUSPTR (status);

    /* allocate memory for array of pointers to events */
    events = NULL;
    events = (InspiralEvent *) LALMalloc (numEvents * sizeof(InspiralEvent));

    /*  Make sure that the allocation was succesful */
    if ( !(events) ){
      ABORT (status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);
    }

    /* copy out events into array */
    i=0;
    for (eventCurrent = *eventHead;
        eventCurrent != NULL;
        eventCurrent = eventCurrent->next
        )
    {
      i++;
      memcpy( (events + i-1), eventCurrent, sizeof(InspiralEvent));
    }

    /* send the array */
    code = MPI_Send( events, numEvents * sizeof(InspiralEvent), MPI_BYTE, 
        exchParams->partnerProcNum, MPIZ, exchParams->mpiComm );
    if ( code != MPI_SUCCESS )
    {
      ABORT( status, FINDCHIRPEXCHH_EMPIE, FINDCHIRPEXCHH_MSGEMPIE);
    }

    /* free the memory */
    LALFree(events);
  }
  else /* I am receiving */
  {

    /* check that this is a new list */
    ASSERT (!*eventHead, status, FINDCHIRPEXCHH_ENNUL, FINDCHIRPEXCHH_MSGENNUL);

    /* receive the number of events */
    LALExchangeUINT4( status->statusPtr, &numEvents, exchParams );
    CHECKSTATUSPTR (status);

    /* allocate memory for array of pointers to events */
    events = NULL;
    events = (InspiralEvent *) LALMalloc (numEvents * sizeof(InspiralEvent));

    /*  Make sure that the allocation was succesful */
    if ( !(events) ){
      ABORT (status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);
    }

    code = MPI_Recv( events, numEvents * sizeof(InspiralEvent), MPI_BYTE, 
        exchParams->partnerProcNum, MPIZ, exchParams->mpiComm, &stat );
    if ( code != MPI_SUCCESS || stat.MPI_ERROR != MPI_SUCCESS )
    {
      ABORT( status, FINDCHIRPEXCHH_EMPIE, FINDCHIRPEXCHH_MSGEMPIE);
    }

    /* recieve the template bank */
    i=0;
    do
    {
      /* create memory for the template */
      eventCurrent = (InspiralEvent *) 
        LALCalloc ( 1, sizeof(InspiralEvent) );

      /* make a note of the first node in the list to return */
      if ( *eventHead == NULL ) *eventHead = eventCurrent;

      /* copy the event from the array */
      memcpy( eventCurrent, (events + i), sizeof(InspiralEvent));
 
      /* point the previous node to this node */
      if ( prev != NULL ) prev->next = eventCurrent;
      prev = eventCurrent;

      /* increment counter */
      ++i;
    } while ( eventCurrent->next != NULL && i< (INT4)numEvents);

    /* Terminate the linked list */
    eventCurrent->next = NULL;

    LALFree(events);
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


#pragma <lalVerbatim file="FindChirpExchCP">
void
LALExchangeTemplateBank (
    LALStatus          *status,
    InspiralTemplate  **tmpltHead,
    ExchParams          *exchParms
    )
#pragma </lalVerbatim>
{
  InspiralTemplate   *tmpltCurrent = NULL;
  InspiralTemplate   *tmpltFine    = NULL;
  InspiralTemplate   *prev         = NULL;

  INITSTATUS (status, "LALExchangeTemplateBank", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  if (exchParms->send)  /* I am sending */
  {

    /* check that we have a bank to send */
    ASSERT (*tmpltHead, status, FINDCHIRPEXCHH_ENULL, FINDCHIRPEXCHH_MSGENULL);

    /* exchange the template bank */
    for (tmpltCurrent = *tmpltHead;
        tmpltCurrent != NULL;
        tmpltCurrent = tmpltCurrent->next
        )
    {
      LALExchangeInspiralTemplate (status->statusPtr, tmpltCurrent, exchParms);
      CHECKSTATUSPTR (status);

      /* exchange the fine grids */
      if (tmpltCurrent->fine != NULL)
      {
        LALExchangeTemplateBank (status->statusPtr, 
            &(tmpltCurrent->fine), exchParms);
        CHECKSTATUSPTR (status);
      }
    }

    ASSERT (!tmpltCurrent, status, 
        FINDCHIRPEXCHH_ENNUL, FINDCHIRPEXCHH_MSGENNUL);

  }
  else /* I am receiving */
  {

    /* check that this is a new list */
    ASSERT (!*tmpltHead, status, FINDCHIRPEXCHH_ENNUL, FINDCHIRPEXCHH_MSGENNUL);

    /* recieve the template bank */
    do
    {
      /* create memory for the template */
      tmpltCurrent = (InspiralTemplate *) 
        LALCalloc ( 1, sizeof(InspiralTemplate) );

      /* make a note of the first node in the list to return */
      if ( *tmpltHead == NULL ) *tmpltHead = tmpltCurrent;

      /* recieve the template */
      LALExchangeInspiralTemplate (status->statusPtr, tmpltCurrent, exchParms);
      CHECKSTATUSPTR (status);

      /* get the fine grids */
      if (tmpltCurrent->fine != NULL)
      {
        tmpltFine = NULL;
        LALExchangeTemplateBank (status->statusPtr, &tmpltFine, exchParms);
        CHECKSTATUSPTR (status);

        /* set the fine pointer */
        tmpltCurrent->fine = tmpltFine;
      }

      /* point the previous node to this node */
      if ( prev != NULL ) prev->next = tmpltCurrent;
      prev = tmpltCurrent;

    } while ( tmpltCurrent->next != NULL );

  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}
