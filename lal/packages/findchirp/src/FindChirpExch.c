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

#include <lal/LALStdlib.h>
#include <lal/DataBuffer.h>
#include <lal/FindChirpExch.h>

NRCSID (FINDCHIRPEXCHC, "$Id$");

void
LALInitializeExchange (
    LALStatus      *status,
    ExchParams **exchParamsOut,
    ExchParams  *exchParamsInp,
    InitExchParams *params
    )
{
  MPIMessage hello; /* initialization message */

  INITSTATUS (status, "LALInitializeExchange", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  ASSERT (exchParamsOut, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (!*exchParamsOut, status, FINDCHIRPEXCH_ENNUL, FINDCHIRPEXCH_MSGENNUL);

  /* now allocate memory for output exchange parameters */
  *exchParamsOut = (ExchParams *) LALMalloc (sizeof(ExchParams));
  ASSERT (*exchParamsOut, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  /* memset( *exchParamsOut, 0, sizeof(ExchParams) ); */

  if (exchParamsInp) /* I am initializing the exchange */
  {
    INT4 dest = exchParamsInp->partnerProcNum;

    /* initialize communications */

    hello.msg    = exchParamsInp->exchObjectType;
    hello.source = params->myProcNum;

    /*
     * just use send field as a means to communicate number of objects:
     * set to negative to indicate that initializer (I) want to receive
     * those objects
     *
     * add one to the number of objects in case it is zero
     */

    ASSERT (exchParamsInp->numObjects >= 0, status,
            FINDCHIRPEXCH_ENOBJ, FINDCHIRPEXCH_MSGENOBJ);
    if (exchParamsInp->send)
    {
      hello.send = exchParamsInp->numObjects + 1;
    }
    else
    {
      hello.send = -(exchParamsInp->numObjects + 1);
    }

    /* send off the communications */
    LALMPISendMsg (status->statusPtr, &hello, dest, params->mpiComm);
    CHECKSTATUSPTR (status);

    /* copy the input structure to the output structure */
    (*exchParamsOut)->exchObjectType = exchParamsInp->exchObjectType;
    (*exchParamsOut)->send           = exchParamsInp->send;
    (*exchParamsOut)->numObjects     = exchParamsInp->numObjects;
    (*exchParamsOut)->partnerProcNum = exchParamsInp->partnerProcNum;

    /* copy the communicator from the parameter structure */
    (*exchParamsOut)->myProcNum      = params->myProcNum;
    (*exchParamsOut)->mpiComm        = params->mpiComm;
  }
  else /* I am waiting for someone else to initialize the exchange */
  {
    /* wait for incoming message */
    LALMPIRecvMsg (status->statusPtr, &hello, params->mpiComm);
    CHECKSTATUSPTR (status);

    /* the message contains all the information needed */
    (*exchParamsOut)->exchObjectType = hello.msg;
    (*exchParamsOut)->send           = (hello.send < 0);
    (*exchParamsOut)->numObjects     = abs(hello.send) - 1;
    (*exchParamsOut)->partnerProcNum = hello.source;

    /* copy the communicator from the parameter structure */
    (*exchParamsOut)->myProcNum      = params->myProcNum;
    (*exchParamsOut)->mpiComm        = params->mpiComm;
  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

void
LALFinalizeExchange (
    LALStatus      *status,
    ExchParams **exchParams
    )
{
  INT2       magic = 0xA505; /* A SOS */
  INT2Vector goodbye;

  INITSTATUS (status, "LALFinalizeExchange", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  ASSERT (exchParams, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (*exchParams, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);

  /* by convension the sending partner initializes the final handshake */

  if ((*exchParams)->send)
  {
    INT4 dest = (*exchParams)->partnerProcNum;

    goodbye.length = 1;
    goodbye.data   = &magic;

    LALMPISendINT2Vector (status->statusPtr, &goodbye, dest, 
        (*exchParams)->mpiComm);
    CHECKSTATUSPTR (status);
  }
  else
  {
    INT4 source  = (*exchParams)->partnerProcNum;
    INT2 myMagic = 0;

    goodbye.length = 1;
    goodbye.data   = &myMagic;

    LALMPIRecvINT2Vector (status->statusPtr, &goodbye, source,
        (*exchParams)->mpiComm);
    CHECKSTATUSPTR (status);

    ASSERT (goodbye.data[0] == magic, status, 
            FINDCHIRPEXCH_EHAND, FINDCHIRPEXCH_MSGEHAND);
  }

  /* empty memory */
  LALFree (*exchParams);
  *exchParams = NULL;

  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
LALExchangeUINT4 (
    LALStatus     *status,
    UINT4         *object,
    ExchParams    *exchParams
    )
{
  CHARVector box; /* a box to hold some bytes of data */

  INITSTATUS (status, "ExchangeDataRequest", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  /* only do a minimal check to see if arguments are somewhat reasonable */
  ASSERT (exchParams, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (object, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);

  /* stuff the event into a box */
  box.length = sizeof (UINT4);
  box.data   = (CHAR *) object;

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


void
LALExchangeDataSegment (
    LALStatus      *status,
    DataSegment *segment,
    ExchParams  *exchParams
    )
{
  CHARVector box; /* a box to hold some bytes of data */

  INITSTATUS (status, "LALExchangeDataSegment", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  /* only do a minimal check to see if arguments are somewhat reasonable */
  ASSERT (exchParams, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (segment, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (segment->data, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (segment->spec, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (segment->resp, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);

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
    LALMPISendINT2TimeSeries (status->statusPtr, segment->data, dest,
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

    /* set relevant fields */
    segment->endOfData = tmpSegment.endOfData;
    segment->newLock   = tmpSegment.newLock;
    segment->newCal    = tmpSegment.newCal;
    segment->number    = tmpSegment.number;

    /* receive remaining fields */

    /* data */
    LALMPIRecvINT2TimeSeries (status->statusPtr, segment->data, source,
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

void
LALExchangeInspiralBankIn (
    LALStatus         *status,
    InspiralBankIn *bankIn,
    ExchParams     *exchParams
    )
{
  CHARVector box; /* a box to hold some bytes of data */

  INITSTATUS (status, "LALExchangeInspiralBankIn", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  /* only do a minimal check to see if arguments are somewhat reasonable */
  ASSERT (exchParams, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (bankIn, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);

  /* stuff the template bank input into a box */
  box.length = sizeof (InspiralBankIn);
  box.data   = (CHAR *) bankIn;

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


void
LALExchangeInspiralTemplate (
    LALStatus           *status,
    InspiralTemplate *tmplt,
    ExchParams       *exchParams
    )
{
  CHARVector box; /* a box to hold some bytes of data */

  INITSTATUS (status, "LALExchangeInspiralTemplate", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  /* only do a minimal check to see if arguments are somewhat reasonable */
  ASSERT (exchParams, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (tmplt, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);

  /* stuff the template into a box */
  box.length = sizeof (InspiralTemplate);
  box.data   = (CHAR *) tmplt;

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


void
LALExchangeInspiralEvent (
    LALStatus        *status,
    InspiralEvent *event,
    ExchParams    *exchParams
    )
{
  CHARVector box; /* a box to hold some bytes of data */

  INITSTATUS (status, "LALExchangeInspiralEvent", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  /* only do a minimal check to see if arguments are somewhat reasonable */
  ASSERT (exchParams, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);
  ASSERT (event, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);

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


void
LALExchangeInspiralEventList (
    LALStatus     *status,
    InspiralEvent **eventHead,
    ExchParams    *exchParams
    )
{
  InspiralEvent   *eventCurrent = NULL;
  InspiralEvent   *prev         = NULL;

  INITSTATUS (status, "LALExchangeInspiralEvent", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  if (exchParams->send)  /* I am sending */
  {

    /* check that we have a bank to send */
    ASSERT (*eventHead, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);

    /* exchange the template bank */
    for (eventCurrent = *eventHead;
        eventCurrent != NULL;
        eventCurrent = eventCurrent->next
        )
    {
      LALExchangeInspiralEvent (status->statusPtr, eventCurrent, exchParams);
      CHECKSTATUSPTR (status);
    }

    ASSERT (!eventCurrent, status, FINDCHIRPEXCH_ENNUL, FINDCHIRPEXCH_MSGENNUL);

  }
  else /* I am receiving */
  {

    /* check that this is a new list */
    ASSERT (!*eventHead, status, FINDCHIRPEXCH_ENNUL, FINDCHIRPEXCH_MSGENNUL);

    /* recieve the template bank */
    do
    {
      /* create memory for the template */
      eventCurrent = (InspiralEvent *) LALMalloc (sizeof(InspiralEvent));

      /* make a note of the first node in the list to return */
      if ( *eventHead == NULL ) *eventHead = eventCurrent;

      /* recieve the template */
      LALExchangeInspiralEvent (status->statusPtr, eventCurrent, exchParams);
      CHECKSTATUSPTR (status);

      /* point the previous node to this node */
      if ( prev != NULL ) prev->next = eventCurrent;
      prev = eventCurrent;

    } while ( eventCurrent->next != NULL );

  }

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

void
LALExchangeTemplateBank (
    LALStatus         *status,
    InspiralTemplate **tmpltHead,
    ExchParams        *exchParms
                 )
{
  InspiralTemplate   *tmpltCurrent = NULL;
  InspiralTemplate   *tmpltFine    = NULL;
  InspiralTemplate   *prev         = NULL;

  INITSTATUS (status, "LALExchangeTemplateBank", FINDCHIRPEXCHC);
  ATTATCHSTATUSPTR (status);

  if (exchParms->send)  /* I am sending */
  {

    /* check that we have a bank to send */
    ASSERT (*tmpltHead, status, FINDCHIRPEXCH_ENULL, FINDCHIRPEXCH_MSGENULL);

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
        LALExchangeTemplateBank (status->statusPtr, &(tmpltCurrent->fine), exchParms);
        CHECKSTATUSPTR (status);
      }
    }

    ASSERT (!tmpltCurrent, status, FINDCHIRPEXCH_ENNUL, FINDCHIRPEXCH_MSGENNUL);

  }
  else /* I am receiving */
  {

    /* check that this is a new list */
    ASSERT (!*tmpltHead, status, FINDCHIRPEXCH_ENNUL, FINDCHIRPEXCH_MSGENNUL);

    /* recieve the template bank */
    do
    {
      /* create memory for the template */
      tmpltCurrent = (InspiralTemplate *) LALMalloc (sizeof(InspiralTemplate));

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
