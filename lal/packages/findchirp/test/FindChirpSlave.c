/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpSlave.c
 *
 * Author: Allen, B., Brown, D., and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/FindChirpExch.h>

NRCSID (FINDCHIRPSLAVEC, "$Id$");

/* global variables */
#include "FindChirpGlobal.h"

void
Slave (LALStatus *status, MPIId id)
{
  InitExchParams    initExchParams;
  ExchParams        exchDataSegments;
  ExchParams        exchInspiralBankIn;
  ExchParams        exchInspiralTemplates;
  ExchParams        exchInspiralEvents;
  ExchParams        exchFinished;
  ExchParams       *thisExch = NULL;
  DataSegment      *segment;
  InspiralBankIn    bankIn;
  InspiralTemplate *tmplt;
  InspiralEvent    *event;
  INT4              i;

  INITSTATUS (status, "Slave", FINDCHIRPSLAVEC);
  ATTATCHSTATUSPTR (status);

  printf ("Slave %d: starting up\n", id.myId);
  fflush (stdout);


  /*
   *
   * set up exchange types
   *
   */

  
  /* init exchane params */
  initExchParams.mpiComm   = MPI_COMM_WORLD;
  initExchParams.myProcNum = id.myId;

  /* exchange data segments */
  exchDataSegments.exchObjectType      = ExchDataSegment;
  exchDataSegments.send                = 0; /* master sends */
  exchDataSegments.numObjects          = numSegments;
  exchDataSegments.partnerProcNum      = 0; /* master */

  /* exchange inspiral template bank input */
  exchInspiralBankIn.exchObjectType    = ExchInspiralBankIn;
  exchInspiralBankIn.send              = 1; /* I send */
  exchInspiralBankIn.numObjects        = 1;
  exchInspiralBankIn.partnerProcNum    = 0; /* master */

  /* perhaps master should initiate reverse-exchange? */
  /* exchange inspiral templates */
  exchInspiralTemplates.exchObjectType = ExchInspiralTemplate;
  exchInspiralTemplates.send           = 0; /* master sends */
  exchInspiralTemplates.numObjects     = maxNumTemplates;
  exchInspiralTemplates.partnerProcNum = 0; /* master */

  /* exchange inspiral events */
  exchInspiralEvents.exchObjectType    = ExchInspiralEvent;
  exchInspiralEvents.send              = 1; /* I send */
  /* actually, this next can be adjusted at will */
  exchInspiralEvents.numObjects        = numPoints; /* I send */
  exchInspiralEvents.partnerProcNum    = 0; /* master */

  /* exchange finished message */
  exchFinished.exchObjectType         = ExchFinished;
  exchFinished.send                   = 0; /* master sends */
  exchFinished.numObjects             = 0; /* irrelevant */
  exchFinished.partnerProcNum         = 0; /* master */


  /*
   *
   * allocate memory
   *
   */


  /* allocate memory for segments */

  segment = (DataSegment *) LALMalloc (numSegments*sizeof(DataSegment));

  for (i = 0; i < numSegments; ++i)
  {
    segment[i].data =
      (INT2TimeSeries *) LALMalloc (sizeof(INT2TimeSeries));
    segment[i].spec =
      (REAL4FrequencySeries *) LALMalloc (sizeof(REAL4FrequencySeries));
    segment[i].resp =
      (COMPLEX8FrequencySeries *) LALMalloc (sizeof(COMPLEX8FrequencySeries));

    segment[i].data->data = NULL;
    segment[i].spec->data = NULL;
    segment[i].resp->data = NULL;

    LALI2CreateVector (status->statusPtr, &segment[i].data->data, numPoints);
    CHECKSTATUSPTR (status);

    LALCreateVector (status->statusPtr, &segment[i].spec->data, numPoints/2 + 1);
    CHECKSTATUSPTR (status);

    LALCCreateVector (status->statusPtr, &segment[i].resp->data, numPoints/2 + 1);
    CHECKSTATUSPTR (status);
  }

  /* allocate (cleared) memory for templates */
  tmplt = (InspiralTemplate *)
    LALCalloc (maxNumTemplates, sizeof(InspiralTemplate));

  /* allocate (cleared) memory for events */
  event = (InspiralEvent *) LALCalloc (numPoints, sizeof(InspiralEvent));


  /*
   *
   * start doing exchanges
   *
   */


  /*
   *
   * exchange template bank input
   *
   */


  LALInitializeExchange (status->statusPtr, &thisExch,
                      &exchInspiralBankIn, &initExchParams);
  CHECKSTATUSPTR (status);

  printf ("Slave %d: sending template band input\n", id.myId);
  
  bankIn.mMin     = 1;
  bankIn.mMax     = 3;
  bankIn.ffCoarse = 0.9;
  bankIn.ffFine   = 0.99;
  bankIn.detector = Caltech40m;
  bankIn.method   = best;

  LALExchangeInspiralBankIn (status->statusPtr, &bankIn, thisExch);

  LALFinalizeExchange (status->statusPtr, &thisExch);
  CHECKSTATUSPTR (status);


  /*
   *
   * exchange some templates
   *
   */


  LALInitializeExchange (status->statusPtr, &thisExch,
                      &exchInspiralTemplates, &initExchParams);
  CHECKSTATUSPTR (status);

  printf ("Slave %d: receiving %d templates\n", id.myId, maxNumTemplates);

  for (i = 0; i < maxNumTemplates; ++i)
  {
    LALExchangeInspiralTemplate (status->statusPtr, tmplt + i, thisExch);
    CHECKSTATUSPTR (status);

    printf ("Slave %d:   got template %d of %d\n", id.myId, i, maxNumTemplates);
    fflush (stdout);

    /* see if template number is negative: signals end of templates */
    if (tmplt[i].number < 0)
    {
      printf ("Slave %d:   template %d indicates end of templates\n",
              id.myId, i);
      fflush (stdout);
      break;
    }
  }

  LALFinalizeExchange (status->statusPtr, &thisExch);
  CHECKSTATUSPTR (status);


  /* 
   *
   * exchange data segments 
   * 
   */


  LALInitializeExchange (status->statusPtr, &thisExch,
                      &exchDataSegments, &initExchParams);
  CHECKSTATUSPTR (status);

  printf ("Slave %d: requesting %d segments\n", id.myId, numSegments);
  fflush (stdout);

  for (i = 0; i < numSegments; ++i)
  {
    LALExchangeDataSegment (status->statusPtr, segment + i, thisExch);
    CHECKSTATUSPTR (status);
    printf ("Slave %d:   got segment %d of %d\n", id.myId, i, numSegments);
    fflush (stdout);
  }

  LALFinalizeExchange (status->statusPtr, &thisExch);
  CHECKSTATUSPTR (status);


  /* 
   *
   * exchange 8 events
   * 
   */


  exchInspiralEvents.numObjects = 8;

  printf ("Slave %d: sending %d events\n",
          id.myId, exchInspiralEvents.numObjects);

  LALInitializeExchange (status->statusPtr, &thisExch,
                      &exchInspiralEvents, &initExchParams);
  CHECKSTATUSPTR (status);

  for (i = 0; i < thisExch->numObjects; ++i)
  {
    printf ("Slave %d:   sent event %d of %d\n",
            id.myId, i, thisExch->numObjects);
    event[i].tmplt               = tmplt[i];
    event[i].snrsq               = 100;
    event[i].chisq               = 50;
    event[i].sigma               = 1;
    LALExchangeInspiralEvent (status->statusPtr, event + i, thisExch);
    CHECKSTATUSPTR (status);
  }

  LALFinalizeExchange (status->statusPtr, &thisExch);
  CHECKSTATUSPTR (status);


  /* 
   *
   * finished message
   * 
   */


  printf ("Slave %d: sending finished message\n", id.myId);
  fflush (stdout);

  LALInitializeExchange (status->statusPtr, &thisExch, &exchFinished, 
      &initExchParams);
  CHECKSTATUSPTR (status);

  LALFinalizeExchange (status->statusPtr, &thisExch);
  CHECKSTATUSPTR (status);

  printf ("Slave %d: shutting down\n", id.myId);
  fflush (stdout);


  /*
   *
   * free memory
   * 
   */


  LALFree (event);
  LALFree (tmplt);

  for (i = 0; i < numSegments; ++i)
  {
    LALI2DestroyVector (status->statusPtr, &segment[i].data->data);
    CHECKSTATUSPTR (status);

    LALDestroyVector (status->statusPtr, &segment[i].spec->data);
    CHECKSTATUSPTR (status);

    LALCDestroyVector (status->statusPtr, &segment[i].resp->data);
    CHECKSTATUSPTR (status);
    
    LALFree (segment[i].data);
    LALFree (segment[i].spec);
    LALFree (segment[i].resp);
  }

  LALFree (segment);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

