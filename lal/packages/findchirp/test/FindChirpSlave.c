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

#ifndef _STDIO_H
#include <stdio.h>
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _FINDCHIRPEXCH_H
#include "FindChirpExch.h"
#ifndef _FINDCHIRPEXCH_H
#define _FINDCHIRPEXCH_H
#endif
#endif

NRCSID (FINDCHIRPSLAVEC, "$Id$");

/* global variables */
#include "FindChirpGlobal.h"

void
Slave (Status *status, MPIId id)
{
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

  INITSTATUS (status, FINDCHIRPSLAVEC);
  ATTATCHSTATUSPTR (status);

  printf ("Slave %d: starting up\n", id.myId);
  fflush (stdout);


  /*
   *
   * set up exchange types
   *
   */


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

    segment[i].data->name        = "anonymous";
    segment[i].data->data        = NULL;
    segment[i].data->sampleUnits = NULL;
    segment[i].spec->name        = "anonymous";
    segment[i].spec->data        = NULL;
    segment[i].spec->sampleUnits = NULL;
    segment[i].resp->name        = "anonymous";
    segment[i].resp->data        = NULL;
    segment[i].resp->sampleUnits = NULL;

    I2CreateVector (status->statusPtr, &segment[i].data->data, numPoints);
    CHECKSTATUSPTR (status);

    CreateVector (status->statusPtr, &segment[i].spec->data, numPoints/2 + 1);
    CHECKSTATUSPTR (status);

    CCreateVector (status->statusPtr, &segment[i].resp->data, numPoints/2 + 1);
    CHECKSTATUSPTR (status);
    
    CHARCreateVector (status->statusPtr, &segment[i].data->sampleUnits, 128);
    CHECKSTATUSPTR (status);

    CHARCreateVector (status->statusPtr, &segment[i].spec->sampleUnits, 128);
    CHECKSTATUSPTR (status);

    CHARCreateVector (status->statusPtr, &segment[i].resp->sampleUnits, 128);
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


  InitializeExchange (status->statusPtr, &thisExch,
                      &exchInspiralBankIn, id.myId);
  CHECKSTATUSPTR (status);

  printf ("Slave %d: sending template band input\n", id.myId);
  
  bankIn.mMin     = 1;
  bankIn.mMax     = 3;
  bankIn.ffCoarse = 0.9;
  bankIn.ffFine   = 0.99;
  bankIn.detector = Caltech40m;
  bankIn.method   = Best;

  ExchangeInspiralBankIn (status->statusPtr, &bankIn, thisExch);

  FinalizeExchange (status->statusPtr, &thisExch);
  CHECKSTATUSPTR (status);


  /*
   *
   * exchange some templates
   *
   */


  InitializeExchange (status->statusPtr, &thisExch,
                      &exchInspiralTemplates, id.myId);
  CHECKSTATUSPTR (status);

  printf ("Slave %d: receiving %d templates\n", id.myId, maxNumTemplates);

  for (i = 0; i < maxNumTemplates; ++i)
  {
    ExchangeInspiralTemplate (status->statusPtr, tmplt + i, thisExch);
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

  FinalizeExchange (status->statusPtr, &thisExch);
  CHECKSTATUSPTR (status);


  /* 
   *
   * exchange data segments 
   * 
   */


  InitializeExchange (status->statusPtr, &thisExch,
                      &exchDataSegments, id.myId);
  CHECKSTATUSPTR (status);

  printf ("Slave %d: requesting %d segments\n", id.myId, numSegments);
  fflush (stdout);

  for (i = 0; i < numSegments; ++i)
  {
    ExchangeDataSegment (status->statusPtr, segment + i, thisExch);
    CHECKSTATUSPTR (status);
    printf ("Slave %d:   got segment %d of %d\n", id.myId, i, numSegments);
    fflush (stdout);
  }

  FinalizeExchange (status->statusPtr, &thisExch);
  CHECKSTATUSPTR (status);


  /* 
   *
   * exchange 8 events
   * 
   */


  exchInspiralEvents.numObjects = 8;

  printf ("Slave %d: sending %d events\n",
          id.myId, exchInspiralEvents.numObjects);

  InitializeExchange (status->statusPtr, &thisExch,
                      &exchInspiralEvents, id.myId);
  CHECKSTATUSPTR (status);

  for (i = 0; i < thisExch->numObjects; ++i)
  {
    printf ("Slave %d:   sent event %d of %d\n",
            id.myId, i, thisExch->numObjects);
    event[i].time.gpsSeconds     = i;
    event[i].time.gpsNanoSeconds = 0;
    event[i].tmplt               = tmplt[i];
    event[i].snrsq               = 100;
    event[i].chisq               = 50;
    event[i].sigma               = 1;
    ExchangeInspiralEvent (status->statusPtr, event + i, thisExch);
    CHECKSTATUSPTR (status);
  }

  FinalizeExchange (status->statusPtr, &thisExch);
  CHECKSTATUSPTR (status);


  /* 
   *
   * finished message
   * 
   */


  printf ("Slave %d: sending finished message\n", id.myId);
  fflush (stdout);

  InitializeExchange (status->statusPtr, &thisExch, &exchFinished, id.myId);
  CHECKSTATUSPTR (status);

  FinalizeExchange (status->statusPtr, &thisExch);
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
    CHARDestroyVector (status->statusPtr, &segment[i].data->sampleUnits);
    CHECKSTATUSPTR (status);

    CHARDestroyVector (status->statusPtr, &segment[i].spec->sampleUnits);
    CHECKSTATUSPTR (status);

    CHARDestroyVector (status->statusPtr, &segment[i].resp->sampleUnits);
    CHECKSTATUSPTR (status);

    I2DestroyVector (status->statusPtr, &segment[i].data->data);
    CHECKSTATUSPTR (status);

    DestroyVector (status->statusPtr, &segment[i].spec->data);
    CHECKSTATUSPTR (status);

    CDestroyVector (status->statusPtr, &segment[i].resp->data);
    CHECKSTATUSPTR (status);
    
    LALFree (segment[i].data);
    LALFree (segment[i].spec);
    LALFree (segment[i].resp);
  }

  LALFree (segment);

  DETATCHSTATUSPTR (status);
  RETURN (status);
}

