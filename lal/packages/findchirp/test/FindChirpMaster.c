/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpMaster.c
 *
 * Author: Allen, B., Brown, D., and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/FindChirpExch.h>

NRCSID (FINDCHIRPMASTERC, "$Id$");

/* global variables */
#include "FindChirpGlobal.h"

static const INT4 numSpec = 8;

void
Master (LALStatus *status, MPIId id)
{
  CHAR                    *framePath;
  InitExchParams           initExchParams;
  DataBuffer              *buffer = NULL;
  DataBufferPar            bufferPar;
  DataSegment             *segment;
  InspiralBankIn           bankIn;
  InspiralTemplate        *tmplt;
  InspiralEvent           *event;
  INT4                     numProcs;
  INT4                     i;

  INITSTATUS (status, "Master", FINDCHIRPMASTERC);
  ATTATCHSTATUSPTR (status);
  
  printf ("Master: starting up\n");

  framePath = getenv ("LAL_FRAME_PATH");


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
   * create data buffer
   *
   */


  bufferPar.numSpec    = numSpec;
  bufferPar.numPoints  = numPoints;
  bufferPar.windowType = Welch;
  bufferPar.plan       = NULL;
  bufferPar.framePath  = framePath;
  LALEstimateFwdRealFFTPlan (status->statusPtr, &bufferPar.plan, numPoints);
  CHECKSTATUSPTR (status);

  LALCreateDataBuffer (status->statusPtr, &buffer, &bufferPar);
  CHECKSTATUSPTR (status);


  /*
   *
   * loop while some slaves are not finished
   *
   */


  numProcs = id.numProcs;
  initExchParams.mpiComm = MPI_COMM_WORLD;

  while (numProcs > 1)
  {
    ExchParams *thisExch = NULL;

    printf ("\nMaster: waiting for request... ");
    initExchParams.myProcNum = id.myId;
    LALInitializeExchange (status->statusPtr, &thisExch, NULL, &initExchParams);
    CHECKSTATUSPTR (status);
    printf ("got request %d from slave %d\n",
            (INT4) thisExch->exchObjectType, thisExch->partnerProcNum);

    switch (thisExch->exchObjectType)
    {

      /*
       *
       * exchange template bank input
       *
       */


      case ExchInspiralBankIn:

        printf ("Master: receiving inspiral bank input params from slave %d\n",
                thisExch->partnerProcNum);

        ASSERT (thisExch->numObjects == 1, status,
                1, "Wrong number objects to exchange");
        ASSERT (thisExch->send == 0, status,
                4, "Master expects to receive");

        LALExchangeInspiralBankIn (status->statusPtr, &bankIn, thisExch);
        CHECKSTATUSPTR (status);

        break;


      /*
       *
       * exchange templates
       *
       */


      case ExchInspiralTemplate:

        printf ("Master: sending 4 templates to slave %d\n",
                thisExch->partnerProcNum);

        ASSERT (thisExch->numObjects == maxNumTemplates, status,
                1, "Wrong number objects to exchange");
        ASSERT (thisExch->send == 1, status,
                2, "Master expects to send");

        /* the four templates */
        for (i = 0; i < 4; ++i)
        {
          printf ("Master:   sent template %d of 4\n", i);
          tmplt[i].number = i;
          tmplt[i].mass1  = 1 + (i + 1)/2;
          tmplt[i].mass2  = 1 + i/2;

          LALExchangeInspiralTemplate (status->statusPtr, tmplt + i, thisExch);
          CHECKSTATUSPTR (status);
        }

        /* one more to say we're done */
        printf ("Master:   sent template 4 -- done\n");
        tmplt[4].number = -1;
        LALExchangeInspiralTemplate (status->statusPtr, tmplt + 4, thisExch);
        CHECKSTATUSPTR (status);

        break;


      /*
       *
       * exchange data segments
       *
       */


      case ExchDataSegment:

        printf ("Master: sending %d segments to slave %d\n",
                numSegments, thisExch->partnerProcNum);

        ASSERT (thisExch->numObjects == numSegments, status,
                1, "Wrong number of objects to exchange");
        ASSERT (thisExch->send, status, 2, "Master expects to send");

        for (i = 0; i < numSegments; ++i)
        {
          LALGetData (status->statusPtr, segment + i, 3*numPoints/4, buffer);
          CHECKSTATUSPTR (status);

          printf ("Master:  sent segment %d of %d\n", i, numSegments);
          LALExchangeDataSegment (status->statusPtr, segment + i, thisExch);
          CHECKSTATUSPTR (status);
        }

        break;


      /*
       *
       * exchange events
       *
       */


      case ExchInspiralEvent:

        printf ("Master: receiving %d events from slave %d\n",
                thisExch->numObjects, thisExch->partnerProcNum);

        ASSERT (thisExch->numObjects <= numPoints, status,
                1, "Wrong number of objects to exchange");
        ASSERT (thisExch->send == 0, status,
                4, "Master expects to receive");

        for (i = 0; i < thisExch->numObjects; ++i)
        {
          printf ("Master:   got event %d of %d\n", i, thisExch->numObjects);
          LALExchangeInspiralEvent (status->statusPtr, event + i, thisExch);
          CHECKSTATUSPTR (status);
        }

        break;


      /*
       *
       * slave is finished
       *
       */


      case ExchFinished:

        printf ("Master: received finished message from slave %d\n",
                thisExch->partnerProcNum);

        --numProcs;
        break;


      /*
       *
       * unrecognized exchange type
       *
       */


      default:
        ASSERT (0, status, 2, "Unrecognized exchange type");

    }

    LALFinalizeExchange (status->statusPtr, &thisExch);
    printf ("Master: done with request\n");
  }


  /*
   *
   * all slaves are finished: shut down
   *
   */

  
  printf ("Master: shutting down\n");


  /*
   *
   * free memory
   *
   */


  LALDestroyRealFFTPlan (status->statusPtr, &bufferPar.plan);
  CHECKSTATUSPTR (status);

  LALDestroyDataBuffer (status->statusPtr, &buffer);
  CHECKSTATUSPTR (status);

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
