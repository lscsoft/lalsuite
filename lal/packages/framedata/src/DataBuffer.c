/*----------------------------------------------------------------------- 
 * 
 * File Name: DataBuffer.c
 *
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include "LALStdlib.h"
#include "SeqFactories.h"
#include "SpecBuffer.h"
#include "DataBuffer.h"

NRCSID (DATABUFFERC, "$Id$");

static void
FillDataBlock (
    Status     *status,
    DataBuffer *buffer,
    INT4        whichBlock
    )
{
  DataBlock    *theBlock;

  INITSTATUS (status, "FillDataBlock", DATABUFFERC);
  ATTATCHSTATUSPTR (status);

  theBlock = buffer->blockArray + whichBlock;

  /* get next block of locked data */
  theBlock->continuous = 1;
  theBlock->anomalous  = 0;
  theBlock->number     = buffer->blocksFilled++;

  GetFrameData (status->statusPtr, theBlock->framedata, buffer->frameData);
  CHECKSTATUSPTR (status);
  if (buffer->frameData->endOfData) /* end of data */
  {
    buffer->endOfData    = 1;
    theBlock->continuous = 0;
    DETATCHSTATUSPTR (status);
    RETURN (status);
  }
  if (buffer->frameData->newCalibration) /* new calib info pending */
  {
    buffer->newCal = 1;
  }

  while (buffer->frameData->newLock)
  {
    INT2TimeSeries dummy;
    INT2Vector     dumvec;
    theBlock->continuous    = 0;
    theBlock->anomalous     = 1;
    buffer->newLock         = 1; /* indicate a new lock is pending */
    dumvec.length           = 3*60/theBlock->framedata->deltaT; /* 3 minutes */
    dumvec.data             = NULL;                             /* seek mode */
    dummy.data              = &dumvec;
    GetFrameData (status->statusPtr, theBlock->framedata, buffer->frameData);
    CHECKSTATUSPTR (status);
    if (buffer->frameData->endOfData) /* end of data */
    {
      buffer->endOfData = 1;
      DETATCHSTATUSPTR (status);
      RETURN (status);
    }
    if (buffer->frameData->newCalibration) /* new calib info pending */
    {
      buffer->newCal = 1;
    }
  }

  /* if the block is not anomalous, add it to the spectrum buffer */
  if (!theBlock->anomalous)
  {
    AddSpectrum (status->statusPtr, buffer->specBuffer, theBlock->framedata);
    CHECKSTATUSPTR (status);
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


static void
FillAllDataBlocks (
    Status     *status,
    DataBuffer *buffer
    )
{
  INT4 block;

  INITSTATUS (status, "FillAllDataBlocks", DATABUFFERC);
  ATTATCHSTATUSPTR (status);

  buffer->newLock      = 0; /* new lock no longer pending */

  for (block = 0; block < buffer->blockArraySize; ++block)
  {
    FillDataBlock (status->statusPtr, buffer, block);
    CHECKSTATUSPTR (status);

    if (buffer->endOfData) /* end of data: immediate exit */
    {
      DETATCHSTATUSPTR (status);
      RETURN (status);
    }

    if (buffer->newLock) /* need to refill from beginning */
    {
      block = -1;
      buffer->newLock = 0;
    }
  }

  buffer->nextData     = 0;
  buffer->nextBlock    = 0;
  buffer->deltaPoints  = 0;
  buffer->newLock      = 0; /* new lock no longer pending */
  
  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
CreateDataBuffer (
    Status         *status,
    DataBuffer    **buffer,
    DataBufferPar  *params
    )
{
  SpectrumBufferPar      specpar;
  CreateVectorSequenceIn vseqin;
  INT4 block;

  INITSTATUS (status, "CreateDataBuffer", DATABUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (buffer, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (params, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);

  /* make sure that the buffer has not already been assigned */
  ASSERT (*buffer == NULL, status, DATABUFFER_ENNUL, DATABUFFER_MSGENNUL);

  /* assign memory for buffer */
  *buffer = (DataBuffer *) LALMalloc (sizeof (DataBuffer));

  /* make sure that the allocation was successful */
  ASSERT (*buffer, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);

  (*buffer)->blockArraySize = params->numSpec;
  vseqin.length             = params->numSpec + 1; /* extra swap space */
  vseqin.vectorLength       = params->numPoints;

  /* create buffer */
  (*buffer)->dataBuffer = NULL;
  I2CreateVectorSequence (status->statusPtr, &(*buffer)->dataBuffer, &vseqin);
  CHECKSTATUSPTR (status);

  /* create block array */
  (*buffer)->blockArray = (DataBlock *)
    LALMalloc ((*buffer)->blockArraySize*sizeof(DataBlock));
  ASSERT ((*buffer)->blockArray, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);

  for (block = 0; block < (*buffer)->blockArraySize; ++block)
  {
    DataBlock *thisBlock = (*buffer)->blockArray + block;

    /* allocate memory for the framedata */
    thisBlock->framedata = (INT2TimeSeries *)
      LALMalloc (sizeof (INT2TimeSeries));
    ASSERT (thisBlock->framedata, status,
            DATABUFFER_ENULL, DATABUFFER_MSGENULL);
    thisBlock->framedata->data = (INT2Vector *)LALMalloc (sizeof(INT2Vector));
    ASSERT (thisBlock->framedata->data, status,
            DATABUFFER_ENULL, DATABUFFER_MSGENULL);
    /* should create units vector ... */
    thisBlock->framedata->data->length = params->numPoints;
    thisBlock->framedata->data->data   = (*buffer)->dataBuffer->data +
      block*params->numPoints;
  }

  /* create spectrum buffer */
  (*buffer)->specBuffer = NULL;
  specpar.numSpec       = params->numSpec;
  specpar.numPoints     = params->numPoints;
  specpar.windowType    = params->windowType;
  specpar.plan          = params->plan;
  CreateSpectrumBuffer (status->statusPtr, &(*buffer)->specBuffer, &specpar);
  CHECKSTATUSPTR (status);

  /* initialize frame data */
  (*buffer)->frameData = NULL;
  InitializeFrameData (
      status->statusPtr,
      &(*buffer)->frameData,
      params->framePath
      );
  CHECKSTATUSPTR (status);

  (*buffer)->blocksFilled = 0;
  (*buffer)->endOfData    = 0;
  (*buffer)->newCal       = 0;
  (*buffer)->numSent      = 0;
  (*buffer)->first        = 1;

  /* fill all data blocks */
  FillAllDataBlocks (status->statusPtr, *buffer);
  CHECKSTATUSPTR (status);

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
DestroyDataBuffer (
    Status      *status,
    DataBuffer **buffer
    )
{
  INT4 block;

  INITSTATUS (status, "DestroyDataBuffer", DATABUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (buffer, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (*buffer, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);

  /* finalize frame data */
  FinalizeFrameData (status->statusPtr, &(*buffer)->frameData);
  CHECKSTATUSPTR (status);

  /* destroy spectrum buffer */
  DestroySpectrumBuffer (status->statusPtr, &(*buffer)->specBuffer);
  CHECKSTATUSPTR (status);

  /* destroy data buffer */
  I2DestroyVectorSequence (status->statusPtr, &(*buffer)->dataBuffer);
  CHECKSTATUSPTR (status);

  for (block = 0; block < (*buffer)->blockArraySize; ++block)
  {
    DataBlock *thisBlock = (*buffer)->blockArray + block;
    LALFree (thisBlock->framedata->data);
    LALFree (thisBlock->framedata);
  }

  LALFree ((*buffer)->blockArray);
  LALFree (*buffer);
  *buffer = NULL;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


void
GetData (
    Status      *status,
    DataSegment *output,
    INT4         advance,
    DataBuffer  *buffer
    )
{
  INT4 numPoints;
  INT4 numBlocks;

  INITSTATUS (status, "GetData", DATABUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (output, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (buffer, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (output->data, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (output->spec, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (output->resp, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (output->data->data, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (output->spec->data, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (output->resp->data, status, DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (output->data->data->data, status,
          DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (output->spec->data->data, status,
          DATABUFFER_ENULL, DATABUFFER_MSGENULL);
  ASSERT (output->resp->data->data, status,
          DATABUFFER_ENULL, DATABUFFER_MSGENULL);

  /* check sizes ... */
  numPoints = output->data->data->length;
  ASSERT (numPoints > 0, status, DATABUFFER_ESIZE, DATABUFFER_MSGESIZE);
  ASSERT (buffer->dataBuffer->vectorLength == numPoints, status,
          DATABUFFER_ESZMM, DATABUFFER_MSGESZMM);
  ASSERT (output->resp->data->length == numPoints/2 + 1, status,
          DATABUFFER_ESZMM, DATABUFFER_MSGESZMM);
  ASSERT (output->spec->data->length == numPoints/2 + 1, status,
          DATABUFFER_ESZMM, DATABUFFER_MSGESZMM);

  /* make sure advance is reasonable */
  ASSERT (advance > 0, status, DATABUFFER_ESIZE, DATABUFFER_MSGESIZE);
  ASSERT (advance < numPoints, status, DATABUFFER_ESIZE, DATABUFFER_MSGESIZE);

  if (buffer->endOfData) /* no more data: immediate exit */
  {
    output->endOfData = 1;
    DETATCHSTATUSPTR (status);
    RETURN (status);
  }

  numBlocks = buffer->blockArraySize;

  if (buffer->first)
  {
    output->newLock = 1;
    buffer->first   = 0;
  }
  else
  {
    output->newLock = 0;
  }
  output->newCal = 0;

  /* 
   * does next block contain continuous data?  (only need to check if
   * buffer->nextData is not at the beginning of a block)
   */
  if (buffer->nextData%numPoints > 0)
  {
    INT4       nextBlockNo = (buffer->nextData/numPoints + 1)%numBlocks;
    DataBlock *nextBlock   = buffer->blockArray + nextBlockNo;
    INT4       continuous  = nextBlock->continuous;

    if (!continuous) /* no: skip into next block */
    {
      output->newLock = 1;

      /* refill buffer */
      FillAllDataBlocks (status->statusPtr, buffer);
      CHECKSTATUSPTR (status);

      if (buffer->endOfData)
      {
        output->endOfData = 1;
        DETATCHSTATUSPTR (status);
        RETURN (status);
      }
    }
  }

  /* 
   * do we need to fill another block of data? (only do so if a
   * buffer->newLock flag is not pending)
   */
  if (buffer->deltaPoints > numBlocks*numPoints/2 && !buffer->newLock)
  {
    FillDataBlock (status->statusPtr, buffer, buffer->nextBlock);
    CHECKSTATUSPTR (status);
    if (buffer->endOfData)
    {
      output->endOfData = 1;
      DETATCHSTATUSPTR (status);
      RETURN (status);
    }

    buffer->deltaPoints -= numPoints;

    if (buffer->nextBlock == 0) /* copy data to end of buffer */
    {
      INT2 *first = buffer->dataBuffer->data;
      INT2 *last  = first + numBlocks*numPoints;
      memcpy (last, first, numPoints*sizeof(INT2));
    }

    buffer->nextBlock = (buffer->nextBlock + 1)%numBlocks;
  }

  /* copy the data to the output vector */
  memcpy (
      output->data->data->data,
      buffer->dataBuffer->data + buffer->nextData,
      numPoints*sizeof(INT2)
      );

  /* set epoch, etc ... */
  {
    INT4            blockno = buffer->nextData/numPoints;
    INT2TimeSeries *tseries = (buffer->blockArray + blockno)->framedata;
    REAL4           time;
    time  = tseries->epoch.gpsSeconds + 1e-9*tseries->epoch.gpsNanoSeconds;
    time += (buffer->nextData%numPoints)*tseries->deltaT;
    
    output->data->epoch.gpsSeconds     = time;
    output->data->epoch.gpsNanoSeconds = fmod (1e9*time, 1e9);
    output->data->deltaT = tseries->deltaT;
    output->data->f0     = tseries->f0;
  }

  /* get average spectrum */
  AverageSpectrum (status->statusPtr, output->spec, buffer->specBuffer);
  CHECKSTATUSPTR (status);

  /* now is the time to deal with new calibration if pending */
  if (buffer->newCal && output->newLock)
  {
    output->newCal = 1;

    /* get new calibration data */
    GetFrameDataResponse (
        status->statusPtr,
        output->resp,
        buffer->frameData
        );
    CHECKSTATUSPTR (status);

    buffer->newCal = 0;
  }

  output->endOfData = 0;
  output->number    = buffer->numSent++;

  /* bookkeeping */
  buffer->deltaPoints += advance;
  buffer->nextData    += advance;

  /* if next data is in end of buffer, move it back to the beginning */
  if (buffer->nextData/numPoints == numBlocks)
    buffer->nextData = buffer->nextData%numPoints;

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}
