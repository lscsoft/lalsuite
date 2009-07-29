/*
*  Copyright (C) 2007 Jolien Creighton
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#if 0  /* autodoc block */

<lalVerbatim file="DataBufferCV">
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{DataBuffer.c}}
\label{ss:DataBuffer.c}

Routines for getting overlapping segments of data along with estimated spectra
and calibration information.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DataBufferCP}
\idx{LALCreateDataBuffer()}
\idx{LALDestroyDataBuffer()}
\idx{LALGetData()}

\subsubsection*{Description}

The routine \texttt{LALCreateDataBuffer()} sets up a buffer for holding
(overlapping) segments of data.  The parameters for this routine are the
number of spectra to average, the number of points in each segment, the type
of window to use in computing the spectra, the forward real FFT plan, and the
directory path of the frame data.  When the user is finished with a data
buffer, it is to be destroyed using \texttt{LALDestroyDataBuffer()}.

The routine \texttt{LALGetData()} acquires the next segment of data, its
spectrum, and its response function, as well as sets flags indicating if the
end of the data has been reached, if a new locked segment is being entered,
and if the response function has been changed.  The data in the buffer is
advanced by an amount \texttt{advance} specified as input.

\subsubsection*{Operating Instructions}

\begin{verbatim}
const  UINT4                   numSpec     = 8;
const  UINT4                   numPoints   = 1024;
CHAR                           framePath[] = "/data/frames"
static LALStatus               status;
static DataBufferPar           params;
static DataBuffer             *buffer;
static DataSegment             segmnt;
static INT2TimeSeries          dmro;
static REAL4FrequencySeries    spec;
static COMPLEX8FrequencySeries resp;

params.numSpec    = numSpec;
params.numPoints  = numPoints;
params.windowType = Welch;
params.framePath  = framePath;
LALCreateForwardRealFFTPlan( &status, &params.plan, numPoints, 0 );
LALCreateDataBuffer( &status, &buffer, &params );
LALI2CreateVector( &status, &dmro.data, numPoints );
LALSCreateVector( &status, &spec.data, numPoints/2 + 1 );
LALCCreateVector( &status, &resp.data, numPoints/2 + 1 );
segmnt.data = &dmro;
segmnt.spec = &spec;
segmnt.resp = &resp;

/* enter infinite loop */
while ( 1 )
{
  /* get next data segment, overlapping by 1/4 of a segment */
  LALGetData( &status, &segmnt, 3*numPoints/4, buffer );

  /* break out of loop when end of data is reached */
  if ( segmnt.endOfData )
  {
    break;
  }

  if ( segmnt.newLock )
  {
    /* new lock acquired */
  }

  if ( segmnt.newCal )
  {
    /* new calibration data */
  }

}

LALCDestroyVector( &status, &resp.data );
LALSDestroyVector( &status, &spec.data );
LALI2DestroyVector( &status, &dmro.data );
LALDestroyDataBuffer( &status, &buffer );
LALDestroyRealFFTPlan( &status, &params.plan );
\end{verbatim}

\subsubsection*{Algorithm}

\subsubsection*{Uses}

\subsubsection*{Notes}
\vfill{\footnotesize\input{DataBufferCV}}

</lalLaTeX>

#endif /* autodoc block */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>
#include <lal/SpecBuffer.h>
#include <lal/DataBuffer.h>

NRCSID (DATABUFFERC, "$Id$");

static void
FillDataBlock (
    LALStatus     *status,
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

  LALGetFrameData (status->statusPtr, theBlock->framedata, buffer->frameData);
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
    LALGetFrameData (status->statusPtr, theBlock->framedata, buffer->frameData);
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
    LALAddSpectrum (status->statusPtr, buffer->specBuffer, theBlock->framedata);
    CHECKSTATUSPTR (status);
  }

  /* normal exit */
  DETATCHSTATUSPTR (status);
  RETURN (status);
}


static void
FillAllDataBlocks (
    LALStatus     *status,
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


/* <lalVerbatim file="DataBufferCP"> */
void
LALCreateDataBuffer (
    LALStatus         *status,
    DataBuffer    **buffer,
    DataBufferPar  *params
    )
{ /* </lalVerbatim> */
  SpectrumBufferPar      specpar;
  CreateVectorSequenceIn vseqin;
  INT4 block;

  INITSTATUS (status, "LALCreateDataBuffer", DATABUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (buffer, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (params, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);

  /* make sure that the buffer has not already been assigned */
  ASSERT (*buffer == NULL, status, DATABUFFERH_ENNUL, DATABUFFERH_MSGENNUL);

  /* assign memory for buffer */
  *buffer = (DataBuffer *) LALMalloc (sizeof (DataBuffer));

  /* make sure that the allocation was successful */
  ASSERT (*buffer, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);

  (*buffer)->blockArraySize = params->numSpec;
  vseqin.length             = params->numSpec + 1; /* extra swap space */
  vseqin.vectorLength       = params->numPoints;

  /* create buffer */
  (*buffer)->dataBuffer = NULL;
  LALI2CreateVectorSequence (status->statusPtr, &(*buffer)->dataBuffer, &vseqin);
  CHECKSTATUSPTR (status);

  /* create block array */
  (*buffer)->blockArray = (DataBlock *)
    LALMalloc ((*buffer)->blockArraySize*sizeof(DataBlock));
  ASSERT ((*buffer)->blockArray, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);

  for (block = 0; block < (*buffer)->blockArraySize; ++block)
  {
    DataBlock *thisBlock = (*buffer)->blockArray + block;

    /* allocate memory for the framedata */
    thisBlock->framedata = (INT2TimeSeries *)
      LALMalloc (sizeof (INT2TimeSeries));
    ASSERT (thisBlock->framedata, status,
            DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
    thisBlock->framedata->data = (INT2Vector *)LALMalloc (sizeof(INT2Vector));
    ASSERT (thisBlock->framedata->data, status,
            DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
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
  LALCreateSpectrumBuffer (status->statusPtr, &(*buffer)->specBuffer, &specpar);
  CHECKSTATUSPTR (status);

  /* initialize frame data */
  (*buffer)->frameData = NULL;
  LALInitializeFrameData (
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


/* <lalVerbatim file="DataBufferCP"> */
void
LALDestroyDataBuffer (
    LALStatus      *status,
    DataBuffer **buffer
    )
{ /* </lalVerbatim> */
  INT4 block;

  INITSTATUS (status, "LALDestroyDataBuffer", DATABUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (buffer, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (*buffer, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);

  /* finalize frame data */
  LALFinalizeFrameData (status->statusPtr, &(*buffer)->frameData);
  CHECKSTATUSPTR (status);

  /* destroy spectrum buffer */
  LALDestroySpectrumBuffer (status->statusPtr, &(*buffer)->specBuffer);
  CHECKSTATUSPTR (status);

  /* destroy data buffer */
  LALI2DestroyVectorSequence (status->statusPtr, &(*buffer)->dataBuffer);
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


/* <lalVerbatim file="DataBufferCP"> */
void
LALGetData (
    LALStatus      *status,
    DataSegment *output,
    INT4         advance,
    DataBuffer  *buffer
    )
{ /* </lalVerbatim> */
  INT4 numPoints;
  INT4 numBlocks;

  INITSTATUS (status, "LALGetData", DATABUFFERC);
  ATTATCHSTATUSPTR (status);

  /* make sure that arguments are not NULL */
  ASSERT (output, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (buffer, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (output->data, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (output->spec, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (output->resp, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (output->data->data, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (output->spec->data, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (output->resp->data, status, DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (output->data->data->data, status,
          DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (output->spec->data->data, status,
          DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);
  ASSERT (output->resp->data->data, status,
          DATABUFFERH_ENULL, DATABUFFERH_MSGENULL);

  /* check sizes ... */
  numPoints = output->data->data->length;
  ASSERT (numPoints > 0, status, DATABUFFERH_ESIZE, DATABUFFERH_MSGESIZE);
  ASSERT ((INT4)buffer->dataBuffer->vectorLength == numPoints, status,
          DATABUFFERH_ESZMM, DATABUFFERH_MSGESZMM);
  ASSERT ((INT4)output->resp->data->length == numPoints/2 + 1, status,
          DATABUFFERH_ESZMM, DATABUFFERH_MSGESZMM);
  ASSERT ((INT4)output->spec->data->length == numPoints/2 + 1, status,
          DATABUFFERH_ESZMM, DATABUFFERH_MSGESZMM);

  /* make sure advance is reasonable */
  ASSERT (advance > 0, status, DATABUFFERH_ESIZE, DATABUFFERH_MSGESIZE);
  ASSERT (advance < numPoints, status, DATABUFFERH_ESIZE, DATABUFFERH_MSGESIZE);

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
    INT8            timeNS;
    timeNS  = 1000000000L*(INT8)(tseries->epoch.gpsSeconds);
    timeNS += tseries->epoch.gpsNanoSeconds;
    timeNS += (INT8)(1e9*(buffer->nextData%numPoints)*tseries->deltaT);
    output->data->epoch.gpsSeconds     = (INT4)(timeNS/1000000000L);
    output->data->epoch.gpsNanoSeconds = (INT4)(timeNS%1000000000L);
    output->data->deltaT = tseries->deltaT;
    output->data->f0     = tseries->f0;
  }

  /* get average spectrum */
  LALAverageSpectrum (status->statusPtr, output->spec, buffer->specBuffer);
  CHECKSTATUSPTR (status);

  /* now is the time to deal with new calibration if pending */
  if (buffer->newCal && output->newLock)
  {
    output->newCal = 1;

    /* get new calibration data */
    LALGetFrameDataResponse (
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
