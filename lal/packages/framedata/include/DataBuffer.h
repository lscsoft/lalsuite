/*----------------------------------------------------------------------- 
 * 
 * File Name: DataBuffer.h
 *
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _DATABUFFER_H
#define _DATABUFFER_H

#include "LALDatatypes.h"
#include "FrameData.h"
#include "SpecBuffer.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (DATABUFFERH, "$Id$");

#define DATABUFFER_ENULL 1
#define DATABUFFER_ENNUL 2
#define DATABUFFER_ESIZE 4
#define DATABUFFER_ESZMM 8

#define DATABUFFER_MSGENULL "Null pointer"
#define DATABUFFER_MSGENNUL "Non-null pointer"
#define DATABUFFER_MSGESIZE "Invalid input size"
#define DATABUFFER_MSGESZMM "Size mismatch"

typedef struct
tagDataBlock
{
  INT4 number;
  INT4 continuous;
  INT4 anomalous;
  INT2TimeSeries *framedata;
}
DataBlock;

typedef struct
tagDataBuffer
{
  /*
   * public data
   */
  INT4 endOfData;
  INT4 newLock;
  INT4 newCal;
  /*
   * private data
   */
  FrameData          *frameData;
  SpectrumBuffer     *specBuffer;
  INT2VectorSequence *dataBuffer;
  DataBlock          *blockArray;
  INT4                blockArraySize;
  INT4                blocksFilled;
  INT4                nextData;
  INT4                nextBlock;
  INT4                deltaPoints;
  INT4                numSent;
  INT4                first;
}
DataBuffer;

typedef struct
tagDataBufferPar
{
  INT4         numSpec;
  INT4         numPoints;
  WindowType   windowType;
  RealFFTPlan *plan;
  CHAR        *framePath;
}
DataBufferPar;

typedef struct
tagDataSegment
{
  INT2TimeSeries          *data;
  REAL4FrequencySeries    *spec;
  COMPLEX8FrequencySeries *resp;
  INT4                     endOfData;
  INT4                     newLock;
  INT4                     newCal;
  INT4                     number;
}
DataSegment;

void
LALCreateDataBuffer (
    LALStatus         *status,
    DataBuffer    **buffer,
    DataBufferPar  *params
    );

void
LALDestroyDataBuffer (
    LALStatus      *status,
    DataBuffer **buffer
    );

void
LALGetData (
    LALStatus      *status,
    DataSegment *output,
    INT4         advance,
    DataBuffer  *buffer
    );


#ifdef  __cplusplus
}
#endif

#endif /* _DATABUFFER_H */
