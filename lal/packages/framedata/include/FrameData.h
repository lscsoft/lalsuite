/*----------------------------------------------------------------------- 
 * 
 * File Name: FrameData.h
 *
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _FRAMEDATA_H
#define _FRAMEDATA_H

#include "LALDatatypes.h"

#ifdef  __cplusplus
extern "C" {
#endif


NRCSID (FRAMEDATAH, "$Id$");

#define FRAMEDATA_ENULL 1
#define FRAMEDATA_ENNUL 2
#define FRAMEDATA_EREAD 4
#define FRAMEDATA_EOPEN 8
#define FRAMEDATA_ENOSS 16
#define FRAMEDATA_EDMRO 32
#define FRAMEDATA_ELOCK 64
#define FRAMEDATA_ELOHI 128
#define FRAMEDATA_ESIZE 256
#define FRAMEDATA_ESSSZ 512

#define FRAMEDATA_MSGENULL "Null pointer"
#define FRAMEDATA_MSGENNUL "Non-null pointer"
#define FRAMEDATA_MSGEREAD "Error reading frame directory"
#define FRAMEDATA_MSGEOPEN "Error opening frame file"
#define FRAMEDATA_MSGENOSS "No sweptsine calibration data in frame"
#define FRAMEDATA_MSGEDMRO "No IFO_DMRO data in frame"
#define FRAMEDATA_MSGELOCK "No IFO_Lock data in frame"
#define FRAMEDATA_MSGELOHI "No locklo/lockhi data in frame"
#define FRAMEDATA_MSGESIZE "Invalid vector length"
#define FRAMEDATA_MSGESSSZ "Bad sweptsine calibration data"

typedef struct
tagFrameData
{
  /*
   *
   * public data:
   *
   */
  INT4 inLock; /* data aquisition mode */
  /* return status flags */
  INT4 endOfData;
  INT4 newLock;
  INT4 newCalibration;
  /*
   *
   * private data:
   *
   */
  INT4                fileOpen;
  INT4                numFiles;
  INT4                fileNum;
  CHARVectorSequence *frameFileNames;
  void               *frameFile;           /* type (struct FrFile *)     */
  void               *frame;               /* type (struct FrameH *)     */
  void               *dmro;                /* type (struct FrAdcData *)  */
  void               *lock;                /* type (struct FrAdcData *)  */
  void               *lockLowHigh;         /* type (struct FrStatData *) */
  void               *calibration;         /* type (struct FrStatData *) */
  LIGOTimeGPS         frameStartTime;
  LIGOTimeGPS         frameFinishTime;
  LIGOTimeGPS         calibrationTime;
  INT4                numDmro;
  INT4                curDmro;
  INT4                numLock;
  INT4                curLock;
  INT4                ratio;
  REAL8               sampleRate;
  INT2                lockLow;
  INT2                lockHigh;
  INT4                dataBreak;
}
FrameData;

void
InitializeFrameData (
    Status     *status,
    FrameData **frameData,
    CHAR       *framePath
    );

void
FinalizeFrameData (
    Status     *status,
    FrameData **frameData
    );

void
GetFrameData (
    Status         *status,
    INT2TimeSeries *data,
    FrameData      *frameData
    );

void
GetFrameDataResponse (
    Status                  *status,
    COMPLEX8FrequencySeries *response,
    FrameData               *frameData
    );


#ifdef  __cplusplus
}
#endif

#endif /* _FRAMEDATA_H */
