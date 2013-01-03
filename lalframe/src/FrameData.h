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
#ifndef _FRAMEDATA_H
#define _FRAMEDATA_H

#include <lal/LALDatatypes.h>

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * \defgroup FrameData_h Header FrameData.h
 * \ingroup pkg_framedata
 *
 * \brief Root finding routines.
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/FrameData.h>
 * \endcode
 *
 * Gets IFO_DMRO data from frame files.
 */
/*@{*/

/** \name Error Codes */
/*@{*/
#define FRAMEDATAH_ENULL 1
#define FRAMEDATAH_ENNUL 2
#define FRAMEDATAH_EREAD 4
#define FRAMEDATAH_EOPEN 8
#define FRAMEDATAH_ENOSS 16
#define FRAMEDATAH_EDMRO 32
#define FRAMEDATAH_ELOCK 64
#define FRAMEDATAH_ELOHI 128
#define FRAMEDATAH_ESIZE 256
#define FRAMEDATAH_ESSSZ 512

#define FRAMEDATAH_MSGENULL "Null pointer"
#define FRAMEDATAH_MSGENNUL "Non-null pointer"
#define FRAMEDATAH_MSGEREAD "Error reading frame directory"
#define FRAMEDATAH_MSGEOPEN "Error opening frame file"
#define FRAMEDATAH_MSGENOSS "No sweptsine calibration data in frame"
#define FRAMEDATAH_MSGEDMRO "No IFO-DMRO data in frame"
#define FRAMEDATAH_MSGELOCK "No IFO-Lock data in frame"
#define FRAMEDATAH_MSGELOHI "No locklo/lockhi data in frame"
#define FRAMEDATAH_MSGESIZE "Invalid vector length"
#define FRAMEDATAH_MSGESSSZ "Bad sweptsine calibration data"
/*@}*/

/**
\code
typedef struct
tagFrameData
{
  INT4 inLock;
  INT4 endOfData;
  INT4 newLock;
  INT4 newCalibration;
  /\* ... private data ... *\/
}
FrameData;
\endcode

This is the frame data parameter structure: think of it as something like a
file pointer to the frame data.  The "public" fields are:

<dl>
<dt>inLock</dt><dd> Boolean that user should set to non-zero if data that
  is not "in lock" according to the IFO_Lock channel is desired.</dd>
<dt>endOfData</dt><dd> Boolean that is non-zero if there is no more data.</dd>
<dt>newLock</dt><dd> Boolean that is non-zero if starting a new locked
  section of data.</dd>
<dt>newCalibration</dt><dd> Boolean that is non-zero if new calibration data
  is available.</dd>
</dl>
*/
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
  /** \cond DONT_DOXYGEN */
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
  /** \endcond */
}
FrameData;









void
LALInitializeFrameData (
    LALStatus     *status,
    FrameData **frameData,
    CHAR       *framePath
    );

void
LALFinalizeFrameData (
    LALStatus     *status,
    FrameData **frameData
    );

void
LALGetFrameData (
    LALStatus         *status,
    INT2TimeSeries *data,
    FrameData      *frameData
    );

void
LALGetFrameDataResponse (
    LALStatus                  *status,
    COMPLEX8FrequencySeries *response,
    FrameData               *frameData
    );

/*@}*/

#ifdef  __cplusplus
}
#endif

#endif /* _FRAMEDATA_H */
