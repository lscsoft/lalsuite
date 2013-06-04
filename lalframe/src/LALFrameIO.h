/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Robert Adam Mercer, Xavier Siemens
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

#ifndef _LALFRAMEIO_H
#define _LALFRAMEIO_H

#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>

#if 1
/* to be deprecated */
#include <lal/LALFrameL.h>
#else
/* new code */
#include <lal/LALFrameU.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* alias some types */
#if 1
/* to be deprecated */
typedef struct FrameH LALFrameH;
typedef struct FrFile LALFrFile;
#else
/* new code */
struct tagLALFrFile;
typedef LALFrameUFrameH LALFrameH;
typedef struct tagLALFrFile LALFrFile;
#endif


#ifdef SWIG // SWIG interface directives
SWIGLAL(EXTERNAL_STRUCT(LALFrameH, XLALFrameFree));
#endif

LALFrFile * XLALFrFileOpenURL( const char *url );
int XLALFrFileCksumValid( const LALFrFile *frfile );

int XLALFrameAddFrHistory( LALFrameH *frame, const char *name, const char *comment );
int XLALFrameAddFrDetector( LALFrameH *frame, const LALDetector *detector );
void XLALFrameFree( LALFrameH *frame );
LALFrameH * XLALFrameNew( LIGOTimeGPS *epoch, double duration,
    const char *project, int run, int frnum, int detectorFlags );

int XLALFrameAddREAL4TimeSeriesAdcData( LALFrameH *frame, REAL4TimeSeries *series );

int XLALFrameAddINT4TimeSeriesProcData( LALFrameH *frame, INT4TimeSeries *series );
int XLALFrameAddREAL4TimeSeriesProcData( LALFrameH *frame, REAL4TimeSeries *series );
int XLALFrameAddREAL8TimeSeriesProcData( LALFrameH *frame, REAL8TimeSeries *series );

int XLALFrameAddREAL4TimeSeriesSimData( LALFrameH *frame, REAL4TimeSeries *series );
int XLALFrameAddREAL8TimeSeriesSimData( LALFrameH *frame, REAL8TimeSeries *series );

/* frame writing function */
int XLALFrameWrite(LALFrameH *frame, const char *fname, int compressLevel);

#ifdef __cplusplus
}
#endif

#endif
