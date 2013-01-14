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
#include <lal/LALCalibration.h>
#include <lal/FrameCache.h>
#include <lal/LALFrameL.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef SWIG // SWIG interface directives
SWIGLAL(EXTERNAL_STRUCT(FrameH, XLALFrameFree));
#endif

struct FrFile * XLALFrOpenURL( const char *url );
int XLALFrFileCheckSum( FrFile *iFile );
FrHistory * XLALFrHistoryAdd( FrameH *frame, const char *name, const char *comment );
FrDetector * XLALFrDetectorNew( int detector );
void XLALFrameFree( FrameH *frame );
FrameH * XLALFrameNew( LIGOTimeGPS *epoch, double duration,
    const char *project, int run, int frnum, int detectorFlags );
FrVect * XLALFrVectINT4TimeSeries( INT4TimeSeries *series );
FrVect * XLALFrVectREAL4TimeSeries( REAL4TimeSeries *series );
FrVect * XLALFrVectREAL8TimeSeries( REAL8TimeSeries *series );
FrVect * XLALFrVectCOMPLEX8TimeSeries( COMPLEX8TimeSeries *series );
FrVect * XLALFrVectCOMPLEX16TimeSeries( COMPLEX16TimeSeries *series );
FrVect * XLALFrVectREAL4FrequencySeries( REAL4FrequencySeries *series );
FrVect * XLALFrVectREAL8FrequencySeries( REAL8FrequencySeries *series );
FrVect * XLALFrVectCOMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series );
FrVect * XLALFrVectCOMPLEX16FrequencySeries( COMPLEX16FrequencySeries *series );

int XLALFrameAddCalRef( FrameH *frame, COMPLEX8FrequencySeries *series, int version, double duration );
int XLALFrameAddREAL8TimeSeriesProcData( FrameH *frame, REAL8TimeSeries *series );
int XLALFrameAddREAL4TimeSeriesProcData( FrameH *frame, REAL4TimeSeries *series );
int XLALFrameAddINT4TimeSeriesProcData( FrameH *frame, INT4TimeSeries *series );
int XLALFrameAddREAL4TimeSeriesSimData( FrameH *frame, REAL4TimeSeries *series );
int XLALFrameAddREAL8TimeSeriesSimData( FrameH *frame, REAL8TimeSeries *series );
int XLALFrameAddREAL4TimeSeriesAdcData( FrameH *frame, REAL4TimeSeries *series );

COMPLEX8FrequencySeries * XLALFrameGetCalRef( LIGOTimeGPS *validUntil, LIGOTimeGPS *epoch, const char *channel, FrameH *frame );

/* int XLALFrameAddCalFac( FrameH *frame, REAL4TimeSeries *series ); */
int XLALFrameAddCalFac( FrameH *frame, REAL4TimeSeries *series, int version );

/* REAL4TimeSeries * XLALFrameGetCalFac( const char *channel, FrameH *frame ); */
REAL4TimeSeries * XLALFrameGetCalFac( LIGOTimeGPS *epoch, const char *channel, FrameH *frame );

/* high-level function */
LALCalData * XLALFrameGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, FrameH *frame );
LALCalData * XLALFrGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, const char *fname );
LALCalData * XLALFrCacheGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, FrCache *cache );

/* frame writing function */
int XLALFrameWrite(FrameH *frame, const char *fname, int compressLevel);

#ifdef __cplusplus
}
#endif

#endif
