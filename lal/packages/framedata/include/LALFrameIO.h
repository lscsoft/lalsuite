#ifndef _LALFRAMEIO_H
#define _LALFRAMEIO_H

#include <FrameL.h>
#include <lal/LALDatatypes.h>
#include <lal/LALCalibration.h>
#include <lal/FrameCache.h>

#ifdef __cplusplus
extern "C" {
#endif

struct FrFile * XLALFrOpenURL( const char *url );
int XLALFrFileCheckSum( FrFile *iFile );
FrHistory * XLALFrHistoryAdd( FrameH *frame, const char *name, const char *comment );
FrDetector * XLALFrDetectorNew( int detector );
FrameH * XLALFrameNew( LIGOTimeGPS *epoch, double duration, 
    const char *project, int run, int frnum, int detectorFlags );
FrVect * XLALFrVectREAL4TimeSeries( REAL4TimeSeries *series );
FrVect * XLALFrVectREAL8TimeSeries( REAL8TimeSeries *series );
FrVect * XLALFrVectCOMPLEX8TimeSeries( COMPLEX8TimeSeries *series );
FrVect * XLALFrVectCOMPLEX16TimeSeries( COMPLEX16TimeSeries *series );
FrVect * XLALFrVectREAL4FrequencySeries( REAL4FrequencySeries *series );
FrVect * XLALFrVectREAL8FrequencySeries( REAL8FrequencySeries *series );
FrVect * XLALFrVectCOMPLEX8FrequencySeries( COMPLEX8FrequencySeries *series );
FrVect * XLALFrVectCOMPLEX16FrequencySeries( COMPLEX16FrequencySeries *series );

int XLALFrameAddCalRef( FrameH *frame, COMPLEX8FrequencySeries *series, int version, double duration );

COMPLEX8FrequencySeries * XLALFrameGetCalRef( LIGOTimeGPS *validUntil, LIGOTimeGPS *epoch, const char *channel, FrameH *frame );

/* int XLALFrameAddCalFac( FrameH *frame, REAL4TimeSeries *series ); */
int XLALFrameAddCalFac( FrameH *frame, REAL4TimeSeries *series, int version );

/* REAL4TimeSeries * XLALFrameGetCalFac( const char *channel, FrameH *frame ); */
REAL4TimeSeries * XLALFrameGetCalFac( LIGOTimeGPS *epoch, const char *channel, FrameH *frame );

/* high-level function */
LALCalData * XLALFrameGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, FrameH *frame );
LALCalData * XLALFrGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, const char *fname );
LALCalData * XLALFrCacheGetCalData( LIGOTimeGPS *epoch, const char *readoutChannel, FrCache *cache );

#ifdef __cplusplus
}
#endif

#endif
