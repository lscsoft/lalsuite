#ifndef GETDATA_H
#define GETDATA_H

#include <lal/LALDatatypes.h>
#include <lal/FrameStream.h>

REAL4TimeSeries * get_simulated_data(
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData,
    REAL8        sampleRate,
    UINT4        simSeed,
    REAL4        simScale
    );

REAL4TimeSeries * get_frame_data(
    const char  *cacheName,
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData
    );

REAL4TimeSeries * get_frame_data_dbl_convert(
    const char  *cacheName,
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData,
    REAL8        dblHighPassFreq,
    REAL8        dblScale
    );

REAL4TimeSeries * fr_get_REAL4TimeSeries( const char *channelName,
    LIGOTimeGPS *epoch, REAL8 duration, FrStream *stream );

REAL8TimeSeries * fr_get_REAL8TimeSeries( const char *channelName,
    LIGOTimeGPS *epoch, REAL8 duration, FrStream *stream );

int resample_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 sampleRate );
int highpass_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 frequency );
int highpass_REAL8TimeSeries( REAL8TimeSeries *series, REAL8 frequency );

#endif /* GETDATA_H */
