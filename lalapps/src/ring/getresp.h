#ifndef GETRESP_H
#define GETRESP_H

#include <lal/LALDatatypes.h>
#include <lal/FrameCache.h>

COMPLEX8FrequencySeries * get_response(
    const char  *cacheName,
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale,
    int          impulseResponse
    );

COMPLEX8FrequencySeries * get_impulse_response(
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale
    );

COMPLEX8FrequencySeries * get_frame_response(
    const char  *cacheName,
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale
    );

COMPLEX8FrequencySeries * get_reference_response_function( FrCache *calCache,
    const char *ifoName );

COMPLEX8FrequencySeries * get_reference_sensing_function( FrCache *calCache,
    const char *ifoName );

COMPLEX8TimeSeries * get_cavity_gain_factor( FrCache *calCache,
    LIGOTimeGPS *epoch, const char *ifoName );

COMPLEX8TimeSeries * get_open_loop_gain_factor( FrCache *calCache,
    LIGOTimeGPS *epoch, const char *ifoName );

#endif /* GETRESP_H */
