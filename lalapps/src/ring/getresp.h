#ifndef GETRESP_H
#define GETRESP_H

/*
 *
 * Routines to construct a response function from information contained in
 * a cache file of calibration frames.
 *
 */

#include <lal/LALDatatypes.h>
#include <lal/FrameCache.h>


/* routine read a response function or construct an impulse response function */
COMPLEX8FrequencySeries * get_response(
    const char  *cacheName,
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale,
    int          impulseResponse
    );


/* routine to construct an impulse (uniform in frequency) response function */
COMPLEX8FrequencySeries * get_impulse_response(
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale
    );


/* routine to construct a response function from calibration frame files */
COMPLEX8FrequencySeries * get_frame_response(
    const char  *cacheName,
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale
    );


/* routine to read the reference response function from a frame file cache */
COMPLEX8FrequencySeries * get_reference_response_function( FrCache *calCache,
    const char *ifoName );


/* routine to read the reference sensing function from a frame file cache */
COMPLEX8FrequencySeries * get_reference_sensing_function( FrCache *calCache,
    const char *ifoName );


/* routine to read one cavity gain factor from a frame file cache */
COMPLEX8TimeSeries * get_cavity_gain_factor( FrCache *calCache,
    LIGOTimeGPS *epoch, const char *ifoName );


/* routine to read one open loop gain factor from a frame file cache */
COMPLEX8TimeSeries * get_open_loop_gain_factor( FrCache *calCache,
    LIGOTimeGPS *epoch, const char *ifoName );

#endif /* GETRESP_H */
