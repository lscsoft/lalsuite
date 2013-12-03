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

#ifndef GETRESP_H
#define GETRESP_H

/*
 *
 * Routines to construct a response function from information contained in
 * a cache file of calibration frames.
 *
 */

#include <lal/LALDatatypes.h>
#include <lal/LALCache.h>


/* routine read a response function or construct an impulse response function */
COMPLEX8FrequencySeries * get_response(
    const char  *cacheName,
    const char  *ifoName,
    LIGOTimeGPS *epoch,
    REAL8        dataDuration,
    REAL8        dataSampleRate,
    REAL4        responseScale,
    int          strainData,
    const char  *channel_name
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
    REAL4        responseScale,
    const char *channel_name
    );


/* routine to read the reference response function from a frame file cache */
COMPLEX8FrequencySeries * get_reference_response_function( LALCache *calCache,
    const char *ifoName );


/* routine to read the reference sensing function from a frame file cache */
COMPLEX8FrequencySeries * get_reference_sensing_function( LALCache *calCache,
    const char *ifoName );


/* routine to read one cavity gain factor from a frame file cache */
COMPLEX8TimeSeries * get_cavity_gain_factor( LALCache *calCache,
    LIGOTimeGPS *epoch, const char *ifoName );


/* routine to read one open loop gain factor from a frame file cache */
COMPLEX8TimeSeries * get_open_loop_gain_factor( LALCache *calCache,
    LIGOTimeGPS *epoch, const char *ifoName );

#endif /* GETRESP_H */
