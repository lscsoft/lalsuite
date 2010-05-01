/*
*  Copyright (C) 2007 Duncan Brown, Jolien Creighton, Lisa M. Goggin
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

#ifndef GETDATA_H
#define GETDATA_H

/*
 *
 * Routines to read frame data.
 * Routines to create simulated data.
 * Routines to condition (resample; highpass filter) data.
 * Routines to convert double-precision data to single precision.
 *
 */

#include <lal/LALDatatypes.h>
#include <lal/FrameStream.h>

/* create simulated data */
REAL4TimeSeries * get_simulated_data(
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData,
    REAL8        sampleRate,
    UINT4        simSeed,
    REAL4        simScale
    );

REAL4TimeSeries * get_zero_data(
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData,
    REAL8        sampleRate
    );


/* read frame data */
REAL4TimeSeries * get_frame_data(
    const char  *cacheName,
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData
    );

/* read double-precision frame data and convert to single-precision data */
REAL4TimeSeries * get_frame_data_dbl_convert(
    const char  *cacheName,
    const char  *channelName,
    LIGOTimeGPS *epoch,
    REAL8        duration,
    int          strainData,
    REAL8        dblHighPassFreq
    );

/* low-level routine to read single-precision frame data */
REAL4TimeSeries * fr_get_REAL4TimeSeries( const char *channelName,
    LIGOTimeGPS *epoch, REAL8 duration, FrStream *stream );

/* low-level routine to read double-precision frame data */
REAL8TimeSeries * fr_get_REAL8TimeSeries( const char *channelName,
    LIGOTimeGPS *epoch, REAL8 duration, FrStream *stream );

/* resample time series */
int resample_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 sampleRate );

/* highpass filter time series */
int highpass_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 frequency );

/* highpass filter double-precision time series */
int highpass_REAL8TimeSeries( REAL8TimeSeries *series, REAL8 frequency );

/* trim padding from data */
int trimpad_REAL4TimeSeries( REAL4TimeSeries *series, REAL8 padData );

#endif /* GETDATA_H */
