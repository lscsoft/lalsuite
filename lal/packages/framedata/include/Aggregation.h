/*
 * Aggregation.h - online frame data aggregation rountines
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 * MA  02111-1307  USA
 *
 * Copyright (C) 2009 Adam Mercer
 */

#ifndef _AGGREGATION_H
#define _AGGREGATION_H

#include <lal/LALDatatypes.h>
#include <lal/FrameCache.h>
#include <lal/FrameStream.h>

#ifdef __cplusplus
extern "C" {
#endif

/* frame channels */
#define ONLINE_STRAIN_CHANNEL "DMT-STRAIN"
#define ONLINE_STATE_VECTOR "DMT-STATE_VECTOR"
#define ONLINE_DQ_VECTOR "DMT-DATA_QUALITY_VECTOR"

/* frame metadata */
#define ONLINE_FRAME_DURATION 16

/* data quality */
#define LAL_DQ_SCIENCE    (1 << 0) /* SV_SCIENCE & LIGHT */
#define LAL_DQ_INJECTION  (1 << 1) /* Injection: same as statevector */
#define LAL_DQ_UP         (1 << 2) /* SV_UP & LIGHT */
#define LAL_DQ_CALIBRATED (1 << 3) /* SV_UP & LIGHT & (not TRANSIENT) */
#define LAL_DQ_BADGAMMA   (1 << 4) /* Bad calibration */
#define LAL_DQ_LIGHT      (1 << 5) /* Light in the arms ok */
#define LAL_DQ_MISSING    (1 << 6) /* Data was dropped in DMT */


/* return frame gps start time for given gps time */
LIGOTimeGPS *XLALAggregationFrameStart(LIGOTimeGPS *gps);

/* return frame type for given ifo */
CHAR *XLALAggregationFrameType(CHAR *ifo);

/* return path for given ifo and gps time */
CHAR *XLALAggregationDirectoryPath(CHAR *ifo,
    LIGOTimeGPS *gps);

/* return frame filename for given ifo and gps time */
CHAR *XLALAggregationFrameFilename(CHAR *ifo,
    LIGOTimeGPS *gps);

/* return full path to frame for given ifo and gps time*/
CHAR *XLALAggregationFramePathFilename(CHAR *ifo,
    LIGOTimeGPS *gps);

/* return url to frame for a given gps time and ifo */
CHAR *XLALAggregationFrameURL(CHAR *ifo,
    LIGOTimeGPS *gps);

/* return gps time of latest frame written to disk */
LIGOTimeGPS *XLALAggregationLatestGPS(CHAR *ifo);

/* return frame cache given ifo, gps time, and duration */
FrCache *XLALAggregationFrameCache(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration);

/* return frame stream for given ifo, gps time, and duration */
FrStream *XLALAggregationFrameStream(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration);

/* return strain data time series for given ifo, gps time, and duration */
REAL8TimeSeries *XLALAggregationStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration);

/* return data quality vector time series for given ifo, gps time and
 * duration */
INT4TimeSeries *XLALAggregationDQVector(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration);

/* return state vector time series for given ifo, gps time, and duration */
INT4TimeSeries *XLALAggregationStateVector(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration);

/* return strain data time series for given ifo, gps time, duration, and
 * data quality vector bitmask */
REAL8TimeSeries *XLALAggregationDQStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration,
    INT4 dq_bitmask);

/* return start position of data gap */
UINT4 XLALAggregationDQGapStart(INT4TimeSeries *series,
    INT4 dq_bitmask);

/* return end position of data gap */
UINT4 XLALAggregationDQGapEnd(INT4TimeSeries *series,
    INT4 dq_bitmask);

/* return end position of data gap - deprecated */
UINT4 XLALAggregationDQGap(INT4TimeSeries *series,
    INT4 dq_bitmask);

/* return strain data time series for given ifo, gps time, duration, and
 * a maximum wait time */
REAL8TimeSeries *XLALAggregationStrainDataWait(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration,
    UINT4 max_wait);

/* return data quality vector time series for given ifo, gps time,
 * duration, and a maximum wait time */
INT4TimeSeries *XLALAggregationDQVectorWait(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration,
    UINT4 max_wait);

/* return state vector time series for given ifo, gps time, duration,
 * and a maximum wait time */
INT4TimeSeries *XLALAggregationStateVectorWait(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration,
    UINT4 max_wait);

/* check that all frames files, for requested data segment, are
 * available */
INT4 XLALAggregationStatFiles(CHAR *ifo,
    LIGOTimeGPS *start,
    REAL8 duration);

#ifdef __cplusplus
}
#endif

#endif
