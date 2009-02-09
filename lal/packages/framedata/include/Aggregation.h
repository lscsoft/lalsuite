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
#define ONLINE_SAMPLE_RATE 16384

/* data quality */
#define DQ_SCIENCE    (1 << 0) /* SV_SCIENCE & LIGHT */
#define DQ_INJECTION  (1 << 1) /* Injection: same as statevector */
#define DQ_UP         (1 << 2) /* SV_UP & LIGHT */
#define DQ_CALIBRATED (1 << 3) /* SV_UP & LIGHT & (not TRANSIENT) */
#define DQ_BADGAMMA   (1 << 4) /* Bad calibration (outside 0.8 < gamma < 1.2) */
#define DQ_LIGHT      (1 << 5) /* Light in the arms ok */
#define DQ_MISSING    (1 << 6) /* Indication that data was dropped in DMT */


/* function prototypes */
LIGOTimeGPS *XLALAggregationFrameStart(LIGOTimeGPS *gps);

CHAR *XLALAggregationDirectoryPath(CHAR *ifo,
    LIGOTimeGPS *gps);

CHAR *XLALAggregationFrameType(CHAR *ifo);

CHAR *XLALAggregationFrameFilename(CHAR *ifo,
    LIGOTimeGPS *gps);

CHAR *XLALAggregationFramePathFilename(CHAR *ifo,
    LIGOTimeGPS *gps);

CHAR *XLALAggregationFrameURL(CHAR *ifo,
    LIGOTimeGPS *gps);

FrCache *XLALAggregationFrameCache(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length);

FrStream *XLALAggregationFrameStream(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length);

REAL8TimeSeries *XLALAggregationStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length);

INT4TimeSeries *XLALAggregationDQVector(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length);

INT4TimeSeries *XLALAggregationStateVector(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length);

REAL8TimeSeries *XLALAggregationDQStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length,
    INT4 dq_bitmask);

REAL8TimeSeries *XLALAggregationStateStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length,
    INT4 state_bitmask);

REAL8TimeSeries *XLALAggregationDQStateStrainData(CHAR *ifo,
    LIGOTimeGPS *start,
    INT4 length,
    INT4 dq_bitmask,
    INT4 state_bitmask);

#ifdef __cplusplus
}
#endif

#endif
