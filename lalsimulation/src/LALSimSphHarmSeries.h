/*
 * Copyright (C) 2013 E. Ochsner, C. Pankow
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

#ifndef _LALSIMSPHHARMSERIES_H
#define _LALSIMSPHHARMSERIES_H

#include <lal/LALDatatypes.h>
#include <lal/TimeSeries.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * @defgroup LALSimSphHarmSeries_h Header LALSimSphHarmSeries.h
 * @ingroup lalsimulation_general
 *
 * @brief Routines for manipulating spherical harmonic time- and
 * frequency-series.
 *
 * @{
 */

/**
 * Structure to carry a collection of spherical harmonic modes in COMPLEX16 
 * time series. Contains convenience getter and setter functions, as well as
 * a convienence "maximum l mode" function. Implemented as a singly forward
 * linked list.
 */
typedef struct tagSphHarmTimeSeries {
    COMPLEX16TimeSeries*            mode; /**< The sequences of sampled data. */
    UINT4                           l; /**< Node mode l  */
    INT4                            m; /**< Node submode m  */
    REAL8Sequence*                  tdata; /**< Timestamp values */
    struct tagSphHarmTimeSeries*    next; /**< next pointer */
} SphHarmTimeSeries;


typedef struct tagSphHarmFrequencySeries {
    COMPLEX16FrequencySeries*       mode; /**< The sequences of sampled data. */
    UINT4                           l; /**< Node mode l  */
    INT4                            m; /**< Node submode m  */
    REAL8Sequence*                  fdata; /**< Frequency values */
    struct tagSphHarmFrequencySeries*    next; /**< next pointer */
} SphHarmFrequencySeries;

/** @} */

SphHarmTimeSeries* XLALSphHarmTimeSeriesAddMode(SphHarmTimeSeries *appended, const COMPLEX16TimeSeries* inmode, UINT4 l, INT4 m);
void XLALSphHarmTimeSeriesSetTData(SphHarmTimeSeries *ts, REAL8Sequence* fdata);
REAL8Sequence* XLALSphHarmTimeSeriesGetTData(SphHarmTimeSeries *ts);
void XLALDestroySphHarmTimeSeries(SphHarmTimeSeries* ts);
UINT4 XLALSphHarmTimeSeriesGetMaxL(SphHarmTimeSeries* ts);

#ifdef SWIG
SWIGLAL(RETURN_OWNED_BY_1ST_ARG(COMPLEX16TimeSeries*, XLALSphHarmTimeSeriesGetMode));
#endif

COMPLEX16TimeSeries* XLALSphHarmTimeSeriesGetMode(SphHarmTimeSeries *ts, UINT4 l, INT4 m);
SphHarmTimeSeries *XLALResizeSphHarmTimeSeries(SphHarmTimeSeries *ts, int first, size_t length);

SphHarmFrequencySeries *XLALSphHarmFrequencySeriesFromSphHarmTimeSeries(SphHarmTimeSeries *hlms_TD);
SphHarmFrequencySeries* XLALSphHarmFrequencySeriesAddMode(SphHarmFrequencySeries *appended, const COMPLEX16FrequencySeries* inmode, UINT4 l, INT4 m);
SphHarmTimeSeries *XLALSphHarmTimeSeriesFromSphHarmFrequencySeriesDataAndPSD(SphHarmFrequencySeries *hlms, COMPLEX16FrequencySeries* data, COMPLEX16FrequencySeries* psd);
void XLALSphHarmFrequencySeriesSetFData(SphHarmFrequencySeries *ts, REAL8Sequence* tdata);
REAL8Sequence* XLALSphHarmFrequencySeriesGetFData(SphHarmFrequencySeries *ts);
void XLALDestroySphHarmFrequencySeries(SphHarmFrequencySeries* ts);
UINT4 XLALSphHarmFrequencySeriesGetMaxL(SphHarmFrequencySeries* ts);

#ifdef SWIG
SWIGLAL(RETURN_OWNED_BY_1ST_ARG(COMPLEX16FrequencySeries*, XLALSphHarmFrequencySeriesGetMode));
#endif

COMPLEX16FrequencySeries* XLALSphHarmFrequencySeriesGetMode(SphHarmFrequencySeries *ts, UINT4 l, INT4 m);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif
