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

#ifndef _LALSIMINSPIRALSPHHARMSERIES_H
#define _LALSIMINSPIRALSPHHARMSERIES_H

#include <lal/LALDatatypes.h>
#include <lal/TimeSeries.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/* 
 * Structure to carry a collection of spherical harmonic modes in COMPLEX16 
 * time series. Contains convenience getter and setter functions, as well as
 * a convienence "maximum l mode" function. Implemented as a singly forward
 * linked list.
 */
typedef struct tagSphHarmTimeSeries {
    COMPLEX16TimeSeries*            mode; /**< The sequences of sampled data. */
    UINT4                           l; /**< Node mode l  */
    INT4                            m; /**< Node submode m  */
	REAL8Sequence*					tdata; /**< Timestamp values */
    struct tagSphHarmTimeSeries*    next; /**< next pointer */
} SphHarmTimeSeries;


typedef struct tagSphHarmFrequencySeries {
    COMPLEX16FrequencySeries*            mode; /**< The sequences of sampled data. */
    UINT4                           l; /**< Node mode l  */
    INT4                            m; /**< Node submode m  */
	REAL8Sequence*					fdata; /**< Frequency values */
    struct tagSphHarmFrequencySeries*    next; /**< next pointer */
} SphHarmFrequencySeries;


/* 
 * Create a SphHarmTimeSeries. If appended is not NULL, this will prepend a new
 * structure to the list by duplicating the mode inmode, mode numbers l, and m, 
 * and then set the next pointer to the appended structure.
 */
SphHarmTimeSeries* XLALSphHarmTimeSeriesAddMode( 
		SphHarmTimeSeries *appended,  /**< List structure to prepend to */
		const COMPLEX16TimeSeries* inmode,  /**< mode series to contain */
		UINT4 l, /**< major mode number */
		INT4 m  /**< minor mode number */
);

/* 
 * Set the tdata pointer to a REAL8Sequence for all members of the 
 * SphHarmTimeSeries linked list. This is mainly intended for use with
 * unevenly sampled time series data
 */
void XLALSphHarmTimeSeriesSetTData( 
		SphHarmTimeSeries *ts,  /**< List structure to set tdata */
		REAL8Sequence* fdata  /**< sequence of timestamps */
);

/* 
 * Get the tdata pointer.
 */
REAL8Sequence* XLALSphHarmTimeSeriesGetTData( 
		SphHarmTimeSeries *ts  /**< List structure to get tdata */
);

/* 
 * Destroy a SphHarmTimeSeries. Note that this will destroy any 
 * COMPLEX16TimeSeries which it has references to.
 */
void XLALDestroySphHarmTimeSeries( SphHarmTimeSeries* ts );

/* 
 * Destroy a SphHarmTimeSeries. Note that this will destroy any 
 * COMPLEX16TimeSeries which it has references to.
 */
UINT4 XLALSphHarmTimeSeriesGetMaxL( SphHarmTimeSeries* ts );

#ifdef SWIG   // SWIG interface directives
SWIGLAL(RETURNS_PROPERTY(COMPLEX16TimeSeries*, XLALSphHarmTimeSeriesGetMode));
#endif

/* 
 * Get the mode-decomposed time series corresponding to l,m.
 */
COMPLEX16TimeSeries* XLALSphHarmTimeSeriesGetMode( 
				SphHarmTimeSeries *ts, 
				UINT4 l, 
				INT4 m 
);

SphHarmTimeSeries *XLALResizeSphHarmTimeSeries(
        SphHarmTimeSeries *ts,
        int first,
        size_t length
        );

SphHarmFrequencySeries *XLALSphHarmFrequencySeriesFromSphHarmTimeSeries(
        SphHarmTimeSeries *hlms_TD
        );

/* 
 * Create a SphHarmFrequencySeries. If appended is not NULL, this will prepend a new
 * structure to the list by duplicating the mode inmode, mode numbers l, and m, 
 * and then set the next pointer to the appended structure.
 */
SphHarmFrequencySeries* XLALSphHarmFrequencySeriesAddMode( 
		SphHarmFrequencySeries *appended,  /**< List structure to prepend to */
		const COMPLEX16FrequencySeries* inmode,  /**< mode series to contain */
		UINT4 l, /**< major mode number */
		INT4 m  /**< minor mode number */
);

SphHarmTimeSeries *XLALSphHarmTimeSeriesFromSphHarmFrequencySeriesDataAndPSD(
                                                                             SphHarmFrequencySeries *hlms, 
                                                                             COMPLEX16FrequencySeries* data,
                                                                             COMPLEX16FrequencySeries* psd
                                                                             );


/* 
 * Set the tdata pointer to a REAL8Sequence for all members of the 
 * SphHarmFrequencySeries linked list. This is mainly intended for use with
 * unevenly sampled time series data
 */
void XLALSphHarmFrequencySeriesSetFData( 
		SphHarmFrequencySeries *ts,  /**< List structure to add tdata */
		REAL8Sequence* tdata  /**< sequence of timestamps */
);

/* 
 * Get the fdata pointer.
 */
REAL8Sequence* XLALSphHarmFrequencySeriesGetFData( 
		SphHarmFrequencySeries *ts  /**< List structure to get fdata */
);

/* 
 * Destroy a SphHarmFrequencySeries. Note that this will destroy any 
 * COMPLEX16TimeSeries which it has references to.
 */
void XLALDestroySphHarmFrequencySeries( SphHarmFrequencySeries* ts );

/* 
 * Destroy a SphHarmFrequencySeries. Note that this will destroy any 
 * COMPLEX16FrequencySeries which it has references to.
 */
UINT4 XLALSphHarmFrequencySeriesGetMaxL( SphHarmFrequencySeries* ts );

#ifdef SWIG   // SWIG interface directives
SWIGLAL(RETURNS_PROPERTY(COMPLEX16FrequencySeries*, XLALSphHarmFrequencySeriesGetMode));
#endif

/* 
 * Get the mode-decomposed frequency series corresponding to l,m.
 */
COMPLEX16FrequencySeries* XLALSphHarmFrequencySeriesGetMode( 
				SphHarmFrequencySeries *ts, 
				UINT4 l, 
				INT4 m 
);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMINSPIRAL_H */
