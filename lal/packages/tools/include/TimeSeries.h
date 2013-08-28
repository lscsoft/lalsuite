/*
 *
 * Copyright (C) 2007  Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


#ifndef _TIMESERIES_H
#define _TIMESERIES_H


#include <stddef.h>
#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup TimeSeriesManipulation
 * \author Kipp Cannon <kipp@gravity.phys.uwm.edu>
 *
 * \brief This is a suite of functions for creating, destroying, and manipulating LAL
 * time series.  One pair of functions (the XLAL version and its LAL
 * counterpart) is available for each method and series type.  For example
 * <tt>XLALCreateREAL4TimeSeries()</tt> is available for creating time series
 * of \c REAL4 data, and the LAL-stype wrapper
 * <tt>LALCreateREAL4TimeSeries()</tt> is provided which is equivalent to the
 * XLAL version in all respects except that it adheres to the LAL calling
 * conventions (eg.\ it takes a \c LALStatus pointer as its first
 * argument, its return type is \c void, etc.).
 *
 */
/*@{*/

/**
 * \name Creation Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/TimeSeries.h>
 *
 * XLALCreate<timeseriestype>()
 * LALCreate<timeseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions create LAL frequency series.  An XLAL function returns a
 * pointer to the newly created series or \c NULL on failure.  The LAL
 * counterpart accepts the address of a pointer which it fills with the
 * address of the newly created series or \c NULL on failure.
 * Additionally, the LAL wrapper provides standard LAL-style error checking
 * via a \c LALStatus pointer.
 */
/*@{*/
COMPLEX8TimeSeries *XLALCreateCOMPLEX8TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
COMPLEX16TimeSeries *XLALCreateCOMPLEX16TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
REAL4TimeSeries *XLALCreateREAL4TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
REAL8TimeSeries *XLALCreateREAL8TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
INT2TimeSeries *XLALCreateINT2TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
INT4TimeSeries *XLALCreateINT4TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
INT8TimeSeries *XLALCreateINT8TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
UINT2TimeSeries *XLALCreateUINT2TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
UINT4TimeSeries *XLALCreateUINT4TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
UINT8TimeSeries *XLALCreateUINT8TimeSeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaT, const LALUnit *sampleUnits, size_t length );
/*@}*/

/**
 * \name Destruction Functions
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/TimeSeries.h>
 *
 * XLALDestroy<timeseriestype>()
 * LALDestroy<timeseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions free all memory associated with a LAL time series.  It is
 * safe to pass \c NULL to these functions.
 *
 */
/*@{*/
void XLALDestroyCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *series );
void XLALDestroyCOMPLEX16TimeSeries ( COMPLEX16TimeSeries *series );
void XLALDestroyREAL4TimeSeries ( REAL4TimeSeries *series );
void XLALDestroyREAL8TimeSeries ( REAL8TimeSeries *series );
void XLALDestroyINT2TimeSeries ( INT2TimeSeries *series );
void XLALDestroyINT4TimeSeries ( INT4TimeSeries *series );
void XLALDestroyINT8TimeSeries ( INT8TimeSeries *series );
void XLALDestroyUINT2TimeSeries ( UINT2TimeSeries *series );
void XLALDestroyUINT4TimeSeries ( UINT4TimeSeries *series );
void XLALDestroyUINT8TimeSeries ( UINT8TimeSeries *series );
/*@}*/


/**
 * \name Cutting Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/TimeSeries.h>
 *
 * XLALCut<timeseriestype>()
 * LALCut<timeseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions create a new time series by extracting a section of an
 * existing time series.
 */
/*@{*/
COMPLEX8TimeSeries *XLALCutCOMPLEX8TimeSeries ( const COMPLEX8TimeSeries *series, size_t first, size_t length );
COMPLEX16TimeSeries *XLALCutCOMPLEX16TimeSeries ( const COMPLEX16TimeSeries *series, size_t first, size_t length );
REAL4TimeSeries *XLALCutREAL4TimeSeries ( const REAL4TimeSeries *series, size_t first, size_t length );
REAL8TimeSeries *XLALCutREAL8TimeSeries ( const REAL8TimeSeries *series, size_t first, size_t length );
INT2TimeSeries *XLALCutINT2TimeSeries ( const INT2TimeSeries *series, size_t first, size_t length );
INT4TimeSeries *XLALCutINT4TimeSeries ( const INT4TimeSeries *series, size_t first, size_t length );
INT8TimeSeries *XLALCutINT8TimeSeries ( const INT8TimeSeries *series, size_t first, size_t length );
UINT2TimeSeries *XLALCutUINT2TimeSeries ( const UINT2TimeSeries *series, size_t first, size_t length );
UINT4TimeSeries *XLALCutUINT4TimeSeries ( const UINT4TimeSeries *series, size_t first, size_t length );
UINT8TimeSeries *XLALCutUINT8TimeSeries ( const UINT8TimeSeries *series, size_t first, size_t length );
/*@}*/

/**
 * \name Resizing Functions
 *
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/TimeSeries.h>
 *
 * XLALResize<timeseriestype>()
 * LALResize<timeseriestype>()
 * XLALShrink<timeseriestype>()
 * LALShrink<timeseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions resize an existing time series.  The new time series will
 * have the given length, and its contents will consist of that part of the
 * original time series that started at sample first.  If first is negative,
 * then the new time series is padded at the start by that many samples.  The
 * time series' epoch is adjusted appropriately.
 */
/*@{*/
COMPLEX8TimeSeries *XLALResizeCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *series, int first, size_t length );
COMPLEX16TimeSeries *XLALResizeCOMPLEX16TimeSeries ( COMPLEX16TimeSeries *series, int first, size_t length );
REAL4TimeSeries *XLALResizeREAL4TimeSeries ( REAL4TimeSeries *series, int first, size_t length );
REAL8TimeSeries *XLALResizeREAL8TimeSeries ( REAL8TimeSeries *series, int first, size_t length );
INT2TimeSeries *XLALResizeINT2TimeSeries ( INT2TimeSeries *series, int first, size_t length );
INT4TimeSeries *XLALResizeINT4TimeSeries ( INT4TimeSeries *series, int first, size_t length );
INT8TimeSeries *XLALResizeINT8TimeSeries ( INT8TimeSeries *series, int first, size_t length );
UINT2TimeSeries *XLALResizeUINT2TimeSeries ( UINT2TimeSeries *series, int first, size_t length );
UINT4TimeSeries *XLALResizeUINT4TimeSeries ( UINT4TimeSeries *series, int first, size_t length );
UINT8TimeSeries *XLALResizeUINT8TimeSeries ( UINT8TimeSeries *series, int first, size_t length );

COMPLEX8TimeSeries *XLALShrinkCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *series, size_t first, size_t length );
COMPLEX16TimeSeries *XLALShrinkCOMPLEX16TimeSeries ( COMPLEX16TimeSeries *series, size_t first, size_t length );
REAL4TimeSeries *XLALShrinkREAL4TimeSeries ( REAL4TimeSeries *series, size_t first, size_t length );
REAL8TimeSeries *XLALShrinkREAL8TimeSeries ( REAL8TimeSeries *series, size_t first, size_t length );
INT2TimeSeries *XLALShrinkINT2TimeSeries ( INT2TimeSeries *series, size_t first, size_t length );
INT4TimeSeries *XLALShrinkINT4TimeSeries ( INT4TimeSeries *series, size_t first, size_t length );
INT8TimeSeries *XLALShrinkINT8TimeSeries ( INT8TimeSeries *series, size_t first, size_t length );
UINT2TimeSeries *XLALShrinkUINT2TimeSeries ( UINT2TimeSeries *series, size_t first, size_t length );
UINT4TimeSeries *XLALShrinkUINT4TimeSeries ( UINT4TimeSeries *series, size_t first, size_t length );
UINT8TimeSeries *XLALShrinkUINT8TimeSeries ( UINT8TimeSeries *series, size_t first, size_t length );
/*@}*/

/**
 * \name Addition Functions
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/TimeSeries.h>
 *
 * XLALAdd<timeseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions add the second argument to the first argument, returning a
 * pointer to the first argument on success or NULL on failure.  The two
 * series must have the same heterodyne frequency and time resolution, and
 * have units that differ only by a dimensionless factor.
 *
 */
/*@{*/
COMPLEX8TimeSeries *XLALAddCOMPLEX8TimeSeries ( COMPLEX8TimeSeries *arg1, const COMPLEX8TimeSeries *arg2 );
COMPLEX16TimeSeries *XLALAddCOMPLEX16TimeSeries ( COMPLEX16TimeSeries *arg1, const COMPLEX16TimeSeries *arg2 );
REAL4TimeSeries *XLALAddREAL4TimeSeries ( REAL4TimeSeries *arg1, const REAL4TimeSeries *arg2 );
REAL8TimeSeries *XLALAddREAL8TimeSeries ( REAL8TimeSeries *arg1, const REAL8TimeSeries *arg2 );
INT2TimeSeries *XLALAddINT2TimeSeries ( INT2TimeSeries *arg1, const INT2TimeSeries *arg2 );
INT4TimeSeries *XLALAddINT4TimeSeries ( INT4TimeSeries *arg1, const INT4TimeSeries *arg2 );
INT8TimeSeries *XLALAddINT8TimeSeries ( INT8TimeSeries *arg1, const INT8TimeSeries *arg2 );
UINT2TimeSeries *XLALAddUINT2TimeSeries ( UINT2TimeSeries *arg1, const UINT2TimeSeries *arg2 );
UINT4TimeSeries *XLALAddUINT4TimeSeries ( UINT4TimeSeries *arg1, const UINT4TimeSeries *arg2 );
UINT8TimeSeries *XLALAddUINT8TimeSeries ( UINT8TimeSeries *arg1, const UINT8TimeSeries *arg2 );
/*@}*/

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif  /* _TIMESERIES_H */
