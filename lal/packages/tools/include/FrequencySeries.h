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


#ifndef _FREQUENCYSERIES_H
#define _FREQUENCYSERIES_H


#include <stddef.h>
#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/**
 * \addtogroup FrequencySeriesManipulation
 * \author Kipp Cannon <kipp@gravity.phys.uwm.edu>
 *
 * \brief This is a suite of functions for creating, destroying, and manipulating LAL
 * frequency series.  One pair of functions (the XLAL version and its LAL
 * counterpart) is available for each method and frequency series type.  For
 * example XLALCreateREAL4FrequencySeries() is available for creating
 * frequency series of \c REAL4 data, and the LAL-stype wrapper
 * LALCreateREAL4FrequencySeries() is provided which is equivalent to
 * the XLAL version in all respects except that it adheres to the LAL calling
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
 * #include <lal/FrequencySeries.h>
 *
 * XLALCreate<frequencyseriestype>()
 * LALCreate<frequencyseriestype>()
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
 *
 */
/*@{*/
COMPLEX8FrequencySeries *XLALCreateCOMPLEX8FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
COMPLEX16FrequencySeries *XLALCreateCOMPLEX16FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
REAL4FrequencySeries *XLALCreateREAL4FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
REAL8FrequencySeries *XLALCreateREAL8FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
INT2FrequencySeries *XLALCreateINT2FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
INT4FrequencySeries *XLALCreateINT4FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
INT8FrequencySeries *XLALCreateINT8FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
UINT2FrequencySeries *XLALCreateUINT2FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
UINT4FrequencySeries *XLALCreateUINT4FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
UINT8FrequencySeries *XLALCreateUINT8FrequencySeries ( const CHAR *name, const LIGOTimeGPS *epoch, REAL8 f0, REAL8 deltaF, const LALUnit *sampleUnits, size_t length );
/*@}*/

/**
 * \name Destruction Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/FrequencySeries.h>
 *
 * XLALDestroy<frequencyseriestype>()
 * LALDestroy<frequencyseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions free all memory associated with a LAL frequency series.  It
 * is safe to pass \c NULL to these functions.
 *
 */
/*@{*/
void XLALDestroyCOMPLEX8FrequencySeries ( COMPLEX8FrequencySeries *series );
void XLALDestroyCOMPLEX16FrequencySeries ( COMPLEX16FrequencySeries *series );
void XLALDestroyREAL4FrequencySeries ( REAL4FrequencySeries *series );
void XLALDestroyREAL8FrequencySeries ( REAL8FrequencySeries *series );
void XLALDestroyINT2FrequencySeries ( INT2FrequencySeries *series );
void XLALDestroyINT4FrequencySeries ( INT4FrequencySeries *series );
void XLALDestroyINT8FrequencySeries ( INT8FrequencySeries *series );
void XLALDestroyUINT2FrequencySeries ( UINT2FrequencySeries *series );
void XLALDestroyUINT4FrequencySeries ( UINT4FrequencySeries *series );
void XLALDestroyUINT8FrequencySeries ( UINT8FrequencySeries *series );
/*@}*/


/**
 * \name Cutting Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/FrequencySeries.h>
 *
 * XLALCut<frequencyseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions create a new frequency series by extracting a section of an
 * existing frequency series.
 *
 */
/*@{*/
COMPLEX8FrequencySeries *XLALCutCOMPLEX8FrequencySeries ( const COMPLEX8FrequencySeries *series, size_t first, size_t length );
COMPLEX16FrequencySeries *XLALCutCOMPLEX16FrequencySeries ( const COMPLEX16FrequencySeries *series, size_t first, size_t length );
REAL4FrequencySeries *XLALCutREAL4FrequencySeries ( const REAL4FrequencySeries *series, size_t first, size_t length );
REAL8FrequencySeries *XLALCutREAL8FrequencySeries ( const REAL8FrequencySeries *series, size_t first, size_t length );
INT2FrequencySeries *XLALCutINT2FrequencySeries ( const INT2FrequencySeries *series, size_t first, size_t length );
INT4FrequencySeries *XLALCutINT4FrequencySeries ( const INT4FrequencySeries *series, size_t first, size_t length );
INT8FrequencySeries *XLALCutINT8FrequencySeries ( const INT8FrequencySeries *series, size_t first, size_t length );
UINT2FrequencySeries *XLALCutUINT2FrequencySeries ( const UINT2FrequencySeries *series, size_t first, size_t length );
UINT4FrequencySeries *XLALCutUINT4FrequencySeries ( const UINT4FrequencySeries *series, size_t first, size_t length );
UINT8FrequencySeries *XLALCutUINT8FrequencySeries ( const UINT8FrequencySeries *series, size_t first, size_t length );
/*@}*/


/**
 * \name Resizing Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/FrequencySeries.h>
 *
 * XLALResize<frequencyseriestype>()
 * XLALShrink<frequencyseriestype>()
 * LALShrink<frequencyseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions resize an existing frequency series.  The new frequency
 * series will have the given length, and its contents will consist of that
 * part of the original time series that started at sample \c first.  If
 * \c first is negative, then the new time series is padded at the start
 * by that many samples.  The frequency series' heterodyne frequency,
 * (f_{0}), is adjusted appropriately.
 *
 * The "Shrink" functions accept non-negative values for the parameter
 * \c first, and are retained only for historical purposes.  New code
 * should use the "Resize" variants.
 *
 */
/*@{*/
COMPLEX8FrequencySeries *XLALResizeCOMPLEX8FrequencySeries ( COMPLEX8FrequencySeries *series, int first, size_t length );
COMPLEX16FrequencySeries *XLALResizeCOMPLEX16FrequencySeries ( COMPLEX16FrequencySeries *series, int first, size_t length );
REAL4FrequencySeries *XLALResizeREAL4FrequencySeries ( REAL4FrequencySeries *series, int first, size_t length );
REAL8FrequencySeries *XLALResizeREAL8FrequencySeries ( REAL8FrequencySeries *series, int first, size_t length );
INT2FrequencySeries *XLALResizeINT2FrequencySeries ( INT2FrequencySeries *series, int first, size_t length );
INT4FrequencySeries *XLALResizeINT4FrequencySeries ( INT4FrequencySeries *series, int first, size_t length );
INT8FrequencySeries *XLALResizeINT8FrequencySeries ( INT8FrequencySeries *series, int first, size_t length );
UINT2FrequencySeries *XLALResizeUINT2FrequencySeries ( UINT2FrequencySeries *series, int first, size_t length );
UINT4FrequencySeries *XLALResizeUINT4FrequencySeries ( UINT4FrequencySeries *series, int first, size_t length );
UINT8FrequencySeries *XLALResizeUINT8FrequencySeries ( UINT8FrequencySeries *series, int first, size_t length );

COMPLEX8FrequencySeries *XLALShrinkCOMPLEX8FrequencySeries ( COMPLEX8FrequencySeries *series, size_t first, size_t length );
COMPLEX16FrequencySeries *XLALShrinkCOMPLEX16FrequencySeries ( COMPLEX16FrequencySeries *series, size_t first, size_t length );
REAL4FrequencySeries *XLALShrinkREAL4FrequencySeries ( REAL4FrequencySeries *series, size_t first, size_t length );
REAL8FrequencySeries *XLALShrinkREAL8FrequencySeries ( REAL8FrequencySeries *series, size_t first, size_t length );
INT2FrequencySeries *XLALShrinkINT2FrequencySeries ( INT2FrequencySeries *series, size_t first, size_t length );
INT4FrequencySeries *XLALShrinkINT4FrequencySeries ( INT4FrequencySeries *series, size_t first, size_t length );
INT8FrequencySeries *XLALShrinkINT8FrequencySeries ( INT8FrequencySeries *series, size_t first, size_t length );
UINT2FrequencySeries *XLALShrinkUINT2FrequencySeries ( UINT2FrequencySeries *series, size_t first, size_t length );
UINT4FrequencySeries *XLALShrinkUINT4FrequencySeries ( UINT4FrequencySeries *series, size_t first, size_t length );
UINT8FrequencySeries *XLALShrinkUINT8FrequencySeries ( UINT8FrequencySeries *series, size_t first, size_t length );
/*@}*/


/**
 * \name Addition Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/FrequencySeries.h>
 *
 * XLALAdd<frequencyseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions add the second argument to the first argument, returning a
 * pointer to the first argument on success or NULL on failure.  The two
 * series must have the same epoch and frequency resolution, and have units
 * that differ only by a dimensionless factor.
 *
 */
/*@{*/
COMPLEX8FrequencySeries *XLALAddCOMPLEX8FrequencySeries ( COMPLEX8FrequencySeries *arg1, const COMPLEX8FrequencySeries *arg2 );
COMPLEX16FrequencySeries *XLALAddCOMPLEX16FrequencySeries ( COMPLEX16FrequencySeries *arg1, const COMPLEX16FrequencySeries *arg2 );
REAL4FrequencySeries *XLALAddREAL4FrequencySeries ( REAL4FrequencySeries *arg1, const REAL4FrequencySeries *arg2 );
REAL8FrequencySeries *XLALAddREAL8FrequencySeries ( REAL8FrequencySeries *arg1, const REAL8FrequencySeries *arg2 );
INT2FrequencySeries *XLALAddINT2FrequencySeries ( INT2FrequencySeries *arg1, const INT2FrequencySeries *arg2 );
INT4FrequencySeries *XLALAddINT4FrequencySeries ( INT4FrequencySeries *arg1, const INT4FrequencySeries *arg2 );
INT8FrequencySeries *XLALAddINT8FrequencySeries ( INT8FrequencySeries *arg1, const INT8FrequencySeries *arg2 );
UINT2FrequencySeries *XLALAddUINT2FrequencySeries ( UINT2FrequencySeries *arg1, const UINT2FrequencySeries *arg2 );
UINT4FrequencySeries *XLALAddUINT4FrequencySeries ( UINT4FrequencySeries *arg1, const UINT4FrequencySeries *arg2 );
UINT8FrequencySeries *XLALAddUINT8FrequencySeries ( UINT8FrequencySeries *arg1, const UINT8FrequencySeries *arg2 );
/*@}*/

/**
 * \name Conjugation Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/FrequencySeries.h>
 *
 * XLALConjugate<frequencyseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions replace a frequency series with its complex conjugate.
 *
 */
/*@{*/
COMPLEX8FrequencySeries *XLALConjugateCOMPLEX8FrequencySeries ( COMPLEX8FrequencySeries *series );
COMPLEX16FrequencySeries *XLALConjugateCOMPLEX16FrequencySeries ( COMPLEX16FrequencySeries *series );
/*@}*/

/**
 * \name Multiplication Functions
 * \heading{Synopsis}
 *
 * \code
 * #include <lal/FrequencySeries.h>
 *
 * XLALMultiply<frequencyseriestype>()
 * \endcode
 *
 * \heading{Description}
 *
 * These functions multiply the first argument by the second argument,
 * returning a pointer to the first argument on success or NULL on failure.
 * The two series must have the same epoch and frequency resolution, and have
 * units that differ only by a dimensionless factor.
 *
 */
/*@{*/
COMPLEX8FrequencySeries *XLALMultiplyCOMPLEX8FrequencySeries ( COMPLEX8FrequencySeries *arg1, const COMPLEX8FrequencySeries *arg2 );
COMPLEX16FrequencySeries *XLALMultiplyCOMPLEX16FrequencySeries ( COMPLEX16FrequencySeries *arg1, const COMPLEX16FrequencySeries *arg2 );
REAL4FrequencySeries *XLALMultiplyREAL4FrequencySeries ( REAL4FrequencySeries *arg1, const REAL4FrequencySeries *arg2 );
REAL8FrequencySeries *XLALMultiplyREAL8FrequencySeries ( REAL8FrequencySeries *arg1, const REAL8FrequencySeries *arg2 );
INT2FrequencySeries *XLALMultiplyINT2FrequencySeries ( INT2FrequencySeries *arg1, const INT2FrequencySeries *arg2 );
INT4FrequencySeries *XLALMultiplyINT4FrequencySeries ( INT4FrequencySeries *arg1, const INT4FrequencySeries *arg2 );
INT8FrequencySeries *XLALMultiplyINT8FrequencySeries ( INT8FrequencySeries *arg1, const INT8FrequencySeries *arg2 );
UINT2FrequencySeries *XLALMultiplyUINT2FrequencySeries ( UINT2FrequencySeries *arg1, const UINT2FrequencySeries *arg2 );
UINT4FrequencySeries *XLALMultiplyUINT4FrequencySeries ( UINT4FrequencySeries *arg1, const UINT4FrequencySeries *arg2 );
UINT8FrequencySeries *XLALMultiplyUINT8FrequencySeries ( UINT8FrequencySeries *arg1, const UINT8FrequencySeries *arg2 );
/*@}*/

/*@}*/

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif  /* _FREQUENCYSERIES_H */
