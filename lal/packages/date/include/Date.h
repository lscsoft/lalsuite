/*
*  Copyright (C) 2012 Karl Wette
*  Copyright (C) 2007 Duncan Brown, David Chin, Jolien Creighton, Kipp Cannon, Reinhard Prix, Stephen Fairhurst
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

#ifndef _DATE_H
#define _DATE_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <lal/LALConstants.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>

#include <lal/DetectorSite.h>

#ifdef  __cplusplus
extern "C"
{
#endif

#define XLAL_BILLION_INT4 1000000000
#define XLAL_BILLION_INT8 LAL_INT8_C( 1000000000 )
#define XLAL_BILLION_REAL8 1e9

/**
 * \addtogroup Date_h
 * \author D.W. Chin, J.D.E. Creighton and Kipp Cannon
 * \brief Provides routines for manipulating date and time information.

\heading{Synopsis}
\code
#include <lal/Date.h>
\endcode

This header covers routines for manipulating date and time
information.  The various time systems are discussed in [\ref esaa1992].

*/
/*@{*/

/** The UNIX time of the GPS origin epoch.
 *
 * 1980 6 JAN 0h UTC is 3657 days after 1970 1 JAN 0h UTC:
 * 8 standard years of 365 days = 2920 days
 * 2 leap years of 366 days = 734 days
 * 5 extra days.
 * Hence 3657*86400=315964800.
 *
 * Note that this deliberately does not account for any leap seconds in the
 * interval.  POSIX:2001 defines the relation between the UNIX time
 * \c time_t \c t and a broken down time \c struct \c tm \c utc as
 * \code
 * t = utc.tm_sec + utc.tm_min*60 + utc.tm_hour*3600
 *     + utc.tm_yday*86400 + (utc.tm_year-70)*31536000
 *     + ((utc.tm_year-69)/4)*86400 - ((utc.tm_year-1)/100)*86400
 *     + ((utc.tm_year+299)/400)*86400;
 * \endcode
 * so if you were to set \c utc.tm_sec=utc.tm_min=utc.tm_hour=0,
 * \c utc.tm_yday=5, and \c utc.tm_year=80, then you get
 * \c t=315964800.  That is what this is.
 */
#define XLAL_EPOCH_UNIX_GPS 315964800

#define XLAL_EPOCH_J2000_0_JD 2451545.0         /**< Julian Day of the J2000.0 epoch (2000 JAN 1 12h UTC). */
#define XLAL_EPOCH_J2000_0_TAI_UTC 32           /**< Leap seconds (TAI-UTC) on the J2000.0 epoch (2000 JAN 1 12h UTC). */
#define XLAL_EPOCH_J2000_0_GPS 630763213        /**< GPS seconds of the J2000.0 epoch (2000 JAN 1 12h UTC). */
#define XLAL_EPOCH_GPS_JD 2444244.5             /**< Julian Day of the GPS epoch (1980 JAN 6 0h UTC) */
#define XLAL_EPOCH_GPS_TAI_UTC 19               /**< Leap seconds (TAI-UTC) on the GPS epoch (1980 JAN 6 0h UTC) */
#define XLAL_MJD_REF 2400000.5                  /**< Reference Julian Day for Mean Julian Day. */
#define XLAL_MODIFIED_JULIEN_DAY(utc) (XLALJulianDay(utc)-XLAL_MJD_REF) /**< Modified Julian Day for specified civil time structure. */

/** This structure stores pointers to a ::LALDetector and a
 * ::LIGOTimeGPS. Its sole purpose is to aggregate these
 * structures for passing to functions.
 */
typedef struct
tagLALPlaceAndGPS
{
    LALDetector *p_detector;   /**< pointer to a detector */
    LIGOTimeGPS *p_gps;        /**< Pointer to a GPS time structure */
}
LALPlaceAndGPS;

/*@}*/

/* ---------- Function prototypes : see respective source.c files for doxygen documentation ---------- */

#ifndef SWIG // exclude from SWIG interface

/* Converts GPS time to nano seconds stored as an INT8. */
INT8 XLALGPSToINT8NS( const LIGOTimeGPS *epoch );

/* Converts nano seconds stored as an INT8 to GPS time. */
LIGOTimeGPS * XLALINT8NSToGPS( LIGOTimeGPS *epoch, INT8 ns );

/* Sets GPS time given GPS integer seconds and residual nanoseconds. */
LIGOTimeGPS * XLALGPSSet( LIGOTimeGPS *epoch, INT4 gpssec, INT8 gpsnan );

/* Sets GPS time given GPS seconds as a REAL8. */
LIGOTimeGPS * XLALGPSSetREAL8( LIGOTimeGPS *epoch, REAL8 t );

/* Returns GPS time as a REAL8. */
REAL8 XLALGPSGetREAL8( const LIGOTimeGPS *epoch );

/** Breaks the GPS time into REAL8 integral and fractional parts,
 * each of which has the same sign as the epoch.  Returns the
 * fractional part, and stores the integral part (as a REAL8)
 * in the object pointed to by iptr.  Like the standard C math
 * library function modf(). */
REAL8 XLALGPSModf( REAL8 *iptr, const LIGOTimeGPS *epoch );

/* Adds dt to a GPS time. */
LIGOTimeGPS * XLALGPSAdd( LIGOTimeGPS *epoch, REAL8 dt );

/* Adds two GPS times. */
LIGOTimeGPS * XLALGPSAddGPS( LIGOTimeGPS *epoch, const LIGOTimeGPS *dt );

/* Difference between two GPS times. */
REAL8 XLALGPSDiff( const LIGOTimeGPS *t1, const LIGOTimeGPS *t0 );

  /* Compares two GPS times. */
int XLALGPSCmp( const LIGOTimeGPS *t0, const LIGOTimeGPS *t1 );

/* Multiply a GPS time by a REAL8 */
LIGOTimeGPS *XLALGPSMultiply( LIGOTimeGPS *gps, REAL8 x );

/* Divide a GPS time by a REAL8 */
LIGOTimeGPS *XLALGPSDivide( LIGOTimeGPS *gps, REAL8 x );

/* Parse an ASCII string into a LIGOTimeGPS structure */
int XLALStrToGPS(LIGOTimeGPS *t, const char *nptr, char **endptr);

/* Return a string containing the ASCII base 10 representation of a LIGOTimeGPS. */
char *XLALGPSToStr(char *, const LIGOTimeGPS *t);

#endif // !SWIG

#ifdef SWIG // SWIG interface directives
SWIGLAL(NEW_EMPTY_ARGUMENT(LIGOTimeGPS*, gpstime));
SWIGLAL(RETURN_VALUE(LIGOTimeGPS*, XLALGPSTimeNow));
#endif

/* This function returns the current GPS time according to the system clock */
LIGOTimeGPS* XLALGPSTimeNow( LIGOTimeGPS *gpstime );

#ifdef SWIG // SWIG interface directives
SWIGLAL_CLEAR(NEW_EMPTY_ARGUMENT(LIGOTimeGPS*, gpstime));
#endif

/* Returns the leap seconds TAI-UTC at a given GPS second. */
int XLALLeapSeconds( INT4 gpssec );

/* Returns the leap seconds GPS-UTC at a given GPS second. */
int XLALGPSLeapSeconds( INT4 gpssec );

/* Returns the leap seconds TAI-UTC for a given UTC broken down time. */
int XLALLeapSecondsUTC( const struct tm *utc );

/* Returns the GPS seconds since the GPS epoch for a specified UTC time structure. */
INT4 XLALUTCToGPS( const struct tm *utc );

#ifdef SWIG // SWIG interface directives
SWIGLAL(EMPTY_ARGUMENT(struct tm*, utc));
SWIGLAL(RETURN_VALUE(struct tm*, XLALGPSToUTC));
#endif

/* Returns a pointer to a tm structure representing the time
 * specified in seconds since the GPS epoch.  */
struct tm* XLALGPSToUTC( struct tm *utc, INT4 gpssec );

#ifdef SWIG // SWIG interface directives
SWIGLAL_CLEAR(EMPTY_ARGUMENT(struct tm*, utc));
#endif

/* Returns the Julian Day (JD) corresponding to the date given in a broken
 * down time structure. */
REAL8 XLALJulianDay( const struct tm *utc );

/* Returns the Modified Julian Day (MJD) corresponding to the date given in a broken down time structure.*/
INT4 XLALModifiedJulianDay( const struct tm *utc );

/* Fill in missing fields of a C 'tm' broken-down time struct. */
int XLALFillBrokenDownTime( struct tm *tm );

/* Returns the Greenwich mean or aparent sideral time in radians. */
REAL8 XLALGreenwichSiderealTime(
        const LIGOTimeGPS *gpstime,
        REAL8 equation_of_equinoxes
);

/* Returns the Greenwich Mean Sidereal Time in RADIANS for a specified GPS time. */
REAL8 XLALGreenwichMeanSiderealTime(
        const LIGOTimeGPS *gpstime
);

/* Returns the GPS time for the given Greenwich mean sidereal time (in radians). */
LIGOTimeGPS *XLALGreenwichMeanSiderealTimeToGPS(
        REAL8 gmst,
        LIGOTimeGPS *gps
);

/* Returns the GPS time for the given Greenwich sidereal time (in radians). */
LIGOTimeGPS *XLALGreenwichSiderealTimeToGPS(
        REAL8 gmst,
        REAL8 equation_of_equinoxes,
        LIGOTimeGPS *gps
);

/* Determines if a given time is playground data. */
int XLALINT8NanoSecIsPlayground (
        INT8 ns
);

#ifdef  __cplusplus
}
#endif

#endif /* _DATE_H */
