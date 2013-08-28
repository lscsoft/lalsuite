/*
*  Copyright (C) 2012 Karl Wette
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton, Kipp Cannon
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

#include <math.h>
#include <time.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>

#include "XLALLeapSeconds.h" /* contains the leap second table */

/**
 * \defgroup XLALCivilTime_c CivilTime
 * \ingroup Date_h
 * \author Chin, D. W. and Creighton, J. D. E.
 *
 * \brief XLAL routines for converting civil time structures to GPS times.
 *
 * Civil time structures, represented in the C library by the \c struct \c tm
 * structure, known as a broken down time structure, gives the time by normal
 * measures (hours, minutes, seconds since midnight, day of the month, month of
 * the year, and year).  LAL uses seconds since the the GPS epoch as its
 * reference.  The GPS epoch is defined to be 0h UTC 6 Jan 1980 in civil time.
 * The process of converting between a civil time and a GPS time is almost done
 * by standard C library functions, though the C library functions use a
 * different reference epoch, which is called the UNIX epoch.
 *
 * The tricky bits are:
 * - What about time zones?
 * - What about leap seconds?
 *
 * This code does not deal with time zones at all.  All civil time structures
 * are taken to be in Coordinated Universal Time or UTC.  The user must convert
 * local times into UTC before using these routines.
 *
 * Leap seconds are accounted for.  But since the addition and subtraction of
 * leap seconds is not deterministic, a leap second table needs to be
 * maintained to account for the number of leap seconds in effect at a
 * particular time.
 *
 * Leap seconds are defined as the difference between UTC and a (yet another)
 * standard of time called the International Atomic Time or TAI.  UTC is
 * ultimately determined by the position of the stars, so it is occationally
 * altered by a second to keep the location of fixed points on the celestial
 * sphere correct to within plus or minus 0.9 seconds.  TAI, like the GPS time
 * used in LAL, is just the number of seconds since some epoch and is not
 * affected by leap seconds.  The difference between TAI and UTC, TAI-UTC,
 * is the number of leap seconds.
 *
 * Note that UTC is fixed to atomic time, though there are an integer number
 * of leap seconds.  The real civil time, as measured in terms of points on
 * the celestial sphere, is UT1.  UTC and UT1 are kept within 0.9 seconds of
 * each other by the introduction of leap seconds.  But if what you want is
 * UT1 note that UTC can be off by up to about a seconds.  For this reason,
 * we assume that accuracy at the second scale is sufficient, so the routines
 * have only second precision.  If you need more accuracy, you'll need to be
 * monitoring UT1.
 *
 * Another way of representing the civil time in in terms of Julian days.
 * There is a routine for converting a UTC time into Julian days.
 * The inverse conversion is not attempted.
 *
 */
/*@{*/

/* change in TAI-UTC from previous second:
 *
 * return values:
 *   -1: TAI-UTC has decreased by one second (this hasn't happened yet).
 *       In this case, UTC will skip a second going from 23:59:58 at
 *       gpssec-1 to 00:00:00 (of the following day) at gpssec.
 *    0: This is not a leap second: UTC has just advanced one second
 *       in going from gpssec-1 to gpssec.
 *   +1: TAI-UTC has increased by one second (this hasn't happened yet)
 *       In this case, UTC will add a second going from 23:59:59 at
 *       gpssec-1 to 23:59:60 (of the same day) at gpssec.
 */
static int delta_tai_utc( INT4 gpssec )
{
  int leap;

  /* assume calling function has already checked this */
  /*
  if ( gpssec <= leaps[0].gpssec )
  {
    fprintf( stderr, "error: don't know leap seconds before gps time %d\n",
        leaps[0].gpssec );
    abort();
  }
  */

  for ( leap = 1; leap < numleaps; ++leap )
    if ( gpssec == leaps[leap].gpssec )
      return leaps[leap].taiutc - leaps[leap-1].taiutc;

  return 0;
}

/** Returns the leap seconds TAI-UTC at a given GPS second. */
int XLALLeapSeconds( INT4 gpssec /**< [In] Seconds relative to GPS epoch.*/ )
{
  int leap;

  if ( gpssec < leaps[0].gpssec )
  {
    XLALPrintError( "XLAL Error - Don't know leap seconds before GPS time %d\n",
        leaps[0].gpssec );
    XLAL_ERROR( XLAL_EDOM );
  }

  /* scan leap second table and locate the appropriate interval */
  for ( leap = 1; leap < numleaps; ++leap )
    if ( gpssec < leaps[leap].gpssec )
      break;

  return leaps[leap-1].taiutc;
}


/** Returns the leap seconds GPS-UTC at a given GPS second. */
int XLALGPSLeapSeconds( INT4 gpssec /**< [In] Seconds relative to GPS epoch.*/ )
{
  int leapTAI;
  int leapGPS;

  leapTAI = XLALLeapSeconds ( gpssec );
  if ( leapTAI < 0 )
    XLAL_ERROR( XLAL_EFUNC );

  leapGPS = leapTAI - 19;	/* subtract 19 seconds to get leap-seconds wrt to GPS epoch */

  return leapGPS;
}


/** Returns the leap seconds TAI-UTC for a given UTC broken down time. */
int XLALLeapSecondsUTC( const struct tm *utc /**< [In] UTC as a broken down time.*/ )
{
  REAL8 jd;
  int leap;

  jd = XLALJulianDay( utc );
  if ( XLAL_IS_REAL8_FAIL_NAN( jd ) )
    XLAL_ERROR( XLAL_EFUNC );

  if ( jd < leaps[0].jd )
  {
    XLALPrintError( "XLAL Error - Don't know leap seconds before Julian Day %9.1f\n", leaps[0].jd );
    XLAL_ERROR( XLAL_EDOM );
  }

  /* scan leap second table and locate the appropriate interval */
  for ( leap = 1; leap < numleaps; ++leap )
    if ( jd < leaps[leap].jd )
      break;

  return leaps[leap-1].taiutc;
}


/**
 * Returns the GPS seconds since the GPS epoch for a
 * specified UTC time structure.
 */
INT4 XLALUTCToGPS( const struct tm *utc /**< [In] UTC time in a broken down time structure. */ )
{
  time_t unixsec;
  INT4 gpssec;
  int leapsec;

  /* compute leap seconds */
  leapsec = XLALLeapSecondsUTC( utc );
  if ( leapsec < 0 )
    XLAL_ERROR( XLAL_EFUNC );
  /* compute unix epoch time: seconds since 1970 JAN 1 0h UTC */
  /* POSIX:2001 definition of seconds since the (UNIX) Epoch */
  unixsec = utc->tm_sec + utc->tm_min*60 + utc->tm_hour*3600
    + utc->tm_yday*86400 + (utc->tm_year-70)*31536000
    + ((utc->tm_year-69)/4)*86400 - ((utc->tm_year-1)/100)*86400
    + ((utc->tm_year+299)/400)*86400;
  gpssec  = unixsec;
  gpssec -= XLAL_EPOCH_UNIX_GPS; /* change to gps epoch */
  gpssec += leapsec - XLAL_EPOCH_GPS_TAI_UTC;
  /* now check to see if this is an additional leap second */
  if ( utc->tm_sec == 60 )
    --gpssec;
  return gpssec;
}


/**
 * Returns a pointer to a tm structure representing the time
 * specified in seconds since the GPS epoch.
 */
struct tm * XLALGPSToUTC(
    struct tm *utc, /**< [Out] Pointer to tm struct where result is stored. */
    INT4 gpssec /**< [In] Seconds since the GPS epoch. */
    )
{
  time_t unixsec;
  int leapsec;
  int delta;
  leapsec = XLALLeapSeconds( gpssec );
  if ( leapsec < 0 )
    XLAL_ERROR_NULL( XLAL_EFUNC );
  unixsec  = gpssec - leapsec + XLAL_EPOCH_GPS_TAI_UTC; /* get rid of leap seconds */
  unixsec += XLAL_EPOCH_UNIX_GPS; /* change to unix epoch */
  memset( utc, 0, sizeof( *utc ) ); /* blank out utc structure */
  utc = gmtime_r( &unixsec, utc ); /* FIXME: use gmtime ?? */
  /* now check to see if we need to add a 60th second to UTC */
  if ( ( delta = delta_tai_utc( gpssec ) ) > 0 )
    utc->tm_sec += 1; /* delta only ever is one, right?? */
  return utc;
}


/**
 * Returns the Julian Day (JD) corresponding to the date given in a broken
 * down time structure.
 *
 * See \ref esaa1992 and \ref green1985 for details.  First, some
 * definitions:
 *
 * Mean Julian Year = 365.25 days
 * Julian Epoch = 1 Jan 4713BCE, 12:00 GMT (4713 BC Jan 01d.5 GMT)
 * Fundamental Epoch J2000.0 = 2001-01-01.5 TDB
 *
 * Julian Date is the amount of time elapsed since the Julian Epoch,
 * measured in days and fractions of a day.  There are a couple of
 * complications arising from the length of a year:  the Tropical Year is
 * 365.2422 days.  First, the Gregorian correction where 10 days
 * (1582-10-05 through 1582-10-14) were eliminated.  Second, leap years:
 * years ending with two zeroes (e.g., 1700, 1800) are leap only if
 * divisible by 400;  so, 400 civil years contain 400 * 365.25 - 3 = 146097
 * days.  So, the Julian Date of J2000.0 is JD 2451545.0, and thus the
 * Julian Epoch = J2000.0 + (JD - 2451545) / 365.25, i.e., number of years
 * elapsed since J2000.0.
 *
 * One algorithm for computing the Julian Day is from \ref vfp1979 based
 * on a formula in \ref esaa1992 where the algorithm is due to
 * \ref fvf1968 and ``compactified'' by P. M. Muller and R. N. Wimberly.
 * The formula is
 *
 * \f[
 * jd = 367 \times y - 7 \times (y + (m + 9)/12)/4 - 3 \times ((y + (m -
 * 9)/7)/100 + 1)/4 + 275 \times m/9 + d + 1721029
 * \f]
 *
 * where jd is the Julian day number, y is the year, m is the month (1-12),
 * and d is the day (1-31).  This formula is valid only for JD > 0, i.e.,
 * after -4713 Nov 23 = 4712 BCE Nov 23.
 *
 * A shorter formula from the same reference, but which only works for
 * dates since 1900 March is:
 *
 * \f[
 * jd = 367 \times y - 7 \times (y + (m + 9)/12)/4 + 275 \times m/9 + d +
 * 1721014
 * \f]
 *
 * We will use this shorter formula since there is unlikely to be any
 * analyzable data from before 1900 March.
 *
 */
REAL8 XLALJulianDay( const struct tm *utc /**< [In] UTC time in a broken down time structure. */ )
{
  const int sec_per_day = 60 * 60 * 24; /* seconds in a day */
  int year, month, day, sec;
  REAL8 jd;

  /* this routine only works for dates after 1900 */
  if ( utc->tm_year <= 0 )
  {
    XLALPrintError( "XLAL Error - Year must be after 1900\n" );
    XLAL_ERROR_REAL8( XLAL_EDOM );
  }

  year  = utc->tm_year + 1900;
  month = utc->tm_mon + 1;     /* month is in range 1-12 */
  day   = utc->tm_mday;        /* day is in range 1-31 */
  sec   = utc->tm_sec + 60*(utc->tm_min + 60*(utc->tm_hour)); /* seconds since midnight */

  jd = 367*year - 7*(year + (month + 9)/12)/4 + 275*month/9 + day + 1721014;
  /* note: Julian days start at noon: subtract half a day */
  jd += (REAL8)sec/(REAL8)sec_per_day - 0.5;
  return jd;
}


/**
 * Returns the Modified Julian Day (MJD) corresponding to the date given
 * in a broken down time structure.
 *
 * Note:
 * - By convention, MJD is an integer.
 * - MJD number starts at midnight rather than noon.
 *
 * If you want a Modified Julian Day that has a fractional part, simply use
 * the macro:
 *
 * \#define XLAL_MODIFIED_JULIAN_DAY(utc) (XLALJulianDay(utc)-XLAL_MJD_REF)
 */
INT4 XLALModifiedJulianDay( const struct tm *utc /**< [In] UTC time in a broken down time structure. */ )
{
  REAL8 jd;
  INT4 mjd;
  jd = XLALJulianDay( utc );
  if ( XLAL_IS_REAL8_FAIL_NAN( jd ) )
    XLAL_ERROR( XLAL_EFUNC );
  mjd = floor( jd - XLAL_MJD_REF );
  return mjd;
}

/**
 * Fill in missing fields of a C 'tm' broken-down time struct.
 *
 * We want to use the C time functions in to fill in values for 'tm_wday' and
 * 'tm_yday', and normalise the ranges of the other members. mktime() does this,
 * but it also assumes local time, so that the 'tm' struct members are adjusted
 * according to the timezone. timegm() would be a more appropriate, but it seems
 * that it is not portable (BSD/Mac, but not standard GNU); neither is using the
 * 'timezone' variable to get the correct offset (works on GNU but not BSD/Mac!)
 * The method used here (idea from somewhere on the internet) should be safe.
 */
int XLALFillBrokenDownTime(struct tm *tm /**< Broken-down time struct. */) {

  /* Check input. */
  XLAL_CHECK( tm, XLAL_EINVAL );

  /* Set timezone. */
  tzset();

  /* Set daylight savings flag to zero, since we want to get the timezone
     difference against UTC. We save its initial value for use later. */
  int isdst = tm->tm_isdst;
  tm->tm_isdst = 0;

  /* Call mktime() to get a time 't1', adjusted for the timezone. */
  time_t t1 = mktime(tm);
  XLAL_CHECK( t1 >= 0, XLAL_ESYS );

  /* If original daylight savings flag was -1 (i.e. daylight savings unknown),
     save the current value of the flag for use later. */
  if (isdst < 0) {
    isdst = tm->tm_isdst;
  }

  /* Convert 't2' back into a 'tm' struct. gmtime() will preserve the timezone. */
  XLAL_CHECK( gmtime_r(&t1, tm) != NULL, XLAL_ESYS );

  /* Now call mktime() again to get time 't2', *twice* adjusted for the timezone. */
  time_t t2 = mktime(tm);
  XLAL_CHECK( t2 >= 0, XLAL_ESYS );

  /* Since 't1' has been adjusted for the timezone once, and 't2' twice, their
     difference is precisely the correct timezone difference! We substract this
     from 't1', which is now the desired time in UTC. */
  t1 -= t2 - t1;

  /* Call gmtime() to convert the desired time 't1' back into a 'tm' struct. */
  XLAL_CHECK( gmtime_r(&t1, tm) != NULL, XLAL_ESYS );

  /* Restore the daylight savings flag. */
  tm->tm_isdst = isdst;

  return XLAL_SUCCESS;

}

/*@}*/
