/** \file
 * \ingroup std
 * \author Chin, D. W. and Creighton, J. D. E.
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
 *   - What about time zones?
 *   - What about leap seconds?
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

#include <math.h>
#include <time.h>
#include <string.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>

#include "XLALLeapSeconds.h" /* contains the leap second table */

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
  static const char *func = "XLALLeapSeconds";
  int leap;

  if ( gpssec < leaps[0].gpssec )
  {
    XLALPrintError( "XLAL Error - Don't know leap seconds before GPS time %d\n",
        leaps[0].gpssec );
    XLAL_ERROR( func, XLAL_EDOM );
  }

  /* scan leap second table and locate the appropriate interval */
  for ( leap = 1; leap < numleaps; ++leap )
    if ( gpssec < leaps[leap].gpssec )
      break;

  return leaps[leap-1].taiutc;
}

/** Returns the leap seconds TAI-UTC for a given UTC broken down time. */
int XLALLeapSecondsUTC( const struct tm *utc /**< [In] UTC as a broken down time.*/ )
{
  static const char *func = "XLALLeapSecondsUTC";
  REAL8 jd;
  int leap;

  jd = XLALJulianDay( utc );
  if ( XLAL_IS_REAL8_FAIL_NAN( jd ) )
    XLAL_ERROR( func, XLAL_EFUNC );

  if ( jd < leaps[0].jd )
  {
    XLALPrintError( "XLAL Error - Don't know leap seconds before Julian Day %9.1f\n", leaps[0].jd );
    XLAL_ERROR( func, XLAL_EDOM );
  }

  /* scan leap second table and locate the appropriate interval */
  for ( leap = 1; leap < numleaps; ++leap )
    if ( jd < leaps[leap].jd )
      break;

  return leaps[leap-1].taiutc;
}


/** Returns the GPS seconds since the GPS epoch for a
 * specified UTC time structure. */
INT4 XLALUTCToGPS( const struct tm *utc /**< [In] UTC time in a broken down time structure. */ )
{
  static const char *func = "XLALUTCToGPS";
  time_t unixsec;
  INT4 gpssec;
  int leapsec;

  /* compute leap seconds */
  leapsec = XLALLeapSecondsUTC( utc );
  if ( leapsec < 0 )
    XLAL_ERROR( func, XLAL_EFUNC );
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


/** Returns a pointer to a tm structure representing the time
 * specified in seconds since the GPS epoch.  */
struct tm * XLALGPSToUTC(
    struct tm *utc, /**< [Out] Pointer to tm struct where result is stored. */
    INT4 gpssec /**< [In] Seconds since the GPS epoch. */
    )
{
  static const char * func = "XLALGPSToUTC";
  time_t unixsec;
  int leapsec;
  int delta;
  leapsec = XLALLeapSeconds( gpssec );
  if ( leapsec < 0 )
    XLAL_ERROR_NULL( func, XLAL_EFUNC );
  unixsec  = gpssec - leapsec + XLAL_EPOCH_GPS_TAI_UTC; /* get rid of leap seconds */
  unixsec += XLAL_EPOCH_UNIX_GPS; /* change to unix epoch */
  memset( utc, 0, sizeof( *utc ) ); /* blank out utc structure */
  utc = gmtime_r( &unixsec, utc ); /* FIXME: use gmtime ?? */
  /* now check to see if we need to add a 60th second to UTC */
  if ( ( delta = delta_tai_utc( gpssec ) ) > 0 )
    utc->tm_sec += 1; /* delta only ever is one, right?? */
  return utc;
}


/** Returns the Julian Day (JD) corresponding to the date given in a broken
 * down time structure. */
REAL8 XLALJulianDay( const struct tm *utc /**< [In] UTC time in a broken down time structure. */ )
{
  static const char * func = "XLALJulianDay";
  const int sec_per_day = 60 * 60 * 24; /* seconds in a day */
  int year, month, day, sec;
  int jd;

  /* this routine only works for dates after 1900 */
  if ( utc->tm_year <= 0 )
  {
    XLALPrintError( "XLAL Error - Year must be after 1900\n" );
    XLAL_ERROR( func, XLAL_EDOM );
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


/** Returns the Modified Julian Day (MJD) corresponding to the date given
 * in a broken down time structure.
 *
 * Note:
 *   - By convention, MJD is an integer.
 *   - MJD number starts at midnight rather than noon.
 *
 * If you want a Modified Julian Day that has a fractional part, simply use
 * the macro:
 *
 * #define XLAL_MODIFIED_JULIAN_DAY(utc) (XLALJulianDay(utc)-XLAL_MJD_REF)
 */
INT4 XLALModifiedJulianDay( const struct tm *utc /**< [In] UTC time in a broken down time structure. */ )
{
  static const char * func = "XLALModifiedJulianDay";
  REAL8 jd;
  INT4 mjd;
  jd = XLALJulianDay( utc );
  if ( XLAL_IS_REAL8_FAIL_NAN( jd ) )
    XLAL_ERROR( func, XLAL_EFUNC );
  mjd = floor( jd - XLAL_MJD_REF );
  return mjd;
}
