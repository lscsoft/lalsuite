//
// Copyright (C) 2004, 2005, 2015 Reinhard Prix
// Copyright (C) 2013 Matt Pitkin
// Copyright (C) 2007 Chris Messenger
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with with program; see the file COPYING. If not, write to the
//  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//  MA  02111-1307  USA
//

#include <math.h>
#include <ctype.h>
#include <errno.h>

#include <lal/StringInput.h>
#include <lal/LALConstants.h>
#include <lal/LALString.h>
#include <lal/UserInputParser.h>

// ---------- global variables ----------
const REAL8 TT_minus_TAI          = 32.184;	///< constant offset between TT and TAI epochs (see e.g. Table 1 in \cite SeidelmannFukushima1992,)
const REAL8 TAI_minus_UTC_at_GPS0 = 19.0;	///< offset between TAI and UTC at GPS epoch, ie number of leap seconds at GPS epoch
const REAL8 TT_minus_UTC_at_GPS0  = 51.184; 	///< = TT_minus_TAI + TAI_minus_UTC_at_GPS0, ie offset between TT and UTC at GPS epoch
const REAL8 GPS0_in_MJDUTC        = 44244.0;	///< GPS epoch [1980 JAN 6 0h UTC] expressed in MJD(UTC)

// ---------- local defines ----------
// these constants are taken from StringConvert.c
#define LAL_sINT4_MAX    LAL_INT8_C(2147483647)
#define LAL_sINT8_MAX    LAL_INT8_C(9223372036854775807)


// ==================== function definitions ====================

/// Parse a string into an INT8
/// This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueToINT8 ( INT8 *valINT8,         ///< [out] return INT8 value
                             const char *valString  ///< [in]  input string value
                             )
{
  XLAL_CHECK ( (valINT8 != NULL) && (valString != NULL ), XLAL_EINVAL );

  errno = 0;
  char *endptr;
  int base10 = 10;
  long long valLLong = strtoll ( valString, &endptr, base10 );
  XLAL_CHECK ( errno == 0, XLAL_EFAILED, "strtoll() failed to convert '%s' into long long!\n", valString );
  XLAL_CHECK ( (*endptr) == '\0', XLAL_EFAILED, "strtoll(): trailing garbage '%s' found after int-conversion of '%s'\n", endptr, valString );

  //  check range and convert long-int into INT8
  if ( sizeof(valLLong) > sizeof(INT8) ) { // avoid warning about trivial check
    XLAL_CHECK ( (valLLong >= -LAL_sINT8_MAX) && (valLLong <= LAL_sINT8_MAX), XLAL_EDOM, "String-conversion '%s' --> '%lli' exceeds INT8 range of +-%"LAL_INT8_FORMAT"\n",
                 valString, valLLong, LAL_sINT8_MAX );
  }

  (*valINT8) = (INT8)valLLong;

  return XLAL_SUCCESS;

} // XLALParseStringValueToINT8()


/// Parse a string into an INT4
/// This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueToINT4 ( INT4 *valINT4,         ///< [out] return INT4 value
                             const char *valString  ///< [in]  input string value
                             )
{
  XLAL_CHECK ( (valINT4 != NULL) && (valString != NULL ), XLAL_EINVAL );

  INT8 valINT8;
  XLAL_CHECK ( XLALParseStringValueToINT8 ( &valINT8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );

  // check range and convert INT8 into INT4
  XLAL_CHECK ( (valINT8 >= -LAL_sINT4_MAX) && (valINT8 <= LAL_sINT4_MAX), XLAL_EDOM, "String-conversion '%s' --> '%"LAL_INT8_FORMAT"' exceeds INT4 range of +-%"LAL_INT8_FORMAT"\n",
               valString, valINT8, LAL_sINT4_MAX );

  (*valINT4) = (INT4)valINT8;

  return XLAL_SUCCESS;

} // XLALParseStringValueToINT4()


/// Parse a string into a REAL8
/// This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueToREAL8 ( REAL8 *valREAL8,         ///< [out] return REAL8 value
                              const char *valString    ///< [in]  input string value
                              )
{
  XLAL_CHECK ( (valREAL8 != NULL) && (valString != NULL ), XLAL_EINVAL );

  errno = 0;
  char *endptr;
  double valDouble = strtod ( valString, &endptr );
  XLAL_CHECK ( errno == 0, XLAL_EFAILED, "strtod() failed to convert '%s' into a double!\n", valString );
  XLAL_CHECK ( (*endptr) == '\0', XLAL_EFAILED, "strtod(): trailing garbage '%s' found after double-conversion of '%s'\n", endptr, valString );

  (*valREAL8) = (REAL8)valDouble;

  return XLAL_SUCCESS;

} // XLALParseStringValueToREAL8()

/// Parse a string into a REAL4.
/// This ignores initial whitespace, but throws an error on _any_ non-converted trailing characters (including whitespace)
int
XLALParseStringValueToREAL4 ( REAL4 *valREAL4,         ///< [out] return REAL4 value
                              const char *valString    ///< [in]  input string value
                              )
{
  XLAL_CHECK ( (valREAL4 != NULL) && (valString != NULL ), XLAL_EINVAL );

  errno = 0;
  char *endptr;
  float valFloat = strtof ( valString, &endptr );
  XLAL_CHECK ( errno == 0, XLAL_EFAILED, "strtof() failed to convert '%s' into a float!\n", valString );
  XLAL_CHECK ( (*endptr) == '\0', XLAL_EFAILED, "strtof(): trailing garbage '%s' found after float-conversion of '%s'\n", endptr, valString );

  (*valREAL4) = (REAL4)valFloat;

  return XLAL_SUCCESS;

} // XLALParseStringValueToREAL4()

///
/// Parse a string containing a floating-point number into integer and fractional part, such that val = valINT + valFrac.
///
/// This is useful for parsing strings representing GPS or MJD times wihout loss of ns accuracy.
///
int
XLALParseStringValueToINT4PlusFrac ( INT4 *valINT4,		///< [out] return INT4 integer part 'xxx'
                                     REAL8 *valFrac,      	///< [out] return fractional part '0.yyyy'
                                     const char *valString	///< [in]  input string value representing a floating-point number "xxx.yyyy"
                                     )
{
  XLAL_CHECK ( (valINT4 != NULL) && (valFrac != NULL) && (valString != NULL), XLAL_EINVAL );
  XLAL_CHECK ( !isspace(valString[0]), XLAL_EINVAL, "No initial whitespace allowed in input string '%s'\n", valString );

  char buf[256];
  strncpy ( buf, valString, sizeof(buf)-1 );
  buf[ sizeof(buf)-1 ] = 0;

  REAL8 sign = 1;
  if ( valString[0] == '-' ) {	 // that's why no initial whitespace is allowed in input string
    sign = -1;
  }
  char *point = strchr ( buf, '.' );
  // is there a fractional part at all? If yes, parse it, if no, set to 0
  if ( point != NULL )
    {
      (*point) = 0;
      char fracString[256] = "0.";
      strcat ( fracString+2, point+1);
      XLAL_CHECK ( XLALParseStringValueToREAL8 ( valFrac, fracString ) == XLAL_SUCCESS, XLAL_EFUNC );
      (*valFrac) *= sign;	// correct sign: must agree with integer part
    }
  else
    {
      (*valFrac) = 0;
    }

  // now parse integer part
  XLAL_CHECK ( XLALParseStringValueToINT4 ( valINT4, buf ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;
} // XLALParseStringValueToINT4PlusFrac()


/// Parse a string into a BOOLEAN
/// Allowed string-values are (case-insensitive):
/// {"yes", "true", "1"} --> TRUE
/// {"no", "false", "0"} --> FALSE
///
/// NOTE: This throws an error on _any_ extraneous leading or trailing characters or whitespace
int
XLALParseStringValueToBOOLEAN ( BOOLEAN *valBOOLEAN,     ///< [out] return BOOLEAN value
                                const char *valString    ///< [in]  input string value
                                )
{
  XLAL_CHECK ( (valBOOLEAN != NULL) && (valString != NULL ), XLAL_EINVAL );

  /* get rid of case ambiguities */
  char *valStringLower;
  XLAL_CHECK ( (valStringLower = XLALMalloc ( strlen(valString) + 1 )) != NULL, XLAL_ENOMEM );
  strcpy ( valStringLower, valString );
  XLALStringToLowerCase ( valStringLower );

  /* parse it as a bool */
  if ( !strcmp(valStringLower, "yes") || !strcmp(valStringLower, "true") || !strcmp(valStringLower, "1") )
    {
      (*valBOOLEAN) = 1;
    }
  else if ( !strcmp(valStringLower, "no") || !strcmp(valStringLower, "false") || !strcmp(valStringLower, "0") )
    {
      (*valBOOLEAN) = 0;
    }
  else
    {
      XLALFree ( valStringLower );
      XLAL_ERROR ( XLAL_EINVAL, "Illegal bool-string '%s', needs to be one of {'yes', 'true', '1'} or {'no', 'false', '0'} (case-insensitive)\n", valString );
    }

  XLALFree ( valStringLower );

  return XLAL_SUCCESS;

} // XLALParseStringValueToBOOLEAN()


///
/// Convert a string representing an angle in the form "degrees:minutues:seconds" into radians.
///
/// It requires that the minutes and seconds values are between 0 to 60. Degrees are allowed to be a
/// positive or negative integer between [-360, 360].
/// An example would be: XLALConvertDMStoRAD ( &radians, "-06:52:16.875" );
///
int
XLALConvertDMStoRAD ( REAL8 *radians, const CHAR *dms )
{
  XLAL_CHECK ( dms != NULL, XLAL_EINVAL, "Angle input string 'dms' is NULL" );
  XLAL_CHECK ( radians != NULL, XLAL_EINVAL );

  XLAL_CHECK ( !isspace(dms[0]), XLAL_EINVAL, "No initial whitespace allowed in input string '%s'\n", dms );

  REAL8 s;
  INT4 d, m;
  int numitems = sscanf(dms, "%d:%d:%lf", &d, &m, &s);

  XLAL_CHECK ( numitems == 3, XLAL_EINVAL, "Angle input string '%s' not in format 'degs:mins:secs'", dms );
  XLAL_CHECK ( d >= -360 && d <= 360, XLAL_EDOM, "Degrees '%d' outside of valid range of [-360,360] deg\n", d );
  XLAL_CHECK ( m >= 0 && m < 60, XLAL_EDOM, "Minutes '%d' outside of the valid range of [0, 59] mins", m );
  XLAL_CHECK ( s >= 0 && s < 60, XLAL_EDOM, "Seconds '%lf' outside of the valid range of [0, 60) secs", s );

  // check if there's a minus sign, and apply to minutes and seconds (degrees would already have it)
  // Note that this is the reason we don't accept initial whitespace in the input string
  REAL8 sig = 1;
  if ( dms[0] == '-' ) {
    sig = -1;
  }

  // now convert the pieces from degrees to radians
  (*radians) =  (LAL_PI/180.0) * ( d + (sig*m / 60.0) + (sig*s / 3600.0) );

  return XLAL_SUCCESS;

} // XLALConvertDMStoRAD()

///
/// Convert a string representing an angle in the form "hours:minutes:seconds" into radians.
///
/// It requires that the hours value to be within [0, 23] hours, and the minutes and seconds values are within [0, 60).
/// An example would be: XLALConvertHMStoRAD( &radians, "12:05:07.765" );
///
int
XLALConvertHMStoRAD ( REAL8 *radians, const CHAR *hms )
{
  XLAL_CHECK_REAL8( hms != NULL, XLAL_EINVAL, "Angle input string 'hms' is NULL" );
  XLAL_CHECK ( radians != NULL, XLAL_EINVAL );

  REAL8 s;
  INT4 h, m;
  int numitems = sscanf(hms, "%d:%d:%lf", &h, &m, &s);

  XLAL_CHECK_REAL8 ( numitems == 3, XLAL_EINVAL, "Angle input string '%s' not in format 'hours:mins:secs'\n", hms );
  XLAL_CHECK_REAL8 ( h >= 0 && h < 24, XLAL_EDOM, "Hours value '%d' must be within [0, 23]\n", h );
  XLAL_CHECK_REAL8 ( m >= 0 && m < 60, XLAL_EDOM, "Minutes value '%d' must be within [0 to 59]\n", m );
  XLAL_CHECK_REAL8 ( s >= 0 && s < 60, XLAL_EDOM, "Seconds value '%lf' must be within [0,60)\n", s );

  /* convert from hh:mm:ss to radians */
  const REAL8 hour2deg = 360./24.;
  const REAL8 deg2rad  = LAL_PI/180.0;
  const REAL8 hour2rad = hour2deg * deg2rad;

  (*radians) = hour2rad * ( h + (m / 60.0) + (s / 3600.0) );

  return XLAL_SUCCESS;

} // XLALConvertHMStoRAD()


///
/// Convert (longitude, right-ascencsion, RA) radians into hours:minutes:seconds (HMS) format, returns allocated string.
///
CHAR *
XLALConvertRADtoHMS ( REAL8 radians )
{
  XLAL_CHECK_NULL ( (radians>=0.0) && (radians < LAL_TWOPI), XLAL_EDOM, "RA %g not in range [0, 2pi) rad\n", radians );

  REAL8 remainderH = radians * 24.0/LAL_TWOPI;
  INT4 hours = (INT4) floor ( remainderH );
  remainderH -= hours;
  INT4 minutes = (INT4) floor ( remainderH * 60.0 );
  remainderH -= minutes / 60.0;
  REAL8 seconds = remainderH * 3600.0;
  INT4 roundedSec = (INT4) round ( seconds * 10000000 ) / 10000000;  // round to 1e-7s accuracy
  if ( roundedSec == 60 )
    {
      seconds = 0;
      minutes ++;
      if ( minutes == 60 )
        {
          minutes = 60;
          hours ++;
          if ( hours == 24 ) {
            hours = 0;
          }
        }
    }
  CHAR hms[256];
  snprintf ( hms, sizeof(hms)-1, "%02d:%02d:%010.7f", hours, minutes, seconds );   // output format taken from tempo2

  return XLALStringDuplicate ( hms );

} // XLALConvertRADtoHMS()

///
/// Convert (latitude, declination, DEC) radians into "sign*degrees:minutes:seconds" (DMS) format, returns allocated string
///
CHAR *
XLALConvertRADtoDMS ( REAL8 radians )
{
  XLAL_CHECK_NULL ( (radians >= -LAL_PI_2) && (radians < LAL_PI_2), XLAL_EDOM, "DEC %g not in range [-pi/2, pi/2) rad\n", radians);

  CHAR sign = (radians < 0) ? '-' : '+';

  REAL8 remainderDeg = fabs ( radians * 360.0/LAL_TWOPI );

  INT4 degrees = (INT4) floor ( remainderDeg );
  remainderDeg -= degrees;
  INT4 arcmins = (INT4) floor ( remainderDeg * 60.0 );
  remainderDeg -= arcmins / 60.0;
  REAL8 arcsecs = remainderDeg * 3600.0;
  INT4 roundedArcsecs = (INT4) round ( arcsecs * 100000 ) / 100000;	// round to 1e-5 arcsec accuracy
  if ( roundedArcsecs == 60 )
    {
      arcsecs = 0;
      arcmins ++;
      if ( arcmins == 60 )
        {
          arcmins = 0;
          degrees ++;
        }
    }

  CHAR dms[256];
  snprintf ( dms, sizeof(dms)-1, "%c%02d:%02d:%08.5f", sign, degrees, arcmins, arcsecs );

  return XLALStringDuplicate ( dms );

} // XLALConvertRADtoDMS()

///
/// convert given MJD(TT) time, mjd = mjdDays + mjdFracDays into LIGOTimeGPS format, preserving full (ns) accuracy.
///
int
XLALConvertMJDTTtoGPS ( LIGOTimeGPS *gps,	///< [out] returned GPS time
                        INT4 mjdDays,		///< [in] input MJD integer days, must be >= 0
                        REAL8 mjdFracDays	///< [in] input MJD fractional days, must be in [0, 1)
                        )
{
  XLAL_CHECK ( gps != NULL, XLAL_EINVAL );
  XLAL_CHECK ( mjdDays >= 0, XLAL_EDOM, "mjdDays = %d must be positive\n", mjdDays );
  XLAL_CHECK ( (mjdFracDays < 1) && (mjdFracDays >= 0), XLAL_EDOM, "mjdFracDays = %g must be within [0, 1) days\n", mjdFracDays );

  INT4 gpsSeconds = (INT4) round ( (mjdDays - GPS0_in_MJDUTC) * 86400.0 );	// use gps0 epoch in mjd(utc) here, [round() just for safety, should be int anyway]
  REAL8 gpsSeconds2 = mjdFracDays * 86400.0 - TT_minus_UTC_at_GPS0;		// correct gps0 epoch from mjd(utc) to mjd(tt)

  REAL8 int2, frac2;
  frac2 = modf ( gpsSeconds2, &int2 );	// get integer and fractional parts

  gpsSeconds += (INT4) int2;
  INT4 gpsNanoSeconds = (INT4) round ( frac2 * 1e9 );
  if ( gpsNanoSeconds >= (INT4)1e9 )
    {
      gpsNanoSeconds = 0;
      gpsSeconds ++;
    }
  if ( gpsNanoSeconds < 0 )
    {
      gpsSeconds --;
      gpsNanoSeconds += (INT4)1e9;
    }

  gps->gpsSeconds = gpsSeconds;
  gps->gpsNanoSeconds = gpsNanoSeconds;

  return XLAL_SUCCESS;

} // XLALConvertMJDTTtoGPS()


///
/// Parse and convert given string representing MJD(TT) time into LIGOTimeGPS gps time, without loss of (ns) accuracy
///
int
XLALConvertStringMJDTTtoGPS ( LIGOTimeGPS *gps,		///< [out] returned GPS time
                              const char *mjdString 	///< [in] input string representing MJD(TT) time
                              )
{
  XLAL_CHECK ( (gps != NULL) && (mjdString != NULL), XLAL_EINVAL );

  INT4 mjdDays;
  REAL8 mjdFracDays;
  XLAL_CHECK ( XLALParseStringValueToINT4PlusFrac ( &mjdDays, &mjdFracDays, mjdString ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALConvertMJDTTtoGPS ( gps, mjdDays, mjdFracDays ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // XLALConvertStringMJDTTtoGPS()
