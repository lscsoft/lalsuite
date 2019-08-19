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
#include <lal/TranslateAngles.h>

// ==================== function definitions ====================

///
/// Translate a string representing an angle in the form "degrees:minutes:seconds" into radians. The string
/// could also just contain "degrees:minutes" or "degrees", where the absent parts are set to zero.
///
/// It requires that the minutes and seconds values are between 0 to 60. Degrees are allowed to be a
/// positive or negative integer between [-360, 360].
/// An example would be: XLALTranslateDMStoRAD ( &radians, "-06:52:16.875" );
///
int
XLALTranslateDMStoRAD ( REAL8 *radians, const CHAR *dms )
{
  XLAL_CHECK ( dms != NULL, XLAL_EINVAL, "Angle input string 'dms' is NULL" );
  XLAL_CHECK ( radians != NULL, XLAL_EINVAL );

  XLAL_CHECK ( !isspace(dms[0]), XLAL_EINVAL, "No initial whitespace allowed in input string '%s'\n", dms );
  XLAL_CHECK ( dms[strlen(dms)-1] != ':', XLAL_EINVAL, "No trailing colons allowed in input string '%s'\n", dms );

  REAL8 s = 0.;
  INT4 d = 0, m = 0;
  int numitems = sscanf(dms, "%d:%d:%lf", &d, &m, &s);

  XLAL_CHECK ( numitems == 3 || numitems == 2 || numitems == 1, XLAL_EINVAL, "Angle input string '%s' not in format 'degs:mins:secs'", dms );
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

} // XLALTranslateDMStoRAD()

///
/// Translate a string representing an angle in the form "hours:minutes:seconds" into radians. The string
/// could also just contain "hours:minutes" or "hours", where the absent parts are set to zero.
///
/// It requires that the hours value to be within [0, 23] hours, and the minutes and seconds values are within [0, 60).
/// An example would be: XLALTranslateHMStoRAD( &radians, "12:05:07.765" );
///
int
XLALTranslateHMStoRAD ( REAL8 *radians, const CHAR *hms )
{
  XLAL_CHECK_REAL8( hms != NULL, XLAL_EINVAL, "Angle input string 'hms' is NULL" );
  XLAL_CHECK ( radians != NULL, XLAL_EINVAL );

  XLAL_CHECK ( hms[strlen(hms)-1] != ':', XLAL_EINVAL, "No trailing colons allowed in input string '%s'\n", hms );

  REAL8 s = 0.;
  INT4 h = 0, m = 0;
  int numitems = sscanf(hms, "%d:%d:%lf", &h, &m, &s);

  XLAL_CHECK_REAL8 ( numitems == 3 || numitems == 2 || numitems == 1, XLAL_EINVAL, "Angle input string '%s' not in format 'hours:mins:secs'\n", hms );
  XLAL_CHECK_REAL8 ( h >= 0 && h < 24, XLAL_EDOM, "Hours value '%d' must be within [0, 23]\n", h );
  XLAL_CHECK_REAL8 ( m >= 0 && m < 60, XLAL_EDOM, "Minutes value '%d' must be within [0 to 59]\n", m );
  XLAL_CHECK_REAL8 ( s >= 0 && s < 60, XLAL_EDOM, "Seconds value '%lf' must be within [0,60)\n", s );

  /* convert from hh:mm:ss to radians */
  const REAL8 hour2deg = 360./24.;
  const REAL8 deg2rad  = LAL_PI/180.0;
  const REAL8 hour2rad = hour2deg * deg2rad;

  (*radians) = hour2rad * ( h + (m / 60.0) + (s / 3600.0) );

  return XLAL_SUCCESS;

} // XLALTranslateHMStoRAD()


///
/// Translate (longitude, right-ascencsion, RA) radians into hours:minutes:seconds (HMS) format, returns allocated string.
///
CHAR *
XLALTranslateRADtoHMS ( REAL8 radians )
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

} // XLALTranslateRADtoHMS()

///
/// Translate (latitude, declination, DEC) radians into "sign*degrees:minutes:seconds" (DMS) format, returns allocated string
///
CHAR *
XLALTranslateRADtoDMS ( REAL8 radians )
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

} // XLALTranslateRADtoDMS()
