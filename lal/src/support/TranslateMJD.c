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

#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/LALConstants.h>
#include <lal/LALString.h>

#include <lal/UserInputParse.h>

#include <lal/TranslateMJD.h>

// ---------- global variables ----------
const REAL8 TT_minus_TAI          = 32.184;	///< constant offset between TT and TAI epochs (see e.g. Table 1 in \cite SeidelmannFukushima1992,)
const REAL8 TAI_minus_UTC_at_GPS0 = 19.0;	///< offset between TAI and UTC at GPS epoch, ie number of leap seconds at GPS epoch
const REAL8 TT_minus_UTC_at_GPS0  = 51.184; 	///< = TT_minus_TAI + TAI_minus_UTC_at_GPS0, ie offset between TT and UTC at GPS epoch
const REAL8 GPS0_in_MJDUTC        = 44244.0;	///< GPS epoch [1980 JAN 6 0h UTC] expressed in MJD(UTC)

// ==================== function definitions ====================

///
/// convert given MJD(TT) time, mjd = mjdDays + mjdFracDays into LIGOTimeGPS format, preserving full (ns) accuracy.
///
/// returns gps input pointer on success, NULL on error.
///
LIGOTimeGPS *
XLALTranslateMJDTTtoGPS ( LIGOTimeGPS *gps,	///< [out] returned GPS time
                          INT4 mjdDays,		///< [in] input MJD integer days, must be >= 0
                          REAL8 mjdFracDays	///< [in] input MJD fractional days, must be in [0, 1)
                          )
{
  XLAL_CHECK_NULL ( gps != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( mjdDays >= 0, XLAL_EDOM, "mjdDays = %d must be positive\n", mjdDays );
  XLAL_CHECK_NULL ( (mjdFracDays < 1) && (mjdFracDays >= 0), XLAL_EDOM, "mjdFracDays = %g must be within [0, 1) days\n", mjdFracDays );

  INT4 gpsSeconds = (INT4) round ( (mjdDays - GPS0_in_MJDUTC) * 86400.0 );	// use gps0 epoch in mjd(utc) here, [round() just for safety, should be int anyway]
  REAL8 gpsSeconds2 = mjdFracDays * 86400.0 - TT_minus_UTC_at_GPS0;		// correct gps0 epoch from mjd(utc) to mjd(tt)

  REAL8 int2, frac2;
  frac2 = modf ( gpsSeconds2, &int2 );	// get integer and fractional parts

  gpsSeconds += (INT4) int2;
  INT4 gpsNanoSeconds = (INT4) round ( frac2 * XLAL_BILLION_REAL8 );
  if ( gpsNanoSeconds >= XLAL_BILLION_INT4 )
    {
      gpsNanoSeconds = 0;
      gpsSeconds ++;
    }
  if ( gpsNanoSeconds < 0 )
    {
      gpsSeconds --;
      gpsNanoSeconds += XLAL_BILLION_INT4;
    }

  gps->gpsSeconds = gpsSeconds;
  gps->gpsNanoSeconds = gpsNanoSeconds;

  return gps;

} // XLALTranslateMJDTTtoGPS()


///
/// Parse and convert given string representing MJD(TT) time into LIGOTimeGPS gps time, without loss of (ns) accuracy
///
/// returns gps input pointer on success, NULL on error.
///
LIGOTimeGPS *
XLALTranslateStringMJDTTtoGPS ( LIGOTimeGPS *gps,		///< [out] returned GPS time
                                const char *mjdString 	///< [in] input string representing MJD(TT) time
                                )
{
  XLAL_CHECK_NULL ( (gps != NULL) && (mjdString != NULL), XLAL_EINVAL );

  INT4 mjdDays;
  REAL8 mjdFracDays;
  XLAL_CHECK_NULL ( XLALParseStringValueAsINT4PlusFrac ( &mjdDays, &mjdFracDays, mjdString ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK_NULL ( XLALTranslateMJDTTtoGPS ( gps, mjdDays, mjdFracDays ) != NULL, XLAL_EFUNC );

  return gps;

} // XLALTranslateStringMJDTTtoGPS()
