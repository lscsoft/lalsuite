//
// Copyright (C) 2015 Reinhard Prix
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

// Tests of the unit conversion functions in the UserInputParser.[ch] module

#include <stdio.h>
#include <math.h>
#include <string.h>

#include <lal/XLALError.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/UserInputParser.h>


// ---------- local prototypes ----------
int test_HMS_RAD ( void );
int test_DMS_RAD ( void );
int test_ParseStringValue ( void );
int test_MJDTT_GPS ( void );

// ==================== function definitions ====================
int main(void)
{

  // ---------- test angle conversions ----------
  XLAL_CHECK_MAIN ( test_HMS_RAD() == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( test_DMS_RAD() == XLAL_SUCCESS, XLAL_EFUNC );

  // ---------- test various string-value parser functions ----------
  XLAL_CHECK_MAIN ( test_ParseStringValue() == XLAL_SUCCESS, XLAL_EFUNC );

  // ---------- test MJD(TT) to GPS conversions ----------
  XLAL_CHECK_MAIN ( test_MJDTT_GPS() == XLAL_SUCCESS, XLAL_EFUNC );

  // check for memory leaks
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

} // main()

///
/// test angle conversions between HMS and RAD: XLALConvertHMStoRAD() and XLALConvertRADtoHMS()
///
int
test_HMS_RAD ( void )
{
  REAL8 diff, tol = 3e-15;
  REAL8 rads, radsRef;
  const char *hmsRef;
  char *hms;

  hmsRef = "06:52:16.8750000";
  radsRef = 1.79891631418447;	// octave>  hms_to_rad ( "6:52:16.875" )
  XLAL_CHECK ( XLALConvertHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (hms = XLALConvertRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( hms, hmsRef ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );
  XLALPrintInfo ("Converted HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );
  XLALFree ( hms );

  hmsRef = "00:52:16.8753234";
  radsRef = 0.228120010907883;	// octave> hms_to_rad ( "00:52:16.8753234" )
  XLAL_CHECK ( XLALConvertHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (hms = XLALConvertRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( hms, hmsRef ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );
  XLALPrintInfo ("Converted HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );
  XLALFree ( hms );

  return XLAL_SUCCESS;
} // test_HMS_RAD()

///
/// test angle conversions between DMS and RAD: XLALConvertDMStoRAD() and XLALConvertRADtoDMS()
///
int
test_DMS_RAD ( void )
{
  REAL8 diff, tol = 3e-15;
  REAL8 rads, radsRef;
  const char *dmsRef;
  char *dms;

  dmsRef = "-06:52:16.87500";
  radsRef = -0.119927754278965;	// octave> dms_to_rad ( "-06:52:16.875" )
  XLAL_CHECK ( XLALConvertDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (dms = XLALConvertRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( dms, dmsRef ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );
  XLALPrintInfo ("Converted DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );
  XLALFree ( dms );

  dmsRef = "+00:52:16.87532";
  radsRef = 0.0152080007107085;	// octave> dms_to_rad ( "00:52:16.87532");
  XLAL_CHECK ( XLALConvertDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (dms = XLALConvertRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( dms, dmsRef ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );
  XLALPrintInfo ("Converted DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );
  XLALFree ( dms );

  return XLAL_SUCCESS;
} // test_DMS_RAD()


///
/// test various string-value parser functions:
/// XLALParseStringValueToINT8(), XLALParseStringValueToINT4(), XLALParseStringValueToREAL8(),
/// XLALParseStringValueToINT4PlusFrac()
///
int
test_ParseStringValue ( void )
{
  const char *valString;

  // ---------- XLALParseStringValueToINT8() ----------
  INT8 valINT8, valINT8Ref;
  valString = "9223372036854775807"; // LAL_INT8_MAX
  valINT8Ref = 9223372036854775807;
  XLAL_CHECK ( XLALParseStringValueToINT8 ( &valINT8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT8 == valINT8Ref, XLAL_ETOL, "XLALParseStringValueToINT8(%s) failed, return = %ld\n", valString, valINT8 );

  valString = "4294967294"; // 2 * LAL_INT4_MAX
  valINT8Ref = 4294967294;
  XLAL_CHECK ( XLALParseStringValueToINT8 ( &valINT8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT8 == valINT8Ref, XLAL_ETOL, "XLALParseStringValueToINT8(%s) failed, return = %ld\n", valString, valINT8 );

  valString = "-4294967294"; // -2 * LAL_INT4_MAX
  valINT8Ref = -4294967294;
  XLAL_CHECK ( XLALParseStringValueToINT8 ( &valINT8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT8 == valINT8Ref, XLAL_ETOL, "XLALParseStringValueToINT8(%s) failed, return = %ld\n", valString, valINT8 );

  // this one needs to fail!
  valString = "18446744073709551616"; // 2 * LAL_INT8_MAX
  XLAL_CHECK ( XLAL_SUCCESS != XLALParseStringValueToINT8 ( &valINT8, valString ), XLAL_EFAILED, "XLALParseStringValueToINT8() failed to catch out-of-range conversion\n" );
  XLALPrintError ("---------- Not to worry, the above failure was on purpose: ----------\n\n");

  // ---------- XLALParseStringValueToINT4() ----------
  INT4 valINT4, valINT4Ref;
  valString = "2147483647"; // LAL_INT4_MAX
  valINT4Ref = 2147483647;
  XLAL_CHECK ( XLALParseStringValueToINT4 ( &valINT4, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT4 == valINT4Ref, XLAL_ETOL, "XLALParseStringValueToINT4(%s) failed, return = %d\n", valString, valINT4 );

  valString = "-1000000";
  valINT4Ref = -1000000;
  XLAL_CHECK ( XLALParseStringValueToINT4 ( &valINT4, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT4 == valINT4Ref, XLAL_ETOL, "XLALParseStringValueToINT4(%s) failed, return = %d\n", valString, valINT4 );

  // this one needs to fail!
  valString = "4294967294"; // 2 * LAL_INT4_MAX
  XLAL_CHECK ( XLAL_SUCCESS != XLALParseStringValueToINT4 ( &valINT4, valString ), XLAL_EFAILED, "XLALParseStringValueToINT4() failed to catch out-of-range conversion\n" );
  XLALPrintError ("---------- Not to worry, the above failure was on purpose: ----------\n\n");

  // ---------- XLALParseStringValueToREAL8() ----------
  REAL8 valREAL8, valREAL8Ref;
  valString = "2147483647";
  valREAL8Ref = 2147483647;
  XLAL_CHECK ( XLALParseStringValueToREAL8 ( &valREAL8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valREAL8 == valREAL8Ref, XLAL_ETOL, "XLALParseStringValueToREAL8(%s) failed, return = %.16g\n", valString, valREAL8 );

  valString = "-1.1234e10";
  valREAL8Ref = -1.1234e10;
  XLAL_CHECK ( XLALParseStringValueToREAL8 ( &valREAL8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( fabs ( (valREAL8 - valREAL8Ref) / valREAL8Ref ) <= LAL_REAL8_EPS, XLAL_ETOL, "XLALParseStringValueToREAL8(%s) failed, return = %.16g\n", valString, valREAL8 );

  // ---------- XLALParseStringValueToREAL4() ----------
  REAL4 valREAL4, valREAL4Ref;
  valString = "2147483647";
  valREAL4Ref = 2147483647;
  XLAL_CHECK ( XLALParseStringValueToREAL4 ( &valREAL4, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valREAL4 == valREAL4Ref, XLAL_ETOL, "XLALParseStringValueToREAL4(%s) failed, return = %.16g\n", valString, valREAL4 );

  valString = "-1.1234e10";
  valREAL4Ref = -1.1234e10;
  XLAL_CHECK ( XLALParseStringValueToREAL4 ( &valREAL4, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( fabs ( (valREAL4 - valREAL4Ref) / valREAL4Ref ) <= LAL_REAL4_EPS, XLAL_ETOL, "XLALParseStringValueToREAL4(%s) failed, return = %.16g\n", valString, valREAL4 );


  // ---------- XLALParseStringValueToINT4PlusFrac() ----------
  INT4 valINT, valINTRef;
  REAL8 valFrac, valFracRef;

  valString = "123456789.12345678912345";
  valINTRef = 123456789;
  valFracRef = 0.12345678912345;
  XLAL_CHECK ( XLALParseStringValueToINT4PlusFrac ( &valINT, &valFrac, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (valINT == valINTRef) && (fabs( (valFrac - valFracRef) / valFracRef ) <= LAL_REAL8_EPS), XLAL_ETOL,
               "XLALParseStringValueToINT4PlusFrac(%s) failed, return = (%d, %.16g)\n", valString, valINT, valFrac );

  valString = "-123456789.12345678912345";
  valINTRef = -123456789;
  valFracRef = -0.12345678912345;
  XLAL_CHECK ( XLALParseStringValueToINT4PlusFrac ( &valINT, &valFrac, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (valINT == valINTRef) && (fabs( (valFrac - valFracRef) / valFracRef ) <= LAL_REAL8_EPS), XLAL_ETOL,
               "XLALParseStringValueToINT4PlusFrac(%s) failed, return = (%d, %.16g)\n", valString, valINT, valFrac );

  return XLAL_SUCCESS;
} // test_ParseStringValue()

///
/// test conversion between MJD(TT) string and GPS value
///
int
test_MJDTT_GPS ( void )
{

  INT4 mjdTTDays;
  REAL8 mjdTTFracDays;
  char mjdTTString[256];
  LIGOTimeGPS gps, gpsRef;

  // ----- example 1: J200 epoch, see https://en.wikipedia.org/wiki/Epoch_%28astronomy%29#Julian_years_and_J2000
  mjdTTDays = 51544;
  mjdTTFracDays = 0.5;
  gpsRef.gpsSeconds = 630763148; // $ lalapps_tconvert "Jan 01 2000 11:58:55 UTC"
  gpsRef.gpsNanoSeconds = 0.816 * 1e9;

  XLAL_CHECK ( XLALConvertMJDTTtoGPS ( &gps, mjdTTDays, mjdTTFracDays ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALConvertMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  sprintf ( mjdTTString, "%d.%014ld", mjdTTDays, (long)round(mjdTTFracDays*1e14) );
  XLAL_CHECK ( XLALConvertStringMJDTTtoGPS ( &gps, mjdTTString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALConvertStringMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  // ----- example 2: Chandra MET http://cxc.cfa.harvard.edu/contrib/arots/time/time_tutorial.html
  mjdTTDays = 50814;
  mjdTTFracDays = 0;
  gpsRef.gpsSeconds = 567647948;
  gpsRef.gpsNanoSeconds = 0.816 * 1e9;

  XLAL_CHECK ( XLALConvertMJDTTtoGPS ( &gps, mjdTTDays, mjdTTFracDays ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALConvertMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  sprintf ( mjdTTString, "%d.%014ld", mjdTTDays, (long)round(mjdTTFracDays*1e14) );
  XLAL_CHECK ( XLALConvertStringMJDTTtoGPS ( &gps, mjdTTString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALConvertStringMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );


  // ----- example 3: RXTE MET https://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html
  mjdTTDays = 49353;
  mjdTTFracDays = 0.000696574074074074;
  gpsRef.gpsSeconds = 441417609;	// $ lalapps_tconvert -g "Jan 1 1994 0:00:00 UTC"
  gpsRef.gpsNanoSeconds = 0;

  XLAL_CHECK ( XLALConvertMJDTTtoGPS ( &gps, mjdTTDays, mjdTTFracDays ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALConvertMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  sprintf ( mjdTTString, "%d.%014ld", mjdTTDays, (long)round(mjdTTFracDays*1e14) );
  XLAL_CHECK ( XLALConvertStringMJDTTtoGPS ( &gps, mjdTTString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALConvertStringMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  return XLAL_SUCCESS;
} // test_MJDTT_GPS()
