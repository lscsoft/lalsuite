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

#include <lal/LALStdio.h>
#include <lal/Date.h>
#include <lal/XLALError.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/ParseStringValue.h>

// ---------- local prototypes ----------
int test_ParseStringValue ( void );

// ==================== function definitions ====================
int main(void)
{
  // ---------- test various string-value parser functions ----------
  XLAL_CHECK_MAIN ( test_ParseStringValue() == XLAL_SUCCESS, XLAL_EFUNC );

  // check for memory leaks
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

} // main()

///
/// test various string-value parser functions:
/// XLALParseStringValueAsINT8(), XLALParseStringValueAsINT4(), XLALParseStringValueAsREAL8(),
/// XLALParseStringValueAsINT4PlusFrac()
///
int
test_ParseStringValue ( void )
{
  const char *valString;

  // ---------- XLALParseStringValueAsINT8() ----------
  INT8 valINT8, valINT8Ref;
  valString = "9223372036854775807"; // LAL_INT8_MAX
  valINT8Ref = 9223372036854775807;
  XLAL_CHECK ( XLALParseStringValueAsINT8 ( &valINT8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT8 == valINT8Ref, XLAL_ETOL, "XLALParseStringValueAsINT8(%s) failed, return = %" LAL_INT8_FORMAT "\n", valString, valINT8 );

  valString = "4294967294"; // 2 * LAL_INT4_MAX
  valINT8Ref = 4294967294;
  XLAL_CHECK ( XLALParseStringValueAsINT8 ( &valINT8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT8 == valINT8Ref, XLAL_ETOL, "XLALParseStringValueAsINT8(%s) failed, return = %" LAL_INT8_FORMAT "\n", valString, valINT8 );

  valString = "-4294967294"; // -2 * LAL_INT4_MAX
  valINT8Ref = -4294967294;
  XLAL_CHECK ( XLALParseStringValueAsINT8 ( &valINT8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT8 == valINT8Ref, XLAL_ETOL, "XLALParseStringValueAsINT8(%s) failed, return = %" LAL_INT8_FORMAT "\n", valString, valINT8 );

  // this one needs to fail!
  //valString = "18446744073709551616"; // 2 * LAL_INT8_MAX
  //XLAL_CHECK ( XLAL_SUCCESS != XLALParseStringValueAsINT8 ( &valINT8, valString ), XLAL_EFAILED, "XLALParseStringValueAsINT8() failed to catch out-of-range conversion\n" );
  //XLALPrintError ("---------- Not to worry, the above failure was on purpose: ----------\n\n");

  // ---------- XLALParseStringValueAsINT4() ----------
  INT4 valINT4, valINT4Ref;
  valString = "2147483647"; // LAL_INT4_MAX
  valINT4Ref = 2147483647;
  XLAL_CHECK ( XLALParseStringValueAsINT4 ( &valINT4, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT4 == valINT4Ref, XLAL_ETOL, "XLALParseStringValueAsINT4(%s) failed, return = %d\n", valString, valINT4 );

  valString = "-1000000";
  valINT4Ref = -1000000;
  XLAL_CHECK ( XLALParseStringValueAsINT4 ( &valINT4, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valINT4 == valINT4Ref, XLAL_ETOL, "XLALParseStringValueAsINT4(%s) failed, return = %d\n", valString, valINT4 );

  // this one needs to fail!
  //valString = "4294967294"; // 2 * LAL_INT4_MAX
  //XLAL_CHECK ( XLAL_SUCCESS != XLALParseStringValueAsINT4 ( &valINT4, valString ), XLAL_EFAILED, "XLALParseStringValueAsINT4() failed to catch out-of-range conversion\n" );
  //XLALPrintError ("---------- Not to worry, the above failure was on purpose: ----------\n\n");

  // ---------- XLALParseStringValueAsREAL8() ----------
  REAL8 valREAL8, valREAL8Ref;
  valString = "2147483647";
  valREAL8Ref = 2147483647;
  XLAL_CHECK ( XLALParseStringValueAsREAL8 ( &valREAL8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valREAL8 == valREAL8Ref, XLAL_ETOL, "XLALParseStringValueAsREAL8(%s) failed, return = %.16g\n", valString, valREAL8 );

  valString = "-1.1234e10";
  valREAL8Ref = -1.1234e10;
  XLAL_CHECK ( XLALParseStringValueAsREAL8 ( &valREAL8, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( fabs ( (valREAL8 - valREAL8Ref) / valREAL8Ref ) <= LAL_REAL8_EPS, XLAL_ETOL, "XLALParseStringValueAsREAL8(%s) failed, return = %.16g\n", valString, valREAL8 );

  // ---------- XLALParseStringValueAsREAL4() ----------
  REAL4 valREAL4, valREAL4Ref;
  valString = "2147483647";
  valREAL4Ref = 2147483647;
  XLAL_CHECK ( XLALParseStringValueAsREAL4 ( &valREAL4, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( valREAL4 == valREAL4Ref, XLAL_ETOL, "XLALParseStringValueAsREAL4(%s) failed, return = %.16g\n", valString, valREAL4 );

  valString = "-1.1234e10";
  valREAL4Ref = -1.1234e10;
  XLAL_CHECK ( XLALParseStringValueAsREAL4 ( &valREAL4, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( fabs ( (valREAL4 - valREAL4Ref) / valREAL4Ref ) <= LAL_REAL4_EPS, XLAL_ETOL, "XLALParseStringValueAsREAL4(%s) failed, return = %.16g\n", valString, valREAL4 );


  // ---------- XLALParseStringValueAsINT4PlusFrac() ----------
  INT4 valINT, valINTRef;
  REAL8 valFrac, valFracRef;

  valString = "123456789.12345678912345";
  valINTRef = 123456789;
  valFracRef = 0.12345678912345;
  XLAL_CHECK ( XLALParseStringValueAsINT4PlusFrac ( &valINT, &valFrac, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (valINT == valINTRef) && (fabs( (valFrac - valFracRef) / valFracRef ) <= LAL_REAL8_EPS), XLAL_ETOL,
               "XLALParseStringValueAsINT4PlusFrac(%s) failed, return = (%d, %.16g)\n", valString, valINT, valFrac );

  valString = "-123456789.12345678912345";
  valINTRef = -123456789;
  valFracRef = -0.12345678912345;
  XLAL_CHECK ( XLALParseStringValueAsINT4PlusFrac ( &valINT, &valFrac, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (valINT == valINTRef) && (fabs( (valFrac - valFracRef) / valFracRef ) <= LAL_REAL8_EPS), XLAL_ETOL,
               "XLALParseStringValueAsINT4PlusFrac(%s) failed, return = (%d, %.16g)\n", valString, valINT, valFrac );

  // ---------- XLALParseStringValueAsGPS() ----------
  LIGOTimeGPS valGPS, valGPSRef = {987654321, 123456789 };

  valString = "987654321.123456789";
  XLAL_CHECK ( XLALParseStringValueAsGPS ( &valGPS, valString ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}\n", valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds );

  // ---------- XLALParseStringValueAsEPOCH() ----------
  XLAL_CHECK ( XLALParseStringValueAsEPOCH ( &valGPS, valString ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}\n", valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds );

  valString = "987654321.123456789GPS";
  XLAL_CHECK ( XLALParseStringValueAsEPOCH ( &valGPS, valString ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}\n", valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds );

  valString = "55675.1848646696387616MJD";
  XLAL_CHECK ( XLALParseStringValueAsEPOCH ( &valGPS, valString ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}\n", valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds );

  valString = "987654321.12345";
  valGPSRef.gpsNanoSeconds = 123450000;
  XLAL_CHECK ( XLALParseStringValueAsEPOCH ( &valGPS, valString ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}, correct = {%d,%d}\n",
               valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds, valGPSRef.gpsSeconds, valGPSRef.gpsNanoSeconds );

  return XLAL_SUCCESS;
} // test_ParseStringValue()
