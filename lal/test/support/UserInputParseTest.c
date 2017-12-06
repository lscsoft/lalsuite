//
// Copyright (C) 2016 Karl Wette
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
#include <lal/StringVector.h>
#include <lal/AVFactories.h>
#include <lal/UserInputPrint.h>

#include <lal/UserInputParse.h>

// ---------- local prototypes ----------
int test_ParseStringValue ( void );
int test_ParseStringVector(void);
int test_ParseREAL8Vector(void);

// ==================== function definitions ====================
int main(void)
{
  // ---------- test various string-value parser functions ----------
  XLAL_CHECK_MAIN ( test_ParseStringValue() == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK_MAIN ( test_ParseStringVector() == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK_MAIN ( test_ParseREAL8Vector() == XLAL_SUCCESS, XLAL_EFUNC );

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
  XLAL_CHECK ( XLALParseStringValueAsGPS ( &valGPS, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}\n", valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds );

  // ---------- XLALParseStringValueAsEPOCH() ----------
  XLAL_CHECK ( XLALParseStringValueAsEPOCH ( &valGPS, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}\n", valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds );

  valString = "987654321.123456789GPS";
  XLAL_CHECK ( XLALParseStringValueAsEPOCH ( &valGPS, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}\n", valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds );

  valString = "55675.1848646696387616MJD";
  XLAL_CHECK ( XLALParseStringValueAsEPOCH ( &valGPS, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}\n", valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds );

  valString = "987654321.12345";
  valGPSRef.gpsNanoSeconds = 123450000;
  XLAL_CHECK ( XLALParseStringValueAsEPOCH ( &valGPS, valString ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALGPSCmp ( &valGPS, &valGPSRef ) == 0, XLAL_ETOL, "XLALParseStringValueAsGPS(%s) failed, return = {%d,%d}, correct = {%d,%d}\n",
               valString, valGPS.gpsSeconds, valGPS.gpsNanoSeconds, valGPSRef.gpsSeconds, valGPSRef.gpsNanoSeconds );

  // ---------- XLALParseStringValueAsREAL8Range() ----------
  REAL8Range real8Range, real8RangeRef;

  valString = "100";
  real8RangeRef[0] = 100; real8RangeRef[1] = 100;
  XLAL_CHECK( XLALParseStringValueAsREAL8Range(&real8Range, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(real8Range[0] - real8RangeRef[0])/real8RangeRef[0] <= LAL_REAL8_EPS && fabs(real8Range[1] - real8RangeRef[1])/real8RangeRef[1] <= LAL_REAL8_EPS && real8Range[0] <= real8Range[1],
              XLAL_ETOL, "XLALParseStringValueAsREAL8Range(%s) failed, return = {%g,%g}\n", valString, real8Range[0], real8Range[1] );

  valString = "150/0.5";
  real8RangeRef[0] = 150; real8RangeRef[1] = 150.5;
  XLAL_CHECK( XLALParseStringValueAsREAL8Range(&real8Range, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(real8Range[0] - real8RangeRef[0])/real8RangeRef[0] <= LAL_REAL8_EPS && fabs(real8Range[1] - real8RangeRef[1])/real8RangeRef[1] <= LAL_REAL8_EPS && real8Range[0] <= real8Range[1],
              XLAL_ETOL, "XLALParseStringValueAsREAL8Range(%s) failed, return = {%g,%g}\n", valString, real8Range[0], real8Range[1] );

  valString = "150/-0.5";
  real8RangeRef[0] = 149.5; real8RangeRef[1] = 150;
  XLAL_CHECK( XLALParseStringValueAsREAL8Range(&real8Range, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(real8Range[0] - real8RangeRef[0])/real8RangeRef[0] <= LAL_REAL8_EPS && fabs(real8Range[1] - real8RangeRef[1])/real8RangeRef[1] <= LAL_REAL8_EPS && real8Range[0] <= real8Range[1],
              XLAL_ETOL, "XLALParseStringValueAsREAL8Range(%s) failed, return = {%g,%g}\n", valString, real8Range[0], real8Range[1] );

  valString = "200,201";
  real8RangeRef[0] = 200; real8RangeRef[1] = 201;
  XLAL_CHECK( XLALParseStringValueAsREAL8Range(&real8Range, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(real8Range[0] - real8RangeRef[0])/real8RangeRef[0] <= LAL_REAL8_EPS && fabs(real8Range[1] - real8RangeRef[1])/real8RangeRef[1] <= LAL_REAL8_EPS && real8Range[0] <= real8Range[1],
              XLAL_ETOL, "XLALParseStringValueAsREAL8Range(%s) failed, return = {%g,%g}\n", valString, real8Range[0], real8Range[1] );

  valString = "203,202";
  real8RangeRef[0] = 202; real8RangeRef[1] = 203;
  XLAL_CHECK( XLALParseStringValueAsREAL8Range(&real8Range, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(real8Range[0] - real8RangeRef[0])/real8RangeRef[0] <= LAL_REAL8_EPS && fabs(real8Range[1] - real8RangeRef[1])/real8RangeRef[1] <= LAL_REAL8_EPS && real8Range[0] <= real8Range[1],
              XLAL_ETOL, "XLALParseStringValueAsREAL8Range(%s) failed, return = {%g,%g}\n", valString, real8Range[0], real8Range[1] );

  valString = "-203,-202";
  real8RangeRef[0] = -203; real8RangeRef[1] = -202;
  XLAL_CHECK( XLALParseStringValueAsREAL8Range(&real8Range, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(real8Range[0] - real8RangeRef[0])/real8RangeRef[0] <= LAL_REAL8_EPS && fabs(real8Range[1] - real8RangeRef[1])/real8RangeRef[1] <= LAL_REAL8_EPS && real8Range[0] <= real8Range[1],
              XLAL_ETOL, "XLALParseStringValueAsREAL8Range(%s) failed, return = {%g,%g}\n", valString, real8Range[0], real8Range[1] );

  valString = "250~5";
  real8RangeRef[0] = 245; real8RangeRef[1] = 255;
  XLAL_CHECK( XLALParseStringValueAsREAL8Range(&real8Range, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(real8Range[0] - real8RangeRef[0])/real8RangeRef[0] <= LAL_REAL8_EPS && fabs(real8Range[1] - real8RangeRef[1])/real8RangeRef[1] <= LAL_REAL8_EPS && real8Range[0] <= real8Range[1],
              XLAL_ETOL, "XLALParseStringValueAsREAL8Range(%s) failed, return = {%g,%g}\n", valString, real8Range[0], real8Range[1] );

  // ---------- XLALParseStringValueAsEPOCHRange() ----------
  LIGOTimeGPSRange gpsRange, gpsRangeRef;

  valString = "100200300.4";
  XLALGPSSet(&gpsRangeRef[0], 100200300, 400000000); XLALGPSSet(&gpsRangeRef[1], 100200300, 400000000);
  XLAL_CHECK( XLALParseStringValueAsEPOCHRange(&gpsRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALGPSCmp(&gpsRange[0], &gpsRangeRef[0]) == 0 && XLALGPSCmp(&gpsRange[1], &gpsRangeRef[1]) == 0 && XLALGPSCmp(&gpsRange[0], &gpsRange[1]) <= 0, XLAL_ETOL,
              "XLALParseStringValueAsEPOCHRange(%s) failed, return = {{%d,%d},{%d,%d}}\n", valString, gpsRange[0].gpsSeconds, gpsRange[0].gpsNanoSeconds, gpsRange[1].gpsSeconds, gpsRange[1].gpsNanoSeconds );

  valString = "100200300.4/800600400.2";
  XLALGPSSet(&gpsRangeRef[0], 100200300, 400000000); XLALGPSSet(&gpsRangeRef[1], 900800700, 600000000);
  XLAL_CHECK( XLALParseStringValueAsEPOCHRange(&gpsRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALGPSCmp(&gpsRange[0], &gpsRangeRef[0]) == 0 && XLALGPSCmp(&gpsRange[1], &gpsRangeRef[1]) == 0 && XLALGPSCmp(&gpsRange[0], &gpsRange[1]) <= 0, XLAL_ETOL,
              "XLALParseStringValueAsEPOCHRange(%s) failed, return = {{%d,%d},{%d,%d}}\n", valString, gpsRange[0].gpsSeconds, gpsRange[0].gpsNanoSeconds, gpsRange[1].gpsSeconds, gpsRange[1].gpsNanoSeconds );

  valString = "100200300.4/-800600400.2";
  XLALGPSSet(&gpsRangeRef[0], -700400099, -800000000); XLALGPSSet(&gpsRangeRef[1], 100200300, 400000000);
  XLAL_CHECK( XLALParseStringValueAsEPOCHRange(&gpsRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALGPSCmp(&gpsRange[0], &gpsRangeRef[0]) == 0 && XLALGPSCmp(&gpsRange[1], &gpsRangeRef[1]) == 0 && XLALGPSCmp(&gpsRange[0], &gpsRange[1]) <= 0, XLAL_ETOL,
              "XLALParseStringValueAsEPOCHRange(%s) failed, return = {{%d,%d},{%d,%d}}\n", valString, gpsRange[0].gpsSeconds, gpsRange[0].gpsNanoSeconds, gpsRange[1].gpsSeconds, gpsRange[1].gpsNanoSeconds );

  valString = "200300400.5,600700800.9";
  XLALGPSSet(&gpsRangeRef[0], 200300400, 500000000); XLALGPSSet(&gpsRangeRef[1], 600700800, 900000000);
  XLAL_CHECK( XLALParseStringValueAsEPOCHRange(&gpsRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALGPSCmp(&gpsRange[0], &gpsRangeRef[0]) == 0 && XLALGPSCmp(&gpsRange[1], &gpsRangeRef[1]) == 0 && XLALGPSCmp(&gpsRange[0], &gpsRange[1]) <= 0, XLAL_ETOL,
              "XLALParseStringValueAsEPOCHRange(%s) failed, return = {{%d,%d},{%d,%d}}\n", valString, gpsRange[0].gpsSeconds, gpsRange[0].gpsNanoSeconds, gpsRange[1].gpsSeconds, gpsRange[1].gpsNanoSeconds );

  valString = "600700800.9,200300400.5";
  XLALGPSSet(&gpsRangeRef[0], 200300400, 500000000); XLALGPSSet(&gpsRangeRef[1], 600700800, 900000000);
  XLAL_CHECK( XLALParseStringValueAsEPOCHRange(&gpsRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALGPSCmp(&gpsRange[0], &gpsRangeRef[0]) == 0 && XLALGPSCmp(&gpsRange[1], &gpsRangeRef[1]) == 0 && XLALGPSCmp(&gpsRange[0], &gpsRange[1]) <= 0, XLAL_ETOL,
              "XLALParseStringValueAsEPOCHRange(%s) failed, return = {{%d,%d},{%d,%d}}\n", valString, gpsRange[0].gpsSeconds, gpsRange[0].gpsNanoSeconds, gpsRange[1].gpsSeconds, gpsRange[1].gpsNanoSeconds );

  valString = "123456789~123456.789";
  XLALGPSSet(&gpsRangeRef[0], 123333332, 211000000); XLALGPSSet(&gpsRangeRef[1], 123580245, 789000000);
  XLAL_CHECK( XLALParseStringValueAsEPOCHRange(&gpsRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALGPSCmp(&gpsRange[0], &gpsRangeRef[0]) == 0 && XLALGPSCmp(&gpsRange[1], &gpsRangeRef[1]) == 0 && XLALGPSCmp(&gpsRange[0], &gpsRange[1]) <= 0, XLAL_ETOL,
              "XLALParseStringValueAsEPOCHRange(%s) failed, return = {{%d,%d},{%d,%d}}\n", valString, gpsRange[0].gpsSeconds, gpsRange[0].gpsNanoSeconds, gpsRange[1].gpsSeconds, gpsRange[1].gpsNanoSeconds );

  // ---------- XLALParseStringValueAsRAJRange() ----------
  REAL8Range rajRange, rajRangeRef;

  valString = "0.1";
  rajRangeRef[0] = 0.1; rajRangeRef[1] = 0.1;
  XLAL_CHECK( XLALParseStringValueAsRAJRange(&rajRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(rajRange[0] - rajRangeRef[0])/rajRangeRef[0] <= LAL_REAL8_EPS && fabs(rajRange[1] - rajRangeRef[1])/rajRangeRef[1] <= LAL_REAL8_EPS && rajRange[0] <= rajRange[1],
              XLAL_ETOL, "XLALParseStringValueAsRAJRange(%s) failed, return = {%g,%g}\n", valString, rajRange[0], rajRange[1] );

  valString = "10:20:30/0.5";
  rajRangeRef[0] = 2.7074420021562036; rajRangeRef[1] = rajRangeRef[0] + 0.5;
  XLAL_CHECK( XLALParseStringValueAsRAJRange(&rajRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(rajRange[0] - rajRangeRef[0])/rajRangeRef[0] <= LAL_REAL8_EPS && fabs(rajRange[1] - rajRangeRef[1])/rajRangeRef[1] <= LAL_REAL8_EPS && rajRange[0] <= rajRange[1],
              XLAL_ETOL, "XLALParseStringValueAsRAJRange(%s) failed, return = {%g,%g}\n", valString, rajRange[0], rajRange[1] );

  valString = "10:20:30/-0.5";
  rajRangeRef[0] = 2.7074420021562036 - 0.5; rajRangeRef[1] = rajRangeRef[0] + 0.5;
  XLAL_CHECK( XLALParseStringValueAsRAJRange(&rajRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(rajRange[0] - rajRangeRef[0])/rajRangeRef[0] <= LAL_REAL8_EPS && fabs(rajRange[1] - rajRangeRef[1])/rajRangeRef[1] <= LAL_REAL8_EPS && rajRange[0] <= rajRange[1],
              XLAL_ETOL, "XLALParseStringValueAsRAJRange(%s) failed, return = {%g,%g}\n", valString, rajRange[0], rajRange[1] );

  valString = "10:20:30,11:22:33";
  rajRangeRef[0] = 2.7074420021562036; rajRangeRef[1] = 2.9781862023718242;
  XLAL_CHECK( XLALParseStringValueAsRAJRange(&rajRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(rajRange[0] - rajRangeRef[0])/rajRangeRef[0] <= LAL_REAL8_EPS && fabs(rajRange[1] - rajRangeRef[1])/rajRangeRef[1] <= LAL_REAL8_EPS && rajRange[0] <= rajRange[1],
              XLAL_ETOL, "XLALParseStringValueAsRAJRange(%s) failed, return = {%g,%g}\n", valString, rajRange[0], rajRange[1] );

  valString = "11:22:33,10:20:30";
  rajRangeRef[0] = 2.7074420021562036; rajRangeRef[1] = 2.9781862023718242;
  XLAL_CHECK( XLALParseStringValueAsRAJRange(&rajRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(rajRange[0] - rajRangeRef[0])/rajRangeRef[0] <= LAL_REAL8_EPS && fabs(rajRange[1] - rajRangeRef[1])/rajRangeRef[1] <= LAL_REAL8_EPS && rajRange[0] <= rajRange[1],
              XLAL_ETOL, "XLALParseStringValueAsRAJRange(%s) failed, return = {%g,%g}\n", valString, rajRange[0], rajRange[1] );

  valString = "11:22:33~00:00:44.55";
  rajRangeRef[0] = 2.9749464349478099; rajRangeRef[1] = 2.9814259697958385;
  XLAL_CHECK( XLALParseStringValueAsRAJRange(&rajRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(rajRange[0] - rajRangeRef[0])/rajRangeRef[0] <= LAL_REAL8_EPS && fabs(rajRange[1] - rajRangeRef[1])/rajRangeRef[1] <= LAL_REAL8_EPS && rajRange[0] <= rajRange[1],
              XLAL_ETOL, "XLALParseStringValueAsRAJRange(%s) failed, return = {%g,%g}\n", valString, rajRange[0], rajRange[1] );

  // ---------- XLALParseStringValueAsDECJRange() ----------
  REAL8Range decjRange, decjRangeRef;

  valString = "0.1";
  decjRangeRef[0] = 0.1; decjRangeRef[1] = 0.1;
  XLAL_CHECK( XLALParseStringValueAsDECJRange(&decjRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(decjRange[0] - decjRangeRef[0])/decjRangeRef[0] <= LAL_REAL8_EPS && fabs(decjRange[1] - decjRangeRef[1])/decjRangeRef[1] <= LAL_REAL8_EPS && decjRange[0] <= decjRange[1],
              XLAL_ETOL, "XLALParseStringValueAsDECJRange(%s) failed, return = {%g,%g}\n", valString, decjRange[0], decjRange[1] );

  valString = "10:20:30/0.5";
  decjRangeRef[0] = (2.7074420021562036/15); decjRangeRef[1] = decjRangeRef[0] + 0.5;
  XLAL_CHECK( XLALParseStringValueAsDECJRange(&decjRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(decjRange[0] - decjRangeRef[0])/decjRangeRef[0] <= LAL_REAL8_EPS && fabs(decjRange[1] - decjRangeRef[1])/decjRangeRef[1] <= LAL_REAL8_EPS && decjRange[0] <= decjRange[1],
              XLAL_ETOL, "XLALParseStringValueAsDECJRange(%s) failed, return = {%g,%g}\n", valString, decjRange[0], decjRange[1] );

  valString = "10:20:30/-0.5";
  decjRangeRef[0] = (2.7074420021562036/15) - 0.5; decjRangeRef[1] = decjRangeRef[0] + 0.5;
  XLAL_CHECK( XLALParseStringValueAsDECJRange(&decjRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(decjRange[0] - decjRangeRef[0])/decjRangeRef[0] <= LAL_REAL8_EPS && fabs(decjRange[1] - decjRangeRef[1])/decjRangeRef[1] <= LAL_REAL8_EPS && decjRange[0] <= decjRange[1],
              XLAL_ETOL, "XLALParseStringValueAsDECJRange(%s) failed, return = {%g,%g}\n", valString, decjRange[0], decjRange[1] );

  valString = "10:20:30,11:22:33";
  decjRangeRef[0] = (2.7074420021562036/15); decjRangeRef[1] = (2.9781862023718242/15);
  XLAL_CHECK( XLALParseStringValueAsDECJRange(&decjRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(decjRange[0] - decjRangeRef[0])/decjRangeRef[0] <= LAL_REAL8_EPS && fabs(decjRange[1] - decjRangeRef[1])/decjRangeRef[1] <= LAL_REAL8_EPS && decjRange[0] <= decjRange[1],
              XLAL_ETOL, "XLALParseStringValueAsDECJRange(%s) failed, return = {%g,%g}\n", valString, decjRange[0], decjRange[1] );

  valString = "11:22:33,10:20:30";
  decjRangeRef[0] = (2.7074420021562036/15); decjRangeRef[1] = (2.9781862023718242/15);
  XLAL_CHECK( XLALParseStringValueAsDECJRange(&decjRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(decjRange[0] - decjRangeRef[0])/decjRangeRef[0] <= LAL_REAL8_EPS && fabs(decjRange[1] - decjRangeRef[1])/decjRangeRef[1] <= LAL_REAL8_EPS && decjRange[0] <= decjRange[1],
              XLAL_ETOL, "XLALParseStringValueAsDECJRange(%s) failed, return = {%g,%g}\n", valString, decjRange[0], decjRange[1] );

  valString = "11:22:33~00:00:44.55";
  decjRangeRef[0] = (2.9749464349478099/15); decjRangeRef[1] = (2.9814259697958385/15);
  XLAL_CHECK( XLALParseStringValueAsDECJRange(&decjRange, valString) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( fabs(decjRange[0] - decjRangeRef[0])/decjRangeRef[0] <= LAL_REAL8_EPS && fabs(decjRange[1] - decjRangeRef[1])/decjRangeRef[1] <= LAL_REAL8_EPS && decjRange[0] <= decjRange[1],
              XLAL_ETOL, "XLALParseStringValueAsDECJRange(%s) failed, return = {%g,%g}\n", valString, decjRange[0], decjRange[1] );

  return XLAL_SUCCESS;
} // test_ParseStringValue()

// test string-vector parsing function XLALParseStringValueAsStringVector()
int
test_ParseStringVector(void)
{
#define STR1 "Hello, world!"
#define STR2 "xyda 3!#4134"
#define STR3 "&\\//.. :: some junk"
#define STR4 "H1"
#define STR5 "H2"
#define STR6 "L1"

  LALStringVector *strVect1;
  XLAL_CHECK ( (strVect1 = XLALCreateStringVector ( STR1, STR2, STR3, STR4, STR5, NULL )) != NULL, XLAL_EFUNC );

  XLAL_CHECK ( (strVect1 = XLALAppendString2Vector ( strVect1, STR6 )) != NULL, XLAL_EFUNC );

  // now 'print' this string-vector as a 'string-value', then re-parse back into a vector:
  CHAR *strValue1 = NULL;
  LALStringVector *strVect2 = NULL;
  XLAL_CHECK ( (strValue1 = XLALPrintStringValueOfSTRINGVector ( &strVect1 )) != NULL, XLAL_EFUNC );
  XLALPrintInfo ("String value of initial string-vector:   %s\n", strValue1 );

  XLAL_CHECK ( XLALParseStringValueAsSTRINGVector ( &strVect2, strValue1 ) == XLAL_SUCCESS, XLAL_EFUNC );
  CHAR *strValue2 = NULL;
  XLAL_CHECK ( (strValue2 = XLALPrintStringValueOfSTRINGVector ( &strVect2 )) != NULL, XLAL_EFUNC );
  XLALPrintInfo ("String value of re-parsed string-vector: %s\n", strValue2 );

  // ----- compare results
  // 1) compare string values
  XLAL_CHECK ( strcmp ( strValue1, strValue2 ) == 0, XLAL_EFAILED, "String values differ:\nstrValue1 = %s\nstrValue2 = %s\n", strValue1, strValue2 );

  // 2) compare string vectors
  UINT4 len1 = strVect1->length;
  UINT4 len2 = strVect2->length;
  XLAL_CHECK ( len1 == len2, XLAL_EFAILED, "String vectors vect1 and vect2 have different lengths (%d != %d )\n", len1, len2 );

  for ( UINT4 i = 0; i < len1; i++ )
    {
      if ( strcmp ( strVect1->data[i], strVect2->data[i] ) != 0 )
        {
          for ( UINT4 j=0; j < len1; j ++ ) {
            XLALPrintError ("j = %d:  s1[j] = %6s, s2[j] = %6s\n", j, strVect1->data[j], strVect2->data[j] );
          }
          XLAL_ERROR ( XLAL_EFAILED, "Printed and re-parsed string-vector differ!\n" );
        } // if s1[i] != s2[i]

    } // for i < len

  // clean up memory
  XLALFree ( strValue1 );
  XLALFree ( strValue2 );

  XLALDestroyStringVector ( strVect1  );
  XLALDestroyStringVector ( strVect2 );

  return XLAL_SUCCESS;

} // test_ParseStringVector()


// test string-vector parsing function XLALParseStringValueAsStringVector()
int
test_ParseREAL8Vector(void)
{
  const char *csvIn = "0.1,5,-5.1e+99,1.23456789e-99,inf,nan";
  REAL8 vals[] = {0.1, 5, -5.1e+99, 1.23456789e-99 }; // only finite values for comparison

  // parse csv string as REAL8Vector:
  REAL8Vector *vect1 = NULL;
  XLAL_CHECK ( XLALParseStringValueAsREAL8Vector ( &vect1, csvIn ) == XLAL_SUCCESS, XLAL_EFUNC );

  // test1: re-print as string, compare strings
  char *csvOut;
  XLAL_CHECK ( (csvOut = XLALPrintStringValueOfREAL8Vector ( &vect1 )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( csvIn, csvOut ) == 0, XLAL_EFAILED, "csvIn != csvOut:\ncsvIn  = %s\ncsvOut = %s\n", csvIn, csvOut );

  // test2: compare finite parsed values:
  for ( UINT4 i = 0; i < 4; i ++ ) {
    XLAL_CHECK ( vect1->data[i] == vals[i], XLAL_EFAILED, "Parsed %d-th value differs from input: %.16g != %.16g\n", i, vect1->data[i], vals[i] );
  }
  // check non-finite values
  XLAL_CHECK ( fpclassify ( vect1->data[4] ) == FP_INFINITE, XLAL_EFAILED, "Failed to parse 'inf'\n");
  XLAL_CHECK ( fpclassify ( vect1->data[5] ) == FP_NAN, XLAL_EFAILED, "Failed to parse 'nan'\n");

  // clean up memory
  XLALFree ( csvOut );
  XLALDestroyREAL8Vector ( vect1  );

  return XLAL_SUCCESS;

} // test_ParseREAL8Vector()
