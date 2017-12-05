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
#include <lal/TranslateAngles.h>

// ---------- local prototypes ----------
int test_HMS_RAD ( void );
int test_DMS_RAD ( void );

// ==================== function definitions ====================
int main(void)
{

  // ---------- test angle conversions ----------
  XLAL_CHECK_MAIN ( test_HMS_RAD() == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( test_DMS_RAD() == XLAL_SUCCESS, XLAL_EFUNC );

  // check for memory leaks
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

} // main()

///
/// test angle conversions between HMS and RAD: XLALTranslateHMStoRAD() and XLALTranslateRADtoHMS()
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
  XLAL_CHECK ( XLALTranslateHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (hms = XLALTranslateRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( hms, hmsRef ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );
  XLALPrintInfo ("Translated HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );
  XLALFree ( hms );

  /* test with form HH:MM */
  hmsRef = "06:52";
  radsRef = 1.797689129554159;  // pulsarpputils.ra_to_rad ( "06:52" )
  XLAL_CHECK ( XLALTranslateHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (hms = XLALTranslateRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( hms, "06:52:00.0000000" ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );
  XLALPrintInfo ("Translated HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );
  XLALFree ( hms );

  /* test with form HH */
  hmsRef = "06";
  radsRef = 1.570796326794897;  // pulsarpputils.ra_to_rad ( "06" )
  XLAL_CHECK ( XLALTranslateHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (hms = XLALTranslateRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( hms, "06:00:00.0000000" ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );
  XLALPrintInfo ("Translated HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );
  XLALFree ( hms );

  hmsRef = "00:52:16.8753234";
  radsRef = 0.228120010907883;	// octave> hms_to_rad ( "00:52:16.8753234" )
  XLAL_CHECK ( XLALTranslateHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (hms = XLALTranslateRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( hms, hmsRef ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );
  XLALPrintInfo ("Translated HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );
  XLALFree ( hms );

  /* test with form HH:MM */
  hmsRef = "00:52";
  radsRef = 0.226892802759263;	// pulsarpputils.ra_to_rad ( "00:52" )
  XLAL_CHECK ( XLALTranslateHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (hms = XLALTranslateRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( hms, "00:52:00.0000000" ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );
  XLALPrintInfo ("Translated HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );
  XLALFree ( hms );

  /* test with form HH */
  hmsRef = "00";
  radsRef = 0.0;	// pulsarpputils.ra_to_rad ( "00" )
  XLAL_CHECK ( XLALTranslateHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (hms = XLALTranslateRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( hms, "00:00:00.0000000" ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );
  XLALPrintInfo ("Translated HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );
  XLALFree ( hms );

  return XLAL_SUCCESS;
} // test_HMS_RAD()

///
/// test angle conversions between DMS and RAD: XLALTranslateDMStoRAD() and XLALTranslateRADtoDMS()
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
  XLAL_CHECK ( XLALTranslateDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (dms = XLALTranslateRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( dms, dmsRef ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );
  XLALPrintInfo ("Translated DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );
  XLALFree ( dms );

  /* test with form DD:MM */
  dmsRef = "-06:52";
  radsRef = -0.119845941970277;	// pulsarpputils.ra_to_rad ( "-06:52" )
  XLAL_CHECK ( XLALTranslateDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (dms = XLALTranslateRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( dms, "-06:52:00.00000" ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );
  XLALPrintInfo ("Translated DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );
  XLALFree ( dms );

  /* test with form DD */
  dmsRef = "-06";
  radsRef = -0.104719755119660;	// pulsarpputils.ra_to_rad ( "-06" )
  XLAL_CHECK ( XLALTranslateDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (dms = XLALTranslateRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( dms, "-06:00:00.00000" ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );
  XLALPrintInfo ("Translated DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );
  XLALFree ( dms );

  dmsRef = "+00:52:16.87532";
  radsRef = 0.0152080007107085;	// octave> dms_to_rad ( "00:52:16.87532");
  XLAL_CHECK ( XLALTranslateDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (dms = XLALTranslateRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( dms, dmsRef ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );
  XLALPrintInfo ("Translated DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );
  XLALFree ( dms );

  /* test with form DD:MM */
  dmsRef = "+00:52";
  radsRef = 0.015126186850618;	// pulsarpputils.ra_to_rad ( "+00:52" )
  XLAL_CHECK ( XLALTranslateDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (dms = XLALTranslateRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( dms, "+00:52:00.00000" ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );
  XLALPrintInfo ("Translated DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );
  XLALFree ( dms );

  /* test with form DD:MM */
  dmsRef = "+00";
  radsRef = 0.0;	// pulsarpputils.ra_to_rad ( "+00" )
  XLAL_CHECK ( XLALTranslateDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALTranslateDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK ( (dms = XLALTranslateRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( strcmp ( dms, "+00:00:00.00000" ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );
  XLALPrintInfo ("Translated DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );
  XLALFree ( dms );

  return XLAL_SUCCESS;
} // test_DMS_RAD()
