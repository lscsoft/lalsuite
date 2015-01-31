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
#include <lal/UserInputParser.h>

int main(void)
{
  REAL8 diff, tol = 3e-15;

  // ---------- test DMS <--> RAD conversion ----------
  REAL8 rads, radsRef;
  const char *dmsRef;
  char *dms;

  dmsRef = "-06:52:16.87500";
  radsRef = -0.119927754278965;	// octave> dms_to_rad ( "-06:52:16.875" )
  XLAL_CHECK_MAIN ( XLALConvertDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK_MAIN ( (dms = XLALConvertRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( strcmp ( dms, dmsRef ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );

  XLALPrintInfo ("Converted DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );

  dmsRef = "+00:52:16.87532";
  radsRef = 0.0152080007107085;	// octave> dms_to_rad ( "00:52:16.87532");
  XLAL_CHECK_MAIN ( XLALConvertDMStoRAD ( &rads, dmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK_MAIN ( (dms = XLALConvertRADtoDMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( strcmp ( dms, dmsRef ) == 0, XLAL_ETOL, "Returned dms string '%s' differs from input '%s'\n", dms, dmsRef );

  XLALPrintInfo ("Converted DMS '%s' into %.16g rad, and back into '%s'\n", dmsRef, rads, dms );
  // ---------- test HMS <--> RAD conversion ----------
  const char *hmsRef;
  char *hms;

  hmsRef = "06:52:16.8750000";
  radsRef = 1.79891631418447;	// octave>  hms_to_rad ( "6:52:16.875" )
  XLAL_CHECK_MAIN ( XLALConvertHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK_MAIN ( (hms = XLALConvertRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( strcmp ( hms, hmsRef ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );

  XLALPrintInfo ("Converted HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );

  hmsRef = "00:52:16.8753234";
  radsRef = 0.228120010907883;	// octave> hms_to_rad ( "00:52:16.8753234" )
  XLAL_CHECK_MAIN ( XLALConvertHMStoRAD ( &rads, hmsRef ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hmsRef, rads, radsRef, diff, tol );

  XLAL_CHECK_MAIN ( (hms = XLALConvertRADtoHMS ( rads )) != NULL, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( strcmp ( hms, hmsRef ) == 0, XLAL_ETOL, "Returned hms string '%s' differs from input '%s'\n", hms, hmsRef );

  XLALPrintInfo ("Converted HMS '%s' into %.16g rad, and back into '%s'\n", hmsRef, rads, hms );
  return EXIT_SUCCESS;

} // main()
