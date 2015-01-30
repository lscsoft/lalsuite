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

#include <lal/XLALError.h>
#include <lal/UserInputParser.h>

int main(void)
{

  // ---------- test DMS --> RAD conversion ----------
  REAL8 rads, radsRef;
  const char *dms;
  REAL8 diff, tol = 3e-15;

  dms = "-06:52:16.875";
  radsRef = -0.119927754278965;	// octave> dms_to_rad ( "-06:52:16.875" )
  XLAL_CHECK_MAIN ( XLALConvertDMStoRAD ( &rads, dms ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dms, rads, radsRef, diff, tol );

  dms = "-0:52:16.875";
  radsRef = -0.0152079991593048;	// octave> dms_to_rad ( " -0:52:16.875 " )
  XLAL_CHECK_MAIN ( XLALConvertDMStoRAD ( &rads, dms ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertDMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", dms, rads, radsRef, diff, tol );

  // ---------- test HMS --> RAD conversion ----------
  const char *hms;

  hms = "6:52:16.875";
  radsRef = 1.79891631418447;	// octave>  hms_to_rad ( "6:52:16.875" )
  XLAL_CHECK_MAIN ( XLALConvertHMStoRAD ( &rads, hms ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hms, rads, radsRef, diff, tol );

  hms = "0:52:16.875";
  radsRef = 0.228119987389571;	// octave> hms_to_rad ( "0:52:16.875" )
  XLAL_CHECK_MAIN ( XLALConvertHMStoRAD ( &rads, hms ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN ( (diff=fabs(rads - radsRef)) < tol, XLAL_ETOL, "Result XLALConvertHMStoRAD(%s)=%.16g differs from reference '%.16g' by %g > tolerance %g\n", hms, rads, radsRef, diff, tol );

  return EXIT_SUCCESS;

} // main()
