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
#include <lal/XLALError.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/TranslateMJD.h>


// ---------- local prototypes ----------
int test_MJDTT_GPS ( void );

// ==================== function definitions ====================
int main(void)
{
  // ---------- test MJD(TT) to GPS conversions ----------
  XLAL_CHECK_MAIN ( test_MJDTT_GPS() == XLAL_SUCCESS, XLAL_EFUNC );

  // check for memory leaks
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

} // main()

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

  XLAL_CHECK ( XLALTranslateMJDTTtoGPS ( &gps, mjdTTDays, mjdTTFracDays ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALTranslateMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  sprintf ( mjdTTString, "%d.%014" LAL_INT8_FORMAT, mjdTTDays, (INT8)round(mjdTTFracDays*1e14) );
  XLAL_CHECK ( XLALTranslateStringMJDTTtoGPS ( &gps, mjdTTString ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALTranslateStringMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  // ----- example 2: Chandra MET http://cxc.cfa.harvard.edu/contrib/arots/time/time_tutorial.html
  mjdTTDays = 50814;
  mjdTTFracDays = 0;
  gpsRef.gpsSeconds = 567647948;
  gpsRef.gpsNanoSeconds = 0.816 * 1e9;

  XLAL_CHECK ( XLALTranslateMJDTTtoGPS ( &gps, mjdTTDays, mjdTTFracDays ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALTranslateMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  sprintf ( mjdTTString, "%d.%014" LAL_INT8_FORMAT, mjdTTDays, (INT8)round(mjdTTFracDays*1e14) );
  XLAL_CHECK ( XLALTranslateStringMJDTTtoGPS ( &gps, mjdTTString ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALTranslateStringMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );


  // ----- example 3: RXTE MET https://heasarc.gsfc.nasa.gov/docs/xte/abc/time_tutorial.html
  mjdTTDays = 49353;
  mjdTTFracDays = 0.000696574074074074;
  gpsRef.gpsSeconds = 441417609;	// $ lalapps_tconvert -g "Jan 1 1994 0:00:00 UTC"
  gpsRef.gpsNanoSeconds = 0;

  XLAL_CHECK ( XLALTranslateMJDTTtoGPS ( &gps, mjdTTDays, mjdTTFracDays ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALTranslateMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  sprintf ( mjdTTString, "%d.%014" LAL_INT8_FORMAT, mjdTTDays, (INT8)round(mjdTTFracDays*1e14) );
  XLAL_CHECK ( XLALTranslateStringMJDTTtoGPS ( &gps, mjdTTString ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( (gps.gpsSeconds == gpsRef.gpsSeconds) && (gps.gpsNanoSeconds == gpsRef.gpsNanoSeconds), XLAL_ETOL,
               "XLALTranslateStringMJDTTtoGPS(%s) = (%d,%d) failed, correct result = (%d,%d)\n",
               mjdTTString, gps.gpsSeconds, gps.gpsNanoSeconds, gpsRef.gpsSeconds, gpsRef.gpsNanoSeconds );

  return XLAL_SUCCESS;
} // test_MJDTT_GPS()
