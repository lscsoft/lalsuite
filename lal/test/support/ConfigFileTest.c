/*
*  Copyright (C) 2010 Larne Pekowsky
*  based on code (C) 2009 Reinhard Prix
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

#include <math.h>

#include <lal/AVFactories.h>
#include <lal/ConfigFile.h>
#include <lal/StringVector.h>
#include <lal/Date.h>

/* Default parameters. */

int main ( int argc, char *argv[])
{
  XLAL_CHECK ( argc == 1, XLAL_EINVAL, "Command-line arguments useless here \n");
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL );

  LALParsedDataFile *cfgdata = NULL;

  BOOLEAN testBool;
  CHAR *string1 = NULL;
  CHAR *string2 = NULL;
  CHAR *string3 = NULL;
  INT4 someint;
  REAL8 float1, float2;
  BOOLEAN wasRead = 0;

  // ---------- TEST 1: make sure the XLAL methods can read config files without sections ----------
  const char *cfgname = TEST_DATA_DIR "ConfigFileSample.cfg";

  XLAL_CHECK ( XLALParseDataFile ( &cfgdata, cfgname ) == XLAL_SUCCESS, XLAL_EFUNC, "Failed to parse config-file '%s'", cfgname );

  // check section-listing function
  LALStringVector *sections;
  XLAL_CHECK ( (sections = XLALListConfigFileSections ( cfgdata )) != NULL, XLAL_EFUNC );
  XLAL_CHECK  ( sections->length == 1, XLAL_EFAILED, "%s: Found %d sections instead of 1\n", cfgname, sections->length );
  XLAL_CHECK ( strcmp ( sections->data[0], "default" ) == 0, XLAL_EFAILED, "%s: found section '%s' instead of 'default'\n", cfgname, sections->data[0] );
  XLALDestroyStringVector ( sections );

  // check config-variable reading
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &float1, cfgdata, 0, "float1", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &float2, cfgdata, 0, "float2", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &string1, cfgdata, 0, "string1", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigINT4Variable ( &someint, cfgdata, 0, "int1", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &string2, cfgdata, 0, "string2", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &string3, cfgdata, 0, "string3", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigBOOLEANVariable ( &testBool, cfgdata, 0, "testBool", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  LIGOTimeGPS epochGPS, epochMJDTT;
  XLAL_CHECK ( XLALReadConfigEPOCHVariable ( &epochGPS, cfgdata, 0, "epochGPS", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigEPOCHVariable ( &epochMJDTT, cfgdata, 0, "epochMJDTT", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  REAL8 longHMS, longRad;
  XLAL_CHECK ( XLALReadConfigRAJVariable ( &longHMS, cfgdata, 0, "longHMS", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigRAJVariable ( &longRad, cfgdata, 0, "longRad", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  REAL8 latDMS, latRad;
  XLAL_CHECK ( XLALReadConfigDECJVariable ( &latDMS, cfgdata, 0, "latDMS", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigDECJVariable ( &latRad, cfgdata, 0, "latRad", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  INT8 longInt;
  XLAL_CHECK ( XLALReadConfigINT8Variable ( &longInt, cfgdata, 0, "longInt", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  UINT4Vector *unread = XLALConfigFileGetUnreadEntries ( cfgdata );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALConfigFileGetUnreadEntries() failed\n");
  XLAL_CHECK ( unread == NULL, XLAL_EFAILED, "Some entries in config-file '%s' have not been parsed!\n", cfgname );

  XLALDestroyParsedDataFile (cfgdata);
  cfgdata = NULL;

  // ----- now check the stuff got read-in correctly
  XLAL_CHECK ( float1 == 1.0, XLAL_EFAILED, "%s: got float1 = %g, expected 1.0\n", cfgname, float1 );

  XLAL_CHECK ( float2 == 2.0, XLAL_EFAILED, "%s: got float2 = %g, expected 2.0\n", cfgname, float2 );

  XLAL_CHECK ( someint == 4, XLAL_EFAILED, "%s: someint = %d, expected 4\n", cfgname, someint );

  const char *string1_ref = "some text. You can also use line-continuation";
  XLAL_CHECK ( strcmp(string1, string1_ref) == 0, XLAL_EFAILED, "%s: got string1 = '%s', expected '%s'\n", cfgname, string1, string1_ref );

  const char *string2_ref = "this is also possible, and # here does nothing; and neither does semi-colon ";
  XLAL_CHECK ( strcmp ( string2, string2_ref ) == 0, XLAL_EFAILED, "%s: got string2 = '%s', expected '%s'\n", cfgname, string2, string2_ref );

  const char *string3_ref = "how about #quotes AND line-continuation?";
  XLAL_CHECK ( strcmp ( string3, string3_ref ) == 0, XLAL_EFAILED, "%s: got string3 = '%s', expected '%s'\n", cfgname, string3, string3_ref );

  XLAL_CHECK ( testBool == 0, XLAL_EFAILED, "%s: got testBool = %d, expected 0\n", cfgname, testBool );

  XLAL_CHECK ( XLALGPSCmp ( &epochGPS, &epochMJDTT ) == 0, XLAL_EFAILED, "GPS epoch {%d,%d} differs from MJD(TT) epoch {%d,%d}\n",
               epochGPS.gpsSeconds, epochGPS.gpsNanoSeconds, epochMJDTT.gpsSeconds, epochMJDTT.gpsNanoSeconds );

  REAL8 diff, tol = 3e-15;
  XLAL_CHECK ( (diff = fabs(longHMS - longRad)) < tol, XLAL_EFAILED, "longitude(HMS) = %.16g differs from longitude(rad) = %.16g by %g > %g\n", longHMS, longRad, diff, tol );
  XLAL_CHECK ( (diff = fabs(latDMS - latRad)) < tol, XLAL_EFAILED, "latitude(HMS) = %.16g differs from latitude(rad) = %.16g by %g > %g\n", latDMS, latRad, diff, tol );

  XLAL_CHECK ( longInt == 4294967294, XLAL_EFAILED, "Failed to read an INT8: longInt = %" LAL_INT8_FORMAT " != 4294967294", longInt );

  XLALFree (string1);
  XLALFree (string2);
  XLALFree (string3);
  string1 = string2 = string3 = NULL;

  // ---------- TEST 2: read some values from different sections ----------

  cfgname = TEST_DATA_DIR "ConfigFileSample2.cfg";
  XLAL_CHECK ( XLALParseDataFile ( &cfgdata, cfgname ) == XLAL_SUCCESS, XLAL_EFUNC, "Failed to parse config-file '%s'", cfgname );

  // check section-listing function
  XLAL_CHECK ( (sections = XLALListConfigFileSections ( cfgdata )) != NULL, XLAL_EFUNC );
  XLAL_CHECK  ( sections->length == 3, XLAL_EFAILED, "%s: found %d sections instead of 3\n", cfgname, sections->length );
  XLAL_CHECK ( strcmp ( sections->data[0], "section1" ) == 0, XLAL_EFAILED, "%s: 1st section found is '%s' instead of 'section1'\n", cfgname, sections->data[0] );
  XLAL_CHECK ( strcmp ( sections->data[1], "section2" ) == 0, XLAL_EFAILED, "%s: 2nd section found is '%s' instead of 'section2'\n", cfgname, sections->data[1] );
  XLAL_CHECK ( strcmp ( sections->data[2], "section3" ) == 0, XLAL_EFAILED, "%s: 3rd section found is '%s' instead of 'section3'\n", cfgname, sections->data[2] );
  XLALDestroyStringVector ( sections );

  // Check for a section we have and one we don't
  XLAL_CHECK ( XLALConfigSectionExists ( cfgdata, "section1" ), XLAL_EFAILED );
  XLAL_CHECK ( !XLALConfigSectionExists ( cfgdata, "section5" ), XLAL_EFAILED );

  // check config-variable reading
  XLAL_CHECK ( XLALReadConfigREAL8Variable  (&float1, cfgdata, "section1", "float1", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGVariable (&string1,   cfgdata, "section1", "string1", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigINT4Variable   (&someint,   cfgdata, "section1", "int1", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigSTRINGVariable(&string2,   cfgdata, "section2", "string2", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGVariable(&string3,   cfgdata, "section3", "string3", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigBOOLEANVariable   (&testBool,  cfgdata, "section3", "testBool", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );


  XLAL_CHECK ( XLALReadConfigEPOCHVariable ( &epochGPS, cfgdata, "section2", "epochGPS", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigEPOCHVariable ( &epochMJDTT, cfgdata, "section3", "epochMJDTT", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  unread = XLALConfigFileGetUnreadEntries ( cfgdata );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALConfigFileGetUnreadEntries() failed\n");
  XLAL_CHECK ( unread == NULL, XLAL_EFAILED, "Some entries in config-file '%s' have not been parsed!\n", cfgname );

  XLALDestroyParsedDataFile (cfgdata);
  cfgdata = NULL;

  // ----- now check the stuff got read-in correctly
  XLAL_CHECK ( float1 == 1.0, XLAL_EFAILED, "%s: got float1 = %g, expected 1.0\n", cfgname, float1 );

  XLAL_CHECK ( float2 == 2.0, XLAL_EFAILED, "%s: got float2 = %g, expected 2.0\n", cfgname, float2 );

  XLAL_CHECK ( someint == 4, XLAL_EFAILED, "%s: someint = %d, expected 4\n", cfgname, someint );

  XLAL_CHECK ( strcmp(string1, string1_ref) == 0, XLAL_EFAILED, "%s: got string1 = '%s', expected '%s'\n", cfgname, string1, string1_ref );

  XLAL_CHECK ( strcmp ( string2, string2_ref ) == 0, XLAL_EFAILED, "%s: got string2 = '%s', expected '%s'\n", cfgname, string2, string2_ref );

  XLAL_CHECK ( strcmp ( string3, string3_ref ) == 0, XLAL_EFAILED, "%s: got string3 = '%s', expected '%s'\n", cfgname, string3, string3_ref );

  XLAL_CHECK ( testBool == 0, XLAL_EFAILED, "%s: got testBool = %d, expected 0\n", cfgname, testBool );

  XLAL_CHECK ( XLALGPSCmp ( &epochGPS, &epochMJDTT ) == 0, XLAL_EFAILED, "GPS epoch {%d,%d} differs from MJD(TT) epoch {%d,%d}\n",
               epochGPS.gpsSeconds, epochGPS.gpsNanoSeconds, epochMJDTT.gpsSeconds, epochMJDTT.gpsNanoSeconds );

  // ---------- TEST 3: check reading of compressed files ----------
  cfgname = TEST_DATA_DIR "ConfigFileSample3.cfg.gz";
  XLAL_CHECK ( XLALParseDataFile ( &cfgdata, cfgname ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( cfgdata->lines->nTokens == 4, XLAL_EFAILED );
  XLAL_CHECK ( strcmp(cfgdata->lines->tokens[0], "a") == 0, XLAL_EFAILED );
  XLAL_CHECK ( strcmp(cfgdata->lines->tokens[1], "b") == 0, XLAL_EFAILED );
  XLAL_CHECK ( strcmp(cfgdata->lines->tokens[2], "c") == 0, XLAL_EFAILED );
  XLAL_CHECK ( strcmp(cfgdata->lines->tokens[3], "d") == 0, XLAL_EFAILED );
  XLALDestroyParsedDataFile (cfgdata);
  cfgdata = NULL;

  XLALFree (string1);
  XLALFree (string2);
  XLALFree (string3);

  // -----
  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()


