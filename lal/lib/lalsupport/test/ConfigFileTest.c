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

#include <lal/AVFactories.h>
#include <lal/ConfigFile.h>

/* Default parameters. */

int main ( int argc, char *argv[])
{
  XLAL_CHECK ( argc == 1, XLAL_EINVAL, "Command-line arguments useless here \n");
  XLAL_CHECK ( argv != NULL, XLAL_EINVAL );

  LALParsedDataFile *cfgdata = NULL;

  BOOLEAN testBool;
  CHAR *string1 = NULL;
  CHARVector *string2 = NULL;
  CHAR *string2b = NULL;
  CHAR *string3 = NULL;
  INT4 someint;
  REAL8 float1, float2;
  BOOLEAN wasRead = 0;

  // ---------- TEST 1: make sure the XLAL methods can read config files without sections ----------
  const char *cfgname = TEST_DATA_DIR "ConfigFileSample.cfg";

  XLAL_CHECK ( XLALParseDataFile ( &cfgdata, cfgname ) == XLAL_SUCCESS, XLAL_EFUNC, "Failed to parse config-file '%s'", cfgname );

  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &float1, cfgdata, 0, "float1", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigREAL8Variable ( &float2, cfgdata, 0, "float2", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &string1, cfgdata, 0, "string1", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigINT4Variable ( &someint, cfgdata, 0, "int1", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  UINT4 trunc_len = 35;
  XLAL_CHECK ( ( string2 = XLALCreateCHARVector (trunc_len) ) != NULL, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGNVariable ( string2, cfgdata, 0, "string2", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &string2b, cfgdata, 0, "string2", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGVariable ( &string3, cfgdata, 0, "string3", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigBOOLVariable ( &testBool, cfgdata, 0, "testBool", &wasRead ) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALCheckConfigReadComplete (cfgdata, CONFIGFILE_ERROR) == XLAL_SUCCESS, XLAL_EFUNC );

  XLALDestroyParsedDataFile (cfgdata);
  cfgdata = NULL;

  // ----- now check the stuff got read-in correctly
  XLAL_CHECK ( float1 == 1.0, XLAL_EFAILED, "%s: got float1 = %g, expected 1.0\n", cfgname, float1 );

  XLAL_CHECK ( float2 == 2.0, XLAL_EFAILED, "%s: got float2 = %g, expected 2.0\n", cfgname, float2 );

  XLAL_CHECK ( someint == 4, XLAL_EFAILED, "%s: someint = %d, expected 4\n", cfgname, someint );

  const char *string1_ref = "some text. You can also use line-continuation";
  XLAL_CHECK ( strcmp(string1, string1_ref) == 0, XLAL_EFAILED, "%s: got string1 = '%s', expected '%s'\n", cfgname, string1, string1_ref );

  const char *string2_ref = "this is also possible, and # here does nothing; and neither does semi-colon ";
  XLAL_CHECK ( string2->length == trunc_len, XLAL_EFAILED, "%s: got string2->length = %d, expected %d\n", cfgname, string2->length, trunc_len );
  XLAL_CHECK ( strncmp ( string2->data, string2_ref, trunc_len-1 ) == 0, XLAL_EFAILED, "%s: got string2[%d] = '%s', expected '%.34s'\n", cfgname, trunc_len, string2->data, string2_ref );
  XLAL_CHECK ( strcmp ( string2b, string2_ref ) == 0, XLAL_EFAILED, "%s: got string2 = '%s', expected '%s'\n", cfgname, string2b, string2_ref );

  const char *string3_ref = "how about #quotes AND line-continuation?";
  XLAL_CHECK ( strcmp ( string3, string3_ref ) == 0, XLAL_EFAILED, "%s: got string3 = '%s', expected '%s'\n", cfgname, string3, string3_ref );

  XLAL_CHECK ( testBool == 0, XLAL_EFAILED, "%s: got testBool = %d, expected 0\n", cfgname, testBool );

  XLALFree (string1);
  XLALFree (string2b);
  XLALFree (string3);
  string1 = string2b = string3 = NULL;

  // ---------- TEST 2: read some values from different sections ----------

  cfgname = TEST_DATA_DIR "ConfigFileSample2.cfg";
  XLAL_CHECK ( XLALParseDataFile ( &cfgdata, cfgname ) == XLAL_SUCCESS, XLAL_EFUNC, "Failed to parse config-file '%s'", cfgname );

  /* Check for a section we have and one we don't */
  XLAL_CHECK ( XLALConfigSectionExists ( cfgdata, "section1" ), XLAL_EFAILED );
  XLAL_CHECK ( !XLALConfigSectionExists ( cfgdata, "section5" ), XLAL_EFAILED );


  XLAL_CHECK ( XLALReadConfigREAL8Variable  (&float1, cfgdata, "section1", "float1", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGVariable (&string1,   cfgdata, "section1", "string1", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigINT4Variable   (&someint,   cfgdata, "section1", "int1", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigSTRINGNVariable(string2,   cfgdata, "section2", "string2", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigSTRINGVariable(&string2b,   cfgdata, "section2", "string2", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALReadConfigSTRINGVariable(&string3,   cfgdata, "section3", "string3", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALReadConfigBOOLVariable   (&testBool,  cfgdata, "section3", "testBool", &wasRead) == XLAL_SUCCESS, XLAL_EFUNC );

  XLAL_CHECK ( XLALCheckConfigReadComplete (cfgdata, CONFIGFILE_ERROR) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALDestroyParsedDataFile (cfgdata);
  cfgdata = NULL;

  // ----- now check the stuff got read-in correctly
  XLAL_CHECK ( float1 == 1.0, XLAL_EFAILED, "%s: got float1 = %g, expected 1.0\n", cfgname, float1 );

  XLAL_CHECK ( float2 == 2.0, XLAL_EFAILED, "%s: got float2 = %g, expected 2.0\n", cfgname, float2 );

  XLAL_CHECK ( someint == 4, XLAL_EFAILED, "%s: someint = %d, expected 4\n", cfgname, someint );

  XLAL_CHECK ( strcmp(string1, string1_ref) == 0, XLAL_EFAILED, "%s: got string1 = '%s', expected '%s'\n", cfgname, string1, string1_ref );

  XLAL_CHECK ( string2->length == trunc_len, XLAL_EFAILED, "%s: got string2->length = %d, expected %d\n", cfgname, string2->length, trunc_len );
  XLAL_CHECK ( strncmp ( string2->data, string2_ref, trunc_len-1 ) == 0, XLAL_EFAILED, "%s: got string2[%d] = '%s', expected '%.34s'\n", cfgname, trunc_len, string2->data, string2_ref );
  XLAL_CHECK ( strcmp ( string2b, string2_ref ) == 0, XLAL_EFAILED, "%s: got string2 = '%s', expected '%s'\n", cfgname, string2b, string2_ref );

  XLAL_CHECK ( strcmp ( string3, string3_ref ) == 0, XLAL_EFAILED, "%s: got string3 = '%s', expected '%s'\n", cfgname, string3, string3_ref );

  XLAL_CHECK ( testBool == 0, XLAL_EFAILED, "%s: got testBool = %d, expected 0\n", cfgname, testBool );

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
  XLALDestroyCHARVector (string2);
  XLALFree (string2b);
  XLALFree (string3);

  // -----
  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()


