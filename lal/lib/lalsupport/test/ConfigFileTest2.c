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

/* Error codes and messages */

#define CONFIGFILETESTC_ENORM 		0
#define CONFIGFILETESTC_EFLOAT 		1
#define CONFIGFILETESTC_EINT 		2
#define CONFIGFILETESTC_EBOOL 		3
#define CONFIGFILETESTC_ESTRING 	4
#define CONFIGFILETESTC_ESUB	 	5
#define CONFIGFILETESTC_EEXISTS     6

#define CONFIGFILETESTC_MSGENORM 	"Normal exit"
#define CONFIGFILETESTC_MSGEFLOAT 	"Read-in REAL8 variable is not what it should be..."
#define CONFIGFILETESTC_MSGEINT 	"Read-in INT4 variable is not what it should be..."
#define CONFIGFILETESTC_MSGEBOOL 	"Read-in BOOL variable is not what it should be..."
#define CONFIGFILETESTC_MSGESTRING 	"Read-in STRING-variable is not what it should be..."
#define CONFIGFILETESTC_MSGESUB	 	"Error occurred in sub-routine"
#define CONFIGFILETESTC_MSGEXISTS   "Error occurrent in sectionExists"


/* Default parameters. */

INT4 lalDebugLevel=3;


/*********************************************************************/
/* Macros for printing errors & testing subroutines (from Creighton) */
/*********************************************************************/

#define ERROR( code, msg, statement )                                \
do {                                                                 \
  if ( lalDebugLevel & LALERROR )                                    \
    XLALPrintError( "Error[0] %d: program %s, file %s, line %d, %s\n" \
                    "        %s %s\n", (code), *argv, __FILE__,       \
                    __LINE__, "$Id$", statement ? statement :  \
                    "", (msg) );                                        \
 } while (0)

/******************************************************************/

#define TRUE (1==1)
#define FALSE (1==0)

int main(int argc, char *argv[]){
  static LALParsedDataFile *cfgdata;

  BOOLEAN testBool;
  CHAR *string1 = NULL;
  CHARVector *string2 = NULL;
  CHAR *string2b = NULL;
  CHAR *string3 = NULL;
  INT4 someint;
  REAL8 somefloat;
  BOOLEAN wasRead = FALSE;

  if ( argc > 1 )
    XLALPrintError ("WARNING: commond-line arguments useless here \n");


  /* First, make sure the XLAL methods can still read config files without sections */
  XLALParseDataFile (&cfgdata, DATADIR "ConfigFileSample.cfg");

  XLALReadConfigREAL8Variable  (&somefloat, cfgdata, 0, "float1", &wasRead);
  XLALReadConfigSTRINGVariable (&string1,   cfgdata, 0, "string1", &wasRead);

  XLALReadConfigINT4Variable   (&someint,   cfgdata, 0, "int1", &wasRead);

  string2 = XLALCreateCHARVector (35);
  XLALReadConfigSTRINGNVariable(string2,   cfgdata, 0, "string2", &wasRead);

  XLALReadConfigSTRINGVariable(&string2b,   cfgdata, 0, "string2", &wasRead);
  XLALReadConfigSTRINGVariable(&string3,   cfgdata, 0, "string3", &wasRead);

  XLALReadConfigBOOLVariable   (&testBool,  cfgdata, 0, "testBool", &wasRead);

  XLALCheckConfigReadComplete (cfgdata, CONFIGFILE_ERROR);
  XLALDestroyParsedDataFile (cfgdata);
  cfgdata = NULL;

  /* now check the stuff got read-in correctly */
  if (somefloat != 1.0) {
    ERROR (CONFIGFILETESTC_EFLOAT, CONFIGFILETESTC_MSGEFLOAT, 0);
    return (CONFIGFILETESTC_EFLOAT);
  }
  if ( strcmp (string1, "some text. You can also use line-continuation") ) {
    XLALPrintError ("read-in: '%s'\n", string1);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }
  if (someint != 4) {
    ERROR (CONFIGFILETESTC_EINT, CONFIGFILETESTC_MSGEINT, 0);
    return (CONFIGFILETESTC_EINT);
  }
  if ( strcmp(string2->data, "this is also possible, and # here ") ) {
    XLALPrintError ("read-in: '%s'\n", string2->data);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }
  if ( strcmp(string2b, "this is also possible, and # here does nothing ")) {
    XLALPrintError ("read-in: '%s'\n", string2b);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }
  if ( strcmp(string3, "how about #quotes AND line-continuation?") ) {
    XLALPrintError ("read-in: '%s'\n", string3);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }


  if ( testBool != 0 ) {
    ERROR (CONFIGFILETESTC_EBOOL, CONFIGFILETESTC_MSGEBOOL, 0);
    return (CONFIGFILETESTC_EBOOL);
  }

  XLALFree (string1);
  XLALDestroyCHARVector (string2);
  XLALFree (string2b);
  XLALFree (string3);


  /* Now, try to read some values from different sections */
  string1 = NULL;
  string2 = NULL;
  string2b = NULL;
  string3 = NULL;

  XLALParseDataFile (&cfgdata, DATADIR "ConfigFileSample2.cfg");

  /* Check for a section we have and one we don't */
  if (XLALConfigSectionExists(cfgdata, "section1") == 0) {
    ERROR (CONFIGFILETESTC_EEXISTS, CONFIGFILETESTC_MSGEXISTS, 0);
    return (CONFIGFILETESTC_EEXISTS);
  }

  if (XLALConfigSectionExists(cfgdata, "section5") == 1) {
    ERROR (CONFIGFILETESTC_EEXISTS, CONFIGFILETESTC_MSGEXISTS, 0);
    return (CONFIGFILETESTC_EEXISTS);
  }

  XLALReadConfigREAL8Variable  (&somefloat, cfgdata, "section1", "float1", &wasRead);
  XLALReadConfigSTRINGVariable (&string1,   cfgdata, "section1", "string1", &wasRead);

  XLALReadConfigINT4Variable   (&someint,   cfgdata, "section1", "int1", &wasRead);

  string2 = XLALCreateCHARVector ( 35 );

  XLALReadConfigSTRINGNVariable(string2,   cfgdata, "section2", "string2", &wasRead);

  XLALReadConfigSTRINGVariable(&string2b,   cfgdata, "section2", "string2", &wasRead);
  XLALReadConfigSTRINGVariable(&string3,   cfgdata, "section3", "string3", &wasRead);

  XLALReadConfigBOOLVariable   (&testBool,  cfgdata, "section3", "testBool", &wasRead);

  XLALCheckConfigReadComplete (cfgdata, CONFIGFILE_ERROR);
  XLALDestroyParsedDataFile (cfgdata);
  cfgdata = NULL;

  /* now check the stuff got read-in correctly */
  if (somefloat != 1.0) {
    ERROR (CONFIGFILETESTC_EFLOAT, CONFIGFILETESTC_MSGEFLOAT, 0);
    return (CONFIGFILETESTC_EFLOAT);
  }
  if ( strcmp (string1, "some text. You can also use line-continuation") ) {
    XLALPrintError ("read-in: '%s'\n", string1);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }
  if (someint != 4) {
    ERROR (CONFIGFILETESTC_EINT, CONFIGFILETESTC_MSGEINT, 0);
    return (CONFIGFILETESTC_EINT);
  }
  if ( strcmp(string2->data, "this is also possible, and # here ") ) {
    XLALPrintError ("read-in: '%s'\n", string2->data);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }
  if ( strcmp(string2b, "this is also possible, and # here does nothing ")) {
    XLALPrintError ("read-in: '%s'\n", string2b);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }
  if ( strcmp(string3, "how about #quotes AND line-continuation?") ) {
    XLALPrintError ("read-in: '%s'\n", string3);
    ERROR (CONFIGFILETESTC_ESTRING, CONFIGFILETESTC_MSGESTRING, 0);
    return (CONFIGFILETESTC_ESTRING);
  }


  if ( testBool != 0 ) {
    ERROR (CONFIGFILETESTC_EBOOL, CONFIGFILETESTC_MSGEBOOL, 0);
    return (CONFIGFILETESTC_EBOOL);
  }

  /* check reading of compressed files */
  XLALParseDataFile (&cfgdata, DATADIR "ConfigFileSample3.cfg.gz");
  XLAL_CHECK_VAL( EXIT_FAILURE, cfgdata->lines->nTokens == 4, XLAL_EFAILED );
  XLAL_CHECK_VAL( EXIT_FAILURE, strcmp(cfgdata->lines->tokens[0], "a") == 0, XLAL_EFAILED );
  XLAL_CHECK_VAL( EXIT_FAILURE, strcmp(cfgdata->lines->tokens[1], "b") == 0, XLAL_EFAILED );
  XLAL_CHECK_VAL( EXIT_FAILURE, strcmp(cfgdata->lines->tokens[2], "c") == 0, XLAL_EFAILED );
  XLAL_CHECK_VAL( EXIT_FAILURE, strcmp(cfgdata->lines->tokens[3], "d") == 0, XLAL_EFAILED );
  XLALDestroyParsedDataFile (cfgdata);
  cfgdata = NULL;

  XLALFree (string1);
  XLALDestroyCHARVector (string2);
  XLALFree (string2b);
  XLALFree (string3);
  LALCheckMemoryLeaks();

  return CONFIGFILETESTC_ENORM;
}


