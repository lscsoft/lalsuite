/*
*  Copyright (C) 2007 Reinhard Prix
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

/**
   \file
   \ingroup ConfigFile_h
   \author Reinhard Prix

   \brief Tests the routines in \ref ConfigFile.h.

   \heading{Usage}
   \code
ConfigFileTest
   \endcode

\heading{Description}

Does some standard-tests for the config-file reading routines. No
extensive error-condition checking is done here, we only check if the
basic functionality works.

*/

/* Error codes and messages */

/**\name Error Codes */ /*@{*/
#define CONFIGFILETESTC_ENORM           0       /**< Normal exit */
#define CONFIGFILETESTC_EFLOAT          1       /**< Read-in REAL8 variable is not what it should be... */
#define CONFIGFILETESTC_EINT            2       /**< Read-in INT4 variable is not what it should be... */
#define CONFIGFILETESTC_EBOOL           3       /**< Read-in BOOL variable is not what it should be... */
#define CONFIGFILETESTC_ESTRING         4       /**< Read-in STRING-variable is not what it should be... */
#define CONFIGFILETESTC_ESUB            5       /**< Error occurred in sub-routine */
/*@}*/

/** \cond DONT_DOXYGEN */
#define CONFIGFILETESTC_MSGENORM 	"Normal exit"
#define CONFIGFILETESTC_MSGEFLOAT 	"Read-in REAL8 variable is not what it should be..."
#define CONFIGFILETESTC_MSGEINT 	"Read-in INT4 variable is not what it should be..."
#define CONFIGFILETESTC_MSGEBOOL 	"Read-in BOOL variable is not what it should be..."
#define CONFIGFILETESTC_MSGESTRING 	"Read-in STRING-variable is not what it should be..."
#define CONFIGFILETESTC_MSGESUB	 	"Error occurred in sub-routine"


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
                   "", (msg) );                                      \
} while (0)

#define INFO( statement )                                            \
do {                                                                 \
  if ( lalDebugLevel & LALINFO )                                     \
    XLALPrintError( "Info[0]: program %s, file %s, line %d, %s\n"     \
                   "        %s\n", *argv, __FILE__, __LINE__,        \
              "$Id$", (statement) );                     \
} while (0)

#define SUB( func, statusptr )                                       \
do {                                                                 \
  if ( (func), (statusptr)->statusCode ) {                           \
    ERROR( CONFIGFILETESTC_ESUB, CONFIGFILETESTC_MSGESUB,      \
           "Function call \"" #func "\" failed:" );                  \
    return CONFIGFILETESTC_ESUB;                                  \
  }                                                                  \
} while (0)
/******************************************************************/

#define TRUE (1==1)
#define FALSE (1==0)

int main(int argc, char *argv[]){
  static LALStatus       status;
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

  SUB (LALParseDataFile (&status, &cfgdata, DATADIR "ConfigFileSample.cfg"), &status);

  SUB (LALReadConfigREAL8Variable  (&status, &somefloat, cfgdata, "float1", &wasRead), &status);
  SUB (LALReadConfigSTRINGVariable (&status, &string1,   cfgdata, "string1", &wasRead), &status);

  SUB (LALReadConfigINT4Variable   (&status, &someint,   cfgdata, "int1", &wasRead), &status);

  SUB (LALCHARCreateVector (&status, &string2, 35), &status);
  SUB (LALReadConfigSTRINGNVariable(&status, string2,   cfgdata, "string2", &wasRead), &status);

  SUB (LALReadConfigSTRINGVariable(&status, &string2b,   cfgdata, "string2", &wasRead), &status);
  SUB (LALReadConfigSTRINGVariable(&status, &string3,   cfgdata, "string3", &wasRead), &status);

  SUB (LALReadConfigBOOLVariable   (&status, &testBool,  cfgdata, "testBool", &wasRead), &status);

  SUB (LALCheckConfigReadComplete (&status, cfgdata, CONFIGFILE_ERROR), &status);
  SUB (LALDestroyParsedDataFile (&status, &cfgdata), &status);

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

  LALFree (string1);
  LALCHARDestroyVector (&status, &string2);
  LALFree (string2b);
  LALFree (string3);

  LALCheckMemoryLeaks();

  return CONFIGFILETESTC_ENORM;
}
/** \endcond */
