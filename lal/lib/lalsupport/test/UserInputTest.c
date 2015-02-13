/*
* Copyright (C) 2010 Reinhard Prix
*  Copyright (C) 2007 Jolien Creighton, Reinhard Prix
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

#include <lal/Date.h>
#include <lal/UserInput.h>

/* Error codes and messages */

#define USERINPUTTESTC_ENORM 		0
#define USERINPUTTESTC_MSGENORM 	"Normal exit"


/* Default parameters. */

typedef struct
{
  REAL8 argNum;
  CHAR * argStr;
  BOOLEAN argBool;
  INT4 argInt;
  BOOLEAN argB2;
  CHAR *string2;	// will be read from config-file
  INT4 dummy;
  LIGOTimeGPS epochGPS;
  LIGOTimeGPS epochMJDTT;
  REAL8 longHMS;
  REAL8 longRad;
  REAL8 latDMS;
  REAL8 latRad;
} UserInput_t;

/**
 * some basic consistency checks of the (XLAL) UserInput module, far from exhaustive,
 * but should be enough to catch big obvious malfunctions
 */
int
main(int argc, char *argv[])
{
  int i, my_argc = 8;
#define CFG_FNAME TEST_DATA_DIR "ConfigFileSample.cfg"
  char **my_argv;
  const char *argv_in[] = { "progname", "--argNum=1", "--argStr=xyz", "--argBool=true", "-a", "1", "-b", "@" CFG_FNAME };
  UserInput_t XLAL_INIT_DECL(my_uvars);

  if ( argc > 1 ) {
    XLALPrintError ("No input arguments allowed.\n");
    XLAL_ERROR ( XLAL_EINVAL );
  }

  my_argv = XLALCalloc ( my_argc, sizeof(char*) );
  for (i=0; i < my_argc; i ++ )
    {
      my_argv[i] = XLALCalloc ( 1, strlen(argv_in[i])+1);
      strcpy ( my_argv[i], argv_in[i] );
    }

  /* ---------- Register all test user-variables ---------- */
  UserInput_t *uvar = &my_uvars;
  if ( XLALregREALUserStruct( argNum, 0, UVAR_REQUIRED, "Testing float argument") != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALregREALUserStruct() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  if ( XLALregSTRINGUserStruct( argStr, 0, UVAR_REQUIRED, "Testing string argument") != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALregSTRINGUserStruct() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  if ( XLALregBOOLUserStruct( argBool, 0, UVAR_REQUIRED, "Testing bool argument") != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALregBOOLUserStruct() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  if ( XLALregINTUserStruct( argInt, 'a', UVAR_REQUIRED, "Testing INT argument") != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALregINTUserStruct() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  if ( XLALregINTUserStruct( dummy,  'c', UVAR_OPTIONAL, "Testing INT argument") != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALregINTUserStruct() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  if ( XLALregBOOLUserStruct( argB2, 'b', UVAR_REQUIRED, "Testing short-option bool argument") != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALregBOOLUserStruct() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  if ( XLALregSTRINGUserStruct( string2, 0, UVAR_REQUIRED, "Testing another string argument") != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALregSTRINGUserStruct() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  if ( XLALregEPOCHUserStruct( epochGPS, 0, UVAR_REQUIRED, "Testing epoch given as GPS time") != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALregEPOCHUserStruct() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  if ( XLALregEPOCHUserStruct( epochMJDTT, 0, UVAR_REQUIRED, "Testing epoch given as MJD(TT) time") != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALregEPOCHUserStruct() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  XLAL_CHECK ( XLALregLONGITUDEUserStruct( longHMS, 0, UVAR_REQUIRED, "Testing LONGITUDE(HMS) argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALregLONGITUDEUserStruct( longRad, 0, UVAR_REQUIRED, "Testing LONGITUDE(rad) argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALregLATITUDEUserStruct( latDMS, 0, UVAR_REQUIRED, "Testing LATITUDE(DMS) argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALregLATITUDEUserStruct( latRad, 0, UVAR_REQUIRED, "Testing LATITUDE(rad) argument") == XLAL_SUCCESS, XLAL_EFUNC );

  /* ---------- now read all input from commandline and config-file ---------- */
  if ( XLALUserVarReadAllInput ( my_argc, my_argv ) != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALUserVarReadAllInput() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* ---------- test help-string generation */
  CHAR *helpstr;
  if ( (helpstr = XLALUserVarHelpString ( argv[0])) == NULL ) {
    XLALPrintError ("%s: XLALUserVarHelpString() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  XLALFree ( helpstr );

  /* ---------- test log-generation */
  CHAR *logstr;
  if ( ( logstr = XLALUserVarGetLog (   UVAR_LOGFMT_CFGFILE )) == NULL ) {
    XLALPrintError ("%s: XLALUserVarGetLog() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }
  XLALFree ( logstr );


  /* ---------- test values were read in correctly ---------- */
  if ( uvar->argNum != 1 ) {
    XLALPrintError ("%s: Failed to read in argNum\n", __func__);
    XLAL_ERROR ( XLAL_EFAILED );
  }
  if ( strcmp ( uvar->argStr, "xyz" ) ) {
    XLALPrintError ("%s: Failed to read in argStr\n", __func__);
    XLAL_ERROR ( XLAL_EFAILED );
  }
  if ( !uvar->argBool ) {
    XLALPrintError ("%s: Failed to read in argBool\n", __func__);
    XLAL_ERROR ( XLAL_EFAILED );
  }
  if ( uvar->argInt != 1 ) {
    XLALPrintError ("%s: Failed to read in argInt\n", __func__);
    XLAL_ERROR ( XLAL_EFAILED );
  }
  if ( !uvar->argB2 ) {
    XLALPrintError ("%s: Failed to read in argB2\n", __func__);
    XLAL_ERROR ( XLAL_EFAILED );
  }
  if ( strcmp ( uvar->string2, "this is also possible, and # here does nothing; and neither does semi-colon " ) ) {
    XLALPrintError ("%s: Failed to read in string2\n", __func__);
    XLAL_ERROR ( XLAL_EFAILED );
  }

  XLAL_CHECK ( XLALGPSCmp ( &uvar->epochGPS, &uvar->epochMJDTT ) == 0, XLAL_EFAILED, "GPS epoch {%d,%d} differs from MJD(TT) epoch {%d,%d}\n",
               uvar->epochGPS.gpsSeconds, uvar->epochGPS.gpsNanoSeconds, uvar->epochMJDTT.gpsSeconds, uvar->epochMJDTT.gpsNanoSeconds );

  REAL8 diff, tol = 3e-15;
  XLAL_CHECK ( (diff = fabs(uvar->longHMS - uvar->longRad)) < tol, XLAL_EFAILED, "longitude(HMS) = %.16g differs from longitude(rad) = %.16g by %g > tolerance\n", uvar->longHMS, uvar->longRad, diff, tol );
  XLAL_CHECK ( (diff = fabs(uvar->latDMS - uvar->latRad)) < tol, XLAL_EFAILED, "latitude(HMS) = %.16g differs from latitude(rad) = %.16g by %g > tolerance\n", uvar->latDMS, uvar->latRad, diff, tol );

  /* ----- cleanup ---------- */
  XLALDestroyUserVars();
  for (i=0; i < my_argc; i ++ )
    XLALFree ( my_argv[i] );
  XLALFree ( my_argv );


  LALCheckMemoryLeaks();

  return (USERINPUTTESTC_ENORM);
}
