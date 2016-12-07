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

#include <lal/LALVCSInfo.h>
#include <lal/Date.h>
#include <lal/LALString.h>

#include <lal/UserInput.h>

// for test print usage and help page
int XLALUserVarPrintUsage ( FILE *file );
int XLALUserVarPrintHelp ( FILE *file );

typedef struct
{
  REAL8 argNum;
  CHAR * argStr;
  BOOLEAN argBool;
  INT4 argInt;
  UINT4 argUInt;
  BOOLEAN argB2;
  CHAR *string2;	// will be read from config-file
  INT4 dummy;
  LIGOTimeGPS epochGPS;
  LIGOTimeGPS epochMJDTT;
  REAL8 longHMS;
  REAL8 longRad;
  REAL8 latDMS;
  REAL8 latRad;
  INT8 longInt;
  int argEnum;
  int argFlag;
  BOOLEAN long_help;
} UserInput_t;

/**
 * some basic consistency checks of the (XLAL) UserInput module, far from exhaustive,
 * but should be enough to catch big obvious malfunctions
 */
int
main(int argc, char *argv[])
{

  const char *argv_in[] = {
    "progname",
    "--argNum=1",
    "--argStr=xyz",
    "--argBool=true",
    "-a", "-1",
    "-A", "7",
    "-b",
    "--argEnum=enum2",
    "--argFlag=flagA,flagC",
    "@" TEST_DATA_DIR "ConfigFileSample.cfg"
  };

  const int my_argc = XLAL_NUM_ELEM( argv_in );
  char **my_argv = XLALCalloc ( my_argc, sizeof(char*) );
  for ( int i=0; i < my_argc; i ++ )
    {
      my_argv[i] = XLALCalloc ( 1, strlen(argv_in[i])+1);
      strcpy ( my_argv[i], argv_in[i] );
    }

  /* ---------- Register all test user-variables ---------- */
  UserInput_t XLAL_INIT_DECL(my_uvars);
  UserInput_t *uvar = &my_uvars;
  uvar->string2 = XLALStringDuplicate ( "this is the default value");
  uvar->argEnum = 0;
  uvar->argFlag = 0;

  const UserChoices enumData = { { -1, "noenum" }, { 1, "enum1" }, { 2, "enum2" }, { 2, "enumB" }, { 0, "enum0" } };
  const UserChoices flagData = { { -1, "noflag" }, { 1, "flagA" }, { 2, "flagB" }, { 4, "flagC" }, { 5, "flagAC" } };

  XLAL_CHECK ( XLALRegisterUvarMember( argNum, REAL8, 0, REQUIRED, "Testing float argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( argStr, STRING, 0, REQUIRED, "Testing string argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( argBool, BOOLEAN, 0, REQUIRED, "Testing bool argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( argInt, INT4, 'a', REQUIRED, "Testing INT4 argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( argUInt, UINT4, 'A', REQUIRED, "Testing UINT4 argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( dummy,  INT4, 'c', OPTIONAL, "Testing INT4 argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( argB2, BOOLEAN, 'b', REQUIRED, "Testing short-option bool argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( string2, STRING, 0, REQUIRED, "Testing another string argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( epochGPS, EPOCH, 0, REQUIRED, "Testing epoch given as GPS time") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( epochMJDTT, EPOCH, 0, REQUIRED, "Testing epoch given as MJD(TT) time") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( longHMS, RAJ, 0, REQUIRED, "Testing RAJ(HMS) argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( longRad, RAJ, 0, REQUIRED, "Testing RAJ(rad) argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( latDMS, DECJ, 0, REQUIRED, "Testing DECJ(DMS) argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( latRad, DECJ, 0, REQUIRED, "Testing DECJ(rad) argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( longInt, INT8, 0, REQUIRED, "Testing INT8 argument") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarAuxDataMember( argEnum, UserEnum, &enumData, 0, REQUIRED, "Testing user enumeration") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarAuxDataMember( argFlag, UserFlag, &flagData, 0, REQUIRED, "Testing user bitflag") == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK ( XLALRegisterUvarMember( long_help, BOOLEAN, 0, NODEFAULT,
                                       "This option is here to test the help page wrapping of long strings. "
                                       "This option is here to test the help page wrapping of long strings. "
                                       "This option is here to test the help page wrapping of long strings. "
                                       "This option is here to test the help page wrapping of long strings. "
                                       "This option is here to test the help page wrapping of long strings. "
                                       "\n"
                                       "This~option~is~here~to~test~the~help~page~wrapping~of~long~strings~without~spaces.~"
                                       "This~option~is~here~to~test~the~help~page~wrapping~of~long~strings~without~spaces.~"
                                       "This~option~is~here~to~test~the~help~page~wrapping~of~long~strings~without~spaces."
                 ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* ---------- now read all input from commandline and config-file ---------- */
  BOOLEAN should_exit = 0;
  if ( argc > 1 ) {
    /* parse actual command-line arguments: this is only useful for debugging and testing, never used by 'make check' */
    XLAL_CHECK ( XLALUserVarReadAllInput ( &should_exit, argc, argv, lalVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( should_exit ) {
      return EXIT_FAILURE;
    }
  } else {
    /* parse constructed command-line arguments: used by 'make check' */
    XLAL_CHECK ( XLALUserVarReadAllInput ( &should_exit, my_argc, my_argv, lalVCSInfoList ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK ( should_exit == 0, XLAL_EFUNC );
  }

  /* ---------- test print usage and help page */
  printf( "=== Begin usage string ===\n" );
  fflush( stdout );
  XLALUserVarPrintUsage( stdout );
  printf( "--- End usage string ---\n" );
  printf( "=== Begin help page ===\n" );
  fflush( stdout );
  XLALUserVarPrintHelp( stdout );
  printf( "--- End help page ---\n" );
  fflush( stdout );

  /* ---------- test log-generation */
  printf( "=== Begin log generation test ===\n" );
  for ( int i = UVAR_LOGFMT_RAWFORM; i < UVAR_LOGFMT_LAST; ++i ) {
    CHAR *logstr;
    XLAL_CHECK ( ( logstr = XLALUserVarGetLog ( i )) != NULL, XLAL_EFUNC );
    printf( "Log format #%i: >>>%s<<<\n", i, logstr );
    XLALFree ( logstr );
  }
  printf( "--- End log generation test ---\n" );

  /* ---------- test values were read in correctly ---------- */
  XLAL_CHECK ( uvar->argNum == 1, XLAL_EFAILED, "Failed to read in argNum\n" );
  XLAL_CHECK ( strcmp ( uvar->argStr, "xyz" ) == 0, XLAL_EFAILED, "Failed to read in argStr\n" );
  XLAL_CHECK ( uvar->argBool, XLAL_EFAILED, "Failed to read in argBool\n" );
  XLAL_CHECK ( uvar->argInt == -1, XLAL_EFAILED, "Failed to read in argInt\n" );
  XLAL_CHECK ( uvar->argUInt == 7, XLAL_EFAILED, "Failed to read in argUInt\n" );
  XLAL_CHECK ( uvar->argB2, XLAL_EFAILED, "Failed to read in argB2\n" );
  XLAL_CHECK ( strcmp ( uvar->string2, "this is also possible, and # here does nothing; and neither does semi-colon " ) == 0, XLAL_EFAILED, "Failed to read in string2\n" );
  XLAL_CHECK ( uvar->argEnum == 2, XLAL_EFAILED, "Failed to read in argEnum\n" );
  XLAL_CHECK ( uvar->argFlag == 5, XLAL_EFAILED, "Failed to read in argFlag\n" );

  char buf1[256], buf2[256];
  XLAL_CHECK ( XLALGPSCmp ( &uvar->epochGPS, &uvar->epochMJDTT ) == 0, XLAL_EFAILED, "GPS epoch %s differs from MJD(TT) epoch %s\n",
               XLALGPSToStr ( buf1, &uvar->epochGPS), XLALGPSToStr ( buf2, &uvar->epochMJDTT ) );

  REAL8 diff, tol = 3e-15;
  XLAL_CHECK ( (diff = fabs(uvar->longHMS - uvar->longRad)) < tol, XLAL_EFAILED, "longitude(HMS) = %.16g differs from longitude(rad) = %.16g by %g > %g\n", uvar->longHMS, uvar->longRad, diff, tol );
  XLAL_CHECK ( (diff = fabs(uvar->latDMS - uvar->latRad)) < tol, XLAL_EFAILED, "latitude(HMS) = %.16g differs from latitude(rad) = %.16g by %g > %g\n", uvar->latDMS, uvar->latRad, diff, tol );

  XLAL_CHECK ( uvar->longInt == 4294967294, XLAL_EFAILED, "Failed to read an INT8: longInt = %" LAL_INT8_FORMAT " != 4294967294", uvar->longInt );

  /* ----- cleanup ---------- */
  XLALDestroyUserVars();
  for ( int i=0; i < my_argc; i ++ ) {
    XLALFree ( my_argv[i] );
  }
  XLALFree ( my_argv );

  LALCheckMemoryLeaks();

  return XLAL_SUCCESS;

} // main()
