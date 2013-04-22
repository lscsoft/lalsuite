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

#include <lal/UserInput.h>

/* Error codes and messages */

#define USERINPUTTESTC_ENORM 		0
#define USERINPUTTESTC_MSGENORM 	"Normal exit"


/* Default parameters. */
INT4 lalDebugLevel=3;

typedef struct
{
  REAL8 argNum;
  CHAR * argStr;
  BOOLEAN argBool;
  INT4 argInt;
  BOOLEAN argB2;
  CHAR *string2;	// will be read from config-file
  INT4 dummy;
} UserInput_t;

UserInput_t empty_UserInput_t;
#define TESTSTRING "this is also possible, and # here does nothing "

const char *cfgfile_content = \
"## Some 'tough' tests for config-file reading routines\n"
"# comment line\n"
"float1 = 1.0    ; ## semi-colon ignored\n"
"\n"
"string1 = some text.\\\n"
"	You can also use\\\n"
"	line-continuation\n"
"\n"
"   int1 = 4      # whatever that means\n"
"\n"
"# Comment before section\n"
"# Comment after section\n"
"string2 = \"" TESTSTRING "\"	# but this is a comment\n"
"\n"
"string3 = \"how about #quotes\\\n"
"	AND line-continuation?\"		# yet another comment\n"
"testBool = False	# use yes/no/0/1/true/false, case INsensitive\n"
"# etc etc.\n"
;



/** some basic consistency checks of the (XLAL) UserInput module, far from exhaustive,
 * but should be enough to catch big obvious malfunctions
 */
int
main(int argc, char *argv[])
{
  int i, my_argc = 8;
  #define CFG_FNAME "ConfigFile.cfg"
  char **my_argv;
  const char *argv_in[] = { "progname", "--argNum=1", "--argStr=xyz", "--argBool=true", "-a", "1", "-b", "@" CFG_FNAME };
  UserInput_t my_uvars = empty_UserInput_t;

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

  /* laldebug level always needs to be read first (before any lal-mallocs!) */
  if ( XLALGetDebugLevel (my_argc, my_argv, 'v') != XLAL_SUCCESS ) {
    XLALPrintError ("%s: XLALGetDebugLevel() failed with code %d\n", __func__, xlalErrno );
    XLAL_ERROR ( XLAL_EFUNC );
  }

  /* ----- dump config-file content into config-file ----- */
  FILE *fid;
  if ( (fid = fopen ( CFG_FNAME, "wb" )) == NULL ) {
    XLALPrintError ("%s: Failed to open configfile '%s' for writing.\n", __func__, CFG_FNAME );
    XLAL_ERROR ( XLAL_ESYS );
  }
  fprintf ( fid, "%s\n", cfgfile_content );
  fclose(fid);

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
  if ( strcmp ( uvar->string2, TESTSTRING ) ) {
    XLALPrintError ("%s: Failed to read in string2\n", __func__);
    XLAL_ERROR ( XLAL_EFAILED );
  }

  /* ----- cleanup ---------- */
  XLALDestroyUserVars();
  for (i=0; i < my_argc; i ++ )
    XLALFree ( my_argv[i] );
  XLALFree ( my_argv );


  LALCheckMemoryLeaks();

  return (USERINPUTTESTC_ENORM);
}
