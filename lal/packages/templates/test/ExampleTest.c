/*
*  Copyright (C) 2007 Bernd Machenschalk, Jolien Creighton
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

/**** <lalVerbatim file="ExampleTestCV">
 * Author: Al A. Lal
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Program \texttt{ExampleTest.c}}
 *
 * %% One-line description of test program.
 *
 * \subsubsection*{Usage}
 *
 * %% Command syntax and options.
 * %% The default mode, with no options specified, is always the
 * %% test that is performed using make check.
 * %%
 * %%   \begin{verbatim}
 * %%   LALExampleTest [options]
 * %%   Options:
 * %%     -d debuglevel  set lalDebugLevel to debuglevel
 * %%     -v             verbose: print lots of output
 * %%     -h             help: print the help message
 * %%   \end{verbatim}
 *
 * \subsubsection*{Description}
 *
 * %% A description of what the test program does.
 *
 * \subsubsection*{Exit codes}
 *
 * %% A table of the possible exit codes and their meanings.
 * %% The program should return exit code 0 on successful completion;
 * %% otherwise it should return some non-zero value.
 * %%
 * %%   \begin{tabular}{|c|l|}\hline Code & Explanation \\
 * %%     \tt 0 & Success, normal exit.       \\
 * %%     \tt 1 & Failure: invalid arguments. \\
 * %%     \tt 2 & Failure: wrong return code. \\
 * %%     \tt 3 & Failure: subroutine failed. \\
 * %%     \tt 4 & Failure: incorrect results. \\
 * %%   \hline\end{tabular}
 *
 * \subsubsection*{Uses}
 *
 * %% List of any external functions called by this function.
 * \begin{verbatim}
 * \end{verbatim}
 *
 * \subsubsection*{Notes}
 *
 * %% Any relevant notes.
 *
 * \vfill{\footnotesize\input{ExampleTestCV}}
 **** </lalLaTeX> */

/**
 ** INCLUDE STANDARD LIBRARY HEADERS
 ** note LALStdLib.h already includes stdio.h and stdarg.h
 **
 **   #include <stdlib.h>
 **   #include <string.h>
 **/

/** INCLUDE ANY LAL HEADERS **/
#include <lal/LALStdlib.h>
#include <lal/Example.h>

NRCSID (EXAMPLETESTC,"$Id$");

/**
 ** DEFINE CONSTANTS, GLOBAL VARIABLES, AND MACROS
 **
 ** useful macros for return codes:
 **
 **   #define FAIL_ARGS 1
 **   #define FAIL_CODE 2
 **   #define FAIL_ESUB 3
 **   #define FAIL_RSLT 4
 **
 **
 ** global variables:
 **
 **   const char *program;
 **   const char *usage = "Usage: %s [options]\nOptions:\n"
 **       "  -d debuglevel  set lalDebugLevel to debuglevel\n"
 **       "  -v             verbose: print lots of output\n"
 **       "  -h             help: print this message\n";
 **   int verbose = 0;
 **/

/**
 ** DECLARE AND SET GLOBAL DEBUG LEVEL
 **
 ** see the section (currently 7.4.1) of the LSD on "Status-reporting objects"
 ** for a list of predefined debug levels
 **/
int lalDebugLevel = LALMSGLVL3;


/** THE PROGRAM **/
int main( int argc, char *argv[] )
{
  /**
   ** variable declarations
   **
   ** the status must be initially blank; this can be done by making it static:
   **
   **   static LALStatus status;
   **
   ** other variables:
   **
   **   ExampleOutput output;
   **   ExampleInput  input;
   **   ExampleParams params;
   **/

  /**
   ** parse arguments, if desired
   **
   ** for example:
   **
   **   program = *argv;
   **   while ( --argc > 0 )
   **   {
   **     ++argv;
   **     if ( ! strcmp( *argv, "-d" ) )
   **     {
   **       --argc;
   **       ++argv;
   **       lalDebugLevel = atoi( *argv );
   **       continue;
   **     }
   **     if ( ! strcmp( *argv, "-v" ) )
   **     {
   **       verbose = 1;
   **       continue;
   **     }
   **     if ( ! strcmp( *argv, "-h" ) )
   **     {
   **       fprintf( stderr, usage, program );
   **       return 0;
   **     }
   **     fprintf( stderr, "no such option %s\n", *argv );
   **     fprintf( stderr, usage, program );
   **     return FAIL_ARGS;
   **   }

  /**
   ** test response to invalid data
   **
   ** these tests must be wrapped so they are not done when debugging is
   ** disabled:
   **
   **   #ifndef LAL_NDEBUG
   **     if ( ! lalNoDebug )
   **     {
   **       LALExample( &status, NULL, &input, &params );
   **       if ( status.statusCode != EXAMPLEH_ENULL
   **            || strcmp( status.statusDescription, EXAMPLEH_MSGENULL ) )
   **       {
   **         fprintf( stderr, "incorrect error code %d and message %s\n",
   **             status.statusCode, status.statusDescription );
   **         fprintf( stderr, "expecting error code %d and message %s\n",
   **             EXAMPLEH_ENULL, EXAMPLEH_MSGENULL );
   **         return FAIL_CODE;
   **       }
   **
   **       fputs( "PASS: Test response to invalid data\n", stderr );
   **     }
   **   #endif
   **/

  /**
   ** test response to valid data
   **
   ** for example:
   **
   **   LALExample( &status, &output, &input, &param );
   **   if ( status.statusCode )
   **   {
   **     fprintf( stderr, "received error code %d and message %s\n",
   **         status.statusCode, status.statusDescription );
   **     return FAIL_ESUB;
   **   }
   **
   ** (perform checks on the contents of output too!)
   **
   **
   ** during default operation, output to the screen should be minimal:
   **
   **   for ( i = 0; i < 1048576; ++i ) printf( "%d\n", i );  !!! BAD !!!
   **
   ** but it is OK to indicate that tests have passed:
   **
   **   fputs( "PASS: Test response to valid data\n", stderr );
   **/

  /**
   ** check for memory leaks and return success:
   **
   **   LALCheckMemoryLeaks();
   **   fputs( "PASS: All tests\n" );
   **/
  return 0;
}
