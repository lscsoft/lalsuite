dnl $Id$
/*************** <lalVerbatim file="ArrayFactoriesTestCV"> ****
$Id$
**************** </lalVerbatim> *******************************/

/* <lalLaTeX>

\subsection{Program \texttt{ArrayFactoriesTest.c}}
\label{ss:ArrayFactoriesTest.c}

A program to test create/destroy array routines.

\subsubsection*{Usage}
\begin{verbatim}
ArrayFactoriesTest [options]
Options:
  -h         print help
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
\end{verbatim}

\subsubsection*{Description}

\subsubsection*{Exit codes}
\begin{tabular}{|c|l|}
\hline
 Code & Explanation                   \\
\hline
\tt 0 & Success, normal exit.         \\
\tt 1 & Subroutine failed.            \\
\hline
\end{tabular}

\subsubsection*{Algorithm}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
TYPECODECreateArray()
TYPECODEDestroyArray()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{ArrayFactoriesTestCV}}

</lalLaTeX> */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "LALConfig.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include "LALStdlib.h"
#include "AVFactories.h"

#define CODES_(x) #x
#define CODES(x) CODES_(x)

NRCSID( MAIN, "$Id$" );

extern char *optarg;
extern int   optind;

int lalDebugLevel = 0;
int verbose    = 0;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

static void
ClearStatus (LALStatus *status);

define(`TYPECODE',`Z')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`C')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`D')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`S')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`I2')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`I4')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`I8')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`U2')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`U4')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`U8')
include(`ArrayFactoriesTestFunction.m4')

define(`TYPECODE',`')
include(`ArrayFactoriesTestFunction.m4')

int main( int argc, char *argv[] )
{
  ParseOptions( argc, argv );

  ArrayFactoriesTest();
  ZArrayFactoriesTest();
  CArrayFactoriesTest();
  DArrayFactoriesTest();
  SArrayFactoriesTest();
  I2ArrayFactoriesTest();
  I4ArrayFactoriesTest();
  I8ArrayFactoriesTest();
  U2ArrayFactoriesTest();
  U4ArrayFactoriesTest();
  U8ArrayFactoriesTest();

  return 0;
}


/*
 * TestStatus ()
 *
 * Routine to check that the status code status->statusCode agrees with one of
 * the codes specified in the space-delimited string ignored; if not,
 * exit to the system with code exitcode.
 *
 */
static void
TestStatus (LALStatus *status, const char *ignored, int exitcode)
{
  char  str[64];
  char *tok;

  if (verbose)
  {
    REPORTSTATUS (status);
  }

  if (strncpy (str, ignored, sizeof (str)))
  {
    if ((tok = strtok (str, " ")))
    {
      do
      {
        if (status->statusCode == atoi (tok))
        {
          return;
        }
      }
      while ((tok = strtok (NULL, " ")));
    }
    else
    {
      if (status->statusCode == atoi (tok))
      {
        return;
      }
    }
  }

  fprintf (stderr, "\nExiting to system with code %d\n", exitcode);
  exit (exitcode);
}


/*
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
void
ClearStatus (LALStatus *status)
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    DETATCHSTATUSPTR (status);
  }
}


/*
 * Usage ()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */
static void
Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h         print this message\n");
  fprintf (stderr, "  -q         quiet: run silently\n");
  fprintf (stderr, "  -v         verbose: print extra information\n");
  fprintf (stderr, "  -d level   set lalDebugLevel to level\n");
  exit (exitcode);
}


/*
 * ParseOptions ()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */
static void
ParseOptions (int argc, char *argv[])
{
  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        freopen ("/dev/null", "w", stderr);
        freopen ("/dev/null", "w", stdout);
        break;

      case 'h':
        Usage (argv[0], 0);
        break;

      default:
        Usage (argv[0], 1);
    }

  }

  if (optind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}

