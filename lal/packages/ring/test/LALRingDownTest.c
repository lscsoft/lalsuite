/****************** <lalVerbatim file="LALRingDownTestCV">
Author: Tibbits, M. M.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{LALRingDownTest.c}}
\label{s:LALRingDownTest.c}

A program to test \texttt{LALDRingDown()}.
- Note only the double precision is tested because both are derived from the same code.
\subsubsection*{Usage}

\begin{verbatim}
./LALRingDownTest [options]
Options:
  -h             print usage message
  -q             quiet: run silently
  -v             verbose: print extra information
  -d level       set lalDebugLevel to level
\end{verbatim}

This program tests the function
\texttt{LALDRingDown()}, which calculates the moment
of a given data set.

First, it tests that the correct error codes 
are generated for the following error conditions (tests in
\textit{italics} are not performed if \verb+LAL_NEDEBUG+ is set, as
the corresponding checks in the code are made using the ASSERT macro):
\begin{itemize}
\item \textit{null pointer to output structure}
\item \textit{null pointer to input structure}
\item \textit{null pointer to data member of input structure}
\item \textit{null pointer to data member of data member of input structure}
\item \textit{zero length}
\end{itemize}

For each successful test
(both of these valid data and the invalid ones described above), it
prints ``\texttt{PASS}'' to standard output; if a test fails, it
prints ``\texttt{FAIL}''.

\subsubsection*{Exit codes}

\subsubsection*{Uses}

\begin{verbatim}
LALDRingDown()
LALSRingDown()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALRingDownTestCV}}

******************************************************* </lalLaTeX> */

#include <lal/LALStdlib.h>

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <config.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALRingDown.h>
#include "CheckStatus.h"


NRCSID (LALRINGDOWNTESTC, "$Id$");


/*  constants  */
#define LALRINGDOWNTESTC_TRUE     1
#define LALRINGDOWNTESTC_FALSE    0

extern char	*optarg;
extern int	optind;


/*  Setting Global debug level  */
int	lalDebugLevel	= 0;


/*  Setting variables to parse command line flags  */
BOOLEAN	optVerbose	= LALRINGDOWNTESTC_FALSE;
UINT4	optLength	= 0;


static void Usage 
(
	const char	*program,
	int		exitflag
);

static void ParseOptions
(
	int		argc,
	char		*argv[]
);


/*************** <lalErrTab > */
#define	LALRINGDOWNTESTC_ENOM	0
#define	LALRINGDOWNTESTC_EARG	1
#define	LALRINGDOWNTESTC_ECHK	2
#define	LALRINGDOWNTESTC_EFLS	3
#define	LALRINGDOWNTESTC_EUSE	4
#define	LALRINGDOWNTESTC_ENULL	5
#define	LALRINGDOWNTESTC_EALOC	6
#define	LALRINGDOWNTESTC_MSGENOM	"Nominal exit"
#define	LALRINGDOWNTESTC_MSGEARG	"Error parsing command-line arguments"
#define	LALRINGDOWNTESTC_MSGECHK	"Error checking failed to catch bad data"
#define	LALRINGDOWNTESTC_MSGEFLS	"Incorrect answer for valid data"
#define	LALRINGDOWNTESTC_MSGEUSE	"Bad user-entered data"
#define	LALRINGDOWNTESTC_MSGENULL	"Null Pointer."
#define LALRINGDOWNTESTC_MSGEALOC	"Memory Allocation Error"
/***************************** </lalErrTab> */



int main( int argc, char *argv[] )
{

	static	LALStatus	status;

	/* Variable declarations */
	COMPLEX16RingDownParams   *nullParams;
	COMPLEX16RingDownParams   params;
	COMPLEX16FrequencySeries  *nullOutput;
	COMPLEX16FrequencySeries  output;
	COMPLEX16Sequence         temp;
	INT4			  code;

	/*  Initialize Variables  */
	nullParams = NULL;
	nullOutput = NULL;

	params.f0 = 3.0;
	params.df = 3.0;
	params.t  = 3.0;
	params.n  = 1;

	ParseOptions( argc, argv );

	printf("\n\nMESG: %s \n",LALRINGDOWNTESTC);

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
	/* test behavior for null pointer to output structure */
	LALDRingDown(&status, nullOutput, &params);
	if ( ( code = CheckStatus(&status, LALRINGDOWNH_ENULL, LALRINGDOWNH_MSGENULL,
					LALRINGDOWNTESTC_ECHK, LALRINGDOWNTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: non-null pointer to output structure results in error:\n");
	printf("       \"%s\"\n", LALRINGDOWNH_MSGENULL);

  }
#endif

	printf("MESG:	More option available from command line.\n");
	printf("MESG:	Type LALRingDownTest -h for options.\n");

	/* normal exit */
	return LALRINGDOWNTESTC_ENOM;
}


/*  CODE LISTED BEYOND THIS POINT WAS TAKEN FROM STOCHASTICCROSSCORRELATIONTEST.C  */

/*
 * Usage ()
 *
 * Prints a usage message for program program and exits with code exitcode.
 *
 */

static void Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h             print this message\n");
  fprintf (stderr, "  -q             quiet: run silently\n");
  fprintf (stderr, "  -v             verbose: print extra information\n");
  fprintf (stderr, "  -d level       set lalDebugLevel to level\n");
  exit (exitcode);
}

/*
 * ParseOptions ()
 *
 * Parses the argc - 1 option strings in argv[].
 *
 */

static void ParseOptions (int argc, char *argv[])
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

      case 'v': /* optVerbose */
        optVerbose = LALRINGDOWNTESTC_TRUE;
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
