/****************************** <lalVerbatim file="LALWienerTestCV">
Author: Dwyer, Steven J.
$Id$
******************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{LALWienerTest.c}}
\label{s:LALWienerTest.c}

A program to test \texttt{LALUnFormattedWiener()} and \texttt{LALFormattedWiener()}.

\subsubsection*{Usage}

\begin{verbatim}
./LALWienerTest [options]
Options:
  -h             print usage message
  -q             quiet: run silently
  -v             verbose: print extra information
  -d level       set lalDebugLevel to level
\end{verbatim}

This program tests the functions
\texttt{LALUnFormattedWiener()} and \texttt{LALFormattedWiener()}, which perform a wiener filter of the data.

First, it tests that the correct error codes 
are generated for the following error conditions (tests in
\textit{italics} are not performed if \verb+LAL_NEDEBUG+ is set, as
the corresponding checks in the code are made using the ASSERT macro):
\begin{itemize}
\item \textit{null pointer to input structure}
\item \textit{null pointer to output structure}
\item \textit{template is longer than source vector}
\end{itemize}

It then verifies that the correct result is obtained from the following casses:
\begin{enumerate}
\item source is longer than template vector.
\item source is same size as template vector.
\end{enumerate}

For each successful test
(both of these valid data and the invalid ones described above), it
prints ``\texttt{PASS}'' to standard output; if a test fails, it
prints ``\texttt{FAIL}''.

\subsubsection*{Exit codes}

\subsubsection*{Uses}

\begin{verbatim}
LALUnFormattedWiener()
LALFormattedWiener()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALWienerTestCV}}

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

#include <lal/LALWiener.h>
#include "CheckStatus.h"

/* Define RCS ID string */
NRCSID(LALWIENERTESTC,"$Id$");

/*  constants  */
#define LALWIENERTESTC_TRUE     1
#define LALWIENERTESTC_FALSE    0

extern char	*optarg;
extern int	optind;

/*  Setting Global debug level  */
int	lalDebugLevel	= 0;

/*  Setting variables to parse command line flags  */
BOOLEAN	optVerbose	= LALWIENERTESTC_FALSE;
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
#define	LALWIENERTESTC_ENOM	0
#define	LALWIENERTESTC_EARG	1
#define	LALWIENERTESTC_ECHK	2
#define	LALWIENERTESTC_EFLS	3
#define	LALWIENERTESTC_EUSE	4
#define	LALWIENERTESTC_ENULL	5
#define	LALWIENERTESTC_EALOC	6
#define	LALWIENERTESTC_MSGENOM	"Nominal exit"
#define	LALWIENERTESTC_MSGEARG	"Error parsing command-line arguments"
#define	LALWIENERTESTC_MSGECHK	"Error checking failed to catch bad data"
#define	LALWIENERTESTC_MSGEFLS	"Incorrect answer for valid data"
#define	LALWIENERTESTC_MSGEUSE	"Bad user-entered data"
#define	LALWIENERTESTC_MSGENULL	"Null Pointer."
#define LALWIENERTESTC_MSGEALOC	"Memory Allocation Error"
/***************************** </lalErrTab> */


int main( int argc, char *argv[] )
{
	static LALStatus stat;

	/* Variable declarations */
	WienerUnFormattedInput	test;
	WienerOutput		result;
	WienerUnFormattedInput	*nullTest;
	WienerOutput		*nullResult;
        INT2			count;
	INT4			code;

	ParseOptions( argc, argv );

	printf("\n\nMESG: %s \n",LALWIENERTESTC);

	/*  Initialize Variables  */
	nullTest   = NULL;
	nullResult = NULL;

	/*  Create input and result vectors  */
	test.h = NULL;
	LALSCreateVector(&stat, &(test.h), 6);

	test.s = NULL;
	LALSCreateVector(&stat, &(test.s), 4);

	result.q = NULL;
	LALSCreateVector(&stat, &(result.q), 10);

	result.z = NULL;
	LALSCreateVector(&stat, &(result.z), 4);

	/*  Initialize Input  */
	(test.h)->data[0] = 0;
	(test.h)->data[1] = 3;
	(test.h)->data[2] = 4;
       	(test.h)->data[3] = 3;
	(test.h)->data[4] = 4;
	(test.h)->data[5] = 0;

	(test.s)->data[0] = 0;
	(test.s)->data[1] = 0.5;
	(test.s)->data[2] = 0.5;
	(test.s)->data[3] = 0;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
	printf("MESG: Test Error Codes I \n");

	/* test behavior for null pointer to input structure */
	LALUnFormattedWiener( &stat, &result, nullTest );
	if ( ( code = CheckStatus(&stat, LALWIENERH_ENULL, LALWIENERH_MSGENULL,
					LALWIENERTESTC_ECHK, LALWIENERTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to input structure results in error:\n");
	printf("       \"%s\"\n", LALWIENERH_MSGENULL);

	/* test behavior for null pointer to output structure */
	LALUnFormattedWiener( &stat, nullResult, &test );
	if ( ( code = CheckStatus(&stat, LALWIENERH_ENULL, LALWIENERH_MSGENULL,
					LALWIENERTESTC_ECHK, LALWIENERTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to output structure results in error:\n");
	printf("       \"%s\"\n", LALWIENERH_MSGENULL);
  }
#endif

	/************ First Test *************/
	printf("MESG: Test LALUnFormattedWiener \n");

	/*  Variables for first test are initialized above to facilitate error checking  */
 
	LALUnFormattedWiener( &stat, &result, &test );
	

	if ( lalDebugLevel > 0 )
	{
             
		printf("MESG: q \n");
		for (count = 0; count < 6; count++)
			printf("MESG: %f \n", result.q->data[count]);

		printf("MESG: z \n");
		for (count = 0; count < 4; count++)
			printf("MESG: %f \n", result.z->data[count]);
	}



	/*  Clean up input and result vectors  */
	LALSDestroyVector(&stat, &(test.h));

	LALSDestroyVector(&stat, &(test.s));    

	LALSDestroyVector(&stat, &(result.q));

	LALSDestroyVector(&stat, &(result.z));


	/************ Second Test *************/
	printf("MESG: Test LALUnFormattedWiener II \n");

	/*  Create input and result vectors  */
	test.h = NULL;
	LALSCreateVector(&stat, &(test.h), 4);
	
	test.s = NULL;
	LALSCreateVector(&stat, &(test.s), 4); 
	        
	result.q = NULL;
	LALSCreateVector(&stat, &(result.q), 4);

	result.z = NULL;
	LALSCreateVector(&stat, &(result.z), 4);

	/*  Initialize Input  */
	(test.h)->data[0] = 0;
	(test.h)->data[1] = 3;
	(test.h)->data[2] = 4;
       	(test.h)->data[3] = 3;

	(test.s)->data[0] = 0;
	(test.s)->data[1] = 0.5;
	(test.s)->data[2] = 0.5;
        (test.s)->data[3] = 0;

	LALUnFormattedWiener( &stat, &result, &test );

	if ( lalDebugLevel > 0 )
	{
		printf("MESG: q \n");
		for (count = 0; count < (INT4)result.q->length; count++)
			printf("MESG: %f \n", result.q->data[count]);

		printf("MESG: z \n");
		for (count = 0; count < (INT4)result.z->length; count++)
			printf("MESG: %f \n", result.z->data[count]);
	}

	/*  Clean up input and result vectors  */
	LALSDestroyVector(&stat, &(test.h));

	LALSDestroyVector(&stat, &(test.s));    
 
	LALSDestroyVector(&stat, &(result.q));

	LALSDestroyVector(&stat, &(result.z));

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
	printf("MESG: Test Error Codes II \n");

	/*  Create input and result vectors  */
	test.h = NULL;
	LALSCreateVector(&stat, &(test.h), 3);

	test.s = NULL;
	LALSCreateVector(&stat, &(test.s), 4); 
        
	result.q = NULL;
	LALSCreateVector(&stat, &(result.q), 7);

	result.z = NULL;
	LALSCreateVector(&stat, &(result.z), 4);


	/*  Initialize Input  */
	(test.h)->data[0] = 0;
	(test.h)->data[1] = 3;
	(test.h)->data[2] = 4;


	(test.s)->data[0] = 0;
	(test.s)->data[1] = 0.5;
	(test.s)->data[1] = 0.5;
        (test.s)->data[3] = 0;	    

	LALUnFormattedWiener( &stat, &result, &test );
	if ( ( code = CheckStatus(&stat, LALWIENERH_ESIZE, LALWIENERH_MSGESIZE,
					LALWIENERTESTC_ECHK, LALWIENERTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: wrong template vector size results in error:\n");
	printf("       \"%s\"\n", LALWIENERH_MSGESIZE);

	/*  Clean up input and result vectors  */
	LALSDestroyVector(&stat, &(test.h));

	LALSDestroyVector(&stat, &(test.s));    
 
	LALSDestroyVector(&stat, &(result.q));

	LALSDestroyVector(&stat, &(result.z));
  }
#endif

	return LALWIENERTESTC_ENOM;
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
        optVerbose = LALWIENERTESTC_TRUE;
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

