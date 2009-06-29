/*
*  Copyright (C) 2007 Jolien Creighton, Matt Tibbits
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

/****************** <lalVerbatim file="LALMomentTestCV">
Author: Tibbits, M. M.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{LALMomentTest.c}}
\label{s:LALMomentTest.c}

A program to test \texttt{LALDMoment()}.
- Note only the double precision is tested because both are derived from the same code.
\subsubsection*{Usage}

\begin{verbatim}
./LALMomentTest [options]
Options:
  -h             print usage message
  -q             quiet: run silently
  -v             verbose: print extra information
  -d level       set lalDebugLevel to level
\end{verbatim}

This program tests the function
\texttt{LALDMoment()}, which calculates the moment
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

It then verifies that the correct moment (value and units) is
generated for each of the following simple test cases:
\begin{enumerate}
\item data set all same value, find moments 2-5.
\item mixed data set, find moments 2-5.
\item evenly distributed data set, find moments 2-5.
\end{enumerate}

For each successful test
(both of these valid data and the invalid ones described above), it
prints ``\texttt{PASS}'' to standard output; if a test fails, it
prints ``\texttt{FAIL}''.

\subsubsection*{Exit codes}

\subsubsection*{Uses}

\begin{verbatim}
LALDMoment()
LALSMoment()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{LALMomentTestCV}}

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

#include <lal/LALMoment.h>
#include "CheckStatus.h"


NRCSID (LALMOMENTTESTC, "$Id$");


/*  constants  */
#define LALMOMENTTESTC_TRUE     1
#define LALMOMENTTESTC_FALSE    0

extern char	*optarg;
extern int	optind;


/*  Setting Global debug level  */
int	lalDebugLevel	= 0;


/*  Setting variables to parse command line flags  */
BOOLEAN	optVerbose	= LALMOMENTTESTC_FALSE;
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
#define	LALMOMENTTESTC_ENOM	0
#define	LALMOMENTTESTC_EARG	1
#define	LALMOMENTTESTC_ECHK	2
#define	LALMOMENTTESTC_EFLS	3
#define	LALMOMENTTESTC_EUSE	4
#define	LALMOMENTTESTC_ENULL	5
#define	LALMOMENTTESTC_EALOC	6
#define	LALMOMENTTESTC_MSGENOM	"Nominal exit"
#define	LALMOMENTTESTC_MSGEARG	"Error parsing command-line arguments"
#define	LALMOMENTTESTC_MSGECHK	"Error checking failed to catch bad data"
#define	LALMOMENTTESTC_MSGEFLS	"Incorrect answer for valid data"
#define	LALMOMENTTESTC_MSGEUSE	"Bad user-entered data"
#define	LALMOMENTTESTC_MSGENULL	"Null Pointer."
#define LALMOMENTTESTC_MSGEALOC	"Memory Allocation Error"
/***************************** </lalErrTab> */



int main( int argc, char *argv[] )
{

	static	LALStatus	status;

	/* Variable declarations */
	REAL8			data[40];
	REAL8			*result;
	REAL8			*nullResult;
	INT4			length;
	INT4			whichMoment;
	INT4			iterator;
	INT4			code;
	REAL8			testOutput[7];
	REAL8Sequence		*sequence;
	REAL8Sequence		*nullSequence;


	const  INT4     	constantData[]	=
				{ 2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,
				2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4,2,4};

	const  INT4     	constantData2[]	=
				{ 10,20,10,20,10,20,10,20,10,20};

	length		=  40;
	nullResult	=  NULL;
	nullSequence	=  NULL;
	whichMoment	=  3;

	result	=  (REAL8*) LALMalloc(sizeof(REAL8));

	if( !result )
	{
		return LALMOMENTTESTC_EALOC;
	}

	sequence =  (REAL8Sequence*) LALMalloc(sizeof(REAL8Sequence));

	if( !sequence )
	{
		return LALMOMENTTESTC_EALOC;
	}

	for ( iterator = 0;  iterator < length;  iterator++ )
	{
		data[iterator] = ((REAL8)(20 - iterator));
	}

	sequence->length = length;

	ParseOptions( argc, argv );

	printf("\n\nMESG: %s \n",LALMOMENTTESTC);

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
	/* test behavior for null pointer to input structure */
	LALDMoment(&status, result, nullSequence, whichMoment);
	if ( ( code = CheckStatus(&status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL,
					LALMOMENTTESTC_ECHK, LALMOMENTTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to input structure results in error:\n");
	printf("       \"%s\"\n", LALMOMENTH_MSGENULL);


	/* test behavior for zero length of input sequnce */
	sequence->length = 0;

	LALDMoment(&status, result, sequence, whichMoment);
	if ( ( code = CheckStatus(&status, LALMOMENTH_ELNTH, LALMOMENTH_MSGELNTH,
					LALMOMENTTESTC_ECHK, LALMOMENTTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to input structure results in error:\n");
	printf("       \"%s\"\n", LALMOMENTH_MSGELNTH);

	/*  Set proper length  */
	sequence->length = length;


	/* test behavior for null pointer to output structure */
	LALDMoment(&status, nullResult, sequence, whichMoment);
	if ( ( code = CheckStatus(&status, LALMOMENTH_ENULL, LALMOMENTH_MSGENULL,
					LALMOMENTTESTC_ECHK, LALMOMENTTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: non-null pointer to output structure results in error:\n");
	printf("       \"%s\"\n", LALMOMENTH_MSGENULL);
  }
#endif


	/*  Setting the parameter data pointer to equal the local data variable  */
	sequence->data = data;

	/*  Setting all data points to 5 for first test  */
	for ( iterator = 0;  iterator < length;  iterator++ )
	{
		data[iterator] = 5.0;
	}

	/**********  First Test  **********/
	for(whichMoment = 2; whichMoment < 6; whichMoment++)
	{
		LALDMoment(&status, result, sequence, whichMoment);
		testOutput[whichMoment] = *result;

		if( lalDebugLevel > 0 )
		{
			printf("MESG:  whichMoment := %d\n",whichMoment);
			printf("MESG:  testOutput[%d] := %f\n\n",whichMoment, testOutput[whichMoment]);
		}
	}

	/*  Data set to equal distribution from 20 to -20  */
	for ( iterator = 0;  iterator < length;  iterator++ )
	{
		data[iterator] = ((REAL8)(20 - iterator));
	}

	/**********  Second Test  **********/
	for(whichMoment = 2; whichMoment < 6; whichMoment++)
	{
		LALDMoment(&status, result, sequence, whichMoment);
		testOutput[whichMoment] = *result;

		if( lalDebugLevel > 0 )
		{
			printf("MESG:  whichMoment := %d\n",whichMoment);
			printf("MESG:  testOutput[%d] := %f\n\n",whichMoment, testOutput[whichMoment]);
		}
	}

	/*  Data set to 50% 2's & 50% 4's  */
	for ( iterator = 0;  iterator < length;  iterator++ )
	{
		data[iterator] = ((REAL8)(constantData[iterator]));
	}

	/**********  Third Test  **********/
	for(whichMoment = 2; whichMoment < 6; whichMoment++)
	{
		LALDMoment(&status, result, sequence, whichMoment);
		testOutput[whichMoment] = *result;

		if( lalDebugLevel > 0 )
		{
			printf("MESG:  whichMoment := %d\n",whichMoment);
			printf("MESG:  testOutput[%d] := %f\n\n",whichMoment, testOutput[whichMoment]);
		}
	}

	/*  Change length  */
	sequence->length = length = 10;

	/*  Data set to 50% 2's & 50% 4's  */
	for ( iterator = 0;  iterator < length;  iterator++ )
	{
		data[iterator] = ((REAL8)(constantData2[iterator]));
	}

	/**********  Fourth Test  **********/
	for(whichMoment = 2; whichMoment < 6; whichMoment++)
	{
		LALDMoment(&status, result, sequence, whichMoment);
		testOutput[whichMoment] = *result;

		if( lalDebugLevel > 0 )
		{
			printf("MESG:  whichMoment := %d\n",whichMoment);
			printf("MESG:  testOutput[%d] := %f\n\n",whichMoment, testOutput[whichMoment]);
		}
	}

	printf("MESG:	More option available from command line.\n");
	printf("MESG:	Type LALMomentTest -h for options.\n");

	/* normal exit */
	return LALMOMENTTESTC_ENOM;
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
        optVerbose = LALMOMENTTESTC_TRUE;
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
