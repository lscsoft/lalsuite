/****************** <lalVerbatim file="NRRandomTestCV">
Author: Tibbits, M. M.
$Id$
********************************* </lalVerbatim> */

/********************************************************** <lalLaTeX>
\subsection{Program \texttt{NRRandomTest.c}}
\label{s:NRRandomTest.c}

A program to test \texttt{NRRan3()}, \texttt{NRRan4()}, and \texttt{NRGasDev()}.

\begin{verbatim}
./NRRandomTest [options]
Options:
  -h             print usage message
  -q             quiet: run silently
  -v             verbose: print extra information
  -d level       set lalDebugLevel to level
\end{verbatim}

This program tests the function

It tests that the correct error codes 
are generated for the following error conditions (tests in
\textit{italics} are not performed if \verb+NR_NEDEBUG+ is set, as
the corresponding checks in the code are made using the ASSERT macro):
\begin{itemize}
\item \textit{null pointer to output structure}
\item \textit{null pointer to input structure}
\item \textit{check out of bounds on params}
\end{itemize}

Note that the check on the NULL pointer to params is only done to the 
\texttt{NRGasDev()} function.

It then tests each of the three functions against the concept that each function
should return the same sequence of numbers twice if the same the seed is used both times.
\\
So for instance, when you pass seed = -248673140 to NRRan3(), you should get the following first five numbers:
\begin{itemize}
\item{0.294567555189132690429687500000}
\item{0.014104902744293212890625000000}
\item{0.089319035410881042480468750000}
\item{0.886621534824371337890625000000}
\item{0.795643210411071777343750000000}
\end{itemize}

For each successful test
(both of these valid data and the invalid ones described above), it
prints ``\texttt{PASS}'' to standard output; if a test fails, it
prints ``\texttt{FAIL}''.

\subsubsection*{Exit codes}

\input{NRRandomTestCE}

\subsubsection*{Uses}

\begin{verbatim}
NRRan3()
NRRan4()
NRGasDev()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{NRRandomTestCV}}

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

#include <lal/NRRandom.h>
#include "CheckStatus.h"


NRCSID (NRRANDOMTESTC, "$Id$");


/*  constants  */
#define NRRANDOMTESTC_TRUE     1
#define NRRANDOMTESTC_FALSE    0

extern char	*optarg;
extern int	optind;


/*  Setting Global debug level  */
int	lalDebugLevel	= 0;


/*  Setting variables to parse command line flags  */
BOOLEAN	optVerbose	= NRRANDOMTESTC_FALSE;
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


/*************** <lalErrTable file="NRRandomTestCE" > */
#define	NRRANDOMTESTC_ENOM      0
#define	NRRANDOMTESTC_EARG      1
#define	NRRANDOMTESTC_ECHK      2
#define	NRRANDOMTESTC_EFLS      3
#define	NRRANDOMTESTC_EUSE      4
#define	NRRANDOMTESTC_ENULL     5
#define	NRRANDOMTESTC_EALOC     6
#define	NRRANDOMTESTC_ESEED     7
#define	NRRANDOMTESTC_MSGENOM   "Nominal exit"
#define	NRRANDOMTESTC_MSGEARG   "Error parsing command-line arguments"
#define	NRRANDOMTESTC_MSGECHK   "Error checking failed to catch bad data"
#define	NRRANDOMTESTC_MSGEFLS   "Incorrect answer for valid data"
#define	NRRANDOMTESTC_MSGEUSE   "Bad user-entered data"
#define	NRRANDOMTESTC_MSGENULL  "Null Pointer."
#define NRRANDOMTESTC_MSGEALOC  "Memory Allocation Error"
#define NRRANDOMTESTC_MSGESEED  "Same seed should logically produce same set."
/***************************** </lalErrTable> */

#define NRRANDOMTESTC_LENGTH 100

int main( int argc, char *argv[] )
{

	static	LALStatus	status;

	/* Variable declarations */
	REAL4Vector	*storeResult1;
	REAL4Vector	*storeResult2;
	REAL4Vector	*temp;
	REAL4		*nullOutput;
	REAL4		output;
	INT4		*nullSeed;
	INT4		iterator;
	INT2		badTest;
	INT4		count;
	INT4		*seed;
	INT4		code;
	INT2		test;

	/*  Initialize Variables  */
	storeResult1	= NULL;
	storeResult2	= NULL;
	nullOutput	= NULL;
	nullSeed	= NULL;
	seed		= NULL;
	temp		= NULL;
	iterator	= 0;
	badTest		= 4;
	count		= 0;
	test		= 0;

	printf("\n\nMESG: %s \n",NRRANDOMTESTC);

	/*  Allocate space for the seed pointer  */
	seed = LALMalloc(sizeof(INT4));
	*seed = -248673140;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {

	if(!(seed))
	{
		return NRRANDOMTESTC_EALOC;
	}

	/********* Error Tests on NRRan3() *********/
	printf("\nMESG: Testing NRRan3()\n");
	/*  test behavior for null pointer to output structure  */
	NRRan3(&status, nullOutput, seed);
	if ( ( code = CheckStatus(&status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL,
					NRRANDOMTESTC_ECHK, NRRANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to output structure results in error:\n");
	printf("       \"%s\"\n", NRRANDOMH_MSGENULL);

	/*  test behavior for null pointer to input structure  */
	NRRan3(&status, &output, nullSeed );
	if ( ( code = CheckStatus(&status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL,
					NRRANDOMTESTC_ECHK, NRRANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to seed structure results in error:\n");
	printf("       \"%s\"\n", NRRANDOMH_MSGENULL);



	/********* Error Tests on NRRan4() *********/
	printf("\nMESG: Testing NRRan4()\n");
	/*  test behavior for null pointer to output structure  */
	NRRan4(&status, nullOutput, seed);
	if ( ( code = CheckStatus(&status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL,
					NRRANDOMTESTC_ECHK, NRRANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to output structure results in error:\n");
	printf("       \"%s\"\n", NRRANDOMH_MSGENULL);

	/*  test behavior for null pointer to input structure  */
	NRRan4(&status, &output, nullSeed );
	if ( ( code = CheckStatus(&status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL,
					NRRANDOMTESTC_ECHK, NRRANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS:   null pointer to seed structure results in error:\n");
	printf("       \"%s\"\n", NRRANDOMH_MSGENULL);



	/********* Error Tests on NRGasDev() *********/
	printf("\nMESG: Testing NRGasDev()\n");
	/*  test behavior for null pointer to output structure  */
	NRGasDev(&status, nullOutput, seed, &test);
	if ( ( code = CheckStatus(&status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL,
					NRRANDOMTESTC_ECHK, NRRANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to output structure results in error:\n");
	printf("       \"%s\"\n", NRRANDOMH_MSGENULL);

	/*  test behavior for null pointer to input structure  */
	NRGasDev(&status, &output, nullSeed, &test );
	if ( ( code = CheckStatus(&status, NRRANDOMH_ENULL, NRRANDOMH_MSGENULL,
					NRRANDOMTESTC_ECHK, NRRANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS:   null pointer to seed structure results in error:\n");
	printf("       \"%s\"\n", NRRANDOMH_MSGENULL);


	/*  test behavior for bad test value  */
	test = 2;
	NRGasDev(&status, &output, seed, &test );
	if ( ( code = CheckStatus(&status, NRRANDOMH_ETEST, NRRANDOMH_MSGETEST,
					NRRANDOMTESTC_ECHK, NRRANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: out of bounds value for the param results in error:\n");
	printf("       \"%s\"\n", NRRANDOMH_MSGETEST);
	test = 0;

  }
#endif

	LALSCreateVector(&status, &storeResult1, NRRANDOMTESTC_LENGTH);
	LALSCreateVector(&status, &storeResult2, NRRANDOMTESTC_LENGTH);
	LALSCreateVector(&status, &temp, NRRANDOMTESTC_LENGTH);

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
	if(!(storeResult1))
	{
		return NRRANDOMTESTC_EALOC;
	}

	if(!(storeResult2))
	{
		return NRRANDOMTESTC_EALOC;
	}

	if(!(temp))
	{
		return NRRANDOMTESTC_EALOC;
	}
  }
#endif
	printf("MESG: Seed Value := %d\n", *seed);
	printf("MESG: Using NRRan3()\n");
	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan3(&status, &(storeResult1->data[iterator]), seed);
	}

	*seed = -658741232;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan3(&status, &(temp->data[iterator]), seed);
	}

	*seed = -248673140;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan3(&status, &(storeResult2->data[iterator]), seed);
	}

	printf("\n");

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; )
	{
		printf(" %1.20f %1.20f", storeResult1->data[iterator], storeResult2->data[iterator]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+1], storeResult2->data[iterator+1]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+2], storeResult2->data[iterator+2]);
		printf(" %1.20f %1.20f\n", storeResult1->data[iterator+3], storeResult2->data[iterator+3]);

		if ( storeResult2->data[iterator] == storeResult1->data[iterator])
			count++;
		if ( storeResult2->data[iterator+1] == storeResult1->data[iterator+1])
			count++;
		if ( storeResult2->data[iterator+2] == storeResult1->data[iterator+2])
			count++;
		if ( storeResult2->data[iterator+3] == storeResult1->data[iterator+3])
			count++;

		iterator +=4;
	}

	printf("\n");
	*seed = -248673140;
	printf("MESG: Seed Value := %d\n", *seed);
	printf("MESG: Using NRRan4()\n");
	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan4(&status, &(storeResult1->data[iterator]), seed);
	}

	*seed = -658741232;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan4(&status, &(temp->data[iterator]), seed);
	}

	*seed = -248673140;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan4(&status, &(storeResult2->data[iterator]), seed);
	}

	printf("\n");

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; )
	{
		printf(" %1.20f %1.20f", storeResult1->data[iterator], storeResult2->data[iterator]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+1], storeResult2->data[iterator+1]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+2], storeResult2->data[iterator+2]);
		printf(" %1.20f %1.20f\n", storeResult1->data[iterator+3], storeResult2->data[iterator+3]);

		if ( storeResult2->data[iterator] == storeResult1->data[iterator])
			count++;
		if ( storeResult2->data[iterator+1] == storeResult1->data[iterator+1])
			count++;
		if ( storeResult2->data[iterator+2] == storeResult1->data[iterator+2])
			count++;
		if ( storeResult2->data[iterator+3] == storeResult1->data[iterator+3])
			count++;

		iterator +=4;
	}

	printf("\n");

	*seed = -248673140;
	printf("MESG: Seed Value := %d\n", *seed);
	printf("MESG: Using NRGasDev()\n");
	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRGasDev(&status, &(storeResult1->data[iterator]), seed, &test);
	}

	*seed = -658741232;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRGasDev(&status, &(temp->data[iterator]), seed, &test);
	}

	*seed = -248673140;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRGasDev(&status, &(storeResult2->data[iterator]), seed, &test);
	}

	printf("\n");

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; )
	{
		printf(" %1.20f %1.20f", storeResult1->data[iterator], storeResult2->data[iterator]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+1], storeResult2->data[iterator+1]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+2], storeResult2->data[iterator+2]);
		printf(" %1.20f %1.20f\n", storeResult1->data[iterator+3], storeResult2->data[iterator+3]);

		if ( storeResult2->data[iterator] == storeResult1->data[iterator])
			count++;
		if ( storeResult2->data[iterator+1] == storeResult1->data[iterator+1])
			count++;
		if ( storeResult2->data[iterator+2] == storeResult1->data[iterator+2])
			count++;
		if ( storeResult2->data[iterator+3] == storeResult1->data[iterator+3])
			count++;

		iterator +=4;
	}

	printf("\n");

	*seed = -322140768;
	printf("MESG: Seed Value := %d\n", *seed);
	printf("MESG: Using NRRan3()\n");
	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan3(&status, &(storeResult1->data[iterator]), seed);
	}

	*seed = -658741232;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan3(&status, &(temp->data[iterator]), seed);
	}

	*seed = -322140768;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan3(&status, &(storeResult2->data[iterator]), seed);
	}

	printf("\n");

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; )
	{
		printf(" %1.20f %1.20f", storeResult1->data[iterator], storeResult2->data[iterator]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+1], storeResult2->data[iterator+1]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+2], storeResult2->data[iterator+2]);
		printf(" %1.20f %1.20f\n", storeResult1->data[iterator+3], storeResult2->data[iterator+3]);

		if ( storeResult2->data[iterator] == storeResult1->data[iterator])
			count++;
		if ( storeResult2->data[iterator+1] == storeResult1->data[iterator+1])
			count++;
		if ( storeResult2->data[iterator+2] == storeResult1->data[iterator+2])
			count++;
		if ( storeResult2->data[iterator+3] == storeResult1->data[iterator+3])
			count++;

		iterator +=4;
	}

	printf("\n");

	*seed = -322140768;
	printf("MESG: Seed Value := %d\n", *seed);
	printf("MESG: Using NRRan4()\n");
	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan4(&status, &(storeResult1->data[iterator]), seed);
	}

	*seed = -658741232;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan4(&status, &(temp->data[iterator]), seed);
	}

	*seed = -322140768;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan4(&status, &(storeResult2->data[iterator]), seed);

	}

	printf("\n");

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; )
	{
		printf(" %1.20f %1.20f", storeResult1->data[iterator], storeResult2->data[iterator]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+1], storeResult2->data[iterator+1]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+2], storeResult2->data[iterator+2]);
		printf(" %1.20f %1.20f\n", storeResult1->data[iterator+3], storeResult2->data[iterator+3]);

		if ( storeResult2->data[iterator] == storeResult1->data[iterator])
			count++;
		if ( storeResult2->data[iterator+1] == storeResult1->data[iterator+1])
			count++;
		if ( storeResult2->data[iterator+2] == storeResult1->data[iterator+2])
			count++;
		if ( storeResult2->data[iterator+3] == storeResult1->data[iterator+3])
			count++;

		iterator +=4;
	}

	printf("\n");

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
	if( count == (NRRANDOMTESTC_LENGTH)*5)
	{
		printf("PASS: Same seed produces same set.\n\n");
	}
	else
	{
		printf("FAIL: Same seed does not produces same set.\n\n");
		return NRRANDOMTESTC_ESEED;
	}

	/*  Reinitialize Counter  */
	count = 0;
  }
#endif

	*seed	= -322140768;
	printf("MESG: Seed Value := %d\n", *seed);
	printf("MESG: Using NRRan3()\n");
	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan3(&status, &(storeResult1->data[iterator]), seed);
	}

	*seed	= -322140768;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan3(&status, &(storeResult2->data[iterator]), seed);
	}

	printf("\n");

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; )
	{
		printf(" %1.20f %1.20f", storeResult1->data[iterator], storeResult2->data[iterator]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+1], storeResult2->data[iterator+1]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+2], storeResult2->data[iterator+2]);
		printf(" %1.20f %1.20f\n", storeResult1->data[iterator+3], storeResult2->data[iterator+3]);

		if ( storeResult2->data[iterator] == storeResult1->data[iterator])
			count++;
		if ( storeResult2->data[iterator+1] == storeResult1->data[iterator+1])
			count++;
		if ( storeResult2->data[iterator+2] == storeResult1->data[iterator+2])
			count++;
		if ( storeResult2->data[iterator+3] == storeResult1->data[iterator+3])
			count++;

		iterator +=4;
	}

	printf("\n");

	*seed	= -322140768;
	printf("MESG: Seed Value := %d\n", *seed);
	printf("MESG: Using NRRan4()\n");
	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan4(&status, &(storeResult1->data[iterator]), seed);
	}

	*seed	= -322140768;

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan4(&status, &(storeResult2->data[iterator]), seed);
	}

	printf("\n");

	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; )
	{
		printf(" %1.20f %1.20f", storeResult1->data[iterator], storeResult2->data[iterator]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+1], storeResult2->data[iterator+1]);
		printf(" %1.20f %1.20f", storeResult1->data[iterator+2], storeResult2->data[iterator+2]);
		printf(" %1.20f %1.20f\n", storeResult1->data[iterator+3], storeResult2->data[iterator+3]);

		if ( storeResult2->data[iterator] == storeResult1->data[iterator])
			count++;
		if ( storeResult2->data[iterator+1] == storeResult1->data[iterator+1])
			count++;
		if ( storeResult2->data[iterator+2] == storeResult1->data[iterator+2])
			count++;
		if ( storeResult2->data[iterator+3] == storeResult1->data[iterator+3])
			count++;

		iterator +=4;
	}

	printf("\n");

	*seed = -248673140;
	printf("MESG: Seed Value := %d\n", *seed);
	printf("MESG: Using NRRan3()\n");
	for( iterator = 0; iterator < NRRANDOMTESTC_LENGTH; iterator++ )
	{
		NRRan3(&status, &(storeResult1->data[iterator]), seed);
	}

	printf("\nvalue[0] = %1.030f\n", storeResult1->data[0]);
	if( storeResult1->data[0] == 0.294567555189132690429687500000)
		count++;
	printf("\nvalue[1] = %1.030f\n", storeResult1->data[1]);
	if( storeResult1->data[1] == 0.014104902744293212890625000000)
		count++;
	printf("\nvalue[2] = %1.030f\n", storeResult1->data[2]);
	if( storeResult1->data[2] == 0.089319035410881042480468750000)
		count++;
	printf("\nvalue[3] = %1.030f\n", storeResult1->data[3]);
	if( storeResult1->data[3] == 0.886621534824371337890625000000)
		count++;
	printf("\nvalue[4] = %1.030f\n", storeResult1->data[4]);
	if( storeResult1->data[4] == 0.795643210411071777343750000000)
		count++;

#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
	if( count == (NRRANDOMTESTC_LENGTH)*2 + 5)
	{
		printf("PASS: Same seed reinitializes function properly.\n\n");
	}
	else
	{
		printf("FAIL: Same seed does not reinitialize function properly.\n\n");
		return NRRANDOMTESTC_ESEED;
	}
  }
#endif

	if ( lalDebugLevel > 0 )
	{
		printf("\n");
		printf("MESG:	More option available from command line.\n");
		printf("MESG:	Type NRRandomTest -h for options.\n");
	}

	LALSDestroyVector(&status, &storeResult1);
	LALSDestroyVector(&status, &storeResult2);
	LALSDestroyVector(&status, &temp);
	/* normal exit */
	return NRRANDOMTESTC_ENOM;
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
        optVerbose = NRRANDOMTESTC_TRUE;
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
