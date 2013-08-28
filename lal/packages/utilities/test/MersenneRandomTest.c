/*
*  Copyright (C) 2007 Jolien Creighton
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

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>

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

#include <lal/Random.h>

/**
 * \author Tibbits, M. M.
 * \file
 * \ingroup Random_h
 *
 * \brief A program to test <tt>LALMersenneRandom()</tt> and <tt>LALMersenneRandomVector()</tt>
 *
 * ### Usage ###
 *
 * \code
 * ./MersenneRandomTest [options]
 * Options:
 * -h             print usage message
 * -q             quiet: run silently
 * -v             verbose: print extra information
 * -d level       set lalDebugLevel to level
 * \endcode
 *
 * This program tests the function
 * <tt>LALMersenneRandom()</tt>,which generates a random
 * number based on the Mersenne Twister algorithm.
 *
 * First, it tests that the correct error codes
 * are generated for the following error conditions passed
 * to the function LALMersenneRandom() (tests in \e italics
 * are not performed if \c LAL_NEDEBUG is set, as the
 * corresponding checks in the code are made using the ASSERT macro):
 *
 * <ul>
 * <li> <em>null pointer to output structure</em></li>
 * <li> <em>null pointer to params structure</em></li>
 * <li> <em>params not initialized</em></li>
 * </ul>
 *
 * Second, it tests that the correct error codes
 * are generated for the following error conditions passed
 * to the function LALMersenneRandomVector() (tests in \e italics
 * are not performed if \c LAL_NEDEBUG is set, as the
 * corresponding checks in the code are made using the ASSERT macro):
 *
 * <ul>
 * <li> <em>null pointer to output structure</em></li>
 * <li> <em>null pointer to params structure</em></li>
 * <li> <em>params not initialized</em></li>
 * <li> <em>outputVector-\f$>\f$length = 0</em></li>
 * </ul>
 *
 * Third, it verifies the output of the generator
 * for each of the following simple test cases:
 * <ol>
 * <li> given a certain seed, does the output match the expected?</li>
 * <li> does calling the function again reinitialize it to the new seed properly?</li>
 * <li> does it create a vector of random numbers correctly?</li>
 * </ol>
 *
 * For each successful test
 * (both of these valid data and the invalid ones described above), it
 * prints "\c PASS" to standard output; if a test fails, it
 * prints "\c FAIL".
 *
 * ### Notes ###
 *
 * <ul>
 * <li>Vector must come in allocated</li>
 * <li>params must be initialized before calls can be made.</li>
 * </ul>
 *
 */

/** \cond DONT_DOXYGEN */

/* bogus type */
struct
tagMTRandomParams
{
        UINT4           seed;
        INT2            initialized;
        void            *priv;
};


/*  constants  */
#define MERSENNERANDOMTESTC_TRUE     1
#define MERSENNERANDOMTESTC_FALSE    0

extern char	*optarg;
extern int	optind;



/*  Setting variables to parse command line flags  */
BOOLEAN	optVerbose	= MERSENNERANDOMTESTC_FALSE;
UINT4	optLength	= 0;

INT4
CheckStatus(LALStatus *status, const INT4 code, const CHAR *message,
            const INT4 exitcode, const CHAR *error);


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


/**\name Error Codes */ /*@{*/
#define	MERSENNERANDOMTESTC_ENOM     0
#define	MERSENNERANDOMTESTC_EARG     1
#define	MERSENNERANDOMTESTC_ECHK     2
#define	MERSENNERANDOMTESTC_EFLS     3
#define	MERSENNERANDOMTESTC_EUSE     4
#define	MERSENNERANDOMTESTC_ENULL    5
#define	MERSENNERANDOMTESTC_EALOC    6
#define	MERSENNERANDOMTESTC_ENMM     7
#define	MERSENNERANDOMTESTC_MSGENOM  "Nominal exit"
#define	MERSENNERANDOMTESTC_MSGEARG  "Error parsing command-line arguments"
#define	MERSENNERANDOMTESTC_MSGECHK  "Error checking failed to catch bad data"
#define	MERSENNERANDOMTESTC_MSGEFLS  "Incorrect answer for valid data"
#define	MERSENNERANDOMTESTC_MSGEUSE  "Bad user-entered data"
#define	MERSENNERANDOMTESTC_MSGENULL "Null Pointer."
#define MERSENNERANDOMTESTC_MSGEALOC "Memory Allocation Error"
#define	MERSENNERANDOMTESTC_MSGENMM  "Randeom Number Mismatch"
/*@}*/



int main( int argc, char *argv[] )
{

	static	LALStatus	status;

	/* Variable declarations */
	REAL8Vector	*outputVector;
	REAL8Vector	*testVector1;
	REAL8Vector	*testVector2;
	REAL8Vector	*testVector3;
	MTRandomParams	*params;
	REAL8		output;
	INT8		iterator;
	INT4		count;

	ParseOptions( argc, argv );

	/*  Initialize Variables  */
	outputVector	= NULL;
	testVector1	= NULL;
	testVector2	= NULL;
	testVector3	= NULL;

	params		= NULL;
	count		= 0;

	LALDCreateVector(&status, &outputVector, 1000);
	LALDCreateVector(&status, &testVector1, 1000);
	LALDCreateVector(&status, &testVector2, 500);
	LALDCreateVector(&status, &testVector3, 500);

	/*  Initialize Parameter structure  */
	LALCreateMTRandomParams(&status, 4357, &params);

	printf("\n\nMESG: %s \n","$Id$");

#ifndef LAL_NDEBUG
	INT4		code;
	REAL8		*nullOutput = NULL;
	MTRandomParams	*nullParams	= NULL;
	REAL8Vector	*nullVector	= NULL;

  if ( ! lalNoDebug )
  {
	/*******************  Test LALMersenneRandom()  *********************/
	/* test behavior for null pointer to output structure */
	LALMersenneRandom(&status, nullOutput, params);
	if ( ( code = CheckStatus(&status, RANDOMH_ENULL, RANDOMH_MSGENULL,
					MERSENNERANDOMTESTC_ECHK, MERSENNERANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to output structure results in error:\n");
	printf("       \"%s\"\n", RANDOMH_MSGENULL);

	/* test behavior for null pointer to params structure */
	LALMersenneRandom(&status, &output, nullParams);
	if ( ( code = CheckStatus(&status, RANDOMH_ENULL, RANDOMH_MSGENULL,
					MERSENNERANDOMTESTC_ECHK, MERSENNERANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to params structure results in error:\n");
	printf("       \"%s\"\n", RANDOMH_MSGENULL);

	/* test behavior for non-initialized params structure  */
	params->initialized = 0;
	LALMersenneRandom(&status, &output, params);
	if ( ( code = CheckStatus(&status, RANDOMH_EINIT, RANDOMH_MSGEINIT,
					MERSENNERANDOMTESTC_ECHK, MERSENNERANDOMTESTC_MSGECHK)) )
	{
		return code;
	}

	printf("\nPASS: non-initialized params structure results in error:\n");
	printf("       \"%s\"\n", RANDOMH_MSGEINIT);
	params->initialized = 1;


	/*******************  Test LALMersenneRandomVector()  *********************/
	/* test behavior for null pointer to output structure */
	LALMersenneRandomVector(&status, nullVector, params);
	if ( ( code = CheckStatus(&status, RANDOMH_ENULL, RANDOMH_MSGENULL,
					MERSENNERANDOMTESTC_ECHK, MERSENNERANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to output structure results in error:\n");
	printf("       \"%s\"\n", RANDOMH_MSGENULL);

	/* test behavior for null pointer to params structure */
	LALMersenneRandomVector(&status, outputVector, nullParams);
	if ( ( code = CheckStatus(&status, RANDOMH_ENULL, RANDOMH_MSGENULL,
					MERSENNERANDOMTESTC_ECHK, MERSENNERANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: null pointer to params structure results in error:\n");
	printf("       \"%s\"\n", RANDOMH_MSGENULL);

	/* test behavior for non-initialized params structure  */
	params->initialized = 0;
	LALMersenneRandomVector(&status, outputVector, params);
	if ( ( code = CheckStatus(&status, RANDOMH_EINIT, RANDOMH_MSGEINIT,
					MERSENNERANDOMTESTC_ECHK, MERSENNERANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: non-initialized params structure results in error:\n");
	printf("       \"%s\"\n", RANDOMH_MSGEINIT);
	params->initialized = 1;

	/* test behavior for outputVector length = 0  */
	UINT8		tempLength = outputVector->length;
	outputVector->length = 0;
	LALMersenneRandomVector(&status, outputVector, params);
	if ( ( code = CheckStatus(&status, RANDOMH_EZERO, RANDOMH_MSGEZERO,
					MERSENNERANDOMTESTC_ECHK, MERSENNERANDOMTESTC_MSGECHK)) )
	{
		return code;
	}
	printf("\nPASS: outputVector length = 0 results in error:\n");
	printf("       \"%s\"\n", RANDOMH_MSGEZERO);
	outputVector->length = tempLength;
  }
#endif

	for (iterator = 0; iterator < 1000; iterator++)
	{
		LALMersenneRandom(&status, &output, params);
		printf("%10.8f ", output);
		if (iterator % 8 == 7)
			printf("\n");
	}

	LALDestroyMTRandomParams(&status, &params);

	LALCreateMTRandomParams(&status, 4357, &params);

	LALMersenneRandomVector(&status, testVector1, params);

	LALDestroyMTRandomParams(&status, &params);

	LALCreateMTRandomParams(&status, 4357, &params);

	LALMersenneRandomVector(&status, testVector2, params);

	for(iterator = 0; iterator < 500; iterator++)
	{
		if (testVector1->data[iterator] == testVector2->data[iterator])
			count++;
	}

	LALMersenneRandomVector(&status, testVector3, params);

	for(; iterator < 1000; iterator++)
	{
		if (testVector1->data[iterator] == testVector3->data[iterator-500])
			count++;
	}

	printf("\ncount := %d\n",count);

	LALDestroyMTRandomParams(&status, &params);

	printf("\n");
	printf("MESG:	More option available from command line.\n");
	printf("MESG:	Type MersenneRandomTest -h for options.\n");

	/* normal exit */

	if ( count == 1000 )
	{
		return MERSENNERANDOMTESTC_ENOM;
	}
	else
	{
		return MERSENNERANDOMTESTC_ENMM;
	}
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
  FILE *fp;

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
        break;

      case 'v': /* optVerbose */
        optVerbose = MERSENNERANDOMTESTC_TRUE;
        break;

      case 'q': /* quiet: run silently (ignore error messages) */
        fp = freopen ("/dev/null", "w", stderr);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
        fp = freopen ("/dev/null", "w", stdout);
        if (fp == NULL)
        {
          fprintf(stderr, "Error: Unable to open /dev/null\n");
          exit(1);
        }
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

INT4
CheckStatus(LALStatus *status, const INT4 code, const CHAR *message,
            const INT4 exitcode, const CHAR *error)
{

  if (optVerbose)
  {
    REPORTSTATUS (status);
  }
  if (status->statusCode!= code)
  {
    if (code) printf ("  FAIL: did not recognize \"%s\"\n", message);
    if (optVerbose) printf("Exiting with error: %s\n", error);
    return(exitcode);
  }
  else if (code && strcmp(message, status->statusDescription))
  {
    printf("  FAIL: incorrect error message \"%s\" not \"%s\"\n",
           status->statusDescription, message);
    if (optVerbose) printf("Exiting with error: %s\n", error);
    return(exitcode);
  }
  return 0;
}

/** \endcond */
