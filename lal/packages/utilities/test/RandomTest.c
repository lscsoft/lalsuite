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
#include <config.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>

/**
   \file
   \ingroup Random_h

   \brief Tests the routines in \ref Random.h.  Outputs random numbers to a file.

\heading{Usage}
\code
RandomTest [options]
Options:
  -h         print this message
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
  -o         output random numbers to files
\endcode

\heading{Exit codes}

<table>
<tr><th>Code</th><th>Explanation</th></tr>
<tr><td>0</td><td>Success, normal exit.</td></tr>
<tr><td>1</td><td>Subroutine failed.</td></tr>
</table>

*/
/** \cond DONT_DOXYGEN */

#define CODES_(x) #x
#define CODES(x) CODES_(x)

extern char *optarg;
extern int   optind;

int output     = 0;
extern int lalDebugLevel;
int verbose    = 0;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

int
main (int argc, char *argv[])
{
  const  INT4          numPoints = 999;
  static LALStatus        status;
  static REAL4Vector  *vector;
  static RandomParams *randpar;
  UINT4                i;

  lalDebugLevel = 0;

  /*
   *
   * Parse the command line options.
   *
   */


  ParseOptions (argc, argv);


  /*
   *
   * Allocate memory.
   *
   */


  if (verbose)
  {
    printf ("\n===== Allocate Memory =====\n");
  }

  LALCreateVector (&status, &vector, numPoints);
  TestStatus (&status, CODES(0), 1);


  /*
   *
   * Create random number parameters with seed drawn from clock.
   *
   */


  if (verbose)
  {
    printf ("\n===== Test Random Routines =====\n");
  }

  LALCreateRandomParams (&status, &randpar, 0);
  TestStatus (&status, CODES(0), 1);


  /*
   *
   * Fill vector with uniform random numbers.
   *
   */


  for (i = 0; i < vector->length; ++i)
  {
    LALUniformDeviate (&status, vector->data + i, randpar);
    if (status.statusCode)
    {
      break;
    }
  }
  TestStatus (&status, CODES(0), 1);

  if (output)
  {
    FILE *fp = fopen ("PrintVector.000", "w");
    for (i = 0; i < vector->length; ++i)
    {
      fprintf (fp, "%e\n", vector->data[i]);
    }
    fclose (fp);
  }


  /*
   *
   * Fill vector with Gaussian random numbers.
   *
   */


  LALNormalDeviates (&status, vector, randpar);
  TestStatus (&status, CODES(0), 1);

  if (output)
  {
    FILE *fp = fopen ("PrintVector.001", "w");
    for (i = 0; i < vector->length; ++i)
    {
      fprintf (fp, "%e\n", vector->data[i]);
    }
    fclose (fp);
  }


  /*
   *
   * Check to make sure that correct error codes are generated.
   *
   */


#ifndef LAL_NDEBUG
  if ( ! lalNoDebug )
  {
    if (verbose || lalDebugLevel)
    {
      printf ("\n===== Check Errors =====\n");
    }

    /* non-null pointer error */

    if (verbose)
    {
      printf ("\n----- Non-Null Pointer Error: Code 2\n");
    }

    LALCreateRandomParams (&status, &randpar, 0);
    TestStatus (&status, CODES(RANDOMH_ENNUL), 1);

    /* null pointer error */

    if (verbose)
    {
      printf ("\n----- Null Pointer Error: Code 1 (8 times)\n");
    }

    LALCreateRandomParams (&status, NULL, 0);
    TestStatus (&status, CODES(RANDOMH_ENULL), 1);

    LALDestroyRandomParams (&status, NULL);
    TestStatus (&status, CODES(RANDOMH_ENULL), 1);

    {
      RandomParams *tmp = NULL;
      LALDestroyRandomParams (&status, &tmp);
    }
    TestStatus (&status, CODES(RANDOMH_ENULL), 1);

    LALUniformDeviate (&status, NULL, randpar);
    TestStatus (&status, CODES(RANDOMH_ENULL), 1);

    LALUniformDeviate (&status, vector->data, NULL);
    TestStatus (&status, CODES(RANDOMH_ENULL), 1);

    LALNormalDeviates (&status, NULL, randpar);
    TestStatus (&status, CODES(RANDOMH_ENULL), 1);

    {
      REAL4Vector tmp;
      tmp.length = 10;
      tmp.data   = NULL;
      LALNormalDeviates (&status, &tmp, randpar);
    }
    TestStatus (&status, CODES(RANDOMH_ENULL), 1);

    LALNormalDeviates (&status, vector, NULL);
    TestStatus (&status, CODES(RANDOMH_ENULL), 1);

    /* vector length error */

    if (verbose)
    {
      printf ("\n----- Invalid Length Error: Code 4\n");
    }

    {
      REAL4Vector tmp;
      tmp.length = 0;
      tmp.data   = (REAL4 *)1;
      LALNormalDeviates (&status, &tmp, randpar);
    }
    TestStatus (&status, CODES(RANDOMH_ESIZE), 1);
  }
#endif


  /*
   *
   * Free memory.
   *
   */


  if (verbose || lalDebugLevel)
  {
    printf ("\n===== Clean up and Exit =====\n");
  }

  LALDestroyRandomParams (&status, &randpar);
  TestStatus (&status, CODES(0), 1);

  LALDestroyVector (&status, &vector);
  TestStatus (&status, CODES(0), 1);

  LALCheckMemoryLeaks ();
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
  fprintf (stderr, "  -o         output random numbers to files\n");
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
  FILE *fp;

  while (1)
  {
    int c = -1;

    c = getopt (argc, argv, "hqvd:""o");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'o': /* set output */
        output = 1;
        break;

      case 'd': /* set debug level */
        lalDebugLevel = atoi (optarg);
        break;

      case 'v': /* verbose */
        ++verbose;
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

/** \endcond */
