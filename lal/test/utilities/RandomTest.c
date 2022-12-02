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
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/
#include <config.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/AVFactories.h>
#include <lal/Random.h>
#include <lal/LALString.h>

/**
 * \file
 * \ingroup Random_h
 *
 * \brief Tests the routines in \ref Random.h.  Outputs random numbers to a file.
 *
 * ### Usage ###
 *
 * \code
 * RandomTest [options]
 * Options:
 * -h         print this message
 * -q         quiet: run silently
 * -v         verbose: print extra information
 * -d level   set lalDebugLevel to level
 * -o         output random numbers to files
 * \endcode
 *
 * ### Exit codes ###
 *
 * <table>
 * <tr><th>Code</th><th>Explanation</th></tr>
 * <tr><td>0</td><td>Success, normal exit.</td></tr>
 * <tr><td>1</td><td>Subroutine failed.</td></tr>
 * </table>
 *
 */
/** \cond DONT_DOXYGEN */

#define CODES_(x) #x
#define CODES(x) CODES_(x)

int output     = 0;
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

  randpar = XLALCreateRandomParams (0);
  if (!randpar)
    exit (1);


  /*
   *
   * Fill vector with uniform random numbers.
   *
   */


  for (i = 0; i < vector->length; ++i)
  {
    vector->data[i] = XLALUniformDeviate (randpar);
    if (XLALIsREAL4FailNaN(vector->data[i]))
    {
      exit (1);
    }
  }

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


  if (XLALNormalDeviates (vector, randpar))
    exit (1);

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


    if (verbose || lalDebugLevel)
    {
      printf ("\n===== Check Errors =====\n");
    }

    /* vector length error */

    if (verbose)
    {
      printf ("\n----- Invalid Length Error: Code 4\n");
    }

    {
      REAL4Vector tmp;
      tmp.length = 0;
      tmp.data   = (REAL4 *)1;
      if (!XLALNormalDeviates (&tmp, randpar))
        exit (1);
      XLALClearErrno();
    }


  /*
   *
   * Free memory.
   *
   */


  if (verbose || lalDebugLevel)
  {
    printf ("\n===== Clean up and Exit =====\n");
  }

  XLALDestroyRandomParams (randpar);

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

  if (XLALStringCopy(str, ignored, sizeof(str)))
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
      if (status->statusCode == atoi (str))
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

    c = LALgetopt (argc, argv, "hqvd:""o");
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

  if (LALoptind < argc)
  {
    Usage (argv[0], 1);
  }

  return;
}

/** \endcond */
