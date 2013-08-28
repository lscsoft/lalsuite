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

/**
 * \file
 * \ingroup Interpolate_h
 *
 * \brief Tests the routines in \ref Interpolate.h.
 *
 * \heading{Usage}
 * \code
 * InterpolateTest [options]
 * Options:
 * -h         print this message
 * -q         quiet: run silently
 * -v         verbose: print extra information
 * -d level   set lalDebugLevel to level
 * \endcode
 *
 * \heading{Exit codes}
 *
 * <table><tr><th>Code</th><th>Explanation</th></tr>
 * <tr><td>0</td><td>Success, normal exit.</td></tr>
 * <tr><td>1</td><td>Subroutine failed.</td></tr>
 * </table>
 *
 */

/** \cond DONT_DOXYGEN */
#include <config.h>

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/Interpolate.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

extern char *optarg;
extern int   optind;

int   verbose    = 0;

static void Usage (const char *program, int exitflag);

static void ParseOptions (int argc, char *argv[]);

static void TestStatus (LALStatus *status, const char *expectCodes, int exitCode);

int main (int argc, char *argv[])
{
  enum {ArraySize = 10};

  static LALStatus       status;
  static REAL4Vector *x;
  static REAL4Vector *y;
  static REAL8Vector *xx;
  static REAL8Vector *yy;

  SInterpolateOut sintout;
  SInterpolatePar sintpar;
  DInterpolateOut dintout;
  DInterpolatePar dintpar;

  int i;

  ParseOptions (argc, argv);

  LALSCreateVector (&status, &x, ArraySize);
  TestStatus (&status, CODES(0), 1);
  LALSCreateVector (&status, &y, ArraySize);
  TestStatus (&status, CODES(0), 1);
  LALDCreateVector (&status, &xx, ArraySize);
  TestStatus (&status, CODES(0), 1);
  LALDCreateVector (&status, &yy, ArraySize);
  TestStatus (&status, CODES(0), 1);

  if ( verbose )
    printf ("Initial data:\n");
  if ( verbose )
    printf ("y = P(x) = 7 - 8x + 2x^2 + 2x^3 - x^4\n");
  if ( verbose )
    printf ("-----------------\n");
  if ( verbose )
    printf ("x\ty\n");
  if ( verbose )
    printf ("=================\n");
  for (i = 0; i < ArraySize; ++i)
  {
    REAL4 xi = x->data[i] = xx->data[i] = i*i;
    y->data[i] = yy->data[i] = 7 + xi*(-8 + xi*(2 + xi*(2 - xi)));
    if ( verbose )
      printf ("%.0f\t%.0f\n", x->data[i], y->data[i]);
  }
  if ( verbose )
    printf ("-----------------\n");

  if ( verbose )
    printf ("\nInterpolate to x = 0.3:\n");
  if ( verbose )
    printf ("---------------------------------\n");
  if ( verbose )
    printf ("order\ty\t\tdy\n");
  if ( verbose )
    printf ("=================================\n");
  sintpar.x = x->data;
  sintpar.y = y->data;
  dintpar.x = xx->data;
  dintpar.y = yy->data;
  for (i = 2; i < ArraySize; ++i)
  {
    sintpar.n = i;
    dintpar.n = i;
    LALSPolynomialInterpolation (&status, &sintout, 0.3, &sintpar);
    TestStatus (&status, CODES(0), 1);
    LALDPolynomialInterpolation (&status, &dintout, 0.3, &dintpar);
    TestStatus (&status, CODES(0), 1);
    if ( verbose )
      printf ("%d\t%f\t%f\t%f\t%f\n", i - 1,
          sintout.y, sintout.dy, dintout.y, dintout.dy);
  }
  if ( verbose )
    printf ("---------------------------------\n");

  if ( verbose )
    printf ("\nExtrapolate to x = -0.3:\n");
  if ( verbose )
    printf ("---------------------------------\n");
  if ( verbose )
    printf ("order\ty\t\tdy\n");
  if ( verbose )
    printf ("=================================\n");
  sintpar.x = x->data;
  sintpar.y = y->data;
  dintpar.x = xx->data;
  dintpar.y = yy->data;
  for (i = 2; i < ArraySize; ++i)
  {
    sintpar.n = i;
    dintpar.n = i;
    LALSPolynomialInterpolation (&status, &sintout, -0.3, &sintpar);
    TestStatus (&status, CODES(0), 1);
    LALDPolynomialInterpolation (&status, &dintout, -0.3, &dintpar);
    TestStatus (&status, CODES(0), 1);
    if ( verbose )
      printf ("%d\t%f\t%f\t%f\t%f\n", i - 1,
          sintout.y, sintout.dy, dintout.y, dintout.dy);
  }
  if ( verbose )
    printf ("---------------------------------\n");


  LALSDestroyVector (&status, &x);
  TestStatus (&status, CODES(0), 1);
  LALSDestroyVector (&status, &y);
  TestStatus (&status, CODES(0), 1);
  LALDDestroyVector (&status, &xx);
  TestStatus (&status, CODES(0), 1);
  LALDDestroyVector (&status, &yy);
  TestStatus (&status, CODES(0), 1);

  LALCheckMemoryLeaks ();

  LALSCreateVector (&status, &x, ArraySize);
  TestStatus (&status, CODES(0), 1);
  LALSCreateVector (&status, &y, ArraySize);
  TestStatus (&status, CODES(0), 1);
  LALDCreateVector (&status, &xx, ArraySize);
  TestStatus (&status, CODES(0), 1);
  LALDCreateVector (&status, &yy, ArraySize);
  TestStatus (&status, CODES(0), 1);


#ifndef LAL_NDEBUG

  if ( ! lalNoDebug )
  {
    if ( verbose )
      printf ("\nCheck error conditions:\n");

    if ( verbose )
      printf ("\nNull pointer:\r");
    LALSPolynomialInterpolation (&status, NULL, -0.3, &sintpar);
    TestStatus (&status, CODES(INTERPOLATEH_ENULL), 1);
    LALDPolynomialInterpolation (&status, NULL, -0.3, &dintpar);
    TestStatus (&status, CODES(INTERPOLATEH_ENULL), 1);
    if ( verbose )
      printf ("Null pointer check passed.\n");

    if ( verbose )
      printf ("\nNull pointer:\r");
    LALSPolynomialInterpolation (&status, &sintout, -0.3, NULL);
    TestStatus (&status, CODES(INTERPOLATEH_ENULL), 1);
    LALDPolynomialInterpolation (&status, &dintout, -0.3, NULL);
    TestStatus (&status, CODES(INTERPOLATEH_ENULL), 1);
    if ( verbose )
      printf ("Null pointer check passed.\n");

    sintpar.x = NULL;
    dintpar.x = NULL;
    if ( verbose )
      printf ("\nNull pointer:\r");
    LALSPolynomialInterpolation (&status, &sintout, -0.3, &sintpar);
    TestStatus (&status, CODES(INTERPOLATEH_ENULL), 1);
    LALDPolynomialInterpolation (&status, &dintout, -0.3, &dintpar);
    TestStatus (&status, CODES(INTERPOLATEH_ENULL), 1);
    if ( verbose )
      printf ("Null pointer check passed.\n");

    sintpar.x = x->data;
    sintpar.y = NULL;
    dintpar.x = xx->data;
    dintpar.y = NULL;
    if ( verbose )
      printf ("\nNull pointer:\r");
    LALSPolynomialInterpolation (&status, &sintout, -0.3, &sintpar);
    TestStatus (&status, CODES(INTERPOLATEH_ENULL), 1);
    LALDPolynomialInterpolation (&status, &dintout, -0.3, &dintpar);
    TestStatus (&status, CODES(INTERPOLATEH_ENULL), 1);
    if ( verbose )
      printf ("Null pointer check passed.\n");

    sintpar.y = y->data;
    sintpar.n = 1;
    dintpar.y = yy->data;
    dintpar.n = 1;
    if ( verbose )
      printf ("\nInvalid size:\r");
    LALSPolynomialInterpolation (&status, &sintout, -0.3, &sintpar);
    TestStatus (&status, CODES(INTERPOLATEH_ESIZE), 1);
    dintout.dy = XLALREAL8PolynomialInterpolation (&(dintout.y), -0.3, dintpar.y, dintpar.x, dintpar.n);
    if (xlalErrno == XLAL_ESIZE)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Invalid size check passed.\n");

    x->data[1]  = x->data[0]  = 2;
    xx->data[1] = xx->data[0] = 2;
    sintpar.n = 3;
    dintpar.n = 3;
    if ( verbose )
      printf ("\nZero divide:\r");
    LALSPolynomialInterpolation (&status, &sintout, -0.3, &sintpar);
    TestStatus (&status, CODES(INTERPOLATEH_EZERO), 1);
    dintout.dy = XLALREAL8PolynomialInterpolation (&(dintout.y), -0.3, dintpar.y, dintpar.x, dintpar.n);
    if (xlalErrno == XLAL_EFPDIV0)
      xlalErrno = 0;
    else
      abort();
    if ( verbose )
      printf ("Zero divide check passed.\n");
  }
#endif


  LALSDestroyVector (&status, &x);
  TestStatus (&status, CODES(0), 1);
  LALSDestroyVector (&status, &y);
  TestStatus (&status, CODES(0), 1);
  LALDDestroyVector (&status, &xx);
  TestStatus (&status, CODES(0), 1);
  LALDDestroyVector (&status, &yy);
  TestStatus (&status, CODES(0), 1);

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

    c = getopt (argc, argv, "hqvd:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'd': /* set debug level */
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently */
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
