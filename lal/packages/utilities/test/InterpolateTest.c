/*----------------------------------------------------------------------- 
 * 
 * File Name: IntegrateTest.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "LALStdlib.h"
#include "AVFactories.h"
#include "Interpolate.h"

#define CODES_(x) #x
#define CODES(x) CODES_(x)

NRCSID (MAIN, "$Id$");

extern char *optarg;
extern int   optind;

int   debuglevel = 0;
int   verbose    = 0;

static void Usage (const char *program, int exitflag);

static void ParseOptions (int argc, char *argv[]);

static void TestStatus (Status *status, const char *expectCodes, int exitCode);

static void ClearStatus (Status *status);

int main (int argc, char *argv[])
{
  enum {ArraySize = 10};

  static Status       status;
  static REAL4Vector *x;
  static REAL4Vector *y;

  InterpolateOut intout;
  InterpolatePar intpar;

  int i;

  ParseOptions (argc, argv);

  CreateVector (&status, &x, ArraySize);
  TestStatus (&status, CODES(0), 1);
  CreateVector (&status, &y, ArraySize);
  TestStatus (&status, CODES(0), 1);

  printf ("Initial data:\n");
  printf ("y = P(x) = 7 - 8x + 2x^2 + 2x^3 - x^4\n");
  printf ("-----------------\n");
  printf ("x\ty\n");
  printf ("=================\n");
  for (i = 0; i < ArraySize; ++i)
  {
    REAL4 xi = x->data[i] = i*i;
    y->data[i] = 7 + xi*(-8 + xi*(2 + xi*(2 - xi)));
    printf ("%.0f\t%.0f\n", x->data[i], y->data[i]);
  }
  printf ("-----------------\n");

  printf ("\nInterpolate to x = 0.3:\n");
  printf ("---------------------------------\n");
  printf ("order\ty\t\tdy\n");
  printf ("=================================\n");
  intpar.x = x->data;
  intpar.y = y->data;
  for (i = 2; i < ArraySize; ++i)
  {
    intpar.n = i;
    PolynomialInterpolation (&status, &intout, 0.3, &intpar);
    TestStatus (&status, CODES(0), 1);
    printf ("%d\t%f\t%f\n", i - 1, intout.y, intout.dy);
  }
  printf ("---------------------------------\n");

  printf ("\nExtrapolate to x = -0.3:\n");
  printf ("---------------------------------\n");
  printf ("order\ty\t\tdy\n");
  printf ("=================================\n");
  intpar.x = x->data;
  intpar.y = y->data;
  for (i = 2; i < ArraySize; ++i)
  {
    intpar.n = i;
    PolynomialInterpolation (&status, &intout, -0.3, &intpar);
    TestStatus (&status, CODES(0), 1);
    printf ("%d\t%f\t%f\n", i - 1, intout.y, intout.dy);
  }
  printf ("---------------------------------\n");


  DestroyVector (&status, &x);
  TestStatus (&status, CODES(0), 1);
  DestroyVector (&status, &y);
  TestStatus (&status, CODES(0), 1);

  LALCheckMemoryLeaks ();

  CreateVector (&status, &x, ArraySize);
  TestStatus (&status, CODES(0), 1);
  CreateVector (&status, &y, ArraySize);
  TestStatus (&status, CODES(0), 1);


  printf ("\nCheck error conditions:\n");

  printf ("\nNull pointer:\r");
  PolynomialInterpolation (&status, NULL, -0.3, &intpar);
  TestStatus (&status, CODES(INTERPOLATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  printf ("\nNull pointer:\r");
  PolynomialInterpolation (&status, &intout, -0.3, NULL);
  TestStatus (&status, CODES(INTERPOLATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  intpar.x = NULL;
  printf ("\nNull pointer:\r");
  PolynomialInterpolation (&status, &intout, -0.3, &intpar);
  TestStatus (&status, CODES(INTERPOLATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  intpar.x = x->data;
  intpar.y = NULL;
  printf ("\nNull pointer:\r");
  PolynomialInterpolation (&status, &intout, -0.3, &intpar);
  TestStatus (&status, CODES(INTERPOLATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  intpar.y = y->data;
  intpar.n = 1;
  printf ("\nInvalid size:\r");
  PolynomialInterpolation (&status, &intout, -0.3, &intpar);
  TestStatus (&status, CODES(INTERPOLATE_ESIZE), 1);
  printf ("Invalid size check passed.\n");

  x->data[1] = x->data[0];
  intpar.n = 3;
  printf ("\nZero divide:\r");
  PolynomialInterpolation (&status, &intout, -0.3, &intpar);
  TestStatus (&status, CODES(INTERPOLATE_EZERO), 1);
  printf ("Zero divide check passed.\n");


  DestroyVector (&status, &x);
  TestStatus (&status, CODES(0), 1);
  DestroyVector (&status, &y);
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
TestStatus (Status *status, const char *ignored, int exitcode)
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
ClearStatus (Status *status)
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
  fprintf (stderr, "  -d level   set debuglevel to level\n");
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
        debuglevel = atoi (optarg);
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'q': /* quiet: run silently */
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


