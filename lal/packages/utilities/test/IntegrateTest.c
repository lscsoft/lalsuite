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
#include "Integrate.h"

#define CODES_(x) #x
#define CODES(x) CODES_(x)

NRCSID (MAIN, "$Id$");


/*
 *
 * These functions are used as integrands for the test integrations.
 *
 */


void f1 (Status *s, REAL4 *y, REAL4 x, void *p)
{
  REAL4  x2 = x*x;
  REAL4  x4 = x2*x2;
  INT4  *n;
  INITSTATUS (s, MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));  
  *y = x4*log(x + sqrt(x2 + 1));
  RETURN (s);
}

void f2 (Status *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/(x*x*x);
  RETURN (s);
}

void f3 (Status *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = exp(-x*x/2);
  RETURN (s);
}

void f4 (Status *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/sqrt(x);
  RETURN (s);
}

void f5 (Status *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = x + 1/sqrt(5 - x);
  RETURN (s);
}

void g (Status *s, REAL4 *z, REAL4 x, void *p)
{
  REAL4 y;
  INITSTATUS (s, MAIN);
  ASSERT (p, s, 1, "Null pointer");
  y  = *((REAL4 *)p);
  *z = exp(-(x*x + y*y)/2);
  RETURN (s);
}

void h (Status *s, REAL4 *z, REAL4 y, void *p)
{
  IntegrateIn intinp;
  INITSTATUS (s, MAIN);
  ATTATCHSTATUSPTR (s);
  ASSERT (!p, s, 2, "Non-null pointer");
  intinp.function = g;
  intinp.xmin     = 0;
  intinp.xmax     = sqrt(y);
  intinp.type     = ClosedInterval;
  RombergIntegrate (s->statusPtr, z, &intinp, &y);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
}


/*
 *
 * This function produces random numbers.  Integration of this function should
 * not converge.  Make this routine fast... no status handling!
 *
 */
void bad (Status *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n = (INT4 *)p;
  *y = *n = 1664525L*(*n) + 1013904223L;
}



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
  const REAL4   epsilon = 1e-6;
  static Status status;
  IntegrateIn   intinp;
  REAL4         result;
  long double   expect;
  INT4          count;

  ParseOptions (argc, argv);


  /*
   *
   * Test 1: Integrate a regular function over a closed interval.
   *
   */


  printf ("Test 1:"
          " Integrate a regular function over a closed interval.\n");

  intinp.function = f1;
  intinp.xmin     = 0;
  intinp.xmax     = 2;
  intinp.type     = ClosedInterval;

  count  = 0;
  expect = 8.153364119811650205;
  RombergIntegrate (&status, &result, &intinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", result);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(result - expect) > epsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  /*
   *
   * Test 2: Integrate to infinity a function with power-law fall-off.
   *
   */


  printf ("\nTest 2:"
          " Integrate to infinity a function with power-law fall-off.\n");
  intinp.function = f2;
  intinp.xmin     = 10;
  intinp.xmax     = 1e30;
  intinp.type     = InfiniteDomainPow;

  count  = 0;
  expect = 1.0/200.0;
  RombergIntegrate (&status, &result, &intinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", result);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(result - expect) > epsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  /*
   *
   * Test 3: Integrate to infinity a function that falls off exponentially.
   *
   */


  printf ("\nTest 3:"
          " Integrate to infinity a function that falls off exponentially.");
  intinp.function = f3;
  intinp.xmin     = 2;
  intinp.xmax     = 1e30;
  intinp.type     = InfiniteDomainExp;

  count  = 0;
  expect = 0.0570261239928920483;
  RombergIntegrate (&status, &result, &intinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", result);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(result - expect) > epsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  /*
   *
   * Test 4: Integrate an integrable singularity at the lower limit.
   *
   */


  printf ("\nTest 4:"
          " Integrate an integrable singularity at the lower limit.\n");
  intinp.function = f4;
  intinp.xmin     = 0;
  intinp.xmax     = 1;
  intinp.type     = SingularLowerLimit;

  count  = 0;
  expect = 2;
  RombergIntegrate (&status, &result, &intinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", result);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(result - expect) > epsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  /*
   *
   * Test 5: Integrate an integrable singularity at the upper limit.
   *
   */


  printf ("\nTest 5:"
          " Integrate an integrable singularity at the upper limit.\n");
  intinp.function = f5;
  intinp.xmin     = 4;
  intinp.xmax     = 5;
  intinp.type     = SingularUpperLimit;

  count  = 0;
  expect = 6.5;
  RombergIntegrate (&status, &result, &intinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", result);
  printf ("expect: %.15Lf\n", expect);
  /* this doesn't work so well: multiply tolerance by factor of three */
  if (fabs(result - expect) > 3*epsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }


  /*
   *
   * Test 6: Two-dimensional integral.
   *
   */


  printf ("\nTest 6: Two-dimensional integral.\n");
  intinp.function = h;
  intinp.xmin     = 0;
  intinp.xmax     = 10;
  intinp.type     = OpenInterval;

  expect = 0.88274109326014810823;
  RombergIntegrate (&status, &result, &intinp, NULL);
  TestStatus (&status, CODES(0), 1);
  printf ("result: %.15f\n", result);
  printf ("expect: %.15Lf\n", expect);
  /* integral isn't very accurate because we needed to use an open interval */
  printf ("error:  %.2f%%\n", 100*fabs(result - expect)/fabs(expect));


  LALCheckMemoryLeaks ();


  /*
   *
   * Check error conditions.
   *
   */


  printf ("\nChecking error conditions:\n");

  printf ("\nNull pointer:\r");
  RombergIntegrate (&status, NULL, &intinp, &count);
  TestStatus (&status, CODES(INTEGRATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  printf ("\nNull pointer:\r");
  RombergIntegrate (&status, &result, NULL, &count);
  TestStatus (&status, CODES(INTEGRATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  printf ("\nNull pointer:\r");
  intinp.function = NULL;
  intinp.xmin     = 0;
  intinp.xmax     = 2;
  intinp.type     = ClosedInterval;
  RombergIntegrate (&status, &result, &intinp, &count);
  TestStatus (&status, CODES(INTEGRATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  printf ("\nInvalid domain:\r");
  intinp.function = f1;
  intinp.xmin     = 0;
  intinp.xmax     = 0;
  intinp.type     = ClosedInterval;
  RombergIntegrate (&status, &result, &intinp, &count);
  TestStatus (&status, CODES(INTEGRATE_EIDOM), 1);
  printf ("Invalid domain check passed.\n");

  printf ("\nUnknown integral type:\r");
  intinp.function = f1;
  intinp.xmin     = 0;
  intinp.xmax     = 2;
  intinp.type     = 999;
  RombergIntegrate (&status, &result, &intinp, &count);
  TestStatus (&status, CODES(INTEGRATE_ETYPE), 1);
  printf ("Unknown integral type check passed.\n");

  printf ("\nMaximum iterations exceeded:\r");
  intinp.function = bad;  /* bad is a quick random number generator */
  intinp.xmin     = 0;
  intinp.xmax     = 2;
  intinp.type     = ClosedInterval;
  count           = 13;   /* count is now used as a random number seed */
  RombergIntegrate (&status, &result, &intinp, &count);
  TestStatus (&status, CODES(INTEGRATE_EMXIT), 1);
  printf ("Maximum iterations exceeded check passed.\n");

  printf ("\nRecursive error:\r");
  intinp.function = f1;
  intinp.xmin     = 0;
  intinp.xmax     = 2;
  intinp.type     = ClosedInterval;
  RombergIntegrate (&status, &result, &intinp, NULL);
  TestStatus (&status, CODES(-1), 1);
  printf ("Recursive error check passed.\n");
  ClearStatus (&status);

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

