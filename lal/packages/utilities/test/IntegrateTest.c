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

#include "LALConfig.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

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
  INITSTATUS (s, "f1", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));  
  *y = x4*log(x + sqrt(x2 + 1));
  RETURN (s);
}

void ff1 (Status *s, REAL8 *y, REAL8 x, void *p)
{
  REAL8  x2 = x*x;
  REAL8  x4 = x2*x2;
  INT4  *n;
  INITSTATUS (s, "ff1", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));  
  *y = x4*log(x + sqrt(x2 + 1));
  RETURN (s);
}

void f2 (Status *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, "f2", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/(x*x*x);
  RETURN (s);
}

void ff2 (Status *s, REAL8 *y, REAL8 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, "ff2", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/(x*x*x);
  RETURN (s);
}

void f3 (Status *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, "f3", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = exp(-x*x/2);
  RETURN (s);
}

void ff3 (Status *s, REAL8 *y, REAL8 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, "ff3", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = exp(-x*x/2);
  RETURN (s);
}

void f4 (Status *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, "f4", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/sqrt(x);
  RETURN (s);
}

void ff4 (Status *s, REAL8 *y, REAL8 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, "ff4", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = 1/sqrt(x);
  RETURN (s);
}

void f5 (Status *s, REAL4 *y, REAL4 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, "f5", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = x + 1/sqrt(5 - x);
  RETURN (s);
}

void ff5 (Status *s, REAL8 *y, REAL8 x, void *p)
{
  INT4 *n;
  INITSTATUS (s, "ff5", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  ++(*(n = (INT4 *)p));
  *y = x + 1/sqrt(5 - x);
  RETURN (s);
}

void g (Status *s, REAL4 *z, REAL4 x, void *p)
{
  REAL4 y;
  INITSTATUS (s, "g", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  y  = *((REAL4 *)p);
  *z = exp(-(x*x + y*y)/2);
  RETURN (s);
}

void gg (Status *s, REAL8 *z, REAL8 x, void *p)
{
  REAL8 y;
  INITSTATUS (s, "gg", MAIN);
  ASSERT (p, s, 1, "Null pointer");
  y  = *((REAL8 *)p);
  *z = exp(-(x*x + y*y)/2);
  RETURN (s);
}

void h (Status *s, REAL4 *z, REAL4 y, void *p)
{
  SIntegrateIn intinp;
  INITSTATUS (s, "h", MAIN);
  ATTATCHSTATUSPTR (s);
  ASSERT (!p, s, 2, "Non-null pointer");
  intinp.function = g;
  intinp.xmin     = 0;
  intinp.xmax     = sqrt(y);
  intinp.type     = ClosedInterval;
  SRombergIntegrate (s->statusPtr, z, &intinp, &y);
  CHECKSTATUSPTR (s);
  DETATCHSTATUSPTR (s);
  RETURN (s);
}

void hh (Status *s, REAL8 *z, REAL8 y, void *p)
{
  DIntegrateIn intinp;
  INITSTATUS (s, "hh", MAIN);
  ATTATCHSTATUSPTR (s);
  ASSERT (!p, s, 2, "Non-null pointer");
  intinp.function = gg;
  intinp.xmin     = 0;
  intinp.xmax     = sqrt(y);
  intinp.type     = ClosedInterval;
  DRombergIntegrate (s->statusPtr, z, &intinp, &y);
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

void bbad (Status *s, REAL8 *y, REAL8 x, void *p)
{
  *y = (REAL8)(++(*(INT4 *)p));
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
  const REAL4   sepsilon = 1e-6;
  const REAL8   depsilon = 1e-13; /* not as good as expected (1e-15) */
  static Status status;
  SIntegrateIn  sintinp;
  DIntegrateIn  dintinp;
  REAL4         sresult;
  REAL8         dresult;
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

  sintinp.function = f1;
  sintinp.xmin     = 0;
  sintinp.xmax     = 2;
  sintinp.type     = ClosedInterval;
  dintinp.function = ff1;
  dintinp.xmin     = 0;
  dintinp.xmax     = 2;
  dintinp.type     = ClosedInterval;

  count  = 0;
  expect = 8.153364119811650205L;
  SRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", sresult);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(sresult - expect) > sepsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }
  count = 0;
  DRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", dresult);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(dresult - expect) > depsilon*fabs(expect))
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
  sintinp.function = f2;
  sintinp.xmin     = 10;
  sintinp.xmax     = 1e30;
  sintinp.type     = InfiniteDomainPow;
  dintinp.function = ff2;
  dintinp.xmin     = 10;
  dintinp.xmax     = 1e300;
  dintinp.type     = InfiniteDomainPow;

  count  = 0;
  expect = 1.0L/200.0L;
  SRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", sresult);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(sresult - expect) > sepsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count = 0;
  DRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", dresult);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(dresult - expect) > depsilon*fabs(expect))
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
  sintinp.function = f3;
  sintinp.xmin     = 2;
  sintinp.xmax     = 1e30;
  sintinp.type     = InfiniteDomainExp;
  dintinp.function = ff3;
  dintinp.xmin     = 2;
  dintinp.xmax     = 1e300;
  dintinp.type     = InfiniteDomainExp;

  count  = 0;
  expect = 0.0570261239928920483L;
  SRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", sresult);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(sresult - expect) > sepsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count = 0;
  DRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", dresult);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(dresult - expect) > depsilon*fabs(expect))
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
  sintinp.function = f4;
  sintinp.xmin     = 0;
  sintinp.xmax     = 1;
  sintinp.type     = SingularLowerLimit;
  dintinp.function = ff4;
  dintinp.xmin     = 0;
  dintinp.xmax     = 1;
  dintinp.type     = SingularLowerLimit;

  count  = 0;
  expect = 2.0L;
  SRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", sresult);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(sresult - expect) > sepsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count  = 0;
  DRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", dresult);
  printf ("expect: %.15Lf\n", expect);
  if (fabs(dresult - expect) > depsilon*fabs(expect))
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
  sintinp.function = f5;
  sintinp.xmin     = 4;
  sintinp.xmax     = 5;
  sintinp.type     = SingularUpperLimit;
  dintinp.function = ff5;
  dintinp.xmin     = 4;
  dintinp.xmax     = 5;
  dintinp.type     = SingularUpperLimit;

  count  = 0;
  expect = 6.5L;
  SRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", sresult);
  printf ("expect: %.15Lf\n", expect);
  /* this doesn't work so well: multiply tolerance by factor of three */
  if (fabs(sresult - expect) > 3*sepsilon*fabs(expect))
  {
    fprintf (stderr, "Integration did not achieve desired accuracy!\n");
    return 1;
  }

  count  = 0;
  DRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(0), 1);
  printf ("number of function calls: %d\n", count);
  printf ("result: %.15f\n", dresult);
  printf ("expect: %.15Lf\n", expect);
  /* this doesn't work so well: multiply tolerance by factor of three */
  if (fabs(dresult - expect) > 3*depsilon*fabs(expect))
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
  sintinp.function = h;
  sintinp.xmin     = 0;
  sintinp.xmax     = 10;
  sintinp.type     = OpenInterval;
  dintinp.function = hh;
  dintinp.xmin     = 0;
  dintinp.xmax     = 10;
  dintinp.type     = OpenInterval;

  expect = 0.88274109326014810823L;
  SRombergIntegrate (&status, &sresult, &sintinp, NULL);
  TestStatus (&status, CODES(0), 1);
  printf ("result: %.15f\n", sresult);
  printf ("expect: %.15Lf\n", expect);
  /* integral isn't very accurate because we needed to use an open interval */
  printf ("error:  %.2f%%\n", 100*fabs(sresult - expect)/fabs(expect));
  /*
   * don't do 2d double-precision: it takes too long!
   *
   * DRombergIntegrate (&status, &dresult, &dintinp, NULL);
   * TestStatus (&status, CODES(0), 1);
   * printf ("result: %.15f\n", dresult);
   * printf ("expect: %.15Lf\n", expect);
   * printf ("error:  %.2f%%\n", 100*fabs(dresult - expect)/fabs(expect));
   *
   */


  LALCheckMemoryLeaks ();


  /*
   *
   * Check error conditions.
   *
   */


  printf ("\nChecking error conditions:\n");

  printf ("\nNull pointer:\r");
  SRombergIntegrate (&status, NULL, &sintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_ENULL), 1);
  DRombergIntegrate (&status, NULL, &dintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  printf ("\nNull pointer:\r");
  SRombergIntegrate (&status, &sresult, NULL, &count);
  TestStatus (&status, CODES(INTEGRATE_ENULL), 1);
  DRombergIntegrate (&status, &dresult, NULL, &count);
  TestStatus (&status, CODES(INTEGRATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  printf ("\nNull pointer:\r");
  sintinp.function = NULL;
  sintinp.xmin     = 0;
  sintinp.xmax     = 2;
  sintinp.type     = ClosedInterval;
  dintinp.function = NULL;
  dintinp.xmin     = 0;
  dintinp.xmax     = 2;
  dintinp.type     = ClosedInterval;
  SRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_ENULL), 1);
  DRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_ENULL), 1);
  printf ("Null pointer check passed.\n");

  printf ("\nInvalid domain:\r");
  sintinp.function = f1;
  sintinp.xmin     = 0;
  sintinp.xmax     = 0;
  sintinp.type     = ClosedInterval;
  dintinp.function = ff1;
  dintinp.xmin     = 0;
  dintinp.xmax     = 0;
  dintinp.type     = ClosedInterval;
  SRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_EIDOM), 1);
  DRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_EIDOM), 1);
  printf ("Invalid domain check passed.\n");

  printf ("\nUnknown integral type:\r");
  sintinp.function = f1;
  sintinp.xmin     = 0;
  sintinp.xmax     = 2;
  sintinp.type     = 999;
  dintinp.function = ff1;
  dintinp.xmin     = 0;
  dintinp.xmax     = 2;
  dintinp.type     = 999;
  SRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_ETYPE), 1);
  DRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_ETYPE), 1);
  printf ("Unknown integral type check passed.\n");

  printf ("\nMaximum iterations exceeded:\r");
  sintinp.function = bad;  /* bad is a quick random number generator */
  sintinp.xmin     = 0;
  sintinp.xmax     = 2;
  sintinp.type     = ClosedInterval;
  dintinp.function = bbad;  /* bbad is a quick random number generator */
  dintinp.xmin     = 0;
  dintinp.xmax     = 2;
  dintinp.type     = ClosedInterval;
  count            = 13;   /* count is now used as a random number seed */
  SRombergIntegrate (&status, &sresult, &sintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_EMXIT), 1);
  count = 1;
  DRombergIntegrate (&status, &dresult, &dintinp, &count);
  TestStatus (&status, CODES(INTEGRATE_EMXIT), 1);
  printf ("Maximum iterations exceeded check passed.\n");

  printf ("\nRecursive error:\r");
  sintinp.function = f1;
  sintinp.xmin     = 0;
  sintinp.xmax     = 2;
  sintinp.type     = ClosedInterval;
  dintinp.function = ff1;
  dintinp.xmin     = 0;
  dintinp.xmax     = 2;
  dintinp.type     = ClosedInterval;
  SRombergIntegrate (&status, &sresult, &sintinp, NULL);
  TestStatus (&status, CODES(-1), 1);
  ClearStatus (&status);
  DRombergIntegrate (&status, &dresult, &dintinp, NULL);
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

