/*----------------------------------------------------------------------- 
 * 
 * File Name: FindRootTest.c
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
#include "FindRoot.h"

#define CODES_(x) #x
#define CODES(x) CODES_(x)

NRCSID (MAIN, "$Id$");

extern char *optarg;
extern int   optind;

int debuglevel = 0;
int verbose    = 0;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (Status *status, const char *expectedCodes, int exitCode);

static void
ClearStatus (Status *status);

/*
 *
 * Function: y = F(x; p) = p + x*x
 *
 */
static void
F (Status *s, REAL4 *y, REAL4 x, void *p)
{
  REAL4 y0;
  INITSTATUS (s, "F", MAIN);
  ASSERT (y, s, 1, "Null pointer");
  ASSERT (p, s, 1, "Null pointer");
  y0 = *(REAL4 *)p;
  *y = y0 + x*x;
  RETURN (s);
}

static void
FF (Status *s, REAL8 *y, REAL8 x, void *p)
{
  REAL8 y0;
  INITSTATUS (s, "FF", MAIN);
  ASSERT (y, s, 1, "Null pointer");
  ASSERT (p, s, 1, "Null pointer");
  y0 = *(REAL8 *)p;
  *y = y0 + x*x;
  RETURN (s);
}


int
main (int argc, char *argv[])
{
  static Status  status;
  SFindRootIn    sinput;
  DFindRootIn    dinput;
  REAL4          y0;
  REAL4          sroot;
  REAL8          yy0;
  REAL8          droot;


  /*
   *
   * Parse the command line options
   *
   */


  ParseOptions (argc, argv);


  /*
   *
   * Set up input structure and function parameter y0.
   *
   */


  y0              = -1;
  sinput.function = F;
  sinput.xmin     = 1e-3;
  sinput.xmax     = 2e-3;
  sinput.xacc     = 1e-6;
  yy0             = -1;
  dinput.function = FF;
  dinput.xmin     = 1e-3;
  dinput.xmax     = 2e-3;
  dinput.xacc     = 1e-15;


  /*
   *
   * Check to see if bracketing and root finding work.
   *
   */


  if (verbose)
  {
    printf ("\n===== Check Root Finding =====\n\n");
  }

  if (verbose)
  {
    printf ("Initial domain: [%e,%e]\n", dinput.xmin, dinput.xmax);
  }

  SBracketRoot (&status, &sinput, &y0);
  TestStatus (&status, CODES(0), 1);

  if (verbose)
  {
    printf ("Bracket domain: [%e,%e]\n", sinput.xmin, sinput.xmax);
  }

  if (sinput.xmin > 1 || sinput.xmax < 1)
  {
    fprintf (stderr, "Root not bracketed correctly\n");
    return 1;
  }

  DBracketRoot (&status, &dinput, &yy0);
  TestStatus (&status, CODES(0), 1);

  if (verbose)
  {
    printf ("Bracket domain: [%e,%e]\n", dinput.xmin, dinput.xmax);
  }

  if (dinput.xmin > 1 || dinput.xmax < 1)
  {
    fprintf (stderr, "Root not bracketed correctly\n");
    return 1;
  }


  SBisectionFindRoot (&status, &sroot, &sinput, &y0);
  TestStatus (&status, CODES(0), 1);

  if (verbose)
  {
    printf ("Root = %e (acc = %e)\n", sroot, sinput.xacc);
  }

  if (fabs(sroot - 1) > sinput.xacc)
  {
    fprintf (stderr, "Root not found to correct accuracy\n");
    return 1;
  }


  DBisectionFindRoot (&status, &droot, &dinput, &yy0);
  TestStatus (&status, CODES(0), 1);

  if (verbose)
  {
    printf ("Root = %.15e (acc = %e)\n", droot, dinput.xacc);
  }

  if (fabs(droot - 1) > dinput.xacc)
  {
    fprintf (stderr, "Root not found to correct accuracy\n");
    return 1;
  }


  /*
   *
   * Check to make sure that correct error codes are generated.
   *
   */


  if (verbose || debuglevel)
  {
    printf ("\n===== Check Errors =====\n");
  }

  /* recursive error from an error occurring in the function */

  if (verbose)
  {
    printf ("\n----- Recursive Error: Code -1 (2 times)\n");
  }

  SBracketRoot (&status, &sinput, NULL);
  TestStatus  (&status, CODES(-1), 1);
  ClearStatus (&status);
  DBracketRoot (&status, &dinput, NULL);
  TestStatus  (&status, CODES(-1), 1);
  ClearStatus (&status);

  /* one of the arguments is a null pointer */

  if (verbose)
  {
    printf ("\n----- Null Pointer Error: Code 1 (10 times)\n");
  }

  SBracketRoot (&status, NULL, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);
  DBracketRoot (&status, NULL, &yy0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  SBisectionFindRoot (&status, NULL, &sinput, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);
  DBisectionFindRoot (&status, NULL, &dinput, &yy0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  SBisectionFindRoot (&status, &sroot, NULL, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);
  DBisectionFindRoot (&status, &droot, NULL, &yy0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  sinput.function = NULL;
  dinput.function = NULL;

  SBracketRoot (&status, &sinput, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);
  DBracketRoot (&status, &dinput, &yy0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  SBisectionFindRoot (&status, &sroot, &sinput, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);
  DBisectionFindRoot (&status, &droot, &dinput, &yy0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  /* invalid initial domain error for BracketRoot() */

  if (verbose)
  {
    printf ("\n----- Invalid Initial Domain: Code 2 (2 times)\n");
  }

  sinput.function = F;
  sinput.xmin     = 5;
  sinput.xmax     = 5;
  dinput.function = FF;
  dinput.xmin     = 5;
  dinput.xmax     = 5;

  SBracketRoot (&status, &sinput, &y0);
  TestStatus (&status, CODES(FINDROOT_EIDOM), 1);
  DBracketRoot (&status, &dinput, &yy0);
  TestStatus (&status, CODES(FINDROOT_EIDOM), 1);

  /* maximum iterations exceeded error */

  if (verbose)
  {
    printf ("\n----- Maximum Iteration Exceeded: Code 4 (4 times)\n");
  }
  
  y0              = 1; /* there is no root when y0 > 0 */
  sinput.xmin     = -1e-18;
  sinput.xmax     = 1e-18;
  yy0             = 1; /* there is no root when y0 > 0 */
  dinput.xmin     = -1e-18;
  dinput.xmax     = 1e-18;

  SBracketRoot (&status, &sinput, &y0);
  TestStatus (&status, CODES(FINDROOT_EMXIT), 1);
  DBracketRoot (&status, &dinput, &yy0);
  TestStatus (&status, CODES(FINDROOT_EMXIT), 1);

  y0              = -1;
  sinput.xmin     = 0;
  sinput.xmax     = 1e19;
  sinput.xacc     = 2e-38;
  yy0             = -1;
  dinput.xmin     = 0;
  dinput.xmax     = 1e19;
  dinput.xacc     = 2e-38;

  SBisectionFindRoot (&status, &sroot, &sinput, &y0);
  TestStatus (&status, CODES(FINDROOT_EMXIT), 1);
  DBisectionFindRoot (&status, &droot, &dinput, &yy0);
  TestStatus (&status, CODES(FINDROOT_EMXIT), 1);

  /* root not bracketed error in BisectionFindRoot() */

  if (verbose)
  {
    printf ("\n----- Root Not Bracketed: Code 8 (2 times)\n");
  }

  sinput.xmin     = -5;
  sinput.xmax     = -3;
  sinput.xacc     = 1e-6;
  dinput.xmin     = -5;
  dinput.xmax     = -3;
  dinput.xacc     = 1e-6;

  SBisectionFindRoot (&status, &sroot, &sinput, &y0);
  TestStatus (&status, CODES(FINDROOT_EBRKT), 1);
  DBisectionFindRoot (&status, &droot, &dinput, &yy0);
  TestStatus (&status, CODES(FINDROOT_EBRKT), 1);

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

