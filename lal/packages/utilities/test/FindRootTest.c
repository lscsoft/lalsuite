/*----------------------------------------------------------------------- 
 * 
 * File Name: FindRootTest.c
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#ifndef _STDIO_H
#include <stdio.h>
#ifndef _STDIO_H
#define _STDIO_H
#endif
#endif

#ifndef _STRING_H
#include <string.h>
#ifndef _STRING_H
#define _STRING_H
#endif
#endif

#ifndef _STDLIB_H
#include <stdlib.h>
#ifndef _STDLIB_H
#define _STDLIB_H
#endif
#endif

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _FINDROOT_H
#include "FindRoot.h"
#ifndef _FINDROOT_H
#define _FINDROOT_H
#endif
#endif

#define _CODES(x) #x
#define CODES(x) _CODES(x)

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

/*
 *
 * Function: y = F(x; p) = p + x*x
 *
 */
static void
F (Status *s, REAL4 *y, REAL4 x, void *p)
{
  REAL4 y0;
  INITSTATUS (s, "Test function");
  ASSERT (y, s, 1, "Null pointer");
  ASSERT (p, s, 1, "Null pointer");
  y0 = *(REAL4 *)p;
  *y = y0 + x*x;
  RETURN (s);
}


int
main (int argc, char *argv[])
{
  static Status  status;
  FindRootIn     input;
  REAL4          y0;
  REAL4          root;


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


  y0             = -1;
  input.function = F;
  input.xmin     = 1e-3;
  input.xmax     = 2e-3;
  input.xacc     = 1e-6;


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
    printf ("Initial domain: [%e,%e]\n", input.xmin, input.xmax);
  }

  BracketRoot (&status, &input, &y0);
  TestStatus (&status, CODES(0), 1);

  if (verbose)
  {
    printf ("Bracket domain: [%e,%e]\n", input.xmin, input.xmax);
  }

  if (input.xmin > 1 || input.xmax < 1)
  {
    fprintf (stderr, "Root not bracketed correctly\n");
    return 1;
  }

  BisectionFindRoot (&status, &root, &input, &y0);
  TestStatus (&status, CODES(0), 1);

  if (verbose)
  {
    printf ("Root = %e (acc = %e)\n", root, input.xacc);
  }

  if (fabs(root - 1) > input.xacc)
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
    printf ("\n----- Recursive Error: Code -1\n");
  }

  BracketRoot (&status, &input, NULL);
  TestStatus (&status, CODES(-1), 1);

  /* one of the arguments is a null pointer */

  if (verbose)
  {
    printf ("\n----- Null Pointer Error: Code 1 (5 times)\n");
  }

  BracketRoot (&status, NULL, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  BisectionFindRoot (&status, NULL, &input, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  BisectionFindRoot (&status, &root, NULL, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  input.function = NULL;

  BracketRoot (&status, &input, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  BisectionFindRoot (&status, &root, &input, &y0);
  TestStatus (&status, CODES(FINDROOT_ENULL), 1);

  /* invalid initial domain error for BracketRoot() */

  if (verbose)
  {
    printf ("\n----- Invalid Initial Domain: Code 2\n");
  }

  input.function = F;
  input.xmin     = 5;
  input.xmax     = 5;

  BracketRoot (&status, &input, &y0);
  TestStatus (&status, CODES(FINDROOT_EIDOM), 1);

  /* maximum iterations exceeded error */

  if (verbose)
  {
    printf ("\n----- Maximum Iteration Exceeded: Code 4 (2 times)\n");
  }
  
  y0             = 1; /* there is no root when y0 > 0 */
  input.xmin     = -5;
  input.xmax     = 5;

  BracketRoot (&status, &input, &y0);
  TestStatus (&status, CODES(FINDROOT_EMXIT), 1);

  y0             = -1;
  input.xmin     = 0;
  input.xmax     = 1e19;
  input.xacc     = 2e-38;

  BisectionFindRoot (&status, &root, &input, &y0);
  TestStatus (&status, CODES(FINDROOT_EMXIT), 1);

  /* root not bracketed error in BisectionFindRoot() */

  if (verbose)
  {
    printf ("\n----- Root Not Bracketed: Code 8\n");
  }

  input.xmin     = -5;
  input.xmax     = -3;
  input.xacc     = 1e-6;

  BisectionFindRoot (&status, &root, &input, &y0);
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

