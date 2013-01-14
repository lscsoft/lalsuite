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
   \file
   \ingroup FindRoot_h

   \brief Tests the routines in \ref FindRoot.h.

\heading{Usage}
\code
FindRootTest [options]
Options:
  -h         print this message
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
\endcode

\heading{Exit codes}

<table>
<tr><th>Code</th><th>Explanation</th></tr>
<tr><td>0</td><td>Success, normal exit.</td></tr>
<tr><td>1</td><td>Subroutine failed.</td></tr>
</table>

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
#include <lal/FindRoot.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

extern char *optarg;
extern int   optind;

extern int lalDebugLevel;
int verbose    = 0;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

#if defined(NDEBUG) || defined(LAL_NDEBUG)
/* debugging is turned off */
#else
static void
ClearStatus (LALStatus *status);
#endif

/*
 *
 * Function: y = F(x; p) = p + x*x
 *
 */
static void
F (LALStatus *s, REAL4 *y, REAL4 x, void *p)
{
  REAL4 y_0;
  INITSTATUS(s);
  ASSERT (y, s, 1, "Null pointer");
  ASSERT (p, s, 1, "Null pointer");
  y_0 = *(REAL4 *)p;
  *y = y_0 + x*x;
  RETURN (s);
}

static void
FF (LALStatus *s, REAL8 *y, REAL8 x, void *p)
{
  REAL8 y_0;
  INITSTATUS(s);
  ASSERT (y, s, 1, "Null pointer");
  ASSERT (p, s, 1, "Null pointer");
  y_0 = *(REAL8 *)p;
  *y = y_0 + x*x;
  RETURN (s);
}


int
main (int argc, char *argv[])
{
  static LALStatus  status;
  SFindRootIn    sinput;
  DFindRootIn    dinput;
  REAL4          y_0;
  REAL4          sroot;
  REAL8          yy0;
  REAL8          droot;

  lalDebugLevel = 0;

  /*
   *
   * Parse the command line options
   *
   */


  ParseOptions (argc, argv);


  /*
   *
   * Set up input structure and function parameter y_0.
   *
   */


  y_0             = -1;
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

  LALSBracketRoot (&status, &sinput, &y_0);
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

  LALDBracketRoot (&status, &dinput, &yy0);
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


  LALSBisectionFindRoot (&status, &sroot, &sinput, &y_0);
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


  LALDBisectionFindRoot (&status, &droot, &dinput, &yy0);
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

#ifndef LAL_NDEBUG

  if ( ! lalNoDebug )
  {

    if (verbose || lalDebugLevel)
    {
      printf ("\n===== Check Errors =====\n");
    }

    /* recursive error from an error occurring in the function */

    if (verbose)
    {
      printf ("\n----- Recursive Error: Code -1 (2 times)\n");
    }

    LALSBracketRoot (&status, &sinput, NULL);
    TestStatus  (&status, CODES(-1), 1);
    ClearStatus (&status);
    LALDBracketRoot (&status, &dinput, NULL);
    TestStatus  (&status, CODES(-1), 1);
    ClearStatus (&status);

    /* one of the arguments is a null pointer */

    if (verbose)
    {
      printf ("\n----- Null Pointer Error: Code 1 (10 times)\n");
    }

    LALSBracketRoot (&status, NULL, &y_0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);
    LALDBracketRoot (&status, NULL, &yy0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);

    LALSBisectionFindRoot (&status, NULL, &sinput, &y_0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);
    LALDBisectionFindRoot (&status, NULL, &dinput, &yy0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);

    LALSBisectionFindRoot (&status, &sroot, NULL, &y_0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);
    LALDBisectionFindRoot (&status, &droot, NULL, &yy0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);

    sinput.function = NULL;
    dinput.function = NULL;

    LALSBracketRoot (&status, &sinput, &y_0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);
    LALDBracketRoot (&status, &dinput, &yy0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);

    LALSBisectionFindRoot (&status, &sroot, &sinput, &y_0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);
    LALDBisectionFindRoot (&status, &droot, &dinput, &yy0);
    TestStatus (&status, CODES(FINDROOTH_ENULL), 1);

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

    LALSBracketRoot (&status, &sinput, &y_0);
    TestStatus (&status, CODES(FINDROOTH_EIDOM), 1);
    LALDBracketRoot (&status, &dinput, &yy0);
    TestStatus (&status, CODES(FINDROOTH_EIDOM), 1);

    /* maximum iterations exceeded error */

    if (verbose)
    {
      printf ("\n----- Maximum Iteration Exceeded: Code 4 (4 times)\n");
    }

    y_0             = 1; /* there is no root when y_0 > 0 */
    sinput.xmin     = -1e-18;
    sinput.xmax     = 1e-18;
    yy0             = 1; /* there is no root when y_0 > 0 */
    dinput.xmin     = -1e-18;
    dinput.xmax     = 1e-18;

    LALSBracketRoot (&status, &sinput, &y_0);
    TestStatus (&status, CODES(FINDROOTH_EMXIT), 1);
    LALDBracketRoot (&status, &dinput, &yy0);
    TestStatus (&status, CODES(FINDROOTH_EMXIT), 1);

    y_0             = -1;
    sinput.xmin     = 0;
    sinput.xmax     = 1e19;
    sinput.xacc     = 2e-38;
    yy0             = -1;
    dinput.xmin     = 0;
    dinput.xmax     = 1e19;
    dinput.xacc     = 2e-38;

    LALSBisectionFindRoot (&status, &sroot, &sinput, &y_0);
    TestStatus (&status, CODES(FINDROOTH_EMXIT), 1);
    LALDBisectionFindRoot (&status, &droot, &dinput, &yy0);
    TestStatus (&status, CODES(FINDROOTH_EMXIT), 1);

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

    LALSBisectionFindRoot (&status, &sroot, &sinput, &y_0);
    TestStatus (&status, CODES(FINDROOTH_EBRKT), 1);
    LALDBisectionFindRoot (&status, &droot, &dinput, &yy0);
    TestStatus (&status, CODES(FINDROOTH_EBRKT), 1);

  }

#endif

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


#if defined(NDEBUG) || defined(LAL_NDEBUG)
/* debugging is turned off */
#else
/*
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
static void
ClearStatus (LALStatus *status)
{
  if (status->statusPtr)
  {
    ClearStatus      (status->statusPtr);
    DETATCHSTATUSPTR (status);
  }
}
#endif

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
