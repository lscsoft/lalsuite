/*----------------------------------------------------------------------- 
 * 
 * File Name: RandomTest.c
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

#ifndef _LALCONFIG_H
#include "LALConfig.h"
#ifndef _LALCONFIG_H
#define _LALCONFIG_H
#endif
#endif

#ifdef HAVE_UNISTD_H
#ifndef _UNISTD_H
#include <unistd.h>
#ifndef _UNISTD_H
#define _UNISTD_H
#endif
#endif
#endif

#ifdef HAVE_GETOPT_H
#ifndef _GETOPT_H
#include <getopt.h>
#ifndef _GETOPT_H
#define _GETOPT_H
#endif
#endif
#endif

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _AVFACTORIES_H
#include "AVFactories.h"
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
#endif
#endif

#ifndef _RANDOM_H
#include "Random.h"
#ifndef _RANDOM_H
#define _RANDOM_H
#endif
#endif

#define _CODES(x) #x
#define CODES(x) _CODES(x)

NRCSID (MAIN, "$Id$");

extern char *optarg;
extern int   optind;

int output     = 0;
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

int
main (int argc, char *argv[])
{
  const  INT4          numPoints = 999;
  static Status        status;
  static REAL4Vector  *vector;
  static RandomParams *randpar;
  INT4                 i;

  
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

  CreateVector (&status, &vector, numPoints);
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

  CreateRandomParams (&status, &randpar, 0);
  TestStatus (&status, CODES(0), 1);


  /*
   *
   * Fill vector with uniform random numbers.
   *
   */


  for (i = 0; i < vector->length; ++i)
  {
    UniformDeviate (&status, vector->data + i, randpar);
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


  NormalDeviates (&status, vector, randpar);
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


  if (verbose || debuglevel)
  {
    printf ("\n===== Check Errors =====\n");
  }

  /* non-null pointer error */

  if (verbose)
  {
    printf ("\n----- Non-Null Pointer Error: Code 2\n");
  }

  CreateRandomParams (&status, &randpar, 0);
  TestStatus (&status, CODES(RANDOM_ENNUL), 1);

  /* null pointer error */

  if (verbose)
  {
    printf ("\n----- Null Pointer Error: Code 1 (8 times)\n");
  }

  CreateRandomParams (&status, NULL, 0);
  TestStatus (&status, CODES(RANDOM_ENULL), 1);

  DestroyRandomParams (&status, NULL);
  TestStatus (&status, CODES(RANDOM_ENULL), 1);

  {
    RandomParams *tmp = NULL;
    DestroyRandomParams (&status, &tmp);
  }
  TestStatus (&status, CODES(RANDOM_ENULL), 1);

  UniformDeviate (&status, NULL, randpar);
  TestStatus (&status, CODES(RANDOM_ENULL), 1);
  
  UniformDeviate (&status, vector->data, NULL);
  TestStatus (&status, CODES(RANDOM_ENULL), 1);

  NormalDeviates (&status, NULL, randpar);
  TestStatus (&status, CODES(RANDOM_ENULL), 1);

  {
    REAL4Vector tmp;
    tmp.length = 10;
    tmp.data   = NULL;
    NormalDeviates (&status, &tmp, randpar);
  }
  TestStatus (&status, CODES(RANDOM_ENULL), 1);
  
  NormalDeviates (&status, vector, NULL);
  TestStatus (&status, CODES(RANDOM_ENULL), 1);

  /* vector length error */

  if (verbose)
  {
    printf ("\n----- Invalid Length Error: Code 4\n");
  }

  {
    REAL4Vector tmp;
    tmp.length = 0;
    tmp.data   = (REAL4 *)1;
    NormalDeviates (&status, &tmp, randpar);
  }
  TestStatus (&status, CODES(RANDOM_ESIZE), 1);
  

  /*
   *
   * Free memory.
   *
   */


  if (verbose || debuglevel)
  {
    printf ("\n===== Clean up and Exit =====\n");
  }

  DestroyRandomParams (&status, &randpar);
  TestStatus (&status, CODES(0), 1);

  DestroyVector (&status, &vector);
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

