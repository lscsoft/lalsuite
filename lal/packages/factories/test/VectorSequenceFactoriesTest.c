/**
\file
\ingroup SeqFactories_h

\brief A program to test create/destroy vector sequence routines.

\heading{Usage}
\code
VectorSequenceFactoriesTest [options]
Options:
  -h         print help
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
\endcode

\heading{Exit codes}
<table><tr><th>Code</th><th>Explanation</th></tr>
<tr><td>0</td><td>Success, normal exit.</td></tr>
<tr><td>1</td><td>Subroutine failed.</td></tr>
</table>

\heading{Algorithm}

\heading{Uses}
\code
lalDebugLevel
<datatype>CreateVectorSequence()
<datatype>DestroyVectorSequence()
\endcode

\heading{Notes}

*/

/** \cond DONT_DOXYGEN */
#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/SeqFactories.h>

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

#define TYPECODE Z
#define TYPE COMPLEX16
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE C
#define TYPE COMPLEX8
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE D
#define TYPE REAL8
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U2
#define TYPE UINT2
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U4
#define TYPE UINT4
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U8
#define TYPE UINT8
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE CHAR
#define TYPE CHAR
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPE REAL4
#include "VectorSequenceFactoriesTest_source.c"
#undef TYPE

int main( int argc, char *argv[] )
{
  lalDebugLevel = 0;

  ParseOptions( argc, argv );

  VectorSequenceFactoriesTest();
  ZVectorSequenceFactoriesTest();
  CVectorSequenceFactoriesTest();
  DVectorSequenceFactoriesTest();
  SVectorSequenceFactoriesTest();
  I2VectorSequenceFactoriesTest();
  I4VectorSequenceFactoriesTest();
  I8VectorSequenceFactoriesTest();
  U2VectorSequenceFactoriesTest();
  U4VectorSequenceFactoriesTest();
  U8VectorSequenceFactoriesTest();
  CHARVectorSequenceFactoriesTest();

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
