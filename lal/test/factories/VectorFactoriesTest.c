/**
 * \file
 * \ingroup AVFactories_h
 *
 * \brief A program to test create/destroy vector routines.
 *
 * ### Usage ###
 *
 * \code
 * VectorFactoriesTest [options]
 * Options:
 * -h         print help
 * -q         quiet: run silently
 * -v         verbose: print extra information
 * -d level   set lalDebugLevel to level
 * \endcode
 *
 * ### Description ###
 *
 *
 * ### Exit codes ###
 *
 * <table><tr><th>Code</th><th>Explanation</th></tr>
 * <tr><td>0</td><td>Success, normal exit.</td></tr>
 * <tr><td>1</td><td>Subroutine failed.</td></tr>
 * </table>
 *
 * ### Algorithm ###
 *
 *
 * ### Uses ###
 *
 * \code
 * lalDebugLevel
 * \<datatype\>CreateVector()
 * \<datatype\>ResizeVector()
 * \<datatype\>DestroyVector()
 * \endcode
 *
 * ### Notes ###
 *
 */
/** \cond DONT_DOXYGEN */
#include <config.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <lal/PrintVector.h>

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/AVFactories.h>
#include <lal/LALString.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

int verbose    = 0;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

#define TYPECODE Z
#define TYPE COMPLEX16
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE C
#define TYPE COMPLEX8
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE D
#define TYPE REAL8
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE S
#define TYPE REAL4
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I2
#define TYPE INT2
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I4
#define TYPE INT4
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE I8
#define TYPE INT8
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U2
#define TYPE UINT2
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U4
#define TYPE UINT4
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE U8
#define TYPE UINT8
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPECODE CHAR
#define TYPE CHAR
#include "VectorFactoriesTest_source.c"
#undef TYPECODE
#undef TYPE

#define TYPE REAL4
#include "VectorFactoriesTest_source.c"
#undef TYPE

int main( int argc, char *argv[] )
{

  ParseOptions( argc, argv );

  VectorFactoriesTest();
  ZVectorFactoriesTest();
  CVectorFactoriesTest();
  DVectorFactoriesTest();
  SVectorFactoriesTest();
  I2VectorFactoriesTest();
  I4VectorFactoriesTest();
  I8VectorFactoriesTest();
  U2VectorFactoriesTest();
  U4VectorFactoriesTest();
  U8VectorFactoriesTest();
  CHARVectorFactoriesTest();

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

    c = LALgetopt (argc, argv, "hqvd:");
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

/** \endcond DONT_DOXYGEN */
