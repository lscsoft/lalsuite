/*----------------------------------------------------------------------- 
 * 
 * File Name: FrameDataTest.c
 *
 * Author: Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "LALConfig.h"

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include "LALStdlib.h"
#include "AVFactories.h"
#include "FrameData.h"

#define CODES_(x) #x
#define CODES(x) CODES_(x)

NRCSID (MAIN, "$Id$");

extern char *optarg;
extern int   optind;

int   debuglevel = 0;
int   verbose    = 0;
int   output     = 0;
char *framePath  = NULL;

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
  const INT4               numPoints = 262144;
  static Status            status;
  FrameData               *frameData = NULL;
  INT2TimeSeries           data;
  COMPLEX8FrequencySeries  response;
  INT4                     seg;

  ParseOptions (argc, argv);

  if (!framePath)
  {
    framePath = getenv ("LAL_FRAME_PATH");
    if (!framePath)
    {
      fprintf (stderr, "error: environment LAL_FRAME_PATH undefined\n");
      return 77;
    }
  }

  data.data = NULL;
  I2CreateVector (&status, &data.data, numPoints);
  TestStatus (&status, "0", 1);

  response.data = NULL;
  CCreateVector (&status, &response.data, numPoints/2 + 1);
  TestStatus (&status, "0", 1);

  InitializeFrameData (&status, &frameData, framePath);
  TestStatus (&status, CODES(0 FRAMEDATA_EREAD), 1);

  for (seg = 0; seg < 99; ++seg)
  {
    GetFrameData (&status, &data, frameData);
    TestStatus (&status, "0", 1);

    fprintf (stderr, "Segment %2d... ", seg);

    if (frameData->endOfData)
    {
      fprintf (stderr, "end of data\n");
      break;
    }

    if (frameData->newLock)
    {
      fprintf (stderr, "starting new locked section, ");
    }
    else
    {
      fprintf (stderr, "continuing locked section, ");
    }

    if (frameData->newCalibration)
    {
      fprintf (stderr, "new calibration data, ");

      GetFrameDataResponse (&status, &response, frameData);
      TestStatus (&status, "0", 1);

      if (output)
      {
        FILE *fp;
        CHAR  fname[64];
        INT4  i;

        sprintf (fname, "Response.%03d", seg);
        fp = fopen (fname, "w");

        for (i = 0; i < response.data->length; ++i)
        {
          fprintf (fp, "%d\t%e\t%e\n", i, response.data->data[i].re,
                   response.data->data[i].im);
        }

        fclose (fp);
      }
    }

    fprintf (stderr, "gps time %d.%09d\n", data.epoch.gpsSeconds,
             data.epoch.gpsNanoSeconds);

    if (output)
    {
      FILE *fp;
      CHAR  fname[64];
      INT4  i;

      sprintf (fname, "Segment.%03d", seg);
      fp = fopen (fname, "w");

      for (i = 0; i < data.data->length; ++i)
      {
        fprintf (fp, "%d\t%d\n", i, data.data->data[i]);
      }

      fclose (fp);
    }

  }


  /* again... just to be annoying if end of data has happened */
  ++seg;

  GetFrameData (&status, &data, frameData);
  TestStatus (&status, "0", 1);

  fprintf (stderr, "Segment %2d... ", seg);
  if (frameData->endOfData)
  {
    fprintf (stderr, "end of data\n");
  }

  /* clean up */

  FinalizeFrameData (&status, &frameData);
  TestStatus (&status, "0", 1);

  CDestroyVector    (&status, &response.data);
  TestStatus (&status, "0", 1);

  I2DestroyVector   (&status, &data.data);
  TestStatus (&status, "0", 1);

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
  fprintf (stderr, "  -o         output framedata to files\n");
  fprintf (stderr, "  -f dir     set frame data path to dir\n");
  fprintf (stderr, "             "
                   "(otherwise use path in environment LAL_FRAME_PATH)\n");
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

    c = getopt (argc, argv, "hqvd:""of:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'f': /* sets frame path */
        framePath = optarg;
        break;

      case 'o': /* sets flag to write output files */
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

