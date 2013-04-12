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
\brief Tests the routines in \ref FrameData.h.

\heading{Usage}
\code
FrameDataTest [options]
Options:
  -h         print this message
  -q         quiet: run silently
  -v         verbose: print extra information
  -d level   set lalDebugLevel to level
  -o         output framedata to files
  -f dir     set frame data path to dir
\endcode

Unless the <tt>-f</tt> option is used, the environment variable
\c LAL_FRAME_PATH must be set to the directory containing the frame files.

\heading{Description}
\heading{Exit codes}
<table>
<tr><th> Code</th><th>Explanation</th></tr>
<tr><td>\c 0</td><td>Success, normal exit.</td></tr>
<tr><td>\c 1</td><td>Subroutine failed.</td></tr>
<tr><td>\c 77</td><td>Ignored failure: \c LAL_FRAME_PATH not set.</td></tr>
</table>

\heading{Uses}
\heading{Notes}

*/
#include <config.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/FrameData.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

extern char *optarg;
extern int   optind;

int   lalDebugLevel = 0;
int   verbose    = 0;
int   output     = 0;
char *framePath  = NULL;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

int
main (int argc, char *argv[])
{
  const INT4               numPoints = 262144;
  static LALStatus            status;
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
  LALI2CreateVector (&status, &data.data, numPoints);
  TestStatus (&status, "0", 1);

  response.data = NULL;
  LALCCreateVector (&status, &response.data, numPoints/2 + 1);
  TestStatus (&status, "0", 1);

  LALInitializeFrameData (&status, &frameData, framePath);
  TestStatus (&status, CODES(0 FRAMEDATAH_EREAD), 1);

  for (seg = 0; seg < 99; ++seg)
  {
    LALGetFrameData (&status, &data, frameData);
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

      LALGetFrameDataResponse (&status, &response, frameData);
      TestStatus (&status, "0", 1);

      if (output)
      {
        FILE  *fp;
        CHAR   fname[64];
        UINT4  i;

        sprintf (fname, "Response.%03d", seg);
        fp = fopen (fname, "w");

        for (i = 0; i < response.data->length; ++i)
        {
          fprintf (fp, "%d\t%e\t%e\n", i, crealf(response.data->data[i]),
                   response.data->data[i].im);
        }

        fclose (fp);
      }
    }

    fprintf (stderr, "gps time %d.%09d\n", data.epoch.gpsSeconds,
             data.epoch.gpsNanoSeconds);

    if (output)
    {
      FILE  *fp;
      CHAR   fname[64];
      UINT4  i;

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

  LALGetFrameData (&status, &data, frameData);
  TestStatus (&status, "0", 1);

  fprintf (stderr, "Segment %2d... ", seg);
  if (frameData->endOfData)
  {
    fprintf (stderr, "end of data\n");
  }

  /* clean up */

  LALFinalizeFrameData (&status, &frameData);
  TestStatus (&status, "0", 1);

  LALCDestroyVector    (&status, &response.data);
  TestStatus (&status, "0", 1);

  LALI2DestroyVector   (&status, &data.data);
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
  FILE *fp;

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

