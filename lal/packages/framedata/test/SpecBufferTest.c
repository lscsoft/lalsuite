/*----------------------------------------------------------------------- 
 * 
 * File Name: SpecBufferTest.c
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
#include <math.h>
#include <lal/LALConfig.h>

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/FrameData.h>
#include <lal/SpecBuffer.h>

#define CODES_(x) #x
#define CODES(x) CODES_(x)

NRCSID (MAIN, "$Id$");

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

static void
ClearStatus (LALStatus *status);

int
main (int argc, char *argv[])
{
  const INT4 numPoints = 65536;
  const INT4 numSpec   = 8;
  const INT4 numSegs   = 10;

  static LALStatus            status;
  FrameData               *frameData = NULL;
  INT2TimeSeries           data;
  COMPLEX8FrequencySeries  response;
  REAL4FrequencySeries     spectrum;
  SpectrumBuffer          *specBuff  = NULL;
  SpectrumBufferPar        specParm;
  INT4                     newLock   = 0;
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
  TestStatus (&status, CODES(0 FRAMEDATA_EREAD), 1);

  spectrum.data = NULL;
  LALCreateVector (&status, &spectrum.data, numPoints/2 + 1);
  TestStatus (&status, "0", 1);

  specParm.numSpec    = numSpec;
  specParm.numPoints  = numPoints;
  specParm.windowType = Welch;
  specParm.plan       = NULL;

  LALEstimateFwdRealFFTPlan (&status, &specParm.plan, numPoints);
  TestStatus (&status, "0", 1);

  LALCreateSpectrumBuffer (&status, &specBuff, &specParm);
  TestStatus (&status, "0", 1);

  for (seg = 0; seg < numSegs; ++seg)
  {
    INT4 newCal = 0;

    fprintf (stderr, "Segment %2d", seg);

    /* get next segment of locked data */
    do
    {
      while (newLock)
      {
        INT2TimeSeries dummy;
        INT2Vector     dumvec;

        dumvec.length   = 3*60/data.deltaT; /* 3 minutes */
        dumvec.data     = NULL;             /* seek mode */
        dummy.data      = &dumvec;

        LALGetFrameData (&status, &dummy, frameData);
        TestStatus (&status, "0", 1);
        if (frameData->endOfData)
        {
          fprintf (stderr, "... end of data\n");
          goto exit;
        }
        fprintf (stderr, "... skipping 3 minutes into new lock");
        newLock = frameData->newLock;
        if (frameData->newCalibration)
        {
          newCal = 1;
        }
      }

      LALGetFrameData (&status, &data, frameData);
      TestStatus (&status, "0", 1);
      if (frameData->endOfData)
      {
        fprintf (stderr, "... end of data\n");
        goto exit;
      }
      newLock = frameData->newLock;
      if (frameData->newCalibration)
      {
        newCal = 1;
      }
    }
    while (newLock);

    if (newCal)
    {
      LALGetFrameDataResponse (&status, &response, frameData);
      TestStatus (&status, "0", 1);

      /* print calibration info */
      if (output)
      {
        FILE *fp;
        CHAR  fname[64];
        INT4  i;

        sprintf (fname, "Response.%03d", seg);
        fp = fopen (fname, "w");

        for (i = 0; i < response.data->length; ++i)
        {
          REAL4 re  = response.data->data[i].re;
          REAL4 im  = response.data->data[i].im;
          REAL4 mod = sqrt (re*re + im*im);
          REAL4 arg = atan2 (im, re);
          fprintf (fp, "%u\t%e\t%e\n", i, mod, arg);
        }

        fclose (fp);
      }

      fprintf (stderr, "... new calibration data");

    }


    /* print data vector */
    if (output)
    {
      FILE *fp;
      CHAR  fname[64];
      INT4  i;

      sprintf (fname, "Segment.%03d", seg);
      fp = fopen (fname, "w");

      for (i = 0; i < data.data->length; ++i)
      {
        fprintf (fp, "%u\t%d\n", i, data.data->data[i]);
      }

      fclose (fp);
    }

    LALAddSpectrum (&status, specBuff, &data);
    TestStatus (&status, "0", 1);

    fprintf (stderr, "\n");
  }

exit:

  LALAverageSpectrum (&status, &spectrum, specBuff);
  TestStatus (&status, CODES(0 SPECBUFFER_ENONE), 1);

  if (output)
  {
    FILE *fp;
    INT4  i;

    fp = fopen ("Spectrum.avg", "w");

    for (i = 0; i < spectrum.data->length; ++i)
      fprintf (fp, "%u\t%e\n", i, spectrum.data->data[i]);

    fclose (fp);
  }

  LALDestroyVector (&status, &spectrum.data);
  TestStatus (&status, "0", 1);

  LALDestroySpectrumBuffer (&status, &specBuff);
  TestStatus (&status, "0", 1);

  LALDestroyRealFFTPlan (&status, &specParm.plan);
  TestStatus (&status, "0", 1);

  LALFinalizeFrameData (&status, &frameData);
  TestStatus (&status, "0", 1);

  LALI2DestroyVector (&status, &data.data);
  TestStatus (&status, "0", 1);

  LALCDestroyVector (&status, &response.data);
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
 *
 * ClearStatus ()
 *
 * Recursively applies DETATCHSTATUSPTR() to status structure to destroy
 * linked list of statuses.
 *
 */
void
ClearStatus (LALStatus *status)
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

