/*
*  Copyright (C) 2007 Duncan Brown
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
 * \ingroup ResampleTimeSeries_h
 * \file
 * \author Brown, D. A.
 *
 * \brief Tests the routines in \ref ResampleTimeSeries.h
 *
 * ### Usage ###
 *
 * \code
 * Usage: ./ResampleTimeSeriesTest [options] defaults shown in brackets
 * Options:
 * -h         print this message
 * -d level   set lalDebugLevel to level
 * -v         verbose: print extra information
 * -n points  number of points in the raw time series (1048576)
 * -i freq    sample rate of input time series (16384)
 * -o freq    sample rate of output time series (4096)
 * -f freq    frequency of sine wave to inject as input (1000.0)
 * -r type    type of filter to use in resampling (ldas)
 * [ldas|butterworth]
 * \endcode
 *
 * ### Description ###
 *
 * Tests the resampling functions by injecting a sine wave of a given
 * frequency into a time series and downsampling it. The raw and output
 * data are returned as frame files for plotting in matlab.
 *
 * ### Sample Results ###
 *
 * Figures.\figref{f_resamp1}-\figref{f_resamp3} show the results of various
 * tests using this program.
 *
 * \image html  resamp_figs1.png "Fig. [f_resamp1]"
 * \image latex resamp_figs1.pdf ""
 *
 * Fig. [f_resamp1]: The left figure shows a 10 Hz sine wave generated at 16384 Hz resampled to
 * 4096 Hz. The right figure shows a 100Hz sine wave generated at 16384 Hz
 * resampled to 4096 Hz. Note that there is no attenuation, time delay or
 * phase shift of the output. FIXME the legend in the right figure is wrong.
 * It should say 100 Hz, not 10 Hz the output.
 *
 * \image html  resamp_figs2.png "Fig. [f_resamp2]"
 * \image latex resamp_figs2.pdf ""
 *
 * Fig. [f_resamp2]: A 100Hz sine wave generated at 16384 Hz resampled to 8192 Hz. The left
 * plot shows the start of the time series and the right plot the end. Note
 * the corruption of points due to the time domain filtering.
 *
 * \image html  resamp_figs3.png "Fig. [f_resamp3]"
 * \image latex resamp_figs3.pdf ""
 *
 * Fig. [f_resamp3]: The left figure shows a 1000 Hz sine wave generated at 16384 Hz resampled
 * to 4096 Hz. The right figure shows a 1000Hz sine wave generated at 16384 Hz
 * resampled to 2048 Hz. Note that there is no attenuation, time delay or
 * phase shift of the output at 4096 Hz, however there is attenuation and
 * phase shift of the output at 2048 Hz. This is due to the fact that the
 * signal is very close to the output Nyquist frequency. Care should be taken
 * to downsample to a suitable rate to avoid this type of attenuation.
 *
 */

/** \cond DONT_DOXYGEN */

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <lal/LALConfig.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Units.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/LALFrStream.h>
#include <lal/AVFactories.h>
#include <lal/LALStdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif

extern char *optarg;
extern int   optind;

int     verbose = 0;
UINT4   numPoints = 1048576;
UINT4   inRate    = 16384;
UINT4   outRate   = 4096;
REAL4   sineFreq  = 1000.0;
ResampleTSFilter filtType = LDASfirLP;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

static void
TestStatus (LALStatus *status, const char *expectedCodes, int exitCode);

int
main (int argc, char *argv[])
{
  static LALStatus      status;
  REAL4TimeSeries       chan;
  ResampleTSParams      resampPars;
  UINT4                 j;
  FrOutPar in_opar = { "F", "IN", ProcDataChannel, 1, 0, 0 };
  FrOutPar out_opar = { "F", "OUT", ProcDataChannel, 1, 0, 0 };

  ParseOptions (argc, argv);

  memset( &chan, 0, sizeof(REAL4TimeSeries) );
  LALSCreateVector( &status, &(chan.data), numPoints );
  TestStatus (&status, "0", 1);

  chan.sampleUnits = lalADCCountUnit;
  chan.deltaT = 1.0 / (REAL8) inRate;
  chan.epoch.gpsSeconds = 100;
  resampPars.deltaT = 1.0 / (REAL8) outRate;
  resampPars.filterType = filtType;

  for ( j = 0; j < chan.data->length; ++j )
  {
    chan.data->data[j] = sin( 2.0 * LAL_PI * sineFreq * j * chan.deltaT );
  }

  snprintf( chan.name, LALNameLength * sizeof(CHAR), "%d_%d_%d_%.2f_%d_in",
      inRate, outRate, numPoints, sineFreq, filtType );
  LALFrWriteREAL4TimeSeries( &status, &chan, &in_opar );
  TestStatus (&status, "0", 1);

  LALResampleREAL4TimeSeries( &status, &chan, &resampPars );
  TestStatus (&status, "0", 1);

  snprintf( chan.name, LALNameLength * sizeof(CHAR), "%d_%d_%d_%.2f_%d_out",
      inRate, outRate, numPoints, sineFreq, filtType );
  LALFrWriteREAL4TimeSeries( &status, &chan, &out_opar );
  TestStatus (&status, "0", 1);

  LALSDestroyVector( &status, &(chan.data) );
  TestStatus (&status, "0", 1);

  LALCheckMemoryLeaks();
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
  fprintf (stderr, "Usage: %s [options] defaults shown in brackets\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h         print this message\n");
  fprintf (stderr, "  -d level   set lalDebugLevel to level\n");
  fprintf (stderr, "  -v         verbose: print extra information\n");
  fprintf (stderr, "  -n points  number of points in the raw time series (1048576)\n");
  fprintf (stderr, "  -i freq    sample rate of input time series (16384)\n");
  fprintf (stderr, "  -o freq    sample rate of output time series (4096)\n");
  fprintf (stderr, "  -f freq    frequency of sine wave to inject as input (1000.0)\n");
  fprintf (stderr, "  -r type    type of filter to use in resampling (ldas)\n");
  fprintf (stderr, "             [ldas|butterworth]\n");
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

    c = getopt (argc, argv, "hvd:n:i:o:f:r:");
    if (c == -1)
    {
      break;
    }

    switch (c)
    {
      case 'h':
        Usage (argv[0], 0);
        break;

      case 'v': /* verbose */
        ++verbose;
        break;

      case 'd': /* set debug level */
        break;

      case 'n': /* sets number of points */
        numPoints = (UINT4) atoi( optarg );
        break;

      case 'i': /* sets number of points */
        inRate = (UINT4) atoi( optarg );
        break;

      case 'o': /* sets number of points */
        outRate = (UINT4) atoi( optarg );
        break;

      case 'f': /* sets number of points */
        sineFreq = (REAL4) atof( optarg );
        break;

      case 'r':
        if ( ! strcmp( "ldas", optarg ) )
        {
          filtType = LDASfirLP;
        }
        else if ( ! strcmp( "butterworth", optarg ) )
        {
          filtType = defaultButterworth;
        }
        else
        {
          Usage (argv[0], 1);
        }
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
