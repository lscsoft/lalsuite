/**** <lalVerbatim file="FrameCalibrationTestCV">
 * Author: Brown, D. A.
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 * \subsection{Program \texttt{FrameCalibrationTest.c}}
 * 
 * Tests the hig-level function to obtain an up-to-date calibration from
 * frame files.
 *
 * \subsubsection*{Usage}
 *
 * \begin{verbatim}
 * FrameCalibrationTest
 * \end{verbatim}
 *
 * \subsection*{Usage}
 * \begin{verbatim}
 * FrameCalibrationTest [options]
 * Options:
 *   -h         print this message
 *   -o         write calibration to file
 *   -v         verbose: print extra information
 * \end{verbatim}
 *
 * \subsubsection*{Description}
 *
 * For each GPS time in the array \verb+calTime+, this program attemps to
 * generate a calibration for that time. It reads a frame cahce file
 * named \verb+L-CAL-gpstimerange.catalog+, where \verb+gpstimerange+ is
 * obtained from the array \verb+cacheTime+. Even elements of the array
 * \verb+calTime+ use the first element of the array \verb+cacheTime+ and odd
 * elements the second. This is done to allow the developer to test a range of
 * GPS times and catalogs quickly.
 *
 * The GPS times for which calibrations are attempted are 714313260,
 * 729925000, 714962262, 729925450 and 800000000. The catalog times are
 * 714240000-715609980 and 729925200-729928800.
 *
 * Since this test requires LIGO science data, the required calibration 
 * frame files are not distributed with LAL, however they can be obtained 
 * via rsync from Caltech. The files required are:
 * \begin{verbatim}
 * L-CAL_REF-715388533-64.gwf
 * L-CAL_FAC-714240000-1369980.gwf
 * L-CAL_REF-729907747-64.gwf
 * L-SenseMonitor_L1_M-729925200-3600.gwf
 * \end{verbatim}
 * In order to run the test, these files should be placed under the framedata
 * test directory before running make check.
 * 
 * If the required calibration data is missing, then the test is not executed.
 *
 * \subsubsection*{Exit codes}
 * \begin{tabular}{|c|l|}
 * \hline
 *  Code & Explanation                   \\
 * \hline
 * \tt 0 & Success, normal exit.         \\
 * \tt 1 & Subroutine failed.            \\
 * \tt77 & Ignored failure: Test frame data not found. \\
 * \hline
 * \end{tabular}
 *
 * \subsubsection*{Uses}
 * \subsubsection*{Notes}
 * 
 **** </lalLaTeX> */


#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif
#include <FrameL.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/FrameStream.h>
#include <lal/FrameCalibration.h>
#include <lal/Calibration.h>
#include <lal/PrintFTSeries.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

#define CHANNEL "L1:LSC-AS_Q"
#define CAL_CATALOG "L-CAL-%s.cache"

INT4 lalDebugLevel =  LALINFO;

extern char *optarg;
extern int   optind;
int verbose = 0;
int output = 0;

static void
Usage (const char *program, int exitflag);

static void
ParseOptions (int argc, char *argv[]);

int main( int argc, char *argv[] )
{
  UINT4 i;

  INT4  calTime[] = { 714241541, 729925000, 714962262, 729925450, 800000000 };
  CHAR  cacheTime[][21] = { "714240000-715609980", "729925200-729928800" };

  static LALStatus      status;
  const CHAR            calCacheName[LALNameLength];
  FrCache              *calCache = NULL;
  UINT4                 numPoints = 262144;
  UINT4                 sampleRate = 4096;
  CHAR                  ifo[3];
  CHAR                  outFile[LALNameLength];

  COMPLEX8FrequencySeries       response;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

  ParseOptions (argc, argv);

  /* clear the response function and create storage for the frequency series */
  memset( &response, 0, sizeof(COMPLEX8FrequencySeries) );
  LALCCreateVector( &status, &(response.data), numPoints / 2 + 1 );
  TESTSTATUS( &status );

  /* set the parameters of the response function to generate */
  response.deltaF = (REAL8) sampleRate / (REAL8) numPoints;
  response.sampleUnits = strainPerCount;
  strncpy( response.name, CHANNEL, LALNameLength * sizeof(CHAR) );

  /* determine the ifo from the channel name */
  strncpy( ifo, response.name, 2 * sizeof(CHAR) );
  ifo[2] = '\0';

  for ( i = 0; i < sizeof(calTime) / sizeof(*calTime); ++i )
  {
    /* set the time of the calibration and the frame cahche file to use */
    LALSnprintf( calCacheName, LALNameLength * sizeof(CHAR), CAL_CATALOG,
        cacheTime[i % 2] );
    response.epoch.gpsSeconds = calTime[i];
    if ( verbose )
    {
      fprintf( stdout, "Calibration for GPS time %d from %s\n",
          response.epoch.gpsSeconds, calCacheName );
      fflush( stdout );
    }

    /* create the response function */
    LALExtractFrameResponse( &status, &response, calCacheName, ifo );
    if ( status.statusCode == FRAMECALIBRATIONH_EOREF )
    {
      if ( verbose )
      {
        LALPrintError( "%s\n", status.statusDescription );
      }
      return 77;
    }
    else if ( status.statusCode == FRAMECALIBRATIONH_ECREF ||
        status.statusCode == FRAMECALIBRATIONH_ECFAC )
    {
      if ( verbose )
      {
        LALPrintError( " %s\n", status.statusDescription );
      }
    }
    else if ( status.statusCode == -1 && status.statusPtr &&
        ( status.statusPtr->statusCode == CALIBRATIONH_ETIME ||
          status.statusPtr->statusCode == CALIBRATIONH_EZERO ||
          status.statusPtr->statusCode == FRAMESTREAMH_EDONE ) )
    {
      if ( verbose )
      {
        LALPrintError( "%s\n", status.statusPtr->statusDescription );
      }
      LALFree( status.statusPtr );
      status.statusPtr = NULL;
    }
    else
    {
      TESTSTATUS( &status );

      /* print out the response function */
      if ( verbose )
      {
        fprintf( stdout, "Calibration updated\n" );
        fflush( stdout );
      }
      if ( output )
      {
        LALSnprintf( outFile, LALNameLength * sizeof(CHAR),
            "Response-%s-%d.txt", ifo, response.epoch.gpsSeconds );
        LALCPrintFrequencySeries( &response, outFile );
      }
    }
  }

  /* free memory */
  LALCDestroyVector( &status, &(response.data) );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}

static void
Usage (const char *program, int exitcode)
{
  fprintf (stderr, "Usage: %s [options]\n", program);
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "  -h         print this message\n");
  fprintf (stderr, "  -o         write calibration to file\n");
  fprintf (stderr, "  -v         verbose: print extra information\n");
  exit (exitcode);
}

static void
ParseOptions (int argc, char *argv[])
{
  while (1)
  {
    int c = -1;
    c = getopt (argc, argv, "hvo");
    if (c == -1)
    {
      break;
    }
    switch (c)
    {
      case 'o': /* sets flag to write output files */
        ++output;
        break;

      case 'v': /* verbose */
        ++verbose;
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
