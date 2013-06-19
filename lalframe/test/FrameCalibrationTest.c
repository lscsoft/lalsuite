/*
*  Copyright (C) 2007 Bernd Machenschalk, Duncan Brown, Stephen Fairhurst
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
 * \author Brown, D. A.
 * \file
 *
 * \brief Tests the hig-level function to obtain an up-to-date calibration from
 * frame files.
 *
 * \heading{Usage}
 *
 * \code
 * FrameCalibrationTest
 * \endcode
 *
 * \heading{Usage}
 * \code
 * FrameCalibrationTest [options]
 * Options:
 *   -h         print this message
 *   -o         write calibration to file
 *   -v         verbose: print extra information
 * \endcode
 *
 * \heading{Description}
 *
 * For each interferometer in the array \c ifoCode the program does the
 * following: For each GPS time in the array \c calTime, an attempt is
 * made to generate a calibration for that time. It reads a frame cache file
 * named <tt>ifo-CAL_V03-gpstimerange.catalog</tt>, where \c ifo is the
 * interferometer obtained from the array \c ifoCode, \c gpstimerange is
 * obtained from the array \c cacheTime. Even elements of the array
 * \c calTime use the first element of the array \c cacheTime and odd
 * elements the second. This is done to allow the developer to test a range of
 * GPS times and catalogs quickly.
 *
 * The GPS times for which calibrations are attempted are 729331981, 729332039,
 * 729332040, 729332041, 729332099 and 800000000. The catalog times are
 * 729273600-734367600 and 729273600-734367600.
 *
 * The test uses that calibration data files
 * \code
 * H-CAL_FAC_V03-729273600-5094000.gwf
 * H-CAL_REF_V03-734073939-64.gwf
 * L-CAL_FAC_V03-729273600-5094000.gwf
 * L-CAL_REF_V03-731488397-64.gwf
 * \endcode
 *
 * FIXME: the following test is not yet implemented: waiting for Steve's
 * change to add the return of the actual values used.
 * The returned values of \f$\alpha\f$ and \f$\alpha\beta\f$ are checked against the
 * values obtained from the files S2 version 3 calibration factor files
 * <tt>H1AlphaBeta.txt</tt> and <tt>L1AlphaBeta.txt</tt> which can be found on the
 * calibration home page.
 *
 * \heading{Exit codes}
 * <table>
 * <tr><th>Code</th><th>Explanation</th></tr>
 * <tr><td>\c 0</td><td>Success, normal exit.</td></tr>
 * <tr><td>\c 1</td><td>Subroutine failed.</td></tr>
 * <tr><td>\c 77</td><td>Ignored failure: Test frame data not found.</td></tr>
 * </table>
 *
 * \heading{Uses}
 * \heading{Notes}
 *
*/


#include <stdio.h>
#include <math.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_GETOPT_H
#include <getopt.h>
#endif
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/PrintFTSeries.h>
#include <lal/LALFrStream.h>
#include <lal/FrameCalibration.h>
#include <lal/LALFrameL.h>
#include <lal/Calibration.h>
#include <lal/PrintFTSeries.h>

#define TESTSTATUS( pstat ) \
  if ( (pstat)->statusCode ) { REPORTSTATUS(pstat); return 1; } else ((void)0)

#define CLEARSTATUS( stat ) \
  LALFree( stat.statusPtr ); \
  stat.statusPtr = NULL; \
  stat.statusCode = 0;

#define CHANNEL "%s:LSC-AS_Q"
#define CAL_CATALOG "%s-CAL-V03-%s.cache"


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
  UINT4 i, j, k;
  const REAL4 tiny = 1e-6;

  CHAR  ifoCode[][3] = { "H1", "L1" };
  INT4  calTime[] = { 729331981, 729332039, 729332040, 729332041, 729332099,
    800000000 };
  CHAR  cacheTime[][21] = { "729273600-734367600", "729273600-734367600" };
  COMPLEX8 H1AlphaBeta[] = { {0.9883124570783, 0}, {0.9883124570783, 0},
    {1.123396433694, 0}, {1.123396433694, 0}, {1.123396433694, 0}, {0, 0} };
  COMPLEX8 L1AlphaBeta[] = { {0, 0}, {0, 0}, {0.6041572088741, 0},
    {0.6041572088741, 0}, {0.6041572088741, 0}, {0, 0} };

  static LALStatus      status;
  const CHAR            calCacheName[LALNameLength];
  LALCache             *calCache = NULL;
  UINT4                 numPoints = 262144;
  UINT4                 sampleRate = 4096;
  CHAR                  outFile[LALNameLength];
  LIGOTimeGPS           duration = {0,0};

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

  /* loop over the three interferometers */
  for ( j = 0; j < sizeof(ifoCode) / sizeof(*ifoCode); ++j )
  {
    snprintf( response.name, LALNameLength * sizeof(CHAR),
        CHANNEL, ifoCode[j] );

    for ( i = 0; i < sizeof(calTime) / sizeof(*calTime); ++i )
    {
      /* set the time of the calibration and the frame cahche file to use */
      snprintf( calCacheName, LALNameLength * sizeof(CHAR), CAL_CATALOG,
          ifoCode[j], cacheTime[i % 2] );
      response.epoch.gpsSeconds = calTime[i];
      if ( verbose )
      {
        fprintf( stdout, "Calibration for GPS time %d from %s\n",
            response.epoch.gpsSeconds, calCacheName );
        fflush( stdout );
      }

      /* create the response function */
      LALExtractFrameResponse( &status, &response, calCacheName, ifoCode[j],
          &duration );
      if ( status.statusCode == -1 && status.statusPtr )
      {
        if ( status.statusPtr->statusCode == FRAMESTREAMH_EDONE &&
            calTime[i] == 800000000 )
        {
          /* no calibration for this time */
          if ( verbose )
          {
            fprintf( stderr, "OK: No calibration for 800000000\n" );
            LALPrintError( "OK: %s\n", status.statusPtr->statusDescription );
          }
          CLEARSTATUS( status );
        }
        else if ( status.statusPtr->statusCode == CALIBRATIONH_EZERO &&
            ! strncmp( ifoCode[j], "L1", 2 * sizeof(CHAR) ) &&
            (calTime[i] == 729331981 || calTime[i] == 729332039) )
        {
          /* alpha is zero at this time */
          if ( verbose )
          {
            fprintf( stderr, "OK: Cal is zero for L1 at 729332039\n" );
            LALPrintError( "OK: %s\n", status.statusPtr->statusDescription );
          }
          CLEARSTATUS( status );
        }
        else
        {
          /* some other error */
          TESTSTATUS( &status );
        }
      }
      else
      {
        TESTSTATUS( &status );

        /* FIXME check that the values match the expected ones          */
        /* H1AlphaBeta[i] and L1AlphaBeta[i] give the correct values    */
        /* for this i, so we need to check that the returned values     */
        /* match the expected ones up to tiny. Need to if on j to get   */
        /* the correct ifo                                              */

        /* test that the response does not contain NaN or Inf */
        for ( k = 0; k < response.data->length; ++k )
        {
          if ( (! finite( response.data->data[k].re )) ||
              (! finite( response.data->data[k].im )) )
          {
            fprintf( stderr, "ERROR: non-finite value found in response "
                "at k = %d (%e,%e)\n",
                k, response.data->data[k].re, response.data->data[k].im );
            exit( 1 );
          }
        }

        /* print out the response function */
        if ( verbose )
        {
          fprintf( stdout, "Calibration updated\n" );
          fflush( stdout );
        }
        if ( output )
        {
          snprintf( outFile, LALNameLength * sizeof(CHAR),
              "Response-%s-%d.txt", ifoCode[j], response.epoch.gpsSeconds );
          LALCPrintFrequencySeries( &response, outFile );
        }
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
        lalDebugLevel |= LALWARNING | LALINFO;
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
