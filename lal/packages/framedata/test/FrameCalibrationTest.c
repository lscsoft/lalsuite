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
 * \subsubsection*{Description}
 *
 * For each GPS time in the array \verb+calTimes+, this program attemps to
 * generate a calibration updated to that time. It attemps to read a frame
 * cahce file names \verb+L-CAL-gpstime.catalog+, where \verb+gpstime+ is 
 * obtained from the array \verb+catalogTime+. Even elements of the arrau
 * \verb+calTimes+ use the first element of the array \verb+catalogTime+
 * and odd elements the second. This is done to allow the developer to test a
 * range of GPS times and catalogs quickly.
 *
 * The GPS times for which calibrations are attempted are 714313260,
 * 729925000, 714962262, 729925450 and 800000000. The catalog times are
 * 715388533 and 729907747.
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
 * 
 * Standard error codes from the calibration function are checked for and
 * handled, and if an updated response is generated it is written to a file.
 * 
 **** </lalLaTeX> */

#include <FrameL.h>

#include <stdio.h>
#include <unistd.h>
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
#define CAL_CATALOG "L-CAL-%d.catalog"

INT4 lalDebugLevel = LALWARNING | LALINFO;

int main( void )
{
  UINT4 i;

  INT4  calTime[] = { 714313260, 729925000, 714962262, 729925450, 800000000 };
  INT4  catalogTime[] = { 715388533, 729907747 };

  static LALStatus      status;
  const CHAR            calCacheName[LALNameLength];
  FrCache              *calCache = NULL;
  UINT4                 numPoints = 262144;
  UINT4                 sampleRate = 4096;
  CHAR                  ifo[3];
  CHAR                  outFile[LALNameLength];

  COMPLEX8FrequencySeries       response;
  const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};

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
        catalogTime[i % 2] );
    response.epoch.gpsSeconds = calTime[i];
    fprintf( stdout, "Getting calibration for GPS time %d from catalog %s\n", 
        response.epoch.gpsSeconds, calCacheName );
    fflush( stdout );

    /* create the response function */
    LALExtractFrameResponse( &status, &response, calCacheName, ifo );
    if ( status.statusCode == FRAMECALIBRATIONH_EMCHE )
    {
      LALPrintError( "No calibration cache file available\n" );
      return 77;
    }
    else if ( status.statusCode == FRAMECALIBRATIONH_ECREF )
    {
      LALPrintError( "Calibration NOT generated: %s\n", 
          status.statusDescription );
    }
    else if ( status.statusCode == FRAMECALIBRATIONH_ECFAC )
    {
      LALPrintError( "Calibration NOT updated: %s\n", 
          status.statusDescription );
    }
    else if ( status.statusCode == -1 && status.statusPtr &&
        ( status.statusPtr->statusCode == CALIBRATIONH_ETIME ||
          status.statusPtr->statusCode == CALIBRATIONH_EZERO ||
          status.statusPtr->statusCode == FRAMESTREAMH_EDONE ) )
    {
      LALPrintError( "Calibration NOT updated: %s\n", 
          status.statusPtr->statusDescription );
      LALFree( status.statusPtr );
      memset( &status, 0, sizeof( status ) );
    }
    else
    {
      TESTSTATUS( &status );

      /* print out the response function */
      fprintf( stdout, "Calibration updated\n" );
      fflush( stdout );
      LALSnprintf( outFile, LALNameLength * sizeof(CHAR),
          "Response-%s-%d.txt", ifo, response.epoch.gpsSeconds );
      LALCPrintFrequencySeries( &response, outFile );
    }
  }

  /* free memory */
  LALCDestroyVector( &status, &(response.data) );
  TESTSTATUS( &status );

  LALCheckMemoryLeaks();
  return 0;
}
