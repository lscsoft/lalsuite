/* <lalVerbatim file="LALTestDetResponseCV">
   $Id$
   </lalVerbatim
*/

/*
<lalLaTeX>

\subsection{Program {\texttt{LALTestDetResponse.c}}}
\label{ss:LALTestDetResponse.c}

\subsubsection*{Usage}

\subsubsection*{Description}

\subsubsection*{Exit codes}

\subsubsection*{Uses}
\begin{verbatim}
lalDebugLevel
LALComputeDetAMResponse()
\end{verbatim}

\subsubsection*{Notes}

</lalLaTeX> 
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <lal/LALConfig.h>

#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALStdlib.h>
#include <lal/DetectorSite.h>
#include <lal/TimeDelay.h>
#include <lal/DetResponse.h>

NRCSID( LALTESTDETRESPONSEC, "$Id$" );

/*
 * Make epsilon larger: GMST1 routine is only good to about 1.2e-5
 */
#define LALDR_REAL4_EPS 1.5e-05

void TestOnePoint(LALStatus              *status,
                  INT4                   *retval,
                  const LALFrDetector    *p_frDetector,
                  const LALSource        *p_source,
                  const LIGOTimeGPS      *p_gps,
                  const LALDetAMResponse *p_expectedResponse,
                  INT4                    use_LHO_p);

static REAL8 deg_to_rad(REAL8 deg);

int lalDebugLevel = 0;
int verbose_p       = 1;
              

int main(int argc, char *argv[])
{
  static LALStatus status;
  LALSource        pulsar;
  LALFrDetector    frdet;
  LIGOTimeGPS      gps;
  LALDate          date;
  LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;
  LALDetAMResponse expectedResponse;
  INT4             passed;
  INT4             testNumber;

  if (argc == 2)
    {
      lalDebugLevel = atoi(argv[1]);
      if (lalDebugLevel >= 8)
        verbose_p = 1;
      else
        verbose_p = 0;
    }

  if (verbose_p > 0)
    printf("LALDR_REAL4_EPS = %2.9e\n", LALDR_REAL4_EPS);

  /* Set source coord system to equatorial */
  pulsar.equatorialCoords.system = COORDINATESYSTEM_EQUATORIAL;

  /*****************
   * Test number 1
   *****************/
  testNumber = 1;

  /* source: RA=0, Dec=0, psi=0 */
  strcpy(pulsar.name, "TEST1: 0,0,0");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = 0.;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 2
   ****************/
  testNumber = 2;

  /* source: RA=0, Dec=0, psi=pi/2 */
  strcpy(pulsar.name, "TEST2: 0,0,pi/2");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = LAL_PI_2;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)-1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }

  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 3
   ****************/
  testNumber = 3;

  /* source: RA=0, Dec=0, psi=pi */
  strcpy(pulsar.name, "TEST3: 0,0,pi");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = LAL_PI;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 4
   ****************/
  testNumber = 4;

  /* source: RA=0, Dec=0, psi=3*pi/2 */
  strcpy(pulsar.name, "TEST4: 0,0,3*pi/2");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = LAL_PI + LAL_PI_2;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)-1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 5
   ****************/
  testNumber = 5;

  /* source: RA=0, Dec=0, psi=pi/4 */
  strcpy(pulsar.name, "TEST5: 0,0,pi/4");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = LAL_PI_4;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)0.;
  expectedResponse.cross  = (REAL4)-1.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 6
   ****************/
  testNumber = 6;

  /* source: RA=0, Dec=0, psi=pi/4 */
  strcpy(pulsar.name, "TEST6: 0,0,3*pi/4");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = 3.*LAL_PI_4;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)0.;
  expectedResponse.cross  = (REAL4)1.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  
  /****************
   * Test number 7
   ****************/
  testNumber = 7;

  /* source: RA=0, Dec=0, psi=5*pi/4 */
  strcpy(pulsar.name, "TEST7: 0,0,5*pi/4");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = 5.*LAL_PI_4;

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST1: (0, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)0.;
  expectedResponse.cross  = (REAL4)-1.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 8
   ****************/
  testNumber = 8;

  /* source: RA=36deg, dec=0, orientation=0 */
  strcpy(pulsar.name, "TEST8: 36deg, 0deg, 0");
  pulsar.equatorialCoords.longitude = deg_to_rad(36.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(0.);
  pulsar.orientation                = deg_to_rad(0.);

  /* detector: (36E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST3: (36, 0), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 36.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 9
   ****************/
  testNumber = 9;

  /* source: RA=0deg, dec=60deg, orientation=0 */
  strcpy(pulsar.name, "TEST9: 0deg, 60deg, 0");
  pulsar.equatorialCoords.longitude = deg_to_rad(0.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(60.);
  pulsar.orientation                = deg_to_rad(0.);

  /* detector: (0E, 60N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST3: (0, 60), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 60.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 10
   ****************/
  testNumber = 10;

  /* source: RA=30deg, dec=60, orientation=0 */
  strcpy(pulsar.name, "TEST10: 30deg, 60deg, 0");
  pulsar.equatorialCoords.longitude = deg_to_rad(30.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(60.);
  pulsar.orientation                = deg_to_rad(0.);

  /* detector: (36E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST4: (30, 60), (0, pi/2)");
  frdet.vertexLongitudeDegrees = 30.;
  frdet.vertexLatitudeDegrees  = 60.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 11
   ****************/
  testNumber = 11;

  /* source: RA=30deg, dec=60deg, orientation=-30deg */
  strcpy(pulsar.name, "TEST11: 30deg, 60deg, -30deg");
  pulsar.equatorialCoords.longitude = deg_to_rad(30.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(60.);
  pulsar.orientation                = deg_to_rad(-30.);

  /* detector: (36E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST5: (30, 60), (30, 120)");
  frdet.vertexLongitudeDegrees = 30.;
  frdet.vertexLatitudeDegrees  = 60.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = LAL_PI / 6.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2 + LAL_PI / 6.;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 12
   ****************/
  testNumber = 12;

  /* source: RA=30deg, dec=60deg, orientation=-75deg */
  strcpy(pulsar.name, "TEST12: 30deg, 60deg, -75deg");
  pulsar.equatorialCoords.longitude = deg_to_rad(30.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(60.);
  pulsar.orientation                = deg_to_rad(-75.);

  /* detector: (36E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST5: (30, 60), (30, 120)");
  frdet.vertexLongitudeDegrees = 30.;
  frdet.vertexLatitudeDegrees  = 60.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = LAL_PI / 6.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2 + LAL_PI / 6.;

  /* GPS time corresponding to 0h GMST1. (by trial and error) */
  /* puts source at zenith of detector */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  /* expected result */
  expectedResponse.plus   = (REAL4)0.;
  expectedResponse.cross  = (REAL4)1.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }


  /****************
   * Test number 13
   ****************/
  testNumber = 13;

  /* source: RA=27deg, dec=53.8deg, orientation=0deg */
  strcpy(pulsar.name, "TEST13: 27deg, -53.8deg, 0.");
  pulsar.equatorialCoords.longitude = deg_to_rad(27.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(-53.8);
  pulsar.orientation                = deg_to_rad(0.);

  /* detector: (27E, 53.8N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST6: (27, 53.8), (0,90)");
  frdet.vertexLongitudeDegrees = 27.;
  frdet.vertexLatitudeDegrees  = 53.8;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* time corresponding to 12h GMST1. (by trial and error) */
  date.unixDate.tm_year  = 94;
  date.unixDate.tm_mon   =  4;
  date.unixDate.tm_mday  = 17;
  date.unixDate.tm_hour  = 20;
  date.unixDate.tm_min   = 18;
  date.unixDate.tm_sec   = 48;
  date.residualNanoSeconds = 790048502;

  LALUTCtoGPS(&status, &gps, &date, &accuracy);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse: error in LALUTCtoGPS, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSEC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }

  /****************
   * Test number 14
   ****************/
  testNumber = 14;

  /* source: RA=27deg, dec=53.8deg, orientation=73.19deg */
  strcpy(pulsar.name, "TEST14: 27deg, -53.8deg, 73.19");
  pulsar.equatorialCoords.longitude = deg_to_rad(27.);
  pulsar.equatorialCoords.latitude  = deg_to_rad(-53.8);
  pulsar.orientation                = deg_to_rad(73.19);

  /* detector: (27E, 53.8N) x-arm: bearing 16.81; y-arm: bearing 286.81 */
  strcpy(frdet.name, "TEST7: (27, 53.8), (16.81,163.19)");
  frdet.vertexLongitudeDegrees = 27.;
  frdet.vertexLatitudeDegrees  = 53.8;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = deg_to_rad(73.19);  /* measured from E */
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2 + deg_to_rad(73.19);

  /* time corresponding to 12h GMST1. (by trial and error) */
  date.unixDate.tm_year  = 94;
  date.unixDate.tm_mon   =  4;
  date.unixDate.tm_mday  = 17;
  date.unixDate.tm_hour  = 20;
  date.unixDate.tm_min   = 18;
  date.unixDate.tm_sec   = 48;
  date.residualNanoSeconds = 790048502;

  LALUTCtoGPS(&status, &gps, &date, &accuracy);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse: error in LALUTCtoGPS, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSEC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /* expected result */
  expectedResponse.plus   = (REAL4)1.;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }
  

  /****************
   * Test number 15
   ****************/
  testNumber = 15;

  /* source */
  strcpy(pulsar.name, "TEST15: 90deg, 0deg, 0");
  pulsar.equatorialCoords.longitude = deg_to_rad(90.);
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = 0.;

  /* detector */
  strcpy(frdet.name, "TEST7");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* time corresponding to 12h GMST1. (by trial and error) */
  date.unixDate.tm_year  = 94;
  date.unixDate.tm_mon   =  4;
  date.unixDate.tm_mday  = 17;
  date.unixDate.tm_hour  = 20;
  date.unixDate.tm_min   = 18;
  date.unixDate.tm_sec   = 48;
  date.residualNanoSeconds = 790048502;

  LALUTCtoGPS(&status, &gps, &date, &accuracy);

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse: error in LALUTCtoGPS, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSEC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /* expected result */
  expectedResponse.plus   = (REAL4).5;
  expectedResponse.cross  = (REAL4)0.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }
  

  /****************
   * Test number 16
   ****************/
  testNumber = 16;

  /* source */
  strcpy(pulsar.name, "TEST16: 0deg, 0deg, -30deg");
  pulsar.equatorialCoords.longitude = 0.;
  pulsar.equatorialCoords.latitude  = 0.;
  pulsar.orientation                = deg_to_rad(-30.);

  /* detector: (0E, 0N) x-arm: bearing 090; y-arm: bearing 000 */
  strcpy(frdet.name, "TEST8");
  frdet.vertexLongitudeDegrees = 0.;
  frdet.vertexLatitudeDegrees  = 0.;
  frdet.vertexElevation        = 0.;
  frdet.xArmAltitudeRadians    = 0.;
  frdet.xArmAzimuthRadians     = 0.;
  frdet.yArmAltitudeRadians    = 0.;
  frdet.yArmAzimuthRadians     = LAL_PI_2;

  /* time */
  gps.gpsSeconds     = 61094;
  gps.gpsNanoSeconds = 640000000;

  if (status.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "LALTestDetResponse: error in LALUTCtoGPS, line %i, %s\n",
              __LINE__, LALTESTDETRESPONSEC);
      REPORTSTATUS(&status);
      return status.statusCode;
    }

  /* expected result */
  expectedResponse.plus   = (REAL4)0.5;
  expectedResponse.cross  = (REAL4)sqrt(3.)/2.;
  expectedResponse.scalar = (REAL4)0.;  /* this is ignored, for now */

  if (verbose_p > 0)
    {
      printf("Test number %d\n", testNumber);
      printf("---------------\n");
    }
  
  TestOnePoint(&status, &passed, &frdet, &pulsar, &gps, &expectedResponse, 0);

  if (passed != 1)
    {
      printf("LALTestDetResponse: ERROR: failed test no. %d\n", testNumber);
      return 1;
    }
  


  
  return 0;
}
              

/* returns 1 if computed values are as expected, 0 otherwise */
/* only tests differential mode IFOs, and only plus and cross components */
void TestOnePoint(LALStatus              *status,
                  INT4                   *retval,
                  const LALFrDetector    *p_frDetector,
                  const LALSource        *p_source,
                  const LIGOTimeGPS      *p_gps,
                  const LALDetAMResponse *p_expectedResponse,
                  INT4                    use_LHO_p)
{
  LALDetector       detector;
  LALDetAndSource   detAndSource;
  LALDetAMResponse  computedResponse;
  char              infostr[1024];

  INITSTATUS(status, "TestOnePoint", LALTESTDETRESPONSEC);
  ATTATCHSTATUSPTR(status);

  if (use_LHO_p == 0)
    {
      TRY( LALCreateDetector(status->statusPtr, &detector, p_frDetector,
                             LALDETECTORTYPE_IFODIFF), status );
      
      detAndSource.pDetector = &detector;
    }
  else
    {
      detAndSource.pDetector = &(lalCachedDetectors[LALDetectorIndexLHODIFF]);
    }
  
  detAndSource.pSource   = p_source;

  if (lalDebugLevel > 0 || verbose_p > 0)
    {
      sprintf(infostr,
              "Detector: (%2.9e deg. E, %2.9e deg. N)\n          (x-azi=%2.9e rad., y-azi=%2.9e rad.)",
              p_frDetector->vertexLongitudeDegrees,
              p_frDetector->vertexLatitudeDegrees,
              p_frDetector->xArmAzimuthRadians,
              p_frDetector->yArmAzimuthRadians);
      LALInfo(status, infostr);

      sprintf(infostr,
              "Source: (%2.9e RA, %2.9e Dec, %2.9e orientation) radians",
              p_source->equatorialCoords.longitude,
              p_source->equatorialCoords.latitude,
              p_source->orientation);
      LALInfo(status, infostr);
    }

  TRY( LALComputeDetAMResponse(status->statusPtr, &computedResponse,
                               &detAndSource, p_gps), status);

  if (lalDebugLevel > 0 || verbose_p > 0)
    {
      printf("expected: plus=%2.9e, cross=%2.9e  +/- %2.9e\n",
             p_expectedResponse->plus, p_expectedResponse->cross,
             LALDR_REAL4_EPS);
      printf("computed: plus=%2.9e, cross=%2.9e\n",
             computedResponse.plus, computedResponse.cross);
      printf("diff.:    plus=%2.9e, cross=%2.9e\n",
             (computedResponse.plus - p_expectedResponse->plus),
             (computedResponse.cross - p_expectedResponse->cross));
    }
  
  if ((REAL4)fabs((double)(computedResponse.plus - p_expectedResponse->plus))
      > LALDR_REAL4_EPS ||
      (REAL4)fabs((double)(computedResponse.cross - p_expectedResponse->cross))
      > LALDR_REAL4_EPS)
    {
      *retval = 0;
    }
  else
    {
      *retval = 1;
    }

  DETATCHSTATUSPTR(status);
  RETURN(status);
}


static REAL8 deg_to_rad(REAL8 deg)
{
  return deg / (REAL8)180. * LAL_PI;
}
