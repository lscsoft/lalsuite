/*
<lalVerbatim file="TestDelayCV">

Author: Chin, David <dwchin@umich.edu> +1-734-730-1274
$Id$

</lalVerbatim>
*/

/*
<lalLaTeX>


\subsection{Program \texttt{TestDelay.c}}
\label{ss:TestDelay.c}

Tests \verb@TimeDelay@ code.


\subsubsection*{Usage}

\begin{verbatim}
TestDelay
\end{verbatim}


\subsubsection*{Description}

This program does zero-th order tests for \texttt{LALTimeDelay()}. 


\subsubsection*{Exit codes}

</lalLaTeX>
*/

/*
<lalErrorTable>
*/

/*
</lalErrorTable>
*/

#include <math.h>
#include <stdlib.h>
#include <errno.h>
/* Darwin doesn't have values.h; the machine constants are defined in
 * float.h */
/* #include <values.h> */
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lal/SkyCoordinates.h>
#include <lal/DetectorSite.h>


NRCSID( LALTESTDELAYC, "$Id$" );

/* This should already be defined as X_EPS in /usr/include/values.h ;
 * in Darwin, it's defined as DBL_EPSILON in /usr/include/float.h */
#define DOUBLE_EPSILON 1.0536712127723507013e-08


int lalDebugLevel = 0;

int main( void )
{
  static LALStatus stat;
  LALFrDetector    frdet1;     /* Framelib detector info */
  LALFrDetector    frdet2;
  LALDetector      detector1;
  LALDetector      detector2;
  LIGOTimeGPS      gps;
  SkyPosition      source;
  REAL8            delay;
  DetTimeAndASource     det1_and_source;
  /* TwoDetsTimeAndASource dets_and_source; */
  LALPlaceAndGPS        det1_and_gps;
  /* LALPlaceAndGPS        det2_and_gps; */

  /*
   * Set up a source that will be used in both LALTimeDelay() and
   * LALTimeDelayFromEarthCenter().
   * Simple source at (RA=0, Dec=0)
   */
  source.longitude = 0.;
  source.latitude  = 0.;
  source.system    = COORDINATESYSTEM_EQUATORIAL;

  /*
   * Now, setup two detectors. One at (0.E, 0.N) and the other at (90.E,
   * 0.N)
   */
  strcpy(frdet1.name, "TEST IFO 1");
  frdet1.vertexLongitudeRadians = 0.;
  frdet1.vertexLatitudeRadians  = 0.;
  frdet1.vertexElevation        = 0.;
  frdet1.xArmAltitudeRadians    = 0.;
  frdet1.xArmAzimuthRadians     = 0.;
  frdet1.yArmAltitudeRadians    = LAL_PI_2;
  frdet1.yArmAzimuthRadians     = 0.;

  LALCreateDetector(&stat, &detector1, &frdet1, LALDETECTORTYPE_IFODIFF);
  if (stat.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestDelay: LALCreateDetector failed, line %i, %s\n",
              __LINE__, LALTESTDELAYC);
      REPORTSTATUS(&stat);
      return stat.statusCode;
    }
  REPORTSTATUS(&stat);

  /*
   * Expect the location vector to be (R, 0, 0): R = radius of Earth
   *                                                 at Equator
   */
  printf("Det #1 location: (%7.4e, %7.4e, %7.4e)\n",
         detector1.location[0], detector1.location[1],
         detector1.location[2]);
    

  strcpy(frdet2.name, "TEST IFO 2");
  frdet2.vertexLongitudeRadians = LAL_PI_2;
  frdet2.vertexLatitudeRadians  = 0.;
  frdet2.vertexElevation        = 0.;
  frdet2.xArmAltitudeRadians    = 0.;
  frdet2.xArmAzimuthRadians     = 0.;
  frdet2.yArmAltitudeRadians    = 0.;
  frdet2.yArmAzimuthRadians     = LAL_PI_2;

  LALCreateDetector(&stat, &detector2, &frdet2, LALDETECTORTYPE_IFODIFF);
  if (stat.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr, "TestDelay: LALCreateDetector failed, line %i, %s\n",
              __LINE__, LALTESTDELAYC);
      REPORTSTATUS(&stat);
      return stat.statusCode;
    }  
  REPORTSTATUS(&stat);

  /*
   * Set a GPS time that's close to 0h GMST1. (Found this by trial and
   * error.) 
   */
  gps.gpsSeconds     = 60858;
  gps.gpsNanoSeconds = 0;

  det1_and_gps.p_detector = &detector1;
  det1_and_gps.p_gps      = &gps;

  det1_and_source.p_det_and_time = &det1_and_gps;
  det1_and_source.p_source       = &source;

  LALTimeDelayFromEarthCenter(&stat, &delay, &det1_and_source);
  if (stat.statusCode && lalDebugLevel > 0)
    {
      fprintf(stderr,
              "TestDelay: LALTimeDelayFromEarthCenter() failed, line %i, %s\n",
              __LINE__, LALTESTDELAYC);
      REPORTSTATUS(&stat);
      return stat.statusCode;
    }
  REPORTSTATUS(&stat);

  /*
   * Expect delay to be roughly c/R, where c=speed of light,
   *                                       R=radius of Earth at
   *                                         Equator
   */
  /*
    printf("Time delay from Earth center = %18.13e sec\n", delay);
    printf("R/c = %18.13e sec\n", (REAL8)LAL_REARTH_SI / (REAL8)LAL_C_SI);

    printf("Diff = %18.13e\n", delay + (REAL8)LAL_REARTH_SI / (REAL8)LAL_C_SI);
    printf("X_EPS = %18.13e\n", (float)X_EPS);
    printf("H_PREC = %18.13e\n", (float)H_PREC);
  */

  printf("delay      = %20.14e\n", delay);
  printf("Rearth / c = %20.14e\n", (REAL8)LAL_REARTH_SI /
         (REAL8)LAL_C_SI);

  if ((fabs(delay) - (REAL8)LAL_REARTH_SI / (REAL8)LAL_C_SI) < DOUBLE_EPSILON)
    return 0;
  else
    return 1;
}
