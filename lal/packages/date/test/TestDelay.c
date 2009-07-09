/*
*  Copyright (C) 2007 David Chin, Jolien Creighton
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

int main(int argc, char **argv)
{
  static LALStatus stat;
  LALFrDetector    frdet1;     /* Framelib detector info */
  LALFrDetector    frdet2;
  LALDetector      detector1;
  LALDetector      detector2;
  LIGOTimeGPS      gps;
  SkyPosition      source;
  REAL8            delay;
  REAL8            difference;
  DetTimeAndASource     det1_and_source;
  /* TwoDetsTimeAndASource dets_and_source; */
  LALPlaceAndGPS        det1_and_gps;
  /* LALPlaceAndGPS        det2_and_gps; */

  if (argc > 1)
    lalDebugLevel = atoi(argv[1]);

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
  if (lalDebugLevel > 2)
    REPORTSTATUS(&stat);

  /*
   * Expect the location vector to be (R, 0, 0): R = radius of Earth
   *                                                 at Equator
   * tolerance, 1.e-04
   */
  if (fabs(detector1.location[0] - LAL_REARTH_SI)/LAL_REARTH_SI > 1.e-4 ||
      detector1.location[1] != 0.                         ||
      detector1.location[2] != 0.)
    {
      fprintf(stderr, "TestDelay: LALCreateDetector output is wrong, line %i, %s\n",
              __LINE__, LALTESTDELAYC);
      fprintf(stderr, "Got Det #1 location: (% 16.8e, % 16.8e, % 16.8e)\n",
              (float)detector1.location[0], (float)detector1.location[1],
              (float)detector1.location[2]);
      fprintf(stderr, "Expected:            (% 16.8e, % 16.8e, % 16.8e)\n",
              (float)LAL_REARTH_SI, 0., 0.);

      return(1);
    }

  if (lalDebugLevel > 2)
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
  if (lalDebugLevel > 2)
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
  if (lalDebugLevel > 2)
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

  if (lalDebugLevel > 2)
    {
      printf("delay      = %20.14e\n", delay);
      printf("Rearth / c = %20.14e\n", (REAL8)LAL_REARTH_SI /
             (REAL8)LAL_C_SI);
    }

  difference = fabs(delay) - (REAL8)LAL_REARTH_SI / (REAL8)LAL_C_SI;

  if (difference < DOUBLE_EPSILON)
    {
      return 0;
    }
  else
    {
      fprintf(stderr, "ERROR: computed delay differs from expected delay by amount greater than DOUBLE_EPSILON (% 14.8e); difference = % 14.8e\n",
              DOUBLE_EPSILON, difference);
      return 1;
    }
}
