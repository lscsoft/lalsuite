/*
<lalVerbatim file="TimeDelayCV">

Author: Chin, David <dwchin@umich.edu> +1-734-730-1274
$Id$

</lalVerbatim>
*/

/*
<lalLaTeX>

\subsection{Module \texttt{TimeDelay.c}}
\label{sec:TimeDelay.c}

Computes difference in arrival time of the same signal at two different
detectors. 


\subsubsection*{Prototypes}
\input{TimeDelayCP}
\index{\texttt{LALTimeDelay()}}


\subsubsection*{Description}

This function computes the difference in time of arrival of a signal at two
detectors from the same source.  The two detectors and the source are
passed in a \verb@TwoDetectorsAndASource@ structure.  The time delay is
defined to be $\delta t = t_2 - t_1$ where $t_1$ is the time the signal
arrives at the first detector and $t_2$ is the time the signal arrives at
the second detector.


\subsubsection*{Algorithm}

TBA. See Anderson, \emph{et al.} in the mean time.

Note that GPS time is passed with both the detectors.  The GPS time of the
second detector is \emph{ignored}, and the GPS time for the first detector
is taken to be the time when the signal arrives at the center of the
Earth.  In practice, this time will be the time of detection of a signal at
the first detector, but, as in Anderson, \emph{et al.}, we make this
approximation as it makes little difference.  This time is used to compute
a GMST which gives us the orientation of the Earth.

\subsubsection*{Uses}


\subsubsection*{Notes}

\vfill{\footnotesize\input{TimeDelayCV}}

</lalLaTeX>
*/

#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/TimeDelay.h>

NRCSID( TIMEDELAYC, "$ID$" );

/* scalar product of two 3-vectors */
static REAL8 dotprod(REAL8 vec1[3], REAL8 vec2[3])
{
  return (vec1[0] * vec2[0] +
          vec1[1] * vec2[1] +
          vec1[2] * vec2[2]);
}

/*
<lalVerbatim file="TimeDelayCP">
*/
void
LALTimeDelay( LALStatus                    *stat,
              REAL8                        *p_time_diff,
              const TwoDetsTimeAndASource  *p_dets_time_and_source )
{ /* </lalVerbatim> */
  /* a and b are WGS-84 ellipsoid parameters (semimajor and semiminor axes) */
  const REAL8 a2 = LAL_AWGS84_SI * LAL_AWGS84_SI;
  const REAL8 b2 = LAL_BWGS84_SI * LAL_BWGS84_SI;

  /* latitude and longitude for each detector */
  REAL8 lat1, lon1;
  REAL8 lat2, lon2;

  /* cos(lat) and sin(lat) for each detector */
  REAL8 cosLat1, sinLat1;
  REAL8 cosLat2, sinLat2;

    /* local radius of curvature at each detector */
  REAL8 R1;
  REAL8 R2;

  /* radius of curvature plus height above ellipsoid for each detector */
  REAL8 Rh1;
  REAL8 Rh2;

  /* GMST of first detector, in radians */
  REAL8 gmst1;
  
  LIGOTimeUnix unixTime; /* Unix time */
  LALDate      date;


  /* NOTE: all source location params are in Earth-fixed frame */
  SkyPosition src_polar;      /* Earth-fixed polar (lon, polar angle) */
  REAL8       sin_pol_angle;  /* sine of src polar angle */
  REAL8       ehat_src[3];    /* unit vector of source location */

  /* location vectors for the two detectors in Earth-fixed frame */
  REAL8 detLoc1[3];
  REAL8 detLoc2[3];
  REAL8 deltaLoc[3];

  /* loop counter */
  INT4  i;


  INITSTATUS( stat, "LALTimeDelay", TIMEDELAYC );
  /* bloody mis-spelling */
  ATTATCHSTATUSPTR( stat );

  /* Ensure non-NULL structures are passed */
  ASSERT( p_time_diff, stat, TIMEDELAYH_ENUL, TIMEDELAYH_MSGENUL );
  ASSERT( p_dets_time_and_source, stat, TIMEDELAYH_ENUL, TIMEDELAYH_MSGENUL );


  /*
   * convert src location in equatorial coordinates to Earth-fixed
   * polar coordinates
   *
   * don't use LALEquatorialToGeographic() since it modifies the input
   * structure
   */
  
  /* need GMST in radians */
  TRY( LALGPStoU( stat->statusPtr, &unixTime,
                  p_dets_time_and_source->det_and_time1->gps ), stat );
  TRY( LALUtime( stat->statusPtr, &date, &unixTime ), stat );
  TRY( LALGMST1( stat->statusPtr, &gmst1, &date, MST_RAD ), stat );

  /* polar angle, theta */
  src_polar.latitude = LAL_PI_2 - p_dets_time_and_source->source->latitude;

  /* azimuthal angle, phi */
  src_polar.longitude = p_dets_time_and_source->source->longitude - gmst1;

  /*
   * compute the unit vector of the source direction
   */
  sin_pol_angle = sin(src_polar.latitude);
  
  ehat_src[0]   = sin_pol_angle * cos(src_polar.longitude);
  ehat_src[1]   = sin_pol_angle * sin(src_polar.longitude);
  ehat_src[2]   = cos(src_polar.latitude);


  /*
   * compute the location vector of detector 1 
   */
  /* latitude of detector 1, in radians */
  lat1 = p_dets_time_and_source->det_and_time1->detector->frDetector.vertexLatitudeDegrees * LAL_PI / (REAL8)180.;
  
  /* longitude of detector 1, in radians */
  lon1 = p_dets_time_and_source->det_and_time1->detector->frDetector.vertexLongitudeDegrees * LAL_PI / (REAL8)180.;

  cosLat1 = cos(lat1);
  sinLat1 = sin(lat1);

  /* local rad. of curv. at detector 2 */
  R1  = a2 / sqrt(a2 * cosLat1 * cosLat1 +
                  b2 * sinLat1 * sinLat1); 
  Rh1 = R1 + p_dets_time_and_source->det_and_time1->detector->frDetector.vertexElevation;

  detLoc1[0] = Rh1 * cosLat1 * cos(lon1);
  detLoc1[1] = Rh1 * cosLat1 * sin(lon1);
  detLoc1[2] = (b2 * R1 / a2) * sinLat1;

  
  /*
   * compute the location vector of detector 2
   */  
  /* latitude of detector 2, in radians */
  lat2 = p_dets_time_and_source->det_and_time2->detector->frDetector.vertexLatitudeDegrees * LAL_PI / (REAL8)180.;
  
  /* longitude of detector 2, in radians */
  lon2 = p_dets_time_and_source->det_and_time2->detector->frDetector.vertexLongitudeDegrees * LAL_PI / (REAL8)180.;

  cosLat2 = cos(lat2);
  sinLat2 = sin(lat2);

  /* local rad. of curv. at detector 2*/
  R2  = a2 / sqrt(a2 * cosLat2 * cosLat2 +
                  b2 * sinLat2 * sinLat2); 
  Rh2 = R2 + p_dets_time_and_source->det_and_time2->detector->frDetector.vertexElevation;

  detLoc2[0] = Rh2 * cosLat2 * cos(lon2);
  detLoc2[1] = Rh2 * cosLat2 * sin(lon2);
  detLoc2[2] = (b2 * R2 / a2) * sinLat2;

  /*
   * displacement of detector 2 from detector 1
   */
  for (i = 0; i < 3; ++i)
    deltaLoc[i] = detLoc2[i] - detLoc1[i];

  /*
   * time difference: time taken for light to travel the intervening
   * distance along the direction to the source
   */
  *p_time_diff = dotprod(ehat_src, deltaLoc) / LAL_C_SI;


  /*
   * House keeping
   */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
}
