/* <lalVerbatim file="TimeDelayFromEarthCenterCV">

Author: David Chin <dwchin@umich.edu> +1-734-730-1274
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{TimeDelayFromEarthCenter.c}}
\label{ss:TimeDelayFromEarthCenter.c}

Computes difference in arrival time of the same signal at detector and
at center of Earth-fixed frame.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TimeDelayFromEarthCenterCP}
\idx{LALTimeDelayFromEarthCenter()}

\subsubsection*{Description}

This function computes the difference in time of arrival of a signal
at a given detector and at the center of the Earth-fixed frame from
the same source.  The detector and the source are passed in a
\verb@DetAndASource@ structure.  The time delay is defined to be
$\delta t = t_d - t_0$, where $t_d$ is the time the signal arrives at
the detector, and $t_0$ is the time the signal arrives at the center
of the Earth.

\subsubsection*{Algorithm}

The vector of the location of the detector is projected onto the line
connecting the center of the Earth and the source.  That projection
gives the distance of the detector from the center of the Earth along
the line connecting the center of the Earth and the source.  The time
delay, then, is the amount of time it takes light to travel that
projected distance.

The GPS time for the detector is taken to be the time when the signal
arrives at the center of the Earth.  As in Anderson, \textit{et al.}
\cite{ABCF:2000}, we make this approximation as this makes little
difference.

\subsubsection*{Uses}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TimeDelayFromEarthCenterCV}}

</lalLaTeX> */

#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/TimeDelay.h>

NRCSID( TIMEDELAYFROMEARTHCENTERC, "$ID$" );

/* scalar product of two 3-vectors */
static REAL8 dotprod(REAL8 vec1[3], REAL8 vec2[3])
{
  return (vec1[0] * vec2[0] +
          vec1[1] * vec2[1] +
          vec1[2] * vec2[2]);
}


/* <lalVerbatim file="TimeDelayFromEarthCenterCP"> */
void
LALTimeDelayFromEarthCenter( LALStatus               *stat,
                             REAL8                   *p_time_diff,
                             const DetTimeAndASource *p_det_time_and_source )
{ /* </lalVerbatim> */
  LALLeapSecAccuracy accuracy = LALLEAPSEC_STRICT;

  /* latitude and longitude for detector */
  REAL8 lat, lon;

  /* cos(lat) and sin(lat) for detector */
  REAL8 cosLat, sinLat;

  /* local radius of curvature at detector */
  REAL8 R;

  /* radius of curvature plus height above ellipsoid for detector */
  REAL8 Rh;

  /* GMST of time of arrival at detector, in radians */
  REAL8 gmst1;
  
  LALDate      date;


  /* NOTE: all source location params are in Earth-fixed frame */
  SkyPosition src_polar;      /* Earth-fixed polar (lon, polar angle) */
  REAL8       sin_pol_angle;  /* sine of src polar angle */
  REAL8       ehat_src[3];    /* unit vector of source location */

  /* location vector for the detector in Earth-fixed frame */
  REAL8 detLoc[3];
  REAL8 *p_detLoc;


  INITSTATUS( stat, "LALTimeDelayFromEarthCenter", TIMEDELAYFROMEARTHCENTERC );
  /* bloody mis-spelling */
  ATTATCHSTATUSPTR( stat );

  /* Ensure non-NULL structures are passed */
  ASSERT( p_time_diff, stat, TIMEDELAYH_ENUL, TIMEDELAYH_MSGENUL );
  ASSERT( p_det_time_and_source, stat, TIMEDELAYH_ENUL, TIMEDELAYH_MSGENUL );

  /*
   * convert src location in equatorial coordinates to Earth-fixed
   * polar coordinates
   *
   * don't use LALEquatorialToGeographic() since it modifies the input
   * structure
   */
  /* need GMST in radians */
  TRY( LALGPStoUTC( stat->statusPtr, &date,
                    p_det_time_and_source->p_det_and_time->p_gps,
                    &accuracy), stat );
  TRY( LALGMST1( stat->statusPtr, &gmst1, &date, MST_RAD ), stat );

  /* polar angle, theta */
  src_polar.latitude = LAL_PI_2 - p_det_time_and_source->p_source->latitude;

  /* azimuthal angle, phi */
  src_polar.longitude = p_det_time_and_source->p_source->longitude - gmst1;


  /*
   * compute the unit vector of the source direction
   */
  sin_pol_angle = sin(src_polar.latitude);
  
  ehat_src[0]   = sin_pol_angle * cos(src_polar.longitude);
  ehat_src[1]   = sin_pol_angle * sin(src_polar.longitude);
  ehat_src[2]   = cos(src_polar.latitude);


  /*
   * location vector of detector 1 
   */
  p_detLoc = p_det_time_and_source->p_det_and_time->p_detector->location;

  /*
   * time difference: time taken for light to travel the distance
   * between Earth center and detector along direction to source.
   * See LALTimeDelay(), and put in Earth-center for detector 1 to see
   * how the -ve sign arises.
   */
  *p_time_diff = -dotprod(ehat_src, p_detLoc) / LAL_C_SI;

  
  /*
   * House keeping
   */
  DETATCHSTATUSPTR( stat );
  RETURN( stat );
} /* END: LALTimeDelayFromEarthCenter */
