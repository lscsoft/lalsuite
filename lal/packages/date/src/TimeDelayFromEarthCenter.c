/* <lalVerbatim file="TimeDelayFromEarthCenterCV">

Author: David Chin <dwchin@umich.edu> +1-734-709-9119
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
\cite{tools:Anderson:2000}, we make this approximation as this makes little
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

NRCSID(TIMEDELAYFROMEARTHCENTERC, "$ID$");

/* scalar product of two 3-vectors */
static double dotprod(const double vec1[3], const double vec2[3])
{
	return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}


/* <lalVerbatim file="TimeDelayFromEarthCenterCP"> */
double XLALTimeDelayFromEarthCenter(
	const double detector_earthfixed_xyz_metres[3],
	double source_declination_radians,
	double source_right_ascension_radians,
	const LIGOTimeGPS *gpstime
)
{/* </lalVerbatim> */
	static const char *func = "XLALTimeDelayFromEarthCenter";
	double ehat_src[3];
	const double colatitude = LAL_PI_2 - source_declination_radians;
	const double greenwich_hour_angle = XLALGreenwichMeanSiderealTime(gpstime) - source_right_ascension_radians;

	if(XLAL_IS_REAL8_FAIL_NAN(greenwich_hour_angle))
		XLAL_ERROR_REAL8(func, XLAL_EFUNC);

	/*
	 * compute the unit vector pointing from the geocenter to the
	 * source
	 */

	ehat_src[0] = sin(colatitude) * cos(greenwich_hour_angle);
	ehat_src[1] = sin(colatitude) * -sin(greenwich_hour_angle);
	ehat_src[2] = cos(colatitude);

	/*
	 * time difference: time taken for light to travel the distance
	 * between Earth center and detector along direction to source.
	 * See LALTimeDelay(), and put in Earth-center for detector 1 to
	 * see how the sign arises.
	 */

	return -dotprod(ehat_src, detector_earthfixed_xyz_metres) / LAL_C_SI;
}


/* <lalVerbatim file="TimeDelayFromEarthCenterCP"> */
void LALTimeDelayFromEarthCenter(
	LALStatus *stat,
	REAL8 *p_time_diff,
	const DetTimeAndASource *p_det_time_and_source
)
{/* </lalVerbatim> */
	INITSTATUS(stat, "LALTimeDelayFromEarthCenter", TIMEDELAYFROMEARTHCENTERC);
	ATTATCHSTATUSPTR(stat);
	ASSERT(p_time_diff, stat, TIMEDELAYH_ENUL, TIMEDELAYH_MSGENUL);
	ASSERT(p_det_time_and_source, stat, TIMEDELAYH_ENUL, TIMEDELAYH_MSGENUL);

	*p_time_diff = XLALTimeDelayFromEarthCenter(p_det_time_and_source->p_det_and_time->p_detector->location, p_det_time_and_source->p_source->latitude, p_det_time_and_source->p_source->longitude, p_det_time_and_source->p_det_and_time->p_gps);

	DETATCHSTATUSPTR(stat);
	RETURN(stat);
}
