/* <lalVerbatim file="TimeDelayCV">

Author: Chin, David <dwchin@umich.edu> +1-734-709-9119, Kipp Cannon <kipp@gravity.phys.uwm.edu>
$Id$

</lalVerbatim> */

/* <lalLaTeX>

\subsection{Module \texttt{TimeDelay.c}}
\label{ss:TimeDelay.c}

Computes difference in arrival time of the same signal at two different
detectors. 

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TimeDelayCP}
\idx{LALTimeDelay()}


\subsubsection*{Description}

The function LALTimeDelay() computes the difference in time of arrival of a
signal at two detectors from the same source.  The two detectors and the source
are passed in a \texttt{TwoDetectorsAndASource} structure.  The time delay is
defined to be $\delta t = t_2 - t_1$ where $t_1$ is the time the signal arrives
at the first detector and $t_2$ is the time the signal arrives at the second
detector.

The function XLALLightTravelTime() computes the light travel time between two detectors and returns the answer in \texttt{INT8} nanoseconds.

\subsubsection*{Algorithm}

TBA. See Anderson, \textit{et al.} \cite{tools:Anderson:2000} in the mean time.

Note that GPS time is passed with both the detectors.  The GPS time of the
second detector is \emph{ignored}, and the GPS time for the first detector
is taken to be the time when the signal arrives at the center of the
Earth.  In practice, this time will be the time of detection of a signal at
the first detector, but, as in Anderson, \textit{et al.}, we make this
approximation as it makes little difference.  This time is used to compute
a GMST which gives us the orientation of the Earth.

\subsubsection*{Uses}


\subsubsection*{Notes}

\vfill{\footnotesize\input{TimeDelayCV}}

</lalLaTeX> */

#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/TimeDelay.h>

NRCSID( TIMEDELAYC, "$Id$" );

/* scalar product of two 3-vectors */
static double dotprod(const double vec1[3], const double vec2[3])
{
	return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}


/* <lalVerbatim file="TimeDelayCP"> */
double
XLALTimeDelay(
	const double detector1_earthfixed_xyz_metres[3],
	const double detector2_earthfixed_xyz_metres[3],
	double source_right_ascension_radians,
	double source_declination_radians,
	const LIGOTimeGPS *gpstime
)
{ /* </lalVerbatim> */
	static const char *func = "XLALTimeDelay";
	double delta_xyz[3];
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
	 * displacement of detector 1 from detector 2
	 */

	delta_xyz[0] = detector1_earthfixed_xyz_metres[0] - detector2_earthfixed_xyz_metres[0];
	delta_xyz[1] = detector1_earthfixed_xyz_metres[1] - detector2_earthfixed_xyz_metres[1];
	delta_xyz[2] = detector1_earthfixed_xyz_metres[1] - detector2_earthfixed_xyz_metres[2];

	/*
	 * Time difference: time taken for light to travel the intervening
	 * distance along the direction to the source.  If (detLoc1 -
	 * detLoc2) \cdot ehat_src is positive, then detector 1 is closer
	 * to the source and hence the time delay is positive.
	 */

	return dotprod(ehat_src, delta_xyz) / LAL_C_SI;
}


/* <lalVerbatim file="TimeDelayCP"> */
void
LALTimeDelay(
	LALStatus *stat,
	REAL8 *p_time_diff,
	const TwoDetsTimeAndASource *p_dets_time_and_source
)
{ /* </lalVerbatim> */
	INITSTATUS(stat, "LALTimeDelay", TIMEDELAYC);
	ATTATCHSTATUSPTR(stat);
	ASSERT(p_time_diff, stat, TIMEDELAYH_ENUL, TIMEDELAYH_MSGENUL);
	ASSERT(p_dets_time_and_source, stat, TIMEDELAYH_ENUL, TIMEDELAYH_MSGENUL);

	*p_time_diff = XLALTimeDelay(p_dets_time_and_source->p_det_and_time1->p_detector->location, p_dets_time_and_source->p_det_and_time2->p_detector->location, p_dets_time_and_source->p_source->longitude, p_dets_time_and_source->p_source->latitude, p_dets_time_and_source->p_det_and_time1->p_gps);

	ASSERT(!XLAL_IS_REAL8_FAIL_NAN(*p_time_diff), stat, DATEH_ERANGEGPSABS, DATEH_MSGERANGEGPSABS);

	DETATCHSTATUSPTR(stat);
	RETURN(stat);
}


/* <lalVerbatim file="TimeDelayCP"> */
INT8
XLALLightTravelTime(
	const LALDetector *aDet,
	const LALDetector *bDet
)
/* </lalVerbatim> */
{
	double deltaLoc[3];

	deltaLoc[0] = aDet->location[0] - bDet->location[0];
	deltaLoc[1] = aDet->location[1] - bDet->location[1];
	deltaLoc[2] = aDet->location[2] - bDet->location[2];

	return (INT8) (1e9 * sqrt(dotprod(deltaLoc, deltaLoc)) / LAL_C_SI);
}
