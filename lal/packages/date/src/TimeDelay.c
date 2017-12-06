/*
 * Copyright (C) 2007  Jolien Creighton, and David Chin, and Steven
 * Fairhurst, and Kipp Cannon, and Alexander Dietz, and Drew Keppel
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */

#include <math.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/SkyCoordinates.h>
#include <lal/TimeDelay.h>
#include <lal/XLALError.h>


/* scalar product of two 3-vectors */
static double dotprod(const double vec1[3], const double vec2[3])
{
	return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
}

/** UNDOCUMENTED */
double
XLALArrivalTimeDiff(
	const double detector1_earthfixed_xyz_metres[3],
	const double detector2_earthfixed_xyz_metres[3],
	const double source_right_ascension_radians,
	const double source_declination_radians,
	const LIGOTimeGPS *gpstime
)
{
	double delta_xyz[3];
	double ehat_src[3];
	const double greenwich_hour_angle = XLALGreenwichMeanSiderealTime(gpstime) - source_right_ascension_radians;

	if(XLAL_IS_REAL8_FAIL_NAN(greenwich_hour_angle))
		XLAL_ERROR_REAL8(XLAL_EFUNC);

	/*
	 * compute the unit vector pointing from the geocenter to the
	 * source
	 */

	ehat_src[0] = cos(source_declination_radians) * cos(greenwich_hour_angle);
	ehat_src[1] = cos(source_declination_radians) * -sin(greenwich_hour_angle);
	ehat_src[2] = sin(source_declination_radians);

	/*
	 * position of detector 2 with respect to detector 1
	 */

	delta_xyz[0] = detector2_earthfixed_xyz_metres[0] - detector1_earthfixed_xyz_metres[0];
	delta_xyz[1] = detector2_earthfixed_xyz_metres[1] - detector1_earthfixed_xyz_metres[1];
	delta_xyz[2] = detector2_earthfixed_xyz_metres[2] - detector1_earthfixed_xyz_metres[2];

	/*
	 * Arrival time at detector 1 - arrival time at detector 2.  This
	 * is positive when the wavefront arrives at detector 1 after
	 * detector 2 (and so t at detector 1 is greater than t at detector
	 * 2).
	 */

	return dotprod(ehat_src, delta_xyz) / LAL_C_SI;
}


/**
 * Compute difference in arrival time of the same signal at detector and at center of Earth-fixed frame.
 */
double XLALTimeDelayFromEarthCenter(
	const double detector_earthfixed_xyz_metres[3],
	double source_right_ascension_radians,
	double source_declination_radians,
	const LIGOTimeGPS *gpstime
)
{
	static const double earth_center[3] = {0.0, 0.0, 0.0};

	/*
	 * This is positive when the wavefront arrives at the detector
	 * after arriving at the geocentre.
	 */

	return XLALArrivalTimeDiff(detector_earthfixed_xyz_metres, earth_center, source_right_ascension_radians, source_declination_radians, gpstime);
}


/**
 * Compute the light travel time between two detectors and returns the answer in \c INT8 nanoseconds.
 */
INT8
XLALLightTravelTime(
	const LALDetector *aDet,
	const LALDetector *bDet
)

{
	double deltaLoc[3];

	deltaLoc[0] = aDet->location[0] - bDet->location[0];
	deltaLoc[1] = aDet->location[1] - bDet->location[1];
	deltaLoc[2] = aDet->location[2] - bDet->location[2];

	return (INT8) (1e9 * sqrt(dotprod(deltaLoc, deltaLoc)) / LAL_C_SI);
}
