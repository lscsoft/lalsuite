/*
*  Copyright (C) 2007 Alexander Dietz, Drew Keppel, Duncan Brown, David Chin, Jolien Creighton, Kipp Cannon, Stephen Fairhurst
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

#ifndef _TIMEDELAY_H
#define _TIMEDELAY_H

#include <lal/Date.h>
#include <lal/DetectorSite.h>

#ifdef __cplusplus
extern "C"
{
#endif


/**
 * \addtogroup TimeDelay_h
 * \author David Chin, Kipp Cannon
 *
 * \brief Provides routines to compute time delay between two detectors
 *
 * ### Synopsis ###
 *
 * \code
 * #include <lal/TimeDelay.h>
 * \endcode
 *
 * ### Description ###
 *
 * The function XLALTimeDelayFromEarthCenter() Computes difference in arrival
 * time of the same signal at detector and at center of Earth-fixed frame.
 *
 * The function XLALLightTravelTime() computes the light travel time between two detectors and returns the answer in \c INT8 nanoseconds.
 *
 * The function XLALPopulateAccuracyParams() creates an instance of ::InspiralAccuracyList populated with
 * the light-travel times between the detectors, using just the previous function.
 * The function XLALPopulateAccuracyParamsExt(), however, creates an instance of ::InspiralAccuracyList
 * populated with the \c real travel time of a putative signal for the given time and the given sky
 * location (in right ascension and declination, both given in degrees).
 *
 * ### Algorithm ###
 *
 * TBA. See Anderson, <em>et al.</em> \cite ABCCRW_2001 in the mean time.
 *
 * Note that GPS time is passed with both the detectors.  The GPS time of the
 * second detector is \e ignored, and the GPS time for the first detector
 * is taken to be the time when the signal arrives at the center of the
 * Earth.  In practice, this time will be the time of detection of a signal at
 * the first detector, but, as in Anderson, <em>et al.</em>, we make this
 * approximation as it makes little difference.  This time is used to compute
 * a GMST which gives us the orientation of the Earth.
 *
 */
/*@{*/


/*
 * Function prototypes
 */

double
XLALArrivalTimeDiff(
	const double detector1_earthfixed_xyz_metres[3],
	const double detector2_earthfixed_xyz_metres[3],
	const double source_right_ascension_radians,
	const double source_declination_radians,
	const LIGOTimeGPS *gpstime
);


INT8
XLALLightTravelTime ( const LALDetector *aDet,
                      const LALDetector *bDet
                     );


REAL8
XLALTimeDelayFromEarthCenter(
	const double detector_earthfixed_xyz_metres[3],
	double source_right_ascension_radians,
	double source_declination_radians,
	const LIGOTimeGPS *gpstime
);

/*@}*/

#ifdef __cplusplus
}
#endif

#endif /* !defined _TIMEDELAY_H */
