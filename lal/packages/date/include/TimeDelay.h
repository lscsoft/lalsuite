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


/**
 * \defgroup TimeDelay_h TimeDelay
 * \ingroup date
 * \brief Provides routines to compute time delay between two detectors
 *
 * \heading{Synopsis}
 * \code
 * #include <lal/TimeDelay.h>
 * \endcode
 */

#ifndef _TIMEDELAY_H
#define _TIMEDELAY_H

#include <lal/Date.h>
#include <lal/DetectorSite.h>

#ifdef __cplusplus
extern "C"
{
#endif

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

#ifdef __cplusplus
}
#endif

#endif /* !defined _TIMEDELAY_H */
