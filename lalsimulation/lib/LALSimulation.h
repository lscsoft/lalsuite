/*
 * Copyright (C) 2008 J. Creighton, T. Creighton
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
 * @author Kipp Cannon, Jolien Creighton, Teviet Creighton
 * @addtogroup LALSimulation_h	Header LALSimulation.h
 * @ingroup lalsimulation_general
 * @brief Routines to calculate detector strain for general waveforms.
 * @details
 * These routines compute the external strain on a detector caused by
 * a general waveform and inject these waveforms into time series.
 *
 * The following example illustrates the basic usage of these routines.
 *
 * @code
 * #include <lal/Date.h>
 * #include <lal/LALSimulation.h>
 * ...
 * const LALDetector *detector;
 * REAL8TimeSeries *data;
 * REAL8TimeSeries *strain;
 * REAL8TimeSeries *hplus;
 * REAL8TimeSeries *hcross;
 * LIGOTimeGPS geocentric_arrival_time;
 * double right_ascension, declination, psi;
 * ...
 * // get detector data
 * // generate waveform hplus and hcross
 * // set geocentric_arrival_time for the injection
 * // set right_ascension, declination, and psi for the injection
 * ...
 * XLALGPSAddGPS(&hplus->epoch, &geocentric_arrival_time);
 * XLALGPSAddGPS(&hcross->epoch, &geocentric_arrival_time);
 * detector = XLALDetectorPrefixToLALDetector(data->name);
 * strain = XLALSimDetectorStrainREAL8TimeSeries(hplus, hcross, right_ascension, declination, psi, detector);
 * XLALSimAddInjectionREAL8TimeSeries(data, strain, NULL);
 * @endcode
 *
 * ### Coordinate Systems
 *
 * The diagram below illustrates the relationship between the wave frame
 * (X,Y,Z) and the Earth's equatorial frame.  The Earth is at the origin
 * of the diagram with its North pole in the direction indicated by N.
 *
 * @anchor lalsimulation_inject
 * @image html lalsimulation_inject.svg "Injection Coordinates"
 *
 * The gravitational wave travels in the Z-direction in this diagram,
 * and the reference direction on the wave-plane (the X-Y-plane) is
 * given by the X-axis.  Note that the direction to the source is in
 * the negative Z direction.
 *
 * The plus- and cross-polarizations of the gravitational waveform are defined
 * in this wave frame.  Specifically, if \f$ h^{ij} \f$ is the metric
 * perturbation, then
 * \f[ h_+ = \frac12 ( \hat{p}_i \hat{p}_j - \hat{q}_i \hat{q}_j ) h^{ij} \f]
 * and
 * \f[ h_\times = \frac12 ( \hat{p}_i \hat{q}_j + \hat{q}_i \hat{p}_j ) h^{ij} \f]
 * where \f$ \hat{p}_i \f$ are the components of the unit vector pointing
 * along the X-axis and \f$ \hat{q}_i \f$ are the components of the unit
 * vector pointing along the Y-axis.
 *
 * The angles relating the wave frame to the equatorial frame are:
 *
 *  * Declination (&delta;).  The angle along the hour circle passing through
 *    the source between the equatorial plane and the source.  It is positive
 *    for a source that is north of the equatorial plane and negative for
 *    a source that is south of the equatorial plane.  The hour circle passing
 *    through the source is the great arc on the celestial sphere that passes
 *    through the north celestial pole (marked as N in the diagram), the
 *    source, and the south celestial pole.
 *
 *    Note that the diagram depicts a source having a negative declination.
 *
 *  * Right ascension (&alpha;).  The angle from the vernal equinox
 *    @htmlonly &#x2648; @endhtmlonly to the hour circle passing through
 *    the source.  The angle is measured counter-clockwise about the axis
 *    pointing toward the north celestial pole.
 *
 *  * Greenwich sidereal time (GST).  The right ascension of the prime
 *    meridian at the time of arrival of the signal.
 *
 *  * Greenwich hour angle (GHA).  The angle along the equatorial plane between
 *    the prime meridian and the hour circle containing the source.  This
 *    angle is measured @e clockwise about the axis pointing toward the
 *    north celestial pole.
 *
 *    The right ascension, Greenwich hour angle, and Greenwich sidereal time
 *    are related by GHA = GST - &alpha;.
 *
 *  * Polarization angle (&psi;).  The angle from the ascending line of
 *    nodes to the X-axis of the wave plane.  The angle is measured counter
 *    clockwise about the Z-axis of the wave plane.
 *
 * @sa
 * A complete description of the coordinate conventions adopted here can be
 * found in
 * > Warren Anderson, Patrick Brady, David Chin, Jolien Creighton, Keith Riles,
 * > and John Whelan,
 * > "Beam Pattern Response Functions and Times of Arrival for Earthbound
 * > Interferometer",
 * > LIGO Technical Document LIGO-T010110-v1 (2009)
 * > https://dcc.ligo.org/LIGO-T010110/public
 *
 * @sa
 * The conventions are also described in Appendix B of
 * > Warren G. Anderson, Patrick R. Brady, Jolien D. E. Creighton, and Eanna E.
 * > Flanagan,
 * > "Excess power statistic for detection of burst sources of gravitational
 * > radiation",
 * > Phys. Rev. D @b 63, 042003 (2001)
 * > http://dx.doi.org/10.1103/PhysRevD.63.042003
 * > http://arxiv.org/abs/gr-qc/0008066
 */

#ifndef _LALSIMULATION_H
#define _LALSIMULATION_H

#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/** @{ */


const LALDetector *XLALDetectorPrefixToLALDetector(const char *string);

REAL8TimeSeries *XLALSimDetectorStrainREAL8TimeSeries(
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross,
	REAL8 right_ascension,
	REAL8 declination,
	REAL8 psi,
	const LALDetector *detector
);

int XLALSimAddInjectionREAL8TimeSeries(
	REAL8TimeSeries *target,
	REAL8TimeSeries *h,
	const COMPLEX16FrequencySeries *response
);

int XLALSimAddInjectionREAL4TimeSeries(
	REAL4TimeSeries *target,
	REAL4TimeSeries *h,
	const COMPLEX8FrequencySeries *response
);

int XLALSimInjectDetectorStrainREAL8TimeSeries(
	REAL8TimeSeries *target,
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX16FrequencySeries *response
);

int XLALSimInjectDetectorStrainREAL4TimeSeries(
	REAL4TimeSeries *target,
	const REAL4TimeSeries *hplus,
	const REAL4TimeSeries *hcross,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX8FrequencySeries *response
);

int XLALSimInjectLWLDetectorStrainREAL8TimeSeries(
	REAL8TimeSeries *target,
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX16FrequencySeries *response
);

int XLALSimInjectLWLDetectorStrainREAL4TimeSeries(
	REAL4TimeSeries *target,
	const REAL4TimeSeries *hplus,
	const REAL4TimeSeries *hcross,
	double ra,
	double dec,
	double psi,
	LALDetector *detector,
	const COMPLEX8FrequencySeries *response
);

/** @} */

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMULATION_H */
