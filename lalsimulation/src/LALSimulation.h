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
 * LALDetector *detector;
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
