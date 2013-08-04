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

#ifndef _LALSIMULATION_H
#define _LALSIMULATION_H

#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

const LALDetector *XLALDetectorPrefixToLALDetector(
	const char *string
);
/* FIXME:  compatibility wrapper.  remove when not needed */
const LALDetector *XLALInstrumentNameToLALDetector(const char *string);


REAL8TimeSeries *XLALSimDetectorStrainREAL8TimeSeries(
	const REAL8TimeSeries *hplus,
	const REAL8TimeSeries *hcross,
	REAL8 right_ascension,
	REAL8 declination,
	REAL8 psi,
	LALDetector *detector
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

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMULATION_H */
