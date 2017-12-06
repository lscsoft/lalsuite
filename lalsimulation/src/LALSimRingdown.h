/*
 * Copyright (C) 2011 J. Clark
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

#ifndef _LALSIMRINGDOWN_H
#define _LALSIMRINGDOWN_H

#include <lal/LALDatatypes.h>

/*#include <lal/LALRCSID.h>*/

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/*NRCSID(LALSIMRINGDOWNH, "$Id$");*/

int XLALSimRingdownFD(
	COMPLEX16FrequencySeries **htilde,
    double f_min,
    double f_max,
	double deltaF,			/**< sampling interval (s) */
	double frequency,		/* ringdown frequency (Hz) */
	double quality,			/* quality factor = pi*decay*frequency*/
    double hrss             /* Root-sum-squared amplitude */
);

/**
 * Computes the waveform for a generic ringdown *
 */
int XLALSimRingdown(
		REAL8TimeSeries **hplus,	/**< plus-polarization waveform [returned] */
		REAL8TimeSeries **hcross,	/**< cross-polarization waveform [returned] */
		double deltaT,			/**< sampling interval (s) */
		double Amp,			/**< initial intrinsic amplitude of ringdown */
        double omega0,		/**< f-mode oscillation frequency (rad) */
        double quality,			/* quality factor = pi*decay*frequency*/
		double phi0,		/**< initial phase of ringdown (rad) */
		double theta,		/**< inclination of source's spin axis (rad) */
		double azimuth,		/**< azimuthal angle */
		int l,		        /**< polar mode number */
        int m		        /**< azimuthal mode number */
		);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMRINGDOWN_H */
