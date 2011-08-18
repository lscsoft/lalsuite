/*
 * Copyright (C) 2011 N. Fotopoulos <nickolas.fotopoulos@ligo.org>
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

#include <lal/LALDatatypes.h>
NRCSID(LALSIMIMRH, "$Id$");

/**
 * Driver routine to compute the non-spinning, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomA in the frequency domain.
 *
 * Reference:
 *   - Waveform: Eq.(4.13) and (4.16) of http://arxiv.org/pdf/0710.2335
 *   - Coefficients: Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and
 *                   Table I of http://arxiv.org/pdf/0712.0343
 *
 * All input parameters should be SI units.
 */
int XLALSimIMRPhenomAGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    LIGOTimeGPS *tRef,                 /**< time at fRef */
    REAL8 phiRef,                      /**< phase at fRef */
    REAL8 fRef,                        /**< reference frequency */
    REAL8 deltaF,                      /**< sampling interval */
    REAL8 m1,                          /**< mass of companion 1 */
    REAL8 m2,                          /**< mass of companion 2 */
    REAL8 f_min,                       /**< start frequency */
    REAL8 f_max,                       /**< end frequency */
    REAL8 distance                     /**< distance of source */
);

/**
 * Driver routine to compute the non-spinning, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomA in the time domain.
 *
 * Reference:
 *   - Waveform: Eq.(4.13) and (4.16) of http://arxiv.org/pdf/0710.2335
 *   - Coefficients: Eq.(4.18) of http://arxiv.org/pdf/0710.2335 and
 *                   Table I of http://arxiv.org/pdf/0712.0343
 *
 * All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomAGenerateTD(
    REAL8TimeSeries **hplus,  /**< +-polarization waveform */
    REAL8TimeSeries **hcross, /**< x-polarization waveform */
    LIGOTimeGPS *tRef,        /**< time at fRef */
    REAL8 phiRef,             /**< phase at fRef */
    REAL8 fRef,               /**< reference frequency */
    REAL8 deltaT,             /**< sampling interval */
    REAL8 m1,                 /**< mass of companion 1 */
    REAL8 m2,                 /**< mass of companion 2 */
    REAL8 f_min,              /**< start frequency */
    REAL8 f_max,              /**< end frequency */
    REAL8 distance,           /**< distance of source */
    REAL8 inclination         /**< inclination of source */
);

/**
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomB in the frequency domain.
 *
 * Reference: http://arxiv.org/pdf/0909.2867
 *   - Waveform: Eq.(1)
 *   - Coefficients: Eq.(2) and Table I
 *
 * All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomBGenerateFD(
    COMPLEX16FrequencySeries **htilde, /**< FD waveform */
    LIGOTimeGPS *tRef,                 /**< time at fRef */
    REAL8 phiRef,                      /**< phase at fRef */
    REAL8 fRef,                        /**< reference frequency */
    REAL8 deltaF,                      /**< sampling interval */
    REAL8 m1,                          /**< mass of companion 1 */
    REAL8 m2,                          /**< mass of companion 2 */
    REAL8 chi,                         /**< mass-weighted aligned-spin parameter */
    REAL8 f_min,                       /**< start frequency */
    REAL8 f_max,                       /**< end frequency */
    REAL8 distance                     /**< distance of source */
);

/**
 * Driver routine to compute the spin-aligned, inspiral-merger-ringdown
 * phenomenological waveform IMRPhenomB in the time domain.
 *
 * Reference: http://arxiv.org/pdf/0909.2867
 *   - Waveform: Eq.(1)
 *   - Coefficients: Eq.(2) and Table I
 *
 * All input parameters should be in SI units. Angles should be in radians.
 */
int XLALSimIMRPhenomBGenerateTD(
    REAL8TimeSeries **hplus,  /**< +-polarization waveform */
    REAL8TimeSeries **hcross, /**< x-polarization waveform */
    LIGOTimeGPS *tRef,        /**< time at fRef */
    REAL8 phiRef,             /**< phase at fRef */
    REAL8 fRef,               /**< reference frequency */
    REAL8 deltaT,             /**< sampling interval */
    REAL8 m1,                 /**< mass of companion 1 */
    REAL8 m2,                 /**< mass of companion 2 */
    REAL8 chi,                /**< mass-weighted aligned-spin parameter */
    REAL8 f_min,              /**< start frequency */
    REAL8 f_max,              /**< end frequency */
    REAL8 distance,           /**< distance of source */
    REAL8 inclination         /**< inclination of source */
);

