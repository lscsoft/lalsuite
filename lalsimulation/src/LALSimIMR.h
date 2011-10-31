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

#ifndef _LALSIMIMR_H
#define _LALSIMIMR_H

#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

NRCSID(LALSIMIMRH, "$Id$");

/**
 * The number of e-folds of ringdown which should be attached for 
 * EOBNR models
 */
#define EOB_RD_EFOLDS 10.0

/**
 * Constant which comes up in some of the EOB models. Its value is
 * (94/3 -41/32*pi*pi)
 */
#define ninty4by3etc 18.687902694437592603

/**
 * Enumerator for choosing the reference frame associated with
 * PSpinInspiralRD waveforms.
 */
typedef enum {
  TotalJ,
  View,
  OrbitalL,
} InputAxis;

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
    REAL8 phi0,                        /**< initial phase */
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
    LIGOTimeGPS *tPeak,       /**< time at peak amplitude */
    REAL8 phiPeak,            /**< phase at peak */
    REAL8 deltaT,             /**< sampling interval */
    REAL8 m1,                 /**< mass of companion 1 */
    REAL8 m2,                 /**< mass of companion 2 */
    REAL8 f_min,              /**< start frequency */
    REAL8 f_max,              /**< end frequency */
    REAL8 distance,           /**< distance of source */
    REAL8 inclination         /**< inclination of source */
);

/**
 * Compute the dimensionless, spin-aligned parameter chi as used in the
 * IMRPhenomB waveform. This is different from chi in SpinTaylorRedSpin!
 * Reference: http://arxiv.org/pdf/0909.2867, paragraph 3.
 */
double XLALSimIMRPhenomBComputeChi(
    const REAL8 m1,                          /**< mass of companion 1 */
    const REAL8 m2,                          /**< mass of companion 2 */
    const REAL8 s1z,                         /**< spin of companion 1 */
    const REAL8 s2z                          /**< spin of companion 2 */
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
    REAL8 phi0,                        /**< initial phase */
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
    LIGOTimeGPS *tPeak,       /**< time at peak amplitude */
    REAL8 phiPeak,            /**< phase at peak */
    REAL8 deltaT,             /**< sampling interval */
    REAL8 m1,                 /**< mass of companion 1 */
    REAL8 m2,                 /**< mass of companion 2 */
    REAL8 chi,                /**< mass-weighted aligned-spin parameter */
    REAL8 f_min,              /**< start frequency */
    REAL8 f_max,              /**< end frequency */
    REAL8 distance,           /**< distance of source */
    REAL8 inclination         /**< inclination of source */
);

/**
 * This function generates the plus and cross polarizations for the dominant
 * (2,2) mode of the EOBNRv2 approximant. This model is defined in Pan et al,
 * arXiv:1106.1021v1 [gr-qc].
 */
int XLALSimIMREOBNRv2DominantMode(
    REAL8TimeSeries **hplus,      /**<< The +-polarization waveform (returned) */
    REAL8TimeSeries **hcross,     /**<< The x-polarization waveform (returned) */
    LIGOTimeGPS      *tC,         /**<< The coalescence time (defined as above) */
    const REAL8       phiC,       /**<< The phase at the coalescence time */
    const REAL8       deltaT,     /**<< Sampling interval (in seconds) */
    const REAL8       m1SI,       /**<< First component mass (in kg) */
    const REAL8       m2SI,       /**<< Second component mass (in kg) */
    const REAL8       fLower,     /**<< Starting frequency (in Hz) */
    const REAL8       distance,   /**<< Distance to source (in metres) */
    const REAL8       inclination /**<< Inclination of the source (in radians) */
);

/**
 * This function generates the plus and cross polarizations for the EOBNRv2 approximant
 * with all available modes included. This model is defined in Pan et al,
 * arXiv:1106.1021v1 [gr-qc].
 */
int XLALSimIMREOBNRv2AllModes(
    REAL8TimeSeries **hplus,      /**<< The +-polarization waveform (returned) */
    REAL8TimeSeries **hcross,     /**<< The x-polarization waveform (returned) */
    LIGOTimeGPS      *tC,         /**<< The coalescence time (defined as above) */
    const REAL8       phiC,       /**<< The phase at the coalescence time */
    const REAL8       deltaT,     /**<< Sampling interval (in seconds) */
    const REAL8       m1SI,       /**<< First component mass (in kg) */
    const REAL8       m2SI,       /**<< Second component mass (in kg) */
    const REAL8       fLower,     /**<< Starting frequency (in Hz) */
    const REAL8       distance,   /**<< Distance to source (in metres) */
    const REAL8       inclination /**<< Inclination of the source (in radians) */
);

/**
 * Driver routine to compute a precessing post-Newtonian inspiral-merger-ringdown waveform
 */
int XLALSimInspiralPSpinInspiralRDGenerator(
    REAL8TimeSeries **hplus,    /**< +-polarization waveform */
    REAL8TimeSeries **hcross,   /**< x-polarization waveform */
    LIGOTimeGPS *t0,            /**< start time */
    REAL8 phi0,                 /**< start phase */
    REAL8 deltaT,               /**< sampling interval */
    REAL8 m1,                   /**< mass of companion 1 */
    REAL8 m2,                   /**< mass of companion 2 */
    REAL8 f_min,                /**< start frequency */
    REAL8 r,                    /**< distance of source */
    REAL8 iota,                 /**< inclination of source (rad) */
    REAL8 spin1[3],             /**< Spin vector on mass1 */
    REAL8 spin2[3],             /**< Spin vector on mass2 */
    int phaseO,                 /**< twice post-Newtonian phase order */
    InputAxis axisChoice       	/**< Choice of axis for input spin params */
    );

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMIMR_H */
