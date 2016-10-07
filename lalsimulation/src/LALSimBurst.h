/*
 * Copyright (C) 2008 J. Creighton
 * Copyright (C) 2008,2015 K. Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */


/**
 * @author Kipp Cannon, Jolien Creighton
 * @addtogroup LALSimBurst_h Header LALSimBurst.h
 * @ingroup lalsimulation_burst
 * @brief Routines to generate burst waveforms.
 * @details
 * These routines generate several burst waveforms used in searches for
 * gravitational waves, including sine-Gaussian waveforms, cosmic string
 * cusp waveforms, and band- and time-limited white-noise burst waveforms.
 * Also included are several general-purpose routines to measure the
 * properties of gravitational wave waveforms like the "hrss" and peak
 * strain.  These are useful for imposing normalizations and other
 * diagnostic activities.
 *
 * \f[
 * \DeclareMathOperator{\order}{O}
 * \newcommand{\Msol}{{M_{\Sol}}}
 * \newcommand{\Sol}{\odot}
 * \newcommand{\aye}{\mathrm{i}}
 * \newcommand{\conj}[1]{#1^{*}}
 * \newcommand{\diff}{\,\mathrm{d}}
 * \newcommand{\ee}{\mathrm{e}}
 * \newcommand{\magnitude}[1]{\left|#1\right|}
 * \newcommand{\mean}[1]{\left\langle#1\right\rangle}
 * \f]
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */

#ifndef _LALSIMBURST_H
#define _LALSIMBURST_H

#include <gsl/gsl_rng.h>
#include <lal/LALDatatypes.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

/*
 * ============================================================================
 *
 *                            Function Prototypes
 *
 * ============================================================================
 */


/** @{ */


int XLALGenerateImpulseBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 hpeak,
	REAL8 delta_t
);


int XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 duration,
	REAL8 frequency,
	REAL8 bandwidth,
	REAL8 eccentricity,
	REAL8 phase,
	REAL8 int_hdot_squared,
	REAL8 delta_t,
	gsl_rng *rng
);


int XLALGenerateStringCusp(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 amplitude,
	REAL8 f_high,
	REAL8 delta_t
);


double XLALSimBurstSineGaussianQ(
	double duration,
	double centre_frequency
);


double XLALSimBurstSineGaussianDuration(
	double Q,
	double centre_frequency
);


int XLALSimBurstSineGaussian(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 Q,
	REAL8 centre_frequency,
	REAL8 hrss,
	REAL8 eccentricity,
	REAL8 phase,
	REAL8 delta_t
);


int XLALSimBurstGaussian(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	REAL8 duration,
	REAL8 hrss,
	REAL8 delta_t
);

int XLALSimBurstImg(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross, 
	REAL8Array *image,
	double dt,
	double df,
	double fstart,
	double hrss,
	double deltaT,
	gsl_rng *rng
);


int XLALSimUnicorn(
	REAL8TimeSeries **hplus,
	REAL8TimeSeries **hcross,
	double f_min,
	double f_max,
	double V,
	double hrss,
	double deltaT,
	gsl_rng *rng
);


COMPLEX16 XLALMeasureHPeak(const REAL8TimeSeries *, const REAL8TimeSeries *, unsigned *);
REAL8 XLALMeasureIntS1S2DT(const REAL8TimeSeries *, const REAL8TimeSeries *);
REAL8 XLALMeasureHrss(const REAL8TimeSeries *, const REAL8TimeSeries *);
REAL8 XLALMeasureIntHDotSquaredDT(const COMPLEX16FrequencySeries *);
REAL8 XLALMeasureEoverRsquared(REAL8TimeSeries *, REAL8TimeSeries *);


/** @} */


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
#endif /*_LALSIMBURST_H */
