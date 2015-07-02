/*
 * Copyright (C) 2008 J. Creighton, K. Cannon
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


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


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
	REAL8 polarization,
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



REAL8 XLALMeasureHPeak(const REAL8TimeSeries *);
REAL8 XLALMeasureIntS1S2DT(const REAL8TimeSeries *, const REAL8TimeSeries *);
REAL8 XLALMeasureHrss(const REAL8TimeSeries *, const REAL8TimeSeries *);
REAL8 XLALMeasureIntHDotSquaredDT(const COMPLEX16FrequencySeries *);
REAL8 XLALMeasureEoverRsquared(REAL8TimeSeries *, REAL8TimeSeries *);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif
