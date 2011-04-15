/*
 * Copyright (C) 2011 Drew Keppel, J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria, Stas Babak, David Churches, B.S. Sathyaprakash, Craig Robinson , Thomas Cokelaer
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

#include <math.h>

#include <lal/LALSimInspiral.h>
#include <lal/LALSimInspiralTaylorT3.h>
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/FindRoot.h>
#include <lal/LALComplex.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALInspiral.h>

#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

NRCSID(LALSIMINSPIRALTAYLORT3C, "$Id$");

/**
 * Estimates the time it takes to go from f_min to f_isco
 * for a given mass pair at a given post-Newtonian order.
 *
 * Estimate is computed by summing the chirp times associated
 * with different post-Newtonian orders.
 */
REAL8 XLALInspiralPNCoalescenceTime(
		REAL8 m1,    /**< mass of companion 1 */
		REAL8 m2,    /**< mass of companion 2 */
		REAL8 f_min, /**< start frequency */
		int O        /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALInspiralPNCoalescenceTime";
	REAL8 tc = 0;
	REAL8 piFl = LAL_PI * f_min;
	REAL8 m = m1 + m2;
	REAL8 nu = m1 * m2 / m / m;
	m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */
	REAL8 nu2 = nu * nu;
	REAL8 nu3 = nu2 * nu;
	REAL8 v = pow(piFl * m, oneby3);
	REAL8 tN = 5.L/256.L / nu * m / pow(v,8.L);

	switch (O) {
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_REAL8(func, XLAL_EINVAL);
		case -1: /* use highest available pN order */
		case 7:
			tc += (-15419335.L / 127008.L - 75703.L / 756.L * nu + 14809.L / 378.L * nu2)
				* LAL_PI * tN * pow(v, 7);
		case 6:
			tc += tN * ((-10052469856691. / 23471078400. + 128. / 3. * LAL_PI * LAL_PI
					+ (15335597827.L / 15240960.L - 451.L / 12.L * LAL_PI * LAL_PI
						+ 352. / 3. * LALINSPIRAL_PNTHETA
						- 2464.L / 9.L * LALINRPIAL_PNLAMBDA) * nu
					+ 6848.L / 105.L * LAL_GAMMA - 15211.L / 1728.L * nu2
					+ 25565.L / 1296.L * nu3)
				+ 6848.L / 105.L * log(4. * v)) * pow(v, 6);
		case 5:
			tc += -5. * (7729. / 252. - 13. / 3. * nu) / (256. * nu * f_min);
		case 4:
			tc += (5.0 / (128.0 * nu * pow(m, oneby3) * pow(piFl, fourby3)))
				* (3058673. / 1016064. + 5429. * nu / 1008. + 617. * nu2 / 144.);
		case 3:
			tc += LAL_PI / (8.0 * nu * pow(m, twoby3) * pow(piFl, fiveby3));
		case 2:
			tc += (3715.0L + (4620.0L * nu)) / (64512.0 * nu * m * pow(piFl, 2.0));
		case 1:
		case 0:
			tc += 5.0L / (256.0L * nu * pow(m, fiveby3) * pow(piFl, eightby3));
	}

	return tc;
}


struct SimInspiralTaylorT3FrequencyParams {
	REAL8 f0;
	REAL8 m1;
	REAL8 m2;
	int O;
};


static REAL8 XLALSimInspiralTaylorT3FrequencyWrapper(REAL8 td, void *pars)
{
	struct SimInspiralTaylorT3FrequencyParams *in = pars;

	return XLALSimInspiralTaylorT3Frequency(td, in->m1, in->m2, in->O) - in->f0;
}


/**
 * Computes the post-Newtonian phase for a given time from coalesence.
 *
 * td = t * nu / (5. * m)
 * where t is the time before coalesnce in units of seconds,
 * nu is the symmetric mass ratio,
 * and m is the total mass in units of kilograms
 */
REAL8 XLALSimInspiralTaylorT3Phasing(
		REAL8 td, /**< post-Newtonian parameter */
	       	REAL8 m1, /**< mass of companion 1 */
	       	REAL8 m2, /**< mass of companion 2 */
	       	int O     /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralTaylorT3Phasing";
	REAL8 m = m1 + m2;
	REAL8 nu = m1 * m2 / m / m;
	m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */
	REAL8 nu2 = nu * nu;
	REAL8 nu3 = nu2 * nu;
	REAL8 thetas[8];
	thetas[0] = 1.;
	thetas[1] = pow(td,-0.125);
	thetas[2] = thetas[1] * thetas[1];
	thetas[5] = thetas[2] * thetas[2] * thetas[1];
	int i;
	REAL8 ans = 0;

	for (i=1; i < O; i++) {
		thetas[i + 1] = thetas[i] * thetas[1];
	}

	switch (O) {
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_REAL8(func, XLAL_EINVAL);
		case -1: /* use highest available pN order */
		case 7:
			ans += thetas[7] * (1.88516689 / 1.73408256 + 488825. / 516096. * nu
				- 141769. / 516096. * nu2) * LAL_PI;
		case 6:
			ans += thetas[6] * ((83.1032450749357 / 5.7682522275840 - 53. / 40. * LAL_PI * LAL_PI
					- 107. / 56. * LAL_GAMMA + (-123.292747421 / 4.161798144
						+ 2.255 / 2.048 * LAL_PI * LAL_PI
						+ 385. / 48. * LALINRPIAL_PNLAMBDA
						- 55. / 16. * LALINSPIRAL_PNTHETA) * nu
					+ 1.54565 / 18.35008 * nu2 - 1.179625 / 1.769472 * nu3)
				+ log(td / 256.) * (107. / 448.));
		case 5:
			/* FIXME: we need to get tn */
			ans += thetas[5] * log(td / 2.977257) * (-(3.8645 / 17.2032 - 65. / 2048. * nu) * LAL_PI);
		case 4:
			ans += thetas[4] * (9.275495 / 14.450688 + 2.84875 / 2.58048 * nu
				+ 1855. / 2048. * nu2);
		case 3:
			ans += thetas[3] * (-0.75) * LAL_PI;
		case 2:
			ans += thetas[2] * (3715. / 8064. + 55. / 96. * nu);
		case 1:
		case 0:
			ans += 1.0;
	}
	ans *= -2. / nu / thetas[5];
	return ans;
}


/**
 * Computes the post-Newtonian frequency for a given time from coalesence.
 *
 * td = t * nu / (5. * m)
 * where t is the time before coalesnce in units of seconds,
 * nu is the symmetric mass ratio,
 * and m is the total mass in units of kilograms
 */
REAL8 XLALSimInspiralTaylorT3Frequency(
		REAL8 td, /**< post-Newtonian parameter */
	       	REAL8 m1, /**< mass of companion 1 */
	       	REAL8 m2, /**< mass of companion 2 */
	       	int O     /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralTaylorT3Frequency";
	REAL8 m = m1 + m2;
	REAL8 nu = m1 * m2 / m / m;
	m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */
	REAL8 nu2 = nu * nu;
	REAL8 nu3 = nu2 * nu;
	REAL8 thetas[8];
	thetas[0] = 1.;
	thetas[1] = pow(td, -0.125);
	thetas[3] = thetas[1] * thetas[1] * thetas[1];
	int i;
	REAL8 ans = 0;

	for (i=1; i < O; i++) {
		thetas[i + 1] = thetas[i] * thetas[1];
	}

	switch (O) {
		default: /* unsupported pN order */
			XLALPrintError("XLAL Error - %s: PN order %d%s not supported\n", func, O/2, O%2?".5":"" );
			XLAL_ERROR_REAL8(func, XLAL_EINVAL);
		case -1: /* use highest available pN order */
		case 7:
			ans += thetas[7] * (-1.88516689 / 4.33520640 - 97765. / 258048. * nu
				+ 141769. / 1290240. * nu2) * LAL_PI;
		case 6:
			ans += thetas[6] * ((-7.20817631400877 / 2.88412611379200
					+ (53. / 200.) * LAL_PI * LAL_PI + 1.07 / 2.80 * LAL_GAMMA
					+ 1.07 / 2.80 * log(2.) + (1.23292747421 / .20808990720
						- 4.51 / 20.48 * LAL_PI * LAL_PI
						- 77. / 48. * LALINRPIAL_PNLAMBDA
						+ 11. / 16. * LALINSPIRAL_PNTHETA) * nu
					- 3.0913 / 183.5008 * nu2 + 2.35925 / 17.69472 * nu3)
				+ log(td) * (-107. / 2240.));
		case 5:
			ans += thetas[5] * (-(7.729 / 21.504 - 13. / 256. * nu) * LAL_PI);
		case 4:
			ans += thetas[4] * (1.855099 / 14.450688 +  5.6975 / 25.8048 * nu
				+ 3.71 / 20.48 * nu2);
		case 3:
			ans += thetas[3] * (-0.3) * LAL_PI;
		case 2:
			ans += thetas[2] * (743. / 2688. + 11. / 32. * nu);
		case 1:
		case 0:
			ans += 1.0;
	}
	ans *= thetas[3] / (8. * LAL_PI * m);
	return ans;
}


/**
 * Computes a post-Newtonian orbit using the Taylor T3 method.
 */
int XLALSimInspiralTaylorT3PNEvolveOrbit(
		REAL8TimeSeries **x,   /**< post-Newtonian parameter [returned] */
	       	REAL8TimeSeries **phi, /**< orbital phase [returned] */
	       	LIGOTimeGPS *tc,       /**< coalescence time */
	       	REAL8 phi0,            /**< initial phase */
	       	REAL8 deltaT,          /**< sampling interval */
		REAL8 m1,              /**< mass of companion 1 */
		REAL8 m2,              /**< mass of companion 2 */
		REAL8 f_min,           /**< start frequency */
		int O                  /**< twice post-Newtonian order */
		)
{
	static const char *func = "XLALSimInspiralPNEvolveOrbitTaylorT3";
	const UINT4 blocklen = 1024;
	const REAL8 xisco = 1.0 / 6.0;
	REAL8 m = m1 + m2;
	REAL8 nu = m1 * m2 / m / m;
	m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */
	REAL8 tC, c1, xmin, xmax, xacc, v, v2, phase, fOld, t, td;
	UINT4 j;
	REAL8 f;
	struct SimInspiralTaylorT3FrequencyParams params;
	void *pars;

	/* allocate memory */

	*x = XLALCreateREAL8TimeSeries("ORBITAL_FREQUENCY_PARAMETER", tc, 0.0, deltaT, &lalDimensionlessUnit,
		blocklen);
	*phi = XLALCreateREAL8TimeSeries("ORBITAL_PHASE", tc, 0.0, deltaT, &lalDimensionlessUnit, blocklen);
	if (!x || !phi)
		XLAL_ERROR(func, XLAL_EFUNC);

	/* compute factor tranforming time to coalecense to pN parameter */

	c1 = nu / (5. * m);

	/* estimate lenth of waveform */

	tC = XLALInspiralPNCoalescenceTime(m1, m2, f_min, O);
	xmin = c1 * tC / 2.;
	xmax = c1 * tC * 3. + 5.;
	xacc = 1.e-6;

	/* fill in parameters for XLALDBracketRoot and XLALDBisectionFindRoot */

	params.m1 = m1;
	params.m2 = m2;
	params.f0 = f_min;
	params.O = O;
	pars = (void*) &params;

	/* bracket root and find root */

	XLALDBracketRoot(XLALSimInspiralTaylorT3FrequencyWrapper, &xmin, &xmax, pars);

	td = XLALDBisectionFindRoot(XLALSimInspiralTaylorT3FrequencyWrapper, xmin, xmax, xacc, pars);
	tC = td / c1;

	/* start waveform generation */

	phase = XLALSimInspiralTaylorT3Phasing(td, m1, m2, O);
	f = f_min;
	v = pow(f * LAL_PI * m, oneby3);
	v2 = v * v;
	(*x)->data->data[0] = v2;
	(*phi)->data->data[0] = phase / 2.;

	j = 0;
	while (1) {

		/* make one step */

		j++;
		fOld = f;
		t = j * deltaT;
		td = c1 * (tC - t);
		phase = XLALSimInspiralTaylorT3Phasing(td, m1, m2, O);
		f = XLALSimInspiralTaylorT3Frequency(td, m1, m2, O);
		v = pow(f * LAL_PI * m, oneby3);
		v2 = v * v;

		/* check termination conditions */

		if (t >= tC) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at coalesence time\n", func);
			break;
		}
		if (v2 >= xisco) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at ISCO\n", func);
			break;
		}
		if (f <= fOld) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated when frequency stalled\n", func);
			break;
		}
	
		/* save current values in vectors but first make sure we don't write past end of vectors */

		if ( j >= (*x)->data->length ) {
			if ( ! XLALResizeREAL8TimeSeries(*x, 0, (*x)->data->length + blocklen) )
				XLAL_ERROR(func, XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*phi, 0, (*phi)->data->length + blocklen) )
				XLAL_ERROR(func, XLAL_EFUNC);
		}
		(*x)->data->data[j] = v2;
		(*phi)->data->data[j] = phase / 2.;
	}

	/* make the correct length */

	if ( ! XLALResizeREAL8TimeSeries(*x, 0, j) )
		XLAL_ERROR(func, XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*phi, 0, j) )
		XLAL_ERROR(func, XLAL_EFUNC);

	/* adjust to correct tc and phic */

	XLALGPSAdd(&(*phi)->epoch, -1.0*j*deltaT);
	XLALGPSAdd(&(*x)->epoch, -1.0*j*deltaT);

	/* adjust so initial phase is phi0 */

	phi0 -= (*phi)->data->data[0];
	for (j = 0; j < (*phi)->data->length; ++j)
		(*phi)->data->data[j] += phi0;

	return (int)(*x)->data->length;
}


/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine allows the user to specify different pN orders
 * for phasing calcuation vs. amplitude calculations.
 */
int XLALSimInspiralTaylorT3PNGenerator(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 x0,                 /**< tail-term gauge choice thing (if you don't know, just set it to zero) */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int amplitudeO,           /**< twice post-Newtonian amplitude order */
	       	int phaseO                /**< twice post-Newtonian phase order */
		)
{
	static const char *func = "XLALSimInspiralPNGenerator";
	REAL8TimeSeries *x;
	REAL8TimeSeries *phi;
	int status;
	int n;
	n = XLALSimInspiralTaylorT3PNEvolveOrbit(&x, &phi, tc, phic, deltaT, m1, m2, f_min, phaseO);
	if ( n < 0 )
		XLAL_ERROR(func, XLAL_EFUNC);
	status = XLALSimInspiralPNPolarizationWaveforms(hplus, hcross, x, phi, x0, m1, m2, r, i, amplitudeO);
	XLALDestroyREAL8TimeSeries(phi);
	XLALDestroyREAL8TimeSeries(x);
	if ( status < 0 )
		XLAL_ERROR(func, XLAL_EFUNC);
	return n;
}


/**
 * Driver routine to compute the post-Newtonian inspiral waveform.
 *
 * This routine uses the same pN order for phasing and amplitude
 * (unless the order is -1 in which case the highest available
 * order is used for both of these -- which might not be the same).
 *
 * Log terms in amplitudes are ignored.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT3PN(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian order */
		)
{
	/* set x0=0 to ignore log terms */
	return XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, f_min, r, i, O, O);
}


/**
 * Driver routine to compute the restricted post-Newtonian inspiral waveform.
 *
 * This routine computes the phasing to the specified order, but
 * only computes the amplitudes to the Newtonian (quadrupole) order.
 *
 * Log terms in amplitudes are ignored.  This is a gauge choice.
 */
int XLALSimInspiralTaylorT3PNRestricted(
		REAL8TimeSeries **hplus,  /**< +-polarization waveform */
	       	REAL8TimeSeries **hcross, /**< x-polarization waveform */
	       	LIGOTimeGPS *tc,          /**< coalescence time */
	       	REAL8 phic,               /**< coalescence phase */
	       	REAL8 deltaT,             /**< sampling interval */
	       	REAL8 m1,                 /**< mass of companion 1 */
	       	REAL8 m2,                 /**< mass of companion 2 */
	       	REAL8 f_min,              /**< start frequency */
	       	REAL8 r,                  /**< distance of source */
	       	REAL8 i,                  /**< inclination of source (rad) */
	       	int O                     /**< twice post-Newtonian phase order */
		)
{
	/* use Newtonian order for amplitude */
	/* set x0=0 to ignore log terms */
	return XLALSimInspiralTaylorT3PNGenerator(hplus, hcross, tc, phic, 0.0, deltaT, m1, m2, f_min, r, i, 0, O);
}
