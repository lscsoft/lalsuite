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
#define LAL_USE_COMPLEX_SHORT_MACROS
#include <lal/FindRoot.h>
#include <lal/LALComplex.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALSimInspiralTaylorT3.h>
#include <lal/LALSimInspiraldEnergyFlux.h>
#include <lal/LALSimInspiralPhasing3.h>
#include <lal/LALSimInspiralFrequency3.h>

#include "check_series_macros.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

NRCSID(LALSIMINSPIRALTAYLORT3C, "$Id$");


typedef struct
{
	REAL8 (*func)(REAL8 tC, expnCoeffsTaylorT3 *ak);
	expnCoeffsTaylorT3 ak;
}
FreqInFromChirptime;

static REAL8 XLALInspiralFrequency3Wrapper(REAL8 tC, void *pars)
{
  static const char *func = "XLALInspiralFrequency3Wrapper";

  FreqInFromChirptime *in;
  REAL8 freq, f;

  in = (FreqInFromChirptime *) pars;
  freq = in->func(tC, &(in->ak));
  if (XLAL_IS_REAL8_FAIL_NAN(freq))
    XLAL_ERROR_REAL8(func, XLAL_EFUNC);
  f = freq - in->ak.f0;

  /*
  fprintf(stderr, "Here freq=%e f=%e tc=%e f0=%e\n", freq, *f, tC, in->ak.f0);
   */

  return f;
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
	const UINT4 blocklen = 1024;
	const REAL8 xisco = 1.0 / 6.0;
	REAL8 m = m1 + m2;
	REAL8 nu = m1 * m2 / m / m;
	m *= LAL_G_SI / pow(LAL_C_SI, 3.0); /* convert m from kilograms to seconds */
	REAL8 mass1 = m1 / LAL_MSUN_SI; /* convert m1 from kilograms to solar masses */
	REAL8 mass2 = m2 / LAL_MSUN_SI; /* convert m2 from kilograms to solar masses */
	REAL8 tmptC, tC, c1, xmin, xmax, xacc, v, v2, phase, fOld, t, td, temp, tempMin = 0, tempMax = 0;
	REAL8 (*freqfunction)(REAL8, void *);
	UINT4 j;
	REAL8 f;
	void *pars;

	expnFuncTaylorT3 expnfunc;
	expnCoeffsTaylorT3 ak;
	FreqInFromChirptime timeIn;

	/* allocate memory */

	*x = XLALCreateREAL8TimeSeries("ORBITAL_FREQUENCY_PARAMETER", tc, 0.0, deltaT, &lalDimensionlessUnit,
		blocklen);
	*phi = XLALCreateREAL8TimeSeries("ORBITAL_PHASE", tc, 0.0, deltaT, &lalDimensionlessUnit, blocklen);
	if (!x || !phi)
		XLAL_ERROR(__func__, XLAL_EFUNC);


	/* initialize expnCoeffsTaylorT3 and expnFuncTaylorT3 structures */
	if (XLALSimInspiralTaylorT3Setup(&ak, &expnfunc, deltaT, mass1, mass2, f_min, O))
		XLAL_ERROR(__func__, XLAL_EFUNC);

	tC = XLALSimInspiralChirpLength(m1, m2, f_min, O);
	c1 = nu/(5.*m);

	/*
	 * In Jan 2003 we realized that the tC determined as a sum of chirp
	 * times is not quite the tC that should enter the definition of Theta
	 * in the expression for the frequency as a function of time (see DIS3,
	 * 2000). This is because chirp times are obtained by inverting t(f).
	 * Rather tC should be obtained by solving the equation f0 - f(tC) = 0.
	 * This is what is implemented below.
	 */

	timeIn.func = expnfunc.frequency3;
	timeIn.ak = ak;
	freqfunction = &XLALInspiralFrequency3Wrapper;
	xmin = c1*tC/2.;
	xmax = c1*tC*2.;
	xacc = 1.e-6;
	pars = (void*) &timeIn;
	/* tc is the instant of coalescence */

	/* we add 5 so that if tC is small then xmax
	 * is always greater than a given value (here 5)*/
	xmax = c1*tC*3 + 5.;

	/* for x in [xmin, xmax], we search the value which gives the max
	 * frequency.  and keep the corresponding rootIn.xmin. */

	for (tmptC = c1*tC/1000.; tmptC < xmax; tmptC+=c1*tC/1000.){
		temp = XLALInspiralFrequency3Wrapper(tmptC , pars);
		if (XLAL_IS_REAL8_FAIL_NAN(temp))
			XLAL_ERROR(__func__, XLAL_EFUNC);
		if (temp > tempMax) {
			xmin = tmptC;
			tempMax = temp;
		}
		if (temp < tempMin) {
			tempMin = temp;
		}
	}

	/* if we have found a value positive then everything should be fine in
	 * the BissectionFindRoot function */
	if (tempMax > 0  &&  tempMin < 0){
		tC = XLALDBisectionFindRoot (freqfunction, xmin, xmax, xacc, pars);
		if (XLAL_IS_REAL8_FAIL_NAN(tC))
			XLAL_ERROR(__func__, XLAL_EFUNC);
	}
	else{
		XLALPrintError("Can't find good bracket for BisectionFindRoot");
		XLAL_ERROR(__func__, XLAL_EMAXITER);
	}

	tC /= c1;

	/* start waveform generation */

	t = 0.;
	td = c1 * (tC - t);
	phase = expnfunc.phasing3(td, &ak);
	if (XLAL_IS_REAL8_FAIL_NAN(phase))
		XLAL_ERROR(__func__, XLAL_EFUNC);
	f = expnfunc.frequency3(td, &ak);
	if (XLAL_IS_REAL8_FAIL_NAN(f))
		XLAL_ERROR(__func__, XLAL_EFUNC);

	v = cbrt(f * LAL_PI * m);
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
		phase = expnfunc.phasing3(td, &ak);
		if (XLAL_IS_REAL8_FAIL_NAN(phase))
			XLAL_ERROR(__func__, XLAL_EFUNC);
		f = expnfunc.frequency3(td, &ak);
		if (XLAL_IS_REAL8_FAIL_NAN(f))
			XLAL_ERROR(__func__, XLAL_EFUNC);
		v = cbrt(f * LAL_PI * m);
		v2 = v * v;

		/* check termination conditions */

		if (t >= tC) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at coalesence time\n", __func__);
			break;
		}
		if (v2 >= xisco) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated at ISCO\n", __func__);
			break;
		}
		if (f <= fOld) {
			XLALPrintInfo("XLAL Info - %s: PN inspiral terminated when frequency stalled\n", __func__);
			break;
		}
	
		/* save current values in vectors but first make sure we don't write past end of vectors */

		if ( j >= (*x)->data->length ) {
			if ( ! XLALResizeREAL8TimeSeries(*x, 0, (*x)->data->length + blocklen) )
				XLAL_ERROR(__func__, XLAL_EFUNC);
			if ( ! XLALResizeREAL8TimeSeries(*phi, 0, (*phi)->data->length + blocklen) )
				XLAL_ERROR(__func__, XLAL_EFUNC);
		}
		(*x)->data->data[j] = v2;
		(*phi)->data->data[j] = phase / 2.;
	}

	/* make the correct length */

	if ( ! XLALResizeREAL8TimeSeries(*x, 0, j) )
		XLAL_ERROR(__func__, XLAL_EFUNC);
	if ( ! XLALResizeREAL8TimeSeries(*phi, 0, j) )
		XLAL_ERROR(__func__, XLAL_EFUNC);

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
