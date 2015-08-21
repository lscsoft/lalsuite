/*
 * Copyright (C) 2015 J. Creighton
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
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/LALSimUtils.h>

#include "check_series_macros.h"

/**
 * Creates a new interpolated psd from an existing psd (if necessary).
 *
 * The returned frequency series will have the specified values of @p f0, @p
 * deltaF, and @p length.  If the original frequency series @p old had the same
 * values, this routine simply returns a pointer to the @p old.  Otherwise,
 * this routine allocates a new frequency series.
 *
 * The interpolation in performed is linear-in-log (since a psd is
 * positive-definite).  Any component of the psd that is 0 is taken to be
 * invalid, and any component of the returned psd that would depend on a 0
 * component of the original psd is also set to 0.  Likewise, any components of
 * the new psd that would require extrapolation of the original psd are set to
 * zero.  The routines defined here adhere to the convention that frequency
 * components in which the psd is 0 are ignored.
 *
 * @note The calling routine must test if the pointer to the psd returned by
 * this routine is the same as the pointer to the psd passed to this routine to
 * determine if the returned psd needs to be freed separately.
 */
static REAL8FrequencySeries * create_interpolated_psd(double f0, double deltaF, size_t length, REAL8FrequencySeries *old)
{
	REAL8FrequencySeries *new;
	size_t k, kmin, kmax;

	/* see if we can get away without interpolating */
	if (fabs(deltaF - old->deltaF) < LAL_REAL8_EPS) {
		int first;

		/* do we need to do anything? */
		if (fabs(f0 - old->f0) < LAL_REAL8_EPS && length == old->data->length)
			return old;

		/* is this just an integer shift / resize? */
		first = round(f0 - old->f0) / deltaF;
		if (fabs(first * deltaF - f0 + old->f0) < LAL_REAL8_EPS) {
			/* is new a subset of old? */
			if (first >=0 && old->data->length >= length + first)
				new = XLALCutREAL8FrequencySeries(old, first, length);
			else {
				/* copy and resize */
				new = XLALCutREAL8FrequencySeries(old, 0, old->data->length);
				new = XLALResizeREAL8FrequencySeries(new, first, length);
			}
			return new;
		}
	}

	/* we have to actually interpolate... */
	new = XLALCreateREAL8FrequencySeries(old->name, &old->epoch, f0, deltaF, &old->sampleUnits, length);

	/* determine the limits of where interpolation can occur */
	if (f0 < old->f0)
		kmin = ceil((old->f0 - f0) / deltaF);
	else
		kmin = 0;

	if (f0 + length * deltaF > old->f0 + old->data->length * old->deltaF)
		kmax = floor((old->f0 + old->data->length * old->deltaF - f0) / deltaF);
	else
		kmax = length;

	/* zero out any invalid regions */
	for (k = 0; k < kmin; ++k)
		new->data->data[k] = 0.0;
	for (k = kmax; k < length; ++k)
		new->data->data[k] = 0.0;

	for (k = kmin; k < kmax; ++k) {
		double x, ix;
		double y0, y1;
		size_t i;

		x = modf((f0 + k * deltaF - old->f0) / old->deltaF, &ix);
		i = (size_t)(ix);

		/* safety checks are probably not necessary */
		y0 = ix < 0.0 ? 0.0 : old->data->data[i];
		y1 = i >= old->data->length - 1 ? 0.0 : old->data->data[i+1];
		
		if (y0 == 0.0 || y1 == 0.0)
			new->data->data[k] = 0.0; /* invalid data */
		else /* linear-in-log interpolation */
			new->data->data[k] = exp((1.0 - x) * log(y0) + x * log(y1));
	}

	return new;
}

/**
 * @brief Computes the sense-monitor range for a binary neutron star standard
 * siren signal for a given one-sided detector noise power spectral density.
 *
 * The "Standard Siren" is a restricted (0 pN in amplitude) gravitational
 * waveform produced by a binary neutron star system comprised of two neutron
 * stars, each having a mass of 1.4 Msun.  A circular inspiral of point
 * particles is also assumed.  The "Sense-Monitor Range" is corresponds to
 * the range measure reported by the LIGO control room monitor @c SenseMonitor.
 * This range is the radius of a sphere that has contains as many sources as
 * the number that would produce a characteristic signal-to-noise ratio 8 in a
 * detector under the assumption of a homogeneous distribution of sources
 * having random orientations.  No cosmological effects are included in this
 * sense-monitor range measure.  The sense-monitor range \f$\cal R\f$ is
 * related to the horizon distance \f$D_{\rm hor}\f$ by
 * \f$D_{\rm hor} \approx 2.26478\cal R\f$.
 *
 * See XLALMeasureStandardSirenHorizonDistance() for a description of the
 * horizon distance.
 *
 * See XLALMeasureSNRFD() for further discussion about characteristic
 * signal-to-noise ratios.
 *
 * @sa Appendix D of
 * Bruce Allen, Warren G. Anderson, Patrick R. Brady, Duncan A. Brown, and
 * Jolien D. E. Creighton, "FINDCHIRP: An algorithm for detection of
 * gravitational waves from inspiraling compact binaries", Phys. Rev. D @b 85,
 * 122006 (2012) http://dx.doi.org/10.1103/PhysRevD.85.122006
 *
 * @param psd The one-sided detector strain noise power spectral density.
 * @param f_min The lower bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to 0 or a negative value for no
 * lower bound.
 * @param f_max The upper bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to a negative value for
 * no upper bound.
 * @returns The sense-monitor range in meters.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 */
double XLALMeasureStandardSirenSenseMonitorRange(REAL8FrequencySeries *psd, double f_min, double f_max)
{
	double horizon_distance = XLALMeasureStandardSirenHorizonDistance(psd, f_min, f_max);
	return horizon_distance / LAL_HORIZON_DISTANCE_OVER_SENSEMON_RANGE;
}

/**
 * @brief Computes the horizon distance for a binary neutron star standard
 * siren signal for a given one-sided detector noise power spectral density.
 *
 * The "Standard Siren" is a restricted (0 pN in amplitude) gravitational
 * waveform produced by a binary neutron star system comprised of two neutron
 * stars, each having a mass of 1.4 Msun.  A circular inspiral of point
 * particles is also assumed.  The horizon distance is the distance at which
 * such a system that is optimally oriented (face on) and located (at an
 * interferometer's zenith) would induce a characteristic signal-to-noise
 * ratio of 8.  No cosmological effects are included in this horizon distance
 * measure.
 *
 * See XLALMeasureSNRFD() for further discussion about characteristic
 * signal-to-noise ratios.
 *
 * @sa Appendix D of
 * Bruce Allen, Warren G. Anderson, Patrick R. Brady, Duncan A. Brown, and
 * Jolien D. E. Creighton, "FINDCHIRP: An algorithm for detection of
 * gravitational waves from inspiraling compact binaries", Phys. Rev. D @b 85,
 * 122006 (2012) http://dx.doi.org/10.1103/PhysRevD.85.122006
 *
 * @param psd The one-sided detector strain noise power spectral density.
 * @param f_min The lower bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to 0 or a negative value for no
 * lower bound.
 * @param f_max The upper bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to a negative value for
 * no upper bound.
 * @returns The horizon distance in meters.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 */
double XLALMeasureStandardSirenHorizonDistance(REAL8FrequencySeries *psd, double f_min, double f_max)
{
	double Mpc = 1e6 * LAL_PC_SI;
	double snr_8 = 8.0;
	double snr;
	snr = XLALMeasureStandardSirenSNR(psd, f_min, f_max);
	return Mpc * snr / snr_8;
}

/**
 * @brief Computes the characteristic signal-to-noise for a binary neutron star
 * standard siren signal located at an effective distance of 1 Mpc for a given
 * one-sided detector noise power spectral density.
 *
 * The "Standard Siren" is a restricted (0 pN in amplitude) gravitational
 * waveform produced by a binary neutron star system comprised of two neutron
 * stars, each having a mass of 1.4 Msun, that is optimally oriented (face on)
 * and located (at an interferometer's zenith) at a distance of 1 Mpc.  A
 * circular inspiral of point particles is also assumed.  No cosmological
 * effects are included in this standard siren.
 *
 * Implements Eq. (D1) of FINDCHIRP for a 1.4 Msun + 1.4 Msun binary neutron
 * star standard siren at an effective distance 1 Mpc.
 *
 * @sa Appendix D of
 * Bruce Allen, Warren G. Anderson, Patrick R. Brady, Duncan A. Brown, and
 * Jolien D. E. Creighton, "FINDCHIRP: An algorithm for detection of
 * gravitational waves from inspiraling compact binaries", Phys. Rev. D @b 85,
 * 122006 (2012) http://dx.doi.org/10.1103/PhysRevD.85.122006
 *
 * @param psd The one-sided detector strain noise power spectral density.
 * @param f_min The lower bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to 0 or a negative value for no
 * lower bound.
 * @param f_max The upper bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to a negative value for
 * no upper bound.
 * @returns The characteristic signal-to-noise ratio of a binary neutron star
 * standard siren at an effective distance of 1 Mpc.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 */
double XLALMeasureStandardSirenSNR(REAL8FrequencySeries *psd, double f_min, double f_max)
{
	size_t k, k_min, k_max;
	double e = 0.0;
	double sum = 0.0;
	double prefac = 1.0;

	/* standard siren constants */
	double Mpc = 1e6 * LAL_PC_SI;
	double m1 = 1.4; /* in solar masses */
	double m2 = 1.4; /* in solar masses */
	double M = m1 + m2;
	double mu = m1 * m2 / M;
	double A_1_Mpc;

	/* find the (inclusive) limits of the summation */
	if (f_min < 0 || f_min <= psd->f0)
		k_min = 1;
	else {
		k_min = round((f_min - psd->f0)/psd->deltaF);
		if (k_min == 0)
			k_min = 1;
	}

	if (f_max < 0 || f_max >= psd->f0 + psd->data->length * psd->deltaF)
		k_max = psd->data->length - 2;
	else {
		k_max = round((f_max - psd->f0)/psd->deltaF);
		if (k_max >= psd->data->length - 2) {
			k_max = psd->data->length - 2;
		}
	}

	/* Kahan's compensated summation algorithm. The summation is done
	 * from lowest to highest frequency under the assumption that high
	 * frequency components tend to add more to the magnitude of the
	 * derivative.  Note that because only half the components of the
	 * Fourier transform are stored a factor of 2 is added after the
	 * sum.  The DC component should only count once, but it does not
	 * contribute anything to the sum so no special case is required to
	 * handle it. */

	for (k = k_min; k <= k_max; ++k) {
		double tmp = sum;
		double f = psd->f0 + k * psd->deltaF;
		double x;

		/* make sure that psd is valid at this point */
		/* if it is 0 then it must be infinity! */
		if (isinf(psd->data->data[k]) || psd->data->data[k] <= 0.0)
			continue; /* no contribution from this component */

		/* what we want to add = f^{-7/3}/S(f) + "error
		 * from last iteration" */
		x = pow(f, -7.0/3.0) / psd->data->data[k] + e;
		/* add */
		sum += x;
		/* negative of what was actually added */
		e = tmp - sum;
		/* what didn't get added, add next time */
		e += x;
	}

	/* because we've only summed the positive frequency components */
	sum *= 2;

	sum *= 2.0 * psd->deltaF * prefac;

	/* from Eq. (3.4b) of FINDCHIRP */
	A_1_Mpc = -sqrt(5.0 / 24.0 / LAL_PI);
	A_1_Mpc *= (LAL_G_SI * LAL_MSUN_SI / LAL_C_SI / LAL_C_SI / Mpc);
	A_1_Mpc *= pow(LAL_PI * LAL_G_SI * LAL_MSUN_SI / LAL_C_SI / LAL_C_SI / LAL_C_SI, -1.0/6.0);
	A_1_Mpc *= sqrt(mu);
	A_1_Mpc *= pow(M, 1.0/3.0);

	return fabs(A_1_Mpc) * sqrt(sum);
}

/**
 * @brief Measures the characteristic signal-to-noise ratio of a gravitational
 * waveform represented in the frequency domain.
 *
 * This routine measures the characteristic signal-to-noise ratio of a signal
 * @p htilde for a detector with a given strain noise power spectral density
 * @p psd.  Only frequency components of the signal between @p f_mim and
 * @p f_max are included in this measurement.  If @p f_min is zero or negative,
 * no lower bound is imposed.  If @p f_max is negative, no upper bound is
 * imposed.
 *
 * The term @e characteristic is used to indicate that this signal-to-noise
 * ratio is the expected value of the signal-to-noise ratio that would be
 * measured for that signal my a perfectly-matched filter in a detector having
 * stationary Gaussian noise with the specified power spectral density.  The
 * signal-to-noise ratio actually recorded by a matched filter is a random
 * variable that depends on the particular instance of the noise.  Thus the
 * @e characteristic signal-to-noise ratio returned by this routine is a
 * property of the signal and the statistical properties of the detector noise.
 *
 * The @e square of the signal-to-noise ratio is given by
 * \f[
 *   (\mbox{signal-to-noise ratio})^2 =
 *   4 \int_{f_{\rm min}}^{f_{\rm max}} \frac{|\tilde{h}(f)|^2}{S_h(f)}\,df
 * \f]
 * where \f$S_h(f)\f$ is the one-sided detector strain noise power spectral
 * density and \f$\tilde{h}(f)\f$ is the Fourier transform of the signal
 * strain time series \f$h(t)\f$
 * \f[
 *   \tilde{h}(f) = \int_{-\infty}^\infty h(t) e^{-2\pi ift}\,dt.
 * \f]
 * The one-sided strain noise power spectral density is defined by
 * \f$\langle\tilde{n}(f)\tilde{n}^\ast(f')\rangle=\frac{1}{2}S_h(|f|)\delta(f-f')\f$
 * where \f$\tilde{n}(f)\f$ is the Fourier transform of the detector noise
 * process \f$n(t)\f$.
 *
 * The discrete versions of these equations are
 * \f[
 *   (\mbox{signal-to-noise ratio})^2 = 4 \Delta f
 *   \sum_{k=k_{\rm min}}^{k_{\rm max}} |\tilde{h}[k]|^2 / S_h[k]
 * \f]
 * where \f$\tilde{h}[k]\f$ for \f$0\le k<N\f$ is the discrete Fourier
 * transform of the discrete time series \f$h[j]=h(j\Delta t)\f$ for
 * \f$0\le j<N\f$:
 * \f[
 *   \tilde{h}[k] = \Delta t \sum_{j=0}^{N-1} h[j] e^{-2\pi ijk/N}
 * \f]
 * (note the factor of \f$\Delta t\f$).  Here \f$\Delta t\f$ is the sampling
 * interval in time, \f$\Delta f=1/(N\Delta t)\f$ is the sampling interval in
 * frequency, and \f$N\f$ is the number of points in the time series.  The
 * discrete one-sided detector strain noise power spectral density is
 * \f$S_h[k]=2\Delta f\langle|\tilde{n}[k]|^2\rangle\f$ where
 * \f$\tilde{n}[k]\f$ is the discrete Fourier transform of the detector noise
 * process \f$n[j]\f$.
 *
 * The limits of the summation are \f$k_{\rm min}=f_{\rm min}/\Delta f\f$ and
 * \f$k_{\rm max}=f_{\rm max}/\Delta f\f$, both rounded to the nearest integer.
 * If \f$k_{\rm min}\f$ is less than 0, it is set to 0.  If \f$k_{\rm max}\f$
 * is negative or greater than \f$N/2\f$ rounded down, it is set to \f$N/2\f$
 * rounded down.  This ends up double-counting the DC and a possible Nyquist
 * component, but it is assumed these terms have negligible contribution to the
 * sum and there is no special case to handle them separately (most likely,
 * the detector will have no sensitivity to those components anyway).  
 *
 * This routine will accept frequency series \f$\tilde{h}\f$ and \f$S_h\f$ that
 * have different frequency resolutions (i.e., different \f$\Delta f\f$).  In
 * this case, a new frequency series \f$S_h\f$ is computed with the same
 * resolution as \f$\tilde{h}\f$ by interpolation.  A simple linear
 * interpolation in the logarithm of the power spectral density is used.  Any
 * points where the power spectral density is zero are considered to be invalid
 * and are omitted from the sum.
 *
 * @param htilde The Fourier transform of the signal strain.
 * @param psd The one-sided detector strain noise power spectral density.
 * @param f_min The lower bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to 0 or a negative value for no
 * lower bound.
 * @param f_max The upper bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to a negative value for
 * no upper bound.
 * @returns The characteristic signal to noise ratio.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 */
double XLALMeasureSNRFD(COMPLEX16FrequencySeries *htilde, REAL8FrequencySeries *psd, double f_min, double f_max)
{
	LALUnit snrUnits;
	REAL8FrequencySeries *S;
	size_t k, k_min, k_max;
	double e = 0.0;
	double sum = 0.0;
	double prefac = 1.0;

	XLAL_CHECK_REAL8(htilde && psd, XLAL_EFAULT);
	XLAL_CHECK_REAL8(htilde->f0 >= 0.0 && psd->f0 >= 0.0, XLAL_EINVAL, "Can only handle non-negative frequencies");

	/* make sure that snr units will be dimensionless,
	 *
	 *     snr ~ deltaF * |htilde|^2 / psd
	 *
	 * and also compute snr unit prefactor */
	XLALUnitSquare(&snrUnits, &htilde->sampleUnits);
	XLALUnitDivide(&snrUnits, &snrUnits, &psd->sampleUnits);
	XLALUnitMultiply(&snrUnits, &snrUnits, &lalHertzUnit);
	if (!XLALUnitIsDimensionless(&snrUnits))
		XLAL_ERROR_REAL8(XLAL_EINVAL, "Incompatible frequency series: incorrect sample units");
	else
		prefac = XLALUnitPrefactor(&snrUnits);

	/* find the (inclusive) limits of the summation */
	if (f_min < 0 || f_min <= htilde->f0)
		k_min = 0;
	else
		k_min = round((f_min - htilde->f0)/htilde->deltaF);

	if (f_max < 0 || f_max >= htilde->f0 + htilde->data->length * htilde->deltaF)
		k_max = htilde->data->length - 1;
	else {
		k_max = round((f_max - htilde->f0)/htilde->deltaF);
		if (k_max >= htilde->data->length - 1) {
			k_max = htilde->data->length - 1;
		}
	}

	/* interpolate psd if necessary */
	S = create_interpolated_psd(htilde->f0, htilde->deltaF, htilde->data->length, psd);
	if (!S)
		XLAL_ERROR_REAL8(XLAL_EFUNC);

	/* Kahan's compensated summation algorithm. The summation is done
	 * from lowest to highest frequency under the assumption that high
	 * frequency components tend to add more to the magnitude of the
	 * derivative.  Note that because only half the components of the
	 * Fourier transform are stored a factor of 2 is added after the
	 * sum.  The DC component should only count once, but it does not
	 * contribute anything to the sum so no special case is required to
	 * handle it. */

	for (k = k_min; k <= k_max; ++k) {
		double tmp = sum;
		double x;

		/* make sure that psd is valid at this point */
		/* if it is 0 then it must be infinity! */
		if (isinf(S->data->data[k]) || S->data->data[k] <= 0.0)
			continue; /* no contribution from this component */

		/* what we want to add = |\tilde{h}(f)|^{2}/S(f) + "error
		 * from last iteration" */
		x = pow(cabs(htilde->data->data[k]), 2) / S->data->data[k] + e;
		/* add */
		sum += x;
		/* negative of what was actually added */
		e = tmp - sum;
		/* what didn't get added, add next time */
		e += x;
	}

	/* because we've only summed the positive frequency components,
	 * multiply sum by two; note that we assume that the DC and Nyquist
	 * components do not contribute significantly to the sum, so there is
	 * no special case to handle these separately */
	sum *= 2;

	sum *= 2.0 * htilde->deltaF * prefac;

	if (S != psd)
		XLALDestroyREAL8FrequencySeries(S);

	return sqrt(sum);
}

/**
 * @brief Measures the characteristic signal-to-noise ratio of a gravitational
 * waveform.
 *
 * This routine measures the characteristic signal-to-noise ratio of a signal
 * @p h for a detector with a given strain noise power spectral density @p psd.
 * Only frequency components of the signal between @p f_mim and @p f_max are
 * included in this measurement.  If @p f_min is zero or negative, no lower
 * bound is imposed.  If @p f_max is negative, no upper bound is imposed.
 *
 * A copy of the strain time series @p h is zero padded up to the next power of
 * two in length and then Fourier transformed; the routine XLALMeasureSNRFD()
 * is then used to compute the characteristic signal-to-noise ratio.  See
 * XLALMeasureSNRFD() for further details.
 *
 * @param h The strain time series of the signal.
 * @param psd The one-sided detector strain noise power spectral density.
 * @param f_min The lower bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to 0 or a negative value for no
 * lower bound.
 * @param f_max The upper bound of the frequency band over which the
 * signal-to-noise ratio will be computed; set to a negative value for
 * no upper bound.
 * @returns The characteristic signal to noise ratio.
 * @retval LAL_REAL8_FAIL_NAN Failure.
 */
double XLALMeasureSNR(REAL8TimeSeries *h, REAL8FrequencySeries *psd, double f_min, double f_max)
{
	REAL8FFTPlan *plan;
	REAL8TimeSeries *hpadded;
	COMPLEX16FrequencySeries *htilde;
	size_t new_length;
	int length_exp;
	double snr;

	XLAL_CHECK_REAL8(h && psd, XLAL_EFAULT);

	/* zero-pad to the next power of two in length */
	frexp(h->data->length, &length_exp);
	new_length = (size_t)ldexp(1.0, length_exp);

	/* create a new time series padded at the end.
	 * first, copy the series, then resize it. */
	hpadded = XLALCutREAL8TimeSeries(h, 0, h->data->length);
	hpadded = XLALResizeREAL8TimeSeries(hpadded, 0, new_length);
	htilde = XLALCreateCOMPLEX16FrequencySeries(NULL, &hpadded->epoch, 0.0, 0.0, &lalDimensionlessUnit, hpadded->data->length / 2 + 1);
	plan = XLALCreateForwardREAL8FFTPlan(hpadded->data->length, 0);
	if (!hpadded || !plan || !htilde) {
		XLALDestroyREAL8FFTPlan(plan);
		XLALDestroyCOMPLEX16FrequencySeries(htilde);
		XLALDestroyREAL8TimeSeries(hpadded);
		XLAL_ERROR(XLAL_EFUNC);
	}

	if (XLALREAL8TimeFreqFFT(htilde, hpadded, plan) < 0) {
		XLALDestroyREAL8FFTPlan(plan);
		XLALDestroyCOMPLEX16FrequencySeries(htilde);
		XLALDestroyREAL8TimeSeries(hpadded);
		XLAL_ERROR(XLAL_EFUNC);
	}
	XLALDestroyREAL8TimeSeries(hpadded);
	XLALDestroyREAL8FFTPlan(plan);
	
	snr = XLALMeasureSNRFD(htilde, psd, f_min, f_max);

	XLALDestroyCOMPLEX16FrequencySeries(htilde);

	if (XLAL_IS_REAL8_FAIL_NAN(snr))
		XLAL_ERROR_REAL8(XLAL_EFUNC);

	return snr;
}
