/*
 * Copyright (C) 2007,2011,2012  Kipp Cannon
 * Copyright (C) 2012,2013  Reinhard Prix, Karl Wette
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <complex.h>
#include <math.h>


#include <lal/EPSearch.h>
#include <lal/FrequencySeries.h>
#include <lal/LALAtomicDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/Sequence.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/Window.h>
#include <lal/XLALError.h>

static double min(double a, double b) { return a < b ? a : b; }
static double max(double a, double b) { return a > b ? a : b; }


/*
 * ============================================================================
 *
 *                         Channel Filter Management
 *
 * ============================================================================
 */


/**
 * Compute the magnitude of the inner product of two arbitrary channel
 * filters.  Note that the sums are done over only the positive frequency
 * components, so this function multiplies by the required factor of 2.
 * The result is the *full* inner product, not the half inner product.  It
 * is safe to pass the same filter as both arguments.  If the PSD is set to
 * NULL then no PSD weighting is applied.  PSD weighting is only used in
 * reconstructing h_rss.
 *
 * The return value is NaN if the input frequency series have incompatible
 * parameters.  Note that the two-point spectral correlation function does
 * not carry enough metadata to determine if it is compatible with the
 * filters or PSD, for example it does not carry a deltaF parameter.  It is
 * left as an excercise for the calling code to ensure the two-point
 * spectral correlation is appropriate.
 */


double XLALExcessPowerFilterInnerProduct(
	const COMPLEX16FrequencySeries *filter1,	/**< frequency-domain filter */
	const COMPLEX16FrequencySeries *filter2,	/**< frequency-domain filter */
	const REAL8Sequence *correlation,		/**< two-point spectral correlation function.  see XLALREAL8WindowTwoPointSpectralCorrelation(). */
	const REAL8FrequencySeries *psd			/**< power spectral density function.  see XLALREAL8AverageSpectrumWelch() and friends. */
)
{
	const int k10 = round(filter1->f0 / filter1->deltaF);
	const int k20 = round(filter2->f0 / filter2->deltaF);
	const double complex *f1data = filter1->data->data;
	const double complex *f2data = filter2->data->data;
	const double *pdata = psd ? psd->data->data - (int) round(psd->f0 / psd->deltaF) : NULL;
	int k1, k2;
	double complex sum = 0;

	/*
	 * check that filters have same frequency resolution, and if a PSD
	 * is provided that it also has the same frequency resolution and
	 * spans the frequencies spanned by the fitlers
	 */

	if(filter1->deltaF != filter2->deltaF || (psd &&
		(psd->deltaF != filter1->deltaF || psd->f0 > min(filter1->f0, filter2->f0) || max(filter1->f0 + filter1->data->length * filter1->deltaF, filter2->f0 + filter2->data->length * filter2->deltaF) > psd->f0 + psd->data->length * psd->deltaF)
	)) {
		XLALPrintError("%s(): filters are incompatible or PSD does not span filters' frequencies", __func__);
		XLAL_ERROR_REAL8(XLAL_EINVAL);
	}

	/*
	 * compute and return inner product
	 */

	for(k1 = 0; k1 < (int) filter1->data->length; k1++) {
		for(k2 = 0; k2 < (int) filter2->data->length; k2++) {
			const unsigned delta_k = abs(k10 + k1 - k20 - k2);
			double sksk = (delta_k & 1 ? -1 : +1) * (delta_k < correlation->length ? correlation->data[delta_k] : 0);

			if(pdata)
				sksk *= sqrt(pdata[k10 + k1] * pdata[k20 + k2]);

			sum += sksk * f1data[k1] * conj(f2data[k2]);
		}
	}

	return 2 * cabs(sum);
}


/**
 * Generate the frequency domain channel filter function.  The filter
 * corresponds to a frequency band [channel_flow, channel_flow +
 * channel_width].  The filter is nominally a Hann window twice the
 * channel's width, centred on the channel's centre frequency.  This makes
 * a sum across channels equivalent to constructing a Tukey window spanning
 * the same frequency band.  This trick is one of the ingredients that
 * allows us to accomplish a multi-resolution tiling using a single
 * frequency channel projection (*).
 *
 * The filter is normalized so that its "magnitude" as defined by the inner
 * product function XLALExcessPowerFilterInnerProduct() is N.  Then the
 * filter is divided by the square root of the PSD frequency series prior
 * to normalilization.  This has the effect of de-emphasizing frequency
 * bins with high noise content, and is called "over whitening".
 *
 * Note:  the number of samples in the window is odd, being one more than
 * the number of frequency bins in twice the channel width.  This gets the
 * Hann windows to super-impose to form a Tukey window.  (you'll have to
 * draw yourself a picture).
 *
 * (*) Really, there's no need for the "effective window" resulting from
 * summing across channels to be something that has a name, any channel
 * filter at all would do, but this way the code's behaviour is more easily
 * understood --- it's easy to say "the channel filter is a Tukey window of
 * variable central width".
 */


COMPLEX16FrequencySeries *XLALCreateExcessPowerFilter(
	double channel_flow,			/**< Hz */
	double channel_width,			/**< Hz */
	const REAL8FrequencySeries *psd,	/**< power spectral density function.  see XLALREAL8AverageSpectrumWelch() and friends. */
	const REAL8Sequence *correlation	/**< two-point spectral correlation function.  see XLALREAL8WindowTwoPointSpectralCorrelation(). */
)
{
	char filter_name[100];
	REAL8Window *hann;
	COMPLEX16FrequencySeries *filter;
	unsigned i;
	double norm;

	/*
	 * create frequency series for filter
	 */

	sprintf(filter_name, "channel %g +/- %g Hz", channel_flow + channel_width / 2, channel_width / 2);
	filter = XLALCreateCOMPLEX16FrequencySeries(filter_name, &psd->epoch, channel_flow - channel_width / 2, psd->deltaF, &lalDimensionlessUnit, 2 * channel_width / psd->deltaF + 1);
	if(!filter)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	if(filter->f0 < 0.0) {
		XLALPrintError("%s(): channel_flow - channel_width / 2 >= 0.0 failed", __func__);
		XLALDestroyCOMPLEX16FrequencySeries(filter);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}

	/*
	 * build real-valued Hann window and copy into filter
	 */

	hann = XLALCreateHannREAL8Window(filter->data->length);
	if(!hann) {
		XLALDestroyCOMPLEX16FrequencySeries(filter);
		XLALDestroyREAL8Window(hann);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	for(i = 0; i < filter->data->length; i++)
		filter->data->data[i] = hann->data->data[i];
	XLALDestroyREAL8Window(hann);

	/*
	 * divide by square root of PSD to whiten
	 */

	if(!XLALWhitenCOMPLEX16FrequencySeries(filter, psd)) {
		XLALDestroyCOMPLEX16FrequencySeries(filter);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/*
	 * normalize the filter.  the filter needs to be normalized so that
	 * it's inner product with itself is (width / delta F), the width
	 * of the filter in bins.
	 */

	norm = XLALExcessPowerFilterInnerProduct(filter, filter, correlation, NULL);
	if(XLAL_IS_REAL8_FAIL_NAN(norm)) {
		XLALDestroyCOMPLEX16FrequencySeries(filter);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}
	norm = sqrt(channel_width / filter->deltaF / norm);
	for(i = 0; i < filter->data->length; i++)
		filter->data->data[i] *= norm;

	/*
	 * success
	 */

	return filter;
}
