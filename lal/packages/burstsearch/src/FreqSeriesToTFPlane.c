/*
 * $Id$
 *
 * Copyright (C) 2007  Kipp Cannon and Flanagan, E.
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


#include <lal/LALRCSID.h>


NRCSID(FREQSERIESTOTFPLANEC, "$Id$");


#include <math.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/Units.h>


/*
 * Generate the frequency domain channel filter function.  The filter is
 * nominally a Hann window twice the channel's width, centred on the
 * channel's centre frequency.  The filter is normalized so that its sum
 * squares is 1.  If the psd parameter is not NULL, then the filter is
 * divided by the square root of this frequency series prior to
 * normalilization.  This has the effect of de-emphasizing frequency bins
 * with high noise content, and is called "over whitening".
 */


static COMPLEX8FrequencySeries *generate_filter(
	const COMPLEX8FrequencySeries *template,
	REAL8 channel_flow,
	REAL8 channel_width,
	const REAL4FrequencySeries *psd
)
{
	static const char func[] = "generate_filter";
	char filter_name[100];
	REAL4Window *hann;
	COMPLEX8FrequencySeries *filter;
	unsigned i;
	REAL4 norm;

	sprintf(filter_name, "channel %g +/- %g Hz", channel_flow + channel_width / 2, channel_width / 2);

	/*
	 * Channel filter is a Hann window twice the channel's width,
	 * centred on the channel's centre frequency.  This makes a sum
	 * across channels equivalent to constructing a Tukey window
	 * spanning the same frequency band.  This trick is one of the
	 * ingredients that allows us to accomplish a multi-resolution
	 * tiling using a single frequency channel projection.  Really,
	 * there's no need for the "effective window" resulting from
	 * summing across channels to be something that has a name, any
	 * channel filter at all would do, but this way the code's
	 * behaviour is more easily understood --- it's easy to say "the
	 * channel filter is a Tukey window of variable centre width".
	 *
	 * Note:  the number of samples in the window is odd, being one
	 * more than the number of frequency bins in twice the channel
	 * width.  This gets the Hann windows to super-impose to form a
	 * Tukey window.  (you'll have to draw yourself a picture).
	 */

	filter = XLALCreateCOMPLEX8FrequencySeries(filter_name, &template->epoch, channel_flow - channel_width / 2, template->deltaF, &lalDimensionlessUnit, 2 * channel_width / template->deltaF + 1);
	hann = XLALCreateHannREAL4Window(filter->data->length);
	if(!filter || !hann) {
		XLALDestroyCOMPLEX8FrequencySeries(filter);
		XLALDestroyREAL4Window(hann);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < filter->data->length; i++) {
		filter->data->data[i].re = hann->data->data[i];
		filter->data->data[i].im = 0.0;
	}
	XLALDestroyREAL4Window(hann);

	/*
	 * divide by PSD if needed
	 */

	if(psd) {
		REAL4 *pdata = psd->data->data + (int) ((filter->f0 - psd->f0) / psd->deltaF);
		for(i = 0; i < filter->data->length; i++) {
			filter->data->data[i].re /= sqrt(pdata[i]);
			filter->data->data[i].im /= sqrt(pdata[i]);
		}
	}

	/*
	 * normalize the filter.  the filter needs to be normalized so that
	 * it's sum squares is 1, but where the sum is done over all
	 * frequency components positive and negative.  since we are
	 * dealing with real data, our frequency series and this filter
	 * contain only the positive frequency components (the negative
	 * frequency components being the complex conjugates of the
	 * positive frequencies).  therefore the sum squares function is
	 * only summing the positive frequencies, and so the value it
	 * returns needs to be doubled;  hence the factor of 2.
	 */

	norm = sqrt(2 * XLALCOMPLEX8SequenceSumSquares(filter->data, 0, filter->data->length));
	for(i = 0; i < filter->data->length; i++) {
		filter->data->data[i].re /= norm;
		filter->data->data[i].im /= norm;
	}

	/*
	 * success
	 */

	return filter;
}


/*
 * Compute the magnitude of the inner product of two channel filters.
 * assumes filter2->f0 >= filter1->f0, and that the highest numbered sample
 * in filter1 is not past the end of filter2.
 */


static REAL8 filter_overlap(
	const COMPLEX8FrequencySeries *filter1,
	const COMPLEX8FrequencySeries *filter2
)
{
	const COMPLEX8 *f1data = filter1->data->data + (int) ((filter2->f0 - filter1->f0) / filter1->deltaF);
	const COMPLEX8 *f1end = filter1->data->data + filter1->data->length;
	const COMPLEX8 *f2data = filter2->data->data;
	COMPLEX16 sum = {0, 0};

	if((filter1->f0 > filter2->f0) || ((f1end - f1data) > (int) filter2->data->length) || (filter1->deltaF != filter2->deltaF))
		return XLAL_REAL8_FAIL_NAN;

	/* sum filter1 * conj(filter2) */
	for(; f1data < f1end; f1data++, f2data++) {
		sum.re += f1data->re * f2data->re + f1data->im * f2data->im;
		sum.im += f1data->im * f2data->re - f1data->re * f2data->im;
	}

	/* return | sum filter1 * conj(filter2) | */
	return sqrt(sum.re * sum.re + sum.im * sum.im);
}


/*
 * Multiply the data by the filter.
 */


static COMPLEX8Sequence *apply_filter(
	COMPLEX8Sequence *outputseq,
	const COMPLEX8FrequencySeries *inputseries,
	const COMPLEX8FrequencySeries *filterseries
)
{
	static const char func[] = "apply_filter";
	int fstart = (filterseries->f0 - inputseries->f0) / filterseries->deltaF;
	COMPLEX8 *output = outputseq->data + (fstart < 0 ? 0 : fstart);
	const COMPLEX8 *input = inputseries->data->data + (fstart < 0 ? 0 : fstart);
	const COMPLEX8 *filter = filterseries->data->data + (fstart < 0 ? -fstart : 0);
	/* an extra 1 is subtracted to ensure the Nyquist is set to 0 */
	size_t fbins = outputseq->length - (fstart < 0 ? 0 : fstart) - 1;
	if(filterseries->data->length - (fstart < 0 ? -fstart : 0) < fbins)
		fbins = filterseries->data->length - (fstart < 0 ? -fstart : 0);

	if(outputseq->length != inputseries->data->length)
		XLAL_ERROR_NULL(func, XLAL_EBADLEN);

	/* zero the product vector */
	memset(outputseq->data, 0, outputseq->length * sizeof(*outputseq->data));

	/* output = inputseries * conj(filter) */
	for(; fbins--; output++, input++, filter++) {
		output->re = input->re * filter->re + input->im * filter->im;
		output->im = input->im * filter->re - input->re * filter->im;
	}

	return outputseq;
}


/*
 * Compute the mean square for a channel from the PSD and the channel's
 * filter.  PSDs computed by LAL obey the convention that for Gaussian
 * noise, the mean square of a frequency bin is psd[k] / (2 deltaF).
 * Therefore, the mean square of a frequency bin after being multiplied by
 * the channel filter, c[k], is psd[k] |c[k]|^2 / (2 deltaF).  The mean
 * square for the channel is the sum of mean squares for the bins within
 * it, if separate frequency bins are statistically independent so that
 * there are no cross terms.  This is true for stationary noise.
 */


static REAL8 channel_mean_square(
	const REAL4FrequencySeries *psd,
	const COMPLEX8FrequencySeries *filter
)
{
	const REAL4 *pdata = psd->data->data + (int) ((filter->f0 - psd->f0) / psd->deltaF);
	const COMPLEX8 *fdata = filter->data->data;
	double sum = 0.0;
	unsigned i;

	for(i = 0; i < filter->data->length; i++, pdata++, fdata++)
		sum += *pdata * (fdata->re * fdata->re + fdata->im * fdata->im);

	return sum / (2 * psd->deltaF);
}


/******** <lalVerbatim file="FreqSeriesToTFPlaneCP"> ********/
int XLALFreqSeriesToTFPlane(
	REAL4TimeFrequencyPlane *plane,
	const COMPLEX8FrequencySeries *fseries,
	const REAL4FrequencySeries *psd,
	const REAL4FFTPlan *reverseplan,
	INT4 enable_over_whitening
)
/******** </lalVerbatim> ********/
{
	static const char func[] = "XLALFreqSeriesToTFPlane";
	COMPLEX8FrequencySeries **filter;
	COMPLEX8Sequence *fcorr;
	unsigned i;

	/* check input parameters */
	if((fmod(plane->deltaF, fseries->deltaF) != 0.0) ||
	   (fmod(plane->flow - fseries->f0, fseries->deltaF) != 0.0))
		XLAL_ERROR(func, XLAL_EINVAL);

	/* make sure the frequency series spans an appropriate band */
	if((plane->flow < fseries->f0) ||
	   (plane->flow + plane->channels * plane->deltaF > fseries->f0 + fseries->data->length * fseries->deltaF))
		XLAL_ERROR(func, XLAL_EDATA);

	/* create temporary vectors */
	filter = XLALMalloc(plane->channels * sizeof(*filter));
	fcorr = XLALCreateCOMPLEX8Sequence(fseries->data->length);
	if(!filter || !fcorr) {
		XLALFree(filter);
		XLALDestroyCOMPLEX8Sequence(fcorr);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	XLALPrintInfo("XLALFreqSeriesToTFPlane(): generating channel filters\n");
	/* generate the frequency domain filter functions */
	for(i = 0; i < plane->channels; i++) {
		filter[i] = generate_filter(fseries, plane->flow + i * plane->deltaF, plane->deltaF, enable_over_whitening ? psd : NULL);
		if(!filter[i]) {
			while(--i)
				XLALDestroyCOMPLEX8FrequencySeries(filter[i]);
			XLALFree(filter);
			XLALDestroyCOMPLEX8Sequence(fcorr);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
	}

	/* compute the channel overlaps */
	for(i = 0; i < plane->channels - 1; i++)
		plane->twice_channel_overlap->data[i] = 2 * filter_overlap(filter[i], filter[i + 1]);

	XLALPrintInfo("XLALFreqSeriesToTFPlane(): projecting data onto time-frequency plane\n");
	/* loop over the time-frequency plane's channels */
	for(i = 0; i < plane->channels; i++) {
		/* cross correlate the input data against the channel
		 * filter by taking their product in the frequency domain
		 * and then inverse transforming to the time domain to
		 * obtain an SNR time series.  Note that
		 * XLALREAL4ReverseFFT() omits the factor of 1/N in the
		 * inverse transform. */
		apply_filter(fcorr, fseries, filter[i]);
		if(XLALREAL4ReverseFFT(plane->channel[i], fcorr, reverseplan)) {
			for(i = 0; i < plane->channels; i++)
				XLALDestroyCOMPLEX8FrequencySeries(filter[i]);
			XLALFree(filter);
			XLALDestroyCOMPLEX8Sequence(fcorr);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
		{
		/* correct for bias resulting from time-domain window used
		 * to put tapers on time series */
		unsigned j;
		for(j = 0; j < plane->channel[i]->length; j++)
			plane->channel[i]->data[j] *= sqrt(plane->window->sumofsquares / plane->window->data->length);
		}

		/* Store the expected root mean square for this channel */
		plane->channel_rms->data[i] = sqrt(channel_mean_square(psd, filter[i]));
	}

	/* clean up */
	for(i = 0; i < plane->channels; i++)
		XLALDestroyCOMPLEX8FrequencySeries(filter[i]);
	XLALFree(filter);
	XLALDestroyCOMPLEX8Sequence(fcorr);

	/* set the name and epoch of the TF plane */
	strncpy(plane->name, fseries->name, LALNameLength);
	plane->epoch = fseries->epoch;

	/* success */
	return 0;
}
