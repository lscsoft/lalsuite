/******** <lalVerbatim file="FreqSeriesToTFPlaneCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

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
 * normalilization.  This has the effect of demphasizing frequency bins
 * with high noise content, and is called "over whitening".
 */


static COMPLEX8FrequencySeries *generate_filter(
	const COMPLEX8FrequencySeries *template,
	REAL8 channel_flow,
	REAL8 channel_width,
	const REAL4FrequencySeries *psd
)
{
	const char func[] = "generate_filter";
	char filter_name[100];
	REAL4Window *hann;
	COMPLEX8FrequencySeries *fdfilter;
	unsigned i;
	REAL4 norm;

	sprintf(filter_name, "channel %g +/- %g Hz", channel_flow + channel_width / 2, channel_width);

	/*
	 * Channel filter is a Hann window twice the channel's width,
	 * centred on the channel's centre frequency.
	 */

	fdfilter = XLALCreateCOMPLEX8FrequencySeries(filter_name, &template->epoch, channel_flow - channel_width / 2, template->deltaF, &lalDimensionlessUnit, 2 * channel_width / template->deltaF);
	hann = XLALCreateHannREAL4Window(fdfilter->data->length);
	if(!fdfilter || !hann) {
		XLALDestroyCOMPLEX8FrequencySeries(fdfilter);
		XLALDestroyREAL4Window(hann);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < fdfilter->data->length; i++) {
		fdfilter->data->data[i].re = hann->data->data[i];
		fdfilter->data->data[i].im = 0.0;
	}
	XLALDestroyREAL4Window(hann);

	/*
	 * divide by PSD if needed
	 */

	if(psd) {
		REAL4 *pdata = psd->data->data + (int) ((fdfilter->f0 - psd->f0) / psd->deltaF);
		for(i = 0; i < fdfilter->data->length; i++) {
			fdfilter->data->data[i].re /= sqrt(pdata[i]);
			fdfilter->data->data[i].im /= sqrt(pdata[i]);
		}
	}

	/*
	 * normalize the filter.  the filter needs to be normalized so that
	 * it's sum squares is 1, but where the sum is done over all
	 * frequency components positive and negative.  since we are
	 * dealing with real data, our frequency series and this filter
	 * contain only the positive frequency components (the negative
	 * frequency components being known as the complex conjugates of
	 * the positive frequencies).  therefore the sum squares function
	 * is only summing the positive frequencies, and so the value
	 * returned needs to be doubled;  hence the factor of 2.
	 */

	norm = sqrt(2 * XLALCOMPLEX8SequenceSumSquares(fdfilter->data, 0, fdfilter->data->length));
	for(i = 0; i < fdfilter->data->length; i++) {
		fdfilter->data->data[i].re /= norm;
		fdfilter->data->data[i].im /= norm;
	}

	/*
	 * success
	 */

	return fdfilter;
}


/*
 * Compute the magnitude of the inner product of neighbouring channel
 * filters.  assumes filter2->f0 >= filter1->f0, and that the highest
 * numbered sample in filter1 is not past the end of filter2.
 */


static REAL8 filter_overlap(
	COMPLEX8FrequencySeries *filter1,
	COMPLEX8FrequencySeries *filter2
)
{
	COMPLEX8 *f1data = filter1->data->data + (int) ((filter2->f0 - filter1->f0) / filter1->deltaF);
	COMPLEX8 *f1end = filter1->data->data + filter1->data->length;
	COMPLEX8 *f2data = filter2->data->data;
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
	const char func[] = "apply_filter";
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
 * filter.  PSD's computed by LAL obey the convention that for Gaussian
 * noise, the mean square of a frequency bin is psd[k] / (2 * deltaF);
 */


static REAL8 channel_mean_square(
	const REAL4FrequencySeries *psd,
	const COMPLEX8FrequencySeries *filter
)
{
	REAL4 *pdata = psd->data->data + (int) ((filter->f0 - psd->f0) / psd->deltaF);
	COMPLEX8 *fdata = filter->data->data;
	double sum = 0.0;
	unsigned i;

	for(i = 0; i < filter->data->length; i++, pdata++, fdata++)
		sum += *pdata * (fdata->re * fdata->re + fdata->im * fdata->im);

	return sum / psd->deltaF;
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
	const char func[] = "XLALFreqSeriesToTFPlane";
	REAL4Sequence *snr;
	COMPLEX8FrequencySeries **filter;
	COMPLEX8Sequence *fcorr;
	unsigned i;
	unsigned j;
	unsigned tstart;
	unsigned tbins_per_sample;
	REAL8 dt;

	/* check input parameters */
	if((fmod(plane->deltaF, fseries->deltaF) != 0.0) ||
	   (fmod(plane->flow - fseries->f0, fseries->deltaF) != 0.0))
		XLAL_ERROR(func, XLAL_EINVAL);

	/* make sure the frequency series spans an appropriate band */
	if((plane->flow < fseries->f0) ||
	   (plane->flow + plane->channels * plane->deltaF > fseries->f0 + fseries->data->length * fseries->deltaF))
		XLAL_ERROR(func, XLAL_EDATA);

	/* create temporary vectors */
	filter = LALMalloc(plane->channels * sizeof(*filter));
	fcorr = XLALCreateCOMPLEX8Sequence(fseries->data->length);
	snr = XLALCreateREAL4Sequence(2 * (fseries->data->length - 1));
	if(!filter || !fcorr || !snr) {
		LALFree(filter);
		XLALDestroyCOMPLEX8Sequence(fcorr);
		XLALDestroyREAL4Sequence(snr);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* sampling rate of time series which gave fseries */
	dt = 1.0 / (snr->length * fseries->deltaF);
	if(fmod(plane->deltaT, dt) != 0.0) {
		LALFree(filter);
		XLALDestroyCOMPLEX8Sequence(fcorr);
		XLALDestroyREAL4Sequence(snr);
		XLAL_ERROR(func, XLAL_EDOM);
	}

	/* number of input time series bins per sample in each of the
	 * time-frequency plane's channels */
	tbins_per_sample = (int) (plane->deltaT / dt);

	/* the time-frequency plane spans less time than the original time
	 * series because the intent is to skip some amount of data at the
	 * start and end to avoid noise from edge effects;  figure out how
	 * many samples at the start of the original time series need to be
	 * skipped in order to centre the time-frequency plane within it.
	 * */
	tstart = (snr->length - plane->channel[0]->length * tbins_per_sample) / 2;

	/* set the name and epoch of the TF plane */
	strncpy(plane->name, fseries->name, LALNameLength);
	plane->epoch = fseries->epoch;
	XLALGPSAdd(&plane->epoch, tstart * dt);

	/* generate the frequency domain filter functions */
	for(i = 0; i < plane->channels; i++) {
		filter[i] = generate_filter(fseries, plane->flow + i * plane->deltaF, plane->deltaF, enable_over_whitening ? psd : NULL);
		if(!filter[i]) {
			while(--i)
				XLALDestroyCOMPLEX8FrequencySeries(filter[i]);
			LALFree(filter);
			XLALDestroyCOMPLEX8Sequence(fcorr);
			XLALDestroyREAL4Sequence(snr);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
	}

	/* compute the channel overlaps */
	for(i = 0; i < plane->channels - 1; i++)
		plane->channel_overlap->data[i] = filter_overlap(filter[i], filter[i + 1]);

	/* loop over time-frequency plane's channels */
	for(i = 0; i < plane->channels; i++) {
		/* cross correlate the input data against the channel
		 * filter by taking their product in the frequency domain
		 * and then inverse transforming to the time domain to
		 * obtain an SNR time series.  Note that
		 * XLALREAL4ReverseFFT() omits the factor of 1/N in the
		 * inverse transform. */
		apply_filter(fcorr, fseries, filter[i]);
		if(XLALREAL4ReverseFFT(snr, fcorr, reverseplan)) {
			for(i = 0; i < plane->channels; i++)
				XLALDestroyCOMPLEX8FrequencySeries(filter[i]);
			LALFree(filter);
			XLALDestroyCOMPLEX8Sequence(fcorr);
			XLALDestroyREAL4Sequence(snr);
			XLAL_ERROR(func, XLAL_EFUNC);
		}

		/* Down-sample the SNR data into the time-frequency plane
		 * storage arrays.  Note that the first tstart samples of
		 * the SNR time series are skipped. */
		for(j = 0; j < plane->channel[i]->length; j++)
			plane->channel[i]->data[j] = snr->data[tstart + j * tbins_per_sample];

		/* Store the predicted mean square for this channel */
		plane->channel_mean_square->data[i] = channel_mean_square(psd, filter[i]);
	}

	/* clean up */
	for(i = 0; i < plane->channels; i++)
		XLALDestroyCOMPLEX8FrequencySeries(filter[i]);
	LALFree(filter);
	XLALDestroyCOMPLEX8Sequence(fcorr);
	XLALDestroyREAL4Sequence(snr);

	/* success */
	return 0;
}
