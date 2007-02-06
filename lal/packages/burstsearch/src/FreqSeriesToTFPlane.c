/******** <lalVerbatim file="FreqSeriesToTFPlaneCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID(FREQSERIESTOTFPLANEC, "$Id$");

#include <math.h>
#include <lal/Date.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/TFTransform.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>


/*
 * Generate the frequency domain filter function by constructing the filter
 * in the time domain and then transforming to the frequency domain.
 */


static COMPLEX8FrequencySeries *generate_filter(const COMPLEX8FrequencySeries *template, REAL8 channel_width, REAL8 dt, INT4 fstart)
{
	const char func[] = "generate_filter";
	REAL4Sequence *tdfilter;
	COMPLEX8FrequencySeries *fdfilter;
	RealFFTPlan *plan;
	REAL8 twopiOverNumpts;
	INT4 fbins_per_channel = channel_width / template->deltaF;
	unsigned firstzero;
	REAL4 norm;
	unsigned j;

	/* keep this many samples from the filter on either side of each
	 * channel */
	const INT4 fwindow = 100;

	tdfilter = XLALCreateREAL4Sequence(2 * (template->data->length - 1));
	fdfilter = XLALCreateCOMPLEX8FrequencySeries(NULL, &template->epoch, 0.0, template->deltaF, &template->sampleUnits, template->data->length);
	plan = XLALCreateForwardREAL4FFTPlan(tdfilter->length, 0);
	if(!tdfilter || !fdfilter || !plan) {
		XLALDestroyREAL4Sequence(tdfilter);
		XLALDestroyCOMPLEX8FrequencySeries(fdfilter);
		XLALDestroyREAL4FFTPlan(plan);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}

	/* zero the time-domain filter */
	memset(tdfilter->data, 0, tdfilter->length * sizeof(*tdfilter->data));

	/* number of points from peak of filter to first zero */
	firstzero = tdfilter->length / fbins_per_channel;

	twopiOverNumpts = 2.0 * LAL_PI / tdfilter->length;
	tdfilter->data[0] = twopiOverNumpts * fbins_per_channel / (LAL_PI * dt);
	for(j = 1; j < firstzero; j++)
		tdfilter->data[j] = tdfilter->data[tdfilter->length - j] = (sin(twopiOverNumpts * j * (fstart + fbins_per_channel)) - sin(twopiOverNumpts * j * fstart)) / (LAL_PI * j * dt);

	/* transform to frequency domain */
	if(XLALREAL4ForwardFFT(fdfilter->data, tdfilter, plan)) {
		XLALDestroyREAL4Sequence(tdfilter);
		XLALDestroyCOMPLEX8FrequencySeries(fdfilter);
		XLALDestroyREAL4FFTPlan(plan);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/* extract the part of the filter to be applied to each channel */
	XLALResizeCOMPLEX8FrequencySeries(fdfilter, fstart - fwindow, fbins_per_channel + 2 * fwindow);

	/* normalize the filter.  the filter needs to be normalized so that
	 * it's sum squares is 1, but where the sum is done over all
	 * frequency components positive and negative.  since we are
	 * dealing with real data, our frequency series and this filter
	 * contain only the positive frequency components (the negative
	 * frequency components being known as the complex conjugates of
	 * the positive frequencies).  therefore the sum squares function
	 * is only summing the positive frequencies, and so the value
	 * returned needs to be doubled;  hence the factor of 2. */
	norm = sqrt(2 * XLALCOMPLEX8SequenceSumSquares(fdfilter->data, 0, fdfilter->data->length));
	for(j = 0; j < fdfilter->data->length; j++) {
		fdfilter->data->data[j].re /= norm;
		fdfilter->data->data[j].im /= norm;
	}

	/* clean up and return filter */
	XLALDestroyREAL4Sequence(tdfilter);
	XLALDestroyREAL4FFTPlan(plan);
	return fdfilter;
}


/*
 * Compute the magnitude of the inner product of neighbouring channel
 * filters.
 */


static REAL8 filter_overlap(
	COMPLEX8FrequencySeries *filter,
	REAL8 channel_spacing
)
{
	COMPLEX8 *fdata = filter->data->data;
	COMPLEX8 sum = {0, 0};
	const unsigned n = filter->data->length - channel_spacing / filter->deltaF;
	unsigned i;

	if(channel_spacing / filter->deltaF > filter->data->length)
		return 0;

	for(i = 0; i < n; i++) {
		sum.re += fdata[n + i].re * fdata[i].re + fdata[n + i].im * fdata[i].im;
		sum.im += fdata[n + i].im * fdata[i].re - fdata[n + i].re * fdata[i].im;
	}

	return sqrt(sum.re * sum.re + sum.im * sum.im);
}


/*
 * Multiply the data by the filter.
 */


static COMPLEX8Sequence *apply_filter(
	COMPLEX8Sequence *outputseq,
	const COMPLEX8FrequencySeries *fseries,
	const COMPLEX8FrequencySeries *filterseries
)
{
	size_t fstart = (int) ((filterseries->f0 - fseries->f0) / filterseries->deltaF);
	COMPLEX8 *output = outputseq->data + fstart;
	const COMPLEX8 *input = fseries->data->data + fstart;
	const COMPLEX8 *filter = filterseries->data->data;
	/* leave Nyquist at 0 */
	size_t fbins = outputseq->length - fstart - 1;
	if(filterseries->data->length < fbins)
		fbins = filterseries->data->length;

	if(outputseq->length != fseries->data->length)
		return NULL;

	/* zero the product vector */
	memset(outputseq->data, 0, outputseq->length * sizeof(*outputseq->data));

	/* output = fseries * conj(filter) */
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


static REAL4 channel_mean_square(
	const REAL4FrequencySeries *psd,
	const COMPLEX8FrequencySeries *filter
)
{
	REAL4 *pdata = psd->data->data + (int) ((filter->f0 - psd->f0) / psd->deltaF);
	COMPLEX8 *fdata = filter->data->data;
	REAL4 sum = 0.0;
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
	const REAL4FFTPlan *reverseplan
)
/******** </lalVerbatim> ********/
{
	const char func[] = "XLALFreqSeriesToTFPlane";
	REAL4Sequence *snr;
	COMPLEX8FrequencySeries *filter;
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
	fcorr = XLALCreateCOMPLEX8Sequence(fseries->data->length);
	snr = XLALCreateREAL4Sequence(2 * (fseries->data->length - 1));
	if(!fcorr || !snr) {
		XLALDestroyCOMPLEX8Sequence(fcorr);
		XLALDestroyREAL4Sequence(snr);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* sampling rate of time series which gave fseries */
	dt = 1.0 / (snr->length * fseries->deltaF);
	if(fmod(plane->deltaT, dt) != 0.0) {
		XLALDestroyCOMPLEX8Sequence(fcorr);
		XLALDestroyREAL4Sequence(snr);
		XLAL_ERROR(func, XLAL_EDOM);
	}

	/* number of input time series bins per sample in each of the
	 * time-frequency plane's channel */
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

	/* generate the frequency domain filter function. */
	filter = generate_filter(fseries, plane->deltaF, dt, (plane->flow - fseries->f0) / fseries->deltaF);
	if(!filter) {
		XLALDestroyCOMPLEX8Sequence(fcorr);
		XLALDestroyREAL4Sequence(snr);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* compute the channel overlap */
	plane->channel_overlap = filter_overlap(filter, plane->deltaF);

	/* loop over time-frequency plane's channels */
	for(i = 0; i < plane->channels; i++, filter->f0 += plane->deltaF) {
		/* cross correlate the input data against the channel
		 * filter by taking their product in the frequency domain
		 * and then inverse transforming to the time domain to
		 * obtain an SNR time series.  Note that
		 * XLALREAL4ReverseFFT() omits the factor of 1/N in the
		 * inverse transform. */
		apply_filter(fcorr, fseries, filter);
		if(XLALREAL4ReverseFFT(snr, fcorr, reverseplan)) {
			XLALDestroyCOMPLEX8FrequencySeries(filter);
			XLALDestroyCOMPLEX8Sequence(fcorr);
			XLALDestroyREAL4Sequence(snr);
			XLAL_ERROR(func, XLAL_EFUNC);
		}

		/* Down-sample the SNR data into the time-frequency plane
		 * storage arrays.  Note that the first tstart samples of
		 * the SNR time series are skipped. */
		for(j = 0; j < plane->channel[i]->length; j++)
			plane->channel[i]->data[j] = snr->data[j * tbins_per_sample + tstart];

		/* Store the predicted mean square for this channel */
		plane->channel_mean_square->data[i] = channel_mean_square(psd, filter);
	}

	/* Get rid of all temporary memory */
	XLALDestroyCOMPLEX8FrequencySeries(filter);
	XLALDestroyCOMPLEX8Sequence(fcorr);
	XLALDestroyREAL4Sequence(snr);

	/* normal exit */
	return 0;
}
