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
 * in the time domain and then transforming to the frequency domain.  This
 * filter should depends on the frequency flow and the fbins_per_channel.
 * It should be non-zero for some amount that depends on fbins_per_channel
 * and the total bandwidth of the input frequency series.
 */

static COMPLEX8FrequencySeries *generate_filter(const COMPLEX8FrequencySeries *template, INT4 fbins_per_channel, REAL8 dt, INT4 fstart)
{
	const char func[] = "generate_filter";
	REAL4Sequence *tdfilter;
	COMPLEX8FrequencySeries *fdfilter;
	RealFFTPlan *plan;
	REAL8 twopiOverNumpts;
	unsigned firstzero;
	REAL4 norm;
	unsigned j;

	/* keep this many samples from the filter on either side of each
	 * channel */
	const INT4 fwindow = 100;

	tdfilter = XLALCreateREAL4Sequence(2 * (template->data->length - 1));
	fdfilter = XLALCreateCOMPLEX8FrequencySeries(NULL, &template->epoch, template->f0, template->deltaF, &template->sampleUnits, template->data->length);
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

	/* normalize the filter */
	norm = sqrt(XLALREAL4SequenceSumSquares(tdfilter, 0, tdfilter->length));
	for(j = 0; j < tdfilter->length; j++)
		tdfilter->data[j] /= norm;

	/* transform to frequency domain */
	if(XLALREAL4ForwardFFT(fdfilter->data, tdfilter, plan)) {
		XLALDestroyREAL4Sequence(tdfilter);
		XLALDestroyCOMPLEX8FrequencySeries(fdfilter);
		XLALDestroyREAL4FFTPlan(plan);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/* extract the part of the filter to be applied to each channel */
	XLALResizeCOMPLEX8FrequencySeries(fdfilter, fstart - fwindow, fbins_per_channel + 2 * fwindow);

	/* clean up and return filter */
	XLALDestroyREAL4Sequence(tdfilter);
	XLALDestroyREAL4FFTPlan(plan);
	return fdfilter;
}


/*
 * Multiply filter by the data.  Don't forget complex conjugate and any other
 * relevant information.
 */

static REAL4 apply_filter(
	COMPLEX8Sequence *outputseq,
	const COMPLEX8FrequencySeries *fseries,
	const COMPLEX8FrequencySeries *filterseries,
	const REAL4FrequencySeries *psdseries
)
{
	size_t fstart = (int) ((filterseries->f0 - fseries->f0) / filterseries->deltaF);
	COMPLEX8 *output = outputseq->data + fstart;
	const COMPLEX8 *fdata = fseries->data->data + fstart;
	const COMPLEX8 *filter = filterseries->data->data;
	const REAL4 *psd = psdseries ? psdseries->data->data + fstart : NULL;
	size_t fbins = filterseries->data->length;
	REAL4 norm = 0.0;

	/* zero the product vector */
	memset(outputseq->data, 0, outputseq->length * sizeof(*outputseq->data));

	/* leave Nyquist at 0 */
	if(fstart + fbins >= outputseq->length)
		fbins = outputseq->length - fstart - 1;

	/* compute:
	 * 	output = fseries * conj(filter) / sqrt(psd)
	 * 	norm = sum(filter * conj(filter) / |psd|)
	 */
	for(; fbins--; output++, fdata++, filter++) {
		REAL4 reFilter = filter->re;
		REAL4 imFilter = filter->im;

		if(psd) {
			reFilter /= sqrt(*psd);
			imFilter /= sqrt(*psd);
			psd++;
		}

		output->re = reFilter * fdata->re + imFilter * fdata->im;
		output->im = reFilter * fdata->im - imFilter * fdata->re;

		norm += reFilter * reFilter + imFilter * imFilter;
	}

	return sqrt(4.0 * norm);
}


/*
 * Evaluate the h_rss estimation factor for each channel
 */

REAL8 *XLALTFPlaneEvalHrssFactor(
	const REAL4TimeFrequencyPlane *plane,
	const COMPLEX8FrequencySeries *response,
	const REAL4FrequencySeries *psd
)
{
	const char func[] = "XLALTFPlaneEvalHrssFactor";
	int fbins_per_channel = plane->deltaF / response->deltaF;
	int fstart = (plane->flow - response->f0) / response->deltaF;
	int tserieslength = 2 * (response->data->length - 1);
	const COMPLEX8 *r = response->data->data + fstart;
	const REAL4 *p = psd->data->data + fstart;
	REAL8 *hrssfactor;
	REAL8 *h;
	int i, j;

	/*
	 * Check that the response is compatible with the spectrum.
	 */

	if((response->f0 != psd->f0) ||
	   (response->deltaF != psd->deltaF) ||
	   (response->data->length != psd->data->length))
	   	XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/*
	 * Allocate memory
	 */

	hrssfactor = LALMalloc(plane->freqBins * sizeof(*hrssfactor));
	if(!hrssfactor)
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);

	/*
	 * Compute the factors
	 */

	h = hrssfactor;
	for(i = 0; i < plane->freqBins; i++, h++) {
		*h = 0.0;
		for(j = 0; j < fbins_per_channel; j++, r++, p++)
			*h += sqrt((r->re * r->re + r->im * r->im) * *p);
		*h *= 0.5 / (fbins_per_channel * tserieslength) / sqrt(psd->deltaF);
	}

	return hrssfactor;
}


/******** <lalVerbatim file="FreqSeriesToTFPlaneCP"> ********/
int XLALFreqSeriesToTFPlane(
	REAL4TimeFrequencyPlane *plane,
	REAL4 *normalisation,
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
	INT4 i;
	INT4 j;
	INT4 nt = plane->timeBins;
	INT4 channels = plane->freqBins;
	INT4 tstart;
	INT4 fstart;
	INT4 fbins_per_channel;
	INT4 tbins_per_sample;
	REAL8 dt;

	/* check input parameters */
	if((nt <= 0) ||
	   (channels <= 0) ||
	   (fseries->deltaF <= 0.0) ||
	   (plane->deltaT <= 0.0) ||
	   (fseries->f0 < 0.0) ||
	   (plane->flow < fseries->f0) ||
	   (fmod(plane->deltaF, fseries->deltaF) != 0.0) ||
	   (fmod(plane->flow - fseries->f0, fseries->deltaF) != 0.0))
		XLAL_ERROR(func, XLAL_EDOM);

	/* number of input frequency series bins per frequency channel in
	 * the time-frequency plane. */
	fbins_per_channel = (int) (plane->deltaF / fseries->deltaF);

	/* time-frequency plane's low frequency cutoff relative to the
	 * input series' heterodyne frequency, in units of frequency bins
	 * */
	fstart = (int) ((plane->flow - fseries->f0) / fseries->deltaF);

	/* make sure we have enough data points in freq series */
	if(fstart + channels * fbins_per_channel > (INT4) fseries->data->length)
		XLAL_ERROR(func, XLAL_EDATA);

	/* create vectors and FFT plans */
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
	tstart = (snr->length - nt * tbins_per_sample) / 2;

	/* set the name and epoch of the TF plane */
	strncpy(plane->name, fseries->name, LALNameLength);
	plane->epoch = fseries->epoch;
	XLALGPSAdd(&plane->epoch, tstart * dt);

	/* generate the frequency domain filter function. */
	filter = generate_filter(fseries, fbins_per_channel, dt, fstart);
	if(!filter) {
		XLALDestroyCOMPLEX8Sequence(fcorr);
		XLALDestroyREAL4Sequence(snr);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* loop over time-frequency plane's channels */
	for(i = 0; i < channels; i++, filter->f0 += fbins_per_channel * filter->deltaF) {
		*normalisation++ = apply_filter(fcorr, fseries, filter, psd);

		/* PRB - Inverse transform the product so that we get a
		 * time series at the full sample rate. */
		if(XLALREAL4ReverseFFT(snr, fcorr, reverseplan)) {
			XLALDestroyCOMPLEX8Sequence(fcorr);
			XLALDestroyREAL4Sequence(snr);
			XLALDestroyCOMPLEX8FrequencySeries(filter);
			XLAL_ERROR(func, XLAL_EFUNC);
		}

		/* PRB - copy the data back into the time series.  In this
		 * process,  one only takes the samples corresponding to the
		 * uncorupted times.   This means that the data should be
		 * longer than the 1/dfmin.   This can mean a change compared
		 * to what has been done before.   We'll have to check that
		 * carefully.  Notice that the originally complex plane is now
		 * real.  I still don't understand why the difference arises
		 * with the Eanns's original implementation.   We'll have to
		 * look at that.  */
		/* Copy the result into appropriate spot in output
		 * structure, skipping the first tstart samples in each
		 * channel */
		for(j = 0; j < nt; j++)
			plane->data[j * channels + i] = snr->data[j * tbins_per_sample + tstart];
	}

	/* Get rid of all temporary memory */
	XLALDestroyCOMPLEX8FrequencySeries(filter);
	XLALDestroyCOMPLEX8Sequence(fcorr);
	XLALDestroyREAL4Sequence(snr);

	/* normal exit */
	return 0;
}
