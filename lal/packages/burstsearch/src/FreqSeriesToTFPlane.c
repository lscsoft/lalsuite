/******** <lalVerbatim file="FreqSeriesToTFPlaneCV"> ********
Author: Flanagan, E. and Cannon, K.
$Id$
********* </lalVerbatim> ********/

#include <lal/LALRCSID.h>

NRCSID(FREQSERIESTOTFPLANEC, "$Id$");

#include <math.h>
#include <lal/Date.h>
#include <lal/LALConstants.h>
#include <lal/RealFFT.h>
#include <lal/TFTransform.h>
#include <lal/Sequence.h>


/*
 * Generate the frequency domain filter function by constructing the filter
 * in the time domain and then transforming to the frequency domain.  This
 * filter should depends on the frequency flow and the fseglength.  It
 * should be non-zero for some amount that depends on fseglength and the
 * total bandwidth of the input frequency series.
 */

static COMPLEX8Sequence *generate_filter(size_t length, INT4 fseglength, REAL8 dt, INT4 flow)
{
	static const char *func = "generate_filter";
	REAL4Sequence *tdfilter;
	COMPLEX8Sequence *fdfilter;
	RealFFTPlan *plan;
	REAL8 twopiOverNumpts;
	INT4 firstzero;
	int j;

	tdfilter = XLALCreateREAL4Sequence(2 * (length - 1));
	fdfilter = XLALCreateCOMPLEX8Sequence(length);
	plan = XLALCreateForwardREAL4FFTPlan(tdfilter->length, 0);
	if(!tdfilter || !fdfilter || !plan) {
		XLALDestroyREAL4Sequence(tdfilter);
		XLALDestroyCOMPLEX8Sequence(fdfilter);
		XLALDestroyREAL4FFTPlan(plan);
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
	}

	/* zero the time-domain filter */
	memset(tdfilter->data, 0, tdfilter->length * sizeof(*tdfilter->data));

	/* number of points from peak of filter to first zero */
	firstzero = tdfilter->length / fseglength;

	twopiOverNumpts = 2.0 * LAL_PI / tdfilter->length;
	tdfilter->data[0] = twopiOverNumpts * fseglength / (LAL_PI * dt);
	for(j = 1; j < firstzero; j++)
		tdfilter->data[j] = tdfilter->data[tdfilter->length - j] = (sin(twopiOverNumpts * j * (flow + fseglength)) - sin(twopiOverNumpts * j * flow)) / (LAL_PI * j * dt);

	if(XLALREAL4ForwardFFT(fdfilter, tdfilter, plan)) {
		XLALDestroyREAL4Sequence(tdfilter);
		XLALDestroyCOMPLEX8Sequence(fdfilter);
		XLALDestroyREAL4FFTPlan(plan);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	XLALDestroyREAL4Sequence(tdfilter);
	XLALDestroyREAL4FFTPlan(plan);
	return(fdfilter);
}


/*
 * Multiply filter by the data.  Don't forget complex conjugate and any other
 * relevant information.
 */

static REAL4 apply_filter(
	COMPLEX8 *output,
	const COMPLEX8 *fseries,
	const COMPLEX8 *filter,
	const REAL4 *psd,
	size_t start,
	size_t length
)
{
	REAL4 norm = 0.0;

	output += start;
	fseries += start;
	if(psd)
		psd += start;

	for(; length--; output++, fseries++, filter++) {
		REAL4 reFilter = filter->re;
		REAL4 imFilter = filter->im;

		if(psd) {
			reFilter /= sqrt(*psd);
			imFilter /= sqrt(*psd);
			psd++;
		}

		output->re = reFilter * fseries->re + imFilter * fseries->im;
		output->im = reFilter * fseries->im - imFilter * fseries->re;

		norm += reFilter * reFilter + imFilter * imFilter;
	}

	return(sqrt(4.0 * norm));
}


/******** <lalVerbatim file="FreqSeriesToTFPlaneCP"> ********/
int XLALFreqSeriesToTFPlane(
	REAL4TimeFrequencyPlane *plane,
	const COMPLEX8FrequencySeries *freqSeries,
	UINT4 windowShift,
	REAL4 *normalisation,
	const REAL4FrequencySeries *psd
)
/******** </lalVerbatim> ********/
{
	static const char *func = "XLALFreqSeriesToTFPlane";
	REAL4Sequence *snr;
	COMPLEX8Sequence *filter;
	COMPLEX8Sequence *fcorr;
	RealFFTPlan *prev;
	INT4 i;
	INT4 j;
	INT4 nt = plane->params.timeBins;
	INT4 nf = plane->params.freqBins;
	INT4 flow;
	INT4 fseglength;
	REAL8 dt;
	/* the frequency window in bins; for each channel contributions
	 * will be zeroed out from fwindow before and after */
	const INT4 fwindow = 100;

	/* check input parameters */
	if((nt <= 0) ||
	   (nf <= 0) ||
	   (freqSeries->deltaF <= 0.0) ||
	   (plane->params.deltaT <= 0.0) ||
	   (freqSeries->f0 < 0.0) ||
	   (plane->params.flow < freqSeries->f0))
		XLAL_ERROR(func, XLAL_EDOM);

	/* number of bins of the frequency series.  Note that the frequency
	 * resolution of the frequency series does not have to be the same
	 * as the time-frequency plane. */
	fseglength = 0.5 + plane->params.deltaF / freqSeries->deltaF;

	/* time-frequency plane's low frequency cutoff relative to the
	 * input series' heterodyne frequency, in units of frequency bins
	 * */
	flow = (plane->params.flow - freqSeries->f0) / freqSeries->deltaF;

	/* make sure have enough data points in freq series */
	if(flow + nf * fseglength > (INT4) freqSeries->data->length)
		XLAL_ERROR(func, XLAL_EDATA);

	/* create vectors and FFT plans */
	fcorr = XLALCreateCOMPLEX8Sequence(freqSeries->data->length);
	snr = XLALCreateREAL4Sequence(2 * (freqSeries->data->length - 1));
	prev = XLALCreateReverseREAL4FFTPlan(snr->length, 0);
	if(!fcorr || !snr || !prev) {
		XLALDestroyCOMPLEX8Sequence(fcorr);
		XLALDestroyREAL4Sequence(snr);
		XLALDestroyREAL4FFTPlan(prev);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	/* sampling rate of time series which gave freqSeries */
	dt = 1.0 / (snr->length * freqSeries->deltaF);

	/* set the epoch of the TF plane */
	plane->epoch = freqSeries->epoch;

	/* generate the frequency domain filter function. */
	filter = generate_filter(freqSeries->data->length, fseglength, dt, flow);
	if(!filter) {
		XLALDestroyCOMPLEX8Sequence(fcorr);
		XLALDestroyREAL4Sequence(snr);
		XLALDestroyREAL4FFTPlan(prev);
		XLAL_ERROR(func, XLAL_EFUNC);
	}
	
	/* heterodyne the filter */
	XLALShiftCOMPLEX8Sequence(filter, -(flow - fwindow));

	/* loop over different basis vectors in the frequency domain */
	for(i = 0; i < nf; i++, normalisation++) {
		/* zero the product vector */
		memset(fcorr->data, 0, fcorr->length * sizeof(*fcorr->data));

		*normalisation = apply_filter(fcorr->data, freqSeries->data->data, filter->data, psd ? psd->data->data : NULL, flow + i * fseglength - fwindow, fseglength + 2 * fwindow);

		/* clean up Nyquist */
		fcorr->data[fcorr->length - 1].re = 0.0;
		fcorr->data[fcorr->length - 1].im = 0.0;

		/* PRB - Inverse transform the product so that we get a
		 * time series at the full sample rate.   Make sure to
		 * check that the sample rate agrees with what you had
		 * before. */
		if(XLALREAL4ReverseFFT(snr, fcorr, prev)) {
			XLALDestroyCOMPLEX8Sequence(fcorr);
			XLALDestroyREAL4Sequence(snr);
			XLALDestroyREAL4FFTPlan(prev);
			XLALDestroyCOMPLEX8Sequence(filter);
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
		/* Copy the result into appropriate spot in output structure */
		for(j = 0; j < nt; j++)
			plane->data[j * nf + i] = snr->data[j * (INT4) (plane->params.deltaT / dt) + windowShift];
	}

	/* Get rid of all temporary memory */
	XLALDestroyREAL4FFTPlan(prev);
	XLALDestroyCOMPLEX8Sequence(filter);
	XLALDestroyCOMPLEX8Sequence(fcorr);
	XLALDestroyREAL4Sequence(snr);

	/* normal exit */
	return(0);
}
