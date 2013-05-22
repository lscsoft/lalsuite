/*
*  Copyright (C) 2007 Jolien Creighton
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

#include <complex.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <lal/Date.h>
#include <lal/LALConstants.h>
#include <lal/LALStdlib.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/TimeSeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/LALSimNoise.h>


/* 
 * This routine generates a single segment of data.  Note that this segment is
 * generated in the frequency domain and is inverse Fourier transformed into
 * the time domain; consequently the data is periodic in the time domain.
 */
static int XLALSimNoiseSegment(REAL8TimeSeries *s, REAL8FrequencySeries *psd, gsl_rng *rng)
{
	size_t k;
	REAL8FFTPlan *plan;
	COMPLEX16FrequencySeries *stilde;

	plan = XLALCreateReverseREAL8FFTPlan(s->data->length, 0);
	if (! plan)
		XLAL_ERROR(XLAL_EFUNC);

	stilde = XLALCreateCOMPLEX16FrequencySeries("STILDE", &s->epoch, 0.0, 1.0/(s->data->length * s->deltaT), &lalSecondUnit, s->data->length/2 + 1);
	if (! stilde) {
		XLALDestroyREAL8FFTPlan(plan);
		XLAL_ERROR(XLAL_EFUNC);
	}

	XLALUnitMultiply(&stilde->sampleUnits, &stilde->sampleUnits, &s->sampleUnits);

	stilde->data->data[0] = 0.0;
	for (k = 0; k < s->data->length/2 + 1; ++k) {
		double sigma = 0.5 * sqrt(psd->data->data[k] / psd->deltaF);
		stilde->data->data[k] = gsl_ran_gaussian_ziggurat(rng, sigma);
		stilde->data->data[k] += I * gsl_ran_gaussian_ziggurat(rng, sigma);
	}

	XLALREAL8FreqTimeFFT(s, stilde, plan);

	XLALDestroyCOMPLEX16FrequencySeries(stilde);
	XLALDestroyREAL8FFTPlan(plan);
	return 0;
}


/**
 * Routine that may be used to generate sequential segments of data with a
 * specified stride from one segment to the next.
 *
 * Calling instructions: for the first call, set stride = 0; subsequent calls
 * should pass the same time series and have non-zero stride.  This routine
 * will advance the time series by an amount given by the stride and will
 * generate new data so that the data is continuous from one segment to the
 * next.  For example: the following routine will output a continuous stream of
 * detector noise with an Initial LIGO spectrum above 40 Hz:
 *
 * \code
 * #include <stdio.h>
 * #include <gsl/gsl_rng.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/FrequencySeries.h>
 * #include <lal/TimeSeries.h>
 * #include <lal/Units.h>
 * #include <lal/LALSimNoise.h>
 * void mkligodata(void)
 * {
 * 	const double flow = 40.0; // 40 Hz low frequency cutoff
 * 	const double duration = 16.0; // 16 second segments
 * 	const double srate = 16384.0; // sampling rate in Hertz
 * 	size_t length = duration * srate; // segment length
 * 	size_t stride = length / 2; // stride between segments
 * 	LIGOTimeGPS epoch = { 0, 0 };
 * 	REAL8FrequencySeries *psd;
 * 	REAL8TimeSeries *seg;
 *	gsl_rng *rng;
 * 	gsl_rng_env_setup();
 * 	rng = gsl_rng_alloc(gsl_rng_default);
 * 	seg = XLALCreateREAL8TimeSeries("STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
 *	psd = XLALCreateREAL8FrequencySeries("LIGO SRD", &epoch, 0.0, 1.0/duration, &lalSecondUnit, length/2 + 1);
 *	XLALSimNoisePSD(psd, flow, XLALSimNoisePSDiLIGOSRD);
 * 	XLALSimNoise(seg, 0, psd, rng); // first time to initialize
 *	while (1) { // infinite loop
 * 		double t0 = XLALGPSGetREAL8(&seg->epoch);
 * 		size_t j;
 *		for (j = 0; j < stride; ++j) // output first stride points
 * 			printf("%.9f\t%e\n", t0 + j*seg->deltaT, seg->data->data[j]);
 *		XLALSimNoise(seg, stride, psd, rng); // make more data
 * 	}
 * }
 * \endcode
 *
 * If only one single segment of data is required, set stride to be the length
 * of the timeseries data vector.  This will make a single segment of data
 * that is *not* periodic (also, in this case it will not advance the epoch of
 * the timeseries).
 *
 * Note:
 *
 *	- If stride = 0, initialize h by generating one (periodic)
 *	realization of noise; subsequent calls should have non-zero
 *	stride.
 *
 *	- If stride = h->data->length then generate one segment of
 *	non-periodic noise by generating two different realizations
 *	and feathering them together.
 *
 * Warning: only the first stride points are valid.
 */
int XLALSimNoise(
	REAL8TimeSeries *s,		/**< [in/out] noise time series */
	size_t stride,			/**< [in] stride (samples) */
	REAL8FrequencySeries *psd,	/**< [in] power spectrum frequency series */
	gsl_rng *rng			/**< [in] GSL random number generator */
)
{
	REAL8Vector *overlap;
	size_t j;

	/* Use a default RNG if a NULL pointer was passed in */
	if (!rng)
		rng = gsl_rng_alloc(gsl_rng_default);

	/* make sure that the resolution of the frequency series is
	 * commensurate with the requested time series */
	if (s->data->length/2 + 1 != psd->data->length
			|| (size_t)floor(0.5 + 1.0/(s->deltaT * psd->deltaF)) != s->data->length)
		XLAL_ERROR(XLAL_EINVAL);

	/* stride cannot be longer than data length */
	if (stride > s->data->length)
		XLAL_ERROR(XLAL_EINVAL);

	if (stride == 0) { /* generate segment with no feathering */
		XLALSimNoiseSegment(s, psd, rng);
		return 0;
	} else if (stride == s->data->length) {
		/* will generate two independent noise realizations
		 * and feather them together with full overlap */
		XLALSimNoiseSegment(s, psd, rng);
		stride = 0;
	}

	overlap = XLALCreateREAL8Sequence(s->data->length - stride);

	/* copy overlap region between the old and the new data to temporary storage */
	memcpy(overlap->data, s->data->data + stride, overlap->length*sizeof(*overlap->data));
	
	/* generate the new data */
	XLALSimNoiseSegment(s, psd, rng);

	/* feather old data in overlap region with new data */
	for (j = 0; j < overlap->length; ++j) {
		double x = cos(LAL_PI*j/(2.0 * overlap->length));
		double y = sin(LAL_PI*j/(2.0 * overlap->length));
		s->data->data[j] = x*overlap->data[j] + y*s->data->data[j];
	}

	XLALDestroyREAL8Sequence(overlap);

	/* advance time */
	XLALGPSAdd(&s->epoch, stride * s->deltaT);
	return 0;
}


/*
 *
 * TEST CODE
 *
 */

#if 0

/* Example routine listed in documentation. */
void mkligodata(void)
{
	const double flow = 40.0; // 40 Hz low frequency cutoff
	const double duration = 16.0; // 16 second segments
	const double srate = 16384.0; // sampling rate in Hertz
	size_t length = duration * srate; // segment length
	size_t stride = length / 2; // stride between segments
	LIGOTimeGPS epoch = { 0, 0 };
	REAL8FrequencySeries *psd;
	REAL8TimeSeries *seg;
	gsl_rng *rng;
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	seg = XLALCreateREAL8TimeSeries("STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
	psd = XLALCreateREAL8FrequencySeries("LIGO SRD", &epoch, 0.0, 1.0/duration, &lalSecondUnit, length/2 + 1);
	XLALSimNoisePSD(psd, flow, XLALSimNoisePSDiLIGOSRD);
	XLALSimNoise(seg, 0, psd, rng); // first time to initialize
	while (1) { // infinite loop
		double t0 = XLALGPSGetREAL8(&seg->epoch);
		size_t j;
		for (j = 0; j < stride; ++j) // output first stride points
			printf("%.9f\t%e\n", t0 + j*seg->deltaT, seg->data->data[j]);
		XLALSimNoise(seg, stride, psd, rng); // make more data
	}
}


/*
 * Test routine that generates 1024 seconds of data in 8 second segments with
 * 50% overlap.  Also produced is the PSD for this data for comparison with the
 * PSD used to generate the data.
 */
int test_noise(void)
{
	const double recdur = 1024.0; // duration of the entire data record in seconds
	const double segdur = 8.0; // duration of a segment in seconds
	const double srate = 4096.0; // sample rate in Hertz
	const double flow = 9.0; // low frequency cutoff in Hertz
	LIGOTimeGPS epoch = {0, 0};
	REAL8TimeSeries *seg;
	REAL8TimeSeries *rec;
	REAL8FrequencySeries *psd;
	REAL8Window *window;
	REAL8FFTPlan *plan;
	gsl_rng *rng;
	double deltaT = 1.0 / srate;
	double deltaF = 1.0 / segdur;
	size_t reclen = recdur * srate;
	size_t seglen = segdur * srate;
	size_t stride = seglen / 2;
	size_t numseg = 1 + (reclen - seglen)/stride;
	size_t klow = flow / deltaF;
	size_t i, j, k;
	FILE *fp;

	if (reclen != (numseg - 1)*stride + seglen) {
		fprintf(stderr, "warning: numseg, seglen, reclen, and stride not commensurate\n");
		reclen = (numseg - 1)*stride + seglen;
	}

	rng = gsl_rng_alloc(gsl_rng_default);
	seg = XLALCreateREAL8TimeSeries("STRAIN", &epoch, 0.0, deltaT, &lalStrainUnit, seglen);
	rec = XLALCreateREAL8TimeSeries("STRAIN", &epoch, 0.0, deltaT, &lalStrainUnit, reclen);
	psd = XLALCreateREAL8FrequencySeries("PSD", &epoch, 0.0, deltaF, &lalSecondUnit, seglen/2 + 1);

	XLALSimNoisePSD(psd, 9.0, XLALSimNoisePSDaLIGOHighFrequency);

	fp = fopen("psd1.dat", "w");
	for (k = klow; k < psd->data->length - 1; ++k)
		fprintf(fp, "%e\t%e\n", k * psd->deltaF, sqrt(psd->data->data[k]));
	fclose(fp);

	for (i = 0; i < numseg; ++i) {
		/* first time, initialize; subsequent times, stride */
		XLALSimNoise(seg, i ? stride : 0, psd, rng);
		/* copy stride amount of data or seglen on last time */
		memcpy(rec->data->data + i*stride, seg->data->data, (i == numseg - 1 ? seglen : stride) * sizeof(*rec->data->data));
	}

	fp = fopen("noise.dat", "w");
	for (j = 0; j < rec->data->length; ++j) {
		epoch = rec->epoch;
		XLALGPSAdd(&epoch, j*rec->deltaT);
		fprintf(fp, "%.9f\t%e\n", XLALGPSGetREAL8(&epoch), rec->data->data[j]);
	}
	fclose(fp);

	plan = XLALCreateForwardREAL8FFTPlan(seglen, 0);
	window = XLALCreateHannREAL8Window(seglen);
	XLALREAL8AverageSpectrumWelch(psd, rec, seglen, stride, window, plan);

	fp = fopen("psd2.dat", "w");
	for (k = klow; k < psd->data->length - 1; ++k)
		fprintf(fp, "%e\t%e\n", k * psd->deltaF, sqrt(psd->data->data[k]));
	fclose(fp);

	XLALDestroyREAL8Window(window);
	XLALDestroyREAL8FFTPlan(plan);
	XLALDestroyREAL8FrequencySeries(psd);
	XLALDestroyREAL8TimeSeries(rec);
	XLALDestroyREAL8TimeSeries(seg);
	gsl_rng_free(rng);
	return 0;
}

int main(void)
{
	XLALSetErrorHandler(XLALAbortErrorHandler);
	gsl_rng_env_setup();
	// mkligodata();
	test_noise();
	LALCheckMemoryLeaks();
	return 0;
}

#endif
