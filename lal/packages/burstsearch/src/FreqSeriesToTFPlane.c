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


#include <math.h>
/* FIXME: decide what to do about this
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
*/


#include <lal/LALRCSID.h>


NRCSID(FREQSERIESTOTFPLANEC, "$Id$");


#include <lal/Date.h>
#include <lal/LALDatatypes.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/Thresholds.h>
#include <lal/XLALError.h>




/*
 * Multiply the data by the filter.  The check that the frequency
 * resolutions are compatible is omitted because it is implied by the
 * calling code.
 */


static double min(double a, double b)
{
	return a < b ? a : b;
}


static double max(double a, double b)
{
	return a > b ? a : b;
}


static COMPLEX16Sequence *apply_filter(
	COMPLEX16Sequence *outputseq,
	const COMPLEX16FrequencySeries *inputseries,
	const COMPLEX16FrequencySeries *filterseries
)
{
	static const char func[] = "apply_filter";
	/* find bounds of common frequencies */
	const double flo = max(filterseries->f0, inputseries->f0);
	const double fhi = min(filterseries->f0 + filterseries->data->length * filterseries->deltaF, inputseries->f0 + inputseries->data->length * inputseries->deltaF);
	COMPLEX16 *output = outputseq->data + (int) ((flo - inputseries->f0) / inputseries->deltaF);
	COMPLEX16 *last = outputseq->data + (int) ((fhi - inputseries->f0) / inputseries->deltaF);
	const COMPLEX16 *input = inputseries->data->data + (int) ((flo - inputseries->f0) / inputseries->deltaF);
	const COMPLEX16 *filter = filterseries->data->data + (int) ((flo - filterseries->f0) / filterseries->deltaF);

	if(outputseq->length != inputseries->data->length)
		XLAL_ERROR_NULL(func, XLAL_EBADLEN);

	/* zero the product vector */
	memset(outputseq->data, 0, outputseq->length * sizeof(*outputseq->data));

	if(fhi < flo)
		/* no op */
		return outputseq;

	/* output = inputseries * conj(filter) */
	for(; output < last; output++, input++, filter++) {
		output->re = input->re * filter->re + input->im * filter->im;
		output->im = input->im * filter->re - input->re * filter->im;
	}

	return outputseq;
}


/******** <lalVerbatim file="FreqSeriesToTFPlaneCP"> ********/
int XLALFreqSeriesToTFPlane(
	REAL8TimeFrequencyPlane *plane,
	const COMPLEX16FrequencySeries *fseries,
	const REAL8FFTPlan *reverseplan
)
/******** </lalVerbatim> ********/
{
	static const char func[] = "XLALFreqSeriesToTFPlane";
	COMPLEX16Sequence *fcorr;
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
	fcorr = XLALCreateCOMPLEX16Sequence(fseries->data->length);
	if(!fcorr) {
		XLALDestroyCOMPLEX16Sequence(fcorr);
		XLAL_ERROR(func, XLAL_EFUNC);
	}

#if 0
	/* diagnostic code to dump data for the s_{k} histogram */
	{
	unsigned k;
	FILE *f = fopen("sk.dat", "a");
	for(k = plane->flow / fseries->deltaF; k < (plane->flow + plane->channels * plane->deltaF) / fseries->deltaF; k++)
		fprintf(f, "%g\n%g\n", fseries->data->data[k].re, fseries->data->data[k].im);
	fclose(f);
	}
#endif
#if 0
	/* diagnostic code to dump data for the s_{k} s^{*}_{k'} histogram
	 * */
	{
	unsigned k, dk;
	FILE *f = fopen("sksk.dat", "a");
	for(dk = 0; dk < 100; dk++) {
		double avg_r = 0;
		double avg_i = 0;
	for(k = plane->flow / fseries->deltaF; k + dk < (plane->flow + plane->channels * plane->deltaF) / fseries->deltaF; k++) {
		double dr = fseries->data->data[k].re;
		double di = fseries->data->data[k].im;
		double dkr = fseries->data->data[k + dk].re;
		double dki = fseries->data->data[k + dk].im;
		avg_r += dr * dkr + di * dki;
		avg_i += di * dkr - dr * dki;
	}
		avg_r /= k - plane->flow / fseries->deltaF;
		avg_i /= k - plane->flow / fseries->deltaF;
		fprintf(f, "%d %g %g\n", dk, avg_r, avg_i);
	}
	fclose(f);
	}
#endif

	/* loop over the time-frequency plane's channels */
	for(i = 0; i < plane->channels; i++) {
		/* cross correlate the input data against the channel
		 * filter by taking their product in the frequency domain
		 * and then inverse transforming to the time domain to
		 * obtain an SNR time series.  Note that
		 * XLALREAL4ReverseFFT() omits the factor of 1 / (N Delta
		 * t) in the inverse transform. */
		apply_filter(fcorr, fseries, plane->filter[i]);
		if(XLALREAL8ReverseFFT(plane->channel[i], fcorr, reverseplan)) {
			XLALDestroyCOMPLEX16Sequence(fcorr);
			XLAL_ERROR(func, XLAL_EFUNC);
		}
	}

	/* clean up */
	XLALDestroyCOMPLEX16Sequence(fcorr);

	/* set the name and epoch of the TF plane */
	strncpy(plane->name, fseries->name, LALNameLength);
	plane->epoch = fseries->epoch;

	/* success */
	return 0;
}


/*
 * ============================================================================
 *
 *                                   Output
 *
 * ============================================================================
 */


/*
 * Convert time-frequency tile info to a SnglBurstTable row.
 */


static SnglBurstTable *XLALTFTileToBurstEvent(
	const REAL8TimeFrequencyPlane *plane,
	unsigned tile_start,
	unsigned tile_length,
	unsigned tile_channel,
	unsigned tile_channels,
	double h_rss,
	double E,
	double d,
	double confidence
)
{
	static const char func[] = "XLALTFTileToBurstEvent";
	SnglBurstTable *event = XLALCalloc(1, sizeof(*event));

	if(!event)
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);

	event->next = NULL;
	strncpy(event->ifo, plane->name, 2);
	event->ifo[2] = '\0';
	strncpy(event->search, "excesspower", LIGOMETA_SEARCH_MAX);
	event->search[LIGOMETA_SEARCH_MAX] = '\0';
	strncpy(event->channel, plane->name, LIGOMETA_CHANNEL_MAX);
	event->channel[LIGOMETA_CHANNEL_MAX] = '\0';
	event->start_time = plane->epoch; 
	XLALGPSAdd(&event->start_time, tile_start * plane->deltaT);
	event->duration = tile_length * plane->deltaT;
	event->peak_time = event->start_time;
	XLALGPSAdd(&event->peak_time, 0.5 * event->duration);
	event->bandwidth = tile_channels * plane->deltaF;
	event->central_freq = plane->flow + tile_channel * plane->deltaF + (0.5 * event->bandwidth);
	/* FIXME: put h_rss into the "hrss" column */
	event->amplitude = h_rss;
	event->snr = E / d - 1;
	/* -ln P(event | stationary Gaussian white noise) */
	event->confidence = confidence;
	/* for safety */
	event->string_cluster_t = XLAL_REAL4_FAIL_NAN;
	/* will be set correctly later */
	event->event_id = 0;

	return(event);
}


/*
 * ============================================================================
 *
 *                         Time-Domain Decorrelation
 *
 * ============================================================================
 */


/* FIXME:  decide what to do about this */
#if 0
static REAL8Sequence *impulse_response(unsigned N, COMPLEX16FrequencySeries **filter, unsigned channel, unsigned channels, REAL8FFTPlan *reverseplan)
{
	COMPLEX16Sequence *f = XLALCreateCOMPLEX16Sequence(N / 2 + 1);
	REAL8Sequence *impulse = XLALCreateREAL8Sequence(N);
	unsigned i, k;

	if(!f || !impulse) {
		XLALDestroyCOMPLEX16Sequence(f);
		XLALDestroyREAL8Sequence(impulse);
		return NULL;
	}

	/* sum the channel filters */
	memset(f->data, 0, f->length * sizeof(*f->data));
	for(i = channel; i < channel + channels; i++) {
		unsigned k0 = filter[i]->f0 / filter[i]->deltaF;
		for(k = 0; k < filter[i]->data->length; k++) {
			f->data[k0 + k].re += filter[i]->data->data[k].re;
			f->data[k0 + k].im += filter[i]->data->data[k].im;
		}
	}

	/* inverse transform */
	if(XLALREAL8ReverseFFT(impulse, f, reverseplan)) {
		XLALDestroyREAL8Sequence(impulse);
		impulse = NULL;
	}

	/* done */
	XLALDestroyCOMPLEX16Sequence(f);
	return impulse;
}


static int impulse_response_matrix_lu(gsl_matrix *m, gsl_permutation *p, unsigned stride, unsigned N, COMPLEX16FrequencySeries **filter, unsigned channel, unsigned channels, REAL8FFTPlan *reverseplan)
{
	REAL8Sequence *impulse = impulse_response(N, filter, channel, channels, reverseplan);
	double norm;
	int ignored[8192];
	int i, j;

	if(!impulse) {
		XLALDestroyREAL8Sequence(impulse);
		return -1;
	}

	/* calculate normalization factor:  sum of squares across each row
	 * of convolution matrix must be 1 */
	norm = 0;
	for(j = 0; (j < (int) m->size2) && (j * stride < impulse->length); j++)
		norm += impulse->data[j * stride] * impulse->data[j * stride];
	norm = sqrt(norm);

	for(i = 0; i < (int) m->size1; i++)
		for(j = 0; j < (int) m->size2; j++) {
			unsigned k = abs(i - j) * stride;
			gsl_matrix_set(m, i, j, k < impulse->length ? impulse->data[k] / norm : 0);
		}

	if(gsl_linalg_LU_decomp(m, p, ignored))
		XLALPrintError("gsl_linalg_LU_decomp() failed!\n");

	XLALDestroyREAL8Sequence(impulse);
	return 0;
}
#endif



/*
 * ============================================================================
 *
 *                               Tile Analysis
 *
 * ============================================================================
 */


/******** <lalVerbatim file="ComputeExcessPowerCP"> ********/
SnglBurstTable *XLALComputeExcessPower(
	const REAL8TimeFrequencyPlane *plane,
	SnglBurstTable *head,
	double confidence_threshold,
	REAL8FFTPlan *reverseplan
)
/******** </lalVerbatim> ********/
{
	static const char func[] = "XLALComputeExcessPower";
	unsigned start;
	unsigned length;
	unsigned end;
	unsigned channel;
	unsigned channels;
	unsigned channel_end;
	double h_rss;
	double confidence;
	/* FIXME: decide what to do about this */
#if 0
	gsl_vector *samples = NULL;
	gsl_matrix *deconv = NULL;
	gsl_permutation *deconvp = NULL;
#endif

	for(channels = plane->tiles.min_channels; channels <= plane->tiles.max_channels; channels *= 2) {
	for(channel_end = (channel = 0) + channels; channel_end <= plane->channels; channel_end = (channel += channels / plane->tiles.inv_fractional_stride) + channels) {
		/* the mean square of the "virtual channel" */
		const double sum_mean_square = channels * plane->deltaF / plane->fseries_deltaF + XLALREAL8SequenceSum(plane->twice_channel_overlap, channel, channels - 1);
		/* the mean square of the "uwsum" quantity computed below,
		 * which is proportional to an approximation of the
		 * unwhitened time series. */
		double uwsum_mean_square;
		/* true unwhitened mean square for this channel.  the ratio
		 * of this to uwsum_mean_square is the correction factor to
		 * be applied to uwsum^{2} to convert it to an
		 * approximation of the square of the unwhitened channel */
		const double unwhitened_mean_square = XLALREAL8SequenceSumSquares(plane->unwhitened_rms, channel, channels) + XLALREAL8SequenceSum(plane->unwhitened_cross, channel, channels - 1);
		unsigned c;

		/* compute uwsum_mean_square */
		uwsum_mean_square = XLALREAL8SequenceSumSquares(plane->unwhitened_rms, channel, channels);
		for(c = channel; c < channel_end - 1; c++)
			uwsum_mean_square += plane->twice_channel_overlap->data[c] * plane->unwhitened_rms->data[c] * plane->unwhitened_rms->data[c + 1] * plane->fseries_deltaF / plane->deltaF;

	/* start with at least 2 degrees of freedom */
	for(length = 2 / (channels * plane->tiles.dof_per_pixel); length <= plane->tiles.max_length; length *= 2) {
		const double tile_dof = (length * channels) * plane->tiles.dof_per_pixel;
		/* number of degrees of freedom in tile = number of
		 * "virtual pixels" in tile.  compute distance between
		 * tile's "virtual pixels" */
		const unsigned stride = length / tile_dof;

		/* LU factorization of channel's normalized time-domain
		 * impulse response convolution matrix */
		/* FIXME: decide what to do about this */
#if 0
		samples = gsl_vector_alloc(tile_dof);
		deconv = gsl_matrix_alloc(tile_dof, tile_dof);
		deconvp = gsl_permutation_alloc(tile_dof);
		impulse_response_matrix_lu(deconv, deconvp, stride, plane->window->data->length, plane->filter, channel, channels, reverseplan);
#endif

	for(end = (start = plane->tiles.tiling_start) + length; end <= plane->tiles.tiling_end; end = (start += length / plane->tiles.inv_fractional_stride) + length) {
		double sumsquares = 0;
		double uwsumsquares = 0;
		unsigned t;

		/* compute sum of squares, and unwhitened sum of squares */
		for(t = start + stride / 2; t < end; t += stride) {
			double sum = 0;
			double uwsum = 0;

			/* compute snr and dewhitened amplitude for this
			 * "virtual pixel". */
			for(c = channel; c < channel_end; c++) {
				sum += plane->channel[c]->data[t];
				uwsum += plane->channel[c]->data[t] * plane->unwhitened_rms->data[c] * sqrt(plane->fseries_deltaF / plane->deltaF);
			}
			/* FIXME: deicide what to do about this */
#if 0
			gsl_vector_set(samples, (t - start - stride / 2) / stride, sum);
#endif

			sumsquares += sum * sum;
			uwsumsquares += uwsum * uwsum;
		}
		/* normalization to give each pixel a mean square of 1.0 */
		sumsquares /= sum_mean_square;

		/* FIXME: decide what to do about this */
#if 0
		gsl_linalg_LU_svx(deconv, deconvp, samples);
#endif

		/* compute statistical confidence, see if it's above
		 * threshold */
		/* FIXME:  the 0.65 is an empirically determined
		 * degree-of-freedom fudge factor.  figure out what its
		 * origin is, and account for it correctly.  it's most
		 * likely due to the time-frequency plane pixels not being
		 * independent of one another as a consequence of a
		 * non-zero inner product of the time-domain impulse
		 * response of the channel filter for adjacent pixels */
		confidence = -XLALlnOneMinusChisqCdf(sumsquares * .65, tile_dof * .65);
		if(XLALIsREAL8FailNaN(confidence))
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		if(confidence >= confidence_threshold) {
			SnglBurstTable *oldhead;

			/* normalization to give each unwhitened pixel a
			 * mean square of 1.0 */
			uwsumsquares /= uwsum_mean_square;

			/* compute h_rss */
			h_rss = sqrt((uwsumsquares - tile_dof) * unwhitened_mean_square * stride * plane->deltaT);

			/* add new event to head of linked list */
			oldhead = head;
			head = XLALTFTileToBurstEvent(plane, start, length, channel, channels, h_rss, sumsquares, tile_dof, confidence);
			if(!head)
				XLAL_ERROR_NULL(func, XLAL_EFUNC);
			head->next = oldhead;

		}
	}
		/* FIXME: decide what to do about this */
#if 0
		gsl_vector_free(samples);
		gsl_matrix_free(deconv);
		gsl_permutation_free(deconvp);
		samples = NULL;
		deconv = NULL;
		deconvp = NULL;
#endif
	}
	}
	}

	/* success */
	return(head);
}
