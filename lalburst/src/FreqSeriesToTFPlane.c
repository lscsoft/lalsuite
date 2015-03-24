/*
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


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <complex.h>
#include <math.h>


#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>
#include <lal/XLALError.h>


/*
 * ============================================================================
 *
 *                      Time-Frequency Plane Projection
 *
 * ============================================================================
 */


/*
 * Multiply the data by the filter.  The check that the frequency
 * resolutions and units are compatible is omitted because it is implied by
 * the calling code.
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
	/* find bounds of common frequencies */
	const double flo = max(filterseries->f0, inputseries->f0);
	const double fhi = min(filterseries->f0 + filterseries->data->length * filterseries->deltaF, inputseries->f0 + inputseries->data->length * inputseries->deltaF);
	COMPLEX16 *output = outputseq->data + (int) round((flo - inputseries->f0) / inputseries->deltaF);
	COMPLEX16 *last = outputseq->data + (int) round((fhi - inputseries->f0) / inputseries->deltaF);
	const COMPLEX16 *input = inputseries->data->data + (int) round((flo - inputseries->f0) / inputseries->deltaF);
	const COMPLEX16 *filter = filterseries->data->data + (int) round((flo - filterseries->f0) / filterseries->deltaF);

	if(outputseq->length != inputseries->data->length)
		XLAL_ERROR_NULL(XLAL_EBADLEN);

	if(((unsigned) (output - outputseq->data) > outputseq->length) || (last - outputseq->data < 0))
		/* inputseries and filterseries don't intersect */
		memset(outputseq->data, 0, outputseq->length * sizeof(*outputseq->data));
	else {
		/* output = inputseries * conj(filter) */
		memset(outputseq->data, 0, (output - outputseq->data) * sizeof(*outputseq->data));
		for(; output < last; output++, input++, filter++)
			*output = *input * conj(*filter);
		memset(last, 0, (outputseq->length - (last - outputseq->data)) * sizeof(*outputseq->data));
	}

	return outputseq;
}


/**
 * Project a frequency series onto the comb of channel filters
 */
int XLALFreqSeriesToTFPlane(
	REAL8TimeFrequencyPlane *plane,
	const LALExcessPowerFilterBank *filter_bank,
	const COMPLEX16FrequencySeries *fseries,
	const REAL8FFTPlan *reverseplan
)
{
	COMPLEX16Sequence *fcorr;
	unsigned i;

	/* check input parameters */
	if((fmod(plane->deltaF, fseries->deltaF) != 0.0) ||
	   (fmod(plane->flow - fseries->f0, fseries->deltaF) != 0.0))
		XLAL_ERROR(XLAL_EINVAL);

	/* make sure the frequency series spans an appropriate band */
	if((plane->flow < fseries->f0) ||
	   (plane->flow + plane->channel_data->size2 * plane->deltaF > fseries->f0 + fseries->data->length * fseries->deltaF))
		XLAL_ERROR(XLAL_EDATA);

	/* create temporary vectors */
	fcorr = XLALCreateCOMPLEX16Sequence(fseries->data->length);
	if(!fcorr)
		XLAL_ERROR(XLAL_EFUNC);

#if 0
	/* diagnostic code to dump data for the \hat{s}_{k} histogram */
	{
	unsigned k;
	FILE *f = fopen("sk.dat", "a");
	for(k = plane->flow / fseries->deltaF; k < (plane->flow + plane->channel_data->size2 * plane->deltaF) / fseries->deltaF; k++)
		fprintf(f, "%g\n%g\n", fseries->data->data[k].re, fseries->data->data[k].im);
	fclose(f);
	}
#endif
#if 0
	/* diagnostic code to dump data for the \hat{s}_{k}
	 * \hat{s}^{*}_{k'} histogram */
	{
	unsigned k, dk;
	FILE *f = fopen("sksk.dat", "a");
	for(dk = 0; dk < 100; dk++) {
		double avg_r = 0;
		double avg_i = 0;
	for(k = plane->flow / fseries->deltaF; k + dk < (plane->flow + plane->channel_data->size2 * plane->deltaF) / fseries->deltaF; k++) {
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
	for(i = 0; i < plane->channel_data->size2; i++) {
		unsigned j;
		/* cross correlate the input data against the channel
		 * filter by taking their product in the frequency domain
		 * and then inverse transforming to the time domain to
		 * obtain an SNR time series.  Note that
		 * XLALREAL8ReverseFFT() omits the factor of 1 / (N Delta
		 * t) in the inverse transform. */
		apply_filter(fcorr, fseries, filter_bank->basis_filters[i].fseries);
		if(XLALREAL8ReverseFFT(plane->channel_buffer, fcorr, reverseplan)) {
			XLALDestroyCOMPLEX16Sequence(fcorr);
			XLAL_ERROR(XLAL_EFUNC);
		}
		/* interleave the result into the channel_data array */
		for(j = 0; j < plane->channel_buffer->length; j++)
			gsl_matrix_set(plane->channel_data, j, i, plane->channel_buffer->data[j]);
	}

	/* clean up */
	XLALDestroyCOMPLEX16Sequence(fcorr);

	/* set the name and epoch of the TF plane */
	strncpy(plane->name, fseries->name, LALNameLength);
	plane->epoch = fseries->epoch;

	/* success */
	return 0;
}
