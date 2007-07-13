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
#include <lal/LALDatatypes.h>
#include <lal/RealFFT.h>
#include <lal/Sequence.h>
#include <lal/TFTransform.h>


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
