#include <math.h>
#include <string.h>

#include <gsl/gsl_rng.h>

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/LALSimBurst.h>

/**
 * \brief Generates a burst injection with a time-frequency structure specified
 * in an image array.
 *
 * \note The image array data is ordered so that the first point corresponds
 * to the upper-left point of the image.  The number of rows in the image
 * is given in NUMROWS = image->dimLength->data[0] and the number of columns in
 * the image is given in NUMCOLUMNS = image->dimLength->data[1].  The data is
 * packed in row-major order so that offset = row * NUMCOLUMNS + column.
 *
 * \note The time-frequency volume of each pixel, dt * df, must be at least
 * 2/pi.
 */
int XLALSimBurstImg(
	REAL8TimeSeries **hplus,	/**< plus-polarization waveform [returned] */
	REAL8TimeSeries **hcross, 	/**< cross-polarization waveform [returned] */
	REAL8Array *image,		/**< image array */
	double dt,			/**< pixel time duration (s) */
	double df,			/**< pixel frequency bandwidth (Hz) */
	double fstart,			/**< start frequency of image (Hz) */
	double hrss,			/**< root-sum-squared value of the waveform (s) */
	double deltaT,			/**< sampleing interval (s) */
	gsl_rng *rng			/**< random number generator */
)
{
	LIGOTimeGPS epoch;
	size_t row, col, nrow, ncol, pad, length;
	double fac, hss;
	size_t j;

	/* sanity check on input values */
	XLAL_CHECK(dt * df > LAL_2_PI, XLAL_EINVAL, "Time-frequency volume dt*df must be greater than 2/pi");
	XLAL_CHECK(image->dimLength->length == 2, XLAL_EINVAL, "Requires a 2-dimensional array");

	nrow = image->dimLength->data[0];
	ncol = image->dimLength->data[1];
	
	/* make a time series that has some padding at start and end */
	pad = floor(15.0 * dt / deltaT);
	length = floor((ncol - 1.) * dt / deltaT) + 2 * pad;
	XLALGPSSetREAL8(&epoch, -1. * pad * deltaT);
	*hplus = XLALCreateREAL8TimeSeries("Image +", &epoch, 0.0, deltaT, &lalStrainUnit, length);
	*hcross = XLALCreateREAL8TimeSeries("Image x", &epoch, 0.0, deltaT, &lalStrainUnit, length);
	if (!*hplus || !*hcross)
		XLAL_ERROR(XLAL_EFUNC);
	memset((*hplus)->data->data, 0, length * sizeof(*(*hplus)->data->data));
	memset((*hcross)->data->data, 0, length * sizeof(*(*hcross)->data->data));

	for (row = 0; row < nrow; ++row) {
		for (col = 0; col < ncol; ++col) {
			REAL8TimeSeries *hp;
			REAL8TimeSeries *hc;
			double t = floor(col * dt / deltaT) * deltaT;
			double f = fstart + (nrow - row) * df;
			double Y = image->data[row * ncol + col];
			if (XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(&hp, &hc, dt, f, df, 0., Y, deltaT, rng) < 0) {
				XLALDestroyREAL8TimeSeries(hc);
				XLALDestroyREAL8TimeSeries(hp);
				XLAL_ERROR(XLAL_EFUNC);
			}
			XLALGPSAdd(&hp->epoch, t);
			XLALGPSAdd(&hc->epoch, t);
			XLALAddREAL8TimeSeries(*hplus, hp);
			XLALAddREAL8TimeSeries(*hcross, hc);
			XLALDestroyREAL8TimeSeries(hc);
			XLALDestroyREAL8TimeSeries(hp);
		}
	}

	/* normalize to desired hrss value */
	hss = 0;
	for (j = 0; j < (*hplus)->data->length; ++j) {
		double h = (*hplus)->data->data[j];
		hss += h * h * (*hplus)->deltaT;
	}
	for (j = 0; j < (*hcross)->data->length; ++j) {
		double h = (*hcross)->data->data[j];
		hss += h * h * (*hcross)->deltaT;
	}
	fac = hrss / sqrt(hss);
	for (j = 0; j < (*hplus)->data->length; ++j)
		(*hplus)->data->data[j] *= fac;
	for (j = 0; j < (*hcross)->data->length; ++j)
		(*hcross)->data->data[j] *= fac;

	return 0;
}
