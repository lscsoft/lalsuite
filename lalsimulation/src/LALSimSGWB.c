/*
*  Copyright (C) 2011 Jolien Creighton
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
#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <lal/LALConstants.h>
#include <lal/LALDetectors.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/Sequence.h>
#include <lal/TimeFreqFFT.h>
#include <lal/Units.h>
#include <lal/LALSimSGWB.h>


/**
 * Creates a frequency series that contains a flat SGWB spectrum with the
 * specified power Omega0 above some low frequency cutoff flow.
 */
REAL8FrequencySeries *XLALSimSGWBOmegaGWFlatSpectrum(
	double Omega0,	/**< [in] sgwb spectrum power (dimensionless) */
	double flow,	/**< [in] low frequncy cutoff of SGWB spectrum (Hz) */
	double deltaF,	/**< [in] frequency bin width (Hz) */
	size_t length	/**< [in] number of frequency bins */
)
{
	REAL8FrequencySeries *OmegaGW;
	LIGOTimeGPS epoch = {0, 0};
	size_t klow = flow / deltaF;
	size_t k;
	OmegaGW = XLALCreateREAL8FrequencySeries("OmegaGW", &epoch, 0.0, deltaF, &lalDimensionlessUnit, length);
	/* zero DC component */
	OmegaGW->data->data[0] = 0.0;
	/* zero up to low frequency cutoff */
	for (k = 1; k < klow; ++k)
		OmegaGW->data->data[k] = 0.0;
	/* set remaining components */
	for (; k < length - 1; ++k)
		OmegaGW->data->data[k] = Omega0;
	/* zero Nyquist component */
	OmegaGW->data->data[length - 1] = 0.0;
	return OmegaGW;
}


/**
 * Creates a frequency series that contains a power law SGWB spectrum with the
 * specified power Omegaref at reference frequency fref and specified power law
 * power alpha above some low frequency cutoff flow.
 */
REAL8FrequencySeries *XLALSimSGWBOmegaGWPowerLawSpectrum(
	double Omegaref,	/**< [in] sgwb spectrum power at reference frequency (dimensionless) */
	double alpha,		/**< [in] sgwb spectrum power law power */
	double fref,		/**< [in] reference frequency (Hz) */
	double flow,		/**< [in] low frequncy cutoff of SGWB spectrum (Hz) */
	double deltaF,		/**< [in] frequency bin width (Hz) */
	size_t length		/**< [in] number of frequency bins */
)
{
	REAL8FrequencySeries *OmegaGW;
	LIGOTimeGPS epoch = {0, 0};
	size_t klow = flow / deltaF;
	size_t k;
	OmegaGW = XLALCreateREAL8FrequencySeries("OmegaGW", &epoch, 0.0, deltaF, &lalDimensionlessUnit, length);
	/* zero DC component */
	OmegaGW->data->data[0] = 0.0;
	/* zero up to low frequency cutoff */
	for (k = 1; k < klow; ++k)
		OmegaGW->data->data[k] = 0.0;
	/* set remaining components */
	for (; k < length - 1; ++k)
		OmegaGW->data->data[k] = Omegaref * pow(k * deltaF / fref, alpha);
	/* zero Nyquist component */
	OmegaGW->data->data[length - 1] = 0.0;
	return OmegaGW;
}


/* 
 * This routine generates a single segment of data.  Note that this segment is
 * generated in the frequency domain and is inverse Fourier transformed into
 * the time domain; consequently the data is periodic in the time domain.
 */
static int XLALSimSGWBSegment(REAL8TimeSeries **h, const LALDetector *detectors, size_t numDetectors, const REAL8FrequencySeries *OmegaGW, double H0, gsl_rng *rng)
{
#	define CLEANUP_AND_RETURN(errnum) do { \
		if (htilde) for (i = 0; i < numDetectors; ++i) XLALDestroyCOMPLEX16FrequencySeries(htilde[i]); \
		XLALFree(htilde); XLALDestroyREAL8FFTPlan(plan); gsl_matrix_free(R); \
		if (errnum) XLAL_ERROR(errnum); else return 0; \
		} while (0)
	REAL8FFTPlan *plan = NULL;
	COMPLEX16FrequencySeries **htilde = NULL;
	gsl_matrix *R = NULL;
	LIGOTimeGPS epoch;
	double psdfac;
	double deltaF;
	size_t length;
	size_t i, j, k;

	epoch = h[0]->epoch;
	length = h[0]->data->length;
	deltaF = 1.0 / (length * h[0]->deltaT);
	psdfac = 0.3 * pow(H0 / LAL_PI, 2.0);

	R = gsl_matrix_alloc(numDetectors, numDetectors);
	if (! R)
		CLEANUP_AND_RETURN(XLAL_ENOMEM);

	plan = XLALCreateReverseREAL8FFTPlan(length, 0);
	if (! plan)
		CLEANUP_AND_RETURN(XLAL_EFUNC);

	/* allocate frequency series for the various detector strains */
	htilde = LALCalloc(numDetectors, sizeof(*htilde));
	if (! htilde)
		CLEANUP_AND_RETURN(XLAL_ENOMEM);
	for (i = 0; i < numDetectors; ++i) {
		htilde[i] = XLALCreateCOMPLEX16FrequencySeries(h[i]->name, &epoch, 0.0, deltaF, &lalSecondUnit, length/2 + 1);
		if (! htilde[i])
			CLEANUP_AND_RETURN(XLAL_EFUNC);

		/* correct units */
		XLALUnitMultiply(&htilde[i]->sampleUnits, &htilde[i]->sampleUnits, &h[i]->sampleUnits);

		/* set data to zero */
		memset(htilde[i]->data->data, 0, htilde[i]->data->length * sizeof(*htilde[i]->data->data));
	}

	/* compute frequencies (excluding DC and Nyquist) */
	for (k = 1; k < length/2; ++k) {
		double f = k * deltaF;
		double sigma = 0.5 * sqrt(psdfac * OmegaGW->data->data[k] * pow(f, -3.0) / deltaF);

		/* construct correlation matrix at this frequency */
		/* diagonal elements of correlation matrix are unity */
		gsl_matrix_set_identity(R);
		/* now do the off-diagonal elements */
		for (i = 0; i < numDetectors; ++i)
			for (j = i + 1; j < numDetectors; ++j) {
				double Rij = XLALSimSGWBOverlapReductionFunction(f, &detectors[i], &detectors[j]);
				/* if the two sites are the same, the overlap reduciton
				 * function will be unity, but this will cause problems
				 * for the cholesky decomposition; a hack is to make it
				 * unity only to single precision */
				if (fabs(Rij - 1.0) < LAL_REAL4_EPS)
					Rij = 1.0 - LAL_REAL4_EPS;

				gsl_matrix_set(R, i, j, Rij);
				gsl_matrix_set(R, j, i, Rij); /* it is symmetric */
			}

		/* perform Cholesky decomposition */
		gsl_linalg_cholesky_decomp(R);

		/* generate numDetector random numbers (both re and im parts) and use
 		 * lower-diagonal part of Cholesky decomposition to create correlations */
		for (j = 0; j < numDetectors; ++j) {
			double re = gsl_ran_gaussian_ziggurat(rng, sigma);
			double im = gsl_ran_gaussian_ziggurat(rng, sigma);
			for (i = j; i < numDetectors; ++i) {
				htilde[i]->data->data[k] += gsl_matrix_get(R, i, j) * re;
				htilde[i]->data->data[k] += I * gsl_matrix_get(R, i, j) * im;
			}
		}
	}

	/* now go back to the time domain */
	for (i = 0; i < numDetectors; ++i)
		XLALREAL8FreqTimeFFT(h[i], htilde[i], plan);

	/* normal exit */
	CLEANUP_AND_RETURN(0);
#	undef CLEANUP_AND_RETURN
}


/**
 * Routine that may be used to generate sequential segments of stochastic
 * background gravitational wave signals for a network of detectors with a
 * specified stride from one segment to the next.
 *
 * The spectrum is specified by the frequency series OmegaGW.
 *
 * Calling instructions: for the first call, set stride = 0; subsequent calls
 * should pass the same time series and have non-zero stride.  This routine
 * will advance the time series by an amount given by the stride and will
 * generate new data so that the data is continuous from one segment to the
 * next.  For example: the following routine will output a continuous stream of
 * stochastic background signals with a "flat" (OmegaGW = const) spectrum for
 * the HLV network.
 *
 * \code
 * #include <stdio.h>
 * #include <gsl/gsl_rng.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/LALDetectors.h>
 * #include <lal/FrequencySeries.h>
 * #include <lal/TimeSeries.h>
 * #include <lal/Units.h>
 * #include <lal/LALSimSGWB.h>
 * int mksgwbdata(void)
 * {
 * 	const double flow = 40.0; // 40 Hz low frequency cutoff
 * 	const double duration = 16.0; // 16 second segments
 * 	const double srate = 16384.0; // sampling rate in Hertz
 * 	const double Omega0 = 1e-6; // fraction of critical energy in GWs
 * 	const double H0 = 0.72 * LAL_H0FAC_SI; // Hubble's constant in seconds
 * 	const LALDetector H1 = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; // Hanford
 * 	const LALDetector L1 = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; // Livingston
 * 	const LALDetector V1 = lalCachedDetectors[LAL_VIRGO_DETECTOR]; // Virgo
 * 	LALDetector detectors[3] = {H1, L1, V1}; // the network of detectors
 * 	size_t length = duration * srate; // segment length
 * 	size_t stride = length / 2; // stride between segments
 * 	LIGOTimeGPS epoch = { 0, 0 };
 * 	REAL8FrequencySeries *OmegaGW; // the spectrum of the SGWB
 * 	REAL8TimeSeries *h[3]; // the strain induced in the network of detectors
 * 	gsl_rng *rng;
 * 	gsl_rng_env_setup();
 * 	rng = gsl_rng_alloc(gsl_rng_default);
 * 	h[0] = XLALCreateREAL8TimeSeries("H1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
 * 	h[1] = XLALCreateREAL8TimeSeries("L1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
 * 	h[2] = XLALCreateREAL8TimeSeries("V1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
	OmegaGW = XLALSimSGWBOmegaGWFlatSpectrum(Omega0, flow, deltaF, seglen/2 + 1);
 * 	XLALSimSGWB(h, detectors, 3, 0, Omega0, flow, H0, rng); // first time to initialize
 * 	while (1) { // infinite loop
 * 		size_t j;
 * 		for (j = 0; j < stride; ++j) // output first stride points
 * 			printf("%.9f\t%e\t%e\t%e\n", XLALGPSGetREAL8(&h[0]->epoch) + j / srate, h[0]->data->data[j], h[1]->data->data[j], h[2]->data->data[j]);
 * 		XLALSimSGWB(h, detectors, 3, stride, Omega0, flow, H0, rng); // make more data
 * 	}
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
int XLALSimSGWB(
	REAL8TimeSeries **h,			/**< [in/out] array of sgwb timeseries for detector network */
	const LALDetector *detectors,		/**< [in] array of detectors in network */
	size_t numDetectors,			/**< [in] number of detectors in network */
	size_t stride,				/**< [in] stride (samples) */
	const REAL8FrequencySeries *OmegaGW,	/**< [in] sgwb spectrum frequeny series */
	double H0,				/**< [in] Hubble's constant (s) */
	gsl_rng *rng				/**< [in] GSL random number generator */
)
{
#	define CLEANUP_AND_RETURN(errnum) do { \
		if (overlap) for (i = 0; i < numDetectors; ++i) XLALDestroyREAL8Sequence(overlap[i]); \
		XLALFree(overlap); \
		if (errnum) XLAL_ERROR(errnum); else return 0; \
		} while (0)
	REAL8Vector **overlap = NULL;
	LIGOTimeGPS epoch;
	size_t length;
	double deltaT;
	size_t i, j;

	length = h[0]->data->length;
	deltaT = h[0]->deltaT;
	epoch = h[0]->epoch;

	/* make sure all the lengths and other metadata are the same */
	for (i = 1; i < numDetectors; ++i)
		if (h[i]->data->length != length
				|| fabs(h[i]->deltaT - deltaT) > LAL_REAL8_EPS
				|| XLALGPSCmp(&epoch, &h[i]->epoch))
			XLAL_ERROR(XLAL_EINVAL);

	/* make sure that the resolution of the frequency series is
	 * commensurate with the requested time series */
	if (length/2 + 1 != OmegaGW->data->length
			|| (size_t)floor(0.5 + 1.0/(deltaT * OmegaGW->deltaF)) != length)
		XLAL_ERROR(XLAL_EINVAL);

	/* stride cannot be longer than data length */
	if (stride > length)
		XLAL_ERROR(XLAL_EINVAL);

	if (stride == 0) { /* generate segment with no feathering */
		XLALSimSGWBSegment(h, detectors, numDetectors, OmegaGW, H0, rng);
		return 0;
	} else if (stride == length) {
		/* will generate two independent noise realizations
		 * and feather them together with full overlap */
		XLALSimSGWBSegment(h, detectors, numDetectors, OmegaGW, H0, rng);
		stride = 0;
	}

	overlap = LALCalloc(numDetectors, sizeof(*overlap));
	if (! overlap)
		XLAL_ERROR(XLAL_ENOMEM);
	for (i = 0; i < numDetectors; ++i) {
		overlap[i] = XLALCreateREAL8Sequence(length - stride);
		if (! overlap[i])
			CLEANUP_AND_RETURN(XLAL_EFUNC);
		/* copy overlap region between the old and the new data to temporary storage */
		memcpy(overlap[i]->data, h[i]->data->data + stride, overlap[i]->length*sizeof(*overlap[i]->data));
	}

	if (XLALSimSGWBSegment(h, detectors, numDetectors, OmegaGW, H0, rng))
		CLEANUP_AND_RETURN(XLAL_EFUNC);

	/* feather old data in overlap region with new data */
	for (j = 0; j < length - stride; ++j) {
		double x = cos(LAL_PI*j/(2.0 * (length - stride)));
		double y = sin(LAL_PI*j/(2.0 * (length - stride)));
		for (i = 0; i < numDetectors; ++i)
			h[i]->data->data[j] = x*overlap[i]->data[j] + y*h[i]->data->data[j];
	}

	/* advance time */
	for (i = 0; i < numDetectors; ++i)
		XLALGPSAdd(&h[i]->epoch, stride * deltaT);

	/* success */
	CLEANUP_AND_RETURN(0);
#	undef CLEANUP_AND_RETURN
}


/**
 * Routine that may be used to generate sequential segments of stochastic
 * background gravitational wave signals for a network of detectors with a
 * specified stride from one segment to the next.
 *
 * The spectrum is flat for frequencies above flow with power given by Omega0,
 * and zero for frequencies below the low frequency cutoff flow.
 *
 * Calling instructions: for the first call, set stride = 0; subsequent calls
 * should pass the same time series and have non-zero stride.  This routine
 * will advance the time series by an amount given by the stride and will
 * generate new data so that the data is continuous from one segment to the
 * next.  For example: the following routine will output a continuous stream of
 * stochastic background signals with a "flat" (OmegaGW = const) spectrum for
 * the HLV network.
 *
 * \code
 * #include <stdio.h>
 * #include <gsl/gsl_rng.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/LALDetectors.h>
 * #include <lal/FrequencySeries.h>
 * #include <lal/TimeSeries.h>
 * #include <lal/Units.h>
 * #include <lal/LALSimSGWB.h>
 * int mkgwbdata_flat(void)
 * {
 * 	const double flow = 40.0; // 40 Hz low frequency cutoff
 * 	const double duration = 16.0; // 16 second segments
 * 	const double srate = 16384.0; // sampling rate in Hertz
 * 	const double Omega0 = 1e-6; // fraction of critical energy in GWs
 * 	const double H0 = 0.72 * LAL_H0FAC_SI; // Hubble's constant in seconds
 * 	const LALDetector H1 = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; // Hanford
 * 	const LALDetector L1 = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; // Livingston
 * 	const LALDetector V1 = lalCachedDetectors[LAL_VIRGO_DETECTOR]; // Virgo
 * 	LALDetector detectors[3] = {H1, L1, V1}; // the network of detectors
 * 	size_t length = duration * srate; // segment length
 * 	size_t stride = length / 2; // stride between segments
 * 	LIGOTimeGPS epoch = { 0, 0 };
 * 	REAL8TimeSeries *h[3]; // the strain induced in the network of detectors
 * 	gsl_rng *rng;
 * 	gsl_rng_env_setup();
 * 	rng = gsl_rng_alloc(gsl_rng_default);
 * 	h[0] = XLALCreateREAL8TimeSeries("H1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
 * 	h[1] = XLALCreateREAL8TimeSeries("L1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
 * 	h[2] = XLALCreateREAL8TimeSeries("V1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
 * 	XLALSimSGWBFlatSpectrum(h, detectors, 3, 0, Omega0, flow, H0, rng); // first time to initialize
 * 	while (1) { // infinite loop
 * 		size_t j;
 * 		for (j = 0; j < stride; ++j) // output first stride points
 * 			printf("%.9f\t%e\t%e\t%e\n", XLALGPSGetREAL8(&h[0]->epoch) + j / srate, h[0]->data->data[j], h[1]->data->data[j], h[2]->data->data[j]);
 * 		XLALSimSGWBFlatSpectrum(h, detectors, 3, stride, Omega0, flow, H0, rng); // make more data
 * 	}
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
int XLALSimSGWBFlatSpectrum(
	REAL8TimeSeries **h,			/**< [in/out] array of sgwb timeseries for detector network */
	const LALDetector *detectors,		/**< [in] array of detectors in network */
	size_t numDetectors,			/**< [in] number of detectors in network */
	size_t stride,				/**< [in] stride (samples) */
	double Omega0,				/**< [in] flat sgwb spectrum power (dimensionless) */
	double flow,				/**< [in] low frequency cutoff (Hz) */
	double H0,				/**< [in] Hubble's constant (s) */
	gsl_rng *rng				/**< [in] GSL random number generator */
)
{
	REAL8FrequencySeries *OmegaGW;
	size_t length;
	double deltaF;
	length = h[0]->data->length;
	deltaF = 1.0/(length * h[0]->deltaT);
	OmegaGW = XLALSimSGWBOmegaGWFlatSpectrum(Omega0, flow, deltaF, length/2 + 1);
	if (! OmegaGW)
		XLAL_ERROR(XLAL_EFUNC);
	if (XLALSimSGWB(h, detectors, numDetectors, stride, OmegaGW, H0, rng)) {
		XLALDestroyREAL8FrequencySeries(OmegaGW);
		XLAL_ERROR(XLAL_EFUNC);
	}
	XLALDestroyREAL8FrequencySeries(OmegaGW);
	return 0;
}


/**
 * Routine that may be used to generate sequential segments of stochastic
 * background gravitational wave signals for a network of detectors with a
 * specified stride from one segment to the next.
 *
 * The spectrum is a power law with power alpha for frequencies above flow with
 * power given by Omegaref at the reference frequency fref, and zero for
 * frequencies below the low frequency cutoff flow.
 *
 * Calling instructions: for the first call, set stride = 0; subsequent calls
 * should pass the same time series and have non-zero stride.  This routine
 * will advance the time series by an amount given by the stride and will
 * generate new data so that the data is continuous from one segment to the
 * next.  For example: the following routine will output a continuous stream of
 * stochastic background signals with a "flat" (OmegaGW = const) spectrum for
 * the HLV network.
 *
 * \code
 * #include <stdio.h>
 * #include <gsl/gsl_rng.h>
 * #include <lal/LALStdlib.h>
 * #include <lal/LALDetectors.h>
 * #include <lal/FrequencySeries.h>
 * #include <lal/TimeSeries.h>
 * #include <lal/Units.h>
 * #include <lal/LALSimSGWB.h>
 * int mksgwbdata_powerlaw(void)
 * {
 * 	const double flow = 40.0; // 40 Hz low frequency cutoff
 * 	const double duration = 16.0; // 16 second segments
 * 	const double srate = 16384.0; // sampling rate in Hertz
 * 	const double Omegaref = 1e-6; // fraction of critical energy in GWs
 * 	const double fref = 100; // reference frequency in Hertz
 * 	const double alpha = 3.0; // sgwb spectrum power law power
 * 	const double H0 = 0.72 * LAL_H0FAC_SI; // Hubble's constant in seconds
 * 	const LALDetector H1 = lalCachedDetectors[LAL_LHO_4K_DETECTOR]; // Hanford
 * 	const LALDetector L1 = lalCachedDetectors[LAL_LLO_4K_DETECTOR]; // Livingston
 * 	const LALDetector V1 = lalCachedDetectors[LAL_VIRGO_DETECTOR]; // Virgo
 * 	LALDetector detectors[3] = {H1, L1, V1}; // the network of detectors
 * 	size_t length = duration * srate; // segment length
 * 	size_t stride = length / 2; // stride between segments
 * 	LIGOTimeGPS epoch = { 0, 0 };
 * 	REAL8TimeSeries *h[3]; // the strain induced in the network of detectors
 * 	gsl_rng *rng;
 * 	gsl_rng_env_setup();
 * 	rng = gsl_rng_alloc(gsl_rng_default);
 * 	h[0] = XLALCreateREAL8TimeSeries("H1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
 * 	h[1] = XLALCreateREAL8TimeSeries("L1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
 * 	h[2] = XLALCreateREAL8TimeSeries("V1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
 * 	XLALSimSGWBPowerLawSpectrum(h, detectors, 3, 0, Omegaref, alpha, fref, flow, H0, rng); // first time to initialize
 * 	while (1) { // infinite loop
 * 		size_t j;
 * 		for (j = 0; j < stride; ++j) // output first stride points
 * 			printf("%.9f\t%e\t%e\t%e\n", XLALGPSGetREAL8(&h[0]->epoch) + j / srate, h[0]->data->data[j], h[1]->data->data[j], h[2]->data->data[j]);
 * 		XLALSimSGWBPowerLawSpectrum(h, detectors, 3, stride, Omegaref, alpha, fref, flow, H0, rng); // make more data
 * 	}
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
int XLALSimSGWBPowerLawSpectrum(
	REAL8TimeSeries **h,			/**< [in/out] array of sgwb timeseries for detector network */
	const LALDetector *detectors,		/**< [in] array of detectors in network */
	size_t numDetectors,			/**< [in] number of detectors in network */
	size_t stride,				/**< [in] stride (samples) */
	double Omegaref,			/**< [in] sgwb spectrum power at reference frequency (dimensionless) */
	double alpha,				/**< [in] sgwb spectrum power power law */
	double fref,				/**< [in] sgwb spectrum reference frequency (Hz) */
	double flow,				/**< [in] low frequency cutoff (Hz) */
	double H0,				/**< [in] Hubble's constant (s) */
	gsl_rng *rng				/**< [in] GSL random number generator */
)
{
	REAL8FrequencySeries *OmegaGW;
	size_t length;
	double deltaF;
	length = h[0]->data->length;
	deltaF = 1.0/(length * h[0]->deltaT);
	OmegaGW = XLALSimSGWBOmegaGWPowerLawSpectrum(Omegaref, alpha, fref, flow, deltaF, length/2 + 1);
	if (! OmegaGW)
		XLAL_ERROR(XLAL_EFUNC);
	if (! XLALSimSGWB(h, detectors, numDetectors, stride, OmegaGW, H0, rng)) {
		XLALDestroyREAL8FrequencySeries(OmegaGW);
		XLAL_ERROR(XLAL_EFUNC);
	}
	XLALDestroyREAL8FrequencySeries(OmegaGW);
	return 0;
}


/*
 *
 * TEST CODE
 *
 */

#if 0

#include <stdio.h>
#include <lal/TimeSeries.h>

/* Example routine listed in documentation. */
int mksgwbdata(void)
{
	const double flow = 40.0; // 40 Hz low frequency cutoff
	const double duration = 16.0; // 16 second segments
	const double srate = 16384.0; // sampling rate in Hertz
	const double Omega0 = 1e-6; // fraction of critical energy in GWs
	const double H0 = 0.72 * LAL_H0FAC_SI; // Hubble's constant in seconds
	const LALDetector H = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
	const LALDetector L = lalCachedDetectors[LAL_LLO_4K_DETECTOR];
	const LALDetector V = lalCachedDetectors[LAL_VIRGO_DETECTOR];
	LALDetector detectors[3] = {H, L, V};
	size_t length = duration * srate; // segment length
	size_t stride = length / 2; // stride between segments
	LIGOTimeGPS epoch = { 0, 0 };
	REAL8TimeSeries *seg[3];
	gsl_rng *rng;
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
	seg[0] = XLALCreateREAL8TimeSeries("H1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
	seg[1] = XLALCreateREAL8TimeSeries("L1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
	seg[2] = XLALCreateREAL8TimeSeries("V1:STRAIN", &epoch, 0.0, 1.0/srate, &lalStrainUnit, length);
	XLALSimSGWBFlatSpectrum(seg, detectors, 3, 0, Omega0, flow, H0, rng); // first time to initialize
	while (1) { // infinite loop
		size_t j;
		for (j = 0; j < stride; ++j) // output first stride points
			printf("%.9f\t%e\t%e\t%e\n", XLALGPSGetREAL8(&seg[0]->epoch) + j / srate, seg[0]->data->data[j], seg[1]->data->data[j], seg[2]->data->data[j]);
		XLALSimSGWBFlatSpectrum(seg, detectors, 3, stride, Omega0, flow, H0, rng); // make more data
	}
}

/*
 * Test routine that generates 1024 seconds of data in 8 second segments with
 * 50% overlap.  Also produced are the PSD for the spectra as well as the
 * normalized cross-spectral densities (to check the overlap reduction
 * function).
 */
int test_sgwb(void)
{
	const double recdur = 1024.0; // duration of the entire data record in seconds
	const double segdur = 8.0; // duration of a segment in seconds
	const double srate = 4096.0; // sample rate in Hertz
	const double flow = 1.0; // low frequency cutoff in Hertz
	const double Omega0 = 1e-6;
	const double H0 = 0.72 * LAL_H0FAC_SI;
	const LIGOTimeGPS epoch = {0, 0};
	enum { numDetectors = 4 };
	const LALDetector H = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
	const LALDetector L = lalCachedDetectors[LAL_LLO_4K_DETECTOR];
	const LALDetector V = lalCachedDetectors[LAL_VIRGO_DETECTOR];
	const LALDetector T = lalCachedDetectors[LAL_TAMA_300_DETECTOR];
	LALDetector detectors[numDetectors] = {H, L, V, H};
	REAL8TimeSeries *seg[numDetectors];
	REAL8TimeSeries *rec[numDetectors];
	REAL8FrequencySeries *psd[numDetectors];
	REAL8FrequencySeries *OmegaGW;
	COMPLEX16FrequencySeries *htilde1;
	COMPLEX16FrequencySeries *htilde2;
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
	size_t i, j, k, l;
	FILE *fp;

	if (reclen != (numseg - 1)*stride + seglen) {
		fprintf(stderr, "warning: numseg, seglen, reclen, and stride not commensurate\n");
		reclen = (numseg - 1)*stride + seglen;
	}

	rng = gsl_rng_alloc(gsl_rng_default);
	for (i = 0; i < numDetectors; ++i) {
		const char *prefix = detectors[i].frDetector.prefix;
		seg[i] = XLALCreateREAL8TimeSeries("STRAIN", &epoch, 0.0, deltaT, &lalStrainUnit, seglen);
		rec[i] = XLALCreateREAL8TimeSeries("STRAIN", &epoch, 0.0, deltaT, &lalStrainUnit, reclen);
		psd[i] = XLALCreateREAL8FrequencySeries("PSD", &epoch, 0.0, deltaF, &lalSecondUnit, seglen/2 + 1);
		snprintf(seg[i]->name, sizeof(seg[i]->name), "%2s:%s", prefix, seg[i]->name);
	}

	OmegaGW = XLALSimSGWBOmegaGWFlatSpectrum(Omega0, flow, deltaF, seglen/2 + 1);

	for (l = 0; l < numseg; ++l) {
		/* first time, initialize; subsequent times, stride */
		XLALSimSGWB(seg, detectors, numDetectors, l ? stride : 0, OmegaGW, H0, rng);
		/* copy stride amount of data or seglen on last time */
		for (i = 0; i < numDetectors; ++i)
			memcpy(rec[i]->data->data + l*stride, seg[i]->data->data, (l == numseg - 1 ? seglen : stride) * sizeof(*rec[i]->data->data));
	}

	fp = fopen("sgwb.dat", "w");
	for (j = 0; j < reclen; ++j) {
		fprintf(fp, "%.9f", j * deltaT);
		for (i = 0; i < numDetectors; ++i)
			fprintf(fp, "\t%e", rec[i]->data->data[j]);
		fprintf(fp, "\n");
	}
	fclose(fp);


	plan = XLALCreateForwardREAL8FFTPlan(seglen, 0);
	window = XLALCreateHannREAL8Window(seglen);
	for (i = 0; i < numDetectors; ++i)
		XLALREAL8AverageSpectrumWelch(psd[i], rec[i], seglen, stride, window, plan);

	/* compute PSDs for each detector */
	fp = fopen("sgwb-psd.dat", "w");
	for (k = klow; k < seglen/2; ++k) {
		double f = k * deltaF;
		double actual = 5.6e-22 * (H0 / LAL_H0FAC_SI) * sqrt(Omega0) * pow(100.0/f, 1.5);
		fprintf(fp, "%.9f\t%e", f, actual);
		for (i = 0; i < numDetectors; ++i)
			fprintf(fp, "\t%e", sqrt(psd[i]->data->data[k]));
		fprintf(fp, "\n");
	}
	fclose(fp);

	/* compute cross spectral densities */
	htilde1 = XLALCreateCOMPLEX16FrequencySeries("htilde1", &epoch, 0.0, deltaF, &lalSecondUnit, seglen/2 + 1);
	htilde2 = XLALCreateCOMPLEX16FrequencySeries("htilde2", &epoch, 0.0, deltaF, &lalSecondUnit, seglen/2 + 1);
	memset(psd[0]->data->data, 0, psd[0]->data->length * sizeof(*psd[0]->data->data));
	for (l = 0; l < numseg; ++l) {
		double fac = seglen / (window->sumofsquares * numseg);
		/* window data */
		for (j = 0; j < seglen; ++j) {
			seg[0]->data->data[j] = window->data->data[j] * rec[0]->data->data[j + l*stride];
			seg[1]->data->data[j] = window->data->data[j] * rec[1]->data->data[j + l*stride];
		}
		
		XLALREAL8TimeFreqFFT(htilde1, seg[0], plan);
		XLALREAL8TimeFreqFFT(htilde2, seg[1], plan);
		for (k = klow; k < seglen/2; ++k) {
			double re = creal(htilde1->data->data[k]);
			double im = cimag(htilde1->data->data[k]);
			psd[0]->data->data[k] += fac * (re * re + im * im);
		}
	}
	fp = fopen("sgwb-orf.dat", "w");
	for (k = klow; k < seglen/2; ++k) {
		double f = k * deltaF;
		double csdfac = 0.15 * pow(H0 / LAL_PI, 2.0) * pow(f, -3.0) * Omega0 * segdur;
		fprintf(fp, "%e\t%e\t%e\n", f, XLALSimSGWBOverlapReductionFunction(f, &detectors[0], &detectors[1]), psd[0]->data->data[k] / csdfac);
	}
	fclose(fp);

	XLALDestroyCOMPLEX16FrequencySeries(htilde2);
	XLALDestroyCOMPLEX16FrequencySeries(htilde1);

	XLALDestroyREAL8Window(window);
	XLALDestroyREAL8FFTPlan(plan);
	XLALDestroyREAL8FrequencySeries(OmegaGW);
	for (i = 0; i < numDetectors; ++i) {
		XLALDestroyREAL8FrequencySeries(psd[i]);
		XLALDestroyREAL8TimeSeries(rec[i]);
		XLALDestroyREAL8TimeSeries(seg[i]);
	}
	gsl_rng_free(rng);

	return 0;
}

int main(void)
{
	lalDebugLevel = 7;
	XLALSetErrorHandler(XLALAbortErrorHandler);
	gsl_rng_env_setup();
	mksgwbdata();
	test_sgwb();
	return 0;
}

#endif
