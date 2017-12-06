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

#ifndef _LALSIMSGWB_H
#define _LALSIMSGWB_H

#include <stddef.h>
#include <gsl/gsl_rng.h>
#include <lal/LALDatatypes.h>
#include <lal/LALDetectors.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif


/*
 *
 *
 * OVERLAP REDUCTION FUNCTION ROUTINE
 *
 *
 */


/**
 * Computes the overlap reduction function between two detectors at a specified
 * frequency.
 *
 * Implements the formulae given in Allen & Romano (1999).
 */
double XLALSimSGWBOverlapReductionFunction(
	double f,			/**< [in] frequency (Hz) */
	const LALDetector *detector1,	/**< [in] 1st detector */
	const LALDetector *detector2	/**< [in] 2nd detector */
);


/*
 *
 *
 * ROUTINES TO GENERATE SGWB SPECTRA
 *
 *
 */


/**
 * Creates a frequency series that contains a flat SGWB spectrum with the
 * specified power Omega0 above some low frequency cutoff flow.
 */
REAL8FrequencySeries *XLALSimSGWBOmegaGWFlatSpectrum(
	double Omega0,	/**< [in] sgwb spectrum power (dimensionless) */
	double flow,	/**< [in] low frequncy cutoff of SGWB spectrum (Hz) */
	double deltaF,	/**< [in] frequency bin width (Hz) */
	size_t length	/**< [in] number of frequency bins */
);


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
);


/*
 *
 * SGWB GENERATION ROUTINES
 *
 */

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
 * - If stride = 0, initialize h by generating one (periodic)
 * realization of noise; subsequent calls should have non-zero
 * stride.
 *
 * - If stride = h->data->length then generate one segment of
 * non-periodic noise by generating two different realizations
 * and feathering them together.
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
);

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
 * - If stride = 0, initialize h by generating one (periodic)
 * realization of noise; subsequent calls should have non-zero
 * stride.
 *
 * - If stride = h->data->length then generate one segment of
 * non-periodic noise by generating two different realizations
 * and feathering them together.
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
);

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
 * - If stride = 0, initialize h by generating one (periodic)
 * realization of noise; subsequent calls should have non-zero
 * stride.
 *
 * - If stride = h->data->length then generate one segment of
 * non-periodic noise by generating two different realizations
 * and feathering them together.
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
);


#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMSGWB_H */
