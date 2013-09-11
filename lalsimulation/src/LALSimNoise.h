/*
 * Copyright (C) 2011 J. Creighton
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

#ifndef _LALSIMNOISE_H
#define _LALSIMNOISE_H

#include <stddef.h>
#include <lal/LALDatatypes.h>
#include <gsl/gsl_rng.h>

#if defined(__cplusplus)
extern "C" {
#elif 0
} /* so that editors will match preceding brace */
#endif

#ifdef SWIG // SWIG interface directives
SWIGLAL(FUNCTION_POINTER(XLALSimNoisePSDiLIGOSRD, XLALSimNoisePSDiLIGOModel,
	XLALSimNoisePSDeLIGOModel, XLALSimNoisePSDVirgo, XLALSimNoisePSDGEO,
	XLALSimNoisePSDTAMA, XLALSimNoisePSDaLIGONoSRMLowPower,
	XLALSimNoisePSDaLIGONoSRMHighPower,
	XLALSimNoisePSDaLIGOZeroDetLowPower,
	XLALSimNoisePSDaLIGOZeroDetHighPower, XLALSimNoisePSDaLIGONSNSOpt,
	XLALSimNoisePSDaLIGOBHBH20Deg, XLALSimNoisePSDaLIGOHighFrequency,
	XLALSimNoisePSDKAGRA, XLALSimNoisePSDAdvVirgo));
#endif

/*
 *
 *
 * PSD GENERATION FUNCTIONS
 *
 *
 */


/*
 *
 * FUNCTIONS TO GENERATE COMPONENT NOISE PSD
 *
 */

/**
 * Provides a rather ad-hoc estimate of the seismic noise power spectral density
 * at a given frequency.
 *
 * This is a crude estimate based on characteristic frequencies for the
 * pendulum and the stack.  What is computed is
 * \f[
 * S_h(f) = L^{-2} S_g(f) (f_{\mathrm{pend}}/f)^4
 * (f_{\mathrm{stack}}/f)^{4n_{\mathrm{stack}}}
 * \f]
 * where
 * \f[
 * S_g(f) = 10^{-18}\,\mathrm{m}^2\,\mathrm{Hz}^{-1}\times(10\,\mathrm{Hz}/f)^4
 * \f]
 * is the displacement power spectrum of ground motion.
 *
 * Warning: the transfer function is only correct at frequencies above the
 * specified characteristic pendulum and stack frequencies.
 */
double XLALSimNoisePSDSeismic(
	double f,		/**< frequency (Hz) */
	double L,		/**< arm length (m) */
	double f_pend,		/**< characteristic frequency of pendulum suspension (Hz) */
	double f_stack,		/**< characteristic frequency of isolation stack (Hz) */
	double n_stack		/**< number of layers of stack */
);

/**
 * Provides a rather ad-hoc estimate of the suspension thermal noise power
 * spectral density at a given frequency.
 *
 * This is a crude estimate based on the characteristic frequency of the
 * pendulum suspension and its quality factor (= 1 / loss angle).  What is
 * computed is
 * \f[
 * S_h(f) = L^{-2} \frac{2 k T}{\pi^3 f_0^3 M Q} \left( \frac{f_0}{f} \right)^5.
 * \f]
 *
 * Warning: this only describes the broadband noise at frequencies above the
 * pendulum frequency; it does not have the correct noise near the resonances.
 */
double XLALSimNoisePSDSuspTherm(
	double f,		/**< frequency (Hz) */
	double L,		/**< arm length (m) */
	double M,		/**< mirror mass (kg) */
	double T,		/**< temperature (K) */
	double f0,		/**< pendulum frequency */
	double Q		/**< pendulum quality */
);

/**
 * Provides a rather ad-hoc estimate of the mirror thermal noise power spectral
 * density at a given frequency.
 *
 * This is a crude estimate based on the characteristic frequency of the
 * mirror/coating internal modes and their quality factor (= 1 / loss angle).
 * What is computed is
 * \f[
 * S_h(f) = L^{-2} \frac{2 k T}{\pi^3 f_0^3 M Q} \frac{f_0}{f}
 * \f]
 *
 * Warning: this only describes the broadband noise at frequencies below the
 * resonance frequency; it does not have the correct noise near the resonances.
 */
double XLALSimNoisePSDMirrorTherm(
	double f,		/**< frequency (Hz) */
	double L,		/**< arm length (m) */
	double M,		/**< mirror mass (kg) */
	double T,		/**< average per mirror power loss */
	double f0,		/**< average per mirror power loss */
	double Q		/**< average per mirror power loss */
);

/**
 * Computes the shot noise in strain-equivalent units using a conventional
 * model appropriate to initial interferometric detectors.
 *
 * Uses the formula for shot noise from
 *
 */
double XLALSimNoisePSDShot(
	double f,		/**< frequency (Hz) */
	double P_BS,		/**< laser power on beamsplitter (W) */
	double lambda,		/**< laser wavelength (m) */
	double L,		/**< arm length (m) */
	double finesse,		/**< arm cavity finesse */
	double eta		/**< effective quantum efficiency of photodiode */
);

/**
 * Computes the quantum noise (shot noise and radiation pressure noise)
 * according to Buonanno and Chen, Phys. Rev. D 64 0402006 (2001).
 *
 * This code is adapted from the GWINC matlab function shotrad.m which includes
 * updated losses by Kirk McKenzie.
 *
 * For simplicity, only losses from the mirrors are included.  Losses from
 * coupling and from the SRC are ignored.  (These could be included as
 * effective losses in A_BS if needed.) A fixed photdiode quantum efficiency of
 * eta = 0.9 is used.
 *
 * Note: this code is adapted from GWINC.
 */
double XLALSimNoisePSDQuantum(
	double f,		/**< frequency (Hz) */
	double I0,		/**< laser power (W) */
	double lambda,		/**< laser wavelength (m) */
	double L,		/**< arm length (m) */
	double M,		/**< mirror mass (kg) */
	double A,		/**< average per mirror power loss */
	double A_BS,		/**< power loss at beam splitter */
	double T_ITM,		/**< transmittance of ITM */
	double T_PRM,		/**< transmittance of PRM */
	double T_SRM,		/**< transmittance of SRM */
	double ds,		/**< detuning phase (rad) */
	double zeta,		/**< demod/detection/homodyne phase */
	double eta		/**< quantum efficiency of photodiode */
);

/*
 *
 * NOISE PSD ROUTINES FOR FIRST GENERATION DETECTORS
 *
 */


/**
 * Provides the noise power spectrum based on a phenomenological fit
 * to the SRD curve for iLIGO.
 *
 * This is a fit to the data provided for the Science Requirements Document
 * (SRD) curve for initial LIGO given, which can be found at
 * http://www.ligo.caltech.edu/~jzweizig/distribution/LSC_Data/
 *
 * The Science Requirements Document is located at
 * http://www.ligo.caltech.edu/docs/E/E950018-02.pdf
 */
double XLALSimNoisePSDiLIGOSRD(double f /**< frequency (Hz) */);

/**
 * Provides the seismic noise power spectrum for iLIGO.
 *
 * Note: only valit for f > 10 Hz.
 * This is mostly a phenomenological fit.
 */
double XLALSimNoisePSDiLIGOSeismic(double f /**< frequency (Hz) */);

/**
 * Provides the thermal noise (suspension + coating) power spectrum for iLIGO.
 *
 * Note: this is a phenomenological fit to the broadband component.
 */
double XLALSimNoisePSDiLIGOThermal(double f /**< frequency (Hz) */);

/**
 * Provides the shot noise power spectrum for iLIGO.
 *
 * Note: the effective quantum efficiency is one-third the actual quantum
 * efficiency owing to the RF readout scheme.  A fiducial value of 250 W
 * of power on the beamsplitter is used.
 */
double XLALSimNoisePSDiLIGOShot(double f /**< frequency (Hz) */);

/**
 * Provides the shot noise power spectrum for eLIGO.
 *
 * Note: A fiducial value of 250 W of power on the beamsplitter is used.
 */
double XLALSimNoisePSDeLIGOShot(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for a model of the iLIGO detector.
 *
 * Warning: not all noise sources are correctly accounted for (in particular,
 * there is no actuation noise modelled) so this noise spectrum does not
 * correspond to the S5 spectrum.
 */
double XLALSimNoisePSDiLIGOModel(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for a model of the eLIGO detector.
 *
 * Warning: not all noise sources are correctly accounted for so this noise
 * spectrum does not correspond to the S6 spectrum.
 */
double XLALSimNoisePSDeLIGOModel(double f /**< frequency (Hz) */);

/**
 * Provides the design noise power spectrum for Virgo based on a
 * phenomenological fit (from the Virgo webiste) that can be approximated by the
 * following:
 * \f{equation}{
 * S_h(f) =
 * s_0 \left ( \frac {7.87f}{f_0} \right )^{-4.8} + \frac{6}{17} \frac{f_0}{f}
 * + \left [1 + \left (\frac {f}{f_0} \right)^2 \right ],
 * \f}
 * where \f$s_0=10.2e-46\f$.
 *
 * Warning: This comes from the deprecated function LALVIRGOPsd in the lal
 * noisemodels package, which comes with no reference to the curve. An updated
 * version of this model, with a reference would be welcomed.
 */
double XLALSimNoisePSDVirgo(double f /**< frequency (Hz) */);

/**
 * Provides a GEO noise power spectrum based on that from Table IV of
 * \cite dis2001.
 *
 * The comes from the deprecated function LALGEOPsd in the lal noisemodels
 * package.
 */
double XLALSimNoisePSDGEO(double f /**< frequency (Hz) */);

/**
 * Provides a GEO-HF noise power spectrum based on a fit to Figure 6
 * from \cite Grote2010.
 *
 * The fit is good between 50Hz to 8kHz and errors between the analytic
 * fit given and the <a href="https://intranet.aei.uni-hannover.de/geo600/geohflogbook.nsf/7e8722dffa24dea0c1256de900406c84/4837a612ac990060c12575ce004e70fd?OpenDocument">estimated curve</a> are less than 1%.
 */
double XLALSimNoisePSDGEOHF(double f /**< frequency (Hz) */);

/**
 * Provides a TAMA300 noise power spectrum based on that from Table IV of
 * \cite dis2001.
 *
 * The comes from the deprecated function LALTAMAPsd in the lal noisemodels
 * package.
 */
double XLALSimNoisePSDTAMA(double f /**< frequency (Hz) */);

/*
 *
 * NOISE PSD ROUTINES FOR SECOND GENERATION DETECTORS
 *
 */


/**
 * Provides the thermal noise (suspension + coating) power spectrum for iLIGO.
 *
 * Note: this is a phenomenological fit to the broadband component.
 */
double XLALSimNoisePSDaLIGOThermal(double f /**< frequency (Hz) */);

/**
 * Provides the quantum noise power spectrum for aLIGO under the low-power
 * no-signal-recycling-mirror configuration.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled No SRM.
 */
double XLALSimNoisePSDaLIGOQuantumNoSRMLowPower(double f /**< frequency (Hz) */);

/**
 * Provides the quantum noise power spectrum for aLIGO under the high-power
 * no-signal-recycling-mirror configuration.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is the same a No SRM but with 125 W laser power.
 */
double XLALSimNoisePSDaLIGOQuantumNoSRMHighPower(double f /**< frequency (Hz) */);

/**
 * Provides the quantum noise power spectrum for aLIGO under the low-power
 * broad-band signal recycling (no detuning of the signal recycling cavity).
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled Zero Detune, Low Power.
 */
double XLALSimNoisePSDaLIGOQuantumZeroDetLowPower(double f /**< frequency (Hz) */);

/**
 * Provides the quantum noise power spectrum for aLIGO under the high-power
 * broad-band signal recycling (no detuning of the signal recycling cavity).
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled Zero Detune, High Power.
 */
double XLALSimNoisePSDaLIGOQuantumZeroDetHighPower(double f /**< frequency (Hz) */);

/**
 * Provides the quantum noise power spectrum for aLIGO under the
 * configuration tuned to optimize sensitivity to NS-NS inspirals.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled NS-NS Opt.
 */
double XLALSimNoisePSDaLIGOQuantumNSNSOpt(double f /**< frequency (Hz) */);

/**
 * Provides the quantum noise power spectrum for aLIGO under the
 * configuration tuned to optimize sensitivity to 30+30 solar mass binary
 * black holes with fixed signal recycling cavity detuning of 20 degrees.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled BHBH 20-degree Detune.
 */
double XLALSimNoisePSDaLIGOQuantumBHBH20Deg(double f /**< frequency (Hz) */);

/**
 * Provides the quantum noise power spectrum for aLIGO under the
 * configuration tuned to narrow-band high-frequency sensitivity around
 * 1 kHz.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled High Freq.
 */
double XLALSimNoisePSDaLIGOQuantumHighFrequency(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for aLIGO under the low-power
 * no-signal-recycling-mirror configuration.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled No SRM.
 *
 * Note: This includes only thermal and quantum noise.
 */
double XLALSimNoisePSDaLIGONoSRMLowPower(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for aLIGO under the high-power
 * no-signal-recycling-mirror configuration.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is the same a No SRM but with 125 W laser power.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGONoSRMHighPower(double f /**< frequency (Hz) */);


/**
 * Provides the noise power spectrum for aLIGO under the low-power
 * broad-band signal recycling (no detuning of the signal recycling cavity).
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled Zero Detune, Low Power.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGOZeroDetLowPower(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for aLIGO under the high-power
 * broad-band signal recycling (no detuning of the signal recycling cavity).
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled Zero Detune, High Power.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGOZeroDetHighPower(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for aLIGO under the
 * configuration tuned to optimize sensitivity to NS-NS inspirals.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled NS-NS Opt.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGONSNSOpt(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for aLIGO under the
 * configuration tuned to optimize sensitivity to 30+30 solar mass binary
 * black holes with fixed signal recycling cavity detuning of 20 degrees.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled BHBH 20-degree Detune.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGOBHBH20Deg(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for aLIGO under the
 * configuration tuned to narrow-band high-frequency sensitivity around
 * 1 kHz.
 *
 * See: LIGO-T0900288-v3 and LIGO-T070247-01.
 * This configuration is labelled High Freq.
 *
 * Warning: This includes only thermal and quantum noise.  It is only valid
 * above around 9 Hz.
 */
double XLALSimNoisePSDaLIGOHighFrequency(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for KAGRA based on that from Eqn 5 of
 * \cite md2012. This is a phenomenological fit to the KAGRA spectrum from
 * http://gwcenter.icrr.u-tokyo.ac.jp/en/researcher/parameter.
 */
double XLALSimNoisePSDKAGRA(double f /**< frequency (Hz) */);

/**
 * Provides the noise power spectrum for AdvVirgo based on that from Eqn 6 of
 * \cite md2012. This is a phenomenological fit to the AdvVirgo spectrum from
 * http://wwwcascina.virgo.infin.it/advirgo.
 */
double XLALSimNoisePSDAdvVirgo(double f /**< frequency (Hz) */);


/*
 *
 * NOISE PSD UTILITY ROUTINES
 *
 */

/**
 * Evaluates a power spectral density function, psdfunc, at the frequencies required
 * to populate the frequency series psd, with a low frequency cutoff flow.
 */
int XLALSimNoisePSD(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow,			/**< low frequency cutoff (Hz) */
	double (*psdfunc)(double)	/**< function that provides the PSD at a specified frequency */
);

/**
 * Reads file fname containing two-column amplitude spectral density data file
 * and interpolates at the frequencies required to populate the frequency
 * series psd, with a low frequency cutoff flow.
 */
int XLALSimNoisePSDFromFile(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow,			/**< low frequency cutoff (Hz) */
	const char *fname		/**< file containing amplitude spectral density data */
);

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "NO_SRM.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGONoSRMLowPowerGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
);

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "ZERO_DET_low_P.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGOZeroDetLowPowerGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
);

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "ZERO_DET_high_P.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGOZeroDetHighPowerGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
);

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "NSNS_Opt.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGONSNSOptGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
);

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "BBH_20deg.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGOBHBH20DegGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
);

/**
 * Returns a frequency series psd with low frequency cutoff flow corresponding
 * to the "High_Freq.txt" data file in LIGO-T0900288.
 */
int XLALSimNoisePSDaLIGOHighFrequencyGWINC(
	REAL8FrequencySeries *psd,	/**< frequency series to be computed */
	double flow 			/**< low frequency cutoff (Hz) */
);

/*
 *
 * NOISE GENERATION ROUTINES
 *
 */


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
int XLALSimNoise(
	REAL8TimeSeries *s,		/**< [in/out] noise time series */
	size_t stride,			/**< [in] stride (samples) */
	REAL8FrequencySeries *psd,	/**< [in] power spectrum frequency series */
	gsl_rng *rng			/**< [in] GSL random number generator */
);

#if 0
{ /* so that editors will match succeeding brace */
#elif defined(__cplusplus)
}
#endif

#endif /* _LALSIMNOISE_H */
