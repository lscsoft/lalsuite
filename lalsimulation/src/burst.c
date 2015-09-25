/*
*  Copyright (C) 2015 Jolien Creighton
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

/**
 * @defgroup lalsim_burst lalsim-burst
 * @ingroup lalsimulation_programs
 *
 * @brief Simulates a generic burst gravitational waveformS
 *
 * ### Synopsis
 *
 *     lalsim-burst [options]
 *
 * ### Description
 *
 * The `lalsim-burst` utility produces a stream of a simulated gravitational
 * waveform for a generic burst.  The output data is the gravitational waveform
 * polarizations in the time domain.  The output is written to standard output
 * as a three-column ascii format.  The first column gives the time
 * corresponding to each sample and the remaining columns give the
 * gravitational waveform values for the two polarizations.
 *
 * ### Options
 * [default values in brackets]
 *
 * <DL>
 * <DT>`-h`, `--help`
 * <DD>print a help message and exit</DD>
 * <DT>`-v`, `--verbose`
 * <DD>verbose output</DD>
 * <DT>`-w` WAVEFM, `--waveform=WAVEFM`</DT>
 * <DD>the waveform to generate</DD>
 * <DT>`-t` DT, `--duration=DT`</DT>
 * <DD>duration (seconds)</DD>
 * <DT>`-f` FREQ, `--frequency=FREQ`</DT>
 * <DD>frequency (Hz)</DD>
 * <DT>`-b` BW, `--bandwidth=BW`</DT>
 * <DD>bandwidth (Hz)</DD>
 * <DT>`-q` Q, `--quality-factor=Q`</DT>
 * <DD>quality factor (dimensionless)</DD>
 * <DT>`-e` ECC, `--eccentricity=ECC`</DT>
 * <DD>eccentricity 0<=ECC<=1 [0]</DD>
 * <DT>`-p` PHI, `--phase=PHI`</DT>
 * <DD>phase (degrees) [0]</DD>
 * <DT>`-A` AMP, `--amplitude=AMP`</DT>
 * <DD>amplitude (dimensionless or \f$ {\rm s}^{-1/3} \f$)</DD>
 * <DT>`-H` HRSS, `--hrss=HRSS`</DT>
 * <DD>root-sum-squared amplitude (\f$ {\rm Hz}^{-1/2} \f$)</DD>
 * <DT>`-F` FLUENCE, `--fluence=`FLUENCE</DT>
 * <DD>isotropic energy fluence (\f$ M_\odot c^2 / {\rm pc}^2 \f$)</DD>
 * <DT>`-R` SRATE, --sample-rate=SRATE</DT>
 * <DD>sample rate (Hz) [16384]</DD>
 * </DL>
 *
 * ### Waveforms Supported
 *
 * <DL>
 * <DT>BLTWNB</DT>
 * <DD>band-limited white-noise burst<BR>
 * required parameters: duration, frequency, bandwidth, eccentricity, fluence
 * </DD>
 * <DT>StringCusp</DT>
 * <DD>cosmic string cusp<BR>
 * required parameters: amplitude, frequency
 * </DD>
 * <DT>SineGaussian</DT>
 * <DD>cosine- or sine-Gaussian<BR>
 * required parameters: quality-factor, frequency, hrss, eccentricity, phase
 * </DD>
 * <DT>Gaussian</DT>
 * <DD>Gaussian<BR>
 * required parameters: duration, hrss
 * </DD>
 * <DT>Impulse</DT>
 * <DD>delta-function impulse<BR>
 * required parameters: amplitude
 * </DD>
 * </DL>
 *
 * ### Environment
 *
 * The `LAL_DEBUG_LEVEL` can used to control the error and warning reporting of
 * `lalsim-burst`.  Common values are: `LAL_DEBUG_LEVEL=0` which suppresses
 * error messages, `LAL_DEBUG_LEVEL=1`  which prints error messages alone,
 * `LAL_DEBUG_LEVEL=3` which prints both error messages and warning messages,
 * and `LAL_DEBUG_LEVEL=7` which additionally prints informational messages.
 *
 * The `GSL_RNG_SEED` and `GSL_RNG_TYPE` environment variables can be used
 * to set the random number generator seed and type respectively when
 * generating band-limited white-noise bursts.
 *
 *
 * ### Exit Status
 *
 * The `lalsim-burst` utility exits 0 on success, and >0 if an error occurs.
 *
 * ### Example
 *
 * The command:
 *
 *     lalsim-burst --waveform BLTWNB --duration 0.1 --frequency 100 --bandwidth 100 --fluence 1e-14
 *
 * produces a three-column ascii output to standard output; the rows are
 * samples (at the default rate of 16384 Hz), and the three columns are 1. the
 * time of each sample, 2. the plus-polarization strain, and 3. the
 * cross-polarization strain.  The waveform produced is a band-limited
 * white-noise burst with equal power in its two polarizations
 * (eccentricity=0, which is the default value), a duration of 0.1 seconds,
 * a frequency of 100 Hz, a bandwidth of 100 Hz, and a fluence of
 * \f$ 10^{-14}\, M_\odot c^2 / {\rm pc}^2 = 0.01\, M_\odot c^2 / {\rm Mpc}^2
 * \simeq 1.9\, {\rm J} / {\rm m}^2 \f$.
 *
 * The command:
 *
 *     lalsim-burst -w SineGaussian -q 9 -e 1 -p 90 -f 150 -H 1e-22
 *
 * produces a Q=9 linearly polarized (since the eccentricity is set to 1)
 * sine-Gaussian (not a cosine-Gaussian since the phase is set to 90 degrees)
 * centered at 150 Hz having root-sum-squared strain
 * \f$ h_{\mathrm{rss}} = 10^{-22}\, {\rm Hz}^{-1/2} \f$.
 */

#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LALgetopt.h>
#include <lal/Date.h>
#include <lal/TimeSeries.h>
#include <lal/LALSimBurst.h>

#define XSTR(x) #x
#define STR(x) XSTR(x)
#define INVALID_DOUBLE (nan(""))
#define IS_INVALID_DOUBLE(x) (isnan(x))
#define DEFAULT_ECCENTRICITY 0 /* circularly polarized / unpolarized */
#define DEFAULT_PHASE 0 /* cosine phase */
#define DEFAULT_SRATE 16384

static int verbose = 0;
#define verbose_output(...) (verbose ? fprintf(stderr, __VA_ARGS__) : 0)

static int waveform_value(const char *waveform);
static double my_atof(const char *s);
static int output(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross);
static int usage(const char *program);
static struct params parseargs(int argc, char **argv);

typedef enum {
	BLTWNB,
	StringCusp,
	SineGaussian,
	Gaussian,
	Impulse,
	NumWaveforms
} waveform_enum;

#define INIT_NAME(a) [a] = #a
static const char *waveform_names[NumWaveforms] = {
	INIT_NAME(BLTWNB),
	INIT_NAME(StringCusp),
	INIT_NAME(SineGaussian),
	INIT_NAME(Gaussian),
	INIT_NAME(Impulse)
};
#undef INIT_NAME

static const char *waveform_long_names[NumWaveforms] = {
	[BLTWNB] = "band-limited white-noise burst",
	[StringCusp] = "cosmic string cusp",
	[SineGaussian] = "cosine- or sine-Gaussian",
	[Gaussian] = "Gaussian",
	[Impulse] = "delta-function impulse"
};

static const char *waveform_parameters[NumWaveforms] = {
	[BLTWNB] = "duration, frequency, bandwidth, eccentricity, fluence",
	[StringCusp] = "amplitude, frequency",
	[SineGaussian] = "quality-factor, frequency, hrss, eccentricity, phase",
	[Gaussian] = "duration, hrss",
	[Impulse] = "amplitude"
};

static int waveform_value(const char *waveform)
{
	int i;
	for (i = 0; i < NumWaveforms; ++i)
		if (strcmp(waveform, waveform_names[i]) == 0)
			return i;
	fprintf(stderr, "error: waveform `%s' not recognized\n", waveform);
	return -1;
}

/* like atof but prints an error message and returns NaN if an error occurs */
static double my_atof(const char *s)
{
	int save_errno = errno;
	char *endptr;
	double val;
	errno = 0;
	val = strtod(s, &endptr);
	if (errno != 0 || *endptr != '\0') {
		save_errno = errno;
		fprintf(stderr, "error: could not parse string `%s' into a double\n", s);
		val = INVALID_DOUBLE;
	}
	errno = save_errno;
	return val;
}

struct params {
	int waveform;
	double duration;
	double frequency;
	double bandwidth;
	double q;
	double phase;
	double eccentricity;
	double amplitude;
	double hrss;
	double fluence;
	double srate;
};

int main(int argc, char **argv)
{
	REAL8TimeSeries *hplus = NULL;
	REAL8TimeSeries *hcross = NULL;
	struct params p;
	int status = 0;
	gsl_rng *rng = NULL;

	XLALSetErrorHandler(XLALAbortErrorHandler);

	p = parseargs(argc, argv);

	if ((int)(p.waveform) < 0 || (int)(p.waveform) >= NumWaveforms) {
		fprintf(stderr, "error: must specify valid waveform\n");
		exit(1);
	}

	if (p.waveform == BLTWNB) { /* set up gsl random number generator */
		gsl_rng_env_setup();
		rng = gsl_rng_alloc(gsl_rng_default);
	}

	if (IS_INVALID_DOUBLE(p.srate) || p.srate <= 0.0) {
		fprintf(stderr, "error: must specify valid sample rate\n");
		status = 1;
	}

	verbose_output("Input parameters:\n");
	verbose_output("%-31s `%s' - %s\n", "waveform:", waveform_names[p.waveform], waveform_long_names[p.waveform]);
	verbose_output("%-31s %g (Hz)\n", "sample rate:", p.srate);

	switch (p.waveform) {

	case BLTWNB:

		/* sanity check relevant parameters */
		if (IS_INVALID_DOUBLE(p.duration) || p.duration < 0.0) {
			fprintf(stderr, "error: must specify valid duration for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.frequency) || p.frequency < 0.0) {
			fprintf(stderr, "error: must specify valid frequency for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.bandwidth) || p.bandwidth < 0.0) {
			fprintf(stderr, "error: must specify valid bandwidth for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.eccentricity) || p.eccentricity < 0.0 || p.eccentricity > 1.0) {
			fprintf(stderr, "error: must specify valid eccentricity in domain [0,1] for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.fluence) || p.fluence < 0.0) {
			fprintf(stderr, "error: must specify valid fluence for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}

		/* detect set but ignored parameters */
		if (!IS_INVALID_DOUBLE(p.q))
			fprintf(stderr, "warning: quality-factor parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (p.phase != DEFAULT_PHASE)
			fprintf(stderr, "warning: phase parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.amplitude))
			fprintf(stderr, "warning: amplitude parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.hrss))
			fprintf(stderr, "warning: hrss parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);

		if (!status) {
			double int_hdot_squared_dt;
			int_hdot_squared_dt = p.fluence * LAL_GMSUN_SI * 4 / LAL_C_SI / LAL_PC_SI / LAL_PC_SI;

			verbose_output("%-31s %g (s)\n", "duration:", p.duration);
			verbose_output("%-31s %g (Hz)\n", "frequency:", p.frequency);
			verbose_output("%-31s %g (Hz)\n", "bandwidth:", p.bandwidth);
			verbose_output("%-31s %g\n", "eccentricity:", p.eccentricity);
			verbose_output("%-31s %g (Msun c^2 pc^-2)\n", "fluence:", p.fluence);
			verbose_output("%-31s %g (s^-1)\n", "integral (dh/dt)^2 dt:", int_hdot_squared_dt);
			verbose_output("%-31s GSL_RNG_TYPE=%s\n", "GSL random number generator:", gsl_rng_name(rng));
			verbose_output("%-31s GSL_RNG_SEED=%lu\n", "GSL random number seed:", gsl_rng_default_seed);
			status = XLALGenerateBandAndTimeLimitedWhiteNoiseBurst(&hplus, &hcross, p.duration, p.frequency, p.bandwidth, p.eccentricity, int_hdot_squared_dt, 1.0/p.srate, rng);
		}
		break;

	case StringCusp:

		/* sanity check relevant parameters */
		if (IS_INVALID_DOUBLE(p.amplitude) || p.amplitude < 0.0) {
			fprintf(stderr, "error: must specify valid amplitude for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.frequency) || p.frequency < 0.0) {
			fprintf(stderr, "error: must specify valid frequency for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}

		/* detect set but ignored parameters */
		if (!IS_INVALID_DOUBLE(p.duration))
			fprintf(stderr, "warning: duration parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.bandwidth))
			fprintf(stderr, "warning: bandwidth parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.q))
			fprintf(stderr, "warning: quality-factor parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (p.eccentricity != DEFAULT_ECCENTRICITY)
			fprintf(stderr, "warning: eccentricity parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (p.phase != DEFAULT_PHASE)
			fprintf(stderr, "warning: phase parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.hrss))
			fprintf(stderr, "warning: hrss parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.fluence))
			fprintf(stderr, "warning: fluence parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);

		if (!status) {
			verbose_output("%-31s %g (s^-1/3)\n", "amplitude:", p.amplitude);
			verbose_output("%-31s %g (Hz)\n", "frequency:", p.frequency);
			status = XLALGenerateStringCusp(&hplus, &hcross, p.amplitude, p.frequency, 1.0/p.srate);
		}
		break;

	case SineGaussian:

		/* sanity check relevant parameters */
		if (IS_INVALID_DOUBLE(p.q) || p.q < 0.0) {
			fprintf(stderr, "error: must specify valid quality factor for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.frequency) || p.frequency < 0.0) {
			fprintf(stderr, "error: must specify valid frequency for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.hrss) || p.hrss < 0.0) {
			fprintf(stderr, "error: must specify valid hrss for waveform `%s\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.eccentricity) || p.eccentricity < 0.0 || p.eccentricity > 1.0) {
			fprintf(stderr, "error: must specify valid eccentricity in domain [0,1] for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.phase)) {
			fprintf(stderr, "error: must specify valid phase for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}

		/* detect set but ignored parameters */
		if (!IS_INVALID_DOUBLE(p.duration))
			fprintf(stderr, "warning: duration parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.bandwidth))
			fprintf(stderr, "warning: bandwidth parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.amplitude))
			fprintf(stderr, "warning: amplitude parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.fluence))
			fprintf(stderr, "warning: fluence parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);

		if (!status) {
			verbose_output("%-31s %g\n", "quality-factor:", p.q);
			verbose_output("%-31s %g (Hz)\n", "frequency:", p.frequency);
			verbose_output("%-31s %g (Hz^-1/2)\n", "root-sum-squared strain:", p.hrss);
			verbose_output("%-31s %g\n", "eccentricity:", p.eccentricity);
			verbose_output("%-31s %g (degrees)\n", "phase:", p.phase / LAL_PI_180);
			status = XLALSimBurstSineGaussian(&hplus, &hcross, p.q, p.frequency, p.hrss, p.eccentricity, p.phase, 1.0/p.srate);
		}
		break;

	case Gaussian:

		/* sanity check relevant parameters */
		if (IS_INVALID_DOUBLE(p.duration) || p.duration < 0.0) {
			fprintf(stderr, "error: must specify valid duration for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		if (IS_INVALID_DOUBLE(p.hrss) || p.hrss < 0.0) {
			fprintf(stderr, "error: must specify valid hrss for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}
		/* detect set but ignored parameters */
		if (!IS_INVALID_DOUBLE(p.frequency))
			fprintf(stderr, "warning: frequency parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.bandwidth))
			fprintf(stderr, "warning: bandwidth parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.q))
			fprintf(stderr, "warning: quality-factor parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (p.eccentricity != DEFAULT_ECCENTRICITY)
			fprintf(stderr, "warning: eccentricity parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (p.phase != DEFAULT_PHASE)
			fprintf(stderr, "warning: phase parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.amplitude))
			fprintf(stderr, "warning: amplitude parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.fluence))
			fprintf(stderr, "warning: fluence parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!status) {
			verbose_output("%-31s %g (s)\n", "duration:", p.duration);
			verbose_output("%-31s %g (Hz^-1/2)\n", "root-sum-squared strain:", p.hrss);
			status = XLALSimBurstGaussian(&hplus, &hcross, p.duration, p.hrss, 1.0/p.srate);
		}

		break;

	case Impulse:

		/* sanity check relevant parameters */
		if (IS_INVALID_DOUBLE(p.amplitude) || p.amplitude < 0.0) {
			fprintf(stderr, "error: must specify valid amplitude for waveform `%s'\n", waveform_names[p.waveform]);
			status = 1;
		}

		/* detect set but ignored parameters */
		if (!IS_INVALID_DOUBLE(p.duration))
			fprintf(stderr, "warning: duration parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.frequency))
			fprintf(stderr, "warning: frequency parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.bandwidth))
			fprintf(stderr, "warning: bandwidth parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.q))
			fprintf(stderr, "warning: quality-factor parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (p.eccentricity != DEFAULT_ECCENTRICITY)
			fprintf(stderr, "warning: eccentricity parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (p.phase != DEFAULT_PHASE)
			fprintf(stderr, "warning: phase parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.hrss))
			fprintf(stderr, "warning: hrss parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);
		if (!IS_INVALID_DOUBLE(p.fluence))
			fprintf(stderr, "warning: fluence parameter is set but ignored for waveform `%s'\n", waveform_names[p.waveform]);

		if (!status) {
			verbose_output("%-31s %g (dimensionless)\n", "amplitude:", p.amplitude);
			status = XLALGenerateImpulseBurst(&hplus, &hcross, p.amplitude, 1.0/p.srate);
		}
		break;

	default:
		fprintf(stderr, "error: unrecognized waveform\n");
		exit(1);
	};

	if (status)
		exit(1);

	if (verbose) {
		char peak_time[32]; // GPS time string - 31 characters is enough
		LIGOTimeGPS tpeak;
		COMPLEX16 hpeak;
		double hrss;
		double fluence;
		unsigned ipeak;
		hpeak = XLALMeasureHPeak(hplus, hcross, &ipeak);
		tpeak = hplus->epoch;
		XLALGPSAdd(&tpeak, ipeak * hplus->deltaT);
		XLALGPSToStr(peak_time, &tpeak);
		hrss = XLALMeasureHrss(hplus, hcross);
		fluence = XLALMeasureEoverRsquared(hplus, hcross);
		verbose_output("Measured parameters:\n");
		verbose_output("%-31s %s (s)\n", "peak time:", peak_time);
		verbose_output("%-31s (h+, hx) = (%g, %g)\n", "peak strain amplitude:", creal(hpeak), cimag(hpeak));
		verbose_output("%-31s abs(h+, hx) = %g\n", "peak strain amplitude:", cabs(hpeak));
		verbose_output("%-31s arg(h+, hx) = %g (rad)\n", "peak strain amplitude:", carg(hpeak));
		verbose_output("%-31s %g (Hz^-1/2)\n", "root-sum-squared strain:", hrss);
		verbose_output("%-31s %g (J m^-2)\n", "isotropic energy fluence:", fluence);
		verbose_output("%-31s %g (Msun c^2 pc^-2)\n", "isotropic energy fluence:", fluence * LAL_PC_SI * LAL_PC_SI / LAL_MSUN_SI / LAL_C_SI / LAL_C_SI);
	}

	output(hplus, hcross);

	XLALDestroyREAL8TimeSeries(hcross);
	XLALDestroyREAL8TimeSeries(hplus);
	LALCheckMemoryLeaks();
	return 0;
}

static int output(REAL8TimeSeries *hplus, REAL8TimeSeries *hcross)
{
	char tstr[32]; // GPS time string - 31 characters is enough
	size_t j;
	fprintf(stdout, "# time (s)\tH_PLUS (strain)\tH_CROSS (strain)\n");
	for (j = 0; j < hplus->data->length; ++j) {
		LIGOTimeGPS t = hplus->epoch;
		fprintf(stdout, "%s\t%e\t%e\n", XLALGPSToStr(tstr, XLALGPSAdd(&t, j * hplus->deltaT)), hplus->data->data[j], hcross->data->data[j]);
	}
	return 0;
}

static int usage(const char *program)
{
	int i;
	const char *options =
/* *INDENT-OFF* */
"[default values in brackets]\n"
"	-h, --help                     print this message and exit\n"
"	-v, --verbose                  verbose output\n"
"	-w WAVEFM, --waveform=WAVEFM   the waveform to generate\n"
"	-t DT, --duration=DT           duration (seconds)\n"
"	-f FREQ, --frequency=FREQ      frequency (Hz)\n"
"	-b BW, --bandwidth=BW          bandwidth (Hz)\n"
"	-q Q, --quality-factor=Q       quality factor (dimensionless)\n"
"	-e ECC, --eccentricity=ECC     eccentricity 0<=ECC<=1 [" STR(DEFAULT_ECCENTRICITY) "]\n"
"	-p PHI, --phase=PHI            phase (degrees) [" STR(DEFAULT_PHASE) "]\n"
"	-A AMP, --amplitude=AMP        amplitude (dimensionless or s^-1/3)\n"
"	-H HRSS, --hrss=HRSS           root-sum-squared amplitude (Hz^-1/2)\n"
"	-F FLUENCE, --fluence=FLUENCE  isotropic energy fluence (Msun c^2 pc^-2)\n"
"	-R SRATE, --sample-rate=SRATE  sample rate (Hz) [" STR(DEFAULT_SRATE) "]\n";
/* *INDENT-ON* */
	fprintf(stderr, "usage: %s [options]\n", program);
	fprintf(stderr, "options:\n%s\n", options);
	fprintf(stderr, "waveforms supported:\n");
	for (i = 0; i < NumWaveforms; ++i) {
		fprintf(stderr, "\n\t`%s' - %s\n", waveform_names[i], waveform_long_names[i]);
		fprintf(stderr, "\t\trequires: %s\n", waveform_parameters[i]);
	}
	return 0;
}

static struct params parseargs(int argc, char **argv)
{
	struct params p = {
		.waveform = -1,				/* invalid */
		.duration = INVALID_DOUBLE,		/* invalid */
		.frequency = INVALID_DOUBLE,	 	/* invalid */
		.bandwidth = INVALID_DOUBLE,		/* invalid */
		.q = INVALID_DOUBLE,			/* invalid */
		.phase = DEFAULT_PHASE,
		.eccentricity = DEFAULT_ECCENTRICITY,
		.amplitude = INVALID_DOUBLE,		/* invalid */
		.hrss = INVALID_DOUBLE,			/* invalid */
		.fluence = INVALID_DOUBLE,		/* invalid */
		.srate = DEFAULT_SRATE
	};
	struct LALoption long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"verbose", no_argument, 0, 'v'},
		{"waveform", required_argument, 0, 'w'},
		{"duration", required_argument, 0, 't'},
		{"frequency", required_argument, 0, 'f'},
		{"bandwidth", required_argument, 0, 'b'},
		{"quality-factor", required_argument, 0, 'q'},
		{"eccentricity", required_argument, 0, 'e'},
		{"phase", required_argument, 0, 'p'},
		{"amplitude", required_argument, 0, 'A'},
		{"hrss", required_argument, 0, 'H'},
		{"fluence", required_argument, 0, 'F'},
		{"sample-rate", required_argument, 0, 'R'},
		{0, 0, 0, 0}
	};
	char args[] = "hvw:t:f:b:q:e:A:H:F:d:R:";

	while (1) {
		int option_index = 0;
		int c;

		c = LALgetopt_long_only(argc, argv, args, long_options, &option_index);
		if (c == -1) /* end of options */
			break;

		switch (c) {
		case 0:		/* if option set a flag, nothing else to do */
			if (long_options[option_index].flag)
				break;
			else {
				fprintf(stderr, "error parsing option %s with argument %s\n", long_options[option_index].name, LALoptarg);
				exit(1);
			}
		case 'h':	/* help */
			usage(argv[0]);
			exit(0);
		case 'v':	/* verbose */
			verbose = 1;
			break;
		case 'w':	/* waveform */
			p.waveform = waveform_value(LALoptarg);
			break;
		case 't':	/* duration */
			p.duration = my_atof(LALoptarg);
			break;
		case 'f':	/* frequency */
			p.frequency = my_atof(LALoptarg);
			break;
		case 'b':	/* bandwidth */
			p.bandwidth = my_atof(LALoptarg);
			break;
		case 'q':	/* quality-factor */
			p.q = my_atof(LALoptarg);
			break;
		case 'e':	/* eccentricity */
			p.eccentricity = my_atof(LALoptarg);
			break;
		case 'p':	/* phase */
			/* convert phase from degrees to radians */
			p.phase = my_atof(LALoptarg) * LAL_PI_180;
			break;
		case 'A':	/* amplitude */
			p.amplitude = my_atof(LALoptarg);
			break;
		case 'H':	/* hrss */
			p.hrss = my_atof(LALoptarg);
			break;
		case 'F':	/* fluence */
			/* convert fluence to sum-squared hdot */
			p.fluence = my_atof(LALoptarg);
			break;
		case 'R':	/* sample-rate */
			p.srate = my_atof(LALoptarg);
			break;
		case '?':
		default:
			fprintf(stderr, "unknown error while parsing options\n");
			exit(1);
		}
	}
	if (LALoptind < argc) {
		fprintf(stderr, "extraneous command line arguments:\n");
		while (LALoptind < argc)
			fprintf(stderr, "%s\n", argv[LALoptind++]);
		exit(1);
	}
	return p;
}
