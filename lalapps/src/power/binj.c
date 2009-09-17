/*
 * Copyright (C) 2007 Kipp Cannon, Lisa M. Goggin, Patrick Brady, Saikat
 * Ray-Majumder, Xavier Siemens
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with with program; see the file COPYING. If not, write to the Free
 * Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
 */


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include <lal/Date.h>
#include <lal/GenerateBurst.h>
#include <lal/LALConstants.h>
#include <lal/LALSimBurst.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/TimeDelay.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>
#include <lal/lalGitID.h>


#include <lalapps.h>
#include <processtable.h>
#include <lalappsGitID.h>


RCSID("$Id$");


#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_binj"


/*
 * ============================================================================
 *
 *                                Command Line
 *
 * ============================================================================
 */


enum population {
	POPULATION_TARGETED,
	POPULATION_ALL_SKY_SINEGAUSSIAN,
	POPULATION_ALL_SKY_BTLWNB,
	POPULATION_STRING_CUSP
};


struct options {
	INT8 gps_start_time;
	INT8 gps_end_time;
	enum population population;
	double ra;
	double dec;
	double maxA;
	double minA;
	double maxbandwidth;
	double minbandwidth;
	double maxduration;
	double minduration;
	double maxEoverr2;
	double minEoverr2;
	double maxf;
	double minf;
	double maxhrss;
	double minhrss;
	double q;
	unsigned long seed;
	double time_step;
	char *user_tag;
};


static struct options options_defaults(void)
{
	struct options defaults;
	
	defaults.gps_start_time = -1;
	defaults.gps_end_time = -1;
	defaults.population = -1;
	defaults.ra = XLAL_REAL8_FAIL_NAN;
	defaults.dec = XLAL_REAL8_FAIL_NAN;
	defaults.maxbandwidth = XLAL_REAL8_FAIL_NAN;
	defaults.minbandwidth = XLAL_REAL8_FAIL_NAN;
	defaults.maxduration = XLAL_REAL8_FAIL_NAN;
	defaults.minduration = XLAL_REAL8_FAIL_NAN;
	defaults.maxEoverr2 = XLAL_REAL8_FAIL_NAN;
	defaults.minEoverr2 = XLAL_REAL8_FAIL_NAN;
	defaults.maxf = XLAL_REAL8_FAIL_NAN;
	defaults.minf = XLAL_REAL8_FAIL_NAN;
	defaults.maxhrss = XLAL_REAL8_FAIL_NAN;
	defaults.minhrss = XLAL_REAL8_FAIL_NAN;
	defaults.minA = XLAL_REAL8_FAIL_NAN;
	defaults.maxA = XLAL_REAL8_FAIL_NAN;
	defaults.q = XLAL_REAL8_FAIL_NAN;
	defaults.seed = 0;
	defaults.time_step = 210.0 / LAL_PI;
	defaults.user_tag = NULL;

	return defaults;
}


static void print_usage(void)
{
	fprintf(stderr, 
"lalapps_binj [options]\n" \
"\n" \
"Options:\n" \
"\n" \
"--gps-end-time seconds\n" \
"--gps-start-time seconds\n" \
"	Bounds of interval in which to synthesize injections.\n" \
"\n" \
"--help\n" \
"	display this message\n" \
"\n" \
"--max-amplitude value\n" \
"--min-amplitude value\n" \
"	Set the bounds of the injection ampltiudes.  These only affect\n" \
"	string cusp injections.\n" \
"\n" \
"--max-bandwidth hertz\n" \
"--min-bandwidth hertz\n" \
"	Set the bounds of the injection bandwidthds.  These only affect\n" \
"	btlwnb waveforms.\n" \
"\n" \
	); fprintf(stderr, 
"--max-duration seconds\n" \
"--min-duration seconds\n" \
"	Set the bounds of the injection durations.  These only affect\n" \
"	btlwnb waveforms.\n" \
"\n" \
"--max-e-over-r2 value\n" \
"--min-e-over-r2 value\n" \
"	Set the bounds of the range of equivalent isotropic radiated\n" \
"	energies of btlwnb waveforms.  The units are M_{sun} / pc^{2} (solar\n" \
"	masses per parsec^{2}).\n" \
"\n" \
	); fprintf(stderr, 
"--max-frequency hertz\n" \
"--min-frequency hertz\n" \
"	Set the bounds of the injection frequencies.  These are the centre\n" \
"	frequencies of sine-Gaussians and btlwnb waveforms, and the\n" \
"	high-frequency cut-offs of string cusp waveforms.\n" \
"\n" \
"--max-hrss value\n" \
"--min-hrss value\n" \
	); fprintf(stderr, 
"	Set the bounds of the injection h_{rss} values.  These only affect\n" \
"	sine-Gaussian injections.  (Actually, these set the bounds of the\n" \
"	product of the waveform's hrss and its duration, which makes the\n" \
"	injections lie along the 50 efficiency curve better.  To convert to\n" \
"	real hrss multiply by sqrt(2) pi f/Q.) \n" \
"\n" \
	); fprintf(stderr, 
"--population name\n" \
"	Select the injection population to synthesize.  Allowed values are\n" \
"	\"targeted\", \"string_cusp\", and \"all_sky_sinegaussian\",\n" \
"	\"all_sky_btlwnb\".\n" \
"\n" \
"--q value\n" \
"	Set the Q for sine-Gaussian injections.\n" \
"\n" \
"--ra-dec ra,dec\n" \
"	Set the right-ascension and declination of the sky location from which\n" \
"	injections should originate when generating a targeted population.\n" \
"	Co-ordinates are in radians.\n" \
"\n" \
	); fprintf(stderr, 
"--seed value\n" \
"	Set the random number generator's seed (0 = off, default = 0).\n" \
"\n" \
"--time-step value\n" \
"	Set the time betwen injections to value/pi seconds (default = 210).\n" \
"\n" \
"--user-tag string\n" \
"	Set the user tag in the process and search summary tables to this.\n"
	);
}


static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param, const ProcessTable *process, const char *type, const char *param, const char *value)
{
	*proc_param = XLALCreateProcessParamsTableRow(process);
	snprintf((*proc_param)->program, sizeof((*proc_param)->program), PROGRAM_NAME);
	snprintf((*proc_param)->type, sizeof((*proc_param)->type), type);
	snprintf((*proc_param)->param, sizeof((*proc_param)->param), "--%s", param);
	snprintf((*proc_param)->value, sizeof((*proc_param)->value), value);

	return(&(*proc_param)->next);
}

	
#define ADD_PROCESS_PARAM(process, type) \
	do { paramaddpoint = add_process_param(paramaddpoint, process, type, long_options[option_index].name, optarg); } while(0)


static struct options parse_command_line(int *argc, char **argv[], const ProcessTable *process, ProcessParamsTable **paramaddpoint)
{
	struct options options = options_defaults();
	int c;
	int option_index;
	struct option long_options[] = {
		{"gps-end-time", required_argument, NULL, 'A'},
		{"gps-start-time", required_argument, NULL, 'B'},
		{"help", no_argument, NULL, 'C'},
		{"max-amplitude", required_argument, NULL, 'D'},
		{"min-amplitude", required_argument, NULL, 'E'},
		{"max-bandwidth", required_argument, NULL, 'F'},
		{"min-bandwidth", required_argument, NULL, 'G'},
		{"max-duration", required_argument, NULL, 'H'},
		{"min-duration", required_argument, NULL, 'I'},
		{"max-e-over-r2", required_argument, NULL, 'S'},
		{"min-e-over-r2", required_argument, NULL, 'T'},
		{"max-frequency", required_argument, NULL, 'J'},
		{"min-frequency", required_argument, NULL, 'K'},
		{"max-hrss", required_argument, NULL, 'L'},
		{"min-hrss", required_argument, NULL, 'M'},
		{"population", required_argument, NULL, 'N'},
		{"q", required_argument, NULL, 'O'},
		{"ra-dec", required_argument, NULL, 'U'},
		{"seed", required_argument, NULL, 'P'},
		{"time-step", required_argument, NULL, 'Q'},
		{"user-tag", required_argument, NULL, 'R'},
		{NULL, 0, NULL, 0}
	};

	do switch(c = getopt_long(*argc, *argv, "", long_options, &option_index)) {
	case 'A':
		XLALClearErrno();
		{
			LIGOTimeGPS tmp;
			XLALStrToGPS(&tmp, optarg, NULL);
			options.gps_end_time = XLALGPSToINT8NS(&tmp);
		}
		if(xlalErrno) {
			fprintf(stderr, "invalid --%s (%s specified)\n", long_options[option_index].name, optarg);
			exit(1);
		}
		ADD_PROCESS_PARAM(process, "lstring");
		break;

	case 'B':
		XLALClearErrno();
		{
			LIGOTimeGPS tmp;
			XLALStrToGPS(&tmp, optarg, NULL);
			options.gps_start_time = XLALGPSToINT8NS(&tmp);
		}
		if(xlalErrno) {
			fprintf(stderr, "invalid --%s (%s specified)\n", long_options[option_index].name, optarg);
			exit(1);
		}
		ADD_PROCESS_PARAM(process, "lstring");
		break;

	case 'C':
		print_usage();
		exit(0);

	case 'D':
		options.maxA = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'E':
		options.minA = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'F':
		options.maxbandwidth = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'G':
		options.minbandwidth = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'H':
		options.maxduration = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'I':
		options.minduration = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'J':
		options.maxf = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'K':
		options.minf = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'L':
		options.maxhrss = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'M':
		options.minhrss = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'N':
		if(!strcmp(optarg, "targeted"))
			options.population = POPULATION_TARGETED;
		else if(!strcmp(optarg, "string_cusp"))
			options.population = POPULATION_STRING_CUSP;
		else if(!strcmp(optarg, "all_sky_sinegaussian"))
			options.population = POPULATION_ALL_SKY_SINEGAUSSIAN;
		else if(!strcmp(optarg, "all_sky_btlwnb"))
			options.population = POPULATION_ALL_SKY_BTLWNB;
		else {
			fprintf(stderr, "error: unrecognized population \"%s\"", optarg);
			exit(1);
		}
		ADD_PROCESS_PARAM(process, "lstring");
		break;

	case 'O':
		options.q = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'P':
		options.seed = atol(optarg);
		ADD_PROCESS_PARAM(process, "int_8u");
		break;

	case 'Q':
		options.time_step = atof(optarg) / LAL_PI;
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'R':
		options.user_tag = optarg;
		ADD_PROCESS_PARAM(process, "lstring");
		break;

	case 'S':
		options.maxEoverr2 = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'T':
		options.minEoverr2 = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;

	case 'U':
		{
			char *end;
			options.ra = strtod(optarg, &end);
			while(isspace(*end))
				end++;
			if(*end != ',') {
				fprintf(stderr, "error: cannot parse --ra-dec \"%s\"\n", optarg);
				exit(1);
			}
			options.dec = strtod(end + 1, &end);
			while(isspace(*end))
				end++;
			if(*end != '\0') {
				fprintf(stderr, "error: cannot parse --ra-dec \"%s\"\n", optarg);
				exit(1);
			}
		}
		ADD_PROCESS_PARAM(process, "lstring");
		break;

	case 0:
		/* option sets a flag */
		break;

	case -1:
		/* end of arguments */
		break;

	case '?':
		/* unrecognized option */
		print_usage();
		exit(1);

	case ':':
		/* missing argument for an option */
		print_usage();
		exit(1);
	} while(c != -1);

	/* check some of the input parameters for consistency */
	if(options.maxA < options.minA) {
		fprintf(stderr, "error: --max-amplitude < --min-amplitude\n");
		exit(1);
	}
	if(options.maxbandwidth < options.minbandwidth) {
		fprintf(stderr, "error: --max-bandwidth < --min-bandwidth\n");
		exit(1);
	}
	if(options.maxduration < options.minduration) {
		fprintf(stderr, "error: --max-duration < --min-duration\n");
		exit(1);
	}
	if(options.maxf < options.minf) {
		fprintf(stderr, "error: --max-frequency < --min-frequency\n");
		exit(1);
	}
	if(options.maxhrss < options.minhrss) {
		fprintf(stderr, "error: --max-hrss < --min-hrss\n");
		exit(1);
	}

	if(options.gps_start_time == -1 || options.gps_end_time == -1) {
		fprintf(stderr, "--gps-start-time and --gps-end-time are both required\n");
		exit(1);
	}
	if(options.gps_end_time < options.gps_start_time) {
		fprintf(stderr, "error: --gps-end-time < --gps-start-time\n");
		exit(1);
	}

	switch(options.population) {
	case POPULATION_TARGETED:
	case POPULATION_ALL_SKY_SINEGAUSSIAN:
	case POPULATION_ALL_SKY_BTLWNB:
	case POPULATION_STRING_CUSP:
		break;

	default:
		fprintf(stderr, "error: --population is required\n");
		exit(1);
	}

	return options;
}


/* 
 * ============================================================================
 *
 *                                 Sequences
 *
 * ============================================================================
 */


/*
 * Repeating arithmetic sequence.
 */


static double sequence_arithmetic_next(double low, double high, double delta)
{
	static int i = 0;
	double x = low + delta * i++;

	if(high < low || high - low < delta)
		return XLAL_REAL8_FAIL_NAN;

	if(x > high) {
		i = 1;
		return low;
	}

	return x;
}


/*
 * Repeating geometric sequence.
 */


static double sequence_geometric_next(double low, double high, double ratio)
{
	static unsigned i = 0;
	double x = low * pow(ratio, i++);

	if(high < low || high / low < ratio)
		return XLAL_REAL8_FAIL_NAN;

	if(x > high) {
		i = 1;
		return low;
	}

	return x;
}


/*
 * Random amplitude presets.
 */


static double sequence_preset_next(gsl_rng *rng)
{
	static const double presets[] = {
		1e-22,
		2e-22,
		5e-22,
		1e-21,
		2e-21,
		5e-21,
		1e-20
	};

	return presets[gsl_rng_uniform_int(rng, sizeof(presets)/sizeof(*presets))];
}


/* 
 * ============================================================================
 *
 *                Logarithmically-Distributed Random Variable
 *
 * ============================================================================
 */


/*
 * Return a random number in the range [a, b), whose logarithm is uniformly
 * distributed
 */


static double ran_flat_log(gsl_rng *rng, double a, double b)
{
	return exp(gsl_ran_flat(rng, log(a), log(b)));
}


/*
 * Logarithmically-distributed random(ish) sequence.  Same as
 * sequence_geometric_next() but where each repetition of the sequence has
 * a random factor applied whose logarithm is uniformly distributed.  The
 * result is a random variable whose logarithm is uniformly distributed on
 * average, but in which neighbouring numbers are separated by a guaranteed
 * minimum ratio.
 */


static double ran_flat_log_discrete(gsl_rng *rng, double a, double b, double ratio)
{
	static double factor = 0.0;
	double x = sequence_geometric_next(a, b / ratio, ratio);

	if(x == a)
		/* sequence has looped.  must happen first time through */
		factor = ran_flat_log(rng, 1.0, ratio);

	return x * factor;
}


/* 
 * ============================================================================
 *
 *                          String Cusp Simulations
 *
 * ============================================================================
 */


/*
 * \theta is the angle the line of sight makes with the cusp.  \theta^{2}
 * is distributed uniformly, and the high frequency cut-off of the
 * injection is \propto \theta^{-3}.  So we first solve for the limits of
 * \theta^{2} from the low- and high bounds of the frequency band of
 * interest, pick a uniformly-distributed number in that range, and convert
 * back to a high frequency cut-off.
 */


static double random_string_cusp_fhigh(double flow, double fhigh, gsl_rng *rng)
{
	const double thetasqmin = pow(fhigh, -2.0 / 3.0);
	const double thetasqmax = pow(flow, -2.0 / 3.0);
	const double thetasq = gsl_ran_flat(rng, thetasqmin, thetasqmax);

	return pow(thetasq, -3.0 / 2.0);
}


/*
 * Uniform sky location, and uniform polarization orientation.
 */


static void random_location_and_polarization(double *ra, double *dec, double *psi, gsl_rng *rng)
{
	*ra = gsl_ran_flat(rng, 0, LAL_TWOPI);
	*dec = asin(gsl_ran_flat(rng, -1, +1));
	*psi = gsl_ran_flat(rng, 0, LAL_TWOPI);
}


/*
 * Pick a random string cusp injection.
 */


static SimBurst *random_string_cusp(double flow, double fhigh, double Alow, double Ahigh, gsl_rng *rng)
{
	SimBurst *sim_burst = XLALCreateSimBurst();

	if(!sim_burst)
		return NULL;

	strcpy(sim_burst->waveform, "StringCusp");

	/* high frequency cut-off and amplitude */

	sim_burst->frequency = random_string_cusp_fhigh(flow, fhigh, rng);
	sim_burst->amplitude = ran_flat_log(rng, Alow, Ahigh);

	/* sky location and wave frame orientation */

	random_location_and_polarization(&sim_burst->ra, &sim_burst->dec, &sim_burst->psi, rng);

	/* string cusp waveform generator makes a linearly polarized
	 * waveform in the + polarization.  it ignores these parameters,
	 * but just for consistency we populate them appropriately */

	sim_burst->pol_ellipse_e = 1.0;
	sim_burst->pol_ellipse_angle = 0.0;

	return sim_burst;
}


/* 
 * ============================================================================
 *
 *                             BTLWNB Injections
 *
 * ============================================================================
 */


/*
 * Pick a random band- and time-limited white noise burst.
 */


static SimBurst *random_directed_btlwnb(double ra, double dec, double psi, double minf, double maxf, double minband, double maxband, double mindur, double maxdur, double minEoverr2, double maxEoverr2, gsl_rng *rng)
{
	REAL8TimeSeries *hplus, *hcross;
	SimBurst *sim_burst = XLALCreateSimBurst();

	if(!sim_burst)
		return NULL;

	strcpy(sim_burst->waveform, "BTLWNB");

	/* sky location and wave frame orientation */

	sim_burst->ra = ra;
	sim_burst->dec = dec;
	sim_burst->psi = psi;

	/* pick a waveform */

	sim_burst->waveform_number = floor(gsl_ran_flat(rng, 0, ULONG_MAX));

	/* centre frequency.  three steps between minf and maxf */

	sim_burst->frequency = ran_flat_log_discrete(rng, minf, maxf, pow(maxf / minf, 1.0 / 3.0));

	/* duration and bandwidth.  keep picking until a valid pair is
	 * obtained (i.e. their product is >= 2 / \pi) */

	do {
		sim_burst->duration = ran_flat_log(rng, mindur, maxdur);
		sim_burst->bandwidth = ran_flat_log(rng, minband, maxband);
	} while(sim_burst->bandwidth * sim_burst->duration < LAL_2_PI);

	/* energy */

	sim_burst->egw_over_rsquared = ran_flat_log(rng, minEoverr2, maxEoverr2) * pow(sim_burst->frequency / 100.0, 4.0);

	/* populate the hrss column for convenience later */
	/* FIXME:  sample rate is hard-coded to 8192 Hz, which is what the
	 * excess power pipeline's .ini file is configured for in CVS, but
	 * because it can be easily changed this is not good */

	XLALGenerateSimBurst(&hplus, &hcross, sim_burst, 1.0 / 8192);
	sim_burst->hrss = XLALMeasureHrss(hplus, hcross);
	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);

	/* not sure if this makes sense.  these parameters are ignored by
	 * the injection code, but post-processing tools sometimes wish to
	 * know with what amplitude an injection should've been seen in an
	 * instrument, and they use these to project the + and x amplitudes
	 * onto the detector.  setting the eccentricity to 0 and the angle
	 * to anything makes such tools believe the amplitude is equally
	 * partitioned between the two polarizations which is true for
	 * these injections. */

	sim_burst->pol_ellipse_e = 0.0;
	sim_burst->pol_ellipse_angle = 0.0;

	/* done */

	return sim_burst;
}


static SimBurst *random_all_sky_btlwnb(double minf, double maxf, double minband, double maxband, double mindur, double maxdur, double minEoverr2, double maxEoverr2, gsl_rng *rng)
{
	double ra, dec, psi;
	SimBurst *sim_burst;

	random_location_and_polarization(&ra, &dec, &psi, rng);

	return random_directed_btlwnb(ra, dec, psi, minf, maxf, minband, maxband, mindur, maxdur, minEoverr2, maxEoverr2, rng);
}


/* 
 * ============================================================================
 *
 *                           All-Sky sine-Gaussians
 *
 * ============================================================================
 */


/*
 * the duration of a sine-Gaussian from its Q and centre frequency
 */


static double duration_from_q_and_f(double Q, double f)
{
	/* compute duration from Q and frequency */
	return Q / (sqrt(2.0) * LAL_PI * f);
}


/*
 * pick a sine-Gaussian
 */


static SimBurst *random_all_sky_sineGaussian(double minf, double maxf, double q, double minhrsst, double maxhrsst, gsl_rng *rng)
{
	SimBurst *sim_burst = XLALCreateSimBurst();

	if(!sim_burst)
		return NULL;

	strcpy(sim_burst->waveform, "SineGaussian");

	/* sky location and wave frame orientation */

	random_location_and_polarization(&sim_burst->ra, &sim_burst->dec, &sim_burst->psi, rng);

	/* q and centre frequency.  three steps between minf and maxf */

	sim_burst->q = q;
	sim_burst->frequency = ran_flat_log_discrete(rng, minf, maxf, pow(maxf / minf, 1.0 / 3.0));

	/* hrss */

	sim_burst->hrss = ran_flat_log(rng, minhrsst, maxhrsst) / duration_from_q_and_f(sim_burst->q, sim_burst->frequency);

	/* hard-code for linearly polarized waveforms in the x
	 * polarization.  induces LAL's sine-Gaussian generator to produce
	 * linearly polarized sine-Gaussians (+ would be a cosine
	 * Gaussian). */

	sim_burst->pol_ellipse_e = 1.0;
	sim_burst->pol_ellipse_angle = LAL_PI_2;

	/* done */

	return sim_burst;
}


/* 
 * ============================================================================
 *
 *                                   Output
 *
 * ============================================================================
 */


static void write_xml(const ProcessTable *process_table_head, const ProcessParamsTable *process_params_table_head, const SearchSummaryTable *search_summary_head, const SimBurst *sim_burst, struct options options)
{
	char fname[256];
	LIGOLwXMLStream *xml;

	if(options.user_tag)
		snprintf(fname, sizeof(fname), "HL-INJECTIONS_%s-%d-%d.xml", options.user_tag, (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));
	else
		snprintf(fname, sizeof(fname), "HL-INJECTIONS-%d-%d.xml", (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));

	xml = XLALOpenLIGOLwXMLFile(fname);

	/* process table */
	if(XLALWriteLIGOLwXMLProcessTable(xml, process_table_head)) {
		/* error occured.  ?? do anything else ?? */
		exit(1);
	}

	/* process params table */
	if(XLALWriteLIGOLwXMLProcessParamsTable(xml, process_params_table_head)) {
		/* error occured.  ?? do anything else ?? */
		exit(1);
	}

	/* search summary table */
	if(XLALWriteLIGOLwXMLSearchSummaryTable(xml, search_summary_head)) {
		/* error occured.  ?? do anything else ?? */
		exit(1);
	}

	/* sim burst table */
	if(XLALWriteLIGOLwXMLSimBurstTable(xml, sim_burst)) {
		/* error occured.  ?? do anything else ?? */
		exit(1);
	}

	XLALCloseLIGOLwXMLFile(xml);
}


/* 
 * ============================================================================
 *
 *                                Entry Point
 *
 * ============================================================================
 */


int main(int argc, char *argv[])
{
	struct options options;
	INT8 tinj;
	gsl_rng *rng;
	ProcessTable *process_table_head;
	ProcessParamsTable *process_params_table_head = NULL;
	SearchSummaryTable *search_summary_head;
	SimBurst *sim_burst_table_head = NULL;
	SimBurst **sim_burst = &sim_burst_table_head;


	/*
	 * Initialize debug handler
	 */


	lal_errhandler = LAL_ERR_EXIT;
	lalDebugLevel = LALINFO | LALWARNING | LALERROR | LALNMEMDBG | LALNMEMPAD | LALNMEMTRK;


	/*
	 * Process table
	 */


	process_table_head = XLALCreateProcessTableRow();
	if (strcmp(CVS_REVISION, "$Revi" "sion$"))
	{
		if(XLALPopulateProcessTable(process_table_head, PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE, 0))
			exit(1);
	}
	else
	{
		if(XLALPopulateProcessTable(process_table_head, PROGRAM_NAME, lalappsGitCommitID, lalappsGitGitStatus, lalappsGitCommitDate, 0))
			exit(1);
	}
	XLALGPSTimeNow(&process_table_head->start_time);
	snprintf(process_table_head->ifos, sizeof(process_table_head->ifos), "H1,H2,L1");


	/*
	 * Command line and process params table.
	 */


	options = parse_command_line(&argc, &argv, process_table_head, &process_params_table_head);
	if(options.user_tag)
		snprintf(process_table_head->comment, sizeof(process_table_head->comment), "%s", options.user_tag);


	/*
	 * Search summary table
	 */


	search_summary_head = XLALCreateSearchSummaryTableRow(process_table_head);
	if(options.user_tag)
		snprintf(search_summary_head->comment, sizeof(search_summary_head->comment), "%s", options.user_tag);
	search_summary_head->nnodes = 1;
	search_summary_head->out_start_time = *XLALINT8NSToGPS(&search_summary_head->in_start_time, options.gps_start_time);
	search_summary_head->out_end_time = *XLALINT8NSToGPS(&search_summary_head->in_end_time, options.gps_end_time);
	snprintf(search_summary_head->ifos, sizeof(search_summary_head->ifos), process_table_head->ifos);


	/*
	 * Initialize random number generator
	 */


	rng = gsl_rng_alloc(gsl_rng_mt19937);
	if(options.seed)
		gsl_rng_set(rng, options.seed);


	/*
	 * Main loop
	 */


	for(tinj = options.gps_start_time; tinj <= options.gps_end_time; tinj += options.time_step * 1e9) {
		/*
		 * Progress bar.
		 */

		XLALPrintInfo("%s: ", argv[0]);
		XLALPrintProgressBar((tinj - options.gps_start_time) / (double) (options.gps_end_time - options.gps_start_time));
		XLALPrintInfo(" complete\n");

		/*
		 * Create an injection
		 */

		switch(options.population) {
		case POPULATION_TARGETED:
			*sim_burst = random_directed_btlwnb(options.ra, options.dec, gsl_ran_flat(rng, 0, LAL_TWOPI), options.minf, options.maxf, options.minbandwidth, options.maxbandwidth, options.minduration, options.maxduration, options.minEoverr2, options.maxEoverr2, rng);
			break;

		case POPULATION_ALL_SKY_SINEGAUSSIAN:
			*sim_burst = random_all_sky_sineGaussian(options.minf, options.maxf, options.q, options.minhrss, options.maxhrss, rng);
			break;

		case POPULATION_ALL_SKY_BTLWNB:
			*sim_burst = random_all_sky_btlwnb(options.minf, options.maxf, options.minbandwidth, options.maxbandwidth, options.minduration, options.maxduration, options.minEoverr2, options.maxEoverr2, rng);
			break;

		case POPULATION_STRING_CUSP:
			*sim_burst = random_string_cusp(options.minf, options.maxf, options.minA, options.maxA, rng);
			break;

		default:
			/* shouldn't get here, command line parsing code
			 * should prevent it */
			XLALPrintError("internal error\n");
			exit(1);
		}

		if(!*sim_burst) {
			XLALPrintError("can't make injection\n");
			exit(1);
		}

		/*
		 * Peak time at geocentre in GPS seconds
		 */

		XLALINT8NSToGPS(&(*sim_burst)->time_geocent_gps, tinj);

		/*
		 * Peak time at geocentre in GMST radians
		 */

		(*sim_burst)->time_geocent_gmst = XLALGreenwichMeanSiderealTime(&(*sim_burst)->time_geocent_gps);

		/*
		 * Move to next injection
		 */

		sim_burst = &(*sim_burst)->next;
	}

	XLALPrintInfo("%s: ", argv[0]);
	XLALPrintProgressBar(1.0);
	XLALPrintInfo(" complete\n");

	/* output */

	XLALGPSTimeNow(&process_table_head->end_time);
	search_summary_head->nevents = XLALSimBurstAssignIDs(sim_burst_table_head, process_table_head->process_id, 0);
	write_xml(process_table_head, process_params_table_head, search_summary_head, sim_burst_table_head, options);

	/* done */

	gsl_rng_free(rng);
	XLALDestroyProcessTable(process_table_head);
	XLALDestroyProcessParamsTable(process_params_table_head);
	XLALDestroySearchSummaryTable(search_summary_head);
	XLALDestroySimBurstTable(sim_burst_table_head);
	exit(0);
}
