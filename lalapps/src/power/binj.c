/*----------------------------------------------------------------------- 
 * 
 * File Name: binj.c
 *
 * Author: Brady, P. R., Brown, D. A., Crieghton, J. D. E., Ray Majumder S,
 * Cannon K. C.
 *
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <regex.h>
#include <time.h>
#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/Random.h>
#include <lal/TimeDelay.h>
#include <lalapps.h>
#include <processtable.h>

int snprintf(char *str, size_t size, const char *format, ...);

RCSID("$Id$");

#define KPC ( 1e3 * LAL_PC_SI )
#define MPC ( 1e6 * LAL_PC_SI )
#define GPC ( 1e9 * LAL_PC_SI )

#define MW_CENTER_J2000_RA_RAD 2.0318570464121519l /* = 27940.04 s = 7 h 45 m 40.04 s */
#define MW_CENTER_J2000_DEC_RAD -0.50628171572274738l /* = -29o 00' 28.1" */


#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "binj"

#define TRUE      1
#define FALSE     0


/*
 * ============================================================================
 *                                Command Line
 * ============================================================================
 */

enum population {
	POPULATION_GALACTIC_CORE,
	POPULATION_UNIFORM_SKY,
	POPULATION_ZENITH
};

enum strain_dist {
	STRAIN_DIST_NONE,
	STRAIN_DIST_CONST_HPEAK,
	STRAIN_DIST_CONST_HRSS,
	STRAIN_DIST_LOG_HPEAK,
	STRAIN_DIST_LOG_HRSS,
	STRAIN_DIST_PRESET_HPEAK
};

enum freq_dist {
	FREQ_DIST_NONE,
	FREQ_DIST_STRING_CUSP,
	FREQ_DIST_MONO_ARITHMETIC,
	FREQ_DIST_MONO_GEOMETRIC,
	FREQ_DIST_RANDOM_GEOMETRIC
};

struct options {
	INT8 gps_start_time;
	INT8 gps_end_time;
	enum population population;
	enum freq_dist freq_dist;
	double flow;
	double fhigh;
	double fratio;
	double deltaf;
	double quality;
	enum strain_dist strain_dist;
	double strain_scale_min;
	double strain_scale_max;
	double log10_max_distance;
	double log10_min_distance;
	int seed;
	INT8 time_step;	/* nanoseconds between injections */
	char *waveform;
	double simwaveform_duration;
	int simwaveform_min_number;
	int simwaveform_max_number;
	double tau;
	char *user_tag;
	int mdc;	/* Need to set this to true if one wants to use MDC signals */
};


static struct options options_defaults(void)
{
	struct options defaults = {
		.gps_start_time = 0,
		.gps_end_time = 0,
		.population = POPULATION_UNIFORM_SKY,
		.freq_dist = FREQ_DIST_NONE,
		.flow = 150.0,
		.fhigh = 1000.0,
		.fratio = 0.0,
		.deltaf = 0.0,
		.quality = -1.0,
		.strain_dist = STRAIN_DIST_NONE,
		.strain_scale_min = 0.0,
		.strain_scale_max = 0.0,
		.log10_max_distance = log10(10000.0),
		.log10_min_distance = log10(100.0),
		.seed = 1,
		.time_step = 210.0 / LAL_PI * 1e9,
		.waveform = "SineGaussian",
		.simwaveform_duration = 0.0,
		.simwaveform_min_number = 0,
		.simwaveform_max_number = 10,
		.tau = 0.1,
		.user_tag = NULL,
		.mdc = FALSE
	};

	return defaults;
}


static void print_usage(const char *prog)
{
	fprintf(stderr, 
"%s [options]\n"\
"\nDefaults are shown in brackets\n\n"\
"\t--help                   display this message\n"\
"\t--gps-start-time SECONDS start injections at GPS time TIME\n"\
"\t--gps-end-time SECONDS   end injections at GPS time TIME\n"\
"\t--time-step STEP         space injections STEP / pi seconds appart (210)\n"\
"\t--population [galactic_core|uniform_sky|zenith] select the population to synthesize\n"\
"\t--freq-dist [monoarithmetic|monogeometric|randgeometric] select the frequency distribution\n"\
"\t--flow FLOW              first frequency of injection (150.0)\n"\
"\t--fhigh FHIGH            only inject frequencies smaller than FHIGH (1000.0)\n"\
"\t--fratio FACTOR          exponential spacing of injection frequencies (0.0)\n"\
"\t--deltaf DELF            linear spacing of injection frequencies (0.0)\n"\
"\t--quality Q              quality factor for SG waveforms TAU=Q/(sqrt(2) pi F)\n"\
"\t--tau TAU                duration of SG waveforms.  Q overrides TAU setting\n"\
"\t--strain-dist [consthpeak|consthrss|hpeakpresets|loghpeak|loghrss] select the strain distribution\n"\
"\t--strain-scale-min VALUE  mininum of the strain distribution\n"\
"\t--strain-scale-max VALUE  maximum of the strain distribution\n"\
"\t--min-distance VALUE     min distance of source in Kpc(default 100Kpc) \n"\
"\t--max-distance VALUE     max distance of source in Kpc(default 10000Kpc) \n"\
"\t--simwaveform-duration   duration of the simulated waveform (Warren/Ott/ZM)\n"\
"\t--simwaveform-min-number min # of the simulated waveform \n"\
"\t--simwaveform-max-number max # of the simulated waveform \n"\
"\t--seed SEED              seed random number generator with SEED (1)\n"\
"\t--waveform NAME          set waveform type to NAME (SineGaussian)\n"\
"\t--user-tag STRING        set the usertag to STRING\n\n", prog);
}


static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param, const char *type, const char *param, const char *value)
{
	*proc_param = LALCalloc(1, sizeof(**proc_param));
	(*proc_param)->next = NULL;
	snprintf((*proc_param)->program, LIGOMETA_PROGRAM_MAX, PROGRAM_NAME);
	snprintf((*proc_param)->type, LIGOMETA_TYPE_MAX, type);
	snprintf((*proc_param)->param, LIGOMETA_PARAM_MAX, "--%s", param);
	snprintf((*proc_param)->value, LIGOMETA_VALUE_MAX, value);

	return(&(*proc_param)->next);
}

	
#define ADD_PROCESS_PARAM(type) \
	do { paramaddpoint = add_process_param(paramaddpoint, type, long_options[option_index].name, optarg); } while(0)


static struct options parse_command_line(int *argc, char **argv[], MetadataTable *procparams)
{
	ProcessParamsTable **paramaddpoint = &procparams->processParamsTable;
	struct options options = options_defaults();
	int c;
	int option_index;
	struct option long_options[] = {
		{"mdc", no_argument, &options.mdc, TRUE},
		{"population", required_argument, NULL, 'A'},
		{"freq-dist", required_argument, NULL, 'F'},
		{"deltaf", required_argument, NULL, 'B'},
		{"fhigh", required_argument, NULL, 'C'},
		{"flow", required_argument, NULL, 'D'},
		{"fratio", required_argument, NULL, 'E'},
		{"gps-end-time", required_argument, NULL, 'G'},
		{"gps-start-time", required_argument, NULL, 'H'},
		{"help", no_argument, NULL, 'I'},
		{"strain-dist", required_argument, NULL, 'J'},
		{"strain-scale-max", required_argument, NULL, 'K'},
		{"strain-scale-min", required_argument, NULL, 'L'},
		{"max-distance", required_argument, NULL, 'M'},
		{"min-distance", required_argument, NULL, 'N'},
		{"quality", required_argument, NULL, 'O'},
		{"seed", required_argument, NULL, 'P'},
		{"simwaveform-duration", required_argument, NULL, 'Q'},
		{"simwaveform-max-number", required_argument, NULL, 'R'},
		{"simwaveform-min-number", required_argument, NULL, 'S'},
		{"tau", required_argument, NULL, 'T'},
		{"time-step", required_argument, NULL, 'U'},
		{"user-tag", required_argument, NULL, 'V'},
		{"waveform", required_argument, NULL, 'W'},
		{NULL, 0, NULL, 0}
	};

	do switch(c = getopt_long(*argc, *argv, "", long_options, &option_index)) {
	case 'A':
		if(!strcmp(optarg, "galactic_core"))
			options.population = POPULATION_GALACTIC_CORE;
		else if(!strcmp(optarg, "uniform_sky"))
			options.population = POPULATION_UNIFORM_SKY;
		else if(!strcmp(optarg, "zenith"))
			options.population = POPULATION_ZENITH;
		else {
			fprintf(stderr, "error: unrecognized population \"%s\"", optarg);
			exit(1);
		}
		ADD_PROCESS_PARAM("string");
		break;

	case 'F':
		if(!strcmp(optarg, "monoarithmetic"))
			options.freq_dist = FREQ_DIST_MONO_ARITHMETIC;
		else if(!strcmp(optarg, "monogeometric"))
			options.freq_dist = FREQ_DIST_MONO_GEOMETRIC;
		else if(!strcmp(optarg, "randgeometric"))
			options.freq_dist = FREQ_DIST_RANDOM_GEOMETRIC;
		else {
			fprintf(stderr, "error: unrecognized frequency distribution \"%s\"", optarg);
			exit(1);
		}
		ADD_PROCESS_PARAM("string");
		break;

	case 'B':
		options.deltaf = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'C':
		options.fhigh = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'D':
		options.flow = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'E':
		options.fratio = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'G':
		XLALClearErrno();
		{
			LIGOTimeGPS tmp;
			XLALStrToGPS(&tmp, optarg, NULL);
			options.gps_end_time = XLALGPSToINT8NS(&tmp);
		}
		if(xlalErrno || (options.gps_end_time < LAL_INT8_C(441417609000000000)) || (options.gps_end_time > LAL_INT8_C(999999999000000000))) {
			fprintf(stderr, "invalid --%s (%s specified)\n", long_options[option_index].name, optarg);
			exit(1);
		}
		ADD_PROCESS_PARAM("string");
		break;

	case 'H':
		XLALClearErrno();
		{
			LIGOTimeGPS tmp;
			XLALStrToGPS(&tmp, optarg, NULL);
			options.gps_start_time = XLALGPSToINT8NS(&tmp);
		}
		if(xlalErrno || (options.gps_start_time < LAL_INT8_C(441417609000000000)) || (options.gps_start_time > LAL_INT8_C(999999999000000000))) {
			fprintf(stderr, "invalid --%s (%s specified)\n", long_options[option_index].name, optarg);
			exit(1);
		}
		ADD_PROCESS_PARAM("string");
		break;

	case 'I':
		print_usage((*argv)[0]);
		exit(0);

	case 'J':
		if(!strcmp(optarg, "consthpeak"))
			options.strain_dist = STRAIN_DIST_CONST_HPEAK;
		else if(!strcmp(optarg, "consthrss"))
			options.strain_dist = STRAIN_DIST_CONST_HRSS;
		else if(!strcmp(optarg, "hpeakpresets"))
			options.strain_dist = STRAIN_DIST_PRESET_HPEAK;
		else if(!strcmp(optarg, "loghpeak"))
			options.strain_dist = STRAIN_DIST_LOG_HPEAK;
		else if(!strcmp(optarg, "loghrss"))
			options.strain_dist = STRAIN_DIST_LOG_HRSS;
		else {
			fprintf(stderr, "error: unrecognized strain distribution \"%s\"", optarg);
			exit(1);
		}
		ADD_PROCESS_PARAM("string");
		break;

	case 'K':
		options.strain_scale_max = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'L':
		options.strain_scale_min = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'M':
		options.log10_max_distance = log10(atof(optarg));
		ADD_PROCESS_PARAM("double");
		break;

	case 'N':
		options.log10_min_distance = log10(atof(optarg));
		ADD_PROCESS_PARAM("double");
		break;

	case 'O':
		options.quality = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'P':
		options.seed = atoi(optarg);
		ADD_PROCESS_PARAM("int");
		break;

	case 'Q':
		options.simwaveform_duration = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'R':
		options.simwaveform_max_number = atoi(optarg);
		ADD_PROCESS_PARAM("int");
		break;

	case 'S':
		options.simwaveform_min_number = atoi(optarg);
		ADD_PROCESS_PARAM("int");
		break;

	case 'T':
		options.tau = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'U':
		options.time_step = atof(optarg) / LAL_PI * 1e9;
		ADD_PROCESS_PARAM("double");
		break;

	case 'V':
		options.user_tag = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 'W':
		options.waveform = optarg;
		if(strlen(options.waveform) > LIGOMETA_WAVEFORM_MAX) {
			fprintf(stderr, "error: --waveform %s exceeds max length of %d characters\n", options.waveform, LIGOMETA_WAVEFORM_MAX);
			exit(1);
		}
		if(!strcmp(options.waveform, "StringCusp"))
			options.freq_dist = FREQ_DIST_STRING_CUSP;
		ADD_PROCESS_PARAM("string");
		break;

	case 0:
		/* option sets a flag */
		break;

	case -1:
		/* end of arguments */
		break;

	case '?':
		/* unrecognized option */
		print_usage((*argv)[0]);
		exit(1);

	case ':':
		/* missing argument for an option */
		print_usage((*argv)[0]);
		exit(1);
	} while(c != -1);


	/* check some of the input parameters for consistency */
	if(options.fhigh == options.flow) {
		/* use monotonic arithmetic distribution for case of fixed
		 * frequency */
		options.freq_dist = FREQ_DIST_MONO_ARITHMETIC;
		options.deltaf = 1.0;
	}
	if(options.fhigh < options.flow) {
		fprintf(stderr, "error: --fhigh < --flow\n");
		exit(1);
	}

	switch(options.freq_dist) {
	case FREQ_DIST_NONE:
		fprintf(stderr, "error: must select a frequency distribution\n");
		exit(1);

	case FREQ_DIST_STRING_CUSP:
		if(options.freq_dist != FREQ_DIST_STRING_CUSP) {
			fprintf(stderr, "error: cannot choose frequency distribution with --waveform=StringCusp\n");
			exit(1);
		}
		break;

	case FREQ_DIST_MONO_ARITHMETIC:
		if((options.deltaf <= 0.0) || (options.fratio != 0.0)) {
			fprintf(stderr, "error: invalid frequency distribution parameters\n");
			exit(1);
		}
		break;

	case FREQ_DIST_MONO_GEOMETRIC:
		if((options.deltaf != 0.0) || (options.fratio <= 0.0)) {
			fprintf(stderr, "error: invalid frequency distribution parameters\n");
			exit(1);
		}
		break;

	case FREQ_DIST_RANDOM_GEOMETRIC:
		if((options.deltaf != 0.0) || (options.fratio <= 0.0) || (options.fhigh / options.flow < sqrt(options.fratio))) {
			fprintf(stderr, "error: invalid frequency distribution parameters\n");
			exit(1);
		}
		break;
	}

	switch(options.strain_dist) {
	case STRAIN_DIST_NONE:
		fprintf(stderr, "error: must select a strain distribution\n");
		exit(1);

	case STRAIN_DIST_CONST_HPEAK:
	case STRAIN_DIST_CONST_HRSS:
		if(options.strain_scale_min == 0.0) {
			fprintf(stderr, "error: must set --strain-scale-min\n");
			exit(1);
		}
		if(options.strain_scale_max != 0.0)
			fprintf(stderr, "warning: --strain-scale-max provided but ignored\n");
		break;

	case STRAIN_DIST_LOG_HPEAK:
		if(options.strain_scale_min == 0.0) {
			fprintf(stderr, "error: must set --strain-scale-min\n");
			exit(1);
		}
		if(options.strain_scale_max == 0.0) {
			fprintf(stderr, "error: must set --strain-scale-max\n");
			exit(1);
		}
		break;

	case STRAIN_DIST_LOG_HRSS:
		if(options.strain_scale_min == 0.0) {
			fprintf(stderr, "error: must set --strain-scale-min\n");
			exit(1);
		}
		if(options.strain_scale_max == 0.0) {
			fprintf(stderr, "error: must set --strain-scale-max\n");
			exit(1);
		}
		break;

	case STRAIN_DIST_PRESET_HPEAK:
		if(options.strain_scale_min != 0.0)
			fprintf(stderr, "warning: --strain-scale-min provided but ignored\n");
		if(options.strain_scale_max != 0.0)
			fprintf(stderr, "warning: --strain-scale-max provided but ignored\n");
		break;
	}

	if(!options.gps_start_time || !options.gps_end_time) {
		fprintf(stderr, "--gps-start-time and --gps-end-time are both required\n");
		exit(1);
	}

	if((options.fratio != 0.0) && (options.deltaf != 0.0)) {
		fprintf(stderr, "error:  cannot specify both --deltaf and --fratio\n");
		exit(1);
	}

	return options;
}


/* 
 * ============================================================================
 *                            Solving For Unknowns
 * ============================================================================
 */

static double tau_from_q_and_f(double Q, double f)
{
	/* compute duration from Q and frequency */
	return Q / (sqrt(2.0) * LAL_PI * f);
}

static double hrss_from_tau_and_hpeak(double tau, double hpeak)
{
	/* compute root-sum-square strain from duration and peak strain */
	return sqrt(sqrt(2.0 * LAL_PI) * tau / 4.0) * hpeak;
}

static double hpeak_from_tau_and_hrss(double tau, double hrss)
{
	/* compute peak strain from duration and root-sum-square strain */
	return hrss / sqrt(sqrt(2.0 * LAL_PI) * tau / 4.0);
}


/* 
 * ============================================================================
 *                          Arrival Time Calculation
 * ============================================================================
 */

static LIGOTimeGPS equatorial_arrival_time(LALDetector detector, LIGOTimeGPS geocent_peak_time, double right_ascension, double declination)
{
	return *XLALGPSAdd(&geocent_peak_time, XLALTimeDelayFromEarthCenter(detector.location, right_ascension, declination, &geocent_peak_time));
}


static LIGOTimeGPS horizon_arrival_time(LALDetector detector, LIGOTimeGPS geocent_peak_time)
{
	const double r = sqrt(detector.location[0] * detector.location[0] + detector.location[1] * detector.location[1] + detector.location[2] * detector.location[2]);

	return *XLALGPSAdd(&geocent_peak_time, -r / LAL_C_SI);
}


/* 
 * ============================================================================
 *                Logarithmically-Distributed Random Variable
 * ============================================================================
 */

static double Log10Deviate(RandomParams *params, double log10_min, double log10_max)
{
	return pow(10.0, log10_min + (log10_max - log10_min) * XLALUniformDeviate(params));
}


/* 
 * ============================================================================
 *                          Frequency Distributions
 * ============================================================================
 */

/*
 * String cusps.  What's distributed uniformly is theta^2, the square of
 * the angle the line of sight makes with direction of the cusp.
 */

static double freq_dist_string_cusp_next(struct options options)
{
	const double thetasqmin = pow(options.fhigh, -2.0 / 3.0);
	const double thetasqmax = pow(options.flow, -2.0 / 3.0);
	const double thetasq = (thetasqmax - thetasqmin) * ((float) rand() / (float) RAND_MAX) + thetasqmin;

	return pow(thetasq, -3.0 / 2.0);
}


/*
 * Monotonically increasing arithmetic sequence.
 */

static double freq_dist_monotonic_arithmetic_next(struct options options)
{
	static int i = 0;
	const double freq = options.flow + options.deltaf * i++;

	if(freq <= options.fhigh)
		return freq;
	i = 1;
	return options.flow;
}


/*
 * Monotonically increasing geometric sequence.
 */

static double freq_dist_monotonic_geometric_next(struct options options)
{
	static int i = 0;
	const double freq = options.flow * pow(options.fratio, i++);

	if(freq <= options.fhigh)
		return freq;
	i = 1;
	return options.flow;
}


/*
 * Logarithmically-distributed random(ish) sequence.
 */

static double freq_dist_random_geometric_next(RandomParams *randparams, struct options options)
{
	const double sqrtratio = sqrt(options.fratio);
	double freq = freq_dist_monotonic_geometric_next(options) * sqrtratio;
	if(freq > options.fhigh)
		freq = freq_dist_monotonic_geometric_next(options) * sqrtratio;

	return Log10Deviate(randparams, log10(freq / sqrtratio), log10(freq * sqrtratio));
}


/* 
 * ============================================================================
 *                               h_peak Presets
 * ============================================================================
 */

static double hpeak_preset_next(RandomParams *params)
{
	static const double presets[] = {
		1e-19,
		2.0e-19,
		3.0e-19,
		3.3e-18,
		5.3e-17
	};
	int i = XLALUniformDeviate(params) * sizeof(presets)/sizeof(*presets);

	return presets[i];
}

/* 
 * ============================================================================
 *                                   Output
 * ============================================================================
 */

/* LIGO LW XML of MDC injections */
static void write_mdc_xml(MetadataTable mdcinjections)
{
	LALStatus status = blank_status;
	CHAR fname[256];
	LIGOLwXMLStream xmlfp;

	memset(&xmlfp, 0, sizeof(xmlfp));

	snprintf(fname, sizeof(fname), "HL-MDCSG10_%d.xml", 1);

	LAL_CALL(LALOpenLIGOLwXMLFile(&status, &xmlfp, fname), &status);
	LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, sim_burst_table), &status);
	LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, mdcinjections, sim_burst_table), &status);
	LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlfp), &status);
	LAL_CALL(LALCloseLIGOLwXMLFile(&status, &xmlfp), &status);
}


/* LIGO LW XML */
static void write_xml(MetadataTable proctable, MetadataTable procparams, MetadataTable injections, struct options options)
{
	LALStatus status = blank_status;
	char fname[256];
	LIGOLwXMLStream xmlfp;

	memset(&xmlfp, 0, sizeof(xmlfp));

	if(options.user_tag)
		snprintf(fname, sizeof(fname), "HL-INJECTIONS_%s-%d-%d.xml", options.user_tag, (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));
	else
		snprintf(fname, sizeof(fname), "HL-INJECTIONS-%d-%d.xml", (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));


	LAL_CALL(LALOpenLIGOLwXMLFile(&status, &xmlfp, fname), &status);

	XLALGPSTimeNow(&proctable.processTable->end_time);
	LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, process_table), &status);
	LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, proctable, process_table), &status);
	LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlfp), &status);

	if(procparams.processParamsTable) {
		LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, process_params_table), &status);
		LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, procparams, process_params_table), &status);
		LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlfp), &status);
	}

	if(injections.simBurstTable) {
		LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, sim_burst_table), &status);
		LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, injections, sim_burst_table), &status);
		LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlfp), &status);
	}

	LAL_CALL(LALCloseLIGOLwXMLFile(&status, &xmlfp), &status);
}


/* 
 * ============================================================================
 *                                Entry Point
 * ============================================================================
 */

int main(int argc, char *argv[])
{
	struct options options;
	INT8 tinj;
	RandomParams *randParams = NULL;
	LALStatus status = blank_status;
	LALDetector lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
	LALDetector llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
	MetadataTable proctable;
	MetadataTable procparams;
	MetadataTable injections;
	MetadataTable mdcinjections;
	SimBurstTable *this_sim_burst = NULL;

	/*
	 * Initialize debug handler
	 */

	lal_errhandler = LAL_ERR_EXIT;
	set_debug_level("LALMSGLVL2");

	/*
	 * Process params table and command line.
	 */

	procparams.processParamsTable = NULL;
	options = parse_command_line(&argc, &argv, &procparams);

	/*
	 * Process table
	 */

	proctable.processTable = calloc(1, sizeof(ProcessTable));
	XLALGPSTimeNow(&proctable.processTable->start_time);
	LAL_CALL(populate_process_table(&status, proctable.processTable, PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE), &status);
	if(options.user_tag)
		snprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, "%s", options.user_tag);

	/*
	 * Initialize random number generator
	 */

	randParams = XLALCreateRandomParams(options.seed);

	/*
	 * Main loop
	 */

	for(tinj = options.gps_start_time; tinj <= options.gps_end_time; tinj += options.time_step) {
		/* allocate the injection */
		if(this_sim_burst) {
			this_sim_burst->next = calloc(1, sizeof(SimBurstTable));
			this_sim_burst = this_sim_burst->next;
		} else
			this_sim_burst = injections.simBurstTable = calloc(1, sizeof(SimBurstTable));

		/* frequency */
		switch(options.freq_dist) {
		case FREQ_DIST_NONE:
			fprintf(stderr, "error: internal error\n");
			exit(1);

		case FREQ_DIST_STRING_CUSP:
			this_sim_burst->freq = freq_dist_string_cusp_next(options);
			break;

		case FREQ_DIST_MONO_ARITHMETIC:
			this_sim_burst->freq = freq_dist_monotonic_arithmetic_next(options);
			break;

		case FREQ_DIST_MONO_GEOMETRIC:
			this_sim_burst->freq = freq_dist_monotonic_geometric_next(options);
			break;

		case FREQ_DIST_RANDOM_GEOMETRIC:
			this_sim_burst->freq = freq_dist_random_geometric_next(randParams, options);
			break;
		}

		/* tau */
		if(options.quality > 0.0)
			this_sim_burst->tau = tau_from_q_and_f(options.quality, this_sim_burst->freq);
		else
			this_sim_burst->tau = options.tau;

		/* waveform */
		snprintf(this_sim_burst->waveform, LIGOMETA_WAVEFORM_MAX, "%s", options.waveform);

		/* peak time at geocentre in GPS seconds */
		XLALINT8NSToGPS(&this_sim_burst->geocent_peak_time, tinj);

		/* peak time at geocentre in Greenwich sidereal hours */
		this_sim_burst->peak_time_gmst = XLALGreenwichMeanSiderealTime(&this_sim_burst->geocent_peak_time) * 12.0 / LAL_PI;

		/* sky location and observatory peak times */
		switch(options.population) {
		case POPULATION_GALACTIC_CORE:
			/* co-ordinate system */
			sprintf(this_sim_burst->coordinates, "EQUATORIAL");

			/* right ascension, declination, and polarization */
			this_sim_burst->longitude = MW_CENTER_J2000_RA_RAD;
			this_sim_burst->latitude = MW_CENTER_J2000_DEC_RAD;
			this_sim_burst->polarization = 2.0 * LAL_PI * XLALUniformDeviate(randParams);

			/* observatory peak times in GPS seconds */
			this_sim_burst->h_peak_time = equatorial_arrival_time(lho, this_sim_burst->geocent_peak_time, this_sim_burst->longitude, this_sim_burst->latitude);
			this_sim_burst->l_peak_time = equatorial_arrival_time(llo, this_sim_burst->geocent_peak_time, this_sim_burst->longitude, this_sim_burst->latitude);
			break;

		case POPULATION_UNIFORM_SKY:
			/* co-ordinate system */
			sprintf(this_sim_burst->coordinates, "EQUATORIAL");

			/* right ascension, declination, and polarization */
			this_sim_burst->longitude = 2.0 * LAL_PI * XLALUniformDeviate(randParams);
			this_sim_burst->latitude = asin(2.0 * XLALUniformDeviate(randParams) - 1.0);
			this_sim_burst->polarization = 2.0 * LAL_PI * XLALUniformDeviate(randParams);

			/* observatory peak times in GPS seconds */
			this_sim_burst->h_peak_time = equatorial_arrival_time(lho, this_sim_burst->geocent_peak_time, this_sim_burst->longitude, this_sim_burst->latitude);
			this_sim_burst->l_peak_time = equatorial_arrival_time(llo, this_sim_burst->geocent_peak_time, this_sim_burst->longitude, this_sim_burst->latitude);
			break;

		case POPULATION_ZENITH:
			/* co-ordinate system */
			sprintf(this_sim_burst->coordinates, "ZENITH");

			/* optimally oriented overhead */
			this_sim_burst->longitude = 0.0;
			this_sim_burst->latitude = 0.0;
			this_sim_burst->polarization = 0.0;

			/* observatory peak times in GPS seconds */
			this_sim_burst->h_peak_time = horizon_arrival_time(lho, this_sim_burst->geocent_peak_time);
			this_sim_burst->l_peak_time = horizon_arrival_time(llo, this_sim_burst->geocent_peak_time);
			break;
		}

		/* strain */
		switch(options.strain_dist) {
		case STRAIN_DIST_NONE:
			fprintf(stderr, "error: internal error\n");
			exit(1);

		case STRAIN_DIST_CONST_HPEAK:
			this_sim_burst->hpeak = options.strain_scale_min;
			this_sim_burst->hrss = hrss_from_tau_and_hpeak(this_sim_burst->tau, this_sim_burst->hpeak);
			break;

		case STRAIN_DIST_CONST_HRSS:
			this_sim_burst->hrss = options.strain_scale_min;
			this_sim_burst->hpeak = hpeak_from_tau_and_hrss(this_sim_burst->tau, this_sim_burst->hrss);
			break;

		case STRAIN_DIST_LOG_HPEAK:
			this_sim_burst->hpeak = Log10Deviate(randParams, options.strain_scale_min, options.strain_scale_max);
			this_sim_burst->hrss = hrss_from_tau_and_hpeak(this_sim_burst->tau, this_sim_burst->hpeak);
			break;

		case STRAIN_DIST_LOG_HRSS:
			this_sim_burst->hrss = Log10Deviate(randParams, options.strain_scale_min, options.strain_scale_max);
			this_sim_burst->hpeak = hpeak_from_tau_and_hrss(this_sim_burst->tau, this_sim_burst->hrss);
			break;

		case STRAIN_DIST_PRESET_HPEAK:
			this_sim_burst->hpeak = hpeak_preset_next(randParams);
			this_sim_burst->hrss = hrss_from_tau_and_hpeak(this_sim_burst->tau, this_sim_burst->hpeak);
			break;
		}

		/* uniform distribution in log(distance) */
		this_sim_burst->distance = Log10Deviate(randParams, options.log10_min_distance, options.log10_max_distance);

		/* dt */
		if(options.simwaveform_duration)
			this_sim_burst->dtplus = this_sim_burst->dtminus = options.simwaveform_duration / 2.0;
		else
			this_sim_burst->dtplus = this_sim_burst->dtminus = 4.0 * this_sim_burst->tau;

		/* set the simulated wavenumber */
		this_sim_burst->zm_number = options.simwaveform_min_number + (options.simwaveform_max_number - options.simwaveform_min_number) * XLALUniformDeviate(randParams);
	}

	/* if using the MDC signal frames for injections:
	 * then read in the ascii file of inj. parameters
	 * and write a sim burst table out of that. However
	 * currently all the fields are not filled.[20040430: SKRM]
	 */
	if(options.mdc) {
		INT4 n = 0;
		INT4 nn, x, rc;
		CHAR mdcSimfile[20];
		FILE *fp;
		float *fHp, *fHx, *fLp, *fLx, *theta, *phi, *psi, *hrss;
		float *gg, *hh, *ii, *jj, *kk, *ll, *mm, *rr;
		int *gps, *tH, *tL, *tauHL;
		int *cc, *dd, *ee, *ff;
		char **waveform;
		char *qq, *y;

		SimBurstTable *this_mdcsim_burst = NULL;

		/*read in the ascii file containing the injection parameters 
		 *the ascii file is assumed to be present in the working dir. 
		 */
		snprintf(mdcSimfile, 20, "mdcsim_%d.dat", 1);
		fp = fopen(mdcSimfile, "r");
		if(fp == NULL) {
			fprintf(stdout, "Error:Must have file mdcsim_1.dat in the working dir.\n");
			exit(1);
		}
		while(fscanf(fp, "%*d %*d %*d %*d %*e %*e %*e %*e %*e %*e %*e %*s %*e") != EOF)
			++n;
		rewind(fp);

		gps = cc = (int *) malloc(n * sizeof(int));
		tL = dd = (int *) malloc(n * sizeof(int));
		tH = ee = (int *) malloc(n * sizeof(int));
		tauHL = ff = (int *) malloc(n * sizeof(int));
		fHp = gg = (float *) malloc(n * sizeof(float));
		fHx = hh = (float *) malloc(n * sizeof(float));
		fLp = ii = (float *) malloc(n * sizeof(float));
		fLx = jj = (float *) malloc(n * sizeof(float));
		theta = kk = (float *) malloc(n * sizeof(float));
		phi = ll = (float *) malloc(n * sizeof(float));
		psi = mm = (float *) malloc(n * sizeof(float));
		hrss = rr = (float *) malloc(n * sizeof(float));
		qq = (char *) malloc((n * 10) * sizeof(char));	/*memory allocation may need to be 
								 *changed if the name of waveform
								 *changes in the mdc file 
								 */
		waveform = (char **) malloc(n * 10 * sizeof(char));

		/* reads in the diff. parameters from the text file */
		nn = 0;
		while((rc = fscanf(fp, "%d %d %d %d %e %e %e %e %e %e %e %s %e", cc++, dd++, ee++, ff++, gg++, hh++, ii++, jj++, kk++, ll++, mm++, qq, rr++)) == 13) {
			waveform[nn++] = qq;	/* scans the name of the waveform */
			qq = qq + 10;
			continue;
		}

		if(rc != EOF)
			printf("ERROR READING FILE\n");

		fclose(fp);

		/* create the first injection */
		this_mdcsim_burst = mdcinjections.simBurstTable = (SimBurstTable *)
		    calloc(1, sizeof(SimBurstTable));
		/*create the corresponding simburst table */
		for(x = 0; x < n; x++) {
			this_mdcsim_burst->geocent_peak_time.gpsSeconds = gps[x];	/* wrong value */
			this_mdcsim_burst->geocent_peak_time.gpsNanoSeconds = 0;	/* wrong value */
			this_mdcsim_burst->h_peak_time.gpsSeconds = gps[x];
			this_mdcsim_burst->h_peak_time.gpsNanoSeconds = tH[x] + (0.5 * 1000000000LL);
			this_mdcsim_burst->l_peak_time.gpsSeconds = gps[x];
			this_mdcsim_burst->l_peak_time.gpsNanoSeconds = tL[x] + (0.5 * 1000000000LL);
			this_mdcsim_burst->longitude = phi[x];
			this_mdcsim_burst->latitude = theta[x];
			this_mdcsim_burst->polarization = psi[x];
			this_mdcsim_burst->hrss = hrss[x];

			/*rip off the frequency from the name of the waveform */
			waveform[x] = waveform[x] + 2;
			y = (char *) malloc(10 * sizeof(char));
			strncpy(y, waveform[x], 3);
			this_mdcsim_burst->freq = atof(y);	/*freq of the SineGaussian */
			free(y);
			this_mdcsim_burst = this_mdcsim_burst->next = (SimBurstTable *)
			    calloc(1, sizeof(SimBurstTable));
		}


		free(gps);
		free(tH);
		free(tL);
		free(tauHL);
		free(fHp);
		free(fHx);
		free(fLp);
		free(fLx);
		free(theta);
		free(phi);
		free(psi);
		free(hrss);
	}

	/* output */
	if(options.mdc) {
		/* this is used when mdc frames are used */
		write_mdc_xml(mdcinjections);
	} else {
		/* non-mdc XML output */
		write_xml(proctable, procparams, injections, options);
	}

	exit(0);
}
