/*
 * Copyright (C) 2007 Denny Mackin, Duncan Brown, Ik Siong Heng, Jolien
 * Creighton, Kipp Cannon, Mark Williamsen, Patrick Brady, Robert Adam
 * Mercer, Saikat Ray-Majumder, Stephen Fairhurst
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


#include <getopt.h>
#include <lalapps.h>
#include <math.h>
#include <processtable.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/EPSearch.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpChisq.h>
#include <lal/FrameCalibration.h>
#include <lal/FrameStream.h>
#include <lal/FrequencySeries.h>
#include <lal/GenerateBurst.h>
#include <lal/Inject.h>
#include <lal/IIRFilter.h>
#include <lal/LALConstants.h>
#include <lal/LALDatatypes.h>
#include <lal/LALError.h>
#include <lal/LALFrameIO.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLArray.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/Random.h>
#include <lal/RealFFT.h>
#include <lal/ResampleTimeSeries.h>
#include <lal/SeqFactories.h>
#include <lal/SimulateCoherentGW.h>
#include <lal/TFTransform.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>
#include <lal/VectorOps.h>
#include <lal/Window.h>


/* ARGH!  allow the code to be C99!  Obsession with C89 will cause bugs */
double trunc(double);
int snprintf(char *str, size_t size, const char *format, ...);
int vsnprintf(char *str, size_t size, const char *format, va_list ap);


NRCSID(POWERC, "power $Id$");
RCSID("power $Id$");


#define PROGRAM_NAME "lalapps_power"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"


#define TRUE       1
#define FALSE      0


/*
 * ============================================================================
 *
 *                               Misc Utilities
 *
 * ============================================================================
 */


/*
 * Return the smaller of two size_t variables.
 */


static size_t min(size_t a, size_t b)
{
	return a < b ? a : b;
}


/*
 * Return TRUE if the given integer is an integer power of 2.  The trick
 * here is that integer powers of 2 (and only integer powers of 2) share
 * exactly 0 bits with the integer 1 less than themselves, so we check to
 * see if that's the case.
 */


static int is_power_of_2(int x)
{
	return !((x - 1) & x);
}


/*
 * Return TRUE if the given double is an integer power of 2.
 */


static int double_is_power_of_2(double x)
{
	if(x <= 0)
		return 0;
	if(x < 1)
		/* might be a -ve power */
		return double_is_power_of_2(1 / x);
	if(x == trunc(x))
		/* is an integer */
		return is_power_of_2(x);
	return 0;
}


/*
 * ============================================================================
 *
 *                       Initialize and parse arguments
 *
 * ============================================================================
 */


/*
 * Parameters from command line
 */


struct options {
	/*
	 * search parameters
	 */

	/* number of samples to use for PSD */
	int psd_length;
	/* number of samples separating PSD starts */
	int psd_shift;
	/* number of samples not tiled at the start of each analysis window
	 * */
	int window_pad;
	/* time-frequency plane stride */
	int window_shift;
	/* random number generator seed */
	unsigned long seed;

	/*
	 * frame data input
	 */

	LIGOTimeGPS gps_start;
	LIGOTimeGPS gps_end;
	/* name of file with frame cache info */
	char *cache_filename;
	/* name of the calibration cache file */
	char *cal_cache_filename;
	/* channel to analyze */
	char *channel_name;
	/* two character interferometer */
	char ifo[3];
	/* don't read a time series longer than this (e.g., because of
	 * available memory) */
	int max_series_length;

	/*
	 * h(t) support
	 */

	/* flag indicating that input is double-precision h(t) */
	int calibrated;

	/*
	 * data conditioning
	 */

	/* sample rate after resampling */
	int resample_rate;
	/* conditioning high pass freq (Hz) */
	double high_pass;
	/* samples corrupted by conditioning */
	int filter_corruption;

	/*
	 * injection support
	 */

	/* Gaussian white noise RMS */
	double noise_rms;
	/* name of mdc signal cache file */
	char *mdc_cache_filename;
	/* mdc signal only channnel info */
	char *mdc_channel_name;
	/* file with list of burst injections */
	char *sim_burst_filename;
	/* file with list of inspiral injections */
	char *sim_inspiral_filename;

	/*
	 * XLALEPSearch() parmaters
	 */

	double confidence_threshold;
	double bandwidth;
	double flow;
	double maxTileBandwidth;
	double maxTileDuration;
	double fractional_stride;
	REAL8Window *window;

	/*
	 * output control
	 */

	/* user comment */
	char *comment;
	/* safety valve (Hz), 0 == disable */
	int max_event_rate;
	char *output_filename;

	/*
	 * diagnostics support
	 */

	/* diagnostics call back for passing to LAL (NULL = disable) */
	struct XLALEPSearchDiagnostics *diagnostics;
};


static struct options *options_new(void)
{
	static char default_comment[] = "";
	struct options *options = malloc(sizeof(*options));

	if(!options)
		return NULL;

	*options = (struct options) {
		.cal_cache_filename = NULL,	/* default == disable */
		.channel_name = NULL,	/* impossible */
		.comment = default_comment,	/* default = "" */
		.filter_corruption = -1,	/* impossible */
		.mdc_cache_filename = NULL,	/* default == disable */
		.mdc_channel_name = NULL,	/* impossible */
		.noise_rms = -1.0,	/* default == disable */
		.diagnostics = NULL,	/* default == disable */
		.psd_length = 0,	/* impossible */
		.psd_shift = 0,	/* impossible */
		.resample_rate = 0,	/* impossible */
		.seed = 0,	/* default == use system clock */
		.max_series_length = 0,	/* default == disable */
		.calibrated = FALSE,	/* default */
		.high_pass = -1.0,	/* impossible */
		.max_event_rate = 0,	/* default == disable */
		.output_filename = NULL,	/* impossible */
		.sim_burst_filename = NULL,	/* default == disable */
		.sim_inspiral_filename = NULL,	/* default == disable */
		.cache_filename = NULL,	/* default == disable */
		.confidence_threshold = XLAL_REAL8_FAIL_NAN,	/* impossible */
		.bandwidth = 0,	/* impossible */
		.flow = -1,	/* impossible */
		.maxTileBandwidth = 0,	/* impossible */
		.maxTileDuration = 0,	/* impossible */
		.fractional_stride = 0,	/* impossible */
		.window = NULL,	/* impossible */
	};

	memset(options->ifo, 0, sizeof(options->ifo));	/* default = "" */
	XLALINT8NSToGPS(&options->gps_start, 0);	/* impossible */
	XLALINT8NSToGPS(&options->gps_end, 0);	/* impossible */

	return options;
}


static void options_free(struct options *options)
{
	if(options) {
		XLALDestroyREAL8Window(options->window);
	}
	free(options);
}


/*
 * Friendly messages.
 */


static void print_usage(const char *program)
{
	fprintf(stderr,
"Usage:  %s [option ...]\n" \
"The following options are recognized.  Options not surrounded in [] are\n" \
"required.\n" \
"	 --bandwidth <Hz>\n" \
"	[--burstinjection-file <XML file name>]\n" \
"	[--calibrated-data]\n" \
"	[--calibration-cache <cache file name>]\n" \
"	 --channel-name <name>\n" \
"	 --confidence-threshold <confidence>\n" \
"	[--debug-level info|warn|error|off]\n" \
"	[--dump-diagnostics <XML file name>]\n" \
"	 --filter-corruption <samples>\n" \
"	 --frame-cache <cache file name>\n" \
"	[--gaussian-noise-rms <RMS amplitude>]\n", program);
	fprintf(stderr,
"	 --gps-end-time <seconds>\n" \
"	 --gps-start-time <seconds>\n" \
"	[--help]\n" \
"	 --high-pass <Hz>\n" \
"	[--inspiralinjection-file <XML file name>]\n" \
"	 --low-freq-cutoff <Hz>\n" \
"	[--max-event-rate <Hz>]\n" \
"	 --max-tile-bandwidth <Hz>\n" \
"	 --max-tile-duration <seconds>\n" \
"	[--mdc-cache <cache file name>]\n" \
"	[--mdc-channel <name>]\n" \
"	[--output <filename>]\n" \
"	 --psd-average-points <samples>\n" \
"	[--ram-limit <MebiBytes>]\n" \
"	 --resample-rate <Hz>\n" \
"	[--seed <seed>]\n");
	fprintf(stderr,
"	 --tile-stride-fraction <fraction>\n" \
"	[--user-tag <comment>]\n" \
"	 --window-length <samples>\n");
}


static void print_bad_argument(const char *prog, const char *arg, const char *msg)
{
	XLALPrintError("%s: error: invalid argument for --%s: %s\n", prog, arg, msg);
}


static void print_missing_argument(const char *prog, const char *arg)
{
	XLALPrintError("%s: error: --%s not specified\n", prog, arg);
}


/*
 * Check for missing command line arguments.
 */


static int all_required_arguments_present(char *prog, struct option *long_options, const struct options *options)
{
	int index;
	int got_all_arguments = TRUE;
	int arg_is_missing;

	for(index = 0; long_options[index].name; index++) {
		switch(long_options[index].val) {
		case 'A':
			arg_is_missing = !options->bandwidth;
			break;

		case 'C':
			arg_is_missing = !options->channel_name;
			break;

		case 'K':
			arg_is_missing = !XLALGPSToINT8NS(&options->gps_end);
			break;

		case 'M':
			arg_is_missing = !XLALGPSToINT8NS(&options->gps_start);
			break;

		case 'Q':
			arg_is_missing = options->flow < 0.0;
			break;

		case 'W':
			arg_is_missing = !options->window;
			break;

		case 'Z':
			arg_is_missing = !options->psd_length;
			break;

		case 'e':
			arg_is_missing = !options->resample_rate;
			break;

		case 'f':
			arg_is_missing = !options->fractional_stride;
			break;

		case 'g':
			arg_is_missing = options->confidence_threshold == XLAL_REAL8_FAIL_NAN;
			break;

		case 'j':
			arg_is_missing = options->filter_corruption < 0;
			break;

		case 'l':
			arg_is_missing = !options->maxTileBandwidth;
			break;

		case 'm':
			arg_is_missing = !options->maxTileDuration;
			break;

		case 'o':
			arg_is_missing = options->high_pass < 0.0;
			break;

		default:
			arg_is_missing = FALSE;
			break;
		}
		if(arg_is_missing) {
			print_missing_argument(prog, long_options[index].name);
			got_all_arguments = FALSE;
		}
	}

	if(!!options->cache_filename + (options->noise_rms > 0.0) != 1) {
		XLALPrintError("%s: must provide exactly one of --frame-cache or --gaussian-noise-rms\n", prog);
		got_all_arguments = FALSE;
	}

	return got_all_arguments;
}


/*
 * Helper function for adding an entry to the process params table.
 */


static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param, const char *type, const char *param, const char *value)
{
	*proc_param = XLALCalloc(1, sizeof(**proc_param));
	(*proc_param)->next = NULL;
	snprintf((*proc_param)->program, LIGOMETA_PROGRAM_MAX, PROGRAM_NAME);
	snprintf((*proc_param)->type, LIGOMETA_TYPE_MAX, type);
	snprintf((*proc_param)->param, LIGOMETA_PARAM_MAX, "--%s", param);
	snprintf((*proc_param)->value, LIGOMETA_VALUE_MAX, "%s", value ? value : "");

	return &(*proc_param)->next;
}


#define ADD_PROCESS_PARAM(type) \
	do { paramaddpoint = add_process_param(paramaddpoint, type, long_options[option_index].name, optarg); } while(0)


/*
 * Parse the --debug-level command line argument.
 */


static int parse_command_line_debug(int argc, char *argv[])
{
	char msg[240];
	int args_are_bad = FALSE;
	int c;
	int option_index;
	struct option long_options[] = {
		{"debug-level", required_argument, NULL, 'D'},
		{NULL, 0, NULL, 0}
	};

	/*
	 * Default == print only error messages.
	 */

	lalDebugLevel = LALERROR | LALNMEMDBG | LALNMEMPAD | LALNMEMTRK;

	/*
	 * Find and parse only the debug level command line options.  Must
	 * jump through this hoop because we cannot edit lalDebugLevel
	 * after any calls to XLALMalloc() and friends.
	 */

	opterr = 0;		/* silence error messages */
	optind = 0;		/* start scanning from argv[0] */
	do switch(c = getopt_long(argc, argv, "-", long_options, &option_index)) {
	case 'D':
		if(!strcmp(optarg, "info"))
			lalDebugLevel = LALINFO | LALWARNING | LALERROR | LALNMEMDBG | LALNMEMPAD | LALNMEMTRK;
		else if(!strcmp(optarg, "warn"))
			lalDebugLevel = LALWARNING | LALERROR | LALNMEMDBG | LALNMEMPAD | LALNMEMTRK;
		else if(!strcmp(optarg, "error"))
			lalDebugLevel = LALERROR | LALNMEMDBG | LALNMEMPAD | LALNMEMTRK;
		else if(!strcmp(optarg, "off"))
			lalDebugLevel = LALNMEMDBG | LALNMEMPAD | LALNMEMTRK;
		else {
			sprintf(msg, "must be one of \"info\", \"warn\", \"error\", or \"off\"");
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		break;

	default:
		break;
	} while(c != -1);

	return args_are_bad ? -1 : 0;
}


/*
 * Parse all the other command line arguments.
 */


static struct options *parse_command_line(int argc, char *argv[], MetadataTable *_process_params_table)
{
	struct options *options;
	char msg[240];
	int args_are_bad = FALSE;
	int c;
	int option_index;
	ProcessParamsTable **paramaddpoint = &_process_params_table->processParamsTable;
	struct option long_options[] = {
		{"bandwidth", required_argument, NULL, 'A'},
		{"burstinjection-file", required_argument, NULL, 'P'},
		{"calibrated-data", no_argument, NULL, 'J'},
		{"calibration-cache", required_argument, NULL, 'B'},
		{"channel-name", required_argument, NULL, 'C'},
		{"confidence-threshold", required_argument, NULL, 'g'},
		{"debug-level", required_argument, NULL, 'D'},
		{"dump-diagnostics", required_argument, NULL, 'X'},
		{"filter-corruption", required_argument, NULL, 'j'},
		{"frame-cache", required_argument, NULL, 'G'},
		{"gaussian-noise-rms", required_argument, NULL, 'V'},
		{"gps-end-time", required_argument, NULL, 'K'},
		{"gps-start-time", required_argument, NULL, 'M'},
		{"help", no_argument, NULL, 'O'},
		{"high-pass", required_argument, NULL, 'o'},
		{"inspiralinjection-file", required_argument, NULL, 'I'},
		{"low-freq-cutoff", required_argument, NULL, 'Q'},
		{"max-event-rate", required_argument, NULL, 'F'},
		{"max-tile-bandwidth", required_argument, NULL, 'l'},
		{"max-tile-duration", required_argument, NULL, 'm'},
		{"mdc-cache", required_argument, NULL, 'R'},
		{"mdc-channel", required_argument, NULL, 'S'},
		{"output", required_argument, NULL, 'b'},
		{"psd-average-points", required_argument, NULL, 'Z'},
		{"ram-limit", required_argument, NULL, 'a'},
		{"resample-rate", required_argument, NULL, 'e'},
		{"seed", required_argument, NULL, 'c'},
		{"tile-stride-fraction", required_argument, NULL, 'f'},
		{"user-tag", required_argument, NULL, 'h'},
		{"window-length", required_argument, NULL, 'W'},
		{NULL, 0, NULL, 0}
	};

	/*
	 * Allocate and initialize options structure
	 */

	options = options_new();
	if(!options)
		return NULL;

	/*
	 * Parse command line.
	 */

	opterr = 1;		/* enable error messages */
	optind = 0;		/* start scanning from argv[0] */
	do switch (c = getopt_long(argc, argv, "", long_options, &option_index)) {
	case 'A':
		options->bandwidth = atof(optarg);
		if(options->bandwidth <= 0 || !double_is_power_of_2(options->bandwidth)) {
			sprintf(msg, "must be greater than 0 and a power of 2 (%g specified)", options->bandwidth);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

	case 'B':
		options->cal_cache_filename = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 'C':
		options->channel_name = optarg;
		memcpy(options->ifo, optarg, sizeof(options->ifo) - 1);
		ADD_PROCESS_PARAM("string");
		break;

	case 'D':
		/* only add --debug-level to params table in this pass */
		ADD_PROCESS_PARAM("string");
		break;

	case 'F':
		options->max_event_rate = atoi(optarg);
		if(options->max_event_rate < 0) {
			sprintf(msg, "must not be negative (%d specified)", options->max_event_rate);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

	case 'G':
		options->cache_filename = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 'I':
		options->sim_inspiral_filename = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 'J':
		options->calibrated = TRUE;
		ADD_PROCESS_PARAM("string");
		break;

	case 'K':
		if(XLALStrToGPS(&options->gps_end, optarg, NULL)) {
			sprintf(msg, "range error parsing \"%s\"", optarg);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("string");
		break;

	case 'M':
		if(XLALStrToGPS(&options->gps_start, optarg, NULL)) {
			sprintf(msg, "range error parsing \"%s\"", optarg);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("string");
		break;

	case 'O':
		print_usage(argv[0]);
		exit(0);
		break;

	case 'P':
		options->sim_burst_filename = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 'Q':
		options->flow = atof(optarg);
		if(options->flow < 0) {
			sprintf(msg, "must not be negative (%f Hz specified)", options->flow);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

	case 'R':
		options->mdc_cache_filename = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 'S':
		options->mdc_channel_name = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 'V':
		options->noise_rms = atof(optarg);
		if(options->noise_rms <= 0.0) {
			sprintf(msg, "must be greater than 0 (%f specified)", options->noise_rms);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

	case 'W':
		{
		int window_length = atoi(optarg);
		if((window_length < 4) || !is_power_of_2(window_length)) {
			sprintf(msg, "must be greater than or equal to 4 and a power of 2 (%i specified)", window_length);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		options->window = XLALCreateHannREAL8Window(window_length);
		if(!options->window) {
			sprintf(msg, "failure generating %d sample Hann window", window_length);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
			/* silence "missing required argument" message */
			options->window = (void *) 1;
		}
		ADD_PROCESS_PARAM("int");
		}
		break;

	case 'X':
		options->diagnostics = malloc(sizeof(*options->diagnostics));
		{
		LALStatus stat;
		memset(&stat, 0, sizeof(stat));
		options->diagnostics->LIGOLwXMLStream = calloc(1, sizeof(*options->diagnostics));
		LALOpenLIGOLwXMLFile(&stat, options->diagnostics->LIGOLwXMLStream, optarg);
		}
		options->diagnostics->XLALWriteLIGOLwXMLArrayREAL8FrequencySeries = XLALWriteLIGOLwXMLArrayREAL8FrequencySeries;
		options->diagnostics->XLALWriteLIGOLwXMLArrayREAL8TimeSeries = XLALWriteLIGOLwXMLArrayREAL8TimeSeries;
		options->diagnostics->XLALWriteLIGOLwXMLArrayCOMPLEX16FrequencySeries = XLALWriteLIGOLwXMLArrayCOMPLEX16FrequencySeries;
		break;

	case 'Z':
		options->psd_length = atoi(optarg);
		ADD_PROCESS_PARAM("int");
		break;

	case 'a':
		/*
		 * Convert the available RAM (in MebiBytes) to a
		 * guestimated limit on the length of a time series to read
		 * in.
		 */
		options->max_series_length = atoi(optarg) * 1024 * 1024 / (8 * sizeof(REAL4));
		if(options->max_series_length <= 0) {
			sprintf(msg, "must be greater than 0 (%i specified)", atoi(optarg));
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

	case 'b':
		options->output_filename = optarg;
		ADD_PROCESS_PARAM("lstring");
		break;

	case 'c':
		options->seed = atol(optarg);
		ADD_PROCESS_PARAM("long");
		break;

	case 'e':
		options->resample_rate = atoi(optarg);
		if(options->resample_rate < 2 || options->resample_rate > 16384 || !is_power_of_2(options->resample_rate)) {
			sprintf(msg, "must be a power of 2 in the rage [2,16384] (%d specified)", options->resample_rate);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

	case 'f':
		options->fractional_stride = atof(optarg);
		if(options->fractional_stride > 1 || !double_is_power_of_2(options->fractional_stride)) {
			sprintf(msg, "must be less than or equal to 1 and a power of 2 (%g specified)", options->fractional_stride);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

	case 'g':
		options->confidence_threshold = atof(optarg);
		if(options->confidence_threshold < 0) {
			sprintf(msg, "must not be negative (%g specified)", options->confidence_threshold);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

	case 'h':
		options->comment = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 'j':
		options->filter_corruption = atoi(optarg);
		if(options->filter_corruption < 0) {
			sprintf(msg, "must not be negative (%d specified)", options->filter_corruption);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

	case 'l':
		options->maxTileBandwidth = atof(optarg);
		if((options->maxTileBandwidth <= 0) || !double_is_power_of_2(options->maxTileBandwidth)) {
			sprintf(msg, "must be greater than 0 and a power of 2 (%f specified)", options->maxTileBandwidth);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

	case 'm':
		options->maxTileDuration = atof(optarg);
		if((options->maxTileDuration <= 0) || !double_is_power_of_2(options->maxTileDuration)) {
			sprintf(msg, "must be greater than 0 and a power of 2 (%f specified)", options->maxTileDuration);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

	case 'o':
		options->high_pass = atof(optarg);
		if(options->high_pass < 0) {
			sprintf(msg, "must not be negative (%f Hz specified)", options->high_pass);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

	/* option sets a flag */
	case 0:
		optarg = NULL;
		ADD_PROCESS_PARAM("string");
		break;

	/* end of arguments */
	case -1:
		break;

	/* unrecognized option */
	case '?':
		print_usage(argv[0]);
		exit(1);
		break;

	/* missing argument for an option */
	case ':':
		print_usage(argv[0]);
		exit(1);
		break;
	} while(c != -1);

	/*
	 * Check for missing command line arguments.
	 */

	if(!all_required_arguments_present(argv[0], long_options, options))
		args_are_bad = TRUE;

	/*
	 * Check the order of the start and stop times.
	 */

	if(XLALGPSCmp(&options->gps_start, &options->gps_end) > 0) {
		XLALPrintError("%s: error: GPS start time > GPS stop time\n", argv[0]);
		args_are_bad = TRUE;
	}

	/*
	 * Warn if calibration cache is not provided that a unit response
	 * will be used for injections and h_rss.
	 */

	if(!options->cal_cache_filename) {
		XLALPrintWarning("warning: no calibration cache is provided:  software injections and hrss will be computed with unit response\n");
	} else if(options->calibrated) {
		XLALPrintError("error: calibration cache provided for use with calibrated data!\n");
		args_are_bad = TRUE;
	}

	/*
	 * Exit if anything was wrong with the command line.
	 */

	if(args_are_bad)
		exit(1);

	/*
	 * Compute timing parameters.
	 */

	if(XLALEPGetTimingParameters(options->window->data->length, options->maxTileDuration * options->resample_rate, options->fractional_stride, &options->psd_length, &options->psd_shift, &options->window_shift, &options->window_pad, NULL) < 0) {
		XLALPrintError("calculation of timing parameters failed\n");
		exit(1);
	}

	/*
	 * Ensure RAM limit is comensurate with the psd_length and its
	 * shift.
	 */

	if(options->max_series_length)
		options->max_series_length = XLALOverlappedSegmentsCommensurate(options->max_series_length, options->psd_length, options->psd_shift);

	/*
	 * Sanitize filter frequencies.
	 */

	if(options->high_pass > options->flow - 10.0)
		XLALPrintWarning("%s: warning: data conditioning high-pass frequency (%f Hz) greater than 10 Hz below TF plane low frequency (%f Hz)\n", argv[0], options->high_pass, options->flow);

	/*
	 * Set output filename to default value if needed.
	 */

	if(!options->output_filename) {
		int length = 1000;
		options->output_filename = calloc(length, sizeof(*options->output_filename));
		if(!options->output_filename) {
			XLALPrintError("memory error");
			exit(1);
		}
		length = snprintf(options->output_filename, length, "%s-POWER_%s-%d-%d.xml", options->ifo, options->comment, options->gps_start.gpsSeconds, options->gps_end.gpsSeconds - options->gps_start.gpsSeconds);
		if(length > 999) {
			XLALPrintError("output filename too long (999 characters max)");
			exit(1);
		}
		options->output_filename = realloc(options->output_filename, (length + 1) * sizeof(*options->output_filename));
		if(!options->output_filename) {
			XLALPrintError("memory error");
			exit(1);
		}
	}

	/*
	 * Miscellaneous chores.
	 */

	XLALPrintInfo("%s: using --psd-average-points %zu\n", argv[0], options->psd_length);
	if(options->max_series_length)
		XLALPrintInfo("%s: available RAM limits analysis to %d samples (%g s)\n", argv[0], options->max_series_length, options->max_series_length / (double) options->resample_rate);

	return options;
}


/*
 * ============================================================================
 *
 *                                   Input
 *
 * ============================================================================
 */


static REAL8TimeSeries *get_time_series(const char *cachefilename, const char *chname, LIGOTimeGPS start, LIGOTimeGPS end, size_t lengthlimit, int getcaltimeseries)
{
	static const char func[] = "get_time_series";
	double duration = XLALGPSDiff(&end, &start);
	FrCache *cache;
	FrStream *stream;
	REAL8TimeSeries *series;

	if(duration < 0)
		XLAL_ERROR_NULL(func, XLAL_EINVAL);

	/*
	 * Open frame stream.
	 */

	cache = XLALFrImportCache(cachefilename);
	if(!cache)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	stream = XLALFrCacheOpen(cache);
	XLALFrDestroyCache(cache);
	if(!stream)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/*
	 * Turn on checking for missing data.
	 */

	stream->mode = LAL_FR_VERBOSE_MODE;

	/*
	 * Get the data.
	 */

	if(getcaltimeseries) {
		series = XLALFrReadREAL8TimeSeries(stream, chname, &start, duration, lengthlimit);
		if(!series)
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		/* FIXME:  ARGH!!!  frame files cannot be trusted to
		 * provide units for their contents!  */
		series->sampleUnits = lalStrainUnit;
	} else {
		/* read data as single-precision, and convert to double */
		REAL4TimeSeries *tmp = XLALFrReadREAL4TimeSeries(stream, chname, &start, duration, lengthlimit);
		unsigned i;
		if(!tmp)
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		series = XLALCreateREAL8TimeSeries(tmp->name, &tmp->epoch, tmp->f0, tmp->deltaT, &tmp->sampleUnits, tmp->data->length);
		if(!series) {
			XLALDestroyREAL4TimeSeries(tmp);
			XLAL_ERROR_NULL(func, XLAL_EFUNC);
		}
		for(i = 0; i < tmp->data->length; i++)
			series->data->data[i] = tmp->data->data[i];
		XLALDestroyREAL4TimeSeries(tmp);
		/* FIXME:  ARGH!!!  frame files cannot be trusted to
		 * provide units for their contents!  */
		series->sampleUnits = lalADCCountUnit;
	}

	/*
	 * Check for missing data.
	 */

	if(stream->state & LAL_FR_GAP) {
		XLALPrintError("get_time_series(): error: gap in data detected between GPS times %d.%09u s and %d.%09u s\n", start.gpsSeconds, start.gpsNanoSeconds, end.gpsSeconds, end.gpsNanoSeconds);
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(func, XLAL_EDATA);
	}

	/*
	 * Close stream.
	 */

	XLALFrClose(stream);

	/*
	 * Check for other failures.
	 */

	if(!series)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/*
	 * Verbosity.
	 */

	XLALPrintInfo("get_time_series(): read %u samples (%.9lf s) at GPS time %u.%09u s\n", series->data->length, series->data->length * series->deltaT, start.gpsSeconds, start.gpsNanoSeconds);

	return series;
}


/*
 * ============================================================================
 *
 *                    Fill a time series with white noise
 *
 * ============================================================================
 */


static void gaussian_noise(REAL8TimeSeries *series, REAL8 rms, gsl_rng *rng)
{
	unsigned i;

	for(i = 0; i < series->data->length; i++)
		series->data->data[i] = gsl_ran_gaussian(rng, rms);
}


/*
 * ============================================================================
 *
 *                             Response function
 *
 * ============================================================================
 */


/*
 * Note, this function can only read new-style calibration data (S5 and
 * later).  Calibration data in earlier formats will need to be converted.
 */


static COMPLEX8FrequencySeries *generate_response(const char *cachefile, const char *channel_name, REAL8 deltaT, LIGOTimeGPS epoch, size_t length)
{
	static const char func[] = "generate_response";
	/* length of time interval spanned by calibration */
	const double duration = length * deltaT;
	/* frequency resolution of response function */
	const double deltaf = 1.0 / (length * deltaT);
	/* number of frequency bins in response function */
	const int n = length / 2 + 1;
	COMPLEX8FrequencySeries *response;

	XLALPrintInfo("generating response function spanning %g s at GPS time %u.%09u s with %g Hz resolution\n", duration, epoch.gpsSeconds, epoch.gpsNanoSeconds, deltaf);

	if(!cachefile) {
		/* generate fake unity response if working with calibrated
		 * data or if there is no calibration information available
		 * */
		const LALUnit strainPerCount = {0, {0, 0, 0, 0, 0, 1, -1}, {0, 0, 0, 0, 0, 0, 0}};
		const COMPLEX8 one = {1.0, 0.0};
		size_t i;

		response = XLALCreateCOMPLEX8FrequencySeries(channel_name, &epoch, 0.0, deltaf, &strainPerCount, n);
		if(!response)
			XLAL_ERROR_NULL(func, XLAL_ENOMEM);

		XLALPrintInfo("generate_response(): generating unit response function\n");

		for(i = 0; i < response->data->length; i++)
			response->data->data[i] = one;

		return response;
	} else {
		FrCache *cache = XLALFrImportCache(cachefile);
		if(cache) {
			LALCalData *caldata = XLALFrCacheGetCalData(&epoch, channel_name, cache);
			XLALFrDestroyCache(cache);
			if(caldata) {
				response = XLALCreateCOMPLEX8Response(&epoch, duration, deltaf, n, caldata);
				XLALDestroyCalData(caldata);
				if(response)
					return response;
			}
		}
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
}


/*
 * ============================================================================
 *
 *                              Burst injections
 *
 * ============================================================================
 */


static int add_burst_injections(REAL8TimeSeries *h, const char *filename)
{
	static const char func[] = "add_burst_injections";
	/* hard-coded speed hack.  only injections whose "times" are within
	 * this many seconds of the time series will be loaded */
	const double longest_injection = 600.0;
	LIGOTimeGPS start = h->epoch;
	LIGOTimeGPS end = h->epoch;
	SimBurst *sim_burst;

	XLALGPSAdd(&start, -longest_injection);
	XLALGPSAdd(&end, h->data->length * h->deltaT + longest_injection);

	XLALPrintInfo("%s(): reading sim_burst table from %s\n", func, filename);

	sim_burst = XLALSimBurstTableFromLIGOLw(filename, &start, &end);
	if(!sim_burst)
		XLAL_ERROR(func, XLAL_EFUNC);

	XLALPrintInfo("%s(): computing injections and adding to input time series\n", func);
	if(XLALBurstInjectSignals(h, sim_burst, NULL)) {
		while(sim_burst) {
			SimBurst *next = sim_burst->next;
			XLALDestroySimBurst(sim_burst);
			sim_burst = next;
		}
		XLAL_ERROR(func, XLAL_EFUNC);
	}

	XLALPrintInfo("%s(): done\n", func);

	while(sim_burst) {
		SimBurst *next = sim_burst->next;
		XLALDestroySimBurst(sim_burst);
		sim_burst = next;
	}

	return 0;
}


/*
 * ============================================================================
 *
 *                              Inspiral injections
 *
 * ============================================================================
 */


/*
 * LAL can only generate single-precision injections so we have to do the
 * injections into a temporary array of zeros, and then add it into the
 * double-precision time series
 */


static REAL8TimeSeries *add_inspiral_injections(LALStatus *stat, char *filename, REAL8TimeSeries *series, COMPLEX8FrequencySeries *response)
{
	static const char func[] = "add_inspiral_injections";
	INT4 start_time = series->epoch.gpsSeconds;
	INT4 stop_time = start_time + series->data->length * series->deltaT;
	SimInspiralTable *injections = NULL;
	REAL4TimeSeries *zeroes = XLALCreateREAL4TimeSeries(series->name, &series->epoch, series->f0, series->deltaT, &series->sampleUnits, series->data->length);
	INT4 numInjections = 0;
	unsigned i;

	if(!zeroes)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	memset(zeroes->data->data, 0, zeroes->data->length * sizeof(*zeroes->data->data));

	if(!response) {
		XLALDestroyREAL4TimeSeries(zeroes);
		XLALPrintError("%s(): error: must supply calibration information for injections\n", func);
		exit(1);
	}

	XLALPrintInfo("%s(): reading sim_inspiral table ...\n", func);

	numInjections = SimInspiralTableFromLIGOLw(&injections, filename, start_time, stop_time);

	if(numInjections < 0) {
		XLALDestroyREAL4TimeSeries(zeroes);
		XLALPrintError("%s(): error: cannot read injection file\n", func);
		exit(1);
	}

	XLALPrintInfo("%s(): computing injections ...\n", func);

	if(numInjections > 0)
		LAL_CALL(LALFindChirpInjectSignals(stat, zeroes, injections, response), stat);

	XLALPrintInfo("%s(): done\n", func);

	while(injections) {
		SimInspiralTable *thisEvent = injections;
		injections = injections->next;
		XLALFree(thisEvent);
	}

	for(i = 0; i < series->data->length; i++)
		series->data->data[i] += zeroes->data->data[i];
	XLALDestroyREAL4TimeSeries(zeroes);

	return series;
}


/*
 * ============================================================================
 *
 *                               MDC injections
 *
 * ============================================================================
 */


static REAL8TimeSeries *add_mdc_injections(const char *mdccachefile, const char *channel_name, REAL8TimeSeries *series)
{
	static const char func[] = "add_mdc_injections";
	LIGOTimeGPS stop;
	REAL8TimeSeries *mdc;
	unsigned i;

	XLALPrintInfo("%s(): mixing data from MDC frames\n", func);

	stop = series->epoch;
	XLALGPSAdd(&stop, series->data->length * series->deltaT);

	mdc = get_time_series(mdccachefile, channel_name, series->epoch, stop, 0, TRUE);
	if(!mdc)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	for(i = 0; i < series->data->length; i++)
		series->data->data[i] += mdc->data->data[i];

	XLALDestroyREAL8TimeSeries(mdc);

	return series;
}


/*
 * ============================================================================
 *
 *                               Analysis Loop
 *
 * ============================================================================
 */


/*
 * Analyze a time series in intervals corresponding to the length of time
 * for which the instrument's noise is stationary.
 */


static SnglBurstTable **analyze_series(SnglBurstTable **addpoint, REAL8TimeSeries *series, int psd_length, int psd_shift, struct options *options)
{
	static const char func[] = "analyze_series";
	unsigned i;

	if((unsigned) psd_length > series->data->length) {
		XLALPrintWarning("%s(): warning: PSD average length exceeds available data, skipping series\n", func);
		return addpoint;
	}

	for(i = 0; i + psd_length < series->data->length + psd_shift; i += psd_shift) {
		int start = min(i, series->data->length - psd_length);
		REAL8TimeSeries *interval = XLALCutREAL8TimeSeries(series, start, psd_length);

		if(!interval)
			XLAL_ERROR_NULL(func, XLAL_EFUNC);

		XLALPrintInfo("%s(): analyzing samples %zu -- %zu (%.9lf s -- %.9lf s)\n", func, start, start + interval->data->length, start * interval->deltaT, (start + interval->data->length) * interval->deltaT);

		*addpoint = XLALEPSearch(
			options->diagnostics,
			interval,
			options->window,
			options->flow,
			options->bandwidth,
			options->confidence_threshold,
			options->fractional_stride,
			options->maxTileBandwidth,
			options->maxTileDuration
		);
		while(*addpoint)
			addpoint = &(*addpoint)->next;

		XLALDestroyREAL8TimeSeries(interval);

		if(xlalErrno)
			XLAL_ERROR_NULL(func, XLAL_EFUNC);

	}

	return addpoint;
}


/*
 * ============================================================================
 *
 *                                   Output
 *
 * ============================================================================
 */


static void output_results(LALStatus *stat, char *file, const char *ifo, MetadataTable *_process_table, MetadataTable *_process_params_table, MetadataTable *_search_summary_table, SnglBurstTable *burstEvent)
{
	LIGOLwXMLStream xml;
	LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
	MetadataTable myTable;

	memset(&xml, 0, sizeof(LIGOLwXMLStream));
	LAL_CALL(LALOpenLIGOLwXMLFile(stat, &xml, file), stat);

	/*
	 * process table
	 */

	LAL_CALL(LALGPSTimeNow(stat, &(_process_table->processTable->end_time), &accuracy), stat);
	LAL_CALL(LALBeginLIGOLwXMLTable(stat, &xml, process_table), stat);
	LAL_CALL(LALWriteLIGOLwXMLTable(stat, &xml, *_process_table, process_table), stat);
	LAL_CALL(LALEndLIGOLwXMLTable(stat, &xml), stat);

	/*
	 * process params table
	 */

	LAL_CALL(LALBeginLIGOLwXMLTable(stat, &xml, process_params_table), stat);
	LAL_CALL(LALWriteLIGOLwXMLTable(stat, &xml, *_process_params_table, process_params_table), stat);
	LAL_CALL(LALEndLIGOLwXMLTable(stat, &xml), stat);

	/*
	 * search summary table
	 */

	snprintf(_search_summary_table->searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
	_search_summary_table->searchSummaryTable->nevents = XLALSnglBurstTableLength(burstEvent);
	LAL_CALL(LALBeginLIGOLwXMLTable(stat, &xml, search_summary_table), stat);
	LAL_CALL(LALWriteLIGOLwXMLTable(stat, &xml, *_search_summary_table, search_summary_table), stat);
	LAL_CALL(LALEndLIGOLwXMLTable(stat, &xml), stat);

	/*
	 * burst table
	 */

	LAL_CALL(LALBeginLIGOLwXMLTable(stat, &xml, sngl_burst_table), stat);
	myTable.snglBurstTable = burstEvent;
	LAL_CALL(LALWriteLIGOLwXMLTable(stat, &xml, myTable, sngl_burst_table), stat);
	LAL_CALL(LALEndLIGOLwXMLTable(stat, &xml), stat);

	LAL_CALL(LALCloseLIGOLwXMLFile(stat, &xml), stat);
}


/*
 * ============================================================================
 *
 *                                Progress Bar
 *
 * ============================================================================
 */


static void print_progress_bar(const char *func, const LIGOTimeGPS *start, const LIGOTimeGPS *end, const LIGOTimeGPS *t)
{
	static const char bar[] = "+++++++++++++++++++++++++++++++++++++++++++++++++)";
	static const char spc[] = "-------------------------------------------------)";
	double fraction = XLALGPSDiff(t, start) / XLALGPSDiff(end, start);
	int l = sizeof(bar)/sizeof(*bar) - 1;
	int offset = fraction * l;

	XLALPrintInfo("%s: [%s%s %.1f%% complete\n", func, bar + l - offset, spc + offset, 100.0 * fraction);
}


/*
 * ============================================================================
 *
 *                                Entry point
 *
 * ============================================================================
 */


int main(int argc, char *argv[])
{
	LALStatus stat;
	LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
	struct options *options;
	LIGOTimeGPS epoch;
	LIGOTimeGPS boundepoch;
	int overlap;
	REAL8TimeSeries *series = NULL;
	SnglBurstTable *burstEvent = NULL;
	SnglBurstTable **EventAddPoint = &burstEvent;
	/* the ugly underscores are because some badger put global symbols
	 * in LAL by exactly these names. it's a mad house, a maad house */
	MetadataTable _process_table;
	MetadataTable _process_params_table;
	MetadataTable _search_summary_table;
	gsl_rng *rng = NULL;

	/*
	 * Command line
	 */

	memset(&stat, 0, sizeof(stat));
	lal_errhandler = LAL_ERR_EXIT;
	if(parse_command_line_debug(argc, argv) < 0)
		exit(1);

	/*
	 * Create the process and process params tables.
	 */

	_process_table.processTable = XLALCalloc(1, sizeof(ProcessTable));
	LAL_CALL(LALGPSTimeNow(&stat, &(_process_table.processTable->start_time), &accuracy), &stat);
	LAL_CALL(populate_process_table(&stat, _process_table.processTable, PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE), &stat);
	_process_params_table.processParamsTable = NULL;

	/*
	 * Parse arguments and fill _process_params_table table.
	 */

	options = parse_command_line(argc, argv, &_process_params_table);
	snprintf(_process_table.processTable->ifos, LIGOMETA_IFOS_MAX, "%s", options->ifo);
	snprintf(_process_table.processTable->comment, LIGOMETA_COMMENT_MAX, "%s", options->comment);

	/*
	 * Create the search summary table.  The number of nodes for a
	 * standalone job is always 1
	 */

	_search_summary_table.searchSummaryTable = XLALCalloc(1, sizeof(SearchSummaryTable));
	snprintf(_search_summary_table.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, "%s", options->comment);
	_search_summary_table.searchSummaryTable->nnodes = 1;
	_search_summary_table.searchSummaryTable->in_start_time = options->gps_start;
	_search_summary_table.searchSummaryTable->in_end_time = options->gps_end;

	/*
	 * determine the input time series post-conditioning overlap, and
	 * set the outer loop's upper bound
	 */

	overlap = options->window->data->length - options->window_shift;
	XLALPrintInfo("%s: time series overlap is %zu samples (%.9lf s)\n", argv[0], overlap, overlap / (double) options->resample_rate);

	boundepoch = options->gps_end;
	XLALGPSAdd(&boundepoch, -(2 * options->filter_corruption + (int) overlap) / (double) options->resample_rate);
	XLALPrintInfo("%s: time series epochs must be less than %u.%09u s\n", argv[0], boundepoch.gpsSeconds, boundepoch.gpsNanoSeconds);


	/*
	 * ====================================================================
	 *
	 *                         Outer Loop
	 *
	 * ====================================================================
	 */

	/*
	 * Split the total length of time to be analyzed into time series
	 * small enough to fit in RAM.
	 */

	for(epoch = options->gps_start; XLALGPSCmp(&epoch, &boundepoch) < 0;) {
		/*
		 * Progress bar.
		 */

		print_progress_bar(argv[0], &options->gps_start, &boundepoch, &epoch);

		/*
		 * Get the data.
		 */

		if(options->cache_filename) {
			/*
			 * Read from frame files
			 */

			series = get_time_series(options->cache_filename, options->channel_name, epoch, options->gps_end, options->max_series_length, options->calibrated);
			if(!series) {
				XLALPrintError("%s: error: failure reading input data\n", argv[0]);
				exit(1);
			}
		} else if(options->noise_rms > 0.0) {
			/*
			 * Synthesize Gaussian white noise.
			 */

			unsigned length = XLALGPSDiff(&options->gps_end, &epoch) * options->resample_rate;
			if(options->max_series_length)
				length = min(options->max_series_length, length);
			series = XLALCreateREAL8TimeSeries(options->channel_name, &epoch, 0.0, (REAL8) 1.0 / options->resample_rate, &lalADCCountUnit, length);
			if(!series) {
				XLALPrintError("%s: error: failure allocating data for Gaussian noise\n", argv[0]);
				exit(1);
			}
			if(!rng) {
				rng = gsl_rng_alloc(gsl_rng_ranlxd1);
				if(options->seed)
					gsl_rng_set(rng, options->seed);
				else {
					/* use time in milliseconds */
					struct timeval t;
					unsigned long seed;
					if(gettimeofday(&t, NULL)) {
						/* failure */
						XLALPrintError("%s: error: cannot get time of day to seed random number generator\n", argv[0]);
						exit(1);
					}
					seed = 1000 * (t.tv_sec + t.tv_usec * 1e-6);
					XLALPrintInfo("%s: using random number seed %lu\n", seed);
					gsl_rng_set(rng, seed);
				}
			}
			gaussian_noise(series, options->noise_rms, rng);
		} else {
			/*
			 * Should never get here.
			 */
			XLALPrintError("%s: error: oops, don't know how to get data\n", argv[0]);
			exit(1);
		}

		/*
		 * Add burst/inspiral injections into the time series if
		 * requested.
		 */

		if(options->sim_burst_filename || options->sim_inspiral_filename) {
			COMPLEX8FrequencySeries *response;

			/* Create the response function (generates unity
			 * response if cache file is NULL). */
			response = generate_response(options->cal_cache_filename, options->channel_name, series->deltaT, series->epoch, series->data->length);
			if(!response)
				exit(1);

			/* perform injections */
			if(options->sim_burst_filename)
				if(add_burst_injections(series, options->sim_burst_filename))
					exit(1);
			if(options->sim_inspiral_filename)
				add_inspiral_injections(&stat, options->sim_inspiral_filename, series, response);

			/* clean up */
			XLALDestroyCOMPLEX8FrequencySeries(response);
		}

		/*
		 * Add MDC injections into the time series if requested.
		 */

		if(options->mdc_cache_filename)
			add_mdc_injections(options->mdc_cache_filename, options->mdc_channel_name, series);

		/*
		 * Condition the time series data.
		 */

		if(XLALEPConditionData(series, options->high_pass, (REAL8) 1.0 / options->resample_rate, options->filter_corruption)) {
			XLALPrintError("%s: XLALEPConditionData() failed.\n", argv[0]);
			exit(1);
		}

		/*
		 * Store the start and end times of the data that actually
		 * gets analyzed.
		 */

		if(!_search_summary_table.searchSummaryTable->out_start_time.gpsSeconds) {
			_search_summary_table.searchSummaryTable->out_start_time = series->epoch;
			XLALGPSAdd(&_search_summary_table.searchSummaryTable->out_start_time, series->deltaT * options->window_pad);
		}
		_search_summary_table.searchSummaryTable->out_end_time = series->epoch;
		XLALGPSAdd(&_search_summary_table.searchSummaryTable->out_end_time, series->deltaT * (series->data->length - options->window_pad));

		/*
		 * Analyze the data
		 */

		EventAddPoint = analyze_series(EventAddPoint, series, options->psd_length, options->psd_shift, options);
		if(!EventAddPoint)
			exit(1);

		/*
		 * Reset for next run
		 *
		 * Advancing the epoch by the post-conditioning series length
		 * provides exactly the overlap needed to account for
		 * conditioning corruption.  The post-conditioning overlap
		 * computed earlier provides the additional overlap that is
		 * needed.
		 */

		XLALGPSAdd(&epoch, (series->data->length - overlap) * series->deltaT);
		XLALDestroyREAL8TimeSeries(series);
	}

	/*
	 * Sort the events, and assign IDs.
	 */

	XLALSortSnglBurst(&burstEvent, XLALCompareSnglBurstByStartTimeAndLowFreq);
	XLALSnglBurstAssignIDs(burstEvent);

	/*
	 * Check event rate limit.
	 */

	if((options->max_event_rate > 0) && (XLALSnglBurstTableLength(burstEvent) > XLALGPSDiff(&_search_summary_table.searchSummaryTable->out_end_time, &_search_summary_table.searchSummaryTable->out_start_time) * options->max_event_rate)) {
		XLALPrintError("%s: event rate limit exceeded!", argv[0]);
		exit(1);
	}

	/*
	 * Output the results.
	 */

	output_results(&stat, options->output_filename, options->ifo, &_process_table, &_process_params_table, &_search_summary_table, burstEvent);

	/*
	 * Final cleanup.
	 */

	if(rng)
		gsl_rng_free(rng);
	XLALFree(_process_table.processTable);
	XLALFree(_search_summary_table.searchSummaryTable);

	while(_process_params_table.processParamsTable) {
		ProcessParamsTable *table = _process_params_table.processParamsTable;
		_process_params_table.processParamsTable = table->next;
		XLALFree(table);
	}

	while(burstEvent) {
		SnglBurstTable *event = burstEvent;
		burstEvent = burstEvent->next;
		XLALFree(event);
	}

	if(options->diagnostics)
		LAL_CALL(LALCloseLIGOLwXMLFile(&stat, options->diagnostics->LIGOLwXMLStream), &stat);
	options_free(options);

	LALCheckMemoryLeaks();
	exit(0);
}
