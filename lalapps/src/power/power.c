/*
*  Copyright (C) 2007 Denny Mackin, Duncan Brown, Ik Siong Heng, Jolien
*  Creighton, Kipp Cannon, Mark Williamsen, Patrick Brady, Robert Adam
*  Mercer, Saikat Ray-Majumder, Stephen Fairhurst
*
*  This program is free software; you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by the
*  Free Software Foundation; either version 2 of the License, or (at your
*  option) any later version.
*
*  This program is distributed in the hope that it will be useful, but
*  WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*  General Public License for more details.
*
*  You should have received a copy of the GNU General Public License along
*  with with program; see the file COPYING. If not, write to the Free
*  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
*  02111-1307  USA
*/

#include <getopt.h>
#include <lalapps.h>
#include <math.h>
#include <processtable.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
 * Round length down so that an integer number of intervals of length
 * block_length, each shifted by block_shift with respect to its
 * neighbours, fits into the result.  This is used to ensure that an
 * integer number of analysis windows fits into the PSD length, and also
 * that an integer number of PSD lengths fits into the RAM limit length.
 */


static size_t block_commensurate(size_t length, size_t block_length, size_t block_shift)
{
	size_t blocks = (length - block_length) / block_shift;

	return length < block_length ? 0 : blocks * block_shift + block_length;
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
		return FALSE;
	if(x - trunc(x) == 0)
		/* is an integer */
		return is_power_of_2(x);
	if(x < 1)
		/* is less than 1 */
		return double_is_power_of_2(1 / x);
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

	int bandwidth;
	unsigned window_length;
	/* number of samples to use for PSD    */
	size_t psd_length;
	/* set non-zero to generate noise      */
	int seed;

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
	/* high pass filter cut-off for double --> single quantization */
	double cal_high_pass;

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
	 * WTF
	 */

	/* file with list of sim injections */
	char *simInjectionFile;
	/* name of the sim waveform cache file */
	char *sim_cache_filename;
	/* Distance at which the sim waveforms have been generated */
	double sim_distance;

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
	LIGOLwXMLStream *diagnostics;
};


static struct options *options_new(void)
{
	static char default_comment[] = "";
	struct options *options = malloc(sizeof(*options));

	if(!options)
		return NULL;

	*options = (struct options) {
		.bandwidth = 0,	/* impossible */
		.cal_cache_filename = NULL,	/* default */
		.channel_name = NULL,	/* impossible */
		.comment = default_comment,	/* default */
		.filter_corruption = -1,	/* impossible */
		.mdc_cache_filename = NULL,	/* default == disable */
		.mdc_channel_name = NULL,	/* default */
		.noise_rms = -1.0,	/* default == disable */
		.diagnostics = NULL,	/* default == disable */
		.psd_length = 0,	/* impossible */
		.resample_rate = 0,	/* impossible */
		.seed = 1,	/* default */
		.max_series_length = 0,	/* default */
		.calibrated = FALSE,	/* default */
		.high_pass = -1.0,	/* impossible */
		.cal_high_pass = -1.0,	/* impossible */
		.max_event_rate = 0,	/* default */
		.output_filename = NULL,	/* impossible */
		.window_length = 0,	/* impossible */
		.sim_distance = 10000.0,	/* default (10 Mpc) */
		.sim_cache_filename = NULL,	/* default */
		.sim_burst_filename = NULL,	/* default == disable */
		.simInjectionFile = NULL,	/* default == disable */
		.sim_inspiral_filename = NULL,	/* default == disable */
		.cache_filename = NULL,	/* default == disable */
	};

	memset(options->ifo, 0, sizeof(options->ifo));	/* empty */
	XLALINT8NSToGPS(&options->gps_start, 0);	/* impossible */
	XLALINT8NSToGPS(&options->gps_end, 0);	/* impossible */

	return options;
}


static void options_free(struct options *options)
{
	free(options);
}


/*
 * Friendly messages.
 */


static void print_usage(const char *program)
{
	fprintf(stderr,
"Usage:  %s <option> [...]\n" \
"The following options are recognized.  Options not surrounded in [] are\n" \
"required.\n" \
"	 --bandwidth <bandwidth>\n" \
"	[--burstinjection-file <file name>]\n" \
"	[--calibrated-data <high pass frequency>]\n" \
"	[--calibration-cache <cache file>]\n" \
"	 --channel-name <string>\n" \
"	 --confidence-threshold <confidence>\n" \
"	[--debug-level info|warn|error|off]\n" \
"	[--dump-diagnostics <xml filename>]\n" \
"	[--enable-over-whitening]\n" \
"	 --filter-corruption <samples>\n" \
"	 --frame-cache <cache file>\n", program);
fprintf(stderr,
"	[--gaussian-noise-rms <rms amplitude>]\n" \
"	 --gps-end-time <seconds>\n" \
"	 --gps-start-time <seconds>\n" \
"	[--help]\n" \
"	 --high-pass <high pass frequency>\n" \
"	[--inspiralinjection-file <file name>]\n" \
"	 --low-freq-cutoff <Hz>\n" \
"	[--max-event-rate <Hz>]\n" \
"	 --max-tileband <Hz>\n" \
"	 --max-tileduration <samples>\n" \
"	[--mdc-cache <cache file>]\n" \
"	[--mdc-channel <channel name>]\n" \
"	[--output <filename>]\n" \
"	 --psd-average-method <method>\n" \
"	 --psd-average-points <samples>\n");
fprintf(stderr,
"	[--ram-limit <MebiBytes>]\n" \
"	 --resample-rate <Hz>\n" \
"	[--seed <seed>]\n" \
"	[--sim-cache <sim cache file>]\n" \
"	[--sim-seconds <sim seconds>]\n" \
"	[--sim-distance <sim distance(Kpc)>]\n" \
"	[--siminjection-file <file name>]\n" \
"	 --tile-stride-fraction <fraction>\n" \
"	[--user-tag <comment>]\n" \
"	 --window-length <samples>\n" \
"	 --window-shift <samples>\n");
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


static int all_required_arguments_present(char *prog, struct option *long_options, const struct options *options, const EPSearchParams *params)
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
			arg_is_missing = params->tf_flow < 0.0;
			break;

		case 'W':
			arg_is_missing = !options->window_length;
			break;

		case 'Y':
			arg_is_missing = params->method == (unsigned) -1;
			break;

		case 'Z':
			arg_is_missing = !options->psd_length;
			break;

		case 'd':
			arg_is_missing = !params->windowShift;
			break;

		case 'e':
			arg_is_missing = !options->resample_rate;
			break;

		case 'f':
			arg_is_missing = !params->fractional_stride;
			break;

		case 'g':
			arg_is_missing = params->confidence_threshold == XLAL_REAL8_FAIL_NAN;
			break;

		case 'j':
			arg_is_missing = options->filter_corruption < 0;
			break;

		case 'l':
			arg_is_missing = !params->maxTileBandwidth;
			break;

		case 'm':
			arg_is_missing = !params->maxTileDuration;
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


static struct options *parse_command_line(int argc, char *argv[], EPSearchParams *params, MetadataTable *_process_params_table)
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
		{"calibrated-data", required_argument, NULL, 'J'},
		{"calibration-cache", required_argument, NULL, 'B'},
		{"channel-name", required_argument, NULL, 'C'},
		{"confidence-threshold", required_argument, NULL, 'g'},
		{"debug-level", required_argument, NULL, 'D'},
		{"dump-diagnostics", required_argument, NULL, 'X'},
		{"enable-over-whitening", no_argument, &params->useOverWhitening, TRUE},
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
		{"max-tileband", required_argument, NULL, 'l'},
		{"max-tileduration", required_argument, NULL, 'm'},
		{"mdc-cache", required_argument, NULL, 'R'},
		{"mdc-channel", required_argument, NULL, 'S'},
		{"output", required_argument, NULL, 'b'},
		{"psd-average-method", required_argument, NULL, 'Y'},
		{"psd-average-points", required_argument, NULL, 'Z'},
		{"ram-limit", required_argument, NULL, 'a'},
		{"resample-rate", required_argument, NULL, 'e'},
		{"seed", required_argument, NULL, 'c'},
		{"sim-cache", required_argument, NULL, 'q'},
		{"sim-distance", required_argument, NULL, 'u'},
		{"siminjection-file", required_argument, NULL, 't'},
		{"tile-stride-fraction", required_argument, NULL, 'f'},
		{"user-tag", required_argument, NULL, 'h'},
		{"window-length", required_argument, NULL, 'W'},
		{"window-shift", required_argument, NULL, 'd'},
		{NULL, 0, NULL, 0}
	};

	/*
	 * Allocate and initialize options structure
	 */

	options = options_new();
	if(!options)
		return NULL;

	/*
	 * Set parameter defaults
	 */

	params->diagnostics = NULL;	/* default == disable */
	params->confidence_threshold = XLAL_REAL8_FAIL_NAN;	/* impossible */
	params->method = -1;	/* impossible */
	params->tf_freqBins = 0;	/* impossible */
	params->tf_deltaF = 0;	/* impossible */
	params->tf_flow = -1.0;	/* impossible */
	params->maxTileBandwidth = 0;	/* impossible */
	params->maxTileDuration = 0;	/* impossible */
	params->fractional_stride = 0;	/* impossible */
	params->windowShift = 0;	/* impossible */
	params->useOverWhitening = FALSE;	/* default */

	/*
	 * Parse command line.
	 */

	opterr = 1;		/* enable error messages */
	optind = 0;		/* start scanning from argv[0] */
	do switch (c = getopt_long(argc, argv, "", long_options, &option_index)) {
	case 'A':
		options->bandwidth = atoi(optarg);
		if(options->bandwidth <= 0 || !is_power_of_2(options->bandwidth)) {
			sprintf(msg, "must be greater than 0 and a power of 2 (%i specified)", options->bandwidth);
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
		options->cal_high_pass = atof(optarg);
		if(options->cal_high_pass < 0) {
			sprintf(msg, "must not be negative (%f Hz specified)", options->cal_high_pass);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("double");
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
		params->tf_flow = atof(optarg);
		if(params->tf_flow < 0) {
			sprintf(msg, "must not be negative (%f Hz specified)", params->tf_flow);
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
		options->window_length = atoi(optarg);
		if((options->window_length < 2) || !is_power_of_2(options->window_length)) {
			sprintf(msg, "must be greater than or equal to 2 and a power of 2 (%i specified)", options->window_length);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

	case 'X':
		{
		LALStatus stat;
		memset(&stat, 0, sizeof(stat));
		options->diagnostics = calloc(1, sizeof(*options->diagnostics));
		LALOpenLIGOLwXMLFile(&stat, options->diagnostics, optarg);
		}
		params->diagnostics = malloc(sizeof(*params->diagnostics));
		params->diagnostics->LIGOLwXMLStream = options->diagnostics;
		params->diagnostics->XLALWriteLIGOLwXMLArrayREAL4FrequencySeries = XLALWriteLIGOLwXMLArrayREAL4FrequencySeries;
		params->diagnostics->XLALWriteLIGOLwXMLArrayREAL4TimeSeries = XLALWriteLIGOLwXMLArrayREAL4TimeSeries;
		params->diagnostics->XLALWriteLIGOLwXMLArrayCOMPLEX8FrequencySeries = XLALWriteLIGOLwXMLArrayCOMPLEX8FrequencySeries;
		break;

	case 'Y':
		if(!strcmp(optarg, "useMean"))
			params->method = useMean;
		else if(!strcmp(optarg, "useMedian"))
			params->method = useMedian;
		else {
			sprintf(msg, "must be \"useMean\", or \"useMedian\"");
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("string");
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
		options->seed = atoi(optarg);
		if(options->seed <= 0) {
			sprintf(msg, "must be greater than 0 (%i specified)", options->seed);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

	case 'd':
		params->windowShift = atoi(optarg);
		ADD_PROCESS_PARAM("int");
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
		params->fractional_stride = atof(optarg);
		if(params->fractional_stride > 1 || !double_is_power_of_2(params->fractional_stride)) {
			sprintf(msg, "must be less than or equal to 1 and a power of 2 (%g specified)", params->fractional_stride);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

	case 'g':
		params->confidence_threshold = atof(optarg);
		if(params->confidence_threshold < 0) {
			sprintf(msg, "must not be negative (%g specified)", params->confidence_threshold);
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
		params->maxTileBandwidth = atof(optarg);
		if((params->maxTileBandwidth <= 0) || !double_is_power_of_2(params->maxTileBandwidth)) {
			sprintf(msg, "must be greater than 0 and a power of 2 (%f specified)", params->maxTileBandwidth);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

	case 'm':
		params->maxTileDuration = atof(optarg);
		if((params->maxTileDuration <= 0) || !double_is_power_of_2(params->maxTileDuration)) {
			sprintf(msg, "must be greater than 0 and a power of 2 (%f specified)", params->maxTileDuration);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		params->tf_deltaF = 1 / (2 * params->maxTileDuration);
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

	case 'q':
		options->sim_cache_filename = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 't':
		options->simInjectionFile = optarg;
		ADD_PROCESS_PARAM("string");
		break;

	case 'u':
		options->sim_distance = atof(optarg);
		if(options->sim_distance <= 0) {
			sprintf(msg, "must be greater than 0 (%f kpc specified), check the specification file for the generation of the sim waveforms to get the distance value", options->sim_distance);
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

	if(!all_required_arguments_present(argv[0], long_options, options, params))
		args_are_bad = TRUE;

	/*
	 * Make sure windowShift and window_length are OK.
	 */

	/* FIXME: check that params->maxTileDuration *
	 * params->fractional_stride yields an allowed delta t */

	if(options->window_length / 2 != block_commensurate(options->window_length / 2, params->maxTileDuration * options->resample_rate, params->maxTileDuration * params->fractional_stride * options->resample_rate)) {
		sprintf(msg, "an integer number of the largest tiles (duration = %g s) must fit into 1/2 of the window length", params->maxTileDuration);
		print_bad_argument(argv[0], "window-length", msg);
		args_are_bad = TRUE;
	}

	if(options->window_length / 2 < params->windowShift) {
		sprintf(msg, "must be >= 2 * --window-shift = %u (%u specified)", 2 * params->windowShift, options->window_length);
		print_bad_argument(argv[0], "window-length", msg);
		args_are_bad = TRUE;
	}

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
	 * Ensure psd_length is comensurate with the analysis window length
	 * and its shift.
	 */

	options->psd_length = block_commensurate(options->psd_length, options->window_length, params->windowShift);

	/*
	 * Ensure RAM limit is comensurate with the psd_length and its
	 * shift.
	 */

	if(options->max_series_length)
		options->max_series_length = block_commensurate(options->max_series_length, options->psd_length, options->psd_length - (options->window_length - params->windowShift));

	/*
	 * Generate time-domain window function.
	 */

	params->window = XLALCreateHannREAL4Window(options->window_length);
	if(!params->window) {
		XLALPrintError("%s: failure generating Hann window\n", argv[0]);
		exit(1);
	}

	/*
	 * Sanitize filter frequencies.
	 */

	if(options->cal_high_pass > params->tf_flow)
		XLALPrintWarning("%s: warning: calibrated data quantization high-pass frequency (%f Hz) greater than TF plane low frequency (%f Hz)\n", argv[0], options->cal_high_pass, params->tf_flow);

	if(options->high_pass > params->tf_flow - 10.0)
		XLALPrintWarning("%s: warning: data conditioning high-pass frequency (%f Hz) greater than 10 Hz below TF plane low frequency (%f Hz)\n", argv[0], options->high_pass, params->tf_flow);

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

	params->tf_freqBins = options->bandwidth / params->tf_deltaF;
	XLALPrintInfo("%s: using --psd-average-points %zu\n", argv[0], options->psd_length);
	if(options->max_series_length)
		XLALPrintInfo("%s: available RAM limits analysis to %d samples\n", argv[0], options->max_series_length);

	return options;
}


/*
 * ============================================================================
 *
 *                                   Input
 *
 * ============================================================================
 */


static REAL4TimeSeries *get_calibrated_data(FrStream *stream, const char *chname, LIGOTimeGPS *start, double duration, size_t lengthlimit, double quantize_high_pass)
{
	static const char func[] = "get_calibrated_data";
	REAL8TimeSeries *calibrated;
	REAL4TimeSeries *series;
	PassBandParamStruc highpassParam;
	unsigned int i;

	/*
	 * retrieve calibrated data as REAL8 time series
	 */

	calibrated = XLALFrReadREAL8TimeSeries(stream, chname, start, duration, lengthlimit);
	if(!calibrated)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/*
	 * high pass filter before casting REAL8 to REAL4
	 */

	highpassParam.nMax = 4;
	highpassParam.f2 = quantize_high_pass;
	highpassParam.f1 = -1.0;
	highpassParam.a2 = 0.9;
	highpassParam.a1 = -1.0;
	if(XLALButterworthREAL8TimeSeries(calibrated, &highpassParam) < 0) {
		XLALDestroyREAL8TimeSeries(calibrated);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/*
	 * copy data into a REAL4 time series
	 */

	series = XLALCreateREAL4TimeSeries(calibrated->name, &calibrated->epoch, calibrated->f0, calibrated->deltaT, &calibrated->sampleUnits, calibrated->data->length);
	if(!series) {
		XLALDestroyREAL8TimeSeries(calibrated);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < series->data->length; i++)
		series->data->data[i] = calibrated->data->data[i];

	XLALDestroyREAL8TimeSeries(calibrated);
	return series;
}


static REAL4TimeSeries *get_time_series(const char *cachefilename, const char *chname, LIGOTimeGPS start, LIGOTimeGPS end, size_t lengthlimit, int getcaltimeseries, double quantize_high_pass)
{
	const char func[] = "get_time_series";
	double duration = XLALGPSDiff(&end, &start);
	FrCache *cache;
	FrStream *stream;
	REAL4TimeSeries *series;

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

	if(getcaltimeseries)
		series = get_calibrated_data(stream, chname, &start, duration, lengthlimit, quantize_high_pass);
	else
		series = XLALFrReadREAL4TimeSeries(stream, chname, &start, duration, lengthlimit);

	/*
	 * Check for missing data.
	 */

	if(stream->state & LAL_FR_GAP) {
		XLALPrintError("get_time_series(): error: gap in data detected between GPS times %d.%09u s and %d.%09u s\n", start.gpsSeconds, start.gpsNanoSeconds, end.gpsSeconds, end.gpsNanoSeconds);
		XLALDestroyREAL4TimeSeries(series);
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


static void gaussian_noise(REAL4TimeSeries *series, REAL4 rms, RandomParams *rparams)
{
	unsigned i;

	XLALNormalDeviates(series->data, rparams);
	for(i = 0; i < series->data->length; i++)
		series->data->data[i] *= rms;
}


/*
 * ============================================================================
 *
 *                             Response function
 *
 * ============================================================================
 */


static COMPLEX8FrequencySeries *generate_response(LALStatus *stat, const char *calcachefile, char *ifo, const char *channel_name, REAL8 deltaT, LIGOTimeGPS epoch, size_t length)
{
	const char func[] = "generate_response";
	COMPLEX8FrequencySeries *response;
	size_t i;
	const LALUnit strainPerCount = { 0, {0, 0, 0, 0, 0, 1, -1}, {0, 0, 0, 0, 0, 0, 0} };
	COMPLEX8 one = { 1.0, 0.0 };
	CalibrationUpdateParams calfacts;
	FrCache *calcache = NULL;

	memset(&calfacts, 0, sizeof(calfacts));
	calfacts.ifo = ifo;

	response = XLALCreateCOMPLEX8FrequencySeries(channel_name, &epoch, 0.0, 1.0 / (length * deltaT), &strainPerCount, length / 2 + 1);
	if(!response) {
		XLAL_ERROR_NULL(func, XLAL_ENOMEM);
		return NULL;
	}

	XLALPrintInfo("generate_response(): working at GPS time %u.%09u s\n", response->epoch.gpsSeconds, response->epoch.gpsNanoSeconds);

	if(!calcachefile) {
		/* generate fake unity response if working with calibrated
		 * data or if there is no calibration information available
		 * */
		XLALPrintInfo("generate_response(): generating unit response function\n");
		for(i = 0; i < response->data->length; i++)
			response->data->data[i] = one;
	} else {
		calcache = XLALFrImportCache(calcachefile);
		LAL_CALL(LALExtractFrameResponse(stat, response, calcache, &calfacts), stat);
		XLALFrDestroyCache(calcache);
	}

	return response;
}


/*
 * ============================================================================
 *
 *                              Burst injections
 *
 * ============================================================================
 */


static REAL4TimeSeries *add_burst_injections(LALStatus *stat, char *filename, REAL4TimeSeries *series, COMPLEX8FrequencySeries *response)
{
	const INT4 start_time = series->epoch.gpsSeconds;
	const INT4 stop_time = start_time + series->data->length * series->deltaT;
	const INT4 calType = 0;
	SimBurstTable *injections = NULL;

	if(!response) {
		XLALPrintError("add_burst_injections(): must supply calibration information for injections\n");
		exit(1);
	}

	XLALPrintInfo("add_burst_injections(): reading in SimBurst Table\n");

	LAL_CALL(LALSimBurstTableFromLIGOLw(stat, &injections, filename, start_time, stop_time), stat);

	XLALPrintInfo("add_burst_injections(): injecting signals into time series\n");

	LAL_CALL(LALBurstInjectSignals(stat, series, injections, response, calType), stat);

	XLALPrintInfo("add_burst_injections(): finished making the injections\n");

	while(injections) {
		SimBurstTable *thisEvent = injections;
		injections = injections->next;
		XLALFree(thisEvent);
	}

	return series;
}


/*
 * ============================================================================
 *
 *                              Inspiral injections
 *
 * ============================================================================
 */


static REAL4TimeSeries *add_inspiral_injections(LALStatus *stat, char *filename, REAL4TimeSeries *series, COMPLEX8FrequencySeries *response)
{
	INT4 start_time = series->epoch.gpsSeconds;
	INT4 stop_time = start_time + series->data->length * series->deltaT;
	SimInspiralTable *injections = NULL;

	INT4 numInjections = 0;

	if(!response) {
		XLALPrintError("add_inspiral_injections(): must supply calibration information for injections\n");
		exit(1);
	}

	XLALPrintInfo("add_inspiral_injections(): reading in SimInspiral Table\n");

	numInjections = SimInspiralTableFromLIGOLw(&injections, filename, start_time, stop_time);

	if(numInjections < 0) {
		XLALPrintError("add_inspiral_injections():error:cannot read injection file\n");
		exit(1);
	}

	XLALPrintInfo("add_inspiral_injections(): injecting signals into time series\n");

	if(numInjections > 0)
		LAL_CALL(LALFindChirpInjectSignals(stat, series, injections, response), stat);

	XLALPrintInfo("add_inspiral_injections(): finished making the injections\n");

	while(injections) {
		SimInspiralTable *thisEvent = injections;
		injections = injections->next;
		XLALFree(thisEvent);
	}

	return series;
}


/*
 * ============================================================================
 *
 *                               MDC injections
 *
 * ============================================================================
 */


static REAL4TimeSeries *add_mdc_injections(const char *mdccachefile, const char *channel_name, REAL4TimeSeries *series, LIGOTimeGPS epoch, LIGOTimeGPS stopepoch, size_t lengthlimit)
{
	const char func[] = "add_mdc_injections";
	REAL4TimeSeries *mdc;
	size_t i;

	XLALPrintInfo("add_mdc_injections(): mixing data from MDC frames\n");

	/*
	 * note: quantization high pass at 40.0 Hz
	 */

	mdc = get_time_series(mdccachefile, channel_name, epoch, stopepoch, lengthlimit, TRUE, 40.0);
	if(!mdc)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/*
	 * add the mdc signal to the given time series
	 */

	for(i = 0; i < series->data->length; i++)
		series->data->data[i] += mdc->data->data[i];

	/*
	 * clean up
	 */

	XLALDestroyREAL4TimeSeries(mdc);

	return series;
}


/*
 * ============================================================================
 *
 *                               Sim injections
 *
 * ============================================================================
 */


#if 0
static void add_sim_injections(LALStatus *stat, REAL4TimeSeries *series, COMPLEX8FrequencySeries *response, size_t lengthlimit)
{
	REAL4TimeSeries *signal;
	DetectorResponse detector;
	LALDetector *tmpDetector = NULL;
	CoherentGW waveform;
	REAL4 *aData;
	LALTimeInterval epochCorrection;
	COMPLEX8FrequencySeries *transfer = NULL;
	COMPLEX8Vector *unity = NULL;

	char pluschan[30];
	char crosschan[30];
	LIGOTimeGPS start, end;
	UINT4 i, n;
	REAL8 simDuration;

	INT4 start_time = series->epoch.gpsSeconds;
	INT4 stop_time = start_time + series->data->length * series->deltaT;
	SimBurstTable *injections = NULL;
	SimBurstTable *simBurst = NULL;
	BurstParamStruc burstParam;

	if(!response) {
		XLALPrintError("add_sim_injections(): must supply calibration information for injections\n");
		exit(1);
	}

	/* allocate memory */
	memset(&detector, 0, sizeof(DetectorResponse));
	transfer = (COMPLEX8FrequencySeries *) XLALCalloc(1, sizeof(COMPLEX8FrequencySeries));
	if(!transfer) {
		XLALPrintError("add_sim_injections(): detector.transfer not allocated\n");
		exit(1);
	}

	memcpy(&(transfer->epoch), &(response->epoch), sizeof(LIGOTimeGPS));
	transfer->f0 = response->f0;
	transfer->deltaF = response->deltaF;

	tmpDetector = detector.site = (LALDetector *) LALMalloc(sizeof(LALDetector));

	/* set the detector site */
	switch(series->name[0]) {
	case 'H':
		*(detector.site) = lalCachedDetectors[LALDetectorIndexLHODIFF];
		LALWarning(stat, "computing waveform for Hanford.");
		break;
	case 'L':
		*(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
		LALWarning(stat, "computing waveform for Livingston.");
		break;
	case 'G':
		*(detector.site) = lalCachedDetectors[LALDetectorIndexGEO600DIFF];
		LALWarning(stat, "computing waveform for GEO600.");
		break;
	default:
		XLALFree(detector.site);
		detector.site = NULL;
		tmpDetector = NULL;
		LALWarning(stat, "Unknown detector site, computing plus mode " "waveform with no time delay");
		break;
	}

	/* set up units for the transfer function */
	{
	RAT4 negOne = { -1, 0 };
	LALUnit unit;
	LALUnitPair pair;
	pair.unitOne = &lalADCCountUnit;
	pair.unitTwo = &lalStrainUnit;
	LAL_CALL(LALUnitRaise(stat, &unit, pair.unitTwo, &negOne), stat);
	pair.unitTwo = &unit;
	LAL_CALL(LALUnitMultiply(stat, &(transfer->sampleUnits), &pair), stat);
	}

	/* invert the response function to get the transfer function */
	LAL_CALL(LALCCreateVector(stat, &(transfer->data), response->data->length), stat);

	LAL_CALL(LALCCreateVector(stat, &unity, response->data->length), stat);
	for(i = 0; i < response->data->length; ++i) {
		unity->data[i].re = 1.0;
		unity->data[i].im = 0.0;
	}

	LAL_CALL(LALCCVectorDivide(stat, transfer->data, unity, response->data), stat);

	LAL_CALL(LALCDestroyVector(stat, &unity), stat);

	/* Set up a time series to hold signal in ADC counts */
	signal = XLALCreateREAL4TimeSeries(series->name, &series->epoch, series->f0, series->deltaT, &series->sampleUnits, series->data->length);

	XLALPrintInfo("add_sim_injections(): reading in SimBurst Table\n");

	LAL_CALL(LALSimBurstTableFromLIGOLw(stat, &injections, options.simInjectionFile, start_time, stop_time), stat);

	simBurst = injections;
	while(simBurst) {
		REAL4TimeSeries *plusseries = NULL;
		REAL4TimeSeries *crossseries = NULL;

		/* set the burst params */
		burstParam.deltaT = series->deltaT;
		if(!(strcmp(simBurst->coordinates, "HORIZON"))) {
			burstParam.system = COORDINATESYSTEM_HORIZON;
		} else if(!(strcmp(simBurst->coordinates, "ZENITH"))) {
			/* set coordinate system for completeness */
			burstParam.system = COORDINATESYSTEM_EQUATORIAL;
			detector.site = NULL;
		} else if(!(strcmp(simBurst->coordinates, "GEOGRAPHIC"))) {
			burstParam.system = COORDINATESYSTEM_GEOGRAPHIC;
		} else if(!(strcmp(simBurst->coordinates, "EQUATORIAL"))) {
			burstParam.system = COORDINATESYSTEM_EQUATORIAL;
		} else if(!(strcmp(simBurst->coordinates, "ECLIPTIC"))) {
			burstParam.system = COORDINATESYSTEM_ECLIPTIC;
		} else if(!(strcmp(simBurst->coordinates, "GALACTIC"))) {
			burstParam.system = COORDINATESYSTEM_GALACTIC;
		} else
			burstParam.system = COORDINATESYSTEM_EQUATORIAL;

		/* Set the channel names */
		snprintf(pluschan, 30, "SIM_plus_%s_%d", simBurst->waveform, simBurst->zm_number);
		snprintf(crosschan, 30, "SIM_cross_%s_%d", simBurst->waveform, simBurst->zm_number);

		/*Set the start and end times of the sim waveforms */
		start.gpsSeconds = 0;
		start.gpsNanoSeconds = 0;
		simDuration = (simBurst->dtplus + simBurst->dtminus);
		n = (INT4) (simDuration / series->deltaT);

		end = start;
		XLALGPSAdd(&end, simDuration);

		/* Get the plus time series */
		plusseries = get_time_series(options.sim_cache_filename, pluschan, start, end, lengthlimit, FALSE, options.cal_high_pass);

		/* Get the cross time series */
		crossseries = get_time_series(options.sim_cache_filename, crosschan, start, end, lengthlimit, FALSE, options.cal_high_pass);

		/* read in the waveform in a CoherentGW struct */
		memset(&waveform, 0, sizeof(CoherentGW));
		/* this step sets the adata,fdata,phidata to 0 */
		LAL_CALL(LALGenerateBurst(stat, &waveform, simBurst, &burstParam), stat);

		/* correct the waveform epoch:
		 * Remember the peak time is always at the center of the frame
		 * Hence the epoch of the waveform is set at 1/2 of the duration
		 * before the geocent_peak_time, since geocent_peak_time should match
		 * with the peak time in the frame
		 */
		simDuration = simDuration / 2.0;
		LAL_CALL(LALFloatToInterval(stat, &epochCorrection, &simDuration), stat);
		LAL_CALL(LALDecrementGPS(stat, &(waveform.a->epoch), &(simBurst->geocent_peak_time), &epochCorrection), stat);

		aData = waveform.a->data->data;

		/* copy the plus and cross data properly scaled for
		 * distance.
		 *
		 * NOTE: options.sim_distance specify the distance at which
		 * the simulated waveforms are produced. Check that the
		 * specified distance is 10 times of the distance as in
		 * parameter file, BBHWaveGen.in. Since in the wave
		 * generation script 1 kpc = 3.086e20 m where as the right
		 * definition is 1 kpc = 3.086e19 m.
		 */
		for(i = 0; i < n; i++) {
			*(aData++) = plusseries->data->data[i] * options.sim_distance / simBurst->distance;
			*(aData++) = crossseries->data->data[i] * options.sim_distance / simBurst->distance;
		}

		/* must set the epoch of signal since it's used by coherent
		 * GW */
		signal->epoch = waveform.a->epoch;
		detector.transfer = NULL;

		/* convert this into an ADC signal */
		LAL_CALL(LALSimulateCoherentGW(stat, signal, &waveform, &detector), stat);
		XLALRespFilt(signal, transfer);

		/* inject the signal into the data channel */
		LAL_CALL(LALSSInjectTimeSeries(stat, series, signal), stat);

		/* Clean up */
		XLALDestroyREAL4TimeSeries(plusseries);
		XLALDestroyREAL4TimeSeries(crossseries);
		LAL_CALL(LALSDestroyVectorSequence(stat, &(waveform.a->data)), stat);
		LAL_CALL(LALSDestroyVector(stat, &(waveform.f->data)), stat);
		LAL_CALL(LALDDestroyVector(stat, &(waveform.phi->data)), stat);
		LALFree(waveform.a);
		waveform.a = NULL;
		LALFree(waveform.f);
		waveform.f = NULL;
		LALFree(waveform.phi);
		waveform.phi = NULL;

		/* reset the detector site information in case it changed */
		detector.site = tmpDetector;

		/* move on to next one */
		simBurst = simBurst->next;
	}

	XLALDestroyREAL4TimeSeries(signal);
	LAL_CALL(LALCDestroyVector(stat, &(transfer->data)), stat);

	if(detector.site)
		XLALFree(detector.site);
	XLALFree(transfer);

	while(injections) {
		SimBurstTable *thisEvent = injections;
		injections = injections->next;
		XLALFree(thisEvent);
	}
}
#endif


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


static SnglBurstTable **analyze_series(SnglBurstTable **addpoint, REAL4TimeSeries *series, size_t psdlength, size_t window_length, EPSearchParams *params)
{
	REAL4TimeSeries *interval;
	size_t i, start;
	size_t overlap = window_length - params->windowShift;

	if(psdlength > series->data->length) {
		psdlength = block_commensurate(series->data->length, window_length, params->windowShift);
		XLALPrintInfo("analyze_series(): warning: PSD average length exceeds available data --- reducing PSD average length to %zu samples\n", psdlength);
		if(!psdlength) {
			XLALPrintInfo("analyze_series(): warning: cowardly refusing to analyze 0 samples... skipping series\n");
			return addpoint;
		}
	}

	for(i = 0; i < series->data->length - overlap; i += psdlength - overlap) {
		start = min(i, series->data->length - psdlength);

		interval = XLALCutREAL4TimeSeries(series, start, psdlength);

		XLALPrintInfo("analyze_series(): analyzing samples %zu -- %zu (%.9lf s -- %.9lf s)\n", start, start + interval->data->length, start * interval->deltaT, (start + interval->data->length) * interval->deltaT);

		*addpoint = XLALEPSearch(interval, params);
		while(*addpoint)
			addpoint = &(*addpoint)->next;
		if(xlalErrno) {
			XLALPrintError("analyze_series(): fatal error: XLALEPSearch() returned failure\n");
			exit(1);
		}

		XLALDestroyREAL4TimeSeries(interval);
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
 *                                Entry point
 *
 * ============================================================================
 */


int main(int argc, char *argv[])
{
	LALStatus stat;
	LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
	struct options *options;
	EPSearchParams params;
	LIGOTimeGPS epoch;
	LIGOTimeGPS boundepoch;
	size_t overlap;
	REAL4TimeSeries *series = NULL;
	SnglBurstTable *burstEvent = NULL;
	SnglBurstTable **EventAddPoint = &burstEvent;
	/* the ugly underscores are because some badger put global symbols
	 * in LAL by exactly these names. it's a mad house, a maad house */
	MetadataTable _process_table;
	MetadataTable _process_params_table;
	MetadataTable _search_summary_table;
	RandomParams *rparams = NULL;

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

	options = parse_command_line(argc, argv, &params, &_process_params_table);
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

	overlap = options->window_length - params.windowShift;
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
		 * Get the data.
		 */

		if(options->cache_filename) {
			/*
			 * Read from frame files
			 */

			series = get_time_series(options->cache_filename, options->channel_name, epoch, options->gps_end, options->max_series_length, options->calibrated, options->cal_high_pass);
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
			series = XLALCreateREAL4TimeSeries(options->channel_name, &epoch, 0.0, (REAL8) 1.0 / options->resample_rate, &lalADCCountUnit, length);
			if(!series) {
				XLALPrintError("%s: error: failure allocating data for Gaussian noise\n", argv[0]);
				exit(1);
			}
			if(!rparams)
				rparams = XLALCreateRandomParams(options->seed);
			gaussian_noise(series, options->noise_rms, rparams);
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

		if(options->sim_burst_filename || options->sim_inspiral_filename || options->sim_cache_filename) {
			COMPLEX8FrequencySeries *response;

			/* Create the response function (generates unity
			 * response if cache file is NULL). */
			response = generate_response(&stat, options->cal_cache_filename, options->ifo, options->channel_name, series->deltaT, series->epoch, series->data->length);
			if(!response)
				exit(1);

			/* perform injections */
			if(options->sim_burst_filename)
				add_burst_injections(&stat, options->sim_burst_filename, series, response);
			if(options->sim_inspiral_filename)
				add_inspiral_injections(&stat, options->sim_inspiral_filename, series, response);
			if(options->sim_cache_filename) {
				fprintf(stderr, "error:  \"sim\" code disabled\n");
				exit(1);
				/*add_sim_injections(&stat, series, response, options->max_series_length);*/
			}

			/*clean up */
			XLALDestroyCOMPLEX8FrequencySeries(response);
		}

		/*
		 * Add MDC injections into the time series if requested.
		 */

		if(options->mdc_cache_filename)
			add_mdc_injections(options->mdc_cache_filename, options->mdc_channel_name, series, epoch, options->gps_end, options->max_series_length);

		/*
		 * Condition the time series data.
		 */

		/* Scale the time series if calibrated data */
		if(options->calibrated) {
			const double scale = 1e10;
			unsigned i;
			for(i = 0; i < series->data->length; i++)
				series->data->data[i] *= scale;
		}

		if(XLALEPConditionData(series, options->high_pass, (REAL8) 1.0 / options->resample_rate, options->filter_corruption)) {
			XLALPrintError("%s: XLALEPConditionData() failed.\n", argv[0]);
			exit(1);
		}

		XLALPrintInfo("%s: %u samples (%.9f s) at GPS time %d.%09u s remain after conditioning\n", argv[0], series->data->length, series->data->length * series->deltaT, series->epoch.gpsSeconds, series->epoch.gpsNanoSeconds);

		/*
		 * Store the start and end times of the data that actually
		 * gets analyzed.
		 */

		if(!_search_summary_table.searchSummaryTable->out_start_time.gpsSeconds) {
			_search_summary_table.searchSummaryTable->out_start_time = series->epoch;
			XLALGPSAdd(&_search_summary_table.searchSummaryTable->out_start_time, series->deltaT * (options->window_length / 2 - params.windowShift));
		}
		_search_summary_table.searchSummaryTable->out_end_time = series->epoch;
		XLALGPSAdd(&_search_summary_table.searchSummaryTable->out_end_time, series->deltaT * (series->data->length - (options->window_length / 2 - params.windowShift)));

		/*
		 * Analyze the data
		 */

		EventAddPoint = analyze_series(EventAddPoint, series, options->psd_length, options->window_length, &params);

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
		XLALDestroyREAL4TimeSeries(series);
	}

	/*
	 * Unscale the h_{rss} estimates if calibrated data
	 */

	if(options->calibrated) {
		const double scale = 1e10;
		SnglBurstTable *event;
		for(event = burstEvent; event; event = event->next)
			event->amplitude /= scale;
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

	XLALDestroyRandomParams(rparams);
	XLALDestroyREAL4Window(params.window);
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
		LAL_CALL(LALCloseLIGOLwXMLFile(&stat, options->diagnostics), &stat);
	options_free(options);

	LALCheckMemoryLeaks();
	exit(0);
}
