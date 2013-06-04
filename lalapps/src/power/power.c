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


/*
 * ============================================================================
 *
 *                                  Preamble
 *
 * ============================================================================
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
#include <lal/BandPassTimeSeries.h>
#include <lal/Date.h>
#include <lal/EPSearch.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpChisq.h>
#include <lal/FrameCalibration.h>
#include <lal/LALFrStream.h>
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
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/LIGOLwXMLInspiralRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataBurstUtils.h>
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

#include <LALAppsVCSInfo.h>



#define PROGRAM_NAME "lalapps_power"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"


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
 * Return non-zero if the given integer is an integer power of 2.  The
 * trick here is that integer powers of 2 (and only integer powers of 2)
 * share exactly 0 bits with the integer 1 less than themselves, so we
 * check to see if that's the case.
 */


static int is_power_of_2(int x)
{
	return !((x - 1) & x);
}


/*
 * Return non-zero if the given double is an integer power of 2.
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
	/* XML file with list(s) of injections */
	char *injection_filename;

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

	options->cal_cache_filename = NULL;	/* default == disable */
	options->channel_name = NULL;	/* impossible */
	options->comment = default_comment;	/* default = "" */
	options->filter_corruption = -1;	/* impossible */
	options->mdc_cache_filename = NULL;	/* default == disable */
	options->mdc_channel_name = NULL;	/* impossible */
	options->noise_rms = -1;	/* default == disable */
	options->diagnostics = NULL;	/* default == disable */
	options->psd_length = 0;	/* impossible */
	options->psd_shift = 0;	/* impossible */
	options->resample_rate = 0;	/* impossible */
	options->seed = 0;	/* default == use system clock */
	options->max_series_length = 0;	/* default == disable */
	options->high_pass = -1;	/* impossible */
	options->max_event_rate = 0;	/* default == disable */
	options->output_filename = NULL;	/* impossible */
	options->injection_filename = NULL;	/* default == disable */
	options->cache_filename = NULL;	/* default == disable */
	options->confidence_threshold = XLAL_REAL8_FAIL_NAN;	/* impossible */
	options->bandwidth = 0;	/* impossible */
	options->flow = -1;	/* impossible */
	options->maxTileBandwidth = 0;	/* impossible */
	options->maxTileDuration = 0;	/* impossible */
	options->fractional_stride = 0;	/* impossible */
	options->window = NULL;	/* impossible */

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
"	[--calibration-cache <cache file name>]\n" \
"	 --channel-name <name>\n" \
"	 --confidence-threshold <confidence>\n" \
"	[--dump-diagnostics <XML file name>]\n" \
"	 --filter-corruption <samples>\n" \
"	 --frame-cache <cache file name>\n" \
"	[--gaussian-noise-rms <RMS amplitude>]\n", program);
	fprintf(stderr,
"	 --gps-end-time <seconds>\n" \
"	 --gps-start-time <seconds>\n" \
"	[--help]\n" \
"	 --high-pass <Hz>\n" \
"	[--injection-file <XML file name>]\n" \
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
	int option_index;
	int got_all_arguments = 1;
	int arg_is_missing;

	for(option_index = 0; long_options[option_index].name; option_index++) {
		switch(long_options[option_index].val) {
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
			arg_is_missing = 0;
			break;
		}
		if(arg_is_missing) {
			print_missing_argument(prog, long_options[option_index].name);
			got_all_arguments = 0;
		}
	}

	if(!!options->cache_filename + (options->noise_rms > 0.0) != 1) {
		XLALPrintError("%s: must provide exactly one of --frame-cache or --gaussian-noise-rms\n", prog);
		got_all_arguments = 0;
	}

	return got_all_arguments;
}


/*
 * Helper function for adding an entry to the process params table.
 */


static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param, const ProcessTable *process, const char *type, const char *param, const char *value)
{
	*proc_param = XLALCreateProcessParamsTableRow(process);
	snprintf((*proc_param)->program, sizeof((*proc_param)->program), PROGRAM_NAME);
	snprintf((*proc_param)->type, sizeof((*proc_param)->type), "%s", type);
	snprintf((*proc_param)->param, sizeof((*proc_param)->param), "--%s", param);
	snprintf((*proc_param)->value, sizeof((*proc_param)->value), "%s", value ? value : "");

	return &(*proc_param)->next;
}


#define ADD_PROCESS_PARAM(process, type) \
	do { paramaddpoint = add_process_param(paramaddpoint, process, type, long_options[option_index].name, optarg); } while(0)


/*
 * Parse the command line arguments.
 */


static struct options *parse_command_line(int argc, char *argv[], const ProcessTable *process, ProcessParamsTable **paramaddpoint)
{
	struct options *options;
	char msg[240];
	int args_are_bad = 0;
	int c;
	int option_index;
	struct option long_options[] = {
		{"bandwidth", required_argument, NULL, 'A'},
		{"injection-file", required_argument, NULL, 'P'},
		{"calibration-cache", required_argument, NULL, 'B'},
		{"channel-name", required_argument, NULL, 'C'},
		{"confidence-threshold", required_argument, NULL, 'g'},
		{"dump-diagnostics", required_argument, NULL, 'X'},
		{"filter-corruption", required_argument, NULL, 'j'},
		{"frame-cache", required_argument, NULL, 'G'},
		{"gaussian-noise-rms", required_argument, NULL, 'V'},
		{"gps-end-time", required_argument, NULL, 'K'},
		{"gps-start-time", required_argument, NULL, 'M'},
		{"help", no_argument, NULL, 'O'},
		{"high-pass", required_argument, NULL, 'o'},
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
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "int");
		break;

	case 'B':
		options->cal_cache_filename = optarg;
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'C':
		options->channel_name = optarg;
		memcpy(options->ifo, optarg, sizeof(options->ifo) - 1);
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'F':
		options->max_event_rate = atoi(optarg);
		if(options->max_event_rate < 0) {
			sprintf(msg, "must not be negative (%d specified)", options->max_event_rate);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "int");
		break;

	case 'G':
		options->cache_filename = optarg;
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'K':
		if(XLALStrToGPS(&options->gps_end, optarg, NULL)) {
			sprintf(msg, "range error parsing \"%s\"", optarg);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'M':
		if(XLALStrToGPS(&options->gps_start, optarg, NULL)) {
			sprintf(msg, "range error parsing \"%s\"", optarg);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'O':
		print_usage(argv[0]);
		exit(0);
		break;

	case 'P':
		options->injection_filename = optarg;
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'Q':
		options->flow = atof(optarg);
		if(options->flow < 0) {
			sprintf(msg, "must not be negative (%f Hz specified)", options->flow);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "float");
		break;

	case 'R':
		options->mdc_cache_filename = optarg;
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'S':
		options->mdc_channel_name = optarg;
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'V':
		options->noise_rms = atof(optarg);
		if(options->noise_rms <= 0.0) {
			sprintf(msg, "must be greater than 0 (%f specified)", options->noise_rms);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "float");
		break;

	case 'W':
		{
		int window_length = atoi(optarg);
		if((window_length < 4) || !is_power_of_2(window_length)) {
			sprintf(msg, "must be greater than or equal to 4 and a power of 2 (%i specified)", window_length);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		options->window = XLALCreateHannREAL8Window(window_length);
		if(!options->window) {
			sprintf(msg, "failure generating %d sample Hann window", window_length);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
			/* silence "missing required argument" message */
			options->window = (void *) 1;
		}
		ADD_PROCESS_PARAM(process, "int");
		}
		break;

	case 'X':
#if 0
		options->diagnostics = malloc(sizeof(*options->diagnostics));
		options->diagnostics->LIGOLwXMLStream = XLALOpenLIGOLwXMLFile(optarg);
		options->diagnostics->XLALWriteLIGOLwXMLArrayREAL8FrequencySeries = XLALWriteLIGOLwXMLArrayREAL8FrequencySeries;
		options->diagnostics->XLALWriteLIGOLwXMLArrayREAL8TimeSeries = XLALWriteLIGOLwXMLArrayREAL8TimeSeries;
		options->diagnostics->XLALWriteLIGOLwXMLArrayCOMPLEX16FrequencySeries = XLALWriteLIGOLwXMLArrayCOMPLEX16FrequencySeries;
#else
		sprintf(msg, "--dump-diagnostics given but diagnostic code not included at compile time");
		args_are_bad = 1;
#endif
		break;

	case 'Z':
		options->psd_length = atoi(optarg);
		ADD_PROCESS_PARAM(process, "int");
		break;

	case 'a':
		/*
		 * Convert the available RAM (in MebiBytes) to a
		 * guestimated limit on the length of a time series to read
		 * in.
		 */
		options->max_series_length = atoi(optarg) * 1024 * 1024 / (8 * sizeof(REAL8));
		if(options->max_series_length <= 0) {
			sprintf(msg, "must be greater than 0 (%i specified)", atoi(optarg));
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "int");
		break;

	case 'b':
		options->output_filename = optarg;
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'c':
		options->seed = atol(optarg);
		ADD_PROCESS_PARAM(process, "long");
		break;

	case 'e':
		options->resample_rate = atoi(optarg);
		if(options->resample_rate < 2 || options->resample_rate > 16384 || !is_power_of_2(options->resample_rate)) {
			sprintf(msg, "must be a power of 2 in the rage [2,16384] (%d specified)", options->resample_rate);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "int");
		break;

	case 'f':
		options->fractional_stride = atof(optarg);
		if(options->fractional_stride > 1 || !double_is_power_of_2(options->fractional_stride)) {
			sprintf(msg, "must be less than or equal to 1 and a power of 2 (%g specified)", options->fractional_stride);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "float");
		break;

	case 'g':
		options->confidence_threshold = atof(optarg);
		if(options->confidence_threshold < 0) {
			sprintf(msg, "must not be negative (%g specified)", options->confidence_threshold);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "float");
		break;

	case 'h':
		options->comment = optarg;
		ADD_PROCESS_PARAM(process, "string");
		break;

	case 'j':
		options->filter_corruption = atoi(optarg);
		if(options->filter_corruption < 0) {
			sprintf(msg, "must not be negative (%d specified)", options->filter_corruption);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "int");
		break;

	case 'l':
		options->maxTileBandwidth = atof(optarg);
		if((options->maxTileBandwidth <= 0) || !double_is_power_of_2(options->maxTileBandwidth)) {
			sprintf(msg, "must be greater than 0 and a power of 2 (%f specified)", options->maxTileBandwidth);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "float");
		break;

	case 'm':
		options->maxTileDuration = atof(optarg);
		if((options->maxTileDuration <= 0) || !double_is_power_of_2(options->maxTileDuration)) {
			sprintf(msg, "must be greater than 0 and a power of 2 (%f specified)", options->maxTileDuration);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "float");
		break;

	case 'o':
		options->high_pass = atof(optarg);
		if(options->high_pass < 0) {
			sprintf(msg, "must not be negative (%f Hz specified)", options->high_pass);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = 1;
		}
		ADD_PROCESS_PARAM(process, "float");
		break;

	/* option sets a flag */
	case 0:
		optarg = NULL;
		ADD_PROCESS_PARAM(process, "string");
		break;

	/* end of arguments */
	case -1:
		break;

	/* unrecognized option */
	case '?':
		print_usage(argv[0]);
		exit(1);

	/* missing argument for an option */
	case ':':
		print_usage(argv[0]);
		exit(1);
	} while(c != -1);

	/*
	 * Check for missing command line arguments.
	 */

	if(!all_required_arguments_present(argv[0], long_options, options))
		args_are_bad = 1;

	/*
	 * Check the order of the start and stop times.
	 */

	if(XLALGPSCmp(&options->gps_start, &options->gps_end) > 0) {
		XLALPrintError("%s: error: GPS start time > GPS stop time\n", argv[0]);
		args_are_bad = 1;
	}

	/*
	 * Warn if calibration cache is not provided that a unit response
	 * will be used for injections and h_rss.
	 */

	if(options->cal_cache_filename)
		XLALPrintWarning("warning: no calibration cache is provided:  software injections and hrss will be computed with unit response\n");

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
		int size = 1024, required_size;
		while(1) {
			options->output_filename = calloc(size, sizeof(*options->output_filename));
			if(!options->output_filename) {
				XLALPrintError("memory error");
				exit(1);
			}
			required_size = snprintf(options->output_filename, size, "%s-POWER_%s-%d-%d.xml", options->ifo, options->comment, options->gps_start.gpsSeconds, options->gps_end.gpsSeconds - options->gps_start.gpsSeconds);
			if(required_size < size)
				break;
			free(options->output_filename);
			size = required_size;
		}
		options->output_filename = realloc(options->output_filename, (required_size + 1) * sizeof(*options->output_filename));
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


static REAL8TimeSeries *get_time_series(const char *cachefilename, const char *chname, LIGOTimeGPS start, LIGOTimeGPS end, size_t lengthlimit)
{
	double duration = XLALGPSDiff(&end, &start);
	LALCache *cache;
	LALFrStream *stream;
	LALTYPECODE series_type;
	REAL8TimeSeries *series;

	if(duration < 0)
		XLAL_ERROR_NULL(XLAL_EINVAL);

	/*
	 * Open frame stream.
	 */

	cache = XLALCacheImport(cachefilename);
	if(!cache)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	stream = XLALFrStreamCacheOpen(cache);
	XLALDestroyCache(cache);
	if(!stream)
		XLAL_ERROR_NULL(XLAL_EFUNC);

	/*
	 * Turn on checking for missing data.
	 */

	XLALFrStreamSetMode(stream, LAL_FR_STREAM_VERBOSE_MODE);

	/*
	 * Get the data.
	 */

	series_type = XLALFrStreamGetTimeSeriesType(chname, stream);
	if((int) series_type < 0) {
		XLALFrStreamClose(stream);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	switch(series_type) {
	case LAL_S_TYPE_CODE: {
		/*
		 * Read data as single-precision and convert to double.
		 */

		REAL4TimeSeries *tmp = XLALFrStreamReadREAL4TimeSeries(stream, chname, &start, duration, lengthlimit);
		unsigned i;
		if(!tmp) {
			XLALFrStreamClose(stream);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
		series = XLALCreateREAL8TimeSeries(tmp->name, &tmp->epoch, tmp->f0, tmp->deltaT, &tmp->sampleUnits, tmp->data->length);
		if(!series) {
			XLALDestroyREAL4TimeSeries(tmp);
			XLALFrStreamClose(stream);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
		for(i = 0; i < tmp->data->length; i++)
			series->data->data[i] = tmp->data->data[i];
		XLALDestroyREAL4TimeSeries(tmp);
		break;
	}

	case LAL_D_TYPE_CODE:
		series = XLALFrStreamReadREAL8TimeSeries(stream, chname, &start, duration, lengthlimit);
		if(!series) {
			XLALFrStreamClose(stream);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}
		break;

	default:
		XLALPrintError("get_time_series(): error: invalid channel data type %d\n", series_type);
		XLALFrStreamClose(stream);
		XLAL_ERROR_NULL(XLAL_EINVAL);
	}

	/*
	 * Check for missing data.
	 */

	if(stream->state & LAL_FR_STREAM_GAP) {
		XLALPrintError("get_time_series(): error: gap in data detected between GPS times %d.%09u s and %d.%09u s\n", start.gpsSeconds, start.gpsNanoSeconds, end.gpsSeconds, end.gpsNanoSeconds);
		XLALDestroyREAL8TimeSeries(series);
		XLAL_ERROR_NULL(XLAL_EDATA);
	}

	/*
	 * Close stream.
	 */

	XLALFrStreamClose(stream);

	/*
	 * Verbosity.
	 */

	XLALPrintInfo("get_time_series(): read %u samples (%.9lf s) at GPS time %u.%09u s\n", series->data->length, series->data->length * series->deltaT, start.gpsSeconds, start.gpsNanoSeconds);

	return series;
}


/*
 * Condition the time series prior to analysis by the power code
 */


static int XLALEPConditionData(
	REAL8TimeSeries *series,
	REAL8 flow,
	REAL8 resampledeltaT,
	INT4 corruption
)
{
	const REAL8 epsilon = 1.0e-3;

	XLALPrintInfo("%s(): conditioning %u samples (%.9f s) at GPS time %d.%09u s ...\n", __func__, series->data->length, series->data->length * series->deltaT, series->epoch.gpsSeconds, series->epoch.gpsNanoSeconds);

	/*
	 * High-pass filter the time series.
	 */

	if(flow > 0.0) {
		PassBandParamStruc highpassParam;
		highpassParam.nMax = 8;
		highpassParam.f2 = flow;
		highpassParam.f1 = -1.0;
		highpassParam.a2 = 0.9;
		highpassParam.a1 = -1.0;
		if(XLALButterworthREAL8TimeSeries(series, &highpassParam))
			XLAL_ERROR(XLAL_EFUNC);
	}

	/*
	 * Resample the time series if necessary
	 */

	if(fabs(resampledeltaT - series->deltaT) / series->deltaT >= epsilon)
		if(XLALResampleREAL8TimeSeries(series, resampledeltaT))
			XLAL_ERROR(XLAL_EFUNC);

	/*
	 * Chop off the ends of the time series.
	 */

	if(!XLALShrinkREAL8TimeSeries(series, corruption, series->data->length - 2 * corruption))
		XLAL_ERROR(XLAL_EFUNC);

	/*
	 * Done.
	 */

	XLALPrintInfo("%s(): %u samples (%.9f s) at GPS time %d.%09u s remain after conditioning\n", __func__, series->data->length, series->data->length * series->deltaT, series->epoch.gpsSeconds, series->epoch.gpsNanoSeconds);

	return(0);
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
		LALUnit strainPerCount;
		unsigned i;

		XLALUnitDivide(&strainPerCount, &lalStrainUnit, &lalADCCountUnit);

		response = XLALCreateCOMPLEX8FrequencySeries(channel_name, &epoch, 0.0, deltaf, &strainPerCount, n);
		if(!response)
			XLAL_ERROR_NULL(XLAL_ENOMEM);

		XLALPrintInfo("generate_response(): generating unit response function\n");

		for(i = 0; i < response->data->length; i++)
			response->data->data[i] = 1.0;

		return response;
	} else {
                XLAL_ERROR_NULL(XLAL_EERR, "Calibration frames no longer supported");
	}
}


/*
 * ============================================================================
 *
 *                                 Injections
 *
 * ============================================================================
 */


/*
 * Load file
 */


struct injection_document {
	ProcessTable *process_table_head;
	ProcessParamsTable *process_params_table_head;
	SearchSummaryTable *search_summary_table_head;
	TimeSlide *time_slide_table_head;
	int has_sim_burst_table;
	SimBurst *sim_burst_table_head;
	int has_sim_inspiral_table;
	SimInspiralTable *sim_inspiral_table_head;
};


static void destroy_injection_document(struct injection_document *doc)
{
	if(doc) {
		XLALDestroyProcessTable(doc->process_table_head);
		XLALDestroyProcessParamsTable(doc->process_params_table_head);
		XLALDestroySearchSummaryTable(doc->search_summary_table_head);
		XLALDestroyTimeSlideTable(doc->time_slide_table_head);
		XLALDestroySimBurstTable(doc->sim_burst_table_head);
		while(doc->sim_inspiral_table_head) {
			SimInspiralTable *next = doc->sim_inspiral_table_head->next;
			XLALFree(doc->sim_inspiral_table_head);
			doc->sim_inspiral_table_head = next;
		}
	}
}


static struct injection_document *load_injection_document(const char *filename, LIGOTimeGPS start, LIGOTimeGPS end)
{
	struct injection_document *new;
	/* hard-coded speed hack.  only injections whose "times" are within
	 * this many seconds of the requested interval will be loaded */
	const double longest_injection = 600.0;

	new = malloc(sizeof(*new));
	if(!new)
		XLAL_ERROR_NULL(XLAL_ENOMEM);

	/*
	 * adjust start and end times
	 */

	XLALGPSAdd(&start, -longest_injection);
	XLALGPSAdd(&end, longest_injection);

	/*
	 * load required tables
	 */

	new->process_table_head = XLALProcessTableFromLIGOLw(filename);
	new->process_params_table_head = XLALProcessParamsTableFromLIGOLw(filename);
	new->search_summary_table_head = XLALSearchSummaryTableFromLIGOLw(filename);
	new->time_slide_table_head = XLALTimeSlideTableFromLIGOLw(filename);

	/*
	 * load optional sim_burst table
	 */

	new->has_sim_burst_table = XLALLIGOLwHasTable(filename, "sim_burst");
	if(new->has_sim_burst_table) {
		new->sim_burst_table_head = XLALSimBurstTableFromLIGOLw(filename, &start, &end);
	} else
		new->sim_burst_table_head = NULL;

	/*
	 * load optional sim_inspiral table
	 */

	new->has_sim_inspiral_table = XLALLIGOLwHasTable(filename, "sim_inspiral");
	if(new->has_sim_inspiral_table) {
		if(SimInspiralTableFromLIGOLw(&new->sim_inspiral_table_head, filename, start.gpsSeconds - 1, end.gpsSeconds + 1) < 0)
			new->sim_inspiral_table_head = NULL;
	} else
		new->sim_inspiral_table_head = NULL;

	/*
	 * did we get it all?
	 */

	if(
		!new->process_table_head ||
		!new->process_params_table_head ||
		!new->search_summary_table_head ||
		!new->time_slide_table_head ||
		(new->has_sim_burst_table && !new->sim_burst_table_head) ||
		(new->has_sim_inspiral_table && !new->sim_inspiral_table_head)
	) {
		destroy_injection_document(new);
		XLAL_ERROR_NULL(XLAL_EFUNC);
	}

	/*
	 * success
	 */

	return new;
}


/*
 * This is the function that gets called.
 */


static int add_xml_injections(REAL8TimeSeries *h, const struct injection_document *injection_document, COMPLEX8FrequencySeries *response)
{
	/*
	 * sim_burst
	 */

	if(injection_document->sim_burst_table_head) {
		XLALPrintInfo("%s(): computing sim_burst injections ...\n", __func__);
		if(XLALBurstInjectSignals(h, injection_document->sim_burst_table_head, injection_document->time_slide_table_head, NULL))
			XLAL_ERROR(XLAL_EFUNC);
		XLALPrintInfo("%s(): done\n", __func__);
	}

	/*
	 * sim_inspiral
	 */

	if(injection_document->sim_inspiral_table_head) {
		LALStatus stat;
		REAL4TimeSeries *mdc;
		unsigned i;

		mdc = XLALCreateREAL4TimeSeries(h->name, &h->epoch, h->f0, h->deltaT, &h->sampleUnits, h->data->length);
		if(!mdc)
			XLAL_ERROR(XLAL_EFUNC);
		memset(mdc->data->data, 0, mdc->data->length * sizeof(*mdc->data->data));
		memset(&stat, 0, sizeof(stat));

		XLALPrintInfo("%s(): computing sim_inspiral injections ...\n", __func__);
		LAL_CALL(LALFindChirpInjectSignals(&stat, mdc, injection_document->sim_inspiral_table_head, response), &stat);
		XLALPrintInfo("%s(): done\n", __func__);

		for(i = 0; i < h->data->length; i++)
			h->data->data[i] += mdc->data->data[i];
		XLALDestroyREAL4TimeSeries(mdc);
	}

	/*
	 * done
	 */

	return 0;
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
	LIGOTimeGPS stop;
	REAL8TimeSeries *mdc;
	unsigned i;

	XLALPrintInfo("%s(): mixing data from MDC frames\n", __func__);

	stop = series->epoch;
	XLALGPSAdd(&stop, series->data->length * series->deltaT);

	mdc = get_time_series(mdccachefile, channel_name, series->epoch, stop, 0);
	if(!mdc)
		XLAL_ERROR_NULL(XLAL_EFUNC);
	/* FIXME:  ARGH!!!  frame files cannot be trusted to provide units
	 * for their contents!  */
	series->sampleUnits = lalStrainUnit;

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


static SnglBurst **analyze_series(SnglBurst **addpoint, REAL8TimeSeries *series, int psd_length, int psd_shift, struct options *options)
{
	unsigned i;

	if((unsigned) psd_length > series->data->length) {
		XLALPrintWarning("%s(): warning: PSD length (%.9lf s) exceeds available data (%.9lf s), skipping series\n", __func__, psd_length * series->deltaT, series->data->length * series->deltaT);
		return addpoint;
	}

	for(i = 0; i + psd_length < series->data->length + psd_shift; i += psd_shift) {
		int errnum;
		int start = min(i, series->data->length - psd_length);
		REAL8TimeSeries *interval = XLALCutREAL8TimeSeries(series, start, psd_length);

		if(!interval)
			XLAL_ERROR_NULL(XLAL_EFUNC);

		XLALPrintInfo("%s(): ", __func__);
		XLALPrintProgressBar(i / (double) (series->data->length + psd_shift - psd_length));
		XLALPrintInfo(" complete\n");
		XLALPrintInfo("%s(): analyzing samples %zu -- %zu (%.9lf s -- %.9lf s)\n", __func__, start, start + interval->data->length, start * interval->deltaT, (start + interval->data->length) * interval->deltaT);

		XLAL_TRY(*addpoint = XLALEPSearch(
			options->diagnostics,
			interval,
			options->window,
			options->flow,
			options->bandwidth,
			options->confidence_threshold,
			options->fractional_stride,
			options->maxTileBandwidth,
			options->maxTileDuration
		), errnum);
		while(*addpoint)
			addpoint = &(*addpoint)->next;

		XLALDestroyREAL8TimeSeries(interval);

		if(errnum) {
			XLALSetErrno(errnum);
			XLAL_ERROR_NULL(XLAL_EFUNC);
		}

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


static void output_results(char *file, const ProcessTable *_process_table, const ProcessParamsTable *_process_params_table, const SearchSummaryTable *_search_summary_table, const SnglBurst *_sngl_burst_table)
{
	LIGOLwXMLStream *xml;

	xml = XLALOpenLIGOLwXMLFile(file);

	/*
	 * process table
	 */

	if(XLALWriteLIGOLwXMLProcessTable(xml, _process_table)) {
		/* FIXME:  error occured. do something smarter */
		exit(1);
	}

	/*
	 * process params table
	 */

	if(XLALWriteLIGOLwXMLProcessParamsTable(xml, _process_params_table)) {
		/* FIXME:  error occured. do something smarter */
		exit(1);
	}

	/*
	 * search summary table
	 */

	if(XLALWriteLIGOLwXMLSearchSummaryTable(xml, _search_summary_table)) {
		/* FIXME:  error occured. do something smarter */
		exit(1);
	}

	/*
	 * burst table
	 */

	if(XLALWriteLIGOLwXMLSnglBurstTable(xml, _sngl_burst_table)) {
		/* FIXME:  error occured. do something smarter */
		exit(1);
	}

	/*
	 * done
	 */

	XLALCloseLIGOLwXMLFile(xml);
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
	struct options *options;
	LIGOTimeGPS epoch;
	LIGOTimeGPS boundepoch;
	int overlap;
	REAL8TimeSeries *series = NULL;
	struct injection_document *injection_document = NULL;
	SnglBurst *_sngl_burst_table = NULL;
	SnglBurst **EventAddPoint = &_sngl_burst_table;
	/* the ugly underscores are because some badger put global symbols
	 * in LAL with exactly these names. it's a mad house, a maad house
	 * */
	ProcessTable *_process_table;
	ProcessParamsTable *_process_params_table = NULL;
	SearchSummaryTable *_search_summary_table;
	gsl_rng *rng = NULL;

	/*
	 * Command line
	 */

	lal_errhandler = LAL_ERR_EXIT;

	/*
	 * Create the process and process params tables.
	 *
	 * FIXME:  hard-coded process ID 9 is used to avoid conflicts with
	 * input injection files.  9 is high-enough for current use cases,
	 * but a better solution would be to set it to 0 and find a way to
	 * merge injection tables with existing process and process params
	 * tables from thsi process.
	 */

	_process_table = XLALCreateProcessTableRow();
	if(XLALPopulateProcessTable(_process_table, PROGRAM_NAME, LALAPPS_VCS_IDENT_ID, LALAPPS_VCS_IDENT_STATUS, LALAPPS_VCS_IDENT_DATE, 9))
		exit(1);

	XLALGPSTimeNow(&_process_table->start_time);

	/*
	 * Parse arguments and fill _process_params_table table.
	 */

	options = parse_command_line(argc, argv, _process_table, &_process_params_table);
	snprintf(_process_table->ifos, sizeof(_process_table->ifos), "%s", options->ifo);
	snprintf(_process_table->comment, sizeof(_process_table->comment), "%s", options->comment);

	/*
	 * Create the search summary table.  The number of nodes for a
	 * standalone job is always 1
	 */

	_search_summary_table = XLALCreateSearchSummaryTableRow(_process_table);
	snprintf(_search_summary_table->comment, sizeof(_search_summary_table->comment), "%s", options->comment);
	snprintf(_search_summary_table->ifos, sizeof(_search_summary_table->ifos), "%s", options->ifo);
	_search_summary_table->nnodes = 1;
	_search_summary_table->in_start_time = options->gps_start;
	_search_summary_table->in_end_time = options->gps_end;

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
	 * load injections if requested
	 */

	if(options->injection_filename) {
		injection_document = load_injection_document(options->injection_filename, options->gps_start, options->gps_end);
		if(!injection_document) {
			XLALPrintError("%s: error: failure reading injections file \"%s\"\n", argv[0], options->injection_filename);
			exit(1);
		}
	}

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

		XLALPrintInfo("%s: ", argv[0]);
		XLALPrintProgressBar(XLALGPSDiff(&epoch, &options->gps_start) / XLALGPSDiff(&boundepoch, &options->gps_start));
		XLALPrintInfo(" complete\n");

		/*
		 * Get the data.
		 */

		if(options->cache_filename) {
			/*
			 * Read from frame files
			 */

			series = get_time_series(options->cache_filename, options->channel_name, epoch, options->gps_end, options->max_series_length);
			if(!series) {
				XLALPrintError("%s: error: failure reading input data\n", argv[0]);
				exit(1);
			}
			/* FIXME:  ARGH!!!  frame files cannot be trusted
			 * to provide units for their contents!  assume ADC
			 * counts if a calibration cache has been provided,
			 * otherwise assume strain. */
			series->sampleUnits = options->cal_cache_filename ? lalADCCountUnit : lalStrainUnit;
		} else if(options->noise_rms > 0.0) {
			/*
			 * Synthesize Gaussian white noise.
			 */

			unsigned length = XLALGPSDiff(&options->gps_end, &epoch) * options->resample_rate;
			if(options->max_series_length)
				length = min(options->max_series_length, length);
			/* for units, assume ADC counts if a calibration
			 * cache has been provided, otherwise assume
			 * strain. */
			series = XLALCreateREAL8TimeSeries(options->channel_name, &epoch, 0.0, (REAL8) 1.0 / options->resample_rate, options->cal_cache_filename ? &lalADCCountUnit : &lalStrainUnit, length);
			if(!series) {
				XLALPrintError("%s: error: failure allocating data for Gaussian noise\n", argv[0]);
				exit(1);
			}
			if(!rng) {
				rng = gsl_rng_alloc(gsl_rng_mt19937);
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
					XLALPrintInfo("%s: using random number seed %lu\n", argv[0], seed);
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
		 * Add XML injections into the time series if requested.
		 */

		if(injection_document) {
			COMPLEX8FrequencySeries *response;

			/* Create the response function (generates unity
			 * response if cache file is NULL). */
			response = generate_response(options->cal_cache_filename, options->channel_name, series->deltaT, series->epoch, series->data->length);
			if(!response)
				exit(1);

			/* perform XML injections */
			if(add_xml_injections(series, injection_document, response))
				exit(1);

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

		if(!_search_summary_table->out_start_time.gpsSeconds) {
			_search_summary_table->out_start_time = series->epoch;
			XLALGPSAdd(&_search_summary_table->out_start_time, series->deltaT * options->window_pad);
		}
		_search_summary_table->out_end_time = series->epoch;
		XLALGPSAdd(&_search_summary_table->out_end_time, series->deltaT * (series->data->length - options->window_pad));

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
		 * computed earlier provides the additional overlap needed
		 * for the time-frequency tiling.
		 */

		XLALGPSAdd(&epoch, (series->data->length - overlap) * series->deltaT);
		XLALDestroyREAL8TimeSeries(series);
	}

	/*
	 * Sort the events, and assign IDs.
	 */

	XLALSortSnglBurst(&_sngl_burst_table, XLALCompareSnglBurstByStartTimeAndLowFreq);
	XLALSnglBurstAssignIDs(_sngl_burst_table, _process_table->process_id, 0);

	/*
	 * Check event rate limit.
	 */

	_search_summary_table->nevents = XLALSnglBurstTableLength(_sngl_burst_table);
	if((options->max_event_rate > 0) && (_search_summary_table->nevents > XLALGPSDiff(&_search_summary_table->out_end_time, &_search_summary_table->out_start_time) * options->max_event_rate)) {
		XLALPrintError("%s: event rate limit exceeded!", argv[0]);
		exit(1);
	}

	/*
	 * Output the results.
	 */

	XLALGPSTimeNow(&_process_table->end_time);
	output_results(options->output_filename, _process_table, _process_params_table, _search_summary_table, _sngl_burst_table);

	/*
	 * Final cleanup.
	 */

	if(rng)
		gsl_rng_free(rng);

	XLALDestroyProcessTable(_process_table);
	XLALDestroyProcessParamsTable(_process_params_table);
	XLALDestroySearchSummaryTable(_search_summary_table);
	XLALDestroySnglBurstTable(_sngl_burst_table);
	destroy_injection_document(injection_document);

	if(options->diagnostics)
		XLALCloseLIGOLwXMLFile(options->diagnostics->LIGOLwXMLStream);
	options_free(options);

	LALCheckMemoryLeaks();
	exit(0);
}
