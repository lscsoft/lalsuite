#include <getopt.h>
#include <lalapps.h>
#include <processtable.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/EPSearch.h>
#include <lal/ExcessPower.h>
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
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/PrintFTSeries.h>
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

int snprintf(char *str, size_t size, const char *format, ...);
int vsnprintf(char *str, size_t size, const char *format, va_list ap);

NRCSID( POWERC, "power $Id$");
RCSID( "power $Id$");

#define PROGRAM_NAME "power"
#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"

#include <config.h>
#ifndef HAVE_LIBLALFRAME
int main(void)
{
	fputs("Disabled: LALApps compiled with non-frame-enabled LAL\n", stderr);
	exit(77);
}
#else

#define TRUE       1
#define FALSE      0
#define SCALE      1e10

/*
 * ============================================================================
 *                                Global Data
 * ============================================================================
 */

/* Parameters from command line */
static struct {
	CHAR *calCacheFile;         /* name of the calibration cache file  */
 	CHAR *simCacheFile;         /* name of the sim waveform cache file */
 	CHAR *simdirname;           /* name of the dir. with the sim. wave */

	int cluster;                /* TRUE == perform clustering          */
	int estimateHrss;           /* TRUE == estimate hrss               */
	CHAR *comment;              /* user comment                        */
	int FilterCorruption;       /* samples corrupted by conditioning   */
	INT4 maxSeriesLength;       /* RAM-limited input length            */
	REAL4 noiseAmpl;            /* gain factor for white noise         */
	INT4 printData;
	size_t PSDAverageLength;    /* number of samples to use for PSD    */
	INT4 ResampleRate;          /* sample rate after resampling        */
	INT4 seed;                  /* set non-zero to generate noise      */
	LIGOTimeGPS startEpoch;     /* gps start time                      */
	LIGOTimeGPS stopEpoch;      /* gps stop time                       */
	int verbose;
	int calibrated;             /* input is double-precision h(t)      */
	int getcaltimeseries;       /* flag to read in double time series  */ 
	REAL8 high_pass;            /* conditioning high pass freq (Hz)    */
	REAL8 cal_high_pass;        /* double->single high pass freq (Hz)  */
	int bandwidth;
	int max_event_rate;         /* safety valve (Hz), 0 == disable     */
	UINT4 windowLength;
	WindowType windowType;
	char *channelName;
	char *mdcchannelName;         /* mdc signal only channnel info       */
} options;
 
/* global variables */
CHAR ifo[3];                        /* two character interferometer        */
CHAR *cachefile;                    /* name of file with frame cache info  */
CHAR *dirname;                      /* name of directory with frames       */

/* data conditioning parameters */
CHAR *burstInjectionFile;           /* file with list of burst injections  */
CHAR *inspiralInjectionFile;        /* file with list of burst injections  */
CHAR *simInjectionFile;             /* file with list of sim injections  */
CHAR *mdcCacheFile;                 /* name of mdc signal cache file       */
CHAR *mdcDirName;                 /* name of mdc signal cache file       */

/*
 * ============================================================================
 *                               Misc Utilities
 * ============================================================================
 */

/*
 * Return the smaller of two size_t variables.
 */

static size_t min(size_t a, size_t b)
{
	return(a < b ? a : b);
}


/*
 * Compare two GPS times.
 */

static int CompareGPS(LALStatus *stat, LIGOTimeGPS *gps1, LIGOTimeGPS *gps2)
{
	LALGPSCompareResult result;
	LAL_CALL(LALCompareGPS(stat, &result, gps1, gps2), stat);
	return(result);
}


/*
 * Return the difference between two GPS times as REAL8.
 */

static REAL8 DeltaGPStoFloat(LALStatus *stat, LIGOTimeGPS *stop, LIGOTimeGPS *start)
{
	REAL8 d;
	LAL_CALL(LALDeltaFloatGPS(stat, &d, stop, start), stat);
	return(d);
}


/*
 * Round length down so that an integer number of intervals of length
 * block_length, each shifted by block_shift with respect to its neighbours,
 * fits into the result.  This is used to ensure that an integer number of
 * analysis windows fits into the PSD length, and also that an integer number
 * of PSD lenghts fits into the RAM limit length.
 */

static size_t block_commensurate(
	size_t length,
	size_t block_length,
	size_t block_shift
)
{
	volatile size_t blocks = (length - block_length) / block_shift;

	return(length < block_length ? 0 : blocks * block_shift + block_length);
}


/*
 * Return TRUE if the given integer is an integer power of 2.  The trick
 * here is that integer powers of 2 (and only integer powers of 2) share
 * exactly 0 bits with the integer 1 less than themselves, so we check to
 * see if that's the case.
 */

static int is_power_of_2(int x)
{
	return(!((x - 1) & x));
}


/*
 * Count the events in a SnglBurstTable.
 */

static int SnglBurstTableLength(SnglBurstTable *list)
{
	int i;

	for(i = 0; list; i++)
		list = list->next;

	return(i);
}


/*
 * ============================================================================
 *                       Initialize and parse arguments
 * ============================================================================
 */

static void print_usage(char *program)
{
	fprintf(stderr,
"Usage:  %s <option> [...]\n" \
"The following options are recognized.  Options not surrounded in [] are\n" \
"required.\n" \
"	 --bandwidth <bandwidth>\n" \
"	[--calibrated-data <high pass frequency>]\n" \
"	[--calibration-cache <cache file>]\n" \
"	 --channel-name <string>\n" \
"	[--cluster]\n" \
"	[--debug-level <level>]\n" \
"	[--estimatehrss]\n" \
"	[--max-event-rate <Hz>]\n" \
"	 --filter-corruption <samples>\n" \
"	 --frame-cache <cache file>\n" \
"	 --frame-dir <directory>\n" \
"	 --gps-end-time <seconds>\n" \
"	 --gps-start-time <seconds>\n" \
"	[--help]\n" \
"	 --high-pass <high pass frequency>\n" \
"	[--burstinjection-file <file name>]\n" \
"	[--inspiralinjection-file <file name>]\n" \
"	 --low-freq-cutoff <Hz>\n" \
"	 --max-tileband <Hz>\n" \
"	 --max-tileduration <samples>\n" \
"	[--mdc-cache <cache file>]\n" \
"	[--mdc-channel <channel name>]\n" \
"	[--noise-amplitude <amplitude>]\n" \
"	[--printData]\n" \
"	[--printSpectrum]\n" \
"	 --psd-average-method <method>\n" \
"	 --psd-average-points <samples>\n" \
"	[--ram-limit <MebiBytes>]\n" \
"	 --resample-rate <Hz>\n" \
"       [--sim-cache <sim cache file>]\n" \
"       [--sim-seconds <sim seconds>]\n" \
"	[--siminjection-file <file name>]\n" \
"	[--seed <seed>]\n" \
"	 --tile-stride-fraction <fraction>\n" \
"	 --lnthreshold <ln threshold>\n" \
"	[--useoverwhitening]\n" \
"	[--user-tag <comment>]\n" \
"	[--verbose]\n" \
"	 --window <window>\n" \
"	 --window-length <samples>\n" \
"	 --window-shift <samples>\n", program);
}

static void print_bad_argument(const char *prog, const char *arg, const char *msg)
{
	fprintf(stderr, "%s: error: invalid argument for --%s: %s\n", prog, arg, msg);
}

static void print_missing_argument(const char *prog, const char *arg)
{
	fprintf(stderr, "%s: error: --%s not specified\n", prog, arg);
}

static void print_alloc_fail(const char *prog, const char *msg)
{
	fprintf(stderr, "%s: error: memory allocation failure %s\n", prog, msg);
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

static int check_for_missing_parameters(char *prog, struct option *long_options, EPSearchParams *params)
{
	int index;
	int got_all_arguments = TRUE;
	int arg_is_missing;

	for(index = 0; long_options[index].name; index++) {
		switch(long_options[index].val) {
			case 'A':
			arg_is_missing = !options.bandwidth;
			break;

			case 'C':
			arg_is_missing = !options.channelName;
			break;

			case 'K':
			arg_is_missing = !XLALGPStoINT8(&options.stopEpoch);
			break;

			case 'M':
			arg_is_missing = !XLALGPStoINT8(&options.startEpoch);
			break;

			case 'Q':
			arg_is_missing = params->tfPlaneParams.flow < 0.0;
			break;

			case 'W':
			arg_is_missing = !options.windowLength;
			break;

			case 'Y':
			arg_is_missing = params->method == (unsigned) -1;
			break;

			case 'Z':
			arg_is_missing = !options.PSDAverageLength;
			break;

			case 'd':
			arg_is_missing = !params->windowShift;
			break;

			case 'e':
			arg_is_missing = !options.ResampleRate;
			break;

			case 'f':
			arg_is_missing = !params->tfTilingInput.inv_fractional_stride;
			break;

			case 'g':
			arg_is_missing = params->lnalphaThreshold == XLAL_REAL8_FAIL_NAN;
			break;

			case 'i':
			arg_is_missing = options.windowType == NumberWindowTypes;
			break;

			case 'j':
			arg_is_missing = options.FilterCorruption < 0;
			break;

			case 'l':
			arg_is_missing = !params->tfTilingInput.maxTileBandwidth;
			break;

			case 'o':
			arg_is_missing = options.high_pass < 0.0;
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

	if(!cachefile && !dirname && (options.noiseAmpl < 0.0)) {
		fprintf(stderr, "%s: must provide at least one of --frame-cache, --frame-dir or --noise-amplitude\n", prog);
		arg_is_missing = TRUE;
	}
	
	return(got_all_arguments);
}


void parse_command_line_debug(
	int argc,
	char *argv[]
)
{
	int c;
	int option_index;
	struct option long_options[] = {
		{"debug-level",         required_argument, NULL,           'D'},
		{NULL, 0, NULL, 0}
	};

	/*
	 * Find and parse only the debug level command line options.  Must jump
	 * through this hoop because we cannot call set_debug_level() after any
	 * calls to LALMalloc() and friends.
	 */

	opterr = 0;	/* silence error messages */
	optind = 0;	/* start scanning from argv[0] */
	do switch(c = getopt_long(argc, argv, "-", long_options, &option_index)) {

		case 'D':
		/* only set the debug level in this pass */
		set_debug_level(optarg);
		default:
		break;
	} while(c != -1);
}



#define ADD_PROCESS_PARAM(type) \
	do { paramaddpoint = add_process_param(paramaddpoint, type, long_options[option_index].name, optarg); } while(0)

void parse_command_line( 
	int argc, 
	char *argv[], 
	EPSearchParams *params,
	MetadataTable *procparams 
)
{ 
	char msg[240];
	int args_are_bad = FALSE;
	int printSpectrum;
	int useoverwhitening;
	int c;
	int option_index;
	int ram;
	ProcessParamsTable **paramaddpoint = &procparams->processParamsTable;
	struct option long_options[] = {
		{"bandwidth",           required_argument, NULL,           'A'},
		{"calibrated-data",	required_argument, NULL,           'J'},
		{"calibration-cache",   required_argument, NULL,           'B'},
		{"channel-name",        required_argument, NULL,           'C'},
		{"cluster",             no_argument, &options.cluster,    TRUE},
		{"estimatehrss",        no_argument, &options.estimateHrss, TRUE},
		{"debug-level",         required_argument, NULL,           'D'},
		{"max-event-rate",      required_argument, NULL,           'F'},
		{"filter-corruption",	required_argument, NULL,           'j'},
		{"frame-cache",         required_argument, NULL,           'G'},
		{"frame-dir",           required_argument, NULL,           'H'},
		{"gps-end-time",        required_argument, NULL,           'K'},
		{"gps-start-time",      required_argument, NULL,           'M'},
		{"help",                no_argument,       NULL,           'O'},
		{"high-pass",           required_argument, NULL,           'o'},
		{"burstinjection-file", required_argument, NULL,           'P'},
		{"inspiralinjection-file", required_argument, NULL,        'I'},
		{"low-freq-cutoff",     required_argument, NULL,           'Q'},
		{"max-tileband",        required_argument, NULL,           'l'},
		{"max-tileduration",    required_argument, NULL,           'm'},
		{"mdc-cache",           required_argument, NULL,           'R'},
		{"mdc-channel",         required_argument, NULL,           'S'},
		{"noise-amplitude",     required_argument, NULL,           'V'},
		{"printData",           no_argument, &options.printData,  TRUE},
		{"printSpectrum",       no_argument, &printSpectrum,      TRUE},
		{"psd-average-method",  required_argument, NULL,           'Y'},
		{"psd-average-points",  required_argument, NULL,           'Z'},
		{"ram-limit",           required_argument, NULL,           'a'},
		{"resample-rate",       required_argument, NULL,           'e'},
		{"sim-cache",           required_argument, NULL,           'q'},
		{"siminjection-file",   required_argument, NULL,           't'},
		{"seed",                required_argument, NULL,           'c'},
		{"tile-stride-fraction", required_argument, NULL,           'f'},
		{"lnthreshold",         required_argument, NULL,           'g'},
		{"useoverwhitening",    no_argument, &useoverwhitening,   TRUE},  
		{"user-tag",            required_argument, NULL,           'h'},
		{"verbose",             no_argument, &options.verbose,    TRUE},
		{"window",              required_argument, NULL,           'i'},
		{"window-length",       required_argument, NULL,           'W'},
		{"window-shift",        required_argument, NULL,           'd'},
		{NULL, 0, NULL, 0}
	};

	/*
	 * Set parameter defaults.
	 */

	params->lnalphaThreshold = XLAL_REAL8_FAIL_NAN;	/* impossible */
	params->method = -1;	/* impossible */
	params->tfPlaneParams.flow = -1.0;	/* impossible */
	params->tfTilingInput.maxTileBandwidth = 0;  /* impossible */
	params->tfTilingInput.inv_fractional_stride = 0;	/* impossible */
	params->windowShift = 0;	/* impossible */

	options.bandwidth = 0;	/* impossible */
	options.calCacheFile = NULL;	/* default */
	options.channelName = NULL;	/* impossible */
	options.cluster = FALSE;	/* default */
	options.estimateHrss = FALSE;	/* default */
	options.comment = "";		/* default */
	options.FilterCorruption = -1;	/* impossible */
	options.mdcchannelName = NULL;	/* default */
	options.noiseAmpl = -1.0;	/* default */
	options.printData = FALSE;	/* default */
	options.PSDAverageLength = 0;	/* impossible */
	options.ResampleRate = 0;	/* impossible */
	options.seed = 1;	        /* default */
	options.verbose = FALSE;	/* default */
	options.calibrated = FALSE;	/* default */
	options.getcaltimeseries = FALSE; /* default */
	options.high_pass = -1.0;	/* impossible */
	options.cal_high_pass = -1.0;	/* impossible */
	options.max_event_rate = 0;	/* default */
	options.windowType = NumberWindowTypes;	/* impossible */
	options.windowLength = 0;	/* impossible */
	XLALINT8toGPS(&options.startEpoch, 0);	/* impossible */
	XLALINT8toGPS(&options.stopEpoch, 0);	/* impossible */

	options.simCacheFile = NULL;	/* default */
	options.simdirname = NULL;      /* default */

	ram = 0;	/* default */

	cachefile = NULL;	        /* default */
	dirname = NULL;	                /* default */
	memset(ifo, 0, sizeof(ifo));	/* default */
	burstInjectionFile = NULL;	/* default */
	simInjectionFile = NULL;	/* default */
	inspiralInjectionFile = NULL;	/* default */
	mdcCacheFile = NULL;	        /* default */
	mdcDirName = NULL;	        /* default */
	printSpectrum = FALSE;	        /* default */
	useoverwhitening = FALSE;       /* default */

	/*
	 * Parse command line.
	 */

	opterr = 1;	/* enable error messages */
	optind = 0;	/* start scanning from argv[0] */
	do switch(c = getopt_long(argc, argv, "", long_options, &option_index)) {
		case 'A':
		options.bandwidth = atoi(optarg);
		if(options.bandwidth <= 0 || !is_power_of_2(options.bandwidth) ) {
			sprintf(msg, "must be > 0 and a power of 2(%i specified)", options.bandwidth);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

		case 'B':
		options.calCacheFile = optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 'C':
		options.channelName = optarg;
		memcpy(ifo, optarg, sizeof(ifo) - 1);
		ADD_PROCESS_PARAM("string");
		break;

		case 'D':
		/* only add --debug-level to params table in this pass */
		ADD_PROCESS_PARAM("string");
		break;

		case 'F':
		options.max_event_rate = atoi(optarg);
		ADD_PROCESS_PARAM("int");
		break;

		case 'G':
		cachefile = optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 'H':
		dirname =  optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 'I':
		inspiralInjectionFile = optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 'J':
		options.calibrated = TRUE;
		options.cal_high_pass = atof(optarg);
		if(options.cal_high_pass < 0.0) {
			sprintf(msg, "must not be negative (%f Hz specified)", options.cal_high_pass);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("double");
		break;

		case 'K':
		if(XLALStrToGPS(&options.stopEpoch, optarg, NULL)) {
			sprintf(msg, "range error parsing \"%s\"", optarg);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		} else if(XLALGPStoINT8(&options.stopEpoch) < LAL_INT8_C(441417609000000000) || XLALGPStoINT8(&options.stopEpoch) > LAL_INT8_C(999999999000000000)) {
			sprintf(msg, "must be in the range [Jan 01 1994 00:00:00 UTC, Sep 14 2011 01:46:26 UTC] (%d.%09d specified)", options.stopEpoch.gpsSeconds, options.stopEpoch.gpsNanoSeconds);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("string");
		break;

		case 'M':
		if(XLALStrToGPS(&options.startEpoch, optarg, NULL)) {
			sprintf(msg, "range error parsing \"%s\"", optarg);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		} else if(XLALGPStoINT8(&options.startEpoch) < LAL_INT8_C(441417609000000000) || XLALGPStoINT8(&options.startEpoch) > LAL_INT8_C(999999999000000000)) {
			sprintf(msg, "must be in the range [Jan 01 1994 00:00:00 UTC, Sep 14 2011 01:46:26 UTC] (%d.%09d specified)", options.startEpoch.gpsSeconds, options.startEpoch.gpsNanoSeconds);
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
		burstInjectionFile = optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 'Q':
		params->tfPlaneParams.flow = atof(optarg);
		if((params->tfPlaneParams.flow < 0.0)) {
			sprintf(msg,"must not be negative (%f Hz specified)", params->tfPlaneParams.flow);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

		case 'R':
		mdcCacheFile = optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 'S':
		options.mdcchannelName = optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 'V':
		options.noiseAmpl = atof(optarg);
		if(options.noiseAmpl <= 0.0) {
			sprintf(msg, "must be > 0.0 (%f specified)", options.noiseAmpl);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

		case 'W':
		options.windowLength = atoi(optarg);
		if((options.windowLength <= 0) || !is_power_of_2(options.windowLength)) {
			sprintf(msg, "must be a power of 2 greater than 0 (%i specified)", options.windowLength);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
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
		options.PSDAverageLength = atoi(optarg);
		ADD_PROCESS_PARAM("int");
		break;

		case 'a':
		ram = atoi(optarg);
		if(ram <= 0) {
			sprintf(msg, "must be > 0 (%i specified)", ram);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

		case 'c':
		options.seed = atoi(optarg);
		if(options.seed <= 0) {
			sprintf(msg, "must be > 0 (%i specified)", options.seed);
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
		options.ResampleRate = (INT4) atoi(optarg);
		if(options.ResampleRate < 2 || options.ResampleRate > 16384 || !is_power_of_2(options.ResampleRate)) {
			sprintf(msg, "must be a power of 2 in the rage [2,16384] (%d specified)", options.ResampleRate);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

		case 'f':
		params->tfTilingInput.inv_fractional_stride = 1.0 / atof(optarg);
		if(params->tfTilingInput.inv_fractional_stride < 0 ||  !is_power_of_2(params->tfTilingInput.inv_fractional_stride)) {
			sprintf(msg, "must be 2^-n, n integer, (%g specified)", 1.0 / params->tfTilingInput.inv_fractional_stride);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

		case 'g':
		params->lnalphaThreshold = atof(optarg);
		ADD_PROCESS_PARAM("float");
		break;

		case 'h':
		options.comment = optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 'i':
		options.windowType = atoi(optarg);
		if(options.windowType >= NumberWindowTypes) {
			sprintf(msg, "must be < %d (%i specified)", NumberWindowTypes, options.windowType);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

		case 'j':
		options.FilterCorruption = atoi(optarg);
		if(options.FilterCorruption < 0) {
			sprintf(msg, "must be >= 0 (%d specified)", options.FilterCorruption);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

		case 'l':
		params->tfTilingInput.maxTileBandwidth = atof(optarg);
		params->tfPlaneParams.deltaT = 1 / (2 * params->tfTilingInput.maxTileBandwidth);
		if((params->tfTilingInput.maxTileBandwidth <= 0) || !is_power_of_2(params->tfTilingInput.maxTileBandwidth)) {
			sprintf(msg,"must be a power of 2 greater than 0 (%f specified)",params->tfTilingInput.maxTileBandwidth);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

		case 'm':
		params->tfTilingInput.maxTileDuration = atof(optarg);
		params->tfPlaneParams.deltaF = 1 / (2 * params->tfTilingInput.maxTileDuration);
		if((params->tfTilingInput.maxTileDuration > 1.0) || !is_power_of_2(1/params->tfTilingInput.maxTileDuration)) {
			sprintf(msg,"must be a power of 2 not greater than 1.0 (%f specified)", params->tfTilingInput.maxTileDuration);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

		case 'o':
		options.high_pass = atof(optarg);
		if(options.high_pass < 0.0) {
			sprintf(msg, "must not be negative (%f Hz specified)", options.high_pass);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

		case 'q':
		options.simCacheFile = optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 't':
		  simInjectionFile = optarg;
		  ADD_PROCESS_PARAM("string");
		  break;

		/* option sets a flag */
		case 0:
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
	 * Make sure windowShift and windowLength are OK.
	 */

	if(options.windowLength < 2 * params->windowShift) {
		sprintf(msg, "must be >= 2 * --window-shift = %u (%u specified)", 2 * params->windowShift, options.windowLength);
		print_bad_argument(argv[0], "window-length", msg);
		args_are_bad = TRUE;
	}

	/*
	 * Check the order of the start and stop times.
	 */

	if(XLALGPStoINT8(&options.startEpoch) > XLALGPStoINT8(&options.stopEpoch)) {
		fprintf(stderr, "%s: error: GPS start time > GPS stop time\n", argv[0]);
		args_are_bad = TRUE;
	}

	/*
	 * Convert the amount of available RAM to a limit on the length of a
	 * time series to read in.
	 */

	options.maxSeriesLength = (ram * 1024 * 1024) / (4 * sizeof(REAL4));

	/*
	 * Check for missing parameters
	 */

	if(!check_for_missing_parameters(argv[0], long_options, params))
		args_are_bad = TRUE;

	/*
	 * Exit if anything was wrong with the command line.
	 */

	if(args_are_bad)
		exit(1);

	/*
	 * Generate time-domain window function.
	 */

	params->window = XLALCreateREAL4Window(options.windowLength, options.windowType, 0.0);
	if(!params->window) {
		fprintf(stderr, "%s: failure generating time-domain window\n", argv[0]);
		exit(1);
	}

	/*
	 * Ensure PSDAverageLength is comensurate with the analysis window
	 * length and its shift.
	 */

	options.PSDAverageLength = block_commensurate(options.PSDAverageLength, options.windowLength, params->windowShift);

	/*
	 * Ensure calibration cache is given when the hrss estimation option is 
	 * turned on
	 */
	if(options.estimateHrss && !options.calCacheFile) {
		fprintf(stderr, "estimate hrss option is turned on but calibration file is missing\n");
		exit(1);
	}

	/*
	 * Ensure RAM limit is comensurate with the PSDAverageLength and
	 * its shift.
	 */

	options.maxSeriesLength = block_commensurate(options.maxSeriesLength, options.PSDAverageLength, options.PSDAverageLength - (options.windowLength - params->windowShift));

	/*
	 * Sanitize filter frequencies.
	 */

	if(options.cal_high_pass > params->tfPlaneParams.flow)
		fprintf(stderr, "%s: warning: calibrated data quantization high-pass frequency (%f Hz) greater than TF plane low frequency (%f Hz)\n", argv[0], options.cal_high_pass, params->tfPlaneParams.flow);

	if(options.high_pass > params->tfPlaneParams.flow - 10.0)
		fprintf(stderr, "%s: warning: data conditioning high-pass frequency (%f Hz) greater than 10 Hz below TF plane low frequency (%f Hz)\n", argv[0], options.high_pass, params->tfPlaneParams.flow);

	/*
	 * Miscellaneous chores.
	 */

	params->tfPlaneParams.timeBins = (options.windowLength / 2) / (options.ResampleRate * params->tfPlaneParams.deltaT);
	params->tfPlaneParams.freqBins = options.bandwidth / params->tfPlaneParams.deltaF;

	params->printSpectrum = printSpectrum;
	params->useOverWhitening = useoverwhitening;

	if(options.verbose) {
		fprintf(stderr, "%s: using --psd-average-points %zu\n", argv[0], options.PSDAverageLength);
		fprintf(stderr, "%s: available RAM limits analysis to %d samples\n", argv[0], options.maxSeriesLength);
		fprintf(stderr, "%s: time-frequency plane has %d time bins by %d frequency bins\n", argv[0], params->tfPlaneParams.timeBins, params->tfPlaneParams.freqBins);
	}
}


/*
 * ============================================================================
 *                                   Input
 * ============================================================================
 */

static REAL4TimeSeries *get_calibrated_data(
	FrStream *stream,
	const char *chname,
	LIGOTimeGPS *start,
	double duration,
	size_t lengthlimit
)
{
	static const char func[] = "get_calibrated_data";
	REAL8TimeSeries *calibrated;
	REAL4TimeSeries *series;
	PassBandParamStruc highpassParam;
	unsigned int i;

	/* retrieve calibrated data as REAL8 time series */
	calibrated = XLALFrReadREAL8TimeSeries(stream, chname, start, duration, lengthlimit);
	if(!calibrated)
		XLAL_ERROR_NULL(func, XLAL_EFUNC);

	/* high pass filter before casting REAL8 to REAL4 */
	highpassParam.nMax = 4;
	highpassParam.f2 = options.cal_high_pass;
	highpassParam.f1 = -1.0;
	highpassParam.a2 = 0.9;
	highpassParam.a1 = -1.0;
	if(XLALButterworthREAL8TimeSeries(calibrated, &highpassParam) < 0) {
		XLALDestroyREAL8TimeSeries(calibrated);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}

	/* copy data into a REAL4 time series */
	series = XLALCreateREAL4TimeSeries(calibrated->name, &calibrated->epoch, calibrated->f0, calibrated->deltaT, &calibrated->sampleUnits, calibrated->data->length);
	if(!series) {
		XLALDestroyREAL8TimeSeries(calibrated);
		XLAL_ERROR_NULL(func, XLAL_EFUNC);
	}
	for(i = 0; i < series->data->length; i++)
		series->data->data[i] = calibrated->data->data[i];

	XLALDestroyREAL8TimeSeries(calibrated);
	return(series);
}

static REAL4TimeSeries *get_time_series(
	LALStatus *stat,
	const char *dirname,
	const char *cachefile,
	const char *chname,
	LIGOTimeGPS start,
	LIGOTimeGPS end,
	size_t lengthlimit,
	int getcaltimeseries
)
{
	REAL4TimeSeries *series;
	FrStream *stream = NULL;
	FrCache *frameCache = NULL;
	double duration = DeltaGPStoFloat(stat, &end, &start);

	/* Open frame stream */
	if(cachefile && dirname && options.verbose)
		fprintf(stderr, "get_time_series(): warning: --frame-cache ignored (using --frame-dir)\n");
	if(dirname)
		LAL_CALL(LALFrOpen(stat, &stream, dirname, "*.gwf"), stat);
	else if(cachefile) {
		LAL_CALL(LALFrCacheImport(stat, &frameCache, cachefile), stat);
		LAL_CALL(LALFrCacheOpen(stat, &stream, frameCache), stat);
		LAL_CALL(LALDestroyFrCache(stat, &frameCache), stat);
	} else
		return(NULL);

	/* Turn on checking for missing data */
	stream->mode = LAL_FR_VERBOSE_MODE;

	/* Get the data */
	if(getcaltimeseries)
		series = get_calibrated_data(stream, chname, &start, duration, lengthlimit);
	else
		series = XLALFrReadREAL4TimeSeries(stream, chname, &start, duration, lengthlimit);


	/* Check for missing data */
	if(stream->state & LAL_FR_GAP) {
		fprintf(stderr, "get_time_series(): error: gap in data detected between GPS times %d.%09d s and %d.%09d s\n", start.gpsSeconds, start.gpsNanoSeconds, end.gpsSeconds, end.gpsNanoSeconds);
		XLALDestroyREAL4TimeSeries(series);
		series = NULL;
	}

	/* verbosity */
	if(options.verbose)
		fprintf(stderr, "get_time_series(): read %u samples (%.9lf s) at GPS time %u.%09u s\n", series->data->length, series->data->length * series->deltaT, start.gpsSeconds, start.gpsNanoSeconds);

	/* Clean up */
	XLALFrClose(stream);

	return(series);
}


/*
 * ============================================================================
 *                    Fill a time series with white noise
 * ============================================================================
 */

static void makeWhiteNoise(
	LALStatus *stat,
	REAL4TimeSeries *series,
	INT4 seed,
	REAL4 amplitude
)
{
	size_t i;
	static RandomParams *rparams = NULL;

	LAL_CALL(LALCreateRandomParams(stat, &rparams, seed), stat);
	LAL_CALL(LALNormalDeviates(stat, series->data, rparams), stat);
	LAL_CALL(LALDestroyRandomParams(stat, &rparams), stat);

	for(i = 0; i < series->data->length; i++)
		series->data->data[i] *= amplitude;
	
	if(options.verbose) {
		REAL4 norm = 0.0;
		for(i = 0; i < series->data->length; i++)
			norm += series->data->data[i] * series->data->data[i];
		fprintf(stderr, "makeWhiteNoise(): the norm is %e\n", sqrt(norm / series->data->length));
	}
}


/*
 * ============================================================================
 *                             Response function
 * ============================================================================
 */

static COMPLEX8FrequencySeries *generate_response(
	LALStatus *stat,
	const char *calcachefile,
	const char *chname,
	REAL8 deltaT,
	LIGOTimeGPS epoch,
	size_t length
)
{
	COMPLEX8FrequencySeries *response;
	size_t i;
	const LALUnit strainPerCount = {0,{0,0,0,0,0,1,-1},{0,0,0,0,0,0,0}};
	COMPLEX8 one = {1.0, 0.0};
	CalibrationUpdateParams calfacts;
	FrCache *calcache = NULL;

	memset(&calfacts, 0, sizeof(calfacts));
	calfacts.ifo = ifo;

	LAL_CALL(LALCreateCOMPLEX8FrequencySeries(stat, &response, chname, epoch, 0.0, 1.0 / (length * deltaT), strainPerCount, length / 2 + 1), stat);

	if(options.verbose) 
		fprintf(stderr, "generate_response(): working at GPS time %u.%09u s\n", response->epoch.gpsSeconds, response->epoch.gpsNanoSeconds );

	/* getting the response is handled differently for calibrated data */
	if(options.calibrated)
		for(i = 0; i < response->data->length; i++)
			response->data->data[i] = one;
	else {
		LAL_CALL(LALFrCacheImport(stat, &calcache, calcachefile), stat);
		LAL_CALL(LALExtractFrameResponse(stat, response, calcache, &calfacts), stat);
		LAL_CALL(LALDestroyFrCache(stat, &calcache), stat);
	}

	return(response);
} 


/*
 * ============================================================================
 *                              Burst injections
 * ============================================================================
 */

static void add_burst_injections(
	LALStatus *stat,
	REAL4TimeSeries *series,
	COMPLEX8FrequencySeries *response
)
{
	INT4 startTime = series->epoch.gpsSeconds;
	INT4 stopTime = startTime + series->data->length * series->deltaT;
	SimBurstTable *injections = NULL;
	INT4 calType=0;

	if(!response) {
		fprintf(stderr, "add_burst_injections(): must supply calibration information for injections\n");
		exit(1);
	}

	if(options.verbose)
		fprintf(stderr, "add_burst_injections(): reading in SimBurst Table\n");

	LAL_CALL(LALSimBurstTableFromLIGOLw(stat, &injections, burstInjectionFile, startTime, stopTime), stat);

	if(options.verbose)
		fprintf(stderr, "add_burst_injections(): injecting signals into time series\n");

	LAL_CALL(LALBurstInjectSignals(stat, series, injections, response, calType), stat); 

	if(options.verbose)
		fprintf(stderr, "add_burst_injections(): finished making the injections\n");

	while(injections) {
		SimBurstTable *thisEvent = injections;
		injections = injections->next;
		LALFree(thisEvent);
	}
}


/*
 * ============================================================================
 *                              Inspiral injections
 * ============================================================================
 */

static void add_inspiral_injections(
	LALStatus *stat,
	REAL4TimeSeries *series,
	COMPLEX8FrequencySeries *response
)
{
	INT4 startTime = series->epoch.gpsSeconds;
	INT4 stopTime = startTime + series->data->length * series->deltaT;
	SimInspiralTable *injections = NULL;

	INT4 numInjections = 0;

	if(!response) {
	  fprintf(stderr, "add_inspiral_injections(): must supply calibration information for injections\n");
	  exit(1);
	}

	if(options.verbose)
	  fprintf(stderr, "add_inspiral_injections(): reading in SimInspiral Table\n");

	numInjections = SimInspiralTableFromLIGOLw(&injections, inspiralInjectionFile, startTime, stopTime);

	if(numInjections < 0){
	  fprintf(stderr,"add_inspiral_injections():error:cannot read injection file\n");
	  exit(1);
	}

	if(options.verbose)
	  fprintf(stderr, "add_inspiral_injections(): injecting signals into time series\n");

	if(numInjections > 0)	
	  LAL_CALL(LALFindChirpInjectSignals(stat, series, injections, response), stat); 

	if(options.verbose)
	  fprintf(stderr, "add_inspiral_injections(): finished making the injections\n");
	
	while(injections) {
	  SimInspiralTable *thisEvent = injections;
	  injections = injections->next;
	  LALFree(thisEvent);
	}
}

/*
 * ============================================================================
 *                               MDC injections
 * ============================================================================
 */

static void add_mdc_injections(
	LALStatus *stat, 
	const char *mdcCacheFile, 
	REAL4TimeSeries *series,
	LIGOTimeGPS epoch,
	LIGOTimeGPS stopepoch,
	size_t lengthlimit
)
{
	REAL4TimeSeries *mdc = NULL;
	size_t i;

	if(options.verbose)
		fprintf(stderr, "add_mdc_injections(): using MDC frames for injections\n");

	{
	REAL8 old_high_pass = options.cal_high_pass;
	options.cal_high_pass = 40.0;
	mdc = get_time_series(stat, mdcDirName, mdcCacheFile, options.mdcchannelName, epoch, stopepoch, lengthlimit, TRUE);
	options.cal_high_pass = old_high_pass;
	}

	/* write diagnostic info to disk */
	if(options.printData)
		LALPrintTimeSeries(mdc, "./timeseriesmdc.dat");

	/* add the mdc signal to the given time series */
	for(i = 0; i < series->data->length; i++)
		series->data->data[i] += mdc->data->data[i];

	/* clean up */
	LAL_CALL(LALDestroyREAL4TimeSeries(stat, mdc), stat);
}

/*
 * ============================================================================
 *                           Sim injections
 * ============================================================================
 */
static void add_sim_injections(
	LALStatus *stat,
	REAL4TimeSeries *series,
	COMPLEX8FrequencySeries *response,
	size_t lengthlimit
	)
{
	REAL4TimeSeries   *signal;
	DetectorResponse   detector;
	LALDetector       *tmpDetector = NULL;
	CoherentGW         waveform;
	REAL4             *aData;
	LALTimeInterval   epochCorrection;
	COMPLEX8FrequencySeries  *transfer = NULL;
	COMPLEX8Vector    *unity = NULL;

	char pluschan[30]; 
	char crosschan[30];
	LIGOTimeGPS start,end;
	FILE *fp = NULL;
	UINT4 i, n;
	REAL8 simDuration;

	INT4 startTime = series->epoch.gpsSeconds;
	INT4 stopTime = startTime + series->data->length * series->deltaT;
	SimBurstTable *injections = NULL;
	SimBurstTable *simBurst = NULL;
	BurstParamStruc    burstParam;
	
	if(!response) {
	  fprintf(stderr, "add_sim_injections(): must supply calibration information for injections\n");
	  exit(1);
	}
	
	/* allocate memory */
	memset( &detector, 0, sizeof( DetectorResponse ) );
	transfer = (COMPLEX8FrequencySeries *)LALCalloc( 1, sizeof(COMPLEX8FrequencySeries) );
	if (!transfer ){
	  fprintf(stderr, "add_sim_injections(): detector.transfer not allocated\n");
	  exit(1);
	}

	memcpy( &(transfer->epoch), &(response->epoch), sizeof(LIGOTimeGPS) );
	transfer->f0 = response->f0;
	transfer->deltaF = response->deltaF;
	
	tmpDetector = detector.site = (LALDetector *) LALMalloc( sizeof(LALDetector) );

	/* set the detector site */
	switch ( series->name[0] )
	  {
	  case 'H':
	    *(detector.site) = lalCachedDetectors[LALDetectorIndexLHODIFF];
	    LALWarning( stat, "computing waveform for Hanford." );
	    break;
	  case 'L':
	    *(detector.site) = lalCachedDetectors[LALDetectorIndexLLODIFF];
	    LALWarning( stat, "computing waveform for Livingston." );
	    break;
	  default:
	    LALFree( detector.site );
	    detector.site = NULL;
	    tmpDetector = NULL;
	    LALWarning( stat, "Unknown detector site, computing plus mode "
			"waveform with no time delay" );
	    break;
	  }

	/* set up units for the transfer function */
	{
	  RAT4 negOne = { -1, 0 };
	  LALUnit unit;
	  LALUnitPair pair;
	  pair.unitOne = &lalADCCountUnit;
	  pair.unitTwo = &lalStrainUnit;
	  LAL_CALL(LALUnitRaise( stat, &unit, pair.unitTwo, &negOne ),stat);
	  pair.unitTwo = &unit;
	  LAL_CALL(LALUnitMultiply( stat, &(transfer->sampleUnits), &pair ),stat);
	}	

	/* invert the response function to get the transfer function */
	LAL_CALL(LALCCreateVector( stat, &( transfer->data ), response->data->length ),stat);
	
	LAL_CALL(LALCCreateVector( stat, &unity, response->data->length ),stat);
	for ( i = 0; i < response->data->length; ++i ) 
	  {
	    unity->data[i].re = 1.0;
	    unity->data[i].im = 0.0;
	  }
	
	LAL_CALL(LALCCVectorDivide( stat, transfer->data, unity, response->data ),stat);
	
	LAL_CALL(LALCDestroyVector( stat, &unity ),stat);
       	
	/* Set up a time series to hold signal in ADC counts */
	LAL_CALL(LALCreateREAL4TimeSeries(stat, &signal, series->name, series->epoch, series->f0, series->deltaT, series->sampleUnits, series->data->length), stat);
 	
	if(options.verbose)
	  fprintf(stderr, "add_sim_injections(): reading in SimBurst Table\n");
	
	LAL_CALL(LALSimBurstTableFromLIGOLw(stat, &injections, simInjectionFile, startTime, stopTime), stat);
	
	simBurst = injections;
	while ( simBurst ){
	  REAL4TimeSeries    *plusseries = NULL;
	  REAL4TimeSeries    *crossseries = NULL;
	  
	  /* set the burst params */
	  burstParam.deltaT = series->deltaT;
	  if( !( strcmp( simBurst->coordinates, "HORIZON" ) ) ){
	    burstParam.system = COORDINATESYSTEM_HORIZON;
	  }
	  else if ( !( strcmp( simBurst->coordinates, "ZENITH" ) ) ){
	    /* set coordinate system for completeness */
	    burstParam.system = COORDINATESYSTEM_EQUATORIAL;
	    detector.site = NULL;
	  }
	  else if ( !( strcmp( simBurst->coordinates, "GEOGRAPHIC" ) ) ){
	    burstParam.system = COORDINATESYSTEM_GEOGRAPHIC;
	  }
	  else if ( !( strcmp( simBurst->coordinates, "EQUATORIAL" ) ) ){
	    burstParam.system = COORDINATESYSTEM_EQUATORIAL;
	  }
	  else if ( !( strcmp( simBurst->coordinates, "ECLIPTIC" ) ) ){
	    burstParam.system = COORDINATESYSTEM_ECLIPTIC;
	  }
	  else if ( !( strcmp( simBurst->coordinates, "GALACTIC" ) ) ){
	    burstParam.system = COORDINATESYSTEM_GALACTIC;
	  }
	  else
	    burstParam.system = COORDINATESYSTEM_EQUATORIAL;

	  /* Set the channel names */
	  snprintf(pluschan,30,"SIM_plus_%s_%d",simBurst->waveform,simBurst->zm_number);	  
	  snprintf(crosschan,30,"SIM_cross_%s_%d",simBurst->waveform,simBurst->zm_number );
	  
	  /*Set the start and end times of the sim waveforms */
	  start.gpsSeconds = 0;
	  start.gpsNanoSeconds = 0;
	  simDuration = (simBurst->dtplus + simBurst->dtminus);
	  n = (INT4) (simDuration / series->deltaT);
 
	  LAL_CALL(LALAddFloatToGPS(stat, &end, &start, simDuration), stat);

	  options.getcaltimeseries = FALSE;
	  /* Get the plus time series */
	  plusseries = get_time_series(stat, options.simdirname, options.simCacheFile, pluschan, start, end, lengthlimit, options.getcaltimeseries );
	  
	  /* Get the cross time series */
	  crossseries = get_time_series(stat, options.simdirname, options.simCacheFile, crosschan, start, end, lengthlimit, options.getcaltimeseries  );
	  
	  /* write diagnostic info to disk */
	  if(options.printData){
	    fp = fopen("timeseriessim.dat","w");
	    for(i = 0;(unsigned)i<plusseries->data->length;i++){
	      fprintf(fp,"%e %e %e\n",i*plusseries->deltaT,plusseries->data->data[i],crossseries->data->data[i] );
	    }
	    fclose(fp);
	  }

	  /* read in the waveform in a CoherentGW struct */
	  memset( &waveform, 0, sizeof(CoherentGW) );
	  /* this step sets the adata,fdata,phidata to 0 */ 
	  LAL_CALL(LALGenerateBurst( stat, &waveform, simBurst, &burstParam ),stat);

	  /* correct the waveform epoch:
	   * Remember the peak time is always at the center of the frame
	   * Hence the epoch of the waveform is set at 1/2 of the duration
	   * before the geocent_peak_time, since geocent_peak_time should match
	   * with the peak time in the frame
	   */
	  simDuration = simDuration / 2.0;
	  LAL_CALL(LALFloatToInterval(stat, &epochCorrection, &simDuration), stat);
	  LAL_CALL(LALDecrementGPS( stat, &(waveform.a->epoch), &(simBurst->geocent_peak_time), &epochCorrection), stat);
	
	  aData = waveform.a->data->data;

	  /* copy the plus and cross data properly scaled for distance.
	   *
	   * NOTE: The waveforms in the frames are always produced at
	   * distance of 10000 Kpc. However the distance in the 
	   * parameter file, BBHWaveGen.in, shu'd be 1000 Kpc., since 
	   * in the wave generation script the definition of kpc 
	   * is 1 kpc = 3.086e20 m . where as the right definition is
	   * 1 kpc = 3.086e19 m.
	   */
	  for( i = 0; i < n; i++){
	    *(aData++) = plusseries->data->data[i] * 10000/ simBurst->distance;
	    *(aData++) = crossseries->data->data[i] * 10000/ simBurst->distance;
	  }

	  /* must set the epoch of signal since it's used by coherent GW */
	  signal->epoch = waveform.a->epoch;
	  detector.transfer=NULL;

	  /* convert this into an ADC signal */
	  LAL_CALL(LALSimulateCoherentGW( stat, signal, &waveform, &detector ),stat);	
	  XLALRespFilt(signal, transfer);

	  if(options.printData){
	    fp = fopen("timeseriessignal.dat","w");
	    for(i = 0;(unsigned)i<signal->data->length;i++){
	      fprintf(fp,"%e %e\n",i*signal->deltaT,signal->data->data[i] );
	    }
	    fclose(fp);
	  }
	  /* inject the signal into the data channel */
	  LAL_CALL(LALSSInjectTimeSeries( stat, series, signal ),stat);

	  /* Clean up */
	  LAL_CALL(LALDestroyREAL4TimeSeries(stat, plusseries), stat);
	  LAL_CALL(LALDestroyREAL4TimeSeries(stat, crossseries), stat);
	  LAL_CALL(LALSDestroyVectorSequence(stat, &( waveform.a->data ) ),stat);
	  LAL_CALL(LALSDestroyVector(stat, &( waveform.f->data ) ),stat);
	  LAL_CALL(LALDDestroyVector(stat, &( waveform.phi->data ) ),stat);
	  LALFree( waveform.a );   waveform.a = NULL;
	  LALFree( waveform.f );   waveform.f = NULL;
	  LALFree( waveform.phi );  waveform.phi = NULL;
	  
	  /* reset the detector site information in case it changed */
	  detector.site = tmpDetector;
	  
	  /* move on to next one */
	  simBurst = simBurst->next;
	}
	
	LAL_CALL(LALDestroyREAL4TimeSeries(stat, signal), stat);
	LAL_CALL(LALCDestroyVector(stat, &( transfer->data )),stat);
	
	if ( detector.site ) 
	  LALFree( detector.site );
	LALFree( transfer );
	
	while(injections) {
	  SimBurstTable *thisEvent = injections;
	  injections = injections->next;
	  LALFree(thisEvent);
	}
	
}


/*
 * ============================================================================
 *                               Analysis Loop
 * ============================================================================
 */

/*
 * Analyze a time series in intervals corresponding to the length of time
 * for which the instrument's noise is stationary.
 */

static SnglBurstTable **analyze_series(
	SnglBurstTable **addpoint,
	COMPLEX8FrequencySeries  *hrssresponse,
	REAL4TimeSeries *series,
	size_t psdlength,
	EPSearchParams *params
)
{
	REAL4TimeSeries *interval;
	size_t i, start;
	size_t overlap = options.windowLength - params->windowShift;

	if(psdlength > series->data->length) {
		psdlength = block_commensurate(series->data->length, options.windowLength, params->windowShift);
		if(options.verbose)
			fprintf(stderr, "analyze_series(): warning: PSD average length exceeds available data --- reducing PSD average length to %zu samples\n", psdlength);
		if(!psdlength) {
			if(options.verbose)
				fprintf(stderr, "analyze_series(): warning: cowardly refusing to analyze 0 samples... skipping series\n");
			return(addpoint);
		}
	}

	for(i = 0; i < series->data->length - overlap; i += psdlength - overlap) {
		start = min(i, series->data->length - psdlength);
		
		interval = XLALCutREAL4TimeSeries(series, start, psdlength);

		if(options.verbose)
			fprintf(stderr, "analyze_series(): analyzing samples %zu -- %zu (%.9lf s -- %.9lf s)\n", start, start + interval->data->length, start * interval->deltaT, (start + interval->data->length) * interval->deltaT);

		*addpoint = XLALEPSearch(hrssresponse, interval, params);
		while(*addpoint)
			addpoint = &(*addpoint)->next;
		if(xlalErrno) {
			fprintf(stderr, "analyze_series(): fatal error: XLALEPSearch() returned failure\n");
			exit(1);
		}

		XLALDestroyREAL4TimeSeries(interval);
	}

	return(addpoint);
}


/*
 * ============================================================================
 *                                   Output
 * ============================================================================
 */

static void output_results(LALStatus *stat, char *file, MetadataTable *procTable, MetadataTable *procparams, MetadataTable *searchsumm, SnglBurstTable *burstEvent)
{
	LIGOLwXMLStream xml;
	LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
	MetadataTable myTable;

	memset(&xml, 0, sizeof(LIGOLwXMLStream));
	LAL_CALL(LALOpenLIGOLwXMLFile(stat, &xml, file), stat);

	/* process table */
	snprintf(procTable->processTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
	LAL_CALL(LALGPSTimeNow(stat, &(procTable->processTable->end_time), &accuracy), stat);
	LAL_CALL(LALBeginLIGOLwXMLTable(stat, &xml, process_table), stat);
	LAL_CALL(LALWriteLIGOLwXMLTable(stat, &xml, *procTable, process_table), stat);
	LAL_CALL(LALEndLIGOLwXMLTable(stat, &xml), stat);

	/* process params table */
	LAL_CALL(LALBeginLIGOLwXMLTable(stat, &xml, process_params_table), stat);
	LAL_CALL(LALWriteLIGOLwXMLTable(stat, &xml, *procparams, process_params_table), stat);
	LAL_CALL(LALEndLIGOLwXMLTable(stat, &xml), stat);

	/* search summary table */
	snprintf(searchsumm->searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, "%s", ifo);
	searchsumm->searchSummaryTable->nevents = XLALCountSnglBurst(burstEvent);
	LAL_CALL(LALBeginLIGOLwXMLTable(stat, &xml, search_summary_table), stat);
	LAL_CALL(LALWriteLIGOLwXMLTable(stat, &xml, *searchsumm, search_summary_table), stat);
	LAL_CALL(LALEndLIGOLwXMLTable(stat, &xml), stat);

	/* burst table */
	LAL_CALL(LALBeginLIGOLwXMLTable(stat, &xml, sngl_burst_table), stat);
	myTable.snglBurstTable = burstEvent;
	LAL_CALL(LALWriteLIGOLwXMLTable(stat, &xml, myTable, sngl_burst_table), stat);
	LAL_CALL(LALEndLIGOLwXMLTable(stat, &xml), stat);

	LAL_CALL(LALCloseLIGOLwXMLFile(stat, &xml), stat);
}


/*
 * ============================================================================
 *                                Entry point
 * ============================================================================
 */

int main( int argc, char *argv[])
{
	LALStatus                 stat;
	LALLeapSecAccuracy        accuracy = LALLEAPSEC_LOOSE;
	EPSearchParams            params;
	LIGOTimeGPS               epoch;
	LIGOTimeGPS               boundepoch;
	size_t                    overlap;
	CHAR                      outfilename[256];
	REAL4TimeSeries          *series = NULL;
	COMPLEX8FrequencySeries  *hrssresponse = NULL;
	SnglBurstTable           *burstEvent = NULL;
	SnglBurstTable          **EventAddPoint = &burstEvent;
	MetadataTable             procTable;
	MetadataTable             procparams;
	MetadataTable             searchsumm;

	/*
	 * Initialize everything
	 */

	memset(&stat, 0, sizeof(stat));
	lal_errhandler = LAL_ERR_EXIT;
	set_debug_level("3");
	parse_command_line_debug(argc, argv);

	/* create the process and process params tables */
	procTable.processTable = LALCalloc(1, sizeof(ProcessTable));
	LAL_CALL(LALGPSTimeNow(&stat, &(procTable.processTable->start_time), &accuracy), &stat);
	LAL_CALL(populate_process_table(&stat, procTable.processTable, PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE), &stat);
	procparams.processParamsTable = NULL;

	/* parse arguments and fill procparams table */
	parse_command_line(argc, argv, &params, &procparams);

	/* create the search summary table */
	searchsumm.searchSummaryTable = LALCalloc(1, sizeof(SearchSummaryTable));

	/* fill the comment */
	snprintf(procTable.processTable->comment, LIGOMETA_COMMENT_MAX, "%s", options.comment);
	snprintf(searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, "%s", options.comment);

	/* the number of nodes for a standalone job is always 1 */
	searchsumm.searchSummaryTable->nnodes = 1;

	/* store the input start and end times */
	searchsumm.searchSummaryTable->in_start_time = options.startEpoch;
	searchsumm.searchSummaryTable->in_end_time = options.stopEpoch;

	/* determine the input time series post-conditioning overlap, and set
	 * the outer loop's upper bound */
	overlap = options.windowLength - params.windowShift;
	if(options.verbose)
		fprintf(stderr, "%s: time series overlap is %zu samples (%.9lf s)\n", argv[0], overlap, overlap / (double) options.ResampleRate);

	LAL_CALL(LALAddFloatToGPS(&stat, &boundepoch, &options.stopEpoch, -(2 * options.FilterCorruption + (int) overlap) / (double) options.ResampleRate), &stat);
	if(options.verbose)
		fprintf(stderr, "%s: time series epochs must be less than %u.%09u s\n", argv[0], boundepoch.gpsSeconds, boundepoch.gpsNanoSeconds);


	/*
	 * ====================================================================
	 *                         Outer Loop
	 * ====================================================================
	 *
	 * Split the total length of time to be analyzed into time series
	 * small enough to fit in RAM.
	 */

	for(epoch = options.startEpoch; CompareGPS(&stat, &epoch, &boundepoch) < 0;) {
		/*
		 * Get the data, if reading calibrated data set the flag
		 */
		if (options.calibrated)
			options.getcaltimeseries = TRUE;

		series = get_time_series(&stat, dirname, cachefile, options.channelName, epoch, options.stopEpoch, options.maxSeriesLength, options.getcaltimeseries);

		/*
		 * If we specified input files but nothing got read, there
		 * was an error.
		 */

		if((dirname || cachefile) && !series) {
			fprintf(stderr, "%s: error: failure reading input data\n", argv[0]);
			exit(1);
		}

		/*
		 * Create an empty series of 1s if we didn't read one (eg.
		 * if no input files were specified), then add white noise
		 * to the series.
		 */

		if(options.noiseAmpl > 0.0) {
			if(!series) {
				size_t i, length;
				length = DeltaGPStoFloat(&stat, &options.stopEpoch, &epoch) * options.ResampleRate;
				if(options.maxSeriesLength)
					length = min(options.maxSeriesLength, length);
				LAL_CALL(LALCreateREAL4TimeSeries(&stat, &series, options.channelName, epoch, 0.0, (REAL8) 1.0 / options.ResampleRate, lalADCCountUnit, length), &stat);
				for(i = 0; i < series->data->length; i++)
					series->data->data[i] = 1.0;
			}
			makeWhiteNoise(&stat, series, options.seed, options.noiseAmpl);
		}

		/*
		 * Write diagnostic info to disk
		 */

		if(options.printData)
			LALPrintTimeSeries(series, "./timeseriesasq.dat");

		/*
		 * Add burst/inspiral injections into the time series if requested.
		 */

		if(burstInjectionFile || inspiralInjectionFile || options.simCacheFile) {
			COMPLEX8FrequencySeries  *response = NULL;

			/* Create the response function */
			if(options.calCacheFile)
				response = generate_response(&stat, options.calCacheFile, options.channelName, series->deltaT, series->epoch, series->data->length);


			if(options.printData)
			  LALCPrintFrequencySeries(response, "./response.dat");

			/* perform injections */
			if(burstInjectionFile)
			  add_burst_injections(&stat, series, response);
			else if(inspiralInjectionFile)
			  add_inspiral_injections(&stat, series, response);
			else if(options.simCacheFile)
			  add_sim_injections(&stat, series, response, options.maxSeriesLength);
		
			if(options.printData)
				LALPrintTimeSeries(series, "./injections.dat");

			/*clean up*/
			LAL_CALL(LALDestroyCOMPLEX8FrequencySeries(&stat, response), &stat);
		}

		/*
		 * Add MDC injections into the time series if requested.
		 */

		if(mdcCacheFile) {
		        add_mdc_injections(&stat, mdcCacheFile, series, epoch, options.stopEpoch, options.maxSeriesLength);
			if(options.printData)
				LALPrintTimeSeries(series, "./timeseriesasqmdc.dat");
		}

		/*
		 * Generate the response function at the right deltaf to be usd for h_rss estimation.
		 */

		if (options.estimateHrss) {
			hrssresponse = generate_response(&stat, options.calCacheFile, options.channelName, series->deltaT, series->epoch, options.windowLength);
			if(options.printData)
			  LALCPrintFrequencySeries(hrssresponse, "./hrssresponse.dat");
		}

		/*
		 * Condition the time series data.
		 */

		/* Scale the time series if calibrated data */
		if (options.calibrated){
		  UINT4 i;
		  for(i = 0;i<series->data->length;i++){
		    series->data->data[i] = SCALE*series->data->data[i];
		  }
		}

		if(XLALEPConditionData(series, options.high_pass, (REAL8) 1.0 / options.ResampleRate, options.FilterCorruption)) {
			fprintf(stderr, "%s: XLALEPConditionData() failed.\n", argv[0]);
			exit(1);
		}

		if(options.printData)
		  LALPrintTimeSeries(series, "./condtimeseries.dat");

		if(options.verbose)
			fprintf(stderr, "%s: %u samples (%.9f s) at GPS time %d.%09d s remain after conditioning\n", argv[0], series->data->length, series->data->length * series->deltaT, series->epoch.gpsSeconds, series->epoch.gpsNanoSeconds);

		/*
		 * Store the start and end times of the data that actually
		 * gets analyzed.
		 */

		if(!searchsumm.searchSummaryTable->out_start_time.gpsSeconds)
			LAL_CALL(LALAddFloatToGPS(&stat, &searchsumm.searchSummaryTable->out_start_time, &series->epoch, series->deltaT * (options.windowLength / 2 - params.windowShift)), &stat);
		LAL_CALL(LALAddFloatToGPS(&stat, &searchsumm.searchSummaryTable->out_end_time, &series->epoch, series->deltaT * (series->data->length - (options.windowLength / 2 - params.windowShift))), &stat);


		/*
		 * Analyze the data
		 */

		EventAddPoint = analyze_series(EventAddPoint, hrssresponse, series, options.PSDAverageLength, &params);

		/*
		 * Reset for next run
		 *
		 * Advancing the epoch by the post-conditioning series length
		 * provides exactly the overlap needed to account for
		 * conditioning corruption.  The post-conditioning overlap
		 * computed earlier provides the additional overlap that is
		 * needed.
		 */

		LAL_CALL(LALAddFloatToGPS(&stat, &epoch, &epoch, (series->data->length - overlap) * series->deltaT), &stat);
		LAL_CALL(LALDestroyREAL4TimeSeries(&stat, series), &stat);
		LAL_CALL(LALDestroyCOMPLEX8FrequencySeries(&stat, hrssresponse), &stat);
	}

	/*
	 * Cluster and sort the events.
	 */

	if(options.cluster)
		XLALClusterSnglBurstTable(&burstEvent, XLALCompareSnglBurstByPeakTime, XLALCompareSnglBurstByPeakTimeAndFreq, XLALSnglBurstCluster);
	XLALSortSnglBurst(&burstEvent, XLALCompareSnglBurstByStartTimeAndLowFreq);

	/*
	 * Check event rate limit.
	 */

	if((options.max_event_rate > 0) && (SnglBurstTableLength(burstEvent) > XLALDeltaFloatGPS(&searchsumm.searchSummaryTable->out_end_time, &searchsumm.searchSummaryTable->out_start_time) * options.max_event_rate)) {
		XLALPrintError("%s: event rate limit exceeded!", argv[0]);
		exit(1);
	}

	/*
	 * Output the results.
	 */

	snprintf(outfilename, sizeof(outfilename)-1, "%s-POWER_%s-%d-%d.xml", ifo, options.comment, options.startEpoch.gpsSeconds, options.stopEpoch.gpsSeconds - options.startEpoch.gpsSeconds);
	outfilename[sizeof(outfilename)-1] = '\0';
	output_results(&stat, outfilename, &procTable, &procparams, &searchsumm, burstEvent);

	/*
	 * Final cleanup.
	 */

	XLALDestroyREAL4Window(params.window);
	LALFree(procTable.processTable);
	LALFree(searchsumm.searchSummaryTable);

	while(procparams.processParamsTable) {
		ProcessParamsTable *table = procparams.processParamsTable;
		procparams.processParamsTable = table->next;
		LALFree(table);
	}

	while(burstEvent) {
		SnglBurstTable *event = burstEvent;
		burstEvent = burstEvent->next;
		LALFree(event);
	}

	LALCheckMemoryLeaks();
	exit(0);
}

#endif
