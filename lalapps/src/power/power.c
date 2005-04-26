#include <getopt.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lal/AVFactories.h>
#include <lal/Date.h>
#include <lal/EPSearch.h>
#include <lal/ExcessPower.h>
#include <lal/FrameCalibration.h>
#include <lal/FrameStream.h>
#include <lal/FrequencySeries.h>
#include <lal/GenerateBurst.h>
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
#include <lal/ResampleTimeSeries.h>
#include <lal/TFTransform.h>
#include <lal/TimeFreqFFT.h>
#include <lal/TimeSeries.h>
#include <lal/Units.h>

#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpTD.h>
#include <lal/FindChirpBCV.h>
#include <lal/FindChirpBCVSpin.h>
#include <lal/FindChirpChisq.h>

#include <lalapps.h>
#include <processtable.h>

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

/*
 * ============================================================================
 *                                Global Data
 * ============================================================================
 */

/* Parameters from command line */
static struct {
	CHAR *calCacheFile;         /* name of the calibration cache file  */
	int cluster;                /* TRUE == perform clustering          */
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
	REAL8 high_pass;            /* conditioning high pass freq (Hz)    */
	REAL8 cal_high_pass;        /* double->single high pass freq (Hz)  */
} options;
 
/* global variables */
FrChanIn mdcchannelIn;              /* mdc signal only channnel info       */
EPSearchParams *mdcparams;          /* mdc search param                    */
CHAR ifo[3];                        /* two character interferometer        */
CHAR *cachefile;                    /* name of file with frame cache info  */
CHAR *dirname;                      /* name of directory with frames       */

/* data conditioning parameters */
CHAR *burstInjectionFile;           /* file with list of burst injections  */
CHAR *inspiralInjectionFile;        /* file with list of burst injections  */
CHAR *mdcCacheFile;                 /* name of mdc signal cache file       */
ResampleTSFilter resampFiltType;


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
"	[--calibration-cache <cache file>]\n" \
"	 --channel-name <string>\n" \
"	[--cluster]\n" \
"	[--debug-level <level>]\n" \
"	 --default-alpha <alpha>\n" \
"	[--event-limit <count>]\n" \
"	 --filter-corruption <samples>\n" \
"	 --frame-cache <cache file>\n" \
"	 --frame-dir <directory>\n" \
"	[--calibrated-data <high pass frequency>]\n" \
"	 --gps-end-time <seconds>\n" \
"	 --gps-end-time-ns <nanoseconds>\n" \
"	 --gps-start-time <seconds>\n" \
"	 --gps-start-time-ns <nanoseconds>\n" \
"	[--help]\n" 
"	[--burstinjection-file <file name>]\n" \
"	[--inspiralinjection-file <file name>]\n" \
"	 --low-freq-cutoff <Hz>\n" \
"	[--mdc-cache <cache file>]\n" \
"	[--mdc-channel <channel name>]\n" \
"	 --max-tileband <Hz>\n" \
"	[--max-tileduration <samples>]\n" \
"	 --min-freq-bin <nfbin>\n" \
"	 --min-time-bin <ntbin>\n" \
"	[--noise-amplitude <amplitude>]\n" \
"	 --nsigma <sigma>\n" \
"	[--printData]\n" \
"	[--printSpectrum]\n" \
"	 --psd-average-method <method>\n" \
"	 --psd-average-points <samples>\n" \
"	[--ram-limit <MebiBytes>]\n" \
"	 --resample-filter <filter type>\n" \
"	 --resample-rate <Hz>\n" \
"	[--seed <seed>]\n" \
"	 --tile-overlap-factor <factor>\n" \
"	 --threshold <threshold>\n" \
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

static int check_for_missing_parameters(LALStatus *stat, char *prog, struct option *long_options, EPSearchParams *params)
{
	int index;
	int got_all_arguments = TRUE;
	int arg_is_missing;
	INT8 gpstmp;

	for(index = 0; long_options[index].name; index++) {
		switch(long_options[index].val) {
			case 'A':
			arg_is_missing = !params->tfTilingInput.length;
			break;

			case 'C':
			arg_is_missing = !params->channelName;
			break;

			case 'E':
			arg_is_missing = params->compEPInput.alphaDefault > 1.0;
			break;

			case 'K':
			LAL_CALL(LALGPStoINT8(stat, &gpstmp, &options.stopEpoch), stat);
			arg_is_missing = !gpstmp;
			break;

			case 'M':
			LAL_CALL(LALGPStoINT8(stat, &gpstmp, &options.startEpoch), stat);
			arg_is_missing = !gpstmp;
			break;

			case 'Q':
			arg_is_missing = params->tfTilingInput.flow < 0.0;
			break;

			case 'T':
			arg_is_missing = params->tfTilingInput.minFreqBins <= 0;
			break;

			case 'U':
			arg_is_missing = params->tfTilingInput.minTimeBins <= 0;
			break;

			case 'W':
			arg_is_missing = !params->windowLength;
			break;

			case 'X':
			arg_is_missing = params->compEPInput.numSigmaMin <= 0.0;
			break;

			case 'Y':
			arg_is_missing = params->method == (unsigned) -1;
			break;

			case 'Z':
			arg_is_missing = !options.PSDAverageLength;
			break;

			case 'b':
			arg_is_missing = resampFiltType == (unsigned) -1;
			break;

			case 'd':
			arg_is_missing = !params->windowShift;
			break;

			case 'e':
			arg_is_missing = !options.ResampleRate;
			break;

			case 'f':
			arg_is_missing = !params->tfTilingInput.overlapFactor;
			break;

			case 'g':
			arg_is_missing = params->alphaThreshold < 0.0;
			break;

			case 'i':
			arg_is_missing = params->windowType == NumberWindowTypes;
			break;

			case 'j':
			arg_is_missing = options.FilterCorruption < 0;
			break;

			case 'l':
			arg_is_missing = !params->tfTilingInput.maxTileBand;
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
	int ram = 0;	/* default */
	INT8 gpstmp;
	INT8 gpsStartTimeNS = 0;	/* impossible */
	INT8 gpsStopTimeNS = 0;	/* impossible */
	ProcessParamsTable **paramaddpoint = &procparams->processParamsTable;
	LALStatus stat = blank_status;
	struct option long_options[] = {
		{"bandwidth",           required_argument, NULL,           'A'},
		{"calibrated-data",	required_argument, NULL,           'J'},
		{"calibration-cache",   required_argument, NULL,           'B'},
		{"channel-name",        required_argument, NULL,           'C'},
		{"cluster",             no_argument, &options.cluster,    TRUE},
		{"debug-level",         required_argument, NULL,           'D'},
		{"default-alpha",       required_argument, NULL,           'E'},
		{"event-limit",         required_argument, NULL,           'F'},
		{"filter-corruption",	required_argument, NULL,           'j'},
		{"frame-cache",         required_argument, NULL,           'G'},
		{"frame-dir",           required_argument, NULL,           'H'},
		{"gps-end-time",        required_argument, NULL,           'K'},
		{"gps-end-time-ns",     required_argument, NULL,           'L'},
		{"gps-start-time",      required_argument, NULL,           'M'},
		{"gps-start-time-ns",   required_argument, NULL,           'N'},
		{"help",                no_argument,       NULL,           'O'},
		{"high-pass",           required_argument, NULL,           'o'},
		{"burstinjection-file", required_argument, NULL,           'P'},
		{"inspiralinjection-file", required_argument, NULL,        'I'},
		{"low-freq-cutoff",     required_argument, NULL,           'Q'},
		{"max-tileband",        required_argument, NULL,           'l'},
		{"max-tileduration",    required_argument, NULL,           'm'},
		{"mdc-cache",           required_argument, NULL,           'R'},
		{"mdc-channel",         required_argument, NULL,           'S'},
		{"min-freq-bin",        required_argument, NULL,           'T'},
		{"min-time-bin",        required_argument, NULL,           'U'},
		{"noise-amplitude",     required_argument, NULL,           'V'},
		{"nsigma",              required_argument, NULL,           'X'},
		{"printData",           no_argument, &options.printData,  TRUE},
		{"printSpectrum",       no_argument, &printSpectrum,      TRUE},
		{"psd-average-method",  required_argument, NULL,           'Y'},
		{"psd-average-points",  required_argument, NULL,           'Z'},
		{"ram-limit",           required_argument, NULL,           'a'},
		{"resample-filter",     required_argument, NULL,           'b'},
		{"resample-rate",       required_argument, NULL,           'e'},
		{"seed",                required_argument, NULL,           'c'},
		{"tile-overlap-factor", required_argument, NULL,           'f'},
		{"threshold",           required_argument, NULL,           'g'},
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

	params->alphaThreshold = -1.0;	/* impossible */
	params->channelName = NULL;	/* impossible */
	params->eventLimit = 999;	/* default */
	params->method = -1;	/* impossible */
	params->compEPInput.alphaDefault = 2.0;	/* impossible */
	params->compEPInput.numSigmaMin = -1.0;	/* impossible */
	params->tfPlaneMethod = useSingleTFPlane ; /* The only method possible */
	params->tfTilingInput.flow = -1.0;	/* impossible */
	params->tfPlaneParams.fhigh = -1.0;	/* impossible */
	params->tfTilingInput.length = 0;	/* impossible */
	params->tfTilingInput.minFreqBins = 0;	/* impossible */
	params->tfTilingInput.minTimeBins = 0;	/* impossible */
	params->tfTilingInput.maxTileBand = 0;  /* impossible */
	/*params->tfTilingInput.maxTileBand = 32.0;	default */
	params->tfTilingInput.overlapFactor = 0;	/* impossible */
	params->windowType = NumberWindowTypes;	/* impossible */
	params->windowLength = 0;	/* impossible */
	params->windowShift = 0;	/* impossible */

	options.calCacheFile = NULL;	/* default */
	options.cluster = FALSE;	/* default */
	options.comment = "";		/* default */
	options.FilterCorruption = -1;	/* impossible */
	options.noiseAmpl = -1.0;	/* default */
	options.printData = FALSE;	/* default */
	options.PSDAverageLength = 0;	/* impossible */
	options.ResampleRate = 0;	/* impossible */
	options.seed = 1;	/* default */
	options.verbose = FALSE;	/* default */
	options.calibrated = FALSE;	/* default */
	options.high_pass = -1.0;	/* impossible */
	options.cal_high_pass = -1.0;	/* impossible */

	cachefile = NULL;	/* default */
	dirname = NULL;	/* default */
	memset(ifo, 0, sizeof(ifo));	/* default */
	burstInjectionFile = NULL;	/* default */
	inspiralInjectionFile = NULL;	/* default */
	mdcCacheFile = NULL;	/* default */
	mdcparams = NULL;	/* default */
	printSpectrum = FALSE;	/* default */
	useoverwhitening = FALSE; /* default */
	resampFiltType = -1;	/* default */

	/*
	 * Parse command line.
	 */

	opterr = 1;	/* enable error messages */
	optind = 0;	/* start scanning from argv[0] */
	do switch(c = getopt_long(argc, argv, "", long_options, &option_index)) {
		case 'A':
		params->tfTilingInput.length = atoi(optarg);
		if(params->tfTilingInput.length <= 0 || !is_power_of_2(params->tfTilingInput.length) ) {
			sprintf(msg, "must be > 0 and a power of 2(%i specified)", params->tfTilingInput.length);
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
		params->channelName = optarg;
		memcpy(ifo, optarg, sizeof(ifo) - 1);
		ADD_PROCESS_PARAM("string");
		break;

		case 'D':
		/* only add --debug-level to params table in this pass */
		ADD_PROCESS_PARAM("string");
		break;

		case 'E':
		params->compEPInput.alphaDefault = atof(optarg);
		if(params->compEPInput.alphaDefault < 0.0 || params->compEPInput.alphaDefault > 1.0) {
			sprintf(msg, "must be in range [0,1] (%f specified)", params->compEPInput.alphaDefault);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

		case 'F':
		params->eventLimit = atoi(optarg);
		if(params->eventLimit < 1 || params->eventLimit > 999) {
			sprintf(msg, "must be in range [1,999] (%i specified)", params->eventLimit);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
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
		gpstmp = atol(optarg);
		if(gpstmp < 441417609 || gpstmp > 999999999) {
			sprintf(msg, "must be in the range [Jan 01 1994 00:00:00 UTC, Sep 14 2011 01:46:26 UTC] (%lld specified)", gpstmp);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		gpsStopTimeNS += gpstmp * 1000000000LL;
		ADD_PROCESS_PARAM("int");
		break;

		case 'L':
		gpstmp = atol(optarg);
		if(gpstmp < 0 || gpstmp > 999999999) {
			sprintf(msg, "must be in the range [0,999999999] (%lld specified)", gpstmp);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		gpsStopTimeNS += gpstmp;
		ADD_PROCESS_PARAM("int");
		break;

		case 'M':
		gpstmp = atol(optarg);
		if(gpstmp < 441417609 || gpstmp > 999999999) {
			sprintf(msg, "must be in the range [Jan 01 1994 00:00:00 UTC, Sep 14 2011 01:46:26 UTC] (%lld specified)", gpstmp);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		gpsStartTimeNS += gpstmp * 1000000000LL;
		ADD_PROCESS_PARAM("int");
		break;

		case 'N':
		gpstmp = atol(optarg);
		if(gpstmp < 0 || gpstmp > 999999999) {
			sprintf(msg, "must be in the range [0,999999999] (%lld specified)", gpstmp);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		gpsStartTimeNS += gpstmp;
		ADD_PROCESS_PARAM("int");
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
		params->tfTilingInput.flow = atof(optarg);
		params->tfPlaneParams.flow = atof(optarg);
		if((params->tfTilingInput.flow < 0.0)) {
			sprintf(msg,"must not be negative (%f Hz specified)", params->tfTilingInput.flow);
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
		mdcparams = LALMalloc(sizeof(*mdcparams));
		if(!mdcparams) {
			print_alloc_fail(argv[0], "for mdcparams");
			args_are_bad = TRUE;
		} else {
			mdcparams->channelName = optarg;
			mdcchannelIn.name = optarg;
			mdcchannelIn.type = ADCDataChannel;
			ADD_PROCESS_PARAM("string");
		}
		break;

		case 'T':
		params->tfTilingInput.minFreqBins = atoi(optarg);
		if(params->tfTilingInput.minFreqBins <= 0) {
			sprintf(msg,"must be > 0 (%i specified)", params->tfTilingInput.minFreqBins);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

		case 'U':
		params->tfTilingInput.minTimeBins = atoi(optarg);
		if(params->tfTilingInput.minTimeBins <= 0) {
			sprintf(msg,"must be > 0 (%i specified)", params->tfTilingInput.minTimeBins);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
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
		params->windowLength = atoi(optarg);
		if(params->windowLength <= 0) {
			sprintf(msg, "must be > 0 (%i specified)", params->windowLength);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

		case 'X':
		params->compEPInput.numSigmaMin = atof(optarg);
		if(params->compEPInput.numSigmaMin <= 0.0) {
			sprintf(msg, "must be > 0.0 (%f specified)", params->compEPInput.numSigmaMin);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

		case 'Y':
		if(!strcmp(optarg, "useMean"))
			params->method = useMean;
		else if(!strcmp(optarg, "useMedian"))
			params->method = useMedian;
		else if (!strcmp(optarg, "useUnity"))
			params->method = useUnity;
		else {
			sprintf(msg, "must be \"useMean\", \"useMedian\", or \"useUnity\"");
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

		case 'b':
		if(!strcmp("ldas", optarg))
			resampFiltType = LDASfirLP;
		else if(!strcmp("butterworth", optarg))
			resampFiltType = defaultButterworth;
		else {
			sprintf(msg, "must be \"ldas\", or \"butterworth\"");
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("string");
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
		params->tfTilingInput.overlapFactor = atoi(optarg);
		if(params->tfTilingInput.overlapFactor <= 0) {
			sprintf(msg, "must be > 0 (%i specified)", params->tfTilingInput.overlapFactor);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("int");
		break;

		case 'g':
		params->alphaThreshold = atof(optarg);
		if(params->alphaThreshold <= 0.0) {
			sprintf(msg, "must be > 0.0 (%f specified)", params->alphaThreshold);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

		case 'h':
		options.comment = optarg;
		ADD_PROCESS_PARAM("string");
		break;

		case 'i':
		params->windowType = atoi(optarg);
		if(params->windowType >= NumberWindowTypes) {
			sprintf(msg, "must be < %d (%i specified)", NumberWindowTypes, params->windowType);
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
		params->tfTilingInput.maxTileBand = atof(optarg);
		params->tfPlaneParams.deltaT = 1 / (2 * params->tfTilingInput.maxTileBand);
		if(params->tfTilingInput.maxTileBand < 0) {
			sprintf(msg,"must be >= 0 (%f specified)",params->tfTilingInput.maxTileBand);
			print_bad_argument(argv[0], long_options[option_index].name, msg);
			args_are_bad = TRUE;
		}
		ADD_PROCESS_PARAM("float");
		break;

		case 'm':
		params->tfTilingInput.maxTileDuration = atof(optarg);
		params->tfPlaneParams.deltaF = 1 / (2 * params->tfTilingInput.maxTileDuration);
		if(params->tfTilingInput.maxTileDuration > 1.0) {
			sprintf(msg,"must be <= 1.0 (%f specified)", params->tfTilingInput.maxTileDuration);
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

	if(params->windowLength < 2 * params->windowShift) {
		sprintf(msg, "must be >= 2 * --window-shift = %u (%u specified)", 2 * params->windowShift, params->windowLength);
		print_bad_argument(argv[0], "window-length", msg);
		args_are_bad = TRUE;
	}

	/*
	 * Convert the start and stop times to LIGOTimeGPS structures.
	 */

	if(gpsStartTimeNS > gpsStopTimeNS) {
		fprintf(stderr, "%s: error: GPS start time > GPS stop time\n", argv[0]);
		args_are_bad = TRUE;
	}
	LAL_CALL(LALINT8toGPS(&stat, &options.startEpoch, &gpsStartTimeNS), &stat);
	LAL_CALL(LALINT8toGPS(&stat, &options.stopEpoch, &gpsStopTimeNS), &stat);

	/*
	 * Convert the amount of available RAM to a limit on the length of a
	 * time series to read in.
	 */

	options.maxSeriesLength = (ram * 1024 * 1024) / (4 * sizeof(REAL4));

	/*
	 * Check for missing parameters
	 */

	if(!check_for_missing_parameters(&stat, argv[0], long_options, params))
		args_are_bad = TRUE;

	/*
	 * Exit if anything was wrong with the command line.
	 */

	if(args_are_bad)
		exit(1);

	/*
	 * Ensure PSDAverageLength is comensurate with the analysis window
	 * length and its shift.
	 */

	options.PSDAverageLength = block_commensurate(options.PSDAverageLength, params->windowLength, params->windowShift);

	/*
	 * Ensure RAM limit is comensurate with the PSDAverageLength and
	 * its shift.
	 */

#if 0
	options.maxSeriesLength = block_commensurate(options.maxSeriesLength, options.PSDAverageLength, options.PSDAverageLength - (params.windowLength - params.windowShift));
#endif

	/*
	 * Sanitize filter frequencies.
	 */

	if(options.cal_high_pass > params->tfTilingInput.flow)
		fprintf(stderr, "%s: warning: calibrated data quantization high-pass frequency (%f Hz) greater than TF plane low frequency (%f Hz)\n", argv[0], options.cal_high_pass, params->tfTilingInput.flow);

	if(options.high_pass > params->tfTilingInput.flow - 10.0)
		fprintf(stderr, "%s: warning: data conditioning high-pass frequency (%f Hz) greater than 10 Hz below TF plane low frequency (%f Hz)\n", argv[0], options.high_pass, params->tfTilingInput.flow);

	/*
	 * Miscellaneous chores.
	 */

	params->tfPlaneParams.fhigh = params->tfPlaneParams.flow + params->tfTilingInput.length;
	params->printSpectrum = printSpectrum;
	params->useOverWhitening = useoverwhitening;

	if(options.verbose) {
		fprintf(stderr, "%s: using --psd-average-points %zu\n", argv[0], options.PSDAverageLength);
		fprintf(stderr, "%s: available RAM limits analysis to %d samples\n", argv[0], options.maxSeriesLength);
	}
}


/*
 * ============================================================================
 *                                   Input
 * ============================================================================
 */

static REAL4TimeSeries *get_calibrated_data(
	LALStatus *stat,
	FrStream *stream,
	const char *chname,
	LIGOTimeGPS *start,
	double duration,
	size_t lengthlimit
)
{
	REAL8TimeSeries *calibrated;
	REAL4TimeSeries *series;
	PassBandParamStruc highpassParam;
	unsigned int i;

	/* retrieve calibrated data as REAL8 time series */
	calibrated = XLALFrReadREAL8TimeSeries(stream, chname, start, duration, lengthlimit);

	/* high pass filter before casting REAL8 to REAL4 */
	highpassParam.nMax = 4;
	highpassParam.f2 = options.cal_high_pass;
	highpassParam.f1 = -1.0;
	highpassParam.a2 = 0.9;
	highpassParam.a1 = -1.0;
	LAL_CALL(LALButterworthREAL8TimeSeries(stat, calibrated, &highpassParam), stat);

	/* copy data into a REAL4 time series */
	LAL_CALL(LALCreateREAL4TimeSeries(stat, &series, calibrated->name, calibrated->epoch, calibrated->f0, calibrated->deltaT, calibrated->sampleUnits, calibrated->data->length), stat);
	for(i = 0; i < series->data->length; i++)
		series->data->data[i] = calibrated->data->data[i];
	LAL_CALL(LALDestroyREAL8TimeSeries(stat, calibrated), stat);

	return(series);
}

static REAL4TimeSeries *get_time_series(
	LALStatus *stat,
	const char *dirname,
	const char *cachefile,
	const char *chname,
	LIGOTimeGPS start,
	LIGOTimeGPS end,
	size_t lengthlimit
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
	if(options.calibrated)
		series = get_calibrated_data(stat, stream, chname, &start, duration, lengthlimit);
	else
		series = XLALFrReadREAL4TimeSeries(stream, chname, &start, duration, lengthlimit);

	/* Check for missing data */
	if(stream->state & LAL_FR_GAP) {
		fprintf(stderr, "get_time_series(): error: gap in data detected between GPS times %d.%09d s and %d.%09d s\n", start.gpsSeconds, start.gpsNanoSeconds, end.gpsSeconds, end.gpsNanoSeconds);
		LAL_CALL(LALDestroyREAL4TimeSeries(stat, series), stat);
		series = NULL;
	}

	/* verbosity */
	if(options.verbose)
		fprintf(stderr, "get_time_series(): read %u samples (%.9lf s) at GPS time %u.%09u s\n", series->data->length, series->data->length * series->deltaT, start.gpsSeconds, start.gpsNanoSeconds);

	/* Clean up */
	LAL_CALL(LALFrClose(stat, &stream), stat);

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

	if(!response) {
		fprintf(stderr, "add_burst_injections(): must supply calibration information for injections\n");
		exit(1);
	}

	if(options.verbose)
		fprintf(stderr, "add_burst_injections(): reading in SimBurst Table\n");

	LAL_CALL(LALSimBurstTableFromLIGOLw(stat, &injections, burstInjectionFile, startTime, stopTime), stat);

	if(options.verbose)
		fprintf(stderr, "add_burst_injections(): injecting signals into time series\n");

	LAL_CALL(LALBurstInjectSignals(stat, series, injections, response), stat); 

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

static void add_mdc_injections(LALStatus *stat, const char *mdcCacheFile, REAL4TimeSeries *series)
{
	REAL4TimeSeries *mdc;
	FrStream *frstream = NULL;
	FrCache *frcache = NULL;
	size_t i;

	if(options.verbose)
		fprintf(stderr, "add_mdc_injections(): using MDC frames for injections\n");

	/* open mdc cache */
	LAL_CALL(LALFrCacheImport(stat, &frcache, mdcCacheFile), stat);
	LAL_CALL(LALFrCacheOpen(stat, &frstream, frcache), stat);
	LAL_CALL(LALDestroyFrCache(stat, &frcache), stat);

	/* create and initialize the mdc time series vector */
	LAL_CALL(LALCreateREAL4TimeSeries(stat, &mdc, mdcparams->channelName, series->epoch, series->f0, series->deltaT, series->sampleUnits, series->data->length), stat);

	/* get the mdc signal data */
	LAL_CALL(LALFrGetREAL4TimeSeries(stat, mdc, &mdcchannelIn, frstream), stat);
	mdc->epoch = series->epoch;
	LAL_CALL(LALFrSeek(stat, &mdc->epoch, frstream), stat);
	LAL_CALL(LALFrGetREAL4TimeSeries(stat, mdc, &mdcchannelIn, frstream), stat);

	/* write diagnostic info to disk */
	if(options.printData)
		LALPrintTimeSeries(mdc, "./timeseriesmdc.dat");

	/* add the mdc signal to the given time series */
	for(i = 0; i < series->data->length; i++)
		series->data->data[i] += mdc->data->data[i];

	/* clean up */
	LAL_CALL(LALDestroyREAL4TimeSeries(stat, mdc), stat);
	LAL_CALL(LALFrClose(stat, &frstream), stat);
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
	LALStatus *stat,
	SnglBurstTable **addpoint,
	REAL4TimeSeries *series,
	size_t psdlength,
	EPSearchParams *params
)
{
	REAL4TimeSeries *interval;
	size_t i, start;
	size_t overlap = params->windowLength - params->windowShift;

	if(psdlength > series->data->length) {
		psdlength = block_commensurate(series->data->length, params->windowLength, params->windowShift);
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
		
		LAL_CALL(LALCutREAL4TimeSeries(stat, &interval, series, start, psdlength), stat);

		if(options.verbose)
			fprintf(stderr, "analyze_series(): analyzing samples %zu -- %zu (%.9lf s)\n", start, start + interval->data->length, interval->data->length * interval->deltaT);

		LAL_CALL(EPSearch(stat, interval, params, addpoint), stat);
		while(*addpoint)
			addpoint = &(*addpoint)->next;

		LAL_CALL(LALDestroyREAL4TimeSeries(stat, interval), stat);
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
	LALTimeInterval           overlapCorrection;
	CHAR                      outfilename[256];
	REAL4TimeSeries          *series = NULL;
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

	/* determine the input time series overlap correction (LAL spec
	 * makes us pass a pointer... sigh), and set the outer loop's upper
	 * bound */
	{
	REAL8 overlap = (REAL8) (params.windowLength - params.windowShift) / options.ResampleRate;

	LAL_CALL(LALFloatToInterval(&stat, &overlapCorrection, &overlap), &stat);
	if(options.verbose)
		fprintf(stderr, "%s: time series overlap correction is %u.%09u s\n", argv[0], overlapCorrection.seconds, overlapCorrection.nanoSeconds);

	LAL_CALL(LALAddFloatToGPS(&stat, &boundepoch, &options.stopEpoch, -2.0 * options.FilterCorruption / options.ResampleRate - overlap), &stat);
	if(options.verbose)
		fprintf(stderr, "%s: time series epochs must be less than %u.%09u s\n", argv[0], boundepoch.gpsSeconds, boundepoch.gpsNanoSeconds);
	}


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
		 * Get the data,
		 */

		series = get_time_series(&stat, dirname, cachefile, params.channelName, epoch, options.stopEpoch, options.maxSeriesLength);

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
				LAL_CALL(LALCreateREAL4TimeSeries(&stat, &series, params.channelName, epoch, 0.0, (REAL8) 1.0 / options.ResampleRate, lalADCCountUnit, length), &stat);
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

		if(burstInjectionFile || inspiralInjectionFile) {
			COMPLEX8FrequencySeries  *response = NULL;

			/* Create the response function */
			if(options.calCacheFile)
				response = generate_response(&stat, options.calCacheFile, params.channelName, series->deltaT, series->epoch, series->data->length);


			if(options.printData)
			  LALCPrintFrequencySeries(response, "./response.dat");

			/* perform injections */
			if(burstInjectionFile)
			  add_burst_injections(&stat, series, response);
			else
			  add_inspiral_injections(&stat, series, response);

			if(options.printData)
				LALPrintTimeSeries(series, "./injections.dat");

			/* clean up */
			LAL_CALL(LALDestroyCOMPLEX8FrequencySeries(&stat, response), &stat);
		}

		/*
		 * Add MDC injections into the time series if requested.
		 */

		if(mdcCacheFile) {
			add_mdc_injections(&stat, mdcCacheFile, series);
			if(options.printData)
				LALPrintTimeSeries(series, "./timeseriesasqmdc.dat");
		}

		/*
		 * Condition the time series data.
		 */

		LAL_CALL(EPConditionData(&stat, series, options.high_pass, (REAL8) 1.0 / options.ResampleRate, resampFiltType, options.FilterCorruption), &stat);

		if(options.verbose)
			fprintf(stderr, "%s: %u samples (%.9f s) at GPS time %d.%09d s remain after conditioning\n", argv[0], series->data->length, series->data->length * series->deltaT, series->epoch.gpsSeconds, series->epoch.gpsNanoSeconds);

		/*
		 * Store the start and end times of the data that actually
		 * gets analyzed.
		 */

		if(!searchsumm.searchSummaryTable->out_start_time.gpsSeconds)
			LAL_CALL(LALAddFloatToGPS(&stat, &searchsumm.searchSummaryTable->out_start_time, &series->epoch, series->deltaT * (params.windowLength / 2 - params.windowShift)), &stat);
		LAL_CALL(LALAddFloatToGPS(&stat, &searchsumm.searchSummaryTable->out_end_time, &series->epoch, series->deltaT * (series->data->length - (params.windowLength / 2 - params.windowShift))), &stat);


		/*
		 * Analyze the data
		 */

		EventAddPoint = analyze_series(&stat, EventAddPoint, series, options.PSDAverageLength, &params);

		/*
		 * Reset for next run
		 *
		 * Advancing the epoch by the post-conditioning series length
		 * provides exactly the overlap needed to account for
		 * conditioning corruption.  overlapCorrection provides any
		 * additional overlap that is needed.
		 */

		LAL_CALL(LALAddFloatToGPS(&stat, &epoch, &epoch, series->data->length * series->deltaT), &stat);
		LAL_CALL(LALDecrementGPS(&stat, &epoch, &epoch, &overlapCorrection), &stat);
		LAL_CALL(LALDestroyREAL4TimeSeries(&stat, series), &stat);
	}

	/*
	 * Cluster and sort the events.
	 */

	if(options.cluster)
		LAL_CALL(LALClusterSnglBurstTable(&stat, &burstEvent, XLALCompareSnglBurstByPeakTime, XLALCompareSnglBurstByPeakTimeAndFreq), &stat);
	LAL_CALL(LALSortSnglBurst(&stat, &burstEvent, XLALCompareSnglBurstByStartTimeAndLowFreq), &stat);

	/*
	 * Output the results.
	 */

	snprintf(outfilename, sizeof(outfilename)-1, "%s-POWER_%s-%d-%d.xml", ifo, options.comment, options.startEpoch.gpsSeconds, options.stopEpoch.gpsSeconds - options.startEpoch.gpsSeconds);
	outfilename[sizeof(outfilename)-1] = '\0';
	output_results(&stat, outfilename, &procTable, &procparams, &searchsumm, burstEvent);

	/*
	 * Final cleanup.
	 */

	LALFree(procTable.processTable);
	LALFree(searchsumm.searchSummaryTable);

	while(procparams.processParamsTable) {
		ProcessParamsTable *table = procparams.processParamsTable;
		procparams.processParamsTable = table->next;
		LALFree(table);
	}

	LALFree(mdcparams);

	while(burstEvent) {
		SnglBurstTable *event = burstEvent;
		burstEvent = burstEvent->next;
		LALFree(event);
	}

	LALCheckMemoryLeaks();
	exit(0);
}

#endif
