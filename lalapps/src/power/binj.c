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

struct options {
	INT8 gps_start_time;
	INT8 gps_end_time;
	char *coordinates;
	double flow;
	double fhigh;
	double fratio;
	double deltaf;
	double quality;
	double log_hpeak_min;
	double log_hpeak_max;
	int use_random_strain;
	double log_max_distance;
	double log_min_distance;
	int seed;
	INT8 time_step;	/* nanoseconds between injections */
	char *waveform;
	double simwaveform_duration;
	int simwaveform_min_number;
	int simwaveform_max_number;
	double tau;
	double freq;
	double hpeak;
	char *user_tag;
	int mdc;	/* Need to set this to true if one wants to use MDC signals */
};


static struct options options_defaults(void)
{
	struct options defaults = {
		.gps_start_time = 0,
		.gps_end_time = 0,
		.coordinates = "EQUATORIAL",
		.flow = 150.0,
		.fhigh = 1000.0,
		.fratio = 0.0,
		.deltaf = 0.0,
		.quality = -1.0,
		.log_hpeak_min = 0.0,
		.log_hpeak_max = 0.0,
		.use_random_strain = 0,
		.log_max_distance = log10(10000.0),
		.log_min_distance = log10(100.0),
		.seed = 1,
		.time_step = 210.0 / LAL_PI * 1e9,
		.waveform = "SineGaussian",
		.simwaveform_duration = 0.0,
		.simwaveform_min_number = 0,
		.simwaveform_max_number = 10,
		.tau = 0.1,
		.freq = 150.0,
		.hpeak = 1.0e-20,
		.user_tag = NULL,
		.mdc = FALSE
	};

	return defaults;
}


static void print_usage(const char *prog)
{
	fprintf(stderr, 
"%s [options]\n"\
"\nDefaults are shown in brackets\n\n" \
"  --help                   display this message\n"\
"  --gps-start-time TIME    start injections at GPS time TIME\n"\
"  --gps-end-time TIME      end injections at GPS time TIME\n"\
"  --time-step STEP         space injections STEP / pi seconds appart (210)\n"\
"  --coordinates COORDS     coordinate system to use for injections\n"\
"  --flow FLOW              first frequency of injection (150.0)\n"\
"  --fhigh FHIGH            only inject frequencies smaller than FHIGH (1000.0)\n"\
"  --fratio FACTOR          exponential spacing of injection frequencies (0.0)\n"\
"  --deltaf DELF            linear spacing of injection frequencies (0.0)\n"\
"  --quality Q              quality factor for SG waveforms TAU=Q/(sqrt(2) pi F)\n"\
"  --tau TAU                duration of SG waveforms.  Q overrides TAU setting\n"\
"  --hpeak HPEAK            amplitude of SG injection in strain units\n"\
"  --log-hpeak-min LOGHMIN  min amplitude of SG injection in strain units\n"\
"  --log-hpeak-max LOGHMAX  max amplitude of SG injection in strain units\n"\
"  --min-distance           min distance of source in Kpc(default 100Kpc) \n"\
"  --max-distance           max distance of source in Kpc(default 10000Kpc) \n"\
"  --simwaveform-duration   duration of th esimulated waveform (Warren/Ott/ZM)\n"\
"  --simwaveform-min-number min # of the simulated waveform \n"\
"  --simwaveform-max-number max # of the simulated waveform \n"\
"  --seed SEED              seed random number generator with SEED (1)\n"\
"  --waveform NAME          set waveform type to NAME (SineGaussian)\n"\
"  --user-tag STRING        set the usertag to STRING\n\n", prog);
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
		{"coordinates", required_argument, NULL, 'A'},
		{"deltaf", required_argument, NULL, 'B'},
		{"fhigh", required_argument, NULL, 'C'},
		{"flow", required_argument, NULL, 'D'},
		{"fratio", required_argument, NULL, 'E'},
		{"freq", required_argument, NULL, 'F'},
		{"gps-end-time", required_argument, NULL, 'G'},
		{"gps-start-time", required_argument, NULL, 'H'},
		{"help", no_argument, NULL, 'I'},
		{"hpeak", required_argument, NULL, 'J'},
		{"log-hpeak-max", required_argument, NULL, 'K'},
		{"log-hpeak-min", required_argument, NULL, 'L'},
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
		options.coordinates = optarg;
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

	case 'F':
		options.freq = atof(optarg);
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
		options.hpeak = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'K':
		options.use_random_strain += 1;
		options.log_hpeak_max = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'L':
		options.use_random_strain += 1;
		options.log_hpeak_min = atof(optarg);
		ADD_PROCESS_PARAM("double");
		break;

	case 'M':
		options.log_max_distance = log10(atof(optarg));
		ADD_PROCESS_PARAM("double");
		break;

	case 'N':
		options.log_min_distance = log10(atof(optarg));
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
	if(strlen(options.coordinates) > LIGOMETA_COORDINATES_MAX) {
		fprintf(stderr, "error: --coordinates %s exceeds max length of %d\n", options.waveform, LIGOMETA_COORDINATES_MAX);
		exit(1);
	}
	if(strlen(options.waveform) > LIGOMETA_WAVEFORM_MAX) {
		fprintf(stderr, "error: --waveform %s exceeds max length of %d\n", options.waveform, LIGOMETA_WAVEFORM_MAX);
		exit(1);
	}
	if(options.use_random_strain && (options.use_random_strain != 2)) {
		fprintf(stderr, "Must supply upper and lower limits when using random strain\n");
		exit(1);
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
 *                          Arrival Time Calculation
 * ============================================================================
 */

static LIGOTimeGPS arrival_time(LALDetector detector, LIGOTimeGPS geocent_peak_time, double latitude, double longitude)
{
	LALStatus status = blank_status;
	LALPlaceAndGPS place_and_gps = {
		.p_detector = &detector,
		.p_gps = &geocent_peak_time
	};
	SkyPosition sky_pos = {
		.longitude = longitude,
		.latitude = latitude,
		.system = COORDINATESYSTEM_EQUATORIAL
	};
	DetTimeAndASource det_time_and_source = {
		.p_det_and_time = &place_and_gps,
		.p_source = &sky_pos
	};
	REAL8 dt;

	LAL_CALL(LALTimeDelayFromEarthCenter(&status, &dt, &det_time_and_source), &status);

	return *XLALGPSAdd(&geocent_peak_time, dt);
}


/* 
 * ============================================================================
 *                                   Output
 * ============================================================================
 */

/* output for LIGO-TAMA simulations */
static void ligo_tama_output(FILE * fpout, SimBurstTable * simBursts)
{
	SimBurstTable *thisEvent = NULL;

	thisEvent = simBursts;
	fprintf(fpout, "# $I" "d$\n");
	fprintf(fpout, "# %s\n", thisEvent->waveform);
	fprintf(fpout, "# geocent_peak_time\tnSec\tdtminus\t\tdtplus\t\tlongitude\tlatitude\tpolarization\tcoordinates\thrss\thpeak\tfreq\ttau\n");
	fflush(fpout);

	while(thisEvent) {
		fprintf(fpout, "%0d\t%0d\t%f\t%f\t%f\t%f\t%f\t%s\t%e\t%e\t%f\t%f\n", thisEvent->geocent_peak_time.gpsSeconds, thisEvent->geocent_peak_time.gpsNanoSeconds, thisEvent->dtminus, thisEvent->dtplus, thisEvent->longitude, thisEvent->latitude, thisEvent->polarization, thisEvent->coordinates, thisEvent->hrss, thisEvent->hpeak, thisEvent->freq, thisEvent->tau);
		thisEvent = thisEvent->next;
	}

	fprintf(fpout, "# $I" "d$\n");
}

static void write_tamma(MetadataTable injections, struct options options)
{
	FILE *fpout;
	char fname[256];

	if(options.user_tag) {
		LALSnprintf(fname, sizeof(fname), "HLT-INJECTIONS_%s-%d-%d.txt", options.user_tag, (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));
	} else {
		LALSnprintf(fname, sizeof(fname), "HLT-INJECTIONS-%d-%d.txt", (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));
	}
	fpout = fopen(fname, "w");
	ligo_tama_output(fpout, injections.simBurstTable);
	fclose(fpout);
}


/* LIGO LW XML of MDC injections */
static void write_mdc_xml(MetadataTable mdcinjections)
{
	LALStatus status = blank_status;
	CHAR fname[256];
	LIGOLwXMLStream xmlfp;

	memset(&xmlfp, 0, sizeof(xmlfp));

	LALSnprintf(fname, sizeof(fname), "HL-MDCSG10_%d.xml", 1);

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
	LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
	LIGOLwXMLStream xmlfp;

	memset(&xmlfp, 0, sizeof(xmlfp));

	if(options.user_tag)
		snprintf(fname, sizeof(fname), "HL-INJECTIONS_%s-%d-%d.xml", options.user_tag, (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));
	else
		snprintf(fname, sizeof(fname), "HL-INJECTIONS-%d-%d.xml", (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));


	LAL_CALL(LALOpenLIGOLwXMLFile(&status, &xmlfp, fname), &status);

	LAL_CALL(LALGPSTimeNow(&status, &proctable.processTable->end_time, &accuracy), &status);
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
	REAL4 deviate = 0.0;
	int iamp = 1;
	/* observatory information */
	LALDetector lho = lalCachedDetectors[LALDetectorIndexLHODIFF];
	LALDetector llo = lalCachedDetectors[LALDetectorIndexLLODIFF];
	/* xml output data */
	LALStatus status = blank_status;
	LALLeapSecAccuracy accuracy = LALLEAPSEC_LOOSE;
	/* tables */
	MetadataTable proctable;
	MetadataTable procparams;
	MetadataTable injections;
	MetadataTable mdcinjections;
	SimBurstTable *this_sim_burst = NULL;

	/* set up inital debugging values */
	lal_errhandler = LAL_ERR_EXIT;
	set_debug_level("LALMSGLVL2");

	/* create the process and process params tables */
	proctable.processTable = calloc(1, sizeof(ProcessTable));
	LAL_CALL(LALGPSTimeNow(&status, &(proctable.processTable->start_time), &accuracy), &status);
	LAL_CALL(populate_process_table(&status, proctable.processTable, PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE), &status);
	LALSnprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, " ");
	procparams.processParamsTable = NULL;

	/* parse command line */
	options = parse_command_line(&argc, &argv, &procparams);
	options.freq = options.flow;

	/* fill the comment */
	if(options.user_tag)
		snprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, "%s", options.user_tag);
	else
		snprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, "");

	/* initialize random number generator */
	LAL_CALL(LALCreateRandomParams(&status, &randParams, options.seed), &status);

	/* if we are doing string cusps; we want to inject the correct
	   population of frequencies; what's distributed uniformly is
	   theta^2 (the square of the angle the line of sight makes
	   with direction of the cusp) */
	if(!strcmp(options.waveform, "StringCusp")) {
		REAL4 thetasqmin = pow(options.fhigh, -2. / 3.), thetasqmax = pow(options.flow, -2. / 3.), thetasq;

		fprintf(stdout, "Generating cusp population\n");

		srand(options.seed);
		thetasq = (thetasqmax - thetasqmin) * ((float) rand() / (float) RAND_MAX) + thetasqmin;
		options.freq = pow(thetasq, -3. / 2.);
	}

	for(tinj = options.gps_start_time; tinj <= options.gps_end_time; tinj += options.time_step) {
		/* compute tau if quality was specified */
		if(options.quality > 0.0)
			options.tau = options.quality / (sqrt(2.0) * LAL_PI * options.freq);

		/* allocate the injection */
		if(!this_sim_burst)
			this_sim_burst = injections.simBurstTable = calloc(1, sizeof(SimBurstTable));
		else {
			this_sim_burst->next = calloc(1, sizeof(SimBurstTable));
			this_sim_burst = this_sim_burst->next;
		}

		/* GPS time of burst */
		XLALINT8NSToGPS(&this_sim_burst->geocent_peak_time, tinj);

		/* GMST of burst in hours */
		this_sim_burst->peak_time_gmst = XLALGreenwichMeanSiderealTime(&this_sim_burst->geocent_peak_time) * 12.0 / LAL_PI;

		/* populate the sim burst table */
		snprintf(this_sim_burst->waveform, LIGOMETA_WAVEFORM_MAX, "%s", options.waveform);
		snprintf(this_sim_burst->coordinates, LIGOMETA_COORDINATES_MAX, "%s", options.coordinates);

		/* sky location and polarizatoin angle */
		if(!strcmp(options.coordinates, "ZENITH")) {
			/* zenith */
			this_sim_burst->longitude = 0.0;
			this_sim_burst->latitude = 0.0;
			this_sim_burst->polarization = 0.0;
		} else {
			LAL_CALL(LALUniformDeviate(&status, &deviate, randParams), &status);
			this_sim_burst->longitude = 2.0 * LAL_PI * deviate;

			LAL_CALL(LALUniformDeviate(&status, &deviate, randParams), &status);
			this_sim_burst->latitude = LAL_PI / 2.0 - acos(2.0 * deviate - 1.0);

			LAL_CALL(LALUniformDeviate(&status, &deviate, randParams), &status);
			this_sim_burst->polarization = 2.0 * LAL_PI * deviate;
		}

		/* compute amplitude information */
		if(options.use_random_strain) {
#if 1
			LAL_CALL(LALUniformDeviate(&status, &deviate, randParams), &status);
			options.hpeak = pow(10, (options.log_hpeak_max - options.log_hpeak_min) * deviate + options.log_hpeak_min);

			/*
			 * Uncomment the lines below and comment the above
			 * if want to inject uniformly spaced amplitude
			 * waveforms
			 */
#else
			if(iamp == 1) {
				hpeak = pow(10, -19);
				iamp++;
			} else if(iamp == 2) {
				hpeak = 2.0 * pow(10, -19);
				iamp++;
			} else if(iamp == 3) {
				hpeak = 3.0 * pow(10, -19);
				iamp++;
			} else if(iamp == 4) {
				hpeak = 3.3 * pow(10, -18);
				iamp++;
			} else if(iamp == 5) {
				hpeak = 5.3 * pow(10, -17);
				iamp = 1;
			}
#endif
		}

		/* uniform distribution in log(distance) */
		LAL_CALL(LALUniformDeviate(&status, &deviate, randParams), &status);
		this_sim_burst->distance = pow(10.0, options.log_min_distance + (options.log_max_distance - options.log_min_distance) * deviate);


		/* deal with the intrinsic signal parameters */
		this_sim_burst->hrss = sqrt(sqrt(2.0 * LAL_PI) * options.tau / 4.0) * options.hpeak;
		this_sim_burst->hpeak = options.hpeak;
		this_sim_burst->freq = options.freq;
		this_sim_burst->tau = options.tau;
		if(options.simwaveform_duration) {
			this_sim_burst->dtplus = options.simwaveform_duration / 2.0;
			this_sim_burst->dtminus = options.simwaveform_duration / 2.0;
		} else {
			this_sim_burst->dtplus = 4.0 * options.tau;
			this_sim_burst->dtminus = 4.0 * options.tau;
		}

		/* set the simulated wavenumber */
		LAL_CALL(LALUniformDeviate(&status, &deviate, randParams), &status);
		this_sim_burst->zm_number = options.simwaveform_min_number + (options.simwaveform_max_number - options.simwaveform_min_number) * deviate;

		/* arrival times */
		this_sim_burst->h_peak_time = arrival_time(lho, this_sim_burst->geocent_peak_time, this_sim_burst->latitude, this_sim_burst->longitude);
		this_sim_burst->l_peak_time = arrival_time(llo, this_sim_burst->geocent_peak_time, this_sim_burst->latitude, this_sim_burst->longitude);

		/* increment to next frequency and test it's still in band */
		if((options.deltaf == 0.0) && (options.fratio != 0.0))
			options.freq *= options.fratio;
		else if((options.deltaf != 0.0) && (options.fratio == 0.0))
			options.freq += options.deltaf;
		else {
			fprintf(stderr, "error: something wrong with --deltaf and --fratio\n");
			exit(1);
		}
		if(options.freq > options.fhigh)
			options.freq = options.flow;

		/* if we are doing string cusps; we want to inject the correct
		   population of frequencies */
		if(!strcmp(options.waveform, "StringCusp")) {
			REAL4 thetasqmin = pow(options.fhigh, -2. / 3.), thetasqmax = pow(options.flow, -2. / 3.), thetasq;

			thetasq = (thetasqmax - thetasqmin) * ((float) rand() / (float) RAND_MAX) + thetasqmin;
			options.freq = pow(thetasq, -3. / 2.);
		}
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
		write_tamma(injections, options);
	}

	exit(0);
}
