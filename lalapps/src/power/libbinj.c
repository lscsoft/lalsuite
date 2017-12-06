/*
 * Copyright (C) 2007 Kipp Cannon, Lisa M. Goggin, Patrick Brady, Saikat
 * Ray-Majumder, Xavier Siemens, Salvatore Vitale
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
#include <lal/LALSimInspiral.h>
#include <lal/LALStdlib.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOLwXMLBurstRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOMetadataBurstUtils.h>
#include <lal/TimeDelay.h>
#include <lal/TimeSeries.h>
#include <lal/XLALError.h>
#include <lal/LALSimNoise.h>
#include <lalapps.h>
#include <processtable.h>
#include <lal/LALSimulation.h>
#include <LALAppsVCSInfo.h>
#include <series.h>
#include <lal/LALDatatypes.h>
#include <lal/FrequencySeries.h>
#include <lal/TimeFreqFFT.h>
#include <lal/DetResponse.h>
#include <lal/Units.h>
#include <lal/LALFrStream.h>
#include <lal/LALFrameIO.h>
#include <lal/LALCache.h>

#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_libbinj"



double q_min=2.0;

/* Variables and functions below are for SNR calculation, if required */
char *snr_ifos      = NULL; 
double single_IFO_SNR_threshold=0.0;
REAL8 ligoStartFreq=-1.0;
REAL8 virgoStartFreq=-1.0;
CHAR *ligoFakePsd=NULL;
CHAR *virgoFakePsd=NULL;
char ** ifonames=NULL;
int numifos=0;
CHAR *ligoPsdFileName   = NULL;
CHAR *virgoPsdFileName  = NULL;
REAL8FrequencySeries *ligoPsd  = NULL;
REAL8FrequencySeries *virgoPsd = NULL;

static REAL8 calculate_SineGaussian_snr(SimBurst *inj, char *IFOname, REAL8FrequencySeries *psd, REAL8 start_freq);
static void get_FakePsdFromString(REAL8FrequencySeries* PsdFreqSeries,char* FakePsdName, REAL8 StartFreq);

struct fvec *interpFromFile(char *filename);
struct fvec {
	REAL8 f;
	REAL8 x;
};
REAL8 interpolate(struct fvec *fvec, REAL8 f);
struct fvec *ligo_interp;
struct fvec *virgo_interp;

/* .. */

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
  POPULATION_ALL_SKY_GAUSSIAN,
	POPULATION_ALL_SKY_SINEGAUSSIAN_F,
	POPULATION_ALL_SKY_BTLWNB,
	POPULATION_STRING_CUSP,
  POPULATION_ALL_SKY_DAMPEDSINUSOID
};

typedef enum tagParDistr {
	
	FIXED,
	UNIFORM,
	UNIFORM_LOG,
	VOLUME,
	GAUSSIAN,
	NUM_ELEMENTS=5
	} ParDistr;
	
static int parse_distr(char* opt){
	
	if( strstr(opt,"constant")|| strstr(opt,"fixed"))
		return 0;
	else if( strstr(opt,"uniform"))
		return 1;
	else if ( strstr(opt,"log"))
		return 2;
	else if ( strstr(opt,"volume"))
		return 3;
	else if( strstr(opt,"gaussian"))
		return 4;
  return 0;

}

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
  ParDistr tau_distr;
  double tau_stdev;
  double tau;
	double maxEoverr2;
	double minEoverr2;
	double maxf;
	double minf;
	double maxhrss;
	double minhrss;
	char *output;
	double q;
	double minq;
	double maxq;
	ParDistr q_distr;
	ParDistr hrss_distr;
	ParDistr f_distr;
	double q_stdev;
	double hrss_stdev;
	double f_stdev;
	unsigned long seed;
	double time_step;
	double jitter;
	char *time_slide_file;
	char *user_tag;
	double f;
	double hrss;
	//ParDistr snr_distr;
	double minsnr;
	double maxsnr;
	char ** ifonames;
	int nIFO;
  ParDistr polee_distr;
  double minpolee;
  double maxpolee;
  double pol_ellipse_e;
  ParDistr polea_distr;
  double minpolea;
  double maxpolea;
  double pol_ellipse_angle;
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
  defaults.tau=XLAL_REAL8_FAIL_NAN;
  defaults.tau_stdev=-1.0;
  defaults.tau_distr=FIXED;
	defaults.maxEoverr2 = XLAL_REAL8_FAIL_NAN;
	defaults.minEoverr2 = XLAL_REAL8_FAIL_NAN;
	defaults.maxf = XLAL_REAL8_FAIL_NAN;
	defaults.minf = XLAL_REAL8_FAIL_NAN;
	defaults.maxhrss = XLAL_REAL8_FAIL_NAN;
	defaults.minhrss = XLAL_REAL8_FAIL_NAN;
	defaults.minA = XLAL_REAL8_FAIL_NAN;
	defaults.maxA = XLAL_REAL8_FAIL_NAN;
	defaults.output = NULL;
	defaults.seed = 0;
	defaults.time_step = 210.0 / LAL_PI;
	defaults.jitter = 0.0;
	defaults.time_slide_file = NULL;
	defaults.user_tag = NULL;
	defaults.q_distr=FIXED;
	defaults.f_distr=UNIFORM_LOG;
	defaults.hrss_distr=NUM_ELEMENTS+1;
	defaults.minq=XLAL_REAL8_FAIL_NAN;
	defaults.maxq=XLAL_REAL8_FAIL_NAN;
	defaults.q_stdev=-1.0;
	defaults.f_stdev=-1.0;
	defaults.hrss_stdev=-1.0;
	defaults.q = XLAL_REAL8_FAIL_NAN;
	defaults.f=XLAL_REAL8_FAIL_NAN;
	defaults.hrss=XLAL_REAL8_FAIL_NAN;
	//defaults.snr_distr=NUM_ELEMENTS+1;
	defaults.minsnr=.0;
	defaults.maxsnr=0.;
	defaults.ifonames=NULL;
	defaults.nIFO=0;
  defaults.polee_distr=NUM_ELEMENTS+1;
  defaults.pol_ellipse_e=XLAL_REAL8_FAIL_NAN;
  defaults.minpolee=XLAL_REAL8_FAIL_NAN;
  defaults.maxpolee=XLAL_REAL8_FAIL_NAN;
  defaults.polea_distr=NUM_ELEMENTS+1;
  defaults.pol_ellipse_angle=XLAL_REAL8_FAIL_NAN;
  defaults.minpolea=XLAL_REAL8_FAIL_NAN;
  defaults.maxpolea=XLAL_REAL8_FAIL_NAN;
	return defaults;
}

static REAL8  calculate_NetSNR(SimBurst *inj,char ** IFOnames, REAL8FrequencySeries **psds,REAL8 *start_freqs, struct options *options, gsl_rng *rng );


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
"	Set the bounds of the injection durations (btlwnb waveforms) or the min and max duration of Gaussians.\n" \
"\n" \
"--max-e-over-r2 value\n" \
"--min-e-over-r2 value\n" \
"	Set the bounds of the range of equivalent isotropic radiated\n" \
"	energies of btlwnb waveforms.  The units are M_{sun} / pc^{2} (solar\n" \
"	masses per parsec^{2}).\n" \
"\n" \
	); fprintf(stderr, 
"[--f-distr constant,uniform,log,gaussian] \n" \
"   Centre frequency distribution of SineGaussian injections:\n" \
"   constant: all the signal have the same frequency (requires --f)\n" \
"   uniform:  uniform distribution (requires --min-frequency and --max-frequency)\n" \
"   log:      log distribution (requires --min-frequency and --max-frequency)\n" \
"   gaussian:  gaussian distribution (requires --f and --f-stdev)\n" \
"\n" \
"--max-frequency hertz\n" \
"--min-frequency hertz\n" \
"	Set the bounds of the injection frequencies (optional for SineGaussian unless f-distr is uniform or log).\n" \
"	These are the centre frequencies of sine-Gaussians and btlwnb waveforms, and the\n" \
"	high-frequency cut-offs of string cusp waveforms.\n" \
"\n" \
"[--f hertz ] For gaussian distribution of centre frequency, mean of the distribution.\n  \t      If --f-distr constant, constant value of centre frequency\n"\
"[--f-stdev  hertz] For gaussian distribution of centre frequency, standard deviation of the distribution.\n" \
"\n" \
"[--hrss-distr constant,uniform,log,gaussian,volume] \n" \
"   hrss distribution of (Sine)Gaussian injections:\n" \
"   constant: all the signal have the same hrss (requires --hrss)\n" \
"   uniform:  uniform distribution (requires --min-hrss and --max-hrss)\n" \
"   log:      log distribution (requires --min-hrss and --max-hrss)\n" \
"   gaussian:  gaussian distribution (requires --hrss and --hrss-stdev)\n" \
"   volume:   distribution uniform in volume (requires --min-hrss and --max-hrss)\n" \
"\n" \
"[--max-hrss value]\n" \
"[--min-hrss value]\n" \
"\n"\
"	Set the bounds of the injection h_{rss} values.  These only affect\n" \
"	(sine)Gaussian injections.\n" \
"\n" \
"[--hrss value ] For gaussian distribution of hrss, mean of the distribution.\n"\
"[--hrss-stdev value ] For gaussian distribution of hrss, standard deviation of the distribution.\n"\
"\n" \
"[--max-hrss value]\n" \
"[--min-hrss value]\n" \
"\n"\
"	Set the bounds of the injection h_{rss} values.  These only affect\n" \
"	(sine)Gaussian injections. Will inherit hrss distribution \n" \
"\n" \
"[--q-distr constant,uniform,log,gaussian] \n" \
"   Q distribution of SineGaussian injections:\n" \
"   constant: all the signal have the same Q (requires --q)\n" \
"   uniform:  uniform distribution (requires --min-q and --max-q)\n" \
"   log:      log distribution (requires --min-q and --max-q)\n" \
"   gaussian:  gaussian distribution (requires --q and --q-stdev)\n" \
"\n" \
"[--max-q] value \n" \
"[--min-q] value \n" \
"	Set the bounds of the injection Qs.\n" \
"\n" \
"[--q value ] For gaussian distribution of Q, mean of the distribution.\n  \t      If --q-distr constant, constant value of Q\n" \
"[--q-stdev  value] For gaussian distribution of Q, standard deviation of the distribution.\n"\
"\n" \
"[--polar-angle-distr constant,uniform] \n"\
"   Distribution of polar ellipse angle. Default: linearly polarized signals\n" \
"   constant: all the signal have the same polar_ellipse_angle (requires --polar-ellipse-angle)\n" \
"   uniform:  uniform distribution (requires --min-polar-angle and --max-polar-angle)\n" \
"\n" \
"[--max-polar-angle] value \n" \
"[--min-polar-angle] value \n" \
"	Set the bounds of the injection ellipse_polar_angle.\n" \
"[--polar-eccentricity-distr constant,uniform] \n"\
"   Distribution of polar ellipse eccentricity. Default: linearly polarized signals\n" \
"   constant: all the signal have the same polar_ellipse_e (requires --polar-ellipse-e)\n" \
"   uniform:  uniform distribution (requires --min-polar-eccentricity and --max-polar-eccentricity)\n" \
"\n" \
"[--max-polar-eccentricity] value \n" \
"[--min-polar-eccentricity] value \n" \
"	Set the bounds of the injection ellipse_polar_eccentricity.\n" \
	); fprintf(stderr, 
"--output filename\n" \
"	Select output name (default is too hard to explain).\n" \
"\n" \
"--population name\n" \
"	Select the injection population to synthesize.  Allowed values are\n" 
"	\"targeted\", \"string_cusp\", \"all_sky_sinegaussian\", \"all_sky_gaussian\",\n" \
"	\"all_sky_sinegaussian_F\", \"all_sky_btlwnb\", \"all_sky_DampedSinusoid\".\n" \
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
"	Set the time betwen injections in seconds (default = 210 / pi).\n" \
"\n" \
"--jitter value\n" \
"	Give the injection time a random offset within a centered window with a\n" \
"	length specified by this parameter. Value is in seconds, default is 0.\n" \
"\n" \
"--time-slide-file filename\n" \
"	Set the name of the LIGO Light-Weight XML file from which to load\n" \
"	the time slide table.  The document must contain exactly 1 time\n" \
"	slide vector, and only the contents of the process, process_params,\n" \
"	search_summary (optional), sim_burst (optional), and time_slide tables\n" \
"	will be copied into " PROGRAM_NAME "'s output.\n" \
"\n" \
"--user-tag string\n" \
"	Set the user tag in the process and search summary tables to this.\n"
	);
}


static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param, const ProcessTable *process, const char *type, const char *param, const char *value)
{
	*proc_param = XLALCreateProcessParamsTableRow(process);
	snprintf((*proc_param)->program, sizeof((*proc_param)->program), "%s", PROGRAM_NAME);
	snprintf((*proc_param)->type, sizeof((*proc_param)->type), "%s", type);
	snprintf((*proc_param)->param, sizeof((*proc_param)->param), "--%s", param);
	snprintf((*proc_param)->value, sizeof((*proc_param)->value), "%s", value);

	return &(*proc_param)->next;
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
		{"output", required_argument, NULL, 'V'},
		{"population", required_argument, NULL, 'N'},
		{"q", required_argument, NULL, 'O'},
		{"min-q", required_argument,NULL,1707},
		{"max-q", required_argument,NULL,1708},
		{"q-distr", required_argument,NULL,1709},
		{"f-distr", required_argument,NULL,1710},
		{"hrss-distr", required_argument,NULL,1711},
		{"q-stdev", required_argument,NULL,1712},
		{"f-stdev", required_argument,NULL,1713},
		{"hrss-stdev", required_argument,NULL,1714},
		{"f", required_argument,NULL,1715},
		{"hrss", required_argument,NULL,1716},
		{"min-snr", required_argument,NULL,1718},
		{"max-snr", required_argument,NULL,1719},
		{"ifos", required_argument, NULL,1720},
		{"ligo-start-freq",required_argument,NULL,1721},
		{"ligo-fake-psd",  required_argument, 0, 1722},
		{"virgo-fake-psd",  required_argument, 0, 1723},
		{"virgo-start-freq", required_argument, 0, 1724},
		{"ligo-psd",  required_argument, 0, 1725},
		{"virgo-psd",  required_argument, 0, 1726},
		{"ra-dec", required_argument, NULL, 'U'},
		{"seed", required_argument, NULL, 'P'},
		{"time-step", required_argument, NULL, 'Q'},
		{"time-slide-file", required_argument, NULL, 'W'},
		{"jitter", required_argument, NULL, 'X'},
		{"user-tag", required_argument, NULL, 'R'},
    {"tau-distr",required_argument,NULL,1727},
    {"tau-stdev",required_argument,NULL,1728},
    {"tau",required_argument,NULL,1729},
    {"min-polar-angle",required_argument,NULL,1731},
    {"max-polar-angle",required_argument,NULL,1732},
    {"polar-angle-distr",required_argument,NULL,1733},
    {"polar-ellipse-angle",required_argument,NULL,1734},
    {"min-polar-eccentricity",required_argument,NULL,1735},
    {"max-polar-eccentricity",required_argument,NULL,1736},
    {"polar-eccentricity-distr",required_argument,NULL,1737},
    {"polar-ellipse-e",required_argument,NULL,1738},
    
		{NULL, 0, NULL, 0}
	};
	int optarg_len=0;
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
		else if(!strcmp(optarg, "all_sky_sinegaussian_F"))
			options.population = POPULATION_ALL_SKY_SINEGAUSSIAN_F;
		else if(!strcmp(optarg, "all_sky_sinegaussian"))
			options.population = POPULATION_ALL_SKY_SINEGAUSSIAN;
    else if(!strcmp(optarg, "all_sky_gaussian"))
			options.population = POPULATION_ALL_SKY_GAUSSIAN;
    else if(!strcmp(optarg, "all_sky_DampedSinusoid"))
      options.population = POPULATION_ALL_SKY_DAMPEDSINUSOID;
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
		options.time_step = atof(optarg);
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

	case 'V':
		options.output = optarg;
		break;

	case 'W':
		options.time_slide_file = optarg;
		ADD_PROCESS_PARAM(process, "lstring");
		break;

	case 'X':
		options.jitter = atof(optarg);
		ADD_PROCESS_PARAM(process, "lstring");
		break;
	case 1707:
		options.minq=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
	case 1708:
		options.maxq=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
	case 1709:
		options.q_distr=(ParDistr) parse_distr(optarg);
		ADD_PROCESS_PARAM(process, "int_4s");
		break;
	case 1710:
		options.f_distr=parse_distr(optarg);
		ADD_PROCESS_PARAM(process, "int_4s");
		break;
	case 1711:
		options.hrss_distr=parse_distr(optarg);
		ADD_PROCESS_PARAM(process, "int_4s");
		break;
  case 1727:
		options.tau_distr=parse_distr(optarg);
		ADD_PROCESS_PARAM(process, "int_4s");
		break;
  case 1728:
		options.tau_stdev=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
  case 1729:
		options.tau=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;    
	case 1712:
		options.q_stdev = atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
	case 1713:
		options.f_stdev=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
	case 1714:
		options.hrss_stdev=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
	case 1715:
		options.f=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
	case 1716:
		options.hrss=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
	/*case 1717:
		options.snr_distr=parse_distr(optarg);
		ADD_PROCESS_PARAM(process, "int_4s");
		break;*/
	case 1718:
		options.minsnr=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
	case 1719:
		options.maxsnr=atof(optarg);
		ADD_PROCESS_PARAM(process, "real_8");
		break;
	case 1720:
		optarg_len = strlen( optarg ) + 1;
    snr_ifos       = calloc( 1, optarg_len * sizeof(char) );
    memcpy( snr_ifos, optarg, optarg_len * sizeof(char) );
    break;
	case 1721:
		ligoStartFreq = (REAL8) atof( optarg );
		break;
	case 1722:
		optarg_len      = strlen( optarg ) + 1;
    ligoFakePsd = calloc( 1, optarg_len * sizeof(char) );
    memcpy( ligoFakePsd, optarg, optarg_len * sizeof(char) );
    break;
  case 1724:
		virgoStartFreq = (REAL8) atof( optarg );
		break;
	case 1723:
		optarg_len      = strlen( optarg ) + 1;
    virgoFakePsd = calloc( 1, optarg_len * sizeof(char) );
    memcpy( virgoFakePsd, optarg, optarg_len * sizeof(char) );
    break;
	case 1725:
    optarg_len      = strlen( optarg ) + 1;
    ligoPsdFileName = calloc( 1, optarg_len * sizeof(char) );
    memcpy( ligoPsdFileName, optarg, optarg_len * sizeof(char) );
    break;
	case 1726:
		optarg_len       = strlen( optarg ) + 1;
    virgoPsdFileName = calloc( 1, optarg_len * sizeof(char) );
    memcpy( virgoPsdFileName, optarg, optarg_len * sizeof(char) );
    break;
  case 1731:
		options.minpolea= (REAL8) atof( optarg );
		ADD_PROCESS_PARAM(process, "real_8");
		break;
  case 1732:
		options.maxpolea= (REAL8) atof( optarg );
		ADD_PROCESS_PARAM(process, "real_8");
		break;
  case 1733:
		options.polea_distr=parse_distr(optarg);
		ADD_PROCESS_PARAM(process, "int_4s");
		break;
  case 1734:
		options.pol_ellipse_angle= (REAL8) atof( optarg );
		ADD_PROCESS_PARAM(process, "real_8");
		break;
  case 1735:
		options.minpolee= (REAL8) atof( optarg );
		ADD_PROCESS_PARAM(process, "real_8");
		break;
  case 1736:
		options.maxpolee= (REAL8) atof( optarg );
		ADD_PROCESS_PARAM(process, "real_8");
		break;
  case 1737:
		options.polee_distr=parse_distr(optarg);
		ADD_PROCESS_PARAM(process, "int_4s");
		break;
  case 1738:
		options.pol_ellipse_e= (REAL8) atof( optarg );
		ADD_PROCESS_PARAM(process, "real_8");
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
	if(options.maxq < options.minq) {
		fprintf(stderr, "error: --max-q < --min-q\n");
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
	if(!options.time_slide_file) {
		fprintf(stderr, "--time-slide-file is required\n");
		exit(1);
	}
	/* Check if gaussian distributions for f, q or hrss are wanted, and, if so, that a value for the standard deviation is provided */
	if (options.q_distr == GAUSSIAN){
		if (options.q_stdev==-1){
			fprintf(stderr,"Must provide a stdev value for Gaussian injections of q with --q-stdev. Exiting...\n");
			exit(1);
		}
		if (options.q_stdev<=0 ){
			fprintf(stderr,"The standard deviation of the Gaussian distribution of q must be larger than 0.\n");
			exit(1);
		}
		if (options.q==XLAL_REAL8_FAIL_NAN){
			fprintf(stderr,"Must provide a value of q for Gaussian injections. Use --q. Exiting...\n");
			exit(1);
		}
	}
	if (options.f_distr == GAUSSIAN){
		if (options.f_stdev==-1){
			fprintf(stderr,"Must provide a stdev value for Gaussian injections of frequency with --f-stdev. Exiting...\n");
			exit(1);
		}
		if (options.f==XLAL_REAL8_FAIL_NAN){
			fprintf(stderr,"Must provide a value of f for Gaussian injections. Use --f. Exiting...\n");
			exit(1);
		}
	}
	if (options.hrss_distr == GAUSSIAN){
		if (options.hrss_stdev==-1){
			fprintf(stderr,"Must provide a stdev value for Gaussian injections of hrss with --hrss-stdev. Exiting...\n");
			exit(1);
		}
		if (options.hrss==XLAL_REAL8_FAIL_NAN){
			fprintf(stderr,"Must provide a value of hrss for Gaussian injections. Use --q. Exiting...\n");
			exit(1);
		}
	}
	if (options.hrss_distr > NUM_ELEMENTS){
		options.hrss_distr=VOLUME;
		fprintf(stderr, "Did not provide hrss-distr. Using default distribution of hrss (uniform in volume) \n");
		
	}
	if (options.minsnr>0. && options.maxsnr>0.){
		
		 /* Check that each ifo has its PSD file or fakePSD */
    char *tmp, *ifo;

    /* Get the number and names of IFOs */
    tmp = LALCalloc(1, strlen(snr_ifos) + 1);
    strcpy(tmp, snr_ifos);
    ifo = strtok (tmp,",");
    while (ifo != NULL)
    {
        numifos += 1;
        ifo= strtok (NULL, ",");
    }
    ifonames=realloc(ifonames,(numifos+2)*sizeof(CHAR **));
    strcpy(tmp, snr_ifos);
    ifo = strtok (tmp,",");
    ifonames[0]=malloc(strlen(ifo)+1);
    sprintf(ifonames[0],"%s",ifo);
    int i=1;
    while (ifo != NULL)
    {
        ifo       = strtok (NULL, ",");
        if (ifo!=NULL){
            ifonames[i]=malloc(strlen(ifo)+1);
            sprintf(ifonames[i],"%s",ifo);
        }
        i++;
    }
    ifonames[numifos]=NULL;
    i=0;
    while (ifonames[i]!= NULL){
        
        if (!strcmp(ifonames[i],"H1") || !strcmp(ifonames[i],"L1")){
            /* Check either PSD file or fakePSD are given */
            if(!(ligoFakePsd || ligoPsdFileName )){
                fprintf( stderr,
                "Must provide PSD file or the name of analytic PSD for LIGO if --snr-distr is given and H1 or L1 are in --ifos. \n" );
                exit( 1 );
            }
            /* Check we didn't give both fake PSD and filename*/
            if ( ligoFakePsd && ligoPsdFileName ){
                fprintf( stderr,"Must provide only one between --ligo-psd and --ligo-fake-psd \n" );
                exit( 1 );
            }
            /* Check flow for SNR calculation was given */
            if (ligoStartFreq < 0) {
                fprintf( stderr, "Must specify --ligo-start-freq together with --ligo-psd.\n");
                exit( 1 );
            }
        }
        else if (!strcmp(ifonames[i],"V1")){
            /* Check either PSD file or fakePSD are given */
            if(!(virgoFakePsd || virgoPsdFileName )){
                fprintf( stderr,
                "Must provide PSD file or the name of analytic PSD for Virgo if --snr-distr is given and V1 is in --ifos. \n" );
                exit( 1 );
            }
            /* Check we didn't give both fake PSD and filename*/
            if ( virgoFakePsd && virgoPsdFileName ){
                fprintf( stderr,"Must provide only one between --virgo-psd and --virgo-fake-psd \n" );
                exit( 1 );
            }
            /* Check flow for SNR calculation was given */
            if (virgoStartFreq < 0) {
              fprintf( stderr,
                "Must specify --virgo-start-freq with --virgo-psd.\n"
                );
              exit( 1 );
            }
        }
        
        i++;
    }
    if (tmp) LALFree(tmp);
    if (ifo) LALFree(ifo);
	if ( options.maxsnr <= options.minsnr )
    {
      fprintf( stderr, "max SNR must be greater than min SNR\n");
      exit( 1 );
    }
    /*if (single_IFO_SNR_threshold<0.0)
      {
        fprintf( stderr,
            "The single IFO SNR threshold must be positive. Exiting...\n" );
        exit( 1 );
      }*/
      
       /* Check custom PSDs */
      if (ligoPsdFileName){
        ligo_interp=interpFromFile(ligoPsdFileName);
            
        
        /* We're done with the filename */
        free(ligoPsdFileName);
        }

    if (virgoPsdFileName) {
      
      virgo_interp=interpFromFile(virgoPsdFileName);
      /* We're done with the filename */
      free(virgoPsdFileName);
    }
  
  
  
    }
	switch(options.population) {
    case POPULATION_TARGETED:
    case POPULATION_ALL_SKY_SINEGAUSSIAN:
    case POPULATION_ALL_SKY_GAUSSIAN:
    case POPULATION_ALL_SKY_BTLWNB:
    case POPULATION_STRING_CUSP:
    case POPULATION_ALL_SKY_SINEGAUSSIAN_F:
    case POPULATION_ALL_SKY_DAMPEDSINUSOID:
      break;
    default:
      fprintf(stderr, "error: --population is required\n");
      exit(1);
	}

	if(!options.output) {
		int max_length = 100;	/* ARGH:  ugly */
		options.output = calloc(max_length + 1, sizeof(*options.output));
		if(options.user_tag)
			snprintf(options.output, max_length, "HL-INJECTIONS_%s-%d-%d.xml", options.user_tag, (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));
		else
			snprintf(options.output, max_length, "HL-INJECTIONS-%d-%d.xml", (int) (options.gps_start_time / LAL_INT8_C(1000000000)), (int) ((options.gps_end_time - options.gps_start_time) / LAL_INT8_C(1000000000)));
	}

	return options;
}


/* 
 * ============================================================================
 *
 *                                XML Handling
 *
 * ============================================================================
 */


#define DEFINELIGOLWTABLEAPPEND(funcroot, rowtype) \
static rowtype *XLAL ## funcroot ## Append(rowtype *head, rowtype *row) \
{ \
	rowtype *tail; \
	if(!head) \
		return row; \
	for(tail = head; tail->next; tail = tail->next); \
	tail->next = row; \
	return head; \
}


DEFINELIGOLWTABLEAPPEND(ProcessTable, ProcessTable)
DEFINELIGOLWTABLEAPPEND(ProcessParamsTable, ProcessParamsTable)
DEFINELIGOLWTABLEAPPEND(TimeSlideTable, TimeSlide)
DEFINELIGOLWTABLEAPPEND(SearchSummaryTable, SearchSummaryTable)
DEFINELIGOLWTABLEAPPEND(SimBurstTable, SimBurst)


static int load_tisl_file_and_merge(const char *filename, ProcessTable **process_table_head, ProcessParamsTable **process_params_table_head, TimeSlide **time_slide_table_head, SearchSummaryTable **search_summary_table_head, SimBurst **sim_burst_table_head)
{
	ProcessTable *tisl_process_table_head, *process_row;
	ProcessParamsTable *tisl_process_params_table_head, *process_params_row;
	SearchSummaryTable *tisl_search_summary_table_head;
	TimeSlide *tisl_time_slide_table_head, *time_slide_row;
	SearchSummaryTable *search_summary_row;
	SimBurst *tisl_sim_burst_table_head, *sim_burst_row;
	long process_id;

	/* load the tables from the time slide document */

	tisl_process_table_head = XLALProcessTableFromLIGOLw(filename);
	if(!tisl_process_table_head)
		return -1;
	tisl_process_params_table_head = XLALProcessParamsTableFromLIGOLw(filename);
	if(!tisl_process_params_table_head)
		return -1;
	tisl_time_slide_table_head = XLALTimeSlideTableFromLIGOLw(filename);
	if(!tisl_time_slide_table_head)
		return -1;
	if(XLALLIGOLwHasTable(filename, "search_summary")) {
		tisl_search_summary_table_head = XLALSearchSummaryTableFromLIGOLw(filename);
		if(!tisl_search_summary_table_head)
			return -1;
	} else
		tisl_search_summary_table_head = NULL;
	if(XLALLIGOLwHasTable(filename, "sim_burst")) {
		tisl_sim_burst_table_head = XLALSimBurstTableFromLIGOLw(filename, NULL, NULL);
		if(!tisl_sim_burst_table_head)
			return -1;
	}
		tisl_sim_burst_table_head = NULL;

	/* check for more than one time slide in the document */

	for(time_slide_row = tisl_time_slide_table_head->next; time_slide_row; time_slide_row = time_slide_row->next)
		if(time_slide_row->time_slide_id != tisl_time_slide_table_head->time_slide_id) {
			fprintf(stderr, "error: time slide file \"%s\" contains more than 1 offset vector", filename);
			return -1;
		}

	/* find the next available process ID and reassign binj's rows to
	 * avoid collisions */

	process_id = XLALProcessTableGetNextID(tisl_process_table_head);
	for(process_row = *process_table_head; process_row; process_row = process_row->next)
		process_row->process_id = process_id;
	for(process_params_row = *process_params_table_head; process_params_row; process_params_row = process_params_row->next)
		process_params_row->process_id = process_id;
	for(search_summary_row = *search_summary_table_head; search_summary_row; search_summary_row = search_summary_row->next)
		search_summary_row->process_id = process_id;
	for(sim_burst_row = *sim_burst_table_head; sim_burst_row; sim_burst_row = sim_burst_row->next)
		sim_burst_row->process_id = process_id;

	/* append our rows to the process and process params tables */

	*process_table_head = XLALProcessTableAppend(tisl_process_table_head, *process_table_head);
	*process_params_table_head = XLALProcessParamsTableAppend(tisl_process_params_table_head, *process_params_table_head);
	*time_slide_table_head = XLALTimeSlideTableAppend(tisl_time_slide_table_head, *time_slide_table_head);
	*search_summary_table_head = XLALSearchSummaryTableAppend(tisl_search_summary_table_head, *search_summary_table_head);
	*sim_burst_table_head = XLALSimBurstTableAppend(tisl_sim_burst_table_head, *sim_burst_table_head);

	/* done */

	return 0;
}


static int qsort_strcmp(char **a, char **b)
{
	return strcmp(*a, *b);
}


static int set_instruments(ProcessTable *process, TimeSlide *time_slide_table_head)
{
	char *ifos = process->ifos;
	char **time_slide_instruments;
	int n_instruments, n;
	TimeSlide *time_slide;

	/* count the rows in the time_slide table */
	for(n_instruments = 0, time_slide = time_slide_table_head; time_slide; n_instruments++, time_slide = time_slide->next);

	/* allocate an array of pointers to the instrument values and sort
	 * in alphabetical order */
	time_slide_instruments = malloc(n_instruments * sizeof(*time_slide_instruments));
	if(!time_slide_instruments)
		return -1;
	for(n = 0, time_slide = time_slide_table_head; time_slide; n++, time_slide = time_slide->next)
		time_slide_instruments[n] = time_slide->instrument;
	qsort(time_slide_instruments, n_instruments, sizeof(*time_slide_instruments), (int (*)(const void *, const void *)) qsort_strcmp);

	/* merge into a comma-delimited string */
	for(n = 0; n < n_instruments; n++) {
		int copied;
		copied = snprintf(ifos, sizeof(process->ifos) - (ifos - process->ifos), "%s%s", time_slide_instruments[n], n < n_instruments - 1 ? "," : "");
		if(copied < 2 + (n < n_instruments - 1 ? 1 : 0)) {
			fprintf(stderr, "error:  too many instruments for process table's ifos column\n");
			free(time_slide_instruments);
			return -1;
		}
		ifos += copied;
	}
	free(time_slide_instruments);
	return 0;
}


/* 
 * ============================================================================
 *
 *                                 Sequences
 *
 * ============================================================================
 */


#if 0
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
#endif


/*
 * Repeating geometric sequence.
 */

/*
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
}*/


#if 0
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
#endif


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

static double draw_uniform(gsl_rng *rng, double a, double b)
{
return a+ (b-a)*gsl_rng_uniform(rng);	
}

static double draw_gaussian(gsl_rng *rng, double mu, double sigma)
{
	return (mu + gsl_ran_gaussian(rng,sigma));
}

static double draw_volume(gsl_rng *rng,double min,double max){
    /* Uniform in volume for hrss (~1/Distance) means p(hrss)~hrss^-4 
     * which implies p(1/hrss^3) = constant 
     */
    REAL8 proposed=0.0;
    
    proposed=1.0/(max*max*max)+(1.0/(min*min*min)- 1.0/(max*max*max))*gsl_rng_uniform(rng);
    proposed=1.0/cbrt(proposed);
    return proposed;
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


static SimBurst *random_directed_btlwnb(double ra, double dec, double psi, double minf, double maxf, double minband, double maxband, double mindur, double maxdur, double minEoverr2, double maxEoverr2, struct options *options, gsl_rng *rng)
{
	REAL8TimeSeries *hplus, *hcross;
	SimBurst *sim_burst = XLALCreateSimBurst();

  (void) minf;
  (void) maxf;
  (void) minEoverr2;
  (void) maxEoverr2;
  
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
  switch (options->f_distr){
		case FIXED:
			sim_burst->frequency = options->f;
			break;
		case UNIFORM:
			sim_burst->frequency=draw_uniform(rng,options->minf,options->maxf);
			break;
		case GAUSSIAN:
			sim_burst->frequency=draw_gaussian(rng,options->f,options->f_stdev);
			break;
		case UNIFORM_LOG:
			sim_burst->frequency=ran_flat_log(rng, options->minf, options->maxf);
			break;
		default:
			fprintf(stderr,"unknown distribution of frequency. Known values are fixed, uniform, gaussian, log.\n");
			exit(1);
	}
	//sim_burst->frequency = ran_flat_log_discrete(rng, minf, maxf, pow(maxf / minf, 1.0 / 3.0));

	/* duration and bandwidth.  keep picking until a valid pair is
 	 * obtained (i.e. their product is >= 2 / \pi) */
	do {
		sim_burst->duration = draw_uniform(rng, mindur, maxdur);
		sim_burst->bandwidth = draw_uniform(rng, minband, maxband);
	} while(sim_burst->bandwidth * sim_burst->duration < LAL_2_PI);

	/* energy -- Since we want to obtain a given distribution on hrss, 
   * instead of drawing an energy and calculate the resulting hrss (as 
   * in binj.c) we draw and hrss and rescale the energy.
   * Start by fixing the energy to an arbitrary value.
   */
	sim_burst->egw_over_rsquared = 1e-10;
  
  REAL8 rand_hrss;
  switch (options->hrss_distr){
			case FIXED:
				rand_hrss = options->hrss;
				break;
			case UNIFORM:
				rand_hrss=draw_uniform(rng,options->minhrss,options->maxhrss);
				break;
			case VOLUME:
				rand_hrss=draw_volume(rng,options->minhrss,options->maxhrss);
				break;
			case GAUSSIAN:
        rand_hrss=draw_gaussian(rng,options->hrss,options->hrss_stdev);
				break;
			case UNIFORM_LOG:
				rand_hrss=ran_flat_log(rng, options->minhrss, options->maxhrss);
				break;
			default:
				fprintf(stderr,"unknown distribution of hrss. Known values are fixed, uniform, gaussian, log, and volume.\n");
				exit(1);	
				
		}
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

	/* populate the hrss column for convenience later */
	/* FIXME:  sample rate is hard-coded to 8192 Hz, which is what the
	 * excess power pipeline's .ini file is configured for in CVS, but
	 * because it can be easily changed this is not good */

	XLALGenerateSimBurst(&hplus, &hcross, sim_burst, 1.0 / 8192);
	if(!hplus || !hcross) {
		XLALDestroyREAL8TimeSeries(hplus);
		XLALDestroyREAL8TimeSeries(hcross);
		XLALDestroySimBurst(sim_burst);
		return NULL;
	}
  
  /* Now fix hrss to the random value and scale the energy */
	sim_burst->hrss = rand_hrss;
  sim_burst->egw_over_rsquared = 1e-10*pow(rand_hrss/XLALMeasureHrss(hplus, hcross),2.);
  
  /* clean up */
	XLALDestroyREAL8TimeSeries(hplus);
	XLALDestroyREAL8TimeSeries(hcross);

	/* done */

	return sim_burst;
}


static SimBurst *random_all_sky_btlwnb(double minf, double maxf, double minband, double maxband, double mindur, double maxdur, double minEoverr2, double maxEoverr2, struct options *options, gsl_rng *rng)
{
	double ra, dec, psi;

	random_location_and_polarization(&ra, &dec, &psi, rng);

	return random_directed_btlwnb(ra, dec, psi, minf, maxf, minband, maxband, mindur, maxdur, minEoverr2, maxEoverr2, options, rng);
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
 


static double duration_from_q_and_f(double Q, double f)
{
	// compute duration from Q and frequency 
	return Q / (sqrt(2.0) * LAL_PI * f);
}*/


/*
 * pick a sine-Gaussian
 */


static SimBurst *random_all_sky_sineGaussian( gsl_rng *rng, struct options *options,REAL8 tinj)
{
	SimBurst *sim_burst = XLALCreateSimBurst();
  REAL8 this_snr;
  int idx=0;
  int max_i=1000;
  if(!sim_burst)
      return NULL;
  XLALINT8NSToGPS(&(sim_burst)->time_geocent_gps, tinj);

  if (options->population==POPULATION_ALL_SKY_SINEGAUSSIAN)
    strcpy(sim_burst->waveform, "SineGaussian");
  else if (options->population==POPULATION_ALL_SKY_SINEGAUSSIAN_F)
    strcpy(sim_burst->waveform, "SineGaussianF");
  else if (options->population==POPULATION_ALL_SKY_GAUSSIAN)
    strcpy(sim_burst->waveform, "Gaussian");
  else if (options->population==POPULATION_ALL_SKY_DAMPEDSINUSOID)
    strcpy(sim_burst->waveform, "DampedSinusoid");
  else{
    fprintf(stderr,"Unrecognized population %d. Exiting\n",options->population);
    exit(1);
  }
  
  while(1){
          
    /* sky location and wave frame orientation */

    random_location_and_polarization(&sim_burst->ra, &sim_burst->dec, &sim_burst->psi, rng);
    
    /* hard-code for linearly polarized waveforms in the x
     * polarization.  induces LAL's sine-Gaussian generator to produce
     * linearly polarized sine-Gaussians (+ would be a cosine
     * Gaussian). */
    /* pick a polar ellipse angle  and eccentricity */
    switch (options->polee_distr){
      case FIXED:
        sim_burst->pol_ellipse_e = options->pol_ellipse_e;
        break;
      case UNIFORM:
        sim_burst->pol_ellipse_e =draw_uniform(rng,options->minpolee,options->maxpolee);
        break;
      default:
        // default are linealry polarized signals
        sim_burst->pol_ellipse_e = 1.0;
    }
    switch (options->polea_distr){
      case FIXED:
        sim_burst->pol_ellipse_angle = options->pol_ellipse_angle;
        break;
      case UNIFORM:
        sim_burst->pol_ellipse_angle =draw_uniform(rng,options->minpolea,options->maxpolea);
        break;
      default:
        // default are linelary polarized signals
        sim_burst->pol_ellipse_angle = LAL_PI_2;
    }
    

    /* q and centre frequency.  three steps between minf and maxf */
    switch (options->q_distr){
      case FIXED:
        sim_burst->q = options->q;
        break;
      case UNIFORM:
        sim_burst->q=draw_uniform(rng,options->minq,options->maxq);
        break;
      case GAUSSIAN:
        do{
        sim_burst->q=draw_gaussian(rng,options->q,options->q_stdev);
        }while(sim_burst->q<=q_min);
        break;
      case UNIFORM_LOG:
        sim_burst->q=ran_flat_log(rng, options->minq, options->maxq);
        break;
      default:
        fprintf(stderr,"unknown distribution of q. Known values are fixed, uniform, gaussian, log.\n");
        exit(1);
    }
      switch (options->f_distr){
      case FIXED:
        sim_burst->frequency = options->f;
        break;
      case UNIFORM:
        sim_burst->frequency=draw_uniform(rng,options->minf,options->maxf);
        break;
      case GAUSSIAN:
        sim_burst->frequency=draw_gaussian(rng,options->f,options->f_stdev);
        break;
      case UNIFORM_LOG:
        sim_burst->frequency=ran_flat_log(rng, options->minf, options->maxf);
        break;
      default:
        fprintf(stderr,"unknown distribution of frequency. Known values are fixed, uniform, gaussian, log.\n");
        exit(1);
    }
    /* duration */
    switch (options->tau_distr){
      case FIXED:
        sim_burst->duration = options->tau;
        break;
      case UNIFORM:
        sim_burst->duration=draw_uniform(rng,options->minduration,options->maxduration);
        break;
      case GAUSSIAN:
        sim_burst->duration=draw_gaussian(rng,options->tau,options->tau_stdev);
        break;
      case UNIFORM_LOG:
        sim_burst->duration=ran_flat_log(rng, options->minduration, options->maxduration);
        break;
      default:
        fprintf(stderr,"unknown distribution of duration tau. Known values are fixed, uniform, gaussian, log.\n");
        exit(1);
    }
    /* hrss */
    if (options->hrss_distr< NUM_ELEMENTS){
      /* Inject using distribution on hrss */
      switch (options->hrss_distr){
        case FIXED:
          sim_burst->hrss = options->hrss;
          break;
        case UNIFORM:
          sim_burst->hrss=draw_uniform(rng,options->minhrss,options->maxhrss);
          break;
        case VOLUME:
          sim_burst->hrss=draw_volume(rng,options->minhrss,options->maxhrss);
          break;
        case GAUSSIAN:
          sim_burst->hrss=draw_gaussian(rng,options->hrss,options->hrss_stdev);
          break;
        case UNIFORM_LOG:
          sim_burst->hrss=ran_flat_log(rng, options->minhrss, options->maxhrss);
          break;
        default:
          fprintf(stderr,"unknown distribution of hrss. Known values are fixed, uniform, gaussian, log, and volume.\n");
          exit(1);	
      }
    }
    if (options->minsnr>0.&& options->maxsnr>0.){
      
      /* Check if SNR is inside range, otherwise redraw al parameters  */
      char *ifo;
      REAL8 *start_freqs;
      REAL8FrequencySeries **psds;
      int i=1;
      UINT4 ui=0;
      /*reset counter */
      ifo=ifonames[0];
      i=0;
      /* Create variables for PSDs and starting frequencies */
      start_freqs = (REAL8 *) LALCalloc(numifos+1, sizeof(REAL8));
      psds        = (REAL8FrequencySeries **) LALCalloc(numifos+1, sizeof(REAL8FrequencySeries *));
      REAL8    srate=8192.0;
      REAL8 segment=8.0;
      size_t seglen=(size_t) segment*srate;
      LIGOTimeGPS ttime;
      memcpy(&ttime,&(sim_burst->time_geocent_gps),sizeof(LIGOTimeGPS));
      
      /* Fill psds and start_freqs */
      /* If the user did not provide files for the PSDs, use XLALSimNoisePSD to fill in ligoPsd and virgoPsd */
      while(ifo !=NULL){
          if(!strcmp("V1",ifo)){
                  start_freqs[i]=virgoStartFreq;
                  if (!virgoPsd){
                              
                      virgoPsd=XLALCreateREAL8FrequencySeries("VPSD",&ttime , 0, 1.0/segment, &lalHertzUnit, seglen/2+1);
                      if (!virgo_interp)
                          get_FakePsdFromString(virgoPsd,virgoFakePsd, virgoStartFreq);
                      else{
                          for (ui=0;ui<virgoPsd->data->length;ui++){
                              virgoPsd->data->data[ui]=interpolate(virgo_interp,ui/segment);
                      
                          }
                      }
                  }
                  if (!virgoPsd) fprintf(stderr,"Failed to produce Virgo PSD series. Exiting...\n");
                  psds[i]=virgoPsd;
          }
          else if (!strcmp("L1",ifo) || !strcmp("H1",ifo)){
        start_freqs[i]=ligoStartFreq;
        if(!ligoPsd){
            ligoPsd=XLALCreateREAL8FrequencySeries("LPSD", &ttime, 0, 1.0/segment, &lalHertzUnit, seglen/2+1);
                  if (! ligo_interp)
                      get_FakePsdFromString(ligoPsd,ligoFakePsd,ligoStartFreq);
                  else{
                      for (ui=0;ui<ligoPsd->data->length;ui++){
                      ligoPsd->data->data[ui]=interpolate(ligo_interp,ui/segment);
                  
                      }
                  }			}   
        if (!ligoPsd) fprintf(stderr,"Failed to produce LIGO PSD series. Exiting...\n");   
        psds[i]=ligoPsd;
          }
          else{
        fprintf(stderr,"Unknown IFO. Allowed IFOs are H1,L1 and V1. Exiting...\n");
        exit(-1);
        }
          i++;
          ifo=ifonames[i];
          }
    
      options->ifonames=ifonames;
      options->nIFO=i+1;
      
      /* If 1 detector is used, turn the single IFO snr check off. */ 
      if (numifos<2){
        fprintf(stdout,"Warning: You are using less than 2 IFOs. Disabling the single IFO SNR threshold check...\n");
        single_IFO_SNR_threshold=0.0;
      }
      
      /* This function takes care of drawing a proposed SNR and set the distance accordingly  */
      this_snr=calculate_NetSNR(sim_burst,ifonames, psds, start_freqs, options, rng);
      
      /* Clean  */
      if (psds) LALFree(psds);
      if (start_freqs) LALFree(start_freqs);
      /* Done */
      
      if (options->minsnr <=this_snr && options->maxsnr>=this_snr){
        break;
      }
      idx+=1;
      if (idx>=max_i){
        fprintf(stdout,"Could not have SNR in the desired range with less than 1000 trials. Keeping last set of parameters.\n");
        break;
      }
    }
    else
      /* No SNR check required, just break*/
      break;
  }
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


static void write_xml(const char *filename, const ProcessTable *process_table_head, const ProcessParamsTable *process_params_table_head, const SearchSummaryTable *search_summary_table_head, const TimeSlide *time_slide_table_head, const SimBurst *sim_burst)
{
	LIGOLwXMLStream *xml;

	xml = XLALOpenLIGOLwXMLFile(filename);

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
	if(XLALWriteLIGOLwXMLSearchSummaryTable(xml, search_summary_table_head)) {
		/* error occured.  ?? do anything else ?? */
		exit(1);
	}

	/* time slide table */
	if(XLALWriteLIGOLwXMLTimeSlideTable(xml, time_slide_table_head)) {
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
	INT8 tinj, jitter;
	gsl_rng *rng;
	ProcessTable *process_table_head = NULL, *process;
	ProcessParamsTable *process_params_table_head = NULL;
	SearchSummaryTable *search_summary_table_head = NULL, *search_summary;
	TimeSlide *time_slide_table_head = NULL;
	SimBurst *sim_burst_table_head = NULL;
	SimBurst **sim_burst = &sim_burst_table_head;


	/*
	 * Initialize debug handler
	 */


	lal_errhandler = LAL_ERR_EXIT;


	/*
	 * Process table
	 */


	process_table_head = process = XLALCreateProcessTableRow();
	if(XLALPopulateProcessTable(process, PROGRAM_NAME, lalAppsVCSIdentId, lalAppsVCSIdentStatus, lalAppsVCSIdentDate, 0))
		exit(1);
	XLALGPSTimeNow(&process->start_time);


	/*
	 * Command line and process params table.
	 */


	options = parse_command_line(&argc, &argv, process, &process_params_table_head);
	if(options.user_tag)
		snprintf(process->comment, sizeof(process->comment), "%s", options.user_tag);


	/*
	 * Search summary table
	 */


	search_summary_table_head = search_summary = XLALCreateSearchSummaryTableRow(process);
	if(options.user_tag)
		snprintf(search_summary->comment, sizeof(search_summary->comment), "%s", options.user_tag);
	search_summary->nnodes = 1;
	search_summary->out_start_time = *XLALINT8NSToGPS(&search_summary->in_start_time, options.gps_start_time);
	search_summary->out_end_time = *XLALINT8NSToGPS(&search_summary->in_end_time, options.gps_end_time);


	/*
	 * Initialize random number generator
	 */


	rng = gsl_rng_alloc(gsl_rng_mt19937);
	if(options.seed)
		gsl_rng_set(rng, options.seed);


	/*
	 * Main loop
	 */

	//DistrOpts *distropt=NULL;
	//distropt=XLALMalloc(sizeof(DistrOpts));
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
			*sim_burst = random_directed_btlwnb(options.ra, options.dec, gsl_ran_flat(rng, 0, LAL_TWOPI), options.minf, options.maxf, options.minbandwidth, options.maxbandwidth, options.minduration, options.maxduration, options.minEoverr2, options.maxEoverr2,&options, rng);
			break;

		case POPULATION_ALL_SKY_SINEGAUSSIAN:
		case POPULATION_ALL_SKY_SINEGAUSSIAN_F:
    case POPULATION_ALL_SKY_GAUSSIAN:
    case POPULATION_ALL_SKY_DAMPEDSINUSOID:
			*sim_burst = random_all_sky_sineGaussian(rng, &options,tinj);
			break;

		case POPULATION_ALL_SKY_BTLWNB:
			*sim_burst = random_all_sky_btlwnb(options.minf, options.maxf, options.minbandwidth, options.maxbandwidth, options.minduration, options.maxduration, options.minEoverr2, options.maxEoverr2, &options, rng);
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

		/* Add "jitter" to the injection if user requests it */
		if(options.jitter > 0) {
			jitter = gsl_ran_flat(rng, -options.jitter/2, options.jitter/2)*1e9;
		} else {
			jitter = 0;
		}

		XLALINT8NSToGPS(&(*sim_burst)->time_geocent_gps, tinj + jitter);

		/*
		 * Peak time at geocentre in GMST radians
		 */

		(*sim_burst)->time_geocent_gmst = XLALGreenwichMeanSiderealTime(&(*sim_burst)->time_geocent_gps);

		/*
		 * Move to next injection
		 */

		sim_burst = &(*sim_burst)->next;
	}
   // XLALFree(distropt);

	XLALPrintInfo("%s: ", argv[0]);
	XLALPrintProgressBar(1.0);
	XLALPrintInfo(" complete\n");

	/* load time slide document and merge our table rows with its */

	if(load_tisl_file_and_merge(options.time_slide_file, &process_table_head, &process_params_table_head, &time_slide_table_head, &search_summary_table_head, &sim_burst_table_head))
		exit(1);
	if(set_instruments(process, time_slide_table_head))
		exit(1);
	snprintf(search_summary->ifos, sizeof(search_summary->ifos), "%s", process->ifos);

	/* output */

	XLALGPSTimeNow(&process->end_time);
	search_summary->nevents = XLALSimBurstAssignIDs(sim_burst_table_head, process->process_id, time_slide_table_head->time_slide_id, 0);
	write_xml(options.output, process_table_head, process_params_table_head, search_summary_table_head, time_slide_table_head, sim_burst_table_head);
	//if (options.write_mdc)
	//	write_mdc(&sim_burst_table_head, time_slide_table_head, &options);
	/* done */

	gsl_rng_free(rng);
	XLALDestroyProcessTable(process_table_head);
	XLALDestroyProcessParamsTable(process_params_table_head);
	XLALDestroyTimeSlideTable(time_slide_table_head);
	XLALDestroySearchSummaryTable(search_summary_table_head);
	XLALDestroySimBurstTable(sim_burst_table_head);
	exit(0);
}


static void get_FakePsdFromString(REAL8FrequencySeries* PsdFreqSeries,char* FakePsdName, REAL8 StartFreq){
    
    /* Call XLALSimNoisePSD to fill the REAL8FrequencySeries PsdFreqSeries (must been already allocated by callers). FakePsdName contains the label of the fake PSD */
       if (!strcmp("LALSimAdVirgo",FakePsdName)){
                XLALSimNoisePSD(PsdFreqSeries,StartFreq,XLALSimNoisePSDAdvVirgo);
        }
        else if (!strcmp("LALSimVirgo",FakePsdName)){
            XLALSimNoisePSD(PsdFreqSeries,StartFreq,XLALSimNoisePSDVirgo);
        }
        else if(!strcmp("LALSimAdLIGO",FakePsdName)){
            XLALSimNoisePSD(PsdFreqSeries,StartFreq,XLALSimNoisePSDaLIGOZeroDetHighPower);
        }
        else if(!strcmp("LALSimLIGO",FakePsdName)){
            XLALSimNoisePSD(PsdFreqSeries,StartFreq,XLALSimNoisePSDiLIGOSRD);
        }
        else{
            fprintf(stderr,"Unknown fake PSD %s. Known types are LALSimLIGO, LALSimAdLIGO, LALSimAdVirgo, LALSimVirgo. Exiting...\n",FakePsdName);
            exit(1);
            }
}

static REAL8  calculate_NetSNR(SimBurst *inj,char ** IFOnames, REAL8FrequencySeries **psds,REAL8 *start_freqs, struct options *options, gsl_rng *rng )
{
    
    //REAL8 local_min=0.0;
    (void) rng;
    (void) options;
    REAL8 net_snr=0.0;
    REAL8 * SNRs=NULL;
    UINT4 j=0;
    UINT4 num_ifos=0;
       
    if (IFOnames ==NULL){
        fprintf(stderr,"Calculate_NetSNR() called with IFOnames=NULL. Exiting...\n");
        exit(1);
        }
    char * ifo=IFOnames[0];
    
    /* Get the number of IFOs from IFOnames*/
    while(ifo !=NULL){
        num_ifos++;
        ifo=IFOnames[num_ifos];
    }
            
    SNRs=calloc(num_ifos+1 ,sizeof(REAL8));
    /* Calculate the single IFO and network SNR for the dummy distance of 100Mpc */
    for (j=0;j<num_ifos;j++){
        SNRs[j]=calculate_SineGaussian_snr(inj,IFOnames[j],psds[j],start_freqs[j]);
        net_snr+=SNRs[j]*SNRs[j];
    }
    net_snr=sqrt(net_snr);
    
    if (SNRs) free(SNRs);
    
    return net_snr;
}

static REAL8 calculate_SineGaussian_snr(SimBurst *inj, char *IFOname, REAL8FrequencySeries *psd, REAL8 start_freq)
{
    /* Calculate and return the single IFO SNR 
     * 
     * Required options:
     * 
     * inj:     SimInspiralTable entry for which the SNR has to be calculated
     * IFOname: The canonical name (e.g. H1, L1, V1) name of the IFO for which the SNR must be calculated
     * PSD:     PSD curve to be used for the overlap integrap
     * start_freq: lower cutoff of the overlap integral
     * 
     * */   
     
    UINT4 j=0;
    /* Fill detector site info */
    LALDetector*  detector=NULL;
    detector=calloc(1,sizeof(LALDetector));
     if(!strcmp(IFOname,"H1")) 			
        memcpy(detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
    if(!strcmp(IFOname,"H2")) 
        memcpy(detector,&lalCachedDetectors[LALDetectorIndexLHODIFF],sizeof(LALDetector));
    if(!strcmp(IFOname,"LLO")||!strcmp(IFOname,"L1")) 
        memcpy(detector,&lalCachedDetectors[LALDetectorIndexLLODIFF],sizeof(LALDetector));
    if(!strcmp(IFOname,"V1")||!strcmp(IFOname,"VIRGO")) 
        memcpy(detector,&lalCachedDetectors[LALDetectorIndexVIRGODIFF],sizeof(LALDetector));
    

    
    REAL8 Q =0.0;
    REAL8 hrss=0.0;
    REAL8 centre_frequency= 0.0;
    REAL8 polar_angle=0.0;
    REAL8 eccentricity=0.0;
    REAL8 latitude=0.0;
    REAL8 polarization=0.0;
    REAL8 injtime=0.0;
    REAL8 longitude;
    LALSimulationDomain modelDomain=LAL_SIM_DOMAIN_FREQUENCY;
    REAL8 f_max;
    Q=inj->q;
    centre_frequency=inj->frequency;
    hrss=inj->hrss;
    polarization=inj->psi;
    polar_angle=inj->pol_ellipse_angle;
    eccentricity=inj->pol_ellipse_e; 
    injtime=inj->time_geocent_gps.gpsSeconds + 1e-9*inj->time_geocent_gps.gpsNanoSeconds;
    latitude=inj->dec;
    longitude=inj->ra;

    const CHAR *WF=inj->waveform;
    BurstApproximant approx=XLALGetBurstApproximantFromString(WF);

    if (XLALSimBurstImplementedFDApproximants(approx))
      modelDomain=LAL_SIM_DOMAIN_FREQUENCY;
    else if (XLALSimBurstImplementedTDApproximants(approx))
      modelDomain=LAL_SIM_DOMAIN_TIME;

    LIGOTimeGPS epoch;  
    
    /* Hardcoded values of srate and segment length. If changed here they must also be changed in inspinj.c */
    REAL8 srate=8192.0;
    REAL8 segment=8.0;
    //REAL8 start=injtime-segment/2.;
    memcpy(&epoch,&inj->time_geocent_gps,sizeof(LIGOTimeGPS));
    XLALGPSAdd(&epoch, -segment/2.0);
    
    f_max=(srate/2.0-(1.0/segment));
    size_t seglen=(size_t) segment*srate;
    REAL8 deltaF=1.0/segment ;
    REAL8 deltaT=1.0/srate;

     
    /* Frequency domain h+ and hx. They are going to be filled either by a FD WF or by the FFT of a TD WF*/ 
    COMPLEX16FrequencySeries *freqHplus;
    COMPLEX16FrequencySeries* freqHcross;
    freqHplus=  XLALCreateCOMPLEX16FrequencySeries("fhplus",
										&epoch,
										0.0,
										deltaF,
										&lalDimensionlessUnit,
										seglen/2+1);
                                        
    freqHcross=XLALCreateCOMPLEX16FrequencySeries("fhcross",
										&epoch,
										0.0,
										deltaF,
										&lalDimensionlessUnit,
										seglen/2+1);
     
    /* If the approximant is on the FD call XLALSimBurstSineGaussianF */
    if(modelDomain == LAL_SIM_DOMAIN_FREQUENCY) {
        COMPLEX16FrequencySeries *hptilde=NULL;
        COMPLEX16FrequencySeries *hctilde=NULL;
      
        XLALSimBurstChooseFDWaveform(&hptilde, &hctilde, deltaF,deltaT,centre_frequency,Q,0.0,0.0,f_max,hrss,polar_angle,eccentricity,NULL,approx);
        
        COMPLEX16 *dataPtr = hptilde->data->data;
        for (j=0; j<(UINT4) freqHplus->data->length; ++j) {
            if(j <= hptilde->data->length){
                freqHplus->data->data[j] = dataPtr[j];
            }else{
                  freqHplus->data->data[j]=0.0+0.0*I;
            }
        }
         dataPtr = hctilde->data->data;
         for (j=0; j<(UINT4) freqHplus->data->length; ++j) {
            if(j <= hctilde->data->length){
                freqHcross->data->data[j] = dataPtr[j];
            }else{
                  freqHcross->data->data[j]=0.0+0.0*I;
            }
        }
    /* Clean */    
    if(hptilde) XLALDestroyCOMPLEX16FrequencySeries(hptilde);
    if(hctilde) XLALDestroyCOMPLEX16FrequencySeries(hctilde);

    }
    else{
        
        REAL8FFTPlan *timeToFreqFFTPlan = XLALCreateForwardREAL8FFTPlan((UINT4) seglen, 0 );
        REAL8TimeSeries *hplus=NULL; 
        REAL8TimeSeries *hcross=NULL; 
        REAL8TimeSeries *timeHplus=NULL;
        REAL8TimeSeries *timeHcross=NULL;
      
        timeHcross=XLALCreateREAL8TimeSeries("timeModelhCross",
																	&epoch,
																	0.0,
																	deltaT,
																	&lalStrainUnit,
																	seglen);
        timeHplus=XLALCreateREAL8TimeSeries("timeModelhplus",
																	&epoch,
																	0.0,
																	deltaT,
																	&lalStrainUnit,
																	seglen);
        XLALGenerateSimBurst(&hplus,&hcross,inj,deltaT);
        memset(timeHplus->data->data, 0, sizeof (REAL8)*timeHplus->data->length);
        memset(timeHcross->data->data, 0, sizeof (REAL8)*timeHcross->data->length);
        XLALGPSSetREAL8(&hcross->epoch, XLALGPSGetREAL8(&inj->time_geocent_gps));
		    XLALGPSSetREAL8(&hplus->epoch, XLALGPSGetREAL8(&inj->time_geocent_gps));
        XLALSimAddInjectionREAL8TimeSeries(timeHplus, hplus, NULL);
        XLALSimAddInjectionREAL8TimeSeries(timeHcross, hcross, NULL);
        
        for (j=0; j<(UINT4) freqHplus->data->length; ++j) {
                freqHplus->data->data[j]=0.0; 
                freqHcross->data->data[j]=0.0;
            }
        REAL8 padding=0.4;
        REAL8Window *window;
        window =XLALCreateTukeyREAL8Window(seglen,(REAL8)2.0*padding*srate/(REAL8)seglen);
        if(!window) {
          XLAL_ERROR(XLAL_EFUNC);
        }
       
        for(j = 0; j < window->data->length; j++) {
          timeHplus ->data->data[j] *= window->data->data[j];
          timeHcross->data->data[j] *= window->data->data[j];
        }
        XLALDestroyREAL8Window(window);   
         
        /* FFT into freqHplus and freqHcross */
        XLALREAL8TimeFreqFFT(freqHplus,timeHplus,timeToFreqFFTPlan);
        XLALREAL8TimeFreqFFT(freqHcross,timeHcross,timeToFreqFFTPlan);

        /* Clean... */
        if ( hplus ) XLALDestroyREAL8TimeSeries(hplus);
        if ( hcross ) XLALDestroyREAL8TimeSeries(hcross);
        if ( timeHplus ) XLALDestroyREAL8TimeSeries(timeHplus);
        if ( timeHcross ) XLALDestroyREAL8TimeSeries(timeHcross);
        if (timeToFreqFFTPlan) LALFree(timeToFreqFFTPlan);
    }

    /* The WF has been generated and is in freqHplus/cross. Now project into the IFO frame */
    double Fplus, Fcross;
    double FplusScaled, FcrossScaled;
    double HSquared;
    double GPSdouble=injtime;
    double gmst;
    LIGOTimeGPS GPSlal;
    XLALGPSSetREAL8(&GPSlal, GPSdouble);
    gmst=XLALGreenwichMeanSiderealTime(&GPSlal);
    /* Fill Fplus and Fcross*/
    XLALComputeDetAMResponse(&Fplus, &Fcross,(const REAL4(*)[3]) detector->response,longitude, latitude, polarization, gmst);
    /* And take the distance into account */
    FplusScaled  = Fplus  ;
    FcrossScaled = Fcross ;
    REAL8 timedelay = XLALTimeDelayFromEarthCenter(detector->location,longitude, latitude, &GPSlal);
    REAL8 timeshift =  timedelay;
    REAL8  twopit    = LAL_TWOPI * timeshift;

    UINT4 lower = (UINT4)ceil(start_freq / deltaF);
    UINT4 upper = (UINT4)floor(f_max / deltaF);
    REAL8 re = cos(twopit*deltaF*lower);
    REAL8  im = -sin(twopit*deltaF*lower);

    /* Incremental values, using cos(theta) - 1 = -2*sin(theta/2)^2 */
    REAL8 dim = -sin(twopit*deltaF);
    REAL8 dre = -2.0*sin(0.5*twopit*deltaF)*sin(0.5*twopit*deltaF);
    REAL8 TwoDeltaToverN = 2.0 *deltaT / ((double) seglen);

    REAL8  plainTemplateReal,  plainTemplateImag,templateReal,templateImag;
    REAL8  newRe, newIm,temp;
    REAL8 this_snr=0.0;

    for (j=lower; j<=(UINT4) upper; ++j){
      /* derive template (involving location/orientation parameters) from given plus/cross waveforms: */
      plainTemplateReal = FplusScaled * creal(freqHplus->data->data[j])  
                          +  FcrossScaled *creal(freqHcross->data->data[j]);
      plainTemplateImag = FplusScaled * cimag(freqHplus->data->data[j])  
                          +  FcrossScaled * cimag(freqHcross->data->data[j]);

      /* do time-shifting...             */
      /* (also un-do 1/deltaT scaling): */
      templateReal = (plainTemplateReal*re - plainTemplateImag*im) / deltaT;
      templateImag = (plainTemplateReal*im + plainTemplateImag*re) / deltaT;
      HSquared  = templateReal*templateReal + templateImag*templateImag ;
      temp = ((TwoDeltaToverN * HSquared) / psd->data->data[j]);
      this_snr  += temp;
      /* Now update re and im for the next iteration. */
      newRe = re + re*dre - im*dim;
      newIm = im + re*dim + im*dre;

      re = newRe;
      im = newIm;
    }
    /* Clean */
    if (freqHcross) XLALDestroyCOMPLEX16FrequencySeries(freqHcross);
    if (freqHplus) XLALDestroyCOMPLEX16FrequencySeries(freqHplus);
    if (detector) free(detector);
    return sqrt(this_snr*2.0);
  
}

struct fvec *interpFromFile(char *filename){
	UINT4 fileLength=0;
	UINT4 i=0;
	UINT4 minLength=100; /* size of initial file buffer, and also size of increment */
	FILE *interpfile=NULL;
	struct fvec *interp=NULL;
	interp=XLALCalloc(minLength,sizeof(struct fvec)); /* Initialise array */
	if(!interp) {printf("Unable to allocate memory buffer for reading interpolation file\n");}
	fileLength=minLength;
	REAL8 f=0.0,x=0.0;
	interpfile = fopen(filename,"r");
	if (interpfile==NULL){
		printf("Unable to open file %s\n",filename);
		exit(1);
	}
	while(2==fscanf(interpfile," %lf %lf ", &f, &x )){
		interp[i].f=f; interp[i].x=x*x;
		i++;
		if(i>fileLength-1){ /* Grow the array */
			interp=XLALRealloc(interp,(fileLength+minLength)*sizeof(struct fvec));
			fileLength+=minLength;
		}
	}
	interp[i].f=0; interp[i].x=0;
	fileLength=i+1;
	interp=XLALRealloc(interp,fileLength*sizeof(struct fvec)); /* Resize array */
	fclose(interpfile);
	printf("Read %i records from %s\n",fileLength-1,filename);
	return interp;
}

REAL8 interpolate(struct fvec *fvec, REAL8 f){
	int i=0;
	REAL8 a=0.0; /* fractional distance between bins */
	REAL8 delta=0.0;
	if(f<fvec[0].f) return(0.0);
	while(fvec[i].f<f && (fvec[i].x!=0.0 )){i++;}; //&& fvec[i].f!=0.0)){i++;};
	if (fvec[i].f==0.0 && fvec[i].x==0.0) /* Frequency above moximum */
	{
		return (fvec[i-1].x);
	}
	a=(fvec[i].f-f)/(fvec[i].f-fvec[i-1].f);
	delta=fvec[i].x-fvec[i-1].x;
	return (fvec[i-1].x + delta*a);
}
