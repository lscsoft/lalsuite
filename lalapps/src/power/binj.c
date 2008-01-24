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
#include <limits.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LIGOLwXML.h>
#include <lal/Date.h>
#include <lal/TimeDelay.h>
#include <lalapps.h>
#include <processtable.h>


int snprintf(char *str, size_t size, const char *format, ...);


RCSID("$Id$");


#define KPC (1e3 * LAL_PC_SI)
#define MPC (1e6 * LAL_PC_SI)
#define GPC (1e9 * LAL_PC_SI)


#define MW_CENTER_J2000_RA_RAD 2.0318570464121519l /* = 27940.04 s = 7 h 45 m 40.04 s */
#define MW_CENTER_J2000_DEC_RAD -0.50628171572274738l /* = -29o 00' 28.1" */


#define CVS_REVISION "$Revision$"
#define CVS_SOURCE "$Source$"
#define CVS_DATE "$Date$"
#define PROGRAM_NAME "lalapps_binj"

#define TRUE      1
#define FALSE     0


/*
 * ============================================================================
 *
 *                                Command Line
 *
 * ============================================================================
 */


enum population {
	POPULATION_GALACTIC_CORE,
	POPULATION_UNIFORM_SKY,
	POPULATION_ZENITH
};


enum strain_dist {
	STRAIN_DIST_NONE,
	STRAIN_DIST_CONST_HRSS,
	STRAIN_DIST_LOG_HRSS,
	STRAIN_DIST_LOG_HRSSTAU,
	STRAIN_DIST_PRESET_HRSS
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
	unsigned long seed;
	double time_step;
	char *waveform;
	double simwaveform_duration;
	unsigned long waveform_number_min;
	unsigned long waveform_number_max;
	char *user_tag;
};


static struct options options_defaults(void)
{
	struct options defaults = {
		.gps_start_time = 0,
		.gps_end_time = 0,
		.population = POPULATION_UNIFORM_SKY,
		.freq_dist = FREQ_DIST_NONE,
		.flow = 70.0,
		.fhigh = 2000.0,
		.fratio = 0.0,
		.deltaf = 0.0,
		.quality = -1.0,
		.strain_dist = STRAIN_DIST_NONE,
		.strain_scale_min = 0.0,
		.strain_scale_max = 0.0,
		.log10_max_distance = log10(10000.0),
		.log10_min_distance = log10(100.0),
		.seed = 0,
		.time_step = 210.0 / LAL_PI,
		.waveform = "SineGaussian",
		.simwaveform_duration = 0.0,
		.waveform_number_min = 0,
		/* FIXME:  this should be ULONG_MAX but metaio can't handle
		 * unsigned ints yet */
		.waveform_number_max = LONG_MAX,
		.user_tag = NULL
	};

	return defaults;
}


static void print_usage(void)
{
	fprintf(stderr, 
"lalapps_binj [options]\n" \
"\n" \
"Options:\n" \
"\n" \
"--deltaf Hertz\n" \
"	linear spacing of injection frequencies (0.0)\n" \
"\n" \
"--fhigh Hertz\n" \
"	only inject frequencies smaller than FHIGH (1000.0)\n" \
"\n" \
"--flow Hertz\n" \
"	first frequency of injection (150.0)\n" \
"\n" \
"--fratio value\n" \
"	exponential spacing of injection frequencies (0.0)\n" \
"\n" \
"--freq-dist [monoarithmetic|monogeometric|randgeometric]\n" \
"	select the frequency distribution\n" \
"\n" \
	); fprintf(stderr, 
"--gps-end-time seconds\n" \
"	end injections at GPS time TIME\n" \
"\n" \
"--gps-start-time seconds\n" \
"	start injections at GPS time TIME\n" \
"\n" \
"--help\n" \
"	display this message\n" \
"\n" \
"--max-distance value\n" \
"	max distance of source in Kpc (default 10000Kpc)\n" \
"\n" \
"--min-distance value\n" \
"	min distance of source in Kpc (default 100Kpc)\n" \
"\n" \
"--population [galactic_core|uniform_sky|zenith]\n" \
"	select the population to synthesize\n" \
"\n" \
	); fprintf(stderr, 
"--quality value\n" \
"	quality factor for SG waveforms\n" \
"\n" \
"--seed value\n" \
"	seed random number generator with this value (0 = off, default = 0)\n" \
"\n" \
"--simwaveform-duration seconds\n" \
"	duration of the simulated waveform (Warren/Ott/ZM)\n" \
"\n" \
"--strain-dist [consthrss|hrsspresets|loghrss|loghrss-t]\n" \
"	select the strain distribution\n" \
"\n" \
"--strain-scale-max value\n" \
"	maximum of the strain distribution\n" \
"\n" \
	); fprintf(stderr, 
"--strain-scale-min value\n" \
"	mininum of the strain distribution\n" \
"\n" \
"--time-step value\n" \
"	space injections value / pi seconds appart (210)\n" \
"\n" \
"--waveform-number-max value\n" \
"--waveform-number-min value\n" \
"	waveform number is uniformly distributed in [min, max)\n" \
"\n" \
"--waveform string\n" \
"	set waveform type to NAME (SineGaussian)\n" \
"\n" \
"--user-tag string\n" \
"	set the usertag to STRING\n"
	);
}


static ProcessParamsTable **add_process_param(ProcessParamsTable **proc_param, const char *type, const char *param, const char *value)
{
	*proc_param = XLALCalloc(1, sizeof(**proc_param));
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
		{"waveform-number-max", required_argument, NULL, 'R'},
		{"waveform-number-min", required_argument, NULL, 'S'},
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
		ADD_PROCESS_PARAM("lstring");
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
		ADD_PROCESS_PARAM("lstring");
		break;

	case 'B':
		options.deltaf = atof(optarg);
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'C':
		options.fhigh = atof(optarg);
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'D':
		options.flow = atof(optarg);
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'E':
		options.fratio = atof(optarg);
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'G':
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
		ADD_PROCESS_PARAM("lstring");
		break;

	case 'H':
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
		ADD_PROCESS_PARAM("lstring");
		break;

	case 'I':
		print_usage();
		exit(0);

	case 'J':
		if(!strcmp(optarg, "consthrss"))
			options.strain_dist = STRAIN_DIST_CONST_HRSS;
		else if(!strcmp(optarg, "hrsspresets"))
			options.strain_dist = STRAIN_DIST_PRESET_HRSS;
		else if(!strcmp(optarg, "loghrss"))
			options.strain_dist = STRAIN_DIST_LOG_HRSS;
		else if(!strcmp(optarg, "loghrss-t"))
			options.strain_dist = STRAIN_DIST_LOG_HRSSTAU;
		else {
			fprintf(stderr, "error: unrecognized strain distribution \"%s\"", optarg);
			exit(1);
		}
		ADD_PROCESS_PARAM("lstring");
		break;

	case 'K':
		options.strain_scale_max = atof(optarg);
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'L':
		options.strain_scale_min = atof(optarg);
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'M':
		options.log10_max_distance = log10(atof(optarg));
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'N':
		options.log10_min_distance = log10(atof(optarg));
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'O':
		options.quality = atof(optarg);
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'P':
		options.seed = atol(optarg);
		ADD_PROCESS_PARAM("int_8u");
		break;

	case 'Q':
		options.simwaveform_duration = atof(optarg);
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'R':
		options.waveform_number_max = atol(optarg);
		ADD_PROCESS_PARAM("int_8u");
		break;

	case 'S':
		options.waveform_number_min = atol(optarg);
		ADD_PROCESS_PARAM("int_8u");
		break;

	case 'U':
		options.time_step = atof(optarg) / LAL_PI;
		ADD_PROCESS_PARAM("real_8");
		break;

	case 'V':
		options.user_tag = optarg;
		ADD_PROCESS_PARAM("lstring");
		break;

	case 'W':
		options.waveform = optarg;
		/* -1 for null terminator */
		if(strlen(options.waveform) > LIGOMETA_WAVEFORM_MAX - 1) {
			fprintf(stderr, "error: --waveform %s exceeds max length of %d characters\n", options.waveform, LIGOMETA_WAVEFORM_MAX - 1);
			exit(1);
		}
		if(!strcmp(options.waveform, "StringCusp"))
			options.freq_dist = FREQ_DIST_STRING_CUSP;
		ADD_PROCESS_PARAM("lstring");
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

	case STRAIN_DIST_CONST_HRSS:
		if(options.strain_scale_min == 0.0) {
			fprintf(stderr, "error: must set --strain-scale-min\n");
			exit(1);
		}
		if(options.strain_scale_max != 0.0)
			fprintf(stderr, "warning: --strain-scale-max provided but ignored\n");
		break;

	case STRAIN_DIST_LOG_HRSS:
	case STRAIN_DIST_LOG_HRSSTAU:
		if(options.strain_scale_min == 0.0) {
			fprintf(stderr, "error: must set --strain-scale-min\n");
			exit(1);
		}
		if(options.strain_scale_max == 0.0) {
			fprintf(stderr, "error: must set --strain-scale-max\n");
			exit(1);
		}
		break;

	case STRAIN_DIST_PRESET_HRSS:
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
 *
 *                            Solving For Unknowns
 *
 * ============================================================================
 */


static double duration_from_q_and_f(double Q, double f)
{
	/* compute duration from Q and frequency */
	return Q / (sqrt(2.0) * LAL_PI * f);
}


/* 
 * ============================================================================
 *
 *                Logarithmically-Distributed Random Variable
 *
 * ============================================================================
 */


static double Log10Deviate(gsl_rng *rng, double log10_min, double log10_max)
{
	/* FIXME:  this distribution must be in GSL */
	return pow(10.0, gsl_ran_flat(rng, log10_min, log10_max));
}


/* 
 * ============================================================================
 *
 *                          Frequency Distributions
 *
 * ============================================================================
 */


/*
 * String cusps.  What's distributed uniformly is theta^2, the square of
 * the angle the line of sight makes with direction of the cusp.
 */


static double freq_dist_string_cusp_next(gsl_rng *rng, struct options options)
{
	const double thetasqmin = pow(options.fhigh, -2.0 / 3.0);
	const double thetasqmax = pow(options.flow, -2.0 / 3.0);
	const double thetasq = gsl_ran_flat(rng, thetasqmin, thetasqmax);

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


static double freq_dist_random_geometric_next(gsl_rng *rng, struct options options)
{
	static double random_factor = 1.0;
	double freq = freq_dist_monotonic_geometric_next(options);
	if(freq * options.fratio > options.fhigh)
		freq = freq_dist_monotonic_geometric_next(options);
	if(freq == options.flow)
		random_factor = Log10Deviate(rng, 0.0, log10(options.fratio));
	return freq * random_factor;
}


/* 
 * ============================================================================
 *
 *                                hrss Presets
 *
 * ============================================================================
 */


static double hrss_preset_next(gsl_rng *rng)
{
	static const double presets[] = {
		1e-19,
		2.0e-19,
		3.0e-19,
		3.3e-18,
		5.3e-17
	};

	return presets[gsl_rng_uniform_int(rng, sizeof(presets)/sizeof(*presets))];
}


/* 
 * ============================================================================
 *
 *                                   Output
 *
 * ============================================================================
 */


static void write_xml(MetadataTable proctable, MetadataTable procparams, MetadataTable searchsumm, const SimBurst *sim_burst, struct options options)
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

	/* process table */
	XLALGPSTimeNow(&proctable.processTable->end_time);
	LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, process_table), &status);
	LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, proctable, process_table), &status);
	LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlfp), &status);

	/* process params table */
	LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, process_params_table), &status);
	LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, procparams, process_params_table), &status);
	LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlfp), &status);

	/* search summary table */
	LAL_CALL(LALBeginLIGOLwXMLTable(&status, &xmlfp, search_summary_table), &status);
	LAL_CALL(LALWriteLIGOLwXMLTable(&status, &xmlfp, searchsumm, search_summary_table), &status);
	LAL_CALL(LALEndLIGOLwXMLTable(&status, &xmlfp), &status);

	/* sim burst table */
	if(XLALWriteLIGOLwXMLSimBurstTable(&xmlfp, sim_burst)) {
		/* error occured.  ?? do anything else ?? */
		return;
	}

	LAL_CALL(LALCloseLIGOLwXMLFile(&status, &xmlfp), &status);
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
	LALStatus status = blank_status;
	MetadataTable proctable;
	MetadataTable procparams;
	MetadataTable searchsumm;
	SimBurst *sim_burst_head = NULL;
	SimBurst *sim_burst = NULL;


	/*
	 * Initialize debug handler
	 */


	lal_errhandler = LAL_ERR_EXIT;
	lalDebugLevel = LALERROR | LALNMEMDBG | LALNMEMPAD | LALNMEMTRK;


	/*
	 * Process params table and command line.
	 */


	procparams.processParamsTable = NULL;
	options = parse_command_line(&argc, &argv, &procparams);


	/*
	 * Process table
	 */


	proctable.processTable = calloc(1, sizeof(*proctable.processTable));
	XLALGPSTimeNow(&proctable.processTable->start_time);
	LAL_CALL(populate_process_table(&status, proctable.processTable, PROGRAM_NAME, CVS_REVISION, CVS_SOURCE, CVS_DATE), &status);
	snprintf(proctable.processTable->ifos, LIGOMETA_IFOS_MAX, "");
	if(options.user_tag)
		snprintf(proctable.processTable->comment, LIGOMETA_COMMENT_MAX, "%s", options.user_tag);


	/*
	 * Search summary table
	 */


	searchsumm.searchSummaryTable = calloc(1, sizeof(*searchsumm.searchSummaryTable));
	if(options.user_tag)
		snprintf(searchsumm.searchSummaryTable->comment, LIGOMETA_COMMENT_MAX, "%s", options.user_tag);
	searchsumm.searchSummaryTable->nnodes = 1;
	searchsumm.searchSummaryTable->out_start_time = *XLALINT8NSToGPS(&searchsumm.searchSummaryTable->in_start_time, options.gps_start_time);
	searchsumm.searchSummaryTable->out_end_time = *XLALINT8NSToGPS(&searchsumm.searchSummaryTable->in_end_time, options.gps_end_time);
	snprintf(searchsumm.searchSummaryTable->ifos, LIGOMETA_IFOS_MAX, proctable.processTable->ifos);
	searchsumm.searchSummaryTable->nevents = 0;


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

		/* allocate the injection */

		if(sim_burst) {
			sim_burst->next = XLALCreateSimBurst();
			sim_burst = sim_burst->next;
		} else
			sim_burst_head = sim_burst = XLALCreateSimBurst();

		/* process and simulation ids */

		sim_burst->process_id = 0;
		sim_burst->simulation_id = searchsumm.searchSummaryTable->nevents++;

		/* frequency */

		switch(options.freq_dist) {
		case FREQ_DIST_NONE:
			fprintf(stderr, "error: internal error\n");
			exit(1);

		case FREQ_DIST_STRING_CUSP:
			sim_burst->frequency = freq_dist_string_cusp_next(rng, options);
			break;

		case FREQ_DIST_MONO_ARITHMETIC:
			sim_burst->frequency = freq_dist_monotonic_arithmetic_next(options);
			break;

		case FREQ_DIST_MONO_GEOMETRIC:
			sim_burst->frequency = freq_dist_monotonic_geometric_next(options);
			break;

		case FREQ_DIST_RANDOM_GEOMETRIC:
			sim_burst->frequency = freq_dist_random_geometric_next(rng, options);
			break;
		}

		/* q, duration, and bandwidth */

		sim_burst->q = options.quality;
		if(options.simwaveform_duration)
			sim_burst->duration = options.simwaveform_duration;
		else
			sim_burst->duration = duration_from_q_and_f(options.quality, sim_burst->frequency);
		sim_burst->bandwidth = LAL_2_PI / sim_burst->duration;

		/* waveform (command line parsing code has already checked
		 * that the length is safe) */

		strcpy(sim_burst->waveform, options.waveform);

		/* peak time at geocentre in GPS seconds */

		XLALINT8NSToGPS(&sim_burst->time_geocent_gps, tinj);

		/* peak time at geocentre in radians */

		sim_burst->time_geocent_gmst = XLALGreenwichMeanSiderealTime(&sim_burst->time_geocent_gps);

		/* sky location and polarization */

		switch(options.population) {
		case POPULATION_GALACTIC_CORE:
			sim_burst->ra = MW_CENTER_J2000_RA_RAD;
			sim_burst->dec = MW_CENTER_J2000_DEC_RAD;
			sim_burst->psi = gsl_ran_flat(rng, 0, LAL_TWOPI);
			break;

		case POPULATION_UNIFORM_SKY:
			sim_burst->ra = gsl_ran_flat(rng, 0, LAL_TWOPI);
			sim_burst->dec = asin(gsl_ran_flat(rng, -1, +1));
			sim_burst->psi = gsl_ran_flat(rng, 0, LAL_TWOPI);
			break;

		case POPULATION_ZENITH:
			/* optimally oriented overhead:  ra and dec are
			 * left undefined (they were initialized to NaNs
			 * when the sim_burst was created) */
			sim_burst->psi = 0.0;
			break;
		}

		/* polarization ellipse for circularly polarized waveforms.
		 * hard-coded to linearly polarized "x" waveforms (makes
		 * LAL's sine-Gaussian generator make linearly-polarized
		 * sine-Gaussians like the old code). */

		sim_burst->pol_ellipse_angle = LAL_PI_2;
		sim_burst->pol_ellipse_e = 1.0;

		/* strain */
		/* FIXME:  must set hrss, amplitude, and egw */

		switch(options.strain_dist) {
		case STRAIN_DIST_NONE:
			fprintf(stderr, "error: internal error\n");
			exit(1);

		case STRAIN_DIST_CONST_HRSS:
			sim_burst->amplitude = sim_burst->hrss = sim_burst->egw_over_rsquared = options.strain_scale_min;
			break;

		case STRAIN_DIST_LOG_HRSS:
			sim_burst->amplitude = sim_burst->hrss = sim_burst->egw_over_rsquared = Log10Deviate(rng, options.strain_scale_min, options.strain_scale_max);
			break;

		case STRAIN_DIST_LOG_HRSSTAU:
			sim_burst->amplitude = sim_burst->hrss = sim_burst->egw_over_rsquared = Log10Deviate(rng, options.strain_scale_min - log10(sim_burst->duration), options.strain_scale_max - log10(sim_burst->duration));
			break;

		case STRAIN_DIST_PRESET_HRSS:
			sim_burst->amplitude = sim_burst->hrss = sim_burst->egw_over_rsquared = hrss_preset_next(rng);
			break;
		}

		/* pick a waveform */

		sim_burst->waveform_number = options.waveform_number_min + (unsigned long) floor(gsl_ran_flat(rng, options.waveform_number_min, options.waveform_number_max));
	}

	/* output */

	write_xml(proctable, procparams, searchsumm, sim_burst_head, options);

	/* done */

	gsl_rng_free(rng);
	exit(0);
}
