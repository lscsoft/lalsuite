/*
 * Copyright (C) 2007 Saikat Ray-Majumder
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/Date.h>
#include <lal/LIGOLwXML.h>
#include <lal/LIGOLwXMLRead.h>
#include <lal/LIGOMetadataTables.h>
#include <lal/LIGOMetadataUtils.h>
#include <lal/LALCache.h>
#include <lal/LALInspiral.h>
#include <lalapps.h>

#define MAXSTR 2048

/* Usage format string. */
#define USAGE "Usage: %s --mass1 mass(Msun) --mass2 mass(Msun)" \
        "--flow low_freq(Hz)" \
	"[--help]\n"

#define MSUNINSEC ( 4.9255*1e-6 )  /* solar mass in sec */



static REAL4
XLALMergerFrequency(
	InspiralTemplate *params
)
{
	REAL4 fmerger;
    
	fmerger = 205 * (20/params->totalMass);

	return(fmerger);
}


static REAL4
XLALQnrFrequency(
	InspiralTemplate *params
)
{
	REAL4 fqnr;
    
	fqnr = 1320 * (20/params->totalMass);

	return(fqnr);
}


static REAL4
XLALMergerEnergyFraction(
	InspiralTemplate *params
)
{
	REAL4 e_merger = 0.0;
    
	e_merger = 0.1 * 16 * params->mu * params->mu /params->totalMass;

	return(e_merger);
}


static void 
XLALMergerDuration(
	InspiralTemplate *params,
	MergerDurationLimits *dur
)
{
	dur->high = 50 * params->totalMass;
	dur->low = 10 * params->totalMass;
}




/*
 * Command-line options
 */


struct options_t {
	int verbose;
	REAL8 mass1;
	REAL8 mass2;
	REAL8 flow;
};


/*
 * Set defaults for command-line options.
 * 
 */

static void set_option_defaults(struct options_t *options)
{
	options->verbose = 0;
	options->mass1 = 0.0;
	options->mass2 = 0.0;
	options->flow = 0.0;	
}


/*
 * Parse command line arguments.
 */

static void parse_command_line(int argc, char **argv, struct options_t *options)
{
	struct option long_options[] = {
		/* these options set a flag */
		{"verbose",         no_argument,        &options->verbose, 1},
		/* parameters which determine the output xml file */
		{"mass1",        required_argument,  NULL,  'a'},
		{"mass2",        required_argument,  NULL,  'b'},
		{"flow",         required_argument,  NULL,  'c'},
		{"help",         no_argument,        NULL,  'o'}, 
		{NULL, 0, NULL, 0}
	};
	int c;
	int option_index;


	do {
		switch(c = getopt_long(argc, argv, "a:b:c:", long_options, &option_index)) {
			case -1:
			case 0:
			break;

			case 'a':
			/*
			 * file containing list of xml files to use
			 */
			options->mass1 = atoi(optarg);
			break;

			case 'b':
			/*
			 * output cache file name
			 */
			options->mass2 = atoi(optarg);
			break;

			case 'c':
			/*
			 * output cache file name
			 */
			options->flow = atof(optarg);
			break;	

			case ':':
			case '?':
			case 'o':
			default:
			/*
			 * print usage
			 */
			LALPrintError(USAGE, *argv);
			exit(1);
		}
	} while(c != -1);

	if(optind < argc) {
		fprintf(stderr, "extraneous command line arguments:\n");
		while(optind < argc)
			fprintf(stderr, "%s\n", argv[optind++]);
		exit(1);
	}

}


/*
 * Entry Point
 */


int main(int argc, char **argv)
{
	static LALStatus  stat;
	struct options_t  options;
	InspiralTemplate  params;
	REAL4 fmerger = 0.0;
	REAL4 fqnr = 0.0;
	REAL4 e_merger = 0.0;
	MergerDurationLimits  dur;

	FILE *fp = NULL;
	/*
	 * Initialize things
	 */

	set_option_defaults(&options);
	parse_command_line(argc, argv, &options);

	params.mass1 = options.mass1;
	params.mass2 = options.mass2;
	params.fLower = options.flow;
	params.massChoice = m1Andm2;
       
	/* generate the template parameters */
	LAL_CALL( LALInspiralParameterCalc(&stat, &params), &stat);

	params.mu = (params.mass1 * params.mass2)/params.totalMass;

	/*merger frequency */
	fmerger = XLALMergerFrequency(&params);

	/*qnr frequency */
	fqnr = XLALQnrFrequency(&params);

	/*merger energy fraction */
	e_merger = XLALMergerEnergyFraction(&params);

	/*merger duration limits in units of M_sun*/
	XLALMergerDuration(&params, &dur);

	fp = fopen("MergerParams.dat","w");
	fprintf(fp,"mass1(Msun) mass2(Msun) totalMass(Msun) mu(Msun) t_insp(seconds) f_merger(Hz) f_qnr(Hz) e_merger(Msun) duration(max)(seconds) duration(minimum)(seconds)\n");
	fprintf(fp,"%f %f %f %f %f %f %f %f %f %f\n",params.mass1, params.mass2, params.totalMass, params.mu, params.t0, fmerger, fqnr, e_merger, dur.high * MSUNINSEC, dur.low * MSUNINSEC); 
	fclose(fp);

	exit(0);
}
