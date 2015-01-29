/*
*  Copyright (C) 2007 Jolien Creighton, Kipp Cannon, Patrick Brady
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*-----------------------------------------------------------------------
 *
 * File Name: ThresholdsTest.c
 *
 * Author: Eanna Flanagan
 *
 *
 *-----------------------------------------------------------------------
 *
 * NAME
 * main()
 *
 * SYNOPSIS
 *
 * DESCRIPTION
 * Test suite for functions in Thresholds.c
 *
 * DIAGNOSTICS
 * Writes PASS or FAIL to stdout as tests are passed or failed.
 *
 * CALLS
 * LALOverlap()
 * LALSCreateVector()
 * LALSDestroyVector()
 * FindRoot()
 *
 * NOTES
 *
 *-----------------------------------------------------------------------
 */


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>

#include <lal/LALStdlib.h>
#include <lal/Thresholds.h>

extern char *optarg;
extern int optind, opterr, optopt;


int verbose = 1;


/*
 * Usage()
 *
 * Prints a usage message for program program and exits with code exitcode.
 */

static void Usage(const char *program, int exitcode)
{
	fprintf(stderr, "Usage: %s [options]\n", program);
	fprintf(stderr, "Options:\n");
	fprintf(stderr, "  -h         print this message\n");
	fprintf(stderr, "  -q         quiet: run silently\n");
	fprintf(stderr, "  -v         verbose: print extra information\n");
	fprintf(stderr, "  -d level   set lalDebugLevel to level\n");
	exit(exitcode);
}


/*
 * ParseOptions()
 *
 * Parses the argc - 1 option strings in argv[].
 */

static void ParseOptions(int argc, char *argv[])
{
	FILE *fp;
	int c;

	while(1) {
		c = getopt(argc, argv, "hqvd:");
		if(c == -1)
			break;
		switch(c) {
		case 'd':
			break;

		case 'v':
			++verbose;
			break;

		case 'q':
			fp = freopen("/dev/null", "w", stderr);
			if (fp == NULL)
			{
				fprintf(stderr, "Error: Unable to open /dev/null\n");
				 exit(1);
			}
			fp = freopen("/dev/null", "w", stdout);
			if (fp == NULL)
			{
				fprintf(stderr, "Error: Unable to open /dev/null\n");
				 exit(1);
			}
			break;

		case 'h':
			Usage(argv[0], 0);
			break;

		default:
			Usage(argv[0], 1);
		}

	}

	if(optind < argc)
		Usage(argv[0], 1);

	return;
}


/*
 * Check the output of functions
 */

#define CHECKOUTPUT(msg, expr, value, acc) { \
	REAL8 result = expr; \
	if(fabs(result - value) > acc) { \
		fprintf(stderr, msg ": expected %.11g, got %.11g\n", value, result); \
		exit(1); \
	} \
	if(XLALGetBaseErrno()) { \
		fprintf(stderr, msg ": returned error\n"); \
		exit(1); \
	} \
};


/*
 * Entry point
 */

int main(int argc, char *argv[])
{
	REAL8 chi2;
	REAL8 dof;

	/*
	 * Parse the command line options
	 */

	ParseOptions(argc, argv);


	/*
	 *  Check to make sure the functions return the correct values.
	 *  First time around
	 */

	chi2 = 2.3l;
	dof = 8.0;

	/* check forward functions */
	CHECKOUTPUT("XLALlnOneMinusChisqCdf(chi2, dof)", XLALlnOneMinusChisqCdf(chi2, dof), -0.030040797757, 1e-9);

	LALCheckMemoryLeaks();

	if(verbose)
		printf("PASS: all tests\n");

	return 0;
}
