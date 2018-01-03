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

#include <lal/LALStdlib.h>
#include <lal/LALgetopt.h>
#include <lal/LALChisq.h>

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
		c = LALgetopt(argc, argv, "hqvd:");
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

	if(LALoptind < argc)
		Usage(argv[0], 1);

	return;
}


/*
 * Check the output of functions
 */

#define CHECKOUTPUT(msg, expr, value, acc) do { \
	double result = expr; \
	if(fabs((result - value) / value) > acc) { \
		fprintf(stderr, "%s:  expected %.17g, got %.17g (fractional error is %g;  upto %g allowed)\n", msg, value, result, fabs((result - value) / value), acc); \
		exit(1); \
	} else if(verbose) { \
		fprintf(stderr, "%s:  expected %.17g, got %.17g (fractional error is %g;  upto %g allowed)\n", msg, value, result, fabs((result - value) / value), acc); \
	} \
	if(XLALGetBaseErrno()) { \
		fprintf(stderr, "%s:  returned error\n", msg); \
		exit(1); \
	} \
} while(0)


#define CHECKXLALLogChisqCCDF(chi2, dof, value, acc) do { \
	char msg[100]; \
	sprintf(msg, "XLALLogChisqCCDF(%.17g, %.17g)", chi2, dof); \
	CHECKOUTPUT(msg, XLALLogChisqCCDF(chi2, dof), value, acc); \
} while(0)


/*
 * Entry point
 */

int main(int argc, char *argv[])
{
	/*
	 * Parse the command line options
	 */

	ParseOptions(argc, argv);

	/*
	 * Check to make sure the functions return the correct values.
	 * "Correct" values obtained with Mathematica, computing
	 * intermediate results to 1000 digits.  E.g., the second result
	 * can be obtained with Mathematica using the expression
	 *
	 * N[Log[1-Q[8 / 2, 0, 2.3 / 2]], 1000]
	 */

	CHECKXLALLogChisqCCDF(2., 64., -1.4417345421413976e-36, 1e-15);
	CHECKXLALLogChisqCCDF(2.3, 8., -0.030040797756978235, 1e-15);
	CHECKXLALLogChisqCCDF(8., 8., -0.83593241162679427, 1e-15);
	CHECKXLALLogChisqCCDF(2.3, 0.5, -2.9095189371057191, 1e-15);
	CHECKXLALLogChisqCCDF(1.2e3, 8., -582.59596635081904, 1e-15);
	CHECKXLALLogChisqCCDF(2e4, 1e4, -1539.4420486763690, 1e-15);

	/*
	 * Done.
	 */

	if(verbose)
		printf("PASS: all tests\n");

	return 0;
}
