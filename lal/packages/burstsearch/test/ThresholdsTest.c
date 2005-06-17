/*-----------------------------------------------------------------------
 *
 * File Name: ThresholdsTest.c
 *
 * Author: Eanna Flanagan
 *
 * Revision: $Id$
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

#include <lal/LALErrno.h>
#include <lal/LALStdlib.h>
#include <lal/Thresholds.h>

NRCSID(MAIN, "$Id$");

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
	int c;

	while(1) {
		c = getopt(argc, argv, "hqvd:");
		if(c == -1)
			break;
		switch(c) {
		case 'd':
			lalDebugLevel = atoi(optarg);
			break;

		case 'v':
			++verbose;
			break;

		case 'q':
			freopen("/dev/null", "w", stderr);
			freopen("/dev/null", "w", stdout);
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


#define CHECKERROR(msg, expr, error) { \
	REAL8 result = expr; \
	if(!XLALIsREAL8FailNaN(result)) { \
		fprintf(stderr, msg ": failed to return REAL8FailNaN\n"); \
		exit(1); \
	} \
	if(XLALGetBaseErrno() != error) { \
		fprintf(stderr, msg ": expected error %d, got %d\n", error, XLALGetBaseErrno()); \
		exit(1); \
	} \
	XLALClearErrno(); \
};


/*
 * Entry point
 */

int main(int argc, char *argv[])
{
	REAL8 chi2;
	REAL8 dof;
	REAL8 rho;
	REAL8 falseAlarm;
	REAL8 falseDismissal;

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
	rho = 2.2;

	/* check forward functions */
	CHECKOUTPUT("1.0 - XLALChisqCdf(chi2, dof)", 1.0 - XLALChisqCdf(chi2, dof), 0.970406, 1e-5);
	CHECKOUTPUT("XLALOneMinusChisqCdf(chi2, dof)", XLALOneMinusChisqCdf(chi2, dof), 0.970406, 1e-5);
	CHECKOUTPUT("XLALNoncChisqCdf(chi2, dof, rho * rho)", XLALNoncChisqCdf(chi2, dof, rho * rho), 0.00439452, 1e-7);
	/* check reverse functions */
	CHECKOUTPUT("XLALChi2Threshold(dof, 1.0 - XLALChisqCdf(chi2, dof))", XLALChi2Threshold(dof, 1.0 - XLALChisqCdf(chi2, dof)), chi2, 1e-5);
#if 0
	/* FIXME */
	CHECKOUTPUT("XLALRhoThreshold(chi2, dof, XLALNoncChisqCdf(chi2, dof, rho * rho))", XLALRhoThreshold(chi2, dof, XLALNoncChisqCdf(chi2, dof, rho * rho)), rho, 1e-5);
#endif


	/* 
	 *  Check to make sure the functions return the correct values.  
	 *  Second time around 
	 */

	chi2 = 12.3l;
	dof = 3.1;
	rho = 2.2;

	/* check forward functions */
	CHECKOUTPUT("1.0 - XLALChisqCdf(chi2, dof)", 1.0 - XLALChisqCdf(chi2, dof), 0.007066, 1e-7);
	CHECKOUTPUT("XLALNoncChisqCdf(chi2, dof, rho * rho)", XLALNoncChisqCdf(chi2, dof, rho * rho), 0.822575, 1e-6);
	/* check reverse functions */
	CHECKOUTPUT("XLALChi2Threshold(dof, 1.0 - XLALChisqCdf(chi2, dof))", XLALChi2Threshold(dof, 1.0 - XLALChisqCdf(chi2, dof)), chi2, 1e-5);
#if 0
	/* FIXME */
	CHECKOUTPUT("XLALRhoThreshold(chi2, dof, XLALNoncChisqCdf(chi2, dof, rho * rho))", XLALRhoThreshold(chi2, dof, XLALNoncChisqCdf(chi2, dof, rho * rho)), rho, 1e-5);
#endif


	/*
	 * Check to make sure that correct error codes are generated.
	 */

#ifndef LAL_NDEBUG
	falseAlarm = 0.970406;
	falseDismissal = 0.00439452;
	if(!lalNoDebug) {
		if(verbose || lalDebugLevel)
			printf("\n===== Check Errors =====\n");

		CHECKERROR("XLALChisqCdf(-chi2, dof)", XLALChisqCdf(-chi2, dof), XLAL_EDOM);
		CHECKERROR("XLALNoncChisqCdf(-chi2, dof, rho * rho)", XLALNoncChisqCdf(-chi2, dof, rho * rho), XLAL_EDOM);
		CHECKERROR("XLALChisqCdf(chi2, -dof)", XLALChisqCdf(chi2, -dof), XLAL_EDOM);
		CHECKERROR("XLALNoncChisqCdf(chi2, -dof, rho * rho)", XLALNoncChisqCdf(chi2, -dof, rho * rho), XLAL_EDOM);
		CHECKERROR("XLALNoncChisqCdf(chi2, dof, -rho * rho)", XLALNoncChisqCdf(chi2, dof, -rho * rho), XLAL_EDOM);
		CHECKERROR("XLALChi2Threshold(-dof, falseAlarm)", XLALChi2Threshold(-dof, falseAlarm), XLAL_EDOM);
		CHECKERROR("XLALRhoThreshold(chi2, -dof, falseDismissal)", XLALRhoThreshold(chi2, -dof, falseDismissal), XLAL_EDOM);
		CHECKERROR("XLALRhoThreshold(-chi2, dof, falseDismissal)", XLALRhoThreshold(-chi2, dof, falseDismissal), XLAL_EDOM);

		/* There is no test here for exceeding the maximum number
		 * of iterations in XLALNoncChisqCdf() since I could not
		 * find a set of parameters which caused this condition to
		 * occur. */

		/* Supplied probabilities must lie between 0 and 1 */
		CHECKERROR("XLALChi2Threshold(dof, -falseAlarm)", XLALChi2Threshold(dof, -falseAlarm), XLAL_EDOM);
		CHECKERROR("XLALChi2Threshold(dof, 2.0)", XLALChi2Threshold(dof, 2.0), XLAL_EDOM);
		CHECKERROR("XLALRhoThreshold(chi2, dof, -falseDismissal)", XLALRhoThreshold(chi2, dof, -falseDismissal), XLAL_EDOM);
		CHECKERROR("XLALRhoThreshold(chi2, dof, -2.0)", XLALRhoThreshold(chi2, dof, -2.0), XLAL_EDOM);
	}
#endif

	LALCheckMemoryLeaks();

	if(verbose)
		printf("PASS: all tests\n");

	return 0;
}
