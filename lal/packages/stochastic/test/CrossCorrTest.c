/*----------------------------------------------------------------------- 
 * 
 * File Name: CrossCorrTest.c
 * 
 * Author: Steve Drasco
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
 * Tests CrossCorr()
 * 
 * DIAGNOSTICS
 * Writes PASS or FAIL to stderr as tests are passed or failed.
 *    
 * CALLS
 * CrossCorr()
 * 
 * NOTES
 *
 *-----------------------------------------------------------------------
 */

#include <string.h>
#include <math.h>
#include "LALStdlib.h"
#include "AVFactories.h"
#include "CrossCorr.h"
#include "PrintVector.h"

NRCSID (CROSSCORRTESTC, "$Id$");

#define	EXPEXTEDOUT1	1.6000e-20	/* expected output for test case 1	*/
#define EXPEXTEDOUT2    1.2937e-20	/* expected output for test case 2	*/
#define EXPEXTEDOUT3	-1.7105e-20	/* expected output for test case 3	*/
#define TOL		1.0e-2		/* allowed % disagreement 		*/

INT4 debuglevel = 0; /* set to 0 if we don't scare anyone with "failed ..." */

int check ( Status*, INT4, CHAR* );

int 
main()
{
	static	Status	status;
	CCIn		in;
	REAL4TimeSeries	g1;
	REAL4TimeSeries g2;
	REAL4		out=0.0;
	INT2		N=100;
	INT2		i;

	/* test for null input structure */
	CrossCorr(&status, &out, NULL);
	if ( check(&status, CROSSCORR_EIN, CROSSCORR_MSGEIN) ) return 1;
	fprintf(stderr,"PASS: CrossCorr test for null input structure\n");

	/* test for null output */
	CrossCorr(&status, NULL, &in);
	if ( check(&status, CROSSCORR_EOUT, CROSSCORR_MSGEOUT) ) return 1;
        fprintf(stderr,"PASS: CrossCorr test for null output\n");

	/* assemble input structure */
	g1.data = g2.data = NULL;
	in.g1 = &g1;
	in.g2 = &g2;
	in.QmaxTilde = NULL;

	/* test (1/2) null detector outputs */
	CrossCorr(&status, &out, &in);
        if ( check(&status, CROSSCORR_ENULL, CROSSCORR_MSGENULL) ) return 1;
        fprintf(stderr,"PASS: CrossCorr test (1/2) for null detector outputs\n");

	/* create vector in->g1->data */
	SCreateVector(&status, &(in.g1->data), N);

	/* test (2/2) null detector outputs */
        CrossCorr(&status, &out, &in);
	if ( check(&status, CROSSCORR_ENULL, CROSSCORR_MSGENULL) ) return 1;
        fprintf(stderr,"PASS: CrossCorr test (2/2) for null detector outputs\n");

	/* create vector in->g2->data (different length) */
	SCreateVector(&status, &(in.g2->data), N-1);

	/* test for detector outputs of unequal length */
	CrossCorr(&status, &out, &in);
        if ( check(&status, CROSSCORR_ESIZE1, CROSSCORR_MSGESIZE1) ) return 1;
        fprintf(stderr,"PASS: CrossCorr test for detector outputs of unequal length\n");

	/* equate lengths */
	SDestroyVector(&status, &(in.g2->data));
	SCreateVector(&status, &(in.g2->data),N);

	/* create kernel vector of wrong length */
	SCreateVector(&status,&in.QmaxTilde,N-1);
	
	/* test for unreasonable length of kernel vector */
	CrossCorr(&status, &out, &in);
        if ( check(&status, CROSSCORR_ESIZE2, CROSSCORR_MSGESIZE2) ) return 1;
        fprintf(stderr,"PASS: CrossCorr test for unreasonable length of kernel vector\n");

	/* correct length of kernel vector */
	SDestroyVector(&status,&in.QmaxTilde);
	SCreateVector(&status,&in.QmaxTilde,N); 

	/* change input to simple case 1 */	
	for (i=0; i < 49; i++) in.g1->data->data[i] = 0.0;
	for (i=49; i < 75; i++) in.g1->data->data[i] = 1.0e-10;
	for (i=75; i < 100;i++) in.g1->data->data[i] = 0.0;
	for (i=0; i < 59; i++) in.g2->data->data[i] = 0.0;
        for (i=59; i < 85; i++) in.g2->data->data[i] = 1.0e-11;
        for (i=85; i < 100; i++) in.g2->data->data[i] = 0.0;
	for (i=0; i< 100; i++) in.QmaxTilde->data[i] = 1.0;

	/* simple test 1/3 */
	CrossCorr(&status, &out, &in);
	fprintf(stderr," RUN: CrossCorr test (1/3) for sample data\n");
	fprintf(stderr, "\t EXPECTED: output = %e\n", EXPEXTEDOUT1);
	fprintf(stderr, "\t    FOUND: output = %e\n", out);
	fprintf(stderr, "\t      TOL: %.2e  (percent)\n", TOL);

	if( fabs( (out-EXPEXTEDOUT1)/EXPEXTEDOUT1 ) > TOL/100.0 ) {
		fprintf(stderr, "FAIL: CrossCorr test (1/3) for sample data\n");
		return 1;
	}
        fprintf(stderr,"PASS: CrossCorr test (1/3) for sample data\n");

        /* change input to simple case 2 */
        for (i=0; i < N; i++)
          in.QmaxTilde->data[i] = i < 10 ? (float) exp((double) (-i*i)/10.0) : 0;


	/* simple test 2/3 */
	CrossCorr(&status, &out, &in);
	fprintf(stderr," RUN: CrossCorr test (2/3) for sample data\n");
	fprintf(stderr, "\t EXPECTED: output = %e\n", EXPEXTEDOUT2);
	fprintf(stderr, "\t    FOUND: output = %e\n", out);
	fprintf(stderr, "\t      TOL: %.2e (percent)\n", TOL);
        if( fabs( (out-EXPEXTEDOUT2)/EXPEXTEDOUT2 ) > TOL/100.0 ) {
                fprintf(stderr, "FAIL: CrossCorr test (2/3) for sample data\n");
		return 1;
        }
        fprintf(stderr,"PASS: CrossCorr test (2/3) for sample data\n");

        /* change input to simple case 3 */
        for (i=0; i < N; i++) in.QmaxTilde->data[i] = (float) -fabs((double) i);

        /* simple test 3/3 */
        CrossCorr(&status, &out, &in);
	fprintf(stderr," RUN: CrossCorr test (3/3) for sample data\n");
        fprintf(stderr, "\t EXPECTED: output = %e\n", EXPEXTEDOUT3);
        fprintf(stderr, "\t    FOUND: output = %e\n", out);
        fprintf(stderr, "\t      TOL: %.2e (percent)\n", TOL);
        if( fabs( (out-EXPEXTEDOUT3)/EXPEXTEDOUT3 ) > TOL/100.0 ) {
                fprintf(stderr, "FAIL: CrossCorr test (3/3) for sample data\n");
                return 1;
        }
        fprintf(stderr,"PASS: CrossCorr test (3/3) for sample data\n");

	/* normal exit */
	return 0;
}

/*
 * check routine borrowed from Dirichlet
 */

int 
check( Status* status, INT4 code, CHAR* message )
{
	if ( status->statusCode != code ) {
		fprintf(stderr, "FAIL: did not recognize \"%s\" \n", message);
		return 1;
	}
	else if ( strcmp( message, status->statusDescription ) ) {
		fprintf(stderr, "FAIL: incorrect warning message \"%s\" not \"%s\" \n",status->statusDescription, message);
		return 1;
	}
	return 0;
}
