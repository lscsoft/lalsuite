/*----------------------------------------------------------------------- 
 * 
 * File Name: OverlapTest.c
 * 
 * Author: J.D. Romano
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
 * Test suite for LALOverlap()
 * 
 * DIAGNOSTICS
 * Writes PASS or FAIL to stdout as tests are passed or failed.
 * Also writes to files the values of the overlap reduction function for 
 * all possible IFO pairs.
 *
 * CALLS
 * LALOverlap()
 * LALSCreateVector()
 * LALSDestroyVector()
 * 
 * NOTES
 *
 *-----------------------------------------------------------------------
 */

#include "LALStdlib.h"
#include <math.h>
#include <string.h>
#include "AVFactories.h"
#include "PrintVector.h"
#include "Overlap.h"

#define PRINT 1 /* set to 1 to write overlap reduction functions to files */ 
INT4 LALDebugLevel = 2; /* set to 2 to get full status information for tests */

int check( LALStatus*, INT4, CHAR* );

NRCSID (MAIN, "$Id$");

int main()
{
  LALStatus               status = {0};
  
  REAL4Vector         *vector = NULL;
  REAL4Vector          dummy;

  OverlapParameters    parameters;

  IFOsite              site1ID, site2ID;

  /* test behavior for null pointer to input parameters */
  LALOverlap( &status, &dummy, NULL );
  if ( check( &status, OVERLAP_ENULLIP, OVERLAP_MSGENULLIP ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGENULLIP);

  /* test behavior for illegitimate site IDs (lower and upper bounds) */
  parameters.site1ID = -1;
  parameters.site2ID = LHO;
  LALOverlap( &status, &dummy, &parameters );
  if ( check( &status, OVERLAP_ESITE, OVERLAP_MSGESITE)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESITE);
  parameters.site1ID = NUMBEROFSITES;
  LALOverlap( &status, &dummy, &parameters );
  if ( check( &status, OVERLAP_ESITE, OVERLAP_MSGESITE)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESITE);

  parameters.site1ID = LHO;
  parameters.site2ID = -1;
  LALOverlap( &status, &dummy, &parameters );
  if ( check( &status, OVERLAP_ESITE, OVERLAP_MSGESITE)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESITE);
  parameters.site2ID = NUMBEROFSITES;
  LALOverlap( &status, &dummy, &parameters );
  if ( check( &status, OVERLAP_ESITE, OVERLAP_MSGESITE)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESITE);

  /* define valid site IDs */
  parameters.site1ID = LHO;
  parameters.site2ID = LLO;

  /* test behavior for specified length of output vector <= 0  */
  parameters.length = -4;
  LALOverlap( &status, &dummy, &parameters);
  if ( check( &status, OVERLAP_ESIZE, OVERLAP_MSGESIZE ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESIZE);

  parameters.length = 0;
  LALOverlap( &status, &dummy, &parameters);
  if ( check( &status, OVERLAP_ESIZE, OVERLAP_MSGESIZE ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESIZE);

  /* define valid vector length */
  parameters.length = 10000;

  /* test behavior for desired frequency spacing <= 0  */
  parameters.deltaF = -1.0;
  LALOverlap( &status, &dummy, &parameters);
  if ( check( &status, OVERLAP_EDELTAF, OVERLAP_MSGEDELTAF ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGEDELTAF);

  parameters.deltaF = 0.0;
  LALOverlap( &status, &dummy, &parameters);
  if ( check( &status, OVERLAP_EDELTAF, OVERLAP_MSGEDELTAF ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGEDELTAF);

  /* define valid frequency spacing */
  parameters.deltaF = 1.0;

  /* test behavior for null pointer to output vector */
  LALOverlap( &status, NULL, &parameters );
  if ( check( &status, OVERLAP_ENULLOP, OVERLAP_MSGENULLOP ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGENULLOP);

  /* test behavior for length of output vector not equal to length */
  /* specified in input parameters */
  dummy.length = 1234;
  LALOverlap( &status, &dummy, &parameters );
  if ( check( &status, OVERLAP_ESIZEMM, OVERLAP_MSGESIZEMM ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESIZEMM);

  /* define valid vector length */
  dummy.length = parameters.length;

  /* test behavior for null pointer to data member of output vector */
  dummy.data = NULL;
  LALOverlap( &status, &dummy, &parameters );
  if ( check( &status, OVERLAP_ENULLD, OVERLAP_MSGENULLD)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGENULLD);

  /* VALID DATA HERE ------------------------------------------------ */

  LALSCreateVector (&status, &vector, 10000);

  /* generate overlap reduction functions for all possible IFO pairs */
  for ( site1ID=LHO; site1ID<NUMBEROFSITES; site1ID++ ) {

    for ( site2ID=LHO; site2ID<site1ID; site2ID++ ) {
      
      parameters.site1ID = site1ID;
      parameters.site2ID = site2ID;

      /* calculate overlap reduction function */
      LALOverlap( &status, vector, &parameters );

      /* write overlap reduction function to file */
      if (PRINT) LALPrintVector( vector );
 
    } /* end site2ID loop */

  } /* end site1ID loop */

  LALSDestroyVector( &status, &vector);

  printf("PASS: all tests\n");

  return 0;
}

/*----------------------------------------------------------------------*/

int check( LALStatus* status, INT4 code, CHAR* message) 
{
  if ( status->statusCode!= code ) {
    printf( "FAIL: did not recognize \"%s\"\n", message );
    return 1;
  }
  else if ( strcmp( message, status->statusDescription ) ) {
    printf( "FAIL: incorrect warning message \"%s\" not \"%s\"\n",
             status->statusDescription, message);
    return 1;
  }

  return 0;
}
