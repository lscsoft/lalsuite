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
 * Test suite for Overlap()
 * 
 * DIAGNOSTICS
 * Writes PASS or FAIL to stdout as tests are passed or failed.
 * Also writes to files the values of the overlap reduction function for 
 * all possible IFO pairs.
 *
 * CALLS
 * Overlap()
 * SCreateVector()
 * SDestroyVector()
 * 
 * NOTES
 *
 *-----------------------------------------------------------------------
 */

#ifndef _LALSTDLIB_H
#include "LALStdlib.h"
#ifndef _LALSTDLIB_H
#define _LALSTDLIB_H
#endif
#endif

#ifndef _MATH_H
#include <math.h>
#ifndef _MATH_H
#define _MATH_H
#endif
#endif

#ifndef _STRING_H
#include <string.h>
#ifndef _STRING_H
#define _STRING_H
#endif
#endif

#ifndef _PRINTVECTOR_H
#include "PrintVector.h"
#ifndef _PRINTVECTOR_H
#define _PRINTVECTOR_H
#endif
#endif

#ifndef _AVFACTORIES_H
#include "AVFactories.h"
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
#endif
#endif

#ifndef _OVERLAP_H
#include "Overlap.h"
#ifndef _OVERLAP_H
#define _OVERLAP_H
#endif
#endif

NRCSID (MAIN, "$Id$");

#define PRINT 1 /* set to 1 to write overlap reduction functions to files */ 
INT4 debuglevel = 2; /* set to 2 to get full status information for tests */

int check( Status*, INT4, CHAR* );

int main()
{
  static Status        status;
  
  REAL4Vector         *vector = NULL;
  REAL4Vector          dummy;

  OverlapParameters    parameters;

  IFOsite              site1ID, site2ID;

  SCreateVector (&status, &vector, 10000);

  /* test behavior for null parameter block */
  Overlap( &status, vector, NULL );
  if ( check( &status, OVERLAP_ENULLP, OVERLAP_MSGENULLP ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGENULLP);

  /* test behavior for null vector */
  Overlap( &status, NULL, &parameters );
  if ( check( &status, OVERLAP_ENULLV, OVERLAP_MSGENULLV ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGENULLV);

  /* test behavior for desired length <=0  */
  parameters.length = 0;
  Overlap( &status, vector, &parameters);
  if ( check( &status, OVERLAP_ESIZE, OVERLAP_MSGESIZE ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESIZE);

  /* define valid vector length */
  parameters.length = 10000;

  /* test behavior for desired frequency spacing <=0  */
  parameters.deltaF = 0;
  Overlap( &status, vector, &parameters);
  if ( check( &status, OVERLAP_EDFREQ, OVERLAP_MSGEDFREQ ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGEDFREQ);

  /* define valid frequency spacing */
  parameters.deltaF = 1.0;

  /* test behavior for undefined site ID (lower and upper bounds) */
  parameters.site1ID = -1;
  parameters.site2ID = LHO;
  Overlap( &status, vector, &parameters );
  if ( check( &status, OVERLAP_ESITE, OVERLAP_MSGESITE)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESITE);
  parameters.site1ID = NUMBEROFSITES;
  Overlap( &status, vector, &parameters );
  if ( check( &status, OVERLAP_ESITE, OVERLAP_MSGESITE)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESITE);

  parameters.site1ID = LHO;
  parameters.site2ID = -1;
  Overlap( &status, vector, &parameters );
  if ( check( &status, OVERLAP_ESITE, OVERLAP_MSGESITE)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESITE);
  parameters.site2ID = NUMBEROFSITES;
  Overlap( &status, vector, &parameters );
  if ( check( &status, OVERLAP_ESITE, OVERLAP_MSGESITE)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESITE);

  /* define valid site IDs */
  parameters.site1ID = LHO;
  parameters.site2ID = LLO;

  /* test behavior for desired vector length not equal to length specified */
  /* by input parameters */
  dummy.length = 1234;
  Overlap( &status, &dummy, &parameters );
  if ( check( &status, OVERLAP_ESZMM, OVERLAP_MSGESZMM ) ) return 1;
  printf("PASS: %s\n", OVERLAP_MSGESZMM);

  /* define valid vector length */
  dummy.length = parameters.length;

  /* test behavior for null vector data area */
  dummy.data = NULL;
  Overlap( &status, &dummy, &parameters );
  if ( check( &status, OVERLAP_ENULLD, OVERLAP_MSGENULLD)) return 1;
  printf("PASS: %s\n", OVERLAP_MSGENULLD);

  /* generate overlap reduction functions for all possible IFO pairs */
  for ( site1ID=LHO; site1ID<NUMBEROFSITES; site1ID++ ) {

    for ( site2ID=LHO; site2ID<site1ID; site2ID++ ) {
      
      parameters.site1ID = site1ID;
      parameters.site2ID = site2ID;

      /* calculate overlap reduction function */
      Overlap( &status, vector, &parameters );

      /* write overlap reduction function to file */
      if (PRINT) PrintVector( vector );
 
    } /* end site2ID loop */

  } /* end site1ID loop */

  SDestroyVector( &status, &vector);

  printf("PASS: all tests\n");

  return 0;
}

/*----------------------------------------------------------------------*/

int check( Status* status, INT4 code, CHAR* message) 
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
