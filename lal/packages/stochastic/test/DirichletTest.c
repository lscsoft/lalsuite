/*----------------------------------------------------------------------- 
 * 
 * File Name: DirichletTest.c
 * 
 * Author: UTB Relativity Group
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
 * Test suite for LALDirichlet()
 * 
 * DIAGNOSTICS
 * Writes PASS or FAIL to stdout as tests are passed or failed.
 * Also writes to a file the values of the LALDirichlet kernel for threee
 * different valid test cases.
 *    
 * CALLS
 * LALDirichlet()
 * LALSCreateVector()
 * LALSDestroyVector()
 * 
 * NOTES
 *
 *-----------------------------------------------------------------------
 */

#include <lal/LALStdlib.h>
#include <math.h>
#include <string.h>
#include <lal/AVFactories.h>
#include <lal/PrintVector.h>
#include <lal/Dirichlet.h>

INT4 lalDebugLevel = 2; /* set to 2 to get full status information for tests */

int check ( LALStatus*, INT4, CHAR* );

NRCSID (MAIN, "$ID: DirichletTest.c$");

int 
main()
{
   
   LALStatus               status  = {0};
   REAL4Vector*         poutput = NULL; 
   REAL4Vector          dummy;
   DirichletParameters  parameters;


   /* test behavior for null pointer to input parameters */
   LALDirichlet( &status, &dummy, NULL );

   if ( check( &status, DIRICHLET_ENULLIP, DIRICHLET_MSGENULLIP ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENULLIP); 

   /* test behavior for LALDirichlet parameter N <= 0  */
   parameters.n = -3;
   LALDirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_ENVALUE, DIRICHLET_MSGENVALUE ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENVALUE);

   parameters.n = 0;
   LALDirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_ENVALUE, DIRICHLET_MSGENVALUE ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENVALUE);

   /* define valid value of N */
   parameters.n = 10 ;

   /* test behavior for specified length of output vector <= 0  */
   parameters.length = -3;
   LALDirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_ESIZE, DIRICHLET_MSGESIZE ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGESIZE);

   parameters.length = 0;
   LALDirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_ESIZE, DIRICHLET_MSGESIZE ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGESIZE);

   /* define valid value for specified length of output vector */
   parameters.length = 11 ;

   /* test behavior for x spacing <= 0 */
   parameters.deltaX = -4;
   LALDirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_EDELTAX, DIRICHLET_MSGEDELTAX ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGEDELTAX );

   parameters.deltaX = 0.0;
   LALDirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_EDELTAX, DIRICHLET_MSGEDELTAX ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGEDELTAX );

   /* define valid delta x */
   parameters.deltaX = 0.1;

   /* test behavior for null pointer to output vector */
   LALDirichlet( &status, NULL, &parameters );
   if ( check( &status, DIRICHLET_ENULLOP, DIRICHLET_MSGENULLOP)) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENULLOP);
      
   /* test behavior for length of output vector not equal to length  */
   /* specified in input parameters */
   dummy.length = 10; 
   LALDirichlet( &status, &dummy, &parameters );
   if ( check( &status, DIRICHLET_ESIZEMM, DIRICHLET_MSGESIZEMM ) ) return 1;
   printf( "PASS: %s\n", DIRICHLET_MSGESIZEMM );

   /* assign valid output vector length */
   dummy.length = parameters.length;

   /* test behavior for null pointer to data member of output vector */
   dummy.data = NULL;
   LALDirichlet( &status, &dummy, &parameters );
   if ( check( &status, DIRICHLET_ENULLD, DIRICHLET_MSGENULLD)) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENULLD);

   /* VALID TEST DATA #1 */
   /* call Dirichet() with valid data (N=even) */
   parameters.n      = 10;
   parameters.length = 101;
   parameters.deltaX = 0.01; 
   LALSCreateVector (&status, &poutput, parameters.length);  

   LALDirichlet( &status, poutput, &parameters );  
   LALPrintVector(poutput); 

   LALSDestroyVector( &status, &poutput );   

   /* VALID TEST DATA #2 */
   /* call Dirichet() with valid data (N=odd) */
   parameters.n      = 11;
   parameters.length = 101;
   parameters.deltaX = 0.01; 
   LALSCreateVector(&status, &poutput, parameters.length);  

   LALDirichlet( &status, poutput, &parameters );  
   LALPrintVector(poutput); 

   LALSDestroyVector( &status, &poutput );   

   /* VALID TEST DATA #3 */
   /* call Dirichet() with valid data (x=0 to 2) */
   parameters.n      = 10;
   parameters.length = 201;
   parameters.deltaX = 0.01; 
   LALSCreateVector(&status, &poutput, parameters.length);  

   LALDirichlet( &status, poutput, &parameters );  
   LALPrintVector(poutput); 

   LALSDestroyVector( &status, &poutput );   

   return 0;
}
/*------------------------------------------------------------------------*/

int 
check( LALStatus* status, INT4 code, CHAR* message )
{
   if ( status->statusCode!= code ) {
     printf ( "FAIL: did not recognize \"%s\"\n", message );
     return 1;
   }
   else if ( strcmp( message, status->statusDescription ) ) {
     printf( "FAIL: incorrect warning message \"%s\" not \"%s\"\n",
	     status->statusDescription, message );

     return 1;
   }

   return 0;
}

