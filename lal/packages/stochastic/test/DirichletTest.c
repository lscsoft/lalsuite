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
 * Test suite for Dirichlet()
 * 
 * DIAGNOSTICS
 * Writes PASS or FAIL to stdout as tests are passed or failed.
 * Also writes to a file the values of the Dirichlet kernel for threee
 * different valid test cases.
 *    
 * CALLS
 * Dirichlet()
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

#ifndef _AVFACTORIES_H
#include "AVFactories.h"
#ifndef _AVFACTORIES_H
#define _AVFACTORIES_H
#endif
#endif

#ifndef _PRINTVECTOR_H
#include "PrintVector.h"
#ifndef _PRINTVECTOR_H
#define _PRINTVECTOR_H
#endif
#endif

#ifndef _DIRICHLET_H
#include "Dirichlet.h"
#ifndef _DIRICHLET_H
#define _DIRICHLET_H
#endif
#endif

INT4 debuglevel = 2; /* set to 2 to get full status information for tests */

int check ( Status*, INT4, CHAR* );

NRCSID (MAIN, "$ID: DirichletTest.c$");

int 
main()
{
   
   Status               status;
   REAL4Vector*         poutput = NULL; 
   REAL4Vector          dummy;
   DirichletParameters  parameters;


   /* test behavior for null pointer to input parameters */
   Dirichlet( &status, &dummy, NULL );

   if ( check( &status, DIRICHLET_ENULLIP, DIRICHLET_MSGENULLIP ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENULLIP); 

   /* test behavior for Dirichlet parameter N <= 0  */
   parameters.n = -3;
   Dirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_ENVALUE, DIRICHLET_MSGENVALUE ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENVALUE);

   parameters.n = 0;
   Dirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_ENVALUE, DIRICHLET_MSGENVALUE ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENVALUE);

   /* define valid value of N */
   parameters.n = 10 ;

   /* test behavior for specified length of output vector <= 0  */
   parameters.length = -3;
   Dirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_ESIZE, DIRICHLET_MSGESIZE ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGESIZE);

   parameters.length = 0;
   Dirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_ESIZE, DIRICHLET_MSGESIZE ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGESIZE);

   /* define valid value for specified length of output vector */
   parameters.length = 11 ;

   /* test behavior for x spacing <= 0 */
   parameters.deltaX = -4;
   Dirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_EDELTAX, DIRICHLET_MSGEDELTAX ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGEDELTAX );

   parameters.deltaX = 0.0;
   Dirichlet( &status, &dummy, &parameters);
   if ( check( &status, DIRICHLET_EDELTAX, DIRICHLET_MSGEDELTAX ) ) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGEDELTAX );

   /* define valid delta x */
   parameters.deltaX = 0.1;

   /* test behavior for null pointer to output vector */
   Dirichlet( &status, NULL, &parameters );
   if ( check( &status, DIRICHLET_ENULLOP, DIRICHLET_MSGENULLOP)) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENULLOP);
      
   /* test behavior for length of output vector not equal to length  */
   /* specified in input parameters */
   dummy.length = 10; 
   Dirichlet( &status, &dummy, &parameters );
   if ( check( &status, DIRICHLET_ESIZEMM, DIRICHLET_MSGESIZEMM ) ) return 1;
   printf( "PASS: %s\n", DIRICHLET_MSGESIZEMM );

   /* assign valid output vector length */
   dummy.length = parameters.length;

   /* test behavior for null pointer to data member of output vector */
   dummy.data = NULL;
   Dirichlet( &status, &dummy, &parameters );
   if ( check( &status, DIRICHLET_ENULLD, DIRICHLET_MSGENULLD)) return 1;
   printf("PASS: %s\n", DIRICHLET_MSGENULLD);

   /* VALID TEST DATA #1 */
   /* call Dirichet() with valid data (N=even) */
   parameters.n      = 10;
   parameters.length = 101;
   parameters.deltaX = 0.01; 
   SCreateVector (&status, &poutput, parameters.length);  

   Dirichlet( &status, poutput, &parameters );  
   PrintVector(poutput); 

   SDestroyVector( &status, &poutput );   

   /* VALID TEST DATA #2 */
   /* call Dirichet() with valid data (N=odd) */
   parameters.n      = 11;
   parameters.length = 101;
   parameters.deltaX = 0.01; 
   SCreateVector(&status, &poutput, parameters.length);  

   Dirichlet( &status, poutput, &parameters );  
   PrintVector(poutput); 

   SDestroyVector( &status, &poutput );   

   /* VALID TEST DATA #3 */
   /* call Dirichet() with valid data (x=0 to 2) */
   parameters.n      = 10;
   parameters.length = 201;
   parameters.deltaX = 0.01; 
   SCreateVector(&status, &poutput, parameters.length);  

   Dirichlet( &status, poutput, &parameters );  
   PrintVector(poutput); 

   SDestroyVector( &status, &poutput );   

   return 0;
}
/*------------------------------------------------------------------------*/

int 
check( Status* status, INT4 code, CHAR* message )
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

