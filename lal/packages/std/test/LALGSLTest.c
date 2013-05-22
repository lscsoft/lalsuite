/*
*  Copyright (C) 2007 Jolien Creighton
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

/**
 * \file
 * \ingroup LALGSL_h
 * \author Creighton, J. D. E.

 This program tests the LAL macros for GSL function calls.  It makes sure
 that a nominal status is returned if the GSL function succeeds, and that
 an error code is returned if the GSL function fails.

*/

/** \cond DONT_DOXYGEN */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <lal/LALStdlib.h>
#include <lal/LALGSL.h>
#include <lal/LALConstants.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/* macro for testing status */
#define TESTSTATUS( pstatus, code ) \
  if ( (pstatus)->statusCode != code ) { REPORTSTATUS(pstatus); exit(1); } \
  else ((void)(0))

/* macro for testing handler */
extern gsl_error_handler_t *gsl_error_handler;
gsl_error_handler_t *original_handler;
#define TESTHANDLER \
  if ( original_handler != gsl_error_handler ) \
    { fprintf( stderr, "Error: handler was not restored!\n" ); exit(2); } \
  else ((void)(0))

/* function for clearing status pointer */
static void ClearStatus( LALStatus *status )
{
  if ( status->statusPtr )
  {
    ClearStatus( status->statusPtr );
    DETATCHSTATUSPTR( status );
  }
  return;
}


/*
 *
 * Basic tests: call the logarithm with both positive and negative values.
 *
 */


/* call the GSL log routine to test the GSL macros */
static void Logarithm( LALStatus *status, REAL8 *output, REAL8 input )
{
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  CALLGSL( *output = gsl_sf_log( input ), status );
  CHECKSTATUSPTR( status );

  /* just for fun, use the TRYGSL macro to do the same thing too */
  TRYGSL( *output = gsl_sf_log( input ), status );

  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/*
 *
 * More complex tests: integrate the Heaviside function, but code the
 * Heaviside function stupidly so that it can cause failures.  This tests
 * that the GSL-macros are reenterant-safe.
 *
 */


/* here is a really stupid Heaviside function: take the log of x;
 * return 1 if no error otherwise return 0 */
static double Heaviside( double x, void UNUSED * params )
{
  LALStatus newStatus;
  double dummy;
  double ans;

  /* blank the new status structure */
  memset( &newStatus, 0, sizeof( newStatus ) );

  /* call the LAL function */
  Logarithm( &newStatus, &dummy, x );
  if ( newStatus.statusCode )
  {
    ClearStatus( &newStatus );
    ans = 0;
  }
  else
    ans = 1;

  return ans;
}

static void Integral( LALStatus *status, REAL8 *y, REAL8 a, REAL8 b, REAL8 eps )
{
  gsl_integration_workspace *work = NULL;
  gsl_function F;
  REAL8  err;
  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );

  F.function = &Heaviside;
  F.params   = NULL;

  TRYGSL( work = gsl_integration_workspace_alloc( 100 ), status );
  CALLGSL( gsl_integration_qags( &F, a, b, eps, eps, 100, work, y, &err ),
      status );
  BEGINFAIL( status )
    TRYGSL( gsl_integration_workspace_free( work ), status );
  ENDFAIL( status );
  TRYGSL( gsl_integration_workspace_free( work ), status );

  /* if eps is set too small the integration is supposed to fail, but on
   * some systems it may happen that it doesn't... in this case, just make
   * this routine fail with a recursive error */
  if ( eps < LAL_REAL8_EPS )
  {
    TRY( Logarithm( status->statusPtr, y, 0 ), status );
  }

  DETATCHSTATUSPTR( status );
  RETURN( status );
}



int main( void )
{
  static LALStatus status;
  double x;
  double y;

  original_handler = gsl_error_handler;

  /* these are simple tests of a LAL routine that calls a GSL function */

  /* this should succeed */
  x = 10;
  Logarithm( &status, &y, x );
  TESTSTATUS( &status, 0 );             /* expect error code 0 */
  TESTHANDLER;

  /* this should fail */
  x = 0;
  Logarithm( &status, &y, x );
  TESTSTATUS( &status, -1 );            /* expect error code -1 */
  TESTHANDLER;
  ClearStatus( &status );

  /* this should succeed again */
  x = 10;
  Logarithm( &status, &y, x );
  TESTSTATUS( &status, 0 );             /* expect error code 0 */
  TESTHANDLER;


  /* these are more complicated tests of a LAL routine that calls a GSL
   * function that ultimately calls another LAL function, which in turn
   * calls a GSL function... to make sure that everything is reenterant */

  /* this should pass */
  Integral( &status, &y, -1, 1, 1e-7 );
  TESTSTATUS( &status, 0 );             /* expect error code 0 */
  TESTHANDLER;

  /* this should fail */
  Integral( &status, &y, -1, 1, LAL_REAL8_MIN );
  TESTSTATUS( &status, -1 );             /* expect error code -1 */
  TESTHANDLER;
  ClearStatus( &status );

  LALCheckMemoryLeaks();

  return 0;
}

/** \endcond */
