/*----------------------------------------------------------------------------
 *
 * File Name: GetErrorMatrixFromSnglInspiral.c
 *
 * Author: Craig Robinson
 *
 * $Id$
 *
 *---------------------------------------------------------------------------*/

#if 0
<lalVerbatim file="LALTrigScanClusterCV">
Author: Craig Robinson
$Id$
</lalVerbatim>
#endif


#include <lal/LALStdlib.h>
#include <lal/LALError.h>
#include <lal/LALGSL.h>
#include <lal/LIGOMetadataTables.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>


/* Function for getting the error matrix from the metric in
 * (tc, tau0, tau3) space.
 */

gsl_matrix * XLALGetErrorMatrixFromSnglInspiral(SnglInspiralTable *event)
{
  static const char *func = "XLALGetErrorMatrixFromSnglInspiral";
  gsl_matrix *shape  = NULL;
  gsl_matrix *fisher = NULL;
  gsl_permutation *p = NULL;

  int signum;

  if (!event)
  {
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  }

  if ( !event->Gamma[0] )
  {
    XLALPrintError( "Metric components are not set.\n" );
    XLAL_ERROR_NULL( func, XLAL_EINVAL );
  }

  /* Allocate memory for the various matrices */
  shape  = gsl_matrix_alloc( 3, 3 );
  fisher = gsl_matrix_alloc( 3, 3 );
  p      = gsl_permutation_alloc( 3 );

  /* Fill in the elements of the fisher matrix */
  gsl_matrix_set( fisher, 0, 0, event->Gamma[0] );
  gsl_matrix_set( fisher, 0, 1, event->Gamma[1] );
  gsl_matrix_set( fisher, 1, 0, event->Gamma[1] );
  gsl_matrix_set( fisher, 0, 2, event->Gamma[2] );
  gsl_matrix_set( fisher, 2, 0, event->Gamma[2] );
  gsl_matrix_set( fisher, 1, 1, event->Gamma[3] );
  gsl_matrix_set( fisher, 1, 2, event->Gamma[4] );
  gsl_matrix_set( fisher, 2, 1, event->Gamma[4] );
  gsl_matrix_set( fisher, 2, 2, event->Gamma[5] );

  gsl_matrix_scale( fisher, 1.0 / 0.3 );

  /* Now invert to get the matrix we need */
  gsl_linalg_LU_decomp( fisher, p, &signum );
  gsl_linalg_LU_invert( fisher, p, shape );

  gsl_matrix_free( fisher );
  gsl_permutation_free( p );
  return shape;
}


/* Returns the position vector in (tc, tau0, tau3) space */
gsl_vector * XLALGetPositionFromSnglInspiral( SnglInspiralTable *table )
{
  static const char *func = "XLALGetPositionFromSnglInspiral";
  gsl_vector *position = NULL;
  REAL8 endTime;

  XLAL_CALLGSL( position = gsl_vector_alloc( 3 ) );
  if ( !position )
    XLAL_ERROR_NULL( func, XLAL_ENOMEM );

  endTime = (REAL8) table->end_time.gpsSeconds +
        (REAL8) table->end_time.gpsNanoSeconds * 1.0e-9;

  gsl_vector_set( position, 0, endTime );
  gsl_vector_set( position, 1, table->tau0 );
  gsl_vector_set( position, 2, table->tau3 );

  return position;
}
