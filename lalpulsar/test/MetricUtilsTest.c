//
// Copyright (C) 2011--2015 Reinhard Prix, Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <lal/MetricUtils.h>
#include <lal/XLALError.h>
#include <lal/UniversalDopplerMetric.h>
#include <lal/LALInitBarycenter.h>

#include <lal/GSLHelpers.h>

int main( void )
{

  // Create a phase metric to use for testing
  gsl_matrix *g_ij = NULL;
  {

    // Load ephemerides
    char earthEphem[] = TEST_DATA_DIR "earth00-19-DE200.dat.gz";
    char sunEphem[]   = TEST_DATA_DIR "sun00-19-DE200.dat.gz";
    EphemerisData *edat = XLALInitBarycenter( earthEphem, sunEphem );
    XLAL_CHECK_MAIN( edat != NULL, XLAL_EFUNC );

    // Initialise UniversalDopplerMetric parameter struct
    LIGOTimeGPS startTime = { 912345678, 0 };
    DopplerMetricParams XLAL_INIT_DECL( par );
    {
      DopplerCoordinateSystem coordSys = { 5, { DOPPLERCOORD_N3X_EQU, DOPPLERCOORD_N3Y_EQU, DOPPLERCOORD_N3Z_EQU, DOPPLERCOORD_FREQ, DOPPLERCOORD_F1DOT } };
      par.coordSys = coordSys;
    }
    par.detMotionType = DETMOTION_SPIN | DETMOTION_ORBIT;
    XLAL_CHECK_MAIN( XLALSegListInitSimpleSegments( &par.segmentList, startTime, 1, 2.5 * LAL_DAYSID_SI ) == XLAL_SUCCESS, XLAL_EFUNC );
    {
      LALStringVector *detNames = XLALCreateStringVector( "H1",  "L1", NULL );
      LALStringVector *sqrtSX   = XLALCreateStringVector( "1.0", "1.0", NULL );
      XLAL_CHECK_MAIN( XLALParseMultiLALDetector( &par.multiIFO, detNames ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLALParseMultiNoiseFloor( &par.multiNoiseFloor, sqrtSX, par.multiIFO.length ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLALDestroyStringVector( detNames );
      XLALDestroyStringVector( sqrtSX );
    }
    par.signalParams.Doppler.refTime = startTime;
    par.signalParams.Doppler.fkdot[0] = 100;
    par.projectCoord = -1;
    par.approxPhase = 1;

    // Compute phase metric
    DopplerPhaseMetric *metric = XLALComputeDopplerPhaseMetric( &par, edat );
    XLAL_CHECK_MAIN( metric != NULL, XLAL_EFUNC );
    g_ij = metric->g_ij;
    metric->g_ij = NULL;

    // Cleanup
    XLALDestroyDopplerPhaseMetric( metric );
    XLALSegListClear( &par.segmentList );
    XLALDestroyEphemerisData( edat );
  }

  // Test XLALMetricEllipseBoundingBox()
  fprintf( stderr, "\n=== Test XLALMetricEllipseBoundingBox() ===\n\n" );
  {
    gsl_vector *bbox = XLALMetricEllipseBoundingBox( g_ij, 0.3 );
    XLAL_CHECK_MAIN( bbox != NULL, XLAL_EFUNC );
    GPVEC( bbox, "%0.6g" );
    const double bbox_0 = 0.148608;
    XLAL_CHECK_MAIN( fabs( gsl_vector_get( bbox, 0 ) - bbox_0 ) <= 1e-5, XLAL_ETOL, "gsl_vector_get( bbox, 0 ) = %0.6g != %0.6g", gsl_vector_get( bbox, 0 ), bbox_0 );
    GFVEC( bbox );
  }

  // Test XLALTransformMetric() and XLALInverseTransformMetric()
  fprintf( stderr, "\n=== Test XLALTransformMetric() and XLALInverseTransformMetric() ===\n\n" );
  {

    // Allocate memory
    gsl_matrix *GAMAT_MAIN( transform, g_ij->size1, g_ij->size2 );
    gsl_matrix *GAMAT_MAIN( gpr_ij, g_ij->size1, g_ij->size2 );

    // Create some transform
    for( size_t i = 0; i < transform->size1; ++i ) {
      for( size_t j = 0; j < transform->size2; ++j ) {
        gsl_matrix_set( transform, i, j, (j >= i) ? pow(2, j - i) : 0.0 );
      }
    }

    // Apply transform then inverse transfrom, should give back same metric
    XLAL_CHECK_MAIN( XLALTransformMetric( &gpr_ij, transform, g_ij ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALInverseTransformMetric( &gpr_ij, transform, gpr_ij ) == XLAL_SUCCESS, XLAL_EFUNC );
    const REAL8 d = XLALCompareMetrics( gpr_ij, g_ij ), d_max = 1e-10;
    XLAL_CHECK_MAIN( d <= d_max, XLAL_ETOL, "XLALCompareMetric( gpr_ij, g_ij ) = %0.2g > %0.2g", d, d_max );

    // Cleanup
    GFMAT( transform, gpr_ij );
  }

  // Test XLALDiagNormalizeMetric()
  fprintf( stderr, "\n=== Test XLALDiagNormalizeMetric() ===\n\n" );
  {
    gsl_matrix *transform = NULL, *gpr_ij = NULL;

    // Diagonally normalize metric
    XLAL_CHECK_MAIN( XLALDiagNormalizeMetric( &gpr_ij, &transform, g_ij ) == XLAL_SUCCESS, XLAL_EFUNC );
    for( size_t i = 0; i < gpr_ij->size1; ++i ) {
      XLAL_CHECK_MAIN( fabs( gsl_matrix_get( gpr_ij, i, i ) - 1.0) <= 1e-5, XLAL_ETOL, "gsl_matrix_get( gpr_ij, %zu, %zu ) = %0.6g != 1.0", i, i, gsl_matrix_get( gpr_ij, i, i ) );
    }

    // Invert transform, should give back same metric
    XLAL_CHECK_MAIN( XLALInverseTransformMetric( &gpr_ij, transform, gpr_ij ) == XLAL_SUCCESS, XLAL_EFUNC );
    const REAL8 d = XLALCompareMetrics( gpr_ij, g_ij ), d_max = 1e-10;
    XLAL_CHECK_MAIN( d <= d_max, XLAL_ETOL, "XLALCompareMetric( gpr_ij, g_ij ) = %0.2g > %0.2g", d, d_max );

    // Cleanup
    GFMAT( transform, gpr_ij );
  }

  // Cleanup
  GFMAT( g_ij );
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
