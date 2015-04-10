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

#include "../src/GSLHelpers.h"

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
    par.metricType = METRIC_TYPE_PHASE;
    par.approxPhase = 1;
    par.nonposEigValThresh = 1;

    // Compute phase metric
    DopplerMetric *metric = XLALDopplerFstatMetric( &par, edat );
    XLAL_CHECK_MAIN( metric != NULL, XLAL_EFUNC );
    g_ij = metric->g_ij;
    metric->g_ij = NULL;

    // Cleanup
    XLALDestroyDopplerMetric( metric );
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

  // Cleanup
  GFMAT( g_ij );
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
