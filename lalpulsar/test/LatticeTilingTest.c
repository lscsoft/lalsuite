//
// Copyright (C) 2014, 2015, 2016 Karl Wette
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

// Tests of the lattice-based template generation code in LatticeTiling.[ch].

#include <config.h>
#include <stdio.h>

#include <lal/LatticeTiling.h>
#include <lal/LALStdlib.h>
#include <lal/Factorial.h>
#include <lal/DopplerFullScan.h>
#include <lal/SuperskyMetrics.h>
#include <lal/LALInitBarycenter.h>

#include <lal/GSLHelpers.h>

#if defined(__GNUC__)
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#define MISM_HIST_BINS 20

// The reference mismatch histograms were generated in Octave
// using the LatticeMismatchHist() function available in OctApps.

const double Z1_mism_hist[MISM_HIST_BINS] = {
  4.531107, 1.870257, 1.430467, 1.202537, 1.057047, 0.953084, 0.875050, 0.813050, 0.762368, 0.719968,
  0.683877, 0.652659, 0.625394, 0.601300, 0.579724, 0.560515, 0.542944, 0.527142, 0.512487, 0.499022
};

const double Z2_mism_hist[MISM_HIST_BINS] = {
  1.570963, 1.571131, 1.571074, 1.571102, 1.570808, 1.570789, 1.570617, 1.570716, 1.570671, 1.570867,
  1.157132, 0.835785, 0.645424, 0.503305, 0.389690, 0.295014, 0.214022, 0.143584, 0.081427, 0.025878
};

const double Z3_mism_hist[MISM_HIST_BINS] = {
  0.608404, 1.112392, 1.440652, 1.705502, 1.934785, 2.139464, 2.296868, 2.071379, 1.748278, 1.443955,
  1.155064, 0.879719, 0.616210, 0.375368, 0.223752, 0.131196, 0.071216, 0.033130, 0.011178, 0.001489
};

#define A1s_mism_hist Z1_mism_hist

const double A2s_mism_hist[MISM_HIST_BINS] = {
  1.210152, 1.210142, 1.209837, 1.209697, 1.209368, 1.209214, 1.209399, 1.209170, 1.208805, 1.208681,
  1.208631, 1.208914, 1.208775, 1.209021, 1.208797, 0.816672, 0.505394, 0.315665, 0.170942, 0.052727
};

const double A3s_mism_hist[MISM_HIST_BINS] = {
  0.327328, 0.598545, 0.774909, 0.917710, 1.040699, 1.150991, 1.250963, 1.344026, 1.431020, 1.512883,
  1.590473, 1.664510, 1.595423, 1.391209, 1.194340, 1.004085, 0.729054, 0.371869, 0.098727, 0.011236
};

const double A4s_mism_hist[MISM_HIST_BINS] = {
  0.088295, 0.264916, 0.441209, 0.617937, 0.794537, 0.971005, 1.147715, 1.324035, 1.500356, 1.677569,
  1.806866, 1.816272, 1.757854, 1.653638, 1.513900, 1.268203, 0.833153, 0.417934, 0.100320, 0.004287
};

static int SerialisationTest(
  const LatticeTiling UNUSED *tiling,
  const UINT8 UNUSED total_ref,
  const int UNUSED total_tol,
  const UINT8 UNUSED total_ckpt_0,
  const UINT8 UNUSED total_ckpt_1,
  const UINT8 UNUSED total_ckpt_2,
  const UINT8 UNUSED total_ckpt_3
  )
{

#if !defined(HAVE_LIBCFITSIO)
  printf( "Skipping serialisation test (CFITSIO library is not available)\n" );
#else // defined(HAVE_LIBCFITSIO)
  printf( "Performing serialisation test ..." );

  const size_t n = XLALTotalLatticeTilingDimensions( tiling );

  const UINT8 total_ckpt[4] = {total_ckpt_0, total_ckpt_1, total_ckpt_2, total_ckpt_3};

  // Create lattice tiling iterator
  LatticeTilingIterator *itr = XLALCreateLatticeTilingIterator( tiling, n );
  XLAL_CHECK( itr != NULL, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingAlternatingIterator( itr, true ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Count number of points
  const UINT8 total = XLALTotalLatticeTilingPoints( itr );
  XLAL_CHECK( total > 0, XLAL_EFUNC );
  XLAL_CHECK( imaxabs( total - total_ref ) <= total_tol, XLAL_EFUNC, "|total - total_ref| = |%" LAL_UINT8_FORMAT " - %" LAL_UINT8_FORMAT "| > %i", total, total_ref, total_tol );

  // Get all points
  gsl_matrix *GAMAT( points, n, total );
  XLAL_CHECK( XLALNextLatticeTilingPoints( itr, &points ) == ( int ) total, XLAL_EFUNC );
  XLAL_CHECK( XLALNextLatticeTilingPoint( itr, NULL ) == 0, XLAL_EFUNC );
  XLAL_CHECK( XLALResetLatticeTilingIterator( itr ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Iterate over all points again, this time with checkpointing
  gsl_vector *GAVEC( point, n );
  size_t k_ckpt = 0;
  for ( UINT8 k = 0; k < total; ++k ) {

    // Check next point for consistency
    XLAL_CHECK( XLALNextLatticeTilingPoint( itr, point ) >= 0, XLAL_EFUNC );
    gsl_vector_const_view points_k_view = gsl_matrix_const_column( points, k );
    gsl_vector_sub( point, &points_k_view.vector );
    double err = gsl_blas_dasum( point ) / n;
    XLAL_CHECK( err < 1e-6, XLAL_EFAILED, "err = %e < 1e-6", err );

    // Checkpoint iterator at certain intervals
    if ( k_ckpt < XLAL_NUM_ELEM( total_ckpt ) && k + 1 >= total_ckpt[k_ckpt] ) {

      // Save iterator to a FITS file
      {
        FITSFile *file = XLALFITSFileOpenWrite( "LatticeTilingTest.fits" );
        XLAL_CHECK( file != NULL, XLAL_EFUNC );
        XLAL_CHECK( XLALSaveLatticeTilingIterator( itr, file, "itr" ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLALFITSFileClose( file );
      }

      // Destroy and recreate lattice tiling iterator
      XLALDestroyLatticeTilingIterator( itr );
      itr = XLALCreateLatticeTilingIterator( tiling, n );
      XLAL_CHECK( itr != NULL, XLAL_EFUNC );
      XLAL_CHECK( XLALSetLatticeTilingAlternatingIterator( itr, true ) == XLAL_SUCCESS, XLAL_EFUNC );

      // Restore iterator to a FITS file
      {
        FITSFile *file = XLALFITSFileOpenRead( "LatticeTilingTest.fits" );
        XLAL_CHECK( file != NULL, XLAL_EFUNC );
        XLAL_CHECK( XLALRestoreLatticeTilingIterator( itr, file, "itr" ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLALFITSFileClose( file );
      }

      printf( " checkpoint at %" LAL_UINT8_FORMAT "/%" LAL_UINT8_FORMAT " ...", k + 1, total );
      ++k_ckpt;

    }

  }
  XLAL_CHECK( XLALNextLatticeTilingPoint( itr, NULL ) == 0, XLAL_EFUNC );

  printf( " done\n" );

  // Cleanup
  XLALDestroyLatticeTilingIterator( itr );
  GFVEC( point );
  GFMAT( points );

#endif // !defined(HAVE_LIBCFITSIO)

  return XLAL_SUCCESS;

}

static int BasicTest(
  const size_t n,
  const int bound_on_0,
  const int bound_on_1,
  const int bound_on_2,
  const int bound_on_3,
  const TilingLattice lattice,
  const UINT8 total_ref_0,
  const UINT8 total_ref_1,
  const UINT8 total_ref_2,
  const UINT8 total_ref_3
  )
{

  const int total_tol = 1;
  const double value_tol = 1000 * LAL_REAL8_EPS;

  const int bound_on[4] = {bound_on_0, bound_on_1, bound_on_2, bound_on_3};
  const UINT8 total_ref[4] = {total_ref_0, total_ref_1, total_ref_2, total_ref_3};

  // Create lattice tiling
  LatticeTiling *tiling = XLALCreateLatticeTiling( n );
  XLAL_CHECK( tiling != NULL, XLAL_EFUNC );

  // Add bounds
  for ( size_t i = 0; i < n; ++i ) {
    XLAL_CHECK( bound_on[i] == 0 || bound_on[i] == 1, XLAL_EFAILED );
    XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, i, 0.0, bound_on[i] * pow( 100.0, 1.0/n ) ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Set metric to the Lehmer matrix
  const double max_mismatch = 0.3;
  {
    gsl_matrix *GAMAT( metric, n, n );
    for ( size_t i = 0; i < n; ++i ) {
      for ( size_t j = 0; j < n; ++j ) {
        const double ii = i+1, jj = j+1;
        gsl_matrix_set( metric, i, j, jj >= ii ? ii/jj : jj/ii );
      }
    }
    XLAL_CHECK( XLALSetTilingLatticeAndMetric( tiling, lattice, metric, max_mismatch ) == XLAL_SUCCESS, XLAL_EFUNC );
    GFMAT( metric );
    printf( "Number of (tiled) dimensions: %zu (%zu)\n", XLALTotalLatticeTilingDimensions( tiling ), XLALTiledLatticeTilingDimensions( tiling ) );
    printf( "  Bounds: %i %i %i %i\n", bound_on_0, bound_on_1, bound_on_2, bound_on_3 );
    printf( "  Lattice type: %i\n", lattice );
  }

  // Check tiled status of lattce tiling dimensions
  for ( size_t i = 0, ti = 0; i < n; ++i ) {
    const int is_tiled_i = XLALIsTiledLatticeTilingDimension( tiling, i );
    XLAL_CHECK( is_tiled_i >= 0, XLAL_EFUNC );
    XLAL_CHECK( !is_tiled_i == !bound_on[i], XLAL_EFAILED, "XLALIsTiledLatticeTilingDimension(tiling, %zu) = %i, should be %i", i, is_tiled_i, bound_on[i] );
    if ( is_tiled_i ) {
      const size_t j = XLALLatticeTilingTiledDimension( tiling, ti );
      XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );
      XLAL_CHECK( i == j, XLAL_EFAILED, "XLALLatticeTilingTiledDimension( tiling, %zu ) = %zu, should be %zu", ti, j, i );
      ++ti;
    }
  }

  // Create lattice tiling locator
  LatticeTilingLocator *loc = XLALCreateLatticeTilingLocator( tiling );
  XLAL_CHECK( loc != NULL, XLAL_EFUNC );
  if ( lalDebugLevel & LALINFOBIT ) {
    printf( "  Index trie:\n" );
    XLAL_CHECK( XLALPrintLatticeTilingIndexTrie( loc, stdout ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  for ( size_t i = 0; i < n; ++i ) {

    // Create lattice tiling iterator over 'i+1' dimensions
    LatticeTilingIterator *itr = XLALCreateLatticeTilingIterator( tiling, i+1 );
    XLAL_CHECK( itr != NULL, XLAL_EFUNC );

    // Count number of points
    const UINT8 total = XLALTotalLatticeTilingPoints( itr );
    XLAL_CHECK( total > 0, XLAL_EFUNC );
    printf( "Number of lattice points in %zu dimensions: %" LAL_UINT8_FORMAT " (vs %" LAL_UINT8_FORMAT ", tolerance = %i)\n", i+1, total, total_ref[i], total_tol );
    XLAL_CHECK( imaxabs( total - total_ref[i] ) <= total_tol, XLAL_EFUNC, "|total - total_ref[%zu]| = |%" LAL_UINT8_FORMAT " - %" LAL_UINT8_FORMAT "| > %i", i, total, total_ref[i], total_tol );
    for ( UINT8 k = 0; XLALNextLatticeTilingPoint( itr, NULL ) > 0; ++k ) {
      const UINT8 itr_index = XLALCurrentLatticeTilingIndex( itr );
      XLAL_CHECK( k == itr_index, XLAL_EFUNC, "k = %" LAL_UINT8_FORMAT " != %" LAL_UINT8_FORMAT " = itr_index", k, itr_index );
    }
    XLAL_CHECK( XLALResetLatticeTilingIterator( itr ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Check tiling statistics
    printf( "  Check tiling statistics ..." );
    for ( size_t j = 0; j < n; ++j ) {
      const LatticeTilingStats *stats = XLALLatticeTilingStatistics( tiling, j );
      XLAL_CHECK( stats != NULL, XLAL_EFUNC );
      XLAL_CHECK( stats->name != NULL, XLAL_EFUNC );
      XLAL_CHECK( imaxabs( stats->total_points - total_ref[j] ) <= total_tol, XLAL_EFAILED, "|total - total_ref[%zu]| = |%" LAL_UINT8_FORMAT " - %" LAL_UINT8_FORMAT "| > %i", j, stats->total_points, total_ref[j], total_tol );
      XLAL_CHECK( stats->min_points <= stats->max_points, XLAL_EFAILED, "min_points = %" LAL_INT4_FORMAT " > %" LAL_INT4_FORMAT " = max_points", stats->min_points, stats->max_points );
      XLAL_CHECK( stats->min_value <= stats->max_value, XLAL_EFAILED, "min_value = %g > %g = max_value", stats->min_value, stats->max_value );
      printf( " %s ...", stats->name );
    }
    printf( " done\n" );

    // Get all points
    gsl_matrix *GAMAT( points, n, total );
    XLAL_CHECK( XLALNextLatticeTilingPoints( itr, &points ) == ( int )total, XLAL_EFUNC );
    XLAL_CHECK( XLALNextLatticeTilingPoint( itr, NULL ) == 0, XLAL_EFUNC );
    for ( UINT8 k = 0; k < total; ++k ) {
      gsl_vector_const_view point_view = gsl_matrix_const_column( points, k );
      const gsl_vector *point = &point_view.vector;
      for ( size_t j = 0; j < n; ++j ) {
        const double point_j = gsl_vector_get( point, j );
        const LatticeTilingStats *stats = XLALLatticeTilingStatistics( tiling, j );
        XLAL_CHECK( point_j >= stats->min_value - value_tol, XLAL_EFAILED, "point_j = %.10g < %.10g = stats[%zu]->min_value", point_j, stats->min_value, j );
        XLAL_CHECK( point_j <= stats->max_value + value_tol, XLAL_EFAILED, "point_j = %.10g > %.10g = stats[%zu]->max_value", point_j, stats->max_value, j );
      }
    }

    // Get nearest points to each template, check for consistency
    printf( "  Testing XLALNearestLatticeTiling{Point|Block}() ..." );
    gsl_vector *GAVEC( nearest, n );
    UINT8Vector *nearest_indexes = XLALCreateUINT8Vector( n );
    XLAL_CHECK( nearest_indexes != NULL, XLAL_ENOMEM );
    for ( UINT8 k = 0; k < total; ++k ) {
      gsl_vector_const_view point_view = gsl_matrix_const_column( points, k );
      const gsl_vector *point = &point_view.vector;
      XLAL_CHECK( XLALNearestLatticeTilingPoint( loc, point, nearest, nearest_indexes ) == XLAL_SUCCESS, XLAL_EFUNC );
      for ( size_t j = 0; j < n; ++j ) {
        const double nearest_j = gsl_vector_get( nearest, j );
        const LatticeTilingStats *stats = XLALLatticeTilingStatistics( tiling, j );
        XLAL_CHECK( nearest_j >= stats->min_value - value_tol, XLAL_EFAILED, "nearest_j = %.10g < %.10g = stats[%zu]->min_value", nearest_j, stats->min_value, j );
        XLAL_CHECK( nearest_j <= stats->max_value + value_tol, XLAL_EFAILED, "nearest_j = %.10g > %.10g = stats[%zu]->max_value", nearest_j, stats->max_value, j );
      }
      gsl_vector_sub( nearest, point );
      double err = gsl_blas_dasum( nearest ) / n;
      XLAL_CHECK( err < 1e-6, XLAL_EFAILED, "err = %e < 1e-6", err );
      XLAL_CHECK( nearest_indexes->data[i] == k, XLAL_EFAILED, "nearest_indexes[%zu] = %" LAL_UINT8_FORMAT " != %" LAL_UINT8_FORMAT "\n", i, nearest_indexes->data[i], k );
      if ( 0 < i ) {
        const LatticeTilingStats *stats = XLALLatticeTilingStatistics( tiling, i );
        UINT8 nearest_index = 0;
        INT4 nearest_left = 0, nearest_right = 0;
        XLAL_CHECK( XLALNearestLatticeTilingBlock( loc, point, i, nearest, &nearest_index, &nearest_left, &nearest_right ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLAL_CHECK( nearest_index == nearest_indexes->data[i-1], XLAL_EFAILED, "nearest_index = %" LAL_UINT8_FORMAT " != %" LAL_UINT8_FORMAT "\n", nearest_index, nearest_indexes->data[i-1] );
        XLAL_CHECK( nearest_left <= nearest_right, XLAL_EFAILED, "invalid [nearest_left, nearest_right] = [%i, %i]\n", nearest_left, nearest_right );
        UINT4 nearest_len = nearest_right - nearest_left + 1;
        XLAL_CHECK( nearest_len <= stats->max_points, XLAL_EFAILED, "nearest_len = %i > %i = stats[%zu]->max_points\n", nearest_len, stats->max_points, i );
      }
      if ( i+1 < n ) {
        const LatticeTilingStats *stats = XLALLatticeTilingStatistics( tiling, i+1 );
        UINT8 nearest_index = 0;
        INT4 nearest_left = 0, nearest_right = 0;
        XLAL_CHECK( XLALNearestLatticeTilingBlock( loc, point, i+1, nearest, &nearest_index, &nearest_left, &nearest_right ) == XLAL_SUCCESS, XLAL_EFUNC );
        XLAL_CHECK( nearest_index == nearest_indexes->data[i], XLAL_EFAILED, "nearest_index = %" LAL_UINT8_FORMAT " != %" LAL_UINT8_FORMAT "\n", nearest_index, nearest_indexes->data[i] );
        XLAL_CHECK( nearest_left <= nearest_right, XLAL_EFAILED, "invalid [nearest_left, nearest_right] = [%i, %i]\n", nearest_left, nearest_right );
        UINT4 nearest_len = nearest_right - nearest_left + 1;
        XLAL_CHECK( nearest_len <= stats->max_points, XLAL_EFAILED, "nearest_len = %i > %i = stats[%zu]->max_points\n", nearest_len, stats->max_points, i+1 );
      }
    }
    printf( " done\n" );

    // Cleanup
    XLALDestroyLatticeTilingIterator( itr );
    GFMAT( points );
    GFVEC( nearest );
    XLALDestroyUINT8Vector( nearest_indexes );

    // Create alternating lattice tiling iterator over 'i+1' dimensions
    printf( "  Testing XLALSetLatticeTilingAlternatingIterator() ..." );
    LatticeTilingIterator *itr_alt = XLALCreateLatticeTilingIterator( tiling, i+1 );
    XLAL_CHECK( itr_alt != NULL, XLAL_EFUNC );
    XLAL_CHECK( XLALSetLatticeTilingAlternatingIterator( itr_alt, true ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Count number of points, check for consistency with non-alternating count
    UINT8 total_alt = 0;
    while ( XLALNextLatticeTilingPoint( itr_alt, NULL ) > 0 ) {
      ++total_alt;
    }
    XLAL_CHECK( imaxabs( total_alt - total_ref[i] ) <= total_tol, XLAL_EFUNC, "alternating |total - total_ref[%zu]| = |%" LAL_UINT8_FORMAT " - %" LAL_UINT8_FORMAT "| > %i", i, total_alt, total_ref[i], total_tol );
    printf( " done\n" );

    // Cleanup
    XLALDestroyLatticeTilingIterator( itr_alt );

  }

  // Perform serialisation test
  XLAL_CHECK( SerialisationTest( tiling, total_ref[n-1], total_tol, total_ref_0, total_ref_1, total_ref_2, total_ref_3 ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Cleanup
  XLALDestroyLatticeTiling( tiling );
  XLALDestroyLatticeTilingLocator( loc );
  LALCheckMemoryLeaks();
  printf( "\n" );

  return XLAL_SUCCESS;

}

static int MismatchTest(
  const LatticeTiling *tiling,
  const gsl_matrix *metric,
  const double max_mismatch,
  const double injs_per_point,
  const double mism_hist_error_tol,
  const double mism_out_of_range_tol,
  const UINT8 total_ref,
  const int total_tol,
  const double mism_hist_ref[MISM_HIST_BINS]
  )
{

  const size_t n = XLALTotalLatticeTilingDimensions( tiling );

  // Create lattice tiling iterator and locator
  LatticeTilingIterator *itr = XLALCreateLatticeTilingIterator( tiling, n );
  XLAL_CHECK( itr != NULL, XLAL_EFUNC );
  LatticeTilingLocator *loc = XLALCreateLatticeTilingLocator( tiling );
  XLAL_CHECK( loc != NULL, XLAL_EFUNC );

  // Count number of points
  const UINT8 total = XLALTotalLatticeTilingPoints( itr );
  XLAL_CHECK( total > 0, XLAL_EFUNC );
  printf( "Number of lattice points: %" LAL_UINT8_FORMAT " (vs %" LAL_UINT8_FORMAT ", tolerance = %i)\n", total, total_ref, total_tol );
  XLAL_CHECK( imaxabs( total - total_ref ) <= total_tol, XLAL_EFUNC, "|total - total_ref| = |%" LAL_UINT8_FORMAT " - %" LAL_UINT8_FORMAT "| > %i", total, total_ref, total_tol );

  // Initialise mismatch histogram counts
  double mism_hist[MISM_HIST_BINS] = {0};
  double mism_hist_total = 0, mism_hist_out_of_range = 0;

  // Perform 'injs_per_point' injections for every template
  printf( "Injections per point: %g\n", injs_per_point );
  {
    UINT8 total_injs = llround(total * injs_per_point);
    const size_t max_injs_per_batch = 100000;

    // Allocate memory
    gsl_matrix *GAMAT( max_injections, n, max_injs_per_batch );
    gsl_matrix *max_nearest = NULL;
    gsl_matrix *GAMAT( max_temp, n, max_injs_per_batch );

    // Allocate random number generator
    RandomParams *rng = XLALCreateRandomParams( total );
    XLAL_CHECK( rng != NULL, XLAL_EFUNC );

    while ( total_injs > 0 ) {
      const size_t injs_per_batch = GSL_MIN( max_injs_per_batch, total_injs );
      total_injs -= injs_per_batch;
      printf("  Injections remaining: %" LAL_UINT8_FORMAT "...\n", total_injs);

      // Create matrix views
      gsl_matrix_view injections_view = gsl_matrix_submatrix( max_injections, 0, 0, n, injs_per_batch );
      gsl_matrix *const injections = &injections_view.matrix;
      gsl_matrix_view temp_view = gsl_matrix_submatrix( max_temp, 0, 0, n, injs_per_batch );
      gsl_matrix *const temp = &temp_view.matrix;

      // Generate random injection points
      XLAL_CHECK( XLALRandomLatticeTilingPoints( tiling, 0.0, rng, injections ) == XLAL_SUCCESS, XLAL_EFUNC );

      // Find nearest lattice template points
      XLAL_CHECK( XLALNearestLatticeTilingPoints( loc, injections, &max_nearest, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      gsl_matrix_view nearest_view = gsl_matrix_submatrix( max_nearest, 0, 0, n, injs_per_batch );
      gsl_matrix *const nearest = &nearest_view.matrix;

      // Compute mismatch between injections
      gsl_matrix_sub( nearest, injections );
      gsl_blas_dsymm( CblasLeft, CblasUpper, 1.0, metric, nearest, 0.0, temp );
      for ( size_t j = 0; j < temp->size2; ++j ) {
        gsl_vector_view temp_j = gsl_matrix_column( temp, j );
        gsl_vector_view nearest_j = gsl_matrix_column( nearest, j );
        double mismatch = 0.0;
        gsl_blas_ddot( &nearest_j.vector, &temp_j.vector, &mismatch );
        mismatch /= max_mismatch;

        // Increment mismatch histogram counts
        ++mism_hist_total;
        if ( mismatch < 0.0 || mismatch > 1.0 ) {
          ++mism_hist_out_of_range;
        } else {
          ++mism_hist[lround( floor( mismatch * MISM_HIST_BINS ) )];
        }

      }

    }

    // Cleanup
    XLALDestroyRandomParams( rng );
    GFMAT( max_injections, max_nearest, max_temp );

  }

  // Normalise histogram
  for ( size_t i = 0; i < MISM_HIST_BINS; ++i ) {
    mism_hist[i] *= MISM_HIST_BINS / mism_hist_total;
  }

  // Print mismatch histogram and its reference
  printf( "Mismatch histogram: " );
  for ( size_t i = 0; i < MISM_HIST_BINS; ++i ) {
    printf( " %0.3f", mism_hist[i] );
  }
  printf( "\n" );
  printf( "Reference histogram:" );
  for ( size_t i = 0; i < MISM_HIST_BINS; ++i ) {
    printf( " %0.3f", mism_hist_ref[i] );
  }
  printf( "\n" );

  // Determine error between mismatch histogram and its reference
  double mism_hist_error = 0.0;
  for ( size_t i = 0; i < MISM_HIST_BINS; ++i ) {
    mism_hist_error += fabs( mism_hist[i] - mism_hist_ref[i] );
  }
  mism_hist_error /= MISM_HIST_BINS;
  printf( "Mismatch histogram error: %0.3e (tolerance %0.3e)\n", mism_hist_error, mism_hist_error_tol );
  if ( mism_hist_error >= mism_hist_error_tol ) {
    XLAL_ERROR( XLAL_EFAILED, "mismatch histogram error exceeds tolerance\n" );
  }

  // Check fraction of injections out of histogram range
  const double mism_out_of_range = mism_hist_out_of_range / mism_hist_total;
  printf( "Fraction of points out of histogram range: %0.3e (tolerance %0.3e)\n", mism_out_of_range, mism_out_of_range_tol );
  if ( mism_out_of_range > mism_out_of_range_tol ) {
    XLAL_ERROR( XLAL_EFAILED, "fraction of points out of histogram range exceeds tolerance\n" );
  }

  // Perform 10 injections outside parameter space
  {
    gsl_matrix *GAMAT( injections, n, 10 );
    gsl_matrix *GAMAT( nearest, n, 10 );
    RandomParams *rng = XLALCreateRandomParams( total );
    XLAL_CHECK( rng != NULL, XLAL_EFUNC );

    // Generate random injection points outside parameter space
    XLAL_CHECK( XLALRandomLatticeTilingPoints( tiling, 5.0, rng, injections ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Find nearest lattice template points
    XLAL_CHECK( XLALNearestLatticeTilingPoints( loc, injections, &nearest, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Cleanup
    GFMAT( injections, nearest );
    XLALDestroyRandomParams( rng );
  }

  // Cleanup
  XLALDestroyLatticeTilingIterator( itr );
  XLALDestroyLatticeTilingLocator( loc );

  return XLAL_SUCCESS;

}

static int MismatchSquareTest(
  const TilingLattice lattice,
  const double freqband,
  const double f1dotband,
  const double f2dotband,
  const UINT8 total_ref,
  const double mism_hist_ref[MISM_HIST_BINS]
  )
{

  const int total_tol = 1;

  // Create lattice tiling
  LatticeTiling *tiling = XLALCreateLatticeTiling( 3 );
  XLAL_CHECK( tiling != NULL, XLAL_EFUNC );

  // Add bounds
  const double fndot[3] = {100, 0, 0};
  const double fndotband[3] = {freqband, f1dotband, f2dotband};
  for ( size_t i = 0; i < 3; ++i ) {
    printf( "Bounds: f%zudot=%0.3g, f%zudotband=%0.3g\n", i, fndot[i], i, fndotband[i] );
    XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, i, fndot[i], fndot[i] + fndotband[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Set metric to the spindown metric
  const double max_mismatch = 0.3;
  gsl_matrix *GAMAT( metric, 3, 3 );
  for ( size_t i = 0; i < metric->size1; ++i ) {
    for ( size_t j = i; j < metric->size2; ++j ) {
      const double Tspan = 432000;
      const double metric_i_j_num = 4.0 * LAL_PI * LAL_PI * pow( Tspan, i + j + 2 ) * ( i + 1 ) * ( j + 1 );
      const double metric_i_j_denom = LAL_FACT[i + 1] * LAL_FACT[j + 1] * ( i + 2 ) * ( j + 2 ) * ( i + j + 3 );
      gsl_matrix_set( metric, i, j, metric_i_j_num / metric_i_j_denom );
      gsl_matrix_set( metric, j, i, gsl_matrix_get( metric, i, j ) );
    }
  }
  printf( "Lattice type: %i\n", lattice );
  XLAL_CHECK( XLALSetTilingLatticeAndMetric( tiling, lattice, metric, max_mismatch ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Perform mismatch test
  XLAL_CHECK( MismatchTest( tiling, metric, max_mismatch, 10, 5e-2, 2e-3, total_ref, total_tol, mism_hist_ref ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Perform serialisation test
  XLAL_CHECK( SerialisationTest( tiling, total_ref, total_tol, 1, 0.2*total_ref, 0.6*total_ref, 0.9*total_ref ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Cleanup
  XLALDestroyLatticeTiling( tiling );
  GFMAT( metric );
  LALCheckMemoryLeaks();
  printf( "\n" );

  return XLAL_SUCCESS;

}

static int MismatchAgeBrakeTest(
  const TilingLattice lattice,
  const double freq,
  const double freqband,
  const UINT8 total_ref,
  const double mism_hist_ref[MISM_HIST_BINS]
  )
{

  const int total_tol = 1;

  // Create lattice tiling
  LatticeTiling *tiling = XLALCreateLatticeTiling( 3 );
  XLAL_CHECK( tiling != NULL, XLAL_EFUNC );

  // Add bounds
  printf( "Bounds: freq=%0.3g, freqband=%0.3g\n", freq, freqband );
  XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, 0, freq, freq + freqband ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingF1DotAgeBrakingBound( tiling, 0, 1, 1e11, 2, 5 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingF2DotBrakingBound( tiling, 0, 1, 2, 2, 5 ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Set metric to the spindown metric
  const double max_mismatch = 0.3;
  gsl_matrix *GAMAT( metric, 3, 3 );
  for ( size_t i = 0; i < metric->size1; ++i ) {
    for ( size_t j = i; j < metric->size2; ++j ) {
      const double Tspan = 1036800;
      const double metric_i_j_num = 4.0 * LAL_PI * LAL_PI * pow( Tspan, i + j + 2 ) * ( i + 1 ) * ( j + 1 );
      const double metric_i_j_denom = LAL_FACT[i + 1] * LAL_FACT[j + 1] * ( i + 2 ) * ( j + 2 ) * ( i + j + 3 );
      gsl_matrix_set( metric, i, j, metric_i_j_num / metric_i_j_denom );
      gsl_matrix_set( metric, j, i, gsl_matrix_get( metric, i, j ) );
    }
  }
  printf( "Lattice type: %i\n", lattice );
  XLAL_CHECK( XLALSetTilingLatticeAndMetric( tiling, lattice, metric, max_mismatch ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Perform mismatch test
  XLAL_CHECK( MismatchTest( tiling, metric, max_mismatch, 10, 5e-2, 2e-3, total_ref, total_tol, mism_hist_ref ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Perform serialisation test
  XLAL_CHECK( SerialisationTest( tiling, total_ref, total_tol, 1, 0.2*total_ref, 0.6*total_ref, 0.9*total_ref ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Cleanup
  XLALDestroyLatticeTiling( tiling );
  GFMAT( metric );
  LALCheckMemoryLeaks();
  printf( "\n" );

  return XLAL_SUCCESS;

}

static int SuperskyTests(
  const UINT8 coh_total_ref_0,
  const UINT8 coh_total_ref_1,
  const UINT8 coh_total_ref_2,
  const UINT8 semi_total_ref
  )
{

  const int total_tol = 15;

  const UINT8 coh_total_ref[3] = {coh_total_ref_0, coh_total_ref_1, coh_total_ref_2};

  printf( "Performing super-sky metric tests ...\n\n" );

  // Compute reduced supersky metrics
  const double Tspan = 90000;
  LIGOTimeGPS ref_time;
  XLALGPSSetREAL8( &ref_time, 900100100 );
  LALSegList segments;
  {
    XLAL_CHECK( XLALSegListInit( &segments ) == XLAL_SUCCESS, XLAL_EFUNC );
    LALSeg segment;
    {
      LIGOTimeGPS start_time = ref_time, end_time = ref_time;
      XLALGPSAdd( &start_time, -100 * Tspan );
      XLALGPSAdd( &end_time, -98 * Tspan );
      XLAL_CHECK( XLALSegSet( &segment, &start_time, &end_time, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALSegListAppend( &segments, &segment ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    {
      LIGOTimeGPS start_time = ref_time, end_time = ref_time;
      XLALGPSAdd( &start_time, -0.5 * Tspan );
      XLALGPSAdd( &end_time, 0.5 * Tspan );
      XLAL_CHECK( XLALSegSet( &segment, &start_time, &end_time, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALSegListAppend( &segments, &segment ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    {
      LIGOTimeGPS start_time = ref_time, end_time = ref_time;
      XLALGPSAdd( &start_time, 92.5 * Tspan );
      XLALGPSAdd( &end_time, 93.5 * Tspan );
      XLAL_CHECK( XLALSegSet( &segment, &start_time, &end_time, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALSegListAppend( &segments, &segment ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    XLAL_CHECK( segments.length == XLAL_NUM_ELEM( coh_total_ref ), XLAL_EFAILED );
  }
  MultiLALDetector detectors = {
    .length = 1,
    .sites = { lalCachedDetectors[LAL_LLO_4K_DETECTOR] }
  };
  EphemerisData *edat = XLALInitBarycenter( TEST_DATA_DIR "earth00-19-DE405.dat.gz",
                                            TEST_DATA_DIR "sun00-19-DE405.dat.gz" );
  XLAL_CHECK( edat != NULL, XLAL_EFUNC );
  const double freq_max = 40.0;
  SuperskyMetrics *metrics = XLALComputeSuperskyMetrics( SUPERSKY_METRIC_TYPE, 1, &ref_time, &segments, freq_max, &detectors, NULL, DETMOTION_SPIN | DETMOTION_PTOLEORBIT, edat );
  XLAL_CHECK( metrics != NULL, XLAL_EFUNC );
  XLAL_CHECK( metrics->num_segments == segments.length, XLAL_EFAILED );

  // Project and rescale semicoherent metric to give equal frequency spacings
  const double coh_max_mismatch = 1.4, semi_max_mismatch = 4.9;
  XLAL_CHECK( XLALEqualizeReducedSuperskyMetricsFreqSpacing( metrics, coh_max_mismatch, semi_max_mismatch ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Create lattice tilings
  LatticeTiling *coh_tiling[metrics->num_segments];
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    coh_tiling[n] = XLALCreateLatticeTiling( 4 );
    XLAL_CHECK( coh_tiling[n] != NULL, XLAL_EFUNC );
  }
  LatticeTiling *semi_tiling = XLALCreateLatticeTiling( 4 );
  XLAL_CHECK( semi_tiling != NULL, XLAL_EFUNC );

  // Add bounds
  const double alpha1 = 0, alpha2 = LAL_PI, delta1 = -LAL_PI_2, delta2 = LAL_PI_2;
  const double freq_min = freq_max - 0.05, f1dot = -3e-9;
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    XLAL_CHECK( XLALSetSuperskyPhysicalSkyBounds( coh_tiling[n], metrics->coh_rssky_metric[n], metrics->coh_rssky_transf[n], alpha1, alpha2, delta1, delta2 ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALSetSuperskyPhysicalSpinBound( coh_tiling[n], metrics->coh_rssky_transf[n], 0, freq_min, freq_max ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALSetSuperskyPhysicalSpinBound( coh_tiling[n], metrics->coh_rssky_transf[n], 1, f1dot, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  XLAL_CHECK( XLALSetSuperskyPhysicalSkyBounds( semi_tiling, metrics->semi_rssky_metric, metrics->semi_rssky_transf, alpha1, alpha2, delta1, delta2 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetSuperskyPhysicalSpinBound( semi_tiling, metrics->semi_rssky_transf, 0, freq_min, freq_max ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetSuperskyPhysicalSpinBound( semi_tiling, metrics->semi_rssky_transf, 1, f1dot, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Set metric
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    XLAL_CHECK( XLALSetTilingLatticeAndMetric( coh_tiling[n], TILING_LATTICE_ANSTAR, metrics->coh_rssky_metric[n], coh_max_mismatch ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  XLAL_CHECK( XLALSetTilingLatticeAndMetric( semi_tiling, TILING_LATTICE_ANSTAR, metrics->semi_rssky_metric, semi_max_mismatch ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Check lattice step sizes in frequency
  const size_t ifreq = 3;
  const double semi_dfreq = XLALLatticeTilingStepSize( semi_tiling, ifreq );
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    const double coh_dfreq = XLALLatticeTilingStepSize( coh_tiling[n], ifreq );
    const double tol = 1e-8;
    XLAL_CHECK( fabs( coh_dfreq - semi_dfreq ) < tol * semi_dfreq, XLAL_EFAILED, "semi_dfreq=%0.15e, coh_dfreq[%zu]=%0.15e, |coh_dfreq - semi_dfreq| >= %0.5g * semi_dfreq", semi_dfreq, n, coh_dfreq, tol );
  }

  // Print information on bounds
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    for ( size_t i = 0; i < XLALTotalLatticeTilingDimensions( coh_tiling[n] ); ++i ) {
      const LatticeTilingStats *stats = XLALLatticeTilingStatistics( coh_tiling[n], i );
      XLAL_CHECK( stats != NULL, XLAL_EFUNC );
      XLAL_CHECK( stats->name != NULL, XLAL_EFUNC );
      printf( "Coherent #%zu  bound #%zu: name=%6s, points=[%4u,%4u]\n", n, i, stats->name, stats->min_points, stats->max_points );
    }
  }
  for ( size_t i = 0; i < XLALTotalLatticeTilingDimensions( semi_tiling ); ++i ) {
    const LatticeTilingStats *stats = XLALLatticeTilingStatistics( semi_tiling, i );
    XLAL_CHECK( stats != NULL, XLAL_EFUNC );
    XLAL_CHECK( stats->name != NULL, XLAL_EFUNC );
    printf( "Semicoherent bound #%zu: name=%6s, points=[%4u,%4u]\n", i, stats->name, stats->min_points, stats->max_points );
  }
  printf( "\n" );

  // Check computation of spindown range for coherent tilings
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    PulsarDopplerParams min_range, max_range;
    XLAL_CHECK( XLALSuperskyLatticePhysicalRange( &min_range, &max_range, coh_tiling[n], metrics->coh_rssky_transf[n] ) == XLAL_SUCCESS, XLAL_EFUNC );
    printf( "Coherent #%zu spindown range: freq=[%0.5g,%0.5g], f1dot=[%0.5g,%0.5g]\n", n, min_range.fkdot[0], max_range.fkdot[0], min_range.fkdot[1], max_range.fkdot[1] );
  }
  printf( "\n" );

  // Perform mismatch test of coherent and semicoherent tilings
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    printf( "Coherent #%zu mismatch tests:\n", n );
    XLAL_CHECK( MismatchTest( coh_tiling[n], metrics->coh_rssky_metric[n], coh_max_mismatch, 0.1, 5e-2, 5e-2, coh_total_ref[n], total_tol, A4s_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );
    printf( "\n" );
  }
  printf( "Semicoherent mismatch tests:\n" );
  XLAL_CHECK( MismatchTest( semi_tiling, metrics->semi_rssky_metric, semi_max_mismatch, 0.0001, 5e-2, 5e-2, semi_total_ref, total_tol, A4s_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );
  printf( "\n" );

  // Cleanup
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    XLALDestroyLatticeTiling( coh_tiling[n] );
  }
  XLALDestroyLatticeTiling( semi_tiling );
  XLALDestroySuperskyMetrics( metrics );
  XLALSegListClear( &segments );
  XLALDestroyEphemerisData( edat );
  LALCheckMemoryLeaks();

  printf( "Finished super-sky metric tests\n" );

  return XLAL_SUCCESS;

}

int main( void )
{

  // Turn off buffering to sync standard output and error printing
  setvbuf( stdout, NULL, _IONBF, 0 );
  setvbuf( stderr, NULL, _IONBF, 0 );

  // Perform basic tests
  XLAL_CHECK_MAIN( BasicTest( 1, 0, 0, 0, 0, TILING_LATTICE_CUBIC,     1,    1,    1,    1 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 1, 1, 1, 1, 1, TILING_LATTICE_ANSTAR,   93,    0,    0,    0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 1, 1, 1, 1, 1, TILING_LATTICE_CUBIC,    93,    0,    0,    0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 2, 0, 0, 0, 0, TILING_LATTICE_ANSTAR,    1,    1,    1,    1 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 2, 1, 1, 1, 1, TILING_LATTICE_ANSTAR,   12,  144,    0,    0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 2, 1, 1, 1, 1, TILING_LATTICE_CUBIC,    13,  190,    0,    0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 3, 0, 0, 0, 0, TILING_LATTICE_CUBIC,     1,    1,    1,    1 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 3, 1, 1, 1, 1, TILING_LATTICE_ANSTAR,    8,   46,  332,    0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 3, 1, 1, 1, 1, TILING_LATTICE_CUBIC,     8,   60,  583,    0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 0, 0, 0, 0, TILING_LATTICE_ANSTAR,    1,    1,    1,    1 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 0, 0, 0, 1, TILING_LATTICE_ANSTAR,    1,    1,    1,    4 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 0, 0, 1, 0, TILING_LATTICE_ANSTAR,    1,    1,    4,    4 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 0, 0, 1, 1, TILING_LATTICE_ANSTAR,    1,    1,    4,   20 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 0, 1, 0, 0, TILING_LATTICE_ANSTAR,    1,    4,    4,    4 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 0, 1, 0, 1, TILING_LATTICE_ANSTAR,    1,    5,    5,   25 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 0, 1, 1, 0, TILING_LATTICE_ANSTAR,    1,    5,   24,   24 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 0, 1, 1, 1, TILING_LATTICE_ANSTAR,    1,    5,   20,  115 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 1, 0, 0, 0, TILING_LATTICE_ANSTAR,    4,    4,    4,    4 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 1, 0, 0, 1, TILING_LATTICE_ANSTAR,    5,    5,    5,   23 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 1, 0, 1, 0, TILING_LATTICE_ANSTAR,    5,    5,   23,   23 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 1, 0, 1, 1, TILING_LATTICE_ANSTAR,    6,    6,   24,  139 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 1, 1, 0, 0, TILING_LATTICE_ANSTAR,    5,   25,   25,   25 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 1, 1, 0, 1, TILING_LATTICE_ANSTAR,    6,   30,   30,  162 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 1, 1, 1, 0, TILING_LATTICE_ANSTAR,    6,   27,  151,  151 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 1, 1, 1, 1, TILING_LATTICE_ANSTAR,    6,   30,  145,  897 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( BasicTest( 4, 1, 1, 1, 1, TILING_LATTICE_CUBIC,     7,   46,  287, 2543 ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Perform mismatch tests with a square parameter space
  XLAL_CHECK_MAIN( MismatchSquareTest( TILING_LATTICE_CUBIC,  0.03,     0,     0, 21460,  Z1_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( MismatchSquareTest( TILING_LATTICE_CUBIC,  2e-4, -2e-9,     0, 23763,  Z2_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( MismatchSquareTest( TILING_LATTICE_CUBIC,  1e-4, -1e-9, 1e-17, 19550,  Z3_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( MismatchSquareTest( TILING_LATTICE_ANSTAR, 0.03,     0,     0, 21460, A1s_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( MismatchSquareTest( TILING_LATTICE_ANSTAR, 2e-4, -2e-9,     0, 18283, A2s_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( MismatchSquareTest( TILING_LATTICE_ANSTAR, 1e-4, -2e-9, 2e-17, 20268, A3s_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Perform mismatch tests with an age--braking index parameter space
  XLAL_CHECK_MAIN( MismatchAgeBrakeTest( TILING_LATTICE_ANSTAR, 100, 4.0e-5, 37868, A3s_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( MismatchAgeBrakeTest( TILING_LATTICE_ANSTAR, 200, 1.5e-5, 37230, A3s_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_MAIN( MismatchAgeBrakeTest( TILING_LATTICE_ANSTAR, 300, 1.0e-5, 37022, A3s_mism_hist ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Perform a variety of tests with the reduced supersky parameter space and metric
  XLAL_CHECK_MAIN( SuperskyTests( 5241516, 770802, 527943, 24586791553 ) == XLAL_SUCCESS, XLAL_EFUNC );

  return EXIT_SUCCESS;

}

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
