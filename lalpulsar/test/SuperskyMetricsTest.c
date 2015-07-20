//
// Copyright (C) 2015 Karl Wette
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

// Tests of the supersky metric code in SuperskyMetrics.[ch].

#include <math.h>
#include <gsl/gsl_blas.h>

#include <lal/SuperskyMetrics.h>
#include <lal/LALStdlib.h>
#include <lal/LALInitBarycenter.h>
#include <lal/MetricUtils.h>

#include "../src/GSLHelpers.h"

#define NUM_POINTS 10

const PulsarDopplerParams phys_points[NUM_POINTS] = {
  { .Alpha = 0.00000000000000, .Delta =  0.000000000000000, .fkdot = {100.0000000000000,  0.00000000000000e-00} },
  { .Alpha = 4.88014010120016, .Delta = -0.954446475246007, .fkdot = { 99.9999983978492,  1.48957780038094e-09} },
  { .Alpha = 0.52587274931672, .Delta =  0.685297319257976, .fkdot = { 99.9999923150006, -1.41365319702693e-09} },
  { .Alpha = 3.53542175437611, .Delta = -1.502778038590950, .fkdot = {100.0000064863180, -1.28748375084384e-09} },
  { .Alpha = 1.36054903961191, .Delta =  0.241343663657163, .fkdot = { 99.9999901679571,  3.37107171004537e-10} },
  { .Alpha = 2.85470536965808, .Delta = 1.1575340928032900, .fkdot = {100.0000074463050,  2.46412240438217e-09} },
  { .Alpha = 1.82755817952460, .Delta =  0.667995269285982, .fkdot = { 99.9999897239871,  1.79900370692270e-10} },
  { .Alpha = 1.70734223243163, .Delta = -1.213787405673430, .fkdot = {100.0000026535270, -1.07122135891104e-09} },
  { .Alpha = 2.30597131157246, .Delta =  0.348657791621429, .fkdot = {100.0000133749770, -5.43309003215614e-10} },
  { .Alpha = 3.31129323970275, .Delta = -1.225892709583030, .fkdot = {100.0000062524320,  8.07713885739405e-10} }
};

const double phys_mismatches[NUM_POINTS][NUM_POINTS] = {
  {  0.0000000000e+00,  2.8764538705e+07,  1.7369679160e+05,  2.8451027863e+07,  3.4148114879e+06,  2.7431945400e+07,  1.3229228239e+07,  2.1651786947e+07,  3.3466651446e+07,  4.2532446390e+07 },
  {  2.8764538705e+07,  0.0000000000e+00,  2.4471879340e+07,  1.7793534479e+03,  1.2362486311e+07,  1.9460081514e+04,  2.9815077312e+06,  5.0457978267e+05,  1.7846091686e+05,  1.3420016910e+06 },
  {  1.7369679160e+05,  2.4471879340e+07,  0.0000000000e+00,  2.4182611494e+07,  2.0496514361e+06,  2.3244047105e+07,  1.0373862422e+07,  1.7950270301e+07,  2.8822526810e+07,  3.7274957328e+07 },
  {  2.8451027863e+07,  1.7793534479e+03,  2.4182611494e+07,  0.0000000000e+00,  1.2159740027e+07,  1.7226320532e+04,  2.8838705301e+06,  4.6363967921e+05,  2.0619676630e+05,  1.4120605271e+06 },
  {  3.4148114879e+06,  1.2362486311e+07,  2.0496514361e+06,  1.2159740027e+07,  0.0000000000e+00,  1.1490795829e+07,  3.2024152317e+06,  7.8746153787e+06,  1.5504363899e+07,  2.1849636064e+07 },
  {  2.7431945400e+07,  1.9460081514e+04,  2.3244047105e+07,  1.7226320532e+04,  1.1490795829e+07,  0.0000000000e+00,  2.5611895406e+06,  3.4663900394e+05,  3.0146103485e+05,  1.6529327758e+06 },
  {  1.3229228239e+07,  2.9815077312e+06,  1.0373862422e+07,  2.8838705301e+06,  3.2024152317e+06,  2.5611895406e+06,  0.0000000000e+00,  1.0352316575e+06,  4.6141065584e+06,  8.3226021460e+06 },
  {  2.1651786947e+07,  5.0457978267e+05,  1.7950270301e+07,  4.6363967921e+05,  7.8746153787e+06,  3.4663900394e+05,  1.0352316575e+06,  0.0000000000e+00,  1.2824844254e+06,  3.4919141748e+06 },
  {  3.3466651446e+07,  1.7846091686e+05,  2.8822526810e+07,  2.0619676630e+05,  1.5504363899e+07,  3.0146103485e+05,  4.6141065584e+06,  1.2824844254e+06,  0.0000000000e+00,  5.4308661488e+05 },
  {  4.2532446390e+07,  1.3420016910e+06,  3.7274957328e+07,  1.4120605271e+06,  2.1849636064e+07,  1.6529327758e+06,  8.3226021460e+06,  3.4919141748e+06,  5.4308661488e+05,  0.0000000000e+00 },
};

const double ussky_metric_ref[5][5] = {
  { 1.895617224814336e+07,  6.731809404038861e+06,  2.921605097922917e+06,  2.046904385910545e+09, -1.656990342031264e+11},
  { 6.731809404038861e+06,  2.391591549452773e+06,  1.037930347967158e+06,  7.269432686412601e+08,  4.092485067498952e+11},
  { 2.921605097922917e+06,  1.037930347967158e+06,  4.504649325079115e+05,  3.154935665162594e+08,  1.815805951797570e+11},
  { 2.046904385910545e+09,  7.269432686412601e+08,  3.154935665162594e+08,  2.210286062098680e+11, -5.119042942971644e-01},
  {-1.656990342031264e+11,  4.092485067498952e+11,  1.815805951797570e+11, -5.119042942971644e-01,  2.474954556318625e+20},
};

const double rssky_metric_ref[4][4] = {
  { 6.568666765075774e+01,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00},
  { 0.000000000000000e+00,  6.236546718525190e+01,  0.000000000000000e+00,  0.000000000000000e+00},
  { 0.000000000000000e+00,  0.000000000000000e+00,  2.474954556329869e+20, -1.208418861582654e+02},
  { 0.000000000000000e+00,  0.000000000000000e+00, -1.208418861582654e+02,  2.210286062098680e+11},
};

const double rssky_transf_ref[5][3] = {
  { 7.928332170641094e-01, -6.094386653165459e-01,  5.600859050599027e-05},
  { 6.094384042541772e-01,  7.928329562280915e-01,  8.572856860540078e-04},
  {-5.668685006887667e-04, -6.455507823947125e-04,  9.999996309620771e-01},
  {-1.538506600219701e-09,  9.036045301253465e-10,  7.329842343015285e-10},
  { 5.337970249571395e-03,  8.252674722491856e-03,  1.420014590206618e-03},
};

#define CHECK_RELERR(A, B, TOL) do { \
    const double lhs = fabs( (A) - (B) ); \
    const double tol = (TOL); \
    const double rhs = GSL_MAX( 1.0, fabs( (A) + (B) ) ); \
    XLALPrintInfo( #A"=%0.5e   "#B"=%0.5e   |"#A" - "#B"|=%0.5e   tol=%0.5e   |"#A" + "#B"|=%0.5e\n", A, B, lhs, tol, rhs ); \
    XLAL_CHECK_MAIN( lhs <= tol * rhs, XLAL_ETOL, "|"#A" - "#B"| = %0.5e > %0.5e = %0.5e * |"#A" + "#B"|", lhs, tol * rhs, tol ); \
  } while(0)

static int CompareDoppler(const PulsarDopplerParams *a, const PulsarDopplerParams *b) {
  CHECK_RELERR( cos(a->Alpha), cos(b->Alpha), 1e-10 );
  CHECK_RELERR( sin(a->Alpha), sin(b->Alpha), 1e-10 );
  CHECK_RELERR( a->Delta, b->Delta, 1e-10 );
  CHECK_RELERR( a->fkdot[0], b->fkdot[0], 1e-10 );
  CHECK_RELERR( a->fkdot[1], b->fkdot[1], 1e-10 );
  return EXIT_SUCCESS;
}

int main( void )
{

  // Compute supersky metrics
  const double Tspan = 3 * 86400;
  LIGOTimeGPS ref_time;
  XLALGPSSetREAL8( &ref_time, 900100100 );
  LALSegList segments;
  {
    XLAL_CHECK_MAIN( XLALSegListInit( &segments ) == XLAL_SUCCESS, XLAL_EFUNC );
    LALSeg segment;
    LIGOTimeGPS start_time = ref_time, end_time = ref_time;
    XLALGPSAdd( &start_time, -0.5 * Tspan );
    XLALGPSAdd( &end_time, 0.5 * Tspan );
    XLAL_CHECK_MAIN( XLALSegSet( &segment, &start_time, &end_time, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_MAIN( XLALSegListAppend( &segments, &segment ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  MultiLALDetector detectors = {
    .length = 1,
    .sites = { lalCachedDetectors[LAL_LLO_4K_DETECTOR] }
  };
  EphemerisData *edat =  XLALInitBarycenter( TEST_DATA_DIR "earth00-19-DE405.dat.gz",
                                             TEST_DATA_DIR "sun00-19-DE405.dat.gz" );
  XLAL_CHECK_MAIN( edat != NULL, XLAL_EFUNC );
  gsl_matrix *rssky_metric = NULL, *rssky_transf = NULL, *ussky_metric = NULL;
  XLAL_CHECK_MAIN( XLALComputeSuperskyMetrics( &rssky_metric, &rssky_transf, &ussky_metric, 1, &ref_time, &segments, 100.0, &detectors, NULL, DETMOTION_SPIN | DETMOTION_PTOLEORBIT, edat ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALSegListClear( &segments );
  XLALDestroyEphemerisData( edat );

  // Check supersky metrics
  {
    gsl_matrix_const_view ussky_metric_ref_view = gsl_matrix_const_view_array( (const double*)ussky_metric_ref, 5, 5 );
    const double err = XLALCompareMetrics( ussky_metric, &ussky_metric_ref_view.matrix ), err_tol = 5e-9;
    XLAL_CHECK_MAIN( err <= err_tol, XLAL_ETOL, "'ussky_metric' check failed: err = %0.3e > %0.3e = err_tol", err, err_tol );
  }
  {
    gsl_matrix_const_view rssky_metric_ref_view = gsl_matrix_const_view_array( (const double*)rssky_metric_ref, 4, 4 );
    const double err = XLALCompareMetrics( rssky_metric, &rssky_metric_ref_view.matrix ), err_tol = 5e-9;
    XLAL_CHECK_MAIN( err <= err_tol, XLAL_ETOL, "'rssky_metric' check failed: err = %0.3e > %0.3e = err_tol", err, err_tol );
  }
  {
    double max_err = 0;
    for( size_t i = 0; i < 5; ++i ) {
      for( size_t j = 0; j < 3; ++j ) {
        const double rssky_transf_ij = gsl_matrix_get( rssky_transf, i, j );
        const double rssky_transf_ref_ij = rssky_transf_ref[i][j];
        const double err_ij = fabs( ( rssky_transf_ij - rssky_transf_ref_ij ) / rssky_transf_ref_ij );
        if( err_ij > max_err ) {
          max_err = err_ij;
        }
      }
    }
    const double err_tol = 1e-6;
    XLAL_CHECK_MAIN( max_err <= err_tol, XLAL_ETOL, "'rssky_transf' check failed: max(err) = %0.3e > %0.3e = err_tol", max_err, err_tol );
  }

  // Check round-trip conversions of each test point
  {
    gsl_vector *GAVEC( ussky_point, 5 );
    for( size_t i = 0; i < NUM_POINTS; ++i ) {
      PulsarDopplerParams XLAL_INIT_DECL(point);
      XLAL_CHECK_MAIN( XLALConvertPhysicalToSupersky( SC_USSKY, ussky_point, &phys_points[i], NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLALConvertSuperskyToPhysical( &point, SC_USSKY, ussky_point, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( CompareDoppler( &phys_points[i], &point ) == EXIT_SUCCESS, XLAL_EFUNC );
    }
    GFVEC( ussky_point );
  }
  {
    gsl_vector *GAVEC( rssky_point, 4 );
    for( size_t i = 0; i < NUM_POINTS; ++i ) {
      PulsarDopplerParams XLAL_INIT_DECL(point);
      XLAL_CHECK_MAIN( XLALConvertPhysicalToSupersky( SC_RSSKY, rssky_point, &phys_points[i], rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLALConvertSuperskyToPhysical( &point, SC_RSSKY, rssky_point, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( CompareDoppler( &phys_points[i], &point ) == EXIT_SUCCESS, XLAL_EFUNC );
    }
    GFVEC( rssky_point );
  }

  // Check mismatches between pairs of points
  {
    gsl_vector *GAVEC( ussky_point_i, 5 );
    gsl_vector *GAVEC( ussky_point_j, 5 );
    gsl_vector *GAVEC( temp, 5 );
    for( size_t i = 0; i < NUM_POINTS; ++i ) {
      XLAL_CHECK_MAIN( XLALConvertPhysicalToSupersky( SC_USSKY, ussky_point_i, &phys_points[i], NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      for( size_t j = 0; j < NUM_POINTS; ++j ) {
        XLAL_CHECK_MAIN( XLALConvertPhysicalToSupersky( SC_USSKY, ussky_point_j, &phys_points[j], NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
        gsl_vector_sub( ussky_point_j, ussky_point_i );
        gsl_blas_dgemv( CblasNoTrans, 1.0, ussky_metric, ussky_point_j, 0.0, temp );
        double mismatch = 0.0;
        gsl_blas_ddot( ussky_point_j, temp, &mismatch );
        CHECK_RELERR( mismatch, phys_mismatches[i][j], 8e-3 );
      }
    }
    GFVEC( ussky_point_i, ussky_point_j, temp );
  }
  {
    gsl_vector *GAVEC( rssky_point_i, 4 );
    gsl_vector *GAVEC( rssky_point_j, 4 );
    gsl_vector *GAVEC( temp, 4 );
    for( size_t i = 0; i < NUM_POINTS; ++i ) {
      XLAL_CHECK_MAIN( XLALConvertPhysicalToSupersky( SC_RSSKY, rssky_point_i, &phys_points[i], rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
      for( size_t j = 0; j < NUM_POINTS; ++j ) {
        XLAL_CHECK_MAIN( XLALConvertPhysicalToSupersky( SC_RSSKY, rssky_point_j, &phys_points[j], rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
        gsl_vector_sub( rssky_point_j, rssky_point_i );
        gsl_blas_dgemv( CblasNoTrans, 1.0, rssky_metric, rssky_point_j, 0.0, temp );
        double mismatch = 0.0;
        gsl_blas_ddot( rssky_point_j, temp, &mismatch );
        CHECK_RELERR( mismatch, phys_mismatches[i][j], 8e-3 );
      }
    }
    GFVEC( rssky_point_i, rssky_point_j, temp );
  }

  // Cleanup
  GFMAT( rssky_metric, rssky_transf, ussky_metric );
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
