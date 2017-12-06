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

// Tests of the supersky metric code in SuperSkyMetrics.[ch].

#include <math.h>
#include <gsl/gsl_blas.h>

#include <lal/SuperSkyMetrics.h>
#include <lal/LALStdlib.h>
#include <lal/LALInitBarycenter.h>

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

const double ssky_mismatches[NUM_POINTS][NUM_POINTS] = {
  {  0.0000000000e+00,  2.8764539641e+07,  1.7339193295e+05,  2.8451029480e+07,  3.4146745049e+06,  2.7431796599e+07,  1.3229138543e+07,  2.1651788499e+07,  3.3466633515e+07,  4.2532448208e+07 },
  {  2.8764539641e+07,  0.0000000000e+00,  2.4471638891e+07,  1.7793497533e+03,  1.2362377917e+07,  1.9343156763e+04,  2.9814372830e+06,  5.0457966981e+05,  1.7844672448e+05,  1.3420018117e+06 },
  {  1.7339193295e+05,  2.4471638891e+07,  0.0000000000e+00,  2.4182449487e+07,  2.0496517155e+06,  2.3244046454e+07,  1.0373862571e+07,  1.7950144640e+07,  2.8822527465e+07,  3.7274828956e+07 },
  {  2.8451029480e+07,  1.7793497533e+03,  2.4182449487e+07,  0.0000000000e+00,  1.2159667072e+07,  1.7147004362e+04,  2.8838230958e+06,  4.6363965844e+05,  2.0618704546e+05,  1.4120604982e+06 },
  {  3.4146745049e+06,  1.2362377917e+07,  2.0496517155e+06,  1.2159667072e+07,  0.0000000000e+00,  1.1490794709e+07,  3.2024149649e+06,  7.8745587587e+06,  1.5504363611e+07,  2.1849578253e+07 },
  {  2.7431796599e+07,  1.9343156763e+04,  2.3244046454e+07,  1.7147004362e+04,  1.1490794709e+07,  0.0000000000e+00,  2.5611892506e+06,  3.4657692079e+05,  3.0146117597e+05,  1.6528701940e+06 },
  {  1.3229138543e+07,  2.9814372830e+06,  1.0373862571e+07,  2.8838230958e+06,  3.2024149649e+06,  2.5611892506e+06,  0.0000000000e+00,  1.0351946794e+06,  4.6141067214e+06,  8.3225648673e+06 },
  {  2.1651788499e+07,  5.0457966981e+05,  1.7950144640e+07,  4.6363965844e+05,  7.8745587587e+06,  3.4657692079e+05,  1.0351946794e+06,  0.0000000000e+00,  1.2824765725e+06,  3.4919140724e+06 },
  {  3.3466633515e+07,  1.7844672448e+05,  2.8822527465e+07,  2.0618704546e+05,  1.5504363611e+07,  3.0146117597e+05,  4.6141067214e+06,  1.2824765725e+06,  0.0000000000e+00,  5.4307920777e+05 },
  {  4.2532448208e+07,  1.3420018117e+06,  3.7274828956e+07,  1.4120604982e+06,  2.1849578253e+07,  1.6528701940e+06,  8.3225648673e+06,  3.4919140724e+06,  5.4307920777e+05,  0.0000000000e+00 }
};

const double rssky_mismatches[NUM_POINTS][NUM_POINTS] = {
  {  0.0000000000e+00,  1.6479386809e+20,  9.9138389085e+19,  2.4635218573e+20,  1.4138082659e+19,  2.0757894095e+20,  9.4957310660e+19,  2.1726843923e+20,  2.8886583299e+19,  2.1920285214e+20 },
  {  1.6479386809e+20,  0.0000000000e+00,  5.1956792459e+20,  8.1704734243e+18,  2.7546939071e+20,  7.4227948167e+20,  5.0993816904e+20,  3.6207318308e+18,  3.3167076696e+20,  3.8741826153e+18 },
  {  9.9138389085e+19,  5.1956792459e+20,  0.0000000000e+00,  6.5804750622e+20,  3.8399860188e+19,  1.9809387909e+19,  4.5038208620e+16,  6.0993463322e+20,  2.0996577377e+19,  6.1317283819e+20 },
  {  2.4635218573e+20,  8.1704734243e+18,  6.5804750622e+20,  0.0000000000e+00,  3.7852327335e+20,  9.0620328287e+20,  6.4720451781e+20,  9.1313682359e+17,  4.4395477816e+20,  7.9229520300e+17 },
  {  1.4138082659e+19,  2.7546939071e+20,  3.8399860188e+19,  3.7852327335e+20,  0.0000000000e+00,  1.1337002134e+20,  3.5814719000e+19,  3.4225342808e+20,  2.6068034507e+18,  3.4468020087e+20 },
  {  2.0757894095e+20,  7.4227948167e+20,  1.9809387909e+19,  9.0620328287e+20,  1.1337002134e+20,  0.0000000000e+00,  2.1743531033e+19,  8.4958425175e+20,  8.1594657153e+19,  8.5340526127e+20 },
  {  9.4957310660e+19,  5.0993816904e+20,  4.5038208620e+16,  6.4720451781e+20,  3.5814719000e+19,  2.1743531033e+19,  0.0000000000e+00,  5.9949723598e+20,  1.9096726661e+19,  6.0270765163e+20 },
  {  2.1726843923e+20,  3.6207318308e+18,  6.0993463322e+20,  9.1313682359e+17,  3.4225342808e+20,  8.4958425175e+20,  5.9949723598e+20,  0.0000000000e+00,  4.0459924448e+20,  4.2866183279e+15 },
  {  2.8886583299e+19,  3.3167076696e+20,  2.0996577377e+19,  4.4395477816e+20,  2.6068034507e+18,  8.1594657153e+19,  1.9096726661e+19,  4.0459924448e+20,  0.0000000000e+00,  4.0723743508e+20 },
  {  2.1920285214e+20,  3.8741826153e+18,  6.1317283819e+20,  7.9229520300e+17,  3.4468020087e+20,  8.5340526127e+20,  6.0270765163e+20,  4.2866183279e+15,  4.0723743508e+20,  0.0000000000e+00 }
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

  // Compute super-sky metrics
  const double Tspan = 3 * 86400;
  gsl_matrix *essky_metric = NULL;
  LIGOTimeGPS ref_time;
  XLALGPSSetREAL8( &ref_time, 900100100 );
  LALSegList segments;
  {
    XLAL_CHECK( XLALSegListInit( &segments ) == XLAL_SUCCESS, XLAL_EFUNC );
    LALSeg segment;
    LIGOTimeGPS start_time = ref_time, end_time = ref_time;
    XLALGPSAdd( &start_time, -0.5 * Tspan );
    XLALGPSAdd( &end_time, 0.5 * Tspan );
    XLAL_CHECK( XLALSegSet( &segment, &start_time, &end_time, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALSegListAppend( &segments, &segment ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  MultiLALDetector detectors = {
    .length = 1,
    .sites = { lalCachedDetectors[LAL_LLO_4K_DETECTOR] }
  };
  EphemerisData *edat =  XLALInitBarycenter( TEST_DATA_DIR "earth00-19-DE405.dat.gz",
                                             TEST_DATA_DIR "sun00-19-DE405.dat.gz" );
  XLAL_CHECK( edat != NULL, XLAL_EFUNC );
  XLAL_CHECK( XLALExpandedSuperSkyMetric( &essky_metric, 1, &ref_time, &segments, 100.0, &detectors, NULL, DETMOTION_SPIN | DETMOTION_PTOLEORBIT, edat ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALSegListClear( &segments );
  XLALDestroyEphemerisData( edat );
  gsl_matrix *ssky_metric = NULL, *rssky_metric = NULL, *rssky_transf = NULL;
  XLAL_CHECK( XLALSuperSkyMetric( &ssky_metric, essky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALReducedSuperSkyMetric( &rssky_metric, &rssky_transf, essky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );
  GFMAT( essky_metric );

  // Check round-trip conversions of each test point
  {
    gsl_vector *GAVEC( ssky_point, 5 );
    for( size_t i = 0; i < NUM_POINTS; ++i ) {
      PulsarDopplerParams XLAL_INIT_DECL(point);
      XLAL_CHECK_MAIN( XLALConvertPhysicalToSuperSky( SSC_SUPER_SKY, ssky_point, &phys_points[i], NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLALConvertSuperSkyToPhysical( &point, SSC_SUPER_SKY, ssky_point, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( CompareDoppler( &phys_points[i], &point ) == EXIT_SUCCESS, XLAL_EFUNC );
    }
    GFVEC( ssky_point );
  }
  {
    gsl_vector *GAVEC( rssky_point, 4 );
    for( size_t i = 0; i < NUM_POINTS; ++i ) {
      PulsarDopplerParams XLAL_INIT_DECL(point);
      XLAL_CHECK_MAIN( XLALConvertPhysicalToSuperSky( SSC_REDUCED_SUPER_SKY, rssky_point, &phys_points[i], rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK_MAIN( XLALConvertSuperSkyToPhysical( &point, SSC_REDUCED_SUPER_SKY, rssky_point, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( CompareDoppler( &phys_points[i], &point ) == EXIT_SUCCESS, XLAL_EFUNC );
    }
    GFVEC( rssky_point );
  }

  // Check mismatches between pairs of points
  {
    gsl_vector *GAVEC( ssky_point_i, 5 );
    gsl_vector *GAVEC( ssky_point_j, 5 );
    gsl_vector *GAVEC( temp, 5 );
    for( size_t i = 0; i < NUM_POINTS; ++i ) {
      XLAL_CHECK_MAIN( XLALConvertPhysicalToSuperSky( SSC_SUPER_SKY, ssky_point_i, &phys_points[i], NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
      for( size_t j = 0; j < NUM_POINTS; ++j ) {
        XLAL_CHECK_MAIN( XLALConvertPhysicalToSuperSky( SSC_SUPER_SKY, ssky_point_j, &phys_points[j], NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
        gsl_vector_sub( ssky_point_j, ssky_point_i );
        gsl_blas_dgemv( CblasNoTrans, 1.0, ssky_metric, ssky_point_j, 0.0, temp );
        double mismatch = 0.0;
        gsl_blas_ddot( ssky_point_j, temp, &mismatch );
        CHECK_RELERR( mismatch, ssky_mismatches[i][j], 1e-5 );
      }
    }
    GFVEC( ssky_point_i, ssky_point_j, temp );
  }
  {
    gsl_vector *GAVEC( rssky_point_i, 4 );
    gsl_vector *GAVEC( rssky_point_j, 4 );
    gsl_vector *GAVEC( temp, 4 );
    for( size_t i = 0; i < NUM_POINTS; ++i ) {
      XLAL_CHECK_MAIN( XLALConvertPhysicalToSuperSky( SSC_SUPER_SKY, rssky_point_i, &phys_points[i], rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
      for( size_t j = 0; j < NUM_POINTS; ++j ) {
        XLAL_CHECK_MAIN( XLALConvertPhysicalToSuperSky( SSC_SUPER_SKY, rssky_point_j, &phys_points[j], rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
        gsl_vector_sub( rssky_point_j, rssky_point_i );
        gsl_blas_dgemv( CblasNoTrans, 1.0, rssky_metric, rssky_point_j, 0.0, temp );
        double mismatch = 0.0;
        gsl_blas_ddot( rssky_point_j, temp, &mismatch );
        CHECK_RELERR( mismatch, rssky_mismatches[i][j], 1e-5 );
      }
    }
    GFVEC( rssky_point_i, rssky_point_j, temp );
  }

  // Cleanup
  GFMAT( ssky_metric, rssky_metric, rssky_transf );
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
