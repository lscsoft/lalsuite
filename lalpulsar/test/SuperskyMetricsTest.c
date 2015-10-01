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

#define REF_TIME { 900100100, 0 }
const LIGOTimeGPS ref_time = REF_TIME;

const PulsarDopplerParams phys_points[NUM_POINTS] = {
  { .refTime = REF_TIME, .Alpha = 0.00000000000000, .Delta =  0.000000000000000, .fkdot = {100.0000000000000,  0.00000000000000e-00} },
  { .refTime = REF_TIME, .Alpha = 4.88014010120016, .Delta = -0.954446475246007, .fkdot = { 99.9999983978492,  1.48957780038094e-09} },
  { .refTime = REF_TIME, .Alpha = 0.52587274931672, .Delta =  0.685297319257976, .fkdot = { 99.9999923150006, -1.41365319702693e-09} },
  { .refTime = REF_TIME, .Alpha = 3.53542175437611, .Delta = -1.502778038590950, .fkdot = {100.0000064863180, -1.28748375084384e-09} },
  { .refTime = REF_TIME, .Alpha = 1.36054903961191, .Delta =  0.241343663657163, .fkdot = { 99.9999901679571,  3.37107171004537e-10} },
  { .refTime = REF_TIME, .Alpha = 2.85470536965808, .Delta = 1.1575340928032900, .fkdot = {100.0000074463050,  2.46412240438217e-09} },
  { .refTime = REF_TIME, .Alpha = 1.82755817952460, .Delta =  0.667995269285982, .fkdot = { 99.9999897239871,  1.79900370692270e-10} },
  { .refTime = REF_TIME, .Alpha = 1.70734223243163, .Delta = -1.213787405673430, .fkdot = {100.0000026535270, -1.07122135891104e-09} },
  { .refTime = REF_TIME, .Alpha = 2.30597131157246, .Delta =  0.348657791621429, .fkdot = {100.0000133749770, -5.43309003215614e-10} },
  { .refTime = REF_TIME, .Alpha = 3.31129323970275, .Delta = -1.225892709583030, .fkdot = {100.0000062524320,  8.07713885739405e-10} }
};

const double Tspan = 3 * 86400;

const double ussky_metric_refs[1][5][5] = {
  {
    { 1.895617224814336e+07,  6.731809404038861e+06,  2.921605097922917e+06,  2.046904385910545e+09, -1.656990342031264e+11},
    { 6.731809404038861e+06,  2.391591549452773e+06,  1.037930347967158e+06,  7.269432686412601e+08,  4.092485067498952e+11},
    { 2.921605097922917e+06,  1.037930347967158e+06,  4.504649325079115e+05,  3.154935665162594e+08,  1.815805951797570e+11},
    { 2.046904385910545e+09,  7.269432686412601e+08,  3.154935665162594e+08,  2.210286062098680e+11, -5.119042942971644e-01},
    {-1.656990342031264e+11,  4.092485067498952e+11,  1.815805951797570e+11, -5.119042942971644e-01,  2.474954556318625e+20},
  }
};

const double rssky_metric_refs[1][4][4] = {
  {
    { 6.568666765075774e+01,  0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00},
    { 0.000000000000000e+00,  6.236546718525190e+01,  0.000000000000000e+00,  0.000000000000000e+00},
    { 0.000000000000000e+00,  0.000000000000000e+00,  2.474954556329869e+20, -1.208418861582654e+02},
    { 0.000000000000000e+00,  0.000000000000000e+00, -1.208418861582654e+02,  2.210286062098680e+11},
  }
};

const double rssky_transf_refs[1][5][3] = {
  {
    { 7.928332170641094e-01, -6.094386653165459e-01,  5.600859050599027e-05},
    { 6.094384042541772e-01,  7.928329562280915e-01,  8.572856860540078e-04},
    {-5.668685006887667e-04, -6.455507823947125e-04,  9.999996309620771e-01},
    {-1.538506600219701e-09,  9.036045301253465e-10,  7.329842343015285e-10},
    { 5.337970249571395e-03,  8.252674722491856e-03,  1.420014590206618e-03},
  }
};

const double phys_mismatches[1][NUM_POINTS][NUM_POINTS] = { 
  {
    {0.000000e+00, 2.876453e+07, 1.736967e+05, 2.845102e+07, 3.414811e+06, 2.743194e+07, 1.322922e+07, 2.165178e+07, 3.346665e+07, 4.253244e+07},
    {2.876453e+07, 0.000000e+00, 2.447187e+07, 1.779353e+03, 1.236248e+07, 1.946008e+04, 2.981507e+06, 5.045797e+05, 1.784609e+05, 1.342001e+06},
    {1.736967e+05, 2.447187e+07, 0.000000e+00, 2.418261e+07, 2.049651e+06, 2.324404e+07, 1.037386e+07, 1.795027e+07, 2.882252e+07, 3.727495e+07},
    {2.845102e+07, 1.779353e+03, 2.418261e+07, 0.000000e+00, 1.215974e+07, 1.722632e+04, 2.883870e+06, 4.636396e+05, 2.061967e+05, 1.412060e+06},
    {3.414811e+06, 1.236248e+07, 2.049651e+06, 1.215974e+07, 0.000000e+00, 1.149079e+07, 3.202415e+06, 7.874615e+06, 1.550436e+07, 2.184963e+07},
    {2.743194e+07, 1.946008e+04, 2.324404e+07, 1.722632e+04, 1.149079e+07, 0.000000e+00, 2.561189e+06, 3.466390e+05, 3.014610e+05, 1.652932e+06},
    {1.322922e+07, 2.981507e+06, 1.037386e+07, 2.883870e+06, 3.202415e+06, 2.561189e+06, 0.000000e+00, 1.035231e+06, 4.614106e+06, 8.322602e+06},
    {2.165178e+07, 5.045797e+05, 1.795027e+07, 4.636396e+05, 7.874615e+06, 3.466390e+05, 1.035231e+06, 0.000000e+00, 1.282484e+06, 3.491914e+06},
    {3.346665e+07, 1.784609e+05, 2.882252e+07, 2.061967e+05, 1.550436e+07, 3.014610e+05, 4.614106e+06, 1.282484e+06, 0.000000e+00, 5.430866e+05},
    {4.253244e+07, 1.342001e+06, 3.727495e+07, 1.412060e+06, 2.184963e+07, 1.652932e+06, 8.322602e+06, 3.491914e+06, 5.430866e+05, 0.000000e+00},
  }
};

#define CHECK_RELERR(A, B, TOL) do { \
    const double lhs = fabs( (A) - (B) ); \
    const double tol = (TOL); \
    const double rhs = GSL_MAX( 1.0, fabs( (A) + (B) ) ); \
    XLALPrintInfo( #A"=%0.5e   "#B"=%0.5e   |"#A" - "#B"|=%0.5e   tol=%0.5e   |"#A" + "#B"|=%0.5e\n", A, B, lhs, tol, rhs ); \
    XLAL_CHECK( lhs <= tol * rhs, XLAL_ETOL, "|"#A" - "#B"| = %0.5e > %0.5e = %0.5e * |"#A" + "#B"|", lhs, tol * rhs, tol ); \
  } while(0)

static int CompareDoppler(const PulsarDopplerParams *a, const PulsarDopplerParams *b)
{
  XLAL_CHECK(XLALGPSCmp(&a->refTime, &b->refTime) == 0, XLAL_ETOL, "Reference time mismatch!");
  CHECK_RELERR(cos(a->Alpha), cos(b->Alpha), 1e-10);
  CHECK_RELERR(sin(a->Alpha), sin(b->Alpha), 1e-10);
  CHECK_RELERR(a->Delta, b->Delta, 1e-10);
  CHECK_RELERR(a->fkdot[0], b->fkdot[0], 1e-10);
  CHECK_RELERR(a->fkdot[1], b->fkdot[1], 1e-10);
  return XLAL_SUCCESS;
}

static int CheckSuperskyMetrics(
  const gsl_matrix *ussky_metric,
  const double ussky_metric_ref[5][5],
  const gsl_matrix *rssky_metric,
  const double rssky_metric_ref[4][4],
  const gsl_matrix *rssky_transf,
  const double rssky_transf_ref[5][3],
  const double phys_mismatch[NUM_POINTS][NUM_POINTS]
  )
{

  // Check supersky metrics
  {
    gsl_matrix_const_view ussky_metric_ref_view = gsl_matrix_const_view_array((const double *)ussky_metric_ref, 5, 5);
    const double err = XLALCompareMetrics(ussky_metric, &ussky_metric_ref_view.matrix), err_tol = 5e-9;
    XLAL_CHECK(err <= err_tol, XLAL_ETOL, "'ussky_metric' check failed: err = %0.3e > %0.3e = err_tol", err, err_tol);
  }
  {
    gsl_matrix_const_view rssky_metric_ref_view = gsl_matrix_const_view_array((const double *)rssky_metric_ref, 4, 4);
    const double err = XLALCompareMetrics(rssky_metric, &rssky_metric_ref_view.matrix), err_tol = 5e-9;
    XLAL_CHECK(err <= err_tol, XLAL_ETOL, "'rssky_metric' check failed: err = %0.3e > %0.3e = err_tol", err, err_tol);
  }
  {
    double max_err = 0;
    for (size_t i = 0; i < 5; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        const double rssky_transf_ij = gsl_matrix_get(rssky_transf, i, j);
        const double rssky_transf_ref_ij = rssky_transf_ref[i][j];
        const double err_ij = fabs((rssky_transf_ij - rssky_transf_ref_ij) / rssky_transf_ref_ij);
        if (err_ij > max_err) {
          max_err = err_ij;
        }
      }
    }
    const double err_tol = 1e-6;
    XLAL_CHECK(max_err <= err_tol, XLAL_ETOL, "'rssky_transf' check failed: max(err) = %0.3e > %0.3e = err_tol", max_err, err_tol);
  }

  // Check round-trip conversions of each test point
  {
    gsl_vector *GAVEC(ussky_point, 5);
    for (size_t i = 0; i < NUM_POINTS; ++i) {
      PulsarDopplerParams XLAL_INIT_DECL(point);
      XLAL_CHECK(XLALConvertPhysicalToSupersky(SC_USSKY, ussky_point, &phys_points[i], NULL, &ref_time) == XLAL_SUCCESS, XLAL_EFUNC);
      XLAL_CHECK(XLALConvertSuperskyToPhysical(&point, SC_USSKY, ussky_point, NULL, &ref_time) == XLAL_SUCCESS, XLAL_EFUNC);
      XLAL_CHECK(CompareDoppler(&phys_points[i], &point) == EXIT_SUCCESS, XLAL_EFUNC);
    }
    GFVEC(ussky_point);
  }
  {
    gsl_vector *GAVEC(rssky_point, 4);
    for (size_t i = 0; i < NUM_POINTS; ++i) {
      PulsarDopplerParams XLAL_INIT_DECL(point);
      XLAL_CHECK(XLALConvertPhysicalToSupersky(SC_RSSKY, rssky_point, &phys_points[i], rssky_transf, &ref_time) == XLAL_SUCCESS, XLAL_EFUNC);
      XLAL_CHECK(XLALConvertSuperskyToPhysical(&point, SC_RSSKY, rssky_point, rssky_transf, &ref_time) == XLAL_SUCCESS, XLAL_EFUNC);
      XLAL_CHECK(CompareDoppler(&phys_points[i], &point) == EXIT_SUCCESS, XLAL_EFUNC);
    }
    GFVEC(rssky_point);
  }

  // Check mismatches between pairs of points
  {
    gsl_vector *GAVEC(ussky_point_i, 5);
    gsl_vector *GAVEC(ussky_point_j, 5);
    gsl_vector *GAVEC(temp, 5);
    for (size_t i = 0; i < NUM_POINTS; ++i) {
      XLAL_CHECK(XLALConvertPhysicalToSupersky(SC_USSKY, ussky_point_i, &phys_points[i], NULL, &ref_time) == XLAL_SUCCESS, XLAL_EFUNC);
      for (size_t j = 0; j < NUM_POINTS; ++j) {
        XLAL_CHECK(XLALConvertPhysicalToSupersky(SC_USSKY, ussky_point_j, &phys_points[j], NULL, &ref_time) == XLAL_SUCCESS, XLAL_EFUNC);
        gsl_vector_sub(ussky_point_j, ussky_point_i);
        gsl_blas_dgemv(CblasNoTrans, 1.0, ussky_metric, ussky_point_j, 0.0, temp);
        double mismatch = 0.0;
        gsl_blas_ddot(ussky_point_j, temp, &mismatch);
        CHECK_RELERR(mismatch, phys_mismatch[i][j], 1e-2);
      }
    }
    GFVEC(ussky_point_i, ussky_point_j, temp);
  }
  {
    gsl_vector *GAVEC(rssky_point_i, 4);
    gsl_vector *GAVEC(rssky_point_j, 4);
    gsl_vector *GAVEC(temp, 4);
    for (size_t i = 0; i < NUM_POINTS; ++i) {
      XLAL_CHECK(XLALConvertPhysicalToSupersky(SC_RSSKY, rssky_point_i, &phys_points[i], rssky_transf, &ref_time) == XLAL_SUCCESS, XLAL_EFUNC);
      for (size_t j = 0; j < NUM_POINTS; ++j) {
        XLAL_CHECK(XLALConvertPhysicalToSupersky(SC_RSSKY, rssky_point_j, &phys_points[j], rssky_transf, &ref_time) == XLAL_SUCCESS, XLAL_EFUNC);
        gsl_vector_sub(rssky_point_j, rssky_point_i);
        gsl_blas_dgemv(CblasNoTrans, 1.0, rssky_metric, rssky_point_j, 0.0, temp);
        double mismatch = 0.0;
        gsl_blas_ddot(rssky_point_j, temp, &mismatch);
        CHECK_RELERR(mismatch, phys_mismatch[i][j], 1e-2);
      }
    }
    GFVEC(rssky_point_i, rssky_point_j, temp);
  }

  return XLAL_SUCCESS;

}

int main(void)
{

  // Compute supersky metrics
  LALSegList segments;
  {
    XLAL_CHECK_MAIN(XLALSegListInit(&segments) == XLAL_SUCCESS, XLAL_EFUNC);
    LALSeg segment;
    LIGOTimeGPS start_time = ref_time, end_time = ref_time;
    XLALGPSAdd(&start_time, -0.5 * Tspan);
    XLALGPSAdd(&end_time, 0.5 * Tspan);
    XLAL_CHECK_MAIN(XLALSegSet(&segment, &start_time, &end_time, 0) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN(XLALSegListAppend(&segments, &segment) == XLAL_SUCCESS, XLAL_EFUNC);
  }
  MultiLALDetector detectors = {
    .length = 1,
    .sites = { lalCachedDetectors[LAL_LLO_4K_DETECTOR] }
  };
  EphemerisData *edat =  XLALInitBarycenter(TEST_DATA_DIR "earth00-19-DE405.dat.gz",
                                            TEST_DATA_DIR "sun00-19-DE405.dat.gz");
  XLAL_CHECK_MAIN(edat != NULL, XLAL_EFUNC);
  gsl_matrix *rssky_metric = NULL, *rssky_transf = NULL, *ussky_metric = NULL;
  XLAL_CHECK_MAIN(XLALComputeSuperskyMetrics(&rssky_metric, &rssky_transf, &ussky_metric, 1, &ref_time, &segments, 100.0, &detectors, NULL, DETMOTION_SPIN | DETMOTION_PTOLEORBIT, edat) == XLAL_SUCCESS, XLAL_EFUNC);
  XLALSegListClear(&segments);
  XLALDestroyEphemerisData(edat);

  // Check supersky metrics
  XLAL_CHECK_MAIN(CheckSuperskyMetrics(
                    ussky_metric, ussky_metric_refs[0],
                    rssky_metric, rssky_metric_refs[0],
                    rssky_transf, rssky_transf_refs[0],
                    phys_mismatches[0]
                    ) == XLAL_SUCCESS, XLAL_EFUNC);

  // Cleanup
  GFMAT(rssky_metric, rssky_transf, ussky_metric);
  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
