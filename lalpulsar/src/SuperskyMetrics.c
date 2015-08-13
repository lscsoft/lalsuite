//
// Copyright (C) 2014, 2015 Karl Wette
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

#include <config.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_roots.h>

#include <lal/SuperskyMetrics.h>
#include <lal/LALConstants.h>
#include <lal/MetricUtils.h>
#include <lal/ExtrapolatePulsarSpins.h>

#include "GSLHelpers.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

// Square of a numbers
#define SQR(x)       ((x)*(x))

// Real part of square root of a number, i.e. zero if x < ~0
#define RE_SQRT(x)   sqrt(GSL_MAX(DBL_EPSILON, (x)))

// 2- and 3-dimensional vector dot products
#define DOT2(x,y)    ((x)[0]*(y)[0] + (x)[1]*(y)[1])
#define DOT3(x,y)    ((x)[0]*(y)[0] + (x)[1]*(y)[1] + (x)[2]*(y)[2])

///
/// Call XLALComputeDopplerPhaseMetric() to compute the phase metric for a given coordinate system.
///
static gsl_matrix *SM_ComputePhaseMetric(
  const DopplerCoordinateSystem *coords,	///< [in] Coordinate system to compute metric for
  const LIGOTimeGPS *ref_time,			///< [in] Reference time for the metric
  const LALSegList *segments,			///< [in] List of segments to average metric over
  const double fiducial_freq,			///< [in] Fiducial frequency for sky-position coordinates
  const MultiLALDetector *detectors,		///< [in] List of detectors to average metric over
  const MultiNoiseFloor *detector_weights,	///< [in] Weights used to combine single-detector metrics (default: unit weights)
  const DetectorMotionType detector_motion,	///< [in] Which detector motion to use
  const EphemerisData *ephemerides		///< [in] Earth/Sun ephemerides
  )
{

  // Check input
  XLAL_CHECK_NULL(coords != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(ref_time != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(segments != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(XLALSegListIsInitialized(segments), XLAL_EINVAL);
  XLAL_CHECK_NULL(segments->length > 0, XLAL_EINVAL);
  XLAL_CHECK_NULL(fiducial_freq > 0, XLAL_EINVAL);
  XLAL_CHECK_NULL(detectors != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(detectors->length > 0, XLAL_EINVAL);
  XLAL_CHECK_NULL(detector_motion > 0, XLAL_EINVAL);
  XLAL_CHECK_NULL(ephemerides != NULL, XLAL_EINVAL);

  // Create parameters struct for XLALComputeDopplerPhaseMetric()
  DopplerMetricParams XLAL_INIT_DECL(par);

  // Set coordinate system
  par.coordSys = *coords;

  // Set detector motion type
  par.detMotionType = detector_motion;

  // Set segment list
  par.segmentList = *segments;

  // Set detectors and detector weights
  par.multiIFO = *detectors;
  if (detector_weights != NULL) {
    par.multiNoiseFloor = *detector_weights;
  } else {
    par.multiNoiseFloor.length = 0;   // Indicates unit weights
  }

  // Set reference time and fiducial frequency
  par.signalParams.Doppler.refTime = *ref_time;
  par.signalParams.Doppler.fkdot[0] = fiducial_freq;

  // Do not project metric
  par.projectCoord = -1;

  // Do not include sky-position-dependent Roemer delay in time variable
  par.approxPhase = 1;

  // Call XLALComputeDopplerPhaseMetric() and check output
  DopplerPhaseMetric *metric = XLALComputeDopplerPhaseMetric(&par, ephemerides);
  XLAL_CHECK_NULL(metric != NULL && metric->g_ij != NULL, XLAL_EFUNC, "XLALComputeDopplerPhaseMetric() failed");

  // Extract metric
  gsl_matrix *g_ij = metric->g_ij;
  metric->g_ij = NULL;

  // Cleanup
  XLALDestroyDopplerPhaseMetric(metric);

  return g_ij;

}

///
/// Find a least-squares linear fit to the orbital X and Y metric elements using the frequency and
/// spindown metric elements. Outputs the intermediate \e fitted supersky metric, and updates the
/// reduced supersky metric coordinate transform data.
///
static int SM_ComputeFittedSuperskyMetric(
  gsl_matrix *fitted_ssky_metric,		///< [out] Fitted supersky metric
  gsl_matrix *rssky_transf,			///< [in,out] Reduced supersky metric coordinate transform data
  const gsl_matrix *ussky_metric,		///< [in] Unrestricted supersky metric
  const gsl_matrix *orbital_metric,		///< [in] Orbital metric in ecliptic coordinates
  const DopplerCoordinateSystem *ocoords,	///< [in] Coordinate system of orbital metric
  const size_t spindowns,			///< [in] Number of frequency+spindown coordinates
  const LIGOTimeGPS *ref_time,			///< [in] Reference time of the metrics
  const LALSegList *segments			///< [in] List of segments metric were averaged over
  )
{

  // Check input
  XLAL_CHECK(fitted_ssky_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(ussky_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(orbital_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(ref_time != NULL, XLAL_EFAULT);

  // Size of the frequency+spindowns block
  const size_t fsize = 1 + spindowns;

  // Allocate memory
  gsl_matrix *GAMAT(tmp, 2 + fsize, 2 + fsize);
  gsl_vector *GAVEC(tmpv, fsize);

  // Compute mid-time of segment list
  LIGOTimeGPS mid_time;
  {
    const LIGOTimeGPS *start_time = &(segments->segs[0].start);
    const LIGOTimeGPS *end_time   = &(segments->segs[segments->length - 1].end);
    const REAL8 time_span = XLALGPSDiff(end_time, start_time);
    mid_time = *start_time;
    XLALGPSAdd(&mid_time, 0.5 * time_span);
  }

  // Internal copy of orbital metric, and various transforms performed on it
  gsl_matrix *orb_metric = NULL, *mid_time_transf = NULL, *diag_norm_transf = NULL;

  // Transform reference time of orbital metric from reference time to segment list mid-time
  const REAL8 Dtau = XLALGPSDiff(&mid_time, ref_time);
  XLAL_CHECK(XLALChangeMetricReferenceTime(&orb_metric, &mid_time_transf, orbital_metric, ocoords, Dtau) == XLAL_SUCCESS, XLAL_EFUNC);

  // Diagonally-normalize orbital metric
  XLAL_CHECK(XLALDiagNormalizeMetric(&orb_metric, &diag_norm_transf, orb_metric) == XLAL_SUCCESS, XLAL_EFUNC);

  // 'fitA' contains the frequency and spindown elements of the orbital metric, used for fitting
  gsl_matrix *GAMAT(fitA, 2 + fsize, fsize);
  {
    gsl_matrix_view orb_metric_fspin = gsl_matrix_submatrix(orb_metric, 0, 2, 2 + fsize, fsize);
    gsl_matrix_memcpy(fitA, &orb_metric_fspin.matrix);
  }

  // Compute 'fitA^T * fitA'
  gsl_matrix *GAMAT(fitAt_fitA, fsize, fsize);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, fitA, fitA, 0.0, fitAt_fitA);

  // Find the singular value decomposition of 'fitA^T * fitA'
  gsl_matrix *GAMAT(svd_U, fsize, fsize);
  gsl_matrix *GAMAT(svd_V, fsize, fsize);
  gsl_vector *GAVEC(svd_S, fsize);
  gsl_matrix_memcpy(svd_U, fitAt_fitA);
  GCALL(gsl_linalg_SV_decomp(svd_U, svd_V, svd_S, tmpv));

  // The columns of 'fitc' contain the least-square fitting coefficients for the orbital X and Y metric elements:
  //    fitc(:,j) = inv(fitA^T * fitA) * fitA^T * orb_metric(:,j)
  // The singular decomposition of fitA^T * fitA is used for the inverse
  gsl_matrix *GAMAT(fitc, fsize, 2);
  for (size_t j = 0; j < 2; ++j) {
    gsl_vector_view orb_metric_j = gsl_matrix_column(orb_metric, j);
    gsl_vector_view fitc_j = gsl_matrix_column(fitc, j);
    gsl_blas_dgemv(CblasTrans, 1.0, fitA, &orb_metric_j.vector, 0.0, tmpv);
    GCALL(gsl_linalg_SV_solve(svd_U, svd_V, svd_S, tmpv, &fitc_j.vector));
  }

  // Construct the matrix 'subtract_orb', which subtracts the least-squares fit of
  // the orbital X and Y metric elements from the orbital metric. Its layout is:
  //
  //   #-------- sky --------#---f/s---#
  //   |                     |         |
  // sky  identity matrix I  |  zeros  |
  //   |                     |         |
  //   |---------------------#---------#
  //   |                     |         |
  // f/s        -fitc        |    I    |
  //   |                     |         |
  //   #---------------------#---------#
  //
  gsl_matrix *GAMAT(subtract_orb, 2 + fsize, 2 + fsize);
  {
    gsl_matrix_set_identity(subtract_orb);
    gsl_matrix_view subtract_orb_fspin_sky = gsl_matrix_submatrix(subtract_orb, 2, 0, fsize, 2);
    gsl_matrix_memcpy(&subtract_orb_fspin_sky.matrix, fitc);
    gsl_matrix_scale(&subtract_orb_fspin_sky.matrix, -1.0);
  }

  // Multiply 'subtract_orb' by the diagonal-normalization and reference time transforms,
  // to obtain the matrix that substracts the fit from the unconstrained supersky metric
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, diag_norm_transf, subtract_orb, 0.0, tmp);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, mid_time_transf, tmp, 0.0, subtract_orb);

  // Construct the matrix 'subtract_ussky', which subtracts the least-squares fit of the orbital X and Y metric
  // elements from the *unconstrained* supersky metric, which is in *equatorial* coordinates. Its layout is:
  //
  //   #------------------------------ sky -------------------------------#---f/s---#
  //   |                                                                  |         |
  // sky                        identity matrix I                         |  zeros  |
  //   |                                                                  |         |
  //   |------------------------------------------------------------------#---------#
  //   |            .                          .                          |         |
  // f/s sub_o(:,0) . sub_o(:,1)*LAL_COSIEARTH . sub_o(:,1)*LAL_SINIEARTH |    I    |
  //   |            .                          .                          |         |
  //   #------------------------------------------------------------------#---------#
  //
  // where 'sub_o' denotes 'subtract_orb'.
  gsl_matrix *GAMAT(subtract_ussky, 3 + fsize, 3 + fsize);
  {
    gsl_matrix_set_identity(subtract_ussky);
    gsl_matrix_view subtract_ussky_fspin_sky = gsl_matrix_submatrix(subtract_ussky, 3, 0, fsize, 3);
    gsl_matrix_view subtract_orb_fspin_sky = gsl_matrix_submatrix(subtract_orb, 2, 0, fsize, 2);
    {
      gsl_vector_view subtract_ussky_fspin_sky_col = gsl_matrix_column(&subtract_ussky_fspin_sky.matrix, 0);
      gsl_vector_view subtract_orb_fspin_sky_col = gsl_matrix_column(&subtract_orb_fspin_sky.matrix, 0);
      gsl_vector_memcpy(&subtract_ussky_fspin_sky_col.vector, &subtract_orb_fspin_sky_col.vector);
    }
    {
      gsl_vector_view subtract_ussky_fspin_sky_col = gsl_matrix_column(&subtract_ussky_fspin_sky.matrix, 1);
      gsl_vector_view subtract_orb_fspin_sky_col = gsl_matrix_column(&subtract_orb_fspin_sky.matrix, 1);
      gsl_vector_memcpy(&subtract_ussky_fspin_sky_col.vector, &subtract_orb_fspin_sky_col.vector);
      gsl_vector_scale(&subtract_ussky_fspin_sky_col.vector, LAL_COSIEARTH);
    }
    {
      gsl_vector_view subtract_ussky_fspin_sky_col = gsl_matrix_column(&subtract_ussky_fspin_sky.matrix, 2);
      gsl_vector_view subtract_orb_fspin_sky_col = gsl_matrix_column(&subtract_orb_fspin_sky.matrix, 1);
      gsl_vector_memcpy(&subtract_ussky_fspin_sky_col.vector, &subtract_orb_fspin_sky_col.vector);
      gsl_vector_scale(&subtract_ussky_fspin_sky_col.vector, LAL_SINIEARTH);
    }
  }

  // Transform the unconstrained supersky metric to the intermediate fitted supersky metric
  XLAL_CHECK(XLALTransformMetric(&fitted_ssky_metric, subtract_ussky, ussky_metric) == XLAL_SUCCESS, XLAL_EFUNC);

  // Extract the sky offset vectors from 'subtract_ussky', and subtract them from the reduced supersky coordinate transform data
  {
    gsl_matrix_view subtract_ussky_fspin_sky = gsl_matrix_submatrix(subtract_ussky, 3, 0, fsize, 3);
    gsl_matrix_view sky_offsets = gsl_matrix_submatrix(rssky_transf, 3, 0, fsize, 3);
    gsl_matrix_sub(&sky_offsets.matrix, &subtract_ussky_fspin_sky.matrix);
  }

  // Cleanup
  GFMAT(diag_norm_transf, fitA, fitAt_fitA, fitc, mid_time_transf, orb_metric, subtract_orb, subtract_ussky, svd_U, svd_V, tmp);
  GFVEC(svd_S, tmpv);

  return XLAL_SUCCESS;

}

///
/// Decouple the sky--sky and freq+spin--freq+spin blocks of the fitted supersky metric. Outputs the
/// intermediate \e decoupled supersky metric, and updates the reduced supersky metric coordinate
/// transform data.
///
static int SM_ComputeDecoupledSuperskyMetric(
  gsl_matrix *decoupled_ssky_metric,		///< [out] Decoupled supersky metric
  gsl_matrix *rssky_transf,			///< [in,out] Reduced supersky metric coordinate transform data
  const gsl_matrix *fitted_ssky_metric,		///< [in] Fitted supersky metric
  const size_t spindowns			///< [in] Number of frequency+spindown coordinates
  )
{

  // Check input
  XLAL_CHECK(decoupled_ssky_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(fitted_ssky_metric != NULL, XLAL_EFAULT);

  // Size of the frequency+spindowns block
  const size_t fsize = 1 + spindowns;

  // Copy fitted metric to decoupled metric
  gsl_matrix_memcpy(decoupled_ssky_metric, fitted_ssky_metric);

  // Create views of the sky--sky, freq+spin--freq+spin, and off-diagonal blocks
  gsl_matrix_view sky_sky     = gsl_matrix_submatrix(decoupled_ssky_metric, 0, 0, 3, 3);
  gsl_matrix_view sky_fspin   = gsl_matrix_submatrix(decoupled_ssky_metric, 0, 3, 3, fsize);
  gsl_matrix_view fspin_sky   = gsl_matrix_submatrix(decoupled_ssky_metric, 3, 0, fsize, 3);
  gsl_matrix_view fspin_fspin = gsl_matrix_submatrix(decoupled_ssky_metric, 3, 3, fsize, fsize);

  // Diagonal-normalise the freq+spin--freq+spin block
  gsl_matrix *fspin_fspin_dnorm = NULL, *fspin_fspin_dnorm_transf = NULL;
  XLAL_CHECK(XLALDiagNormalizeMetric(&fspin_fspin_dnorm, &fspin_fspin_dnorm_transf, &fspin_fspin.matrix) == XLAL_SUCCESS, XLAL_EFUNC);

  // Invert the freq+spin--freq+spin block
  gsl_matrix *GAMAT(fspin_fspin_dnorm_LU, fsize, fsize);
  gsl_matrix *GAMAT(fspin_fspin_dnorm_inv, fsize, fsize);
  gsl_permutation *GAPERM(fspin_fspin_dnorm_LU_perm, fsize);
  int fspin_fspin_dnorm_LU_sign = 0;
  gsl_matrix_memcpy(fspin_fspin_dnorm_LU, fspin_fspin_dnorm);
  GCALL(gsl_linalg_LU_decomp(fspin_fspin_dnorm_LU, fspin_fspin_dnorm_LU_perm, &fspin_fspin_dnorm_LU_sign));
  GCALL(gsl_linalg_LU_invert(fspin_fspin_dnorm_LU, fspin_fspin_dnorm_LU_perm, fspin_fspin_dnorm_inv));

  // Compute the additional sky offsets required to decouple the sky--sky and frequency blocks:
  //   decouple_sky_offsets = fspin_fspin_dnorm_transf * inv(fspin_fspin_dnorm) * fspin_fspin_dnorm_transf * fspin_sky
  // Uses fspin_sky as a temporary matrix, since it will be zeroed out anyway
  gsl_matrix *GAMAT(decouple_sky_offsets, fsize, 3);
  gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, fspin_fspin_dnorm_transf, &fspin_sky.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, fspin_fspin_dnorm_inv, &fspin_sky.matrix, 0.0, decouple_sky_offsets);
  gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, fspin_fspin_dnorm_transf, decouple_sky_offsets);

  // Add the additional sky offsets to the reduced supersky coordinate transform data
  gsl_matrix_view sky_offsets = gsl_matrix_submatrix(rssky_transf, 3, 0, fsize, 3);
  gsl_matrix_add(&sky_offsets.matrix, decouple_sky_offsets);

  // Apply the decoupling transform to the sky--sky block:
  //   sky_sky = sky_sky - sky_fspin * decouplp_sky_offsets
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &sky_fspin.matrix, decouple_sky_offsets, 1.0, &sky_sky.matrix);

  // Zero out the off-diagonal blocks
  gsl_matrix_set_zero(&sky_fspin.matrix);
  gsl_matrix_set_zero(&fspin_sky.matrix);

  // Cleanup
  GFPERM(fspin_fspin_dnorm_LU_perm);
  GFMAT(fspin_fspin_dnorm, fspin_fspin_dnorm_LU, fspin_fspin_dnorm_inv, fspin_fspin_dnorm_transf, decouple_sky_offsets);

  return XLAL_SUCCESS;

}

///
/// Align the sky--sky block of the decoupled supersky metric with its eigenvalues by means of a
/// rotation. Outputs the intermediate \e aligned supersky metric, and updates the reduced supersky
/// metric coordinate transform data.
///
static int SM_ComputeAlignedSuperskyMetric(
  gsl_matrix *aligned_ssky_metric,		///< [out] Aligned supersky metric
  gsl_matrix *rssky_transf,			///< [in,out] Reduced supersky metric coordinate transform data
  const gsl_matrix *decoupled_ssky_metric,	///< [in] Decoupled supersky metric
  const size_t spindowns			///< [in] Number of frequency+spindown coordinates
  )
{

  // Check input
  XLAL_CHECK(aligned_ssky_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(decoupled_ssky_metric != NULL, XLAL_EFAULT);

  // Size of the frequency+spindowns block
  const size_t fsize = 1 + spindowns;

  // Allocate memory
  gsl_matrix *GAMAT(tmp, fsize, 3);

  // Copy decoupled metric to aligned metric
  gsl_matrix_memcpy(aligned_ssky_metric, decoupled_ssky_metric);

  // Compute the eigenvalues/vectors of the sky--sky block
  gsl_vector *GAVEC(sky_eval, 3);
  gsl_matrix *GAMAT(sky_evec, 3, 3);
  gsl_eigen_symmv_workspace *GALLOC(wksp, gsl_eigen_symmv_alloc(3));
  gsl_matrix_view sky_sky = gsl_matrix_submatrix(aligned_ssky_metric, 0, 0, 3, 3);
  GCALL(gsl_eigen_symmv(&sky_sky.matrix, sky_eval, sky_evec, wksp));

  // Sort the eigenvalues/vectors by descending absolute eigenvalue
  GCALL(gsl_eigen_symmv_sort(sky_eval, sky_evec, GSL_EIGEN_SORT_ABS_DESC));

  // Set the sky--sky block to the diagonal matrix of eigenvalues
  gsl_matrix_set_zero(&sky_sky.matrix);
  gsl_vector_view sky_sky_diag = gsl_matrix_diagonal(&sky_sky.matrix);
  gsl_vector_memcpy(&sky_sky_diag.vector, sky_eval);

  // Ensure that the matrix of eigenvalues has a positive diagonal; this and
  // the determinant constraints ensures fully constraints the eigenvector signs
  for (size_t j = 0; j < 3; ++j) {
    gsl_vector_view col = gsl_matrix_column(sky_evec, j);
    if (gsl_vector_get(&col.vector, j) < 0.0) {
      gsl_vector_scale(&col.vector, -1.0);
    }
  }

  // Store the alignment transform in the reduced supersky coordinate transform data
  gsl_matrix_view align_sky = gsl_matrix_submatrix(rssky_transf, 0, 0, 3, 3);
  gsl_matrix_transpose_memcpy(&align_sky.matrix, sky_evec);

  // Ensure that the alignment transform has a positive determinant,
  // to ensure that that it represents a rotation
  gsl_permutation *GAPERM(LU_perm, 3);
  int LU_sign = 0;
  GCALL(gsl_linalg_LU_decomp(sky_evec, LU_perm, &LU_sign));
  if (gsl_linalg_LU_det(sky_evec, LU_sign) < 0.0) {
    gsl_vector_view col = gsl_matrix_column(&align_sky.matrix, 2);
    gsl_vector_scale(&col.vector, -1.0);
  }

  // Multiply the sky offsets by the alignment transform to transform to aligned sky coordinates:
  //   aligned_sky_off = sky_offsets * alignsky^T;
  gsl_matrix_view aligned_sky_offsets = gsl_matrix_submatrix(rssky_transf, 3, 0, fsize, 3);
  gsl_matrix_memcpy(tmp, &aligned_sky_offsets.matrix);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp, &align_sky.matrix, 0.0, &aligned_sky_offsets.matrix);

  // Cleanup
  gsl_eigen_symmv_free(wksp);
  GFMAT(sky_evec, tmp);
  GFPERM(LU_perm);
  GFVEC(sky_eval);

  return XLAL_SUCCESS;

}

///
/// Extract the reduced supersky metric from the aligned supersky metric.
///
static int SM_ExtractReducedSuperskyMetric(
  gsl_matrix *rssky_metric,			///< [out] Reduced supersky metric
  gsl_matrix *rssky_transf,			///< [in,out] Reduced supersky metric coordinate transform data
  const gsl_matrix *aligned_ssky_metric		///< [in] Aligned supersky metric
  )
{

  // Check input
  XLAL_CHECK(rssky_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(aligned_ssky_metric != NULL, XLAL_EFAULT);
  const size_t n = aligned_ssky_metric->size1;
  const size_t m = rssky_metric->size1;

  // Internal copy of aligned supersky metric
  gsl_matrix *GAMAT(aln_metric, n, n);
  gsl_matrix_memcpy(aln_metric, aligned_ssky_metric);

  // Move the 3rd row/column of 'aln_metric', which is the 'n_c' sky coordinate with
  // the smallest eigenvalue, to the last row/column, so it can be easily dropped
  for (size_t i = 2; i + 1 < n; ++i) {
    gsl_matrix_swap_rows(aln_metric, i, i + 1);
    gsl_matrix_swap_columns(aln_metric, i, i + 1);
  }

  // Move the 3rd row/column of 'aln_metric', which is *now* the frequency,
  // to the second-to-last row/column, i.e. still before 'n_c'
  for (size_t i = 2; i + 2 < n; ++i) {
    gsl_matrix_swap_rows(aln_metric, i, i + 1);
    gsl_matrix_swap_columns(aln_metric, i, i + 1);
  }

  // Copy the first 'm' dimensions of 'aln_metric' to 'rssky_metric', dropping 'n_c'
  {
    gsl_matrix_view aln_metric_nm1_nm1 = gsl_matrix_submatrix(aln_metric, 0, 0, m, m);
    gsl_matrix_memcpy(rssky_metric, &aln_metric_nm1_nm1.matrix);
  }

  // Move the 4th row of 'rssky_transf', which is the coordinate
  // transform data for frequency, to the last row
  for (size_t i = 3; i + 1 < n; ++i) {
    gsl_matrix_swap_rows(rssky_transf, i, i + 1);
  }

  // Ensure reduced supersky metric is symmetric
  for (size_t i = 0; i < m; ++i) {
    for (size_t j = i + 1; j < m; ++j) {
      const double gij = gsl_matrix_get(rssky_metric, i, j);
      const double gji = gsl_matrix_get(rssky_metric, j, i);
      const double g = 0.5 * (gij + gji);
      gsl_matrix_set(rssky_metric, i, j, g);
      gsl_matrix_set(rssky_metric, j, i, g);
    }
  }

  // Ensure reduced supersky metric is positive definite
  for (size_t s = 1; s <= rssky_metric->size1; ++s) {
    gsl_matrix_view rssky_metric_s = gsl_matrix_submatrix(rssky_metric, 0, 0, s, s);
    const double det_s = XLALMetricDeterminant(&rssky_metric_s.matrix);
    XLAL_CHECK(det_s > 0, XLAL_EFAILED, "Reduced supersky metric is not positive definite (s=%zu, det_s=%0.3e)", s, det_s);
  }

  return XLAL_SUCCESS;

}

int XLALComputeSuperskyMetrics(
  gsl_matrix **p_rssky_metric,
  gsl_matrix **p_rssky_transf,
  gsl_matrix **p_ussky_metric,
  const size_t spindowns,
  const LIGOTimeGPS *ref_time,
  const LALSegList *segments,
  const double fiducial_freq,
  const MultiLALDetector *detectors,
  const MultiNoiseFloor *detector_weights,
  const DetectorMotionType detector_motion,
  const EphemerisData *ephemerides
  )
{

  // Check input
  XLAL_CHECK(p_rssky_metric == NULL || *p_rssky_metric == NULL, XLAL_EINVAL, "'*p_rssky_metric' must be NULL of 'p_rssky_metric' is non-NULL");
  XLAL_CHECK(p_rssky_transf == NULL || *p_rssky_transf == NULL, XLAL_EINVAL, "'*p_rssky_transf' must be NULL of 'p_rssky_transf' is non-NULL");
  XLAL_CHECK(p_ussky_metric == NULL || *p_ussky_metric == NULL, XLAL_EINVAL, "'*p_ussky_metric' must be NULL of 'p_ussky_metric' is non-NULL");
  XLAL_CHECK((p_rssky_metric != NULL) == (p_rssky_transf != NULL), XLAL_EINVAL, "Both 'p_rssky_metric' and 'p_rssky_transf' must be either NULL or non-NULL");
  XLAL_CHECK((p_rssky_metric != NULL) || (p_ussky_metric != NULL), XLAL_EINVAL, "At least one of 'p_rssky_metric' or 'p_ussky_metric' must be non-NULL");
  XLAL_CHECK(spindowns <= 3, XLAL_EINVAL);
  XLAL_CHECK(ref_time != NULL, XLAL_EFAULT);
  XLAL_CHECK(segments != NULL, XLAL_EFAULT);
  XLAL_CHECK(XLALSegListIsInitialized(segments), XLAL_EINVAL);
  XLAL_CHECK(segments->length > 0, XLAL_EINVAL);
  XLAL_CHECK(fiducial_freq > 0, XLAL_EINVAL);
  XLAL_CHECK(detectors != NULL, XLAL_EFAULT);
  XLAL_CHECK(detectors->length > 0, XLAL_EINVAL);
  XLAL_CHECK(detector_motion > 0, XLAL_EINVAL);
  XLAL_CHECK(ephemerides != NULL, XLAL_EINVAL);

  // Size of the frequency+spindowns block
  const size_t fsize = 1 + spindowns;

  // Fiducial frequency at which to numerically calculate metrics, which
  // are then rescaled to input 'fiducial_freq' based on known scalings
  const double fiducial_calc_freq = 100.0;

  // Build coordinate system for the unrestricted supersky metric and orbital metric
  DopplerCoordinateSystem XLAL_INIT_DECL(ucoords);
  DopplerCoordinateSystem XLAL_INIT_DECL(ocoords);
  {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_N3X_EQU;
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_N3Y_EQU;
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_N3Z_EQU;
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_FREQ;
  }
  {
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_N3OX_ECL;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_N3OY_ECL;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_FREQ;
  }
  if (spindowns >= 1) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F1DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F1DOT;
  }
  if (spindowns >= 2) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F2DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F2DOT;
  }
  if (spindowns >= 3) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F3DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F3DOT;
  }

  // Compute the unrestricted supersky metric
  gsl_matrix *ussky_metric = SM_ComputePhaseMetric(&ucoords, ref_time, segments, fiducial_calc_freq, detectors, detector_weights, detector_motion, ephemerides);
  XLAL_CHECK(ussky_metric != NULL, XLAL_EFUNC);

  // Compute the reduced supersky metric and coordinate transform data
  if (p_rssky_metric != NULL) {

    // Allocate memory
    GAMAT(*p_rssky_metric, 2 + fsize, 2 + fsize);
    GAMAT(*p_rssky_transf, 3 + fsize, 3);
    gsl_matrix *GAMAT(interm_ssky_metric, 3 + fsize, 3 + fsize);

    // Compute the orbital metric in ecliptic coordinates
    gsl_matrix *orbital_metric = SM_ComputePhaseMetric(&ocoords, ref_time, segments, fiducial_calc_freq, detectors, detector_weights, detector_motion, ephemerides);
    XLAL_CHECK(orbital_metric != NULL, XLAL_EFUNC);

    // Compute the reduced supersky metric from the unrestricted supersky metric and the orbital metric
    XLAL_CHECK(SM_ComputeFittedSuperskyMetric(interm_ssky_metric, *p_rssky_transf, ussky_metric, orbital_metric, &ocoords, spindowns, ref_time, segments) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK(SM_ComputeDecoupledSuperskyMetric(interm_ssky_metric, *p_rssky_transf, interm_ssky_metric, spindowns) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK(SM_ComputeAlignedSuperskyMetric(interm_ssky_metric, *p_rssky_transf, interm_ssky_metric, spindowns) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK(SM_ExtractReducedSuperskyMetric(*p_rssky_metric, *p_rssky_transf, interm_ssky_metric) == XLAL_SUCCESS, XLAL_EFUNC);

    // Cleanup
    GFMAT(orbital_metric, interm_ssky_metric);

  }

  // Return or free unrestricted supersky metric
  if (p_ussky_metric != NULL) {
    *p_ussky_metric = ussky_metric;
  } else {
    GFMAT(ussky_metric);
  }

  // Rescale metrics to input 'fiducial_freq' based on known scalings
  const double fiducial_scale = fiducial_freq / fiducial_calc_freq;
  if (p_ussky_metric != NULL) {
    gsl_matrix_view sky_sky = gsl_matrix_submatrix(*p_ussky_metric, 0, 0, 3, 3);
    gsl_matrix_scale(&sky_sky.matrix, SQR(fiducial_scale));
    gsl_matrix_view sky_freq = gsl_matrix_submatrix(*p_ussky_metric, 0, 3, 3, fsize);
    gsl_matrix_scale(&sky_freq.matrix, fiducial_scale);
    gsl_matrix_view freq_sky = gsl_matrix_submatrix(*p_ussky_metric, 3, 0, fsize, 3);
    gsl_matrix_scale(&freq_sky.matrix, fiducial_scale);
  }
  if (p_rssky_metric != NULL) {
    gsl_matrix_view sky_sky = gsl_matrix_submatrix(*p_rssky_metric, 0, 0, 2, 2);
    gsl_matrix_scale(&sky_sky.matrix, SQR(fiducial_scale));
    gsl_matrix_view sky_offsets = gsl_matrix_submatrix(*p_rssky_transf, 3, 0, fsize, 3);
    gsl_matrix_scale(&sky_offsets.matrix, fiducial_scale);
  }

  return XLAL_SUCCESS;

}

/**
 * Convert from 2-dimensional reduced supersky coordinates to 3-dimensional aligned sky
 * coordinates.
 *
 * The 2-dimensional reduced supersky coordinates \c rss = (\c A, \c B) encode the two
 * hemispheres of the sky as two neighbouring unit disks. The conversion to 3-dimensional
 * aligned sky coordinates is illustrated in the following diagram:
 *
 \verbatim
 | as[1] =    __________________________________________
 |     B = 1_|         _____             _____         |
 |           |      .-'     '-.       .-'     '-.      |
 |           |    .'           '.   .'           '.    |
 |           |   /               \ /               \   |
 |           |  ;                 ;                 ;  |
 |         0-|  |                 |                 |  |
 |           |  ;                 ;                 ;  |
 |           |   \               / \               /   |
 |           |    '.           .'   '.           .'    |
 |          _|      '-._____.-'       '-._____.-'      |
 |        -1 |__.________.________.________.________.__|
 |              '        '        '        '        '
 |         A = -2       -1        0        1        2
 |     as[0] =  1        0       -1        0        1
 |     as[2] =  0       -1        0        1        0
 \endverbatim
 *
 * Points outside the unit disks are moved radially onto their boundaries.
 */
static void SM_ReducedToAligned(
  double as[3],					///< [out] 3-dimensional aligned sky coordinates
  const gsl_vector *rss				///< [in] 2-dimensional reduced supersky coordinates
  )
{
  const double A = gsl_vector_get(rss, 0);
  const double B = gsl_vector_get(rss, 1);
  const double dA = fabs(A) - 1.0;
  const double R = sqrt(SQR(dA) + SQR(B));
  const double Rmax = GSL_MAX(1.0, R);
  as[0] = dA / Rmax;
  as[1] = B / Rmax;
  as[2] = GSL_SIGN(A) * RE_SQRT(1.0 - DOT2(as, as));
}

///
/// Convert from 3-dimensional aligned sky coordinates to 2-dimensional reduced supersky
/// coordinates.
///
static void SM_AlignedToReduced(
  gsl_vector *rss,				///< [out] 2-dimensional reduced supersky coordinates
  const double as[3]				///< [in] 3-dimensional aligned sky coordinates
  )
{
  const double r = sqrt(DOT3(as, as));
  const double A = GSL_SIGN(as[2]) * ((as[0] / r) + 1.0);
  const double B = as[1] / r;
  gsl_vector_set(rss, 0, A);
  gsl_vector_set(rss, 1, B);
}

int XLALConvertSuperskyCoordinates(
  const SuperskyCoordinates out,
  gsl_matrix **out_points,
  const SuperskyCoordinates in,
  const gsl_matrix *in_points,
  const gsl_matrix *rssky_transf
  )
{

  // Check input
  XLAL_CHECK(out < SC_MAX, XLAL_EINVAL);
  XLAL_CHECK(in < SC_MAX, XLAL_EINVAL);
  XLAL_CHECK(out_points != NULL, XLAL_EFAULT);
  XLAL_CHECK(in_points != NULL, XLAL_EINVAL);
  XLAL_CHECK(rssky_transf != NULL || (out != SC_RSSKY && in != SC_RSSKY), XLAL_EINVAL);
  XLAL_CHECK(rssky_transf == NULL || rssky_transf->size1 > 3, XLAL_ESIZE);
  XLAL_CHECK(rssky_transf == NULL || rssky_transf->size2 == 3, XLAL_ESIZE);

  // Deduce number of input sky coordinates, and frequency/spindown coordinates
  const size_t in_ssize = (in == SC_USSKY) ? 3 : 2;
  XLAL_CHECK(in_points->size1 > in_ssize, XLAL_EINVAL);
  const size_t fsize = in_points->size1 - in_ssize;

  // Resize or allocate output points matrix, if required
  const size_t out_ssize = (out == SC_USSKY) ? 3 : 2;
  const size_t out_rows = fsize + out_ssize;
  if (*out_points != NULL) {
    if ((*out_points)->size1 != out_rows || (*out_points)->size2 != in_points->size2) {
      GFMAT(*out_points);
      *out_points = NULL;
    }
  }
  if (*out_points == NULL) {
    GAMAT(*out_points, out_rows, in_points->size2);
  }

  // If input and output coordinate systems are the same, copy input matrix and exit
  if (in == out) {
    gsl_matrix_memcpy(*out_points, in_points);
    return XLAL_SUCCESS;
  }

  // Iterate over input points
  for (size_t j = 0; j < in_points->size2; ++j) {

    // Create array for point in intermediate coordinates
    double tmp[3 + fsize];
    gsl_vector_view tmp_sky = gsl_vector_view_array(&tmp[0], 3);
    gsl_vector_view tmp_fspin = gsl_vector_view_array(&tmp[3], fsize);

    // Copy input point to intermediate point
    for (size_t i = 0; i < in_ssize; ++i) {
      tmp[i] = gsl_matrix_get(in_points, i, j);
    }
    for (size_t i = 0; i < fsize; ++i) {
      tmp[3 + i] = gsl_matrix_get(in_points, in_ssize + i, j);
    }

    // Initialise current coordinate system
    SuperskyCoordinates curr = in;

    // Convert physical coordinates to supersky coordinates
    if (curr == SC_PHYS && out > curr) {

      // Convert right ascension and declination to supersky position
      const double alpha = tmp[0];
      const double delta = tmp[1];
      const double cos_delta = cos(delta);
      tmp[0] = cos(alpha) * cos_delta;
      tmp[1] = sin(alpha) * cos_delta;
      tmp[2] = sin(delta);

      // Update current coordinate system
      curr = SC_USSKY;

    }

    // Convert supersky coordinates to reduced supersky coordinates
    if (curr == SC_USSKY && out > curr) {

      // Move frequency to after spindowns
      const double freq = tmp[3];
      memmove(&tmp[3], &tmp[4], (fsize - 1) * sizeof(tmp[0]));
      tmp[2 + fsize] = freq;

      // Create views of the sky alignment transform and sky offset vectors
      gsl_matrix_const_view align_sky = gsl_matrix_const_submatrix(rssky_transf, 0, 0, 3, 3);
      gsl_matrix_const_view sky_offsets = gsl_matrix_const_submatrix(rssky_transf, 3, 0, fsize, 3);

      // Apply the alignment transform to the supersky position to produced the aligned sky position:
      //   asky = align_sky * ssky
      double asky[3];
      gsl_vector_view asky_v = gsl_vector_view_array(asky, 3);
      gsl_blas_dgemv(CblasNoTrans, 1.0, &align_sky.matrix, &tmp_sky.vector, 0.0, &asky_v.vector);

      // Add the inner product of the sky offsets with the aligned sky position
      // to the supersky spins and frequency to get the reduced supersky quantities:
      //   rssky_fspin[i] = ussky_fspin[i] + dot(sky_offsets[i], asky)
      gsl_blas_dgemv(CblasNoTrans, 1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &tmp_fspin.vector);

      // Convert from 3-dimensional aligned sky coordinates to 2-dimensional reduced supersky coordinates
      SM_AlignedToReduced(&tmp_sky.vector, asky);

      // Update current coordinate system
      curr = SC_RSSKY;

    }

    // Convert reduced supersky coordinates to supersky coordinates
    if (curr == SC_RSSKY && out < curr) {

      // Create views of the sky alignment transform and sky offset vectors
      gsl_matrix_const_view align_sky = gsl_matrix_const_submatrix(rssky_transf, 0, 0, 3, 3);
      gsl_matrix_const_view sky_offsets = gsl_matrix_const_submatrix(rssky_transf, 3, 0, fsize, 3);

      // Convert from 2-dimensional reduced supersky coordinates to 3-dimensional aligned sky coordinates
      double asky[3];
      SM_ReducedToAligned(asky, &tmp_sky.vector);
      gsl_vector_view asky_v = gsl_vector_view_array(asky, 3);

      // Subtract the inner product of the sky offsets with the aligned sky position
      // from the reduced supersky spins and frequency to get the supersky quantities:
      //   ussky_fspin[i] = rssky_fspin[i] - dot(sky_offsets[i], asky)
      gsl_blas_dgemv(CblasNoTrans, -1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &tmp_fspin.vector);

      // Apply the inverse alignment transform to the aligned sky position to produced the supersky position:
      //   ssky = align_sky^T * asky
      gsl_blas_dgemv(CblasTrans, 1.0, &align_sky.matrix, &asky_v.vector, 0.0, &tmp_sky.vector);

      // Move frequency to before spindowns
      const double freq = tmp[2 + fsize];
      memmove(&tmp[4], &tmp[3], (fsize - 1) * sizeof(tmp[0]));
      tmp[3] = freq;

      // Update current coordinate system
      curr = SC_USSKY;

    }

    // Convert supersky coordinates to physical coordinates
    if (curr == SC_USSKY && out < curr) {

      // Convert supersky position to right ascension and declination
      const double nx = tmp[0];
      const double ny = tmp[1];
      const double nz = tmp[2];
      tmp[0] = atan2(ny, nx);
      tmp[1] = atan2(nz, sqrt(SQR(nx) + SQR(ny)));
      XLALNormalizeSkyPosition(&tmp[0], &tmp[1]);

      // Update current coordinate system
      curr = SC_PHYS;

    }

    // Check that correct coordinate system has been converted to
    XLAL_CHECK(curr == out, XLAL_EFAILED);

    // Copy intermediate point to output point
    for (size_t i = 0; i < out_ssize; ++i) {
      gsl_matrix_set(*out_points, i, j, tmp[i]);
    }
    for (size_t i = 0; i < fsize; ++i) {
      gsl_matrix_set(*out_points, out_ssize + i, j, tmp[3 + i]);
    }

  }

  return XLAL_SUCCESS;

}

int XLALConvertPhysicalToSupersky(
  const SuperskyCoordinates out,
  gsl_vector *out_point,
  const PulsarDopplerParams *in_phys,
  const gsl_matrix *rssky_transf,
  const LIGOTimeGPS *ref_time
  )
{

  // Check input
  XLAL_CHECK(SC_PHYS < out && out < SC_MAX, XLAL_EINVAL);
  XLAL_CHECK(out_point != NULL, XLAL_EFAULT);
  XLAL_CHECK(in_phys != NULL, XLAL_EFAULT);
  XLAL_CHECK(ref_time != NULL, XLAL_EFAULT);

  // Deduce number of sky coordinates, and frequency/spindown coordinates
  const size_t ssize = (out == SC_USSKY) ? 3 : 2;
  XLAL_CHECK(out_point->size > ssize, XLAL_EINVAL);
  const size_t fsize = out_point->size - ssize;
  XLAL_CHECK(fsize <= PULSAR_MAX_SPINS, XLAL_EFAILED);

  // Transform input physical point to reference time of coordinate transform data
  PulsarDopplerParams in_phys_ref = *in_phys;
  {
    const REAL8 dtau = XLALGPSDiff(ref_time, &in_phys_ref.refTime);
    XLAL_CHECK(XLALExtrapolatePulsarSpins(in_phys_ref.fkdot, in_phys_ref.fkdot, dtau) == XLAL_SUCCESS, XLAL_EFUNC);
  }

  // Copy input physical point to array
  double in_point[2 + fsize];
  in_point[0] = in_phys_ref.Alpha;
  in_point[1] = in_phys_ref.Delta;
  memcpy(&in_point[2], in_phys_ref.fkdot, fsize * sizeof(in_point[0]));

  // Convert input physical point to output supersky coordinate point
  gsl_matrix_view out_point_view = gsl_matrix_view_vector(out_point, out_point->size, 1);
  gsl_matrix_const_view in_point_view = gsl_matrix_const_view_array(in_point, 2 + fsize, 1);
  gsl_matrix *out_point_view_ptr = &out_point_view.matrix;
  XLAL_CHECK(XLALConvertSuperskyCoordinates(out, &out_point_view_ptr, SC_PHYS, &in_point_view.matrix, rssky_transf) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK(out_point_view_ptr == &out_point_view.matrix, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

int XLALConvertSuperskyToPhysical(
  PulsarDopplerParams *out_phys,
  const SuperskyCoordinates in,
  const gsl_vector *in_point,
  const gsl_matrix *rssky_transf,
  const LIGOTimeGPS *ref_time
  )
{

  // Check input
  XLAL_CHECK(out_phys != NULL, XLAL_EFAULT);
  XLAL_CHECK(SC_PHYS < in && in < SC_MAX, XLAL_EINVAL);
  XLAL_CHECK(in_point != NULL, XLAL_EFAULT);
  XLAL_CHECK(ref_time != NULL, XLAL_EFAULT);

  // Deduce number of sky coordinates, and frequency/spindown coordinates
  const size_t ssize = (in == SC_USSKY) ? 3 : 2;
  XLAL_CHECK(in_point->size > ssize, XLAL_EINVAL);
  const size_t fsize = in_point->size - ssize;
  XLAL_CHECK(fsize <= PULSAR_MAX_SPINS, XLAL_EFAILED);

  // Convert input supersky coordinate point to output physical point
  double out_point[2 + fsize];
  gsl_matrix_view out_point_view = gsl_matrix_view_array(out_point, 2 + fsize, 1);
  gsl_matrix_const_view in_point_view = gsl_matrix_const_view_vector(in_point, in_point->size, 1);
  gsl_matrix *out_point_view_ptr = &out_point_view.matrix;
  XLAL_CHECK(XLALConvertSuperskyCoordinates(SC_PHYS, &out_point_view_ptr, in, &in_point_view.matrix, rssky_transf) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK(out_point_view_ptr == &out_point_view.matrix, XLAL_EFAILED);

  // Copy output physical point from array
  out_phys->Alpha = out_point[0];
  out_phys->Delta = out_point[1];
  memcpy(out_phys->fkdot, &out_point[2], fsize * sizeof(out_point[0]));

  // Set output physical point reference time to that of of coordinate transform data
  out_phys->refTime = *ref_time;

  return XLAL_SUCCESS;

}

// Structure which stores bounds data for PhysicalSkyBound()
#define MAX_PHYS_SKY_BOUND_DATA 6
typedef struct {
  double max_A;
  int type;
  double r;
  double na0;
  bool altroot;
  double angle0;
  double C;
  double S;
  double Z;
} PhysicalSkyBoundData;

static double PhysicalSkyBound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_vector *point
  )
{

  // Get bounds data
  const PhysicalSkyBoundData *bdata = (const PhysicalSkyBoundData *) data;

  // Decode the reduced supersky coordinates to get
  //   na = as[0] = Q_na . n
  const double A = gsl_vector_get(point, 0);
  const double rssky[2] = { A, 0 };
  gsl_vector_const_view rssky_view = gsl_vector_const_view_array(rssky, 2);
  double as[3];
  SM_ReducedToAligned(as, &rssky_view.vector);
  const double na = as[0];

  // Absolute limiting bound on 'nb = +/- sqrt(1 - na^2)'
  const double limit = RE_SQRT(1 - SQR(na));

  // Loop over bounds to find bound which currently applies, based on 'max_A'
  double bound = GSL_NAN;
  for (size_t i = 0; i < MAX_PHYS_SKY_BOUND_DATA; ++i) {
    const PhysicalSkyBoundData b = bdata[i];
    if (A <= b.max_A) {

      if (b.type != 0) {

        // Set bound 'nb = +/- sqrt(1 - na^2)'
        bound = b.type * limit;

      } else {

        // Set bound 'nb' to either a constant right ascension or constant declination bound,
        // depending on bounds data set by XLALSetSuperskyLatticeTilingPhysicalSkyBounds()
        const double c = (na - b.na0) / b.r;
        double angle = asin(GSL_MAX(-1, GSL_MIN(c, 1)));
        if (b.altroot) {
          angle = LAL_PI - angle;
        }
        angle -= b.angle0;
        bound = b.C*cos(angle) + b.S*sin(angle) + b.Z;

      }

      break;

    }
  }

  return GSL_MAX(-limit, GSL_MIN(bound, limit));

}

int XLALSetSuperskyLatticeTilingPhysicalSkyBounds(
  LatticeTiling *tiling,
  gsl_matrix *rssky_metric,
  gsl_matrix *rssky_transf,
  const double alpha1,
  const double alpha2,
  const double delta1,
  const double delta2
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(gsl_matrix_get(rssky_metric, 0, 1) == 0, XLAL_EINVAL);
  XLAL_CHECK(gsl_matrix_get(rssky_metric, 1, 0) == 0, XLAL_EINVAL);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_metric->size1 + 1 == rssky_transf->size1, XLAL_ESIZE);
  XLAL_CHECK(rssky_transf->size2 == 3, XLAL_ESIZE);
  XLAL_CHECK(fabs(alpha1 - alpha2) <= LAL_PI || (fabs(alpha1 - alpha2) >= LAL_TWOPI && fabs(delta1) >= LAL_PI_2 && fabs(delta2) >= LAL_PI_2), XLAL_EINVAL);
  XLAL_CHECK(-LAL_PI_2 <= delta1 && delta1 <= LAL_PI_2, XLAL_EINVAL);
  XLAL_CHECK(-LAL_PI_2 <= delta2 && delta2 <= LAL_PI_2, XLAL_EINVAL);

  // If parameter space is a single point:
  if (alpha1 == alpha2 && delta1 == delta2) {

    // Convert physical point to reduced supersky coordinates A and B
    double rssky_point[rssky_metric->size1];
    gsl_vector_view rssky_point_view = gsl_vector_view_array(rssky_point, rssky_metric->size1);
    PulsarDopplerParams XLAL_INIT_DECL(phys_point);
    phys_point.Alpha = alpha1;
    phys_point.Delta = delta1;
    XLAL_CHECK(XLALConvertPhysicalToSupersky(SC_RSSKY, &rssky_point_view.vector, &phys_point, rssky_transf, &phys_point.refTime) == XLAL_SUCCESS, XLAL_EFUNC);

    // Set the parameter-space bounds on reduced supersky sky coordinates A and B
    for (size_t dim = 0; dim < 2; ++dim) {
      XLAL_CHECK(XLALSetLatticeTilingConstantBound(tiling, dim, rssky_point[dim], rssky_point[dim]) == XLAL_SUCCESS, XLAL_EFUNC);
    }

    return XLAL_SUCCESS;

  }

  // Allocate and initialise bounds data
  const size_t data_len = MAX_PHYS_SKY_BOUND_DATA * sizeof(PhysicalSkyBoundData);
  PhysicalSkyBoundData *data_lower = XLALCalloc(1, data_len);
  XLAL_CHECK(data_lower != NULL, XLAL_ENOMEM);
  PhysicalSkyBoundData *data_upper = XLALCalloc(1, data_len);
  XLAL_CHECK(data_upper != NULL, XLAL_ENOMEM);
  for (size_t i = 0; i < MAX_PHYS_SKY_BOUND_DATA; ++i) {
    data_lower[i].max_A = data_upper[i].max_A = GSL_NEGINF;
  }

  // Special bounds data representing the lower/upper circular bounds on reduced supersky coordinate B
  const PhysicalSkyBoundData lower_circular = { .type = -1 };
  const PhysicalSkyBoundData upper_circular = { .type = 1 };

  // If parameter space is the entire sky:
  if (fabs(alpha1 - alpha2) >= LAL_TWOPI) {

    // Set bounds data to lower/upper circular bounds
    data_lower[0] = lower_circular;
    data_lower[0].max_A = GSL_POSINF;
    data_upper[0] = upper_circular;
    data_upper[0].max_A = GSL_POSINF;

    // Set the parameter-space bounds on reduced supersky sky coordinates A and B
    XLAL_CHECK(XLALSetLatticeTilingConstantBound(tiling, 0, -2, 2) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 1, PhysicalSkyBound, data_len, data_lower, data_upper) == XLAL_SUCCESS, XLAL_EFUNC);

    return XLAL_SUCCESS;

  }

  // Determine the right-handed angle of rotation 'phi' in the reduced supersky
  // plane 'nc = 0' to apply to the matrix
  //   Q^T = [Q_na; Q_nb; Q_nc] = rssky_transf[0:2, 0:2]
  // such that the physical sky point 'alpha = min(alpha1,alpha2), delta = 0' is
  // mapped to the reduced supersky point 'A = B = 0'. This will transform the
  // physical sky region '[alpha1,alpha2] x [delta1,delta2]' such that moving from
  // 'delta1' to 'delta2' will run through the 'A' coordinate, and moving from
  // 'alpha1' to 'alpha2' will roughly run through the 'B' coordinate:
  //          __________________________________________
  //   B = 1_|         _____             _____         |
  //         |      .-'     '-.       .-'     '-.      |
  //         |    .'           '.   .'           '.    |
  //         |   /               \ /               \   |
  //         |  ;                 ;                 ;  |     \.
  //       0-|  |                 |                 |  |      |---> ~direction of delta
  //         |  ;                 ;                 ;  |  <__/
  //         |   \               / \               /   |  ~direction of alpha
  //         |    '.           .'   '.           .'    |
  //        _|      '-._____.-'       '-._____.-'      |
  //      -1 |__.________.________.________.________.__|
  //            '        '        '        '        '
  //       A = -2       -1        0        1        2
  // where '#' indicates the area being tiled.
  double phi = 0;
  {
    const double alpha = GSL_MIN(alpha1, alpha2);
    const double cos_alpha = cos(alpha);
    const double sin_alpha = sin(alpha);
    const double Q_na_0 = gsl_matrix_get(rssky_transf, 0, 0);
    const double Q_na_1 = gsl_matrix_get(rssky_transf, 0, 1);
    const double Q_nb_0 = gsl_matrix_get(rssky_transf, 1, 0);
    const double Q_nb_1 = gsl_matrix_get(rssky_transf, 1, 1);
    const double na = Q_na_0*cos_alpha + Q_na_1*sin_alpha;
    const double nb = Q_nb_0*cos_alpha + Q_nb_1*sin_alpha;
    phi = atan2(-nb, na);
    const double na_rot = na*cos(phi) + nb*sin(phi);
    if (na_rot < 0) {
      phi -= LAL_PI;
    } else if (na_rot > 0) {
      phi += LAL_PI;
    }
  }
  const double cos_phi = cos(phi);
  const double sin_phi = sin(phi);

  // Apply the right-handed rotation matrix
  //   R = [cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1]
  // to the reduced supersky coordinate transform data
  //   rssky_transf = [Q^T; Delta^s]
  // where 'Q' is the sky alignment matrix and 'Delta^s' are the sky offset vectors.
  // The correct transformation to apply is:
  //   Q^T ==> R * Q^T, Delta^s ==> Delta^s . R^T
  for (size_t j = 0; j < 3; ++j) {
    const double Q_na_j = gsl_matrix_get(rssky_transf, 0, j);
    const double Q_nb_j = gsl_matrix_get(rssky_transf, 1, j);
    const double Q_na_j_rot = Q_na_j*cos_phi - Q_nb_j*sin_phi;
    const double Q_nb_j_rot = Q_nb_j*cos_phi + Q_na_j*sin_phi;
    gsl_matrix_set(rssky_transf, 0, j, Q_na_j_rot);
    gsl_matrix_set(rssky_transf, 1, j, Q_nb_j_rot);
  }
  for (size_t i = 3; i < rssky_transf->size1; ++i) {
    const double Delta_0 = gsl_matrix_get(rssky_transf, i, 0);
    const double Delta_1 = gsl_matrix_get(rssky_transf, i, 1);
    const double Delta_0_rot = Delta_0*cos_phi - Delta_1*sin_phi;
    const double Delta_1_rot = Delta_1*cos_phi + Delta_0*sin_phi;
    gsl_matrix_set(rssky_transf, i, 0, Delta_0_rot);
    gsl_matrix_set(rssky_transf, i, 1, Delta_1_rot);
  }

  // Apply 'R' to the sky-sky block of the reduced supersky metric
  //   g_nn = [g_na_na, 0; 0, g_nb_nb]
  // to compensate for the changes to the coordinate transform data.
  // The correct transformation to apply is:
  //   Lambda ==> R * Lambda * R^T
  // where 'Lambda' re-introduces the supressed 'nc' coordinate dimension
  //   Lambda = [g_na_na, 0, 0; 0, g_nb_nb, 0; 0, 0, 0]
  // but since 'R' is a rotation in the 'nc = 0' plane, the metric in
  // this dimension need not be known, and thus can be assumed to be zero.
  {
    const double g_na_na = gsl_matrix_get(rssky_metric, 0, 0);
    const double g_nb_nb = gsl_matrix_get(rssky_metric, 1, 1);
    const double g_na_na_rot = g_na_na*SQR(cos_phi) + g_nb_nb*SQR(sin_phi);
    const double g_na_nb_rot = (g_na_na - g_nb_nb)*cos_phi*sin_phi;
    const double g_nb_nb_rot = g_na_na*SQR(sin_phi) + g_nb_nb*SQR(cos_phi);
    gsl_matrix_set(rssky_metric, 0, 0, g_na_na_rot);
    gsl_matrix_set(rssky_metric, 0, 1, g_na_nb_rot);
    gsl_matrix_set(rssky_metric, 1, 0, g_na_nb_rot);
    gsl_matrix_set(rssky_metric, 1, 1, g_nb_nb_rot);
  }

  // Get components of the vectors 'Q_na', 'Q_nb', and 'Q_nc' from coordinate transform data
  const double Q_na[3] = { gsl_matrix_get(rssky_transf, 0, 0), gsl_matrix_get(rssky_transf, 0, 1), gsl_matrix_get(rssky_transf, 0, 2) };
  const double Q_nb[3] = { gsl_matrix_get(rssky_transf, 1, 0), gsl_matrix_get(rssky_transf, 1, 1), gsl_matrix_get(rssky_transf, 1, 2) };
  const double Q_nc[3] = { gsl_matrix_get(rssky_transf, 2, 0), gsl_matrix_get(rssky_transf, 2, 1), gsl_matrix_get(rssky_transf, 2, 2) };

  // Determine the minimum and maximum right ascension and declination
  const double alphas[2] = { GSL_MIN(alpha1, alpha2), GSL_MAX(alpha1, alpha2) };
  const double deltas[2] = { GSL_MIN(delta1, delta2), GSL_MAX(delta1, delta2) };

  // Create bound data for declination bounds, at constant minimum/maximum right ascension:
  // Given known 'na' from previous bound, and known minimum/maximum 'alpha', solve
  //   na = ( Q_na[0]*cos(alpha) + Q_na[1]*sin(alpha) )*cos(delta) + Q_na[2]*sin(delta)
  // for
  //   delta = asin( ( na - na0 ) / r ) - angle0
  // and compute
  //   nb = C*cos(delta) + S*sin(delta) + Z
  PhysicalSkyBoundData const_alpha[2];
  for (size_t i = 0; i < 2; ++i) {
    XLAL_INIT_MEM(const_alpha[i]);
    const_alpha[i].type = 0;
    const double cos_alpha = cos(alphas[i]);
    const double sin_alpha = sin(alphas[i]);
    const double x = Q_na[2];
    const double y = Q_na[0]*cos_alpha + Q_na[1]*sin_alpha;
    const_alpha[i].na0 = 0;
    const_alpha[i].r = sqrt(SQR(x) + SQR(y));
    const_alpha[i].angle0 = atan2(y, x);
    const_alpha[i].C = Q_nb[0]*cos_alpha + Q_nb[1]*sin_alpha;
    const_alpha[i].S = Q_nb[2];
    const_alpha[i].Z = 0;
  }

  // Create bound data for right ascension bounds, at constant minimum/maximum declination:
  // Given known 'na' from previous bound, and known minimum/maximum 'delta', solve
  //   na = ( Q_na[0]*cos(alpha) + Q_na[1]*sin(alpha) )*cos(delta) + Q_na[2]*sin(delta)
  // for
  //   alpha = asin( ( na - na0 ) / r ) - angle0
  // and compute
  //   nb = C*cos(alpha) + S*sin(alpha) + Z
  PhysicalSkyBoundData const_delta[2];
  for (size_t j = 0; j < 2; ++j) {
    XLAL_INIT_MEM(const_delta[j]);
    const_delta[j].type = 0;
    const double cos_delta = cos(deltas[j]);
    const double sin_delta = sin(deltas[j]);
    const double x = Q_na[1]*cos_delta;
    const double y = Q_na[0]*cos_delta;
    const_delta[j].na0 = Q_na[2]*sin_delta;
    const_delta[j].r = sqrt(SQR(x) + SQR(y));
    const_delta[j].angle0 = atan2(y, x);
    const_delta[j].C = Q_nb[0]*cos_delta;
    const_delta[j].S = Q_nb[1]*cos_delta;
    const_delta[j].Z = Q_nb[2]*sin_delta;
  }

  // Determine corner points in reduced supersky coordinate A of the region
  //   '[min(alpha),max(alpha)] x [min(delta),max(delta)]'
  double corner_A[2][2];
  for (size_t i = 0; i < 2; ++i) {
    for (size_t j = 0; j < 2; ++j) {
      double rssky_point[rssky_metric->size1];
      gsl_vector_view rssky_point_view = gsl_vector_view_array(rssky_point, rssky_metric->size1);
      PulsarDopplerParams XLAL_INIT_DECL(phys_point);
      phys_point.Alpha = alphas[i];
      phys_point.Delta = deltas[j];
      XLAL_CHECK(XLALConvertPhysicalToSupersky(SC_RSSKY, &rssky_point_view.vector, &phys_point, rssky_transf, &phys_point.refTime) == XLAL_SUCCESS, XLAL_EFUNC);
      corner_A[i][j] = rssky_point[0];
    }
  }

  // Use corner points to classify parameter space into different shapes and set bounds data
  double min_A = GSL_NEGINF, max_A = GSL_POSINF;
  if (corner_A[1][0] < 0 && corner_A[1][1] <= 0) {

    if (corner_A[1][1] < corner_A[1][0]) {

      //          __________________________________________  Lower bound(s) on B:
      //   B = 1_|         _____             _____         |  0 = right ascension bound at max(delta) until max_A
      //         |      .-'     '-.       .-'     '-.      |
      //         |    .'           '.   .'           '.    |  Upper bound(s) on B:
      //         |   /               \ /               \   |  0 = declination bound at max(alpha) until corner_A[1][0]
      //         |  ;          __2_   ;                 ;  |  1 = right ascension bound at min(delta) until corner_A[0][0]
      //       0-|  |         /####|  |                 |  |  2 = declination bound at min(alpha) until max_A
      //         |  ;      _1-#####;  ;                 ;  |
      //         |   \   0/#######/  / \               /   |
      //         |    '. /###0##.' .'   '.           .'    |
      //        _|      '-._____.-'       '-._____.-'      |
      //      -1 |__.________.________.________.________.__|
      //            '        '        '        '        '
      //       A = -2       -1        0        1        2
      min_A = corner_A[1][1];
      max_A = corner_A[0][1];
      data_lower[0] = const_delta[1];
      data_lower[0].max_A = GSL_POSINF;
      data_upper[0] = const_alpha[1];
      data_upper[0].max_A = corner_A[1][0];
      data_upper[1] = const_delta[0];
      data_upper[1].max_A = corner_A[0][0];
      data_upper[2] = const_alpha[0];
      data_upper[2].max_A = GSL_POSINF;
      data_upper[2].altroot = true;

    } else {

      //          __________________________________________  Lower bound(s) on B:
      //   B = 1_|         _____             _____         |  0 = declination bound at max(alpha) until corner_A[1][1]
      //         |      .-'     '-.       .-'     '-.      |  1 = right ascension bound at max(delta) until max_A
      //         |    .'           '.   .'           '.    |
      //         |   /               \ /               \   |  Upper bound(s) on B:
      //         |  ;          __1_   ;                 ;  |  0 = right ascension bound at min(delta) until corner_A[0][0]
      //       0-|  |        0/####|  |                 |  |  1 = declination bound at min(alpha) until max_A
      //         |  ;        -#####;  ;                 ;  |
      //         |   \      0\####/  / \               /   |
      //         |    '.      \1.' .'   '.           .'    |
      //        _|      '-._____.-'       '-._____.-'      |
      //      -1 |__.________.________.________.________.__|
      //            '        '        '        '        '
      //       A = -2       -1        0        1        2
      min_A = corner_A[1][0];
      max_A = corner_A[0][1];
      data_lower[0] = const_alpha[1];
      data_lower[0].max_A = corner_A[1][1];
      data_lower[0].altroot = true;
      data_lower[1] = const_delta[1];
      data_lower[1].max_A = GSL_POSINF;
      data_upper[0] = const_delta[0];
      data_upper[0].max_A = corner_A[0][0];
      data_upper[1] = const_alpha[0];
      data_upper[1].max_A = GSL_POSINF;
      data_upper[1].altroot = true;

    }

  } else if (0 <= corner_A[1][0] && 0 < corner_A[1][1]) {

    if (corner_A[1][1] < corner_A[1][0]) {

      //          __________________________________________  Lower bound(s) on B:
      //   B = 1_|         _____             _____         |  0 = right ascension bound at min(delta) until max_A
      //         |      .-'     '-.       .-'     '-.      |
      //         |    .'           '.   .'           '.    |  Upper bound(s) on B:
      //         |   /               \ /               \   |  0 = declination bound at min(alpha) until corner_A[0][1]
      //         |  ;                 ;   _0__          ;  |  1 = right ascension bound at max(delta) until corner_A[1][1x]
      //       0-|  |                 |  |####\         |  |  2 = declination bound at max(alpha) until max_A
      //         |  ;                 ;  ;#####-1_      ;  |
      //         |   \               / \  \#######\2   /   |
      //         |    '.           .'   '. '.##0###\ .'    |
      //        _|      '-._____.-'       '-._____.-'      |
      //      -1 |__.________.________.________.________.__|
      //            '        '        '        '        '
      //       A = -2       -1        0        1        2
      min_A = corner_A[0][0];
      max_A = corner_A[1][0];
      data_lower[0] = const_delta[0];
      data_lower[0].max_A = GSL_POSINF;
      data_upper[0] = const_alpha[0];
      data_upper[0].max_A = corner_A[0][1];
      data_upper[1] = const_delta[1];
      data_upper[1].max_A = corner_A[1][1];
      data_upper[2] = const_alpha[1];
      data_upper[2].max_A = GSL_POSINF;
      data_upper[2].altroot = true;

    } else {

      //          __________________________________________  Lower bound(s) on B:
      //   B = 1_|         _____             _____         |  0 = right ascension bound at min(delta) until corner_A[1][0]
      //         |      .-'     '-.       .-'     '-.      |  1 = declination bound at max(alpha) until max_A
      //         |    .'           '.   .'           '.    |
      //         |   /               \ /               \   |  Upper bound(s) on B:
      //         |  ;                 ;   _0__          ;  |  0 = declination bound at min(alpha) until corner_A[0][1]
      //       0-|  |                 |  |####\1        |  |  1 = right ascension bound at max(delta) until max_A
      //         |  ;                 ;  ;#####-        ;  |
      //         |   \               / \  \####/1      /   |
      //         |    '.           .'   '. \0.'      .'    |
      //        _|      '-._____.-'       '-._____.-'      |
      //      -1 |__.________.________.________.________.__|
      //            '        '        '        '        '
      //       A = -2       -1        0        1        2
      min_A = corner_A[0][0];
      max_A = corner_A[1][1];
      data_lower[0] = const_delta[0];
      data_lower[0].max_A = corner_A[1][0];
      data_lower[1] = const_alpha[1];
      data_lower[1].max_A = GSL_POSINF;
      data_upper[0] = const_alpha[0];
      data_upper[0].max_A = corner_A[0][1];
      data_upper[1] = const_delta[1];
      data_upper[1].max_A = GSL_POSINF;

    }

  } else {

    // This parameter space straddles both reduced supersky hemispheres
    // Find the value of 'na' where this occurs at max(alpha) by solving
    //   nc = ( Q_nc[0]*cos(max(alpha)) + Q_nc[1]*sin(max(alpha)) )*cos(delta) + Q_nc[2]*sin(delta)
    // for delta, then computing
    //   na = ( Q_na[0]*cos(alpha) + Q_na[1]*sin(alpha) )*cos(delta) + Q_na[2]*sin(delta)
    const double cos_alpha_split = cos(alphas[1]);
    const double sin_alpha_split = sin(alphas[1]);
    const double delta_split = -1 * atan2(Q_nc[0]*cos_alpha_split + Q_nc[1]*sin_alpha_split, Q_nc[2]);
    const double na_split = (Q_na[0]*cos_alpha_split + Q_na[1]*sin_alpha_split)*cos(delta_split) + Q_na[2]*sin(delta_split);
    const double split_A[2] = { -1 - na_split, 1 + na_split };

    if (split_A[0] < corner_A[1][0]) {
      if (corner_A[1][1] < split_A[1]) {

        //          __________________________________________  Lower bound(s) on B:
        //   B = 1_|         _____             _____         |  0 = lower circular bound until max_A
        //         |      .-'     '-.       .-'     '-.      |
        //         |    .'           '.   .'           '.    |  Upper bound(s) on B:
        //         |   /               \ /               \   |  0 = declination bound at max(alpha) until corner_A[1][0]
        //         |  ;            __2__;___3___          ;  |  1 = right ascension bound at min(delta) until corner_A[0][0]
        //       0-|  |           /#####|#######\         |  |  2 = declination bound at min(alpha) until A = 0
        //         |  ;      ___1-######;########-4_      ;  |  3 = declination bound at min(alpha) until corner_A[0][1]
        //         |   \    /##########/ \##########\    /   |  4 = right ascension bound at max(delta) until corner_A[1][1]
        //         |    '.0/#########.'   '.#########\5.'    |  5 = declination bound at max(alpha) until max_A
        //        _|      '-.##0##.-'       '-.##0##.-'      |
        //      -1 |__.________.________.________.________.__|
        //            '        '        '        '        '
        //       A = -2       -1        0        1        2
        min_A = split_A[0];
        max_A = split_A[1];
        data_lower[0] = lower_circular;
        data_lower[0].max_A = GSL_POSINF;
        data_upper[0] = const_alpha[1];
        data_upper[0].max_A = corner_A[1][0];
        data_upper[1] = const_delta[0];
        data_upper[1].max_A = corner_A[0][0];
        data_upper[2] = const_alpha[0];
        data_upper[2].max_A = 0;
        data_upper[2].altroot = true;
        data_upper[3] = const_alpha[0];
        data_upper[3].max_A = corner_A[0][1];
        data_upper[4] = const_delta[1];
        data_upper[4].max_A = corner_A[1][1];
        data_upper[5] = const_alpha[1];
        data_upper[5].max_A = GSL_POSINF;
        data_upper[5].altroot = true;

      } else {

        //          __________________________________________  Lower bound(s) on B:
        //   B = 1_|         _____             _____         |  0 = lower circular bound until split_A[1]
        //         |      .-'     '-.       .-'     '-.      |  1 = declination bound at max(alpha) until max_A
        //         |    .'           '.   .'           '.    |
        //         |   /               \ /               \   |  Upper bound(s) on B:
        //         |  ;            __2__;___3___          ;  |  0 = declination bound at max(alpha) until corner_A[1][0]
        //       0-|  |           /#####|#######\         |  |  1 = right ascension bound at min(delta) until corner_A[0][0]
        //         |  ;      ___1-######;########-4_      ;  |  2 = declination bound at min(alpha) until A = 0
        //         |   \    /##########/ \##########/    /   |  3 = declination bound at min(alpha) until corner_A[0][1]
        //         |    '.0/#########.'   '.#######/1  .'    |  4 = right ascension bound at max(delta) until max_A
        //        _|      '-.##0##.-'       '-.#0#/_.-'      |
        //      -1 |__.________.________.________.________.__|
        //            '        '        '        '        '
        //       A = -2       -1        0        1        2
        min_A = split_A[0];
        max_A = corner_A[1][1];
        data_lower[0] = lower_circular;
        data_lower[0].max_A = split_A[1];
        data_lower[1] = const_alpha[1];
        data_lower[1].max_A = GSL_POSINF;
        data_upper[0] = const_alpha[1];
        data_upper[0].max_A = corner_A[1][0];
        data_upper[1] = const_delta[0];
        data_upper[1].max_A = corner_A[0][0];
        data_upper[2] = const_alpha[0];
        data_upper[2].max_A = 0;
        data_upper[2].altroot = true;
        data_upper[3] = const_alpha[0];
        data_upper[3].max_A = corner_A[0][1];
        data_upper[4] = const_delta[1];
        data_upper[4].max_A = GSL_POSINF;

      }

    } else {
      if (corner_A[1][1] < split_A[1]) {

        //          __________________________________________  Lower bound(s) on B:
        //   B = 1_|         _____             _____         |  0 = declination bound at max(alpha) until split_A[0]
        //         |      .-'     '-.       .-'     '-.      |  1 = lower circular bound until split_A[1]
        //         |    .'           '.   .'           '.    |
        //         |   /               \ /               \   |  Upper bound(s) on B:
        //         |  ;            __1__;___2___          ;  |  0 = right ascension bound at min(delta) until corner_A[0][0]
        //       0-|  |           /#####|#######\         |  |  1 = declination bound at min(alpha) until A = 0
        //         |  ;      ___0-######;########-3_      ;  |  2 = declination bound at min(alpha) until corner_A[0][1]
        //         |   \    \##########/ \##########\    /   |  3 = right ascension bound at max(delta) until corner_A[1][1]
        //         |    '.  0\#######.'   '.#########\4.'    |  4 = declination bound at max(alpha) until max_A
        //        _|      '-._'.#1.-'       '-.##1##.-'      |
        //      -1 |__.________.________.________.________.__|
        //            '        '        '        '        '
        //       A = -2       -1        0        1        2
        min_A = corner_A[1][0];
        max_A = split_A[1];
        data_lower[0] = const_alpha[1];
        data_lower[0].max_A = split_A[0];
        data_lower[0].altroot = true;
        data_lower[1] = lower_circular;
        data_lower[1].max_A = GSL_POSINF;
        data_upper[0] = const_delta[0];
        data_upper[0].max_A = corner_A[0][0];
        data_upper[1] = const_alpha[0];
        data_upper[1].max_A = 0;
        data_upper[1].altroot = true;
        data_upper[2] = const_alpha[0];
        data_upper[2].max_A = corner_A[0][1];
        data_upper[3] = const_delta[1];
        data_upper[3].max_A = corner_A[1][1];
        data_upper[4] = const_alpha[1];
        data_upper[4].max_A = GSL_POSINF;
        data_upper[4].altroot = true;

      } else {

        //          __________________________________________  Lower bound(s) on B:
        //   B = 1_|         _____             _____         |  0 = declination bound at max(alpha) until split_A[0]
        //         |      .-'     '-.       .-'     '-.      |  1 = lower circular bound until split_A[1]
        //         |    .'           '.   .'           '.    |  2 = declination bound at max(alpha) until max_A
        //         |   /               \ /               \   |
        //         |  ;            __1__;___2___          ;  |  Upper bound(s) on B:
        //       0-|  |           /#####|#######\         |  |  0 = right ascension bound at min(delta) until corner_A[0][0]
        //         |  ;      ___0-######;########-3_      ;  |  1 = declination bound at min(alpha) until A = 0
        //         |   \    \##########/ \##########/    /   |  2 = declination bound at min(alpha) until corner_A[0][1]
        //         |    '.  0\#######.'   '.#######/2  .'    |  3 = right ascension bound at max(delta) until max_A
        //        _|      '-._'.#1.-'       '-.#1.'_.-'      |
        //      -1 |__.________.________.________.________.__|
        //            '        '        '        '        '
        //       A = -2       -1        0        1        2
        min_A = corner_A[1][0];
        max_A = corner_A[1][1];
        data_lower[0] = const_alpha[1];
        data_lower[0].max_A = split_A[0];
        data_lower[0].altroot = true;
        data_lower[1] = lower_circular;
        data_lower[1].max_A = split_A[1];
        data_lower[2] = const_alpha[1];
        data_lower[2].max_A = GSL_POSINF;
        data_upper[0] = const_delta[0];
        data_upper[0].max_A = corner_A[0][0];
        data_upper[1] = const_alpha[0];
        data_upper[1].max_A = 0;
        data_upper[1].altroot = true;
        data_upper[2] = const_alpha[0];
        data_upper[2].max_A = corner_A[0][1];
        data_upper[3] = const_delta[1];
        data_upper[3].max_A = GSL_POSINF;

      }
    }

  }


  // Set the parameter-space bounds on reduced supersky sky coordinates A and B
  XLAL_CHECK(XLALSetLatticeTilingConstantBound(tiling, 0, min_A, max_A) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, 1, PhysicalSkyBound, data_len, data_lower, data_upper) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

int XLALSetSuperskyLatticeTilingPhysicalSkyPatch(
  LatticeTiling *tiling,
  gsl_matrix *rssky_metric,
  gsl_matrix *rssky_transf,
  const UINT4 patch_count,
  const UINT4 patch_index
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(gsl_matrix_get(rssky_metric, 0, 1) == 0, XLAL_EINVAL);
  XLAL_CHECK(gsl_matrix_get(rssky_metric, 1, 0) == 0, XLAL_EINVAL);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_metric->size1 + 1 == rssky_transf->size1, XLAL_ESIZE);
  XLAL_CHECK(rssky_transf->size2 == 3, XLAL_ESIZE);
  XLAL_CHECK(patch_count > 0, XLAL_EINVAL);
  XLAL_CHECK(patch_index < patch_count, XLAL_EINVAL);

  // Number of patch divisions in 'alpha'; for less than 4 patches, divide only in 'alpha' to prevent
  // 'alpha' range in [pi,2*pi], which XLALSetSuperskyLatticeTilingPhysicalSkyBounds() cannot handle
  const UINT4 alpha_count = (patch_count < 4) ? patch_count : ((UINT4) ceil(sqrt(patch_count)));

  // Mininum number of patch divisions in 'sin(delta)'; note integer division equivalent to floor()
  const UINT4 min_sdelta_count = patch_count / alpha_count;

  // Excess number of patches, which must be added on to get 'patch_count'
  INT4 patch_excess = patch_count - alpha_count * min_sdelta_count;
  XLAL_CHECK(patch_excess >= 0, XLAL_EFAILED);

  // Initialise number of patch divisions in 'sin(delta)'; if there are excess patches, add an extra patch
  UINT4 sdelta_count = min_sdelta_count;
  if (patch_excess > 0) {
    ++sdelta_count;
  }

  // Calculate range of indices in 'alpha', and number of patch divisions and index in 'sin(delta)'.
  // The divisions in 'alpha' are set in proportion to the range of 'alpha_index', i.e. the number of
  // divisions in 'sin(delta)' for that range of 'alpha_index'. This is so that, if 'patch_excess' is
  // not zero, and therefore the number of divisions in 'sin(delta)' is not constant, patch areas should
  // still be equal. Example:
  //   patch_count=7 patch_index=0 | alpha_index=0--3 sdelta_count=3 sdelta_index=0
  //   patch_count=7 patch_index=1 | alpha_index=0--3 sdelta_count=3 sdelta_index=1
  //   patch_count=7 patch_index=2 | alpha_index=0--3 sdelta_count=3 sdelta_index=2
  //   patch_count=7 patch_index=3 | alpha_index=3--5 sdelta_count=2 sdelta_index=0
  //   patch_count=7 patch_index=4 | alpha_index=3--5 sdelta_count=2 sdelta_index=1
  //   patch_count=7 patch_index=5 | alpha_index=5--7 sdelta_count=2 sdelta_index=0
  //   patch_count=7 patch_index=6 | alpha_index=5--7 sdelta_count=2 sdelta_index=1
  UINT4 alpha_index1 = 0, alpha_index2 = sdelta_count, sdelta_index = patch_index;
  while (sdelta_index >= sdelta_count) {

    // Decrease index in 'sin(delta)'; we are done when 'sdelta_index' < 'sdelta_count'
    sdelta_index -= sdelta_count;

    // Decrease number of excess patches; if zero, subtract extra patch from patch divisions in 'sin(delta)'
    --patch_excess;
    if (patch_excess == 0) {
      --sdelta_count;
    }

    // Store the current last 'alpha' index in 'alpha_index1', and increase
    // 'alpha_index2' by the current number of patch divisions in 'sin(delta)'
    alpha_index1 = alpha_index2;
    alpha_index2 += sdelta_count;

  }

  // Compute range of 'alpha' to bound
  const double alpha1 = LAL_TWOPI * ((double) alpha_index1) / ((double) patch_count);
  const double alpha2 = LAL_TWOPI * ((double) alpha_index2) / ((double) patch_count);

  // Compute range of 'sin(delta)' to bound
  const double sdelta1 = -1 + 2 * ((double) sdelta_index) / ((double) sdelta_count);
  const double sdelta2 = -1 + 2 * ((double) sdelta_index + 1) / ((double) sdelta_count);

  // Set the parameter-space bounds on physical sky position 'alpha' and 'delta'
  XLAL_CHECK(XLALSetSuperskyLatticeTilingPhysicalSkyBounds(tiling, rssky_metric, rssky_transf, alpha1, alpha2, asin(sdelta1), asin(sdelta2)) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

static double PhysicalSpinBound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_vector *point
  )
{

  // Get bounds data
  const double *sky_offsets = ((const double *) data);
  double bound = ((const double *) data)[3];

  // Add the inner product of the sky offsets with the aligned sky
  // position to the physical bound to get the reduced supersky bound
  double as[3];
  SM_ReducedToAligned(as, point);
  bound += DOT3(sky_offsets, as);

  return bound;

}

int XLALSetSuperskyLatticeTilingPhysicalSpinBound(
  LatticeTiling *tiling,
  const gsl_matrix *rssky_transf,
  const size_t s,
  const double bound1,
  const double bound2
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf->size1 > 3, XLAL_ESIZE);
  XLAL_CHECK(rssky_transf->size2 == 3, XLAL_ESIZE);
  XLAL_CHECK(isfinite(bound1), XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound2), XLAL_EINVAL);
  const size_t smax = rssky_transf->size1 - 4;
  XLAL_CHECK(s <= smax, XLAL_ESIZE);
  const size_t dim = (s == 0) ? (2 + smax) : (1 + s);

  // Allocate memory
  const size_t data_len = 4 * sizeof(double);
  double *data_lower = XLALMalloc(data_len);
  XLAL_CHECK(data_lower != NULL, XLAL_ENOMEM);
  double *data_upper = XLALMalloc(data_len);
  XLAL_CHECK(data_upper != NULL, XLAL_ENOMEM);

  // Copy the sky offset vector to bounds data
  for (size_t j = 0; j < 3; ++j) {
    data_lower[j] = data_upper[j] = gsl_matrix_get(rssky_transf, dim + 1, j);
  }

  // Set the parameter-space bound on physical frequency/spindown coordinate
  data_lower[3] = GSL_MIN(bound1, bound2);
  data_upper[3] = GSL_MAX(bound1, bound2);
  XLAL_CHECK(XLALSetLatticeTilingBound(tiling, dim, PhysicalSpinBound, data_len, data_lower, data_upper) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

int XLALSetSuperskyLatticeTilingCoordinateSpinBound(
  LatticeTiling *tiling,
  const gsl_matrix *rssky_transf,
  const size_t s,
  const double bound1,
  const double bound2
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf->size1 > 3, XLAL_ESIZE);
  XLAL_CHECK(rssky_transf->size2 == 3, XLAL_ESIZE);
  XLAL_CHECK(isfinite(bound1), XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound2), XLAL_EINVAL);
  const size_t smax = rssky_transf->size1 - 4;
  XLAL_CHECK(s <= smax, XLAL_ESIZE);
  const size_t dim = (s == 0) ? (2 + smax) : (1 + s);

  // Set the parameter-space bound on reduced supersky frequency/spindown coordinate
  XLAL_CHECK(XLALSetLatticeTilingConstantBound(tiling, dim, bound1, bound2) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

int XLALSuperskyLatticePulsarSpinRange(
  PulsarSpinRange *spin_range,
  LatticeTiling *tiling,
  const gsl_matrix *rssky_transf,
  const LIGOTimeGPS *ref_time
  )
{

  // Check input
  XLAL_CHECK(spin_range != NULL, XLAL_EFAULT);
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf->size1 > 3, XLAL_ESIZE);
  XLAL_CHECK(rssky_transf->size2 == 3, XLAL_ESIZE);
  XLAL_CHECK(ref_time != NULL, XLAL_EFAULT);

  // Get rectange containing range reduced supersky coordinates
  double skyA_rect[4], skyB_rect[4];
  {
    const LatticeTilingStats *stats = XLALLatticeTilingStatistics(tiling, 0);
    XLAL_CHECK(stats != NULL, XLAL_EFUNC);
    skyA_rect[0] = stats->min_value_pass;
    skyA_rect[1] = stats->max_value_pass;
    skyA_rect[2] = stats->max_value_pass;
    skyA_rect[3] = stats->min_value_pass;
  }
  {
    const LatticeTilingStats *stats = XLALLatticeTilingStatistics(tiling, 1);
    XLAL_CHECK(stats != NULL, XLAL_EFUNC);
    skyB_rect[0] = stats->min_value_pass;
    skyB_rect[1] = stats->min_value_pass;
    skyB_rect[2] = stats->max_value_pass;
    skyB_rect[3] = stats->max_value_pass;
  }

  // Get range of physical frequency/spindowns
  const size_t smax = rssky_transf->size1 - 4;
  for (size_t i = 0; i < 4; ++i) {

    // Construct reduced supersky point
    double in_rssky[3 + smax];
    in_rssky[0] = skyA_rect[i];
    in_rssky[1] = skyB_rect[i];
    for (size_t s = 0; s <= smax; ++s) {
      const size_t dim = (s == 0) ? (2 + smax) : (1 + s);
      const LatticeTilingStats *stats = XLALLatticeTilingStatistics(tiling, dim);
      XLAL_CHECK(stats != NULL, XLAL_EFUNC);
      in_rssky[dim] = stats->min_value_pass;
    }

    // Convert reduced supersky point to physical coordinates
    gsl_vector_view in_rssky_view = gsl_vector_view_array(in_rssky, 3 + smax);
    PulsarDopplerParams XLAL_INIT_DECL(out_phys);
    XLAL_CHECK(XLALConvertSuperskyToPhysical(&out_phys, SC_RSSKY, &in_rssky_view.vector, rssky_transf, ref_time) == XLAL_SUCCESS, XLAL_EFUNC);

    // Store minimum/maximum physical frequency/spindown in 'spin_range'
    for (size_t s = 0; s <= smax; ++s) {
      if (i == 0 || out_phys.fkdot[s] < spin_range->fkdot[s]) {
        spin_range->fkdot[s] = out_phys.fkdot[s];
      }
      if (i == 0 || out_phys.fkdot[s] > spin_range->fkdotBand[s]) {
        spin_range->fkdotBand[s] = out_phys.fkdot[s];
      }
    }

  }
  for (size_t s = 0; s <= smax; ++s) {
    spin_range->fkdotBand[s] -= spin_range->fkdot[s];
  }

  // Adjust 'spin_range' bands to include width of supersky frequency/spindown parameter space
  for (size_t s = 0; s <= smax; ++s) {
    const size_t dim = (s == 0) ? (2 + smax) : (1 + s);
    const LatticeTilingStats *stats = XLALLatticeTilingStatistics(tiling, dim);
    XLAL_CHECK(stats != NULL, XLAL_EFUNC);
    spin_range->fkdotBand[s] += stats->max_value_pass - stats->min_value_pass;
  }

  // Set reference time of 'spin_range' to that of coordinate transform data
  spin_range->refTime = *ref_time;

  return XLAL_SUCCESS;

}
