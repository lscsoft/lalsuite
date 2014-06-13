//
// Copyright (C) 2014 Karl Wette
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

#include <lal/SuperSkyMetrics.h>
#include <lal/LALConstants.h>

#include "GSLHelpers.h"

///
/// Ensure a matrix is exactly symmetric, by ensuring it is equal to its transpose.
///
static void SSM_MakeSymmetric(
  gsl_matrix* matrix				///< [in/out] Matrix to make symmetric
  )
{
  for (size_t i = 0; i < matrix->size1; ++i) {
    for (size_t j = i + 1; j < matrix->size2; ++j) {
      const double mij = gsl_matrix_get(matrix, i, j);
      const double mji = gsl_matrix_get(matrix, j, i);
      const double m = 0.5 * (mij + mji);
      gsl_matrix_set(matrix, i, j, m);
      gsl_matrix_set(matrix, j, i, m);
    }
  }
}

///
/// Build matrix to reconstruct the super-sky metric in equatorial coordinates from the expanded super-sky metric.
///
static void SSM_ReconstructionMatrix(
  gsl_matrix* reconstruct			///< [in/out] Reconstruction matrix
  )
{
  const size_t fsize = reconstruct->size1 - 5;
  gsl_matrix_set_zero(reconstruct);
  const double reconstruct_ss_val[5][3] = {
    {1, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, LAL_COSIEARTH, LAL_SINIEARTH}, {0, -LAL_SINIEARTH, LAL_COSIEARTH}
  };
  gsl_matrix_const_view reconstruct_ss_val_view = gsl_matrix_const_view_array(&reconstruct_ss_val[0][0], 5, 3);
  gsl_matrix_view reconstruct_ss = gsl_matrix_submatrix(reconstruct, 0, 0, 5, 3);
  gsl_matrix_memcpy(&reconstruct_ss.matrix, &reconstruct_ss_val_view.matrix);
  gsl_matrix_view reconstruct_ff = gsl_matrix_submatrix(reconstruct, 5, 3, fsize, fsize);
  gsl_matrix_set_identity(&reconstruct_ff.matrix);
}

///
/// Diagonal-normalise the given matrix, and return the normalisation transform.
///
static void SSM_DiagonalNormalise(
  gsl_matrix* matrix,				///< [in/out] Matrix to diagonal-normalise
  gsl_matrix* transf				///< [in/out] Normalisation transform
  )
{
  for (size_t i = 0; i < matrix->size1; ++i) {
    const double dii = sqrt(fabs(gsl_matrix_get(matrix, i, i)));
    gsl_matrix_set(transf, i, i, dii > 0 ? 1.0/dii : 1.0);
  }
  for (size_t i = 0; i < matrix->size1; ++i) {
    const double sii = gsl_matrix_get(transf, i, i);
    for (size_t j = 0; j < matrix->size2; ++j) {
      const double sjj = gsl_matrix_get(transf, j, j);
      const double mij = gsl_matrix_get(matrix, i, j);
      gsl_matrix_set(matrix, i, j, mij * sii * sjj);
    }
  }
}

int XLALExpandedSuperSkyMetric(
  gsl_matrix **essky_metric,
  const size_t spindowns,
  const LIGOTimeGPS* ref_time,
  const LALSegList* segments,
  const double fiducial_freq,
  const MultiLALDetector* detectors,
  const MultiNoiseFloor* detector_weights,
  const DetectorMotionType detector_motion,
  const EphemerisData* ephemerides
  )
{

  // Check input
  XLAL_CHECK(essky_metric != NULL && *essky_metric == NULL, XLAL_EFAULT);
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

  // Create parameters struct for XLALDopplerFstatMetric()
  DopplerMetricParams XLAL_INIT_DECL(par);

  // Set expanded super-sky coordinate system: x,y spin motion in equatorial coordinates,
  // then X,Y,Z orbital motion in ecliptic coordinates, then spindowns, then frequency
  {
    size_t i = 0;
    par.coordSys.coordIDs[i++] = DOPPLERCOORD_N3SX_EQU;
    par.coordSys.coordIDs[i++] = DOPPLERCOORD_N3SY_EQU;
    par.coordSys.coordIDs[i++] = DOPPLERCOORD_N3OX_ECL;
    par.coordSys.coordIDs[i++] = DOPPLERCOORD_N3OY_ECL;
    par.coordSys.coordIDs[i++] = DOPPLERCOORD_N3OZ_ECL;
    if (spindowns >= 1) {
      par.coordSys.coordIDs[i++] = DOPPLERCOORD_F1DOT;
    }
    if (spindowns >= 2) {
      par.coordSys.coordIDs[i++] = DOPPLERCOORD_F2DOT;
    }
    if (spindowns >= 3) {
      par.coordSys.coordIDs[i++] = DOPPLERCOORD_F3DOT;
    }
    par.coordSys.coordIDs[i++] = DOPPLERCOORD_FREQ;
    par.coordSys.dim = i;
  }

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

  // Only compute the phase metric
  par.metricType = METRIC_TYPE_PHASE;

  // Do not include sky-position-dependent Roemer delay in time variable
  par.approxPhase = 1;

  // Allow metric to have at most 1 non-positive eigenvalue
  par.nonposEigValThresh = 2;

  // Call XLALDopplerFstatMetric() and extract phase metric
  DopplerMetric* retn = XLALDopplerFstatMetric(&par, ephemerides);
  XLAL_CHECK(retn != NULL, XLAL_EFUNC, "XLALDopplerFstatMetric() failed");
  *essky_metric = retn->g_ij;
  retn->g_ij = NULL;
  XLALDestroyDopplerMetric(retn);

  // Ensured returned metric is symmetric
  SSM_MakeSymmetric(*essky_metric);

  return XLAL_SUCCESS;

}

int XLALSuperSkyMetric(
  gsl_matrix **ssky_metric,
  const gsl_matrix* essky_metric
  )
{

  // Check input
  XLAL_CHECK(ssky_metric != NULL && *ssky_metric == NULL, XLAL_EFAULT);
  XLAL_CHECK(essky_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(essky_metric->size1 == essky_metric->size2, XLAL_EINVAL);
  XLAL_CHECK(essky_metric->size1 > 5, XLAL_EINVAL);

  // Size of the frequency+spindowns block of the expanded super-sky metric
  const size_t fsize = essky_metric->size1 - 5;

  // Allocate memory
  GAMAT(*ssky_metric, 3 + fsize, 3 + fsize);
  gsl_matrix* GAMAT(reconstruct, 5 + fsize, 3 + fsize);
  gsl_matrix* GAMAT(tmp, 5 + fsize, 3 + fsize);

  // Build matrix to reconstruct the super-sky metric from the expanded super-sky metric
  SSM_ReconstructionMatrix(reconstruct);

  // Reconstruct the super-sky metric from the expanded super-sky metric:
  //   *ssky_metric = reconstruct^T * essky_metric * reconstruct
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, essky_metric, reconstruct, 0.0, tmp);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, reconstruct, tmp, 0.0, *ssky_metric);

  // Ensured returned metric is symmetric
  SSM_MakeSymmetric(*ssky_metric);

  // Cleanup
  GFMAT(reconstruct, tmp);

  return XLAL_SUCCESS;

}

int XLALReducedSuperSkyMetric(
  gsl_matrix **rssky_metric,
  gsl_matrix **rssky_transf,
  const gsl_matrix* essky_metric
  )
{

  // Check input
  XLAL_CHECK(rssky_metric != NULL && *rssky_metric == NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL && *rssky_transf == NULL, XLAL_EFAULT);
  XLAL_CHECK(essky_metric != NULL, XLAL_EFAULT);
  XLAL_CHECK(essky_metric->size1 == essky_metric->size2, XLAL_EINVAL);
  XLAL_CHECK(essky_metric->size1 > 5, XLAL_EINVAL);

  // Size of the frequency+spindowns block of the expanded super-sky metric
  const size_t fsize = essky_metric->size1 - 5;

  // Allocate memory
  GAMAT(*rssky_metric, 2 + fsize, 2 + fsize);
  GAMAT(*rssky_transf, 3 + fsize, 3);
  gsl_matrix* GAMAT(intm_ssky_metric, 3 + fsize, 3 + fsize);

  //
  // Find least-squares linear fit to orbital X and Y by frequency and spindowns
  // Produces the intermediate "residual" super-sky metric
  //

  {

    // Allocate memory
    gsl_matrix* GAMAT(essky_dnorm, 5 + fsize, 5 + fsize);
    gsl_matrix* GAMAT(essky_dnorm_transf, 5 + fsize, 5 + fsize);
    gsl_matrix* GAMAT(reconstruct, 5 + fsize, 3 + fsize);
    gsl_matrix* GAMAT(residual, 5 + fsize, 3 + fsize);
    gsl_matrix* GAMAT(subtractfit, 5 + fsize, 5 + fsize);
    gsl_matrix* GAMAT(tmp, 5 + fsize, 3 + fsize);

    // Build matrix to reconstruct the super-sky metric from the expanded super-sky metric
    SSM_ReconstructionMatrix(reconstruct);

    // Diagonal-normalise the expanded super-sky metric
    gsl_matrix_memcpy(essky_dnorm, essky_metric);
    SSM_DiagonalNormalise(essky_dnorm, essky_dnorm_transf);

    gsl_matrix* fitA = NULL;
    gsl_matrix* fitAt_fitA = NULL;
    gsl_matrix* fitAt_fitA_svd_AU = NULL;
    gsl_matrix* fitAt_fitA_svd_V = NULL;
    gsl_vector* fitAt_fitA_svd_S = NULL;
    gsl_vector* tmpv = NULL;
    for (size_t n = fsize; n > 0; --n) {

      // Allocate memory
      GAMAT(fitA, 5 + fsize, n);
      GAMAT(fitAt_fitA, n, n);
      GAMAT(fitAt_fitA_svd_AU, n, n);
      GAMAT(fitAt_fitA_svd_V, n, n);
      GAVEC(fitAt_fitA_svd_S, n);
      GAVEC(tmpv, n);

      // Copy to fitA the fitting columns of the expanded super-sky metric. The order is always:
      //   f1dot, ... f(n-1)dot, freq
      gsl_vector_view freq = gsl_matrix_column(essky_dnorm, 5 + fsize - 1);
      gsl_vector_view fitA_freq = gsl_matrix_column(fitA, n - 1);
      gsl_vector_memcpy(&fitA_freq.vector, &freq.vector);
      if (n > 1) {
        gsl_matrix_view nspins = gsl_matrix_submatrix(essky_dnorm, 0, 5, 5 + fsize, n - 1);
        gsl_matrix_view fitA_nspins = gsl_matrix_submatrix(fitA, 0, 0, 5 + fsize, n - 1);
        gsl_matrix_memcpy(&fitA_nspins.matrix, &nspins.matrix);
      }

      // Compute fitA^T * fitA
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, fitA, fitA, 0.0, fitAt_fitA);

      // Find the singular value decomposition of fitA^T * fitA
      gsl_matrix_memcpy(fitAt_fitA_svd_AU, fitAt_fitA);
      GCALL(gsl_linalg_SV_decomp(fitAt_fitA_svd_AU, fitAt_fitA_svd_V, fitAt_fitA_svd_S, tmpv));

      // Determine if fitA^T * fitA is rank deficient, by looking at whether singular values
      // are below a given tolerance - taken from the Octave rank() function
      bool rank_deficient = false;
      const double tolerance = fsize * gsl_vector_get(fitAt_fitA_svd_S, 0) * DBL_EPSILON;
      for (size_t i = 0; i < fitAt_fitA_svd_S->size; ++i) {
        if (gsl_vector_get(fitAt_fitA_svd_S, i) <= tolerance) {
          rank_deficient = true;
        }
      }

      // If fitA^T * fitA is not rank deficient, we can use it for fitting
      if (!rank_deficient) {
        break;
      }

      // Otherwise cleanup for next loop, where we will reduce the rank of the fitting matrix
      GFMAT(fitA, fitAt_fitA, fitAt_fitA_svd_AU, fitAt_fitA_svd_V);
      GFVEC(fitAt_fitA_svd_S, tmpv);

    }

    // Get the rank of the fitting matrix, which is guaranteed to be at least 1
    const size_t n = fitA->size2;

    // Allocate memory
    gsl_vector* GAVEC(fitcoeff, n);

    // Construct subtraction matrix, which subtracts the least-squares linear fit by
    // frequency and spindowns from the orbital X and Y, using the singular value
    // decomposition of fitA^T * fitA to perform its inverse:
    //   fit_coeff = inv(fitAt_fitA) * fitA^T * fity);
    gsl_matrix_set_identity(subtractfit);
    for (size_t j = 2; j < 4; ++j) {

      // Construct a view of the column of the expanded super-sky metric to be fitted:
      // these are the X,Y orbital motion in ecliptic coordinates
      gsl_vector_view fity = gsl_matrix_column(essky_dnorm, j);

      // Compute the least-square fitting coefficients, using SVD for the inverse:
      //   fitcoeff = inv(fitA^T * fitA) * fitA^T * fity
      gsl_blas_dgemv(CblasTrans, 1.0, fitA, &fity.vector, 0.0, tmpv);
      GCALL(gsl_linalg_SV_solve(fitAt_fitA_svd_AU, fitAt_fitA_svd_V, fitAt_fitA_svd_S, tmpv, fitcoeff));

      // To subtract the fit, the fitting coefficients must be negated
      gsl_vector_scale(fitcoeff, -1.0);

      // Copy the fitting coefficents to the subtraction matrix
      gsl_matrix_set(subtractfit, 5 + fsize - 1, j, gsl_vector_get(fitcoeff, n - 1));
      if (n > 1) {
        gsl_matrix_view subtractfit_fitcoeff_nspins = gsl_matrix_submatrix(subtractfit, 5, j, n - 1, 1);
        gsl_matrix_view fitcoeff_nspins = gsl_matrix_view_vector(fitcoeff, n - 1, 1);
        gsl_matrix_memcpy(&subtractfit_fitcoeff_nspins.matrix, &fitcoeff_nspins.matrix);
      }

    }

    // Build matrix which constructs the intermediate "residual" super-sky metric:
    //   residual = essky_dnorm_transf * subtractfit * inv(essky_dnorm_transf) * reconstruct
    gsl_matrix_memcpy(tmp, reconstruct);
    gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, essky_dnorm_transf, tmp);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, subtractfit, tmp, 0.0, residual);
    gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, essky_dnorm_transf, residual);

    // Construct the intermediate "residual" super-sky metric from the fitted expanded super-sky metric:
    //    interm_ssky_metric = residual^T * essky_metric * residual
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, essky_metric, residual, 0.0, tmp);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, residual, tmp, 0.0, intm_ssky_metric);

    // Extract the sky offset vectors from the residual transform matrix,
    // and subtract them from the reduced super-sky coordinate transform data
    gsl_matrix_view residual_sky_offsets = gsl_matrix_submatrix(residual, 5, 0, fsize, 3);
    gsl_matrix_view sky_offsets = gsl_matrix_submatrix(*rssky_transf, 3, 0, fsize, 3);
    gsl_matrix_sub(&sky_offsets.matrix, &residual_sky_offsets.matrix);

    // Cleanup
    GFMAT(essky_dnorm, essky_dnorm_transf, fitA, fitAt_fitA, fitAt_fitA_svd_AU, fitAt_fitA_svd_V, reconstruct, residual, subtractfit, tmp);
    GFVEC(fitAt_fitA_svd_S, fitcoeff, tmpv);

    // Ensured intermediate metric is symmetric
    SSM_MakeSymmetric(intm_ssky_metric);

  }

  //
  // De-couple the sky-sky and frequency-frequency blocks
  // Produces the intermediate "decoupled" super-sky metric
  //

  {

    // Allocate memory
    gsl_matrix* GAMAT(decouple_sky_offsets, fsize, 3);
    gsl_matrix* GAMAT(freq_freq_dnorm, fsize, fsize);
    gsl_matrix* GAMAT(freq_freq_dnorm_transf, fsize, fsize);

    // Create views of the sky-sky, frequency-frequency, and off-diagonal blocks
    gsl_matrix_view sky_sky   = gsl_matrix_submatrix(intm_ssky_metric, 0, 0, 3, 3);
    gsl_matrix_view sky_freq  = gsl_matrix_submatrix(intm_ssky_metric, 0, 3, 3, fsize);
    gsl_matrix_view freq_sky  = gsl_matrix_submatrix(intm_ssky_metric, 3, 0, fsize, 3);
    gsl_matrix_view freq_freq = gsl_matrix_submatrix(intm_ssky_metric, 3, 3, fsize, fsize);

    // Diagonal-normalise the frequency-frequency block
    gsl_matrix_memcpy(freq_freq_dnorm, &freq_freq.matrix);
    SSM_DiagonalNormalise(freq_freq_dnorm, freq_freq_dnorm_transf);

    // Compute the additional sky offsets required to decouple the sky-sky and frequency blocks:
    //   decouple_sky_offsets = freq_freq_dnorm_transf * inv(freq_freq_dnorm) * freq_freq_dnorm_transf * freq_sky
    // Uses freq_sky as a temporary matrix, since it will be zeroed out anyway
    gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, freq_freq_dnorm_transf, &freq_sky.matrix);
    GCALL(gsl_linalg_cholesky_decomp(freq_freq_dnorm));
    GCALL(gsl_linalg_cholesky_invert(freq_freq_dnorm));
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, freq_freq_dnorm, &freq_sky.matrix, 0.0, decouple_sky_offsets);
    gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, freq_freq_dnorm_transf, decouple_sky_offsets);

    // Add the additional sky offsets to the reduced super-sky coordinate transform data
    gsl_matrix_view sky_offsets = gsl_matrix_submatrix(*rssky_transf, 3, 0, fsize, 3);
    gsl_matrix_add(&sky_offsets.matrix, decouple_sky_offsets);

    // Apply the decoupling transform to the sky-sky block:
    //   sky_sky = sky_sky - sky_freq * decoupl_sky_offsets
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &sky_freq.matrix, decouple_sky_offsets, 1.0, &sky_sky.matrix);

    // Zero out the off-diagonal blocks
    gsl_matrix_set_zero(&sky_freq.matrix);
    gsl_matrix_set_zero(&freq_sky.matrix);

    // Cleanup
    GFMAT(freq_freq_dnorm, freq_freq_dnorm_transf, decouple_sky_offsets);

    // Ensured intermediate metric is symmetric
    SSM_MakeSymmetric(intm_ssky_metric);

  }

  //
  // Align the sky-sky metric with its eigenvalues by means of a rotation
  // Produces the intermediate "aligned" super-sky metric
  //

  {

    // Allocate memory
    gsl_eigen_symmv_workspace* GALLOC(wksp, gsl_eigen_symmv_alloc(3));
    gsl_matrix* GAMAT(sky_evec, 3, 3);
    gsl_matrix* GAMAT(tmp, fsize, 3);
    gsl_permutation* GAPERM(luperm, 3);
    gsl_vector* GAVEC(sky_eval, 3);

    // Compute the eigenvalues/vectors of the sky-sky block
    gsl_matrix_view sky_sky = gsl_matrix_submatrix(intm_ssky_metric, 0, 0, 3, 3);
    GCALL(gsl_eigen_symmv(&sky_sky.matrix, sky_eval, sky_evec, wksp));

    // Sort the eigenvalues/vectors by descending absolute eigenvalue
    GCALL(gsl_eigen_symmv_sort(sky_eval, sky_evec, GSL_EIGEN_SORT_ABS_DESC));

    // Check that the first two eigenvalues are strictly positive
    for (size_t i = 0; i < 2; ++i) {
      XLAL_CHECK(gsl_vector_get(sky_eval, i) > 0.0, XLAL_ELOSS,
                 "Negative eigenvalue #%zu = %0.6e", i, gsl_vector_get(sky_eval, i));
    }

    // Set the sky-sky block to the diagonal matrix of eigenvalues
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

    // Store the alignment transform in the reduced super-sky coordinate transform data
    gsl_matrix_view align_sky = gsl_matrix_submatrix(*rssky_transf, 0, 0, 3, 3);
    gsl_matrix_transpose_memcpy(&align_sky.matrix, sky_evec);

    // Ensure that the alignment transform has a positive determinant,
    // to ensure that that it represents a rotation
    int lusign = 0;
    GCALL(gsl_linalg_LU_decomp(sky_evec, luperm, &lusign));
    if (gsl_linalg_LU_det(sky_evec, lusign) < 0.0) {
      gsl_vector_view col = gsl_matrix_column(&align_sky.matrix, 2);
      gsl_vector_scale(&col.vector, -1.0);
    }

    // Multiply the sky offsets by the alignment transform to transform to aligned sky coordinates:
    //   aligned_sky_off = sky_offsets * alignsky^T;
    gsl_matrix_view aligned_sky_offsets = gsl_matrix_submatrix(*rssky_transf, 3, 0, fsize, 3);
    gsl_matrix_memcpy(tmp, &aligned_sky_offsets.matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, tmp, &align_sky.matrix, 0.0, &aligned_sky_offsets.matrix);

    // Cleanup
    gsl_eigen_symmv_free(wksp);
    GFMAT(sky_evec, tmp);
    GFPERM(luperm);
    GFVEC(sky_eval);

    // Ensured intermediate metric is symmetric
    SSM_MakeSymmetric(intm_ssky_metric);

  }

  //
  // Drop the 3rd row/column of the metric, corresponding to the smallest sky eigenvalue
  // Produces the final reduced super-sky metric
  //

  {
    gsl_matrix_view aligned_sky2_sky2 = gsl_matrix_submatrix(intm_ssky_metric, 0, 0, 2, 2);
    gsl_matrix_view reduced_sky_sky = gsl_matrix_submatrix(*rssky_metric, 0, 0, 2, 2);
    gsl_matrix_memcpy(&reduced_sky_sky.matrix, &aligned_sky2_sky2.matrix);
    gsl_matrix_view aligned_freq_freq = gsl_matrix_submatrix(intm_ssky_metric, 3, 3, fsize, fsize);
    gsl_matrix_view reduced_freq_freq = gsl_matrix_submatrix(*rssky_metric, 2, 2, fsize, fsize);
    gsl_matrix_memcpy(&reduced_freq_freq.matrix, &aligned_freq_freq.matrix);
  }

  // Cleanup
  GFMAT(intm_ssky_metric);

  return XLAL_SUCCESS;

}
