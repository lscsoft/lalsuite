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

///
/// Convert from 2-dimensional reduced super-sky coordinates to 3-dimensional aligned sky coordinates
///
/// The 2-dimensional reduced super-sky coordinates \c rss = (\c A, \c B) encode the two hemispheres
/// of the sky as two neighbouring unit disks. The conversion to 3-dimensional aligned sky coordinates
/// is illustrated in the following diagram:
///
/** \verbatim
as[1] = B =___________________
        1_|   _____   _____   |
          |  /     \ /     \  |
        0-| |       |       | |
         _|  \_____/ \_____/  |
       -1 |_.___.___.___.___._|
            '   '   '   '   '
       A = -2  -1   0   1   2
   as[0] =  1   0  -1   0   1
   as[2] =  0  -1   0   1   0
\endverbatim */
///
/// Points outside the unit disks are moved radially onto their boundaries.
///
static void SSM_ReducedToAligned(
  double as[3],					///< [out] 3-dimensional aligned sky coordinates
  const gsl_vector* rss				///< [in] 2-dimensional reduced super-sky coordinates
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
/// Convert from 3-dimensional aligned sky coordinates to 2-dimensional reduced super-sky coordinates
///
static void SSM_AlignedToReduced(
  gsl_vector* rss,				///< [out] 2-dimensional reduced super-sky coordinates
  const double as[3]				///< [in] 3-dimensional aligned sky coordinates
  )
{
  const double r = sqrt(DOT3(as, as));
  const double A = GSL_SIGN(as[2]) * ((as[0] / r) + 1.0);
  const double B = as[1] / r;
  gsl_vector_set(rss, 0, A);
  gsl_vector_set(rss, 1, B);
}

int XLALConvertSuperSkyCoordinates(
  const SuperSkyCoordinates out,
  gsl_matrix** out_points,
  const SuperSkyCoordinates in,
  const gsl_matrix* in_points,
  const gsl_matrix* rssky_transf
  )
{

  // Check input
  XLAL_CHECK(out < SSC_MAX, XLAL_EINVAL);
  XLAL_CHECK(in < SSC_MAX, XLAL_EINVAL);
  XLAL_CHECK(out_points != NULL, XLAL_EFAULT);
  XLAL_CHECK(in_points != NULL, XLAL_EINVAL);
  XLAL_CHECK(rssky_transf != NULL || (out != SSC_REDUCED_SUPER_SKY && in != SSC_REDUCED_SUPER_SKY), XLAL_EINVAL);

  // Deduce number of input sky coordinates, and frequency/spindown coordinates
  const size_t in_ssize = (in == SSC_SUPER_SKY) ? 3 : 2;
  XLAL_CHECK(in_points->size1 > in_ssize, XLAL_EINVAL);
  const size_t fsize = in_points->size1 - in_ssize;

  // (Re)Allocate output points matrix
  const size_t out_ssize = (out == SSC_SUPER_SKY) ? 3 : 2;
  const size_t out_rows = fsize + out_ssize;
  if (*out_points != NULL && ((*out_points)->size1 != out_rows || (*out_points)->size2 != in_points->size2)) {
    gsl_matrix_free(*out_points);
    *out_points = NULL;
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
    SuperSkyCoordinates curr = in;

    // Convert physical coordinates to super-sky coordinates
    if (curr == SSC_PHYSICAL && out > curr) {

      // Convert right ascension and declination to super-sky position
      const double alpha = tmp[0];
      const double delta = tmp[1];
      const double cos_delta = cos(delta);
      tmp[0] = cos(alpha) * cos_delta;
      tmp[1] = sin(alpha) * cos_delta;
      tmp[2] = sin(delta);

      // Move frequency to last dimension
      const double freq = tmp[3];
      memmove(&tmp[3], &tmp[4], (fsize - 1) * sizeof(tmp[0]));
      tmp[2 + fsize] = freq;

      // Update current coordinate system
      curr = SSC_SUPER_SKY;

    }

    // Convert super-sky coordinates to reduced super-sky coordinates
    if (curr == SSC_SUPER_SKY && out > curr) {

      // Create views of the sky alignment transform and sky offset vectors
      gsl_matrix_const_view align_sky = gsl_matrix_const_submatrix(rssky_transf, 0, 0, 3, 3);
      gsl_matrix_const_view sky_offsets = gsl_matrix_const_submatrix(rssky_transf, 3, 0, fsize, 3);

      // Apply the alignment transform to the super-sky position to produced the aligned sky position:
      //   asky = align_sky * ssky
      double asky[3];
      gsl_vector_view asky_v = gsl_vector_view_array(asky, 3);
      gsl_blas_dgemv(CblasNoTrans, 1.0, &align_sky.matrix, &tmp_sky.vector, 0.0, &asky_v.vector);

      // Add the inner product of the sky offsets with the aligned sky position
      // to the super-sky spins and frequency to get the reduced super-sky quantities:
      //   rssky_fspin[i] = ssky_fspin[i] + dot(sky_offsets[i], asky)
      gsl_blas_dgemv(CblasNoTrans, 1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &tmp_fspin.vector);

      // Convert from 3-dimensional aligned sky coordinates to 2-dimensional reduced super-sky coordinates
      SSM_AlignedToReduced(&tmp_sky.vector, asky);

      // Update current coordinate system
      curr = SSC_REDUCED_SUPER_SKY;

    }

    // Convert reduced super-sky coordinates to super-sky coordinates
    if (curr == SSC_REDUCED_SUPER_SKY && out < curr) {

      // Create views of the sky alignment transform and sky offset vectors
      gsl_matrix_const_view align_sky = gsl_matrix_const_submatrix(rssky_transf, 0, 0, 3, 3);
      gsl_matrix_const_view sky_offsets = gsl_matrix_const_submatrix(rssky_transf, 3, 0, fsize, 3);

      // Convert from 2-dimensional reduced super-sky coordinates to 3-dimensional aligned sky coordinates
      double asky[3];
      SSM_ReducedToAligned(asky, &tmp_sky.vector);
      gsl_vector_view asky_v = gsl_vector_view_array(asky, 3);

      // Subtract the inner product of the sky offsets with the aligned sky position
      // from the reduced super-sky spins and frequency to get the super-sky quantities:
      //   ssky_fspin[i] = rssky_fspin[i] - dot(sky_offsets[i], asky)
      gsl_blas_dgemv(CblasNoTrans, -1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &tmp_fspin.vector);

      // Apply the inverse alignment transform to the aligned sky position to produced the super-sky position:
      //   ssky = align_sky^T * asky
      gsl_blas_dgemv(CblasTrans, 1.0, &align_sky.matrix, &asky_v.vector, 0.0, &tmp_sky.vector);

      // Update current coordinate system
      curr = SSC_SUPER_SKY;

    }

    // Convert super-sky coordinates to physical coordinates
    if (curr == SSC_SUPER_SKY && out < curr) {

      // Convert super-sky position to right ascension and declination
      const double nx = tmp[0];
      const double ny = tmp[1];
      const double nz = tmp[2];
      tmp[0] = atan2(ny, nx);
      tmp[1] = atan2(nz, sqrt(SQR(nx) + SQR(ny)));
      XLALNormalizeSkyPosition(&tmp[0], &tmp[1]);

      // Move frequency to first dimension
      const double freq = tmp[2 + fsize];
      memmove(&tmp[4], &tmp[3], (fsize - 1) * sizeof(tmp[0]));
      tmp[3] = freq;

      // Update current coordinate system
      curr = SSC_PHYSICAL;

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

int XLALConvertPhysicalToSuperSky(
  const SuperSkyCoordinates out,
  gsl_vector* out_point,
  const PulsarDopplerParams* in_phys,
  const gsl_matrix* rssky_transf
  )
{

  // Check input
  XLAL_CHECK(SSC_PHYSICAL < out && out < SSC_MAX, XLAL_EINVAL);
  XLAL_CHECK(out_point != NULL, XLAL_EFAULT);
  XLAL_CHECK(in_phys != NULL, XLAL_EFAULT);

  // Deduce number of sky coordinates, and frequency/spindown coordinates
  const size_t ssize = (out == SSC_SUPER_SKY) ? 3 : 2;
  XLAL_CHECK(out_point->size > ssize, XLAL_EINVAL);
  const size_t fsize = out_point->size - ssize;
  XLAL_CHECK(fsize <= PULSAR_MAX_SPINS, XLAL_EFAILED);

  // Copy input physical point to array
  double in_point[2 + fsize];
  in_point[0] = in_phys->Alpha;
  in_point[1] = in_phys->Delta;
  memcpy(&in_point[2], in_phys->fkdot, fsize * sizeof(in_point[0]));

  // Convert input physical point to output super-sky coordinate point
  gsl_matrix_view out_point_view = gsl_matrix_view_vector(out_point, out_point->size, 1);
  gsl_matrix_const_view in_point_view = gsl_matrix_const_view_array(in_point, 2 + fsize, 1);
  gsl_matrix* out_point_view_ptr = &out_point_view.matrix;
  XLAL_CHECK(XLALConvertSuperSkyCoordinates(out, &out_point_view_ptr, SSC_PHYSICAL, &in_point_view.matrix, rssky_transf) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK(out_point_view_ptr == &out_point_view.matrix, XLAL_EFAILED);

  return XLAL_SUCCESS;

}

int XLALConvertSuperSkyToPhysical(
  PulsarDopplerParams* out_phys,
  const SuperSkyCoordinates in,
  const gsl_vector* in_point,
  const gsl_matrix* rssky_transf
  )
{

  // Check input
  XLAL_CHECK(out_phys != NULL, XLAL_EFAULT);
  XLAL_CHECK(SSC_PHYSICAL < in && in < SSC_MAX, XLAL_EINVAL);
  XLAL_CHECK(in_point != NULL, XLAL_EFAULT);

  // Deduce number of sky coordinates, and frequency/spindown coordinates
  const size_t ssize = (in == SSC_SUPER_SKY) ? 3 : 2;
  XLAL_CHECK(in_point->size > ssize, XLAL_EINVAL);
  const size_t fsize = in_point->size - ssize;
  XLAL_CHECK(fsize <= PULSAR_MAX_SPINS, XLAL_EFAILED);

  // Convert input super-sky coordinate point to output physical point
  double out_point[2 + fsize];
  gsl_matrix_view out_point_view = gsl_matrix_view_array(out_point, 2 + fsize, 1);
  gsl_matrix_const_view in_point_view = gsl_matrix_const_view_vector(in_point, in_point->size, 1);
  gsl_matrix* out_point_view_ptr = &out_point_view.matrix;
  XLAL_CHECK(XLALConvertSuperSkyCoordinates(SSC_PHYSICAL, &out_point_view_ptr, in, &in_point_view.matrix, rssky_transf) == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK(out_point_view_ptr == &out_point_view.matrix, XLAL_EFAILED);

  // Copy output physical point from array
  out_phys->Alpha = out_point[0];
  out_phys->Delta = out_point[1];
  memcpy(out_phys->fkdot, &out_point[2], fsize * sizeof(out_point[0]));

  return XLAL_SUCCESS;

}

static void ReducedSuperSkyBCoordBound(
  const size_t dimension UNUSED,
  const gsl_vector* point,
  const gsl_vector* bbox,
  const void* data,
  double* bound,
  double* padding
  )
{

  // Get bounds data
  const double sgn = ((const double*) data)[0];

  // Get 2-dimensional reduced super-sky A coordinate
  const double A = gsl_vector_get(point, 0);
  const double dA = fabs(A) - 1.0;

  // Set bound on 2-dimensional reduced super-sky B coordinate
  const double maxB = RE_SQRT(1.0 - SQR(dA));
  *bound = sgn * maxB;

  // Add sufficient extra padding on B, such that the bounding box of the
  // boundary templates will not intersect the circular A-B parameter spaces
  const double hbboxA = 0.5 * gsl_vector_get(bbox, 0);
  const double absdA = fabs(dA);
  if (absdA <= hbboxA) {
    *padding += 1.0 - maxB;
  } else if (absdA <= 1.0 + hbboxA) {
    const double ddA = (dA < 0.0) ? dA + hbboxA : dA - hbboxA;
    *padding += RE_SQRT(1.0 - ddA * ddA) - maxB;
  }

}

int XLALSetLatticeReducedSuperSkyBounds(
  LatticeTiling* tiling
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);

  // Set the parameter-space bound on 2-dimensional reduced super-sky A coordinate
  XLAL_CHECK(XLALSetLatticeConstantBound(tiling, 0, -2.0, 2.0) == XLAL_SUCCESS, XLAL_EFUNC);

  // Allocate memory
  const size_t data_len = sizeof(double);
  double* data_lower = XLALMalloc(data_len);
  XLAL_CHECK(data_lower != NULL, XLAL_ENOMEM);
  double* data_upper = XLALMalloc(data_len);
  XLAL_CHECK(data_upper != NULL, XLAL_ENOMEM);

  // Set the parameter-space bound on 2-dimensional reduced super-sky B coordinate
  data_lower[0] = -1.0;
  data_upper[0] = +1.0;
  XLAL_CHECK(XLALSetLatticeBound(tiling, 1, ReducedSuperSkyBCoordBound, data_len, data_lower, data_upper) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

int XLALSetLatticeReducedSuperSkyPointBounds(
  LatticeTiling* tiling,
  const gsl_matrix* rssky_transf,
  const double alpha,
  const double delta
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);

  // Allocate memory
  gsl_vector* GAVEC(rssky_point, rssky_transf->size1 - 1);

  // Convert right ascension and declination to reduced super-sky coordinates
  PulsarDopplerParams XLAL_INIT_DECL(doppler);
  doppler.Alpha = alpha;
  doppler.Delta = delta;
  XLAL_CHECK(XLALConvertPhysicalToSuperSky(SSC_REDUCED_SUPER_SKY, rssky_point, &doppler, rssky_transf) == XLAL_SUCCESS, XLAL_EFUNC);

  // Set the parameter-space bounds on 2-dimensional reduced super-sky A and B coordinates
  for (size_t i = 0; i < 2; ++i) {
    const double rssky_point_i = gsl_vector_get(rssky_point, i);
    XLAL_CHECK(XLALSetLatticeConstantBound(tiling, i, rssky_point_i, rssky_point_i) == XLAL_SUCCESS, XLAL_EFUNC);
  }

  // Cleanup
  GFVEC(rssky_point);

  return XLAL_SUCCESS;

}

static void SuperSkySpinBound(
  const size_t dimension UNUSED,
  const gsl_vector* point,
  const gsl_vector* bbox UNUSED,
  const void* data,
  double* bound,
  double* padding UNUSED
  )
{

  // Get bounds data
  const double* sky_offsets = ((const double*) data);
  const double padding_scale = ((const double*) data)[3];
  *bound = ((const double*) data)[4];

  // Add the inner product of the sky offsets with the aligned sky
  // position to the super-sky bound to get the reduced super-sky bound
  double as[3];
  SSM_ReducedToAligned(as, point);
  *bound += DOT3(sky_offsets, as);

  // Add sufficient extra padding to the reduced super-sky bound to cover any
  // mismatch in sky position, within the bounding box around each sky point
  const double bboxA = gsl_vector_get(bbox, 0);
  const double bboxB = gsl_vector_get(bbox, 1);
  *padding += padding_scale * (fabs(sky_offsets[0]) * bboxA + fabs(sky_offsets[1]) * bboxB);

}

int XLALSetLatticeSuperSkySpinBound(
  LatticeTiling* tiling,
  const gsl_matrix* rssky_transf,
  const size_t s,
  const double bound1,
  const double bound2,
  bool extra_padding
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf->size1 > 3, XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound1), XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound2), XLAL_EINVAL);
  const size_t smax = rssky_transf->size1 - 4;
  XLAL_CHECK(s <= smax, XLAL_ESIZE);
  const size_t dimension = (s == 0) ? (2 + smax) : (1 + s);
  XLAL_CHECK(dimension < XLALLatticeDimensions(tiling), XLAL_ESIZE);

  // Allocate memory
  const size_t data_len = 5 * sizeof(double);
  double* data_lower = XLALMalloc(data_len);
  XLAL_CHECK(data_lower != NULL, XLAL_ENOMEM);
  double* data_upper = XLALMalloc(data_len);
  XLAL_CHECK(data_upper != NULL, XLAL_ENOMEM);

  // Copy the sky offset vector to bounds data
  for (size_t j = 0; j < 3; ++j) {
    data_lower[j] = data_upper[j] = gsl_matrix_get(rssky_transf, dimension + 1, j);
  }

  // Indicate whether to use extra padding to cover sky position mismatch
  data_lower[3] = data_upper[3] = (extra_padding ? 1.0 : 0.0);

  // Set the parameter-space bound on (untransformed) super-sky frequency/spindown coordinate
  data_lower[4] = GSL_MIN(bound1, bound2);
  data_upper[4] = GSL_MAX(bound1, bound2);
  XLAL_CHECK(XLALSetLatticeBound(tiling, dimension, SuperSkySpinBound, data_len, data_lower, data_upper) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}

int XLALSetLatticeReducedSuperSkySpinBound(
  LatticeTiling* tiling,
  const gsl_matrix* rssky_transf,
  const size_t s,
  const double bound1,
  const double bound2
  )
{

  // Check input
  XLAL_CHECK(tiling != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf != NULL, XLAL_EFAULT);
  XLAL_CHECK(rssky_transf->size1 > 3, XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound1), XLAL_EINVAL);
  XLAL_CHECK(isfinite(bound2), XLAL_EINVAL);
  const size_t smax = rssky_transf->size1 - 4;
  XLAL_CHECK(s <= smax, XLAL_ESIZE);
  const size_t dimension = (s == 0) ? (2 + smax) : (1 + s);
  XLAL_CHECK(dimension < XLALLatticeDimensions(tiling), XLAL_ESIZE);

  // Set the parameter-space bound on reduced super-sky frequency/spindown coordinate
  XLAL_CHECK(XLALSetLatticeConstantBound(tiling, dimension, bound1, bound2) == XLAL_SUCCESS, XLAL_EFUNC);

  return XLAL_SUCCESS;

}
