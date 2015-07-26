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
  XLAL_CHECK_NULL( coords != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( segments != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( XLALSegListIsInitialized( segments ), XLAL_EINVAL );
  XLAL_CHECK_NULL( segments->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( fiducial_freq > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( detectors != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( detectors->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( detector_motion > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( ephemerides != NULL, XLAL_EINVAL );

  // Create parameters struct for XLALComputeDopplerPhaseMetric()
  DopplerMetricParams XLAL_INIT_DECL( par );

  // Set coordinate system
  par.coordSys = *coords;

  // Set detector motion type
  par.detMotionType = detector_motion;

  // Set segment list
  par.segmentList = *segments;

  // Set detectors and detector weights
  par.multiIFO = *detectors;
  if( detector_weights != NULL ) {
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
  DopplerPhaseMetric *metric = XLALComputeDopplerPhaseMetric( &par, ephemerides );
  XLAL_CHECK_NULL( metric != NULL && metric->g_ij != NULL, XLAL_EFUNC, "XLALComputeDopplerPhaseMetric() failed" );

  // Extract metric
  gsl_matrix *g_ij = metric->g_ij;
  metric->g_ij = NULL;

  // Cleanup
  XLALDestroyDopplerPhaseMetric( metric );

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
  XLAL_CHECK( fitted_ssky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf != NULL, XLAL_EFAULT );
  XLAL_CHECK( ussky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( orbital_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );

  // Size of the frequency+spindowns block
  const size_t fsize = 1 + spindowns;

  // Allocate memory
  gsl_matrix *GAMAT( tmp, 2 + fsize, 2 + fsize );
  gsl_vector *GAVEC( tmpv, fsize );

  // Compute mid-time of segment list
  LIGOTimeGPS mid_time;
  {
    const LIGOTimeGPS *start_time = &(segments->segs[0].start);
    const LIGOTimeGPS *end_time   = &(segments->segs[segments->length - 1].end);
    const REAL8 time_span = XLALGPSDiff( end_time, start_time );
    mid_time = *start_time;
    XLALGPSAdd( &mid_time, 0.5 * time_span );
  }

  // Internal copy of orbital metric, and various transforms performed on it
  gsl_matrix *orb_metric = NULL, *mid_time_transf = NULL, *diag_norm_transf = NULL;

  // Transform reference time of orbital metric from reference time to segment list mid-time
  const REAL8 Dtau = XLALGPSDiff( &mid_time, ref_time );
  XLAL_CHECK( XLALChangeMetricReferenceTime( &orb_metric, &mid_time_transf, orbital_metric, ocoords, Dtau ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Diagonally-normalize orbital metric
  XLAL_CHECK( XLALDiagNormalizeMetric( &orb_metric, &diag_norm_transf, orb_metric ) == XLAL_SUCCESS, XLAL_EFUNC );

  // 'fitA' contains the frequency and spindown elements of the orbital metric, used for fitting
  gsl_matrix *GAMAT( fitA, 2 + fsize, fsize );
  {
    gsl_matrix_view orb_metric_fspin = gsl_matrix_submatrix( orb_metric, 0, 2, 2 + fsize, fsize );
    gsl_matrix_memcpy( fitA, &orb_metric_fspin.matrix );
  }

  // Compute 'fitA^T * fitA'
  gsl_matrix *GAMAT( fitAt_fitA, fsize, fsize );
  gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, fitA, fitA, 0.0, fitAt_fitA );

  // Find the singular value decomposition of 'fitA^T * fitA'
  gsl_matrix *GAMAT( svd_U, fsize, fsize );
  gsl_matrix *GAMAT( svd_V, fsize, fsize );
  gsl_vector *GAVEC( svd_S, fsize );
  gsl_matrix_memcpy( svd_U, fitAt_fitA );
  GCALL( gsl_linalg_SV_decomp( svd_U, svd_V, svd_S, tmpv ) );

  // The columns of 'fitc' contain the least-square fitting coefficients for the orbital X and Y metric elements:
  //    fitc(:,j) = inv(fitA^T * fitA) * fitA^T * orb_metric(:,j)
  // The singular decomposition of fitA^T * fitA is used for the inverse
  gsl_matrix *GAMAT( fitc, fsize, 2 );
  for( size_t j = 0; j < 2; ++j ) {
    gsl_vector_view orb_metric_j = gsl_matrix_column( orb_metric, j );
    gsl_vector_view fitc_j = gsl_matrix_column( fitc, j );
    gsl_blas_dgemv( CblasTrans, 1.0, fitA, &orb_metric_j.vector, 0.0, tmpv );
    GCALL( gsl_linalg_SV_solve( svd_U, svd_V, svd_S, tmpv, &fitc_j.vector ) );
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
  gsl_matrix *GAMAT( subtract_orb, 2 + fsize, 2 + fsize );
  {
    gsl_matrix_set_identity( subtract_orb );
    gsl_matrix_view subtract_orb_fspin_sky = gsl_matrix_submatrix( subtract_orb, 2, 0, fsize, 2 );
    gsl_matrix_memcpy( &subtract_orb_fspin_sky.matrix, fitc );
    gsl_matrix_scale( &subtract_orb_fspin_sky.matrix, -1.0 );
  }

  // Multiply 'subtract_orb' by the diagonal-normalization and reference time transforms,
  // to obtain the matrix that substracts the fit from the unconstrained supersky metric
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, diag_norm_transf, subtract_orb, 0.0, tmp );
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, mid_time_transf, tmp, 0.0, subtract_orb );

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
  gsl_matrix *GAMAT( subtract_ussky, 3 + fsize, 3 + fsize );
  {
    gsl_matrix_set_identity( subtract_ussky );
    gsl_matrix_view subtract_ussky_fspin_sky = gsl_matrix_submatrix( subtract_ussky, 3, 0, fsize, 3 );
    gsl_matrix_view subtract_orb_fspin_sky = gsl_matrix_submatrix( subtract_orb, 2, 0, fsize, 2 );
    {
      gsl_vector_view subtract_ussky_fspin_sky_col = gsl_matrix_column( &subtract_ussky_fspin_sky.matrix, 0 );
      gsl_vector_view subtract_orb_fspin_sky_col = gsl_matrix_column( &subtract_orb_fspin_sky.matrix, 0 );
      gsl_vector_memcpy( &subtract_ussky_fspin_sky_col.vector, &subtract_orb_fspin_sky_col.vector );
    }
    {
      gsl_vector_view subtract_ussky_fspin_sky_col = gsl_matrix_column( &subtract_ussky_fspin_sky.matrix, 1 );
      gsl_vector_view subtract_orb_fspin_sky_col = gsl_matrix_column( &subtract_orb_fspin_sky.matrix, 1 );
      gsl_vector_memcpy( &subtract_ussky_fspin_sky_col.vector, &subtract_orb_fspin_sky_col.vector );
      gsl_vector_scale( &subtract_ussky_fspin_sky_col.vector, LAL_COSIEARTH );
    }
    {
      gsl_vector_view subtract_ussky_fspin_sky_col = gsl_matrix_column( &subtract_ussky_fspin_sky.matrix, 2 );
      gsl_vector_view subtract_orb_fspin_sky_col = gsl_matrix_column( &subtract_orb_fspin_sky.matrix, 1 );
      gsl_vector_memcpy( &subtract_ussky_fspin_sky_col.vector, &subtract_orb_fspin_sky_col.vector );
      gsl_vector_scale( &subtract_ussky_fspin_sky_col.vector, LAL_SINIEARTH );
    }
  }

  // Transform the unconstrained supersky metric to the intermediate fitted supersky metric
  XLAL_CHECK( XLALTransformMetric( &fitted_ssky_metric, subtract_ussky, ussky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Extract the sky offset vectors from 'subtract_ussky', and subtract them from the reduced supersky coordinate transform data
  {
    gsl_matrix_view subtract_ussky_fspin_sky = gsl_matrix_submatrix( subtract_ussky, 3, 0, fsize, 3 );
    gsl_matrix_view sky_offsets = gsl_matrix_submatrix( rssky_transf, 3, 0, fsize, 3 );
    gsl_matrix_sub( &sky_offsets.matrix, &subtract_ussky_fspin_sky.matrix );
  }

  // Cleanup
  GFMAT( diag_norm_transf, fitA, fitAt_fitA, fitc, mid_time_transf, orb_metric, subtract_orb, subtract_ussky, svd_U, svd_V, tmp );
  GFVEC( svd_S, tmpv );

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
  XLAL_CHECK( decoupled_ssky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf != NULL, XLAL_EFAULT );
  XLAL_CHECK( fitted_ssky_metric != NULL, XLAL_EFAULT );

  // Size of the frequency+spindowns block
  const size_t fsize = 1 + spindowns;

  // Copy fitted metric to decoupled metric
  gsl_matrix_memcpy(decoupled_ssky_metric, fitted_ssky_metric);

  // Create views of the sky--sky, freq+spin--freq+spin, and off-diagonal blocks
  gsl_matrix_view sky_sky     = gsl_matrix_submatrix( decoupled_ssky_metric, 0, 0, 3, 3 );
  gsl_matrix_view sky_fspin   = gsl_matrix_submatrix( decoupled_ssky_metric, 0, 3, 3, fsize );
  gsl_matrix_view fspin_sky   = gsl_matrix_submatrix( decoupled_ssky_metric, 3, 0, fsize, 3 );
  gsl_matrix_view fspin_fspin = gsl_matrix_submatrix( decoupled_ssky_metric, 3, 3, fsize, fsize );

  // Diagonal-normalise the freq+spin--freq+spin block
  gsl_matrix *fspin_fspin_dnorm = NULL, *fspin_fspin_dnorm_transf = NULL;
  XLAL_CHECK( XLALDiagNormalizeMetric( &fspin_fspin_dnorm, &fspin_fspin_dnorm_transf, &fspin_fspin.matrix ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Invert the freq+spin--freq+spin block
  gsl_matrix *GAMAT( fspin_fspin_dnorm_LU, fsize, fsize );
  gsl_matrix *GAMAT( fspin_fspin_dnorm_inv, fsize, fsize );
  gsl_permutation *GAPERM( fspin_fspin_dnorm_LU_perm, fsize );
  int fspin_fspin_dnorm_LU_sign = 0;
  gsl_matrix_memcpy( fspin_fspin_dnorm_LU, fspin_fspin_dnorm );
  GCALL( gsl_linalg_LU_decomp( fspin_fspin_dnorm_LU, fspin_fspin_dnorm_LU_perm, &fspin_fspin_dnorm_LU_sign ) );
  GCALL( gsl_linalg_LU_invert( fspin_fspin_dnorm_LU, fspin_fspin_dnorm_LU_perm, fspin_fspin_dnorm_inv ) );

  // Compute the additional sky offsets required to decouple the sky--sky and frequency blocks:
  //   decouple_sky_offsets = fspin_fspin_dnorm_transf * inv(fspin_fspin_dnorm) * fspin_fspin_dnorm_transf * fspin_sky
  // Uses fspin_sky as a temporary matrix, since it will be zeroed out anyway
  gsl_matrix *GAMAT( decouple_sky_offsets, fsize, 3 );
  gsl_blas_dtrmm( CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, fspin_fspin_dnorm_transf, &fspin_sky.matrix );
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, fspin_fspin_dnorm_inv, &fspin_sky.matrix, 0.0, decouple_sky_offsets );
  gsl_blas_dtrmm( CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, fspin_fspin_dnorm_transf, decouple_sky_offsets );

  // Add the additional sky offsets to the reduced supersky coordinate transform data
  gsl_matrix_view sky_offsets = gsl_matrix_submatrix( rssky_transf, 3, 0, fsize, 3 );
  gsl_matrix_add( &sky_offsets.matrix, decouple_sky_offsets );

  // Apply the decoupling transform to the sky--sky block:
  //   sky_sky = sky_sky - sky_fspin * decouplp_sky_offsets
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, -1.0, &sky_fspin.matrix, decouple_sky_offsets, 1.0, &sky_sky.matrix );

  // Zero out the off-diagonal blocks
  gsl_matrix_set_zero( &sky_fspin.matrix );
  gsl_matrix_set_zero( &fspin_sky.matrix );

  // Cleanup
  GFPERM( fspin_fspin_dnorm_LU_perm );
  GFMAT( fspin_fspin_dnorm, fspin_fspin_dnorm_LU, fspin_fspin_dnorm_inv, fspin_fspin_dnorm_transf, decouple_sky_offsets );

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
  XLAL_CHECK( aligned_ssky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf != NULL, XLAL_EFAULT );
  XLAL_CHECK( decoupled_ssky_metric != NULL, XLAL_EFAULT );

  // Size of the frequency+spindowns block
  const size_t fsize = 1 + spindowns;

  // Allocate memory
  gsl_matrix *GAMAT( tmp, fsize, 3 );

  // Copy decoupled metric to aligned metric
  gsl_matrix_memcpy(aligned_ssky_metric, decoupled_ssky_metric);

  // Compute the eigenvalues/vectors of the sky--sky block
  gsl_vector *GAVEC( sky_eval, 3 );
  gsl_matrix *GAMAT( sky_evec, 3, 3 );
  gsl_eigen_symmv_workspace *GALLOC( wksp, gsl_eigen_symmv_alloc( 3 ) );
  gsl_matrix_view sky_sky = gsl_matrix_submatrix( aligned_ssky_metric, 0, 0, 3, 3 );
  GCALL( gsl_eigen_symmv( &sky_sky.matrix, sky_eval, sky_evec, wksp ) );

  // Sort the eigenvalues/vectors by descending absolute eigenvalue
  GCALL( gsl_eigen_symmv_sort( sky_eval, sky_evec, GSL_EIGEN_SORT_ABS_DESC ) );

  // Set the sky--sky block to the diagonal matrix of eigenvalues
  gsl_matrix_set_zero( &sky_sky.matrix );
  gsl_vector_view sky_sky_diag = gsl_matrix_diagonal( &sky_sky.matrix );
  gsl_vector_memcpy( &sky_sky_diag.vector, sky_eval );

  // Ensure that the matrix of eigenvalues has a positive diagonal; this and
  // the determinant constraints ensures fully constraints the eigenvector signs
  for( size_t j = 0; j < 3; ++j ) {
    gsl_vector_view col = gsl_matrix_column( sky_evec, j );
    if( gsl_vector_get( &col.vector, j ) < 0.0 ) {
      gsl_vector_scale( &col.vector, -1.0 );
    }
  }

  // Store the alignment transform in the reduced supersky coordinate transform data
  gsl_matrix_view align_sky = gsl_matrix_submatrix( rssky_transf, 0, 0, 3, 3 );
  gsl_matrix_transpose_memcpy( &align_sky.matrix, sky_evec );

  // Ensure that the alignment transform has a positive determinant,
  // to ensure that that it represents a rotation
  gsl_permutation *GAPERM( LU_perm, 3 );
  int LU_sign = 0;
  GCALL( gsl_linalg_LU_decomp( sky_evec, LU_perm, &LU_sign ) );
  if( gsl_linalg_LU_det( sky_evec, LU_sign ) < 0.0 ) {
    gsl_vector_view col = gsl_matrix_column( &align_sky.matrix, 2 );
    gsl_vector_scale( &col.vector, -1.0 );
  }

  // Multiply the sky offsets by the alignment transform to transform to aligned sky coordinates:
  //   aligned_sky_off = sky_offsets * alignsky^T;
  gsl_matrix_view aligned_sky_offsets = gsl_matrix_submatrix( rssky_transf, 3, 0, fsize, 3 );
  gsl_matrix_memcpy( tmp, &aligned_sky_offsets.matrix );
  gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, tmp, &align_sky.matrix, 0.0, &aligned_sky_offsets.matrix );

  // Cleanup
  gsl_eigen_symmv_free( wksp );
  GFMAT( sky_evec, tmp );
  GFPERM( LU_perm );
  GFVEC( sky_eval );

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
  XLAL_CHECK( rssky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf != NULL, XLAL_EFAULT );
  XLAL_CHECK( aligned_ssky_metric != NULL, XLAL_EFAULT );
  const size_t n = aligned_ssky_metric->size1;
  const size_t m = rssky_metric->size1;

  // Internal copy of aligned supersky metric
  gsl_matrix *GAMAT( aln_metric, n, n );
  gsl_matrix_memcpy( aln_metric, aligned_ssky_metric );

  // Move the 3rd row/column of 'aln_metric', which is the 'n_c' sky coordinate with
  // the smallest eigenvalue, to the last row/column, so it can be easily dropped
  for( size_t i = 2; i + 1 < n; ++i ) {
    gsl_matrix_swap_rows( aln_metric, i, i + 1 );
    gsl_matrix_swap_columns( aln_metric, i, i + 1 );
  }

  // Move the 3rd row/column of 'aln_metric', which is *now* the frequency,
  // to the second-to-last row/column, i.e. still before 'n_c'
  for( size_t i = 2; i + 2 < n; ++i ) {
    gsl_matrix_swap_rows( aln_metric, i, i + 1 );
    gsl_matrix_swap_columns( aln_metric, i, i + 1 );
  }

  // Copy the first 'm' dimensions of 'aln_metric' to 'rssky_metric', dropping 'n_c'
  {
    gsl_matrix_view aln_metric_nm1_nm1 = gsl_matrix_submatrix( aln_metric, 0, 0, m, m );
    gsl_matrix_memcpy( rssky_metric, &aln_metric_nm1_nm1.matrix );
  }

  // Move the 4th row of 'rssky_transf', which is the coordinate
  // transform data for frequency, to the last row
  for( size_t i = 3; i + 1 < n; ++i ) {
    gsl_matrix_swap_rows( rssky_transf, i, i + 1 );
  }

  // Ensure reduced supersky metric is symmetric
  for( size_t i = 0; i < m; ++i ) {
    for( size_t j = i + 1; j < m; ++j ) {
      const double gij = gsl_matrix_get( rssky_metric, i, j );
      const double gji = gsl_matrix_get( rssky_metric, j, i );
      const double g = 0.5 * ( gij + gji );
      gsl_matrix_set( rssky_metric, i, j, g );
      gsl_matrix_set( rssky_metric, j, i, g );
    }
  }

  // Ensure reduced supersky metric is positive definite
  for ( size_t s = 1; s <= rssky_metric->size1; ++s ) {
    gsl_matrix_view rssky_metric_s = gsl_matrix_submatrix( rssky_metric, 0, 0, s, s );
    const double det_s = XLALMetricDeterminant( &rssky_metric_s.matrix );
    XLAL_CHECK( det_s > 0, XLAL_EFAILED, "Reduced supersky metric is not positive definite (s=%zu, det_s=%0.3e)", s, det_s );
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
  XLAL_CHECK( p_rssky_metric == NULL || *p_rssky_metric == NULL, XLAL_EINVAL, "'*p_rssky_metric' must be NULL of 'p_rssky_metric' is non-NULL" );
  XLAL_CHECK( p_rssky_transf == NULL || *p_rssky_transf == NULL, XLAL_EINVAL, "'*p_rssky_transf' must be NULL of 'p_rssky_transf' is non-NULL" );
  XLAL_CHECK( p_ussky_metric == NULL || *p_ussky_metric == NULL, XLAL_EINVAL, "'*p_ussky_metric' must be NULL of 'p_ussky_metric' is non-NULL" );
  XLAL_CHECK( (p_rssky_metric != NULL) == (p_rssky_transf != NULL), XLAL_EINVAL, "Both 'p_rssky_metric' and 'p_rssky_transf' must be either NULL or non-NULL" );
  XLAL_CHECK( (p_rssky_metric != NULL) || (p_ussky_metric != NULL), XLAL_EINVAL, "At least one of 'p_rssky_metric' or 'p_ussky_metric' must be non-NULL" );
  XLAL_CHECK( spindowns <= 3, XLAL_EINVAL );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK( segments != NULL, XLAL_EFAULT );
  XLAL_CHECK( XLALSegListIsInitialized( segments ), XLAL_EINVAL );
  XLAL_CHECK( segments->length > 0, XLAL_EINVAL );
  XLAL_CHECK( fiducial_freq > 0, XLAL_EINVAL );
  XLAL_CHECK( detectors != NULL, XLAL_EFAULT );
  XLAL_CHECK( detectors->length > 0, XLAL_EINVAL );
  XLAL_CHECK( detector_motion > 0, XLAL_EINVAL );
  XLAL_CHECK( ephemerides != NULL, XLAL_EINVAL );

  // Size of the frequency+spindowns block
  const size_t fsize = 1 + spindowns;

  // Build coordinate system for the unrestricted supersky metric and orbital metric
  DopplerCoordinateSystem XLAL_INIT_DECL( ucoords );
  DopplerCoordinateSystem XLAL_INIT_DECL( ocoords );
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
  if( spindowns >= 1 ) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F1DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F1DOT;
  }
  if( spindowns >= 2 ) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F2DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F2DOT;
  }
  if( spindowns >= 3 ) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F3DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F3DOT;
  }

  // Compute the unrestricted supersky metric
  gsl_matrix *ussky_metric = SM_ComputePhaseMetric( &ucoords, ref_time, segments, fiducial_freq, detectors, detector_weights, detector_motion, ephemerides );
  XLAL_CHECK( ussky_metric != NULL, XLAL_EFUNC );

  // Compute the reduced supersky metric and coordinate transform data
  if( p_rssky_metric != NULL ) {

    // Allocate memory
    GAMAT( *p_rssky_metric, 2 + fsize, 2 + fsize );
    GAMAT( *p_rssky_transf, 3 + fsize, 3 );
    gsl_matrix *GAMAT( interm_ssky_metric, 3 + fsize, 3 + fsize );

    // Compute the orbital metric in ecliptic coordinates
    gsl_matrix *orbital_metric = SM_ComputePhaseMetric( &ocoords, ref_time, segments, fiducial_freq, detectors, detector_weights, detector_motion, ephemerides );
    XLAL_CHECK( orbital_metric != NULL, XLAL_EFUNC );

    // Compute the reduced supersky metric from the unrestricted supersky metric and the orbital metric
    XLAL_CHECK( SM_ComputeFittedSuperskyMetric( interm_ssky_metric, *p_rssky_transf, ussky_metric, orbital_metric, &ocoords, spindowns, ref_time, segments ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( SM_ComputeDecoupledSuperskyMetric( interm_ssky_metric, *p_rssky_transf, interm_ssky_metric, spindowns ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( SM_ComputeAlignedSuperskyMetric( interm_ssky_metric, *p_rssky_transf, interm_ssky_metric, spindowns ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( SM_ExtractReducedSuperskyMetric( *p_rssky_metric, *p_rssky_transf, interm_ssky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Cleanup
    GFMAT( orbital_metric, interm_ssky_metric );

  }

  // Return or free unrestricted supersky metric
  if( p_ussky_metric != NULL ) {
    *p_ussky_metric = ussky_metric;
  } else {
    GFMAT( ussky_metric );
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
 | as[1] =    ___________________
 |     B = 1_|   _____   _____   |
 |           |  /     \ /     \  |
 |         0-| |       |       | |
 |          _|  \_____/ \_____/  |
 |        -1 |_.___.___.___.___._|
 |             '   '   '   '   '
 |        A = -2  -1   0   1   2
 |    as[0] =  1   0  -1   0   1
 |    as[2] =  0  -1   0   1   0
 \endverbatim
 *
 * Points outside the unit disks are moved radially onto their boundaries.
 */
static void SM_ReducedToAligned(
  double as[3],					///< [out] 3-dimensional aligned sky coordinates
  const gsl_vector *rss				///< [in] 2-dimensional reduced supersky coordinates
  )
{
  const double A = gsl_vector_get( rss, 0 );
  const double B = gsl_vector_get( rss, 1 );
  const double dA = fabs( A ) - 1.0;
  const double R = sqrt( SQR( dA ) + SQR( B ) );
  const double Rmax = GSL_MAX( 1.0, R );
  as[0] = dA / Rmax;
  as[1] = B / Rmax;
  as[2] = GSL_SIGN( A ) * RE_SQRT( 1.0 - DOT2( as, as ) );
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
  const double r = sqrt( DOT3( as, as ) );
  const double A = GSL_SIGN( as[2] ) * ( ( as[0] / r ) + 1.0 );
  const double B = as[1] / r;
  gsl_vector_set( rss, 0, A );
  gsl_vector_set( rss, 1, B );
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
  XLAL_CHECK( out < SC_MAX, XLAL_EINVAL );
  XLAL_CHECK( in < SC_MAX, XLAL_EINVAL );
  XLAL_CHECK( out_points != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_points != NULL, XLAL_EINVAL );
  XLAL_CHECK( rssky_transf != NULL || ( out != SC_RSSKY && in != SC_RSSKY ), XLAL_EINVAL );

  // Deduce number of input sky coordinates, and frequency/spindown coordinates
  const size_t in_ssize = ( in == SC_USSKY ) ? 3 : 2;
  XLAL_CHECK( in_points->size1 > in_ssize, XLAL_EINVAL );
  const size_t fsize = in_points->size1 - in_ssize;

  // Resize or allocate output points matrix, if required
  const size_t out_ssize = ( out == SC_USSKY ) ? 3 : 2;
  const size_t out_rows = fsize + out_ssize;
  if( *out_points != NULL ) {
    if( ( *out_points )->size1 != out_rows || ( *out_points )->size2 != in_points->size2 ) {
      gsl_matrix_free( *out_points );
      *out_points = NULL;
    }
  }
  if( *out_points == NULL ) {
    GAMAT( *out_points, out_rows, in_points->size2 );
  }

  // If input and output coordinate systems are the same, copy input matrix and exit
  if( in == out ) {
    gsl_matrix_memcpy( *out_points, in_points );
    return XLAL_SUCCESS;
  }

  // Iterate over input points
  for( size_t j = 0; j < in_points->size2; ++j ) {

    // Create array for point in intermediate coordinates
    double tmp[3 + fsize];
    gsl_vector_view tmp_sky = gsl_vector_view_array( &tmp[0], 3 );
    gsl_vector_view tmp_fspin = gsl_vector_view_array( &tmp[3], fsize );

    // Copy input point to intermediate point
    for( size_t i = 0; i < in_ssize; ++i ) {
      tmp[i] = gsl_matrix_get( in_points, i, j );
    }
    for( size_t i = 0; i < fsize; ++i ) {
      tmp[3 + i] = gsl_matrix_get( in_points, in_ssize + i, j );
    }

    // Initialise current coordinate system
    SuperskyCoordinates curr = in;

    // Convert physical coordinates to supersky coordinates
    if( curr == SC_PHYS && out > curr ) {

      // Convert right ascension and declination to supersky position
      const double alpha = tmp[0];
      const double delta = tmp[1];
      const double cos_delta = cos( delta );
      tmp[0] = cos( alpha ) * cos_delta;
      tmp[1] = sin( alpha ) * cos_delta;
      tmp[2] = sin( delta );

      // Update current coordinate system
      curr = SC_USSKY;

    }

    // Convert supersky coordinates to reduced supersky coordinates
    if( curr == SC_USSKY && out > curr ) {

      // Move frequency to after spindowns
      const double freq = tmp[3];
      memmove( &tmp[3], &tmp[4], ( fsize - 1 ) * sizeof( tmp[0] ) );
      tmp[2 + fsize] = freq;

      // Create views of the sky alignment transform and sky offset vectors
      gsl_matrix_const_view align_sky = gsl_matrix_const_submatrix( rssky_transf, 0, 0, 3, 3 );
      gsl_matrix_const_view sky_offsets = gsl_matrix_const_submatrix( rssky_transf, 3, 0, fsize, 3 );

      // Apply the alignment transform to the supersky position to produced the aligned sky position:
      //   asky = align_sky * ssky
      double asky[3];
      gsl_vector_view asky_v = gsl_vector_view_array( asky, 3 );
      gsl_blas_dgemv( CblasNoTrans, 1.0, &align_sky.matrix, &tmp_sky.vector, 0.0, &asky_v.vector );

      // Add the inner product of the sky offsets with the aligned sky position
      // to the supersky spins and frequency to get the reduced supersky quantities:
      //   rssky_fspin[i] = ussky_fspin[i] + dot(sky_offsets[i], asky)
      gsl_blas_dgemv( CblasNoTrans, 1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &tmp_fspin.vector );

      // Convert from 3-dimensional aligned sky coordinates to 2-dimensional reduced supersky coordinates
      SM_AlignedToReduced( &tmp_sky.vector, asky );

      // Update current coordinate system
      curr = SC_RSSKY;

    }

    // Convert reduced supersky coordinates to supersky coordinates
    if( curr == SC_RSSKY && out < curr ) {

      // Create views of the sky alignment transform and sky offset vectors
      gsl_matrix_const_view align_sky = gsl_matrix_const_submatrix( rssky_transf, 0, 0, 3, 3 );
      gsl_matrix_const_view sky_offsets = gsl_matrix_const_submatrix( rssky_transf, 3, 0, fsize, 3 );

      // Convert from 2-dimensional reduced supersky coordinates to 3-dimensional aligned sky coordinates
      double asky[3];
      SM_ReducedToAligned( asky, &tmp_sky.vector );
      gsl_vector_view asky_v = gsl_vector_view_array( asky, 3 );

      // Subtract the inner product of the sky offsets with the aligned sky position
      // from the reduced supersky spins and frequency to get the supersky quantities:
      //   ussky_fspin[i] = rssky_fspin[i] - dot(sky_offsets[i], asky)
      gsl_blas_dgemv( CblasNoTrans, -1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &tmp_fspin.vector );

      // Apply the inverse alignment transform to the aligned sky position to produced the supersky position:
      //   ssky = align_sky^T * asky
      gsl_blas_dgemv( CblasTrans, 1.0, &align_sky.matrix, &asky_v.vector, 0.0, &tmp_sky.vector );

      // Move frequency to before spindowns
      const double freq = tmp[2 + fsize];
      memmove( &tmp[4], &tmp[3], ( fsize - 1 ) * sizeof( tmp[0] ) );
      tmp[3] = freq;

      // Update current coordinate system
      curr = SC_USSKY;

    }

    // Convert supersky coordinates to physical coordinates
    if( curr == SC_USSKY && out < curr ) {

      // Convert supersky position to right ascension and declination
      const double nx = tmp[0];
      const double ny = tmp[1];
      const double nz = tmp[2];
      tmp[0] = atan2( ny, nx );
      tmp[1] = atan2( nz, sqrt( SQR( nx ) + SQR( ny ) ) );
      XLALNormalizeSkyPosition( &tmp[0], &tmp[1] );

      // Update current coordinate system
      curr = SC_PHYS;

    }

    // Check that correct coordinate system has been converted to
    XLAL_CHECK( curr == out, XLAL_EFAILED );

    // Copy intermediate point to output point
    for( size_t i = 0; i < out_ssize; ++i ) {
      gsl_matrix_set( *out_points, i, j, tmp[i] );
    }
    for( size_t i = 0; i < fsize; ++i ) {
      gsl_matrix_set( *out_points, out_ssize + i, j, tmp[3 + i] );
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
  XLAL_CHECK( SC_PHYS < out && out < SC_MAX, XLAL_EINVAL );
  XLAL_CHECK( out_point != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );

  // Deduce number of sky coordinates, and frequency/spindown coordinates
  const size_t ssize = ( out == SC_USSKY ) ? 3 : 2;
  XLAL_CHECK( out_point->size > ssize, XLAL_EINVAL );
  const size_t fsize = out_point->size - ssize;
  XLAL_CHECK( fsize <= PULSAR_MAX_SPINS, XLAL_EFAILED );

  // Transform input physical point to reference time of coordinate transform data
  PulsarDopplerParams in_phys_ref = *in_phys;
  {
    const REAL8 dtau = XLALGPSDiff ( ref_time, &in_phys_ref.refTime );
    XLAL_CHECK ( XLALExtrapolatePulsarSpins ( in_phys_ref.fkdot, in_phys_ref.fkdot, dtau ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Copy input physical point to array
  double in_point[2 + fsize];
  in_point[0] = in_phys_ref.Alpha;
  in_point[1] = in_phys_ref.Delta;
  memcpy( &in_point[2], in_phys_ref.fkdot, fsize * sizeof( in_point[0] ) );

  // Convert input physical point to output supersky coordinate point
  gsl_matrix_view out_point_view = gsl_matrix_view_vector( out_point, out_point->size, 1 );
  gsl_matrix_const_view in_point_view = gsl_matrix_const_view_array( in_point, 2 + fsize, 1 );
  gsl_matrix *out_point_view_ptr = &out_point_view.matrix;
  XLAL_CHECK( XLALConvertSuperskyCoordinates( out, &out_point_view_ptr, SC_PHYS, &in_point_view.matrix, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( out_point_view_ptr == &out_point_view.matrix, XLAL_EFAILED );

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
  XLAL_CHECK( out_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( SC_PHYS < in && in < SC_MAX, XLAL_EINVAL );
  XLAL_CHECK( in_point != NULL, XLAL_EFAULT );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );

  // Deduce number of sky coordinates, and frequency/spindown coordinates
  const size_t ssize = ( in == SC_USSKY ) ? 3 : 2;
  XLAL_CHECK( in_point->size > ssize, XLAL_EINVAL );
  const size_t fsize = in_point->size - ssize;
  XLAL_CHECK( fsize <= PULSAR_MAX_SPINS, XLAL_EFAILED );

  // Convert input supersky coordinate point to output physical point
  double out_point[2 + fsize];
  gsl_matrix_view out_point_view = gsl_matrix_view_array( out_point, 2 + fsize, 1 );
  gsl_matrix_const_view in_point_view = gsl_matrix_const_view_vector( in_point, in_point->size, 1 );
  gsl_matrix *out_point_view_ptr = &out_point_view.matrix;
  XLAL_CHECK( XLALConvertSuperskyCoordinates( SC_PHYS, &out_point_view_ptr, in, &in_point_view.matrix, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( out_point_view_ptr == &out_point_view.matrix, XLAL_EFAILED );

  // Copy output physical point from array
  out_phys->Alpha = out_point[0];
  out_phys->Delta = out_point[1];
  memcpy( out_phys->fkdot, &out_point[2], fsize * sizeof( out_point[0] ) );

  // Set output physical point reference time to that of of coordinate transform data
  out_phys->refTime = *ref_time;

  return XLAL_SUCCESS;

}

static double SuperskyAHemiFracRoot(
  double x,
  void *params
  )
{

  // Fractional area of reduced supersky hemisphere to left of line of constant 'x'
  const double frac = LAL_PI + x * sqrt( 1 - x*x ) - acos( x );

  // Target fractional area we are trying to find 'x' to satisfy
  const double target_frac = *( ( double * ) params );

  return frac - target_frac;

}

static double SuperskyBCoordBound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_vector *point
  )
{

  // Get bounds data
  const double semiB = *( ( const double * ) data );

  // Get 2-dimensional reduced supersky A coordinate
  const double A = gsl_vector_get( point, 0 );
  const double dA = fabs( A ) - 1.0;

  // Set bound on 2-dimensional reduced supersky B coordinate
  const double bound = semiB * RE_SQRT( 1.0 - SQR( dA ) );

  return bound;

}

int XLALSetSuperskyLatticeTilingAllSkyBounds(
  LatticeTiling *tiling,
  const double patch_B_extent,
  const UINT8 patch_count,
  const UINT8 patch_index
  )
{

  // Check input
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( patch_B_extent > 0.0, XLAL_EINVAL );
  XLAL_CHECK( patch_count > 0, XLAL_EINVAL );
  XLAL_CHECK( patch_index < patch_count, XLAL_EINVAL );

  // Calculate patch indexes and counts in reduced supersky coordinate A and B directions.
  // 'max_patch_B_extent' is used to calculate a minimum number of patches in the B coordinate
  // direction, 'min_patch_count_B', which with 'patch_count' is used to calculate the number
  // of patches in the A coordinate direction, 'patch_count_A'. The produce of these two counts
  // will never exceed 'patch_count'; any excess patches are then added to 'min_patch_count_B'
  // as required to make up 'patch_count'. For example, given
  //   'patch_count' = 14, 'patch_count_A' = 4, 'min_patch_count_B' = 3,
  // the number of patches is partitioned as [4, 4, 3, 3], i.e.
  //   'patch_index_A' = 0, 'patch_count_B' = 4, 'patch_index_B' = 0 to 3
  //   'patch_index_A' = 1, 'patch_count_B' = 4, 'patch_index_B' = 0 to 3
  //   'patch_index_A' = 2, 'patch_count_B' = 3, 'patch_index_B' = 0 to 2
  //   'patch_index_A' = 3, 'patch_count_B' = 3, 'patch_index_B' = 0 to 2

  // Approximate minimum number of patches in B coordinate
  const double approx_min_patch_count_B = 2.0 / patch_B_extent;

  // Number of patches in A coordinate
  const UINT8 patch_count_A = GSL_MAX( 1, lround( floor( patch_count / approx_min_patch_count_B ) ) );

  // Actual minimum number of patches in B coordinate; note integer division equivalent to floor()
  const UINT8 min_patch_count_B = patch_count / patch_count_A;

  // Excess number of patches which must be added on to get 'patch_count'
  INT8 patch_excess = patch_count - patch_count_A * min_patch_count_B;
  XLAL_CHECK( patch_excess >= 0, XLAL_EFAILED );

  // Initialise number of patches in B coordinate; if there are excess patches, add an extra patch
  UINT8 patch_count_B = min_patch_count_B;
  if( patch_excess > 0) {
    ++patch_count_B;
  }

  // Initialise patch indexes in A and B coordinates
  UINT8 patch_index_A = 0, patch_index_B = patch_index;

  while( patch_index_B >= patch_count_B ) {

    // Increase patch index in A coordinate, substract patch count in B coordinate from patch index
    ++patch_index_A;
    patch_index_B -= patch_count_B;

    // Decrease number of excess patches; if zero, subtract extra patch from patch count in B coordinate
    --patch_excess;
    if( patch_excess == 0 ) {
      --patch_count_B;
    }

  }

  // The reduced supersky A coordinate is bounded from 'skyA_bounds[0]' to 'skyA_bounds[1]'. The
  // value of 'skyA_bounds[i]' is chosen such that the fractional area of the sky (shaded # in the
  // diagram) to the left of the line A = 'skyA_bounds[i]' (dotted : vertical line in the diagram)
  // is equal to
  //   2*pi * ('patch_index_A' + 'i') / 'patch_count_A'.
  // This divides the sky into equal-area bands at constant 'patch_index_A'.
  //
  //        ______:____________
  // B = 1_|   ___:_   _____   |
  //       |  /###: \ /     \  |
  //     0-| |####:  |       | |
  //      _|  \###:_/ \_____/  |
  //    -1 |_.____:__._______._|
  //         '    :  '       '
  //    A = -2       0       2
  //
  // To find 'skyA_bounds[i]', we numerically solve for 'x' in [-1.0, 1.0] such that
  //   int_{-1}^{x} 2*sqrt(1 - x'^2) dx' = pi + x*sqrt(1 - x^2) - acos(x) = 'skyA_frac'
  // where 'skyA_frac' is the fractional area of one sky hemisphere. The bound is then given by
  //   'skyA_bounds[i]' = -1.0 + 2.0*'skyA_hemi' + 'x'
  // where 'skyA_hemi' is 0.0 for the left hemisphere and 1.0 for the right hemisphere.

  // Allocate GSL root solver for finding reduced supersky A coordinate bounds
  gsl_root_fsolver *GALLOC( skyA_bounds_fsolver, gsl_root_fsolver_alloc( gsl_root_fsolver_brent ) );

  // Compute the lower and upper bounds on reduced supersky A coordinate
  double skyA_bounds[2] = {0, 0};
  for( size_t i = 0; i < 2; ++i ) {
    const UINT8 iA = patch_index_A + i;

    // Treat special value of 'iA' where root-finding is not required separately
    if( iA == 0 ) {
      skyA_bounds[i] = -2.0;
      continue;
    }
    if( patch_count_A % 2 == 0 && iA == patch_count_A/2 ) {
      skyA_bounds[i] = 0.0;
      continue;
    }
    if( iA == patch_count_A ) {
      skyA_bounds[i] = 2.0;
      continue;
    }

    // Compute which hemisphere 'skyA_bounds[i]' is in ('skyA_hemi'), and
    // the fractional area of the hemisphere ('skyA_frac')
    double skyA_hemi = 0;
    double skyA_frac = LAL_PI * modf( 2.0 * ( ( double ) iA ) / patch_count_A, &skyA_hemi );

    // Initialise GSL root finder function and set parameter to 'skyA_frac'
    gsl_function skyA_bounds_fsolver_F = { .function = &SuperskyAHemiFracRoot, .params = &skyA_frac };

    // Set GSL root finder bounds on 'x' to [-1.0, 1.0] and solve for 'x'
    GCALL( gsl_root_fsolver_set( skyA_bounds_fsolver, &skyA_bounds_fsolver_F, -1.0, 1.0 ) );
    double x = -1.0;
    const double epsabs = 1e-4;
    int status = GSL_CONTINUE;
    int iterations = 100;
    while( status == GSL_CONTINUE && iterations-- > 0 ) {
      GCALL( gsl_root_fsolver_iterate( skyA_bounds_fsolver ) );
      x = gsl_root_fsolver_root( skyA_bounds_fsolver );
      const double x_lower = gsl_root_fsolver_x_lower( skyA_bounds_fsolver );
      const double x_upper = gsl_root_fsolver_x_upper( skyA_bounds_fsolver );
      status = gsl_root_test_interval( x_lower, x_upper, epsabs, 0.0 );
    }
    XLAL_CHECK( status == GSL_SUCCESS, XLAL_EMAXITER, "GSL root solver failed to converge: x=%0.6g, f(x)=%0.6g, epsabs=%0.6g", x, GSL_FN_EVAL(&skyA_bounds_fsolver_F, x), epsabs );

    // Set 'skyA_bounds[i]' from 'skyA_hemi' and 'x'
    skyA_bounds[i] = -1.0 + 2.0*skyA_hemi + x;

  }

  // Set the parameter-space bound on reduced supersky A coordinate
  XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, 0, skyA_bounds[0], skyA_bounds[1] ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingBoundPadding( tiling, 0, patch_index_A == 0, patch_index_A == patch_count_A - 1 ) == XLAL_SUCCESS, XLAL_EFUNC );

  // The reduced supersky B coordinate is bounded from below and above by ellipses, with semi-major
  // axes (in A) of 1.0, and semi-minor axes (in B) given by 'semiB', which has values:
  //   lower 'semiB' = -1.0 + 2.0 *   'patch_index_B'       / 'patch_count_B'
  //   upper 'semiB' = -1.0 + 2.0 * ( 'patch_index_B' + 1 ) / 'patch_count_B'
  // This divides the sky into equal-area bands at constant 'patch_index_B'.

  // Allocate memory
  const size_t data_len = sizeof( double );
  double *data_lower = XLALMalloc( data_len );
  XLAL_CHECK( data_lower != NULL, XLAL_ENOMEM );
  double *data_upper = XLALMalloc( data_len );
  XLAL_CHECK( data_upper != NULL, XLAL_ENOMEM );

  // Set the parameter-space bound on reduced supersky B coordinate
  data_lower[0] = -1.0 + 2.0 * ( ( double ) patch_index_B + 0 ) / patch_count_B;
  data_upper[0] = -1.0 + 2.0 * ( ( double ) patch_index_B + 1 ) / patch_count_B;
  XLAL_CHECK( XLALSetLatticeTilingBound( tiling, 1, SuperskyBCoordBound, data_len, data_lower, data_upper ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingBoundPadding( tiling, 1, patch_index_B == 0, patch_index_B == patch_count_B - 1 ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Cleanup
  gsl_root_fsolver_free( skyA_bounds_fsolver );

  return XLAL_SUCCESS;

}

int XLALSetSuperskyLatticeTilingSkyPointBounds(
  LatticeTiling *tiling,
  const gsl_matrix *rssky_transf,
  const double alpha,
  const double delta
  )
{

  // Check input
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf != NULL, XLAL_EFAULT );

  // Allocate memory
  gsl_vector *GAVEC( rssky_point, rssky_transf->size1 - 1 );

  // Convert right ascension and declination to reduced supersky coordinates
  PulsarDopplerParams XLAL_INIT_DECL( doppler );
  doppler.Alpha = alpha;
  doppler.Delta = delta;
  XLAL_CHECK( XLALConvertPhysicalToSupersky( SC_RSSKY, rssky_point, &doppler, rssky_transf, &doppler.refTime ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Set the parameter-space bounds on 2-dimensional reduced supersky A and B coordinates
  for( size_t i = 0; i < 2; ++i ) {
    const double rssky_point_i = gsl_vector_get( rssky_point, i );
    XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, i, rssky_point_i, rssky_point_i ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Cleanup
  GFVEC( rssky_point );

  return XLAL_SUCCESS;

}

static double PhysicalSpinBound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_vector *point
  )
{

  // Get bounds data
  const double *sky_offsets = ( ( const double * ) data );
  double bound = ( ( const double * ) data )[3];

  // Add the inner product of the sky offsets with the aligned sky
  // position to the physical bound to get the reduced supersky bound
  double as[3];
  SM_ReducedToAligned( as, point );
  bound += DOT3( sky_offsets, as );

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
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf->size1 > 3, XLAL_EINVAL );
  XLAL_CHECK( isfinite( bound1 ), XLAL_EINVAL );
  XLAL_CHECK( isfinite( bound2 ), XLAL_EINVAL );
  const size_t smax = rssky_transf->size1 - 4;
  XLAL_CHECK( s <= smax, XLAL_ESIZE );
  const size_t dim = ( s == 0 ) ? ( 2 + smax ) : ( 1 + s );

  // Allocate memory
  const size_t data_len = 4 * sizeof( double );
  double *data_lower = XLALMalloc( data_len );
  XLAL_CHECK( data_lower != NULL, XLAL_ENOMEM );
  double *data_upper = XLALMalloc( data_len );
  XLAL_CHECK( data_upper != NULL, XLAL_ENOMEM );

  // Copy the sky offset vector to bounds data
  for( size_t j = 0; j < 3; ++j ) {
    data_lower[j] = data_upper[j] = gsl_matrix_get( rssky_transf, dim + 1, j );
  }

  // Set the parameter-space bound on physical frequency/spindown coordinate
  data_lower[3] = GSL_MIN( bound1, bound2 );
  data_upper[3] = GSL_MAX( bound1, bound2 );
  XLAL_CHECK( XLALSetLatticeTilingBound( tiling, dim, PhysicalSpinBound, data_len, data_lower, data_upper ) == XLAL_SUCCESS, XLAL_EFUNC );

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
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf->size1 > 3, XLAL_EINVAL );
  XLAL_CHECK( isfinite( bound1 ), XLAL_EINVAL );
  XLAL_CHECK( isfinite( bound2 ), XLAL_EINVAL );
  const size_t smax = rssky_transf->size1 - 4;
  XLAL_CHECK( s <= smax, XLAL_ESIZE );
  const size_t dim = ( s == 0 ) ? ( 2 + smax ) : ( 1 + s );

  // Set the parameter-space bound on reduced supersky frequency/spindown coordinate
  XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, dim, bound1, bound2 ) == XLAL_SUCCESS, XLAL_EFUNC );

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
  XLAL_CHECK( spin_range != NULL, XLAL_EFAULT );
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf != NULL, XLAL_EFAULT );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );

  // Get rectange containing range reduced supersky coordinates
  double skyA_rect[4], skyB_rect[4];
  {
    const LatticeTilingStats *stats = XLALLatticeTilingStatistics( tiling, 0 );
    XLAL_CHECK( stats != NULL, XLAL_EFUNC );
    skyA_rect[0] = stats->min_value_pass;
    skyA_rect[1] = stats->max_value_pass;
    skyA_rect[2] = stats->max_value_pass;
    skyA_rect[3] = stats->min_value_pass;
  }
  {
    const LatticeTilingStats *stats = XLALLatticeTilingStatistics( tiling, 1 );
    XLAL_CHECK( stats != NULL, XLAL_EFUNC );
    skyB_rect[0] = stats->min_value_pass;
    skyB_rect[1] = stats->min_value_pass;
    skyB_rect[2] = stats->max_value_pass;
    skyB_rect[3] = stats->max_value_pass;
  }

  // Get range of physical frequency/spindowns
  const size_t smax = rssky_transf->size1 - 4;
  for( size_t i = 0; i < 4; ++i ) {

    // Construct reduced supersky point
    double in_rssky[3 + smax];
    in_rssky[0] = skyA_rect[i];
    in_rssky[1] = skyB_rect[i];
    for( size_t s = 0; s <= smax; ++s ) {
      const size_t dim = ( s == 0 ) ? ( 2 + smax ) : ( 1 + s );
      const LatticeTilingStats *stats = XLALLatticeTilingStatistics( tiling, dim );
      XLAL_CHECK( stats != NULL, XLAL_EFUNC );
      in_rssky[dim] = stats->min_value_pass;
    }

    // Convert reduced supersky point to physical coordinates
    gsl_vector_view in_rssky_view = gsl_vector_view_array( in_rssky, 3 + smax );
    PulsarDopplerParams XLAL_INIT_DECL( out_phys );
    XLAL_CHECK( XLALConvertSuperskyToPhysical( &out_phys, SC_RSSKY, &in_rssky_view.vector, rssky_transf, ref_time ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Store minimum/maximum physical frequency/spindown in 'spin_range'
    for( size_t s = 0; s <= smax; ++s ) {
      if( i == 0 || out_phys.fkdot[s] < spin_range->fkdot[s] ) {
        spin_range->fkdot[s] = out_phys.fkdot[s];
      }
      if( i == 0 || out_phys.fkdot[s] > spin_range->fkdotBand[s] ) {
        spin_range->fkdotBand[s] = out_phys.fkdot[s];
      }
    }

  }
  for( size_t s = 0; s <= smax; ++s ) {
    spin_range->fkdotBand[s] -= spin_range->fkdot[s];
  }

  // Adjust 'spin_range' bands to include width of supersky frequency/spindown parameter space
  for( size_t s = 0; s <= smax; ++s ) {
    const size_t dim = ( s == 0 ) ? ( 2 + smax ) : ( 1 + s );
    const LatticeTilingStats *stats = XLALLatticeTilingStatistics( tiling, dim );
    XLAL_CHECK( stats != NULL, XLAL_EFUNC );
    spin_range->fkdotBand[s] += stats->max_value_pass - stats->min_value_pass;
  }

  // Set reference time of 'spin_range' to that of coordinate transform data
  spin_range->refTime = *ref_time;

  return XLAL_SUCCESS;

}
