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
#include <lal/LogPrintf.h>
#include <lal/MetricUtils.h>
#include <lal/ExtrapolatePulsarSpins.h>

#include <lal/GSLHelpers.h>

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

// Dimensions/indexes into reduced supersky coordinate transform data
#define _RSSKY_TRANSF_XTRADIM   2
#define _RSSKY_TRANSF_REFTIME   0, 0
#define _RSSKY_TRANSF_FIDFREQ   0, 1

// Check reduced supersky coordinate metric and/or transform data
#define CHECK_RSSKY_METRIC_TRANSF(M, RT) \
  ((M) != NULL && CHECK_RSSKY_TRANSF(RT) && (M)->size1 == (RT)->size1 - _RSSKY_TRANSF_XTRADIM)
#define CHECK_RSSKY_TRANSF(RT) \
  ((RT) != NULL && (RT)->size1 > (2 + _RSSKY_TRANSF_XTRADIM) && (RT)->size2 == 3)

// Decompose reduced supersky coordinate transform data
#define _DECOMPOSE_RSSKY_TRANSF(RT, GSLMAT) \
  UNUSED const size_t ndim = (RT)->size1 - _RSSKY_TRANSF_XTRADIM; \
  UNUSED const size_t smax = ndim - 3; \
  LIGOTimeGPS _decomposed_ref_time; \
  UNUSED const LIGOTimeGPS *ref_time = XLALGPSSetREAL8(&_decomposed_ref_time, gsl_matrix_get((RT), _RSSKY_TRANSF_REFTIME)); \
  UNUSED const double fiducial_freq = gsl_matrix_get((RT), _RSSKY_TRANSF_FIDFREQ); \
  UNUSED GSLMAT##_view align_sky = GSLMAT##_submatrix((RT), 1, 0, 3, 3); \
  UNUSED GSLMAT##_view sky_offsets = GSLMAT##_submatrix((RT), 2 + _RSSKY_TRANSF_XTRADIM, 0, 1 + smax, 3); \
  do { } while(0)
#define DECOMPOSE_RSSKY_TRANSF(RT)         _DECOMPOSE_RSSKY_TRANSF(RT, gsl_matrix)
#define DECOMPOSE_CONST_RSSKY_TRANSF(RT)   _DECOMPOSE_RSSKY_TRANSF(RT, gsl_matrix_const)

// Determine which dimension stores the reduced supersky frequency/spindown of order 's'
#define RSSKY_FKDOT_OFFSET(s)   (((s) == 0) ? (smax) : ((size_t)((s) - 1)))
#define RSSKY_FKDOT_DIM(s)      (2 + RSSKY_FKDOT_OFFSET(s))

///
/// Fiducial frequency at which to numerically calculate metrics, which
/// are then rescaled to user-requested frequency based on known scalings
///
const double fiducial_calc_freq = 100.0;

///
/// Call XLALComputeDopplerPhaseMetric() to compute the phase metric for a given coordinate system.
///
static gsl_matrix *SM_ComputePhaseMetric(
  const DopplerCoordinateSystem *coords,        ///< [in] Coordinate system to compute metric for
  const LIGOTimeGPS *ref_time,                  ///< [in] Reference time of the metric
  const LIGOTimeGPS *start_time,                ///< [in] Start time of the metric
  const LIGOTimeGPS *end_time,                  ///< [in] End time of the metric
  const MultiLALDetector *detectors,            ///< [in] List of detectors to average metric over
  const MultiNoiseFloor *detector_weights,      ///< [in] Weights used to combine single-detector metrics (default: unit weights)
  const DetectorMotionType detector_motion,     ///< [in] Which detector motion to use
  const EphemerisData *ephemerides              ///< [in] Earth/Sun ephemerides
  )
{

  // Check input
  XLAL_CHECK_NULL( coords != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( start_time != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( end_time != NULL, XLAL_EFAULT );
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
  LALSegList segments;
  XLAL_CHECK_NULL( XLALSegListInit( &segments ) == XLAL_SUCCESS, XLAL_EFUNC );
  LALSeg segment;
  XLAL_CHECK_NULL( XLALSegSet( &segment, start_time, end_time, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_NULL( XLALSegListAppend( &segments, &segment ) == XLAL_SUCCESS, XLAL_EFUNC );
  par.segmentList = segments;

  // Set detectors and detector weights
  par.multiIFO = *detectors;
  if ( detector_weights != NULL ) {
    par.multiNoiseFloor = *detector_weights;
  } else {
    par.multiNoiseFloor.length = 0;   // Indicates unit weights
  }

  // Set reference time and fiducial frequency
  par.signalParams.Doppler.refTime = *ref_time;
  par.signalParams.Doppler.fkdot[0] = fiducial_calc_freq;

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
  XLALSegListClear( &segments );

  return g_ij;

}

///
/// Find a least-squares linear fit to the orbital X and Y metric elements using the frequency and
/// spindown metric elements. Outputs the intermediate \e fitted supersky metric, and updates the
/// reduced supersky metric coordinate transform data.
///
static int SM_ComputeFittedSuperskyMetric(
  gsl_matrix *fitted_ssky_metric,               ///< [out] Fitted supersky metric
  gsl_matrix *rssky_transf,                     ///< [in,out] Reduced supersky metric coordinate transform data
  const gsl_matrix *ussky_metric,               ///< [in] Unrestricted supersky metric
  const gsl_matrix *orbital_metric,             ///< [in] Orbital metric in ecliptic coordinates
  const DopplerCoordinateSystem *ocoords,       ///< [in] Coordinate system of orbital metric
  const LIGOTimeGPS *start_time,                ///< [in] Start time of the metrics
  const LIGOTimeGPS *end_time                   ///< [in] End time of the metrics
  )
{

  // Check input
  XLAL_CHECK( fitted_ssky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( ussky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( orbital_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( ocoords != NULL, XLAL_EFAULT );
  XLAL_CHECK( start_time != NULL, XLAL_EFAULT );
  XLAL_CHECK( end_time != NULL, XLAL_EFAULT );

  // Decompose coordinate transform data
  DECOMPOSE_RSSKY_TRANSF( rssky_transf );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );

  // Allocate memory
  gsl_matrix *GAMAT( tmp, 3 + smax, 3 + smax );
  gsl_vector *GAVEC( tmpv, 1 + smax );

  // Compute mid-time of segment list
  LIGOTimeGPS mid_time = *start_time;
  XLALGPSAdd( &mid_time, 0.5 * XLALGPSDiff( end_time, start_time ) );

  // Internal copy of orbital metric, and various transforms performed on it
  gsl_matrix *orb_metric = NULL, *mid_time_transf = NULL, *diag_norm_transf = NULL;

  // Transform reference time of orbital metric from reference time to segment list mid-time
  const REAL8 Dtau = XLALGPSDiff( &mid_time, ref_time );
  XLAL_CHECK( XLALChangeMetricReferenceTime( &orb_metric, &mid_time_transf, orbital_metric, ocoords, Dtau ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Diagonally-normalize orbital metric
  XLAL_CHECK( XLALDiagNormalizeMetric( &orb_metric, &diag_norm_transf, orb_metric ) == XLAL_SUCCESS, XLAL_EFUNC );

  // 'fitA' contains the frequency and spindown elements of the orbital metric, used for fitting
  gsl_matrix *GAMAT( fitA, 3 + smax, 1 + smax );
  {
    gsl_matrix_view orb_metric_fspin = gsl_matrix_submatrix( orb_metric, 0, 2, 3 + smax, 1 + smax );
    gsl_matrix_memcpy( fitA, &orb_metric_fspin.matrix );
  }

  // Compute 'fitA^T * fitA'
  gsl_matrix *GAMAT( fitAt_fitA, 1 + smax, 1 + smax );
  gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, fitA, fitA, 0.0, fitAt_fitA );

  // Find the singular value decomposition of 'fitA^T * fitA'
  gsl_matrix *GAMAT( svd_U, 1 + smax, 1 + smax );
  gsl_matrix *GAMAT( svd_V, 1 + smax, 1 + smax );
  gsl_vector *GAVEC( svd_S, 1 + smax );
  gsl_matrix_memcpy( svd_U, fitAt_fitA );
  XLAL_CHECK( gsl_linalg_SV_decomp( svd_U, svd_V, svd_S, tmpv ) == 0, XLAL_EFAILED );

  // The columns of 'fitc' contain the least-square fitting coefficients for the orbital X and Y metric elements:
  //    fitc(:,j) = inv(fitA^T * fitA) * fitA^T * orb_metric(:,j)
  // The singular decomposition of fitA^T * fitA is used for the inverse
  gsl_matrix *GAMAT( fitc, 1 + smax, 2 );
  for ( size_t j = 0; j < 2; ++j ) {
    gsl_vector_view orb_metric_j = gsl_matrix_column( orb_metric, j );
    gsl_vector_view fitc_j = gsl_matrix_column( fitc, j );
    gsl_blas_dgemv( CblasTrans, 1.0, fitA, &orb_metric_j.vector, 0.0, tmpv );
    XLAL_CHECK( gsl_linalg_SV_solve( svd_U, svd_V, svd_S, tmpv, &fitc_j.vector ) == 0, XLAL_EFAILED );
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
  gsl_matrix *GAMAT( subtract_orb, 3 + smax, 3 + smax );
  {
    gsl_matrix_set_identity( subtract_orb );
    gsl_matrix_view subtract_orb_fspin_sky = gsl_matrix_submatrix( subtract_orb, 2, 0, 1 + smax, 2 );
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
  gsl_matrix *GAMAT( subtract_ussky, 4 + smax, 4 + smax );
  {
    gsl_matrix_set_identity( subtract_ussky );
    gsl_matrix_view subtract_ussky_fspin_sky = gsl_matrix_submatrix( subtract_ussky, 3, 0, 1 + smax, 3 );
    gsl_matrix_view subtract_orb_fspin_sky = gsl_matrix_submatrix( subtract_orb, 2, 0, 1 + smax, 2 );
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
    gsl_matrix_view subtract_ussky_fspin_sky = gsl_matrix_submatrix( subtract_ussky, 3, 0, 1 + smax, 3 );
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
  gsl_matrix *decoupled_ssky_metric,            ///< [out] Decoupled supersky metric
  gsl_matrix *rssky_transf,                     ///< [in,out] Reduced supersky metric coordinate transform data
  const gsl_matrix *fitted_ssky_metric          ///< [in] Fitted supersky metric
  )
{

  // Check input
  XLAL_CHECK( decoupled_ssky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( fitted_ssky_metric != NULL, XLAL_EFAULT );

  // Decompose coordinate transform data
  DECOMPOSE_RSSKY_TRANSF( rssky_transf );

  // Copy fitted metric to decoupled metric
  gsl_matrix_memcpy( decoupled_ssky_metric, fitted_ssky_metric );

  // Create views of the sky--sky, freq+spin--freq+spin, and off-diagonal blocks
  gsl_matrix_view sky_sky     = gsl_matrix_submatrix( decoupled_ssky_metric, 0, 0, 3, 3 );
  gsl_matrix_view sky_fspin   = gsl_matrix_submatrix( decoupled_ssky_metric, 0, 3, 3, 1 + smax );
  gsl_matrix_view fspin_sky   = gsl_matrix_submatrix( decoupled_ssky_metric, 3, 0, 1 + smax, 3 );
  gsl_matrix_view fspin_fspin = gsl_matrix_submatrix( decoupled_ssky_metric, 3, 3, 1 + smax, 1 + smax );

  // Diagonal-normalise the freq+spin--freq+spin block
  gsl_matrix *fspin_fspin_dnorm = NULL, *fspin_fspin_dnorm_transf = NULL;
  XLAL_CHECK( XLALDiagNormalizeMetric( &fspin_fspin_dnorm, &fspin_fspin_dnorm_transf, &fspin_fspin.matrix ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Invert the freq+spin--freq+spin block
  gsl_matrix *GAMAT( fspin_fspin_dnorm_LU, 1 + smax, 1 + smax );
  gsl_matrix *GAMAT( fspin_fspin_dnorm_inv, 1 + smax, 1 + smax );
  gsl_permutation *GAPERM( fspin_fspin_dnorm_LU_perm, 1 + smax );
  int fspin_fspin_dnorm_LU_sign = 0;
  gsl_matrix_memcpy( fspin_fspin_dnorm_LU, fspin_fspin_dnorm );
  XLAL_CHECK( gsl_linalg_LU_decomp( fspin_fspin_dnorm_LU, fspin_fspin_dnorm_LU_perm, &fspin_fspin_dnorm_LU_sign ) == 0, XLAL_EFAILED );
  XLAL_CHECK( gsl_linalg_LU_invert( fspin_fspin_dnorm_LU, fspin_fspin_dnorm_LU_perm, fspin_fspin_dnorm_inv ) == 0, XLAL_EFAILED );

  // Compute the additional sky offsets required to decouple the sky--sky and frequency blocks:
  //   decouple_sky_offsets = fspin_fspin_dnorm_transf * inv(fspin_fspin_dnorm) * fspin_fspin_dnorm_transf * fspin_sky
  // Uses fspin_sky as a temporary matrix, since it will be zeroed out anyway
  gsl_matrix *GAMAT( decouple_sky_offsets, 1 + smax, 3 );
  gsl_blas_dtrmm( CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, fspin_fspin_dnorm_transf, &fspin_sky.matrix );
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, fspin_fspin_dnorm_inv, &fspin_sky.matrix, 0.0, decouple_sky_offsets );
  gsl_blas_dtrmm( CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, fspin_fspin_dnorm_transf, decouple_sky_offsets );

  // Add the additional sky offsets to the reduced supersky coordinate transform data
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
  gsl_matrix *aligned_ssky_metric,              ///< [out] Aligned supersky metric
  gsl_matrix *rssky_transf,                     ///< [in,out] Reduced supersky metric coordinate transform data
  const gsl_matrix *decoupled_ssky_metric       ///< [in] Decoupled supersky metric
  )
{

  // Check input
  XLAL_CHECK( aligned_ssky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( decoupled_ssky_metric != NULL, XLAL_EFAULT );

  // Decompose coordinate transform data
  DECOMPOSE_RSSKY_TRANSF( rssky_transf );

  // Allocate memory
  gsl_matrix *GAMAT( tmp, 1 + smax, 3 );

  // Copy decoupled metric to aligned metric
  gsl_matrix_memcpy( aligned_ssky_metric, decoupled_ssky_metric );

  // Compute the eigenvalues/vectors of the sky--sky block
  gsl_vector *GAVEC( sky_eval, 3 );
  gsl_matrix *GAMAT( sky_evec, 3, 3 );
  gsl_eigen_symmv_workspace *GALLOC( wksp, gsl_eigen_symmv_alloc( 3 ) );
  gsl_matrix_view sky_sky = gsl_matrix_submatrix( aligned_ssky_metric, 0, 0, 3, 3 );
  XLAL_CHECK( gsl_eigen_symmv( &sky_sky.matrix, sky_eval, sky_evec, wksp ) == 0, XLAL_EFAILED );

  // Sort the eigenvalues/vectors by descending absolute eigenvalue
  XLAL_CHECK( gsl_eigen_symmv_sort( sky_eval, sky_evec, GSL_EIGEN_SORT_ABS_DESC ) == 0, XLAL_EFAILED );

  // Set the sky--sky block to the diagonal matrix of eigenvalues
  gsl_matrix_set_zero( &sky_sky.matrix );
  gsl_vector_view sky_sky_diag = gsl_matrix_diagonal( &sky_sky.matrix );
  gsl_vector_memcpy( &sky_sky_diag.vector, sky_eval );

  // Ensure that the matrix of eigenvalues has a positive diagonal; this and
  // the determinant constraints ensures fully constraints the eigenvector signs
  for ( size_t j = 0; j < 3; ++j ) {
    gsl_vector_view col = gsl_matrix_column( sky_evec, j );
    if ( gsl_vector_get( &col.vector, j ) < 0.0 ) {
      gsl_vector_scale( &col.vector, -1.0 );
    }
  }

  // Store the alignment transform in the reduced supersky coordinate transform data
  gsl_matrix_transpose_memcpy( &align_sky.matrix, sky_evec );

  // Ensure that the alignment transform has a positive determinant,
  // to ensure that that it represents a rotation
  gsl_permutation *GAPERM( LU_perm, 3 );
  int LU_sign = 0;
  XLAL_CHECK( gsl_linalg_LU_decomp( sky_evec, LU_perm, &LU_sign ) == 0, XLAL_EFAILED );
  if ( gsl_linalg_LU_det( sky_evec, LU_sign ) < 0.0 ) {
    gsl_vector_view col = gsl_matrix_column( &align_sky.matrix, 2 );
    gsl_vector_scale( &col.vector, -1.0 );
  }

  // Multiply the sky offsets by the alignment transform to transform to aligned sky coordinates:
  //   aligned_sky_off = sky_offsets * alignsky^T;
  gsl_matrix_memcpy( tmp, &sky_offsets.matrix );
  gsl_blas_dgemm( CblasNoTrans, CblasTrans, 1.0, tmp, &align_sky.matrix, 0.0, &sky_offsets.matrix );

  // Cleanup
  gsl_eigen_symmv_free( wksp );
  GFMAT( sky_evec, tmp );
  GFPERM( LU_perm );
  GFVEC( sky_eval );

  return XLAL_SUCCESS;

}

///
/// Compute the reduced supersky metric
///
static int SM_ComputeReducedSuperskyMetric(
  gsl_matrix **rssky_metric,                    ///< [out] Reduced supersky metric
  gsl_matrix **rssky_transf,                    ///< [out] Reduced supersky metric coordinate transform data
  const gsl_matrix *ussky_metric,               ///< [in] Unrestricted supersky metric
  const DopplerCoordinateSystem *ucoords,       ///< [in] Coordinate system of unrestricted supersky metric
  const gsl_matrix *orbital_metric,             ///< [in] Orbital metric in ecliptic coordinates
  const DopplerCoordinateSystem *ocoords,       ///< [in] Coordinate system of orbital metric
  const LIGOTimeGPS *ref_time,                  ///< [in] Reference time of the metrics
  const LIGOTimeGPS *start_time,                ///< [in] Start time of the metrics
  const LIGOTimeGPS *end_time                   ///< [in] End time of the metrics
  )
{

  // Check input
  XLAL_CHECK( rssky_metric != NULL && *rssky_metric == NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_transf != NULL && *rssky_transf == NULL, XLAL_EFAULT );
  XLAL_CHECK( ussky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( ussky_metric->size1 == ussky_metric->size2, XLAL_ESIZE );
  XLAL_CHECK( ucoords != NULL, XLAL_EFAULT );
  XLAL_CHECK( orbital_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( orbital_metric->size1 == orbital_metric->size2, XLAL_ESIZE );
  XLAL_CHECK( ocoords != NULL, XLAL_EFAULT );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK( start_time != NULL, XLAL_EFAULT );
  XLAL_CHECK( end_time != NULL, XLAL_EFAULT );

  // Allocate memory
  GAMAT( *rssky_metric, orbital_metric->size1, orbital_metric->size1 );
  GAMAT( *rssky_transf, _RSSKY_TRANSF_XTRADIM + orbital_metric->size1, 3 );
  gsl_matrix *GAMAT( assky_metric, 1 + orbital_metric->size1, 1 + orbital_metric->size1 );

  // Set coordinate transform metadata
  gsl_matrix_set( *rssky_transf, _RSSKY_TRANSF_REFTIME, XLALGPSGetREAL8( ref_time ) );
  gsl_matrix_set( *rssky_transf, _RSSKY_TRANSF_FIDFREQ, fiducial_calc_freq );

  // Compute the aligned supersky metric from the unrestricted supersky metric and the orbital metric
  XLAL_CHECK( SM_ComputeFittedSuperskyMetric( assky_metric, *rssky_transf, ussky_metric, orbital_metric, ocoords, start_time, end_time ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( SM_ComputeDecoupledSuperskyMetric( assky_metric, *rssky_transf, assky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( SM_ComputeAlignedSuperskyMetric( assky_metric, *rssky_transf, assky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );

  {
    // Move the row/column of the aligned supersky metric corresponding to the frequency to the last row/column
    // Move the row of the coordinate transform data corresponding to the frequency to the last row
    const int ifreq = XLALFindDopplerCoordinateInSystem( ucoords, DOPPLERCOORD_FREQ );
    XLAL_CHECK( 0 <= ifreq, XLAL_EFAILED );
    for ( size_t i = ifreq; i + 1 < assky_metric->size1; ++i ) {
      gsl_matrix_swap_rows( assky_metric, i, i + 1 );
      gsl_matrix_swap_columns( assky_metric, i, i + 1 );
      gsl_matrix_swap_rows( *rssky_transf, i + _RSSKY_TRANSF_XTRADIM - 1, i + _RSSKY_TRANSF_XTRADIM );
    }

    {
      // Move the row/column of the aligned supersky metric corresponding to the 'n_c' sky coordinate to the last row/column, so it can be easily dropped
      const int inc = XLALFindDopplerCoordinateInSystem( ucoords, DOPPLERCOORD_N3Z_EQU );
      XLAL_CHECK( 0 <= inc && inc < ifreq, XLAL_EFAILED );
      for ( size_t i = inc; i + 1 < assky_metric->size1; ++i ) {
        gsl_matrix_swap_rows( assky_metric, i, i + 1 );
        gsl_matrix_swap_columns( assky_metric, i, i + 1 );
      }
    }
  }

  // Copy all but the last row/column to the aligned supersky metric to the reduced supersky metric, dropping 'n_c'
  gsl_matrix_view extract_rssky_metric = gsl_matrix_submatrix( assky_metric, 0, 0, ( *rssky_metric )->size1, ( *rssky_metric )->size1 );
  gsl_matrix_memcpy( *rssky_metric, &extract_rssky_metric.matrix );

  // Ensure reduced supersky metric is symmetric
  for ( size_t i = 0; i < ( *rssky_metric )->size1; ++i ) {
    for ( size_t j = i + 1; j < ( *rssky_metric )->size1; ++j ) {
      const double gij = gsl_matrix_get( *rssky_metric, i, j );
      const double gji = gsl_matrix_get( *rssky_metric, j, i );
      const double g = 0.5 * ( gij + gji );
      gsl_matrix_set( *rssky_metric, i, j, g );
      gsl_matrix_set( *rssky_metric, j, i, g );
    }
  }

  // Ensure reduced supersky metric is positive definite
  for ( size_t s = 1; s <= ( *rssky_metric )->size1; ++s ) {
    gsl_matrix_view rssky_metric_s = gsl_matrix_submatrix( *rssky_metric, 0, 0, s, s );
    const double det_s = XLALMetricDeterminant( &rssky_metric_s.matrix );
    XLAL_CHECK( det_s > 0, XLAL_EFAILED, "Reduced supersky metric is not positive definite (s=%zu, det_s=%0.3e)", s, det_s );
  }

  // Cleanup
  GFMAT( assky_metric );

  return XLAL_SUCCESS;

}

SuperskyMetrics *XLALComputeSuperskyMetrics(
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
  XLAL_CHECK_NULL( spindowns <= 3, XLAL_EINVAL );
  XLAL_CHECK_NULL( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( segments != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( XLALSegListIsInitialized( segments ), XLAL_EINVAL );
  XLAL_CHECK_NULL( segments->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( fiducial_freq > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( detectors != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( detectors->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( detector_motion > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( ephemerides != NULL, XLAL_EINVAL );

  // Build coordinate system for the unrestricted supersky metric and orbital metric
  DopplerCoordinateSystem XLAL_INIT_DECL( ucoords );
  DopplerCoordinateSystem XLAL_INIT_DECL( ocoords );
  {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_N3X_EQU;
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_N3Y_EQU;
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_N3Z_EQU;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_N3OX_ECL;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_N3OY_ECL;
  }
  {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_FREQ;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_FREQ;
  }
  if ( spindowns >= 1 ) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F1DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F1DOT;
  }
  if ( spindowns >= 2 ) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F2DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F2DOT;
  }
  if ( spindowns >= 3 ) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F3DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F3DOT;
  }

  // Allocate memory for output struct
  SuperskyMetrics *metrics = XLALCalloc( 1, sizeof( *metrics ) );
  XLAL_CHECK_NULL( metrics != NULL, XLAL_ENOMEM );
  metrics->num_segments = segments->length;

  // Allocate memory for arrays of coherent metrics
  metrics->coh_rssky_metric = XLALCalloc( metrics->num_segments, sizeof( *metrics->coh_rssky_metric ) );
  XLAL_CHECK_NULL( metrics->coh_rssky_metric != NULL, XLAL_ENOMEM );
  metrics->coh_rssky_transf = XLALCalloc( metrics->num_segments, sizeof( *metrics->coh_rssky_transf ) );
  XLAL_CHECK_NULL( metrics->coh_rssky_transf != NULL, XLAL_ENOMEM );

  // Allocate memory for averaged metrics
  gsl_matrix *GAMAT_NULL( ussky_metric_avg, 4 + spindowns, 4 + spindowns );
  gsl_matrix *GAMAT_NULL( orbital_metric_avg, 3 + spindowns, 3 + spindowns );

  // Compute the coherent supersky metrics for each segment
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    const LIGOTimeGPS *start_time_seg = &segments->segs[n].start;
    const LIGOTimeGPS *end_time_seg = &segments->segs[n].end;

    // Compute the unrestricted supersky metric
    gsl_matrix *ussky_metric_seg = SM_ComputePhaseMetric( &ucoords, ref_time, start_time_seg, end_time_seg, detectors, detector_weights, detector_motion, ephemerides );
    XLAL_CHECK_NULL( ussky_metric_seg != NULL, XLAL_EFUNC );
    gsl_matrix_add( ussky_metric_avg, ussky_metric_seg );

    // Compute the orbital metric in ecliptic coordinates
    gsl_matrix *orbital_metric_seg = SM_ComputePhaseMetric( &ocoords, ref_time, start_time_seg, end_time_seg, detectors, detector_weights, detector_motion, ephemerides );
    XLAL_CHECK_NULL( orbital_metric_seg != NULL, XLAL_EFUNC );
    gsl_matrix_add( orbital_metric_avg, orbital_metric_seg );

    // Compute the coherent reduced supersky metric
    XLAL_CHECK_NULL( SM_ComputeReducedSuperskyMetric( &metrics->coh_rssky_metric[n], &metrics->coh_rssky_transf[n], ussky_metric_seg, &ucoords, orbital_metric_seg, &ocoords, ref_time, start_time_seg, end_time_seg ) == XLAL_SUCCESS, XLAL_EFUNC );
    LogPrintf( LOG_DEBUG, "Computed coherent reduced supersky metric for segment %zu/%zu\n", n, metrics->num_segments );

    // Cleanup
    GFMAT( ussky_metric_seg, orbital_metric_seg );

  }

  // Normalise averaged metrics by number of segments
  gsl_matrix_scale( ussky_metric_avg, 1.0 / metrics->num_segments );
  gsl_matrix_scale( orbital_metric_avg, 1.0 / metrics->num_segments );

  // Compute the semicoherent supersky metric for all segments
  const LIGOTimeGPS *start_time_avg = &segments->segs[0].start;
  const LIGOTimeGPS *end_time_avg = &segments->segs[segments->length - 1].end;
  XLAL_CHECK_NULL( SM_ComputeReducedSuperskyMetric( &metrics->semi_rssky_metric, &metrics->semi_rssky_transf, ussky_metric_avg, &ucoords, orbital_metric_avg, &ocoords, ref_time, start_time_avg, end_time_avg ) == XLAL_SUCCESS, XLAL_EFUNC );
  LogPrintf( LOG_DEBUG, "Computed semicoherent reduced supersky metric for %zu segments\n", metrics->num_segments );

  // Rescale metrics to input fiducial frequency
  XLALScaleSuperskyMetricsFiducialFreq( metrics, fiducial_freq );

  // Cleanup
  GFMAT( ussky_metric_avg, orbital_metric_avg );

  return metrics;

}

void XLALDestroySuperskyMetrics(
  SuperskyMetrics *metrics                      /// [in] Supersky metrics struct
  )
{
  if ( metrics != NULL ) {
    for ( size_t n = 0; n < metrics->num_segments; ++n ) {
      GFMAT( metrics->coh_rssky_metric[n], metrics->coh_rssky_transf[n] );
    }
    GFMAT( metrics->semi_rssky_metric, metrics->semi_rssky_transf );
    XLALFree( metrics->coh_rssky_metric );
    XLALFree( metrics->coh_rssky_transf );
    XLALFree( metrics );
  }
}

int XLALSuperskyMetricsDimensions(
  const SuperskyMetrics *metrics,
  size_t *spindowns
  )
{


  // Check input
  XLAL_CHECK( metrics != NULL, XLAL_EFAULT );
  XLAL_CHECK( metrics->num_segments > 0, XLAL_EINVAL );
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    XLAL_CHECK( CHECK_RSSKY_METRIC_TRANSF( metrics->coh_rssky_metric[n], metrics->coh_rssky_transf[n] ), XLAL_EINVAL );
  }
  XLAL_CHECK( CHECK_RSSKY_METRIC_TRANSF( metrics->semi_rssky_metric, metrics->semi_rssky_transf ), XLAL_EINVAL );

  // Decompose coordinate transform data
  DECOMPOSE_RSSKY_TRANSF( metrics->semi_rssky_transf );

  // Return dimensions
  if ( spindowns != NULL ) {
    *spindowns = smax;
  }

  return XLAL_SUCCESS;

}

int XLALScaleSuperskyMetricFiducialFreq(
  gsl_matrix *rssky_metric,
  gsl_matrix *rssky_transf,
  const double new_fiducial_freq
  )
{

  // Check input
  XLAL_CHECK( CHECK_RSSKY_METRIC_TRANSF( rssky_metric, rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( new_fiducial_freq > 0, XLAL_EINVAL );

  // Decompose coordinate transform data
  DECOMPOSE_RSSKY_TRANSF( rssky_transf );
  XLAL_CHECK( fiducial_freq > 0, XLAL_EINVAL );

  // Rescale metrics to 'new_fiducial_freq' based on known scalings
  const double fiducial_scale = new_fiducial_freq / fiducial_freq;
  gsl_matrix_view sky_sky = gsl_matrix_submatrix( rssky_metric, 0, 0, 2, 2 );
  gsl_matrix_scale( &sky_sky.matrix, SQR( fiducial_scale ) );
  gsl_matrix_scale( &sky_offsets.matrix, fiducial_scale );

  // Set new fiducial frequency
  gsl_matrix_set( rssky_transf, _RSSKY_TRANSF_FIDFREQ, new_fiducial_freq );

  return XLAL_SUCCESS;

}

int XLALScaleSuperskyMetricsFiducialFreq(
  SuperskyMetrics *metrics,
  const double new_fiducial_freq
  )
{

  // Check input
  XLAL_CHECK( metrics != NULL, XLAL_EFAULT );
  XLAL_CHECK( metrics->num_segments > 0, XLAL_EINVAL );
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    XLAL_CHECK( CHECK_RSSKY_METRIC_TRANSF( metrics->coh_rssky_metric[n], metrics->coh_rssky_transf[n] ), XLAL_EINVAL );
  }
  XLAL_CHECK( CHECK_RSSKY_METRIC_TRANSF( metrics->semi_rssky_metric, metrics->semi_rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( new_fiducial_freq > 0, XLAL_EINVAL );

  // Rescale all metrics to 'new_fiducial_freq'
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    XLAL_CHECK( XLALScaleSuperskyMetricFiducialFreq( metrics->coh_rssky_metric[n], metrics->coh_rssky_transf[n], new_fiducial_freq ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  XLAL_CHECK( XLALScaleSuperskyMetricFiducialFreq( metrics->semi_rssky_metric, metrics->semi_rssky_transf, new_fiducial_freq ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

int XLALEqualizeReducedSuperskyMetricsFreqSpacing(
  SuperskyMetrics *metrics,
  const double coh_max_mismatch,
  const double semi_max_mismatch
  )
{

  // Check input
  XLAL_CHECK( metrics != NULL, XLAL_EFAULT );
  XLAL_CHECK( metrics->num_segments > 0, XLAL_EINVAL );
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    XLAL_CHECK( CHECK_RSSKY_METRIC_TRANSF( metrics->coh_rssky_metric[n], metrics->coh_rssky_transf[n] ), XLAL_EINVAL );
  }
  XLAL_CHECK( CHECK_RSSKY_METRIC_TRANSF( metrics->semi_rssky_metric, metrics->semi_rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( coh_max_mismatch > 0, XLAL_EINVAL );
  XLAL_CHECK( semi_max_mismatch > 0, XLAL_EINVAL );

  // Decompose coordinate transform data
  DECOMPOSE_RSSKY_TRANSF( metrics->semi_rssky_transf );
  const size_t ifreq = RSSKY_FKDOT_DIM( 0 );

  // Find the maximum, over both coherent and semicoherent metrics, of 'g_{ff} / mu',
  // where 'g_{ff}' is the frequency-frequency metric element, and mu is the mismatch
  double max_rssky_metric_ff_d_mu = 0;
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    const double coh_rssky_metric_ff_d_mu = gsl_matrix_get( metrics->coh_rssky_metric[n], ifreq, ifreq ) / coh_max_mismatch;
    max_rssky_metric_ff_d_mu = GSL_MAX( max_rssky_metric_ff_d_mu, coh_rssky_metric_ff_d_mu );
  }
  {
    const double semi_rssky_metric_ff_d_mu = gsl_matrix_get( metrics->semi_rssky_metric, ifreq, ifreq ) / semi_max_mismatch;
    max_rssky_metric_ff_d_mu = GSL_MAX( max_rssky_metric_ff_d_mu, semi_rssky_metric_ff_d_mu );
  }

  // Project all metrics in the frequency dimension, and set frequency-frequency
  // metric element to 'max_rssky_metric_ff_d_mu' * 'mu', where mu is the mismatch
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    XLAL_CHECK( XLALProjectMetric( &metrics->coh_rssky_metric[n], metrics->coh_rssky_metric[n], ifreq ) == XLAL_SUCCESS, XLAL_EFUNC );
    gsl_matrix_set( metrics->coh_rssky_metric[n], ifreq, ifreq, max_rssky_metric_ff_d_mu * coh_max_mismatch );
  }
  {
    XLAL_CHECK( XLALProjectMetric( &metrics->semi_rssky_metric, metrics->semi_rssky_metric, ifreq ) == XLAL_SUCCESS, XLAL_EFUNC );
    gsl_matrix_set( metrics->semi_rssky_metric, ifreq, ifreq, max_rssky_metric_ff_d_mu * semi_max_mismatch );
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
  double as[3],                                 ///< [out] 3-dimensional aligned sky coordinates
  const gsl_vector *rss                         ///< [in] 2-dimensional reduced supersky coordinates
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
  gsl_vector *rss,                              ///< [out] 2-dimensional reduced supersky coordinates
  const double as[3]                            ///< [in] 3-dimensional aligned sky coordinates
  )
{
  const double r = sqrt( DOT3( as, as ) );
  const double A = GSL_SIGN( as[2] ) * ( ( as[0] / r ) + 1.0 );
  const double B = as[1] / r;
  gsl_vector_set( rss, 0, A );
  gsl_vector_set( rss, 1, B );
}

int XLALConvertPhysicalToSuperskyPoint(
  gsl_vector *out_rssky,
  const PulsarDopplerParams *in_phys,
  const gsl_matrix *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EFAULT );

  // Decompose coordinate transform data
  DECOMPOSE_CONST_RSSKY_TRANSF( rssky_transf );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK( out_rssky->size == ndim, XLAL_ESIZE );

  // Transform input physical point to reference time of coordinate transform data
  PulsarDopplerParams in_phys_ref = *in_phys;
  {
    const REAL8 dtau = XLALGPSDiff( ref_time, &in_phys_ref.refTime );
    XLAL_CHECK( XLALExtrapolatePulsarSpins( in_phys_ref.fkdot, in_phys_ref.fkdot, dtau ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Create array for intermediate coordinates
  double intm[4 + smax];
  gsl_vector_view intm_sky2 = gsl_vector_view_array( &intm[0], 2 );
  gsl_vector_view intm_sky3 = gsl_vector_view_array( &intm[0], 3 );
  gsl_vector_view intm_fspin = gsl_vector_view_array( &intm[3], 1 + smax );

  // Convert right ascension and declination to equatorial coordinates
  {
    const double cos_Delta = cos( in_phys_ref.Delta );
    intm[0] = cos( in_phys_ref.Alpha ) * cos_Delta;
    intm[1] = sin( in_phys_ref.Alpha ) * cos_Delta;
    intm[2] = sin( in_phys_ref.Delta );
  }

  // Copy frequency/spindowns to intermediate array; frequency goes last
  intm[3 + smax] = in_phys_ref.fkdot[0];
  for ( size_t s = 1; s <= smax; ++s ) {
    intm[2 + s] = in_phys_ref.fkdot[s];
  }

  // Apply the alignment transform to the supersky position to produced the aligned sky position:
  //   asky = align_sky * ssky
  double asky[3];
  gsl_vector_view asky_v = gsl_vector_view_array( asky, 3 );
  gsl_blas_dgemv( CblasNoTrans, 1.0, &align_sky.matrix, &intm_sky3.vector, 0.0, &asky_v.vector );

  // Add the inner product of the sky offsets with the aligned sky position
  // to the supersky spins and frequency to get the reduced supersky quantities:
  //   rssky_fspin[i] = ussky_fspin[i] + dot(sky_offsets[i], asky)
  gsl_blas_dgemv( CblasNoTrans, 1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &intm_fspin.vector );

  // Convert from 3-dimensional aligned sky coordinates to 2-dimensional reduced supersky coordinates
  SM_AlignedToReduced( &intm_sky3.vector, asky );

  // Copy intermediate array to output point
  {
    gsl_vector_view out_sky2 = gsl_vector_subvector( out_rssky, 0, 2 );
    gsl_vector_view out_fspin = gsl_vector_subvector( out_rssky, 2, 1 + smax );
    gsl_vector_memcpy( &out_sky2.vector, &intm_sky2.vector );
    gsl_vector_memcpy( &out_fspin.vector, &intm_fspin.vector );
  }

  return XLAL_SUCCESS;

}

int XLALConvertSuperskyToPhysicalPoint(
  PulsarDopplerParams *out_phys,
  const gsl_vector *in_rssky,
  const gsl_matrix *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EFAULT );

  // Decompose coordinate transform data
  DECOMPOSE_CONST_RSSKY_TRANSF( rssky_transf );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_rssky->size == ndim, XLAL_ESIZE );

  // Set output physical point reference time to that of of coordinate transform data
  out_phys->refTime = *ref_time;

  // Create array for intermediate coordinates
  double intm[4 + smax];
  gsl_vector_view intm_sky2 = gsl_vector_view_array( &intm[0], 2 );
  gsl_vector_view intm_sky3 = gsl_vector_view_array( &intm[0], 3 );
  gsl_vector_view intm_fspin = gsl_vector_view_array( &intm[3], 1 + smax );

  // Copy input point to intermediate array
  {
    gsl_vector_const_view in_sky2 = gsl_vector_const_subvector( in_rssky, 0, 2 );
    gsl_vector_const_view in_fspin = gsl_vector_const_subvector( in_rssky, 2, 1 + smax );
    gsl_vector_memcpy( &intm_sky2.vector, &in_sky2.vector );
    gsl_vector_memcpy( &intm_fspin.vector, &in_fspin.vector );
  }

  // Convert from 2-dimensional reduced supersky coordinates to 3-dimensional aligned sky coordinates
  double asky[3];
  SM_ReducedToAligned( asky, &intm_sky3.vector );
  gsl_vector_view asky_v = gsl_vector_view_array( asky, 3 );

  // Subtract the inner product of the sky offsets with the aligned sky position
  // from the reduced supersky spins and frequency to get the supersky quantities:
  //   ussky_fspin[i] = rssky_fspin[i] - dot(sky_offsets[i], asky)
  gsl_blas_dgemv( CblasNoTrans, -1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &intm_fspin.vector );

  // Apply the inverse alignment transform to the aligned sky position to produced the supersky position:
  //   ssky = align_sky^T * asky
  gsl_blas_dgemv( CblasTrans, 1.0, &align_sky.matrix, &asky_v.vector, 0.0, &intm_sky3.vector );

  // Copy frequency/spindowns to output physical point; frequency goes first
  out_phys->fkdot[0] = intm[3 + smax];
  for ( size_t s = 1; s <= smax; ++s ) {
    out_phys->fkdot[s] = intm[2 + s];
  }

  // Convert supersky position in equatorial coordinates to right ascension and declination
  out_phys->Alpha = atan2( intm[1], intm[0] );
  out_phys->Delta = atan2( intm[2], sqrt( SQR( intm[0] ) + SQR( intm[1] ) ) );
  XLALNormalizeSkyPosition( &out_phys->Alpha, &out_phys->Delta );

  return XLAL_SUCCESS;

}

int XLALConvertSuperskyToSuperskyPoint(
  gsl_vector *out_rssky,
  const gsl_matrix *out_rssky_transf,
  const gsl_vector *in_rssky,
  const gsl_matrix *in_rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( out_rssky_transf ), XLAL_EFAULT );
  XLAL_CHECK( in_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( in_rssky_transf ), XLAL_EFAULT );

  // Convert input reduced supersky point to physical coordinates
  PulsarDopplerParams XLAL_INIT_DECL( phys );
  XLAL_CHECK( XLALConvertSuperskyToPhysicalPoint( &phys, in_rssky, in_rssky_transf ) == XLAL_SUCCESS, XLAL_EINVAL );

  // Convert physical point to output reduced supersky coordinates
  XLAL_CHECK( XLALConvertPhysicalToSuperskyPoint( out_rssky, &phys, out_rssky_transf ) == XLAL_SUCCESS, XLAL_EINVAL );

  return XLAL_SUCCESS;

}

int XLALConvertPhysicalToSuperskyPoints(
  gsl_matrix **out_rssky,
  const gsl_matrix *in_phys,
  const gsl_matrix *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EFAULT );

  // Decompose coordinate transform data
  DECOMPOSE_CONST_RSSKY_TRANSF( rssky_transf );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_phys->size1 == ndim, XLAL_ESIZE );

  // Resize or allocate output matrix, if required
  if ( *out_rssky != NULL ) {
    if ( ( *out_rssky )->size1 != in_phys->size1 || ( *out_rssky )->size2 != in_phys->size2 ) {
      GFMAT( *out_rssky );
      *out_rssky = NULL;
    }
  }
  if ( *out_rssky == NULL ) {
    GAMAT( *out_rssky, in_phys->size1, in_phys->size2 );
  }

  // Loop over all input points
  for ( size_t j = 0; j < in_phys->size2; ++j ) {

    // Fill PulsarDopplerParams struct from input point
    PulsarDopplerParams XLAL_INIT_DECL( in_phys_j );
    in_phys_j.refTime = *ref_time;
    in_phys_j.Alpha = gsl_matrix_get( in_phys, 0, j );
    in_phys_j.Delta = gsl_matrix_get( in_phys, 1, j );
    for ( size_t s = 0; s <= smax; ++s ) {
      in_phys_j.fkdot[s] = gsl_matrix_get( in_phys, 2 + s, j );
    }

    // Convert point from physical to supersky coordinates
    gsl_vector_view out_rssky_j = gsl_matrix_column( *out_rssky, j );
    XLAL_CHECK( XLALConvertPhysicalToSuperskyPoint( &out_rssky_j.vector, &in_phys_j, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );

  }

  return XLAL_SUCCESS;

}

int XLALConvertSuperskyToPhysicalPoints(
  gsl_matrix **out_phys,
  const gsl_matrix *in_rssky,
  const gsl_matrix *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EFAULT );

  // Decompose coordinate transform data
  DECOMPOSE_CONST_RSSKY_TRANSF( rssky_transf );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_rssky->size1 == ndim, XLAL_ESIZE );

  // Resize or allocate output matrix, if required
  if ( *out_phys != NULL ) {
    if ( ( *out_phys )->size1 != in_rssky->size1 || ( *out_phys )->size2 != in_rssky->size2 ) {
      GFMAT( *out_phys );
      *out_phys = NULL;
    }
  }
  if ( *out_phys == NULL ) {
    GAMAT( *out_phys, in_rssky->size1, in_rssky->size2 );
  }

  // Loop over all input points
  for ( size_t j = 0; j < in_rssky->size2; ++j ) {

    // Convert point from supersky to physical coordinates
    gsl_vector_const_view in_rssky_j = gsl_matrix_const_column( in_rssky, j );
    PulsarDopplerParams XLAL_INIT_DECL( out_phys_j );
    XLAL_CHECK( XLALConvertSuperskyToPhysicalPoint( &out_phys_j, &in_rssky_j.vector, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Fill output point from PulsarDopplerParams struct
    gsl_matrix_set( *out_phys, 0, j, out_phys_j.Alpha );
    gsl_matrix_set( *out_phys, 1, j, out_phys_j.Delta );
    for ( size_t s = 0; s <= smax; ++s ) {
      gsl_matrix_set( *out_phys, 2 + s, j, out_phys_j.fkdot[s] );
    }

  }

  return XLAL_SUCCESS;

}

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
} PhysicalSkyBoundPiece;
typedef struct {
  double min_A;
  double min_A_bound;
  double max_A;
  double max_A_bound;
  PhysicalSkyBoundPiece pieces[6];
} PhysicalSkyBoundData;

static double PhysicalSkyBound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_vector *point
  )
{

  // Get bounds data
  const PhysicalSkyBoundData *psbd = ( const PhysicalSkyBoundData * ) data;

  // Decode the reduced supersky coordinates to get
  //   na = as[0] = Q_na . n
  const double A = gsl_vector_get( point, 0 );
  const double rssky[2] = { A, 0 };
  gsl_vector_const_view rssky_view = gsl_vector_const_view_array( rssky, 2 );
  double as[3];
  SM_ReducedToAligned( as, &rssky_view.vector );
  const double na = as[0];

  // Absolute limiting bound on 'nb = +/- sqrt(1 - na^2)'
  const double limit = RE_SQRT( 1 - SQR( na ) );

  // If 'A' is outside range '(min_A, max_A)', set bound to 'min_A_bound' or 'max_A_bound'
  double bound = GSL_NAN;
  if ( A <= psbd->min_A ) {
    bound = psbd->min_A_bound;
  } else if ( A >= psbd->max_A ) {
    bound = psbd->max_A_bound;
  } else {

    // Loop over bound pieces to find which one currently applies, based on 'max_A'
    for ( size_t i = 0; i < XLAL_NUM_ELEM( psbd->pieces ); ++i ) {
      const PhysicalSkyBoundPiece p = psbd->pieces[i];
      if ( A <= p.max_A ) {

        if ( p.type != 0 ) {

          // Set bound 'nb = +/- sqrt(1 - na^2)'
          bound = p.type * limit;

        } else {

          // Set bound 'nb' to either a constant right ascension or constant declination bound,
          // depending on bounds data set by XLALSetSuperskyPhysicalSkyBounds()
          const double c = ( na - p.na0 ) / p.r;
          double angle = asin( GSL_MAX( -1, GSL_MIN( c, 1 ) ) );
          if ( p.altroot ) {
            angle = LAL_PI - angle;
          }
          angle -= p.angle0;
          bound = p.C*cos( angle ) + p.S*sin( angle ) + p.Z;

        }

        break;

      }
    }

  }

  return GSL_MAX( -limit, GSL_MIN( bound, limit ) );

}

int XLALSetSuperskyPhysicalSkyBounds(
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
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_METRIC_TRANSF( rssky_metric, rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( gsl_matrix_get( rssky_metric, 0, 1 ) == 0, XLAL_EINVAL );
  XLAL_CHECK( gsl_matrix_get( rssky_metric, 1, 0 ) == 0, XLAL_EINVAL );

  // Decompose coordinate transform data
  DECOMPOSE_RSSKY_TRANSF( rssky_transf );

  // Check sky bound ranges
  XLAL_CHECK( ( fabs( alpha1 - alpha2 ) > 0 ) == ( fabs( delta1 - delta2 ) > 0 ), XLAL_EINVAL,
              "|alpha1 - alpha2| and |delta1 - delta2| must either be both zero, or both nonzero" );
  if ( fabs( alpha1 - alpha2 ) >= LAL_TWOPI ) {
    XLAL_CHECK( fabs( delta1 ) >= LAL_PI_2 && fabs( delta2 ) >= LAL_PI_2, XLAL_EINVAL,
                "If |alpha1 - alpha2| covers the whole sky, then |delta1| == |delta2| == PI/2 is required" );
  } else {
    XLAL_CHECK( fabs( alpha1 - alpha2 ) <= LAL_PI, XLAL_EINVAL,
                "If |alpha1 - alpha2| does not cover the whole sky, then |alpha1 - alpha2| <= PI is required" );
    XLAL_CHECK( -LAL_PI_2 <= delta1 && delta1 <= LAL_PI_2, XLAL_EINVAL,
                "If |alpha1 - alpha2| does not cover the whole sky, then -PI/2 <= delta1 <= PI/2 is required" );
    XLAL_CHECK( -LAL_PI_2 <= delta2 && delta2 <= LAL_PI_2, XLAL_EINVAL,
                "If |alpha1 - alpha2| does not cover the whole sky, then -PI/2 <= delta2 <= PI/2 is required" );
  }

  // If parameter space is a single point:
  if ( alpha1 == alpha2 && delta1 == delta2 ) {

    // Convert physical point to reduced supersky coordinates A and B
    double rssky_point[rssky_metric->size1];
    gsl_vector_view rssky_point_view = gsl_vector_view_array( rssky_point, rssky_metric->size1 );
    PulsarDopplerParams XLAL_INIT_DECL( phys_point );
    phys_point.Alpha = alpha1;
    phys_point.Delta = delta1;
    XLAL_CHECK( XLALConvertPhysicalToSuperskyPoint( &rssky_point_view.vector, &phys_point, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Set the parameter-space bounds on reduced supersky sky coordinates A and B
    for ( size_t dim = 0; dim < 2; ++dim ) {
      XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, dim, rssky_point[dim], rssky_point[dim] ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

    return XLAL_SUCCESS;

  }

  // Initialise bounds data
  PhysicalSkyBoundData XLAL_INIT_DECL( data_lower );
  PhysicalSkyBoundData XLAL_INIT_DECL( data_upper );
  for ( size_t i = 0; i < XLAL_NUM_ELEM( data_lower.pieces ); ++i ) {
    data_lower.pieces[i].max_A = data_upper.pieces[i].max_A = GSL_NEGINF;
  }

  // Special bounds data representing the lower/upper circular bounds on reduced supersky coordinate B
  const PhysicalSkyBoundPiece lower_circular = { .type = -1 };
  const PhysicalSkyBoundPiece upper_circular = { .type = 1 };

  // If parameter space is the entire sky:
  if ( fabs( alpha1 - alpha2 ) >= LAL_TWOPI ) {

    // Set bounds data to lower/upper circular bounds
    data_lower.min_A = data_upper.min_A = -2;
    data_lower.min_A_bound = data_upper.min_A_bound = 0;
    data_lower.max_A = data_upper.max_A = 2;
    data_lower.max_A_bound = data_upper.max_A_bound = 0;
    data_lower.pieces[0] = lower_circular;
    data_lower.pieces[0].max_A = GSL_POSINF;
    data_upper.pieces[0] = upper_circular;
    data_upper.pieces[0].max_A = GSL_POSINF;

    // Set the parameter-space bounds on reduced supersky sky coordinates A and B
    XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, 0, data_lower.min_A, data_lower.max_A ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALSetLatticeTilingBound( tiling, 1, PhysicalSkyBound, sizeof( data_lower ), &data_lower, &data_upper ) == XLAL_SUCCESS, XLAL_EFUNC );

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
    const double alpha = GSL_MIN( alpha1, alpha2 );
    const double cos_alpha = cos( alpha );
    const double sin_alpha = sin( alpha );
    const double Q_na_0 = gsl_matrix_get( &align_sky.matrix, 0, 0 );
    const double Q_na_1 = gsl_matrix_get( &align_sky.matrix, 0, 1 );
    const double Q_nb_0 = gsl_matrix_get( &align_sky.matrix, 1, 0 );
    const double Q_nb_1 = gsl_matrix_get( &align_sky.matrix, 1, 1 );
    const double na = Q_na_0*cos_alpha + Q_na_1*sin_alpha;
    const double nb = Q_nb_0*cos_alpha + Q_nb_1*sin_alpha;
    phi = atan2( -nb, na );
    const double na_rot = na*cos( phi ) + nb*sin( phi );
    if ( na_rot < 0 ) {
      phi -= LAL_PI;
    } else if ( na_rot > 0 ) {
      phi += LAL_PI;
    }
  }
  const double cos_phi = cos( phi );
  const double sin_phi = sin( phi );

  // Apply the right-handed rotation matrix
  //   R = [cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1]
  // to the reduced supersky coordinate transform data
  //   rssky_transf = [Q^T; Delta^s]
  // where 'Q' is the sky alignment matrix and 'Delta^s' are the sky offset vectors.
  // The correct transformation to apply is:
  //   Q^T ==> R * Q^T, Delta^s ==> Delta^s . R^T
  for ( size_t j = 0; j < 3; ++j ) {
    const double Q_na_j = gsl_matrix_get( &align_sky.matrix, 0, j );
    const double Q_nb_j = gsl_matrix_get( &align_sky.matrix, 1, j );
    const double Q_na_j_rot = Q_na_j*cos_phi - Q_nb_j*sin_phi;
    const double Q_nb_j_rot = Q_nb_j*cos_phi + Q_na_j*sin_phi;
    gsl_matrix_set( &align_sky.matrix, 0, j, Q_na_j_rot );
    gsl_matrix_set( &align_sky.matrix, 1, j, Q_nb_j_rot );
  }
  for ( size_t i = 0; i < sky_offsets.matrix.size1; ++i ) {
    const double Delta_0 = gsl_matrix_get( &sky_offsets.matrix, i, 0 );
    const double Delta_1 = gsl_matrix_get( &sky_offsets.matrix, i, 1 );
    const double Delta_0_rot = Delta_0*cos_phi - Delta_1*sin_phi;
    const double Delta_1_rot = Delta_1*cos_phi + Delta_0*sin_phi;
    gsl_matrix_set( &sky_offsets.matrix, i, 0, Delta_0_rot );
    gsl_matrix_set( &sky_offsets.matrix, i, 1, Delta_1_rot );
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
    const double g_na_na = gsl_matrix_get( rssky_metric, 0, 0 );
    const double g_nb_nb = gsl_matrix_get( rssky_metric, 1, 1 );
    const double g_na_na_rot = g_na_na*SQR( cos_phi ) + g_nb_nb*SQR( sin_phi );
    const double g_na_nb_rot = ( g_na_na - g_nb_nb )*cos_phi*sin_phi;
    const double g_nb_nb_rot = g_na_na*SQR( sin_phi ) + g_nb_nb*SQR( cos_phi );
    gsl_matrix_set( rssky_metric, 0, 0, g_na_na_rot );
    gsl_matrix_set( rssky_metric, 0, 1, g_na_nb_rot );
    gsl_matrix_set( rssky_metric, 1, 0, g_na_nb_rot );
    gsl_matrix_set( rssky_metric, 1, 1, g_nb_nb_rot );
  }

  // Get components of the vectors 'Q_na', 'Q_nb', and 'Q_nc' from coordinate transform data
  const double Q_na[3] = { gsl_matrix_get( &align_sky.matrix, 0, 0 ), gsl_matrix_get( &align_sky.matrix, 0, 1 ), gsl_matrix_get( &align_sky.matrix, 0, 2 ) };
  const double Q_nb[3] = { gsl_matrix_get( &align_sky.matrix, 1, 0 ), gsl_matrix_get( &align_sky.matrix, 1, 1 ), gsl_matrix_get( &align_sky.matrix, 1, 2 ) };
  const double Q_nc[3] = { gsl_matrix_get( &align_sky.matrix, 2, 0 ), gsl_matrix_get( &align_sky.matrix, 2, 1 ), gsl_matrix_get( &align_sky.matrix, 2, 2 ) };

  // Determine the minimum and maximum right ascension and declination
  const double alphas[2] = { GSL_MIN( alpha1, alpha2 ), GSL_MAX( alpha1, alpha2 ) };
  const double deltas[2] = { GSL_MIN( delta1, delta2 ), GSL_MAX( delta1, delta2 ) };

  // Create bound data for declination bounds, at constant minimum/maximum right ascension:
  // Given known 'na' from previous bound, and known minimum/maximum 'alpha', solve
  //   na = ( Q_na[0]*cos(alpha) + Q_na[1]*sin(alpha) )*cos(delta) + Q_na[2]*sin(delta)
  // for
  //   delta = asin( ( na - na0 ) / r ) - angle0
  // and compute
  //   nb = C*cos(delta) + S*sin(delta) + Z
  PhysicalSkyBoundPiece const_alpha[2];
  for ( size_t i = 0; i < 2; ++i ) {
    XLAL_INIT_MEM( const_alpha[i] );
    const_alpha[i].type = 0;
    const double cos_alpha = cos( alphas[i] );
    const double sin_alpha = sin( alphas[i] );
    const double x = Q_na[2];
    const double y = Q_na[0]*cos_alpha + Q_na[1]*sin_alpha;
    const_alpha[i].na0 = 0;
    const_alpha[i].r = sqrt( SQR( x ) + SQR( y ) );
    const_alpha[i].angle0 = atan2( y, x );
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
  PhysicalSkyBoundPiece const_delta[2];
  for ( size_t j = 0; j < 2; ++j ) {
    XLAL_INIT_MEM( const_delta[j] );
    const_delta[j].type = 0;
    const double cos_delta = cos( deltas[j] );
    const double sin_delta = sin( deltas[j] );
    const double x = Q_na[1]*cos_delta;
    const double y = Q_na[0]*cos_delta;
    const_delta[j].na0 = Q_na[2]*sin_delta;
    const_delta[j].r = sqrt( SQR( x ) + SQR( y ) );
    const_delta[j].angle0 = atan2( y, x );
    const_delta[j].C = Q_nb[0]*cos_delta;
    const_delta[j].S = Q_nb[1]*cos_delta;
    const_delta[j].Z = Q_nb[2]*sin_delta;
  }

  // Determine corner points in reduced supersky coordinate A of the region
  //   '[min(alpha),max(alpha)] x [min(delta),max(delta)]'
  double corner_A[2][2], corner_B[2][2];
  for ( size_t i = 0; i < 2; ++i ) {
    for ( size_t j = 0; j < 2; ++j ) {
      double rssky_point[rssky_metric->size1];
      gsl_vector_view rssky_point_view = gsl_vector_view_array( rssky_point, rssky_metric->size1 );
      PulsarDopplerParams XLAL_INIT_DECL( phys_point );
      phys_point.Alpha = alphas[i];
      phys_point.Delta = deltas[j];
      XLAL_CHECK( XLALConvertPhysicalToSuperskyPoint( &rssky_point_view.vector, &phys_point, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
      corner_A[i][j] = rssky_point[0];
      corner_B[i][j] = rssky_point[1];
    }
  }

  // Use corner points to classify parameter space into different shapes and set bounds data
  data_lower.min_A = data_upper.min_A = GSL_NEGINF;
  data_lower.max_A = data_upper.max_A = GSL_POSINF;
  if ( corner_A[1][0] < 0 && corner_A[1][1] <= 0 ) {

    if ( corner_A[1][1] < corner_A[1][0] ) {

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
      data_lower.min_A = data_upper.min_A = corner_A[1][1];
      data_lower.min_A_bound = data_upper.min_A_bound = corner_B[1][1];
      data_lower.max_A = data_upper.max_A = corner_A[0][1];
      data_lower.max_A_bound = data_upper.max_A_bound = corner_B[0][1];
      data_lower.pieces[0] = const_delta[1];
      data_lower.pieces[0].max_A = GSL_POSINF;
      data_upper.pieces[0] = const_alpha[1];
      data_upper.pieces[0].max_A = corner_A[1][0];
      data_upper.pieces[1] = const_delta[0];
      data_upper.pieces[1].max_A = corner_A[0][0];
      data_upper.pieces[2] = const_alpha[0];
      data_upper.pieces[2].max_A = GSL_POSINF;
      data_upper.pieces[2].altroot = true;

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
      data_lower.min_A = data_upper.min_A = corner_A[1][0];
      data_lower.min_A_bound = data_upper.min_A_bound = corner_B[1][0];
      data_lower.max_A = data_upper.max_A = corner_A[0][1];
      data_lower.max_A_bound = data_upper.max_A_bound = corner_B[0][1];
      data_lower.pieces[0] = const_alpha[1];
      data_lower.pieces[0].max_A = corner_A[1][1];
      data_lower.pieces[0].altroot = true;
      data_lower.pieces[1] = const_delta[1];
      data_lower.pieces[1].max_A = GSL_POSINF;
      data_upper.pieces[0] = const_delta[0];
      data_upper.pieces[0].max_A = corner_A[0][0];
      data_upper.pieces[1] = const_alpha[0];
      data_upper.pieces[1].max_A = GSL_POSINF;
      data_upper.pieces[1].altroot = true;

    }

  } else if ( 0 <= corner_A[1][0] && 0 < corner_A[1][1] ) {

    if ( corner_A[1][1] < corner_A[1][0] ) {

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
      data_lower.min_A = data_upper.min_A = corner_A[0][0];
      data_lower.min_A_bound = data_upper.min_A_bound = corner_B[0][0];
      data_lower.max_A = data_upper.max_A = corner_A[1][0];
      data_lower.max_A_bound = data_upper.max_A_bound = corner_B[1][0];
      data_lower.pieces[0] = const_delta[0];
      data_lower.pieces[0].max_A = GSL_POSINF;
      data_upper.pieces[0] = const_alpha[0];
      data_upper.pieces[0].max_A = corner_A[0][1];
      data_upper.pieces[1] = const_delta[1];
      data_upper.pieces[1].max_A = corner_A[1][1];
      data_upper.pieces[2] = const_alpha[1];
      data_upper.pieces[2].max_A = GSL_POSINF;
      data_upper.pieces[2].altroot = true;

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
      data_lower.min_A = data_upper.min_A = corner_A[0][0];
      data_lower.min_A_bound = data_upper.min_A_bound = corner_B[0][0];
      data_lower.max_A = data_upper.max_A = corner_A[1][1];
      data_lower.max_A_bound = data_upper.max_A_bound = corner_B[1][1];
      data_lower.pieces[0] = const_delta[0];
      data_lower.pieces[0].max_A = corner_A[1][0];
      data_lower.pieces[1] = const_alpha[1];
      data_lower.pieces[1].max_A = GSL_POSINF;
      data_upper.pieces[0] = const_alpha[0];
      data_upper.pieces[0].max_A = corner_A[0][1];
      data_upper.pieces[1] = const_delta[1];
      data_upper.pieces[1].max_A = GSL_POSINF;

    }

  } else {

    // This parameter space straddles both reduced supersky hemispheres
    // Find the value of 'na' where this occurs at max(alpha) by solving
    //   nc = ( Q_nc[0]*cos(max(alpha)) + Q_nc[1]*sin(max(alpha)) )*cos(delta) + Q_nc[2]*sin(delta)
    // for delta, then computing
    //   na = ( Q_na[0]*cos(alpha) + Q_na[1]*sin(alpha) )*cos(delta) + Q_na[2]*sin(delta)
    const double cos_alpha_split = cos( alphas[1] );
    const double sin_alpha_split = sin( alphas[1] );
    const double delta_split = -1 * atan2( Q_nc[0]*cos_alpha_split + Q_nc[1]*sin_alpha_split, Q_nc[2] );
    const double na_split = ( Q_na[0]*cos_alpha_split + Q_na[1]*sin_alpha_split )*cos( delta_split ) + Q_na[2]*sin( delta_split );
    const double split_A[2] = { -1 - na_split, 1 + na_split };
    const double split_B = -1 * RE_SQRT( 1 - SQR( na_split ) );

    if ( split_A[0] < corner_A[1][0] ) {
      if ( corner_A[1][1] < split_A[1] ) {

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
        data_lower.min_A = data_upper.min_A = split_A[0];
        data_lower.min_A_bound = data_upper.min_A_bound = split_B;
        data_lower.max_A = data_upper.max_A = split_A[1];
        data_lower.max_A_bound = data_upper.max_A_bound = split_B;
        data_lower.pieces[0] = lower_circular;
        data_lower.pieces[0].max_A = GSL_POSINF;
        data_upper.pieces[0] = const_alpha[1];
        data_upper.pieces[0].max_A = corner_A[1][0];
        data_upper.pieces[1] = const_delta[0];
        data_upper.pieces[1].max_A = corner_A[0][0];
        data_upper.pieces[2] = const_alpha[0];
        data_upper.pieces[2].max_A = 0;
        data_upper.pieces[2].altroot = true;
        data_upper.pieces[3] = const_alpha[0];
        data_upper.pieces[3].max_A = corner_A[0][1];
        data_upper.pieces[4] = const_delta[1];
        data_upper.pieces[4].max_A = corner_A[1][1];
        data_upper.pieces[5] = const_alpha[1];
        data_upper.pieces[5].max_A = GSL_POSINF;
        data_upper.pieces[5].altroot = true;

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
        data_lower.min_A = data_upper.min_A = split_A[0];
        data_lower.min_A_bound = data_upper.min_A_bound = split_B;
        data_lower.max_A = data_upper.max_A = corner_A[1][1];
        data_lower.max_A_bound = data_upper.max_A_bound = corner_B[1][1];
        data_lower.pieces[0] = lower_circular;
        data_lower.pieces[0].max_A = split_A[1];
        data_lower.pieces[1] = const_alpha[1];
        data_lower.pieces[1].max_A = GSL_POSINF;
        data_upper.pieces[0] = const_alpha[1];
        data_upper.pieces[0].max_A = corner_A[1][0];
        data_upper.pieces[1] = const_delta[0];
        data_upper.pieces[1].max_A = corner_A[0][0];
        data_upper.pieces[2] = const_alpha[0];
        data_upper.pieces[2].max_A = 0;
        data_upper.pieces[2].altroot = true;
        data_upper.pieces[3] = const_alpha[0];
        data_upper.pieces[3].max_A = corner_A[0][1];
        data_upper.pieces[4] = const_delta[1];
        data_upper.pieces[4].max_A = GSL_POSINF;

      }

    } else {
      if ( corner_A[1][1] < split_A[1] ) {

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
        data_lower.min_A = data_upper.min_A = corner_A[1][0];
        data_lower.min_A_bound = data_upper.min_A_bound = corner_B[1][0];
        data_lower.max_A = data_upper.max_A = split_A[1];
        data_lower.max_A_bound = data_upper.max_A_bound = split_B;
        data_lower.pieces[0] = const_alpha[1];
        data_lower.pieces[0].max_A = split_A[0];
        data_lower.pieces[0].altroot = true;
        data_lower.pieces[1] = lower_circular;
        data_lower.pieces[1].max_A = GSL_POSINF;
        data_upper.pieces[0] = const_delta[0];
        data_upper.pieces[0].max_A = corner_A[0][0];
        data_upper.pieces[1] = const_alpha[0];
        data_upper.pieces[1].max_A = 0;
        data_upper.pieces[1].altroot = true;
        data_upper.pieces[2] = const_alpha[0];
        data_upper.pieces[2].max_A = corner_A[0][1];
        data_upper.pieces[3] = const_delta[1];
        data_upper.pieces[3].max_A = corner_A[1][1];
        data_upper.pieces[4] = const_alpha[1];
        data_upper.pieces[4].max_A = GSL_POSINF;
        data_upper.pieces[4].altroot = true;

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
        data_lower.min_A = data_upper.min_A = corner_A[1][0];
        data_lower.min_A_bound = data_upper.min_A_bound = corner_B[1][0];
        data_lower.max_A = data_upper.max_A = corner_A[1][1];
        data_lower.max_A_bound = data_upper.max_A_bound = corner_B[1][1];
        data_lower.pieces[0] = const_alpha[1];
        data_lower.pieces[0].max_A = split_A[0];
        data_lower.pieces[0].altroot = true;
        data_lower.pieces[1] = lower_circular;
        data_lower.pieces[1].max_A = split_A[1];
        data_lower.pieces[2] = const_alpha[1];
        data_lower.pieces[2].max_A = GSL_POSINF;
        data_upper.pieces[0] = const_delta[0];
        data_upper.pieces[0].max_A = corner_A[0][0];
        data_upper.pieces[1] = const_alpha[0];
        data_upper.pieces[1].max_A = 0;
        data_upper.pieces[1].altroot = true;
        data_upper.pieces[2] = const_alpha[0];
        data_upper.pieces[2].max_A = corner_A[0][1];
        data_upper.pieces[3] = const_delta[1];
        data_upper.pieces[3].max_A = GSL_POSINF;

      }
    }

  }

  // Set the parameter-space bounds on reduced supersky sky coordinates A and B
  XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, 0, data_lower.min_A, data_lower.max_A ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingBound( tiling, 1, PhysicalSkyBound, sizeof( data_lower ), &data_lower, &data_upper ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

int XLALComputePhysicalSkyEqualAreaPatch(
  double *alpha1,
  double *alpha2,
  double *delta1,
  double *delta2,
  const UINT4 patch_count,
  const UINT4 patch_index
  )
{

  // Check input
  XLAL_CHECK( alpha1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( alpha2 != NULL, XLAL_EFAULT );
  XLAL_CHECK( delta1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( delta2 != NULL, XLAL_EFAULT );
  XLAL_CHECK( patch_count > 0, XLAL_EINVAL );
  XLAL_CHECK( patch_index < patch_count, XLAL_EINVAL );

  // Number of patch divisions in 'alpha'; for less than 4 patches, divide only in 'alpha' to prevent
  // 'alpha' range in [pi,2*pi], which XLALSetSuperskyPhysicalSkyBounds() cannot handle
  const UINT4 alpha_count = ( patch_count < 4 ) ? patch_count : ( ( UINT4 ) ceil( sqrt( patch_count ) ) );

  // Mininum number of patch divisions in 'sin(delta)'; note integer division equivalent to floor()
  const UINT4 min_sdelta_count = patch_count / alpha_count;

  // Excess number of patches, which must be added on to get 'patch_count'
  INT4 patch_excess = patch_count - alpha_count * min_sdelta_count;
  XLAL_CHECK( patch_excess >= 0, XLAL_EFAILED );

  // Initialise number of patch divisions in 'sin(delta)'; if there are excess patches, add an extra patch
  UINT4 sdelta_count = min_sdelta_count;
  if ( patch_excess > 0 ) {
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
  while ( sdelta_index >= sdelta_count ) {

    // Decrease index in 'sin(delta)'; we are done when 'sdelta_index' < 'sdelta_count'
    sdelta_index -= sdelta_count;

    // Decrease number of excess patches; if zero, subtract extra patch from patch divisions in 'sin(delta)'
    --patch_excess;
    if ( patch_excess == 0 ) {
      --sdelta_count;
    }

    // Store the current last 'alpha' index in 'alpha_index1', and increase
    // 'alpha_index2' by the current number of patch divisions in 'sin(delta)'
    alpha_index1 = alpha_index2;
    alpha_index2 += sdelta_count;

  }

  // Compute range of 'alpha' to bound
  *alpha1 = LAL_TWOPI * ( ( double ) alpha_index1 ) / ( ( double ) patch_count );
  *alpha2 = LAL_TWOPI * ( ( double ) alpha_index2 ) / ( ( double ) patch_count );

  // Compute range of 'delta' to bound
  *delta1 = asin( -1 + 2 * ( ( double ) sdelta_index ) / ( ( double ) sdelta_count ) );
  *delta2 = asin( -1 + 2 * ( ( double ) sdelta_index + 1 ) / ( ( double ) sdelta_count ) );

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

int XLALSetSuperskyPhysicalSpinBound(
  LatticeTiling *tiling,
  const gsl_matrix *rssky_transf,
  const size_t s,
  const double bound1,
  const double bound2
  )
{

  // Check input
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( isfinite( bound1 ), XLAL_EINVAL );
  XLAL_CHECK( isfinite( bound2 ), XLAL_EINVAL );

  // Decompose coordinate transform data
  DECOMPOSE_CONST_RSSKY_TRANSF( rssky_transf );
  XLAL_CHECK( s <= smax, XLAL_ESIZE );

  // Copy the sky offset vector to bounds data
  double data_lower[4], data_upper[4];
  for ( size_t j = 0; j < 3; ++j ) {
    data_lower[j] = data_upper[j] = gsl_matrix_get( &sky_offsets.matrix, RSSKY_FKDOT_OFFSET( s ), j );
  }
  data_lower[3] = GSL_MIN( bound1, bound2 );
  data_upper[3] = GSL_MAX( bound1, bound2 );

  // Set the parameter-space bound on physical frequency/spindown coordinate
  XLAL_CHECK( XLALSetLatticeTilingBound( tiling, RSSKY_FKDOT_DIM( s ), PhysicalSpinBound, sizeof( data_lower ), &data_lower, &data_upper ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

int XLALSetSuperskyCoordinateSpinBound(
  LatticeTiling *tiling,
  const gsl_matrix *rssky_transf,
  const size_t s,
  const double bound1,
  const double bound2
  )
{

  // Check input
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( isfinite( bound1 ), XLAL_EINVAL );
  XLAL_CHECK( isfinite( bound2 ), XLAL_EINVAL );

  // Decompose coordinate transform data
  DECOMPOSE_CONST_RSSKY_TRANSF( rssky_transf );
  XLAL_CHECK( s <= smax, XLAL_ESIZE );

  // Set the parameter-space bound on reduced supersky frequency/spindown coordinate
  XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, RSSKY_FKDOT_DIM( s ), bound1, bound2 ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

int XLALSuperskyLatticePulsarSpinRange(
  PulsarSpinRange *spin_range,
  LatticeTiling *tiling,
  const gsl_matrix *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( spin_range != NULL, XLAL_EFAULT );
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );

  // Decompose coordinate transform data
  DECOMPOSE_CONST_RSSKY_TRANSF( rssky_transf );
  XLAL_CHECK( ref_time != NULL, XLAL_EFAULT );

  // Create arrays for minimum/maximum physical frequency/spindowns
  PulsarSpins fkdotMin, fkdotMax;
  for ( size_t s = 0; s <= smax; ++s ) {
    fkdotMin[s] = GSL_POSINF;
    fkdotMax[s] = GSL_NEGINF;
  }

  // Get frequency step size
  const double dfreq = XLALLatticeTilingStepSizes( tiling, ndim - 1 );

  // Create iterator over reduced supersky coordinates
  LatticeTilingIterator *itr = XLALCreateLatticeTilingIterator( tiling, ndim - 1 );
  XLAL_CHECK( itr != NULL, XLAL_EFUNC );

  // Iterate over reduced supersky coordinates
  double in_rssky_array[ndim];
  gsl_vector_view in_rssky_view = gsl_vector_view_array( in_rssky_array, ndim );
  gsl_vector *const in_rssky = &in_rssky_view.vector;
  PulsarDopplerParams XLAL_INIT_DECL( out_phys );
  while ( XLALNextLatticeTilingPoint( itr, in_rssky ) > 0 ) {

    // Convert reduced supersky point to physical coordinates
    XLAL_CHECK( XLALConvertSuperskyToPhysicalPoint( &out_phys, in_rssky, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Get indexes of left/right-most point in current frequency block
    INT4 left = 0, right = 0;
    XLAL_CHECK( XLALCurrentLatticeTilingBlock( itr, ndim - 1, &left, &right ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Store minimum/maximum physical frequency
    fkdotMin[0] = GSL_MIN( fkdotMin[0], out_phys.fkdot[0] + dfreq * left );
    fkdotMax[0] = GSL_MAX( fkdotMax[0], out_phys.fkdot[0] + dfreq * right );

    // Store minimum/maximum physical spindowns
    for ( size_t s = 0; s <= smax; ++s ) {
      fkdotMin[s] = GSL_MIN( fkdotMin[s], out_phys.fkdot[s] );
      fkdotMax[s] = GSL_MAX( fkdotMax[s], out_phys.fkdot[s] );
    }

  }

  // Initialise 'spin_range' to zero
  XLAL_INIT_MEM( *spin_range );

  // Set reference time of 'spin_range' to that of coordinate transform data
  spin_range->refTime = *ref_time;

  // Set spindown range
  for ( size_t s = 0; s <= smax; ++s ) {
    spin_range->fkdot[s] = fkdotMin[s];
    spin_range->fkdotBand[s] = fkdotMax[s] - fkdotMin[s];
  }

  // Cleanup
  XLALDestroyLatticeTilingIterator( itr );

  return XLAL_SUCCESS;

}
