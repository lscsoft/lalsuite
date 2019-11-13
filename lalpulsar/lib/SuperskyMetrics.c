//
// Copyright (C) 2014--2017 Karl Wette
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

// Maximum number of sky offsets required
#define MAX_SKY_OFFSETS PULSAR_MAX_SPINS

// FIXME: replace 'SMAX' with either 'nspins', 'nsky_offsets - 1', or something else...
#define SMAX nspins

// Internal definition of reduced supersky metric coordinate transform data
struct tagSuperskyTransformData {
  UINT4 ndim;                                   ///< Dimensions of the corresponding metric
  LIGOTimeGPS ref_time;                         ///< Reference time of the metric
  REAL8 fiducial_freq;                          ///< Fiducial frequency of metric
  UINT4 nspins;                                 ///< Number of spindown dimensions
  REAL8 align_sky[3][3];                        ///< Alignment transform of the supersky metric
  UINT4 nsky_offsets;                           ///< Number of sky offsets
  REAL8 sky_offsets[MAX_SKY_OFFSETS][3];        ///< Sky offsets of the supersky metric
};

// Check reduced supersky coordinate metric and/or transform data
#define CHECK_RSSKY_METRIC_TRANSF(M, RT) \
  ((M) != NULL && CHECK_RSSKY_TRANSF(RT) && (M)->size1 == (RT)->ndim && (M)->size2 == (RT)->ndim)
#define CHECK_RSSKY_TRANSF(RT) \
  ((RT) != NULL && (RT)->ndim > 0 && (RT)->fiducial_freq > 0 && (RT)->nsky_offsets == 1 + (RT)->nspins)

// Determine which dimension stores the reduced supersky frequency/spindown of order 's'
#define RSSKY_FKDOT_OFFSET(RT, S)   (((S) == 0) ? ((RT)->SMAX) : ((size_t)((S) - 1)))
#define RSSKY_FKDOT_DIM(RT, S)      (2 + RSSKY_FKDOT_OFFSET(RT, S))

///
/// Fiducial frequency at which to numerically calculate metrics, which
/// are then rescaled to user-requested frequency based on known scalings
///
const double fiducial_calc_freq = 100.0;

// Structures for callback functions
typedef struct tagSM_CallbackParam {
  const SuperskyTransformData *rssky_transf;    ///< Reduced supersky coordinate transform data
  const SuperskyTransformData *rssky2_transf;   ///< Other reduced supersky coordinate transform data
} SM_CallbackParam;
typedef struct tagSM_CallbackOut {
  PulsarDopplerParams min_phys;                 ///< Minimum physical range
  PulsarDopplerParams max_phys;                 ///< Maximum physical range
  double min_rssky2_array[32];                  ///< Minimum range of other reduced supersky coordinates
  gsl_vector_view min_rssky2_view;
  double max_rssky2_array[32];                  ///< Maximum range of other reduced supersky coordinates
  gsl_vector_view max_rssky2_view;
} SM_CallbackOut;

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

  // Supersky metric cannot (reliably) be computed for segment lengths <= ~24 hours
  XLAL_CHECK_NULL( XLALGPSDiff( end_time, start_time ) >= 81000, XLAL_ERANGE, "Supersky metric cannot be computed for segment lengths <= ~24 hours" );

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

  // Set reference time to segment mid-time, to improve stability of numerical calculation of metric
  par.signalParams.Doppler.refTime = *start_time;
  XLALGPSAdd( &par.signalParams.Doppler.refTime, 0.5 * XLALGPSDiff( end_time, start_time ) );

  // Set fiducial frequency
  par.signalParams.Doppler.fkdot[0] = fiducial_calc_freq;

  // Do not project metric
  par.projectCoord = -1;

  // Do not include sky-position-dependent Roemer delay in time variable
  par.approxPhase = 1;

  // Call XLALComputeDopplerPhaseMetric() and check output
  DopplerPhaseMetric *metric = XLALComputeDopplerPhaseMetric( &par, ephemerides );
  XLAL_CHECK_NULL( metric != NULL && metric->g_ij != NULL, XLAL_EFUNC, "XLALComputeDopplerPhaseMetric() failed" );

  // Extract metric, while transforming its reference time from segment mid-time to time specified by 'ref_time'
  gsl_matrix *g_ij = NULL;
  const REAL8 Dtau = XLALGPSDiff( ref_time, &par.signalParams.Doppler.refTime );
  XLAL_CHECK_NULL( XLALChangeMetricReferenceTime( &g_ij, NULL, metric->g_ij, coords, Dtau ) == XLAL_SUCCESS, XLAL_EFUNC );

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
  SuperskyTransformData *rssky_transf,          ///< [in,out] Reduced supersky metric coordinate transform data
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

  // Allocate memory
  gsl_matrix *GAMAT( tmp, 3 + rssky_transf->SMAX, 3 + rssky_transf->SMAX );
  gsl_vector *GAVEC( tmpv, 1 + rssky_transf->SMAX );

  // Compute mid-time of segment list
  LIGOTimeGPS mid_time = *start_time;
  XLALGPSAdd( &mid_time, 0.5 * XLALGPSDiff( end_time, start_time ) );

  // Internal copy of orbital metric, and various transforms performed on it
  gsl_matrix *orb_metric = NULL, *mid_time_transf = NULL, *diag_norm_transf = NULL;

  // Transform reference time of orbital metric from reference time to segment list mid-time
  const REAL8 Dtau = XLALGPSDiff( &mid_time, &rssky_transf->ref_time );
  XLAL_CHECK( XLALChangeMetricReferenceTime( &orb_metric, &mid_time_transf, orbital_metric, ocoords, Dtau ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Diagonally-normalize orbital metric
  XLAL_CHECK( XLALDiagNormalizeMetric( &orb_metric, &diag_norm_transf, orb_metric ) == XLAL_SUCCESS, XLAL_EFUNC );

  // 'fitA' contains the frequency and spindown elements of the orbital metric, used for fitting
  gsl_matrix *GAMAT( fitA, 3 + rssky_transf->SMAX, 1 + rssky_transf->SMAX );
  {
    gsl_matrix_view orb_metric_fspin = gsl_matrix_submatrix( orb_metric, 0, 2, 3 + rssky_transf->SMAX, 1 + rssky_transf->SMAX );
    gsl_matrix_memcpy( fitA, &orb_metric_fspin.matrix );
  }

  // Compute 'fitA^T * fitA'
  gsl_matrix *GAMAT( fitAt_fitA, 1 + rssky_transf->SMAX, 1 + rssky_transf->SMAX );
  gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, fitA, fitA, 0.0, fitAt_fitA );

  // Find the singular value decomposition of 'fitA^T * fitA'
  gsl_matrix *GAMAT( svd_U, 1 + rssky_transf->SMAX, 1 + rssky_transf->SMAX );
  gsl_matrix *GAMAT( svd_V, 1 + rssky_transf->SMAX, 1 + rssky_transf->SMAX );
  gsl_vector *GAVEC( svd_S, 1 + rssky_transf->SMAX );
  gsl_matrix_memcpy( svd_U, fitAt_fitA );
  XLAL_CHECK( gsl_linalg_SV_decomp( svd_U, svd_V, svd_S, tmpv ) == 0, XLAL_EFAILED );

  // The columns of 'fitc' contain the least-square fitting coefficients for the orbital X and Y metric elements:
  //    fitc(:,j) = inv(fitA^T * fitA) * fitA^T * orb_metric(:,j)
  // The singular decomposition of fitA^T * fitA is used for the inverse
  gsl_matrix *GAMAT( fitc, 1 + rssky_transf->SMAX, 2 );
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
  gsl_matrix *GAMAT( subtract_orb, 3 + rssky_transf->SMAX, 3 + rssky_transf->SMAX );
  {
    gsl_matrix_set_identity( subtract_orb );
    gsl_matrix_view subtract_orb_fspin_sky = gsl_matrix_submatrix( subtract_orb, 2, 0, 1 + rssky_transf->SMAX, 2 );
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
  gsl_matrix *GAMAT( subtract_ussky, 4 + rssky_transf->SMAX, 4 + rssky_transf->SMAX );
  {
    gsl_matrix_set_identity( subtract_ussky );
    gsl_matrix_view subtract_ussky_fspin_sky = gsl_matrix_submatrix( subtract_ussky, 3, 0, 1 + rssky_transf->SMAX, 3 );
    gsl_matrix_view subtract_orb_fspin_sky = gsl_matrix_submatrix( subtract_orb, 2, 0, 1 + rssky_transf->SMAX, 2 );
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
    gsl_matrix_view sky_offsets = gsl_matrix_view_array( &rssky_transf->sky_offsets[0][0], rssky_transf->nsky_offsets, 3 );
    gsl_matrix_view subtract_ussky_fspin_sky = gsl_matrix_submatrix( subtract_ussky, 3, 0, 1 + rssky_transf->SMAX, 3 );
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
  SuperskyTransformData *rssky_transf,          ///< [in,out] Reduced supersky metric coordinate transform data
  const gsl_matrix *fitted_ssky_metric          ///< [in] Fitted supersky metric
  )
{

  // Check input
  XLAL_CHECK( decoupled_ssky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( fitted_ssky_metric != NULL, XLAL_EFAULT );

  // Copy fitted metric to decoupled metric
  gsl_matrix_memcpy( decoupled_ssky_metric, fitted_ssky_metric );

  // Create views of the sky--sky, freq+spin--freq+spin, and off-diagonal blocks
  gsl_matrix_view sky_sky     = gsl_matrix_submatrix( decoupled_ssky_metric, 0, 0, 3, 3 );
  gsl_matrix_view sky_fspin   = gsl_matrix_submatrix( decoupled_ssky_metric, 0, 3, 3, 1 + rssky_transf->SMAX );
  gsl_matrix_view fspin_sky   = gsl_matrix_submatrix( decoupled_ssky_metric, 3, 0, 1 + rssky_transf->SMAX, 3 );
  gsl_matrix_view fspin_fspin = gsl_matrix_submatrix( decoupled_ssky_metric, 3, 3, 1 + rssky_transf->SMAX, 1 + rssky_transf->SMAX );

  // Diagonal-normalise the freq+spin--freq+spin block
  gsl_matrix *fspin_fspin_dnorm = NULL, *fspin_fspin_dnorm_transf = NULL;
  XLAL_CHECK( XLALDiagNormalizeMetric( &fspin_fspin_dnorm, &fspin_fspin_dnorm_transf, &fspin_fspin.matrix ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Invert the freq+spin--freq+spin block
  gsl_matrix *GAMAT( fspin_fspin_dnorm_LU, 1 + rssky_transf->SMAX, 1 + rssky_transf->SMAX );
  gsl_matrix *GAMAT( fspin_fspin_dnorm_inv, 1 + rssky_transf->SMAX, 1 + rssky_transf->SMAX );
  gsl_permutation *GAPERM( fspin_fspin_dnorm_LU_perm, 1 + rssky_transf->SMAX );
  int fspin_fspin_dnorm_LU_sign = 0;
  gsl_matrix_memcpy( fspin_fspin_dnorm_LU, fspin_fspin_dnorm );
  XLAL_CHECK( gsl_linalg_LU_decomp( fspin_fspin_dnorm_LU, fspin_fspin_dnorm_LU_perm, &fspin_fspin_dnorm_LU_sign ) == 0, XLAL_EFAILED );
  XLAL_CHECK( gsl_linalg_LU_invert( fspin_fspin_dnorm_LU, fspin_fspin_dnorm_LU_perm, fspin_fspin_dnorm_inv ) == 0, XLAL_EFAILED );

  // Compute the additional sky offsets required to decouple the sky--sky and frequency blocks:
  //   decouple_sky_offsets = fspin_fspin_dnorm_transf * inv(fspin_fspin_dnorm) * fspin_fspin_dnorm_transf * fspin_sky
  // Uses fspin_sky as a temporary matrix, since it will be zeroed out anyway
  gsl_matrix *GAMAT( decouple_sky_offsets, 1 + rssky_transf->SMAX, 3 );
  gsl_blas_dtrmm( CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, fspin_fspin_dnorm_transf, &fspin_sky.matrix );
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, fspin_fspin_dnorm_inv, &fspin_sky.matrix, 0.0, decouple_sky_offsets );
  gsl_blas_dtrmm( CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, fspin_fspin_dnorm_transf, decouple_sky_offsets );

  // Add the additional sky offsets to the reduced supersky coordinate transform data
  gsl_matrix_view sky_offsets = gsl_matrix_view_array( &rssky_transf->sky_offsets[0][0], rssky_transf->nsky_offsets, 3 );
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
  SuperskyTransformData *rssky_transf,          ///< [in,out] Reduced supersky metric coordinate transform data
  const gsl_matrix *decoupled_ssky_metric       ///< [in] Decoupled supersky metric
  )
{

  // Check input
  XLAL_CHECK( aligned_ssky_metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( decoupled_ssky_metric != NULL, XLAL_EFAULT );

  // Allocate memory
  gsl_matrix *GAMAT( tmp, 1 + rssky_transf->SMAX, 3 );

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
  gsl_matrix_view align_sky = gsl_matrix_view_array( &rssky_transf->align_sky[0][0], 3, 3 );
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
  gsl_matrix_view sky_offsets = gsl_matrix_view_array( &rssky_transf->sky_offsets[0][0], rssky_transf->nsky_offsets, 3 );
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
  SuperskyTransformData **rssky_transf,         ///< [out] Reduced supersky metric coordinate transform data
  const size_t spindowns,                       ///< [in] Number of frequency+spindown coordinates
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
  *rssky_transf = XLALCalloc( 1, sizeof( **rssky_transf ) );
  XLAL_CHECK( *rssky_transf != NULL, XLAL_ENOMEM );
  gsl_matrix *GAMAT( assky_metric, 1 + orbital_metric->size1, 1 + orbital_metric->size1 );

  // Set coordinate transform data
  ( *rssky_transf )->ndim = ( *rssky_metric )->size1;
  ( *rssky_transf )->ref_time = *ref_time;
  ( *rssky_transf )->fiducial_freq = fiducial_calc_freq;
  ( *rssky_transf )->nspins = spindowns;
  ( *rssky_transf )->nsky_offsets = 1 + spindowns;

  // Compute the aligned supersky metric from the unrestricted supersky metric and the orbital metric
  XLAL_CHECK( SM_ComputeFittedSuperskyMetric( assky_metric, *rssky_transf, ussky_metric, orbital_metric, ocoords, start_time, end_time ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( SM_ComputeDecoupledSuperskyMetric( assky_metric, *rssky_transf, assky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( SM_ComputeAlignedSuperskyMetric( assky_metric, *rssky_transf, assky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );

  {
    // Move the row/column of the aligned supersky metric corresponding to the frequency to the last row/column
    const int ifreq = XLALFindDopplerCoordinateInSystem( ucoords, DOPPLERCOORD_FREQ );
    XLAL_CHECK( 3 <= ifreq, XLAL_EFAILED );   // 'ucoords' has 3 sky positions
    for ( size_t i = ifreq; i + 1 < assky_metric->size1; ++i ) {
      gsl_matrix_swap_rows( assky_metric, i, i + 1 );
      gsl_matrix_swap_columns( assky_metric, i, i + 1 );
    }

    // Move the row of the coordinate transform data corresponding to the frequency to the last row
    const int isky_offset_freq = ifreq - 3;   // 'ucoords' has 3 sky positions
    gsl_matrix_view sky_offsets = gsl_matrix_view_array(&( *rssky_transf )->sky_offsets[0][0], ( *rssky_transf )->nsky_offsets, 3);
    for ( size_t i = isky_offset_freq; i + 1 < sky_offsets.matrix.size1; ++i ) {
      gsl_matrix_swap_rows( &sky_offsets.matrix, i, i + 1 );
    }

    // Move the row/column of the aligned supersky metric corresponding to the 'n_c' sky coordinate to the last row/column, so it can be easily dropped
    const int inc = XLALFindDopplerCoordinateInSystem( ucoords, DOPPLERCOORD_N3Z_EQU );
    XLAL_CHECK( 0 <= inc && inc < ifreq, XLAL_EFAILED );
    for ( size_t i = inc; i + 1 < assky_metric->size1; ++i ) {
      gsl_matrix_swap_rows( assky_metric, i, i + 1 );
      gsl_matrix_swap_columns( assky_metric, i, i + 1 );
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
  const SuperskyMetricType type,
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
  XLAL_CHECK_NULL( type < MAX_METRIC_TYPE, XLAL_EINVAL );
  XLAL_CHECK_NULL( spindowns <= 4, XLAL_EINVAL );
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
  if ( spindowns >= 4 ) {
    ucoords.coordIDs[ucoords.dim++] = DOPPLERCOORD_F4DOT;
    ocoords.coordIDs[ocoords.dim++] = DOPPLERCOORD_F4DOT;
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
    XLAL_CHECK_NULL( SM_ComputeReducedSuperskyMetric( &metrics->coh_rssky_metric[n], &metrics->coh_rssky_transf[n], spindowns, ussky_metric_seg, &ucoords, orbital_metric_seg, &ocoords, ref_time, start_time_seg, end_time_seg ) == XLAL_SUCCESS, XLAL_EFUNC );
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
  XLAL_CHECK_NULL( SM_ComputeReducedSuperskyMetric( &metrics->semi_rssky_metric, &metrics->semi_rssky_transf, spindowns, ussky_metric_avg, &ucoords, orbital_metric_avg, &ocoords, ref_time, start_time_avg, end_time_avg ) == XLAL_SUCCESS, XLAL_EFUNC );
  LogPrintf( LOG_DEBUG, "Computed semicoherent reduced supersky metric for %zu segments\n", metrics->num_segments );

  // Rescale metrics to input fiducial frequency
  XLALScaleSuperskyMetricsFiducialFreq( metrics, fiducial_freq );

  // Cleanup
  GFMAT( ussky_metric_avg, orbital_metric_avg );

  return metrics;

}

SuperskyMetrics *XLALCopySuperskyMetrics(
  const SuperskyMetrics *metrics
  )
{

  // Check input
  XLAL_CHECK_NULL( metrics != NULL, XLAL_EFAULT );

  // Allocate memory for copied struct
  SuperskyMetrics *copy_metrics = XLALCalloc( 1, sizeof( *copy_metrics ) );
  XLAL_CHECK_NULL( copy_metrics != NULL, XLAL_ENOMEM );
  copy_metrics->num_segments = metrics->num_segments;

  // Allocate memory for copies of arrays of coherent metrics
  copy_metrics->coh_rssky_metric = XLALCalloc( copy_metrics->num_segments, sizeof( *copy_metrics->coh_rssky_metric ) );
  XLAL_CHECK_NULL( copy_metrics->coh_rssky_metric != NULL, XLAL_ENOMEM );
  copy_metrics->coh_rssky_transf = XLALCalloc( copy_metrics->num_segments, sizeof( *copy_metrics->coh_rssky_transf ) );
  XLAL_CHECK_NULL( copy_metrics->coh_rssky_transf != NULL, XLAL_ENOMEM );

  // Copy coherent metrics and transform data
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    GAMAT_NULL( copy_metrics->coh_rssky_metric[n], metrics->coh_rssky_metric[n]->size1, metrics->coh_rssky_metric[n]->size2 );
    gsl_matrix_memcpy( copy_metrics->coh_rssky_metric[n], metrics->coh_rssky_metric[n] );
    copy_metrics->coh_rssky_transf[n] = XLALCalloc( 1, sizeof( *copy_metrics->coh_rssky_transf[n] ) );
    XLAL_CHECK_NULL( copy_metrics->coh_rssky_transf[n] != NULL, XLAL_ENOMEM );
    memcpy( copy_metrics->coh_rssky_transf[n], metrics->coh_rssky_transf[n], sizeof( *copy_metrics->coh_rssky_transf[n] ) );
  }

  // Copy semocherent metrics and transform data
  GAMAT_NULL( copy_metrics->semi_rssky_metric, metrics->semi_rssky_metric->size1, metrics->semi_rssky_metric->size2 );
  gsl_matrix_memcpy( copy_metrics->semi_rssky_metric, metrics->semi_rssky_metric );
  copy_metrics->semi_rssky_transf = XLALCalloc( 1, sizeof( *copy_metrics->semi_rssky_transf ) );
  XLAL_CHECK_NULL( copy_metrics->semi_rssky_transf != NULL, XLAL_ENOMEM );
  memcpy( copy_metrics->semi_rssky_transf, metrics->semi_rssky_transf, sizeof( *copy_metrics->semi_rssky_transf ) );

  return copy_metrics;

}

void XLALDestroySuperskyMetrics(
  SuperskyMetrics *metrics
  )
{
  if ( metrics != NULL ) {
    for ( size_t n = 0; n < metrics->num_segments; ++n ) {
      GFMAT( metrics->coh_rssky_metric[n] );
      XLALFree( metrics->coh_rssky_transf[n] );
    }
    XLALFree( metrics->coh_rssky_metric );
    XLALFree( metrics->coh_rssky_transf );
    GFMAT( metrics->semi_rssky_metric );
    XLALFree( metrics->semi_rssky_transf );
    XLALFree( metrics );
  }
}

///
/// Initialise a FITS table for writing/reading a table of SuperskyTransformData entries
///
static int fits_table_init_SuperskyTransformData(
  FITSFile *file
  )
{
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_FITS_TABLE_COLUMN_BEGIN( SuperskyTransformData );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, UINT4, ndim ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, ref_time ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, fiducial_freq ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, UINT4, nspins ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY2( file, REAL8, align_sky ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, UINT4, nsky_offsets ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_ARRAY2( file, REAL8, sky_offsets ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;
}

int XLALFITSWriteSuperskyMetrics(
  FITSFile *file,
  const SuperskyMetrics *metrics
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( metrics != NULL, XLAL_EFAULT );

  // Write coherent metrics to a FITS array
  {
    size_t dims[3] = {
      metrics->coh_rssky_metric[0]->size1,
      metrics->coh_rssky_metric[0]->size2,
      metrics->num_segments
    };
    XLAL_CHECK( XLALFITSArrayOpenWrite( file, "coh_rssky_metric", 3, dims, "coherent supersky metrics" ) == XLAL_SUCCESS, XLAL_EFUNC );
    size_t idx[3];
    for ( idx[2] = 0; idx[2] < dims[2]; ++idx[2] ) {
      XLAL_CHECK( XLALFITSArrayWriteGSLMatrix( file, idx, metrics->coh_rssky_metric[idx[2]] ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Write coherent metric transform data to a FITS table
  {
    XLAL_CHECK( XLALFITSTableOpenWrite( file, "coh_rssky_transf", "coherent supersky metric transform data" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( fits_table_init_SuperskyTransformData( file ) == XLAL_SUCCESS, XLAL_EFUNC );
    for ( size_t i = 0; i < metrics->num_segments; ++i ) {
      XLAL_CHECK( XLALFITSTableWriteRow( file, metrics->coh_rssky_transf[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Write semicoherent metric to a FITS array
  {
    XLAL_CHECK( XLALFITSArrayOpenWrite2( file, "semi_rssky_metric", metrics->semi_rssky_metric->size1, metrics->semi_rssky_metric->size2, "semicoherent supersky metric" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSArrayWriteGSLMatrix( file, NULL, metrics->semi_rssky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write semicoherent metric transform data to a FITS table
  {
    XLAL_CHECK( XLALFITSTableOpenWrite( file, "semi_rssky_transf", "semicoherent supersky metric transform data" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( fits_table_init_SuperskyTransformData( file ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSTableWriteRow( file, metrics->semi_rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}

int XLALFITSReadSuperskyMetrics(
  FITSFile *file,
  SuperskyMetrics **metrics
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( metrics != NULL && *metrics == NULL, XLAL_EFAULT );

  // Allocate memory
  *metrics = XLALCalloc( 1, sizeof( **metrics ) );
  XLAL_CHECK( *metrics != NULL, XLAL_ENOMEM );

  // Read coherent metrics from a FITS array
  {
    size_t ndim, dims[FFIO_MAX];
    XLAL_CHECK( XLALFITSArrayOpenRead( file, "coh_rssky_metric", &ndim, dims ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( ndim == 3, XLAL_EIO );
    ( *metrics )->num_segments = dims[2];
    ( *metrics )->coh_rssky_metric = XLALCalloc( ( *metrics )->num_segments, sizeof( *( *metrics )->coh_rssky_metric ) );
    XLAL_CHECK( ( *metrics )->coh_rssky_metric != NULL, XLAL_ENOMEM );
    size_t idx[3];
    for ( idx[2] = 0; idx[2] < dims[2]; ++idx[2] ) {
      XLAL_CHECK( XLALFITSArrayReadGSLMatrix( file, idx, &( *metrics )->coh_rssky_metric[idx[2]] ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Read coherent metric transform data from a FITS table
  {
    UINT8 nrows = 0;
    XLAL_CHECK( XLALFITSTableOpenRead( file, "coh_rssky_transf", &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( nrows == ( *metrics )->num_segments, XLAL_EIO );
    XLAL_CHECK( fits_table_init_SuperskyTransformData( file ) == XLAL_SUCCESS, XLAL_EFUNC );
    ( *metrics )->coh_rssky_transf = XLALCalloc( ( *metrics )->num_segments, sizeof( *( *metrics )->coh_rssky_transf ) );
    XLAL_CHECK( ( *metrics )->coh_rssky_transf != NULL, XLAL_ENOMEM );
    for ( size_t i = 0; i < ( *metrics )->num_segments; ++i ) {
      ( *metrics )->coh_rssky_transf[i] = XLALCalloc( 1, sizeof( *( *metrics )->coh_rssky_transf[i] ) );
      XLAL_CHECK( ( *metrics )->coh_rssky_transf[i] != NULL, XLAL_ENOMEM );
      XLAL_CHECK( XLALFITSTableReadRow( file, ( *metrics )->coh_rssky_transf[i], &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Read semicoherent metric from a FITS array
  {
    XLAL_CHECK( XLALFITSArrayOpenRead2( file, "semi_rssky_metric", NULL, NULL ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSArrayReadGSLMatrix( file, NULL, &( *metrics )->semi_rssky_metric ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Read semicoherent metric transform data from a FITS table
  {
    UINT8 nrows = 0;
    XLAL_CHECK( XLALFITSTableOpenRead( file, "semi_rssky_transf", &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( nrows == 1, XLAL_EIO );
    XLAL_CHECK( fits_table_init_SuperskyTransformData( file ) == XLAL_SUCCESS, XLAL_EFUNC );
    ( *metrics )->semi_rssky_transf = XLALCalloc( 1, sizeof( *( *metrics )->semi_rssky_transf ) );
    XLAL_CHECK( ( *metrics )->semi_rssky_transf != NULL, XLAL_ENOMEM );
    XLAL_CHECK( XLALFITSTableReadRow( file, ( *metrics )->semi_rssky_transf, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

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

  // Return dimensions
  if ( spindowns != NULL ) {
    *spindowns = metrics->semi_rssky_transf->nspins;
  }

  return XLAL_SUCCESS;

}

///
/// Scale a given supersky metric and its coordinate transform data to a new fiducial frequency.
///
static int SM_ScaleSuperskyMetricFiducialFreq(
  gsl_matrix *rssky_metric,                     ///< [in] Reduced supersky metric
  SuperskyTransformData *rssky_transf,          ///< [in] Reduced supersky metric coordinate transform data
  const double new_fiducial_freq                ///< [in] New fiducial frequency
  )
{

  // Check input
  XLAL_CHECK( CHECK_RSSKY_METRIC_TRANSF( rssky_metric, rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( new_fiducial_freq > 0, XLAL_EINVAL );

  // Rescale metrics to 'new_fiducial_freq' based on known scalings
  const double fiducial_scale = new_fiducial_freq / rssky_transf->fiducial_freq;
  gsl_matrix_view sky_sky = gsl_matrix_submatrix( rssky_metric, 0, 0, 2, 2 );
  gsl_matrix_view sky_offsets = gsl_matrix_view_array( &rssky_transf->sky_offsets[0][0], rssky_transf->nsky_offsets, 3 );
  gsl_matrix_scale( &sky_sky.matrix, SQR( fiducial_scale ) );
  gsl_matrix_scale( &sky_offsets.matrix, fiducial_scale );

  // Set new fiducial frequency
  rssky_transf->fiducial_freq = new_fiducial_freq;

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
    XLAL_CHECK( SM_ScaleSuperskyMetricFiducialFreq( metrics->coh_rssky_metric[n], metrics->coh_rssky_transf[n], new_fiducial_freq ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  XLAL_CHECK( SM_ScaleSuperskyMetricFiducialFreq( metrics->semi_rssky_metric, metrics->semi_rssky_transf, new_fiducial_freq ) == XLAL_SUCCESS, XLAL_EFUNC );

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

  // Find the maximum, over both semicoherent and coherent metrics, of 'g_{ff} / mu',
  // where 'g_{ff}' is the frequency-frequency metric element, and mu is the mismatch
  const size_t ifreq = RSSKY_FKDOT_DIM( metrics->semi_rssky_transf, 0 );
  double max_gff_d_mu = gsl_matrix_get( metrics->semi_rssky_metric, ifreq, ifreq ) / semi_max_mismatch;
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    const double coh_gff_d_mu = gsl_matrix_get( metrics->coh_rssky_metric[n], ifreq, ifreq ) / coh_max_mismatch;
    if ( max_gff_d_mu < coh_gff_d_mu ) {
      max_gff_d_mu = coh_gff_d_mu;
    }
  }

  // Rescale rows 'g_{f*}' and columns 'g_{*f}' of both semicoherent and coherent metrics by
  // 'sqrt( max{ g_{ff} / mu } / ( g_{ff} / mu ) )'; this scales the frequency coordinate such
  // that 'g_{ff}' will give the same frequency spacing for all metrics, while also scaling the
  // off-diagonal frequency terms (e.g. frequency-spindown correlations) correctly.
  {
    const double semi_gff_d_mu = gsl_matrix_get( metrics->semi_rssky_metric, ifreq, ifreq ) / semi_max_mismatch;
    const double scale = sqrt( max_gff_d_mu / semi_gff_d_mu );
    XLAL_CHECK( scale >= 1, XLAL_EFAILED );
    {
      gsl_vector_view rssky_metric_f_row = gsl_matrix_row( metrics->semi_rssky_metric, ifreq );
      gsl_vector_scale( &rssky_metric_f_row.vector, scale );
    } {
      gsl_vector_view rssky_metric_f_col = gsl_matrix_column( metrics->semi_rssky_metric, ifreq );
      gsl_vector_scale( &rssky_metric_f_col.vector, scale );
    }
  }
  for ( size_t n = 0; n < metrics->num_segments; ++n ) {
    const double coh_gff_d_mu = gsl_matrix_get( metrics->coh_rssky_metric[n], ifreq, ifreq ) / coh_max_mismatch;
    const double scale = sqrt( max_gff_d_mu / coh_gff_d_mu );
    XLAL_CHECK( scale >= 1, XLAL_EFAILED );
    {
      gsl_vector_view rssky_metric_f_row = gsl_matrix_row( metrics->coh_rssky_metric[n], ifreq );
      gsl_vector_scale( &rssky_metric_f_row.vector, scale );
    } {
      gsl_vector_view rssky_metric_f_col = gsl_matrix_column( metrics->coh_rssky_metric[n], ifreq );
      gsl_vector_scale( &rssky_metric_f_col.vector, scale );
    }
  }

  return XLAL_SUCCESS;

}

int XLALSetPhysicalPointSuperskyRefTime(
  PulsarDopplerParams *out_phys,
  const SuperskyTransformData *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EFAULT );

  // Set output physical point reference time to that of of coordinate transform data
  out_phys->refTime = rssky_transf->ref_time;

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
  const gsl_vector *rss,                        ///< [in] 2-dimensional reduced supersky coordinates
  const double hemi                             ///< [in] Supersky coordinate hemisphere; only sign is used
  )
{
  const double A = GSL_SIGN( hemi ) * gsl_vector_get( rss, 0 ) - 1;
  const double B = gsl_vector_get( rss, 1 );
  const double R = sqrt( SQR( A ) + SQR( B ) );
  const double Rmax = GSL_MAX( 1.0, R );
  as[0] = A / Rmax;
  as[1] = B / Rmax;
  as[2] = GSL_SIGN( hemi ) * RE_SQRT( 1.0 - DOT2( as, as ) );
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
  const SuperskyTransformData *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EFAULT );
  XLAL_CHECK( out_rssky->size == rssky_transf->ndim, XLAL_ESIZE );

  // Transform input physical point to reference time of coordinate transform data
  PulsarDopplerParams in_phys_ref = *in_phys;
  {
    const REAL8 dtau = XLALGPSDiff( &rssky_transf->ref_time, &in_phys_ref.refTime );
    XLAL_CHECK( XLALExtrapolatePulsarSpins( in_phys_ref.fkdot, in_phys_ref.fkdot, dtau ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Create array for intermediate coordinates
  double intm[4 + rssky_transf->SMAX];
  gsl_vector_view intm_sky2 = gsl_vector_view_array( &intm[0], 2 );
  gsl_vector_view intm_sky3 = gsl_vector_view_array( &intm[0], 3 );
  gsl_vector_view intm_fspin = gsl_vector_view_array( &intm[3], 1 + rssky_transf->SMAX );

  // Convert right ascension and declination to equatorial coordinates
  {
    const double cos_Delta = cos( in_phys_ref.Delta );
    intm[0] = cos( in_phys_ref.Alpha ) * cos_Delta;
    intm[1] = sin( in_phys_ref.Alpha ) * cos_Delta;
    intm[2] = sin( in_phys_ref.Delta );
  }

  // Copy frequency/spindowns to intermediate array; frequency goes last
  intm[3 + rssky_transf->SMAX] = in_phys_ref.fkdot[0];
  for ( size_t s = 1; s <= rssky_transf->SMAX; ++s ) {
    intm[2 + s] = in_phys_ref.fkdot[s];
  }

  // Apply the alignment transform to the supersky position to produced the aligned sky position:
  //   asky = align_sky * ssky
  double asky[3];
  gsl_vector_view asky_v = gsl_vector_view_array( asky, 3 );
  gsl_matrix_const_view align_sky = gsl_matrix_const_view_array( &rssky_transf->align_sky[0][0], 3, 3 );
  gsl_blas_dgemv( CblasNoTrans, 1.0, &align_sky.matrix, &intm_sky3.vector, 0.0, &asky_v.vector );

  // Add the inner product of the sky offsets with the aligned sky position
  // to the supersky spins and frequency to get the reduced supersky quantities:
  //   rssky_fspin[i] = ussky_fspin[i] + dot(sky_offsets[i], asky)
  gsl_matrix_const_view sky_offsets = gsl_matrix_const_view_array( &rssky_transf->sky_offsets[0][0], rssky_transf->nsky_offsets, 3 );
  gsl_blas_dgemv( CblasNoTrans, 1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &intm_fspin.vector );

  // Convert from 3-dimensional aligned sky coordinates to 2-dimensional reduced supersky coordinates
  SM_AlignedToReduced( &intm_sky3.vector, asky );

  // Copy intermediate array to output point
  {
    gsl_vector_view out_sky2 = gsl_vector_subvector( out_rssky, 0, 2 );
    gsl_vector_view out_fspin = gsl_vector_subvector( out_rssky, 2, 1 + rssky_transf->SMAX );
    gsl_vector_memcpy( &out_sky2.vector, &intm_sky2.vector );
    gsl_vector_memcpy( &out_fspin.vector, &intm_fspin.vector );
  }

  return XLAL_SUCCESS;

}

int XLALConvertSuperskyToPhysicalPoint(
  PulsarDopplerParams *out_phys,
  const gsl_vector *in_rssky,
  const gsl_vector *ref_rssky,
  const SuperskyTransformData *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EFAULT );
  XLAL_CHECK( in_rssky->size == rssky_transf->ndim, XLAL_ESIZE );
  XLAL_CHECK( ref_rssky == NULL || in_rssky->size == rssky_transf->ndim, XLAL_ESIZE );

  // Set output physical point reference time to that of of coordinate transform data
  out_phys->refTime = rssky_transf->ref_time;

  // Create array for intermediate coordinates
  double intm[4 + rssky_transf->SMAX];
  gsl_vector_view intm_sky2 = gsl_vector_view_array( &intm[0], 2 );
  gsl_vector_view intm_sky3 = gsl_vector_view_array( &intm[0], 3 );
  gsl_vector_view intm_fspin = gsl_vector_view_array( &intm[3], 1 + rssky_transf->SMAX );

  // Copy input point to intermediate array
  {
    gsl_vector_const_view in_sky2 = gsl_vector_const_subvector( in_rssky, 0, 2 );
    gsl_vector_const_view in_fspin = gsl_vector_const_subvector( in_rssky, 2, 1 + rssky_transf->SMAX );
    gsl_vector_memcpy( &intm_sky2.vector, &in_sky2.vector );
    gsl_vector_memcpy( &intm_fspin.vector, &in_fspin.vector );
  }

  // Convert from 2-dimensional reduced supersky coordinates to 3-dimensional aligned sky coordinates
  const double A = gsl_vector_get( ref_rssky != NULL ? ref_rssky : in_rssky, 0 );
  double asky[3];
  SM_ReducedToAligned( asky, &intm_sky3.vector, A );
  gsl_vector_view asky_v = gsl_vector_view_array( asky, 3 );

  // Subtract the inner product of the sky offsets with the aligned sky position
  // from the reduced supersky spins and frequency to get the supersky quantities:
  //   ussky_fspin[i] = rssky_fspin[i] - dot(sky_offsets[i], asky)
  gsl_matrix_const_view sky_offsets = gsl_matrix_const_view_array( &rssky_transf->sky_offsets[0][0], rssky_transf->nsky_offsets, 3 );
  gsl_blas_dgemv( CblasNoTrans, -1.0, &sky_offsets.matrix, &asky_v.vector, 1.0, &intm_fspin.vector );

  // Apply the inverse alignment transform to the aligned sky position to produced the supersky position:
  //   ssky = align_sky^T * asky
  gsl_matrix_const_view align_sky = gsl_matrix_const_view_array( &rssky_transf->align_sky[0][0], 3, 3 );
  gsl_blas_dgemv( CblasTrans, 1.0, &align_sky.matrix, &asky_v.vector, 0.0, &intm_sky3.vector );

  // Copy frequency/spindowns to output physical point; frequency goes first
  out_phys->fkdot[0] = intm[3 + rssky_transf->SMAX];
  for ( size_t s = 1; s <= rssky_transf->SMAX; ++s ) {
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
  const SuperskyTransformData *out_rssky_transf,
  const gsl_vector *in_rssky,
  const gsl_vector *ref_rssky,
  const SuperskyTransformData *in_rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( out_rssky_transf ), XLAL_EFAULT );
  XLAL_CHECK( in_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( in_rssky_transf ), XLAL_EFAULT );

  // Convert input reduced supersky point to physical coordinates
  PulsarDopplerParams XLAL_INIT_DECL( phys );
  XLAL_CHECK( XLALConvertSuperskyToPhysicalPoint( &phys, in_rssky, ref_rssky, in_rssky_transf ) == XLAL_SUCCESS, XLAL_EINVAL );

  // Convert physical point to output reduced supersky coordinates
  XLAL_CHECK( XLALConvertPhysicalToSuperskyPoint( out_rssky, &phys, out_rssky_transf ) == XLAL_SUCCESS, XLAL_EINVAL );

  return XLAL_SUCCESS;

}

int XLALConvertPhysicalToSuperskyPoints(
  gsl_matrix **out_rssky,
  const gsl_matrix *in_phys,
  const SuperskyTransformData *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EFAULT );
  XLAL_CHECK( in_phys->size1 == rssky_transf->ndim, XLAL_ESIZE );

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
    in_phys_j.refTime = rssky_transf->ref_time;
    in_phys_j.Alpha = gsl_matrix_get( in_phys, 0, j );
    in_phys_j.Delta = gsl_matrix_get( in_phys, 1, j );
    for ( size_t s = 0; s <= rssky_transf->SMAX; ++s ) {
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
  const SuperskyTransformData *rssky_transf
  )
{

  // Check input
  XLAL_CHECK( out_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( in_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EFAULT );
  XLAL_CHECK( in_rssky->size1 == rssky_transf->ndim, XLAL_ESIZE );

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
    XLAL_CHECK( XLALConvertSuperskyToPhysicalPoint( &out_phys_j, &in_rssky_j.vector, NULL, rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Fill output point from PulsarDopplerParams struct
    gsl_matrix_set( *out_phys, 0, j, out_phys_j.Alpha );
    gsl_matrix_set( *out_phys, 1, j, out_phys_j.Delta );
    for ( size_t s = 0; s <= rssky_transf->SMAX; ++s ) {
      gsl_matrix_set( *out_phys, 2 + s, j, out_phys_j.fkdot[s] );
    }

  }

  return XLAL_SUCCESS;

}

static void SkyBoundCache(
  const size_t dim UNUSED,
  const gsl_vector *point,
  gsl_vector* cache
  )
{

  // Convert from 2-dimensional reduced supersky coordinates to 3-dimensional aligned sky coordinates
  const double A = gsl_vector_get( point, 0 );
  double as[3];
  SM_ReducedToAligned( as, point, A );

  // Store aligned sky position in cache
  for ( size_t i = 0; i < 3; ++i ) {
    gsl_vector_set( cache, i, as[i] );
  }

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
  const gsl_matrix *cache UNUSED,
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
  SM_ReducedToAligned( as, &rssky_view.vector, A );
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
  SuperskyTransformData *rssky_transf,
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
  gsl_matrix_view align_sky = gsl_matrix_view_array( &rssky_transf->align_sky[0][0], 3, 3 );
  gsl_matrix_view sky_offsets = gsl_matrix_view_array( &rssky_transf->sky_offsets[0][0], rssky_transf->nsky_offsets, 3 );

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

  // Set parameter-space bound names
  XLAL_CHECK( XLALSetLatticeTilingBoundName( tiling, 0, "sskyA" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingBoundName( tiling, 1, "sskyB" ) == XLAL_SUCCESS, XLAL_EFUNC );

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
    XLAL_CHECK( XLALSetLatticeTilingBoundCacheFunction( tiling, 1, SkyBoundCache ) == XLAL_SUCCESS, XLAL_EFUNC );

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
    XLAL_CHECK( XLALSetLatticeTilingBoundCacheFunction( tiling, 1, SkyBoundCache ) == XLAL_SUCCESS, XLAL_EFUNC );

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
  XLAL_CHECK( XLALSetLatticeTilingBoundCacheFunction( tiling, 1, SkyBoundCache ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Set the parameter-space origin on reduced supersky sky coordinates A and B
  XLAL_CHECK( XLALSetLatticeTilingOrigin( tiling, 0, 0.0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingOrigin( tiling, 1, 0.0 ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

static double ConstantBoundB(
  const void *data,
  const size_t dim UNUSED,
  const gsl_matrix *cache UNUSED,
  const gsl_vector *point
  )
{

  // Get bounds data
  const double bound = *( ( const double * ) data );

  // Decode the reduced supersky coordinates to get
  //   na = as[0] = Q_na . n
  const double A = gsl_vector_get( point, 0 );
  const double rssky[2] = { A, 0 };
  gsl_vector_const_view rssky_view = gsl_vector_const_view_array( rssky, 2 );
  double as[3];
  SM_ReducedToAligned( as, &rssky_view.vector, A );
  const double na = as[0];

  // Absolute limiting bound on 'nb = +/- sqrt(1 - na^2)'
  const double limit = RE_SQRT( 1 - SQR( na ) );

  return GSL_MAX( -limit, GSL_MIN( bound, limit ) );

}

static double EqualAreaSkyBoundSolverA(
  double A1,
  void *params
  )
{

  // Get parameters
  const double target_area = ( ( double* ) params )[0];

  // Compute area of unit disk to the left of 'A1'
  const double area = LAL_PI_2 + A1*RE_SQRT( 1 - SQR( A1 ) ) + asin( A1 );

  return area - target_area;

}

static double EqualAreaSkyBoundSolverB(
  double B1,
  void *params
  )
{

  // Get parameters
  const double target_area = ( ( double* ) params )[0];
  double A0 = ( ( double* ) params )[1];
  double A1 = ( ( double* ) params )[2];

  // Compute area of unit disk between 'A0' and 'A1'
  const double max_area = A1*RE_SQRT( 1 - SQR( A1 ) ) - A0*RE_SQRT( 1 - SQR( A0 ) ) + asin( A1 ) - asin( A0 );

  // Work out where '-|B1| = const' line intersects unit disk boundary
  const double Ai = RE_SQRT( 1 - SQR( B1 ) );

  // Restrict range '[A0, A1]' if '-|B1| = const' line intersects unit disk boundary within it
  if ( A0 < -Ai ) {
    A0 = -Ai;
  }
  if ( A1 > Ai ) {
    A1 = Ai;
  }

  // Compute area of unit disk between 'A0' and 'A1' and below '-|B1|'
  double area = -fabs( B1 )*( A1 - A0 ) + 0.5*( A1*RE_SQRT( 1 - SQR( A1 ) ) - A0*RE_SQRT( 1 - SQR( A0 ) ) + asin( A1 ) - asin( A0 ) );

  // For positive B, substract 'area' from 'max_area'
  if ( B1 > 0 ) {
    area = max_area - area;
  }

  return area - target_area;

}

int XLALSetSuperskyEqualAreaSkyBounds(
  LatticeTiling *tiling,
  const gsl_matrix *rssky_metric,
  const double max_mismatch,
  const UINT4 patch_count,
  const UINT4 patch_index
  )
{

  // Check input
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( rssky_metric->size1 == rssky_metric->size2, XLAL_EINVAL );
  XLAL_CHECK( rssky_metric->size1 >= 2, XLAL_EINVAL );
  XLAL_CHECK( gsl_matrix_get( rssky_metric, 0, 1 ) == 0, XLAL_EINVAL );
  XLAL_CHECK( gsl_matrix_get( rssky_metric, 1, 0 ) == 0, XLAL_EINVAL );
  XLAL_CHECK( max_mismatch > 0, XLAL_EINVAL );
  XLAL_CHECK( patch_count > 0, XLAL_EINVAL );
  XLAL_CHECK( patch_count == 1 || patch_count % 2 == 0, XLAL_EINVAL, "'patch_count' must be either one or even" );
  XLAL_CHECK( patch_index < patch_count, XLAL_EINVAL );

  // Set parameter-space bound names
  XLAL_CHECK( XLALSetLatticeTilingBoundName( tiling, 0, "sskyA" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingBoundName( tiling, 1, "sskyB" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Parameter-space bounds on reduced supersky sky coordinates A and B
  double A_bound[2] = {GSL_NAN, GSL_NAN};
  double B_bound[2] = {-1, 1};
  LatticeTilingPaddingFlags A_padf = LATTICE_TILING_PAD_LHBBX | LATTICE_TILING_PAD_UHBBX;
  LatticeTilingPaddingFlags B_padf = LATTICE_TILING_PAD_LHBBX | LATTICE_TILING_PAD_UHBBX;

  // Handle special cases of 1 and 2 patches
  if ( patch_count <= 2 ) {
    A_bound[0] = -2 + 4 * ( ( double ) patch_index ) / ( ( double ) patch_count );
    A_bound[1] = A_bound[0] + 4 / ( ( double ) patch_count );
  } else {

    // Number of patches per hemisphere for odd number of patches
    const UINT4 hemi_patch_count = patch_count / 2;

    // Index of patch within hemisphere; increases to 'hemi_patch_count' then decreases to zero
    const UINT4 hemi_patch_index = patch_index < hemi_patch_count ? patch_index : patch_count - patch_index - 1;

    // Minimum number of patch divisions in supersky coordinate A within one hemisphere
    // - Prevents patches from being too squashed in either A or B, which lead to overcoverage (i.e. the
    //   sum of points over all sky patches being much larger than the number of points in the whole sky)
    const UINT4 min_A_count = ceil( sqrt( hemi_patch_count ) );

    // Maximum extent of metric ellipse bounding box in supersky coordinate B
    const double max_B_bbox = 2 * sqrt( max_mismatch / gsl_matrix_get( rssky_metric, 1, 1 ) );

    // Target number of patch divisions in supersky coordinate B
    const double target_B_count = 1.0 / ( 1.0 * max_B_bbox );

    // Number of patch divisions in supersky coordinate A within one hemisphere
    const UINT4 A_count = GSL_MIN( min_A_count, ceil( hemi_patch_count / target_B_count ) );

    // Minimum number of patch divisions in supersky coordinate B
    const UINT4 min_B_count = hemi_patch_count / A_count;

    // Excess number of patches, which must be added on to get 'hemi_patch_count'
    INT4 patch_excess = hemi_patch_count - A_count * min_B_count;
    XLAL_CHECK( patch_excess >= 0, XLAL_EFAILED );

    // Initialise number of patch divisions in 'B'; if there are excess patches, add an extra patch
    UINT4 B_count = min_B_count;
    if ( patch_excess > 0 ) {
      ++B_count;
    }

    // Calculate range of indices in 'A', and number of patch divisions and index in 'B'.
    // - The divisions in 'A' are set in proportion to the range of 'A_index', i.e. the number of
    //   divisions in 'B' for that range of 'A_index'. This is so that, if 'patch_excess' is
    //   not zero, and therefore the number of divisions in 'B' is not constant, patch areas should
    //   still be equal. Example:
    //     hemi_patch_count=7 hemi_patch_index=0 | A_index=0--3 B_count=3 B_index=0
    //     hemi_patch_count=7 hemi_patch_index=1 | A_index=0--3 B_count=3 B_index=1
    //     hemi_patch_count=7 hemi_patch_index=2 | A_index=0--3 B_count=3 B_index=2
    //     hemi_patch_count=7 hemi_patch_index=3 | A_index=3--5 B_count=2 B_index=0
    //     hemi_patch_count=7 hemi_patch_index=4 | A_index=3--5 B_count=2 B_index=1
    //     hemi_patch_count=7 hemi_patch_index=5 | A_index=5--7 B_count=2 B_index=0
    //     hemi_patch_count=7 hemi_patch_index=6 | A_index=5--7 B_count=2 B_index=1
    UINT4 A_index[2] = {0, B_count};
    UINT4 B_index = hemi_patch_index;
    while ( B_index >= B_count ) {

      // Decrease index in 'B'; we are done when 'B_index' < 'B_count'
      B_index -= B_count;

      // Decrease number of excess patches; if zero, subtract extra patch from patch divisions in 'B'
      --patch_excess;
      if ( patch_excess == 0 ) {
        --B_count;
      }

      // Store the current last 'A' index in 'A_index[0]', and increase
      // 'A_index[1]' by the current number of patch divisions in 'B'
      A_index[0] = A_index[1];
      A_index[1] += B_count;

    }

    // Decide which patches to add padding to
    A_padf = ( ( A_index[0] == 0 && patch_index <  hemi_patch_count ) ? LATTICE_TILING_PAD_LHBBX : 0 )
      |      ( ( A_index[0] == 0 && patch_index >= hemi_patch_count ) ? LATTICE_TILING_PAD_UHBBX : 0 );
    B_padf = ( ( B_index     == 0       ) ? LATTICE_TILING_PAD_LHBBX : 0 )
      |      ( ( B_index + 1 == B_count ) ? LATTICE_TILING_PAD_UHBBX : 0 );

    // Allocate a GSL root solver
    gsl_root_fsolver *fs = gsl_root_fsolver_alloc( gsl_root_fsolver_brent );
    XLAL_CHECK( fs != NULL, XLAL_ENOMEM );

    // Find bounds on 'A' corresponding to the computed indexes
    for ( size_t i = 0; i < 2; ++i ) {
      const UINT4 A_index_i = A_index[i];

      // Handle boundaries as special cases
      if ( A_index_i == 0 ) {
        A_bound[i] = -1;
      } else if ( A_index_i == hemi_patch_count ) {
        A_bound[i] = 1;
      } else {

        // Calculate the target area of unit disk
        const double target_area = LAL_PI * ( ( double ) A_index_i ) / ( ( ( double ) patch_count ) / 2 );

        // Set up GSL root solver
        double params[] = { target_area };
        gsl_function F = { .function = EqualAreaSkyBoundSolverA, .params = params };
        double A_lower = -1, A_upper = 1;
        XLAL_CHECK( gsl_root_fsolver_set( fs, &F, A_lower, A_upper ) == 0, XLAL_EFAILED );

        // Try to find root
        int status = 0, iter = 0;
        do {
          XLAL_CHECK( gsl_root_fsolver_iterate( fs ) == 0, XLAL_EFAILED );
          A_lower = gsl_root_fsolver_x_lower( fs );
          A_upper = gsl_root_fsolver_x_upper( fs );
          status = gsl_root_test_interval( A_lower, A_upper, 1e-5, 1e-5 );
        } while (status == GSL_CONTINUE && ++iter < 1000);
        XLAL_CHECK( status == GSL_SUCCESS, XLAL_EMAXITER, "Could not find bound for A_index[%zu]=%i; best guess [%g, %g]", i, A_index_i, A_lower, A_upper );

        // Store bound
        A_bound[i] = gsl_root_fsolver_root( fs );

      }

    }

    // Find bounds on 'B' corresponding to the computed indexes
    for ( size_t i = 0; i < 2; ++i ) {
      const UINT4 B_index_i = B_index + i;

      // Maximum possible value of 'B' within the region bound by 'A_bound'
      // - For a single patch in 'A', 'B' must span the entire unit disk
      const double B_max = ( A_count == 1 ) ? 1 : RE_SQRT( 1 - GSL_MIN( SQR( A_bound[0] ), SQR( A_bound[1] ) ) );

      // Handle boundaries as special cases
      if ( B_index_i == 0 ) {
        B_bound[i] = -B_max;
      } else if ( B_index_i == B_count ) {
        B_bound[i] = B_max;
      } else {

        // Calculate the target area of unit disk
        const double A0 = A_bound[0];
        const double A1 = A_bound[1];
        const double target_area = ( A1*RE_SQRT( 1 - SQR( A1 ) ) - A0*RE_SQRT( 1 - SQR( A0 ) ) + asin( A1 ) - asin( A0 ) ) * ( ( double ) B_index_i ) / ( ( double ) B_count );

        // Set up GSL root solver
        double params[] = { target_area, A0, A1 };
        gsl_function F = { .function = EqualAreaSkyBoundSolverB, .params = params };
        double B_lower = -B_max, B_upper = B_max;
        XLAL_CHECK( gsl_root_fsolver_set( fs, &F, B_lower, B_upper ) == 0, XLAL_EFAILED );

        // Try to find root
        int status = 0, iter = 0;
        do {
          XLAL_CHECK( gsl_root_fsolver_iterate( fs ) == 0, XLAL_EFAILED );
          B_lower = gsl_root_fsolver_x_lower( fs );
          B_upper = gsl_root_fsolver_x_upper( fs );
          status = gsl_root_test_interval( B_lower, B_upper, 1e-5, 1e-5 );
        } while (status == GSL_CONTINUE && ++iter < 1000);
        XLAL_CHECK( status == GSL_SUCCESS, XLAL_EMAXITER, "Could not find bound for B_index[%zu]=%i; best guess [%g, %g]", i, B_index_i, B_lower, B_upper );

        // Store bound
        B_bound[i] = gsl_root_fsolver_root( fs );

      }

    }

    // Restrict range 'A' if 'B = const' bounds intersect unit disk boundary within it
    // - Only start to do this when there are 3 patches in 'B' direction
    // - Do not do this for the middle 'B' patch which straddles 'B = 0'
    if ( B_count >= 3 && B_index != (B_count-1)/2 ) {
      const double Ai = RE_SQRT( 1 - GSL_MIN( SQR( B_bound[0] ), SQR( B_bound[1] ) ) );
      if ( A_bound[0] < -Ai ) {
        A_bound[0] = -Ai;
      }
      if ( A_bound[1] > Ai ) {
        A_bound[1] = Ai;
      }
    }

    // Shift bounds on 'A' into the correct hemisphere, and put in correct order
    if ( patch_index < hemi_patch_count ) {    // This patch is in the left hemisphere
      A_bound[0] = A_bound[0] - 1;
      A_bound[1] = A_bound[1] - 1;
    } else {                                          // This patch is in the left hemisphere
      A_bound[0] = -A_bound[0] + 1;
      A_bound[1] = -A_bound[1] + 1;
    }

    // Cleanup
    gsl_root_fsolver_free( fs );

  }

  // Set the parameter-space bounds on reduced supersky sky coordinates A and B
  XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, 0, A_bound[0], A_bound[1] ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingBound( tiling, 1, ConstantBoundB, sizeof( B_bound[0] ), &B_bound[0], &B_bound[1] ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingBoundCacheFunction( tiling, 1, SkyBoundCache ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Set the parameter-space origin on reduced supersky sky coordinates A and B
  XLAL_CHECK( XLALSetLatticeTilingOrigin( tiling, 0, 0.0 ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingOrigin( tiling, 1, 0.0 ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Set the parameter-space padding control flags for reduced supersky sky coordinates A and B
  XLAL_CHECK( XLALSetLatticeTilingPaddingFlags( tiling, 0, A_padf ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingPaddingFlags( tiling, 1, B_padf ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

static double PhysicalSpinBound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_matrix *cache,
  const gsl_vector *point UNUSED
  )
{

  // Get bounds data
  const double *sky_offsets = ( ( const double * ) data );
  double bound = ( ( const double * ) data )[3];

  // Add the inner product of the sky offsets with the aligned sky
  // position to the physical bound to get the reduced supersky bound
  for ( size_t i = 0; i < 3; ++i ) {
    bound += sky_offsets[i] * gsl_matrix_get( cache, 1, i );
  }

  return bound;

}

int XLALSetSuperskyPhysicalSpinBound(
  LatticeTiling *tiling,
  const SuperskyTransformData *rssky_transf,
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
  XLAL_CHECK( s <= rssky_transf->SMAX, XLAL_ESIZE );

  // Decompose coordinate transform data
  gsl_matrix_const_view sky_offsets = gsl_matrix_const_view_array( &rssky_transf->sky_offsets[0][0], rssky_transf->nsky_offsets, 3 );

  // Set parameter-space bound name
  XLAL_CHECK( XLALSetLatticeTilingBoundName( tiling, RSSKY_FKDOT_DIM( rssky_transf, s ), "nu%zudot", s ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Copy the sky offset vector to bounds data
  double data_lower[4], data_upper[4];
  for ( size_t j = 0; j < 3; ++j ) {
    data_lower[j] = data_upper[j] = gsl_matrix_get( &sky_offsets.matrix, RSSKY_FKDOT_OFFSET( rssky_transf, s ), j );
  }
  data_lower[3] = GSL_MIN( bound1, bound2 );
  data_upper[3] = GSL_MAX( bound1, bound2 );

  // Set the parameter-space bound on physical frequency/spindown coordinate
  XLAL_CHECK( XLALSetLatticeTilingBound( tiling, RSSKY_FKDOT_DIM( rssky_transf, s ), PhysicalSpinBound, sizeof( data_lower ), &data_lower, &data_upper ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Supersky metric requires searching over spindown if also searching over sky
  const int tiled_sskyA = XLALIsTiledLatticeTilingDimension( tiling, 0 );
  XLAL_CHECK( tiled_sskyA != XLAL_FAILURE, XLAL_EFUNC );
  const int tiled_sskyB = XLALIsTiledLatticeTilingDimension( tiling, 1 );
  XLAL_CHECK( tiled_sskyB != XLAL_FAILURE, XLAL_EFUNC );
  const int tiled_fkdot = XLALIsTiledLatticeTilingDimension( tiling, RSSKY_FKDOT_DIM( rssky_transf, s ) );
  XLAL_CHECK( tiled_fkdot != XLAL_FAILURE, XLAL_EFUNC );
  XLAL_CHECK( !( tiled_sskyA || tiled_sskyB ) || tiled_fkdot, XLAL_EINVAL, "Must search over %zu-order spindown if also searching over sky", s );

  return XLAL_SUCCESS;

}

static int SM_LatticePhysicalRangeCallback(
  const bool first_call,
  const LatticeTiling *tiling,
  const LatticeTilingIterator *itr,
  const gsl_vector *point,
  const size_t changed_i,
  const void *param,
  void *out
  )
{

  // Only care about changes in sky position
  if ( 2 <= changed_i ) {
    return XLAL_SUCCESS;
  }

  // Get callback data
  const SM_CallbackParam *cparam = ( ( const SM_CallbackParam * ) param );
  SM_CallbackOut *cout = ( ( SM_CallbackOut * ) out );

  // Initialise translation data
  if ( first_call ) {
    XLAL_INIT_MEM( cout->min_phys );
    XLAL_INIT_MEM( cout->max_phys );
    cout->min_phys.refTime = cparam->rssky_transf->ref_time;
    cout->max_phys.refTime = cparam->rssky_transf->ref_time;
    cout->min_phys.Alpha = GSL_POSINF;
    cout->max_phys.Alpha = GSL_NEGINF;
    cout->min_phys.Delta = GSL_POSINF;
    cout->max_phys.Delta = GSL_NEGINF;
    for ( size_t s = 0; s <= cparam->rssky_transf->SMAX; ++s ) {
      cout->min_phys.fkdot[s] = GSL_POSINF;
      cout->max_phys.fkdot[s] = GSL_NEGINF;
    }
  }

  // Convert point from reduced supersky coordinates to physical coordinates
  PulsarDopplerParams XLAL_INIT_DECL( phys );
  XLAL_CHECK( XLALConvertSuperskyToPhysicalPoint( &phys, point, NULL, cparam->rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Store minimum/maximum values of physical sky position
  cout->min_phys.Alpha = GSL_MIN( cout->min_phys.Alpha, phys.Alpha );
  cout->max_phys.Alpha = GSL_MAX( cout->max_phys.Alpha, phys.Alpha );
  cout->min_phys.Delta = GSL_MIN( cout->min_phys.Delta, phys.Delta );
  cout->max_phys.Delta = GSL_MAX( cout->max_phys.Delta, phys.Delta );

  for ( size_t s = 0; s <= cparam->rssky_transf->SMAX; ++s ) {

    // Get indexes of left/right-most point in current frequency block
    INT4 left = 0, right = 0;
    XLAL_CHECK( XLALCurrentLatticeTilingBlock( itr, RSSKY_FKDOT_DIM( cparam->rssky_transf, s ), &left, &right ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Get step size
    const double step = XLALLatticeTilingStepSize( tiling, RSSKY_FKDOT_DIM( cparam->rssky_transf, s ) );

    // Store minimum/maximum values of physical frequency and spindowns
    if ( XLALIsTiledLatticeTilingDimension( tiling, RSSKY_FKDOT_DIM( cparam->rssky_transf, s ) ) ) {
      cout->min_phys.fkdot[s] = GSL_MIN( cout->min_phys.fkdot[s], phys.fkdot[s] + step*(left - 1) );
      cout->max_phys.fkdot[s] = GSL_MAX( cout->max_phys.fkdot[s], phys.fkdot[s] + step*(right + 1) );
    } else {
      cout->min_phys.fkdot[s] = cout->max_phys.fkdot[s] = phys.fkdot[s];
    }

  }

  return XLAL_SUCCESS;

}

int XLALRegisterSuperskyLatticePhysicalRangeCallback(
  LatticeTiling *tiling,
  const SuperskyTransformData *rssky_transf,
  const PulsarDopplerParams **min_phys,
  const PulsarDopplerParams **max_phys
  )
{

  // Check input
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( min_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( max_phys != NULL, XLAL_EFAULT );

  // Register callback function
  const SM_CallbackParam param = {
    .rssky_transf = rssky_transf,
  };
  const SM_CallbackOut *out = XLALRegisterLatticeTilingCallback( tiling, SM_LatticePhysicalRangeCallback, sizeof( param ), &param, sizeof( *out ) );
  XLAL_CHECK( out != NULL, XLAL_EFUNC );

  // Set output parameters
  *min_phys = &out->min_phys;
  *max_phys = &out->max_phys;

  return XLAL_SUCCESS;

}

static int SM_LatticeSuperskyRangeCallback(
  const bool first_call,
  const LatticeTiling *tiling,
  const LatticeTilingIterator *itr,
  const gsl_vector *point,
  const size_t changed_i,
  const void *param,
  void *out
  )
{

  // Only care about changes in sky position
  if ( 2 <= changed_i ) {
    return XLAL_SUCCESS;
  }

  // Get callback data
  const SM_CallbackParam *cparam = ( ( const SM_CallbackParam * ) param );
  SM_CallbackOut *cout = ( ( SM_CallbackOut * ) out );

  // Initialise translation data
  if ( first_call ) {
    cout->min_rssky2_view = gsl_vector_view_array( cout->min_rssky2_array, cparam->rssky2_transf->ndim );
    cout->max_rssky2_view = gsl_vector_view_array( cout->max_rssky2_array, cparam->rssky2_transf->ndim );
    gsl_vector_set_all( &cout->min_rssky2_view.vector, GSL_POSINF );
    gsl_vector_set_all( &cout->max_rssky2_view.vector, GSL_NEGINF );
  }

  // Convert point from reduced supersky coordinates to other reduced supersky coordinates
  double rssky2_array[cparam->rssky_transf->ndim];
  gsl_vector_view rssky2_view = gsl_vector_view_array( rssky2_array, cparam->rssky2_transf->ndim );
  XLAL_CHECK( XLALConvertSuperskyToSuperskyPoint( &rssky2_view.vector, cparam->rssky2_transf, point, point, cparam->rssky_transf ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Store minimum/maximum values of other reduced supersky sky coordinates
  for ( size_t i = 0; i < 2; ++i ) {
    cout->min_rssky2_array[i] = GSL_MIN( cout->min_rssky2_array[i], rssky2_array[i] );
    cout->max_rssky2_array[i] = GSL_MAX( cout->max_rssky2_array[i], rssky2_array[i] );
  }

  for ( size_t i = 2; i < cparam->rssky2_transf->ndim; ++i ) {

    // Get indexes of left/right-most point in current frequency block
    INT4 left = 0, right = 0;
    XLAL_CHECK( XLALCurrentLatticeTilingBlock( itr, i, &left, &right ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Get step size
    const double step = XLALLatticeTilingStepSize( tiling, i );

    // Store minimum/maximum values of other reduced supersky frequency and spindown coordinates
    if ( XLALIsTiledLatticeTilingDimension( tiling, i ) ) {
      cout->min_rssky2_array[i] = GSL_MIN( cout->min_rssky2_array[i], rssky2_array[i] + step*(left - 1) );
      cout->max_rssky2_array[i] = GSL_MAX( cout->max_rssky2_array[i], rssky2_array[i] + step*(right + 1) );
    } else {
      cout->min_rssky2_array[i] = cout->max_rssky2_array[i] = rssky2_array[i];
    }

  }

  return XLAL_SUCCESS;

}

int XLALRegisterSuperskyLatticeSuperskyRangeCallback(
  LatticeTiling *tiling,
  const SuperskyTransformData *rssky_transf,
  const SuperskyTransformData *rssky2_transf,
  const gsl_vector **min_rssky2,
  const gsl_vector **max_rssky2
  )
{

  // Check input
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky_transf ), XLAL_EINVAL );
  XLAL_CHECK( CHECK_RSSKY_TRANSF( rssky2_transf ), XLAL_EINVAL );
  XLAL_CHECK( min_rssky2 != NULL, XLAL_EFAULT );
  XLAL_CHECK( max_rssky2 != NULL, XLAL_EFAULT );

  // Register callback function
  const SM_CallbackParam param = {
    .rssky_transf = rssky_transf,
    .rssky2_transf = rssky2_transf,
  };
  const SM_CallbackOut *out = XLALRegisterLatticeTilingCallback( tiling, SM_LatticeSuperskyRangeCallback, sizeof( param ), &param, sizeof( *out ) );
  XLAL_CHECK( out != NULL, XLAL_EFUNC );
  XLAL_CHECK( rssky2_transf->ndim <= XLAL_NUM_ELEM( out->min_rssky2_array ), XLAL_EFAILED );
  XLAL_CHECK( rssky2_transf->ndim <= XLAL_NUM_ELEM( out->max_rssky2_array ), XLAL_EFAILED );

  // Set output parameters
  *min_rssky2 = &out->min_rssky2_view.vector;
  *max_rssky2 = &out->max_rssky2_view.vector;

  return XLAL_SUCCESS;

}

int XLALSetSuperskyRangeBounds(
  LatticeTiling *tiling,
  const gsl_vector *min_rssky,
  const gsl_vector *max_rssky
  )
{

  // Check input
  XLAL_CHECK( tiling != NULL, XLAL_EFAULT );
  const size_t n = XLALTotalLatticeTilingDimensions( tiling );
  XLAL_CHECK( min_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( min_rssky->size == n, XLAL_EINVAL );
  XLAL_CHECK( max_rssky != NULL, XLAL_EFAULT );
  XLAL_CHECK( max_rssky->size == n, XLAL_EINVAL );

  // Set the parameter-space bounds on reduced supersky sky coordinates A and B
  double A_bound[2] = { gsl_vector_get( min_rssky, 0 ), gsl_vector_get( max_rssky, 0 ) };
  double B_bound[2] = { gsl_vector_get( min_rssky, 1 ), gsl_vector_get( max_rssky, 1 ) };
  XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, 0, A_bound[0], A_bound[1] ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingBound( tiling, 1, ConstantBoundB, sizeof( B_bound[0] ), &B_bound[0], &B_bound[1] ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALSetLatticeTilingBoundCacheFunction( tiling, 1, SkyBoundCache ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Set the parameter-space bounds on all other coordinates
  for ( size_t j = 2; j < n; ++j ) {
    XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, j, gsl_vector_get( min_rssky, j ), gsl_vector_get( max_rssky, j ) ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
