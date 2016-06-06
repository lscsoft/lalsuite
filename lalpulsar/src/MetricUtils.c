//
// Copyright (C) 2011--2015 Reinhard Prix, Karl Wette
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <lal/MetricUtils.h>
#include <lal/XLALError.h>
#include <lal/LogPrintf.h>
#include <lal/Factorial.h>

#include <lal/GSLHelpers.h>

// Shortcuts for integer powers
#define POW2(a)  ( (a) * (a) )
#define POW3(a)  ( (a) * (a) * (a) )
#define POW4(a)  ( (a) * (a) * (a) * (a) )
#define POW5(a)  ( (a) * (a) * (a) * (a) * (a) )

///
/// Flexible comparison function between two metrics \f$ g^1_{ij} \f$ and \f$ g^2_{ij} \f$.
///
/// Returns maximal relative deviation \f$ \max_{i,j} \delta g_{ij} \f$, measured in terms of
/// diagonal entries, i.e. \f$ \delta g_{ij} = ( g^1_{ij} - g^2_{ij} ) / \sqrt{ g^1_{ii} g^1_{jj} } \f$.
/// This should always be well-defined, as we deal with positive-definite square matrices.
///
REAL8
XLALCompareMetrics(
  const gsl_matrix *g1_ij,			///< [in] Metric to compare \f$ g^1_{ij} \f$.
  const gsl_matrix *g2_ij			///< [in] Metric to compare \f$ g^2_{ij} \f$.
  )
{

  // Check input
  XLAL_CHECK_REAL8 ( g1_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK_REAL8 ( g2_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK_REAL8 ( g1_ij->size1 == g1_ij->size2, XLAL_ESIZE );
  XLAL_CHECK_REAL8 ( g2_ij->size1 == g2_ij->size2, XLAL_ESIZE );
  XLAL_CHECK_REAL8 ( g1_ij->size1 == g2_ij->size1, XLAL_ESIZE );

  if (lalDebugLevel & LALINFOBIT) {
    fprintf(stderr, "%s(): comparing this metric:", __func__);
    XLALfprintfGSLmatrix ( stderr, "%.15e", g1_ij );
    fprintf(stderr, "%s(): to this metric:", __func__);
    XLALfprintfGSLmatrix ( stderr, "%.15e", g2_ij );
  }

  REAL8 errmax = 0;
  UINT4 dim = g1_ij->size1;
  for ( UINT4 i = 0; i < dim; i ++ )
    {
      for ( UINT4 j = 0; j < dim; j ++ )
        {
          REAL8 norm = sqrt ( gsl_matrix_get ( g1_ij, i, i ) * gsl_matrix_get ( g1_ij, j, j ) );
          REAL8 e1 = gsl_matrix_get ( g1_ij, i, j ) / norm;
          REAL8 e2 = gsl_matrix_get ( g2_ij, i, j ) / norm;
          REAL8 base;
          if ( e2 == 0 )
            base = 1;
          else
            base = e2;
          REAL8 reldiff = fabs ( e1 - e2 ) / base;

          errmax = fmax ( errmax, reldiff );
        } // for j < dim
    } // for i < dim

  return errmax;

} // XLALCompareMetrics()

///
/// Compute the determinant of a metric \f$ g_{ij} \f$.
///
double XLALMetricDeterminant(
  const gsl_matrix *g_ij		///< [in] Parameter-space metric \f$ g_{ij} \f$
  )
{

  // Check input
  XLAL_CHECK_REAL8( g_ij != NULL, XLAL_EFAULT );

  const size_t n = g_ij->size1;

  // Allocate memory
  gsl_matrix *GAMAT_REAL8( LU_decomp, n, n );
  gsl_permutation *GAPERM_REAL8( LU_perm, n );

  // Compute LU decomposition
  gsl_matrix_memcpy( LU_decomp, g_ij );
  int LU_sign = 0;
  XLAL_CHECK_REAL8( gsl_linalg_LU_decomp( LU_decomp, LU_perm, &LU_sign ) == 0, XLAL_EFAILED, "'g_ij' cannot be LU-decomposed" );

  // Compute determinant
  const double determinant = gsl_linalg_LU_det( LU_decomp, LU_sign );

  // Cleanup
  GFMAT( LU_decomp );
  GFPERM( LU_perm );

  return determinant;

}

///
/// Compute the extent of the bounding box of the mismatch ellipse of a metric \f$ g_{ij} \f$.
///
gsl_vector *XLALMetricEllipseBoundingBox(
  const gsl_matrix *g_ij,		///< [in] Parameter-space metric \f$ g_{ij} \f$
  const double max_mismatch		///< [in] Maximum prescribed mismatch
  )
{

  // Check input
  XLAL_CHECK_NULL( g_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( g_ij->size1 == g_ij->size2, XLAL_ESIZE );

  const size_t n = g_ij->size1;

  // Allocate memory
  gsl_matrix *GAMAT_NULL( LU_decomp, n, n );
  gsl_permutation *GAPERM_NULL( LU_perm, n );
  gsl_matrix *GAMAT_NULL( diag_norm, n, n );
  gsl_matrix *GAMAT_NULL( inverse, n, n );
  gsl_vector *GAVEC_NULL( bounding_box, n );

  // Diagonally normalize metric
  XLAL_CHECK_NULL( XLALDiagNormalizeMetric( &LU_decomp, &diag_norm, g_ij ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Compute metric inverse
  int LU_sign = 0;
  XLAL_CHECK_NULL( gsl_linalg_LU_decomp( LU_decomp, LU_perm, &LU_sign ) == 0, XLAL_EFAILED, "'g_ij' cannot be LU-decomposed" );
  XLAL_CHECK_NULL( gsl_linalg_LU_invert( LU_decomp, LU_perm, inverse ) == 0, XLAL_EFAILED, "'g_ij' cannot be inverted" );

  // Compute bounding box, and invert diagonal normalization
  for( size_t i = 0; i < n; ++i ) {
    const double diag_norm_i = gsl_matrix_get( diag_norm, i, i );
    const double bounding_box_i = 2.0 * sqrt( max_mismatch * gsl_matrix_get( inverse, i, i ) ) * diag_norm_i;
    gsl_vector_set( bounding_box, i, bounding_box_i );
  }

  // Cleanup
  GFMAT( LU_decomp, diag_norm, inverse );
  GFPERM( LU_perm );

  return bounding_box;

} // XLALMetricEllipseBoundingBox()

///
/// Calculate the projected metric onto the subspace orthogonal to coordinate-axis \c c, namely
/// \f$ g^{\prime}_{ij} = g_{ij} - ( g_{ic} g_{jc} / g_{cc} ) \f$, where \f$c\f$ is the index of
/// the projected coordinate.
///
/// \note \c *p_gpr_ij will be allocated if \c NULL. \c *p_gpr_ij and \c g_ij may point to the same
/// matrix.
///
int
XLALProjectMetric(
  gsl_matrix **p_gpr_ij,			///< [in,out] Pointer to projected matrix \f$g^{\prime}_{ij}\f$
  const gsl_matrix *g_ij,			///< [in] Matrix to project \f$g_{ij}\f$
  const UINT4 c					///< [in] Index of projected coordinate
  )
{

  // Check input
  XLAL_CHECK( g_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK( g_ij->size1 == g_ij->size2, XLAL_ESIZE );
  XLAL_CHECK( p_gpr_ij != NULL, XLAL_EFAULT );
  if ( *p_gpr_ij != NULL ) {
    XLAL_CHECK( (*p_gpr_ij)->size1 == (*p_gpr_ij)->size2, XLAL_ESIZE );
    XLAL_CHECK( (*p_gpr_ij)->size1 == g_ij->size1, XLAL_ESIZE );
  } else {
    GAMAT( *p_gpr_ij, g_ij->size1, g_ij->size2 );
  }
  XLAL_CHECK( c < g_ij->size1, XLAL_EINVAL );

  // Allocate temporary matrix
  gsl_matrix *GAMAT( ret_ij, g_ij->size1, g_ij->size2 );

  // Compute projected matrix
  for ( UINT4 i=0; i < g_ij->size1; i++)
    {
    for ( UINT4 j=0; j < g_ij->size2; j++ )
      {
        if ( i==c || j==c )
          {
            gsl_matrix_set ( ret_ij, i, j, 0.0 );
          }
        else
          {
            double proj = gsl_matrix_get(g_ij, i, j) - (gsl_matrix_get(g_ij, i, c) * gsl_matrix_get(g_ij, j, c) / gsl_matrix_get(g_ij, c, c));
            gsl_matrix_set ( ret_ij, i, j, proj );
          }
      } /* for j < dim2 */

    } /* for i < dim1 */
  gsl_matrix_memcpy( *p_gpr_ij, ret_ij );

  // Cleanup
  GFMAT( ret_ij );

  return XLAL_SUCCESS;

} // XLALProjectMetric()

///
/// Decompose a metric \f$\mathbf{G}\f$ as \f$ \mathbf{G} = \mathbf{L} \mathbf{D}
/// \mathbf{L}^{\mathrm{T}} \f$, where \f$\mathbf{L}\f$ is a lower-triangular matrix
/// with unit diagonal, and \f$\mathbf{D}\f$ is a diagonal matrix. This decomposition
/// may be useful if the metric cannot yet be guaranteed to be positive definite.
///
int XLALCholeskyLDLTDecompMetric(
  gsl_matrix **p_cholesky,			///< [in,out] Pointer to decomposition; stores \f$\mathbf{L}\f$ in lower triangular part \f$\mathbf{D}\f$ on diagonal
  const gsl_matrix *g_ij			///< [in] Matrix to decompose \f$\mathbf{G}\f$
  )
{


  // Check input
  XLAL_CHECK( g_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK( g_ij->size1 == g_ij->size2, XLAL_ESIZE );
  XLAL_CHECK( p_cholesky != NULL, XLAL_EFAULT );
  if ( *p_cholesky != NULL ) {
    XLAL_CHECK( (*p_cholesky)->size1 == g_ij->size1, XLAL_ESIZE );
    XLAL_CHECK( (*p_cholesky)->size2 == g_ij->size2, XLAL_ESIZE );
  } else {
    GAMAT( *p_cholesky, g_ij->size1, g_ij->size2 );
  }

  // Straightforward implementation of Cholesky–Banachiewicz algorithm
  gsl_matrix_set_zero( *p_cholesky );
#define A(i,j) *gsl_matrix_const_ptr( g_ij, i, j )
#define D(j)   *gsl_matrix_ptr( *p_cholesky, j, j )
#define L(i,j) *gsl_matrix_ptr( *p_cholesky, i, j )
  for ( size_t i = 0; i < g_ij->size1; ++i ) {
    for ( size_t j = 0; j <= i; ++j ) {
      if ( i == j ) {
        D(j) = A(j, j);
        for ( size_t k = 0; k < j; ++k ) {
          D(j) -= L(j, k) * L(j, k) * D(k);
        }
      } else {
        L(i, j) = A(i, j);
        for ( size_t k = 0; k < j; ++k ) {
          L(i, j) -= L(i, k) * L(j, k) * D(k);
        }
        L(i, j) /= D(j);
      }
    }
  }
#undef A
#undef D
#undef L

  return XLAL_SUCCESS;

}

///
/// Apply a transform \f$\mathbf{A} = (a_{ij})\f$ to a metric \f$\mathbf{G} = (g_{ij})\f$ such that
/// \f$ \mathbf{G} \rightarrow \mathbf{G}^{\prime} = \mathbf{A}^{\mathrm{T}} \mathbf{G} \mathbf{A}
/// \f$, or equivalently \f$ g_{ij} \rightarrow g^{\prime}_{kl} = g_{ij} a_{ik} a_{jl} \f$.
///
/// \note \c *p_gpr_ij will be allocated if \c NULL. \c *p_gpr_ij and \c g_ij may point to the same
/// matrix.
///
int XLALTransformMetric(
  gsl_matrix **p_gpr_ij,			///< [in,out] Pointer to transformed matrix \f$\mathbf{G}^{\prime}\f$
  const gsl_matrix *transform,			///< [in] Transform to apply \f$\mathbf{A}\f$
  const gsl_matrix *g_ij			///< [in] Matrix to transform \f$\mathbf{G}\f$
  )
{

  // Check input
  XLAL_CHECK( g_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK( transform != NULL, XLAL_EFAULT );
  XLAL_CHECK( g_ij->size1 == g_ij->size2, XLAL_ESIZE );
  XLAL_CHECK( g_ij->size2 == transform->size1, XLAL_ESIZE );
  XLAL_CHECK( p_gpr_ij != NULL, XLAL_EFAULT );
  if ( *p_gpr_ij != NULL ) {
    XLAL_CHECK( (*p_gpr_ij)->size1 == transform->size2, XLAL_ESIZE );
    XLAL_CHECK( (*p_gpr_ij)->size2 == transform->size2, XLAL_ESIZE );
  } else {
    GAMAT( *p_gpr_ij, transform->size2, transform->size2 );
  }

  // Allocate temporary matrix
  gsl_matrix *tmp = gsl_matrix_alloc( g_ij->size1, transform->size2 );
  XLAL_CHECK( tmp != NULL, XLAL_ENOMEM );

  // Perform transform
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, g_ij, transform, 0.0, tmp );
  gsl_blas_dgemm( CblasTrans, CblasNoTrans, 1.0, transform, tmp, 0.0, *p_gpr_ij );

  // Ensure transformed g_ij is exactly symmetric
  for( size_t i = 0; i < (*p_gpr_ij)->size1; ++i ) {
    for( size_t j = i + 1; j < (*p_gpr_ij)->size2; ++j ) {
      const double gij = gsl_matrix_get( *p_gpr_ij, i, j );
      const double gji = gsl_matrix_get( *p_gpr_ij, j, i );
      const double g = 0.5 * ( gij + gji );
      gsl_matrix_set( *p_gpr_ij, i, j, g );
      gsl_matrix_set( *p_gpr_ij, j, i, g );
    }
  }

  // Cleanup
  gsl_matrix_free( tmp );

  return XLAL_SUCCESS;

} // XLALTransformMetric()

///
/// Apply the inverse of a transform \f$\mathbf{A}^{-1} = \mathbf{B} = (b_{ij})\f$ to a metric
/// \f$\mathbf{G} = (g_{ij})\f$ such that \f$ \mathbf{G} \rightarrow \mathbf{G}^{\prime} =
/// \mathbf{B}^{\mathrm{T}} \mathbf{G} \mathbf{B} \f$, or equivalently \f$ g_{ij} \rightarrow
/// g^{\prime}_{kl} = g_{ij} b_{ik} b_{jl} \f$.
///
/// \note \c *p_gpr_ij will be allocated if \c NULL. \c *p_gpr_ij and \c g_ij may point to the same
/// matrix.
///
int XLALInverseTransformMetric(
  gsl_matrix **p_gpr_ij,			///< [in,out] Pointer to transformed matrix \f$\mathbf{G}^{\prime}\f$
  const gsl_matrix *transform,			///< [in] Transform \f$\mathbf{A}\f$, the inverse of which to apply
  const gsl_matrix *g_ij			///< [in] Matrix to transform \f$\mathbf{G}\f$
  )
{

  // Check input
  XLAL_CHECK( g_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK( transform != NULL, XLAL_EFAULT );
  XLAL_CHECK( g_ij->size1 == g_ij->size2, XLAL_ESIZE );
  XLAL_CHECK( g_ij->size2 == transform->size1, XLAL_ESIZE );
  XLAL_CHECK( transform->size1 == transform->size2, XLAL_ESIZE );
  XLAL_CHECK( p_gpr_ij != NULL, XLAL_EFAULT );

  // Allocate memory
  gsl_matrix *GAMAT( LU_decomp, transform->size1, transform->size2 );
  gsl_permutation *GAPERM( LU_perm, transform->size1 );
  gsl_matrix *GAMAT( inverse, transform->size1, transform->size2 );

  // Compute inverse of transform
  int LU_sign = 0;
  gsl_matrix_memcpy( LU_decomp, transform );
  XLAL_CHECK( gsl_linalg_LU_decomp( LU_decomp, LU_perm, &LU_sign ) == 0, XLAL_EFAILED, "'transform' cannot be LU-decomposed" );
  XLAL_CHECK( gsl_linalg_LU_invert( LU_decomp, LU_perm, inverse ) == 0, XLAL_EFAILED, "'transform' cannot be inverted" );

  // Apply transform
  XLAL_CHECK( XLALTransformMetric( p_gpr_ij, inverse, g_ij ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Cleanup
  GFMAT( LU_decomp, inverse );
  GFPERM( LU_perm );

  return XLAL_SUCCESS;

} // XLALInverseTransformMetric()

///
/// Diagonally-normalize a metric \f$ g_{ij} \f$.  <i>Diagonally-normalize</i> means normalize
/// metric by its diagonal, namely apply the transformation \f$ g_{ij} \rightarrow g^{\prime}_{ij} =
/// g_{ij} / \sqrt{g_{ii} g_{jj}} \f$, resulting in a lower condition number and unit diagonal
/// elements.
///
/// If \c p_transform is non-NULL, return the diagonal-normalization transform in \c *p_transform.
/// If \c p_gpr_ij is non-NULL, apply the transform to the metric \c *p_gpr_ij.
///
/// \note \c *p_gpr_ij and \c *p_transform will be allocated if \c NULL. \c *p_gpr_ij and \c g_ij
/// may point to the same matrix.
///
int
XLALDiagNormalizeMetric(
  gsl_matrix **p_gpr_ij,			///< [in,out,optional] Pointer to transformed matrix \f$g^{\prime}_{ij}\f$
  gsl_matrix **p_transform,			///< [in,out,optional] Pointer to diagonal-normalization transform
  const gsl_matrix *g_ij			///< [in] Matrix to transform \f$g_{ij}\f$
  )
{

  // Check input
  XLAL_CHECK( g_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK( g_ij->size1 == g_ij->size2, XLAL_ESIZE );
  if ( p_transform != NULL ) {
    if ( *p_transform != NULL ) {
      XLAL_CHECK( (*p_transform)->size1 == g_ij->size1, XLAL_ESIZE );
      XLAL_CHECK( (*p_transform)->size2 == g_ij->size2, XLAL_ESIZE );
    } else {
      GAMAT( *p_transform, g_ij->size1, g_ij->size2 );
    }
  }

  // Create diagonal normalization transform
  gsl_matrix *GAMAT( transform, g_ij->size1, g_ij->size2 );
  gsl_matrix_set_zero( transform );
  for ( size_t i = 0; i < g_ij->size1; ++i ) {
    const double gii = gsl_matrix_get(g_ij, i, i);
    XLAL_CHECK( gii > 0, XLAL_EINVAL, "Diagonal normalization not defined for non-positive diagonal elements! g_ii(%zu,%zu) = %g\n", i, i, gii );
    gsl_matrix_set( transform, i, i, 1.0 / sqrt( gii ) );
  }

  // Apply transform
  if ( p_gpr_ij != NULL ) {
    XLAL_CHECK( XLALTransformMetric( p_gpr_ij, transform, g_ij ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Return transform
  if ( p_transform != NULL ) {
    gsl_matrix_memcpy( *p_transform, transform );
  }

  // Cleanup
  GFMAT( transform );

  return XLAL_SUCCESS;

} // XLALDiagNormalizeMetric()

///
/// Return a metric in <i>naturalized</i> coordinates.
/// Frequency coordinates of spindown order \f$s\f$ are scaled by
/// \f[ \frac{2\pi}{(s+1)!} \left(\frac{\overline{\Delta T}}{2}\right)^{s+1} \f]
/// where \f$\overline{\Delta T}\equiv\sum_{k}^{N} \Delta T_k\f$ is the average segment-length
/// over all \f$N\f$ segments.
///
/// Sky coordinates are scaled by Holgers' units, see Eq.(44) in PRD82,042002(2010),
/// without the equatorial rotation in alpha:
/// \f[ \frac{2\pi f R_{E}  \cos(\delta_{D})}{c} \f]
/// where \f$f\f$ is the frequency and
/// \f$R_{E}\f$ the Earth radius, and \f$\delta_{D}\f$ is the detectors latitude.
///
/// If \c p_transform is non-NULL, return the naturalization transform in \c *p_transform.
/// If \c p_gpr_ij is non-NULL, apply the transform to the metric \c *p_gpr_ij.
///
/// \note \c *p_gpr_ij and \c *p_transform will be allocated if \c NULL. \c *p_gpr_ij and \c g_ij
/// may point to the same matrix.
///
int
XLALNaturalizeMetric(
  gsl_matrix **p_gpr_ij,			///< [in,out,optional] Pointer to transformed matrix \f$g^{\prime}_{ij}\f$
  gsl_matrix **p_transform,			///< [in,out,optional] Pointer to naturalization transform
  const gsl_matrix *g_ij,			///< [in] Matrix to transform \f$g_{ij}\f$
  const DopplerMetricParams *metricParams	///< [in] Input parameters used to calculate naturalization transform
  )
{

  // Check input
  XLAL_CHECK( g_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK( g_ij->size1 == g_ij->size2, XLAL_ESIZE );
  if ( p_transform != NULL ) {
    if ( *p_transform != NULL ) {
      XLAL_CHECK( (*p_transform)->size1 == g_ij->size1, XLAL_ESIZE );
      XLAL_CHECK( (*p_transform)->size2 == g_ij->size2, XLAL_ESIZE );
    } else {
      GAMAT( *p_transform, g_ij->size1, g_ij->size2 );
    }
  }
  XLAL_CHECK( metricParams != NULL, XLAL_EFAULT );

  // Compute average segment length
  XLAL_CHECK ( XLALSegListIsInitialized ( &(metricParams->segmentList) ), XLAL_EINVAL, "Passed un-initialzied segment list 'metricParams->segmentList'\n");
  UINT4 Nseg = metricParams->segmentList.length;
  REAL8 sumTseg = 0;
  for ( UINT4 k=0; k < Nseg; k ++ )
    {
      LIGOTimeGPS *startTime_k = &(metricParams->segmentList.segs[k].start);
      LIGOTimeGPS *endTime_k   = &(metricParams->segmentList.segs[k].end);
      sumTseg += XLALGPSDiff( endTime_k, startTime_k );
    }
  REAL8 avgTseg = sumTseg / Nseg;

  // Compute naturalization transform
  gsl_matrix *GAMAT( transform, g_ij->size1, g_ij->size2 );
  gsl_matrix_set_zero( transform );
  for (size_t i = 0; i < g_ij->size1; ++i) {
    const DopplerCoordinateID coordID = metricParams->coordSys.coordIDs[i];
    const double Freq = metricParams->signalParams.Doppler.fkdot[0];
    const double T = avgTseg;
    double scale = 0;
    switch (coordID) {
    case DOPPLERCOORD_NONE:
      scale = 1;
      break;

    case DOPPLERCOORD_FREQ:
    case DOPPLERCOORD_GC_NU0:
      scale = LAL_TWOPI * LAL_FACT_INV[1] * (0.5 * T);
      break;

    case DOPPLERCOORD_F1DOT:
    case DOPPLERCOORD_GC_NU1:
      scale = LAL_TWOPI * LAL_FACT_INV[2] * POW2(0.5 * T);
      break;

    case DOPPLERCOORD_F2DOT:
    case DOPPLERCOORD_GC_NU2:
      scale = LAL_TWOPI * LAL_FACT_INV[3] * POW3(0.5 * T);
      break;

    case DOPPLERCOORD_F3DOT:
    case DOPPLERCOORD_GC_NU3:
      scale = LAL_TWOPI * LAL_FACT_INV[4] * POW4(0.5 * T);
      break;

    case DOPPLERCOORD_ALPHA:
    case DOPPLERCOORD_DELTA:
    case DOPPLERCOORD_N2X_EQU:
    case DOPPLERCOORD_N2Y_EQU:
    case DOPPLERCOORD_N2X_ECL:
    case DOPPLERCOORD_N2Y_ECL:
    case DOPPLERCOORD_N3X_EQU:
    case DOPPLERCOORD_N3Y_EQU:
    case DOPPLERCOORD_N3Z_EQU:
    case DOPPLERCOORD_N3X_ECL:
    case DOPPLERCOORD_N3Y_ECL:
    case DOPPLERCOORD_N3Z_ECL:
    case DOPPLERCOORD_N3SX_EQU:
    case DOPPLERCOORD_N3SY_EQU:
    case DOPPLERCOORD_N3OX_ECL:
    case DOPPLERCOORD_N3OY_ECL:
      {
        const REAL8 rEarth_c = (LAL_REARTH_SI / LAL_C_SI);
        REAL8 cosdD = cos ( metricParams->multiIFO.sites[0].frDetector.vertexLatitudeRadians );
        scale = LAL_TWOPI * Freq * rEarth_c * cosdD;
      }
      break;

    default:
      scale = 1;
      break;
    } // switch(coordID)
    gsl_matrix_set( transform, i, i, 1.0 / scale );
  }

  // Apply transform
  if ( p_gpr_ij != NULL ) {
    XLAL_CHECK( XLALTransformMetric( p_gpr_ij, transform, g_ij ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Return transform
  if ( p_transform != NULL ) {
    gsl_matrix_memcpy( *p_transform, transform );
  }

  // Cleanup
  GFMAT( transform );

  return XLAL_SUCCESS;

} // XLALNaturalizeMetric()

///
/// Compute the transform which changes the metric reference time \f$ \tau_0 \rightarrow \tau_1 =
/// \tau_0 + \Delta\tau \f$.
///
/// If \c p_transform is non-NULL, return the reference-time transform in \c *p_transform.
/// If \c p_gpr_ij is non-NULL, apply the transform to the metric \c *p_gpr_ij.
///
/// \note \c *p_gpr_ij and \c *p_transform will be allocated if \c NULL. \c *p_gpr_ij and \c g_ij
/// may point to the same matrix.
///
int
XLALChangeMetricReferenceTime(
  gsl_matrix **p_gpr_ij,			///< [in,out,optional] Pointer to transformed matrix \f$g^{\prime}_{ij}\f$
  gsl_matrix **p_transform,			///< [in,out,optional] Pointer to reference time transform
  const gsl_matrix *g_ij,			///< [in] Matrix to transform \f$g_{ij}\f$
  const DopplerCoordinateSystem *coordSys,	///< [in] Coordinate system of metric
  const double Dtau				///< [in] Difference between new and old reference times \f$\Delta\tau\f$
  )
{

  // Check input
  XLAL_CHECK( g_ij != NULL, XLAL_EFAULT );
  XLAL_CHECK( g_ij->size1 == g_ij->size2, XLAL_ESIZE );
  if ( p_transform != NULL ) {
    if ( *p_transform != NULL ) {
      XLAL_CHECK( (*p_transform)->size1 == g_ij->size1, XLAL_ESIZE );
      XLAL_CHECK( (*p_transform)->size2 == g_ij->size2, XLAL_ESIZE );
    } else {
      GAMAT( *p_transform, g_ij->size1, g_ij->size2 );
    }
  }
  XLAL_CHECK( coordSys != NULL, XLAL_EFAULT );
  XLAL_CHECK( coordSys->dim == g_ij->size1, XLAL_ESIZE );
  for( size_t i = 0; i < coordSys->dim; ++i ) {
    switch( coordSys->coordIDs[i] ) {
    case DOPPLERCOORD_GC_NU0:
    case DOPPLERCOORD_GC_NU1:
    case DOPPLERCOORD_GC_NU2:
    case DOPPLERCOORD_GC_NU3:
      XLAL_ERROR( XLAL_EINVAL, "GCT coordinates are not supported" );
    default:
      break;
    }
  }

  // Compute transformation matrix for frequency-spindown coordinates
  // from Prix: "Frequency metric for CW searches" (2014-08-17), p. 4
  gsl_matrix *GAMAT( transform, g_ij->size1, g_ij->size2 );
  gsl_matrix_set_identity( transform );
  if( Dtau != 0.0 ) {
    for( size_t i = 0; i < coordSys->dim; ++i ) {
      const DopplerCoordinateID icoord = coordSys->coordIDs[i];
      if (DOPPLERCOORD_FREQ <= icoord && icoord <= DOPPLERCOORD_LASTFDOT) {
        const size_t ispinorder = icoord - DOPPLERCOORD_FREQ;

        for( size_t j = i + 1; j < coordSys->dim; ++j ) {
          const DopplerCoordinateID jcoord = coordSys->coordIDs[j];
          if (DOPPLERCOORD_FREQ <= jcoord && jcoord <= DOPPLERCOORD_LASTFDOT) {
            const size_t jspinorder = jcoord - DOPPLERCOORD_FREQ;

            if( jspinorder > ispinorder ) {
              const double tij = LAL_FACT_INV[jspinorder - ispinorder] * pow( -1 * Dtau, jspinorder - ispinorder );
              gsl_matrix_set( transform, i, j, tij );
            }

          }
        }

      }
    }
  }

  // Apply transform
  if ( p_gpr_ij != NULL ) {
    XLAL_CHECK( XLALTransformMetric( p_gpr_ij, transform, g_ij ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Return transform
  if ( p_transform != NULL ) {
    gsl_matrix_memcpy( *p_transform, transform );
  }

  // Cleanup
  GFMAT( transform );

  return XLAL_SUCCESS;

} // XLALPhaseMetricRefTimeTransform()
