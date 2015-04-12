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

#include "GSLHelpers.h"

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
  gsl_matrix *GAMAT_NULL( inverse, n, n );
  gsl_vector *GAVEC_NULL( bounding_box, n );

  // Copy metric, and ensure it is diagonally normalised
  for( size_t i = 0; i < n; ++i ) {
    const double norm_i = gsl_matrix_get( g_ij, i, i );
    for( size_t j = 0; j < n; ++j ) {
      const double norm_j = gsl_matrix_get( g_ij, j, j );
      gsl_matrix_set( LU_decomp, i, j, gsl_matrix_get( g_ij, i, j ) / sqrt( norm_i * norm_j ) );
    }
  }

  // Compute metric inverse
  int LU_sign = 0;
  GCALL_NULL( gsl_linalg_LU_decomp( LU_decomp, LU_perm, &LU_sign ), "'g_ij' cannot be LU-decomposed" );
  GCALL_NULL( gsl_linalg_LU_invert( LU_decomp, LU_perm, inverse ), "'g_ij' cannot be inverted" );

  // Compute bounding box, and reverse diagonal scaling
  for( size_t i = 0; i < n; ++i ) {
    const double norm_i = gsl_matrix_get( g_ij, i, i );
    const double bounding_box_i = 2.0 * sqrt( max_mismatch * gsl_matrix_get( inverse, i, i ) / norm_i );
    gsl_vector_set( bounding_box, i, bounding_box_i );
  }

  // Cleanup
  GFMAT( LU_decomp, inverse );
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
