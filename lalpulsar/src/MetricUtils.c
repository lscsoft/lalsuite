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
