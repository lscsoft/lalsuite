//
// Copyright (C) 2019--2023 Benjamin Grace
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

#ifndef _PIECEWISEMODEL_H
#define _PIECEWISEMODEL_H

#include <stdio.h>
#include <math.h>
#include <lal/LatticeTiling.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Sets the bounds for the piecewise model
///
int XLALSetLatticeTilingPiecewiseBounds(
  LatticeTiling* tiling,
  const double fmin,            /// Minimum spin frequency to search over
  const double fmax,            /// Maximum spin frequency to search over
  const double fmaxtrue,        /// Maximum spin frequency used to calculate k and knots (useful for computing tiles in parrallel and fmax != fmaxtrue)
  const double nmin,            /// Minimum braking index
  const double nmax,            /// Maximum braking index
  const double nmin0,           /// Minimum braking index for the first knot. Useful if you want to brake up a search into partitions separated by templates with braking indices within a certain range
  const double nmax0,           /// Maximum braking index for the first knot. Useful if you want to brake up a search into partitions separated by templates with braking indices within a certain range
  const double ntol,            /// Tolerance (percentage per second) between braking indices on adjacent knots
  const double taumin,          /// Minimum spin half life when n = nmax, f0 = fmaxtrue
  const double taumax,          /// Maximum spin half life when n = nmax, f0 = fmaxtrue
  const double ktol,            /// Tolerance (percentage per second) between k values on adjacent knots
  const gsl_vector* knots,      /// List of knots
  const int finalknot,          /// The number of the final knot
  const gsl_vector* bboxpad,    /// Vector containing fractional bounding box padding
  const gsl_vector_int* intpad, ///< Vector containing number of integer points to use for padding
  const int reset               ///< 1 to use the resetting methods, -1 to not. If the bboxpad and intpad vectors contain non-zero elements, then this will automatically be set to 1
  );

///
/// Sets the bounds for the piecewise model when we are using 2 spin down parameters for each knot
///
int XLALSetLatticeTilingPiecewiseBoundsS2(
  LatticeTiling* tiling,
  const double fmin,            /// Minimum spin frequency to search over
  const double fmax,            /// Maximum spin frequency to search over
  const double fmaxtrue,        /// Maximum spin frequency used to calculate k and knots (useful for computing tiles in parrallel and fmax != fmaxtrue)
  const double nmin,            /// Minimum braking index
  const double nmax,            /// Maximum braking index
  const double nmin0,           /// Minimum braking index for the first knot. Useful if you want to brake up a search into partitions separated by templates with braking indices within a certain range
  const double nmax0,           /// Maximum braking index for the first knot. Useful if you want to brake up a search into partitions separated by templates with braking indices within a certain range
  const double ntol,            /// Tolerance (percentage per second) between braking indices on adjacent knots
  const double taumin,          /// Minimum spin half life when n = nmax, f0 = fmaxtrue
  const double taumax,          /// Maximum spin half life when n = nmax, f0 = fmaxtrue
  const double ktol,            /// Tolerance (percentage per second) between k values on adjacent knots
  const gsl_vector* knots,      /// List of knots
  const int finalknot,          /// The number of the final knot
  const gsl_vector* bboxpad,    /// Vector containing fractional bounding box padding
  const gsl_vector_int* intpad  /// Vector containing number of integer points to use for padding
  );


///
/// Returns the upper or lower bound for the dimension 'dim' given values for all previous dimensions, point_up_to_dim, as well as the relevant parameter space information
///
double XLALPiecewiseParameterBounds(
  const size_t      dim,              /// The dimension of the parameter we wish to calculate the bounds for
  const gsl_vector* point_up_to_dim,  /// The point/vector containing all parameter values up to the given dimension
  const int         upperlower,       /// +1 to return upper bound, -1 to return lower bound
  const double      fmin,             /// Global maximum frequency
  const double      fmax,             /// Global minimum frequency
  const double      nmin,             /// Minimum braking index
  const double      nmax,             /// Maximum braking index
  const double      ntol,             /// Tolerance (percentage per second) between braking indices on adjacent knots
  const double      kmin,             /// Minimum k value
  const double      kmax,             /// Maximum k value
  const double      ktol,             /// Tolerance (percentage per second) between k values on adjacent knots
  const double      segment_length,   /// The length of the segment. The time between the previous knot and the knot that the parameter 'dim' resides on
  const gsl_vector* knots,            /// The knots of the piecewise model
  const double      reset             /// 1 to use the resetting methods, -1 to not.
  );


#ifdef __cplusplus
}
#endif

#endif //_PIECEWISEMODEL_H

