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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301 USA
//

#include <stdio.h>
#include <math.h>
#include <lal/LatticeTiling.h>
#include <lal/LogPrintf.h>
#include <lal/PiecewiseModel.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

///
/// A struct containing the relevant information for calculating upper and lower bounds on a specific piecewise segment
///
typedef struct tagPiecewiseBoundInfo {
  size_t s;                   ///< Number of frequency/spindown parameters per knot
  double n;                   ///< Braking index
  double k;                   ///< k value
  double dt;                  ///< Time between knots
} PiecewiseBoundInfo;

///
/// The general torque equation and its first two derivatives
///
static double GTEAndDerivs(
  double f0,               ///< Initial frequency
  double n,                ///< Braking index
  double k,                ///< k value
  double t,                ///< Time at which to evaluate the GTE
  int d                    ///< Derivative order (d <= 2)
)
{
  double base = 1 + ( n - 1 ) * k * pow( f0, n - 1 ) * t;
  double power = 1 / ( 1 - n );

  if ( d == 0 ) {
    return f0 * pow( base, power );
  } else if ( d == 1 ) {
    double factor = -1 * pow( f0, n ) * k;
    return factor * pow( base, power - 1 );
  } else if ( d == 2 ) {
    double factor = pow( f0, 2 * n - 1 ) * k * k * n;
    double secondderiv = factor * pow( base, power - 2 );
    return secondderiv;
  }

  return NAN;
}

///
/// Sets the bound on the frequency parameter
///
static double F0Bound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_matrix *cache UNUSED,
  const gsl_vector *point
)
{
  const PiecewiseBoundInfo *info = ( const PiecewiseBoundInfo * )data;

  // Frequency at previous knot
  const double f0 = gsl_vector_get( point, dim - info->s );

  // Frequency at this knot, propagated by dt
  const double f = GTEAndDerivs( f0, info->n, info->k, info->dt, 0 );

  if ( f >= f0 ) {
    printf( "f >= f0, %g >= %g", f, f0 );
  }

  return f;
}

///
/// Sets the bound on the first derivative frequency parameter
///
static double F1Bound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_matrix *cache UNUSED,
  const gsl_vector *point
)
{
  const PiecewiseBoundInfo *info = ( const PiecewiseBoundInfo * )data;

  // Frequency at this knot
  const double f0 = gsl_vector_get( point, dim - 1 );

  // First derivative of the general torque equation
  // - Evaluated at dt=0, i.e. at this knot
  const double fd = GTEAndDerivs( f0, info->n, info->k, 0, 1 );

  return fd;
}

///
/// Sets the bounds on the second derivative frequency parameter
///
static double F2Bound(
  const void *data,
  const size_t dim UNUSED,
  const gsl_matrix *cache UNUSED,
  const gsl_vector *point
)
{
  const PiecewiseBoundInfo *info = ( const PiecewiseBoundInfo * )data;

  // Frequency at this knot
  const double f0 = gsl_vector_get( point, dim - 2 );

  // Second derivative of the general torque equation
  // - Evaluated at dt=0, i.e. at this knot
  const double fdd = GTEAndDerivs( f0, info->n, info->k, 0, 2 );

  return fdd;
}

///
/// Sets the bounds for the piecewise model
///
int XLALSetLatticeTilingPiecewiseBounds(
  LatticeTiling *tiling,        ///< Lattice tiling
  const size_t s,               ///< Number of frequency/spindown parameters per knot
  const double fmin,            ///< Minimum initial frequency
  const double fmax,            ///< Maximum initial frequency
  const double nmin,            ///< Minimum braking index
  const double nmax,            ///< Maximum braking index
  const double kmin,            ///< Minimum k value
  const double kmax,            ///< Maximum k value
  const gsl_vector *knots,      ///< List of knots
  const gsl_vector *bboxpad,    ///< Vector containing fractional bounding box padding
  const gsl_vector_int *intpad  ///< Vector containing number of integer points to use for padding
)
{

  XLAL_CHECK( tiling != NULL,   XLAL_EINVAL );
  XLAL_CHECK( s == 2 || s == 3, XLAL_EINVAL );
  XLAL_CHECK( fmin <= fmax,     XLAL_EINVAL, "fmin greater than fmax: [%g, %g]", fmin, fmax );
  XLAL_CHECK( nmin <= nmax,     XLAL_EINVAL, "nmin greater than nmax: [%g, %g]", nmin, nmax );
  XLAL_CHECK( kmin <= kmax,     XLAL_EINVAL, "kmin greater than kmax: [%g, %g]", kmin, kmax );

  PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_lower );
  PiecewiseBoundInfo XLAL_INIT_DECL( info_knot_upper );

  info_knot_lower.s = info_knot_upper.s = s;

  info_knot_lower.n = nmax;
  info_knot_lower.k = kmax;

  info_knot_upper.n = nmin;
  info_knot_upper.k = kmin;

  // Setting the first knot bounds
  // - dt is not used
  XLAL_CHECK( XLALSetLatticeTilingConstantBound( tiling, 0, fmin, fmax ) == XLAL_SUCCESS, XLAL_EFAILED );
  XLAL_CHECK( XLALSetLatticeTilingBound( tiling, 1, F1Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper ) == XLAL_SUCCESS, XLAL_EFAILED );

  if ( s == 3 ) {
    XLAL_CHECK( XLALSetLatticeTilingBound( tiling, 2, F2Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper ) == XLAL_SUCCESS, XLAL_EFAILED );
  }

  // Setting the bounds for all following knots
  for ( size_t knot = 1; knot < knots->size; ++knot ) {
    info_knot_lower.dt = info_knot_upper.dt = gsl_vector_get( knots, knot ) - gsl_vector_get( knots, knot - 1 );
    size_t dimindex = s * knot;
    XLAL_CHECK( XLALSetLatticeTilingBound( tiling, dimindex, F0Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper ) == XLAL_SUCCESS, XLAL_EFAILED );
    XLAL_CHECK( XLALSetLatticeTilingBound( tiling, dimindex + 1, F1Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper ) == XLAL_SUCCESS, XLAL_EFAILED );
    if ( s == 3 ) {
      XLAL_CHECK( XLALSetLatticeTilingBound( tiling, dimindex + 2, F2Bound, sizeof( info_knot_lower ), &info_knot_lower, &info_knot_upper ) == XLAL_SUCCESS, XLAL_EFAILED );
    }
  }

  // Disabling using extrema for parameter space bounds
  /*
  for ( size_t dim = 0; dim < s * knots->size; ++dim ) {
    XLAL_CHECK( XLALSetLatticeTilingPadding( tiling, dim, -1, -1, -1, -1, false ) == XLAL_SUCCESS, XLAL_EFAILED );
  }
  */

  for ( size_t dim = 0; dim < s * knots->size; ++dim ) {
    XLAL_CHECK( XLALSetLatticeTilingPadding( tiling, dim, gsl_vector_get( bboxpad, dim ), gsl_vector_get( bboxpad, dim ), gsl_vector_int_get( intpad, dim ), gsl_vector_int_get( intpad, dim ), false ) == XLAL_SUCCESS, XLAL_EFAILED );
  }

  return XLAL_SUCCESS;
}
