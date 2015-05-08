//
// Copyright (C) 2007, 2008, 2012, 2014, 2015 Karl Wette
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

#ifndef _LATTICETILING_H
#define _LATTICETILING_H

#include <stdbool.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <lal/LALStdlib.h>
#include <lal/Random.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// \defgroup LatticeTiling_h Header LatticeTiling.h
/// \ingroup lalpulsar_templbank
/// \author Karl Wette
/// \brief Lattice-based template generation for constant-metric parameter spaces, described in
/// \cite Wette2009a and \cite Wette2014a .
///
/// @{
///

///
/// Describes a lattice tiling parameter-space bounds and metric.
///
typedef struct tagLatticeTiling LatticeTiling;

///
/// Iterates over all points in a lattice tiling.
///
typedef struct tagLatticeTilingIterator LatticeTilingIterator;

///
/// Locates the nearest point in a lattice tiling.
///
typedef struct tagLatticeTilingLocator LatticeTilingLocator;

///
/// Type of lattice to generate tiling with.
///
typedef enum tagTilingLattice {
  TILING_LATTICE_CUBIC,			///< Cubic (\f$Z_n\f$) lattice
  TILING_LATTICE_ANSTAR,		///< An-star (\f$A_n^*\f$) lattice
  TILING_LATTICE_MAX
} TilingLattice;

///
/// Order in which to iterate over lattice tiling points with a lattice tiling iterator.
///
typedef enum tagTilingOrder {
  TILING_ORDER_POSITIVE,		///< Iterate in positive order (i.e. lower bound to upper bound)
  TILING_ORDER_ALTERNATING,		///< Alternate between positive and negative order (i.e. upper bound to lower bound) after every pass over each dimension
  TILING_ORDER_MAX
} TilingOrder;

///
/// Function which returns a bound on a dimension of the lattice tiling.
///
typedef double( *LatticeTilingBound )(
  const void *data,			///< [in] Arbitrary data describing parameter space bound
  const size_t dim,			///< [in] Dimension on which bound applies
  const gsl_vector *point		///< [in] Point at which to find bound
  );

///
/// Create a new lattice tiling.
///
LatticeTiling *XLALCreateLatticeTiling(
  const size_t ndim			///< [in] Number of parameter-space dimensions
  );

///
/// Destroy a lattice tiling.
///
void XLALDestroyLatticeTiling(
  LatticeTiling *tiling			///< [in] Lattice tiling
  );

///
/// Set a parameter-space bound on a dimension of the lattice tiling.  The bound is described by a
/// function \c func, and two data of length \c data_len, \c data_lower and \c data_upper,
/// describing the lower and upper parameter space bounds respectively. If \c data_lower and \c
/// data_upper are identical, this parameter-space dimension will be treated as a single point, and
/// will not be tiled.
///
int XLALSetLatticeTilingBound(
  LatticeTiling *tiling,		///< [in] Lattice tiling
  const size_t dim,			///< [in] Dimension on which bound applies
  const LatticeTilingBound func,	///< [in] Parameter space bound function
  const size_t data_len,		///< [in] Length of arbitrary data describing parameter space bounds
  void *data_lower,			///< [in] Arbitrary data describing lower parameter space bound
  void *data_upper			///< [in] Arbitrary data describing upper parameter space bound
  );

///
/// Set a constant lattice tiling parameter-space bound, given by the minimum and maximum of the two
/// supplied bounds, on a dimension of the lattice tiling.
///
int XLALSetLatticeTilingConstantBound(
  LatticeTiling *tiling,		///< [in] Lattice tiling
  const size_t dim,			///< [in] Dimension on which bound applies
  const double bound1,			///< [in] First bound on dimension
  const double bound2			///< [in] Second bound on dimension
  );

///
/// Set the tiling lattice, parameter-space metric, and maximum prescribed mismatch.  The lattice
/// tiling \c tiling is now fully initialised, and can be used to create tiling iterators [via
/// XLALCreateLatticeTilingIterator()] and locators [via XLALCreateLatticeTilingLocator()].
///
int XLALSetTilingLatticeAndMetric(
  LatticeTiling *tiling,		///< [in] Lattice tiling
  const TilingLattice lattice,		///< [in] Type of lattice to generate tiling with
  const gsl_matrix *metric,		///< [in] Parameter-space metric
  const double max_mismatch		///< [in] Maximum prescribed mismatch
  );

///
/// Return the total number of dimensions of the lattice tiling.
///
size_t XLALNumberOfLatticeTilingDimensions(
  const LatticeTiling *tiling		///< [in] Lattice tiling
  );

///
/// Return the number of tiled dimensions of the lattice tiling.
///
size_t XLALNumberOfTiledLatticeDimensions(
  const LatticeTiling *tiling		///< [in] Lattice tiling
  );

///
/// Return the step size of the lattice tiling in a given dimension, or 0 for non-tiled dimensions.
///
REAL8 XLALLatticeTilingStepSizes(
  const LatticeTiling *tiling,		///< [in] Lattice tiling
  const size_t dim			///< [in] Dimension of which to return step size
  );

///
/// Generate random points within the parameter space of the lattice tiling.  Points can be scaled
/// to fill the parameter space exactly (<tt>scale == 0</tt>), fill a subset of the parameter space
/// (<tt>-1 < scale < 0</tt>), or fill outside the parameter space (<tt>scale > 0</tt>).
///
int XLALRandomLatticeTilingPoints(
  const LatticeTiling *tiling,		///< [in] Lattice tiling
  const double scale,			///< [in] Scale of random points
  RandomParams *rng,			///< [in] Random number generator
  gsl_matrix *random_points		///< [out] Matrix whose columns are the random points
  );

///
/// Create a new lattice tiling iterator.
///
LatticeTilingIterator *XLALCreateLatticeTilingIterator(
  const LatticeTiling *tiling,		///< [in] Lattice tiling
  const size_t itr_ndim,		///< [in] Number of parameter-space dimensions to iterate over
  const TilingOrder order		///< [in] Order in which to iterate over lattice tiling points
  );

///
/// Destroy a lattice tiling iterator.
///
void XLALDestroyLatticeTilingIterator(
  LatticeTilingIterator *itr		///< [in] Lattice tiling iterator
  );

///
/// Reset an iterator to the beginning of a lattice tiling.
///
int XLALResetLatticeTilingIterator(
  LatticeTilingIterator *itr		///< [in] Lattice tiling iterator
  );

///
/// Advance lattice tiling iterator, and optionally return the current point in \c point. Returns >0
/// if there are points remaining, 0 if there are no more points, and XLAL_FAILURE on error.
///
int XLALNextLatticeTilingPoint(
  LatticeTilingIterator *itr,		///< [in] Lattice tiling iterator
  gsl_vector *point			///< [out] Current point
  );

///
/// Return the total number of points in the lattice tiling.
///
UINT8 XLALNumberOfLatticeTilingPoints(
  LatticeTilingIterator *itr		///< [in] Lattice tiling iterator
  );

///
/// Return ranges of the number of lattice tiling points per dimension.
///
int XLALRangesOfLatticeTilingPoints(
  LatticeTilingIterator *itr,		///< [in] Lattice tiling iterator
  const size_t dim,			///< [in] Dimension in which to return ranges
  long *min_points,			///< [out] Minimum number of points in this dimension
  double *avg_points,			///< [out] Average number of points in this dimension
  long *max_points			///< [out] Maximum number of points in this dimension
  );

///
/// Create a new lattice tiling locator. If <tt>bound_ndim > 0</tt>, an index trie is internally built.
///
LatticeTilingLocator *XLALCreateLatticeTilingLocator(
  const LatticeTiling *tiling,		///< [in] Lattice tiling
  const size_t bound_ndim		///< [in] Number of parameter-space dimensions to enforce bounds over
  );

///
/// Destroy a lattice tiling locator.
///
void XLALDestroyLatticeTilingLocator(
  LatticeTilingLocator *loc		///< [in] Lattice tiling locator
  );

///
/// Locate the nearest point in a lattice tiling to a given point. Return the nearest point in
/// \c nearest_point, and optionally: its generating integers in \c nearest_int_point, and its
/// unique index <i>over the bound-enforced dimensions</i> in \c nearest_index (requires \c loc
/// to have been created with an index trie).
///
#ifdef SWIG // SWIG interface directives
SWIGLAL( INOUT_STRUCTS( gsl_vector **, nearest_point ) );
SWIGLAL( INOUT_STRUCTS( gsl_vector_long **, nearest_int_point ) );
SWIGLAL( INOUT_SCALARS( UINT8 *, nearest_index ) );
#endif
int XLALNearestLatticeTilingPoint(
  const LatticeTilingLocator *loc,	///< [in] Lattice tiling locator
  const gsl_vector *point,		///< [in] Point for which to find nearest point
  gsl_vector **nearest_point,		///< [out] Corresponding nearest point
  gsl_vector_long **nearest_int_point,	///< [out] Generating integers of nearest point
  UINT8 *nearest_index			///< [out] Unique tiling index of nearest point
  );

///
/// Multi-point version of XLALNearestLatticeTilingPoint(): locate the nearest points in a lattice
/// tiling to a given a set of points.
///
#ifdef SWIG // SWIG interface directives
SWIGLAL( INOUT_STRUCTS( gsl_matrix **, nearest_points ) );
SWIGLAL( INOUT_STRUCTS( gsl_matrix_long **, nearest_int_points ) );
SWIGLAL( INOUT_STRUCTS( UINT8Vector **, nearest_indexes ) );
#endif
int XLALNearestLatticeTilingPoints(
  const LatticeTilingLocator *loc,	///< [in] Lattice tiling locator
  const gsl_matrix *points,		///< [in] Columns are set of points for which to find nearest points
  gsl_matrix **nearest_points,		///< [out] Columns are the corresponding nearest points
  gsl_matrix_long **nearest_int_points,	///< [out] Columns are the generating integers of the nearest points
  UINT8Vector **nearest_indexes		///< [out] Unique tiling indexes of nearest points
  );

///
/// Print the internal index trie of a lattice tiling locator to the given file pointer.
///
int XLALPrintLatticeTilingIndexTrie(
  const LatticeTilingLocator *loc,	///< [in] Lattice tiling locator
  FILE *file				///< [in] File pointer to print trie to
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _LATTICETILING_H
