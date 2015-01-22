//
// Copyright (C) 2007, 2008, 2012, 2014 Karl Wette
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
/// \brief Lattice-based template generation for constant-metric parameter
/// spaces, described in \cite Wette2009a and \cite Wette2014a .
///
/// @{
///

///
/// Lattice tiling parameter space.
///
typedef struct tagLatticeTilingSpace LatticeTilingSpace;

///
/// Lattice tiling state structure.
///
typedef struct tagLatticeTiling LatticeTiling;

///
/// Lattice tiling iterator.
///
typedef struct tagLatticeTilingIterator LatticeTilingIterator;

///
/// Type of lattice to generate tiling with.
///
typedef enum tagTilingLattice {
  /// Cubic (\f$Z_n\f$) lattice.
  TILING_LATTICE_CUBIC,
  /// An-star (\f$A_n^*\f$) lattice.
  TILING_LATTICE_ANSTAR,
  /// \cond DONT_DOXYGEN
  TILING_LATTICE_MAX
  /// \endcond
} TilingLattice;

///
/// Function which returns a bound on a dimension of the lattice tiling
/// parameter space.
///
typedef double (*LatticeTilingBound)(
  /// [in] Arbitrary data describing parameter space bound.
  const void* data,
  /// [in] Dimension on which bound applies.
  const size_t dim,
  /// [in] Point at which to find bound.
  const gsl_vector* point
  );

///
/// Compute the extent of the bounding box of the mismatch ellipse of a metric.
///
gsl_vector* XLALMetricEllipseBoundingBox(
  /// [in] Parameter-space metric.
  const gsl_matrix* metric,
  /// [in] Maximum mismatch.
  const double max_mismatch
  );

///
/// Create a new lattice tiling parameter space.
///
LatticeTilingSpace* XLALCreateLatticeTilingSpace(
  /// [in] Number of parameter-space dimensions.
  const size_t ndim
  );

///
/// Destroy a lattice tiling space.
///
void XLALDestroyLatticeTilingSpace(
  /// [in] Lattice tiling parameter space.
  LatticeTilingSpace* space
  );

///
/// Set a parameter-space bound on a dimension of the lattice tiling.
/// The bound is described by a function \c func, and two data of length
/// \c data_len, \c data_lower and \c data_upper, describing the lower
/// and upper parameter space bounds respectively. If \c data_lower and
/// \c data_upper are identical, this parameter-space dimension will be
/// treated as a single point, and will not be tiled.
///
int XLALSetLatticeTilingBound(
  /// [in] Lattice tiling parameter space.
  LatticeTilingSpace* space,
  /// [in] Dimension on which bound applies.
  const size_t dim,
  /// [in] Parameter space bound function.
  const LatticeTilingBound func,
  /// [in] Length of arbitrary data describing parameter space bounds.
  const size_t data_len,
  /// [in] Arbitrary data describing lower parameter space bound.
  void* data_lower,
  /// [in] Arbitrary data describing upper parameter space bound.
  void* data_upper
  );

///
/// Set a constant lattice tiling parameter-space bound, given by the
/// minimum and maximum of the two supplied bounds, on a dimension of
/// the lattice tiling.
///
int XLALSetLatticeTilingConstantBound(
  /// [in] Lattice tiling parameter space.
  LatticeTilingSpace* space,
  /// [in] Dimension on which bound applies.
  const size_t dim,
  /// [in] First bound on dimension.
  const double bound1,
  /// [in] Second bound on dimension.
  const double bound2
  );

///
/// Generate random points within the lattice tiling parameter space.
///
int XLALRandomLatticeTilingPoints(
  /// [in] Lattice tiling parameter space.
  const LatticeTilingSpace* space,
  /// [in] Random number generator.
  RandomParams* rng,
  /// [out] Matrix whose columns are the random points.
  gsl_matrix* random_points
  );

///
/// Create a lattice tiling given a parameter space, lattice type,
/// metric, and maximum mismatch.
///
LatticeTiling* XLALCreateLatticeTiling(
  /// [in] Lattice tiling parameter space.
  const LatticeTilingSpace* space,
  /// [in] Type of lattice to generate tiling with.
  const TilingLattice lattice,
  /// [in] Parameter-space metric.
  const gsl_matrix* metric,
  /// [in] Maximum mismatch.
  const double max_mismatch
  );

///
/// Destroy a lattice tiling state structure.
///
void XLALDestroyLatticeTiling(
  /// [in] Lattice tiling state structure.
  LatticeTiling* tiling
  );

///
/// Return the total number of dimensions of the lattice tiling.
///
size_t XLALTotalLatticeTilingDimensions(
  /// [in] Lattice tiling state structure.
  const LatticeTiling* tiling
  );

///
/// Return the number of tiled dimensions of the lattice tiling.
///
size_t XLALTiledLatticeTilingDimensions(
  /// [in] Lattice tiling state structure.
  const LatticeTiling* tiling
  );

///
/// Return the total number of points in the lattice tiling.
///
UINT8 XLALLatticeTilingTotalPointCount(
  /// [in] Lattice tiling state structure.
  const LatticeTiling* tiling
  );

///
/// Return the number of unique points in the lattice tiling up to the
/// specified dimension.
///
UINT8 XLALLatticeTilingUniquePointCount(
  /// [in] Lattice tiling state structure.
  const LatticeTiling* tiling,
  /// [in] Dimension up to which to count points.
  const size_t dim
  );

///
/// Return the unique points in the lattice tiling up to the specified dimension.
///
gsl_matrix* XLALLatticeTilingUniquePoints(
  /// [in] Lattice tiling state structure.
  const LatticeTiling* tiling,
  /// [in] Dimension up to which to count points.
  const size_t dim
  );

///
/// Create an iterator over all unique points in the lattice tiling up
/// to the specified dimension.
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(RETURNS_PROPERTY(LatticeTilingIterator*, XLALCreateLatticeTilingIterator));
#endif
LatticeTilingIterator* XLALCreateLatticeTilingIterator(
  /// [in] Lattice tiling state structure.
  const LatticeTiling* tiling,
  /// [in] Dimension up to which to iterate over points.
  const size_t dim
  );

///
/// Destroy an lattice tiling iterator.
///
void XLALDestroyLatticeTilingIterator(
  /// [in] Lattice tiling iterator.
  LatticeTilingIterator* itr
  );

///
/// Advance iterator and optionally return the current point in \c
/// point. Returns >0 if there are points remaining, 0 if there are no
/// more points, and XLAL_FAILURE on error.
///
int XLALNextLatticeTilingPoint(
  /// [in] Lattice tiling iterator.
  LatticeTilingIterator* itr,
  /// [out] Current point up to iterator dimension.
  gsl_vector* point
  );

///
/// Given a set of points, find their nearest points in the lattice
/// tiling, and optionally their indices.
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(gsl_matrix**, nearest_points));
SWIGLAL(INOUT_STRUCTS(UINT8Vector**, nearest_indices));
#endif
int XLALNearestLatticeTilingPoints(
  /// [in] Lattice tiling state structure.
  const LatticeTiling* tiling,
  /// [in] Columns are set of points for which to find nearest points.
  const gsl_matrix* points,
  /// [out] Columns are the corresponding nearest points.
  gsl_matrix** nearest_points,
  /// [out] Corresponding tiling indices of nearest points.
  UINT8Vector** nearest_indices
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _LATTICETILING_H
