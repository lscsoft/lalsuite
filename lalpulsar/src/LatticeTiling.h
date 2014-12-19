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
/// \file
/// \author Karl Wette
/// \brief Lattice-based template generation for constant-metric parameter spaces
///

///
/// Lattice tiling parameter space
///
typedef struct tagLatticeTilingSpace LatticeTilingSpace;

///
/// Lattice tiling state structure
///
typedef struct tagLatticeTiling LatticeTiling;

///
/// Lattice tiling band iterator
///
typedef struct tagLatticeTilingBandIterator LatticeTilingBandIterator;

///
/// Lattice tiling iterator
///
typedef struct tagLatticeTilingIterator LatticeTilingIterator;

///
/// Type of lattice to generate tiling with
///
typedef enum tagTilingLattice {
  TILING_LATTICE_CUBIC,			///< Cubic (\f$Z_n\f$) lattice
  TILING_LATTICE_ANSTAR,		///< An-star (\f$A_n^*\f$) lattice
  TILING_LATTICE_MAX
} TilingLattice;

///
/// Function which returns a bound on a dimension of the lattice tiling parameter space
///
typedef double (*LatticeTilingBound)(
  const void* data,			///< [in] Arbitrary data describing parameter space bound
  const size_t dim,			///< [in] Dimension on which bound applies
  const gsl_vector* point		///< [in] Point at which to find bound
  );

///
/// Create a new lattice tiling parameter space
///
LatticeTilingSpace* XLALCreateLatticeTilingSpace(
  const size_t n			///< [in] Number of parameter-space dimensions
  );

///
/// Destroy a lattice tiling space
///
void XLALDestroyLatticeTilingSpace(
  LatticeTilingSpace* space		///< [in] Lattice tiling parameter space
  );

///
/// Set a parameter-space bound on a dimension of the lattice tiling.  The bound is described by a function
/// \c func, and two data of length \c data_len, \c data_lower and \c data_upper, describing the lower and
/// upper parameter space bounds respectively. If \c data_lower and \c data_upper are identical, this
/// parameter-space dimension will be treated as a single point, and will not be tiled.
///
int XLALSetLatticeTilingBound(
  LatticeTilingSpace* space,		///< [in] Lattice tiling parameter space
  const size_t dim,			///< [in] Dimension on which bound applies
  const LatticeTilingBound func,	///< [in] Parameter space bound function
  const size_t data_len,		///< [in] Length of arbitrary data describing parameter space bounds
  void* data_lower,			///< [in] Arbitrary data describing lower parameter space bound
  void* data_upper			///< [in] Arbitrary data describing upper parameter space bound
  );

///
/// Set a constant lattice tiling parameter-space bound, given by the minimum and maximum of the two supplied
/// bounds, on a dimension of the lattice tiling.
///
int XLALSetLatticeTilingConstantBound(
  LatticeTilingSpace* space,		///< [in] Lattice tiling parameter space
  const size_t dim,			///< [in] Dimension on which bound applies
  const double bound1,			///< [in] First bound on dimension
  const double bound2			///< [in] Second bound on dimension
  );

///
/// Generate random points within the lattice tiling parameter space
///
int XLALRandomLatticeTilingPoints(
  const LatticeTilingSpace* space,	///< [in] Lattice tiling parameter space
  RandomParams* rng,			///< [in] Random number generator
  gsl_matrix* random_points		///< [out] Matrix whose columns are the random points
  );

///
/// Create a lattice tiling given a parameter space, lattice type, metric, and maximum mismatch
///
LatticeTiling* XLALCreateLatticeTiling(
  const LatticeTilingSpace* space,	///< [in] Lattice tiling parameter space
  const TilingLattice lattice,		///< [in] Type of lattice to generate tiling with
  const gsl_matrix* metric,		///< [in] Parameter-space metric
  const double max_mismatch		///< [in] Maximum mismatch
  );

///
/// Destroy a lattice tiling state structure
///
void XLALDestroyLatticeTiling(
  LatticeTiling* tiling			///< [in] Lattice tiling state structure
  );

///
/// Return the total number of dimensions of the lattice tiling
///
size_t XLALTotalLatticeTilingDimensions(
  const LatticeTiling* tiling		///< [in] Lattice tiling state structure
  );

///
/// Return the number of tiled dimensions of the lattice tiling
///
size_t XLALTiledLatticeTilingDimensions(
  const LatticeTiling* tiling		///< [in] Lattice tiling state structure
  );

///
/// Return the total number of points in the lattice tiling
///
UINT8 XLALLatticeTilingTotalPointCount(
  const LatticeTiling* tiling		///< [in] Lattice tiling state structure
  );

///
/// Return the number of unique points in the lattice tiling up to the specified dimension
///
UINT8 XLALLatticeTilingUniquePointCount(
  const LatticeTiling* tiling,		///< [in] Lattice tiling state structure
  const size_t dim			///< [in] Dimension up to which to count points
  );

///
/// Return the unique points in the lattice tiling up to the specified dimension
///
gsl_matrix* XLALLatticeTilingUniquePoints(
  const LatticeTiling* tiling,		///< [in] Lattice tiling state structure
  const size_t dim			///< [in] Dimension up to which to count points
  );

///
/// Create an iterator over all unique points in the lattice tiling up to the specified dimension
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(RETURNS_PROPERTY(LatticeTilingIterator*, XLALCreateLatticeTilingIterator));
#endif
LatticeTilingIterator* XLALCreateLatticeTilingIterator(
  const LatticeTiling* tiling,		///< [in] Lattice tiling state structure
  const size_t dim			///< [in] Dimension up to which to iterate over points
  );

///
/// Destroy an lattice tiling iterator
///
void XLALDestroyLatticeTilingIterator(
  LatticeTilingIterator* itr		///< [in] Lattice tiling iterator
  );

///
/// Advance iterator and optionally return the current point in \c point. Returns >0 if there are
/// points remaining, 0 if there are no more points, and XLAL_FAILURE on error.
///
int XLALNextLatticeTilingPoint(
  LatticeTilingIterator* itr,		///< [in] Lattice tiling iterator
  gsl_vector* point			///< [out] Current point up to iterator dimension
  );

///
/// Given a set of points, find their nearest points in the lattice tiling, and optionally their indices
///
#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(gsl_matrix**, nearest_points));
SWIGLAL(INOUT_STRUCTS(UINT8Vector**, nearest_indices));
#endif
int XLALNearestLatticeTilingPoints(
  const LatticeTiling* tiling,		///< [in] Lattice tiling state structure
  const gsl_matrix* points,		///< [in] Columns are set of points for which to find nearest points
  gsl_matrix** nearest_points,		///< [out] Columns are the corresponding nearest points
  UINT8Vector** nearest_indices		///< [out] Corresponding tiling indices of nearest points
  );

#ifdef __cplusplus
}
#endif

#endif
