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

#ifndef _FLATLATTICETILING_H
#define _FLATLATTICETILING_H

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
/// \brief Lattice-based template generation for flat metric parameter spaces
///

///
/// Function which defines bounds on the lattice tiling parameter space
///
typedef void (*LatticeBoundFunction)(
  const size_t dimension,			///< [in] Dimension on which bound applies
  const gsl_vector* point,			///< [in] Point on which to find bound
  const gsl_vector* bbox,			///< [in] Metric ellipse bounding box
  const void* data,				///< [in] Arbitrary data describing parameter space bound
  double* bound,				///< [out] Bound on point in dimension
  double* padding				///< [out,optional] Padding on parameter space bound (ignored for non-tiled dimensions)
  );

///
/// Type of lattice to generate flat tiling with
///
typedef enum tagLatticeType {
  LATTICE_TYPE_CUBIC,				///< Cubic (\f$Z_n\f$) lattice
  LATTICE_TYPE_ANSTAR,				///< An-star (\f$A_n^*\f$) lattice
  LATTICE_TYPE_MAX
} LatticeType;

///
/// Lattice tiling state structure
///
typedef struct tagLatticeTiling LatticeTiling;

///
/// Find the bounding box of the mismatch ellipses of a metric
///
gsl_vector* XLALMetricEllipseBoundingBox(
  const gsl_matrix* metric,			///< [in] Metric to bound
  const double max_mismatch			///< [in] Maximum mismatch with respect to metric
  );

///
/// Compute a lower triangular basis matrix whose columns are orthonormal with respect to a given metric
///
gsl_matrix* XLALComputeMetricOrthoBasis(
  const gsl_matrix* metric			///< [in] Metric to orthonormalise with respect to
  );

///
/// Compute a lower triangular generator matrix for a given lattice type and mismatch
///
gsl_matrix* XLALComputeLatticeGenerator(
  const size_t dimensions,			///< [in] Number of dimensions
  const LatticeType lattice,			///< [in] Lattice type
  const double max_mismatch			///< [in] Maximum prescribed mismatch
  );

///
/// Create a new lattice tiling state structure
///
LatticeTiling* XLALCreateLatticeTiling(
  const size_t dimensions			///< [in] Number of parameter-space dimensions
  );

///
/// Destroy a lattice tiling state structure
///
void XLALDestroyLatticeTiling(
  LatticeTiling* tiling				///< [in] Tiling state
  );

///
/// Return the total number of dimensions of the lattice tiling
///
size_t XLALGetLatticeTotalDimensions(
  const LatticeTiling* tiling			///< [in] Tiling state
  );

///
/// Return the number of tiled dimensions of the lattice,
/// i.e. exclusing dimensions containing only a single point
///
size_t XLALGetLatticeTiledDimensions(
  const LatticeTiling* tiling			///< [in] Tiling state
  );

///
/// Return the current number of lattice tiling parameter-space points
///
uint64_t XLALGetLatticePointCount(
  const LatticeTiling* tiling			///< [in] Tiling state
  );

///
/// Calculate the total number of lattice tiling parameter-space points
///
uint64_t XLALCountLatticePoints(
  LatticeTiling* tiling				///< [in] Tiling state
  );

///
/// Set a parameter-space bound on a dimension of the lattice tiling.
/// The bound is described by a function \c func, and two data of length
/// \c data_len, \c data_lower and \c data_upper, describing the lower
/// and upper parameter space bounds respectively. If \c data_lower and
/// \c data_upper are identical, this parameter-space dimension will be
/// treated as a single point, and will not be tiled.
///
int XLALSetLatticeBound(
  LatticeTiling* tiling,			///< [in] Tiling state
  const size_t dimension,			///< [in] Dimension on which bound applies
  const LatticeBoundFunction func,		///< [in] Parameter space bound function
  const size_t data_len,			///< [in] Length of arbitrary data describing parameter space bounds
  void* data_lower,				///< [in] Arbitrary data describing lower parameter space bound
  void* data_upper				///< [in] Arbitrary data describing upper parameter space bound
  );

///
/// Set a constant parameter-space bound, given by the minimum and maximum
/// of the two supplied bounds, on a dimension of the lattice tiling
///
int XLALSetLatticeConstantBound(
  LatticeTiling* tiling,			///< [in] Tiling state
  const size_t dimension,			///< [in] Dimension on which bound applies
  const double bound1,				///< [in] First bound on dimension
  const double bound2				///< [in] Second bound on dimension
  );

///
/// Set the flat tiling lattice type, metric and maximum mismatch
///
int XLALSetLatticeTypeAndMetric(
  LatticeTiling* tiling,			///< [in] Tiling state
  const LatticeType lattice,			///< [in] Lattice type
  const gsl_matrix* metric,			///< [in] parameter-space metric
  const double max_mismatch			///< [in] Maximum prescribed mismatch
  );

///
/// Move to the next point in the lattice tiling parameter space. Returns a
/// positive number while the template bank is being generated, or -1 once the
/// template bank is exhausted. Optionally, return the current lattice point.
///
int XLALNextLatticePoint(
  LatticeTiling* tiling,			///< [in] Tiling state
  gsl_vector* curr_point			///< [in/out] Current lattice point
  );

///
/// Fast-forward the lattice tiling through the highest tiled dimension of
/// the parameter space, so that then calling XLALNextLatticePoint() will
/// advance the next highest tiled dimension. Optionally, return the count of
/// and spacing between the points fast-forwarded over.
///
int XLALFastForwardLatticeTiling(
  LatticeTiling* tiling,			///< [in] Tiling state
  uint32_t *point_count,			///< [out] Count of points fast-forwarded over
  double *point_spacing				///< [out] Spacing between points fast-forwarded over
  );

///
/// Return to the beginning of a lattice tiling
///
int XLALRestartLatticeTiling(
  LatticeTiling* tiling				///< [in] Tiling state
  );

///
/// Generate random points within the lattice tiling parameter space
///
int XLALRandomLatticePoints(
  const LatticeTiling* tiling,			///< [in] Tiling state
  RandomParams* rng,				///< [in] Random number generator
  gsl_matrix* random_points			///< [in/out] Matrix whose columns are the random points
  );

///
/// Build a lookup of the index to each point in the lattice tiling; required by
/// XLALNearestLatticePoints() when finding the indices of nearest lattice points
///
int XLALBuildLatticeIndexLookup(
  LatticeTiling* tiling				///< [in] Tiling state
  );

///
/// Print the lookup generated by XLALBuildLatticeIndexLookup() to the given file pointer
///
int XLALPrintLatticeIndexLookup(
  const LatticeTiling* tiling,			///< [in] Tiling state
  FILE* file					///< [in] File pointer to print to
  );

///
/// Find the nearest point in the lattice tiling, and optionally its index, to each given
///
int XLALNearestLatticePoints(
  const LatticeTiling* tiling,			///< [in] Tiling state
  const gsl_matrix* points,			///< [in] Matrix whose columns are the given points
  gsl_matrix* nearest_points,			///< [in/out] Matrix whose columns are the nearest points
  UINT8Vector* nearest_indices			///< [in/out] Vector of indices of nearest points
  );

#ifdef __cplusplus
}
#endif

#endif
