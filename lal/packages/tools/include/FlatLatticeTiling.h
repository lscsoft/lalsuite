//
// Copyright (C) 2007, 2008, 2012 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

///
/// \addtogroup FlatLatticeTiling_h
/// \author Karl Wette
/// \brief Lattice-based template generation for flat metric parameter spaces
///

/// @{

#ifndef _FLATLATTICETILING_H
#define _FLATLATTICETILING_H

#include <stdbool.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALDatatypes.h>
#include <lal/Random.h>
#include <lal/GSLSupport.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Flat lattice tiling bound function
///
typedef void (*FlatLatticeBound)(
  const size_t dimension,	///< [in] Dimension on which bound applies
  const gsl_vector_uint* bound,	///< [in] Indices of current bounds
  const gsl_vector* point,	///< [in] Point on which to find bounds
  const void* data,		///< [in] Arbitrary data describing parameter space
  gsl_vector* lower,		///< [out] Lower bounds on point in dimension
  gsl_vector* upper		///< [out] Upper bounds on point in dimension
  );

///
/// Flat tiling lattice generator function
///
typedef int (*FlatLatticeGenerator)(
  const size_t dimensions,	///< [in] Number of dimensions
  gsl_matrix** generator,	///< [out] Generator matrix
  double* norm_thickness	///< [out] Normalised thickness
  );

///
/// Flat lattice tiling state structure
///
typedef struct tagFlatLatticeTiling FlatLatticeTiling;

///
/// Create a new flat lattice tiling state structure
///
FlatLatticeTiling* XLALCreateFlatLatticeTiling(
  const size_t dimensions	///< [in] Number of parameter space dimensions
  );

///
/// Destroy a flat lattice tiling state structure
///
void XLALDestroyFlatLatticeTiling(
  FlatLatticeTiling* tiling	///< [in] Tiling state
  );

///
/// Return the number of dimensions being tiled
///
size_t XLALGetFlatLatticeDimensions(
  FlatLatticeTiling* tiling	///< [in] Tiling state
  );

///
/// Return the current number of flat lattice tiling parameter space points
///
unsigned long XLALGetFlatLatticePointCount(
  FlatLatticeTiling* tiling	///< [in] Tiling state
  );

///
/// Add a parameter space bound to the flat lattice tiling
///
int XLALSetFlatLatticeBound(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const size_t dimension,		///< [in] Dimension on which bound applies
  const bool singular,			///< [in] Is bound composed of single points?
  const FlatLatticeBound func,		///< [in] Parameter space bound function
  void* data				///< [in] Arbitrary data describing parameter space
  );

///
/// Set the flat tiling lattice generator
///
int XLALSetFlatLatticeGenerator(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const FlatLatticeGenerator generator	///< [in] Lattice generator function
  );

///
/// Set the flat lattice tiling metric and maximum mismatch
///
int XLALSetFlatLatticeMetric(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const gsl_matrix* metric,		///< [in] Parameter space metric
  const double max_mismatch		///< [in] Maximum prescribed mismatch
  );

#ifdef SWIG // SWIG interface directives
SWIGLAL(NO_NEW_OBJECT(XLALNextFlatLatticePoint));
#endif
///
/// Return the next point in the flat lattice tiling parameter space
///
gsl_vector* XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling		///< [in] Tiling state
  );

///
/// Return to the beginning of a flat lattice tiling
///
int XLALRestartFlatLatticeTiling(
  FlatLatticeTiling* tiling		///< [in] Tiling state
  );

///
/// Calculate the total number of flat lattice tiling parameter space points
///
unsigned long XLALCountTotalFlatLatticePoints(
  FlatLatticeTiling* tiling		///< [in] Tiling state
  );

#ifdef SWIG // SWIG interface directives
SWIGLAL(FUNCTION_POINTER(XLALCubicLatticeGenerator));
#endif
///
/// Calculate the generator matrix for a cubic (\f$Z_n\f$) lattice
///
int XLALCubicLatticeGenerator(
  const size_t dimensions,	///< [in] Number of dimensions
  gsl_matrix** generator,	///< [out] Generator matrix
  double* norm_thickness	///< [out] Normalised thickness
  );

#ifdef SWIG // SWIG interface directives
SWIGLAL(FUNCTION_POINTER(XLALAnstarLatticeGenerator));
#endif
///
/// Calculate the generator matrix for a \f$A_n^*\f$ lattice
///
int XLALAnstarLatticeGenerator(
  const size_t dimensions,	///< [in] Number of dimensions
  gsl_matrix** generator,	///< [out] Generator matrix
  double* norm_thickness	///< [out] Normalised thickness
  );

///
/// Add a constant parameter space bound to the flat lattice tiling
///
int XLALSetFlatLatticeConstantBound(
  FlatLatticeTiling* tiling,	///< [in] Tiling state
  const size_t dimension,	///< [in] Dimension on which bound applies
  const double lower,		///< [in] Lower bound on dimension
  const double upper		///< [in] Upper bound on dimension
  );

///
/// Find the bounding box of the mismatch ellipses of a metric
///
gsl_vector* XLALMetricEllipseBoundingBox(
  gsl_matrix* metric,		///< [in] Metric to bound
  const double max_mismatch	///< [in] Maximum mismatch with respect to metric
  );

///
/// Orthonormalise the columns of a matrix with respect to a metric (matrix is lower triangular)
///
int XLALOrthonormaliseWRTMetric(
  gsl_matrix* matrix,		///< [in] Matrix of columns to orthonormalise
  const gsl_matrix* metric	///< [in] Metric to orthonormalise with respect to
  );

///
/// Transform a lattice generator to a square lower triangular form
///
gsl_matrix* XLALSquareLowerTriangularLatticeGenerator(
  gsl_matrix* generator		///< [in] Generator matrix of lattice
  );

///
/// Normalise a lattice generator matrix to have a specified covering radius
///
int XLALNormaliseLatticeGenerator(
  gsl_matrix* generator,	///< [in] Generator matrix of lattice
  const double norm_thickness,	///< [in] Normalised thickness of lattice
  const double covering_radius	///< [in] Desired covering radius
  );

#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(gsl_matrix**, random_points, nearest_points, workspace));
SWIGLAL(INOUT_STRUCTS(gsl_vector**, nearest_distances));
SWIGLAL(INOUT_STRUCTS(gsl_vector_ulong**, nearest_indices));
#endif
///
/// Generate random points within the flat lattice tiling parameter space,
/// then calculate the nearest flat lattice point to each random point
///
int XLALNearestFlatLatticePointToRandomPoints(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  RandomParams* rng,			///< [in] Random number generator
  const size_t num_random_points,	///< [in] Number of random points to generate
  gsl_matrix** random_points,		///< [in/out] Pointer to matrix of random points
  gsl_matrix** nearest_points,		///< [in/out] Pointer to matrix of nearest lattice points to each random point
  gsl_vector_ulong** nearest_indices,	///< [in/out] Pointer to vector of indices of nearest lattice point
  gsl_vector** nearest_distances,	///< [in/out] Pointer to vector of distances to nearest lattice point
  gsl_matrix** workspace		///< [in/out] Pointer to workspace matrix for computing distances
  );

#ifdef __cplusplus
}
#endif

#endif

/// @}
