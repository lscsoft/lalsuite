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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
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
/// Flat lattice tiling bound function
///
typedef void (*FlatLatticeBound)(
  const size_t dimension,                       ///< [in] Dimension on which bound applies
  const gsl_vector_uint* bound,                 ///< [in] Indices of current bounds
  const gsl_vector* point,                      ///< [in] Point on which to find bounds
  const void* data,                             ///< [in] Arbitrary data describing parameter space
  const gsl_vector* incr,                       ///< [in] Increments of the lattice tiling generator
  const gsl_vector* bbox,                       ///< [in] Metric ellipse bounding box extents
  gsl_vector* lower,                            ///< [out] Lower bounds on point in dimension
  gsl_vector* upper,                            ///< [out] Upper bounds on point in dimension
  double* lower_pad,                            ///< [out] Padding of lower parameter space bounds
  double* upper_pad                             ///< [out] Padding of upper parameter space bounds
  );

///
/// Type of lattice to generate flat tiling with
///
typedef enum tagFlatLatticeType {
  FLAT_LATTICE_TYPE_CUBIC,                      ///< Cubic (\f$Z_n\f$) lattice
  FLAT_LATTICE_TYPE_ANSTAR,                     ///< An-star (\f$A_n^*\f$) lattice
  FLAT_LATTICE_TYPE_MAX
} FlatLatticeType;

///
/// Flat lattice tiling state structure
///
typedef struct tagFlatLatticeTiling FlatLatticeTiling;

///
/// Create a new flat lattice tiling state structure
///
FlatLatticeTiling* XLALCreateFlatLatticeTiling(
  const size_t dimensions                       ///< [in] Number of parameter space dimensions
  );

///
/// Destroy a flat lattice tiling state structure
///
void XLALDestroyFlatLatticeTiling(
  FlatLatticeTiling* tiling                     ///< [in] Tiling state
  );

///
/// Return the number of dimensions being tiled
///
size_t XLALGetFlatLatticeDimensions(
  FlatLatticeTiling* tiling                     ///< [in] Tiling state
  );

#ifdef SWIG // SWIG interface directives
SWIGLAL(RETURNS_PROPERTY(XLALGetFlatLatticePoint));
#endif
///
/// Return the current lattice tiling parameter space point
///
const gsl_vector* XLALGetFlatLatticePoint(
  FlatLatticeTiling* tiling                     ///< [in] Tiling state
  );

///
/// Return the current number of flat lattice tiling parameter space points
///
unsigned long XLALGetFlatLatticePointCount(
  FlatLatticeTiling* tiling                     ///< [in] Tiling state
  );

///
/// Return the increment vectors which are used to generate the lattice.
///
gsl_matrix* XLALGetFlatLatticeIncrements(
  FlatLatticeTiling* tiling                     ///< [in] Tiling state
  );

///
/// Set a parameter space bound on the flat lattice tiling
///
int XLALSetFlatLatticeBound(
  FlatLatticeTiling* tiling,                    ///< [in] Tiling state
  const size_t dimension,                       ///< [in] Dimension on which bound applies
  const bool singular,                          ///< [in] Is bound composed of single points?
  const FlatLatticeBound func,                  ///< [in] Parameter space bound function
  void* data                                    ///< [in] Arbitrary data describing parameter space
  );

///
/// Set the flat tiling lattice type, metric and maximum mismatch
///
int XLALSetFlatLatticeTypeAndMetric(
  FlatLatticeTiling* tiling,                    ///< [in] Tiling state
  const FlatLatticeType lattice,                ///< [in] Lattice type
  const gsl_matrix* metric,                     ///< [in] Parameter space metric
  const double max_mismatch                     ///< [in] Maximum prescribed mismatch
  );

///
/// Move to the next point in the flat lattice tiling parameter space.
/// Returns the index of the lowest dimension where the point has changed,
/// or a negative number when the template bank is exhausted.
///
int XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling                     ///< [in] Tiling state
  );

///
/// Return to the beginning of a flat lattice tiling
///
int XLALRestartFlatLatticeTiling(
  FlatLatticeTiling* tiling                     ///< [in] Tiling state
  );

///
/// Calculate the total number of flat lattice tiling parameter space points
///
unsigned long XLALCountTotalFlatLatticePoints(
  FlatLatticeTiling* tiling                     ///< [in] Tiling state
  );

#ifdef SWIG // SWIG interface directives
SWIGLAL(INOUT_STRUCTS(gsl_matrix**, random_points, workspace));
SWIGLAL(INOUT_STRUCTS(gsl_vector**, nearest_distances));
SWIGLAL(INOUT_STRUCTS(gsl_vector_ulong**, nearest_indices));
#endif
///
/// Generate random points within the flat lattice tiling parameter space,
/// then calculate the nearest flat lattice point to each random point
///
int XLALNearestFlatLatticePointToRandomPoints(
  FlatLatticeTiling* tiling,                    ///< [in] Tiling state
  RandomParams* rng,                            ///< [in] Random number generator
  const size_t num_random_points,               ///< [in] Number of random points to generate
  gsl_matrix** random_points,                   ///< [in/out] Pointer to matrix of random points
  gsl_vector_ulong** nearest_indices,           ///< [in/out] Pointer to vector of indices of nearest lattice point
  gsl_vector** nearest_distances,               ///< [in/out] Pointer to vector of distances to nearest lattice point
  gsl_matrix** workspace                        ///< [in/out] Pointer to workspace matrix for computing distances
  );

///
/// Set a constant parameter space bound, given by the minimum and
/// maximum of the two supplied bounds, on the flat lattice tiling
///
int XLALSetFlatLatticeConstantBound(
  FlatLatticeTiling* tiling,                    ///< [in] Tiling state
  const size_t dimension,                       ///< [in] Dimension on which bound applies
  const double bound1,                          ///< [in] First bound on dimension
  const double bound2                           ///< [in] Second bound on dimension
  );

///
/// Set elliptical bounds in two dimensions on the flat lattice tiling
///
int XLALSetFlatLatticeEllipticalBounds(
  FlatLatticeTiling* tiling,                    ///< [in] Tiling state
  const size_t dimension,                       ///< [in] Dimension of X bound (Y bound is one higher)
  const double x_centre,                        ///< [in] X centre of ellipse
  const double y_centre,                        ///< [in] Y centre of ellipse
  const double x_semi,                          ///< [in] Length of X semi-diameter
  const double y_semi                           ///< [in] Length of Y semi-diameter
  );

#ifdef __cplusplus
}
#endif

#endif
