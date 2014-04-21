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
  const size_t dimensions                       ///< [in] Number of parameter-space dimensions
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
  const FlatLatticeTiling* tiling               ///< [in] Tiling state
  );

///
/// Return the current number of flat lattice tiling parameter-space points
///
UINT8 XLALGetFlatLatticePointCount(
  const FlatLatticeTiling* tiling               ///< [in] Tiling state
  );

///
/// Calculate the total number of flat lattice tiling parameter-space points
///
UINT8 XLALCountFlatLatticePoints(
  FlatLatticeTiling* tiling                     ///< [in] Tiling state
  );

///
/// Set the parameter-space bounds on dimension \f$n\f$ of the flat lattice tiling.
///
/// The lower/upper bound \f$X_n\f$ is specified by
/// \f[
/// X_n - a_n = \sum_{k=0}^{P-1} c_{N,k} \left[ \sum_{j=0}^{N-1} c_{j,k} \prod_{i=0}^{n-1} (x_i - a_i)^{m_{nj+i,k}} \right]^{m_{nN,k}}
/// \f]
/// where \f$x = (x_0,\cdots,x_{n-1})\f$ is the current parameter-space point, \f$a\f$ is a vector of \f$n + 1 > 0\f$ offsets,
/// \f$c\f$ is a matrix of \f$N + 1 \times P\f$ coefficients, and \f$m\f$ is a matrix of \f$nM + 1 \times P\f$ exponents.
/// Note that \f$x^m \rightarrow 0\f$ if the result is not finite and real-valued.
///
int XLALSetFlatLatticeBound(
  FlatLatticeTiling* tiling,                    ///< [in] Tiling state
  const size_t dimension,                       ///< [in] Dimension on which bound applies (\f$n\f$)
  const gsl_vector* a,                          ///< [in] Vector of offsets (\f$a\f$)
  const gsl_matrix* c_lower,                    ///< [in] Matrix of coefficients (\f$c\f$) for the lower bound
  const gsl_matrix* m_lower,                    ///< [in] Matrix of exponents (\f$m\f$) for the lower bound
  const gsl_matrix* c_upper,                    ///< [in] Matrix of coefficients (\f$c\f$) for the upper bound
  const gsl_matrix* m_upper                     ///< [in] Matrix of exponents (\f$m\f$) for the upper bound
  );

///
/// Set a constant parameter-space bound, given by the minimum and
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

///
/// Set the flat tiling lattice type, metric and maximum mismatch
///
int XLALSetFlatLatticeTypeAndMetric(
  FlatLatticeTiling* tiling,                    ///< [in] Tiling state
  const FlatLatticeType lattice,                ///< [in] Lattice type
  const gsl_matrix* metric,                     ///< [in] parameter-space metric
  const double max_mismatch                     ///< [in] Maximum prescribed mismatch
  );

///
/// Move to the next point in the flat lattice tiling parameter space.
/// Returns the index of the lowest dimension where the point has changed,
/// or a negative number when the template bank is exhausted.
/// Optionally, return the current flat lattice point.
///
int XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling,                    ///< [in] Tiling state
  gsl_vector* curr_point                        ///< [in/out] Current flat lattice point
  );

///
/// Return to the beginning of a flat lattice tiling
///
int XLALRestartFlatLatticeTiling(
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
/* int XLALNearestFlatLatticePointToRandomPoints( */
/*   FlatLatticeTiling* tiling,                    ///< [in] Tiling state */
/*   RandomParams* rng,                            ///< [in] Random number generator */
/*   const size_t num_random_points,               ///< [in] Number of random points to generate */
/*   gsl_matrix** random_points,                   ///< [in/out] Pointer to matrix of random points */
/*   gsl_vector_ulong** nearest_indices,           ///< [in/out] Pointer to vector of indices of nearest lattice point */
/*   gsl_vector** nearest_distances,               ///< [in/out] Pointer to vector of distances to nearest lattice point */
/*   gsl_matrix** workspace                        ///< [in/out] Pointer to workspace matrix for computing distances */
/*   ); */

#ifdef __cplusplus
}
#endif

#endif
