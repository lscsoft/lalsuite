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
/// Type of lattice to generate flat tiling with
///
typedef enum tagLatticeType {
  LATTICE_TYPE_CUBIC,				///< Cubic (\f$Z_n\f$) lattice
  LATTICE_TYPE_ANSTAR,				///< An-star (\f$A_n^*\f$) lattice
  LATTICE_TYPE_MAX
} LatticeType;

///
/// Flat lattice tiling state structure
///
typedef struct tagLatticeTiling LatticeTiling;

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
/// Set the parameter-space bounds on dimension \f$n\f$ of the lattice tiling.
///
/// The lower/upper bound \f$X_n\f$ is specified by
/// \f[
/// X_n - a_n = \sum_{k=0}^{P-1} c_{N,k} \left[ \sum_{j=0}^{N-1} c_{j,k} \prod_{i=0}^{n-1} (x_i - a_i)^{m_{nj+i,k}} \right]^{m_{nN,k}}
/// \f]
/// where \f$x = (x_0,\cdots,x_{n-1})\f$ is the current parameter-space point, \f$a\f$ is a vector of \f$n + 1\f$ offsets,
/// \f$c\f$ is a matrix of \f$N + 1 \times P\f$ coefficients, and \f$m\f$ is a matrix of \f$nM + 1 \times P\f$ exponents.
/// Note that \f$x^m \rightarrow 0\f$ if the result is not finite and real-valued.
///
int XLALSetLatticeBound(
  LatticeTiling* tiling,			///< [in] Tiling state
  const size_t dimension,			///< [in] Dimension on which bound applies (\f$n\f$)
  const gsl_vector* a,				///< [in] Vector of offsets (\f$a\f$)
  const gsl_matrix* c_lower,			///< [in] Matrix of coefficients (\f$c\f$) for the lower bound
  const gsl_matrix* m_lower,			///< [in] Matrix of exponents (\f$m\f$) for the lower bound
  const gsl_matrix* c_upper,			///< [in] Matrix of coefficients (\f$c\f$) for the upper bound
  const gsl_matrix* m_upper			///< [in] Matrix of exponents (\f$m\f$) for the upper bound
  );

///
/// Set a constant parameter-space bound, given by the minimum and
/// maximum of the two supplied bounds, on the lattice tiling
///
int XLALSetLatticeConstantBound(
  LatticeTiling* tiling,			///< [in] Tiling state
  const size_t dimension,			///< [in] Dimension on which bound applies
  const double bound1,				///< [in] First bound on dimension
  const double bound2				///< [in] Second bound on dimension
  );

///
/// Set elliptical bounds in two dimensions on the lattice tiling
///
int XLALSetLatticeEllipticalBounds(
  LatticeTiling* tiling,			///< [in] Tiling state
  const size_t dimension,			///< [in] Dimension of X bound (Y bound is one higher)
  const double x_centre,			///< [in] X centre of ellipse
  const double y_centre,			///< [in] Y centre of ellipse
  const double x_semi,				///< [in] Length of X semi-diameter
  const double y_semi				///< [in] Length of Y semi-diameter
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
/// the parameter space, so that then calling XLALNextFlatticePoint() will
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

#ifdef __cplusplus
}
#endif

#endif
