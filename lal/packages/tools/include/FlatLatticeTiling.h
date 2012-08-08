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

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALDatatypes.h>
#include <lal/Random.h>
#include <lal/GSLSupport.h>

#ifdef __cplusplus
extern "C" <% // {
#endif

///
/// Flat lattice tiling bound function
///
typedef void (*FlatLatticeTilingBound)(
  double* lower,	///< [out] Lower bound on point in dimension
  double* upper,	///< [out] Upper bound on point in dimension
  gsl_vector* point,	///< [in] Point on which to find bounds
  gsl_vector* data	///< [in] Arbitrary data describing parameter space
  );

///
/// Flat tiling lattice generator function
///
typedef int (*FlatTilingLatticeGenerator)(
  size_t dimensions,		///< [in] Number of dimensions
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
  size_t dimensions		///< [in] Number of parameter space dimensions
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
size_t XLALGetFlatLatticeTilingDimensions(
  FlatLatticeTiling* tiling	///< [in] Tiling state
  );

///
/// Return the flat lattice tiling metric
///
gsl_matrix* XLALGetFlatLatticeTilingMetric(
  FlatLatticeTiling* tiling	///< [in] Tiling state
  );

///
/// Return the current number of flat lattice tiling parameter space points
///
uint64_t XLALGetFlatLatticePointCount(
  FlatLatticeTiling* tiling	///< [in] Tiling state
  );

///
/// Add a parameter space bound to the flat lattice tiling
///
int XLALAddFlatLatticeTilingBound(
  FlatLatticeTiling* tiling,	///< [in] Tiling state
  FlatLatticeTilingBound func,	///< [in] Parameter space bound function
  gsl_vector* data		///< [in] Arbitrary data describing parameter space
  );


///
/// Set the flat lattice tiling metric and maximum mismatch
///
int XLALSetFlatLatticeTilingMetric(
  FlatLatticeTiling* tiling,	///< [in] Tiling state
  gsl_matrix* metric,		///< [in] Parameter space metric
  double max_mismatch,		///< [in] Maximum prescribed mismatch
  gsl_vector* real_scale	///< [in] Multiply to get real metric, may be NULL
  );

///
/// Set the flat tiling lattice generator
///
int XLALSetFlatTilingLatticeGenerator(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  FlatTilingLatticeGenerator generator	///< [in] Lattice generator function
  );


///
/// Return the next point in the flat lattice tiling parameter space
///
gsl_vector* XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling		///< [in] Tiling state
  );

///
/// Calculate the total number of flat lattice tiling parameter space points
///
uint64_t XLALCountTotalFlatLatticePoints(
  FlatLatticeTiling* tiling		///< [in] Tiling state
  );

int XLALRandomPointInFlatLatticeParamSpace(FlatLatticeTiling*, RandomParams*, gsl_vector*, gsl_vector*, double*);

///
/// Calculate the generator matrix for a cubic (\f$Z_n\f$) lattice
///
int XLALCubicLatticeGenerator(
  size_t dimensions,		///< [in] Number of dimensions
  gsl_matrix** generator,	///< [out] Generator matrix
  double* norm_thickness	///< [out] Normalised thickness
  );

///
/// Calculate the generator matrix for a \f$A_n^*\f$ lattice
///
int XLALAnstarLatticeGenerator(
  size_t dimensions,		///< [in] Number of dimensions
  gsl_matrix** generator,	///< [out] Generator matrix
  double* norm_thickness	///< [out] Normalised thickness
  );

///
/// Add a constant parameter space bound to the flat lattice tiling
///
int XLALAddFlatLatticeTilingConstantBound(
  FlatLatticeTiling* tiling,	///< [in] Tiling state
  double lower,			///< [in] Lower bound on dimension
  double upper			///< [in] Upper bound on dimension
  );

#ifdef __cplusplus
%> // }
#endif

#endif

/// @}
