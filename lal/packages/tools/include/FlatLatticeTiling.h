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
  double* lower,		///< [out] Lower bound on point in dimension
  double* upper,		///< [out] Upper bound on point in dimension
  const gsl_vector* point,	///< [in] Point on which to find bounds
  const void* data		///< [in] Arbitrary data describing parameter space
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

#ifdef SWIG // SWIG interface directives
SWIGLAL(NO_NEW_OBJECT(XLALNextFlatLatticePoint));
SWIGLAL(FUNCTION_POINTER(XLALCubicLatticeGenerator));
SWIGLAL(FUNCTION_POINTER(XLALAnstarLatticeGenerator));
#endif

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
uint64_t XLALGetFlatLatticePointCount(
  FlatLatticeTiling* tiling	///< [in] Tiling state
  );

///
/// Add a parameter space bound to the flat lattice tiling
///
int XLALSetFlatLatticeBound(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const size_t dimension,		///< [in] Dimension on which bound applies
  const FlatLatticeBound func,		///< [in] Parameter space bound function
  void* data				///< [in] Arbitrary data describing parameter space
  );

///
/// Set the flat lattice tiling metric and maximum mismatch
///
int XLALSetFlatLatticeMetric(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const gsl_matrix* metric,		///< [in] Parameter space metric
  const double max_mismatch		///< [in] Maximum prescribed mismatch
  );

///
/// Set the flat tiling lattice generator
///
int XLALSetFlatLatticeGenerator(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  const FlatLatticeGenerator generator	///< [in] Lattice generator function
  );

///
/// Return the next point in the flat lattice tiling parameter space
///
gsl_vector* XLALNextFlatLatticePoint(
  FlatLatticeTiling* tiling		///< [in] Tiling state
  );

///
/// Return a set of points in the flat lattice tiling parameter space
///
size_t XLALNextFlatLatticePoints(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  gsl_matrix* points,			///< [in] Flat lattice tiling points
  const bool fill_last			///< [in] If not enought points to fill 'points', whether to fill in using the last tiling point
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
uint64_t XLALCountTotalFlatLatticePoints(
  FlatLatticeTiling* tiling		///< [in] Tiling state
  );

///
/// Generate random points within the flat lattice tiling parameter space
///
int XLALGenerateRandomFlatLatticePoints(
  FlatLatticeTiling* tiling,		///< [in] Tiling state
  RandomParams* randpar,		///< [in] Random number generator state
  gsl_matrix* randpoints		///< [in] Random points (column-wise)
  );

///
/// Calculate the generator matrix for a cubic (\f$Z_n\f$) lattice
///
int XLALCubicLatticeGenerator(
  const size_t dimensions,	///< [in] Number of dimensions
  gsl_matrix** generator,	///< [out] Generator matrix
  double* norm_thickness	///< [out] Normalised thickness
  );

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
/// Workspace for computing the nearest template to a set of injections
///
typedef struct tagNearestTemplateWorkspace NearestTemplateWorkspace;

///
/// Create a new workspace for computing the nearest template to a set of injections
///
NearestTemplateWorkspace* XLALCreateNearestTemplateWorkspace(
  const gsl_matrix* metric,		///< [in] Parameter space metric
  const size_t num_templates,		///< [in] Number of templates to compare at once
  const size_t num_injections		///< [in] Nunber of injections to compare at once
  );

///
/// Destroy a nearest template workspace
///
void XLALDestroyNearestTemplateWorkspace(
  NearestTemplateWorkspace* wksp	///< [in] Nearest template workspace
  );

///
/// Update the templates used to compute distances from, and reset nearest template index
///
int XLALUpdateWorkspaceTemplates(
  NearestTemplateWorkspace* wksp,	///< [in] Nearest template workspace
  const gsl_matrix* templates,		///< [in] Template bank
  gsl_vector_uint* nearest_template	///< [in] Index of nearest template
  );

///
/// Update the injections used to compute distances to, and reset minimum distances
///
int XLALUpdateWorkspaceInjections(
  NearestTemplateWorkspace* wksp,	///< [in] Nearest template workspace
  const gsl_matrix* injections,		///< [in] Injection set
  gsl_vector* min_distance		///< [in] Distance from injection to nearest template
  );

int XLALUpdateNearestTemplateToInjections(
  NearestTemplateWorkspace* wksp,	///< [in] Nearest template workspace
  gsl_vector* min_distance,		///< [in] Distance from injection to nearest template
  gsl_vector_uint* nearest_template	///< [in] Index of nearest template
  );

#ifdef __cplusplus
}
#endif

#endif

/// @}
