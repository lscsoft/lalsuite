/*
 *  Copyright (C) 2007, 2008, 2012 Karl Wette
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/**
 * \author Karl Wette
 * \ingroup pulsarTODO
 * \brief Lattice-based template generation for flat metric parameter spaces
 */

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

/**
 * Flat lattice tiling bound
 */
typedef BOOLEAN (*FlatLatticeTilingBoundFunc)(
  void *data,        /**< Arbitrary data describing parameter space */
  INT4 dimension,    /**< Dimension on which bound applies */
  gsl_vector *point, /**< Point on which to find bounds */
  REAL8 *lower,      /**< Lower bound on point in dimension */
  REAL8 *upper       /**< Upper bound on point in dimension */
  );
typedef void (*FlatLatticeTilingBoundFree)(
  void *data /**< Arbitrary data describing parameter space */
  );

/**
 * Flat tiling lattice generator
 */
typedef int (*FlatTilingLatticeGenerator)(
  INT4 dimensions,        /**< Number of dimensions */
  gsl_matrix** generator, /**< Generator matrix */
  REAL8* norm_thickness   /**< Normalised thickness */
  );

/**
 * Flat lattice tiling state/input structure
 */
typedef struct tagFlatLatticeTiling FlatLatticeTiling;

/**
 * Core functions
 */
FlatLatticeTiling* XLALCreateFlatLatticeTiling(INT4);
INT4 XLALFlatLatticeTilingDimension(FlatLatticeTiling*);
gsl_matrix* XLALFlatLatticeTilingMetric(FlatLatticeTiling*);
void XLALDestroyFlatLatticeTiling(FlatLatticeTiling*);
int XLALAddFlatLatticeTilingBound(FlatLatticeTiling*, UINT8, FlatLatticeTilingBoundFunc, void*, FlatLatticeTilingBoundFree);
int XLALSetFlatLatticeTilingMetric(FlatLatticeTiling*, gsl_matrix*, REAL8, gsl_vector*);
int XLALSetFlatTilingLattice(FlatLatticeTiling*, FlatTilingLatticeGenerator);
int XLALNextFlatLatticePoint(FlatLatticeTiling*);
gsl_vector* XLALCurrentFlatLatticePoint(FlatLatticeTiling*);
UINT4 XLALTotalFlatLatticePointCount(FlatLatticeTiling*);
int XLALRandomPointInFlatLatticeParamSpace(FlatLatticeTiling*, RandomParams*, gsl_vector*, gsl_vector*, REAL8*);

/**
 * Support functions
 */
gsl_matrix* XLALMetricEllipsePrincipalAxes(gsl_matrix*, REAL8);
gsl_vector* XLALMetricEllipseBoundingBox(gsl_matrix*, REAL8);
int XLALOrthonormaliseWRTMetric(gsl_matrix*, gsl_matrix*);
gsl_matrix* XLALSquareLowerTriangularLatticeGenerator(gsl_matrix*);
int XLALNormaliseLatticeGenerator(gsl_matrix*, REAL8, REAL8);

/**
 * Specific lattices and parameter spaces
 */
int XLALSetFlatTilingCubicLattice(FlatLatticeTiling*);
int XLALSetFlatTilingAnstarLattice(FlatLatticeTiling*);
int XLALAddFlatLatticeTilingConstantBound(FlatLatticeTiling*, INT4, REAL8, REAL8);

#ifdef __cplusplus
%> // }
#endif

#endif
