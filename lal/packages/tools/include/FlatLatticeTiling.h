/*
 *  Copyright (C) 2007, 2008 Karl Wette
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
 * \file
 * \brief Lattice-based template generation for flat metric parameter spaces
 */

#ifndef _FLATLATTICETILING_H
#define _FLATLATTICETILING_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALRCSID.h>
#include <lal/LALDatatypes.h>
#include <lal/Random.h>
#include <lal/GSLSupport.h>

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(FLATLATTICETILINGH, "$Id$");

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
  typedef struct tagFlatLatticeTilingBound {

    /* Number of bound dimensions */
    INT4 dimensions;

    /* Dimensions which are bound */
    UINT8 is_bound;

    /* Parameter space bound function */
    FlatLatticeTilingBoundFunc func;

    /* Arbitrary data describing parameter space */
    void *data;

    /* Cleanup function */
    FlatLatticeTilingBoundFree free;

  } FlatLatticeTilingBound;

  /**
   * Flat lattice tiling subspace
   */
  typedef struct tagFlatLatticeTilingSubspace {

    /* Total number of tiled (non-flat) dimensions */
    INT4 dimensions;

    /* Dimensions which are tiled (non-flat) */
    UINT8 is_tiled;

    /* Padding of bounds along each dimension */
    gsl_vector *padding;

    /* Increment vectors of the lattice tiling generator */
    gsl_matrix *increment;

  } FlatLatticeTilingSubspace;

  /**
   * Flat tiling lattice generator
   */
  typedef int (*FlatTilingLatticeGenerator)(
					    INT4 dimensions,        /**< Number of dimensions */
					    gsl_matrix** generator, /**< Generator matrix */
					    REAL8* norm_thickness   /**< Normalised thickness */
					    );

  /**
   * State of the flat lattice tiling algorithm
   */
  typedef enum tagFlatLatticeTilingState {
    FLT_S_NotInitialised,
    FLT_S_NotStarted,
    FLT_S_InProgress,
    FLT_S_Finished
  } FlatLatticeTilingState;
  typedef struct tagFlatLatticeTiling {

    /* Dimension of the parameter space */
    INT4 dimensions;

    /* Parameter space bounds */
    INT4 num_bounds;
    FlatLatticeTilingBound **bounds;
    gsl_vector_int *bound_map;
    gsl_vector *bound_point;

    /* Metric of the parameter space in normalised coordinates */
    gsl_matrix *metric;

    /* Normalised to real parameter coordinates scaling and offset */
    gsl_vector *real_scale;
    gsl_vector *real_offset;

    /* Maximum metric mismatch between the templates */
    REAL8 max_mismatch;

    /* Flat tiling lattice generator */
    FlatTilingLatticeGenerator generator;

    /* Cache of generated tiling subspaces */
    INT4 num_subspaces;
    FlatLatticeTilingSubspace **subspaces;

    /* Scaling of the padding of bounds (for testing) */
    REAL8 scale_padding;

    /* Current dimensions which are tiled (non-flat) */
    UINT8 curr_is_tiled;

    /* Current tiling subspace */
    FlatLatticeTilingSubspace *curr_subspace;

    /* Current lattice point */
    gsl_vector *curr_point;

    /* Bounds on current point */
    gsl_vector *curr_lower;
    gsl_vector *curr_upper;

    /* Current template */
    gsl_vector *current;

    /* Total number of points generated so far */
    UINT4 count;

    /* State of the tiling */
    FlatLatticeTilingState state;

  } FlatLatticeTiling;

  /**
   * Core functions
   */
  FlatLatticeTiling* XLALCreateFlatLatticeTiling(INT4);
  void XLALFreeFlatLatticeTiling(FlatLatticeTiling*);
  int XLALAddFlatLatticeTilingBound(FlatLatticeTiling*, UINT8, FlatLatticeTilingBoundFunc, void*, FlatLatticeTilingBoundFree);
  int XLALSetFlatLatticeTilingMetric(FlatLatticeTiling*, gsl_matrix*, REAL8, gsl_vector*);
  int XLALSetFlatTilingLattice(FlatLatticeTiling*, FlatTilingLatticeGenerator);
  int XLALNextFlatLatticePoint(FlatLatticeTiling*);
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
}
#endif

#endif
