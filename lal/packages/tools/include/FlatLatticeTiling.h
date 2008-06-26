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
#include <lal/GSLSupport.h>

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(FLATLATTICETILINGH, "$Id$");

  /**
   * Types of parameter space bounds
   */
  typedef enum tagFlatLatticeTilingBoundType {
    FLT_BT_Undefined,
    FLT_BT_Singular,
    FLT_BT_Polynomial
  } FlatLatticeTilingBoundType;

  /**
   * Description of parameter space bounds
   */
  typedef struct tagFlatLatticeTilingBound {

    /* Dimension on which bound applies */
    INT4 dimension;

    /* Zone within dimension on which bound applies */
    INT4 zone;

    /* Type of bound */
    FlatLatticeTilingBoundType type;

    /* Singular bound variables */
    REAL8 singular_value;

    /* Polynomial bound variables */
    gsl_vector *poly_lower_const;
    gsl_matrix_int *poly_lower_exp;
    gsl_vector *poly_upper_const;
    gsl_matrix_int *poly_upper_exp;

  } FlatLatticeTilingBound;
  
  /**
   * State of the tiling
   */
  typedef enum tagFlatLatticeTilingState {
    FLT_S_NotInitialised,
    FLT_S_NotStarted,
    FLT_S_InProgress,
    FLT_S_Finished
  } FlatLatticeTilingState;

  /**
   * Information for the flat lattice tiling algorithm 
   */
  typedef struct tagFlatLatticeTiling {

    /* Dimension of the parameter space */
    INT4 dimensions;

    /* Parameter space bounds */
    INT4 num_bounds;
    FlatLatticeTilingBound *bounds;

    /* Index of current bound zone */
    gsl_vector_int *bound_zone;

    /* Maximum number of bound zones */
    gsl_vector_int *max_bound_zone;

    /* Metric of the parameter space in normalised coordinates */
    gsl_matrix *norm_metric;

    /* Conversion from normalised to real parameter coordinates */
    gsl_vector *norm_to_real;

    /* Maximum metric mismatch between the templates */
    REAL8 max_mismatch;

    /* Reduced dimension (singular dimensions excluded) */
    INT4 reduced_dims;

    /* Dimension map to/from reduced and full dimensions */
    gsl_vector_int *reduced_map;
    gsl_vector_int *dimension_map;

    /* Transformation to/from lattice coordinates and normalised space */
    gsl_matrix *latt_to_norm;
    gsl_matrix *norm_to_latt;

    /* Metric of the parameter space in lattice coordinates */
    gsl_matrix *latt_metric;

    /* Current point */
    gsl_vector_int *latt_current;
    gsl_vector *norm_current;

    /* Bounds on current point */
    gsl_vector *norm_lower;
    gsl_vector *norm_upper;

    /* Padding of bounds along each dimension */
    gsl_vector *padding;

    /* Current template point */
    gsl_vector *current;

    /* Total number of points generated so far */
    INT8 count;

    /* State of the tiling */
    FlatLatticeTilingState state;

  } FlatLatticeTiling;

  /**
   * Core functions 
   */
  FlatLatticeTiling* XLALCreateFlatLatticeTiling(INT4);
  void XLALFreeFlatLatticeTiling(FlatLatticeTiling*);
  int XLALSetFlatLatticeTilingMetric(FlatLatticeTiling*, gsl_matrix*, REAL8, gsl_vector*);
  int XLALAddFlatLatticeTilingSingularBound(FlatLatticeTiling*, INT4, INT4, REAL8);
  int XLALAddFlatLatticeTilingPolynomialBound(FlatLatticeTiling*, INT4, INT4, gsl_vector*, gsl_matrix_int*, gsl_vector*, gsl_matrix_int*);
  int XLALFinaliseFlatLatticeTilingBounds(FlatLatticeTiling*);
  int XLALSetFlatTilingLattice(FlatLatticeTiling*, gsl_matrix*, REAL8);
  int XLALIsPointInFlatLatticeParamSpace(FlatLatticeTiling*, gsl_vector*);
  int XLALNextFlatLatticePoint(FlatLatticeTiling*);
  REAL8 XLALTotalFlatLatticePointCount(FlatLatticeTiling*);
  
  /**
   * Supplementary function
   */
  gsl_vector* XLALMetricEllipseBoundingBox(gsl_matrix*, REAL8);
  int XLALOrthonormaliseWRTMetric(gsl_matrix*, gsl_matrix*);
  gsl_matrix* XLALSquareLowerTriangularLatticeGenerator(gsl_matrix*);
  int XLALNormaliseLatticeGenerator(gsl_matrix*, REAL8, REAL8);
  
  /**
   * Specific lattices and parameter spaces
   */
  int XLALSetFlatTilingCubicLattice(FlatLatticeTiling*);
  int XLALSetFlatTilingAnstarLattice(FlatLatticeTiling*);
  int XLALFlatLatticeTilingSquareParamSpace(FlatLatticeTiling**, gsl_vector*);
  
#ifdef __cplusplus
}
#endif

#endif
