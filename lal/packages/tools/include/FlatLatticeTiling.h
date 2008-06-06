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
#include <lal/FlatLatticeTilingSupport.h>

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(FLATLATTICETILINGH, "$Id$");

  /**
   * Types of parameter space bounds
   */
  enum tagFlatLatticeTilingParamSpaceBoundType {
    FLT_PSBT_Undefined,
    FLT_PSBT_Singular,
    FLT_PSBT_Linear,
    FLT_PSBT_Quadratic
  };
  
  /**
   * Information for the flat lattice tiling algorithm 
   */
  typedef struct tagFlatLatticeTiling {

    /* Dimension of the parameter space */
    INT4 dimension;

    /* Type of parameter space bound on each dimension */
    gsl_vector_int *bound_type;

    /* Singular bounds */
    gsl_vector *singular_bound;

    /* Linear bounds */
    gsl_matrix *linear_bound_A;
    gsl_vector *linear_bound_b;
    gsl_matrix *linear_vertices;
    gsl_vector_int *linear_map;
    gsl_vector *linear_min_vertex;
    gsl_vector *linear_max_vertex;
    SimplexMethodTableau *linear_tableau;
    gsl_matrix *linear_metric;
    gsl_vector *linear_current;
    gsl_vector *linear_transl_b;
    gsl_vector *linear_optimal;
    gsl_vector *linear_temp;

    /* Quadratic bounds */
    gsl_matrix *quadratic_bound_Q;
    gsl_vector_int *quadratic_map;

    /* Inverse lookup for bound index maps */
    gsl_vector_int *bound_inverse_map;

    /* Carefully check bounds only on these dimensions */
    gsl_vector_int *careful_linear;
    gsl_vector_int *careful_quadratic;

    /* Metric of the parameter space in normalised coordinates */
    gsl_matrix *norm_metric;

    /* Conversion from normalised to real parameter coordinates */
    gsl_vector *norm_to_real_mul;
    gsl_vector *norm_to_real_add;

    /* Maximum metric mismatch between the templates */
    REAL8 max_mismatch;

    /* Reduced dimension (singular dimensions excluded) */
    INT4 reduced_dim;

    /* Dimension map from reduced to full dimensions */
    gsl_vector_int *reduced_map;

    /* Transformation from lattice coordinates to normalised space */
    gsl_matrix *latt_to_norm;

    /* Padding along each dimension */
    gsl_vector *norm_padding;

    /* Current lattice point */
    gsl_vector_int *latt_current;
    gsl_vector *norm_current;

    /* Bounds on current point */
    gsl_vector *norm_lower;
    gsl_vector *norm_upper;

    /* Pointer to current point */
    gsl_vector *current_tile;

    /* Total number of points generated so far */
    INT8 tile_count;

    /* State of the tiling */
    int state;

  } FlatLatticeTiling;

  /**
   * Functions 
   */
  FlatLatticeTiling *XLALCreateFlatLatticeTiling(INT4);
  void XLALFreeFlatLatticeTiling(FlatLatticeTiling*);
  int XLALSetFlatLatticeTilingSingularBound(FlatLatticeTiling*, INT4, REAL8);
  int XLALAddFlatLatticeTilingLinearBound(FlatLatticeTiling*, INT4, gsl_vector*, REAL8, BOOLEAN);
  int XLALSetFlatLatticeTilingMetric(FlatLatticeTiling*, gsl_matrix*, REAL8, gsl_vector*);
  int XLALSetCubicTilingLattice(FlatLatticeTiling*);
  int XLALSetAnstarTilingLattice(FlatLatticeTiling*);
  int XLALNextFlatLatticeTile(FlatLatticeTiling*);
  REAL8 XLALTotalFlatLatticeTileCount(FlatLatticeTiling*);

#ifdef __cplusplus
}
#endif

#endif
