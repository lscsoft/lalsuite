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

#include <math.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_permutation.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_nan.h>

#include <lal/LALRCSID.h>
#include <lal/LALStdlib.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <lal/GSLSupport.h>
#include <lal/FlatLatticeTiling.h>

NRCSID(FLATLATTICETILINGC, "$Id$");

static const BOOLEAN TRUE  = (1==1);
static const BOOLEAN FALSE = (0==1);

/**
 * Create a new flat lattice tiling structure
 */
FlatLatticeTiling *XLALCreateFlatLatticeTiling(
	 				       INT4 dimensions /**< Dimensionality of the parameter space */
					       )
{

  FlatLatticeTiling *tiling = NULL;

  /* Check dimension */
  if (dimensions < 1)
    XLAL_ERROR_NULL("'dimensions' must be non-zero", XLAL_EINVAL);

  /* Allocate memory */
  if ((tiling = (FlatLatticeTiling*)LALMalloc(sizeof(FlatLatticeTiling))) == NULL)
    XLAL_ERROR_NULL("Could not allocate 'tiling'", XLAL_ENOMEM);

  /* Initialise structure */
  tiling->dimensions = dimensions;
  tiling->num_bounds = 0;
  tiling->bounds = NULL;
  tiling->bound_zone = NULL;
  tiling->max_bound_zone = NULL;
  tiling->norm_metric = NULL;
  tiling->norm_to_real = NULL;
  tiling->max_mismatch = 0.0;
  tiling->reduced_dims = -1;
  tiling->reduced_map = NULL;
  tiling->dimension_map = NULL;
  tiling->latt_to_norm = NULL;
  tiling->norm_to_latt = NULL;
  tiling->latt_metric = NULL;
  tiling->latt_current = NULL;
  tiling->norm_current = NULL;
  tiling->norm_lower = NULL;
  tiling->norm_upper = NULL;
  tiling->padding = NULL;
  tiling->current = NULL;
  tiling->count = 0;
  tiling->state = FLT_S_NotInitialised;

  return tiling;

}

/**
 * Free a flat lattice tiling structure
 */
void XLALFreeFlatLatticeTiling(
			       FlatLatticeTiling *tiling /**< Tiling structure */
			       )
{

  int k;
  
  if (tiling) {

    for (k = 0; k < tiling->num_bounds; ++k) {

      FREE_GSL_VECTOR(tiling->bounds[k].poly_lower_const);
      FREE_GSL_VECTOR(tiling->bounds[k].poly_upper_const);
      FREE_GSL_MATRIX_INT(tiling->bounds[k].poly_lower_exp);
      FREE_GSL_MATRIX_INT(tiling->bounds[k].poly_upper_exp);
      
    }
    LALFree(tiling->bounds);

    FREE_GSL_VECTOR_INT(tiling->bound_zone);
    FREE_GSL_VECTOR_INT(tiling->max_bound_zone);
    FREE_GSL_MATRIX(tiling->norm_metric);
    FREE_GSL_VECTOR(tiling->norm_to_real);
    FREE_GSL_VECTOR_INT(tiling->reduced_map);
    FREE_GSL_VECTOR_INT(tiling->dimension_map);
    FREE_GSL_MATRIX(tiling->latt_to_norm);
    FREE_GSL_MATRIX(tiling->norm_to_latt);
    FREE_GSL_MATRIX(tiling->latt_metric);
    FREE_GSL_VECTOR_INT(tiling->latt_current);
    FREE_GSL_VECTOR(tiling->norm_current);
    FREE_GSL_VECTOR(tiling->norm_lower);
    FREE_GSL_VECTOR(tiling->norm_upper);
    FREE_GSL_VECTOR(tiling->padding);
    FREE_GSL_VECTOR(tiling->current);

    LALFree(tiling);
  
  }
  
}

/**
 * Set the flat lattice tiling metric and maximum mismatch
 */
int XLALSetFlatLatticeTilingMetric(
				   FlatLatticeTiling *tiling, /**< Tiling structure */
				   gsl_matrix *metric,        /**< Metric */
				   REAL8 max_mismatch,        /**< Maximum mismatch */
				   gsl_vector *norm_to_real   /**< Multiply to get real metric, may be NULL */
				   )
{

  const int n = tiling->dimensions;
  
  int i, j;
  
  /* Check tiling state */
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR("'tiling' has already been initialised", XLAL_EFAILED);

  /* Check metric and mismatch */
  if (tiling->norm_metric)
    XLAL_ERROR("'tiling->norm_metric' has already been set", XLAL_EFAILED);
  if (!metric)
    XLAL_ERROR("'metric' must be allocated", XLAL_EINVAL);
  if (metric->size1 != metric->size2)
    XLAL_ERROR("'metric' must be square", XLAL_EINVAL);
  if (n != (int)metric->size1)
    XLAL_ERROR("'metric' size must match tiling dimension", XLAL_EINVAL);
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      if (gsl_matrix_get(metric, i, j) != gsl_matrix_get(metric, j, i))
	XLAL_ERROR("'metric' must be symmetric", XLAL_EINVAL);
  if (max_mismatch < 0.0)
    XLAL_ERROR("'max_mismatch' must be strictly positive", XLAL_EINVAL);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(tiling->norm_metric,     n, n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->norm_to_real,       n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->norm_current,       n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->norm_lower,         n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->norm_upper,         n, XLAL_ERROR);
  ALLOC_GSL_VECTOR(tiling->current,            n, XLAL_ERROR);
  ALLOC_GSL_VECTOR_INT(tiling->bound_zone,     n, XLAL_ERROR);
  ALLOC_GSL_VECTOR_INT(tiling->max_bound_zone, n, XLAL_ERROR);

  /* Copy metric, normalised to real conversion, and mismatch */
  gsl_matrix_memcpy(tiling->norm_metric, metric);
  if (norm_to_real)
    gsl_vector_memcpy(tiling->norm_to_real, norm_to_real);
  else
    gsl_vector_set_all(tiling->norm_to_real, 1.0);
  tiling->max_mismatch = max_mismatch;

  /* Initialise maximum bound zones */
  gsl_vector_int_set_all(tiling->max_bound_zone, -1);

  return XLAL_SUCCESS;

}

/**
 * Add a bound to the parameter space bounds
 */
static FlatLatticeTilingBound *AddFlatLatticeTilingBound(
							 FlatLatticeTiling *tiling,      /**< Tiling structure */
							 INT4 dimension,                 /**< Dimension on which bound applies */
							 INT4 zone,                      /**< Zone within dimension on which bound applies */
							 FlatLatticeTilingBoundType type /**< Type of bound */
							 )
{

  const int n = tiling->dimensions;

  FlatLatticeTilingBound *last = NULL;

  /* Check tiling state */
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR_NULL("'tiling' has already been initialised", XLAL_EFAILED);
  if (!tiling->norm_to_real)
    XLAL_ERROR_NULL("'XLALSetFlatLatticeTilingMetric' has not been called", XLAL_EFAILED);

  /* Check input */
  if (dimension < 0 || n <= dimension)
    XLAL_ERROR_NULL("'dimension' is out of range", XLAL_EFAILED);
  if (zone < 0)
    XLAL_ERROR_NULL("'zone' is out of range", XLAL_EFAILED);
  if (zone - gsl_vector_int_get(tiling->max_bound_zone, dimension) > 1)
    XLAL_ERROR_NULL("'zone' is discontinuous", XLAL_EFAILED);  

  /* Increment number of bounds */
  ++tiling->num_bounds;

  /* (Re)Allocate memory */
  if ((tiling->bounds = (FlatLatticeTilingBound*)LALRealloc(tiling->bounds, tiling->num_bounds * sizeof(FlatLatticeTilingBound))) == NULL)
    XLAL_ERROR_NULL("Could not (re)allocate 'tiling->bounds'", XLAL_ENOMEM);

  /* Initialise last element */
  last = &tiling->bounds[tiling->num_bounds - 1];
  last->dimension = dimension;
  last->zone = zone;
  last->type = type;
  last->singular_value = 0.0;
  last->poly_lower_const = NULL;
  last->poly_lower_exp = NULL;
  last->poly_upper_const = NULL;
  last->poly_upper_exp = NULL;

  /* Initialise maximum bound zones */
  gsl_vector_int_set(tiling->max_bound_zone, dimension, zone);

  return last;

}

/**
 * Add a singular bound to the parameter space
 */
int XLALAddFlatLatticeTilingSingularBound(
					  FlatLatticeTiling *tiling, /**< Tiling structure */
					  INT4 dimension,            /**< Dimension on which bound applies */
					  INT4 zone,                 /**< Zone within dimension on which bound applies */
					  REAL8 value                /**< Value of the singular bound */
					  )
{

  FlatLatticeTilingBound *bound = NULL;

  /* Add bound to tiling */
  if ((bound = AddFlatLatticeTilingBound(tiling, dimension, zone, FLT_BT_Singular)) == NULL)
    XLAL_ERROR("AddFlatLatticeTilingBound failed", XLAL_EFAILED);

  /* Copy singular bound */
  bound->singular_value = value;

  /* Normalise singular bound */
  bound->singular_value /= gsl_vector_get(tiling->norm_to_real, dimension);

  return XLAL_SUCCESS;

}

/**
 * Add a polynomial bound to the parameter space
 */
int XLALAddFlatLatticeTilingPolynomialBound(
					    FlatLatticeTiling *tiling, /**< Tiling structure */
					    INT4 dimension,            /**< Dimension on which bound applies */
					    INT4 zone,                 /**< Zone within dimension on which bound applies */
					    gsl_vector *lower_const,   /**< Constants for each term of lower polynomial */
					    gsl_matrix_int *lower_exp, /**< Exponents (columns) for each term of lower polynomial */
					    gsl_vector *upper_const,   /**< Constants for each term of upper polynomial */
					    gsl_matrix_int *upper_exp  /**< Exponents (columns) for each term of upper polynomial */
					    )
{

  FlatLatticeTilingBound *bound = NULL;

  /* Check input */
  if (!lower_const)
    XLAL_ERROR("'lower_const' must be allocated", XLAL_EFAILED);
  if (!upper_const)
    XLAL_ERROR("'upper_const' must be allocated", XLAL_EFAILED);
  if (lower_const->size != upper_const->size)
    XLAL_ERROR("'lower_const' and 'upper_const' are different sizes", XLAL_EFAILED);
  if (lower_exp) {
    if (lower_exp->size2 != lower_const->size)
      XLAL_ERROR("'lower_exp' and 'lower_const' are inconsistent sizes", XLAL_EFAILED);
    if (lower_exp->size1 != (size_t)dimension)
      XLAL_ERROR("'lower_exp' has incorrect number of rows", XLAL_EFAILED);
  }
  if (upper_exp) {
    if (upper_exp->size2 != upper_const->size)
      XLAL_ERROR("'upper_exp' and 'upper_const' are inconsistent sizes", XLAL_EFAILED);
    if (upper_exp->size1 != (size_t)dimension)
      XLAL_ERROR("'lower_exp' has incorrect number of rows", XLAL_EFAILED);
  }

  /* Add bound to tiling */
  if ((bound = AddFlatLatticeTilingBound(tiling, dimension, zone, FLT_BT_Polynomial)) == NULL)
    XLAL_ERROR("AddFlatLatticeTilingBound failed", XLAL_EFAILED);

  /* Copy polynomial bound */
  ALLOC_GSL_VECTOR(bound->poly_lower_const, lower_const->size, XLAL_ERROR);
  gsl_vector_memcpy(bound->poly_lower_const, lower_const);
  if (lower_exp) {
    ALLOC_GSL_MATRIX_INT(bound->poly_lower_exp, lower_exp->size1, lower_exp->size1, XLAL_ERROR);
    gsl_matrix_int_memcpy(bound->poly_lower_exp, lower_exp);
  }
  ALLOC_GSL_VECTOR(bound->poly_upper_const, upper_const->size, XLAL_ERROR);
  gsl_vector_memcpy(bound->poly_upper_const, upper_const);
  if (upper_exp) {
    ALLOC_GSL_MATRIX_INT(bound->poly_upper_exp, upper_exp->size1, upper_exp->size1, XLAL_ERROR);
    gsl_matrix_int_memcpy(bound->poly_upper_exp, upper_exp);
  }

  /* Normalise polynomial bound */
  gsl_vector_scale(bound->poly_lower_const, 1.0 / gsl_vector_get(tiling->norm_to_real, dimension));
  gsl_vector_scale(bound->poly_upper_const, 1.0 / gsl_vector_get(tiling->norm_to_real, dimension));

  return XLAL_SUCCESS;

}

/**
 * Check the parameter space bounds and get the reduced dimension
 */
int XLALFinaliseFlatLatticeTilingBounds(
					FlatLatticeTiling *tiling /**< Tiling structure */
					)
{
  
  const int N = tiling->num_bounds;
  const int n = tiling->dimensions;
  
  int in, k;
  gsl_vector_int *is_singular = NULL;
  
  /* Check tiling state */
  if (tiling->state != FLT_S_NotInitialised)
    XLAL_ERROR("'tiling' has already been initialised", XLAL_EFAILED);
  if (!tiling->norm_to_real)
    XLAL_ERROR("'XLALSetFlatLatticeTilingMetric' has not been called", XLAL_EFAILED);

  /* Allocate memory */
  ALLOC_GSL_VECTOR_INT(is_singular, n, XLAL_ERROR);
  
  /* Check parameter space bounds */
  if (tiling->num_bounds == 0 || !tiling->bounds)
    XLAL_ERROR("No parameter space bounds have been added", XLAL_EFAILED);
  for (in = 0; in < n; ++in)
    if (gsl_vector_int_get(tiling->max_bound_zone, in) < 0)
      XLAL_ERROR("No all dimensions have been bounded", XLAL_EFAILED);
  
  /* Check for singular dimensions */
  gsl_vector_int_set_all(is_singular, TRUE);
  for (k = 0; k < N; ++k)
    if (tiling->bounds[k].type != FLT_BT_Singular)
      gsl_vector_int_set(is_singular, tiling->bounds[k].dimension, FALSE);
  tiling->reduced_dims = n;
  for (in = 0; in < n; ++in)
    if (gsl_vector_int_get(is_singular, in))
      --tiling->reduced_dims;
  
  {
    const int r = tiling->reduced_dims;

    int ir;

    /* Allocate memory */
    if (r > 0)
      ALLOC_GSL_VECTOR_INT(tiling->reduced_map,   r, XLAL_ERROR);
    ALLOC_GSL_VECTOR_INT(tiling->dimension_map,   n, XLAL_ERROR);
    
    /* Make maps to/from reduced and full dimensions */
    gsl_vector_int_set_all(tiling->dimension_map, -1);
    for (ir = 0, in = 0; ir < r; ++ir, ++in) {
      while (gsl_vector_int_get(is_singular, in))
	++in;
      if (r > 0)
	gsl_vector_int_set(tiling->reduced_map, ir, in);
      gsl_vector_int_set(tiling->dimension_map, in, ir);
    }
    
  }
  
  /* Cleanup */
  FREE_GSL_VECTOR_INT(is_singular);
  
  return XLAL_SUCCESS;

}

/**
 * Calculate the polynomial bound
 */
static REAL8 CalculateBoundPolynomial(
				      gsl_vector *point,       /**< Point */
				      gsl_vector *poly_const,  /**< Constants */
				      gsl_matrix_int *poly_exp /**< Exponents */
				      )
{
  
  int i, j;  
  double bound = 0.0;

  for (j = 0; j < (int)poly_const->size; ++j) {
    double x = gsl_vector_get(poly_const, j);
    if (point && poly_exp)
      for (i = 0; i < (int)poly_exp->size1; ++i)
	x *= pow(gsl_vector_get(point, i),
		 gsl_matrix_int_get(poly_exp, i, j));
    bound += x;
  }

  return bound;

}

/**
 * Get bounds of the specified dimension/zone
 */
static int GetBounds(
		     FlatLatticeTiling *tiling, /**< Tiling structure */
		     INT4 dimension,            /**< Dimension on which bound applies */
		     INT4 zone,                 /**< Zone within dimension on which bound applies */
		     gsl_vector *point,         /**< Point on which to find bounds */
		     REAL8 *lower,              /**< Lower bound on dimension */
		     REAL8 *upper               /**< Upper bound on dimension */
		     )
{

  const int N = tiling->num_bounds;

  int k;
  gsl_vector_view partial_view;
  gsl_vector *partial = NULL;

  /* Check dimension and zone */
  if (dimension < 0 || tiling->dimensions <= dimension)
    XLAL_ERROR("'dimension' is out of bounds", XLAL_EFAILED);
  if (zone < 0 || gsl_vector_int_get(tiling->max_bound_zone, dimension) < zone)
    XLAL_ERROR("'zone' is out of bounds", XLAL_EFAILED);

  /* Get partial view of point (should only use this part) */
  if (dimension > 0) {
    partial_view = gsl_vector_subvector(point, 0, dimension);
    partial = &partial_view.vector;
  }
  else
    partial = NULL;

  /* Iterate over all bounds */
  *lower = GSL_NEGINF;
  *upper = GSL_POSINF;
  for (k = 0; k < N; ++k) {

    /* Get current bound */
    const FlatLatticeTilingBound *this = &(tiling->bounds[k]);

    REAL8 lower_k, upper_k;

    /* Filter on dimension and zone */
    if (this->dimension != dimension || this->zone != zone)
      continue;

    /* Switch on boundary type */
    switch (this->type) {
      
    case FLT_BT_Singular:
      
      /* Singular boundary */
      lower_k = upper_k = this->singular_value;
      
      break;
      
    case FLT_BT_Polynomial:
      
      /* Calculate lower and upper bound polynomial */
      lower_k = CalculateBoundPolynomial(partial, this->poly_lower_const, this->poly_lower_exp);
      upper_k = CalculateBoundPolynomial(partial, this->poly_upper_const, this->poly_upper_exp);

      break;
      
    default:
      XLAL_ERROR("Unknown boundary type", XLAL_EFAILED);
      
    }
    
    /* Take the most exclusive bound */
    *lower = GSL_MAX_DBL(*lower, lower_k);
    *upper = GSL_MIN_DBL(*upper, upper_k);

  }

  /* Return if bounds were found and are sensible */
  if (!gsl_finite(*lower) || !gsl_finite(*upper))
    XLAL_ERROR("Bounds are nonsensical", XLAL_EFAILED);

  return XLAL_SUCCESS;

}

/**
 * Set the tiling lattice generator
 */
int XLALSetFlatTilingLattice(
			     FlatLatticeTiling *tiling, /**< Tiling structure */
			     gsl_matrix *generator,     /**< Lattice generator */
			     REAL8 norm_thickness       /**< Normalised thickness of lattice */
			     )
{




  
  const int r = tiling->reduced_dims;
  
  /* Check for non-singular dimensions */
  if (r > 0) {

    int i, j;  
    gsl_matrix *reduced_metric = NULL;
    gsl_matrix *ortho_directions = NULL;
    gsl_matrix *sq_lwtri_generator = NULL;
    gsl_matrix *temp = NULL;
    gsl_permutation *perm = NULL;
    int sign = 0;
    
    /* Allocate memory */
    ALLOC_GSL_MATRIX(tiling->latt_to_norm,     r, r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(tiling->norm_to_latt,     r, r, XLAL_ERROR); 
    ALLOC_GSL_MATRIX(tiling->latt_metric ,     r, r, XLAL_ERROR);
    ALLOC_GSL_VECTOR_INT(tiling->latt_current,    r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(reduced_metric,           r, r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(ortho_directions,         r, r, XLAL_ERROR);
    ALLOC_GSL_MATRIX(temp,                     r, r, XLAL_ERROR);
    ALLOC_GSL_PERMUTATION(perm,                   r, XLAL_ERROR);

    /* Extract reduced metric from normalised metric */
    for (i = 0; i < r; ++i) {
      const int in = gsl_vector_int_get(tiling->reduced_map, i);
      for (j = 0; j < r; ++j) {
	const int jn = gsl_vector_int_get(tiling->reduced_map, j);
	gsl_matrix_set(reduced_metric, i, j, gsl_matrix_get(tiling->norm_metric, in, jn));
      }
    }

    /* Use reduced metric ellipse bounding box as padding */
    if (NULL == (tiling->padding = XLALMetricEllipseBoundingBox(reduced_metric, tiling->max_mismatch)))
      XLAL_ERROR("XLALMetricEllipseBoundingBox failed", XLAL_EFAILED);

    /* Find orthonormal directions of reduced metric */
    gsl_matrix_set_identity(ortho_directions);
    if (XLAL_SUCCESS != XLALOrthonormaliseWRTMetric(ortho_directions, reduced_metric))
      XLAL_ERROR("XLALOrthonormaliseWRTMetric failed", XLAL_EFAILED);

    /* Transform lattice generator to square lower triangular */
    if (NULL == (sq_lwtri_generator = XLALSquareLowerTriangularLatticeGenerator(generator)))
      XLAL_ERROR("XLALSquareLowerTriangularLatticeGenerator failed", XLAL_EFAILED);
    
    /* Normalise lattice generator so covering radius is sqrt(mismatch) */
    if (XLAL_SUCCESS != XLALNormaliseLatticeGenerator(sq_lwtri_generator, norm_thickness, sqrt(tiling->max_mismatch)))
      XLAL_ERROR("XLALNormaliseLatticeGenerator failed", XLAL_EFAILED);
    
    /* Compute transformation from lattice coordinates to normalised parameter space */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, ortho_directions, sq_lwtri_generator, 0.0, tiling->latt_to_norm);

    /* Compute the inverse transformation for use in injections */
    gsl_matrix_memcpy(temp, tiling->latt_to_norm);
    gsl_linalg_LU_decomp(temp, perm, &sign);
    gsl_linalg_LU_invert(temp, perm, tiling->norm_to_latt);

    /* Compute the metric in lattice coordinates */
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, reduced_metric,       tiling->latt_to_norm, 0.0, temp);
    gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, tiling->latt_to_norm, temp,                 0.0, tiling->latt_metric);
    
    /* Cleanup */
    FREE_GSL_MATRIX(reduced_metric);
    FREE_GSL_MATRIX(ortho_directions);
    FREE_GSL_MATRIX(sq_lwtri_generator);
    FREE_GSL_MATRIX(temp);
    FREE_GSL_PERMUTATION(perm);
        
  }

  /* Tiling is now fully initialised */
  tiling->state = FLT_S_NotStarted;

  return XLAL_SUCCESS;

}

/**
 * Update the current point in normalised parameter space from 
 * the current point in lattice coordinates
 */
static void UpdateNormCurrent(
			      FlatLatticeTiling *tiling /**< Tiling structure */
			      )
{
  
  const int r = tiling->reduced_dims;

  int i, j;
  
  /* Compute the matrix transform for non-singular dimensions */
  gsl_vector_memcpy(tiling->norm_current, tiling->norm_lower);
  for (i = 0; i < r; ++i) {
    const int in = gsl_vector_int_get(tiling->reduced_map, i);
    double x = gsl_vector_get(tiling->norm_current, in);
    for (j = 0; j < r; ++j) {
      x += gsl_matrix_get(tiling->latt_to_norm, i, j) * 
	gsl_vector_int_get(tiling->latt_current, j);
    }
    gsl_vector_set(tiling->norm_current, in, x);
  }

}

/**
 * Return to the first point in specified and higher dimensions
 */
static int ReturnToFirstPoint(
			      FlatLatticeTiling *tiling, /**< Tiling structure */
			      INT4 dimension             /**< Starting dimension */
			      )
{
  
  const int n = tiling->dimensions;

  int in;
  REAL8 lower, upper;

  for (in = dimension; in < n; ++in) {
    
    /* Get reduced dimension */
    const int ir = gsl_vector_int_get(tiling->dimension_map, in);

    /* Get and store bounds */
    if (XLAL_SUCCESS != GetBounds(tiling, in, gsl_vector_int_get(tiling->bound_zone, in), tiling->norm_current, &lower, &upper))
      XLAL_ERROR("GetBounds failed", XLAL_EFAILED);
    gsl_vector_set(tiling->norm_lower, in, lower);
    gsl_vector_set(tiling->norm_upper, in, upper);
    
    /* Move the current point to the lower bound minus the padding */
    if (ir >= 0) {
      const int latt = gsl_vector_int_get(tiling->latt_current, ir);
      const double dnorm = (gsl_vector_get(tiling->norm_lower, in) -
			    gsl_vector_get(tiling->padding, ir) - 
			    gsl_vector_get(tiling->norm_current, in)) /
	gsl_matrix_get(tiling->latt_to_norm, ir, ir);
      const int dlatt = (int)(dnorm > 0.0 ? ceil(dnorm) : floor(dnorm));
      gsl_vector_int_set(tiling->latt_current, ir, latt + dlatt);
      UpdateNormCurrent(tiling);
    }
    
  }
  
  return XLAL_SUCCESS;
  
}

/**
 * Check whether the metric ellipse about the
 * current point intersects the parameter space.
 */
int XLALIsPointInFlatLatticeParamSpace(
				       FlatLatticeTiling *tiling, /**< Tiling structure */
				       gsl_vector *norm_point     /**< Normalised point */
				       )
{

  const int r = tiling->reduced_dims;

  int i;

  return XLAL_SUCCESS;

  /* Check whether point is in parameter space */
  for (i = 0; i < r; ++i) {
    const int in = gsl_vector_int_get(tiling->reduced_map, i);
    const double lower = gsl_vector_get(tiling->norm_lower, in);
    const double point = gsl_vector_get(norm_point, in);
    const double upper = gsl_vector_get(tiling->norm_upper, in);
    if (point < lower || upper < point) {
      return XLAL_FAILURE;
    }
  }

  return XLAL_SUCCESS;

}

/**
 * Move to the next point in the flat lattice tiling
 */
int XLALNextFlatLatticePoint(
			     FlatLatticeTiling *tiling /**< Tiling structure */
			     )
{
  
  const int n = tiling->dimensions;
  const int r = tiling->reduced_dims;

  /* Switch on tiling state */
  switch (tiling->state) {
    
  case FLT_S_NotInitialised:
    
    /* Fail if uninitialised */
    XLAL_ERROR("'tiling' has not been fully initialised", XLAL_EINVAL);
    
  case FLT_S_Finished:
    
    /* Exit if finished */
    return XLAL_FAILURE;
    
  case FLT_S_NotStarted:
    
    /* Initialise tiling to first point */
    tiling->count = 0;
    gsl_vector_int_set_zero(tiling->bound_zone);
    if (XLAL_SUCCESS != ReturnToFirstPoint(tiling, 0))
      XLAL_ERROR("ReturnToFirstPoint failed", XLAL_EFAILED);
    break;

  default:
    break;
    
  }
  
  /* Loop until a point is in the parameter space */
  do {
    
    int in = n - 1;

    /* Loop until a point is found */
    do {
      
/*       printf("# in: %i\n", in); */
/*       { */
/* 	int i; */
/* 	printf("# latt_current:"); */
/* 	for (i = 0; i < r; ++i) */
/* 	  printf(" %i", gsl_vector_int_get(tiling->latt_current, i)); */
/* 	printf("\n"); */
/* 	printf("# norm_current: "); */
/* 	for (i = 0; i < n; ++i) */
/* 	  printf(" %g", gsl_vector_get(tiling->norm_current, i)); */
/* 	printf("\n"); */
/* 	printf("# bound_zone/max_bound_zone: "); */
/* 	for (i = 0; i < n; ++i) */
/* 	  printf(" %i/%i", gsl_vector_int_get(tiling->bound_zone, i), gsl_vector_int_get(tiling->max_bound_zone, i)); */
/* 	printf("\n"); */
/*       } */
      
      /* If first point, test in loop condition */
      if (tiling->state == FLT_S_NotStarted) {
	tiling->state = FLT_S_InProgress;
	break;
      }
      
      /* Check for non-singular dimensions */
      if (r > 0) {
	
	/* Get reduced dimension */
	const int ir = gsl_vector_int_get(tiling->dimension_map, in);
	
	/* If this is a non-singular dimension */
	if (ir >= 0) {
	  
	  /* Move to next lattice point along this dimension */
	  {
	    const int latt = gsl_vector_int_get(tiling->latt_current, ir);
	    gsl_vector_int_set(tiling->latt_current, ir, latt + 1);
	    UpdateNormCurrent(tiling);
	  }
	  
	  /* Move point back to lower bound in higher dimensions */
	  if (XLAL_SUCCESS != ReturnToFirstPoint(tiling, in + 1))
	    XLAL_ERROR("ReturnToFirstPoint failed", XLAL_EFAILED);

	  /* If point is still in bounds */
	  {
	    const double point = gsl_vector_get(tiling->norm_current, in);
	    const double upper = gsl_vector_get(tiling->norm_upper, in);
	    const double padding = gsl_vector_get(tiling->padding, ir);
/* 	    printf("# point: %g upper: %g padding: %g\n", point, upper, padding); */
	    if (point <= upper + padding) {
	      
	      /* Test point in loop condition */
	      break;
	      
	    }
	  }
	  
	}
      }
      
      /* Move to next zone in this dimension */
      {
	const int zone = gsl_vector_int_get(tiling->bound_zone, in);
	gsl_vector_int_set(tiling->bound_zone, in, zone + 1);
      }
      
      /* If point is still in bounds */
      {
	const int zone = gsl_vector_int_get(tiling->bound_zone, in);
	const int max_zone = gsl_vector_int_get(tiling->max_bound_zone, in);
	if (zone <= max_zone) {
	  
	  /* Move point back to lower bound in higher dimensions */
	  if (XLAL_SUCCESS != ReturnToFirstPoint(tiling, in))
	    XLAL_ERROR("ReturnToFirstPoint failed", XLAL_EFAILED);
	  
	  /* Test point in loop condition */
	  break;
	  
	}
      }
      
      /* Reset zone to zero in this and higher dimensions */
      {
	gsl_vector_int_view v = gsl_vector_int_subvector(tiling->bound_zone, in, n - in);
	gsl_vector_int_set_zero(&v.vector);
      }

      /* Move to lower dimension */
      --in;
      
    } while (in >= 0);
	
    /* If this is the lowest dimension, we're done! */
    if (in < 0)
      return XLAL_FAILURE;

  } while (XLAL_SUCCESS != XLALIsPointInFlatLatticeParamSpace(tiling, tiling->norm_current));

  /* Found a point, increment count */
  ++tiling->count;

  /* Copy template and convert to real coordinates */
  gsl_vector_memcpy(tiling->current, tiling->norm_current);
  gsl_vector_mul(tiling->current, tiling->norm_to_real);
  
  return XLAL_SUCCESS;

}

/**
 * Return the count of the total number of flat lattice points
 */
REAL8 XLALTotalFlatLatticePointCount(
				     FlatLatticeTiling *tiling /**< Tiling structure */
				     )
{
  
  /* Switch on tiling state */
  switch (tiling->state) {
    
  case FLT_S_NotInitialised:

    /* Fail if uninitialised */
    XLAL_ERROR("'tiling' has not been fully initialised", XLAL_EINVAL);
    
  case FLT_S_NotStarted:
    {
      
      int retn;

      /* Iterate through all templates */
      while ((retn = XLALNextFlatLatticePoint(tiling)) == XLAL_SUCCESS);
      if (retn != XLAL_FAILURE)
	XLAL_ERROR_REAL8("XLALNextFlatLatticePoint failed", XLAL_EFAILED);
      
      /* Reset tiling */
      tiling->state = FLT_S_NotStarted;
      
    }
    break;

  default:
    break;

  }
    
  /* Return the current/final number of templates */
  return (REAL8) tiling->count;
  
}

/**
 * Find the bounding box of the mismatch ellipses of a metric
 */
gsl_vector *XLALMetricEllipseBoundingBox(
					 gsl_matrix *metric, /**< Metric to bound */
					 REAL8 max_mismatch  /**< Maximum mismatch w.r.t metric */
					 )
{

  const int n = metric->size1;

  int i;
  gsl_matrix *LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  gsl_matrix *inverse = NULL;
  int LU_sign = 0;
  gsl_vector *bound_box = NULL;

  /* Check input */
  if ((int)metric->size1 != n || (int)metric->size2 != n)
    XLAL_ERROR_NULL("'metric' is not square", XLAL_ESIZE);

  /* Allocate memory */
  ALLOC_GSL_MATRIX(LU_decomp, n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_PERMUTATION(LU_perm, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(inverse, n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_VECTOR(bound_box, n, XLAL_ERROR_NULL);

  /* Compute metric inverse */
  gsl_matrix_memcpy(LU_decomp, metric);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);
  gsl_linalg_LU_invert(LU_decomp, LU_perm, inverse);

  /* Compute bounding box */
  for (i = 0; i < n; ++i) {
    gsl_vector_set(bound_box, i, sqrt(max_mismatch * gsl_matrix_get(inverse, i ,i)));
  }

  /* Cleanup */
  FREE_GSL_MATRIX(LU_decomp);
  FREE_GSL_PERMUTATION(LU_perm);
  FREE_GSL_MATRIX(inverse);

  return bound_box;

}

/**
 * Orthonormalise the columns of a matrix with respect to a metric (matrix is lower triangular)
 */
int XLALOrthonormaliseWRTMetric(
				gsl_matrix *matrix, /**< Matrix of columns to orthonormalise */
				gsl_matrix *metric  /**< Metric to orthonormalise with respect to */
				)
{
  
  const int n = metric->size1;

  int i, j;
  gsl_vector *temp = NULL;
  double inner_prod = 0.0;

  /* Check input */
  if ((int)metric->size1 != n || (int)metric->size2 != n)
    XLAL_ERROR("'metric' is not square", XLAL_ESIZE);
  if (metric->size1 != matrix->size2 || metric->size2 != matrix->size2)
    XLAL_ERROR("'matrix' is not the same size as 'metric'", XLAL_ESIZE);

  /* Allocate */
  ALLOC_GSL_VECTOR(temp, n, XLAL_ERROR);

  /* Orthonormalise the columns of the matrix using numerically stabilised Gram-Schmidt */
  for (i = n - 1; i >= 0; --i) {
    gsl_vector_view col_i = gsl_matrix_column(matrix, i);
    for (j = n - 1; j > i; --j) {
      gsl_vector_view col_j = gsl_matrix_column(matrix, j);

      /* Compute inner product of jth and ith columns with the metric */
      gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &col_j.vector, 0.0, temp);
      gsl_blas_ddot(&col_i.vector, temp, &inner_prod);

      /* Subtract component of jth column from ith column */
      gsl_vector_memcpy(temp, &col_j.vector);
      gsl_vector_scale(temp, inner_prod);
      gsl_vector_sub(&col_i.vector, temp);
      
    }

    /* Compute inner product of ith column with itself */
    gsl_blas_dgemv(CblasNoTrans, 1.0, metric, &col_i.vector, 0.0, temp);
    gsl_blas_ddot(&col_i.vector, temp, &inner_prod);

    /* Normalise ith column */
    gsl_vector_scale(&col_i.vector, 1.0 / sqrt(inner_prod));    

  }

  /* Cleanup */
  FREE_GSL_VECTOR(temp);

  return XLAL_SUCCESS;

}

/**
 * Transform a lattice generator to a square lower triangular form
 */
gsl_matrix *XLALSquareLowerTriangularLatticeGenerator(
						      gsl_matrix *generator /**< Generator matrix of lattice */
						      )
{

  const int m = generator->size1;
  const int n = generator->size2;

  int i, j;
  gsl_matrix *QR_decomp = NULL;
  gsl_vector *QR_tau = NULL;
  gsl_matrix *Q = NULL;
  gsl_matrix *R = NULL;
  gsl_matrix *perm_sign = NULL;
  gsl_matrix *left = NULL;
  gsl_matrix *right = NULL;
  gsl_matrix *temp = NULL;
  gsl_matrix *result = NULL;
  double x = 0.0;

  /* Check input */
  if (m < n)
    XLAL_ERROR_NULL("'generator' must have number of rows >= number of columns", XLAL_ESIZE);
  
  /* Allocate memory */
  ALLOC_GSL_MATRIX(QR_decomp, m, n, XLAL_ERROR_NULL);
  ALLOC_GSL_VECTOR(QR_tau,       n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(Q,         m, m, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(R,         m, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(perm_sign, n, m, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(left,      n, m, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(right,     n, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(temp,      m, n, XLAL_ERROR_NULL);
  ALLOC_GSL_MATRIX(result,    n, n, XLAL_ERROR_NULL);

  /* Find the QR decomposition of the generator */
  gsl_matrix_memcpy(QR_decomp, generator);
  gsl_linalg_QR_decomp(QR_decomp, QR_tau);
  gsl_linalg_QR_unpack(QR_decomp, QR_tau, Q, R);

  /* Build matrix to permute column order and make signs to diagonal positive */
  gsl_matrix_set_zero(perm_sign);
  for (i = 0; i < n; ++i) {
    for (j = 0; j < m; ++j) {
      if (i + j == n - 1) {
	x = gsl_matrix_get(R, j, j);
	gsl_matrix_set(perm_sign, i, j, x < 0 ? -1.0 : (x > 0 ? 1.0 : 0.0));
      }
    }
  }

  /* Calculate left side of transform (Q is transposed to get inverse) */
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, perm_sign, Q, 0.0, left);

  /* Build right side of transform */
  gsl_matrix_set_zero(right);
  for (i = 0; i < n; ++i)
    for (j = 0; j < n; ++j)
      if (i + j == n - 1)
	gsl_matrix_set(right, i, j, 1.0);

  /* Transform generator */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, generator, right, 0.0, temp);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, left, temp, 0.0, result);

  /* Generator will be lower triangular, so zero out upper triangle */
  for (i = 0; i < n; ++i)
    for (j = i + 1; j < n; ++j)
      gsl_matrix_set(result, i, j, 0.0);

  /* Cleanup */
  FREE_GSL_MATRIX(QR_decomp);
  FREE_GSL_VECTOR(QR_tau);
  FREE_GSL_MATRIX(Q);
  FREE_GSL_MATRIX(R);
  FREE_GSL_MATRIX(perm_sign);
  FREE_GSL_MATRIX(left);
  FREE_GSL_MATRIX(right);
  FREE_GSL_MATRIX(temp);
  
  return result;

}

/**
 * Normalise a lattice generator matrix to have a specified covering radius
 */
int XLALNormaliseLatticeGenerator(
				  gsl_matrix *generator, /**< Generator matrix of lattice */
				  REAL8 norm_thickness,  /**< Normalised thickness of lattice */
				  REAL8 covering_radius  /**< Desired covering radius */
				  )
{
  
  const int m = generator->size1;
  const int n = generator->size2;

  gsl_matrix *LU_decomp = NULL;
  gsl_permutation *LU_perm = NULL;
  int LU_sign = 0;
  double generator_determinant = 0.0;
  double generator_covering_radius = 0.0;

  /* Allocate memory */
  ALLOC_GSL_MATRIX(LU_decomp, m, n, XLAL_ERROR);
  ALLOC_GSL_PERMUTATION(LU_perm, m, XLAL_ERROR);

  /* Compute generator LU decomposition */
  gsl_matrix_memcpy(LU_decomp, generator);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);

  /* Compute generator determinant */
  generator_determinant = gsl_linalg_LU_det(LU_decomp, LU_sign);

  /* Compute covering radius */
  generator_covering_radius = pow(norm_thickness * generator_determinant, 1.0 / n);

  /* Normalise so covering spheres have specified covering radius */
  gsl_matrix_scale(generator, covering_radius / generator_covering_radius);

  /* Cleanup */
  FREE_GSL_MATRIX(LU_decomp);
  FREE_GSL_PERMUTATION(LU_perm);

  return XLAL_SUCCESS;
  
}

/**
 * Set the tiling to a cubic lattice generator
 */
int XLALSetFlatTilingCubicLattice(
				  FlatLatticeTiling *tiling /**< Tiling structure */
				  )
{
  
  const int r = tiling->reduced_dims;
  
  gsl_matrix *generator = NULL;
  REAL8 norm_thickness = 0.0;

  /* Check for non-singular dimensions */
  if (r > 0) {

    /* Allocate memory */
    ALLOC_GSL_MATRIX(generator, r, r, XLAL_ERROR);
    
    /* Create generator */
    gsl_matrix_set_identity(generator);
    
    /* Calculate normalised thickness */
    norm_thickness = pow(sqrt(r)/2, r);

  }

  /* Set lattice generator */
  if (XLALSetFlatTilingLattice(tiling, generator, norm_thickness) != XLAL_SUCCESS)
    XLAL_ERROR("SetTilingLattice failed", XLAL_EFAILED);

  /* Cleanup */
  FREE_GSL_MATRIX(generator);

  return XLAL_SUCCESS;

}

/**
 * Set the tiling to an \fA_n^*\f lattice generator
 */
int XLALSetFlatTilingAnstarLattice(
				   FlatLatticeTiling *tiling /**< Tiling structure */
				   )
{
  
  const int r = tiling->reduced_dims;

  gsl_matrix *generator = NULL;
  REAL8 norm_thickness = 0.0;

  /* Check for non-singular dimensions */
  if (r > 0) {

    /* Allocate memory */
    ALLOC_GSL_MATRIX(generator, r + 1, r, XLAL_ERROR);
    
    /* Create generator in (r + 1) space */
    gsl_matrix_set_all(generator, 0.0);
    {
      gsl_vector_view first_row = gsl_matrix_row(generator, 0);
      gsl_vector_view sub_diag = gsl_matrix_subdiagonal(generator, 1);
      gsl_vector_view last_col = gsl_matrix_column(generator, r - 1);
      gsl_vector_set_all(&first_row.vector, 1.0);
      gsl_vector_set_all(&sub_diag.vector, -1.0);
      gsl_vector_set_all(&last_col.vector, 1.0 / (r + 1.0));
      gsl_vector_set(&last_col.vector, 0, -r * 1.0 / (r + 1.0));
    }
    
    /* Calculate normalised thickness */
    norm_thickness = sqrt(r + 1.0)*pow((1.0*r*(r + 2))/(12.0*(r + 1)), 0.5*r);

  }

  /* Set lattice generator */
  if (XLALSetFlatTilingLattice(tiling, generator, norm_thickness) != XLAL_SUCCESS)
    XLAL_ERROR("SetTilingLattice failed", XLAL_EFAILED);

  /* Cleanup */
  FREE_GSL_MATRIX(generator);

  return XLAL_SUCCESS;

}

/**
 * Create a flat lattice tiling with a simple square parameter space
 */
int XLALFlatLatticeTilingSquareParamSpace(
					  FlatLatticeTiling **tiling, /**< Tiling structure */
					  gsl_vector *square_bounds   /**< Bounds on each dimension: start,width,start,width,... */
					  )
{

  /* Create tiling if it hasn't been created */
  if (!*tiling) {

    /* Check there is an even number of bounds */
    if (!square_bounds || GSL_IS_ODD(square_bounds->size))
      XLAL_ERROR("'square_bounds' must have an even number of elements", XLAL_EFAILED);
    
    /* Create a flat lattice tiling */
    if ((*tiling = XLALCreateFlatLatticeTiling(square_bounds->size / 2)) == NULL)
      XLAL_ERROR("XLALCreateFlatLatticeTiling failed", XLAL_EFAILED);

  }

  /* Otherwise create square parameter space */
  else {

    const int n = (*tiling)->dimensions;

    int i;
    gsl_vector *lower_const = NULL;
    gsl_vector *upper_const = NULL;
    
    /* Allocate memory */
    ALLOC_GSL_VECTOR(lower_const, 1, XLAL_ERROR);
    ALLOC_GSL_VECTOR(upper_const, 1, XLAL_ERROR);

    /* Iterate over start/width pairs */
    for (i = 0; i < n; ++i) {

      /* Get start/width pair */
      const double start = gsl_vector_get(square_bounds, 2*i    );
      const double width = gsl_vector_get(square_bounds, 2*i + 1);

      /* If width is zero, set singular bound */
      if (width == 0.0) {
	if (XLAL_SUCCESS != XLALAddFlatLatticeTilingSingularBound(*tiling, i, 0, start))
	  XLAL_ERROR("XLALAddFlatLatticeTilingSingularBound failed", XLAL_EFAILED);
      }

      /* Otherwise set polynomial bound with only constant term */
      else {
	gsl_vector_set(lower_const, 0, start);
	gsl_vector_set(upper_const, 0, start + width);
	if (XLAL_SUCCESS != XLALAddFlatLatticeTilingPolynomialBound(*tiling, i, 0, lower_const, NULL, upper_const, NULL))
	  XLAL_ERROR("XLALAddFlatLatticeTilingPolynomialBound failed", XLAL_EFAILED);
      }

    }

    /* Finalise bounds */
    if (XLAL_SUCCESS != XLALFinaliseFlatLatticeTilingBounds(*tiling))
      XLAL_ERROR("XLALFinaliseFlatLatticeTilingBounds failed", XLAL_EFAILED);
    
    /* Cleanup */
    FREE_GSL_VECTOR(lower_const);
    FREE_GSL_VECTOR(upper_const);
    
  }
  
  return XLAL_SUCCESS;

}	     
