/*********************************************************************************
 *  \author K. Wette
 *  \file
 *  \ingroup templateBanks
 *  \brief
 *  Flat lattice tiling over multi-dimensioned parameter spaces
 *********************************************************************************/

/******** Includes ********/

#include "FlatLatticeTiling.h"

#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_nan.h>

#include <lal/LALRCSID.h>
#include <lal/LALError.h>
#include <lal/XLALError.h>
#include <lal/LALMalloc.h>
#include <lal/LALConstants.h>
#include <lal/LALError.h>
#include <lal/LALStdlib.h>
#include <lal/LogPrintf.h>

/******** Constants ********/

static const int TRUE  = (1==1);
static const int FALSE = (1==0);

/******** Declarations ********/

NRCSID(FLATLATTICETILINGC, "$Id$");

/******** Functions ********/

/*
 *  Create a FlatLatticeTiling structure to hold information
 *  for the flat lattice tiling routines.
 */
FlatLatticeTiling *XLALCreateFlatLatticeTiling(
					       UINT4 dimension, /**< [in] Dimension of the parameter space */
					       REAL8 mismatch   /**< [in] Maximum mismatch of the parameter space templates */
					       )
{

  FlatLatticeTiling *tiling = NULL;

  if (dimension < 1) {
    LALPrintError("%s\nERROR: Tiling dimension must be non-zero\n", FLATLATTICETILINGC);
    XLAL_ERROR_NULL("XLALCreateFlatLatticeTiling", XLAL_EINVAL);
  }

  if (mismatch <= 0.0) {
    LALPrintError("%s\nERROR: Tiling mismatch must be non-zero\n", FLATLATTICETILINGC);
    XLAL_ERROR_NULL("XLALCreateFlatLatticeTiling", XLAL_EINVAL);
  }

  if ((tiling = (FlatLatticeTiling*) LALMalloc(sizeof(FlatLatticeTiling))) == NULL) {
    LALPrintError("%s\nERROR: Could not allocate a FlatLatticeTiling\n", FLATLATTICETILINGC);
    XLAL_ERROR_NULL("XLALCreateFlatLatticeTiling", XLAL_ENOMEM);
  }

  tiling->dimension = dimension;

  tiling->is_degenerate = gsl_vector_int_alloc(dimension);
  gsl_vector_int_set_all(tiling->is_degenerate, FALSE);
  tiling->first_degenerate = -1;

  tiling->metric = gsl_matrix_calloc(dimension, dimension);
  tiling->metric_congruent_factor = gsl_matrix_alloc(dimension, dimension);
  gsl_matrix_set_identity(tiling->metric_congruent_factor);

  tiling->mismatch = mismatch;

  tiling->generator = gsl_matrix_calloc(dimension, dimension);
  tiling->generator_norm_thickness = 0.0;

  tiling->num_bounds = 0;

  tiling->param_bound_normal = NULL;
  tiling->param_bound_origin = NULL;

  tiling->latt_bound_normal = NULL;
  tiling->latt_bound_dot = NULL;

  tiling->param_lower = gsl_vector_calloc(dimension);
  tiling->param_upper = gsl_vector_calloc(dimension);

  tiling->latt_to_param = gsl_matrix_alloc(dimension, dimension);

  tiling->latt_lower = gsl_vector_long_calloc(dimension);
  tiling->latt_upper = gsl_vector_long_calloc(dimension);

  tiling->line_index = 0;

  tiling->iter_order = gsl_vector_int_alloc(dimension);

  tiling->latt_current = NULL;
  tiling->current = NULL;

  tiling->line_latt_lower = 0;
  tiling->line_latt_upper = 0;

  tiling->line_start = gsl_vector_alloc(dimension);
  tiling->line_dir = gsl_vector_alloc(dimension);

  tiling->line_start_dot = NULL;
  tiling->line_dir_dot = NULL;

  tiling->line_intersect = NULL;

  tiling->in_param_space = NULL;
  tiling->in_param_space_args = NULL;

  return tiling;

}

/*
 *  Destroy a FlatLatticeTiling structure
 */
void XLALDestroyFlatLatticeTiling(
				  FlatLatticeTiling *tiling /**< [in] Structure */
				  )
{

  if (tiling) {

    if (tiling->is_degenerate) gsl_vector_int_free(tiling->is_degenerate);

    if (tiling->metric) gsl_matrix_free(tiling->metric);
    if (tiling->metric_congruent_factor) gsl_matrix_free(tiling->metric_congruent_factor);

    if (tiling->generator) gsl_matrix_free(tiling->generator);

    if (tiling->param_bound_normal) gsl_matrix_free(tiling->param_bound_normal);
    if (tiling->param_bound_origin) gsl_matrix_free(tiling->param_bound_origin);

    if (tiling->latt_bound_normal) gsl_matrix_free(tiling->latt_bound_normal);
    if (tiling->latt_bound_dot) gsl_vector_free(tiling->latt_bound_dot);

    if (tiling->param_lower) gsl_vector_free(tiling->param_lower);
    if (tiling->param_upper) gsl_vector_free(tiling->param_upper);

    if (tiling->latt_to_param) gsl_matrix_free(tiling->latt_to_param);

    if (tiling->latt_lower) gsl_vector_long_free(tiling->latt_lower);
    if (tiling->latt_upper) gsl_vector_long_free(tiling->latt_upper);

    if (tiling->iter_order) gsl_vector_int_free(tiling->iter_order);

    if (tiling->latt_current) gsl_vector_long_free(tiling->latt_current);
    if (tiling->current) gsl_vector_free(tiling->current);

    if (tiling->line_start) gsl_vector_free(tiling->line_start);
    if (tiling->line_dir) gsl_vector_free(tiling->line_dir);

    if (tiling->line_start_dot) gsl_vector_free(tiling->line_start_dot);
    if (tiling->line_dir_dot) gsl_vector_free(tiling->line_dir_dot);

    if (tiling->line_intersect) gsl_vector_free(tiling->line_intersect);

    if (tiling->in_param_space_args) LALFree(tiling->in_param_space_args);

    LALFree(tiling);
    tiling = NULL;

  }

}

/*
 *  Add a parameter space bound to a FlatLatticeTiling structure
 */
int XLALAddParameterSpaceBound(
			       FlatLatticeTiling* tiling, /**< [in] Structure */
			       UINT4 num_bounds,          /**< [in] Number of bounds */
			       UINT4 bound_index,         /**< [in] Index of this bound */
			       gsl_vector *normal,        /**< [in] Normal vector of bound */
			       gsl_vector *origin         /**< [in] Origin of bound */
			       )
{

  UINT4 i;

  if (num_bounds < 1) {
    LALPrintError("%s\nERROR: Number of bounds must be non-zero\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALAddParameterSpaceBound", XLAL_EINVAL);
  }
  
  if (tiling->num_bounds == 0) {

    tiling->num_bounds = num_bounds;
  
    tiling->param_bound_normal = gsl_matrix_calloc(tiling->dimension, tiling->num_bounds);
    tiling->param_bound_origin = gsl_matrix_calloc(tiling->dimension, tiling->num_bounds);

    tiling->latt_bound_normal = gsl_matrix_calloc(tiling->dimension, tiling->num_bounds);
    tiling->latt_bound_dot = gsl_vector_calloc(tiling->num_bounds);

  }

  for (i = 0; i < tiling->dimension; ++i) {
    gsl_matrix_set(tiling->param_bound_normal, i, bound_index, gsl_vector_get(normal, i));
    gsl_matrix_set(tiling->param_bound_origin, i, bound_index, gsl_vector_get(origin, i));
  }

  return XLAL_SUCCESS;

}

/*
 *  Creates a square parameter space for the flat lattice tiling
 */
int XLALSquareParameterSpace(
			     FlatLatticeTiling *tiling, /**< [in] Structure */
			     gsl_vector *lower,         /**< [in] Lower corner of space */
			     gsl_vector *upper          /**< [in] Upper corner of space */
			     )
{

  UINT4 i;

  gsl_vector *normal;
  gsl_vector *origin;

  normal = gsl_vector_alloc(tiling->dimension);
  origin = gsl_vector_alloc(tiling->dimension);

  /* Check dimensions agree */
  if (lower->size != upper->size || upper->size != tiling->dimension) {
    LALPrintError("%s\nERROR: Lengths of lower and upper must be equal to dimension\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALSquareParameterSpace", XLAL_EINVAL);
  }

  /* Copy lower and upper bounds */
  gsl_vector_memcpy(tiling->param_lower, lower);
  gsl_vector_memcpy(tiling->param_upper, upper);

  /* Iterate over each dimension */
  for (i = 0; i < tiling->dimension; ++i) {

    /* Zero vectors */
    gsl_vector_set_zero(normal);
    gsl_vector_set_zero(origin);

    /* Create lower bound */
    gsl_vector_set(normal, i, -1.0);
    gsl_vector_set(origin, i, gsl_vector_get(lower, i));
    if (XLALAddParameterSpaceBound(tiling, 2*tiling->dimension, 2*i, normal, origin) != XLAL_SUCCESS) {
      LALPrintError("%s\nERROR: XLALAddParameterSpaceBound failed\n", FLATLATTICETILINGC);
      XLAL_ERROR("XLALSquareParameterSpace", XLAL_EFUNC);
    }

    /* Create upper bound */
    gsl_vector_set(normal, i, 1.0);
    gsl_vector_set(origin, i, gsl_vector_get(upper, i));
    if (XLALAddParameterSpaceBound(tiling, 2*tiling->dimension, 2*i + 1, normal, origin) != XLAL_SUCCESS) {
      LALPrintError("%s\nERROR: XLALAddParameterSpaceBound failed\n", FLATLATTICETILINGC);
      XLAL_ERROR("XLALSquareParameterSpace", XLAL_EFUNC);
    }

  }

  /* Cleanup */
  gsl_vector_free(normal);
  gsl_vector_free(origin);

  return XLAL_SUCCESS;

}

/*
 *  Returns the generator matrix for a $Z_n$ or cubic lattice,
 *  normalised so that the covering radius is unity.
 */
int XLALCubicLatticeGenerator(
			      gsl_matrix* generator, /**< [in/out] Generator matrix of lattice */
			      REAL8 *norm_thickness  /**< [out] Normalised thickness of lattice */
			      )
{

  double n = 0.0;

  /* Check matrix is square */
  if (generator->size1 != generator->size2) {
    LALPrintError("%s\nERROR: Generator must be square\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALCubicLatticeGenerator", XLAL_EINVAL);
  }
  n = generator->size1;

  /* Generator */
  gsl_matrix_set_identity(generator);

  /* Normalised thickness */
  *norm_thickness = pow(sqrt(n)/2, n);

  return XLAL_SUCCESS;

}

/*
 *  Returns the generator matrix for a $A_n^*$ lattice, which are the optimal known
 *  covering lattice for $n <= 23$. For example, $A_3^*$ is the body-centered cubic
 *  lattice. See \ref CS99, p115.
 */
int XLALAnstarLatticeGenerator(
			       gsl_matrix* generator, /**< [in/out] Generator matrix of lattice */
			       REAL8 *norm_thickness  /**< [out] Normalised thickness of lattice */
			       )
{

  UINT4 i;
  double n = 0.0;

  /* Check matrix is square */
  if (generator->size1 != generator->size2) {
    LALPrintError("%s\nERROR: Generator must be square\n", FLATLATTICETILINGC);
    XLAL_ERROR("XLALAnStarLatticeGenerator", XLAL_EINVAL);
  }
  n = generator->size1;

  /* Generator */
  gsl_matrix_set_identity(generator);
  for (i = 0; i < generator->size1; ++i) {
    gsl_matrix_set(generator, i, 0, 0.5);
  }
  if (generator->size1 == 2) {
    /* Special case for planar hexagonal */
    gsl_matrix_set(generator, 0, 0, sqrt(3.0) / 2.0);
  }

  /* Normalised thickness */
  *norm_thickness = sqrt(n + 1)*pow((n*(n + 2))/(12*(n + 1)), n/2);

  return XLAL_SUCCESS;

}

/*
 *  Complete initialisation of a FlatLatticeTiling structure,
 *  setting it up for use by XLALNextFlatLatticeTile.
 */
int XLALSetupFlatLatticeTiling(
			       FlatLatticeTiling *tiling /**< [in/out] Structure */
			       )
{

  UINT4 i, j;

  UINT4 n = tiling->dimension;

  gsl_eigen_symmv_workspace *eigen_wksp = NULL;
  gsl_permutation *LU_perm = NULL;
  gsl_matrix *LU_decomp = NULL;
  gsl_matrix *eigen_vec = NULL;
  gsl_matrix *real_metric = NULL;
  gsl_matrix *temp_1 = NULL;
  gsl_matrix *temp_2 = NULL;
  gsl_vector *eigen_val = NULL;
  gsl_vector *latt_lower = NULL;
  gsl_vector *latt_upper = NULL;
  gsl_vector *bound_normal = NULL;
  gsl_vector *bound_origin = NULL;
  gsl_vector *temp_3 = NULL;
  gsl_vector_int *index = NULL;

  double eigen_rhs = 0.0;
  double mismatch_factor = 0.0;
  double eigen_val_j = 0.0;
  double generator_determinant = 0.0;
  double covering_radius = 0.0;
  double lower = 0.0;
  double upper = 0.0;
  double bound_dot = 0.0;
  double bound_norm = 0.0;
  long range = 0;
  long max_range = 0;
  int LU_sign = 0;

  eigen_wksp = gsl_eigen_symmv_alloc(n);
  LU_perm =    gsl_permutation_alloc(n);
  LU_decomp =       gsl_matrix_alloc(n, n);
  eigen_vec =       gsl_matrix_alloc(n, n);
  real_metric =     gsl_matrix_alloc(n, n);
  temp_1 =          gsl_matrix_alloc(n, n);
  temp_2 =          gsl_matrix_alloc(n, n);
  eigen_val =       gsl_vector_alloc(n);
  latt_lower =      gsl_vector_alloc(n);
  latt_upper =      gsl_vector_alloc(n);
  bound_normal =    gsl_vector_alloc(n);
  bound_origin =    gsl_vector_alloc(n);
  temp_3 =          gsl_vector_alloc(n);
  index =           gsl_vector_int_alloc(n);
  
  /* Print debugging */
  LogPrintf(LOG_DETAIL, "Setting up flat lattice tiling...\n");

  /* Check that the parameter space bounds are the right way round */
  for (i = 0; i < n; ++i) {
    if (gsl_vector_get(tiling->param_lower, i) > gsl_vector_get(tiling->param_upper, i)) {
      LALPrintError("%s\nERROR: lower parameter space bound must be less than upper bound\n", FLATLATTICETILINGC);
      XLAL_ERROR("XLALSetupFlatLatticeTiling", XLAL_EINVAL);
    }
  }

  /*
   *  Compute the eigenvectors and values of the metric
   */

  gsl_matrix_memcpy(temp_1, tiling->metric);
  gsl_eigen_symmv(temp_1, eigen_val, eigen_vec, eigen_wksp);
  gsl_eigen_symmv_sort(eigen_val, eigen_vec, GSL_EIGEN_SORT_ABS_ASC);

  /* Print debugging */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tiling->metric, eigen_vec, 0.0, temp_2);
  LogPrintfVerbatim(LOG_DETAIL, "Error in calculation of metric eigenvectors/values:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      eigen_rhs = gsl_matrix_get(eigen_vec, i, j) * gsl_vector_get(eigen_val, j);
      LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", (gsl_matrix_get(temp_2, i, j) - eigen_rhs) /  eigen_rhs);
    }
    LogPrintfVerbatim(LOG_DETAIL, ";\n");
  }
  
  /* Calculate quadratic product and print debugging */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tiling->metric, eigen_vec, 0.0, temp_1);
  gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, eigen_vec, temp_1, 0.0, temp_2);
  for (j = 0; j < n; ++j) {
    gsl_matrix_set(temp_2, j, j, gsl_matrix_get(temp_2, j, j) - mismatch_factor);
  }
  LogPrintfVerbatim(LOG_DETAIL, "Error in quadratic product of metric with re-scaled eigenvectors:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", fabs(gsl_matrix_get(temp_2, i, j)));
    }
    LogPrintfVerbatim(LOG_DETAIL, ";\n");
  }
  
  /*
   *  Divide the generator matrix through by the covering radius (worked out from the
   *  normalised thickness and the determinant) so that the covering spheres are unity.
   */

  /* Compute generator LU decomposition */
  gsl_matrix_memcpy(LU_decomp, tiling->generator);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);

  /* Compute generator determinant */
  generator_determinant = gsl_linalg_LU_det(LU_decomp, LU_sign);

  /* Compute covering radius */
  covering_radius = pow(tiling->generator_norm_thickness * generator_determinant, 1.0 / n);

  /* Scale generator */
  gsl_matrix_scale(tiling->generator, 1.0 / covering_radius);

  /* Print debugging */
  LogPrintfVerbatim(LOG_DETAIL, "Generator matrix:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", gsl_matrix_get(tiling->generator, i, j));
    }
    LogPrintfVerbatim(LOG_DETAIL, ";\n");
  }

  /*  WRONG!!!
   *  Re-scale each eigenvector $v$ by $v' = v \sqrt{\mu / \lambda}$, where
   *  $\mu$ is the mismatch, $\lambda$ is the eigenvalue, so that:
   *  $g_{ij} v'^i v'^j = V'^{tr} G V'                        $
   *  $                 = \mu / \lambda * (V^{tr} G V)$
   *  $                 = \mu / \lambda * \lambda     $
   *  $                 = \mu                         $
   */

  mismatch_factor = tiling->mismatch;
  for (j = 0; j < n; ++j) {
    eigen_val_j = gsl_vector_get(eigen_val, j);
    for (i = 0; i < n; ++i) {
      gsl_matrix_set(eigen_vec, i, j, gsl_matrix_get(eigen_vec, i, j) *
		     (eigen_val_j > 0 ? 1 : -1) *
		     sqrt(mismatch_factor / fabs(eigen_val_j)));
    }
  }

  /*
   *  Calculate the transformation from the lattice to the parameter space:
   *  $p = T \eta = P^{-tr} V G \eta$, where
   *  $\eta \el Z_n$ is a vector of integers describing the position of a lattice point
   *  relative to the origin (e.g. how many points across? up? down?), $G$ is the generator
   *  matrix of the lattice, which transforms $\eta$ into the lattice points with unit 
   *  covering radius, $V$ are the rescaled eigenvectors, which scale and rotate $G \eta$
   *  appropriately, and $P$ is the congruent factor of the real metric $M$ (such that
   *  real_metric = $P^{tr}$ tiling->metric $P$) which converts $V$ into the eigenvectors
   *  of the real matrix $V'$, such that $V'^{tr} M V'$ is equal to the appropriate maximum
   *  mismatch factors.
   */

  /* Compute congruent factor LU decomposition */
  gsl_matrix_memcpy(LU_decomp, tiling->metric_congruent_factor);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);

  /* Compute congruent factor inverse */
  gsl_linalg_LU_invert(LU_decomp, LU_perm, temp_1);

  /* Compute part of transformation matrix */
  gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, temp_1, eigen_vec, 0.0, tiling->latt_to_param);

  /* Calculate real metric and print debugging */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, tiling->metric, tiling->metric_congruent_factor, 0.0, temp_1);
  gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, tiling->metric_congruent_factor, temp_1, 0.0, real_metric);
  LogPrintfVerbatim(LOG_DETAIL, "Real parameter space metric:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", gsl_matrix_get(real_metric, i, j));
    }
    LogPrintfVerbatim(LOG_DETAIL, ";\n");
  }
  
  /* Calculate quadratic product and print debugging */
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, real_metric, tiling->latt_to_param, 0.0, temp_1);
  gsl_blas_dgemm(CblasTrans,   CblasNoTrans, 1.0, tiling->latt_to_param, temp_1, 0.0, temp_2);
  for (j = 0; j < n; ++j) {
    gsl_matrix_set(temp_2, j, j, gsl_matrix_get(temp_2, j, j) - mismatch_factor);
  }
  LogPrintfVerbatim(LOG_DETAIL, "Error in quadratic product of real metric with transformation matrix (excluding generator):\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", fabs(gsl_matrix_get(temp_2, i, j)));
    }
    LogPrintfVerbatim(LOG_DETAIL, ";\n");
  }
  
  /* Compute rest of transformation matrix */
  gsl_matrix_memcpy(temp_1, tiling->latt_to_param);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, temp_1, tiling->generator, 0.0, tiling->latt_to_param);

  /* Print debugging */
  LogPrintfVerbatim(LOG_DETAIL, "Transformation matrix:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", gsl_matrix_get(tiling->latt_to_param, i, j));
    }
    LogPrintfVerbatim(LOG_DETAIL, ";\n");
  }

  /*
   *  Check that the size of the templating step (diagonal element of the transformation matrix)
   *  are not larger than the size of the parameter space in that dimension. If it is, a full space
   *  of templates will not be generated for that dimension - in fact, if it's small enough, no templates
   *  will be generated, even though there are still templates in other dimensions. So, we mark the 
   *  dimension as degenerate. Also, find the store the first such dimension.
   */

  tiling->first_degenerate = -1;
  for (i = 0; i < n; ++i) {
    lower = gsl_vector_get(tiling->param_lower, i);
    upper = gsl_vector_get(tiling->param_upper, i);
    if (fabs((upper - lower) / gsl_matrix_get(tiling->latt_to_param, i, i)) < 1.0) {
      gsl_vector_int_set(tiling->is_degenerate, i, TRUE);
      if (tiling->first_degenerate < 0) {
	tiling->first_degenerate = i;
      }
    }
  }
  
  /* Print debugging */
  LogPrintfVerbatim(LOG_DETAIL, "Degenerate dimensions:\n");
  for (i = 0; i < n; ++i) {
    LogPrintfVerbatim(LOG_DETAIL, "%s ", (gsl_vector_int_get(tiling->is_degenerate, i) ? "true" : "false"));
  }
  LogPrintfVerbatim(LOG_DETAIL, "\n");

  /*
   *  Find the lower and upper bounds on the parameter space in lattice coordinates.
   *  The points are first translated by the lower bound of the parameter space, so:
   *  $\eta = T^{-1} (p - p_{lb})$
   *  $p    = T    \eta + p_{lb} $
   *  The bounds are found from all possible combination of the bounds in parameter
   *  space coordinates, i.e. the $2^n$ corners of the bounding cube in $R^n$. The
   *  bounds therefore form a maximal bounding box around the parameter space.
   */

  /* Compute LU decomposition of transformation matrix */
  gsl_matrix_memcpy(LU_decomp, tiling->latt_to_param);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);

  /* Initialise bounds and index */
  for (i = 0; i < n; ++i) {
    gsl_vector_set(latt_lower, i, GSL_POSINF);
    gsl_vector_set(latt_upper, i, GSL_NEGINF);
    gsl_vector_int_set(index, i, 0);
  }

  /* Iterate through all corners of the parameter space */
  while (gsl_vector_int_get(index, n-1) <= 1) {

    /* Construct lattice point */
    for (i = 0; i < n; ++i) {
      gsl_vector_set(temp_3, i, gsl_vector_int_get(index, i) > 0 ? 
		     gsl_vector_get(tiling->param_upper, i) - gsl_vector_get(tiling->param_lower, i) : 0.0);
    }

    /* Transform from parameter to lattice space */
    gsl_linalg_LU_svx(LU_decomp, LU_perm, temp_3);

    /* Update bounds */
    for (i = 0; i < n; ++i) {
      gsl_vector_set(latt_lower, i, GSL_MIN_DBL(gsl_vector_get(latt_lower, i), gsl_vector_get(temp_3, i)));
      gsl_vector_set(latt_upper, i, GSL_MAX_DBL(gsl_vector_get(latt_upper, i), gsl_vector_get(temp_3, i)));
    }

    /* Find the next corner */
    for (i = 0; i < n; ++i) {
      gsl_vector_int_set(index, i, gsl_vector_int_get(index, i) + 1);
      if (gsl_vector_int_get(index, i) <= 1) {
	break;
      }
      else if (i < n-1) {
	gsl_vector_int_set(index, i, 0);
      }
    }

  }

  /* Convert the (double) lattice space bounds to the (long) bounds */
  for (i = 0; i < n; ++i) {
    gsl_vector_long_set(tiling->latt_lower, i, floor(gsl_vector_get(latt_lower, i)));
    gsl_vector_long_set(tiling->latt_upper, i,  ceil(gsl_vector_get(latt_upper, i)));
  }

  /* Print debugging */
  LogPrintfVerbatim(LOG_DETAIL, "Lower -> upper bound on parameter space in parameter (lattice) coordinates:\n");
  for (i = 0; i < n; ++i) {
    LogPrintfVerbatim(LOG_DETAIL, "% 0.8e -> % 0.8e (% ld -> % ld)\n",
		      gsl_vector_get(tiling->param_lower, i), gsl_vector_get(tiling->param_upper, i),
		      gsl_vector_long_get(tiling->latt_lower, i), gsl_vector_long_get(tiling->latt_upper, i));
  }

  /*
   *  Convert the normal and origin of the parameter space bounds from parameter to lattice space.
   *  The normal $n$ and origin $x_0$ describe a hyperspace containing points $x$ with the relation
   *  $n . (x - x_0) = n^{T} (x - x_0) = 0$ or $n . x = c$
   *  where $n, x, x_0$ are in real parameter coordinates. If we transform to lattice coordinates with
   *  $x = T x' + p_{lb}, x_0 = T x_0' + p_{lb}, n = T^{-T} n'$
   *  where $n', x', x_0'$ are in lattice coordinates, then the inner product is preserved:
   *  $n . (x - x_0) = n^{T} (x - x_0)          $
   *  $              = n' T^{-1} (T x' - T x_0')$
   *  $              = n' . (x' - x_0')         $
   *  The lattice coordinates are therefore:
   *  $x_0' = T^{-1} (x_0 - p_{lb}), n' = T^{T} n$
   */

  /* Compute LU decomposition of transformation matrix */
  gsl_matrix_memcpy(LU_decomp, tiling->latt_to_param);
  gsl_linalg_LU_decomp(LU_decomp, LU_perm, &LU_sign);

  /* Iterate over the normal and origin vectors */
  for (j = 0; j < tiling->num_bounds; ++j) {

    /* Copy normal and origin vectors */
    for (i = 0; i < n; ++i) {
      gsl_vector_set(bound_normal, i, gsl_matrix_get(tiling->param_bound_normal, i, j));
      gsl_vector_set(bound_origin, i, gsl_matrix_get(tiling->param_bound_origin, i, j) - gsl_vector_get(tiling->param_lower, i));
    }

    /* Transform normal vector */
    gsl_vector_memcpy(temp_3, bound_normal);
    gsl_blas_dgemv(CblasTrans, 1.0, tiling->latt_to_param, temp_3, 0.0, bound_normal);

    /* Transform origin vector */
    gsl_linalg_LU_svx(LU_decomp, LU_perm, bound_origin);

    /* Normalise normal vector */
    bound_norm = gsl_blas_dnrm2(bound_normal);
    gsl_vector_scale(bound_normal, 1.0 / bound_norm);

    /* Calculate dot product */
    gsl_blas_ddot(bound_normal, bound_origin, &bound_dot);

    /* Store transformed normal vector and dot product */
    for (i = 0; i < n; ++i) {
      gsl_matrix_set(tiling->latt_bound_normal, i, j, gsl_vector_get(bound_normal, i));
    }
    gsl_vector_set(tiling->latt_bound_dot, j, bound_dot);

  }

  /* Print debugging */
  LogPrintfVerbatim(LOG_DETAIL, "Parameter space bound normal vectors in real parameter coordinates:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < tiling->num_bounds; ++j) {
      LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", gsl_matrix_get(tiling->param_bound_normal, i, j));
    }
    LogPrintfVerbatim(LOG_DETAIL, ";\n");
  }
  LogPrintfVerbatim(LOG_DETAIL, "Parameter space bound origin vectors in real parameter coordinates:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < tiling->num_bounds; ++j) {
      LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", gsl_matrix_get(tiling->param_bound_origin, i, j));
    }
    LogPrintfVerbatim(LOG_DETAIL, ";\n");
  }
  LogPrintfVerbatim(LOG_DETAIL, "Parameter space bound normalised normal vectors in lattice coordinates:\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < tiling->num_bounds; ++j) {
      LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", gsl_matrix_get(tiling->latt_bound_normal, i, j));
    }
    LogPrintfVerbatim(LOG_DETAIL, ";\n");
  }
  LogPrintfVerbatim(LOG_DETAIL, "Parameter space bound constant in lattice coordinates:\n");
  for (j = 0; j < tiling->num_bounds; ++j) {
    LogPrintfVerbatim(LOG_DETAIL, "% 0.8e ", gsl_vector_get(tiling->latt_bound_dot, j));
  }
  LogPrintfVerbatim(LOG_DETAIL, ";\n");

  /*
   *  Find the largest dimension in lattice space coordinates. This dimension will be reduced
   *  by finding where a line along it intersects with the parameter space bounds. Save the bounds
   *  on this dimension, and find the resulting order of dimensions to be searched over, 
   *  with this dimension first.
   */

  /* Find largest dimension */
  max_range = 0;
  tiling->line_index = -1;
  for (i = 0; i < n; ++i) {
    range = gsl_vector_long_get(tiling->latt_upper, i) - gsl_vector_long_get(tiling->latt_lower, i);
    if (range > max_range) {
      tiling->line_index = i;
      max_range = range;
    }
  }

  /* Save bounds on this dimension */
  tiling->line_latt_lower = gsl_vector_long_get(tiling->latt_lower, tiling->line_index);
  tiling->line_latt_upper = gsl_vector_long_get(tiling->latt_upper, tiling->line_index);

  /* Find resulting order */
  gsl_vector_int_set(tiling->iter_order, 0, tiling->line_index);
  j = 0;
  for (i = 1; i < n; ++j) {
    if (j != tiling->line_index) {
      gsl_vector_int_set(tiling->iter_order, i++, j);
    }
  }

  /*
   *  Cleanup
   */
  gsl_eigen_symmv_free(eigen_wksp);
  gsl_permutation_free(LU_perm);
  gsl_matrix_free(LU_decomp);
  gsl_matrix_free(eigen_vec);
  gsl_matrix_free(real_metric);
  gsl_matrix_free(temp_1);
  gsl_matrix_free(temp_2);
  gsl_vector_free(eigen_val);
  gsl_vector_free(latt_lower);
  gsl_vector_free(latt_upper);
  gsl_vector_free(bound_normal);
  gsl_vector_free(bound_origin);
  gsl_vector_free(temp_3);
  gsl_vector_int_free(index);

  return XLAL_SUCCESS;

}

/*
 *  Finds the next template point in the flat lattice tiling described by the structure
 */
int XLALNextFlatLatticePoint(
			       FlatLatticeTiling *tiling /**< [in] Structure */
			       )

{

  UINT4 i, j;

  UINT4 n = tiling->dimension;

  BOOLEAN in_bounds = 0;
  long latt_curr = 0;
  BOOLEAN find_bounds = 0;
  double line_lower = 0.0;
  double line_upper = 0.0;
  double point = 0.0;
  double dir = 0.0;

  /* Keep searching until a point is found */
  in_bounds = FALSE;
  while (!in_bounds) {

    /* Assume point is with bounds already set */
    in_bounds = TRUE;
    find_bounds = FALSE;

    /* If first point, initialise point vectors */
    if (tiling->latt_current == NULL) {

      tiling->latt_current = gsl_vector_long_alloc(n);
      tiling->current = gsl_vector_alloc(n);

      tiling->line_start_dot = gsl_vector_alloc(tiling->num_bounds);
      tiling->line_dir_dot = gsl_vector_alloc(tiling->num_bounds);

      tiling->line_intersect = gsl_vector_alloc(tiling->num_bounds);

      gsl_vector_long_memcpy(tiling->latt_current, tiling->latt_lower);

      gsl_vector_set_zero(tiling->line_dir);
      gsl_vector_set(tiling->line_start, tiling->line_index, tiling->line_latt_lower);
      gsl_vector_set(tiling->line_dir, tiling->line_index, tiling->line_latt_upper - tiling->line_latt_lower);

      tiling->template_count = 0;

      find_bounds = TRUE;

    }

    /* Otherwise, advance point */
    else {
      for (i = 0; i < n; ++i) {

	/* Find current dimension from iteration order */
	j = gsl_vector_int_get(tiling->iter_order, i);

	/* Advance point along this dimension and check if still in bounds */
	latt_curr = gsl_vector_long_get(tiling->latt_current, j);
	++latt_curr;
	if (latt_curr <= gsl_vector_long_get(tiling->latt_upper, j)) {
	  gsl_vector_long_set(tiling->latt_current, j, latt_curr);
	  break;
	}
	else {

	  /* If last dimension, we are done */
	  if (i == n-1) {
	    return XLAL_FAILURE;
	  }

	  /* Otherwise reset this dimension and advance next */
	  else {
	    gsl_vector_long_set(tiling->latt_current, j, gsl_vector_long_get(tiling->latt_lower, j));
	    find_bounds = TRUE;
	  }

	}

      }
    }

    /*
     *  If requested, find the intersection of a line, intersecting the current point and parallel
     *  to the line_index'th dimension, with each of the parameter space bounds. For each intersecting
     *  point, check that it is within (in the sense of the direction of the normal vector) each of
     *  the other bounds, then consider it as a bound on the parameter space.
     */
    if (find_bounds) {

      /* Initialise line bounds */
      line_lower = GSL_NEGINF;
      line_upper = GSL_POSINF;
/*       line_lower = tiling->line_latt_upper; */
/*       line_upper = tiling->line_latt_lower; */

      /* Initialise start and direction of line */
      for (i = 0; i < n; ++i) {
	if (i != tiling->line_index) {
	  gsl_vector_set(tiling->line_start, i, gsl_vector_long_get(tiling->latt_current, i));
	}
      }
      
      /*
       *  Points on the boundary of the parameter space are defined by
       *  $n . x = c$
       *  where $n$ is the normal vector and $c$ a constant. A line is of the form
       *  $x = u + v t$
       *  for $u, d$ vectors and t is the parameter, which we solve for:
       *  $n . (u + v t) = n.u + n.v t = c$
       *  $ => t = (c - n.u) / n.v        $
       *  which substituted back gives the point of intersection 
       *  $ => x = u + (c - n.u) / n.v * v$.
       */
      
      /* Compute dot products */
      gsl_blas_dgemv(CblasTrans, 1.0, tiling->latt_bound_normal, tiling->line_start, 0.0, tiling->line_start_dot);
      gsl_blas_dgemv(CblasTrans, 1.0, tiling->latt_bound_normal, tiling->line_dir, 0.0, tiling->line_dir_dot);

      /* Calculate the intersection of the line and bound */
      gsl_vector_memcpy(tiling->line_intersect, tiling->latt_bound_dot);
      gsl_vector_sub(tiling->line_intersect, tiling->line_start_dot);
      gsl_vector_div(tiling->line_intersect, tiling->line_dir_dot);
      gsl_vector_scale(tiling->line_intersect, gsl_vector_get(tiling->line_dir, tiling->line_index));
      gsl_vector_add_constant(tiling->line_intersect, gsl_vector_get(tiling->line_start, tiling->line_index));

      /* Iterate over the normal and origin vectors */
      for (j = 0; j < tiling->num_bounds; ++j) {

	/* 
	 *  Check the projection of the normal vector onto the line direction:
	 *    - if it is less than zero, the intersection is a lower bound;
	 *    - if it is greater than zero, the intersection is an upper bound;
	 *    - if it is zero, ignore it, as the line as parallel to the bound
	 */
	dir = gsl_vector_get(tiling->line_dir_dot, j);
	point = gsl_vector_get(tiling->line_intersect, j);
	if (dir < 0.0) {
 	  line_lower = GSL_MAX_DBL(line_lower, point);
	}
	else if (dir > 0.0) {
 	  line_upper = GSL_MIN_DBL(line_upper, point);
	}

      }

      /* Complain if bounds are still infinite */
      if (!(gsl_finite(line_lower) && gsl_finite(line_upper))) {
	LALPrintError("%s\nERROR: parameter space bounds appear to be insufficient\n", FLATLATTICETILINGC);
	XLAL_ERROR("XLALNextFlatLatticePoint", XLAL_EINVAL);
      }
      in_bounds = (line_lower <= line_upper);

      /* Set line bounds */
      gsl_vector_long_set(tiling->latt_lower, tiling->line_index, ceil(line_lower));
      gsl_vector_long_set(tiling->latt_upper, tiling->line_index, floor(line_upper)); 
      
      /* Move current point to lower line bound */
      gsl_vector_long_set(tiling->latt_current, tiling->line_index, gsl_vector_long_get(tiling->latt_lower, tiling->line_index));

    }

    /* Transform point from lattice to parameter space */
    for (i = 0; i < n; ++i) {
      point = gsl_vector_get(tiling->param_lower, i);
      for (j = 0; j < tiling->dimension; ++j) {
	point += gsl_matrix_get(tiling->latt_to_param, i, j) * gsl_vector_long_get(tiling->latt_current, j);
      }
      gsl_vector_set(tiling->current, i, point);
    }

    /* Final check to see if point is in parameter space (used for curved spaces) */
    if (tiling->in_param_space != NULL) {
      in_bounds = (tiling->in_param_space)(tiling);
    }

  }
  
  ++tiling->template_count;
  
  return XLAL_SUCCESS;
  
}
