/*********************************************************************************
 *  \author K. Wette
 *  \file
 *  \ingroup templateBanks
 *  \brief
 *  Frequency-spindown metric and parameter spaces for FlatLatticeTiling
 *********************************************************************************/

/******** Includes ********/

#include "DirectedSearch.h"

#include <math.h>

#include <lal/LALError.h>
#include <lal/XLALError.h>
#include <lal/LALConstants.h>
#include <lal/LALMalloc.h>

/******** Constants ********/

static const int TRUE  = (1==1);
static const int FALSE = (1==0);

/******** Declarations ********/

NRCSID(DIRECTEDSEARCHC, "$Id$");

REAL8 factorial(INT4);

BOOLEAN IsInBrakingIndexParameterSpace(FlatLatticeTiling*);

/******** Functions ********/

/* Private factorial function */
REAL8 factorial(INT4 n) {
  REAL8 f = 1.0;
  while (n > 0) {
    f *= n--;
  }
  return f;
}

/*
 *  Returns the frequency-spindown metric, derived in the same manner as PtoleMetric,
 *  (see that code for further background), except in terms of the frequency time
 *  derivatives expected by ComputeFstat, and scaled so that the mismatch can be
 *  calculated as $\mu = g_{ij} \Delta f^{(i)} \Delta f^{(j)}$. The time span-dependent
 *  components are included in the congruent factor, if included.
 */
int XLALSpindownMetric(
		       gsl_matrix* metric   , /**< [in/out] Matrix containing the congruent metric */
		       gsl_matrix* congruent, /**< [in/out] Optional matrix of the congruent factor */
		       REAL8 Tspan            /**< [in] Time span of the data set */
		       )
{

  UINT4 i, j;
  double v;

  /* Check matrix is square */
  if (metric->size1 != metric->size2) {
    LALPrintError("%s\nERROR: Metric must be square\n", DIRECTEDSEARCHC);
    XLAL_ERROR("XLALSpindownMetric", XLAL_EINVAL);
  }

  /* Check time span is non-zero */
  if (Tspan == 0.0) {
    LALPrintError("%s\nERROR: Time span must be greater than zero\n", DIRECTEDSEARCHC);
    XLAL_ERROR("XLALSpindownMetric", XLAL_EINVAL);
  }

  /* Generate metric */
  for (i = 0; i < metric->size1; ++i) {

    for (j = 0; j < metric->size2; ++j) {

      if (congruent != NULL) {
	gsl_matrix_set(congruent, i, j, (i == j ? pow(Tspan, 1 + i) : 0.0));
      }

      v = pow(LAL_PI, 2) / (factorial(i) * factorial(j) * (2 + i) * (2 + j) * (3 + i + j));
      if (congruent == NULL) {
	v *= pow(Tspan, 2 + i + j);
      }
      gsl_matrix_set(metric, i, j, v);

    }
  }

  return XLAL_SUCCESS;

}

/*
 *  Creates a braking index-based parameter space for FlatLatticeTiling
 */
int XLALBrakingIndexParameterSpace(
				   FlatLatticeTiling* tiling, /**< [in] Structure */
				   REAL8 freq,                /**< [in] Starting frequency */
				   REAL8 freq_band,           /**< [in] Frequency band */
				   REAL8 age,                 /**< [in] Characteristic age */
				   INT4 minBrakingIndex,      /**< [in] Minimum braking index */
				   INT4 maxBrakingIndex       /**< [in] Maximum braking index */
				   )
{

  gsl_vector *lower = NULL;
  gsl_vector *upper = NULL;

  /* Only up to 3 dimensions are currently supported */
  if (tiling->dimension > 3) {
    LALPrintError("%s\nERROR: Only up to 3 dimensions currently supported!\n", DIRECTEDSEARCHC);
    XLAL_ERROR("XLALBrakingIndexParameterSpace", XLAL_EINVAL);
  }

  /* Allocate variables */
  lower = gsl_vector_calloc(tiling->dimension);
  upper = gsl_vector_calloc(tiling->dimension);

  /* Set frequency bounds */
  gsl_vector_set(lower, 0, freq);
  gsl_vector_set(upper, 0, freq + freq_band);

  /* Set first spindown bounds */
  if (tiling->dimension > 1) {
    gsl_vector_set(lower, 1, -gsl_vector_get(upper, 0) / (minBrakingIndex - 1) / age);
    gsl_vector_set(upper, 1, -gsl_vector_get(lower, 0) / (maxBrakingIndex - 1) / age);
  }

  /* Set second spindown bounds */
  if (tiling->dimension > 2) {
    gsl_vector_set(lower, 2, minBrakingIndex * pow(gsl_vector_get(upper, 1), 2) / gsl_vector_get(upper, 0));
    gsl_vector_set(upper, 2, maxBrakingIndex * pow(gsl_vector_get(lower, 1), 2) / gsl_vector_get(lower, 0));
  }

  /* Create parameter space */
  XLALSquareParameterSpace(tiling, lower, upper);

  /* Add function to do curved bits */
  tiling->in_param_space = &IsInBrakingIndexParameterSpace;

  /* Extra data for function */
  if ((tiling->in_param_space_args = (double*) LALCalloc(3, sizeof(double))) == NULL) {
    LALPrintError("%s\nERROR: Could not allocate tiling->in_param_space_args\n", DIRECTEDSEARCHC);
    XLAL_ERROR("XLALBrakingIndexParameterSpace", XLAL_ENOMEM);
  }
  tiling->in_param_space_args[0] = age;
  tiling->in_param_space_args[1] = (double) minBrakingIndex;
  tiling->in_param_space_args[2] = (double) maxBrakingIndex;

  /* Cleanup */
  gsl_vector_free(lower);
  gsl_vector_free(upper);

  return XLAL_SUCCESS;

}

/*
 *  Checks that the point is within the braking index parameter space
 */
BOOLEAN IsInBrakingIndexParameterSpace(FlatLatticeTiling *tiling) {

  double *age = NULL;
  double *minBrakingIndex = NULL;
  double *maxBrakingIndex = NULL;

  double point = 0.0;
  double lower = 0.0;
  double upper = 0.0;

  /* Initialise pointers */
  age             = &(tiling->in_param_space_args[0]);
  minBrakingIndex = &(tiling->in_param_space_args[1]);
  maxBrakingIndex = &(tiling->in_param_space_args[2]);

  /* Check bounds on first spindown */
  if (tiling->dimension > 1) {
    point = gsl_vector_get(tiling->current, 1);
    lower = -gsl_vector_get(tiling->current, 0) / (*minBrakingIndex - 1) / *age;
    upper = -gsl_vector_get(tiling->current, 0) / (*maxBrakingIndex - 1) / *age;
    if (point < lower || upper < point) return FALSE;
  }

  /* Check bounds on second spindown */
  if (tiling->dimension > 2) {
    point = gsl_vector_get(tiling->current, 2);
    lower = *minBrakingIndex * pow(gsl_vector_get(tiling->current, 1), 2) / gsl_vector_get(tiling->current, 0);
    upper = *maxBrakingIndex * pow(gsl_vector_get(tiling->current, 1), 2) / gsl_vector_get(tiling->current, 0);
    if (point < lower || upper < point) return FALSE;
  }

  /* Point is in bounds */
  return TRUE;

}
