/**
 * \author K. Wette
 * \file
 * \brief Frequency-spindown metric
 */

#include "DirectedSearch.h"

#include <math.h>

#include <lal/LALError.h>
#include <lal/XLALError.h>
#include <lal/LALConstants.h>
#include <lal/LALMalloc.h>

NRCSID(DIRECTEDSEARCHC, "$Id$");

static REAL8 factrl(INT4 n) {
  REAL8 f = 1.0;
  while (n > 0) {
    f *= n--;
  }
  return f;
}

/**
 * Frequency and frequency derivative components of the metric, suitable for a directed
 * search with only one fixed sky position. The units are those expected by ComputeFstat.
 */
gsl_matrix *XLALSpindownMetric(
			       UINT4 dimension, /**< [in] Dimension of the lattice */
 			       REAL8 Tspan      /**< [in] Time span of the data set */
			       )
{

  UINT4 i, j;
  gsl_matrix *metric = NULL;

  /* Allocate metric */
  metric = gsl_matrix_alloc(dimension, dimension);

  /* Calculate metric */
  for (i = 0; i < metric->size1; ++i) {
    for (j = 0; j < metric->size2; ++j) {

      gsl_matrix_set(metric, i, j, 4 * pow(LAL_PI, 2) * pow(Tspan, 2 + i + j) / 
		     (factrl(i) * factrl(j) * (2 + i) * (2 + j) * (3 + i + j)));

    }
  }

  return metric;

}
