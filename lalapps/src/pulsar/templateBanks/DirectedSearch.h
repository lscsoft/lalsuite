/**
 * \author K. Wette
 * \file
 * \brief Frequency-spindown metric
 */

#ifndef _DIRECTED_H
#define _DIRECTED_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALRCSID.h>
#include <lal/LALDatatypes.h>

#include "FlatLatticeTiling.h"

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(DIRECTEDSEARCHH, "$Id$");

  gsl_matrix *XLALSpindownMetric(UINT4, REAL8);

#ifdef __cplusplus
}
#endif

#endif
