/*********************************************************************************
 *  \author K. Wette
 *  \file
 *  \ingroup templateBanks
 *  \brief
 *  Frequency-spindown metric and parameter spaces for FlatLatticeTiling
 *********************************************************************************/

#ifndef _DIRECTED_H
#define _DIRECTED_H

/******** Includes ********/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALRCSID.h>
#include <lal/LALDatatypes.h>

#include "FlatLatticeTiling.h"

/******** Declarations ********/

#ifdef __cplusplus
extern "C" {
#endif

  NRCSID(DIRECTEDSEARCHH, "$Id$");

  int XLALSpindownMetric(gsl_matrix*, gsl_vector*, REAL8);

  int XLALBrakingIndexParameterSpace(FlatLatticeTiling*, REAL8, REAL8, REAL8, INT4, INT4);

#ifdef __cplusplus
}
#endif

#endif
