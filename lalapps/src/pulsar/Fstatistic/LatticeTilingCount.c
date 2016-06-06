//
// Copyright (C) 2007, 2008, 2016 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA  02111-1307  USA
//

///
/// \file
/// \ingroup lalapps_pulsar_Fstatistic
/// \author Karl Wette
/// \brief Count number of templates in a given lattice tiling
///

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>

#include <lal/LALStdlib.h>
#include <lal/LALString.h>
#include <lal/XLALError.h>
#include <lal/UserInput.h>
#include <lal/LatticeTiling.h>
#include <lal/DopplerFullScan.h>
#include <lal/PtoleMetric.h>
#include <lal/GSLHelpers.h>

#include <LALAppsVCSInfo.h>

typedef struct {
  REAL8 time_span;
  REAL8Vector *square;
  REAL8Vector *age_braking;
  REAL8 max_mismatch;
  CHAR *lattice;
  CHAR *metric;
} UserVariables;

int main(int argc, char *argv[])
{

  // Check VCS information
  XLALAppsVCSInfoCheck();

  // Initialise user variables
  UserVariables uvar_struct = {
    .lattice = XLALStringDuplicate("an-star"),
    .metric = XLALStringDuplicate("spindown"),
  };
  UserVariables *const uvar = &uvar_struct;

  // Register user variables
  XLAL_CHECK_MAIN(XLALRegisterUvarMember(time_span, REAL8, 'T', REQUIRED, "Time-span of the data set (in seconds)") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN(XLALRegisterUvarMember(square, REAL8Vector, 0, OPTIONAL, "Square parameter space: start,width,...") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN(XLALRegisterUvarMember(age_braking, REAL8Vector, 0, OPTIONAL, "Age/braking index parameter space: alpha,delta,freq,freqband,age,minbrake,maxbrake") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN(XLALRegisterUvarMember(max_mismatch, REAL8, 'X', REQUIRED, "Maximum allowed mismatch between the templates") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN(XLALRegisterUvarMember(lattice, STRING, 'L', REQUIRED, "Lattice: 'an-star' or 'cubic'") == XLAL_SUCCESS, XLAL_EFUNC);
  XLAL_CHECK_MAIN(XLALRegisterUvarMember(metric, STRING, 'M', OPTIONAL, "Metric: 'spindown' or 'eye'") == XLAL_SUCCESS, XLAL_EFUNC);

  // Parse user input
  BOOLEAN should_exit = 0;
  XLAL_CHECK( XLALUserVarReadAllInput( &should_exit, argc, argv ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Check user input
  XLALUserVarCheck( &should_exit, UVAR_SET2(square, age_braking) == 1, "Exactly one of " UVAR_STR2AND(square, age_braking) " must be specified" );

  // Exit if required
  if ( should_exit ) {
    return EXIT_FAILURE;
  }

  LatticeTiling *tiling = NULL;
  if (UVAR_SET(square)) {

    // Create square parameter space
    XLAL_CHECK_MAIN(GSL_IS_EVEN(uvar->square->length), XLAL_EINVAL, "'square' must have an even number of arguments");
    const size_t n = uvar->square->length/2;
    tiling = XLALCreateLatticeTiling(n);
    XLAL_CHECK_MAIN(tiling != NULL, XLAL_EFUNC);
    for (size_t i = 0; i < n; ++i) {
      XLAL_CHECK_MAIN(XLALSetLatticeTilingConstantBound(tiling, i, uvar->square->data[2*i], uvar->square->data[2*i] + uvar->square->data[2*i + 1]) == XLAL_SUCCESS, XLAL_EFUNC);
    }

  } else {

    // Create age--braking index parameter space
    XLAL_CHECK_MAIN(uvar->age_braking->length == 7, XLAL_EINVAL, "'age-braking' must have exactly 7 arguments");
    const double alpha       = uvar->age_braking->data[0];
    const double delta       = uvar->age_braking->data[1];
    const double freq        = uvar->age_braking->data[2];
    const double freq_band   = uvar->age_braking->data[3];
    const double age         = uvar->age_braking->data[4];
    const double min_braking = uvar->age_braking->data[5];
    const double max_braking = uvar->age_braking->data[6];
    tiling = XLALCreateLatticeTiling(5);
    XLAL_CHECK_MAIN(tiling != NULL, XLAL_EFUNC);
    XLAL_CHECK_MAIN(XLALSetLatticeTilingConstantBound(tiling, 0, alpha, alpha) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN(XLALSetLatticeTilingConstantBound(tiling, 1, delta, delta) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN(XLALSetLatticeTilingConstantBound(tiling, 2, freq, freq + freq_band) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN(XLALSetLatticeTilingF1DotAgeBrakingBound(tiling, 2, 3, age, min_braking, max_braking) == XLAL_SUCCESS, XLAL_EFUNC);
    XLAL_CHECK_MAIN(XLALSetLatticeTilingF2DotBrakingBound(tiling, 2, 3, 4, min_braking, max_braking) == XLAL_SUCCESS, XLAL_EFUNC);

  }
  const size_t n = XLALTotalLatticeTilingDimensions(tiling);

  // Set lattice and metric
  gsl_matrix *metric = NULL;
  if (XLALStringCaseCompare(uvar->metric, "spindown") == 0) {
    GAMAT(metric, n, n);
    gsl_matrix_set_identity(metric);
    gsl_matrix_view spin_metric = gsl_matrix_submatrix(metric, 2, 2, n - 2, n - 2);
    XLAL_CHECK_MAIN(XLALSpindownMetric(&spin_metric.matrix, uvar->time_span) == XLAL_SUCCESS, XLAL_EFUNC);
  } else if (XLALStringCaseCompare(uvar->metric, "eye") == 0) {
    GAMAT(metric, n, n);
    gsl_matrix_set_identity(metric);
  } else {
    XLAL_ERROR_MAIN(XLAL_EINVAL, "Invalid value '%s' for 'metric'", uvar->metric);
  }
  XLAL_CHECK_MAIN(XLALSetTilingLatticeAndMetric(tiling, uvar->lattice, metric, uvar->max_mismatch)  == XLAL_SUCCESS, XLAL_EFUNC);
  gsl_matrix_free(metric);

  // Create a lattice iterator
  LatticeTilingIterator *itr = XLALCreateLatticeTilingIterator(tiling, n);
  XLAL_CHECK_MAIN(itr != NULL, XLAL_EFUNC);

  // Print number of templates
  UINT8 ntemplates = XLALTotalLatticeTilingPoints(itr);
  XLAL_CHECK_MAIN(ntemplates > 0, XLAL_EFUNC);
  printf("%" LAL_UINT8_FORMAT "\n", ntemplates);

  // Cleanup
  XLALDestroyLatticeTilingIterator(itr);
  XLALDestroyLatticeTiling(tiling);
  XLALDestroyUserVars();

  LALCheckMemoryLeaks();

  return EXIT_SUCCESS;

}
