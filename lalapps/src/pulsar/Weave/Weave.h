//
// Copyright (C) 2016, 2017 Karl Wette
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with with program; see the file COPYING. If not, write to the
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#ifndef _WEAVE_H
#define _WEAVE_H

///
/// \defgroup lalapps_pulsar_Weave Weave Search Application
/// \ingroup lalapps_pulsar_Apps
/// \author Karl Wette
///

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALStdlib.h>
#include <lal/LALStdio.h>
#include <lal/LALString.h>
#include <lal/PulsarDataTypes.h>
#include <lal/GSLHelpers.h>
#include <lal/FITSFileIO.h>

#include <LALAppsVCSInfo.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#ifdef __cplusplus
extern "C" {
#endif

///
/// Bitflags representing search simulation levels
///
typedef enum {
  /// Simulate search (implicitly with full memory allocation)
  WEAVE_SIMULATE                        = 0001,
  /// Simulate search with minimal memory allocation
  WEAVE_SIMULATE_MIN_MEM                = 0002,
} WeaveSimulationLevel;

#ifdef __cplusplus
}
#endif

#endif // _WEAVE_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
