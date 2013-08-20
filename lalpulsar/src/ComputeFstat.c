//
// Copyright (C) 2012, 2013 Karl Wette
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

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include <lal/ComputeFstat.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

static const LALStatus empty_status;

///// Internal function prototypes /////

// Performs some common input sanity checks, and allocate a FstatInputData* struct
static FstatInputData*
SetupFstat_Common(
  MultiSFTVector **multiSFTs,
  MultiNoiseWeights **multiWeights,
  const EphemerisData *edat,
  const SSBprecision SSBprec
  );

///// Internal struct definitions /////

// Common input data for F-statistic algorithms
typedef struct {
  UINT4 numDetectors;					// Number of detectors
  CHAR detectorNames[PULSAR_MAX_DETECTORS][3];		// Names of detectors
  REAL8 Fnorm;						// F-statistic normalisation factor
  MultiNoiseWeights *multiWeights;			// Multi-detector noise weights
} FstatInputData_Common;

// Input data specific to F-statistic algorithms
typedef struct tagFstatInputData_Demod FstatInputData_Demod;
typedef struct tagFstatInputData_Resamp FstatInputData_Resamp;

// Internal definition of input data structure
struct tagFstatInputData {
  FstatInputData_Common common;				// Common input data
  FstatInputData_Demod* demod;				// Demodulation input data
  FstatInputData_Resamp* resamp;			// Resampling input data
};

///// Include F-statistic algorithm implementations /////

#include "ComputeFstat_Demod.c"
#include "ComputeFstat_Resamp.c"

///// Function definitions /////

static FstatInputData*
SetupFstat_Common(
  UNUSED MultiSFTVector **multiSFTs,
  UNUSED MultiNoiseWeights **multiWeights,
  UNUSED const EphemerisData *edat,
  UNUSED const SSBprecision SSBprec
  )
{
  XLAL_ERROR_NULL( XLAL_EFAILED, "Unimplemented!" );
}

FstatInputDataVector*
XLALCreateFstatInputDataVector(
  const UINT4 length
  )
{

  // Allocate and initialise vector container
  FstatInputDataVector* inputs = XLALCalloc(1, sizeof(*inputs));
  XLAL_CHECK_NULL(inputs != NULL, XLAL_ENOMEM);
  inputs->length = length;

  // Allocate and initialise vector data
  if (inputs->length > 0) {
    inputs->data = XLALCalloc(inputs->length, sizeof(inputs->data[0]));
    XLAL_CHECK_NULL(inputs->data != NULL, XLAL_ENOMEM);
  }

  return inputs;

}

void
XLALDestroyFstatInputDataVector(
  FstatInputDataVector* inputs
  )
{
  if (inputs != NULL) {
    for (UINT4 i = 0; i < inputs->length; ++i) {
      XLALDestroyFstatInputData(inputs->data[i]);
    }
    XLALFree(inputs->data);
    XLALFree(inputs);
  }
}

int
XLALComputeFstat(
  UNUSED FstatResults **Fstats,
  UNUSED FstatInputData *input,
  UNUSED const PulsarDopplerParams *doppler,
  UNUSED const REAL8 dFreq,
  UNUSED const UINT4 numFreqBins,
  UNUSED const FstatQuantities whatToCompute
  )
{
  XLAL_ERROR( XLAL_EFAILED, "Unimplemented!" );
}

void
XLALDestroyFstatInputData(
  FstatInputData* input
  )
{
  if (input != NULL) {
    XLALDestroyMultiNoiseWeights(input->common.multiWeights);
    if (input->demod != NULL) {
      DestroyFstatInputData_Demod(input->demod);
    } else if (input->resamp != NULL) {
      DestroyFstatInputData_Resamp(input->resamp);
    }
    XLALFree(input);
  }
}

void
XLALDestroyFstatResults(
  UNUSED FstatResults* Fstats
  )
{
  XLAL_ERROR_VOID( XLAL_EFAILED, "Unimplemented!" );
}
