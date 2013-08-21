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
  MultiSFTVector **multiSFTs,
  MultiNoiseWeights **multiWeights,
  const EphemerisData *edat,
  const SSBprecision SSBprec
  )
{

  // Check input
  XLAL_CHECK_NULL(multiSFTs != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(*multiSFTs != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(multiWeights != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(edat != NULL, XLAL_EFAULT);
  XLAL_CHECK_NULL(SSBprec < SSBPREC_LAST, XLAL_EINVAL);

  // Check number of SFTs
  XLAL_CHECK_NULL((*multiSFTs)->length > 0, XLAL_EINVAL, "Found no SFTs!");
  XLAL_CHECK_NULL((*multiSFTs)->length <= PULSAR_MAX_DETECTORS, XLAL_EINVAL, "Supports only up to PULSAR_MAX_DETECTORS=%u detectors", PULSAR_MAX_DETECTORS);
  for (UINT4 X = 0; X < (*multiSFTs)->length; ++X) {

    // Check number of SFTs for each detector
    XLAL_CHECK_NULL((*multiSFTs)->data[X]->length > 0, XLAL_EINVAL, "Found no SFTs from detector %u", X);
    for (UINT4 alpha = 0; alpha < (*multiSFTs)->data[X]->length; ++alpha) {

      // Check length of SFTs
      XLAL_CHECK_NULL((*multiSFTs)->data[X]->data[alpha].data->length > 0, XLAL_EINVAL,
                      "Found zero-length SFT from detector %u, position %u", X, alpha);

    }

  }

  // If noise weights were supplied ...
  if (*multiWeights != NULL) {

    // Check numbers of noise weights match SFTs
    XLAL_CHECK_NULL((*multiWeights)->length == (*multiSFTs)->length, XLAL_EINVAL,
                    "Number of noise weight detectors does not match SFTS: %u != %u",
                    (*multiWeights)->length, (*multiSFTs)->length);
    for (UINT4 X = 0; X < (*multiSFTs)->length; ++X) {

      // Check number of noise weights for each detector
      XLAL_CHECK_NULL((*multiWeights)->data[X]->length == (*multiSFTs)->data[X]->length, XLAL_EINVAL,
                      "Number of noise weights from detector %u does not match SFTS: %u != %u",
                      X, (*multiWeights)->data[X]->length, (*multiSFTs)->data[X]->length);

    }

  }

  // Allocate input data struct
  FstatInputData* input = XLALCalloc(1, sizeof(*input));
  XLAL_CHECK_NULL(input != NULL, XLAL_ENOMEM);

  // Save number of detectors, and copy name of each detector
  input->common.numDetectors = (*multiSFTs)->length;
  for (UINT4 X = 0; X < (*multiSFTs)->length; ++X) {
    strncpy(input->common.detectorNames[X], (*multiSFTs)->data[X]->data[0].name, 2);
  }

  // If no noise weights were supplied ...
  if (*multiWeights == NULL) {

    // Correction to F-statistic quantities computed by XLALComputeFstat() without noise weights
    const REAL8 Tsft = 1.0 / (*multiSFTs)->data[0]->data[0].deltaF;
    input->common.Fnorm = 1.0 / sqrt( 0.5 * Tsft );

  } else {
    input->common.Fnorm = 0.0;
  }

  // Save pointer to input noise weights, set supplied pointer to NULL
  input->common.multiWeights = *multiWeights;
  *multiWeights = NULL;

  return input;

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
  FstatResults **Fstats,
  FstatInputData *input,
  const PulsarDopplerParams *doppler,
  const REAL8 dFreq,
  const UINT4 numFreqBins,
  const FstatQuantities whatToCompute
  )
{

  // Check input
  XLAL_CHECK(Fstats != NULL, XLAL_EFAULT);
  XLAL_CHECK(input != NULL, XLAL_EFAULT);
  XLAL_CHECK(doppler != NULL, XLAL_EFAULT);
  XLAL_CHECK(doppler->orbit == NULL, XLAL_EINVAL, "Binary parameters are currently not supported!");
  XLAL_CHECK(dFreq > 0 || (numFreqBins == 1 && dFreq >= 0), XLAL_EINVAL);
  XLAL_CHECK(numFreqBins > 0, XLAL_EINVAL);
  XLAL_CHECK(0 < whatToCompute && whatToCompute < FSTATQ_LAST, XLAL_EINVAL);

  // Allocate results struct, if needed
  if (*Fstats == NULL) {
    *Fstats = XLALCalloc(1, sizeof(**Fstats));
    XLAL_CHECK(*Fstats != NULL, XLAL_ENOMEM);
  }

  // Get constant pointer to common input data
  const FstatInputData_Common *common = &(input->common);
  const UINT4 numDetectors = common->numDetectors;

  // Enlarge result arrays if they are too small
  const bool moreFreqBins = (numFreqBins > (*Fstats)->internalalloclen);
  const bool moreDetectors = (numDetectors > (*Fstats)->numDetectors);
  if (moreFreqBins || moreDetectors) {

    // Enlarge multi-detector 2F array
    if ((whatToCompute & FSTATQ_2F) && moreFreqBins) {
      (*Fstats)->twoF = XLALRealloc((*Fstats)->twoF, numFreqBins*sizeof((*Fstats)->twoF[0]));
      XLAL_CHECK((*Fstats)->twoF != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->twoF to length %u", numFreqBins);
    }

    // Enlarge multi-detector Fa & Fb array
    if ((whatToCompute & FSTATQ_FAFB) && moreFreqBins) {
      (*Fstats)->FaFb = XLALRealloc((*Fstats)->FaFb, numFreqBins*sizeof((*Fstats)->FaFb[0]));
      XLAL_CHECK((*Fstats)->FaFb != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->FaFb to length %u", numFreqBins);
    }

    // Enlarge 2F per detector arrays
    if ((whatToCompute & FSTATQ_2F_PER_DET) && (moreFreqBins || moreDetectors)) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        (*Fstats)->twoFPerDet[X] = XLALRealloc((*Fstats)->twoFPerDet[X], numFreqBins*sizeof((*Fstats)->twoFPerDet[X][0]));
        XLAL_CHECK((*Fstats)->twoFPerDet[X] != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->twoFPerDet[%u] to length %u", X, numFreqBins);
      }
    }

    // Enlarge Fa & Fb per detector arrays
    if ((whatToCompute & FSTATQ_FAFB_PER_DET) && (moreFreqBins || moreDetectors)) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        (*Fstats)->FaFbPerDet[X] = XLALRealloc((*Fstats)->FaFbPerDet[X], numFreqBins*sizeof((*Fstats)->FaFbPerDet[X][0]));
        XLAL_CHECK((*Fstats)->FaFbPerDet[X] != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->FaFbPerDet[%u] to length %u", X, numFreqBins);
      }
    }

    // Enlarge F-atoms per detector arrays, and initialise to NULL
    if ((whatToCompute & FSTATQ_ATOMS_PER_DET) && (moreFreqBins || moreDetectors)) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        (*Fstats)->multiFatoms = XLALRealloc((*Fstats)->multiFatoms, numFreqBins*sizeof((*Fstats)->multiFatoms[0]));
        XLAL_CHECK((*Fstats)->multiFatoms != NULL, XLAL_EINVAL, "Failed to (re)allocate (*Fstats)->multiFatoms to length %u", numFreqBins);

        // If more detectors are needed, destroy multi-F-atom vectors so they can be re-allocated later
        if (moreDetectors) {
          for (UINT4 k = 0; k < numFreqBins; ++k) {
            XLALDestroyMultiFstatAtomVector((*Fstats)->multiFatoms[k]);
            (*Fstats)->multiFatoms[k] = NULL;
          }
        } else {
          for (UINT4 k = (*Fstats)->internalalloclen; k < numFreqBins; ++k) {
            (*Fstats)->multiFatoms[k] = NULL;
          }
        }

      }
    }

    // Update allocated length of arrays
    (*Fstats)->internalalloclen = numFreqBins;

  } // if (moreFreqBins || moreDetectors)

  // Initialise result struct parameters
  (*Fstats)->doppler = *doppler;
  (*Fstats)->dFreq = dFreq;
  (*Fstats)->numFreqBins = numFreqBins;
  (*Fstats)->numDetectors = numDetectors;
  memcpy((*Fstats)->detectorNames, common->detectorNames, sizeof(common->detectorNames));
  (*Fstats)->whatWasComputed = whatToCompute;

  // Call the appropriate algorithm function to compute the F-statistic
  if (input->demod != NULL) {
    XLAL_CHECK(ComputeFstat_Demod(*Fstats, common, input->demod) == XLAL_SUCCESS, XLAL_EFUNC);
  } else if (input->resamp != NULL) {
    XLAL_CHECK(ComputeFstat_Resamp(*Fstats, common, input->resamp) == XLAL_SUCCESS, XLAL_EFUNC);
  } else {
    XLAL_ERROR(XLAL_EFAILED, "Invalid FstatInputData struct passed to %s()", __func__);
  }

  // Correct F-statistic quantities when no noise weights are given
  if (common->Fnorm != 0.0) {
    const REAL8 Fnorm = common->Fnorm;
    const REAL8 Fnorm_sqr = Fnorm * Fnorm;

    // Correct antenna pattern matrix
    (*Fstats)->Mmunu.Sinv_Tsft = 2.0 / Fnorm_sqr;   // equivalent to Tsft

    // Correct multi-detector 2F array
    if (whatToCompute & FSTATQ_2F) {
      for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
        (*Fstats)->twoF[k] *= Fnorm_sqr;
        (*Fstats)->twoF[k] += 4;
      }
    }

    // Correct multi-detector F-parts array
    if (whatToCompute & FSTATQ_FAFB) {
      for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
        (*Fstats)->FaFb[k].Fa *= Fnorm;
        (*Fstats)->FaFb[k].Fb *= Fnorm;
      }
    }

    // Correct 2F per detector arrays
    if (whatToCompute & FSTATQ_2F_PER_DET) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
          (*Fstats)->twoFPerDet[X][k] *= Fnorm_sqr;
          (*Fstats)->twoFPerDet[X][k] += 4;
        }
      }
    }

    // Correct F-parts per detector arrays
    if (whatToCompute & FSTATQ_FAFB_PER_DET) {
      for (UINT4 X = 0; X < numDetectors; ++X) {
        for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
          (*Fstats)->FaFbPerDet[X][k].Fa *= Fnorm;
          (*Fstats)->FaFbPerDet[X][k].Fb *= Fnorm;
        }
      }
    }

    // Correct F-atoms per detector arrays, and initialise to NULL
    if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
      for (UINT4 k = 0; k < (*Fstats)->numFreqBins; ++k) {
        for (UINT4 X = 0; X < numDetectors; ++X) {
          FstatAtomVector *atomX = (*Fstats)->multiFatoms[k]->data[X];
          for (UINT4 alpha = 0; alpha < atomX->length; ++alpha) {
            atomX->data[alpha].Fa_alpha *= Fnorm;
            atomX->data[alpha].Fb_alpha *= Fnorm;
          }
        }
      }
    }

  } // if (common->Fnorm != 0.0)

  return XLAL_SUCCESS;

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
  FstatResults* Fstats
  )
{
  if (Fstats != NULL) {
    XLALFree(Fstats->twoF);
    XLALFree(Fstats->FaFb);
    for (UINT4 X = 0; X < PULSAR_MAX_DETECTORS; ++X) {
      XLALFree(Fstats->twoFPerDet[X]);
      XLALFree(Fstats->FaFbPerDet[X]);
      if (Fstats->multiFatoms != NULL) {
        for (UINT4 n = 0; n < Fstats->internalalloclen; ++n) {
          XLALDestroyMultiFstatAtomVector(Fstats->multiFatoms[n]);
        }
        XLALFree(Fstats->multiFatoms);
      }
    }
    XLALFree(Fstats);
  }
}
