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

// This file implements the F-statistic resampling algorithm. It is not compiled directly, but
// included from ComputeFstat.c

struct tagFstatInputData_Resamp {
  MultiSFTVector *multiSFTs;			// Input multi-detector SFTs
  ComputeFParams params;			// Additional parameters for ComputeFStatFreqBand_RS()
  REAL4 *Fout;					// Output array of *F* values passed to ComputeFStatFreqBand_RS()
};

static inline void
DestroyFstatInputData_Resamp(
  FstatInputData_Resamp* resamp
  )
{
  XLALDestroyMultiSFTVector(resamp->multiSFTs);
  XLALEmptyComputeFBuffer_RS(resamp->params.buffer);
  XLALFree(resamp->params.buffer);
  XLALFree(resamp->Fout);
  XLALFree(resamp);
}

FstatInputData*
XLALSetupFstat_Resamp(
  MultiSFTVector **multiSFTs,
  MultiNoiseWeights **multiWeights,
  const EphemerisData *edat,
  const SSBprecision SSBprec
  )
{

  // Check common input and allocate input data struct
  FstatInputData* input = SetupFstat_Common(multiSFTs, multiWeights, edat, SSBprec);
  XLAL_CHECK_NULL(input != NULL, XLAL_EFUNC);

  // Allocate resampling input data struct
  FstatInputData_Resamp *resamp = XLALCalloc(1, sizeof(FstatInputData_Resamp));
  XLAL_CHECK_NULL(resamp != NULL, XLAL_ENOMEM);

  // Save pointer to input SFTs, set supplied pointer to NULL
  resamp->multiSFTs = *multiSFTs;
  multiSFTs = NULL;

  // Set parameters to pass to ComputeFStatFreqBand_RS()
  resamp->params.Dterms = 0;
  resamp->params.SSBprec = SSBprec;
  resamp->params.buffer = NULL;
  resamp->params.bufferedRAA = 0;
  resamp->params.edat = edat;
  resamp->params.upsampling = 1;
  resamp->params.useRAA = 0;

  // Initialise output array of *F* values to NULL
  resamp->Fout = NULL;

  // Save pointer to resampling input data
  input->resamp = resamp;

  return input;

}

static int
ComputeFstat_Resamp(
  FstatResults* Fstats,
  const FstatInputData_Common *common,
  FstatInputData_Resamp* resamp
  )
{

  // Check input
  XLAL_CHECK(Fstats != NULL, XLAL_EFAULT);
  XLAL_CHECK(common != NULL, XLAL_EFAULT);
  XLAL_CHECK(resamp != NULL, XLAL_EFAULT);

  // Get which F-statistic quantities to compute
  const FstatQuantities whatToCompute = Fstats->whatWasComputed;

  // Check which quantities can be computed
  XLAL_CHECK(!(whatToCompute & FSTATQ_FAFB), XLAL_EINVAL, "Resamping does not currently support Fa & Fb");
  XLAL_CHECK(!(whatToCompute & FSTATQ_FAFB_PER_DET), XLAL_EINVAL, "Resamping does not currently support Fa & Fb per detector");
  XLAL_CHECK(!(whatToCompute & FSTATQ_ATOMS_PER_DET), XLAL_EINVAL, "Resamping does not currently support atoms per detector");

  // Set parameters to pass to ComputeFStatFreqBand_RS()
  resamp->params.returnSingleF = whatToCompute & FSTATQ_2F_PER_DET;
  resamp->params.returnAtoms = 0;

  // Save local copy of doppler point
  PulsarDopplerParams thisPoint = Fstats->doppler;

  // (Re)allocate output array of *F* values
  const UINT4 FoutN = resamp->params.returnSingleF ? (Fstats->numDetectors + 1) : 1;
  resamp->Fout = XLALRealloc(resamp->Fout, Fstats->numFreqBins * FoutN * sizeof(resamp->Fout[0]));
  XLAL_CHECK(resamp->Fout != NULL, XLAL_ENOMEM);

  // Create REAL4FrequencySeries to receive 2F values
  static const REAL4Sequence empty_REAL4Sequence;
  REAL4Sequence CFSFB_RS_data = empty_REAL4Sequence;
  CFSFB_RS_data.length = Fstats->numFreqBins * FoutN;
  CFSFB_RS_data.data = resamp->Fout;
  static const REAL4FrequencySeries empty_REAL4FrequencySeries;
  REAL4FrequencySeries CFSFB_RS = empty_REAL4FrequencySeries;
  CFSFB_RS.deltaF = Fstats->dFreq;
  CFSFB_RS.data = &CFSFB_RS_data;

  // Call ComputeFStatFreqBand_RS()
  {
    LALStatus status = empty_status;
    ComputeFStatFreqBand_RS(&status, &CFSFB_RS, &thisPoint, resamp->multiSFTs, common->multiWeights, &resamp->params);
    if (status.statusCode) {
      XLAL_ERROR(XLAL_EFAILED, "ComputeFStatFreqBand_RS() failed: %s (statusCode=%i)", status.statusDescription, status.statusCode);
    }
  }
  for (UINT4 k = 0; k < Fstats->numFreqBins; ++k) {

    // Return multi-detector 2F
    if (whatToCompute & FSTATQ_2F) {
      Fstats->twoF[k] = 2.0 * resamp->Fout[k];   // *** Return value of 2F ***
    }

    // Return multi-detector Fa & Fb
    if (whatToCompute & FSTATQ_FAFB) {
      XLAL_ERROR(XLAL_EFAILED, "Unimplemented!");
    }

    // Return 2F per detector
    if (whatToCompute & FSTATQ_2F_PER_DET) {
      for (UINT4 X = 0; X < Fstats->numDetectors; ++X) {
        Fstats->twoFPerDet[X][k] = 2.0 * resamp->Fout[(X+1)*Fstats->numFreqBins + k];   // *** Return value of 2F ***
      }
    }

    // Return Fa & Fb per detector
    if (whatToCompute & FSTATQ_FAFB_PER_DET) {
      XLAL_ERROR(XLAL_EFAILED, "Unimplemented!");
    }

    // Return F-atoms per detector
    if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
      XLAL_ERROR(XLAL_EFAILED, "Unimplemented!");
    }

  } // for k < Fstats->numFreqBins

  // Resampling cannot currently return amplitude modulation coefficients
  Fstats->Mmunu.Ad = NAN;
  Fstats->Mmunu.Bd = NAN;
  Fstats->Mmunu.Cd = NAN;
  Fstats->Mmunu.Ed = NAN;
  Fstats->Mmunu.Dd = NAN;
  Fstats->Mmunu.Sinv_Tsft = NAN;

  return XLAL_SUCCESS;

}
