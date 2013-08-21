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

// This file implements the F-statistic demodulation algorithm. It is not compiled directly, but
// included from ComputeFstat.c

struct tagFstatInputData_Demod {
  MultiSFTVector *multiSFTs;			// Input multi-detector SFTs
  MultiDetectorStateSeries *multiDetStates;	// Multi-detector state series
  ComputeFParams params;			// Additional parameters for ComputeFStat()
  ComputeFBuffer buffer;			// Internal buffer for ComputeFStat()
};

static inline void
DestroyFstatInputData_Demod(
  FstatInputData_Demod* demod
  )
{
  XLALDestroyMultiSFTVector(demod->multiSFTs);
  XLALDestroyMultiDetectorStateSeries(demod->multiDetStates);
  XLALEmptyComputeFBuffer(&demod->buffer);
  XLALFree(demod);
}

FstatInputData*
XLALSetupFstat_Demod(
  MultiSFTVector **multiSFTs,
  MultiNoiseWeights **multiWeights,
  const EphemerisData *edat,
  const SSBprecision SSBprec,
  const DemodAMType demodAM,
  const UINT4 Dterms
  )
{

  // Check non-common input
  XLAL_CHECK_NULL(demodAM < DEMODAM_LAST, XLAL_EINVAL);
  XLAL_CHECK_NULL(Dterms > 0, XLAL_EINVAL);

  // Check common input and allocate input data struct
  FstatInputData* input = SetupFstat_Common(multiSFTs, multiWeights, edat, SSBprec);
  XLAL_CHECK_NULL(input != NULL, XLAL_EFUNC);

  // Allocate demodulation input data struct
  FstatInputData_Demod *demod = XLALCalloc(1, sizeof(FstatInputData_Demod));
  XLAL_CHECK_NULL(demod != NULL, XLAL_ENOMEM);

  // Save pointer to input SFTs, set supplied pointer to NULL
  demod->multiSFTs = *multiSFTs;
  *multiSFTs = NULL;

  // Calculate the detector states from the SFTs
  {
    LALStatus status = empty_status;
    LALGetMultiDetectorStates(&status, &demod->multiDetStates, demod->multiSFTs, edat);
    if (status.statusCode) {
      XLAL_ERROR_NULL(XLAL_EFAILED, "LALGetMultiDetectorStates() failed: %s (statusCode=%i)", status.statusDescription, status.statusCode);
    }
  }

  // Set parameters to pass to ComputeFStat()
  demod->params.Dterms = Dterms;
  demod->params.SSBprec = SSBprec;
  demod->params.buffer = NULL;
  demod->params.bufferedRAA = (demodAM & DEMODAM_BUFFERED_RIGID_ADIABATIC);
  demod->params.edat = edat;
  demod->params.upsampling = 1;
  demod->params.useRAA = (demodAM & DEMODAM_RIGID_ADIABATIC) || (demodAM & DEMODAM_BUFFERED_RIGID_ADIABATIC);

  // Save pointer to demodulation input data
  input->demod = demod;

  return input;

}

static int
ComputeFstat_Demod(
  FstatResults* Fstats,
  const FstatInputData_Common *common,
  FstatInputData_Demod* demod
  )
{

  // Check input
  XLAL_CHECK(Fstats != NULL, XLAL_EFAULT);
  XLAL_CHECK(common != NULL, XLAL_EFAULT);
  XLAL_CHECK(demod != NULL, XLAL_EFAULT);

  // Get which F-statistic quantities to compute
  const FstatQuantities whatToCompute = Fstats->whatWasComputed;

  // Check which quantities can be computed
  XLAL_CHECK(!(whatToCompute & FSTATQ_FAFB_PER_DET), XLAL_EINVAL, "Demodulation does not currently support Fa & Fb per detector");

  // Set parameters to pass to ComputeFStat()
  demod->params.returnSingleF = (whatToCompute & FSTATQ_2F_PER_DET);
  demod->params.returnAtoms = (whatToCompute & FSTATQ_ATOMS_PER_DET);

  // Save local copy of doppler point and starting frequency
  PulsarDopplerParams thisPoint = Fstats->doppler;
  const REAL8 fStart = thisPoint.fkdot[0];

  // Call ComputeFStat() for each frequency bin
  for (UINT4 k = 0; k < Fstats->numFreqBins; ++k) {

    // Set frequency to search at
    thisPoint.fkdot[0] = fStart + k * Fstats->dFreq;

    // Call ComputeFStat()
    Fcomponents Fcomp;
    {
      LALStatus status = empty_status;
      ComputeFStat(&status, &Fcomp, &thisPoint, demod->multiSFTs, common->multiWeights, demod->multiDetStates, &demod->params, &demod->buffer);
      if (status.statusCode) {
        XLAL_ERROR(XLAL_EFAILED, "ComputeFStat() failed: %s (statusCode=%i)", status.statusDescription, status.statusCode);
      }
    }

    // Return multi-detector 2F
    if (whatToCompute & FSTATQ_2F) {
      Fstats->twoF[k] = 2.0 * Fcomp.F;   // *** Return value of 2F ***
    }

    // Return multi-detector Fa & Fb
    if (whatToCompute & FSTATQ_FAFB) {
      Fstats->FaFb[k].Fa = Fcomp.Fa;
      Fstats->FaFb[k].Fb = Fcomp.Fb;
    }

    // Return 2F per detector
    if (whatToCompute & FSTATQ_2F_PER_DET) {
      for (UINT4 X = 0; X < Fstats->numDetectors; ++X) {
        Fstats->twoFPerDet[X][k] = 2.0 * Fcomp.FX[X];   // *** Return value of 2F ***
      }
    }

    // Return Fa & Fb per detector
    if (whatToCompute & FSTATQ_FAFB_PER_DET) {
      XLAL_ERROR(XLAL_EFAILED, "Unimplemented!");
    }

    // Return F-atoms per detector
    if (whatToCompute & FSTATQ_ATOMS_PER_DET) {
      XLALDestroyMultiFstatAtomVector(Fstats->multiFatoms[k]);
      Fstats->multiFatoms[k] = Fcomp.multiFstatAtoms;
    }

  } // for k < Fstats->numFreqBins

  // Return amplitude modulation coefficients
  if (demod->buffer.multiCmplxAMcoef != NULL) {
    Fstats->Mmunu = demod->buffer.multiCmplxAMcoef->Mmunu;
  } else if (demod->buffer.multiAMcoef != NULL) {
    Fstats->Mmunu.Ad = demod->buffer.multiAMcoef->Mmunu.Ad;
    Fstats->Mmunu.Bd = demod->buffer.multiAMcoef->Mmunu.Bd;
    Fstats->Mmunu.Cd = demod->buffer.multiAMcoef->Mmunu.Cd;
    Fstats->Mmunu.Ed = 0;
    Fstats->Mmunu.Dd = demod->buffer.multiAMcoef->Mmunu.Dd;
    Fstats->Mmunu.Sinv_Tsft = demod->buffer.multiAMcoef->Mmunu.Sinv_Tsft;
  }

  return XLAL_SUCCESS;

}
