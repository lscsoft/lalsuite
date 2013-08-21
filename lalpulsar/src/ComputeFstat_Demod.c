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
};

static inline void
DestroyFstatInputData_Demod(
  UNUSED FstatInputData_Demod* demod
  )
{
  XLAL_ERROR_VOID( XLAL_EFAILED, "Unimplemented!" );
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
  UNUSED FstatInputData* input = SetupFstat_Common(multiSFTs, multiWeights, edat, SSBprec);

  XLAL_ERROR_NULL( XLAL_EFAILED, "Unimplemented!" );

}

static int
ComputeFstat_Demod(
  UNUSED FstatResults* Fstats,
  UNUSED const FstatInputData_Common *common,
  UNUSED FstatInputData_Demod* demod
  )
{
  XLAL_ERROR( XLAL_EFAILED, "Unimplemented!" );
}
