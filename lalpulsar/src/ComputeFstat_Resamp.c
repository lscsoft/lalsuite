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
};

static inline void
DestroyFstatInputData_Resamp(
  UNUSED FstatInputData_Resamp* resamp
  )
{
  XLAL_ERROR_VOID( XLAL_EFAILED, "Unimplemented!" );
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
  UNUSED FstatInputData* input = SetupFstat_Common(multiSFTs, multiWeights, edat, SSBprec);

  XLAL_ERROR_NULL( XLAL_EFAILED, "Unimplemented!" );

}
