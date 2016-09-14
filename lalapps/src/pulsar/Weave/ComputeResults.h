//
// Copyright (C) 2016 Karl Wette
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

#ifndef _COMPUTE_RESULTS_H
#define _COMPUTE_RESULTS_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"

#include <lal/ComputeFstat.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Input data required for computing coherent results
///
typedef struct tagWeaveCohInput WeaveCohInput;

///
/// Results of a coherent computation on a single segment
///
typedef struct tagWeaveCohResults WeaveCohResults;

///
/// Results of a semicoherent computation over many segments
///
typedef struct tagWeaveSemiResults WeaveSemiResults;

///
/// \name Routines which compute coherent and semicoherent results
///
/// @{

WeaveCohInput *XLALWeaveCohInputCreate(
  FstatInput *Fstat_input,
  const LALStringVector *per_detectors
  );
void XLALWeaveCohInputDestroy(
  WeaveCohInput *coh_input
  );
int XLALWeaveCohResultsCompute(
  WeaveCohResults **coh_res,
  WeaveCohInput *coh_input,
  const PulsarDopplerParams *coh_phys,
  const UINT4 coh_nfreqs
  );
void XLALWeaveCohResultsDestroy(
  WeaveCohResults *coh_res
  );
WeaveSemiResults *XLALWeaveSemiResultsCreate(
  const LALStringVector *per_detectors,
  const UINT4 per_nsegments,
  const double dfreq
  );
void XLALWeaveSemiResultsDestroy(
  WeaveSemiResults *semi_res
  );
int XLALWeaveSemiResultsInit(
  WeaveSemiResults *semi_res,
  const PulsarDopplerParams *semi_phys,
  const UINT4 semi_nfreqs
  );
int XLALWeaveSemiResultsAdd(
  WeaveSemiResults *semi_res,
  const WeaveCohResults *coh_res,
  const UINT4 coh_offset
  );
int XLALWeaveFillOutputResultItem(
  WeaveOutputResultItem **item,
  BOOLEAN *full_init,
  const WeaveSemiResults *semi_res,
  const size_t freq_idx
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _COMPUTE_RESULTS_H
