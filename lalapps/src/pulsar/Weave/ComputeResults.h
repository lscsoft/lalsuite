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
/// \name Routines which compute coherent and semicoherent results
///
/// @{

///
/// Input data required for computing coherent results
///
typedef struct tagWeaveCohInput WeaveCohInput;

///
/// Results of a coherent computation on a single segment
///
typedef struct tagWeaveCohResults WeaveCohResults;

///
/// Partial results of a semicoherent computation in progress
///
typedef struct tagWeaveSemiPartials WeaveSemiPartials;

///
/// Final results of a semicoherent computation over many segments
///
typedef struct tagWeaveSemiResults WeaveSemiResults;

WeaveCohInput *XLALWeaveCohInputCreate(
  const BOOLEAN shortcut_compute,
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
int XLALWeaveSemiPartialsInit(
  WeaveSemiPartials **semi_parts,
  const BOOLEAN shortcut_compute,
  const LALStringVector *per_detectors,
  const UINT4 per_nsegments,
  const PulsarDopplerParams *semi_phys,
  const double dfreq,
  const UINT4 semi_nfreqs
  );
void XLALWeaveSemiPartialsDestroy(
  WeaveSemiPartials *semi_parts
  );
int XLALWeaveSemiPartialsAdd(
  WeaveSemiPartials *semi_parts,
  const WeaveCohResults *coh_res,
  const UINT4 coh_offset
  );
int XLALWeaveSemiResultsCompute(
  WeaveSemiResults **semi_res,
  const WeaveSemiPartials *semi_parts
  );
void XLALWeaveSemiResultsDestroy(
  WeaveSemiResults *semi_res
  );
int XLALWeaveFillOutputResultItem(
  WeaveOutputResultItem **item,
  BOOLEAN *full_init,
  const WeaveSemiResults *semi_res,
  const size_t nspins,
  const size_t freq_idx
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _COMPUTE_RESULTS_H
