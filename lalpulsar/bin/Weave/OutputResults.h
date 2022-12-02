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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301 USA
//

#ifndef _OUTPUT_RESULTS_H
#define _OUTPUT_RESULTS_H

///
/// \file
/// \ingroup lalpulsar_bin_Weave
/// \brief Module which handles the output results
///

#include "Weave.h"
#include "SetupData.h"
#include "ComputeResults.h"

#include <lal/UserInputParse.h>
#include <lal/LFTandTSutils.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Extra toplist output fields
///
enum tagWeaveToplistExtraOutputs {
  WEAVE_TOPLIST_EXTRA_OUTPUTS
};

WeaveOutputResults *XLALWeaveOutputResultsCreate(
  const LIGOTimeGPS *ref_time,
  const size_t nspins,
  WeaveStatisticsParams *statistics_params,
  const UINT4 toplist_limit,
  const BOOLEAN mean2F_hgrm
  );
void XLALWeaveOutputResultsDestroy(
  WeaveOutputResults *out
  );
int XLALWeaveOutputResultsAdd(
  WeaveOutputResults *out,
  const WeaveSemiResults *semi_res,
  const UINT4 semi_nfreqs
  );
int XLALWeaveOutputResultsCompletionLoop(
  WeaveOutputResults *out
  );
int XLALWeaveOutputResultsWrite(
  FITSFile *file,
  const WeaveOutputResults *out
  );
int XLALWeaveOutputResultsReadAppend(
  FITSFile *file,
  WeaveOutputResults **out,
  UINT4 toplist_limit
  );
int XLALWeaveOutputResultsCompare(
  BOOLEAN *equal,
  const WeaveSetupData *setup,
  const BOOLEAN sort_by_semi_phys,
  const UINT4 round_param_to_dp,
  const UINT4 round_param_to_sf,
  const REAL8 unmatched_item_tol,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const UINT4 toplist_compare_limit,
  const WeaveOutputResults *out_1,
  const WeaveOutputResults *out_2
  );

#ifdef __cplusplus
}
#endif

#endif // _OUTPUT_RESULTS_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
