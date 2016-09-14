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

#ifndef _RESULTS_BASKET_H
#define _RESULTS_BASKET_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"
#include "SetupData.h"
#include "ComputeResults.h"
#include "OutputResults.h"

#include <lal/LALHeap.h>
#include <lal/LFTandTSutils.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Basket of output results
///
typedef struct tagWeaveResultsBasket WeaveResultsBasket;

///
/// \name Routines which handle baskets of results
///
/// @{

WeaveResultsBasket *XLALWeaveResultsBasketCreate(
  const WeaveOutputParams *par,
  const char *stat_name,
  const char *stat_desc,
  const int toplist_limit,
  LALHeapCmpFcn toplist_result_item_compare_fcn
  );
void XLALWeaveResultsBasketDestroy(
  WeaveResultsBasket *basket
  );
int XLALWeaveResultsBasketAdd(
  WeaveResultsBasket *basket,
  const WeaveSemiResults *semi_res,
  const UINT4 semi_nfreqs
  );
int XLALWeaveResultsBasketWrite(
  FITSFile *file,
  const WeaveResultsBasket *basket
  );
int XLALWeaveResultsBasketReadAppend(
  FITSFile *file,
  WeaveResultsBasket *basket
  );
int XLALWeaveResultsBasketCompare(
  BOOLEAN *equal,
  const WeaveSetupData *setup,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const WeaveResultsBasket *basket_1,
  const WeaveResultsBasket *basket_2
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _RESULTS_BASKET_H
