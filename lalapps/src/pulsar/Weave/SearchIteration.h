//
// Copyright (C) 2017 Karl Wette
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

#ifndef _SEARCH_ITERATION_H
#define _SEARCH_ITERATION_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
/// \brief Module which implements iterators over search parameter spaces
///

#include "Weave.h"

#include <lal/LatticeTiling.h>

#ifdef __cplusplus
extern "C" {
#endif

WeaveSearchIterator *XLALWeaveMainLoopSearchIteratorCreate(
  const LatticeTiling *semi_tiling,
  const UINT4 freq_partitions,
  const UINT4 f1dot_partitions
  );
void XLALWeaveSearchIteratorDestroy(
  WeaveSearchIterator *itr
  );
int XLALWeaveSearchIteratorSave(
  const WeaveSearchIterator *itr,
  FITSFile *file
  );
int XLALWeaveSearchIteratorRestore(
  WeaveSearchIterator *itr,
  FITSFile *file
  );
int XLALWeaveSearchIteratorNext(
  WeaveSearchIterator *itr,
  BOOLEAN *iteration_complete,
  BOOLEAN *expire_cache,
  UINT8 *semi_index,
  const gsl_vector **semi_rssky,
  INT4 *semi_left,
  INT4 *semi_right,
  UINT4 *repetition_index
  );
REAL8 XLALWeaveSearchIteratorProgress(
  const WeaveSearchIterator *itr
  );
REAL8 XLALWeaveSearchIteratorRemainingTime(
  const WeaveSearchIterator *itr,
  const REAL8 elapsed_time
  );

#ifdef __cplusplus
}
#endif

#endif // _SEARCH_ITERATION_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
