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
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

#ifndef _ITERATION_H
#define _ITERATION_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"

#include <lal/LatticeTiling.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// \name Functions which implement iterators over search parameter spaces
///
/// @{

///
/// Iterator over a search parameter space
///
typedef struct tagWeaveIterator WeaveIterator;

WeaveIterator *XLALWeaveMainLoopIteratorCreate(
  const LatticeTiling *semi_tiling,
  const UINT4 freq_partitions
  );
void XLALWeaveIteratorDestroy(
  WeaveIterator *itr
  );
int XLALWeaveIteratorSave(
  const WeaveIterator *itr,
  FITSFile *file
  );
int XLALWeaveIteratorRestore(
  WeaveIterator *itr,
  FITSFile *file
  );
int XLALWeaveIteratorNext(
  WeaveIterator *itr,
  BOOLEAN *iteration_complete,
  BOOLEAN *expire_cache,
  UINT8 *semi_index,
  const gsl_vector **semi_rssky,
  INT4 *semi_left,
  INT4 *semi_right,
  UINT4 *repetition_index
  );
REAL8 XLALWeaveIteratorProgress(
  const WeaveIterator *itr
  );
REAL8 XLALWeaveIteratorRemainingTime(
  const WeaveIterator *itr,
  const REAL8 elapsed_time
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _ITERATION_H

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
