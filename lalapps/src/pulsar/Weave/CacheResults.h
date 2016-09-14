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

#ifndef _CACHE_RESULTS_H
#define _CACHE_RESULTS_H

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "Weave.h"
#include "ComputeResults.h"

#include <lal/LatticeTiling.h>

#ifdef __cplusplus
extern "C" {
#endif

///
/// Cache used to store coherent results
///
typedef struct tagWeaveCache WeaveCache;

///
/// Storage for a series of cache queries
///
typedef struct tagWeaveCacheQueries WeaveCacheQueries;

///
/// \name Routines which cache computed coherent results
///
/// @{

WeaveCache *XLALWeaveCacheCreate(
  const LatticeTiling *coh_tiling,
  const BOOLEAN interpolation,
  const WeavePhysicalToLattice phys_to_latt,
  const WeaveLatticeToPhysical latt_to_phys,
  const void *coh_transf_data,
  const void *semi_transf_data,
  WeaveCohInput *coh_input,
  const size_t max_size,
  const size_t gc_limit,
  const BOOLEAN per_seg_info
  );
void XLALWeaveCacheDestroy(
  WeaveCache *cache
  );
WeaveCacheQueries *XLALWeaveCacheQueriesCreate(
  const LatticeTiling *semi_tiling,
  const WeavePhysicalToLattice phys_to_latt,
  const WeaveLatticeToPhysical latt_to_phys,
  const void *semi_transf_data,
  const UINT4 nqueries,
  const UINT4 npartitions
  );
void XLALWeaveCacheQueriesDestroy(
  WeaveCacheQueries *queries
  );
int XLALWeaveCacheQueriesInit(
  WeaveCacheQueries *queries,
  const LatticeTilingIterator *semi_itr,
  const gsl_vector *semi_point
  );
int XLALWeaveCacheQuery(
  const WeaveCache *cache,
  const UINT8 semi_index,
  WeaveCacheQueries *queries,
  const UINT4 query_index
  );
int XLALWeaveCacheQueriesFinal(
  WeaveCacheQueries *queries,
  const UINT4 partition_index,
  PulsarDopplerParams *semi_phys,
  const double dfreq,
  UINT4 *semi_nfreqs
  );
int XLALWeaveCacheRetrieve(
  WeaveCache *cache,
  WeaveCacheQueries *queries,
  const UINT4 query_index,
  const WeaveCohResults **coh_res,
  UINT4 *coh_offset,
  WeaveOutputPerSegInfo *per_seg_info
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _CACHE_RESULTS_H
