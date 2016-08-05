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

#ifndef _WEAVE_H
#define _WEAVE_H

///
/// \defgroup lalapps_pulsar_Weave Weave Search Application
/// \ingroup lalapps_pulsar_Apps
/// \author Karl Wette
///

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>

#include <lal/LALStdlib.h>
#include <lal/ComputeFstat.h>
#include <lal/FITSFileIO.h>
#include <lal/LALBarycenter.h>
#include <lal/Segments.h>
#include <lal/LatticeTiling.h>
#include <lal/SuperskyMetrics.h>
#include <lal/VectorMath.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#ifdef __cplusplus
extern "C" {
#endif

///
/// Compare two quantities, and return a sort order value if they are unequal
///
#define WEAVE_COMPARE_BY( x, y ) do { \
    if ( (x) < (y) ) return -1; \
    if ( (x) > (y) ) return +1; \
  } while(0)

///
/// Function which transforms a point from physical coordinates to lattice tiling coordinates
///
typedef int ( *WeavePhysicalToLattice )( gsl_vector *out_latt, const PulsarDopplerParams *in_phys, const void *transf_data );

///
/// Function which transforms a point from lattice tiling coordinates to physical coordinates
///
typedef int ( *WeaveLatticeToPhysical )( PulsarDopplerParams *out_phys, const gsl_vector *in_latt, const void *transf_data );

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
/// Cache used to store coherent results
///
typedef struct tagWeaveCache WeaveCache;

///
/// Storage for a series of cache queries
///
typedef struct tagWeaveCacheQueries WeaveCacheQueries;

///
/// Output data from a search
///
typedef struct tagWeaveOutput WeaveOutput;

///
/// Setup data which is computed only once for a given search setup
///
typedef struct {
  /// Physical to lattice coordinate transform
  WeavePhysicalToLattice phys_to_latt;
  /// Lattice to physical coordinate transform
  WeaveLatticeToPhysical latt_to_phys;
  /// Reference time at which search is conducted
  LIGOTimeGPS ref_time;
  /// List of detector names for which metrics were computed
  LALStringVector *detectors;
  /// Segment list for which metrics were computed
  LALSegList *segments;
  /// Reduced supersky parameter-space metrics
  SuperskyMetrics *metrics;
  /// Ephemeris data over time-span of segments
  EphemerisData *ephemerides;
} WeaveSetup;

///
/// Output toplist per-segment item
///
typedef struct tagWeaveOutputToplistPerSegItem {
  /// Physical coordinates of coherent template
  PulsarDopplerParams coh_phys;
  /// Coherent multi-detector F-statistic
  REAL4 twoF;
  /// Coherent per-detector F-statistic
  REAL4 twoF_per_det[PULSAR_MAX_DETECTORS];
} WeaveOutputToplistPerSegItem;

///
/// Output toplist item
///
typedef struct tagWeaveOutputToplistItem {
  /// Per-segment items (optional)
  WeaveOutputToplistPerSegItem *per_seg;
  /// Physical coordinates of semicoherent template
  PulsarDopplerParams semi_phys;
  /// Mean multi-detector F-statistic
  REAL4 mean_twoF;
  /// Mean per-detector F-statistic
  REAL4 mean_twoF_per_det[PULSAR_MAX_DETECTORS];
} WeaveOutputToplistItem;

///
/// Output various information for each segment
///
typedef struct tagWeaveOutputPerSegInfo {
  /// Start time of segment
  LIGOTimeGPS segment_start;
  /// End time of segment
  LIGOTimeGPS segment_end;
  /// Timestamp of first SFT from each detector
  LIGOTimeGPS sft_first[PULSAR_MAX_DETECTORS];
  /// Timestamp of last SFT from each detector
  LIGOTimeGPS sft_last[PULSAR_MAX_DETECTORS];
  /// Number of SFTs from each detector
  INT4 sft_count[PULSAR_MAX_DETECTORS];
  /// Minimum of frequency covering range
  REAL8 min_cover_freq;
  /// Maximum of frequency covering range
  REAL8 max_cover_freq;
  /// Total number of computed coherent results
  INT4 coh_total;
  /// Total number of recomputed coherent results
  INT4 coh_total_recomp;
} WeaveOutputPerSegInfo;

///
/// \name Routines which handle the setup data
///
/// @{

void XLALWeaveSetupClear(
  WeaveSetup *setup
  );
int XLALWeaveSetupWrite(
  FITSFile *file,
  const WeaveSetup *setup
  );
int XLALWeaveSetupRead(
  FITSFile *file,
  WeaveSetup *setup
  );

/// @}

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
int XLALWeaveFillOutputToplistItem(
  WeaveOutputToplistItem **item,
  BOOLEAN *full_init,
  const WeaveSemiResults *semi_res,
  const size_t freq_idx
  );

/// @}

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

///
/// \name Routines which handle the output data
///
/// @{

WeaveOutput *XLALWeaveOutputCreate(
  const LIGOTimeGPS *ref_time,
  const int toplist_limit,
  const size_t nspins,
  const LALStringVector *per_detectors,
  const UINT4 per_nsegments
  );
void XLALWeaveOutputDestroy(
  WeaveOutput *out
  );
int XLALWeaveOutputAdd(
  WeaveOutput *out,
  const WeaveSemiResults *semi_res,
  const UINT4 semi_nfreqs
  );
int XLALWeaveOutputWrite(
  FITSFile *file,
  const WeaveOutput *out
  );
int XLALWeaveOutputWriteExtra(
  FITSFile *file,
  const LALStringVector *detectors,
  const size_t nsegments,
  const WeaveOutputPerSegInfo *per_seg_info
  );

/// @}

#ifdef __cplusplus
}
#endif

#endif // _WEAVE_H
