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

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "ComputeResults.h"

#include <lal/UserInputPrint.h>
#include <lal/ExtrapolatePulsarSpins.h>

///
/// Aligned arrays use maximum required alignment, i.e.\ 32 bytes for AVX
///
const UINT4 alignment = 32;

///
/// Input data segment info
///
typedef struct {
  /// Start time of segment
  LIGOTimeGPS segment_start;
  /// End time of segment
  LIGOTimeGPS segment_end;
  /// Timestamp of first SFT from each detector
  LIGOTimeGPS sft_first[PULSAR_MAX_DETECTORS];
  /// Timestamp of last SFT from each detector
  LIGOTimeGPS sft_last[PULSAR_MAX_DETECTORS];
  /// Number of SFTs from each detector
  UINT4 sft_count[PULSAR_MAX_DETECTORS];
  /// Minimum of frequency range loaded from SFTs
  REAL8 sft_min_freq;
  /// Maximum of frequency range loaded from SFTs
  REAL8 sft_max_freq;
} segment_info;

///
/// Input data required for computing coherent results
///
struct tagWeaveCohInput {
  /// List of detector names from setup file
  const LALStringVector *setup_detectors;
  /// Bitflag representing search simulation level
  WeaveSimulationLevel simulation_level;
  /// Input data segment info
  segment_info seg_info;
  /// Whether input data segment info contains SFT info
  BOOLEAN seg_info_have_sft_info;
  /// F-statistic input data
  FstatInput *Fstat_input;
  /// What F-statistic quantities to compute
  FstatQuantities Fstat_what_to_compute;
  /// Number of detectors in F-statistic data
  UINT4 Fstat_ndetectors;
  /// Map detectors in F-statistic data in this segment to 'global' detector list across segments in coherent results
  size_t Fstat_res_idx[PULSAR_MAX_DETECTORS];
  /// Whether F-statistic timing info is being collected
  BOOLEAN Fstat_collect_timing;
};

///
/// Results of a coherent computation on a single segment
///
struct tagWeaveCohResults {
  /// Coherent template parameters of the first frequency bin
  PulsarDopplerParams coh_phys;
  /// Number of frequencies
  UINT4 nfreqs;
  /// Multi-detector F-statistics per frequency
  REAL4Vector *coh2F;
  /// Per-detector F-statistics per frequency
  REAL4Vector *coh2F_det[PULSAR_MAX_DETECTORS];
};

///
/// \name Internal functions
///
/// @{

static int semi_res_sum_2F( UINT4 *nsum, REAL4 *sum2F, const REAL4 *coh2F, const UINT4 nfreqs );
static int semi_res_max_2F( UINT4 *nmax, REAL4 *max2F, const REAL4 *coh2F, const UINT4 nfreqs );

/// @}

///
/// Create coherent input data
///
WeaveCohInput *XLALWeaveCohInputCreate(
  const LALStringVector *setup_detectors,
  const WeaveSimulationLevel simulation_level,
  const SFTCatalog *sft_catalog,
  const UINT4 segment_index,
  const LALSeg *segment,
  const PulsarDopplerParams *min_phys,
  const PulsarDopplerParams *max_phys,
  const double dfreq,
  const EphemerisData *ephemerides,
  const LALStringVector *sft_noise_sqrtSX,
  const LALStringVector *Fstat_assume_sqrtSX,
  FstatOptionalArgs *Fstat_opt_args,
  const WeaveStatisticsParams *statistics_params,
  BOOLEAN recalc_stage
  )
{

  // Check input
  XLAL_CHECK_NULL( setup_detectors != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( ( simulation_level & WEAVE_SIMULATE_MIN_MEM ) || ( sft_catalog != NULL ), XLAL_EFAULT );
  XLAL_CHECK_NULL( segment != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( min_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( max_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( dfreq >= 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( ephemerides != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( Fstat_opt_args != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( statistics_params != NULL, XLAL_EFAULT );

  // Allocate memory
  WeaveCohInput *coh_input = XLALCalloc( 1, sizeof( *coh_input ) );
  XLAL_CHECK_NULL( coh_input != NULL, XLAL_ENOMEM );

  // Set fields
  coh_input->setup_detectors = setup_detectors;
  coh_input->simulation_level = simulation_level;
  coh_input->seg_info_have_sft_info = ( sft_catalog != NULL );
  coh_input->Fstat_collect_timing = Fstat_opt_args->collectTiming;

  // Record information from segment
  coh_input->seg_info.segment_start = segment->start;
  coh_input->seg_info.segment_end = segment->end;

  // Decide what F-statistic quantities to compute
  WeaveStatisticType requested_stats = ( recalc_stage ) ? statistics_params->completionloop_statistics[1] : statistics_params->mainloop_statistics;
  if ( requested_stats & WEAVE_STATISTIC_COH2F ) {
    coh_input->Fstat_what_to_compute |= FSTATQ_2F;
  }
  if ( requested_stats & WEAVE_STATISTIC_COH2F_DET ) {
    coh_input->Fstat_what_to_compute |= FSTATQ_2F_PER_DET;
  }

  // Return now if simulating search with minimal memory allocation
  if ( coh_input->simulation_level & WEAVE_SIMULATE_MIN_MEM ) {
    return coh_input;
  }

  // Get a timeslice of SFT catalog restricted to the given segment
  SFTCatalog XLAL_INIT_DECL( sft_catalog_seg );
  XLAL_CHECK_NULL( XLALSFTCatalogTimeslice( &sft_catalog_seg, sft_catalog, &segment->start, &segment->end ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK_NULL( sft_catalog_seg.length > 0, XLAL_EINVAL, "No SFTs found in segment %u", segment_index );

  // Check that the number of SFTs in each segment matches the number provided by the segment list in the setup file, if nonzero
  const UINT4 sft_count = ( UINT4 ) segment->id;
  XLAL_CHECK_NULL( sft_count == 0 || sft_catalog_seg.length == sft_count, XLAL_EFAILED, "Number of SFTs found for segment %u (%u) is inconsistent with expected number of SFTs given by segment list (%u)", segment_index, sft_catalog_seg.length, sft_count );

  // Get list of detectors of SFT catalog in this segment
  LALStringVector *sft_catalog_seg_detectors = XLALListIFOsInCatalog( &sft_catalog_seg );
  XLAL_CHECK_NULL( sft_catalog_seg_detectors != NULL, XLAL_EFUNC );

  // Compute frequency range covered by spindown range over in the given segment
  LIGOTimeGPS sft_start = sft_catalog_seg.data[0].header.epoch;
  LIGOTimeGPS sft_end = sft_catalog_seg.data[sft_catalog_seg.length - 1].header.epoch;
  const double sft_end_timebase = 1.0 / sft_catalog_seg.data[sft_catalog_seg.length - 1].header.deltaF;
  XLALGPSAdd( &sft_end, sft_end_timebase );
  PulsarSpinRange XLAL_INIT_DECL( spin_range );
  XLAL_CHECK_NULL( XLALInitPulsarSpinRangeFromSpins( &spin_range, &min_phys->refTime, min_phys->fkdot, max_phys->fkdot ) == XLAL_SUCCESS, XLAL_EFUNC );
  double sft_min_cover_freq = 0, sft_max_cover_freq = 0;
  XLAL_CHECK_NULL( XLALCWSignalCoveringBand( &sft_min_cover_freq, &sft_max_cover_freq, &sft_start, &sft_end, &spin_range, 0, 0, 0 ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Parse SFT noise sqrt(Sh) string vector for detectors in this segment
  // - This is important when segments contain data from a subset of detectors
  MultiNoiseFloor Fstat_injectSqrtSX;
  if ( sft_noise_sqrtSX != NULL ) {
    XLAL_CHECK_NULL( XLALParseMultiNoiseFloorMapped( &Fstat_injectSqrtSX, sft_catalog_seg_detectors, sft_noise_sqrtSX, setup_detectors ) == XLAL_SUCCESS, XLAL_EFUNC );
    Fstat_opt_args->injectSqrtSX = &Fstat_injectSqrtSX;
  }

  // Parse F-statistic assumed sqrt(Sh) string vector for detectors in this segment
  // - This is important when segments contain data from a subset of detectors
  MultiNoiseFloor Fstat_assumeSqrtSX;
  if ( Fstat_assume_sqrtSX != NULL ) {
    XLAL_CHECK_NULL( XLALParseMultiNoiseFloorMapped( &Fstat_assumeSqrtSX, sft_catalog_seg_detectors, Fstat_assume_sqrtSX, setup_detectors ) == XLAL_SUCCESS, XLAL_EFUNC );
    Fstat_opt_args->assumeSqrtSX = &Fstat_assumeSqrtSX;
  }

  // Load F-statistic input data
  coh_input->Fstat_input = XLALCreateFstatInput( &sft_catalog_seg, sft_min_cover_freq, sft_max_cover_freq, dfreq, ephemerides, Fstat_opt_args );
  XLAL_CHECK_NULL( coh_input->Fstat_input != NULL, XLAL_EFUNC );
  Fstat_opt_args->prevInput = coh_input->Fstat_input;

  // Map detectors in F-statistic data in the given segment to their index in the coherent results
  // - This is important when segments contain data from a subset of detectors
  // - Map entry 'i' in 'Fstat_detector_info' (F-statistic data) to entry 'idx' in 'detectors' (coherent results)
  if ( coh_input->Fstat_what_to_compute & FSTATQ_2F_PER_DET ) {
    const MultiLALDetector *Fstat_detector_info = XLALGetFstatInputDetectors( coh_input->Fstat_input );
    coh_input->Fstat_ndetectors = Fstat_detector_info->length;
    char *statistics_detectors_string = XLALConcatStringVector( statistics_params->detectors, "," );
    for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
      const char *prefix = Fstat_detector_info->sites[i].frDetector.prefix;
      const int idx = XLALFindStringInVector( prefix, statistics_params->detectors );
      XLAL_CHECK_NULL( idx >= 0, XLAL_EFAILED, "Detector '%s' from F-statistic data not found in list of detectors '%s'", prefix, statistics_detectors_string );
      coh_input->Fstat_res_idx[i] = idx;
    }
    XLALFree( statistics_detectors_string );
  }

  // Record information from SFTs in the given segment
  {
    MultiSFTCatalogView *sft_catalog_seg_view = XLALGetMultiSFTCatalogView( &sft_catalog_seg );
    XLAL_CHECK_NULL( sft_catalog_seg_view != NULL, XLAL_EINVAL );
    for ( size_t j = 0; j < sft_catalog_seg_view->length; ++j ) {
      XLAL_CHECK_NULL( sft_catalog_seg_view->data[j].length > 0, XLAL_EINVAL );
      char *detector_name = XLALGetChannelPrefix( sft_catalog_seg_view->data[j].data[0].header.name );
      XLAL_CHECK_NULL( detector_name != NULL, XLAL_EFUNC );
      const int k = XLALFindStringInVector( detector_name, sft_catalog_seg_detectors );
      if ( k >= 0 ) {
        const UINT4 length = sft_catalog_seg_view->data[j].length;
        coh_input->seg_info.sft_first[k] = sft_catalog_seg_view->data[j].data[0].header.epoch;
        coh_input->seg_info.sft_last[k] = sft_catalog_seg_view->data[j].data[length - 1].header.epoch;
        coh_input->seg_info.sft_count[k] = length;
      }
      XLALFree( detector_name );
    }
    XLALDestroyMultiSFTCatalogView( sft_catalog_seg_view );
  }
  XLAL_CHECK_NULL( XLALGetFstatInputSFTBand( coh_input->Fstat_input, &coh_input->seg_info.sft_min_freq, &coh_input->seg_info.sft_max_freq ) == XLAL_SUCCESS, XLAL_EFUNC );


  // Cleanup
  XLALDestroyStringVector( sft_catalog_seg_detectors );

  return coh_input;

}

///
/// Destroy coherent input data
///
void XLALWeaveCohInputDestroy(
  WeaveCohInput *coh_input
  )
{
  if ( coh_input != NULL ) {
    XLALDestroyFstatInput( coh_input->Fstat_input );
    XLALFree( coh_input );
  }
}

///
/// Write various information from coherent input data to a FITS file
///
int XLALWeaveCohInputWriteInfo(
  FITSFile *file,
  const size_t ncoh_input,
  WeaveCohInput *const *coh_input
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( ncoh_input > 0, XLAL_ESIZE );
  XLAL_CHECK( coh_input != NULL, XLAL_EFAULT );

  // Write total number of SFTs used by search
  {
    UINT4 nsfts = 0;
    for ( size_t i = 0; i < ncoh_input; ++i ) {
      for ( size_t j = 0; j < PULSAR_MAX_DETECTORS; ++j ) {
        nsfts += coh_input[i]->seg_info.sft_count[j];
      }
    }
    XLAL_CHECK_MAIN( XLALFITSHeaderWriteUINT4( file, "nsfts", nsfts, "number of SFTs used by search" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Return now if simulating search
  if ( coh_input[0]->simulation_level & WEAVE_SIMULATE ) {
    return XLAL_SUCCESS;
  }

  // Write F-statistic method name
  {
    const char *method_name = XLALGetFstatInputMethodName( coh_input[0]->Fstat_input );
    XLAL_CHECK( method_name != NULL, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteString( file, "fstat method", method_name, "name of F-statistic method" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write F-statistic timing information
  if ( coh_input[0]->Fstat_collect_timing ) {

    // Get timing information from F-statistic input data
    REAL4 tauF_eff = 0, tauF_core = 0, tauF_buffer = 0, NCalls = 0, NBufferMisses = 0;
    UINT4 nmodel = 0;
    const char *XLAL_INIT_DECL( model_names, [TIMING_MODEL_MAX_VARS] );
    REAL4 XLAL_INIT_DECL( model_values, [TIMING_MODEL_MAX_VARS] );
    for ( size_t i = 0; i < ncoh_input; ++i ) {

      // Get timing data
      FstatTimingGeneric XLAL_INIT_DECL( timing_generic );
      FstatTimingModel XLAL_INIT_DECL( timing_model );
      XLAL_CHECK( XLALGetFstatTiming( coh_input[i]->Fstat_input, &timing_generic, &timing_model ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( timing_generic.NCalls > 0, XLAL_EFAILED );

      // Accumulate generic timing constants
      tauF_eff += timing_generic.tauF_eff;
      tauF_core += timing_generic.tauF_core;
      tauF_buffer += timing_generic.tauF_buffer;
      NCalls += timing_generic.NCalls;
      NBufferMisses += timing_generic.NBufferMisses;

      // Get names of method-specific timing constants
      if ( nmodel == 0 ) {
        nmodel = timing_model.numVariables;
        for ( size_t j = 0; j < nmodel; ++j ) {
          model_names[j] = timing_model.names[j];
        }
      } else {
        XLAL_CHECK( nmodel == timing_model.numVariables, XLAL_EFAILED );
        for ( size_t j = 0; j < nmodel; ++j ) {
          XLAL_CHECK( strcmp( model_names[j], timing_model.names[j] ) == 0, XLAL_EFAILED );
        }
      }

      // Accumulate method-specific timing constants
      for ( size_t j = 0; j < nmodel; ++j ) {
        model_values[j] += timing_model.values[j];
      }

    }

    // Write generic timing constants
    XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, "fstat tauF_eff", tauF_eff / ncoh_input, "F-statistic generic timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, "fstat tauF_core", tauF_core / ncoh_input, "F-statistic generic timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, "fstat tauF_buffer", tauF_buffer / ncoh_input, "F-statistic generic timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, "fstat b", NBufferMisses / NCalls, "F-statistic generic timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Write method-specific timing constants
    for ( size_t j = 0; j < nmodel; ++j ) {
      char keyword[64];
      snprintf( keyword, sizeof( keyword ), "fstat %s", model_names[j] );
      XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, keyword, model_values[j] / ncoh_input, "F-statistic method-specific timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  }

  return XLAL_SUCCESS;

}

///
/// Write various segment information from coherent input data to a FITS file
///
int XLALWeaveCohInputWriteSegInfo(
  FITSFile *file,
  const size_t ncoh_input,
  WeaveCohInput *const *coh_input
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( ncoh_input > 0, XLAL_ESIZE );
  XLAL_CHECK( coh_input != NULL, XLAL_EFAULT );

  // Begin FITS table
  XLAL_CHECK( XLALFITSTableOpenWrite( file, "segment_info", "segment information" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Describe FITS table
  char col_name[64];
  XLAL_FITS_TABLE_COLUMN_BEGIN( segment_info );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, segment_start ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, segment_end ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( coh_input[0]->seg_info_have_sft_info ) {
    for ( size_t i = 0; i < coh_input[0]->setup_detectors->length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "sft_first_%s", coh_input[0]->setup_detectors->data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, GPSTime, sft_first[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      snprintf( col_name, sizeof( col_name ), "sft_last_%s", coh_input[0]->setup_detectors->data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, GPSTime, sft_last[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      snprintf( col_name, sizeof( col_name ), "sft_count_%s", coh_input[0]->setup_detectors->data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, UINT4, sft_count[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, sft_min_freq ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, sft_max_freq ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write FITS table
  for ( size_t i = 0; i < ncoh_input; ++i ) {
    XLAL_CHECK( XLALFITSTableWriteRow( file, &coh_input[i]->seg_info ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}

///
/// Create and compute coherent results
///
int XLALWeaveCohResultsCompute(
  WeaveCohResults **coh_res,
  WeaveCohInput *coh_input,
  const PulsarDopplerParams *coh_phys,
  const UINT4 coh_nfreqs,
  WeaveSearchTiming *tim
  )
{

  // Check input
  XLAL_CHECK( coh_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( coh_input != NULL, XLAL_EFAULT );
  XLAL_CHECK( coh_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( coh_nfreqs > 0, XLAL_EINVAL );

  // Allocate results struct if required
  if ( *coh_res == NULL ) {
    *coh_res = XLALCalloc( 1, sizeof( **coh_res ) );
    XLAL_CHECK( *coh_res != NULL, XLAL_ENOMEM );
  }

  // Set fields
  ( *coh_res )->coh_phys = *coh_phys;
  ( *coh_res )->nfreqs = coh_nfreqs;

  // Return now if simulating search with minimal memory allocation
  if ( coh_input->simulation_level & WEAVE_SIMULATE_MIN_MEM ) {
    return XLAL_SUCCESS;
  }

  // Reallocate vectors of multi- and per-detector F-statistics per frequency
  if ( ( *coh_res )->coh2F == NULL || ( *coh_res )->coh2F->length < ( *coh_res )->nfreqs ) {
    ( *coh_res )->coh2F = XLALResizeREAL4Vector( ( *coh_res )->coh2F, ( *coh_res )->nfreqs );
    XLAL_CHECK( ( *coh_res )->coh2F != NULL, XLAL_ENOMEM );
  }
  if ( coh_input->Fstat_what_to_compute & FSTATQ_2F_PER_DET ) {
    for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
      const size_t idx = coh_input->Fstat_res_idx[i];
      if ( ( *coh_res )->coh2F_det[idx] == NULL || ( *coh_res )->coh2F_det[idx]->length < ( *coh_res )->nfreqs ) {
        ( *coh_res )->coh2F_det[idx] = XLALResizeREAL4Vector( ( *coh_res )->coh2F_det[idx], ( *coh_res )->nfreqs );
        XLAL_CHECK( ( *coh_res )->coh2F_det[idx] != NULL, XLAL_ENOMEM );
      }
    }
  }

  // Return now if simulating search
  if ( coh_input->simulation_level & WEAVE_SIMULATE ) {
    return XLAL_SUCCESS;
  }

  // Start timing of coherent results
  if ( tim != NULL ) {
    XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_NONE, WEAVE_STATISTIC_COH2F ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Use a local F-statistic results structure, since we supply our own memory
  // - The 'internalalloclen' field stores the memory size in elements of the results arrays (e.g. 'coh2F'),
  //   as opposed to 'numFreqBins' which stores how many elements of the result arrays are in use.
  FstatResults XLAL_INIT_DECL( Fstat_res_struct );
  FstatResults *Fstat_res = &Fstat_res_struct;
  Fstat_res->internalalloclen = ( *coh_res )->nfreqs;
  if ( coh_input->Fstat_what_to_compute & FSTATQ_2F_PER_DET ) {
    Fstat_res->numDetectors = coh_input->Fstat_ndetectors;
  }
  Fstat_res->twoF = ( *coh_res )->coh2F->data;
  for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
    const size_t idx = coh_input->Fstat_res_idx[i];
    Fstat_res->twoFPerDet[i] = ( *coh_res )->coh2F_det[idx]->data;
  }

  // Compute the F-statistic starting at the point 'coh_phys', with 'nfreqs' frequency bins
  XLAL_CHECK( XLALComputeFstat( &Fstat_res, coh_input->Fstat_input, coh_phys, ( *coh_res )->nfreqs, coh_input->Fstat_what_to_compute ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Sanity check the F-statistic results structure
  XLAL_CHECK( Fstat_res->internalalloclen == ( *coh_res )->nfreqs, XLAL_EFAILED );
  XLAL_CHECK( Fstat_res->numFreqBins == ( *coh_res )->nfreqs, XLAL_EFAILED );
  XLAL_CHECK( !( coh_input->Fstat_what_to_compute & FSTATQ_2F_PER_DET ) || ( Fstat_res->numDetectors == coh_input->Fstat_ndetectors ), XLAL_EFAILED );
  XLAL_CHECK( Fstat_res->twoF == ( *coh_res )->coh2F->data, XLAL_EFAILED );
  for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
    const size_t idx = coh_input->Fstat_res_idx[i];
    XLAL_CHECK( Fstat_res->twoFPerDet[i] == ( *coh_res )->coh2F_det[idx]->data, XLAL_EFAILED, "%p vs %p", Fstat_res->twoFPerDet[i], ( *coh_res )->coh2F_det[idx]->data );
  }

  // Stop timing of coherent results
  if ( tim != NULL ) {
    XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_COH2F, WEAVE_STATISTIC_NONE ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}

///
/// Destroy coherent results
///
void XLALWeaveCohResultsDestroy(
  WeaveCohResults *coh_res
  )
{
  if ( coh_res != NULL ) {
    XLALDestroyREAL4Vector( coh_res->coh2F );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALDestroyREAL4Vector( coh_res->coh2F_det[i] );
    }
    XLALFree( coh_res );
  }
}

///
/// Create and initialise semicoherent results
///
int XLALWeaveSemiResultsInit(
  WeaveSemiResults **semi_res,
  const WeaveSimulationLevel simulation_level,
  const UINT4 ndetectors,
  const UINT4 nsegments,
  const UINT8 semi_index,
  const PulsarDopplerParams *semi_phys,
  const double dfreq,
  const UINT4 semi_nfreqs,
  const WeaveStatisticsParams *statistics_params
  )
{

  // Check input
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( dfreq >= 0, XLAL_EINVAL );
  XLAL_CHECK( semi_nfreqs > 0, XLAL_EINVAL );
  XLAL_CHECK( statistics_params != NULL, XLAL_EFAULT );

  WeaveStatisticType mainloop_stats = statistics_params->mainloop_statistics;

  const WeaveStatisticType supported_mainloop = (
    0
    | WEAVE_STATISTIC_COH2F
    | WEAVE_STATISTIC_COH2F_DET
    | WEAVE_STATISTIC_MAX2F
    | WEAVE_STATISTIC_MAX2F_DET
    | WEAVE_STATISTIC_SUM2F
    | WEAVE_STATISTIC_SUM2F_DET
    | WEAVE_STATISTIC_MEAN2F
    | WEAVE_STATISTIC_BSGL
    | WEAVE_STATISTIC_BSGLtL
    | WEAVE_STATISTIC_BtSGLtL
    );

  WeaveStatisticType unsupported = ( mainloop_stats & ~supported_mainloop );
  if ( unsupported != 0 ) {
    char *unsupported_names = XLALPrintStringValueOfUserFlag( ( const int * )&unsupported, &WeaveStatisticChoices );
    XLALPrintError( "BUG: unsupported main-loop statistics requested: %s\n", unsupported_names );
    XLALFree( unsupported_names );
    XLAL_ERROR( XLAL_EERR );
  }

  // Allocate results struct if required
  if ( *semi_res == NULL ) {
    *semi_res = XLALCalloc( 1, sizeof( **semi_res ) );
    XLAL_CHECK( *semi_res != NULL, XLAL_ENOMEM );

    // Allocate array for per-segment index
    ( *semi_res )->coh_index = XLALCalloc( nsegments, sizeof( *( *semi_res )->coh_index ) );
    XLAL_CHECK( ( *semi_res )->coh_index != NULL, XLAL_ENOMEM );

    // Allocate array for per-segment coordinates
    ( *semi_res )->coh_phys = XLALCalloc( nsegments, sizeof( *( *semi_res )->coh_phys ) );
    XLAL_CHECK( ( *semi_res )->coh_phys != NULL, XLAL_ENOMEM );

    // If we need coh2F in "main loop": allocate vector
    if ( mainloop_stats & WEAVE_STATISTIC_COH2F ) {
      ( *semi_res )->coh2F = XLALCalloc( nsegments, sizeof( *( *semi_res )->coh2F ) );
      XLAL_CHECK( ( *semi_res )->coh2F != NULL, XLAL_ENOMEM );
    }

    // If we need coh2F_det in "main loop": allocate vector
    if ( mainloop_stats & WEAVE_STATISTIC_COH2F_DET ) {
      for ( size_t i = 0; i < ndetectors; ++i ) {
        ( *semi_res )->coh2F_det[i] = XLALCalloc( nsegments, sizeof( *( *semi_res )->coh2F_det[i] ) );
        XLAL_CHECK( ( *semi_res )->coh2F_det[i] != NULL, XLAL_ENOMEM );
      }
    }

  }

  // Set fields
  ( *semi_res )->simulation_level = simulation_level;
  ( *semi_res )->ndetectors = ndetectors;
  ( *semi_res )->nsegments = nsegments;
  ( *semi_res )->dfreq = dfreq;
  ( *semi_res )->nfreqs = semi_nfreqs;
  ( *semi_res )->semi_index = semi_index;
  ( *semi_res )->semi_phys = *semi_phys;
  ( *semi_res )->statistics_params = statistics_params;

  // Initialise number of coherent results
  ( *semi_res )->ncoh_res = 0;

  // Initialise number of max-over-segments to multi- and per-detector F-statistics
  ( *semi_res )->nmax2F = 0;
  XLAL_INIT_MEM( ( *semi_res )->nmax2F_det );
  // If we need max2F in "main loop": Reallocate vector of max-over-segments of multi-detector F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_MAX2F ) {
    if ( ( *semi_res )->max2F == NULL || ( *semi_res )->max2F->length < ( *semi_res )->nfreqs ) {
      ( *semi_res )->max2F = XLALResizeREAL4VectorAligned( ( *semi_res )->max2F, ( *semi_res )->nfreqs, alignment );
      XLAL_CHECK( ( *semi_res )->max2F != NULL, XLAL_ENOMEM );
    }
  }
  // If we need max2F_det in "main loop": Reallocate vector of max-over-segments of per-detector F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_MAX2F_DET ) {
    for ( size_t i = 0; i < ( *semi_res )->ndetectors; ++i ) {
      if ( ( *semi_res )->max2F_det[i] == NULL || ( *semi_res )->max2F_det[i]->length < ( *semi_res )->nfreqs ) {
        ( *semi_res )->max2F_det[i] = XLALResizeREAL4VectorAligned( ( *semi_res )->max2F_det[i], ( *semi_res )->nfreqs, alignment );
        XLAL_CHECK( ( *semi_res )->max2F_det[i] != NULL, XLAL_ENOMEM );
      }
    }
  }

  // Initialise number of additions to multi- and per-detector F-statistics
  ( *semi_res )->nsum2F = 0;
  XLAL_INIT_MEM( ( *semi_res )->nsum2F_det );
  // If we need sum2F in "main loop": Reallocate vector of sum of multi-detector F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_SUM2F ) {
    if ( ( *semi_res )->sum2F == NULL || ( *semi_res )->sum2F->length < ( *semi_res )->nfreqs ) {
      ( *semi_res )->sum2F = XLALResizeREAL4VectorAligned( ( *semi_res )->sum2F, ( *semi_res )->nfreqs, alignment );
      XLAL_CHECK( ( *semi_res )->sum2F != NULL, XLAL_ENOMEM );
    }
  }
  // If we need sum2F_det in "main loop": Reallocate vector of sum of per-detector F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_SUM2F_DET ) {
    for ( size_t i = 0; i < ( *semi_res )->ndetectors; ++i ) {
      if ( ( *semi_res )->sum2F_det[i] == NULL || ( *semi_res )->sum2F_det[i]->length < ( *semi_res )->nfreqs ) {
        ( *semi_res )->sum2F_det[i] = XLALResizeREAL4VectorAligned( ( *semi_res )->sum2F_det[i], ( *semi_res )->nfreqs, alignment );
        XLAL_CHECK( ( *semi_res )->sum2F_det[i] != NULL, XLAL_ENOMEM );
      }
    }
  }

  // If we compute mean2F in "main loop": Reallocate vector of mean multi-F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_MEAN2F ) {
    if ( ( *semi_res )->mean2F == NULL || ( *semi_res )->mean2F->length < ( *semi_res )->nfreqs ) {
      ( *semi_res )->mean2F = XLALResizeREAL4VectorAligned( ( *semi_res )->mean2F, ( *semi_res )->nfreqs, alignment );
      XLAL_CHECK( ( *semi_res )->mean2F != NULL, XLAL_ENOMEM );
    }
  }

  // (Re-)allocate vectors of line-robust log10(B_S/GL) statistic IFF used as a toplist statistic
  if ( mainloop_stats & WEAVE_STATISTIC_BSGL ) {
    if ( ( *semi_res )->log10BSGL == NULL || ( *semi_res )->log10BSGL->length < ( *semi_res )->nfreqs ) {
      ( *semi_res )->log10BSGL = XLALResizeREAL4VectorAligned( ( *semi_res )->log10BSGL, ( *semi_res )->nfreqs, alignment );
      XLAL_CHECK( ( *semi_res )->log10BSGL != NULL, XLAL_ENOMEM );
    }
  }

  // (Re-)allocate vectors of transient-line-robust log10(B_S/GLtL) statistic IFF used as a toplist statistic
  if ( mainloop_stats & WEAVE_STATISTIC_BSGLtL ) {
    if ( ( *semi_res )->log10BSGLtL == NULL || ( *semi_res )->log10BSGLtL->length < ( *semi_res )->nfreqs ) {
      ( *semi_res )->log10BSGLtL = XLALResizeREAL4VectorAligned( ( *semi_res )->log10BSGLtL, ( *semi_res )->nfreqs, alignment );
      XLAL_CHECK( ( *semi_res )->log10BSGLtL != NULL, XLAL_ENOMEM );
    }
  }

  // (Re-)allocate vectors of transient-signal line-robust log10(B_tS/GLtL) statistic IFF used as a toplist statistic
  if ( mainloop_stats & WEAVE_STATISTIC_BtSGLtL ) {
    if ( ( *semi_res )->log10BtSGLtL == NULL || ( *semi_res )->log10BtSGLtL->length < ( *semi_res )->nfreqs ) {
      ( *semi_res )->log10BtSGLtL = XLALResizeREAL4VectorAligned( ( *semi_res )->log10BtSGLtL, ( *semi_res )->nfreqs, alignment );
      XLAL_CHECK( ( *semi_res )->log10BtSGLtL != NULL, XLAL_ENOMEM );
    }
  }

  return XLAL_SUCCESS;

}

///
/// Add F-statistic array 'coh2F' to summed array 'sum2F', and keep track of the number of summations 'nsum'
///
static int semi_res_sum_2F(
  UINT4 *nsum,
  REAL4 *sum2F,
  const REAL4 *coh2F,
  const UINT4 nfreqs
  )
{
  if ( ( *nsum )++ == 0 ) {
    // If this is the first summation, just use memcpy()
    memcpy( sum2F, coh2F, sizeof( *sum2F ) * nfreqs );
    return XLAL_SUCCESS;
  }
  return XLALVectorAddREAL4( sum2F, sum2F, coh2F, nfreqs );
}

///
/// Track maximum between F-statistic array 'coh2F' and 'max2F', and keep track of the number of maximum-comparisons 'nmax'
///
static int semi_res_max_2F(
  UINT4 *nmax,
  REAL4 *max2F,
  const REAL4 *coh2F,
  const UINT4 nfreqs
  )
{
  if ( ( *nmax )++ == 0 ) {
    // If this is the first max-comparison, just use memcpy()
    memcpy( max2F, coh2F, sizeof( *max2F ) * nfreqs );
    return XLAL_SUCCESS;
  }
  return XLALVectorMaxREAL4( max2F, max2F, coh2F, nfreqs );
}

///
/// Add a new set of coherent results to the semicoherent results
///
int XLALWeaveSemiResultsAdd(
  WeaveSemiResults *semi_res,
  const WeaveCohResults *coh_res,
  const UINT8 coh_index,
  const UINT4 coh_offset,
  WeaveSearchTiming *tim
  )
{

  // Check input
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res->ncoh_res < semi_res->nsegments, XLAL_EINVAL );
  XLAL_CHECK( coh_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res->statistics_params != NULL, XLAL_EFAULT );
  XLAL_CHECK( tim != NULL, XLAL_EFAULT );

  // Check that offset does not overrun coherent results arrays
  XLAL_CHECK( coh_offset + semi_res->nfreqs <= coh_res->nfreqs, XLAL_EFAILED, "Coherent offset (%u) + number of semicoherent frequency bins (%u) > number of coherent frequency bins (%u)", coh_offset, semi_res->nfreqs, coh_res->nfreqs );

  WeaveStatisticType mainloop_stats = semi_res->statistics_params->mainloop_statistics;

  // Increment number of processed coherent results
  const size_t j = semi_res->ncoh_res;
  ++semi_res->ncoh_res;

  // Store per-segment coherent template index
  semi_res->coh_index[j] = coh_index;

  // Store per-segment coherent template parameters
  semi_res->coh_phys[j] = coh_res->coh_phys;
  semi_res->coh_phys[j].fkdot[0] += semi_res->dfreq * coh_offset;

  // Return now if simulating search
  if ( semi_res->simulation_level & WEAVE_SIMULATE ) {
    return XLAL_SUCCESS;
  }

  // Store per-segment F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_COH2F ) {
    semi_res->coh2F[j] = coh_res->coh2F->data + coh_offset;
  }

  // Store per-segment per-detector F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_COH2F_DET ) {
    for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
      semi_res->coh2F_det[i][j] = ( coh_res->coh2F_det[i] != NULL ) ? coh_res->coh2F_det[i]->data + coh_offset : NULL;
    }
  }

  // Start timing of semicoherent results
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_NONE, WEAVE_STATISTIC_MAX2F ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Add to max-over-segments multi-detector F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_MAX2F ) {
    XLAL_CHECK( semi_res_max_2F( &semi_res->nmax2F, semi_res->max2F->data, coh_res->coh2F->data + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Switch timed statistic
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_MAX2F, WEAVE_STATISTIC_MAX2F_DET ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Add to max-over-segments per-detector F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_MAX2F_DET ) {
    for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
      if ( coh_res->coh2F_det[i] != NULL ) {
        XLAL_CHECK( semi_res_max_2F( &semi_res->nmax2F_det[i], semi_res->max2F_det[i]->data, coh_res->coh2F_det[i]->data + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }
  }

  // Switch timed statistic
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_MAX2F_DET, WEAVE_STATISTIC_SUM2F ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Add to summed multi-detector F-statistics per frequency, and increment number of additions thus far
  if ( mainloop_stats & WEAVE_STATISTIC_SUM2F ) {
    XLAL_CHECK( semi_res_sum_2F( &semi_res->nsum2F, semi_res->sum2F->data, coh_res->coh2F->data + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
  } else {
    semi_res->nsum2F ++;             // Even if not summing here: count number of 2F summands for (potential) completion-loop usage
  }

  // Switch timed statistic
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_SUM2F, WEAVE_STATISTIC_SUM2F_DET ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Add to summed per-detector F-statistics per frequency, and increment number of additions thus far
  for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
    if ( coh_res->coh2F_det[i] != NULL ) {
      if ( mainloop_stats & WEAVE_STATISTIC_SUM2F_DET ) {
        XLAL_CHECK( semi_res_sum_2F( &semi_res->nsum2F_det[i], semi_res->sum2F_det[i]->data, coh_res->coh2F_det[i]->data + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
      } else {
        semi_res->nsum2F_det[i] ++;  // Even if not summing here: count number of per-detector 2F summands for (potential) completion-loop usage
      }
    }
  }

  // Stop timing of semicoherent results
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_SUM2F_DET, WEAVE_STATISTIC_NONE ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Compute all remaining *toplist-ranking* semicoherent statistics (ie 'mainloop-statistics').
/// For efficiency reasons any statistics not needed here will be computed later in the
/// "completion loop" on the final toplist.
///
int XLALWeaveSemiResultsComputeMain(
  WeaveSemiResults *semi_res,
  WeaveSearchTiming *tim
  )
{
  // Check input
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res->ncoh_res == semi_res->nsegments, XLAL_EINVAL );
  XLAL_CHECK( semi_res->statistics_params != NULL, XLAL_EFAULT );
  XLAL_CHECK( tim != NULL, XLAL_EFAULT );

  WeaveStatisticType mainloop_stats = semi_res->statistics_params->mainloop_statistics;

  // Return now if simulating search
  if ( semi_res->simulation_level & WEAVE_SIMULATE ) {
    return XLAL_SUCCESS;
  }

  // Store results from 'semi_res' in some convenience variables
  const REAL4 *sum2F = semi_res->sum2F != NULL ? semi_res->sum2F->data : NULL;
  const REAL4 *sum2F_det[PULSAR_MAX_DETECTORS];
  for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
    sum2F_det[i] = semi_res->sum2F_det[i] != NULL ? semi_res->sum2F_det[i]->data : NULL;
  }
  const REAL4 *max2F = semi_res->max2F != NULL ? semi_res->max2F->data : NULL;
  const REAL4 *max2F_det[PULSAR_MAX_DETECTORS];
  for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
    max2F_det[i] = semi_res->max2F_det[i] != NULL ? semi_res->max2F_det[i]->data : NULL;
  }

  //
  // Compute any remaining (ie that don't directly depend on coh_2F or coh2F_det) toplist ranking statistics
  //

  // Start timing of semicoherent results
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_NONE, WEAVE_STATISTIC_MEAN2F ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Compute mean multi-detector F-statistics per frequency:
  if ( mainloop_stats & WEAVE_STATISTIC_MEAN2F ) {
    XLAL_CHECK( XLALVectorScaleREAL4( semi_res->mean2F->data, 1.0 / semi_res->nsum2F, sum2F, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Switch timed statistic
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_MEAN2F, WEAVE_STATISTIC_BSGL ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Compute line-robust log10(B_S/GL) statistic per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_BSGL ) {
    XLAL_CHECK( XLALVectorComputeBSGL( semi_res->log10BSGL->data, sum2F, sum2F_det, semi_res->nfreqs, semi_res->statistics_params->BSGL_setup ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Switch timed statistic
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_BSGL, WEAVE_STATISTIC_BSGLtL ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Compute transient-line-robust log10(B_S/GL) statistic per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_BSGLtL ) {
    XLAL_CHECK( XLALVectorComputeBSGLtL( semi_res->log10BSGLtL->data, sum2F, sum2F_det, max2F_det, semi_res->nfreqs, semi_res->statistics_params->BSGL_setup ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Switch timed statistic
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_BSGLtL, WEAVE_STATISTIC_BtSGLtL ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Compute transient-signal line-robust log10(B_tS/GL) statistic per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_BtSGLtL ) {
    XLAL_CHECK( XLALVectorComputeBtSGLtL( semi_res->log10BtSGLtL->data, max2F, sum2F_det, max2F_det, semi_res->nfreqs, semi_res->statistics_params->BSGL_setup ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Stop timing of semicoherent results
  XLAL_CHECK( XLALWeaveSearchTimingStatistic( tim, WEAVE_STATISTIC_BtSGLtL, WEAVE_STATISTIC_NONE ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Destroy final semicoherent results
///
void XLALWeaveSemiResultsDestroy(
  WeaveSemiResults *semi_res
  )
{
  if ( semi_res == NULL ) {
    return;
  }

  XLALFree( semi_res->coh_index );
  XLALFree( semi_res->coh_phys );
  XLALFree( semi_res->coh2F );

  XLALDestroyREAL4VectorAligned( semi_res->max2F );
  XLALDestroyREAL4VectorAligned( semi_res->sum2F );
  XLALDestroyREAL4VectorAligned( semi_res->mean2F );

  for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
    XLALFree( semi_res->coh2F_det[i] );
    XLALDestroyREAL4VectorAligned( semi_res->max2F_det[i] );
    XLALDestroyREAL4VectorAligned( semi_res->sum2F_det[i] );
  }
  XLALDestroyREAL4VectorAligned( semi_res->log10BSGL );
  XLALDestroyREAL4VectorAligned( semi_res->log10BSGLtL );
  XLALDestroyREAL4VectorAligned( semi_res->log10BtSGLtL );

  XLALFree( semi_res );
  return;

}


/// Simple API function to extract pointers to 2F results from WeaveCohResults
int XLALWeaveCohResultsExtract(
  REAL4Vector **coh2F,
  REAL4Vector *coh2F_det[PULSAR_MAX_DETECTORS],
  BOOLEAN *have_coh2F_det,
  WeaveCohResults *coh_res,
  const WeaveCohInput *coh_input
  )
{
  XLAL_CHECK( coh_input != NULL, XLAL_EINVAL );
  XLAL_CHECK( coh_res != NULL, XLAL_EINVAL );
  XLAL_CHECK( coh_res->nfreqs >= 1, XLAL_EINVAL );
  XLAL_CHECK( coh2F != NULL && coh2F_det != NULL, XLAL_EINVAL );
  XLAL_CHECK( have_coh2F_det != NULL, XLAL_EINVAL );

  ( *coh2F ) = coh_res->coh2F;
  ( *have_coh2F_det ) = 0;
  if ( coh_input->Fstat_what_to_compute & FSTATQ_2F_PER_DET ) {
    ( *have_coh2F_det ) = 1;
    // Set all pointers to NULL first, the copy only the results from 'active' IFOs
    memset( coh2F_det, 0, PULSAR_MAX_DETECTORS * sizeof( coh2F_det[0] ) );
    for ( UINT4 i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
      const size_t idx = coh_input->Fstat_res_idx[i];
      coh2F_det[idx] = coh_res->coh2F_det[idx];
    }
  }

  return XLAL_SUCCESS;
}


// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
