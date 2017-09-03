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
// Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307 USA
//

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "ComputeResults.h"

#include <lal/UserInputPrint.h>

///
/// Aligned arrays use maximum required alignment, i.e.\ 32 bytes for AVX
///
const UINT4 alignment = 32;

///
/// Internal definition of input data required for computing coherent results
///
struct tagWeaveCohInput {
  /// Bitflag representing search simulation level
  WeaveSimulationLevel simulation_level;
  /// F-statistic input data
  FstatInput *Fstat_input;
  /// What F-statistic quantities to compute
  FstatQuantities what_to_compute;
  /// Number of detectors in F-statistic data
  UINT4 Fstat_ndetectors;
  /// Map detectors in F-statistic data in this segment to 'global' detector list across segments in coherent results
  size_t Fstat_res_idx[PULSAR_MAX_DETECTORS];
};

///
/// Internal definition of results of a coherent computation on a single segment
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
  const WeaveSimulationLevel simulation_level,
  FstatInput *Fstat_input,
  const WeaveStatisticsParams *statistics_params
  )
{

  // Allocate memory
  WeaveCohInput *coh_input = XLALCalloc( 1, sizeof( *coh_input ) );
  XLAL_CHECK_NULL( coh_input != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL( ( simulation_level & WEAVE_SIMULATE_MIN_MEM ) || ( Fstat_input != NULL ), XLAL_EFAULT );
  XLAL_CHECK_NULL( statistics_params != NULL, XLAL_EFAULT );

  // Set fields
  coh_input->simulation_level = simulation_level;
  coh_input->Fstat_input = Fstat_input;

  WeaveStatisticType mainloop_stats = statistics_params -> mainloop_statistics;

  // Decide what F-statistic quantities to compute
  if ( mainloop_stats & WEAVE_STATISTIC_COH2F ) {
    coh_input->what_to_compute |= FSTATQ_2F;
  }
  if ( mainloop_stats & WEAVE_STATISTIC_COH2F_DET ) {
    coh_input->what_to_compute |= FSTATQ_2F_PER_DET;
  }

  // Map detectors in F-statistic data in this segment to their index in the coherent results.
  // This is important when segments contain data from a subset of detectors.
  if ( !(simulation_level & WEAVE_SIMULATE_MIN_MEM) ) {

    // Get detectors in F-statistic data
    const MultiLALDetector *detector_info = XLALGetFstatInputDetectors( Fstat_input );
    coh_input->Fstat_ndetectors = detector_info->length;

    // Map entry 'i' in 'detector_info' (F-statistic data) to entry 'idx' in 'detectors' (coherent results)
    char *detectors_string = XLALConcatStringVector( statistics_params->detectors, "," );
    for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
      const char *prefix = detector_info->sites[i].frDetector.prefix;
      const int idx = XLALFindStringInVector( prefix, statistics_params->detectors );
      XLAL_CHECK_NULL( idx >= 0, XLAL_EFAILED, "Detector '%s' from F-statistic data not found in list of detectors '%s'", prefix, detectors_string );
      coh_input->Fstat_res_idx[i] = idx;
    }
    XLALFree( detectors_string );

  }

  return coh_input;

} // XLALWeaveCohInputCreate()

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
} // XLALWeaveCohInputDestroy()

///
/// Write F-statistic method name from coherent input data to a FITS file
///
int XLALWeaveCohInputWriteFstatMethod(
  FITSFile *file,
  const WeaveCohInput *coh_input
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( coh_input != NULL, XLAL_EFAULT );

  // Get method name
  const char *method_name = XLALGetFstatInputMethodName( coh_input->Fstat_input );
  XLAL_CHECK( method_name != NULL, XLAL_EFUNC );

  // Write method name
  XLAL_CHECK( XLALFITSHeaderWriteString( file, "fmethod", method_name, "name of F-statistic method" ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Write F-statistic timing information from coherent input data to a FITS file
///
int XLALWeaveCohInputWriteFstatTiming(
  FITSFile *file,
  const size_t ncoh_input,
  WeaveCohInput *const *coh_input
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( ncoh_input > 0, XLAL_ESIZE );
  XLAL_CHECK( coh_input != NULL, XLAL_EFAULT );

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
  XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, "ftiming tauF_eff", tauF_eff / ncoh_input, "F-statistic generic timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, "ftiming tauF_core", tauF_core / ncoh_input, "F-statistic generic timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, "ftiming tauF_buffer", tauF_buffer / ncoh_input, "F-statistic generic timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, "ftiming b", NBufferMisses / NCalls, "F-statistic generic timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write method-specific timing constants
  for ( size_t j = 0; j < nmodel; ++j ) {
    char keyword[64];
    snprintf( keyword, sizeof( keyword ), "ftiming %s", model_names[j] );
    XLAL_CHECK( XLALFITSHeaderWriteREAL4( file, keyword, model_values[j] / ncoh_input, "F-statistic method-specific timing constant" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} // XLALWeaveCohInputWriteTiming()

///
/// Create and compute coherent results
///
int XLALWeaveCohResultsCompute(
  WeaveCohResults **coh_res,
  WeaveCohInput *coh_input,
  const PulsarDopplerParams *coh_phys,
  const UINT4 coh_nfreqs
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
  for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
    const size_t idx = coh_input->Fstat_res_idx[i];
    if ( ( *coh_res )->coh2F_det[idx] == NULL || ( *coh_res )->coh2F_det[idx]->length < ( *coh_res )->nfreqs ) {
      ( *coh_res )->coh2F_det[idx] = XLALResizeREAL4Vector( ( *coh_res )->coh2F_det[idx], ( *coh_res )->nfreqs );
      XLAL_CHECK( ( *coh_res )->coh2F_det[idx] != NULL, XLAL_ENOMEM );
    }
  }

  // Return now if simulating search
  if ( coh_input->simulation_level & WEAVE_SIMULATE ) {
    return XLAL_SUCCESS;
  }

  // Use a local F-statistic results structure, since we supply our own memory
  // - The 'internalalloclen' field stores the memory size in elements of the results arrays (e.g. 'coh2F'),
  //   as opposed to 'numFreqBins' which stores how many elements of the result arrays are in use.
  FstatResults XLAL_INIT_DECL( Fstat_res_struct );
  FstatResults *Fstat_res = &Fstat_res_struct;
  Fstat_res->internalalloclen = ( *coh_res )->nfreqs;
  if ( coh_input->what_to_compute & FSTATQ_2F_PER_DET ) {
    Fstat_res->numDetectors = coh_input->Fstat_ndetectors;
  }
  Fstat_res->twoF = ( *coh_res )->coh2F->data;
  for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
    const size_t idx = coh_input->Fstat_res_idx[i];
    Fstat_res->twoFPerDet[i] = ( *coh_res )->coh2F_det[idx]->data;
  }

  // Compute the F-statistic starting at the point 'coh_phys', with 'nfreqs' frequency bins
  XLAL_CHECK( XLALComputeFstat( &Fstat_res, coh_input->Fstat_input, coh_phys, ( *coh_res )->nfreqs, coh_input->what_to_compute ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Sanity check the F-statistic results structure
  XLAL_CHECK( Fstat_res->internalalloclen == ( *coh_res )->nfreqs, XLAL_EFAILED );
  XLAL_CHECK( Fstat_res->numFreqBins == ( *coh_res )->nfreqs, XLAL_EFAILED );
  XLAL_CHECK( !( coh_input->what_to_compute & FSTATQ_2F_PER_DET ) || ( Fstat_res->numDetectors == coh_input->Fstat_ndetectors ), XLAL_EFAILED );
  XLAL_CHECK( Fstat_res->twoF == ( *coh_res )->coh2F->data, XLAL_EFAILED );
  for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
    const size_t idx = coh_input->Fstat_res_idx[i];
    XLAL_CHECK( Fstat_res->twoFPerDet[i] == ( *coh_res )->coh2F_det[idx]->data, XLAL_EFAILED, "%p vs %p", Fstat_res->twoFPerDet[i], ( *coh_res )->coh2F_det[idx]->data );
  }

  return XLAL_SUCCESS;

} // XLALWeaveCohResultsCompute()

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
} // XLALWeaveCohResultsDestroy()

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

  WeaveStatisticType unsupported = (mainloop_stats & ~supported_mainloop);
  if ( unsupported != 0 ) {
    char *unsupported_names = XLALPrintStringValueOfUserFlag( (const int*)&unsupported, &WeaveStatisticChoices );
    XLALPrintError ( "BUG: unsupported main-loop statistics requested: %s\n", unsupported_names );
    XLALFree ( unsupported_names );
    XLAL_ERROR ( XLAL_EERR );
  }

  // Allocate results struct if required
  if ( *semi_res == NULL ) {
    *semi_res = XLALCalloc( 1, sizeof( **semi_res ) );
    XLAL_CHECK( *semi_res != NULL, XLAL_ENOMEM );

    // Allocate array for per-segment index
    ( *semi_res )->coh_index = XLALCalloc( nsegments, sizeof( *( *semi_res )->coh_index ) );
    XLAL_CHECK( ( *semi_res )->coh_index != NULL, XLAL_ENOMEM );

    // allocate array for per-segment coordinates
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

  } // if *semi_res == NULL

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

} // XLALWeaveSemiResultsInit()

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
  if ( (*nsum)++ == 0 ) {
    // If this is the first summation, just use memcpy()
    memcpy( sum2F, coh2F, sizeof( *sum2F ) * nfreqs );
    return XLAL_SUCCESS;
  }
  return XLALVectorAddREAL4( sum2F, sum2F, coh2F, nfreqs );
} // semi_res_sum_2F()

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
  if ( (*nmax)++ == 0 ) {
    // If this is the first max-comparison, just use memcpy()
    memcpy( max2F, coh2F, sizeof( *max2F ) * nfreqs );
    return XLAL_SUCCESS;
  }
  return XLALVectorMaxREAL4( max2F, max2F, coh2F, nfreqs );
} // semi_res_max_2F()

///
/// Add a new set of coherent results to the semicoherent results
///
int XLALWeaveSemiResultsAdd(
  WeaveSemiResults *semi_res,
  const WeaveCohResults *coh_res,
  const UINT8 coh_index,
  const UINT4 coh_offset
  )
{

  // Check input
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res->ncoh_res < semi_res->nsegments, XLAL_EINVAL );
  XLAL_CHECK( coh_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res->statistics_params != NULL, XLAL_EFAULT );

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
    } // for i < ndetectors
  }

  // Add to max-over-segments multi-detector F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_MAX2F ) {
    XLAL_CHECK( semi_res_max_2F( &semi_res->nmax2F, semi_res->max2F->data, coh_res->coh2F->data + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  // Add to max-over-segments per-detector F-statistics per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_MAX2F_DET ) {
    for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
      if ( coh_res->coh2F_det[i] != NULL ) {
        XLAL_CHECK( semi_res_max_2F( &semi_res->nmax2F_det[i], semi_res->max2F_det[i]->data, coh_res->coh2F_det[i]->data + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }
  }

  // Add to summed multi-detector F-statistics per frequency, and increment number of additions thus far
  if ( mainloop_stats & WEAVE_STATISTIC_SUM2F ) {
    XLAL_CHECK( semi_res_sum_2F( &semi_res->nsum2F, semi_res->sum2F->data, coh_res->coh2F->data + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
  } else {
    semi_res->nsum2F ++;             // even if not summing here: count number of 2F summands for (potential) completion-loop usage
  }
  // Add to summed per-detector F-statistics per frequency, and increment number of additions thus far
  for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
    if ( coh_res->coh2F_det[i] != NULL ) {
      if ( mainloop_stats & WEAVE_STATISTIC_SUM2F_DET ) {
        XLAL_CHECK( semi_res_sum_2F( &semi_res->nsum2F_det[i], semi_res->sum2F_det[i]->data, coh_res->coh2F_det[i]->data + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
      } else {
        semi_res->nsum2F_det[i] ++;  // even if not summing here: count number of per-detector 2F summands for (potential) completion-loop usage
      }
    }
  }

  return XLAL_SUCCESS;

} // XLALWeaveSemiResultsAdd()

///
/// Compute all remaining *toplist-ranking* semicoherent statistics (ie 'mainloop-statistics').
/// For efficiency reasons any statistics not needed here will be computed later in the
/// "completion loop" on the final toplist.
///
int XLALWeaveSemiResultsComputeMain(
  WeaveSemiResults *semi_res
  )
{
  // Check input
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res->ncoh_res == semi_res->nsegments, XLAL_EINVAL );
  XLAL_CHECK( semi_res->statistics_params != NULL, XLAL_EFAULT );

  // Return now if simulating search
  if ( semi_res->simulation_level & WEAVE_SIMULATE ) {
    return XLAL_SUCCESS;
  }

  //
  // Compute any remaining (ie that don't directly depend on coh_2F or coh2F_det) toplist ranking statistics
  //
  WeaveStatisticType mainloop_stats = semi_res -> statistics_params -> mainloop_statistics;

  // mean multi-detector F-statistics per frequency:
  if ( mainloop_stats & WEAVE_STATISTIC_MEAN2F ) {
    XLAL_CHECK( XLALVectorScaleREAL4( semi_res->mean2F->data, 1.0 / semi_res->nsum2F, semi_res->sum2F->data, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // various line-robust statistics
  const REAL4 *sum2F_det[PULSAR_MAX_DETECTORS];
  const REAL4 *max2F_det[PULSAR_MAX_DETECTORS];
  if ( mainloop_stats & ( WEAVE_STATISTIC_BSGL | WEAVE_STATISTIC_BSGLtL |  WEAVE_STATISTIC_BtSGLtL ) ) {
    for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
      sum2F_det[i] = semi_res->sum2F_det[i]->data;
    }
  }
  if ( mainloop_stats & (WEAVE_STATISTIC_BSGLtL | WEAVE_STATISTIC_BtSGLtL) ) {
    for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
      max2F_det[i] = semi_res->max2F_det[i]->data;
    }
  }

  // line-robust log10(B_S/GL) statistic per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_BSGL ) {
    XLAL_CHECK ( XLALVectorComputeBSGL ( semi_res->log10BSGL->data, semi_res->sum2F->data, sum2F_det, semi_res->nfreqs, semi_res->statistics_params->BSGL_setup ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // transient-line-robust log10(B_S/GL) statistic per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_BSGLtL ) {
    XLAL_CHECK ( XLALVectorComputeBSGLtL ( semi_res->log10BSGLtL->data, semi_res->sum2F->data, sum2F_det, max2F_det, semi_res->nfreqs, semi_res->statistics_params->BSGL_setup ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // transient-signal line-robust log10(B_tS/GL) statistic per frequency
  if ( mainloop_stats & WEAVE_STATISTIC_BtSGLtL ) {
    XLAL_CHECK ( XLALVectorComputeBtSGLtL ( semi_res->log10BtSGLtL->data, semi_res->max2F->data, sum2F_det, max2F_det, semi_res->nfreqs, semi_res->statistics_params->BSGL_setup ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} // XLALWeaveSemiResultsComputeMain()

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

} // XLALWeaveSemiResultsDestroy()

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
