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

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "ComputeResults.h"

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
  /// If computing per-detector quantities, list of detectors
  const LALStringVector *per_detectors;
  /// What F-statistic quantities to compute
  FstatQuantities what_to_compute;
  /// Number of detectors in F-statistic data
  UINT4 Fstat_ndetectors;
  /// Map detectors in F-statistic data to coherent results
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
/// Internal definition of partial results of a semicoherent computation in progress
///
struct tagWeaveSemiPartials {
  /// Bitflag representing search simulation level
  WeaveSimulationLevel simulation_level;
  /// Number of detectors
  UINT4 ndetectors;
  /// Number of segments
  UINT4 nsegments;
  /// Frequency spacing for semicoherent results
  double dfreq;
  /// Number of frequencies
  UINT4 nfreqs;
  /// Per-segment coherent template parameters of the first frequency bin (optional)
  PulsarDopplerParams *coh_phys;
  /// Per-segment multi-detector F-statistics per frequency (optional)
  const REAL4 **coh2F;
  /// Per-segment per-detector F-statistics per frequency (optional)
  const REAL4 **coh2F_det[PULSAR_MAX_DETECTORS];
  /// Number of coherent results processed thus far
  UINT4 ncoh_res;
  /// Semicoherent template parameters of the first frequency bin
  PulsarDopplerParams semi_phys;
  /// Summed multi-detector F-statistics per frequency
  REAL4VectorAligned *sum2F;
  /// Number of additions to multi-detector F-statistics thus far
  UINT4 nsum2F;
  /// Summed per-detector F-statistics per frequency
  REAL4VectorAligned *sum2F_det[PULSAR_MAX_DETECTORS];
  /// Number of additions to per-detector F-statistics thus far
  UINT4 nsum2F_det[PULSAR_MAX_DETECTORS];
};

///
/// \name Internal routines
///
/// @{

static int semi_parts_sum_2F( UINT4 *nsum, REAL4 *sum2F, const REAL4 *coh2F, const UINT4 nfreqs );

/// @}

///
/// Create coherent input data
///
WeaveCohInput *XLALWeaveCohInputCreate(
  const WeaveSimulationLevel simulation_level,
  FstatInput *Fstat_input,
  const LALStringVector *per_detectors
  )
{

  // Allocate memory
  WeaveCohInput *coh_input = XLALCalloc( 1, sizeof( *coh_input ) );
  XLAL_CHECK_NULL( coh_input != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL( ( simulation_level & WEAVE_SIMULATE_MIN_MEM ) || ( Fstat_input != NULL ), XLAL_EFAULT );

  // Set fields
  coh_input->simulation_level = simulation_level;
  coh_input->Fstat_input = Fstat_input;
  coh_input->per_detectors = per_detectors;

  // Decide what F-statistic quantities to compute
  coh_input->what_to_compute = FSTATQ_2F;
  if ( per_detectors != NULL ) {
    coh_input->what_to_compute |= FSTATQ_2F_PER_DET;
  }

  // If computing per-detector quantities, map detectors in F-statistic data to their index in the
  // coherent results. This is important when segments contain data from a subset of detectors.
  if ( per_detectors != NULL ) {

    // Get detectors in F-statistic data
    const MultiLALDetector *detector_info = XLALGetFstatInputDetectors( Fstat_input );
    coh_input->Fstat_ndetectors = detector_info->length;

    // Map entry 'i' in 'detector_info' (F-statistic data) to entry 'idx' in 'detectors' (coherent results)
    char *per_detectors_string = XLALConcatStringVector( per_detectors, "," );
    for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
      const char *prefix = detector_info->sites[i].frDetector.prefix;
      const int idx = XLALFindStringInVector( prefix, per_detectors );
      XLAL_CHECK_NULL( idx >= 0, XLAL_EFAILED, "Detector '%s' from F-statistic data not found in list of detectors '%s'", prefix, per_detectors_string );
      coh_input->Fstat_res_idx[i] = idx;
    }
    XLALFree( per_detectors_string );

  }

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

  // Double-check that F-statistic results detectors correctly match up with coherent results detectors
  for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
    const size_t idx = coh_input->Fstat_res_idx[i];
    const char *Fstat_ndetector = Fstat_res->detectorNames[i];
    const char *res_detector = coh_input->per_detectors->data[idx];
    XLAL_CHECK( strcmp( Fstat_ndetector, res_detector ) == 0, XLAL_EFAILED, "Detector #%zu in F-statistic results '%s' does not match detector #%zu in coherent results '%s'", i, Fstat_ndetector, idx, res_detector );
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
/// Create and initialise partial semicoherent results
///
int XLALWeaveSemiPartialsInit(
  WeaveSemiPartials **semi_parts,
  const WeaveSimulationLevel simulation_level,
  const UINT4 ndetectors,
  const UINT4 nsegments,
  const PulsarDopplerParams *semi_phys,
  const double dfreq,
  const UINT4 semi_nfreqs
  )
{

  // Check input
  XLAL_CHECK( semi_parts != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( dfreq >= 0, XLAL_EINVAL );
  XLAL_CHECK( semi_nfreqs > 0, XLAL_EINVAL );

  // Allocate results struct if required
  if ( *semi_parts == NULL ) {
    *semi_parts = XLALCalloc( 1, sizeof( **semi_parts ) );
    XLAL_CHECK( *semi_parts != NULL, XLAL_ENOMEM );
    ( *semi_parts )->coh_phys = XLALCalloc( nsegments, sizeof( *( *semi_parts )->coh_phys ) );
    XLAL_CHECK( ( *semi_parts )->coh_phys != NULL, XLAL_ENOMEM );
    ( *semi_parts )->coh2F = XLALCalloc( nsegments, sizeof( *( *semi_parts )->coh2F ) );
    XLAL_CHECK( ( *semi_parts )->coh2F != NULL, XLAL_ENOMEM );
    for ( size_t i = 0; i < ndetectors; ++i ) {
      ( *semi_parts )->coh2F_det[i] = XLALCalloc( nsegments, sizeof( *( *semi_parts )->coh2F_det[i] ) );
      XLAL_CHECK( ( *semi_parts )->coh2F_det[i] != NULL, XLAL_ENOMEM );
    }
  }

  // Set fields
  ( *semi_parts )->simulation_level = simulation_level;
  ( *semi_parts )->ndetectors = ndetectors;
  ( *semi_parts )->nsegments = nsegments;
  ( *semi_parts )->dfreq = dfreq;
  ( *semi_parts )->nfreqs = semi_nfreqs;
  ( *semi_parts )->semi_phys = *semi_phys;

  // Initialise number of coherent results
  ( *semi_parts )->ncoh_res = 0;

  // Initialise number of additions to multi- and per-detector F-statistics
  ( *semi_parts )->nsum2F = 0;
  XLAL_INIT_MEM( ( *semi_parts )->nsum2F_det );

  // Reallocate vectors of summed multi- and per-detector F-statistics per frequency
  if ( ( *semi_parts )->sum2F == NULL || ( *semi_parts )->sum2F->length < ( *semi_parts )->nfreqs ) {
    XLALDestroyREAL4VectorAligned( ( *semi_parts )->sum2F );
    ( *semi_parts )->sum2F = XLALCreateREAL4VectorAligned( ( *semi_parts )->nfreqs, alignment );
    XLAL_CHECK( ( *semi_parts )->sum2F != NULL, XLAL_ENOMEM );
  }
  for ( size_t i = 0; i < ( *semi_parts )->ndetectors; ++i ) {
    if ( ( *semi_parts )->sum2F_det[i] == NULL || ( *semi_parts )->sum2F_det[i]->length < ( *semi_parts )->nfreqs ) {
      XLALDestroyREAL4VectorAligned( ( *semi_parts )->sum2F_det[i] );
      ( *semi_parts )->sum2F_det[i] = XLALCreateREAL4VectorAligned( ( *semi_parts )->nfreqs, alignment );
      XLAL_CHECK( ( *semi_parts )->sum2F_det[i] != NULL, XLAL_ENOMEM );
    }
  }

  return XLAL_SUCCESS;

}

///
/// Destroy partial semicoherent results
///
void XLALWeaveSemiPartialsDestroy(
  WeaveSemiPartials *semi_parts
  )
{
  if ( semi_parts != NULL ) {
    XLALFree( semi_parts->coh_phys );
    XLALFree( semi_parts->coh2F );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALFree( semi_parts->coh2F_det[i] );
    }
    XLALDestroyREAL4VectorAligned( semi_parts->sum2F );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALDestroyREAL4VectorAligned( semi_parts->sum2F_det[i] );
    }
    XLALFree( semi_parts );
  }
}

///
/// Add F-statistic array 'coh2F' to summed array 'sum2F', and keep track of the number of summations 'nsum'
///
int semi_parts_sum_2F(
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
}

///
/// Add a new set of coherent results to the partial semicoherent results
///
int XLALWeaveSemiPartialsAdd(
  WeaveSemiPartials *semi_parts,
  const WeaveCohResults *coh_res,
  const UINT4 coh_offset
  )
{

  // Check input
  XLAL_CHECK( semi_parts != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_parts->ncoh_res < semi_parts->nsegments, XLAL_EINVAL );
  XLAL_CHECK( coh_res != NULL, XLAL_EFAULT );

  // Check that offset does not overrun coherent results arrays
  XLAL_CHECK( coh_offset + semi_parts->nfreqs <= coh_res->nfreqs, XLAL_EFAILED, "Coherent offset (%u) + number of semicoherent frequency bins (%u) > number of coherent frequency bins (%u)", coh_offset, semi_parts->nfreqs, coh_res->nfreqs );

  // Increment number of processed coherent results
  const size_t j = semi_parts->ncoh_res;
  ++semi_parts->ncoh_res;

  // Store per-segment coherent template parameters
  semi_parts->coh_phys[j] = coh_res->coh_phys;
  semi_parts->coh_phys[j].fkdot[0] += semi_parts->dfreq * coh_offset;

  // Return now if simulating search
  if ( semi_parts->simulation_level & WEAVE_SIMULATE ) {
    return XLAL_SUCCESS;
  }

  // Store per-segment F-statistics per frequency
  semi_parts->coh2F[j] = coh_res->coh2F->data + coh_offset;
  for ( size_t i = 0; i < semi_parts->ndetectors; ++i ) {
    if ( coh_res->coh2F_det[i] != NULL ) {
      semi_parts->coh2F_det[i][j] = coh_res->coh2F_det[i]->data + coh_offset;
    } else {
      semi_parts->coh2F_det[i][j] = NULL;
    }
  }

  // Add to summed multi-detector F-statistics per frequency, and increment number of additions thus far
  XLAL_CHECK( semi_parts_sum_2F( &semi_parts->nsum2F, semi_parts->sum2F->data, semi_parts->coh2F[j], semi_parts->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Add to summed per-detector F-statistics per frequency, and increment number of additions thus far
  for ( size_t i = 0; i < semi_parts->ndetectors; ++i ) {
    if ( semi_parts->coh2F_det[i][j] != NULL ) {
      XLAL_CHECK( semi_parts_sum_2F( &semi_parts->nsum2F_det[i], semi_parts->sum2F_det[i]->data, semi_parts->coh2F_det[i][j], semi_parts->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  return XLAL_SUCCESS;

}

///
/// Create and compute final semicoherent results
///
int XLALWeaveSemiResultsCompute(
  WeaveSemiResults **semi_res,
  const WeaveSemiPartials *semi_parts
  )
{

  // Check input
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_parts != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_parts->ncoh_res == semi_parts->nsegments, XLAL_EINVAL );

  // Allocate results struct if required
  if ( *semi_res == NULL ) {
    *semi_res = XLALCalloc( 1, sizeof( **semi_res ) );
    XLAL_CHECK( *semi_res != NULL, XLAL_ENOMEM );
  }

  // Set fields
  ( *semi_res )->simulation_level = semi_parts->simulation_level;
  ( *semi_res )->ndetectors = semi_parts->ndetectors;
  ( *semi_res )->nsegments = semi_parts->nsegments;
  ( *semi_res )->dfreq = semi_parts->dfreq;
  ( *semi_res )->nfreqs = semi_parts->nfreqs;
  ( *semi_res )->coh_phys = semi_parts->coh_phys;
  ( *semi_res )->coh2F = semi_parts->coh2F;
  memcpy( ( *semi_res )->coh2F_det, semi_parts->coh2F_det, sizeof( ( *semi_res )->coh2F_det ) );
  ( *semi_res )->semi_phys = semi_parts->semi_phys;

  // Reallocate vectors of mean multi- and per-detector F-statistics per frequency
  if ( ( *semi_res )->mean2F == NULL || ( *semi_res )->mean2F->length < semi_parts->nfreqs ) {
    XLALDestroyREAL4VectorAligned( ( *semi_res )->mean2F );
    ( *semi_res )->mean2F = XLALCreateREAL4VectorAligned( semi_parts->nfreqs, alignment );
    XLAL_CHECK( ( *semi_res )->mean2F != NULL, XLAL_ENOMEM );
  }
  for ( size_t i = 0; i < semi_parts->ndetectors; ++i ) {
    if ( ( *semi_res )->mean2F_det[i] == NULL || ( *semi_res )->mean2F_det[i]->length < semi_parts->nfreqs ) {
      XLALDestroyREAL4VectorAligned( ( *semi_res )->mean2F_det[i] );
      ( *semi_res )->mean2F_det[i] = XLALCreateREAL4VectorAligned( semi_parts->nfreqs, alignment );
      XLAL_CHECK( ( *semi_res )->mean2F_det[i] != NULL, XLAL_ENOMEM );
    }
  }

  // Return now if simulating search
  if ( ( *semi_res )->simulation_level & WEAVE_SIMULATE ) {
    return XLAL_SUCCESS;
  }

  // Calculate mean multi-detector F-statistics per frequency
  XLAL_CHECK( XLALVectorScaleREAL4( ( *semi_res )->mean2F->data, 1.0 / semi_parts->nsum2F, semi_parts->sum2F->data, ( *semi_res )->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Calculate mean per-detector F-statistics per frequency
  for ( size_t i = 0; i < ( *semi_res )->ndetectors; ++i ) {
    XLAL_CHECK( XLALVectorScaleREAL4( ( *semi_res )->mean2F_det[i]->data, 1.0 / semi_parts->nsum2F_det[i], semi_parts->sum2F_det[i]->data, ( *semi_res )->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}

///
/// Destroy final semicoherent results
///
void XLALWeaveSemiResultsDestroy(
  WeaveSemiResults *semi_res
  )
{
  if ( semi_res != NULL ) {
    XLALDestroyREAL4VectorAligned( semi_res->mean2F );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALDestroyREAL4VectorAligned( semi_res->mean2F_det[i] );
    }
    XLALFree( semi_res );
  }
}

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
