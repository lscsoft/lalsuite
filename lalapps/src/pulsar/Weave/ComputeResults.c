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

#include <lal/VectorMath.h>

// Aligned arrays use maximum required alignment, i.e. 32 bytes for AVX
const UINT4 alignment = 32;

///
/// Internal definition of input data required for computing coherent results
///
struct tagWeaveCohInput {
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
  REAL4 *twoF;
  /// Per-detector F-statistics per frequency
  REAL4 *twoF_per_det[PULSAR_MAX_DETECTORS];
};

///
/// Internal definition of results of a semicoherent computation over many segments
///
struct tagWeaveSemiResults {
  /// Per-segment coherent results (optional)
  WeaveCohResults *coh_res;
  /// Semicoherent template parameters of the first frequency bin
  PulsarDopplerParams semi_phys;
  /// Frequency spacing for semicoherent results
  double dfreq;
  /// Number of frequencies
  UINT4 nfreqs;
  /// Number of segments added to multi-detector results thus far
  UINT4 nsegments;
  /// Summed multi-detector F-statistics per frequency
  REAL4VectorAligned *sum_twoF;
  /// Number of detectors for per-detector results
  UINT4 ndetectors;
  /// Number of segments added to per-detector results thus far
  UINT4 nsegments_per_det[PULSAR_MAX_DETECTORS];
  /// Summed per-detector F-statistics per frequency
  REAL4VectorAligned *sum_twoF_per_det[PULSAR_MAX_DETECTORS];
};

///
/// \name Internal routines
///
/// @{

static int semi_results_add_REAL4( const size_t nsum, REAL4 *sum, const REAL4 *x, const size_t nbins );

/// @}

///
/// Create coherent input data
///
WeaveCohInput *XLALWeaveCohInputCreate(
  FstatInput *Fstat_input,
  const LALStringVector *per_detectors
  )
{

  // Check input
  XLAL_CHECK_NULL( Fstat_input != NULL, XLAL_EFAULT );

  // Allocate memory
  WeaveCohInput *coh_input = XLALCalloc( 1, sizeof( *coh_input ) );
  XLAL_CHECK_NULL( coh_input != NULL, XLAL_ENOMEM );

  // Set fields
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

  // Store coherent template parameters
  ( *coh_res )->coh_phys = *coh_phys;

  // Store number of coherent frequency bins
  const UINT4 old_coh_res_nfreqs = ( *coh_res )->nfreqs;
  ( *coh_res )->nfreqs = coh_nfreqs;

  // Reallocate arrays if required
  if ( old_coh_res_nfreqs < coh_nfreqs ) {

    // Reallocate arrays of multi- and per-detector F-statistics per frequency
    ( *coh_res )->twoF = XLALRealloc( ( *coh_res )->twoF, coh_nfreqs * sizeof( ( *coh_res )->twoF[0] ) );
    XLAL_CHECK( ( *coh_res )->twoF != NULL, XLAL_ENOMEM );
    for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
      const size_t idx = coh_input->Fstat_res_idx[i];
      ( *coh_res )->twoF_per_det[idx] = XLALRealloc( ( *coh_res )->twoF_per_det[idx], coh_nfreqs * sizeof( ( *coh_res )->twoF_per_det[idx][0] ) );
      XLAL_CHECK( ( *coh_res )->twoF_per_det[idx] != NULL, XLAL_ENOMEM );
    }

  }

  // Use a local F-statistic results structure, since we supply our own memory
  // - The 'internalalloclen' field stores the memory size in elements of the results arrays (e.g. 'twoF'),
  //   as opposed to 'numFreqBins' which stores how many elements of the result arrays are in use.
  FstatResults XLAL_INIT_DECL( Fstat_res_struct );
  FstatResults *Fstat_res = &Fstat_res_struct;
  Fstat_res->internalalloclen = ( *coh_res )->nfreqs;
  Fstat_res->twoF = ( *coh_res )->twoF;
  for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
    const size_t idx = coh_input->Fstat_res_idx[i];
    Fstat_res->twoFPerDet[i] = ( *coh_res )->twoF_per_det[idx];
  }

  // Compute the F-statistic starting at the point 'coh_phys', with 'nfreqs' frequency bins
  XLAL_CHECK( XLALComputeFstat( &Fstat_res, coh_input->Fstat_input, coh_phys, coh_nfreqs, coh_input->what_to_compute ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Sanity check the F-statistic results structure
  XLAL_CHECK( Fstat_res->internalalloclen == ( *coh_res )->nfreqs, XLAL_EFAILED );
  XLAL_CHECK( Fstat_res->numFreqBins == ( *coh_res )->nfreqs, XLAL_EFAILED );
  XLAL_CHECK( Fstat_res->twoF == ( *coh_res )->twoF, XLAL_EFAILED );
  for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
    const size_t idx = coh_input->Fstat_res_idx[i];
    XLAL_CHECK( Fstat_res->twoFPerDet[i] == ( *coh_res )->twoF_per_det[idx], XLAL_EFAILED );
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
    XLALFree( coh_res->twoF );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALFree( coh_res->twoF_per_det[i] );
    }
    XLALFree( coh_res );
  }
}

///
/// Create storage for semicoherent results
///
WeaveSemiResults *XLALWeaveSemiResultsCreate(
  const LALStringVector *per_detectors,
  const UINT4 per_nsegments,
  const double dfreq
  )
{

  // Check input
  XLAL_CHECK_NULL( dfreq >= 0, XLAL_EINVAL );

  // Allocate memory
  WeaveSemiResults *semi_res = XLALCalloc( 1, sizeof( *semi_res ) );
  XLAL_CHECK_NULL( semi_res != NULL, XLAL_ENOMEM );
  if ( per_nsegments > 0 ) {
    semi_res->coh_res = XLALCalloc( per_nsegments, sizeof( *semi_res->coh_res ) );
    XLAL_CHECK_NULL( semi_res->coh_res != NULL, XLAL_ENOMEM );
  }

  // Set fields
  semi_res->dfreq = dfreq;

  // Initialise number of detectors for per-detector results
  semi_res->ndetectors = ( per_detectors != NULL ) ? per_detectors->length : 0;

  return semi_res;


}

///
/// Destroy storage for semicoherent results
///
void XLALWeaveSemiResultsDestroy(
  WeaveSemiResults *semi_res
  )
{
  if ( semi_res != NULL ) {
    XLALFree( semi_res->coh_res );
    XLALDestroyREAL4VectorAligned( semi_res->sum_twoF );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALDestroyREAL4VectorAligned( semi_res->sum_twoF_per_det[i] );
    }
    XLALFree( semi_res );
  }
}

///
/// Create and initialise semicoherent results
///
int XLALWeaveSemiResultsInit(
  WeaveSemiResults *semi_res,
  const PulsarDopplerParams *semi_phys,
  const UINT4 semi_nfreqs
  )
{

  // Check input
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_nfreqs > 0, XLAL_EINVAL );

  // Initialise number of segments for which results have been added
  semi_res->nsegments = 0;
  XLAL_INIT_MEM( semi_res->nsegments_per_det );

  // Copy start point of the current semicoherent frequency block and number of frequency points
  semi_res->semi_phys = *semi_phys;
  semi_res->nfreqs = semi_nfreqs;

  // Reallocate arrays of summed multi- and per-detector F-statistics per frequency
  if ( semi_res->sum_twoF == NULL || semi_res->sum_twoF->length < semi_res->nfreqs ) {
    XLALDestroyREAL4VectorAligned( semi_res->sum_twoF );
    semi_res->sum_twoF = XLALCreateREAL4VectorAligned( semi_res->nfreqs, alignment );
    XLAL_CHECK( semi_res->sum_twoF != NULL, XLAL_ENOMEM );
  }
  for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
    if ( semi_res->sum_twoF_per_det[i] == NULL || semi_res->sum_twoF_per_det[i]->length < semi_res->nfreqs ) {
      XLALDestroyREAL4VectorAligned( semi_res->sum_twoF_per_det[i] );
      semi_res->sum_twoF_per_det[i] = XLALCreateREAL4VectorAligned( semi_res->nfreqs, alignment );
      XLAL_CHECK( semi_res->sum_twoF_per_det[i] != NULL, XLAL_ENOMEM );
    }
  }

  return XLAL_SUCCESS;

}

///
/// Add 'nbins' of REAL4 array 'x' to 'sum', and keep track of the number of summations 'nsum'
///
int semi_results_add_REAL4( const size_t nsum, REAL4 *sum, const REAL4 *x, const size_t nbins )
{
  if ( nsum == 0 ) { // If this is the first summation, just use memcpy()
    memcpy( sum, x, sizeof( *sum ) * nbins );
    return XLAL_SUCCESS;
  }
  return XLALVectorAddREAL4( sum, sum, x, nbins );
}

///
/// Add a new set of coherent results to the semicoherent results
///
int XLALWeaveSemiResultsAdd(
  WeaveSemiResults *semi_res,
  const WeaveCohResults *coh_res,
  const UINT4 coh_offset
  )
{

  // Check input
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( coh_res != NULL, XLAL_EFAULT );

  // Check that offset does not overrun coherent results arrays
  XLAL_CHECK( coh_offset + semi_res->nfreqs <= coh_res->nfreqs, XLAL_EFAILED, "Coherent offset (%u) + number of semicoherent frequency bins (%u) > number of coherent frequency bins (%u)", coh_offset, semi_res->nfreqs, coh_res->nfreqs );

  // Store per-segment coherent template parameters and F-statistics per frequency
  if ( semi_res->coh_res != NULL ) {
    WeaveCohResults *semi_res_coh_res = &semi_res->coh_res[semi_res->nsegments];
    *semi_res_coh_res = *coh_res;

    // Offset coherent template frequency
    semi_res_coh_res->coh_phys.fkdot[0] += semi_res->dfreq * coh_offset;

    // Offset arrays of coherent F-statistics per frequency
    semi_res_coh_res->twoF += coh_offset;
    for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
      if ( semi_res_coh_res->twoF_per_det[i] != NULL ) {
        semi_res_coh_res->twoF_per_det[i] += coh_offset;
      }
    }

  }

  // Add to summed multi-detector F-statistics per frequency
  XLAL_CHECK( semi_results_add_REAL4( semi_res->nsegments, semi_res->sum_twoF->data, coh_res->twoF + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Increment number of segments added to multi-detector results thus far
  ++semi_res->nsegments;

  for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
    if ( coh_res->twoF_per_det[i] != NULL ) {

      // Add to summed per-detector F-statistics per frequency
      XLAL_CHECK( semi_results_add_REAL4( semi_res->nsegments_per_det[i], semi_res->sum_twoF_per_det[i]->data, coh_res->twoF_per_det[i] + coh_offset, semi_res->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

      // Increment number of segments added to this detector results thus far
      ++semi_res->nsegments_per_det[i];

    }
  }

  return XLAL_SUCCESS;

}

///
/// Fill a output result item, creating a new one if needed
///
int XLALWeaveFillOutputResultItem(
  WeaveOutputResultItem **item,
  BOOLEAN *full_init,
  const WeaveSemiResults *semi_res,
  const size_t freq_idx
  )
{

  // Check input
  XLAL_CHECK( item != NULL, XLAL_EFAULT );
  XLAL_CHECK( full_init != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( freq_idx < semi_res->nfreqs, XLAL_EINVAL );

  const UINT4 per_nsegments = ( semi_res->coh_res != NULL ) ? semi_res->nsegments : 0;
  const UINT4 per_ndetectors = semi_res->ndetectors;

  // Fully initialise all output result item field, if requested
  // - Otherwise only output result item fields that change with 'freq_idx' are updated
  if ( *full_init ) {

    // Set all semicoherent and coherent template parameters
    ( *item )->semi_phys = semi_res->semi_phys;
    for ( size_t s = 0; s < per_nsegments; ++s ) {
      ( *item )->per_seg[s].coh_phys = semi_res->coh_res[s].coh_phys;
    }

    // Next time, only output result item fields that change with 'freq_idx' should need updating
    *full_init = 0;

  }

  // Update semicoherent and coherent template frequency
  const double foffset = freq_idx * semi_res->dfreq;
  ( *item )->semi_phys.fkdot[0] = semi_res->semi_phys.fkdot[0] + foffset;
  for ( size_t s = 0; s < per_nsegments; ++s ) {
    ( *item )->per_seg[s].coh_phys.fkdot[0] = semi_res->coh_res[s].coh_phys.fkdot[0] + foffset;
  }

  // Update multi-detector F-statistics
  ( *item )->mean_twoF = semi_res->sum_twoF->data[freq_idx] / semi_res->nsegments;
  for ( size_t s = 0; s < per_nsegments; ++s ) {
    ( *item )->per_seg[s].twoF = semi_res->coh_res[s].twoF[freq_idx];
  }

  // Update per-detector F-statistics
  for ( size_t i = 0; i < per_ndetectors; ++i ) {
    ( *item )->mean_twoF_per_det[i] = semi_res->sum_twoF_per_det[i]->data[freq_idx] / semi_res->nsegments_per_det[i];
    for ( size_t s = 0; s < per_nsegments; ++s ) {
      if ( semi_res->coh_res[s].twoF_per_det[i] != NULL ) {
        ( *item )->per_seg[s].twoF_per_det[i] = semi_res->coh_res[s].twoF_per_det[i][freq_idx];
      } else {
        // There is not per-detector F-statistic for this segment, usually because this segment contains
        // no data from this detector. In this case we output a clearly invalid F-statistic value.
        ( *item )->per_seg[s].twoF_per_det[i] = NAN;
      }
    }
  }

  return XLAL_SUCCESS;

}
