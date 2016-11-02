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

///
/// Aligned arrays use maximum required alignment, i.e.\ 32 bytes for AVX
///
const UINT4 alignment = 32;

///
/// Internal definition of input data required for computing coherent results
///
struct tagWeaveCohInput {
  /// Whether to shortcut computing coherent results
  BOOLEAN shortcut_compute;
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
/// Internal definition of partial results of a semicoherent computation in progress
///
struct tagWeaveSemiPartials {
  /// Whether to shortcut computing semicoherent results
  BOOLEAN shortcut_compute;
  /// Per-segment coherent results (optional)
  WeaveCohResults *coh_res;
  /// Number of per-segment coherent results added thus far
  UINT4 ncoh_res;
  /// Semicoherent template parameters of the first frequency bin
  PulsarDopplerParams semi_phys;
  /// Frequency spacing for semicoherent results
  double dfreq;
  /// Number of frequencies
  UINT4 nfreqs;
  /// Number of detectors for per-detector results
  UINT4 ndetectors;
  /// Summed multi-detector F-statistics per frequency
  REAL4VectorAligned *sum_twoF;
  /// Number of additions to multi-detector F-statistics thus far
  UINT4 nsum_twoF;
  /// Summed per-detector F-statistics per frequency
  REAL4VectorAligned *sum_twoF_per_det[PULSAR_MAX_DETECTORS];
  /// Number of additions to per-detector F-statistics thus far
  UINT4 nsum_twoF_per_det[PULSAR_MAX_DETECTORS];
};

///
/// Internal definition of final results of a semicoherent computation over many segments
///
struct tagWeaveSemiResults {
  /// Per-segment coherent results (optional)
  const WeaveCohResults *coh_res;
  /// Number of per-segment coherent results (optional)
  UINT4 ncoh_res;
  /// Semicoherent template parameters of the first frequency bin
  PulsarDopplerParams semi_phys;
  /// Frequency spacing for semicoherent results
  double dfreq;
  /// Number of frequencies
  UINT4 nfreqs;
  /// Number of detectors for per-detector results
  UINT4 ndetectors;
  /// Mean multi-detector F-statistics per frequency
  REAL4VectorAligned *mean_twoF;
  /// Mean per-detector F-statistics per frequency
  REAL4VectorAligned *mean_twoF_per_det[PULSAR_MAX_DETECTORS];
};

///
/// \name Internal routines
///
/// @{

static int semi_partials_add_REAL4( UINT4 *nsum, REAL4 *sum, const REAL4 *x, const UINT4 nbins );

/// @}

///
/// Create coherent input data
///
WeaveCohInput *XLALWeaveCohInputCreate(
  const BOOLEAN shortcut_compute,
  FstatInput *Fstat_input,
  const LALStringVector *per_detectors
  )
{

  // Allocate memory
  WeaveCohInput *coh_input = XLALCalloc( 1, sizeof( *coh_input ) );
  XLAL_CHECK_NULL( coh_input != NULL, XLAL_ENOMEM );

  // Set fields
  coh_input->shortcut_compute = shortcut_compute;
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

  // Reallocate arrays of multi- and per-detector F-statistics per frequency
  if ( ( *coh_res )->nfreqs < coh_nfreqs ) {
    ( *coh_res )->twoF = XLALRealloc( ( *coh_res )->twoF, coh_nfreqs * sizeof( ( *coh_res )->twoF[0] ) );
    XLAL_CHECK( ( *coh_res )->twoF != NULL, XLAL_ENOMEM );
    for ( size_t i = 0; i < coh_input->Fstat_ndetectors; ++i ) {
      const size_t idx = coh_input->Fstat_res_idx[i];
      ( *coh_res )->twoF_per_det[idx] = XLALRealloc( ( *coh_res )->twoF_per_det[idx], coh_nfreqs * sizeof( ( *coh_res )->twoF_per_det[idx][0] ) );
      XLAL_CHECK( ( *coh_res )->twoF_per_det[idx] != NULL, XLAL_ENOMEM );
    }
  }
  ( *coh_res )->nfreqs = coh_nfreqs;

  // Return now if computation shortcut is in place
  if ( coh_input->shortcut_compute ) {
    return XLAL_SUCCESS;
  }

  // Use a local F-statistic results structure, since we supply our own memory
  // - The 'internalalloclen' field stores the memory size in elements of the results arrays (e.g. 'twoF'),
  //   as opposed to 'numFreqBins' which stores how many elements of the result arrays are in use.
  FstatResults XLAL_INIT_DECL( Fstat_res_struct );
  FstatResults *Fstat_res = &Fstat_res_struct;
  Fstat_res->internalalloclen = ( *coh_res )->nfreqs;
  if ( coh_input->what_to_compute & FSTATQ_2F_PER_DET ) {
    Fstat_res->numDetectors = coh_input->Fstat_ndetectors;
  }
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
  XLAL_CHECK( !( coh_input->what_to_compute & FSTATQ_2F_PER_DET ) || ( Fstat_res->numDetectors == coh_input->Fstat_ndetectors ), XLAL_EFAILED );
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
/// Create and initialise partial semicoherent results
///
int XLALWeaveSemiPartialsInit(
  WeaveSemiPartials **semi_parts,
  const BOOLEAN shortcut_compute,
  const LALStringVector *per_detectors,
  const UINT4 per_nsegments,
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
    if ( per_nsegments > 0 ) {
      ( *semi_parts )->coh_res = XLALCalloc( per_nsegments, sizeof( *( *semi_parts )->coh_res ) );
      XLAL_CHECK( ( *semi_parts )->coh_res != NULL, XLAL_ENOMEM );
    }
  }

  // Set fields
  ( *semi_parts )->shortcut_compute = shortcut_compute;
  ( *semi_parts )->semi_phys = *semi_phys;
  ( *semi_parts )->dfreq = dfreq;
  ( *semi_parts )->nfreqs = semi_nfreqs;
  ( *semi_parts )->ndetectors = ( per_detectors != NULL ) ? per_detectors->length : 0;

  // Initialise number of per-segment coherent results
  ( *semi_parts )->ncoh_res = 0;

  // Initialise number of additions to multi- and per-detector F-statistics
  ( *semi_parts )->nsum_twoF = 0;
  XLAL_INIT_MEM( ( *semi_parts )->nsum_twoF_per_det );

  // Reallocate arrays of summed multi- and per-detector F-statistics per frequency
  if ( ( *semi_parts )->sum_twoF == NULL || ( *semi_parts )->sum_twoF->length < ( *semi_parts )->nfreqs ) {
    XLALDestroyREAL4VectorAligned( ( *semi_parts )->sum_twoF );
    ( *semi_parts )->sum_twoF = XLALCreateREAL4VectorAligned( ( *semi_parts )->nfreqs, alignment );
    XLAL_CHECK( ( *semi_parts )->sum_twoF != NULL, XLAL_ENOMEM );
  }
  for ( size_t i = 0; i < ( *semi_parts )->ndetectors; ++i ) {
    if ( ( *semi_parts )->sum_twoF_per_det[i] == NULL || ( *semi_parts )->sum_twoF_per_det[i]->length < ( *semi_parts )->nfreqs ) {
      XLALDestroyREAL4VectorAligned( ( *semi_parts )->sum_twoF_per_det[i] );
      ( *semi_parts )->sum_twoF_per_det[i] = XLALCreateREAL4VectorAligned( ( *semi_parts )->nfreqs, alignment );
      XLAL_CHECK( ( *semi_parts )->sum_twoF_per_det[i] != NULL, XLAL_ENOMEM );
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
    XLALFree( semi_parts->coh_res );
    XLALDestroyREAL4VectorAligned( semi_parts->sum_twoF );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALDestroyREAL4VectorAligned( semi_parts->sum_twoF_per_det[i] );
    }
    XLALFree( semi_parts );
  }
}

///
/// Add 'nbins' of REAL4 array 'x' to 'sum', and keep track of the number of summations 'nsum'
///
int semi_partials_add_REAL4(
  UINT4 *nsum,
  REAL4 *sum,
  const REAL4 *x,
  const UINT4 nbins
  )
{
  if ( (*nsum)++ == 0 ) {
    // If this is the first summation, just use memcpy()
    memcpy( sum, x, sizeof( *sum ) * nbins );
    return XLAL_SUCCESS;
  }
  return XLALVectorAddREAL4( sum, sum, x, nbins );
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
  XLAL_CHECK( coh_res != NULL, XLAL_EFAULT );

  // Check that offset does not overrun coherent results arrays
  XLAL_CHECK( coh_offset + semi_parts->nfreqs <= coh_res->nfreqs, XLAL_EFAILED, "Coherent offset (%u) + number of semicoherent frequency bins (%u) > number of coherent frequency bins (%u)", coh_offset, semi_parts->nfreqs, coh_res->nfreqs );

  // Store per-segment coherent template parameters and F-statistics per frequency
  if ( semi_parts->coh_res != NULL ) {
    WeaveCohResults *semi_parts_coh_res = &semi_parts->coh_res[semi_parts->ncoh_res++];
    *semi_parts_coh_res = *coh_res;

    // Offset coherent template frequency
    semi_parts_coh_res->coh_phys.fkdot[0] += semi_parts->dfreq * coh_offset;

    // Offset arrays of coherent F-statistics per frequency
    semi_parts_coh_res->twoF += coh_offset;
    for ( size_t i = 0; i < semi_parts->ndetectors; ++i ) {
      if ( semi_parts_coh_res->twoF_per_det[i] != NULL ) {
        semi_parts_coh_res->twoF_per_det[i] += coh_offset;
      }
    }

  }

  // Return now if computation shortcut is in place
  if ( semi_parts->shortcut_compute ) {
    return XLAL_SUCCESS;
  }

  // Add to summed multi-detector F-statistics per frequency, and increment number of additions thus far
  XLAL_CHECK( semi_partials_add_REAL4( &semi_parts->nsum_twoF, semi_parts->sum_twoF->data, coh_res->twoF + coh_offset, semi_parts->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Add to summed per-detector F-statistics per frequency, and increment number of additions thus far
  for ( size_t i = 0; i < semi_parts->ndetectors; ++i ) {
    if ( coh_res->twoF_per_det[i] != NULL ) {
      XLAL_CHECK( semi_partials_add_REAL4( &semi_parts->nsum_twoF_per_det[i], semi_parts->sum_twoF_per_det[i]->data, coh_res->twoF_per_det[i] + coh_offset, semi_parts->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
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

  // Allocate results struct if required
  if ( *semi_res == NULL ) {
    *semi_res = XLALCalloc( 1, sizeof( **semi_res ) );
    XLAL_CHECK( *semi_res != NULL, XLAL_ENOMEM );
  }

  // Set fields
  ( *semi_res )->coh_res = semi_parts->coh_res;
  ( *semi_res )->ncoh_res = semi_parts->ncoh_res;
  ( *semi_res )->semi_phys = semi_parts->semi_phys;
  ( *semi_res )->dfreq = semi_parts->dfreq;
  ( *semi_res )->nfreqs = semi_parts->nfreqs;
  ( *semi_res )->ndetectors = semi_parts->ndetectors;

  // Reallocate arrays of mean multi- and per-detector F-statistics per frequency
  if ( ( *semi_res )->mean_twoF == NULL || ( *semi_res )->mean_twoF->length < semi_parts->nfreqs ) {
    XLALDestroyREAL4VectorAligned( ( *semi_res )->mean_twoF );
    ( *semi_res )->mean_twoF = XLALCreateREAL4VectorAligned( semi_parts->nfreqs, alignment );
    XLAL_CHECK( ( *semi_res )->mean_twoF != NULL, XLAL_ENOMEM );
  }
  for ( size_t i = 0; i < semi_parts->ndetectors; ++i ) {
    if ( ( *semi_res )->mean_twoF_per_det[i] == NULL || ( *semi_res )->mean_twoF_per_det[i]->length < semi_parts->nfreqs ) {
      XLALDestroyREAL4VectorAligned( ( *semi_res )->mean_twoF_per_det[i] );
      ( *semi_res )->mean_twoF_per_det[i] = XLALCreateREAL4VectorAligned( semi_parts->nfreqs, alignment );
      XLAL_CHECK( ( *semi_res )->mean_twoF_per_det[i] != NULL, XLAL_ENOMEM );
    }
  }

  // Return now if computation shortcut is in place
  if ( semi_parts->shortcut_compute ) {
    return XLAL_SUCCESS;
  }

  // Calculate mean multi-detector F-statistics per frequency
  XLAL_CHECK( XLALVectorScaleREAL4( ( *semi_res )->mean_twoF->data, 1.0 / semi_parts->nsum_twoF, semi_parts->sum_twoF->data, ( *semi_res )->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Calculate mean per-detector F-statistics per frequency
  for ( size_t i = 0; i < ( *semi_res )->ndetectors; ++i ) {
    XLAL_CHECK( XLALVectorScaleREAL4( ( *semi_res )->mean_twoF_per_det[i]->data, 1.0 / semi_parts->nsum_twoF_per_det[i], semi_parts->sum_twoF_per_det[i]->data, ( *semi_res )->nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
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
    XLALDestroyREAL4VectorAligned( semi_res->mean_twoF );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALDestroyREAL4VectorAligned( semi_res->mean_twoF_per_det[i] );
    }
    XLALFree( semi_res );
  }
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

  // Fully initialise all output result item field, if requested
  // - Otherwise only output result item fields that change with 'freq_idx' are updated
  if ( *full_init ) {

    // Set all semicoherent and coherent template parameters
    ( *item )->semi_phys = semi_res->semi_phys;
    for ( size_t j = 0; j < semi_res->ncoh_res; ++j ) {
      ( *item )->per_seg[j].coh_phys = semi_res->coh_res[j].coh_phys;
    }

    // Next time, only output result item fields that change with 'freq_idx' should need updating
    *full_init = 0;

  }

  // Update semicoherent and coherent template frequency
  const double foffset = freq_idx * semi_res->dfreq;
  ( *item )->semi_phys.fkdot[0] = semi_res->semi_phys.fkdot[0] + foffset;
  for ( size_t j = 0; j < semi_res->ncoh_res; ++j ) {
    ( *item )->per_seg[j].coh_phys.fkdot[0] = semi_res->coh_res[j].coh_phys.fkdot[0] + foffset;
  }

  // Update multi-detector F-statistics
  ( *item )->mean_twoF = semi_res->mean_twoF->data[freq_idx];
  for ( size_t j = 0; j < semi_res->ncoh_res; ++j ) {
    ( *item )->per_seg[j].twoF = semi_res->coh_res[j].twoF[freq_idx];
  }

  // Update per-detector F-statistics
  for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
    ( *item )->mean_twoF_per_det[i] = semi_res->mean_twoF_per_det[i]->data[freq_idx];
    for ( size_t j = 0; j < semi_res->ncoh_res; ++j ) {
      if ( semi_res->coh_res[j].twoF_per_det[i] != NULL ) {
        ( *item )->per_seg[j].twoF_per_det[i] = semi_res->coh_res[j].twoF_per_det[i][freq_idx];
      } else {
        // There is not per-detector F-statistic for this segment, usually because this segment contains
        // no data from this detector. In this case we output a clearly invalid F-statistic value.
        ( *item )->per_seg[j].twoF_per_det[i] = NAN;
      }
    }
  }

  return XLAL_SUCCESS;

}
