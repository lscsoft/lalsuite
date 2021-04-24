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

///
/// \file
/// \ingroup lalapps_pulsar_Weave
///

#include "SearchIteration.h"

///
/// Iterator over a search parameter space
///
struct tagWeaveSearchIterator {
  /// Number of parameter-space dimensions
  size_t ndim;
  /// Iterator over semicoherent parameter space
  LatticeTilingIterator *semi_itr;
  /// Current lattice tiling point
  gsl_vector *semi_rssky;
  /// Dimension in which to partition iterator output
  int partition_dim;
  /// Maximum number of points in the partitioned dimension
  UINT4 partition_max;
  /// Number of points in each partition
  UINT4 partition_size;
  /// Number of partitions of iterator output
  UINT4 partition_count;
  /// Index of the current iterator output partition
  UINT4 partition_index;
  /// Number of times to repeat iteration over semicoherent parameter space
  UINT4 repetition_count;
  /// Index of the current repetition
  UINT4 repetition_index;
  /// Progress count for iteration
  UINT8 prog_count;
  /// Progress index for iteration
  UINT8 prog_index;
};

///
/// Create iterator over the main loop search parameter space
///
WeaveSearchIterator *XLALWeaveMainLoopSearchIteratorCreate(
  const LatticeTiling *semi_tiling,
  const UINT4 freq_partitions,
  const UINT4 f1dot_partitions
  )
{

  // Check input
  XLAL_CHECK_NULL( semi_tiling != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL( freq_partitions > 0, XLAL_EINVAL );

  // Allocate memory
  WeaveSearchIterator *itr = XLALCalloc( 1, sizeof( *itr ) );
  XLAL_CHECK_NULL( itr != NULL, XLAL_ENOMEM );

  // Get number of parameter-space dimensions
  itr->ndim = XLALTotalLatticeTilingDimensions( semi_tiling );

  // Create iterator over semicoherent tiling
  // - The last parameter-space dimension is always frequency and is not iterated over, since we
  //   always operate over a block of frequencies at once. Since the frequency spacing is always
  //   equal over all tilings due to XLALEqualizeReducedSuperskyMetricsFreqSpacing(), operations
  //   such as nearest point finding can be performed once per frequency block instead of per bin.
  itr->semi_itr = XLALCreateLatticeTilingIterator( semi_tiling, itr->ndim - 1 );

  // Allocate current lattice tiling point
  GAVEC_NULL( itr->semi_rssky, itr->ndim );

  // Partition iterator output in reduced supersky spindown dimension
  {
    itr->partition_count = f1dot_partitions;
    itr->partition_index = 0;
    itr->partition_dim = XLALLatticeTilingDimensionByName( semi_tiling, "nu1dot" );
    XLAL_CHECK_NULL( itr->partition_dim >= 0, XLAL_EFUNC );
    const LatticeTilingStats *stats = XLALLatticeTilingStatistics( semi_tiling, itr->partition_dim );
    XLAL_CHECK_NULL( stats != NULL, XLAL_EFUNC );
    XLAL_CHECK_NULL( stats->max_points > 0, XLAL_ESIZE );
    itr->partition_max = stats->max_points;
    itr->partition_size = 1 + ( itr->partition_max - 1 ) / itr->partition_count;
  }

  // Set repetition count and index
  itr->repetition_count = freq_partitions;
  itr->repetition_index = 0;

  // Set progress count and index for iteration
  itr->prog_count = itr->repetition_count * XLALTotalLatticeTilingPoints( itr->semi_itr );
  XLAL_CHECK_NULL( itr->prog_count > 0, XLAL_EFUNC );
  itr->prog_index = 0;

  return itr;

}

///
/// Destroy iterator
///
void XLALWeaveSearchIteratorDestroy(
  WeaveSearchIterator *itr
  )
{
  if ( itr != NULL ) {
    XLALDestroyLatticeTilingIterator( itr->semi_itr );
    GFVEC( itr->semi_rssky );
    XLALFree( itr );
  }
}

///
/// Save state of iterator to a FITS file
///
int XLALWeaveSearchIteratorSave(
  const WeaveSearchIterator *itr,
  FITSFile *file
  )
{

  // Check input
  XLAL_CHECK( itr != NULL, XLAL_EFAULT );
  XLAL_CHECK( file != NULL, XLAL_EFAULT );

  // Write state of iterator over semicoherent parameter space
  XLAL_CHECK( XLALSaveLatticeTilingIterator( itr->semi_itr, file, "itrstate" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write partition index
  XLAL_CHECK( XLALFITSHeaderWriteUINT4( file, "partidx", itr->partition_index, "partition index" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write repetition index
  XLAL_CHECK( XLALFITSHeaderWriteUINT4( file, "reptidx", itr->repetition_index, "repetition index" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write progress index
  XLAL_CHECK( XLALFITSHeaderWriteUINT8( file, "progidx", itr->prog_index, "progress index" ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Restore state of iterator from a FITS file
///
int XLALWeaveSearchIteratorRestore(
  WeaveSearchIterator *itr,
  FITSFile *file
  )
{

  // Check input
  XLAL_CHECK( itr != NULL, XLAL_EFAULT );
  XLAL_CHECK( file != NULL, XLAL_EFAULT );

  // Read state of iterator over semicoherent parameter space
  XLAL_CHECK( XLALRestoreLatticeTilingIterator( itr->semi_itr, file, "itrstate" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Read partition index
  XLAL_CHECK( XLALFITSHeaderReadUINT4( file, "partidx", &itr->partition_index ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( itr->partition_index < itr->partition_count, XLAL_EIO );

  // Read repetition index
  XLAL_CHECK( XLALFITSHeaderReadUINT4( file, "reptidx", &itr->repetition_index ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( itr->repetition_index < itr->repetition_count, XLAL_EIO );

  // Read progress index
  XLAL_CHECK( XLALFITSHeaderReadUINT8( file, "progidx", &itr->prog_index ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( itr->prog_index < itr->prog_count, XLAL_EIO );

  return XLAL_SUCCESS;

}

///
/// Advance to next state of iterator
///
int XLALWeaveSearchIteratorNext(
  WeaveSearchIterator *itr,
  BOOLEAN *iteration_complete,
  BOOLEAN *expire_cache,
  UINT8 *semi_index,
  const gsl_vector **semi_rssky,
  INT4 *semi_left,
  INT4 *semi_right,
  UINT4 *repetition_index
  )
{

  // Check input
  XLAL_CHECK( itr != NULL, XLAL_EFAULT );
  XLAL_CHECK( iteration_complete != NULL, XLAL_EFAULT );
  XLAL_CHECK( expire_cache != NULL, XLAL_EFAULT );

  // Initialise output flags
  *iteration_complete = 0;
  *expire_cache = 0;

  while ( 1 ) {

    // Get the next frequency block in the semicoherent tiling iterator
    // - XLALNextLatticeTilingPoint() returns mid-point in non-iterated dimensions
    const int itr_retn = XLALNextLatticeTilingPoint( itr->semi_itr, itr->semi_rssky );
    XLAL_CHECK( itr_retn >= 0, XLAL_EFUNC );
    if ( itr_retn == 0 ) {

      // Move to the next partition
      ++itr->partition_index;
      if ( itr->partition_index == itr->partition_count ) {
        itr->partition_index = 0;

        // Move to the next repetition
        ++itr->repetition_index;
        if ( itr->repetition_index == itr->repetition_count ) {

          // Iteration is complete
          *iteration_complete = 1;
          return XLAL_SUCCESS;

        }

      }

      // Expire cache items from previous partitions/repetitions
      *expire_cache = 1;

      // Reset iterator over semicoherent tiling
      XLAL_CHECK( XLALResetLatticeTilingIterator( itr->semi_itr ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLALNextLatticeTilingPoint( itr->semi_itr, itr->semi_rssky ) > 0, XLAL_EFUNC );

    }

    // Partition iterator output, if requested
    if ( itr->partition_count > 1 ) {
      INT4 current_left = 0, current_right = 0;
      XLAL_CHECK( XLALCurrentLatticeTilingBlock( itr->semi_itr, itr->partition_dim, &current_left, &current_right ) == XLAL_SUCCESS, XLAL_EINVAL );
      const UINT4 current_max = current_right - current_left + 1;
      XLAL_CHECK( current_max <= itr->partition_max, XLAL_ESIZE );
      const UINT4 current_index = -current_left / itr->partition_size;
      XLAL_CHECK( current_index < itr->partition_count, XLAL_ESIZE );
      if ( current_index != itr->partition_index ) {
        continue;
      }
    }

    // Next frequency block has been found
    break;

  }

  // Update iteration progress
  ++itr->prog_index;

  // Return iterator state
  *semi_index = XLALCurrentLatticeTilingIndex( itr->semi_itr );
  *semi_rssky = itr->semi_rssky;
  XLAL_CHECK( XLALCurrentLatticeTilingBlock( itr->semi_itr, itr->ndim - 1, semi_left, semi_right ) == XLAL_SUCCESS, XLAL_EFUNC );
  *repetition_index = itr->repetition_index;

  return XLAL_SUCCESS;

}

///
/// Return progress of iterator as a percentage
///
REAL8 XLALWeaveSearchIteratorProgress(
  const WeaveSearchIterator *itr
  )
{

  // Check input
  XLAL_CHECK_REAL8( itr != NULL, XLAL_EFAULT );

  // Return iteration progress
  return GSL_MAX( 0.0, GSL_MIN( ( ( double ) itr->prog_index ) / itr->prog_count, 1.0 ) ) * 100.0;

}

///
/// Return estimate of time remaining for iteration to complete,
/// assuming a equal dstribution in computation cost over time
///
REAL8 XLALWeaveSearchIteratorRemainingTime(
  const WeaveSearchIterator *itr,
  const REAL8 elapsed_time
  )
{

  // Check input
  XLAL_CHECK_REAL8( itr != NULL, XLAL_EFAULT );

  // Return estimate of time remaining for iteration to complete,
  return elapsed_time * ( itr->prog_count - itr->prog_index ) / itr->prog_index;

}

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
