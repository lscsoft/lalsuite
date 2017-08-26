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

#include "CacheResults.h"

#include <lal/LALHeap.h>
#include <lal/LALHashTbl.h>
#include <lal/LALBitset.h>

// Compare two quantities, and return a sort order value if they are unequal
#define COMPARE_BY( x, y ) do { if ( (x) < (y) ) return -1; if ( (x) > (y) ) return +1; } while(0)

///
/// Internal definition of an item stored in the cache
///
typedef struct {
  /// Frequency partition index, used both to find items
  /// in cache and to decide how long to keep items
  UINT4 partition_index;
  /// Relevance, used to decide how long to keep items
  REAL4 relevance;
  /// Coherent locator index, used to find items in cache
  UINT8 coh_index;
  /// Results of a coherent computation on a single segment
  WeaveCohResults *coh_res;
} cache_item;

///
/// Internal definition of cache
///
struct tagWeaveCache {
  /// Number of parameter-space dimensions
  size_t ndim;
  /// Coherent lattice tiling bounding box in dimension 0
  double coh_bound_box_0;
  /// Physical to lattice coordinate transform
  WeavePhysicalToLattice phys_to_latt;
  /// Lattice to physical coordinate transform
  WeaveLatticeToPhysical latt_to_phys;
  /// Coordinate transform data for coherent lattice
  const void *coh_transf_data;
  /// Coordinate transform data for semicoherent lattice
  const void *semi_transf_data;
  /// Input data required for computing coherent results
  WeaveCohInput *coh_input;
  /// Coherent parameter-space tiling locator
  LatticeTilingLocator *coh_locator;
  /// Heap which ranks cache items by relevance
  LALHeap *relevance_heap;
  /// Hash table which looks up cache items by index
  LALHashTbl *coh_index_hash;
  /// Bitset which records whether an item has ever been computed
  LALBitset *coh_computed_bitset;
  /// Partition index for which 'coh_computed_bitset' applies
  UINT4 coh_computed_bitset_partition_index;
  /// Save an no-longer-used cache item for re-use
  cache_item *saved_item;
  /// Cache garbage collection limit
  size_t gc_limit;
};

///
/// Internal definition of results from a series of cache queries
///
struct tagWeaveCacheQueries {
  /// Number of parameter-space dimensions
  size_t ndim;
  /// Semicoherent lattice tiling bounding box in dimension 0
  double semi_bound_box_0;
  /// Physical to lattice coordinate transform
  WeavePhysicalToLattice phys_to_latt;
  /// Lattice to physical coordinate transform
  WeaveLatticeToPhysical latt_to_phys;
  /// Coordinate transform data for semicoherent lattice
  const void *semi_transf_data;
  /// Number of queries for which space is allocated
  UINT4 nqueries;
  /// Sequential indexes for each queried coherent frequency block
  UINT8 *coh_index;
  /// Physical points of each queried coherent frequency block
  PulsarDopplerParams *coh_phys;
  /// Indexes of left-most point in queried coherent frequency block
  INT4 *coh_left;
  /// Indexes of right-most point in queried coherent frequency block
  INT4 *coh_right;
  /// Relevance of each queried coherent frequency block
  REAL4 *coh_relevance;
  /// Number of partitions to divide semicoherent frequency block into
  UINT4 npartitions;
  /// Index to current partition of semicoherent frequency block
  UINT4 partition_index;
  /// Physical coordinates of the current semicoherent frequency block
  PulsarDopplerParams semi_phys;
  /// Index of left-most point in current semicoherent frequency block
  INT4 semi_left;
  /// Index of right-most point in current semicoherent frequency block
  INT4 semi_right;
  /// Relevance of the current semicoherent frequency block
  REAL4 semi_relevance;
  /// Maximum number of bins per partition in semicoherent frequency block
  UINT4 max_semi_part_nfreqs;
};

///
/// \name Internal functions
///
/// @{

static UINT8 cache_item_hash( const void *x );
static int cache_item_compare_by_coh_index( const void *x, const void *y );
static int cache_item_compare_by_relevance( const void *x, const void *y );
static void cache_item_destroy( void *x );

/// @}

///
/// Destroy a cache item
///
void cache_item_destroy(
  void *x
  )
{
  if ( x != NULL ) {
    cache_item *ix = ( cache_item * ) x;
    XLALWeaveCohResultsDestroy( ix->coh_res );
    XLALFree( ix );
  }
} // cache_item_destroy()

///
/// Compare cache items by partition index, then relevance
///
int cache_item_compare_by_relevance(
  const void *x,
  const void *y
  )
{
  const cache_item *ix = ( const cache_item * ) x;
  const cache_item *iy = ( const cache_item * ) y;
  COMPARE_BY( ix->partition_index, iy->partition_index );   // Compare in ascending order
  COMPARE_BY( ix->relevance, iy->relevance );   // Compare in ascending order
  return 0;
} // cache_item_compare_by_relevance()

///
/// Compare cache items by partition index, then locator index
///
int cache_item_compare_by_coh_index(
  const void *x,
  const void *y
  )
{
  const cache_item *ix = ( const cache_item * ) x;
  const cache_item *iy = ( const cache_item * ) y;
  COMPARE_BY( ix->partition_index, iy->partition_index );   // Compare in ascending order
  COMPARE_BY( ix->coh_index, iy->coh_index );   // Compare in ascending order
  return 0;
} // cache_item_compare_by_coh_index()

///
/// Hash cache items by partition index and locator index
///
UINT8 cache_item_hash(
  const void *x
  )
{
  const cache_item *ix = ( const cache_item * ) x;
  UINT4 hval = 0;
  XLALPearsonHash( &hval, sizeof( hval ), &ix->partition_index, sizeof( ix->partition_index ) );
  XLALPearsonHash( &hval, sizeof( hval ), &ix->coh_index, sizeof( ix->coh_index ) );
  return hval;
} // cache_item_hash()

///
/// Create a cache
///
WeaveCache *XLALWeaveCacheCreate(
  const LatticeTiling *coh_tiling,
  const BOOLEAN interpolation,
  const WeavePhysicalToLattice phys_to_latt,
  const WeaveLatticeToPhysical latt_to_phys,
  const void *coh_transf_data,
  const void *semi_transf_data,
  WeaveCohInput *coh_input,
  const size_t max_size,
  const size_t gc_limit
  )
{

  // Check input
  XLAL_CHECK_NULL( coh_tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( phys_to_latt != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( latt_to_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( coh_input != NULL, XLAL_EFAULT );

  // Allocate memory
  WeaveCache *cache = XLALCalloc( 1, sizeof( *cache ) );
  XLAL_CHECK_NULL( cache != NULL, XLAL_ENOMEM );

  // Set fields
  cache->phys_to_latt = phys_to_latt;
  cache->latt_to_phys = latt_to_phys;
  cache->coh_transf_data = coh_transf_data;
  cache->semi_transf_data = semi_transf_data;
  cache->coh_input = coh_input;
  cache->gc_limit = gc_limit;

  // Get number of parameter-space dimensions
  cache->ndim = XLALTotalLatticeTilingDimensions( coh_tiling );
  XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );

  // Get coherent lattice tiling bounding box in dimension 0
  cache->coh_bound_box_0 = XLALLatticeTilingBoundingBox( coh_tiling, 0 );
  XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );

  // If this is an interpolating search, create a lattice tiling locator
  if ( interpolation ) {
    cache->coh_locator = XLALCreateLatticeTilingLocator( coh_tiling );
    XLAL_CHECK_NULL( cache->coh_locator != NULL, XLAL_EFUNC );
  }

  // Create a heap which sorts items by "relevance", a quantity which determines how long
  // cache items are kept. Consider the following scenario:
  //
  //   +-----> parameter-space dimension 0
  //   |
  //   V other parameter-space dimensions
  //
  //        :
  //        : R[S1] = relevance of semicoherent point 'S1'
  //        :
  //        +-----+
  //        | /`\ |
  //        || S1||
  //        | \,/ |  :
  //        +-----+  : R[S2] = relevance of semicoherent point 'S2'
  //   +-------+     :
  //   | /```\ |     +-----+
  //   |/     \|     | /`\ |
  //   ||  C  ||     || S2||
  //   |\     /|     | \,/ |
  //   | \,,,/ |     +-----+
  //   +-------+
  //           :
  //           : R[C] = relevance of coherent point 'C'
  //           :
  //
  // The relevance R[C] of the coherent point 'C' is given by the coordinate in dimension 0 of
  // the *rightmost* edge of the bounding box surrounding its metric ellipse. The relevances of
  // two semicoherent points S1 and S2, R[S1] and R[S2], are given by the *leftmost* edges of the
  // bounding box surrounding their metric ellipses.
  //
  // Note that iteration over the parameter space is ordered such that dimension 0 is the slowest
  // coordinate, i.e. dimension 0 is passed over only once, therefore coordinates in this dimension
  // are monotonically increasing, therefore relevances are also monotonically increasing.
  //
  // Suppose S1 is the current point in the semicoherent parameter-space tiling; note that R[C] > R[S1].
  // It is clear that, as iteration progresses, some future semicoherent points will overlap with C.
  // Therefore C cannot be discarded from the cache, since it will be the closest point for future
  // semicoherent points. Now suppose that S2 is the current point; note that R[C] < R[S1]. It is clear
  // that neither S2, nor any future point in the semicoherent parameter-space tiling, can ever overlap
  // with C. Therefore C can never be the closest point for any future semicoherent points, and therefore
  // it can safely be discarded from the cache.
  //
  // In short, an item in the cache can be discarded once its relevance falls below the threshold set
  // by the current point semicoherent in the semicoherent parameter-space tiling.
  //
  // Items removed from the heap are destroyed by calling cache_item_destroy().
  cache->relevance_heap = XLALHeapCreate( cache_item_destroy, max_size, -1, cache_item_compare_by_relevance );
  XLAL_CHECK_NULL( cache->relevance_heap != NULL, XLAL_EFUNC );

  // Create a hash table which looks up cache items by partition and locator index. Items removed
  // from the hash table are NOT destroyed, since items are shared with 'relevance_heap'.
  cache->coh_index_hash = XLALHashTblCreate( NULL, cache_item_hash, cache_item_compare_by_coh_index );
  XLAL_CHECK_NULL( cache->coh_index_hash != NULL, XLAL_EFUNC );

  // Create a bitset which records which cache items have ever been computed.
  cache->coh_computed_bitset = XLALBitsetCreate();
  XLAL_CHECK_NULL( cache->coh_computed_bitset != NULL, XLAL_EFUNC );
  cache->coh_computed_bitset_partition_index = LAL_UINT4_MAX;

  return cache;

} // XLALWeaveCacheCreate()

///
/// Destroy a cache
///
void XLALWeaveCacheDestroy(
  WeaveCache *cache
  )
{
  if ( cache != NULL ) {
    XLALDestroyLatticeTilingLocator( cache->coh_locator );
    XLALHeapDestroy( cache->relevance_heap );
    XLALHashTblDestroy( cache->coh_index_hash );
    cache_item_destroy( cache->saved_item );
    XLALBitsetDestroy( cache->coh_computed_bitset );
    XLALFree( cache );
  }
} // XLALWeaveCacheDestroy()

///
/// Create storage for a series of cache queries
///
WeaveCacheQueries *XLALWeaveCacheQueriesCreate(
  const LatticeTiling *semi_tiling,
  const WeavePhysicalToLattice phys_to_latt,
  const WeaveLatticeToPhysical latt_to_phys,
  const void *semi_transf_data,
  const UINT4 nqueries,
  const UINT4 npartitions
  )
{

  // Check input
  XLAL_CHECK_NULL( semi_tiling != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( phys_to_latt != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( latt_to_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( nqueries > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( npartitions > 0, XLAL_EINVAL );

  // Allocate memory
  WeaveCacheQueries *queries = XLALCalloc( 1, sizeof( *queries ) );
  XLAL_CHECK_NULL( queries != NULL, XLAL_ENOMEM );
  queries->coh_index = XLALCalloc( nqueries, sizeof( *queries->coh_index ) );
  XLAL_CHECK_NULL( queries->coh_index != NULL, XLAL_ENOMEM );
  queries->coh_phys = XLALCalloc( nqueries, sizeof( *queries->coh_phys ) );
  XLAL_CHECK_NULL( queries->coh_phys != NULL, XLAL_ENOMEM );
  queries->coh_left = XLALCalloc( nqueries, sizeof( *queries->coh_left ) );
  XLAL_CHECK_NULL( queries->coh_left != NULL, XLAL_ENOMEM );
  queries->coh_right = XLALCalloc( nqueries, sizeof( *queries->coh_right ) );
  XLAL_CHECK_NULL( queries->coh_right != NULL, XLAL_ENOMEM );
  queries->coh_relevance = XLALCalloc( nqueries, sizeof( *queries->coh_relevance ) );
  XLAL_CHECK_NULL( queries->coh_relevance != NULL, XLAL_ENOMEM );

  // Set fields
  queries->phys_to_latt = phys_to_latt;
  queries->latt_to_phys = latt_to_phys;
  queries->semi_transf_data = semi_transf_data;
  queries->nqueries = nqueries;
  queries->npartitions = npartitions;

  // Get number of parameter-space dimensions
  queries->ndim = XLALTotalLatticeTilingDimensions( semi_tiling );
  XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );

  // Get semicoherent lattice tiling bounding box in dimension 0
  queries->semi_bound_box_0 = XLALLatticeTilingBoundingBox( semi_tiling, 0 );
  XLAL_CHECK_NULL( xlalErrno == 0, XLAL_EFUNC );

  // Get maximum number of bins (per partition) in semicoherent frequency block
  {
    const LatticeTilingStats *semi_freq_stats = XLALLatticeTilingStatistics( semi_tiling, queries->ndim - 1 );
    XLAL_CHECK_NULL( semi_freq_stats != NULL, XLAL_EFUNC );
    XLAL_CHECK_NULL( semi_freq_stats->max_points > 0, XLAL_EFAILED );
    queries->max_semi_part_nfreqs = 1 + ( ( semi_freq_stats->max_points - 1 ) / queries->npartitions );
  }

  return queries;

} // XLALWeaveCacheQueriesCreate()

///
/// Destroy storage for a series of cache queries
///
void XLALWeaveCacheQueriesDestroy(
  WeaveCacheQueries *queries
  )
{
  if ( queries != NULL ) {
    XLALFree( queries->coh_index );
    XLALFree( queries->coh_phys );
    XLALFree( queries->coh_left );
    XLALFree( queries->coh_right );
    XLALFree( queries->coh_relevance );
    XLALFree( queries );
  }
} // XLALWeaveCacheQueriesDestroy()

///
/// Initialise a series of cache queries
///
int XLALWeaveCacheQueriesInit(
  WeaveCacheQueries *queries,
  const LatticeTilingIterator *semi_itr,
  const gsl_vector *semi_point
  )
{

  // Check input
  XLAL_CHECK( queries != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_itr != NULL, XLAL_EFAULT );

  // Reset sequential indexes to zero, indicating no query
  memset( queries->coh_index, 0, queries->nqueries * sizeof( queries->coh_index[0] ) );

  // Compute the relevance of the current semicoherent frequency block
  // - Subtract half the semicoherent lattice tiling bounding box in dimension 0
  queries->semi_relevance = gsl_vector_get( semi_point, 0 ) - 0.5 * queries->semi_bound_box_0;

  // Convert semicoherent point to physical coordinates
  XLAL_CHECK( ( queries->latt_to_phys )( &queries->semi_phys, semi_point, NULL, queries->semi_transf_data ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Get indexes of left/right-most point in current semicoherent frequency block
  XLAL_CHECK( XLALCurrentLatticeTilingBlock( semi_itr, queries->ndim - 1, &queries->semi_left, &queries->semi_right ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // XLALWeaveCacheQueriesInit()

///
/// Query a cache for the results nearest to a given coherent point
///
int XLALWeaveCacheQuery(
  const WeaveCache *cache,
  const UINT8 semi_index,
  WeaveCacheQueries *queries,
  const UINT4 query_index
  )
{

  // Check input
  XLAL_CHECK( cache != NULL, XLAL_EFAULT );
  XLAL_CHECK( queries != NULL, XLAL_EFAULT );
  XLAL_CHECK( query_index < queries->nqueries, XLAL_EINVAL );

  // Convert semicoherent physical point to coherent lattice tiling coordinates
  double coh_point_array[cache->ndim];
  gsl_vector_view coh_point_view = gsl_vector_view_array( coh_point_array, cache->ndim );
  XLAL_CHECK( ( cache->phys_to_latt )( &coh_point_view.vector, &queries->semi_phys, cache->coh_transf_data ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Initialise nearest point to semicoherent point in coherent lattice tiling coordinates
  double coh_near_point_array[cache->ndim];
  gsl_vector_view coh_near_point_view = gsl_vector_view_array( coh_near_point_array, cache->ndim );
  gsl_vector_memcpy( &coh_near_point_view.vector, &coh_point_view.vector );

  // Initialise values for a non-interpolating search
  queries->coh_index[query_index] = semi_index;
  queries->coh_left[query_index] = queries->semi_left;
  queries->coh_right[query_index] = queries->semi_right;

  // If performing an interpolating search, find the nearest coherent frequency pass in this segment
  // - 'coh_near_point' is set to the nearest point to the mid-point of the semicoherent frequency block
  // - 'coh_index' is set to the index of this coherent frequency block, used for cache lookup
  // - 'coh_left' and 'coh_right' are set of number of points to the left/right of 'coh_point'
  if ( cache->coh_locator != NULL ) {
    XLAL_CHECK( XLALNearestLatticeTilingBlock( cache->coh_locator, &coh_point_view.vector, cache->ndim - 1, &coh_near_point_view.vector, &queries->coh_index[query_index], &queries->coh_left[query_index], &queries->coh_right[query_index] ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( queries->coh_left[query_index] <= queries->coh_right[query_index], XLAL_EINVAL );
  }

  // Make 'coh_index' a 1-based index, so zero can be used to indicate a missing query
  ++queries->coh_index[query_index];

  // Convert nearest coherent point to physical coordinates
  XLAL_INIT_MEM( queries->coh_phys[query_index] );
  XLAL_CHECK( ( cache->latt_to_phys )( &queries->coh_phys[query_index], &coh_near_point_view.vector, &coh_point_view.vector, cache->coh_transf_data ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Compute the relevance of the current coherent frequency block
  // - Add half the coherent lattice tiling bounding box in dimension 0
  // - Convert coherent point to semicoherent lattice tiling coordinates
  // - Relevance is semicoherent point coordinate in dimension 0
  {
    double tmp_point_array[cache->ndim];
    gsl_vector_view tmp_point_view = gsl_vector_view_array( tmp_point_array, cache->ndim );
    gsl_vector_memcpy( &tmp_point_view.vector, &coh_near_point_view.vector );
    *gsl_vector_ptr( &tmp_point_view.vector, 0 ) += 0.5 * cache->coh_bound_box_0;
    PulsarDopplerParams XLAL_INIT_DECL( tmp_phys );
    XLAL_CHECK( ( cache->latt_to_phys )( &tmp_phys, &tmp_point_view.vector, NULL, cache->coh_transf_data ) == XLAL_SUCCESS, XLAL_EINVAL );
    XLAL_CHECK( ( cache->phys_to_latt )( &tmp_point_view.vector, &tmp_phys, cache->semi_transf_data ) == XLAL_SUCCESS, XLAL_EFUNC );
    queries->coh_relevance[query_index] = gsl_vector_get( &tmp_point_view.vector, 0 );
  }

  return XLAL_SUCCESS;

} // XLALWeaveCacheQuery()

///
/// Finalise a series of cache queries
///
int XLALWeaveCacheQueriesFinal(
  WeaveCacheQueries *queries,
  const UINT4 partition_index,
  PulsarDopplerParams *semi_phys,
  const double dfreq,
  UINT4 *semi_nfreqs
  )
{

  // Check input
  XLAL_CHECK( queries != NULL, XLAL_EFAULT );
  XLAL_CHECK( partition_index < queries->npartitions, XLAL_EINVAL );
  XLAL_CHECK( semi_phys != NULL, XLAL_EFAULT );
  XLAL_CHECK( dfreq >= 0, XLAL_EINVAL );
  XLAL_CHECK( semi_nfreqs != NULL, XLAL_EFAULT );

  // Save index to current partition of semicoherent frequency block
  queries->partition_index = partition_index;

  for ( size_t i = 0; i < queries->nqueries; ++i ) {

    // Check that 'coh_index' is at least 1, i.e. a query was made
    XLAL_CHECK( queries->coh_index[i] > 0, XLAL_EINVAL, "Missing query at index %zu", i );

    // Check that the semicoherent left/right-most indexes do not exceed the corresponding
    // coherent left/right-most indexes, i.e. that the semicoherent parameter space is a
    // subset of the coherent parameter spaces. This check is for safety as this should be
    // taken care of by additional padding of the lattice tiling parameter spaces.
    queries->semi_left = GSL_MAX( queries->semi_left, queries->coh_left[i] );
    queries->semi_right = GSL_MIN( queries->semi_right, queries->coh_right[i] );

  }

  // Compute the total number of points in the semicoherent frequency block
  *semi_nfreqs = queries->semi_right - queries->semi_left + 1;

  // Adjust semicoherent left/right-most indexes to select the given partition.
  // Shrinks the semicoherent frequency block from '*semi_nfreqs' points to
  // (at most) 'queries->max_semi_part_nfreqs' points.
  const INT4 semi_left_offset = queries->max_semi_part_nfreqs * queries->partition_index;
  const INT4 semi_right_offset = GSL_MIN( 0, semi_left_offset + ( ( INT4 ) queries->max_semi_part_nfreqs ) - ( ( INT4 ) *semi_nfreqs ) );
  queries->semi_left += semi_left_offset;
  queries->semi_right += semi_right_offset;

  // Recompute the number of points in the selected semicoherent frequency block partition.
  // This could now potentially be zero, as some semicoherent frequency blocks may not
  // have enough points for all partitions - in which case skip computing results for the
  // selected partition.
  *semi_nfreqs = GSL_MAX( 0, queries->semi_right - queries->semi_left + 1 );
  if ( *semi_nfreqs == 0 ) {
    return XLAL_SUCCESS;
  }

  // Adjust coherent left/right-most indexes to enclose the given partition.
  // Shrinks each coherent frequency block by an amount independent of the current
  // semicoherent frequency block, since coherent frequency blocks are likely to be
  // queried by multiple semicoherent frequency blocks and therefore need to have
  // enough points to satisfy each query.
  const INT4 coh_left_offset = semi_left_offset;
  const INT4 coh_right_offset = queries->max_semi_part_nfreqs * ( queries->npartitions - queries->partition_index - 1 );
  XLAL_CHECK( coh_left_offset <= semi_left_offset, XLAL_EFAILED, "Adjustment to coherent left-most index (%i) exceeds adjustment to semicoherent left-must index (%i)", coh_left_offset, semi_left_offset );
  XLAL_CHECK( semi_right_offset <= coh_right_offset, XLAL_EFAILED, "Adjustment to semicoherent right-most index (%i) exceeds adjustment to coherent right-must index (%i)", semi_right_offset, coh_right_offset );
  for ( size_t i = 0; i < queries->nqueries; ++i ) {
    queries->coh_left[i] += coh_left_offset;
    queries->coh_right[i] += coh_right_offset;
  }

  // Shift physical frequencies to first point in coherent/semicoherent frequency blocks
  for ( size_t i = 0; i < queries->nqueries; ++i ) {
    queries->coh_phys[i].fkdot[0] += dfreq * queries->coh_left[i];
  }
  queries->semi_phys.fkdot[0] += dfreq * queries->semi_left;

  // Return first point in semicoherent frequency block
  *semi_phys = queries->semi_phys;

  return XLAL_SUCCESS;

} // XLALWeaveCacheQueriesFinal()

///
/// Retrieve coherent results for a given query, or compute new coherent results if not found
///
int XLALWeaveCacheRetrieve(
  WeaveCache *cache,
  WeaveCacheQueries *queries,
  const UINT4 query_index,
  const WeaveCohResults **coh_res,
  UINT8 *coh_index,
  UINT4 *coh_offset,
  UINT8 *coh_nres,
  UINT8 *coh_ntmpl
  )
{

  // Check input
  XLAL_CHECK( cache != NULL, XLAL_EFAULT );
  XLAL_CHECK( queries != NULL, XLAL_EFAULT );
  XLAL_CHECK( queries != NULL, XLAL_EFAULT );
  XLAL_CHECK( query_index < queries->nqueries, XLAL_EINVAL );
  XLAL_CHECK( coh_res != NULL, XLAL_EFAULT );
  XLAL_CHECK( coh_offset != NULL, XLAL_EFAULT );
  XLAL_CHECK( coh_nres != NULL, XLAL_EFAULT );
  XLAL_CHECK( coh_ntmpl != NULL, XLAL_EFAULT );

  // Reset bitset if partition index has changed
  if ( cache->coh_computed_bitset_partition_index != queries->partition_index ) {
    XLAL_CHECK( XLALBitsetClear( cache->coh_computed_bitset ) == XLAL_SUCCESS, XLAL_EFUNC );
    cache->coh_computed_bitset_partition_index = queries->partition_index;
  }

  // See if coherent results are already cached
  const cache_item find_key = { .partition_index = queries->partition_index, .coh_index = queries->coh_index[query_index] };
  const cache_item *find_item = NULL;
  XLAL_CHECK( XLALHashTblFind( cache->coh_index_hash, &find_key, ( const void ** ) &find_item ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( find_item == NULL ) {

    // Reuse 'saved_item' if possible, otherwise allocate memory for a new cache item
    if ( cache->saved_item == NULL ) {
      cache->saved_item = XLALCalloc( 1, sizeof( *cache->saved_item ) );
      XLAL_CHECK( cache->saved_item != NULL, XLAL_EINVAL );
    }
    cache_item *new_item = cache->saved_item;
    find_item = new_item;

    // Set the key of the new cache item for future lookups
    new_item->partition_index = find_key.partition_index;
    new_item->coh_index = find_key.coh_index;

    // Set the relevance of the coherent frequency block associated with the new cache item
    new_item->relevance = queries->coh_relevance[query_index];

    // Determine the number of points in the coherent frequency block
    const UINT4 coh_nfreqs = queries->coh_right[query_index] - queries->coh_left[query_index] + 1;

    // Compute coherent results for the new cache item
    XLAL_CHECK( XLALWeaveCohResultsCompute( &new_item->coh_res, cache->coh_input, &queries->coh_phys[query_index], coh_nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Add new cache item to the index hash table
    XLAL_CHECK( XLALHashTblAdd( cache->coh_index_hash, new_item ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Get the item in the cache with the smallest relevance
    const cache_item *least_relevant_item = ( const cache_item * ) XLALHeapRoot( cache->relevance_heap );
    XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

    // Create a 'fake' item specifying thresholds for cache item relevance, for comparison with 'least_relevant_item'
    const cache_item relevance_threshold = { .partition_index = queries->partition_index, .relevance = queries->semi_relevance };

    // If item's relevance has fallen below the threshold relevance, it can be removed from the cache
    if ( least_relevant_item != NULL && least_relevant_item != new_item && cache_item_compare_by_relevance( least_relevant_item, &relevance_threshold ) < 0 ) {

      // Remove least relevant item from index hash table
      XLAL_CHECK( XLALHashTblRemove( cache->coh_index_hash, least_relevant_item ) == XLAL_SUCCESS, XLAL_EFUNC );

      // Exchange 'saved_item' with the least relevant item in the relevance heap
      XLAL_CHECK( XLALHeapExchangeRoot( cache->relevance_heap, ( void ** ) &cache->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );

      // Remove any addition items that are no longer relevant, up to a maximum of 'gc_limit'
      for ( size_t i = 0; i < cache->gc_limit; ++i ) {

        // Get the item in the cache with the smallest relevance
        least_relevant_item = ( const cache_item * ) XLALHeapRoot( cache->relevance_heap );
        XLAL_CHECK( xlalErrno == 0, XLAL_EFUNC );

        // If item's relevance has fallen below the threshold relevance, it can be removed from the cache
        if ( least_relevant_item != NULL && least_relevant_item != new_item && cache_item_compare_by_relevance( least_relevant_item, &relevance_threshold ) < 0 ) {

          // Remove least relevant item from index hash table
          XLAL_CHECK( XLALHashTblRemove( cache->coh_index_hash, least_relevant_item ) == XLAL_SUCCESS, XLAL_EFUNC );

          // Remove and destroy least relevant item from the relevance heap
          XLAL_CHECK( XLALHeapRemoveRoot( cache->relevance_heap ) == XLAL_SUCCESS, XLAL_EINVAL );

        } else {

          // All cache items are still relevant
          break;

        }

      }

    } else {

      // Add new cache item to the relevance heap; 'saved_item' many now contains an item removed from the heap
      XLAL_CHECK( XLALHeapAdd( cache->relevance_heap, ( void ** ) &cache->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );

      // If 'saved_item' contains an item removed from the heap, also remove it from the index hash table
      if ( cache->saved_item != NULL ) {
        XLAL_CHECK( XLALHashTblRemove( cache->coh_index_hash, cache->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );
      }

    }

    // Increment number of computed coherent results
    *coh_nres += coh_nfreqs;

    // Check if coherent results have been computed previously
    BOOLEAN computed = 0;
    XLAL_CHECK( XLALBitsetGet( cache->coh_computed_bitset, find_key.coh_index, &computed ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( !computed ) {

      // Coherent results have not been computed before: increment the number of coherent templates
      *coh_ntmpl += coh_nfreqs;

      // This coherent result has now been computed
      XLAL_CHECK( XLALBitsetSet( cache->coh_computed_bitset, find_key.coh_index, 1 ) == XLAL_SUCCESS, XLAL_EFUNC );

    }

  }

  // Return coherent results from cache
  *coh_res = find_item->coh_res;

  // Return index of coherent result
  *coh_index = find_item->coh_index;

  // Return offset at which coherent results should be combined with semicoherent results
  *coh_offset = queries->semi_left - queries->coh_left[query_index];

  return XLAL_SUCCESS;

} // XLALWeaveCacheRetrieve()

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
