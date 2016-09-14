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

#include "OutputResults.h"

#include <lal/LALHeap.h>

///
/// Internal definition of output results from a search
///
struct tagWeaveOutputResults {
  /// Reference time at which search is conducted
  LIGOTimeGPS ref_time;
  /// Number of spindown parameters to output
  size_t nspins;
  /// If outputting per-detector quantities, list of detectors
  LALStringVector *per_detectors;
  /// Number of per-segment items to output (may be zero)
  UINT4 per_nsegments;
  /// Total number of semicoherent results added to output
  INT8 semi_total;
  /// Toplist ranked by mean multi-detector F-statistic
  LALHeap *mean_twoF_toplist;
  /// Save a no-longer-used toplist item for re-use
  WeaveOutputToplistItem *saved_item;
};

///
/// \name Internal routines
///
/// @{

static LALHeap *toplist_create( int toplist_limit, LALHeapCmpFcn toplist_item_compare_fcn );
static WeaveOutputToplistItem *toplist_item_create( const LIGOTimeGPS *ref_time, const UINT4 per_nsegments );
static int toplist_compare( BOOLEAN *equal, const WeaveSetupData *setup, const REAL8 param_tol_mism, const VectorComparison *result_tol, const LALStringVector *detectors, const size_t nsegments, const LALHeap *toplist_1, const LALHeap *toplist_2 );
static int toplist_compare_results( BOOLEAN *equal, const VectorComparison *result_tol, const REAL4Vector *res_1, const REAL4Vector *res_2 );
static int toplist_compare_templates( BOOLEAN *equal, const char *loc_str, const char *tmpl_str, const REAL8 param_tol_mism, const WeavePhysicalToLattice phys_to_latt, const gsl_matrix *metric, const void *transf_data, const PulsarDopplerParams *phys_1, const PulsarDopplerParams *phys_2 );
static int toplist_fits_table_init( FITSFile *file, const size_t nspins, const LALStringVector *per_detectors, const UINT4 per_nsegments );
static int toplist_fits_table_read( FITSFile *file, const char *name, WeaveOutputResults *out, LALHeap **toplist, LALHeapCmpFcn toplist_item_compare_fcn );
static int toplist_fits_table_write( FITSFile *file, const char *name, const char *comment, const WeaveOutputResults *out, LALHeap *toplist );
static int toplist_fits_table_write_visitor( void *param, const void *x );
static int toplist_item_add( BOOLEAN *full_init, WeaveOutputResults *out, LALHeap *toplist, const WeaveSemiResults *semi_res, const size_t freq_idx );
static int toplist_item_compare_by_mean_twoF( const void *x, const void *y );
static int toplist_item_sort_by_semi_phys( const void *x, const void *y );
static void toplist_item_destroy( void *x );

/// @}

///
/// Create a toplist
///
LALHeap *toplist_create(
  int toplist_limit,
  LALHeapCmpFcn toplist_item_compare_fcn
  )
{
  LALHeap *toplist = XLALHeapCreate( toplist_item_destroy, toplist_limit, +1, toplist_item_compare_fcn );
  XLAL_CHECK_NULL( toplist != NULL, XLAL_EFUNC );
  return toplist;
}

///
/// Initialise a FITS table for writing/reading a toplist
///
int toplist_fits_table_init(
  FITSFile *file,
  const size_t nspins,
  const LALStringVector *per_detectors,
  const UINT4 per_nsegments
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( nspins > 0, XLAL_EINVAL );

  char col_name[32];

  // Begin FITS table description
  XLAL_FITS_TABLE_COLUMN_BEGIN( WeaveOutputToplistItem );

  // Add columns for semicoherent template parameters
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_phys.Alpha, "alpha [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_phys.Delta, "delta [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_phys.fkdot[0], "freq [Hz]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  for ( size_t k = 1; k <= nspins; ++k ) {
    snprintf( col_name, sizeof( col_name ), "f%zudot [Hz/s^%zu]", k, k );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_phys.fkdot[k], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Add columns for mean multi- and per-detector F-statistic
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, mean_twoF ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( per_detectors != NULL ) {
    for ( size_t i = 0; i < per_detectors->length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "mean_twoF_%s", per_detectors->data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL4, mean_twoF_per_det[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Begin FITS table description for per-segment items (optional)
  if ( per_nsegments > 0 ) {
    XLAL_FITS_TABLE_COLUMN_PTR_BEGIN( per_seg, WeaveOutputToplistPerSegItem, per_nsegments );
    for ( size_t s = 0; s < per_nsegments; ++s ) {

      // Add columns for coherent template parameters
      snprintf( col_name, sizeof( col_name ), "seg%zu_alpha [rad]", s + 1 );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_phys.Alpha, col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      snprintf( col_name, sizeof( col_name ), "seg%zu_delta [rad]", s + 1 );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_phys.Delta, col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      snprintf( col_name, sizeof( col_name ), "seg%zu_freq [Hz]", s + 1 );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_phys.fkdot[0], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      for ( size_t k = 1; k <= nspins; ++k ) {
        snprintf( col_name, sizeof( col_name ), "seg%zu_f%zudot [Hz/s^%zu]", s + 1, k, k );
        XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_phys.fkdot[k], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      // Add columns for coherent multi- and per-detector F-statistic
      snprintf( col_name, sizeof( col_name ), "seg%zu_twoF", s + 1 );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL4, twoF, col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      if ( per_detectors != NULL ) {
        for ( size_t i = 0; i < per_detectors->length; ++i ) {
          snprintf( col_name, sizeof( col_name ), "seg%zu_twoF_%s", s + 1, per_detectors->data[i] );
          XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL4, twoF_per_det[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
      }

    }
  }

  return XLAL_SUCCESS;

}

///
/// Visitor function for writing a toplist to a FITS table
///
int toplist_fits_table_write_visitor(
  void *param,
  const void *x
  )
{
  FITSFile *file = ( FITSFile * ) param;
  XLAL_CHECK( XLALFITSTableWriteRow( file, x ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;
}

///
/// Write toplist items to a FITS table
///
int toplist_fits_table_write(
  FITSFile *file,
  const char *name,
  const char *comment,
  const WeaveOutputResults *out,
  LALHeap *toplist
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( name != NULL, XLAL_EFAULT );
  XLAL_CHECK( comment != NULL, XLAL_EFAULT );
  XLAL_CHECK( out != NULL, XLAL_EFAULT );
  XLAL_CHECK( toplist != NULL, XLAL_EFAULT );

  // Open FITS table for writing and initialise
  XLAL_CHECK( XLALFITSTableOpenWrite( file, name, comment ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( toplist_fits_table_init( file, out->nspins, out->per_detectors, out->per_nsegments ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write all items to FITS table
  XLAL_CHECK( XLALHeapVisit( toplist, toplist_fits_table_write_visitor, file ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write maximum size of toplist to FITS header
  XLAL_CHECK( XLALFITSHeaderWriteINT8( file, "toplimit", XLALHeapMaxSize( toplist ), "maximum size of toplist" ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Read items from a FITS table and either create, or append to existing, toplist
///
int toplist_fits_table_read(
  FITSFile *file,
  const char *name,
  WeaveOutputResults *out,
  LALHeap **toplist,
  LALHeapCmpFcn toplist_item_compare_fcn
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( name != NULL, XLAL_EFAULT );
  XLAL_CHECK( out != NULL, XLAL_EFAULT );
  XLAL_CHECK( toplist != NULL, XLAL_EFAULT );

  // Decide whether to create, or append to existing, toplist
  const BOOLEAN create = ( *toplist == NULL );

  // Open FITS table for reading and initialise
  UINT8 nrows = 0;
  XLAL_CHECK( XLALFITSTableOpenRead( file, name, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( toplist_fits_table_init( file, out->nspins, out->per_detectors, out->per_nsegments ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Read maximum size of toplist to FITS header
  INT8 toplist_limit = 0;
  XLAL_CHECK( XLALFITSHeaderReadINT8( file, "toplimit", &toplist_limit ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Create a new toplist if required
  if ( create ) {
    *toplist = toplist_create( toplist_limit, toplist_item_compare_fcn );
    XLAL_CHECK( *toplist != NULL, XLAL_EFUNC );
  }

  // Read all items from FITS table
  while ( nrows > 0 ) {

    // Create a new toplist item if needed
    if ( out->saved_item == NULL ) {
      out->saved_item = toplist_item_create( &out->ref_time, out->per_nsegments );
      XLAL_CHECK( out->saved_item != NULL, XLAL_ENOMEM );
    }

    // Read item from FITS table
    XLAL_CHECK( XLALFITSTableReadRow( file, out->saved_item, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Add item to toplist
    XLAL_CHECK( XLALHeapAdd( *toplist, ( void ** ) &out->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );

  }

  return XLAL_SUCCESS;

}

///
/// Compute two template parameters
///
int toplist_compare_templates(
  BOOLEAN *equal,
  const char *loc_str,
  const char *tmpl_str,
  const REAL8 param_tol_mism,
  const WeavePhysicalToLattice phys_to_latt,
  const gsl_matrix *metric,
  const void *transf_data,
  const PulsarDopplerParams *phys_1,
  const PulsarDopplerParams *phys_2
  )
{

  // Check input
  XLAL_CHECK( equal != NULL, XLAL_EINVAL );
  XLAL_CHECK( loc_str != NULL, XLAL_EINVAL );
  XLAL_CHECK( tmpl_str != NULL, XLAL_EINVAL );
  XLAL_CHECK( param_tol_mism > 0, XLAL_EINVAL );
  XLAL_CHECK( phys_to_latt != NULL, XLAL_EFAULT );
  XLAL_CHECK( metric != NULL, XLAL_EFAULT );
  XLAL_CHECK( transf_data != NULL, XLAL_EFAULT );
  XLAL_CHECK( phys_1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( phys_2 != NULL, XLAL_EFAULT );

  // Transform physical point to lattice coordinates
  double latt_1_array[metric->size1];
  gsl_vector_view latt_1_view = gsl_vector_view_array( latt_1_array, metric->size1 );
  gsl_vector *const latt_1 = &latt_1_view.vector;
  XLAL_CHECK( ( phys_to_latt )( latt_1, phys_1, transf_data ) == XLAL_SUCCESS, XLAL_EFUNC );
  double latt_2_array[metric->size1];
  gsl_vector_view latt_2_view = gsl_vector_view_array( latt_2_array, metric->size1 );
  gsl_vector *const latt_2 = &latt_2_view.vector;
  XLAL_CHECK( ( phys_to_latt )( latt_2, phys_2, transf_data ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Store difference between lattice coordinates in 'u'
  double u_array[metric->size1];
  gsl_vector_view u_view = gsl_vector_view_array( u_array, metric->size1 );
  gsl_vector *const u = &u_view.vector;
  gsl_vector_memcpy( u, latt_1 );
  gsl_vector_sub( u, latt_2 );

  // Multiply 'u' by metric, storing result in 'v'
  double v_array[metric->size1];
  gsl_vector_view v_view = gsl_vector_view_array( v_array, metric->size1 );
  gsl_vector *const v = &v_view.vector;
  gsl_blas_dsymv( CblasUpper, 1.0, metric, u, 0.0, v );

  // Compute mismatch and compare to tolerance
  REAL8 mism = 0;
  gsl_blas_ddot( u, v, &mism );

  // If mismatch is above tolerance, print error message
  if ( mism > param_tol_mism ) {
    *equal = 0;
    XLALPrintInfo( "%s: at %s, mismatch between %s template parameters exceeds tolerance: %g > %g\n", __func__, loc_str, tmpl_str, mism, param_tol_mism );
    const PulsarDopplerParams *phys[2] = { phys_1, phys_2 };
    gsl_vector *latt[2] = { latt_1, latt_2 };
    for ( size_t i = 0; i < 2; ++i ) {
      XLALPrintInfo( "%s:     physical %zu = {%.15g,%.15g,%.15g,%.15g}\n", __func__, i, phys[i]->Alpha, phys[i]->Delta, phys[i]->fkdot[0], phys[i]->fkdot[1] );
    }
    for ( size_t i = 0; i < 2; ++i ) {
      XLALPrintInfo( "%s:     lattice %zu = ", __func__, i );
      for ( size_t j = 0; j < latt[i]->size; ++j ) {
        XLALPrintInfo( "%c%.15g", j == 0 ? '{' : ',', gsl_vector_get( latt[i], j ) );
      }
      XLALPrintInfo( "}\n" );
    }
    XLALPrintInfo( "%s:     lattice diff = ", __func__ );
    for ( size_t j = 0; j < u->size; ++j ) {
      XLALPrintInfo( "%c%.15g", j == 0 ? '{' : ',', gsl_vector_get( u, j ) );
    }
    XLALPrintInfo( "}\n" );
    XLALPrintInfo( "%s:     metric dot = ", __func__ );
    for ( size_t j = 0; j < u->size; ++j ) {
      XLALPrintInfo( "%c%.15g", j == 0 ? '{' : ',', gsl_vector_get( u, j ) * gsl_vector_get( v, j ) );
    }
    XLALPrintInfo( "}\n" );
  }

  return XLAL_SUCCESS;

}

///
/// Compare vectors of results from two toplists
///
int toplist_compare_results(
  BOOLEAN *equal,
  const VectorComparison *result_tol,
  const REAL4Vector *res_1,
  const REAL4Vector *res_2
  )
{
  VectorComparison XLAL_INIT_DECL( result_diff );
  int errnum = 0;
  XLAL_TRY( XLALCompareREAL4Vectors( &result_diff, res_1, res_2, result_tol ), errnum );
  if ( errnum == XLAL_ETOL ) {
    *equal = 0;
  } else if ( errnum != 0 ) {
    XLAL_ERROR( XLAL_EFUNC );
  }
  return XLAL_SUCCESS;
}

///
/// Compare two toplists and return whether they are equal
///
int toplist_compare(
  BOOLEAN *equal,
  const WeaveSetupData *setup,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const LALStringVector *detectors,
  const size_t nsegments,
  const LALHeap *toplist_1,
  const LALHeap *toplist_2
  )
{

  // Check input
  XLAL_CHECK( equal != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup != NULL, XLAL_EFAULT );
  XLAL_CHECK( param_tol_mism > 0, XLAL_EINVAL );
  XLAL_CHECK( result_tol != NULL, XLAL_EFAULT );
  XLAL_CHECK( toplist_1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( toplist_2 != NULL, XLAL_EFAULT );

  // Compare lengths of toplists
  const size_t n = XLALHeapSize( toplist_1 );
  {
    const size_t n_2 = XLALHeapSize( toplist_2 );
    if ( n != n_2 ) {
      *equal = 0;
      XLALPrintInfo( "%s: unequal toplist sizes: %zu != %zu\n", __func__, n, n_2 );
      return XLAL_SUCCESS;
    }
  }

  // Get lists of toplist items
  const WeaveOutputToplistItem **items_1 = ( const WeaveOutputToplistItem ** ) XLALHeapElements( toplist_1 );
  XLAL_CHECK( items_1 != NULL, XLAL_EFUNC );
  XLAL_CHECK( nsegments == 0 || items_1[0]->per_seg != NULL, XLAL_EINVAL );
  const WeaveOutputToplistItem **items_2 = ( const WeaveOutputToplistItem ** ) XLALHeapElements( toplist_2 );
  XLAL_CHECK( items_2 != NULL, XLAL_EFUNC );
  XLAL_CHECK( nsegments == 0 || items_2[0]->per_seg != NULL, XLAL_EINVAL );

  // Sort toplist items by physical coordinates of semicoherent template
  // - Template coordinates are less likely to suffer from numerical differences
  //   than result values, and therefore provide more stable sort values to ensure
  //   that equivalent items in both templates match up with each other.
  // - Ideally one would compare toplist items with possess the minimum mismatch
  //   in template parameters with respect to each other, but that would require
  //   of order 'n^2' mismatch calculations, which may be too expensive
  qsort( items_1, n, sizeof( *items_1 ), toplist_item_sort_by_semi_phys );
  qsort( items_2, n, sizeof( *items_2 ), toplist_item_sort_by_semi_phys );

  // Allocate vectors for storing results for comparison with toplist_compare_results()
  REAL4Vector *res_1 = XLALCreateREAL4Vector( n );
  XLAL_CHECK( res_1 != NULL, XLAL_EFUNC );
  REAL4Vector *res_2 = XLALCreateREAL4Vector( n );
  XLAL_CHECK( res_2 != NULL, XLAL_EFUNC );

  // Compare toplist semicoherent and coherent template parameters
  for ( size_t i = 0; i < n; ++i ) {
    char loc_str[256];

    // Compare semicoherent template parameters
    snprintf( loc_str, sizeof(loc_str), "toplist item %zu", i );
    XLAL_CHECK( toplist_compare_templates( equal, loc_str, "semicoherent", param_tol_mism, setup->phys_to_latt, setup->metrics->semi_rssky_metric, setup->metrics->semi_rssky_transf, &items_1[i]->semi_phys, &items_2[i]->semi_phys ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Compare coherent template parameters
    for ( size_t j = 0; j < nsegments; ++j ) {
      snprintf( loc_str, sizeof(loc_str), "toplist item %zu, segment %zu", i, j );
      XLAL_CHECK( toplist_compare_templates( equal, loc_str, "coherent", param_tol_mism, setup->phys_to_latt, setup->metrics->coh_rssky_metric[j], setup->metrics->coh_rssky_transf[j], &items_1[i]->per_seg[j].coh_phys, &items_2[i]->per_seg[j].coh_phys ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  }
  if ( !*equal ) {
    goto CLEANUP;
  }

  // Compare mean multi-detector F-statistics
  XLALPrintInfo( "%s: comparing mean multi-detector F-statistics ...\n", __func__ );
  for ( size_t i = 0; i < n; ++i ) {
    res_1->data[i] = items_1[i]->mean_twoF;
    res_2->data[i] = items_2[i]->mean_twoF;
  }
  XLAL_CHECK( toplist_compare_results( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( !*equal ) {
    goto CLEANUP;
  }

  // Compare mean per-detector F-statistic
  if ( detectors != NULL ) {
    for ( size_t k = 0; k < detectors->length; ++k ) {
      XLALPrintInfo( "%s: comparing mean per-detector F-statistics for detector '%s'...\n", __func__, detectors->data[k] );
      for ( size_t i = 0; i < n; ++i ) {
        res_1->data[i] = items_1[i]->mean_twoF_per_det[k];
        res_2->data[i] = items_2[i]->mean_twoF_per_det[k];
      }
      XLAL_CHECK( toplist_compare_results( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    if ( !*equal ) {
      goto CLEANUP;
    }
  }

  // Compare per-segment coherent multi-detector F-statistics
  for ( size_t j = 0; j < nsegments; ++j ) {
    XLALPrintInfo( "%s: comparing coherent multi-detector F-statistics for segment %zu...\n", __func__, j );
    for ( size_t i = 0; i < n; ++i ) {
      res_1->data[i] = items_1[i]->per_seg[j].twoF;
      res_2->data[i] = items_2[i]->per_seg[j].twoF;
    }
    XLAL_CHECK( toplist_compare_results( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  if ( !*equal ) {
    goto CLEANUP;
  }

  // Compare per-segment per-detector F-statistics
  if ( detectors != NULL ) {
    for ( size_t j = 0; j < nsegments; ++j ) {
      for ( size_t k = 0; k < detectors->length; ++k ) {
        if ( isfinite( items_1[0]->per_seg[j].twoF_per_det[k] ) || isfinite( items_2[0]->per_seg[j].twoF_per_det[k] ) ) {
          XLALPrintInfo( "%s: comparing per-segment per-detector F-statistics for segment %zu, detector '%s'...\n", __func__, j, detectors->data[k] );
          for ( size_t i = 0; i < n; ++i ) {
            res_1->data[i] = items_1[i]->per_seg[j].twoF_per_det[k];
            res_2->data[i] = items_2[i]->per_seg[j].twoF_per_det[k];
          }
          XLAL_CHECK( toplist_compare_results( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        } else {
          XLALPrintInfo( "%s: no per-segment per-detector F-statistics for segment %zu, detector '%s'; skipping comparison\n", __func__, j, detectors->data[k] );
        }
      }
    }
  }
  if ( !*equal ) {
    goto CLEANUP;
  }

CLEANUP:

  // Cleanup
  XLALFree( items_1 );
  XLALFree( items_2 );
  XLALDestroyREAL4Vector( res_1 );
  XLALDestroyREAL4Vector( res_2 );

  return XLAL_SUCCESS;

}

///
/// Create a toplist item
///
WeaveOutputToplistItem *toplist_item_create(
  const LIGOTimeGPS *ref_time,
  const UINT4 per_nsegments
  )
{

  // Allocate memory
  WeaveOutputToplistItem *item = XLALCalloc( 1, sizeof( *item ) );
  XLAL_CHECK_NULL( item != NULL, XLAL_ENOMEM );
  if ( per_nsegments > 0 ) {
    item->per_seg = XLALCalloc( per_nsegments, sizeof( *item->per_seg ) );
    XLAL_CHECK_NULL( item->per_seg != NULL, XLAL_ENOMEM );
  }

  // Set reference time of physical coordinates
  item->semi_phys.refTime = *ref_time;
  for ( size_t j = 0; j < per_nsegments; ++j ) {
    item->per_seg[j].coh_phys.refTime = *ref_time;
  }

  return item;

}

///
/// Destroy a toplist item
///
void toplist_item_destroy(
  void *x
  )
{
  if ( x != NULL ) {
    WeaveOutputToplistItem *ix = ( WeaveOutputToplistItem * ) x;
    XLALFree( ix->per_seg );
    XLALFree( ix );
  }
}

///
/// Sort toplist items by physical coordinates of semicoherent template
///
int toplist_item_sort_by_semi_phys(
  const void *x,
  const void *y
  )
{
  const WeaveOutputToplistItem *ix = *( const WeaveOutputToplistItem *const * ) x;
  const WeaveOutputToplistItem *iy = *( const WeaveOutputToplistItem *const * ) y;
  WEAVE_COMPARE_BY( ix->semi_phys.Alpha, iy->semi_phys.Alpha );   // Compare in ascending order
  WEAVE_COMPARE_BY( ix->semi_phys.Delta, iy->semi_phys.Delta );   // Compare in ascending order
  for ( size_t s = 0; s < XLAL_NUM_ELEM( ix->semi_phys.fkdot ); ++s ) {
    WEAVE_COMPARE_BY( ix->semi_phys.fkdot[s], iy->semi_phys.fkdot[s] );   // Compare in ascending order
  }
  return 0;
}

///
/// Compare toplist items by mean multi-detector F-statistic
///
int toplist_item_compare_by_mean_twoF(
  const void *x,
  const void *y
  )
{
  const WeaveOutputToplistItem *ix = ( const WeaveOutputToplistItem * ) x;
  const WeaveOutputToplistItem *iy = ( const WeaveOutputToplistItem * ) y;
  WEAVE_COMPARE_BY( iy->mean_twoF, ix->mean_twoF );   // Compare in descending order
  return 0;
}

///
/// Fill a toplist item, creating a new one if needed
///
int toplist_item_add(
  BOOLEAN *full_init,
  WeaveOutputResults *out,
  LALHeap *toplist,
  const WeaveSemiResults *semi_res,
  const size_t freq_idx
  )
{

  // Check input
  XLAL_CHECK( full_init != NULL, XLAL_EFAULT );
  XLAL_CHECK( out != NULL, XLAL_EFAULT );
  XLAL_CHECK( toplist != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );

  // Create a new toplist item if needed
  if ( out->saved_item == NULL ) {
    out->saved_item = toplist_item_create( &out->ref_time, out->per_nsegments );
    XLAL_CHECK( out->saved_item != NULL, XLAL_ENOMEM );

    // Toplist item must be fully initialised
    *full_init = 1;

  }

  // Fill toplist item, creating a new one if needed
  XLAL_CHECK( XLALWeaveFillOutputToplistItem( &out->saved_item, full_init, semi_res, freq_idx ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Add item to toplist
  const WeaveOutputToplistItem *prev_item = out->saved_item;
  XLAL_CHECK( XLALHeapAdd( toplist, ( void ** ) &out->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );

  // If toplist added item and returned a different, no-unused item,
  // that item will have to be fully initialised at the next call
  *full_init = ( out->saved_item != prev_item );

  return XLAL_SUCCESS;

}

///
/// Create output results
///
WeaveOutputResults *XLALWeaveOutputResultsCreate(
  const LIGOTimeGPS *ref_time,
  const int toplist_limit,
  const size_t nspins,
  const LALStringVector *per_detectors,
  const UINT4 per_nsegments
  )
{

  // Check input
  XLAL_CHECK_NULL( nspins > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( toplist_limit >= 0, XLAL_EINVAL );

  // Allocate memory
  WeaveOutputResults *out = XLALCalloc( 1, sizeof( *out ) );
  XLAL_CHECK_NULL( out != NULL, XLAL_ENOMEM );

  // Set fields
  out->ref_time = *ref_time;
  out->nspins = nspins;
  out->per_nsegments = per_nsegments;
  out->semi_total = 0;

  // Copy list of detectors
  if ( per_detectors != NULL ) {
    out->per_detectors = XLALCopyStringVector( per_detectors );
    XLAL_CHECK_NULL( out->per_detectors != NULL, XLAL_EFUNC );
  }

  // Create a toplist ranked by mean multi-detector F-statistic
  out->mean_twoF_toplist = toplist_create( toplist_limit, toplist_item_compare_by_mean_twoF );
  XLAL_CHECK_NULL( out->mean_twoF_toplist != NULL, XLAL_EFUNC );

  return out;

}

///
/// Free output results
///
void XLALWeaveOutputResultsDestroy(
  WeaveOutputResults *out
  )
{
  if ( out != NULL ) {
    toplist_item_destroy( out->saved_item );
    XLALHeapDestroy( out->mean_twoF_toplist );
    XLALDestroyStringVector( out->per_detectors );
    XLALFree( out );
  }
}

///
/// Add semicoherent results to output
///
int XLALWeaveOutputResultsAdd(
  WeaveOutputResults *out,
  const WeaveSemiResults *semi_res,
  const UINT4 semi_nfreqs
  )
{

  // Check input
  XLAL_CHECK( out != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );

  // Must initialise all toplist item fields the first time
  BOOLEAN full_init = 1;

  // Iterate over the frequency bins of the semicoherent results
  for ( size_t i = 0; i < semi_nfreqs; ++i ) {

    // Add item to toplist ranked by mean multi-detector F-statistic
    XLAL_CHECK( toplist_item_add( &full_init, out, out->mean_twoF_toplist, semi_res, i ) == XLAL_SUCCESS, XLAL_EFUNC );

  }

  // Increment total number of semicoherent results
  out->semi_total += semi_nfreqs;

  return XLAL_SUCCESS;

}

///
/// Write output results to a FITS file
///
int XLALWeaveOutputResultsWrite(
  FITSFile *file,
  const WeaveOutputResults *out
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( out != NULL, XLAL_EFAULT );

  // Write reference time
  XLAL_CHECK( XLALFITSHeaderWriteGPSTime( file, "date-obs", &out->ref_time, "reference time" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write number of spindowns
  XLAL_CHECK( XLALFITSHeaderWriteINT4( file, "nspins", out->nspins, "number of spindowns" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write if outputting per-detector quantities
  XLAL_CHECK( XLALFITSHeaderWriteBOOLEAN( file, "perdet", out->per_detectors != NULL, "output per detector?" ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( out->per_detectors != NULL ) {
    XLAL_CHECK( XLALFITSHeaderWriteStringVector( file, "detect", out->per_detectors, "setup detectors" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write if outputting per-segment quantities
  XLAL_CHECK( XLALFITSHeaderWriteBOOLEAN( file, "perseg", out->per_nsegments > 0, "output per segment?" ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( out->per_nsegments > 0 ) {
    XLAL_CHECK( XLALFITSHeaderWriteINT4( file, "nsegment", out->per_nsegments, "number of segments" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write total number of semicoherent results added to output
  XLAL_CHECK( XLALFITSHeaderWriteINT8( file, "semitot", out->semi_total, "total semicoherent templates searched" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write toplist ranked by mean multi-detector F-statistic
  XLAL_CHECK( toplist_fits_table_write( file, "mean_twoF_toplist", "toplist ranked by mean multi-detector F-statistic", out, out->mean_twoF_toplist ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Read output results from a FITS file and either create, or append to existing, output results
///
int XLALWeaveOutputResultsReadAppend(
  FITSFile *file,
  WeaveOutputResults **out
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( out != NULL, XLAL_EFAULT );

  // Decide whether to create, or append to existing, output results
  const BOOLEAN create = ( *out == NULL );

  // Allocate memory if required
  if ( create ) {
    *out = XLALCalloc( 1, sizeof( **out ) );
    XLAL_CHECK( *out != NULL, XLAL_ENOMEM );
  }

  // Read and either set, or check, reference time
  {
    LIGOTimeGPS ref_time;
    XLAL_CHECK( XLALFITSHeaderReadGPSTime( file, "date-obs", &ref_time ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( create ) {
      ( *out )->ref_time = ref_time;
    } else {
      XLAL_CHECK( XLALGPSCmp( &ref_time, &( *out )->ref_time ) == 0, XLAL_EIO, "Inconsistent reference time: %" LAL_GPS_FORMAT " != %" LAL_GPS_FORMAT, LAL_GPS_PRINT( ref_time ), LAL_GPS_PRINT( ( *out )->ref_time ) );
    }
  }

  // Read and either set, or check, number of spindowns
  {
    INT4 nspins = 0;
    XLAL_CHECK( XLALFITSHeaderReadINT4( file, "nspins", &nspins ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( create ) {
      ( *out )->nspins = nspins;
    } else {
      XLAL_CHECK( (size_t) nspins == ( *out )->nspins, XLAL_EIO, "Inconsistent number of spindowns: %i != %zu", nspins, ( *out )->nspins );
    }
  }

  // Read and either set, or check, if outputting per-detector quantities, and list of detectors
  {
    BOOLEAN perdet = 0;
    XLAL_CHECK( XLALFITSHeaderReadBOOLEAN( file, "perdet", &perdet ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( perdet ) {
      LALStringVector *per_detectors = NULL;
      XLAL_CHECK( XLALFITSHeaderReadStringVector( file, "detect", &per_detectors ) == XLAL_SUCCESS, XLAL_EFUNC );
      if ( create ) {
        ( *out )->per_detectors = per_detectors;
      } else {
        XLAL_CHECK( ( *out )->per_detectors != NULL, XLAL_EIO, "Inconsistent output per detector?" );
        XLAL_CHECK( per_detectors->length == ( *out )->per_detectors->length, XLAL_EIO, "Inconsistent number of detectors: %u != %u", per_detectors->length, ( *out )->per_detectors->length );
        for ( size_t i = 0; i < per_detectors->length; ++i ) {
          XLAL_CHECK( strcmp( per_detectors->data[i], ( *out )->per_detectors->data[i] ) == 0, XLAL_EIO, "Inconsistent detectors: %s != %s", per_detectors->data[i], ( *out )->per_detectors->data[i] );
        }
        XLALDestroyStringVector( per_detectors );
      }
    } else {
      if ( create ) {
        ( *out )->per_detectors = NULL;
      } else {
        XLAL_CHECK( ( *out )->per_detectors == NULL, XLAL_EIO, "Inconsistent output per detector?" );
      }
    }
  }

  // Read and either set, or check, if outputting per-segment quantities, and number of per-segment items
  {
    BOOLEAN perseg = 0;
    XLAL_CHECK( XLALFITSHeaderReadBOOLEAN( file, "perseg", &perseg ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( perseg ) {
      INT4 per_nsegments = 0;
      XLAL_CHECK( XLALFITSHeaderReadINT4( file, "nsegment", &per_nsegments ) == XLAL_SUCCESS, XLAL_EFUNC );
      if ( create ) {
        ( *out )->per_nsegments = per_nsegments;
      } else {
        XLAL_CHECK( (size_t) per_nsegments == ( *out )->per_nsegments, XLAL_EIO, "Inconsistent number of segments: %i != %u", per_nsegments, ( *out )->per_nsegments );
      }
    } else {
      if ( create ) {
        ( *out )->per_nsegments = 0;
      } else {
        XLAL_CHECK( ( *out )->per_nsegments == 0, XLAL_EIO, "Inconsistent output per segment?" );
      }
    }
  }

  // Read and increment total number of semicoherent results added to output
  {
    INT8 semi_total = 0;
    XLAL_CHECK( XLALFITSHeaderReadINT8( file, "semitot", &semi_total ) == XLAL_SUCCESS, XLAL_EFUNC );
    ( *out )->semi_total += semi_total;
  }

  // Read and append to toplist ranked by mean multi-detector F-statistic
  XLAL_CHECK( toplist_fits_table_read( file, "mean_twoF_toplist", *out, &( *out )->mean_twoF_toplist, toplist_item_compare_by_mean_twoF ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Compare two output results and return whether they are equal
///
int XLALWeaveOutputResultsCompare(
  BOOLEAN *equal,
  const WeaveSetupData *setup,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const WeaveOutputResults *out_1,
  const WeaveOutputResults *out_2
  )
{

  // Check input
  XLAL_CHECK( equal != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup != NULL, XLAL_EFAULT );
  XLAL_CHECK( param_tol_mism > 0, XLAL_EINVAL );
  XLAL_CHECK( result_tol != NULL, XLAL_EFAULT );
  XLAL_CHECK( out_1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( out_2 != NULL, XLAL_EFAULT );

  // Output results are assumed equal until we find otherwise
  *equal = 1;

  // Compare reference times
  if ( XLALGPSCmp( &out_1->ref_time, &out_2->ref_time ) != 0 ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal reference times: %" LAL_GPS_FORMAT " != %" LAL_GPS_FORMAT "\n", __func__, LAL_GPS_PRINT( out_1->ref_time ), LAL_GPS_PRINT( out_2->ref_time ) );
    return XLAL_SUCCESS;
  }

  // Compare number of spindowns
  if ( out_1->nspins != out_2->nspins ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal number of spindowns: %zu != %zu\n", __func__, out_1->nspins, out_2->nspins );
    return XLAL_SUCCESS;
  }

  // Compare if outputting per-detector quantities, and list of detectors
  {
    const BOOLEAN perdet_1 = ( out_1->per_detectors != NULL ), perdet_2 = ( out_2->per_detectors != NULL );
    if ( perdet_1 != perdet_2 ) {
      *equal = 0;
      XLALPrintInfo( "%s: unequal output per detector?: %i != %i\n", __func__, perdet_1, perdet_2 );
      return XLAL_SUCCESS;
    }
    if ( perdet_1 ) {
      if ( out_1->per_detectors->length != out_2->per_detectors->length ) {
        *equal = 0;
        XLALPrintInfo( "%s: unequal number of detectors: %u != %u\n", __func__, out_1->per_detectors->length, out_2->per_detectors->length );
        return XLAL_SUCCESS;
      }
      for ( size_t i = 0; i < out_1->per_detectors->length; ++i ) {
        if ( strcmp( out_1->per_detectors->data[i], out_2->per_detectors->data[i] ) != 0 ) {
          *equal = 0;
          XLALPrintInfo( "%s: unequal detectors: %s != %s\n", __func__, out_1->per_detectors->data[i], out_2->per_detectors->data[i] );
          return XLAL_SUCCESS;
        }
      }
    }
  }
  const LALStringVector *detectors = out_1->per_detectors;

  // Compare number of per-segment items
  if ( out_1->per_nsegments != out_2->per_nsegments ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal number of segments: %u != %u\n", __func__, out_1->per_nsegments, out_2->per_nsegments );
    return XLAL_SUCCESS;
  }
  const size_t nsegments = out_1->per_nsegments;

  // Compare total number of semicoherent results
  if ( out_1->semi_total != out_2->semi_total ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal total number of semicoherent results: %" LAL_INT8_FORMAT " != %" LAL_INT8_FORMAT "\n", __func__, out_1->semi_total, out_2->semi_total );
    return XLAL_SUCCESS;
  }

  // Compare toplists ranked by mean multi-detector F-statistic
  XLALPrintInfo( "%s: comparing toplists ranked by mean multi-detector F-statistic ...\n", __func__ );
  XLAL_CHECK( toplist_compare( equal, setup, param_tol_mism, result_tol, detectors, nsegments, out_1->mean_twoF_toplist, out_2->mean_twoF_toplist ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( !*equal ) {
    return XLAL_SUCCESS;
  }

  return XLAL_SUCCESS;

}
