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

#include "ResultsBasket.h"

///
/// Internal definition of basket of output results
///
struct tagWeaveResultsBasket {
  /// Various varameters required to output results
  const WeaveOutputParams *par;
  /// Name of ranking statistic
  const char *stat_name;
  /// Description of ranking statistic
  const char *stat_desc;
  /// Toplist which ranks output results by a particular statistic
  LALHeap *toplist;
  /// Save a no-longer-used output result item for re-use
  WeaveOutputResultItem *saved_item;
};

///
/// \name Internal routines
///
/// @{

static int compare_templates( BOOLEAN *equal, const char *loc_str, const char *tmpl_str, const REAL8 param_tol_mism, const WeavePhysicalToLattice phys_to_latt, const gsl_matrix *metric, const void *transf_data, const PulsarDopplerParams *phys_1, const PulsarDopplerParams *phys_2 );
static int compare_vectors( BOOLEAN *equal, const VectorComparison *result_tol, const REAL4Vector *res_1, const REAL4Vector *res_2 );
static int result_item_sort_by_semi_phys( const void *x, const void *y );
static int toplist_fits_table_init( FITSFile *file, const WeaveOutputParams *par );
static int toplist_fits_table_write_visitor( void *param, const void *x );

/// @}

///
/// Initialise a FITS table for writing/reading a toplist
///
int toplist_fits_table_init(
  FITSFile *file,
  const WeaveOutputParams *par
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );

  char col_name[32];

  // Begin FITS table description
  XLAL_FITS_TABLE_COLUMN_BEGIN( WeaveOutputResultItem );

  // Add columns for semicoherent template parameters
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_phys.Alpha, "alpha [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_phys.Delta, "delta [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_phys.fkdot[0], "freq [Hz]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  for ( size_t k = 1; k <= par->nspins; ++k ) {
    snprintf( col_name, sizeof( col_name ), "f%zudot [Hz/s^%zu]", k, k );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_phys.fkdot[k], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Add columns for mean multi- and per-detector F-statistic
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, mean_twoF ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( par->per_detectors != NULL ) {
    for ( size_t i = 0; i < par->per_detectors->length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "mean_twoF_%s", par->per_detectors->data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL4, mean_twoF_per_det[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Begin FITS table description for per-segment items (optional)
  if ( par->per_nsegments > 0 ) {
    XLAL_FITS_TABLE_COLUMN_PTR_BEGIN( per_seg, WeaveOutputPerSegResultItem, par->per_nsegments );
    for ( size_t s = 0; s < par->per_nsegments; ++s ) {

      // Add columns for coherent template parameters
      snprintf( col_name, sizeof( col_name ), "seg%zu_alpha [rad]", s + 1 );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_phys.Alpha, col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      snprintf( col_name, sizeof( col_name ), "seg%zu_delta [rad]", s + 1 );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_phys.Delta, col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      snprintf( col_name, sizeof( col_name ), "seg%zu_freq [Hz]", s + 1 );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_phys.fkdot[0], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      for ( size_t k = 1; k <= par->nspins; ++k ) {
        snprintf( col_name, sizeof( col_name ), "seg%zu_f%zudot [Hz/s^%zu]", s + 1, k, k );
        XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_phys.fkdot[k], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      }

      // Add columns for coherent multi- and per-detector F-statistic
      snprintf( col_name, sizeof( col_name ), "seg%zu_twoF", s + 1 );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL4, twoF, col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      if ( par->per_detectors != NULL ) {
        for ( size_t i = 0; i < par->per_detectors->length; ++i ) {
          snprintf( col_name, sizeof( col_name ), "seg%zu_twoF_%s", s + 1, par->per_detectors->data[i] );
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
/// Sort output result items by physical coordinates of semicoherent template
///
int result_item_sort_by_semi_phys(
  const void *x,
  const void *y
  )
{
  const WeaveOutputResultItem *ix = *( const WeaveOutputResultItem *const * ) x;
  const WeaveOutputResultItem *iy = *( const WeaveOutputResultItem *const * ) y;
  WEAVE_COMPARE_BY( ix->semi_phys.Alpha, iy->semi_phys.Alpha );   // Compare in ascending order
  WEAVE_COMPARE_BY( ix->semi_phys.Delta, iy->semi_phys.Delta );   // Compare in ascending order
  for ( size_t s = 0; s < XLAL_NUM_ELEM( ix->semi_phys.fkdot ); ++s ) {
    WEAVE_COMPARE_BY( ix->semi_phys.fkdot[s], iy->semi_phys.fkdot[s] );   // Compare in ascending order
  }
  return 0;
}

///
/// Compute two template parameters
///
int compare_templates(
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
/// Compare two vectors of results
///
int compare_vectors(
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
/// Create results basket
///
WeaveResultsBasket *XLALWeaveResultsBasketCreate(
  const WeaveOutputParams *par,
  const char *stat_name,
  const char *stat_desc,
  const int toplist_limit,
  LALHeapCmpFcn toplist_result_item_compare_fcn
  )
{

  // Check input
  XLAL_CHECK_NULL( par != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( stat_name != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( stat_desc != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( toplist_limit >= 0, XLAL_EINVAL );

  // Allocate memory
  WeaveResultsBasket *basket = XLALCalloc( 1, sizeof( *basket ) );
  XLAL_CHECK_NULL( basket != NULL, XLAL_ENOMEM );

  // Set fields
  basket->par = par;
  basket->stat_name = stat_name;
  basket->stat_desc = stat_desc;

  // Create toplist which ranks output results using the given comparison function
  basket->toplist = XLALHeapCreate( ( LALHeapDtorFcn ) XLALWeaveOutputResultItemDestroy, toplist_limit, +1, toplist_result_item_compare_fcn );
  XLAL_CHECK_NULL( basket->toplist != NULL, XLAL_EFUNC );

  return basket;

}

///
/// Free results basket
///
void XLALWeaveResultsBasketDestroy(
  WeaveResultsBasket *basket
  )
{
  if ( basket != NULL ) {
    XLALHeapDestroy( basket->toplist );
    XLALWeaveOutputResultItemDestroy( basket->saved_item );
    XLALFree( basket );
  }
}

///
/// Add semicoherent results to basket
///
int XLALWeaveResultsBasketAdd(
  WeaveResultsBasket *basket,
  const WeaveSemiResults *semi_res,
  const UINT4 semi_nfreqs
  )
{

  // Check input
  XLAL_CHECK( basket != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );

  // Must initialise all output result item fields the first time
  BOOLEAN full_init = 1;

  // Iterate over the frequency bins of the semicoherent results
  for ( size_t i = 0; i < semi_nfreqs; ++i ) {

    // Create a new output result item if needed
    if ( basket->saved_item == NULL ) {
      basket->saved_item = XLALWeaveOutputResultItemCreate( basket->par );
      XLAL_CHECK( basket->saved_item != NULL, XLAL_ENOMEM );

      // Output result item must be fully initialised
      full_init = 1;

    }

    // Fill output result item, creating a new one if needed
    XLAL_CHECK( XLALWeaveFillOutputResultItem( &basket->saved_item, &full_init, semi_res, i ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Add item to toplist
    const WeaveOutputResultItem *prev_item = basket->saved_item;
    XLAL_CHECK( XLALHeapAdd( basket->toplist, ( void ** ) &basket->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );

    // If toplist added item and returned a different, now-unused item,
    // that item will have to be fully initialised at the next call
    full_init = ( basket->saved_item != prev_item );

  }

  return XLAL_SUCCESS;

}

///
/// Write results basket to a FITS file
///
int XLALWeaveResultsBasketWrite(
  FITSFile *file,
  const WeaveResultsBasket *basket
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( basket != NULL, XLAL_EFAULT );

  // Write toplist
  {

    // Format name and description of statistic
    char name[256];
    snprintf( name, sizeof( name ), "%s_toplist", basket->stat_name );
    char desc[256];
    snprintf( desc, sizeof( desc ), "toplist ranked by %s", basket->stat_desc );

    // Open FITS table for writing and initialise
    XLAL_CHECK( XLALFITSTableOpenWrite( file, name, desc ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( toplist_fits_table_init( file, basket->par ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Write all toplist items to FITS table
    XLAL_CHECK( XLALHeapVisit( basket->toplist, toplist_fits_table_write_visitor, file ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Write maximum size of toplist to FITS header
    XLAL_CHECK( XLALFITSHeaderWriteINT8( file, "toplimit", XLALHeapMaxSize( basket->toplist ), "maximum size of toplist" ) == XLAL_SUCCESS, XLAL_EFUNC );

  }

  return XLAL_SUCCESS;

}

///
/// Read results from a FITS file and append to existing results basket
///
int XLALWeaveResultsBasketReadAppend(
  FITSFile *file,
  WeaveResultsBasket *basket
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( basket != NULL, XLAL_EFAULT );

  // Read and append to toplist
  {

    // Format name of statistic
    char name[256];
    snprintf( name, sizeof( name ), "%s_toplist", basket->stat_name );

    // Open FITS table for reading and initialise
    UINT8 nrows = 0;
    XLAL_CHECK( XLALFITSTableOpenRead( file, name, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( toplist_fits_table_init( file, basket->par ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Read maximum size of toplist from FITS header
    INT8 toplist_limit = 0;
    XLAL_CHECK( XLALFITSHeaderReadINT8( file, "toplimit", &toplist_limit ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Maximize size of toplist
    if ( toplist_limit > XLALHeapSize( basket->toplist ) ) {
      XLAL_CHECK( XLALHeapResize( basket->toplist, toplist_limit ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

    // Read all items from FITS table
    while ( nrows > 0 ) {

      // Create a new output result item if needed
      if ( basket->saved_item == NULL ) {
        basket->saved_item = XLALWeaveOutputResultItemCreate( basket->par );
        XLAL_CHECK( basket->saved_item != NULL, XLAL_ENOMEM );
      }

      // Read item from FITS table
      XLAL_CHECK( XLALFITSTableReadRow( file, basket->saved_item, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );

      // Add item to toplist
      XLAL_CHECK( XLALHeapAdd( basket->toplist, ( void ** ) &basket->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );

    }

  }

  return XLAL_SUCCESS;

}

///
/// Compare two results baskets and return whether they are equal
///
int XLALWeaveResultsBasketCompare(
  BOOLEAN *equal,
  const WeaveSetupData *setup,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const WeaveResultsBasket *basket_1,
  const WeaveResultsBasket *basket_2
  )
{

  // Check input
  XLAL_CHECK( equal != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup != NULL, XLAL_EFAULT );
  XLAL_CHECK( param_tol_mism > 0, XLAL_EINVAL );
  XLAL_CHECK( result_tol != NULL, XLAL_EFAULT );
  XLAL_CHECK( basket_1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( basket_2 != NULL, XLAL_EFAULT );
  XLAL_CHECK( strcmp( basket_1->stat_name, basket_2->stat_name ) == 0, XLAL_EINVAL );
  XLAL_CHECK( strcmp( basket_1->stat_desc, basket_2->stat_desc ) == 0, XLAL_EINVAL );

  const WeaveOutputParams *par = basket_1->par;
  const char *stat_desc = basket_1->stat_desc;

  // Results baskets are assumed equal until we find otherwise
  *equal = 1;

  // Compare toplists
  XLALPrintInfo( "%s: comparing toplists ranked by %s ...\n", __func__, stat_desc );
  {

    // Compare lengths of toplists
    const size_t n = XLALHeapSize( basket_1->toplist );
    {
      const size_t n_2 = XLALHeapSize( basket_2->toplist );
      if ( n != n_2 ) {
        *equal = 0;
        XLALPrintInfo( "%s: unequal size %s toplists: %zu != %zu\n", __func__, stat_desc, n, n_2 );
        return XLAL_SUCCESS;
      }
    }

    // Get lists of output result items
    const WeaveOutputResultItem **items_1 = ( const WeaveOutputResultItem ** ) XLALHeapElements( basket_1->toplist );
    XLAL_CHECK( items_1 != NULL, XLAL_EFUNC );
    XLAL_CHECK( par->per_nsegments == 0 || items_1[0]->per_seg != NULL, XLAL_EINVAL );
    const WeaveOutputResultItem **items_2 = ( const WeaveOutputResultItem ** ) XLALHeapElements( basket_2->toplist );
    XLAL_CHECK( items_2 != NULL, XLAL_EFUNC );
    XLAL_CHECK( par->per_nsegments == 0 || items_2[0]->per_seg != NULL, XLAL_EINVAL );

    // Sort output result items by physical coordinates of semicoherent template
    // - Template coordinates are less likely to suffer from numerical differences
    //   than result values, and therefore provide more stable sort values to ensure
    //   that equivalent items in both templates match up with each other.
    // - Ideally one would compare output result items with possess the minimum mismatch
    //   in template parameters with respect to each other, but that would require
    //   of order 'n^2' mismatch calculations, which may be too expensive
    qsort( items_1, n, sizeof( *items_1 ), result_item_sort_by_semi_phys );
    qsort( items_2, n, sizeof( *items_2 ), result_item_sort_by_semi_phys );

    // Allocate vectors for storing results for comparison with compare_vectors()
    REAL4Vector *res_1 = XLALCreateREAL4Vector( n );
    XLAL_CHECK( res_1 != NULL, XLAL_EFUNC );
    REAL4Vector *res_2 = XLALCreateREAL4Vector( n );
    XLAL_CHECK( res_2 != NULL, XLAL_EFUNC );

    do {   // So we can use 'break' to skip comparisons on failure

      // Compare semicoherent and coherent template parameters
      for ( size_t i = 0; i < n; ++i ) {
        char loc_str[256];

        // Compare semicoherent template parameters
        snprintf( loc_str, sizeof(loc_str), "output result item %zu", i );
        XLAL_CHECK( compare_templates( equal, loc_str, "semicoherent", param_tol_mism, setup->phys_to_latt, setup->metrics->semi_rssky_metric, setup->metrics->semi_rssky_transf, &items_1[i]->semi_phys, &items_2[i]->semi_phys ) == XLAL_SUCCESS, XLAL_EFUNC );

        // Compare coherent template parameters
        for ( size_t j = 0; j < par->per_nsegments; ++j ) {
          snprintf( loc_str, sizeof(loc_str), "output result item %zu, segment %zu", i, j );
          XLAL_CHECK( compare_templates( equal, loc_str, "coherent", param_tol_mism, setup->phys_to_latt, setup->metrics->coh_rssky_metric[j], setup->metrics->coh_rssky_transf[j], &items_1[i]->per_seg[j].coh_phys, &items_2[i]->per_seg[j].coh_phys ) == XLAL_SUCCESS, XLAL_EFUNC );
        }

      }
      if ( !*equal ) {
        break;
      }

      // Compare mean multi-detector F-statistics
      XLALPrintInfo( "%s: comparing mean multi-detector F-statistics ...\n", __func__ );
      for ( size_t i = 0; i < n; ++i ) {
        res_1->data[i] = items_1[i]->mean_twoF;
        res_2->data[i] = items_2[i]->mean_twoF;
      }
      XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
      if ( !*equal ) {
        break;
      }

      // Compare mean per-detector F-statistic
      if ( par->per_detectors != NULL ) {
        for ( size_t k = 0; k < par->per_detectors->length; ++k ) {
          XLALPrintInfo( "%s: comparing mean per-detector F-statistics for detector '%s'...\n", __func__, par->per_detectors->data[k] );
          for ( size_t i = 0; i < n; ++i ) {
            res_1->data[i] = items_1[i]->mean_twoF_per_det[k];
            res_2->data[i] = items_2[i]->mean_twoF_per_det[k];
          }
          XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
        if ( !*equal ) {
          break;
        }
      }

      // Compare per-segment coherent multi-detector F-statistics
      for ( size_t j = 0; j < par->per_nsegments; ++j ) {
        XLALPrintInfo( "%s: comparing coherent multi-detector F-statistics for segment %zu...\n", __func__, j );
        for ( size_t i = 0; i < n; ++i ) {
          res_1->data[i] = items_1[i]->per_seg[j].twoF;
          res_2->data[i] = items_2[i]->per_seg[j].twoF;
        }
        XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
      if ( !*equal ) {
        break;
      }

      // Compare per-segment per-detector F-statistics
      if ( par->per_detectors != NULL ) {
        for ( size_t j = 0; j < par->per_nsegments; ++j ) {
          for ( size_t k = 0; k < par->per_detectors->length; ++k ) {
            if ( isfinite( items_1[0]->per_seg[j].twoF_per_det[k] ) || isfinite( items_2[0]->per_seg[j].twoF_per_det[k] ) ) {
              XLALPrintInfo( "%s: comparing per-segment per-detector F-statistics for segment %zu, detector '%s'...\n", __func__, j, par->per_detectors->data[k] );
              for ( size_t i = 0; i < n; ++i ) {
                res_1->data[i] = items_1[i]->per_seg[j].twoF_per_det[k];
                res_2->data[i] = items_2[i]->per_seg[j].twoF_per_det[k];
              }
              XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
            } else {
              XLALPrintInfo( "%s: no per-segment per-detector F-statistics for segment %zu, detector '%s'; skipping comparison\n", __func__, j, par->per_detectors->data[k] );
            }
          }
        }
      }
      if ( !*equal ) {
        break;
      }

    } while (0);

    // Cleanup
    XLALFree( items_1 );
    XLALFree( items_2 );
    XLALDestroyREAL4Vector( res_1 );
    XLALDestroyREAL4Vector( res_2 );

    if ( !*equal ) {
      return XLAL_SUCCESS;
    }

  }

  return XLAL_SUCCESS;

}
