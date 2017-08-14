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

#include "ResultsToplist.h"

#include <lal/VectorMath.h>
#include <lal/UserInputPrint.h>

///
/// Internal definition of toplist of output results
///
struct tagWeaveResultsToplist {
  /// Struct holding all parameters for which statistics to output and compute, when, and how
  /// NOTE: this is only a reference-pointer to WeaveStatisticsParams, while WeaveSemiResults is the *owner*
  const WeaveStatisticsParams *statistics_params;
  /// Number of spindown parameters to output
  size_t nspins;
  /// Name of ranking statistic
  const char *stat_name;
  /// Description of ranking statistic
  const char *stat_desc;
  /// Function which returns pointer to array of statistics by which toplist items are ranked
  WeaveResultsToplistRankingStats rank_stats_fcn;
  /// Function which returns the value of the statistic by which toplist items are ranked
  WeaveResultsToplistItemGetRankStat item_get_rank_stat_fcn;
  /// Function which sets the value of the statistic by which toplist items are ranked
  WeaveResultsToplistItemSetRankStat item_set_rank_stat_fcn;
  /// Vector of indexes of toplist results which should be considered for addition
  UINT4Vector *maybe_add_freq_idxs;
  /// Heap which ranks toplist items
  LALHeap *heap;
  /// Save a no-longer-used toplist item for re-use
  WeaveResultsToplistItem *saved_item;
};

///
/// \name Internal functions
///
/// @{

static WeaveResultsToplistItem *toplist_item_create( const WeaveResultsToplist *toplist );
static int compare_templates( BOOLEAN *equal, const char *loc_str, const char *tmpl_str, const REAL8 param_tol_mism, const WeavePhysicalToLattice phys_to_latt, const gsl_matrix *metric, const void *transf_data, const PulsarDopplerParams *phys_1, const PulsarDopplerParams *phys_2 );
static int compare_vectors( BOOLEAN *equal, const VectorComparison *result_tol, const REAL4Vector *res_1, const REAL4Vector *res_2 );
static int toplist_fits_table_init( FITSFile *file, const WeaveResultsToplist *toplist );
static int toplist_fits_table_write_visitor( void *param, const void *x );
static int toplist_item_sort_by_semi_phys( const void *x, const void *y );
static void toplist_item_destroy( WeaveResultsToplistItem *item );
static int toplist_item_compare( void *param, const void *x, const void *y );
static int toplist_fill_completionloop_stats( void *param, void *x );

/// @}

///
/// Create a toplist item
///
WeaveResultsToplistItem *toplist_item_create(
  const WeaveResultsToplist *toplist
  )
{

  // Check input
  XLAL_CHECK_NULL( toplist != NULL, XLAL_EFAULT );

  // Allocate memory for item
  WeaveResultsToplistItem *item = XLALCalloc( 1, sizeof( *item ) );
  XLAL_CHECK_NULL( item != NULL, XLAL_ENOMEM );

  const WeaveStatisticsParams *params = toplist-> statistics_params;
  WeaveStatisticType store_per_segment_stats = (params -> statistics_to_output & (WEAVE_STATISTIC_COH2F | WEAVE_STATISTIC_COH2F_DET));
  WeaveStatisticType store_coh2F             = (params -> mainloop_statistics_to_keep & WEAVE_STATISTIC_COH2F);
  WeaveStatisticType store_coh2F_det         = (params -> mainloop_statistics_to_keep & WEAVE_STATISTIC_COH2F_DET);

  // Allocate memory for per-segment output results
  if ( store_per_segment_stats ) {
    item->coh_alpha = XLALCalloc( params -> nsegments, sizeof( *item->coh_alpha ) );
    XLAL_CHECK_NULL( item->coh_alpha != NULL, XLAL_ENOMEM );
    item->coh_delta = XLALCalloc( params -> nsegments, sizeof( *item->coh_delta ) );
    XLAL_CHECK_NULL( item->coh_delta != NULL, XLAL_ENOMEM );
    for ( size_t k = 0; k <= toplist->nspins; ++k ) {
      item->coh_fkdot[k] = XLALCalloc( params -> nsegments, sizeof( *item->coh_fkdot[k] ) );
      XLAL_CHECK_NULL( item->coh_fkdot[k] != NULL, XLAL_ENOMEM );
    }
  }

  // Allocate memory for per-segment coh2F and coh2F_det statistics if requested
  if ( store_coh2F ) {
    item->coh2F = XLALCalloc( params -> nsegments, sizeof( *item->coh2F ) );
    XLAL_CHECK_NULL( item->coh2F != NULL, XLAL_ENOMEM );
  }
  if ( store_coh2F_det ) {
    for ( size_t i = 0; i < params -> detectors->length; ++i ) {
      item->coh2F_det[i] = XLALCalloc( params -> nsegments, sizeof( *item->coh2F_det[i] ) );
      XLAL_CHECK_NULL( item->coh2F_det[i] != NULL, XLAL_ENOMEM );
    }
  }

  return item;

} // toplist_item_create()

///
/// Destroy a toplist item
///
void toplist_item_destroy(
  WeaveResultsToplistItem *item
  )
{
  if ( item != NULL ) {
    XLALFree( item->coh_alpha );
    XLALFree( item->coh_delta );
    for ( size_t k = 0; k < PULSAR_MAX_SPINS; ++k ) {
      XLALFree( item->coh_fkdot[k] );
    }
    XLALFree( item->coh2F );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALFree( item->coh2F_det[i] );
    }
    XLALFree( item );
  }
} // toplist_item_destroy()

///
/// Compare toplist items
///
int toplist_item_compare(
  void *param,
  const void *x,
  const void *y
  )
{
  WeaveResultsToplistItemGetRankStat item_get_rank_stat_fcn = ( WeaveResultsToplistItemGetRankStat ) param;
  const WeaveResultsToplistItem *ix = ( const WeaveResultsToplistItem * ) x;
  const WeaveResultsToplistItem *iy = ( const WeaveResultsToplistItem * ) y;
  WEAVE_COMPARE_BY( item_get_rank_stat_fcn( iy ), item_get_rank_stat_fcn( ix ) );   // Compare in descending order
  return 0;
} // toplist_item_compare()

///
/// Initialise a FITS table for writing/reading a toplist
///
int toplist_fits_table_init(
  FITSFile *file,
  const WeaveResultsToplist *toplist
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );

  char col_name[64];

  // Begin FITS table description
  XLAL_FITS_TABLE_COLUMN_BEGIN( WeaveResultsToplistItem );

  // Add columns for semicoherent template parameters
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_alpha, "alpha [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_delta, "delta [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_fkdot[0], "freq [Hz]" ) == XLAL_SUCCESS, XLAL_EFUNC );
  for ( size_t k = 1; k <= toplist->nspins; ++k ) {
    snprintf( col_name, sizeof( col_name ), "f%zudot [Hz/s^%zu]", k, k );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_fkdot[k], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // which statistics have been requested for output?
  const WeaveStatisticsParams *params = toplist -> statistics_params;
  WeaveStatisticType statistics_to_output = params -> statistics_to_output;

  // Add column for mean multi-detector F-statistic
  if ( statistics_to_output & WEAVE_STATISTIC_MEAN2F ) {
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, mean2F ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  // Add columns for mean per-detector F-statistic
  if ( statistics_to_output & WEAVE_STATISTIC_MEAN2F_DET ) {
    for ( size_t i = 0; i < params -> detectors -> length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "mean2F_%s", params -> detectors -> data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL4, mean2F_det[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // We tie the output of per-segment coordinates to the output of any per-segment statistics
  WeaveStatisticType per_segment_stats = (statistics_to_output & (WEAVE_STATISTIC_COH2F | WEAVE_STATISTIC_COH2F_DET));
  // Add columns for coherent template parameters
  if ( per_segment_stats ) {
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_PTR_ARRAY_NAMED( file, REAL8, params -> nsegments, coh_alpha, "alpha_seg [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_PTR_ARRAY_NAMED( file, REAL8, params -> nsegments, coh_delta, "delta_seg [rad]" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_PTR_ARRAY_NAMED( file, REAL8, params -> nsegments, coh_fkdot[0], "freq_seg [Hz]" ) == XLAL_SUCCESS, XLAL_EFUNC );
    for ( size_t k = 1; k <= toplist->nspins; ++k ) {
      snprintf( col_name, sizeof( col_name ), "f%zudot_seg [Hz/s^%zu]", k, k );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_PTR_ARRAY_NAMED( file, REAL8, params -> nsegments, coh_fkdot[k], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }
  // Add column for per-segment coherent multi-detector 2F statistics
  if ( statistics_to_output & WEAVE_STATISTIC_COH2F ) {
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_PTR_ARRAY_NAMED( file, REAL4, params -> nsegments, coh2F, "coh2F" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  // Add columns for per-segment coherent per-detector 2F statistics
  if ( statistics_to_output & WEAVE_STATISTIC_COH2F_DET ) {
    for ( size_t i = 0; i < params -> detectors->length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "coh2F_%s", params -> detectors->data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_PTR_ARRAY_NAMED( file, REAL4, params -> nsegments, coh2F_det[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Add column for multi-detector max2F statistic
  if ( statistics_to_output & WEAVE_STATISTIC_MAX2F ) {
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, max2F ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  // Add columns for per-detector max2F_det statistic
  if ( statistics_to_output & WEAVE_STATISTIC_MAX2F_DET ) {
    for ( size_t i = 0; i < params -> detectors->length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "max2F_%s", params -> detectors->data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL4, max2F_det[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Add column for multi-detector sum2F statistic
  if ( statistics_to_output & WEAVE_STATISTIC_SUM2F ) {
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, sum2F ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  // Add columns for per-detector sum2F statistic
  if ( statistics_to_output & WEAVE_STATISTIC_SUM2F_DET ) {
    for ( size_t i = 0; i < params -> detectors->length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "sum2F_%s", params -> detectors->data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL4, sum2F_det[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  // Add column for BSGL statistic
  if ( statistics_to_output & WEAVE_STATISTIC_BSGL ) {
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, log10BSGL ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Add column for BSGLtL statistic
  if ( statistics_to_output & WEAVE_STATISTIC_BSGLtL ) {
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, log10BSGLtL ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Add column for BtSGLtL statistic
  if ( statistics_to_output & WEAVE_STATISTIC_BtSGLtL ) {
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, log10BtSGLtL ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Add column for multi-detector number-count statistic
  if ( statistics_to_output & WEAVE_STATISTIC_NCOUNT ) {
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, ncount ) == XLAL_SUCCESS, XLAL_EFUNC );
  }
  // Add column for per-detector number-count statistics
  if ( statistics_to_output & WEAVE_STATISTIC_NCOUNT_DET ) {
    for ( size_t i = 0; i < params -> detectors->length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "ncount_%s", params -> detectors->data[i] );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL4, ncount_det[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  return XLAL_SUCCESS;

} // toplist_fits_table_init()

///
/// Function to update given toplist item with missing 'completion loop' statistics
///
int toplist_fill_completionloop_stats(
  void *param,
  void *x
  )
{
  XLAL_CHECK ( param != NULL, XLAL_EFAULT );
  XLAL_CHECK ( x != NULL, XLAL_EFAULT );

  WeaveResultsToplistItem *item = ( WeaveResultsToplistItem * ) x;
  WeaveStatisticsParams *stats_params = (WeaveStatisticsParams *)param;
  WeaveStatisticType completionloop_stats = stats_params -> completionloop_statistics;
  UINT4 ndetectors = stats_params -> detectors -> length;
  UINT4 nsegments  = stats_params -> nsegments;

  if ( completionloop_stats & WEAVE_STATISTIC_MAX2F ) {
    item -> max2F = 0;
    for ( size_t l = 0; l < nsegments; ++l ) {
      item -> max2F = fmaxf( item -> max2F, item->coh2F[l] );
    }
  }

  if ( completionloop_stats & WEAVE_STATISTIC_MAX2F_DET ) {
    for ( size_t X = 0; X < ndetectors; ++X ) {
      item -> max2F_det[X] = 0;
      for ( size_t l = 0; l < nsegments; ++l ) {
        REAL4 item_Xl = item->coh2F_det[X][l];
        item -> max2F_det[X] = isnan(item_Xl) ? item -> max2F_det[X] : fmaxf ( item -> max2F_det[X], item_Xl );
      }
    }
  }

  if ( completionloop_stats & WEAVE_STATISTIC_SUM2F ) {
    item -> sum2F = 0;
    for ( size_t l = 0; l < nsegments; ++l ) {
      item -> sum2F += item->coh2F[l];
    }
  }

  if ( completionloop_stats & WEAVE_STATISTIC_SUM2F_DET ) {
    for ( size_t X = 0; X < ndetectors; ++X ) {
      item -> sum2F_det[X] = 0;
      for ( size_t l = 0; l < nsegments; ++l ) {
        REAL4 item_Xl = item->coh2F_det[X][l];
        item -> sum2F_det[X] += isnan(item_Xl) ? 0 : item_Xl;
      }
    }
  }

  if ( completionloop_stats & WEAVE_STATISTIC_MEAN2F ) {
    item -> mean2F = item -> sum2F / stats_params->nsum2F;
  }

  if ( completionloop_stats & WEAVE_STATISTIC_MEAN2F_DET ) {
    for ( size_t X = 0; X < ndetectors; ++X ) {
      item -> mean2F_det[X] = item -> sum2F_det[X] / stats_params->nsum2F_det[X];
    }
  }

  if ( completionloop_stats & WEAVE_STATISTIC_BSGL ) {
    item -> log10BSGL = XLALComputeBSGL ( item -> sum2F, item -> sum2F_det, stats_params -> BSGL_setup );
  }

  if ( completionloop_stats & WEAVE_STATISTIC_BSGLtL ) {
    item -> log10BSGLtL = XLALComputeBSGLtL ( item -> sum2F, item -> sum2F_det, item -> max2F_det, stats_params -> BSGL_setup );
  }

  if ( completionloop_stats & WEAVE_STATISTIC_BtSGLtL ) {
    item -> log10BtSGLtL = XLALComputeBtSGLtL ( item -> max2F, item -> sum2F_det, item -> max2F_det, stats_params -> BSGL_setup );
  }

  if ( completionloop_stats & WEAVE_STATISTIC_NCOUNT ) {
    item -> ncount = 0;
    for ( size_t l = 0; l < nsegments; ++l ) {
      item -> ncount += ( item->coh2F[l] > stats_params->nc_2Fth ) ? 1 : 0;
    }
  }

  if ( completionloop_stats & WEAVE_STATISTIC_NCOUNT_DET ) {
    for ( size_t X = 0; X < ndetectors; ++X ) {
      item -> ncount_det[X] = 0;
      for ( size_t l = 0; l < nsegments; ++l ) {
        REAL4 item_Xl = item->coh2F_det[X][l];
        if ( isnan(item_Xl) ) { continue; }
        item -> ncount_det[X] += ( item_Xl > stats_params->nc_2Fth ) ? 1 : 0;
      }
    }
  }

  return XLAL_SUCCESS;

} // toplist_fill_completionloop_stats()

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
} // toplist_fits_table_write_visitor()

///
/// Sort toplist items by physical coordinates of semicoherent template.
///
/// For stable comparisons, the order of parameter comparisons should be the same
/// as the order in which parameters are generated by the search lattice tiling.
///
int toplist_item_sort_by_semi_phys(
  const void *x,
  const void *y
  )
{
  const WeaveResultsToplistItem *ix = *( const WeaveResultsToplistItem *const * ) x;
  const WeaveResultsToplistItem *iy = *( const WeaveResultsToplistItem *const * ) y;
  WEAVE_COMPARE_BY( ix->semi_alpha, iy->semi_alpha );   // Compare in ascending order
  WEAVE_COMPARE_BY( ix->semi_delta, iy->semi_delta );   // Compare in ascending order
  for ( size_t s = 1; s < XLAL_NUM_ELEM( ix->semi_fkdot ); ++s ) {
    WEAVE_COMPARE_BY( ix->semi_fkdot[s], iy->semi_fkdot[s] );   // Compare in ascending order
  }
  WEAVE_COMPARE_BY( ix->semi_fkdot[0], iy->semi_fkdot[0] );   // Compare in ascending order
  return 0;
} // toplist_item_sort_by_semi_phys()

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
  XLAL_CHECK( param_tol_mism >= 0, XLAL_EINVAL );
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

} // compare_templates()

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
} // compare_vectors()

///
/// Create results toplist
///
WeaveResultsToplist *XLALWeaveResultsToplistCreate(
  const size_t nspins,
  const WeaveStatisticsParams *statistics_params,
  const char *stat_name,
  const char *stat_desc,
  const UINT4 toplist_limit,
  WeaveResultsToplistRankingStats toplist_rank_stats_fcn,
  WeaveResultsToplistItemGetRankStat toplist_item_get_rank_stat_fcn,
  WeaveResultsToplistItemSetRankStat toplist_item_set_rank_stat_fcn
  )
{

  // Check input
  XLAL_CHECK_NULL( stat_name != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( stat_desc != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( statistics_params != NULL, XLAL_EFAULT );

  // Allocate memory
  WeaveResultsToplist *toplist = XLALCalloc( 1, sizeof( *toplist ) );
  XLAL_CHECK_NULL( toplist != NULL, XLAL_ENOMEM );

  // Set fields
  toplist->nspins = nspins;
  toplist->stat_name = stat_name;
  toplist->stat_desc = stat_desc;
  toplist->rank_stats_fcn = toplist_rank_stats_fcn;
  toplist->item_get_rank_stat_fcn = toplist_item_get_rank_stat_fcn;
  toplist->item_set_rank_stat_fcn = toplist_item_set_rank_stat_fcn;
  toplist->statistics_params = statistics_params;

  // Create heap which ranks toplist items
  toplist->heap = XLALHeapCreate2( ( LALHeapDtorFcn ) toplist_item_destroy, toplist_limit, +1, toplist_item_compare, toplist_item_get_rank_stat_fcn );
  XLAL_CHECK_NULL( toplist->heap != NULL, XLAL_EFUNC );

  return toplist;

} // XLALWeaveResultsToplistCreate()

///
/// Free results toplist
///
void XLALWeaveResultsToplistDestroy(
  WeaveResultsToplist *toplist
  )
{
  if ( toplist != NULL ) {
    XLALDestroyUINT4Vector( toplist->maybe_add_freq_idxs );
    XLALHeapDestroy( toplist->heap );
    toplist_item_destroy( toplist->saved_item );
    XLALFree( toplist );
  }
} // XLALWeaveResultsToplistDestroy()

///
/// Add semicoherent results to toplist
///
int XLALWeaveResultsToplistAdd(
  WeaveResultsToplist *toplist,
  const WeaveSemiResults *semi_res,
  const UINT4 semi_nfreqs
  )
{
  // Check input
  XLAL_CHECK( toplist != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );

  // Reallocate vector of indexes of toplist results which should be considered for addition
  if ( toplist->maybe_add_freq_idxs == NULL || toplist->maybe_add_freq_idxs->length < semi_nfreqs ) {
    toplist->maybe_add_freq_idxs = XLALResizeUINT4Vector( toplist->maybe_add_freq_idxs, semi_nfreqs );
    XLAL_CHECK( toplist->maybe_add_freq_idxs != NULL, XLAL_EFUNC );
  }

  // Get pointer to array of ranking statistics
  const REAL4 *toplist_rank_stats = toplist->rank_stats_fcn( semi_res );

  // Get ranking statistic of heap root (or -infinity if heap is not yet full)
  const int heap_full = XLALHeapIsFull( toplist->heap );
  XLAL_CHECK( heap_full >= 0, XLAL_EFUNC );
  const REAL4 heap_root_rank_stat = heap_full ? toplist->item_get_rank_stat_fcn( XLALHeapRoot( toplist->heap ) ) : GSL_NEGINF;

  // Find the indexes of the semicoherent results whose ranking statistic equals or exceeds
  // that of the heap root; only select these results for possible insertion into the toplist
  UINT4 n_maybe_add = 0;
  XLAL_CHECK( XLALVectorFindScalarLessEqualREAL4( &n_maybe_add, toplist->maybe_add_freq_idxs->data, heap_root_rank_stat, toplist_rank_stats, semi_nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

  // whether we output per-segment template coordinates is currently tied to output of any per-segment statistics
  const WeaveStatisticsParams *params = toplist->statistics_params;
  WeaveStatisticType per_seg_coords = params -> statistics_to_output & (WEAVE_STATISTIC_COH2F | WEAVE_STATISTIC_COH2F_DET);

  // Iterate over semicoherent results which have been selected for possible toplist insertion
  for ( UINT4 idx = 0; idx < n_maybe_add; ++idx ) {
    const UINT4 freq_idx = toplist->maybe_add_freq_idxs->data[idx];

    // Create a new toplist item if needed
    if ( toplist->saved_item == NULL ) {
      toplist->saved_item = toplist_item_create( toplist );
      XLAL_CHECK( toplist->saved_item != NULL, XLAL_ENOMEM );
    }
    WeaveResultsToplistItem *item = toplist->saved_item;

    // Set ranking statistic of toplist item
    toplist->item_set_rank_stat_fcn( item, toplist_rank_stats[freq_idx] );

    // Possibly add toplist item to heap
    XLAL_CHECK( XLALHeapAdd( toplist->heap, ( void ** ) &toplist->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Skip remainder of loop if toplist item was not added to heap
    if ( item == toplist->saved_item ) {
      continue;
    }

    // Set all semicoherent template parameters
    item->semi_alpha = semi_res->semi_phys.Alpha;
    item->semi_delta = semi_res->semi_phys.Delta;
    for ( size_t k = 0; k <= toplist->nspins; ++k ) {
      item->semi_fkdot[k] = semi_res->semi_phys.fkdot[k];
    }

    // Set all coherent template parameters if outputting per-segment statistics
    if ( per_seg_coords ) {
      for ( size_t j = 0; j < semi_res->nsegments; ++j ) {
        item->coh_alpha[j] = semi_res->coh_phys[j].Alpha;
        item->coh_delta[j] = semi_res->coh_phys[j].Delta;
        for ( size_t k = 0; k <= toplist->nspins; ++k ) {
          item->coh_fkdot[k][j] = semi_res->coh_phys[j].fkdot[k];
        }
      }
    }

    // Update semicoherent and coherent template frequency
    item->semi_fkdot[0] = semi_res->semi_phys.fkdot[0] + freq_idx * semi_res->dfreq;
    if ( per_seg_coords ) {
      for ( size_t j = 0; j < semi_res->nsegments; ++j ) {
        item->coh_fkdot[0][j] = semi_res->coh_phys[j].fkdot[0] + freq_idx * semi_res->dfreq;
      }
    }

    // Skip remainder of loop if simulating search
    if ( semi_res->simulation_level & WEAVE_SIMULATE ) {
      continue;
    }

    //
    // copy all 'mainloop_statistics_to_keep' statistic values, as they will be needed either 1) for output or 2) computing remaining completion-loop statistics
    //
    WeaveStatisticType stats_to_keep  = params -> mainloop_statistics_to_keep;

    if ( stats_to_keep & WEAVE_STATISTIC_COH2F ) {
      for ( size_t j = 0; j < semi_res->nsegments; ++j ) {
        item->coh2F[j] = (semi_res->coh2F[j] != NULL) ? semi_res->coh2F[j][freq_idx] : NAN;
      }
    }
    if ( stats_to_keep & WEAVE_STATISTIC_COH2F_DET ) {
      for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
        for ( size_t j = 0; j < semi_res->nsegments; ++j ) {
          if ( semi_res->coh2F_det[i][j] != NULL ) {
            item->coh2F_det[i][j] = semi_res->coh2F_det[i][j][freq_idx];
          } else {
            // There is not per-detector F-statistic for this segment, usually because this segment contains
            // no data from this detector. In this case we output a clearly invalid F-statistic value.
            item->coh2F_det[i][j] = NAN;
          }
        }
      }
    }

    if ( stats_to_keep & WEAVE_STATISTIC_MAX2F ) {
      item->max2F = semi_res->max2F->data[freq_idx];
    }
    if ( stats_to_keep & WEAVE_STATISTIC_MAX2F_DET ) {
      for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
        item->max2F_det[i]  = semi_res->max2F_det[i]->data[freq_idx];
      }
    }

    if ( stats_to_keep & WEAVE_STATISTIC_SUM2F ) {
      item->sum2F = semi_res->sum2F->data[freq_idx];
    }
    if ( stats_to_keep & WEAVE_STATISTIC_SUM2F_DET ) {
      for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
        item->sum2F_det[i]  = semi_res->sum2F_det[i]->data[freq_idx];
      }
    }

    if ( stats_to_keep & WEAVE_STATISTIC_MEAN2F ) {
      item->mean2F = semi_res->mean2F->data[freq_idx];
    }
    if ( stats_to_keep & WEAVE_STATISTIC_MEAN2F_DET ) {
      for ( size_t i = 0; i < semi_res->ndetectors; ++i ) {
        item->mean2F_det[i]  = semi_res->mean2F_det[i]->data[freq_idx];
      }
    }

    if ( stats_to_keep & WEAVE_STATISTIC_BSGL ) {
      item->log10BSGL = semi_res->log10BSGL->data[freq_idx];
    }

    if ( stats_to_keep & WEAVE_STATISTIC_BSGLtL ) {
      item->log10BSGLtL = semi_res->log10BSGLtL->data[freq_idx];
    }

    if ( stats_to_keep & WEAVE_STATISTIC_BtSGLtL ) {
      item->log10BtSGLtL = semi_res->log10BtSGLtL->data[freq_idx];
    }

  } // for i < n_maybe_add

  return XLAL_SUCCESS;

} // XLALWeaveResultsToplistAdd()

///
/// Compute all missing 'extra' (non-toplist-ranking) statistics for all toplist entries
///
int XLALWeaveResultsToplistCompletionLoop(
  WeaveResultsToplist *toplist
  )
{
  // Check input
  XLAL_CHECK( toplist != NULL, XLAL_EFAULT );

  // compute all completion-loop statistics on toplist items
  WeaveStatisticsParams params = *toplist->statistics_params;	// local copy
  XLAL_CHECK( XLALHeapModify( toplist->heap, toplist_fill_completionloop_stats, &params ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // XLALWeaveResultsToplistCompletionLoop()

///
/// Write results toplist to a FITS file
///
int XLALWeaveResultsToplistWrite(
  FITSFile *file,
  const WeaveResultsToplist *toplist
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( toplist != NULL, XLAL_EFAULT );

  // Format name and description of statistic
  char name[256];
  snprintf( name, sizeof( name ), "%s_toplist", toplist->stat_name );
  char desc[256];
  snprintf( desc, sizeof( desc ), "toplist ranked by %s", toplist->stat_desc );

  // Open FITS table for writing and initialise
  XLAL_CHECK( XLALFITSTableOpenWrite( file, name, desc ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( toplist_fits_table_init( file, toplist ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write all heap items to FITS table
  XLAL_CHECK( XLALHeapVisit( toplist->heap, toplist_fits_table_write_visitor, file ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // XLALWeaveResultsToplistWrite

///
/// Read results from a FITS file and append to existing results toplist
///
int XLALWeaveResultsToplistReadAppend(
  FITSFile *file,
  WeaveResultsToplist *toplist
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( toplist != NULL, XLAL_EFAULT );

  // Format name of statistic
  char name[256];
  snprintf( name, sizeof( name ), "%s_toplist", toplist->stat_name );

  // Open FITS table for reading and initialise
  UINT8 nrows = 0;
  XLAL_CHECK( XLALFITSTableOpenRead( file, name, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( toplist_fits_table_init( file, toplist ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Read all items from FITS table
  while ( nrows > 0 ) {

    // Create a new toplist item if needed
    if ( toplist->saved_item == NULL ) {
      toplist->saved_item = toplist_item_create( toplist );
      XLAL_CHECK( toplist->saved_item != NULL, XLAL_ENOMEM );
    }

    // Read item from FITS table
    XLAL_CHECK( XLALFITSTableReadRow( file, toplist->saved_item, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Add item to heap
    XLAL_CHECK( XLALHeapAdd( toplist->heap, ( void ** ) &toplist->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );

  }

  return XLAL_SUCCESS;

} // XLALWeaveResultsToplistReadAppend()

///
/// Compare two results toplists and return whether they are equal
///
int XLALWeaveResultsToplistCompare(
  BOOLEAN *equal,
  const WeaveSetupData *setup,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const WeaveResultsToplist *toplist_1,
  const WeaveResultsToplist *toplist_2
  )
{

  // Check input
  XLAL_CHECK( equal != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup != NULL, XLAL_EFAULT );
  XLAL_CHECK( param_tol_mism >= 0, XLAL_EINVAL );
  XLAL_CHECK( result_tol != NULL, XLAL_EFAULT );
  XLAL_CHECK( toplist_1 != NULL, XLAL_EFAULT );
  XLAL_CHECK( toplist_2 != NULL, XLAL_EFAULT );
  XLAL_CHECK( strcmp( toplist_1->stat_name, toplist_2->stat_name ) == 0, XLAL_EINVAL );
  XLAL_CHECK( strcmp( toplist_1->stat_desc, toplist_2->stat_desc ) == 0, XLAL_EINVAL );

  const WeaveResultsToplist *toplist = toplist_1;
  const WeaveStatisticsParams *params = toplist -> statistics_params;
  WeaveStatisticType stats_to_output = params->statistics_to_output;

  // Results toplists are assumed equal until we find otherwise
  *equal = 1;

  // Compare toplists
  XLALPrintInfo( "%s: comparing toplists ranked by %s ...\n", __func__, toplist->stat_desc );
  {

    // Compare lengths of heaps
    const size_t n = XLALHeapSize( toplist_1->heap );
    {
      const size_t n_2 = XLALHeapSize( toplist_2->heap );
      if ( n != n_2 ) {
        *equal = 0;
        XLALPrintInfo( "%s: unequal size %s toplists: %zu != %zu\n", __func__, toplist->stat_desc, n, n_2 );
        return XLAL_SUCCESS;
      }
    }

    // Get lists of toplist items
    const WeaveResultsToplistItem **items_1 = ( const WeaveResultsToplistItem ** ) XLALHeapElements( toplist_1->heap );
    XLAL_CHECK( items_1 != NULL, XLAL_EFUNC );
    const WeaveResultsToplistItem **items_2 = ( const WeaveResultsToplistItem ** ) XLALHeapElements( toplist_2->heap );
    XLAL_CHECK( items_2 != NULL, XLAL_EFUNC );

    // Sort toplist items by physical coordinates of semicoherent template
    // - Template coordinates are less likely to suffer from numerical differences
    //   than result values, and therefore provide more stable sort values to ensure
    //   that equivalent items in both templates match up with each other.
    // - Ideally one would compare toplist items with possess the minimum mismatch
    //   in template parameters with respect to each other, but that would require
    //   of order 'n^2' mismatch calculations, which may be too expensive
    qsort( items_1, n, sizeof( *items_1 ), toplist_item_sort_by_semi_phys );
    qsort( items_2, n, sizeof( *items_2 ), toplist_item_sort_by_semi_phys );

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
        {
          snprintf( loc_str, sizeof(loc_str), "toplist item %zu", i );
          PulsarDopplerParams XLAL_INIT_DECL( semi_phys_1 );
          PulsarDopplerParams XLAL_INIT_DECL( semi_phys_2 );
          semi_phys_1.Alpha = items_1[i]->semi_alpha;
          semi_phys_2.Alpha = items_2[i]->semi_alpha;
          semi_phys_1.Delta = items_1[i]->semi_delta;
          semi_phys_2.Delta = items_2[i]->semi_delta;
          for ( size_t k = 0; k <= toplist->nspins; ++k ) {
            semi_phys_1.fkdot[k] = items_1[i]->semi_fkdot[k];
            semi_phys_2.fkdot[k] = items_2[i]->semi_fkdot[k];
          };
          XLAL_CHECK( compare_templates( equal, loc_str, "semicoherent", param_tol_mism, setup->phys_to_latt, setup->metrics->semi_rssky_metric, setup->metrics->semi_rssky_transf, &semi_phys_1, &semi_phys_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        }

        // Compare coherent template parameters
        if ( stats_to_output & (WEAVE_STATISTIC_COH2F|WEAVE_STATISTIC_COH2F_DET) ) {
          for ( size_t j = 0; j < params -> nsegments; ++j ) {
            snprintf( loc_str, sizeof(loc_str), "toplist item %zu, segment %zu", i, j );
            PulsarDopplerParams XLAL_INIT_DECL( coh_phys_1 );
            PulsarDopplerParams XLAL_INIT_DECL( coh_phys_2 );
            coh_phys_1.Alpha = items_1[i]->coh_alpha[j];
            coh_phys_2.Alpha = items_2[i]->coh_alpha[j];
            coh_phys_1.Delta = items_1[i]->coh_delta[j];
            coh_phys_2.Delta = items_2[i]->coh_delta[j];
            for ( size_t k = 0; k <= toplist->nspins; ++k ) {
              coh_phys_1.fkdot[k] = items_1[i]->coh_fkdot[k][j];
              coh_phys_2.fkdot[k] = items_2[i]->coh_fkdot[k][j];
            };
            XLAL_CHECK( compare_templates( equal, loc_str, "coherent", param_tol_mism, setup->phys_to_latt, setup->metrics->coh_rssky_metric[j], setup->metrics->coh_rssky_transf[j], &coh_phys_1, &coh_phys_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
          }
        } // if output per-segment coordinates

      }
      if ( !*equal ) {
        break;
      }

      // Compare mean multi-detector F-statistics
      if ( stats_to_output & WEAVE_STATISTIC_MEAN2F ) {
        XLALPrintInfo( "%s: comparing mean multi-detector F-statistics ...\n", __func__ );
        for ( size_t i = 0; i < n; ++i ) {
          res_1->data[i] = items_1[i]->mean2F;
          res_2->data[i] = items_2[i]->mean2F;
        }
        XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        if ( !*equal ) {
          break;
        }
      }

      // Compare mean per-detector F-statistic
      if ( stats_to_output & WEAVE_STATISTIC_MEAN2F_DET ) {
        for ( size_t k = 0; k < params -> detectors->length; ++k ) {
          XLALPrintInfo( "%s: comparing mean per-detector F-statistics for detector '%s'...\n", __func__, params -> detectors->data[k] );
          for ( size_t i = 0; i < n; ++i ) {
            res_1->data[i] = items_1[i]->mean2F_det[k];
            res_2->data[i] = items_2[i]->mean2F_det[k];
          }
          XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
        if ( !*equal ) {
          break;
        }
      }

      // Compare per-segment coherent multi-detector F-statistics
      if ( stats_to_output & WEAVE_STATISTIC_COH2F ) {
        for ( size_t j = 0; j < params -> nsegments; ++j ) {
          XLALPrintInfo( "%s: comparing coherent multi-detector F-statistics for segment %zu...\n", __func__, j );
          for ( size_t i = 0; i < n; ++i ) {
            res_1->data[i] = items_1[i]->coh2F[j];
            res_2->data[i] = items_2[i]->coh2F[j];
          }
          XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
        if ( !*equal ) {
          break;
        }
      }

      // Compare per-segment per-detector F-statistics
      if ( stats_to_output & WEAVE_STATISTIC_COH2F_DET ) {
        for ( size_t j = 0; j < params -> nsegments; ++j ) {
          for ( size_t k = 0; k < params -> detectors->length; ++k ) {
            if ( isfinite( items_1[0]->coh2F_det[k][j] ) || isfinite( items_2[0]->coh2F_det[k][j] ) ) {
              XLALPrintInfo( "%s: comparing per-segment per-detector F-statistics for segment %zu, detector '%s'...\n", __func__, j, params -> detectors->data[k] );
              for ( size_t i = 0; i < n; ++i ) {
                res_1->data[i] = items_1[i]->coh2F_det[k][j];
                res_2->data[i] = items_2[i]->coh2F_det[k][j];
              }
              XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
            } else {
              XLALPrintInfo( "%s: no per-segment per-detector F-statistics for segment %zu, detector '%s'; skipping comparison\n", __func__, j, params -> detectors->data[k] );
            }
          }
        }
      }
      if ( !*equal ) {
        break;
      }

      // Compare segment-max multi-detector F-statistics
      if ( stats_to_output & WEAVE_STATISTIC_MAX2F ) {
        XLALPrintInfo( "%s: comparing max multi-detector F-statistics ...\n", __func__ );
        for ( size_t i = 0; i < n; ++i ) {
          res_1->data[i] = items_1[i]->max2F;
          res_2->data[i] = items_2[i]->max2F;
        }
        XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        if ( !*equal ) {
          break;
        }
      }

      // Compare segment-max per-detector F-statistic
      if ( stats_to_output & WEAVE_STATISTIC_MAX2F_DET ) {
        for ( size_t k = 0; k < params -> detectors->length; ++k ) {
          XLALPrintInfo( "%s: comparing max per-detector F-statistics for detector '%s'...\n", __func__, params -> detectors->data[k] );
          for ( size_t i = 0; i < n; ++i ) {
            res_1->data[i] = items_1[i]->max2F_det[k];
            res_2->data[i] = items_2[i]->max2F_det[k];
          }
          XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
        if ( !*equal ) {
          break;
        }
      }

      // Compare summed multi-detector F-statistics
      if ( stats_to_output & WEAVE_STATISTIC_SUM2F ) {
        XLALPrintInfo( "%s: comparing sum multi-detector F-statistics ...\n", __func__ );
        for ( size_t i = 0; i < n; ++i ) {
          res_1->data[i] = items_1[i]->sum2F;
          res_2->data[i] = items_2[i]->sum2F;
        }
        XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        if ( !*equal ) {
          break;
        }
      }

      // Compare sum per-detector F-statistic
      if ( stats_to_output & WEAVE_STATISTIC_SUM2F_DET ) {
        for ( size_t k = 0; k < params -> detectors->length; ++k ) {
          XLALPrintInfo( "%s: comparing sum per-detector F-statistics for detector '%s'...\n", __func__, params -> detectors->data[k] );
          for ( size_t i = 0; i < n; ++i ) {
            res_1->data[i] = items_1[i]->sum2F_det[k];
            res_2->data[i] = items_2[i]->sum2F_det[k];
          }
          XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        }
        if ( !*equal ) {
          break;
        }
      }

      // Compare line-robust BSGL statistic
      if ( stats_to_output & WEAVE_STATISTIC_BSGL ) {
        XLALPrintInfo( "%s: comparing line-robust B_S/GL statistic ...\n", __func__ );
        for ( size_t i = 0; i < n; ++i ) {
          res_1->data[i] = items_1[i]->log10BSGL;
          res_2->data[i] = items_2[i]->log10BSGL;
        }
        XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        if ( !*equal ) {
          break;
        }
      }

      // Compare transient line-robust BSGLtL statistic
      if ( stats_to_output & WEAVE_STATISTIC_BSGLtL ) {
        XLALPrintInfo( "%s: comparing transient line-robust B_S/GLtL statistic ...\n", __func__ );
        for ( size_t i = 0; i < n; ++i ) {
          res_1->data[i] = items_1[i]->log10BSGLtL;
          res_2->data[i] = items_2[i]->log10BSGLtL;
        }
        XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        if ( !*equal ) {
          break;
        }
      }

      // Compare transient signal line-robust BtSGLtL statistic
      if ( stats_to_output & WEAVE_STATISTIC_BtSGLtL ) {
        XLALPrintInfo( "%s: comparing transient signal line-robust B_tS/GLtL statistic ...\n", __func__ );
        for ( size_t i = 0; i < n; ++i ) {
          res_1->data[i] = items_1[i]->log10BtSGLtL;
          res_2->data[i] = items_2[i]->log10BtSGLtL;
        }
        XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        if ( !*equal ) {
          break;
        }
      }

      // Compare 'Hough' multi-detector line statistics
      if ( stats_to_output & WEAVE_STATISTIC_NCOUNT ) {
        XLALPrintInfo( "%s: comparing 'Hough' multi-detector number count statistic ...\n", __func__ );
        for ( size_t i = 0; i < n; ++i ) {
          res_1->data[i] = items_1[i]->ncount;
          res_2->data[i] = items_2[i]->ncount;
        }
        XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
        if ( !*equal ) {
          break;
        }
      }
      // Compare 'Hough' per-detector line statistics
      if ( stats_to_output & WEAVE_STATISTIC_NCOUNT_DET ) {
        for ( size_t k = 0; k < params -> detectors->length; ++k ) {
          XLALPrintInfo( "%s: comparing 'Hough' per-detector number-count statistic for detector '%s'...\n", __func__, params -> detectors->data[k] );
          for ( size_t i = 0; i < n; ++i ) {
            res_1->data[i] = items_1[i]->ncount_det[k];
            res_2->data[i] = items_2[i]->ncount_det[k];
          }
          XLAL_CHECK( compare_vectors( equal, result_tol, res_1, res_2 ) == XLAL_SUCCESS, XLAL_EFUNC );
          if ( !*equal ) {
            break;
          }
        }
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

} // XLALWeaveResultsToplistCompare()

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
