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
#include "ResultsToplist.h"

///
/// Internal definition of output results from a search
///
struct tagWeaveOutputResults {
  /// Various varameters required to output results
  WeaveOutputParams par;
  /// Number of computed semicoherent results
  INT8 semi_ncomp;
  /// Toplist of output results ranked by mean multi-detector F-statistic
  WeaveResultsToplist *mean_twoF_toplist;
};

///
/// \name Internal routines
///
/// @{

static int result_item_compare_by_mean_twoF( const void *x, const void *y );

/// @}

///
/// Compare output result items by mean multi-detector F-statistic
///
int result_item_compare_by_mean_twoF(
  const void *x,
  const void *y
  )
{
  const WeaveOutputResultItem *ix = ( const WeaveOutputResultItem * ) x;
  const WeaveOutputResultItem *iy = ( const WeaveOutputResultItem * ) y;
  WEAVE_COMPARE_BY( iy->mean_twoF, ix->mean_twoF );   // Compare in descending order
  return 0;
}

///
/// Create a output result item
///
WeaveOutputResultItem *XLALWeaveOutputResultItemCreate(
  const WeaveOutputParams *par
  )
{

  // Check input
  XLAL_CHECK_NULL( par != NULL, XLAL_EFAULT );

  // Allocate memory for item
  WeaveOutputResultItem *item = XLALCalloc( 1, sizeof( *item ) );
  XLAL_CHECK_NULL( item != NULL, XLAL_ENOMEM );

  // Allocate memory for per-segment output results
  if ( par->per_nsegments > 0 ) {
    item->coh_alpha = XLALCalloc( par->per_nsegments, sizeof( *item->coh_alpha ) );
    XLAL_CHECK_NULL( item->coh_alpha != NULL, XLAL_ENOMEM );
    item->coh_delta = XLALCalloc( par->per_nsegments, sizeof( *item->coh_delta ) );
    XLAL_CHECK_NULL( item->coh_delta != NULL, XLAL_ENOMEM );
    for ( size_t k = 0; k <= par->nspins; ++k ) {
      item->coh_fkdot[k] = XLALCalloc( par->per_nsegments, sizeof( *item->coh_fkdot[k] ) );
      XLAL_CHECK_NULL( item->coh_fkdot[k] != NULL, XLAL_ENOMEM );
    }
    item->twoF = XLALCalloc( par->per_nsegments, sizeof( *item->twoF ) );
    XLAL_CHECK_NULL( item->twoF != NULL, XLAL_ENOMEM );
  }

  // Allocate memory for per-detector and per-segment output results
  if ( par->per_detectors != NULL && par->per_nsegments > 0 ) {
    for ( size_t i = 0; i < par->per_detectors->length; ++i ) {
      item->twoF_per_det[i] = XLALCalloc( par->per_nsegments, sizeof( *item->twoF_per_det[i] ) );
      XLAL_CHECK_NULL( item->twoF_per_det[i] != NULL, XLAL_ENOMEM );
    }
  }

  return item;

}

///
/// Destroy a output result item
///
void XLALWeaveOutputResultItemDestroy(
  WeaveOutputResultItem *item
  )
{
  if ( item != NULL ) {
    XLALFree( item->coh_alpha );
    XLALFree( item->coh_delta );
    for ( size_t k = 0; k < PULSAR_MAX_SPINS; ++k ) {
      XLALFree( item->coh_fkdot[k] );
    }
    XLALFree( item->twoF );
    for ( size_t i = 0; i < PULSAR_MAX_DETECTORS; ++i ) {
      XLALFree( item->twoF_per_det[i] );
    }
    XLALFree( item );
  }
}

///
/// Create output results
///
WeaveOutputResults *XLALWeaveOutputResultsCreate(
  const LIGOTimeGPS *ref_time,
  const size_t nspins,
  const LALStringVector *per_detectors,
  const UINT4 per_nsegments,
  const int toplist_limit
  )
{

  // Check input
  XLAL_CHECK_NULL( ref_time != NULL, XLAL_EFAULT );

  // Allocate memory
  WeaveOutputResults *out = XLALCalloc( 1, sizeof( *out ) );
  XLAL_CHECK_NULL( out != NULL, XLAL_ENOMEM );

  // Set fields
  out->par.ref_time = *ref_time;
  out->par.nspins = nspins;
  out->par.per_nsegments = per_nsegments;

  // Copy list of detectors
  if ( per_detectors != NULL ) {
    out->par.per_detectors = XLALCopyStringVector( per_detectors );
    XLAL_CHECK_NULL( out->par.per_detectors != NULL, XLAL_EFUNC );
  }

  // Create a toplist which ranks results by mean multi-detector F-statistic
  out->mean_twoF_toplist = XLALWeaveResultsToplistCreate( &out->par, "mean_twoF", "mean multi-detector F-statistic", toplist_limit, result_item_compare_by_mean_twoF );
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
    XLALDestroyStringVector( out->par.per_detectors );
    XLALWeaveResultsToplistDestroy( out->mean_twoF_toplist );
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

  // Increment number of computed semicoherent results
  out->semi_ncomp += semi_nfreqs;

  // Add results to toplist ranked by mean multi-detector F-statistic
  XLAL_CHECK( XLALWeaveResultsToplistAdd( out->mean_twoF_toplist, semi_res, semi_nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

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
  XLAL_CHECK( XLALFITSHeaderWriteGPSTime( file, "date-obs", &out->par.ref_time, "reference time" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write number of spindowns
  XLAL_CHECK( XLALFITSHeaderWriteUINT4( file, "nspins", out->par.nspins, "number of spindowns" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write if outputting per-detector quantities
  XLAL_CHECK( XLALFITSHeaderWriteBOOLEAN( file, "perdet", out->par.per_detectors != NULL, "output per detector?" ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( out->par.per_detectors != NULL ) {
    XLAL_CHECK( XLALFITSHeaderWriteStringVector( file, "detect", out->par.per_detectors, "setup detectors" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write if outputting per-segment quantities
  XLAL_CHECK( XLALFITSHeaderWriteBOOLEAN( file, "perseg", out->par.per_nsegments > 0, "output per segment?" ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( out->par.per_nsegments > 0 ) {
    XLAL_CHECK( XLALFITSHeaderWriteUINT4( file, "nsegment", out->par.per_nsegments, "number of segments" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write number of computed semicoherent results
  XLAL_CHECK( XLALFITSHeaderWriteUINT8( file, "semicomp", out->semi_ncomp, "number of computed semicoherent results" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write toplist ranked by mean multi-detector F-statistic
  XLAL_CHECK( XLALWeaveResultsToplistWrite( file, out->mean_twoF_toplist ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Read results from a FITS file and append to new/existing output results
///
int XLALWeaveOutputResultsReadAppend(
  FITSFile *file,
  WeaveOutputResults **out
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( out != NULL, XLAL_EFAULT );

  // Read reference time
  LIGOTimeGPS ref_time;
  XLAL_CHECK( XLALFITSHeaderReadGPSTime( file, "date-obs", &ref_time ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Read number of spindowns
  UINT4 nspins = 0;
  XLAL_CHECK( XLALFITSHeaderReadUINT4( file, "nspins", &nspins ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Read if outputting per-detector quantities, and list of detectors
  BOOLEAN perdet = 0;
  LALStringVector *per_detectors = NULL;
  XLAL_CHECK( XLALFITSHeaderReadBOOLEAN( file, "perdet", &perdet ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( perdet ) {
    XLAL_CHECK( XLALFITSHeaderReadStringVector( file, "detect", &per_detectors ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Read if outputting per-segment quantities, and number of per-segment items
  BOOLEAN perseg = 0;
  UINT4 per_nsegments = 0;
  XLAL_CHECK( XLALFITSHeaderReadBOOLEAN( file, "perseg", &perseg ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( perseg ) {
    XLAL_CHECK( XLALFITSHeaderReadUINT4( file, "nsegment", &per_nsegments ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  if ( *out == NULL ) {

    // Create new output results
    *out = XLALWeaveOutputResultsCreate( &ref_time, nspins, per_detectors, per_nsegments, 0 );
    XLAL_CHECK( *out != NULL, XLAL_EFUNC );

  } else {

    // Check reference time
    XLAL_CHECK( XLALGPSCmp( &ref_time, &( *out )->par.ref_time ) == 0, XLAL_EIO, "Inconsistent reference time: %" LAL_GPS_FORMAT " != %" LAL_GPS_FORMAT, LAL_GPS_PRINT( ref_time ), LAL_GPS_PRINT( ( *out )->par.ref_time ) );

    // Check number of spindowns
    XLAL_CHECK( (size_t) nspins == ( *out )->par.nspins, XLAL_EIO, "Inconsistent number of spindowns: %i != %zu", nspins, ( *out )->par.nspins );

    // Check if outputting per-detector quantities, and list of detectors
    XLAL_CHECK( ( !perdet ) == ( ( *out )->par.per_detectors == NULL ), XLAL_EIO, "Inconsistent output per detector flag vs list of detectors: %i != %i", perseg, ( *out )->par.per_detectors != NULL );
    if ( perdet ) {
      XLAL_CHECK( per_detectors->length == ( *out )->par.per_detectors->length, XLAL_EIO, "Inconsistent number of detectors: %u != %u", per_detectors->length, ( *out )->par.per_detectors->length );
      for ( size_t i = 0; i < per_detectors->length; ++i ) {
        XLAL_CHECK( strcmp( per_detectors->data[i], ( *out )->par.per_detectors->data[i] ) == 0, XLAL_EIO, "Inconsistent detectors: %s != %s", per_detectors->data[i], ( *out )->par.per_detectors->data[i] );
      }
    }

    // Check if outputting per-segment quantities, and number of per-segment items
    XLAL_CHECK( ( !perseg ) == ( ( *out )->par.per_nsegments == 0 ), XLAL_EIO, "Inconsistent output per segment flag vs number of segments: %i != %i", perseg, ( *out )->par.per_nsegments > 0 );
    if ( perseg ) {
      XLAL_CHECK( (size_t) per_nsegments == ( *out )->par.per_nsegments, XLAL_EIO, "Inconsistent number of segments: %i != %u", per_nsegments, ( *out )->par.per_nsegments );
    }

  }

  // Read and increment number of computed semicoherent results
  {
    UINT8 semi_ncomp = 0;
    XLAL_CHECK( XLALFITSHeaderReadUINT8( file, "semicomp", &semi_ncomp ) == XLAL_SUCCESS, XLAL_EFUNC );
    ( *out )->semi_ncomp += semi_ncomp;
  }

  // Read and append to toplist ranked by mean multi-detector F-statistic
  XLAL_CHECK( XLALWeaveResultsToplistReadAppend( file, ( *out )->mean_twoF_toplist ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Cleanup
  XLALDestroyStringVector( per_detectors );

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
  if ( XLALGPSCmp( &out_1->par.ref_time, &out_2->par.ref_time ) != 0 ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal reference times: %" LAL_GPS_FORMAT " != %" LAL_GPS_FORMAT "\n", __func__, LAL_GPS_PRINT( out_1->par.ref_time ), LAL_GPS_PRINT( out_2->par.ref_time ) );
    return XLAL_SUCCESS;
  }

  // Compare number of spindowns
  if ( out_1->par.nspins != out_2->par.nspins ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal number of spindowns: %zu != %zu\n", __func__, out_1->par.nspins, out_2->par.nspins );
    return XLAL_SUCCESS;
  }

  // Compare if outputting per-detector quantities, and list of detectors
  {
    const BOOLEAN perdet_1 = ( out_1->par.per_detectors != NULL ), perdet_2 = ( out_2->par.per_detectors != NULL );
    if ( perdet_1 != perdet_2 ) {
      *equal = 0;
      XLALPrintInfo( "%s: unequal output per detector?: %i != %i\n", __func__, perdet_1, perdet_2 );
      return XLAL_SUCCESS;
    }
    if ( perdet_1 ) {
      if ( out_1->par.per_detectors->length != out_2->par.per_detectors->length ) {
        *equal = 0;
        XLALPrintInfo( "%s: unequal number of detectors: %u != %u\n", __func__, out_1->par.per_detectors->length, out_2->par.per_detectors->length );
        return XLAL_SUCCESS;
      }
      for ( size_t i = 0; i < out_1->par.per_detectors->length; ++i ) {
        if ( strcmp( out_1->par.per_detectors->data[i], out_2->par.per_detectors->data[i] ) != 0 ) {
          *equal = 0;
          XLALPrintInfo( "%s: unequal detectors: %s != %s\n", __func__, out_1->par.per_detectors->data[i], out_2->par.per_detectors->data[i] );
          return XLAL_SUCCESS;
        }
      }
    }
  }

  // Compare number of per-segment items
  if ( out_1->par.per_nsegments != out_2->par.per_nsegments ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal number of segments: %u != %u\n", __func__, out_1->par.per_nsegments, out_2->par.per_nsegments );
    return XLAL_SUCCESS;
  }

  // Compare number of computed semicoherent results
  if ( out_1->semi_ncomp != out_2->semi_ncomp ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal number of computed semicoherent results: %" LAL_INT8_FORMAT " != %" LAL_INT8_FORMAT "\n", __func__, out_1->semi_ncomp, out_2->semi_ncomp );
    return XLAL_SUCCESS;
  }

  // Compare toplists ranked by mean multi-detector F-statistic
  XLAL_CHECK( XLALWeaveResultsToplistCompare( equal, setup, param_tol_mism, result_tol, out_1->mean_twoF_toplist, out_2->mean_twoF_toplist ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( !*equal ) {
    return XLAL_SUCCESS;
  }

  return XLAL_SUCCESS;

}

///
/// Write miscellaneous per-segment information to a FITS file
///
int XLALWeaveOutputMiscPerSegInfoWrite(
  FITSFile *file,
  const WeaveSetupData *setup,
  const BOOLEAN write_SFT_info,
  const UINT4 nsegments,
  const WeaveOutputMiscPerSegInfo *per_seg_info
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup != NULL, XLAL_EFAULT );
  XLAL_CHECK( nsegments > 0, XLAL_EINVAL );
  XLAL_CHECK( per_seg_info != NULL, XLAL_EFAULT );

  // Begin FITS table
  XLAL_CHECK( XLALFITSTableOpenWrite( file, "per_seg_info", "miscellaneous per-segment information" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Describe FITS table
  {
    char col_name[32];
    XLAL_FITS_TABLE_COLUMN_BEGIN( WeaveOutputMiscPerSegInfo );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, segment_start ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, segment_end ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( write_SFT_info ) {
      for ( size_t i = 0; i < setup->detectors->length; ++i ) {
        snprintf( col_name, sizeof( col_name ), "sft_first_%s", setup->detectors->data[i] );
        XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, GPSTime, sft_first[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
        snprintf( col_name, sizeof( col_name ), "sft_last_%s", setup->detectors->data[i] );
        XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, GPSTime, sft_last[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
        snprintf( col_name, sizeof( col_name ), "sft_count_%s", setup->detectors->data[i] );
        XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, INT4, sft_count[i], col_name ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, sft_min_cover_freq ) == XLAL_SUCCESS, XLAL_EFUNC );
      XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, sft_max_cover_freq ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, INT4, coh_n1comp ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, INT4, coh_nrecomp ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write FITS table
  for ( size_t i = 0; i < nsegments; ++i ) {
    XLAL_CHECK( XLALFITSTableWriteRow( file, &per_seg_info[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}
