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
#include "ResultsBasket.h"

///
/// Internal definition of output results from a search
///
struct tagWeaveOutputResults {
  /// Various varameters required to output results
  WeaveOutputParams par;
  /// Total number of semicoherent results added to output
  INT8 semi_total;
  /// Basket of output results ranked by mean multi-detector F-statistic
  WeaveResultsBasket *mean_twoF_basket;
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

  // Allocate memory
  WeaveOutputResultItem *item = XLALCalloc( 1, sizeof( *item ) );
  XLAL_CHECK_NULL( item != NULL, XLAL_ENOMEM );
  if ( par->per_nsegments > 0 ) {
    item->per_seg = XLALCalloc( par->per_nsegments, sizeof( *item->per_seg ) );
    XLAL_CHECK_NULL( item->per_seg != NULL, XLAL_ENOMEM );
  }

  // Set reference time of physical coordinates
  item->semi_phys.refTime = par->ref_time;
  for ( size_t j = 0; j < par->per_nsegments; ++j ) {
    item->per_seg[j].coh_phys.refTime = par->ref_time;
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
    XLALFree( item->per_seg );
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

  // Create a basket which ranks results by mean multi-detector F-statistic
  out->mean_twoF_basket = XLALWeaveResultsBasketCreate( &out->par, "mean_twoF", "mean multi-detector F-statistic", toplist_limit, result_item_compare_by_mean_twoF );
  XLAL_CHECK_NULL( out->mean_twoF_basket != NULL, XLAL_EFUNC );

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
    XLALWeaveResultsBasketDestroy( out->mean_twoF_basket );
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

  // Increment total number of semicoherent results
  out->semi_total += semi_nfreqs;

  // Add results to basket ranked by mean multi-detector F-statistic
  XLAL_CHECK( XLALWeaveResultsBasketAdd( out->mean_twoF_basket, semi_res, semi_nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );

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
  XLAL_CHECK( XLALFITSHeaderWriteINT4( file, "nspins", out->par.nspins, "number of spindowns" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write if outputting per-detector quantities
  XLAL_CHECK( XLALFITSHeaderWriteBOOLEAN( file, "perdet", out->par.per_detectors != NULL, "output per detector?" ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( out->par.per_detectors != NULL ) {
    XLAL_CHECK( XLALFITSHeaderWriteStringVector( file, "detect", out->par.per_detectors, "setup detectors" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write if outputting per-segment quantities
  XLAL_CHECK( XLALFITSHeaderWriteBOOLEAN( file, "perseg", out->par.per_nsegments > 0, "output per segment?" ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( out->par.per_nsegments > 0 ) {
    XLAL_CHECK( XLALFITSHeaderWriteINT4( file, "nsegment", out->par.per_nsegments, "number of segments" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write total number of semicoherent results added to output
  XLAL_CHECK( XLALFITSHeaderWriteINT8( file, "semitot", out->semi_total, "total semicoherent templates searched" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write basket ranked by mean multi-detector F-statistic
  XLAL_CHECK( XLALWeaveResultsBasketWrite( file, out->mean_twoF_basket ) == XLAL_SUCCESS, XLAL_EFUNC );

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
  INT4 nspins = 0;
  XLAL_CHECK( XLALFITSHeaderReadINT4( file, "nspins", &nspins ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Read if outputting per-detector quantities, and list of detectors
  BOOLEAN perdet = 0;
  LALStringVector *per_detectors = NULL;
  XLAL_CHECK( XLALFITSHeaderReadBOOLEAN( file, "perdet", &perdet ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( perdet ) {
    XLAL_CHECK( XLALFITSHeaderReadStringVector( file, "detect", &per_detectors ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Read if outputting per-segment quantities, and number of per-segment items
  BOOLEAN perseg = 0;
  INT4 per_nsegments = 0;
  XLAL_CHECK( XLALFITSHeaderReadBOOLEAN( file, "perseg", &perseg ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( perseg ) {
    XLAL_CHECK( XLALFITSHeaderReadINT4( file, "nsegment", &per_nsegments ) == XLAL_SUCCESS, XLAL_EFUNC );
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
    XLAL_CHECK( !perdet == ( ( *out )->par.per_detectors == NULL ), XLAL_EIO, "Inconsistent output per detector flag vs list of detectors: %i != %i", perseg, ( *out )->par.per_detectors != NULL );
    if ( perdet ) {
      XLAL_CHECK( per_detectors->length == ( *out )->par.per_detectors->length, XLAL_EIO, "Inconsistent number of detectors: %u != %u", per_detectors->length, ( *out )->par.per_detectors->length );
      for ( size_t i = 0; i < per_detectors->length; ++i ) {
        XLAL_CHECK( strcmp( per_detectors->data[i], ( *out )->par.per_detectors->data[i] ) == 0, XLAL_EIO, "Inconsistent detectors: %s != %s", per_detectors->data[i], ( *out )->par.per_detectors->data[i] );
      }
    }

    // Check if outputting per-segment quantities, and number of per-segment items
    XLAL_CHECK( !perseg == ( ( *out )->par.per_nsegments == 0 ), XLAL_EIO, "Inconsistent output per segment flag vs number of segments: %i != %i", perseg, ( *out )->par.per_nsegments > 0 );
    if ( perseg ) {
      XLAL_CHECK( (size_t) per_nsegments == ( *out )->par.per_nsegments, XLAL_EIO, "Inconsistent number of segments: %i != %u", per_nsegments, ( *out )->par.per_nsegments );
    }

  }

  // Read and increment total number of semicoherent results added to output
  {
    INT8 semi_total = 0;
    XLAL_CHECK( XLALFITSHeaderReadINT8( file, "semitot", &semi_total ) == XLAL_SUCCESS, XLAL_EFUNC );
    ( *out )->semi_total += semi_total;
  }

  // Read and append to basket ranked by mean multi-detector F-statistic
  XLAL_CHECK( XLALWeaveResultsBasketReadAppend( file, ( *out )->mean_twoF_basket ) == XLAL_SUCCESS, XLAL_EFUNC );

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

  // Compare total number of semicoherent results
  if ( out_1->semi_total != out_2->semi_total ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal total number of semicoherent results: %" LAL_INT8_FORMAT " != %" LAL_INT8_FORMAT "\n", __func__, out_1->semi_total, out_2->semi_total );
    return XLAL_SUCCESS;
  }

  // Compare baskets ranked by mean multi-detector F-statistic
  XLAL_CHECK( XLALWeaveResultsBasketCompare( equal, setup, param_tol_mism, result_tol, out_1->mean_twoF_basket, out_2->mean_twoF_basket ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( !*equal ) {
    return XLAL_SUCCESS;
  }

  return XLAL_SUCCESS;

}
