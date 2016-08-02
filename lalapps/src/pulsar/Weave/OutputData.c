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

#include "Weave.h"

#include <lal/LALStdio.h>
#include <lal/LALString.h>
#include <lal/LALHeap.h>
#include <lal/UserInput.h>

///
/// Internal definition of output data from a search
///
struct tagWeaveOutput {
  /// Reference time at which search is conducted
  LIGOTimeGPS ref_time;
  /// Number of spindown parameters to output
  size_t nspins;
  /// If outputting per-detector quantities, list of detectors
  const LALStringVector *per_detectors;
  /// Number of per-segment items to output (may be zero)
  UINT4 per_nsegments;
  /// Total number of semicoherent results added to output
  INT8 semi_total;
  /// Toplist ranked by mean multi-detector F-statistic
  LALHeap *toplist_mean_twoF;
  /// Save a no-longer-used output data item for re-use
  WeaveOutputDataItem *saved_item;
};

///
/// Initialise a FITS table for writing/reading a toplist
///
static int toplist_fits_table_init(
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
  XLAL_FITS_TABLE_COLUMN_BEGIN( WeaveOutputDataItem );

  // Add columns for semicoherent template parameters
  XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_par.alpha, "alpha [rad]" );
  XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_par.delta, "delta [rad]" );
  XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_par.fkdot[0], "freq [Hz]" );
  for ( size_t k = 1; k <= nspins; ++k ) {
    snprintf( col_name, sizeof( col_name ), "f%zudot [Hz/s^%zu]", k, k );
    XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL8, semi_par.fkdot[k], col_name );
  }

  // Add columns for mean multi- and per-detector F-statistic
  XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, mean_twoF );
  if ( per_detectors != NULL ) {
    for ( size_t i = 0; i < per_detectors->length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "mean_twoF_%s", per_detectors->data[i] );
      XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, REAL4, mean_twoF_per_det[i], col_name );
    }
  }

  // Begin FITS table description for per-segment items (optional)
  if ( per_nsegments > 0 ) {
    XLAL_FITS_TABLE_COLUMN_PTR_BEGIN( per_seg, WeaveOutputDataPerSegItem, per_nsegments );
    for ( size_t s = 0; s < per_nsegments; ++s ) {

      // Add columns for coherent template parameters
      snprintf( col_name, sizeof( col_name ), "seg%zu_alpha [rad]", s + 1 );
      XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_par.alpha, col_name );
      snprintf( col_name, sizeof( col_name ), "seg%zu_delta [rad]", s + 1 );
      XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_par.delta, col_name );
      snprintf( col_name, sizeof( col_name ), "seg%zu_freq [Hz]", s + 1 );
      XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_par.fkdot[0], col_name );
      for ( size_t k = 1; k <= nspins; ++k ) {
        snprintf( col_name, sizeof( col_name ), "seg%zu_f%zudot [Hz/s^%zu]", s + 1, k, k );
        XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL8, coh_par.fkdot[k], col_name );
      }

      // Add columns for coherent multi- and per-detector F-statistic
      snprintf( col_name, sizeof( col_name ), "seg%zu_twoF", s + 1 );
      XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL4, twoF, col_name );
      if ( per_detectors != NULL ) {
        for ( size_t i = 0; i < per_detectors->length; ++i ) {
          snprintf( col_name, sizeof( col_name ), "seg%zu_twoF_%s", s + 1, per_detectors->data[i] );
          XLAL_FITS_TABLE_COLUMN_PTR_ADD_NAMED( file, s, REAL4, twoF_per_det[i], col_name );
        }
      }

    }
  }

  return XLAL_SUCCESS;

}

///
/// Visitor function for writing a toplist to a FITS table
///
static int toplist_fits_write_visitor(
  void *param,
  const void *x
  )
{
  FITSFile *file = ( FITSFile * ) param;
  XLAL_CHECK( XLALFITSTableWriteRow( file, x ) == XLAL_SUCCESS, XLAL_EFUNC );
  return XLAL_SUCCESS;
}

///
/// Destroy an output data item
///
static void output_data_item_destroy(
  void *x
  )
{
  if ( x != NULL ) {
    WeaveOutputDataItem *ix = ( WeaveOutputDataItem * ) x;
    XLALFree( ix->per_seg );
    XLALFree( ix );
  }
}

///
/// Compare output data items by mean multi-detector F-statistic
///
static int output_data_item_compare_by_mean_twoF(
  const void *x,
  const void *y
  )
{
  const WeaveOutputDataItem *ix = ( const WeaveOutputDataItem * ) x;
  const WeaveOutputDataItem *iy = ( const WeaveOutputDataItem * ) y;
  if ( ix->mean_twoF > iy->mean_twoF ) {
    return -1;
  }
  if ( ix->mean_twoF < iy->mean_twoF ) {
    return +1;
  }
  return 0;
}

///
/// Fill a output data item, creating a new one if needed
///
static int output_data_item_add(
  BOOLEAN *full_init,
  WeaveOutput *out,
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

  // Fill output data item, creating a new one if needed
  XLAL_CHECK( XLALWeaveFillOutputDataItem( &out->saved_item, full_init, semi_res, freq_idx ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Add item to toplist
  const WeaveOutputDataItem *prev_item = out->saved_item;
  XLAL_CHECK( XLALHeapAdd( toplist, ( void ** ) &out->saved_item ) == XLAL_SUCCESS, XLAL_EFUNC );

  // If toplist added item and returned a different, no-unused item,
  // that item will have to be fully initialised at the next call
  *full_init = ( out->saved_item != prev_item );

  return XLAL_SUCCESS;

}

///
/// Create output data
///
WeaveOutput *XLALWeaveOutputCreate(
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
  WeaveOutput *out = XLALCalloc( 1, sizeof( *out ) );
  XLAL_CHECK_NULL( out != NULL, XLAL_ENOMEM );

  // Set fields
  out->ref_time = *ref_time;
  out->nspins = nspins;
  out->per_detectors = per_detectors;
  out->per_nsegments = per_nsegments;
  out->semi_total = 0;

  // Create a toplist ranked by mean multi-detector F-statistic
  out->toplist_mean_twoF = XLALHeapCreate( output_data_item_destroy, toplist_limit, +1, output_data_item_compare_by_mean_twoF );
  XLAL_CHECK_NULL( out->toplist_mean_twoF != NULL, XLAL_EFUNC );

  return out;

}

///
/// Free output data
///
void XLALWeaveOutputDestroy(
  WeaveOutput *out
  )
{
  if ( out != NULL ) {
    output_data_item_destroy( out->saved_item );
    XLALHeapDestroy( out->toplist_mean_twoF );
    XLALFree( out );
  }
}

///
/// Add semicoherent results to output
///
int XLALWeaveOutputAdd(
  WeaveOutput *out,
  const WeaveSemiResults *semi_res,
  const UINT4 semi_nfreqs
  )
{

  // Check input
  XLAL_CHECK( out != NULL, XLAL_EFAULT );
  XLAL_CHECK( semi_res != NULL, XLAL_EFAULT );

  // Must initialise all output data item fields the first time
  BOOLEAN full_init = 1;

  // Iterate over the frequency bins of the semicoherent results
  for ( size_t i = 0; i < semi_nfreqs; ++i ) {

    // Add item to toplist ranked by mean multi-detector F-statistic
    XLAL_CHECK( output_data_item_add( &full_init, out, out->toplist_mean_twoF, semi_res, i ) == XLAL_SUCCESS, XLAL_EFUNC );

  }

  // Increment total number of semicoherent results
  out->semi_total += semi_nfreqs;

  return XLAL_SUCCESS;

}

///
/// Write output data to a FITS file
///
int XLALWeaveOutputWrite(
  FITSFile *file,
  const WeaveOutput *out
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( out != NULL, XLAL_EFAULT );

  // Write command line
  {
    char *command_line = XLALUserVarGetLog( UVAR_LOGFMT_CMDLINE );
    XLAL_CHECK( command_line != NULL, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteHistory( file, "Command line: %s", command_line ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALFree( command_line );
  }

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
  XLAL_CHECK( XLALFITSHeaderWriteBOOLEAN( file, "perdet", out->per_nsegments > 0, "output per segment?" ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( out->per_nsegments > 0 ) {
    XLAL_CHECK( XLALFITSHeaderWriteINT4( file, "nsegment", out->per_nsegments, "number of segments" ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write total number of semicoherent results added to output
  XLAL_CHECK( XLALFITSHeaderWriteINT8( file, "semitot", out->semi_total, "total semicoherent templates searched" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write toplist ranked by mean multi-detector F-statistic
  XLAL_CHECK( XLALFITSTableOpenWrite( file, "toplist_mean_twoF", "toplist ranked by mean multi-detector F-statistic" ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( toplist_fits_table_init( file, out->nspins, out->per_detectors, out->per_nsegments ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALHeapVisit( out->toplist_mean_twoF, toplist_fits_write_visitor, file ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

}

///
/// Write output detailed information per segment to a FITS file
///
int XLALWeaveOutputWriteDetails(
  FITSFile *file,
  const LALStringVector *detectors,
  const size_t nsegments,
  WeaveOutputDetails details[]
  )
{

  // Check input
  XLAL_CHECK( file != NULL, XLAL_EFAULT );
  XLAL_CHECK( detectors != NULL, XLAL_EFAULT );
  XLAL_CHECK( nsegments > 0, XLAL_EINVAL );
  XLAL_CHECK( details != NULL, XLAL_EFAULT );

  // Begin FITS table
  XLAL_CHECK( XLALFITSTableOpenWrite( file, "details", "details segment information" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Describe FITS table
  {
    char col_name[32];
    XLAL_FITS_TABLE_COLUMN_BEGIN( WeaveOutputDetails );
    XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, segment_start );
    XLAL_FITS_TABLE_COLUMN_ADD( file, GPSTime, segment_end );
    for ( size_t i = 0; i < detectors->length; ++i ) {
      snprintf( col_name, sizeof( col_name ), "sft_first_%s", detectors->data[i] );
      XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, GPSTime, sft_first[i], col_name );
      snprintf( col_name, sizeof( col_name ), "sft_last_%s", detectors->data[i] );
      XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, GPSTime, sft_last[i], col_name );
      snprintf( col_name, sizeof( col_name ), "sft_count_%s", detectors->data[i] );
      XLAL_FITS_TABLE_COLUMN_ADD_NAMED( file, INT4, sft_count[i], col_name );
    }
    XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, min_cover_freq );
    XLAL_FITS_TABLE_COLUMN_ADD( file, REAL8, max_cover_freq );
    XLAL_FITS_TABLE_COLUMN_ADD( file, INT4, coh_total );
    XLAL_FITS_TABLE_COLUMN_ADD( file, INT4, coh_total_recomp );
  }

  // Write output details for each segment
  for ( size_t i = 0; i < nsegments; ++i ) {
    XLAL_CHECK( XLALFITSTableWriteRow( file, &details[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

}
