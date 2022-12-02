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
// Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
// MA 02110-1301 USA
//

///
/// \file
/// \ingroup lalpulsar_bin_Weave
///

#include "OutputResults.h"
#include "ResultsToplist.h"
#include "Statistics.h"

#include <lal/UserInputPrint.h>

///
/// Histogram bin width of mean multi-F-statistics
///
const REAL4 mean2F_hgrm_bin_width = 0.1;

///
/// Histogram bin of mean multi-F-statistics
///
typedef struct {
  // Lower bin boundary
  REAL4 lower;
  // Upper bin boundary
  REAL4 upper;
  // Bin count
  UINT8 count;
} WeaveMean2FHistogramBin;

///
/// Output results from a search
///
struct tagWeaveOutputResults {
  /// Struct holding all parameters for which statistics to output and compute, when, and how
  /// NOTE: this is the *owner* of WeaveStatisticsParams, which is where it will be freed at the end
  /// while toplists will simply hold a reference-pointer
  WeaveStatisticsParams *statistics_params;
  /// Reference time at which search is conducted
  LIGOTimeGPS ref_time;
  /// Number of spindown parameters to output
  size_t nspins;
  // Maximum size of toplists
  UINT4 toplist_limit;
  /// Number of output results toplists
  size_t ntoplists;
  /// Output result toplists
  WeaveResultsToplist *toplists[8];
  // Vector to store histogram of mean multi-F-statistics
  UINT8Vector* mean2F_hgrm_bins;
  // Number of mean multi-F-statistics below range of histogram
  UINT8 mean2F_hgrm_underflow;
  // Number of mean multi-F-statistics above range of histogram
  UINT8 mean2F_hgrm_overflow;
  // Temporary REAL4 vector for generating histogram of mean multi-F-statistics
  REAL4Vector* mean2F_hgrm_tmp_REAL4;
  // Temporary INT4 vector for generating histogram of mean multi-F-statistics
  INT4Vector* mean2F_hgrm_tmp_INT4;
};

///
/// \name Internal functions
///
/// @{

/// @}

///
/// \name Functions for results toplist ranked by summed multi-detector F-statistic
///
/// @{

static const REAL4 *toplist_results_sum2F( const WeaveSemiResults *semi_res )
{
  return semi_res->sum2F->data;
}
static REAL4 toplist_item_get_sum2F( const WeaveResultsToplistItem *item )
{
  return item->stage[0].sum2F;
}
static void toplist_item_set_sum2F( WeaveResultsToplistItem *item, const REAL4 value )
{
  item->stage[0].sum2F = value;
}

/// @}
///
/// \name Functions for results toplist ranked by mean multi-detector F-statistic
///
/// @{

static const REAL4 *toplist_results_mean2F( const WeaveSemiResults *semi_res )
{
  return semi_res->mean2F->data;
}
static REAL4 toplist_item_get_mean2F( const WeaveResultsToplistItem *item )
{
  return item->stage[0].mean2F;
}
static void toplist_item_set_mean2F( WeaveResultsToplistItem *item, const REAL4 value )
{
  item->stage[0].mean2F = value;
}

/// @}
///
/// \name Functions for results toplist ranked by line-robust log10(B_S/GL) statistic
///
/// @{

static const REAL4 *toplist_results_log10BSGL( const WeaveSemiResults *semi_res )
{
  return semi_res->log10BSGL->data;
}
static REAL4 toplist_item_get_log10BSGL( const WeaveResultsToplistItem *item )
{
  return item->stage[0].log10BSGL;
}
static void toplist_item_set_log10BSGL( WeaveResultsToplistItem *item, const REAL4 value )
{
  item->stage[0].log10BSGL = value;
}

/// @}
///
/// \name Functions for results toplist ranked by transient-line-robust log10(B_S/GLtL) statistic
///
/// @{

static const REAL4 *toplist_results_log10BSGLtL( const WeaveSemiResults *semi_res )
{
  return semi_res->log10BSGLtL->data;
}
static REAL4 toplist_item_get_log10BSGLtL( const WeaveResultsToplistItem *item )
{
  return item->stage[0].log10BSGLtL;
}
static void toplist_item_set_log10BSGLtL( WeaveResultsToplistItem *item, const REAL4 value )
{
  item->stage[0].log10BSGLtL = value;
}

/// @}
///
/// \name Functions for results toplist ranked by transient-signal line-robust log10(B_tS/GLtL) statistic
///
/// @{

static const REAL4 *toplist_results_log10BtSGLtL( const WeaveSemiResults *semi_res )
{
  return semi_res->log10BtSGLtL->data;
}
static REAL4 toplist_item_get_log10BtSGLtL( const WeaveResultsToplistItem *item )
{
  return item->stage[0].log10BtSGLtL;
}
static void toplist_item_set_log10BtSGLtL( WeaveResultsToplistItem *item, const REAL4 value )
{
  item->stage[0].log10BtSGLtL = value;
}

/// @}

///
/// Create output results
///
WeaveOutputResults *XLALWeaveOutputResultsCreate(
  const LIGOTimeGPS *ref_time,
  const size_t nspins,
  WeaveStatisticsParams *statistics_params,
  const UINT4 toplist_limit,
  const BOOLEAN mean2F_hgrm
  )
{
  // Check input
  XLAL_CHECK_NULL( ref_time != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( statistics_params != NULL, XLAL_EFAULT );

  // Allocate memory
  WeaveOutputResults *out = XLALCalloc( 1, sizeof( *out ) );
  XLAL_CHECK_NULL( out != NULL, XLAL_ENOMEM );

  // Set fields
  out->ref_time = *ref_time;
  out->nspins = nspins;
  out->toplist_limit = toplist_limit;
  out->statistics_params = statistics_params;

  WeaveStatisticType toplist_statistics = statistics_params->toplist_statistics;

  // Initialise number of toplists
  out->ntoplists = 0;

  // Create a toplist which ranks results by mean multi-detector F-statistic
  if ( toplist_statistics & WEAVE_STATISTIC_MEAN2F ) {
    out->toplists[out->ntoplists] = XLALWeaveResultsToplistCreate( nspins, statistics_params, WEAVE_STATISTIC_NAME( WEAVE_STATISTIC_MEAN2F ), "average multi-detector F-statistic", toplist_limit, toplist_results_mean2F, toplist_item_get_mean2F, toplist_item_set_mean2F );
    XLAL_CHECK_NULL( out->toplists[out->ntoplists] != NULL, XLAL_EFUNC );
    XLAL_CHECK_NULL( out->ntoplists < XLAL_NUM_ELEM( out->toplists ), XLAL_EFAILED );
    out->ntoplists++;
  }

  // Create a toplist which ranks results by summed multi-detector F-statistic
  if ( toplist_statistics & WEAVE_STATISTIC_SUM2F ) {
    out->toplists[out->ntoplists] = XLALWeaveResultsToplistCreate( nspins, statistics_params, WEAVE_STATISTIC_NAME( WEAVE_STATISTIC_SUM2F ), "summed multi-detector F-statistic", toplist_limit, toplist_results_sum2F, toplist_item_get_sum2F, toplist_item_set_sum2F );
    XLAL_CHECK_NULL( out->toplists[out->ntoplists] != NULL, XLAL_EFUNC );
    XLAL_CHECK_NULL( out->ntoplists < XLAL_NUM_ELEM( out->toplists ), XLAL_EFAILED );
    out->ntoplists++;
  }

  // Create a toplist which ranks results by line-robust log10(B_S/GL) statistic
  if ( toplist_statistics & WEAVE_STATISTIC_BSGL ) {
    out->toplists[out->ntoplists] = XLALWeaveResultsToplistCreate( nspins, statistics_params, WEAVE_STATISTIC_NAME( WEAVE_STATISTIC_BSGL ), "line-robust log10BSGL statistic", toplist_limit, toplist_results_log10BSGL, toplist_item_get_log10BSGL, toplist_item_set_log10BSGL );
    XLAL_CHECK_NULL( out->toplists[out->ntoplists] != NULL, XLAL_EFUNC );
    XLAL_CHECK_NULL( out->ntoplists < XLAL_NUM_ELEM( out->toplists ), XLAL_EFAILED );
    out->ntoplists++;
  }

  // Create a toplist which ranks results by transient-line-robust log10(B_S/GLtL) statistic
  if ( toplist_statistics & WEAVE_STATISTIC_BSGLtL ) {
    out->toplists[out->ntoplists] = XLALWeaveResultsToplistCreate( nspins, statistics_params, WEAVE_STATISTIC_NAME( WEAVE_STATISTIC_BSGLtL ), "transient line-robust log10BSGLtL statistic", toplist_limit, toplist_results_log10BSGLtL, toplist_item_get_log10BSGLtL, toplist_item_set_log10BSGLtL );
    XLAL_CHECK_NULL( out->toplists[out->ntoplists] != NULL, XLAL_EFUNC );
    XLAL_CHECK_NULL( out->ntoplists < XLAL_NUM_ELEM( out->toplists ), XLAL_EFAILED );
    out->ntoplists++;
  }

  // Create a toplist which ranks results by transient-signal line-robust log10(B_tS/GLtL) statistic
  if ( toplist_statistics & WEAVE_STATISTIC_BtSGLtL ) {
    out->toplists[out->ntoplists] = XLALWeaveResultsToplistCreate( nspins, statistics_params, WEAVE_STATISTIC_NAME( WEAVE_STATISTIC_BtSGLtL ), "transient signal line-robust log10BtSGLtL statistic", toplist_limit, toplist_results_log10BtSGLtL, toplist_item_get_log10BtSGLtL, toplist_item_set_log10BtSGLtL );
    XLAL_CHECK_NULL( out->toplists[out->ntoplists] != NULL, XLAL_EFUNC );
    XLAL_CHECK_NULL( out->ntoplists < XLAL_NUM_ELEM( out->toplists ), XLAL_EFAILED );
    out->ntoplists++;
  }

  // Comnsistency check on number of toplists
  XLAL_CHECK_NULL( out->ntoplists == statistics_params->ntoplists, XLAL_EFAILED );

  // Create histogram of mean multi-F-statistic
  if ( mean2F_hgrm ) {
    out->mean2F_hgrm_bins = XLALCreateUINT8Vector( 10000 );
    XLAL_CHECK_NULL( out->mean2F_hgrm_bins != NULL, XLAL_EFUNC );
    memset( out->mean2F_hgrm_bins->data, 0, sizeof( out->mean2F_hgrm_bins->data[0] ) * out->mean2F_hgrm_bins->length );
  }

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
    XLALWeaveStatisticsParamsDestroy( out->statistics_params );
    for ( size_t i = 0; i < out->ntoplists; ++i ) {
      XLALWeaveResultsToplistDestroy( out->toplists[i] );
    }
    XLALDestroyUINT8Vector( out->mean2F_hgrm_bins );
    XLALDestroyREAL4Vector( out->mean2F_hgrm_tmp_REAL4 );
    XLALDestroyINT4Vector( out->mean2F_hgrm_tmp_INT4 );
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

  // Add results to toplists
  for ( size_t i = 0; i < out->ntoplists; ++i ) {
    XLAL_CHECK( XLALWeaveResultsToplistAdd( out->toplists[i], semi_res, semi_nfreqs ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Add to histogram of mean multi-F-statistics
  if ( out->mean2F_hgrm_bins != NULL ) {

    // Check input
    XLAL_CHECK( semi_res->mean2F != NULL, XLAL_EINVAL );

    // Allocate memory
    if ( out->mean2F_hgrm_tmp_REAL4 == NULL || out->mean2F_hgrm_tmp_REAL4->length < semi_res->mean2F->length ) {
      out->mean2F_hgrm_tmp_REAL4 = XLALResizeREAL4Vector( out->mean2F_hgrm_tmp_REAL4, semi_res->mean2F->length );
      XLAL_CHECK( out->mean2F_hgrm_tmp_REAL4 != NULL, XLAL_EFUNC );
      out->mean2F_hgrm_tmp_INT4 = XLALResizeINT4Vector( out->mean2F_hgrm_tmp_INT4, semi_res->mean2F->length );
      XLAL_CHECK( out->mean2F_hgrm_tmp_INT4 != NULL, XLAL_EFUNC );
    }

    // Put mean multi-F-statistics into bins
    XLAL_CHECK( XLALVectorScaleREAL4( out->mean2F_hgrm_tmp_REAL4->data, 1.0 / mean2F_hgrm_bin_width, semi_res->mean2F->data, semi_res->mean2F->length ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALVectorINT4FromREAL4( out->mean2F_hgrm_tmp_INT4->data, out->mean2F_hgrm_tmp_REAL4->data, semi_res->mean2F->length ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Add to histogram bins
    for ( size_t j = 0; j < semi_res->mean2F->length; ++j ) {
      INT4 bin = out->mean2F_hgrm_tmp_INT4->data[j];
      const INT4 bin_max = out->mean2F_hgrm_bins->length;
      if ( bin < 0 ) {
        ++out->mean2F_hgrm_underflow;
      } else if ( bin >= bin_max ) {
        ++out->mean2F_hgrm_overflow;
      } else {
        ++out->mean2F_hgrm_bins->data[bin];
      }
    }

  }

  return XLAL_SUCCESS;

}

///
/// Compute all the missing 'completion-loop' statistics for all toplist entries
///
int XLALWeaveOutputResultsCompletionLoop(
  WeaveOutputResults *out
  )
{
  // Check input
  XLAL_CHECK( out != NULL, XLAL_EFAULT );

  // Iterate over all toplists
  for ( size_t i = 0; i < out->ntoplists; ++i ) {
    XLAL_CHECK( XLALWeaveResultsToplistCompletionLoop( out->toplists[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

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
  XLAL_CHECK( XLALFITSHeaderWriteUINT4( file, "nspins", out->nspins, "number of spindowns" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write list of detectors (if outputting per-detector quantities
  XLAL_CHECK( XLALFITSHeaderWriteStringVector( file, "detect", out->statistics_params->detectors, "list of detectors" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write number of segments
  XLAL_CHECK( XLALFITSHeaderWriteUINT4( file, "segments", out->statistics_params->nsegments, "number of segments" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write names of selected toplist types (ie ranking statistics)
  {
    char *toplist_statistics = XLALPrintStringValueOfUserFlag( ( const int * )&( out->statistics_params->toplist_statistics ), &WeaveToplistChoices );
    XLAL_CHECK( toplist_statistics != NULL, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteString( file, "toplists", toplist_statistics, "names of selected toplist statistics" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALFree( toplist_statistics );
  }

  // Write names of all selected 'extra' output statistics
  {
    WeaveStatisticType extras = ( out->statistics_params->statistics_to_output[0] & ~out->statistics_params->toplist_statistics );
    if ( extras == 0 ) {
      extras = WEAVE_STATISTIC_NONE;
    }
    char *extras_names = XLALPrintStringValueOfUserFlag( ( const int * )&extras, &WeaveStatisticChoices );
    XLAL_CHECK( extras_names != NULL, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteString( file, "extras", extras_names, "names of additional selected output statistics" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALFree( extras_names );
  }
  // Write names of all requested 'recalc' statistics
  {
    char *recalc_names = XLALPrintStringValueOfUserFlag( ( const int * )&out->statistics_params->statistics_to_output[1], &WeaveStatisticChoices );
    XLAL_CHECK( recalc_names != NULL, XLAL_EFUNC );
    XLAL_CHECK( XLALFITSHeaderWriteString( file, "recalc", recalc_names, "names of selected recalc statistics" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALFree( recalc_names );
  }

  // Write maximum size of toplists
  XLAL_CHECK( XLALFITSHeaderWriteUINT4( file, "toplimit", out->toplist_limit, "maximum size of toplists" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write whether a histogram of mean multi-F-statistics will be written
  XLAL_CHECK( XLALFITSHeaderWriteBOOLEAN( file, "m2Fhgrm", out->mean2F_hgrm_bins != NULL, "mean 2F histogram?" ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Write toplists
  for ( size_t i = 0; i < out->ntoplists; ++i ) {
    XLAL_CHECK( XLALWeaveResultsToplistWrite( file, out->toplists[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Write histogram of mean multi-F-statistics
  if ( out->mean2F_hgrm_bins != NULL ) {

    // Open and describe FITS table for writing histogram bins
    XLAL_CHECK( XLALFITSTableOpenWrite( file, "mean2F_hgrm", "histogram of mean multi-F-statistics" ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_FITS_TABLE_COLUMN_BEGIN( WeaveMean2FHistogramBin );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, lower ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, upper ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, UINT8, count ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Write histogram bins
    WeaveMean2FHistogramBin bin;
    if ( out->mean2F_hgrm_underflow > 0 ) {
      bin.lower = GSL_NEGINF;
      bin.upper = 0;
      bin.count = out->mean2F_hgrm_underflow;
      XLAL_CHECK( XLALFITSTableWriteRow( file, &bin ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
    for ( size_t j = 0; j < out->mean2F_hgrm_bins->length; ++j ) {
      const UINT4 bin_count = out->mean2F_hgrm_bins->data[j];
      if ( bin_count > 0 ) {
        bin.lower = mean2F_hgrm_bin_width * j;
        bin.upper = mean2F_hgrm_bin_width * (j + 1);
        bin.count = bin_count;
        XLAL_CHECK( XLALFITSTableWriteRow( file, &bin ) == XLAL_SUCCESS, XLAL_EFUNC );
      }
    }
    if ( out->mean2F_hgrm_underflow > 0 ) {
      bin.lower = mean2F_hgrm_bin_width * out->mean2F_hgrm_bins->length;
      bin.upper = GSL_POSINF;
      bin.count = out->mean2F_hgrm_underflow;
      XLAL_CHECK( XLALFITSTableWriteRow( file, &bin ) == XLAL_SUCCESS, XLAL_EFUNC );
    }

  }

  return XLAL_SUCCESS;

}

///
/// Read results from a FITS file and append to new/existing output results
///
int XLALWeaveOutputResultsReadAppend(
  FITSFile *file,
  WeaveOutputResults **out,
  UINT4 toplist_limit
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

  // ----- read elements for 'statistics_params' struct ----------
  WeaveStatisticsParams *statistics_params = XLALCalloc( 1, sizeof( *statistics_params ) );
  XLAL_CHECK( statistics_params != NULL, XLAL_ENOMEM );

  // Read list of detectors
  XLAL_CHECK( XLALFITSHeaderReadStringVector( file, "detect", &( statistics_params->detectors ) ) == XLAL_SUCCESS, XLAL_EFUNC );

  /// Number of segments
  XLAL_CHECK( XLALFITSHeaderReadUINT4( file, "segments", &( statistics_params->nsegments ) ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Read names of selected toplist statistics
  char *toplist_stats_names = NULL;
  int toplist_stats = 0;
  XLAL_CHECK( XLALFITSHeaderReadString( file, "toplists", &toplist_stats_names ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALParseStringValueAsUserFlag( &toplist_stats, &WeaveToplistChoices, toplist_stats_names ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALFree( toplist_stats_names );

  // Read names of selected extra output stats
  char *extras_names = NULL;
  int extra_stats = 0;
  XLAL_CHECK( XLALFITSHeaderReadString( file, "extras", &extras_names ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLAL_CHECK( XLALParseStringValueAsUserFlag( &extra_stats, &WeaveStatisticChoices, extras_names ) == XLAL_SUCCESS, XLAL_EFUNC );
  XLALFree( extras_names );

  // Read names of selected recalc stats
  int recalc_stats = 0;
  BOOLEAN exists = 0;
  XLAL_CHECK( XLALFITSHeaderQueryKeyExists( file, "recalc" , &exists ) == XLAL_SUCCESS, XLAL_EFUNC );
  if ( exists ) {
    char *recalc_names = NULL;
    XLAL_CHECK( XLALFITSHeaderReadString( file, "recalc", &recalc_names ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLALParseStringValueAsUserFlag( &recalc_stats, &WeaveStatisticChoices, recalc_names ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLALFree( recalc_names );
  }

  // Compute and fill the full stats-dependency map
  XLAL_CHECK( XLALWeaveStatisticsParamsSetDependencyMap( statistics_params, toplist_stats, extra_stats, recalc_stats ) == XLAL_SUCCESS, XLAL_EFUNC );

  // Read whether a histogram of mean multi-F-statistics will be written
  BOOLEAN mean2F_hgrm = 0;
  {
    exists = 0;
    XLAL_CHECK( XLALFITSHeaderQueryKeyExists( file, "m2Fhgrm", &exists ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( exists ) {
      XLAL_CHECK( XLALFITSHeaderReadBOOLEAN( file, "m2Fhgrm", &mean2F_hgrm ) == XLAL_SUCCESS, XLAL_EFUNC );
    }
  }

  if ( *out == NULL ) {

    // Create new output results
    *out = XLALWeaveOutputResultsCreate( &ref_time, nspins, statistics_params, toplist_limit, mean2F_hgrm );
    XLAL_CHECK( *out != NULL, XLAL_EFUNC );

  } else {

    // Check reference time
    XLAL_CHECK( XLALGPSCmp( &ref_time, &( *out )->ref_time ) == 0, XLAL_EIO, "Inconsistent reference time: %" LAL_GPS_FORMAT " != %" LAL_GPS_FORMAT, LAL_GPS_PRINT( ref_time ), LAL_GPS_PRINT( ( *out )->ref_time ) );

    // Check number of spindowns
    XLAL_CHECK( ( size_t ) nspins == ( *out )->nspins, XLAL_EIO, "Inconsistent number of spindowns: %i != %zu", nspins, ( *out )->nspins );

    // Check if list of detectors agrees
    if ( statistics_params->detectors != NULL ) {
      XLAL_CHECK( statistics_params->detectors->length == ( *out )->statistics_params->detectors->length, XLAL_EIO, "Inconsistent number of detectors: %u != %u", statistics_params->detectors->length, ( *out )->statistics_params->detectors->length );
      for ( size_t i = 0; i < statistics_params->detectors->length; ++i ) {
        XLAL_CHECK( strcmp( statistics_params->detectors->data[i], ( *out )->statistics_params->detectors->data[i] ) == 0, XLAL_EIO, "Inconsistent detectors: %s != %s", statistics_params->detectors->data[i], ( *out )->statistics_params->detectors->data[i] );
      }
    }

    // Check if number of segments agrees
    XLAL_CHECK( statistics_params->nsegments == ( *out )->statistics_params->nsegments, XLAL_EIO, "Inconsistent output per segment?: %i != %u", statistics_params->nsegments, ( *out )->statistics_params->nsegments );

    // Check list of selected toplist statistics
    if ( statistics_params->toplist_statistics != ( *out )->statistics_params->toplist_statistics ) {
      char *toplists1, *toplists2;
      toplists1 = XLALPrintStringValueOfUserFlag( ( const int * )&( statistics_params->toplist_statistics ), &WeaveToplistChoices );
      XLAL_CHECK( toplists1 != NULL, XLAL_EFUNC );
      toplists2 = XLALPrintStringValueOfUserFlag( ( const int * )&( ( *out )->statistics_params->toplist_statistics ), &WeaveToplistChoices );
      XLAL_CHECK( toplists2 != NULL, XLAL_EFUNC );
      XLALPrintError( "Inconsistent set of toplist statistics: %s != %s\n", toplists1, toplists2 );
      XLALFree( toplists1 );
      XLALFree( toplists2 );
      XLAL_ERROR( XLAL_EIO );
    }
    // Check list of selected output statistics
    for ( UINT4 istage = 0; istage < 2; ++ istage ) {
      if ( statistics_params->statistics_to_output[istage] != ( *out )->statistics_params->statistics_to_output[istage] ) {
        char *output1, *output2;
        output1 = XLALPrintStringValueOfUserFlag( ( const int * )&( statistics_params->statistics_to_output[istage] ), &WeaveStatisticChoices );
        XLAL_CHECK( output1 != NULL, XLAL_EFUNC );
        output2 = XLALPrintStringValueOfUserFlag( ( const int * )&( ( *out )->statistics_params->statistics_to_output[istage] ), &WeaveStatisticChoices );
        XLAL_CHECK( output2 != NULL, XLAL_EFUNC );
        XLALPrintError( "Inconsistent set of stage-%d output statistics: {%s} != {%s}\n", istage, output1, output2 );
        XLALFree( output1 );
        XLALFree( output2 );
        XLAL_ERROR( XLAL_EIO );
      }
    }

    XLALWeaveStatisticsParamsDestroy( statistics_params );  // Not creating a new output, so we need to free this

    // Check whether a histogram of mean multi-F-statistics will be written
    XLAL_CHECK( !mean2F_hgrm == !(( *out )->mean2F_hgrm_bins != NULL), XLAL_EIO, "Inconsistent mean 2F histogram? %i != %i", mean2F_hgrm, ( *out )->mean2F_hgrm_bins != NULL );

  }

  // Read and append to toplists
  for ( size_t i = 0; i < ( *out )->ntoplists; ++i ) {
    XLAL_CHECK( XLALWeaveResultsToplistReadAppend( file, ( *out )->toplists[i] ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  // Read and append histogram of mean multi-F-statistics
  if ( ( *out )->mean2F_hgrm_bins != NULL ) {

    // Open and describe FITS table for writing histogram bins
    UINT8 nrows = 0;
    XLAL_CHECK( XLALFITSTableOpenRead( file, "mean2F_hgrm", &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_FITS_TABLE_COLUMN_BEGIN( WeaveMean2FHistogramBin );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, lower ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, REAL4, upper ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK( XLAL_FITS_TABLE_COLUMN_ADD( file, UINT8, count ) == XLAL_SUCCESS, XLAL_EFUNC );

    // Read histogram bins
    WeaveMean2FHistogramBin bin;
    const REAL4 bin_upper_max = mean2F_hgrm_bin_width * ( *out )->mean2F_hgrm_bins->length;
    while ( nrows > 0 ) {
      XLAL_CHECK( XLALFITSTableReadRow( file, &bin, &nrows ) == XLAL_SUCCESS, XLAL_EFUNC );
      if ( bin.lower < 0 ) {
        ( *out )->mean2F_hgrm_underflow += bin.count;
      } else if ( bin.upper >= bin_upper_max ) {
        ( *out )->mean2F_hgrm_overflow += bin.count;
      } else {
        const size_t j = bin.lower / mean2F_hgrm_bin_width;
        ( *out )->mean2F_hgrm_bins->data[j] += bin.count;
      }
    }

  }

  return XLAL_SUCCESS;

}

///
/// Compare two output results and return whether they are equal
///
int XLALWeaveOutputResultsCompare(
  BOOLEAN *equal,
  const WeaveSetupData *setup,
  const BOOLEAN sort_by_semi_phys,
  const UINT4 round_param_to_dp,
  const UINT4 round_param_to_sf,
  const REAL8 unmatched_item_tol,
  const REAL8 param_tol_mism,
  const VectorComparison *result_tol,
  const UINT4 toplist_compare_limit,
  const WeaveOutputResults *out_1,
  const WeaveOutputResults *out_2
  )
{

  // Check input
  XLAL_CHECK( equal != NULL, XLAL_EFAULT );
  XLAL_CHECK( setup != NULL, XLAL_EFAULT );
  XLAL_CHECK( unmatched_item_tol >= 0, XLAL_EINVAL );
  XLAL_CHECK( param_tol_mism >= 0, XLAL_EINVAL );
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

  // Compare list of detectors
  {
    const BOOLEAN have_det_1 = ( out_1->statistics_params->detectors != NULL ), have_det_2 = ( out_2->statistics_params->detectors != NULL );
    if ( have_det_1 != have_det_2 ) {
      *equal = 0;
      XLALPrintInfo( "%s: unequal presence of detector lists: %i != %i\n", __func__, have_det_1, have_det_2 );
      return XLAL_SUCCESS;
    }
    if ( have_det_1 ) {
      if ( out_1->statistics_params->detectors->length != out_2->statistics_params->detectors->length ) {
        *equal = 0;
        XLALPrintInfo( "%s: unequal number of detectors: %u != %u\n", __func__, out_1->statistics_params->detectors->length, out_2->statistics_params->detectors->length );
        return XLAL_SUCCESS;
      }
      for ( size_t i = 0; i < out_1->statistics_params->detectors->length; ++i ) {
        if ( strcmp( out_1->statistics_params->detectors->data[i], out_2->statistics_params->detectors->data[i] ) != 0 ) {
          *equal = 0;
          XLALPrintInfo( "%s: unequal detectors: %s != %s\n", __func__, out_1->statistics_params->detectors->data[i], out_2->statistics_params->detectors->data[i] );
          return XLAL_SUCCESS;
        }
      }
    }
  }

  // Compare number of per-segment items
  if ( out_1->statistics_params->nsegments != out_2->statistics_params->nsegments ) {
    *equal = 0;
    XLALPrintInfo( "%s: unequal number of segments: %u != %u\n", __func__, out_1->statistics_params->nsegments, out_2->statistics_params->nsegments );
    return XLAL_SUCCESS;
  }

  // Compare toplist statistics
  if ( out_1->statistics_params->toplist_statistics != out_2->statistics_params->toplist_statistics ) {
    *equal = 0;
    char *toplists1, *toplists2;
    toplists1 = XLALPrintStringValueOfUserFlag( ( const int * )&( out_1->statistics_params->toplist_statistics ), &WeaveToplistChoices );
    XLAL_CHECK( toplists1 != NULL, XLAL_EFUNC );
    toplists2 = XLALPrintStringValueOfUserFlag( ( const int * )&( out_2->statistics_params->toplist_statistics ), &WeaveToplistChoices );
    XLAL_CHECK( toplists2 != NULL, XLAL_EFUNC );
    XLALPrintError( "%s: Inconsistent set of toplist statistics: {%s} != {%s}\n", __func__, toplists1, toplists2 );
    XLALFree( toplists1 );
    XLALFree( toplists2 );
    return XLAL_SUCCESS;
  }

  // Compare statistics_to_output
  for ( UINT4 istage = 0; istage < 2; ++ istage ) {
    if ( out_1->statistics_params->statistics_to_output[istage] != out_2->statistics_params->statistics_to_output[istage] ) {
      *equal = 0;
      char *outputs1, *outputs2;
      outputs1 = XLALPrintStringValueOfUserFlag( ( const int * )&( out_1->statistics_params->statistics_to_output[istage] ), &WeaveStatisticChoices );
      XLAL_CHECK( outputs1 != NULL, XLAL_EFUNC );
      outputs2 = XLALPrintStringValueOfUserFlag( ( const int * )&( out_2->statistics_params->statistics_to_output[istage] ), &WeaveStatisticChoices );
      XLAL_CHECK( outputs2 != NULL, XLAL_EFUNC );
      XLALPrintError( "%s: Inconsistent set of stage-%d ouput statistics: {%s} != {%s}\n", __func__, istage, outputs1, outputs2 );
      XLALFree( outputs1 );
      XLALFree( outputs2 );
      return XLAL_SUCCESS;
    }
  }

  // Compare toplists
  for ( size_t i = 0; i < out_1->ntoplists; ++i ) {
    XLAL_CHECK( XLALWeaveResultsToplistCompare( equal,
                                                setup, sort_by_semi_phys,
                                                round_param_to_dp, round_param_to_sf, unmatched_item_tol, param_tol_mism, result_tol, toplist_compare_limit,
                                                out_1->toplists[i], out_2->toplists[i]
                  ) == XLAL_SUCCESS, XLAL_EFUNC );
    if ( !*equal ) {
      return XLAL_SUCCESS;
    }
  }


  return XLAL_SUCCESS;

}

// Local Variables:
// c-file-style: "linux"
// c-basic-offset: 2
// End:
