/*
 * Copyright (C) 2020 David Keitel
 * Copyright (C) 2014, 2015 Evan Goetz
 * Copyright (C) 2014 Matthew Pitkin
 * Copyright (C) 2010, 2012--2017, 2020, 2022 Karl Wette
 * Copyright (C) 2010 Chris Messenger
 * Copyright (C) 2009, 2011, 2015 Adam Mercer
 * Copyright (C) 2004--2006, 2009-- 2011, 2013, 2015--2018 Reinhard Prix
 * Copyright (C) 2004--2006 Bernd Machenschalk
 * Copyright (C) 2004, 2005 Alicia Sintes
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with with program; see the file COPYING. If not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA  02110-1301  USA
 */

/*---------- includes ----------*/

#include <lal/ConfigFile.h>

#include "SFTinternal.h"

/*---------- global variables ----------*/

// XLALReadSegmentsFromFile(): applications which still must support
// the deprecated 4-column format should set this variable to non-zero
int XLALReadSegmentsFromFile_support_4column_format = 0;

/*========== function definitions ==========*/

/// \addtogroup SFTfileIO_h
/// @{

/** Allocate a LIGOTimeGPSVector */
LIGOTimeGPSVector *
XLALCreateTimestampVector( UINT4 length )
{
  int len;
  LIGOTimeGPSVector *out = XLALCalloc( 1, len = sizeof( LIGOTimeGPSVector ) );
  if ( out == NULL ) {
    XLAL_ERROR_NULL( XLAL_ENOMEM, "Failed to allocate XLALCalloc(1,%d)\n", len );
  }

  out->length = length;
  out->data = XLALCalloc( 1, len = length * sizeof( LIGOTimeGPS ) );
  if ( out->data == NULL ) {
    XLALFree( out );
    XLAL_ERROR_NULL( XLAL_ENOMEM, "Failed to allocate XLALCalloc(1,%d)\n", len );
  }

  return out;

} /* XLALCreateTimestampVector() */


/** De-allocate a LIGOTimeGPSVector */
void
XLALDestroyTimestampVector( LIGOTimeGPSVector *vect )
{
  if ( !vect ) {
    return;
  }

  XLALFree( vect->data );
  XLALFree( vect );

  return;

} /* XLALDestroyTimestampVector() */


/** Resize a LIGOTimeGPSVector */
LIGOTimeGPSVector *
XLALResizeTimestampVector( LIGOTimeGPSVector *vector, UINT4 length )
{
  if ( ! vector ) {
    return XLALCreateTimestampVector( length );
  }
  if ( ! length ) {
    XLALDestroyTimestampVector( vector );
    return NULL;
  }

  vector->data = XLALRealloc( vector->data, length * sizeof( LIGOTimeGPS ) );

  if ( ! vector->data ) {
    vector->length = 0;
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }
  vector->length = length;
  return vector;
}


/**
 * Simple creator function for MultiLIGOTimeGPSVector with numDetectors entries
 */
MultiLIGOTimeGPSVector *
XLALCreateMultiLIGOTimeGPSVector( UINT4 numDetectors )
{
  MultiLIGOTimeGPSVector *ret;

  if ( ( ret = XLALMalloc( sizeof( *ret ) ) ) == NULL ) {
    XLALPrintError( "%s: XLALMalloc(%zu) failed.\n", __func__, sizeof( *ret ) );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  ret->length = numDetectors;
  if ( ( ret->data = XLALCalloc( numDetectors, sizeof( *ret->data ) ) ) == NULL ) {
    XLALPrintError( "%s: XLALCalloc(%d, %zu) failed.\n", __func__, numDetectors, sizeof( *ret->data ) );
    XLALFree( ret );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  return ret;

} /* XLALCreateMultiLIGOTimeGPSVector() */


/**
 * Destroy a MultiLIGOTimeGPSVector timestamps vector
 */
void
XLALDestroyMultiTimestamps( MultiLIGOTimeGPSVector *multiTS )
{
  UINT4 numIFOs, X;

  if ( !multiTS ) {
    return;
  }

  numIFOs = multiTS->length;
  for ( X = 0; X < numIFOs; X ++ ) {
    XLALDestroyTimestampVector( multiTS->data[X] );
  }

  XLALFree( multiTS->data );
  XLALFree( multiTS );

  return;

} /* XLALDestroyMultiTimestamps() */


// Find index values of first and last timestamp within given timeslice range XLALCWGPSinRange(minStartGPS, maxStartGPS)
int
XLALFindTimesliceBounds( UINT4 *iStart,
                         UINT4 *iEnd,
                         const LIGOTimeGPSVector *timestamps,
                         const LIGOTimeGPS *minStartGPS,
                         const LIGOTimeGPS *maxStartGPS
                       )
{
  XLAL_CHECK( ( iStart != NULL ) && ( iEnd != NULL ) && ( timestamps != NULL ) && ( minStartGPS != NULL ) && ( maxStartGPS != NULL ), XLAL_EINVAL );
  XLAL_CHECK( XLALGPSCmp( minStartGPS, maxStartGPS ) < 1, XLAL_EINVAL, "minStartGPS (%"LAL_GPS_FORMAT") is greater than maxStartGPS (%"LAL_GPS_FORMAT")\n",
              LAL_GPS_PRINT( *minStartGPS ), LAL_GPS_PRINT( *maxStartGPS ) );

  UINT4 N = timestamps->length;
  ( *iStart ) = 0;
  ( *iEnd )   = N - 1;

  // check if there's any timestamps falling into the requested timeslice at all
  if ( ( ( XLALCWGPSinRange( timestamps->data[0], minStartGPS, maxStartGPS ) == 1 ) || ( XLALCWGPSinRange( timestamps->data[N - 1], minStartGPS, maxStartGPS ) == -1 ) ) ) {
    // if not: set an emtpy index interval in this case
    ( *iStart ) = 1;
    ( *iEnd )   = 0;
    XLALPrintInfo( "Returning empty timeslice: Timestamps span [%"LAL_GPS_FORMAT", %"LAL_GPS_FORMAT "]"
                   " has no overlap with requested timeslice range [%"LAL_GPS_FORMAT", %"LAL_GPS_FORMAT").\n",
                   LAL_GPS_PRINT( timestamps->data[0] ), LAL_GPS_PRINT( timestamps->data[N - 1] ),
                   LAL_GPS_PRINT( *minStartGPS ), LAL_GPS_PRINT( *maxStartGPS )
                 );
    return XLAL_SUCCESS;
  }

  while ( ( *iStart ) <= ( *iEnd ) && XLALCWGPSinRange( timestamps->data[( *iStart ) ], minStartGPS, maxStartGPS ) < 0 ) {
    ++ ( *iStart );
  }
  while ( ( *iStart ) <= ( *iEnd ) && XLALCWGPSinRange( timestamps->data[( *iEnd ) ], minStartGPS, maxStartGPS ) > 0 ) {
    -- ( *iEnd );
  }
  // note: *iStart >=0, *iEnd >= 0 is now guaranteed due to previous range overlap-check

  // check if there is any timestamps found witin the interval, ie if iStart <= iEnd
  if ( ( *iStart ) > ( *iEnd ) ) {
    XLALPrintInfo( "Returning empty timeslice: no sfttimes fall within given GPS range [%"LAL_GPS_FORMAT", %"LAL_GPS_FORMAT"). "
                   "Closest timestamps are: %"LAL_GPS_FORMAT" and %"LAL_GPS_FORMAT"\n",
                   LAL_GPS_PRINT( *minStartGPS ), LAL_GPS_PRINT( *maxStartGPS ),
                   LAL_GPS_PRINT( timestamps->data[( *iEnd ) ] ), LAL_GPS_PRINT( timestamps->data[( *iStart ) ] )
                 );
    ( *iStart ) = 1;
    ( *iEnd )   = 0;
  }

  return XLAL_SUCCESS;

} // XLALFindTimesliceBounds()


/**
 * Given a start-time, Tspan, Tsft and Toverlap, returns a list of timestamps
 * covering this time-stretch (allowing for overlapping SFTs).
 *
 * NOTE: boundary-handling: the returned list of timestamps fall within the
 * interval [tStart, tStart+Tspan),
 * consistent with the convention defined in XLALCWGPSinRange().
 * Assuming each SFT covers a stretch of data of length 'Tsft',
 * the returned timestamps correspond to a set of SFTs that is guaranteed
 * to *cover* all data between 'tStart' and 'tStart+Tspan'.
 * This implies that, while the last timestamp returned will always be
 * 'ret->data[numSFTs-1] < tStart+Tspan',
 * the actual data coverage can extend up to 'Tsft' beyond 'tStart+duration'.
 */
LIGOTimeGPSVector *
XLALMakeTimestamps( LIGOTimeGPS tStart,        /**< GPS start-time */
                    REAL8 Tspan,               /**< total duration to cover, in seconds */
                    REAL8 Tsft,                /**< length of the SFT corresponding to each timestamp, in seconds */
                    REAL8 Toverlap             /**< time to overlap successive SFTs by, in seconds */
                  )
{
  XLAL_CHECK_NULL( Tspan > 0, XLAL_EDOM );
  XLAL_CHECK_NULL( Tsft  > 0, XLAL_EDOM );
  XLAL_CHECK_NULL( Toverlap  >= 0, XLAL_EDOM );
  XLAL_CHECK_NULL( Toverlap < Tsft, XLAL_EDOM );        // we must actually advance

  REAL8 Tstep = Tsft - Toverlap;        // guaranteed > 0
  UINT4 numSFTsMax = ceil( Tspan * fudge_down / Tstep );                        /* >= 1 !*/
  // now we might be covering the end-time several times, if using overlapping SFTs, so
  // let's trim this back down so that end-time is covered exactly once
  UINT4 numSFTs = numSFTsMax;
  while ( ( numSFTs >= 2 ) && ( ( numSFTs - 2 ) * Tstep + Tsft >= Tspan ) ) {
    numSFTs --;
  }

  LIGOTimeGPSVector *ret;
  XLAL_CHECK_NULL( ( ret = XLALCreateTimestampVector( numSFTs ) ) != NULL, XLAL_EFUNC );

  ret->deltaT = Tsft;

  LIGOTimeGPS tt = tStart;      /* initialize to start-time */
  for ( UINT4 i = 0; i < numSFTs; i++ ) {
    ret->data[i] = tt;
    /* get next time-stamp */
    /* NOTE: we add the interval tStep successively (rounded correctly to ns each time!)
     * instead of using iSFT*Tsft, in order to avoid possible ns-rounding problems
     * with REAL8 intervals, which becomes critial from about 100days on...
     */
    XLAL_CHECK_NULL( XLALGPSAdd( &tt, Tstep ) != NULL, XLAL_EFUNC );

  } /* for i < numSFTs */

  return ret;

} /* XLALMakeTimestamps() */


/**
 * Same as XLALMakeTimestamps() just for several detectors,
 * additionally specify the number of detectors.
 */
MultiLIGOTimeGPSVector *
XLALMakeMultiTimestamps( LIGOTimeGPS tStart,   /**< GPS start-time */
                         REAL8 Tspan,          /**< total duration to cover, in seconds */
                         REAL8 Tsft,           /**< Tsft: SFT length of each timestamp, in seconds */
                         REAL8 Toverlap,       /**< time to overlap successive SFTs by, in seconds */
                         UINT4 numDet          /**< number of timestamps-vectors to generate */
                       )
{
  XLAL_CHECK_NULL( numDet >= 1, XLAL_EINVAL );

  MultiLIGOTimeGPSVector *ret;
  XLAL_CHECK_NULL( ( ret = XLALCalloc( 1, sizeof( *ret ) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL( ( ret->data = XLALCalloc( numDet, sizeof( ret->data[0] ) ) ) != NULL, XLAL_ENOMEM );
  ret->length = numDet;

  for ( UINT4 X = 0; X < numDet; X ++ ) {
    XLAL_CHECK_NULL( ( ret->data[X] = XLALMakeTimestamps( tStart, Tspan, Tsft, Toverlap ) ) != NULL, XLAL_EFUNC );
  } // for X < numDet

  return ret;

} /* XLALMakeMultiTimestamps() */


/// backwards compatible wrapper to XLALReadTimestampsFileConstrained() without GPS-time constraints
LIGOTimeGPSVector *
XLALReadTimestampsFile( const CHAR *fname )
{
  return XLALReadTimestampsFileConstrained( fname, NULL, NULL );
}


/// backwards compatible wrapper to XLALReadMultiTimestampsFilesConstrained() without GPS-time constraints
MultiLIGOTimeGPSVector *
XLALReadMultiTimestampsFiles( const LALStringVector *fnames )
{
  return XLALReadMultiTimestampsFilesConstrained( fnames, NULL, NULL );
}


/**
 * Load timestamps file 'fname' into LIGOTimeGPSVector struct, allocated here.
 *
 * The timestamps file must contain one GPS time per line, allowing for '%#' as comments, which are ignored.
 * The constraints 'minGPS', 'maxGPS' are applied by returning only timestamps that fall within
 * the range defined by XLALCWGPSinRange(gps, minGPS, maxGPS) == 0.
 */
LIGOTimeGPSVector *
XLALReadTimestampsFileConstrained( const CHAR *fname, const LIGOTimeGPS *minGPS, const LIGOTimeGPS *maxGPS )
{
  /** check input consistency */
  XLAL_CHECK_NULL( fname != NULL, XLAL_EINVAL );

  /* read and parse timestamps-list file contents*/
  LALParsedDataFile *flines = NULL;
  XLAL_CHECK_NULL( XLALParseDataFile( &flines, fname ) == XLAL_SUCCESS, XLAL_EFUNC );

  UINT4 numTS = flines->lines->nTokens;
  UINT4 numTSinRange = 0;

  /* allocate and initialized segment list */
  LIGOTimeGPSVector *timestamps = NULL;
  XLAL_CHECK_NULL( ( timestamps = XLALCreateTimestampVector( numTS ) ) != NULL, XLAL_EFUNC );

  char buf[256];
  for ( UINT4 iTS = 0; iTS < numTS; iTS ++ ) {
    LIGOTimeGPS gps;
    INT4 secs, ns;
    char junk[11] = "";

    // first check if obsolete <sec ns> format is found on this line
    if ( sscanf( flines->lines->tokens[iTS], "%" LAL_INT4_FORMAT " %" LAL_INT4_FORMAT "%10s\n", &secs, &ns, junk ) > 1 ) {
      gps.gpsSeconds = secs;
      gps.gpsNanoSeconds = ns;

      XLALPrintWarning( "Line %d: found obsolete 'sec ns' timestamps format '%s', use 'xx.yy[GPS|MJD]' instead: %s\n", iTS, flines->lines->tokens[iTS], XLALGPSToStr( buf, &gps ) );
      if ( junk[0] != 0 ) {
        XLALDestroyTimestampVector( timestamps );
        XLALDestroyParsedDataFile( flines );
        XLAL_ERROR_NULL( XLAL_EINVAL, "Unconverted trailing junk '%s' found: invalid\n", junk );
      }
    } // end: if old-style format 'sec ns' is found
    else {
      if ( XLALParseStringValueAsEPOCH( &gps, flines->lines->tokens[iTS] ) != XLAL_SUCCESS ) {
        XLALDestroyTimestampVector( timestamps );
        XLALDestroyParsedDataFile( flines );
        XLAL_ERROR_NULL( XLAL_EINVAL, "Failed to parse line %d into epoch: '%s'\n", iTS, flines->lines->tokens[iTS] );
      }
    }

    if ( XLALCWGPSinRange( gps, minGPS, maxGPS ) == 0 ) {
      timestamps->data[numTSinRange ++] = gps;
    }

  } /* for iTS < numTS */

  /* free parsed segment file contents */
  XLALDestroyParsedDataFile( flines );

  // adjust size of timestamps vector to those found within range
  timestamps->length = numTSinRange;
  timestamps->data = XLALRealloc( timestamps->data, numTSinRange * sizeof( timestamps->data[0] ) );

  return timestamps;

} /* XLALReadTimestampsFileConstrained() */


/**
 * Load several timestamps files, return a MultiLIGOTimeGPSVector struct, allocated here.
 *
 * The timestamps files must contain one GPS time per line, allowing for '%#' as comments, which are ignored.
 * The constraints 'minGPS', 'maxGPS' are applied by returning only timestamps that fall within
 * the range defined by XLALCWGPSinRange(gps, minGPS, maxGPS) == 0.
 *
 */
MultiLIGOTimeGPSVector *
XLALReadMultiTimestampsFilesConstrained( const LALStringVector *fnames, const LIGOTimeGPS *minGPS, const LIGOTimeGPS *maxGPS )
{
  XLAL_CHECK_NULL( fnames != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL( fnames->data != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL( fnames->length > 0, XLAL_EDOM );

  UINT4 numDet = fnames->length;

  // ----- prepare output container
  MultiLIGOTimeGPSVector *multiTS;
  XLAL_CHECK_NULL( ( multiTS = XLALCalloc( 1, sizeof( *multiTS ) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL( ( multiTS->data = XLALCalloc( numDet, sizeof( multiTS->data[0] ) ) ) != NULL, XLAL_ENOMEM );
  multiTS->length = numDet;

  for ( UINT4 X = 0; X < numDet; X ++ ) {
    XLAL_CHECK_NULL( fnames->data[X] != NULL, XLAL_EINVAL );
    XLAL_CHECK_NULL( ( multiTS->data[X] = XLALReadTimestampsFileConstrained( fnames->data[X], minGPS, maxGPS ) ) != NULL, XLAL_EFUNC );
  } // for X < numDet

  return multiTS;

} // XLALReadMultiTimestampsFilesConstrained()


/**
 * Extract timstamps-vector from the given SFTVector
 */
LIGOTimeGPSVector *
XLALExtractTimestampsFromSFTs( const SFTVector *sfts )          /**< [in] input SFT-vector  */
{
  /* check input consistency */
  if ( !sfts ) {
    XLALPrintError( "%s: invalid NULL input 'sfts'\n", __func__ );
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }

  UINT4 numSFTs = sfts->length;
  /* create output vector */
  LIGOTimeGPSVector *ret = NULL;
  if ( ( ret = XLALCreateTimestampVector( numSFTs ) ) == NULL ) {
    XLALPrintError( "%s: XLALCreateTimestampVector(%d) failed.\n", __func__, numSFTs );
    XLAL_ERROR_NULL( XLAL_EFUNC );
  }
  ret->deltaT = TSFTfromDFreq( sfts->data[0].deltaF );

  UINT4 i;
  for ( i = 0; i < numSFTs; i ++ ) {
    ret->data[i] = sfts->data[i].epoch;
  }

  /* done: return Ts-vector */
  return ret;

} /* XLALExtractTimestampsFromSFTs() */


/**
 * Given a multi-SFT vector, return a MultiLIGOTimeGPSVector holding the
 * SFT timestamps
 */
MultiLIGOTimeGPSVector *
XLALExtractMultiTimestampsFromSFTs( const MultiSFTVector *multiSFTs )
{
  /* check input consistency */
  if ( !multiSFTs || multiSFTs->length == 0 ) {
    XLALPrintError( "%s: illegal NULL or empty input 'multiSFTs'.\n", __func__ );
    XLAL_ERROR_NULL( XLAL_EINVAL );
  }
  UINT4 numIFOs = multiSFTs->length;

  /* create output vector */
  MultiLIGOTimeGPSVector *ret = NULL;
  if ( ( ret = XLALCalloc( 1, sizeof( *ret ) ) ) == NULL ) {
    XLALPrintError( "%s: failed to XLALCalloc ( 1, %zu ).\n", __func__, sizeof( *ret ) );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  if ( ( ret->data = XLALCalloc( numIFOs, sizeof( *ret->data ) ) ) == NULL ) {
    XLALPrintError( "%s: failed to XLALCalloc ( %d, %zu ).\n", __func__, numIFOs, sizeof( ret->data[0] ) );
    XLALFree( ret );
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }
  ret->length = numIFOs;

  /* now extract timestamps vector from each SFT-vector */
  UINT4 X;
  for ( X = 0; X < numIFOs; X ++ ) {
    if ( ( ret->data[X] = XLALExtractTimestampsFromSFTs( multiSFTs->data[X] ) ) == NULL ) {
      XLALPrintError( "%s: XLALExtractTimestampsFromSFTs() failed for X=%d\n", __func__, X );
      XLALDestroyMultiTimestamps( ret );
      XLAL_ERROR_NULL( XLAL_EFUNC );
    }

  } /* for X < numIFOs */

  return ret;

} /* XLALExtractMultiTimestampsFromSFTs() */


/**
 * Extract timestamps-vector of *unique* timestamps from the given SFTCatalog
 *
 * NOTE: when dealing with catalogs of frequency-slided SFTs, each timestamp will appear in the
 * catalog multiple times, depending on how many frequency slices have been read in.
 * In such cases this function will return the list of *unique* timestamps.
 *
 * NOTE 2: This function will also enfore the multiplicity of each timestamp to be the
 * same through the whole catalog, corresponding to the case of 'frequency-sliced' SFTs,
 * while non-constant multiplicities would indicate a potential problem somewhere.
 */
LIGOTimeGPSVector *
XLALTimestampsFromSFTCatalog( const SFTCatalog *catalog )               /**< [in] input SFT-catalog  */
{
  // check input consistency
  XLAL_CHECK_NULL( catalog != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL( catalog->length > 0, XLAL_EINVAL );

  UINT4 numEntries = catalog->length;

  // create output vector, assuming maximal length, realloc at the end
  LIGOTimeGPSVector *ret;
  XLAL_CHECK_NULL( ( ret = XLALCreateTimestampVector( numEntries ) ) != NULL, XLAL_EFUNC );

  REAL8 Tsft0 = 1.0 / catalog->data[0].header.deltaF;
  if ( fabs( ( Tsft0 - round( Tsft0 ) ) ) / Tsft0 < 10 * LAL_REAL8_EPS ) { // 10-eps 'snap' to closest integer
    ret->deltaT = round( Tsft0 );
  } else {
    ret->deltaT = Tsft0;
  }

  // For dealing with SFTCatalogs corresponding to frequency-sliced input SFTs:
  // Given the guaranteed GPS-ordering of XLALSFTDataFind(), we can rely on duplicate
  // timestamps to all be found next to each other, and therefore can easily skip them
  ret->data[0] = catalog->data[0].header.epoch;
  UINT4 numUnique = 1;
  UINT4 stride = 0;
  for ( UINT4 i = 1; i < numEntries; i ++ ) {
    UINT4 thisStride = 1;
    const LIGOTimeGPS *ti   = &( catalog->data[i].header.epoch );
    const LIGOTimeGPS *tim1 = &( catalog->data[i - 1].header.epoch );
    if ( XLALGPSCmp( ti, tim1 ) == 0 ) {
      thisStride ++;
      continue;       // skip duplicates
    }
    ret->data[numUnique] = catalog->data[i].header.epoch;
    numUnique ++;

    // keep track of stride, ensure that it's the same for every unique timestamp
    if ( stride == 0 ) {
      stride = thisStride;
    } else {
      XLAL_CHECK_NULL( stride == thisStride, XLAL_EINVAL, "Suspicious SFT Catalog with non-constant timestamps multiplicities '%u != %u'\n", stride, thisStride );
    }
  } // for i < numEntries


  // now truncate output vector to actual length of unique timestamps
  ret->length = numUnique;
  XLAL_CHECK_NULL( ( ret->data = XLALRealloc( ret->data, numUnique * sizeof( ( *ret->data ) ) ) ) != NULL, XLAL_ENOMEM );

  // done: return Ts-vector
  return ret;

} /* XLALTimestampsFromSFTCatalog() */


/**
 * Given a multi-SFTCatalogView, return a MultiLIGOTimeGPSVector holding the
 * SFT timestamps
 */
MultiLIGOTimeGPSVector *
XLALTimestampsFromMultiSFTCatalogView( const MultiSFTCatalogView *multiView )
{
  /* check input consistency */
  XLAL_CHECK_NULL( multiView != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL( multiView->length > 0, XLAL_EINVAL );

  UINT4 numIFOs = multiView->length;

  /* create output vector */
  MultiLIGOTimeGPSVector *ret;
  XLAL_CHECK_NULL( ( ret = XLALCalloc( 1, sizeof( *ret ) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL( ( ret->data = XLALCalloc( numIFOs, sizeof( *( ret->data ) ) ) ) != NULL, XLAL_ENOMEM );
  ret->length = numIFOs;

  /* now extract timestamps vector from each IFO's SFT-Catalog */
  for ( UINT4 X = 0; X < numIFOs; X ++ ) {
    XLAL_CHECK_NULL( ( ret->data[X] = XLALTimestampsFromSFTCatalog( &( multiView->data[X] ) ) ) != NULL, XLAL_EFUNC );
  } /* for X < numIFOs */

  return ret;

} /* XLALTimestampsFromMultiSFTCatalogView() */


/**
 * Function to read a segment list from given filename, returns a *sorted* LALSegList
 *
 * The segment list file format is repeated lines (excluding comment lines beginning with
 * <tt>\%</tt> or <tt>#</tt>) of one of the following forms:
 * - <tt>startGPS endGPS</tt>
 * - <tt>startGPS endGPS NumSFTs</tt> (NumSFTs must be a positive integer)
 * - <tt>startGPS endGPS duration NumSFTs</tt> (\b DEPRECATED, duration is ignored)
 *
 * \note We (ab)use the integer \p id field in LALSeg to carry the total number of SFTs
 * contained in that segment if <tt>NumSFTs</tt> was provided in the segment file.
 * This can be used as a consistency check when loading SFTs for these segments.
 */
LALSegList *
XLALReadSegmentsFromFile( const char *fname     /**< name of file containing segment list */
                        )
{
  LALSegList *segList = NULL;

  /* check input consistency */
  XLAL_CHECK_NULL( fname != NULL, XLAL_EFAULT );

  /* read and parse segment-list file contents*/
  LALParsedDataFile *flines = NULL;
  XLAL_CHECK_NULL( XLALParseDataFile( &flines, fname ) == XLAL_SUCCESS, XLAL_EFUNC );
  const UINT4 numSegments = flines->lines->nTokens;
  XLAL_CHECK_NULL( numSegments > 0, XLAL_EINVAL, "%s: segment file '%s' does not contain any segments", __func__, fname );

  /* allocate and initialized segment list */
  XLAL_CHECK_NULL( ( segList = XLALCalloc( 1, sizeof( *segList ) ) ) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL( XLALSegListInit( segList ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* determine number of columns */
  int ncol = 0;
  {
    REAL8 col[4];
    ncol = sscanf( flines->lines->tokens[0], "%lf %lf %lf %lf", &col[0], &col[1], &col[2], &col[3] );
    switch ( ncol ) {
    case 2:
    case 3:
      break;
    case 4:
      if ( XLALReadSegmentsFromFile_support_4column_format ) {
        XLALPrintError( "\n%s: WARNING: segment file '%s' is in DEPRECATED 4-column format (startGPS endGPS duration NumSFTs, duration is ignored)\n", __func__, fname );
      } else {
        XLAL_ERROR_NULL( XLAL_EIO, "%s: segment file '%s' is in DEPRECATED 4-column format (startGPS endGPS duration NumSFTs)\n", __func__, fname );
      }
      break;
    default:
      XLAL_ERROR_NULL( XLAL_EIO, "%s: segment file '%s' contains an unknown %i-column format", __func__, fname, ncol );
    }
  }

  /* parse segment list */
  for ( UINT4 iSeg = 0; iSeg < numSegments; iSeg ++ ) {

    /* parse line of segment file, depending on determined number of columns */
    REAL8 start = 0, end = 0, duration = 0;
    INT4 NumSFTs = 0;
    int ret;
    switch ( ncol ) {
    case 2:
      ret = sscanf( flines->lines->tokens[iSeg], "%lf %lf", &start, &end );
      XLAL_CHECK_NULL( ret == 2, XLAL_EIO, "%s: number of columns in segment file '%s' is inconsistent (line 1: %i, line %u: %i)", __func__, fname, ncol, iSeg + 1, ret );
      break;
    case 3:
      ret = sscanf( flines->lines->tokens[iSeg], "%lf %lf %i", &start, &end, &NumSFTs );
      XLAL_CHECK_NULL( ret == 3, XLAL_EIO, "%s: number of columns in segment file '%s' is inconsistent (line 1: %i, line %u: %i)", __func__, fname, ncol, iSeg + 1, ret );
      XLAL_CHECK_NULL( NumSFTs > 0, XLAL_EIO, "%s: number of SFTs (3rd column) in segment file '%s' must be a positive integer if given (line %u: %i)", __func__, fname, iSeg + 1, NumSFTs );
      break;
    case 4:
      ret = sscanf( flines->lines->tokens[iSeg], "%lf %lf %lf %i", &start, &end, &duration, &NumSFTs );
      XLAL_CHECK_NULL( ret == 4, XLAL_EIO, "%s: number of columns in segment file '%s' is inconsistent (line 1 = %i, line %u = %i)", __func__, fname, ncol, iSeg + 1, ret );
      break;
    default:
      XLAL_ERROR_NULL( XLAL_EFAILED, "Unexpected error!" );
    }

    /* set GPS start and end times */
    LIGOTimeGPS startGPS, endGPS;
    XLALGPSSetREAL8( &startGPS, start );
    XLALGPSSetREAL8( &endGPS, end );

    /* create segment and append to list
       - we set number of SFTs as 'id' field, as we have no other use for it */
    LALSeg thisSeg;
    XLAL_CHECK_NULL( XLALSegSet( &thisSeg, &startGPS, &endGPS, NumSFTs ) == XLAL_SUCCESS, XLAL_EFUNC );
    XLAL_CHECK_NULL( XLALSegListAppend( segList, &thisSeg ) == XLAL_SUCCESS, XLAL_EFUNC );

  } /* for iSeg < numSegments */

  /* sort final segment list in increasing GPS start-times */
  XLAL_CHECK_NULL( XLALSegListSort( segList ) == XLAL_SUCCESS, XLAL_EFUNC );

  /* free parsed segment file contents */
  XLALDestroyParsedDataFile( flines );

  return segList;

} /* XLALReadSegmentsFromFile() */


/**
 * Extract timestamps-vector from a segment file, with functionality based on MakeSFTDAG
 * The filename should point to a file containing \<GPSstart GPSend\> of segments or \<GPSstart GPSend segLength numSFTs\> where segLength is in hours.
 * adjustSegExtraTime is used in MakeSFTDAG to maximize the number of SFTs in each segement by placing the SFTs in the middle of the segment.
 * synchronize is used to force the start times of the SFTs to be integer multiples of Toverlap from the start time of the first SFT.
 * adjustSegExtraTime and synchronize cannot be used concurrently (synchronize will be preferred if both values are non-zero).
 */
LIGOTimeGPSVector *
XLALTimestampsFromSegmentFile( const char *filename,       //!< filename: Input filename
                               REAL8 Tsft,                 //!< Tsft: SFT length of each timestamp, in seconds
                               REAL8 Toverlap,             //!< Toverlap: time to overlap successive SFTs by, in seconds
                               BOOLEAN adjustSegExtraTime, //!< adjustSegExtraTime: remove the unused time from beginning and end of the segments (see MakeSFTDAG)
                               BOOLEAN synchronize         //!< synchronize: synchronize SFT start times according to the start time of the first SFT. Start time of first SFT is shifted to next higher integer value of Tsft
                             )
{
  XLAL_CHECK_NULL( filename != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL( !( adjustSegExtraTime && synchronize ), XLAL_EINVAL, "Must specify only one of adjustSegExtraTime or synchronize" );

  LALSegList *list = NULL;
  XLAL_CHECK_NULL( ( list = XLALReadSegmentsFromFile( filename ) ) != NULL, XLAL_EFUNC );

  //Need to know the number of SFTs before creating the timestamps vector, so we have to do the same loop twice
  INT4 numSFTs = 0;
  REAL8 firstSFTstartTime = 0.0, overlapFraction = Toverlap / Tsft;

  for ( UINT4 i = 0; i < list->length; i++ ) {
    INT4 numThisSeg = 0;
    REAL8 analysisStartTime, analysisEndTime;
    if ( adjustSegExtraTime && !synchronize ) {
      REAL8 segStartTime = XLALGPSGetREAL8( &( list->segs[i].start ) );
      REAL8 segEndTime = XLALGPSGetREAL8( &( list->segs[i].end ) );
      REAL8 segExtraTime = fmod( ( segEndTime - segStartTime ), Tsft );
      if ( overlapFraction != 0.0 ) {
        if ( ( segEndTime - segStartTime ) > Tsft ) {
          segExtraTime = fmod( ( segEndTime - segStartTime - Tsft ), ( ( 1.0 - overlapFraction ) * Tsft ) );
        }
      } else {
        segExtraTime = fmod( ( segEndTime - segStartTime ), Tsft );
      }
      REAL8 segExtraStart =  segExtraTime / 2;
      REAL8 segExtraEnd = segExtraTime - segExtraStart;
      analysisStartTime = segStartTime + segExtraStart;
      if ( analysisStartTime > segEndTime ) {
        analysisStartTime = segEndTime;
      }
      analysisEndTime = segEndTime - segExtraEnd;
      if ( analysisEndTime < segStartTime ) {
        analysisEndTime = segStartTime;
      }
    } // if adjustSegExtraTime && !synchronize
    else if ( synchronize ) {
      REAL8 segStartTime = XLALGPSGetREAL8( &( list->segs[i].start ) );
      REAL8 segEndTime = XLALGPSGetREAL8( &( list->segs[i].end ) );
      if ( firstSFTstartTime == 0.0 ) {
        firstSFTstartTime = ceil( segStartTime / Tsft ) * Tsft;
      }
      analysisStartTime = round( ceil( ( segStartTime - firstSFTstartTime ) / ( ( 1.0 - overlapFraction ) * Tsft ) ) * ( 1.0 - overlapFraction ) * Tsft ) + firstSFTstartTime;
      if ( analysisStartTime > segEndTime ) {
        analysisStartTime = segEndTime;
      }
      analysisEndTime = round( floor( ( segEndTime - analysisStartTime - Tsft ) / ( ( 1.0 - overlapFraction ) * Tsft ) ) * ( 1.0 - overlapFraction ) * Tsft ) + Tsft + analysisStartTime;
      if ( analysisEndTime < segStartTime ) {
        analysisEndTime = segStartTime;
      }
    } else {
      analysisStartTime = XLALGPSGetREAL8( &( list->segs[i].start ) );
      analysisEndTime = XLALGPSGetREAL8( &( list->segs[i].end ) );
    }

    REAL8 endTime = analysisStartTime;
    while ( endTime < analysisEndTime ) {
      if ( numThisSeg == 0 ) {
        endTime += Tsft;
      } else {
        endTime += ( 1.0 - overlapFraction ) * Tsft;
      }
      if ( endTime <= analysisEndTime ) {
        numThisSeg++;
      }
    } // while endTime < analysisEndTime
    numSFTs += numThisSeg;
  } // for i < length

  LIGOTimeGPSVector *ret;
  XLAL_CHECK_NULL( ( ret = XLALCreateTimestampVector( numSFTs ) ) != NULL, XLAL_EFUNC );

  ret->deltaT = Tsft;

  //Second time doing the same thing, but now we can set the times of the SFTs in the timestamps vector
  firstSFTstartTime = 0.0;
  UINT4 j = 0;
  for ( UINT4 i = 0; i < list->length; i++ ) {
    INT4 numThisSeg = 0;
    REAL8 analysisStartTime, analysisEndTime;
    if ( adjustSegExtraTime && !synchronize ) {
      REAL8 segStartTime = XLALGPSGetREAL8( &( list->segs[i].start ) );
      REAL8 segEndTime = XLALGPSGetREAL8( &( list->segs[i].end ) );
      REAL8 segExtraTime = fmod( ( segEndTime - segStartTime ), Tsft );
      if ( overlapFraction != 0.0 ) {
        if ( ( segEndTime - segStartTime ) > Tsft ) {
          segExtraTime = fmod( ( segEndTime - segStartTime - Tsft ), ( ( 1.0 - overlapFraction ) * Tsft ) );
        }
      } else {
        segExtraTime = fmod( ( segEndTime - segStartTime ), Tsft );
      }
      REAL8 segExtraStart =  segExtraTime / 2;
      REAL8 segExtraEnd = segExtraTime - segExtraStart;
      analysisStartTime = segStartTime + segExtraStart;
      if ( analysisStartTime > segEndTime ) {
        analysisStartTime = segEndTime;
      }
      analysisEndTime = segEndTime - segExtraEnd;
      if ( analysisEndTime < segStartTime ) {
        analysisEndTime = segStartTime;
      }
    } // if adjustSegExtraTime && !synchronize
    else if ( synchronize ) {
      REAL8 segStartTime = XLALGPSGetREAL8( &( list->segs[i].start ) );
      REAL8 segEndTime = XLALGPSGetREAL8( &( list->segs[i].end ) );
      if ( firstSFTstartTime == 0.0 ) {
        firstSFTstartTime = ceil( segStartTime / Tsft ) * Tsft;
      }
      analysisStartTime = round( ceil( ( segStartTime - firstSFTstartTime ) / ( ( 1.0 - overlapFraction ) * Tsft ) ) * ( 1.0 - overlapFraction ) * Tsft ) + firstSFTstartTime;
      if ( analysisStartTime > segEndTime ) {
        analysisStartTime = segEndTime;
      }
      analysisEndTime = round( floor( ( segEndTime - analysisStartTime - Tsft ) / ( ( 1.0 - overlapFraction ) * Tsft ) ) * ( 1.0 - overlapFraction ) * Tsft ) + Tsft + analysisStartTime;
      if ( analysisEndTime < segStartTime ) {
        analysisEndTime = segStartTime;
      }
    } else {
      analysisStartTime = XLALGPSGetREAL8( &( list->segs[i].start ) );
      analysisEndTime = XLALGPSGetREAL8( &( list->segs[i].end ) );
    }

    REAL8 endTime = analysisStartTime;
    while ( endTime < analysisEndTime ) {
      if ( numThisSeg == 0 ) {
        endTime += Tsft;
      } else {
        endTime += ( 1.0 - overlapFraction ) * Tsft;
      }
      if ( endTime <= analysisEndTime ) {
        numThisSeg++;
        LIGOTimeGPS sftStart;
        XLALGPSSetREAL8( &sftStart, endTime - Tsft );
        ret->data[j] = sftStart;
        j++;
      }
    } // while ( endTime < analysisEndTime )
  } // for i < length

  XLALSegListFree( list );

  /* done: return Ts-vector */
  return ret;

} // XLALTimestampsFromSegmentFile()

/// @}
