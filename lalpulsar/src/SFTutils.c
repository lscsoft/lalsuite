/*
 * Copyright (C) 2005 Reinhard Prix
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

/*---------- INCLUDES ----------*/
#include <stdarg.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sort_double.h>

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FrequencySeries.h>
#include <lal/NormalizeSFTRngMed.h>
#include <lal/LISAspecifics.h>
#include <lal/Date.h>
#include <lal/Units.h>

#include <lal/SFTutils.h>

/*---------- DEFINES ----------*/
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

/*----- SWITCHES -----*/

/*---------- internal types ----------*/

/*---------- Global variables ----------*/
static REAL8 fudge_up   = 1 + 10 * LAL_REAL8_EPS;	// about ~1 + 2e-15
static REAL8 fudge_down = 1 - 10 * LAL_REAL8_EPS;	// about ~1 - 2e-15

/* empty struct initializers */
const PSDVector empty_PSDVector;
const MultiPSDVector empty_MultiPSDVector;
const MultiNoiseWeights empty_MultiNoiseWeights;

/*---------- internal prototypes ----------*/

/*==================== FUNCTION DEFINITIONS ====================*/

// ---------- obsolete LAL-API was moved into external file
#include "SFTutils-LAL.c"
// ------------------------------

/**
 * XLAL function to create one SFT-struct.
 * Note: Allows for numBins == 0, in which case only the header is
 * allocated, with a NULL data pointer.
 */
SFTtype *
XLALCreateSFT ( UINT4 numBins )
{
  SFTtype *sft;

  if ( (sft = XLALCalloc (1, sizeof(*sft) )) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "XLALCalloc (1, %d) failed.\n", sizeof(*sft) );

  if ( numBins )
    {
      if ( (sft->data = XLALCreateCOMPLEX8Vector ( numBins )) == NULL )
        {
          XLALFree ( sft );
          XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALCreateCOMPLEX8Vector ( %s ) failed. xlalErrno = %d\n", numBins, xlalErrno );
        }
    }
  else
    sft->data = NULL;	/* no data, just header */

  return sft;

} /* XLALCreateSFT() */


/** Destructor for one SFT */
void
XLALDestroySFT ( SFTtype *sft )
{
  if ( !sft )
    return;

  if ( sft->data )
    XLALDestroyCOMPLEX8Vector ( sft->data );

  XLALFree ( sft );

  return;

} /* XLALDestroySFT() */


/**
 * XLAL function to create an SFTVector of \c numSFT SFTs with \c SFTlen frequency-bins
 */
SFTVector *
XLALCreateSFTVector (UINT4 numSFTs, 	/**< number of SFTs */
		     UINT4 numBins	/**< number of frequency-bins per SFT */
		     )
{
  UINT4 iSFT;
  SFTVector *vect;

  if ( (vect = XLALCalloc ( 1, sizeof(*vect) )) == NULL ) {
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  vect->length = numSFTs;
  if ( (vect->data = XLALCalloc (1, numSFTs * sizeof ( *vect->data ) )) == NULL ) {
    XLALFree (vect);
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  for ( iSFT = 0; iSFT < numSFTs; iSFT ++)
    {
      COMPLEX8Vector *data = NULL;

      /* allow SFTs with 0 bins: only header */
      if ( numBins )
	{
	  if ( (data = XLALCreateCOMPLEX8Vector ( numBins )) == NULL )
	    {
	      UINT4 j;
	      for ( j = 0; j < iSFT; j++ )
		XLALDestroyCOMPLEX8Vector ( vect->data[j].data );
	      XLALFree (vect->data);
	      XLALFree (vect);
	      XLAL_ERROR_NULL( XLAL_ENOMEM );
	    }
	}

      vect->data[iSFT].data = data;

    } /* for iSFT < numSFTs */

  return vect;

} /* XLALCreateSFTVector() */


/**
 * XLAL interface to destroy an SFTVector
 */
void
XLALDestroySFTVector ( SFTVector *vect )
{
  if ( !vect )
    return;

  for ( UINT4 i=0; i < vect->length; i++ )
    {
      SFTtype *sft = &( vect->data[i] );
      if ( sft->data )
	{
	  if ( sft->data->data )
	    XLALFree ( sft->data->data );
	  XLALFree ( sft->data );
	}
    } // for i < numSFTs

  XLALFree ( vect->data );
  XLALFree ( vect );

  return;

} /* XLALDestroySFTVector() */


/**
 * Destroy a PSD-vector
 */
void
XLALDestroyPSDVector ( PSDVector *vect )	/**< the PSD-vector to free */
{
  if ( vect == NULL )	/* nothing to be done */
    return;

  for ( UINT4 i=0; i < vect->length; i++ )
    {
      REAL8FrequencySeries *psd = &( vect->data[i] );
      if ( psd->data )
	{
	  if ( psd->data->data )
	    XLALFree ( psd->data->data );
	  XLALFree ( psd->data );
	}
    } // for i < numPSDs

  XLALFree ( vect->data );
  XLALFree ( vect );

  return;

} /* XLALDestroyPSDVector() */


/**
 * Destroy a multi SFT-vector
 */
void
XLALDestroyMultiSFTVector ( MultiSFTVector *multvect )	/**< the SFT-vector to free */
{
  if ( multvect == NULL )	/* nothing to be done */
    return;

  for ( UINT4 i = 0; i < multvect->length; i++ )
    XLALDestroySFTVector ( multvect->data[i] );

  XLALFree( multvect->data );
  XLALFree( multvect );

  return;

} /* XLALDestroyMultiSFTVector() */


/**
 * Destroy a multi PSD-vector
 */
void
XLALDestroyMultiPSDVector ( MultiPSDVector *multvect )	/**< the SFT-vector to free */
{
  if ( multvect == NULL )
    return;

  for ( UINT4 i = 0; i < multvect->length; i++ )
    XLALDestroyPSDVector ( multvect->data[i] );

  XLALFree( multvect->data );
  XLALFree( multvect );

  return;

} /* XLALDestroyMultiPSDVector() */


/** Allocate a LIGOTimeGPSVector */
LIGOTimeGPSVector *
XLALCreateTimestampVector ( UINT4 length )
{
  int len;
  LIGOTimeGPSVector *out = XLALCalloc (1, len = sizeof(LIGOTimeGPSVector));
  if (out == NULL)
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to allocate XLALCalloc(1,%d)\n", len );

  out->length = length;
  out->data = XLALCalloc (1, len = length * sizeof(LIGOTimeGPS));
  if (out->data == NULL) {
    XLALFree (out);
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to allocate XLALCalloc(1,%d)\n", len );
  }

  return out;

} /* XLALCreateTimestampVector() */


/** De-allocate a LIGOTimeGPSVector */
void
XLALDestroyTimestampVector ( LIGOTimeGPSVector *vect)
{
  if ( !vect )
    return;

  XLALFree ( vect->data );
  XLALFree ( vect );

  return;

} /* XLALDestroyTimestampVector() */


/**
 * Given a start-time, Tspan, Tsft and Toverlap, returns a list of timestamps
 * covering this time-stretch (allowing for overlapping SFT timestamps).
 *
 * NOTE: boundary-handling: the returned list of timestamps are guaranteed to *cover* the
 * interval [tStart, tStart+duration), assuming a each timestamp covers a length of 'Tsft'
 * This implies that the actual timestamps-coverage can extend up to 'Tsft' beyond 'tStart+duration'.
 */
LIGOTimeGPSVector *
XLALMakeTimestamps ( LIGOTimeGPS tStart,	/**< GPS start-time */
                     REAL8 Tspan, 		/**< total duration to cover, in seconds */
                     REAL8 Tsft,		/**< Tsft: SFT length of each timestamp, in seconds */
                     REAL8 Toverlap		/**< time to overlap successive SFTs by, in seconds */
                     )
{
  XLAL_CHECK_NULL ( Tspan > 0, XLAL_EDOM );
  XLAL_CHECK_NULL ( Tsft  > 0, XLAL_EDOM );
  XLAL_CHECK_NULL ( Toverlap  >= 0, XLAL_EDOM );
  XLAL_CHECK_NULL ( Toverlap < Tsft, XLAL_EDOM );	// we must actually advance

  REAL8 Tstep = Tsft - Toverlap;	// guaranteed > 0
  UINT4 numSFTsMax = ceil ( Tspan * fudge_down / Tstep );			/* >= 1 !*/
  // now we might be covering the end-time several times, if using overlapping SFTs, so
  // let's trim this back down so that end-time is covered exactly once
  UINT4 numSFTs = numSFTsMax;
  while ( (numSFTs >= 2) && ( (numSFTs - 1) * Tstep + Tsft > Tspan) ) {
    numSFTs --;
  }

  LIGOTimeGPSVector *ret;
  XLAL_CHECK_NULL ( (ret = XLALCreateTimestampVector ( numSFTs )) != NULL, XLAL_EFUNC );

  ret->deltaT = Tsft;

  LIGOTimeGPS tt = tStart;	/* initialize to start-time */
  for ( UINT4 i = 0; i < numSFTs; i++ )
    {
      ret->data[i] = tt;
      /* get next time-stamp */
      /* NOTE: we add the interval tStep successively (rounded correctly to ns each time!)
       * instead of using iSFT*Tsft, in order to avoid possible ns-rounding problems
       * with REAL8 intervals, which becomes critial from about 100days on...
       */
      XLAL_CHECK_NULL ( XLALGPSAdd ( &tt, Tstep ) != NULL, XLAL_EFUNC );

    } /* for i < numSFTs */

  return ret;

} /* XLALMakeTimestamps() */


/**
 * Same as XLALMakeTimestamps() just for several detectors,
 * additionally specify the number of detectors.
 */
MultiLIGOTimeGPSVector *
XLALMakeMultiTimestamps ( LIGOTimeGPS tStart,	/**< GPS start-time */
                          REAL8 Tspan, 		/**< total duration to cover, in seconds */
                          REAL8 Tsft,		/**< Tsft: SFT length of each timestamp, in seconds */
                          REAL8 Toverlap,	/**< time to overlap successive SFTs by, in seconds */
                          UINT4 numDet		/**< number of timestamps-vectors to generate */
                          )
{
  XLAL_CHECK_NULL ( numDet >= 1, XLAL_EINVAL );

  MultiLIGOTimeGPSVector *ret;
  XLAL_CHECK_NULL ( ( ret = XLALCalloc ( 1, sizeof(*ret))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( ( ret->data = XLALCalloc ( numDet, sizeof(ret->data[0]) )) != NULL, XLAL_ENOMEM );
  ret->length = numDet;

  for ( UINT4 X=0; X < numDet; X ++ )
    {
      XLAL_CHECK_NULL ( (ret->data[X] = XLALMakeTimestamps ( tStart, Tspan, Tsft, Toverlap ) ) != NULL, XLAL_EFUNC );
    } // for X < numDet

  return ret;

} /* XLALMakeMultiTimestamps() */


/**
 * Extract timstamps-vector from the given SFTVector
 */
LIGOTimeGPSVector *
XLALExtractTimestampsFromSFTs ( const SFTVector *sfts )		/**< [in] input SFT-vector  */
{
  /* check input consistency */
  if ( !sfts ) {
    XLALPrintError ("%s: invalid NULL input 'sfts'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  UINT4 numSFTs = sfts->length;
  /* create output vector */
  LIGOTimeGPSVector *ret = NULL;
  if ( ( ret = XLALCreateTimestampVector ( numSFTs )) == NULL ) {
    XLALPrintError ("%s: XLALCreateTimestampVector(%d) failed.\n", __func__, numSFTs );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  REAL8 Tsft = 1.0 / sfts->data[0].deltaF;
  ret->deltaT = Tsft;

  UINT4 i;
  for ( i=0; i < numSFTs; i ++ )
    ret->data[i] = sfts->data[i].epoch;

  /* done: return Ts-vector */
  return ret;

} /* XLALExtractTimestampsFromSFTs() */

/**
 * Given a multi-SFT vector, return a MultiLIGOTimeGPSVector holding the
 * SFT timestamps
 */
MultiLIGOTimeGPSVector *
XLALExtractMultiTimestampsFromSFTs ( const MultiSFTVector *multiSFTs )
{
  /* check input consistency */
  if ( !multiSFTs || multiSFTs->length == 0 ) {
    XLALPrintError ("%s: illegal NULL or empty input 'multiSFTs'.\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  UINT4 numIFOs = multiSFTs->length;

  /* create output vector */
  MultiLIGOTimeGPSVector *ret = NULL;
  if ( (ret = XLALCalloc ( 1, sizeof(*ret) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( 1, %d ).\n", __func__, sizeof(*ret));
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  if ( (ret->data = XLALCalloc ( numIFOs, sizeof(*ret->data) )) == NULL ) {
    XLALPrintError ("%s: failed to XLALCalloc ( %d, %d ).\n", __func__, numIFOs, sizeof(ret->data[0]) );
    XLALFree (ret);
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }
  ret->length = numIFOs;

  /* now extract timestamps vector from each SFT-vector */
  UINT4 X;
  for ( X=0; X < numIFOs; X ++ )
    {
      if ( (ret->data[X] = XLALExtractTimestampsFromSFTs ( multiSFTs->data[X] )) == NULL ) {
        XLALPrintError ("%s: XLALExtractTimestampsFromSFTs() failed for X=%d\n", __func__, X );
        XLALDestroyMultiTimestamps ( ret );
        XLAL_ERROR_NULL ( XLAL_EFUNC );
      }

    } /* for X < numIFOs */

  return ret;

} /* XLALExtractMultiTimestampsFromSFTs() */


/**
 * Extract timstamps-vector from the given SFTVector
 */
LIGOTimeGPSVector *
XLALTimestampsFromSFTCatalog ( const SFTCatalog *catalog )		/**< [in] input SFT-catalog  */
{
  /* check input consistency */
  XLAL_CHECK_NULL ( catalog != NULL, XLAL_EINVAL );

  UINT4 numSFTs = catalog->length;

  /* create output vector */
  LIGOTimeGPSVector *ret;
  XLAL_CHECK_NULL ( ( ret = XLALCreateTimestampVector ( numSFTs )) != NULL, XLAL_EINVAL, "Failed to XLALCreateTimestampVector ( %d )\n", numSFTs );

  REAL8 Tsft = 1.0 / catalog->data[0].header.deltaF;
  ret->deltaT = Tsft;

  for ( UINT4 i = 0; i < numSFTs; i ++ ) {
    ret->data[i] = catalog->data[i].header.epoch;
  }

  /* done: return Ts-vector */
  return ret;

} /* XLALTimestampsFromSFTCatalog() */


/**
 * Given a multi-SFTCatalogView, return a MultiLIGOTimeGPSVector holding the
 * SFT timestamps
 */
MultiLIGOTimeGPSVector *
XLALTimestampsFromMultiSFTCatalogView ( const MultiSFTCatalogView *multiView )
{
  /* check input consistency */
  XLAL_CHECK_NULL ( multiView != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( multiView->length > 0, XLAL_EINVAL );

  UINT4 numIFOs = multiView->length;

  /* create output vector */
  MultiLIGOTimeGPSVector *ret;
  XLAL_CHECK_NULL ( (ret = XLALCalloc ( 1, sizeof(*ret) )) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( (ret->data = XLALCalloc ( numIFOs, sizeof(*(ret->data)) )) != NULL, XLAL_ENOMEM );
  ret->length = numIFOs;

  /* now extract timestamps vector from each IFO's SFT-Catalog */
  for ( UINT4 X=0; X < numIFOs; X ++ )
    {
      XLAL_CHECK_NULL ( (ret->data[X] = XLALTimestampsFromSFTCatalog ( &(multiView->data[X]) )) != NULL, XLAL_EFUNC );
    } /* for X < numIFOs */

  return ret;

} /* XLALTimestampsFromMultiSFTCatalogView() */


/**
 * Destroy a MultiLIGOTimeGPSVector timestamps vector
 */
void
XLALDestroyMultiTimestamps ( MultiLIGOTimeGPSVector *multiTS )
{
  UINT4 numIFOs, X;

  if ( !multiTS )
    return;

  numIFOs = multiTS->length;
  for ( X=0; X < numIFOs; X ++ ) {
    XLALDestroyTimestampVector ( multiTS->data[X] );
  }

  XLALFree ( multiTS->data );
  XLALFree ( multiTS );

  return;

} /* XLALDestroyMultiTimestamps() */


/**
 * Extract/construct the unique 2-character "channel prefix" from the given
 * "detector-name", which unfortunately will not always follow any of the
 * official detector-naming conventions given in the Frames-Spec LIGO-T970130-F-E
 * This function therefore sometime has to do some creative guessing:
 *
 * NOTE: in case the channel-number can not be deduced from the name,
 * it is set to '1', and a warning will be printed if lalDebugLevel > 0.
 *
 * NOTE2: the returned string is allocated here!
 *
 * Note3: if more than one valid detector-string is found in the input, an error is returned
 *
 */
CHAR *
XLALGetChannelPrefix ( const CHAR *name )
{
  CHAR *channel = XLALCalloc( 3, sizeof(CHAR) );  /* 2 chars + \0 */

#define CHECK_UNIQUE do { if ( channel[0] != 0 ) XLAL_ERROR_NULL ( XLAL_EINVAL, "More than one matching detector name found in '%s'", name ); } while(0)

  if ( !channel ) {
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "Failed to calloc(3)!\n" );
  }
  if ( !name ) {
    XLALFree ( channel );
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Invalid NULL input 'name'" );
  }

  /* first handle (currently) unambiguous ones */
  if ( strstr( name, "ALLEGRO") || strstr ( name, "A1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "A1");
  }
  if ( strstr(name, "NIOBE") || strstr( name, "B1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "B1");
  }
  if ( strstr(name, "EXPLORER") || strstr( name, "E1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "E1");
  }
  if ( strstr(name, "GEO") || strstr(name, "G1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "G1" );
  }
  if ( strstr(name, "ACIGA") || strstr (name, "K1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "K1" );
  }
  if ( strstr(name, "LLO") || strstr(name, "Livingston") || strstr(name, "L1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "L1" );
  }
  if ( strstr(name, "Nautilus") || strstr(name, "N1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "N1" );
  }
  if ( strstr(name, "AURIGA") || strstr(name,"O1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "O1" );
  }
  if ( strstr(name, "CIT_40") || strstr(name, "Caltech-40") || strstr(name, "P1") ) {
    CHECK_UNIQUE;
    strcpy ( channel, "P1" );
  }
  if ( strstr(name, "TAMA") || strstr(name, "T1") ) {
    CHECK_UNIQUE;
    strcpy (channel, "T1" );
  }
  /* currently the only real ambiguity arises with H1 vs H2 */
  if ( strstr(name, "LHO") || strstr(name, "Hanford") || strstr(name, "H1") || strstr(name, "H2") ) {
    if ( strstr(name, "LHO_2k") ||  strstr(name, "H2") )
      {
        CHECK_UNIQUE;
        strcpy ( channel, "H2" );
      }
    if ( strstr(name, "LHO_4k") ||  strstr(name, "H1") )
      {
        CHECK_UNIQUE;
        strcpy ( channel, "H1" );
      }
    /* otherwise: guess */
    if ( channel[0] == 0 )
      {
        strcpy ( channel, "H1" );
        XLALPrintWarning("Detector-name '%s' ambiguous, guessing '%s'\n", name, channel );
      }
  } /* if LHO */
  /* LISA channel names are simply left unchanged */
  if ( strstr(name, "Z1") || strstr(name, "Z2") || strstr(name, "Z3")
       || strstr(name, "Z4") || strstr(name, "Z5") || strstr(name, "Z6")
       || strstr(name, "Z7") || strstr(name, "Z8") || strstr(name, "Z9") )
    {
      CHECK_UNIQUE;
      strncpy ( channel, name, 2);
      channel[2] = 0;
    }
  if ( strstr(name, "Virgo") || strstr(name, "VIRGO") || strstr(name, "V1") || strstr(name, "V2") )
    {
      if ( strstr(name, "Virgo_CITF") || strstr(name, "V2") )
        {
          CHECK_UNIQUE;
          strcpy ( channel, "V2" );
        }
      if ( strstr(name, "Virgo") || strstr(name, "VIRGO") || strstr(name, "V1") )
        {
          CHECK_UNIQUE;
          strcpy ( channel, "V1" );
        }
    } /* if Virgo */

  /* Did we fail to find any matches? */
  if ( channel[0] == 0 )
    XLAL_ERROR_NULL ( XLAL_EINVAL, "Unknown detector-name '%s'", name );
  else
    return channel;

} /* XLALGetChannelPrefix() */


/**
 * Find the site geometry-information 'LALDetector' (mis-nomer!) given a detector-name.
 * The LALDetector struct is allocated here.
 */
LALDetector *
XLALGetSiteInfo ( const CHAR *name )
{
  CHAR *channel;
  LALDetector *site;

  /* first turn the free-form 'detector-name' into a well-defined channel-prefix */
  if ( ( channel = XLALGetChannelPrefix ( name ) ) == NULL ) {
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }

  if ( ( site = XLALCalloc ( 1, sizeof( *site) )) == NULL ) {
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  }

  switch ( channel[0] )
    {
    case 'T':
      (*site) = lalCachedDetectors[LAL_TAMA_300_DETECTOR];
      break;
    case 'V':
      (*site) = lalCachedDetectors[LAL_VIRGO_DETECTOR];
      break;
    case 'G':
      (*site) = lalCachedDetectors[LAL_GEO_600_DETECTOR];
      break;
    case 'H':
      if ( channel[1] == '1' )
	(*site) = lalCachedDetectors[LAL_LHO_4K_DETECTOR];
      else
	(*site) = lalCachedDetectors[LAL_LHO_2K_DETECTOR];
      break;
    case 'L':
      (*site) = lalCachedDetectors[LAL_LLO_4K_DETECTOR];
      break;
    case 'P':
      (*site) = lalCachedDetectors[LAL_CIT_40_DETECTOR];
      break;
    case 'N':
      (*site) = lalCachedDetectors[LAL_NAUTILUS_DETECTOR];
      break;

    case 'Z':       /* create dummy-sites for LISA  */
      if ( XLALcreateLISA ( site, channel[1] ) != 0 )
	{
	  XLALPrintError("\nFailed to created LISA detector '%d'\n\n", channel[1]);
	  XLALFree ( site );
	  XLALFree ( channel );
	  XLAL_ERROR_NULL ( XLAL_EFUNC );
	}
      break;

    default:
      XLALPrintError ( "\nSorry, I don't have the site-info for '%c%c'\n\n", channel[0], channel[1]);
      XLALFree(site);
      XLALFree(channel);
      XLAL_ERROR_NULL ( XLAL_EINVAL );
      break;
    } /* switch channel[0] */

  XLALFree ( channel );

  return site;

} /* XLALGetSiteInfo() */


/**
 * Computes weight factors arising from MultiSFTs with different noise
 * floors
 */
MultiNoiseWeights *
XLALComputeMultiNoiseWeights ( const MultiPSDVector  *rngmed,
			       UINT4                 blocksRngMed,
			       UINT4                 excludePercentile)
{
  XLAL_CHECK_NULL ( rngmed != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( rngmed->data != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( rngmed->length != 0, XLAL_EINVAL );

  UINT4 numIFOs = rngmed->length;
  REAL8 Tsft = 1.0 / rngmed->data[0]->data[0].deltaF;

  /* create multi noise weights for output */
  MultiNoiseWeights *multiWeights = NULL;
  XLAL_CHECK_NULL ( (multiWeights = XLALCalloc(1, sizeof(*multiWeights))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( (multiWeights->data = XLALCalloc ( numIFOs, sizeof(*multiWeights->data))) != NULL, XLAL_ENOMEM );
  multiWeights->length = numIFOs;

  UINT4 numSFTsTot = 0;
  REAL8 sumWeights = 0;

  for ( UINT4 X = 0; X < numIFOs; X++)
    {
      UINT4 numSFTs = rngmed->data[X]->length;
      numSFTsTot += numSFTs;

      /* create k^th weights vector */
      if( ( multiWeights->data[X] = XLALCreateREAL8Vector ( numSFTs ) ) == NULL )
        {
          /* free weights vectors created previously in loop */
          XLALDestroyMultiNoiseWeights ( multiWeights );
          XLAL_ERROR_NULL ( XLAL_EFUNC, "Failed to allocate noiseweights for IFO X = %d\n", X );
        } /* if XLALCreateREAL8Vector() failed */

      /* loop over rngmeds and calculate weights -- one for each sft */
      for ( UINT4 alpha = 0; alpha < numSFTs; alpha++)
	{
	  UINT4 halfBlock = blocksRngMed/2;

	  REAL8FrequencySeries *thisrm = &(rngmed->data[X]->data[alpha]);

	  UINT4 lengthsft = thisrm->data->length;

	  XLAL_CHECK_NULL ( lengthsft >= blocksRngMed, XLAL_EINVAL );

	  UINT4 length = lengthsft - blocksRngMed + 1;
	  UINT4 halfLength = length/2;

	  /* calculate index in power medians vector from which to calculate mean */
	  UINT4 excludeIndex =  excludePercentile * halfLength ; /* integer arithmetic */
	  excludeIndex /= 100; /* integer arithmetic */

	  REAL8 Tsft_avgS2 = 0.0;	// 'S2' refers to double-sided PSD
	  for ( UINT4 k = halfBlock + excludeIndex; k < lengthsft - halfBlock - excludeIndex; k++)
	    {
	      Tsft_avgS2 += thisrm->data->data[k];
	    }
	  Tsft_avgS2 /= lengthsft - 2*halfBlock - 2*excludeIndex;

          REAL8 wXa = 1.0/Tsft_avgS2;	// unnormalized weight
	  multiWeights->data[X]->data[alpha] = wXa;

	  sumWeights += wXa;	// sum the weights to normalize this at the end
	} /* end loop over sfts for each ifo */

    } /* end loop over ifos */

  /* overall noise-normalization factor Sinv = 1/Nsft sum_Xa Sinv_Xa,
   * see Eq.(60) in CFSv2 notes:
   * https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=1665&version=3
   */
  REAL8 TsftS2_inv = sumWeights / numSFTsTot;	// this is double-sided PSD 'S2'

  /* make weights of order unity by normalizing with TsftS2_inv, see Eq.(58) in CFSv2 notes (v3) */
  for ( UINT4 X = 0; X < numIFOs; X ++) {
    UINT4 numSFTs = multiWeights->data[X]->length;
    for ( UINT4 alpha = 0; alpha < numSFTs; alpha ++)
      {
	multiWeights->data[X]->data[alpha] /= TsftS2_inv;
      }
  }

  multiWeights->Sinv_Tsft = 0.5 * Tsft*Tsft * TsftS2_inv;		/* 'Sinv * Tsft' refers to single-sided PSD!! Eq.(60) in CFSv2 notes (v3)*/

  return multiWeights;

} /* XLALComputeMultiNoiseWeights() */

/** Destroy a MultiNoiseWeights object */
void
XLALDestroyMultiNoiseWeights ( MultiNoiseWeights *weights )
{
  if ( weights == NULL)
    return;

  for ( UINT4 k = 0; k < weights->length; k++ )
    XLALDestroyREAL8Vector ( weights->data[k] );

  XLALFree ( weights->data );
  XLALFree ( weights );

  return;

} /* XLALDestroyMultiNoiseWeights() */


/**
 * Interpolate frequency-series to newLen frequency-bins.
 * This is using DFT-interpolation (derived from zero-padding).
 */
COMPLEX8Vector *
XLALrefineCOMPLEX8Vector (const COMPLEX8Vector *in,
			  UINT4 refineby,
			  UINT4 Dterms)
{
  UINT4 newLen, oldLen, l;
  COMPLEX8Vector *ret = NULL;

  if ( !in )
    return NULL;

  oldLen = in->length;
  newLen = oldLen * refineby;

  /* the following are used to speed things up in the innermost loop */
  if ( (ret = XLALCreateCOMPLEX8Vector ( newLen )) == NULL )
    return NULL;

  for (l=0; l < newLen; l++)
    {

      REAL8 kappa_l_k;
      REAL8 remain, kstarREAL;
      UINT4 kstar, kmin, kmax, k;
      REAL8 sink, coskm1;
      REAL8 Yk_re, Yk_im, Xd_re, Xd_im;

      kstarREAL = 1.0 * l  / refineby;
      kstar = (INT4)( kstarREAL + 0.5);	/* round to closest bin */
      kstar = MIN ( kstar, oldLen - 1 );	/* stay within the old SFT index-bounds */
      remain = kstarREAL - kstar;

      /* boundaries for innermost loop */
      kmin = MAX( 0, (INT4)kstar - (INT4)Dterms );
      kmax = MIN( oldLen, kstar + Dterms );

      Yk_re = Yk_im = 0;
      if ( fabs(remain) > 1e-5 )	/* denominater doens't vanish */
	{
	  /* Optimization: sin(2pi*kappa(l,k)) = sin(2pi*kappa(l,0) and idem for cos */
	  sink = sin ( LAL_TWOPI * remain );
	  coskm1 = cos ( LAL_TWOPI * remain ) - 1.0;

	  /* ---------- innermost loop: k over 2*Dterms around kstar ---------- */
	  for (k = kmin; k < kmax; k++)
	    {
	      REAL8 Plk_re, Plk_im;

	      Xd_re = crealf(in->data[k]);
	      Xd_im = cimagf(in->data[k]);

	      kappa_l_k = kstarREAL - k;

	      Plk_re = sink / kappa_l_k;
	      Plk_im = coskm1 / kappa_l_k;

	      Yk_re += Plk_re * Xd_re - Plk_im * Xd_im;
	      Yk_im += Plk_re * Xd_im + Plk_im * Xd_re;

	    } /* hotloop over Dterms */
	}
      else	/* kappa -> 0: Plk = 2pi delta(k, l) */
	{
	  Yk_re = LAL_TWOPI * crealf(in->data[kstar]);
	  Yk_im = LAL_TWOPI * cimagf(in->data[kstar]);
	}

      const REAL8 OOTWOPI = (1.0 / LAL_TWOPI );
      ret->data[l].realf_FIXME = OOTWOPI* Yk_re;
      ret->data[l].imagf_FIXME = OOTWOPI * Yk_im;

    }  /* for l < newlen */

  return ret;

} /* XLALrefineCOMPLEX8Vector() */


/**
 * Function to read a segment list from given filename, returns a *sorted* SegmentList
 *
 * The segment-list format parse here is consistent with Xavie's segment lists used previously
 * and follows the format <repeated lines of form "startGPS endGPS duration[h] NumSFTs">,
 * allowed comment-characters are '%' and '#'
 *
 * \note we (ab)use the integer 'id' field in LALSeg to carry the total number of SFTs
 * contained in that segment. This can be used as a consistency check when loading SFTs for these segments.
 *
 */
LALSegList *
XLALReadSegmentsFromFile ( const char *fname	/**< name of file containing segment list */
                           )
{
  LALSegList *segList = NULL;

  /** check input consistency */
  if ( !fname ) {
    XLALPrintError ( "%s: NULL input 'fname'", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  /* read and parse segment-list file contents*/
  LALParsedDataFile *flines = NULL;
  if ( XLALParseDataFile ( &flines, fname ) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  UINT4 numSegments = flines->lines->nTokens;
  /* allocate and initialized segment list */
  if ( (segList = XLALCalloc ( 1, sizeof(*segList) )) == NULL )
    XLAL_ERROR_NULL ( XLAL_ENOMEM );
  if ( XLALSegListInit ( segList ) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );


  UINT4 iSeg;
  for ( iSeg = 0; iSeg < numSegments; iSeg ++ )
    {
      REAL8 t0, t1, TspanHours;
      INT4 NSFT;
      LALSeg thisSeg;
      int ret;
      ret = sscanf ( flines->lines->tokens[iSeg], "%lf %lf %lf %d", &t0, &t1, &TspanHours, &NSFT );
      if ( ret != 4 ) {
        XLALPrintError ("%s: failed to parse data-line %d (%d) in segment-list %s: '%s'\n", __func__, iSeg, ret, fname, flines->lines->tokens[iSeg] );
        XLALSegListClear ( segList );
        XLALFree ( segList );
        XLALDestroyParsedDataFile ( flines );
        XLAL_ERROR_NULL ( XLAL_ESYS );
      }
      /* check internal consistency of these numbers */
      REAL8 hours = 3600.0;
      if ( fabs ( t1 - t0 - TspanHours * hours ) >= 1.0 ) {
        XLALPrintError ("%s: Inconsistent segment list, in line %d: t0 = %f, t1 = %f, Tspan = %f != t1 - t0 (to within 1s)\n", __func__, iSeg, t0, t1, TspanHours );
        XLAL_ERROR_NULL ( XLAL_EDOM );
      }

      LIGOTimeGPS start, end;
      XLALGPSSetREAL8( &start, t0 );
      XLALGPSSetREAL8( &end,   t1 );

      /* we set number of SFTs as 'id' field, as we have no other use for it */
      if ( XLALSegSet ( &thisSeg, &start, &end, NSFT ) != XLAL_SUCCESS )
        XLAL_ERROR_NULL ( XLAL_EFUNC );

      if ( XLALSegListAppend ( segList, &thisSeg ) != XLAL_SUCCESS )
        XLAL_ERROR_NULL ( XLAL_EFUNC );

    } /* for iSeg < numSegments */

  /* sort final segment list in increasing GPS start-times */
  if ( XLALSegListSort( segList ) != XLAL_SUCCESS )
    XLAL_ERROR_NULL ( XLAL_EFUNC );

  /* free parsed segment file contents */
  XLALDestroyParsedDataFile ( flines );

  return segList;

} /* XLALReadSegmentsFromFile() */

/**
 * Return a vector of SFTs containing only the bins in [fMin, fMin+Band].
 * Note: the output SFT is guaranteed to "cover" the input boundaries 'fMin'
 * and 'fMin+Band', ie if necessary the output SFT contains one additional
 * bin on either end of the interval.
 *
 * This uses the conventions in XLALFindCoveringSFTBins() to determine
 * the 'effective' frequency-band to extract.
 *
 */
SFTVector *
XLALExtractBandFromSFTVector ( const SFTVector *inSFTs, REAL8 fMin, REAL8 Band )
{
  XLAL_CHECK_NULL ( inSFTs != NULL, XLAL_EINVAL, "Invalid NULL input SFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( inSFTs->length > 0, XLAL_EINVAL, "Invalid zero-length input SFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( fMin >= 0, XLAL_EDOM, "Invalid negative frequency fMin = %g\n", fMin );
  XLAL_CHECK_NULL ( Band > 0, XLAL_EDOM, "Invalid non-positive Band = %g\n", Band );

  UINT4 numSFTs = inSFTs->length;

  REAL8 df      = inSFTs->data[0].deltaF;
  REAL8 Tsft    = 1.0 / df;

  REAL8 fMinSFT    = inSFTs->data[0].f0;
  UINT4 numBinsSFT = inSFTs->data[0].data->length;
  REAL8 BandSFT    = df * numBinsSFT;
  UINT4 firstBinSFT= round ( fMinSFT / df );	// round to closest bin
  UINT4 lastBinSFT = firstBinSFT + ( numBinsSFT - 1 );

  // find 'covering' SFT-band to extract
  UINT4 firstBinExt, numBinsExt;
  XLAL_CHECK_NULL ( XLALFindCoveringSFTBins ( &firstBinExt, &numBinsExt, fMin, Band, Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );
  UINT4 lastBinExt = firstBinExt + ( numBinsExt - 1 );

  XLAL_CHECK_NULL ( firstBinExt >= firstBinSFT && (lastBinExt <= lastBinSFT), XLAL_EINVAL,
                    "Requested frequency-bins [%f,%f]Hz = [%d, %d] not contained within SFT's [%f, %f]Hz = [%d,%d].\n",
                    fMin, fMin + Band, firstBinExt, lastBinExt, fMinSFT, fMinSFT + BandSFT, firstBinSFT, lastBinSFT );

  INT4 firstBinOffset = firstBinExt - firstBinSFT;

  SFTVector *ret;
  XLAL_CHECK_NULL ( (ret = XLALCreateSFTVector ( numSFTs, numBinsExt )) != NULL, XLAL_EFUNC );

  for ( UINT4 i = 0; i < numSFTs; i ++ )
    {
      SFTtype *dest = &(ret->data[i]);
      SFTtype *src =  &(inSFTs->data[i]);
      COMPLEX8Vector *ptr = dest->data;

      /* copy complete header first */
      memcpy ( dest, src, sizeof(*dest) );
      /* restore data-pointer */
      dest->data = ptr;
      /* set correct fMin */
      dest->f0 = firstBinExt * df ;

      /* copy the relevant part of the data */
      memcpy ( dest->data->data, src->data->data + firstBinOffset, numBinsExt * sizeof( dest->data->data[0] ) );

    } /* for i < numSFTs */

  /* return final SFT-vector */
  return ret;

} /* XLALExtractBandFromSFTVector() */


/**
 * Adds SFT-data from MultiSFTvector 'b' to elements of MultiSFTVector 'a'
 *
 * NOTE: the inputs 'a' and 'b' must have consistent number of IFO, number of SFTs,
 * IFO-names, start-frequency, frequency-spacing, timestamps, units and number of bins.
 *
 * The 'name' field of input/output SFTs in 'a' is not modified!
 */
int
XLALMultiSFTVectorAdd ( MultiSFTVector *a,	/**< [in/out] MultiSFTVector to be added to */
                        const MultiSFTVector *b	/**< [in] MultiSFTVector data to be added */
                        )
{
  XLAL_CHECK ( a != NULL, XLAL_EINVAL );
  XLAL_CHECK ( b != NULL, XLAL_EINVAL );

  XLAL_CHECK ( a->length == b->length, XLAL_EINVAL );
  UINT4 numIFOs = a->length;

  for ( UINT4 X = 0; X < numIFOs; X ++ )
    {
      SFTVector *vect1 = a->data[X];
      SFTVector *vect2 = b->data[X];

      XLAL_CHECK ( XLALSFTVectorAdd ( vect1, vect2 ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALSFTVectorAdd() failed for SFTVector %d out of %d\n", X+1, numIFOs );

    } // for X < numIFOs

  return XLAL_SUCCESS;

} /* XLALMultiSFTVectorAdd() */


/**
 * Adds SFT-data from SFTvector 'b' to elements of SFTVector 'a'
 *
 * NOTE: the inputs 'a' and 'b' must have consistent number of SFTs, IFO-names,
 * start-frequency, frequency-spacing, timestamps, units and number of bins.
 *
 * The 'name' field of input/output SFTs in 'a' is not modified!
 */
int
XLALSFTVectorAdd ( SFTVector *a,	/**< [in/out] SFTVector to be added to */
                   const SFTVector *b	/**< [in] SFTVector data to be added */
                   )
{
  XLAL_CHECK ( a != NULL, XLAL_EINVAL );
  XLAL_CHECK ( b != NULL, XLAL_EINVAL );

  XLAL_CHECK ( a->length == b->length, XLAL_EINVAL );
  UINT4 numSFTs = a->length;

  for ( UINT4 k = 0; k < numSFTs; k ++ )
    {
      SFTtype *sft1 = &(a->data[k]);
      SFTtype *sft2 = &(b->data[k]);

      XLAL_CHECK ( XLALSFTAdd ( sft1, sft2 ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALSFTAdd() failed for SFTs k = %d out of %d SFTs\n", k, numSFTs );

    } // for k < numSFTs

  return XLAL_SUCCESS;
} /* XLALSFTVectorAdd() */


/**
 * Adds SFT-data from SFT 'b' to SFT 'a'
 *
 * NOTE: the inputs 'a' and 'b' must have consistent IFO-names,
 * start-frequency, frequency-spacing, timestamps, units and number of bins.
 *
 * The 'name' field of input/output SFTs in 'a' is not modified!
 */
int
XLALSFTAdd ( SFTtype *a,		/**< [in/out] SFT to be added to */
             const SFTtype *b	/**< [in] SFT data to be added */
             )
{
  XLAL_CHECK ( a != NULL, XLAL_EINVAL );
  XLAL_CHECK ( b != NULL, XLAL_EINVAL );
  XLAL_CHECK ( a->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( b->data != NULL, XLAL_EINVAL );

  XLAL_CHECK ( strncmp ( a->name, b->name, 2 ) == 0, XLAL_EINVAL, "SFT detectors differ '%c%c' != '%c%c'\n", a->name[0], a->name[1], b->name[0], b->name[1] );
  XLAL_CHECK ( XLALGPSDiff ( &(a->epoch), &(b->epoch) ) == 0, XLAL_EINVAL, "SFT epochs differ %ld != %ld ns\n", XLALGPSToINT8NS ( &(a->epoch) ), XLALGPSToINT8NS ( &(b->epoch) ) );

  REAL8 tol = 10 * LAL_REAL8_EPS;	// generously allow up to 10*eps tolerance
  XLAL_CHECK ( gsl_fcmp ( a->f0, b->f0, tol ) == 0, XLAL_ETOL, "SFT frequencies relative deviation exceeds %g: %.16g != %.16g\n", tol, a->f0, b->f0 );
  XLAL_CHECK ( gsl_fcmp ( a->deltaF, b->deltaF, tol ) == 0, XLAL_ETOL, "SFT frequency-steps relative deviation exceeds %g: %.16g != %.16g\n", tol, a->deltaF, b->deltaF );
  XLAL_CHECK ( XLALUnitCompare ( &(a->sampleUnits), &(b->sampleUnits) ) == 0, XLAL_EINVAL, "SFT sample units differ\n" );
  XLAL_CHECK ( a->data->length == b->data->length, XLAL_EINVAL, "SFT lengths differ: %d != %d\n", a->data->length, b->data->length );

  UINT4 numBins = a->data->length;
  for ( UINT4 k = 0; k < numBins; k ++ )
    {
      a->data->data[k].realf_FIXME += b->data->data[k].realf_FIXME;
      a->data->data[k].imagf_FIXME += b->data->data[k].imagf_FIXME;
    }

  return XLAL_SUCCESS;

} /* XLALSFTAdd() */

/**
 * Return the 'effective' frequency-band [fMinEff, fMaxEff] = [firstBin, lastBin] * 1/Tsft,
 * with numBins = lastBin - firstBin + 1
 * which is the smallest band of SFT-bins that fully covers a given band [fMin, fMin+Band]
 *
 * ==> calculate "effective" fMinEff by rounding down from fMin to closest (firstBin/Tsft)
 * and rounds up in the same way to fMaxEff = (lastBin/Tsft).
 *
 * The 'fudge region' allowing for numerical noise is eps= 10*LAL_REAL8_EPS ~2e-15
 * relative deviation: ie if the SFT contains a bin at 'fi', then we consider for example
 * "fMin == fi" if  fabs(fi - fMin)/fi < eps.
 */
int
XLALFindCoveringSFTBins ( UINT4 *firstBin,	///< [out] effective lower frequency-bin fMinEff = firstBin/Tsft
                          UINT4 *numBins,	///< [out] effective Band of SFT-bins, such that BandEff = (numBins-1)/Tsft
                          REAL8 fMinIn,		///< [in] input lower frequency
                          REAL8 BandIn,		///< [in] input frequency band
                          REAL8 Tsft		///< [in] SFT duration 'Tsft'
                          )
{
  XLAL_CHECK ( firstBin != NULL, XLAL_EINVAL );
  XLAL_CHECK ( numBins  != NULL, XLAL_EINVAL );
  XLAL_CHECK ( fMinIn >= 0, XLAL_EDOM );
  XLAL_CHECK ( BandIn >= 0, XLAL_EDOM );
  XLAL_CHECK ( Tsft > 0, XLAL_EDOM );

  volatile REAL8 dFreq = 1.0 / Tsft;
  volatile REAL8 tmp;
  // NOTE: don't "simplify" this: we try to make sure
  // the result of this will be guaranteed to be IEEE-compliant,
  // and identical to other locations, such as in SFT-IO

  // ----- lower effective frequency
  tmp = fMinIn / dFreq;
  UINT4 imin = (UINT4) floor( tmp * fudge_up );	// round *down*, allowing for eps 'fudge'

  // ----- upper effective frequency
  REAL8 fMaxIn = fMinIn + BandIn;
  tmp = fMaxIn / dFreq;
  UINT4 imax = (UINT4) ceil ( tmp * fudge_down );  // round *up*, allowing for eps fudge

  // ----- effective band
  UINT4 num_bins = (UINT4) (imax - imin + 1);

  // ----- return these
  (*firstBin) = imin;
  (*numBins)  = num_bins;

  return XLAL_SUCCESS;

} // XLALFindCoveringSFTBins()
