/*
 * Copyright (C) 2019, 2020 David Keitel
 * Copyright (C) 2019 Pep Covas
 * Copyright (C) 2010, 2013--2016, 2022 Karl Wette
 * Copyright (C) 2010 Chris Messenger
 * Copyright (C) 2004--2009, 2013, 2014 Reinhard Prix
 * Copyright (C) 2004--2006, 2010 Bernd Machenschalk
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

#include <gsl/gsl_math.h>

#include <lal/Units.h>
#include <lal/FrequencySeries.h>

#include "SFTinternal.h"

/*========== function definitions ==========*/

/// \addtogroup SFTfileIO_h
/// @{

/**
 * Defines the official CW convention for whether a GPS time is 'within' a given range, defined
 * as the half-open interval [minGPS, maxGPS)
 *
 * This function should be used when dealing with SFTs, segment lists, etc. It returns:
 *
 * - -1 if \f$ \mathtt{gps} < \mathtt{minGPS} \f$ , i.e. GPS time is 'below' the range;
 * -  0 if \f$ \mathtt{minGPS} \le \mathtt{gps} < \mathtt{maxGPS} \f$ , i.e. GPS time is 'within' the range;
 * -  1 if \f$ \mathtt{maxGPS} \le \mathtt{gos} \f$ , i.e. GPS time is 'above' the range;
 *
 * If either \c minGPS or \c maxGPS are \c NULL, there are treated as \f$ -\infty \f$ or \f$ +\infty \f$ respectively.
 */
int XLALCWGPSinRange( const LIGOTimeGPS gps, const LIGOTimeGPS* minGPS, const LIGOTimeGPS* maxGPS )
{
  if (minGPS != NULL && GPS2REAL8(gps) < GPS2REAL8(*minGPS)) {
    return -1;
  }
  if (maxGPS != NULL && GPS2REAL8(gps) >= GPS2REAL8(*maxGPS)) {
    return 1;
  }
  return 0;
} /* XLALCWGPSinRange() */


/**
 * Round a REAL8 frequency down to the nearest integer SFT bin number.
 *
 * This function provides an official rounding convention,
 * including a "fudge" factor.
 */
UINT4 XLALRoundFrequencyDownToSFTBin( const REAL8 freq, const REAL8 df )
{
  return (UINT4) floor (freq / df * fudge_up);
} /* XLALRoundFrequencyDownToSFTBin() */


/**
 * Round a REAL8 frequency up to the nearest integer SFT bin number.
 *
 * This function provides an official rounding convention,
 * including a "fudge" factor.
 */
UINT4 XLALRoundFrequencyUpToSFTBin( const REAL8 freq, const REAL8 df )
{
  return (UINT4) ceil (freq / df * fudge_down);
} /* XLALRoundFrequencyUpToSFTBin() */


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
 *
 * Note: this function is most useful for internal operations, e.g. where a generated
 * time series needs to be over-sampled to cover the SFT frequency band of interest.
 * Ultimately SFTs covering a half-open interval [fMinIn,BandIn) should be returned to
 * the user using XLALExtractStrictBandFromSFTVector().
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
    XLAL_ERROR_NULL ( XLAL_ENOMEM, "XLALCalloc (1, %zu) failed.\n", sizeof(*sft) );

  if ( numBins )
    {
      if ( (sft->data = XLALCreateCOMPLEX8Vector ( numBins )) == NULL )
        {
          XLALFree ( sft );
          XLAL_ERROR_NULL ( XLAL_EFUNC, "XLALCreateCOMPLEX8Vector ( %d ) failed. xlalErrno = %d\n", numBins, xlalErrno );
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
 * Copy an entire SFT-type into another.
 * We require the destination-SFT to have a NULL data-entry, as the
 * corresponding data-vector will be allocated here and copied into
 *
 * Note: the source-SFT is allowed to have a NULL data-entry,
 * in which case only the header is copied.
 */
int
XLALCopySFT ( SFTtype *dest, 		/**< [out] copied SFT (needs to be allocated already) */
              const SFTtype *src	/**< input-SFT to be copied */
              )
{
  // check input sanity
  XLAL_CHECK ( dest != NULL, XLAL_EINVAL );
  XLAL_CHECK ( dest->data == NULL, XLAL_EINVAL );
  XLAL_CHECK ( src != NULL, XLAL_EINVAL );

  /* copy complete head (including data-pointer, but this will be separately alloc'ed and copied in the next step) */
  (*dest) = (*src);	// struct copy

  /* copy data (if there's any )*/
  if ( src->data )
    {
      UINT4 numBins = src->data->length;
      XLAL_CHECK ( (dest->data = XLALCreateCOMPLEX8Vector ( numBins )) != NULL, XLAL_EFUNC );
      memcpy ( dest->data->data, src->data->data, numBins * sizeof (dest->data->data[0]));
    }

  return XLAL_SUCCESS;

} // XLALCopySFT()


/**
 * XLAL function to create an SFTVector of \c numSFT SFTs with \c SFTlen frequency-bins (which will be allocated too).
 */
SFTVector *
XLALCreateSFTVector ( UINT4 numSFTs, 	/**< number of SFTs */
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
 * XLAL function to create an SFTVector of \c numSFT SFTs (which are not allocated).
 */
SFTVector *
XLALCreateEmptySFTVector ( UINT4 numSFTs 	/**< number of SFTs */
                      )
{
  SFTVector *vect;

  if ( (vect = XLALCalloc ( 1, sizeof(*vect) )) == NULL ) {
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  vect->length = numSFTs;
  if ( (vect->data = XLALCalloc (1, numSFTs * sizeof ( *vect->data ) )) == NULL ) {
    XLALFree (vect);
    XLAL_ERROR_NULL( XLAL_ENOMEM );
  }

  return vect;

} /* XLALCreateEmptySFTVector() */


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
 * Create a complete copy of an SFT vector
 */
SFTVector *
XLALDuplicateSFTVector ( const SFTVector *sftsIn )
{
  XLAL_CHECK_NULL ( (sftsIn != NULL) && ( sftsIn->length > 0), XLAL_EINVAL );

  UINT4 numSFTs = sftsIn->length;
  UINT4 numBins = sftsIn->data[0].data->length;

  SFTVector *sftsOut;
  XLAL_CHECK_NULL ( (sftsOut = XLALCreateSFTVector ( numSFTs, numBins )) != NULL, XLAL_EFUNC );

  for ( UINT4 alpha=0; alpha < numSFTs; alpha++ )
    {
      SFTtype *thisSFTIn = &sftsIn->data[alpha];
      SFTtype *thisSFTOut = &sftsOut->data[alpha];

      COMPLEX8Vector *tmp = thisSFTOut->data;
      memcpy ( thisSFTOut, thisSFTIn, sizeof(*thisSFTOut) );
      thisSFTOut->data = tmp;
      thisSFTOut->data->length = numBins;
      memcpy ( thisSFTOut->data->data, thisSFTIn->data->data, numBins * sizeof(thisSFTOut->data->data[0]) );

    } // for alpha < numSFTs

  return sftsOut;

} // XLALDuplicateSFTVector()


/**
 * Create a multi-IFO SFT vector with a given number of bins per SFT and number of SFTs per IFO (which will be allocated too).
 *
 * Note that the input argument "length" refers to the number of frequency bins in each SFT.
 * The length of the returned MultiSFTVector (i.e. the number of IFOs)
 * is set from the length of the input numsft vector instead.
 */
MultiSFTVector *XLALCreateMultiSFTVector (
  UINT4 length,          /**< number of SFT data points (frequency bins) */
  UINT4Vector *numsft    /**< number of SFTs in each per-detector SFTVector */
  )
{

  XLAL_CHECK_NULL( length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( numsft != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( numsft->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( numsft->data != NULL, XLAL_EFAULT );

  MultiSFTVector *multSFTVec = NULL;

  XLAL_CHECK_NULL( ( multSFTVec = XLALCalloc( 1, sizeof(*multSFTVec) ) ) != NULL, XLAL_ENOMEM );

  const UINT4 numifo = numsft->length;
  multSFTVec->length = numifo;

  XLAL_CHECK_NULL( ( multSFTVec->data = XLALCalloc( numifo, sizeof(*multSFTVec->data) ) ) != NULL, XLAL_ENOMEM );

  for ( UINT4 k = 0; k < numifo; k++) {
    XLAL_CHECK_NULL( ( multSFTVec->data[k] = XLALCreateSFTVector( numsft->data[k], length ) ) != NULL, XLAL_ENOMEM );
  } /* loop over ifos */

  return multSFTVec;

} /* XLALCreateMultiSFTVector() */


/**
 * Create an empty multi-IFO SFT vector with a given number of SFTs per IFO (which are not allocated).
 */
MultiSFTVector *XLALCreateEmptyMultiSFTVector (
  UINT4Vector *numsft    /**< number of SFTs in each per-detector SFTVector */
  )
{
  XLAL_CHECK_NULL( numsft != NULL, XLAL_EFAULT );
  XLAL_CHECK_NULL( numsft->length > 0, XLAL_EINVAL );
  XLAL_CHECK_NULL( numsft->data != NULL, XLAL_EFAULT );

  MultiSFTVector *multSFTVec = NULL;

  XLAL_CHECK_NULL( ( multSFTVec = XLALCalloc( 1, sizeof(*multSFTVec) ) ) != NULL, XLAL_ENOMEM );

  const UINT4 numifo = numsft->length;
  multSFTVec->length = numifo;

  XLAL_CHECK_NULL( ( multSFTVec->data = XLALCalloc( numifo, sizeof(*multSFTVec->data) ) ) != NULL, XLAL_ENOMEM );

  for ( UINT4 k = 0; k < numifo; k++) {
    XLAL_CHECK_NULL( ( multSFTVec->data[k] = XLALCreateEmptySFTVector( numsft->data[k] ) ) != NULL, XLAL_ENOMEM );
  } /* loop over ifos */

  return multSFTVec;

} /* XLALCreateEmptyMultiSFTVector() */


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
 * Return an SFTs containing only the bins in [fMin, fMin+Band].
 * Note: the output SFT is guaranteed to "cover" the input boundaries 'fMin'
 * and 'fMin+Band', ie if necessary the output SFT contains one additional
 * bin on either end of the interval.
 *
 * This uses the conventions in XLALFindCoveringSFTBins() to determine
 * the 'effective' frequency-band to extract.
 *
 * \warning This convention is deprecated. Please use either
 * XLALExtractStrictBandFromSFT(), or else XLALSFTResizeBand()
 * if you really need a covering frequency band.
 */
int
XLALExtractBandFromSFT ( SFTtype **outSFT,	///< [out] output SFT (alloc'ed or re-alloced as required)
                         const SFTtype *inSFT,	///< [in] input SFT
                         REAL8 fMin,		///< [in] lower end of frequency interval to return
                         REAL8 Band		///< [in] band width of frequency interval to return
                         )
{
  XLAL_PRINT_DEPRECATION_WARNING("XLALExtractStrictBandFromSFT");

  XLAL_CHECK ( outSFT != NULL, XLAL_EINVAL );
  XLAL_CHECK ( inSFT != NULL, XLAL_EINVAL );
  XLAL_CHECK ( inSFT->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( fMin >= 0, XLAL_EDOM, "Invalid negative frequency fMin = %g\n", fMin );
  XLAL_CHECK ( Band > 0, XLAL_EDOM, "Invalid non-positive Band = %g\n", Band );

  REAL8 df      = inSFT->deltaF;
  REAL8 Tsft    = TSFTfromDFreq ( df );

  REAL8 fMinSFT    = inSFT->f0;
  UINT4 numBinsSFT = inSFT->data->length;
  UINT4 firstBinSFT= round ( fMinSFT / df );	// round to closest bin
  UINT4 lastBinSFT = firstBinSFT + ( numBinsSFT - 1 );

  // find 'covering' SFT-band to extract
  UINT4 firstBinExt, numBinsExt;
  XLAL_CHECK ( XLALFindCoveringSFTBins ( &firstBinExt, &numBinsExt, fMin, Band, Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );
  UINT4 lastBinExt = firstBinExt + ( numBinsExt - 1 );

  XLAL_CHECK ( firstBinExt >= firstBinSFT && (lastBinExt <= lastBinSFT), XLAL_EINVAL,
               "Requested frequency-bins [%f,%f]Hz = [%d, %d] not contained within SFT's [%f, %f]Hz = [%d,%d].\n",
               fMin, fMin + Band, firstBinExt, lastBinExt, fMinSFT, fMinSFT + (numBinsSFT-1) * df, firstBinSFT, lastBinSFT );

  INT4 firstBinOffset = firstBinExt - firstBinSFT;

  if ( (*outSFT) == NULL ) {
    XLAL_CHECK ( ((*outSFT) = XLALCalloc(1, sizeof(*(*outSFT)))) != NULL, XLAL_ENOMEM );
  }
  if ( (*outSFT)->data == NULL ) {
    XLAL_CHECK ( ((*outSFT)->data = XLALCreateCOMPLEX8Vector ( numBinsExt )) != NULL, XLAL_EFUNC );
  }
  if ( (*outSFT)->data->length != numBinsExt ) {
    XLAL_CHECK ( ((*outSFT)->data->data = XLALRealloc ( (*outSFT)->data->data, numBinsExt * sizeof((*outSFT)->data->data[0]))) != NULL, XLAL_ENOMEM );
    (*outSFT)->data->length = numBinsExt;
  }

  COMPLEX8Vector *ptr = (*outSFT)->data;	// keep copy to data-pointer
  (*(*outSFT)) = (*inSFT);			// copy complete header
  (*outSFT)->data = ptr;	  		// restore data-pointer
  (*outSFT)->f0 = firstBinExt * df ;	  	// set correct new fMin

  /* copy the relevant part of the data */
  memcpy ( (*outSFT)->data->data, inSFT->data->data + firstBinOffset, numBinsExt * sizeof( (*outSFT)->data->data[0] ) );

  return XLAL_SUCCESS;

} // XLALExtractBandFromSFT()


/**
 * Return a vector of SFTs containing only the bins in [fMin, fMin+Band].
 * Note: the output SFT is guaranteed to "cover" the input boundaries 'fMin'
 * and 'fMin+Band', ie if necessary the output SFT contains one additional
 * bin on either end of the interval.
 *
 * This uses the conventions in XLALFindCoveringSFTBins() to determine
 * the 'effective' frequency-band to extract.
 *
 * \warning This convention is deprecated. Please use either
 * XLALExtractStrictBandFromSFTVector(), or else XLALSFTVectorResizeBand()
 * if you really need a covering frequency band.
 */
SFTVector *
XLALExtractBandFromSFTVector ( const SFTVector *inSFTs,	///< [in] input SFTs
                               REAL8 fMin,		///< [in] lower end of frequency interval to return
                               REAL8 Band		///< [in] band width of frequency interval to return
                               )
{
  XLAL_PRINT_DEPRECATION_WARNING("XLALExtractStrictBandFromSFTVector");

  XLAL_CHECK_NULL ( inSFTs != NULL, XLAL_EINVAL, "Invalid NULL input SFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( inSFTs->length > 0, XLAL_EINVAL, "Invalid zero-length input SFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( fMin >= 0, XLAL_EDOM, "Invalid negative frequency fMin = %g\n", fMin );
  XLAL_CHECK_NULL ( Band > 0, XLAL_EDOM, "Invalid non-positive Band = %g\n", Band );

  UINT4 numSFTs = inSFTs->length;

  SFTVector *ret;
  XLAL_CHECK_NULL ( (ret = XLALCreateSFTVector ( numSFTs, 0 )) != NULL, XLAL_EFUNC );

  for ( UINT4 i = 0; i < numSFTs; i ++ )
    {
      SFTtype *dest = &(ret->data[i]);
      SFTtype *src =  &(inSFTs->data[i]);

      XLAL_CHECK_NULL ( XLALExtractBandFromSFT ( &dest, src, fMin, Band ) == XLAL_SUCCESS, XLAL_EFUNC );

    } /* for i < numSFTs */

  /* return final SFT-vector */
  return ret;

} /* XLALExtractBandFromSFTVector() */


/**
 * Return a MultiSFT vector containing only the bins in [fMin, fMin+Band].
 * Note: the output MultiSFT is guaranteed to "cover" the input boundaries 'fMin'
 * and 'fMin+Band', ie if necessary the output SFT contains one additional
 * bin on either end of the interval.
 *
 * This uses the conventions in XLALFindCoveringSFTBins() to determine
 * the 'effective' frequency-band to extract.
 *
 * \warning This convention is deprecated. Please use either
 * XLALExtractStrictBandFromMultiSFTVector(), or else XLALMultiSFTVectorResizeBand()
 * if you really need a covering frequency band.
 */
MultiSFTVector *
XLALExtractBandFromMultiSFTVector ( const MultiSFTVector *inSFTs,      ///< [in] input MultiSFTs
                                    REAL8 fMin,                                ///< [in] lower end of frequency interval to return
                                    REAL8 Band                         ///< [in] band width of frequency interval to return
                                    )
{
  XLAL_PRINT_DEPRECATION_WARNING("XLALExtractStrictBandFromMultiSFTVector");

  XLAL_CHECK_NULL ( inSFTs != NULL, XLAL_EINVAL, "Invalid NULL input MultiSFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( inSFTs->length > 0, XLAL_EINVAL, "Invalid zero-length input MultiSFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( fMin >= 0, XLAL_EDOM, "Invalid negative frequency fMin = %g\n", fMin );
  XLAL_CHECK_NULL ( Band > 0, XLAL_EDOM, "Invalid non-positive Band = %g\n", Band );

  MultiSFTVector *ret = NULL;
  XLAL_CHECK_NULL ( (ret = XLALCalloc(1, sizeof(*ret))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( (ret->data = XLALCalloc(inSFTs->length, sizeof(ret->data[0]))) != NULL, XLAL_ENOMEM );
  ret->length = inSFTs->length;

  for (UINT4 X = 0; X < inSFTs->length; X++) {
     XLAL_CHECK_NULL( (ret->data[X] = XLALExtractBandFromSFTVector(inSFTs->data[X], fMin, Band)) != NULL, XLAL_EFUNC );
  }

  return ret;
} //XLALExtractBandFromMultiSFTVector()


/**
 * Return a copy of an SFT containing only the bins in [fMin, fMin+Band).
 */
int
XLALExtractStrictBandFromSFT ( SFTtype **outSFT,	///< [out] output SFT (alloc'ed or re-alloced as required)
                               const SFTtype *inSFT,	///< [in] input SFT
                               REAL8 fMin,		///< [in] lower end of frequency interval to return
                               REAL8 Band		///< [in] band width of frequency interval to return
  )
{
  XLAL_CHECK ( outSFT != NULL, XLAL_EINVAL );
  XLAL_CHECK ( inSFT != NULL, XLAL_EINVAL );
  XLAL_CHECK ( inSFT->data != NULL, XLAL_EINVAL );
  XLAL_CHECK ( fMin >= 0, XLAL_EDOM, "Invalid negative frequency fMin = %g\n", fMin );
  XLAL_CHECK ( Band > 0, XLAL_EDOM, "Invalid non-positive Band = %g\n", Band );

  REAL8 df = inSFT->deltaF;

  REAL8 fMinSFT    = inSFT->f0;
  UINT4 numBinsSFT = inSFT->data->length;
  UINT4 firstBinSFT= round ( fMinSFT / df );	// round to closest bin
  UINT4 lastBinSFT = firstBinSFT + ( numBinsSFT - 1 );

  UINT4 firstBinExt = XLALRoundFrequencyDownToSFTBin( fMin, df );
  UINT4 lastBinExt = XLALRoundFrequencyUpToSFTBin( fMin + Band, df ) - 1;
  UINT4 numBinsExt = lastBinExt - firstBinExt + 1;

  XLAL_CHECK ( firstBinExt >= firstBinSFT && (lastBinExt <= lastBinSFT), XLAL_EINVAL,
               "Requested frequency-bins [%f,%f)Hz = [%d, %d] not contained within SFT's [%f, %f)Hz = [%d,%d].\n",
               fMin, fMin + Band, firstBinExt, lastBinExt, fMinSFT, fMinSFT + (numBinsSFT-1) * df, firstBinSFT, lastBinSFT );

  INT4 firstBinOffset = firstBinExt - firstBinSFT;

  if ( (*outSFT) == NULL ) {
    XLAL_CHECK ( ((*outSFT) = XLALCalloc(1, sizeof(*(*outSFT)))) != NULL, XLAL_ENOMEM );
  }
  if ( (*outSFT)->data == NULL ) {
    XLAL_CHECK ( ((*outSFT)->data = XLALCreateCOMPLEX8Vector ( numBinsExt )) != NULL, XLAL_EFUNC );
  }
  if ( (*outSFT)->data->length != numBinsExt ) {
    XLAL_CHECK ( ((*outSFT)->data->data = XLALRealloc ( (*outSFT)->data->data, numBinsExt * sizeof((*outSFT)->data->data[0]))) != NULL, XLAL_ENOMEM );
    (*outSFT)->data->length = numBinsExt;
  }

  COMPLEX8Vector *ptr = (*outSFT)->data;	// keep copy to data-pointer
  (*(*outSFT)) = (*inSFT);			// copy complete header
  (*outSFT)->data = ptr;	  		// restore data-pointer
  (*outSFT)->f0 = firstBinExt * df ;	  	// set correct new fMin

  /* copy the relevant part of the data */
  memcpy ( (*outSFT)->data->data, inSFT->data->data + firstBinOffset, numBinsExt * sizeof( (*outSFT)->data->data[0] ) );

  return XLAL_SUCCESS;

} // XLALExtractStrictBandFromSFT()


/**
 * Return a copy of a vector of SFTs containing only the bins in [fMin, fMin+Band).
 */
SFTVector *
XLALExtractStrictBandFromSFTVector ( const SFTVector *inSFTs,	///< [in] input SFTs
                                     REAL8 fMin,		///< [in] lower end of frequency interval to return
                                     REAL8 Band		///< [in] band width of frequency interval to return
  )
{
  XLAL_CHECK_NULL ( inSFTs != NULL, XLAL_EINVAL, "Invalid NULL input SFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( inSFTs->length > 0, XLAL_EINVAL, "Invalid zero-length input SFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( fMin >= 0, XLAL_EDOM, "Invalid negative frequency fMin = %g\n", fMin );
  XLAL_CHECK_NULL ( Band > 0, XLAL_EDOM, "Invalid non-positive Band = %g\n", Band );

  UINT4 numSFTs = inSFTs->length;

  SFTVector *ret;
  XLAL_CHECK_NULL ( (ret = XLALCreateSFTVector ( numSFTs, 0 )) != NULL, XLAL_EFUNC );

  for ( UINT4 i = 0; i < numSFTs; i ++ )
    {
      SFTtype *dest = &(ret->data[i]);
      SFTtype *src =  &(inSFTs->data[i]);

      XLAL_CHECK_NULL ( XLALExtractStrictBandFromSFT ( &dest, src, fMin, Band ) == XLAL_SUCCESS, XLAL_EFUNC );

    } /* for i < numSFTs */

  /* return final SFT-vector */
  return ret;

} /* XLALExtractStrictBandFromSFTVector() */


/**
 * Return a copy of a MultiSFT vector containing only the bins in [fMin, fMin+Band).
 */
MultiSFTVector *
XLALExtractStrictBandFromMultiSFTVector ( const MultiSFTVector *inSFTs,      ///< [in] input MultiSFTs
                                          REAL8 fMin,                        ///< [in] lower end of frequency interval to return
                                          REAL8 Band                         ///< [in] band width of frequency interval to return
  )
{
  XLAL_CHECK_NULL ( inSFTs != NULL, XLAL_EINVAL, "Invalid NULL input MultiSFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( inSFTs->length > 0, XLAL_EINVAL, "Invalid zero-length input MultiSFT vector 'inSFTs'\n");
  XLAL_CHECK_NULL ( fMin >= 0, XLAL_EDOM, "Invalid negative frequency fMin = %g\n", fMin );
  XLAL_CHECK_NULL ( Band > 0, XLAL_EDOM, "Invalid non-positive Band = %g\n", Band );

  MultiSFTVector *ret = NULL;
  XLAL_CHECK_NULL ( (ret = XLALCalloc(1, sizeof(*ret))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL ( (ret->data = XLALCalloc(inSFTs->length, sizeof(ret->data[0]))) != NULL, XLAL_ENOMEM );
  ret->length = inSFTs->length;

  for (UINT4 X = 0; X < inSFTs->length; X++) {
     XLAL_CHECK_NULL( (ret->data[X] = XLALExtractStrictBandFromSFTVector(inSFTs->data[X], fMin, Band)) != NULL, XLAL_EFUNC );
  }

  return ret;
} //XLALExtractStrictBandFromMultiSFTVector()


/** Append the given SFTtype to the SFT-vector (no SFT-specific checks are done!) */
int XLALAppendSFT2Vector (SFTVector *vect,		/**< destinatino SFTVector to append to */
                          const SFTtype *sft            /**< the SFT to append */
                          )
{
  UINT4 oldlen = vect->length;

  if ( (vect->data = LALRealloc ( vect->data, (oldlen + 1)*sizeof( *vect->data ) )) == NULL ) {
     XLAL_ERROR(XLAL_ENOMEM);
  }
  memset ( &(vect->data[oldlen]), 0, sizeof( vect->data[0] ) );
  vect->length ++;

  XLALCopySFT(&vect->data[oldlen], sft );

  return XLAL_SUCCESS;

} /* XLALAppendSFT2Vector() */


/**
 * Reorder the MultiSFTVector with specified list of IFOs
 */
int XLALReorderMultiSFTVector( MultiSFTVector *multiSFTs, const LALStringVector *IFOs)
{
  XLAL_CHECK( multiSFTs!=NULL && IFOs!=NULL && multiSFTs->length==IFOs->length && multiSFTs->length<=PULSAR_MAX_DETECTORS, XLAL_EINVAL );

  // Initialize array of reordered SFTVector pointers
  SFTVector *reordered[PULSAR_MAX_DETECTORS];
  XLAL_INIT_MEM(reordered);

  // Loop through IFO list and reorder if necessary
  for (UINT4 i=0; i < IFOs->length; i ++ )
    {
      UINT4 j=0;
      while ( (j < IFOs->length) && (strncmp ( IFOs->data[i], multiSFTs->data[j]->data[0].name, 2 ) != 0) ) {
        j++;
      }
      XLAL_CHECK ( j < IFOs->length, XLAL_EINVAL, "IFO %c%c not found", IFOs->data[i][0], IFOs->data[i][1] );
      reordered[i] = multiSFTs->data[j]; // copy the SFTVector pointer
    }

  // Replace the old pointers with the new values
  for ( UINT4 i=0; i < multiSFTs->length; i ++ )
    {
      multiSFTs->data[i] = reordered[i];
    }

  return XLAL_SUCCESS;

} // XLALReorderMultiSFTVector()


/**
 * Adds SFT-data from SFT 'b' to SFT 'a'
 *
 * NOTE: the inputs 'a' and 'b' must have consistent
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

  XLAL_CHECK ( XLALGPSDiff ( &(a->epoch), &(b->epoch) ) == 0, XLAL_EINVAL, "SFT epochs differ %"LAL_INT8_FORMAT" != %"LAL_INT8_FORMAT" ns\n", XLALGPSToINT8NS ( &(a->epoch) ), XLALGPSToINT8NS ( &(b->epoch) ) );

  REAL8 tol = 10 * LAL_REAL8_EPS;	// generously allow up to 10*eps tolerance
  XLAL_CHECK ( gsl_fcmp ( a->f0, b->f0, tol ) == 0, XLAL_ETOL, "SFT frequencies relative deviation exceeds %g: %.16g != %.16g\n", tol, a->f0, b->f0 );
  XLAL_CHECK ( gsl_fcmp ( a->deltaF, b->deltaF, tol ) == 0, XLAL_ETOL, "SFT frequency-steps relative deviation exceeds %g: %.16g != %.16g\n", tol, a->deltaF, b->deltaF );
  XLAL_CHECK ( XLALUnitCompare ( &(a->sampleUnits), &(b->sampleUnits) ) == 0, XLAL_EINVAL, "SFT sample units differ\n" );
  XLAL_CHECK ( a->data->length == b->data->length, XLAL_EINVAL, "SFT lengths differ: %d != %d\n", a->data->length, b->data->length );

  UINT4 numBins = a->data->length;
  for ( UINT4 k = 0; k < numBins; k ++ )
    {
      a->data->data[k] += b->data->data[k];
    }

  return XLAL_SUCCESS;

} /* XLALSFTAdd() */


/**
 * Adds SFT-data from SFTvector 'b' to elements of SFTVector 'a'
 *
 * NOTE: the inputs 'a' and 'b' must have consistent number of SFTs,
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

      XLAL_CHECK ( strncmp ( vect1->data[0].name, vect2->data[0].name, 2 ) == 0, XLAL_EINVAL, "SFT detectors differ '%c%c' != '%c%c'\n", vect1->data[0].name[0], vect1->data[0].name[1], vect2->data[0].name[0], vect2->data[0].name[1] );

      XLAL_CHECK ( XLALSFTVectorAdd ( vect1, vect2 ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALSFTVectorAdd() failed for SFTVector %d out of %d\n", X+1, numIFOs );

    } // for X < numIFOs

  return XLAL_SUCCESS;

} /* XLALMultiSFTVectorAdd() */


/**
 * Resize the frequency-band of a given SFT to [f0, f0+Band].
 *
 * NOTE: If the frequency band is extended in any direction, the corresponding bins
 * will be set to zero
 *
 * NOTE2: This uses the conventions in XLALFindCoveringSFTBins() to determine
 * the 'effective' frequency-band to resize to, in order to coincide with SFT frequency bins.
 *
 */
int
XLALSFTResizeBand ( SFTtype *SFT,	///< [in/out] SFT to resize
                    REAL8 f0,		///< [in] new start frequency
                    REAL8 Band		///< [in] new frequency Band
                    )
{
  XLAL_CHECK ( SFT != NULL, XLAL_EINVAL );
  XLAL_CHECK ( f0 >= 0, XLAL_EINVAL );
  XLAL_CHECK ( Band >= 0, XLAL_EINVAL );


  REAL8 Tsft = TSFTfromDFreq ( SFT->deltaF );
  REAL8 f0In = SFT->f0;

  UINT4 firstBinIn = (UINT4) lround ( f0In / SFT->deltaF );

  UINT4 firstBinOut;
  UINT4 numBinsOut;
  XLAL_CHECK ( XLALFindCoveringSFTBins ( &firstBinOut, &numBinsOut, f0, Band, Tsft ) == XLAL_SUCCESS, XLAL_EFUNC );

  int firstRelative = firstBinOut - firstBinIn;

  XLAL_CHECK ( (SFT = XLALResizeCOMPLEX8FrequencySeries ( SFT, firstRelative, numBinsOut )) != NULL, XLAL_EFUNC );

  return XLAL_SUCCESS;

} // XLALSFTResizeBand()


/**
 * Resize the frequency-band of a given SFT vector to [f0, f0+Band].
 *
 * NOTE: If the frequency band is extended in any direction, the corresponding bins
 * will be set to zero
 *
 * NOTE2: This uses the conventions in XLALFindCoveringSFTBins() to determine
 * the 'effective' frequency-band to resize to, in order to coincide with SFT frequency bins.
 *
 */
int
XLALSFTVectorResizeBand ( SFTVector *SFTs,	///< [in/out] SFT vector to resize
                          REAL8 f0,		///< [in] new start frequency
                          REAL8 Band		///< [in] new frequency Band
                          )
{
  XLAL_CHECK ( SFTs != NULL, XLAL_EINVAL );
  XLAL_CHECK ( f0 >= 0, XLAL_EINVAL );
  XLAL_CHECK ( Band >= 0, XLAL_EINVAL );

  for ( UINT4 alpha = 0; alpha < SFTs->length; alpha ++ ) {
    XLAL_CHECK ( XLALSFTResizeBand ( &(SFTs->data[alpha]), f0, Band ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} // XLALSFTVectorResizeBand()


/**
 * Resize the frequency-band of a given multi-SFT vector to [f0, f0+Band].
 *
 * NOTE: If the frequency band is extended in any direction, the corresponding bins
 * will be set to zero
 *
 * NOTE2: This uses the conventions in XLALFindCoveringSFTBins() to determine
 * the 'effective' frequency-band to resize to, in order to coincide with SFT frequency bins.
 *
 */
int
XLALMultiSFTVectorResizeBand ( MultiSFTVector *multiSFTs,	///< [in/out] multi-SFT vector to resize
                               REAL8 f0,			///< [in] new start frequency
                               REAL8 Band			///< [in] new frequency Band
                               )
{
  XLAL_CHECK ( multiSFTs != NULL, XLAL_EINVAL );
  XLAL_CHECK ( f0 >= 0, XLAL_EINVAL );
  XLAL_CHECK ( Band >= 0, XLAL_EINVAL );

  for ( UINT4 X = 0; X < multiSFTs->length; X ++ ) {
    XLAL_CHECK ( XLALSFTVectorResizeBand ( multiSFTs->data[X], f0, Band ) == XLAL_SUCCESS, XLAL_EFUNC );
  }

  return XLAL_SUCCESS;

} // XLALMultiSFTVectorResizeBand()


/**
 * Finds the earliest timestamp in a multi-SFT data structure
 */
int
XLALEarliestMultiSFTsample ( LIGOTimeGPS *out,              /**< [out] earliest GPS time */
                             const MultiSFTVector *multisfts      /**< [in] multi SFT vector */
                             )
{
  UINT4 i,j;

  /* check sanity of input */
  if ( !multisfts || (multisfts->length == 0) )
    {
      XLALPrintError ("%s: empty multiSFT input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  for (i=0;i<multisfts->length;i++)
    {
      if ( !multisfts->data[i] || (multisfts->data[i]->length == 0) )
        {
          XLALPrintError ("%s: empty multiSFT->data[%d] input!\n", __func__,i );
          XLAL_ERROR (XLAL_EINVAL);
        }
    }

  /* initialise the earliest sample value */
  out->gpsSeconds = multisfts->data[0]->data[0].epoch.gpsSeconds;
  out->gpsNanoSeconds = multisfts->data[0]->data[0].epoch.gpsNanoSeconds;

  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {

    /* loop over all SFTs to determine the earliest SFT epoch */
    for (j=0;j<multisfts->data[i]->length;j++) {

      /* compare current SFT epoch with current earliest */
      if ( (XLALGPSCmp(out,&multisfts->data[i]->data[j].epoch) == 1 ) ) {
        out->gpsSeconds = multisfts->data[i]->data[j].epoch.gpsSeconds;
        out->gpsNanoSeconds = multisfts->data[i]->data[j].epoch.gpsNanoSeconds;
      }

    }

  }

  /* success */
  return XLAL_SUCCESS;

} /* XLALEarliestMultiSFTsample() */


/**
 * Find the time of the end of the latest SFT in a multi-SFT data structure
 */
int
XLALLatestMultiSFTsample ( LIGOTimeGPS *out,              /**< [out] latest GPS time */
                           const MultiSFTVector *multisfts      /**< [in] multi SFT vector */
                           )
{
  UINT4 i,j;
  SFTtype *firstSFT;

  /* check sanity of input */
  if ( !multisfts || (multisfts->length == 0) )
    {
      XLALPrintError ("%s: empty multiSFT input!\n", __func__ );
      XLAL_ERROR (XLAL_EINVAL);
    }
  for (i=0;i<multisfts->length;i++)
    {
      if ( !multisfts->data[i] || (multisfts->data[i]->length == 0) )
        {
          XLALPrintError ("%s: empty multiSFT->data[%d] input!\n", __func__,i );
          XLAL_ERROR (XLAL_EINVAL);
        }
    }

  /* define some useful quantities */
  firstSFT = (multisfts->data[0]->data);        /* a pointer to the first SFT of the first detector */
  REAL8 Tsft = TSFTfromDFreq ( firstSFT->deltaF );    /* the length of the SFTs in seconds assuming 1/T freq resolution */

  /* initialise the latest sample value */
  out->gpsSeconds = firstSFT->epoch.gpsSeconds;
  out->gpsNanoSeconds = firstSFT->epoch.gpsNanoSeconds;

  /* loop over detectors */
  for (i=0;i<multisfts->length;i++) {

    /* loop over all SFTs to determine the earliest SFT midpoint of the input data in the SSB frame */
    for (j=0;j<multisfts->data[i]->length;j++) {

      /* compare current SFT epoch with current earliest */
      if ( (XLALGPSCmp(out,&multisfts->data[i]->data[j].epoch) == -1 ) ) {
        out->gpsSeconds = multisfts->data[i]->data[j].epoch.gpsSeconds;
        out->gpsNanoSeconds = multisfts->data[i]->data[j].epoch.gpsNanoSeconds;
      }
    }

  }

  /* add length of SFT to the result so that we output the end of the SFT */
  if ( XLALGPSAdd(out,Tsft) == NULL )
    {
      XLALPrintError ("%s: NULL pointer returned from XLALGPSAdd()!\n", __func__ );
      XLAL_ERROR (XLAL_EFAULT);
    }

  /* success */
  return XLAL_SUCCESS;

} /* XLALLatestMultiSFTsample() */


/**
 * Extract an SFTVector from another SFTVector but only those timestamps matching
 *
 * Timestamps must be a subset of those sfts in the SFTVector or an error occurs
 */
SFTVector *
XLALExtractSFTVectorWithTimestamps ( const SFTVector *sfts,                 /**< input SFTs */
                                     const LIGOTimeGPSVector *timestamps    /**< timestamps */
                                     )
{
  // check input sanity
  XLAL_CHECK_NULL( sfts != NULL, XLAL_EINVAL);
  XLAL_CHECK_NULL( timestamps != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL( sfts->length >= timestamps->length, XLAL_EINVAL );

  SFTVector *ret = NULL;
  XLAL_CHECK_NULL( (ret = XLALCreateSFTVector(timestamps->length, 0)) != NULL, XLAL_EFUNC );

  UINT4 indexOfInputSFTVector = 0;
  UINT4 numberOfSFTsLoadedIntoOutputVector = 0;
  for (UINT4 ii=0; ii<timestamps->length; ii++)
    {
      XLAL_CHECK_NULL( indexOfInputSFTVector < sfts->length, XLAL_FAILURE, "At least one timestamp is not in the range specified by the SFT vector" );

      for (UINT4 jj=indexOfInputSFTVector; jj<sfts->length; jj++)
        {
        if ( XLALGPSCmp(&(sfts->data[jj].epoch), &(timestamps->data[ii])) == 0 )
          {
            indexOfInputSFTVector = jj+1;
            XLAL_CHECK_NULL( XLALCopySFT(&(ret->data[ii]), &(sfts->data[jj])) == XLAL_SUCCESS, XLAL_EFUNC );
            numberOfSFTsLoadedIntoOutputVector++;
            break;
          } // if SFT epoch matches timestamp epoch
        } // for jj < sfts->length
    } // for ii < timestamps->length

  XLAL_CHECK_NULL( numberOfSFTsLoadedIntoOutputVector == ret->length, XLAL_FAILURE, "Not all timestamps were found in the input SFT vector" );

  return ret;

} // XLALExtractSFTVectorWithTimestamps


/**
 * Extract a MultiSFTVector from another MultiSFTVector but only those timestamps matching
 *
 * Timestamps in each LIGOTimeGPSVector must be a subset of those sfts in each SFTVector or an error occurs
 */
MultiSFTVector *
XLALExtractMultiSFTVectorWithMultiTimestamps ( const MultiSFTVector *multiSFTs,                 /**< input SFTs */
                                               const MultiLIGOTimeGPSVector *multiTimestamps    /**< timestamps */
                                               )
{
  // check input sanity
  XLAL_CHECK_NULL( multiSFTs != NULL, XLAL_EINVAL);
  XLAL_CHECK_NULL( multiTimestamps != NULL, XLAL_EINVAL );

  MultiSFTVector *ret = NULL;
  XLAL_CHECK_NULL( (ret = XLALCalloc(1, sizeof(*ret))) != NULL, XLAL_ENOMEM );
  XLAL_CHECK_NULL( (ret->data = XLALCalloc(multiSFTs->length, sizeof(*ret->data))) != NULL, XLAL_ENOMEM );
  ret->length = multiSFTs->length;

  for (UINT4 X=0; X<multiSFTs->length; X++)
    {
       XLAL_CHECK_NULL( (ret->data[X] = XLALExtractSFTVectorWithTimestamps(multiSFTs->data[X], multiTimestamps->data[X])) != NULL, XLAL_EFUNC );
    }

  return ret;

} // XLALExtractMultiSFTVectorWithMultiTimestamps


// Function to 'safely' invert Tsft=1/dFreq to avoid roundoff error for close-but-not-quite integers after inversion
// policy: if (1/dFreq) is within 10 * eps of an integer, round, otherwise leave as a fractional number
// comment: unfortunately Tsft is allowed by spec to be a double, but really should be limited to integer seconds,
// however with this function we should be able to safely work with integer-valued Tsft without leaving spec (yet)
REAL8
TSFTfromDFreq ( REAL8 dFreq )
{
  REAL8 Tsft0 = 1.0 / dFreq;
  REAL8 Tsft;
  if ( fabs ( (Tsft0 - round(Tsft0)) ) / Tsft0 < 10 * LAL_REAL8_EPS ) {
    Tsft = round(Tsft0);
  } else {
    Tsft = Tsft0;
  }

  return Tsft;

} // TSFTfromDFreq()


/* compare two SFT-descriptors by their GPS-epoch, then starting frequency */
int
compareSFTdesc(const void *ptr1, const void *ptr2)
{
  const SFTDescriptor *desc1 = ptr1;
  const SFTDescriptor *desc2 = ptr2;

  if      ( GPS2REAL8( desc1->header.epoch ) < GPS2REAL8 ( desc2->header.epoch ) )
    return -1;
  else if ( GPS2REAL8( desc1->header.epoch ) > GPS2REAL8 ( desc2->header.epoch ) )
    return 1;
  else if ( desc1->header.f0 < desc2->header.f0 )
    return -1;
  else if ( desc1->header.f0 > desc2->header.f0 )
    return 1;
  else
    return 0;
} /* compareSFTdesc() */


/* compare two SFT-descriptors by their locator (f0, file, position) */
int
compareSFTloc(const void *ptr1, const void *ptr2)
{
  const SFTDescriptor *desc1 = ptr1;
  const SFTDescriptor *desc2 = ptr2;
  int s;
  if ( desc1->header.f0 < desc2->header.f0 )
    return -1;
  else if ( desc1->header.f0 > desc2->header.f0 )
    return 1;
  s = strcmp(desc1->locator->fname, desc2->locator->fname);
  if(!s) {
    if (desc1->locator->offset < desc2->locator->offset)
      return(-1);
    else if (desc1->locator->offset > desc2->locator->offset)
      return(1);
    else
      return(0);
  }
  return(s);
} /* compareSFTloc() */


/* compare two SFT-catalog by detector name in alphabetic order */
int
compareDetNameCatalogs ( const void *ptr1, const void *ptr2 )
{
  SFTCatalog const* cat1 = (SFTCatalog const*)ptr1;
  SFTCatalog const* cat2 = (SFTCatalog const*)ptr2;
  const char *name1 = cat1->data[0].header.name;
  const char *name2 = cat2->data[0].header.name;

  if ( name1[0] < name2[0] )
    return -1;
  else if ( name1[0] > name2[0] )
    return 1;
  else if ( name1[1] < name2[1] )
    return -1;
  else if ( name1[1] > name2[1] )
    return 1;
  else
    return 0;

} /* compareDetNameCatalogs() */


/* compare two SFT-descriptors by their GPS-epoch, then starting frequency */
int compareSFTepoch(const void *ptr1, const void *ptr2)
 {
   const SFTtype *desc1 = ptr1;
   const SFTtype *desc2 = ptr2;

   if      ( XLALGPSGetREAL8( &desc1->epoch ) < XLALGPSGetREAL8 ( &desc2->epoch ) )
     return -1;
   else if ( XLALGPSGetREAL8( &desc1->epoch ) > XLALGPSGetREAL8 ( &desc2->epoch ) )
     return 1;
   else if ( desc1->f0 < desc2->f0 )
     return -1;
   else if ( desc1->f0 > desc2->f0 )
     return 1;
   else
     return 0;
} /* compareSFTepoch() */

/// @}
