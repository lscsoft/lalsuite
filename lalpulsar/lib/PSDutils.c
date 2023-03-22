/*
 * Copyright (C) 2013, 2020 David Keitel
 * Copyright (C) 2010, 2013, 2014, 2016, 2017, 2022 Karl Wette
 * Copyright (C) 2010 Chris Messenger
 * Copyright (C) 2007, 2013 John T. Whelan
 * Copyright (C) 2006 Badri Krishnan
 * Copyright (C) 2004--2006, 2008--2013, 2016 Reinhard Prix
 * Copyright (C) 2004, 2005 Bernd Machenschalk
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
#include <gsl/gsl_sort_double.h>

#include <lal/Units.h>
#include <lal/FrequencySeries.h>
#include <lal/NormalizeSFTRngMed.h>

#include <lal/PSDutils.h>
#include "SFTinternal.h"

/*---------- global variables ----------*/

// for backwards compatibility, we currently allow calling the same internal
// enum entries by either human-friendly names or by the old-style numerical
// arguments. The latter may get deprecated at some point in the future.
const UserChoices MathOpTypeChoices = {
  { MATH_OP_ARITHMETIC_SUM,    "arithsum" },
  { MATH_OP_ARITHMETIC_MEAN,   "arithmean" },
  { MATH_OP_ARITHMETIC_MEDIAN, "arithmedian" },
  { MATH_OP_HARMONIC_SUM,      "harmsum" },
  { MATH_OP_HARMONIC_MEAN,     "harmmean" },
  { MATH_OP_POWERMINUS2_SUM,   "powerminus2sum" },
  { MATH_OP_POWERMINUS2_MEAN,  "powerminus2mean" },
  { MATH_OP_MINIMUM,           "min" },
  { MATH_OP_MAXIMUM,           "max" },
  { MATH_OP_ARITHMETIC_SUM,    "0" },
  { MATH_OP_ARITHMETIC_MEAN,   "1" },
  { MATH_OP_ARITHMETIC_MEDIAN, "2" },
  { MATH_OP_HARMONIC_SUM,      "3" },
  { MATH_OP_HARMONIC_MEAN,     "4" },
  { MATH_OP_POWERMINUS2_SUM,   "5" },
  { MATH_OP_POWERMINUS2_MEAN,  "6" },
  { MATH_OP_MINIMUM,           "7" },
  { MATH_OP_MAXIMUM,           "8" },
};

/*========== function definitions ==========*/

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

///
/// Create a \c MultiNoiseWeights of the given length.
///
/// NOTE: this does not allocate the individual data entries.
///
MultiNoiseWeights*
XLALCreateMultiNoiseWeights ( const UINT4 length )  ///< [in] Length of the \c MultiNoiseWeights.
{

  MultiNoiseWeights *multiWeights = NULL;
  XLAL_CHECK_NULL ( (multiWeights = XLALCalloc(1, sizeof(*multiWeights))) != NULL, XLAL_ENOMEM );
  multiWeights->length = length;
  if ( length > 0 ) {
    XLAL_CHECK_NULL ( (multiWeights->data = XLALCalloc ( length, sizeof(*multiWeights->data))) != NULL, XLAL_ENOMEM );
  }

  return multiWeights;

} // XLALCreateMultiNoiseWeights()

///
/// Create a full copy of a \c MultiNoiseWeights object.
///
MultiNoiseWeights*
XLALCopyMultiNoiseWeights ( const MultiNoiseWeights *multiWeights )  ///< [in] Length of the \c MultiNoiseWeights.
{

  // check input sanity
  XLAL_CHECK_NULL ( multiWeights != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( multiWeights->data != NULL, XLAL_EINVAL );
  XLAL_CHECK_NULL ( multiWeights->length > 0, XLAL_EINVAL );

  MultiNoiseWeights *copiedWeights = NULL;
  XLAL_CHECK_NULL ( ( copiedWeights = XLALCreateMultiNoiseWeights ( multiWeights->length ) ) != NULL, XLAL_EFUNC );
  copiedWeights->length = multiWeights->length;
  copiedWeights->Sinv_Tsft = multiWeights->Sinv_Tsft;
  copiedWeights->isNotNormalized = multiWeights->isNotNormalized;

  for ( UINT4 X = 0; X < multiWeights->length; X++) {
    /* create k^th weights vector */
    XLAL_CHECK_NULL ( ( copiedWeights->data[X] = XLALCreateREAL8Vector ( multiWeights->data[X]->length ) ) != NULL, XLAL_EFUNC );
    for ( UINT4 alpha = 0; alpha < multiWeights->data[X]->length; alpha++) {
        copiedWeights->data[X]->data[alpha] = multiWeights->data[X]->data[alpha];
    }
  }

  return copiedWeights;

} // XLALCopyMultiNoiseWeights()


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
 * Function that *truncates the PSD in place* to the requested frequency-bin interval [firstBin, lastBin] for the given multiPSDVector.
 * Now also truncates the original SFT vector, as necessary for correct computation of normalized SFT power.
 */
int
XLALCropMultiPSDandSFTVectors ( MultiPSDVector *multiPSDVect,
                         MultiSFTVector *multiSFTVect,
                         UINT4 firstBin,
                         UINT4 lastBin
                         )
{
  /* check user input */
  if ( !multiPSDVect ) {
    XLALPrintError ("%s: invalid NULL input 'multiPSDVect'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( !multiSFTVect ) {
    XLALPrintError ("%s: invalid NULL input 'multiSFTVect'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( lastBin < firstBin ) {
    XLALPrintError ("%s: empty bin interval requested [%d, %d]\n", __func__, firstBin, lastBin );
    XLAL_ERROR ( XLAL_EDOM );
  }

  UINT4 numIFOs = multiPSDVect->length;
  UINT4 numBins = multiPSDVect->data[0]->data[0].data->length;

  if ( numIFOs != multiSFTVect->length ) {
    XLALPrintError ("%s: inconsistent number of IFOs between PSD (%d) and SFT (%d) vectors.\n", __func__, numIFOs, multiSFTVect->length );
    XLAL_ERROR ( XLAL_EDOM );
  }

  if ( numBins != multiSFTVect->data[0]->data[0].data->length ) {
    XLALPrintError ("%s: inconsistent number of bins between PSD (%d bins) and SFT (%d bins) vectors.\n", __func__, numBins, multiSFTVect->data[0]->data[0].data->length );
    XLAL_ERROR ( XLAL_EDOM );
  }

  if ( (firstBin >= numBins) || (lastBin >= numBins ) ) {
    XLALPrintError ("%s: requested bin-interval [%d, %d] outside of PSD bins [0, %d]\n", __func__, firstBin, lastBin, numBins - 1 );
    XLAL_ERROR ( XLAL_EDOM );
  }

  /* ----- check if there's anything to do at all? ----- */
  if ( (firstBin == 0)  && (lastBin == numBins - 1) )
    return XLAL_SUCCESS;

  /* ----- loop over detectors, timestamps, then crop each PSD ----- */
  UINT4 X;
  for ( X=0; X < numIFOs; X ++ )
    {
      PSDVector *thisPSDVect = multiPSDVect->data[X];
      SFTVector *thisSFTVect = multiSFTVect->data[X];
      UINT4 numTS   = thisPSDVect->length;

      UINT4 iTS;
      for ( iTS = 0; iTS < numTS; iTS ++ )
        {
          REAL8FrequencySeries *thisPSD = &thisPSDVect->data[iTS];
          COMPLEX8FrequencySeries *thisSFT = &thisSFTVect->data[iTS];

          if ( numBins != thisPSD->data->length ) {
            XLALPrintError ("%s: inconsistent number of frequency-bins across multiPSDVector: X=%d, iTS=%d: numBins = %d != %d\n",
                            __func__, X, iTS, numBins, thisPSD->data->length );
            XLAL_ERROR ( XLAL_EDOM );
          }

          if ( numBins != thisSFT->data->length ) {
            XLALPrintError ("%s: inconsistent number of frequency-bins across multiSFTVector: X=%d, iTS=%d: numBins = %d != %d\n",
                            __func__, X, iTS, numBins, thisSFT->data->length );
            XLAL_ERROR ( XLAL_EDOM );
          }

          UINT4 numNewBins = lastBin - firstBin + 1;

          /* crop PSD */
          XLAL_CHECK(XLALResizeREAL8FrequencySeries(thisPSD, firstBin, numNewBins) != NULL, XLAL_EFUNC);

          /* crop SFT */
          XLAL_CHECK(XLALResizeCOMPLEX8FrequencySeries(thisSFT, firstBin, numNewBins) != NULL, XLAL_EFUNC);

        } /* for iTS < numTS */

    } /* for X < numIFOs */

  /* that should be all ... */
  return XLAL_SUCCESS;

} /* XLALCropMultiPSDandSFTVectors() */


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
  REAL8 Tsft = TSFTfromDFreq ( rngmed->data[0]->data[0].deltaF );

  /* create multi noise weights for output */
  MultiNoiseWeights *multiWeights = NULL;
  XLAL_CHECK_NULL ( ( multiWeights = XLALCreateMultiNoiseWeights ( numIFOs ) ) != NULL, XLAL_EFUNC );

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


/**
 * Compute the "data-quality factor" \f$ \mathcal{Q}(f) = \sum_X \frac{\epsilon_X}{\mathcal{S}_X(f)} \f$ over the given SFTs.
 * The input \a multiPSD is a pre-computed PSD map \f$ \mathcal{S}_{X,i}(f) \f$ , over IFOs \f$ X \f$ , SFTs \f$ i \f$ 
 * and frequency \f$ f \f$ .
 *
 * \return the output is a vector \f$ \mathcal{Q}(f) \f$ .
 *
 */
REAL8FrequencySeries *
XLALComputeSegmentDataQ ( const MultiPSDVector *multiPSDVect, 	/**< input PSD map over IFOs, SFTs, and frequencies */
                          LALSeg segment		  	/**< segment to compute Q for */
                          )
{
  /* check input consistency */
  if ( multiPSDVect == NULL ) {
    XLALPrintError ("%s: NULL input 'multiPSDVect'\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }
  if ( multiPSDVect->length == 0 || multiPSDVect->data==0 ) {
    XLALPrintError ("%s: invalid multiPSDVect input (length=0 or data=NULL)\n", __func__ );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  REAL8 Tseg = XLALGPSDiff ( &segment.end, &segment.start );
  if ( Tseg <= 0 ) {
    XLALPrintError ("%s: negative segment-duration '%g'\n", __func__, Tseg );
    XLAL_ERROR_NULL ( XLAL_EINVAL );
  }

  REAL8 Tsft 	 = 1.0 / multiPSDVect->data[0]->data[0].deltaF;
  REAL8 f0       = multiPSDVect->data[0]->data[0].f0;
  REAL8 dFreq    = multiPSDVect->data[0]->data[0].deltaF;
  UINT4 numFreqs = multiPSDVect->data[0]->data[0].data->length;

  REAL8FrequencySeries *Q, *SXinv;
  if ( (Q = XLALCreateREAL8FrequencySeries ( "Qfactor", &segment.start, f0, dFreq, &lalHertzUnit, numFreqs )) == NULL ) {
    XLALPrintError ("%s: Q = XLALCreateREAL8FrequencySeries(numFreqs=%d) failed with xlalErrno = %d\n", __func__, numFreqs, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  if ( (SXinv = XLALCreateREAL8FrequencySeries ( "SXinv", &segment.start, f0, dFreq, &lalHertzUnit, numFreqs )) == NULL ) {
    XLALPrintError ("%s: SXinv = XLALCreateREAL8FrequencySeries(numFreqs=%d) failed with xlalErrno = %d\n", __func__, numFreqs, xlalErrno );
    XLAL_ERROR_NULL ( XLAL_EFUNC );
  }
  /* initialize Q-array to zero, as we'll keep adding to it */
  memset ( Q->data->data, 0, Q->data->length * sizeof(Q->data->data[0]) );

  /* ----- loop over IFOs ----- */
  UINT4 numIFOs = multiPSDVect->length;
  UINT4 X;
  for ( X = 0; X < numIFOs; X ++ )
    {
      PSDVector *thisPSDVect = multiPSDVect->data[X];

      /* initialize SXinv-array to zero, as we'll keep adding to it */
      memset ( SXinv->data->data, 0, SXinv->data->length * sizeof(SXinv->data->data[0]) );
      UINT4 numSFTsInSeg = 0;	/* reset counter of SFTs within this segment */

      /* ----- loop over all timestamps ----- */
      /* find SFTs inside segment, count them and combine their PSDs */
      UINT4 numTS = thisPSDVect->length;
      UINT4 iTS;
      for ( iTS = 0; iTS < numTS; iTS++ )
        {
          REAL8FrequencySeries *thisPSD = &thisPSDVect->data[iTS];

          /* some internal consistency/paranoia checks */
          if ( ( f0 != thisPSD->f0) || ( dFreq != thisPSD->deltaF ) || (numFreqs != thisPSD->data->length ) ) {
            XLALPrintError ("%s: %d-th timestamp %f: inconsistent PSDVector: f0 = %g : %g,  dFreq = %g : %g, numFreqs = %d : %d \n",
                            __func__, iTS, XLALGPSGetREAL8( &thisPSD->epoch ), f0, thisPSD->f0, dFreq, thisPSD->deltaF, numFreqs, thisPSD->data->length );
            XLAL_ERROR_NULL ( XLAL_EDOM );
          }

          int cmp = XLALCWGPSinRange( thisPSD->epoch, &segment.start, &segment.end );

          if ( cmp < 0 )	/* SFT-end before segment => advance to the next one */
            continue;
          if ( cmp > 0 )	/* SFT-start past end of segment: ==> terminate loop */
            break;

          if ( cmp == 0 )	/* this SFT is inside segment */
            {
              numSFTsInSeg ++;
              /* add SXinv(f) += 1/SX_i(f) over all frequencies */
              UINT4 iFreq;
              for ( iFreq = 0; iFreq < numFreqs; iFreq++ )
                SXinv->data->data[iFreq] += 1.0 / thisPSD->data->data[iFreq] ;

            }	/* if SFT inside segment */

        } /* for iTS < numTS */

      /* compute duty-cycle eps_X = nX * Tsft / Tseg for this IFO */
      REAL8 duty_X = numSFTsInSeg * Tsft / Tseg;
      /* sanity check: eps in [0, 1]*/
      if ( (duty_X < 0) && (duty_X > 1 ) ) {
        XLALPrintError ("%s: something is WRONG: duty-cyle = %g not within [0,1]!\n", __func__, duty_X );
        XLAL_ERROR_NULL ( XLAL_EFAILED );
      }

      /* add duty_X-weighted SXinv to Q */
      UINT4 iFreq;
      for ( iFreq = 0; iFreq < numFreqs; iFreq ++ )
        Q->data->data[iFreq] += duty_X * SXinv->data->data[iFreq] / numSFTsInSeg;

    } /* for X < numIFOs */

  /* clean up, free memory */
  XLALDestroyREAL8FrequencySeries ( SXinv );


  return Q;

} /* XLALComputeSegmentDataQ() */


/**
 * Compute various types of "math operations" over the entries of an array.
 *
 * The supported optypes (e.g. sums and averages) are defined in \c MathOpType.
 *
 * This can be used e.g. for the different established conventions of combining
 * SFTs for a PSD estimate.
 *
 */
REAL8 XLALMathOpOverArray ( const REAL8* data,      /**< input data array */
                            const size_t length,    /**< length of the input data array */
                            const MathOpType optype /**< type of operation */
                          ) {

  UINT4 i;
  REAL8 res = 0.0;

  switch (optype) {

  case MATH_OP_ARITHMETIC_SUM: /* sum(data) */

    for (i = 0; i < length; ++i) res += *(data++);

    break;

  case MATH_OP_ARITHMETIC_MEAN: /* sum(data)/length  */

    for (i = 0; i < length; ++i) res += *(data++);
    res /= (REAL8)length;

    break;

  case MATH_OP_HARMONIC_SUM: /* 1 / sum(1 / data) */

    for (i = 0; i < length; ++i) res += 1.0 / *(data++);
    res = 1.0 / res;

    break;

  case MATH_OP_HARMONIC_MEAN: /* length / sum(1 / data) */

    for (i = 0; i < length; ++i) res += 1.0 / *(data++);
    res = (REAL8)length / res;

    break;

  case MATH_OP_POWERMINUS2_SUM: /*   1 / sqrt ( sum(1 / data/data) )*/

    for (i = 0; i < length; ++i) res += 1.0 / (data[i]*data[i]);
    res = 1.0 / sqrt(res);

    break;

   case MATH_OP_POWERMINUS2_MEAN: /*   1 / sqrt ( sum(1/data/data) / length )*/

    for (i = 0; i < length; ++i) res += 1.0 / (data[i]*data[i]);
    res = 1.0 / sqrt(res / (REAL8)length);

    break;

   /* these cases require sorted data;
    * we use gsl_sort_index() to avoid in-place modification of the input data
    * by gsl_sort()
    */
   case MATH_OP_ARITHMETIC_MEDIAN:
   case MATH_OP_MINIMUM:
   case MATH_OP_MAXIMUM:
    ; /* empty statement because declaration cannot be first line in a switch case */
    size_t * sortidx;
    if ( ( sortidx = XLALMalloc ( length * sizeof(size_t) )) == NULL ) {
        XLALPrintError ("XLALMalloc(%ld) failed.\n", length );
        XLAL_ERROR ( XLAL_ENOMEM);
    }
    gsl_sort_index(sortidx, data, 1, length);

    switch (optype) {

    case MATH_OP_ARITHMETIC_MEDIAN: /* middle element of sorted data */

      if (length/2 == (length+1)/2) /* length is even */ {
        res = (data[sortidx[length/2-1]] + data[sortidx[length/2]])/2;
      }
      else /* length is odd */ {
        res = data[sortidx[length/2]];
      }

      break;

    case MATH_OP_MINIMUM: /* first element of sorted data */

      res = data[sortidx[0]];
      break;

    case MATH_OP_MAXIMUM: /* last element of sorted data */

      res = data[sortidx[length-1]];
      break;

    default:

      XLAL_ERROR_REAL8(XLAL_EINVAL, "For optype=%i this block should not have been reached.", optype);

    }

    XLALFree ( sortidx );
    break;

  default:

    XLAL_ERROR_REAL8(XLAL_EINVAL, "'%i' is not an implemented math. operation", optype);

  } /* switch (optype) */

  return res;

} /* XLALMathOpOverArray() */


/**
 * Compute an additional normalization factor from the total number of SFTs to be applied after a loop of mathops over SFTs and IFOs.
 *
 * Currently only implemented for:
 * MATH_OP_HARMONIC_SUM (to emulate MATH_OP_HARMONIC_MEAN over the combined array)
 * MATH_OP_POWERMINUS2_SUM (to emulate MATH_OP_POWERMINUS2_MEAN over the combined array)
 *
 * This function exists only as a simple workaround for when we want to compute
 * some types of mathops, e.g. the harmonic or power2 mean,
 * over all SFTs from all detectors together:
 * then call XLALComputePSDandNormSFTPower with PSDmthopSFTs=PSDmthopIFOs=MATH_OP_[HARMONIC/POWERMINUS2]_SUM
 * and then this factor is applied at the end,
 * which gives equivalent results to if we had rewritten the loop order in XLALComputePSDandNormSFTPower
 * to allow for calling MATH_OP_[HARMONIC/POWERMINUS2]_MEAN over the combined array.
 */
REAL8 XLALGetMathOpNormalizationFactorFromTotalNumberOfSFTs (
                            const UINT4 totalNumSFTs, /**< total number of SFTs from all IFOs */
                            const MathOpType optypeSFTs /**< type of operations that was done over SFTs */
                          ) {

  REAL8 factor;

  switch (optypeSFTs) {

  case MATH_OP_ARITHMETIC_SUM: /* sum(data) */
  case MATH_OP_ARITHMETIC_MEAN: /* sum(data)/length  */
  case MATH_OP_HARMONIC_MEAN:
  case MATH_OP_ARITHMETIC_MEDIAN:
  case MATH_OP_MINIMUM:
  case MATH_OP_MAXIMUM:
  case MATH_OP_POWERMINUS2_MEAN:
    XLAL_ERROR_REAL8(XLAL_EINVAL, "For math operation '%i' it makes no sense to call this function!", optypeSFTs);
    break;

  case MATH_OP_HARMONIC_SUM: /* length / sum(1 / data) */
    factor = totalNumSFTs;
    break;

  case MATH_OP_POWERMINUS2_SUM: /*   1 / sqrt ( sum(1 / data/data) )*/
    factor = sqrt(totalNumSFTs);
    break;

  default:

    XLAL_ERROR_REAL8(XLAL_EINVAL, "'%i' is not an implemented math. operation", optypeSFTs);

  } /* switch (optypeSFTs) */

  return factor;

} /* XLALGetMathOpNormalizationFactorFromTotalNumberOfSFTs() */


/**
 * Compute the PSD (power spectral density) and the "normalized power" \f$ P_\mathrm{SFT} \f$ over a MultiSFTVector.
 *
 * \return: a REAL8Vector of the final normalized and averaged/summed PSD;
 * plus the full multiPSDVector array,
 * and optionally a REAL8Vector of normalized SFT power (only if not passed as NULL)
 *
 * If normalizeSFTsInPlace is TRUE, then the inputSFTs will be modified by XLALNormalizeMultiSFTVect().
 * If it is FALSE, an internal copy will be used and the inputSFTs will be returned unmodified.
 */
int
XLALComputePSDandNormSFTPower ( REAL8Vector **finalPSD, /* [out] final PSD averaged over all IFOs and timestamps */
                                MultiPSDVector **multiPSDVector, /* [out] complete MultiPSDVector over IFOs, timestamps and freq-bins */
                                REAL8Vector **normSFT, /* [out] normalised SFT power (optional) */
                                MultiSFTVector *inputSFTs, /* [in] the multi-IFO SFT data */
                                const BOOLEAN returnMultiPSDVector, /* [in] whether to return multiPSDVector */
                                const BOOLEAN returnNormSFT, /* [in] whether to return normSFT vector */
                                const UINT4 blocksRngMed, /* [in] running Median window size */
                                const MathOpType PSDmthopSFTs, /* [in] math operation over SFTs for each IFO (PSD) */
                                const MathOpType PSDmthopIFOs, /* [in] math operation over IFOs (PSD) */
                                const MathOpType nSFTmthopSFTs, /* [in] math operation over SFTs for each IFO (normSFT) */
                                const MathOpType nSFTmthopIFOs, /* [in] math operation over IFOs (normSFT) */
                                const BOOLEAN normalizeByTotalNumSFTs, /* [in] whether to include a final normalization factor derived from the total number of SFTs (over all IFOs); only useful for some mthops */
                                const REAL8 FreqMin, /* [in] starting frequency -> first output bin (if -1: use full SFT data including rngmed bins, else must be >=0) */
                                const REAL8 FreqBand, /* [in] frequency band to cover with output bins (must be >=0) */
                                const BOOLEAN normalizeSFTsInPlace /* [in] if FALSE, a copy of inputSFTs will be used internally and the original will not be modified */
                              )
{

  /* sanity checks */
  XLAL_CHECK ( finalPSD, XLAL_EINVAL );
  XLAL_CHECK ( multiPSDVector, XLAL_EINVAL );
  XLAL_CHECK ( normSFT, XLAL_EINVAL );
  XLAL_CHECK ( inputSFTs && inputSFTs->data && inputSFTs->length>0, XLAL_EINVAL, "inputSFTs must be pre-allocated." );
  XLAL_CHECK ( ( FreqMin>=0 && FreqBand>=0 ) || ( FreqMin==-1 && FreqBand==0 ), XLAL_EINVAL, "Need either both Freqmin>=0 && FreqBand>=0 (truncate PSD band to this range) or FreqMin=-1 && FreqBand=0 (use full band as loaded from SFTs, including rngmed-sidebands.");

  /* get power running-median rngmed[ |data|^2 ] from SFTs *
   * as the output of XLALNormalizeMultiSFTVect(),
   * multiPSD will be a rng-median smoothed periodogram over the normalized SFTs.
   * The inputSFTs themselves will also be normalized in this call
   * unless overruled by normalizeSFTsInPlace option.
  */
  UINT4Vector *numSFTsVec = NULL;
  MultiSFTVector *SFTs = NULL;
  UINT4 numIFOs = inputSFTs->length;
  if ( normalizeSFTsInPlace ) {
      SFTs = inputSFTs;
  }
  else {
      numSFTsVec = XLALCreateUINT4Vector ( numIFOs );
      for (UINT4 X = 0; X < numIFOs; ++X) {
        numSFTsVec->data[X] = inputSFTs->data[X]->length;
      }
      SFTs = XLALCreateEmptyMultiSFTVector ( numSFTsVec );
      for (UINT4 X = 0; X < numIFOs; ++X) {
        for (UINT4 k = 0; k < numSFTsVec->data[X]; ++k) {
          XLAL_CHECK ( XLALCopySFT(&(SFTs->data[X]->data[k]), &(inputSFTs->data[X]->data[k])) == XLAL_SUCCESS, XLAL_EFUNC );
        }
      }
  }
  XLAL_CHECK ( (*multiPSDVector = XLALNormalizeMultiSFTVect ( SFTs, blocksRngMed, NULL )) != NULL, XLAL_EFUNC);
  UINT4 numBinsSFTs = SFTs->data[0]->data[0].data->length;
  REAL8 dFreq = (*multiPSDVector)->data[0]->data[0].deltaF;

  /* restrict to just the "physical" band if requested */
  /* first, figure out the effective PSD bin-boundaries for user */
  UINT4 firstBin, lastBin;
  if ( FreqMin>0 ) { /* user requested a constrained band */
    REAL8 fminSFTs = SFTs->data[0]->data[0].f0;
    REAL8 fmaxSFTs = fminSFTs + numBinsSFTs*dFreq;
    XLAL_CHECK ( ( FreqMin >= fminSFTs ) && ( FreqMin+FreqBand <= fmaxSFTs ), XLAL_EDOM, "Requested frequency range [%f,%f]Hz not covered by available SFT band [%f,%f]Hz", FreqMin, FreqMin+FreqBand, fminSFTs, fmaxSFTs );
    /* check band wide enough for wings */
    UINT4 rngmedSideBandBins = blocksRngMed / 2 + 1; /* truncates down plus add one bin extra safety! */
    REAL8 rngmedSideBand = rngmedSideBandBins * dFreq;
    /* set the actual output range in bins */
    UINT4 minBinPSD = XLALRoundFrequencyDownToSFTBin(FreqMin, dFreq);
    UINT4 maxBinPSD = XLALRoundFrequencyUpToSFTBin(FreqMin + FreqBand, dFreq) - 1;
    UINT4 minBinSFTs = XLALRoundFrequencyDownToSFTBin(fminSFTs, dFreq);
    UINT4 maxBinSFTs = XLALRoundFrequencyUpToSFTBin(fmaxSFTs, dFreq) - 1;
    if ( ( minBinPSD < minBinSFTs ) || ( maxBinPSD > maxBinSFTs ) ) {
      XLAL_ERROR( XLAL_ERANGE, "Requested frequency range [%f,%f)Hz means rngmedSideBand=%fHz (%d bins) would spill over available SFT band [%f,%f)Hz",
                  FreqMin, FreqMin+FreqBand, rngmedSideBand, rngmedSideBandBins, fminSFTs, fmaxSFTs );
    }
    UINT4 binsBand = maxBinPSD - minBinPSD + 1;
    firstBin = minBinPSD - minBinSFTs;
    lastBin = firstBin + binsBand - 1;
  }
  else { /* output all bins loaded from SFTs (includes rngmed-sidebands) */
    firstBin = 0;
    lastBin = numBinsSFTs - 1;
  }
  XLAL_PRINT_INFO("loaded SFTs have %d bins, effective PSD output band is [%d, %d]", numBinsSFTs, firstBin, lastBin);
  /* now truncate */
  XLAL_CHECK ( XLALCropMultiPSDandSFTVectors ( *multiPSDVector, SFTs, firstBin, lastBin ) == XLAL_SUCCESS, XLAL_EFUNC, "Failed call to XLALCropMultiPSDandSFTVectors (multiPSDVector, SFTs, firstBin=%d, lastBin=%d)", firstBin, lastBin );

  /* number of raw bins in final PSD */
  UINT4 numBins = (*multiPSDVector)->data[0]->data[0].data->length;

  /* allocate main output struct, using cropped length */
  XLAL_CHECK ( (*finalPSD = XLALCreateREAL8Vector ( numBins )) != NULL, XLAL_ENOMEM, "Failed to create REAL8Vector for finalPSD.");

  /* number of IFOs */
  REAL8Vector *overIFOs = NULL; /* one frequency bin over IFOs */
  XLAL_CHECK ( (overIFOs = XLALCreateREAL8Vector ( numIFOs )) != NULL, XLAL_ENOMEM, "Failed to create REAL8Vector for overIFOs array.");

  /* maximum number of SFTs */
  UINT4 maxNumSFTs = 0;
  for (UINT4 X = 0; X < numIFOs; ++X) {
    maxNumSFTs = GSL_MAX(maxNumSFTs, (*multiPSDVector)->data[X]->length);
  }

  REAL8Vector *overSFTs = NULL; /* one frequency bin over SFTs */
  XLAL_CHECK ( (overSFTs = XLALCreateREAL8Vector ( maxNumSFTs )) != NULL, XLAL_ENOMEM, "Failed to create REAL8Vector for overSFTs array.");

  /* normalize rngmd(power) to get proper *single-sided* PSD: Sn = (2/Tsft) rngmed[|data|^2]] */
  REAL8 normPSD = 2.0 * dFreq;

  /* optionally include a normalizing factor for when we want to compute e.g. the harmonic or power2 mean over all SFTs from all detectors together:
   * call with PSDmthopSFTs=PSDmthopIFOs=MATH_OP_[HARMONIC/POWERMINUS2]_SUM
   * and then this factor is applied at the end,
   * which gives equivalent results to if we had rewritten the loop order
   * to allow for calling MATH_OP_[HARMONIC/POWERMINUS2]_MEAN over the combined array.
   */
  REAL8 totalNumSFTsNormalizingFactor = 1.;
  if ( normalizeByTotalNumSFTs ) {
    UINT4 totalNumSFTs = 0;
    for (UINT4 X = 0; X < numIFOs; ++X) {
      totalNumSFTs += (*multiPSDVector)->data[X]->length;
    }
    totalNumSFTsNormalizingFactor = XLALGetMathOpNormalizationFactorFromTotalNumberOfSFTs ( totalNumSFTs, PSDmthopSFTs );
  }

  XLAL_PRINT_INFO("Computing spectrogram and PSD ...");
  /* loop over frequency bins in final PSD */
  for (UINT4 k = 0; k < numBins; ++k) {

    /* loop over IFOs */
    for (UINT4 X = 0; X < numIFOs; ++X) {

      /* number of SFTs for this IFO */
      UINT4 numSFTs = (*multiPSDVector)->data[X]->length;

      /* copy PSD frequency bins and normalise multiPSDVector for later use */
      for (UINT4 alpha = 0; alpha < numSFTs; ++alpha) {
        (*multiPSDVector)->data[X]->data[alpha].data->data[k] *= normPSD;
        overSFTs->data[alpha] = (*multiPSDVector)->data[X]->data[alpha].data->data[k];
      }

      /* compute math. operation over SFTs for this IFO */
      overIFOs->data[X] = XLALMathOpOverArray(overSFTs->data, numSFTs, PSDmthopSFTs);
      XLAL_CHECK ( !XLAL_IS_REAL8_FAIL_NAN(overIFOs->data[X]), XLAL_EFUNC, "XLALMathOpOverArray() returned NAN for overIFOs->data[X=%d]", X );

    } /* for IFOs X */

    /* compute math. operation over IFOs for this frequency */
    (*finalPSD)->data[k] = XLALMathOpOverArray(overIFOs->data, numIFOs, PSDmthopIFOs);
    XLAL_CHECK ( !XLAL_IS_REAL8_FAIL_NAN((*finalPSD)->data[k]), XLAL_EFUNC, "XLALMathOpOverArray() returned NAN for finalPSD->data[k=%d]", k );

    if ( normalizeByTotalNumSFTs ) {
       (*finalPSD)->data[k] *= totalNumSFTsNormalizingFactor;
    }

  } /* for freq bins k */
  XLAL_PRINT_INFO("done.");

  /* compute normalised SFT power */
  if ( returnNormSFT ) {
    XLAL_PRINT_INFO("Computing normalised SFT power ...");
    XLAL_CHECK ( (*normSFT = XLALCreateREAL8Vector ( numBins )) != NULL, XLAL_ENOMEM, "Failed to create REAL8Vector for normSFT.");

    /* loop over frequency bins in SFTs */
    for (UINT4 k = 0; k < numBins; ++k) {

      /* loop over IFOs */
      for (UINT4 X = 0; X < numIFOs; ++X) {

        /* number of SFTs for this IFO */
        UINT4 numSFTs = SFTs->data[X]->length;

        /* compute SFT power */
        for (UINT4 alpha = 0; alpha < numSFTs; ++alpha) {
          COMPLEX8 bin = SFTs->data[X]->data[alpha].data->data[k];
          overSFTs->data[alpha] = crealf(bin)*crealf(bin) + cimagf(bin)*cimagf(bin);
        }

        /* compute math. operation over SFTs for this IFO */
        overIFOs->data[X] = XLALMathOpOverArray(overSFTs->data, numSFTs, nSFTmthopSFTs);
        XLAL_CHECK ( !XLAL_IS_REAL8_FAIL_NAN(overIFOs->data[X]), XLAL_EFUNC, "XLALMathOpOverArray() returned NAN for overIFOs->data[X=%d]", X );

      } /* over IFOs */

      /* compute math. operation over IFOs for this frequency */
      (*normSFT)->data[k] = XLALMathOpOverArray(overIFOs->data, numIFOs, nSFTmthopIFOs);
      XLAL_CHECK ( !XLAL_IS_REAL8_FAIL_NAN((*normSFT)->data[k]), XLAL_EFUNC, "XLALMathOpOverArray() returned NAN for normSFT->data[k=%d]", k );

    } /* over freq bins */
    XLAL_PRINT_INFO("done.");
  } /* if returnNormSFT */

  XLALDestroyREAL8Vector ( overSFTs );
  XLALDestroyREAL8Vector ( overIFOs );

  if ( !returnMultiPSDVector ) {
    XLALDestroyMultiPSDVector ( *multiPSDVector );
    *multiPSDVector = NULL;
  }

  if ( !normalizeSFTsInPlace ) {
    XLALDestroyUINT4Vector ( numSFTsVec );
    XLALDestroyMultiSFTVector ( SFTs );
  }

  return XLAL_SUCCESS;

} /* XLALComputePSDandNormSFTPower() */


/**
 * Compute the PSD (power spectral density) over a MultiSFTVector.
 * This is just a convenience wrapper around XLALComputePSDandNormSFTPower()
 * without the extra options for normSFTpower computation.
 *
 * \return a REAL8Vector of the final normalized and averaged/summed PSD.
 *
 * NOTE: contrary to earlier versions of this function, it no longer modifies
 * the inputSFTs in-place, so now it can safely be called multiple times.
 *
 */
int
XLALComputePSDfromSFTs ( REAL8Vector **finalPSD, /* [out] final PSD averaged over all IFOs and timestamps */
                         MultiSFTVector *inputSFTs, /* [in] the multi-IFO SFT data */
                         const UINT4 blocksRngMed, /* [in] running Median window size */
                         const MathOpType PSDmthopSFTs, /* [in] math operation over SFTs for each IFO (PSD) */
                         const MathOpType PSDmthopIFOs, /* [in] math operation over IFOs (PSD) */
                         const BOOLEAN normalizeByTotalNumSFTs, /* [in] whether to include a final normalization factor derived from the total number of SFTs (over all IFOs); only useful for some mthops */
                         const REAL8 FreqMin, /* [in] starting frequency -> first output bin (if -1: use full SFT data including rngmed bins, else must be >=0) */
                         const REAL8 FreqBand /* [in] frequency band to cover with output bins (must be >=0) */
                       )
{

  MultiPSDVector *multiPSDVector = NULL;
  REAL8Vector *normSFT = NULL;
  XLAL_CHECK ( XLALComputePSDandNormSFTPower ( finalPSD,
                      &multiPSDVector,
                      &normSFT,
                      inputSFTs,
                      0, // returnMultiPSDVector
                      0, // returnNormSFT
                      blocksRngMed,
                      PSDmthopSFTs,
                      PSDmthopIFOs,
                      0, // nSFTmthopSFTs
                      0, // nSFTmthopIFOs
                      normalizeByTotalNumSFTs,
                      FreqMin,
                      FreqBand,
                      FALSE // normalizeSFTsInPlace
  ) == XLAL_SUCCESS, XLAL_EFUNC );

  return XLAL_SUCCESS;

} /* XLALComputePSDfromSFTs() */


/**
 * Dump complete multi-PSDVector over IFOs, timestamps and frequency-bins into
 * per-IFO ASCII output-files 'outbname-IFO'
 *
 */
int
XLALDumpMultiPSDVector ( const CHAR *outbname,			/**< output basename 'outbname' */
                         const MultiPSDVector *multiPSDVect	/**< multi-psd vector to output */
                  )
{
  /* check input consistency */
  if ( outbname == NULL ) {
    XLALPrintError ("%s: NULL input 'outbname'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( multiPSDVect == NULL ) {
    XLALPrintError ("%s: NULL input 'multiPSDVect'\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }
  if ( multiPSDVect->length == 0 || multiPSDVect->data==0 ) {
    XLALPrintError ("%s: invalid multiPSDVect input (length=0 or data=NULL)\n", __func__ );
    XLAL_ERROR ( XLAL_EINVAL );
  }

  CHAR *fname;
  FILE *fp;

  UINT4 len = strlen ( outbname ) + 4;
  if ( ( fname = XLALMalloc ( len * sizeof(CHAR) )) == NULL ) {
    XLALPrintError ("%s: XLALMalloc(%d) failed.\n", __func__, len );
    XLAL_ERROR ( XLAL_ENOMEM);
  }

  UINT4 numIFOs = multiPSDVect->length;
  UINT4 X;
  for ( X = 0; X < numIFOs; X ++ )
    {
      PSDVector *thisPSDVect = multiPSDVect->data[X];
      char buf[100];

      sprintf ( fname, "%s-%c%c", outbname, thisPSDVect->data[0].name[0], thisPSDVect->data[0].name[1] );

      if ( ( fp = fopen( fname, "wb" ))  == NULL ) {
        XLALPrintError ("%s: Failed to open PSDperSFT file '%s' for writing!\n", __func__, fname );
        XLALFree ( fname );
        XLAL_ERROR ( XLAL_ESYS );
      }

      REAL8 f0       = thisPSDVect->data[0].f0;
      REAL8 dFreq    = thisPSDVect->data[0].deltaF;
      UINT4 numFreqs = thisPSDVect->data[0].data->length;
      UINT4 iFreq;

      /* write comment header line into this output file */
      /* FIXME: output code-version/cmdline/history info */
      fprintf(fp,"%%%% first line holds frequencies [Hz] of PSD-Columns\n");
      /* loop over frequency and output comnment-header markers */
      fprintf(fp,"%%%%%-17s", "dummy");
      for ( iFreq = 0; iFreq < numFreqs; iFreq ++ )
        {
          sprintf (buf, "f%d [Hz]", iFreq + 1 );
          fprintf (fp, " %-21s", buf );
        }
      fprintf (fp, "\n");

      /* write parseable header-line giving bin frequencies for PSDs */
      fprintf (fp, "%-19d", -1 );
      for (iFreq = 0; iFreq < numFreqs; iFreq++ )
        fprintf (fp, " %-21.16g", f0 + iFreq * dFreq );
      fprintf (fp, "\n\n\n");

      /* output another header line describing the following format "ti[GPS] PSD(f1) ... " */
      fprintf(fp,"%%%%%-17s", "ti[GPS]");
      for ( iFreq = 0; iFreq < numFreqs; iFreq ++ )
        {
          sprintf (buf, "PSD(f%d)", iFreq + 1 );
          fprintf (fp, " %-21s", buf );
        }
      fprintf (fp, "\n");

      /* loop over timestamps: dump all PSDs over frequencies into one line */
      UINT4 numTS = thisPSDVect->length;
      UINT4 iTS;
      for ( iTS = 0; iTS < numTS; iTS++ )
        {
          REAL8FrequencySeries *thisPSD = &thisPSDVect->data[iTS];

          /* first output timestamp GPS time for this line */
          REAL8 tGPS = XLALGPSGetREAL8( &thisPSD->epoch );
          fprintf (fp, "%-19.18g", tGPS );

          /* some internal consistency/paranoia checks */
          if ( ( f0 != thisPSD->f0) || ( dFreq != thisPSD->deltaF ) || (numFreqs != thisPSD->data->length ) ) {
            XLALPrintError ("%s: %d-th timestamp %f: inconsistent PSDVector: f0 = %g : %g,  dFreq = %g : %g, numFreqs = %d : %d \n",
                            __func__, iTS, tGPS, f0, thisPSD->f0, dFreq, thisPSD->deltaF, numFreqs, thisPSD->data->length );
            XLALFree ( fname );
            fclose ( fp );
            XLAL_ERROR ( XLAL_EDOM );
          }

          /* loop over all frequencies and dump PSD-value */
          for ( iFreq = 0; iFreq < numFreqs; iFreq ++ )
            fprintf (fp, " %-21.16g", thisPSD->data->data[iFreq] );

          fprintf (fp, "\n");

        } /* for iTS < numTS */

      fclose ( fp );

    } /* for X < numIFOs */

  XLALFree ( fname );

  return XLAL_SUCCESS;

} /* XLALDumpMultiPSDVector() */


/**
 * Write a PSD vector as an ASCII table to an open output file pointer.
 *
 * Adds a column header string, but version or commandline info
 * needs to be printed by the caller.
 *
 * Standard and (optional) columns: FreqBinStart (FreqBinEnd) PSD (normSFTpower)
 */
int
XLALWritePSDtoFilePointer ( FILE *fpOut,              /**< output file pointer */
                            REAL8Vector *PSDVect,     /**< required: PSD vector */
                            REAL8Vector *normSFTVect, /**< optional: normSFT vector (can be NULL if outputNormSFT==False) */
                            BOOLEAN outputNormSFT,    /**< output a column for normSFTVect? */
                            BOOLEAN outFreqBinEnd,    /**< output a column for the end frequency of each bin? */
                            INT4 PSDmthopBins,        /**< math operation for binning of the PSD */
                            INT4 nSFTmthopBins,       /**< math operation for binning of the normSFT */
                            INT4 binSize,             /**< output bin size (in number of input bins, 1 to keep input binning) */
                            INT4 binStep,             /**< output bin step (in number of input bins, 1 to keep input binning) */
                            REAL8 Freq0,              /**< starting frequency of inputs */
                            REAL8 dFreq               /**< frequency step of inputs */
                          ) {

    /* input sanity checks */
    XLAL_CHECK ( ( PSDVect != NULL ) && ( PSDVect->data != NULL ) && ( PSDVect->length > 0 ), XLAL_EINVAL, "PSDVect must be allocated and not zero length." );
    if (outputNormSFT) {
      XLAL_CHECK ( ( normSFTVect != NULL ) && ( normSFTVect->data != NULL ), XLAL_EINVAL, "If requesting outputNormSFT, normSFTVect must be allocated" );
      XLAL_CHECK ( normSFTVect->length == PSDVect->length, XLAL_EINVAL, "Lengths of PSDVect and normSFTVect do not match (%d!=%d)", PSDVect->length, normSFTVect->length );
    }
    XLAL_CHECK ( binSize > 0, XLAL_EDOM, "output bin size must be >0");
    XLAL_CHECK ( binStep > 0, XLAL_EDOM, "output bin step must be >0");
    XLAL_CHECK ( Freq0 >= 0., XLAL_EDOM, "starting frequency must be nonnegative");
    XLAL_CHECK ( dFreq >= 0., XLAL_EDOM, "frequency step must be nonnegative");

    /* work out total number of bins */
    UINT4 numBins = (UINT4)floor((PSDVect->length - binSize) / binStep) + 1;

    /* write column headings */
    fprintf(fpOut,"%%%% columns:\n%%%% FreqBinStart");
    if (outFreqBinEnd)
      fprintf(fpOut," FreqBinEnd");
    fprintf(fpOut," PSD");
    if (outputNormSFT)
      fprintf(fpOut," normSFTpower");
    fprintf(fpOut,"\n");

    for (UINT4 k = 0; k < numBins; ++k) {
      UINT4 b = k * binStep;

      REAL8 f0 = Freq0 + b * dFreq;
      REAL8 f1 = f0 + binStep * dFreq;
      fprintf(fpOut, "%f", f0);
      if (outFreqBinEnd)
        fprintf(fpOut, "   %f", f1);

      REAL8 psd = XLALMathOpOverArray(&(PSDVect->data[b]), binSize, PSDmthopBins);
      XLAL_CHECK ( !XLAL_IS_REAL8_FAIL_NAN(psd), XLAL_EFUNC, "XLALMathOpOverArray() returned NAN for psd[k=%d]", k );

      fprintf(fpOut, "   %e", psd);

      if (outputNormSFT) {
        REAL8 nsft = XLALMathOpOverArray(&(normSFTVect->data[b]), binSize, nSFTmthopBins);
        XLAL_CHECK ( !XLAL_IS_REAL8_FAIL_NAN(nsft), XLAL_EFUNC, "XLALMathOpOverArray() returned NAN for nsft[k=%d]", k );

        fprintf(fpOut, "   %f", nsft);
      }

      fprintf(fpOut, "\n");
    } // k < numBins

    return XLAL_SUCCESS;

} /* XLALWritePSDtoFile() */
