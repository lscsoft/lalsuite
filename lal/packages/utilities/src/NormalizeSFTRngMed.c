/*
*  Copyright (C) 2007 Badri Krishnan, Reinhard Prix
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

/**
 * \file
 * \ingroup NormalizeSFTRndMed.c
 * \author Badri Krishnan and Alicia Sintes
 * \date $Date$
 * \brief Normalizes SFTs based on their noise floor calculated using the running median
 *
 * $Id$
 *
 * History: Created by B. Krishnan Aug, 2004
 *       Taken from SFTbin.c and PeakSelect.c from hough dir in lalapps
 *
 *

 \par Description

This module contains functions for normalizing SFTs.  Currently two normalizations
are supported.  Given SFT data \f$\tilde{x}_k \f$ where \f$ k\f$ labels a frequency bin,
the normalized SFT is either \f$ \tilde{x}_k/\sqrt{ < |\tilde{x}_k|^2 >} \f$ or
\f$ \sqrt{2N} \tilde{x}_k/\sqrt{ < |\tilde{x}_k|^2 >} \f$, where \f$ N \f$ is the number
of frequency bins in the SFT.   The first normalization
ensures that the SFT power follows an exponential distribution with unit mean
(if the SFT data is distributed normally), while the second normalization is appropriate
in the time domain.  In either case, the mean of \f$ |\tilde{x}_k|^2 \f$ is
estimated using the median, suitably normalized assuming that the power is
distributed is exponentially.


\par Uses
\code
LALSFTtoPeriodogram ()
LALPeriodoToRngmed ()
LALNormalizeSFT ()
LALNormalizeSFTVect ()
LALNormalizeMultiSFTVect ()
\endcode


The function LALNormalizeSFTVect() takes as input a vector of SFTs and normalizes
them.  This function calls the functions LALNormalizeSFT() which normalizes a
single SFT, LALSFTtoPeriodogram() which calculates the \f$ |\tilde{x}|^2 \f$ and
LALPeriodoToRngmed () which applies the running median algorithm to find a vector
of medians.  The function LALNormalizeMultiSFTVect() normalizes a multi-IFO collection
of SFT vectors and also returns a collection of power-estimates for these vectors using
the Running median method.

*/

#include <lal/NormalizeSFTRngMed.h>


NRCSID (NORMALIZESFTRNGMEDC, "$Id$");


/** Calculate the "periodogram" of an SFT, ie the modulus-squares of the SFT-data.
 */
void
LALSFTtoPeriodogram (LALStatus    *status,
		     REAL8FrequencySeries    *periodo,	/**< [out] mod squares of SFT data */
		     const COMPLEX8FrequencySeries *SFT	/**< [in] input SFT */
		     )
{
  UINT4     length, j;
  REAL8    *out;
  COMPLEX8 *in;

  INITSTATUS (status, "LALSFTtoPeriodogram", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL */
  ASSERT (periodo, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (periodo->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (periodo->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (periodo->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (SFT, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (SFT->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (SFT->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (SFT->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  /* copy values from SFT */
  strcpy ( periodo->name, SFT->name );
  periodo->epoch.gpsSeconds = SFT->epoch.gpsSeconds;
  periodo->epoch.gpsNanoSeconds = SFT->epoch.gpsNanoSeconds;
  periodo->f0 = SFT->f0;
  periodo->deltaF = SFT->deltaF;

  /* check lengths are same */
  length = SFT->data->length;
  ASSERT (length == periodo->data->length, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);

  out = periodo->data->data;
  in = SFT->data->data;

  for (j=0; j<length; j++) {
    /* extra-paranoia: make absolutely sure that the calculation below is in REAL8
     * in order to avoid underflow-problems (data 'in' can be of order ~ 1e-20 )
     */
    *out = ((REAL8)in->re)*((REAL8)in->re) + ((REAL8)in->im)*((REAL8)in->im);
    ++out;
    ++in;
  }

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* end LALSFTtoPeriodogram() */


/** Calculates running median over a single periodogram.
*/
void
LALPeriodoToRngmed (LALStatus  *status,
		    REAL8FrequencySeries  *rngmed,		/**< [out] resulting 'smoothed' periodogram */
		    const REAL8FrequencySeries  *periodo,	/**< [in] input periodogram */
		    UINT4 blockSize			/**< Running median block size */
		    )
{
  UINT4 blocks2;
  UINT4 j;
  UINT4 length;
  LALRunningMedianPar rngMedPar;
  REAL8Sequence mediansV, inputV;
  REAL8 medianBias;

  INITSTATUS (status, "LALPeriodoToRngmed", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL */
  ASSERT (periodo, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (periodo->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (periodo->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (periodo->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (rngmed, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (rngmed->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (rngmed->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (rngmed->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (blockSize > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);


  /* copy values from the periodogram */
  strcpy ( rngmed->name, periodo->name );
  rngmed->epoch.gpsSeconds = periodo->epoch.gpsSeconds;
  rngmed->epoch.gpsNanoSeconds = periodo->epoch.gpsNanoSeconds;
  rngmed->f0 = periodo->f0;
  rngmed->deltaF = periodo->deltaF;

  /* check lengths are same */
  length = periodo->data->length;
  ASSERT (length == rngmed->data->length, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);

  if ( length <= blockSize ) {
    XLALPrintError ("Need at least %d bins in SFT (have %d) to perform running median!\n", blockSize + 1, length );
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  }

  blocks2 = blockSize/2 - 1; /* integer division */

  rngMedPar.blocksize = (UINT4)blockSize;
  inputV.length = length;
  inputV.data = periodo->data->data;
  mediansV.length= length - blockSize + 1;
  mediansV.data = rngmed->data->data + blocks2;

  TRY( LALDRunningMedian2(status->statusPtr, &mediansV, &inputV, rngMedPar), status);

  /* copy values in the wings */
  for (j=0; j<blocks2; j++)
    rngmed->data->data[j] = rngmed->data->data[blocks2];

  for (j=blocks2+length-blockSize+1; j<length; j++)
    rngmed->data->data[j] = rngmed->data->data[blocks2 + length-blockSize];

  /* get the bias factor -- for estimating the mean from the median */
  TRY ( LALRngMedBias( status->statusPtr, &medianBias, blockSize ), status);

  /* normalize by the bias factor */
  for (j=0; j<length; j++)
    rngmed->data->data[j] /= medianBias;

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* LALPeriodoToRngmed() */


/** Calculates a smoothed (running-median) periodogram for the given SFT.
 */
void
LALSFTtoRngmed (LALStatus  *status,
		REAL8FrequencySeries  *rngmed,		/**< [out] running-median smoothed periodo [must be allocated!] */
		const COMPLEX8FrequencySeries *sft,	/**< [in]  input SFT */
		UINT4 blockSize				/**< Running median block size */
		)
{
  REAL8FrequencySeries periodo;
  INT4 length;

  INITSTATUS (status, "LALSFTtoRngmed", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL */
  ASSERT (sft, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sft->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sft->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (sft->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (rngmed, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (rngmed->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (rngmed->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (rngmed->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  length = sft->data->length;

  periodo.data = NULL;
  if ( (periodo.data = (REAL8Sequence *)LALCalloc(1, sizeof(REAL8Sequence))) == NULL) {
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }
  periodo.data->length = length;
  if ( (periodo.data->data = (REAL8 *)LALCalloc( length,  sizeof(REAL8))) == NULL) {
    LALFree(periodo.data);
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  /* calculate the periodogram */
  LALSFTtoPeriodogram (status->statusPtr, &periodo, sft);
  BEGINFAIL (status) {
    LALFree (periodo.data->data);
    LALFree (periodo.data);
  } ENDFAIL (status);

  /* calculate the rngmed */
  if ( blockSize > 0 )
    {
      LALPeriodoToRngmed (status->statusPtr, rngmed, &periodo, blockSize);
      BEGINFAIL (status) {
	LALFree (periodo.data->data);
	LALFree (periodo.data);
      } ENDFAIL (status);

      /* free memory */
      LALFree(periodo.data->data);
      LALFree(periodo.data);
    }
  else
    (*rngmed) = periodo;	/* struct-copy */

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* LALSFTtoRngmed() */



/** Normalize an sft based on RngMed estimated PSD.
 */
void
LALNormalizeSFT (LALStatus           *status,
		 REAL8FrequencySeries *rngmed, 	/**< [out] rng-median smoothed periodogram over SFT (Tsft*Sn/2) */
		 SFTtype              *sft,     /**< SFT to be normalized */
		 UINT4                blockSize)/**< Running median block size for rngmed calculation */
{
  UINT4 j;
  REAL8 Tsft_Sn_b2;	/* Wiener-Kinchine: E[|data|^2] = Tsft * Sn / 2 */

  INITSTATUS (status, "LALNormalizeSFT", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL and other sanity checks*/
  ASSERT (sft, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sft->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sft->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (sft->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  ASSERT (rngmed, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (rngmed->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (rngmed->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (rngmed->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  /* make sure there is no size mismatch */
  if ( rngmed->data->length != sft->data->length ) {
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  }

  /* calculate the rngmed */
  TRY (LALSFTtoRngmed (status->statusPtr, rngmed, sft, blockSize), status);

  /* loop over sft and normalize */
  for (j = 0; j < sft->data->length; j++) {

    Tsft_Sn_b2 = rngmed->data->data[j];

    /* frequency domain normalization */
    sft->data->data[j].re /= sqrt(Tsft_Sn_b2);
    sft->data->data[j].im /= sqrt(Tsft_Sn_b2);
  }

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* end LALNormalizeSFT() */


/** Function for normalizing a vector of SFTs.
 */
void
LALNormalizeSFTVect (LALStatus  *status,
		     SFTVector  *sftVect,	/**< [in/out] pointer to a vector of SFTs which will be normalized */
		     UINT4     blockSize	/**< Running median window size */
		     )
{
  UINT4 j, lengthsft;
  REAL8FrequencySeries *rngmed = NULL;

  INITSTATUS (status, "LALNormalizeSFT", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL and other sanity checks*/
  ASSERT (sftVect, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sftVect->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sftVect->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);

  /* memory allocation of rngmed using length of first sft
     -- assume all sfts have the same length*/
  lengthsft = sftVect->data->data->length;

  /* allocate memory for a single rngmed */
  if ( (rngmed = (REAL8FrequencySeries *)LALCalloc(1, sizeof(REAL8FrequencySeries))) == NULL){
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  rngmed->data = NULL;
  if ( (rngmed->data = (REAL8Sequence *)LALCalloc(1, sizeof(REAL8Sequence))) == NULL) {
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  rngmed->data->length = lengthsft;
  if ( (rngmed->data->data = (REAL8 *)LALCalloc( lengthsft, sizeof(REAL8))) == NULL) {
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  /* loop over sfts and normalize them */
  for (j = 0; j < sftVect->length; j++) {
    SFTtype *sft;

    if ( (sft = sftVect->data + j) == NULL) {
      ABORT( status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
    }

    if (sft->data == NULL) {
      ABORT( status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
    }

    /* check there is no mismatch in length of sft */
    if (sft->data->length != lengthsft) {
      ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
    }

    if (sft->data->data == NULL) {
      ABORT (status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
    }

    /* call sft normalization function */
    LALNormalizeSFT (status->statusPtr, rngmed, sft, blockSize);
    BEGINFAIL (status) {
      LALFree (rngmed->data->data);
      LALFree (rngmed->data);
      LALFree (rngmed);
    } ENDFAIL (status);

  } /* for loop over sfts */

  /* free memory for psd */
  LALFree(rngmed->data->data);
  LALFree(rngmed->data);
  LALFree(rngmed);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* end LALNormalizeSFTVect() */



/** Function for normalizing a multi vector of SFTs in a multi IFO search and also
 * returns the running-median estimates of the power.
 */
void
LALNormalizeMultiSFTVect (LALStatus      *status,
			  MultiPSDVector **multiRngmed,	/**< [out] multi running-median power estimates of input SFTs */
			  MultiSFTVector *multsft,	/**< [in/out] multi-vector of SFTs which will be normalized */
			  UINT4          blockSize	/**< Running median window size */
			  )
{

  UINT4 k, j; /* k loops over IFOs and j over SFTs for each IFO */
  UINT4 jCleanUp, kCleanUp; /* indices used in clean-up loops */
  UINT4 numifo, numsft;
  MultiPSDVector *ret = NULL;

  INITSTATUS (status, "LALNormalizeMultiSFT", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL and other sanity checks*/
  ASSERT (multsft, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (multsft->length, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (multsft->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  ASSERT (multiRngmed, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT ( *multiRngmed == NULL, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  /* first memory allocation for multipsd structure */
  if ( (ret = (MultiPSDVector *)LALCalloc(1, sizeof(MultiPSDVector))) == NULL) {
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  ret->length = numifo = multsft->length;
  if ( (ret->data = (PSDVector **)LALCalloc( numifo, sizeof(PSDVector *))) == NULL) {
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  /* loop over ifos */
  for ( k = 0; k < numifo; k++) {

    /* second memory allocation for psd vector */
    if ( (ret->data[k] = (PSDVector *)LALCalloc(1, sizeof(PSDVector))) == NULL) {
      ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
    }

    ret->data[k]->length = numsft = multsft->data[k]->length;
    if ( (ret->data[k]->data = (REAL8FrequencySeries *)LALCalloc(numsft, sizeof(REAL8FrequencySeries))) == NULL) {
      ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
    }

    /* loop over sfts for each ofo */
    for (j = 0; j < numsft; j++) {

      SFTtype *sft;
      UINT4 lengthsft;

      if ( (sft = multsft->data[k]->data + j) == NULL) {
	ABORT( status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
      }

      if (sft->data == NULL) {
	ABORT( status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
      }

      if (sft->data->length == 0) {
	ABORT( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
      }

      if (sft->data->data == NULL) {
	ABORT( status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
      }

      /* final memory allocation for psd */
      ret->data[k]->data[j].data = NULL;
      if ( (ret->data[k]->data[j].data = (REAL8Sequence *)LALCalloc(1, sizeof(REAL8Sequence))) == NULL) {
	ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
      }

      ret->data[k]->data[j].data->length = lengthsft = sft->data->length;
      if ( (ret->data[k]->data[j].data->data = (REAL8 *)LALCalloc( lengthsft, sizeof(REAL8))) == NULL) {
	ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
      }

      LALNormalizeSFT (status->statusPtr, ret->data[k]->data + j, sft, blockSize);
      BEGINFAIL (status) {
        /* clean up for this value of k -- note that j and k have not been incremented at this stage*/
	for ( jCleanUp = 0; jCleanUp < j+1; jCleanUp++) {
	  LALFree( ret->data[k]->data[jCleanUp].data->data);
	  LALFree( ret->data[k]->data[jCleanUp].data);
	}
	LALFree( ret->data[k]->data);
	LALFree( ret->data[k]);
	/* clean up for previous values of k */
	for ( kCleanUp = 0; kCleanUp < k; kCleanUp++) {
	  for ( jCleanUp = 0; jCleanUp < multsft->data[kCleanUp]->length; jCleanUp++) {
	    LALFree( ret->data[kCleanUp]->data[jCleanUp].data->data);
	    LALFree( ret->data[kCleanUp]->data[jCleanUp].data);
	  }
	  LALFree( ret->data[kCleanUp]->data);
	  LALFree( ret->data[kCleanUp]);
	}
	/* clean up memory allocated outside loop */
	LALFree(ret->data);
	LALFree(ret);
      } ENDFAIL (status);
    } /* loop over sfts ++j */
  } /* loop over ifos ++k */

  (*multiRngmed) = ret;

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
} /* LALNormalizeMultiSFTVect() */


/** Calculate the cross-correlation periodogram from 2 SFTs.
 */
void
LALSFTstoCrossPeriodogram (LALStatus    *status,
			   REAL8FrequencySeries *periodo,	/**< [out] modulus square of SFT data  */
			   const COMPLEX8FrequencySeries *sft1,	/**< [in] pointer to first SFT  */
			   const COMPLEX8FrequencySeries *sft2	/**< [in] pointer to second SFT */
			   )
{
  UINT4     length, j;
  REAL8    *out;
  COMPLEX8 *in1, *in2;

  INITSTATUS (status, "LALSFTstoCrossPeriodogram", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL */
  ASSERT (periodo, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (periodo->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (periodo->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (periodo->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sft1, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sft1->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sft1->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (sft1->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  ASSERT (sft2, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sft2->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);
  ASSERT (sft2->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (sft2->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  /* make sure both sfts are consistent in frequency and freq. band
     -- time stamps need not be consistent */
  ASSERT (sft2->data->length == sft1->data->length, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (sft2->f0 == sft1->f0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (sft2->deltaF == sft1->deltaF, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);

  /* copy values from SFT */
  /*   periodo->epoch.gpsSeconds = sft1->epoch.gpsSeconds; */
  /*   periodo->epoch.gpsNanoSeconds = sft1->epoch.gpsNanoSeconds; */
  periodo->f0 = sft1->f0;
  periodo->deltaF = sft1->deltaF;

  /* check lengths are same */
  length = sft1->data->length;
  ASSERT (length == periodo->data->length, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);

  out = periodo->data->data;
  in1 = sft1->data->data;
  in2 = sft2->data->data;

  for (j=0; j<length; j++) {
    /* extra-paranoia: make absolutely sure that the calculation below is in REAL8
     * in order to avoid underflow-problems (data 'in' can be of order ~ 1e-20 )
     */
    *out = ((REAL8)in1->re)*((REAL8)in2->re) + ((REAL8)in1->im)*((REAL8)in2->im);
    ++out;
    ++in1;
    ++in2;
  }

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* end LALSFTtoPeriodogram() */
