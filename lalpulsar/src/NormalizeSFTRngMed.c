/*
 * (C) 2012 Reinhard Prix
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

#include <lal/NormalizeSFTRngMed.h>

static const LALStatus empty_LALStatus;

/**
 * \addtogroup NormalizeSFTRngMed_h
 * \author Badri Krishnan and Alicia Sintes
 * \brief Normalizes SFTs based on their noise floor calculated using the running median
 *
 * History: Created by B. Krishnan Aug, 2004
 * Taken from SFTbin.c and PeakSelect.c from hough dir in lalapps
 *
 * ### Description ###
 *
 * This module contains functions for normalizing SFTs.  Currently two normalizations
 * are supported.  Given SFT data \f$\tilde{x}_k \f$ where \f$ k\f$ labels a frequency bin,
 * the normalized SFT is either \f$ \tilde{x}_k/\sqrt{ < |\tilde{x}_k|^2 >} \f$ or
 * \f$ \sqrt{2N} \tilde{x}_k/\sqrt{ < |\tilde{x}_k|^2 >} \f$, where \f$ N \f$ is the number
 * of frequency bins in the SFT.   The first normalization
 * ensures that the SFT power follows an exponential distribution with unit mean
 * (if the SFT data is distributed normally), while the second normalization is appropriate
 * in the time domain.  In either case, the mean of \f$ |\tilde{x}_k|^2 \f$ is
 * estimated using the median, suitably normalized assuming that the power is
 * distributed is exponentially.
 *
 * ### Uses ###
 *
 * \code
 * LALSFTtoPeriodogram ()
 * LALPeriodoToRngmed ()
 * LALNormalizeSFT ()
 * LALNormalizeSFTVect ()
 * LALNormalizeMultiSFTVect ()
 * \endcode
 *
 * The function LALNormalizeSFTVect() takes as input a vector of SFTs and normalizes
 * them.  This function calls the functions LALNormalizeSFT() which normalizes a
 * single SFT, LALSFTtoPeriodogram() which calculates the \f$ |\tilde{x}|^2 \f$ and
 * LALPeriodoToRngmed () which applies the running median algorithm to find a vector
 * of medians.  The function LALNormalizeMultiSFTVect() normalizes a multi-IFO collection
 * of SFT vectors and also returns a collection of power-estimates for these vectors using
 * the Running median method.
 *
 */

/**
 * Normalize an sft based on RngMed estimated PSD, and returns running-median.
 */
int
XLALNormalizeSFT ( REAL8FrequencySeries *rngmed, 	/**< [out] rng-median smoothed periodogram over SFT (Tsft*Sn/2) (must be allocated) */
                   SFTtype              *sft,		/**< SFT to be normalized */
                   UINT4                blockSize	/**< Running median block size for rngmed calculation */
                   )
{
  /* check input argments */
  XLAL_CHECK (sft && sft->data && sft->data->data && sft->data->length > 0,
              XLAL_EINVAL, "Invalid NULL or zero-length input in 'sft'" );

  XLAL_CHECK ( rngmed && rngmed->data && rngmed->data->data && rngmed->data->length > 0,
               XLAL_EINVAL, "Invalid NULL or zero-length input in 'rngmed'" );
  /* make sure there is no size mismatch */
  UINT4 length = sft->data->length;
  XLAL_CHECK ( length == rngmed->data->length, XLAL_EINVAL, "SFT length (%d) differs from rngmed length (%d)", length, rngmed->data->length );

  /* calculate the rngmed */
  XLAL_CHECK ( XLALSFTtoRngmed (rngmed, sft, blockSize) == XLAL_SUCCESS, XLAL_EFUNC, "XLALSFTtoRngmed() failed" );

  /* loop over sft and normalize */
  for (UINT4 j = 0; j < length; j++)
    {
      REAL8 Tsft_Sn_b2 = rngmed->data->data[j];		/* Wiener-Kinchine: E[|data|^2] = Tsft * Sn / 2 */
      REAL8 norm = 1.0 / sqrt(Tsft_Sn_b2);
      /* frequency domain normalization */
      sft->data->data[j] *= ((REAL4) norm);
    } // for j < length

  return XLAL_SUCCESS;

} /* XLALNormalizeSFT() */


/**
 * Function for normalizing a vector of SFTs.
 */
int
XLALNormalizeSFTVect ( SFTVector  *sftVect,	/**< [in/out] pointer to a vector of SFTs which will be normalized */
                       UINT4     blockSize	/**< Running median window size */
                       )
{
  /* check input argments */
  XLAL_CHECK ( sftVect && sftVect->data && sftVect->length > 0,
               XLAL_EINVAL, "Invalid NULL or zero-length input in 'sftVect'");

  /* memory allocation of rngmed using length of first sft
     -- assume all sfts have the same length*/
  UINT4 lengthsft = sftVect->data->data->length;

  /* allocate memory for a single rngmed */
  REAL8FrequencySeries *rngmed;
  XLAL_CHECK ( ( rngmed = XLALCalloc(1, sizeof(*rngmed))) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1,%d)", sizeof(*rngmed) );
  XLAL_CHECK ( ( rngmed->data = XLALCreateREAL8Vector ( lengthsft ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector ( %d ) failed.", lengthsft );

  /* loop over sfts and normalize them */
  for (UINT4 j = 0; j < sftVect->length; j++)
    {
      SFTtype *sft = &sftVect->data[j];

      /* call sft normalization function */
      XLAL_CHECK ( XLALNormalizeSFT ( rngmed, sft, blockSize ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALNormalizeSFT() failed." );

    } /* for j < sftVect->length */

  /* free memory for psd */
  XLALDestroyREAL8Vector ( rngmed->data );
  XLALFree(rngmed);

  return XLAL_SUCCESS;

} /* XLALNormalizeSFTVect() */


/**
 * Function for normalizing a multi vector of SFTs in a multi IFO search and
 * returns the running-median estimates of the power.
 */
MultiPSDVector *
XLALNormalizeMultiSFTVect ( MultiSFTVector *multsft,		/**< [in/out] multi-vector of SFTs which will be normalized */
                            UINT4          blockSize		/**< Running median window size */
                            )
{
  /* check input argments */
  XLAL_CHECK_NULL (multsft && multsft->data && multsft->length > 0,
              XLAL_EINVAL, "Invalid NULL or zero-length input 'multsft'");

  /* allocate multipsd structure */
  MultiPSDVector *multiPSD;
  XLAL_CHECK_NULL ( ( multiPSD = XLALCalloc (1, sizeof(*multiPSD))) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1, sizeof(*multiPSD))");

  UINT4 numifo = multsft->length;
  multiPSD->length = numifo;
  XLAL_CHECK_NULL ( ( multiPSD->data = XLALCalloc ( numifo, sizeof(*multiPSD->data))) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( %d, %d)", numifo, sizeof(*multiPSD->data) );

  /* loop over ifos */
  for ( UINT4 X = 0; X < numifo; X++)
    {
      UINT4 numsft = multsft->data[X]->length;

      /* allocation of psd vector over SFTs for this detector X */
      XLAL_CHECK_NULL ( (multiPSD->data[X] = XLALCalloc(1, sizeof(*multiPSD->data[X]))) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc(1, %d)", sizeof(*multiPSD->data[X]));

      multiPSD->data[X]->length = numsft;
      XLAL_CHECK_NULL ( (multiPSD->data[X]->data = XLALCalloc ( numsft, sizeof(*(multiPSD->data[X]->data)))) != NULL, XLAL_ENOMEM, "Failed to XLALCalloc ( %d, size)", numsft );

      /* loop over sfts for this IFO X */
      for (UINT4 j = 0; j < numsft; j++)
        {
          SFTtype *sft = &multsft->data[X]->data[j];

          /* memory allocation of psd vector for this SFT */
          UINT4 lengthsft = sft->data->length;
          XLAL_CHECK_NULL ( (multiPSD->data[X]->data[j].data = XLALCreateREAL8Vector ( lengthsft ) ) != NULL, XLAL_EFUNC, "XLALCreateREAL8Vector(%d) failed.", lengthsft );

          XLAL_CHECK_NULL( XLALNormalizeSFT ( &multiPSD->data[X]->data[j], sft, blockSize ) == XLAL_SUCCESS, XLAL_EFUNC, "XLALNormalizeSFT() failed");

        } /* for j < numsft */

    } /* for X < numifo */

  return multiPSD;

} /* XLALNormalizeMultiSFTVect() */


/**
 * Calculates a smoothed (running-median) periodogram for the given SFT.
 */
int
XLALSFTtoRngmed ( REAL8FrequencySeries *rngmed,	/**< [out] running-median smoothed periodo [must be allocated!] */
                  const SFTtype *sft,		/**< [in]  input SFT */
                  UINT4 blockSize		/**< Running median block size */
                  )
{
  /* check argments */
  XLAL_CHECK ( sft != NULL, XLAL_EINVAL, "Invalid NULL pointer passed in 'sft'" );
  XLAL_CHECK ( sft->data != NULL, XLAL_EINVAL, "Invalid NULL pointer in sft->data" );
  XLAL_CHECK ( sft->data->length > 0, XLAL_EINVAL, "Zero-length SFT passed (sft->data->length=0)" );
  XLAL_CHECK ( sft->data->data != NULL, XLAL_EINVAL, "Invalid NULL data in sft->data->data" );
  XLAL_CHECK ( rngmed != NULL, XLAL_EINVAL, "Invalid NULL pointer passed in 'rngmed'" );
  XLAL_CHECK ( rngmed->data != NULL, XLAL_EINVAL, "Invalid NULL pointer in rngmed->data" );
  XLAL_CHECK ( rngmed->data->length == sft->data->length, XLAL_EINVAL, "Allocated rngmed data-vector has to have same length (%d) as the SFT (%d)",
               rngmed->data->length, sft->data->length );
  XLAL_CHECK ( rngmed->data->data != NULL, XLAL_EINVAL, "Invalid NULL pointer in rngmed->data->data" );

  UINT4 length = sft->data->length;

  REAL8FrequencySeries periodo;
  XLAL_CHECK ( (periodo.data = XLALCreateREAL8Vector ( length )) != NULL, XLAL_EFUNC, "Failed to allocate periodo.data of length %d", length);

  /* calculate the periodogram */
  XLAL_CHECK ( XLALSFTtoPeriodogram ( &periodo, sft ) == XLAL_SUCCESS, XLAL_EFUNC, "Call to XLALSFTtoPeriodogram() failed.\n");

  /* calculate the rngmed */
  if ( blockSize > 0 )
    {
      XLAL_CHECK ( XLALPeriodoToRngmed ( rngmed, &periodo, blockSize ) == XLAL_SUCCESS, XLAL_EFUNC, "Call to XLALPeriodoToRngmed() failed." );
    }
  else	// blockSize==0 means don't use any running-median, just *copy* the periodogram contents into the output
    {
      strcpy ( rngmed->name, periodo.name );
      rngmed->epoch 	= periodo.epoch;
      rngmed->f0 	= periodo.f0;
      rngmed->deltaF 	= periodo.deltaF;
      rngmed->sampleUnits=periodo.sampleUnits;
      memcpy ( rngmed->data->data, periodo.data->data, periodo.data->length * sizeof(periodo.data->data[0]) );
    }

  /* free memory */
  XLALFree ( periodo.data->data );
  XLALFree ( periodo.data );

  return XLAL_SUCCESS;

} /* XLALSFTtoRngmed() */

/**
 * Calculate the "periodogram" of an SFT, ie the modulus-squares of the SFT-data.
 */
int
XLALSFTtoPeriodogram ( REAL8FrequencySeries    *periodo,	/**< [out] mod squares of SFT data (has to be allocated) */
                       const COMPLEX8FrequencySeries *SFT	/**< [in] input SFT */
                       )
{
  /* check input argments */
  XLAL_CHECK ( SFT != NULL, XLAL_EINVAL, "Invalid NULL pointer passed in 'SFT'" );
  XLAL_CHECK ( SFT->data != NULL, XLAL_EINVAL, "Invalid NULL pointer in SFT->data" );
  XLAL_CHECK ( SFT->data->length > 0, XLAL_EINVAL, "Zero-length SFT passed (SFT->data->length=0)" );
  XLAL_CHECK ( SFT->data->data != NULL, XLAL_EINVAL, "Invalid NULL data in SFT->data->data" );
  XLAL_CHECK ( periodo != NULL, XLAL_EINVAL, "Invalid NULL pointer passed in 'periodo'" );
  XLAL_CHECK ( periodo->data != NULL, XLAL_EINVAL, "Invalid NULL pointer in periodo->data" );
  XLAL_CHECK ( periodo->data->length == SFT->data->length, XLAL_EINVAL, "Allocated periodo data-vector has to have same length (%d) as the SFT (%d)",
               periodo->data->length, SFT->data->length );
  XLAL_CHECK ( periodo->data->data != NULL, XLAL_EINVAL, "Invalid NULL pointer in periodo->data->data" );

  /* copy SFT header */
  strcpy ( periodo->name, SFT->name );
  periodo->epoch = SFT->epoch;
  periodo->f0 = SFT->f0;
  periodo->deltaF = SFT->deltaF;

  /* check lengths are same */
  UINT4 length = SFT->data->length;
  REAL8 *out = periodo->data->data;
  COMPLEX8 *in = SFT->data->data;

  for (UINT4 j=0; j<length; j++)
    {
      /* extra-paranoia: make absolutely sure that the calculation below is in REAL8
       * in order to avoid underflow-problems (data 'in' can be of order ~ 1e-20 )
       */
      *out = ((REAL8)crealf(*in))*((REAL8)crealf(*in)) + ((REAL8)cimagf(*in))*((REAL8)cimagf(*in));
      ++out;
      ++in;
    } // for j<length

  return XLAL_SUCCESS;

} /* XLALSFTtoPeriodogram() */

/**
 * Calculates running median over a single periodogram.
 */
int
XLALPeriodoToRngmed ( REAL8FrequencySeries  *rngmed,		/**< [out] resulting 'smoothed' periodogram (must be allocated) */
                      const REAL8FrequencySeries  *periodo,	/**< [in] input periodogram */
                      UINT4 blockSize				/**< Running median block size */
                      )
{
  /* check input argments are not NULL */
  XLAL_CHECK ( periodo != NULL && periodo->data != NULL && periodo->data->data && periodo->data->length > 0,
               XLAL_EINVAL, "Invalid input 'periodo': needs to be allocated and non-zero length" );
  XLAL_CHECK ( rngmed != NULL && rngmed->data != NULL &&  rngmed->data->data != NULL && rngmed->data->length > 0,
               XLAL_EINVAL, "Invalid input 'rngmend': needs to be allocated and non-zero length" );
  UINT4 length = periodo->data->length;
  XLAL_CHECK ( length == rngmed->data->length,
               XLAL_EINVAL, "'periodo' vector must be same length (%d) as 'rngmed' vector (%d)", periodo->data->length, rngmed->data->length );
  XLAL_CHECK ( blockSize > 0, XLAL_EINVAL, "'blockSize = %d' must be > 0", blockSize );
  XLAL_CHECK( length >= blockSize, XLAL_EINVAL, "Need at least %d bins in SFT (have %d) to perform running median!\n", blockSize, length );

  /* copy periodogram header */
  strcpy ( rngmed->name, periodo->name );
  rngmed->epoch = periodo->epoch;
  rngmed->f0 = periodo->f0;
  rngmed->deltaF = periodo->deltaF;

  UINT4 blocks2 = blockSize/2; /* integer division, round down */

  LALRunningMedianPar rngMedPar;
  rngMedPar.blocksize = blockSize;
  REAL8Sequence mediansV, inputV;
  inputV.length = length;
  inputV.data = periodo->data->data;
  UINT4 medianVLength = length - blockSize + 1;
  mediansV.length = medianVLength;
  mediansV.data = rngmed->data->data + blocks2;

  LALStatus status = empty_LALStatus;
  LALDRunningMedian2 ( &status, &mediansV, &inputV, rngMedPar);
  XLAL_CHECK ( status.statusCode == 0, XLAL_EFAILED, "LALDRunningMedian2() failed with statusCode = %d", status.statusCode );

  /* copy values in the wings */
  for ( UINT4 j=0; j<blocks2; j++)
    rngmed->data->data[j] = rngmed->data->data [ blocks2 ];

  for (UINT4 j=blocks2 + medianVLength; j<length; j++)
    rngmed->data->data[j] = rngmed->data->data [ blocks2 + medianVLength - 1 ];

  /* get the bias factor -- for estimating the mean from the median */
  REAL8 medianBias = XLALRngMedBias ( blockSize );
  XLAL_CHECK ( xlalErrno == 0, XLAL_EFUNC, "XLALRngMedBias() failed");

  /* normalize by the bias factor */
  REAL8 medianBiasInv = 1.0 / medianBias;
  for (UINT4 j=0; j<length; j++)
    rngmed->data->data[j] *= medianBiasInv;

  return XLAL_SUCCESS;

} /* XLALPeriodoToRngmed() */


/**
 * Calculate the cross-correlation periodogram from 2 SFTs.
 */
int
XLALSFTstoCrossPeriodogram ( REAL8FrequencySeries *periodo,		/**< [out] modulus square of SFT data (must be allocated) */
                             const COMPLEX8FrequencySeries *sft1,	/**< [in] pointer to first SFT  */
                             const COMPLEX8FrequencySeries *sft2	/**< [in] pointer to second SFT */
			   )
{
  /* check input argments */
  XLAL_CHECK ( periodo && periodo->data && periodo->data->data && periodo->data->length > 0,
               XLAL_EINVAL, "Invalid NULL or zero-length input 'periodo'");

  XLAL_CHECK ( sft1 && sft1->data && sft1->data->data && sft1->data->length > 0,
               XLAL_EINVAL, "Invalud NULL or zero-length input 'sft1'");

  XLAL_CHECK ( sft2 && sft2->data && sft2->data->data && sft2->data->length > 0,
               XLAL_EINVAL, "Invalud NULL or zero-length input 'sft2'");

  /* make sure both sfts are consistent in frequency and freq. band
     -- time stamps need not be consistent */
  XLAL_CHECK ( sft2->data->length == sft1->data->length, XLAL_EINVAL, "SFT lengths differ len1 = %d, len2 = %d", sft1->data->length, sft2->data->length );
  XLAL_CHECK ( sft2->f0 == sft1->f0, XLAL_EINVAL, "SFT start-frequencies differ f0_1 = %g, f0_2 = %g", sft1->f0, sft2->f0 );
  XLAL_CHECK ( sft2->deltaF == sft1->deltaF, XLAL_EINVAL, "SFT frequency spacings differ deltaF_1 = %g, deltaF_2 = %g", sft1->deltaF, sft2->deltaF );

  /* copy values from SFT */
  /*   periodo->epoch.gpsSeconds = sft1->epoch.gpsSeconds; */
  /*   periodo->epoch.gpsNanoSeconds = sft1->epoch.gpsNanoSeconds; */
  periodo->f0 = sft1->f0;
  periodo->deltaF = sft1->deltaF;

  /* check lengths are same */
  UINT4 length = sft1->data->length;
  XLAL_CHECK ( length == periodo->data->length, XLAL_EINVAL, "SFT length (%d) differs from periodo length (%d)", length, periodo->data->length );

  REAL8    *out = periodo->data->data;
  COMPLEX8 *in1 = sft1->data->data;
  COMPLEX8 *in2 = sft2->data->data;

  for (UINT4 j=0; j<length; j++)
    {
      /* extra-paranoia: make absolutely sure that the calculation below is in REAL8
       * in order to avoid underflow-problems (data 'in' can be of order ~ 1e-20 )
       */
      *out = ((REAL8)crealf(*in1))*((REAL8)crealf(*in2)) + ((REAL8)cimagf(*in1))*((REAL8)cimagf(*in2));
      ++out;
      ++in1;
      ++in2;
    } // for j < length

  return XLAL_SUCCESS;

} /* XLALSFTstoCrossPeriodogram() */


// ****************************** OBSOLETE + DEPRECATED LAL-INTERFACE FUNCTIONS ******************************

/**
 * \deprecated use XLALSFTtoPeriodogram() instead
 * Calculate the "periodogram" of an SFT, ie the modulus-squares of the SFT-data.
 */
void
LALSFTtoPeriodogram (LALStatus    *status,		/**< pointer to LALStatus structure */
		     REAL8FrequencySeries    *periodo,	/**< [out] mod squares of SFT data */
		     const COMPLEX8FrequencySeries *SFT	/**< [in] input SFT */
		     )
{
  INITSTATUS(status);

  if ( XLALSFTtoPeriodogram ( periodo, SFT ) != XLAL_SUCCESS )
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL );

  /* normal exit */
  RETURN (status);

} /* end LALSFTtoPeriodogram() */

/**
 * \deprecated use XLALPeriodoToRngmed() instead
 */
void
LALPeriodoToRngmed (LALStatus  *status,				/**< pointer to LALStatus structure */
		    REAL8FrequencySeries  *rngmed,		/**< [out] resulting 'smoothed' periodogram */
		    const REAL8FrequencySeries  *periodo,	/**< [in] input periodogram */
		    UINT4 blockSize				/**< Running median block size */
		    )
{
  INITSTATUS(status);

  if ( XLALPeriodoToRngmed ( rngmed, periodo, blockSize ) != XLAL_SUCCESS )
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL );

  /* normal exit */
  RETURN (status);

} /* LALPeriodoToRngmed() */


/**
 * \deprecated use XLALSFTtoRngmed() instead
 */
void
LALSFTtoRngmed (LALStatus  *status,		/**< pointer to LALStatus structure */
		REAL8FrequencySeries  *rngmed,		/**< [out] running-median smoothed periodo [must be allocated!] */
		const COMPLEX8FrequencySeries *sft,	/**< [in]  input SFT */
		UINT4 blockSize				/**< Running median block size */
		)
{
  INITSTATUS(status);

  if ( XLALSFTtoRngmed ( rngmed, sft, blockSize ) != XLAL_SUCCESS )
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL );

  /* normal exit */
  RETURN (status);

} /* LALSFTtoRngmed() */


/**
 * \deprecated use XLALNormalizeSFT() instead
 */
void
LALNormalizeSFT (LALStatus           *status,		/**< pointer to LALStatus structure */
		 REAL8FrequencySeries *rngmed, 	/**< [out] rng-median smoothed periodogram over SFT (Tsft*Sn/2) */
		 SFTtype              *sft,     /**< SFT to be normalized */
		 UINT4                blockSize)/**< Running median block size for rngmed calculation */
{
  INITSTATUS(status);

  if ( XLALNormalizeSFT ( rngmed, sft, blockSize) != XLAL_SUCCESS )
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL );

  /* normal exit */
  RETURN (status);

} /* LALNormalizeSFT() */


/**
 * Function for normalizing a vector of SFTs.
 */
void
LALNormalizeSFTVect (LALStatus  *status,		/**< pointer to LALStatus structure */
		     SFTVector  *sftVect,	/**< [in/out] pointer to a vector of SFTs which will be normalized */
		     UINT4     blockSize	/**< Running median window size */
		     )
{
  INITSTATUS(status);

  if ( XLALNormalizeSFTVect ( sftVect, blockSize ) != XLAL_SUCCESS )
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL );

  /* normal exit */
  RETURN (status);

} /* LALNormalizeSFTVect() */


/**
 * Function for normalizing a multi vector of SFTs in a multi IFO search and also
 * returns the running-median estimates of the power.
 */
void
LALNormalizeMultiSFTVect (LALStatus      *status,		/**< pointer to LALStatus structure */
			  MultiPSDVector **multiRngmed,	/**< [out] multi running-median power estimates of input SFTs */
			  MultiSFTVector *multsft,	/**< [in/out] multi-vector of SFTs which will be normalized */
			  UINT4          blockSize	/**< Running median window size */
			  )
{
  INITSTATUS(status);

  if ( multiRngmed == NULL ) ABORT ( status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL );
  if ( (*multiRngmed) != NULL ) ABORT ( status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  MultiPSDVector *ret;
  if ( (ret = XLALNormalizeMultiSFTVect ( multsft, blockSize )) == NULL )
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL );

  (*multiRngmed) = ret;

  /* normal exit */
  RETURN (status);

} /* LALNormalizeMultiSFTVect() */


/**
 * Calculate the cross-correlation periodogram from 2 SFTs.
 */
void
LALSFTstoCrossPeriodogram (LALStatus    *status,		/**< pointer to LALStatus structure */
			   REAL8FrequencySeries *periodo,	/**< [out] modulus square of SFT data  */
			   const COMPLEX8FrequencySeries *sft1,	/**< [in] pointer to first SFT  */
			   const COMPLEX8FrequencySeries *sft2	/**< [in] pointer to second SFT */
			   )
{
  INITSTATUS(status);

  if ( XLALSFTstoCrossPeriodogram ( periodo, sft1, sft2 ) != XLAL_SUCCESS )
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL );

  /* normal exit */
  RETURN (status);

} /* LALSFTstoCrossPeriodogram() */
