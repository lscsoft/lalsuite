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
LALPeriodoToPSDRngMed ()
LALNormalizeSFT ()
LALNormalizeSFTVect ()
LALNormalizeMultiSFTVect ()
\endcode				      


The function LALNormalizeSFTVect() takes as input a vector of SFTs and normalizes
them.  This function calls the functions LALNormalizeSFT() which normalizes a
single SFT, LALSFTtoPeriodogram() which calculates the \f$ |\tilde{x}|^2 \f$ and 
LALPeriodoToPSDRngMed() which applies the running median algorithm to find a vector
of medians.  The function LALNormalizeMultiSFTVect() normalizes a multi-ifo collection 
of SFT vectors and also returns a collection of PSD estimates for these vectors using
the Running median method. 

*/


/*----------- laldoc documentation -------------- */

/************************************ <lalVerbatim file="NormalizeSFTRngMedCV">
Author: Krishnan, B.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{Module \texttt{NormalizeSFTRngMed.c}}
\label{ss:NormalizeSFTRngMed.c}


Normalizing SFTs using the Running median

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{NormalizeSFTRngMedD}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

Normalizing SFTs using the Running Median


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALHO()
\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{NormalizeSFTRngMedCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*********************************************** </lalLaTeX> */


#include <lal/NormalizeSFTRngMed.h>



NRCSID (NORMALIZESFTRNGMEDC, "$Id$");


/********************************* <lalVerbatim file="NormalizeSFTRngMedD"> */
/** \brief Calculate the modulus square of a single SFT 
    \param *SFT : pointer to a SFT
    \param *periodo : pointer to REAL8FrequencySeries containing modulus square of SFT data 
*/
void LALSFTtoPeriodogram (LALStatus    *status,
			  REAL8FrequencySeries    *periodo,
			  const COMPLEX8FrequencySeries *SFT)
{/*   *********************************************  </lalVerbatim> */


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




/* *******************************  <lalVerbatim file="NormalizeSFTRngMedD"> */
/** \brief  Calculates running psd for a single periodogram using the running median
    \param  *periodo : pointer to REAL8FrequencySeries containing square modulus of a SFT
    \param blockSize : Running median block size
    \param *psd : pointer to REAL8FrequencySeries containing psd estimate

*/
void LALPeriodoToPSDRngMed (LALStatus  *status,
			    REAL8FrequencySeries  *psd,
			    const REAL8FrequencySeries  *periodo,
			    UINT4                  blockSize)
{/*   *********************************************  </lalVerbatim> */
  UINT4 blocks2;
  UINT4 j;
  UINT4 length;
  LALRunningMedianPar rngMedPar;
  REAL8Sequence mediansV, inputV;
  REAL8 medianBias;

  INITSTATUS (status, "LALPeriodoToPSDRngMed", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL */
  ASSERT (periodo, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (periodo->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (periodo->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 
  ASSERT (periodo->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (psd, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);  
  ASSERT (psd->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (psd->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);  
  ASSERT (psd->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (blockSize > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 


  /* copy values from the periodogram */
  psd->epoch.gpsSeconds = periodo->epoch.gpsSeconds;
  psd->epoch.gpsNanoSeconds = periodo->epoch.gpsNanoSeconds;
  psd->f0 = periodo->f0;
  psd->deltaF = periodo->deltaF;

  /* check lengths are same */
  length = periodo->data->length;
  ASSERT (length == psd->data->length, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);  
  ASSERT (length > blockSize, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);  

  blocks2 = blockSize/2; /* integer division */

  rngMedPar.blocksize = (UINT4)blockSize;
  inputV.length = length;
  inputV.data = periodo->data->data;
  mediansV.length= length - blockSize + 1;
  mediansV.data = psd->data->data + blocks2;    

  TRY( LALDRunningMedian2(status->statusPtr, &mediansV, &inputV, rngMedPar), status);

  /* copy values in the wings */
  for (j=0; j<blocks2; j++)
    psd->data->data[j] = psd->data->data[blocks2];

  for (j=blocks2+length-blockSize+1; j<length; j++)
    psd->data->data[j] = psd->data->data[blocks2 + length-blockSize];

  /* get the bias factor -- for estimating the mean from the median */
  TRY ( LALRngMedBias( status->statusPtr, &medianBias, blockSize ), status);

  /* normalize by the bias factor */
  for (j=0; j<length; j++)
    psd->data->data[j] /= medianBias;

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* end LALPeriodoToPSDRngMed() */




/* *******************************  <lalVerbatim file="NormalizeSFTRngMedD"> */
/** \brief  Calculates running psd for a single sft using the running median
    \param  *sft : pointer to COMPLEX8FrequencySeries (the SFT data)
    \param blockSize : Running median block size
    \param *psd : pointer to correctly allocated REAL8FrequencySeries containing psd estimate 
*/
void LALSFTtoPSDRngMed (LALStatus  *status,
			REAL8FrequencySeries  *psd,
			const COMPLEX8FrequencySeries *sft,
			UINT4                  blockSize)
{/*   *********************************************  </lalVerbatim> */

  REAL8FrequencySeries periodo;
  INT4 length;

  INITSTATUS (status, "LALSFTtoPSDRngMed", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL */
  ASSERT (sft, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (sft->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (sft->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 
  ASSERT (sft->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (psd, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);  
  ASSERT (psd->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (psd->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);  
  ASSERT (psd->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (blockSize > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 

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

  /* calculate the psd */
  LALPeriodoToPSDRngMed (status->statusPtr, psd, &periodo, blockSize);
  BEGINFAIL (status) { 
    LALFree (periodo.data->data);
    LALFree (periodo.data);
  } ENDFAIL (status);
  
  /* free memory */
  LALFree(periodo.data->data);
  LALFree(periodo.data);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* end LALSFTtoPSDRngMed() */




/* *******************************  <lalVerbatim file="NormalizeSFTRngMedD"> */
/** \brief Normalizes a sft based on RngMed 
    \param [out] psd estimate from sft using running median
    \param *sft : pointer to a SFT which will be normalized
    \param blockSize : Running median block size
*/
void LALNormalizeSFT (LALStatus            *status,
 		      REAL8FrequencySeries *psd,     /**< [out] psd estimate from sft using running median */
		      SFTtype              *sft,      /**< SFT to be normalized */
		      UINT4                blockSize) /**< Running median block size for psd calculation */ 
{/*   *********************************************  </lalVerbatim> */

  UINT4 j;
  REAL8 Sn;

  INITSTATUS (status, "LALNormalizeSFT", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL and other sanity checks*/
  ASSERT (sft, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (sft->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (sft->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 
  ASSERT (sft->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 

  ASSERT (psd, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (psd->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (psd->data->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 
  ASSERT (psd->data->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 

  ASSERT (blockSize > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 

  /* make sure there is no size mismatch */
  if ( psd->data->length != sft->data->length ) {
    ABORT ( status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 
  }
  
  /* calculate the psd */
  TRY (LALSFTtoPSDRngMed (status->statusPtr, psd, sft, blockSize), status);

  /* loop over sft and normalize */
  for (j = 0; j < sft->data->length; j++) {
  
    Sn = psd->data->data[j]; 
    
    /* frequency domain normalization */
    sft->data->data[j].re /= sqrt(Sn);
    sft->data->data[j].im /= sqrt(Sn);
  }


  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
} /* end LALNormalizeSFT() */




/* *******************************  <lalVerbatim file="NormalizeSFTRngMedD"> */
/** \brief Function for normalizing a vector of SFTs
    \param *sftVect  pointer to a vector of SFTs which will be normalized
    \param blockSize : Running median window size
*/
void LALNormalizeSFTVect (LALStatus  *status,
			  SFTVector  *sftVect,
			  UINT4     blockSize)
{/*   *********************************************  </lalVerbatim> */
  /* normalizes a sft vector using RngMed */
  INT4 j, lengthsft;
  REAL8FrequencySeries *psd = NULL;

  INITSTATUS (status, "LALNormalizeSFT", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL and other sanity checks*/
  ASSERT (sftVect, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (sftVect->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (sftVect->length > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 
  ASSERT (blockSize > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 

  /* memory allocation of psd using length of first sft 
     -- assume all sfts have the same length*/
  lengthsft = sftVect->data->data->length;

  /* allocate memory for a single psd */
  if ( (psd = (REAL8FrequencySeries *)LALCalloc(1, sizeof(REAL8FrequencySeries))) == NULL){
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  psd->data = NULL;
  if ( (psd->data = (REAL8Sequence *)LALCalloc(1, sizeof(REAL8Sequence))) == NULL) {
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  psd->data->length = lengthsft;
  if ( (psd->data->data = (REAL8 *)LALCalloc( lengthsft, sizeof(REAL8))) == NULL) {
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
    LALNormalizeSFT (status->statusPtr, psd, sft, blockSize);
    BEGINFAIL (status) { 
      LALFree (psd->data->data);
      LALFree (psd->data);
      LALFree (psd);
    } ENDFAIL (status);

  } /* for loop over sfts */

  /* free memory for psd */
  LALFree(psd->data->data);
  LALFree(psd->data); 
  LALFree(psd); 

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);

} /* end LALNormalizeSFTVect() */



/** \brief Function for normalizing a multi vector of SFTs in a multi IFO search and also  calculates the PSDs
    \param [out] **psdvect multi vector of PSD estimates of input SFTs using the running median
    \param *MultiSFTVect  pointer to a vector of SFTs which will be normalized
    \param [in] blockSize : Running median window size
*/
void LALNormalizeMultiSFTVect (LALStatus      *status,
			       MultiPSDVector **out,
			       MultiSFTVector *multsft,
			       UINT4          blockSize)
{

  UINT4 k, j; /* k loops over IFOs and j over SFTs for each IFO */
  UINT4 jCleanUp, kCleanUp; /* indices used in clean-up loops */
  UINT4 numifo, numsft;
  MultiPSDVector *multpsd = NULL;

  INITSTATUS (status, "LALNormalizeMultiSFT", NORMALIZESFTRNGMEDC);
  ATTATCHSTATUSPTR (status);

  /* check argments are not NULL and other sanity checks*/
  ASSERT (multsft, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT (multsft->length, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL);
  ASSERT (multsft->data, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL);

  ASSERT (out, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 
  ASSERT ( *out == NULL, status, NORMALIZESFTRNGMEDH_ENULL, NORMALIZESFTRNGMEDH_MSGENULL); 

  ASSERT (blockSize > 0, status, NORMALIZESFTRNGMEDH_EVAL, NORMALIZESFTRNGMEDH_MSGEVAL); 

  /* first memory allocation for multipsd structure */
  if ( (multpsd = (MultiPSDVector *)LALCalloc(1, sizeof(MultiPSDVector))) == NULL) {
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  multpsd->length = numifo = multsft->length;
  if ( (multpsd->data = (PSDVector **)LALCalloc( numifo, sizeof(PSDVector *))) == NULL) {
    ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
  }

  /* loop over ifos */
  for ( k = 0; k < numifo; k++) {
   
    /* second memory allocation for psd vector */
    if ( (multpsd->data[k] = (PSDVector *)LALCalloc(1, sizeof(PSDVector))) == NULL) {
      ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
    }

    multpsd->data[k]->length = numsft = multsft->data[k]->length;
    if ( (multpsd->data[k]->data = (REAL8FrequencySeries *)LALCalloc(numsft, sizeof(REAL8FrequencySeries))) == NULL) {
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
      multpsd->data[k]->data[j].data = NULL;
      if ( (multpsd->data[k]->data[j].data = (REAL8Sequence *)LALCalloc(1, sizeof(REAL8Sequence))) == NULL) {
	ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
      }

      multpsd->data[k]->data[j].data->length = lengthsft = sft->data->length;
      if ( (multpsd->data[k]->data[j].data->data = (REAL8 *)LALCalloc( lengthsft, sizeof(REAL8))) == NULL) {
	ABORT( status, NORMALIZESFTRNGMEDH_EMEM, NORMALIZESFTRNGMEDH_MSGEMEM);
      }

      LALNormalizeSFT (status->statusPtr, multpsd->data[k]->data + j, sft, blockSize);
      BEGINFAIL (status) { 
        /* clean up for this value of k */
	for ( jCleanUp = 0; jCleanUp < j; jCleanUp++) {
	  LALFree( multpsd->data[k]->data[jCleanUp].data->data);
	  LALFree( multpsd->data[k]->data[jCleanUp].data);
	}
	LALFree( multpsd->data[k]->data);
	LALFree( multpsd->data[k]);
	/* clean up for previous values of k */
	for ( kCleanUp = 0; kCleanUp < k-1; kCleanUp++) {
	  for ( jCleanUp = 0; jCleanUp < multsft->data[kCleanUp]->length; jCleanUp++) {
	    LALFree( multpsd->data[kCleanUp]->data[jCleanUp].data->data);
	    LALFree( multpsd->data[kCleanUp]->data[jCleanUp].data);
	  }
	  LALFree( multpsd->data[kCleanUp]->data);
	  LALFree( multpsd->data[kCleanUp]);
	}
	/* clean up memory allocated outside loop */
	LALFree(multpsd->data);
	LALFree(multpsd);
      } ENDFAIL (status);          
    }
  }

  *out = multpsd;

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
} /* LALNormalizeMultiSFTVect() */















