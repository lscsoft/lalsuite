/*
 *  Copyright (C) 2005 Badri Krishnan, Alicia Sintes  
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
 * \file PeakSelect.c
 * \author Sintes, A.M. and Krishnan, B.
 * \brief Functions for selecting peaks from SFTs.
 *
 * History:   Created by Sintes May 21, 2003
 * Modified by Krishnan Oct, 2005
 *
 * \heading{\ref PeakSelect.c}
 * Routines for reading SFT binary files
 *
 * \heading{Prototypes}
 * <tt>PeakSelectHeader1()</tt>
 *
 * \heading{Description}
 * the output of the periodogram should be properly normalized !!!
 *
 * \heading{Uses}
 * \code
 * LALHO()
 * \endcode
 */

#include "./PeakSelect.h"

/*
 * The functions that make up the guts of this module
 */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALComputeMeanPower (LALStatus  *status,
			  REAL8                *mean,
			  REAL8Periodogram1    *peri)
{ 
  
  INT4       length;
  REAL8      sum;
  
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (peri,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (mean,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  
  length = peri->length;
  sum = 0.0;
  
  if (length > 0){
    REAL8      *in1;
    INT4      n;
    ASSERT (peri->data, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    in1= peri->data;
    n= length;
    while (n-- >0){
      sum += *in1;
      ++in1;
    }
    sum = sum/length;
  }
  
  *mean = sum;
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALSelectPeakWhiteNoise(LALStatus  *status,
			     UCHARPeakGram        *pg,
			     REAL8                *thr,
			     REAL8Periodogram1    *peri)
{ 
  

  INT4  length;
  INT4  nPeaks;
  
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (peri, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (pg,   status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (thr,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);

  pg->epoch.gpsSeconds     = peri->epoch.gpsSeconds;
  pg->epoch.gpsNanoSeconds = peri->epoch.gpsNanoSeconds;
  pg->timeBase     = peri->timeBase;
  pg->fminBinIndex = peri->fminBinIndex;
  
  length = peri->length;
  nPeaks = 0;
  
  ASSERT (length==pg->length,  status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  
  if (length > 0){
    UCHAR  *out;
    REAL8  *in1;
    REAL8  threshold;
    INT4      n;

    ASSERT (peri->data, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    ASSERT (pg->data,   status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    out = pg->data;
    in1 = peri->data;
    threshold = *thr;
    n = length;
    
    while (n-- >0){
      if ( *in1 > threshold){
        *out = 1 ;
	nPeaks++ ;
      } else {
        *out = 0;
      }     
      ++out;
      ++in1;    
    }         
  }
  
  pg->nPeaks = nPeaks;
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALUCHAR2HOUGHPeak(LALStatus  *status,
			HOUGHPeakGram        *pgOut,
			UCHARPeakGram        *pgIn)
{ 
  
  INT4  length;
  UINT4  nPeaks;
  
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (pgIn,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (pgOut, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);

  length = pgIn->length;
  pgOut->deltaF = 1./pgIn->timeBase;
  pgOut->fBinIni = pgIn->fminBinIndex;
  pgOut->fBinFin = pgOut->fBinIni + length -1;
  
  nPeaks = pgIn->nPeaks;
  ASSERT (nPeaks==pgOut->length,  status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  
  if (nPeaks >0){
    INT4    *peak;
    UCHAR   *in1;
    INT4    i;
    
    ASSERT (pgOut->peak,   status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    ASSERT (pgIn->data,   status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    peak = pgOut->peak;
    in1 = pgIn->data;
    
    for (i=0; i < length; i++){
      if (*in1){
        *peak = i;
	++peak;
      }
      ++in1;
    }
  }
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */


/******************************************************************/
/*  estimate psd from periodogram using runing-median  */
/*  this is just a wrapper to Mohanty's rngmed */
/******************************************************************/

void LALPeriodo2PSDrng (LALStatus  *status,
                     REAL8Periodogram1    *psd,
                     REAL8Periodogram1    *peri,
                     INT4                *blocksRNG)
{ 

  INT4   length, j;
  UINT2  blockSize, nblocks2;
  REAL8  *medians;
  LALRunningMedianPar rngmedPar;
  REAL8Sequence mediansV;
  REAL8Sequence inputV;
  
  
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);

  /*   Make sure the arguments are not NULL: */
  ASSERT (psd,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (peri, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (blocksRNG, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);

  psd->epoch.gpsSeconds     = peri->epoch.gpsSeconds;
  psd->epoch.gpsNanoSeconds = peri->epoch.gpsNanoSeconds;
  psd->timeBase     = peri->timeBase;
  psd->fminBinIndex = peri->fminBinIndex;

  length = peri->length;
  blockSize = *blocksRNG;

  ASSERT (length==psd->length,  status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  ASSERT (blockSize,  status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  ASSERT (blockSize<=length, status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);

  if (length > 0){
    ASSERT (peri->data,   status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    ASSERT (psd->data,   status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);

    nblocks2 = floor(blockSize/2.0);
    
    medians = psd->data+nblocks2; 

    rngmedPar.blocksize = (UINT4)blockSize;      
    inputV.length = length;
    inputV.data = peri->data;
    mediansV.length= length-blockSize+1;
    mediansV.data = psd->data+nblocks2;    
   
    /* calling Mohanty's function. It is not LAL compliant */
    /*rngmed(data, length, blockSize, medians); */
    /* this now calls Bernd's functions since other one is buggy */
    TRY( LALDRunningMedian2(status->statusPtr,
    			 &mediansV, &inputV, rngmedPar), status);
    

    for (j=0; j<nblocks2 ; ++j){psd->data[j] = medians[0]; }
    for (j=nblocks2+length-blockSize+1; j<length ; ++j){
      psd->data[j] = medians[length-blockSize];
    }
  }

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */

void LALSelectPeakColorNoise(LALStatus  *status,
			     UCHARPeakGram        *pg,
			     REAL8                *thr,
			     REAL8PeriodoPSD      *in)
{  
  
  INT4  length;
  INT4  nPeaks;
  
  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (pg,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (thr, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (in,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  
  /* check that the psd and periodogram are consistent */
  ASSERT(in->psd.timeBase==in->periodogram.timeBase,status,  PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  ASSERT(in->psd.fminBinIndex==in->periodogram.fminBinIndex,status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  ASSERT(in->psd.length==in->periodogram.length, status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  ASSERT(in->psd.epoch.gpsSeconds==in->periodogram.epoch.gpsSeconds, status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  ASSERT(in->psd.epoch.gpsNanoSeconds==in->periodogram.epoch.gpsNanoSeconds, status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);


  pg->epoch.gpsSeconds     = in->psd.epoch.gpsSeconds; 
  pg->epoch.gpsNanoSeconds = in->psd.epoch.gpsNanoSeconds;
  pg->timeBase     = in->psd.timeBase;
  pg->fminBinIndex = in->psd.fminBinIndex;
  
  length = in->psd.length;

  /* make sure that the length of the peakgram is consistent */
  ASSERT (length==pg->length,  status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);

  nPeaks = 0;

  if (length > 0){
    UCHAR  *out;
    REAL8  *psd; 
    REAL8  *peri;
    REAL8  threshold;
    INT4   n;
    
    ASSERT (in->periodogram.data, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    ASSERT (in->psd.data, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    ASSERT (pg->data,   status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    
    out = pg->data;
    peri = in->periodogram.data;
    psd  = in->psd.data;
    
    threshold = *thr;
    n = length;
    
    while (n-- >0){
      if ( *peri > threshold * (*psd) ){
        *out = 1 ;
	nPeaks++ ;
      } else {
        *out = 0;
      }     
      ++out;
      ++psd;
      ++peri;
    }
    
  }
  
  pg->nPeaks = nPeaks;
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}
/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */



/**
 * Convert a normalized sft into a peakgram by selecting bins in
 * which power exceeds a threshold.
 */

void SFTtoUCHARPeakGram(LALStatus        *status,
			UCHARPeakGram    *pg,
			const SFTtype    *sft,   
			REAL8            thr)
{  
  
  INT4  length;
  INT4  nPeaks;
  REAL8 f0, deltaF;
  REAL8FrequencySeries periodo;

  /* --------------------------------------------- */
  INITSTATUS(status);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (pg,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (sft, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  

  pg->epoch.gpsSeconds     = sft->epoch.gpsSeconds; 
  pg->epoch.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  deltaF = sft->deltaF;
  f0 = sft->f0;
  pg->timeBase     = 1.0/deltaF;
  pg->fminBinIndex = floor( f0 /deltaF + 0.5);
  length = sft->data->length;

  /* make sure that the length of the peakgram is consistent */
  ASSERT (length==pg->length,  status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);

  /* allocate memory for normalized periodogram */
  periodo.data = NULL;
  periodo.data = (REAL8Sequence *)LALMalloc(sizeof(REAL8Sequence));
  periodo.data->length = length;
  periodo.data->data = (REAL8 *)LALMalloc( length * sizeof(REAL8));

  TRY( LALSFTtoPeriodogram( status->statusPtr, &periodo, sft), status);

  nPeaks = 0;
  if (length > 0){
    UCHAR  *out;
    REAL8  *peri; 
    INT4   n;
    
    ASSERT (pg->data,   status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
    
    out = pg->data;
    peri = periodo.data->data;
    
    n = length;
    
    while (n-- >0){
      if ( *peri > thr ){
        *out = 1 ;
	nPeaks++ ;
      } else {
        *out = 0;
      }     
      ++out;
      ++peri;
    } /* end loop over n */
  }/* end if statement */ 
  
  pg->nPeaks = nPeaks;

  LALFree(periodo.data->data);
  LALFree(periodo.data);
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}
