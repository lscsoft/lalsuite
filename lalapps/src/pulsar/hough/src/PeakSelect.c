/*-----------------------------------------------------------------------
 *
 * File Name: PeakSelect.c
 *
 * Authors: Sintes, A.M.,  
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes May 21, 2003
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */

/************************************ <lalVerbatim file="PeakSelectCV">
Author: Sintes, A.M.,
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{MOdule \texttt{PeakSelect.c}}
\label{ss:PeakSelect.c}
Routines for reading SFT binary files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{PeakSelectD}
\index{\verb&PeakSelectHeader1()&}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

the output of the periodogram should be properly normalized !!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALHO()
\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{PeakSelectCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*********************************************** </lalLaTeX> */


#include "./PeakSelect.h"

NRCSID (PEAKSELECTC, "$Id$");

/*
 * The functions that make up the guts of this module
 */

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="PeakSelectD"> */
void LALComputeMeanPower (LALStatus  *status,
			  REAL8                *mean,
			  REAL8Periodogram1    *peri)
{ /*   *********************************************  </lalVerbatim> */
  
  INT4       length;
  REAL8      sum;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "LALComputeMeanPower", PEAKSELECTC);
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
/* *******************************  <lalVerbatim file="PeakSelectD"> */
void LALSelectPeakWhiteNoise(LALStatus  *status,
			     UCHARPeakGram        *pg,
			     REAL8                *thr,
			     REAL8Periodogram1    *peri)
{ /*   *********************************************  </lalVerbatim> */
  

  INT4  length;
  INT4  nPeaks;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "LALSelectPeakWhiteNoise", PEAKSELECTC);
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
/* *******************************  <lalVerbatim file="PeakSelectD"> */
void LALUCHAR2HOUGHPeak(LALStatus  *status,
			HOUGHPeakGram        *pgOut,
			UCHARPeakGram        *pgIn)
{ /*   *********************************************  </lalVerbatim> */
  
  INT4  length;
  INT4  nPeaks;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "LALUCHAR2HOUGHPeak", PEAKSELECTC);
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
/* *******************************  <lalVerbatim file="PeakSelectD"> */
void LALPeriodo2PSDrng (LALStatus  *status,
                     REAL8Periodogram1    *psd,
                     REAL8Periodogram1    *peri,
                     UINT2                *blocksRNG)
{ /*   *********************************************  </lalVerbatim> */

  INT4   length, j;
  UINT2  blockSize, nblocks2;
  REAL8  *data;
  REAL8  *medians;
  LALRunningMedianPar rngmedPar;
  REAL8Sequence mediansV;
  REAL8Sequence inputV;
  
  
  /* --------------------------------------------- */
  INITSTATUS (status, "LALPeriodo2PSDrng", PEAKSELECTC);
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
    data=peri->data;

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
/* *******************************  <lalVerbatim file="PeakSelectD"> */
void LALSelectPeakColorNoise(LALStatus  *status,
			     UCHARPeakGram        *pg,
			     REAL8                *thr,
			     REAL8PeriodoPSD      *in)
{ /*   *********************************************  </lalVerbatim> */ 
  
  INT4  length;
  INT4  nPeaks;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "LALSelectPeakColorNoise", PEAKSELECTC);
  ATTATCHSTATUSPTR (status);
  
  /*   Make sure the arguments are not NULL: */
  ASSERT (pg,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (thr, status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  ASSERT (in,  status, PEAKSELECTH_ENULL, PEAKSELECTH_MSGENULL);
  
  pg->epoch.gpsSeconds     = in->psd.epoch.gpsSeconds;
  pg->epoch.gpsNanoSeconds = in->psd.epoch.gpsNanoSeconds;
  pg->timeBase     = in->psd.timeBase;
  pg->fminBinIndex = in->psd.fminBinIndex;
  
  length = in->psd.length;
  nPeaks = 0;

  ASSERT (length==pg->length,  status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  ASSERT (length==in->periodogram.length,  status, PEAKSELECTH_EVAL, PEAKSELECTH_MSGEVAL);
  
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

