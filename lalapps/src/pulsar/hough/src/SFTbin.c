/*-----------------------------------------------------------------------
 *
 * File Name: SFTbin.c
 *
 * Authors: Sintes, A.M.,  Krishnan, B. & inspired from Siemens, X.
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes May 21, 2003
 *            Modified by Krishnan on Feb 22, 2004
 *
 *-----------------------------------------------------------------------
 */

/** OBSOLETE -- Use LAL functions from SFTfileIO.c instead */

/************************************ <lalVerbatim file="SFTbinCV">
Author: Sintes, A.M., Krishnan, B.
$Id$
************************************* </lalVerbatim> */

/* <lalLaTeX>  *******************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\section{MOdule \texttt{SFTbin.c}}
\label{ss:SFTbin.c}
Routines for reading SFT binary files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{SFTbinD}
\index{\verb&ReadSFTbinHeader1()&}
\index{\verb&ReadCOMPLEX8SFTbinData1()&}
\index{\verb&ReadCOMPLEX16SFTbinData1()&}
\index{\verb&COMPLEX8SFT2Periodogram1()&}
\index{\verb&COMPLEX16SFT2Periodogram1()&}
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
\vfill{\footnotesize\input{SFTbinCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

*********************************************** </lalLaTeX> */


#include "./SFTbin.h"



NRCSID (SFTBINC, "$Id$");

/*
 * The functions that make up the guts of this module
 */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTbinD"> */
void COMPLEX8SFT2Periodogram1 (LALStatus  *status,
                   REAL8Periodogram1     *peri,
		   COMPLEX8SFTData1      *sft)
{ /*   *********************************************  </lalVerbatim> */

  INT4       length;
  
  /* --------------------------------------------- */
  INITSTATUS (status, "COMPLEX8SFT2Periodogram1", SFTBINC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (peri,  status, SFTBINH_ENULL, SFTBINH_MSGENULL);

  peri->epoch.gpsSeconds = sft->epoch.gpsSeconds;
  peri->epoch.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  peri->timeBase = sft->timeBase;
  peri->fminBinIndex = sft->fminBinIndex;
  
  length = sft->length;
  ASSERT (length==peri->length,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);

  if (length > 0){
    REAL8      *out;
    COMPLEX8   *in1;
    REAL8      re,im, factor;
    INT4      n;

    ASSERT (sft->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    ASSERT (peri->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    out = peri->data;
    in1 = sft->data;
    n= length;

    /* if data was properly normalized */
    factor = 1./peri->timeBase;
    
    /* if data is not normalized... this factor needs to be clarified  */
    /* note inconsistency with the previous function and quantity *factor*  there */
    
    /* if bin zero is included should be treated properly because of factor 2 */
    
    while (n-- >0){
      re = in1->re;
      im = in1->im;
      *out = (re*re + im*im) *factor; /* factor 2 still missing if one-sided*/
      ++out;
      ++in1;
    }
  }
  
  

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}
 



/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="SFTbinD"> */
void SFT2Periodogram (LALStatus  *status,
                   REAL8Periodogram1     *peri,
		   SFTtype      *sft)
{ /*   *********************************************  </lalVerbatim> */

  INT4       length;
  REAL8      f0, deltaF;
  /* --------------------------------------------- */
  INITSTATUS (status, "COMPLEX8SFT2Periodogram1", SFTBINC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (sft,   status, SFTBINH_ENULL, SFTBINH_MSGENULL);
  ASSERT (peri,  status, SFTBINH_ENULL, SFTBINH_MSGENULL);

  f0 = sft->f0;
  deltaF = sft->deltaF;  
  length = sft->data->length;

  peri->epoch.gpsSeconds = sft->epoch.gpsSeconds;
  peri->epoch.gpsNanoSeconds = sft->epoch.gpsNanoSeconds;
  peri->timeBase = 1.0/deltaF;
  peri->fminBinIndex = floor(f0/deltaF + 0.5);
  
  ASSERT (length==peri->length,  status, SFTBINH_EVAL, SFTBINH_MSGEVAL);

  if (length > 0){
    REAL8      *out;
    COMPLEX8   *in1;
    REAL8      re,im, factor;
    INT4      n;

    ASSERT (sft->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    ASSERT (peri->data, status, SFTBINH_ENULL,  SFTBINH_MSGENULL);
    out = peri->data;
    in1 = sft->data->data;
    n= length;

    /* if data was properly normalized */
    factor = 1./peri->timeBase;
    
    /* if data is not normalized... this factor needs to be clarified  */
    /* note inconsistency with the previous function and quantity *factor*  there */
    
    /* if bin zero is included should be treated properly because of factor 2 */
    
    while (n-- >0){
      re = in1->re;
      im = in1->im;
      *out = (re*re + im*im) *factor; /* factor 2 still missing if one-sided*/
      ++out;
      ++in1;
    }
  }
  
  

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}
 








