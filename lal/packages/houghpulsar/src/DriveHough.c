/*-----------------------------------------------------------------------
 *
 * File Name: DriveHough.c
 *
 * Authors: Sintes, A.M., 
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes August 3, 2001
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */

/************************************ <lalVerbatim file="DriveHoughCV">
Author: Sintes, A. M. 
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Module \texttt{DriveHough.c}}
\label{ss:DriveHough.c}

Routines for building and updating the space of partial Hough map derivatives
({\sc phmd}), 
and related functions needed for the construction of  total Hough maps at
different frequencies and possible residual spin down parameters.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{DriveHoughD}
\index{\verb&LALHOUGHConstructSpacePHMD()&}
\index{\verb&LALHOUGHupdateSpacePHMDup()&}
\index{\verb&LALHOUGHupdateSpacePHMDdn()&}
\index{\verb&LALHOUGHConstructHMT()&}
\index{\verb&LALHOUGHComputeFBinMap()&}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

${}$

The function \verb&LALHOUGHConstructSpacePHMD()& constructs the space of {\sc
phmd} \verb&PHMDVectorSequence *phmdVS&, given a 
\verb&HOUGHPeakGramVector *pgV& and \verb&HOUGHptfLUTVector *lutV&.
The minimum frequency bin present corresponds to \verb&phmdVS->fBinMin& and the
total number of different frequencies is 
\verb&phmdVS->nfSize&. At this moment the \verb@fBinMin@ line corresponds to the
first row of the cylinder and \verb&phmdVS->breakLine& is set to zero. 
\verb&phmdVS->breakLine&
 $\in\, [0,$\verb&nfSize&) is {\it the pointer} which
 identifies the position of the  \verb@fBinMin@ row in the circular-cylinder
 buffer. \\
 
 
 The function \verb&LALHOUGHupdateSpacePHMDup()& updates the space of {\sc
phmd} increasing the frequency \verb&phmdVS->fBinMin& by one.\\

  The function \verb&LALHOUGHupdateSpacePHMDdn()& updates the space of {\sc
phmd} decreasing the frequency \verb&phmdVS->fBinMin& by one.\\

Given \verb@PHMDVectorSequence *phmdVS@, the space of {\sc phmd}, and 
\verb&UINT8FrequencyIndexVector *freqInd&, a structure containing the frequency
indices  of the   {\sc phmd} at different time stamps that have to be combined
to form a Hough map, the function \verb&LALHOUGHConstructHMT()& produces the
total Hough map.\\
  
  The function \verb&LALHOUGHComputeFBinMap()& computes the corresponding frequency bin of
  a  {\sc phmd} \verb&UINT8 *fBinMap& for a given  intrinsic
  frequency bin of a source \verb&UINT8 *f0Bin&, and information regarding the
  time and the 
  residual spin down parameters \verb&HOUGHResidualSpinPar *rs&.
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
\begin{verbatim}
LALHOUGHPeak2PHMD()
LALHOUGHAddPHMD2HD()
LALHOUGHIntegrHD2HT()
\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{DriveHoughCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
</lalLaTeX> */




#include <lal/LALHough.h>

NRCSID (DRIVEHOUGHC, "$Id$");

#define rint(x) floor((x)+0.5)

/*
 * The functions that make up the guts of this module
 */


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="DriveHoughD"> */
void LALHOUGHConstructSpacePHMD  (LALStatus            *status, 
				  PHMDVectorSequence   *phmdVS,
				  HOUGHPeakGramVector  *pgV, 
				  HOUGHptfLUTVector    *lutV) 
{ /*   *********************************************  </lalVerbatim> */

  UINT4    k,j;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  UINT8    fBinMin;   /* present minimum frequency bin */ 
  UINT8    fBin;      /* present frequency bin */
  REAL8    deltaF;    /* frequency resolution */


  /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHConstructSpacePHMD", DRIVEHOUGHC);
  ATTATCHSTATUSPTR (status); 


  /*   Make sure the arguments are not NULL: */ 
  ASSERT (phmdVS, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV,   status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV->pg,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV->lut,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (pgV->length == lutV->length, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (pgV->length == phmdVS->length, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);

  /* at the  beggining, the  fBinMin line corresponds to the first row */
  phmdVS->breakLine = 0; /* mark [0,nfSize) (of the circular buffer)
			    pointing to the starting of the fBinMin line */

  length = phmdVS->length;
  nfSize = phmdVS->nfSize; 
  fBinMin = phmdVS->fBinMin;
  deltaF =  phmdVS->deltaF = lutV->lut[0].deltaF;


  for ( k=0; k<length; ++k ){ 

    /* make sure all deltaF are consistent */
    ASSERT (deltaF == lutV->lut[k].deltaF, 
	    status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);

    fBin = fBinMin;

    for ( j=0; j<  nfSize; ++j ){ 
      phmdVS->phmd[ j*length+k ].fBin = fBin;

      TRY( LALHOUGHPeak2PHMD(status->statusPtr,
			     &(phmdVS->phmd[ j*length+k ]),
			     &(lutV->lut[k]), &(pgV->pg[k]) ), status);
      ++fBin;
    }
  }
 

  DETATCHSTATUSPTR (status);
   /* normal exit */
  RETURN (status);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="DriveHoughD"> */
void LALHOUGHupdateSpacePHMDup  (LALStatus            *status, 
				  PHMDVectorSequence   *phmdVS,
				  HOUGHPeakGramVector  *pgV, 
				  HOUGHptfLUTVector    *lutV) 
{ /*   *********************************************  </lalVerbatim> */
  UINT4    k,breakLine;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  UINT8    fBinMin;   /* minimum frequency bin */ 
  UINT8    fBin;      /* present frequency bin */
  REAL8    deltaF;    /* frequency resolution */


  /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHupdateSpacePHMDup", DRIVEHOUGHC);
  ATTATCHSTATUSPTR (status); 


  /*   Make sure the arguments are not NULL: */ 
  ASSERT (phmdVS, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV,   status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV->pg,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV->lut,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (pgV->length == lutV->length, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (pgV->length == phmdVS->length, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
 /* -------------------------------------------   */


  length = phmdVS->length;
  nfSize = phmdVS->nfSize; 
  deltaF =  phmdVS->deltaF;
  
  breakLine = phmdVS->breakLine; /* old Break Line */
  fBinMin = phmdVS->fBinMin; /* initial frequency value  */

  /* Make sure initial breakLine is in [0,nfSize)  */
  ASSERT ( breakLine < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
  ASSERT ( breakLine >= 0,     status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
  /* -------------------------------------------   */

  /* Updating the space of PHMD increasing frequency */ 
  
  fBin = fBinMin + nfSize;

  for ( k=0; k<length; ++k ){ 
    /* make sure all deltaF are consistent */
    ASSERT (deltaF == lutV->lut[k].deltaF, 
	    status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
   
    phmdVS->phmd[ breakLine*length+k ].fBin = fBin;
    TRY( LALHOUGHPeak2PHMD(status->statusPtr,
			   &(phmdVS->phmd[ breakLine*length+k ]),
			   &(lutV->lut[k]), &(pgV->pg[k]) ), status);
  }

  /* Shift fBinMin and its mark */
  ++phmdVS->fBinMin;
  
  phmdVS->breakLine = (breakLine +1) % nfSize; 
  /* mark [0,nfSize) (of the circular buffer, modulus nfSize)
     pointing to the starting of the new fBinMin line */
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="DriveHoughD"> */
void LALHOUGHupdateSpacePHMDdn  (LALStatus            *status, 
				 PHMDVectorSequence   *phmdVS,
				 HOUGHPeakGramVector  *pgV, 
				 HOUGHptfLUTVector    *lutV) 
{ /*   *********************************************  </lalVerbatim> */
  UINT4    k,breakLine;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  UINT8    fBinMin;   /* minimum frequency bin */ 
  UINT8    fBin;      /* present frequency bin */
  REAL8    deltaF;    /* frequency resolution */


  /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHupdateSpacePHMDdn", DRIVEHOUGHC);
  ATTATCHSTATUSPTR (status); 


  /*   Make sure the arguments are not NULL: */ 
  ASSERT (phmdVS, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV,   status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (pgV->pg,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (lutV->lut,    status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (pgV->length == lutV->length, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (pgV->length == phmdVS->length, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements to be computed*/
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  /* -------------------------------------------   */

  length = phmdVS->length;
  nfSize = phmdVS->nfSize; 
  deltaF =  phmdVS->deltaF;
  
  breakLine = phmdVS->breakLine; /* old Break Line */

  /* Make sure initial breakLine is in [0,nfSize)  */
  ASSERT ( breakLine < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
  ASSERT ( breakLine >= 0,     status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
  /* -------------------------------------------   */

  /* Updating the space of PHMD decreasing frequency */ 
  
  /* Shift fBinMin and its mark */
  fBin =  --phmdVS->fBinMin; /* initial frequency value  */
  
  phmdVS->breakLine = (breakLine + nfSize- 1) % nfSize; 
  /* mark [0,nfSize) (of the circular buffer, modulus nfSize)
     pointing to the starting of the new fBinMin line */
  
  breakLine = phmdVS->breakLine; /* the new Break Line */

  for ( k=0; k<length; ++k ){ 
    /* make sure all deltaF are consistent */
    ASSERT (deltaF == lutV->lut[k].deltaF, 
	    status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
   
    phmdVS->phmd[ breakLine*length+k ].fBin = fBin;
    TRY( LALHOUGHPeak2PHMD(status->statusPtr,
			   &(phmdVS->phmd[ breakLine*length+k ]),
			   &(lutV->lut[k]), &(pgV->pg[k]) ), status);
  }

 
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}

/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="DriveHoughD"> */
void LALHOUGHConstructHMT  (LALStatus                  *status, 
			    HOUGHMapTotal              *ht,
			    UINT8FrequencyIndexVector  *freqInd,
			    PHMDVectorSequence         *phmdVS)
{ /*   *********************************************  </lalVerbatim> */


  UINT4    k,j;
  UINT4    breakLine;
  UINT4    nfSize;    /* number of different frequencies */
  UINT4    length;    /* number of elements for each frequency */
  UINT8    fBinMin;   /* present minimum frequency bin */ 
  UINT8    fBin;      /* present frequency bin */
  UINT2    xSide,ySide;
 
  HOUGHMapDeriv hd; /* the Hough map derivative */

  /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHConstructHMT", DRIVEHOUGHC);
  ATTATCHSTATUSPTR (status); 

  /*   Make sure the arguments are not NULL: */ 
  ASSERT (phmdVS,  status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (ht,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (freqInd, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */

  ASSERT (phmdVS->phmd,  status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (freqInd->data, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
 /* -------------------------------------------   */

  /* Make sure there is no size mismatch */
  ASSERT (freqInd->length == phmdVS->length, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  ASSERT (freqInd->deltaF == phmdVS->deltaF, status, 
	  LALHOUGHH_ESZMM, LALHOUGHH_MSGESZMM);
  /* -------------------------------------------   */

  /* Make sure there are elements  */
  ASSERT (phmdVS->length, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (phmdVS->nfSize, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  /* -------------------------------------------   */
  
   /* Make sure the ht map contains some pixels */
  ASSERT (ht->xSide, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);
  ASSERT (ht->ySide, status, LALHOUGHH_ESIZE, LALHOUGHH_MSGESIZE);

  length = phmdVS->length;
  nfSize = phmdVS->nfSize; 
  
  fBinMin = phmdVS->fBinMin; /* initial frequency value  od the cilinder*/
  
  breakLine = phmdVS->breakLine;

  /* number of physical pixels */
  xSide = ht->xSide;
  ySide = ht->ySide;
  
  /* Make sure initial breakLine is in [0,nfSize)  */
  ASSERT ( breakLine < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
  ASSERT ( breakLine >= 0,     status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
  /* -------------------------------------------   */
  
  /* Initializing  hd map and memory allocation */
  hd.xSide = xSide;
  hd.ySide = ySide;
  hd.map = (HoughDT *)LALMalloc(ySide*(xSide+1)*sizeof(HoughDT));
  /* -------------------------------------------   */

 
  TRY( LALHOUGHInitializeHD(status->statusPtr, &hd), status);
  for ( k=0; k<length; ++k ){ 
    /* read the frequency index and make sure is in the proper interval*/
    fBin =freqInd->data[k] -fBinMin;

    ASSERT ( fBin < nfSize, status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
    ASSERT ( fBin >= 0,     status, LALHOUGHH_EVAL, LALHOUGHH_MSGEVAL);
 
    /* find index */
    j = (fBin + breakLine) % nfSize;

    /* Add the corresponding PHMD to HD */
    TRY( LALHOUGHAddPHMD2HD(status->statusPtr,
			    &hd, &(phmdVS->phmd[j*length+k]) ), status);
  }

  TRY( LALHOUGHIntegrHD2HT(status->statusPtr, ht, &hd), status);
  
  /* Free memory and exit */
  LALFree(hd.map);

  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}


/* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< */
/* *******************************  <lalVerbatim file="DriveHoughD"> */
void LALHOUGHComputeFBinMap (LALStatus             *status, 
			     UINT8                 *fBinMap, 
			     UINT8                 *f0Bin,
			     HOUGHResidualSpinPar  *rs)
{ /*   *********************************************  </lalVerbatim> */

  UINT4    i;
  INT4    shiftFBin;
  REAL8   shiftF;

  UINT4   spinOrder;
  REAL8   *spinF;
  REAL8   timeDiff;    /*  T(t)-T(t0) */
  REAL8   timeDiffProd;  
  REAL8   deltaF;  /*  df=1/TCOH  */
  /* --------------------------------------------- */
  INITSTATUS (status, "LALHOUGHComputeFBinMap", DRIVEHOUGHC);
  ATTATCHSTATUSPTR (status); 
  
  /*   Make sure the arguments are not NULL: */ 
  ASSERT (fBinMap, status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (f0Bin,   status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  ASSERT (rs,      status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
  /* -------------------------------------------   */
  
  /*   Make sure the Input/Output pointers are not the same */
  ASSERT (fBinMap != f0Bin, status, LALHOUGHH_ESAME, LALHOUGHH_MSGESAME);
  
  shiftFBin = 0;
  shiftF = 0.0;
  
  spinOrder = rs->spinRes.length;
  
  if(spinOrder){
    ASSERT (rs->spinRes.data , status, LALHOUGHH_ENULL, LALHOUGHH_MSGENULL);
    timeDiff = rs->timeDiff;
    timeDiffProd =  timeDiff;
    
    deltaF = rs->deltaF; 
    spinF = rs->spinRes.data;

    for (i=0; i<spinOrder; ++i ){
      shiftF += spinF[i] * timeDiffProd;
      timeDiffProd *= timeDiff;
    }
    shiftFBin = rint( shiftF/deltaF) ; /* positive or negative */
  }
  
  *fBinMap = *f0Bin + shiftFBin;
  
  DETATCHSTATUSPTR (status);
  /* normal exit */
  RETURN (status);
}


