/*-----------------------------------------------------------------------
 *
 * File Name: Peak2PHMD.c
 *
 * Authors: Sintes, A.M., 
 *
 * Revision: $Id$
 *
 * History:   Created by Sintes June 20, 2001
 *            Modified...
 *
 *-----------------------------------------------------------------------
 */

/************************************ <lalVerbatim file="Peak2PHMDCV">
Author: Sintes, A. M. 
$Id$
************************************* </lalVerbatim> */


/* <lalLaTeX>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsection{Module \texttt{Peak2PHMD.c}}
\label{ss:Peak2PHMD.c}
Construction of Partial-Hough-Map-Derivatives ({\sc phmd}) given a peak-gram
and the look-up-table.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Prototypes}
\vspace{0.1in}
\input{Peak2PHMDD}
\index{\verb&LALHOUGHPeak2PHMD()&}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Description}

This routine produces a {\sc phmd} at a certain frequency for a given  peak-gram and
look-up-table.\\

The inputs are:

\verb@phmd->fBin@: The frequency bin of this {\sc phmd}.

\verb@*lut@: The look-up-table (of type  \verb@HOUGHptfLUT@)

\verb@*pg@: The peak-gram  (of type  \verb@HOUGHPeakGram@)\\


The function \verb@LALHOUGHPeak2PHMD@ makes sure that the  {\sc lut}, the
peak-gram and also the frequency of the {\sc phmd}
are compatible.\\

The output \verb@HOUGHphmd  *phmd@ is  a structure
containing the frequency bin of this {\sc phmd},
the total number of borders of each type ({\it Left and Right}) to be
marked, the pointers to the borders in the corresponding
look-up-table, plus  {\it border} effects of clipping  on a finite
patch. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Uses}
%\begin{verbatim}
%LALZDestroyVector()
%\end{verbatim}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\subsubsection*{Notes}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\vfill{\footnotesize\input{Peak2PHMDCV}}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
</lalLaTeX> */

#include <lal/PHMD.h>

NRCSID (PEAK2PHMDC, "$Id$");

/* *******************************  <lalVerbatim file="Peak2PHMDD"> */
void LALHOUGHPeak2PHMD (LALStatus    *status,
			HOUGHphmd    *phmd, /* partial Hough map derivative */
			HOUGHptfLUT  *lut, /* Look up table */
			HOUGHPeakGram *pg)  /* peakgram */
{ /* *********************************************  </lalVerbatim> */

  INT2    i,j;
  UCHAR   *column1P;
  INT4    fBinDif,shiftPeak,minPeakBin,maxPeakBin,thisPeak;
  INT2    relatIndex;
  UINT4   lengthLeft,lengthRight,n, searchIndex;
  UINT8   firstBin,lastBin,pgI,pgF;
  /* --------------------------------------------- */

  INITSTATUS (status, "LALHOUGHPeak2PHMD", PEAK2PHMDC);
  ATTATCHSTATUSPTR (status); 

  /*   Make sure the arguments are not NULL: */ 
  ASSERT (phmd, status, PHMDH_ENULL, PHMDH_MSGENULL);
  ASSERT (lut,  status, PHMDH_ENULL, PHMDH_MSGENULL);
  ASSERT (pg,   status, PHMDH_ENULL, PHMDH_MSGENULL);
  ASSERT (phmd->firstColumn, status, PHMDH_ENULL, PHMDH_MSGENULL);
  
  /*  Make sure there are elements in firstColumn */
  ASSERT (phmd->ySide, status, PHMDH_ESIZE, PHMDH_MSGESIZE);


  /* Make sure peakgram and lut have same frequency discretization */
  ASSERT (lut->deltaF == pg->deltaF,  status, PHMDH_EVAL, PHMDH_MSGEVAL);

  /* Make sure phmd.fBin and lut are compatible */
  fBinDif = abs( (phmd->fBin) - (lut->f0Bin) );
  ASSERT (fBinDif < lut->nFreqValid, status, PHMDH_EFREQ, PHMDH_MSGEFREQ);

  pgI = pg->fBinIni;
  pgF = pg->fBinFin;
  /* bounds of interval to look at in the peakgram */
  firstBin = (phmd->fBin)+(lut->iniBin);
  lastBin  = firstBin + (lut->nBin)-1;

 /* Make sure peakgram f-interval and phmd.fBin+lut are compatible */
  ASSERT ( pgI <= firstBin, status, PHMDH_EINT, PHMDH_MSGEINT);
  ASSERT ( pgF >= lastBin,  status, PHMDH_EINT, PHMDH_MSGEINT);

  /* -------------------------------------------   */
  /* initialization */
  n= pg->length;

  lengthLeft = 0;
  lengthRight= 0;

  column1P = &(phmd->firstColumn[0]);

  for(i=0; i< phmd->ySide; ++i ){
    *column1P = 0;
    ++column1P;
  }
  
  minPeakBin = firstBin - pgI;
  maxPeakBin =  lastBin - pgI;
    
  /* -------------------------------------------   */
  if(n){ /* only if there are peaks present */
    INT2  lb1,rb1,lb2,rb2;
    INT2  max1,min1,max2,min2;
    INT2  nBinPos;

    nBinPos = (lut->iniBin) + (lut->nBin) -1;

    /* searching for the initial peak to look at */
    searchIndex = (n * minPeakBin ) /(pgF-pgI+1);
    
    while( pg->peak[searchIndex -1] >=  minPeakBin )
      --searchIndex;

    while( pg->peak[searchIndex] <  minPeakBin )
      ++searchIndex;
    
    shiftPeak = pgI - (phmd->fBin);

    /* for all the interesting peaks (or none) */
    while( (thisPeak = pg->peak[searchIndex] ) <=  maxPeakBin ){
      /* relative Index */
      relatIndex = thisPeak + shiftPeak;

      i = relatIndex;
      if( relatIndex < 0 )  i =  nBinPos - relatIndex;

      /* Reading the bin information */
      lb1 = lut->bin[i].leftB1;
      rb1 = lut->bin[i].rightB1;
      lb2 = lut->bin[i].leftB2;
      rb2 = lut->bin[i].rightB2;
      
      max1 = lut->bin[i].piece1max;
      min1 = lut->bin[i].piece1min; 
      max2 = lut->bin[i].piece2max;
      min2 = lut->bin[i].piece2min; 

      /* border selection from lut */

      if(lb1){
	phmd->leftBorderP[lengthLeft] = &( lut->border[lb1] );
	++lengthLeft;
      }

      if(lb2){
	phmd->leftBorderP[lengthLeft] = &( lut->border[lb2] );
	++lengthLeft;
      }

      if(rb1){
	phmd->rightBorderP[lengthRight] = &( lut->border[rb1] );
	++lengthRight;
      }

      if(rb2){
	phmd->rightBorderP[lengthRight] = &( lut->border[rb2] );
	++lengthRight;
      }

      /* correcting 1st column */
      for(j=min1; j<=max1; ++j) { phmd->firstColumn[j] = 1; }
      for(j=min2; j<=max2; ++j) { phmd->firstColumn[j] = 1; }


      ++searchIndex;
    }

  }

  /* -------------------------------------------   */
 
  phmd->lengthLeft = lengthLeft;
  phmd->lengthRight= lengthRight;
  /* -------------------------------------------   */
  
  DETATCHSTATUSPTR (status);
  
  /* normal exit */
  RETURN (status);
}
