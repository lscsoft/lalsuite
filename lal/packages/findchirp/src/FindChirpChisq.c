/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpChisq.c
 *
 * Author: Anderson, W. G., and Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpChisqCV">
Author: Anderson, W. G., and Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpChisq.c}}
\label{ss:FindChirpChisq.c}

Module to implement the $\chi^2$ veto for the stationary phase chirp.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpChisqCP}
\idx{LALFindChirpChisqVeto()}

\subsubsection*{Description}

The function \texttt{LALFindChirpChisqVeto()} perfoms a $\chi^2$ veto on 
an entire data segment using the algorithm described below. On exit the
vector \texttt{chisqVec} contains the value $\chi^2(t_j)$ for the data
segment.

\subsubsection*{Algorithm}

chisq algorithm here

\subsubsection*{Uses}
\begin{verbatim}
LALCreateReverseComplexFFTPlan()
LALDestroyComplexFFTPlan()
LALCCreateVector()
LALCDestroyVector()
LALCOMPLEX8VectorFFT()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpChisqCV}}
</lalLaTeX>
#endif

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>

NRCSID (FINDCHIRPCHISQC, "$Id$");


/* <lalVerbatim file="FindChirpChisqCP"> */
void
LALFindChirpComputeChisqBins(
    LALStatus                  *status,
    UINT4Vector                *chisqBinVec,
    FindChirpSegment           *fcSeg,
    UINT4                       kmax
    )
/* </lalVerbatim> */
{
  UINT4         k, incIdx;
  REAL4        *tmpltPower;
  UINT4        *chisqBin = NULL;
  UINT4         numChisqBins;
  UINT4         chisqPt;
  REAL4         increment;
  REAL4         nextBin;
  REAL4         partSum;

  INITSTATUS( status, "LALFindChirpComputeChisqBins", FINDCHIRPCHISQC );
  ATTATCHSTATUSPTR( status );


  ASSERT( chisqBinVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( chisqBinVec->length > 1, status,
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );
  ASSERT( ! chisqBinVec->data, status,
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );
  
  ASSERT( fcSeg, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( fcSeg->segNorm, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( fcSeg->segNorm->data, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( fcSeg->segNorm->length, status, 
      FINDCHIRPCHISQH_ENUMZ, FINDCHIRPCHISQH_MSGENUMZ );
  ASSERT( fcSeg->tmpltPowerVec, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( fcSeg->tmpltPowerVec->length, status, 
      FINDCHIRPCHISQH_ENUMZ, FINDCHIRPCHISQH_MSGENUMZ );
  ASSERT( fcSeg->tmpltPowerVec->data, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );


  /*
   *
   * calculate the chisq bins for the segment and template
   *
   */


  /* number of chisq bins is one less than the number of bin boundaries */
  numChisqBins = chisqBinVec->length - 1;

  /* vector containing k^(-7/3) / S_h(f) */
  tmpltPower = fcSeg->tmpltPowerVec->data;

  /* compute the amount of power in a chisq bin */
  incIdx = kmax > fcSeg->segNorm->length-1 ? fcSeg->segNorm->length-1 : kmax;
  increment = fcSeg->segNorm->data[incIdx] / (REAL4) numChisqBins;

  /* allocate memory for the bin boundaries */
  chisqBin = chisqBinVec->data = (UINT4 *)
    LALCalloc( chisqBinVec->length, sizeof(UINT4) );
  if ( ! chisqBinVec->data )
  {
    ABORT( status, FINDCHIRPCHISQH_EALOC, FINDCHIRPCHISQH_MSGEALOC );
  }

  /* initalize the first bin boundary */
  nextBin   = increment;
  chisqPt   = 0;
  partSum   = 0.0;
  chisqBin[chisqPt++] = 0;

  /* calculate the frequencies of the chi-squared bin boundaries */
  for ( k = 1; k < fcSeg->tmpltPowerVec->length; ++k ) 
  {
    partSum += tmpltPower[k];
    if ( partSum >= nextBin ) 
    {
      chisqBin[chisqPt++] = k;
      nextBin += increment;
      if ( chisqPt == numChisqBins ) break;
    }
  }

  /* check that we have sucessfully allocated all the bins */
  if ( k = fcSeg->tmpltPowerVec->length && chisqPt != numChisqBins )
  {
    /* if we have reaced the end of the template power vec and not */
    /* allocated all the bin boundaries then there is a problem    */
    ABORT( status, FINDCHIRPCHISQH_EBINS, FINDCHIRPCHISQH_MSGEBINS );
  }

  /* the last bin boundary is at can be at Nyquist since   */
  /* qtilde is zero above the ISCO of the current template */
  chisqBin[numChisqBins] = fcSeg->data->data->length;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpChisqCP"> */
void
LALFindChirpChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqParams       *params
    )
/* </lalVerbatim> */
{
  UINT4                 j, l;
  UINT4                 numPoints;

  REAL4                *chisq;

  COMPLEX8             *q;
  COMPLEX8             *qtilde;

  UINT4                 numChisqPts;
  UINT4                 numChisqBins;
  UINT4                *chisqBin;
  REAL4                 chisqNorm;
  REAL4                 rhosq;

  COMPLEX8             *qtildeBin;

  INITSTATUS( status, "LALFindChirpChisqVeto", FINDCHIRPCHISQC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* check that the output pointer is non-null and has room to store data */ 
  ASSERT( chisqVec, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( chisqVec->data, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the chisq bin vector is reasonable */
  ASSERT( params->chisqBinVec, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVec->data, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVec->length > 0, status, 
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );

  /* check that the fft plan exists */
  ASSERT( params->plan, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the input exists */
  ASSERT( input, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the input contains some data */
  ASSERT( input->qVec, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( input->qVec->data, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( input->qtildeVec, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( input->qtildeVec->data, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT( params->qtildeBinVec, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVec->data, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVec->length > 0, status, 
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );

  /* check that we are using the correct approximant */
  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorF2:
    case GeneratePPN:
    case FindChirpSP:
      break;
    default:
      ABORT( status, FINDCHIRPCHISQH_EIAPX, FINDCHIRPCHISQH_MSGEIAPX );
      break;
  }


  /*
   *
   * point local pointers to structure pointers
   *
   */


  chisq     = chisqVec->data;

  numPoints = input->qVec->length;
  q         = input->qVec->data;
  qtilde    = input->qtildeVec->data;

  numChisqPts  = params->chisqBinVec->length;
  numChisqBins = numChisqPts - 1;
  chisqBin     = params->chisqBinVec->data;
  chisqNorm    = sqrt( params->norm );

  qtildeBin = params->qtildeBinVec->data;


  /* 
   *
   * fill the numBins time series vectors for the chi-squared statistic
   *
   */


  for ( l = 0; l < numChisqBins; ++l ) 
  {
    memset( qtildeBin, 0, numPoints * sizeof(COMPLEX8) );

    memcpy( qtildeBin + chisqBin[l], qtilde + chisqBin[l], 
        (chisqBin[l+1] - chisqBin[l]) * sizeof(COMPLEX8) );

    LALCOMPLEX8VectorFFT( status->statusPtr, params->qBinVecPtr[l], 
        params->qtildeBinVec, params->plan );
    CHECKSTATUSPTR( status );
  }


  /* 
   *
   * calculate the chi-squared value at each time
   *
   */


  memset( chisq, 0, numPoints * sizeof(REAL4) );

  for ( j = 0; j < numPoints; ++j ) 
  {
    for ( l = 0; l < numChisqBins; ++l ) 
    {
      REAL4 Xl = params->qBinVecPtr[l]->data[j].re;
      REAL4 Yl = params->qBinVecPtr[l]->data[j].im;
      REAL4 deltaXl = chisqNorm * Xl -
        (chisqNorm * q[j].re / (REAL4) (numChisqBins));
      REAL4 deltaYl = chisqNorm * Yl -
        (chisqNorm * q[j].im / (REAL4) (numChisqBins));

      chisq[j] += deltaXl * deltaXl + deltaYl * deltaYl;
    }
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
