/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpDTF1Chisq.c
 *
 * Author: Anderson, W. G., and Brown, D. A. and Messaritaki E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpDTF1ChisqCV">
Author: Anderson, W. G., and Brown D. A., DTF1-Modifications: Messaritaki E.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpDTF1Chisq.c}}
\label{ss:FindChirpDTF1Chisq.c}

Module to implement the $\chi^2$ veto for the DTF1 templates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpDTF1ChisqCP}
\idx{LALFindChirpDTF1ChisqVeto()}

\subsubsection*{Description}

The function \texttt{LALFindChirpDTF1ChisqVeto()} perfoms a $\chi^2$ veto on an
entire data segment using the corresponding algorithm for the DTF1 templates,
described below. On exit the vector \texttt{chisqVec} contains the value
$\chi^2(t_j)$ for the data segment.


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

\vfill{\footnotesize\input{FindChirpDTF1ChisqCV}}
</lalLaTeX>
#endif

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>
#include <lal/FindChirpDTF1.h>

NRCSID (FINDCHIRPDTF1CHISQC, "$Id$");


/* <lalVerbatim file="FindChirpDTF1ChisqCP"> */
void
LALFindChirpDTF1ChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqInput        *inputDTF1,
    FindChirpChisqParams       *params
    )
/* </lalVerbatim> */
{
  UINT4                 j, l;
  UINT4                 numPoints;

  REAL4                *chisq;

  COMPLEX8             *q;
  COMPLEX8             *qDTF1;
  COMPLEX8             *qtilde;
  COMPLEX8             *qtildeDTF1;

  UINT4                 numChisqPts;
  UINT4                 numChisqBins;
  UINT4                *chisqBin;
  UINT4                *chisqBinDTF1;
  REAL4                 chisqNorm;
  /* REAL4                 rhosq; */

  COMPLEX8             *qtildeBin;
  COMPLEX8             *qtildeBinDTF1;

  INITSTATUS( status, "LALFindChirpDTF1ChisqVeto", FINDCHIRPDTF1CHISQC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* check that the output pointer is non-null and has room to store data */
  ASSERT(chisqVec,status,FINDCHIRPCHISQH_ENULL,FINDCHIRPCHISQH_MSGENULL );
  ASSERT(chisqVec->data,status,FINDCHIRPCHISQH_ENULL,FINDCHIRPCHISQH_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the chisq bin vectors is reasonable */
  ASSERT( params->chisqBinVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVec->length > 0, status,
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );
  ASSERT( params->chisqBinVecDTF1, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVecDTF1->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVecDTF1->length > 0, status,
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );

  /* check that the fft plan exists */
  ASSERT( params->plan,status,FINDCHIRPCHISQH_ENULL,FINDCHIRPCHISQH_MSGENULL );

  /* check that the input exists */
  ASSERT( input, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputDTF1, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the input contains some data */
  ASSERT( input->qVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( input->qVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( input->qtildeVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( input->qtildeVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputDTF1->qVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputDTF1->qVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputDTF1->qtildeVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputDTF1->qtildeVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT( params->qtildeBinVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVec->length > 0, status,
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );
  ASSERT( params->qtildeBinVecDTF1, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVecDTF1->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVecDTF1->length > 0, status,
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );

  /* check that we are using the correct approximant */
  if ( params->approximant != DTF1 )
  {
    ABORT( status, FINDCHIRPCHISQH_EIAPX, FINDCHIRPCHISQH_MSGEIAPX );
  }


  /*
   *
   * point local pointers to structure pointers
   *
   */


  chisq        = chisqVec->data;

  numPoints    = input->qVec->length;
  q            = input->qVec->data;
  qDTF1         = inputDTF1->qVec->data;
  qtilde       = input->qtildeVec->data;
  qtildeDTF1    = inputDTF1->qtildeVec->data;

  numChisqPts     = params->chisqBinVec->length;
  numChisqBins    = numChisqPts - 1;
  chisqBin     = params->chisqBinVec->data;
  chisqBinDTF1  = params->chisqBinVecDTF1->data;
  chisqNorm    = params->norm ; /* norm squared */

  qtildeBinDTF1 = params->qtildeBinVecDTF1->data;
  qtildeBin    = params->qtildeBinVec->data;


  /* 
   *
   * fill the numBins time series vectors for the chi-squared statistic
   *
   */


  for ( l = 0; l < numChisqBins; ++l )
  {
    memset( qtildeBin, 0, numPoints * sizeof(COMPLEX8) );
    memset( qtildeBinDTF1, 0, numPoints * sizeof(COMPLEX8) );

    memcpy( qtildeBin + chisqBin[l], qtilde + chisqBin[l],
        (chisqBin[l+1] - chisqBin[l]) * sizeof(COMPLEX8) );
    memcpy( qtildeBinDTF1 + chisqBinDTF1[l], qtildeDTF1 + chisqBinDTF1[l],
        (chisqBinDTF1[l+1] - chisqBinDTF1[l]) * sizeof(COMPLEX8) );

    LALCOMPLEX8VectorFFT( status->statusPtr, params->qBinVecPtr[l],
        params->qtildeBinVec, params->plan );
    CHECKSTATUSPTR( status );
    LALCOMPLEX8VectorFFT( status->statusPtr, params->qBinVecPtrDTF1[l],
        params->qtildeBinVecDTF1, params->plan );
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
      REAL4 X1 = params->qBinVecPtr[l]->data[j].re;  
      REAL4 Y1 = params->qBinVecPtr[l]->data[j].im;

      REAL4 mod1 = ( ( X1 - q[j].re / (REAL4) (numChisqBins) ) * 
          ( X1 - q[j].re / (REAL4) (numChisqBins) ) + 
          ( Y1 - q[j].im / (REAL4) (numChisqBins) ) * 
          ( Y1 - q[j].im / (REAL4) (numChisqBins) ) ); 

      chisq[j] += chisqNorm * mod1 ;
    }
  }  

  for ( j = 0; j < numPoints; ++j )
  {
    for ( l = 0; l < numChisqBins; ++l )
    {

      REAL4 X2 = params->qBinVecPtrDTF1[l]->data[j].re;
      REAL4 Y2 = params->qBinVecPtrDTF1[l]->data[j].im;

      REAL4 mod2 = ( ( X2 - qDTF1[j].re / (REAL4) (numChisqBins) ) * 
          ( X2 - qDTF1[j].re / (REAL4) (numChisqBins) ) +
          ( Y2 - qDTF1[j].im / (REAL4) (numChisqBins) ) * 
          ( Y2 - qDTF1[j].im / (REAL4) (numChisqBins) ) );

      chisq[j] += chisqNorm * mod2;

    }
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
