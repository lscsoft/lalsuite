/*----------------------------------------------------------------------- 
 * 
 * File Name: TwoInterfFindChirpChisq.c
 *
 * Author: Noel, J. S., Anderson, W. G., and Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="TwoInterfFindChirpChisqCV">
Author: Noel, J. S., Anderson, W. G., and Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{TwoInterfFindChirpChisq.c}}
\label{ss:TwoInterfFindChirpChisq.c}

Module to implement the $\chi^2$ veto.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{TwoInterfFindChirpChisqCP}
\idx{LALTwoInterfFindChirpChisqVetoInit()}
\idx{LALTwoInterfFindChirpChisqVetoFinalize()}
\idx{LALTwoInterfFindChirpChisqVeto()}

\subsubsection*{Uses}
\begin{verbatim}
LALCreateReverseComplexFFTPlan()
LALDestroyComplexFFTPlan()
LALCCreateVector()
LALCDestroyVector()
LALCOMPLEX8VectorFFT()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{TwoInterfFindChirpChisqCV}}
</lalLaTeX>
#endif

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>

NRCSID (TWOINTERFFINDCHIRPCHISQC, "$Id$");

/* <lalVerbatim file="TwoInterfFindChirpChisqCP"> */
void
LALTwoInterfFindChirpChisqVetoInit (
    LALStatus                  *status,
    FindChirpChisqParams       *params,
    UINT4                       numChisqBins,
    UINT4                       numPoints
    )
/* </lalVerbatim> */
{
  UINT4                         l, m;

  INITSTATUS( status, "TwoInterfFindChirpChisqVetoInit", TWOINTERFFINDCHIRPCHISQC );
  ATTATCHSTATUSPTR( status );

  ASSERT( params, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  ASSERT( numPoints > 0, status,
      FINDCHIRPCHISQH_ENUMZ, FINDCHIRPCHISQH_MSGENUMZ );

  ASSERT( ! params->plan, status, 
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );
  ASSERT( ! params->qtildeBinVec, status, 
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );
  ASSERT( ! params->qtildeBinVec, status, 
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );


  /*
   *
   * if numChisqBins is zero, don't initialize anything
   *
   */


  if ( numChisqBins == 0 )
  {
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }



  /*
   *
   * create storage
   *
   */


  /* create plan for chisq filter */
  LALCreateReverseComplexFFTPlan( status->statusPtr, 
      &(params->plan), numPoints, 0 );
  CHECKSTATUSPTR( status );

  /* create one vector for the fourier domain data */
  LALCCreateVector( status->statusPtr, 
      &(params->qtildeBinVec), numPoints );
  BEGINFAIL( status )
  {
    TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
          &(params->plan) ), status );
  }
  ENDFAIL( status );

  /* create numBins vectors for the time domain data */
  params->qBinVecPtr = (COMPLEX8Vector **) 
    LALCalloc( 1, numChisqBins * sizeof(COMPLEX8Vector*) );
  if ( ! params->qBinVecPtr )
  {
    TRY( LALCDestroyVector( status->statusPtr, 
           &(params->qtildeBinVec) ), status );
    TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
          &(params->plan) ), status );
    ABORT( status, FINDCHIRPCHISQH_EALOC, FINDCHIRPCHISQH_MSGEALOC );
  }

  /* Create the chisqBin vector*/
  /* The same vector is reused when numChisqBin > 1 */

  l = 0;
  {
    LALCCreateVector( status->statusPtr, params->qBinVecPtr + l, numPoints );
    BEGINFAIL( status )
    {
      for ( m = 0; m < l ; ++m )
      {
        TRY( LALCDestroyVector( status->statusPtr, 
              params->qBinVecPtr + m ), status );
      }
      LALFree( params->qBinVecPtr );
      TRY( LALCDestroyVector( status->statusPtr, 
            &(params->qtildeBinVec) ), status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
            &(params->plan) ), status );
    }
    ENDFAIL( status );
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="TwoInterfFindChirpChisqCP"> */
void
LALTwoInterfFindChirpChisqVetoFinalize (
    LALStatus                  *status,
    FindChirpChisqParams       *params,
    UINT4                       numChisqBins
    )
/* </lalVerbatim> */
{
  UINT4                         l;

  INITSTATUS( status, "TwoInterfFindChirpChisqVetoInit", TWOINTERFFINDCHIRPCHISQC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * if numChisqBins is zero, don't finalize anything
   *
   */


  if ( numChisqBins == 0 )
  {
    DETATCHSTATUSPTR( status );
    RETURN( status );
  }


  /*
   *
   * check arguments
   *
   */


  ASSERT( params, status, 
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->plan, status, 
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );
  ASSERT( params->qtildeBinVec, status, 
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );
  ASSERT( params->qtildeBinVec, status, 
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );

  /*
   *
   * destroy storage
   *
   */
  
  l = 0;
  {
    LALCDestroyVector( status->statusPtr, (params->qBinVecPtr + l) );
    CHECKSTATUSPTR( status );
  }
  
  LALFree( params->qBinVecPtr );

  LALCDestroyVector( status->statusPtr, &(params->qtildeBinVec) );
  CHECKSTATUSPTR( status );

  /* destroy plan for chisq filter */
  LALDestroyComplexFFTPlan( status->statusPtr, &(params->plan) );
  CHECKSTATUSPTR( status );


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="TwoInterfFindChirpChisqCP"> */
void
LALTwoInterfFindChirpChisqVeto (
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

#if 0
  REAL4                 rhosq;
  REAL4                 mismatch;
#endif

  COMPLEX8             *qtildeBin;

  INITSTATUS( status, "LALTwoInterfFindChirpChisqVeto", TWOINTERFFINDCHIRPCHISQC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* check that the output pointer is non-null and has room to store data */ 
  ASSERT( chisqVec, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( chisqVec->data, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

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
  ASSERT( params->plan, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

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

#if 0
  /* check that the bank match has been set */
  if ( params->bankMatch < 0 || params->bankMatch >= 1 )
  {
    ABORT( status, FINDCHIRPCHISQH_EMTCH, FINDCHIRPCHISQH_MSGEMTCH );
  }
#endif

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
#if 0
  mismatch     = 1.0 - params->bankMatch;
#endif

  qtildeBin = params->qtildeBinVec->data;


  /* 
   *
   * fill the numBins time series (reused) vector for the chi-squared statistic
   *
   */
  
  memset( chisq, 0, numPoints * sizeof(REAL4) );
  
  for ( l = 0; l < numChisqBins; ++l ) {
    memset( qtildeBin, 0, numPoints * sizeof(COMPLEX8) );
    
    memcpy( qtildeBin + chisqBin[l], qtilde + chisqBin[l],  
	    (chisqBin[l+1] - chisqBin[l]) * sizeof(COMPLEX8) ); 
    
    LALCOMPLEX8VectorFFT( status->statusPtr, params->qBinVecPtr[0],  
			  params->qtildeBinVec, params->plan ); 
    CHECKSTATUSPTR( status );
    
    /* 
     *
     * calculate the chi-squared value at each time
     *
     */
    
    for ( j = 0; j < numPoints; ++j ) {
      REAL4 Xl = params->qBinVecPtr[0]->data[j].re; 
      REAL4 Yl = params->qBinVecPtr[0]->data[j].im;
      REAL4 deltaXl = chisqNorm * Xl - 
	(chisqNorm * q[j].re / (REAL4) (numChisqBins));
      REAL4 deltaYl = chisqNorm * Yl - 
	(chisqNorm * q[j].im / (REAL4) (numChisqBins)); 
      
      chisq[j] += deltaXl * deltaXl + deltaYl * deltaYl; 
    }
  }
  
#if 0
  /* now modify the value to compute the new veto */
  for ( j = 0; j < numPoints; ++j ) 
    {
      rhosq = params->norm * (q[j].re * q[j].re + q[j].im * q[j].im);
      chisq[j] /= 1.0 + rhosq * mismatch * mismatch / (REAL4) numChisqBins;
    }
#endif 
  
  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
