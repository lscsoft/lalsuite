/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpChisq.c
 *
 * Author: Anderson, W. G., and Brown, D. A., BCV-Modifications: Messaritaki E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpChisqCV">
Author: Anderson, W. G., and Brown D. A., BCV-Modifications: Messaritaki E.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpChisq.c}}
\label{ss:FindChirpChisq.c}

Module to implement the $\chi^2$ veto.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpChisqCP}
\idx{LALFindChirpChisqVetoInit()}
\idx{LALFindChirpChisqVetoFinalize()}
\idx{LALFindChirpChisqVeto()}

\subsubsection*{Description}

The function \texttt{LALFindChirpChisqVetoInit()} takes as input the number of
bins required to contruct the $\chi^2$ veto and the number of points a data
segment as a parameter. The pointer \texttt{*params} must contain the
address of a structure of type \texttt{FindChirpChisqParams} for which storage
has already been allocated.  On exit this structure will be populated with the
correct values for execution of the function \texttt{LALFindChirpChisqVeto()}.
The workspace arrays and the inverse FFTW plan used by the veto will be
created.

The function \texttt{LALFindChirpChisqVetoFinalize()} takes the address of a
structure of type \texttt{FindChirpChisqParams} which has been populated by
\texttt{LALFindChirpChisqVetoInit()} as input. It takes the number of bins
required to contruct the $\chi^2$ veto and as a parameter. On exit all memory
allocated by the \texttt{LALFindChirpChisqVetoInit()} will be freed.

The function \texttt{LALFindChirpChisqVeto()} perfoms a $\chi^2$ veto on 
an entire data segment using the algorithm described below. On exit the
vector \texttt{chisqVec} contains the value $\chi^2(t_j)$ for the data
segment.

The function \texttt{LALFindChirpBCVChisqVeto()} perfoms a $\chi^2$ veto on 
an entire data segment using the corresponding algorithm for the BCV templates,
described below. On exit the
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
LALFindChirpChisqVetoInit (
    LALStatus                  *status,
    FindChirpChisqParams       *params,
    UINT4                       numChisqBins,
    UINT4                       numPoints
    )
/* </lalVerbatim> */
{
  UINT4                         l, m;

  INITSTATUS( status, "FindChirpChisqVetoInit", FINDCHIRPCHISQC );
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
  ASSERT( ! params->qtildeBinVecBCV, status,
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );
  ASSERT( ! params->qtildeBinVecBCV, status,
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


  /* create one vector for the additional BCV fourier domain data */
  LALCCreateVector( status->statusPtr,
      &(params->qtildeBinVecBCV), numPoints );
  BEGINFAIL( status )
  {
    TRY( LALDestroyComplexFFTPlan( status->statusPtr,
          &(params->plan) ), status );
    TRY( LALCDestroyVector( status->statusPtr,
	  &(params->qtildeBinVec) ), status );
  }
  ENDFAIL( status );


  /* create numBins vectors for the time domain data */
  params->qBinVecPtr = (COMPLEX8Vector **) 
    LALCalloc( 1, numChisqBins * sizeof(COMPLEX8Vector*) );
  if ( ! params->qBinVecPtr )
  {
    TRY( LALCDestroyVector( status->statusPtr, 
           &(params->qtildeBinVec) ), status );
    TRY( LALCDestroyVector( status->statusPtr,
           &(params->qtildeBinVecBCV) ), status );
    TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
          &(params->plan) ), status );
    ABORT( status, FINDCHIRPCHISQH_EALOC, FINDCHIRPCHISQH_MSGEALOC );
  }

  for ( l = 0; l < numChisqBins; ++l )
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
      TRY( LALCDestroyVector( status->statusPtr,
            &(params->qtildeBinVecBCV) ), status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
            &(params->plan) ), status );
    }
    ENDFAIL( status );
  }


  /* create additional numBins vectors for the BCV time domain data */
  params->qBinVecPtrBCV = (COMPLEX8Vector **)
    LALCalloc( 1, numChisqBins * sizeof(COMPLEX8Vector*) );
  if ( ! params->qBinVecPtrBCV )
  {
    TRY( LALCDestroyVector( status->statusPtr,
           &(params->qtildeBinVec) ), status );	    
    TRY( LALCDestroyVector( status->statusPtr,
           &(params->qtildeBinVecBCV) ), status );
    TRY( LALDestroyComplexFFTPlan( status->statusPtr,
           &(params->plan) ), status );
    ABORT( status, FINDCHIRPCHISQH_EALOC, FINDCHIRPCHISQH_MSGEALOC );
  }

  for ( l = 0; l < numChisqBins; ++l )
  {
    LALCCreateVector( status->statusPtr, params->qBinVecPtrBCV + l, numPoints );
    BEGINFAIL( status )
    {
      for ( m = 0; m < l ; ++m )
      {
        TRY( LALCDestroyVector( status->statusPtr,
              params->qBinVecPtrBCV + m ), status );
      }
      for ( m = 0; m < l ; ++m )
      {
        TRY( LALCDestroyVector( status->statusPtr,
              params->qBinVecPtr + m ), status );
      }
      LALFree( params->qBinVecPtrBCV );
      TRY( LALCDestroyVector( status->statusPtr,
            &(params->qtildeBinVec) ), status );
      TRY( LALCDestroyVector( status->statusPtr,
            &(params->qtildeBinVecBCV) ), status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(params->plan) ), status );
    }
    ENDFAIL( status );
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



/* <lalVerbatim file="FindChirpChisqCP"> */
void
LALFindChirpChisqVetoFinalize (
    LALStatus                  *status,
    FindChirpChisqParams       *params,
    UINT4                       numChisqBins
    )
/* </lalVerbatim> */
{
  UINT4                         l;

  INITSTATUS( status, "FindChirpChisqVetoInit", FINDCHIRPCHISQC );
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
  ASSERT( params->qtildeBinVecBCV, status,
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );
  ASSERT( params->qtildeBinVecBCV, status,
      FINDCHIRPCHISQH_ENNUL, FINDCHIRPCHISQH_MSGENNUL );


  /*
   *
   * destroy storage
   *
   */


  for ( l = 0; l < numChisqBins; ++l )
  {
    LALCDestroyVector( status->statusPtr, (params->qBinVecPtr + l) );
    CHECKSTATUSPTR( status );
  }

  for ( l = 0; l < numChisqBins; ++l )
  {
    LALCDestroyVector( status->statusPtr, (params->qBinVecPtrBCV + l) );
    CHECKSTATUSPTR( status );
  }

  LALFree( params->qBinVecPtr );
  LALFree( params->qBinVecPtrBCV );

  LALCDestroyVector( status->statusPtr, &(params->qtildeBinVec) );
  CHECKSTATUSPTR( status );
  LALCDestroyVector( status->statusPtr, &(params->qtildeBinVecBCV) );
  CHECKSTATUSPTR( status );


  /* destroy plan for chisq filter */
  LALDestroyComplexFFTPlan( status->statusPtr, &(params->plan) );
  CHECKSTATUSPTR( status );


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
#if 0
  REAL4                 mismatch;
#endif

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



/* <lalVerbatim file="FindChirpChisqCP"> */
void
LALFindChirpBCVChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqInput        *inputBCV,
    FindChirpChisqParams       *params
    )
/* </lalVerbatim> */
{
  UINT4                 j, l;
  UINT4                 numPoints;

  REAL4                *chisq;

  COMPLEX8             *q;
  COMPLEX8             *qBCV;
  COMPLEX8             *qtilde;
  COMPLEX8             *qtildeBCV;

  UINT4                 numChisqPts;
  UINT4                 numChisqBins;
  UINT4                *chisqBin;
  UINT4                *chisqBinBCV;
  REAL4                 chisqNorm;
  REAL4                 rhosq;
#if 0
  REAL4                 mismatch;
#endif

  COMPLEX8             *qtildeBin;
  COMPLEX8             *qtildeBinBCV;

  INITSTATUS( status, "LALFindChirpBCVChisqVeto", FINDCHIRPCHISQC );
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
  ASSERT( params->chisqBinVecBCV, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVecBCV->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->chisqBinVecBCV->length > 0, status,
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );


  /* check that the fft plan exists */
  ASSERT( params->plan,status,FINDCHIRPCHISQH_ENULL,FINDCHIRPCHISQH_MSGENULL );

  /* check that the input exists */
  ASSERT( input, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputBCV, status, FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );

  /* check that the input contains some data */
  ASSERT( input->qVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( input->qVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( input->qtildeVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( input->qtildeVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputBCV->qVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputBCV->qVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputBCV->qtildeVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( inputBCV->qtildeVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );


  /* check that the workspace vectors exist */
  ASSERT( params->qtildeBinVec, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVec->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVec->length > 0, status,
      FINDCHIRPCHISQH_ECHIZ, FINDCHIRPCHISQH_MSGECHIZ );
  ASSERT( params->qtildeBinVecBCV, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVecBCV->data, status,
      FINDCHIRPCHISQH_ENULL, FINDCHIRPCHISQH_MSGENULL );
  ASSERT( params->qtildeBinVecBCV->length > 0, status,
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


  chisq        = chisqVec->data;

  numPoints    = input->qVec->length;
  q            = input->qVec->data;
  qBCV         = inputBCV->qVec->data;
  qtilde       = input->qtildeVec->data;
  qtildeBCV    = inputBCV->qtildeVec->data;

  numChisqPts  = params->chisqBinVec->length;
  numChisqBins = numChisqPts - 1;
  chisqBin     = params->chisqBinVec->data;
  chisqBinBCV  = params->chisqBinVecBCV->data;
  chisqNorm    = sqrt( params->norm );
#if 0
  mismatch     = 1.0 - params->bankMatch;
#endif


    qtildeBinBCV = params->qtildeBinVecBCV->data;
    qtildeBin    = params->qtildeBinVec->data;


  /* 
   *
   * fill the numBins time series vectors for the chi-squared statistic
   *
   */


  for ( l = 0; l < numChisqBins; ++l )
  {
    memset( qtildeBin, 0, numPoints * sizeof(COMPLEX8) );
    memset( qtildeBinBCV, 0, numPoints * sizeof(COMPLEX8) );

    memcpy( qtildeBin + chisqBin[l], qtilde + chisqBin[l],
        (chisqBin[l+1] - chisqBin[l]) * sizeof(COMPLEX8) );
    memcpy( qtildeBinBCV + chisqBinBCV[l], qtildeBCV + chisqBinBCV[l],
        (chisqBinBCV[l+1] - chisqBinBCV[l]) * sizeof(COMPLEX8) );

    LALCOMPLEX8VectorFFT( status->statusPtr, params->qBinVecPtr[l],
        params->qtildeBinVec, params->plan );
    CHECKSTATUSPTR( status );
    LALCOMPLEX8VectorFFT( status->statusPtr, params->qBinVecPtrBCV[l],
	params->qtildeBinVecBCV, params->plan );
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
      REAL4 X2 = params->qBinVecPtrBCV[l]->data[j].re; 
      REAL4 Y2 = params->qBinVecPtrBCV[l]->data[j].im;

#if 0
      REAL4 deltaXl = chisqNorm * Xl -
        (chisqNorm * q[j].re / (REAL4) (numChisqBins));
      REAL4 deltaYl = chisqNorm * Yl -
        (chisqNorm * q[j].im / (REAL4) (numChisqBins));

      chisq[j] += deltaXl * deltaXl + deltaYl * deltaYl;
#endif

      REAL4 mod1 = ( ( X1 - q[j].re / (REAL4) (numChisqBins) ) * 
		     ( X1 - q[j].re / (REAL4) (numChisqBins) ) + 
	      	     ( Y1 - q[j].im / (REAL4) (numChisqBins) ) * 
		     ( Y1 - q[j].im / (REAL4) (numChisqBins) ) ); 

      REAL4 mod2 = ( ( X2 - qBCV[j].re / (REAL4) (numChisqBins) ) * 
		     ( X2 - qBCV[j].re / (REAL4) (numChisqBins) ) +
		     ( Y2 - qBCV[j].im / (REAL4) (numChisqBins) ) * 
		     ( Y2 - qBCV[j].im / (REAL4) (numChisqBins) ) );

      chisq[j] += ( chisqNorm * chisqNorm * ( mod1 + mod2 ) ) ;	    

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

