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

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>

NRCSID (FINDCHIRPCHISQC, "$Id$");

void
LALFindChirpChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqParams       *params
                      )
{
  UINT4                 j, l;
  UINT4                 numPoints;

  REAL4                *chisq;

  COMPLEX8             *q;
  COMPLEX8             *qtilde;

  UINT4                 numChisqPts;
  UINT4                 numChisqBins;
  UINT4                *chisqBin;

  COMPLEX8Vector       *qtildeBinVec;
  COMPLEX8             *qtildeBin;

  COMPLEX8Vector      **qBinVecPtr;

  INITSTATUS( status, "LALFindChirpChisqVeto", FINDCHIRPCHISQC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* check that the output pointer is non-null and has room to store data */ 
  ASSERT( chisqVec, status, FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );
  ASSERT( chisqVec->data, status, FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );

  /* check that the chisq bin vector is reasonable */
  ASSERT( params->chisqBinVec, status, 
      FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );
  ASSERT( params->chisqBinVec->data, status, 
      FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );
  ASSERT( params->chisqBinVec->length > 0, status, 
      FINDCHIRPCHISQ_ECHIZ, FINDCHIRPCHISQ_MSGECHIZ );

  /* check that the fft plan exists */
  ASSERT( params->plan, status, FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );

  /* check that the input exists */
  ASSERT( input, status, FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );

  /* check that the input contains some data */
  ASSERT( input->qVec, status, 
      FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );
  ASSERT( input->qVec->data, status, 
      FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );
  ASSERT( input->qtildeVec, status, 
      FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );
  ASSERT( input->qtildeVec->data, status, 
      FINDCHIRPCHISQ_ENULL, FINDCHIRPCHISQ_MSGENULL );


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


  /*
   *
   * create storage (this should be moved to an init)
   *
   */


  /* create one vector for the fourier domain data */
  qtildeBinVec = NULL;

  LALCCreateVector( status->statusPtr, &qtildeBinVec, numPoints );
  CHECKSTATUSPTR( status );

  qtildeBin = qtildeBinVec->data;

  /* create numBins vectors for the time domain data */
  qBinVecPtr = (COMPLEX8Vector **) 
    LALMalloc( numChisqBins * sizeof(COMPLEX8Vector*) );

  for ( l = 0; l < numChisqBins; ++l )
  {
    qBinVecPtr[l] = NULL;

    LALCCreateVector( status->statusPtr, qBinVecPtr + l, numPoints );
    CHECKSTATUSPTR( status );
  }


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

    LALCOMPLEX8VectorFFT( status->statusPtr, qBinVecPtr[l], qtildeBinVec,
        params->plan );
    CHECKSTATUSPTR( status );
  }


  /* 
   *
   * calculate the chi-squared value at each time
   *
   */


#if 0
  /* this is for debugging */
  {
    memset( chisq, 0, numPoints * sizeof(REAL4) );

    for ( j = 0; j < numPoints; ++j ) 
    {
      REAL4 sumXl = 0.0;
      REAL4 sumYl = 0.0;

      for ( l = 0; l < numChisqBins; ++l ) 
      {
        sumXl += params->chisqNorm * qBinVecPtr[l]->data[j].re;
        sumYl += params->chisqNorm * qBinVecPtr[l]->data[j].im;
      }
      chisq[j] = sumXl * sumXl + sumYl * sumYl;
    }
  }
#endif

  memset( chisq, 0, numPoints * sizeof(REAL4) );

  for ( j = 0; j < numPoints; ++j ) 
  {
    for ( l = 0; l < numChisqBins; ++l ) 
    {
      REAL4 Xl = qBinVecPtr[l]->data[j].re;
      REAL4 Yl = qBinVecPtr[l]->data[j].im;
      REAL4 deltaXl = params->chisqNorm * Xl -
        (params->chisqNorm * q[j].re / (REAL4) (numChisqBins));
      REAL4 deltaYl = params->chisqNorm * Yl -
        (params->chisqNorm * q[j].im / (REAL4) (numChisqBins));

      chisq[j] += deltaXl * deltaXl + deltaYl * deltaYl;
    }
  }


  /*
   *
   * destroy storage (this should be moved to a finalize)
   *
   */


  for ( l = 0; l < numChisqBins; ++l )
  {
    LALCDestroyVector( status->statusPtr, (qBinVecPtr + l) );
    CHECKSTATUSPTR( status );
  }

  LALFree( qBinVecPtr );

  LALCDestroyVector( status->statusPtr, &qtildeBinVec );
  CHECKSTATUSPTR( status );


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
