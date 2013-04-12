/*
*  Copyright (C) 2007 Duncan Brown, Eirini Messaritaki, Jolien Creighton
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

/*-----------------------------------------------------------------------
 *
 * File Name: FindChirpBCVChisq.c
 *
 * Author: Anderson, W. G., and Brown, D. A. and Messaritaki E.
 *
 *-----------------------------------------------------------------------
 */

/**

\author Anderson, W. G., and Brown D. A., BCV-Modifications: Messaritaki E.
\file
\ingroup FindChirpBCV_h

\brief Module to implement the \f$\chi^2\f$ veto for the BCV templates.

\heading{Description}

The function <tt>LALFindChirpBCVChisqVeto()</tt> perfoms a \f$\chi^2\f$ veto on an
entire data segment using the corresponding algorithm for the BCV templates,
described below. On exit the vector \c chisqVec contains the value
\f$\chi^2(t_j)\f$ for the data segment.


\heading{Algorithm}

chisq algorithm here

\heading{Uses}
\code
LALCreateReverseComplexFFTPlan()
LALDestroyComplexFFTPlan()
LALCCreateVector()
LALCDestroyVector()
LALCOMPLEX8VectorFFT()
\endcode

*/

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>
#include <lal/FindChirpBCV.h>

void
LALFindChirpBCVChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqInput        *inputBCV,
    FindChirpChisqParams       *params
    )

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
  /* REAL4                 rhosq; */

  COMPLEX8             *qtildeBin;
  COMPLEX8             *qtildeBinBCV;

  INITSTATUS(status);
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

  /* check that we are using the correct approximant */
  if ( params->approximant != BCV )
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
  qBCV         = inputBCV->qVec->data;
  qtilde       = input->qtildeVec->data;
  qtildeBCV    = inputBCV->qtildeVec->data;

  numChisqPts     = params->chisqBinVec->length;
  numChisqBins    = numChisqPts - 1;
  chisqBin     = params->chisqBinVec->data;
  chisqBinBCV  = params->chisqBinVecBCV->data;
  chisqNorm    = params->norm ; /* norm squared */

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
      REAL4 X1 = crealf(params->qBinVecPtr[l]->data[j]);
      REAL4 Y1 = cimagf(params->qBinVecPtr[l]->data[j]);

      REAL4 mod1 = ( ( X1 - crealf(q[j]) / (REAL4) (numChisqBins) ) *
          ( X1 - crealf(q[j]) / (REAL4) (numChisqBins) ) +
          ( Y1 - cimagf(q[j]) / (REAL4) (numChisqBins) ) *
          ( Y1 - cimagf(q[j]) / (REAL4) (numChisqBins) ) );

      chisq[j] += chisqNorm * mod1 ;
    }
  }

  for ( j = 0; j < numPoints; ++j )
  {
    for ( l = 0; l < numChisqBins; ++l )
    {

      REAL4 X2 = crealf(params->qBinVecPtrBCV[l]->data[j]);
      REAL4 Y2 = cimagf(params->qBinVecPtrBCV[l]->data[j]);

      REAL4 mod2 = ( ( X2 - crealf(qBCV[j]) / (REAL4) (numChisqBins) ) *
          ( X2 - crealf(qBCV[j]) / (REAL4) (numChisqBins) ) +
          ( Y2 - cimagf(qBCV[j]) / (REAL4) (numChisqBins) ) *
          ( Y2 - cimagf(qBCV[j]) / (REAL4) (numChisqBins) ) );

      chisq[j] += chisqNorm * mod2;

    }
  }

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
