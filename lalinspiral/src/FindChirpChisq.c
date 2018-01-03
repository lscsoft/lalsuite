/*
*  Copyright (C) 2010 Karsten Wiesner, Stas Babak, Duncan Brown, Eirini Messaritaki, Jolien Creighton, Reinhard Prix, Craig Robinson
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
 * File Name: FindChirpChisq.c
 *
 * Author: Wiesner, K., Anderson, W. G., and Brown, D. A.
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Anderson, W. G., and Brown D. A.
 * \file
 * \ingroup FindChirpChisq_h
 *
 * \brief Module to implement the \f$\chi^2\f$ veto for the stationary phase chirp.
 *
 * ### Description ###
 *
 * The function <tt>LALFindChirpChisqVeto()</tt> perfoms a \f$\chi^2\f$ veto on
 * an entire data segment using the algorithm described below. On exit the
 * vector \c chisqVec contains the value \f$\chi^2(t_j)\f$ for the data
 * segment.
 *
 * ### Algorithm ###
 *
 * chisq algorithm here
 *
 * ### Uses ###
 *
 * \code
 * LALCreateReverseComplexFFTPlan()
 * LALDestroyComplexFFTPlan()
 * LALCCreateVector()
 * LALCDestroyVector()
 * LALCOMPLEX8VectorFFT()
 * \endcode
 *
 * ### Notes ###
 *
 */

#include <stdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/ComplexFFT.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpChisq.h>
#include <lal/Chisq_GPU.h>
#include <lal/Chisq_CPU.h>
#include <lal/LALInspiralConfig.h>


void
LALFindChirpComputeChisqBins(
    LALStatus                  *status,
    UINT4Vector                *chisqBinVec,
    FindChirpSegment           *fcSeg,
    UINT4                       kmax
    )

{
  UINT4         k, incIdx;
  REAL4        *tmpltPower;
  UINT4        *chisqBin = NULL;
  UINT4         numChisqBins;
  UINT4         chisqPt;
  REAL4         increment;
  REAL4         nextBin;
  REAL4         partSum;

  INITSTATUS(status);
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
  for ( k = 1; k < incIdx; ++k )
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
  if ( k == fcSeg->tmpltPowerVec->length && chisqPt != numChisqBins )
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



void
LALFindChirpChisqVeto (
    LALStatus                  *status,
    REAL4Vector                *chisqVec,
    FindChirpChisqInput        *input,
    FindChirpChisqParams       *params
    )

{

  INITSTATUS(status);
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
    case PadeT1:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case FindChirpSP:
    case IMRPhenomB:
      break;
    default:
      ABORT( status, FINDCHIRPCHISQH_EIAPX, FINDCHIRPCHISQH_MSGEIAPX );
      break;
  }


#ifdef LALINSPIRAL_CUDA_ENABLED
  Chisq_GPU(chisqVec->data, input->qVec->data, input->qtildeVec->data, params->chisqBinVec->data,
      input->qVec->length, params->chisqBinVec->length - 1, sqrt( params->norm ));
#else
  Chisq_CPU(chisqVec->data, input->qVec->data, input->qtildeVec->data, params,
      input->qVec->length, params->chisqBinVec->length - 1, sqrt( params->norm ), status);
#endif

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
