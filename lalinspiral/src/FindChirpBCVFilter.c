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
 * File Name: FindChirpBCVFilter.c
 *
 * Author: Brown, D. A. and Messaritaki E.
 *
 *-----------------------------------------------------------------------
 */

/**

\author Brown D. A. and Messaritaki E.
\file
\ingroup FindChirpBCV_h

\brief This module provides the core of the matched filter for binary inspiral
chirps for BCV templates.

*/

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <math.h>
#include <lal/LALErrno.h>
#include <lal/XLALError.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpBCV.h>

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

void
LALFindChirpBCVFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )

{
  UINT4                 j, k, kFinal/*, kOpt*/;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex = 0;
  UINT4                 ignoreIndex;
  REAL4                 UNUSED myfmin;
  REAL4                 deltaT, deltaF;
  REAL4                 norm;
  REAL4                 modqsqThresh;
  REAL4                 rhosqThresh;
  REAL4                 mismatch;
  REAL4                 chisqThreshFac = 0;
  /* REAL4                 modChisqThresh; */
  UINT4                 numChisqBins = 0;
  UINT4                 eventStartIdx = 0;
  UINT4                *chisqBin      = NULL;
  UINT4                *chisqBinBCV   = NULL;
  UINT4                 chisqPt;
  REAL4                 chirpTime     = 0;
  REAL4                 Power         = 0.0;
  REAL4                 PowerBCV      = 0.0;
  REAL4                 increment, nextBin, partSum;
  REAL4                *tmpltPower    = NULL;
  REAL4                *tmpltPowerBCV = NULL;
  BOOLEAN               haveChisq     = 0;
  COMPLEX8             *qtilde        = NULL;
  COMPLEX8             *qtildeBCV     = NULL;
  COMPLEX8             *q             = NULL;
  COMPLEX8             *qBCV          = NULL;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *inputDataBCV  = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;
  REAL4                 a1 = 0.0;
  REAL4                 b1 = 0.0;
  REAL4                 b2 = 0.0;
  REAL4                 UNUSED templateNorm;
  /* REAL4                 modqsq; */
  REAL4                 Num1, Num2, Den1, Den2;
  REAL4                 omega, InvTan1, InvTan2;
  REAL4                 m, fFinal;
  REAL4                 UNUSED psi0, UNUSED psi3;
  /* CHAR                  infomsg[256]; */


  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh >= 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );
  ASSERT( params->chisqThresh >= 0, status,
      FINDCHIRPH_ECHIT, FINDCHIRPH_MSGECHIT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT(params->qVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec->data,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);
  ASSERT(params->qVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVecBCV->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV->data,status, FINDCHIRPH_ENULL,
      FINDCHIRPH_MSGENULL);

  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInput,   status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInputBCV,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec )
  {
    ASSERT( params->rhosqVec->data->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* if we are doing a chisq, check we can store the data */
  if ( input->segment->chisqBinVec->length )
  {
    ASSERT( params->chisqVec, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->chisqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure the filter has been initialized for the correct approximant */
  if ( params->approximant != BCV )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }

  /* make sure that the template and the segment are both BCV */
  ASSERT( input->fcTmplt->tmplt.approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );


  /*
   *
   * point local pointers to input and output pointers
   *
   */


  /* workspace vectors */
  q    = params->qVec->data;
  qBCV = params->qVecBCV->data;
  qtilde    = params->qtildeVec->data;
  qtildeBCV = params->qtildeVecBCV->data;


  /* template and data */
  inputData     = input->segment->data->data->data;
  inputDataBCV  = input->segment->dataBCV->data->data;
  tmpltSignal   = input->fcTmplt->data->data;
  templateNorm  = input->fcTmplt->tmpltNorm;
  deltaT        = params->deltaT;


  /* the length of the chisq bin vec is the number of bin   */
  /* _boundaries_ so the number of chisq bins is length - 1 */

  if ( input->segment->chisqBinVec->length )
  {
    /*
     * at this point, numChisqBins is only used as a parameter
     * on the basis of which we decide whether we will do a chisq test or not.
     * the actual number of chisq bins is:
     */
    numChisqBins = input->segment->chisqBinVec->length - 1;
    chisqBin    = input->segment->chisqBinVec->data;
    chisqBinBCV = input->segment->chisqBinVecBCV->data;
    tmpltPower    = input->segment->tmpltPowerVec->data;
    tmpltPowerBCV = input->segment->tmpltPowerVecBCV->data;
  }


  /* number of points in a segment */
  if (params->qVec->length != params->qVecBCV->length)
  {
    ABORT(status, FINDCHIRPBCVH_EQLEN, FINDCHIRPBCVH_MSGEQLEN);
  }
  numPoints = params->qVec->length;

  /*
   * template parameters, since FindChirpBCVFilterSegment is run
   * for every template
   */
  psi0 = input->fcTmplt->tmplt.psi0;
  psi3 = input->fcTmplt->tmplt.psi3;
  fFinal = input->fcTmplt->tmplt.fFinal;

  /*
   *
   * compute viable search regions in the snrsq vector
   *
   */


  {
    /* length of the chirp:                                            */
    /* calculated according to chirplen, using the values of M and eta */
    /* for the BCV tempaltes, as calculated using psi0 and psi3        */

#if 0
    /* REAL4 eta = input->fcTmplt->tmplt.eta; */
    /* m1 and m2 are currently equal to 0 in the BCV template bank */
    REAL4 m1 = input->fcTmplt->tmplt.mass1;
    REAL4 m2 = input->fcTmplt->tmplt.mass2;
    REAL4 myfmin = input->segment->fLow;
    /* total mass and eta, for use in chirpTime calculation */
    REAL4 m =  fabs(psi3) / (16.0 * LAL_MTSUN_SI * LAL_PI * LAL_PI * psi0) ;
    REAL4 eta = 3.0 / (128.0*psi0 * pow( (m*LAL_MTSUN_SI*LAL_PI), (5.0/3.0)) );
    REAL4 c0 = 5*m*LAL_MTSUN_SI/(256*eta);
    REAL4 c2 = 743.0/252.0 + eta*11.0/3.0;
    REAL4 c3 = -32*LAL_PI/3;
    REAL4 c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
    REAL4 x  = pow(LAL_PI*m*LAL_MTSUN_SI*myfmin, 1.0/3.0);
    REAL4 x2 = x*x;
    REAL4 x3 = x*x2;
    REAL4 x4 = x2*x2;
    REAL4 x8 = x4*x4;
    chirpTime = fabs(c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8);
#endif

#if 0
    /* temporarily set chirpTime equal to 0.5 seconds */
    chirpTime = 0.5;
    deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );
#endif


    /* Calculate deltaEventIndex : the only acceptable clustering */
    /* is "window" method, for BCV                                */
    if ( params->clusterMethod == FindChirpClustering_window )
    {
      deltaEventIndex=(UINT4) rint((params->clusterWindow/params->deltaT)+1.0);
    }
    else if ( params->clusterMethod == FindChirpClustering_tmplt )
    {
      ABORT( status, FINDCHIRPBCVH_ECLUW, FINDCHIRPBCVH_MSGECLUW );
    }


    /* ignore corrupted data at start and end */
    ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

    if ( lalDebugLevel & LALINFO )
    {
      CHAR newinfomsg[256];

#if 0
      snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
          "m = %e eta = %e => %e seconds => %d points\n"
          "invSpecTrunc = %d => ignoreIndex = %d\n",
          m, eta, chirpTime, deltaEventIndex,
          input->segment->invSpecTrunc, ignoreIndex );
#endif

      snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
          "chirp time = %e seconds => %d points\n"
          "invSpecTrunc = %d => ignoreIndex = %d\n",
          chirpTime, deltaEventIndex,
          input->segment->invSpecTrunc, ignoreIndex );
      LALInfo( status, newinfomsg );
    }

    /* XXX check that we are not filtering corrupted data XXX */
    /* XXX this is hardwired to 1/4 segment length        XXX */
    if ( ignoreIndex > numPoints / 4 )
    {
      ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
    }
    /* XXX reset ignoreIndex to one quarter of a segment XXX */
    ignoreIndex = numPoints / 4;
  }

  if ( lalDebugLevel & LALINFO )
  {
    CHAR newinfomsg[256];

    snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, newinfomsg );
  }

  /* k that corresponds to fFinal */
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  kFinal = fFinal / deltaF < numPoints/2 ? floor(fFinal / deltaF)
    : floor(numPoints/2);
  myfmin = input->segment->fLow;

  /* assign the values to a1, b1 and b2 */
  a1 = input->segment->a1->data[kFinal];
  b1 = input->segment->b1->data[kFinal];
  b2 = input->segment->b2->data[kFinal];


  /*
   * if one or more of a1, b1 and b2 are 0, output a message
   */

  if ( !fabs(a1) || !fabs(b1) || !fabs(b2) )
  {
    if ( lalDebugLevel & LALINFO )
    {
       CHAR newinfomsg[256];
       snprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
              "a1 = %e b1 = %e b2 = %e\n"
              "fFinal = %e deltaF = %e numPoints = %d => kFinal = %d\n",
               a1, b1, b2, fFinal, deltaF, numPoints, kFinal );
       LALInfo( status, newinfomsg );
    }
  }


  /*
   *
   * compute qtilde, qtildeBCV, and q, qBCV
   * using the correct combination of inputData and inputDataBCV
   *
   */


  memset( qtilde,    0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCV, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r    = a1 * crealf(inputData[k]);
    REAL4 s    = a1 * cimagf(inputData[k]);
    REAL4 rBCV = b1 * crealf(inputData[k]) + b2 * crealf(inputDataBCV[k]);
    REAL4 sBCV = b1 * cimagf(inputData[k]) + b2 * cimagf(inputDataBCV[k]);
    REAL4 x = crealf(tmpltSignal[k]);
    REAL4 y = 0.0 - cimagf(tmpltSignal[k]); /* note complex conjugate */

    qtilde[k].realf_FIXME = r * x - s * y ;
    qtilde[k].imagf_FIXME = r * y + s * x ;
    qtildeBCV[k].realf_FIXME = rBCV * x - sBCV * y ;
    qtildeBCV[k].imagf_FIXME = rBCV * y + sBCV * x ;
  }

  /* inverse fft to get q, and qBCV */
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec,
      params->qtildeVec, params->invPlan );
  CHECKSTATUSPTR( status );
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVecBCV,
      params->qtildeVecBCV, params->invPlan );
  CHECKSTATUSPTR( status );



  /*
   *
   * calculate signal to noise squared
   *
   */



  /* if full snrsq vector is required, set it to zero */
  if ( params->rhosqVec )
    memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );

  rhosqThresh = params->rhosqThresh;

  norm = deltaT / ((REAL4) numPoints) ;
  /* notice difference from corresponding factor in the sp templates: */
  /* no factor of 4 (taken care of in inputData and inputDataBCV      */
  /* and no segnorm, since we already multiplied by a1, b1 and b2.    */


  /* normalized snr threhold */
  modqsqThresh = rhosqThresh / norm ;

  /* we threshold on the "modified" chisq threshold computed from       */
  /*   chisqThreshFac = delta^2 * norm / p                              */
  /*   rho^2 = norm * modqsq                                            */
  /*                                                                    */
  /* So we actually threshold on                                        */
  /*                                                                    */
  /*    r^2 < chisqThresh * ( 1 + modqsq * chisqThreshFac )             */
  /*                                                                    */
  /* which is the same as thresholding on                               */
  /*    r^2 < chisqThresh * ( 1 + rho^2 * delta^2 / p )                 */
  /* and since                                                          */
  /*    chisq = p r^2                                                   */
  /* this is equivalent to thresholding on                              */
  /*    chisq < chisqThresh * ( p + rho^2 delta^2 )                     */
  /*                                                                    */
  /* The raw chisq is stored in the database. this quantity is chisq    */
  /* distributed with 2p-2 degrees of freedom.                          */

  mismatch = 1.0 - input->fcTmplt->tmplt.minMatch;
  if ( !numChisqBins )
  {
    chisqThreshFac = norm * mismatch * mismatch / (REAL4) numChisqBins;
  }

  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec )
  {
    memcpy( params->rhosqVec->name, input->segment->data->name,
        LALNameLength * sizeof(CHAR) );
    memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch),
        sizeof(LIGOTimeGPS) );
    params->rhosqVec->deltaT = input->segment->deltaT;

    for ( j = 0; j < numPoints; ++j )
    {
      REAL4 modqsqSP  = crealf(q[j]) * crealf(q[j]) + cimagf(q[j]) * cimagf(q[j]) ;
      REAL4 modqsqBCV = crealf(qBCV[j]) * crealf(qBCV[j]) + cimagf(qBCV[j]) * cimagf(qBCV[j]) ;
      REAL4 ImProd = 2.0 * ( - crealf(q[j]) * cimagf(qBCV[j]) + crealf(qBCV[j]) * cimagf(q[j]) ) ;

      REAL4 newmodqsq = ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) *
        ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) ;

      params->rhosqVec->data->data[j] = norm * newmodqsq;
    }
  }


  /*
   * calculation of the chisq bin boundaries
   * only if we are going to do a chisq veto
   */

  numChisqBins = 0;
  /* the following happens only if we are doing a BCV chisq veto */
  if( input->segment->chisqBinVec->length )
  {
    /*
     * Decide whether we do a chisq veto and what the
     * correct number of chisq bins is, based on fFinal
     */
    if ( fFinal <= 200.0 )
    {
      numChisqBins = 0; /* no chisq veto for these frequencies */
    }
    else if ( fFinal > 200.0 && fFinal <= 450.0 )
    {
      numChisqBins = 4;
    }
    else if ( fFinal > 450.0 && fFinal <= 700.0 )
    {
      numChisqBins = 6;
    }
    else if ( fFinal > 700.0 )
    {
      numChisqBins = 8;
    }

    /* if we do a chisq veto, calculate the bin boundaries */
    if ( numChisqBins )
    {

      /* sum up the template power */
      for ( k = 1; k < kFinal; ++k )
      {
         Power    += 4.0 * a1 * a1 * tmpltPower[k] * tmpltPower[k];
         PowerBCV += 4.0 * ( b1 * tmpltPower[k] + b2 * tmpltPowerBCV[k] )
                         * ( b1 * tmpltPower[k] + b2 * tmpltPowerBCV[k] );
      }


      /* First set of chisq bins */
      increment = Power / (REAL4) numChisqBins ;
      nextBin   = increment;
      chisqPt   = 0;
      partSum   = 0.0;

      /* calculate the frequencies of the chi-squared bin boundaries */
      chisqBin[chisqPt++] = 0;

      for ( k = 1; k < kFinal; ++k )
      {
        partSum += 4.0 * a1 * a1 * tmpltPower[k] * tmpltPower[k];
        if ( partSum >= nextBin )
        {
          chisqBin[chisqPt++] = k;
          nextBin += increment;
          if ( chisqPt == numChisqBins ) break;
        }
      }
      chisqBin[numChisqBins] = input->segment->data->data->length;

      /* Second set of chisq bins */
      increment = PowerBCV / (REAL4) numChisqBins;
      nextBin   = increment;
      chisqPt   = 0;
      partSum   = 0.0;

      /* calculate the frequencies of the chi-squared bin boundaries */
      chisqBinBCV[chisqPt++] = 0;

      for ( k = 1; k < kFinal; ++k )
      {
        partSum += 4.0 * ( b1 * tmpltPower[k] + b2 * tmpltPowerBCV[k] )
                       * ( b1 * tmpltPower[k] + b2 * tmpltPowerBCV[k] );
        if ( partSum >= nextBin )
        {
          chisqBinBCV[chisqPt++] = k;
          nextBin += increment;
          if ( chisqPt == numChisqBins ) break;
        }
      }
      chisqBinBCV[numChisqBins] = input->segment->dataBCV->data->length;

    } /* end if ( numChisqBins ) */

  } /* end: if( input->segment->chisqBinVec->length ) */


  /* look for an event in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
    REAL4 modqsqSP  = crealf(q[j]) * crealf(q[j]) + cimagf(q[j]) * cimagf(q[j]) ;
    REAL4 modqsqBCV = crealf(qBCV[j]) * crealf(qBCV[j]) + cimagf(qBCV[j]) * cimagf(qBCV[j]) ;
    REAL4 ImProd = 2.0 * ( - crealf(q[j]) * cimagf(qBCV[j]) + crealf(qBCV[j]) * cimagf(q[j]) ) ;

    REAL4 newmodqsq = ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
        0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) *
      ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
        0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) ;


    /* if snrsq exceeds threshold at any point */
    if ( newmodqsq > modqsqThresh )
    {

      /* compute chisq vector if it does not exist and we want it */
      if ( ! haveChisq  && numChisqBins )
      {
        memset( params->chisqVec->data, 0,
            params->chisqVec->length * sizeof(REAL4) );

        /* pointers to chisq input */
        params->chisqInput->qtildeVec    = params->qtildeVec;
        params->chisqInput->qVec         = params->qVec;
        params->chisqInputBCV->qtildeVec = params->qtildeVecBCV;
        params->chisqInputBCV->qVec      = params->qVecBCV;

        /* pointers to the chisq bin vectors in the segment */
        params->chisqParams->chisqBinVec    = input->segment->chisqBinVec;
        params->chisqParams->chisqBinVecBCV = input->segment->chisqBinVecBCV;

        /* pointers to chisq normalization */
        params->chisqParams->norm           = norm;
        params->chisqParams->a1             = a1 ;
        params->chisqParams->b1             = b1 ;
        params->chisqParams->b2             = b2 ;
#if 0
        params->chisqParams->bankMatch   = input->fcTmplt->tmplt.minMatch;
#endif

        /* compute the chisq vector: this is slow! */
        LALFindChirpBCVChisqVeto( status->statusPtr, params->chisqVec,
            params->chisqInput, params->chisqInputBCV, params->chisqParams );
        CHECKSTATUSPTR (status);

        haveChisq = 1;
      }

      if ( ! numChisqBins )
      {
        memset( params->chisqVec->data, 0,
            params->chisqVec->length * sizeof(REAL4) );
      }


      /*
       * when we decide to impose a cut on alphaF, the calculation of
       * alphaF must be done at this point, and the check on alphaF
       * should be in the if statement that follows
       */
#if 0
      /* calculate alphaF */
      Num1 = qBCV[j].re + q[j].im ;
      Num2 = qBCV[j].re - q[j].im ;
      Den1 = q[j].re - qBCV[j].im ;
      Den2 = q[j].re + qBCV[j].im ;

      InvTan1 = (REAL4) atan2(Num1, Den1);
      InvTan2 = (REAL4) atan2(Num2, Den2);

      omega = 0.5 * InvTan1 + 0.5 * InvTan2 ;
      alpha = - b2 * tan(omega) / ( a1 + b1*tan(omega) );
      alpha *= pow(params->deltaT, 2.0/3.0);
      alphaF = alpha * pow(fFinal, 2.0/3.0);
      if ( (alphaF >= 0.0 && alphaF <= 2.0) &&
           ( ! numChisqBins || params->chisqVec->data[j] <
           (params->chisqThresh * ( 1.0 + newmodqsq * chisqThreshFac )) ) )
#endif

      /*
       * if we don't have a chisq or the chisq drops below the
       * modified chisq threshold, start processing events
       */
      if ( ! numChisqBins ||
          params->chisqVec->data[j] <
          (params->chisqThresh * ( 1.0 + newmodqsq * chisqThreshFac )) )
      {
        if ( ! *eventList )
        {
          /* store the start of the crossing */
          eventStartIdx = j;

          /* if this is the first event, start the list */
          thisEvent = *eventList = (SnglInspiralTable *)
            LALCalloc( 1, sizeof(SnglInspiralTable) );
          if ( ! thisEvent )
          {
            ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
          }

          /* record the data that we need for the clustering algorithm */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = newmodqsq;
        }
        else if ( ! params->clusterMethod == FindChirpClustering_none &&
            j <= thisEvent->end_time.gpsSeconds + deltaEventIndex &&
            newmodqsq > thisEvent->snr )
        {
          /* if this is the same event, update the maximum */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = newmodqsq;
        }
        else if (j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
              params->clusterMethod == FindChirpClustering_none )
        {
          /* clean up this event */
          SnglInspiralTable *lastEvent;
          INT8               timeNS;
          INT4               timeIndex = thisEvent->end_time.gpsSeconds;

          /* set the event LIGO GPS time of the event */
          timeNS = 1000000000L *
            (INT8) (input->segment->data->epoch.gpsSeconds);
          timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
          timeNS += (INT8) (1e9 * timeIndex * deltaT);
          thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
          thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
          thisEvent->end_time_gmst = fmod(XLALGreenwichMeanSiderealTime(
	      &thisEvent->end_time), LAL_TWOPI) * 24.0 / LAL_TWOPI;	/* hours */
          ASSERT( !XLAL_IS_REAL8_FAIL_NAN(thisEvent->end_time_gmst), status, LAL_FAIL_ERR, LAL_FAIL_MSG );

          /* set the impuse time for the event */
          thisEvent->template_duration = (REAL8) chirpTime;

          /* record the ifo and channel name for the event */
          strncpy( thisEvent->ifo, input->segment->data->name,
              2 * sizeof(CHAR) );
          strncpy( thisEvent->channel, input->segment->data->name + 3,
              (LALNameLength - 3) * sizeof(CHAR) );
          thisEvent->impulse_time = thisEvent->end_time;

          /* record coalescence phase and alpha */

          /* calculate the numerators and the denominators */
          Num1 = crealf(qBCV[timeIndex]) + cimagf(q[timeIndex]) ;
          Num2 = crealf(qBCV[timeIndex]) - cimagf(q[timeIndex]) ;
          Den1 = crealf(q[timeIndex]) - cimagf(qBCV[timeIndex]) ;
          Den2 = crealf(q[timeIndex]) + cimagf(qBCV[timeIndex]) ;

          InvTan1 = (REAL4) atan2(Num1, Den1);
          InvTan2 = (REAL4) atan2(Num2, Den2);

          thisEvent->coa_phase = 0.5 * InvTan1 - 0.5 * InvTan2 ;
          omega = 0.5 * InvTan1 + 0.5 * InvTan2 ;
          thisEvent->alpha = - b2 * tan(omega)
            / ( a1 + b1 * tan(omega) );
          thisEvent->alpha *= pow(params->deltaT, 2.0/3.0);

          /* copy the template into the event */
          thisEvent->psi0   = (REAL4) input->fcTmplt->tmplt.psi0;
          thisEvent->psi3   = (REAL4) input->fcTmplt->tmplt.psi3;
          /* chirp mass in units of M_sun */
          thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
            pow( 3.0 / 128.0 / input->fcTmplt->tmplt.psi0 , 3.0/5.0 );
          m =  fabs(thisEvent->psi3) /
            (16.0 * LAL_MTSUN_SI * LAL_PI * LAL_PI * thisEvent->psi0) ;
          thisEvent->eta = 3.0 / (128.0*thisEvent->psi0 *
              pow( (m*LAL_MTSUN_SI*LAL_PI), (5.0/3.0)) );
          thisEvent->f_final  = (REAL4) input->fcTmplt->tmplt.fFinal ;

          /* set the type of the template used in the analysis */
          snprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
              "FindChirpBCV" );

          /* set snrsq, chisq, sigma and effDist for this event */
          if ( input->segment->chisqBinVec->length )
          {
            /* we store chisq distributed with 2p - 2 degrees of freedom */
            /* in the database. params->chisqVec->data = r^2 = chisq / p */
            /* so we multiply r^2 by p here to get chisq                 */
            thisEvent->chisq =
              params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
            thisEvent->chisq_dof = 2 * numChisqBins; /* double for BCV */
          }
          else
          {
            thisEvent->chisq     = 0;
            thisEvent->chisq_dof = 0;
          }
          thisEvent->sigmasq = sqrt( norm / a1 );
          thisEvent->eff_distance =
            input->fcTmplt->tmpltNorm / norm / thisEvent->snr;
          thisEvent->eff_distance = sqrt( thisEvent->eff_distance ) /
            pow(params->deltaT, 1.0/6.0);

          thisEvent->snr *= norm;
          thisEvent->snr = sqrt( thisEvent->snr );

          /* compute the time since the snr crossing */
          thisEvent->event_duration =
            (REAL8) timeIndex - (REAL8) eventStartIdx;
          thisEvent->event_duration *= (REAL8) deltaT;

          /* store the start of the crossing */
          eventStartIdx = j;

          /* allocate memory for the newEvent */
          lastEvent = thisEvent;

          lastEvent->next = thisEvent = (SnglInspiralTable *)
            LALCalloc( 1, sizeof(SnglInspiralTable) );
          if ( ! lastEvent->next )
          {
            ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
          }

          /* stick minimal data into the event */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = newmodqsq;
        }
      }
    }
  }


  /*
   *
   * clean up the last event if there is one
   *
   */


  if ( thisEvent )
  {
    INT8           timeNS;
    INT4           timeIndex = thisEvent->end_time.gpsSeconds;

    /* set the event LIGO GPS time of the event */
    timeNS = 1000000000L *
      (INT8) (input->segment->data->epoch.gpsSeconds);
    timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
    timeNS += (INT8) (1e9 * timeIndex * deltaT);
    thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
    thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
    thisEvent->end_time_gmst = fmod(XLALGreenwichMeanSiderealTime(
        &thisEvent->end_time), LAL_TWOPI) * 24.0 / LAL_TWOPI;	/* hours */
    ASSERT( !XLAL_IS_REAL8_FAIL_NAN(thisEvent->end_time_gmst), status, LAL_FAIL_ERR, LAL_FAIL_MSG );

    /* set the impuse time for the event */
    thisEvent->template_duration = (REAL8) chirpTime;

    /* record the ifo name for the event */
    strncpy( thisEvent->ifo, input->segment->data->name,
        2 * sizeof(CHAR) );
    strncpy( thisEvent->channel, input->segment->data->name + 3,
        (LALNameLength - 3) * sizeof(CHAR) );
    thisEvent->impulse_time = thisEvent->end_time;

    /* record coalescence phase and alpha */

    /* calculate the numerators and the denominators */
    Num1 = crealf(qBCV[timeIndex]) + cimagf(q[timeIndex]) ;
    Num2 = crealf(qBCV[timeIndex]) - cimagf(q[timeIndex]) ;
    Den1 = crealf(q[timeIndex]) - cimagf(qBCV[timeIndex]) ;
    Den2 = crealf(q[timeIndex]) + cimagf(qBCV[timeIndex]) ;

    InvTan1 = (REAL4) atan2(Num1, Den1);
    InvTan2 = (REAL4) atan2(Num2, Den2 );


    thisEvent->coa_phase = 0.5 * InvTan1 - 0.5 * InvTan2 ;
    omega = 0.5 * InvTan1 + 0.5 * InvTan2 ;
    thisEvent->alpha = - b2 * tan(omega)
      / ( a1 + b1 * tan(omega) );
    thisEvent->alpha *= pow(params->deltaT, 2.0/3.0);


    /* copy the template into the event */
    thisEvent->psi0   = (REAL4) input->fcTmplt->tmplt.psi0;
    thisEvent->psi3   = (REAL4) input->fcTmplt->tmplt.psi3;
    /* chirp mass in units of M_sun */
    thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
      pow( 3.0 / 128.0 / input->fcTmplt->tmplt.psi0, 3.0/5.0 );
    thisEvent->f_final  = (REAL4) input->fcTmplt->tmplt.fFinal;
    m =  fabs(thisEvent->psi3) /
          (16.0 * LAL_MTSUN_SI * LAL_PI * LAL_PI * thisEvent->psi0) ;
    thisEvent->eta = 3.0 / (128.0*thisEvent->psi0 *
          pow( (m*LAL_MTSUN_SI*LAL_PI), (5.0/3.0)) );



    /* set the type of the template used in the analysis */
    snprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
        "FindChirpBCV" );

    /* set snrsq, chisq, sigma and effDist for this event */
    if ( input->segment->chisqBinVec->length )
    {
      /* we store chisq distributed with 2p - 2 degrees of freedom */
      /* in the database. params->chisqVec->data = r^2 = chisq / p */
      /* so we multiply r^2 by p here to get chisq                 */
      thisEvent->chisq =
        params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
      thisEvent->chisq_dof =  2 * numChisqBins; /* double for BCV */
    }
    else
    {
      thisEvent->chisq     = 0;
      thisEvent->chisq_dof = 0;
    }
    thisEvent->sigmasq = sqrt( norm / a1 );
    thisEvent->eff_distance = input->fcTmplt->tmpltNorm / norm /
      thisEvent->snr;
    thisEvent->eff_distance = sqrt( thisEvent->eff_distance ) /
      pow(params->deltaT,1.0/6.0);

    thisEvent->snr *=  norm ;
    thisEvent->snr = sqrt( thisEvent->snr );

    /* compute the time since the snr crossing */
    thisEvent->event_duration = (REAL8) timeIndex - (REAL8) eventStartIdx;
    thisEvent->event_duration *= (REAL8) deltaT;

  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
