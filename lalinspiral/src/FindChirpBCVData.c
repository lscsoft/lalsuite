/*
*  Copyright (C) 2007 Duncan Brown, Eirini Messaritaki, Jolien Creighton, Patrick Brady
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
 * File Name: FindChirpBCVData.c
 *
 * Author: Brown D. A. and Messaritaki E.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBCVDataCV">
Author: Brown, D. A. and Messaritaki E.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpBCVData.c}}
\label{ss:FindChirpBCVData.c}

\input{FindChirpBCVDataCDoc}

\vfill{\footnotesize\input{FindChirpBCVDataCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <lal/FindChirpBCV.h>


NRCSID (FINDCHIRPBCVDATAC, "$Id$");

/* <lalVerbatim file="FindChirpBCVDataCP"> */
void
LALFindChirpBCVData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    )
/* </lalVerbatim> */
{
  UINT4                 i, k;
  UINT4                 cut;
  CHAR                  infoMsg[512];

  REAL4                *w;
  REAL4                *amp;
  REAL4                *ampBCV;
  COMPLEX8             *wtilde;
  REAL4                *tmpltPower;
  REAL4                *tmpltPowerBCV;

  REAL4Vector          *dataVec;
  REAL4                *spec;
  COMPLEX8             *resp;

  COMPLEX8             *outputData;
  COMPLEX8             *outputDataBCV;

  UINT4                *chisqBin    = NULL;
  UINT4                *chisqBinBCV = NULL;
  UINT4                 numChisqBins;
#if 0
  UINT4                 chisqPt;
  REAL4                 increment;
  REAL4                 nextBin;
  REAL4                 partSum;
#endif
#if 0
  REAL4                 Power  = 0.0 ;
  REAL4                 PowerBCV = 0.0 ;
#endif
  REAL4                 I73 = 0.0;
  REAL4                 I53 = 0.0;
  REAL4                 I1 = 0.0;
  REAL4                 segNormSum;

  FindChirpSegment     *fcSeg;
  DataSegment          *dataSeg;


  INITSTATUS( status, "LALFindChirpBCVData", FINDCHIRPBCVDATAC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * make sure that the arguments are reasonable
   *
   */


  /* check that the output exists */
  ASSERT( fcSegVec, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": fcSegVec" );
  ASSERT( fcSegVec->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": fcSegVec->data" );
  ASSERT( fcSegVec->data->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": fcSegVec->data->data" );
  ASSERT( fcSegVec->data->data->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": fcSegVec->data->data->data" );
  ASSERT( fcSegVec->data->dataBCV, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": fcSegVec->data->dataBCV" );
  ASSERT( fcSegVec->data->dataBCV->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": fcSegVec->data->dataBCV->data" );

  /* check that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPBCVH_ENULL,
      FINDCHIRPBCVH_MSGENULL ": params" );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  ASSERT( params->approximant == BCV, status,
      FINDCHIRPBCVH_EMAPX, FINDCHIRPBCVH_MSGEMAPX );

  /* check that the workspace vectors exist */
  ASSERT( params->ampVec, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );
  ASSERT( params->ampVec->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );
  ASSERT( params->ampVecBCV, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );
  ASSERT( params->ampVecBCV->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );

  ASSERT( params->wVec, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );
  ASSERT( params->wVec->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );

  ASSERT( params->wtildeVec, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );
  ASSERT( params->wtildeVec->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );

  ASSERT( params->tmpltPowerVec, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );
  ASSERT( params->tmpltPowerVec->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );
  ASSERT( params->tmpltPowerVecBCV, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );
  ASSERT( params->tmpltPowerVecBCV->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );


  /* check that the fft plans exist */
  ASSERT( params->fwdPlan, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );
  ASSERT( params->invPlan, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL );

  /* check that the parameter values are reasonable */
  ASSERT( params->fLow >= 0, status,
      FINDCHIRPBCVH_EFLOW, FINDCHIRPBCVH_MSGEFLOW );
  ASSERT( params->dynRange > 0, status,
      FINDCHIRPBCVH_EDYNR, FINDCHIRPBCVH_MSGEDYNR );

  /* check that the input exists */
  ASSERT( dataSegVec, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": dataSegVec" );
  ASSERT( dataSegVec->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": dataSegVec->data" );
  ASSERT( dataSegVec->data->chan, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": dataSegVec->data->chan" );
  ASSERT( dataSegVec->data->chan->data, status,
      FINDCHIRPBCVH_ENULL, FINDCHIRPBCVH_MSGENULL
      ": dataSegVec->data->chan->data" );


  /*
   *
   * set up local segment independent pointers
   *
   */


  w             = params->wVec->data;
  amp           = params->ampVec->data;
  ampBCV        = params->ampVecBCV->data;
  wtilde        = params->wtildeVec->data;
  tmpltPower    = params->tmpltPowerVec->data;
  tmpltPowerBCV = params->tmpltPowerVecBCV->data;


  /*
   *
   * loop over data segments
   *
   */


  for ( i = 0; i < dataSegVec->length; ++i )
  {


    /*
     *
     * set up segment dependent pointers
     *
     */


    dataSeg      = &(dataSegVec->data[i]);
    fcSeg        = &(fcSegVec->data[i]);

    dataVec      = dataSeg->chan->data;
    spec         = dataSeg->spec->data->data;
    resp         = dataSeg->resp->data->data;

    outputData    = fcSeg->data->data->data;
    outputDataBCV = fcSeg->dataBCV->data->data;

    if ( fcSeg->chisqBinVec->length )
    {
      chisqBin     = fcSeg->chisqBinVec->data;
      chisqBinBCV  = fcSeg->chisqBinVecBCV->data;
      numChisqBins = fcSeg->chisqBinVec->length - 1;
    }
    else
    {
      numChisqBins = 0;
    }

    ASSERT( params->wtildeVec->length == fcSeg->data->data->length, status,
        FINDCHIRPBCVH_EMISM, FINDCHIRPBCVH_MSGEMISM );
    ASSERT( params->wtildeVec->length == fcSeg->dataBCV->data->length, status,
        FINDCHIRPBCVH_EMISM, FINDCHIRPBCVH_MSGEMISM );

    /* store the waveform approximant in the data segment */
    fcSeg->approximant = BCV;


    /*
     *
     * compute htilde and store in fcSeg
     *
     */


    LALForwardRealFFT( status->statusPtr, fcSeg->data->data,
        dataVec, params->fwdPlan );
    CHECKSTATUSPTR( status );
    LALForwardRealFFT( status->statusPtr, fcSeg->dataBCV->data,
        dataVec, params->fwdPlan );
    CHECKSTATUSPTR( status );

    /* compute strain */
    for ( k = 0; k < fcSeg->data->data->length; ++k )
    {
      REAL4 p = outputData[k].re;
      REAL4 q = outputData[k].im;
      REAL4 pBCV = outputDataBCV[k].re;
      REAL4 qBCV = outputDataBCV[k].im;
      REAL4 x = resp[k].re * params->dynRange;
      REAL4 y = resp[k].im * params->dynRange;

      outputData[k].re =  p*x - q*y;
      outputData[k].im =  p*y + q*x;
      outputDataBCV[k].re =  pBCV*x - qBCV*y;
      outputDataBCV[k].im =  pBCV*y + qBCV*x;
    }


    /*
     *
     * compute inverse power spectrum
     *
     */


    /* set low frequency cutoff inverse power spectrum */
    cut = params->fLow / dataSeg->spec->deltaF > 1 ?
      params->fLow / dataSeg->spec->deltaF : 1;
    snprintf( infoMsg, sizeof(infoMsg)/sizeof(*infoMsg),
        "low frequency cut off index = %d\n", cut );
    LALInfo( status, infoMsg );

    /* set inverse power spectrum to zero */
    memset( wtilde, 0, params->wtildeVec->length * sizeof(COMPLEX8) );

    /* compute inverse of S_v */
    for ( k = cut; k < params->wtildeVec->length; ++k )
    {
      if ( spec[k] == 0 )
      {
        ABORT( status, FINDCHIRPBCVH_EDIVZ, FINDCHIRPBCVH_MSGEDIVZ );
      }
      wtilde[k].re = 1.0 / spec[k];
    }

    /*
     *
     * truncate inverse power spectrum in time domain if required
     *
     */


    if ( params->invSpecTrunc )
    {
      /* compute square root of inverse power spectrum */
      for ( k = cut; k < params->wtildeVec->length; ++k )
      {
        wtilde[k].re = sqrt( wtilde[k].re );
      }

      /* set nyquist and dc to zero */
      wtilde[params->wtildeVec->length - 1].re = 0.0;
      wtilde[0].re                             = 0.0;

      /* transform to time domain */
      LALReverseRealFFT( status->statusPtr, params->wVec, params->wtildeVec,
          params->invPlan );
      CHECKSTATUSPTR (status);

      /* truncate in time domain */
      memset( w + params->invSpecTrunc/2, 0,
          (params->wVec->length - params->invSpecTrunc) * sizeof(REAL4) );

      /* transform to frequency domain */
      LALForwardRealFFT( status->statusPtr, params->wtildeVec, params->wVec,
          params->fwdPlan );
      CHECKSTATUSPTR (status);

      /* normalise fourier transform and square */
      {
        REAL4 norm = 1.0 / (REAL4) params->wVec->length;
        for ( k = cut; k < params->wtildeVec->length; ++k )
        {
          wtilde[k].re *= norm;
          wtilde[k].re *= wtilde[k].re;
          wtilde[k].im = 0.0;
        }
      }

      /* set nyquist and dc to zero */
      wtilde[params->wtildeVec->length - 1].re = 0.0;
      wtilde[0].re                             = 0.0;
    }

    /* set inverse power spectrum below cut to zero */
    memset( wtilde, 0, cut * sizeof(COMPLEX8) );

    /* convert from S_v to S_h */
    for ( k = cut; k < params->wtildeVec->length; ++k )
    {
      REAL4 respRe = resp[k].re * params->dynRange;
      REAL4 respIm = resp[k].im * params->dynRange;
      REAL4 modsqResp = (respRe * respRe + respIm * respIm);
      REAL4 invmodsqResp;
      if ( modsqResp == 0 )
      {
        ABORT( status, FINDCHIRPBCVH_EDIVZ, FINDCHIRPBCVH_MSGEDIVZ );
      }
      invmodsqResp = 1.0 / modsqResp;
      wtilde[k].re *= invmodsqResp;
    }


    /*
     *
     * compute BCV normalisation parameters a1, b1 and b2,
     * segment normalization, outputData, point fcSeg at data segment
     *
     */


#if 0
    Power    = 0.0;
    PowerBCV = 0.0;
#endif
    I73 = 0.0;
    I53 = 0.0;
    I1 = 0.0;


    for ( k = 0; k < cut; ++k )
    {
      outputData[k].re = 0.0;
      outputData[k].im = 0.0;
      outputDataBCV[k].re = 0.0;
      outputDataBCV[k].im = 0.0;
    }

    memset( tmpltPower, 0, params->tmpltPowerVec->length * sizeof(REAL4) );
    memset( tmpltPowerBCV,0, params->tmpltPowerVecBCV->length * sizeof(REAL4) );
    memset( fcSeg->a1->data, 0, fcSeg->a1->length * sizeof(REAL4) );
    memset( fcSeg->b1->data, 0, fcSeg->b1->length * sizeof(REAL4) );
    memset( fcSeg->b2->data, 0, fcSeg->b2->length * sizeof(REAL4) );
    memset( fcSeg->segNorm->data, 0, fcSeg->segNorm->length * sizeof(REAL4) );

    segNormSum = 0.0;
    /*
     * moments necessary for the calculation of
     * the BCV normalization parameters
     */

    for ( k = 1; k < fcSeg->data->data->length; ++k )
    {
      I73 += 4.0 * amp[k] * amp[k] * wtilde[k].re ;
      I53 += 4.0 * amp[k] *  ampBCV[k] * wtilde[k].re ;
      I1 += 4.0 * ampBCV[k] * ampBCV[k] * wtilde[k].re;

      segNormSum += amp[k] * amp[k] * wtilde[k].re;
      fcSeg->segNorm->data[k] = segNormSum;

      /* calculation of a1, b1 and b2 for each ending frequency */
      /* making sure to avoid division by 0                     */
      if ( (I1 - I53*I53/I73) > 0 )
      {
        fcSeg->b2->data[k] = 1.0 / sqrt( I1 - I53*I53/I73 ) ;
      }
      if ( I73 > 0 )
      {
        fcSeg->a1->data[k] = 1.0 / sqrt(I73) ;
        fcSeg->b1->data[k] = - I53 * fcSeg->b2->data[k] / I73 ;
      }
    }


    /* the following is only necessary if we are doing a chisq veto... */
    for ( k = 1; k < fcSeg->data->data->length; ++k )
    {
      tmpltPower[k]    = amp[k] * sqrt(wtilde[k].re);
      /* Power += tmpltPower[k]; */
      tmpltPowerBCV[k] = ampBCV[k] * sqrt(wtilde[k].re);
      /* PowerBCV += tmpltPowerBCV[k] ; */
    }

    for ( k = cut; k < fcSeg->data->data->length; ++k )
    {
      outputData[k].re  *= 4.0 * amp[k] * wtilde[k].re ;
      outputData[k].im  *= 4.0 * amp[k] * wtilde[k].re ;
      outputDataBCV[k].re *= 4.0 * ampBCV[k] * wtilde[k].re ;
      outputDataBCV[k].im *= 4.0 * ampBCV[k] * wtilde[k].re ;
    }

    /* set output frequency series parameters */
    strncpy( fcSeg->data->name, dataSeg->chan->name, LALNameLength );
    strncpy( fcSeg->dataBCV->name, dataSeg->chan->name, LALNameLength );

    fcSeg->data->epoch.gpsSeconds      = dataSeg->chan->epoch.gpsSeconds;
    fcSeg->data->epoch.gpsNanoSeconds  = dataSeg->chan->epoch.gpsNanoSeconds;
    fcSeg->dataBCV->epoch.gpsSeconds     = dataSeg->chan->epoch.gpsSeconds;
    fcSeg->dataBCV->epoch.gpsNanoSeconds = dataSeg->chan->epoch.gpsNanoSeconds;

    fcSeg->data->f0     = dataSeg->chan->f0;
    fcSeg->data->deltaF = 1.0 /
      ( (REAL8) dataSeg->chan->data->length * dataSeg->chan->deltaT ) ;
    fcSeg->dataBCV->f0     = dataSeg->chan->f0;
    fcSeg->dataBCV->deltaF = 1.0 /
      ( (REAL8) dataSeg->chan->data->length * dataSeg->chan->deltaT ) ;

    fcSeg->deltaT       = dataSeg->chan->deltaT;
    fcSeg->number       = dataSeg->number;
    fcSeg->analyzeSegment = dataSeg->analyzeSegment;

    /* store low frequency cutoff and invSpecTrunc in segment */
    fcSeg->fLow         = params->fLow;
    fcSeg->invSpecTrunc = params->invSpecTrunc;


#if 0
    /*
     *
     * calculate the chisq bins for the segment and template
     *
     */


    /* First set of chisq bins */
    if ( numChisqBins )
    {
      increment = Power / (REAL4) numChisqBins;
      nextBin   = increment;
      chisqPt   = 0;
      partSum   = 0.0;

      /* calculate the frequencies of the chi-squared bin boundaries */
      chisqBin[chisqPt++] = 0;

      for ( k = 1; k < fcSeg->data->data->length; ++k )
      {
        partSum += tmpltPower[k] ;
        if ( partSum >= nextBin )
        {
          chisqBin[chisqPt++] = k;
          nextBin += increment;
          if ( chisqPt == numChisqBins ) break;
        }
      }
      chisqBin[numChisqBins] = fcSeg->data->data->length;
    }

    /* Second set of chisq bins */
    if ( numChisqBins )
    {
      increment = PowerBCV / (REAL4) numChisqBins;
      nextBin   = increment;
      chisqPt   = 0;
      partSum   = 0.0;

      /* calculate the frequencies of the chi-squared bin boundaries */
      chisqBinBCV[chisqPt++] = 0;

      for ( k = 1; k < fcSeg->dataBCV->data->length; ++k )
      {
        partSum += tmpltPowerBCV[k] ;
        if ( partSum >= nextBin )
        {
          chisqBinBCV[chisqPt++] = k;
          nextBin += increment;
          if ( chisqPt == numChisqBins ) break;
        }
      }
      chisqBinBCV[numChisqBins] = fcSeg->dataBCV->data->length;
    }
#endif

  } /* end loop over data segments */


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
