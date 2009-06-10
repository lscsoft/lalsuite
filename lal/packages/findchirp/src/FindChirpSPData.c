/*
*  Copyright (C) 2007 Duncan Brown, Eirini Messaritaki, Gareth Jones, Jolien Creighton, Patrick Brady
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
 * File Name: FindChirpSPData.c
 *
 * Author: Brown D. A.
 *
 * Revision: $Id$
 *
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpSPDataCV">
Author: Brown, D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\subsection{Module \texttt{FindChirpSPData.c}}
\label{ss:FindChirpSPData.c}

\input{FindChirpSPDataCDoc}

\vfill{\footnotesize\input{FindChirpSPDataCV}}
</lalLaTeX>
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <math.h>

NRCSID (FINDCHIRPSPDATAC, "$Id$");


/* <lalVerbatim file="FindChirpSPDataCP"> */
void
LALFindChirpSPData (
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
  COMPLEX8             *wtilde;
  REAL4                *tmpltPower;

  REAL4Vector          *dataVec;
  REAL4                *spec;
  COMPLEX8             *resp;

  COMPLEX8             *outputData;
  REAL4                 segNormSum;

  /* stuff added for continous chisq test */
  REAL4Vector          *dataPower = NULL;
  REAL4		        PSDsum = 0;
  INT4 			startIX = 0;
  INT4			endIX = 0;
  INT4			sortFlag = 0;
  COMPLEX8Vector       *fftVec = NULL;
  FindChirpSegment     *fcSeg;
  DataSegment          *dataSeg;
  INITSTATUS( status, "LALFindChirpSPData", FINDCHIRPSPDATAC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * make sure that the arguments are reasonable
   *
   */


  /* check that the output exists */
  ASSERT( fcSegVec, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": fcSegVec" );
  ASSERT( fcSegVec->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": fcSegVec->data" );
  ASSERT( fcSegVec->data->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": fcSegVec->data->dat" );
  ASSERT( fcSegVec->data->data->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": fcSegVec->data->data->data" );

  /* check that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPSPH_ENULL,
      FINDCHIRPSPH_MSGENULL ": params" );

  /* check that the workspace vectors exist */
  ASSERT( params->ampVec, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->ampVec->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  ASSERT( params->wVec, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->wVec->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  ASSERT( params->wtildeVec, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->wtildeVec->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  ASSERT( params->tmpltPowerVec, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->tmpltPowerVec->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the fft plans exist */
  ASSERT( params->fwdPlan, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->invPlan, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the parameter values are reasonable */
  ASSERT( params->fLow >= 0, status,
      FINDCHIRPSPH_EFLOW, FINDCHIRPSPH_MSGEFLOW );
  ASSERT( params->dynRange > 0, status,
      FINDCHIRPSPH_EDYNR, FINDCHIRPSPH_MSGEDYNR );

  /* check that the input exists */
  ASSERT( dataSegVec, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": dataSegVec" );
  ASSERT( dataSegVec->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": dataSegVec->data" );
  ASSERT( dataSegVec->data->chan, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": dataSegVec->data->chan" );
  ASSERT( dataSegVec->data->chan->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": dataSegVec->data->chan->data" );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  if ( params->approximant != FindChirpSP )
  {
    ABORT( status, FINDCHIRPSPH_EMAPX, FINDCHIRPSPH_MSGEMAPX );
  }


  /*
   *
   * set up local segment independent pointers
   *
   */


  w          = params->wVec->data;
  amp        = params->ampVec->data;
  wtilde     = params->wtildeVec->data;
  tmpltPower = params->tmpltPowerVec->data;

  /* allocate memory to store some temporary info for the
     continous chisq test */
  fcSeg        = &(fcSegVec->data[0]);
  fftVec = XLALCreateCOMPLEX8Vector( fcSeg->data->data->length );

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

    outputData   = fcSeg->data->data->data;
    dataPower    = fcSeg->dataPower->data;

    ASSERT( params->wtildeVec->length == fcSeg->data->data->length, status,
        FINDCHIRPSPH_EMISM, FINDCHIRPSPH_MSGEMISM );


    /* store the waveform approximant in the data segment */
    fcSeg->approximant = params->approximant;


    /*
     *
     * compute htilde and store in fcSeg
     *
     */


    LALForwardRealFFT( status->statusPtr, fcSeg->data->data,
        dataVec, params->fwdPlan );
    CHECKSTATUSPTR( status );

    /* compute strain */
    for ( k = 0; k < fcSeg->data->data->length; ++k )
    {
      REAL4 p = outputData[k].re;
      REAL4 q = outputData[k].im;
      REAL4 x = resp[k].re * params->dynRange;
      REAL4 y = resp[k].im * params->dynRange;

      outputData[k].re =  p*x - q*y;
      outputData[k].im =  p*y + q*x;
    }


    /*
     *
     * compute inverse power spectrum
     *
     */


    /* set low frequency cutoff inverse power spectrum */
    cut = params->fLow / dataSeg->spec->deltaF > 1 ?
      params->fLow / dataSeg->spec->deltaF : 1;
    LALSnprintf( infoMsg, sizeof(infoMsg)/sizeof(*infoMsg),
        "low frequency cut off index = %d\n", cut );
    LALInfo( status, infoMsg );

    /* set inverse power spectrum to zero */
    memset( wtilde, 0, params->wtildeVec->length * sizeof(COMPLEX8) );

    /* compute inverse of S_v */
    for ( k = cut; k < params->wtildeVec->length; ++k )
    {
      if ( spec[k] == 0 )
      {

        ABORT( status, FINDCHIRPSPH_EDIVZ, FINDCHIRPSPH_MSGEDIVZ );
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
        ABORT( status, FINDCHIRPSPH_EDIVZ, FINDCHIRPSPH_MSGEDIVZ );
      }
      invmodsqResp = 1.0 / modsqResp;
      wtilde[k].re *= invmodsqResp;
    }


    /*
     *
     * compute segment normalisation, outputData, point fcSeg at data segment
     *
     */



    for ( k = 0; k < cut; ++k )
    {
      outputData[k].re = 0.0;
      outputData[k].im = 0.0;
    }

    for ( k = 0; k < cut; ++k )
    {
      fftVec->data[k].re = 0.0;
      fftVec->data[k].im = 0.0;
    }


    memset( tmpltPower, 0, params->tmpltPowerVec->length * sizeof(REAL4) );
    memset( fcSeg->segNorm->data, 0, fcSeg->segNorm->length * sizeof(REAL4) );

    fcSeg->tmpltPowerVec = params->tmpltPowerVec;

    segNormSum = 0.0;
    for ( k = 1; k < fcSeg->data->data->length; ++k )
    {
      tmpltPower[k] = amp[k] * amp[k] * wtilde[k].re;
      segNormSum += tmpltPower[k];
      fcSeg->segNorm->data[k] = segNormSum;
    }

    /*  Compute whitened data for continous chisq test */
    for ( k = 0; k < fcSeg->data->data->length; ++k )
    {
      fftVec->data[k].re  = outputData[k].re * sqrt( wtilde[k].re );
      fftVec->data[k].im  = outputData[k].im * sqrt( wtilde[k].re );
    }

    /* get the whitened time series */
    LALReverseRealFFT( status->statusPtr, dataPower, fftVec,
          params->invPlan );
    dataPower->data[0] = 0;

    /* compute the cumulative power used for the continous
       chisq test */
    for ( k = 1; k < dataPower->length; k++ )
      {
      dataPower->data[k] =
        dataPower->data[k-1] +
        dataPower->data[k] * dataPower->data[k];
      }

    /* hard wired to quarter segment !! */
    startIX = floor(1.0/4.0 * (REAL4) dataPower->length + 0.5);
    endIX = floor(3.0/4.0 * (REAL4) dataPower->length + 0.5);
    /* compute the total power in the uncorrupted data */
    dataPower->data[dataPower->length - 1 ] = 2.0 *
      (dataPower->data[endIX] - dataPower->data[startIX]);
    for ( k = cut; k < fcSeg->data->data->length; ++k )
    {
      outputData[k].re  *= wtilde[k].re * amp[k];
      outputData[k].im  *= wtilde[k].re * amp[k];
    }

    /* set output frequency series parameters */
    strncpy( fcSeg->data->name, dataSeg->chan->name, LALNameLength );

    fcSeg->data->epoch.gpsSeconds      = dataSeg->chan->epoch.gpsSeconds;
    fcSeg->data->epoch.gpsNanoSeconds  = dataSeg->chan->epoch.gpsNanoSeconds;

    fcSeg->data->f0     = dataSeg->chan->f0;
    fcSeg->data->deltaF = 1.0 /
      ( (REAL8) dataSeg->chan->data->length * dataSeg->chan->deltaT ) ;

    fcSeg->deltaT       = dataSeg->chan->deltaT;
    fcSeg->number       = dataSeg->number;
    fcSeg->analyzeSegment = dataSeg->analyzeSegment;

    /* store low frequency cutoff and invSpecTrunc in segment */
    fcSeg->fLow         = params->fLow;
    fcSeg->invSpecTrunc = params->invSpecTrunc;


  } /* end loop over data segments */


  /* Find the min power from the whitened time series */
  /* For the continuous chisq test */
  fcSeg = &(fcSegVec->data[0]);
  PSDsum = fcSeg->dataPower->data->data[fcSeg->dataPower->data->length - 1 ];

  for ( i = 1; i < dataSegVec->length; ++i )
  {
    fcSeg = &(fcSegVec->data[i]);
    if
    (
    (fcSeg->dataPower->data->data[fcSeg->dataPower->data->length - 1 ] < PSDsum)    &&
    (fcSeg->dataPower->data->data[fcSeg->dataPower->data->length - 1 ] > 0)
    ||
    PSDsum == 0
    )
    {
      PSDsum = fcSeg->dataPower->data->data[fcSeg->dataPower->data->length - 1 ];
    }

  }
  /* reset each dataPower's last element to the min power */
/*  sortFlag = rint(dataSegVec->length / 2 + 1);*/
  sortFlag = 0;
  for ( i = 0; i < dataSegVec->length; ++i )
  {
    fcSeg = &(fcSegVec->data[i]);
    fcSeg->dataPower->data->data[fcSeg->dataPower->data->length - 1 ] = PSDsum;

  }

  /* clean up the data used for the continous chisq test */
  XLALDestroyCOMPLEX8Vector( fftVec );

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
