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
 *-----------------------------------------------------------------------
 */

/**

\author Brown, D. A.
\file
\ingroup FindChirpSP_h

\brief Provides functions to condition the input data from the interferometer
to a form that can be used by the <tt>FindChirpFilter()</tt> function.

At the present time this also includes the template independent part of the
stationary phase filter.

\section sec_fcd_desc Description

<dl>
<dt><tt> LALFindChirpSPDataInit()</tt></dt><dd> takes as input the address
of a structure of type \c FindChirpInitParams containing the correct
values to intialize a search. It creates a structure of type
\c FindChirpSPDataParams as described above and returns its address.</dd>

<dt><tt>LALFindChirpSPDataFinalize()</tt></dt><dd> takes as the address
of a structure of type \c FindChirpSPDataParams destroys this
structure and sets the address to NULL.</dd>

<dt><tt>LALFindChirpSPData()</tt></dt><dd> conditions the interferometer data
as described by the algorithm below.</dd>

<dt><tt>LALFindChirpBCVData()</tt></dt><dd> conditions the interferometer data
as described by the algorithm below.</dd>
</dl>

\subsection sec_fcd_alg Algorithm for SP templates

The <tt>LALFindChirpSPData()</tt> function takes as input three quantities
<ol>
<li> An uncallibrated input data channel \f$v_j\f$ (typically LSC-AS_Q).</li>
<li> A one sided power spectral density of the input data \f$S_v(|f_k|)\f$.</li>
<li> The frequency domain response function of the input data channel
\f$R(f_k)\f$ which is used to convert from an uncallibrated input data into
strain.</li>
</ol>

The input parameters also contain a dynmaic range scaling parameter
\c dynRange. This is used to keep the quantities in the range of
\c REAL4 and cancells from the filter output. It is typically \f$2^{69}\f$
for the LIGO channel \c AS_Q.

The discrete low frequency cutoff is computed from \c fLow by
\f{equation}{
k_{\mathrm{low}} = \frac{\mathtt{fLow}}{\Delta f}.
\f}
\f$k_{\mathrm{low}}\f$ is set to unity if this quantity is less than one.

\heading{Computation of strain and inverse power spectrum}

The uncallibrated input data channel \f$v_j\f$ is Fourier transformed into the
frequency domain to obtain \f$\tilde{v}_k\f$. This is then converted into strain
\f$\tilde{h}_k\f$ by computing
\f{equation}{
\tilde{h}_k = \mathtt{dynRange} \times R(f_k) \tilde{v}_k.
\f}
The inverse power spectrum \f$S^{-1}_v(|f_k|)\f$ is computed between the low
frequency cutoff \f$f_{\mathrm{low}}\f$ and the Nyquist frequency. Below
\f$f_{\mathrm{low}}\f$, the inverse power spectrum is set to zero. If the low
frequency cutoff is set to \f$0\f$, then the DC component of the spectrum is
set to zero.

\heading{Truncation of the Inverse Power Spectrum in the Time Domain}

Recall that the FFT we use to compute the match filter treats the data as
being periodic and that we had to ignore part of the filter output that was
corrupted due to wraparound of the filter output from the chirp signal.

As well as the chirp, we are also filtering the data against the inverse power
spectrum. The chirp has a duration that is typically much less than then
length of the data segment. It only corrupts a region that is the length of
the chirp at the start of the data segment. However, in the time domain the
inverse power spectrum is the same length of the data segmemt. This will cause
the filter output to be corrupted for all times, as it is non-zero over the
entire data segment.

To prevent this, we truncate the inverse power spectrum to a specified length
in the time domain. This has the effect of smoothing out the high \f$Q\f$ features
which and restructing the corruption of the filter to the part of the data
segment where the power spectrum is non-zero. These regions can then be
ignored when searching for chirps in the filter output.

The parameter that controls the duration of the power spectrum in the time
domain is \c invSpecTrunc and the algorithm used to perform the
truncation is as follows:
<ol>
<li> Compute the square root of the inverse power spectrum,
\f$\sqrt{S^{-1}_v(|f_k|)}\f$.</li>
<li> Set the Nyquist, \f$k = N/2\f$ and DC \f$k = 0\f$ components of this to zero.</li>
<li> Inverse FFT to to obtain the time domain inverse PSD of length \f$N\f$
points.</li>
<li> Zero the spectrum between the points \f$j = \mathtt{invSpecTrunc}/2\f$ and
\f$j = N - \mathtt{invSpecTrunc}/2\f$. This sets the length of the inverse
spectrum in the time domain to be \c invSpecTrunc points.</li>
<li> FFT the time domain quantity back to the frequency domain.</li>
<li> Divide by \f$N\f$ so that to recover the quantity before the inverse FFT.</li>
<li> Square this quantity to recover \f$S^{-1}_v(|f_k|)\f$.</li>
<li> Set the Nyqist and DC frequencies to zero and zero the inverse power
spectrum below \f$f_{\mathrm{low}}\f$.</li>
</ol>

The strain inverse power spectral density is then computed by
\f{equation}{
S^{-1}_h(|f_k|) = \frac{1}{\left|\mathtt{dynRange} \times R(f_k)\right|^2}
 S^{-1}_v(|f_k|).
\f}

\heading{Output Data}
The quantity \c segNorm is computed by
\f{equation}{
\mathtt{segNorm} =
   \sum_{k = k_{\mathrm{low}}}^{N/2} \frac{k^{-\frac{7}{3}}}
   {\left|\mathtt{dynRange}\times R(f_k)\right|^2 S_v(|f_k|)}.
\f}

The output data is given by
\f{equation}{
\mathtt{outputData[k]} =
\frac{k^{-7/6} \times \mathtt{dynRange} \times R(f_k)\tilde{v}_k}
{\left|\mathtt{dynRange} \times R(f_k)\right|^2 S_v(|f_k|)}
\f}
and is stored in the \c FindChirpSegmentVector structure. Note the
quantity \f$k^{-\frac{7}{6}}\f$ which is specific to the stationary phase chirps
used.

\heading{Calculation of the \f$\chi^2\f$ Bins}

If a \f$\chi^2\f$ veto is requested, the bin boundaries for the veto are computed
at this point. The indices \f$k\f$ in the frequency domain that divide the power
in the quantity
\f{equation}{
\frac{k^{-\frac{7}{3}}}
{\left|\mathtt{dynRange}\times R(f_k)\right|^2 S_v(|f_k|)}.
\f}
into equal intervals are stored in the array \c chisqBin.

\subsection sec_fcd_bcv Algorithm for BCV templates

The <tt>LALFindChirpBCVData()</tt> function takes as input...

\heading{Uses}
\code
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
LALCreateForwardRealFFTPlan()
LALDestroyRealFFTPlan()
LALCCreateVector()
LALCDestroyVector()
LALForwardRealFFT()
LALReverseRealFFT()
\endcode

\heading{Notes}

*/

#define LAL_USE_OLD_COMPLEX_STRUCTS
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>
#include <math.h>

void
LALFindChirpSPData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    )

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
  COMPLEX8Vector       *fftVec = NULL;
  FindChirpSegment     *fcSeg;
  DataSegment          *dataSeg;
  INITSTATUS(status);
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
      REAL4 p = crealf(outputData[k]);
      REAL4 q = cimagf(outputData[k]);
      REAL4 x = crealf(resp[k]) * params->dynRange;
      REAL4 y = cimagf(resp[k]) * params->dynRange;

      outputData[k].realf_FIXME =  p*x - q*y;
      outputData[k].imagf_FIXME =  p*y + q*x;
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

        ABORT( status, FINDCHIRPSPH_EDIVZ, FINDCHIRPSPH_MSGEDIVZ );
      }
      wtilde[k].realf_FIXME = 1.0 / spec[k];
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
        wtilde[k].realf_FIXME = sqrt( crealf(wtilde[k]) );
      }

      /* set nyquist and dc to zero */
      wtilde[params->wtildeVec->length - 1].realf_FIXME = 0.0;
      wtilde[0].realf_FIXME                             = 0.0;

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
          wtilde[k].realf_FIXME *= norm;
          wtilde[k].realf_FIXME *= crealf(wtilde[k]);
          wtilde[k].imagf_FIXME = 0.0;
        }
      }

      /* set nyquist and dc to zero */
      wtilde[params->wtildeVec->length - 1].realf_FIXME = 0.0;
      wtilde[0].realf_FIXME                             = 0.0;
    }

    /* set inverse power spectrum below cut to zero */
    memset( wtilde, 0, cut * sizeof(COMPLEX8) );

    /* convert from S_v to S_h */
    for ( k = cut; k < params->wtildeVec->length; ++k )
    {
      REAL4 respRe = crealf(resp[k]) * params->dynRange;
      REAL4 respIm = cimagf(resp[k]) * params->dynRange;
      REAL4 modsqResp = (respRe * respRe + respIm * respIm);
      REAL4 invmodsqResp;
      if ( modsqResp == 0 )
      {
        ABORT( status, FINDCHIRPSPH_EDIVZ, FINDCHIRPSPH_MSGEDIVZ );
      }
      invmodsqResp = 1.0 / modsqResp;
      wtilde[k].realf_FIXME *= invmodsqResp;
    }


    /*
     *
     * compute segment normalisation, outputData, point fcSeg at data segment
     *
     */



    for ( k = 0; k < cut; ++k )
    {
      outputData[k].realf_FIXME = 0.0;
      outputData[k].imagf_FIXME = 0.0;
    }

    for ( k = 0; k < cut; ++k )
    {
      fftVec->data[k].realf_FIXME = 0.0;
      fftVec->data[k].imagf_FIXME = 0.0;
    }


    memset( tmpltPower, 0, params->tmpltPowerVec->length * sizeof(REAL4) );
    memset( fcSeg->segNorm->data, 0, fcSeg->segNorm->length * sizeof(REAL4) );

    fcSeg->tmpltPowerVec = params->tmpltPowerVec;

    segNormSum = 0.0;
    for ( k = 1; k < fcSeg->data->data->length; ++k )
    {
      tmpltPower[k] = amp[k] * amp[k] * crealf(wtilde[k]);
      segNormSum += tmpltPower[k];
      fcSeg->segNorm->data[k] = segNormSum;
    }

    /*  Compute whitened data for continous chisq test */
    for ( k = 0; k < fcSeg->data->data->length; ++k )
    {
      fftVec->data[k].realf_FIXME  = crealf(outputData[k]) * sqrt( crealf(wtilde[k]) );
      fftVec->data[k].imagf_FIXME  = cimagf(outputData[k]) * sqrt( crealf(wtilde[k]) );
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
      outputData[k].realf_FIXME  *= crealf(wtilde[k]) * amp[k];
      outputData[k].imagf_FIXME  *= crealf(wtilde[k]) * amp[k];
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
    ((fcSeg->dataPower->data->data[fcSeg->dataPower->data->length - 1 ] < PSDsum)    &&
    (fcSeg->dataPower->data->data[fcSeg->dataPower->data->length - 1 ] > 0))
    ||
    PSDsum == 0
    )
    {
      PSDsum = fcSeg->dataPower->data->data[fcSeg->dataPower->data->length - 1 ];
    }

  }
  /* reset each dataPower's last element to the min power */
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
