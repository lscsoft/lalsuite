/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpBCVSpinData.c
 *
 * Author: Brown D. A., Spinning BCV-Modifications: Jones, G
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpBCVSpinDataCV">
Author: Brown, D. A., Spinning BCV-Modifications: Jones, G.
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpBCVSpinData.c}}
\label{ss:FindChirpBCVSpinData.c}

Provides functions to condition data prior to filtering with spinning BCV
detection templates.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpBCVDataCP}
\idx{LALFindChirpBCVData()}

The function \texttt{LALFindChirpBCVSpinData()} constions the data
as described by the algorithm below.

\subsubsection*{Algorithm}

Blah.

\subsubsection*{Uses}
\begin{verbatim}
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\end{verbatim}

\subsubsection*{Notes}

\vfill{\footnotesize\input{FindChirpBCVSpinDataCV}}
</lalLaTeX> 
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpBCVSpin.h>


NRCSID (FINDCHIRPBCVSPINDATAC, "$Id$");

/* <lalVerbatim file="FindChirpBCVSpinDataCP"> */
void
LALFindChirpBCVSpinData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpDataParams        *params
    )
/* </lalVerbatim> */
{
  UINT4                 i, k;
  UINT4                 cut;

  REAL4                *w;
  REAL4                *amp;
  REAL4                *ampBCV; 
  COMPLEX8             *wtilde;
  REAL4		       *tmpltPower;
  REAL4		       *tmpltPowerBCV;

  REAL4Vector          *dataVec;
  REAL4                *spec;
  COMPLEX8             *resp;

  COMPLEX8             *outputData;
  COMPLEX8             *outputDataBCV;

  UINT4                *chisqBin    = NULL;
  UINT4                *chisqBinBCV = NULL;
  UINT4                 numChisqBins;
  UINT4                 chisqPt;
  REAL4                 increment;
  REAL4                 nextBin;
  REAL4                 partSum;
  REAL4                 Power  = 0.0 ;
  REAL4                 PowerBCV = 0.0 ;

  FindChirpSegment     *fcSeg;
  DataSegment          *dataSeg;

  FILE                 *fpdataVec;
  FILE                 *fpspec;
  FILE		       *fprespRe;
  FILE                 *fprespIm;
  FILE                 *fpStrainRe;
  FILE                 *fpStrainIm;

  FILE                 *fpDataIn = NULL;


  /*declaration*/
  INITSTATUS( status, "LALFindChirpBCVSpinData", FINDCHIRPBCVSPINDATAC );
  ATTATCHSTATUSPTR( status );

  /*
   *
   * make sure that the arguments are reasonable
   * Identical to Eirini's except for aproximant term
   *
   */

/*  fprintf ( stdout, "just before assert statements of LALFindChirpBCVSpinData \n");*/

  /* check that the output exists */
  ASSERT( fcSegVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec" );
  ASSERT( fcSegVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec->data" );
  ASSERT( fcSegVec->data->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec->data->data" );
  ASSERT( fcSegVec->data->data->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec->data->data->data" );
 /* ASSERT( fcSegVec->data->dataBCV, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec->data->dataBCV" );
  ASSERT( fcSegVec->data->dataBCV->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": fcSegVec->data->dataBCV->data" );*/

  /* check that the parameter structure exists */
  ASSERT( params, status, 
	  FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL 
	  ": params" );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  ASSERT( params->approximant == BCVSpin, status, 
      FINDCHIRPBCVSPINH_EMAPX, FINDCHIRPBCVSPINH_MSGEMAPX );
  
  /* check that the workspace vectors exist */
  ASSERT( params->ampVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->ampVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
 /*ASSERT( params->ampVecBCV, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL ); 
   ASSERT( params->ampVecBCV->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );*/
ASSERT( params->ampVecBCVSpin1, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
   ASSERT( params->ampVecBCVSpin1->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
ASSERT( params->ampVecBCVSpin2, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
   ASSERT( params->ampVecBCVSpin2->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  ASSERT( params->wVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->wVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  ASSERT( params->wtildeVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->wtildeVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  ASSERT( params->tmpltPowerVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->tmpltPowerVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
/*  ASSERT( params->tmpltPowerVecBCV, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->tmpltPowerVecBCV->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );*/


  /* check that the fft plans exist */
  ASSERT( params->fwdPlan, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );
  ASSERT( params->invPlan, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL );

  /* check that the parameter values are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPBCVSPINH_EDELT, FINDCHIRPBCVSPINH_MSGEDELT );
  ASSERT( params->fLow >= 0, status,
      FINDCHIRPBCVSPINH_EFLOW, FINDCHIRPBCVSPINH_MSGEFLOW );
  ASSERT( params->dynRange > 0, status,
      FINDCHIRPBCVSPINH_EDYNR, FINDCHIRPBCVSPINH_MSGEDYNR );

  /* check that the input exists */
  ASSERT( dataSegVec, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": dataSegVec" );
  ASSERT( dataSegVec->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": dataSegVec->data" );
  ASSERT( dataSegVec->data->chan, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": dataSegVec->data->chan" );
  ASSERT( dataSegVec->data->chan->data, status,
      FINDCHIRPBCVSPINH_ENULL, FINDCHIRPBCVSPINH_MSGENULL
      ": dataSegVec->data->chan->data" );
  
  fpdataVec = fopen ("dataVec.dat","w");
  fpspec    = fopen ("spec.dat","w");
  fprespRe  = fopen ("respRe.dat","w");
  fprespIm  = fopen ("respIm.dat","w");
  fpStrainRe= fopen ("StrainRe.dat","w");
  fpStrainIm= fopen ("StrainIm.dat","w");

  fpDataIn = fopen ("DataIn.dat","w");





 /*fprintf ( stdout, "just after assert statements of LALFindChirpBCVSpinData \n");*/

  /*
   *
   * set up local segment independent pointers
   * will I need the last two? maybe for chisq
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

 /*fprintf ( stdout, "just before loop over data segments in LALFindChirpBCVSpinData \n");*/

fprintf (stdout, "dataSegVec->length %d\n", dataSegVec->length);


  for ( i = 0; i < dataSegVec->length; ++i )
  {

/*fprintf ( stdout, "just inside loop over data segments in LALFindChirpBCVSpinData \n");*/

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


/*    outputDataBCV = fcSeg->dataBCV->data->data;*/

/*fprintf ( stdout, "just after sep. dep. pointers in LALFindChirpBCVSpinData \n");*/


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
        FINDCHIRPBCVSPINH_EMISM, FINDCHIRPBCVSPINH_MSGEMISM );

/*fprintf ( stdout, "just after chisq in LALFindChirpBCVSpinData \n");*/


    /* store the waveform approximant in the data segment */
    fcSeg->approximant = BCVSpin;


  for ( i = 0; i < dataVec->length; ++i )
       {
        fprintf (fpDataIn, "%d\t%e\n",  i, dataVec->data[i]);
       }


/*fprintf ( stdout, "just after approx = BCVSpin in LALFindChirpBCVSpinData \n");*/
 

    /*
     *
     * compute htilde and store in fcSeg
     *
     */


    LALForwardRealFFT( status->statusPtr, fcSeg->data->data,
        dataVec, params->fwdPlan );
    CHECKSTATUSPTR( status );
/*    LALForwardRealFFT( status->statusPtr, fcSeg->dataBCV->data, 
        dataVec, params->fwdPlan );                             
    CHECKSTATUSPTR( status );*/                                   

fprintf ( stdout, "fcSeg->data->data->length  %d\n", fcSeg->data->data->length);


  /* compute strain */
  for ( k = 0; k < fcSeg->data->data->length; ++k )
  {
    REAL4 p = outputData[k].re;
    REAL4 q = outputData[k].im;
    REAL4 x = resp[k].re * params->dynRange;
    REAL4 y = resp[k].im * params->dynRange;

   
    fprintf (fprespRe, "%d\t%e\n",k,resp[k].re);
    fprintf (fprespIm, "%d\t%e\n",k,resp[k].im);


    outputData[k].re =  p*x - q*y;
    outputData[k].im =  p*y + q*x;
  
    fprintf (fpStrainRe, "%d\t%e\n",k,outputData[k].re);
    fprintf (fpStrainIm, "%d\t%e\n",k,outputData[k].im);


  }

 /*fprintf ( stdout, "just after htilde etc. calc. of LALFindChirpBCVSpinData \n");*/


    /*
     *
     * compute inverse power spectrum
     *
     */


    /* set low frequency cutoff inverse power spectrum */
    cut = params->fLow / dataSeg->spec->deltaF > 1 ?
      params->fLow / dataSeg->spec->deltaF : 1;

    /* set inverse power spectrum to zero */
    memset( wtilde, 0, params->wtildeVec->length * sizeof(COMPLEX8) );

    fprintf (stdout, "cut %d\n", cut);
    fprintf (stdout, "params->wtildeVec->length %d\n", params->wtildeVec->length);

    /* compute inverse of S_v */
    for ( k = cut; k < params->wtildeVec->length; ++k )
    {
    
	 fprintf (fpspec, "%d\t%e\n",k,spec[k]);


      if ( spec[k] == 0 )
      {
        ABORT( status, FINDCHIRPBCVSPINH_EDIVZ, FINDCHIRPBCVSPINH_MSGEDIVZ );
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

      fprintf(stdout, "truncating wtilde!");
     
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
        ABORT( status, FINDCHIRPBCVSPINH_EDIVZ, FINDCHIRPBCVSPINH_MSGEDIVZ );
      }
      invmodsqResp = 1.0 / modsqResp;
   
   /*  fprintf(stdout, "invmodsqResp %e\n", invmodsqResp);*/

      wtilde[k].re *= invmodsqResp;
    }

/*fprintf ( stdout, "just after wtilde etc. calc. of LALFindChirpBCVSpinData \n");*/

  /*
     *
     * compute BCV normalisation parameters a1, b1 and b2, 
     * segment normalization, outputData, point fcSeg at data segment
     *
     */

/*
    fcSeg->a1 = 0.0;
    fcSeg->b1 = 0.0;
    fcSeg->b2 = 0.0;
    fcSeg->segNorm = 0.0;

    for ( k = 0; k < cut; ++k )
    {
      outputData[k].re = 0.0;
      outputData[k].im = 0.0;
      outputDataBCV[k].re = 0.0;
      outputDataBCV[k].im = 0.0;
    }

    memset( tmpltPower, 0, params->tmpltPowerVec->length * sizeof(REAL4) );
    memset( tmpltPowerBCV,0, params->tmpltPowerVecBCV->length * sizeof(REAL4) );
*/
    /* 
     * moments necessary for the calculation of
     * the BCVSpin normalization parameters
     */
    
#if 0
    for ( k = 1; k < fcSeg->data->data->length; ++k )
    {
      fcSeg->segNorm += amp[k] * amp[k] * wtilde[k].re ; /* for std-candle */
    }
#endif

    /*
     *
     * tmplt power calcs require moments etc, not calc yet
     *
     */

  

    /* set output frequency series parameters */
    strncpy( fcSeg->data->name, dataSeg->chan->name, LALNameLength );
   /* strncpy( fcSeg->dataBCV->name, dataSeg->chan->name, LALNameLength );*/

    fcSeg->data->epoch.gpsSeconds      = dataSeg->chan->epoch.gpsSeconds;
    fcSeg->data->epoch.gpsNanoSeconds  = dataSeg->chan->epoch.gpsNanoSeconds;
 /*   fcSeg->dataBCV->epoch.gpsSeconds     = dataSeg->chan->epoch.gpsSeconds;
    fcSeg->dataBCV->epoch.gpsNanoSeconds = dataSeg->chan->epoch.gpsNanoSeconds;*/

    fcSeg->data->f0     = dataSeg->chan->f0;
    fcSeg->data->deltaF = 1.0 /
      ( (REAL8) dataSeg->chan->data->length * dataSeg->chan->deltaT ) ;
  /*  fcSeg->dataBCV->f0     = dataSeg->chan->f0;
    fcSeg->dataBCV->deltaF = 1.0 /
      ( (REAL8) dataSeg->chan->data->length * dataSeg->chan->deltaT ) ;*/

    fcSeg->deltaT       = dataSeg->chan->deltaT;
    fcSeg->number       = dataSeg->number;

    /* store low frequency cutoff and invSpecTrunc in segment */
    fcSeg->fLow         = params->fLow;
    fcSeg->invSpecTrunc = params->invSpecTrunc;


  /* code */

  } /* end of loop over data segments */
 
fclose(fprespRe);
fclose(fprespIm);
fclose(fpspec);
fclose(fpdataVec);
fclose(fpStrainRe);
fclose(fpStrainIm);
fclose (fpDataIn);


  DETATCHSTATUSPTR( status );
  RETURN( status );

}
