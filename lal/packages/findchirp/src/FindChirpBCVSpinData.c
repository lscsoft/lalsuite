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

\input{FindChirpBCVSpinDataCDoc}

\vfill{\footnotesize\input{FindChirpBCVSpinDataCV}}
</lalLaTeX> 
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>


NRCSID (FINDCHIRPBCVSPINDATAC, "$Id$");

/*documenation later*/
void
LALFindChirpBCVSpinData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpSPDataParams      *params
    )

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



  /*declaration*/
  INITSTATUS( status, "LALFindChirpBCVSpinData", FINDCHIRPBCVSPINDATAC );
  ATTATCHSTATUSPTR( status );

  /*
   *
   * make sure that the arguments are reasonable
   * Identical to Eirini's except for aproximant term
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
      ": fcSegVec->data->data" );
  ASSERT( fcSegVec->data->data->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": fcSegVec->data->data->data" );
  ASSERT( fcSegVec->data->dataBCV, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": fcSegVec->data->dataBCV" );
  ASSERT( fcSegVec->data->dataBCV->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL
      ": fcSegVec->data->dataBCV->data" );

  /* check that the parameter structure exists */
  ASSERT( params, status, 
	  FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL 
	  ": params" );

  /* check that the parameter structure is set */
  /* to the correct waveform approximant       */
  ASSERT( params->approximant == BCVSpin, status, 
      FINDCHIRPSPH_EMAPX, FINDCHIRPSPH_MSGEMAPX );
  
  /* check that the workspace vectors exist */
  ASSERT( params->ampVec, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->ampVec->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->ampVecBCV, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL ); 
  ASSERT( params->ampVecBCV->data, status,
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
  ASSERT( params->tmpltPowerVecBCV, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->tmpltPowerVecBCV->data, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );


  /* check that the fft plans exist */
  ASSERT( params->fwdPlan, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );
  ASSERT( params->invPlan, status,
      FINDCHIRPSPH_ENULL, FINDCHIRPSPH_MSGENULL );

  /* check that the parameter values are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPSPH_EDELT, FINDCHIRPSPH_MSGEDELT );
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
  /*
   *
   * set up local segment independent pointers
   * will I need the last two
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
        FINDCHIRPSPH_EMISM, FINDCHIRPSPH_MSGEMISM );

    /* store the waveform approximant in the data segment */
    fcSeg->approximant = BCVSpin;


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
    REAL4 x = resp[k].re * params->dynRange;
    REAL4 y = resp[k].im * params->dynRange;

    outputData[k].re =  p*x - q*y;
    outputData[k].im =  p*y + q*x;
  }

  /* code */

  }
 
  DETATCHSTATUSPTR( status );
  RETURN( status );

}
