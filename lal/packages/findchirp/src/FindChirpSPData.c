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

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpSP.h>


NRCSID (FINDCHIRPSPDATAC, "$Id$");

void
LALFindChirpSPDataInit (
    LALStatus                  *status,
    FindChirpSPDataParams     **output,
    FindChirpInitParams        *params
    )
{
  FindChirpSPDataParams        *dataParamPtr;
  REAL4                        *amp;
  const REAL4                   exponent = -7.0/6.0;
  UINT4                         k;

  INITSTATUS( status, "LALFindChirpSPDataInit", FINDCHIRPSPDATAC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( output, status, FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( !*output, status, FINDCHIRPSP_ENNUL, FINDCHIRPSP_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT (params, status, FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL);

  /* make sure that the number of points in a segment is positive */
  ASSERT (params->numPoints > 0, status, 
      FINDCHIRPSP_ENUMZ, FINDCHIRPSP_MSGENUMZ);


  /*
   *
   * allocate memory for the FindChirpSPDataParams
   *
   */


  /* create the output structure */
  dataParamPtr = *output = (FindChirpSPDataParams *)
    LALMalloc( sizeof(FindChirpSPDataParams) );
  ASSERT( dataParamPtr, status, FINDCHIRPSP_EALOC, FINDCHIRPSP_MSGEALOC );
  memset( dataParamPtr, 0, sizeof(FindChirpSPDataParams) );

  /* should need this because of the memset above */
  /* set contents to reasonable values */
  /* dataParamPtr->ampVec        = NULL; */
  /* dataParamPtr->fwdPlan       = NULL; */
  /* dataParamPtr->invPlan       = NULL; */
  /* dataParamPtr->vVec          = NULL; */
  /* dataParamPtr->wVec          = NULL; */
  /* dataParamPtr->wtildeVec     = NULL; */
  /* dataParamPtr->tmpltPowerVec = NULL; */
  /* dataParamPtr->deltaT        = 0; */
  /* dataParamPtr->fLow          = 0; */
  /* dataParamPtr->dynRange      = 0; */
  /* dataParamPtr->invSpecTrunc  = 0; */
  

  /*
   *
   * allocate and fill vector for exponent of amplitude
   *
   */


  LALCreateVector( status->statusPtr, &dataParamPtr->ampVec, 
        (params->numPoints)/2 + 1 );
  BEGINFAIL( status )
  {
    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  amp = dataParamPtr->ampVec->data;
  amp[0] = 0.0;

  for ( k = 1; k < dataParamPtr->ampVec->length; ++k )
    amp[k] = pow( ((REAL4) k / (REAL4)params->numPoints), exponent );


  /*
   *
   * create fft plans and workspace vectors
   *
   */


  /* foward fft plan */
  LALEstimateFwdRealFFTPlan (status->statusPtr, &dataParamPtr->fwdPlan, 
        params->numPoints);
  BEGINFAIL( status )
  {
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), status );

    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* inverse fft plan */
  LALEstimateInvRealFFTPlan( status->statusPtr, &dataParamPtr->invPlan, 
        params->numPoints );
  BEGINFAIL( status )
  {
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), status );

    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* workspace vector v */
  LALCreateVector( status->statusPtr, &dataParamPtr->vVec, 
      params->numPoints );
  BEGINFAIL( status )
  {
    TRY( LALDestroyRealFFTPlan (status->statusPtr, &dataParamPtr->invPlan ), status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), status );

    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* workspace vector w: time domain */
  LALCreateVector( status->statusPtr, &dataParamPtr->wVec, 
      params->numPoints );
  BEGINFAIL( status )
  {
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->vVec ), status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->invPlan ), status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), status );

    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* workspace vector w: freq domain */
  LALCCreateVector( status->statusPtr, &dataParamPtr->wtildeVec, 
      params->numPoints/2 + 1 );
  BEGINFAIL( status )
  {
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->wVec ), status ); 
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->vVec ), status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->invPlan ), status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), status );

    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );
  CHECKSTATUSPTR (status);
  
  /* template power vector */
  LALCreateVector( status->statusPtr, &dataParamPtr->tmpltPowerVec, 
      params->numPoints/2 + 1 );
  BEGINFAIL( status )
  {
    TRY( LALCDestroyVector( status->statusPtr, &dataParamPtr->wtildeVec), status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->wVec ), status ); 
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->vVec ), status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->invPlan ), status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), status );

    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );
  

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



void
LALFindChirpSPDataFinalize (
    LALStatus                  *status,
    FindChirpSPDataParams     **output
    )
{
  FindChirpSPDataParams        *dataParamPtr;

  INITSTATUS( status, "LALFindChirpSPDataFinalize", FINDCHIRPSPDATAC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  /* make sure handle is non-null and points to a non-null pointer */
  ASSERT( output, status, FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( *output, status, FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );

  
  /*
   *
   * destroy fft plans and workspace vectors
   *
   */


  /* local pointer to structure */
  dataParamPtr = *output;
  
  LALDestroyRealFFTPlan (status->statusPtr, &dataParamPtr->fwdPlan);
  CHECKSTATUSPTR (status);

  LALDestroyRealFFTPlan (status->statusPtr, &dataParamPtr->invPlan);
  CHECKSTATUSPTR (status);

  LALDestroyVector (status->statusPtr, &dataParamPtr->vVec);
  CHECKSTATUSPTR (status);

  LALDestroyVector (status->statusPtr, &dataParamPtr->wVec);
  CHECKSTATUSPTR (status);

  LALCDestroyVector (status->statusPtr, &dataParamPtr->wtildeVec);
  CHECKSTATUSPTR (status);

  LALDestroyVector (status->statusPtr, &dataParamPtr->tmpltPowerVec);
  CHECKSTATUSPTR (status);


  /*
   *
   * destroy vector for exponent of amplitude
   *
   */


  LALDestroyVector (status->statusPtr, &dataParamPtr->ampVec);
  CHECKSTATUSPTR (status);


  /*
   *
   * free memory for the FindChirpSPDataParams
   *
   */


  LALFree (dataParamPtr);
  *output = NULL;


  DETATCHSTATUSPTR (status);
  RETURN (status);
}



void
LALFindChirpSPData (
    LALStatus                  *status,
    FindChirpSegmentVector     *fcSegVec,
    DataSegmentVector          *dataSegVec,
    FindChirpSPDataParams      *params
    )
{
  UINT4                 i, j, k; 
  UINT4                 cut;

  REAL4                *v;
  REAL4                *w;
  REAL4                *amp;
  COMPLEX8             *wtilde;
  REAL4                *tmpltPower;
  
  REAL4                *data;
  REAL4                *spec;
  COMPLEX8             *resp;
  
  COMPLEX8             *outputData;

  UINT4                *chisqBin;
  UINT4                 numChisqBins;
  UINT4                 chisqPt;
  REAL4                 increment;
  REAL4                 nextBin;
  REAL4                 partSum;
  
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
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( fcSegVec->data, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( fcSegVec->data->data, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( fcSegVec->data->data->data, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  
  /* check that the workspace vectors exist */
  ASSERT( params->ampVec, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( params->ampVec->data, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );

  ASSERT( params->vVec, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( params->vVec->data, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );

  ASSERT( params->wVec, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( params->wVec->data, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );

  ASSERT( params->wtildeVec, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( params->wtildeVec->data, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );

  ASSERT( params->tmpltPowerVec, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( params->tmpltPowerVec->data, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );

  /* check that the fft plans exist */
  ASSERT( params->fwdPlan, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( params->invPlan, status, 
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );

  /* check that the parameter values are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPSP_EDELT, FINDCHIRPSP_MSGEDELT );
  ASSERT( params->fLow >= 0, status,
      FINDCHIRPSP_EFLOW, FINDCHIRPSP_MSGEFLOW );
  ASSERT( params->dynRange > 0, status,
      FINDCHIRPSP_EDYNR, FINDCHIRPSP_MSGEDYNR );
  
  /* check that the input exists */
  ASSERT( dataSegVec, status,
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( dataSegVec->data, status,
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( dataSegVec->data->real4Data, status,
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );
  ASSERT( dataSegVec->data->real4Data->data, status,
      FINDCHIRPSP_ENULL, FINDCHIRPSP_MSGENULL );


  /*
   *
   * set up local segment independent pointers
   *
   */


  v          = params->vVec->data;
  w          = params->wVec->data;
  amp        = params->ampVec->data;
  wtilde     = params->wtildeVec->data;
  tmpltPower = params->tmpltPowerVec->data;


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

    data         = dataSeg->real4Data->data->data;
    spec         = dataSeg->spec->data->data;
    resp         = dataSeg->resp->data->data;

    outputData   = fcSeg->data->data->data;

    chisqBin     = fcSeg->chisqBinVec->data;
    numChisqBins = fcSeg->chisqBinVec->length - 1;

    ASSERT( params->wtildeVec->length == fcSeg->data->data->length, status,
        FINDCHIRPSP_EMISM, FINDCHIRPSP_MSGEMISM );


    /*
     *
     * compute htilde and store in wtilde
     *
     */


    /* fft the IFO data */
    for ( j = 0; j < params->vVec->length; ++j )
    {
      v[j] = (REAL4) data[j];
    }

    LALFwdRealFFT( status->statusPtr, fcSeg->data->data, 
        params->vVec, params->fwdPlan );
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

    /* set inverse power spectrum to zero */
    memset( wtilde, 0, params->wtildeVec->length * sizeof(COMPLEX8) );

    /* compute inverse power spectrum above fLow */
    for ( k = cut; k < params->wtildeVec->length; ++k )
    {
      REAL4 respRe = resp[k].re * params->dynRange;
      REAL4 respIm = resp[k].im * params->dynRange;
      REAL4 modsqResp = respRe * respRe + respIm * respIm;
      wtilde[k].re = spec[k] * modsqResp;
    }

    for ( k = cut; k < params->wtildeVec->length; ++k )
    {
      wtilde[k].re = 1.0 / wtilde[k].re;
    }


    /*
     *
     * compute segment normalisation, outputData, point fcSeg at data segment
     *
     */


    fcSeg->segNorm = 0.0;

    for ( k = 0; k < cut; ++k )
    {
      outputData[k].re = 0.0;
      outputData[k].im = 0.0;
    }

    memset( tmpltPower, 0, params->tmpltPowerVec->length * sizeof(REAL4) );

    for ( k = 1; k < fcSeg->data->data->length; ++k )
    {
      tmpltPower[k]   = amp[k] * amp[k] * wtilde[k].re;
      fcSeg->segNorm += tmpltPower[k];
    }

    for ( k = cut; k < fcSeg->data->data->length; ++k )
    {
      wtilde[k].re      *= amp[k];
      outputData[k].re  *= wtilde[k].re;
      outputData[k].im  *= wtilde[k].re;
    }

    /* set output frequency series parameters */
    fcSeg->data->epoch  = dataSeg->data->epoch;
    fcSeg->data->f0     = dataSeg->data->f0;
    fcSeg->data->deltaF = 1.0 / 
      ( (REAL8) dataSeg->data->data->length * dataSeg->data->deltaT ) ;
    fcSeg->deltaT       = dataSeg->data->deltaT;
    fcSeg->number       = dataSeg->number;;


    /*
     *
     * calculate the chisq bins for the segment and template
     *
     */


    increment = fcSeg->segNorm / (REAL4) numChisqBins;
    nextBin   = increment;
    chisqPt   = 0;
    partSum   = 0.0;

    /* calulate the frequencies of the chi-squared bin boundaries */
    chisqBin[chisqPt++] = 0;

    for ( k = 1; k < fcSeg->data->data->length; ++k ) 
    {
      partSum += tmpltPower[k];
      if ( partSum >= nextBin ) 
      {
        chisqBin[chisqPt++] = k;
        nextBin += increment;
        if ( chisqPt == numChisqBins ) break;
      }
    }
    chisqBin[numChisqBins] = fcSeg->data->data->length;

  } /* end loop over data segments */


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
