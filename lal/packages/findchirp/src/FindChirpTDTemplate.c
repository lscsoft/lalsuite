/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpTDTemplate.c
 *
 * Author: Brown D. A., and Creighton, J. D. E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpTDTemplateCV">
Author: Brown, D. A., and Creighton, J. D. E.
$Id$
</lalVerbatim> 

<lalLaTeX>
\subsection{Module \texttt{FindChirpTDTemplate.c}}
\label{ss:FindChirpTDTemplate.c}

Provides functions to create time domain inspiral templates in a
form that can be used by the \texttt{FindChirpFilter()} function.

\subsubsection*{Prototypes}
\vspace{0.1in}
\input{FindChirpTDTemplateCP}
\idx{LALFindChirpTDTemplate()}
\idx{LALFindChirpTDNormalize()}

The function \texttt{LALFindChirpTDTemplate()} creates a time domain template
template using the inspiral package.

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

\vfill{\footnotesize\input{FindChirpTDTemplateCV}}
</lalLaTeX> 
#endif

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpTD.h>


NRCSID (FINDCHIRPTDTEMPLATEC, "$Id$");

/* <lalVerbatim file="FindChirpTDTemplateCP"> */
void
LALFindChirpTDTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params
    )
/* </lalVerbatim> */
{
  UINT4         j;
  UINT4         waveLength;
  UINT4         numPoints;
  REAL4        *xfac;
  REAL8         deltaF;
  REAL8         deltaT;
  REAL8         sampleRate;
  const REAL4   cannonDist = 1.0; /* Mpc */
  CHAR          infomsg[512];
  PPNParamStruc ppnParams;
  CoherentGW    waveform;


  INITSTATUS( status, "LALFindChirpTDTemplate", FINDCHIRPTDTEMPLATEC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */
  

  /* check that the output structures exist */
  ASSERT( fcTmplt, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( fcTmplt->data, status,
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( fcTmplt->data->data, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check that the parameter structure exists */
  ASSERT( params, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->xfacVec, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->xfacVec->data, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check we have an fft plan for the template */
  ASSERT( params->fwdPlan, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check that the timestep is positive */
  ASSERT( params->deltaT > 0, status, 
      FINDCHIRPTDH_EDELT, FINDCHIRPTDH_MSGEDELT );

  /* check that the input exists */
  ASSERT( tmplt, status, FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check that the parameter structure is set to a time domain approximant */
  switch ( params->approximant )
  {
    case TaylorT1:
      LALInfo( status, "Generating template using TaylorT1" );
      break;
      
    case TaylorT2:
      LALInfo( status, "Generating template using TaylorT2" );
      break;
      
    case TaylorT3:
      LALInfo( status, "Generating template using TaylorT3" );
      break;
      
    case GeneratePPN:
      LALInfo( status, "Generating template using GeneratePPN" );
      break;

    default:
      ABORT( status, FINDCHIRPTDH_EMAPX, FINDCHIRPTDH_MSGEMAPX );
      break;
  }
  
  /* store deltaT and zero out the time domain waveform vector */
  deltaT = params->deltaT;
  xfac = params->xfacVec->data;
  numPoints =  params->xfacVec->length;
  memset( xfac, 0, numPoints * sizeof(REAL4) );

  ASSERT( numPoints == (2 * (fcTmplt->data->length - 1)), status,
      FINDCHIRPTDH_EMISM, FINDCHIRPTDH_MSGEMISM );
  

  /* choose the time domain template */
  if ( params->approximant == GeneratePPN )
  {


    /*
     *
     * generate the waveform using LALGeneratePPNInspiral() from inject
     *
     */



    /* input parameters */
    memset( &ppnParams, 0, sizeof(PPNParamStruc) );
    ppnParams.deltaT = deltaT;
    ppnParams.mTot = tmplt->mass1 + tmplt->mass2;
    ppnParams.eta = tmplt->mass1 * tmplt->mass2 /
      ( ppnParams.mTot * ppnParams.mTot );
    ppnParams.d = 1.0;
    ppnParams.fStartIn = params->fLow;
    ppnParams.fStopIn = -1.0 / 
      (6.0 * sqrt(6.0) * LAL_PI * ppnParams.mTot * LAL_MTSUN_SI);

    /* generate waveform amplitude and phase */
    memset( &waveform, 0, sizeof(CoherentGW) );
    LALGeneratePPNInspiral( status->statusPtr, &waveform, &ppnParams );
    CHECKSTATUSPTR( status );

    /* print termination information and check sampling */
    LALInfo( status, ppnParams.termDescription );
    if ( ppnParams.dfdt > 2.0 )
    {
      ABORT( status, FINDCHIRPTDH_ESMPL, FINDCHIRPTDH_MSGESMPL );
    }
    if ( waveform.a->data->length > numPoints )
    {
      ABORT( status, FINDCHIRPTDH_ELONG, FINDCHIRPTDH_MSGELONG );
    }

    /* compute h(t) */
    for ( j = 0; j < waveform.a->data->length; ++j )
    {
      xfac[j] = 
        waveform.a->data->data[2*j] * cos( waveform.phi->data->data[j] );
    }

    /* free the memory allocated by LALGeneratePPNInspiral() */
    LALSDestroyVectorSequence( status->statusPtr, &(waveform.a->data) );
    CHECKSTATUSPTR( status );

    LALSDestroyVector( status->statusPtr, &(waveform.f->data) );
    CHECKSTATUSPTR( status );

    LALDDestroyVector( status->statusPtr, &(waveform.phi->data) );
    CHECKSTATUSPTR( status );

    LALFree( waveform.a );
    LALFree( waveform.f );
    LALFree( waveform.phi );
    
    /* waveform parameters needed for findchirp filter */
    tmplt->approximant = params->approximant;
    tmplt->tC = ppnParams.tc;
    tmplt->fFinal = ppnParams.fStop;
    
    fcTmplt->tmpltNorm = params->dynRange / ( cannonDist * 1.0e6 * LAL_PC_SI );
    fcTmplt->tmpltNorm *= fcTmplt->tmpltNorm;
  }
  else
  {


    /*
     *
     * generate the waveform by calling LALInspiralWave() from inspiral
     *
     */


    /* set up additional template parameters */
    deltaF = 1.0 / ((REAL8) numPoints * deltaT);
    sampleRate = 1.0 / deltaT;
    tmplt->ieta            = 1;
    tmplt->approximant     = params->approximant;
    tmplt->order           = twoPN;
    tmplt->massChoice      = m1Andm2;
    tmplt->tSampling       = sampleRate;
    tmplt->fLower          = params->fLow;
    tmplt->fCutoff         = sampleRate / 2.0 - deltaF;
    tmplt->signalAmplitude = 1.0;

    /* compute the tau parameters from the input template */
    LALInspiralParameterCalc( status->statusPtr, tmplt );
    CHECKSTATUSPTR( status );

    /* determine the length of the chirp in sample points */
    LALInspiralWaveLength( status->statusPtr, &waveLength, *tmplt );
    CHECKSTATUSPTR( status );

    if ( waveLength > numPoints )
    {
      ABORT( status, FINDCHIRPTDH_ELONG, FINDCHIRPTDH_MSGELONG );
    }

    /* generate the chirp in the time domain */
    LALInspiralWave( status->statusPtr, params->xfacVec, tmplt );
    CHECKSTATUSPTR( status );

    /* template dependent normalization */
    fcTmplt->tmpltNorm  = 2 * tmplt->mu;
    fcTmplt->tmpltNorm *= 2 * LAL_MRSUN_SI / ( cannonDist * 1.0e6 * LAL_PC_SI );
    fcTmplt->tmpltNorm *= params->dynRange;
    fcTmplt->tmpltNorm *= fcTmplt->tmpltNorm;
  }


  /*
   *
   * create the frequency domain findchirp template
   *
   */


  /* shift chirp to end of vector so it is the correct place for the filter */
  j = numPoints - 1;
  while ( xfac[j] == 0 )
  {
    /* search for the end of the chirp but don't fall off the array */
    if ( --j == 0 )
    {
      ABORT( status, FINDCHIRPTDH_EEMTY, FINDCHIRPTDH_MSGEEMTY );
    }
  }
  ++j;
  memmove( xfac + numPoints - j, xfac, j * sizeof( *xfac ) );
  memset( xfac, 0, ( numPoints - j ) * sizeof( *xfac ) );

  /* fft chirp */
  LALForwardRealFFT( status->statusPtr, fcTmplt->data, params->xfacVec, 
      params->fwdPlan );
  CHECKSTATUSPTR( status );

  /* copy the template parameters to the findchirp template structure */
  memcpy( &(fcTmplt->tmplt), tmplt, sizeof(InspiralTemplate) );

  /* print the template normalization constant */
  if ( lalDebugLevel & LALINFO )
  {
    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg), 
        "tmpltNorm = %e\n", fcTmplt->tmpltNorm );
    LALInfo( status, infomsg );
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpTDTemplateCP"> */
void
LALFindChirpTDNormalize(
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    FindChirpSegment           *fcSeg,
    FindChirpDataParams        *params
    )
/* </lalVerbatim> */
{
  UINT4         k;
  REAL4        *tmpltPower;
  COMPLEX8     *wtilde;
  REAL4        *segNorm;
  REAL4         segNormSum;

  INITSTATUS( status, "LALFindChirpTDNormalize", FINDCHIRPTDTEMPLATEC );

  /* check the required input exists */
  ASSERT( fcTmplt, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( fcSeg, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  ASSERT( params, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  ASSERT( params->wtildeVec, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->wtildeVec->data, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  ASSERT( params->tmpltPowerVec, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );
  ASSERT( params->tmpltPowerVec->data, status, 
      FINDCHIRPTDH_ENULL, FINDCHIRPTDH_MSGENULL );

  /* check that the parameter structure is set to a time domain approximant */
  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case GeneratePPN:
      break;
    default:
      ABORT( status, FINDCHIRPTDH_EMAPX, FINDCHIRPTDH_MSGEMAPX );
      break;
  }

  tmpltPower = params->tmpltPowerVec->data;
  wtilde     = params->wtildeVec->data;
  segNorm    = fcSeg->segNorm->data;

  memset( tmpltPower, 0, params->tmpltPowerVec->length * sizeof(REAL4) );
  memset( segNorm, 0, fcSeg->segNorm->length * sizeof(REAL4) );
  
  /* re-compute data normalization using template power */
  segNormSum = 0;
  for ( k = 1; k < fcTmplt->data->length; ++k )
  {
    REAL4 re = fcTmplt->data->data[k].re;
    REAL4 im = fcTmplt->data->data[k].im;
    REAL4 power = re * re + im * im;
    tmpltPower[k] = power * wtilde[k].re;
    segNormSum += tmpltPower[k];
    segNorm[k] = segNormSum;
  }

  /* normal exit */
  RETURN( status );
}
