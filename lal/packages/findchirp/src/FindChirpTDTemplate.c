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
  REAL8         deltaF;
  REAL8         sampleRate;
  const REAL4   cannonDist = 1.0; /* Mpc */

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
    case TaylorT2:
    case TaylorT3:
      break;
    default:
      ABORT( status, FINDCHIRPTDH_EMAPX, FINDCHIRPTDH_MSGEMAPX );
      break;
  }
  

  /*
   *
   * compute the time domain template
   *
   */

  
  /* set up additional template parameters */
  numPoints = 2 * (fcTmplt->data->length + 1);
  deltaF = 1.0 / ((REAL8) numPoints * params->deltaT);
  sampleRate = 1.0 / params->deltaT;
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

  if ( waveLength > params->xfacVec->length )
  {
    ABORT( status, FINDCHIRPTDH_ELONG, FINDCHIRPTDH_MSGELONG );
  }

  /* generate the chirp in the time domain */
  LALInspiralWave( status->statusPtr, params->xfacVec, tmplt );
  CHECKSTATUSPTR( status );

  /* shift chirp to end of vector so it is the correct place for the filter */
  j = params->xfacVec->length;
  while ( params->xfacVec->data[--j] == 0 )
  {
    /* search for the end of the chirp */
  }
  ++j;
  memmove( params->xfacVec->data + params->xfacVec->length - j, 
      params->xfacVec->data,
      j * sizeof( *params->xfacVec->data ) );
  memset( params->xfacVec->data, 0, 
      ( params->xfacVec->length - j ) * sizeof( *params->xfacVec->data ) );

  /* fft chirp */
  LALForwardRealFFT( status->statusPtr, fcTmplt->data, params->xfacVec, 
      params->fwdPlan );
  CHECKSTATUSPTR( status );

  /* template dependent normalization */
  fcTmplt->tmpltNorm  = 2 * tmplt->mu;
  fcTmplt->tmpltNorm *= 2 * LAL_MRSUN_SI / ( cannonDist * 1.0e6 * LAL_PC_SI );
  fcTmplt->tmpltNorm *= params->dynRange;
  fcTmplt->tmpltNorm *= fcTmplt->tmpltNorm;

  /* copy the template parameters to the finchirp template structure */
  memcpy( &(fcTmplt->tmplt), tmplt, sizeof(InspiralTemplate) );

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
