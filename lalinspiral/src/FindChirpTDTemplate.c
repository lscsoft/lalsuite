/*
*  Copyright (C) 2007 Stas Babak, Duncan Brown, Jolien Creighton, Craig Robinson
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
 * File Name: FindChirpTDTemplate.c
 *
 * Author: Brown D. A., and Creighton, J. D. E.
 *
 */

/**
\author Brown, D. A., and Creighton, J. D. E.
\file
\ingroup FindChirpTD_h

\brief Provides functions to create time domain inspiral templates in a
form that can be used by the <tt>FindChirpFilter()</tt> function.

\heading{Prototypes}

The function <tt>LALFindChirpTDTemplate()</tt> creates a time domain template
template using the inspiral package.

\heading{Algorithm}

Blah.

\heading{Uses}
\code
LALCalloc()
LALFree()
LALCreateVector()
LALDestroyVector()
\endcode

\heading{Notes}

*/

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpTD.h>

void
LALFindChirpTDTemplate (
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    InspiralTemplate           *tmplt,
    FindChirpTmpltParams       *params
    )

{
  UINT4         j;
  UINT4         shift;
  UINT4         waveLength;
  UINT4         numPoints;
  REAL4        *xfac;
  REAL8         deltaF;
  REAL8         deltaT;
  REAL8         sampleRate;
  const REAL4   cannonDist = 1.0; /* Mpc */
  /*CHAR          infomsg[512];*/
  PPNParamStruc ppnParams;
  CoherentGW    waveform;

  REAL4Vector  *tmpxfac = NULL; /* Used for band-passing */


  INITSTATUS(status);
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
    case TaylorT4:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case EOBNR:
    case FindChirpPTF:
    case EOBNRv2:
    case IMRPhenomB:
    case IMRPhenomC:
      break;

    default:
      ABORT( status, FINDCHIRPTDH_EMAPX, FINDCHIRPTDH_MSGEMAPX );
      break;
  }

  /* store deltaT and zero out the time domain waveform vector */
  deltaT = params->deltaT;
  sampleRate = 1.0 / deltaT;
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
    tmplt->ieta            = 1;
    tmplt->approximant     = params->approximant;
    tmplt->order           = params->order;
    tmplt->massChoice      = m1Andm2;
    tmplt->tSampling       = sampleRate;
    tmplt->fLower          = params->fLow;
    tmplt->fCutoff         = sampleRate / 2.0 - deltaF;
    /* get the template norm right */
    if ( params->approximant==EOBNR )
    {
     /* lalinspiral EOBNR code produces correct norm when 
        fed unit signalAmplitude and non-physical distance */
      tmplt->signalAmplitude = 1.0;
      tmplt->distance      = -1.0;
    }
    else if ( params->approximant==EOBNRv2 )
    {
     /* this formula sets the ampl0 variable to 1.0 
      * within the lalsimulation EOBNRv2 waveform engine 
      * which again produces a correct template norm     */
      tmplt->distance      = tmplt->totalMass*LAL_MRSUN_SI;
    }
    else if ( (params->approximant==IMRPhenomB) || (params->approximant==IMRPhenomC) )
    {
      /* 1Mpc standard distance - not clear if this produces correct norm */
      tmplt->distance      = 1.0;
      tmplt->spin1[2]      = 2 * tmplt->chi/(1. + sqrt(1.-4.*tmplt->eta));
    }

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


  /* Taper the waveform if required */
  if ( params->taperTmplt != LAL_SIM_INSPIRAL_TAPER_NONE )
  {
    if ( XLALSimInspiralREAL4WaveTaper( params->xfacVec, params->taperTmplt )
           == XLAL_FAILURE )
    {
      ABORTXLAL( status );
    }
  }

  /* Find the end of the chirp */
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

  /* Band pass the template if required */
  if ( params->bandPassTmplt )
  {
    REAL4Vector bpVector; /*Used to save time */

    /* We want to shift the template to the middle of the vector so */
    /* that band-passing will work properly */
    shift = ( numPoints - j ) / 2;
    memmove( xfac + shift, xfac, j * sizeof( *xfac ) );
    memset( xfac, 0, shift * sizeof( *xfac ) );
    memset( xfac + ( numPoints + j ) / 2, 0,
         ( numPoints - ( numPoints + j ) / 2 ) * sizeof( *xfac ) );


    /* Select an appropriate part of the vector to band pass. */
    /* band passing the whole thing takes a lot of time */
    if ( j > 2 * sampleRate && 2 * j <= numPoints )
    {
      bpVector.length = 2 * j;
      bpVector.data   = params->xfacVec->data + numPoints / 2 - j;
    }
    else if ( j <= 2 * sampleRate && j + 2 * sampleRate <= numPoints )
    {
      bpVector.length = j + 2 * sampleRate;
      bpVector.data   = params->xfacVec->data
                   + ( numPoints - j ) / 2 - (INT4)sampleRate;
    }
    else
    {
      bpVector.length = params->xfacVec->length;
      bpVector.data   = params->xfacVec->data;
    }

    if ( XLALBandPassInspiralTemplate( &bpVector, 0.95 * tmplt->fLower,
                 1.02 * tmplt->fFinal, sampleRate ) == XLAL_FAILURE )
    {
      ABORTXLAL( status );
    }

    /* Now we need to do the shift to the end. */
    /* Use a temporary vector to avoid mishaps */
    if ( ( tmpxfac = XLALCreateREAL4Vector( numPoints ) ) == NULL )
    {
      ABORTXLAL( status );
    }

    if ( params->approximant == EOBNR || params->approximant == EOBNRv2
        || params->approximant == IMRPhenomB || params->approximant == IMRPhenomC )
    {
      /* We need to do something slightly different for EOBNR */
      UINT4 endIndx = (UINT4) (tmplt->tC * sampleRate);

      memcpy( tmpxfac->data, xfac + ( numPoints - j ) / 2 + endIndx,
          ( numPoints - ( numPoints - j ) / 2 - endIndx ) * sizeof( *xfac ) );

      memcpy( tmpxfac->data + numPoints - ( numPoints - j ) / 2 - endIndx,
                  xfac, ( ( numPoints - j ) / 2 + endIndx ) * sizeof( *xfac ) );
    }
    else
    {
      memcpy( tmpxfac->data, xfac + ( numPoints + j ) / 2,
          ( numPoints - ( numPoints + j ) / 2 ) * sizeof( *xfac ) );

      memcpy( tmpxfac->data + numPoints - ( numPoints + j ) / 2,
                    xfac, ( numPoints + j ) / 2 * sizeof( *xfac ) );
    }

    memcpy( xfac, tmpxfac->data, numPoints * sizeof( *xfac ) );

    XLALDestroyREAL4Vector( tmpxfac );
    tmpxfac = NULL;
  }
  else if ( params->approximant == EOBNR || params->approximant == EOBNRv2
      || params->approximant == IMRPhenomB || params->approximant == IMRPhenomC )
  {
    /* For EOBNR we shift so that tC is at the end of the vector */
    if ( ( tmpxfac = XLALCreateREAL4Vector( numPoints ) ) == NULL )
    {
      ABORTXLAL( status );
    }

    /* Set the coalescence index depending on tC */
    j = (UINT4) (tmplt->tC * sampleRate);
    memcpy( tmpxfac->data + numPoints - j, xfac, j * sizeof( *xfac ) );
    memcpy( tmpxfac->data, xfac + j, ( numPoints - j ) * sizeof( *xfac ) );
    memcpy( xfac, tmpxfac->data, numPoints * sizeof( *xfac ) );
    XLALDestroyREAL4Vector( tmpxfac );
    tmpxfac = NULL;
  }
  else
  {
    /* No need for so much shifting around if not band passing */
    /* shift chirp to end of vector so it is the correct place for the filter */
    memmove( xfac + numPoints - j, xfac, j * sizeof( *xfac ) );
    memset( xfac, 0, ( numPoints - j ) * sizeof( *xfac ) );
  }

  /*
   *
   * create the frequency domain findchirp template
   *
   */

  /* fft chirp */
  if ( XLALREAL4ForwardFFT( fcTmplt->data, params->xfacVec,
      params->fwdPlan ) == XLAL_FAILURE )
  {
    ABORTXLAL( status );
  }

  /* copy the template parameters to the findchirp template structure */
  memcpy( &(fcTmplt->tmplt), tmplt, sizeof(InspiralTemplate) );

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



void
LALFindChirpTDNormalize(
    LALStatus                  *status,
    FindChirpTemplate          *fcTmplt,
    FindChirpSegment           *fcSeg,
    FindChirpDataParams        *params
    )

{
  UINT4         k;
  REAL4        *tmpltPower;
  COMPLEX8     *wtilde;
  REAL4        *segNorm;
  REAL4         segNormSum;

  INITSTATUS(status);

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
    case TaylorT4:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case EOBNR:
    case FindChirpPTF:
    case EOBNRv2:
    case IMRPhenomB:
    case IMRPhenomC:
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
    REAL4 re = crealf(fcTmplt->data->data[k]);
    REAL4 im = cimagf(fcTmplt->data->data[k]);
    REAL4 power = re * re + im * im;
    tmpltPower[k] = power * crealf(wtilde[k]);
    segNormSum += tmpltPower[k];
    segNorm[k] = segNormSum;
  }

  /* normal exit */
  RETURN( status );
}
