/*
*  Copyright (C) 2007 Diego Fazi, Duncan Brown, Craig Robinson
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
 * File Name: FindChirpTemplate.c
 *
 * Author: Brown D. A.
 *
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown, D. A.
 * \file
 * \ingroup FindChirp_h
 *
 * \brief Provides functions to initialize template creation routines.
 *
 * \heading{Prototypes}
 *
 * The function <tt>LALFindChirpTemplateInit()</tt> takes as input the address
 * of a structure of type \c FindChirpInitParams containing the correct
 * values to intialize a search. It creates a structure of type
 * \c FindChirpTmpltParams as described above and returns its address.
 *
 * The function <tt>LALFindChirpTemplateFinalize()</tt> takes as the address
 * of a structure of type \c FindChirpTmpltParams destroys this
 * structure and sets the address to NULL.
 *
 * \heading{Algorithm}
 *
 * Blah.
 *
 * \heading{Uses}
 * \code
 * LALCalloc()
 * LALFree()
 * LALCreateVector()
 * LALDestroyVector()
 * \endcode
 *
 * \heading{Notes}
 *
 */

#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpACTD.h>

void
LALFindChirpTemplateInit (
    LALStatus                  *status,
    FindChirpTmpltParams      **output,
    FindChirpInitParams        *params
    )

{
  UINT4                         k;
  FindChirpTmpltParams         *outputPtr;
  REAL4                        *xfac = NULL;
  const REAL4                   exponent = -1.0/3.0;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( output, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*output, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the number of points is positive */
  ASSERT( params->numPoints > 0, status,
      FINDCHIRPH_ENUMZ, FINDCHIRPH_MSGENUMZ );


  /*
   *
   * create tmplt generation parameters structure
   *
   */


  /* create the output structure */
  outputPtr = *output = (FindChirpTmpltParams *)
    LALCalloc( 1, sizeof(FindChirpTmpltParams) );
  if ( ! outputPtr )
  {
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* store the waveform approximant */
  outputPtr->approximant = params->approximant;

  switch ( params->approximant )
  {
    case FindChirpPTF:
      /* create workspace for the dynamical variables needed to */
      /* compute the PTF Q(t) vectors                           */
      outputPtr->PTFphi = XLALCreateVector( params->numPoints );
      if ( ! outputPtr->PTFphi )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }
      outputPtr->PTFomega_2_3 = XLALCreateVector( params->numPoints );
      if ( ! outputPtr->PTFomega_2_3 )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }
      outputPtr->PTFe1 = XLALCreateVectorSequence( 3, params->numPoints );
      if ( ! outputPtr->PTFe1 )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }
      outputPtr->PTFe2 = XLALCreateVectorSequence( 3, params->numPoints );
      if ( ! outputPtr->PTFe2 )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }

    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case IMRPhenomB:
      /* time domain waveforms use xfac to store the time domain waveform */
      LALCreateVector( status->statusPtr, &(outputPtr->xfacVec),
          params->numPoints );
      BEGINFAIL( status )
      {
        LALFree( outputPtr );
        *output = NULL;
      }
      ENDFAIL( status );

      memset( outputPtr->xfacVec->data, 0,
          outputPtr->xfacVec->length * sizeof(REAL4) );

      /* create an fft plan for the time domain waveform */
      LALCreateForwardRealFFTPlan( status->statusPtr, &(outputPtr->fwdPlan),
          params->numPoints, 1 );
      BEGINFAIL( status )
      {
        TRY( LALDestroyVector( status->statusPtr, &(outputPtr->xfacVec) ),
            status );

        LALFree( outputPtr );
        *output = NULL;
      }
      ENDFAIL( status );
      break;

    case FindChirpSP:
    case BCV:
    case BCVSpin:
      /* freq domain waveforms need xfac vector containing k^(-1/3) */
      LALCreateVector( status->statusPtr, &(outputPtr->xfacVec),
          params->numPoints/2 + 1 );
      BEGINFAIL( status )
      {
        LALFree( outputPtr );
        *output = NULL;
      }
      ENDFAIL( status );

      xfac = outputPtr->xfacVec->data;

      xfac[0] = 0;
      for (k = 1; k < outputPtr->xfacVec->length; ++k)
        xfac[k] = pow( (REAL4) k, exponent );
      break;

    case AmpCorPPN:
      /* create workspace memory for the time-domain Q vectors */
      outputPtr->ACTDVecs =
                      XLALCreateVectorSequence( NACTDVECS, params->numPoints );
      if ( ! outputPtr->ACTDVecs )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }

      /* create a forward FFT plan */
      outputPtr->fwdPlan =
        XLALCreateForwardREAL4FFTPlan( params->numPoints, 1 );
      if ( ! outputPtr->fwdPlan )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }

      break;

    default:
      /* unknown approximant type */
      LALFree( outputPtr );
      *output = NULL;
      ABORT( status, FINDCHIRPH_EUAPX, FINDCHIRPH_MSGEUAPX );
      break;
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN (status);
}




void
LALFindChirpTemplateFinalize (
    LALStatus                  *status,
    FindChirpTmpltParams      **output
    )

{
  FindChirpTmpltParams         *outputPtr;

  INITSTATUS(status);
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure handle is non-null and points to a non-null pointer */
  ASSERT( output, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( *output, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );


  /*
   *
   * destroy tmplt generation parameters structure
   *
   */


  /* local pointer to output */
  outputPtr = *output;

  /* destroy the fft plan if it exists */
  if ( outputPtr->fwdPlan )
  {
    LALDestroyRealFFTPlan( status->statusPtr, &(outputPtr->fwdPlan) );
    CHECKSTATUSPTR( status );
  }

  /* destroy the vector used to store part of the template if it exists */
  if ( outputPtr->xfacVec )
  {
  LALDestroyVector( status->statusPtr, &(outputPtr->xfacVec) );
  CHECKSTATUSPTR( status );
  }

  /* destroy the vectors used for the PTF template if they exist */
  if ( outputPtr->PTFphi )
  {
    XLALDestroyVector( outputPtr->PTFphi );
  }
  if ( outputPtr->PTFomega_2_3 )
  {
    XLALDestroyVector( outputPtr->PTFomega_2_3 );
  }
  if ( outputPtr->PTFe1 )
  {
    XLALDestroyVectorSequence( outputPtr->PTFe1 );
  }
  if ( outputPtr->PTFe2 )
  {
    XLALDestroyVectorSequence( outputPtr->PTFe2 );
  }

  /* destroy ACTD vectors if they exist */
  if ( outputPtr->ACTDVecs )
  {
    XLALDestroyVectorSequence( outputPtr->ACTDVecs );
  }


  /* free the structure */
  LALFree( outputPtr );
  *output = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
