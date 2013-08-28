/*
*  Copyright (C) 2007 Sukanta Bose, Diego Fazi, Duncan Brown, Gareth Jones, Craig Robinson , Sean Seader
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
 * File Name: FindChirpFilterInit.c
 *
 * Author: Brown, D. A., BCV-Modifications by Messaritaki E.
 *
 *
 *-----------------------------------------------------------------------
 */

/**
 * \author Brown D. A., BCV-Modifications by Messaritaki E.
 * \file
 *
 */

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/SeqFactories.h>
#include <lal/FindChirp.h>
#include <lal/LALDatatypes.h>
#include <lal/FindChirpACTD.h>

void
LALCreateFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output,
    FindChirpInitParams        *params
    )

{
  FindChirpFilterInput         *outputPtr;

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

  /* check that the approximants is of a known type */
  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorF2:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case FindChirpSP:
    case FindChirpPTF:
    case BCV:
    case BCVSpin:
    case AmpCorPPN:
    case IMRPhenomB:
      break;
    default:
      ABORT( status, FINDCHIRPH_EUAPX, FINDCHIRPH_MSGEUAPX );
      break;
  }


  /*
   *
   * create the findchirp filter input structure
   *
   */


  /* create the output structure */
  outputPtr = *output = (FindChirpFilterInput *)
    LALCalloc( 1, sizeof(FindChirpFilterInput) );
  if ( ! outputPtr )
  {
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* create memory for the chirp template structure */
  outputPtr->fcTmplt = (FindChirpTemplate *)
    LALCalloc( 1, sizeof(FindChirpTemplate) );
  if ( !outputPtr->fcTmplt )
  {
    LALFree( *output );
    *output = NULL;
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  if( params->approximant == AmpCorPPN )
  {
    outputPtr->fcTmplt->ACTDtilde = 
     XLALCreateCOMPLEX8VectorSequence( NACTDVECS, params->numPoints / 2 + 1 );
    if ( ! outputPtr->fcTmplt->ACTDtilde )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }
  }
  else 
  {  
    if ( params->approximant == FindChirpPTF )
    {
      /* create memory for the PTF template data */

      outputPtr->fcTmplt->PTFQ = XLALCreateVectorSequence( 5, params->numPoints );
      if ( ! outputPtr->fcTmplt->PTFQ )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }

      outputPtr->fcTmplt->PTFQtilde = 
        XLALCreateCOMPLEX8VectorSequence( 5, params->numPoints / 2 + 1 );
      if ( ! outputPtr->fcTmplt->PTFQtilde )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }

      outputPtr->fcTmplt->PTFBinverse = XLALCreateArrayL( 2, 5, 5 );
      if ( ! outputPtr->fcTmplt->PTFBinverse )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }

      outputPtr->fcTmplt->PTFB = XLALCreateArrayL( 2, 5, 5 );
      if ( ! outputPtr->fcTmplt->PTFB )
      {
        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
      }
    }
    /* create memory for the chirp template data */
    LALCCreateVector (status->statusPtr, &(outputPtr->fcTmplt->data),
        params->numPoints / 2 + 1 );
    BEGINFAIL( status )
    {
      LALFree( outputPtr->fcTmplt );
      outputPtr->fcTmplt = NULL;
      LALFree( *output );
      *output = NULL;
    }
    ENDFAIL( status );
  }

  /* create memory for BCVSpin orthonormalised amplitude vectors */
  if ( params->approximant == BCVSpin )
  {
    LALDCreateVector (status->statusPtr, &(outputPtr->fcTmplt->A1BCVSpin),
        ((params->numPoints)/2)+1 );
    BEGINFAIL( status )
    {
      LALFree( outputPtr->fcTmplt );
      outputPtr->fcTmplt = NULL;
      LALFree( *output );
      *output = NULL;
    }
    ENDFAIL( status );

    LALDCreateVector (status->statusPtr, &(outputPtr->fcTmplt->A2BCVSpin),
        ((params->numPoints)/2)+1 );
    BEGINFAIL( status )
    {
      LALFree( outputPtr->fcTmplt );
      outputPtr->fcTmplt = NULL;
      LALFree( *output );
      *output = NULL;
    }
    ENDFAIL( status );

    LALDCreateVector (status->statusPtr, &(outputPtr->fcTmplt->A3BCVSpin),
        ((params->numPoints)/2)+1 );
    BEGINFAIL( status )
    {
      LALFree( outputPtr->fcTmplt );
      outputPtr->fcTmplt = NULL;
      LALFree( *output );
      *output = NULL;
    }
    ENDFAIL( status );
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}




void
LALDestroyFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output
    )

{
  FindChirpFilterInput         *outputPtr;

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
   * destroy the findchirp input structure
   *
   */

  /* local pointer to output */
  outputPtr = *output;

  /* destroy the chirp template data storage */
  if ( outputPtr->fcTmplt->data )
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->fcTmplt->data) );
    CHECKSTATUSPTR( status );
  }

  /* destroy BCVSpin amplitude vectors */
  if (outputPtr->fcTmplt->A1BCVSpin)
  {
    LALDDestroyVector( status->statusPtr, &(outputPtr->fcTmplt->A1BCVSpin) );
    CHECKSTATUSPTR( status );
  }

  if (outputPtr->fcTmplt->A2BCVSpin)
  {
    LALDDestroyVector( status->statusPtr, &(outputPtr->fcTmplt->A2BCVSpin) );
    CHECKSTATUSPTR( status );
  }

  if (outputPtr->fcTmplt->A3BCVSpin)
  {
    LALDDestroyVector( status->statusPtr, &(outputPtr->fcTmplt->A3BCVSpin) );
    CHECKSTATUSPTR( status );
  }

  /* destroy the PTF work space */
  if ( outputPtr->fcTmplt->PTFQ )
  {
    XLALDestroyVectorSequence( outputPtr->fcTmplt->PTFQ );
  }
  if ( outputPtr->fcTmplt->PTFQtilde )
  {
    XLALDestroyCOMPLEX8VectorSequence( outputPtr->fcTmplt->PTFQtilde );
  }
  if ( outputPtr->fcTmplt->PTFBinverse )
  {
    XLALDestroyArray( outputPtr->fcTmplt->PTFBinverse );
  }
  if ( outputPtr->fcTmplt->PTFB )
  {
    XLALDestroyArray( outputPtr->fcTmplt->PTFB );
  }

  /* destroy the ACTDtilde vector sequence */
  if ( outputPtr->fcTmplt->ACTDtilde )
  {
    XLALDestroyCOMPLEX8VectorSequence( outputPtr->fcTmplt->ACTDtilde );
  }



  /* destroy the chirp template structure */
  LALFree( outputPtr->fcTmplt );

  /* destroy the filter input structure */
  LALFree( outputPtr );
  *output = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



void
LALFindChirpFilterInit (
    LALStatus                  *status,
    FindChirpFilterParams     **output,
    FindChirpInitParams        *params
    )

{
  FindChirpFilterParams        *outputPtr;

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

  /* make sure that the number of points in a segment is positive */
  ASSERT( params->numPoints > 0,  status,
      FINDCHIRPH_ENUMZ, FINDCHIRPH_MSGENUMZ );

  /* check that the user has given a known approximant */
  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorF2:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case EOBNR:
    case EOBNRv2:
    case FindChirpSP:
    case FindChirpPTF:
    case BCV:
    case BCVSpin:
    case IMRPhenomB:
    case AmpCorPPN:
      break;
    default:
      ABORT( status, FINDCHIRPH_EUAPX, FINDCHIRPH_MSGEUAPX );
      break;
  }


  /*
   *
   * allocate memory for the FindChirpFilterParams
   *
   */


  /* create the output structure */
  outputPtr = *output = (FindChirpFilterParams *)
    LALCalloc( 1, sizeof(FindChirpFilterParams) );
  if ( ! outputPtr )
  {
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* create memory for the chisq parameters */
  outputPtr->chisqParams = (FindChirpChisqParams *)
    LALCalloc( 1, sizeof(FindChirpChisqParams) );
  if ( !outputPtr->chisqParams )
  {
    LALFree( outputPtr );
    *output = NULL;
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* store the filter approximant in the filter and chisq params */
  outputPtr->approximant = outputPtr->chisqParams->approximant =
    params->approximant;

  outputPtr->order = params->order;

  /* create memory for the chisq input */
  outputPtr->chisqInput = (FindChirpChisqInput *)
    LALCalloc( 1, sizeof(FindChirpChisqInput) );
  if ( !outputPtr->chisqInput )
  {
    LALFree( outputPtr->chisqParams );
    LALFree( outputPtr );
    *output = NULL;
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* create memory for the additional BCV chisq input */
  if ( params->approximant == BCV )
  {
    outputPtr->chisqInputBCV = (FindChirpChisqInput *)
      LALCalloc( 1, sizeof(FindChirpChisqInput) );
    if ( !outputPtr->chisqInputBCV )
    {
      LALFree( outputPtr->chisqInput );
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }
  }


  /*
   *
   * create fft plans and workspace vectors
   *
   */


  /* create plan for optimal filter */
  LALCreateReverseComplexFFTPlan( status->statusPtr,
      &(outputPtr->invPlan), params->numPoints, 1 );
  BEGINFAIL( status )
  {
    LALFree( outputPtr->chisqInput );
    if ( outputPtr->chisqInputBCV )
    {
      LALFree( outputPtr->chisqInputBCV );
    }
    LALFree( outputPtr->chisqParams );
    LALFree( outputPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  if ( params->approximant == FindChirpPTF )
  {
    /* create workspace vector for PTF filter: time domain*/
    outputPtr->PTFqVec =
      XLALCreateCOMPLEX8VectorSequence ( 5, params->numPoints );
    if ( ! outputPtr->PTFqVec )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }

    outputPtr->PTFA = XLALCreateArrayL( 2, 5, 5 );
    if ( ! outputPtr->PTFA )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }

    outputPtr->PTFMatrix = XLALCreateArrayL( 2, 5, 5 );
    if ( ! outputPtr->PTFMatrix )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }

    outputPtr->PTFsnrVec = XLALCreateCOMPLEX8Vector( params->numPoints );
    if ( ! outputPtr->PTFsnrVec )
    {
      ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
    }
  }
  else if ( params->approximant == AmpCorPPN )
  {
    /* workspace vector for optimal AmpCorPPN filter: time-domain */
    INT4 i;

    outputPtr->qVecACTD =
                       LALCalloc( NACTDVECS, sizeof( *outputPtr->qVecACTD ) );

    for( i=0; i < NACTDVECS; ++i )
    {
      outputPtr->qVecACTD[i] = XLALCreateCOMPLEX8Vector( params->numPoints );
      BEGINFAIL( status )
      {
        INT4 j;

        for( j = 0; j < i; ++j )
        {
          LALCDestroyVector( status->statusPtr, &( outputPtr->qVecACTD[i] ) );
          CHECKSTATUSPTR( status );
        }
        LALFree( outputPtr->qVecACTD );
        LALFree( outputPtr->chisqInput );
        LALFree( outputPtr->chisqParams );
        LALFree( outputPtr );
        *output = NULL;
      }
      ENDFAIL( status );
    }
  }
  else
  {
    /* create workspace vector for optimal filter: time domain */
    LALCCreateVector( status->statusPtr, &(outputPtr->qVec),
        params->numPoints );
    BEGINFAIL( status )
    {
      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );

      LALFree( outputPtr->chisqInput );
      if ( outputPtr->chisqInputBCV )
      {
        LALFree( outputPtr->chisqInputBCV );
      }
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }

  /* create additional workspace vector for optimal BCV filter: time domain */
  if ( params->approximant == BCV )
  {
    LALCCreateVector( status->statusPtr, &(outputPtr->qVecBCV),
        params->numPoints );
    BEGINFAIL( status )
    {
      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
          status );

      LALFree( outputPtr->chisqInput );
      LALFree( outputPtr->chisqInputBCV );
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }
  else if ( params->approximant == BCVSpin )
  {
    /* workspace vector for optimal BCVSpin filter: time domain */
    LALCCreateVector( status->statusPtr, &(outputPtr->qVecBCVSpin1),
        params->numPoints );
    BEGINFAIL( status )
    {
      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
          status );

      LALFree( outputPtr->chisqInput );
      if ( outputPtr->chisqInputBCV )
      {
        LALFree( outputPtr->chisqInputBCV );
      }
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );

    LALCCreateVector( status->statusPtr, &(outputPtr->qVecBCVSpin2),
        params->numPoints );
    BEGINFAIL( status )
    {
      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
          status );
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin1) ),
          status );

      LALFree( outputPtr->chisqInput );
      if ( outputPtr->chisqInputBCV )
      {
        LALFree( outputPtr->chisqInputBCV );
      }
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }


  if ( params->approximant == AmpCorPPN )
  {
    /* workspace vector for optimal AmpCorPPN filter: freq-domain */
    INT4 i;

    outputPtr->qtildeVecACTD =
                   LALCalloc( NACTDVECS, sizeof( *outputPtr->qtildeVecACTD ) );

    for( i=0; i < NACTDVECS; ++i )
    {
      outputPtr->qtildeVecACTD[i] =
                                 XLALCreateCOMPLEX8Vector( params->numPoints );
      BEGINFAIL( status )
      {
        INT4 j;

        for( j = 0; j < i; ++j )
        {
          LALCDestroyVector( status->statusPtr,
													                 &( outputPtr->qtildeVecACTD[i] ) );
          CHECKSTATUSPTR( status );
        }
        LALFree( outputPtr->qtildeVecACTD );
        LALFree( outputPtr->chisqInput );
        LALFree( outputPtr->chisqParams );
        LALFree( outputPtr );
        *output = NULL;
      }
      ENDFAIL( status );
    }

  }
  else
  {
    /* create workspace vector for optimal filter: freq domain */
    LALCCreateVector( status->statusPtr, &(outputPtr->qtildeVec),
        params->numPoints );
    BEGINFAIL( status )
    {
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
          status );
      if ( outputPtr->qVecBCV )
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCV) ),
            status );
      }
      if ( outputPtr->qVecBCVSpin1 )
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin1) ),
            status );
      }
      if ( outputPtr->qVecBCVSpin2 )
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin2) ),
            status );
      }

      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );

      LALFree( outputPtr->chisqInput );
      if ( outputPtr->chisqInputBCV )
      {
        LALFree( outputPtr->chisqInputBCV );
      }
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }

  /* create workspace vector for optimal filter: freq domain */
  if ( params->approximant == BCV )
  {
    LALCCreateVector( status->statusPtr, &(outputPtr->qtildeVecBCV),
        params->numPoints );
    BEGINFAIL( status )
    {
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
          status );
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCV) ),
          status );
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ),
          status );

      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );

      LALFree( outputPtr->chisqInput );
      LALFree( outputPtr->chisqInputBCV );
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }
  else if ( params->approximant == BCVSpin )
  {
    /* create workspace vector for optimal filter: freq domain */
    LALCCreateVector( status->statusPtr, &(outputPtr->qtildeVecBCVSpin1),
        params->numPoints );
    BEGINFAIL( status )
    {
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
          status );
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ),
          status );
      TRY( LALCDestroyVector(status->statusPtr,&(outputPtr->qVecBCVSpin1)),
          status );
      TRY( LALCDestroyVector(status->statusPtr,&(outputPtr->qVecBCVSpin2)),
          status );

      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );

      LALFree( outputPtr->chisqInput );
      if ( outputPtr->chisqInputBCV )
      {
        LALFree( outputPtr->chisqInputBCV );
      }
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );

    LALCCreateVector( status->statusPtr, &(outputPtr->qtildeVecBCVSpin2),
        params->numPoints );
    BEGINFAIL( status )
    {
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
          status );
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ),
          status );
      TRY( LALCDestroyVector(status->statusPtr,&(outputPtr->qVecBCVSpin1)),
          status );
      TRY( LALCDestroyVector(status->statusPtr,&(outputPtr->qVecBCVSpin2)),
          status );
      TRY( LALCDestroyVector(status->statusPtr,&(outputPtr->qtildeVecBCVSpin1)),
          status );

      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );

      LALFree( outputPtr->chisqInput );
      if ( outputPtr->chisqInputBCV )
      {
        LALFree( outputPtr->chisqInputBCV );
      }
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }

  if ( params->numChisqBins )
  {
    /* create workspace vector for chisq filter */
    LALCreateVector (status->statusPtr, &(outputPtr->chisqVec),
        params->numPoints);
    BEGINFAIL( status )
    {
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ),
          status );
      if ( outputPtr->qtildeVecBCV)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVecBCV) ),
            status );
      }
      if ( outputPtr->qtildeVecBCVSpin1)
      {
        TRY( LALCDestroyVector( status->statusPtr,
              &(outputPtr->qtildeVecBCVSpin1) ),
            status );
      }
      if ( outputPtr->qtildeVecBCVSpin2)
      {
        TRY( LALCDestroyVector( status->statusPtr,
              &(outputPtr->qtildeVecBCVSpin2) ),
            status );
      }
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
          status );
      if ( outputPtr->qVecBCV)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCV) ),
            status );
      }
      if ( outputPtr->qVecBCVSpin1)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin1) ),
            status );
      }
      if ( outputPtr->qVecBCVSpin2)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin2) ),
            status );
      }

      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );

      LALFree( outputPtr->chisqInput );
      if ( outputPtr->chisqInputBCV )
      {
        LALFree( outputPtr->chisqInputBCV );
      }
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }


  /*
   *
   * create vector to store snrsq, if required
   *
   */


  if ( params->createRhosqVec )
  {
    outputPtr->rhosqVec = (REAL4TimeSeries *)
      LALCalloc( 1, sizeof(REAL4TimeSeries) );
    LALCreateVector (status->statusPtr, &(outputPtr->rhosqVec->data),
        params->numPoints);
    BEGINFAIL( status )
    {
      if ( outputPtr->chisqVec )
      {
        TRY( LALDestroyVector( status->statusPtr, &(outputPtr->chisqVec) ),
            status );
      }
      if ( outputPtr->qtildeVec )
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ),
            status );
      }
      if ( outputPtr->qtildeVecBCV )
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVecBCV)),
            status );
      }
      if ( outputPtr->qtildeVecBCVSpin1 )
      {
        TRY( LALCDestroyVector( status->statusPtr,
              &(outputPtr->qtildeVecBCVSpin1)), status );
      }
      if ( outputPtr->qtildeVecBCVSpin2 )
      {
        TRY( LALCDestroyVector( status->statusPtr,
              &(outputPtr->qtildeVecBCVSpin2)), status );
      }
      if ( outputPtr->qVec )
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
            status );
      }
      if ( outputPtr->qVecBCV)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCV) ),
            status );
      }
      if ( outputPtr->qVecBCVSpin1)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin1) ),
            status );
      }
      if ( outputPtr->qVecBCVSpin2)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin2) ),
            status );
      }

      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );

      LALFree( outputPtr->chisqInput );
      if ( outputPtr->chisqInputBCV )
      {
        LALFree( outputPtr->chisqInputBCV );
      }
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }


  /*
   *
   * create vector to store c data, if required
   *
   */


  if ( params->createCVec )
  {
    outputPtr->cVec = (COMPLEX8TimeSeries *)
      LALCalloc( 1, sizeof(COMPLEX8TimeSeries) );
    LALCCreateVector (status->statusPtr, &(outputPtr->cVec->data),
        params->numPoints);
    BEGINFAIL( status )
    {
      if ( outputPtr->chisqVec )
      {
        TRY( LALDestroyVector( status->statusPtr, &(outputPtr->chisqVec) ),
            status );
      }
      if ( outputPtr->qtildeVec )
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ),
            status );
      }
      if ( outputPtr->rhosqVec )
      {
        TRY( LALDestroyVector( status->statusPtr, &(outputPtr->rhosqVec->data)),
            status );
      }
      if ( outputPtr->qtildeVecBCV )
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVecBCV)),
            status );
      }
      if ( outputPtr->qtildeVecBCVSpin1 )
      {
        TRY( LALCDestroyVector( status->statusPtr,
              &(outputPtr->qtildeVecBCVSpin1)), status );
      }
      if ( outputPtr->qtildeVecBCVSpin2 )
      {
        TRY( LALCDestroyVector( status->statusPtr,
              &(outputPtr->qtildeVecBCVSpin2)), status );
      }
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ),
          status );
      if ( outputPtr->qVecBCV)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCV) ),
            status );
      }
      if ( outputPtr->qVecBCVSpin1)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin1) ),
            status );
      }
      if ( outputPtr->qVecBCVSpin2)
      {
        TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin2) ),
            status );
      }

      TRY( LALDestroyComplexFFTPlan( status->statusPtr,
            &(outputPtr->invPlan) ), status );

      LALFree( outputPtr->chisqInput );
      if ( outputPtr->rhosqVec )
      {
        LALFree( outputPtr->rhosqVec );
      }
      if ( outputPtr->chisqInputBCV )
      {
        LALFree( outputPtr->chisqInputBCV );
      }
      LALFree( outputPtr->chisqParams );
      LALFree( outputPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}




void
LALFindChirpFilterFinalize (
    LALStatus                  *status,
    FindChirpFilterParams     **output
    )

{
  FindChirpFilterParams        *outputPtr;
  /*UINT4 i,j;*/

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
   * destroy the filter parameter structure
   *
   */


  /* local pointer to output structure */
  outputPtr = *output;


  /*
   *
   * destroy fft plans and workspace vectors
   *
   */

  /* destroy plan for optimal filter */
  LALDestroyComplexFFTPlan( status->statusPtr, &(outputPtr->invPlan) );
  CHECKSTATUSPTR( status );

  /* destroy workspace vector for optimal filter: freq domain */
  if ( outputPtr->qVec )
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) );
    CHECKSTATUSPTR( status );
  }


  /* destroy the PTF memory */
  if ( outputPtr->PTFqVec )
  {
    XLALDestroyCOMPLEX8VectorSequence( outputPtr->PTFqVec );
  }
  if ( outputPtr->PTFsnrVec )
  {
    XLALDestroyCOMPLEX8Vector( outputPtr->PTFsnrVec );
  }
  if ( outputPtr->PTFA )
  {
    XLALDestroyArray( outputPtr->PTFA );
  }
  if ( outputPtr->PTFMatrix )
  {
    XLALDestroyArray( outputPtr->PTFMatrix );
  }

  /* BCV */
  if (outputPtr->qVecBCV)
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCV) );
    CHECKSTATUSPTR( status );
  }

  if (outputPtr->qVecBCVSpin1)
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin1) );
    CHECKSTATUSPTR( status );
  }

  if (outputPtr->qVecBCVSpin2)
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->qVecBCVSpin2) );
    CHECKSTATUSPTR( status );
  }

  /* AmpCorPPN */
  if( outputPtr->qVecACTD )
  {
    INT4 k;
    for( k = 0 ; k < NACTDVECS; ++k )
    {
      LALCDestroyVector( status->statusPtr, &( outputPtr->qVecACTD[k] ) );
      CHECKSTATUSPTR( status );
    }
    LALFree( outputPtr->qVecACTD );
  }

  /* destroy workspace vector for optimal filter: freq domain */
  if (outputPtr->qtildeVec)
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ) ;
    CHECKSTATUSPTR( status );
  }
  if  (outputPtr->qtildeVecBCV)
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVecBCV) );
    CHECKSTATUSPTR( status );
  }

  if (outputPtr->qtildeVecBCVSpin1)
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVecBCVSpin1) );
    CHECKSTATUSPTR( status );
  }

  if (outputPtr->qtildeVecBCVSpin2)
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVecBCVSpin2) );
    CHECKSTATUSPTR( status );
  }

  if ( outputPtr->qtildeVecACTD )
  {
    INT4 k;
    for( k = 0 ; k < NACTDVECS; ++k )
    {
      LALCDestroyVector( status->statusPtr, &( outputPtr->qtildeVecACTD[k] ) );
      CHECKSTATUSPTR( status );
    }
    LALFree( outputPtr->qtildeVecACTD );
  }

  /* destroy workspace vector for chisq filter */
  if ( outputPtr->chisqVec )
  {
    LALDestroyVector( status->statusPtr, &(outputPtr->chisqVec) );
    CHECKSTATUSPTR( status );
  }

  /*
   *
   * free the chisq structures
   *
   */


  /* parameter structure */
  LALFree( outputPtr->chisqParams );

  /* input structures */
  LALFree( outputPtr->chisqInput );
  if ( outputPtr->chisqInputBCV )
  {
    LALFree( outputPtr->chisqInputBCV );
  }


  /*
   *
   * destroy vector to store snrsq, if it exists
   *
   */


  if ( outputPtr->rhosqVec )
  {
    LALDestroyVector( status->statusPtr, &(outputPtr->rhosqVec->data) );
    CHECKSTATUSPTR( status );

    LALFree( outputPtr->rhosqVec );
  }

  if ( outputPtr->cVec )
  {
    LALCDestroyVector( status->statusPtr, &(outputPtr->cVec->data) );
    CHECKSTATUSPTR( status );

    LALFree( outputPtr->cVec );
  }

  /*
   *
   * free memory for the FindChirpFilterParams
   *
   */


  LALFree( outputPtr );
  *output = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
