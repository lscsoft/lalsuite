/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpData.c
 *
 * Author: Brown D. A., BCV-Modifications: Messaritaki E., BCV-Spin: Jones G.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0 
<lalVerbatim file="FindChirpDataCV">
Author: Brown, D. A., BCV-Modifications: Messaritaki E.
$Id$
</lalVerbatim> 

<lalLaTeX>
\input{FindChirpDataCDoc}

\vfill{\footnotesize\input{FindChirpDataCV}}
</lalLaTeX> 
#endif

#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/AVFactories.h>
#include <lal/DataBuffer.h>
#include <lal/LALInspiral.h>
#include <lal/FindChirp.h>


NRCSID (FINDCHIRPDATAC, "$Id$");

/* <lalVerbatim file="FindChirpDataCP"> */
void
LALFindChirpDataInit (
    LALStatus                  *status,
    FindChirpDataParams       **output,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
{
  FindChirpDataParams          *dataParamPtr;
  REAL4                        *amp;
  REAL4                        *ampBCV;
  REAL8                        *ampBCVSpin1; 
  REAL8                        *ampBCVSpin2; 
  const REAL4                   exponent    = -7.0/6.0;
  const REAL4                   exponentBCV = -1.0/2.0;
  const REAL8                   exponentBCVSpin1 = -7.0/6.0;
  const REAL8                   exponentBCVSpin2 = -2.0/3.0;
  UINT4                         k;

  INITSTATUS( status, "LALFindChirpDataInit", FINDCHIRPDATAC );
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
  ASSERT (params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);

  /* make sure that the number of points in a segment is positive */
  ASSERT (params->numPoints > 0, status, 
      FINDCHIRPH_ENUMZ, FINDCHIRPH_MSGENUMZ);

  /* check that the approximant is of a known type */
  switch ( params->approximant )
  {
    case TaylorT1:
    case TaylorT2:
    case TaylorT3:
    case TaylorF2:
    case GeneratePPN:
    case PadeT1:
    case EOB:
    case FindChirpSP:
    case BCV:
    case BCVSpin:
      break;
    default:
      ABORT( status, FINDCHIRPH_EUAPX, FINDCHIRPH_MSGEUAPX );
      break;
  }


  /*
   *
   * allocate memory for the FindChirpDataParams
   *
   */


  /* create the output structure */
  dataParamPtr = *output = (FindChirpDataParams *)
    LALCalloc( 1, sizeof(FindChirpDataParams) );
  if ( ! dataParamPtr )
  {
    ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
  }

  /* record the type of waveform that we are initializing for */
  dataParamPtr->approximant = params->approximant;


  /*
   *
   * allocate and fill vector for f^(-7/6)
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
  {
    amp[k] = pow( ((REAL4) k / (REAL4)params->numPoints), exponent );
  }


  /*
   *
   * for the BCV templates, allocate and fill vector for f^(-1/2)
   *
   */


  /* only do this if a BCV waveform is requested */
  if ( params->approximant == BCV )
  {
    LALCreateVector( status->statusPtr, &dataParamPtr->ampVecBCV,
        (params->numPoints)/2 + 1 );
    BEGINFAIL( status )
    {
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), 
          status );
      LALFree( dataParamPtr );
      *output = NULL;
    }
    ENDFAIL( status );

    ampBCV = dataParamPtr->ampVecBCV->data;
    ampBCV[0] = 0.0;

    for ( k = 1; k < dataParamPtr->ampVecBCV->length; ++k )
    {
      ampBCV[k] = pow( ((REAL4) k / (REAL4)params->numPoints), exponentBCV );
    }
  }


  /*
   *
   * for the BCVSpin templates, allocate and fill vector for f^(-5/3) 
   *
   */


  /* only do this if a BCVSpin waveform is requested */
  if ( params->approximant == BCVSpin )
  {
    LALDCreateVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin1,
        (params->numPoints)/2 + 1 );
    BEGINFAIL( status )
    {
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ),
          status );
      LALFree( dataParamPtr );
      *output = NULL;
    }
    ENDFAIL( status );

    ampBCVSpin1 = dataParamPtr->ampVecBCVSpin1->data;
    ampBCVSpin1[0] = 0.0;

    for ( k = 1; k < dataParamPtr->ampVecBCVSpin1->length; ++k )
    {
      ampBCVSpin1[k] = pow( ((REAL4) k / (REAL4)params->numPoints), 
          exponentBCVSpin1 );
    }
  }


  /*
   *
   * for the BCVSpin templates, allocate and fill vector for f^(-2/3) 
   *
   */


  /* only do this if a BCVSpin waveform is requested */
  if ( params->approximant == BCVSpin )
  {
    LALDCreateVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin2,
        (params->numPoints)/2 + 1 );
    BEGINFAIL( status )
    {
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ),
          status );
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin1 ),
          status );
      LALFree( dataParamPtr );
      *output = NULL;
    }
    ENDFAIL( status );

    ampBCVSpin2 = dataParamPtr->ampVecBCVSpin2->data;
    ampBCVSpin2[0] = 0.0;

    for ( k = 1; k < dataParamPtr->ampVecBCVSpin2->length; ++k )
    {
      ampBCVSpin2[k] = pow( ((REAL4) k / (REAL4)params->numPoints), 
          exponentBCVSpin2 );
    }
  }


  /*
   *
   * create fft plans and workspace vectors
   *
   */


  /* foward fft plan */
  LALCreateForwardRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan, 
      params->numPoints, 0 );
  BEGINFAIL( status )
  {
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), 
        status );
    if ( dataParamPtr->ampVecBCV )
    {
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCV ), 
          status);
    }
    if ( dataParamPtr->ampVecBCVSpin1 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin1 ),
          status);
    }
    if ( dataParamPtr->ampVecBCVSpin2 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin2 ),
          status);
    }
    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* inverse fft plan */
  LALCreateReverseRealFFTPlan( status->statusPtr, &dataParamPtr->invPlan, 
      params->numPoints, 0 );
  BEGINFAIL( status )
  {
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), 
        status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), 
        status );
    if ( dataParamPtr->ampVecBCV )
    {
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCV ), 
          status );
    }
    if ( dataParamPtr->ampVecBCVSpin1 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin1 ),
          status);
    }
    if ( dataParamPtr->ampVecBCVSpin2 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin2 ),
          status);
    }
    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* workspace vector w: time domain */
  LALCreateVector( status->statusPtr, &dataParamPtr->wVec, 
      params->numPoints );
  BEGINFAIL( status )
  {
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->invPlan ), 
        status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), 
        status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), 
        status );
    if ( dataParamPtr->ampVecBCV )
    {
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCV ),
          status );
    }
    if ( dataParamPtr->ampVecBCVSpin1 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin1 ),
          status);
    }
    if ( dataParamPtr->ampVecBCVSpin2 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin2 ),
          status);
    }
    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* workspace vector w: freq domain */
  LALCCreateVector( status->statusPtr, &dataParamPtr->wtildeVec, 
      params->numPoints/2 + 1 );
  BEGINFAIL( status )
  {
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->wVec ), 
        status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->invPlan ), 
        status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), 
        status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), 
        status );
    if ( dataParamPtr->ampVecBCV )
    {
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCV ), 
          status );
    }
    if ( dataParamPtr->ampVecBCVSpin1 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin1 ),
          status);
    }
    if ( dataParamPtr->ampVecBCVSpin2 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin2 ),
          status);
    }
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
    TRY( LALCDestroyVector( status->statusPtr, &dataParamPtr->wtildeVec), 
        status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->wVec ), 
        status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->invPlan ), 
        status ); 
    TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), 
        status );
    TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), 
        status );
    if ( dataParamPtr->ampVecBCV )
    {
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCV ), 
          status);
    }
    if ( dataParamPtr->ampVecBCVSpin1 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin1 ),
          status);
    }
    if ( dataParamPtr->ampVecBCVSpin2 )
    {
      TRY( LALDDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCVSpin2 ),
          status);
    }
    LALFree( dataParamPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* additional BCV template power vector */
  if ( params->approximant == BCV )
  {
    LALCreateVector( status->statusPtr, &dataParamPtr->tmpltPowerVecBCV,
        params->numPoints/2 + 1 );
    BEGINFAIL( status )
    {
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->tmpltPowerVec ),
          status );
      TRY( LALCDestroyVector( status->statusPtr, &dataParamPtr->wtildeVec ), 
          status );
      TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->wVec ), 
          status );
      TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->invPlan ), 
          status ); 
      TRY( LALDestroyRealFFTPlan( status->statusPtr, &dataParamPtr->fwdPlan ), 
          status );
      TRY(LALDestroyVector( status->statusPtr, &dataParamPtr->ampVec ), 
          status ); 
      if ( dataParamPtr->ampVecBCV )
      {
        TRY( LALDestroyVector( status->statusPtr, &dataParamPtr->ampVecBCV ), 
            status );
      }
      if ( dataParamPtr->ampVecBCVSpin1 )
      {
        TRY( LALDDestroyVector( status->statusPtr, 
              &dataParamPtr->ampVecBCVSpin1 ), status );
      }
      if ( dataParamPtr->ampVecBCVSpin2 )
      {
        TRY( LALDDestroyVector( status->statusPtr, 
              &dataParamPtr->ampVecBCVSpin2 ), status );
      }
      LALFree( dataParamPtr );
      *output = NULL;
    }
    ENDFAIL( status );
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



/* <lalVerbatim file="FindChirpDataCP"> */
void
LALFindChirpDataFinalize (
    LALStatus                  *status,
    FindChirpDataParams       **output
    )
/* </lalVerbatim> */
{
  FindChirpDataParams          *dataParamPtr;

  INITSTATUS( status, "LALFindChirpDataFinalize", FINDCHIRPDATAC );
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
   * destroy fft plans and workspace vectors
   *
   */


  /* local pointer to structure */
  dataParamPtr = *output;

  LALDestroyRealFFTPlan (status->statusPtr, &dataParamPtr->fwdPlan);
  CHECKSTATUSPTR (status);

  LALDestroyRealFFTPlan (status->statusPtr, &dataParamPtr->invPlan);
  CHECKSTATUSPTR (status);

  LALDestroyVector (status->statusPtr, &dataParamPtr->wVec);
  CHECKSTATUSPTR (status);

  LALCDestroyVector (status->statusPtr, &dataParamPtr->wtildeVec);
  CHECKSTATUSPTR (status);

  LALDestroyVector (status->statusPtr, &dataParamPtr->tmpltPowerVec);
  CHECKSTATUSPTR (status);

  if ( dataParamPtr->tmpltPowerVecBCV )
  {
    LALDestroyVector (status->statusPtr, &dataParamPtr->tmpltPowerVecBCV);
    CHECKSTATUSPTR (status);
  }


  /*
   *
   * destroy vectors for exponent of amplitude
   *
   */


  LALDestroyVector (status->statusPtr, &dataParamPtr->ampVec);
  CHECKSTATUSPTR (status);

  if ( dataParamPtr->ampVecBCV )
  {
    LALDestroyVector (status->statusPtr, &dataParamPtr->ampVecBCV);
    CHECKSTATUSPTR (status);
  }

  if ( dataParamPtr->ampVecBCVSpin1 )
  {
    LALDDestroyVector (status->statusPtr, &dataParamPtr->ampVecBCVSpin1);
    CHECKSTATUSPTR (status);
  }

  if ( dataParamPtr->ampVecBCVSpin2 )
  {
    LALDDestroyVector (status->statusPtr, &dataParamPtr->ampVecBCVSpin2);
    CHECKSTATUSPTR (status);
  }


  /*
   *
   * free memory for the FindChirpDataParams
   *
   */


  LALFree (dataParamPtr);
  *output = NULL;


  DETATCHSTATUSPTR (status);
  RETURN (status);
}
