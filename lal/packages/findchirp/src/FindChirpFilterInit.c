/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilterInit.c
 *
 * Author: Brown, D. A., BCV-Modifications by Messaritaki E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpFilterInitCV">
Author: Brown D. A., BCV-Modifications by Messaritaki E.
$Id$
</lalVerbatim>

<lalLaTeX>
\input{FindChirpFilterInitCDoc}

\vfill{\footnotesize\input{FindChirpFilterInitCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/LALDatatypes.h>

double rint(double x);

NRCSID (FINDCHIRPFILTERINITC, "$Id$");


/* <lalVerbatim file="FindChirpFilterInitCP"> */
void
LALCreateFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
{
  FindChirpFilterInput         *outputPtr;

  INITSTATUS( status, "LALCreateFindChirpFilterInput", FINDCHIRPFILTERINITC );
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

  /* create memory for the chirp template data */
  LALCCreateVector (status->statusPtr, &(outputPtr->fcTmplt->data), 
      params->numPoints);
  BEGINFAIL( status )
  {
    LALFree( outputPtr->fcTmplt );
    outputPtr->fcTmplt = NULL;
    LALFree( *output );
    *output = NULL;
  }
  ENDFAIL( status );

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



/* <lalVerbatim file="FindChirpFilterInitCP"> */
void
LALDestroyFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output
    )
/* </lalVerbatim> */
{
  FindChirpFilterInput         *outputPtr;

  INITSTATUS( status, "LALDestroyFindChirpFilterInput", FINDCHIRPFILTERINITC );
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
  LALCDestroyVector( status->statusPtr, &(outputPtr->fcTmplt->data) );
  CHECKSTATUSPTR( status );

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
    
  
  /* destroy the chirp template structure */
  LALFree( outputPtr->fcTmplt );

  /* destroy the filter input structure */
  LALFree( outputPtr );
  *output = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpFilterInitCP"> */
void
LALFindChirpFilterInit (
    LALStatus                  *status,
    FindChirpFilterParams     **output,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
{
  FindChirpFilterParams        *outputPtr;

  INITSTATUS( status, "LALFindChirpFilterInit", FINDCHIRPFILTERINITC );
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
  if ( params->approximant != TaylorF2 && params->approximant != BCV &&
      params->approximant != BCVSpin )
  {
    ABORT( status, FINDCHIRPH_EUAPX, FINDCHIRPH_MSGEUAPX );
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
      &(outputPtr->invPlan), params->numPoints, 0 );
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
      TRY( LALDestroyVector( status->statusPtr, &(outputPtr->chisqVec) ), 
          status ); 
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ),
          status );
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
      TRY( LALDestroyVector( status->statusPtr, &(outputPtr->chisqVec) ), 
          status ); 
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ),
          status );
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



/* <lalVerbatim file="FindChirpFilterInitCP"> */
void
LALFindChirpFilterFinalize (
    LALStatus                  *status,
    FindChirpFilterParams     **output
    )
/* </lalVerbatim> */
{
  FindChirpFilterParams        *outputPtr;

  INITSTATUS( status, "LALFindChirpFilterFinalize", FINDCHIRPFILTERINITC );
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
  LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) );
  CHECKSTATUSPTR( status );

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

  /* destroy workspace vector for optimal filter: freq domain */
  LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) );
  CHECKSTATUSPTR( status );

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

  /* destroy workspace vector for chisq filter */
  LALDestroyVector( status->statusPtr, &(outputPtr->chisqVec) );
  CHECKSTATUSPTR( status );


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
