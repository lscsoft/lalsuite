/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilter.c
 *
 * Author: Brown, D. A., BCV-Modifications by Messaritaki E.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpFilterCV">
Author: Brown D. A., BCV-Modifications by Messaritaki E.
$Id$
</lalVerbatim>

<lalLaTeX>
\input{FindChirpFilterCDoc}

\vfill{\footnotesize\input{FindChirpFilterCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>

double rint(double x);

NRCSID (FINDCHIRPFILTERC, "$Id$");


/* <lalVerbatim file="FindChirpFilterCP"> */
void
LALCreateFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
{
  FindChirpFilterInput         *outputPtr;

  INITSTATUS( status, "LALCreateFindChirpFilterInput", FINDCHIRPFILTERC );
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


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



/* <lalVerbatim file="FindChirpFilterCP"> */
void
LALDestroyFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output
    )
/* </lalVerbatim> */
{
  FindChirpFilterInput         *outputPtr;

  INITSTATUS( status, "LALDestroyFindChirpFilterInput", FINDCHIRPFILTERC );
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

  /* destroy the chirp template structure */
  LALFree( outputPtr->fcTmplt );

  /* destroy the filter input structure */
  LALFree( outputPtr );
  *output = NULL;


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}


/* <lalVerbatim file="FindChirpFilterCP"> */
void
LALFindChirpFilterInit (
    LALStatus                  *status,
    FindChirpFilterParams     **output,
    FindChirpInitParams        *params
    )
/* </lalVerbatim> */
{
  FindChirpFilterParams        *outputPtr;

  INITSTATUS( status, "LALFindChirpFilterInit", FINDCHIRPFILTERC );
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
      TRY( LALDestroyVector (status->statusPtr, &(outputPtr->chisqVec) ), 
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


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



/* <lalVerbatim file="FindChirpFilterCP"> */
void
LALFindChirpFilterFinalize (
    LALStatus                  *status,
    FindChirpFilterParams     **output
    )
/* </lalVerbatim> */
{
  FindChirpFilterParams        *outputPtr;

  INITSTATUS( status, "LALFindChirpFilterFinalize", FINDCHIRPFILTERC );
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



/* <lalVerbatim file="FindChirpFilterCP"> */
void
LALFindChirpFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 j, k;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  REAL4                 deltaT;
  REAL4                 norm;
  REAL4                 modqsqThresh;
  REAL4                 mismatch;
  REAL4                 chisqThreshFac;
  REAL4                 modChisqThresh;
  UINT4                 numChisqBins;
  UINT4                 eventStartIdx = 0;
  REAL4                 chirpTime     = 0;
  BOOLEAN               haveChisq     = 0;
  COMPLEX8             *qtilde        = NULL;
  COMPLEX8             *q             = NULL;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;
  LALMSTUnitsAndAcc     gmstUnits;

  INITSTATUS( status, "LALFindChirpFilter", FINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh > 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );
  ASSERT( params->chisqThresh > 0, status,
      FINDCHIRPH_ECHIT, FINDCHIRPH_MSGECHIT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT( params->qVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->qVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->qtildeVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInput, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec ) 
  {
    ASSERT( params->rhosqVec->data->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status, 
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* if a chisqVec vector has been created, check we can store data in it */
  if ( params->chisqVec ) 
  {
    ASSERT( params->chisqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->tmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure the filter has been initialized for the correct approximant */
  if ( params->approximant != TaylorF2 )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }

  /* make sure that the template and the segment are both stationary phase */
  ASSERT( input->fcTmplt->approximant == TaylorF2, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == TaylorF2, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );


  /*
   *
   * point local pointers to input and output pointers
   *
   */


  /* workspace vectors */
  q = params->qVec->data;
  qtilde = params->qtildeVec->data;

  /* template and data */
  inputData = input->segment->data->data->data;
  tmpltSignal = input->fcTmplt->data->data;
  deltaT = params->deltaT;

  /* the length of the chisq bin vec is the number of bin   */
  /* _boundaries_ so the number of chisq bins is length - 1 */
  numChisqBins = input->segment->chisqBinVec->length - 1;

  /* number of points in a segment */
  numPoints = params->qVec->length;

  /* set the gmst units and strictness */
  gmstUnits.units = MST_HRS;
  gmstUnits.accuracy = LALLEAPSEC_STRICT;


  /*
   *
   * compute viable search regions in the snrsq vector
   *
   */


  {
    /* calculate the length of the chirp */
    REAL4 eta = input->tmplt->eta;
    REAL4 m1 = input->tmplt->mass1;
    REAL4 m2 = input->tmplt->mass2;
    REAL4 fmin = input->segment->fLow;
    REAL4 m = m1 + m2;
    REAL4 c0 = 5*m*LAL_MTSUN_SI/(256*eta);
    REAL4 c2 = 743.0/252.0 + eta*11.0/3.0;
    REAL4 c3 = -32*LAL_PI/3;
    REAL4 c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
    REAL4 x  = pow(LAL_PI*m*LAL_MTSUN_SI*fmin, 1.0/3.0);
    REAL4 x2 = x*x;
    REAL4 x3 = x*x2;
    REAL4 x4 = x2*x2;
    REAL4 x8 = x4*x4;
    chirpTime = c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8;

    deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );

    /* ignore corrupted data at start and end */
    ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

    if ( lalDebugLevel | LALINFO )
    {
      CHAR infomsg[256];

      LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
          "m1 = %e m2 = %e => %e seconds => %d points\n"
          "invSpecTrunc = %d => ignoreIndex = %d\n", 
          m1, m2, chirpTime, deltaEventIndex, 
          input->segment->invSpecTrunc, ignoreIndex );
      LALInfo( status, infomsg );
    }

    /* XXX check that we are not filtering corrupted data XXX */
    /* XXX this is hardwired to 1/4 segment length        XXX */
    if ( ignoreIndex > numPoints / 4 )
    {
      ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
    }
    /* XXX reset ignoreIndex to one quarter of a segment XXX */
    ignoreIndex = numPoints / 4;
  }

  if ( lalDebugLevel | LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg), 
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, infomsg );
  }

  /*
   *
   * compute qtilde and q
   *
   */


  memset( qtilde, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r = inputData[k].re;
    REAL4 s = inputData[k].im;
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = 0 - tmpltSignal[k].im;       /* note complex conjugate */

    qtilde[k].re = r*x - s*y;
    qtilde[k].im = r*y + s*x;
  }

  /* qtilde negative frequency only: not DC or nyquist */
  if ( params->computeNegFreq )
  {
    for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
    {
      REAL4 r = inputData[k].re;
      REAL4 s = inputData[k].im;
      REAL4 x = tmpltSignal[k].re;
      REAL4 y = 0 - tmpltSignal[k].im;     /* note complex conjugate */

      qtilde[k].re = r*x - s*y;
      qtilde[k].im = r*y + s*x;
    }
  }

  /* inverse fft to get q */
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec, params->qtildeVec, 
      params->invPlan );
  CHECKSTATUSPTR( status );


  /* 
   *
   * calculate signal to noise squared 
   *
   */


  /* if full snrsq vector is required, set it to zero */
  if ( params->rhosqVec )
    memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );

  /* normalisation */
  params->norm = norm = 
    4.0 * (deltaT / (REAL4)numPoints) / input->segment->segNorm;

  /* normalised snr threhold */
  modqsqThresh = params->rhosqThresh / norm;

  /* we threshold on the "modified" chisq threshold computed from       */
  /*   chisqThreshFac = delta^2 * norm / p                              */
  /*   rho^2 = norm * modqsq                                            */
  /*                                                                    */
  /* So we actually threshold on                                        */
  /*                                                                    */
  /*    r^2 < chisqThresh * ( 1 + modqsq * chisqThreshFac )             */
  /*                                                                    */
  /* which is the same as thresholding on                               */
  /*    r^2 < chisqThresh * ( 1 + rho^2 * delta^2 / p )                 */
  /* and since                                                          */
  /*    chisq = p r^2                                                   */
  /* this is equivalent to thresholding on                              */
  /*    chisq < chisqThresh * ( p + rho^2 delta^2 )                     */
  /*                                                                    */
  /* The raw chisq is stored in the database. this quantity is chisq    */
  /* distributed with 2p-2 degrees of freedom.                          */
  mismatch = 1.0 - input->tmplt->minMatch;
  chisqThreshFac = norm * mismatch * mismatch / (REAL4) numChisqBins;

  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec ) 
  {
    memcpy( params->rhosqVec->name, input->segment->data->name,
        LALNameLength * sizeof(CHAR) );
    memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch), 
        sizeof(LIGOTimeGPS) );
    params->rhosqVec->deltaT = input->segment->deltaT;

    for ( j = 0; j < numPoints; ++j )
    {
      REAL4 modqsq = q[j].re * q[j].re + q[j].im * q[j].im;
      params->rhosqVec->data->data[j] = norm * modqsq;
    }
  }

  /* look for an events in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
    REAL4 modqsq = q[j].re * q[j].re + q[j].im * q[j].im;

    /* if snrsq exceeds threshold at any point */
    if ( modqsq > modqsqThresh )
    {
      /* compute chisq vector if it does not exist and we want it */
      if ( ! haveChisq  && input->segment->chisqBinVec->length )
      {
        memset( params->chisqVec->data, 0, 
            params->chisqVec->length * sizeof(REAL4) );

        /* pointers to chisq input */
        params->chisqInput->qtildeVec = params->qtildeVec;
        params->chisqInput->qVec      = params->qVec;

        /* pointer to the chisq bin vector in the segment */
        params->chisqParams->chisqBinVec = input->segment->chisqBinVec;
        params->chisqParams->norm        = norm;
#if 0
        params->chisqParams->bankMatch   = input->tmplt->minMatch;
#endif

        /* compute the chisq threshold: this is slow! */
        LALFindChirpChisqVeto( status->statusPtr, params->chisqVec, 
            params->chisqInput, params->chisqParams );
        CHECKSTATUSPTR (status); 

        haveChisq = 1;
      }

      /* if we have don't have a chisq or the chisq drops below the       */
      /* modified chisq threshold, start processing events                */
      if ( ! input->segment->chisqBinVec->length ||
          params->chisqVec->data[j] < 
          (params->chisqThresh * ( 1.0 + modqsq * chisqThreshFac )) )
      {
        if ( ! *eventList )
        {
          /* store the start of the crossing */
          eventStartIdx = j;

          /* if this is the first event, start the list */
          thisEvent = *eventList = (SnglInspiralTable *) 
            LALCalloc( 1, sizeof(SnglInspiralTable) );
          if ( ! thisEvent )
          {
            ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
          }

          /* record the data that we need for the clustering algorithm */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = modqsq;
        }
        else if ( params->maximiseOverChirp &&
            j <= thisEvent->end_time.gpsSeconds + deltaEventIndex &&
            modqsq > thisEvent->snr )
        {
          /* if this is the same event, update the maximum */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = modqsq;
        }
        else if ( j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
            ! params->maximiseOverChirp )
        {
          /* clean up this event */
          SnglInspiralTable *lastEvent;
          INT8               timeNS;
          INT4               timeIndex = thisEvent->end_time.gpsSeconds;

          /* set the event LIGO GPS time of the event */
          timeNS = 1000000000L * 
            (INT8) (input->segment->data->epoch.gpsSeconds);
          timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
          timeNS += (INT8) (1e9 * timeIndex * deltaT);
          thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
          thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
          LALGPStoGMST1( status->statusPtr, &(thisEvent->end_time_gmst),
              &(thisEvent->end_time), &gmstUnits );
          CHECKSTATUSPTR( status );

          /* set the impuse time for the event */
          thisEvent->template_duration = (REAL8) chirpTime;

          /* record the ifo and channel name for the event */
          strncpy( thisEvent->ifo, input->segment->data->name, 
              2 * sizeof(CHAR) );
          strncpy( thisEvent->channel, input->segment->data->name + 3,
              (LALNameLength - 3) * sizeof(CHAR) );
          thisEvent->impulse_time = thisEvent->end_time;

          /* record the coalescence phase of the chirp */
          thisEvent->coa_phase = (REAL4) 
            atan2( q[timeIndex].im, q[timeIndex].re ); 

          /* copy the template into the event */
          thisEvent->mass1  = (REAL4) input->tmplt->mass1;
          thisEvent->mass2  = (REAL4) input->tmplt->mass2;
          thisEvent->mchirp = (REAL4) input->tmplt->chirpMass;
          thisEvent->eta    = (REAL4) input->tmplt->eta;
          thisEvent->tau0   = (REAL4) input->tmplt->t0;
          thisEvent->tau2   = (REAL4) input->tmplt->t2;
          thisEvent->tau3   = (REAL4) input->tmplt->t3;
          thisEvent->tau4   = (REAL4) input->tmplt->t4;
          thisEvent->tau5   = (REAL4) input->tmplt->t5;
          thisEvent->ttotal = (REAL4) input->tmplt->tC;

          /* set the type of the template used in the analysis */
          LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
              "FindChirpSPtwoPN" );

          /* set snrsq, chisq, sigma and effDist for this event */
          if ( input->segment->chisqBinVec->length )
          {
            /* we store chisq distributed with 2p - 2 degrees of freedom */
            /* in the database. params->chisqVec->data = r^2 = chisq / p */
            /* so we multiply r^2 by p here to get chisq                 */
            thisEvent->chisq = 
              params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
            thisEvent->chisq_dof = numChisqBins;
          }
          else
          {
            thisEvent->chisq     = 0;
            thisEvent->chisq_dof = 0;
          }
          thisEvent->sigmasq = norm * input->segment->segNorm * 
            input->segment->segNorm * input->fcTmplt->tmpltNorm;
          thisEvent->eff_distance = 
            (input->fcTmplt->tmpltNorm * input->segment->segNorm * 
             input->segment->segNorm) / thisEvent->snr;
          thisEvent->eff_distance = sqrt( thisEvent->eff_distance );

          thisEvent->snr *= norm;
          thisEvent->snr = sqrt( thisEvent->snr );

          /* compute the time since the snr crossing */
          thisEvent->event_duration = (REAL8) timeIndex - (REAL8) eventStartIdx;
          thisEvent->event_duration *= (REAL8) deltaT;

          /* store the start of the crossing */
          eventStartIdx = j;

          /* allocate memory for the newEvent */
          lastEvent = thisEvent;

          lastEvent->next = thisEvent = (SnglInspiralTable *) 
            LALCalloc( 1, sizeof(SnglInspiralTable) );
          if ( ! lastEvent->next )
          {
            ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
          }

          /* stick minimal data into the event */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = modqsq;
        }
      }
    }
  }


  /* 
   *
   * clean up the last event if there is one
   *
   */


  if ( thisEvent )
  {
    INT8           timeNS;
    INT4           timeIndex = thisEvent->end_time.gpsSeconds;

    /* set the event LIGO GPS time of the event */
    timeNS = 1000000000L * 
      (INT8) (input->segment->data->epoch.gpsSeconds);
    timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
    timeNS += (INT8) (1e9 * timeIndex * deltaT);
    thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
    thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
    LALGPStoGMST1( status->statusPtr, &(thisEvent->end_time_gmst),
        &(thisEvent->end_time), &gmstUnits );
    CHECKSTATUSPTR( status );

    /* set the impuse time for the event */
    thisEvent->template_duration = (REAL8) chirpTime;

    /* record the ifo name for the event */
    strncpy( thisEvent->ifo, input->segment->data->name, 
        2 * sizeof(CHAR) );
    strncpy( thisEvent->channel, input->segment->data->name + 3,
        (LALNameLength - 3) * sizeof(CHAR) );
    thisEvent->impulse_time = thisEvent->end_time;

    /* record the coalescence phase of the chirp */
    if ( q[timeIndex].re == 0 )
    {
      if ( q[timeIndex].im >= 0 )
      {
        thisEvent->coa_phase = LAL_PI / 2.0;
      }
      else
      {
        thisEvent->coa_phase = - LAL_PI / 2.0;
      }
    }
    else
    {
      thisEvent->coa_phase = (REAL4) 
        atan( q[timeIndex].im / q[timeIndex].re );
    }

    /* copy the template into the event */
    thisEvent->mass1  = (REAL4) input->tmplt->mass1;
    thisEvent->mass2  = (REAL4) input->tmplt->mass2;
    thisEvent->mchirp = (REAL4) input->tmplt->chirpMass;
    thisEvent->eta    = (REAL4) input->tmplt->eta;
    thisEvent->tau0   = (REAL4) input->tmplt->t0;
    thisEvent->tau2   = (REAL4) input->tmplt->t2;
    thisEvent->tau3   = (REAL4) input->tmplt->t3;
    thisEvent->tau4   = (REAL4) input->tmplt->t4;
    thisEvent->tau5   = (REAL4) input->tmplt->t5;
    thisEvent->ttotal = (REAL4) input->tmplt->tC;

    /* set the type of the template used in the analysis */
    LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
        "FindChirpSPtwoPN" );

    /* set snrsq, chisq, sigma and effDist for this event */
    if ( input->segment->chisqBinVec->length )
    {
      /* we store chisq distributed with 2p - 2 degrees of freedom */
      /* in the database. params->chisqVec->data = r^2 = chisq / p */
      /* so we multiply r^2 by p here to get chisq                 */
      thisEvent->chisq = 
        params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
      thisEvent->chisq_dof = numChisqBins;
    }
    else
    {
      thisEvent->chisq     = 0;
      thisEvent->chisq_dof = 0;
    }
    thisEvent->sigmasq = norm * input->segment->segNorm * 
      input->segment->segNorm * input->fcTmplt->tmpltNorm;
    thisEvent->eff_distance = 
      (input->fcTmplt->tmpltNorm * input->segment->segNorm * 
       input->segment->segNorm) / thisEvent->snr;
    thisEvent->eff_distance = sqrt( thisEvent->eff_distance );
    thisEvent->snr *= norm;
    thisEvent->snr = sqrt( thisEvent->snr );

    /* compute the time since the snr crossing */
    thisEvent->event_duration = (REAL8) timeIndex - (REAL8) eventStartIdx;
    thisEvent->event_duration *= (REAL8) deltaT;
  }    


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}



/* <lalVerbatim file="FindChirpFilterCP"> */
void
LALFindChirpBCVFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 j, k;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex;
  UINT4                 ignoreIndex;
  REAL4                 deltaT;
  REAL4                 norm;
  REAL4                 modqsqThresh;
  REAL4                 rhosqThresh;
  REAL4                 mismatch;
  REAL4                 chisqThreshFac;
  REAL4                 modChisqThresh;
  UINT4                 numChisqBins;
  UINT4                 eventStartIdx = 0;
  REAL4                 chirpTime     = 0;
  BOOLEAN               haveChisq     = 0;
  COMPLEX8             *qtilde        = NULL; 
  COMPLEX8             *qtildeBCV     = NULL; 
  COMPLEX8             *q             = NULL; 
  COMPLEX8             *qBCV          = NULL;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *inputDataBCV  = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;
  LALMSTUnitsAndAcc     gmstUnits;
  REAL4                 a1;
  REAL4                 b1;                  
  REAL4                 b2;                  
  REAL4                 templateNorm;
  REAL4                 modqsq;
  REAL4                 omega;
  REAL4                 Num1;
  REAL4                 Num2;
  REAL4                 Den1;
  REAL4                 Den2;
  REAL4                 InvTan1;
  REAL4                 InvTan2;

  /*
     FindChirpChisqInput  *chisqInput;
     FindChirpChisqInput  *chisqInputBCV;
   */

  INITSTATUS( status, "LALFindChirpBCVFilter", FINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *    
   * check that the arguments are reasonable
   *          
   */

  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRPH_ENNUL, FINDCHIRPH_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRPH_EDTZO, FINDCHIRPH_MSGEDTZO );
  ASSERT( params->rhosqThresh > 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );
  ASSERT( params->chisqThresh > 0, status,
      FINDCHIRPH_ECHIT, FINDCHIRPH_MSGECHIT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT(params->qVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec->data,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);
  ASSERT(params->qVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVecBCV->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVecBCV->data,status, FINDCHIRPH_ENULL, 
      FINDCHIRPH_MSGENULL);

  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInput,   status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( params->chisqInputBCV,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec )
  {
    ASSERT( params->rhosqVec->data->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* if a chisqVec vector has been created, check we can store data in it */
  if ( params->chisqVec )
  {
    ASSERT( params->chisqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->tmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure the filter has been initialized for the correct approximant */
  if ( params->approximant != BCV )
  {
    ABORT( status, FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  }

  /* make sure that the template and the segment are both BCV */
  ASSERT( input->fcTmplt->approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );


  /*
   *
   * point local pointers to input and output pointers
   *
   */


  /* workspace vectors */
  q    = params->qVec->data;
  qBCV = params->qVecBCV->data;
  qtilde    = params->qtildeVec->data;
  qtildeBCV = params->qtildeVecBCV->data;


  /* template and data */
  inputData     = input->segment->data->data->data;
  inputDataBCV  = input->segment->dataBCV->data->data;
  tmpltSignal   = input->fcTmplt->data->data;
  templateNorm  = input->fcTmplt->tmpltNorm;
  deltaT        = params->deltaT;

  /* the length of the chisq bin vec is the number of bin   */
  /* _boundaries_ so the number of chisq bins is length - 1 */
  numChisqBins = input->segment->chisqBinVec->length - 1;

  /* number of points in a segment */
  numPoints = params->qVec->length;

  /* set the gmst units and strictness */
  gmstUnits.units = MST_HRS;
  gmstUnits.accuracy = LALLEAPSEC_STRICT;


  /*
   *
   * compute viable search regions in the snrsq vector
   *
   */


  {
    /* length of the chirp:                                            */
    /* calculated according to chirplen, using the values of M and eta */
    /* for the BCV tempaltes, as calculated using psi0 and psi3        */ 
    REAL4 psi0 = input->tmplt->psi0;
    REAL4 psi3 = input->tmplt->psi3;
    REAL4 fFinal = input->tmplt->fFinal;

    /* REAL4 eta = input->tmplt->eta; */

    REAL4 m1 = input->tmplt->mass1;
    REAL4 m2 = input->tmplt->mass2;

    REAL4 fmin = input->segment->fLow;
    /* total mass, in seconds */
    REAL4 m =  abs(psi3) / (16.0 * LAL_PI * LAL_PI * psi0) ;
    REAL4 eta = 3.0 / (128.0 * psi0 * pow( (m*LAL_PI), (5.0/3.0)) );
    /*REAL4 c0 = 5*m*LAL_MTSUN_SI/(256*eta);*/
    REAL4 c0 = 5*m/(256*eta);
    REAL4 c2 = 743.0/252.0 + eta*11.0/3.0;
    REAL4 c3 = -32*LAL_PI/3;
    REAL4 c4 = 3058673.0/508032.0 + eta*(5429.0/504.0 + eta*617.0/72.0);
    REAL4 x  = pow(LAL_PI*m*fmin, 1.0/3.0);
    REAL4 x2 = x*x;
    REAL4 x3 = x*x2;
    REAL4 x4 = x2*x2;
    REAL4 x8 = x4*x4;
    REAL4 chirpTime = abs(c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8);

    /* k that corresponds to fFinal, currently not used      */
    /* UINT4 kFinal = floor( numPoints * deltaT * fFinal ); */  
    /* BCV normalization parameters */
    REAL4 a1 = input->segment->a1;
    REAL4 b1 = input->segment->b1;  
    REAL4 b2 = input->segment->b2; 

    deltaEventIndex = (UINT4) rint( (chirpTime / deltaT) + 1.0 );

    /* ignore corrupted data at start and end */
    ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

    if ( lalDebugLevel | LALINFO )
    {
      CHAR infomsg[256];

      LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
          "m1 = %e m2 = %e => %e seconds => %d points\n"
          "invSpecTrunc = %d => ignoreIndex = %d\n",
          m1, m2, chirpTime, deltaEventIndex,
          input->segment->invSpecTrunc, ignoreIndex );
      LALInfo( status, infomsg );
    }

    /* XXX check that we are not filtering corrupted data XXX */
    /* XXX this is hardwired to 1/4 segment length        XXX */
    if ( ignoreIndex > numPoints / 4 )
    {
      ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
    }
    /* XXX reset ignoreIndex to one quarter of a segment XXX */
    ignoreIndex = numPoints / 4;
  }

  if ( lalDebugLevel | LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, infomsg );
  }


  /*
   *
   * compute qtilde, qtildeBCV, and q, qBCV 
   *
   */


  memset( qtilde,    0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCV, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r    = inputData[k].re;
    REAL4 s    = 0.0 - inputData[k].im;    /* note complex conjugate */
    REAL4 rBCV = inputDataBCV[k].re;
    REAL4 sBCV = 0.0 - inputDataBCV[k].im; /* note complex conjugate */
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = tmpltSignal[k].im;     

    qtilde[k].re = r * x - s * y ;
    qtilde[k].im = r * y + s * x ;
    qtildeBCV[k].re = rBCV * x - sBCV * y ;
    qtildeBCV[k].im = rBCV * y + sBCV * x ;
  }

  /* qtilde negative frequency only: not DC or nyquist */
  if ( params->computeNegFreq )
  {
    for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
    {
      REAL4 r = inputData[k].re;
      REAL4 s = 0.0 - inputData[k].im;    /* note complex conjugate */
      REAL4 rBCV = inputDataBCV[k].re;
      REAL4 sBCV = 0.0 - inputDataBCV[k].im; /* note complex conjugate */
      REAL4 x = tmpltSignal[k].re;
      REAL4 y = tmpltSignal[k].im;

      qtilde[k].re = r * x - s * y ;
      qtilde[k].im = r * y + s * x ;
      qtildeBCV[k].re = rBCV * x - sBCV * y ;
      qtildeBCV[k].im = rBCV * y + sBCV * x ;
    }
  }

  /* inverse fft to get q, and qBCV */
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec, 
      params->qtildeVec, params->invPlan );
  CHECKSTATUSPTR( status );
  LALCOMPLEX8VectorFFT( status->statusPtr, params->qVecBCV, 
      params->qtildeVecBCV, params->invPlan );
  CHECKSTATUSPTR( status );




  /* 
   *
   * calculate signal to noise squared 
   *
   */



  /* if full snrsq vector is required, set it to zero */
  if ( params->rhosqVec )
    memset( params->rhosqVec->data->data, 0, numPoints * sizeof( REAL4 ) );

  /* normalisation */
  /*  For the BCV templates, templateNorm is equal to                      */
  /*  (5*eta/96) * (M/pi^2)^(2/3) (Tsun/Dt)^(-1/3) (2*Msun(L)*d/1Mpc)^2    */
  /* which is the square of the normalization factor that multiplies the   */
  /* template                                                              */

  rhosqThresh = params->rhosqThresh;

  params->norm = norm = deltaT / ((REAL4) numPoints) ;


  /* normalized snr threhold */
  modqsqThresh = rhosqThresh / norm ;

  /* we threshold on the "modified" chisq threshold computed from       */
  /*   chisqThreshFac = delta^2 * norm / p                              */
  /*   rho^2 = norm * modqsq                                            */
  /*                                                                    */
  /* So we actually threshold on                                        */
  /*                                                                    */
  /*    r^2 < chisqThresh * ( 1 + modqsq * chisqThreshFac )             */
  /*                                                                    */
  /* which is the same as thresholding on                               */
  /*    r^2 < chisqThresh * ( 1 + rho^2 * delta^2 / p )                 */
  /* and since                                                          */
  /*    chisq = p r^2                                                   */
  /* this is equivalent to thresholding on                              */
  /*    chisq < chisqThresh * ( p + rho^2 delta^2 )                     */
  /*                                                                    */
  /* The raw chisq is stored in the database. this quantity is chisq    */
  /* distributed with 2p-2 degrees of freedom.                          */

  mismatch = 1.0 - input->tmplt->minMatch;
  chisqThreshFac = norm * mismatch * mismatch / (REAL4) numChisqBins;

  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec )
  {
    memcpy( params->rhosqVec->name, input->segment->data->name,
        LALNameLength * sizeof(CHAR) );
    memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch),
        sizeof(LIGOTimeGPS) );
    params->rhosqVec->deltaT = input->segment->deltaT;

    for ( j = 0; j < numPoints; ++j )
    {
      REAL4 modqsqSP  = q[j].re * q[j].re + q[j].im * q[j].im ;
      REAL4 modqsqBCV = qBCV[j].re * qBCV[j].re + qBCV[j].im * qBCV[j].im ;
      REAL4 ImProd = 2.0 * ( - q[j].re * qBCV[j].im + qBCV[j].re * q[j].im ) ;

      REAL4 modqsq = ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) *
        ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) ;

      params->rhosqVec->data->data[j] = norm * modqsq;   
    }
  }


  /* look for an events in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
    REAL4 modqsqSP  = q[j].re * q[j].re + q[j].im * q[j].im ;
    REAL4 modqsqBCV = qBCV[j].re * qBCV[j].re + qBCV[j].im * qBCV[j].im ;
    REAL4 ImProd = 2.0 * ( - q[j].re * qBCV[j].im + qBCV[j].re * q[j].im ) ;

    REAL4 modqsq = ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
        0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) *
      ( 0.5 * sqrt( modqsqSP + modqsqBCV + ImProd ) +
        0.5 * sqrt( modqsqSP + modqsqBCV - ImProd ) ) ;


    /* if snrsq exceeds threshold at any point */
    if ( modqsq > modqsqThresh )                  
    {
      /* compute chisq vector if it does not exist and we want it */
      if ( ! haveChisq  && input->segment->chisqBinVec->length )
      {
        memset( params->chisqVec->data, 0,
            params->chisqVec->length * sizeof(REAL4) );

        /* pointers to chisq input */
        params->chisqInput->qtildeVec    = params->qtildeVec;
        params->chisqInput->qVec         = params->qVec;
        params->chisqInputBCV->qtildeVec = params->qtildeVecBCV;
        params->chisqInputBCV->qVec      = params->qVecBCV;

        /* pointer to the chisq bin vector in the segment */
        params->chisqParams->chisqBinVec    = input->segment->chisqBinVec;
        params->chisqParams->chisqBinVecBCV = input->segment->chisqBinVecBCV;
        params->chisqParams->norm           = norm;
        params->chisqParams->a1             = input->segment->a1 ;
        params->chisqParams->b1             = input->segment->b1 ;
        params->chisqParams->b2             = input->segment->b2 ;
#if 0
        params->chisqParams->bankMatch   = input->tmplt->minMatch;
#endif

        /* compute the chisq threshold: this is slow! */
        LALFindChirpBCVChisqVeto( status->statusPtr, params->chisqVec,
            params->chisqInput, params->chisqInputBCV, params->chisqParams );
        CHECKSTATUSPTR (status);

        haveChisq = 1;
      }


      /* if we don't have a chisq or the chisq drops below the            */
      /* modified chisq threshold, start processing events                */
      if ( ! input->segment->chisqBinVec->length ||
          params->chisqVec->data[j] <
          (params->chisqThresh * ( 1.0 + modqsq * chisqThreshFac )) )
      {
        if ( ! *eventList )
        {
          /* store the start of the crossing */
          eventStartIdx = j;

          /* if this is the first event, start the list */
          thisEvent = *eventList = (SnglInspiralTable *)
            LALCalloc( 1, sizeof(SnglInspiralTable) );
          if ( ! thisEvent )
          {
            ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
          }

          /* record the data that we need for the clustering algorithm */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = modqsq;
        }
        else if ( params->maximiseOverChirp &&
            j <= thisEvent->end_time.gpsSeconds + deltaEventIndex &&
            modqsq > thisEvent->snr )
        {
          /* if this is the same event, update the maximum */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = modqsq;
        }
        else if ( (j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
              ! params->maximiseOverChirp) )
        {
          /* clean up this event */
          SnglInspiralTable *lastEvent;
          INT8               timeNS;
          INT4               timeIndex = thisEvent->end_time.gpsSeconds;

          /* set the event LIGO GPS time of the event */
          timeNS = 1000000000L *
            (INT8) (input->segment->data->epoch.gpsSeconds);
          timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
          timeNS += (INT8) (1e9 * timeIndex * deltaT);
          thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
          thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
          LALGPStoGMST1( status->statusPtr, &(thisEvent->end_time_gmst),
              &(thisEvent->end_time), &gmstUnits );
          CHECKSTATUSPTR( status );

          /* set the impuse time for the event */
          thisEvent->template_duration = (REAL8) chirpTime;

          /* record the ifo and channel name for the event */
          strncpy( thisEvent->ifo, input->segment->data->name,
              2 * sizeof(CHAR) );
          strncpy( thisEvent->channel, input->segment->data->name + 3,
              (LALNameLength - 3) * sizeof(CHAR) );
          thisEvent->impulse_time = thisEvent->end_time;

          /* record coalescence phase and alpha */

          /* calculate the numerators and the denominators */
          Num1 = qBCV[timeIndex].re + q[timeIndex].im ;
          Num2 = qBCV[timeIndex].re - q[timeIndex].re ;
          Den1 = q[timeIndex].re - qBCV[timeIndex].im ;
          Den2 = q[timeIndex].re + qBCV[timeIndex].im ;

          InvTan1 = (REAL4) atan2(Num1, Den1) ;
          InvTan2 = (REAL4) atan2(Num2, Den2);

          thisEvent->coa_phase = 0.5 * InvTan1 - 0.5 * InvTan2 ;
          omega = 0.5 * InvTan1 + 0.5 * InvTan2 ;
          thisEvent->alpha = - input->segment->b2 * tan(omega) / 
            ( input->segment->a1 + input->segment->b1*tan(omega) );
          /* actually record alpha * fcut^(2/3) which must be b/t 0 and 1 */
          thisEvent->alpha *= pow((input->tmplt->fFinal) , 2.0/3.0);   

          /* copy the template into the event */
          thisEvent->psi0   = (REAL4) input->tmplt->psi0; 
          thisEvent->psi3   = (REAL4) input->tmplt->psi3;
          /* chirp mass in units of M_sun */
          thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
            pow( 3.0 / 128.0 / input->tmplt->psi0 , 3.0/5.0 );
          thisEvent->f_final  = (REAL4) input->tmplt->fFinal ;

          /* set the type of the template used in the analysis */
          LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
              "FindChirpBCV" );

          /* set snrsq, chisq, sigma and effDist for this event */
          if ( input->segment->chisqBinVec->length )
          {
            /* we store chisq distributed with 2p - 2 degrees of freedom */
            /* in the database. params->chisqVec->data = r^2 = chisq / p */
            /* so we multiply r^2 by p here to get chisq                 */
            thisEvent->chisq =
              params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
            thisEvent->chisq_dof = numChisqBins;
          }
          else
          {
            thisEvent->chisq     = 0;
            thisEvent->chisq_dof = 0;
          }
#if 0
          thisEvent->sigmasq = norm * input->segment->segNorm *
            input->segment->segNorm * input->fcTmplt->tmpltNorm;
          thisEvent->eff_distance =
            (input->fcTmplt->tmpltNorm * input->segment->segNorm *
             input->segment->segNorm) / thisEvent->snr;
          thisEvent->eff_distance = sqrt( thisEvent->eff_distance );
#endif

          thisEvent->snr *= norm;      
          thisEvent->snr = sqrt( thisEvent->snr );

          /* compute the time since the snr crossing */
          thisEvent->event_duration = 
            (REAL8) timeIndex - (REAL8) eventStartIdx;
          thisEvent->event_duration *= (REAL8) deltaT;

          /* store the start of the crossing */
          eventStartIdx = j;

          /* allocate memory for the newEvent */
          lastEvent = thisEvent;

          lastEvent->next = thisEvent = (SnglInspiralTable *)
            LALCalloc( 1, sizeof(SnglInspiralTable) );
          if ( ! lastEvent->next )
          {
            ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
          }

          /* stick minimal data into the event */
          thisEvent->end_time.gpsSeconds = j;
          thisEvent->snr = modqsq;
        }
      }
    }
  }


  /* 
   *
   * clean up the last event if there is one
   *
   */


  if ( thisEvent )
  {
    INT8           timeNS;
    INT4           timeIndex = thisEvent->end_time.gpsSeconds;

    /* set the event LIGO GPS time of the event */
    timeNS = 1000000000L *
      (INT8) (input->segment->data->epoch.gpsSeconds);
    timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
    timeNS += (INT8) (1e9 * timeIndex * deltaT);
    thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
    thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
    LALGPStoGMST1( status->statusPtr, &(thisEvent->end_time_gmst),
        &(thisEvent->end_time), &gmstUnits );
    CHECKSTATUSPTR( status );

    /* set the impuse time for the event */
    thisEvent->template_duration = (REAL8) chirpTime; 

    /* record the ifo name for the event */
    strncpy( thisEvent->ifo, input->segment->data->name,
        2 * sizeof(CHAR) );
    strncpy( thisEvent->channel, input->segment->data->name + 3,
        (LALNameLength - 3) * sizeof(CHAR) );
    thisEvent->impulse_time = thisEvent->end_time;

    /* record coalescence phase and alpha */

    /* calculate the numerators and the denominators */
    Num1 = qBCV[timeIndex].re + q[timeIndex].im ;
    Num2 = qBCV[timeIndex].re - q[timeIndex].re ;
    Den1 = q[timeIndex].re - qBCV[timeIndex].im ;
    Den2 = q[timeIndex].re + qBCV[timeIndex].im ;

    InvTan1 = (REAL4) atan2(Num1, Den1);
    InvTan2 = (REAL4) atan2(Num2, Den2 );


    thisEvent->coa_phase = 0.5 * InvTan1 - 0.5 * InvTan2 ;
    omega = 0.5 * InvTan1 + 0.5 * InvTan2 ;
    thisEvent->alpha = - input->segment->b2 * tan(omega) /
      ( input->segment->a1 + input->segment->b1*tan(omega) );
    /* actually record alpha * fcut^(2/3) which must be b/t 0 and 1 */
    thisEvent->alpha *= pow( (input->tmplt->fFinal), 2.0/3.0);   


    /* copy the template into the event */
    thisEvent->psi0   = (REAL4) input->tmplt->psi0;   
    thisEvent->psi3   = (REAL4) input->tmplt->psi3;  
    /* chirp mass in units of M_sun */
    thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
      pow( 3.0 / 128.0 / input->tmplt->psi0, 3.0/5.0 );
    thisEvent->f_final  = (REAL4) input->tmplt->fFinal;


    /* set the type of the template used in the analysis */
    LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
        "FindChirpBCV" );

    /* set snrsq, chisq, sigma and effDist for this event */
    if ( input->segment->chisqBinVec->length )
    {
      /* we store chisq distributed with 2p - 2 degrees of freedom */
      /* in the database. params->chisqVec->data = r^2 = chisq / p */
      /* so we multiply r^2 by p here to get chisq                 */
      thisEvent->chisq =
        params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
      thisEvent->chisq_dof = numChisqBins;
    }
    else
    {
      thisEvent->chisq     = 0;
      thisEvent->chisq_dof = 0;
    }
#if 0    
    thisEvent->sigmasq = norm * input->segment->segNorm *
      input->segment->segNorm * input->fcTmplt->tmpltNorm;
    thisEvent->eff_distance =
      (input->fcTmplt->tmpltNorm * input->segment->segNorm *
       input->segment->segNorm) / thisEvent->snr;
    thisEvent->eff_distance = sqrt( thisEvent->eff_distance ); 
#endif      

    thisEvent->snr *=  norm ;   
    thisEvent->snr = sqrt( thisEvent->snr );

    /* compute the time since the snr crossing */
    thisEvent->event_duration = (REAL8) timeIndex - (REAL8) eventStartIdx;
    thisEvent->event_duration *= (REAL8) deltaT;
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
