/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilter.c
 *
 * Author: Brown, D. A.
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/FindChirp.h>

NRCSID (FINDCHIRPFILTERC, "$Id$");


void
LALCreateFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output,
    FindChirpInitParams        *params
    )
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
  ASSERT( output, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );  
  ASSERT( !*output, status, FINDCHIRP_ENNUL, FINDCHIRP_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* make sure that the number of points is positive */
  ASSERT( params->numPoints > 0, status, 
      FINDCHIRP_ENUMZ, FINDCHIRP_MSGENUMZ );


  /*
   *
   * create the findchirp filter input structure
   *
   */


  /* create the output structure */
  outputPtr = *output = (FindChirpFilterInput *)
    LALMalloc( sizeof(FindChirpFilterInput) );
  ASSERT( outputPtr, status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
  memset( outputPtr, 0, sizeof(FindChirpFilterInput) );

  /* don't need this because of the memset above */
  /* outputPtr->tmplt   = NULL; */
  /* outputPtr->segment = NULL; */

  /* create memory for the chirp template structure */
  outputPtr->fcTmplt = (FindChirpTemplate *)
    LALMalloc( sizeof(FindChirpTemplate) );
  if ( !outputPtr->fcTmplt )
  {
    LALFree( *output );
    *output = NULL;
    ABORT( status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
  }
  memset( outputPtr->fcTmplt, 0, sizeof(FindChirpTemplate) );

  /* don't need this because of the memset above */
  /* outputPtr->fcTmplt->tmpltNorm = 0.0; */
  /* outputPtr->fcTmplt->data = NULL; */

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



void
LALDestroyFindChirpInput (
    LALStatus                  *status,
    FindChirpFilterInput      **output
    )
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
  ASSERT( output, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( *output, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );


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
    

void
LALFindChirpFilterInit (
    LALStatus                  *status,
    FindChirpFilterParams     **output,
    FindChirpInitParams        *params
    )
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
  ASSERT( output, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( !*output, status, FINDCHIRP_ENNUL, FINDCHIRP_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* make sure that the number of points in a segment is positive */
  ASSERT( params->numPoints > 0,  status, 
      FINDCHIRP_ENUMZ, FINDCHIRP_MSGENUMZ );


  /*
   *
   * allocate memory for the FindChirpFilterParams
   *
   */


  /* create the output structure */
  outputPtr = *output = (FindChirpFilterParams *)
    LALMalloc( sizeof(FindChirpFilterParams) );
  ASSERT( outputPtr, status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
  memset( outputPtr, 0, sizeof(FindChirpFilterParams) );

  /* don't need this because of the memset above */
  /* set contents to reasonable values */
  /* outputPtr->deltaT         = 0; */
  /* outputPtr->rhosqThresh    = 0; */
  /* outputPtr->computeNegFreq = 0; */
  /* outputPtr->qVec           = NULL; */
  /* outputPtr->qtildeVec      = NULL; */
  /* outputPtr->invPlan        = NULL; */
  /* outputPtr->chisqVec       = NULL; */
  /* outputPtr->rhosqVec       = NULL; */

  /* create memory for the chisq parameters */
  outputPtr->chisqParams = (FindChirpChisqParams *)
    LALMalloc( sizeof(FindChirpChisqParams) );
  if ( !outputPtr->chisqParams )
  {
    LALFree( outputPtr );
    *output = NULL;
    ABORT( status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
  }
  memset( outputPtr->chisqParams, 0, sizeof(FindChirpChisqParams) );

  /* create memory for the chisq parameters */
  outputPtr->chisqInput = (FindChirpChisqInput *)
    LALMalloc( sizeof(FindChirpChisqInput) );
  if ( !outputPtr->chisqInput )
  {
    LALFree( outputPtr->chisqParams );
    LALFree( outputPtr );
    *output = NULL;
    ABORT( status, FINDCHIRP_EALOC, FINDCHIRP_MSGEALOC );
  }
  memset( outputPtr->chisqInput, 0, sizeof(FindChirpChisqInput) );
    
  /* don't need this because of the memset above */
  /* outputPtr->chisqParams->chisqBinVec = NULL; */
  /* outputPtr->chisqParams->plan        = NULL; */
  /* outputPtr->chisqInput->qtildeVec = NULL; */
  /* outputPtr->chisqInput->qVec      = NULL; */


  /*
   *
   * create fft plans and workspace vectors
   *
   */


  /* create plan for optimal filter */
  LALEstimateInvComplexFFTPlan( status->statusPtr, 
      &(outputPtr->invPlan), params->numPoints );
  BEGINFAIL( status )
  {
    LALFree( outputPtr->chisqInput );
    LALFree( outputPtr->chisqParams );
    LALFree( outputPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* create plan for chisq filter */
  LALEstimateInvComplexFFTPlan( status->statusPtr, 
      &(outputPtr->chisqParams->plan), params->numPoints );
  BEGINFAIL( status )
  {
    TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
        &(outputPtr->invPlan) ), status );

    LALFree( outputPtr->chisqInput );
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
        &(outputPtr->chisqParams->plan) ), status );
    TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
        &(outputPtr->invPlan) ), status );
    
    LALFree( outputPtr->chisqInput );
    LALFree( outputPtr->chisqParams );
    LALFree( outputPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* create workspace vector for optimal filter: freq domain */
  LALCCreateVector( status->statusPtr, &(outputPtr->qtildeVec), 
      params->numPoints );
  BEGINFAIL( status )
  {
    TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ), status );

    TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
        &(outputPtr->chisqParams->plan) ), status );
    TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
        &(outputPtr->invPlan) ), status );
    
    LALFree( outputPtr->chisqInput );
    LALFree( outputPtr->chisqParams );
    LALFree( outputPtr );
    *output = NULL;
  }
  ENDFAIL( status );

  /* create workspace vector for chisq filter */
  LALCreateVector (status->statusPtr, &(outputPtr->chisqVec), 
      params->numPoints);
  BEGINFAIL( status )
  {
    TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ), status );
    TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ), status );

    TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
        &(outputPtr->chisqParams->plan) ), status );
    TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
        &(outputPtr->invPlan) ), status );
    
    LALFree( outputPtr->chisqInput );
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
    LALCreateVector (status->statusPtr, &(outputPtr->rhosqVec), 
        params->numPoints);
    BEGINFAIL( status )
    {
      TRY( LALDestroyVector (status->statusPtr, &(outputPtr->chisqVec) ), status ); 
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ), status );
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ), status );

      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
          &(outputPtr->chisqParams->plan) ), status );
      TRY( LALDestroyComplexFFTPlan( status->statusPtr, 
          &(outputPtr->invPlan) ), status );

      LALFree( outputPtr->chisqInput );
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

  INITSTATUS( status, "LALFindChirpFilterFinalize", FINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure handle is non-null and points to a non-null pointer */
  ASSERT( output, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( *output, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );


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

  /* destroy plan for chisq filter */
  LALDestroyComplexFFTPlan( status->statusPtr, &(outputPtr->chisqParams->plan) );
  CHECKSTATUSPTR( status );

  /* destroy workspace vector for optimal filter: freq domain */
  LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) );
  CHECKSTATUSPTR( status );

  /* destroy workspace vector for optimal filter: freq domain */
  LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) );
  CHECKSTATUSPTR( status );

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
  
  /* input structure */
  LALFree( outputPtr->chisqInput );


  /*
   *
   * destroy vector to store snrsq, if it exists
   *
   */


  if ( outputPtr->rhosqVec )
  {
    LALDestroyVector( status->statusPtr, &(outputPtr->rhosqVec) );
    CHECKSTATUSPTR( status );
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



void
LALFindChirpFilterSegment (
    LALStatus                  *status,
    InspiralEvent             **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
{
  UINT4                 j, k;
  UINT4                 numPoints;
  REAL4                 rhosq;
  REAL4                 deltaT;
  REAL4                 norm;
  BOOLEAN               haveChisq   = 0;
  COMPLEX8             *qtilde      = NULL;
  COMPLEX8             *q           = NULL;
  COMPLEX8             *inputData   = NULL;
  COMPLEX8             *tmpltSignal = NULL;

  REAL4                 distsq      = 0.0;

  INITSTATUS( status, "LALFindChirpFilter", FINDCHIRPFILTERC );
  ATTATCHSTATUSPTR( status );


  /*
   *
   * check that the arguments are reasonable
   *
   */


  /* make sure the output handle exists, but points to a null pointer */
  ASSERT( eventList, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( !*eventList, status, FINDCHIRP_ENNUL, FINDCHIRP_MSGENNUL );

  /* make sure that the parameter structure exists */
  ASSERT( params, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* check that the filter parameters are reasonable */
  ASSERT( params->deltaT > 0, status,
      FINDCHIRP_EDTZO, FINDCHIRP_MSGEDTZO );
  ASSERT( params->rhosqThresh > 0, status,
      FINDCHIRP_ERHOT, FINDCHIRP_MSGERHOT );
  ASSERT( params->chisqThresh > 0, status,
      FINDCHIRP_ECHIT, FINDCHIRP_MSGECHIT );

  /* check that the fft plan exists */
  ASSERT( params->invPlan, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* check that the workspace vectors exist */
  ASSERT( params->qVec, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( params->qVec->data, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( params->qtildeVec, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( params->qtildeVec->data, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* check that the chisq parameter and input structures exist */
  ASSERT( params->chisqParams, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( params->chisqInput, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec ) 
  {
    ASSERT( params->rhosqVec->data, status, 
        FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  }
  
  /* if a chisqVec vector has been created, check we can store data in it */
  if ( params->chisqVec ) 
  {
    ASSERT( params->chisqVec->data, status,
        FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->tmplt, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( input->fcTmplt, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRP_ENULL, FINDCHIRP_MSGENULL );


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

  /* number of points in a segment */
  numPoints = params->qVec->length;


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
    memset( params->rhosqVec->data, 0, numPoints * sizeof( REAL4 ) );

  /* normalisation */
  norm = 4.0 * (deltaT / (REAL4)numPoints) / input->segment->segNorm;

  for ( j = 0; j < numPoints; ++j )
  {
    REAL4 modqsq = q[j].re * q[j].re + q[j].im * q[j].im;

    /* signal to noise squared */
    rhosq = norm * modqsq;

    /* if full snrsq vector is required, store the snrsq */
    if ( params->rhosqVec ) params->rhosqVec->data[j] = rhosq;

    /* test for rhosq threshold */
    if ( rhosq > params->rhosqThresh )
    {
      /* compute chisqVector */
      if ( ! haveChisq )
      {
        memset( params->chisqVec->data, 0, 
            params->chisqVec->length * sizeof(REAL4) );

        /* pointers to chisq input */
        params->chisqInput->qtildeVec = params->qtildeVec;
        params->chisqInput->qVec      = params->qVec;

        /* pointer to the chisq bin vector in the segment */
        params->chisqParams->chisqBinVec = input->segment->chisqBinVec;
        params->chisqParams->chisqNorm   = sqrt( norm );

        LALFindChirpChisqVeto( status->statusPtr, params->chisqVec, 
            params->chisqInput, params->chisqParams );
        CHECKSTATUSPTR (status); 

        haveChisq = 1;
      }
    }

    /* this needs updating: it is just temporary code for the mdc */
    if ( params->chisqVec->data[j] < params->chisqThresh &&
        rhosq > params->rhosqThresh )
    {
      if ( ! *eventList )
      {
        *eventList = (InspiralEvent *) LALMalloc( sizeof(InspiralEvent) );
        memset( *eventList, 0, sizeof(InspiralEvent) );

        (*eventList)->time.gpsSeconds     = 25;
        (*eventList)->time.gpsNanoSeconds = 0176;
        
        memcpy( &((*eventList)->tmplt), input->tmplt , sizeof(InspiralTemplate) );
      }

      if ( rhosq > (*eventList)->snrsq )
      {
        (*eventList)->snrsq = rhosq;
        (*eventList)->chisq = params->chisqVec->data[j];

        /* effective distance */
        distsq = (input->fcTmplt->tmpltNorm * 
            input->segment->segNorm * input->segment->segNorm) / modqsq;

        (*eventList)->sigma = sqrt( distsq );
      }
    }
  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
