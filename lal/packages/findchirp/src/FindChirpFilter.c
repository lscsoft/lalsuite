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

#if 0
<lalVerbatim file="FindChirpFilterCV">
Author: Brown D. A.
$Id$
</lalVerbatim>

<lalLaTeX>
\input{FindChirpFilterCDoc}

\vfill{\footnotesize\input{FindChirpFilterCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
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
    TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ), 
        status );

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
    TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qtildeVec) ), 
        status );
    TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ), 
        status );

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
      TRY( LALCDestroyVector( status->statusPtr, &(outputPtr->qVec) ), 
          status );

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
  UINT4                 eventId = 0;
  REAL4                 deltaT;
  REAL4                 norm;
  REAL4                 modqsqThresh;
  UINT4                 eventStartIdx = 0;
  REAL4                 chirpTime     = 0;
  BOOLEAN               haveChisq     = 0;
  COMPLEX8             *qtilde        = NULL;
  COMPLEX8             *q             = NULL;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;

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
   * compute viable search regions in the snrsq vector
   *
   */


  {
    /* calculate the length of the chirp */
    REAL4 eta = input->tmplt->eta;
    REAL4 m1 = input->tmplt->mass1;
    REAL4 m2 = input->tmplt->mass2;
    REAL4 fmin = input->segment->fLow;
    REAL4 m = 2 * ( m1 > m2 ? m2 : m1 );
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

#if 0
    fprintf( stdout, "m1 = %e m2 = %e => %e seconds => %d points\n"
        "invSpecTrunc = %d => ignoreIndex = %d\n", 
        m1, m2, chirpTime, deltaEventIndex, 
        input->segment->invSpecTrunc, ignoreIndex );
    fflush( stdout );
#endif

    /* XXX check that we are not filtering corrupted data XXX */
    /* XXX this is hardwired to 1/4 segment length        XXX */
    if ( ignoreIndex > numPoints / 4 )
    {
      ABORT( status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP );
    }
    /* XXX reset ignoreIndex to one quarter of a segment XXX */
    ignoreIndex = numPoints / 4;
  }

#if 0
    fprintf( stdout, "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    fflush( stdout );
#endif

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

  /* normalised threhold */
  modqsqThresh = params->rhosqThresh / norm;

  /* if full snrsq vector is required, store the snrsq */
  if ( params->rhosqVec ) 
  {
    snprintf( params->rhosqVec->name, LALNameLength * sizeof(CHAR),
        "%s:rhosq:output", input->segment->data->name );
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
        params->chisqParams->bankMatch   = input->tmplt->minMatch;

        /* compute the chisq threshold: this is slow! */
        LALFindChirpChisqVeto( status->statusPtr, params->chisqVec, 
            params->chisqInput, params->chisqParams );
        CHECKSTATUSPTR (status); 

        haveChisq = 1;
      }

      /* if we have don't have a chisq or the chisq drops below threshold */
      /* start processing events                                          */
      if ( ! input->segment->chisqBinVec->length ||
          params->chisqVec->data[j] < params->chisqThresh )
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

          /* set the impuse time for the event */
          thisEvent->template_duration = (REAL8) chirpTime;
          
          /* record the ifo and channel name for the event */
          strncpy( thisEvent->ifo, input->segment->data->name, 
              2 * sizeof(CHAR) );
          strncpy( thisEvent->channel, input->segment->data->name + 3,
              (LALNameLength - 3) * sizeof(CHAR) );

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

          /* set snrsq, chisq, sigma and effDist for this event */
          if ( input->segment->chisqBinVec->length )
          {
            thisEvent->chisq     = params->chisqVec->data[timeIndex];
            thisEvent->chisq_dof = input->segment->chisqBinVec->length - 1;
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
          thisEvent->event_duration /= (REAL8) deltaT;
          
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

    /* set the impuse time for the event */
    thisEvent->template_duration = (REAL8) chirpTime;

    /* record the ifo name for the event */
    strncpy( thisEvent->ifo, input->segment->data->name, 
        2 * sizeof(CHAR) );
    strncpy( thisEvent->channel, input->segment->data->name + 3,
        (LALNameLength - 3) * sizeof(CHAR) );

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

    /* set snrsq, chisq, sigma and effDist for this event */
    if ( input->segment->chisqBinVec->length )
    {
      thisEvent->chisq = params->chisqVec->data[timeIndex];
      thisEvent->chisq_dof = input->segment->chisqBinVec->length - 1;
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
    thisEvent->event_duration /= (REAL8) deltaT;
  }    


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
