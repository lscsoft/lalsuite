/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpFilter.c
 *
 * Author: Brown, D. A., Modified by Messaritaki, E.
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
    thisEvent->event_duration /= (REAL8) deltaT;
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
  UINT4                 kendBCV, kTop;        
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
  COMPLEX8Vector       *qtildeVec1    = NULL; 
  COMPLEX8Vector       *qtildeVec2    = NULL; 
  COMPLEX8             *qtilde1       = NULL; 
  COMPLEX8             *qtilde2       = NULL; 
  COMPLEX8             *q             = NULL;
  COMPLEX8Vector       *qVec1         = NULL; 
  COMPLEX8Vector       *qVec2         = NULL; 
  COMPLEX8             *q1            = NULL; 
  COMPLEX8             *q2            = NULL; 
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;
  LALMSTUnitsAndAcc     gmstUnits;
  REAL4                 normFFT;              
  REAL4                 b1;                  
  REAL4                 a2;                  
  REAL4                 templateNorm;

  FindChirpChisqInput  *chisqInput1;
  FindChirpChisqInput  *chisqInput2;

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
  ASSERT(params->qVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qVec->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT(params->qtildeVec->data,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);
  
ASSERT( qVec1, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( qVec1->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( qtildeVec1, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( qtildeVec1->data,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);
  ASSERT( qVec2, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( qVec2->data, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( qtildeVec2, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( qtildeVec2->data,status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL);


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

  q1 = qVec1->data; 
  q2 = qVec2->data; 
  qtilde1 = qtildeVec1->data; 
  qtilde2 = qtildeVec2->data; 
  
  
  /* template and data */
  inputData    = input->segment->data->data->data;
  tmpltSignal  = input->fcTmplt->data->data;
  templateNorm = input->fcTmplt->tmpltNorm;
  deltaT       = params->deltaT;

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
    /* length of the chirp:                      */
    /* For BCV, the chirpTime is not calculated  */
    /* but set to 1.5 seconds.                   */
    /* chirpTime is not recorded                 */
    /* just used to maximize over chirp.         */
#if 0
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
    REAL4 chirpTime = c0*(1 + c2*x2 + c3*x3 + c4*x4)/x8; 
#endif
    REAL4 m1 = input->tmplt->mass1;
    REAL4 m2 = input->tmplt->mass2;
    REAL4 m = m1 + m2;
    REAL4 fmin = input->segment->fLow;
    REAL4 x  = pow(LAL_PI*m*LAL_MTSUN_SI*fmin, 1.0/3.0);
    REAL4 chirpTime = 1.5;

    /* template parameters */
    REAL4 psi0 = input->tmplt->psi0;                        
    REAL4 psi3 = input->tmplt->psi3;                        
    REAL4 fendBCV = input->tmplt->fendBCV;                  
    /* k that corresponds to fendBCV  */
    UINT4 kendBCV = floor( numPoints * deltaT * fendBCV );  
    /* BCV normalization */
    REAL4 b1 = input->segment->b1;  
    REAL4 a2 = input->segment->a2; 

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
   * compute qtilde, qtilde1, qtilde2 and q, q1, q2 
   *
   */


  memset( qtilde, 0, numPoints * sizeof(COMPLEX8) );
  memset( qtilde1, 0, numPoints * sizeof(COMPLEX8) );
  memset( qtilde2, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  kTop = ( (numPoints/2) < kendBCV ? numPoints/2 : kendBCV );
  for ( k = 1; k < kTop; ++k )
  {
    REAL4 r = inputData[k].re;
    REAL4 s = 0 - inputData[k].im;   /* note complex conjugate */
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = tmpltSignal[k].im;     

    qtilde[k].re = r*x - s*y;
    qtilde[k].im = r*y + s*x;
    qtilde1[k].re = qtilde[k].re; 
    qtilde1[k].im = qtilde[k].im; 
    qtilde2[k].re = qtilde1[k].re * pow(k, (2/3)); 
    qtilde2[k].im = qtilde1[k].im * pow(k, (2/3)); 
  }

  /* qtilde negative frequency only: not DC or nyquist */
  if ( params->computeNegFreq )
  {
    kTop = ( (numPoints-1) < kendBCV ? (numPoints-1) : kendBCV );	  
    for ( k = numPoints/2 + 2; k < numPoints - 1; ++k )
    {
      REAL4 r = inputData[k].re;
      REAL4 s = 0 - inputData[k].im; /* note complex conjugate */
      REAL4 x = tmpltSignal[k].re;
      REAL4 y = tmpltSignal[k].im;   

      qtilde[k].re = r*x - s*y;
      qtilde[k].im = r*y + s*x;
      qtilde1[k].re = qtilde[k].re; 
      qtilde1[k].im = qtilde[k].im; 
      qtilde2[k].re = qtilde1[k].re * pow(k, (2/3));
      qtilde2[k].im = qtilde1[k].im * pow(k, (2/3)); 

    }
   }

   /* inverse fft to get q, q1 and q2 */
   LALCOMPLEX8VectorFFT( status->statusPtr, params->qVec, params->qtildeVec,
       params->invPlan ); 
   CHECKSTATUSPTR( status );
   LALCOMPLEX8VectorFFT( status->statusPtr, qVec1, qtildeVec1,
       params->invPlan ); 
   CHECKSTATUSPTR( status ); 
   LALCOMPLEX8VectorFFT( status->statusPtr, qVec2, qtildeVec2,
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
  /*  For the BCV templates, templateNorm is equal to                    */
  /*  (5*m/96) * (M/pi^2)^(2/3) (Tsun/Dt)^(-1/3) (2*Msun(L)*d/1Mpc)^2    */

  params->norm = norm =
     (deltaT / (REAL4)numPoints) * (deltaT / (REAL4)numPoints) * templateNorm;

  /* normalised snr threhold */
  modqsqThresh = params->rhosqThresh / norm ;

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
      REAL4 modq1sq = 
	      ( 2.0 * b1 * q1[j].re - 2.0 * a2 * q2[j].re ) * 
	      ( 2.0 * b1 * q1[j].re - 2.0 * a2 * q2[j].re )    
	    + ( 2.0 * b1 * q1[j].im - 2.0 * a2 * q2[j].im ) * 
	      ( 2.0 * b1 * q1[j].im - 2.0 * a2 * q2[j].im );    
      REAL4 modq2sq = 
              ( 2.0 * b1 * q1[j].re + 2.0 * a2 * q2[j].re ) *
              ( 2.0 * b1 * q1[j].re + 2.0 * a2 * q2[j].re )      
            + ( 2.0 * b1 * q1[j].im + 2.0 * a2 * q2[j].im ) *
              ( 2.0 * b1 * q1[j].im + 2.0 * a2 * q2[j].im );   
      REAL4 modqsq = modq1sq + modq2sq + 
	      2.0 * sqrt( modq1sq ) * sqrt( modq2sq ) ;       
      params->rhosqVec->data->data[j] = norm * modqsq;   
    }
  }


  /* look for an events in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
  {
      REAL4 modq1sq = 
	      ( 2.0 * b1 * q1[j].re - 2.0 * a2 * q2[j].re ) *
	      ( 2.0 * b1 * q1[j].re - 2.0 * a2 * q2[j].re )  
	    + ( 2.0 * b1 * q1[j].im - 2.0 * a2 * q2[j].im ) *
	      ( 2.0 * b1 * q1[j].im - 2.0 * a2 * q2[j].im );    
      REAL4 modq2sq =      
              ( 2.0 * b1 * q1[j].re + 2.0 * a2 * q2[j].re ) *
              ( 2.0 * b1 * q1[j].re + 2.0 * a2 * q2[j].re )         
            + ( 2.0 * b1 * q1[j].im + 2.0 * a2 * q2[j].im ) *
              ( 2.0 * b1 * q1[j].im + 2.0 * a2 * q2[j].im );   
      REAL4 modqsq = modq1sq + modq2sq + 
              2.0 * sqrt( modq1sq ) * sqrt( modq2sq ) ;       


    /* if snrsq exceeds threshold at any point */
    if ( modqsq > modqsqThresh )                  
    {
      /* compute chisq vector if it does not exist and we want it */
      if ( ! haveChisq  && input->segment->chisqBinVec->length )
      {
        memset( params->chisqVec->data, 0,
            params->chisqVec->length * sizeof(REAL4) );

        /* pointers to chisq input */
#if 0
        params->chisqInput->qtildeVec = params->qtildeVec;
        params->chisqInput->qVec      = params->qVec;
#endif
	chisqInput1->qtildeVec = qtildeVec1;
	chisqInput1->qVec = qVec1;
	chisqInput2->qtildeVec = qtildeVec2;
	chisqInput2->qVec = qVec2;
	
        /* pointer to the chisq bin vector in the segment */
        params->chisqParams->chisqBinVec = input->segment->chisqBinVec;
        params->chisqParams->norm        = norm;
	params->chisqParams->b1          = b1 ;
	params->chisqParams->a2          = a2 ;
#if 0
        params->chisqParams->bankMatch   = input->tmplt->minMatch;
#endif

        /* compute the chisq threshold: this is slow! */
        LALFindChirpBCVChisqVeto( status->statusPtr, params->chisqVec,
            chisqInput1, chisqInput2, params->chisqParams );
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
           /* thisEvent->template_duration = (REAL8) chirpTime; */

	   /* record the ifo and channel name for the event */
	   strncpy( thisEvent->ifo, input->segment->data->name,
	       2 * sizeof(CHAR) );
	   strncpy( thisEvent->channel, input->segment->data->name + 3,
	       (LALNameLength - 3) * sizeof(CHAR) );
	   thisEvent->impulse_time = thisEvent->end_time;

           /* copy the template into the event */

#if 0      
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
#endif	   
	   /* temporarily */
	   thisEvent->tau0   = (REAL4) input->tmplt->psi0; 
	   thisEvent->tau3   = (REAL4) input->tmplt->psi3;
	   thisEvent->mchirp = LAL_1_PI *
		   pow( 3.0 / 128.0 / input->tmplt->psi0 , 3.0/5.0 );

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
           thisEvent->event_duration= (REAL8) timeIndex - (REAL8) eventStartIdx;
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
    LALGPStoGMST1( status->statusPtr, &(thisEvent->end_time_gmst),
        &(thisEvent->end_time), &gmstUnits );
    CHECKSTATUSPTR( status );

    /* set the impuse time for the event */
    /* thisEvent->template_duration = (REAL8) chirpTime; */

    /* record the ifo name for the event */
    strncpy( thisEvent->ifo, input->segment->data->name,
        2 * sizeof(CHAR) );
    strncpy( thisEvent->channel, input->segment->data->name + 3,
        (LALNameLength - 3) * sizeof(CHAR) );
    thisEvent->impulse_time = thisEvent->end_time;

    /* copy the template into the event */

#if 0    
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
#endif    
    /* temporarily */
    thisEvent->tau0   = (REAL4) input->tmplt->psi0;   
    thisEvent->tau3   = (REAL4) input->tmplt->psi3;  
    thisEvent->mchirp = LAL_1_PI *
	    pow( 3.0 / 128.0 / input->tmplt->psi0, 3.0/5.0 );


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
      thisEvent->event_duration /= (REAL8) deltaT;
  }


   /* normal exit */
   DETATCHSTATUSPTR( status );
   RETURN( status );
}
