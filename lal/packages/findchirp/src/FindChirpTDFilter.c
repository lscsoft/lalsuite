/* -----------------------------------
 *     File Name: FindChirpTDFilter.c
 *     Author: S. Babak
 *------------------------------------
 */

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirpTDTemplate.h>

NRCSID (FINDCHIRPTDFILTERC, "$Id$");


void LALFindChirpTDFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
{
   UINT4	j,k;
   UINT4	numPoints;
   UINT4	numChisqBins;
   UINT4	deltaEventIndex;
   UINT4	ignoreIndex;
/*   BOOLEAN	haveChisq = 0;*/
   REAL4	deltaT;
   REAL4	modqsqThresh;
   REAL4	lowerThresh;
   REAL4   	norm;
   REAL4 	chirpTime;
   REAL4 	mismatch;
   REAL4	chisqThreshFac;
 
   REAL4        distNorm;
  const REAL4   cannonDist = 1.0; /* Mpc */
   REAL4 	tmpltNorm;
   REAL4        segNorm;

   REAL4	m1;
   REAL4	m2;
   REAL4 	mu;
   REAL4 	m;
   REAL4 	eta;

   BOOLEAN	haveChisq=0;
   UINT4	eventStartIdx=0;

   COMPLEX8	*qtilde=NULL;
   COMPLEX8	*q=NULL;
   COMPLEX8	*inputData=NULL;
   COMPLEX8	*tmpltSignal=NULL;
   
   SnglInspiralTable	*thisEvent=NULL;
   LALMSTUnitsAndAcc	gmstUnits;

   INITSTATUS(status, "LALFindChirpTDFilter", FINDCHIRPTDFILTERC);
   ATTATCHSTATUSPTR(status);

   /*   Check arguments    */

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

  /* make sure that the template and the segment are both stationary phase 
  ASSERT( input->segment->approximant == TaylorF2, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  */

  ASSERT( input->fcTmplt->approximant == TaylorT1 || input->fcTmplt->approximant == TaylorT3
        || input->fcTmplt->approximant == PadeT1 || input->fcTmplt->approximant == EOB,
      status, FINDCHIRPTDT_ETAPX, FINDCHIRPTDT_MSGETAPX );

  ASSERT( input->fcTmplt->approximant == input->segment->approximant, status,
      FINDCHIRPTDT_ETAPX, FINDCHIRPTDT_MSGETAPX );
   

    
   /*   assign local pointers */

   q = params->qVec->data;
   qtilde = params->qtildeVec->data;

   /*  template and data */

   inputData = input->segment->data->data->data;
   tmpltSignal = input->fcTmplt->data->data;
   deltaT = params->deltaT;

   /* number of bins */

   numChisqBins = input->segment->chisqBinVec->length ?
            input->segment->chisqBinVec->length - 1 : 0;

   /*  number of points in a segment */

   numPoints = params->qVec->length;

   /* set the gmst units  */

   gmstUnits.units = MST_HRS;
   gmstUnits.accuracy = LALLEAPSEC_STRICT;


   chirpTime = input->tmplt->tC;
   if ( chirpTime <= 0.0 )
    {
      ABORT( status, FINDCHIRPH_ECHTZ, FINDCHIRPH_MSGECHTZ );
    }



   /* ignore corrupted data at start and end */

   deltaEventIndex  = (UINT4) floor(chirpTime/deltaT + 1.0);
   ignoreIndex = input->segment->invSpecTrunc/2 + deltaEventIndex;

   m1 = input->tmplt->mass1;
   m2 = input->tmplt->mass2;
   if(lalDebugLevel & LALINFO){
      CHAR infomsg[256];


       LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
          "m1 = %e m2 = %e => %e seconds => %d points\n"
          "invSpecTrunc = %d => ignoreIndex = %d\n", 
          m1, m2, chirpTime, deltaEventIndex, 
          input->segment->invSpecTrunc, ignoreIndex );
       LALInfo( status, infomsg );
   }

   if(ignoreIndex >numPoints/4){
       ABORT(status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP);
   }
   ignoreIndex = numPoints/4;

   if ( lalDebugLevel & LALINFO )
  {
    CHAR infomsg[256];

    LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, infomsg );
  }



   /*  compute qtilde and q  */

   memset(qtilde, 0, numPoints*sizeof(COMPLEX8));

    
   for(k=1; k < numPoints/2; ++k){
       REAL4 r = inputData[k].re;
       REAL4 s = inputData[k].im;
       REAL4 x = tmpltSignal[k].im;
       REAL4 y = 0 - tmpltSignal[k].im; /* note complex conjugate */
 
       qtilde[k].re = r*x - s*y;
       qtilde[k].im = r*y + s*x;
   }

   /* not implemented for negative freuencies */

   /*  inverse fft to get q */
       
   LALCOMPLEX8VectorFFT(status->statusPtr, params->qVec, params->qtildeVec, \
   			params->invPlan);
   CHECKSTATUSPTR(status);
   
   /* compute SNR */

   if(params->rhosqVec)
      memset(params->rhosqVec->data->data, 0, numPoints*sizeof(REAL4));
      
      /* segNorm is stored in fcTmplt->tmpltNorm */
   params->norm = norm = 4.0*(deltaT/(REAL4)numPoints)/input->fcTmplt->tmpltNorm;   

      /* normalised snr threhold */
   modqsqThresh = params->rhosqThresh / norm;

  /*  preparation (modifying amplitude of chi^2) for \chi^2 */
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


  /* We need to compute tmpltNorm */
  /* template dependent normalisation */
 /*------------------------------------------------*/
      mu = (REAL4) input->tmplt->mu;
      m = (REAL4) input->tmplt->totalMass;
      eta = (REAL4) input->tmplt->eta;
      distNorm = 2.0 * LAL_MRSUN_SI / (cannonDist * 1.0e6 * LAL_PC_SI);
      distNorm *= input->segment->segNorm->data[1]; /* multiplication with dynRange */
      tmpltNorm = sqrt( (5.0*mu) / 96.0 ) *
                  pow( m / (LAL_PI*LAL_PI) , 1.0/3.0 ) *
                  pow( LAL_MTSUN_SI / deltaT, -1.0/6.0 );

      tmpltNorm *= tmpltNorm;
      tmpltNorm *= distNorm * distNorm;
 /*****************************************************/
      segNorm = input->fcTmplt->tmpltNorm;


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

        /* compute the chisq bin boundaries for this template */
        if ( ! params->chisqParams->chisqBinVec->data )
        {
           /*  I will implement chi^2 later */
      /*    LALFindChirpComputeChisqBins( status->statusPtr,
              params->chisqParams->chisqBinVec, input->segment, kmax );
          CHECKSTATUSPTR( status ); */
        }

        /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             CHI^2 is not implemented yet
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

        /* compute the chisq threshold: this is slow! */
/*        LALFindChirpChisqVeto( status->statusPtr, params->chisqVec,
            params->chisqInput, params->chisqParams );
          CHECKSTATUSPTR (status);

        haveChisq = 1;
*/
    }


    if ( ! input->segment->chisqBinVec->length ||
          params->chisqVec->data[j] <
          (params->chisqThresh * ( 1.0 + modqsq * chisqThreshFac )) )
     {
        if ( ! *eventList )
        {
          /* store the stttt of the crossing */
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
            /* so we multiply r^2 by p here to get chisq                 
            thisEvent->chisq =
              params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
            thisEvent->chisq_dof = numChisqBins;
            */
          }
          else
          {
            thisEvent->chisq     = 0;
            thisEvent->chisq_dof = 0;
          }



          thisEvent->sigmasq = norm * segNorm *segNorm* tmpltNorm;
          thisEvent->eff_distance = (tmpltNorm * segNorm *segNorm) / thisEvent->snr;
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
            /* so we multiply r^2 by p here to get chisq                 
            thisEvent->chisq =
              params->chisqVec->data[timeIndex] * (REAL4) numChisqBins;
            thisEvent->chisq_dof = numChisqBins;
          */
      }
    else
    {
      thisEvent->chisq     = 0;
      thisEvent->chisq_dof = 0;
    }
    thisEvent->sigmasq = norm * segNorm * segNorm * tmpltNorm;
    thisEvent->eff_distance = (tmpltNorm * segNorm * segNorm ) / thisEvent->snr;
    thisEvent->eff_distance = sqrt( thisEvent->eff_distance );

    thisEvent->snr *= norm;
    thisEvent->snr = sqrt( thisEvent->snr );

     /* compute the time since the snr crossing */
    thisEvent->event_duration = (REAL8) timeIndex - (REAL8) eventStartIdx;
    thisEvent->event_duration *= (REAL8) deltaT;

  }

/*
   *
   * check the events pass the filter output veto
   *
   *


  if ( params->filterOutputVetoParams )
  {
    LALFindChirpFilterOutputVeto( status->statusPtr, eventList, params->qVec,
        norm, params->filterOutputVetoParams );
    CHECKSTATUSPTR( status );
  }
*/

  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
















    


   

