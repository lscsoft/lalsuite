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


void LALFindChirpTDTFilterSegment (
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
   BOOLEAN	haveChisq = 0;
   REAL4	deltaT;
   REAL4	modqsqThresh;
   REAL4	lowerThresh;
   REAL4   	norm;

   COMPLEX8	*qtilde=NULL;
   COMPLEX8	*q=NULL;
   COMPLEX8	*inputData=NULL;
   COMPLEX8	*tmpltSignal=NULL;
   
   SnglInspiralTable	*thisEvent=NULL;
   LALMSTUnitsAndAcc	gmstUnits;

   INITSTATUS(status, "LALFindChirpTDFilter", FINDCHIRPTDFILTERC);
   ATTACHSTATUSPTR(status);

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
  ASSERT( input->fcTmplt->approximant == TaylorF2, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == TaylorF2, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  */
   

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      make sure that approximant is known and PN order is reasonable some where 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   */
    
   /*   assign local pointers */

   q = params->qVec->data;
   qtilde = params->qtildeVec->data;

   /*  template and data */

   inputData = input->segment->data->data->data;
   tmpltSignal = input->fcTmplt->data->data;
   deltaT = params->deltaT;

   /* number of bins */

   numChisqBins = input->segment->chisqBinVec->length -1;

   /*  number of points in a segment */

   numPoints = params->qVec->length;

   /* set the gmst units  */

   gmstUnits.units = MST_HRS;
   gmstUnits.accuracy = LALLEAPSEC_STRICT;

   /* ignore corrupted data at start and end */

   deltaEventIndex  = (UINT4) floor(input->tmplt->tC/deltaT + 1.0);
   ignoreIndex = input->segment->invSpecTrunc/2 + deltaEventIndex;

   CHAR infomsg[256];

   REAL4 m1 = input->tmplt->mass1;
   REAL4 m2 = input->tmplt->mass2;
   REAL4 tmlptDuration = input->tmplt->tC;

   LALSnprintf( infomsg, sizeof(infomsg) / sizeof(*infomsg),
          "m1 = %e m2 = %e => %e seconds => %d points\n"
          "invSpecTrunc = %d => ignoreIndex = %d\n", 
          m1, m2, tmpltDuration, deltaEventIndex, 
          input->segment->invSpecTrunc, ignoreIndex );
   LALInfo( status, infomsg );

   if(ignoreIndex >numPoints/4){
       ABORT(status, FINDCHIRPH_ECRUP, FINDCHIRPH_MSGECRUP);
   }
   ignoreIndex = numPoints/4;

   /*  compute qtilde and q  */

   memset(qtilde, 0, numPoints*sizeof(COMPLEX8));

   REAL4
    
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
      
   norm = 1.0/input->fcTmplt->tmpltNorm;   
   /*******!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        I NEED TO APPLY PROPER NORMALIZATION  AND MODIFY CHI^2 FACTOR
		ACCORDING TO UWM DEFINITION

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   **************/

   if ( params->rhosqVec ){
       memcpy( params->rhosqVec->name, input->segment->data->name,
              LALNameLength * sizeof(CHAR) );
       memcpy( &(params->rhosqVec->epoch), &(input->segment->data->epoch), 
              sizeof(LIGOTimeGPS) );
       params->rhosqVec->deltaT = input->segment->deltaT;
       for ( j = 0; j < numPoints; ++j ){
          REAL4 modqsq = q[j].re * q[j].re + q[j].im * q[j].im;
          params->rhosqVec->data->data[j] = norm * modqsq;
       }
   }

   /* look for an event in the filter output using double threshold method */ 
   

   UINT4 calm = 1;
   UINT4 max_cross = 0;
   UINT4 bin = 0;
   UINT4 startbin=0;
   modqsqThresh = params->rhosqThresh;
   lowerThresh  = modqsqThresh/16.0;
   REAL4 snr = 0.0;

   for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j ){
       REAL4 z = (q[j].re*q[j].re + q[j].im*q[j].im)*norm;
       if(z >= modqsqThresh){
           startbin = j;
           max_cross = 1;        /* we found event */
	   calm = 0;
	   if(z > snr){
	      snr = z;
	      bin = j
	   }
       }else{
           max_cross = 0;
	   if(z <= lowerThresh){
	       if(!calm){

	       /** we found an event at bin/samplingRate 
	           with SNR^2 = snr. We need to compute chi^2
		   and if event passes chi^2 threshold add it to event
		   list 
	        */
		    if( input->segment->chisqBinVec->length ){
			memset(params->chisqVec->data, 0,
				params->chisqVec->length * sizeof(REAL4));
			params->chisqInput->qtildeVec = params->qtildeVec;
			params-chisqInput->qVec	= params->qVEc;
			params->chisqParams->chisqBinVec = input->segment->chisqBinVec;
			params->chisqParams->norm = norm;

		/*	LALFindChirpTDTChisqVeto(status->statusPtr,
			   params->chisqVec, params->chisqInput, params->chisqParams);
			CHECKSTATUSPTR(status); */
		    }   
		    if( !input->segment->chisqBinVec->length 
		    /* || !!!!! here I need to set threshold on chi^2 !!!!!*/){
		       if( ! *eventList){
			  /* start the list if it is first event */
			  thisEvent = *eventList = (SndlInspiralTable *)
			  	LALCalloc(1, sizeof(SnglInspiralTable));
			  if ( ! thisEvent ){
			        ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
          		  }
		       }
		       SnglInspiralTable *lastEvent;
		       INT8 timeNS;

		       timeNS =  1000000000L * (INT8) (input->segment->data->epoch.gpsSeconds);
		       timeNS += (INT8) (input->segment->data->epoch.gpsNanoSeconds);
		       timeNS += (INT8) (1e9 * bin * deltaT);
		       thisEvent->impulse_time.gpsSeconds = (INT4) (timeNS/1000000000L);
          	       thisEvent->impulse_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
		       thisEvent->template_duration = (REAL8)input->tmplt->tC;
		       timeNS += (INT8) (1e9 * input->tmplt-tC);
		       thisEvent->end_time.gpsSeconds = (INT4) (timeNS/1000000000L);
          	       thisEvent->end_time.gpsNanoSeconds = (INT4) (timeNS%1000000000L);
		       LALGPStoGMST1( status->statusPtr, &(thisEvent->end_time_gmst),
              			&(thisEvent->end_time), &gmstUnits );
		       CHECKSTATUSPTR( status );
		        /* record the coalescence phase of the chirp */
		       if ( q[timeIndex].re == 0 ){
		         if ( q[timeIndex].im >= 0 ){
		            thisEvent->coa_phase = LAL_PI / 2.0;
            	         }else{
		            thisEvent->coa_phase = - LAL_PI / 2.0;
            	         }
          	       }else{
		         thisEvent->coa_phase = (REAL4)atan( q[timeIndex].im / q[timeIndex].re );
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
		       if(input->segment->chisqBinVec->length){
        		  /*  !!!! record chi^2 HERE !!!!!!!!!!!!!!!!! */
		       }else{
		          thisEvent->chisq = 0.0;
		          thisEvent->chisq_dof = 0;
		       }
		       thisEvent->snr = sqrt(snr);
		       thisEvent->event_duration = ((REAL8)(bin - startbin))*deltaT;
		       /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		           I need to compute sigmasq, edd_distance
			  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			*/
		       lastEvent = thisEvent;
		       lastEvent->next = thisEvent = (SnglInspiralTable *)
		       		LALCalloc(1, sizeof(SnglInspiralTable));
		       if ( ! lastEvent->next ){
		            ABORT( status, FINDCHIRPH_EALOC, FINDCHIRPH_MSGEALOC );
          	       }

	               /* stick minimal data into the event */
	              thisEvent->end_time.gpsSeconds =bin;
   	              thisEvent->snr = snr;
		     }

		    }
	            calm = 1;
		    snr = 0.0;
		    bin =0;
		    startbin=0;
	       }
	   }
       }
      /* normal exit */
     DETATCHSTATUSPTR( status );
     RETURN( status );
    
   }
         
    /*  max_cross is not used currently, but might be usefull if event 
      spans across segments. We need to take this into account. We might miss
      event if lower threshold is not crossed
     */

    


   

