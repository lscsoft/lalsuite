#include <lal/CoherentEstimation.h>
#include <lal/DetResponse.h>
#include <lal/TimeDelay.h>
#include <lal/AVFactories.h>
#include <strings.h>
#include <math.h>

#define EPS 1.0e-7

NRCSID( COHERENTESTIMATIONC, "$Id$" );



void
LALCoherentEstimation ( LALStatus          *stat,
			COMPLEX8TimeSeries *output,
			CoherentEstimation *params,
			DetectorsData      *in ) {
  /* 
     NOTES:
     o destroys input (in)
     o order of in must be same as order of params->filters
     o output time wrt center of Earth
     o output.re = plus, output.im = cross
  */

  INT4 i, j, /* counters */
    iPad, ePad, /* used for padding */
    del;  /* integer time delay */
  REAL4 y; /* dummy */
  REAL8 p1,p2; /* weights in interpolation */
  REAL8 *tDelays; /* time delays wrt Earth center */
  REAL8 mDelay; /* max of time delays */
  REAL8 tmid; /* mid time point */
  LALDetAMResponse *F; /* response functions */
  DetTimeAndASource dtS; /* for time delays */
  LALPlaceAndGPS pGPS; /* for time delays */
  LALDetAndSource dAs; /* for responses */
  LALSource source; /* for responses */
  LIGOTimeGPS t0; /* start time of data */
  REAL8 normPlus, normCross; /* for normalization */

  /***********************************************************************/
  /* initialize status & validate input                                  */
  /***********************************************************************/
  INITSTATUS( stat, "LALCoherentEstimation", COHERENTESTIMATIONC);
  ATTATCHSTATUSPTR( stat );

  ASSERT ( in, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( in->data, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( in->Ndetectors > 0 && in->Ndetectors < 10, stat, COHERENTESTIMATIONH_E0DEC, COHERENTESTIMATIONH_MSGE0DEC );

  for(i=0;i<in->Ndetectors;i++) {
    ASSERT ( in->data[i].data, stat, COHERENTESTIMATIONH_E0DEC, COHERENTESTIMATIONH_MSGE0DEC );
  }

  t0 = in->data[0].epoch;

  for(i=1;i<in->Ndetectors;i++) {
    ASSERT ( in->data[i].epoch.gpsSeconds == t0.gpsSeconds &&
	     in->data[i].epoch.gpsNanoSeconds == t0.gpsNanoSeconds, stat, COHERENTESTIMATIONH_EDST, COHERENTESTIMATIONH_MSGEDST );
  }

  ASSERT ( params, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( params->detectors, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( params->filters, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( params->position, stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  ASSERT ( params->Ndetectors == in->Ndetectors, stat, COHERENTESTIMATIONH_EICE, COHERENTESTIMATIONH_MSGEICE );

  for(i=0;i<in->Ndetectors;i++) {
    ASSERT ( params->filters[i], stat, COHERENTESTIMATIONH_EICE, COHERENTESTIMATIONH_MSGEICE );
  }

  /***********************************************************************/
  /* if params->preProcessed is 0, need to pre-processed by applying the */
  /* IIR filters twice                                                   */
  /***********************************************************************/
  if(!(params->preProcessed) && params->nPreProcessed) {
  
    REAL8Vector *tmpR8 = NULL;
    UINT4 jind;

    /* loop over all detectors, and filter twice */
    for(i=0; i<params->Ndetectors; i++) {

      TRY ( LALDCreateVector( stat->statusPtr,
			      &tmpR8,
			      in->data[i].data->length ), stat );

      for(jind = 0; jind<in->data[i].data->length; jind++) {
	tmpR8->data[jind] = (REAL8)(in->data[i].data->data[jind]);
      }
     
      for(j=0; j<params->nPreProcessed; j++) { 
        TRY ( LALIIRFilterREAL8Vector( stat->statusPtr,
	       			       tmpR8,
				       params->filters[i]), stat );
      }

      for(jind = 0; jind<in->data[i].data->length; jind++) {
	in->data[i].data->data[jind] = (REAL4)(tmpR8->data[jind]);
      }

      TRY ( LALDDestroyVector (  stat->statusPtr,
				 &tmpR8 ) , stat );

    }

    params->preProcessed = 1;

  }


  /***********************************************************************/
  /* Compute time delays and response functions                          */
  /***********************************************************************/
  if(!(tDelays = (REAL8 *)LALMalloc(params->Ndetectors * sizeof(REAL8)))) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }

  dtS.p_source = params->position;

  /* delays are computed wrt to center of data stretch */
  dtS.p_det_and_time = &pGPS;
  pGPS.p_gps = &t0;
  tmid = 0.5 * (REAL8)(in->data[0].data->length) * in->data[0].deltaT;
  pGPS.p_gps->gpsSeconds += (INT4)floor(tmid);
  pGPS.p_gps->gpsNanoSeconds += (INT4)floor(1E9*(tmid-floor(tmid)));
  if(pGPS.p_gps->gpsNanoSeconds >= 1000000000) {
    pGPS.p_gps->gpsSeconds += (INT4)floor((REAL8)(pGPS.p_gps->gpsNanoSeconds) / 1E9);
    pGPS.p_gps->gpsNanoSeconds -= 1000000000 * (INT4)floor((REAL8)(pGPS.p_gps->gpsNanoSeconds) / 1E9);
  }


  if(!(F = (LALDetAMResponse *)LALMalloc(params->Ndetectors * sizeof(LALDetAMResponse)))) {
    ABORT ( stat, COHERENTESTIMATIONH_EMEM, COHERENTESTIMATIONH_MSGEMEM );
  }

  dAs.pSource = &source;
  source.equatorialCoords = *(params->position);
  source.orientation = params->polAngle;

  for(i=0; i<params->Ndetectors; i++) {

    pGPS.p_detector = dAs.pDetector = params->detectors + i;

    /* tDelays = arrival time at detector - arrival time a center of Earth */
    TRY ( LALTimeDelayFromEarthCenter( stat->statusPtr, 
				       tDelays + i, 
				       &dtS ), stat );


    if(isnan(tDelays[i])) {
      ABORT ( stat, COHERENTESTIMATIONH_ENUM, COHERENTESTIMATIONH_MSGENUM );
    }

    TRY ( LALComputeDetAMResponse ( stat->statusPtr, 
				    F + i,
				    &dAs,
				    pGPS.p_gps ), stat );

    if(isnan(F[i].cross) || isnan(F[i].plus)) {
      ABORT ( stat, COHERENTESTIMATIONH_ENUM, COHERENTESTIMATIONH_MSGENUM );
    }
  }

  /***********************************************************************/
  /* Compute and store estimated data                                    */
  /***********************************************************************/
  /* update output parameters */
  output->epoch = in->data[0].epoch;
  output->deltaT = in->data[0].deltaT;
  output->f0 = in->data[0].f0;
  memcpy(&(output->sampleUnits), &(in->data[0].sampleUnits), sizeof(LALUnit));

  /* make sure output is zero */
  bzero(output->data->data, output->data->length * sizeof(COMPLEX8));

  /* set time origine on detector with largest delay */
  mDelay = tDelays[0];
  for(i=1; i<params->Ndetectors; i++) {
	if(tDelays[i] > mDelay) {
		mDelay = tDelays[i];
	}
  }

  for(i=0; i<params->Ndetectors; i++) {
    tDelays[i] -= mDelay;
  }

  /* response function */
  for(i=0; i<params->Ndetectors; i++) {

    F[i].cross = F[i].cross / (fabs(F[i].cross) * (params->across + pow(fabs(F[i].cross), params->bcross)));
    F[i].plus = F[i].plus / (fabs(F[i].plus) * ( params->aplus + pow(fabs(F[i].plus), params->bplus) ) );

  }


  normPlus = normCross = 0.0;
  for(i=0; i<params->Ndetectors; i++) {
    normCross += fabs(F[i].cross);
    normPlus += F[i].plus*F[i].plus;
  }

  normPlus = sqrt(normPlus);


  for(i=0; i<params->Ndetectors; i++) {
    F[i].plus /= normPlus;
    F[i].cross /= normCross;
  }

  /* loop */
  for(i=0; i<params->Ndetectors; i++) {

    /* setup padding and weights */
    if(tDelays[i] < 0.0) { 
      /* need padding at beginning */
      iPad = (INT4)floor(-tDelays[i]/output->deltaT);
      ePad = 0;

      /* set integer delay (for p1 weight) */
      del = -iPad;

      /* set weights */
      p1 = ceil(tDelays[i] / output->deltaT) - tDelays[i] / output->deltaT;
      p2 = 1.0 - p1;
    } else { 
      /* need padding at end */ 
      iPad = 0;
      ePad = (INT4)ceil(tDelays[i]/output->deltaT);

      /* integer delay */
      del = ePad;

      /* weights */
      p1 = ceil(tDelays[i] / output->deltaT) - tDelays[i] / output->deltaT;
      p2 = 1.0 - p1;
    }


    /* interpolate using time delays */
    for(j=iPad+1; j<output->data->length - ePad; j++) {

      y = p1 * in->data[i].data->data[del+j-1] + p2 * in->data[i].data->data[del+j];

      output->data->data[j].re += y * F[i].plus;
      output->data->data[j].im += y * F[i].cross;

    }

  }

  /***********************************************************************/
  /* clean up and normal return                                          */
  /***********************************************************************/
  LALFree(tDelays);
  LALFree(F);

  DETATCHSTATUSPTR( stat );
  RETURN( stat );

}


void 
LALClearCoherentData (
		      LALStatus     *stat,
		      DetectorsData *dat
		      ) {

  UINT4 i;
  
  INITSTATUS( stat, "LALClearCoherentData", COHERENTESTIMATIONC);
  ATTATCHSTATUSPTR( stat );

  if(!dat) {
    ABORT ( stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  }

  for(i=0; i<dat->Ndetectors; i++) {
    if(dat->data[i].data) {
      TRY ( LALDestroyVector(stat->statusPtr, &(dat->data[i].data)), stat );
    }
  }

  LALFree(dat->data);

  DETATCHSTATUSPTR( stat );
  RETURN( stat );

}



void 
LALClearCoherentInfo (
		      LALStatus     *stat,
		      CoherentEstimation *dat
		      ) {

  UINT4 i;
  
  INITSTATUS( stat, "LALClearCoherentInfo", COHERENTESTIMATIONC);
  ATTATCHSTATUSPTR( stat );

  if(!dat) {
    ABORT ( stat, COHERENTESTIMATIONH_ENULL, COHERENTESTIMATIONH_MSGENULL );
  }

  LALFree(dat->detectors);

  for(i=0; i<dat->Ndetectors; i++) {
    if(dat->filters[i]) {
      TRY ( LALDestroyREAL8IIRFilter(stat->statusPtr, dat->filters + i), stat );
    }
  }

  LALFree(dat->filters);

  DETATCHSTATUSPTR( stat );
  RETURN( stat );

}
