/*----------------------------------------------------------------------- 
 * 
 * File Name: FindChirpBCVCFilter.c
 *
 * Author: 
 * 
 * Revision: $Id$
 * 
 *-----------------------------------------------------------------------
 */

#if 0
<lalVerbatim file="FindChirpBCVCFilterCV">
Author: 
$Id$
</lalVerbatim>

<lalLaTeX>
\input{FindChirpBCVCFilterCDoc}

\vfill{\footnotesize\input{FindChirpBCVCFilterCV}}
</lalLaTeX>
#endif

#include <math.h>
#include <lal/LALStdio.h>
#include <lal/LALStdlib.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/AVFactories.h>
#include <lal/FindChirp.h>
#include <lal/FindChirpBCV.h>

double rint(double x);

NRCSID (FINDCHIRPBCVCFILTERC, "$Id$");


/* <lalVerbatim file="FindChirpBCVCFilterCP"> */
void
LALFindChirpBCVCFilterSegment (
    LALStatus                  *status,
    SnglInspiralTable         **eventList,
    FindChirpFilterInput       *input,
    FindChirpFilterParams      *params
    )
/* </lalVerbatim> */
{
  UINT4                 j, k, kFinal/*, kOpt*/;
  UINT4                 numPoints;
  UINT4                 deltaEventIndex = 0;
  UINT4                 ignoreIndex;
  REAL4                 myfmin;
  REAL4                 deltaT, deltaF;
  REAL4                 norm;
  REAL4                 modqsqThresh;
  REAL4                 rhosqThresh;
  REAL4                 mismatch;

  UINT4                 eventStartIdx = 0;
  REAL4                 chirpTime     = 0;
  REAL4                 increment, nextBin, partSum;

  COMPLEX8             *qtilde        = NULL; 
  COMPLEX8             *qtildeBCVC    = NULL; 
  COMPLEX8             *q             = NULL; 
  COMPLEX8             *qBCVC         = NULL;
  COMPLEX8             *inputData     = NULL;
  COMPLEX8             *inputDataBCVC = NULL;
  COMPLEX8             *tmpltSignal   = NULL;
  SnglInspiralTable    *thisEvent     = NULL;
  LALMSTUnitsAndAcc     gmstUnits;
  REAL4                 a1 = 0.0;
  REAL4                 b1 = 0.0;                  
  REAL4                 b2 = 0.0;                  
  REAL4                 templateNorm;

  REAL4                 Num1, Num2, Den1, Den2;
  REAL4                 omega, InvTan1, InvTan2;
  REAL4                 m, psi0, psi3, fFinal;

  
  /* for BCVC */
  REAL4 alphaMax, thetab, thetav, temp, w, thetac, alpha2Store, alphaC2Store;
  REAL4 rhoConstraint = 0.0;
  REAL4 rhoUnconstraint = 0.0;
  REAL4 deltaTPower2by3 ;
  REAL4 deltaTPower1by6 ;
  REAL4 ThreeByFive = 3./5.;
  REAL4 FiveByThree = 5./3.;
  INITSTATUS( status, "LALFindChirpBCVCFilter", FINDCHIRPBCVCFILTERC );
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
  ASSERT( params->rhosqThresh >= 0, status,
      FINDCHIRPH_ERHOT, FINDCHIRPH_MSGERHOT );
  ASSERT( params->chisqThresh >= 0, status,
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

  /* if a rhosqVec vector has been created, check we can store data in it */
  if ( params->rhosqVec )
  {
    ASSERT( params->rhosqVec->data->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
    ASSERT( params->rhosqVec->data, status,
        FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  }

  /* make sure that the input structure exists */
  ASSERT( input, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the input structure contains some input */
  ASSERT( input->fcTmplt, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );
  ASSERT( input->segment, status, FINDCHIRPH_ENULL, FINDCHIRPH_MSGENULL );

  /* make sure that the template and the segment are both BCVC */
  ASSERT( input->fcTmplt->tmplt.approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );
  ASSERT( input->segment->approximant == BCV, status,
      FINDCHIRPH_EAPRX, FINDCHIRPH_MSGEAPRX );



  /* workspace vectors */
  q    = params->qVec->data;
  qBCVC = params->qVecBCV->data;
  qtilde    = params->qtildeVec->data;
  qtildeBCVC = params->qtildeVecBCV->data;


  /* template and data */
  inputData     = input->segment->data->data->data;
  inputDataBCVC  = input->segment->dataBCV->data->data;
  tmpltSignal   = input->fcTmplt->data->data;
  templateNorm  = input->fcTmplt->tmpltNorm;
  deltaT        = params->deltaT;
  deltaTPower2by3 = pow(params->deltaT, 2./3.);  
  deltaTPower1by6 = pow(params->deltaT, 1./6.);  



  /* number of points in a segment */
  if (params->qVec->length != params->qVecBCV->length)
  {
    ABORT(status, FINDCHIRPBCVH_EQLEN, FINDCHIRPBCVH_MSGEQLEN);
  }
  numPoints = params->qVec->length;

  /* set the gmst units and strictness */
  gmstUnits.units = MST_HRS;
  gmstUnits.accuracy = LALLEAPSEC_STRICT;

  /* 
   * template parameters, since FindChirpBCVCFilterSegment is run
   * for every template
   */
  psi0 = input->fcTmplt->tmplt.psi0;  
  psi3 = input->fcTmplt->tmplt.psi3;
  fFinal = input->fcTmplt->tmplt.fFinal;

  {
  /* Calculate deltaEventIndex : the only acceptable clustering */
  /* is "window" method, for BCVC                                */
    if ( params->clusterMethod == window )
    {
      deltaEventIndex=(UINT4) rint((params->clusterWindow/params->deltaT)+1.0);
    }
    else if ( params->clusterMethod == tmplt ) 
    {
      ABORT( status, FINDCHIRPBCVH_ECLUW, FINDCHIRPBCVH_MSGECLUW );
    }
       

    /* ignore corrupted data at start and end */
    ignoreIndex = ( input->segment->invSpecTrunc / 2 ) + deltaEventIndex;

    if ( lalDebugLevel & LALINFO )
    {
      CHAR newinfomsg[256];


      LALSnprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
          "chirp time = %e seconds => %d points\n"
          "invSpecTrunc = %d => ignoreIndex = %d\n",
          chirpTime, deltaEventIndex,
          input->segment->invSpecTrunc, ignoreIndex );
      LALInfo( status, newinfomsg );
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

  if ( lalDebugLevel & LALINFO )
  {
    CHAR newinfomsg[256];

    LALSnprintf( newinfomsg, sizeof(newinfomsg) / sizeof(*newinfomsg),
        "filtering from %d to %d\n",
        ignoreIndex, numPoints - ignoreIndex );
    LALInfo( status, newinfomsg );
  }

  /* k that corresponds to fFinal */
  deltaF = 1.0 / ( (REAL4) params->deltaT * (REAL4) numPoints );
  kFinal = fFinal / deltaF < numPoints/2 ? floor(fFinal / deltaF) 
    : floor(numPoints/2); 
  myfmin = input->segment->fLow;

  /* assign the values to a1, b1 and b2 */
  a1 = input->segment->a1->data[kFinal];
  b1 = input->segment->b1->data[kFinal];
  b2 = input->segment->b2->data[kFinal];


  /*
   *
   * compute qtilde, qtildeBCVC, and q, qBCVC 
   * using the correct combination of inputData and inputDataBCVC
   *
   */


  memset( qtilde,    0, numPoints * sizeof(COMPLEX8) );
  memset( qtildeBCVC, 0, numPoints * sizeof(COMPLEX8) );

  /* qtilde positive frequency, not DC or nyquist */
  for ( k = 1; k < numPoints/2; ++k )
  {
    REAL4 r    = a1 * inputData[k].re;
    REAL4 s    = a1 * inputData[k].im;    
    REAL4 rBCVC = b1 * inputData[k].re + b2 * inputDataBCVC[k].re;
    REAL4 sBCVC = b1 * inputData[k].im + b2 * inputDataBCVC[k].im; 
    REAL4 x = tmpltSignal[k].re;
    REAL4 y = 0.0 - tmpltSignal[k].im; /* note complex conjugate */     

    qtilde[k].re = r * x - s * y ;
    qtilde[k].im = r * y + s * x ;
    qtildeBCVC[k].re = rBCVC * x - sBCVC * y ;
    qtildeBCVC[k].im = rBCVC * y + sBCVC * x ;
  }

  /* inverse fft to get q, and qBCVC */
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

  rhosqThresh = params->rhosqThresh;

  params->norm = norm = deltaT / ((REAL4) numPoints) ;
  /* notice difference from corresponding factor in the sp templates: */
  /* no factor of 4 (taken care of in inputData and inputDataBCVC      */
  /* and no segnorm, since we already multiplied by a1, b1 and b2.    */


  /* normalized snr threhold */
  modqsqThresh = rhosqThresh / norm ;

  mismatch = 1.0 - input->fcTmplt->tmplt.minMatch;


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
      REAL4 modqsqBCVC = qBCVC[j].re * qBCVC[j].re + qBCVC[j].im * qBCVC[j].im ;
      REAL4 ImProd = 2.0 * ( - q[j].re * qBCVC[j].im + qBCVC[j].re * q[j].im ) ;

      REAL4 newmodqsq = ( 0.5 * sqrt( modqsqSP + modqsqBCVC + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCVC - ImProd ) ) *
        ( 0.5 * sqrt( modqsqSP + modqsqBCVC + ImProd ) +
          0.5 * sqrt( modqsqSP + modqsqBCVC - ImProd ) ) ;

      params->rhosqVec->data->data[j] = norm * newmodqsq;   
    }
  }



#if 0 
 { 

   FILE *K1File, *K2File, *K3File, *K4File, *V0File, *V1File,*V2File,
     *alphaFile,  *rho1File, *rho2File, *rho3File,
     *phaseFile, *phaseFileE, *alphaFileE, *omegaFile, *rhoFile, *omegaFileE;

   thetac = atan(-a1/b1);
   
   fprintf(stderr, "a1=%e b1=%e b2 =%e fFinal=%e\n", a1, b1, b2, fFinal);
   
   alphaMax = pow(fFinal, -2./3.)*pow(params->deltaT,-2./3.);
   
   
   thetab   = fabs(-(a1 * alphaMax)/(b2+b1*alphaMax));
   thetab   = atan(thetab);
   fprintf(stderr, "alphaMax -->thetab = %e\n", thetab);
   fprintf(stderr, "alphaMax -->thetac = %e\n", thetac);
   
   
   K1File=fopen("K1File.dat","w");/*ok*/
   K2File=fopen("K2File.dat","w");/*ok*/
   K3File=fopen("K3File.dat","w");/*ok*/
   K4File=fopen("K4File.dat","w");/*ok*/
   V0File=fopen("V0File.dat","w");/*ok*/
   V1File=fopen("V1File.dat","w");/*ok*/
   V2File=fopen("V2File.dat","w");/*ok*/
   rho1File=fopen("rho1File.dat","w");
   rho2File=fopen("rho2File.dat","w");
   rho3File=fopen("rho3File.dat","w");
   alphaFile=fopen("alphaFile.dat","w");
   omegaFile=fopen("omegaFile.dat","w");
   phaseFile=fopen("phaseFile.dat","w");

   phaseFileE=fopen("phaseFileE.dat","w");
   alphaFileE=fopen("alphaFileE.dat","w");
   omegaFileE=fopen("omegaFileE.dat","w");
   rhoFile=fopen("rhoFile.dat","w");
   
   
   for ( j = 0; j < numPoints; ++j )
     {      
       REAL4 K1 = q[j].re ;
       REAL4 K3 = q[j].im ;
       REAL4 K2 = qBCVC[j].re ;
       REAL4 K4 = qBCVC[j].im ;
       
       REAL4 V0 = K1*K1 + K2*K2 + K3*K3 + K4*K4;
       REAL4 V1 = K1*K1 + K3*K3 - K2*K2 - K4*K4;
       REAL4 V2 = 2.*(K1*K2 + K3*K4); /* signe - ??*/
       REAL4 modqsq = 0.5 * (V0+sqrt(V1*V1+V2*V2));
       
       thetav = atan2(V2,V1);
       thetav *= .5;
       
       Num1 = K2 + K3 ;
       Num2 = K2 - K3 ;
       Den1 = K1 - K4;
       Den2 = K1 + K4;
       
       InvTan1 = (REAL4) atan2(Num1, Den1);
       InvTan2 = (REAL4) atan2(Num2, Den2);
       
       fprintf(K1File,"%e\n",K1);
       fprintf(K2File,"%e\n",K2);
       fprintf(K3File,"%e\n",K3);
       fprintf(K4File,"%e\n",K4);
       fprintf(V0File,"%e\n",V0);
       fprintf(V1File,"%e\n",V1);
       fprintf(V2File,"%e\n",V2);
       
       fprintf(omegaFile,"%e\n",thetav);
       
       fprintf(alphaFile,"%e\n", -(b2 * tan(thetav))
	       / (a1 + b1* tan(thetav))*pow(params->deltaT*fFinal, 2.0/3.0));
       
       omega = 0.5 * InvTan1 + 0.5 * InvTan2 ;
       
       fprintf(phaseFileE,"%e\n",0.5 * InvTan1 - 0.5 * InvTan2);

       fprintf(phaseFile,"%e\n",0.5 * InvTan1 - 0.5 * InvTan2);


       fprintf(omegaFileE,"%e\n",omega);
       
       
       fprintf(alphaFileE,"%e\n", (- b2 * tan(omega) 
				   / ( a1 + b1 * tan(omega)))
	       *pow(params->deltaT*fFinal, 2.0/3.0) );
       
       
       fprintf(rho1File,"%e\n",sqrt(modqsq*norm));
       fprintf(rho2File,"%e\n",sqrt((V0 + V1)/2.*norm));
       fprintf(rho3File,"%e\n",sqrt((V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.*norm));
       
       
      if ((-thetab<thetav && thetav<0)){	  
	fprintf(rhoFile,"%e\n",sqrt(norm*modqsq));	
      }
      else if ((0 <= thetav && thetav < thetac)){
	
	fprintf(rhoFile, "%e\n", sqrt((V0 + V1)/2.*norm));	  
	
      }
      else if( (-LAL_PI_2-1e-6 < thetav && thetav<=-thetab  )||(thetac <= thetav && thetav< LAL_PI_2+1e-6)){
	
	fprintf(rhoFile,"%e\n",sqrt(norm*(V0+V1*cos(2*thetab)+V2*sin(2*thetab))/2.));
	
      }
      else 
	{
	  fprintf(stderr,"must never enter here\n");
	  exit(0);
	}
     }
   fclose(rhoFile);
   fclose(alphaFileE);

   fclose(K1File);
   fclose(K2File);
   fclose(K3File);
   fclose(K4File);

   fclose(V0File);
   fclose(V1File);
   fclose(V2File);
   
   fclose(alphaFile);

   fclose(phaseFile);
   fclose(omegaFile);
   fclose(phaseFileE);
   fclose(omegaFileE);

   fclose(rho1File);
   fclose(rho2File);
   fclose(rho3File);
 }
#endif






  thetac   = atan(-a1/b1);
  alphaMax = pow(fFinal, -2./3.)*pow(params->deltaT,-2./3.);
  thetab   = fabs(-(a1 * alphaMax)/(b2+b1*alphaMax));
  thetab   = atan(thetab);
  thetav = 0;
  w = 0;
  /* look for an event in the filter output */
  for ( j = ignoreIndex; j < numPoints - ignoreIndex; ++j )
    {
#if 0 
      REAL4 modqsqSP  = q[j].re * q[j].re + q[j].im * q[j].im ;
      REAL4 modqsqBCVC = qBCVC[j].re * qBCVC[j].re + qBCVC[j].im * qBCVC[j].im ;
      REAL4 ImProd = 2.0 * ( - q[j].re * qBCVC[j].im + qBCVC[j].re * q[j].im ) ;
      
      REAL4 newmodqsq = ( 0.5 * sqrt( modqsqSP + modqsqBCVC + ImProd ) +
			  0.5 * sqrt( modqsqSP + modqsqBCVC - ImProd ) ) *
	( 0.5 * sqrt( modqsqSP + modqsqBCVC + ImProd ) +
	  0.5 * sqrt( modqsqSP + modqsqBCVC - ImProd ) ) ;
#endif
      /*
	REAL4 V0 = K1*K1 + K2*K2 + K3*K3 + K4*K4;
	REAL4 V1 = K1*K1 + K3*K3 - K2*K2 - K4*K4;
	REAL4 V2 = 2.*(K1*K2 + K3*K4); 
	REAL4 modqsq = 0.5 * (V0+sqrt(V1*V1+V2*V2));
      */
      REAL4 V0 = q[j].re * q[j].re  
	+ qBCVC[j].re * qBCVC[j].re  
	+ q[j].im * q[j].im 
	+ qBCVC[j].im * qBCVC[j].im ;
      
      REAL4 V1 = q[j].re * q[j].re  
	+ q[j].im * q[j].im  
	- qBCVC[j].im * qBCVC[j].im 
	-  qBCVC[j].re * qBCVC[j].re;
      
      REAL4 V2 =  2 * ( q[j].re * qBCVC[j].re + qBCVC [j].im * q[j].im); 
      REAL4 rhoUnconstraint = .5*( V0 + sqrt(V1*V1 + V2*V2));


      /* <-- That piece of code is the constraint BCV code */
      thetav = 0.5 * atan2(V2,V1);
      /* In the future, we can get rid of that since it concerns only the unconstraint code*/
      alpha2Store  =     - b2 * tan(thetav)  / ( a1 + b1 * tan(thetav) );
      
      if (-thetab < thetav && thetav < 0){	  
  	rhoConstraint = rhoUnconstraint;
  	w = thetav;
	/*alphaC2Store =  - b2 * tan(thetav)  / ( a1 + b1 * tan(thetav) );*/
      }
      else if (0 <= thetav && thetav < thetac){
	rhoConstraint = 0.5 * (V0 + V1);	  
  	w = 0.;	
/*	alphaC2Store  = 0.;*/
      }
      else if( (-LAL_PI_2-1e-6 < thetav && thetav<=-thetab )||(thetac <= thetav && thetav< LAL_PI_2+1e-6)){
  	rhoConstraint = 0.5 * (V0 + V1*cos(thetab * 2) + V2*sin( 2 * thetab));
  	w = -thetab; 
/*	alphaC2Store =  - b2 * tan(-thetab)  / ( a1 + b1 * tan(-thetab) );*/
      }
      else 
	{
	  /* TODO ABORT();*/
	  fprintf(stderr,"we must never enter here \n");
	  exit(0);
	}
      /* --> this is the end of the constraint BCV code */
      /* Now we maximise only over constraint BCV. So it is difficult to compare with
	 unconstraint snr now. */
      
      if ( rhoConstraint > modqsqThresh )                  
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
          thisEvent->snr   = rhoConstraint;
/*	  thisEvent->tau5  = alpha2Store * deltaTPower2by3;*/
	  thisEvent->tau4  = rhoUnconstraint;
	  thisEvent->alpha = alpha2Store * deltaTPower2by3; 
        }
        else if ( ! params->clusterMethod == noClustering &&
		  j <= thisEvent->end_time.gpsSeconds + deltaEventIndex &&
		  rhoConstraint > thisEvent->snr )
	  {
	    /* if this is the same event, update the maximum */
	    thisEvent->end_time.gpsSeconds = j;
	    thisEvent->snr = rhoConstraint;
/*	    thisEvent->tau5 = alpha2Store * deltaTPower2by3;*/
	    thisEvent->tau4 = rhoUnconstraint;
	    thisEvent->alpha = alpha2Store * deltaTPower2by3; 
	  }
      else if (j > thisEvent->end_time.gpsSeconds + deltaEventIndex ||
	       params->clusterMethod == noClustering )
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


	  /* I have decided to not keep the coa phase for the time being since
	   it uses two atan2 and is not used.*/
	  /*  
	  Num1 = qBCVC[timeIndex].re + q[timeIndex].im ;
          Num2 = qBCVC[timeIndex].re - q[timeIndex].im ;
          Den1 = q[timeIndex].re - qBCVC[timeIndex].im ;
          Den2 = q[timeIndex].re + qBCVC[timeIndex].im ;

          InvTan1 = (REAL4) atan2(Num1, Den1);
          InvTan2 = (REAL4) atan2(Num2, Den2);

          thisEvent->coa_phase = 0.5 * InvTan1 - 0.5 * InvTan2 ;
	  */

          /* copy the template into the event */
          thisEvent->psi0   = (REAL4) input->fcTmplt->tmplt.psi0; 
          thisEvent->psi3   = (REAL4) input->fcTmplt->tmplt.psi3;
	  
          /* chirp mass in units of M_sun */
          thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
            pow( 3.0 / 128.0 / input->fcTmplt->tmplt.psi0 , ThreeByFive );
          m =  fabs(thisEvent->psi3) / 
            (16.0 * LAL_MTSUN_SI * LAL_PI * LAL_PI * thisEvent->psi0) ;
          thisEvent->eta = 3.0 / (128.0*thisEvent->psi0 * 
              pow( (m*LAL_MTSUN_SI*LAL_PI), FiveByThree) );
          thisEvent->f_final  = (REAL4) input->fcTmplt->tmplt.fFinal ;

          /* set the type of the template used in the analysis */
          LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
              "BCVC" );
	  
	  thisEvent->chisq     = 0;
	  thisEvent->chisq_dof = 0;
          
          thisEvent->sigmasq = sqrt( norm / a1 );
          thisEvent->eff_distance =
            input->fcTmplt->tmpltNorm / norm / thisEvent->snr;
          thisEvent->eff_distance = sqrt( thisEvent->eff_distance ) / deltaTPower1by6;

          thisEvent->snr *= norm;      
          thisEvent->snr = sqrt( thisEvent->snr );

          thisEvent->tau4 *= norm;      
          thisEvent->tau4 *= sqrt(thisEvent->tau4);      
          /* compute the time since the snr crossing */
          thisEvent->event_duration =  (REAL8) timeIndex - (REAL8) eventStartIdx;
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
          thisEvent->snr = rhoConstraint;
/*	  thisEvent->tau5  = alpha2Store * deltaTPower2by3;*/
	  thisEvent->tau4 = rhoUnconstraint;
	  thisEvent->alpha = alpha2Store  * deltaTPower2by3; 
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
    /*    
    Num1 = qBCVC[timeIndex].re + q[timeIndex].im ;
    Num2 = qBCVC[timeIndex].re - q[timeIndex].im ;
    Den1 = q[timeIndex].re - qBCVC[timeIndex].im ;
    Den2 = q[timeIndex].re + qBCVC[timeIndex].im ;
    
    InvTan1 = (REAL4) atan2(Num1, Den1);
    InvTan2 = (REAL4) atan2(Num2, Den2 );
    
    thisEvent->coa_phase = 0.5 * InvTan1 - 0.5 * InvTan2 ;
    */

    

    /*thisEvent->tau5  = alpha2Store * deltaTPower2by3;*/
    thisEvent->alpha = alpha2Store  * deltaTPower2by3; 
    
    /* copy the template into the event */
    thisEvent->psi0   = (REAL4) input->fcTmplt->tmplt.psi0;   
    thisEvent->psi3   = (REAL4) input->fcTmplt->tmplt.psi3;  
    /* chirp mass in units of M_sun */
    thisEvent->mchirp = (1.0 / LAL_MTSUN_SI) * LAL_1_PI *
      pow( 3.0 / 128.0 / input->fcTmplt->tmplt.psi0, ThreeByFive);
    thisEvent->f_final  = (REAL4) input->fcTmplt->tmplt.fFinal;
    m =  fabs(thisEvent->psi3) /
          (16.0 * LAL_MTSUN_SI * LAL_PI * LAL_PI * thisEvent->psi0) ;
    thisEvent->eta = 3.0 / (128.0*thisEvent->psi0 *
          pow( (m*LAL_MTSUN_SI*LAL_PI), FiveByThree) );



    /* set the type of the template used in the analysis */
    LALSnprintf( thisEvent->search, LIGOMETA_SEARCH_MAX * sizeof(CHAR),
        "BCVC" );


    thisEvent->chisq     = 0;
    thisEvent->chisq_dof = 0;

    thisEvent->sigmasq = sqrt( norm / a1 );
    thisEvent->eff_distance = input->fcTmplt->tmpltNorm / norm /
      thisEvent->snr;
    thisEvent->eff_distance = sqrt( thisEvent->eff_distance ) /  deltaTPower1by6;

    thisEvent->snr *=  norm ;   
    thisEvent->snr = sqrt( thisEvent->snr );
    
    thisEvent->tau4 *= norm;      
    thisEvent->tau4 *= sqrt(thisEvent->tau4);      
 
    /* compute the time since the snr crossing */
    thisEvent->event_duration = (REAL8) timeIndex - (REAL8) eventStartIdx;
    thisEvent->event_duration *= (REAL8) deltaT;

  }


  /* normal exit */
  DETATCHSTATUSPTR( status );
  RETURN( status );
}
