/* Write function definitions for reading in the Crab ephemeris and */
/* computing the coefficients of the spin down model for f3 and f4  */ 

/* Matt Pitkin 10/03/04 */
/* LALTimingNoiseHeterodyne function has now been changed so that it 
   doesn't have to heterodyne at the SSB, cos the data has already been
	 transformed to the SSB - soooo much easier! */
	 
/* Matt Pitkin 26/03/04 */
/* The function's inputs params have be swithced about to conform to 
	 the LAL Spec */
/*
$Id$
*/

#include <math.h>
#include <stdio.h>

#include <lal/LALConstants.h>
#include "HeterodyneCrabPulsar.h"

#define LAL_DAYS_PER_SECOND 1.1574074074074074074074074074074074e-5L /* 1 second as a fraction of a day */

/* DEFINE RCS ID STRING */
NRCSID( HETERODYNECRABPULSARC, "$Id$"); 

void 
LALGetCrabEphemeris	( LALStatus			*status,
			  CrabSpindownParamsInput	*output,
				GetCrabEphemerisInput *input )
{
  UINT4 i=0, j=0;
  
  REAL8 MJD;
  REAL8 MJDVec[1000];
  
  /* array of the time of the pulse from the ephemeris  */
  /* minus the time of the first ephemeris data point */
  REAL8 GPStArr;
  REAL8 tArrVec[1000];
  
  REAL8 f0;
  REAL8 f1;

  FILE *fp;
  
  INITSTATUS( status, "LALGetCrabEphemeris", HETERODYNECRABPULSARC);
  ATTATCHSTATUSPTR (status);
  
  /* Check Input Arguments */
  ASSERT(output != (CrabSpindownParamsInput *)NULL, status, HETERODYNECRABPULSARH_ENULLOUTPUT,
  HETERODYNECRABPULSARH_MSGENULLOUTPUT);
  
  ASSERT(input->filename != NULL, status, HETERODYNECRABPULSARH_EEPHEMERISFILENAME,
  HETERODYNECRABPULSARH_MSGEEPHEMERISFILENAME);
  
  if((fp = fopen(input->filename,"r")) == NULL){
    fprintf(stderr,"Unable to open file!");
    exit(1);
  }
  
  /* read in ephermeris data */
  while(!feof(fp)){
    fscanf(fp,"%lf%lf%lf%lf",&MJD,&GPStArr,&f0,&f1);
    output->f0->data[j] = f0;
    output->f1->data[j] = f1;
    MJDVec[j] = MJD;
    tArrVec[j] = GPStArr;
    j++;
  }
  
  fclose(fp);
  
	ASSERT(j != 0, status, HETERODYNECRABPULSARH_ENUMEPHEMERISDATA,
  HETERODYNECRABPULSARH_MSGENUMEPHEMERISDATA);
  
	output->numOfData = j;
	
	/* check that the number of times in the ephemeris file (N2) is equal to the */
  /* number of lines read in (i-1) if not give error and exit                */
  /*if ( N2 != i-1 ){
    fprintf(stderr,"Error in ephemeris file, check or update\n");
    exit(1);
  }*/
  
  /* convert time in MJD to GPS time */
  for(i=0;i<j;i++){
    /* add leap seconds */
    if (MJDVec[i] >= 44786) MJDVec[i] += LAL_DAYS_PER_SECOND;
    
    if (MJDVec[i] >= 45151) MJDVec[i] += LAL_DAYS_PER_SECOND;

    if (MJDVec[i] >= 45516) MJDVec[i] += LAL_DAYS_PER_SECOND;

    if (MJDVec[i] >= 46247) MJDVec[i] += LAL_DAYS_PER_SECOND;

    if (MJDVec[i] >= 47161) MJDVec[i] += LAL_DAYS_PER_SECOND;

    if (MJDVec[i] >= 47892) MJDVec[i] += LAL_DAYS_PER_SECOND;

    if (MJDVec[i] >= 48257) MJDVec[i] += LAL_DAYS_PER_SECOND;

    if (MJDVec[i] >= 48804) MJDVec[i] += LAL_DAYS_PER_SECOND;

    if (MJDVec[i] >= 49169) MJDVec[i] += LAL_DAYS_PER_SECOND;

    if (MJDVec[i] >= 49534) MJDVec[i] += LAL_DAYS_PER_SECOND;
  
    if (MJDVec[i] >= 50083) MJDVec[i] += LAL_DAYS_PER_SECOND;
                                   
    if (MJDVec[i] >= 50630) MJDVec[i] += LAL_DAYS_PER_SECOND;
                                     
    if (MJDVec[i] >= 51179) MJDVec[i] += LAL_DAYS_PER_SECOND;
    /* need to add extra as leap seconds are added */
    /*printf("%lf\n",fmod(MJDVec[i],50000));*/ 
    output->tArr->data[i] = (MJDVec[i] - 44244.0)*24.0*60.0*60.0 + tArrVec[i];   
  }
  
  /* LALCheckMemoryLeaks();*/
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALComputeFreqDerivatives	( LALStatus			*status,
				  CrabSpindownParamsOutput	*output,
					CrabSpindownParamsInput	*input )
{
  REAL8 t, t1; /* temp time values to calculate freq derivs */
  REAL8 a, b, c; /* values used to calculate frequency derivs */
  REAL8 tempf2, tempf3; /* temp values of f2 and f3 */
  UINT4 i;
  REAL8 f0;
  REAL8 f1;
  REAL8 f2;  /* second deriv of nu */
  REAL8 f3;  /* third deriv of nu */
  REAL8 f4;  /* fourth deriv of nu */
  REAL8 tArr;	/* arrival time of pulse (GPS) */
  REAL8 phase0, phase; /* phase of signal */
  REAL8 phaseNew;
  
  INITSTATUS( status, "LALComputeFreqDerivatives", HETERODYNECRABPULSARC);
  ATTATCHSTATUSPTR (status);
  
  /* Check Input Arguments */
  ASSERT(input != (CrabSpindownParamsInput *)NULL, status, HETERODYNECRABPULSARH_ENULLINPUT,
  HETERODYNECRABPULSARH_MSGENULLINPUT);
  
  ASSERT(output != (CrabSpindownParamsOutput *)NULL, status, HETERODYNECRABPULSARH_ENULLOUTPUT,
  HETERODYNECRABPULSARH_MSGENULLOUTPUT);
  
	ASSERT(input->numOfData != 0, status, HETERODYNECRABPULSARH_ENUMEPHEMERISDATA,
  HETERODYNECRABPULSARH_MSGENUMEPHEMERISDATA);
  
	input->f1->data[0] = input->f1->data[0]*1.e-15;  /* convert f1 to correct units */ 
  phase0 = 0.0;
  phaseNew = 0.0;
  
  for (i=1;i<input->numOfData;i++){
    tArr = input->tArr->data[i-1];
    t = input->tArr->data[i]-input->tArr->data[i-1]; 
    t1 = t/2.0;
    input->f1->data[i] = input->f1->data[i]*1.0e-15;

    /* calculate values for spline */ 
    b = input->f0->data[i-1] + input->f1->data[i-1]*t1 - input->f0->data[i] +
    input->f1->data[i]*t1;
    c = input->f1->data[i-1] - input->f1->data[i];
    
    /* calculate second and third derivatives of f with constrained f0 and f1 */
    
    tempf2 = -(3.0/(2.0*t1*t1))*b - (1.0/(2.0*t1))*c;
    tempf3 = (3.0/(2.0*t1*t1*t1))*b;

    /* calculate the phase */
    phase = phase0 +input->f0->data[i-1]*t + (input->f1->data[i-1]/2.0)*t*t +
    (tempf2/6.0)*t*t*t + (tempf3/24.0)*t*t*t*t;

    /* round phase to nearest integer */
		/* (phase is tracked well enough using the values of f0, f1, f2 and f3 
		calculated above to be able to constrain the phase after time t to the nearest
		integer) */
    if((ceil(phase)-phase<0.5) && (ceil(phase)-phase>=0.0))
      phaseNew = ceil(phase);
    else
      phaseNew = floor(phase);
    /* fprintf(stdout,"%lf\t%lf\n",phase,phaseNew); */
    
    /* recalculate spindown params with constrained phase */
    a = phase0 + input->f0->data[i-1]*t1 + (input->f1->data[i-1]/2.0)*t1*t1 -
    phaseNew + input->f0->data[i]*t1 - (input->f1->data[i]/2.0)*t1*t1;
    
    /* calculate the second, third and fourth derivatives of frequency */ 
    f0 = input->f0->data[i-1];
    f1 = input->f1->data[i-1];
    f2 = -(15.0/(2.0*t1*t1*t1))*a - (3.0/(2.0*t1*t1))*b + (3.0/(4.0*t1))*c;
    f3 = (45.0/(2.0*t1*t1*t1*t1))*a + (3.0/(2.0*t1*t1*t1))*b - (15.0/(4.0*t1*t1))*c;
    f4 = -(45.0/(2.0*t1*t1*t1*t1*t1))*a + (15.0/(4.0*t1*t1*t1))*c;
    
    output->f0->data[i-1] = f0;
    output->f1->data[i-1] = f1;
    output->f2->data[i-1] = f2;
    output->f3->data[i-1] = f3;
    output->f4->data[i-1] = f4;
    output->tArr->data[i-1] = tArr;
		output->numOfData = input->numOfData;
  }
  
  /* LALCheckMemoryLeaks(); */
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

void
LALSetSpindownParams	( LALStatus			*status,
			  ParamsForHeterodyne		*output,
				CrabSpindownParamsOutput	*input,
				LIGOTimeGPS			epoch)
{
  UINT4 i = 0;
  REAL8 dataEpoch;
  /*REAL8 dt;*/
  
  INITSTATUS( status, "LALSetSpindownParams", HETERODYNECRABPULSARC);
  ATTATCHSTATUSPTR (status);
  
  /* Check Input Arguments */ 
  ASSERT(input != (CrabSpindownParamsOutput *)NULL, status, HETERODYNECRABPULSARH_ENULLINPUT,
  HETERODYNECRABPULSARH_MSGENULLINPUT);
  
	ASSERT(input->numOfData != 0, status, HETERODYNECRABPULSARH_ENUMEPHEMERISDATA,
  HETERODYNECRABPULSARH_MSGENUMEPHEMERISDATA);
	
  dataEpoch = (REAL8)epoch.gpsSeconds + ((REAL8)epoch.gpsNanoSeconds/1.0e9);
    
  for(i=0;i<input->numOfData;i++){
    if((dataEpoch<input->tArr->data[i+1])&&(dataEpoch>=input->tArr->data[i])){
      output->f0 = input->f0->data[i];
			output->f1 = input->f1->data[i];
			output->f2 = input->f2->data[i];
			output->f3 = input->f3->data[i];
			output->f4 = input->f4->data[i];
	
			output->epoch = input->tArr->data[i];
			
			break;
    }
  }
  
  DETATCHSTATUSPTR(status);
  RETURN(status);
}

/* function that performs the timing noise heterodyne but no longer at SSB */
/* also it now does one point at a time */
void
LALTimingNoiseHeterodyne	( LALStatus		*status,
				  TNHeterodyneOutput	*output,
					TNHeterodyneInput	*input,
				  ParamsForHeterodyne	*params )
{
	UINT4 i;
  REAL8 DPhi;	/* phase difference to be removed */
  REAL8 t1;		 
  REAL8 t;
  REAL8 phi, phi0; /* phases with and without timing noise */
  INT8 length;
  
  INITSTATUS( status, "LALTimingNoiseHeterodyne", HETERODYNECRABPULSARC);
  ATTATCHSTATUSPTR (status);
  
  /* Check Input Arguments */ 
  ASSERT(input != (TNHeterodyneInput *)NULL, status, HETERODYNECRABPULSARH_ENULLINPUT,
  HETERODYNECRABPULSARH_MSGENULLINPUT);
  
  ASSERT(output != (TNHeterodyneOutput *)NULL, status, HETERODYNECRABPULSARH_ENULLOUTPUT,
  HETERODYNECRABPULSARH_MSGENULLOUTPUT);
  
  ASSERT(params->f0 > 0, status, HETERODYNECRABPULSARH_EINVALIDF0,
  HETERODYNECRABPULSARH_MSGEINVALIDF0);

  t = (REAL8)input->epoch.gpsSeconds+(REAL8)input->epoch.gpsNanoSeconds/1.e9;

  /* Dt between time of data point and the epoch of the FIRST data point */
  /* epoch of the FIRST data point is always the GPS time of that point  */
  t -= (REAL8)input->t0;    
    
  /* Dt between time of data point and the epoch of that data point */
  t1 = ((REAL8)input->epoch.gpsSeconds+(REAL8)input->epoch.gpsNanoSeconds/1.e9)-(REAL8)params->epoch;
		
  /* calculate phase difference between signal with timing noise and one without */
  /* phi - phi0 */
  /* phi0 = (input->f0*t + 0.5*input->f1*t*t + (1.0/6.0)*input->f2*t*t*t);/*-(input->f0*t0 + 0.5*input->f1*t0*t0);*/
  /* input f0 is already at GW freq (i.e. spin freq * 2) */
	phi0 = input->f0*t + 0.5*input->f1*t*t + (1.0/6.0)*input->f2*t*t*t;
    
  phi = 2.0*(params->f0*t1 + 0.5*params->f1*t1*t1 + (params->f2/6.0)*t1*t1*t1 +
    (params->f3/24.0)*t1*t1*t1*t1 +
    (params->f4/120.0)*t1*t1*t1*t1*t1);
    
  DPhi = phi-phi0;
  	 
  DPhi = 2.0*(REAL8)LAL_PI*fmod(DPhi,1.0);
	
	output->Dphase = DPhi;
	    
  /* Heterodyne to remove timing noise */
  output->Vh.re = (REAL8)input->Vh.re*cos(-DPhi)
    			    -(REAL8)input->Vh.im*sin(-DPhi);
  output->Vh.im = (REAL8)input->Vh.re*sin(-DPhi)
              +(REAL8)input->Vh.im*cos(-DPhi);
    
  DETATCHSTATUSPTR(status);
  RETURN(status);
}
