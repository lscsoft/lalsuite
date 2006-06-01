/* This driver code performs the timing noise heterodyne for 
   the Crab pulsar on completely heterodyned data.
	 No calibration or noise estimation is needed or performed */
	 
/* Matt Pitkin (26/03/04) CrabTNHeterodyne.c v0.1 */

/*
$Id$
*/

/* headers */
#include <stdio.h>
#include <string.h>
#include <math.h>

/* lal headers */
#include <lal/LALStdlib.h>
#include <lal/AVFactories.h>
#include <lal/LALConstants.h>
#include <lal/BinaryPulsarTiming.h>

/* Crab heterodyne specific headers */
#include "HeterodyneCrabPulsar.h"

INT4 lalDebugLevel = 1;

#define EPHEMFILE "crab_ephemeris.txt"
#define MAXLENGTH 200000
#define NUM 1000 /* max num of lines in crab ephem file (ok for next 61 years
										file contains 266 lines as of 15 Jan 2004) */
#define EFILE "/home/matthew/lscsoft/lal/packages/pulsar/test/earth05-09.dat"
#define SFILE "/home/matthew/lscsoft/lal/packages/pulsar/test/sun05-09.dat"

int main(int argc, char *argv[]){
	static LALStatus status;
	FILE *fpin=NULL, *fpout=NULL;
	CHAR inputFile[256];
	CHAR outputFile[256];
	
	FILE *phifp=NULL;
	
	CHAR psrInput[256]; /* pulsar input file containing params f0, f1 etc. */
	FILE *psrfp=NULL;
	REAL8 f0;
	REAL8 f1;
	REAL8 f2;
	REAL8 fepoch;
	REAL8 val;
	CHAR txt[10];
	
	UINT4 i=0, j; /* counter */
	REAL8Vector *time=NULL;
	COMPLEX16Vector *B=NULL;
	COMPLEX16Vector *output=NULL;
	
	/* vars for Crab */
	GetCrabEphemerisInput input;
	CrabSpindownParamsInput crabEphemerisData;
  CrabSpindownParamsOutput crabOutput;
  ParamsForHeterodyne hetParams;
	TNHeterodyneInput TNInput;
  TNHeterodyneOutput TNOutput;
	
	LIGOTimeGPS dataEpoch;
	
	EphemerisData *edat = NULL;
	CHAR det[5];
	
	/* check command line inputs */
	if(argc!=4){
		fprintf(stderr, "Wrong number of input params:\n\tinputfile outputfile detector\n");
		return 0;
	}
	
	sprintf(inputFile, "%s", argv[1]);
	sprintf(outputFile, "%s", argv[2]); 
	sprintf(det, "%s", argv[3]);
	
	/* set detector */
	if (!strcmp(det,"H1") || !strcmp(det,"H2")) 
    hetParams.detector = lalCachedDetectors[LALDetectorIndexLHODIFF];      
  else if (   !strcmp(det,"L1"))
   	hetParams.detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if (   !strcmp(det,"GEO"))
   	hetParams.detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];             
  else
  {
    fprintf(stderr,"Error DETECTOR must be either H1, H1, L1, or GEO\n");
    return 2;
  }
	 
	/* read in Crab pulsar parameters used in heterodyning */
	sprintf(psrInput, "J0534+2200");
	psrfp = fopen(psrInput, "r");
	
	while (2==fscanf(psrfp,"%s %lf", &txt, &val))
  {
    if( !strcmp(txt,"f0") || !strcmp(txt,"F0")) {
      f0 = val;
    }
    else if( !strcmp(txt,"f1") || !strcmp(txt,"F1")) {
      f1 = val;
    }
    else if( !strcmp(txt,"f2") || !strcmp(txt,"F2")) {
      f2 = val;
    }
    else if( !strcmp(txt,"fepoch") || !strcmp(txt,"FEPOCH")) {
      fepoch = LALTDBMJDtoGPS(val);
    }
		else if (!strcmp(txt,"RA") || !strcmp(txt,"ra")){
			hetParams.alpha = val;
		}
		else if (!strcmp(txt,"DEC") || !strcmp(txt,"dec")){
			hetParams.delta = val;
		}
  }
	
	fclose(psrfp);
	
	/* allocate memory for crab ephemeris */
  crabEphemerisData.f1 = NULL;
  LALDCreateVector( &status, &crabEphemerisData.f1, NUM);
  
  crabEphemerisData.f0 = NULL;
  LALDCreateVector( &status, &crabEphemerisData.f0, NUM);
  
  crabEphemerisData.tArr = NULL;
  LALDCreateVector( &status, &crabEphemerisData.tArr, NUM);
  
  crabOutput.tArr = NULL;
  LALDCreateVector( &status, &crabOutput.tArr, NUM);
      
  crabOutput.f0 = NULL;
  LALDCreateVector( &status, &crabOutput.f0, NUM);
  
  crabOutput.f1 = NULL;
  LALDCreateVector( &status, &crabOutput.f1, NUM);
  
  crabOutput.f2 = NULL;
  LALDCreateVector( &status, &crabOutput.f2, NUM);
  
  crabOutput.f3 = NULL;
  LALDCreateVector( &status, &crabOutput.f3, NUM);
  
  crabOutput.f4 = NULL;
  LALDCreateVector( &status, &crabOutput.f4, NUM);

	input.filename = EPHEMFILE;

	/* read Crab ephemeris data */
	LALGetCrabEphemeris( &status, &crabEphemerisData, &input );
	
	/* compute frequency derivatives */
	LALComputeFreqDerivatives ( &status, &crabOutput, &crabEphemerisData );
	
	/* allocate mem for input data */
	LALDCreateVector(&status, &time, MAXLENGTH);
	LALZCreateVector(&status, &B, MAXLENGTH);
	
	/* read in input data */
	fpin = fopen(inputFile, "r");
	
	i=0;
	while(!feof(fpin)){
		if(i >= MAXLENGTH){
			fprintf(stderr, "MAXLENGTH at %d is not big enough.\n", MAXLENGTH);
			return 0;
		}
		
		fscanf(fpin, "%lf%lf%lf", &time->data[i], &B->data[i].re, &B->data[i].im);
		i++;
	}
	
	fclose(fpin);
	
	/* set input params for timing noise heterodyne (2*pulsar spin freq) */
	TNInput.f0 = 2.0*f0;
	TNInput.f1 = 2.0*f1;
	TNInput.f2 = 2.0*f2;
	TNInput.t0 = (fepoch-44244.)*86400.; /* convert from time in MJD to secs in TDB */
	
	fpout = fopen(outputFile, "w");
	phifp = fopen("DPhase.txt", "w");
	
	/* set up LALBarycenter */
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));    
  (*edat).leap = 13;
  (*edat).ephiles.earthEphemeris = EFILE;
  (*edat).ephiles.sunEphemeris = SFILE;	 
  LALInitBarycenter(&status, edat);
  hetParams.edat = edat;
	
	/* perform timing noise heterodyne */
	for(j=0;j<i-1;j++){
		dataEpoch.gpsSeconds = (INT8)floor(time->data[j]);
		dataEpoch.gpsNanoSeconds = 0;
		
		TNInput.epoch.gpsSeconds = dataEpoch.gpsSeconds;
		TNInput.epoch.gpsNanoSeconds = 0;
		
		TNInput.Vh.re = B->data[j].re;
		TNInput.Vh.im = B->data[j].im;
		
		/* set values for timing noise removal */
  	LALSetSpindownParams( &status, &hetParams, &crabOutput, dataEpoch );
		
		if(j==0){
			REAL8 dt;
			dt = time->data[j]-hetParams.epoch;
			fprintf(stderr, "f0 = %f,\nf1 = %e,\nf2 = %e,\nf3 = %e,\nf4 = %e.\n", 
			hetParams.f0, hetParams.f1, hetParams.f2, hetParams.f3, hetParams.f4);
			fprintf(stderr, "dt = %f, \n(1/120)*f4*dt^5 = %e.\n", dt, 
			(1.0/120.0)*hetParams.f4*dt*dt*dt*dt*dt);
		}
		
		/* perform TN heterodyne */
		LALTimingNoiseHeterodyne( &status, &TNOutput, &TNInput, &hetParams);
		
		fprintf(fpout, "%lf\t%e\t%e\n", time->data[j], TNOutput.Vh.re,
		TNOutput.Vh.im);
		
		fprintf(phifp, "%f\t%f\t%f\n", TNOutput.phi0, TNOutput.phi1, TNOutput.Dphase);
	}
	
	fclose(fpout);
	fclose(phifp);
	
	LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
	
	/* destroy vectors */
	LALZDestroyVector(&status, &B);
	LALDDestroyVector(&status, &time);
	LALDDestroyVector(&status, &crabEphemerisData.f1);
	LALDDestroyVector(&status, &crabEphemerisData.f0);
	LALDDestroyVector(&status, &crabEphemerisData.tArr);
	LALDDestroyVector(&status, &crabOutput.tArr);
	LALDDestroyVector(&status, &crabOutput.f0);
	LALDDestroyVector(&status, &crabOutput.f1);
	LALDDestroyVector(&status, &crabOutput.f2);
	LALDDestroyVector(&status, &crabOutput.f3);
	LALDDestroyVector(&status, &crabOutput.f4);
	
	return 0;
}
