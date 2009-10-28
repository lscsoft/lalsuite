/*
*  Copyright (C) 2007 Matt Pitkin
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/*********************************************************************/
/* 	Modified 24/06/04 to read binary info from a table, which has been
		precreated containing all relevant pulsar params - like for isolated
		pulsars. Also to be modified so frame data is only read in once
		rather than re-read for each pulsar.

		Added stuff to include binary pulsars - Matt Pitkin 05/05/04      */
/*********************************************************************/

/*********************************************************************************/
/*                    Time domain heterodyning code for known pulsars (S2/S3)    */
/*                                                                               */
/*			               Réjean Dupuis                             */
/*                                                                               */
/*                  University of Glasgow - last modified 26/03/2004             */
/*********************************************************************************/
/*
$Id$
*/
#include "TargetedPulsars.h"
/* include binary pulsar timing header file */
#include "BinaryPulsarTiming.h"

#define CONDOR 1 
#define DBUG 1 /* set to 1 to run in debug mode with more print outs */

typedef struct
tagFilters{
  REAL8IIRFilter *iirFilter1Re;
  REAL8IIRFilter *iirFilter1Im;
  REAL8IIRFilter *iirFilter2Re;
  REAL8IIRFilter *iirFilter2Im;
  REAL8IIRFilter *iirFilter3Re;
  REAL8IIRFilter *iirFilter3Im;
} Filters;

/* define a function to make a short list of frame files needed for processing
- this will mean that FrFileINew does not have to open an extremely large list
of files just a few. The function takes in the file name of a file containing 
a long list of frames, it takes in the length of the data you want to look for
and it takes in the start times from which you want */
char *MakeFrameFileList(char *framelist, int length, int tstart);

int
main( int argc, char *argv[])
{
  /***** declare variables *****/
  static LALStatus status;
  FILE *psrfp=NULL,*fpout3=NULL;
  FILE *fpout[MAXPSRS]; /* pointer to array of file pointers for output files */
	FILE *fpolg, *fpsf, *fpab;
  FrFile *frfile;         
  FrVect *frvect;    
  INT4 starts[MAXSTART],nstarts=0, count=0, tbase=INT_MAX, k=0;
  INT4 opencount=0, flag = 0, Nab, jobnum, i, irun, iflag;
  INT4 pnum, dec_deg, dec_min, ra_hrs, ra_min, counter = 0, i_pulsar = 0, j;
  char psrinput[256], psrout[256][MAXPSRS], chname[256];
  char framelist[256],earthfile[256],sunfile[256], outfile[256];
  char olgfile[256],sffile[256], abfile[256];
  char jnkstr1[256],jnkstr2[256],jnkstr3[256], jnkstr4[256],jnkstr5[256];   
  char pname[MAXPSRS][256], ppmra[30], ppmdec[30], pposepoch[16], ppepoch[16];  
  COMPLEX16 B; 
  REAL4 wc;
  REAL8Vector *alpha=NULL, *beta=NULL, *tab=NULL;
  REAL8 tt;
  REAL8 RA[MAXPSRS], DEC[MAXPSRS], f0[MAXPSRS], f1[MAXPSRS], f2[MAXPSRS];
	INT4 jtemp=0;
	
	/************* binary specific variables ********************/
	REAL8 T0[MAXPSRS], Tasc[MAXPSRS], Pb[MAXPSRS], Pbdot[MAXPSRS], x[MAXPSRS];
	REAL8 xdot[MAXPSRS], w0[MAXPSRS], wdot[MAXPSRS], e[MAXPSRS], edot[MAXPSRS];
	REAL8 eps1[MAXPSRS], eps1dot[MAXPSRS], eps2[MAXPSRS], eps2dot[MAXPSRS];
	REAL8 xpbdot[MAXPSRS], gamma[MAXPSRS], s[MAXPSRS], M[MAXPSRS], m2[MAXPSRS];
	REAL8 dr[MAXPSRS], dth[MAXPSRS], a0[MAXPSRS], b0[MAXPSRS];
	
	char model[MAXPSRS][5];
	
	REAL8 pT0, pTasc, pPb, pPbdot, px, pxdot, pw0, pwdot, pe, pedot, peps1;
	REAL8 peps1dot, peps2, peps2dot, pxpbdot, pgamma, ps, pM, pm2, pdr, pdth;
	REAL8 pa0, pb0;
	
	char pmodel[5]="";
  /***********************************************************/
	
	REAL8 a, b, ab, tg, f=0, jnknum1;
  REAL8Vector *H0pul=NULL, *phaseH0pul=NULL;
	REAL8Vector *C0pul=NULL, *phaseC0pul=NULL;
	
	/* calibration variables only used if H2 during S2 run */
	REAL8Vector *H0pul2=NULL, *phaseH0pul2=NULL;
	REAL8Vector *C0pul2=NULL, *phaseC0pul2=NULL;
	/*******************************************************/
	
	REAL8 ra_sec, dec_sec, pf0, pmra, pmdec,ffepoch, posepoch, fepoch_gps;
  REAL8 pf1, pf2, s_ra, s_dec;
  LIGOTimeGPS fepoch[MAXPSRS],epoch; 
  UINT4 npts;
  HeterodyneInput input;    HeterodyneOutput output;    HeterodyneParams params; 
  EphemerisData *edat = NULL;
  
	COMPLEX16ZPGFilter *zpgFilter1 = NULL;
	COMPLEX16ZPGFilter *zpgFilter2 = NULL;
	COMPLEX16ZPGFilter *zpgFilter3 = NULL;
  
	/* make filters an array, so there is one for each pulsar */
	/*REAL8IIRFilter *iirFilter1Re = NULL;   REAL8IIRFilter *iirFilter1Im = NULL;
	REAL8IIRFilter *iirFilter2Re = NULL;   REAL8IIRFilter *iirFilter2Im = NULL;
  REAL8IIRFilter *iirFilter3Re = NULL;   REAL8IIRFilter *iirFilter3Im = NULL; */
	Filters iirFilters[MAXPSRS];

	CalibrateFineParams calParams;    CalibrateFineInput calInput;  CalibrateFineOutput calOutput; 
  REAL8 REF_MJD, REF_GPS_SECS;
	
  /* get command line arguments and print error if there is wrong number of arguments.  if 
  the variable condor is defined, then read in the first argument jobnum from stdin */

  #ifndef CONDOR
  if (argc !=7 || (jobnum=atoi(argv[1]))<0 || jobnum>99999){
    int a;
    fprintf(stderr,"Syntax:\n\t%s N DIR1 DIR2 DETECTOR RUN InjectionsFlag where 0<=N<=99999.\n",argv[0]);
    fprintf(stderr,"Files used are jobdata.N, jobtimes.N, and pulsar.list\n");
    fprintf(stderr,"Input files are taken from DIR1\n");
    fprintf(stderr,"Output files are put into DIR2\n");
    fprintf(stderr,"DETECTOR is either H1, H2, or L1\n");
    fprintf(stderr,"Run - 2 (S2) or 3 (S3)\n");
    fprintf(stderr,"injections? 0 for no, 1 for yes\n");
    fprintf(stderr,"There were argc=%d arguments:\n",argc);
    fflush(stderr);
    for (a=0;a<argc;a++)
			fprintf(stderr,"arg=%d: %s\n",a,argv[a]);
    fflush(stderr);
    return 2;
  }
	#else
  scanf("%d",&jobnum);
  fprintf(stderr,"condor job no %d\n", jobnum);
	#endif 

  irun = atoi(argv[5]);
  if (irun != 2 && irun != 3 && irun!=4){
    fprintf(stderr, "irun = ??? Should be 2 (S2) or 3 (S3)!\n");
    return(1);
  }
   
  iflag = atoi(argv[6]);
  if (iflag != 0 && iflag != 1){
    fprintf(stderr, "iflag has to be 0 or 1. Either regular or hardware injections\n");
    return(0);
  }
      
	/************** BEGIN READING PULSAR PARAMETERS ******************/     

	/* read input file with pulsar information */
  sprintf(psrinput,"%s/pulsar.list",argv[2]);
  psrfp=tryopen(psrinput,"r");
       
  while (!feof(psrfp)){
    /* if iflag = 1 then we are analysing data from a hardware injection.  the format of the input 
		file is slighly different for the hardware injections */
  
  	if (iflag == 1)
    	fscanf(psrfp, "%d %s %lf %lf %s %s %s %lf %lf %lf %s", &pnum, &pname[counter][0], &ra_sec,&dec_sec,
    	&ppmra[0], &ppmdec[0], &pposepoch[0], &pf0, &pf1, &pf2,  &ppepoch[0]); 
  	else{ 
     	fscanf(psrfp, "%d %s %lf %lf %s %s %s %lf %lf %lf %s %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
     	&pnum, &pname[counter][0], &ra_sec, &dec_sec,
     	&ppmra[0], &ppmdec[0], &pposepoch[0], &pf0, &pf1, &pf2,  &ppepoch[0],
		 	&pmodel[0], &pT0, &pTasc, &pPb, &pPbdot, &px, &pxdot, &pw0, &pwdot, &pe,
		 	&pedot, &peps1, &peps1dot, &peps2, &peps2dot, &pxpbdot, &pgamma, &ps, &pM,
		 	&pm2, &pdr, &pdth, &pa0, &pb0); 
  	}
				
		if(feof(psrfp) && counter > 0)
			break;

		fprintf(stderr, "Pulsar name is %s and model is %s.\n", pname[counter], pmodel);

    f0[counter] = pf0;
    f1[counter] = pf1;
    f2[counter] = pf2;
		T0[counter] = pT0;
		Tasc[counter] = pTasc;
		Pb[counter] = pPb;
		Pbdot[counter] = pPbdot;
		x[counter] = px;
		xdot[counter] = pxdot;
		w0[counter] = pw0;
		wdot[counter] = pwdot;
		e[counter] = pe;
		edot[counter] = pedot;
		eps1[counter] = peps1;
		eps1dot[counter] = peps1dot;
		eps2[counter] = peps2;
		eps2dot[counter] = peps2dot;
		xpbdot[counter] = pxpbdot;
		gamma[counter] = pgamma;
		s[counter] = ps;
		M[counter] = pM;
		m2[counter] = pm2;
		dr[counter] = pdr;
		dth[counter] = pdth;
		a0[counter] = pa0;
		b0[counter] = pb0;
		sprintf(model[counter], "%s", pmodel);
		
		/* use LALTDBtoGPS to convert T0 and Tasc to GPS times */
		if(strcmp(model[counter],"*")!=0){
			if(T0[counter] > 0.0){
				T0[counter] = LALTDBtoGPS(T0[counter]);
				fprintf(stderr, "T0 = %.9f\n", T0[counter]); 
			}
			if(Tasc[counter] > 0.0){
				Tasc[counter] = LALTDBtoGPS(Tasc[counter]);
				fprintf(stderr, "Tacs = %.9f\n", Tasc[counter]);
			}
		}
    
		/* if this is a hardware injection then divide the spin parameters provided by 2.0 since
    the rest of the code is especting the pulsar rotation frequency as opposed to the gw frequency */
    if (iflag == 1){
      f0[counter] /= 2.0;
      f1[counter] /= 2.0;
      f2[counter] /= 2.0;
    }
    if ((2.0*f0[counter]) > 1500.0 || (2.0*f0[counter]) < 10.0){
      fprintf(stderr, "f0 outside range for truncated calibration info!\n");
      return(6);
    }

    /* RA and Dec now read in in Radians as conversion is done in preprocessing
		of the table */
    
    RA[counter] = ra_sec;
    DEC[counter] = dec_sec;
    
		/* check if either the position epoch or period epoch are zero and set accordingly */
		posepoch = atof(pposepoch);
		ffepoch = atof(ppepoch);
		
		if(posepoch == 0.0 && ffepoch != 0.0){
			posepoch = ffepoch;
		}
		else if(ffepoch == 0.0 && posepoch != 0.0){
			ffepoch = posepoch;
		}
		else if(ffepoch == 0.0 && posepoch == 0.0){
			fprintf(stderr, "Neither the frequency or position epoch are set for PSR %s.\n", pname[counter]);
			continue;
		}
		
    if (strcmp(ppmra, "*")) { pmra = atof(ppmra);}
      else {pmra = 0.0;}
    if (strcmp(ppmdec, "*")) { pmdec = atof(ppmdec);}
      else {pmdec = 0.0;}
    if (strcmp(pposepoch, "*")) { posepoch = atof(pposepoch);}
      else {posepoch = 0.0;}
    if (strcmp(ppepoch, "*")) { ffepoch = atof(ppepoch);}
      else {ffepoch = 0.0;}
  
    /*adjust RA and DEC to epoch at beg of run */
    if (irun==2){
      REF_MJD = S2_MJD;
      REF_GPS_SECS = S2_GPS_SECS;
    }
    else if (irun==3){
      REF_MJD = S3_MJD;
      REF_GPS_SECS = S3_GPS_SECS;
    }
		else if (irun==4){
    	REF_MJD = S4_MJD;
			REF_GPS_SECS = S4_GPS_SECS;
		}
		else{
      fprintf(stderr, "irun must be either 2, 3 or 4!\n");
      return(1);
    }   
    
    /*if (iflag != 1) { 
      RA[counter] +=  ((REF_MJD - posepoch)/365.26) * pmra;
      DEC[counter] += ((REF_MJD - posepoch)/365.26) * pmdec;}*/
		/* pmra and pmdec are in rad/s */
		if (iflag != 1) { 
      RA[counter] +=  ((REF_MJD - posepoch)*86400.0) * pmra;
      DEC[counter] += ((REF_MJD - posepoch)*86400.0) * pmdec;
		}
 
    posepoch = REF_MJD;
    /*if (ffepoch != 0.0) fepoch_gps = (-REF_MJD + ffepoch)*86400.0 + REF_GPS_SECS;*/
		if (ffepoch != 0.0) fepoch_gps = LALTDBtoGPS(ffepoch);
    else fepoch_gps = REF_GPS_SECS;
   
    fepoch[counter].gpsSeconds = (INT4)floor(fepoch_gps);
    fepoch[counter].gpsNanoSeconds = (INT4)floor((fmod(fepoch_gps, 1.0)*1.e9)); 

		/* if S2 hardware injections then hardwire the frequency epochs here
		- due to fact that I did not have these epochs in MJD format so easier to 
		just put them in manually here */	    
    if (iflag==1 && irun==2)
    {      
      if ((f0[counter]*2.0) < 1285.0) /* thus signal 1*/
      {
        fepoch[counter].gpsSeconds =733967667;
        fepoch[counter].gpsNanoSeconds = 26112310;
      }
      else if ((f0[counter]*2.0) > 1285.0) /* thus signal 2*/
     	{
        fepoch[counter].gpsSeconds =733967751;
        fepoch[counter].gpsNanoSeconds =  522490380;
     	}
     	else
     	{fprintf(stderr, "frequency of signal should not be exactly 1285Hz!\n"); return(0);}    
    }
  
    if (DBUG) fprintf(stderr, "%s\tfepoch.gpsSeconds=%d\nfepoch.gpsNanoSeconds=%d\n", 
        pname[counter], fepoch[counter].gpsSeconds, fepoch[counter].gpsNanoSeconds);

    /* print out pulsar information read to output directory, this file will be used later
    in the ComputePDF code */
    sprintf(psrout[counter],"%s/%s/%s/%s",argv[3],argv[4], pname[counter],pname[counter]);
    fpout3=tryopen(psrout[counter],"w");
    fprintf(fpout3, "ra\t%f\ndec\t%f\npmra\t%f\npmdec\t%f\npepoch\t%f\nf0\t%f\nf1\t%e\nf2\t%e\nfepoch\t%f",
		     RA[counter], DEC[counter], pmra, pmdec, posepoch, f0[counter],
		     f1[counter], f2[counter], ffepoch); 
    fflush(fpout3);
    fclose(fpout3);
    
		counter++;   
  }  /* loop until all the lines are read from the pulsar.list file */
  fclose(psrfp);
  /************** END READING PULSAR PARAMETERS ******************/   
  
  /************** INITIALIZE 'STATIC' VARIABLES BEFORE GOING INTO LOOP OVER PULSARS *******/
  if (!strcmp(argv[4],"H1") || !strcmp(argv[4],"H2")) 
    params.detector = lalCachedDetectors[LALDetectorIndexLHODIFF];      
  else if (   !strcmp(argv[4],"L1"))
   	params.detector = lalCachedDetectors[LALDetectorIndexLLODIFF];
  else if (   !strcmp(argv[4],"GEO"))
   	params.detector = lalCachedDetectors[LALDetectorIndexGEO600DIFF];             
  else
  {
    fprintf(stderr,"Error DETECTOR must be either H1, H1, L1, or GEO\n");
    return 2;
  }  
          
  /* set up three Butterworth 3 order filters before looping over pulsars 
      cut off frequencies fc1, fc2, and fc3 are defined in TargetedPulsars.h */
  wc = tan(LAL_PI * fc1 / SRATE);
  LALCreateCOMPLEX16ZPGFilter( &status, &zpgFilter1, 0, 3 );
  zpgFilter1->poles->data[0].re = wc*SQRT3_2;
  zpgFilter1->poles->data[0].im = wc*0.5;
  zpgFilter1->poles->data[1].re = 0.0;
  zpgFilter1->poles->data[1].im = wc;
  zpgFilter1->poles->data[2].re = -wc*SQRT3_2;
  zpgFilter1->poles->data[2].im = wc*0.5;
  zpgFilter1->gain.re = 0.0;
  zpgFilter1->gain.im = wc*wc*wc;
  LALWToZCOMPLEX16ZPGFilter( &status, zpgFilter1 );
       
  wc = tan(LAL_PI * fc2 / SRATE);
  LALCreateCOMPLEX16ZPGFilter( &status, &zpgFilter2, 0, 3 );
  zpgFilter2->poles->data[0].re = wc*SQRT3_2;
  zpgFilter2->poles->data[0].im = wc*0.5;
  zpgFilter2->poles->data[1].re = 0.0;
  zpgFilter2->poles->data[1].im = wc;
  zpgFilter2->poles->data[2].re = -wc*SQRT3_2;
  zpgFilter2->poles->data[2].im = wc*0.5;
  zpgFilter2->gain.re = 0.0;
  zpgFilter2->gain.im = wc*wc*wc;
  LALWToZCOMPLEX16ZPGFilter( &status, zpgFilter2 );
      
  wc = tan(LAL_PI * fc3 / SRATE);
  LALCreateCOMPLEX16ZPGFilter( &status, &zpgFilter3, 0, 3 );
  zpgFilter3->poles->data[0].re = wc*SQRT3_2;
  zpgFilter3->poles->data[0].im = wc*0.5;
  zpgFilter3->poles->data[1].re = 0.0;
  zpgFilter3->poles->data[1].im = wc;
  zpgFilter3->poles->data[2].re = -wc*SQRT3_2;
  zpgFilter3->poles->data[2].im = wc*0.5;
  zpgFilter3->gain.re = 0.0;
  zpgFilter3->gain.im = wc*wc*wc;
  LALWToZCOMPLEX16ZPGFilter( &status, zpgFilter3 );
     
	for(i_pulsar=0;i_pulsar<counter;i_pulsar++){
	  /* create array of filters, one for each pulsar */
	  iirFilters[i_pulsar].iirFilter1Re = NULL;
	  iirFilters[i_pulsar].iirFilter1Im = NULL;
	  LALCreateREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter1Re, zpgFilter1 );
   	TESTSTATUS(&status);
		LALCreateREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter1Im, zpgFilter1 );
	  TESTSTATUS(&status);
		
	  iirFilters[i_pulsar].iirFilter2Re = NULL;
	  iirFilters[i_pulsar].iirFilter2Im = NULL;
	  LALCreateREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter2Re, zpgFilter2 );
	  TESTSTATUS(&status);
		LALCreateREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter2Im, zpgFilter2 );
		TESTSTATUS(&status);

	  iirFilters[i_pulsar].iirFilter3Re = NULL;
	  iirFilters[i_pulsar].iirFilter3Im = NULL;		 
	  LALCreateREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter3Re, zpgFilter3 );
	  TESTSTATUS(&status);
		LALCreateREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter3Im, zpgFilter3 );
	  TESTSTATUS(&status);
	}
	
	LALDestroyCOMPLEX16ZPGFilter( &status, &zpgFilter1 );
	LALDestroyCOMPLEX16ZPGFilter( &status, &zpgFilter2 );
	LALDestroyCOMPLEX16ZPGFilter( &status, &zpgFilter3 );
	 
	fprintf(stderr, "I've create the array of filters.\n");
			
	/* begin reading calibration information for this pulsar */   
  if (strcmp(argv[4],"GEO")){
    REAL8Vector *H0=NULL, *phaseH0=NULL, *C0=NULL, *phaseC0=NULL;
		REAL8Vector *freqOLG=NULL, *freqSen=NULL;
		 
		/* vectors if det is H2 druing S2 run */
		REAL8Vector *H02=NULL, *phaseH02=NULL, *C02=NULL, *phaseC02=NULL;
		REAL8Vector *freqOLG2=NULL, *freqSen2=NULL;
		/**************************************/
		 
		REAL8 H0temp, phaseH0temp;
		REAL8 C0temp, phaseC0temp;
		
		/* allocate memory for ASQ sensing function and open loop gain */
		LALDCreateVector(&status, &H0pul, MAXPSRS);
		LALDCreateVector(&status, &phaseH0pul, MAXPSRS);
		LALDCreateVector(&status, &C0pul, MAXPSRS);
		LALDCreateVector(&status, &phaseC0pul, MAXPSRS);
		 
		LALDCreateVector(&status, &H0pul2, MAXPSRS);
		LALDCreateVector(&status, &phaseH0pul2, MAXPSRS);
		LALDCreateVector(&status, &C0pul2, MAXPSRS);
		LALDCreateVector(&status, &phaseC0pul2, MAXPSRS);
		
		/* allocate memory to store the alpha/beta calibration parameters */   
    LALDCreateVector(&status, &alpha, MAXLINES);
		LALDCreateVector(&status, &beta, MAXLINES);
   
    /* allocate memory for time stamps for each alpha/beta */    
    LALDCreateVector(&status, &tab, MAXLINES);
			
		/* allocate memory for calibration info */
	  LALDCreateVector(&status, &freqOLG, MAXLINES);
	  LALDCreateVector(&status, &freqSen, MAXLINES);
	  LALDCreateVector(&status, &H0, MAXLINES);
	  LALDCreateVector(&status, &phaseH0, MAXLINES);
	  LALDCreateVector(&status, &C0, MAXLINES);
	  LALDCreateVector(&status, &phaseC0, MAXLINES);
		
		sprintf(olgfile,"%s/calibration/%sopenloopgain.txt",argv[2],argv[4]);    
    fpolg = tryopen(olgfile,"r");
      
		i=0;
		while(!feof(fpolg)){
			fscanf(fpolg, "%lf%lf%lf", &f, &H0temp, &phaseH0temp);
			freqOLG->data[i] = f;
			H0->data[i] = H0temp;
			phaseH0->data[i] = phaseH0temp;
			i++;
			if (i > MAXLINES){ fprintf(stderr,"MAXLINES not big enough!\n");return(4);}
		}
			
		fclose(fpolg);
			
		sprintf(sffile,"%s/calibration/%sASQsensing.txt",argv[2],argv[4]);   
    fpsf = tryopen(sffile,"r");
			
		i=0;
		while(!feof(fpsf)){
			fscanf(fpsf, "%lf%lf%lf", &f, &C0temp, &phaseC0temp);
			freqSen->data[i] = f;
			C0->data[i] = C0temp;
			phaseC0->data[i] = phaseC0temp;
			i++;
			if (i > MAXLINES){ fprintf(stderr,"MAXLINES not big enough!\n");return(4);}
		}
			
		fclose(fpsf);
			
    /* read in alpha beta values */
		sprintf(abfile,"%s/calibration/%salphabeta.txt",argv[2],argv[4]);
  
    fpab = tryopen(abfile,"r");
    fscanf(fpab, "%s%s%s%s%s", &jnkstr1[0], &jnkstr2[0], &jnkstr3[0],
		&jnkstr4[0], &jnkstr5[0]);
    
    i = 0;
    while (!feof(fpab)){
      fscanf(fpab, "%lf%lf%lf%lf", &tg, &ab, &a, &b);
      if (a>=ALPHAMIN && a<=ALPHAMAX){
        tab->data[i] = tg;
        alpha->data[i] = a;
        beta->data[i] = b;
        i++;
  	    if (i > MAXLINES){ fprintf(stderr,"MAXLINES not big enough!\n");return(4);}
      }
    }
    fclose(fpab);
    Nab = i;
    
		fprintf(stderr, "I've read in the calibration info.\n");
    /* end reading calibration information for this pulsar */
	   
	  LALDCreateVector(&status, &H02, MAXLINES);
	  LALDCreateVector(&status, &phaseH02, MAXLINES);
	  LALDCreateVector(&status, &C02, MAXLINES);
	  LALDCreateVector(&status, &phaseC02, MAXLINES);
		LALDCreateVector(&status, &freqOLG2, MAXLINES);
	  LALDCreateVector(&status, &freqSen2, MAXLINES);
		
		/* if det is H2 during S2 run then read in values for second epoch */
		if (!strcmp(argv[4],"H2") && flag == 0 && irun == 2){
		  flag = 1;
	    if (DBUG) fprintf(stderr, "reading second epoch of S2 calibration for H2\n");
       
			sprintf(olgfile,"%s/calibration/%sopenloopgainb.txt",argv[2],argv[4]);
      fpolg = tryopen(olgfile,"r");
       
			i=0;
      while(!feof(fpolg)){  
        fscanf(fpolg,"%lf%lf%lf",&f,&H0temp,&phaseH0temp);
				freqOLG2->data[i] = f;
				H02->data[i] = H0temp;
				phaseH02->data[i] = phaseH0temp;
				i++;
			}
        
      fclose(fpolg);    
        
			sprintf(sffile,"%s/calibration/%sASQsensingb.txt",argv[2],argv[4]);  
      fpsf = tryopen(sffile,"r");    
       
			i = 0;
      while(f<2.0*f0[i_pulsar]){
        fscanf(fpsf,"%lf%lf%lf",&f,&C0temp,&phaseC0temp);
				freqSen2->data[i] = f;
				C02->data[i] = C0temp;
				phaseC02->data[i] = phaseC0temp;
				i++;
      } 
      fclose(fpsf); 
		}
		 
	  /* get calibration values for each pulsar */
	 	for(i_pulsar=0;i_pulsar<counter;i_pulsar++){
	 	  i=0;
     
		  while(freqOLG->data[i]<2.0*f0[i_pulsar]){  
        H0pul->data[i_pulsar] = H0->data[i];
			  phaseH0pul->data[i_pulsar] = phaseH0->data[i];
			  i++;
		  }       
       
      /*if (fabs(f-2.0*f0[i_pulsar]) > 0.1) 
        {fprintf(stderr,"Open loop gain file does not have frequency bin near f0\n");  return(4);}
      */       
      
		  i=0;
      while(freqSen->data[i]<2.0*f0[i_pulsar]){
        C0pul->data[i_pulsar] = C0->data[i];
			  phaseC0pul->data[i_pulsar] = phaseC0->data[i];
			  i++;
		  }
    
		  /* if det is H2 during S2 */
			if (!strcmp(argv[4],"H2") && flag == 0 && irun == 2){
			  while(freqOLG2->data[i]<2.0*f0[i_pulsar]){  
          H0pul2->data[i_pulsar] = H02->data[i];
			    phaseH0pul2->data[i_pulsar] = phaseH02->data[i];
			    i++;
		    }             
      
		    i=0;
        while(freqSen2->data[i]<2.0*f0[i_pulsar]){
          C0pul2->data[i_pulsar] = C02->data[i];
			    phaseC0pul2->data[i_pulsar] = phaseC02->data[i];
			    i++;
		    }
			}
		
      /* if (fabs(f-2.0*f0[i_pulsar]) > 0.1) 
        {fprintf(stderr,"Sensing function file does not have frequency bin near f0\n"); return(4);}
      */
		}
		 
		LALDDestroyVector(&status, &freqOLG);
		LALDDestroyVector(&status, &freqSen);
		LALDDestroyVector(&status, &H0);
		LALDDestroyVector(&status, &phaseH0);
		LALDDestroyVector(&status, &C0);
		LALDDestroyVector(&status, &phaseC0);
		 
		LALDDestroyVector(&status, &freqOLG2);
		LALDDestroyVector(&status, &freqSen2);
		LALDDestroyVector(&status, &H02);
		LALDDestroyVector(&status, &phaseH02);
		LALDDestroyVector(&status, &C02);
		LALDDestroyVector(&status, &phaseC02);        
	}
	
	fprintf(stderr, "I've set the calibration values for each pulsar.\n");
	
	/* set up LALBarycenter */
  edat = (EphemerisData *)LALMalloc(sizeof(EphemerisData));    
  sprintf(earthfile,"%s/earth03-06.dat",argv[2]);
  sprintf(sunfile,"%s/sun03-06.dat",argv[2]);
  (*edat).ephiles.earthEphemeris = earthfile;
  (*edat).ephiles.sunEphemeris = sunfile;	 
  LALInitBarycenter(&status, edat);
  if(status.statusCode){
    fprintf(stderr, "Unexpectedly got error code %d and message %s\n", 
    status.statusCode, status.statusDescription);
    return 0;
  } 
  params.edat = edat;
	
	/* construct frame channel name */   
  if (!strcmp(argv[4],"GEO")){
   	/* sprintf(chname,"G1:DER_DATA_HP"); channel name when looking at injection
		from June 04 */
		sprintf(chname,"G1:DER_DATA_H");
  }
	else 
    sprintf(chname,"%s:LSC-AS_Q",argv[4]);
   	
	/* set up the frame reading structures */   
  {
    char timelist[256];
    INT4 tbasei;
    FILE *stfp=NULL;
    sprintf(timelist,"%s/data%s/jobtimes.%05d",argv[2],argv[4],jobnum);
    stfp=tryopen(timelist,"r"); 
    
		while (2==fscanf(stfp,"%d %d",starts+nstarts,&tbasei)){
      if (nstarts==0)
        tbase=tbasei;
      else if (tbasei!=tbase){
        fprintf(stderr,"File %s contains inconsistent SFT time baselines %d != %d\n",timelist,tbasei,tbase);
	      fflush(stderr);
	      return 3;
      }
      nstarts++;
      if (nstarts>=MAXSTART){
	      fprintf(stderr,"More than MAXSTART=%d lines in file: %s.  Increase MAXSTART\n",
	      MAXSTART,timelist);
	      fflush(stderr);
	      return 3;
      }
    }
    
		if (!nstarts){
      fprintf(stderr,"File %s didn't contain any valid lines!\n",timelist);
      fflush(stderr);
      return 3;
    }
    fclose(stfp);
  }
	  
  tbase = TBASE;
  npts = SRATE*tbase;
	 
  sprintf(framelist,"%s/data%s/jobdata.%05d.ffl",argv[2],argv[4],jobnum);
  /*opencount=0;
  while (!(frfile = FrFileINew(framelist))){
    fprintf(stderr, "Couldnt open frame file %s\n", framelist);
    fflush(stderr);
    if (opencount++<10)
      sleep(10);
    else
      return(3);
  } */
     
	fprintf(stderr, "I've set up the frame reading input structures.\n");
	/* i've set up the frame reading structures */
	
	/* set up output files for each pulsar - only once */
	for(i_pulsar=0;i_pulsar<counter;i_pulsar++){
		sprintf(outfile,"%s/%s/%s/outfine.%s_%s.%d",argv[3],argv[4],pname[i_pulsar],
pname[i_pulsar],argv[4],starts[1]);
		fpout[i_pulsar] = tryopen(outfile, "w");
	}
	
	k=0;
	/* loop to read in data */
	for (count=0;count<nstarts;count++){
    INT4 startEpoch;
		char *shortlist=NULL; /* short list of frame files */
		
		epoch.gpsSeconds=starts[count];
    epoch.gpsNanoSeconds=0;    
    
		/* allocate memory for 60 second chucks of raw data*/
  	input.V.data = NULL;
  	LALDCreateVector( &status, &input.V.data, npts);  
  	input.V.deltaT = 1.0/SRATE;
  	input.V.data->length = npts;
		
		input.V.epoch = epoch;
    
		if(k==1)
			startEpoch = epoch.gpsSeconds;
 
		/* make short list of frames files from framelist */
		shortlist = MakeFrameFileList(framelist, tbase ,starts[count]);
 		/* open frame file */
		
		/* in case my MakeFileFrameList code returns null then open the whole list*/
		if(shortlist != NULL){
			if(!(frfile = FrFileINew(shortlist))){
				fprintf(stderr, "Couldnt open frame file %s\n", shortlist);
				return(3);
			}
		}
		else if(shortlist == NULL){
			if(!(frfile = FrFileINew(framelist))){
				fprintf(stderr, "Couldnt open frame file %s\n", framelist);
				return(3);
			}
		}
 
    frvect = FrFileIGetVAdc(frfile, chname, (double)epoch.gpsSeconds, (double)tbase, 0);
    if (frvect == NULL){
      fprintf(stderr, "Data between times %d and %d not present\n",epoch.gpsSeconds,epoch.gpsSeconds+tbase);
      fflush(stderr);
			
			/* destroy input data vector */
			LALDDestroyVector(&status, &input.V.data);
			FrFileIEnd(frfile);
      continue;
    }
    if (strcmp(argv[4],"GEO")){     
      if ( frvect->type != FR_VECT_4R ){
        fprintf(stderr, "Wrong data type (not FR_VECT_4R) found in frame!\n" );
        fflush(stderr);
        return(5);
      }   
    } 
    else /* if GEO */{
      if ( frvect->type != FR_VECT_8R ){
        fprintf(stderr, "Wrong data type (not FR_VECT_8R) found in frame!\n" );
        fflush(stderr);
        return(5);
      }   
    }
    
    if (frvect->next){
      fprintf(stderr, "Data between times %d and %d had a gap\n",epoch.gpsSeconds,epoch.gpsSeconds+tbase);
      fflush(stderr);
      FrVectFree(frvect);
      frvect=NULL;
      
			/* destroy input data vector */
			LALDDestroyVector(&status, &input.V.data);
			FrFileIEnd(frfile);
			continue;
    }    
   
    if (strcmp(argv[4],"GEO"))  /* if LIGO */{     
      for (i=0;i<npts;i++) 
        input.V.data->data[i] = frvect->dataF[i];
    }
    else /* if GEO */{
      for (i=0;i<npts;i++){ 
        input.V.data->data[i] = frvect->dataD[i];
			}
    }
    fprintf(stderr, "In loop %d of %d.\n", count+1, nstarts);
    /* have read in frame data */
	 
    FrVectFree(frvect);
    frvect=NULL;
		
 		/*********** begin loop over pulsars ***************/ 
  	for (i_pulsar=0;i_pulsar<counter;i_pulsar++){
  
 			/* set up input paramters for LALHeterodyneToPulsar for pulsar i_pulsar */ 
   		input.V.f0 = 0.0; /* base frequency */
   		input.f0 = 2.0*f0[i_pulsar];
   		input.f1 = 2.0*f1[i_pulsar];
   		input.f2 = 2.0*f2[i_pulsar]; 
   		input.source.longitude = RA[i_pulsar];
   		input.source.latitude = DEC[i_pulsar];
   		input.fEpochGPS = (double)fepoch[i_pulsar].gpsSeconds + (double)fepoch[i_pulsar].gpsNanoSeconds*1.0e-9;
   		input.V.deltaT = 1.0/SRATE; 
			
	 		/* set correct filters for pulsar */
	 		params.iirFilter1Re = iirFilters[i_pulsar].iirFilter1Re;
   		params.iirFilter1Im = iirFilters[i_pulsar].iirFilter1Im;
	 
	 		params.iirFilter2Re = iirFilters[i_pulsar].iirFilter2Re;
   		params.iirFilter2Im = iirFilters[i_pulsar].iirFilter2Im;
	 
	 		params.iirFilter3Re = iirFilters[i_pulsar].iirFilter3Re;
   		params.iirFilter3Im = iirFilters[i_pulsar].iirFilter3Im;
			
	 		/* set binary params */
	 		params.binaryInput.tbflag = "GPS";
	 		params.binaryParams.T0 = T0[i_pulsar];
	 		params.binaryParams.Tasc = Tasc[i_pulsar];
	 		params.binaryParams.Pb = Pb[i_pulsar];
	 		params.binaryParams.Pbdot = Pbdot[i_pulsar];
	 		params.binaryParams.x = x[i_pulsar];
	 		params.binaryParams.xdot = xdot[i_pulsar];
	 		params.binaryParams.w0 = w0[i_pulsar];
	 		params.binaryParams.wdot = wdot[i_pulsar];
	 		params.binaryParams.e = e[i_pulsar];
	 		params.binaryParams.edot = edot[i_pulsar];
	 		params.binaryParams.eps1 = eps1[i_pulsar];
	 		params.binaryParams.eps1dot = eps1dot[i_pulsar];
	 		params.binaryParams.eps2 = eps2[i_pulsar];
	 		params.binaryParams.eps2dot = eps2dot[i_pulsar];
	 		params.binaryParams.xpbdot = xpbdot[i_pulsar];
	 		params.binaryParams.gamma = gamma[i_pulsar];
	 		params.binaryParams.s = s[i_pulsar];
	 		params.binaryParams.M = M[i_pulsar];
	 		params.binaryParams.m2 = m2[i_pulsar];
	 		params.binaryParams.dr = dr[i_pulsar];
	 		params.binaryParams.dth = dth[i_pulsar];
	 		params.binaryParams.a0 = a0[i_pulsar];
	 		params.binaryParams.b0 = b0[i_pulsar];
			params.binaryParams.model = model[i_pulsar]; 
    	params.binaryParams.nEll = 0;
			
			/* if the binary model is ELL1 then check whether we have eps derivs or w
			time derivs */
			if(strstr(params.binaryParams.model, "ELL1") != NULL){
				if(params.binaryParams.edot != 0.0 || params.binaryParams.wdot != 0.0){
					params.binaryParams.nEll = 1;
				}
			} 
			 
	 		/* set calibration data */
	 		if(strcmp(argv[4],"GEO")){
	   		calParams.H0 = H0pul->data[i_pulsar];
		 		calParams.phaseH0 = phaseH0pul->data[i_pulsar];
		 		calParams.C0 = C0pul->data[i_pulsar];
		 		calParams.phaseC0 = phaseC0pul->data[i_pulsar];
	 		}

			/******* CALL LALHeterodyneToPulsar ************/     
   		LALHeterodyneToPulsar( &status, &output, &input, &params ); 
   		if(status.statusCode){
     		fprintf(stderr, "Unexpectedly got error code %d and message %s\n", 
     		status.statusCode, status.statusDescription);
     		return 0;
   		}      

   		if (k>0)/* due to  to iir step response do not write first minute after start */
   		{   
       	B.re =  output.B.re;
       	B.im =  output.B.im; 
       	tt = (double)epoch.gpsSeconds + (double)epoch.gpsNanoSeconds*1.0e-9;  
   		
   			/**** begin calibrating Bk's and writing Bk's to file****/   
 				/* if (DBUG) fprintf(stderr, "opening output file: epoch %d, psrname %s\n", epoch.gpsSeconds, pname[i_pulsar]); */
 
 				/* create output file name and try to open (will fail if output director doesn't exist) */
 			
/*sprintf(outfile,"%s/%s/%s/outfine.%s_%s.%d",argv[3],argv[4],pname[i_pulsar],
pname[i_pulsar],argv[4],startEpoch);*/
 				
				/* if this is first output just create output file (or overwrite 
					 an existing file, else just amend to the, already created output */
				
				/* fpout=tryopen(outfile,"a"); */   
 
 				flag=0;
 				/* only calibrate if we have LIGO data */
 				if (strcmp(argv[4],"GEO")){    
      		/* calibration information has already been read. however, there are two reference 
      	 	epochs for H2 S2.  so if this is H2 S2 and we need to use the second epoch, then
	 				read the second reference ASQ sensing and open loop gain */ 
      		if (!strcmp(argv[4],"H2") && tt > 731849043.0 && flag == 0 && irun == 2)
      		{								        
        		flag =1 ;
						calParams.H0 = H0pul2->data[i_pulsar]; 
        		calParams.phaseH0 = phaseH0pul2->data[i_pulsar];
        		calParams.C0 = C0pul2->data[i_pulsar]; 
        		calParams.phaseC0 = phaseC0pul2->data[i_pulsar];       
      		} 
      		/* find first alpha/beta time stamp that is less than 60 seconds away from Bk time stamp */    
      		for (j=jtemp;j<Nab;j++) {
        		if(fabs(tt - tab->data[j])<=30){
  	  				calInput.alpha = alpha->data[j];
	  					calInput.beta = beta->data[j];
							calInput.B.re = B.re;
							calInput.B.im = B.im;
		
          		/* calibrate Bk's and write to file */
	  					LALCalibrateFineHeterodyne(&status, &calOutput, &calInput, &calParams);
	  					fprintf(fpout[i_pulsar],"%f\t%e\t%e\n", tt, calOutput.B.re,
calOutput.B.im);
							jtemp = j;
          		break;
        		}
      		} /* loop over alpha/beta */
  			} /* if LIGO data then calibrate */
  
  			else /* if GEO no need to calibrate just write to file */
  			{  
       		fprintf(fpout[i_pulsar],"%f\t%e\t%e\n", tt, B.re, B.im);              
  			} 
				/* fclose(fpout); */
 
			} /* if (k>0) close */
  	      
 			/**** end calibrating Bk's and writing Bk's to file****/   
		} /*********** end loop over pulsars ***************/
		
		/* destroy input data vector */
		LALDDestroyVector(&status, &input.V.data);
		FrFileIEnd(frfile);
		k++;
	} /***** end extracting from frvect and heterodyning the data *****/
  	
	LALFree(edat->ephemE);
  LALFree(edat->ephemS);
  LALFree(edat);
	
  /* free memory and check for memory leaks */
  
	if(strcmp(argv[4],"GEO")){
		LALDDestroyVector(&status, &tab);   
		LALDDestroyVector(&status, &alpha);   
		LALDDestroyVector(&status, &beta);
		LALDDestroyVector(&status, &H0pul);
		LALDDestroyVector(&status, &phaseH0pul);
		LALDDestroyVector(&status, &C0pul);
		LALDDestroyVector(&status, &phaseC0pul);
		LALDDestroyVector(&status, &H0pul2);
		LALDDestroyVector(&status, &phaseH0pul2);
		LALDDestroyVector(&status, &C0pul2);
		LALDDestroyVector(&status, &phaseC0pul2);
	}
	
	for(i_pulsar=0;i_pulsar<counter;i_pulsar++){	   
  	LALDestroyREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter1Re);
  	LALDestroyREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter1Im ); 
  	LALDestroyREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter2Re );
  	LALDestroyREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter2Im );   
  	LALDestroyREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter3Re );
  	LALDestroyREAL8IIRFilter( &status, &iirFilters[i_pulsar].iirFilter3Im );
		
		/* close pulsar output files */
		fclose(fpout[i_pulsar]);
  }

	fprintf(stderr, "I've destroyed the filters.\n");
	 
  TESTSTATUS( &status );
  LALCheckMemoryLeaks();
  return 0;
} /* main */

char *MakeFrameFileList(char *framelist, int length, int tstart){
	FILE *fp;
	char *shortlist=NULL;
	char FrFileNames[10000][256];
	int startTimes[10000], dur[10000]; /* start time and duration of the frame
files */
	int i=0, j=0, k=0;
	
	/* open up frame list */
	if((fp = fopen(framelist, "r"))==NULL){
		fprintf(stderr, "Error opening frame list file!\n");
		return NULL;
	}
	
	i=0;
	/* read in all the frame filenames */ 
	while(!feof(fp)){
		fscanf(fp, "%s", &FrFileNames[i]);
		
		/* file names should end with GPS-length.gwf - *********-***.gwf, so we can
			 locate the first '-' before the end of the file name */
		/* locate the start time and length of the data */
		for (j=strlen(FrFileNames[i])-1; j>=0; j--)
			if (FrFileNames[i][j]=='-') break;
		
		/* get start times and durations (assumes the GPS time is 9 digits long) */
		startTimes[i]=atoi(FrFileNames[i]+j-9);
    dur[i]=atoi(FrFileNames[i]+j+1);
			
		if(feof(fp)) break;
		i++;
	}
	
	fclose(fp);
	
	i--;
	
	/* find frame file that matches our tstart */
	if(tstart < startTimes[0] || tstart >= startTimes[i]){
		fprintf(stderr, "The start time is before that of the first frame or \
after that of the last frame!\n");
		return NULL;
	}
	
	for(j=0;j<i;j++){
		if(tstart >= startTimes[j] && tstart < startTimes[j]+dur[j]){
			shortlist = FrFileNames[j];
			break;
		}
	}
	
	k=1;
	/* add all the files needed to shortlist to cover length */
	while(startTimes[j+k] <= tstart + length){
		sprintf(shortlist, "%s %s", shortlist, FrFileNames[j+k]);
		k++;
	}
	
	return shortlist;
}
